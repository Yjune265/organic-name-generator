/**
 * stereo.js - CIP priorities, R/S and E/Z.
 *
 * The molecule graph carries raw geometry:
 *   atom.stereo = {order: [n0, n1, n2, n3], clockwise: bool}
 *     "looking from n0 towards the centre, n1 -> n2 -> n3 turn clockwise".
 *     An implicit hydrogen is written as -1 in `order`.
 *   bond.stereo = {refA, refB, same: bool}
 *     refA (on bond.a) and refB (on bond.b) lie on the same side of the
 *     double bond when `same` is true.
 * This module turns that into the descriptors a name needs.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};
var MAX_DEPTH = 9;

function atomicNumber(mol, idx) {
  if (idx < 0) return 1;                       // implicit hydrogen
  var info = ONG.ELEMENTS[mol.atoms[idx].element];
  return info ? info.z : 0;
}

/* ------------------------------------------------------------------ *
 * Hierarchical digraph
 * ------------------------------------------------------------------ */

function node(idx, from, dup, z) {
  return {idx: idx, from: from, dup: dup, z: z};
}

var PHANTOM = {idx: -2, from: -2, dup: true, z: 0};

/**
 * Children of a digraph node. Multiple bonds contribute duplicate atoms
 * (which themselves have no children), hydrogens are explicit leaves.
 */
function childrenOf(mol, n) {
  if (n.dup || n.idx < 0) return [];
  var out = [], links = mol.links(n.idx), k;
  for (k = 0; k < links.length; k++) {
    var l = links[k], w = l.atom;
    if (mol.atoms[w].element !== 'H' && w !== n.from) {
      out.push(node(w, n.idx, false, atomicNumber(mol, w)));
    } else if (mol.atoms[w].element === 'H' && w !== n.from) {
      out.push(node(w, n.idx, false, 1));
    }
    for (var d = 1; d < l.bond.order; d++) {
      out.push(node(w, n.idx, true, atomicNumber(mol, w)));   // duplicate atom
    }
  }
  for (var h = 0; h < mol.implicitH(n.idx); h++) out.push(node(-1, n.idx, false, 1));
  return out;
}

/**
 * CIP comparison of two branches. Returns 1 when `a` outranks `b`.
 * Rule 1a (atomic number) and 1b (duplicate atoms) are applied; ties that
 * survive are reported as 0 so the caller can refuse to guess.
 */
function compareNodes(mol, a, b, depth, memo) {
  if (a.z !== b.z) return a.z > b.z ? 1 : -1;
  if (depth <= 0) return 0;
  if (a.dup && b.dup) return 0;
  if (a.dup !== b.dup) return a.dup ? -1 : 1;     // a real atom outranks its duplicate
  if (a.idx === b.idx && a.from === b.from) return 0;

  var key = a.idx + ',' + a.from + ',' + b.idx + ',' + b.from + ',' + depth;
  if (memo[key] !== undefined) return memo[key];
  memo[key] = 0;                                   // guards cyclic re-entry

  var ca = sortedChildren(mol, a, depth, memo);
  var cb = sortedChildren(mol, b, depth, memo);
  var result = 0;
  var n = Math.max(ca.length, cb.length);
  for (var i = 0; i < n && !result; i++) {
    var x = ca[i] || PHANTOM, y = cb[i] || PHANTOM;
    if (x.z !== y.z) { result = x.z > y.z ? 1 : -1; break; }
  }
  if (!result) {
    for (var j = 0; j < n && !result; j++) {
      result = compareNodes(mol, ca[j] || PHANTOM, cb[j] || PHANTOM, depth - 1, memo);
    }
  }
  memo[key] = result;
  return result;
}

function sortedChildren(mol, n, depth, memo) {
  var kids = childrenOf(mol, n);
  kids.sort(function (x, y) { return compareNodes(mol, y, x, depth - 1, memo); });
  return kids;
}

/**
 * Rank an atom's neighbours by CIP priority.
 * Returns {order: [neighbour indices, highest first], resolved: bool}.
 */
function rankNeighbors(mol, centre, neighbors) {
  var memo = {};
  var nodes = neighbors.map(function (j) {
    return j < 0 ? node(-1, centre, false, 1) : node(j, centre, false, atomicNumber(mol, j));
  });
  var sorted = nodes.slice().sort(function (x, y) {
    return compareNodes(mol, y, x, MAX_DEPTH, memo);
  });
  var resolved = true;
  for (var i = 1; i < sorted.length; i++) {
    if (compareNodes(mol, sorted[i - 1], sorted[i], MAX_DEPTH, memo) === 0) resolved = false;
  }
  return {order: sorted.map(function (n) { return n.idx; }), resolved: resolved};
}

/* ------------------------------------------------------------------ *
 * R / S
 * ------------------------------------------------------------------ */

/** Neighbours of a potential stereocentre, implicit H written as -1. */
function centreNeighbors(mol, i) {
  var list = mol.neighbors(i).slice();
  for (var h = 0; h < mol.implicitH(i); h++) list.push(-1);
  return list;
}

/** Parity of the permutation taking `from` to `to` (both same members). */
function permutationSign(from, to) {
  var arr = from.slice(), sign = 1;
  for (var i = 0; i < to.length; i++) {
    if (arr[i] === to[i]) continue;
    var j = arr.indexOf(to[i], i);
    if (j < 0) return 0;
    var tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
    sign = -sign;
  }
  return sign;
}

/** True when the four neighbours are all different by CIP. */
function isStereocentre(mol, i) {
  var a = mol.atoms[i];
  if (a.element !== 'C' && a.element !== 'Si' && a.element !== 'S' && a.element !== 'N') return false;
  var neighbors = centreNeighbors(mol, i);
  if (neighbors.length !== 4) return false;
  if (mol.bondOrderSum(i) !== mol.degree(i)) return false;   // must be saturated
  if (a.element === 'N' && !a.charge) return false;          // amine nitrogen inverts
  var ranked = rankNeighbors(mol, i, neighbors);
  return ranked.resolved;
}

/** 'R', 'S' or null. */
function descriptorFor(mol, i) {
  var a = mol.atoms[i];
  if (!a.stereo) return null;
  var neighbors = centreNeighbors(mol, i);
  if (neighbors.length !== 4) return null;
  var ranked = rankNeighbors(mol, i, neighbors);
  if (!ranked.resolved) return null;
  var sign = permutationSign(a.stereo.order, ranked.order);
  if (!sign) return null;
  var clockwise = a.stereo.clockwise ? 1 : -1;
  // sign' = +1 means: seen from the top-ranked neighbour, 2 -> 3 -> 4 runs
  // clockwise, which is exactly R once the lowest one points away.
  return (clockwise * sign) > 0 ? 'R' : 'S';
}

/* ------------------------------------------------------------------ *
 * E / Z
 * ------------------------------------------------------------------ */

/** True when both ends of the double bond carry two different groups. */
function isStereoDoubleBond(mol, bond) {
  if (bond.order !== 2) return false;
  if (mol.isRingBond(bond)) {
    var ring = smallestRingWith(mol, bond);
    if (ring && ring.length < 8) return false;    // small rings cannot be E
  }
  return endIsDefined(mol, bond.a, bond.b) && endIsDefined(mol, bond.b, bond.a);
}

function smallestRingWith(mol, bond) {
  var rings = mol.rings(), best = null;
  rings.forEach(function (ring) {
    for (var i = 0; i < ring.length; i++) {
      if (mol.bondBetween(ring[i], ring[(i + 1) % ring.length]) === bond) {
        if (!best || ring.length < best.length) best = ring;
      }
    }
  });
  return best;
}

function endIsDefined(mol, end, other) {
  var subs = mol.neighbors(end).filter(function (j) { return j !== other; });
  for (var h = 0; h < mol.implicitH(end); h++) subs.push(-1);
  if (subs.length !== 2) return false;
  var ranked = rankNeighbors(mol, end, subs);
  return ranked.resolved;
}

/** 'E', 'Z' or null. */
function bondDescriptor(mol, bond) {
  if (!bond.stereo || bond.order !== 2) return null;
  var aSubs = mol.neighbors(bond.a).filter(function (j) { return j !== bond.b; });
  var bSubs = mol.neighbors(bond.b).filter(function (j) { return j !== bond.a; });
  for (var h = 0; h < mol.implicitH(bond.a); h++) aSubs.push(-1);
  for (var g = 0; g < mol.implicitH(bond.b); g++) bSubs.push(-1);
  if (aSubs.length !== 2 || bSubs.length !== 2) return null;

  var rankedA = rankNeighbors(mol, bond.a, aSubs);
  var rankedB = rankNeighbors(mol, bond.b, bSubs);
  if (!rankedA.resolved || !rankedB.resolved) return null;

  var topA = rankedA.order[0], topB = rankedB.order[0];
  // Translate the stored reference pair into "are the top two on one side?"
  var same = bond.stereo.same;
  if (bond.stereo.refA !== topA) same = !same;
  if (bond.stereo.refB !== topB) same = !same;
  return same ? 'Z' : 'E';
}

/* ------------------------------------------------------------------ *
 * Whole-molecule summary
 * ------------------------------------------------------------------ */

/**
 * Returns {atoms: {idx: 'R'|'S'}, bonds: {bondIdx: 'E'|'Z'},
 *          undefinedCentres: [idx], undefinedBonds: [bondIdx], any: bool}
 */
function analyse(mol) {
  var result = {atoms: {}, bonds: {}, undefinedCentres: [], undefinedBonds: [], any: false};
  var i;
  for (i = 0; i < mol.atoms.length; i++) {
    if (mol.atoms[i].stereo) {
      var d = descriptorFor(mol, i);
      if (d) { result.atoms[i] = d; result.any = true; continue; }
    }
    if (isStereocentre(mol, i)) result.undefinedCentres.push(i);
  }
  for (i = 0; i < mol.bonds.length; i++) {
    var bond = mol.bonds[i];
    if (bond.stereo) {
      var e = bondDescriptor(mol, bond);
      if (e) { result.bonds[i] = e; result.any = true; continue; }
    }
    if (isStereoDoubleBond(mol, bond)) result.undefinedBonds.push(i);
  }
  return result;
}

ONG.stereo = {
  analyse: analyse,
  descriptorFor: descriptorFor,
  bondDescriptor: bondDescriptor,
  isStereocentre: isStereocentre,
  isStereoDoubleBond: isStereoDoubleBond,
  rankNeighbors: rankNeighbors,
  centreNeighbors: centreNeighbors,
  permutationSign: permutationSign
};

if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
