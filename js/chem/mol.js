/**
 * mol.js - molecule graph model.
 *
 * Everything downstream (SMILES parsing, ring perception, IUPAC naming,
 * structure lookup) is built on the plain graph defined here.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

/* ------------------------------------------------------------------ *
 * Element data
 * ------------------------------------------------------------------ */

// val: allowed valences, smallest first. Empty -> never carries implicit H.
var ELEMENTS = {
  H:  {z: 1,  w: 1.008,   val: [1]},
  B:  {z: 5,  w: 10.81,   val: [3]},
  C:  {z: 6,  w: 12.011,  val: [4]},
  N:  {z: 7,  w: 14.007,  val: [3]},
  O:  {z: 8,  w: 15.999,  val: [2]},
  F:  {z: 9,  w: 18.998,  val: [1]},
  Si: {z: 14, w: 28.085,  val: [4]},
  P:  {z: 15, w: 30.974,  val: [3, 5]},
  S:  {z: 16, w: 32.06,   val: [2, 4, 6]},
  Cl: {z: 17, w: 35.45,   val: [1]},
  Se: {z: 34, w: 78.971,  val: [2, 4, 6]},
  Br: {z: 35, w: 79.904,  val: [1]},
  I:  {z: 53, w: 126.904, val: [1]},
  Li: {z: 3,  w: 6.94,    val: []},
  Na: {z: 11, w: 22.990,  val: []},
  Mg: {z: 12, w: 24.305,  val: []},
  K:  {z: 19, w: 39.098,  val: []},
  Ca: {z: 20, w: 40.078,  val: []},
  Fe: {z: 26, w: 55.845,  val: []},
  Zn: {z: 30, w: 65.38,   val: []}
};

var HALOGENS = {F: 1, Cl: 1, Br: 1, I: 1};

/** Valence correction for a charged atom. */
function chargeAdjust(symbol, charge) {
  if (!charge) return 0;
  if (symbol === 'C' || symbol === 'Si') return -Math.abs(charge); // no lone pair to donate
  if (symbol === 'B') return -charge;                              // electron deficient
  return charge;                                                   // N+, O-, ...
}

/* ------------------------------------------------------------------ *
 * Mol
 * ------------------------------------------------------------------ */

/**
 * atom: {element, charge, explicitH|null, aromatic, idx}
 *   explicitH === null means "work it out from the valence rules".
 * bond: {a, b, order}   order 1|2|3 (aromatic input is kekulized on the way in)
 */
function Mol() {
  this.atoms = [];
  this.bonds = [];
  this._adj = [];      // atom idx -> [{atom, bond}]
  this._cache = {};
}

Mol.prototype.addAtom = function (element, opts) {
  opts = opts || {};
  var a = {
    idx: this.atoms.length,
    element: element,
    charge: opts.charge || 0,
    explicitH: (opts.explicitH === undefined || opts.explicitH === null) ? null : opts.explicitH,
    isotope: opts.isotope || null,
    aromatic: false
  };
  this.atoms.push(a);
  this._adj.push([]);
  this._cache = {};
  return a.idx;
};

Mol.prototype.addBond = function (a, b, order) {
  if (a === b) return -1;
  var existing = this.bondBetween(a, b);
  if (existing) { existing.order = Math.max(existing.order, order || 1); return existing.idx; }
  var bond = {idx: this.bonds.length, a: a, b: b, order: order || 1, aromatic: false};
  this.bonds.push(bond);
  this._adj[a].push({atom: b, bond: bond});
  this._adj[b].push({atom: a, bond: bond});
  this._cache = {};
  return bond.idx;
};

Mol.prototype.neighbors = function (i) {
  return this._adj[i].map(function (e) { return e.atom; });
};

Mol.prototype.links = function (i) { return this._adj[i]; };

Mol.prototype.degree = function (i) { return this._adj[i].length; };

Mol.prototype.bondBetween = function (a, b) {
  var links = this._adj[a];
  for (var i = 0; i < links.length; i++) if (links[i].atom === b) return links[i].bond;
  return null;
};

Mol.prototype.bondOrder = function (a, b) {
  var bond = this.bondBetween(a, b);
  return bond ? bond.order : 0;
};

Mol.prototype.other = function (bond, i) { return bond.a === i ? bond.b : bond.a; };

Mol.prototype.bondOrderSum = function (i) {
  var sum = 0, links = this._adj[i];
  for (var k = 0; k < links.length; k++) sum += links[k].bond.order;
  return sum;
};

/** Hydrogens that are not drawn as their own atom. */
Mol.prototype.implicitH = function (i) {
  var a = this.atoms[i];
  if (a.explicitH !== null) return a.explicitH;
  var info = ELEMENTS[a.element];
  if (!info || !info.val.length) return 0;
  var sum = this.bondOrderSum(i);
  var adj = chargeAdjust(a.element, a.charge);
  for (var k = 0; k < info.val.length; k++) {
    var v = info.val[k] + adj;
    if (v >= sum) return v - sum;
  }
  return 0;
};

/** Implicit hydrogens plus any drawn as explicit H atoms. */
Mol.prototype.totalH = function (i) {
  var n = this.implicitH(i), nb = this.neighbors(i);
  for (var k = 0; k < nb.length; k++) if (this.atoms[nb[k]].element === 'H') n++;
  return n;
};

Mol.prototype.isHalogen = function (i) { return !!HALOGENS[this.atoms[i].element]; };

/** Neighbours excluding plain hydrogen atoms. */
Mol.prototype.heavyNeighbors = function (i) {
  var self = this;
  return this.neighbors(i).filter(function (j) { return self.atoms[j].element !== 'H'; });
};

Mol.prototype.formula = function () {
  var counts = {}, i;
  for (i = 0; i < this.atoms.length; i++) {
    var el = this.atoms[i].element;
    counts[el] = (counts[el] || 0) + 1;
    if (el !== 'H') counts.H = (counts.H || 0) + this.implicitH(i);
  }
  // Hill order: C, H, then alphabetical
  var keys = Object.keys(counts).filter(function (k) { return k !== 'C' && k !== 'H'; }).sort();
  var order = [];
  if (counts.C) order.push('C');
  if (counts.H) order.push('H');
  order = order.concat(keys);
  return order.map(function (k) {
    return k + (counts[k] > 1 ? counts[k] : '');
  }).join('');
};

Mol.prototype.molecularWeight = function () {
  var w = 0;
  for (var i = 0; i < this.atoms.length; i++) {
    var info = ELEMENTS[this.atoms[i].element];
    if (!info) continue;
    w += info.w;
    if (this.atoms[i].element !== 'H') w += this.implicitH(i) * ELEMENTS.H.w;
  }
  return w;
};

/** Connected components, as arrays of atom indices. */
Mol.prototype.components = function () {
  if (this._cache.components) return this._cache.components;
  var seen = [], out = [], i;
  for (i = 0; i < this.atoms.length; i++) {
    if (seen[i]) continue;
    var stack = [i], comp = [];
    seen[i] = true;
    while (stack.length) {
      var v = stack.pop();
      comp.push(v);
      var nb = this.neighbors(v);
      for (var k = 0; k < nb.length; k++) if (!seen[nb[k]]) { seen[nb[k]] = true; stack.push(nb[k]); }
    }
    comp.sort(function (x, y) { return x - y; });
    out.push(comp);
  }
  this._cache.components = out;
  return out;
};

/** Sub-molecule holding the given atoms; returns {mol, map, rmap}. */
Mol.prototype.subMol = function (atomIdxs) {
  var sub = new Mol(), map = {}, i;
  for (i = 0; i < atomIdxs.length; i++) {
    var a = this.atoms[atomIdxs[i]];
    map[atomIdxs[i]] = sub.addAtom(a.element, {charge: a.charge, explicitH: a.explicitH});
  }
  // Stereochemistry refers to atoms by index, so it has to be remapped.
  function remap(j) { return j < 0 ? j : (map[j] === undefined ? null : map[j]); }
  for (i = 0; i < atomIdxs.length; i++) {
    var src = this.atoms[atomIdxs[i]];
    if (!src.stereo) continue;
    var order = src.stereo.order.map(remap);
    if (order.some(function (x) { return x === null; })) continue;
    sub.atoms[map[atomIdxs[i]]].stereo = {order: order, clockwise: src.stereo.clockwise};
  }
  for (i = 0; i < this.bonds.length; i++) {
    var b = this.bonds[i];
    if (map[b.a] === undefined || map[b.b] === undefined) continue;
    var idx = sub.addBond(map[b.a], map[b.b], b.order);
    if (b.stereo && idx >= 0) {
      var refA = remap(b.stereo.refA), refB = remap(b.stereo.refB);
      if (refA !== null && refB !== null) {
        sub.bonds[idx].stereo = {refA: refA, refB: refB, same: b.stereo.same};
      }
    }
  }
  return {mol: sub, map: map, rmap: atomIdxs.slice()};
};

/* ------------------------------------------------------------------ *
 * Ring perception
 * ------------------------------------------------------------------ */

/** Shortest path between two atoms, optionally forbidding one bond. */
Mol.prototype.shortestPath = function (from, to, skipBond) {
  if (from === to) return [from];
  var prev = {}, queue = [from], seen = {};
  seen[from] = true;
  while (queue.length) {
    var v = queue.shift();
    var links = this._adj[v];
    for (var k = 0; k < links.length; k++) {
      if (links[k].bond === skipBond) continue;
      var w = links[k].atom;
      if (seen[w]) continue;
      seen[w] = true;
      prev[w] = v;
      if (w === to) {
        var path = [to], cur = to;
        while (cur !== from) { cur = prev[cur]; path.push(cur); }
        return path.reverse();
      }
      queue.push(w);
    }
  }
  return null;
};

/**
 * Smallest set of smallest rings. Each ring is a cycle-ordered atom list.
 */
Mol.prototype.rings = function () {
  if (this._cache.rings) return this._cache.rings;
  var nRings = this.bonds.length - this.atoms.length + this.components().length;
  var out = [];
  if (nRings > 0) {
    var cands = [], i;
    for (i = 0; i < this.bonds.length; i++) {
      var bond = this.bonds[i];
      var path = this.shortestPath(bond.a, bond.b, bond);
      if (path) cands.push(path);
    }
    cands.sort(function (x, y) { return x.length - y.length; });
    var basis = [];
    for (i = 0; i < cands.length && out.length < nRings; i++) {
      var vec = this._ringBondVector(cands[i]);
      var reduced = vec.slice(), independent = false, j;
      for (j = 0; j < basis.length; j++) {
        var pivot = basis[j].pivot;
        if (reduced[pivot]) for (var k = 0; k < reduced.length; k++) reduced[k] ^= basis[j].vec[k];
      }
      for (j = 0; j < reduced.length; j++) if (reduced[j]) { independent = true; break; }
      if (independent) {
        var p = reduced.indexOf(1);
        basis.push({vec: reduced, pivot: p});
        out.push(cands[i]);
      }
    }
  }
  this._cache.rings = out;
  return out;
};

Mol.prototype._ringBondVector = function (ring) {
  var vec = new Array(this.bonds.length);
  for (var i = 0; i < this.bonds.length; i++) vec[i] = 0;
  for (i = 0; i < ring.length; i++) {
    var bond = this.bondBetween(ring[i], ring[(i + 1) % ring.length]);
    if (bond) vec[bond.idx] = 1;
  }
  return vec;
};

Mol.prototype.isInRing = function (i) {
  var rs = this.rings();
  for (var k = 0; k < rs.length; k++) if (rs[k].indexOf(i) >= 0) return true;
  return false;
};

Mol.prototype.isRingBond = function (bond) {
  var rs = this.rings();
  for (var k = 0; k < rs.length; k++) {
    var ring = rs[k];
    for (var i = 0; i < ring.length; i++) {
      if (this.bondBetween(ring[i], ring[(i + 1) % ring.length]) === bond) return true;
    }
  }
  return false;
};

/** Number of SSSR rings each atom belongs to. */
Mol.prototype.ringMembership = function (i) {
  var rs = this.rings(), n = 0;
  for (var k = 0; k < rs.length; k++) if (rs[k].indexOf(i) >= 0) n++;
  return n;
};

/**
 * Fused ring systems: rings sharing at least one bond are grouped together.
 * Returns [{atoms: [...], rings: [ring, ...]}]
 */
Mol.prototype.ringSystems = function () {
  if (this._cache.ringSystems) return this._cache.ringSystems;
  var rings = this.rings(), groups = [], i, k;
  for (i = 0; i < rings.length; i++) {
    var merged = {atoms: rings[i].slice(), rings: [rings[i]]};
    var rest = [];
    for (k = 0; k < groups.length; k++) {
      if (shareTwoAtoms(groups[k].atoms, merged.atoms)) {
        merged.atoms = union(merged.atoms, groups[k].atoms);
        merged.rings = merged.rings.concat(groups[k].rings);
      } else rest.push(groups[k]);
    }
    rest.push(merged);
    groups = rest;
  }
  groups.forEach(function (g) { g.atoms.sort(function (a, b) { return a - b; }); });
  this._cache.ringSystems = groups;
  return groups;
};

Mol.prototype.ringSystemOf = function (i) {
  var systems = this.ringSystems();
  for (var k = 0; k < systems.length; k++) if (systems[k].atoms.indexOf(i) >= 0) return systems[k];
  return null;
};

function shareTwoAtoms(a, b) {
  var n = 0;
  for (var i = 0; i < a.length; i++) if (b.indexOf(a[i]) >= 0) n++;
  return n >= 2;
}

function union(a, b) {
  var out = a.slice();
  for (var i = 0; i < b.length; i++) if (out.indexOf(b[i]) < 0) out.push(b[i]);
  return out;
}

/* ------------------------------------------------------------------ *
 * Aromaticity
 * ------------------------------------------------------------------ */

/**
 * Hueckel-style perception over SSSR rings: every ring atom must be able to
 * contribute a p orbital, and the ring must hold 4n+2 pi electrons.
 * Sets atom.aromatic / bond.aromatic and returns the aromatic rings.
 */
Mol.prototype.perceiveAromaticity = function () {
  if (this._cache.aromatic) return this._cache.aromatic;
  var self = this, rings = this.rings(), aromaticRings = [];
  var i;
  for (i = 0; i < this.atoms.length; i++) this.atoms[i].aromatic = false;
  for (i = 0; i < this.bonds.length; i++) this.bonds[i].aromatic = false;

  rings.forEach(function (ring) {
    if (ring.length < 5 || ring.length > 7) return;
    var pi = 0, ok = true;
    ring.forEach(function (idx) {
      if (!ok) return;
      var e = contribution(self, idx, ring);
      if (e === null) ok = false; else pi += e;
    });
    if (!ok) return;
    if ((pi - 2) % 4 !== 0 || pi <= 0) return;
    aromaticRings.push(ring);
  });

  aromaticRings.forEach(function (ring) {
    ring.forEach(function (idx) { self.atoms[idx].aromatic = true; });
    for (var k = 0; k < ring.length; k++) {
      var bond = self.bondBetween(ring[k], ring[(k + 1) % ring.length]);
      if (bond) bond.aromatic = true;
    }
  });
  this._cache.aromatic = aromaticRings;
  return aromaticRings;
};

/** pi electrons an atom brings to a ring, or null if it cannot be aromatic. */
function contribution(mol, idx, ring) {
  var a = mol.atoms[idx];
  var el = a.element;
  if (el !== 'C' && el !== 'N' && el !== 'O' && el !== 'S' && el !== 'Se' && el !== 'P' && el !== 'B') return null;
  var links = mol.links(idx), inRingDouble = false, exoDouble = null;
  for (var k = 0; k < links.length; k++) {
    var l = links[k];
    if (l.bond.order === 3) return null;
    if (l.bond.order === 2) {
      if (ring.indexOf(l.atom) >= 0) inRingDouble = true;
      else exoDouble = mol.atoms[l.atom];
    }
  }
  if (inRingDouble) return 1;
  if (exoDouble) {
    // A double bond shared with a fused neighbour ring still donates one electron.
    var sys = mol.ringSystemOf(idx);
    if (sys && sys.atoms.indexOf(exoDouble.idx) >= 0) return 1;
    // C=O / C=S in the ring plane keeps sp2 but donates nothing.
    var el2 = exoDouble.element;
    return (el2 === 'O' || el2 === 'S' || el2 === 'N') ? 0 : null;
  }
  if (el === 'C') {
    if (a.charge === -1) return 2;
    if (a.charge === 1) return 0;
    return null;                       // sp3 carbon breaks the ring
  }
  if (el === 'B') return 0;
  if (a.charge === 1 && (el === 'N' || el === 'P')) return null;
  return 2;                            // N, O, S lone pair
}

Mol.prototype.isAromaticAtom = function (i) {
  this.perceiveAromaticity();
  return !!this.atoms[i].aromatic;
};

/* ------------------------------------------------------------------ *
 * Morgan / Weisfeiler-Lehman colouring - used for structure matching
 * ------------------------------------------------------------------ */

function hashString(s) {
  var h = 2166136261;
  for (var i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i);
    h = (h * 16777619) >>> 0;
  }
  return h >>> 0;
}

/** Number of double / triple bonds at an atom (Kekule-invariant per atom). */
Mol.prototype.piCounts = function (i) {
  var d = 0, t = 0, links = this.links(i);
  for (var k = 0; k < links.length; k++) {
    if (links[k].bond.order === 2) d++;
    else if (links[k].bond.order === 3) t++;
  }
  return {d: d, t: t};
};

/**
 * Atom label used for structure comparison. Deliberately built from
 * per-atom pi counts rather than individual bond orders, so that two
 * Kekule drawings of the same aromatic ring compare equal.
 */
Mol.prototype.atomLabel = function (i) {
  var a = this.atoms[i], pi = this.piCounts(i);
  return [a.element, a.charge, this.totalH(i), pi.d, pi.t].join('|');
};

/**
 * Refined colours, invariant under atom renumbering (Weisfeiler-Lehman).
 * Bonds are compared by adjacency only - the unsaturation lives in the
 * atom labels - which keeps the colouring resonance-independent.
 */
Mol.prototype.colors = function () {
  if (this._cache.colors) return this._cache.colors;
  var n = this.atoms.length, colors = new Array(n), i;
  for (i = 0; i < n; i++) colors[i] = hashString(this.atomLabel(i) + '|' + this.degree(i));
  for (var round = 0; round < Math.min(n, 12); round++) {
    var next = new Array(n);
    for (i = 0; i < n; i++) {
      var parts = this.neighbors(i).map(function (j) { return String(colors[j]); });
      parts.sort();
      next[i] = hashString(colors[i] + '<' + parts.join(',') + '>');
    }
    colors = next;
  }
  this._cache.colors = colors;
  return colors;
};

/** Order independent fingerprint of the whole structure. */
Mol.prototype.structureKey = function () {
  if (this._cache.key) return this._cache.key;
  var colors = this.colors().slice().sort(function (a, b) { return a - b; });
  this._cache.key = this.formula() + '/' + hashString(colors.join(','));
  return this._cache.key;
};

ONG.ELEMENTS = ELEMENTS;
ONG.HALOGENS = HALOGENS;
ONG.Mol = Mol;
ONG.hashString = hashString;

if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
