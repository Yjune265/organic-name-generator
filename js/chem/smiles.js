/**
 * smiles.js - a compact SMILES reader.
 *
 * Supports: organic subset, bracket atoms, charges, explicit H counts,
 * bond symbols, ring-closure bonds, branches, dot-disconnected parts and
 * lowercase aromatic atoms (kekulized on the way in).
 * Stereo descriptors (@ / @@ / cis-trans slashes) are parsed and ignored -
 * this project names constitution, not configuration.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};
var Mol = ONG.Mol;

var ORGANIC = ['Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I'];
var AROMATIC_LOWER = {c: 'C', n: 'N', o: 'O', s: 'S', p: 'P', b: 'B', se: 'Se'};

function parseSmiles(text) {
  if (!text || !text.trim()) throw new Error('빈 SMILES 입니다.');
  var s = text.trim();
  var mol = new Mol();
  var stack = [];         // branch stack of previous atom indices
  var prev = null;        // previous atom index
  var pendingBond = null; // bond order waiting for the next atom
  var pendingAromaticBond = false;
  var pendingDirection = null;  // '/' or '\\' waiting for the next atom
  var ringBonds = {};     // ring-closure digit -> {atom, order, aromatic}
  var stereoOrder = {};   // atom -> neighbours in written order (-1 = implicit H)
  var chirality = {};     // atom -> '@' | '@@'
  var directional = [];   // {bond, from, to, up} for / and \\
  var aromaticAtoms = [];
  var aromaticBonds = [];
  var i = 0;

  function readBracketAtom() {
    var end = s.indexOf(']', i);
    if (end < 0) throw new Error('SMILES 대괄호가 닫히지 않았습니다.');
    var body = s.slice(i + 1, end);
    i = end + 1;
    var m = /^(\d+)?([A-Za-z][a-z]?)(@{1,2})?(H\d*)?((?:[+-]\d*|\++|-+))?(?::\d+)?$/.exec(body);
    if (!m) throw new Error('해석할 수 없는 원자 표기: [' + body + ']');
    var symbol = m[2], aromatic = false;
    if (AROMATIC_LOWER[symbol]) { aromatic = true; symbol = AROMATIC_LOWER[symbol]; }
    else if (symbol[0] === symbol[0].toLowerCase() && symbol !== 'H') {
      symbol = symbol[0].toUpperCase() + symbol.slice(1);
    }
    var explicitH = 0;
    if (m[4]) explicitH = m[4].length > 1 ? parseInt(m[4].slice(1), 10) : 1;
    var charge = 0;
    if (m[5]) {
      var c = m[5];
      if (/^[+-]\d+$/.test(c)) charge = parseInt(c, 10);
      else charge = (c[0] === '+' ? 1 : -1) * c.length;
    }
    var idx = mol.addAtom(symbol, {
      charge: charge, explicitH: explicitH,
      isotope: m[1] ? parseInt(m[1], 10) : null
    });
    if (aromatic) aromaticAtoms.push(idx);
    if (m[3]) {
      chirality[idx] = m[3];
      // The implicit hydrogen of a chiral bracket atom counts at the position
      // right after the preceding atom (or first, when there is none).
      stereoOrder[idx] = [];
      if (prev !== null) stereoOrder[idx].push(prev);
      if (explicitH > 0) stereoOrder[idx].push(-1);
    }
    return idx;
  }

  function readOrganicAtom() {
    for (var k = 0; k < ORGANIC.length; k++) {
      var sym = ORGANIC[k];
      if (s.substr(i, sym.length) === sym) {
        i += sym.length;
        return {idx: mol.addAtom(sym), aromatic: false};
      }
    }
    var two = s.substr(i, 2).toLowerCase();
    if (AROMATIC_LOWER[two]) { i += 2; return {idx: mol.addAtom(AROMATIC_LOWER[two]), aromatic: true}; }
    var one = s[i];
    if (AROMATIC_LOWER[one]) { i += 1; return {idx: mol.addAtom(AROMATIC_LOWER[one]), aromatic: true}; }
    return null;
  }

  function noteNeighbor(atom, neighbor) {
    if (stereoOrder[atom]) stereoOrder[atom].push(neighbor);
  }

  function attach(idx, aromatic) {
    if (prev !== null) {
      noteNeighbor(prev, idx);
      if (stereoOrder[idx] && stereoOrder[idx].indexOf(prev) < 0) stereoOrder[idx].unshift(prev);
      var order = pendingBond;
      var isAromaticBond = pendingAromaticBond ||
        (order === null && aromatic && mol.atoms[prev].__aromatic);
      if (order === null) order = 1;
      var bond = mol.bonds[mol.addBond(prev, idx, order)];
      if (isAromaticBond && bond) aromaticBonds.push(bond);
      if (pendingDirection && bond) {
        directional.push({bond: bond, from: prev, to: idx, up: pendingDirection === '/'});
      }
      pendingDirection = null;
    }
    mol.atoms[idx].__aromatic = !!aromatic;
    prev = idx;
    pendingBond = null;
    pendingAromaticBond = false;
  }

  while (i < s.length) {
    var ch = s[i];

    if (ch === '[') { var bidx = readBracketAtom(); attach(bidx, aromaticAtoms.indexOf(bidx) >= 0); continue; }
    if (ch === '(') { stack.push(prev); i++; continue; }
    if (ch === ')') {
      if (!stack.length) throw new Error('괄호 짝이 맞지 않습니다.');
      prev = stack.pop(); i++; continue;
    }
    if (ch === '.') { prev = null; pendingBond = null; i++; continue; }
    if (ch === '-') { pendingBond = 1; i++; continue; }
    if (ch === '=') { pendingBond = 2; i++; continue; }
    if (ch === '#') { pendingBond = 3; i++; continue; }
    if (ch === ':') { pendingBond = 1; pendingAromaticBond = true; i++; continue; }
    if (ch === '/' || ch === '\\') { pendingDirection = ch; i++; continue; }
    if (ch === '%' || /\d/.test(ch)) {
      var label;
      if (ch === '%') { label = s.substr(i + 1, 2); i += 3; }
      else { label = ch; i += 1; }
      if (prev === null) throw new Error('고리 닫기 번호의 위치가 잘못되었습니다.');
      if (stereoOrder[prev]) stereoOrder[prev].push({ring: label});
      if (ringBonds[label]) {
        var open = ringBonds[label];
        delete ringBonds[label];
        resolveRingPlace(open.atom, label, prev);
        resolveRingPlace(prev, label, open.atom);
        var order = pendingBond !== null ? pendingBond : (open.order !== null ? open.order : 1);
        var rb = mol.bonds[mol.addBond(open.atom, prev, order)];
        var bothAromatic = mol.atoms[open.atom].__aromatic && mol.atoms[prev].__aromatic;
        if (rb && (open.aromatic || pendingAromaticBond || (pendingBond === null && bothAromatic))) {
          aromaticBonds.push(rb);
        }
      } else {
        ringBonds[label] = {atom: prev, order: pendingBond, aromatic: pendingAromaticBond};
      }
      pendingBond = null;
      pendingAromaticBond = false;
      continue;
    }
    if (/[A-Za-z]/.test(ch)) {
      var atom = readOrganicAtom();
      if (!atom) throw new Error('알 수 없는 문자: "' + ch + '"');
      if (atom.aromatic) aromaticAtoms.push(atom.idx);
      attach(atom.idx, atom.aromatic);
      continue;
    }
    throw new Error('알 수 없는 문자: "' + ch + '"');
  }

  if (Object.keys(ringBonds).length) throw new Error('닫히지 않은 고리 번호가 있습니다.');

  function resolveRingPlace(atom, label, partner) {
    var list = stereoOrder[atom];
    if (!list) return;
    for (var k = 0; k < list.length; k++) {
      if (list[k] && list[k].ring === label) { list[k] = partner; return; }
    }
  }

  function applyStereo() {
    Object.keys(chirality).forEach(function (key) {
      var idx = parseInt(key, 10);
      var order = (stereoOrder[idx] || []).filter(function (x) { return typeof x === 'number'; });
      var expected = mol.neighbors(idx).length + mol.implicitH(idx);
      if (order.length !== 4 || expected !== 4) return;   // only tetrahedral centres
      mol.atoms[idx].stereo = {order: order, clockwise: chirality[key] === '@@'};
    });

    // "/" and "\" describe where a substituent sits relative to its double bond.
    mol.bonds.forEach(function (bond) {
      if (bond.order !== 2) return;
      var a = sideOf(bond.a), b = sideOf(bond.b);
      if (!a || !b) return;
      bond.stereo = {refA: a.atom, refB: b.atom, same: a.label === b.label};
    });

    function sideOf(end) {
      for (var k = 0; k < directional.length; k++) {
        var d = directional[k];
        if (d.from === end) return {atom: d.to, label: d.up ? 'U' : 'D'};
        if (d.to === end) return {atom: d.from, label: d.up ? 'D' : 'U'};
      }
      return null;
    }
  }

  applyStereo();

  if (aromaticAtoms.length) kekulize(mol, aromaticAtoms);
  mol.atoms.forEach(function (a) { delete a.__aromatic; });
  mol._cache = {};
  return mol;
}

/**
 * Turn a lowercase-aromatic input into an explicit Kekule structure by
 * finding a perfect matching over the aromatic atoms that still need a
 * double bond. Pyrrole-type N/O/S and carbonyl carbons are left alone.
 */
function kekulize(mol, aromaticAtoms) {
  var need = [];
  aromaticAtoms.forEach(function (idx) {
    var a = mol.atoms[idx];
    var sum = mol.bondOrderSum(idx) + (a.explicitH || 0);
    var hasDouble = mol.links(idx).some(function (l) { return l.bond.order >= 2; });
    if (hasDouble) return;
    if (a.element === 'C') {
      // aromatic carbon with three heavy neighbours and no H is already saturated
      var open = 4 - sum - (a.explicitH === null ? mol.implicitH(idx) : 0);
      if (a.explicitH !== null) { if (4 - sum <= 0) return; }
      else if (mol.degree(idx) >= 4) return;
      need.push(idx);
    } else if (a.element === 'N' || a.element === 'P') {
      var maxVal = 3 + (a.charge > 0 ? 1 : 0);
      if (sum >= maxVal) return;            // pyrrole-type N already full
      if (a.explicitH) return;
      need.push(idx);
    }
    // O, S in aromatic rings never take a double bond in the ring
  });

  var pairs = {};
  var ok = match(0);
  function match(k) {
    while (k < need.length && pairs[need[k]] !== undefined) k++;
    if (k >= need.length) return true;
    var idx = need[k];
    var links = mol.links(idx);
    for (var m = 0; m < links.length; m++) {
      var j = links[m].atom;
      if (pairs[j] !== undefined) continue;
      if (need.indexOf(j) < 0) continue;
      if (links[m].bond.order !== 1) continue;
      pairs[idx] = j; pairs[j] = idx;
      links[m].bond.order = 2;
      if (match(k + 1)) return true;
      links[m].bond.order = 1;
      delete pairs[idx]; delete pairs[j];
    }
    return false;
  }
  if (!ok) {
    // Leave the structure as drawn; aromaticity perception may still cope.
    return false;
  }
  return true;
}

ONG.parseSmiles = parseSmiles;
ONG.kekulize = kekulize;
if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
