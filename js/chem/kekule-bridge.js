/**
 * kekule-bridge.js - turn what the user drew in the Kekule editor into the
 * molecule graph the naming engine works on.
 *
 * Every call into Kekule is feature-detected: different builds of the
 * library expose slightly different helpers, and a missing one should
 * degrade to a readable message rather than a stack trace.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

/** Molecules currently on the editor's canvas. */
function moleculesOf(composer, Kekule) {
  var out = [];
  try {
    if (composer.exportObjs && Kekule && Kekule.Molecule) {
      out = composer.exportObjs(Kekule.Molecule) || [];
    }
  } catch (e) { out = []; }
  if (out.length) return out;

  // Fall back to walking the edited chem space.
  var root = composer.getChemObj && composer.getChemObj();
  if (!root) return [];
  var stack = [root], seen = [];
  while (stack.length) {
    var obj = stack.pop();
    if (!obj || seen.indexOf(obj) >= 0) continue;
    seen.push(obj);
    if (obj.getNodeCount && obj.getNodeCount() > 0) { out.push(obj); continue; }
    if (obj.getChildCount) {
      for (var i = 0; i < obj.getChildCount(); i++) stack.push(obj.getChildAt(i));
    }
  }
  return out;
}

/** Expand abbreviations (Ph, Et ...) when the build supports it. */
function flatten(fragment) {
  try {
    if (fragment.getFlattenedShadowFragment) {
      var flat = fragment.getFlattenedShadowFragment();
      if (flat && flat.getNodeCount && flat.getNodeCount() > 0) return flat;
    }
  } catch (e) { /* keep the original fragment */ }
  return fragment;
}

/** 2D position of a drawn node, or null. */
function coordOf(node) {
  try {
    var c = node.getAbsCoord2D ? node.getAbsCoord2D() : (node.getCoord2D ? node.getCoord2D() : null);
    if (c && typeof c.x === 'number' && typeof c.y === 'number') return {x: c.x, y: c.y};
  } catch (e) { /* no coordinates */ }
  return null;
}

/**
 * Wedge/hash bonds and 2D positions become the stereo the namer needs.
 *   atom.stereo : which side of the drawing plane each neighbour is on
 *   bond.stereo : which side of a double bond its substituents sit on
 */
function applyDrawnStereo(mol, coords, wedges) {
  var i;

  /* --- tetrahedral centres ----------------------------------------- */
  for (i = 0; i < mol.atoms.length; i++) {
    var neighbors = mol.neighbors(i);
    var implicit = mol.implicitH(i);
    if (neighbors.length + implicit !== 4 || neighbors.length < 3) continue;
    if (mol.bondOrderSum(i) !== mol.degree(i)) continue;          // must be saturated
    if (!coords[i]) continue;
    var hasWedge = false, vectors = [], order = [], ok = true;

    neighbors.forEach(function (j) {
      if (!ok) return;
      if (!coords[j]) { ok = false; return; }
      var z = 0, w = wedges[bondKey(i, j)];
      if (w) {
        // Only a wedge whose narrow end sits on this atom says anything.
        if (w.narrow === i) { z = w.up ? 0.6 : -0.6; hasWedge = true; }
      }
      order.push(j);
      vectors.push({x: coords[j].x - coords[i].x, y: coords[j].y - coords[i].y, z: z});
    });
    if (!ok || !hasWedge) continue;

    if (order.length === 3) {
      // The implicit hydrogen sits opposite the three drawn bonds.
      var sum = vectors.reduce(function (acc, v) {
        return {x: acc.x + v.x, y: acc.y + v.y, z: acc.z + v.z};
      }, {x: 0, y: 0, z: 0});
      vectors.push({x: -sum.x, y: -sum.y, z: -sum.z || -0.6});
      order.push(-1);
    }
    var volume = signedVolume(vectors[0], vectors[1], vectors[2], vectors[3]);
    if (Math.abs(volume) < 1e-6) continue;                        // degenerate drawing
    mol.atoms[i].stereo = {order: order, clockwise: volume < 0};
  }

  /* --- double bond geometry ---------------------------------------- */
  mol.bonds.forEach(function (bond) {
    if (bond.order !== 2 || bond.stereoUnspecified) return;
    var ca = coords[bond.a], cb = coords[bond.b];
    if (!ca || !cb) return;
    var axis = {x: cb.x - ca.x, y: cb.y - ca.y};
    var refA = pickSide(mol, bond.a, bond.b, ca, axis);
    var refB = pickSide(mol, bond.b, bond.a, cb, axis);
    if (!refA || !refB) return;
    bond.stereo = {refA: refA.atom, refB: refB.atom, same: refA.side === refB.side};
  });

  function pickSide(mol, end, other, origin, axis) {
    var subs = mol.neighbors(end).filter(function (j) { return j !== other; });
    for (var k = 0; k < subs.length; k++) {
      var c = coords[subs[k]];
      if (!c) continue;
      var cross = axis.x * (c.y - origin.y) - axis.y * (c.x - origin.x);
      if (Math.abs(cross) > 1e-6) return {atom: subs[k], side: cross > 0 ? 1 : -1};
    }
    return null;
  }
}

function signedVolume(v1, v2, v3, v4) {
  var a = {x: v1.x - v4.x, y: v1.y - v4.y, z: v1.z - v4.z};
  var b = {x: v2.x - v4.x, y: v2.y - v4.y, z: v2.z - v4.z};
  var c = {x: v3.x - v4.x, y: v3.y - v4.y, z: v3.z - v4.z};
  var cross = {
    x: b.y * c.z - b.z * c.y,
    y: b.z * c.x - b.x * c.z,
    z: b.x * c.y - b.y * c.x
  };
  return a.x * cross.x + a.y * cross.y + a.z * cross.z;
}

function bondKey(a, b) { return a < b ? a + '-' + b : b + '-' + a; }

function symbolOf(node) {
  try {
    if (node.getSymbol) {
      var s = node.getSymbol();
      if (typeof s === 'string' && s) return s;
    }
  } catch (e) { /* fall through */ }
  try {
    if (node.getIsotope && node.getIsotope() && node.getIsotope().getSymbol) {
      return node.getIsotope().getSymbol();
    }
  } catch (e) { /* fall through */ }
  return null;
}

/**
 * Build a Mol from the editor contents.
 * Returns {mol, warnings} - mol is null when nothing usable was drawn.
 */
function fromComposer(composer, Kekule) {
  var fragments = moleculesOf(composer, Kekule);
  var warnings = [];
  var mol = new ONG.Mol();
  var map = [];           // [kekuleNode, molIndex]
  var aromaticAtoms = [], aromaticBonds = [];
  var coords = {}, wedges = {};
  var total = 0;
  var BondStereo = (Kekule && Kekule.BondStereo) || {};

  fragments.forEach(function (raw) {
    var frag = flatten(raw);
    if (!frag.getNodeCount) return;
    var n = frag.getNodeCount();
    for (var i = 0; i < n; i++) {
      var node = frag.getNodeAt(i);
      total++;
      if (node.getNodeCount && node.getNodeCount() > 0) {
        warnings.push('구조식 안의 축약기(subgroup)를 펼치지 못했습니다. 전체 구조를 원자 단위로 그려 주세요.');
        continue;
      }
      var symbol = symbolOf(node);
      if (!symbol || !ONG.ELEMENTS[symbol]) {
        warnings.push('이름을 붙일 수 없는 원자가 있습니다' + (symbol ? ' (' + symbol + ')' : '') +
          '. R기나 가변 원자는 지원하지 않습니다.');
        continue;
      }
      var charge = 0;
      try { charge = node.getCharge ? (node.getCharge() || 0) : 0; } catch (e) { charge = 0; }
      var explicitH = null;
      try {
        var h = node.getExplicitHydrogenCount ? node.getExplicitHydrogenCount() : null;
        if (typeof h === 'number' && h > 0) explicitH = h;
      } catch (e) { explicitH = null; }
      var idx = mol.addAtom(symbol, {charge: Math.round(charge), explicitH: explicitH});
      map.push([node, idx]);
      var xy = coordOf(node);
      if (xy) coords[idx] = xy;
    }

    var m = frag.getConnectorCount ? frag.getConnectorCount() : 0;
    for (var k = 0; k < m; k++) {
      var conn = frag.getConnectorAt(k);
      var objs = [];
      try { objs = conn.getConnectedObjs ? conn.getConnectedObjs() : []; } catch (e) { objs = []; }
      var ends = objs.map(indexOfNode).filter(function (x) { return x >= 0; });
      if (ends.length !== 2) continue;
      var order = 1;
      try { order = conn.getBondOrder ? conn.getBondOrder() : 1; } catch (e) { order = 1; }
      var aromatic = false;
      if (!order || order < 1) order = 1;
      if (order > 3) { aromatic = true; order = 1; }   // Kekule marks aromatic bonds as order 10
      var bondIdx = mol.addBond(ends[0], ends[1], order);
      if (bondIdx >= 0) recordWedge(conn, ends, mol.bonds[bondIdx]);
      if (aromatic && bondIdx >= 0) {
        aromaticBonds.push(mol.bonds[bondIdx]);
        if (aromaticAtoms.indexOf(ends[0]) < 0) aromaticAtoms.push(ends[0]);
        if (aromaticAtoms.indexOf(ends[1]) < 0) aromaticAtoms.push(ends[1]);
      }
    }
  });

  /** Wedge ("up") and hash ("down") bonds, with their narrow end. */
  function recordWedge(conn, ends, bond) {
    var stereo = null;
    try { stereo = conn.getStereo ? conn.getStereo() : null; } catch (e) { stereo = null; }
    if (stereo === null || stereo === undefined) return;
    var UP = BondStereo.UP === undefined ? 1 : BondStereo.UP;
    var UP_INV = BondStereo.UP_INVERTED === undefined ? 2 : BondStereo.UP_INVERTED;
    var DOWN = BondStereo.DOWN === undefined ? 3 : BondStereo.DOWN;
    var DOWN_INV = BondStereo.DOWN_INVERTED === undefined ? 4 : BondStereo.DOWN_INVERTED;
    var EITHER = [BondStereo.UP_OR_DOWN, BondStereo.UP_OR_DOWN_INVERTED, BondStereo.CIS_OR_TRANS];
    if (EITHER.indexOf(stereo) >= 0) { bond.stereoUnspecified = true; return; }
    var up = (stereo === UP || stereo === UP_INV);
    var down = (stereo === DOWN || stereo === DOWN_INV);
    if (!up && !down) return;
    var inverted = (stereo === UP_INV || stereo === DOWN_INV);
    var narrow = inverted ? ends[1] : ends[0];
    var wide = inverted ? ends[0] : ends[1];
    wedges[bondKey(ends[0], ends[1])] = {narrow: narrow, wide: wide, up: up};
  }

  function indexOfNode(node) {
    for (var i = 0; i < map.length; i++) if (map[i][0] === node) return map[i][1];
    return -1;
  }

  if (aromaticAtoms.length && ONG.kekulize) {
    ONG.kekulize(mol, aromaticAtoms);
  }

  if (!mol.atoms.length) return {mol: null, warnings: warnings, drawn: total};
  mol._cache = {};
  try {
    applyDrawnStereo(mol, coords, wedges);
  } catch (e) {
    warnings.push('입체 배치(쐐기 결합)를 읽는 중 문제가 있었습니다: ' + e.message);
  }
  return {mol: mol, warnings: unique(warnings), drawn: total};
}

function unique(list) {
  var out = [];
  list.forEach(function (x) { if (out.indexOf(x) < 0) out.push(x); });
  return out;
}

ONG.fromComposer = fromComposer;
if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
