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
  var total = 0;

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
      if (aromatic && bondIdx >= 0) {
        aromaticBonds.push(mol.bonds[bondIdx]);
        if (aromaticAtoms.indexOf(ends[0]) < 0) aromaticAtoms.push(ends[0]);
        if (aromaticAtoms.indexOf(ends[1]) < 0) aromaticAtoms.push(ends[1]);
      }
    }
  });

  function indexOfNode(node) {
    for (var i = 0; i < map.length; i++) if (map[i][0] === node) return map[i][1];
    return -1;
  }

  if (aromaticAtoms.length && ONG.kekulize) {
    ONG.kekulize(mol, aromaticAtoms);
  }

  if (!mol.atoms.length) return {mol: null, warnings: warnings, drawn: total};
  mol._cache = {};
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
