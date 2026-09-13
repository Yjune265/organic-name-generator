/**
 * variants.js - the other names the same structure legitimately goes by.
 *
 * Students meet a compound under several names: the preferred IUPAC name,
 * the pre-1993 locant style (2-butanol), retained substituent prefixes
 * (isopropyl, tert-butyl, neopentyl) and functional class names
 * (isopropyl alcohol, ethyl methyl ether). All of them are generated here.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

/** Retained names for whole molecules, keyed by the systematic name. */
var RETAINED_PARENTS = {
  '2-methylpropane': 'isobutane',
  '2-methylbutane': 'isopentane',
  '2,2-dimethylpropane': 'neopentane',
  '2,2,4-trimethylpentane': 'isooctane',
  'ethene': 'ethylene',
  'prop-1-ene': 'propylene',
  'ethyne': 'acetylene',
  '2-methylbuta-1,3-diene': 'isoprene',
  'methylbenzene': 'toluene',
  'ethenylbenzene': 'styrene',
  'methoxybenzene': 'anisole',
  'ethanoic acid': 'acetic acid',
  'methanoic acid': 'formic acid',
  'propanoic acid': 'propionic acid',
  'butanoic acid': 'butyric acid',
  'ethanamide': 'acetamide',
  'ethanal': 'acetaldehyde',
  'methanal': 'formaldehyde',
  'propan-2-one': 'acetone',
  'ethanenitrile': 'acetonitrile',
  'trichloromethane': 'chloroform',
  'ethanoic anhydride': 'acetic anhydride'
};

/** Word-level retained forms that may appear inside a larger name. */
var RETAINED_FRAGMENTS = [
  [/\bethanoic acid\b/g, 'acetic acid'],
  [/\bethanamide\b/g, 'acetamide'],
  [/\bethanoate\b/g, 'acetate'],
  [/\bethanoyl\b/g, 'acetyl'],
  [/\bmethanoic acid\b/g, 'formic acid'],
  [/\bmethanoate\b/g, 'formate'],
  [/\bmethanamide\b/g, 'formamide']
];

/**
 * Pre-1993 locant style: the parent's locant moves in front of its stem.
 *   butan-2-ol          -> 2-butanol
 *   hex-5-en-2-ol       -> 5-hexen-2-ol
 *   cyclohex-2-en-1-one -> 2-cyclohexen-1-one
 * Only the parent is rewritten; substituent prefixes keep their locants.
 */
function classicLocants(parts) {
  if (!parts || !parts.parent) return null;
  var m = /^([a-z]+)-(\d[\d,]*)-/.exec(parts.parent);
  if (!m) return null;
  var parent = m[2] + '-' + m[1] + parts.parent.slice(m[0].length);
  var name = parts.prefix ? parts.prefix + '-' + parent : parent;
  if (parts.esterAlkyl) name = parts.esterAlkyl + ' ' + name;
  return name + (parts.tail || '');
}

/** Functional class names: "isopropyl alcohol", "ethyl methyl ether". */
var HALIDE_WORD = {F: 'fluoride', Cl: 'chloride', Br: 'bromide', I: 'iodide'};

function functionalClassName(mol) {
  var groups = ONG.perceiveGroups(mol);
  var all = groups.all;
  if (all.length !== 1) return null;
  var g = all[0];
  var style = {retained: true};

  // Only simple groups get a functional class name; anything needing a
  // locant of its own (but-2-en-1-yl ...) reads better systematically.
  function substituentOf(start, blocked) {
    var atoms = ONG.naming.branchAtoms(mol, start, blocked, ONG.naming.allSet(mol));
    var c = new ONG.naming.Ctx(mol, style);
    var unit = ONG.naming.nameUnit(c, atoms, start);
    if (!unit || c.warnings.length) return null;
    return /[\d(]/.test(unit.name) ? null : unit.name;
  }

  function twoSided(centre, skipElements) {
    var sides = mol.heavyNeighbors(centre).filter(function (j) {
      return !skipElements || skipElements.indexOf(mol.atoms[j].element) < 0;
    });
    if (sides.length !== 2) return null;
    var blocked = {};
    blocked[centre] = true;
    var a = substituentOf(sides[0], blocked), b = substituentOf(sides[1], blocked);
    if (!a || !b) return null;
    return a === b ? ['di' + a] : [a, b].sort();
  }

  var blocked = {}, r, pair;
  switch (g.type) {
    case 'ether':
      pair = twoSided(g.oxygen);
      return pair ? pair.join(' ') + ' ether' : null;
    case 'sulfide':
      pair = twoSided(g.sulfur);
      return pair ? pair.join(' ') + ' sulfide' : null;
    case 'ketone':
      pair = twoSided(g.carbon, ['O']);
      return pair ? pair.join(' ') + ' ketone' : null;
    case 'alcohol':
      blocked[g.oxygen] = true;
      r = substituentOf(g.carrier, blocked);
      return r ? r + ' alcohol' : null;
    case 'thiol':
      blocked[g.sulfur] = true;
      r = substituentOf(g.carrier, blocked);
      return r ? r + ' mercaptan' : null;
    case 'halide': {
      var carrier = mol.heavyNeighbors(g.atom)[0];
      if (carrier === undefined) return null;
      blocked[g.atom] = true;
      r = substituentOf(carrier, blocked);
      return r ? r + ' ' + HALIDE_WORD[mol.atoms[g.atom].element] : null;
    }
    case 'amine': {
      var names = [], ok = true;
      blocked[g.nitrogen] = true;
      mol.heavyNeighbors(g.nitrogen).forEach(function (j) {
        var n = substituentOf(j, blocked);
        if (!n) ok = false; else names.push(n);
      });
      if (!ok || !names.length) return null;
      var counts = {};
      names.forEach(function (n) { counts[n] = (counts[n] || 0) + 1; });
      return Object.keys(counts).sort().map(function (n) {
        return ({1: '', 2: 'di', 3: 'tri'}[counts[n]] || '') + n;
      }).join('') + 'amine';
    }
    default:
      return null;
  }
}

/**
 * Every acceptable alternative for one molecule.
 * Returns [{name, kind, note}] with no duplicates and without the primary.
 */
function alternatives(mol, primary, primaryParts) {
  var out = [], seen = {};
  if (primary) seen[primary] = true;

  function add(name, kind, note) {
    if (!name || seen[name]) return;
    seen[name] = true;
    out.push({name: name, kind: kind, note: note});
  }

  // 1. retained substituent prefixes (isopropyl, tert-butyl, neopentyl ...)
  var retained = null, retainedParts = null;
  try {
    var ctx = new ONG.naming.Ctx(mol, {retained: true});
    var unit = ONG.naming.nameUnit(ctx, mol.atoms.map(function (a, i) { return i; }), null);
    if (unit && !ctx.warnings.length) { retained = unit.name; retainedParts = unit.parts; }
  } catch (e) { retained = null; }
  if (retained && retained !== primary) {
    add(retained, 'retained', '관용 접두사(iso-/sec-/tert-/neo-)를 사용한 이름');
  }

  // 2. retained name for the whole molecule (isobutane, acetone, toluene ...)
  [primary, retained].forEach(function (base) {
    if (base && RETAINED_PARENTS[base]) add(RETAINED_PARENTS[base], 'trivial', '관용명');
  });

  // 3. retained acyl/acid stems inside a longer name (acetic acid, acetyl ...)
  [primary, retained].forEach(function (base) {
    if (!base) return;
    var swapped = base;
    RETAINED_FRAGMENTS.forEach(function (pair) { swapped = swapped.replace(pair[0], pair[1]); });
    if (swapped !== base) add(swapped, 'retained', '관용 어간(acetic/acetyl 등)을 사용한 이름');
  });

  // 4. pre-1993 locant placement (2-butanol, 4-penten-2-ol)
  [primaryParts, retainedParts].forEach(function (parts) {
    add(classicLocants(parts), 'classic', 'locant을 앞에 쓰는 예전(1993년 이전) 방식');
  });

  // 5. functional class name (isopropyl alcohol, ethyl methyl ether)
  add(functionalClassName(mol), 'functional', '작용기 분류명(functional class name)');

  return out;
}

ONG.variants = {
  alternatives: alternatives,
  classicLocants: classicLocants,
  functionalClassName: functionalClassName,
  RETAINED_PARENTS: RETAINED_PARENTS
};

if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
