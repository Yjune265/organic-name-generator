/**
 * rings.js - recognition and IUPAC numbering of ring systems.
 *
 * Named ring systems are stored as templates whose atoms are listed in
 * IUPAC numbering order. Matching a template onto a drawn ring system by
 * graph isomorphism therefore yields both the name and every numbering the
 * fusion rules allow; the naming code then picks the numbering that gives
 * the lowest locants.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};
var Mol = ONG.Mol;

/* ------------------------------------------------------------------ *
 * Template table
 * ------------------------------------------------------------------ */

// atoms: "El:label ..."   bonds: "label-label" (single), "=" double, "#" triple
function tpl(name, atoms, bonds, opts) {
  var t = Object.assign({name: name, atomSpec: atoms, bondSpec: bonds}, opts || {});
  return t;
}

var TEMPLATES = [
  /* --- benzenoid --- */
  tpl('naphthalene',
      'C:1 C:2 C:3 C:4 C:4a C:5 C:6 C:7 C:8 C:8a',
      '1=2 2-3 3=4 4-4a 4a-5 5=6 6-7 7=8 8-8a 8a-1 4a=8a', {aromatic: true}),
  tpl('anthracene',
      'C:1 C:2 C:3 C:4 C:4a C:5 C:6 C:7 C:8 C:8a C:9 C:9a C:10 C:10a',
      '1=2 2-3 3=4 4-4a 4a-9a 9a-1 4a=10 10-10a 10a=8a 10a-5 8a-9 9=9a 5=6 6-7 7=8 8-8a',
      {aromatic: true}),
  tpl('phenanthrene',
      'C:1 C:2 C:3 C:4 C:4a C:4b C:5 C:6 C:7 C:8 C:8a C:9 C:10 C:10a',
      '1=2 2-3 3=4 4-4a 4a=10a 10a-1 4a-4b 4b=8a 8a-9 9=10 10-10a 4b-5 5=6 6-7 7=8 8-8a',
      {aromatic: true}),
  tpl('1H-indene',
      'C:1 C:2 C:3 C:3a C:4 C:5 C:6 C:7 C:7a',
      '1-2 2=3 3-3a 3a-4 4=5 5-6 6=7 7-7a 7a-1 3a=7a'),
  tpl('2,3-dihydro-1H-indene',
      'C:1 C:2 C:3 C:3a C:4 C:5 C:6 C:7 C:7a',
      '1-2 2-3 3-3a 3a-4 4=5 5-6 6=7 7-7a 7a-1 3a=7a', {trivial: 'indane'}),
  tpl('9H-fluorene',
      'C:1 C:2 C:3 C:4 C:4a C:4b C:5 C:6 C:7 C:8 C:8a C:9 C:9a',
      '1=2 2-3 3=4 4-4a 4a=9a 9a-1 4a-4b 4b-5 5=6 6-7 7=8 8-8a 8a=4b 8a-9 9-9a'),

  /* --- six-membered heteroaromatics --- */
  tpl('pyridine',    'N:1 C:2 C:3 C:4 C:5 C:6', '1=2 2-3 3=4 4-5 5=6 6-1', {aromatic: true}),
  tpl('pyridazine',  'N:1 N:2 C:3 C:4 C:5 C:6', '1=2 2-3 3=4 4-5 5=6 6-1', {aromatic: true}),
  tpl('pyrimidine',  'N:1 C:2 N:3 C:4 C:5 C:6', '1=2 2-3 3=4 4-5 5=6 6-1', {aromatic: true}),
  tpl('pyrazine',    'N:1 C:2 C:3 N:4 C:5 C:6', '1=2 2-3 3=4 4-5 5=6 6-1', {aromatic: true}),
  tpl('1,3,5-triazine', 'N:1 C:2 N:3 C:4 N:5 C:6', '1=2 2-3 3=4 4-5 5=6 6-1', {aromatic: true}),

  /* --- five-membered heteroaromatics --- */
  tpl('1H-pyrrole',  'N:1 C:2 C:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true}),
  tpl('furan',       'O:1 C:2 C:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true}),
  tpl('thiophene',   'S:1 C:2 C:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true}),
  tpl('1H-imidazole','N:1 C:2 N:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true}),
  tpl('1H-pyrazole', 'N:1 N:2 C:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true}),
  tpl('1,3-oxazole', 'O:1 C:2 N:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true, trivial: 'oxazole'}),
  tpl('1,2-oxazole', 'O:1 N:2 C:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true, trivial: 'isoxazole'}),
  tpl('1,3-thiazole','S:1 C:2 N:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true, trivial: 'thiazole'}),
  tpl('1,2-thiazole','S:1 N:2 C:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true, trivial: 'isothiazole'}),
  tpl('1H-1,2,3-triazole', 'N:1 N:2 N:3 C:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true}),
  tpl('1H-1,2,4-triazole', 'N:1 N:2 C:3 N:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true}),
  tpl('1H-tetrazole',      'N:1 N:2 N:3 N:4 C:5', '1-2 2=3 3-4 4=5 5-1', {aromatic: true}),

  /* --- fused heteroaromatics --- */
  tpl('quinoline',   'N:1 C:2 C:3 C:4 C:4a C:5 C:6 C:7 C:8 C:8a',
      '1=2 2-3 3=4 4-4a 4a-5 5=6 6-7 7=8 8-8a 8a-1 4a=8a', {aromatic: true}),
  tpl('isoquinoline','C:1 N:2 C:3 C:4 C:4a C:5 C:6 C:7 C:8 C:8a',
      '1=2 2-3 3=4 4-4a 4a-5 5=6 6-7 7=8 8-8a 8a-1 4a=8a', {aromatic: true}),
  tpl('quinazoline', 'N:1 C:2 N:3 C:4 C:4a C:5 C:6 C:7 C:8 C:8a',
      '1=2 2-3 3=4 4-4a 4a-5 5=6 6-7 7=8 8-8a 8a-1 4a=8a', {aromatic: true}),
  tpl('quinoxaline', 'N:1 C:2 C:3 N:4 C:4a C:5 C:6 C:7 C:8 C:8a',
      '1=2 2-3 3=4 4-4a 4a-5 5=6 6-7 7=8 8-8a 8a-1 4a=8a', {aromatic: true}),
  tpl('1H-indole',   'N:1 C:2 C:3 C:3a C:4 C:5 C:6 C:7 C:7a',
      '1-2 2=3 3-3a 3a-4 4=5 5-6 6=7 7-7a 7a-1 3a=7a', {aromatic: true}),
  tpl('2H-isoindole','C:1 N:2 C:3 C:3a C:4 C:5 C:6 C:7 C:7a',
      '1-2 2-3 3=3a 3a-4 4=5 5-6 6=7 7-7a 7a=1 3a-7a', {aromatic: true}),
  tpl('1-benzofuran','O:1 C:2 C:3 C:3a C:4 C:5 C:6 C:7 C:7a',
      '1-2 2=3 3-3a 3a-4 4=5 5-6 6=7 7-7a 7a-1 3a=7a', {aromatic: true, trivial: 'benzofuran'}),
  tpl('1-benzothiophene','S:1 C:2 C:3 C:3a C:4 C:5 C:6 C:7 C:7a',
      '1-2 2=3 3-3a 3a-4 4=5 5-6 6=7 7-7a 7a-1 3a=7a', {aromatic: true}),
  tpl('1H-benzimidazole','N:1 C:2 N:3 C:3a C:4 C:5 C:6 C:7 C:7a',
      '1-2 2=3 3-3a 3a-4 4=5 5-6 6=7 7-7a 7a-1 3a=7a', {aromatic: true}),
  tpl('9H-purine',   'N:1 C:2 N:3 C:4 C:5 C:6 N:7 C:8 N:9',
      '1=2 2-3 3=4 4-5 5=6 6-1 5-7 7=8 8-9 9-4', {aromatic: true}),
  tpl('7H-purine',   'N:1 C:2 N:3 C:4 C:5 C:6 N:7 C:8 N:9',
      '1=2 2-3 3=4 4-5 5=6 6-1 5-7 7-8 8=9 9-4', {aromatic: true}),
  tpl('9H-carbazole','C:1 C:2 C:3 C:4 C:4a C:4b C:5 C:6 C:7 C:8 C:8a N:9 C:9a',
      '1=2 2-3 3=4 4-4a 4a=9a 9a-1 4a-4b 4b-5 5=6 6-7 7=8 8-8a 8a=4b 8a-9 9-9a',
      {aromatic: true}),

  /* --- saturated / partly saturated rings --- */
  tpl('oxirane',   'O:1 C:2 C:3', '1-2 2-3 3-1'),
  tpl('aziridine', 'N:1 C:2 C:3', '1-2 2-3 3-1'),
  tpl('thiirane',  'S:1 C:2 C:3', '1-2 2-3 3-1'),
  tpl('oxetane',   'O:1 C:2 C:3 C:4', '1-2 2-3 3-4 4-1'),
  tpl('azetidine', 'N:1 C:2 C:3 C:4', '1-2 2-3 3-4 4-1'),
  tpl('oxolane',   'O:1 C:2 C:3 C:4 C:5', '1-2 2-3 3-4 4-5 5-1', {trivial: 'tetrahydrofuran (THF)'}),
  tpl('pyrrolidine','N:1 C:2 C:3 C:4 C:5', '1-2 2-3 3-4 4-5 5-1'),
  tpl('thiolane',  'S:1 C:2 C:3 C:4 C:5', '1-2 2-3 3-4 4-5 5-1', {trivial: 'tetrahydrothiophene'}),
  tpl('1,3-dioxolane','O:1 C:2 O:3 C:4 C:5', '1-2 2-3 3-4 4-5 5-1'),
  tpl('oxane',     'O:1 C:2 C:3 C:4 C:5 C:6', '1-2 2-3 3-4 4-5 5-6 6-1', {trivial: 'tetrahydropyran (THP)'}),
  tpl('piperidine','N:1 C:2 C:3 C:4 C:5 C:6', '1-2 2-3 3-4 4-5 5-6 6-1'),
  tpl('piperazine','N:1 C:2 C:3 N:4 C:5 C:6', '1-2 2-3 3-4 4-5 5-6 6-1'),
  tpl('morpholine','O:1 C:2 C:3 N:4 C:5 C:6', '1-2 2-3 3-4 4-5 5-6 6-1'),
  tpl('thiane',    'S:1 C:2 C:3 C:4 C:5 C:6', '1-2 2-3 3-4 4-5 5-6 6-1'),
  tpl('1,4-dioxane','O:1 C:2 C:3 O:4 C:5 C:6', '1-2 2-3 3-4 4-5 5-6 6-1'),
  tpl('1,3-dioxane','O:1 C:2 O:3 C:4 C:5 C:6', '1-2 2-3 3-4 4-5 5-6 6-1'),
  tpl('2H-pyran',  'O:1 C:2 C:3 C:4 C:5 C:6', '1-2 2-3 3=4 4-5 5=6 6-1'),
  tpl('4H-pyran',  'O:1 C:2 C:3 C:4 C:5 C:6', '1-2 2=3 3-4 4-5 5=6 6-1'),
  tpl('2H-1-benzopyran', 'O:1 C:2 C:3 C:4 C:4a C:5 C:6 C:7 C:8 C:8a',
      '1-2 2-3 3=4 4-4a 4a-5 5=6 6-7 7=8 8-8a 8a-1 4a=8a', {trivial: '2H-chromene'}),
  tpl('3,4-dihydro-2H-1-benzopyran', 'O:1 C:2 C:3 C:4 C:4a C:5 C:6 C:7 C:8 C:8a',
      '1-2 2-3 3-4 4-4a 4a-5 5=6 6-7 7=8 8-8a 8a-1 4a=8a', {trivial: 'chromane'}),
  tpl('2,3-dihydro-1H-indole', 'N:1 C:2 C:3 C:3a C:4 C:5 C:6 C:7 C:7a',
      '1-2 2-3 3-3a 3a-4 4=5 5-6 6=7 7-7a 7a-1 3a=7a', {trivial: 'indoline'}),
  tpl('1,2,3,4-tetrahydronaphthalene', 'C:1 C:2 C:3 C:4 C:4a C:5 C:6 C:7 C:8 C:8a',
      '1-2 2-3 3-4 4-4a 4a-5 5=6 6-7 7=8 8-8a 8a-1 4a=8a', {trivial: 'tetralin'}),
  tpl('1,3-thiazolidine', 'S:1 C:2 N:3 C:4 C:5', '1-2 2-3 3-4 4-5 5-1'),
  tpl('1,3-oxazolidine',  'O:1 C:2 N:3 C:4 C:5', '1-2 2-3 3-4 4-5 5-1')
];

/** Build the Mol for a template lazily and cache it. */
function templateMol(t) {
  if (t._mol) return t._mol;
  var mol = new Mol(), index = {};
  t.labels = [];
  t.atomSpec.trim().split(/\s+/).forEach(function (item) {
    var parts = item.split(':');
    index[parts[1]] = mol.addAtom(parts[0]);
    t.labels.push(parts[1]);
  });
  t.bondSpec.trim().split(/\s+/).forEach(function (item) {
    var m = /^(.+?)([-=#])(.+)$/.exec(item);
    if (!m) throw new Error('bad template bond: ' + item + ' in ' + t.name);
    var order = m[2] === '=' ? 2 : (m[2] === '#' ? 3 : 1);
    if (index[m[1]] === undefined || index[m[3]] === undefined) {
      throw new Error('unknown template atom in ' + t.name + ': ' + item);
    }
    mol.addBond(index[m[1]], index[m[3]], order);
  });
  t._mol = mol;
  return mol;
}

/* ------------------------------------------------------------------ *
 * Ring-system matching
 * ------------------------------------------------------------------ */

/**
 * Comparison label for ring-system matching.
 *
 * Uses the number of pi bonds an atom has *inside* the ring system, which
 * is the same for every Kekule structure of that system, and deliberately
 * ignores hydrogens and exocyclic double bonds so that, for example,
 * cyclohexanone still matches the cyclohexane skeleton.
 */
function ringLabel(mol, i) {
  var a = mol.atoms[i], piIn = 0;
  var links = mol.links(i);
  for (var k = 0; k < links.length; k++) {
    if (links[k].bond.order >= 2) piIn++;
  }
  return [a.element, a.charge, piIn, mol.degree(i)].join('|');
}

var MATCH_OPTS = {atomLabel: ringLabel, ignoreColors: true, ignoreBondOrder: true, limit: 64};

/**
 * Identify a ring system.
 * Returns null when the skeleton is not one this program can name, else
 *   {name, kind, atoms, numberings, labels, template}
 * where each numbering is a list of {label, atom} in numbering order.
 */
function describeRingSystem(mol, system) {
  var atoms = system.atoms.slice();
  var sub = mol.subMol(atoms);
  var subMol = sub.mol;
  var toMol = sub.rmap;                          // sub index -> mol index
  var i;

  for (i = 0; i < TEMPLATES.length; i++) {
    var t = TEMPLATES[i];
    var tMol = templateMol(t);
    if (tMol.atoms.length !== subMol.atoms.length) continue;
    if (tMol.bonds.length !== subMol.bonds.length) continue;
    var maps = ONG.isomorphisms(tMol, subMol, MATCH_OPTS);
    if (!maps.length) continue;
    var numberings = maps.map(function (map) {
      return t.labels.map(function (label, pos) {
        return {label: label, atom: toMol[map[pos]]};
      });
    });
    return {
      kind: 'template', name: t.name, trivial: t.trivial || null,
      atoms: atoms, numberings: numberings, template: t
    };
  }

  // Plain monocyclic carbocycle: named from its size and unsaturation.
  if (system.rings.length === 1 && atoms.every(function (a) { return mol.atoms[a].element === 'C'; })) {
    var ring = system.rings[0];
    var n = ring.length;
    var numberings = [];
    for (var start = 0; start < n; start++) {
      for (var dir = 1; dir >= -1; dir -= 2) {
        var seq = [];
        for (var k = 0; k < n; k++) {
          seq.push({label: String(k + 1), atom: ring[(start + dir * k + 2 * n) % n]});
        }
        numberings.push(seq);
      }
    }
    // benzene is a retained name, not "cyclohexa-1,3,5-triene"
    var doubleCount = 0;
    for (k = 0; k < n; k++) {
      if (mol.bondOrder(ring[k], ring[(k + 1) % n]) === 2) doubleCount++;
    }
    if (n === 6 && doubleCount === 3) {
      return {kind: 'benzene', name: 'benzene', atoms: atoms, numberings: numberings};
    }
    return {kind: 'carbocycle', size: n, atoms: atoms, numberings: numberings};
  }

  return null;
}

ONG.rings = {
  TEMPLATES: TEMPLATES,
  templateMol: templateMol,
  ringLabel: ringLabel,
  describeRingSystem: describeRingSystem
};

if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
