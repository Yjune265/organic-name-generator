/** Alternative-name and compound-database tests. */
var ONG = require('./lib.js');

var pass = 0, fail = [];

function check(label, ok, detail) {
  if (ok) pass++;
  else fail.push(label + (detail ? ' — ' + detail : ''));
}

/* --- the database is self-consistent ------------------------------- */
// Stereoisomers share a constitution key on purpose, so entries may only
// collide when the engine can still tell them apart by name.
var keys = {};
ONG.COMPOUNDS.forEach(function (entry) {
  var mol;
  try {
    mol = ONG.parseSmiles(entry.smiles);
  } catch (e) {
    check('parse ' + entry.smiles, false, e.message);
    return;
  }
  var generated = ONG.generatedNameFor(mol);
  var key = mol.structureKey() + '|' + generated;
  check('distinguishable ' + entry.smiles, !keys[key], 'indistinguishable from ' + keys[key]);
  keys[key] = entry.smiles;
  var found = ONG.lookupCompound(mol, generated);
  check('lookup ' + entry.smiles, found && found.smiles === entry.smiles,
        'got ' + (found && found.smiles));
  check('has a name ' + entry.smiles, !!(entry.iupac || (entry.common && entry.common.length)));
});

/* --- alternative names --------------------------------------------- */
var VARIANTS = [
  ['CC(C)O', ['2-propanol', 'isopropyl alcohol']],
  ['CC(C)C', ['isobutane']],
  ['CC(C)(C)C', ['neopentane']],
  ['CC(C)c1ccccc1', ['isopropylbenzene']],
  ['CC(C)(C)c1ccccc1', ['tert-butylbenzene']],
  ['CCC(C)O', ['2-butanol', 'sec-butyl alcohol']],
  ['CCOCC', ['diethyl ether']],
  ['COCC', ['ethyl methyl ether']],
  ['CCNCC', ['diethylamine']],
  ['CC(=O)OCC', ['ethyl acetate']],
  ['CC(C)CC(C)(C)C', ['isooctane']],
  ['C=CCO', ['allyl alcohol', '2-propen-1-ol']],
  ['CC(C)Cc1ccc(cc1)C(C)C(=O)O', ['2-(4-isobutylphenyl)propanoic acid']],
  ['CC(C)CO', ['isobutyl alcohol', '2-methyl-1-propanol']],
  ['CC(C)Cl', ['isopropyl chloride']]
];

VARIANTS.forEach(function (c) {
  var mol = ONG.parseSmiles(c[0]);
  var result = ONG.nameStructure(mol);
  var names = result.alternatives.map(function (a) { return a.name; })
    .concat(result.common || []);
  c[1].forEach(function (want) {
    check('variant ' + c[0] + ' -> ' + want, names.indexOf(want) >= 0, 'got: ' + names.join(', '));
  });
});

/* --- database-backed answers --------------------------------------- */
var DB_CASES = [
  ['CN1C=NC2=C1C(=O)N(C)C(=O)N2C', '1,3,7-trimethylpurine-2,6-dione', 'caffeine'],
  ['CC(=O)Oc1ccccc1C(=O)O', '2-(acetyloxy)benzoic acid', 'aspirin'],
  ['CC1(C)C2CCC1(C)C(=O)C2', '1,7,7-trimethylbicyclo[2.2.1]heptan-2-one', 'camphor'],
  ['CN1CCCC1c1cccnc1', '3-(1-methylpyrrolidin-2-yl)pyridine', 'nicotine']
];
DB_CASES.forEach(function (c) {
  var r = ONG.nameStructure(ONG.parseSmiles(c[0]));
  check('db name ' + c[2], r.name === c[1], 'got ' + r.name);
  check('db common ' + c[2], (r.common || []).indexOf(c[2]) >= 0, 'got ' + (r.common || []).join(','));
  // Only "cannot name this skeleton" complaints are noise once a curated
  // name is used; an undefined-stereocentre note is still worth showing.
  var noisy = r.warnings.filter(function (w) { return w.indexOf('고리 골격') >= 0; });
  check('db quiet ' + c[2], noisy.length === 0, 'warnings: ' + r.warnings.join('|'));
});

/* --- structures drawn differently still match ----------------------- */
var SAME = [
  ['c1ccccc1', 'C1=CC=CC=C1'],
  ['CN1C=NC2=C1C(=O)N(C)C(=O)N2C', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'],
  ['OC(=O)c1ccccc1', 'c1ccccc1C(O)=O']
];
SAME.forEach(function (pair) {
  var a = ONG.nameStructure(ONG.parseSmiles(pair[0]));
  var b = ONG.nameStructure(ONG.parseSmiles(pair[1]));
  check('same name for ' + pair[0] + ' / ' + pair[1], a.name === b.name, a.name + ' vs ' + b.name);
});

console.log('variants + database: ' + pass + '/' + (pass + fail.length) + ' passed');
fail.forEach(function (f) { console.log('  ✗ ' + f); });
if (fail.length) process.exitCode = 1;
