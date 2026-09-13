/**
 * Stereochemistry (R/S, E/Z) and the newer functional groups.
 */
var ONG = require('./lib.js');

var pass = 0, fail = [];
function check(label, got, want) {
  if (got === want) pass++;
  else fail.push(label + '\n      want: ' + want + '\n      got : ' + got);
}

function generated(smiles) {
  return ONG.generatedNameFor(ONG.parseSmiles(smiles));
}

/* --- R/S ------------------------------------------------------------ */
// Reference structures whose configuration is not in doubt.
var CHIRAL = [
  ['N[C@@H](C)C(=O)O', '(2S)-2-aminopropanoic acid'],              // L-alanine
  ['N[C@H](C)C(=O)O', '(2R)-2-aminopropanoic acid'],               // D-alanine
  ['C([C@@H](C(=O)O)N)O', '(2S)-2-amino-3-hydroxypropanoic acid'], // L-serine
  ['C([C@H](C=O)O)O', '(2R)-2,3-dihydroxypropanal'],               // D-glyceraldehyde
  ['SC[C@H](N)C(=O)O', '(2R)-2-amino-3-sulfanylpropanoic acid'],   // L-cysteine is (R)
  ['CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O',
   '(1R,2S,5R)-5-methyl-2-(propan-2-yl)cyclohexan-1-ol'],          // (-)-menthol
  ['OC(=O)[C@H](O)[C@@H](O)C(=O)O', '(2R,3R)-2,3-dihydroxybutanedioic acid'],
  ['OC(=O)[C@H](O)[C@H](O)C(=O)O', '(2R,3S)-2,3-dihydroxybutanedioic acid'],  // meso
  ['CC[C@H](C)O', '(2S)-butan-2-ol'],
  ['CC[C@@H](C)O', '(2R)-butan-2-ol'],
  ['OC(=O)c1ccccc1[C@H](C)O', '2-[(1S)-1-hydroxyethyl]benzoic acid'],
  ['COC(=O)[C@H](O)C', 'methyl (2R)-2-hydroxypropanoate']
];
CHIRAL.forEach(function (c) { check('R/S ' + c[0], generated(c[0]), c[1]); });

/* --- E/Z ------------------------------------------------------------ */
var GEOMETRIC = [
  ['C/C=C/C', '(2E)-but-2-ene'],
  ['C/C=C\\C', '(2Z)-but-2-ene'],
  ['OC(=O)/C=C/C(=O)O', '(2E)-but-2-enedioic acid'],           // fumaric
  ['OC(=O)/C=C\\C(=O)O', '(2Z)-but-2-enedioic acid'],          // maleic
  ['CCCCCCCC/C=C\\CCCCCCCC(=O)O', '(9Z)-octadec-9-enoic acid'],  // oleic
  ['O=C/C=C/c1ccccc1', '(2E)-3-phenylprop-2-enal'],            // cinnamaldehyde
  ['C/C=C/C=C/C', '(2E,4E)-hexa-2,4-diene']
];
GEOMETRIC.forEach(function (c) { check('E/Z ' + c[0], generated(c[0]), c[1]); });

/* --- a flat drawing gets no descriptors, but is flagged ------------- */
check('no stereo drawn', generated('CC(C)C(O)CC'), 'hexan-3-ol'.replace('hexan-3-ol', '2-methylpentan-3-ol'));
var flat = ONG.nameStructure(ONG.parseSmiles('NC(C)C(=O)O'));
check('undefined centre is reported',
  flat.warnings.some(function (w) { return w.indexOf('입체중심') >= 0; }), true);
var drawn = ONG.nameStructure(ONG.parseSmiles('N[C@@H](C)C(=O)O'));
check('defined centre is not reported',
  drawn.warnings.some(function (w) { return w.indexOf('입체중심') >= 0; }), false);

/* --- rings and terminal double bonds are not stereogenic ------------ */
check('cyclohexene has no E/Z',
  ONG.nameStructure(ONG.parseSmiles('C1=CCCCC1')).warnings.length, 0);
check('1,1-disubstituted alkene has no E/Z',
  ONG.nameStructure(ONG.parseSmiles('CC(C)=CC')).warnings.length, 0);

/* --- newer functional groups ---------------------------------------- */
var GROUPS = [
  ['OS(=O)(=O)c1ccccc1', 'benzenesulfonic acid'],
  ['CS(=O)(=O)O', 'methanesulfonic acid'],
  ['Cc1ccc(cc1)S(=O)(=O)Cl', '4-methylbenzene-1-sulfonyl chloride'],
  ['Nc1ccc(cc1)S(=O)(=O)N', '4-aminobenzene-1-sulfonamide'],
  ['COS(C)(=O)=O', 'methyl methanesulfonate'],
  ['NCCS(=O)(=O)O', '2-aminoethane-1-sulfonic acid'],
  ['OC(=O)c1ccc(cc1)S(=O)(=O)O', '4-sulfobenzoic acid'],
  ['CS(C)=O', '(methanesulfinyl)methane'],
  ['CS(C)(=O)=O', '(methanesulfonyl)methane'],
  ['CSc1ccccc1', '(methylsulfanyl)benzene'],
  ['CC(=O)OC(C)=O', 'ethanoic anhydride'],
  ['CC(=O)OC(=O)CC', 'ethanoic propanoic anhydride'],
  ['O=C1CCC(=O)O1', 'oxolane-2,5-dione'],
  ['CC=N', 'ethanimine'],
  ['CC=NC', 'N-methylethanimine'],
  ['CC(C)=S', 'propane-2-thione'],
  ['C[N+](C)(C)C', 'N,N,N-trimethylmethanaminium'],
  ['[N-]=[N+]=Nc1ccccc1', 'azidobenzene'],
  ['O=C=Nc1ccccc1', 'isocyanatobenzene'],
  ['S=C=Nc1ccccc1', 'isothiocyanatobenzene'],
  ['O=Nc1ccccc1', 'nitrosobenzene'],
  ['N#CSc1ccccc1', 'thiocyanatobenzene']
];
GROUPS.forEach(function (c) { check('group ' + c[0], generated(c[0]), c[1]); });

console.log('stereo + groups: ' + pass + '/' + (pass + fail.length) + ' passed');
fail.forEach(function (f) { console.log('  ✗ ' + f); });
if (fail.length) process.exitCode = 1;
