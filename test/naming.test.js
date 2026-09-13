/**
 * Expected-name regression tests. Names follow the 2013 IUPAC
 * recommendations (preferred IUPAC names) unless noted.
 */
var ONG = require('./lib.js');

var CASES = [
  // --- alkanes -----------------------------------------------------
  ['CCCC', 'butane'],
  ['CC(C)C', '2-methylpropane'],
  ['CC(C)(C)C', '2,2-dimethylpropane'],
  ['CCC(C)C', '2-methylbutane'],
  ['CC(C)C(C)C', '2,3-dimethylbutane'],
  ['CCCCCCC', 'heptane'],
  ['CCC(CC)CCC', '3-ethylhexane'],
  ['CC(C)CC(C)(C)C', '2,2,4-trimethylpentane'],
  // --- unsaturation ------------------------------------------------
  ['C=C', 'ethene'],
  ['CC=C', 'prop-1-ene'],
  ['CC=CC', 'but-2-ene'],
  ['C=CC=C', 'buta-1,3-diene'],
  ['C#C', 'ethyne'],
  ['CC#CC', 'but-2-yne'],
  ['C=CC#C', 'but-1-en-3-yne'],
  ['CC(C)=C', '2-methylprop-1-ene'],
  // --- halides -----------------------------------------------------
  ['ClCCl', 'dichloromethane'],
  ['ClC(Cl)(Cl)Cl', 'tetrachloromethane'],
  ['CCCl', 'chloroethane'],
  ['CCCBr', '1-bromopropane'],
  ['CC(C)(C)Cl', '2-chloro-2-methylpropane'],
  ['FC(F)(F)C(Cl)Br', '2-bromo-2-chloro-1,1,1-trifluoroethane'],
  // --- alcohols ----------------------------------------------------
  ['CO', 'methanol'],
  ['CCO', 'ethanol'],
  ['CCCO', 'propan-1-ol'],
  ['CC(C)O', 'propan-2-ol'],
  ['CC(C)(C)O', '2-methylpropan-2-ol'],
  ['OCCO', 'ethane-1,2-diol'],
  ['OCC(O)CO', 'propane-1,2,3-triol'],
  ['CC(C)CCO', '3-methylbutan-1-ol'],
  ['C1CCCCC1O', 'cyclohexanol'],
  // --- carbonyls ---------------------------------------------------
  ['CC=O', 'ethanal'],
  ['CCC=O', 'propanal'],
  ['CC(C)=O', 'propan-2-one'],
  ['CCC(C)=O', 'butan-2-one'],
  ['O=C1CCCCC1', 'cyclohexanone'],
  ['O=CCCC=O', 'butanedial'],
  ['CC(=O)CC(C)=O', 'pentane-2,4-dione'],
  // --- acids and derivatives ---------------------------------------
  ['CC(=O)O', 'ethanoic acid'],
  ['CCC(=O)O', 'propanoic acid'],
  ['OC(=O)CCC(=O)O', 'butanedioic acid'],
  ['CC(O)C(=O)O', '2-hydroxypropanoic acid'],
  ['CC(=O)OCC', 'ethyl ethanoate'],
  ['CCC(=O)OC', 'methyl propanoate'],
  ['CCC(=O)Cl', 'propanoyl chloride'],
  ['CC(=O)N', 'ethanamide'],
  ['CCC#N', 'propanenitrile'],
  ['OC(=O)CC(O)(CC(=O)O)C(=O)O', '2-hydroxypropane-1,2,3-tricarboxylic acid'],
  // --- amines and ethers -------------------------------------------
  ['CCN', 'ethanamine'],
  ['CCNCC', 'N-ethylethanamine'],
  ['CN(C)C', 'N,N-dimethylmethanamine'],
  ['CCOCC', 'ethoxyethane'],
  ['COC', 'methoxymethane'],
  ['CCOC', 'methoxyethane'],
  ['CS', 'methanethiol'],
  // --- rings -------------------------------------------------------
  ['C1CCCCC1', 'cyclohexane'],
  ['C1CC1', 'cyclopropane'],
  ['C1=CCCCC1', 'cyclohexene'],
  ['CC1CCCCC1', 'methylcyclohexane'],
  ['c1ccccc1', 'benzene'],
  ['Cc1ccccc1', 'methylbenzene'],
  ['Cc1ccc(C)cc1', '1,4-dimethylbenzene'],
  ['Clc1ccccc1Cl', '1,2-dichlorobenzene'],
  ['Oc1ccccc1', 'phenol'],
  ['Nc1ccccc1', 'aniline'],
  ['O=Cc1ccccc1', 'benzaldehyde'],
  ['OC(=O)c1ccccc1', 'benzoic acid'],
  ['Oc1ccc(Cl)cc1', '4-chlorophenol'],
  ['C=Cc1ccccc1', 'ethenylbenzene'],
  ['CC(=O)c1ccccc1', '1-phenylethan-1-one'],
  ['OCCc1ccccc1', '2-phenylethan-1-ol'],
  ['OC(=O)C1CCCCC1', 'cyclohexanecarboxylic acid'],
  ['c1ccc2ccccc2c1', 'naphthalene'],
  ['Cc1cccc2ccccc12', '1-methylnaphthalene'],
  ['c1ccncc1', 'pyridine'],
  ['OC(=O)c1cccnc1', 'pyridine-3-carboxylic acid'],
  ['C1CCOC1', 'oxolane'],
  ['C1COCCN1', 'morpholine'],
  ['c1ccc2[nH]ccc2c1', '1H-indole'],
  // --- multi-group / real drugs ------------------------------------
  ['NCC(=O)O', '2-aminoethanoic acid'],
  ['CC(=O)Oc1ccccc1C(=O)O', '2-(acetyloxy)benzoic acid'],
  ['CC(C)Cc1ccc(cc1)C(C)C(=O)O', '2-[4-(2-methylpropyl)phenyl]propanoic acid'],
  ['CC(=O)Nc1ccc(O)cc1', 'N-(4-hydroxyphenyl)ethanamide'],
  ['COc1cc(C=O)ccc1O', '4-hydroxy-3-methoxybenzaldehyde'],
  ['Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', '2-methyl-1,3,5-trinitrobenzene'],
  ['[O-][N+](=O)c1ccccc1', 'nitrobenzene'],
  ['OCC(O)C(O)C(O)C(O)C=O', '2,3,4,5,6-pentahydroxyhexanal']
];

var pass = 0, fail = [];
CASES.forEach(function (c) {
  var got;
  try {
    got = ONG.nameMolecule(ONG.parseSmiles(c[0])).name;
  } catch (e) {
    got = 'ERROR: ' + e.message;
  }
  if (got === c[1]) pass++;
  else fail.push({smiles: c[0], want: c[1], got: got});
});

console.log('naming: ' + pass + '/' + CASES.length + ' passed');
fail.forEach(function (f) {
  console.log('  ✗ ' + f.smiles + '\n      want: ' + f.want + '\n      got : ' + f.got);
});
if (fail.length) process.exitCode = 1;
