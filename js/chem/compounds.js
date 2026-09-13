/**
 * compounds.js - well known compounds.
 *
 * Each entry is looked up by exact structure match (constitution only -
 * stereochemistry is not compared), so a drawing is recognised however it
 * was drawn. `iupac` is a curated name; it is trusted over the generated
 * one, which lets the app name things the engine cannot derive by itself
 * (fused polycyclics, steroids, bridged rings...).
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

var COMPOUNDS = [
  /* --- hydrocarbons ------------------------------------------------ */
  {smiles: 'C', iupac: 'methane', common: ['메테인 (메탄)', 'marsh gas']},
  {smiles: 'CC', iupac: 'ethane', common: ['에테인 (에탄)']},
  {smiles: 'CCC', iupac: 'propane', common: ['프로페인 (프로판)']},
  {smiles: 'CC(C)C', iupac: '2-methylpropane', common: ['isobutane', '아이소뷰테인']},
  {smiles: 'CC(C)(C)C', iupac: '2,2-dimethylpropane', common: ['neopentane', '네오펜테인']},
  {smiles: 'CCC(C)C', iupac: '2-methylbutane', common: ['isopentane', 'isoamylhydride']},
  {smiles: 'CC(C)CC(C)(C)C', iupac: '2,2,4-trimethylpentane', common: ['isooctane', '아이소옥테인']},
  {smiles: 'C=C', iupac: 'ethene', common: ['ethylene', '에틸렌']},
  {smiles: 'CC=C', iupac: 'prop-1-ene', common: ['propylene', '프로필렌']},
  {smiles: 'C#C', iupac: 'ethyne', common: ['acetylene', '아세틸렌']},
  {smiles: 'CC(=C)C=C', iupac: '2-methylbuta-1,3-diene', common: ['isoprene', '아이소프렌']},
  {smiles: 'c1ccccc1', iupac: 'benzene', common: ['벤젠']},
  {smiles: 'Cc1ccccc1', iupac: 'methylbenzene', common: ['toluene', '톨루엔']},
  {smiles: 'Cc1ccccc1C', iupac: '1,2-dimethylbenzene', common: ['o-xylene', '오쏘-자일렌']},
  {smiles: 'Cc1cccc(C)c1', iupac: '1,3-dimethylbenzene', common: ['m-xylene', '메타-자일렌']},
  {smiles: 'Cc1ccc(C)cc1', iupac: '1,4-dimethylbenzene', common: ['p-xylene', '파라-자일렌']},
  {smiles: 'C=Cc1ccccc1', iupac: 'ethenylbenzene', common: ['styrene', '스타이렌', 'vinylbenzene']},
  {smiles: 'CC(C)c1ccccc1', iupac: '(propan-2-yl)benzene', common: ['cumene', '큐멘', 'isopropylbenzene']},
  {smiles: 'Cc1ccc(C(C)C)cc1', iupac: '1-methyl-4-(propan-2-yl)benzene', common: ['p-cymene']},
  {smiles: 'c1ccc2ccccc2c1', iupac: 'naphthalene', common: ['나프탈렌']},
  {smiles: 'c1ccc2cc3ccccc3cc2c1', iupac: 'anthracene', common: ['안트라센']},
  {smiles: 'c1ccc2c(c1)ccc1ccccc12', iupac: 'phenanthrene', common: ['페난트렌']},
  {smiles: 'C1CCCCC1', iupac: 'cyclohexane', common: ['사이클로헥세인']},

  /* --- alcohols, phenols, ethers ----------------------------------- */
  {smiles: 'CO', iupac: 'methanol', common: ['methyl alcohol', '메탄올', 'wood alcohol']},
  {smiles: 'CCO', iupac: 'ethanol', common: ['ethyl alcohol', '에탄올', 'grain alcohol']},
  {smiles: 'CCCO', iupac: 'propan-1-ol', common: ['n-propyl alcohol']},
  {smiles: 'CC(C)O', iupac: 'propan-2-ol', common: ['isopropyl alcohol', 'IPA', '아이소프로판올']},
  {smiles: 'CC(C)(C)O', iupac: '2-methylpropan-2-ol', common: ['tert-butyl alcohol', 'tert-butanol']},
  {smiles: 'OCCO', iupac: 'ethane-1,2-diol', common: ['ethylene glycol', '에틸렌 글라이콜']},
  {smiles: 'OCC(O)CO', iupac: 'propane-1,2,3-triol', common: ['glycerol', 'glycerin', '글리세롤']},
  {smiles: 'Oc1ccccc1', iupac: 'phenol', common: ['carbolic acid', '페놀']},
  {smiles: 'Oc1ccccc1O', iupac: 'benzene-1,2-diol', common: ['catechol', 'pyrocatechol']},
  {smiles: 'Oc1cccc(O)c1', iupac: 'benzene-1,3-diol', common: ['resorcinol']},
  {smiles: 'Oc1ccc(O)cc1', iupac: 'benzene-1,4-diol', common: ['hydroquinone', '하이드로퀴논']},
  {smiles: 'Cc1ccc(O)cc1', iupac: '4-methylphenol', common: ['p-cresol', '파라-크레졸']},
  {smiles: 'COc1ccccc1', iupac: 'methoxybenzene', common: ['anisole', '아니솔']},
  {smiles: 'CCOCC', iupac: 'ethoxyethane', common: ['diethyl ether', '다이에틸 에터', 'ether']},
  {smiles: 'C1CCOC1', iupac: 'oxolane', common: ['tetrahydrofuran', 'THF']},
  {smiles: 'C1CO1', iupac: 'oxirane', common: ['ethylene oxide', '에틸렌 옥사이드']},
  {smiles: 'C1COCCO1', iupac: '1,4-dioxane', common: ['dioxane']},

  /* --- carbonyl compounds ------------------------------------------ */
  {smiles: 'C=O', iupac: 'methanal', common: ['formaldehyde', '폼알데하이드']},
  {smiles: 'CC=O', iupac: 'ethanal', common: ['acetaldehyde', '아세트알데하이드']},
  {smiles: 'O=Cc1ccccc1', iupac: 'benzaldehyde', common: ['벤즈알데하이드', 'oil of bitter almond']},
  {smiles: 'CC(C)=O', iupac: 'propan-2-one', common: ['acetone', '아세톤', 'dimethyl ketone']},
  {smiles: 'CCC(C)=O', iupac: 'butan-2-one', common: ['methyl ethyl ketone', 'MEK']},
  {smiles: 'O=C1CCCCC1', iupac: 'cyclohexanone', common: ['사이클로헥산온']},
  {smiles: 'CC(=O)c1ccccc1', iupac: '1-phenylethan-1-one', common: ['acetophenone', '아세토페논']},
  {smiles: 'O=C(c1ccccc1)c1ccccc1', iupac: 'diphenylmethanone', common: ['benzophenone']},
  {smiles: 'O=C1C=CC(=O)C=C1', iupac: 'cyclohexa-2,5-diene-1,4-dione', common: ['1,4-benzoquinone', 'quinone']},

  /* --- acids and derivatives --------------------------------------- */
  {smiles: 'OC=O', iupac: 'methanoic acid', common: ['formic acid', '폼산 (개미산)']},
  {smiles: 'CC(=O)O', iupac: 'ethanoic acid', common: ['acetic acid', '아세트산 (초산)']},
  {smiles: 'CCC(=O)O', iupac: 'propanoic acid', common: ['propionic acid', '프로피온산']},
  {smiles: 'CCCC(=O)O', iupac: 'butanoic acid', common: ['butyric acid', '뷰티르산']},
  {smiles: 'CCCCCCCCCCCCCCCCCC(=O)O', iupac: 'octadecanoic acid', common: ['stearic acid', '스테아르산']},
  {smiles: 'CCCCCCCCCCCCCCCC(=O)O', iupac: 'hexadecanoic acid', common: ['palmitic acid', '팔미트산']},
  {smiles: 'OC(=O)C(=O)O', iupac: 'ethanedioic acid', common: ['oxalic acid', '옥살산']},
  {smiles: 'OC(=O)CC(=O)O', iupac: 'propanedioic acid', common: ['malonic acid']},
  {smiles: 'OC(=O)CCC(=O)O', iupac: 'butanedioic acid', common: ['succinic acid', '숙신산']},
  {smiles: 'OC(=O)CCCC(=O)O', iupac: 'pentanedioic acid', common: ['glutaric acid']},
  {smiles: 'OC(=O)CCCCC(=O)O', iupac: 'hexanedioic acid', common: ['adipic acid', '아디프산']},
  {smiles: 'OC(=O)C=CC(=O)O', iupac: 'but-2-enedioic acid',
   common: ['maleic acid (cis)', 'fumaric acid (trans)'],
   note: '이 프로그램은 cis/trans(E/Z)를 구별하지 않아 두 이름을 모두 보여줍니다.'},
  {smiles: 'CC(O)C(=O)O', iupac: '2-hydroxypropanoic acid', common: ['lactic acid', '젖산']},
  {smiles: 'OC(=O)CC(O)(CC(=O)O)C(=O)O', iupac: '2-hydroxypropane-1,2,3-tricarboxylic acid',
   common: ['citric acid', '시트르산 (구연산)']},
  {smiles: 'OC(=O)C(O)C(O)C(=O)O', iupac: '2,3-dihydroxybutanedioic acid', common: ['tartaric acid', '타타르산']},
  {smiles: 'OC(=O)c1ccccc1', iupac: 'benzoic acid', common: ['벤조산']},
  {smiles: 'OC(=O)c1ccccc1O', iupac: '2-hydroxybenzoic acid', common: ['salicylic acid', '살리실산']},
  {smiles: 'CC(=O)Oc1ccccc1C(=O)O', iupac: '2-(acetyloxy)benzoic acid',
   common: ['aspirin', '아스피린', 'acetylsalicylic acid']},
  {smiles: 'OC(=O)c1ccccc1C(=O)O', iupac: 'benzene-1,2-dicarboxylic acid', common: ['phthalic acid']},
  {smiles: 'CC(=O)OCC', iupac: 'ethyl ethanoate', common: ['ethyl acetate', '에틸 아세테이트']},
  {smiles: 'CC(=O)N', iupac: 'ethanamide', common: ['acetamide', '아세트아마이드']},
  {smiles: 'NC(N)=O', iupac: 'urea', common: ['carbamide', '요소']},
  {smiles: 'CC#N', iupac: 'ethanenitrile', common: ['acetonitrile', '아세토나이트릴']},
  {smiles: 'CN(C)C=O', iupac: 'N,N-dimethylmethanamide', common: ['N,N-dimethylformamide', 'DMF']},

  /* --- nitrogen and sulfur compounds -------------------------------- */
  {smiles: 'CN', iupac: 'methanamine', common: ['methylamine', '메틸아민']},
  {smiles: 'Nc1ccccc1', iupac: 'aniline', common: ['아닐린', 'phenylamine']},
  {smiles: '[O-][N+](=O)c1ccccc1', iupac: 'nitrobenzene', common: ['나이트로벤젠']},
  {smiles: 'Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]',
   iupac: '2-methyl-1,3,5-trinitrobenzene', common: ['TNT', 'trinitrotoluene']},
  {smiles: 'CS(C)=O', iupac: '(methanesulfinyl)methane', common: ['dimethyl sulfoxide', 'DMSO']},
  {smiles: 'CS', iupac: 'methanethiol', common: ['methyl mercaptan']},
  {smiles: 'c1ccncc1', iupac: 'pyridine', common: ['피리딘']},
  {smiles: 'c1cc[nH]c1', iupac: '1H-pyrrole', common: ['피롤']},
  {smiles: 'c1ccoc1', iupac: 'furan', common: ['퓨란']},
  {smiles: 'c1ccsc1', iupac: 'thiophene', common: ['싸이오펜']},
  {smiles: 'c1c[nH]cn1', iupac: '1H-imidazole', common: ['이미다졸']},
  {smiles: 'c1ccc2[nH]ccc2c1', iupac: '1H-indole', common: ['인돌']},

  /* --- halogen compounds -------------------------------------------- */
  {smiles: 'ClCCl', iupac: 'dichloromethane', common: ['methylene chloride', 'DCM']},
  {smiles: 'ClC(Cl)Cl', iupac: 'trichloromethane', common: ['chloroform', '클로로폼']},
  {smiles: 'ClC(Cl)(Cl)Cl', iupac: 'tetrachloromethane', common: ['carbon tetrachloride', '사염화탄소']},
  {smiles: 'ClC(Cl)(F)F', iupac: 'dichlorodifluoromethane', common: ['Freon-12', 'CFC-12']},
  {smiles: 'FC(F)(F)C(Cl)Br', iupac: '2-bromo-2-chloro-1,1,1-trifluoroethane', common: ['halothane']},
  {smiles: 'C=CCl', iupac: 'chloroethene', common: ['vinyl chloride', '염화 바이닐']},

  /* --- amino acids --------------------------------------------------- */
  {smiles: 'NCC(=O)O', iupac: '2-aminoethanoic acid', common: ['glycine', 'Gly', '글라이신']},
  {smiles: 'CC(N)C(=O)O', iupac: '2-aminopropanoic acid', common: ['alanine', 'Ala', '알라닌']},
  {smiles: 'CC(C)C(N)C(=O)O', iupac: '2-amino-3-methylbutanoic acid', common: ['valine', 'Val']},
  {smiles: 'CC(C)CC(N)C(=O)O', iupac: '2-amino-4-methylpentanoic acid', common: ['leucine', 'Leu']},
  {smiles: 'CCC(C)C(N)C(=O)O', iupac: '2-amino-3-methylpentanoic acid', common: ['isoleucine', 'Ile']},
  {smiles: 'OCC(N)C(=O)O', iupac: '2-amino-3-hydroxypropanoic acid', common: ['serine', 'Ser']},
  {smiles: 'SCC(N)C(=O)O', iupac: '2-amino-3-sulfanylpropanoic acid', common: ['cysteine', 'Cys']},
  {smiles: 'CC(O)C(N)C(=O)O', iupac: '2-amino-3-hydroxybutanoic acid', common: ['threonine', 'Thr']},
  {smiles: 'CSCCC(N)C(=O)O', iupac: '2-amino-4-(methylsulfanyl)butanoic acid', common: ['methionine', 'Met']},
  {smiles: 'OC(=O)C(N)Cc1ccccc1', iupac: '2-amino-3-phenylpropanoic acid', common: ['phenylalanine', 'Phe']},
  {smiles: 'OC(=O)C(N)Cc1ccc(O)cc1', iupac: '2-amino-3-(4-hydroxyphenyl)propanoic acid', common: ['tyrosine', 'Tyr']},
  {smiles: 'OC(=O)C(N)Cc1c[nH]c2ccccc12', iupac: '2-amino-3-(1H-indol-3-yl)propanoic acid', common: ['tryptophan', 'Trp']},
  {smiles: 'OC(=O)CC(N)C(=O)O', iupac: '2-aminobutanedioic acid', common: ['aspartic acid', 'Asp']},
  {smiles: 'OC(=O)CCC(N)C(=O)O', iupac: '2-aminopentanedioic acid', common: ['glutamic acid', 'Glu', '글루탐산']},
  {smiles: 'NCCCCC(N)C(=O)O', iupac: '2,6-diaminohexanoic acid', common: ['lysine', 'Lys']},
  {smiles: 'OC(=O)C1CCCN1', iupac: 'pyrrolidine-2-carboxylic acid', common: ['proline', 'Pro']},
  {smiles: 'OC(=O)C(N)Cc1c[nH]cn1', iupac: '2-amino-3-(1H-imidazol-4-yl)propanoic acid', common: ['histidine', 'His']},

  /* --- natural products and drugs ------------------------------------ */
  {smiles: 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C', iupac: '1,3,7-trimethylpurine-2,6-dione',
   common: ['caffeine', '카페인'],
   note: '완전한 형태: 1,3,7-trimethyl-3,7-dihydro-1H-purine-2,6-dione'},
  {smiles: 'CN1C=NC2=C1C(=O)NC(=O)N2C', iupac: '3,7-dimethylpurine-2,6-dione',
   common: ['theobromine', '테오브로민']},
  {smiles: 'CN1C(=O)N(C)c2[nH]cnc2C1=O', iupac: '1,3-dimethylpurine-2,6-dione',
   common: ['theophylline', '테오필린']},
  {smiles: 'CN1CCCC1c1cccnc1', iupac: '3-(1-methylpyrrolidin-2-yl)pyridine',
   common: ['nicotine', '니코틴']},
  {smiles: 'CC(=O)Nc1ccc(O)cc1', iupac: 'N-(4-hydroxyphenyl)ethanamide',
   common: ['paracetamol', 'acetaminophen', '아세트아미노펜', 'Tylenol']},
  {smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', iupac: '2-[4-(2-methylpropyl)phenyl]propanoic acid',
   common: ['ibuprofen', '이부프로펜']},
  {smiles: 'COc1cc(C=O)ccc1O', iupac: '4-hydroxy-3-methoxybenzaldehyde',
   common: ['vanillin', '바닐린']},
  {smiles: 'C=CCc1ccc(O)c(OC)c1', iupac: '2-methoxy-4-(prop-2-en-1-yl)phenol',
   common: ['eugenol', '오이게놀']},
  {smiles: 'CC(C)c1ccc(C)cc1O', iupac: '5-methyl-2-(propan-2-yl)phenol', common: ['carvacrol']},
  {smiles: 'CC(C)c1ccc(C)c(O)c1', iupac: '2-methyl-5-(propan-2-yl)phenol', common: ['thymol', '타이몰']},
  {smiles: 'CC(C)C1CCC(C)CC1O', iupac: '5-methyl-2-(propan-2-yl)cyclohexan-1-ol',
   common: ['menthol', '멘톨']},
  {smiles: 'CC(=C)C1CCC(C)=CC1', iupac: '1-methyl-4-(prop-1-en-2-yl)cyclohex-1-ene',
   common: ['limonene', '리모넨']},
  {smiles: 'CC1(C)C2CCC1(C)C(=O)C2', iupac: '1,7,7-trimethylbicyclo[2.2.1]heptan-2-one',
   common: ['camphor', '캠퍼 (장뇌)']},
  {smiles: 'O=Cc1ccccc1C=C', iupac: '2-ethenylbenzaldehyde', common: []},
  {smiles: 'NCCc1ccc(O)c(O)c1', iupac: '4-(2-aminoethyl)benzene-1,2-diol',
   common: ['dopamine', '도파민']},
  {smiles: 'CNCC(O)c1ccc(O)c(O)c1', iupac: '4-[1-hydroxy-2-(methylamino)ethyl]benzene-1,2-diol',
   common: ['adrenaline', 'epinephrine', '아드레날린']},
  {smiles: 'NCCc1c[nH]c2ccc(O)cc12', iupac: '3-(2-aminoethyl)-1H-indol-5-ol',
   common: ['serotonin', '세로토닌']},
  {smiles: 'COc1ccc2[nH]cc(CCNC(C)=O)c2c1', iupac: 'N-[2-(5-methoxy-1H-indol-3-yl)ethyl]acetamide',
   common: ['melatonin', '멜라토닌']},
  {smiles: 'NCCc1c[nH]cn1', iupac: '2-(1H-imidazol-4-yl)ethan-1-amine',
   common: ['histamine', '히스타민']},
  {smiles: 'OCC(O)C(O)C(O)C(O)C=O', iupac: '2,3,4,5,6-pentahydroxyhexanal',
   common: ['glucose', '포도당', 'dextrose'],
   note: '입체(D/L, α/β)는 구분하지 않습니다. D-글루코스의 완전한 이름은 (2R,3S,4R,5R)-2,3,4,5,6-pentahydroxyhexanal 입니다.'},
  {smiles: 'OCC(O)C(O)C(O)C(=O)CO', iupac: '1,3,4,5,6-pentahydroxyhexan-2-one',
   common: ['fructose', '과당'],
   note: '입체는 구분하지 않습니다.'},
  {smiles: 'CC(=O)OC1=CC=CC=C1C(=O)OC', iupac: 'methyl 2-(acetyloxy)benzoate', common: []},
  {smiles: 'OC(=O)CC(C)(O)CC(=O)O', iupac: '3-hydroxy-3-methylpentanedioic acid', common: []},
  {smiles: 'CCCCCCCCCCCCCCCCCC(=O)OCC(O)CO', iupac: '2,3-dihydroxypropyl octadecanoate',
   common: ['glyceryl monostearate']},

  /* --- sulfur and nitrogen functional groups -------------------------- */
  {smiles: 'OS(=O)(=O)c1ccccc1', iupac: 'benzenesulfonic acid', common: ['벤젠술폰산']},
  {smiles: 'CS(=O)(=O)O', iupac: 'methanesulfonic acid', common: ['메실산', 'MsOH']},
  {smiles: 'Cc1ccc(cc1)S(=O)(=O)O', iupac: '4-methylbenzene-1-sulfonic acid',
   common: ['p-toluenesulfonic acid', 'TsOH', '토실산']},
  {smiles: 'Cc1ccc(cc1)S(=O)(=O)Cl', iupac: '4-methylbenzene-1-sulfonyl chloride',
   common: ['tosyl chloride', 'TsCl']},
  {smiles: 'Nc1ccc(cc1)S(=O)(=O)N', iupac: '4-aminobenzene-1-sulfonamide',
   common: ['sulfanilamide', '설파닐아마이드']},
  {smiles: 'NCCS(=O)(=O)O', iupac: '2-aminoethane-1-sulfonic acid', common: ['taurine', '타우린']},
  {smiles: 'CS(C)(=O)=O', iupac: '(methanesulfonyl)methane', common: ['dimethyl sulfone', 'MSM']},
  {smiles: 'CSc1ccccc1', iupac: '(methylsulfanyl)benzene', common: ['thioanisole']},
  {smiles: 'CC(=O)OC(C)=O', iupac: 'ethanoic anhydride', common: ['acetic anhydride', '무수 아세트산']},
  {smiles: 'O=C1CCC(=O)O1', iupac: 'oxolane-2,5-dione', common: ['succinic anhydride', '무수 숙신산']},
  {smiles: 'O=C1OC(=O)c2ccccc12', iupac: '2-benzofuran-1,3-dione',
   common: ['phthalic anhydride', '무수 프탈산']},
  {smiles: 'O=C=Nc1ccccc1', iupac: 'isocyanatobenzene', common: ['phenyl isocyanate']},
  {smiles: 'O=Nc1ccccc1', iupac: 'nitrosobenzene', common: ['나이트로소벤젠']},
  {smiles: 'C[N+](C)(C)C', iupac: 'N,N,N-trimethylmethanaminium',
   common: ['tetramethylammonium', '테트라메틸암모늄']},

  /* --- stereochemistry: the isomer decides the common name ------------ */
  {smiles: 'N[C@@H](C)C(=O)O', iupac: '(2S)-2-aminopropanoic acid',
   common: ['L-alanine', 'L-Ala', 'L-알라닌']},
  {smiles: 'N[C@H](C)C(=O)O', iupac: '(2R)-2-aminopropanoic acid',
   common: ['D-alanine', 'D-알라닌']},
  {smiles: 'OC[C@H](N)C(=O)O', iupac: '(2S)-2-amino-3-hydroxypropanoic acid',
   common: ['L-serine', 'L-Ser']},
  {smiles: 'SC[C@H](N)C(=O)O', iupac: '(2R)-2-amino-3-sulfanylpropanoic acid',
   common: ['L-cysteine', 'L-Cys'],
   note: 'L-시스테인은 황의 우선순위 때문에 다른 L-아미노산과 달리 (R)로 표기됩니다.'},
  {smiles: 'CC(C)[C@H](N)C(=O)O', iupac: '(2S)-2-amino-3-methylbutanoic acid',
   common: ['L-valine', 'L-Val']},
  {smiles: 'CC(C)C[C@H](N)C(=O)O', iupac: '(2S)-2-amino-4-methylpentanoic acid',
   common: ['L-leucine', 'L-Leu']},
  {smiles: 'OC(=O)[C@@H](N)Cc1ccccc1', iupac: '(2S)-2-amino-3-phenylpropanoic acid',
   common: ['L-phenylalanine', 'L-Phe']},
  {smiles: 'C[C@@H](O)[C@H](N)C(=O)O', iupac: '(2S,3R)-2-amino-3-hydroxybutanoic acid',
   common: ['L-threonine', 'L-Thr']},
  {smiles: 'OC(=O)[C@@H]1CCCN1', iupac: '(2S)-pyrrolidine-2-carboxylic acid',
   common: ['L-proline', 'L-Pro']},
  {smiles: 'C[C@H](O)C(=O)O', iupac: '(2S)-2-hydroxypropanoic acid',
   common: ['L-lactic acid', 'L-젖산']},
  {smiles: 'OC[C@@H](O)C=O', iupac: '(2R)-2,3-dihydroxypropanal',
   common: ['D-glyceraldehyde', 'D-글리세르알데하이드']},
  {smiles: 'OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C=O',
   iupac: '(2R,3S,4R,5R)-2,3,4,5,6-pentahydroxyhexanal',
   common: ['D-glucose', 'D-포도당', 'dextrose'],
   note: '사슬형(open-chain) 구조 기준입니다. 고리형(α/β-피라노스)은 따로 그려야 합니다.'},
  {smiles: 'OC[C@@H](O)[C@@H](O)[C@@H](O)C=O', iupac: '(2R,3R,4R)-2,3,4,5-tetrahydroxypentanal',
   common: ['D-ribose', 'D-리보스']},
  {smiles: 'CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O',
   iupac: '(1R,2S,5R)-5-methyl-2-(propan-2-yl)cyclohexan-1-ol',
   common: ['L-menthol', '(−)-멘톨']},
  {smiles: 'OC(=O)[C@H](O)[C@@H](O)C(=O)O', iupac: '(2R,3R)-2,3-dihydroxybutanedioic acid',
   common: ['L-(+)-tartaric acid', 'L-타타르산']},
  {smiles: 'OC(=O)[C@H](O)[C@H](O)C(=O)O', iupac: '(2R,3S)-2,3-dihydroxybutanedioic acid',
   common: ['meso-tartaric acid'],
   note: '두 입체중심이 서로 상쇄되어 광학 활성이 없는 meso 화합물입니다.'},
  {smiles: 'CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O',
   iupac: '(2S)-2-[4-(2-methylpropyl)phenyl]propanoic acid',
   common: ['dexibuprofen', '(S)-이부프로펜'],
   note: '이부프로펜의 두 거울상 중 실제 약효를 내는 쪽입니다.'},
  {smiles: 'COc1ccc2cc([C@H](C)C(=O)O)ccc2c1',
   iupac: '(2S)-2-(6-methoxynaphthalen-2-yl)propanoic acid', common: ['naproxen', '나프록센']},
  {smiles: 'CNC[C@H](O)c1ccc(O)c(O)c1',
   iupac: '4-[(1R)-1-hydroxy-2-(methylamino)ethyl]benzene-1,2-diol',
   common: ['(R)-adrenaline', '(R)-에피네프린']},
  {smiles: 'CC(=C)[C@@H]1CCC(C)=CC1', iupac: '(4R)-1-methyl-4-(prop-1-en-2-yl)cyclohex-1-ene',
   common: ['(R)-limonene', 'D-리모넨'], note: '오렌지 향이 나는 쪽 거울상입니다.'},
  {smiles: 'OC(=O)/C=C/C(=O)O', iupac: '(2E)-but-2-enedioic acid',
   common: ['fumaric acid', '푸마르산']},
  {smiles: 'OC(=O)/C=C\\C(=O)O', iupac: '(2Z)-but-2-enedioic acid',
   common: ['maleic acid', '말레산']},
  {smiles: 'CCCCCCCC/C=C\\CCCCCCCC(=O)O', iupac: '(9Z)-octadec-9-enoic acid',
   common: ['oleic acid', '올레산']},
  {smiles: 'CCCCCCCC/C=C/CCCCCCCC(=O)O', iupac: '(9E)-octadec-9-enoic acid',
   common: ['elaidic acid', '엘라이드산']},
  {smiles: 'O=C/C=C/c1ccccc1', iupac: '(2E)-3-phenylprop-2-enal',
   common: ['cinnamaldehyde', '계피 알데하이드']}
];

var index = null;

function hasStereoMarkup(smiles) {
  return /@|\/|\\/.test(smiles);
}

function buildIndex() {
  if (index) return index;
  index = {};
  COMPOUNDS.forEach(function (entry) {
    try {
      entry.mol = ONG.parseSmiles(entry.smiles);
      entry.hasStereo = hasStereoMarkup(entry.smiles);
      var key = entry.mol.structureKey();
      if (!index[key]) index[key] = [];
      index[key].push(entry);
    } catch (e) {
      if (typeof console !== 'undefined') console.warn('compound DB: ' + entry.smiles + ': ' + e.message);
    }
  });
  return index;
}

/** The engine's name for a database entry, cached. */
function entryName(entry) {
  if (entry._name === undefined) entry._name = ONG.generatedNameFor(entry.mol);
  return entry._name;
}

/**
 * Exact structure lookup. Constitution decides which entries are candidates;
 * when the drawing carries stereochemistry, the matching stereoisomer wins,
 * so L-alanine and D-alanine do not answer for each other.
 */
function lookup(mol, generatedName) {
  var idx = buildIndex();
  var bucket = idx[mol.structureKey()];
  if (!bucket) return null;
  var matches = bucket.filter(function (entry) {
    return ONG.isIsomorphic(mol, entry.mol, {limit: 1, ignoreBondOrder: true});
  });
  if (!matches.length) return null;

  if (generatedName) {
    for (var i = 0; i < matches.length; i++) {
      if (matches[i].hasStereo && entryName(matches[i]) === generatedName) return matches[i];
    }
  }
  for (var k = 0; k < matches.length; k++) {
    if (!matches[k].hasStereo) return matches[k];
  }
  // Only stereo-specific entries exist and none of them matches: the drawing
  // is a different stereoisomer of a known compound.
  return Object.assign({}, matches[0], {stereoMismatch: true});
}

ONG.COMPOUNDS = COMPOUNDS;
ONG.lookupCompound = lookup;
if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
