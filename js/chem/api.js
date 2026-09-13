/**
 * api.js - what the page calls.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

var GROUP_KO = {
  carboxylicAcid: '카복실산 (-oic acid / -carboxylic acid)',
  ester: '에스터 (-oate)',
  acylHalide: '아실 할라이드 (-oyl halide)',
  amide: '아마이드 (-amide)',
  nitrile: '나이트릴 (-nitrile)',
  aldehyde: '알데하이드 (-al)',
  ketone: '케톤 (-one)',
  alcohol: '알코올 (-ol)',
  thiol: '싸이올 (-thiol)',
  amine: '아민 (-amine)',
  ether: '에터 (alkoxy- 접두사)',
  sulfide: '설파이드 (alkylsulfanyl- 접두사)',
  halide: '할로젠 (halo- 접두사)',
  nitro: '나이트로기 (nitro- 접두사)'
};

/** Human readable account of how the name was built. */
function explain(unit) {
  if (!unit || !unit.evaluation) return [];
  var ev = unit.evaluation, steps = [], numbering = ev.numbering;

  var parentDesc;
  if (ev.cand.kind === 'ring') {
    parentDesc = '고리 ' + ev.cand.atoms.length + '원자 → ' + ev.cand.info.name +
      (ev.cand.info.trivial ? ' (' + ev.cand.info.trivial + ')' : '');
  } else {
    parentDesc = '가장 긴 사슬의 탄소 ' + ev.cand.atoms.length + '개';
  }
  steps.push({title: '① 모체(parent) 선택', body: parentDesc});

  if (ev.served.length) {
    var g = ev.served[0].group;
    steps.push({
      title: '② 주작용기 결정',
      body: (GROUP_KO[g.type] || g.type) + ' — 가장 서열이 높은 작용기이므로 접미사로 붙입니다.' +
        (ev.served.length > 1 ? ' (' + ev.served.length + '개)' : '')
    });
  } else {
    steps.push({title: '② 주작용기 결정', body: '접미사로 쓸 작용기가 없어 모체 이름이 그대로 쓰입니다.'});
  }

  var locantParts = [];
  if (ev.served.length) {
    locantParts.push('주작용기 ' + numbering.suffixLocants.join(',') + '번');
  }
  if (numbering.ene.length) locantParts.push('이중결합 ' + numbering.ene.join(',') + '번');
  if (numbering.yne.length) locantParts.push('삼중결합 ' + numbering.yne.join(',') + '번');
  if (ev.branches.length) {
    locantParts.push('치환기 ' + ONG.naming.sortLocants(ev.branches.map(function (b) {
      return numbering.locantOf[b.parentAtom];
    })).join(',') + '번');
  }
  steps.push({
    title: '③ 번호 매기기',
    body: locantParts.length
      ? locantParts.join(' / ') + ' — 주작용기 → 다중결합 → 치환기 순으로 낮은 번호를 줍니다.'
      : '번호를 구분할 필요가 없습니다.'
  });

  if (ev.branches.length) {
    var listed = ev.branches.map(function (b) {
      var loc = numbering.locantOf[b.parentAtom];
      return (loc ? loc + '번 ' : '') + b.name;
    });
    steps.push({title: '④ 치환기 이름 붙이기', body: listed.join(', ') + ' — 알파벳 순으로 나열합니다.'});
  }
  return steps;
}

/**
 * Name a structure and gather everything the UI shows.
 */
function nameStructure(mol) {
  var components = mol.components();
  var warnings = [], generated = [], steps = [], parts = null;

  components.forEach(function (atoms) {
    var sub = mol.subMol(atoms);
    var ctx = new ONG.naming.Ctx(sub.mol);
    var unit = null;
    try {
      unit = ONG.naming.nameUnit(ctx, sub.mol.atoms.map(function (a, i) { return i; }), null);
    } catch (e) {
      ctx.warn('이름을 만드는 중 오류가 발생했습니다: ' + e.message);
    }
    ctx.warnings.forEach(function (w) { if (warnings.indexOf(w) < 0) warnings.push(w); });
    generated.push(unit && !ctx.warnings.length ? unit.name : (unit ? unit.name : null));
    if (components.length === 1 && unit) {
      steps = explain(unit);
      parts = unit.parts;
    }
  });

  var generatedName = generated.every(function (n) { return n; }) ? generated.join(' ; ') : null;
  var entry = ONG.lookupCompound ? ONG.lookupCompound(mol, generatedName) : null;

  var stereoInfo = ONG.stereo.analyse(mol);
  var stereoDefined = stereoInfo.any;
  if (stereoInfo.undefinedCentres.length) {
    warnings.push('입체중심(카이랄 탄소)이 ' + stereoInfo.undefinedCentres.length +
      '개 있지만 배치가 지정되어 있지 않아 R/S를 붙이지 않았습니다. ' +
      '편집기에서 쐐기(굵은 선)·점선 결합으로 표시하거나 SMILES에 @/@@를 쓰면 됩니다.');
  }
  if (stereoInfo.undefinedBonds.length) {
    warnings.push('이중결합의 기하 배치(E/Z)가 지정되어 있지 않습니다. ' +
      '구조를 그릴 때 치환기 위치를 분명히 하거나 SMILES에 /, \\ 를 쓰면 E/Z를 붙여 드립니다.');
  }

  var alternatives = [];
  if (components.length === 1) {
    try {
      alternatives = ONG.variants.alternatives(mol, generatedName, parts);
    } catch (e) { /* alternatives are a bonus, never fatal */ }
  }

  // A curated name covers skeletons the generator cannot derive, but a
  // generated name that carries stereodescriptors says more, so it wins.
  var curated = entry && entry.iupac ? entry.iupac : null;
  var preferGenerated = stereoDefined && generatedName &&
    curated && !/\([\dRSEZ,]+\)-/.test(curated);
  var primary = preferGenerated ? generatedName : (curated || generatedName);
  var source = (primary === curated) ? 'database' : (generatedName ? 'generated' : null);

  if (curated && primary !== curated) {
    alternatives.unshift({
      name: curated, kind: 'database',
      note: '입체 배치를 빼고 부르는 이름'
    });
  } else if (curated && generatedName && generatedName !== curated) {
    alternatives.unshift({
      name: generatedName, kind: 'generated',
      note: '이 프로그램이 규칙대로 만들어 낸 이름'
    });
  }
  if (entry && entry.stereoMismatch) {
    warnings.push('그린 구조는 ' + (entry.common && entry.common[0] ? entry.common[0] : '수록된 화합물') +
      '의 다른 입체 이성질체입니다. 관용명은 참고로만 보세요.');
  }
  // Only say something when the curated name really does ignore the stereo.
  var curatedHasStereo = curated && /\([\dRSEZ,]+\)-/.test(curated);
  if (stereoDefined && entry && !entry.hasStereo && !curatedHasStereo) {
    warnings.push('관용명은 입체 이성질체를 구분하지 않고 부르는 이름입니다.');
  }
  alternatives = alternatives.filter(function (a) { return a.name !== primary; });

  // A curated name already answers the question, so the generator's
  // complaints about skeletons it cannot handle are just noise.
  if (entry && entry.iupac) {
    warnings = warnings.filter(function (w) { return w.indexOf('고리 골격') < 0 && w.indexOf('치환기가 있습니다') < 0; });
  }

  return {
    formula: mol.formula(),
    mass: mol.molecularWeight(),
    atomCount: mol.atoms.length,
    name: primary,
    generatedName: generatedName,
    stereo: stereoInfo,
    source: source,
    alternatives: alternatives,
    common: entry ? (entry.common || []) : [],
    note: entry ? entry.note : null,
    steps: steps,
    warnings: warnings,
    multipart: components.length > 1
  };
}

/** Convenience wrapper for SMILES input. */
function nameSmiles(smiles) {
  return nameStructure(ONG.parseSmiles(smiles));
}

ONG.nameStructure = nameStructure;
ONG.nameSmiles = nameSmiles;
// Kept for the tests, which only care about the generated name.
ONG.nameMolecule = function (mol) {
  var r = nameStructure(mol);
  return {
    name: r.generatedName,
    formula: r.formula, mass: r.mass, warnings: r.warnings,
    components: mol.components(), full: r
  };
};

if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
