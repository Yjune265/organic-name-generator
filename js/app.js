/**
 * app.js - wiring between the Kekule editor and the naming engine.
 */
(function () {
'use strict';

var ONG = window.ONG;
var composer = null;

var EXAMPLES = [
  ['카페인', 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C'],
  ['아스피린', 'CC(=O)Oc1ccccc1C(=O)O'],
  ['이부프로펜', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'],
  ['아세톤', 'CC(C)=O'],
  ['네오펜테인', 'CC(C)(C)C'],
  ['아이소뷰틸 알코올', 'CC(C)CO'],
  ['바닐린', 'COc1cc(C=O)ccc1O'],
  ['글루코스', 'OCC(O)C(O)C(O)C(O)C=O'],
  ['시트르산', 'OC(=O)CC(O)(CC(=O)O)C(=O)O'],
  ['나프탈렌', 'c1ccc2ccccc2c1']
];

/* ------------------------------------------------------------------ *
 * start-up
 * ------------------------------------------------------------------ */

function whenKekuleReady(callback) {
  var waited = 0;
  (function poll() {
    if (window.Kekule && Kekule.Editor && Kekule.Editor.Composer && document.getElementById('editor')) {
      callback(true);
      return;
    }
    waited += 120;
    if (waited > 15000) { callback(false); return; }
    setTimeout(poll, 120);
  })();
}

function init() {
  document.getElementById('btn-name').addEventListener('click', nameFromEditor);
  document.getElementById('btn-clear').addEventListener('click', clearEditor);
  document.getElementById('btn-smiles').addEventListener('click', nameFromSmiles);
  document.getElementById('smiles-input').addEventListener('keydown', function (e) {
    if (e.key === 'Enter') nameFromSmiles();
  });

  var box = document.getElementById('examples');
  EXAMPLES.forEach(function (ex) {
    var b = document.createElement('button');
    b.textContent = ex[0];
    b.addEventListener('click', function () {
      document.getElementById('smiles-input').value = ex[1];
      loadIntoEditor(ex[1]);
      nameFromSmiles();
    });
    box.appendChild(b);
  });

  whenKekuleReady(function (ok) {
    var holder = document.getElementById('editor');
    if (!ok) {
      holder.innerHTML = '<div id="editor-fallback">구조 편집기를 불러오지 못했습니다.<br>' +
        '아래 SMILES 입력칸은 그대로 사용할 수 있습니다.</div>';
      return;
    }
    holder.innerHTML = '';
    try {
      composer = new Kekule.Editor.Composer(holder);
      if (composer.setPredefinedSetting) composer.setPredefinedSetting('fullFunc');
      if (composer.setEnableCreateNewChild) composer.setEnableCreateNewChild(true);
      if (composer.resized) composer.resized();
    } catch (e) {
      holder.innerHTML = '<div id="editor-fallback">구조 편집기 초기화에 실패했습니다: ' +
        escapeHtml(e.message) + '</div>';
      composer = null;
    }
  });
}

/* ------------------------------------------------------------------ *
 * actions
 * ------------------------------------------------------------------ */

function nameFromEditor() {
  if (!composer) {
    showError('구조 편집기가 준비되지 않았습니다. SMILES 입력칸을 사용해 주세요.');
    return;
  }
  var extracted;
  try {
    extracted = ONG.fromComposer(composer, window.Kekule);
  } catch (e) {
    showError('그린 구조를 읽지 못했습니다: ' + e.message);
    return;
  }
  if (!extracted.mol) {
    showError(extracted.warnings.length ? extracted.warnings.join(' ') : '먼저 구조를 그려 주세요.');
    return;
  }
  render(ONG.nameStructure(extracted.mol), extracted.warnings);
}

function nameFromSmiles() {
  var text = document.getElementById('smiles-input').value.trim();
  if (!text) { showError('SMILES를 입력해 주세요.'); return; }
  var mol;
  try {
    mol = ONG.parseSmiles(text);
  } catch (e) {
    showError('SMILES를 해석하지 못했습니다: ' + e.message);
    return;
  }
  render(ONG.nameStructure(mol), []);
}

function loadIntoEditor(smiles) {
  if (!composer || !window.Kekule || !Kekule.IO || !Kekule.IO.loadFormatData) return;
  try {
    var obj = Kekule.IO.loadFormatData(smiles, 'smi');
    if (obj) composer.setChemObj(obj);
  } catch (e) { /* the editor simply keeps whatever is on the canvas */ }
}

function clearEditor() {
  if (composer && composer.newDoc) {
    try { composer.newDoc(); } catch (e) { /* ignore */ }
  }
  document.getElementById('smiles-input').value = '';
  document.getElementById('result-name').className = 'result-name empty';
  document.getElementById('result-name').innerHTML =
    '구조를 그리고 <b>이름 생성하기</b>를 누르세요.';
  document.getElementById('result-body').innerHTML = '';
}

/* ------------------------------------------------------------------ *
 * rendering
 * ------------------------------------------------------------------ */

function showError(message) {
  var el = document.getElementById('result-name');
  el.className = 'result-name empty';
  el.textContent = message;
  document.getElementById('result-body').innerHTML = '';
}

function render(result, extraWarnings) {
  var nameEl = document.getElementById('result-name');
  if (result.name) {
    nameEl.className = 'result-name';
    nameEl.textContent = result.name;
  } else {
    nameEl.className = 'result-name empty';
    nameEl.textContent = '이 구조의 이름을 만들지 못했습니다. 아래 설명을 확인해 주세요.';
  }

  var html = '';
  if (result.name) {
    html += '<span class="badge">' +
      (result.source === 'database' ? '잘 알려진 화합물 — 검증된 이름' : 'IUPAC 규칙으로 생성') +
      '</span>';
  }

  html += '<div class="facts">' +
    '<div>분자식<b>' + formatFormula(result.formula) + '</b></div>' +
    '<div>분자량<b>' + result.mass.toFixed(2) + '</b></div>' +
    '<div>원자 수(H 제외)<b>' + result.atomCount + '</b></div>' +
    '</div>';

  if (result.common && result.common.length) {
    html += section('관용명 (common name)',
      '<div class="common">' + result.common.map(function (c) {
        return '<span>' + escapeHtml(c) + '</span>';
      }).join('') + '</div>');
  }

  if (result.alternatives && result.alternatives.length) {
    html += section('함께 쓰이는 다른 이름',
      '<ul class="names">' + result.alternatives.map(function (a) {
        return '<li><span class="n">' + escapeHtml(a.name) + '</span>' +
          (a.note ? '<span class="why">' + escapeHtml(a.note) + '</span>' : '') + '</li>';
      }).join('') + '</ul>');
  }

  if (result.note) {
    html += section('참고', '<p class="hint">' + escapeHtml(result.note) + '</p>');
  }

  if (result.steps && result.steps.length) {
    html += section('이름을 붙인 과정',
      '<ol class="steps">' + result.steps.map(function (s) {
        return '<li><b>' + escapeHtml(s.title) + '</b><span>' + escapeHtml(s.body) + '</span></li>';
      }).join('') + '</ol>');
  }

  var warnings = (result.warnings || []).concat(extraWarnings || []);
  if (warnings.length) {
    html += '<div class="warnings"><b>알려 드릴 점</b><ul>' +
      warnings.map(function (w) { return '<li>' + escapeHtml(w) + '</li>'; }).join('') +
      '</ul></div>';
  }

  document.getElementById('result-body').innerHTML = html;
}

function section(title, body) {
  return '<div class="section"><h3>' + escapeHtml(title) + '</h3>' + body + '</div>';
}

/** C8H10N4O2 with subscripted digits. */
function formatFormula(formula) {
  return escapeHtml(formula).replace(/(\d+)/g, '<sub>$1</sub>');
}

function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, function (c) {
    return {'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c];
  });
}

if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
else init();

})();
