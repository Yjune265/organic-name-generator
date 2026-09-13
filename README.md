# Organic Name Generator

구조식을 그리면 **IUPAC 이름**을 알려 주는 웹 페이지입니다.
접두사 차이로 생기는 여러 이름 후보(iso-/sec-/tert-/neo-, 예전 locant 표기, 작용기 분류명)를 함께
보여 주고, 잘 알려진 화합물이면 **관용명**(caffeine, aspirin, 아세톤 …)도 같이 표시합니다.

서버가 필요 없습니다. 모든 계산이 브라우저 안에서 일어나므로 `index.html`을 열거나
GitHub Pages에 그대로 올리면 됩니다.

## 쓰는 법

1. `index.html`을 브라우저로 엽니다. (또는 `python3 -m http.server` 로 띄운 뒤 접속)
2. 왼쪽 [Kekule.js](https://partridgejiang.github.io/Kekule.js/) 편집기에서 구조를 그립니다.
3. **이름 생성하기**를 누릅니다.
4. SMILES 문자열이 있다면 아래 입력칸에 바로 넣어도 됩니다.

## 무엇을 보여 주나

| 항목 | 설명 |
| --- | --- |
| IUPAC 이름 | 2013 IUPAC 권고(치환 명명법)에 따라 생성 |
| 함께 쓰이는 다른 이름 | 관용 접두사(isopropyl, tert-butyl, neopentyl …), 1993년 이전 locant 표기(2-butanol), 작용기 분류명(isopropyl alcohol, diethyl ether) |
| 관용명 | 수록된 화합물 130여 종의 통용 이름 (한글 이름 포함) |
| 이름을 붙인 과정 | 모체 선택 → 주작용기 → 번호 매기기 → 치환기 순서 설명 |

## 구조

```
index.html            페이지
css/app.css           페이지 스타일 (css/kekule.css 는 편집기 스타일)
js/app.js             편집기 ↔ 명명 엔진 연결, 화면 표시
js/chem/mol.js        분자 그래프, 원자가·고리·방향족성 인식
js/chem/smiles.js     SMILES 파서 (방향족 소문자 표기 Kekule 화 포함)
js/chem/match.js      그래프 동형 사상 (구조 검색·고리 템플릿 대조)
js/chem/rings.js      이름이 있는 고리 골격 59종과 IUPAC 번호 매기기
js/chem/groups.js     작용기 인식과 서열
js/chem/naming.js     모체 선택·번호 매기기·이름 조립 (핵심)
js/chem/branches.js   치환기 접두사 이름
js/chem/variants.js   다른 이름 후보 생성
js/chem/compounds.js  잘 알려진 화합물 목록
js/chem/api.js        페이지가 부르는 진입점
js/chem/kekule-bridge.js  그린 구조 → 분자 그래프
test/                 Node 로 도는 회귀 테스트
```

## 테스트

```bash
node test/run.js
```

이름 87종 + 관용명·후보 이름 430종을 검사합니다.

## 명명 규칙 적용 범위

지원: 사슬·고리 알케인/알켄/알카인, 카복실산 · 에스터 · 아실할라이드 · 아마이드 · 나이트릴 ·
알데하이드 · 케톤 · 알코올(페놀) · 싸이올 · 아민의 서열과 접미사, ether/halide/nitro 등 접두사,
N-치환기, 벤젠 관용 모체(phenol, aniline, benzoic acid, benzaldehyde), 고리 59종
(naphthalene, pyridine, indole, quinoline, morpholine, oxolane …)의 고정 번호 매기기,
-carboxylic acid 형 접미사(시트르산 등), 알파벳 순 배열과 괄호 겹침 처리.

아직 못 하는 것:

- **입체화학** (E/Z, R/S) — 구성만 보고 이름을 만듭니다.
- **다리 걸친 고리·스피로 고리** (bicyclo/spiro) — 수록된 화합물이면 검증된 이름이 나옵니다(예: 캠퍼).
- 표에 없는 융합 고리 골격, 설폰산 등 일부 작용기.

이런 구조는 이름 대신 그 이유를 화면에 적어 줍니다. 수록된 화합물 목록에 있으면
그 검증된 이름이 우선 표시됩니다(카페인처럼 엔진이 직접 유도하지 못하는 이름도 정확히 나옵니다).
