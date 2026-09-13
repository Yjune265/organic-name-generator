/**
 * naming.js - substitutive IUPAC nomenclature.
 *
 * The same routine (nameUnit) names the parent hydride and, recursively,
 * every substituent: pick candidate parents, score every allowed numbering,
 * then assemble prefixes + parent + suffix.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

/* ------------------------------------------------------------------ *
 * Vocabulary
 * ------------------------------------------------------------------ */

var STEMS = ['', 'meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non', 'dec',
  'undec', 'dodec', 'tridec', 'tetradec', 'pentadec', 'hexadec', 'heptadec', 'octadec',
  'nonadec', 'icos', 'henicos', 'docos', 'tricos', 'tetracos', 'pentacos', 'hexacos',
  'heptacos', 'octacos', 'nonacos', 'triacont'];

var MULT = ['', '', 'di', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'nona', 'deca',
  'undeca', 'dodeca', 'trideca', 'tetradeca'];

var MULT_COMPLEX = ['', '', 'bis', 'tris', 'tetrakis', 'pentakis', 'hexakis', 'heptakis',
  'octakis', 'nonakis', 'decakis'];

var SUFFIX = {
  ammonium:       {inChain: 'aminium'},
  carboxylicAcid: {inChain: 'oic acid', attached: 'carboxylic acid'},
  sulfonicAcid:   {attached: 'sulfonic acid'},
  sulfinicAcid:   {attached: 'sulfinic acid'},
  sulfonateEster: {attached: 'sulfonate'},
  sulfonylHalide: {attached: 'sulfonyl'},
  sulfonamide:    {attached: 'sulfonamide'},
  thione:         {inChain: 'thione'},
  imine:          {inChain: 'imine'},
  ester:          {inChain: 'oate',     attached: 'carboxylate'},
  acylHalide:     {inChain: 'oyl',      attached: 'carbonyl'},
  amide:          {inChain: 'amide',    attached: 'carboxamide'},
  nitrile:        {inChain: 'nitrile',  attached: 'carbonitrile'},
  aldehyde:       {inChain: 'al',       attached: 'carbaldehyde'},
  ketone:         {inChain: 'one'},
  alcohol:        {inChain: 'ol'},
  thiol:          {inChain: 'thiol'},
  amine:          {inChain: 'amine'}
};

var VOWELS = 'aeiouy';

// Retained substituent prefixes that many textbooks still use. They are
// only emitted when the caller asks for the "retained" style.
var RETAINED_SUBSTITUENTS = {
  'propan-2-yl': 'isopropyl',
  '2-methylpropyl': 'isobutyl',
  'butan-2-yl': 'sec-butyl',
  '2-methylpropan-2-yl': 'tert-butyl',
  '2,2-dimethylpropyl': 'neopentyl',
  '3-methylbutyl': 'isopentyl',
  '2-methylbutan-2-yl': 'tert-pentyl',
  'prop-1-en-2-yl': 'isopropenyl',
  'ethenyl': 'vinyl',
  'ethynyl': 'ethynyl',
  'prop-2-en-1-yl': 'allyl',
  'phenylmethyl': 'benzyl',
  '2-phenylethyl': 'phenethyl',
  'diphenylmethyl': 'benzhydryl'
};

// Suffixes that can only occupy a terminal skeletal atom of a chain.
var TERMINAL_SUFFIX = {
  carboxylicAcid: true, ester: true, acylHalide: true,
  amide: true, nitrile: true, aldehyde: true
};

function stemFor(n) { return STEMS[n] || ('C' + n); }
function multiplier(n, complex) {
  var table = complex ? MULT_COMPLEX : MULT;
  if (n < table.length) return table[n];
  return n + (complex ? 'kis' : '-');
}

/* ------------------------------------------------------------------ *
 * Locant helpers
 * ------------------------------------------------------------------ */

/** '4a' -> 4.01, so fused-ring locants sort the way IUPAC reads them. */
function locantValue(label) {
  if (label === 'N') return 0.5;
  var m = /^(\d+)([a-z]*)$/.exec(String(label));
  if (!m) return 999;
  var v = parseInt(m[1], 10);
  if (m[2]) v += 0.01 * (m[2].charCodeAt(0) - 96);
  return v;
}

function sortLocants(list) {
  return list.slice().sort(function (a, b) { return locantValue(a) - locantValue(b); });
}

/** First point of difference; a shorter list wins when otherwise equal. */
function compareLocantLists(a, b) {
  var n = Math.min(a.length, b.length);
  for (var i = 0; i < n; i++) {
    var d = locantValue(a[i]) - locantValue(b[i]);
    if (Math.abs(d) > 1e-9) return d < 0 ? -1 : 1;
  }
  return a.length - b.length;
}

function compareVectors(a, b) {
  for (var i = 0; i < Math.max(a.length, b.length); i++) {
    var x = a[i], y = b[i];
    if (x === undefined) return -1;
    if (y === undefined) return 1;
    var c;
    if (Array.isArray(x)) c = compareLocantLists(x, y);
    else if (typeof x === 'string') c = x < y ? -1 : (x > y ? 1 : 0);
    else c = x - y;
    if (c) return c < 0 ? -1 : 1;
  }
  return 0;
}

/* ------------------------------------------------------------------ *
 * Naming context
 * ------------------------------------------------------------------ */

function Ctx(mol, style) {
  this.mol = mol;
  this.style = style || {};
  this.stereoInfo = (ONG.stereo && !(this.style || {}).noStereo)
    ? ONG.stereo.analyse(mol) : {atoms: {}, bonds: {}, undefinedCentres: [], undefinedBonds: []};
  this.groups = ONG.perceiveGroups(mol);
  this.warnings = [];
  this.steps = [];
}

Ctx.prototype.warn = function (message) {
  if (this.warnings.indexOf(message) < 0) this.warnings.push(message);
};

/* ------------------------------------------------------------------ *
 * Candidate parents
 * ------------------------------------------------------------------ */

/** Every maximal carbon chain inside `allowed` (ring atoms excluded). */
function chainCandidates(ctx, allowedSet, attachment, excluded) {
  var mol = ctx.mol, pool = [];
  for (var i = 0; i < mol.atoms.length; i++) {
    if (!allowedSet[i]) continue;
    if (mol.atoms[i].element !== 'C') continue;
    if (mol.isInRing(i)) continue;
    if (excluded && excluded[i]) continue;
    pool.push(i);
  }
  if (!pool.length) return [];
  var inPool = {};
  pool.forEach(function (i) { inPool[i] = true; });

  function poolNeighbors(i) {
    return mol.heavyNeighbors(i).filter(function (j) { return inPool[j]; });
  }

  var chains = [], seen = {};
  function addChain(path) {
    var key = path.join(',');
    var rkey = path.slice().reverse().join(',');
    if (seen[key] || seen[rkey]) return;
    seen[key] = seen[rkey] = true;
    chains.push(path);
  }

  var wantsAttachment = (attachment !== null && attachment !== undefined);
  if (wantsAttachment && !inPool[attachment]) return [];

  // The pool is a forest, so the path between two chain ends is unique.
  // A substituent's parent is the longest chain *containing* the attachment
  // atom - the free valence then takes the lowest locant it can
  // (propan-2-yl, not 1-methylethyl).
  var ends = pool.filter(function (i) { return poolNeighbors(i).length <= 1; });
  if (!ends.length) ends = pool.slice();
  for (var x = 0; x < ends.length; x++) {
    for (var y = x; y < ends.length; y++) {
      var path = pathBetween(ends[x], ends[y]);
      if (!path) continue;
      if (wantsAttachment && path.indexOf(attachment) < 0) continue;
      addChain(path);
    }
  }
  return chains;

  function pathBetween(from, to) {
    if (from === to) return [from];
    var prev = {}, queue = [from], seenA = {};
    seenA[from] = true;
    while (queue.length) {
      var v = queue.shift();
      var nb = poolNeighbors(v);
      for (var k = 0; k < nb.length; k++) {
        var w = nb[k];
        if (seenA[w]) continue;
        seenA[w] = true; prev[w] = v;
        if (w === to) {
          var out = [to], cur = to;
          while (cur !== from) { cur = prev[cur]; out.push(cur); }
          return out.reverse();
        }
        queue.push(w);
      }
    }
    return null;
  }
}

function ringCandidates(ctx, allowedSet, attachment) {
  var mol = ctx.mol, out = [];
  mol.ringSystems().forEach(function (system) {
    if (!system.atoms.every(function (a) { return allowedSet[a]; })) return;
    if (attachment !== null && attachment !== undefined && system.atoms.indexOf(attachment) < 0) return;
    var info = ONG.rings.describeRingSystem(mol, system);
    if (!info) {
      ctx.warn('이 프로그램이 아직 이름을 붙이지 못하는 고리 골격이 있습니다 (다리걸친 고리·스피로 고리 등).');
      return;
    }
    info.system = system;
    out.push(info);
  });
  return out;
}

/* ------------------------------------------------------------------ *
 * nameUnit - the core
 * ------------------------------------------------------------------ */

/**
 * Name one parent unit.
 *  allowed     : array of atom indices this unit may use
 *  attachment  : atom bonded to the outside (substituent), or null for the
 *                molecule's own parent
 * Returns {name, sortName, complex, locantOfAttachment, ...} or null.
 */
function nameUnit(ctx, allowed, attachment) {
  var mol = ctx.mol;
  var allowedSet = {};
  allowed.forEach(function (i) { allowedSet[i] = true; });
  var isSubstituent = (attachment !== null && attachment !== undefined);

  var principal = isSubstituent ? [] : ctx.groups.principal().filter(function (g) {
    return allowed.indexOf(anchorOf(ctx, g)) >= 0 || true;
  });

  var candidates = [];
  function addChains(excluded) {
    chainCandidates(ctx, allowedSet, isSubstituent ? attachment : null, excluded).forEach(function (path) {
      candidates.push({kind: 'chain', atoms: path, numberings: [path, path.slice().reverse()]});
    });
  }
  addChains(null);
  // Groups such as -COOH may also be expressed as "-carboxylic acid" on a
  // shorter parent (citric acid), so offer chains that leave them out.
  var carboCarbons = {};
  principal.forEach(function (g) {
    var spec = SUFFIX[g.type];
    if (spec && spec.attached && g.carbon !== undefined) carboCarbons[g.carbon] = true;
  });
  if (Object.keys(carboCarbons).length) addChains(carboCarbons);
  ringCandidates(ctx, allowedSet, isSubstituent ? attachment : null).forEach(function (info) {
    candidates.push({
      kind: 'ring', ringKind: info.kind, atoms: info.atoms, info: info,
      numberings: info.numberings
    });
  });
  if (!candidates.length) {
    if (!isSubstituent) {
      var hasCarbon = allowed.some(function (i) { return mol.atoms[i].element === 'C'; });
      ctx.warn(hasCarbon
        ? '이 골격에서는 모체(parent)를 정하지 못했습니다.'
        : '탄소가 없어 유기 화합물 명명법의 대상이 아닙니다 (물, 암모니아, 무기 이온 등).');
    }
    return null;
  }

  // Candidates that lose should not leave their complaints behind: keep the
  // warnings raised while naming the parent we actually chose.
  var best = null, before = ctx.warnings.slice(), everything = before.slice();
  candidates.forEach(function (cand) {
    ctx.warnings = before.slice();
    var evaluated = evaluateCandidate(ctx, cand, allowed, allowedSet, principal, attachment);
    ctx.warnings.forEach(function (w) { if (everything.indexOf(w) < 0) everything.push(w); });
    if (!evaluated) return;
    evaluated.warnings = ctx.warnings.slice();
    if (!best || compareVectors(evaluated.score, best.score) < 0) best = evaluated;
  });
  ctx.warnings = best ? best.warnings : everything;
  if (!best) return null;
  return assemble(ctx, best, isSubstituent);
}

/** Atom whose locant a group is cited at, given the parent it belongs to. */
function anchorOf(ctx, group) {
  switch (group.type) {
    case 'alcohol': return group.carrier;
    case 'thiol': return group.carrier;
    case 'amine': return group.nitrogen;
    case 'ammonium': return group.nitrogen;
    // Only the sulfur oxo-acids are anchored on sulfur; a thione is not.
    default: return (group.sulfur !== undefined && group.carrier !== undefined)
      ? group.sulfur : group.carbon;
  }
}

function evaluateCandidate(ctx, cand, allowed, allowedSet, principal, attachment) {
  var mol = ctx.mol;
  var inParent = {};
  cand.atoms.forEach(function (i) { inParent[i] = true; });

  /* -- which principal groups does this parent express? --------------- */
  var servedInChain = [], servedAttached = [];
  principal.forEach(function (g) {
    var placement = placeGroup(ctx, g, inParent);
    if (!placement) return;
    (placement.mode === 'inChain' ? servedInChain : servedAttached).push({group: g, placement: placement});
  });
  // A name cites one suffix form, so the two modes cannot be mixed.
  var served = servedInChain.length >= servedAttached.length ? servedInChain : servedAttached;
  if (principal.length && !served.length && cand.kind === 'chain') {
    // A chain that expresses none of the principal groups is never the parent.
    return null;
  }

  /* -- branches hanging off the parent -------------------------------- */
  var branches = collectBranches(ctx, cand.atoms, inParent, allowedSet, served);
  if (branches === null) return null;

  /* -- skeletal unsaturation ------------------------------------------ */
  var isCarbocycle = cand.kind === 'ring' && (cand.ringKind === 'carbocycle');
  var isChain = cand.kind === 'chain';

  var bestNumbering = null;
  cand.numberings.forEach(function (numbering) {
    var locantOf = {};
    numbering.forEach(function (item, pos) {
      if (typeof item === 'number') locantOf[item] = String(pos + 1);
      else locantOf[item.atom] = item.label;
    });
    // Substituents on the nitrogen of a principal amine/amide are cited as N-.
    served.forEach(function (s) {
      if (s.group.nitrogen !== undefined) locantOf[s.group.nitrogen] = 'N';
    });
    var atomsInOrder = numbering.map(function (item) {
      return typeof item === 'number' ? item : item.atom;
    });

    var suffixLocants = served.map(function (s) {
      return locantOf[s.placement.locantAtom];
    }).filter(function (x) { return x !== undefined; });
    suffixLocants = sortLocants(suffixLocants);

    var ene = [], yne = [];
    if (isChain || isCarbocycle) {
      var n = atomsInOrder.length;
      var limit = isChain ? n - 1 : n;
      for (var k = 0; k < limit; k++) {
        var a = atomsInOrder[k], b = atomsInOrder[(k + 1) % n];
        var order = mol.bondOrder(a, b);
        if (order === 2) ene.push(String(k + 1));
        else if (order === 3) yne.push(String(k + 1));
      }
    }

    var prefixLocants = sortLocants(branches.map(function (br) { return locantOf[br.parentAtom]; }));
    var attachmentLocant = attachment !== undefined && attachment !== null ? locantOf[attachment] : null;

    var score = [];
    if (attachmentLocant !== null) score.push(locantValue(attachmentLocant));
    score.push(suffixLocants, sortLocants(ene.concat(yne)), ene, prefixLocants);

    var entry = {
      numbering: numbering, locantOf: locantOf, atomsInOrder: atomsInOrder,
      suffixLocants: suffixLocants, ene: ene, yne: yne,
      prefixLocants: prefixLocants, attachmentLocant: attachmentLocant, score: score
    };
    if (!bestNumbering || compareVectors(entry.score, bestNumbering.score) < 0) bestNumbering = entry;
  });
  if (!bestNumbering) return null;

  // Are all skeletal positions equivalent? (benzene: yes, naphthalene: no)
  var firstPositions = {};
  cand.numberings.forEach(function (numbering) {
    var head = numbering[0];
    firstPositions[typeof head === 'number' ? head : head.atom] = true;
  });
  var allEquivalent = cand.atoms.every(function (a) { return firstPositions[a]; });

  var nRings = cand.kind === 'ring' ? cand.info.system.rings.length : 0;
  var nMultiple = bestNumbering.ene.length + bestNumbering.yne.length;

  var score = [
    -served.length,
    cand.kind === 'ring' ? 0 : 1,
    -nRings,
    -cand.atoms.length,
    -nMultiple,
    -bestNumbering.ene.length,
    bestNumbering.suffixLocants,
    sortLocants(bestNumbering.ene.concat(bestNumbering.yne)),
    -branches.length,
    bestNumbering.prefixLocants,
    branches.map(function (b) { return b.sortName; }).sort().join(',')
  ];

  return {
    cand: cand, served: served, branches: branches, numbering: bestNumbering,
    allEquivalent: allEquivalent, parentSize: cand.atoms.length,
    // N-substituents say nothing about where the suffix sits on the skeleton,
    // so they do not force skeletal locants to be cited.
    citations: served.length + (attachment === null || attachment === undefined ? 0 : 1) +
      branches.filter(function (b) {
        return bestNumbering.locantOf[b.parentAtom] !== 'N';
      }).length,
    isChain: isChain, isCarbocycle: isCarbocycle, score: score, attachment: attachment
  };
}

/**
 * Where a characteristic group sits relative to a candidate parent.
 * Returns null when the parent cannot express it.
 */
function placeGroup(ctx, group, inParent) {
  var mol = ctx.mol;
  var anchor = anchorOf(ctx, group);
  if (group.type === 'alcohol' || group.type === 'thiol') {
    return inParent[group.carrier] ? {mode: 'inChain', locantAtom: group.carrier} : null;
  }
  // -SO3H, -SO2NH2 and friends hang off a carbon: that carbon takes the locant.
  if (group.sulfur !== undefined && group.carrier !== undefined) {
    return inParent[group.carrier] ? {mode: 'attached', locantAtom: group.carrier} : null;
  }
  if (group.type === 'amine' || group.type === 'ammonium') {
    var carriers = mol.heavyNeighbors(group.nitrogen).filter(function (j) { return inParent[j]; });
    return carriers.length ? {mode: 'inChain', locantAtom: carriers[0]} : null;
  }
  if (inParent[anchor]) return {mode: 'inChain', locantAtom: anchor};
  // suffix carbon outside the parent: -carboxylic acid / -carbaldehyde style
  if (!SUFFIX[group.type] || !SUFFIX[group.type].attached) return null;
  var links = mol.heavyNeighbors(anchor).filter(function (j) { return inParent[j]; });
  return links.length ? {mode: 'attached', locantAtom: links[0]} : null;
}

/** Everything hanging off the parent, already named as prefixes. */
function collectBranches(ctx, parentAtoms, inParent, allowedSet, served) {
  var mol = ctx.mol;
  var servedAtoms = {};
  served.forEach(function (s) {
    (s.group.atoms || []).forEach(function (a) { if (!inParent[a]) servedAtoms[a] = true; });
    if (s.placement.mode === 'attached') servedAtoms[s.group.carbon] = true;
  });

  var visited = {}, branches = [];
  var roots = parentAtoms.slice();
  // A principal amine / amide nitrogen carries its own substituents (N-...).
  served.forEach(function (s) {
    if (s.group.nitrogen !== undefined && !inParent[s.group.nitrogen]) {
      roots.push(s.group.nitrogen);
      delete servedAtoms[s.group.nitrogen];
      servedAtoms[s.group.nitrogen] = 'root';
    }
  });
  var rootSet = {};
  roots.forEach(function (r) { rootSet[r] = true; });

  for (var p = 0; p < roots.length; p++) {
    var atom = roots[p];
    var links = mol.links(atom);
    for (var k = 0; k < links.length; k++) {
      var nb = links[k].atom;
      if (inParent[nb] || rootSet[nb] || !allowedSet[nb]) continue;
      if (mol.atoms[nb].element === 'H') continue;
      if (servedAtoms[nb] || visited[nb]) continue;
      var blocked = {};
      Object.keys(inParent).forEach(function (key) { blocked[key] = true; });
      roots.forEach(function (r) { blocked[r] = true; });
      var atoms = branchAtoms(mol, nb, blocked, allowedSet);
      atoms.forEach(function (a) { visited[a] = true; });
      // A branch touching the parent twice would need bridged nomenclature.
      var contacts = 0;
      atoms.forEach(function (a) {
        mol.heavyNeighbors(a).forEach(function (j) { if (inParent[j] || rootSet[j]) contacts++; });
      });
      if (contacts > 1) {
        ctx.warn('다리 걸친(bridged) 구조는 아직 지원하지 않습니다.');
        return null;
      }
      var named = ONG.nameBranch(ctx, atoms, nb, atom, links[k].bond.order);
      if (!named) return null;
      named.parentAtom = atom;
      branches.push(named);
    }
  }
  return branches;
}

function branchAtoms(mol, start, inParent, allowedSet) {
  var stack = [start], seen = {}, out = [];
  seen[start] = true;
  while (stack.length) {
    var v = stack.pop();
    out.push(v);
    mol.heavyNeighbors(v).forEach(function (j) {
      if (seen[j] || inParent[j] || !allowedSet[j]) return;
      seen[j] = true; stack.push(j);
    });
  }
  return out;
}

/* ------------------------------------------------------------------ *
 * Assembly
 * ------------------------------------------------------------------ */

/** Drop the parent's final "e" when the suffix starts with a vowel. */
function elide(stem, suffixCore) {
  var letter = (suffixCore || '')[0];
  if (VOWELS.indexOf(letter) >= 0 && stem[stem.length - 1] === 'e') return stem.slice(0, -1);
  return stem;
}

/** parent hydride: stem + ene/yne endings */
function parentHydride(ev) {
  var cand = ev.cand, numbering = ev.numbering;
  // "ethene", "cyclohexene": one multiple bond, nothing else cited, all
  // skeletal positions alike - then the locant carries no information.
  var nUnsat = numbering.ene.length + numbering.yne.length;
  var omit = nUnsat === 1 &&
    (ev.parentSize <= 2 || (ev.citations === 0 && cand.kind === 'ring' && ev.allEquivalent));
  if (cand.kind === 'ring') {
    if (cand.ringKind === 'benzene') return 'benzene';
    if (cand.ringKind === 'template') return cand.info.name;
    return unsaturate('cyclo' + stemFor(cand.atoms.length), numbering.ene, numbering.yne, omit);
  }
  return unsaturate(stemFor(cand.atoms.length), numbering.ene, numbering.yne, omit);
}

function unsaturate(stem, ene, yne, omitLocants) {
  if (!ene.length && !yne.length) return stem + 'ane';
  var s = stem;
  if (ene.length) {
    s += (ene.length > 1 ? 'a' : '') + (omitLocants && ene.length === 1 ? '' : '-' + ene.join(',') + '-') +
         multiplier(ene.length) + 'ene';
  }
  if (yne.length) {
    if (ene.length) s = s.replace(/e$/, '');
    else s += (yne.length > 1 ? 'a' : '');
    s += (omitLocants && yne.length === 1 && !ene.length ? '' : '-' + yne.join(',') + '-') +
         multiplier(yne.length) + 'yne';
  }
  return s;
}

/** Retained names for a benzene parent carrying one senior group. */
var BENZENE_RETAINED = {
  'alcohol/inChain': 'phenol',
  'amine/inChain': 'aniline',
  'carboxylicAcid/attached': 'benzoic acid',
  'aldehyde/attached': 'benzaldehyde',
  'ester/attached': 'benzoate',
  'amide/attached': 'benzamide',
  'nitrile/attached': 'benzonitrile',
  'acylHalide/attached': 'benzoyl'
};

function assemble(ctx, ev, isSubstituent) {
  var mol = ctx.mol;
  var numbering = ev.numbering;
  var parts = [];

  /* --- prefixes ---------------------------------------------------- */
  var prefixes = ev.branches.map(function (br) {
    return {
      name: br.name, sortName: br.sortName, complex: br.complex,
      locant: numbering.locantOf[br.parentAtom]
    };
  });
  var merged = {};
  prefixes.forEach(function (p) {
    var key = p.name;
    if (!merged[key]) merged[key] = {name: p.name, sortName: p.sortName, complex: p.complex, locants: []};
    merged[key].locants.push(p.locant);
  });
  // Locants carry no information on a one-atom parent, or when a lone
  // citation sits on a skeleton whose positions are all equivalent.
  var omitLocants = ev.parentSize === 1 ||
    (ev.citations === 1 && (ev.parentSize <= 2 || ev.allEquivalent));

  var prefixList = Object.keys(merged).map(function (k) { return merged[k]; });
  prefixList.sort(function (a, b) {
    if (a.sortName === b.sortName) return 0;
    return a.sortName < b.sortName ? -1 : 1;
  });
  var prefixText = prefixList.map(function (p) {
    var locants = sortLocants(p.locants);
    var name = p.name;
    var needsMarks = p.complex || (/[\s\d-]/.test(name) && !/^(tert|sec|iso|neo)-?/.test(name));
    if (needsMarks) name = enclose(name);
    var mult = locants.length > 1 ? multiplier(locants.length, needsMarks) : '';
    var isN = locants.some(function (l) { return l === 'N'; });
    var show = locants[0] !== undefined && (isN || !omitLocants);
    return (show ? locants.join(',') + '-' : '') + mult + name;
  }).reduce(function (acc, piece) {
    if (!acc) return piece;
    return acc + (/^[\dN]/.test(piece) ? '-' : '') + piece;
  }, '');

  /* --- substituent: parent + yl ------------------------------------ */
  if (isSubstituent) {
    var hydride = parentHydride(ev);
    var loc = numbering.attachmentLocant;
    var name;
    if (ev.cand.kind === 'ring' && ev.cand.ringKind === 'benzene') {
      name = 'phenyl';
    } else if (ev.cand.kind === 'ring') {
      var base = hydride.replace(/e$/, '');
      name = (ev.isCarbocycle && loc === '1') ? base + 'yl' : base + '-' + loc + '-yl';
    } else {
      var hasUnsat = numbering.ene.length || numbering.yne.length;
      if (loc === '1' && !hasUnsat) name = hydride.replace(/ane$/, 'yl');
      else if (loc === '1' && ev.parentSize <= 2) name = hydride.replace(/e$/, '') + 'yl';
      else name = hydride.replace(/e$/, '') + '-' + loc + '-yl';
    }
    var full = joinPrefix(prefixText, name);
    var subStereo = stereoDescriptors(ctx, ev, omitLocants);
    if (subStereo) {
      return {
        name: subStereo + full, sortName: sortKeyOf(full),
        complex: true, evaluation: ev
      };
    }
    if (ctx.style && ctx.style.retained && RETAINED_SUBSTITUENTS[full]) {
      full = RETAINED_SUBSTITUENTS[full];
      return {name: full, sortName: sortKeyOf(full), complex: false, evaluation: ev};
    }
    return {
      name: full,
      sortName: sortKeyOf(full),
      complex: !!prefixText || /[-\d]/.test(name),
      evaluation: ev
    };
  }

  /* --- suffix ------------------------------------------------------- */
  var hydride = parentHydride(ev);
  var served = ev.served;
  var type = served.length ? served[0].group.type : null;
  var mode = served.length ? served[0].placement.mode : null;
  var suffixText = '', esterAlkyl = null, halideWord = null, tail = '';

  if (served.length) {
    var spec = SUFFIX[type];
    // Classes with no suffix form of their own (anhydrides) are named elsewhere.
    if (!spec) return null;
    var word = mode === 'attached' ? spec.attached : spec.inChain;
    if (!word) return null;
    var locants = sortLocants(served.map(function (s) {
      return numbering.locantOf[s.placement.locantAtom];
    }));
    // -oic acid / -al / -nitrile and friends can only sit at a chain end,
    // so their locants are never cited.
    var terminalOnChain = TERMINAL_SUFFIX[type] && mode === 'inChain' && ev.isChain;
    var omit = terminalOnChain || ev.parentSize === 1 ||
      (ev.citations === 1 && (ev.parentSize <= 2 || ev.allEquivalent));
    if (mode === 'attached' && ev.cand.kind === 'ring' && ev.cand.ringKind === 'benzene' &&
        served.length === 1 && BENZENE_RETAINED[type + '/attached']) {
      hydride = BENZENE_RETAINED[type + '/attached'];
      word = null;
    } else if (mode === 'inChain' && ev.cand.kind === 'ring' && ev.cand.ringKind === 'benzene' &&
        served.length === 1 && BENZENE_RETAINED[type + '/inChain']) {
      hydride = BENZENE_RETAINED[type + '/inChain'];
      word = null;
    }
    if (word) {
      var core = multiplier(served.length) + word;
      hydride = elide(hydride, core) + (omit ? '' : '-' + locants.join(',') + '-') + core;
    }
    if (type === 'ester' || type === 'sulfonateEster') {
      var alkylNames = served.map(function (s) {
        var anchorAtom = s.group.carbon !== undefined ? s.group.carbon : s.group.sulfur;
        var alkylStart = mol.heavyNeighbors(s.group.esterO).filter(function (j) {
          return j !== anchorAtom;
        })[0];
        var sub = ONG.nameBranch(ctx, branchAtoms(mol, alkylStart, indexSet(ev.cand.atoms.concat([s.group.esterO])), allSet(mol)), alkylStart, s.group.esterO, 1);
        return sub ? sub.name : '?';
      });
      esterAlkyl = uniqueJoin(alkylNames);
    }
    if (type === 'acylHalide' || type === 'sulfonylHalide') {
      halideWord = served.map(function (s) { return s.group.halide; })[0];
      tail = ' ' + halideWord;
    }
  }

  var full = joinPrefix(prefixText, hydride) + tail;
  var stereo = stereoDescriptors(ctx, ev, omitLocants && !ev.served.length);
  if (stereo) full = stereo + full;
  // "methyl (2R)-2-hydroxypropanoate": the descriptor belongs to the acid part.
  if (esterAlkyl) full = esterAlkyl + ' ' + full;

  return {
    name: full,
    sortName: sortKeyOf(full),
    complex: true,
    evaluation: ev,
    // Keeping the pieces lets variants.js rebuild the name in other styles
    // without having to parse it back apart.
    parts: {prefix: prefixText, parent: hydride, tail: tail, esterAlkyl: esterAlkyl}
  };
}

/**
 * Stereodescriptors for one parent unit: "(2R,3S)", "(E)" ...
 * Locants are cited unless the name omits them everywhere else.
 */
function stereoDescriptors(ctx, ev, omitLocants) {
  var info = ctx.stereoInfo;
  if (!info) return '';
  var mol = ctx.mol, numbering = ev.numbering, items = [];

  ev.cand.atoms.forEach(function (atom) {
    var d = info.atoms[atom];
    if (d) items.push({locant: numbering.locantOf[atom], descriptor: d});
  });

  var atomsInOrder = numbering.atomsInOrder;
  var n = atomsInOrder.length;
  var limit = ev.isChain ? n - 1 : n;
  for (var k = 0; k < limit && n > 1; k++) {
    var a = atomsInOrder[k], b = atomsInOrder[(k + 1) % n];
    var bond = mol.bondBetween(a, b);
    if (!bond || bond.order !== 2) continue;
    var e = info.bonds[bond.idx];
    if (e) items.push({locant: String(k + 1), descriptor: e});
  }
  if (!items.length) return '';

  items.sort(function (x, y) { return locantValue(x.locant) - locantValue(y.locant); });
  var text = items.map(function (it) {
    return (omitLocants || it.locant === undefined ? '' : it.locant) + it.descriptor;
  }).join(',');
  return '(' + text + ')-';
}

/** Prefixes and parent need a hyphen when the parent starts with a locant. */
function joinPrefix(prefixText, rest) {
  if (!prefixText) return rest;
  return prefixText + (/^[\dN]/.test(rest) ? '-' : '') + rest;
}

/** Enclosing marks nest as ( ) -> [ ] -> { }. */
function enclose(name) {
  if (name.indexOf('[') >= 0 || name.indexOf('{') >= 0) return '{' + name + '}';
  if (name.indexOf('(') >= 0) return '[' + name + ']';
  return '(' + name + ')';
}

function indexSet(list) {
  var s = {};
  list.forEach(function (i) { s[i] = true; });
  return s;
}

function allSet(mol) {
  var s = {};
  for (var i = 0; i < mol.atoms.length; i++) s[i] = true;
  return s;
}

function uniqueJoin(names) {
  var counts = {};
  names.forEach(function (n) { counts[n] = (counts[n] || 0) + 1; });
  return Object.keys(counts).map(function (n) {
    return (counts[n] > 1 ? multiplier(counts[n], /[-\d]/.test(n)) : '') + n;
  }).join(' ');
}

/** Alphabetisation key: ignore locants, enclosing marks and italic prefixes. */
function sortKeyOf(name) {
  return String(name)
    .replace(/\(|\)|\[|\]/g, '')
    .replace(/^[\d,'\-]+/, '')
    .replace(/\b(tert|sec|N|O|o|m|p)-/g, '')
    .replace(/[\d,\-]/g, '')
    .toLowerCase();
}

/**
 * R-CO-O-CO-R' is named from its two acid halves: "ethanoic anhydride",
 * "ethanoic propanoic anhydride".
 */
function anhydrideName(mol) {
  var groups = ONG.perceiveGroups(mol);
  var anhydrides = groups.ofType('anhydride');
  if (anhydrides.length !== 1) return null;
  if (groups.principalRank !== anhydrides[0].rank) return null;
  var g = anhydrides[0];
  var blocked = {};
  blocked[g.bridge] = true;

  var halves = [g.carbon, g.partner].map(function (acyl) {
    var atoms = branchAtoms(mol, acyl, blocked, allSet(mol));
    var sub = mol.subMol(atoms);
    sub.mol.addBond(sub.map[acyl], sub.mol.addAtom('O'), 1);   // close it into an acid
    var name = generatedNameFor(sub.mol);
    return name ? name.replace(/\s*acid$/, '').trim() : null;
  });
  if (halves.some(function (h) { return !h; })) return null;
  if (halves[0] === halves[1]) return halves[0] + ' anhydride';
  return halves.sort().join(' ') + ' anhydride';
}

/** The engine's own name for a molecule, or null. Never throws. */
function generatedNameFor(mol, style) {
  try {
    var ctx = new Ctx(mol, style);
    var unit = nameUnit(ctx, mol.atoms.map(function (a, i) { return i; }), null);
    if (unit && !ctx.warnings.length) return unit.name;
  } catch (e) { /* fall through to the special cases */ }
  try {
    return anhydrideName(mol);
  } catch (e) {
    return null;
  }
}

ONG.generatedNameFor = generatedNameFor;
ONG.anhydrideName = anhydrideName;

ONG.naming = {
  STEMS: STEMS, MULT: MULT, MULT_COMPLEX: MULT_COMPLEX, SUFFIX: SUFFIX,
  stemFor: stemFor, multiplier: multiplier, locantValue: locantValue,
  sortLocants: sortLocants, compareLocantLists: compareLocantLists,
  compareVectors: compareVectors, Ctx: Ctx, nameUnit: nameUnit,
  RETAINED_SUBSTITUENTS: RETAINED_SUBSTITUENTS,
  branchAtoms: branchAtoms, sortKeyOf: sortKeyOf, elide: elide, enclose: enclose,
  joinPrefix: joinPrefix, stereoDescriptors: stereoDescriptors,
  anchorOf: anchorOf, unsaturate: unsaturate, allSet: allSet
};

if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
