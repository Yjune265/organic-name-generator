/**
 * branches.js - naming everything that hangs off a parent as a prefix.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

var CONTRACTED_ALKOXY = {
  methyl: 'methoxy', ethyl: 'ethoxy', propyl: 'propoxy', butyl: 'butoxy', phenyl: 'phenoxy'
};

var RETAINED_ACYL = {methanoyl: 'formyl', ethanoyl: 'acetyl', benzenecarbonyl: 'benzoyl'};

/**
 * Name one branch as a substituent prefix.
 *  atoms      : the branch's atoms
 *  start      : the branch atom bonded to the parent
 *  parentAtom : the parent atom it is bonded to
 *  bondOrder  : order of that bond
 */
function nameBranch(ctx, atoms, start, parentAtom, bondOrder) {
  var mol = ctx.mol;
  var group = ctx.groups.groupOwning(start);
  var el = mol.atoms[start].element;

  /* --- doubly bonded single atoms: =O, =S, =NH ---------------------- */
  if (bondOrder === 2 && atoms.length === 1) {
    if (el === 'O') return plain('oxo');
    if (el === 'S') return plain('sulfanylidene');
    if (el === 'N') return plain('imino');
  }
  if (bondOrder === 2 && el === 'N') {
    var rest = restOf(atoms, [start]);
    var sub = rest.length ? nameSub(ctx, atoms, start, mol.heavyNeighbors(start).filter(function (j) {
      return j !== parentAtom;
    })[0]) : null;
    return sub ? complexName(sub.name + 'imino') : plain('imino');
  }

  if (group) {
    // Groups that only ever appear as a prefix carry their own word.
    if (group.prefix && group.atoms.indexOf(start) === 0) return plain(group.prefix);
    switch (group.type) {
      case 'halide':
        return plain(group.prefix);
      case 'sulfonicAcid':
        return plain('sulfo');
      case 'sulfonamide':
        return plain('sulfamoyl');
      case 'nitro':
        return plain('nitro');
      case 'alcohol':
        return plain('hydroxy');
      case 'thiol':
        return plain('sulfanyl');
      case 'nitrile':
        if (atoms.length === 2) return plain('cyano');
        break;
      case 'carboxylicAcid':
        if (atoms.length === 3) return plain('carboxy');
        break;
      case 'aldehyde':
        if (atoms.length === 2) return plain('formyl');
        break;
      case 'ether': {
        var other = mol.heavyNeighbors(start).filter(function (j) { return j !== parentAtom; })[0];
        if (other === undefined) break;
        var named = nameSub(ctx, restOf(atoms, [start]), other, other);
        if (!named) break;
        if (CONTRACTED_ALKOXY[named.name]) return plain(CONTRACTED_ALKOXY[named.name]);
        return complexName(wrap(named) + 'oxy');
      }
      case 'sulfide': {
        var otherS = mol.heavyNeighbors(start).filter(function (j) { return j !== parentAtom; })[0];
        if (otherS === undefined) break;
        var namedS = nameSub(ctx, restOf(atoms, [start]), otherS, otherS);
        if (!namedS) break;
        return complexName(wrap(namedS) + 'sulfanyl');
      }
      case 'amine': {
        var subs = nSubstituents(ctx, atoms, start, parentAtom);
        if (!subs.length) return plain('amino');
        return complexName(joinMultiplied(subs) + 'amino');
      }
      case 'amide': {
        // parent is on the carbonyl side: -C(=O)NR2
        if (start === group.carbon) {
          var nSubs = nSubstituents(ctx, atoms, group.nitrogen, group.carbon);
          if (!nSubs.length) return plain('carbamoyl');
          return complexName(nSubs.map(function (n) { return 'N-' + n; }).join('') + 'carbamoyl');
        }
        break;
      }
      case 'ester': {
        if (start === group.carbon) {           // parent bears the acyl carbon
          var alkyl = mol.heavyNeighbors(group.esterO).filter(function (j) { return j !== group.carbon; })[0];
          var alkylName = alkyl === undefined ? null : nameSub(ctx, restOf(atoms, [group.carbon, group.esterO]), alkyl, alkyl);
          if (!alkylName) break;
          var alkoxy = CONTRACTED_ALKOXY[alkylName.name] || (wrap(alkylName) + 'oxy');
          return complexName(alkoxy + 'carbonyl');
        }
        if (start === group.esterO) {           // parent bears the ester oxygen
          var acyl = acylName(ctx, atoms, group);
          if (acyl) return complexName(acyl + 'oxy');
        }
        break;
      }
      case 'ketone':
        break;                                   // handled inside the branch chain
      default:
        break;
    }
  }

  /* --- sulfoxide / sulfone bridges ---------------------------------- */
  if (el === 'S') {
    var oxo = mol.heavyNeighbors(start).filter(function (j) {
      return mol.atoms[j].element === 'O' && mol.bondOrder(start, j) === 2;
    });
    var carbon = mol.heavyNeighbors(start).filter(function (j) {
      return j !== parentAtom && oxo.indexOf(j) < 0;
    })[0];
    if (oxo.length && carbon !== undefined) {
      // "(methanesulfinyl)methane", "benzenesulfonyl": these prefixes are
      // built on the parent hydride name, not on the -yl substituent name.
      var word = oxo.length === 1 ? 'sulfinyl' : 'sulfonyl';
      var side = nameSub(ctx, restOf(atoms, [start].concat(oxo)), carbon, carbon);
      var hydride = side ? hydrideNameOf(side) : null;
      if (hydride) return complexName(hydride + word);
    }
  }

  /* --- bare heteroatom fallbacks ------------------------------------ */
  if (atoms.length === 1 && mol.totalH(start) > 0) {
    if (el === 'N') return plain('amino');
    if (el === 'O') return plain('hydroxy');
    if (el === 'S') return plain('sulfanyl');
  }

  /* --- ordinary carbon / ring substituent --------------------------- */
  var unit = ONG.naming.nameUnit(ctx, atoms, start);
  if (!unit) {
    ctx.warn('이름을 붙일 수 없는 치환기가 있습니다.');
    return null;
  }
  return unit;

  function plain(name) { return {name: name, sortName: ONG.naming.sortKeyOf(name), complex: false}; }
  function complexName(name) { return {name: name, sortName: ONG.naming.sortKeyOf(name), complex: true}; }
}

function wrap(named) {
  return /[-\d\s]/.test(named.name) ? '(' + named.name + ')' : named.name;
}

function restOf(atoms, remove) {
  return atoms.filter(function (a) { return remove.indexOf(a) < 0; });
}

/** Substituent names on an amine / amide nitrogen. */
function nSubstituents(ctx, atoms, nitrogen, cameFrom) {
  var mol = ctx.mol, out = [];
  mol.heavyNeighbors(nitrogen).forEach(function (j) {
    if (j === cameFrom) return;
    if (atoms.indexOf(j) < 0) return;
    var sub = nameSub(ctx, ONG.naming.branchAtoms(mol, j, indexOf([nitrogen, cameFrom]), indexOf(atoms.concat([nitrogen]))), j, j);
    if (sub) out.push(wrap(sub));
  });
  return out.sort();
}

function indexOf(list) {
  var s = {};
  list.forEach(function (i) { s[i] = true; });
  return s;
}

function joinMultiplied(names) {
  var counts = {};
  names.forEach(function (n) { counts[n] = (counts[n] || 0) + 1; });
  return Object.keys(counts).sort().map(function (n) {
    return (counts[n] > 1 ? ONG.naming.multiplier(counts[n], /[-\d(]/.test(n)) : '') + n;
  }).join('');
}

/** The parent hydride behind a substituent name: methyl -> methane. */
function hydrideNameOf(unit) {
  var cand = unit.evaluation && unit.evaluation.cand;
  if (!cand) return null;
  if (cand.kind === 'ring') {
    if (cand.ringKind === 'benzene') return 'benzene';
    if (cand.ringKind === 'template') return cand.info.name;
    return 'cyclo' + ONG.naming.stemFor(cand.atoms.length) + 'ane';
  }
  if (/[-\d(]/.test(unit.name)) return null;      // substituted: keep it simple
  return ONG.naming.stemFor(cand.atoms.length) + 'ane';
}

/** "acetyl", "propanoyl", "benzoyl" for an acyl group. */
function acylName(ctx, atoms, group) {
  var mol = ctx.mol;
  var acylAtoms = atoms.filter(function (a) { return a !== group.esterO; });
  var unit = ONG.naming.nameUnit(ctx, acylAtoms, group.carbon);
  if (!unit) return null;
  // nameUnit gives us e.g. "methyl"/"ethyl"; rebuild the acyl form from the chain length.
  var chain = unit.evaluation && unit.evaluation.cand;
  if (!chain) return null;
  var base;
  if (chain.kind === 'ring') base = ONG.naming.elide(unit.name.replace(/yl$/, 'e'), 'c') + 'carbonyl';
  else base = ONG.naming.stemFor(chain.atoms.length) + 'anoyl';
  if (chain.kind === 'ring' && chain.ringKind === 'benzene') base = 'benzenecarbonyl';
  return RETAINED_ACYL[base] || base;
}

function nameSub(ctx, atoms, start) {
  if (!atoms || !atoms.length) return null;
  return ONG.naming.nameUnit(ctx, atoms, start);
}

ONG.nameBranch = nameBranch;
if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
