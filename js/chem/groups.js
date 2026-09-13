/**
 * groups.js - functional group perception.
 *
 * Ring heteroatoms belong to the parent hydride, so they are never reported
 * as ethers/amines; a lactone or lactam therefore falls out as a ring
 * ketone (oxolan-2-one, piperidin-2-one), which is what IUPAC asks for.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

// Seniority of the classes that can be expressed as a suffix (low = senior).
var SUFFIX_RANK = {
  ammonium: 1,
  carboxylicAcid: 2,
  sulfonicAcid: 3,
  sulfinicAcid: 4,
  anhydride: 5,
  ester: 6,
  sulfonateEster: 7,
  acylHalide: 8,
  sulfonylHalide: 9,
  amide: 10,
  sulfonamide: 11,
  nitrile: 12,
  aldehyde: 13,
  ketone: 14,
  thione: 15,
  alcohol: 16,
  thiol: 17,
  amine: 18,
  imine: 19
};

var HALIDE_WORD = {F: 'fluoride', Cl: 'chloride', Br: 'bromide', I: 'iodide'};
var HALO_PREFIX = {F: 'fluoro', Cl: 'chloro', Br: 'bromo', I: 'iodo'};

function perceiveGroups(mol) {
  mol.perceiveAromaticity();
  var groups = [], owner = {}, i;

  function claim(group) {
    group.rank = SUFFIX_RANK[group.type] || 99;
    group.idx = groups.length;
    groups.push(group);
    (group.atoms || []).forEach(function (a) { if (owner[a] === undefined) owner[a] = group.idx; });
    return group;
  }

  var inRing = mol.atoms.map(function (a, idx) { return mol.isInRing(idx); });

  /* --- carbonyl and nitrile carbons ---------------------------------- */
  for (i = 0; i < mol.atoms.length; i++) {
    var a = mol.atoms[i];
    if (a.element !== 'C') continue;
    var links = mol.links(i);
    var dblO = null, dblS = null, dblN = null, singleOH = null, singleOR = null,
        singleN = null, halo = null, tripleN = null;
    links.forEach(function (l) {
      var nb = mol.atoms[l.atom];
      if (l.bond.order === 2 && nb.element === 'O') dblO = l.atom;
      else if (l.bond.order === 2 && nb.element === 'S') dblS = l.atom;
      else if (l.bond.order === 2 && nb.element === 'N' && !inRing[l.atom]) dblN = l.atom;
      else if (l.bond.order === 3 && nb.element === 'N') tripleN = l.atom;
      else if (l.bond.order === 1 && nb.element === 'O') {
        if (mol.totalH(l.atom) > 0 && mol.heavyNeighbors(l.atom).length === 1) singleOH = l.atom;
        else if (!inRing[l.atom] && mol.heavyNeighbors(l.atom).length === 2) singleOR = l.atom;
      } else if (l.bond.order === 1 && nb.element === 'N' && !inRing[l.atom]) {
        if (isPlainNitrogen(mol, l.atom, i)) singleN = l.atom;
      } else if (l.bond.order === 1 && mol.isHalogen(l.atom)) halo = l.atom;
    });

    if (tripleN !== null && mol.heavyNeighbors(tripleN).length === 1) {
      // R-S-C#N is a thiocyanate and R-O-C#N a cyanate: a real nitrile
      // carbon is attached to carbon (or to nothing else at all).
      var attached = mol.heavyNeighbors(i).filter(function (j) { return j !== tripleN; });
      var onCarbon = attached.every(function (j) { return mol.atoms[j].element === 'C'; });
      if (onCarbon) {
        claim({type: 'nitrile', carbon: i, nitrogen: tripleN, atoms: [i, tripleN]});
        continue;
      }
    }
    // Two cumulated heteroatom double bonds (O=C=N, S=C=N) belong to an
    // isocyanate / isothiocyanate, not to a carbonyl or a thione.
    var heteroDoubles = [dblO, dblS, dblN].filter(function (x) { return x !== null; });
    if (heteroDoubles.length >= 2) continue;

    if (dblO === null) {
      // C=S is a thione, C=N an imine (both rank below the carbonyl classes).
      if (dblS !== null && mol.heavyNeighbors(dblS).length === 1) {
        claim({type: 'thione', carbon: i, sulfur: dblS, atoms: [i, dblS]});
      } else if (dblN !== null && isImineNitrogen(mol, dblN, i)) {
        claim({type: 'imine', carbon: i, nitrogen: dblN, atoms: [i, dblN]});
      }
      continue;
    }

    // R-CO-O-CO-R is an anhydride, not two esters.
    if (singleOR !== null && bridgesTwoAcyls(mol, singleOR, i)) {
      var partner = mol.heavyNeighbors(singleOR).filter(function (j) { return j !== i; })[0];
      if (i < partner) {
        claim({type: 'anhydride', carbon: i, partner: partner, bridge: singleOR,
               atoms: [i, dblO, singleOR]});
      }
      continue;
    }

    if (singleOH !== null) {
      claim({type: 'carboxylicAcid', carbon: i, atoms: [i, dblO, singleOH]});
    } else if (singleOR !== null) {
      claim({type: 'ester', carbon: i, esterO: singleOR, atoms: [i, dblO, singleOR]});
    } else if (halo !== null) {
      claim({type: 'acylHalide', carbon: i, halogen: halo,
             halide: HALIDE_WORD[mol.atoms[halo].element], atoms: [i, dblO, halo]});
    } else if (singleN !== null) {
      claim({type: 'amide', carbon: i, nitrogen: singleN, atoms: [i, dblO, singleN]});
    } else if (mol.totalH(i) >= 1) {
      claim({type: 'aldehyde', carbon: i, atoms: [i, dblO]});
    } else {
      claim({type: 'ketone', carbon: i, atoms: [i, dblO]});
    }
  }

  /* --- sulfonic / sulfinic acids and their derivatives ---------------- */
  for (i = 0; i < mol.atoms.length; i++) {
    if (mol.atoms[i].element !== 'S' || owner[i] !== undefined || inRing[i]) continue;
    var sLinks = mol.links(i);
    var oxo = [], hydroxy = null, alkoxy = null, carbon = null, nitrogen = null, sHalo = null;
    sLinks.forEach(function (l) {
      var nb = mol.atoms[l.atom];
      if (nb.element === 'O' && l.bond.order === 2) oxo.push(l.atom);
      else if (nb.element === 'O' && mol.totalH(l.atom) > 0 && mol.heavyNeighbors(l.atom).length === 1) {
        hydroxy = l.atom;
      } else if (nb.element === 'O' && mol.heavyNeighbors(l.atom).length === 2) alkoxy = l.atom;
      else if (nb.element === 'C') carbon = l.atom;
      else if (nb.element === 'N' && isSulfonamideNitrogen(mol, l.atom)) nitrogen = l.atom;
      else if (mol.isHalogen(l.atom)) sHalo = l.atom;
    });
    if (carbon === null || !oxo.length) continue;
    var base = [i].concat(oxo);
    if (oxo.length === 2 && hydroxy !== null) {
      claim({type: 'sulfonicAcid', sulfur: i, carrier: carbon, atoms: base.concat([hydroxy])});
    } else if (oxo.length === 1 && hydroxy !== null) {
      claim({type: 'sulfinicAcid', sulfur: i, carrier: carbon, atoms: base.concat([hydroxy])});
    } else if (oxo.length === 2 && alkoxy !== null) {
      claim({type: 'sulfonateEster', sulfur: i, carrier: carbon, esterO: alkoxy,
             atoms: base.concat([alkoxy])});
    } else if (oxo.length === 2 && sHalo !== null) {
      claim({type: 'sulfonylHalide', sulfur: i, carrier: carbon, halogen: sHalo,
             halide: HALIDE_WORD[mol.atoms[sHalo].element], atoms: base.concat([sHalo])});
    } else if (oxo.length === 2 && nitrogen !== null) {
      claim({type: 'sulfonamide', sulfur: i, carrier: carbon, nitrogen: nitrogen,
             atoms: base.concat([nitrogen])});
    }
  }

  /* --- quaternary ammonium ------------------------------------------- */
  for (i = 0; i < mol.atoms.length; i++) {
    var atom = mol.atoms[i];
    if (atom.element !== 'N' || atom.charge !== 1 || owner[i] !== undefined) continue;
    if (mol.heavyNeighbors(i).length !== 4) continue;
    if (mol.bondOrderSum(i) !== 4) continue;
    claim({type: 'ammonium', nitrogen: i, atoms: [i]});
  }

  /* --- azide, isocyanate, nitroso and friends ------------------------ */
  for (i = 0; i < mol.atoms.length; i++) {
    if (owner[i] !== undefined || inRing[i]) continue;
    var el0 = mol.atoms[i].element;
    if (el0 === 'N') {
      var azide = azideChain(mol, i);
      if (azide) { claim({type: 'azide', nitrogen: i, prefix: 'azido', atoms: azide}); continue; }
      var cumulated = cumulatedCarbon(mol, i);
      if (cumulated) {
        claim({type: cumulated.type, nitrogen: i, prefix: cumulated.prefix, atoms: cumulated.atoms});
        continue;
      }
      var nitrosoO = mol.links(i).filter(function (l) {
        return l.bond.order === 2 && mol.atoms[l.atom].element === 'O' &&
          mol.heavyNeighbors(l.atom).length === 1;
      });
      if (nitrosoO.length === 1 && mol.heavyNeighbors(i).length === 2) {
        claim({type: 'nitroso', nitrogen: i, prefix: 'nitroso', atoms: [i, nitrosoO[0].atom]});
        continue;
      }
    }
    if (el0 === 'S') {
      var thiocyanate = thiocyanateChain(mol, i);
      if (thiocyanate) {
        claim({type: 'thiocyanate', sulfur: i, prefix: 'thiocyanato', atoms: thiocyanate});
      }
    }
  }

  /* --- nitro --------------------------------------------------------- */
  for (i = 0; i < mol.atoms.length; i++) {
    if (mol.atoms[i].element !== 'N' || owner[i] !== undefined || inRing[i]) continue;
    var oxygens = mol.neighbors(i).filter(function (j) {
      return mol.atoms[j].element === 'O' && mol.heavyNeighbors(j).length === 1;
    });
    if (oxygens.length === 2 && mol.heavyNeighbors(i).length === 3) {
      claim({type: 'nitro', nitrogen: i, atoms: [i].concat(oxygens)});
    }
  }

  /* --- hydroxy, thiol, ether, sulfide, amine ------------------------- */
  for (i = 0; i < mol.atoms.length; i++) {
    var el = mol.atoms[i].element;
    if (owner[i] !== undefined) continue;
    var heavy = mol.heavyNeighbors(i);

    if (el === 'O' && !inRing[i]) {
      if (mol.totalH(i) > 0 && heavy.length === 1) {
        claim({type: 'alcohol', oxygen: i, carrier: heavy[0],
               phenol: mol.atoms[heavy[0]].aromatic, atoms: [i]});
      } else if (heavy.length === 2) {
        claim({type: 'ether', oxygen: i, atoms: [i]});
      }
    } else if (el === 'S' && !inRing[i]) {
      if (mol.totalH(i) > 0 && heavy.length === 1) {
        claim({type: 'thiol', sulfur: i, carrier: heavy[0], atoms: [i]});
      } else if (heavy.length === 2) {
        claim({type: 'sulfide', sulfur: i, atoms: [i]});
      }
    } else if (el === 'N' && !inRing[i]) {
      if (mol.atoms[i].charge) continue;
      if (!isPlainNitrogen(mol, i, -1)) continue;
      claim({type: 'amine', nitrogen: i, atoms: [i]});
    }
  }

  /* --- halogens as prefixes ------------------------------------------ */
  for (i = 0; i < mol.atoms.length; i++) {
    if (owner[i] === undefined && mol.isHalogen(i)) {
      claim({type: 'halide', atom: i, prefix: HALO_PREFIX[mol.atoms[i].element], atoms: [i]});
    }
  }

  var principal = null;
  groups.forEach(function (g) {
    if (!SUFFIX_RANK[g.type]) return;
    if (!principal || g.rank < principal) principal = g.rank;
  });

  return {
    all: groups,
    owner: owner,
    principalRank: principal,
    ofType: function (type) { return groups.filter(function (g) { return g.type === type; }); },
    principal: function () {
      return principal === null ? [] : groups.filter(function (g) { return g.rank === principal; });
    },
    groupOwning: function (atom) {
      return owner[atom] === undefined ? null : groups[owner[atom]];
    }
  };
}

/** The oxygen of R-CO-O-CO-R bridges two acyl carbons. */
function bridgesTwoAcyls(mol, oxygen, carbon) {
  var others = mol.heavyNeighbors(oxygen).filter(function (j) { return j !== carbon; });
  return others.length === 1 && isCarbonyl(mol, others[0]);
}

/** C=N-R / C=N-H, but not an amide, oxime ether or ring member. */
function isImineNitrogen(mol, n, carbon) {
  if (mol.atoms[n].charge) return false;
  var links = mol.links(n);
  if (links.length > 2) return false;
  for (var k = 0; k < links.length; k++) {
    var l = links[k];
    if (l.atom === carbon) continue;
    var el = mol.atoms[l.atom].element;
    if (l.bond.order !== 1) return false;
    if (el !== 'C' && el !== 'H') return false;      // =N-OH is an oxime, not an imine
  }
  return true;
}

/** -N=[N+]=[N-] or -N-N#N */
function azideChain(mol, n1) {
  var heavy = mol.heavyNeighbors(n1);
  if (heavy.length !== 2) return null;
  var carbon = heavy.filter(function (j) { return mol.atoms[j].element === 'C'; });
  var n2 = heavy.filter(function (j) { return mol.atoms[j].element === 'N'; });
  if (carbon.length !== 1 || n2.length !== 1) return null;
  var rest = mol.heavyNeighbors(n2[0]).filter(function (j) { return j !== n1; });
  if (rest.length !== 1 || mol.atoms[rest[0]].element !== 'N') return null;
  if (mol.heavyNeighbors(rest[0]).length !== 1) return null;
  return [n1, n2[0], rest[0]];
}

/** -N=C=O and -N=C=S */
function cumulatedCarbon(mol, n) {
  var heavy = mol.heavyNeighbors(n);
  if (heavy.length !== 2) return null;
  var carbons = heavy.filter(function (j) {
    return mol.atoms[j].element === 'C' && mol.bondOrder(n, j) === 2;
  });
  if (carbons.length !== 1) return null;
  var c = carbons[0];
  var tail = mol.heavyNeighbors(c).filter(function (j) { return j !== n; });
  if (tail.length !== 1 || mol.bondOrder(c, tail[0]) !== 2) return null;
  var el = mol.atoms[tail[0]].element;
  if (el === 'O') return {type: 'isocyanate', prefix: 'isocyanato', atoms: [n, c, tail[0]]};
  if (el === 'S') return {type: 'isothiocyanate', prefix: 'isothiocyanato', atoms: [n, c, tail[0]]};
  return null;
}

/** -S-C#N */
function thiocyanateChain(mol, s) {
  var heavy = mol.heavyNeighbors(s);
  if (heavy.length !== 2) return null;
  for (var k = 0; k < heavy.length; k++) {
    var c = heavy[k];
    if (mol.atoms[c].element !== 'C' || mol.bondOrder(s, c) !== 1) continue;
    var tail = mol.heavyNeighbors(c).filter(function (j) { return j !== s; });
    if (tail.length !== 1 || mol.bondOrder(c, tail[0]) !== 3) continue;
    if (mol.atoms[tail[0]].element !== 'N') continue;
    return [s, c, tail[0]];
  }
  return null;
}

/** N of a sulfonamide: single bonds only, to the sulfur plus C/H. */
function isSulfonamideNitrogen(mol, n) {
  if (mol.atoms[n].charge) return false;
  var links = mol.links(n), sulfurs = 0;
  for (var k = 0; k < links.length; k++) {
    if (links[k].bond.order !== 1) return false;
    var el = mol.atoms[links[k].atom].element;
    if (el === 'S') sulfurs++;
    else if (el !== 'C' && el !== 'H') return false;
  }
  return sulfurs === 1;
}

/** True when N carries only single bonds to C/H (an amine-type nitrogen). */
function isPlainNitrogen(mol, n, exceptCarbon) {
  var links = mol.links(n), ok = true;
  links.forEach(function (l) {
    if (l.bond.order !== 1) ok = false;
    var el = mol.atoms[l.atom].element;
    if (el !== 'C' && el !== 'H' && el !== 'N') ok = false;
    // A nitrogen already bonded to another carbonyl carbon is an amide N.
    if (el === 'C' && l.atom !== exceptCarbon && isCarbonyl(mol, l.atom)) ok = false;
  });
  return ok;
}

function isCarbonyl(mol, c) {
  if (mol.atoms[c].element !== 'C') return false;
  return mol.links(c).some(function (l) {
    return l.bond.order === 2 && mol.atoms[l.atom].element === 'O';
  });
}

ONG.perceiveGroups = perceiveGroups;
ONG.SUFFIX_RANK = SUFFIX_RANK;
ONG.HALO_PREFIX = HALO_PREFIX;
if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
