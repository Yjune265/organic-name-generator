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
  carboxylicAcid: 1,
  ester: 2,
  acylHalide: 3,
  amide: 4,
  nitrile: 5,
  aldehyde: 6,
  ketone: 7,
  alcohol: 8,
  thiol: 9,
  amine: 10
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
    var dblO = null, dblS = null, singleOH = null, singleOR = null, singleN = null, halo = null, tripleN = null;
    links.forEach(function (l) {
      var nb = mol.atoms[l.atom];
      if (l.bond.order === 2 && nb.element === 'O') dblO = l.atom;
      else if (l.bond.order === 2 && nb.element === 'S') dblS = l.atom;
      else if (l.bond.order === 3 && nb.element === 'N') tripleN = l.atom;
      else if (l.bond.order === 1 && nb.element === 'O') {
        if (mol.totalH(l.atom) > 0 && mol.heavyNeighbors(l.atom).length === 1) singleOH = l.atom;
        else if (!inRing[l.atom] && mol.heavyNeighbors(l.atom).length === 2) singleOR = l.atom;
      } else if (l.bond.order === 1 && nb.element === 'N' && !inRing[l.atom]) {
        if (isPlainNitrogen(mol, l.atom, i)) singleN = l.atom;
      } else if (l.bond.order === 1 && mol.isHalogen(l.atom)) halo = l.atom;
    });

    if (tripleN !== null && mol.heavyNeighbors(tripleN).length === 1) {
      claim({type: 'nitrile', carbon: i, nitrogen: tripleN, atoms: [i, tripleN]});
      continue;
    }
    if (dblO === null) continue;

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
