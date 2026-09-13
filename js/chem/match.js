/**
 * match.js - graph isomorphism used for compound lookup and for mapping
 * named ring templates onto a drawn ring system.
 */
(function (global) {
'use strict';

var ONG = global.ONG = global.ONG || {};

/**
 * All isomorphisms from molA onto molB.
 *
 * opts.atomLabel(mol, idx) -> string   comparison label (default: Mol#atomLabel)
 * opts.limit                           stop after N mappings (default 1)
 *
 * Returns an array of mappings; mapping[aIdx] = bIdx.
 */
function isomorphisms(molA, molB, opts) {
  opts = opts || {};
  var limit = opts.limit || 1;
  var label = opts.atomLabel || function (mol, i) { return mol.atomLabel(i); };
  var n = molA.atoms.length;
  if (n !== molB.atoms.length || molA.bonds.length !== molB.bonds.length) return [];

  var labelsA = [], labelsB = [], i;
  for (i = 0; i < n; i++) { labelsA.push(label(molA, i)); labelsB.push(label(molB, i)); }
  if (labelsA.slice().sort().join(';') !== labelsB.slice().sort().join(';')) return [];

  var colorsA = molA.colors(), colorsB = molB.colors();
  var histA = {}, histB = {};
  colorsA.forEach(function (c) { histA[c] = (histA[c] || 0) + 1; });
  colorsB.forEach(function (c) { histB[c] = (histB[c] || 0) + 1; });
  var sameColorHistogram = Object.keys(histA).length === Object.keys(histB).length &&
    Object.keys(histA).every(function (k) { return histA[k] === histB[k]; });
  // Colours are only a pruning aid; callers using custom labels opt out.
  var useColors = opts.ignoreColors ? false : sameColorHistogram;
  if (!opts.ignoreColors && !sameColorHistogram) return [];

  // Most constrained atoms first shrinks the search dramatically.
  var order = [];
  for (i = 0; i < n; i++) order.push(i);
  order.sort(function (x, y) {
    var d = molA.degree(y) - molA.degree(x);
    if (d) return d;
    return (histA[colorsA[y]] || 0) - (histA[colorsA[x]] || 0);
  });
  // Keep the search connected: each new atom should touch an already mapped one.
  var seq = [], used = {};
  while (seq.length < n) {
    var pick = -1;
    for (i = 0; i < order.length; i++) {
      var cand = order[i];
      if (used[cand]) continue;
      if (!seq.length) { pick = cand; break; }
      if (molA.neighbors(cand).some(function (nb) { return used[nb]; })) { pick = cand; break; }
    }
    if (pick < 0) { for (i = 0; i < order.length; i++) if (!used[order[i]]) { pick = order[i]; break; } }
    used[pick] = true;
    seq.push(pick);
  }

  var mapping = new Array(n), taken = new Array(n), results = [];
  for (i = 0; i < n; i++) { mapping[i] = -1; taken[i] = false; }

  function search(k) {
    if (results.length >= limit) return;
    if (k >= n) { results.push(mapping.slice()); return; }
    var a = seq[k];
    for (var b = 0; b < n; b++) {
      if (taken[b]) continue;
      if (labelsA[a] !== labelsB[b]) continue;
      if (useColors && colorsA[a] !== colorsB[b]) continue;
      if (molA.degree(a) !== molB.degree(b)) continue;
      var ok = true, nbs = molA.links(a);
      for (var m = 0; m < nbs.length && ok; m++) {
        var mapped = mapping[nbs[m].atom];
        if (mapped < 0) continue;
        var other = molB.bondBetween(b, mapped);
        if (!other) ok = false;
        else if (!opts.ignoreBondOrder && other.order !== nbs[m].bond.order) ok = false;
      }
      // Nothing already mapped may bond to b unless it bonds to a as well.
      if (ok) {
        var nbsB = molB.links(b);
        for (var q = 0; q < nbsB.length && ok; q++) {
          var back = mapping.indexOf(nbsB[q].atom);
          if (back >= 0 && !molA.bondBetween(a, back)) ok = false;
        }
      }
      if (!ok) continue;
      mapping[a] = b; taken[b] = true;
      search(k + 1);
      mapping[a] = -1; taken[b] = false;
      if (results.length >= limit) return;
    }
  }
  search(0);
  return results;
}

function isIsomorphic(molA, molB, opts) {
  return isomorphisms(molA, molB, opts).length > 0;
}

ONG.isomorphisms = isomorphisms;
ONG.isIsomorphic = isIsomorphic;
if (typeof module !== 'undefined' && module.exports) module.exports = ONG;

})(typeof globalThis !== 'undefined' ? globalThis : this);
