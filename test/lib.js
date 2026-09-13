// Shared loader for the Node test scripts.
['mol', 'smiles', 'match', 'stereo', 'rings', 'groups', 'naming', 'branches', 'compounds', 'variants', 'api']
  .forEach(function (f) { require('../js/chem/' + f + '.js'); });
module.exports = globalThis.ONG;
