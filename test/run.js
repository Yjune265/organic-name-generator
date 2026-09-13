// Runs every test file. Usage: node test/run.js
var path = require('path');
var files = ['naming.test.js', 'stereo.test.js', 'variants.test.js'];
var failed = 0;
files.forEach(function (f) {
  delete require.cache[require.resolve('./' + f)];
  process.exitCode = 0;
  require('./' + f);
  if (process.exitCode) failed++;
});
process.exitCode = failed ? 1 : 0;
console.log(failed ? '\n' + failed + ' test file(s) failed' : '\nall test files passed');
