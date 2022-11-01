
function res = run_QUPS_tests(tag)
% RUN_QUPS_TESTS - Run tests for QUPS
%
% RUN_QUPS_TESTS Runs the full suite of tests for QUPS. These
% include default initialization, simulation, beamforming, and some
% computational kernel test.
%
% RUN_QUPS_TESTS('Github') runs a subset of the tests that complete within
% a more rasonable amount of time.
%
% res = RUN_QUPS_TESTS(...) returns the test results array res.

arguments, tag (1,1) string {mustBeMember(tag, ["full", "Github"])} = "full", end

% get all of the test
base_dir = fileparts(mfilename('fullpath'));
suite = matlab.unittest.TestSuite.fromFolder(fullfile(base_dir, 'test'));
% tags = string(unique([suite.Tags])) % show all the tags

% filter by the tag
suite = suite.selectIf(matlab.unittest.selectors.HasTag(tag));

% create a test runner
runner = matlab.unittest.TestRunner.withTextOutput();

% run all the tests
res = runner.run(suite);

