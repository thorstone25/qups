
function res = runProjectTests(tag, parenv, kwargs)
% runProjectTests - Run tests for QUPS
%
% runProjectTests Runs the full suite of tests for QUPS. These
% include default initialization, simulation, beamforming, and some
% computational kernel tests.
%
% runProjectTests('Github') runs a subset of the tests that complete within
% a more rasonable amount of time.
%
% runProjectTests('full') runs the full suite.
% 
% runProjectTests(tag, parenv) uses the given parallel environment.
% 
% runProjectTests(..., 0, 'report', true) creates a coverage report or 
% runProjectTests(..., 0, 'cobertura', true) creates a cobertura format
% coverage.xml file. These options and running in parallel together are
% mutually exclusive.
% 
% res = runProjectTests(...) returns the test results array res.

arguments
    tag (1,1) string {mustBeMember(tag, ["full", "Github"])} = "full"
    parenv {mustBeScalarOrEmpty, mustBeA(parenv, ["double","parallel.Pool"])} = gcp('nocreate')
    kwargs.cobertura (1,1) logical = false
    kwargs.report (1,1) logical = false
end

% create a test runner
runner = matlab.unittest.TestRunner.withTextOutput();

% add coverage reporting
flds = ["kern", "src", "utils"];
if kwargs.cobertura
    fmt = matlab.unittest.plugins.codecoverage.CoberturaFormat('coverage.xml');
    runner.addPlugin(matlab.unittest.plugins.CodeCoveragePlugin.forFolder(flds, "Producing", fmt));
end
if kwargs.report
    fmt = matlab.unittest.plugins.codecoverage.CoverageReport('CoverageReport');
    runner.addPlugin(matlab.unittest.plugins.CodeCoveragePlugin.forFolder(flds, "Producing", fmt));
end

% run all the tests
if isempty(parenv) || (parenv == 0)
    suite = matlab.unittest.TestSuite.fromProject(fileparts(which("Qups.prj")), "Tag", tag);
    res = runner.run(suite); % no parallel
elseif isnumeric(parenv) % make tmp pool
    hcp = parpool(parenv);
    suite = matlab.unittest.TestSuite.fromProject(fileparts(which("Qups.prj")), "Tag", tag);
    res = runner.runInParallel(suite); % proc
    delete(hcp);

elseif isa(parenv, "parallel.ProcessPool")
    suite = matlab.unittest.TestSuite.fromProject(fileparts(which("Qups.prj")), "Tag", tag);
    res = runner.runInParallel(suite); % proc

elseif isa(parenv, "parallel.ThreadPool")   
    prj = openProject(fileparts(which("Qups.prj")));
    fls = [prj.Files.Path]; % all paths
    clnms = argn(2, @fileparts, fls(endsWith(fls, '.m'))); % filenames
    clnms = clnms(logical(arrayfun(@(n) exist(n,'class'), clnms))); % class files
    clnms = clnms(arrayfun(@(c) isa(meta.class.fromName(c), "matlab.unittest.meta.class"), clnms)); % test class files
    clss = arrayfun(@meta.class.fromName, clnms);
    suite = arrayfun(@(c)matlab.unittest.TestSuite.fromClass(c, "Tag", tag), clss, 'UniformOutput',false);
    suite = [suite{:}];

    res = runner.runInParallel(suite);
else
    error('Ambiguous pool option - set parenv directly');
end

