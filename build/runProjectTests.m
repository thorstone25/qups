
function [res, suite] = runProjectTests(tag, parenv, kwargs)
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
% runProjectTests(tag, parenv) uses the given parallel environment parenv.
% 
% runProjectTests(..., 0, 'cobertura', true) creates a cobertura format
% coverage.xml file. 
% 
% runProjectTests(..., 0, 'report', true) creates a coverage report. This
% option and running in parallel are mutually exclusive.
% 
% runProjectTests(..., 0, 'verbosity', level) sets the verbosity level.
% 
% res = runProjectTests(...) returns the test results array res.
%
% [res, suite] = runProjectTests(...) also returns tbe TestSuite suite.
% 
% See also matlab.unittest.Verbosity matlab.unittest.TestRunner matlab.unittest.plugins.codecoverage

arguments
    tag (1,1) string {mustBeMember(tag, ["full", "Github", "benchmark"])} = "full"
    parenv {mustBeScalarOrEmpty, mustBeA(parenv, ["double","parallel.Pool"])} = gcp('nocreate')
    kwargs.cobertura (1,1) logical = false
    kwargs.report (1,1) logical = false
    kwargs.verbosity (1,1) matlab.unittest.Verbosity = matlab.unittest.Verbosity.Concise
end

% create a test runner
runner = matlab.unittest.TestRunner.withNoPlugins();
plugs = { ...
    matlab.unittest.plugins.TestRunProgressPlugin.withVerbosity(kwargs.verbosity), ...
    matlab.unittest.plugins.DiagnosticsRecordingPlugin("LoggingLevel", kwargs.verbosity, "OutputDetail", kwargs.verbosity), ...
    matlab.unittest.plugins.LoggingPlugin.withVerbosity(kwargs.verbosity), ...
    };

% add coverage reporting
flds = ["kern", "src", "utils"];
if kwargs.cobertura
    fmt = matlab.unittest.plugins.codecoverage.CoberturaFormat('build/coverage.xml');
    plugs{end+1} = (matlab.unittest.plugins.CodeCoveragePlugin.forFolder(flds, "Producing", fmt));
end
if kwargs.report
    fmt = matlab.unittest.plugins.codecoverage.CoverageReport('build/CoverageReport');
    plugs{end+1} = (matlab.unittest.plugins.CodeCoveragePlugin.forFolder(flds, "Producing", fmt));
end

% filter plugs if running in parallel
if ~isempty(parenv) || (isnumeric(parenv) && parenv > 0) % delete non-parallel plug-ins
    i = ~cellfun(@(p)p.supportsParallel, plugs); % is paralle
    if any(i), warning("Deleting plug-ins that lack parallel support."); end
    plugs(i) = [];
end 
cellfun(@runner.addPlugin, plugs); % add each plugin

% run all the tests
if isempty(parenv) || (isnumeric(parenv) && (parenv == 0))
    suite = matlab.unittest.TestSuite.fromProject(fileparts(which("Qups.prj")), "Tag", tag);
    res = runner.run(suite); % no parallel
elseif isnumeric(parenv) % make tmp pool
    hcp = parpool(parenv);
    suite = matlab.unittest.TestSuite.fromProject(fileparts(which("Qups.prj")), "Tag", tag);
    res = runner.runInParallel(suite); % proc
    delete(hcp);

elseif isa(parenv, "parallel.ProcessPool")
    prj = openProject(fileparts(which("Qups.prj"))); % load project
    suite = matlab.unittest.TestSuite.fromProject(prj.RootFolder, "Tag", tag);
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
    error("Unrecognized pool of type '" + class(parenv) + "'.");
end

