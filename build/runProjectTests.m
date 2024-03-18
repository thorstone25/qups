
function [res, suite, runner] = runProjectTests(tag, parenv, kwargs)
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
% runProjectTests(..., 'report', reports) adds plugins to create
% coverage or test reports in the corresponding format. These options and
% running in parallel may be mutually exclusive.
% 
% runProjectTests(..., 'verbosity', level) sets the verbosity level.
% 
% runProjectTests(..., 'coverage', cov) adds a CoverageResult plugin cov to
% the tests. This can later be used to generate further reports.
% 
% runProjectTests(..., 'coverage', cov, 'MetricLevel', lvl) records results
% at the MetricLevel lvl. The default is "statement".
% 
% res = runProjectTests(...) returns the test results array res.
%
% [res, suite] = runProjectTests(...) also returns the TestSuite suite.
% 
% [res, suite, runner] = runProjectTests(...) also returns the TestRunner
% runner.
% 
% [res, suite, runner] = runProjectTests(..., 'dryrun', true) sets up the
% TestSuite suite and TestRunner runner but does not run the tests. The
% runner can be configured with further plugins and run with: 
%     `res = runner.run(suite);` 
% or
%     `res = runner.runInParallel(suite);`
% when ready.
% 
% Example:
% % Generate reports for the "syntax" tag
% parpool local; % setup acceleration (optional)
% res = runProjectTests("syntax", "report",  ["test-html", "cov-report"]);
% 
%
% % Example 2 (requires R2023a+ & MATLAB Test Toolbox): 
% % Generate more detailed results 
% cov = matlab.unittest.plugins.codecoverage.CoverageResult;
% fld = "build/Coverage-Report-Example"; % output report folder
% runProjectTests("full", coverage=cov, MetricLevel="condition");
% filename = generateHTMLReport(cov.Result, fld);
% open(filename); % open the report
% 
% See also matlab.unittest.Verbosity matlab.unittest.TestRunner
% matlab.unittest.plugins.codecoverage matlab.coverage.Result

arguments
    tag (1,1) string {mustBeMember(tag, ["syntax", "full", "Github", "build", "benchmark"])} = "syntax"
    parenv {mustBeScalarOrEmpty, mustBeA(parenv, ["double","parallel.Pool"])} = gcp('nocreate')
    kwargs.report (1,:) string {mustBeMember(kwargs.report,["test-html", "test-pdf", "cov-xml", "cov-report"])} = string.empty
    kwargs.verbosity (1,1) matlab.unittest.Verbosity = matlab.unittest.Verbosity.Concise
    kwargs.dryrun (1,1) logical = false
    kwargs.coverage (1,1) matlab.unittest.plugins.codecoverage.CoverageResult
    kwargs.MetricLevel (1,1) string {mustBeMember(kwargs.MetricLevel, ["statement","decision","condition","mcdc"])} = "statement"
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
if isfield(kwargs, 'coverage')
    fmt = kwargs.coverage; % returned by handle
    plugs{end+1} = (matlab.unittest.plugins.CodeCoveragePlugin.forFolder(flds, "Producing", fmt, "MetricLevel",kwargs.MetricLevel));
end
if ismember("cov-xml", kwargs.report)
    fmt = matlab.unittest.plugins.codecoverage.CoberturaFormat(fullfile('build','coverage-'+tag+'.xml'));
    plugs{end+1} = (matlab.unittest.plugins.CodeCoveragePlugin.forFolder(flds, "Producing", fmt));
end
if ismember("cov-report", kwargs.report)
    fmt = matlab.unittest.plugins.codecoverage.CoverageReport(fullfile('build','CoverageReport'+"-"+tag));
    plugs{end+1} = (matlab.unittest.plugins.CodeCoveragePlugin.forFolder(flds, "Producing", fmt));
end
if ismember("test-html", kwargs.report)
    plugs{end+1} = matlab.unittest.plugins.TestReportPlugin.producingHTML( ...
        fullfile('build',"test-results-"+tag), 'LoggingLevel',kwargs.verbosity);
end
if ismember("test-pdf", kwargs.report)
    plugs{end+1} = matlab.unittest.plugins.TestReportPlugin.producingPDF( ...
        fullfile('build',"test-results-"+tag+".pdf"), 'LoggingLevel',kwargs.verbosity);
end


% filter plugs if running in parallel
if ~isempty(parenv) || (isnumeric(parenv) && parenv > 0) % delete non-parallel plug-ins
    i = ~cellfun(@(p)p.supportsParallel, plugs); % is paralle
    if any(i), warning("Deleting plug-ins that lack parallel support."); end
    plugs(i) = [];
end 
cellfun(@runner.addPlugin, plugs); % add each plugin


% create the test suite
isparallel = ~(isempty(parenv) || (isnumeric(parenv) && (parenv == 0)));
if ~isparallel
    suite = matlab.unittest.TestSuite.fromProject(fileparts(which("Qups.prj")), "Tag", tag);
elseif isnumeric(parenv) % make tmp pool
    hcp = parpool(parenv);
    suite = matlab.unittest.TestSuite.fromProject(fileparts(which("Qups.prj")), "Tag", tag);

elseif isa(parenv, "parallel.ProcessPool")
    prj = openProject(fileparts(which("Qups.prj"))); % load project
    suite = matlab.unittest.TestSuite.fromProject(prj.RootFolder, "Tag", tag);

elseif isa(parenv, "parallel.ThreadPool")
    prj = openProject(fileparts(which("Qups.prj")));
    fls = [prj.Files.Path]; % all paths
    clnms = argn(2, @fileparts, fls(endsWith(fls, '.m'))); % filenames
    clnms = clnms(logical(arrayfun(@(n) exist(n,'class'), clnms))); % class files
    clnms = clnms(arrayfun(@(c) isa(meta.class.fromName(c), "matlab.unittest.meta.class"), clnms)); % test class files
    clss = arrayfun(@meta.class.fromName, clnms);
    suite = arrayfun(@(c)matlab.unittest.TestSuite.fromClass(c, "Tag", tag), clss, 'UniformOutput',false);
    suite = [suite{:}];

else
    error("Unrecognized pool of type '" + class(parenv) + "'.");
end

% run ... or not
if kwargs.dryrun
    res = matlab.unittest.TestResult.empty;
elseif ~isparallel
    res = runner.run(suite); % no parallel
else
    res = runner.runInParallel(suite); % in parallel   
    if exist('hcp','var'), delete(hcp); end % cleanup
end
