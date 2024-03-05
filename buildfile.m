function plan = buildfile
import matlab.buildtool.tasks.CleanTask
import matlab.buildtool.tasks.CodeIssuesTask
import matlab.buildtool.tasks.TestTask

% Create a plan from task functions
plan = buildplan(localfunctions);
base = plan.RootFolder; % base directory

%% Installing Extensions
% extension names
ext_nms = ["FieldII" "kWave" "OpenCL", "MUST", "USTB"];

plan(  "install_FieldII").Outputs = fullfile(base, "../FieldII");
plan(  "install_kWave"  ).Outputs = fullfile(base, "../k-wave");
plan(  "install_OpenCL" ).Outputs = fullfile(base, "../Matlab-OpenCL");
plan(  "install_MUST"   ).Outputs = fullfile(base, "../MUST");
plan(  "install_USTB"   ).Outputs = fullfile(base, "../ustb");

plan("uninstall_FieldII").Inputs  = plan("install_FieldII").Outputs;
plan("uninstall_kWave"  ).Inputs  = plan("install_kWave"  ).Outputs;
plan("uninstall_OpenCL" ).Inputs  = plan("install_OpenCL" ).Outputs;
plan("uninstall_MUST"   ).Inputs  = plan("install_MUST"   ).Outputs;
plan("uninstall_USTB"   ).Inputs  = plan("install_USTB"   ).Outputs;

for nm = ext_nms, plan(  "install_"+nm).Description = "Download and install " + nm; end
for nm = ext_nms, plan("uninstall_"+nm).Description = "Uninstall " + nm; end

plan("install") = matlab.buildtool.Task( ...
    "Description", "Download and install all extensions", ...
    "Dependencies", "install_"+ext_nms ...
    );

plan("uninstall") = matlab.buildtool.Task( ...
    "Description", "Uninstall all extensions", ...
    "Dependencies", "uninstall_"+ext_nms ...
    );

%% Testing and Artifacts
% source files (folders)
src = ["kern","src","utils"] + filesep;

% Add a task to identify code issues
plan("check") = CodeIssuesTask([src, "examples"+filesep]);

% Add a task to remove output files and projects
plan("clean") = CleanTask("Dependencies","uninstall");

% Add a task to run tests
plan("coverage") = TestTask( ...
    "Tag","full","SourceFiles",src ...
    ,CodeCoverageResults=fullfile("build","code-coverage",["coverage.xml" "html/index.html"]) ...
    ,TestResults="build/test-results/report.html" ...
    );

% Add a task to run brief tests
plan("test") = TestTask( ...
    "Tag","Github","SourceFiles",src ...
    ,CodeCoverageResults=fullfile("build","code-coverage-github",["coverage.xml" "html/index.html"]) ...
    ,TestResults="build/test-results-github/report.html" ...
    );

%% Add compilation dependencies

% Make the "compile" task dependent on the "check and "test" tasks
% plan("compile").Dependencies = ["test", "check"];0

% make dummy compile mex task
mfls = dir(fullfile(base, "src", "**", "msfm*.c"));
[mfls, nms] = deal(string(fullfile({mfls.folder}, {mfls.name})), string({mfls.name}));
ofls = replace(fullfile(base, "bin", nms), '.c', "."+mexext());
tnm = "compile_mex_"+extractBefore(nms,'.c');
for i = 1:numel(mfls)
    plan(tnm(i)) = matlab.buildtool.tasks.MexTask(mfls(i), fileparts(ofls(i)), ...
        "Options",join("-I"+fullfile(base,"src","FMM","functions")) ...
        );
    plan(tnm(i)).Outputs = ofls(i);
end

plan("compile_mex") = matlab.buildtool.Task( ...
    "Description","Compile mex kernels.", ...
    "Dependencies", tnm ...
    );

% get CUDA files
fls = dir(fullfile(base, "src", "**", "*.cu")); % inputs
[fls, nms] = deal(string(fullfile({fls.folder}, {fls.name})), string({fls.name}));
ofls = replace(fullfile(base, "bin", nms), '.cu', '.ptx'); % outputs
ofls(endsWith(ofls, "sizes.ptx")) = []; % delete - no matching ptx output

% mark expected inputs/outputs
plan("compile_CUDA").Inputs  = fullfile(fls);
plan("compile_CUDA").Outputs = fullfile(ofls);

% make dependent tasks
plan("compile") = matlab.buildtool.Task( ...
    "Description", "Compile all kernels", ...
    "Dependencies", "compile_"+["CUDA", "mex"] ...
    );

%% Final
% Mark the default task in the plan
plan.DefaultTasks = "compile";

% Full release, without fresh install
plan("release") = matlab.buildtool.Task( ...
    "Description", "Check code, compile kernels, and report coverage", ...
    "Dependencies", ["check", "compile", "coverage"] ...
    );

end

function archiveTask(~)
% Create ZIP file
filename = "source_" + ...
    string(datetime("now",Format="yyyyMMdd'T'HHmmss"));
zip(filename,"*")
end

function compile_CUDATask(context, arch)
arguments
    context
    arch = "compute_" + [90 89 87 86 80 75 72 70 60];% 52 50]; % CC values we support
end
% Compile CUDA kernels
setup CUDA no-path; % add CUDA and US to path

% compile
defs = UltrasoundSystem.genCUDAdefs(); % definition structs
us = UltrasoundSystem('recompile', false);
us.recompileCUDA(defs, arch(end)); % compile

% copy files to bin
fls = fullfile(us.tmp_folder, "*.ptx");
ofl = fullfile(context.Plan.RootFolder,"bin"); % most back-compatible
copyfile(fls, ofl);

end

function   install_FieldIITask(context),   installer(context, "FieldII"); end % Download and install FieldII
function   install_OpenCLTask( context),   installer(context, "OpenCL" ); end % Download and install Matlab-OpenCL
function   install_kWaveTask(  context),   installer(context, "kWave"  ); end % Download and install k-Wave
function   install_MUSTTask(   context),   installer(context, "MUST"   ); end % Download and install MUST
function   install_USTBTask(   context),   installer(context, "USTB"   ); end % Download and install USTB

function uninstall_FieldIITask(context), uninstaller(context, "FieldII"); end % Uninstall FieldII
function uninstall_OpenCLTask( context), uninstaller(context, "OpenCL" ); end % Uninstall Matlab-OpenCL
function uninstall_kWaveTask(  context), uninstaller(context, "kWave"  ); end % Uninstall k-Wave
function uninstall_MUSTTask(   context), uninstaller(context, "MUST"   ); end % Uninstall MUST
function uninstall_USTBTask(   context), uninstaller(context, "USTB"   ); end % Uninstall USTB


function uninstaller(context, ext)
arguments
    context
    ext (1,1) string
end

prj = matlab.project.rootProject; % get current project
epth = context.Task.Inputs(1).Path; % extension path
try % try to remove - hard to tell if it's already referenced via path alone
    prj.removeReference(epth);
catch ME
    if ME.identifier == "MATLAB:project:api:IsNotProjectReference"
        warning("Cannot find " + ext + ": it may already be uninstalled.");
    else, rethrow(ME)
    end
    disp(ME);
end % try to remove via path
end

function installer(context, ext, proj)
arguments
    context 
    ext (1,1) string {mustBeMember(ext, ["FieldII", "kWave", "MUST", "USTB", "OpenCL"])}
    proj (1,1) logical = true; % with a project
end

og = context.Plan.RootFolder; % base directory
fld = context.Task.Outputs(1).Path; % output folder for the project (relative)

% try to fetch the git repo if it's already there
try repo = gitrepo(fld); dwnld = false; catch ME, dwnld = true; end

% if unable to initialize repo, create one
if dwnld
    switch ext
        case "FieldII"
            % download into local folder
            % fld = "../FieldII";
            if ispc,       os = "windows";
            elseif ismac,  os = "mac";
            elseif isunix, os = "linux";
            else, error("Cannot identify system type for FieldII download.");
            end

            url = "https://www.field-ii.dk/program_code/matlab_2021/Field_II_ver_3_30_"+os+".tar.gz";
            untar(url, fld);

            % create a git repo
            vnm = replace(extractBetween(url, "Field_II_ver_", "_"+os),"_","."); % version name
            repo = gitinit(fld, "InitialBranch", "main");
            repo.add(fld); % add all files
            repo.commit(Message="Initial commit - " + vnm);

        case "MUST"
            url = "https://www.biomecardio.com/MUST/functions/MUST.zip";
            unzip(url, fld);

            % create a git repo
            repo = gitinit(fld, "InitialBranch", "main");
            repo.add(fld); % add all files
            repo.commit(Message="Initial commit - " + string(datetime()))

        case "kWave"
            url = "https://github.com/ucl-bug/k-wave.git";
            repo = gitclone(url, fld);

        case "OpenCL"
            url = "https://github.com/thorstone25/Matlab-OpenCL.git";
            repo = gitclone(url, fld, RecurseSubmodules=true);

        case "USTB"
            url = "https://bitbucket.org/ustb/ustb.git";
            repo = gitclone(url, fld);

            % move the USTB.prj file to USTB-pck.prj
            % HACK: git-mv not available yet: git-rm, mv, git-add (should be same effect)
            repo.rm( fullfile(fld,"USTB.prj")); % delete old location
            movefile(fullfile(fld,"USTB.prj"), fullfile(fld,"USTB-pkg.prj")); % move
            repo.add(fullfile(fld,"USTB-pkg.prj")); % add updated location

        otherwise
            error("Project " + ext + " not implemented ");
    end
end
% pth - path to folder containing project

% make a project
if ~proj, return; end % short circuit if not making it a project dependency

if isempty(dir(fullfile(fld, "*"+ext+"*.prj"))) ... % proper project file and folders
        || ~exist(fullfile(fld,'resources'), 'dir') % resource folders

    % project folders to add to path, w.r.t. base path (`pth`)
    switch ext
        case "kWave", pflds = fullfile("k-Wave",[".", "binaries"]);
        otherwise, pflds = ["."]; % only the base directory
    end

    % make a project
    eprj = matlab.project.createProject("Folder",pth,"Name",ext);
    arrayfun(@(p) eprj.addPath(fullfile(pth, p)), pflds); % add project folders

else % project should already exist
    eprj = matlab.project.loadProject(pth);
end

% add as a QUPS reference
prj = openProject(fullfile(og, "Qups.prj"));
if ~any(eprj == [prj.ProjectReferences.Project]) % check if already installed
    addReference(prj, eprj);
else
    warning(ext + " already installed.");
end
end
