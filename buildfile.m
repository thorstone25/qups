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

plan("patch_MUST"  ).Inputs = plan("install_MUST"   ).Outputs;
plan("patch_kWave" ).Inputs = plan("install_kWave"  ).Outputs;
plan("patch_OpenCL").Inputs = plan("install_OpenCL" ).Outputs;

plan("benchmark").Outputs   = fullfile(base, "build", "benchmark.log");


for nm = ext_nms, plan(  "install_"+nm).Description = "Download and install " + nm; end
for nm = ext_nms, plan("uninstall_"+nm).Description = "Uninstall " + nm; end

% aggregate
plan("install") = matlab.buildtool.Task( ...
    "Description", "Download and install all extensions", ...
    "Dependencies", "install_"+ext_nms ...
    );

plan("uninstall") = matlab.buildtool.Task( ...
    "Description", "Uninstall all extensions", ...
    "Dependencies", "uninstall_"+ext_nms ...
    );

plan("patch") = matlab.buildtool.Task( ...
    "Description", "Patch extensions", ...
    "Dependencies", "patch_"+[ext_nms([2 4 3])] ...
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
    "Tag","syntax","SourceFiles",src ...
    ,CodeCoverageResults=fullfile("build","code-coverage-brief",["coverage.xml" "html/index.html"]) ...
    ,TestResults="build/test-results-brief/report.html" ...
    ,Dependencies="compatability" ...
    ,OutputDetail="Concise" ...
    );

%% Add compilation dependencies

% Make the "compile" task dependent on the "check and "test" tasks
% plan("compile").Dependencies = ["test", "check"];

mfls = dir(fullfile(base, "src", "**", "msfm*.c"));
[mfls, nms] = deal(string(fullfile({mfls.folder}, {mfls.name})), string({mfls.name}));
ofls = fullfile(base, "bin", replace(nms, ".c" + lineBoundary("end"), "."+mexext()));
tnm = "compile_mex_"+extractBefore(nms,'.c');
for i = 1:numel(mfls)
    plan(tnm(i)) = matlab.buildtool.tasks.MexTask(mfls(i), fileparts(ofls(i)), ...
        "Options",join("-I"+fullfile(base,"src","FMM","functions")) ...
        );
    plan(tnm(i)).Outputs = ofls(i);
end

% make dummy compile mex task
plan("compile_mex") = matlab.buildtool.Task( ...
    "Description","Compile MEX kernels.", ...
    "Dependencies", tnm ...
    );

% get CUDA files
fls = dir(fullfile(base, "src", "**", "*.cu")); % inputs
fls(endsWith({fls.name}, "sizes.cu")) = []; % delete - no matching ptx output
[fls, nms] = deal(string(fullfile({fls.folder}, {fls.name})), string({fls.name}));
ofls = replace(fullfile(base, "bin", nms), '.cu', '.ptx'); % outputs

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
    "Description", "Check code, compile kernels, syntax test and report coverage", ...
    "Dependencies", ["check", "compile", "test"] ...
    );

% all the things (not the default!)
plan("all") = matlab.buildtool.Task( ...
    "Description", "Install and patch all extensions, check code, compile kernels, full test and report coverage", ...
    "Dependencies", ["install", "patch", "release"] ...
    );

end

function archiveTask(~)
% Create ZIP file
filename = "source_" + ...
    string(datetime("now",Format="yyyyMMdd'T'HHmmss"));
zip(filename,"*")
end

function compile_CUDATask(context, arch)
% compile CUDA kernels
arguments
    context matlab.buildtool.TaskContext
    arch (1,1) string = "compute_"+60; % + ;
end
supp = "compute_"+[90 89 87 86 80 75 72 70 60];  % supported CC values
if ~ismember(arch, supp)
    warning("QUPS:build:compileCUDA:unsupportedArchitecture", ...
        "Expected the architecture to be one of: "+newline+join("'"+supp+"'",newline));
end

odirs = fileparts(string({context.Task.Outputs.Path}));
ifls  = string({context.Task.Inputs.Path});
defs = UltrasoundSystem.genCUDAdefs(); % matching definition structs

try % via mexcuda
    for i = numel(defs):- 1:1
        d = defs(i);
        ifl = ifls(endsWith(ifls, d.Source));
        com{i} = cat(1,...
            ifl, ...
            "-ptx",...
            "-outdir", odirs(i), ...
            ...join("--" + d.CompileOptions),...
            join("-I" + d.IncludePath), ...
            join("-L" + d.Libraries), ...
            ...join("-W" + d.Warnings), ...
            join("-D" + d.DefinedMacros)...
            );
    end

    args = cellfun(@cellstr, com, 'UniformOutput', false);
    cellfun(@(args) mexcuda(args{:}), args);

catch ME % via nvcc
    warning("QUPS:build:"+ME.identifier, ...
        join(["Unable to compile via mexcuda:",...
        ME.message, ...
        "Attempting command line compilation instead."...
        ], newline));

    % Compile CUDA kernels
    setup CUDA no-path; % add CUDA and US to path

    % compile
    us = UltrasoundSystem('recompile', false);
    us.recompileCUDA(defs, arch(end)); % compile

    % copy files to bin
    fls = fullfile(us.tmp_folder, "*.ptx");
    ofl = fullfile(context.Plan.RootFolder,"bin"); % most back-compatible
    copyfile(fls, ofl);
end

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

            isM1 = ismac() && computer() == "MACA64"; % M1 chips
            if isM1
                url = "https://www.field-ii.dk/program_code/matlab_2023/Field_II_ver_4_11_"+os+".zip";
                unzip(url, fld);
                movefile(fullfile(fld, "m_files","*"), fld); % move sub-folder files to main folder
            else
                url = "https://www.field-ii.dk/program_code/matlab_2021/Field_II_ver_3_30_"+os+".tar.gz";
                untar(string(gunzip(url, fld)), fld);
            end

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
            repo.commit(Message="Initial commit - " + string(datetime()));

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

try eprj = matlab.project.loadProject(fld); % if project already exists
catch ME % else probably no project, so make one
    if ~ismember(ME.identifier, "MATLAB:project:api:LoadFail"), rethrow(ME); end % something else?
    % project folders to add to path, w.r.t. base path (`pth`)
    switch ext
        case "kWave", pflds = fullfile("k-Wave",[".", "binaries"]);
        otherwise, pflds = "."; % only the base directory
    end

    % make a project
    eprj = matlab.project.createProject("Folder",fld,"Name",ext);
    arrayfun(@(p) eprj.addPath(fullfile(fld, p)), pflds); % add project folders

end

% add as a QUPS reference
prj = openProject(fullfile(og, "Qups.prj"));
prj.addReference(eprj); % save to add if already exists
end

function patch_MUSTTask(context)
% Patch MUST to remove GUI calls and enable ThreadPools
fld = context.Task.Inputs.Path; % output folder for the project (relative)
fls = (dir(fullfile(fld, '**', 'pfield*.m')));
fls = string(fullfile({fls.folder}, {fls.name}));

% shadow problematic functions
pnc_dsp = whitespacePattern(0,Inf) + ("," | "%" | lineBoundary("end")); % implicit display call punctuation
cll = optionalPattern("("+wildcardPattern+")") + whitespacePattern(0,Inf) + (";" | "," | "%" | lineBoundary("end"));
shdw = ["AdMessage", "MUSTstat"]; % to shadow
for fl = fls
    munlock(fl); clear(fl); % refresh
    txt = string(readlines(fl));
    for s = shdw
        j = find(contains(txt,   s + cll) & ~contains(txt, "function")); % find the line with the call
        fcn    = extract(txt(j), s + cll); % extract function call
        fcn    = replace(fcn, pnc_dsp, ";"); % silence the output
        txt(j) = replace(txt(j), s + cll, s+"=@(varargin)deal([]); "+fcn); % replace call
    end
    writelines(txt,fl); % save modified file
end

end
function patch_OpenCLTask(context)
% Patch MatCL to silence mexPrintf debug statements
fld = context.Task.Inputs.Path; % output folder for the project (relative)
fls = (dir(fullfile(fld,"sub","MatCL", '**', 'cl_*.*pp')));
fls = string(fullfile({fls.folder}, {fls.name}));

% shadow problematic functions
pat = lineBoundary('start') + whitespacePattern + "mexPrintf";
str = "// mexPrintf";
for fl = fls
    txt = string(readlines(fl));
    txt = replace(txt, pat, str);
    writelines(txt,fl); % save modified file
end

end

function patch_kWaveTask(context)
% Patch k-Wave to avoid race conditions and enable parallel simulations
fld = context.Task.Inputs.Path; % output folder for the project (relative)
fl = (dir(fullfile(fld, '**', 'kspaceFirstOrder3DC.m')));
fl = string(fullfile({fl.folder}, {fl.name}));

% replace tempdir with new tempname folder, and truncate before deletion
txt = string(readlines(fl));
pat = ["data_path = tempdir;", "delete(input_filename);", "delete(output_filename);"];
rep = ["data_path = tempname;" + newline + " mkdir(data_path);", ...
    "if isunix, system(""truncate -s 0 "" +  input_filename); end; " + pat(2) ...
    "if isunix, system(""truncate -s 0 "" + output_filename); end; " + pat(3) + newline + " rmdir(data_path);" ...    
    ];
for i = 1:length(pat)
    txt = replace(txt, pat(i), rep(i)); % find & replace
end
writelines(txt, fl); % save modified file

fls = (dir(fullfile(fld, '**', 'kspaceFirstOrder_*nput*.m'))); % private input checking files
fls = string(fullfile({fls.folder}, {fls.name}));

% replace calls to 'ver' or 'verLessThan' - QUPS requires parallel
% computing anyway, and this prevents ThreadPools
pat = ["v = ver;" "ismember('Parallel Computing Toolbox',"+wildcardPattern()+")", "verLessThan("+wildcardPattern()+")"];
rep = ["v = false;","true(1)", "false(1)"];
for fl = fls
    txt = string(readlines(fl));
    for i = 1:length(pat)
        txt = replace(txt, pat(i), rep(i)); % find & replace
    end
    writelines(txt, fl); % save modified file
end

% if recording pressure vs. time on gpu, keep storage on CPU when
% allocating memory
fl = dir(fullfile(fld, '**', "kspaceFirstOrder_createStorageVariables.m"));
fl = string(fullfile({fl.folder}, {fl.name}));
pat = "        sensor_data.p = castZeros([num_sensor_points, num_recorded_time_points]);";
rep = join([
    "        sensor_data.p = castZeros([num_sensor_points, 1]);", ...
    "        if isa(sensor_data.p, 'gpuArray'), sensor_data.p = gather(sensor_data.p); end % always retain on CPU", ...
    "        sensor_data.p = repmat(sensor_data.p, [1, num_recorded_time_points]);", ...
    ], newline);
txt = string(readlines(fl));
txt = replace(txt, pat, rep); % find & replace
writelines(txt, fl); % save modified file

end

function benchmarkTask(context)
% Benchmark simulation and beamforming routines
arguments
    context matlab.buildtool.TaskContext
end

ofl  = context.Task.Outputs.Path;
res = runProjectTests("benchmark", 0);
txt = string({cat(2,cat(2,res.Details).DiagnosticRecord).Report})';
writelines(txt, ofl);

end

function compatabilityTask(context)
% Check for required toolboxes
req = ["Signal Processing", "Parallel Computing"];
rec = ["Image Processing", "Statistics and Machine Learning"];
sug = reshape(string.empty,1,0);
pck = {req,rec,sug}; % packages
sev = ["required", "recommended", "suggested"]; % severity
rep = "The following packages are " + sev ... 
    + ": " + sprintf('\t') + cellfun(@(s)join(s,", ",2),pck);

disp(newline); arrayfun(@disp, rep(~ismissing(rep))); disp(newline);
vn = string({ver().Name});

for i = 1:numel(sev)
    if isempty(pck{i}), continue; end
    j = ismember(pck{i} + " Toolbox", vn);
    if all(j)
        disp("All "+sev(i)+" toolboxes are installed.");
    else
        msg = ("The following toolboxes are "+sev(i)+", but not installed: "...
            +newline+join(pck{i}(~j)+ " Toolbox", newline)+newline);
        switch sev(i)
            case "required",    error(  msg);
            case "recommended", warning(msg);
            otherwise,          disp(   msg);
        end
    end
end
disp(newline);

end
