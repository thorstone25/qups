function setup(opts)
% SETUP - Setup the workspace for QUPS
%
% SETUP - by itself adds relevant paths for QUPS classes and functions.
%
% SETUP parallel - additionally initiates a parallel pool that is used by 
% parfor loops.
%
% SETUP cache - recompiles CUDA and mex files and creates a bin folder to 
% be used to cache the files. 
%
% SETUP CUDA - adds the default CUDA installation paths to the system
% environmental variable PATH so that nvcc can be called. The function will
% attempt to find an installation of CUDA with nvcc and add it to the
% system PATH. On Windows, it will also attempt to find an installation of
% MSVC C/C++. If the 'MW_NVCC_PATH' environment variable is set, it will
% use that version of CUDA first. On Windows, if the 'VCToolsInstallDir'
% environmental variable is set, it will use that C/C++ compiler.
%
% Note: In MATLAB R2023a+ ptx compilation in provided via mexcuda, which
% will also use the 'MW_NVCC_PATH' environment variable if set.
%
% SETUP disable-gpu - disables gpu support by shadowing `gpuDeviceCount` so
% that it always returns 0. This prevents implicit gpu usage by some
% functions.
% 
% SETUP enable-gpu undoes the effects of the above.
% 
% SETUP disable-ocl - disables gpu support by shadowing `oclDeviceCount` so
% that it always returns 0, similar to the disable-gpu option. This 
% prevents implicit oclDevice and oclKernel usage and checks, which may be
% necessary to run some function with a parallel.ThreadPool.
% 
% SETUP enable-ocl undoes the effects of the above.
% 
% See also TEARDOWN

arguments(Repeating)
    opts (1,1) string {mustBeMember(opts, ["CUDA", "cache", "parallel", "disable-gpu", "disable-ocl", "enable-gpu", "enable-ocl", "no-path"])}
end

base_path = string(fileparts(mfilename('fullpath')));
prj = matlab.project.rootProject;
nms = recursiveProjectName(prj); % get all open project names
if ~any(cellfun(@(o) o=="no-path", opts)) % don't modify paths if asked not to
    if isempty(prj)
        openProject(fullfile(base_path, 'Qups.prj')); % open the project
    elseif ismember(nms, "qups") % qups already open
        % if no argumnets, this is likely a re-initialization
        if isempty(opts), warning("QUPS:setup:AlreadyInitialized","QUPS is already initialized here: '" + prj.RootFolder + "'."), end
    else % addpaths manually, so as not to disturb open projects
        rel_paths = [".", "kern", "src", "utils", "bin"];
        paths = fullfile(base_path, rel_paths);
        paths = paths(7 == arrayfun(@exist, paths));
        addpath(paths{:});
    end
end

i = 1;
while i <= nargin % go through arguments sequentially
    switch opts{i}
        case 'parallel'
            hcp = gcp('nocreate');
            if isvalid(hcp)
                warning('A valid parallel pool already exists.');
            else
                if verLessThan('matlab', '9.8.0') % R2020a
                    parpool('local', 'SpmdEnabled', false); % backwards compatibles
                else
                    parpool('threads'); % thread-based pool - no copying temp variables
                end
            end
            
        case 'cache'
            us = UltrasoundSystem('recompile', false); % needs paths to have been added already
            us.recompile(); % attempt to recompile code
            copyfile(us.tmp_folder, fullfile(base_path, "bin")); % copy it
            % addpath(fullfile(base_path, 'bin')); % add new cache path
        
        case 'CUDA' % add CUDA to the path
            
            % get the nvcc executable path
            if isunix
                p = getenv('MW_NVCC_PATH'); % use env. var if set
                if ~isfolder(p) % bin should have nvcc
                    p = "/usr/local/cuda/bin"; % linux default nvcc path
                end                
                if ~exist(fullfile(p, 'nvcc')), warning("nvcc not found at " + p); end
            
            elseif ispc
                % get all the windows drives, from A to Z
                isdrive = @(c) logical(exist([c ':\'], "dir"));
                wdrvs = char(double('A'):double('Z'));
                wdrvs = wdrvs(arrayfun(isdrive, wdrvs));

                % get the nvcc path
                % find the windows (default) nvcc paths
                p1 = getenv('MW_NVCC_PATH'); % use env. var if set
                if ~isfolder(p1) % should have nvcc
                    l = arrayfun(@(d) {dir(fullfile(d + ":\Program Files\NVIDIA GPU Computing Toolkit\CUDA","**","nvcc*"))}, wdrvs); % search for nvcc
                    l = cat(1, l{:});
                    if ~isempty(l) % we found at least one
                        p1 = l(1).folder;  % grab the 1st folder - TODO: grab the best instead
                    else % nothing found
                        p1 = repmat("", [1, 0]); % empty string
                    end
                end

                if isempty(p1),                          warning("nvcc not found.")
                elseif ~exist(fullfile(p1, 'nvcc.exe')), warning("nvcc not found at " + p1); 
                end
                p1 = string(p1); % enforce string type for casting/sizing
                
                % MSVC find cl.exe
                % search for cl.exe within MSVC default install directories
                p2 = getenv('VCToolsInstallDir');
                if isfolder(p2)
                    p2 = fullfile(p2, 'bin', 'Hostx86','x64');
                else
                    l = arrayfun(@(d) {dir(fullfile(d + ":\Program Files*\Microsoft Visual Studio*",'**','Hostx86','x64', 'cl.exe'))}, wdrvs);
                    l = cat(1, l{:});
                    if ~isempty(l) % we found at least one
                        p2 = l(1).folder;  % grab the 1st folder - TODO: grab the best instead
                    else % nothing found
                        p2 = repmat("", [1, 0]); % empty string
                    end
                end
                
                if isempty(p2),                        warning("cl not found.");
                elseif ~exist(fullfile(p2, 'cl.exe')), warning("cl not found at " + p2); end
                p2 = string(p2); % enforce string type for casting/sizing
                
                % join nvcc and CUDA paths (if they exist)
                p = strjoin([p1, p2],pathsep);
                
            else 
                warning('CUDA compilation paths undefined if not unix or pc.');
                i = i + 1;
                continue;
            end
            
            % add them to the system path
            % TODO: keep track of this so that we can remove it during
            % teardown.m
            setenv('PATH', join([p, string(getenv('PATH'))],pathsep));
            
        case {"disable-gpu", "disable-ocl", "enable-gpu", "enable-ocl"}
            
            dev = extractAfter(opts{i}, "able-"); % dev type
            fl = fullfile(base_path,"utils",dev+"DeviceCount.m"); % file
            switch extractBefore(opts{i}, "-")
                case "disable"
                    % shadow gpu/ocl support by setting the device count to 0
                    writelines("function n = "+dev+"DeviceCount(), n = 0;", fl);
                case "enable"
                    % delete the shadowing file (if it exists)
                    delete(fl);                    
            end
        case "no-path" % pass
        otherwise
            error("Unrecognized option " + opts{i} + ".");
    end
    i = i + 1; % move to next argument
end

function nms = recursiveProjectName(prj)
arguments 
    prj matlab.project.Project
end
if isempty(prj)
    nms = string.empty;
else
    nms = cellfun(@recursiveProjectName, {prj.ProjectReferences.Project}, 'UniformOutput', false);
    nms = unique([prj.Name, nms{:}]);
end


