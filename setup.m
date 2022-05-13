function setup(varargin)
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
% environmental variable PATH so that nvcc can be called
%
% SETUP CUDA NVCC_PATH - specifies the path for the nvcc executable. On
% Linux the default is '/usr/local/cuda/bin'. On Windows it is the first
% installation under
% 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\**\bin'
%
% SETUP CUDA NVCC_PATH MSVC_CL_PATH - specifies the path for the C/C++
% compilier under MSVC. By default this is the first installation under 
% 'C:\Program Files (x86)\Microsoft Visual Studio *\**\cl.exe'
%
% See also TEARDOWN

base_path = fileparts(mfilename('fullpath'));
rel_paths = {'.', 'fun', 'src', 'utils'};
paths = cellfun(@(p)fullfile(base_path, p), rel_paths, 'UniformOutput', false);
paths = paths(7 == cellfun(@exist, paths));
addpath(paths{:});

i = 1;
while i <= nargin % go through arguments sequentially
    switch varargin{i}
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
            us = UltrasoundSystem(); % needs paths to have been added already
            us.recompile(); % attempt to recompile code
            copyfile(us.tmp_folder, fullfile(base_path, "bin")); % copy it
            addpath(fullfile(base_path, 'bin')); % add new cache path
        
        case 'CUDA' % add CUDA to the path
            
            % get the nvcc executable path
            if isunix
                if i+1 <= nargin && isfolder(varargin{i+1}) % if next arg is user provided path
                    p = varargin{i+1}; % store user path
                else 
                    p = "/usr/local/cuda/bin"; % linux default nvcc path
                end                
                if ~exist(fullfile(p, 'nvcc')), warning("nvcc not found at " + p); end
                p = pathsep + p; % prepend path separator here
            elseif ispc
                % get the nvcc path
                if  i+1 <= nargin && isfolder(varargin{i+1}) % if next arg is user provided path
                    p1 = varargin{i+1}; % store user path
                    i = i + 1; % move to next argument
                else % find the windows (default) nvcc paths
                    p1 = getenv('CUDA_PATH'); % use env. var if set
                    if isfolder(p1)
                        p1 = fullfile(p1, 'bin'); % bin should have nvcc
                    else
                        l = dir(fullfile("C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA","**","nvcc*")); % search for nvcc
                        if ~isempty(l) % we found at least one
                            p1 = l(1).folder;  % grab the 1st folder - TODO: grab the best instead
                        else % nothing found
                            p1 = repmat("", [0, 1]); % empty string
                        end
                    end
                end
                if ~isempty(p1) || ~exist(fullfile(p1, 'nvcc.exe')), warning("nvcc not found at " + p1); end
                
                % MSVC find cl.exe
                if i+1 <= nargin && isfolder(varargin{i+1}) % if next arg is user provided path
                    p2 = varargin{i+1}; % store user path
                    i = i + 1; % move to next argument
                else % search for cl.exe within MSVC default install directories
                    l = dir(fullfile('C:\Program Files (x86)\Microsoft Visual Studio *','**', 'cl.exe'));                
                    if ~isempty(l) % we found at least one
                        p2 = l(1).folder;  % grab the 1st folder - TODO: grab the best instead
                    else % nothing found
                        p2 = repmat("", [0, 1]); % empty string
                    end
                end
                if ~isempty(p2) || ~exist(fullfile(p2, 'cl.exe')), warning("cl not found at " + p2); end
                
                % join nvcc and CUDA paths (if they exist)
                p = strjoin([(pathsep + p1); (pathsep + p2)],'');                
                
            else 
                error('CUDA compilation paths undefined if not unix or pc.');
            end
            
            % add them to the system path
            % TODO: keep track of this so that we can remove it during
            % teardown.m
            setenv('PATH', fullfile(getenv('PATH'), p));
            
        otherwise
            error("Unrecognized option " + varargin{i} + ".");
    end
    i = i + 1; % move to next argument
end

