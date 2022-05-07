function setup(varargin)
base_path = fileparts(mfilename('fullpath'));
rel_paths = {'.', 'fun', 'src', 'utils'};
paths = cellfun(@(p)fullfile(base_path, p), rel_paths, 'UniformOutput', false);
paths = paths(7 == cellfun(@exist, paths));
addpath(paths{:});

for i = 1:nargin
    switch varargin{i}
        case 'parallel'
            hcp = gcp('nocreate');
            if isvalid(hcp)
                warning('A valid parallel pool already exists.');
            else 
                if verLessThan('matlab', '9.8.0') % R2020a
                    parpool('local', 'SpmdEnabled', false);
                else
                    parpool('threads');
                end
            end
        case 'cache'
            us = UltrasoundSystem(); % needs paths to have been added already
            us.recompile();
            copyfile(us.tmp_folder, fullfile(base_path, "bin"));
            addpath(fullfile(base_path, 'bin')); % added new cache path
        otherwise 
            error('Unrecognized option'); 
    end
end

