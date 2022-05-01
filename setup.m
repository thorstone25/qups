function setup(varargin)
base_path = fileparts(mfilename('fullpath'));
rel_paths = {'.', 'bin', 'fun', 'src', 'utils'};
paths = cellfun(@(p)fullfile(base_path, p), rel_paths, 'UniformOutput', false);

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
            us = UltrasoundSystem();
            us.recompile();
            copyfile(us.tmp_folder, fullfile(base_path, "bin"));
        otherwise 
            error('Unrecognized option'); 
    end
end

% add (existent) paths after (optionally) adding bin folder
paths = paths(7 == cellfun(@exist, paths));
addpath(paths{:});
