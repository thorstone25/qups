function teardown
base_path = fileparts(mfilename('fullpath'));
rel_paths = {'.', 'bin', 'fun', 'src', 'utils'};
paths = cellfun(@(p)fullfile(base_path, p), rel_paths, 'UniformOutput', false);
paths = paths(7 == cellfun(@exist, paths)); % existing paths only
rmpath(paths{:});

