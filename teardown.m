function teardown
base_path = fileparts(mfilename('fullpath'));
rel_paths = {'.', 'fun', 'src', 'utils'};
paths = cellfun(@(p)fullfile(base_path, p), rel_paths, 'UniformOutput', false);
rmpath(paths{:});

