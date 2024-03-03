function compileCUDABinaries()
% compileCUDABinaries - Compile CUDA binaries for the oldest supported CC.

% architectures (last is oldest support)
arch = "compute_" + [90 89 87 86 80 75 72 70 60];% 52 50];

% set paths
prj = matlab.project.rootProject;
if isempty(prj) || prj.Name ~= "qups"
    prj = openProject(which('Qups.prj'));
end
setup CUDA no-path;

% compile
us = UltrasoundSystem();
defs = UltrasoundSystem.genCUDAdefs(); % definition structs
com = us.recompileCUDA(defs, arch(end)); % compile

% copy to bin folder
copyfile(fullfile(us.tmp_folder, "*.ptx"), fullfile(fileparts(mfilename('fullpath')),"..","bin"+filesep));

end