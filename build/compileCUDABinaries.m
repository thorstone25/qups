function compileCUDABinaries(ofl)
arguments
    ofl (1,1) string = fullfile(fileparts(mfilename('fullpath')),"..","bin"+filesep)
end
% compileCUDABinaries - Compile CUDA binaries for the oldest supported CC.

% architectures (last is oldest support)
arch = "compute_" + [90 89 87 86 80 75 72 70 60];% 52 50];

% compile
us = UltrasoundSystem('recompile', false);
defs = UltrasoundSystem.genCUDAdefs(); % definition structs
com = us.recompileCUDA(defs, arch(end)); % compile
fls = fullfile(us.tmp_folder, "*.ptx");

% copy to bin folder
copyfile(fls, ofl);

end