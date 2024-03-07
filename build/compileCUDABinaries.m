function [s, com] = compileCUDABinaries(ofl)
% compileCUDABinaries - Compile CUDA binaries for the oldest supported CC.
%
% s = compileCUDABinaries(outdir) compiles the binaries and places them in
% directory outdir and returns the success status s for each of the
% compialtion defintions. The default directory is the 'bin' folder.
%
% [s, com] = compileCUDABinaries(...) additionally returns the compilation
% command stirngs com.
%
% See also UltrasoundSystem.genCUDAdefs UltrasoundSystem.getSrcFolder
arguments
    ofl (1,1) string = fullfile(fileparts(mfilename('fullpath')),"..","bin"+filesep)
end

% architectures (last is oldest support)
arch = "compute_" + [90 89 87 86 80 75 72 70 60];% 52 50];

%% compile
% us = UltrasoundSystem('recompile', false);
defs = UltrasoundSystem.genCUDAdefs(); % definition structs
src = UltrasoundSystem.getSrcFolder(); % source folder
com = arrayfun(@(d)makeDefCommand(d, arch(end), src, ofl), defs); % compile each

% try to compile
try s = arrayfun(@system, com); return;
catch ME % not sure why ...
    rethrow(ME)
end

function com = makeDefCommand(d, arch, src_folder, bin_folder)
arguments
    d (1,1) struct
    arch (1,1) string
    src_folder {mustBeFolder} = UltrasoundSystem.getSrcFolder();
    bin_folder {mustBeFolder} = src_folder
end
% trun def into a full command
com = join(cat(1,...
    "nvcc ", ...
    "-arch=" + arch + " ", ... compile for active gpu
    "-o " + fullfile(bin_folder, strrep(d.Source, '.cu', '.ptx')), ...
    join("--" + d.CompileOptions),...
    join("-I" + d.IncludePath), ...
    join("-L" + d.Libraries), ...
    join("-W" + d.Warnings), ...
    join("-D" + d.DefinedMacros),...
    "--ptx " + fullfile(src_folder, d.Source) ...
    ));
