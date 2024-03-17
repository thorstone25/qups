function err = teardown()
% TEARDOWN - Remove QUPS from the path
%
% teardown() closes the QUPS projects or removes QUPS folders from
% the path directly if the QUPS project is not the current project.
%
% err = teardown() additionally returns true if the paths cannot be removed
% due to project references or false otherwise.
%
% Note: this function may be removed in a future release.
% 
% Example:
% setup CUDA cache; % one-time setup
% ... QUPS code here ...
% teardown; % remove all references to QUPS
% 
% See also SETUP
prj = matlab.project.rootProject;
nms = recursiveProjectName(prj);
err = false;

if isempty(prj) % no current project
    rmproj(); % remove paths directly
elseif prj.Name == "qups" % is the current project
    close(prj);
elseif ismember(nms, "qups") % must be referenced by current project
    warning( ...
        "QUPS:teardown:ReferencedProject", ...
        "Cannot teardown QUPS because it is refernced by the project " + prj.Name + "." ...
        );
    err = true;
else % not referenced by currnet project
    rmproj();
end

function rmproj()
base_path = fileparts(mfilename("fullpath"));
rel_paths = [".", "bin", "kern", "src", "utils"];
paths = fullfile(base_path, rel_paths);
paths = paths(7 == arrayfun(@exist, paths)); % existing paths only
rmpath(paths{:});


function nms = recursiveProjectName(prj)
arguments 
    prj matlab.project.Project
end
if isempty(prj)
    nms = string.empty;
else
    nms = cellfun(@recursiveProjectName, {prj.ProjectReferences.Project});
    nms = unique([prj.Name, nms]);
end


