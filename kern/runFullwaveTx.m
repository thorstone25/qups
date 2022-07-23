
function runFullwaveTx(icmat, simdir, outdir)
% RUNFULLWAVETX Run a Fullwave simulation for a single transmit
%
% RUNFULLWAVETX(icmat) runs a fullwave simulation for the transmit icmat
%
% RUNFULLWAVETX(icmat, simdir) uses the data files in the directory simdir
% and outputs the transmit to a random filename in simdir
%
% RUNFULLWAVETX(icmat, simdir, outdir) outputs the data files to the
% directory simdir
%
% See also READFULLWAVESIM

% TODO: use a time based filename, not a random one.

% TODO: report success or something from the executable

% defaults
if nargin < 2, simdir = pwd; end
if nargin < 3, outdir = tempname(simdir); end
ogdir = pwd;

% get full path of simulation directory
simdir = getfield(what(simdir), 'path');

% move to (temporary) simulation directory on the local machine
tmp = tempname(); 
mkdir(tmp);

try % protect against not cleaning-up errors
    % ensure that the output directory is in fact a directory
    copyfile(fullfile(simdir, 'fullwave2_executable'), tmp); % copy over the executable

    % make the output directory
    mkdir(fullfile(outdir));

    % link to all input (.dat) files
    listing = dir(simdir);
    dfiles = string({listing.name});
    dfiles = dfiles(endsWith(dfiles, '.dat')); % all .dat files
    dfilenames = fullfile(simdir, dfiles); % full paths
    err = system(join("ln -s " + dfilenames + " " + tmp + " ; ")); % symbolic link
    if err % link didn't work: copy everything?
        warning('Unable to link files: copying instead.');
        arrayfun(@copyfile, dfilenames); % copy everything
    end

    % move to the tmp directory
    cd(tmp);

    % write the data to file in the tmp directory
    fid = fopen('icmat.dat','w+');
    fwrite(fid, real(icmat), 'float');
    fclose(fid);

    % call exectuable (blocking)
    disp("Running from " + pwd());
    [err, res] = system('./fullwave2_executable'); % TODO: not sure how to call from Windows/Mac?

    disp("Received code " + err)
    disp(res)

    % return to original directory
    cd(ogdir);

    % move output files to simulation output directory 
    copyfile(fullfile(tmp, 'genout.dat'), fullfile(outdir));

    % TODO: get memory map for the data if an output was requested?

catch ME, getReport(ME) % show the error
    
end

% ensure back in the original folder
cd(ogdir);

% force delete the simulation folder
suc = rmdir(tmp, 's'); % cleanup: recursive deletion
