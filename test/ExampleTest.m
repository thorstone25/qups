classdef ExampleTest < matlab.unittest.TestCase

    % files
    properties
        run_file (1,1) logical  = true; % set false to only write to file
        delete_file (1,1) logical = true; % whether to delete the file after running
        meta_path (1,1) string = "test/meta"; % path for example scripts (relative to proj)
        base_dir  (1,1) string = fullfile(fileparts(mfilename('fullpath')), '..')
    end
    properties(Dependent)
        meta_dir (1,1) string
    end
    methods
        function m = get.meta_dir(test), m = fullfile(test.base_dir, test.meta_path); end
    end
    
    % blacklists
    properties
        bl_file = [ ...
    "import_verasonics_data", "run_QUPS_tests", "cheat_sheet", ...
    ("msfm" + ["","2d","3d"]), ...
    "setup", "teardown", ...
    ];
        bl_var = [ ...
    "my_VSX_data.mat", "Trans", "RcvData", "Resource", ... Verasonics
    ] + (pnc() | whitespacePattern);

        bl_fcn = [ ...
            ... "bfAdjoint", ... QUPS (RAM)
    "getFullwaveTransducer", "fullwaveSim", "fullwaveConf", "fullwaveJob", "mapToCoords", ... % fullwave
    "kspaceFirstOrder"+optionalPattern(digitsPattern(1)+"D"), ... k-Wave
    "QUPS2USTB", "uff." + alphanumericsPattern + wildcardPattern, ... USTB
    "calc_scat"+["","_all", "_multi"], ... FieldII
    "simus", ... MUST
    "recompile" + ["","Mex","CUDA"], ... compilation (optional, requires CUDA)
    ] + "(";
        bl_gpu = [
            "wbilerpg", ... CUDAKernel or oclKernel support required
            ]
    end
    
    % ------------------------------------------------------------ %
    
    methods (TestClassSetup)
        function setup_metadir(test)
            % folder for temp files
            if ~exist(test.meta_dir,'dir'), mkdir(test.meta_dir); end
            addpath(test.meta_dir);
        end
    end
    methods (TestClassTeardown)
        function teardown_metadir(test)
            % delete folder for temp files
            rmpath(test.meta_dir);
            if test.delete_file, try rmdir(test.meta_dir); end, end %#ok<TRYNC> % don't force deletion
        end
    end

    % ---------------------------------------------- % 
    
    methods (TestMethodSetup)
        % Setup for each test
    end
    methods(TestMethodTeardown)
        function cleanup_test(test)
            if test.run_file, try, gpuDevice([]); end, end % clear the gpu if there
            close all; % close all figures
        end
    end

    properties (TestParameter)
        fls (1,:) cell 
    end

    methods(Static, TestParameterDefinition)
        function fls = get_proj_files()
            % project files
            % if ~isempty(matlab.project.rootProject) && (currentProject().Name == "qups")
            %     prj = currentProject();
            % elseif exist("Qups.prj", "file")
            %     prj = openProject(which("Qups.prj"));
            % else
            %     prj = struct.empty;
            % end
            % TODO: access via project fixture?
            prj = struct.empty;

            % get all the files in the repo
            if isempty(prj)
                base_dir = fullfile(fileparts(mfilename('fullpath')), '..');
                list = dir(fullfile(base_dir, '**','*.m'));
                fls = fullfile(string({list.folder}), string({list.name}));
            else
                fls = [prj.Files.Path];
                fls = fls(endsWith(fls, '.m')); % m-files only
            end

            % filter through the file blacklist
            fls = fls(~contains(fls, filesep + "test" + filesep + "*.m")); % exclude test folder m-files
            fls = cellstr(fls);
        end
    end

    % ---------------------------------------------- % 
   
    methods (Test, TestTags=["Github", "full"])
        % Test methods
        function run_examples(test, fls)
            % arguments
            %     test matlab.unittest.TestCase
            %     fls (1,1) string
            % end
            fls = string(fls);

            % assume it's not on the blacklist or the meta_path
            blf = [test.bl_file, test.meta_path + wildcardPattern];
            if(endsWith(fls,blf + optionalPattern(".m"))) % pass any blacklisted files
                warning("[blacklist] Skipping " + fls); return; 
            end 

            % extract any and all examples from each file
            ws0 = whitespacePattern(0,Inf); % alias
            
            % start/end pattern
            spat = optionalPattern("%")+ws0+"Example"+ws0+digitsPattern(0,Inf)+ws0+optionalPattern(":") + ws0 + lineBoundary; % starting pat
            epat = "See also" + optionalPattern(":") + ws0; % (optional) ending pat
            pt = [spat epat]; 

            % comment block pattern
            cpat = ws0 + "%" + wildcardPattern;% + lineBoundary('end');

            % name of the file
            [~, n] = fileparts(fls);

            % files (kernel functions) that require CUDAKernel or oclKernel support
            test.assumeTrue(~ismember(n, test.bl_gpu) || gpuDeviceCount || (exist('oclDeviceCount', 'file') && oclDeviceCount));

            % read in code
            txt = readlines(fls);

            % split into continuous comment blocks
            cb = matches(txt, cpat);
            s = cb; s(2:end  ) = ~s(1:end-1) &  s(2:end); s = find(s);
            e = cb; e(1:end-1) =  e(1:end-1) & ~e(2:end); e = find(e);
            i = arrayfun(@colon, s, e, 'UniformOutput',false);
            i = i(cellfun(@length, i) > 1); % must be more than a one-line comment at least ...

            % for each comment block
            h = false;
            for blk = i(:)'
                code = txt(blk{:}); % block of text

                % extract an example if it exists
                ord = ["first", "last"];
                for j = 1:2
                    bdln = find(contains(code, pt(j),'IgnoreCase',false));
                    if isempty(bdln)
                        bdln = nan;
                    elseif ~isscalar(bdln)
                        warning("Multiple lines match " + string(pt(j)) + " in file " + fls + ": choosing the "+ord(j)+".");
                        switch j, case 1, bdln = bdln(1); case 2, bdln = bdln(end); end
                    end
                    k(j) = bdln; %#ok<AGROW>
                end
                k = k + [1, -1]; % go after/before matching start/end
                if isnan(k(1))
                    % warning("No example found for lines "+join(string(blk{:}([1 end])),"-")+" in file "+f+".");
                    continue;
                end
                if isnan(k(2)), k(2) = length(code); end % can't find an ending - assume it's okay
                code = code(k(1):k(2));

                % strip the (first) comment character
                code = extractAfter(code, lineBoundary('start')+ws0+"%");

                % check for blacklisted variables or functions
                l = contains(code, test.bl_var);
                m = contains(code, test.bl_fcn);
                if any(l) || any(m)
                    warning("[blacklist] Filtering lines "+join(string(blk{:}([1 end])),"-")+" in file "+fls+"."...
                        +newline+" Matching line(s):" + newline + join(code(l | m), newline));
                    continue;
                end

                % TODO: pass or filter based on Tag
                % assume less than 25 scatterers with the greens function
                S = check_num_scat(code); % number of scatterers
                L = 1e4;
                test.assumeTrue(S <= L, "Example uses " + S + " scatterers, exceeding the limit of " + L + ".");

                % make into a function
                fnm = join([n, blk{1}(k)],"_"); % function/file name
                header = "function " + fnm + "()";

                % copy into an (output) file
                ofl = fullfile(test.meta_dir, fnm +".m");
                lwritelines([header; code], ofl);

                % delete on cleanup
                if test.delete_file, test.addTeardown(@delete, ofl); end

                % assert a clean run
                test.assertWarningFree(str2func("@"+fnm), "Example "+fnm+" did not complete without a warning!");

                % mark at least 1 help file
                h = true;
            end
            if ~h, warning("No examples found in file " + fls + "."); end
        end
    end
end

function p = pnc(), p = ("." | "(" | ")" | "{" | "}" | "," | "=" | "*" | "+"); end % punctuation

function S = check_num_scat(code)
% check that less than N scatterers
pat = "greens(" + alphanumericsPattern() + ",";
i = find(contains(code,pat)); % lines with call to greens

if isempty(i), S = 0; return; end % 0 if greens not found
if ~isscalar(i), i = i(1); warning("Multiple calls to greens: choosing the first."); end
v = strip(extractBetween(code(i), pat, ("," | ")"))); % scat variable

% write code to a file
fl = tempname + ".m";
if exist('writelines','file')
    lwritelines(code(1:i-1), fl);
else
end

% run the code up to right before calling greens (assumes no loops/branches)
try run(fl);
    % evaluate number of scatterers
    S = eval(v+".numScat");
    delete(fl); % cleanup

catch ME % failed to run
    disp(ME);
    delete(fl); % cleanup
    S = Inf; %
end
end

function fid = lwritelines(txt, fl)
if exist('writelines','file') % 2022a+
    writelines(txt, fl);
else
    fid = fopen(fl, 'w+');
    fwrite(fid, join(txt,newline));
    fclose(fid);
end
end
