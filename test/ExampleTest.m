classdef ExampleTest < matlab.unittest.TestCase

    % files
    properties
        run_file (1,1) logical  = true; % set false to only write to file
        delete_file (1,1) logical = true; % whether to delete the file after running
        scat_lim (1,1) double = 1e1; % limit of scatterers
        pulse_lim (1,1) double = 1; % max number of pulses to sim in k-Wave
    end
    
    % blacklists
    properties
        bl_file = ExampleTest.generateBlacklistFiles()
        bl_var = [ ...
    "my_VSX_data.mat", "Trans", "RcvData", "Resource", ... Verasonics
    ] + (pnc() | whitespacePattern);

        % TODO: filter if package unavailable
        bl_fcn = ExampleTest.generateBlacklistFunctions()
    end

    % ------------------------------------------------------------ %
    methods(Static)
        function fls = get_proj_files()
            prj = matlab.project.rootProject; % assume the project is loaded
            % project files
            % if ~isempty(matlab.project.rootProject) && (currentProject().Name == "qups")
            %     prj = currentProject();
            % elseif exist("Qups.prj", "file")
            %     prj = openProject(which("Qups.prj"));
            % else
            %     prj = struct.empty;
            % end
            % TODO: access via project fixture?

            % get all the files in the repo
            if isempty(prj) || prj.Name ~= "qups"
                base_dir = fullfile(fileparts(mfilename('fullpath')), '..');
                list = dir(fullfile(base_dir, '**','*.m'));
                fls = fullfile(string({list.folder}), string({list.name}));
            else
                fls = [prj.Files.Path];
                fls = fls(endsWith(fls, '.m')); % m-files only
            end

            % filter through the file blacklist
            fls = fls(~contains(fls, filesep + ("test"|"build")  + filesep + wildcardPattern() + ".m")); % exclude test,build folder m-files
            fls = fls(~endsWith(fls,ExampleTest.generateBlacklistFiles() + optionalPattern(".m"))); % exclude blacklist
            fls = cellstr(fls);
        end
        function bl_fcn = generateBlacklistFunctions()
            bl_fcn = reshape({ ...
                ... "bfAdjoint", ... QUPS (RAM)
                'kWave',    "kspaceFirstOrder"+optionalPattern(digitsPattern(1)+"D"), ... k-Wave (still too large)
                'USTB',     ["QUPS2USTB", "uff." + alphanumericsPattern + wildcardPattern], ... USTB
                'FieldII',  "calc_scat"+["","_all", "_multi"], ... FieldII
                'MUST',     "simus", ... MUST
                'fullwave', ["getFullwaveTransducer", "fullwaveSim", "fullwaveConf", "fullwaveJob", "mapToCoords"], ... % fullwave
                'mex',      ["bfEikonal", "msfm"+["","2d","3d"]], ... mex binaries required ... compilation (optional, requires CUDA)
                'gpu',            "wbilerpg", ... CUDAKernel or oclKernel support required
                'comp',     ("recompile" + ["","Mex","CUDA"]) ... compilations setup required
                }, 2, [])'; % + "(";
            
            % filter by installed projects
            prj = matlab.project.rootProject;
            i = ismember(string(bl_fcn(:,1)), [cat(2,prj.ProjectReferences.Project).Name]);
            bl_fcn(i,:) = []; % Remove from blacklist if the project is installed.

            % filter by available binaries
            kerns = ("wbilerp")+".ptx"; % require GPU binaries
            if all(arrayfun(@(k) exist(fullfile(prj.RootFolder, "bin", k), 'file'), kerns))
                bl_fcn(bl_fcn(:,1) == "gpu",:) = []; 
            end
            kerns = ["msfm2d", "msfm3d"]+"."+mexext; % require mex binaries
            if all(arrayfun(@(k) exist(fullfile(prj.RootFolder, "bin", k), 'file'), kerns))
                bl_fcn(bl_fcn(:,1) == "mex",:) = []; 
            end

            % filter by gpu compiler availability
            % (not quite working - seems the environment changes?)
            if ~isempty(argn(2, @system, "which nvcc"))
                % bl_fcn(bl_fcn(:,1) == "comp",:) = [];
            end
        end
        function bl_file = generateBlacklistFiles()
            bl_file = [ ...
                "import_verasonics_data", "cheat_sheet", ...
                ("msfm" + ["","2d","3d"]), ...
                "setup", "teardown", ...
                ];
        end

    end
    methods
        function filterBlacklist(test, code, blk, fls)
                % check for blacklisted variables or functions
                for i = length(test.bl_fcn(:,1)):-1:1
                    m(:,i) = contains(code, test.bl_fcn{i,2} + "("); % lines x pcks
                end
                m   = [m    , contains(code, test.bl_var)]; % lines x (pcks+1)
                pkg = [string(test.bl_fcn(:,1))', "var"]; % packages
                l = any(m,2); % matching lines
                p = any(m,1); % matching packages
                test.assumeFalse(any(m,'all'), ...
                    join(rmmissing([ ...
                    "[blacklist] Filtering lines " + join(string(blk([1 end])),"-") + " in file " + fls + "." ...
                        ,"Matching line(s):" + newline + join(code(l), newline) ...
                        ,"Matching packages: " + join(string(pkg(p)), ", ") ...
                    ]), newline));
        end
    end
    methods (TestClassSetup)
        function move2TmpFolder(testCase)
            % folder for temp files
            import matlab.unittest.fixtures.TemporaryFolderFixture;
            import matlab.unittest.fixtures.CurrentFolderFixture;

            % Create a temporary folder and make it the current working folder.
            tempFolder = testCase.applyFixture(TemporaryFolderFixture);
            testCase.applyFixture(CurrentFolderFixture(tempFolder.Folder));
        end
        function hideFigures(testCase)
            % turn off figures during these tests
            state = get(0, 'defaultFigureVisible');
            testCase.addTeardown(@()set(0, 'defaultFigureVisible', state))
            set(0, 'defaultFigureVisible', 'off');
        end
        function silenceAcceptableWarnings(testCase)
            W = warning(); % get current state
            testCase.addTeardown(@() warning(W)); % restore on exit
            warning('off', 'MATLAB:ver:ProductNameDeprecated'); % in k-Wave
            warning('off', "QUPS:kspaceFirstOrder:upsampling"); % in US
        end
        function limitSims(test)
            if gpuDeviceCount()
                test.pulse_lim = 25; else, test.pulse_lim = 3; 
            end
            if gpuDeviceCount() || (exist('oclDeviceCount', 'file') && oclDeviceCount())
                test.scat_lim = 1e4; else, test.scat_lim = 1e2;
            end
        end
    end
    methods (TestClassTeardown)
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
        lns (1,:) cell
        % blk (1,:) cell
    end

    methods(Static, TestParameterDefinition)    
        function [fls, lns] = findExamples()
            % et files
            fls_ = string(ExampleTest.get_proj_files());

            % extract any and all examples from each file
            ws0 = whitespacePattern(0,Inf); % alias

            % start/end pattern
            spat = optionalPattern("%")+ws0+"Example"+ws0+digitsPattern(0,Inf)+ws0+optionalPattern(":") + ws0 + lineBoundary; % starting pat
            epat = "See also" + optionalPattern(":") + ws0; % (optional) ending pat
            pt = [spat epat];

            % comment block pattern
            cpat = ws0 + "%" + wildcardPattern;% + lineBoundary('end');

            % init
            [blk, lns, fls] = deal(cell(1,0));

            for m = 1 : numel(fls_)
                h = false;

                % read in code
                txt = readlines(fls_(m));

                % split into continuous comment blocks
                cb = matches(txt, cpat);
                s = cb; s(2:end  ) = ~s(1:end-1) &  s(2:end); s = find(s);
                e = cb; e(1:end-1) =  e(1:end-1) & ~e(2:end); e = find(e);
                i = arrayfun(@colon, s, e, 'UniformOutput',false);
                i = i(cellfun(@length, i) > 1); % must be more than a one-line comment at least ...

                % for each comment block
                for n = 1:numel(i)
                    code_ = txt(i{n}); % block of text

                    % extract an example if it exists
                    ord = ["first", "last"];
                    for j = 1:2
                        bdln = find(contains(code_, pt(j),'IgnoreCase',false));
                        if isempty(bdln)
                            bdln = nan;
                        elseif ~isscalar(bdln)
                            warning("Multiple lines match " + string(pt(j)) + " in file " + fls_(m) + ": choosing the "+ord(j)+".");
                            switch j, case 1, bdln = bdln(1); case 2, bdln = bdln(end); end
                        end
                        k(j) = bdln; %#ok<AGROW>
                    end
                    k = k + [1, -1]; % go after/before matching start/end
                    if isnan(k(1))
                        % warning("No example found for lines "+join(string(blk{:}([1 end])),"-")+" in file "+f+".");
                        continue;
                    end
                    if isnan(k(2)), k(2) = length(code_); end % can't find an ending - assume it's okay

                    % save
                    % code(end+1) = {code_(k(1):k(2))};
                    fls(end+1) = cellstr(fls_(m));
                    lns(end+1) = {join(string(i{n}(k)),"-")}; % for labelling
                    blk(end+1) = {i{n}(k(1):k(end))};
                    h = true;
                end

                if ~h,
                    fls(end+1) = cellstr(fls_(m));
                    lns(end+1) = {join(string([1 0]),"-")}; % for labelling
                    blk(end+1) = {[1 0]};
                end

            end
        end
    end

    % ---------------------------------------------- % 

    methods (Test, ParameterCombination="sequential", TestTags=["Github", "full", "build"])
        % Test methods
        function run_examples(test, fls, lns)
            % arguments
            %     test matlab.unittest.TestCase
            %     fls (1,1) string
            % end

            % name of the file
            fls = string(fls);
            [~, n] = fileparts(fls);

            % get the code snippet
            lns = double(split(lns,"-"));
            blk = lns(1):lns(end);
            txt = readlines(fls);
            code = txt(blk);

            % if it's empty, there's nothing there ...
            test.assumeNotEmpty(code, ...
                "No examples found in file " + fls + "." ...
                );

            % strip the (first) comment character
            ws0 = whitespacePattern(0,Inf); % alias
            code = extractAfter(code, lineBoundary('start')+ws0+"%");

            % check for blacklisted variables or functions
            filterBlacklist(test, code, blk, fls);

            % assume less than a limit with simulator functions
            args = { 2, "numScat","greens";  1, "seq.numPulse","kspaceFirstOrder"};
            L = [test.scat_lim, test.pulse_lim];
            for l = 1:size(args,1)
                S = check_num_prop(code, args{l,:}); % size
                test.assumeLessThanOrEqual(S, L(l), ...
                    "Skipping example " + n + " <" ...
                    + blk((1)) + "-" + blk((end)) + "> : " +args{l,2}+"=" ...
                    + S + " exceeds the limit of " ...
                    + L(l) + "." ...
                    );
            end

            % make into a function
            fnm = join([n, blk([1 end])],"_"); % function/file name
            header = "function " + fnm + "()";

            % copy into an (output) file
            ofl = fnm + ".m";
            localwritelines([header; code], ofl);

            % delete on cleanup
            if test.delete_file, test.addTeardown(@delete, fullfile(pwd, ofl)); end

            % assert a clean run
            test.assertWarningFree(str2func("@"+fnm), "Example "+fnm+" did not complete without a warning!");

        end
    end
    methods
    end
end

function p = pnc(), p = ("." | "(" | ")" | "{" | "}" | "," | "=" | "*" | "+"); end % punctuation

function V = check_num_prop(code, arg, prp, fnm)
% check that less than N scatterers
% arg = 1; % is us
% prp = "numPulse"; % property
% fnm = "kspaceFirstOrder";
pat = fnm + "(" + asManyOfPattern(alphanumericsPattern() + ",", arg(1)-1,arg(end)-1);
i = find(contains(code,pat)); % lines with call to kspaceFirstOrder

if isempty(i), V = 0; return; end % 0 if greens not found
if ~isscalar(i), i = i(1); warning("Multiple calls to "+fnm+"(): choosing the first."); end
v = strip(extractBetween(code(i), pat, ("," | ")"))); % us variable

% write code to a file
fl = tempname + ".m";
localwritelines(code(1:i-1), fl);

% run the code up to right before calling greens (assumes no loops/branches)
try run(fl);
    % evaluate number of scatterers
    V = eval(v+"."+prp);
    delete(fl); % cleanup

catch ME % failed to run
    disp(ME);
    delete(fl); % cleanup
    V = Inf; %
end
end


function fid = localwritelines(txt, fl)
if exist('writelines','file') % 2022a+
    writelines(txt, fl);
else
    fid = fopen(fl, 'w+');
    fwrite(fid, join(txt,newline));
    fclose(fid);
end
end
