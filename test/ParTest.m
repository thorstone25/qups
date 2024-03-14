classdef ParTest < matlab.unittest.TestCase
    properties
        us
        sct
        md
        grd
        dif
        cd
        sq
        tlimit = 60; % if expected to exceed this many seconds, don't continue
    end

    methods(TestClassSetup)
        function setupData(test)
            % setup data

            % system
            us_ = UltrasoundSystem( ...
                'seq', SequenceRadial('angles',[-10 0 10], 'c0', 1500), ...
                'xdc', TransducerArray('numel',9), ...
                'recompile', false);
            us_.fs = single(4*us_.xdc.fc);
            us_ = scale(us_, 'dist',1e3, 'time',1e6);

            % Media
            [db_spek, db_scat] = deal(-40, 20);
            [c0, rho0] = deal(us_.seq.c0, us_.seq.c0 / 1.5);
            pb = us_.xdc.bounds();
            grid = ScanCartesian('x', 10*pb(1,:), 'y',0, 'z', [1, 51]*us_.lambda);
            [grid.dx, grid.dz] = deal(us_.lambda / 2.5);
            scat_1 = Scatterers('pos', mean(cat(1,grid.xb,grid.yb,grid.zb),2), 'c0', c0); % center of grid
            % vsz = prod(range([grid.xb; grid.zb],2) ./ us_.lambda); % voxel size (resolution)
            S = 2 .^ (0 : 3 : 27); % round(max(1,100*vsz))], % scale up to N scats per voxel
            scat_d = arrayfun(@(S) Scatterers.Diffuse(grid, S, db_scat, 'c0', c0), S); 
            ONE = ones(grid.size);
            c   =   c0 * ONE;
            rho = rho0 * ONE + db2mag(db_spek) * randn(grid.size);
            med = Medium.Sampled(grid, c, rho, 'c0', c0, 'rho0', rho0);

            % Imaging
            us_.scan = copy(grid);
            [us_.scan.dx, us_.scan.dz] = deal(us_.lambda / 4);

            % ChannelData
            chd = singleT(gather(greens(us_, scat_1)));

            % k-Wave parallelized sequence
            try M = str2double(argn(2, @system, 'nproc')) / 2;
            catch ME, M = 16; end % guess at how many threads to test
            seq = SequenceRadial('angles',linspace(-10, 10, M), 'c0', us_.seq.c0);

            % save
            test.us = us_;
            test.md = med;
            test.grd = grid;
            test.sct = scat_1; % single
            test.dif = scat_d; % diffuse
            test.cd = chd; % data
            test.sq = seq; % alternate sequence
        end

        function configurePenv(test)
            hcp = gcp('nocreate'); % get current pool type

            % teardown to restore pool
            if isempty(hcp)
                % pass
            elseif isa(hcp, "parallel.ThreadPool")
                test.addTeardown(@() parpool("Threads"));
            elseif isa(hcp, "parallel.ProcessPool")
                test.addTeardown(@() parpool(hcp.Cluster));
            elseif isa(hcp, "parallel.ClusterPool")
                test.addTeardown(@() parpool(hcp.Cluster));
            else
                warning("Unable to restore pool on completion");
            end

            % teardown to delete current pool (if it exists then)
            test.addTeardown(@() delete(gcp('nocreate')));

            % actually delete the current pool
            delete(hcp);
        end
    end

    % ------------------------------------------------------------ %
    % files
    properties(TestParameter)
        par_prof
        dev
    end

    methods(Static, TestParameterDefinition)
        function dev = getdevs()
            dev = {0};
            if gpuDeviceCount
                dev = [dev, arrayfun(@gpuDevice, 1:gpuDeviceCount(), 'UniformOutput', false)];
            end
            if exist('oclDeviceCount', 'file') && oclDeviceCount()
                dev = [dev, arrayfun(@oclDevice, 1:oclDeviceCount(), 'UniformOutput', false)];
            end
            [~, i] = unique(cellfun(@(d) string(d.Name)+"-"+class(d), dev(2:end)), 'stable');
            dev = [dev(1), dev(1+i)]; % uniquely named devices only
        end
        function par_prof = getParProf()
            par_prof = ["none", parallel.clusterProfiles];
            if exist('parallel.BackgroundPool', 'class'), par_prof(end+1) = "background"; end
            par_prof = cellstr(par_prof);
        end
    end
    % ---------------------------------------------- %

    methods
        function setDev(test, dev)
            import matlab.unittest.fixtures.TemporaryFolderFixture;
            import matlab.unittest.fixtures.CurrentFolderFixture;

            % Create a temporary folder and make it the current working folder.
            tempFolder = test.applyFixture(TemporaryFolderFixture);
            test.applyFixture(CurrentFolderFixture(tempFolder.Folder));

            % The test can now write files to the current working folder.
            switch class(dev)
                case "double",                  dtype = "none"; 
                case "parallel.gpu.CUDADevice", dtype = "gpu";  gpuDevice(dev.Index);
                case "oclDevice",               dtype = "ocl";  oclDevice(dev.Index);
            end
            odevs = setdiff(["gpu", "ocl"], dtype);
            for odev = odevs
                localwritelines("function n = "+odev+"DeviceCount(), n = 0;", odev+"DeviceCount.m");
                test.addTeardown(@()delete(fullfile(pwd, odev+"DeviceCount.m")));
            end
        end

        function setupPenv(test, par_prof)
            % TODO: allow testing for each parallel.clusterProfiles
            par_prof = string(par_prof);
            hcp = gcp('nocreate'); % active pool
            switch par_prof
                case "none",        oldT = "";
                case "background",  oldT = "";
                case "local",       oldT = "parallel.ProcessPool";
                case "Processes",   oldT = "parallel.ProcessPool";
                case "Threads",     oldT = "parallel.ThreadPool";
                otherwise,          oldT = "";
            end

            if (~isempty(hcp) && ~isa(hcp, oldT)) % active pool matches
            elseif isempty(hcp) && ~strlength(oldT) % no active pool, none needed
            else % mismatch
                disp("Switching from " + class(hcp) + " to a " + par_prof + " profile.");
                if ~isempty(hcp), delete(hcp); end
                switch par_prof
                    case "none"
                    case "background"
                    case "local",     parpool("local"       );
                    case "Processes", parpool("Processes"   );
                    case "Threads",   parpool("Threads"     );
                    otherwise,        % parpool(parcluster()  ); % TODO: flag for parpool vs. parcluster input
                end
            end
        end
    end

    % ---------------------------------------------- %
    methods
        function penv = getPenv(test, par_prof)
            % get the penv
            par_prof = string(par_prof);
            switch par_prof
                case "none",        penv = 0;
                case "background",  penv = backgroundPool();
                case "local",       penv = gcp();
                case "Processes",   penv = gcp();
                case "Threads",     penv = gcp();
                otherwise,          penv = parcluster(par_prof);
            end
        end

        function logTestCheck(tst, chd, tt)
            tst.log("Completed in " + 1e3*toc(tt) + " milliseconds.");
            tst.assertFalse( anynan(chd.data) , "NaN values were detected." )
            tst.assertFalse(~logical(nnz(chd.data)), "All values are 0.");
        end
    end

    methods (Test, TestTags=["full","build"], ParameterCombination='sequential')
        % test simulators can run on different devices
        function greens_das_dev(tst, dev)
            % test that greens works
            
            tt = tic(); chd = greens(tst.us, tst.sct, 'device', -1*(dev ~= 0));
            tst.logTestCheck(chd, tt);

            % test that DAS works
            tt = tic(); chd.data = DAS(tst.us, chd, 'device', -1*(dev ~= 0));
            tst.logTestCheck(chd, tt);
        end
        function greens_penv(tst, par_prof) % test greens
            tst.setup(par_prof);
            tt = tic(); chd = greens(tst.us, tst.sct, 'parenv', tst.getPenv(par_prof));
            tst.logTestCheck(chd, tt);
        end
        function fieldII_penv(tst, par_prof) % test FieldII
            tst.setup(par_prof);
            penv = tst.getPenv(par_prof); % get penv for this profile
            if contains(par_prof, ["Threads", "background"]), return; end % incompatible
            tt = tic(); chd = calc_scat_multi(tst.us, tst.sct, 'parenv', penv);
            tst.logTestCheck(chd, tt);
        end
        function kwave_penv(tst, par_prof) % test k-Wave
            tst.setup(par_prof);
            tt = tic(); chd = kspaceFirstOrder(tst.us, tst.md, tst.grd, 'parenv', tst.getPenv(par_prof));
            tst.logTestCheck(chd, tt);
        end
        function simus_penv(tst, par_prof) % test simus
            tst.setup(par_prof);
            [us_, sct_] = dealfun(@(x) scale(x,'dist',1e-3,'time',1e-6), tst.us, tst.sct); % return to SI units - MUST is hard-coded
            tt = tic(); chd = simus(us_, sct_(1), 'parenv', tst.getPenv(par_prof), 'dims',2);
            tst.logTestCheck(chd, tt);
        end
    end

    methods(Test, TestTags="benchmark",  ParameterCombination='sequential')
        function greens_benchmark(tst, dev)
            if ~isnumeric(dev) || dev ~= 0, tst.setDev(dev); else, return; end % pass on no dev
            f = @greens;
            dur = nan(1,numel(tst.dif));
            S   = [tst.dif.numScat, Inf];
            for i = 1:numel(tst.dif)
                tst.log("Benchmarking " + func2str(f) + " method on device " + dev.Name + " with " + tst.dif(i).numScat + " scatterers.");
                tt = tic(); chd = f(tst.us, tst.dif(i), 'device', -1); tst.logTestCheck(chd, tt); dur(i) = toc(tt); toc(tt)
                if dur(i) * S(i+1)/S(i) > tst.tlimit, break; end % don't push it on the scaling ...
            end
            msg = "Method "+func2str(f)+" on "+string(class(dev))+": " + 1e6*dur(i)/S(i) + " microseconds / scatterer ("+S(i)+").";
            disp(msg); tst.log(1, join(msg,newline));
        end
        function fieldII_benchmark(tst, par_prof)
            if ~contains(par_prof, ["Threads", "background"]), tst.setupPenv(par_prof); else, return; end
            f = @calc_scat_multi;
            dur = nan(1,numel(tst.dif));
            S   = [tst.dif.numScat, Inf];
            for i = 1:numel(tst.dif) % last is dev only for now
                tst.log("Benchmarking " + func2str(f) + " method on profile " + par_prof + " with " + tst.dif(i).numScat + " scatterers.");
                tt = tic(); chd = f(tst.us, tst.dif(i), 'parenv', tst.getPenv(par_prof)); tst.logTestCheck(chd, tt); dur(i) = toc(tt); toc(tt)
                if dur(i) * S(i+1)/S(i) > tst.tlimit, break; end % don't push it on the scaling ...
            end            
            msg = "Method "+ func2str(f) +" on profile "+par_prof+": " + 1e6*dur(i)/S(i) + " microseconds / scatterer ("+S(i)+").";
            disp(msg); tst.log(1, msg);
        end
        function kWave_benchmark(tst, par_prof)
            if ~contains(par_prof, "background"), tst.setupPenv(par_prof); else, return; end
            f = @kspaceFirstOrder;
            us_ = copy(tst.us); us_.seq = tst.sq; 
            t0 = 2*us_.xdc.impulse.t0 + us_.seq.pulse.t0 + min(us_.seq.delays(us_.xdc),[],'all');
            T = 1; dur = 0; stp = Inf; grid = copy(tst.grd); PML = sort(floor(max(0, 225 - grid.size(1:2)) / 2)) + [-1 1];
            scl = 2^3; % amount to increase each simulation
            tst.log("Benchmarking " + func2str(f) + " method on profile " + par_prof + " with grid size [" + join(string(grid.size),",") + "].")
            while(dur*scl < tst.tlimit)
                tt = tic(); chd = f(us_, tst.md, grid, 'PML', PML, 'T', t0 + (T / us_.fs), 'parenv', tst.getPenv(par_prof)); 
                tst.logTestCheck(chd, tt); dur = toc(tt); toc(tt)
                stp = min(stp, dur); % log setup time as shortest time to run
                dur = max(eps, dur - stp); % subtract pre-processing time
                T = T * min(scl, tst.tlimit ./ dur); % increase time steps
            end
            msg = "Method "+ func2str(f) +" on profile "+par_prof+": " + 1e3*dur*us_.seq.numPulse/chd.T ...
                + " milliseconds / time step (grid = ["...
                + join(string([225 225 1]),",") + ...
            "]).";
            disp(msg); tst.log(1, msg);
        end
    end
end


function localwritelines(txt, fl)
if exist('writelines','file') % 2022a+
    writelines(txt, fl);
else
    fid = fopen(fl, 'w+');
    fwrite(fid, join(txt,newline));
    fclose(fid);
end
end
