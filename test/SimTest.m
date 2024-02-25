classdef SimTest < matlab.unittest.TestCase
    % SimTest - Simulation test class
    %
    % This class test that all simulation routines properly simulate an
    % array of point targets
    
    %#ok<*PROPLC> object properties can have the same name as local variables
    %#ok<*INUSD> unused inputs
    properties
        us % UltrasoundSystem
        scat % Scatterers
        scanc % Image Scan (Cartesian)
        tscan % Medium Scan
        med % Medium
        clu % processing cluster
    end

    properties(ClassSetupParameter)
        gpu   = getgpu()
        clp = poolfilter({'pool', 'default','none' , 'threads', 'background'}); % cluster: local has to create a temp cluster each time
        xdc_seq_name = struct(...
            'linfsa', string({"L11-5V", 'FSA'}),...
            'linpw', string({"L11-5V",'Plane-wave'}), ...
            'linvs', string({"L11-5V",'Focused'}), ...
            'crvfsa', string({"C5-2V",'FSA'}), ...
            'crvvs', string({"C5-2V",'sector'}), ...
            'matpw', string({"PO192O", 'Plane-wave'}) ...
            );
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
        % Shared setup for the entire test class
        function setupQUPS(test, gpu, clp, xdc_seq_name)
            cd(SimTest.proj_folder); % call setup relative to here
            setup(gpu{:}); % compile what we can
            if ~exist('bin', 'dir'), setup cache; end % recompile and make a cache

            % ensure we can actually start these pools
            switch clp
                case "background", test.assumeTrue(logical(exist('backgroundPool','builtin')), ...
                        'No backgroundPool available on this platform.');
            end
            

            % setup the parcluster/parpool
            hcp = gcp('nocreate'); % find the active pool
            if ~isempty(hcp)
                switch clp
                    case {"none", "default"}, del = isvalid(hcp); % delete if anything valid/active
                    case "threads", del = ~isa(hcp, 'parallel.ThreadPool'); % delete if not a threadpool
                    case "background", del = ~isa(hcp, 'parallel.BackgroundPool'); % delete if not a backgroundPool
                    case "pool", del = ~isa(hcp, 'parallel.ProcessPool'); % delete if not a parpool (process pool)
                end
                % delete only if we gotta
                if del, delete(hcp); end 
            end
            
            % choose the new cluster, creating it if appropriate
            ecp = isempty(gcp('nocreate'));
            switch clp 
                case "none",                test.clu = 0;
                case "default",             test.clu = parcluster();
                case "threads",     if ecp, test.clu = parpool('threads'); end
                case "background",  if ecp, test.clu = backgroundPool(); end
                case "pool",        if ecp, test.clu = parpool('local', 'SpmdEnabled',true); end % default
            end
        end

        function setupQUPSdata(test, xdc_seq_name)
            % create point target data with the configuration

            % simple point target 30 mm depth
            target_depth = 15e-3;
            scat = Scatterers('pos', [0;0;target_depth], 'c0', 1500);
            rho0 = 1000; % ambient density
            rho_scat = 2; % make density scatterers at 2x the density

            % Choose a transducer
            xdc_name = xdc_seq_name(1);
            switch xdc_name
                case 'L11-5V', xdc = TransducerArray.L11_5v(); % linear array
                case 'L12-3V', xdc = TransducerArray.L12_3v(); % another linear array
                case 'C5-2V' , xdc = TransducerConvex.C5_2v(); % convex array
                case 'PO192O', xdc = TransducerMatrix.PO192O(); % matrix array
            end
            
            % reduce number of elements - still an array, but we only check the center
            switch xdc_name
                case 'PO192O', xdc.numd(:) = 3;
                otherwise, xdc.numel = 3;
            end

            % Choose a transmit sequence
            seq_name = xdc_seq_name(2);
            if isa(xdc, 'TransducerArray')
                switch seq_name
                    case 'FSA', seq = Sequence('type', 'FSA', 'c0', scat.c0, 'numPulse', xdc.numel); % set a Full Synthetic Aperture (FSA) sequence
                    case 'Plane-wave'
                        [amin, amax, Na] = deal( -25 ,  25 , 3 );
                        seq = SequenceRadial('type', 'PW', 'c0', scat.c0, ...
                            'ranges', 1, 'angles',  linspace(amin, amax, Na)); % Plane Wave (PW) sequence
                    case 'Focused'
                        [xmin, xmax, Nx] = deal( -10 ,  10 , 3 );
                        seq = Sequence('type', 'FC', 'c0', scat.c0, ...
                            'focus', [1;0;0] .* 1e-3*linspace(xmin, xmax, Nx) + scat.pos ... % translating aperture: depth of 30mm, lateral stride of 2mm
                            ...'focus', [1;0;0] .* 1e-3*(-10 : 0.2 : 10) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 0.2 mm
                            );
                end

            elseif isa(xdc, 'TransducerConvex')
                switch seq_name
                    case 'sector'
                        [amin, amax, Na] = deal( -40 ,  40 , 3 );
                        seq = SequenceRadial('type', 'FC', 'c0', scat.c0, ...
                            'angles', linspace(amin, amax, Na), ...
                            'ranges', norm(xdc.center - scat.pos), 'apex', xdc.center ...
                            ); % sector scan sequence
                    case 'FSA'
                        seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'c0', scat.c0 ...
                            ); % FSA - convex sequence
                end
            elseif isa(xdc, 'TransducerMatrix')
                switch seq_name
                    case 'Plane-wave'
                        [amin, amax, Na] = deal( -25 ,  25 , 3 );
                        seq = SequenceRadial('type', 'PW', 'c0', scat.c0, ...
                            'ranges', 1, 'angles',  linspace(amin, amax, Na)); % Plane Wave (PW) sequence
                end
            end

            % make a cartesian scan
            % set the scan at the edge of the transducer
            pn = xdc.positions(); % element positions
            xb = pn(1,[1,end]); % x-limits are the edge of the aperture
            zb = [-10e-3, 10e-3] + [min(scat.pos(3,:)), max(scat.pos(3,:))]; % z-limits surround the point target

            Npixd = 2^7;
            scanc = ScanCartesian(...
                'x', linspace(xb(1), xb(end), Npixd+1), ...
                'z', linspace(zb(1), zb(end), Npixd+1) ...
                ); % X x Z scan
            scanc.x(end) = [];
            scanc.z(end) = [];

            % choose the scan
            if isa(xdc, 'TransducerConvex') && seq_name == "sector" % sector scan 
                scan = ScanPolar('origin', xdc.center, 'a', -40:0.5:40, ...
                    'r', norm(xdc.center) + linspace(zb(1), zb(end), Npixd+1)...
                    ); % R x A scan
                scan.r(end) = [];
            else
                scan = scanc; % use the cartesian one
            end
            
            % Choose the simulation region (eikonal)
            switch xdc_name
                case "C5-2V",   tscan = ScanCartesian('x',linspace(-5e-3, 5e-3, 1+10*2^3), 'z', linspace(-2e-3, 23e-3, 20*2^3));
                case "PO192O",  tscan = ScanCartesian('x',linspace(-2e-3, 2e-3, 1+04*2^3), 'z', linspace(-2e-3, 18e-3, 1+20*2^3));
                                tscan.y = tscan.x;
                otherwise,      tscan = ScanCartesian('x',linspace(-5e-3, 5e-3, 1+10*2^4), 'z', linspace(-5e-3, 20e-3, 1+25*2^4));
            end

            % create a distributed medium based on the point scatterers
            % (this logic will later be implemented in a class)
            s_rad = max([tscan.dx, tscan.dz]); %, 260e-6); % scatterer radius
            nextdim = @(p) ndims(p) + 1;
            ifun = @(p) any(vecnorm(p - swapdim(scat.pos,2,nextdim(p)),2,1) < s_rad, nextdim(p));
            med = Medium('c0', scat.c0, 'rho0', rho0, 'pertreg', {{ifun, [scat.c0, rho0 * rho_scat]}});

            % Construct an UltrasoundSystem object, combining all of these properties
            us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'scan', scan, 'fs', 40e6);

            % Make the transducer impulse tighter
            % us.xdc.fc = us.fs / 8; % 5e6, with numeric precision w.r.t. us.fs
            us.xdc.bw_frac = 0.2; % SIMUS errors above 0.2
            us.xdc.impulse = us.xdc.xdcImpulse(); 
            % us.xdc.impulse = Waveform.Delta(); 
            % us.seq.pulse = Waveform.Delta(); % adding this may make the
            % length of the signal go to 0, causing problems for interpolators

            % point per wavelength - aim for >2 for a simulation
            ppw = scat.c0/xdc.fc/min([tscan.dx, tscan.dz]);
            if ppw < 2, warning("Fewer than 2 points per wavelength (" + ppw + ")."); end

            % save QUPS objects for this test case
            test.us = us; 
            test.scat = scat;
            test.tscan = tscan; 
            test.med = med;

        end
    end
    methods(TestClassTeardown)
    end

    % some of these options aren't fully supported yet.
    properties(TestParameter)
        terp = {'nearest', 'linear', 'cubic'};  % interpolation
        sim_name = {'Greens', 'FieldII', 'FieldII_multi', 'SIMUS', 'kWave'} % k-Wave just takes a while, SIMUS has difficult to predict phase errors
    end
    methods(TestMethodSetup)
        function resetGPU(test)
            if gpuDeviceCount()
                reselectgpu();
                if ~isempty(gcp('nocreate'))
                    wait(parfevalOnAll(gcp(), @reselectgpu, 0));
                end
            end

            % helper
            function reselectgpu()
                id = gpuDevice().Index; gpuDevice([]); gpuDevice(id); 
            end
        end
    end

    % Github test routine
    methods(Test, ParameterCombination = 'sequential', TestTags={'Github'})
        function github_pscat(test)
            % only test Green's function
            % switch terp, case {'nearest','linear','cubic'}, otherwise, return; end
            % terp = 'cubic';

            % skip local pools - they take too long;
            if isa(test.clu, 'parallel.Cluster') || isa(test.clu, 'parallel.ProcessPool'), return; end
            
            % forward remaining
            pscat(test, 'Greens', 'cubic');
        end
    end

    % Full test routine
    methods(Test, ParameterCombination = 'exhaustive', TestTags={'full'})
        function full_pscat(test, sim_name, terp), pscat(test, sim_name, terp); end % forward all
        function full_greens_devs(test, sim_name, terp), greens_devs(test, sim_name, terp); end % forward all
    end

    
    % test method implementations
    methods
        function pscat(test, sim_name, terp)
            % PSF - Test the PSF
            % 
            % Test that the PSF for a point at a reasonable distance 
            % properly beamforms for all transmit sequences, transducer 
            % types and compute device.

            % unpack
            [us, scat, med, tscan, clu] = deal(...
                test.us, test.scat, test.med, test.tscan, test.clu ...
                );

            % check if we can even call the sim
            switch sim_name
                case {'FieldII', 'FieldII_multi'} 
                                test.assumeTrue(logical(exist('field_init', 'file')));
                case {'SIMUS'}, test.assumeTrue(logical(exist('pfield3'   , 'file')));
                case {'kWave'}, test.assumeTrue(logical(exist('kWaveGrid' , 'file')));
            end

            % k-Wave/calc_scat_multi don't use interpolation: pass on all but one option
            if ismember(sim_name,["FieldII_multi", "kWave"]) && terp ~= "nearest", return; end
            
            % SIMUS restricts the fractional bandwidth to < 0.2 and 
            % requires sampling frequency be a multiple of 4 x the central
            % frequency
            if sim_name == "SIMUS", us.xdc.bw_frac = 0.2; end

            
            % simulate based on the simulation routine
            opts = {'parenv', clu, 'interp', terp, };
            switch sim_name
                case 'FieldII',       chd = calc_scat_all   (us, scat, opts{:}); % use FieldII,
                case 'FieldII_multi', chd = calc_scat_multi (us, scat, opts{[1:2]}); % use FieldII,
                case 'SIMUS'  ,       chd = simus           (us, scat, 'periods', 1, 'dims', 3, opts{[1:2]}); % use MUST: note that we have to use a tone burst or LFM chirp, not seq.pulse
                case 'Greens' ,       chd = greens          (us, scat, opts{[3:4]});
                case 'kWave',         if(gpuDeviceCount && (clu == 0 || isa(clu, 'parallel.Cluster'))), dtype = 'gpuArray-single'; else, dtype = 'single'; end % data type for k-Wave
                                      T = 2.1 * max(vecnorm(scat.pos - us.xdc.positions,2,1)) / med.c0 + us.seq.pulse.tend + 2 * us.xdc.impulse.tend; % signal end time
                                      chd = kspaceFirstOrder(us, med, tscan, 'CFL_max', 0.5, 'PML', [8 64], 'parenv', clu, 'PlotSim', false, 'DataCast', dtype, "T", T); % run locally, and use an FFT friendly PML size
                otherwise, warning('Simulator not recognized'); return;
            end

            % peak should be ~near~ 40us at the center elements for
            % FSA and PW, ~near~ 20us (at the focus) for FC

            % truncate the data if we observed the transmit pulse
            if sim_name == "kWave",
                wind = 2*vecnorm(scat.pos) / scat.c0; % 2-way propagation time - use as a search window size
                buf = 2*us.xdc.impulse.duration ... 2-way impulse
                    + range(us.seq.delays(us.tx), 'all') ... delays
                    + vecnorm(range(us.xdc.positions,2),2,1) ./ scat.c0 ... cross-talk
                    ... + 10e-6 ... % heuristic for things to calm down
                    ;
                tfilt = chd.t0 + buf + wind/2 < chd.time & chd.time < chd.t0 + buf + 3/2*wind; % 40us window
                tfilt = find(tfilt,1,'first'):find(tfilt,1,'last');
                chd.data = sub(chd.data, tfilt, 1);
                chd.t0   = chd.t0 + (tfilt(1)-1)/chd.fs;
            end
            
            % hilbert needed for more accuracy with spectral methods
            if ismember(sim_name, ["kWave", "SIMUS"])
                chd = hilbert(zeropad(chd, 0,2^(1+nextpow2(chd.T))-chd.T));
            end
            
            % get the peak of the center element / center transmit
            n = median(1:chd.N); % must be odd integer size
            m = median(1:chd.M); % must be odd integer size
            ip = argmax(sub(chd.data, {n,m,1}, [2,3,4]), [], 1); % slice rx/tx/frames
            tau = gather(double(chd.time(ip,median(1:size(chd.time,2)), median(1:size(chd.time,3)), 1)));

            % true peak time - exactly 10us travel time
            switch us.seq.type
                case {'PW', 'FSA'}, t0 = 20e-6;
                case {'FC'},        t0 = 10e-6;
                otherwise, error("Unrecognized sequence type " + us.seq.type + ".");
            end
            switch sim_name
                case "kWave", tol = double(10*(tscan.dz / scat.c0)); % within 10 samples of the true location
                case "SIMUS", tol = double(1/us.xdc.fc); % SIMUS does not have calibrated phase: be within 1 wavelength of the true location
                otherwise,    tol = double(1.1/chd.fs); % must be accurate down to the sample
            end
            
            % test
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;

            % temporal offset
            test.assertThat(tau, IsEqualTo(t0, 'Within', AbsoluteTolerance(tol)), sprintf(...
                "Peak of the data (" + tau + "us) is offset by more than " + tol + "from the true peak (" + t0 + "us)." ...
                ));
        end
    
        function greens_devs(test, sim_name, terp)
            % test that device options produce (roughly) the same result

            if sim_name ~= "Greens", return; end % pass if not Greens
           
            
            % gpu-cubic is a different algorithm, so we get very different results
            switch terp, case "cubic", return; end

            % unpack
            [us, scat, clu] = deal(test.us, test.scat, test.clu);

            % test definitions
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.AbsoluteTolerance;

            % simulate based on the simulation routine
            opts = {'interp', terp, 'parenv', clu};
            switch sim_name
                case 'Greens' 
                    us.fs = single(us.fs); % implicit cast to single type
                    chd0 = gather(greens(us, scat, [1,1], opts{1:2}, 'device', 0 , 'tall', false)); % reference
                    [xo, to] = deal(double(chd0.data), double(chd0.t0));
                    for usetall = [true, false]
                        for dev = [0 -1]
                            chd = gather(greens(us, scat, [1,1], opts{1:2}, 'device', dev , 'tall', usetall));
                            [x, t] = deal(double(gather(chd.data)), gather(double(chd.t0)));
                            test.assertEqual(x, xo, 'AbsTol', 1e-3, 'RelTol', 1e-3, sprintf(...
                                "The data is different on device " + dev + " and tall set to " + usetall + " for a " + us.seq.type + " sequence."  ...
                                ));

                            % accurate to 10 nanoseconds
                            test.assertEqual(t, to, 'AbsTol', 1e-9, 'RelTol', 1e-9, sprintf(...
                                "The time axis is different on device " + dev + " and tall set to " + usetall + " for a " + us.seq.type + " sequence."  ...
                                ));
                        end
                    end
            end
        end
    end

    methods(Static)
        % PROJ_FOLDER - Identifies the base folder for the project
        function f = proj_folder(), f = fullfile(fileparts(mfilename('fullpath')), '..'); end

        function addcache() 
            % ADDCACHE - add the cached bin folder
            % adds the cached bin folder - this prevents recompilation for
            % each and every test, which can be a majority of the
            % computational effort / time 
            try addpath(fullfile(SimTest.proj_folder,'bin')); end %#ok<TRYNC>
        end
    end
end

function s = getgpu()
s.no_ptx = {};
if gpuDeviceCount, s.ptx = {'CUDA'}; end
end

function pnms = poolfilter(pnms)
for i = numel(pnms):-1:1
    tf(i) = true;
    switch pnms{i}
        case "none",
        case "default",     try test.clu = parcluster(); catch, tf(i) = false; end
        case "threads",     try test.clu = parpool('threads'); close(gcp()); catch,  tf(i) = false; end
        case "background",  try test.clu = backgroundPool(); close(gcp()); catch, tf(i) = false; end
        case "pool",        try test.clu = parpool('local', 'SpmdEnabled',true); close(gcp()); catch, tf(i) = false; end
    end
end
pnms = pnms(tf); % filter the pools that fail to launch
end