classdef SimTest < matlab.unittest.TestCase
    % SimTest - Simulation test class
    %
    % This class test that all simulation routines properly simulate an
    % array of point targets
    
    %#ok<*PROPLC> object properties can have the same name as local variables
    %#ok<*INUSD> unused inputs
    properties
        us % UltrasoundSystem
        targ % Target
        scanc % Image Scan (Cartesian)
        tscan % Medium Scan
        med % Medium
        clu % processing cluster
    end

    properties(ClassSetupParameter)
        gpu   = getgpu()
        clp = poolfilter({'none', 'threads', 'pool', 'background', 'local'}); % cluster: local has to create a temp cluster each time
        xdc_seq_name = struct(...
            'linfsa', string({"L11-5V", 'FSA'}),...
            'linpw', string({"L11-5V",'Plane-wave'}), ...
            'linvs', string({"L11-5V",'Focused'}), ...
            'crvfsa', string({"C5-2V",'FSA'}), ...
            'crvvs', string({"C5-2V",'sector'}) ...
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
                case {"local", "pool"}, test.assumeTrue(ismember('local', parallel.clusterProfiles()), ...
                        'No local cluster profile available on this platform.');
            end
            

            % setup the parcluster/parpool
            hcp = gcp('nocreate'); % find the active pool
            if ~isempty(hcp)
                switch clp
                    case {"none", "local"}, del = isvalid(hcp); % delete if anything valid/active
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
                case "local",               test.clu = parcluster('local');
                case "threads",     if ecp, test.clu = parpool('threads'); end
                case "background",  if ecp, test.clu = backgroundPool(); end
                case "pool",        if ecp, test.clu = parpool('local', 'SpmdEnabled',true); end % default
            end
        end

        function setupQUPSdata(test, xdc_seq_name)
            % create point target data with the configuration

            % simple point target 30 mm depth
            target_depth = 30e-3;
            targ = Target('pos', [0;0;target_depth], 'c0', 1500); 
            targ.rho_scat = 2; % make density scatterers at 2x the density
            targ.scat_mode = 'ratio';

            % Choose a transducer
            xdc_name = xdc_seq_name(1);
            switch xdc_name
                case 'L11-5V', xdc = TransducerArray.L11_5V(); % linear array
                case 'L12-3V', xdc = TransducerArray.L12_3V(); % another linear array
                case 'C5-2V' , xdc = TransducerConvex.C5_2V(); % convex array
            end
            xdc.numel = 3; % low number of elements - still an array, but we only check the center

            % Choose a transmit sequence
            seq_name = xdc_seq_name(2);
            if isa(xdc, 'TransducerArray')
                switch seq_name
                    case 'FSA', seq = Sequence('type', 'FSA', 'c0', targ.c0, 'numPulse', xdc.numel); % set a Full Synthetic Aperture (FSA) sequence
                    case 'Plane-wave'
                        [amin, amax, Na] = deal( -25 ,  25 , 3 );
                        seq = SequenceRadial('type', 'PW', ...
                            'ranges', 1, 'angles',  linspace(amin, amax, Na), 'c0', targ.c0); % Plane Wave (PW) sequence
                    case 'Focused'
                        [xmin, xmax, Nx] = deal( -10 ,  10 , 3 );
                        seq = Sequence('type', 'VS', 'c0', targ.c0, ...
                            'focus', [1;0;0] .* 1e-3*linspace(xmin, xmax, Nx) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 2mm
                            ...'focus', [1;0;0] .* 1e-3*(-10 : 0.2 : 10) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 0.2 mm
                            );
                end

            elseif isa(xdc, 'TransducerConvex')
                switch seq_name
                    case 'sector'
                        [amin, amax, Na] = deal( -40 ,  40 , 3 );
                        seq = SequenceRadial('type', 'VS', 'c0', targ.c0, ...
                            'angles', linspace(amin, amax, Na), ...
                            'ranges', norm(xdc.center) + target_depth, 'apex', xdc.center ...
                            ); % sector scan sequence
                    case 'FSA'
                        seq = SequenceRadial('type', 'FSA', 'numPulse', xdc.numel, 'c0', targ.c0 ...
                            , 'apex', xdc.center ... % for convex transducer only!
                            ); % FSA - convex sequence
                end
            end

            % make a cartesian scan
            % set the scan at the edge of the transducer
            pn = xdc.positions(); % element positions
            xb = pn(1,[1,end]); % x-limits are the edge of the aperture
            zb = [-10e-3, 10e-3] + [min(targ.pos(3,:)), max(targ.pos(3,:))]; % z-limits surround the point target

            Npixd = 2^8;
            scanc = ScanCartesian(...
                'x', linspace(xb(1), xb(end), Npixd+1), ...
                'z', linspace(zb(1), zb(end), Npixd+1) ...
                ); % X x Z scan
            scanc.x(end) = [];
            scanc.z(end) = [];

            % For linear transducers only!
            if isa(xdc, 'TransducerArray'),
                scan = scanc; % use the cartesian one
            elseif isa(xdc, 'TransducerConvex') % use with a SequenceRadial!
                scan = ScanPolar('origin', seq.apex, 'a', -40:0.5:40, ...
                    'r', norm(seq.apex) + linspace(zb(1), zb(end), Npixd+1)...
                    ); % R x A scan
                scan.r(end) = [];
            end
            
            % Choose the simulation region (eikonal)
            switch xdc_name
                case "C5-2V", tscan = ScanCartesian('x',linspace(-50e-3, 50e-3, 1+100*2^3), 'z', linspace(-30e-3, 60e-3, 1+90*2^3));
                otherwise,    tscan = ScanCartesian('x',linspace(-20e-3, 20e-3, 1+40 *2^3), 'z', linspace(  -05e-3, 55e-3, 1+60*2^3));
            end

            % create a distributed medium based on the point scatterers
            % (this logic will later be implemented in a class)
            s_rad = max([tscan.dx, tscan.dz]); %, 260e-6); % scatterer radius
            nextdim = @(p) ndims(p) + 1;
            ifun = @(p) any(vecnorm(p - swapdim(targ.pos,2,nextdim(p)),2,1) < s_rad, nextdim(p));
            med = Medium('c0', targ.c0, 'rho0', targ.rho0, 'pertreg', {{ifun, [targ.c0*targ.c_scat, targ.rho0*targ.rho_scat]}});

            % Construct an UltrasoundSystem object, combining all of these properties
            us = UltrasoundSystem('xdc', xdc, 'sequence', seq, 'scan', scan, 'fs', 40e6);

            % Make the transducer impulse tighter
            us.xdc.fc = us.fs / 8; % 5e6, with numeric precision w.r.t. us.fs
            us.xdc.bw_frac = 0.2; % SIMUS errors above 0.2
            us.xdc.impulse = us.xdc.ultrasoundTransducerImpulse(); 
            % us.xdc.impulse = Waveform.Delta(); 
            % us.sequence.pulse = Waveform.Delta(); % adding this may make the
            % length of the signal go to 0, causing problems for interpolators

            % point per wavelength - aim for >2 for a simulation
            ppw = targ.c0/xdc.fc/min([tscan.dx, tscan.dz]);
            if ppw < 2, warning("Fewer than 2 points per wavelength (" + ppw + ")."); end

            % save QUPS objects for this test case
            test.us = us; 
            test.targ = targ;
            test.tscan = tscan; 
            test.med = med;

        end
    end
    methods(TestClassTeardown)
    end

    % some of these options aren't fully supported yet.
    properties(TestParameter)
        terp = {'nearest', 'linear', 'cubic'};  % interpolation
        sim_name = {'Greens', 'FieldII', 'FieldII_multi', 'kWave', 'SIMUS'}% k-Wave just takes a while, SIMUS has seemingly random phase errors
    end
    methods(TestMethodSetup)
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
            [us, targ, med, tscan, clu] = deal(...
                test.us, test.targ, test.med, test.tscan, test.clu ...
                );

            % k-Wave doesn't use interpolation: pass on all but one option
            if sim_name == "kWave" && terp ~= "nearest", return; end
            
            % SIMUS restricts the fractional bandwidth to < 0.2 and 
            % requires sampling frequency be a multiple of 4 x the central
            % frequency
            if sim_name == "SIMUS", us.xdc.bw_frac = 0.2; end

            
            % simulate based on the simulation routine
            opts = {'interp', terp, 'parcluster', clu};
            switch sim_name
                case 'FieldII',       chd = calc_scat_all   (us, targ, [1,1], opts{:}); % use FieldII,
                case 'FieldII_multi', chd = calc_scat_multi (us, targ, [1,1], opts{:}); % use FieldII,
                case 'SIMUS'  ,       chd = simus           (us, targ, 'periods', 1, 'dims', 3, opts{:}); % use MUST: note that we have to use a tone burst or LFM chirp, not seq.pulse
                case 'Greens' ,       chd = greens          (us, targ, [1,1], opts{:});
                case 'kWave',         if(gpuDeviceCount) && (clu == 0 || isa(clu, 'parallel.Cluster')), dtype = 'gpuArray-double'; else, dtype = 'double'; end % data type for k-Wave
                                      chd = kspaceFirstOrder(us, med, tscan, 'CFL_max', 0.5, 'PML', [64 128], 'parcluster', clu, 'PlotSim', false, 'DataCast', dtype); % run locally, and use an FFT friendly PML size
                otherwise, warning('Simulator not recognized'); return;
            end

            % peak should be ~near~ 40us at the center elements for
            % FSA and PW, ~near~ 20us (at the focus) for VS

            % truncate the data if we observed the transmit pulse
            if sim_name == "kWave",
                wind = 2*vecnorm(targ.pos) / targ.c0; % 2-way propagation time - use as a search window size
                buf = 2*us.xdc.impulse.duration ... 2-way impulse
                    + range(us.sequence.delays(us.tx), 'all') ... delays
                    + vecnorm(range(us.xdc.positions,2),2,1) ./ targ.c0 ... cross-talk
                    ... + 10e-6 ... % heuristic for things to calm down
                    ;
                tfilt = chd.t0 + buf + wind/2 < chd.time & chd.time < chd.t0 + buf + 3/2*wind; % 40us window
                tfilt = find(tfilt,1,'first'):find(tfilt,1,'last');
                chd.data = sub(chd.data, tfilt, 1);
                chd.t0   = chd.t0 + (tfilt(1)-1)/chd.fs;
            end
            
            % hilbert needed for more accuracy with spectral methods
            if ismember(sim_name, ["kWave", "SIMUS"])
                chd = hilbert(zeropad(chd, 0,2^nextpow2(chd.T)-chd.T));
            end
            
            % get the peak of the center element / center transmit
            n = median(1:chd.N); % must be odd integer size
            m = median(1:chd.M); % must be odd integer size
            ip = argmax(sub(chd.data, {n,m,1}, [2,3,4]), [], 1); % slice rx/tx/frames
            tau = gather(double(chd.time(ip,median(1:size(chd.time,2)), median(1:size(chd.time,3)), 1)));

            % true peak time - exactly 20us travel time
            switch us.sequence.type
                case {'PW', 'FSA'}, t0 = 40e-6;
                case {'VS'},        t0 = 20e-6;
                otherwise, error("Unrecognized sequence type " + us.sequence.type + ".");
            end
            switch sim_name
                case "kWave", tol = 10*(tscan.dz / targ.c0); % within 10 samples of the true location
                case "SIMUS", tol = 1/us.xdc.fc; % SIMUS does not have calibrated phase: be within 1 wavelength of the true location
                otherwise,    tol = 1/chd.fs; % must be accurate down to the sample
            end
            
            % test
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;

            % temporal offset
            test.assertThat(tau, IsEqualTo(t0, 'Within', AbsoluteTolerance(tol)), sprintf(...
                "Peak of the data (" + tau + "us) is offset from the true peak (" + t0 + "us)." ...
                ));
        end
    
        function greens_devs(test, sim_name, terp)
            % test that device options produce (roughly) the same result

            if sim_name ~= "Greens", return; end % pass if not Greens
           
            
            % gpu-cubic is a different algorithm, so we get very different results
            switch terp, case "cubic", return; end

            % unpack
            [us, targ, clu] = deal(test.us, test.targ, test.clu);

            % test definitions
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.AbsoluteTolerance;

            % simulate based on the simulation routine
            opts = {'interp', terp, 'parcluster', clu};
            switch sim_name
                case 'Greens' ,
                    us.sequence.pulse.fun = @(t) single(t==0); % implicit cast to single type
                    chd0 = gather(greens(us, targ, [1,1], opts{:}, 'device', 0 , 'tall', false)); % reference
                    [xo, to] = deal(chd0.data, chd0.t0);
                    for usetall = [true, false]
                        for dev = [0 -1]
                            chd = doubleT(gather(greens(us, targ, [1,1], opts{:}, 'device', dev , 'tall', usetall)));
                            [x, t] = deal(gather(chd.data), gather(chd.t0));
                            test.assertEqual(x, xo, 'AbsTol', 1e-3, 'RelTol', 1e-3, sprintf(...
                                "The data is different on device " + dev + " and tall set to " + usetall + " for a " + us.sequence.type + " sequence."  ...
                                ));

                            % accurate to 10 nanoseconds
                            test.assertEqual(t, to, 'AbsTol', 1e-9, 'RelTol', 1e-9, sprintf(...
                                "The time axis is different on device " + dev + " and tall set to " + usetall + " for a " + us.sequence.type + " sequence."  ...
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
        case "local",       try test.clu = parcluster('local'); catch, tf(i) = false; end
        case "threads",     try test.clu = parpool('threads'); close(gcp()); catch,  tf(i) = false; end
        case "background",  try test.clu = backgroundPool(); close(gcp()); catch, tf(i) = false; end
        case "pool",        try test.clu = parpool('local', 'SpmdEnabled',true); close(gcp()); catch tf(i) = false; end
    end
end
pnms = pnms(tf); % filter the pools that fail to launch
end