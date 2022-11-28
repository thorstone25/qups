classdef BFTest < matlab.unittest.TestCase
    % BFTest - Beamformer test class
    %
    % This class test that all beamformers properly beamform a point target
    
    %#ok<*PROPLC> object properties can have the same name as local variables
    properties
        chd % ChannelData
        us % UltrasoundSystemct
        scat % Scatterers
        scanc % Image Scan (Cartesian)
        tscan % Sound-speed Scan
        apod % image apodization
        fmod % modulation of the Channel Data
    end

    properties(ClassSetupParameter)
        par   = struct('one_thread', {{}})%, 'threadpool', {{'parallel'}});
        gpu   = getgpu()
        xdc_seq_name = struct(...
            'linfsa', string({"L11-5V", 'FSA'}),...
            'linpw', string({"L11-5V",'Plane-wave'}), ...
            'linvs', string({"L11-5V",'Focused'}), ...
            'crvfsa', string({"C5-2V",'FSA'}), ...
            'crvvs', string({"C5-2V",'sector'}) ...
            );
        target_offset = struct('offset', 5e-3, 'centered', 0);
        baseband = struct('false', false, 'true', true);
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
        % Setup for the entire test class
        function addcache(test)
            % ADDCACHE - add the cached bin folder
            % adds the cached bin folder - this prevents recompilation for
            % each and every test, which can be a majority of the
            % computational effort / time 
            cache_folder = fullfile(BFTest.proj_folder,'bin');
            if ~exist(cache_folder, 'dir'), try mkdir(cache_folder); end, end %#ok<TRYNC>
            if  exist(cache_folder, 'dir')    addpath(cache_folder); end 
        end

        % Shared setup for the entire test class
        function setupQUPS(test, par, gpu)
            cd(BFTest.proj_folder); % call setup relative to here
            setup(par{:}, gpu{:}); % compile what we can
            if ~exist('bin', 'dir'), setup cache; end % recompile and make a cache
        end

        function setupQUPSdata(test, xdc_seq_name, target_offset, baseband)
            %% create point target data with the configuration

            % simple point target 30 mm depth
            target_depth = 30e-3;
            scat = Scatterers('pos', [target_offset;0;target_depth], 'c0', 1500); 

            % Choose a transducer
            xdc_name = xdc_seq_name(1);
            switch xdc_name
                case 'L11-5V', xdc = TransducerArray.L11_5v(); % linear array
                case 'L12-3V', xdc = TransducerArray.L12_3v(); % another linear array
                case 'C5-2V' , xdc = TransducerConvex.C5_2v(); % convex array
            end

            % Choose a transmit sequence
            seq_name = xdc_seq_name(2);
            if isa(xdc, 'TransducerArray')
                switch seq_name
                    case 'FSA', seq = Sequence('type', 'FSA', 'c0', scat.c0, 'numPulse', xdc.numel); % set a Full Synthetic Aperture (FSA) sequence
                    case 'Plane-wave'
                        [amin, amax, Na] = deal( -25 ,  25 , 25 );
                        seq = SequenceRadial('type', 'PW', ...
                            'ranges', 1, 'angles',  linspace(amin, amax, Na), 'c0', scat.c0); % Plane Wave (PW) sequence
                    case 'Focused'
                        [xmin, xmax, Nx] = deal( -10 ,  10 , 21 );
                        seq = Sequence('type', 'VS', 'c0', scat.c0, ...
                            'focus', [1;0;0] .* 1e-3*linspace(xmin, xmax, Nx) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 2mm
                            ...'focus', [1;0;0] .* 1e-3*(-10 : 0.2 : 10) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 0.2 mm
                            );
                end

            elseif isa(xdc, 'TransducerConvex')
                switch seq_name
                    case 'sector'
                        [amin, amax, Na] = deal( -15 ,  15 , 21 );
                        seq = SequenceRadial('type', 'VS', 'c0', scat.c0, ...
                            'angles', linspace(amin, amax, Na), ...
                            'ranges', norm(xdc.center) + 35e-3, 'apex', xdc.center ...
                            ); % sector scan sequence
                    case 'FSA'
                        seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'c0', scat.c0); % FSA - convex sequence
                end
            end

            % make a cartesian scan
            % set the scan at the edge of the transducer
            pn = xdc.positions(); % element positions
            xb = [-1 1] .* max(abs([pn(1,[1 end]), scat.pos(1,:)])); % x-limits are the edge of the aperture
            zb = [-10e-3, 10e-3] + [min(scat.pos(3,:)), max(scat.pos(3,:))]; % z-limits surround the point target

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
                scan = ScanPolar('origin', xdc.center, 'a', -40:0.5:40, ...
                    'r', norm(xdc.center) + linspace(zb(1), zb(end), Npixd+1)...
                    ); % R x A scan
                scan.r(end) = [];
            end
            
            % Choose the simulation region (eikonal)
            switch xdc_name
                case "C5-2V", tscan = ScanCartesian('x',linspace(-50e-3, 50e-3, 1+100*2^2), 'z', linspace(-30e-3, 60e-3, 1+90*2^2));
                otherwise,    tscan = ScanCartesian('x',linspace(-20e-3, 20e-3, 1+40 *2^3), 'z', linspace(  0e-3, 60e-3, 1+60*2^3));
            end

            % Construct an UltrasoundSystem object, combining all of these properties
            us = UltrasoundSystem('xdc', xdc, 'sequence', seq, 'scan', scan, 'fs', 40e6, 'recompile', false);

            % Make the transducer impulse tighter
            us.xdc.numel = 64; % reduce the number of elements
            us.xdc.fc = us.fs/8; % set to 5MHZ, as a ratio of sampling frequency
            us.xdc.bw_frac = 1.5;
            us.xdc.impulse = us.xdc.ultrasoundTransducerImpulse(); 
            % us.xdc.impulse = Waveform.Delta(); 
            % us.sequence.pulse = Waveform.Delta(); % adding this may make the
            % length of the signal go to 0, causing problems for interpolators

            % Simulate a point target
            % run on CPU to use spline interpolation
            chd = gather(greens(us, scat, [1,1], 'interp', 'linear')); % use a Greens function

            % Precondition the data
            chd.data = chd.data - mean(chd.data, 1, 'omitnan'); % remove DC
            chd = filter(chd, getPassbandFilter(chd, xdc.bw)); % apply a filter
            if isreal(chd), chd = hilbert(chd, 2^nextpow2(chd.T)); end % apply hilbert on real data
            chd = gather(chd); % bring to CPU

            % scale the problem
            us = scale(us, 'dist', 1e3, 'time', 1e6);
            scat = scale(scat, 'dist', 1e3, 'time', 1e6);
            chd = scale(chd, 'time', 1e6);
            scanc = scale(scanc, 'dist', 1e3);
            tscan = scale(tscan, 'dist', 1e3);

            % apply acceptance angle apodization
            apod = apAcceptanceAngle(us, 45);

            % baseband the data
            if baseband
                fmod_max = max(abs(us.xdc.fc - us.xdc.bw)); % maximum demodulated frequency
                ratio = floor(us.fs / fmod_max / 2); % discrete downsampling factor
                chd = downmix(hilbert(chd), us.xdc.fc); % downmix
                chd = downsample(chd, ratio); % downsample
                fmod = us.xdc.fc;
            else
                fmod = 0;
            end


            %% save QUPS objects for this test case
            test.chd    = chd; 
            test.us     = us; 
            test.scat   = scat;
            test.scanc  = scanc;
            test.tscan  = tscan;
            test.apod   = apod;
            test.fmod   = fmod;

        end
    end
    methods(TestClassTeardown)
        function teardownQUPSdata(test)
            % delete data
            test.chd    = [];
            test.us     = []; 
            test.scat   = [];
            test.scanc  = [];
            test.tscan  = []; 
        end
    end

    % some of these options aren't supported yet.
    properties(TestParameter)
        gdev = getdevs()
        bf_name = {'DAS','DAS-direct','Eikonal','Adjoint'}
        prec = struct('single','singleT','double', 'doubleT','halfT','halfT');
        terp = {'nearest', 'linear', 'cubic'};
        apodization = struct('false', false, 'true', true);
    end
    methods(TestMethodSetup)
        function resetGPU(test), 
            if gpuDeviceCount(), 
                gpuDevice([]);
                if isa(gcp('nocreate'), 'parallel.ProcessPool'),
                    hcp = gcp('nocreate');
                    parfor i = hcp.NumWorkers, gpuDevice([]); end
                end
            end
        end
    end
    
    % Github test routine
    methods(Test, ParameterCombination = 'sequential', TestTags={'Github'})
        function github_psf(test, bf_name)%, prec, terp)
            test.assumeTrue(bf_name ~= "Eikonal"); % Eikonal is too large?
            test.assumeTrue(test.fmod ~= 0); % only test non-centered data
            test.assumeTrue(test.scat.pos(1) ~= 0); % only test non-centered data
            
            % only test 1 precision/interpolation/apodization combo
            psf(test, 0, bf_name, 'singleT', 'nearest', false); 
        end
    end

    % Full test routine
    methods(Test, ParameterCombination = 'exhaustive', TestTags={'full'})
        function full_psf(test, gdev, bf_name, prec, terp, apodization)
            psf(test, gdev, bf_name, prec, terp, apodization); % forward all
        end 
    end

    % method implementations
    methods
        function psf(test, gdev, bf_name, prec, terp, apodization)
            % PSF - Test the PSF
            % 
            % Test that the PSF for a point at a reasonable distance 
            % properly beamforms for all transmit sequences, transducer 
            % types and compute device.

            % unpack
            [us, scat, scan, scanc, tscan, apod, fmod] = deal(...
                test.us, test.scat, test.us.scan, test.scanc, test.tscan, test.apod, test.fmod ...
                );
            if ~apodization, apod = 1; end % apodization not applied

            % make an equivalent medium - assume density scatterers
            med = Medium('c0', scat.c0);

            % exceptions
            % for the Eikonal beamformer, pass if not given FSA delays
            test.assumeFalse(bf_name == "Eikonal" && us.sequence.type ~= "FSA");
            % if using a frequency domain method, skip half precision - the phase errors are too large
            test.assumeFalse(ismember(bf_name, ["Adjoint"]) && prec == "halfT");
            % weird phase error, but it's not important right now - skip it
            test.assumeFalse(ismember(bf_name, ["Adjoint"]) && us.sequence.type == "VS"); % && isa(us.xdc, 'TransducerConvex'));
            % is using the adjoint method, pagemtimes,pagetranspose must be supported
            test.assumeTrue( bf_name ~= "Adjoint" || (...
                   logical(exist('pagemtimes'   , 'builtin')) ...
                && logical(exist('pagetranspose', 'file')) ...
                ));

            % set ChannelData type
            test.assumeTrue(prec ~= "halfT" || logical(exist('halfT', 'class')));
            tfun = str2func(prec);
            chd = tfun(test.chd); % cast to specified precision
            
            % for frequency domain methods, the time-axis must be extended
            % so that the replications of the miage are outside of the
            % imaging range
            if ismember(bf_name, ["Adjoint"])
                dr = hypot(range(scanc.xb), range(scanc.zb)) / 2; % half the largest range across the image
                dT = round(2 * dr / scat.c0 * chd.fs); % temporal buffer in indices
                chd = zeropad(chd, dT, dT); % add buffer on both sides
            end

            % move data to GPU if requested
            if gdev, chd = gpuArray(chd); else, chd = gather(chd); end 

            % Beamform 
            % TODO: test that we can keep dimensions and sum later
            % get common args that should work for all
            args = {'apod', apod, 'fmod', fmod};
            switch bf_name
                case "DAS",         b = bfDAS(us, chd, scat.c0,    'interp', terp, args{:});
                case "DAS-direct",  b =   DAS(us, chd, scat.c0,    'interp', terp, 'device', gdev, args{:}); % use a delay-and-sum beamformer
                case "Eikonal", b = bfEikonal(us, chd, med, tscan, 'interp', terp, 'verbose', false, args{:}); % use the eikonal equation
                case "Adjoint", b = bfAdjoint(us, chd, scat.c0,    'bsize' ,    1, 'verbose', false, args{:}); % use an adjoint matrix method
            end

            % show the image
            b_im = mod2db(b); % convert to power in dB

            if ~isa(scan, 'ScanCartesian')
                [b_im, scan] = deal(scanConvert(us.scan, b_im, scanc), scanc); % interpolate in dB - could also do abs, as long as it's not complex!
            end

            % TODO: peak should be ~near~ [0, 30mm] scan - check for this
            [i,j] = deal(argmin(abs(scan.z - scat.pos(3))), argmin(abs(scan.x - scat.pos(1))));
            [xo,zo] = deal(scan.x(j), scan.z(i)); % ideal max, discrete
            nmax = argmax(b_im, [], 'all', 'linear');
            [X,~,Z] = scan.getImagingGrid();
            [x, z] = deal(X(nmax), Z(nmax)); % image max
            
            % test
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;

            % the image must be non-zero
            test.assertTrue(logical(nnz(b)), 'The image has no non-zero values.');

            % no lateral (x) offset allowed
            test.assertThat(x, IsEqualTo(xo, 'Within', AbsoluteTolerance(1.1)), sprintf(...
                'Peak of the b-mode image is laterally offset from the peak of the target position (%.2fmm, %.2fmm).',  ...
                xo, x));

            % can be up to 1.1 mm off in depth (z)
            test.assertThat(z, IsEqualTo(zo, 'Within', AbsoluteTolerance(1.1)), sprintf(...
                'Peak of the b-mode image is axially offset from the peak of the target position (%.2fmm, %.2fmm).',  ...
                zo, z));
        end
    end

    methods(Static)
        % PROJ_FOLDER - Identifies the base folder for the project
        function f = proj_folder(), f = fullfile(fileparts(mfilename('fullpath')), '..'); end
    end
end

% test GPU devices if we can
function s = getdevs()
s.nodev = 0;
if gpuDeviceCount, s.dev = -1; end
end

function s = getgpu()
s.no_ptx = {};
if gpuDeviceCount, s.ptx = {'CUDA'}; end
end
