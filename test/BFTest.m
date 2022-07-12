classdef BFTest < matlab.unittest.TestCase
    % BFTest - Beamformer test class
    %
    % This class test that all beamformers properly beamform a point target
    
    %#ok<*PROPLC> object properties can have the same name as local variables
    properties
        chd % ChannelData
        us % UltrasoundSystemct
        targ % Target
        scanc % Image Scan (Cartesian)
        tscan % Sound-speed Scan
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
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
        % Setup for the entire test class
        function addcache(test)
            % ADDCACHE - add the cached bin folder
            % adds the cached bin folder - this prevents recompilation for
            % each and every test, which can be a majority of the
            % computational effort / time 
            try mkdir  (fullfile(BFTest.proj_folder,'bin')); end %#ok<TRYNC>
            try addpath(fullfile(BFTest.proj_folder,'bin')); end %#ok<TRYNC>
        end

        % Shared setup for the entire test class
        function setupQUPS(test, par, gpu)
            cd(BFTest.proj_folder); % call setup relative to here
            setup(par{:}, gpu{:}); % compile what we can
            if ~exist('bin', 'dir'), setup cache; end % recompile and make a cache
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

            % Choose a transmit sequence
            seq_name = xdc_seq_name(2);
            if isa(xdc, 'TransducerArray')
                switch seq_name
                    case 'FSA', seq = Sequence('type', 'FSA', 'c0', targ.c0, 'numPulse', xdc.numel); % set a Full Synthetic Aperture (FSA) sequence
                    case 'Plane-wave'
                        [amin, amax, Na] = deal( -25 ,  25 , 25 );
                        seq = SequenceRadial('type', 'PW', ...
                            'ranges', 1, 'angles',  linspace(amin, amax, Na), 'c0', targ.c0); % Plane Wave (PW) sequence
                    case 'Focused'
                        [xmin, xmax, Nx] = deal( -10 ,  10 , 21 );
                        seq = Sequence('type', 'VS', 'c0', targ.c0, ...
                            'focus', [1;0;0] .* 1e-3*linspace(xmin, xmax, Nx) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 2mm
                            ...'focus', [1;0;0] .* 1e-3*(-10 : 0.2 : 10) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 0.2 mm
                            );
                end

            elseif isa(xdc, 'TransducerConvex')
                switch seq_name
                    case 'sector'
                        [amin, amax, Na] = deal( -15 ,  15 , 21 );
                        seq = SequenceRadial('type', 'VS', 'c0', targ.c0, ...
                            'angles', linspace(amin, amax, Na), ...
                            'ranges', norm(xdc.center) + 35e-3, 'apex', xdc.center ...
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
                case "C5-2V", tscan = ScanCartesian('x',linspace(-50e-3, 50e-3, 1+100*2^2), 'z', linspace(-30e-3, 60e-3, 1+90*2^2));
                otherwise,    tscan = ScanCartesian('x',linspace(-20e-3, 20e-3, 1+40 *2^3), 'z', linspace(  0e-3, 60e-3, 1+60*2^3));
            end

            % Construct an UltrasoundSystem object, combining all of these properties
            us = UltrasoundSystem('xdc', xdc, 'sequence', seq, 'scan', scan, 'fs', 40e6);

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
            chd = gather(greens(us, targ, [1,1], 'interp', 'linear')); % use a Greens function

            % Precondition the data
            chd.data = chd.data - mean(chd.data, 1, 'omitnan'); % remove DC
            chd = filter(chd, getPassbandFilter(chd, xdc.bw)); % apply a filter
            if isreal(chd), chd = hilbert(chd, 2^nextpow2(chd.T)); end % apply hilbert on real data

            % save QUPS objects for this test case
            test.chd    = chd; 
            test.us     = us; 
            test.targ   = targ;
            test.scanc  = scanc;
            test.tscan  = tscan; 

        end
    end
    methods(TestClassTeardown)
        function teardownQUPSdata(test)
            % delete data
            test.chd    = [];
            test.us     = []; 
            test.targ   = [];
            test.scanc  = [];
            test.tscan  = []; 
        end
    end

    % some of these options aren't supported yet.
    properties(TestParameter)
        gdev = getdevs()
        bf_name = {'DAS','DAS-direct','Eikonal','Adjoint'}
        prec = struct('double', 'doubleT','single','singleT', 'halfT','halfT');
        terp = {'nearest', 'linear', 'cubic'};
    end
    methods(TestMethodSetup)
    end
    
    % Github test routine
    methods(Test, ParameterCombination = 'sequential', TestTags={'Github'})
        function github_psf(test, bf_name)%, prec, terp)
            switch bf_name, case {'Eikonal', 'Adjoint'}, return; end % Adjoint not supported, Eikonal too large?
            
            % only test 1 precision/interpolation
            psf(test, 0, bf_name, 'singleT', 'nearest'); 
        end
    end

    % Full test routine
    methods(Test, ParameterCombination = 'exhaustive', TestTags={'full'})
        function full_psf(test, gdev, bf_name, prec, terp)
            psf(test, gdev, bf_name, prec, terp); % forward all
        end 
    end

    % method implementations
    methods
        function psf(test, gdev, bf_name, prec, terp)
            % PSF - Test the PSF
            % 
            % Test that the PSF for a point at a reasonable distance 
            % properly beamforms for all transmit sequences, transducer 
            % types and compute device.

            % unpack
            [us, targ, scan, scanc, tscan] = deal(...
                test.us, test.targ, test.us.scan, test.scanc, test.tscan...
                );

            % for the Eikonal beamformer, pass if not given FSA delays
            if bf_name == "Eikonal" && us.sequence.type ~= "FSA", return; end

            % set ChannelData type
            tfun = str2func(prec);
            chd = tfun(test.chd); % cast to specified precision
            if gdev, chd = gpuArray(chd); else, chd = gather(chd); end % move data to GPU if requested

            % Beamform 
            switch bf_name
                case "DAS",         b = bfDAS(us, chd, targ.c0,     'interp', terp);
                case "DAS-direct",  b =   DAS(us, chd, targ.c0,     'interp', terp, 'device', gdev); % use a delay-and-sum beamformer
                case "Eikonal", b = bfEikonal(us, chd, targ, tscan, 'interp', terp); % use the eikonal equation
                case "Adjoint", b = bfAdjoint(us, chd, targ.c0                    ); % use an adjoint matrix method
            end

            % show the image
            b_im = mod2db(b); % convert to power in dB

            if ~isa(scan, 'ScanCartesian')
                [b_im, scan] = deal(scanConvert(scan, b_im, scanc), scanc); % interpolate in dB - could also do abs, as long as it's not complex!
            end

            % TODO: peak should be ~near~ [0, 30mm] scan - check for this
            [i,j] = deal(argmin(abs(scan.z - targ.pos(3))), argmin(abs(scan.x - targ.pos(1))));
            [xo,zo] = deal(scan.x(j), scan.z(i)); % ideal max, discrete
            nmax = argmax(b_im, [], 'all', 'linear');
            [X,~,Z] = scan.getImagingGrid();
            [x, z] = deal(X(nmax), Z(nmax)); % image max
            
            % test
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;

            % no lateral (x) offset allowed
            test.assertThat(x, IsEqualTo(xo, 'Within', AbsoluteTolerance(1.1e-3)), sprintf(...
                'Peak of the b-mode image is laterally offset from the peak of the target position (%.2fmm, %.2fmm).',  ...
                xo, x));

            % can be up to 1.1 mm off in depth (z)
            test.assertThat(z, IsEqualTo(zo, 'Within', AbsoluteTolerance(1.1e-3)), sprintf(...
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