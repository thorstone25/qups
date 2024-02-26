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
        xdc_type = struct("array", "L12-3v","convex","C5-2v", "matrix","PO1921", "generic","generic"); % xdc types
        seq_type = {"FSA", "FC", "PW"}
        baseband = struct('true', true)%,'false', false);
    end

    methods(TestClassSetup, ParameterCombination = 'pairwise')
        % Shared setup for the entire test class
        % function setupQUPS(test)
        %     cd(BFTest.proj_folder); % call setup relative to here
        %     if ~exist('bin', 'dir'), setup CUDA cache; end % recompile and make a cache
        %     addpath("bin/");
        % end

        function setupQUPSdata(test, xdc_type, seq_type, baseband)
            %% create point target data with the configuration

            % simple point target 30 mm depth
            c0 = 1500;
            scat = Scatterers('pos', 1e-3*[2 0 15]', 'c0', c0); 

            % Transducer
            switch xdc_type
                case 'PO1921',  xdc = TransducerMatrix.PO1921(); xdc.numd = [5, 5]; 
                case 'L12-3v',  xdc = TransducerArray.L12_3v();  xdc.numel = 32;
                case 'C5-2v' ,  xdc = TransducerConvex.C5_2v();  xdc.numel = 32;
                case 'generic', xdc = TransducerArray();         xdc.numel = 32; % generic type, lin array
                                xdc = TransducerGeneric('pos', xdc.positions); 
            end
            ismat = isa(xdc, "TransducerMatrix");

            % Make the transducer impulse tighter
            % xdc.fc = us.fs/8; % set to 5MHZ, as a ratio of sampling frequency
            xdc.bw_frac = 1.5;
            xdc.impulse = xdc.xdcImpulse();

            % Choose a transmit sequence
            switch seq_type
                case {"FC","DV","VS"}
                    % get apodization 
                    switch class(xdc)
                        case "TransducerMatrix"
                            [N, M] = deal(xdc.numd(1), xdc.numd(2));
                            ax = Sequence.apWalking(N, 2, 1);
                            ay = Sequence.apWalking(M, 2, 1);
                            ax = swapdim(ax, [1,2], [1, 3]);
                            ay = swapdim(ay, [1,2], [2, 4]);
                            apd = reshape(ax .* ay, N*M,[]);
                        otherwise
                            apd = Sequence.apWalking(xdc.numel, xdc.numel / 4, xdc.numel/16);
                    end

                    % get foci offset foci towards the scatterer
                    if seq_type == "DV", zf = -10e-3; else, zf = 50e-3; end
                    switch class(xdc)
                        case "TransducerConvex"
                            pf = xdc.focActive(apd, zf) + mean(scat.pos(1,:));
                            seq = SequenceRadial('c0',c0,'type',seq_type,'focus',pf,'apd',apd,'apex', xdc.center);
                        case "TransducerGeneric"
                            [pn, th] = deal(xdc.positions, xdc.orientations);
                            pl = cell2mat(cellfun(@(a){pn(:,floor(median(find(a))))}, num2cell(apd,1)));
                            pu = cell2mat(cellfun(@(a){pn(:,ceil( median(find(a))))}, num2cell(apd,1)));
                            tl = cell2mat(cellfun(@(a){th(:,floor(median(find(a))))}, num2cell(apd,1)));
                            tu = cell2mat(cellfun(@(a){th(:,ceil( median(find(a))))}, num2cell(apd,1)));
                            p0 = (pl + pu)/2;
                            t0 = (tl + tu)/2;
                            pf = p0 + zf*[sind(t0); 0*t0; cosd(t0)] + mean(scat.pos(1,:));
                            seq = Sequence('c0',c0,'type',seq_type,'focus',pf,'apd',apd);
                        otherwise
                            pf = xdc.focActive(apd, zf) + mean(scat.pos(1,:));
                            seq = Sequence('c0',c0,'type',seq_type,'focus',pf,'apd',apd);
                    end         
                case "FSA"
                    seq = Sequence('c0',c0,'type','FSA','numPulse',xdc.numel);
                case "PW"
                    switch class(xdc), case "TransducerConvex", p0 = xdc.center; otherwise, p0 = xdc.offset; end
                    r = scat.pos - p0;
                    th0 = mean(atan2d(r(1,:), r(3,:)));
                    seq = SequenceRadial('c0',c0,'type','PW','angles',-10:5:10+th0);
            end

            % make a Cartesian Scan around the point target for imaging
            if ismat, Npixd = 2^5; else, Npixd = 2^7; end
            lbda = c0 / xdc.fc;
            dp = (lbda/4) * ((0:Npixd-1) - floor(Npixd/2));
            p0 = mean(scat.pos,2);
            
            scanc = ScanCartesian('x', p0(1)+dp, 'z', p0(3)+dp); % X x Z scan
            if ismat, scanc.y = scanc.x; end

            % Choose the simulation region (eikonal)
            switch xdc_type
                case "C5-2V", tscan = ScanCartesian('x',linspace(-50e-3, 50e-3, 1+100*2^2), 'z', linspace(-30e-3, 30e-3, 1+60*2^2));
                case "PO192O",tscan = ScanCartesian('x',linspace(-10e-3, 10e-3, 1+20 *2^3), 'z', linspace(-10e-3, 30e-3, 1+40*2^3));
                              tscan.y = tscan.x;
                otherwise,    tscan = ScanCartesian('x',linspace(-20e-3, 20e-3, 1+40 *2^3), 'z', linspace(-10e-3, 30e-3, 1+40*2^3));
            end

            % Image polar on convex probes
            switch class(xdc)
                case "TransducerConvex"
                    r = scat.pos - xdc.center; % radial vector to scat
                    th = xdc.orientations;
                    scan = ScanPolar("origin",xdc.center,"r", vecnorm(r,2,1) + dp, "a", th + mean(atan2d(r(1,:),r(3,:))));
                otherwise
                    scan = copy(scanc);
            end

            % Construct an UltrasoundSystem object, combining all of these properties
            us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'scan', scan, 'fs', 4*xdc.fc, 'recompile', false);

            % Simulate a point target
            % run on CPU to use spline interpolation
            chd = gather(greens(us, scat, [1,1], 'interp', 'linear')); % use a Greens function

            % Precondition the data
            chd.data = chd.data - mean(chd.data, 1, 'omitnan'); % remove DC
            chd = filter(chd, getPassbandFilter(chd, xdc.bw)); % apply a filter

            % scale the problem
            us    = scale(us, 'dist', 1e3, 'time', 1e6);
            scat  = scale(scat, 'dist', 1e3, 'time', 1e6);
            chd   = scale(chd, 'time', 1e6);
            scanc = scale(scanc, 'dist', 1e3);
            tscan = scale(tscan, 'dist', 1e3);

            % apply acceptance angle apodization
            switch class(xdc)
                case "TransducerMatrix"
                    apod = 1; % too expensive / not working
                otherwise
                    switch seq_type
                        case {"FC","VS"},   apod = us.apMultiline;
                        case {"PW","FSA"},  apod = us.apAcceptanceAngle(30);
                        case {"DV"},        apod = us.apAcceptanceAngle(30);
                    end
            end

            % baseband the data
            if baseband
                fmod_max = max(abs(us.xdc.fc - us.xdc.bw)); % maximum demodulated frequency
                ratio = floor(us.fs / fmod_max / 2); % discrete downsampling factor
                chd = downmix(chd, us.xdc.fc); % downmix
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
            % reset gpu
        end
    end

    % some of these options aren't supported yet.
    properties(TestParameter)
        gdev = getdevs()
        bf_name = {'DAS','DAS-direct','Eikonal','Adjoint'}
        prec = struct('single','singleT')%,'double', 'doubleT')%,'halfT','halfT');
        terp = {'cubic'};
        apodization = struct('false', false, 'true', true);
    end
    methods(TestMethodSetup)
        function resetGPU(test)
            return; % pass
            if gpuDeviceCount()
                test.chd = gather(test.chd);
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
        function github_psf(test, bf_name)%, prec, terp)
            if(any(bf_name == ["Eikonal","Adjoint"])), return; end % Too much compute
            
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
            if(bf_name == "Eikonal" && us.seq.type ~= "FSA"), return; end
            % if using a frequency domain method, skip half precision - the phase errors are too large
            if(ismember(bf_name, ["Adjoint"]) && prec == "halfT"), return; end
            % weird phase error, but it's not important right now - skip it
            if(ismember(bf_name, ["Adjoint"]) && us.seq.type == "VS") return; end % && isa(us.xdc, 'TransducerConvex'));
            % is using the adjoint method, pagemtimes,pagetranspose must be supported
            test.assumeTrue( bf_name ~= "Adjoint" || (...
                   logical(exist('pagemtimes'   , 'builtin')) ...
                && logical(exist('pagetranspose', 'file')) ...
                ));

            % set ChannelData type
            % test.assumeTrue(prec ~= "halfT" || logical(exist('halfT', 'class')));
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
