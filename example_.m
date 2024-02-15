%% QUPS - Quick Ultrasound Processing & Simulation
%% What?
% QUPS is designed to be an accessible, versatile, shareable, lightweight codebase 
% designed to make developing and running ultrasound algorithms quick and easy! 
% It offers a standardization of data formatting that serves most common pulse-echo 
% ultrasound systems and supports scatterer simulators such as FieldII and MUST 
% as well as finite difference simulation programs such as k-Wave.
%% Why?
% This package seeks to lower the barrier to entry for doing research in ultrasound 
% by offering a lightweight, easy to use package. Most of the underlying implementations 
% are here so that users can get to their work more efficiently! 
% 
% Currently you can:
%% 
% * Create linear, curvilinear (convex), and matrix transducers, as well as 
% custom transducers
% * Define full-synthetic-aperture (FSA), focused (FC), diverging (DV), plane-wave 
% (PW), or arbitrary pulse sequences
% * Simulate point targets natively or via MUST or FieldII
% * Simulate distributed media via k-Wave
% * Filter, downmix (demodulate), downsample (decimate), and resample channel 
% data
% * Define arbitrary apodization schemes across transmits, receives, and pixels, 
% jointly
% * Beamform using traditional delay-and-sum, eikonal delays (via the FMM toolbox), 
% an adjoint matrix method, or Stolt's migration
% * Focus FSA data to some pulse sequence, or retrospectively refocus data back 
% to FSA
% * Compute coherence based images with techniques including SLSC, DMAS, and 
% coherence factor
%% 
% Please submit issues, feature requests or documentation requests via <https://github.com/thorstone25/qups/issues 
% github>!
% 
% 
% 
% 
% Setup the workspace

 
%#ok<*UNRCH> ignore unreachable code due to constant values
%#ok<*BDLGI> ignore casting numbers to logical values
gpu = logical(gpuDeviceCount); % set to false to remain on the cpu
if ~exist('UltrasoundSystem.m','file')
    % this setup function adds the proper paths to use QUPS: it only needs
    % to be run once.

    % setup;  % add all the necessary paths
    % setup parallel; % start a parpool for faster CPU processing
    setup CUDA cache; % setup CUDA paths & recompile and cache local mex/CUDA binaries (ignore any warnings)
end

% if Matlab-OpenCL is installed, you can select an OpenCL device to
% accelerate some functions
if exist('oclDeviceCount', 'file') && oclDeviceCount()
    oclDevice(1); % select the first OpenCL device
end 
%% Create a simple simulation
% Create some Scatterers

switch "grid"
    case 'single'  
        target_depth = 1e-3 * 30;
        scat = Scatterers('pos', [0;0;target_depth], 'c0', 1500); % simple point target
    case 'grid'  
        [~, px, py, pz] = ndgrid(0, -10:5:10, 0, 20:10:40); % grid in each axes
        ps = [px; py; pz]; % form a grid of point targets
        scat = Scatterers('pos', 1e-3 * ps(:,:), 'c0', 1500); % targets every 5mm
    case 'diffuse', N = 1000; % number of random scatterers
        sz = 1e-3*[40;10;60]; % rectangular region size
        cen = [true;true;false]; % uniform distribution offset
        scat = Scatterers('pos', sz.*(rand([3,N]) - 0.5*cen), 'amp', rand([1,N]), 'c0', 1500); % diffuse scattering
end
% Choose a Transducer

switch "L11-5v"
    case 'L11-5v', xdc = TransducerArray.L11_5v();  % linear array
    case 'L12-3v', xdc = TransducerArray.L12_3v();  % another linear array
    case 'P4-2v',  xdc = TransducerArray.P4_2v();   % a phased array 
    case 'C5-2v' , xdc = TransducerConvex.C5_2v();  % convex array
    case 'PO192O', xdc = TransducerMatrix.PO192O(); % matrix array
end
xdc.impulse = xdc.ultrasoundTransducerImpulse(); % set the impulse response function
% Define the Simulation Region

if isa(xdc, 'TransducerMatrix')
    tscan = ScanCartesian('x', 1e-3*(-10 : 1/16 : 10), 'z', 1e-3*(-5 : 1/16 : 50)); 
    tscan.y = tscan.x; % 3D grid
else
    tscan = ScanCartesian('x', 1e-3*(-50 : 1/16 : 50), 'z', 1e-3*(-20 : 1/16 : 60)); % 2D grid
    tscan.y = 0; % 2D grid
end

% point per wavelength - aim for >2 for a simulation
ppw = scat.c0 / xdc.fc / min([tscan.dx, tscan.dy, tscan.dz], [], 'omitnan'); 

% Define the Transmit  Pulse Sequence

switch "FSA"
    case 'FSA'
        seq = Sequence('type', 'FSA', 'numPulse', xdc.numel); % set a Full Synthetic Aperture (FSA) sequence
    case 'Plane-wave'

%% 
% Plane waves:

        th = linspace( -15 ,  15 , 31 );
        seq = SequenceRadial('type', 'PW', 'angles', th); % Plane Wave (PW) sequence
    case 'Focused'
        
%% 
% Focused pulses:

        zf = 60 ; xf = linspace( -10 ,  10 , 11 ); % Focused (VS) Sequence  
        pf = 1e-3*([1,0,0]'.*xf + [0,0,zf]');
        seq = Sequence('type', 'FC', 'focus', pf);
    case 'Sector'
        
%% 
% Sector scans:

        r = 60 ; th = linspace( -40 ,  40 , 41 ); % Sector Scan
        if isa(xdc, 'TransducerConvex'), cen = xdc.center;
        else,                            cen = [0;0;0]; 
        end
        seq = SequenceRadial( ...
            'type', 'FC', ...
            'angles', th, ...
            'ranges', norm(cen) + 1e-3*r, ...
            'apex', cen ...
            ); % sector scan sequence
end

% set the sound speed to match the scatterers
seq.c0 = scat.c0;

% Define the imaging region

% set the spatial resolution
dr = deal(scat.c0 / xdc.fc / 4);

% For linear transducers only!
if isa(xdc, 'TransducerArray') 
    % set the scan boundaires at the edge of the transducer with a 5mm
    % buffer
    pb = xdc.bounds(); % element positions
    xb = pb(1,[1,end]) + [-5e-3, 5e-3]; % x-limits are the edge of the aperture
    zb = [-15e-3, 15e-3] + [min(scat.pos(3,:)), max(scat.pos(3,:))]; % z-limits surround the point target

    scan = ScanCartesian('x', xb, 'z', zb);
    [scan.dx, scan.dz] = deal(dr);

    % For convex transducers only!
elseif isa(xdc, 'TransducerConvex') 
    % ranges
    r = 0 : dr : max(vecnorm(scat.pos,2,1)) + 15e-3;

    % use with a SequenceRadial for best results
    scan = ScanPolar('origin', xdc.center, ... set the cartesian origin of the axes
        'a', -40:0.5:40, ... angles
        'r', norm(xdc.center) + r...
        ); % R x A scan

end

% create a distributed medium based on the point scatterers 
% (this logic will later be implemented in a class)

% ambient density
rho0 = 1000; 

% make a scatterer, if not diffuse
if scat.numScat < 500
    s_rad = max([tscan.dx,tscan.dy,tscan.dz],[],'omitnan'); % scatterer radius
    ifun = @(p) any(vecnorm(p - swapdim(scat.pos,2,5),2,1) < s_rad, 5); % finds all points within scatterer radius
    med = Medium('c0', scat.c0, 'rho0', rho0, 'pertreg', {{ifun, [scat.c0, rho0*2]'}});
else
    med = Medium('c0', scat.c0, 'rho0', rho0);
end
%% 
% 

% Show the transducer's impulse response
if false, figure; plot(xdc.impulse, '.-'); title('Element Impulse Response'); end

% Show the transmit signal - a single point means it's a delta function
if false, figure; plot(seq.pulse, '.-'); title('Transmit signal'); end
%  
% Plot configuration of the simulation

figure; hold on; title('Geometry');

% plot the medium
imagesc(med, tscan, 'props', "c"); colorbar; % show the background medium for the simulation/imaging region

% Construct an UltrasoundSystem object, combining all of these properties
fs = single(ceil(2.5*xdc.bw(end)/1e6)*1e6);
us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'scan', scan, 'fs', fs, 'recompile', false);

% plot the system
hold on;
hs = plot(us);

% plot the point scatterers
if scat.numScat < 500
    hs(end+1) = plot(scat, 'k', 'LineStyle', 'none', 'Marker', 'diamond', 'MarkerSize', 5, 'DisplayName', 'Scatterers'); % point scatterers
else
    disp('INFO: Too many scatterers to display.');
end
hl = legend(gca, 'Location','bestoutside');
%% Simulate the Point Scatterer(s)

 
%% 
% 


% Simulate a point target
tic;
switch "Greens"
    case 'Greens' , chd0 = greens(us, scat); % use a Greens function with a GPU if available!su-vpn.stanford.edu
    case 'FieldII', chd0 = calc_scat_all(us, scat); % use FieldII to simulate FSA, then focus in QUPS
    case 'FieldII-multi', 
                    chd0 = calc_scat_multi(us, scat); % use FieldII to simulate sequence directly
    case 'SIMUS',   us.fs = 4 * us.fc; % to address a bug in early versions of MUST where fs must be a ~factor~ of 4 * us.fc
                    chd0 = simus(us, scat, 'periods', 1, 'dims', 3); % use MUST: note that we have to use a tone burst or LFM chirp, not seq.pulse
    case 'kWave',   chd0 = kspaceFirstOrder(us, med, tscan, 'CFL_max', 0.5, 'PML', [64 128], 'parenv', 0, 'PlotSim', true); % run locally, and use an FFT friendly PML size
end
toc; 
chd0
% Show the Channel Data

 
% display the channel data across the transmits
figure; 
h = imagesc(chd0); 
colormap jet; colorbar; 
cmax = gather(mod2db(max(chd0.data(:)))); % maximum power
caxis([-80 0] + cmax); % plot up to 80 dB down
% animate(chd0.data, h, 'fs', 10, 'loop', false); % show all transmits
for i = 1:chd0.M, h.CData(:) = mod2db(chd0.data(:,:,i)); drawnow limitrate; pause(1/10); end % implement above manually for live editor
%% Create a B-mode Image

 
%% 
% 

% Precondition the data
chd = singleT(chd0); % convert to single precision to use less data
if gpu, chd = gpuArray(chd); end % move data onto a GPU if one is available

% optionally apply a passband filter to retain only the bandwidth of the transducer
D = chd.getPassbandFilter(xdc.bw, 25); % get a passband filter for the transducer bandwidth
chd = filter(chd, D); % apply passband filter for transducer bandwidth

% apply hilbert on real data to get the complex analog
if isreal(chd.data), chd = hilbert(chd, 2^nextpow2(chd.T)); end 

% optionally demodulate and downsample the data
demod = false;
if demod
    demod_thresh_db = 80; % db threshold for demodulation
    
    % find the largest frequency more the 80dB down from the peak
    f = chd.fftaxis();
    fpow = mod2db(max(fft(hilbert(chd)).data, [], 2:3)); % max power per frequency
    if_max = find(fpow > max(fpow) - demod_thresh_db, 1, 'last'); % dB threshold
    fmod = f(if_max); % maximum frequency past the threshold
    
    chd = downmix(chd, fmod); % downmix - this reduces the central frequency of the data    
    chd = downsample(chd, floor(chd.fs / fmod)); % now we can downsample without losing information
else
    fmod = 0; % no modulation
end

% Choose how to scale apodization, laterally or angularly
switch class(xdc)
    case 'TransducerArray' , scl = xdc.pitch;           % Definitions in elements
    case 'TransducerConvex', scl = xdc.angular_pitch;   % Definitions in degrees
    case 'TransducerMatrix', scl = min(abs(xdc.pitch)); % Definitions in elements
    otherwise,               scl = us.lambda;           % unknown - default to 1 wavelength
        warning("Unknown scaling for a "+class(xdc)+": defaulting to 1 wavelength.");
end

% Choose the apodization (beamforming weights)

%% 
% Choose apodization scheme(s) suitable for the transmit sequence.
% 
% Scanline / Multiline <- Use either of these with focused or diverging sequences 
% (FC / DV) to emulate a traditional "scan-line" beamformer
% 
% Translating Aperture <- Use this with focused or diverging sequences (FC / 
% DV) to restrict the number of active elements per transmit
% 
% Aperture Growth / Acceptance Angle <- Use either of these with any sequence 
% to improve SNR throughtou the image

apod = 1;
if false, apod = apod .* apMultiline(us); end
if false, apod = apod .* apScanline(us); end
if false, apod = apod .* apTranslatingAperture(us, 32*scl); end
if false, apod = apod .* apApertureGrowth(us, 2); end
if false, apod = apod .* apAcceptanceAngle(us, 50); end

%% 
% Choose a beamforming method

bf_args = {'apod', apod, 'fmod', fmod}; % arguments for all beamformers
bscan = us.scan; % default b-mode image scan
switch "DAS"
    case "DAS-direct"
        b = DAS(      us, chd, bf_args{:}); % use a specialized delay-and-sum beamformer
    case "DAS"
        b = bfDAS(    us, chd, bf_args{:}); % use a generic delay-and-sum beamformer
    case "DASLUT"
        [~, tau_rx, tau_tx] = bfDAS(us, chd, 'delay_only', true); % pre-compute delays
        b = bfDASLUT( us, chd, tau_rx, tau_tx, bf_args{:}); % use a look-up table with pre-computed delays
    case "Adjoint"
        b = bfAdjoint(us, chd, bf_args{:}, 'fthresh', -20); % use an adjoint matrix method, top 20dB frequencies only
    case "Eikonal"
        b = bfEikonal(us, chd, med, tscan, bf_args{:}); % use the eikonal equation
    case "Migration"
        % Stolt's f-k Migrartion (no apodization accepted)
        % NOTE: this function works best with small-angle (< 10 deg) plane waves
        [b, bscan] = bfMigration(us, chd, 'Nfft', [2*chd.T, 4*chd.N], 'fmod', fmod); 
end

% show the image
dbR = 40;
figure; imagesc(bscan, b); % display
colormap gray; colorbar;
caxis(max(caxis) + [-dbR 0]); % 60dB dynamic range
%% Scan Convert for a Sector Scan

if ~isa(bscan, 'ScanCartesian')
    scanc = ScanCartesian(bscan); % make a ScanCartesian from the ScanPolar
    [scanc.dx, scanc.dz] = deal(1/32e3); % (optional) increase resolution
    bc_im = (scanConvert(bscan, mod2db(b), scanc)); % interpolate in dB - could also do abs, as long as it's not complex!
    
    imagesc(scanc, bc_im); 
    caxis([-dbR, 0] + max(caxis)); % set 60dB dynamic range
    colormap gray; colorbar;
end