%% QUPS - Quick Ultrasound Processing & Simulation
%% What?
% QUPS is designed to be an accessible, verastile, shareable, lightweight codebase 
% designed to make developing and running ultrasound algorithms quick and easy! 
% It offers a standardization of data formatting that serves most common pulse-echo 
% ultrasound systems and supports homogeneous sound speed simulators such as FieldII 
% and MUST as well as heterogeneous simulation programs such as k-Wave.
%% Why?
% This package seeks to lower the barrier to entry for doing research in ultrasound 
% by offering a lightweight, easy to use package. Most of the underlying implementations 
% are here so that users can get to their work more efficiently! 
% 
% Currently you can:
%% 
% * Create linear and curvilinear (convex) transducers
% * Define full-synthetic-aperture, focused, diverging, or plane-wave transmits
% * Simulate point targets (via MUST or FieldII) or distributed media (via k-Wave)
% * Filter, demodulate, and downsample channel data
% * Define arbitrary apodization schemes over transmits, receives, and pixels, 
% jointly
% * Beamform using traditional delay-and-sum, the eikonal equation, or an adjoint 
% matrix method to create a b-mode image
%% 
% Please submit issues, feature requests or documentation requests via <https://github.com/thorstone25/qups/issues 
% github>!
% 
% 
% 
% 
% Setup the workspace
% (Run once)

 
%#ok<*UNRCH> ignore unreachable code due to constant values
%#ok<*BDLGI> ignore casting numbers to logical values
dev = -logical(gpuDeviceCount); % select 0 for cpu, -1 for gpu if you have one
% setup;  % add all the necessary paths
% setup parallel; % start a parpool for faster CPU processing
setup CUDA cache;  % setup CUDA paths & recompile and cache local mex/CUDA binaries
%% Create a simple simulation
% Choose a Target

target_depth = 1e-3 * 30;
switch "single"
    case 'single' , scat = Scatterers('pos', [0;0;target_depth], 'c0', 1500); % simple point target
    case 'array'  , 
        ps_x =  [1e-3;0;0] * (-10:5:10);
        ps = [ps_x + [0;0;30e-3], ps_x + [0;0;20e-3], ps_x + [0;0;40e-3]];
        scat = Scatterers('pos', ps, 'c0', 1500); % targets every 5mm
    case 'diffuse', N = 2500; % number of random scatterers
                    scat = Scatterers('pos', 1e-3*[40;10;60].*([1;0;1].*rand(3,N)-[0.5;0;0]), 'amp', rand(1,N), 'c0', 1500); % diffuse scattering
end
% Choose a transducer

switch "L11-5V"
    case 'L11-5V', xdc = TransducerArray.L11_5V(); % linear array
    case 'L12-3V', xdc = TransducerArray.L12_3V(); % another linear array
    case 'C5-2V' , xdc = TransducerConvex.C5_2V(); % convex array
end
% Choose the simulation region
% Make sure it covers the target and the transducer!

switch "small"
    case 'small', tscan = ScanCartesian('x',linspace(-20e-3, 20e-3, 40*2^4), 'z', linspace(  -4e-3, 36e-3, 40*2^4));
    case 'large', tscan = ScanCartesian('x',linspace(-50e-3, 50e-3, 100*2^4), 'z', linspace(-30e-3, 70e-3, 100*2^4));
end

% point per wavelength - aim for >2 for a simulation
ppw = scat.c0/xdc.fc/min([tscan.dx, tscan.dz]); 

% Choose a transmit sequence
% For linear transducers only!

if isa(xdc, 'TransducerArray') 
    switch "FSA"
        case 'FSA', seq = Sequence('type', 'FSA', 'c0', scat.c0, 'numPulse', xdc.numel); % set a Full Synthetic Aperture (FSA) sequence
        case 'Plane-wave', 
            [amin, amax, Na] = deal( -25 ,  25 , 11 );
            seq = SequenceRadial('type', 'PW', ...
                'ranges', 1, 'angles',  linspace(amin, amax, Na), 'c0', scat.c0); % Plane Wave (PW) sequence
        case 'Focused', 
            [xmin, xmax, Nx] = deal( -10 ,  10 , 11 );
            seq = Sequence('type', 'VS', 'c0', scat.c0, ...
        'focus', [1;0;0] .* 1e-3*linspace(xmin, xmax, Nx) + [0;0;50e-3] ... % translating aperture: depth of 30mm, lateral stride of 2mm
        ...'focus', [1;0;0] .* 1e-3*(-10 : 0.2 : 10) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 0.2 mm
        ); 
    end

    
%% 
% For convex transducers only!

elseif isa(xdc, 'TransducerConvex') 
    switch "sector"
        case 'sector'
            [amin, amax, Na] = deal( -40 ,  40 , 41 );
            seq = SequenceRadial('type', 'VS', 'c0', scat.c0, ...
                'angles', linspace(amin, amax, Na), ...
                'ranges', norm(xdc.center) + 60e-3, 'apex', xdc.center ...
                ); % sector scan sequence
        case 'FSA'
            seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'c0', scat.c0); % FSA - convex sequence
    end
end
% Choose an imaging region

% For linear transducers only!
if isa(xdc, 'TransducerArray') 
    % set the scan at the edge of the transducer
    pn = xdc.positions(); % element positions
    xb = pn(1,[1,end]); % x-limits are the edge of the aperture
    zb = [-10e-3, 10e-3] + [min(scat.pos(3,:)), max(scat.pos(3,:))]; % z-limits surround the point target
    switch "low"
        case "high"
            scan = ScanCartesian(...
                'x', linspace(xb(1), xb(end), 2^9), ...
                'z', linspace(zb(1), zb(end), 2^9) ...
                ); % X x Z scan
        case "low"
            scan = ScanCartesian(...
                'x', linspace(xb(1), xb(end), 2^6), ...
                'z', linspace(zb(1), zb(end), 2^6) ...
                ); % X x Z scan
    end
% For convex transducers only!
elseif isa(xdc, 'TransducerConvex') 

    % use with a SequenceRadial!
    scan = ScanPolar('origin', xdc.center, 'a', -40:0.5:40, ...
        'r', norm(xdc.center) + linspace(0, 2*max(vecnorm(scat.pos,2,1)), 2^9)...
        ); % R x A scan

end
%% 
% 

% create a distributed medium based on the point scatterers 
% (this logic will later be implemented in a class)
s_rad = max([tscan.dx, tscan.dz]); %, 260e-6); % scatterer radius
nextdim = @(p) ndims(p) + 1;
rho0 = 1000; % ambient density
if scat.numScat < 500
    ifun = @(p) any(vecnorm(p - swapdim(scat.pos,2,nextdim(p)),2,1) < s_rad, nextdim(p));
    med = Medium('c0', scat.c0, 'rho0', rho0, 'pertreg', {{ifun, [scat.c0, rho0*2]}});
else
    med = Medium('c0', scat.c0, 'rho0', rho0);
end
%% 
% 

% Show the transducer's impulse response
% figure; plot(xdc.impulse, '.-'); title('Element Impulse Response'); 

% Show the transmit signal - a single point means it's a delta function
% figure; plot(seq.pulse, '.-'); title('Transmit signal');
%  Plot configuration of the simulation

figure; hold on; title('Geometry');
switch "sound-speed"
    case 'sound-speed', imagesc(med, tscan, 'props', 'c'  ); colorbar; % show the background medium for the simulation/imaging region
    case 'density',     imagesc(med, tscan, 'props', 'rho'); colorbar; % show the background medium for the simulation/imaging region
end
plot(xdc, 'r+', 'DisplayName', 'Elements'); % elements
hps = plot(scan, 'w.', 'MarkerSize', 0.5, 'DisplayName', 'Image'); % the imaging points
switch seq.type % show the transmit sequence
    case 'PW', plot(seq, 3e-2, 'k.', 'DisplayName', 'Tx Sequence'); % scale the vectors for the plot
    otherwise, plot(seq, 'k.', 'DisplayName', 'Tx Sequence'); % plot focal points, if they exist
end
if scat.numScat < 500
    plot(scat, 'k', 'LineStyle', 'none', 'Marker', 'diamond', 'MarkerSize', 5, 'DisplayName', 'Scatterers'); % point scatterers
else
    disp('Info: Too many scatterers to display.');
end
hl = legend(gca, 'Location','bestoutside'); 
set(gca, 'YDir', 'reverse'); % set transducer at the top of the image
%% Simulate the Point Scatterer(s)

 
%% 
% 

% Construct an UltrasoundSystem object, combining all of these properties
us = UltrasoundSystem('xdc', xdc, 'sequence', seq, 'scan', scan, 'fs', single(40e6), 'recompile', false);

% Simulate a point target
tic;
switch "Greens"
    case 'Greens' , chd0 = greens(us, scat); % use a Greens function with a GPU if available!su-vpn.stanford.edu
    case 'FieldII', chd0 = calc_scat_all(us, scat); % use FieldII, 
    case 'FieldII-multi', chd0 = calc_scat_multi(us, scat); % use FieldII, 
    case 'SIMUS'  , us.fs = 4 * us.fc; % to address a bug in MUST where fs must be a ~factor~ of 4 * us.fc
                    chd0 = simus(us, scat, 'periods', 1, 'dims', 3); % use MUST: note that we have to use a tone burst or LFM chirp, not seq.pulse
    case 'kWave', chd0 = kspaceFirstOrder(us, med, tscan, 'CFL_max', 0.5, 'PML', [64 128], 'parenv', 0, 'PlotSim', true); % run locally, and use an FFT friendly PML size
end
toc; 
chd0
%%

% display the channel data across the transmits
chd = mod2db(chd0); % == 20*log10(abs(x)) -> the power of the signl in decibels
figure; h = imagesc(chd, 1); colormap jet; colorbar; caxis(gather([-80 0] + (max(chd.data(chd.data < inf)))))
xlabel('Channel'); ylabel('Time (s)'); ylim([min(chd.time(:)), max(chd.time(:))]);
for m = 1:size(chd.data,3), if isvalid(h), h.CData(:) = chd.data(:,:,m); h.YData(:) = chd.time(:,:,min(m,size(chd.time,3))); drawnow limitrate; title(h.Parent, "Tx " + m); pause(1/10); end, end
%% Create a B-mode Image

 
%% 
% 

% Precondition the data
chd = chd0;
chd = singleT(chd); % use less data

% apply a passband filter to retain only the badwidth of the transducer
D = chd.getPassbandFilter(xdc.bw, 25); % get a passband filter for the transducer bandwidth
chd = filter(chd, D); % apply passband filter for transducer bandwidth

if isreal(chd.data), chd = hilbert(chd, 2^nextpow2(chd.T)); end % apply hilbert on real data

% optionally demodulate and downsample the data
demod = true;
demod_thresh_db = 80; % db threshold for demodulation
if demod
    fpow = max(fft(chd).data, [], setdiff(1:ndims(chd.data), chd.tdim)); % max power per frequency
    if_max = gather(find(mod2db(fpow) > max(mod2db(fpow)) - demod_thresh_db, 1, 'last')); % dB threshold
    fmod = (if_max - 1) * chd.fs / chd.T; % maximum frequency past the threshold
    chd = downmix(chd, fmod); % downmix - this reduces the central frequency of the data    
    chd = downsample(chd, floor(chd.fs / fmod)); % now we can downsample without losing information

end % demodulate and downsample (by any whole number)
if dev, chd = gpuArray(chd); end % move data to GPU

% Run a simple DAS algorithm
switch class(xdc)
    case 'TransducerArray' , scale = xdc.pitch; % Definitions in elements
    case 'TransducerConvex', scale = xdc.angular_pitch; % Definitions in degrees
end
switch seq.type
    case "VS", 

%% 
% Choose an apodization method for Virtual Source (VS) transmit sequences

        switch "none"
            case 'multiline', apod = multilineApodization(us.scan, us.sequence);
            case 'scanline', apod = scanlineApodization(us.scan, us.sequence);
            case 'translating', apod = translatingApertureApodization(us.scan, us.sequence, us.rx, 32*scale);
            case 'aperture-growth', apod = apertureGrowthApodization(us.scan, us.sequence, us.rx, 1.8);
            case 'accept', apod = acceptanceAngleApodization(us.scan, us.sequence, us.rx, 55); 
            case 'none', apod = 1;
        end

%% 
% Choose an apodization method for Full Synthetic Aperture (FSA) transmit sequences

    case "FSA"
        us.sequence.focus = us.tx.positions(); % set the sequence foci to be the location of the transmitters for these profiles
        switch "none"
            case 'translating', apod = translatingApertureApodization(us.scan, us.sequence, us.rx, 32*scale);
            case 'aperture-growth', apod = apertureGrowthApodization(us.scan, us.sequence, us.rx, 1.8);
            case 'accept', apod = acceptanceAngleApodization(us.scan, us.sequence, us.rx, 55); 
            case 'none', apod = 1;
        end

%% 
% Choose an apodization method for Plane Wave (PW) transmit sequences

    case "PW"
        switch "none"
            case 'aperture-growth', apod = apertureGrowthApodization(us.scan, us.sequence, us.rx, 2);
            case 'accept', apod = acceptanceAngleApodization(us.scan, us.sequence, us.rx, 20); 
            case 'none', apod = 1;
        end
        
        
    otherwise, apod = 1; % apodization profiles for plane-wave not implemented :(
end
%% 
% Choose a beamforming method

bf_args = {'apod', apod, 'fmod', fmod}; % arguments for all beamformers
switch "DAS"
    case "DAS"
        b = bfDAS(us, chd, scat.c0, bf_args{:}); % use a vanilla delay-and-sum beamformer
    case "Adjoint"
        b = bfAdjoint(us, chd, scat.c0, 'fthresh', -20, bf_args{:}); % use an adjoint matrix method
    case "Eikonal"
        b = bfEikonal(us, chd, med, tscan, bf_args{:}); % use the eikonal equation
    case "DAS-direct"
        b = DAS(us, chd, scat.c0, bf_args{:}); % use a specialized delay-and-sum beamformer
end

% show the image
b_im = mod2db(b); % convert to power in dB
figure; imagesc(scan, b_im(:,:,1), [-60, 0] + max(b_im(:))); % display with 80dB dynamic range
colormap gray; colorbar;
%% Scan Convert for a Sector Scan

if ~isa(scan, 'ScanCartesian')
    scanc = scanCartesian(scan); % mind the caps! this gives us a ScanCartesian
    [scanc.nx, scanc.nz] = deal(2^9); % don't change boundaries, but increase resolution
    bc_im = (scanConvert(scan, mod2db(b), scanc)); % interpolate in dB - could also do abs, as long as it's not complex!

    imagesc(scanc, bc_im, [-60, 0] + max(bc_im(bc_im < Inf)));
    colormap gray; colorbar;
end