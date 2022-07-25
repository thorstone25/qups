%% QUPS - Quick Ultrasound Processing & Simulation
%% What?
% QUPS is designed to be an accessible, verastile, shareable, lightweight codebase 
% designed to make developing and running ultrasound algorithms quick and easy! 
% It offers a standardization of data formatting that serves most common pulse-echo 
% ultrasound systems and supports simulation programs such as k-Wave or Fullwave.
%% Why?
% This package seeks to lower the barrier to entry for doing research in ultrasound 
% by offering a lightweight, easy to use package. Most of the underlying implementations 
% are here so that users can get to their work more efficiently! 
% 
% Currently you can:
%% 
% * Create linear and curvilinear (convex) transducers
% * Define full-synthetic-aperture, focused, diverging, or plane-wave transmits
% * Beamform using traditional delay-and-sum, the eikonal equation, or the adjoint 
% method to create a b-mode image
% * Simulate point targets (via MUST or FieldII) or distributed media (via k-Wave)
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
%% 
%% Create a simple simulation
% Choose a Target

target_depth = 1e-3 * 30;
switch "single"
    case 'single' , targ = Target('pos', [0;0;target_depth], 'c0', 1500); % simple point target
    case 'array'  , targ = Target('pos', [0;0;1] * 1e-3*(10:5:50), 'c0', 1500); % targets every 5mm
    case 'diffuse', N = 2500; % a random number of scatterers
                    targ = Target('pos', 1e-3*[40;10;60].*([1;0;1].*rand(3,N)-[0.5;0;0]), 'amp', rand(1,N), 'c0', 1500); % diffuse scattering
end
% make scatterers twice the density
targ.scat_mode = 'ratio';
targ.rho_scat = 2; 
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
ppw = targ.c0/xdc.fc/min([tscan.dx, tscan.dz]); 

% Choose a transmit sequence
% For linear transducers only!

if isa(xdc, 'TransducerArray') 
    switch "Focused"
        case 'FSA', seq = Sequence('type', 'FSA', 'c0', targ.c0, 'numPulse', xdc.numel); % set a Full Synthetic Aperture (FSA) sequence
        case 'Plane-wave', 
            [amin, amax, Na] = deal( -25 ,  25 , 11 );
            seq = SequenceRadial('type', 'PW', ...
                'ranges', 1, 'angles',  linspace(amin, amax, Na), 'c0', targ.c0); % Plane Wave (PW) sequence
        case 'Focused', 
            [xmin, xmax, Nx] = deal( -10 ,  10 , 11 );
            seq = Sequence('type', 'VS', 'c0', targ.c0, ...
        'focus', [1;0;0] .* 1e-3*linspace(xmin, xmax, Nx) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 2mm
        ...'focus', [1;0;0] .* 1e-3*(-10 : 0.2 : 10) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 0.2 mm
        ); 
    end

    
%% 
% For convex transducers only!

elseif isa(xdc, 'TransducerConvex') 
    switch "sector"
        case 'sector'
            [amin, amax, Na] = deal( -40 ,  40 , 41 );
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
% Choose an imaging region

% For linear transducers only!
if isa(xdc, 'TransducerArray') 
    % set the scan at the edge of the transducer
    pn = xdc.positions(); % element positions
    xb = pn(1,[1,end]); % x-limits are the edge of the aperture
    zb = [-10e-3, 10e-3] + [min(targ.pos(3,:)), max(targ.pos(3,:))]; % z-limits surround the point target
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
    scan = ScanPolar('origin', seq.apex, 'a', -40:0.5:40, ...
        'r', norm(seq.apex) + linspace(0, 2*max(vecnorm(targ.pos,2,1)), 2^9)...
        ); % R x A scan

end
%% 
% 

% create a distributed medium based on the point scatterers 
% (this logic will later be implemented in a class)
s_rad = max([tscan.dx, tscan.dz]); %, 260e-6); % scatterer radius
nextdim = @(p) ndims(p) + 1;
ifun = @(p) any(vecnorm(p - swapdim(targ.pos,2,nextdim(p)),2,1) < s_rad, nextdim(p));
med = Medium('c0', targ.c0, 'rho0', targ.rho0, 'pertreg', {{ifun, [targ.c0*targ.c_scat, targ.rho0*targ.rho_scat]}});
%% 
% 

% Show the transducer's impulse response
% figure; plot(xdc.impulse, '.-'); title('Element Impulse Response'); 

% Show the transmit signal - a single point means it's a delta function
% figure; plot(seq.pulse, '.-'); title('Transmit signal');
%  Plot configuration of the simulation

figure; hold on; title('Geometry');
switch "sound-speed"
    case 'sound-speed', imagesc(med, tscan, 'c'  ); colorbar; % show the background medium for the simulation/imaging region
    case 'density',     imagesc(med, tscan, 'rho'); colorbar; % show the background medium for the simulation/imaging region
end
plot(xdc, 'r+', 'DisplayName', 'Elements'); % elements
hps = plot(scan, 'w.', 'MarkerSize', 0.5, 'DisplayName', 'Image'); % the imaging points
switch seq.type % show the transmit sequence
    case 'PW', plot(seq, 3e-2, 'k.', 'DisplayName', 'Tx Sequence'); % scale the vectors for the plot
    otherwise, plot(seq, 'k.', 'DisplayName', 'Tx Sequence'); % plot focal points, if they exist
end
if targ.numScat < 500
    plot(targ, 'k', 'LineStyle', 'none', 'Marker', 'diamond', 'MarkerSize', 5, 'DisplayName', 'Scatterers'); % point scatterers
else
    disp('Info: Too many scatterers to display.');
end
hl = legend(gca, 'Location','bestoutside'); 
set(gca, 'YDir', 'reverse'); % set transducer at the top of the image
%% Simulate the Point Scatterer(s)

 
%% 
% 

% Construct an UltrasoundSystem object, combining all of these properties
us = UltrasoundSystem('xdc', xdc, 'sequence', seq, 'scan', scan, 'fs', 40e6);

% Simulate a point target
% run on CPU to use spline interpolation
switch "Greens"
    case 'FieldII', chd0 = calc_scat_all(us, targ, [1,1], 'device', dev, 'interp', 'cubic'); % use FieldII, 
    case 'FieldII-multi', chd0 = calc_scat_multi(us, targ, [1,1]); % use FieldII, 
    case 'SIMUS'  , us.fs = 4 * us.fc; % to address a bug in MUST where fs must be a ~factor~ of 4 * us.fc
                    chd0 = simus(us, targ, 'periods', 1, 'dims', 3, 'interp', 'cubic'); % use MUST: note that we have to use a tone burst or LFM chirp, not seq.pulse
    case 'Greens' , chd0 = greens(us, targ, [1,1], 'device', dev, 'interp', 'cubic'); % use a Greens function with a GPU if available!
    case 'kWave', chd0 = kspaceFirstOrder(us, med, tscan, 'CFL_max', 0.5, 'PML', [64 128], 'parcluster', 0, 'PlotSim', true); % run locally, and use an FFT friendly PML size
end
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
chd.data = chd.data - mean(chd.data, 1, 'omitnan'); % remove DC 
D = chd.getPassbandFilter(xdc.bw, 25); % get a passband filter for the transducer bandwidth
chd = filter(chd, D); % apply passband filter for transducer bandwidth
if dev, chd = gpuArray(chd); end % move data to GPU
if isreal(chd.data), chd = hilbert(chd, 2^nextpow2(chd.T)); end % apply hilbert on real data

% Run a simple DAS algorithm
switch class(xdc)
    case 'TransducerArray' , scale = 1e-3; % Definitions in millimeters
    case 'TransducerConvex', scale = 1; % Definitions in degrees
end
switch seq.type
    case "VS", 

%% 
% Choose an apodization method for Virtual Source (VS) transmit sequences

        switch "none"
            case 'multiline', apod = multilineApodization(us.scan, us.sequence);
            case 'scanline', apod = scanlineApodization(us.scan, us.sequence);
            case 'translating', apod = translatingApertureApodization(us.scan, us.sequence, us.rx, 19.8*scale);
            case 'aperture-growth', apod = apertureGrowthApodization(us.scan, us.sequence, us.rx, 1.8);
            case 'accept', apod = acceptanceAngleApodization(us.scan, us.sequence, us.rx, 55); 
            case 'none', apod = 1;
        end

%% 
% Choose an apodization method for Full Synthetic Aperture (FSA) transmit sequences

    case "FSA"
        us.sequence.focus = us.tx.positions(); % set the sequence foci to be the location of the transmitters for these profiles
        switch "none"
            case 'translating', apod = translatingApertureApodization(us.scan, us.sequence, us.rx, 19.8*scale);
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

switch "DAS"
    case "DAS"
        b = bfDAS(us, chd, targ.c0, 'apod', apod, 'interp', 'cubic'); % use a vanilla delay-and-sum beamformer
    case "Adjoint"
        b = bfAdjoint(us, chd, targ.c0, 'apod', apod, 'fthresh', -20); % use an adjoint matrix method
    case "Eikonal"
        b = bfEikonal(us, chd, targ, tscan, 'apod', apod, 'interp', 'cubic'); % use the eikonal equation
    case "DAS-direct"
        b = DAS(us, chd, targ.c0, 'apod', apod, 'interp', 'cubic', 'device', dev); % use a specialized delay-and-sum beamformer
end

% show the image
b_im = mod2db(b); % convert to power in dB
figure; imagesc(scan, b_im(:,:,1), [-80, 0] + max(b_im(:))); % display with 80dB dynamic range
colormap gray; colorbar;
%% Scan Convert for a Sector Scan

if ~isa(scan, 'ScanCartesian')
    scanc = scanCartesian(scan); % mind the caps! this gives us a ScanCartesian
    [scanc.nx, scanc.nz] = deal(2^9); % don't change boundaries, but increase resolution
    bc_im = (scanConvert(scan, mod2db(b), scanc)); % interpolate in dB - could also do abs, as long as it's not complex!

    imagesc(scanc, bc_im, [-80, 0] + max(bc_im(bc_im < Inf)));
    colormap gray; colorbar;
end