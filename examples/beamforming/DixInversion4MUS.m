%% Download Full-Synthetic Aperture Dataset
% This example loads data from 
% https://github.com/rehmanali1994/DixInversion4MedicalUltrasound
fls  = fullfile("data","MultistaticDataset"+[1520 1540 1570]'+".mat"); % files to download (~68MB each)
repo = "https://github.com/rehmanali1994/DixInversion4MedicalUltrasound";
urls = repo+"/releases/download/v1.0.0/"+extract(fls,"MultistaticDataset"+digitsPattern(4)+".mat"); % urls
i = find(~arrayfun(@exist, fls)); % check if not there already 
[~,~] = arrayfun(@mkdir, unique(fileparts(fls)), 'uni', 0); % make output folder(s)
arrayfun(@websave, fls(i), urls(i), 'uni', 0); % download if not there

%% Options
kwargs.verbose = true;
kwargs.update = true;
c = 1420 : 1e0 : 1620; % sound speeds for imaging
dc = 0.01; % sound speed upsampling interval
el_str = 4; % element stride - use every el_str transmit elements

%% Convert to QUPS
% load data: set index to on of {1, 2, 3} for each sound speed
dat = load(fls(1)); % load

% Transducer
pitch = mean(diff(dat.rxAptPos(:,1))); % interelement spacing
xdc = TransducerArray('pitch', pitch, 'numel', size(dat.rxAptPos,1));

% Scan (imaging domain)
pb = xdc.bounds(); % boundaries of the Transducer
scan = ScanCartesian('x',pb(1,[1 end]), 'z', [2 35]*1e-3);
[scan.nx, scan.nz] = deal(100, 600);

% Channel Data
fs = 1 / mean(diff(dat.time)); % sampling frequency
chd = ChannelData('data', dat.scat, 't0', dat.time(1), 'fs', fs);

% Medium
grd = ScanCartesian('x', dat.x, 'z', dat.z);
med = Medium.Sampled(grd, dat.C, 'c0', 1540);

% Sequence
tx_elmts = 1 : el_str : xdc.numel; % choose the set of transmits
pn = xdc.positions(); % element positions (foci)
apd = eye(xdc.numel); % element transmit weights (for FSA)

% truncate transmits
if ~isequal(tx_elmts, 1:xdc.numel) % we are downsampling
    pn  = pn(    :, tx_elmts          ); % truncate transmits
    apd = apd(   :, tx_elmts          ); % truncate transmits
    chd = subD(chd, tx_elmts, chd.mdim); % truncate transmits
end

% construct the Sequence
seq = Sequence('type','DV', 'focus',pn, "apd",apd, "c0",med.c0);

% UltrasoundSystem (synthesize)
us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'scan', scan, 'fs', chd.fs);

%% Focused Transmit Image Reconstruction
% Display the B-mode and Coherence Factor images
bmax = gather(mod2db(max(chd.data,[],'all') * us.xdc.numel)); % estimate max image power
figure; clf();
h    = imagesc(us.scan, zeros(us.scan.size), nexttile()); dbr b-mode;
h(2) = imagesc(us.scan, zeros(us.scan.size), nexttile()); dbr corr; clim([0 1]);
clim(h(1).Parent, [-60 0] + bmax); % set axis
colormap(h(2).Parent, 'hot');
ttls = ["B-mode Image", "Coherence Factor"];
title(h(1).Parent, ttls(1));
title(h(2).Parent, ttls(2));

if canUseGPU(), chd = gpuArray(chd); end
chd = hilbert(singleT(chd)); % pre-processing

%% 
if kwargs.verbose, fprintf("\nComputing the coherence factor ... %2.0f%%", 0); tt = tic; end
[cf, bm] = deal(cell(1, numel(c))); % init

% for each sound speed
for i = 1 : numel(c)
    % set the sound speed
    us.seq.c0 = c(i); 

    % image at this sound speed, per rx
    b = us.DAS(chd, 'keep_rx', true);
    rxdim = ndims(b);

    % compute CF image (@cohfac) and B-mode image (@sum)
    cf{i} = gather(cohfac(b, rxdim)); % compute coherence factor over rx
    bm{i} = gather(sum(   b, rxdim)); % compute b-mode image

    % update display
    if kwargs.update
        h(1).CData(:) = mod2db(sum(b, rxdim));
        h(2).CData(:) = abs(cf{i});
        arrayfun(@title, [h.Parent], ttls + " (" + c(i) + " m/s)");
        drawnow limitrate;
    end

    % update progress
    if kwargs.verbose, fprintf("\b\b\b%2.0f%%", floor(100 * i / numel(c))); end
end

% combine sound speed frames into dim 4
bm = cat(4, bm{:});
cf = cat(4, cf{:});
if kwargs.verbose, fprintf(" Done!\n"); toc(tt); end

%% Display all B-modes and Coherence Factor imagesc
figure;
ttlsm = ttls + newline + "Sound Speed : " + c' + " (m/s)";
animate({bm, cf}, h, 'title', ttlsm, 'loop', false);

%%
% Assemble Coherence Factor Curves for Each Depth
cf_avg = mean(cf, 2); % Laterally Average CF
cf_avg = permute(cf_avg, [4 1 2 3]); % -> (C x Z x 1 x 1)

% Estimate Average Sound Speed Upto Depth
c_up = c(1) : dc : c(end); % upsampled sound speed axis
cf_up = interp1(c, cf_avg, c_up); % upsample average coherence factor over sound speeds
c_avg = c_up(argmax(cf_up))'; % choose maximum over sound speeds (Z x 1)

% Smoothing Regularization
reg = 1.0e8; % Regularization Value
I = eye(us.scan.nz); % Identity Matrix
F = I(1:end-3,:) - 3*I(2:end-2,:) + 3*I(3:end-1,:) - I(4:end,:);
c_avg_smoothed = (I+reg*(F'*F)) \ c_avg; % (Z x 1)

% Calculate Local Sound Speed
[z, dz] = deal(us.scan.z', us.scan.dz); % (Z x 1), (1 x 1)
c_local = dz ./ diff( z ./ c_avg_smoothed ); % ((Z - 1) x 1)

% Plot Average Sound Speed vs Depth
figure; imagesc(z, c_up, cf_up);
title('Coherence-Based Average Sound Speed Estimates');
cbar = colorbar(); ylabel(cbar, 'Coherence Factor'); 
hold on;
plot(z, c_avg         , 'k*', 'Linewidth', 2);
plot(z, c_avg_smoothed, 'r', 'Linewidth', 2);
xlabel('Imaging Depth [m]'); 
ylabel('Average Sound Speed [m/s]'); 
axis([min(z),max(z),min(c),max(c)]);
legend('Measured', 'Smoothed', 'Location', 'Northwest');
set(gca, 'YDir', 'normal');

% Plot Local Sound Speed vs Depth
figure; hold on; 
plot((z(2:end)+z(1:end-1))/2, c_local, 'Linewidth', 2);
plot(dat.z, mean(dat.C,2), 'Linewidth', 2);
grid on; grid minor;
title('Local Sound Speed Estimates')
xlabel('Imaging Depth [m]'); 
ylabel('Local Sound Speed [m/s]'); 
axis([min(z),max(z), min(c),max(c)]);
legend('Estimated', 'True', 'Location', 'Northwest');
