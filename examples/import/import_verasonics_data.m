% Example script for importing Verasonics Vantage data into QUPS
%
% This example walks through how to import verasonics data into a QUPS
% compatible format. The example data contains multiple transmit foci and
% involves receive aperture multiplexing captured on an L12-3v transducer
% on a Vantage UTA-260-D platform. 
% 
% There is 1 TX per 1 Receive. Since every TX is duplicated, only the first
% of each duplicate pair is used to define the Sequence in QUPS' framework.
% The joining of the receive apertures must be handled outside of the QUPS.
% 
% Verasonics' programming guide is confidential - therefore, the data,
% generation, and usage of the relevant properties are not described.
% 

%% 0) Load the data
fn = 'DATA_L12_3v_MultiFocal.mat'; % filename
dat = load(which(fn)); % req'd: Trans, TX, Receive, RcvData | opt: TW, PData, Resource

%% 1) Construct the UltrasoundSystem piece by piece

% reference sound speed
if isfield(dat, 'Resource'), c0 = dat.Resource.Parameters.speedOfSound;
else,                        c0 = 1540; % default
end

% Transducer (req'd: Trans, c0)
xdc = Transducer.Verasonics(dat.Trans, c0);
lbda = c0 / xdc.fc; % wavelengths

% Scan (req'd: PData, units)
if isfield(dat, 'PData') % import
    switch dat.Trans.units, case "mm", scl = 1; case "wavelengths", scl = lbda; end
    scan = Scan.Verasonics(dat.PData, scl);
else % declare
    pn = xdc.positions;
    scan = ScanCartesian('xb', pn(1,[1,end]), 'zb', [0 40e-3]);
end

% Sequence (req'd: TX, Trans, c0 | opt: TW)
j = 1:2:numel(dat.TX); % data has multiplexing, so we only need every other TX
k = unique([dat.TX(j).waveform]); % select the appropriate transmit waveform
[seq, t0q] = Sequence.Verasonics(dat.TX(j), dat.Trans, dat.TW(k), 'c0', c0); % import

% UltrasoundSystem
us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'scan', scan);

%% 2) Construct the ChannelData
% import 1 frame of data (req'd: RcvData, Receive | opt: Trans)
[chd, fmod] = ChannelData.Verasonics(dat.RcvData, dat.Receive, dat.Trans, 'frames', 1); 

% Fix the time axes (assumes Sequence import was valid)
tlens = - 2*d.Trans.lensCorrection / xdc.fc; % lens correction
chd.t0 = t0q + seq.pulse.t0 + tlens; % beamforming delay corrections

%% 3) Tweak the output to handle any multiplexing 
% We need to handle the rx-aperture multiplexing manually. This particular
% dataset multiplexes the receive aperture for an L12-3v with 192 elements
% connected to a UTA-260-D with 128 channels. The middle 64 elements
% overlap.
%
% We can either truncate the data to a 'left' and 'right' aperture of 96
% elements each or average the 64 overlapping elements.

switch "trunc"
    case "trunc" % truncation method
        [l, r] = deal(1:96, 97:192); % left and right halves
        chd.data = cat(3, chd.data(:,1:2:end,l,:), chd.data(:,2:2:end,r,:)); % truncate and combine

    case "avg" % overlap averaging method
        chd.data = chd.data(:,1:2:end,:,:) + chd.data(:,2:2:end,:,:); % sum
        i = [true(1,64), false(1,64), true(1,64)]; % non-overlap region
        chd.data(:,i,:,:) = 2 * chd.data(:,i,:,:); % double the non-overlapped region (equivalent to halving the overlap)
end








