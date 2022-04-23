function y = cbeamform(fun, Pi, Pr, Pv, cgrid, x, c, t0, fs, varargin)

% CBEAMFORM - Sound speed aware beamforming
%
% y = CBEAMFORM(fun, Pi, Pr, Pv, cgrid, x, c, t0, fs)
%
% beamforms the data optionally summing across the aperture.
%
% Given a set of pixels, (virtual or plane wave) transmitter locations,
% receiver locations, as well as a datacube equipped with a time,
% transmitter and receiver axis, beamforming the data without summation.
% The data is linearly interpolated at the sample time.
%
% All positions are in vector coordinates.
%
% If the virtual transmitter normal has a fourth component that is 0, this
% indicates that the transmission should be treated as a plane wave
% transmission instead of a virtual source (focused) transmission.
%
% The value of t = 0 must be the time when the peak of the wavefront
% reaches the virtual source location. Because this time must be the same
% for all transmits, the datacube must be stitched together in such a way
% that for all transmits, the same time axis is used.
%
% If the function is 'delays', the output will be the time delay for each
% transmitter-receiver pair. The values for 'x', 't0', and 'fs' will be
% ignored.
%
% Inputs:
%   fun:         algorithm for aperture summation {'DAS'|'SYN'|'BF'|'delays'}
%   Pi:          pixel positions (3 x I)
%   Pr:          receiver positions (3 x N)
%   Pv:          (virtual) transmitter positions (3 x M)
%   cgrid:       sound speed grid structure with
%     cgrid.origin (3 x 1)  origin of the grid (smallest value) in x/y/z
%     cgrid.step   (3 x 1)  step size between grid points in x/y/z.
%                           Use inf/NaN for a sliced dimension.
%     cgrid.size   (3 x 1)  Grid size in x/y/z. Use 1 for a sliced
%                           dimension.
%   x:           datacube of complex sample values (T x N x M)
%   t0:          initial time for the data
%   fs:          sampling frequency of the data
%   c:           speed of sound map used for beamforming. Must match
%                cgrid.size
%
% Outputs:
%   y:           complex pixel values per transmit/channel (I [x N [x M]])
%       I -> pixels, M -> transmitters, N -> receivers, T -> time samples
%
%

% default parameters
VS = true;
odataPrototype = complex(zeros(0, 'like', x));
interp_type = 'linear';
isType = @(c,T) isa(c, T) || isa(c, 'gpuArray') && strcmp(classUnderlying(c), T);
if isType(x, 'single')
    idataType = 'single';
else
    idataType = 'double';
end
if any(cellfun(@(c) isType(c, 'single'), {Pi, Pr, Pv}))
    posType = 'single';
else
    posType = 'double';
end
if gpuDeviceCount && any(cellfun(@(v)isa(v, 'gpuArray'), {Pi, Pr, Pv, x}))
    device = -1; % default GPU
else
    device = 0; % CPU
end

% optional inputs
nargs = int64(numel(varargin));
n = int64(1);
while (n <= nargs)
    switch varargin{n}
        case 'plane-waves'
            VS = false;
        case 'virtual-source'
            VS = true;
        case 'output-prototype'
            n = n + 1;
            odataPrototype = varargin{n};
        case 'output-precision'
            n = n + 1;
            odataPrototype = complex(zeros(0, varargin{n}));
        case 'position-precision'
            n = n + 1;
            posType = varargin{n};
        case 'input-precision'
            n = n + 1;
            idataType = varargin{n};
        case 'device'
            n = n + 1;
            device = int64(varargin{n});
        case 'interp'
            n = n + 1;
            interp_type = varargin{n};
        otherwise
            error('Unrecognized option');
    end
    n = n + 1;
end

% defaults
if nargin < 9 || isempty(c), c = 1540; end
if (nargin < 8 || isempty(fs))
    if isvector(t0)
        fs = mean(diff(t0)); % find the sampling frequency
        t0 = min(t0); % extract start time
    elseif strcmp(fun, 'delays')
        fs = [];
    else
        error('Undefined sampling rate.');
    end
end

% parse inputs
parse_inputs();

% store this option input for convenience
per_cell = {'UniformOutput', false};


% cast constant data on CPU
dprototype = complex(zeros(0, idataType));
pprototype = zeros(0, posType);

x = cast(x, 'like', dprototype);
[tmp] = cellfun(@(p)cast(p, 'like', pprototype), ...
    {Pi, Pr, Pv}, per_cell{:});
[Pi, Pr, Pv] = deal(tmp{:});

% expand all inputs to 3D
expand_inputs();

% data sizing
[T, M, N] = size(x);
[~, I] = size(Pi);

% cast to integers
tmp = cellfun(@uint64, {T,M,N,I}, per_cell{:});
[T, M, N, I] = deal(tmp{:});

% get output data sizing
switch fun
    case 'DAS'
        osize = {I, uint64(1)}; % delay and sum (tx/rx)
    case 'SYN'
        osize = {I, N}; % beamform and sum the transmit aperture only
    case {'BF', 'delays'}
        osize = {I, N, M}; % beamform only | times only
end

% time vector
t = t0 + cast(colon(0,T-1).', 'like', t0) ./ fs;

% get sound speed map dimension info
nsdims = cgrid.size ~= 1; % non-singleton dimensions

% ensure grid sizing is the same
cstep = unique(cgrid.step(nsdims));
assert(numel(cstep) == 1, ...
    'The cgrid must have equally sized steps in all non-singleton dimensions.'...
    );

% convert positions to sound speed (1-based) grid coordinates
Pvc = (Pv - cgrid.origin) ./ cgrid.step + 1; % 3 x M
Prc = (Pr - cgrid.origin) ./ cgrid.step + 1; % 3 x N
Pic = (Pi - cgrid.origin) ./ cgrid.step + 1; % 3 x I

% trim singleton dimension(s)
[Pvc, Prc, Pic] = dealfun(@(n) n(nsdims,:), Pvc, Prc, Pic);

% send const data to all workers maybe
hcp = gcp('nocreate');
if(isempty(hcp))
    pvc.Value = Pvc ;
    prc.Value = Prc ;
    pic.Value = Pic'; % transpose here for convenience
    cnorm.Value = c ./ cstep;
else
    pvc   = parallel.pool.Constant(Pvc );
    prc   = parallel.pool.Constant(Prc );
    pic   = parallel.pool.Constant(Pic'); % transpose here for convenience
    cnorm = parallel.pool.Constant(c ./ cstep);
end

% get one-way delays within the field then generate samplers, using
% reduced dimensions ([M|N] x {Cx x Cy x Cz})
[tx_samp, rx_samp] = deal(cell([M,1]),cell([N,1]));
gi_opts = {'linear', 'none'};
tt = tic; fprintf('\n');
parfor m = 1:M
    fprintf('tx %i\n', m);
    [tau_map_tx] = msfm(squeeze(cnorm.Value), double(pvc.Value(:,m))); %#ok<PFBNS>
    tx_samp{m} = griddedInterpolant(tau_map_tx,gi_opts{:}); %#ok<PFBNS>
end
parfor n = 1:N
    fprintf('rx %i\n', n);
    [tau_map_rx] = msfm(squeeze(cnorm.Value), double(prc.Value(:,n))); %#ok<PFBNS>
    rx_samp{n} = griddedInterpolant(tau_map_rx,gi_opts{:}); %#ok<PFBNS>
end
fprintf('\nEikonal time delays completed in %0.3f seconds.\n', toc(tt));

if device

    % currently selected gpu device handle
    g = gpuDevice();

    % reselect gpu device and copy over inputs if necessary
    if device > 0 && g.Index ~= device
        inps = {Pi, Pr, Pv, x, t0, fs, c, cgrid};
        oclass = cellfun(@class, inps, per_cell{:});
        tmp = cellfun(@gather, inps, per_cell{:});
        g = gpuDevice(device);
        tmp = cellfun(@(inp, cls)cast(inp, cls), tmp, oclass, per_cell{:});
        [Pi, Pr, Pv, x, t0, fs, c, cgrid] = deal(tmp{:});
    end

    % move data to gpu
    [x, t] = deal(gpuArray(x), gpuArray(t));

    y = zeros([osize{:}], 'like', odataPrototype);

    % create the image
    tt = tic;
    switch fun
        case 'DAS'
            for m = 1:M
                tau_tx = gpuArray(tx_samp{m}(pic.Value)); % get tx delay (I x 1 x 1)
                xm = x(:,:,m); % slice tx data
                tt2 = tic;
                for n = 1:N
                    tau_rx = (rx_samp{n}(pic.Value)); % get rx delay (I x 1 x 1)
                    xmn = xm(:,n); % slice tx data
                    ys = interpn(t, xmn, tau_tx + tau_rx, interp_type, 0); % sample
                    y(:) = y + ys; % combine for output (I x 1 x 1)
                end
                fprintf('tx %i - %0.6f\n', m, toc(tt2));
            end

        case 'SYN'
            for m = 1:M
                tau_tx = gpuArray(tx_samp{m}(pic.Value)); % get tx delay (I x 1 x 1)
                xm = x(:,:,m); % slice tx data
                for n = 1:N
                    fprintf('tx/rx %i/%i\n', m, n);
                    tau_rx = (rx_samp{n}(pic.Value)); % get rx delay (I x 1 x 1)
                    xmn = xm(:,n); % slice rx data
                    ys = interpn(t, xmn, tau_tx + tau_rx, interp_type, 0); % sample
                    y(:,n) = y(:,n) + ys; % combine for output (I x N)
                end
            end

        case 'BF'
            for m = 1:M
                tau_tx = gpuArray(tx_samp{m}(pic.Value)); % get tx delay (I x 1 x 1)
                xm = x(:,:,m); % slice tx data
                for n = 1:N
                    fprintf('tx/rx %i/%i\n', m, n);
                    tau_rx = (rx_samp{n}(pic.Value)); % get rx delay (I x 1 x 1)
                    xmn = xm(:,n); % slice rx data
                    ys = interpn(t, xmn, tau_tx + tau_rx, interp_type, 0); % sample
                    y(:,n,m) = ys; % combine for output (I x N x M)
                end
            end

        case 'delays'
            for m = 1:M
                tau_tx = gpuArray(tx_samp{m}(pic.Value));
                for n = 1:N
                    fprintf('tx/rx %i/%i\n', m, n);
                    tau = tau_tx + (rx_samp{n}(pic.Value)); % time delay I x 1 x 1
                    y(:,n,m) = tau; % combine for output (I x N x M)
                end
            end
    end
    fprintf('\nBeamforming completed in %0.3f seconds.\n',toc(tt));
else

    % accumulate over transmitters / receivers
    y = zeros([osize{:}], 'like', odataPrototype);

    % create the image
    tt = tic;
    switch fun
        case 'DAS'
            hcp = gcp('nocreate');
            for m = 1:M
                tau_tx = tx_samp{m}(pic.Value); % get tx delay (I x 1 x 1)
                xm = x(:,:,m); % slice tx data
                hcp.ticBytes()
                tt2 = tic;
                parfor n = 1:N
                    tau_rx = rx_samp{n}(pic.Value); %#ok<PFBNS> % get rx delay (I x 1 x 1)
                    xmn = xm(:,n); % slice tx data
                    ys = interpn(t, xmn, tau_tx + tau_rx, interp_type, 0); % sample
                    y = y + cast(ys, 'like', odataPrototype); % combine for output (I x 1 x 1)
                end
                fprintf('tx %i - %0.6f\n', m, toc(tt2));
                hcp.tocBytes()
            end

        case 'SYN'
            for m = 1:M
                tau_tx = tx_samp{m}(pic.Value); % get tx delay (I x 1 x 1)
                xm = x(:,:,m); % slice tx data
                parfor n = 1:N
                    fprintf('tx/rx %i/%i\n', m, n);
                    tau_rx = rx_samp{n}(pic.Value); %#ok<PFBNS>% get rx delay (I x 1 x 1)
                    xmn = xm(:,n); % slice rx data
                    ys = interpn(t, xmn, tau_tx + tau_rx, interp_type, 0); % sample
                    y(:,n) = y(:,n) + cast(ys, 'like', odataPrototype); % combine for output (I x N)
                end
            end

        case 'BF'
            for m = 1:M
                tau_tx = tx_samp{m}(pic.Value); % get tx delay (I x 1 x 1)
                xm = x(:,:,m); % slice tx data
                parfor n = 1:N
                    fprintf('tx/rx %i/%i\n', m, n);
                    tau_rx = rx_samp{n}(pic.Value); %#ok<PFBNS>% get rx delay (I x 1 x 1)
                    xmn = xm(:,n); % slice rx data
                    ys = interpn(t, xmn, tau_tx + tau_rx, interp_type, 0); % sample
                    y(:,n,m) = cast(ys,'like', odataPrototype); % combine for output (I x N x M)
                end
            end

        case 'delays'
            for m = 1:M
                tau_tx = tx_samp{m}(pic.Value);
                parfor n = 1:N
                    fprintf('tx/rx %i/%i\n', m, n);
                    tau = tau_tx + rx_samp{n}(pic.Value); %#ok<PFBNS>% time delay I x 1 x 1
                    y(:,n,m) = cast(tau, 'like', odataPrototype); % combine for output (I x N x M)
                end
            end
    end
    fprintf('\nBeamforming completed in %0.3f seconds.\n',toc(tt));
end

    function parse_inputs()
        % * Inputs:
        % *  y:           complex pixel values per transmit/channel (M x N x I)
        % *  Pi:          pixel positions (3 x I)
        % *  Pr:          receiver positions (3 x N)
        % *  Pv:          (virtual) transmitter positions (3 x M)
        % *  Nv:          (virtual) transmitter normal (3 x M)
        % *  x:           datacube of complex sample values (T x M x N)
        % *  t0:          initial time for the data
        % *  fs:          sampling frequency of the data
        % *  cinv:        inverse of the speed of sound used for beamforming
        % *
        % * I -> pixels, M -> transmitters, N -> receivers, T -> time samples

        % check the function call is valid
        switch fun
            case {'DAS', 'SYN', 'BF', 'delays'}
            otherwise, error('Invalid beamformer.');
        end

        % check shape
        Pi = modSize(Pi);
        Pv = modSize(Pv);
        Pr = modSize(Pr);

        % ensure coordinates in the 1st dimension
        function x = modSize(x)
            if size(x,1) <= 4
                % seems legit
            elseif size(x,2) <= 4 && (size(x,1) == 1 || size(x,1) > 4)
                x = x.';
            else % ambiguous ... we'll see what happens
            end
        end
    end

    function expand_inputs()

        % data sizing
        if strcmp(fun, 'delays')
            [T, M, N] = deal(0);
            [~, I] = size(Pi);
            [~, Mv] = size(Pv);
            [~, Nr] = size(Pr);

            % expand vectors
            M = max(Mv, Mnv);
            N = Nr;
            if Mv  == 1, Pv = repmat(Pv,1,M); Mv  = M; end
        else
            [T, M, N] = size(x);
            [~, I] = size(Pi);
            [~, Mv] = size(Pv);
            [~, Nr] = size(Pr);

            % expand vectors
            if Mv  == 1 , Pv = repmat(Pv,1,M); Mv  = M; end
            if Nr  == 1 , Pr = repmat(Pr,1,N); Nr  = N; end
        end

        % check sizing
        assert((M == Mv), 'Inconsistent transmitter data size.');
        assert(N == Nr, 'Inconsistent receiver data size.');

        % expand to 3D: 1D -> x, 2D -> x,z, 4D -> x/w, y/w, z/w
        Pi = modDim(Pi);
        Pv = modDim(Pv);
        Pr = modDim(Pr);

        % move to 3D
        function x = modDim(x)
            if size(x, 1) == 1
                x = [x(1,:); zeros([2, size(x,2)], 'like', x)];
            elseif size(x,1) == 2
                x = [x(1,:); zeros([1, size(x,2)], 'like', x); x(2,:)];
            elseif size(x,1) == 3
            elseif size(x,1) == 4
                x = x(1:3,:) ./ x(4,:); % for projective coordinates, project
            else
                error('Improper coordinate dimension.');
            end
        end
    end
end

%#ok<*ASGSL>
%#ok<*NBRAK>



