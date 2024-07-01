function [y, k, pre_args, post_args] = das_spec(fun, Pi, Pr, Pv, Nv, x, t0, fs, c, varargin)
% DAS_SPEC - Specialized delay-and-sum beamformer
% 
% y = DAS_SPEC(fun, Pi, Pr, Pv, Nv, x, t0, fs, c) beamforms the given data
% by computing and applying beamforming delays based on an assumed sound
% speed, receiver positions, and (virtual) positions transmit positions. 
% Optionally sum across the receive aperture and/or transmit aperture 
% (default).
% 
% Inputs:
%  fun -         algorithm for aperture summation {'DAS'|'SYN'|'BF'|'delays'}
%  Pi -          pixel positions (3 x [I]) == (3 x I1 x I2 x I3)
%  Pr -          receiver positions (3 x N)
%  Pv -          (virtual) transmit foci (3 x M)
%  Nv -          plane-wave transmit normal (3 x M)
%  x -           datacube of complex sample values (T x N x M x F x ...)
%  t0 -          initial time for the data (scalar)
%  fs -          sampling frequency of the data (scalar)
%  c -           sound speed used for beamforming (scalar)
%
% Outputs: 
%  y -           beamformed data (I1 x I2 x I3 x [1|N] x [1|M] x F x ...)
% 
% Where T -> time samples, N -> receivers, M -> transmits, F -> data frames
% I = [I1 x I2 x I3] -> pixels
% 
% fun must be one of {'DAS'*|'SYN'|'BF'|'delays'}:
%   DAS - sum across both apertures
%   SYN - sum across the transmi apeture only
%   BF  - do not sum
%   delays - return time delays
%
% If the function is 'delays', the output will be the time delay for each
% transmitter-receiver pair. The values for 'x', 't0', and 'fs' will be
% ignored.
% 
% y = DAS_SPEC(..., 'plane-waves', ...) uses a plane-wave delay model 
% instead of a virtual-source delay model (default).
% 
% y = DAS_SPEC(..., 'diverging-waves', ...) specifies a diverging-wave 
% delay model, in which the time delay is always positive rather than
% negative prior to reaching the virtual source in the virtual-source delay
% model (default).
% 
% In a plane-wave model, time t == 0 is when the wave passes through the 
% origin of the coordinate system. Pv is the origin of the coordinate
% system (typically [0;0;0]) and Nv is the normal vector of the plane wave.
% 
% In a virtual-source model, time t == 0 is when the wavefront (in theory) 
% passes through the focus. For a full-synthetic-aperture (FSA) 
% acquisition, Pv is the position of the transmitting element. In a focused
% or diverging wave transmit, Pv is the focus. For a focused wave, Nv is
% used to determine whether the theoretical wavefront has reached the focal
% point. For a diverging wave, Nv is ignored. 
% 
% y = DAS_SPEC(..., 'apod', apod, ...) applies apodization across the image
% and data dimensions. apod must be able to broadcast to size 
% I1 x I2 x I3 x N x M
% 
% y = DAS_SPEC(..., 'device', device, ...) forces selection of a gpuDevice.
% device can be one of {0, -1, +n}:
%    0 - use native MATLAB calls to interp1
%   -1 - use a CUDAKernel on the current GPU
%   +n - reset gpuDevice n and use a CUDAKernel (caution: this will clear
%        all of your gpuArray variables!
% 
% y = DAS_SPEC(..., 'interp', method, ...) uses method for the underlying
% intepolation. On the gpu, method can be one of 
% {"nearest","linear"*,"cubic","lanczos3"}. In MATLAB, support is
% determined by the interp1 function. 
% 
% [y, k, PRE_ARGS, POST_ARGS] = DAS_SPEC(...) when the CUDA ptx is used returns
% the parallel.gpu.CUDAKernel k as well as the arguments for calling the
% data PRE_ARGS and POST_ARGS. The kernel can then be called per frame f as 
%
%     y{f} = k.feval(PRE_ARGS{:}, x{f}, POST_ARGS{:}); 
%
% The data x{f} must be in dimensions T x N x M. If x{f} is a gpuArray, 
% it must have the same type as was used to create the
% parallel.gpu.CUDAKernel k. This is useful for processing many identical
% frames with minimal overhead.
%
% NOTE: if the input data is smaller than was used to create the
% parallel.gpu.CUDAKernel k, an illegal address error may occur, requiring
% MATLAB to be restarted!
%
% See also ULTRASOUNDSYSTEM/DAS CHANNELDATA/SAMPLE INTERP1

% TODO: switch to kwargs struct to support arguments block
% default parameters
VS = true; % whither plane wave
DV = false; % whither diverging wave
interp_type = 'linear'; 
apod = cell(1,0);
isType = @(c,T) isa(c, T) || isa(c, 'gpuArray') && strcmp(classUnderlying(c), T);
if isType(x, 'single')
    idataType = 'single';
elseif isType(x, 'double')
    idataType = 'double';
elseif isType(x, 'halfT')
    idataType = 'halfT'; 
else
    idataType = 'double'; % default
end
if gpuDeviceCount && any(cellfun(@(v)isa(v, 'gpuArray'), {Pi, Pr, Pv, Nv, x}))
    device = -1; % GPU
else
    device = 0; % CPU
end
fmod = 0; % modulation frequency

% optional inputs
nargs = int64(numel(varargin));
n = int64(1);
while (n <= nargs)
    switch varargin{n}
        case 'plane-waves'
            VS = false;
        case 'virtual-source'
            VS = true;
        case 'diverging-waves'
            DV = true;
        case 'focused-waves'
            DV = false;
        case 'input-precision'
            n = n + 1;
            idataType = varargin{n};
        case 'device'
            n = n + 1;
            device = int64(varargin{n});
        case 'interp'
            n = n + 1;
            interp_type = varargin{n};
        case 'apod'
            n = n + 1; 
            apod{end+1} = varargin{n};
        case 'modulation'
            n = n + 1;
            fmod = varargin{n};
        otherwise
            error('Unrecognized option');
    end
    n = n + 1;
end

% defaults
if nargin < 9 || isempty(c), c = 1540; end
if (nargin < 8 || isempty(fs))
    if isvector(t0)
        fs = mean(diff(t0,1,1)); % find the sampling frequency
        t0 = min(t0,[],1); % extract start time
    elseif strcmp(fun, 'delays')
        fs = [];
    else
        error('Undefined sampling rate.');
    end
end

% ensure non-empty
if isempty(apod), apod = {1}; end

% init
[k, pre_args, post_args] = deal({});

% get speed of sound inverse
cinv = 1./c;

% parse inputs
parse_inputs();

% store this option input for convenience
per_cell = {'UniformOutput', false};

% move replicating dimensions of the data to dims 6+
% TODO: support doing this with non-scalar t0 definition
x = permute(x, [1:3,(max(3,ndims(x))+[1,2]),4:ndims(x)]); % (T x N x M x 1 x 1 x F x ...)
fsz = size(x, 6:max(6,ndims(x))); % frame size: starts at dim 6

% get devices: must have package, device, and kernel
gdev = (exist('gpuDeviceCount','file') && gpuDeviceCount() && device <= gpuDeviceCount() && logical(exist('bf.ptx', 'file'))); ... % PTX track available
odev = (exist('oclDeviceCount','file') && oclDeviceCount() && device <= oclDeviceCount() && logical(exist('bf.cl' , 'file'))); ... OpenCL kernel available
if odev % check for data type support
    d = oclDevice();
    odev = ~isempty(d) ... % must have a device selected
        && (((idataType == "double") && d.SupportsDouble) ... double support 
        ||   (idataType == "single") ...                      single always supported)
        ||  ((idataType == "halfT" ) && d.SupportsHalf  )); % half support 
end 

% dispatch
if device && (gdev || odev)

    % ocl requires exactly this many separable apodization matrices
    AMAX = 1; % can be modified by modifying the kernel
    if ~gdev && numel(apod) > AMAX, error("QUPS:das_spec:OpenCLApodizationSupport","OpenCL only supports up to "+AMAX+" apodization matrices."); end
    if ~gdev, [apod{end+1:AMAX}] = deal(1); end
    
    % warn if non-linear interp was requested
    switch interp_type
        case "nearest", flagnum = 0;
        case "linear",  flagnum = 1;
        case "cubic",   flagnum = 2;
        case "lanczos3",flagnum = 3;
        otherwise
            error('QUPS:das_spec:UnrecognizedInput', "Unrecognized interpolation of type " ...
                + interp_type ...
                + ": must be one of {'nearest', 'linear', 'cubic', 'lanczos3'}.");
    end

    % modify flag to choose whether to store rx / tx
    keep_rx = ismember(fun, {'SYN', 'BF'});
    keep_tx = ismember(fun, {'MUL', 'BF'});
    flagnum = flagnum + 8 * keep_rx + 16 * keep_tx;

    % load the kernel
    src_folder = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    
    switch idataType
        case 'halfT',  postfix = 'h';
        case 'single', postfix = 'f';
        case 'double', postfix = '';
    end
    
    % reselect gpu device and copy over inputs if necessary
    if device > 0 && gdev && getfield(gpuDevice(), 'Index') ~= device
        inps = {Pi, Pr, Pv, Nv, x, t0, fs, cinv};
        oclass = cellfun(@class, inps, per_cell{:}); 
        tmp = cellfun(@gather, inps, per_cell{:});
        % g = gpuDevice(device); 
        tmp = cellfun(@(inp, cls)cast(inp, cls), tmp, oclass, per_cell{:});
        [Pi, Pr, Pv, Nv, x, t0, fs, cinv] = deal(tmp{:});
    end

    % send constant data to GPU
    ptypefun = str2func(idataType); % function to cast types
    if gdev, ptypefun = @(x) gpuArray(ptypefun(x)); end
    dtypefun = @(x) complex(ptypefun(x));
    if idataType == "halfT"
        if gdev, ptypefun = @(x) gpuArray(single(x)); end
        if odev, ptypefun = @(x)          single(x) ; end
        dtypefun = @(x) complex(halfT(x));
    end
    [x, apod{:}] = dealfun(dtypefun, x, apod{:});
    [Pi, Pr, Pv, Nv, cinv] = dealfun(ptypefun, Pi, Pr, Pv, Nv, cinv);

    % expand all inputs to 3D
    expand_inputs();
    
    % data sizing
    [T, N, M] = size(x, 1:3);
    Isz = size(Pi, 2:4); % I1 x I2 x I3 == I
    I = prod(Isz);
    S = numel(apod);

    % get stride for apodization and sound speed
    Icsz = size(cinv,1:5)';
    Iasz = cell2mat(cellfun(@(a){size(a,1:5)'}, apod));
    cstride = [[        1; cumprod(Icsz(1:end-1  ),1)] .* (Icsz ~= 1); 0                                    ];
    astride = [[ones(1,S); cumprod(Iasz(1:end-1,:),1)] .* (Iasz ~= 1); sub(cumsum([0 prod(Iasz,1)],2),1:S,2)];
    
    % get kernel and frame sizing
    switch fun
        case 'DAS', osize = {1,1}; % delay and sum over tx/rx
        case 'SYN', osize = {N,1}; % beamform and sum the transmit aperture only
        case 'MUL', osize = {1,M}; % beamform and sum the receive aperture only
        case 'BF',  osize = {N,M}; % beamform only 
        case 'delays', osize = {N,M}; % delays only
    end
    
    % get output data type
    switch fun
        case {'DAS', 'SYN', 'BF', 'MUL'}
            obufproto = x;
        case {'delays'}
            obufproto = real(x([]));
    end
    
    if gdev
        % reference selected device 
        g = gpuDevice();

        % load kernel
        k = parallel.gpu.CUDAKernel(...
            'bf.ptx',...
            fullfile(src_folder, 'bf.cu'),...
            ['DAS', postfix]); % use the same kernel, just modify the flag here.

        % constant arg type casting
        tmp = cellfun(@uint64, {T,M,N,I,Isz,S}, per_cell{:});
        [T, M, N, I, Isz, S] = deal(tmp{:});

        % set constant args
        k.setConstantMemory('QUPS_I', I); % gauranteed
        try k.setConstantMemory('QUPS_T', T); end %#ok<TRYNC> % if not const compiled with ChannelData
        try k.setConstantMemory('QUPS_M', M, 'QUPS_N', N, 'QUPS_S', S, 'QUPS_VS', VS, 'QUPS_DV', DV, 'QUPS_I1', Isz(1), 'QUPS_I2', Isz(2), 'QUPS_I3', Isz(3)); end %#ok<TRYNC> % if not const compiled

        % kernel sizes
        k.ThreadBlockSize(1) = k.MaxThreadsPerBlock; % threads per block
        k.GridSize(1) = min([...
            g.MaxThreadBlockSize(1), ... max cause device reqs
            ceil(I / k.ThreadBlockSize(1))... max cause number of pixels
            ceil(g.AvailableMemory / (2^8) / (prod(osize{:}) * k.ThreadBlockSize(1))),... max cause GPU memory reqs (empirical)
            ]); % blocks per frame

    elseif odev
        
        % load kernel refernece
        k = oclKernel(which('bf.cl'), 'DAS');
        
        % select device - no risk of resetting anything
        if device > 0, oclDevice(device); end 
        k.Device = oclDevice(); 

        % set precision
        switch idataType
            case "single", prc = 32;
            case "double", prc = 64;
            otherwise, error("Not implemented :(");
        end
        k.defineTypes(repmat(idataType,[1,3]), ["V","T","U"]); % time / data / position

        % set constant args
        k.macros = [k.macros, ("QUPS_" + ["I","T","M","N","S","VS","DV","I1","I2","I3"] + "=" + [I,T,M,N,S,VS,DV,Isz])];
        k.macros = [k.macros, ("QUPS_" + ["F"] + "=" + [1])]; % unused ...
        k.macros = [k.macros, ("QUPS_BF_" + ["FLAG"] + "=" + [flagnum])];
        k.macros(end+1) = "QUPS_PRECISION="+prc;
        k.opts = ["-cl-mad-enable"];%, "-cl-fp32-correctly-rounded-divide-sqrt", "-cl-opt-disable"];

        % compile kernel
        k.build();
    
        % set kernel size
        k.ThreadBlockSize(1) = min(I, k.MaxThreadsPerBlock);
        k.GridSize(1) = ceil(I / k.ThreadBlockSize(1));

        % expand inputs to 4D - in  OpenCL, all 3D vectors interpreted are 4D underlying
        [Pi(4,:), Pr(4,:), Pv(4,:), Nv(4,:)] = deal(0);
    end

    % vectorize (combine) all apodization matrices for CUDA
    if gdev
        [apod{:}] = dealfun(@(x) x(:), apod{:});
        apod = {dtypefun(cat(1, apod{:}))};
    end

    
    % allocate output data buffer
    osize = cat(2, {I}, osize);
    osize = cellfun(@uint64, osize, per_cell{:});
    yg = zeros([osize{:}], 'like', obufproto);
    
    % for half types, alias the weights/data, recast the positions as
    % single
    if idataType == "halfT"
        [Pi, Pr, Pv, Nv] = dealfun(@single, Pi, Pr, Pv, Nv);
        [yg, apod{:}, x] = dealfun(@(x)getfield(alias(x),'val'), yg, apod{:}, x);
    end

    % combine timing info with the transmit positions
    Pv(4,:) = t0;
    
    % partition data per frame
    F = prod(fsz);

    % for each data frame, run the kernel
    switch fun
        case {'DAS','SYN','BF','MUL'}
            if ~isreal(obufproto), yg = complex(yg); end % force complex if x is
            % I [x N [x M]] x 1 x 1 x {F x ...}
            for f = F:-1:1 % beamform each data frame
                y{f} = k.feval(yg, Pi, Pr, Pv, Nv, apod{:}, cinv, [cstride, astride], x(:,:,:,f), flagnum, [fs, fmod]);
            end
            % unpack frames
            y = cat(6, y{:});
        case {'delays'}
            y = k.feval(yg, Pi, Pr, Pv, Nv, cinv(1)); ... TODO: fix the bf.cu kernel
    end
    
    % reshape output and truncate garbage
    y = reshape(y, [Isz, size(y,2), size(y,3), fsz]); % I1 x I2 x I3 x [1|N] x [1|M] x F x ...
    
    % if it's a half type, make an aliased halfT
    if idataType == "halfT", y_ = alias(halfT([])); y_.val = y; y = y_; end

    % save the pre/post arguments if requested
    if nargout > 1
        pre_args = [{yg, Pi, Pr, Pv, Nv}, apod, {cinv, [cstride, astride]}];
        post_args = {flagnum, [fs, fmod]};
    end
else
    
    % cast constant data on CPU
    typefun = str2func(idataType);
    ptypefun = @(x) typefun(x);
    dtypefun = @(x) complex(ptypefun(x));
    [x] = dealfun(dtypefun, x);
    [Pi, Pr, Pv, Nv] = dealfun(ptypefun, Pi, Pr, Pv, Nv);
    
    % expand all inputs to 3D
    expand_inputs();
    
    % data sizing
    [T, N, M] = size(x, 1:3);
    Isz = size(Pi,2:4);
    I = prod(Isz);
    
    % cast to integers
    tmp = cellfun(@uint64, {T,N,M,I}, per_cell{:});
    [T, N, M, I] = deal(tmp{:});
    
    % permute to compatible casting dimensions
    Pr = swapdim(Pr, 2, 5); % 3 x 1 x 1 x 1 x N
    Pv = swapdim(Pv, 2, 6); % 3 x 1 x 1 x 1 x 1 x M
    Nv = swapdim(Nv, 2, 6); % 3 x 1 x 1 x 1 x 1 x M
    t0 = swapdim(t0(:), 1, 5); % 1 x 1 x 1 x 1 x M
    t0 = repmat(t0,[ones(1,4), M / size(t0,5)]); % explicit broadcast
    
    % transmit sensing vector
    rv = Pi - Pv; % 3 x I1 x I2 x I3 x 1 x M
    if VS % virtual source delays
        if DV, s = 1; else, s = sign(sum(rv .* Nv,1)); end % diverging or focused
        dv = vecnorm(rv, 2, 1) .* s;
    else % plane-wave delays
        dv = sum(rv .* Nv, 1);
    end % 1 x I1 x I2 x I3 x 1 x M
    
    % receive sensing vector
    dr = vecnorm(Pi - Pr, 2, 1); % 1 x I1 x I2 x I3 x N x 1
    
    % bring to I1 x I2 x I3 x N x M == [I] x N x M
    dv = reshape(dv, size(dv, [2:ndims(dv),1])); 
    dr = reshape(dr, size(dr, [2:ndims(dr),1]));
    
    % apply modulation frequency
    if fmod       
        t = t0 + double((0:T-1).') ./ fs; % time vector
        x = x .* exp(2i*pi*fmod.*t);
    end
    
    % temporal packaging function
    pck = @(x) num2cell(x, [1:3,6:ndims(x)]); % pack for I, F
    osz = [Isz,1,1,fsz]; % interp1 output size (unsqueeze)
    S = numel(apod);

    switch fun
        case 'delays'
            y = cinv .* (dv + dr);
            
        case 'DAS'
            y = dtypefun(zeros([Isz, 1, 1]));
            dvm = num2cell(dv, [1:3]); % ({I} x 1 x M)
            drn = num2cell(dr, [1:3]); % ({I} x N x 1)
            xmn = shiftdim(num2cell(x,  [1,6:ndims(x)]),1); % ({T} x N x M x 1 x 1 x {F x ...})
            cinvmn = shiftdim(num2cell(cinv, [1:3]),3); % ({I} x [1|N] x [1|M])
            for s = S:-1:1
                amn = swapdim(num2cell(apod{s}, [1:3]),[4 5]); % ({I} x [1|N] x [1|M])
                asnm(s,:,:) = swapdim(repmat(amn, [double([N M]) ./ size(amn)]),[1,2],[2,3]); % broadcast to 1 x [1|N] x M
            end
            parfor m = 1:M
                yn = dtypefun(zeros([Isz, 1, 1]));
                asn = asnm(:,:,m);
                if size(cinvmn,5) == 1, cinvn = cinvmn; else, cinvn = cinvmn(:,m); end % ({I} x [1|N] x 1)
                for n = 1:N
                    if isscalar(cinvn), cinv_ = cinvn{1}; else, cinv_ = cinvn{n}; end % ({I} x 1 x 1)
                    % time delay (I x 1 x 1)
                    tau = cinv_ .* (dvm{m} + drn{n}) - t0(m); %#ok<PFBNS>

                    % extract apodization
                    % a = an{min(n,Na)}; % amn(n,m) || amn(1,m)
                    a = asn{end,n}; for s=1:S-1, a=a.*asn{s,n}; end % reduce apodizations
                    
                    % sample and output (I x 1 x 1)
                    yn = yn + a .* (...
                        reshape(interp1(xmn{n,m}, 1 + tau * fs, interp_type, 0), osz) ...
                       );
                end
                y = y + yn;
            end
            y = reshape(y, [Isz, 1, 1, fsz]);
        case 'MUL'
            y = dtypefun(zeros([Isz, 1, M]));
            drn = num2cell(dr, [1:3]); % ({I} x  N  x 1)
            xn  = num2cell(swapdim(x,3,5),  [1,5,6:ndims(x)]); % ({T} x N {x M x F x ...)
            cinvn = num2cell(cinv, [1:3 5]); % ({I} x N x {M})
            for s = S:-1:1
                an = num2cell(apod{s}, [1:3, 5]);
                asn(s,:) = repmat(an, double([ones(1,3),N]) ./ [ones(1,3), size(an,4)]); % broadcast to ({I} x N x {M})
            end
            parfor n = 1:N
                % time delay (I x 1 x M)
                if isscalar(cinvn), cinv_ = cinvn{1}; else, cinv_ = cinvn{n}; end % ([1|I] x 1 x [1|M])
                tau = cinv_ .* (dv + drn{n}) - t0; 

                % extract apodization
                % am = repmat(an{n}, double([ones(1,4), M]) ./ [ones(1,4), size(an{n},5)]);
                as = asn(:,n); % init
                a = as{end}; for s = 1:S-1, a = a .* as{s}; end % accumulate
                a = repmat(a, [ones(1,4), double(M ./ size(a,5))]); % broadcast

                % sample and output (I x 1 x M)
                ym = (cellfun(...
                    @(x, tau, a) ...
                    a .* reshape(interp1(x,  1 + tau * fs, interp_type, 0),osz), ...
                    num2cell(xn{n},[1,6:ndims(xn{n})]), pck(tau), pck(a), per_cell{:})) ...
                    ; %#ok<PFBNS>
                ym = reshape(cat(4,ym{:}), [Isz, size(ym,4:ndims(ym)), fsz]);
                y = y + ym;
            end
            y = reshape(y, [Isz, 1, M, fsz]);
        case 'SYN'
            y = dtypefun(zeros([Isz, N, 1]));
            dvm = num2cell(dv, [1:3]); % ({I} x  1  x M)
            xm  = num2cell(swapdim(x,2,4),  [1,4,6:ndims(x)]); % ({T x N} x M x {F x ...)
            cinvm = num2cell(cinv, [1:4]); % ({I x N} x M)
            for s = S:-1:1
                am = num2cell(apod{s}, [1:4]);
                asm(s,:) = repmat(am, double([ones(1,4),M]) ./ [ones(1,4), size(am,5)]); % broadcast to ({I x N} x M)
            end
            parfor m = 1:M
                % time delay (I x N x 1)
                if isscalar(cinvm), cinv_ = cinvm{1}; else, cinv_ = cinvm{m}; end % ([1|I] x [1|N] x 1)
                tau = cinv_ .* (dvm{m} + dr) - t0(m); 

                % extract apodization
                % an = repmat(am{m}, double([ones(1,3), N]) ./ [ones(1,3), size(am{m},4)]);
                as = asm(:,m); % init
                a = as{end}; for s = 1:S-1, a = a .* as{s}; end  % accumulate
                a = repmat(a, [ones(1,3), double(N ./ size(a,4)), 1]); % broadcast

                % sample and output (I x N x 1)
                ym = (cellfun(...
                    @(x, tau, a) ...
                    a .* reshape(interp1(x,  1 + tau * fs, interp_type, 0),osz), ...
                    num2cell(xm{m},[1,6:ndims(xm{m})]), pck(tau), pck(a), per_cell{:})) ...
                    ; %#ok<PFBNS>
                ym = reshape(cat(4,ym{:}), [Isz, size(ym,4:ndims(ym)), fsz]);
                y = y + ym;
            end
            y = reshape(y, [Isz, N, 1, fsz]);
        case 'BF'
            % time delay ([I] x N x M)
            tau = cinv .* (dv + dr) - t0;

            % set size of x
            xmn = permute(x, [1,4,5,2,3,6:ndims(x)]); % (T x 1 x 1 x N x M x {F x ...})
            
            % sample and output ([I] x N x M x F x ...)
            y = cellfun(... 
                @(x, tau, a) ...
                reshape(interp1(x,  1 + tau * fs, interp_type, 0), osz), ...
                pck(xmn), pck(tau), per_cell{:} ...
                );
            y = cat(4, y{:}); % unpack
            y = reshape(y, [Isz N M fsz]); % resize
            for s = 1:S, y = y .* apod{s}; end % apply weights
    end
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
            case {'DAS', 'SYN', 'BF', 'MUL', 'delays'}
            otherwise, error('Invalid beamformer.');
        end
        
        % check shape
        Pi = modSize(Pi);
        Pv = modSize(Pv);
        Nv = modSize(Nv);
        Pr = modSize(Pr);
        
        % TODO: avoid handling ambiguous inputs /reorgnanize this code
        % ensure coordinates in the 1st dimension
        function x = modSize(x)
            if size(x,1) <= 4
                % seems legit
            elseif size(x,2) <= 4 && (size(x,1) == 1 || size(x,1) > 4)
                x = x.';
            else % ambiguous ... we'll see what happens
                warning('Input data size is ambiguous.');
            end
        end
    end

    function expand_inputs()
        
        % data sizing
        if strcmp(fun, 'delays')
            [T, N, M] = deal(0);
            [Mv] = size(Pv,2);
            [Mnv] = size(Nv,2);
            [Nr] = size(Pr,2);
            Isz = size(Pi, 2:4);
            I = prod(Isz);

            % expand vectors
            M = max(Mv, Mnv);
            N = Nr;
            if Mv  == 1, Pv   = repmat(Pv,1,M);     Mv  = M; end
            if Mnv == 1, Nv   = repmat(Nv,1,M);     Mnv = M; end
        else
            [T, N, M] = size(x, 1:3);
            [Mv] = size(Pv,2);
            [Mnv] = size(Nv,2);
            [Nr] = size(Pr,2);
            Isz = size(Pi, 2:4);
            Iasz = cell2mat(cellfun(@(a){size(a, 1:5)}, apod(:)));
            Icsz = size(cinv, 1:5);
            I = prod(Isz);
            S  = numel(apod);

            % expand vectors
            if Mv  == 1,   Pv = repmat(Pv,1,M);     Mv  = M; end
            if Mnv == 1,   Nv = repmat(Nv,1,M);     Mnv = M; end
            if Nr  == 1,   Pr = repmat(Pr,1,N);     Nr  = N; end

            % check sizing
            assert(all(Iasz(:,1:3) == 1 | Iasz(:,1:3) == Isz,'all'), 'Apodization data size inconsistent with pixel data size');
            assert(all(Iasz(:, 4 ) == 1 | Iasz(:, 4 ) == N), 'Apodization data size inconsistent with receiver data size');
            assert(all(Iasz(:, 5 ) == 1 | Iasz(:, 5 ) == M), 'Apodization data size inconsistent with transmit data size');
            assert(all(Icsz(1:3) == 1 | Icsz(1:3) == Isz), 'Sound speed data size inconsistent with pixel data size');
            assert(Icsz(4) == 1 || Icsz(4) == N, 'Sound speed data size inconsistent with receiver data size');
            assert(Icsz(5) == 1 || Icsz(5) == M, 'Sound speed data size inconsistent with transmit data size');

        end
        
        % check sizing
        assert(Mv == Mnv || (M == Mv && M == Mnv), 'Inconsistent transmitter data size.');
        assert(N == Nr, 'Inconsistent receiver data size.');
                
        % expand to 3D: 1D -> x, 2D -> x,z, 4D -> x/w, y/w, z/w
        Pi = modDim(Pi);
        Pv = modDim(Pv);
        Nv = modDim(Nv);
        Pr = modDim(Pr);
        
        % move to 3D
        function x = modDim(x)
            xsz = size(x); % size of new x
            xsz(1) = 1; % set first dim to 1
            if size(x, 1) == 1 % assume data is in x
                x = cat(1, x, repmat(cast(zeros(xsz), 'like', x), [2,1]));
            elseif size(x,1) == 2 % assume data is in (x,z)
                x = cat(1, sub(x,1,1), cast(zeros(xsz),'like', x), sub(x,2,1));
            elseif size(x,1) == 3 % assume data is in (x,y,z)
            elseif size(x,1) == 4 % assume data is in (x,y,z,w)
                x = sub(x,1:3,1) ./ sub(x,4,1); % for projective coordinates, project
            else
                error('Improper coordinate dimension.');
            end
        end
    end
end

%#ok<*ASGSL>
%#ok<*NBRAK>



