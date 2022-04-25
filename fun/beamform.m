function y = beamform(fun, Pi, Pr, Pv, Nv, x, t0, fs, c, varargin)

%
% function y = beamform(fun, Pi, Pr, Pv, Nv, x, t0, fs, c, varargin)
%
% Beamform the data at the given pixels. Optionally sum across apertures.
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
%  fun -         algorithm for aperture summation {'DAS'|'SYN'|'BF'|'delays'}
%  Pi -          pixel positions (3 x I[1 x I2 x I3])
%  Pr -          receiver positions (3 x N)
%  Pv -          (virtual) transmitter positions (3 x M)
%  Nv -          (virtual) transmitter normal (3 x M)
%  x -           datacube of complex sample values (T x N x M)
%  t0 -          initial time for the data
%  fs -          sampling frequency of the data
%  c -           sound speed used for beamforming
% 
% Outputs:
%   y -          complex pixel values per transmit/channel (I[ x N x M])
%                I -> pixels, M -> transmitters, N -> receivers, T -> time samples
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
if any(cellfun(@(c) isType(c, 'single'), {Pi, Pr, Pv, Nv}))
    posType = 'single';
else
    posType = 'double';
end
if gpuDeviceCount && any(cellfun(@(v)isa(v, 'gpuArray'), {Pi, Pr, Pv, Nv, x}))
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

% get speed of sound inverse
cinv = 1./c;

% parse inputs
parse_inputs();

% store this option input for convenience
per_cell = {'UniformOutput', false};
    
if device

    % match position precision to inputData precision
    if ~strcmp(idataType, posType) && ~strcmp(fun, 'delays'), 
        warning('Position data type will match the input data type.');
        idataType = posType; 
    end

    % warn if non-linear interp was requested
    if ~strcmp(interp_type, 'linear')
        warning('Only linear interpolation supported on the GPU. Set ''device'' to 0 to use an alternate interpolation method.');
    end

    % load the kernel
    src_folder = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    
    postfix = '';
    switch idataType
        case 'single'
            postfix = [postfix 'f'];
        case 'double'
            postfix = postfix;
    end
    
    % currently selected gpu device handle
    g = gpuDevice();
    
    % reselect gpu device and copy over inputs if necessary
    if device > 0 && g.Index ~= device
        inps = {Pi, Pr, Pv, Nv, x, t0, fs, cinv};
        oclass = cellfun(@class, inps, per_cell{:}); 
        tmp = cellfun(@gather, inps, per_cell{:});
        g = gpuDevice(device); 
        tmp = cellfun(@(inp, cls)cast(inp, cls), tmp, oclass, per_cell{:});
        [Pi, Pr, Pv, Nv, x, t0, fs, cinv] = deal(tmp{:});
    end
    
    % If bf.ptx isn't there, ask the user to generate it
    if ~exist('bf.ptx', 'file')
        error('The source code is not compiled for MATLAB on this system. Try using UltrasoundSystem.recompileCUDA() to compile it.')
    end
    
    k = parallel.gpu.CUDAKernel(...
        'bf.ptx',...
        fullfile(src_folder, 'bf.cu'),...
        [fun, postfix]);
    
    % send constant data to GPU
    dprototype = complex(gpuArray(zeros(0, idataType)));
    pprototype = gpuArray(zeros(0, idataType));
    
    x = cast(x, 'like', dprototype);
    [tmp] = cellfun(@(p)cast(p, 'like', pprototype), {Pi, Pr, Pv, Nv}, per_cell{:});
    [Pi, Pr, Pv, Nv,] = deal(tmp{:});
    
    % expand all inputs to 3D
    expand_inputs();
    
    % data sizing
    [T, N, M] = size(x);
    I = numel(Pi) / 3; % guaranteed 3D
    Isz = size(Pi, 2:4); % I1 x I2 x I3 == I
    
    % get kernel and frame sizing
    switch fun
        case 'DAS'
            osize = {1}; % delay and sum (tx/rx)
        case 'SYN'
            osize = {N}; % beamform and sum the transmit aperture only
        case {'BF', 'delays'}
            osize = {N, M}; % beamform only | times only
    end
    
    % get output data type
    switch fun
        case {'DAS', 'SYN', 'BF'}
            obufprototype = complex(dprototype);
        case {'delays'}
            obufprototype = real(dprototype);
    end
   
    % kernel sizes
    nThreads = k.MaxThreadsPerBlock; % threads per block
    nBlocks = min([...
        g.MaxThreadBlockSize(1), ... max cause device reqs
        ceil(I / nThreads)... max cause number of pixels
        ceil(g.AvailableMemory / (2^8) / (prod(osize{:}) * nThreads)),... max cause GPU memory reqs (empirical)
        ]); % blocks per frame
    kI = nThreads*nBlocks; % pixels per frame
    nF = ceil(I ./ kI); % number of frames
    
    % constant arg type casting
    tmp = cellfun(@uint64, {T,M,N,kI,nF}, per_cell{:});
    [T, M, N, kI, nF] = deal(tmp{:});
    
    % set constant args
    k.setConstantMemory('T', T, 'M', M, 'N', N, 'I', kI, 'VS', VS);
    
    % set kernel size
    k.ThreadBlockSize = nThreads;
    k.GridSize = nBlocks;
    
    % allocate output data buffer
    osize = cat(2, {kI}, osize);
    osize = cellfun(@uint64, osize, per_cell{:});
    yg = zeros([osize{:}], 'like', obufprototype);
    
    % partition input pixels per frame
    Pif = NaN([3, kI, nF], 'like', Pi); % initialize value
    Pif(1:3, 1:I) = Pi(:,:); % place valid pixel positions
    Pif = num2cell(Pif, [1,2]); % pack in cells per frame
    
    % for each frame, run the kernel
    switch fun
        case {'DAS','SYN','BF'}
            y = cellfun(@(pi) cast(...
                k.feval(yg, pi, Pr, Pv, Nv, x, t0, fs, cinv),...
                'like', odataPrototype), Pif, per_cell{:});
        case {'delays'}
            y = cellfun(@(pi) cast(...
                k.feval(yg, pi, Pr, Pv, Nv, cinv), ...
                'like', odataPrototype), Pif, per_cell{:});
    end
    
    % reshape output and truncate garbage
    y = cat(1, y{:}); % I' [x N [x M]]
    y = y(1:I,:,:); % trim the junk
    y = reshape(y, [Isz, size(y,2), size(y,3)]); % I1 x I2 x I3[ x N [x M]]
    
else
    
    % cast constant data on CPU
    dprototype = complex(zeros(0, idataType));
    pprototype = zeros(0, posType);
    
    [x, t0, fs] = dealfun(@(x) cast(x, 'like', dprototype), x, t0, fs);
    [tmp] = cellfun(@(p)cast(p, 'like', pprototype), ... 
        {Pi, Pr, Pv, Nv}, per_cell{:});
    [Pi, Pr, Pv, Nv] = deal(tmp{:});
    
    % expand all inputs to 3D
    expand_inputs();
    
    % data sizing
    [T, N, M] = size(x);
    I = numel(Pi) / 3;
    Isz = size(Pi,2:4);
    
    % cast to integers
    tmp = cellfun(@uint64, {T,N,M,I}, per_cell{:});
    [T, N, M, I] = deal(tmp{:});
    
    % permute to compatible casting dimensions
    Pr = swapdim(Pr, 2, 5); % 3 x 1 x 1 x 1 x N
    Pv = swapdim(Pv, 2, 6); % 3 x 1 x 1 x 1 x 1 x M
    Nv = swapdim(Nv, 2, 6); % 3 x 1 x 1 x 1 x 1 x M
    
    % transmit sensing vector
    rv = Pi - Pv; % 3 x I1 x I2 x I3 x 1 x M
    if VS % virtual source delays
        dv = vecnorm(rv, 2, 1) .* sign(sum(rv .* Nv,1));
    else % plane-wave delays
        dv = sum(rv .* Nv, 1);
    end % 1 x I1 x I2 x I3 x 1 x M
    
    % receive sensing vector
    dr = vecnorm(Pi - Pr, 2, 1); % 1 x I1 x I2 x I3 x N x 1
    
    % bring to I1 x I2 x I3 x N x M == I x N x M
    dv = shiftdim(dv, 1); 
    dr = shiftdim(dr, 1);
    
    % time vector
    t = t0 + cast(colon(0,T-1).', 'like', t0) ./ fs;
    
    % temporal packaging function
    pck = @(x) num2cell(x, [1:3]); % pack for I
    
    switch fun
        case 'delays'
            y = cast(cinv .* (dv + dr), 'like', odataPrototype);
            
        case 'DAS'
            y = zeros([Isz, 1, 1], 'like', odataPrototype);
            dvm = num2cell(dv, [1:3]); % ({I} x 1 x M)
            xmn = num2cell(x,  [1]); % ({T} x N x M)
            parfor m = 1:M
                yn = zeros([Isz, 1, 1], 'like', odataPrototype);
                drn = num2cell(dr, [1]); % ({I} x N x 1)
                for n = 1:N
                    % time delay (I x 1 x 1)
                    tau = cinv .* (dvm{m} + drn{n});
                    
                    % sample and output (I x 1 x 1)
                    yn = yn + cast(...
                        interpn(t, xmn{1,n,m}, tau, interp_type, 0), ...
                        'like', odataPrototype);
                end
                y = y + yn;
            end
            
        case 'SYN'
            y = zeros([Isz, N, 1], 'like', odataPrototype);
            dvm = num2cell(dv, [1:3]); % ({I} x  1  x M)
            xm  = num2cell(x,  [1,2]); % ({T} x {N} x M)
            parfor m = 1:M
                % time delay (I x N x 1)
                tau = cinv .* (dvm{m} + dr); 
                
                % sample and output (I x N x 1)
                y = y + cast(cell2mat(cellfun(...
                    @(x, tau) ...
                    interpn(t, x, tau, interp_type, 0), ...
                    num2cell(xm{m},1), pck(tau), per_cell{:})), ...
                    'like', odataPrototype); %#ok<PFBNS>
            end
            
        case 'BF'
            % time delay (I x N x M)
            tau = cinv .* (dv + dr);
            
            % sample and output (I x N x M)
            y = cell2mat(cellfun(... 
                @(x, tau) cast(...
                interpn(t, x, tau, interp_type, 0), ...
                'like', odataPrototype), ...
                num2cell(x,1), pck(tau), per_cell{:}));
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
            case {'DAS', 'SYN', 'BF', 'delays'}
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
            [~, Mv] = size(Pv);
            [~, Mnv] = size(Nv);
            [~, Nr] = size(Pr);
            Isz = size(Pi, 2:ndims(Pi));
            I = prod(Isz);

            % expand vectors
            M = max(Mv, Mnv);
            N = Nr;
            if Mv  == 1, Pv   = repmat(Pv,1,M);     Mv  = M; end
            if Mnv == 1, Nv   = repmat(Nv,1,M);     Mnv = M; end
        else
            [T, N, M] = size(x);
            [~, Mv] = size(Pv);
            [~, Mnv] = size(Nv);
            [~, Nr] = size(Pr);
            Isz = size(Pi, 2:ndims(Pi));
            I = prod(Isz);
            Ia = I;

            % expand vectors
            if Mv  == 1,   Pv = repmat(Pv,1,M);     Mv  = M; end
            if Mnv == 1,   Nv = repmat(Nv,1,M);     Mnv = M; end
            if Nr  == 1,   Pr = repmat(Pr,1,N);     Nr  = N; end
        end
        
        % check sizing
        assert(Mv == Mnv || (M == Mv && M == Mnv), 'Inconsistent transmitter data size.');
        assert(N == Nr, 'Inconsistent receiver data size.');
        assert(Ia == 1 || Ia == I, 'Inconsistent apodization data size');
        
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
                x = cat(1, x, repmat(zeros(xsz, 'like', x), [2,1]));
            elseif size(x,1) == 2 % assume data is in (x,z)
                x = cat(1, sub(x,1,1), zeros(xsz,'like', x), sub(x,2,1));
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



