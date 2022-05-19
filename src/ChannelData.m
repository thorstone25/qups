% CHANNELDATA - Store and process channel data
%
% The ChannelData class stores an N-dimensional datacube and it's temporal 
% axes and provides overloaded methods for manipulating and plotting the 
% data. The ChannelData must have a t0 value consistent with the
% definition of t0 in QUPS to be used with other beamforming algorithms in
% QUPS. Most methods that affect the time axes, such as zeropad or filter, 
% will shift the time axes accordingly.
%
% The underlying datacube can be N-dimensional as long as the first
% dimension is time. The second and third dimensions should be receivers 
% and transmits respectively to be compatible with QUPS. All data must
% share the same sampling frequency fs, but the start time t0 may vary
% across any dimension(s) except for the first and second dimensions. For
% example, if each transmit has a different t0, this can be represented by
% an array of size [1,1,M].
%
% The underlying numeric type of the data can be cast by appending 'T' to
% the type (e.g. singleT(chd) produces a ChannelData) whereas the datacube 
% itself can be cast using the numeric type constructor (e.g. single(chd) 
% produces an array). Overloaded methods are also provided to cast to a 
% gpuArray or tall type. The time axes is also cast to the corresponding
% type. This enables MATLABian casting rules to apply to the object, which
% can be used by other functions.
% 
% See also SEQUENCE TRANSDUCER

classdef ChannelData < matlab.mixin.Copyable

    properties
        data    % channel data (T x N x M x F x ...)
        t0 = 0  % start time (1 x 1 x [1|M] x [1|F] x ...)
        fs = 1  % sampling freuqency (scalar)
    end
    properties(SetAccess=protected,GetAccess=public)
        ord = 'TNM'; % data order: T: time, N: receive, M: transmit
    end
    properties (Dependent)
        time    % time axis (T x 1 x [1|M] x [1|F] x ...)
    end
    properties(Hidden, Dependent)
        T       % number of time samples
        N       % number of receiver channels
        M       % number of transmits
        rxs     % receives vector
        txs     % transmits vector
    end

    % constructor/destructor
    methods
        function self = ChannelData(varargin)
            % CHANNELDATA/CHANNELDATA - Construct a channel data object
            %
            % ChannelData(Name1, Value1, ...) constructs a channel data
            % object via name/value pairs.
            %
            % 

            % set each property by name-value pairs
            for i = 1:2:nargin, self.(lower(varargin{i})) = varargin{i+1}; end
        end
    end
    
    % copyable overloads
    methods(Access=protected)
        function chd = copyElement(self)
            chd = ChannelData('data',self.data,'fs', self.fs,'t0',self.t0);
        end
    end

    % conversion functions
    methods
        function channel_data = getUSTBChannelData(self, sequence, xdc)
            % GETUSTBCHANNELDATA - Create a USTB channel data object
            % 
            % channel_data = getUSTBChannelData(self, sequence, xdc) 
            % creates a USTB compatible channel data object from the QUPS 
            % channel data. USTB must be on the path.
            %
            % 
            
            channel_data = uff.channel_data(...
                'sampling_frequency', self.fs, ...
                'sound_speed', sequence.c0, ...
                'sequence', sequence.getUSTBSequence(xdc, self.t0), ...
                'probe', xdc.getUSTBProbe(), ...
                'data', self.data(:,:,:,:) ... limit to 4 dimensions
                );

        end
    end

    % helper functions
    methods(Hidden)
        function chd = applyFun2Props(chd, fun), 
            chd = copy(chd);
            [chd.t0, chd.fs, chd.data] = deal(fun(chd.t0), fun(chd.fs), fun(chd.data));
        end
        function chd = applyFun2Data(chd, fun), chd = copy(chd); chd.data = fun(chd.data); end
    end

    % data type overloads
    methods
        function chd = gather(chd)  , chd = applyFun2Props(chd, @gather); end
        % gather the underlying data
        function chd = gpuArray(chd), chd = applyFun2Props(chd, @gpuArray); end
        % cast underlying type to gpuArray
        function chd = tall(chd)    , chd = applyFun2Data (chd, @tall); end
        % cast underlying type to tall
        function chd = sparse(chd)  , chd = applyFun2Data (chd, @sparse); end
        % cast underlying type to sparse
        function chd = doubleT(chd) , chd = applyFun2Props(chd, @double); end
        % cast underlying type to double
        function chd = singleT(chd) , chd = applyFun2Props(chd, @single); end
        % cast underlying type to single
        function chd =   halfT(chd) , chd = applyFun2Data (chd, @half); end
        % cast underlying type to half
        function chd =  int64T(chd) , chd = applyFun2Data (chd, @int64); end
        % cast underlying type to int64
        function chd = uint64T(chd) , chd = applyFun2Data (chd, @uint64); end
        % cast underlying type to uint64
        function chd =  int32T(chd) , chd = applyFun2Data (chd, @int32); end
        % cast underlying type to int32
        function chd = uint32T(chd) , chd = applyFun2Data (chd, @uint32); end
        % cast underlying type to uint32
        function chd =  int16T(chd) , chd = applyFun2Data (chd, @int16); end
        % cast underlying type to int16
        function chd = uint16T(chd) , chd = applyFun2Data (chd, @uint16); end
        % cast underlying type to uint16
        function chd =   int8T(chd) , chd = applyFun2Data (chd, @int8); end
        % cast underlying type to int8
        function chd =  uint8T(chd) , chd = applyFun2Data (chd, @uint8); end
        % cast underlying type to uint8
        function T = classUnderlying(self), try T = classUnderlying(self.data); catch, T = class(self.data); end, end % revert to class if undefined
        % underlying class of the data or class of the data
        function T = underlyingType(self), try T = underlyingType(self.data); catch, T = class(self.data); end, end % R2020b+ overload
        % underlying type of the data or class of the data
        function tf = isreal(self), tf = isreal(self.data); end
    end
    
    % implicit casting: functions that request a numeric type may call
    % these functions
    methods
        function x = double(chd), x = double(chd.data); end
        % convert to a double array
        function x = single(chd), x = single(chd.data); end
        % convert to a single array
        function x =   half(chd), x =   half(chd.data); end
        % convert to a half array
        function x =  int64(chd), x =  int64(chd.data); end
        % convert to a int64 array
        function x = uint64(chd), x = uint64(chd.data); end
        % convert to a uint64 array
        function x =  int32(chd), x =  int32(chd.data); end
        % convert to a int32 array
        function x = uint32(chd), x = uint32(chd.data); end
        % convert to a uint32 array
        function x =  int16(chd), x =  int16(chd.data); end
        % convert to a int16 array
        function x = uint16(chd), x = uint16(chd.data); end
        % convert to a uint16 array
        function x =   int8(chd), x =   int8(chd.data); end
        % convert to a int8 array
        function x =  uint8(chd), x =  uint8(chd.data); end
        % convert to a uint8 array
    end

    % DSP overloads 
    methods
        function D = getPassbandFilter(chd, bw, N)
            % GETPASSBANDFILTER Get a passband filter
            %
            % D = GETPASSBANDFILTER(chd, bw) creates a FIR bandpass 
            % digitalFilter object D with a passband between bw(1) and 
            % bw(end). It can be used to filter the ChannelData.
            %
            % D = GETPASSBANDFILTER(chd, bw, N) uses N coefficients.
            %
            % See also DESIGNFILT DIGITALFILTER CHANNELDATA/FILTER
            % CHANNELDATA/FILTFILT CHANNELDATA/FFTFILT


            % defaults
            if nargin < 3, N = 25; end

            % make a
            D = designfilt('bandpassfir', ...
                'SampleRate',chd.fs, ...
                'FilterOrder', N, ...
                'CutoffFrequency1', bw(1), ...
                'CutoffFrequency2', bw(end), ...
                'DesignMethod', 'window' ...
                );
        end
        function chd = filter(chd, D)
            % FILTER Filter data with a digitalFilter
            %
            % chd = FILTER(chd, D) filters the channel data with the
            % digitalFilter D. Use DESIGNFILT to design a digital filter
            %
            % See also DESIGNFILT DIGITALFILTER FILTER

            % data must be in dim 1! but we don't check for this, just
            % require it

            % hard error if we aren't given a digitalFilter
            assert(isa(D, 'digitalFilter'), "Expected a 'digitalFilter' but got a " + class(D) + " instead.");

            % copy semantics
            chd = copy(chd);

            % filter: always applied in dim 1
            chd.data = filter(D, chd.data);

            % adjust time axes
            chd.t0 = chd.t0 - (D.FilterOrder-1)/2/chd.fs;
        end
        function chd = filtfilt(chd, D)
            % FILTFILT Filter data with a digitalFilter
            %
            % chd = FILTFILT(chd, D) filters the channel data with the
            % digitalFilter D. Use DESIGNFILT to design a digital filter
            %
            % See also DESIGNFILT DIGITALFILTER FILTFILT

            % data must be in dim 1! but we don't check for this, just
            % require it

            % hard error if we aren't given a digitalFilter
            assert(isa(D, 'digitalFilter'), "Expected a 'digitalFilter' but got a " + class(D) + " instead.");

            % copy semantics
            chd = copy(chd);

            % filter: always applied in dim 1
            chd.data = filtfilt(D, chd.data);
        end
        function chd = fftfilt(chd, D)
            % FFTFILT Filter data with a digitalFilter
            %
            % chd = FFTFILT(chd, D) filters the channel data with the
            % digitalFilter D. Use DESIGNFILT to design a digital filter
            %
            % See also DESIGNFILT DIGITALFILTER FFTFILT

            % data must be in dim 1! but we don't check for this, just
            % require it

            % hard error if we aren't given a digitalFilter
            assert(isa(D, 'digitalFilter'), "Expected a 'digitalFilter' but got a " + class(D) + " instead.");

            % copy semantics
            chd = copy(chd);

            % filter: always applied in dim 1
            chd.data(:,:) = fftfilt(D, chd.data(:,:));

            % adjust time axes
            chd.t0 = chd.t0 - (D.FilterOrder-1)/2/chd.fs;
        end
        function chd = hilbert(chd, varargin)
            % HILBERT - overloads the hilbert function
            %
            % chd = hilbert(chd) applies the hilbert function to the data.
            %
            % chd = hilbert(chd, ...) forwards arguments to MATLAB's 
            % hilbert function
            %
            % See also HILBERT
            chd = copy(chd);
            chd.data = hilbert(chd.data, varargin{:});
        end
        function chd = fft(chd, N, dim)
            % FFT - overload of fft
            %
            % chd = FFT(chd) computes the fft of the channel data along the 
            % time axis. The time axes is unchanged.
            %
            % chd = FFT(chd, N) computes the N-point fft.
            %
            % chd = FFT(chd, N, dim) or FFT(chd, [], dim) operates along
            % dimension dim.
            %
            % See also FFT CHANNELDATA/FFTSHIFT

            % defaults
            if nargin < 3, dim = 1; end
            if nargin < 2, N = []; end
            chd = copy(chd);
            chd.data = fft(chd.data, N, dim); % take the fourier transform
        end
        function chd = fftshift(chd, dim)
            % FFTSHIFT - overload of fftshift
            %
            % chd = FFTSHIFT(chd) swaps the left and right halves of the 
            % data along the time dimension. The time and frequency axes  
            % are unchanged.
            %
            % chd = FFTSHIFT(chd, dim) operates along dimension dim.
            %
            % See also FFTSHIFT CHANNELDATA/FFT 

            if nargin < 2, dim = 1; end
            chd = copy(chd);
            chd.data = fftshift(chd.data, dim);
        end
        function chd = ifft(chd, N, dim)
            % IFFT - overload of fft
            %
            % chd = IFFT(chd) computes the inverse fft of the channel data 
            % along the time axis. The time and frequency axes are 
            % unchanged.
            %
            % chd = IFFT(chd, N) computes the N-point inverse fft.
            %
            % chd = IFFT(chd, N, dim) or IFFT(chd, [], dim) operates along
            % dimension dim.
            %
            % See also IFFT CHANNELDATA/IFFTSHIFT


            % defaults
            if nargin < 3, dim = 1; end
            if nargin < 2, N = []; end
            chd = copy(chd);
            chd.data = ifft(chd.data, N, dim); % take the fourier transform
        end
        function chd = ifftshift(chd, dim)
            % IFFTSHIFT - overload of fftshift
            %
            % chd = IFFTSHIFT(chd) swaps the left and right halves of the 
            % data along the time dimension. The time and frequency axes  
            % are unchanged.
            %
            % chd = IFFTSHIFT(chd, dim) operates along dimension dim.
            %
            % IFFTSHIFT undoes the effects of fftshift
            % 
            % See also IFFTSHIFT CHANNELDATA/IFFT 
            if nargin < 2, dim = 1; end
            chd = copy(chd);
            chd.data = ifftshift(chd.data, dim);
        end
        function chd = resample(chd, fs, varargin)
            % RESAMPLE - Resample the data in time
            %
            % chd = RESAMPLE(chd, fs) resamples the data at sampling
            % frequency fs. And returns a new ChannelData object.
            %
            % chd = RESAMPLE(chd, fs, ..., METHOD) specifies the method of 
            % interpolation. The default is linear interpolation.  
            % Available methods are:
            %   'linear' - linear interpolation
            %   'pchip'  - shape-preserving piecewise cubic interpolation
            %   'spline' - piecewise cubic spline interpolation
            %
            % chd = RESAMPLE(chd, fs, ..., arg1, arg2, ...) forwards 
            % arguments to MATLAB's RESAMPLE function
            %
            % See also RESAMPLE

            % save original data prototypes
            [Tt, Tf, Td] = deal(chd.t0, chd.fs, zeros([0,0], 'like', chd.data));
            
            % Make new ChannelData (avoid modifying the original)
            chd = copy(chd);

            % ensure numeric args are non-sparse, double
            chd = (doubleT(chd)); % data is type double
            fs = (double(fs)); % new frequency is type double
            inum = cellfun(@isnumeric, varargin); % all numeric inputs are type double
            varargin(inum) = cellfun(@double, varargin(inum), 'UniformOutput', false);

            % resample in time - no support for other dims: fs is required arg
            [y, ty] = resample(chd.data, chd.time, fs, varargin{:}, 'Dimension', 1);
            [chd.fs, chd.t0, chd.data] = deal(fs, ty(1), y);            

            % cast back to original type
            tmp = cellfun(@(x,T) cast(x, 'like', T), {chd.t0, chd.fs, chd.data}, {Tt, Tf, Td}, 'UniformOutput', false);
            [chd.t0, chd.fs, chd.data] = tmp{:};

        end
        function chd = real(chd)    , chd = applyFun2Data (chd, @real); end
        function chd = imag(chd)    , chd = applyFun2Data (chd, @imag); end
        function chd = abs(chd)     , chd = applyFun2Data (chd, @abs); end
        function chd = angle(chd)   , chd = applyFun2Data (chd, @angle); end
        function chd = mag2db(chd)  , chd = applyFun2Data (chd, @mag2db); end
        function chd = mod2db(chd)  , chd = applyFun2Data (chd, @mod2db); end
    end

    % DSP helpers
    methods
        function chd = zeropad(chd, B, A)
            % ZEROPAD - Zero pad the data in time
            %
            % chd = ZEROPAD(chd, B) prepends B zeros to the ChannelData 
            % data in time
            %
            % chd = ZEROPAD(chd, B, A) also appends A zeros to the 
            % ChannelData data in time
            %
            % When using this function, the time axis is adjusted.
            % 
            % See also CIRCSHIFT

            if nargin < 2 || isempty(B), B = 0; end
            if nargin < 3 || isempty(A), A = 0; end
            assert(A >= 0 && B >= 0, 'Data append or prepend size must be positive.');

            chd = copy(chd); % copy semantics
            chd.data(end+(B+A),:) = 0; % append A + B zeros in time to the data
            chd.data = circshift(chd.data, B, 1); % shift B of the zeros to the front
            chd.t0 = chd.t0 - B ./ chd.fs; % shift start time for B of the zeros
        end
    
        function fc = estfc(chd)
            % ESTFC - Estimate the central frequency
            %
            % fc = ESTFC(chd) estimates the central frequency of the
            % ChannelData chd by choosing the mode of the maximum frequency
            % across all channels. This method should be updated with a
            % better heuristic. Alternatively, the central frequency may be
            % defined by the Waveform that was transmitted, the impulse
            % response of the Transducer(s), or the Transducer itself
            %
            % See also TRANSDUCER WAVEFORM SEQUENCE
            % TRANSDUCER/ULTRASOUNDTRANSDUCERIMPULSE


            f = chd.fs * ((0:chd.T-1)' ./ chd.T); % compute frequency axis
            y = fft(chd, [], 1); % compute fft
            fc = f(mode(argmax(abs(y.data), [], 1), 'all')); % select peak over receives/transmits
        end
    
        function chd = rectifyt0(chd, interp, t0_)
            % RECTIFYT0 - Collapse t0 to a scalar
            %
            % chd = RECTIFYT0(chd) returns a ChannelData object with a
            % single value of t0 by resampling all channels onto a single
            % time axis. 
            %
            % chd = RECTIFYT0(chd, interp) specifices the interpolation
            % method as recognized by the sample function.
            %
            % See also CHANNELDATA/SAMPLE 


            if nargin < 2, interp = 'cubic'; end % default interpolation
            if nargin < 3, t0_ = min(chd.t0, [], 'all'); end % get global start time
            toff = chd.t0 - t0_; % get offset across upper dimensions
            npad = ceil(max(toff,[],'all') * chd.fs); % furthest sample containing data
            chd = zeropad(chd,0,npad); % extend time-axes
            tau = chd.time + toff; % get delays to resample all traces
            y = chd.sample(tau, interp); % resample
            chd = ChannelData('data', y, 't0', t0_, 'fs', chd.fs); % make new object

        end
    
        function y = sample(chd, tau, interp)
            % SAMPLE Sample the channel data in time
            %
            % y = SAMPLE(chd, tau) samples the ChannelData chd at the times
            % given by the delays tau. 
            %
            % tau must have broadcastable sizing in the non-temporal
            % dimensions. In other words, in all dimensions d, either of
            % the following must hold
            %   1)   size(tau,d)    ==   size(chd.data,d) 
            %   2)  (size(tau,d) == 1 || size(chd.data,d) == 1)
            %
            % The underlying routines are optimized for compute
            % performance. Consider using a for-loop if memory is a
            % concern.
            % 
            % y = SAMPLE(chd, tau, interp) specifies the interpolation
            % method. Interpolation is handled by the built-in interp1 
            % function. The available methods are:
            %
            %   'linear'   - (default) linear interpolation **
            %   'nearest'  - nearest neighbor interpolation **
            %   'next'     - next neighbor interpolation
            %   'previous' - previous neighbor interpolation
            %   'spline'   - piecewise cubic spline interpolation 
            %   'pchip'    - shape-preserving piecewise cubic interpolation
            %   'cubic'    - cubic convolution interpolation for ***
            %                uniformly-spaced data **
            %   'v5cubic'  - same as 'cubic' ***
            %   'makima'   - modified Akima cubic interpolation
            %   'freq'     - frequency domain sinc interpolation ****
            %   'lanczos3' - lanczos kernel with a = 3 **
            % 
            %    **   GPU support is enabled via interpd
            %    ***  GPU support is enabled via interp1
            %    **** GPU support is native
            %
            % See also INTERP1 INTERPD CHANNELDATA/RECTIFYT0

            % defaults
            if nargin < 3, interp = 'linear'; end

            % check condition that we can implicitly broadcast
            for d = 2:ndims(tau)
                assert(any(size(tau,d) == size(chd.data,d)) || ...
                      any([size(tau,d)  , size(chd.data,d)] == 1), ...
                    'Delay size must match the data size (%i) or be singleton in dimension %i.',...
                    size(chd.data, d), d ...
                    );
            end

            % dispatch
            assert(chd.tdim == 1, 'Time must be in the first dimension of the data.'); % TODO: generalize this restriction if necessary
            if (isa(chd.data, 'gpuArray') ... 
                    && ismember(interp, ["nearest", "linear", "cubic", "lanczos3"])) ... 
                    && logical(exist('interpd.ptx', 'file'))
                    % interpolate on the gpu via ptx if we can
                    y = interpd(chd.data, (tau - chd.t0) * chd.fs, 1, interp, 0); 

            elseif interp == "freq"
                % extend data if necessary
                ntau = (tau - chd.t0) * chd.fs; % sample delays (I x [1|N] x [1|M] x [1|F] x ...)
                nwrap = min(0,floor(min(ntau,[],'all'))); % amount the signal is wrapped around from negative time
                next  = max(0, ceil(max(ntau,[],'all')) - chd.T-1); % amount the signal is extended in positive time
                chd = zeropad(chd, -nwrap, next); % extend signal to ensure we sample zeros

                x = chd.data; % reference data (T x N x M x [1|F] x ...)
                L = chd.T; % fft length
                % l = (0:L-1)'; % make new time vector in sample dimension
                d = max(ndims(x), ndims(ntau)); % find max dimension
                ntau = swapdim(ntau, d+1, 1); % move sampling to next one, a free dimension
                
                % apply phase shifts and sum (code written serially to 
                % request in-place operation from MATLAB)
                x = fft(x, L, 1); % put data in freq domain (L x N x M x [1|F] x ... x I)
                wL = exp(2i*pi*ntau./L); % sampling steering vector (1 x [1|N] x [1|M] x [1|F] x ... x I)
                xl = num2cell(x, 2:ndims(x)); % splice data in freq domain (L x {N x M x [1|F] x ... x I})
                y = 0; % initialize accumulator
                if isa(x, 'gpuArray') || isa(wL, 'gpuArray'), clu = 0; % avoid parfor on gpuArray
                else, clu = gcp('nocreate'); if isempty(clu), clu = 0; end % otherwise, use current parpool   
                end
                parfor (l = (1:L), clu), y = y + wL.^(l-1) .* xl{l}; end % apply phase shift and sum over freq (1 x N x M x [1|F] x ... x I)
                y = swapdim(y, 1, d+1); % move samples back to first dim (I x N x 1 x M' x F x ...)
                                    
            else % use interp1, iterating over matched dimensions, broadcasting otherwise
                    % convert to index based coordinates
                    ntau = (tau - chd.t0) * chd.fs; % get full size sample delays
                    x = chd.data; % reference data

                    % get dims to pack versus loop
                    mxdim = max(ndims(x), ndims(tau)); % maximum number of dimensions
                    xsing = 1+find(size(x   ,2:mxdim) == 1); % x    singular
                    nsing = 1+find(size(ntau,2:mxdim) == 1); % ntau singular
                    bdim  = setxor(xsing, nsing); % broadcast dimensions
                    pdim  = union(1,bdim); % pack dim1 & all broadcast dimensions
                    ldim  = setxor(1:mxdim,pdim); % loop over matching dims (the compliment)
                    Dn    = max(pdim(size(ntau,pdim) ~= 1)); if isempty(Dn), Dn = 1; end % dimensions
                    Dx    = max(pdim(size(x   ,pdim) ~= 1)); if isempty(Dx), Dx = 1; end % dimensions


                    % sample: output is [size(ntau,1), size(ntau,2:Dn), size(x,2:Dx)].
                    xc = num2cell(x   ,pdim); % pack/splice data
                    nc = num2cell(ntau,pdim); % pack/splice data
                    if isa(x, 'gpuArray') || isa(ntau, 'gpuArray'), clu = 0; % avoid parfor on gpuArray
                    elseif numel(xc) == 1, clu = 0; % don't parallelize over 1 thing
                    else, clu = gcp('nocreate'); if isempty(clu), clu = 0; end % otherwise, use current parpool
                    end
                    parfor (i = 1:numel(xc), clu), y{i} = interp1(xc{i}, nc{i}, interp, 0); end
                    
                    % check my logic ...
                    assert(...
                        isempty(2:Dn) ||  all(cellfun(@(y) all(size(y,2:Dn) == size(ntau,(2:Dn)) | ~ismember(2:Dn, pdim)), y)), ...
                        'Internal error. Please either check your sizing, check your MATLAB version, submit an issue, or give up.' ...
                        );

                    assert(...
                        isempty(2:Dx) || all(cellfun(@(y) all(size(y,Dn-1+(2:Dx)) == size(x,(2:Dx)) | ~ismember(2:Dx, pdim)), y)), ...
                        'Internal error. Please either check your sizing, check your MATLAB version, submit an issue, or give up.' ...
                        );

                    % output is [Tp, size(ntau,2:Dn), size(x,2:Dx)] - we want to identify
                    % the broadcast dimensions of x and pull them back into
                    % their original dimensions
                    yord = [1, 2:Dn, Dn-1+(2:Dx)]; % full dimension size
                    yord(nsing) = Dn+nsing-1; % swap out entries of ntau that are singular
                    yord(Dn+nsing-1) = nsing; % corresponding to where x is non-singular
                    if isequal(yord, [1]), yord(2) = 2; end %#ok<NBRAK2> % special case: [1,] -> [1,] not accepted by MATLAB
                    y = cellfun(@(y) {permute(y, yord)}, y); % and permute it down

                    % output is now (Tp x N x M x F x ...) in principal,
                    % but with cell arrays over the looped dimensions
                    % restore output sizing using trailing singleton dimensions
                    y = cat(mxdim+1, y{:}); % (Tp x [1|N] x [1|M] x F x ... x [N|1] x [M|1]
                    if ~isempty(ldim), lsz = size(ntau,ldim); else, lsz = []; end % forward empty
                    y = reshape(y, [size(y,1:mxdim), lsz]); % restore data size in upper dimension
                    y = swapdim(y,ldim,mxdim+(1:numel(ldim))); % fold upper dimensions back into original dimensions
            end
        end
        
        function chd = alignInt(chd, interp)
            % ALIGNINT - Align data to integer sampling
            %
            % chd = ALIGNINT(chd, interp) returns a ChannelData object
            % resampled to have an integer time axis. 
            % 
            % This can be used to compress the data to an integer type, but
            % may erase information about the true sampling frequency and
            % initial time which is used by beamforming algorithms. It can 
            % be useful if you want to store many ChannelData objects with 
            % less precision and only cast them to a floating-point type 
            % when they are being processed.
            % 
            % See also CHANNELDATA/RECTIFYT0 CHANNELDATA/SAMPLE 

            if nargin < 2, interp = 'cubic'; end
            n0_ = floor(min(chd.t0, [], 'all') * chd.fs); % get minimum integer time
            chd = rectifyt0(chd, interp, n0_  / chd.fs); % set data on the same time axis
            [chd.t0, chd.fs] = deal(n0_, 1);
        end    
    end

    % plotting and display
    methods
        function h = imagesc(self, m, varargin)
            % IMAGESC - Overload of imagesc function
            %
            % h = IMAGESC(self, m) displays transmit m of the channel data 
            % and returns the handle h.
            %
            % h = IMAGESC(self, m, 'YData', yax) uses the yax for the
            % y-axis instead of the time domain.
            % 
            % h = IMAGESC(self, m, Name, Value, ...) passes the followin
            % arguments to the imagesc function.
            %
            % Example:
            %
            % % Show the data in the frequency domain in MHz
            % figure; 
            % imagesc(fft(chd), 1, 'YData', 1e-6*chd.fs*(0:chd.T-1)/chd.T)
            %
            % See also IMAGESC

            % parse inputs
            if nargin < 2, m = floor((self.M+1)/2); end
            if nargin >= 3 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                ax = varargin{1}; varargin(1) = [];
            else
                ax = gca;
            end

            % get full data sizing
            dims = max([ndims(self.data), ndims(self.time), 3]);
            dsz = size(self.data, 3:dims);
            tsz = size(self.time, 3:dims);

            % we index the data linearly, but the time axes may be
            % implicitly broadcasted: we can use ind2sub to recover it's
            % sizing
            i = cell([dims-3+1, 1]); 
            [i{:}] = ind2sub(dsz, m); % get cell array of indices
            i = min([i{:}], tsz); % restrict to size of self.time

            % select the transmit/frame - we can only show first 2
            % dimensions
            d = self.data(:,:,m); 

            % choose to show real part or dB magnitude
            if isnumeric(d), d = double(d); end % cast integer types for abs, half types for imagesc
            if ~isreal(d), d = mod2db(d); end

            % get the time axes for this channel
            t = double(sub(self.time, num2cell(i), [3:dims])); %#ok<NBRAK> 

            % choose which dimensions to show
            axes_args = {'XData', 1:self.N, 'YData', t};

            % show the data
            h = imagesc(ax, d, axes_args{:}, varargin{:});
        end

        function h = animate(self, varargin)
            % ANIMATE Show the data across transmits
            %
            % h = ANIMATE(self, ...) iteratively calls imagesc to quickly
            % display the data. All trailing arguments are passed to
            % ChannelData/imagesc.
            %
            % See also CHANNELDATA/IMAGESC IMAGESC

            if nargin >= 2 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                ax = varargin{1}; varargin(1) = [];
            else
                ax = gca;
            end

            % now use the handle only
            for f = 1:prod(size(self.data,3:ndims(self.data)))
                if ~isvalid(ax), break; end % quit if handle destroyed
                h = imagesc(self, f, ax, varargin{:});
                drawnow limitrate; pause(1/20);
            end
        end
    end

    % dependent methods
    methods
        function set.time(self, t)
            % set.time - set the time axis
            %
            % time must be of dimensions (T x 1 x [1|M]) where M is the
            % number of transmits.

            % TODO: validate size
            % TODO: warn or error if time axis is not regularly spaced

            %  get the possible sampling freuqencies
            fs_ = unique(diff(t,1,1));

            % choose the sampling frequency with the best adherence to the
            % data
            t0_ = sub(t, 1, 1);
            m = arrayfun(@(fs) sum(abs((t0_ + (0 : numel(t) - 1)' ./ fs ) - t), 'all'), fs_);
            fs_ = fs_(argmin(m));

            % set the time axis
            self.t0 = t0_;
            self.fs = fs_;
        end
        function t = get.time(self), t = cast(self.t0 + shiftdim((0 : gather(self.T) - 1)', self.tdim-1) ./ self.fs, 'like', real(self.data)); end % match data type, except always real
        function T = get.T(self), T = size(self.data,self.tdim); end
        function N = get.N(self), N = size(self.data,self.ndim); end
        function M = get.M(self), M = size(self.data,self.mdim); end
        function n = get.rxs(self), n=cast(shiftdim(1:self.N,0 ), 'like', real(self.data)); end
        function m = get.txs(self), m=cast(shiftdim(1:self.M,-1), 'like', real(self.data)); end
    end

    % convenience functions (used internally for now)
    properties(Hidden, Dependent)
        tdim
        ndim
        mdim
    end
    methods
        function d = get.tdim(self), d = find(self.ord == 'T'); end
        function d = get.ndim(self), d = find(self.ord == 'N'); end
        function d = get.mdim(self), d = find(self.ord == 'M'); end
    end
end

% TODO: make a default interpolation method property