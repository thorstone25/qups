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
    properties(Access=public)
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
            % CHANNELDATA - Construct a ChannelData object
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
            chd = ChannelData('data',self.data,'fs', self.fs,'t0',self.t0,'ord',self.ord);
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
            self = rectifyDims(self); % make sure it's in order 'TNM' first
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
        function chd = applyFun2Dim(chd, fun, dim, varargin),
            chd = copy(chd); % copy semantics
            chd.data = matlab.tall.transform(@dimfun, chd.data, varargin{:}); % apply function in dim 1; % set output data

            % dim1 mapping function: dim d gets sent to dim 1 and back.
            function x = dimfun(x, varargin)
                x = swapdim(x, 1, dim); % send dim d to dim 1
                x = fun(x, varargin{:}); % operate in dim 1
                x = swapdim(x, 1, dim); % send dim d back
            end
        end
    end

    % data type overloads
    methods
        function chd = gather(chd)  , chd = applyFun2Data(chd, @gather); end
        % gather the underlying data
        function chd = gpuArray(chd), chd = applyFun2Data (chd, @gpuArray); end
        % cast underlying type to gpuArray
        function chd = tall(chd)    , chd = applyFun2Data (chd, @tall); end
        % cast underlying type to tall
        function chd = sparse(chd)  , chd = applyFun2Data (chd, @sparse); end
        % cast underlying type to sparse
        function chd = doubleT(chd) , chd = applyFun2Data (chd, @double); end
        % cast underlying type to double
        function chd = singleT(chd) , chd = applyFun2Data (chd, @single); end
        % cast underlying type to single
        function chd =   halfT(chd) , chd = applyFun2Data (chd, @halfT); end
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
        % whether the underlying data is real
        function tf = istall(self), tf = istall(self.data); end
        % whether the underlying data is tall
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
        function chd = filter(chd, D, dim)
            % FILTER Filter data with a digitalFilter
            %
            % chd = FILTER(chd, D) filters the channel data with the
            % digitalFilter D. Use DESIGNFILT to design a digital filter.
            %
            % chd = FILTER(chd, D, dim) applies the filter in dimension
            % dim. The default is the time dimension.
            %
            % See also DESIGNFILT DIGITALFILTER FILTER

            % hard error if we aren't given a digitalFilter
            assert(isa(D, 'digitalFilter'), "Expected a 'digitalFilter' but got a " + class(D) + " instead.");

            % defaults
            if nargin < 3, dim = chd.tdim; end

            % filter: always applied in dim 1
            chd = applyFun2Dim(chd, @(x) filter(D, x), dim);

            % adjust time axes
            if dim == chd.tdim
                chd.t0 = chd.t0 - (D.FilterOrder-1)/2/chd.fs;
            end
        end
        function chd = filtfilt(chd, D, dim)
            % FILTFILT Filter data with a digitalFilter
            %
            % chd = FILTFILT(chd, D) filters the channel data with the
            % digitalFilter D. Use DESIGNFILT to design a digital filter
            %
            % chd = FILTFILT(chd, D, dim) applies the filter in dimension
            % dim. The default is the time dimension.
            %
            % See also DESIGNFILT DIGITALFILTER FILTFILT

            % hard error if we aren't given a digitalFilter
            assert(isa(D, 'digitalFilter'), "Expected a 'digitalFilter' but got a " + class(D) + " instead.");

            % defaults
            if nargin < 3, dim = chd.tdim; end

            % filter: always applied in dim 1
            chd = applyFun2Dim(chd, @(x) cast(filtfilt(D, double(x)), 'like', x), dim);
        end
        function chd = fftfilt(chd, D, dim)
            % FFTFILT Filter data with a digitalFilter
            %
            % chd = FFTFILT(chd, D) filters the channel data with the
            % digitalFilter D. Use DESIGNFILT to design a digital filter
            %
            % chd = FFTFILT(chd, D, dim) applies the filter in dimension
            % dim. The default is the time dimension.
            %
            % See also DESIGNFILT DIGITALFILTER FFTFILT

            % hard error if we aren't given a digitalFilter
            assert(isa(D, 'digitalFilter'), "Expected a 'digitalFilter' but got a " + class(D) + " instead.");
            
            % defaults
            if nargin < 3, dim = chd.tdim; end

            % filter: always applied in dim 1
            chd = applyFun2Dim(chd, @(x) reshape(cast(fftfilt(D, double(x(:,:))), 'like', x), size(x)), dim);

            % adjust time axes
            if dim == chd.tdim
                chd.t0 = chd.t0 - (D.FilterOrder-1)/2/chd.fs;
            end
        end
        function chd = hilbert(chd, N, dim)
            % HILBERT - overloads the hilbert function
            %
            % chd = HILBERT(chd) applies the hilbert function to the data
            % in the time dimension.
            %
            % chd = HILBERT(chd, N) computes the N-point Hilbert transform. 
            % The data is padded with zeros if it has less than N points, 
            % and truncated if it has more.
            %
            % chd = HILBERT(chd, N, dim) operates in dimension dim. The
            % default is the time dimension.
            %
            % See also HILBERT

            % parse inputs
            if nargin < 3, dim = chd.tdim; end
            if nargin < 2 || isempty(N), N = size(chd.data, dim); end

            % use MATLAB's optimized implementation
            % chd = applyFun2Dim(chd, @hilbert, chd.tdim, varargin{:});
            
            % apply natively to support half type
            chd = fft(chd, N, chd.tdim); % send to freq domain
            Nd2 = floor(N/2); % number of postive/negative frequencies
            w = [1; 2*ones([Nd2-1,1]); 1+mod(N,2); zeros([N-Nd2-1,1])]; % hilbert weight vector
            chd.data = chd.data .* shiftdim(w, 1-chd.tdim); % apply 
            chd = ifft(chd, N, chd.tdim);
        end
        function chd = fft(chd, N, dim)
            % FFT - overload of fft
            %
            % chd = FFT(chd) computes the fft of the channel data along the 
            % time axis. The time and frequency axes are unchanged.
            %
            % chd = FFT(chd, N) computes the N-point fft.
            %
            % chd = FFT(chd, N, dim) or FFT(chd, [], dim) operates along
            % dimension dim.
            %
            % See also FFT CHANNELDATA/FFTSHIFT

            % defaults
            if nargin < 3, dim = chd.tdim; end
            if nargin < 2 || isempty(N), N = size(chd.data, dim); end
            if istall(chd) && dim == 1, error('Cannot compute fft in the tall dimension.'); end
            chd = copy(chd);
            chd.data = matlab.tall.transform(@fft, chd.data, N, dim); % take the fourier transform
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

            if nargin < 2, dim = chd.tdim; end
            if istall(chd) && dim == 1, error('Cannot compute fftshift in the tall dimension.'); end
            chd = copy(chd);
            chd.data = matlab.tall.transform(@fftshift, chd.data, dim);
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
            if nargin < 3, dim = chd.tdim; end
            if nargin < 2 || isempty(N), N = size(chd.data, dim); end
            if istall(chd) && dim == 1, error('Cannot compute ifft in the tall dimension.'); end
            chd = copy(chd);
            chd.data = matlab.tall.transform(@ifft, chd.data, N, dim); % take the fourier transform
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
            if nargin < 2, dim = chd.tdim; end
            if istall(chd) && dim == 1, error('Cannot compute ifftshift in the tall dimension.'); end
            chd = copy(chd);
            chd.data = matlab.tall.transform(@ifftshift, chd.data, dim);
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
            [Tt, Tf, Td] = deal(chd.t0, chd.fs, cast(zeros([0,0]), 'like', chd.data));
            
            % Make new ChannelData (avoid modifying the original)
            chd = copy(chd);

            % ensure numeric args are non-sparse, double
            chd = (doubleT(chd)); % data is type double
            fs = (double(fs)); % new frequency is type double
            inum = cellfun(@isnumeric, varargin); % all numeric inputs are type double
            varargin(inum) = cellfun(@double, varargin(inum), 'UniformOutput', false);

            % resample in time - no support for other dims: fs is required arg
            % [y, ty] = resample(chd.data, chd.time, fs, varargin{:}, 'Dimension', chd.tdim);
            % [chd.fs, chd.t0, chd.data] = deal(fs, ty(1), y);
            y = matlab.tall.transform(@resample, chd.data, chd.time, fs, varargin{:}, 'Dimension', chd.tdim);
            [chd.fs, chd.data] = deal(fs, y);

            % cast back to original type
            tmp = cellfun(@(x,T) cast(x, 'like', T), {chd.t0, chd.fs, chd.data}, {Tt, Tf, Td}, 'UniformOutput', false);
            [chd.t0, chd.fs, chd.data] = tmp{:};

        end
        function chd = real(chd)    , chd = applyFun2Data (chd, @real); end
        function chd = imag(chd)    , chd = applyFun2Data (chd, @imag); end
        function chd = abs(chd)     , chd = applyFun2Data (chd, @abs); end
        function chd = angle(chd)   , chd = applyFun2Data (chd, @angle); end
        function chd = rad2deg(chd) , chd = applyFun2Data (chd, @rad2deg); end
        function chd = deg2rad(chd) , chd = applyFun2Data (chd, @deg2rad); end
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
            % chd.data(end+(B+A),:) = 0; % append A + B zeros in time to the data
            s = repmat({':'}, [1,gather(ndims(chd.data))]); % splice in all other dimensions
            s{chd.tdim} = chd.T + (1:(B+A)); % expand by B+A in time dimension
            chd.data = subsasgn(chd.data,substruct('()',s),0); % call assignment - set to zero
            chd.data = matlab.tall.transform(@circshift, chd.data, B, chd.tdim); % shift B of the zeros to the front
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
            y = fft(chd, [], chd.tdim); % compute fft
            z = argmax(abs(y.data), [], chd.tdim); % get peak over frequencies
            z = matlab.tall.transform(@mode, z, 2:gather(ndims(y.data))); % mode over non-tall dims
            fc = f(median(gather(z))); % select median over tall dims
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

            if isscalar(chd.t0), chd = copy(chd); end % short-circuit
            if nargin < 2, interp = 'cubic'; end % default interpolation
            if nargin < 3, t0_ = min(chd.t0, [], 'all'); end % get global start time
            toff = chd.t0 - t0_; % get offset across upper dimensions
            npad = ceil(max(toff,[],'all') * chd.fs); % furthest sample containing data
            chd = zeropad(chd,0,npad); % extend time-axes
            tau = chd.time + toff; % get delays to resample all traces
            y = chd.sample(tau, interp); % resample
            chd.data = y; % make new object
            chd.t0 = t0_; % make new object
        end
    
        function y = sample(chd, tau, interp, w, sdim)
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
            % y = SAMPLE(x, t, interp, w) applies the weighting array via 
            % point-wise  multiplication after sampling the the data. The 
            % dimensions must be compatible with the sampling array t in 
            % the sampling dimension dim. The default is 1.
            % 
            % y = SAMPLE(x, t, interp, w, sdim) sums the data in the 
            % dimension(s) sdim after the weighting matrix has been applied.
            % The default is [] (no dimensions).
            % 
            % See also INTERP1 INTERPD INTERPF WSINTERPD CHANNELDATA/RECTIFYT0

            % defaults
            if nargin < 5, sdim = [];           end
            if nargin < 4, w = 1;               end
            if nargin < 3, interp = 'linear';   end

            % check condition that we can implicitly broadcast
            for d = setdiff(1:gather(ndims(tau)), chd.tdim) % all dims except time must match
                assert(any(size(tau,d) == size(chd.data,d)) || ...
                      any([size(tau,d)  , size(chd.data,d)] == 1), ...
                    'Delay size must match the data size (%i) or be singleton in dimension %i.',...
                    size(chd.data, d), d ...
                    );
            end

            % compute the integer delays (I x [1|N] x [1|M] x [1|F] x ...) (default order)
            ntau = (tau - chd.t0) * chd.fs;

            % dispatch
            if interp ~= "freq" 
                if istall(ntau) || istall(chd.data)
                    y = matlab.tall.transform(@wsinterpd, chd.data, ntau, chd.tdim, w, sdim, interp, 0);
                else
                    y = wsinterpd(chd.data, ntau, chd.tdim, w, sdim, interp, 0);
                end
            elseif interp == "freq"
                % extend data if necessary
                nwrap = min(0,floor(min(ntau,[],'all'))); % amount the signal is wrapped around from negative time
                next  = max(0, ceil(max(ntau,[],'all')) - chd.T-1); % amount the signal is extended in positive time
                chd = zeropad(chd, -nwrap, next); % extend signal to ensure we sample zeros
                y = interpf(chd.data, ntau, chd.tdim, w, sdim); % sample in the frequency domain
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
            % h = IMAGESC(self, m, ax) uses the axes ax instead of the axes
            % returned by gca. 
            %
            % h = IMAGESC(..., 'YData', yax) uses the yax for the
            % y-axis instead of the time domain.
            % 
            % h = IMAGESC(..., Name, Value, ...) passes the following
            % arguments to the imagesc function.
            %
            % Example:
            %
            % % Show the data in the frequency domain in MHz
            % f = chd.fs*(0:chd.T-1)/chd.T; % frequency axis
            % figure; 
            % imagesc(fft(chd), 1, 'YData', 1e-6*f);
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
            dims = gather(max(3, [ndims(self.data)])); % (minimum) number of dimensions
            idims = [self.tdim, self.ndim]; % image dimensions
            fdims = setdiff(1:dims, idims); % frame dimensions
            dsz = gather(size(self.data, fdims)); % full size - data dimensions
            tsz = gather(size(self.time, fdims));

            % we index the data linearly, but the time axes may be
            % implicitly broadcasted: we can use ind2sub to recover it's
            % sizing
            ix = cell([numel(fdims), 1]); % indices of the frame for the data
            [ix{:}] = ind2sub(dsz, gather(m)); % get cell array of indices
            it = gather(min([ix{:}], tsz)); % restrict to size of self.time
            ix = cellfun(@gather, ix, 'UniformOutput', false); % enforce on CPU

            % select the transmit/frame - we can only show first 2
            % dimensions, time x rx
            d = gather(sub(self.data, ix, fdims));
            d = permute(d, [idims, fdims]); % move image dimensions down 

            % choose to show real part or dB magnitude
            if isnumeric(d), d = double(d); end % cast integer types for abs, half types for imagesc
            d = squeeze(double(d)); % make friendly to imagesc
            if ~isreal(d), d = mod2db(d); end % assume complex data should be shown in dB

            % get the time axes for this frame
            t = gather(double(sub(self.time, num2cell(it), fdims)));

            % choose which dimensions to show
            axes_args = {'XData', 1:self.N, 'YData', t}; % ndim, tdim labels

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
            for f = 1:gather(prod(size(self.data,3:max(3,gather(ndims(self.data))))))
                if ~isvalid(ax), break; end % quit if handle destroyed
                h = imagesc(self, f, ax, varargin{:});
                drawnow limitrate; pause(1/20);
            end
        end

        function gif(chd, filename, h, varargin)
            % GIF - Write the ChannelData to a GIF file
            %
            % GIF(chd, filename) writes the ChannelData chd to the file
            % filename.
            %
            % GIF(chd, filename, h) updates the image handle h rather than
            % creating a new image. Use imagesc to create an image handle.
            % You can then format the figure prior to calling this
            % function.
            %
            % GIF(..., Name, Value, ...) forwards Name/Value pairs to
            % imwrite. 
            %
            % Example:
            % sz = [2^8, 2^6, 2^6];
            % x = complex(rand(sz), rand(sz)) - (0.5 + 0.5i);
            % chd = angle(ChannelData('data', x));
            % figure;
            % h = imagesc(chd, 1);
            % colormap hsv; 
            % colorbar;
            % gif(chd, 'random_phase.gif', h, 'DelayTime', 1/20);
            %
            % See also IMAGESC ANIMATE

            % defaults
            kwargs = struct('LoopCount', Inf, 'DelayTime', 1/15);

            % parse inputs
            for i = 1:2:numel(varargin), kwargs.(varargin{1}) = varargin{i+1}; end

            % if no image handle, create a new image
            if nargin < 3, h = imagesc(chd, 1); end
            
            chd = rectifyDims(chd); % put in proper order
            x = chd.data; % all data
            M_ = prod(size(x, 3:min(3,ndims(x)))); % all slices

            % get image frames
            % TODO: there's some weird bug where the size is randomly off 
            % by 10 pixels here? Can I work around it?
            for m = M_:-1:1, h.CData(:) = x(:,:,m); fr{m} = getframe(h.Parent.Parent); end
            
            % get color space for the image
            [~, map] = rgb2ind(fr{1}.cdata, 256, 'nodither');

            % get all image data
            im = cellfun(@(fr) {rgb2ind(fr.cdata,map,'nodither')}, fr);
            im = cat(4, im{:});

            % forward to imwrite
            nvkwargs = struct2nvpair(kwargs);
            imwrite(im, map, filename, nvkwargs{:});
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
        function t = get.time(self), t = cast(self.t0 + shiftdim((0 : self.T - 1)', 1-self.tdim) ./ self.fs, 'like', real(self.data)); end % match data type, except always real
        function T = get.T(self), T = gather(size(self.data,self.tdim)); end
        function N = get.N(self), N = gather(size(self.data,self.ndim)); end
        function M = get.M(self), M = gather(size(self.data,self.mdim)); end
        function n = get.rxs(self), n=cast(shiftdim((1:self.N)',1-self.ndim), 'like', real(self.data)); end
        function m = get.txs(self), m=cast(shiftdim((1:self.M)',1-self.mdim), 'like', real(self.data)); end
    end

    % sizing functions (used to control tall behaviour and reshaping)
    properties(Hidden, Dependent)
        tdim
        ndim
        mdim
    end

    % data indexing functions
    methods
        function chd = join(chds, dim)
            % JOIN - Merge an array of ChannelData objects
            %
            % chd = JOIN(chds, dim) merges the ChannelData array chds into
            % a single ChannelData objects chd in dimension dim of the
            % data. 
            % 
            % The temporal dimension is expanded as needed. All other
            % dimensions must be an identical size or empty. The sampling 
            % frequency must be identical.
            %
            % The start time t0 must be either full in dimension dim or 
            % identical in dimension dim for all ChannelData objects.
            %
            % See also CHANNELDATA/SPLICE

            assert(isalmostn([chds.fs], repmat(median([chds.fs]), [1,numel(chds)]))); % sampling frequency must be identical
            assert(all(string(chds(1).ord) == {chds.ord})); % data orders are identical

            T_ = max([chds.T]); % maximum length of data
            chds = arrayfun(@(chds) zeropad(chds,0,T_ - chds.T), chds); % make all the same length
            chd = ChannelData('t0', cat(dim,chds.t0), 'fs', median([chds.fs]), 'data', cat(dim, chds.data)); % combine
            if all(chd.t0 == sub(chd.t0,1,dim),'all'), chd.t0 = sub(chd.t0,1,dim); end % simplify for identical t0

        end
        function [chds, ix] = splice(chd, dim, bsize)
            % SPLICE - Split the ChannelData into an array of ChannelDatas
            % 
            % chds = SPLICE(chd, dim) returns an array of ChannelData
            % objects chds where each element contains a slice of the data
            % in dimensions dim. It is useful for iterating over
            % ChannelData without explicitly indexing the data and time
            % dimensions manually.
            %
            % chds = SPLICE(chd, dim, bsize) uses a maximum block size of 
            % bsize to partition the ChannelData objects. The default is 1.
            %
            % Example:
            % 
            % % Create data
            % [T, N, M] = deal(2^10, 8, 8);
            % t = (0 : T - 1)'; % time axis
            % f = shiftdim((M/2 : 0.5 : M - 1/2) / M / 2, -1); % frequencies per transmit (normalized)
            % x = sinpi(2 .* f .* t); % sinusoids
            % w = randn([T, N, M]); % white noise
            % chd = ChannelData('data', x + w, 't0', t(1), 'fs', 1);
            % 
            % % Estimate the central frequency for each transmit
            % chds = splice(chd, chd.mdim); % splice per transmit
            % for m = 1:chd.M, fc(m) = estfc(chds(m)); end % estimate per transmit
            % figure; plot(1:chd.M, [f(:), fc(:)], '.-'); 
            % title('Estimated frequency per transmit');
            % legend({'True', 'Estimated'});
            % 
            % See also CHANNELDATA/JOIN SUB
            
            assert(isscalar(dim), 'Dimension must be scalar!'); 
            if nargin < 3, bsize = 1; end

            % S = gather(size(chd.data, dim)); % slices
            St = gather(size(chd.t0  ,dim)); % size in t
            Sx = gather(size(chd.data,dim)); % size in x
            it = num2cell((0:bsize:St-1) + (1:bsize)',1); % index up to maximum t
            ix = num2cell((0:bsize:Sx-1) + (1:bsize)',1); % index up to maximum x
            it{end}(it{end} > St) = []; % restrict to actual max t
            ix{end}(ix{end} > Sx) = []; % restrict to actual max x

            % splice data and time axes
            t = cellfun(@(i) sub(chd.t0  , i, dim), it, 'UniformOutput',false);
            x = cellfun(@(i) sub(chd.data, i, dim), ix, 'UniformOutput',false);
            
            % make array of new ChannelData objects
            chds = repmat(ChannelData('fs', chd.fs, 'ord', chd.ord), [max(numel(t),numel(x)),1]); % new ChannelData objects
            chds = arrayfun(@copy, shiftdim(chds, 1-dim)); % make unique and move to dimension dim
            [chds.t0  ] = deal(t{:}); % set start time(s)
            [chds.data] = deal(x{:}); % set data

        end
        function chd = sub(chd, ind, dim)
            if ~iscell(ind), ind = {ind}; end % enforce cell syntax
            tind = ind; % separate copy for the time indices
            tind(size(chd.time,dim) == 1) = {1}; % set singleton
            has_tdim = any(dim == chd.tdim);
            if has_tdim, 
                n = tind{dim == chd.tdim}; % get the time indices
                assert(issorted(n, 'strictascend') && all(n == round(n)), ... % check the index in time is sorted, continous
                    'The temporal index must be a strictly increasing set of indices.');
            end
            if has_tdim, t = chd.time; else, t = chd.t0; end % index in time if necessary
            t0_   = sub(t, tind, dim); % extract
            data_ = sub(chd.data, ind, dim); % extract
            if has_tdim, t0_ = sub(t0_, 1, chd.tdim); end % extract only the first value in time
            
            chd = copy(chd); % copy semantics
            chd.t0 = t0_; % assign
            chd.data = data_; % assign

        end
        function chd = setOrder(chd, cord)
            assert(...
                length(cord) >= ndims(chd.data),...
                'Number of dimension labels must be greater than or equal to the number of data dimensions.'...
                ); 
            assert(...
                all(ismember('TMN', cord)), ...
                "Dimension labels must contain 'T', 'N', and 'M'."...
                );
            chd.ord = cord; 
        end
        function chd = expandDims(chd, d)
            chd = copy(chd); % copy semantics
            nc = numel(chd.ord); % number of dimension labels
            ccand = setdiff(char(double('F') + (0 : 2*(d-nc))), chd.ord); % get unique labels, starting with 'F'
            chd.ord(nc+1:d) = ccand(1:d-nc);
        end
        function chd = truncateDims(chd)
            nd = numel(chd.ord); % number of (labelled) dimensions
            [~,o] = ismember('TMN', chd.ord); % position of necessary params
            sdims = find(size(chd.data,1:nd) ~= 1); % singleton dims
            kdims = sort(union(o, sdims)); % dims to keep: necessary or non-singleton
            rdims = setdiff(1:nd, kdims); % dims to remove: all others
            chd = permute(chd, [kdims, rdims]); % squeeze data down
            chd.ord = chd.ord(kdims); % remove unnecesary dimensions 
        end
        function chd = swapdim(chd, i, o)
            dord = 1:max([i,o,ndims(chd.data)]); % max possible dim
            dord(i) = o; % swap
            dord(o) = i; % swap
            chd = permute(chd, dord); % move data dimensions
        end
        function chd = permute(chd, dord)
            chd = expandDims(chd, max(dord)); % expand to have enough dim labels
            chd = copy(chd); % copy semantics
            chd.data = permute(chd.data, dord); % change data dimensions
            chd.ord(1:numel(dord)) = chd.ord(dord); % change data order
        end
        function [chd, dord] = rectifyDims(chd)
            % RECTIFYDIMS - Set dimensions to default order
            %
            % chd = RECTIFYDIMS(chd) sets the dimension of the data to
            % their default order of time x receive x transmit x ...
            %

            D = gather(max(numel(chd.ord), ndims(chd.data))); % number of dimensions
            dord = arrayfun(@(o) find(chd.ord == o), 'TNM'); % want this order to start
            dord = [dord, setdiff(1:D, dord)]; % make sure we have all dimensions accounted for
            chd = permute(chd, dord); % set dims to match in lower dimensions
            % chd = truncateDims(chd); % remove unnecessary dimensions
        end
    end
    methods
        function d = get.tdim(self), d = find(self.ord == 'T'); end
        function d = get.ndim(self), d = find(self.ord == 'N'); end
        function d = get.mdim(self), d = find(self.ord == 'M'); end
    end
end

% TODO: make a default interpolation method property