% CHANNELDATA - Store channel data and associated properties
%
% The ChannelData class stores a datacube and it's temporal axes and
% provides overloaded methods for convenience in manipulating and plotting 
% the data. The ChannelData must have a t0 value consistent with the
% definition of t0 in QUPS to be used with higher beamforming algorithms in
% QUPS.
%
% All data is assumed to be sampled as the same sampling freuqency, but the
% start time t0 may vary across transmits.
%
% See also SEQUENCE TRANSDUCER

classdef ChannelData < matlab.mixin.Copyable

    properties
        data    % channel data (T x N x M)
        t0 = 0  % start time (T x 1 x [1|M])
        fs      % sampling freuqency
    end
    properties(SetAccess=protected,GetAccess=public)
        ord = 'TNM'; % data order
    end
    properties (Dependent)
        time    % time axis (T x 1 x [1|M])
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
    
    % copy
    methods(Access=protected)
        function chd = copyElement(self)
            chd = ChannelData('fs', self.fs,'t0',self.t0,'data',self.data);
        end
    end

    % conversion functions
    methods
        function channel_data = getUSTBChannelData(self, sequence, xdc)
            % GETUSTBCHANNELDATA - Create a USTB channel data object
            % 
            % channel_data = getUSTBChannelData(self) creates a USTB channel data
            % object from the given data parameters
            
            channel_data = uff.channel_data(...
                'sampling_frequency', self.fs, ...
                'sound_speed', sequence.c0, ...
                'sequence', sequence.getUSTBSequence(xdc, self.t0), ...
                'probe', xdc.getUSTBProbe(), ...
                'data', self.data ...
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
        function chd = gpuArray(chd), chd = applyFun2Props(chd, @gpuArray); end
        function chd = tall(chd)    , chd = applyFun2Data (chd, @tall); end
        function chd = sparse(chd)  , chd = applyFun2Data (chd, @sparse); end
        function chd = doubleT(chd) , chd = applyFun2Props(chd, @double); end
        function chd = singleT(chd) , chd = applyFun2Props(chd, @single); end
        function chd =  int64T(chd) , chd = applyFun2Data (chd, @int64); end
        function chd = uint64T(chd) , chd = applyFun2Data (chd, @uint64); end
        function chd =  int32T(chd) , chd = applyFun2Data (chd, @int32); end
        function chd = uint32T(chd) , chd = applyFun2Data (chd, @uint32); end
        function chd =  int16T(chd) , chd = applyFun2Data (chd, @int16); end
        function chd = uint16T(chd) , chd = applyFun2Data (chd, @uint16); end
        function chd =  int8T(chd)  , chd = applyFun2Data (chd, @int8); end
        function chd = uint8T(chd)  , chd = applyFun2Data (chd, @uint8); end
        function T = classUnderlying(self), try T = classUnderlying(self.data); catch, T = class(self.data); end, end % revert to class if undefined
        function T = underlyingType(self), try T = underlyingType(self.data); catch, T = class(self.data); end, end % R2020b+ overload
        function tf = isreal(self), tf = isreal(self.data); end
    end
    
    % implicit casting: functions that request a numeric type may call
    % these functions
    methods
        function x = double(chd), x = double(chd.data); end
        function x = single(chd), x = single(chd.data); end
        function x =  int64(chd), x =  int64(chd.data); end
        function x = uint64(chd), x = uint64(chd.data); end
        function x =  int32(chd), x =  int32(chd.data); end
        function x = uint32(chd), x = uint32(chd.data); end
        function x =  int16(chd), x =  int16(chd.data); end
        function x = uint16(chd), x = uint16(chd.data); end
        function x =   int8(chd), x =   int8(chd.data); end
        function x =  uint8(chd), x =  uint8(chd.data); end
    end

    % DSP overloads 
    methods
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
            % CHANNELDATA/FFT - overload of fft
            %
            % chd = FFT(chd) computes the fft of the channel data along the 
            % time axis. 
            %
            % See FFT

            % defaults
            if nargin < 3, dim = 1; end
            if nargin < 2, N = []; end
            chd = copy(chd);
            chd.data = fft(chd.data, N, dim); % take the fourier transform
        end
        function chd = fftshift(chd, dim)
            % CHANNELDATA/FFTSHIFT - overload of fftshift
            %
            % See FFTSHIFT
            if nargin < 2, dim = 1; end
            chd = copy(chd);
            chd.data = fftshift(chd.data, dim);
        end
        function chd = ifft(chd, N, dim)
            % CHANNELDATA/IFFT - overload of ifft
            %
            % See IFFT

            % defaults
            if nargin < 3, dim = 1; end
            if nargin < 2, N = []; end
            chd = copy(chd);
            chd.data = ifft(chd.data, N, dim); % take the fourier transform
        end
        function chd = ifftshift(chd, dim)
            % CHANNELDATA/IFFTSHIFT - overload of ifftshift
            %
            % See IFFTSHIFT
            if nargin < 2, dim = 1; end
            chd = copy(chd);
            chd.data = ifftshift(chd.data, dim);
        end
        function chd = resample(chd, fs, varargin)
            % CHANNELDATA/RESAMPLE - Overload of resample
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
            % CHANNELDATA/ZEROPAD - Zero pad the data in time
            %
            % chd = ZEROPAD(chd, B) prepends B zeros to the ChannelData 
            % data in time
            %
            % chd = ZEROPAD(chd, B, A) also appends A zeros to the 
            % ChannelData data in time
            %
            % See also CIRCSHIFT

            if nargin < 2 || isempty(B), B = 0; end
            if nargin < 3 || isempty(A), A = 0; end
            assert(A >= 0 && B >= 0);

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
            % y = SAMPLE(chd, tau) samples the data at the times given by 
            % the delays tau. 
            %
            % y = SAMPLE(chd, tau, interp) specifies the interpolation
            % method. For a gpu, currently ["linear", "nearest", "cubic"] 
            % are supported. On a cpu, support is provided by interpn
            %
            % tau must have broadcastable sizing in the non-temporal
            % dimensions.
            %
            % See also INTERPN CHANNELDATA/RECTIFYT0

            if nargin < 3, interp = 'linear'; end

            % get sizing in first 3 dims (matching dimensions)
            [T_, N_, M_] = size(tau, 1:3); % T' x [1|N] x [1|M] x ...

            % check implicit broadcasting compatability
            assert(N_ == 1 || N_ == chd.N);
            assert(M_ == 1 || M_ == chd.M);

            % pre-allocate
            y = zeros(size(tau), 'like', chd.data); 
            
            % dispatch
            assert(chd.tdim == 1); % TODO: generalize this restriction if necessary
            
            % TODO: use {inner|outer|matching} approach to allow broadcast sampling
            % for the native implementation, iterate over matching 
            % dimensions and use interp1 to sample in time (the sole inner
            % dimension). sub can now be used to index multiple matching
            % dimensions. For the GPU kernel, the variables of iteration
            % must be created for having the outer dimensions in both the
            % data and the sampling times.

            % TODO: implement frequency domain sinc interpolation as in
            % UltrasoundSystem/FocusTX to be used as an interp option
            switch class(chd.data)
                case "gpuArray"
                    ntau = (tau - chd.t0) * chd.fs;
                    mxdim = max(ndims(chd.data), ndims(ntau)); % max dims
                    if all(size(chd.data,2:mxdim) <= size(tau, 2:mxdim)) % if time delays defined across all dimensions of the data
                        y(:) = interpd(chd.data, ntau, chd.tdim, interp); % sample all at once
                    else % TODO: implement implicit broadcast across receives/transmits for tau within interpd
                        for m = chd.M:-1:1
                            for n = chd.N:-1:1
                                y(:,n,m,:) = interpd(chd.data(:,n,m), ntau(:,min(n,N_),min(m,M_),:), chd.tdim, interp); % easy iteration
                            end
                        end
                    end

                otherwise % interpolate, iterate over n,m
                    % TODO: switch to 2D interpolation for MATLAB optimization
                    [t, x] = deal(chd.time, chd.data); % splice
                    [Nt, Mt] = deal(size(chd.t0,2), size(chd.t0,3)); % sizing
                    for m = 1:chd.M
                        for n = 1:chd.N
                            y(:,n,m,:) = interpn(t(:,min(n,Nt),min(m,Mt)), x(:,n,m), tau(:,min(n,N_),min(m,M_),:), interp, 0); 
                        end
                    end
            end
        end
        
        function chd = alignInt(chd, interp)
            % ALIGNINT - Align data to integer sampling
            %
            % chd = ALIGNINT(chd, interp) returns a ChannelData object
            % resampled to have use an integer time axis. This can be used
            % to compress the data to an integer type, but erases
            % information about the true sampling frequency. 
            % It can be useful if you want to store many ChannelData 
            % objects and only cast them to a floating-point type when they
            % are being processed.
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
            % CHANNELDATA/IMAGESC - Overload of imagesc function
            %
            % h = imagesc(self, m) displays transmit m of the channel data 
            % and returns the handle h.
            %
            % h = imagesc(self, m, 'YData', yax) uses the yax for the
            % y-axis instead of the time domain.
            % 
            % h = imagesc(self, m, Name, Value, ...) passes the followin
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
            if isinteger(d), d = double(d); end % cast integer types for abs
            if ~isreal(d), d = mod2db(d); end

            % get the time axes for this channel
            t = sub(self.time, num2cell(i), [3:dims]); %#ok<NBRAK> 

            % choose which dimensions to show
            axes_args = {'XData', 1:self.N, 'YData', t};

            % show the data
            h = imagesc(ax, d, axes_args{:}, varargin{:});
        end

        function h = animate(self, varargin)
            % no halp :(

            if nargin >= 2 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                ax = varargin{1}; varargin(1) = [];
            else
                ax = gca;
            end

            % now use the handle only
            for f = 1:prod(size(self.data,3:ndims(self.data)))
                if ~isvalid(ax), break; end % quit if handle destroyed
                h = imagesc(self, f, ax, varargin{:});
                drawnow limitrate; pause(1/5);
            end
        end
    end

    % dependent methods
    methods
        function set.time(self, t)
            % set.time - set the time axis
            %
            % time must be (T x 1 x [1|M])

            % TODO: warn or error if time axis is not regularly spaced

            %  get the possible sampling freuqencies
            fs_ = unique(diff(t,1));

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
