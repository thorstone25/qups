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
        function chd = single(chd)  , chd = applyFun2Props(chd, @single); end
        function chd = double(chd)  , chd = applyFun2Props(chd, @double); end
        function chd = gpuArray(chd), chd = applyFun2Props(chd, @gpuArray); end
        function chd = tall(chd)    , chd = applyFun2Data (chd, @tall); end
        function chd = sparse(chd)  , chd = applyFun2Data (chd, @sparse); end
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

            % Make new ChannelData (avoid modifying the original)
            chd = copy(chd);

            % save original type
            [Tt, Tf, Td] = deal(chd.t0, chd.fs, chd(1));
            
            % ensure numeric args are non-sparse, double
            chd = double(chd);
            fs = double(fs);
            inum = cellfun(@isnumeric, varargin);
            varargin(inum) = cellfun(@double, varargin(inum), 'UniformOutput', false);

            % resample in time - no support for other dims: fs required
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

            if nargin < 2, B = 0; end
            if nargin < 3, A = 0; end
            assert(A >= 0 && B >= 0);

            % append A + B zeros to the data
            chd.data(end+(B+A),:,:) = 0;
            % shift B of the zeros to the front
            chd.data = circshift(chd.data, B, 1);
            % shift start time
            chd.t0 = chd.t0 - B ./ chd.fs;
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


            f = chd.fs * ((0:chd.T-1) ./ chd.T); % compute frequency axis
            y = fft(chd, [], 1); % compute fft
            fc = f(mode(argmax(abs(y.data), [], 1), 'all')); % select peak over receives/transmits
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

            % select the channel
            d = sub(self.data, m, 3);

            % choose to show real part or dB magnitude
            if ~isreal(self.data), d = mag2db(abs(d)); end

            % get the time axes for this channel
            t = sub(self.time, min(m, size(self.time,3)), 3);

            % choose which dimensions to show
            axes_args = {'XData', 1:self.N, 'YData', t};

            % show the data
            h = imagesc(ax, d, axes_args{:}, varargin{:});
        end
    end

    % dependent methods
    methods
        function set.time(self, t)
            % set.time - set the time axis
            %
            % time must be (T x 1 x [1|M])

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
        function t = get.time(self), t = cast(self.t0 + (0 : gather(self.T) - 1)' ./ self.fs, 'like', real(self.data(1))); end % match data type, except always real
        function T = get.T(self), T = size(self.data,1); end
        function N = get.N(self), N = size(self.data,2); end
        function M = get.M(self), M = size(self.data,3); end
        function n = get.rxs(self), n=cast(shiftdim(1:self.N,0 ), 'like', real(self.data(1))); end
        function m = get.txs(self), m=cast(shiftdim(1:self.M,-1), 'like', real(self.data(1))); end
    end    
end
