% WAVEFORM Arbitrary Waveform definition class
%
% The Waveform class stores the function form of a signal and relevant time
% domain so that it can be sampled and resampled as needed by relevant 
% methods. All waveforms must have a start time and end time after which 
% the function goes to 0.
%
% When a Waveform has a sampling frequency, the 'time' and 'samples'
% properties will be available for discrete sampling. The time axis is
% guaranteed to include t == 0 if t0 <= 0 and tend >= 0. 
%
% See also TRANSDUCER SEQUENCE
classdef Waveform < matlab.mixin.Copyable
    
    properties
        % WAVEFORM/FUN - Function defining the waveform
        %
        % FUN is a function_handle to a function that accepts an ND-array 
        % of time samples and returns a corresponding ND-array of samples 
        % at those points in time.
        %
        % See also FUNCTION_HANDLE
        fun {mustBeA(fun, {'function_handle', 'griddedInterpolant'})} = @(t) sin(2*pi*5e6*t)  % functional form of the waveform
        % WAVEFORM/T0 - Start time of the Waveform
        %
        % T0 defines the start of the Waveform. All samples prior to this
        % point are 0.
        %
        % See also WAVEFORM.TEND WAVEFORM.DURATION
        t0 (1,1) {mustBeNumeric} = -1 / 5e6         % start time
        % WAVEFORM/TEND - End time of the Waveform
        %
        % T0 defines the start of the Waveform. All samples after to this
        % point are 0.
        %
        % See also WAVEFORM.T0 WAVEFORM.DURATION
        tend (1,1) {mustBeNumeric} = 1 / 5e6        % end time
        % WAVEFORM/FS - Sampling frequency for the Waveform
        %
        % FS defines a sampling frequency for sampling the Waveform. When
        % defined, the Waveform/dt, Waveform/time, and Waveform/samples
        % properties are available.
        %
        % See also WAVEFORM.DT WAVEFORM.TIME WAVEFORM.SAMPLES
        fs {mustBePositive, mustBeNumeric, mustBeScalarOrEmpty} % sampling frequency
    end
    
    properties(Dependent)
        % WAVEFORM/DURATION - duration of the signal
        %
        % The duration of the signal is portion of the signal that may be 
        % non-zero, defined as Waveform/tend - Waveform/t0.
        %
        % See also WAVEFORM/FS WAVEFORM/TIME WAVEFORM/TEND WAVEFORM/T0

        duration (1,1) % non-zero duration of the signal
        % WAVEFORM/DT - sampling interval
        %
        % The sampling interval is the inverse of the sampling frequency 
        % i.e. DT = 1 / Waveform/fs
        %
        % See also WAVEFORM/FS WAVEFORM/TIME

        dt (1,1) {mustBePositive, mustBeNumeric, mustBeScalarOrEmpty} % sampling interval
        % WAVEFORM/SAMPLES - Waveform samples
        %
        % When the sampling frequency Waveform/fs is set, the 'samples' 
        % property is the set of samples of the Waveform at the times given
        % by Waveform/time.
        % See also WAVEFORM/TIME WAVEFORM/FS

        samples (:,1) {mustBeNumeric} % Waveform samples

        % WAVEFORM/TIME - Waveform time axis
        %
        % When the sampling frequency Waveform/fs is set, the 'time' 
        % property provides a time axis for the non-zero portion of the 
        % signal. The time axis is guaranteed to pass through t == 0.
        %
        % See also WAVEFORM/FS

        time (:,1) {mustBeNumeric}    % Waveform sample times
    end
    properties(Hidden)
        tscale (1,1) double = 1; % scaling factor
        % TLIM - Time axis szie limit
        %
        % TLIM is the maximum size of the time axis. A larger sized time
        % axis will produce an error when calling the time function. This
        % is meant to prevent dependent properties from producing
        % unreasonably large temporary results, such as if for example the 
        % sampling frequency is temporarily set to be 1e6 times larger than
        % it's final state.
        Tlim (1,1) double = 2^13; % limit on size of time axes 
    end
    properties(Dependent, Hidden)
        T % size of the time axis
    end
    
    % constructor / actions
    methods
        function wv = Waveform(kwargs)
            % WAVEFORM - Waveform constructor
            % 
            % wv = Waveform('t0', t0, 'tend', tend, 'fun', fun) creates a
            % Waveform defined from time t0 to time tend with a sample
            % generating function fun.
            % 
            % wv = Waveform('t', t_axis, 'samples', s) creates a Waveform 
            % defined on the regularaly spaced time axis t_axis with 
            % samples s.
            %
            % wv = Waveform(..., 'fs', fs) sets the output sampling 
            % frequency to be fs, and implicitly defines the Waveform/time,
            % Waveform/samples and Waveform/dt properties.
            % 
            % wv = Waveform(..., 'dt', dt) sets the output sampling
            % interval dt = 1/fs. It has the same effect as setting the
            % sampling frequency.
            %
            % Example:
            % 
            % % Create a 5Hz sinusoid
            % wv1 = Waveform('t0', -1.0, 'tend', 1.0, 'fun', @(t)  cospi(10*t ));
            % 
            % % Create a symmetric comb window from sampled data
            % wv2 = Waveform('t', -0.3:0.1:0.3, 'samples', [1 0 1 0 1 0 1])
            %
            % % Set the sampling frequency of the signal to 25Hz
            % [wv1.fs, wv2.fs] = deal(25);
            %
            % % Convolve with an intermediate upsampling frequency of 100Hz
            % wvc = conv(wv1, wv2, 100);
            %
            % % Plot the results
            % figure;
            % hold on;
            % plot(wv1, '.-', 'DisplayName', 'Sinusoid');
            % plot(wv2, '.-', 'DisplayName', 'Comb');
            % plot(wvc, '.-', 'DisplayName', 'Convolution');
            % 
            % See also SEQUENCE TRANSDUCER
            arguments
                kwargs.fun {mustBeA(kwargs.fun, {'function_handle', 'griddedInterpolant'})}
                kwargs.t0 (1,1) {mustBeNumeric}
                kwargs.tend (1,1) {mustBeNumeric}
                kwargs.t (:,1) {mustBeNumeric}
                kwargs.samples (:,1) {mustBeNumeric}
                kwargs.dt {mustBeNumeric, mustBeScalarOrEmpty}
                kwargs.fs {mustBeNumeric, mustBeScalarOrEmpty}
            end

            % set time and functions
            for f = string(fieldnames(kwargs))'
                switch f
                    case {'fun', 't0', 'tend', 'dt', 'fs'}
                        wv.(f)  = kwargs.(f);
                    case 't'
                        t = kwargs.(f);
                        wv.t0   = t(1);
                        wv.tend = t(end);
                        wv.dt   = mode(diff(t));
                end
            end
                    
            % make a function out of the samples if given
            if isfield(kwargs, 'samples')
                if isfield(kwargs, 't')
                    wv.fun = griddedInterpolant(kwargs.t, gather(kwargs.samples), 'cubic', 'none');
                elseif all(isfield(kwargs, {'t0', 'dt', 'tend'}))
                    t = kwargs.t0 : kwargs.dt : kwargs.tend;
                    wv.fun = griddedInterpolant(t, gather(kwargs.samples), 'cubic', 'none');
                elseif all(isfield(kwargs, {'t0', 'fs', 'tend'}))
                    t = kwargs.t0 : 1/kwargs.fs : kwargs.tend;
                    wv.fun = griddedInterpolant(t, gather(kwargs.samples), 'cubic', 'none');
                else
                    wv.fun = griddedInterpolant(wv.time, gather(kwargs.samples), 'cubic', 'none');
                end
            end
        end
        
        function varargout = plot(wv, varargin, kwargs, plot_args)
            % WAVEFORM/PLOT - plot the Waveform
            % 
            % PLOT(wv) plots the Waveform
            %
            % PLOT(wv, ax) uses axes ax instead of the current axes
            %
            % PLOT(..., 'freqs', true) toggles whether to also plot the
            % frequency domain of the waveform
            % 
            % PLOT(..., Name, Value, ...) passes name-value pairs to the
            % built-in plot function.
            %
            % h = PLOT(...) returns a handle to the plot of the signal in
            % the time-domain.
            %
            % [ht, hf] = PLOT(..., 'freqs', true, ...) also returns a
            % handle to the plot of the signal in the frequency domain.
            % 
            % Example:
            % 
            % % Create a 5Hz sinusoid
            % wv1 = Waveform('t0', -1.0, 'tend', 1.0, 'fun', @(t)  cospi(10*t ));
            % % Create a 0.2 second rect window
            % wv2 = Waveform('t0', -0.1, 'tend', 0.1, 'fun', @(t) ones(size(t)));
            %
            % % Set the sampling frequency of the signal to 25Hz
            % [wv1.fs, wv2.fs] = deal(25);
            %
            % % Convolve with an intermediate upsampling frequency of 100Hz
            % wvc = conv(wv1, wv2, 100);
            %
            % % Plot the results
            % figure;
            % plot(wvc, '.-');
            % 
            % See also PLOT

            % parse inputs
            arguments
                wv Waveform
            end
            arguments(Repeating)
                varargin
            end
            arguments
                kwargs.freqs (1,1) logical = false
                plot_args.?matlab.graphics.chart.primitive.Line
            end

            % extract axis and other non-Name/Value pair arguments
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes')
                axs = varargin{1}; varargin(1) = [];
            else, axs = gca;
            end

            % get data
            if isempty(wv.fs)
                wv = copy(wv); % don't modify the original
                N = 2^nextpow2(1e4);
                fs_ = (N - 1)/(wv.duration);
                if isinf(fs_), fs_ = 1e9; end % handle divide by zero
            else
                N = 2^nextpow2(max(1e4,wv.T)); % use at least 1e4-point DFT
                fs_ = wv.fs; % use original sampling frequency
            end
            wv = copy(wv); 
            wv.Tlim = N+1;
            wv.fs = fs_;
            t = wv.time;
            v = wv.samples;
            w = abs(fftshift(fft(v,N,1),1));
            k = (((0 : (N - 1)) / N ) - 1 / 2) * fs_;
            
            % show complex components
            if isreal(v) 
                [v, lbl] = deal(v(:), {'real'});
            else
                v = [real(v(:)), imag(v(:)), abs(v(:))];
                lbl = {'real', 'imag', 'abs'};
            end
            t = t(:);
            
            % plot with arguments
            plot_args = struct2nvpair(plot_args);
            if kwargs.freqs
                h_time = subplot(2,1,1);
                hp = plot(h_time, t, v, varargin{:}, plot_args{:});
                xlabel('Time');
                ylabel('Amplitude');
                grid on;
                
                h_freq = subplot(2,1,2);
                hp2 = plot(h_freq, k, w, varargin{:}, plot_args{:});
                xlabel('Frequency');
                ylabel('Amplitude');
                grid on;
            else
                hp = plot(axs, t, v, varargin{:}, plot_args{:});
                xlabel(axs, 'Time');
                ylabel(axs, 'Amplitude');
                grid on;
            end
            
            legend(hp, lbl);
            
            if nargout > 0, varargout{1} = hp; varargout{2} = hp2; end
        end
        
        function wv = scale(wv, kwargs)
            % SCALE - Scale the units of the Waveform
            %
            % wv = SCALE(wv, 'time', tscale) scales the values of
            % time by tscale, and values of (temporal) frequency by the
            % inverse of tscale.
            %
            % Example:
            % % Define a system in meters, seconds
            % wv = Waveform('t0', -1/5e6, 'tend', 1/5e6, 'fun', @(t)cos(2i*pi*5e6*t)); % defined in meters, seconds
            %
            % % Scale the values to millimiters, microseconds
            % wv = scale(wv, 'time', 1e6);
            % 
            %
            arguments
                wv Waveform
                kwargs.time (1,1) double
            end
            wv = copy(wv);
            [wv.t0, wv.tend] = dealfun(@(x) kwargs.time * x, wv.t0, wv.tend);
            if ~isempty(wv.dt), wv.dt = wv.dt * kwargs.time; end
            wv.tscale = wv.tscale ./ kwargs.time;
        end
        
        function s = sample(wv, t)
            % SAMPLE - Sample the Waveform
            %
            % s = SAMPLE(wv, t) samples the Waveform wv at times t.
            %
            % To ensure proper sampling, set the sampling frequency i.e. 
            % the 'fs' property and use the 'samples' property.
            %
            % Example:
            % 
            % % Create a waveform (10 periods @ 5Hz)
            % wv = Waveform('t0', -1, 'tend', 1, 'fun', @(t)cos(2i*pi*5*t));
            %
            % % Re-sample a continuous waveform with offsets
            % wv.fs = 100; % 100 Hz
            % toff = (rand([wv.T, 1]) - 0.5) / 5 / 20; % pi/10 phase noise
            % sig = wv.sample(wv.time + toff); % sample with phase noise
            % 
            % % Plot
            % figure;
            % plot(wv.time, wv.samples, '.-', 'DisplayName', 'Noiseless');
            % plot(wv.time, sig       , '.-', 'DisplayName', 'Noisy');
            %
            % 
           arguments
               wv Waveform
               t {mustBeNumeric}
           end
           
            % Get waveform samples for the times specified in the vector t.
            tu = unique(t); % unique values of t - we only need to sample these
            [~, i] = ismember(t, tu); % indices where t matches tu 
            n = wv.t0 <= tu & tu <= wv.tend; % non-OOB indices
            f = zeros(size(tu), 'like', tu);
            f(n) = wv.fun(wv.tscale .* tu(n)); % compute for unique, non-OOB
            s = reshape(f(i),size(t)); % map output
        end
        
        function wv = conv(this, that, fs)
            % CONV Waveform convolution
            %
            % wv = CONV(this, that) convoles 2 Waveform objects with
            % matching sampling modes to form another Waveform
            %
            % wv = CONV(this, that, fs) uses an intermediate sampling
            % frequency of fs when sampling the waveform to perform the 
            % convolution. The default is the maximum sampling frequency.
            %
            % Example:
            % 
            % % Create a 5Hz sinusoid
            % wv1 = Waveform('t0', -1.0, 'tend', 1.0, 'fun', @(t)  cospi(10*t ));
            % % Create a 0.2 second rect window
            % wv2 = Waveform('t0', -0.1, 'tend', 0.1, 'fun', @(t) ones(size(t)));
            %
            % % Set the sampling frequency of the signal to 25Hz
            % [wv1.fs, wv2.fs] = deal(25);
            %
            % % Convolve with an intermediate upsampling frequency of 100Hz
            % wvc = conv(wv1, wv2, 100);
            %
            % See also WAVEFORM/SAMPLE

            arguments
                this (1,1) Waveform
                that (1,1) Waveform
                fs double = max([this.fs, that.fs])
            end

            % make a convolution function, using the sampling freuqency fs
            if this.t0 == this.tend
                % this is a delta function: use the identity
                f = @(t) this.sample(0) * that.sample(t);
            elseif that.t0 == that.tend
                % that is a delta function: use the identity
                f = @(t) that.sample(0) * this.sample(t);
            else % this and that are not deltas: resample and convolve via matrix multiplication
                k = (floor((this.t0+that.t0)*fs) : ceil((this.tend+that.tend)*fs)) / fs; % 1 x K
                f = @(t) reshape(this.sample(t(:) - k) * that.sample(k'), size(t));
            end
            wv = Waveform('fun', f, 't0', this.t0 + that.t0, 'tend', this.tend + that.tend, 'fs', max([this.fs, that.fs]));
        end
    end
    
    methods
        function s = get.samples(wv), s = wv.sample(wv.time); end
        function set.samples(wv, s), wv.fun = griddedInterpolant(wv.time, s, 'cubic', 'none'); end
        function dur = get.duration(wv), dur = wv.tend - wv.t0; end
        function t = get.time(wv)
            if wv.T > wv.Tlim, error("Time axis too large to compute. If more than " + wv.Tlim + " values are required, increase the 'Tlim' property."); end
            t = (floor(wv.t0 * wv.fs) : 1 : ceil(wv.tend * wv.fs))' .* wv.dt;            
        end
        function T_ = get.T(wv), T_ = 1 + ceil(wv.tend * wv.fs) - floor(wv.t0 * wv.fs);  end 
        function dt = get.dt(wv), dt = 1./wv.fs; end
        function set.dt(wv, dt), wv.fs = 1./dt; end
    end

    methods(Static)
        function wv = Delta()
            % DELTA - Delta function Waveform constructor
            %
            % wv = WAVEFORM.DELTA() constructs a Delta function waveform
            %
            % 
            wv = Waveform('t0', 0, 'tend', 0, 'fun', @(t) t == 0);
        end
    end
end

