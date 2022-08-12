% WAVEFORM Arbitrary Waveform definition class
%
% 
classdef Waveform < handle
    
    properties(GetAccess=public, SetAccess=public)
        fun function_handle = @(t) sin(2*pi*5e6*t)  % functional form of the waveform
        t0 (1,1) {mustBeNumeric} = -1 / 5e6         % start time
        tend (1,1) {mustBeNumeric} = 1 / 5e6        % end time
    end
    
    properties(GetAccess=public, SetAccess=public)
        samples = []         % sampled form of the waveform
        mode (1,1) string {mustBeMember(mode, ["fun", "samp"])} = 'fun'         % echo mode {*fun | samp}
        dt {mustBeNumeric, mustBeScalarOrEmpty} = []              % time resolution in sampled mode
    end
    
    properties(Dependent)
       time
       duration
    end
    
    % constructor / actions
    methods(Access=public)
        function self = Waveform(varargin)
            % WAVEFORM/WAVEFORM - Waveform constructor
            % 
            % Waveform(Name,Value, ...) constructs a Waveform using
            % name-value pairs to create a Waveform object by specifying a 
            % start/end time t0/tend or time samples t and specifying 
            % either a function valid in the range t0 <= t <= tend or 
            % samples defined on t0 <= t <= tend. Pass arguments as 
            % name/value pairs. The Waveform object operates in *either* 
            % function or sampled mode
            %
            % See also SEQUENCE
            for i = 1:2:nargin
                switch varargin{i}
                    case {'fun', 't0', 'tend'}
                        self.(varargin{i})  = varargin{i+1};
                    case 't'
                        t = varargin{i+1};
                        self.t0           = t(1);
                        self.tend         = t(end);
                        self.dt = mode(diff(t));
                    case 'samples'
                        % check if dt given after
                        if i + 2 <= nargin ...
                                && strcmpi(varargin{i + 2}, 'dt')
                            self.setSampled(varargin{i+1}, 1 / varargin{i+3});
                        else
                            self.setSampled(varargin{i+1});
                        end
                        
                end
            end
        end
        
        function s = toSampled(self, fs)
            % samples the waveform over the valid time period at sampling
            % frequency fs
            
            % if already in sampled mode, do nothing
            if strcmp(self.mode, 'samp'), return, end
            
            % get sample times
            t = self.getSampleTimes(fs);
            
            % sample it
            samps = self.sample(t);
            
            % set the mode (clears the function)
            self.setSampled(samps, fs);
            
            s = self.samples;
        end
        
        function varargout = plot(self, varargin)
            % WAVEFORM/PLOT - plot the Waveform
            % 
            % PLOT(self) plots the Waveform
            %
            % PLOT(self, ax) uses axes ax instead of the current axes
            %
            % PLOT(..., 'freqs', true) toggles whether to also plot the
            % frequency domain of the waveform
            % 
            % PLOT(..., Name, Value, ...) passes name-value pairs to the
            % built-in plot function.
            % 
            % arguments for plot can be passed as with the plot function
            % set the property 'freqs' to true to plot frequencies
            % returns a handle to the plot or to 2 plots if 'freqs' is an
            % input

            % parse inputs
            if nargin >= 2 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                axs = varargin{1}; varargin(1) = [];
            else
                axs = gca;
            end
            
            % search for a 'freqs' argument
            plot_freqs = false;
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'freqs'
                        l = i;
                        plot_freqs = logical(varargin{i+1});
                end
            end
            
            % remove axes and freqs arguments
            if exist('l','var'), varargin(l:l+1) = []; end
            
            % get data
            switch self.mode
                case 'fun'
                    p = nextpow2(1e4);
                    N = 2^p;
                    fs = (N - 1)/(self.tend - self.t0);
                    if isinf(fs), fs = 1e9; end % handle divide by zero
                    t = self.getSampleTimes(fs);
                    v = self.sample(t);
                    w = abs(fftshift(fft(v,N,1),1));
                    % dt = mode(diff(t)); %#ok<PROPLC,CPROPLC>
                    % fs = 1 / dt; %#ok<PROPLC>
                    k = (((0 : (N - 1)) / N ) - 1 / 2) * fs;
                case 'samp'
                    t = self.getSampleTimes();
                    v = self.samples;
                    M = numel(t);
                    % p = nextpow2(M);
                    %echo N = 2^p;
                    N = M;
                    w = abs(fftshift(fft(v,N,1),1));
                    fs = 1 / self.dt;
                    k = (((0 : (N - 1)).' / N ) - 1 / 2) * fs;
            end
            
            % show complex components
            if isreal(v), 
                [v, lbl] = deal(v(:), {'real'});
            else, 
                v = [real(v(:)), imag(v(:)), abs(v(:))];
                lbl = {'real', 'imag', 'abs'};
            end
            t = t(:);
            
            % plot with arguments
            if plot_freqs
                h_time = subplot(2,1,1);
                hp = plot(h_time, t, v, varargin{:});
                xlabel('Time');
                ylabel('Amplitude');
                grid on;
                
                h_freq = subplot(2,1,2);
                hp2 = plot(h_freq, k, w, varargin{:});
                xlabel('Frequency');
                ylabel('Amplitude');
                grid on;
            else
                hp = plot(axs, t, v, varargin{:});
                xlabel(axs, 'Time');
                ylabel(axs, 'Amplitude');
                grid on;
            end
            
            legend(hp, lbl);
            
            if nargout > 0, varargout{1} = hp; varargout{2} = hp2; end
        end
        
        function w = copy(self)
            params = {'t0', self.t0, 'tend', self.tend};        
            switch self.mode
                case 'fun'
                    params = [params(:)', {'fun'}, {self.fun}];
                case 'samp'
                    params = [params(:)', {'samples'}, {self.samples}];
            end
            w = Waveform(params{:});
        end

        function wv = scale(self, kwargs)
            arguments
                self Waveform
                kwargs.time (1,1) double
            end
            % create a new waveform object with the properties scaled
            w = kwargs.time; % factor
            f = self.fun; % function
            if self.mode == "fun"
                wv = Waveform( ...
                    't0', w * self.t0, 'tend', w * self.tend, ...
                    'fun', @(t) f(t./w) ...
                    );
            elseif self.mode == "samp"
                wv = Waveform( ...
                    't0', w * self.t0, 'tend', w * self.tend, ...
                    'dt', w * self.dt ...
                    );
            else 
                error('Undefined state.')
            end
        end
        
        function s = sample(self, t)
            % Get waveform samples for the times specified in the vector t.
            % If in sampled mode, t must match the sample time axis
            switch self.mode
                case 'fun'
                    s = zeros(size(t), 'like', t);
                    n = self.t0 <= t & t <= self.tend;
                    s(n) = self.fun(t(n));
                case 'samp'
                    % make sure all time points are actually sampled
                    t_s = self.getSampleTimes();
                    [m, i] = ismembertol(t, t_s, self.dt / 16);
                    
                    % if all points are sampled, return sampled points
                    if all(m)
                        s = self.samples(i);
                    else
                        error('The requested times do not match the sampled times for this Waveform.');
                    end
            end
        end
        
        function wv = conv(self, other, fs)
            % CONV Waveform convolution
            %
            % wv = CONV(self, other) convoles 2 Waveform objects with
            % matching sampling modes to form another Waveform
            %
            % wv = CONV(self, other, fs) uses an intermediate sampling
            % frequency of fs when sampling the waveform to perform the 
            % convolution.
            %
            % Note: this function needs to undergo more testing
            %
            % See also WAVEFORM/SAMPLE

            assert(isequal(self.mode, other.mode), "Both waveforms must have matching sample type.")
            if self.mode == "samp" % sampled
                if self.dt == other.dt % matching sampling frequency
                wv = Waveform(...
                    'samples', conv(self.samples, other.samples), ...
                    'dt', self.dt, 't0', self.t0 + other.t0, ...
                    'tend', self.tend + other.tend...
                    );
                else
                    error('Cannot convolve waveforms with different sampling freuqencies (not implemented)');
                end
            else
                % make a convolution function, using the sampling freuqency fs
                if self.t0 == self.tend
                    % this is a delta function: use the identity
                    f = @(t) self.sample(0) * other.sample(t);
                elseif other.t0 == other.tend
                    % that is a delta function: use the identity
                    f = @(t) other.sample(0) * self.sample(t);
                else
                    k = (floor((self.t0+other.t0)*fs) : ceil((self.tend+other.tend)*fs)) / fs; % 1 x K
                    f = @(t) reshape(self.sample(t(:) - k) * other.sample(k'), size(t));
                end
                wv = Waveform('fun', f, 't0', self.t0 + other.t0, 'tend', self.tend + other.tend);
            end
        end
    end
    
    % get/set methods
    methods(Access=public)
        function t = getSampleTimes(self, fs)
            switch self.mode
                case 'samp'
                    if nargin ~= 1
                        error('No arguments expected in sampled mode');
                    end
                    t = self.t0 : self.dt : self.tend;
                case 'fun'
                    if nargin ~= 2
                        error('1 argument expected in function mode');
                    end
                    dt = 1/fs; %#ok<PROPLC>
                    t = (ceil(self.t0 / dt) : 1 : floor(self.tend / dt)) * dt; %#ok<PROPLC>
                otherwise
                    error('Waveform in an undefined state.');
            end
            t = t(:);
        end
        
        function setSampled(self, samples, fs)
            if nargin < 3
                self.dt = abs(self.tend - self.t0) / numel(samples);
            else
                self.dt = 1/fs;
            end
            self.samples = samples;
            self.mode = 'samp';
            self.fun = @(t) interp1(self.samples, 1 + (t - self.t0) * fs, "cubic", 0);
        end        
    end    
    
    methods
        function dur = get.duration(self), dur = self.tend - self.t0; end
        function t = get.time(self), t = self.getSampleTimes(); end
    end

    methods(Static)
        function wv = Delta()
            wv = Waveform('t0', 0, 'tend', 0, 'fun', @(t) t == 0);
        end
    end
end

