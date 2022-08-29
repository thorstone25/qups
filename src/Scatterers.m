% SCATTERERS - Point scatterer definition class
%
% The Scatterers class stores definitions and provides convenience methods
% representing a collection of point scatterers within a homogeneous
% medium.
%
% The Scatterers definition can be used to simulate ChannelData from an
% UltrasoundSystem using a point scatterer simulation call, such as
% the UltrasoundSystem/greens method.
%
% See also MEDIUM ULTRASOUNDSYSTEM/GREENS CHANNELDATA
classdef Scatterers < matlab.mixin.Copyable
   
    properties
        % POS - positions of the scatterers
        %
        % POS defines the positions of the scatterers as a 3 x S array
        % in 3D cartesian coordinates. The dimensions are in the order
        % {'x', 'y', 'z'}.
        %
        % See also SCATTERERS/AMP SCATTERERS/C0 SCATTERERS/BOUNDS
        pos (3,:) double = [0;0;30e-3] % positions for scatterers
        % AMP - amplitude of the scatterers
        %
        % AMP defines the amplitude of the reflection from each scatterer.
        % It must be a vector with the same number of scatterers as are
        % defined by the Scatterers/pos property.
        %
        % See also SCATTERERS/POS SCATTERERS/C0
        amp (1,:) double = 1 % amplitude for the scatterers
        % C0 - ambient sound speed
        %
        % C0 defines the speed of sound of the homogeneous medium on which
        % the scatterers are defined.
        %
        % See also SCATTERERS/POS SCATTERERS/AMP SCATTERERS/ALPHA0
        c0  (1,1) double = 1540; % sound speed of the Medium

        % ALPHA0 - ambient attenuation
        %
        % ALPHA0 defines an attenuation constant for simulation routines 
        % that support one.
        %
        % See also SCATTERERS/C0 SCATTERERS/POS SCATTERERS/AMP
        alpha0 (1,1) double = nan; % sound speed of the Medium
    end

    properties(Dependent, Hidden)
        % XB - x boundaries of the scatterer positions
        %
        % XB gives the minimum and maximum x-value for the
        % Scatterers' positions.
        %
        % See also SCATTERERS/BOUNDS SCATTERERS/YB SCATTERERS/ZB
        xb (1,2) double % x bounds
        % YB - y boundaries of the scatterer positions
        %
        % YB gives the minimum and maximum y-value for the
        % Scatterers' positions.
        %
        % See also SCATTERERS/BOUNDS SCATTERERS/XB SCATTERERS/ZB
        yb (1,2) double % y bounds
        % ZB - z boundaries of the scatterer positions
        %
        % ZB gives the minimum and maximum z-value for the
        % Scatterers' positions.
        %
        % See also SCATTERERS/BOUNDS SCATTERERS/YB SCATTERERS/XB
        zb (1,2) double % z bounds
    end
    properties(Dependent)
        % NUMSCAT - number of scatterers
        %
        % NUMSCAT gives the number of points represented by the Scatterers
        % object. If the positions and amplitudes are not properly
        % initialized, this property will error.
        %
        % See also SCATTERERS/POS SCATTERERS/AMP
        numScat (1,1) double % number of scatterers
        % BOUNDS - boundaries of the scatterer positions
        %
        % BOUNDS gives the minimum and maximum coordinates for the
        % Scatterers' positions in each dimension as a 3 x 2 array.
        %
        % See also SCATTERERS/XB SCATTERERS/YB SCATTERERS/ZB
        bounds (3,2) double % 3D bounds of the scatterers
    end
    
    % constructor and get/set methods
    methods
        function self = Scatterers(kwargs)
            % SCATTERERS - Construct a Scatterers object
            %
            % scat = SCATTERERS('pos', pos) constructs a Scatterers
            % object scat with scatterers located at the positions pos. The
            % amplitudes are initialized to 1.
            %
            % scat = SCATTERERS('amp', amp) constructs a Scatterers
            % object scat with scatterers having the amplitudes amp. The 
            % positions of all scatterers are initialized to [0;0;0].
            %
            % scat = SCATTERERS('pos', pos, 'amp', amp) constructs a 
            % Scatterers object scat with scatterers located at the
            % positions pos and having the amplitudes amp.
            %
            % scat = SCATTERERS(..., 'c0', c0) sets the homogeneous sound
            % speed to be c0.
            %
            % Example:
            %
            % % Create a scatterer with 30(mm) depth
            % scat1 = Scatterers('pos', [0;0;30e-3]); 
            %
            % % Create an array of scatterers with different amplitudes
            % scat2 = Scatterers('pos', [0;0;1e-3] .* (10:10:50), 'amp', 1:5);
            %
            % % Create a scatterer in a Medium with a slower sound speed
            % scat3 = Scatterers('pos', [0;0;30e-3], 'c0', 1400);
            %
            % See also MEDIUM

            arguments
                kwargs.pos (3,:) double
                kwargs.amp (1,:) double
                kwargs.c0  (1,1) double
            end
            
            % Name/Value initialization
            for f = string(fieldnames(kwargs))', self.(f) = kwargs.(f); end

            % input data sizing
            [Sa, Sp] = deal(size(self.amp,2), size(self.pos,2)); % == 0 if not set

            % check/default amplitudes / positions
            if      Sp &&  Sa % both non-zero
                if Sp ~= Sa, error('Number of scatterers must match!'); end 
            elseif  Sp && ~Sa % default amplitude
                self.amp = ones([1,Sp]);
            elseif ~Sp &&  Sa % default positions
                self.pos = zeros([3,Sa]);
            else % ~Sp && ~Sa % no info - do nothing                
            end
        end
        
        function self = scale(self, kwargs)
            % SCALE - Scale units
            %
            % scat = SCALE(scat, 'dist', factor) scales the distance of the
            % properties by factor. This can be used to convert from meters
            % to millimeters for example.
            %
            % scat = SCALE(scat, 'time', factor) scales the temporal
            % properties by factor. This can be used to convert from
            % seconds to microseconds and hertz to megahertz.
            %
            % Example:
            %
            % % Create a Scatterers
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', 1500) % m, s
            %
            % % convert from meters to millimeters, seconds to microsecond
            % scat = scale(scat, 'dist', 1e3, 'time', 1e6); % mm, us
            % scat.pos
            % scat.c0
            %
            %

            arguments
                self Scatterers
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            self = copy(self); % copy semantics
            cscale = 1; % speed scaling
            if isfield(kwargs, 'dist') % scale distance (e.g. m -> mm)
                self.pos = kwargs.dist .* self.pos;
                cscale = cscale * kwargs.dist; % scale speed
            end
            if isfield(kwargs, 'time')
                cscale = cscale / kwargs.time; % scale speed
            end
            self.c0 = self.c0 * cscale; % scale speed
        end
    end

    % plot methods
    methods 
        function h = plot(self, varargin, plot_args)
            % PLOT - overload the plot function
            % 
            % PLOT(self) plots the locations of the Scatterers self.
            %
            % PLOT(self, ax) uses the axes ax instead of the current axes.
            % 
            % PLOT(..., Name, Value, ...) passes name-value pair arguments
            % to the built-in plot function so that name value pairs that 
            % are valid for plot are valid here. 
            % 
            % h = PLOT(...) returns the handle to the plot.
            % 
            % Plots only the x-z slice.
            % 
            % See also MEDIUM/IMAGESC
            arguments
                self (1,1) Scatterers
            end 
            arguments(Repeating)
                varargin
            end
            arguments
                plot_args.?matlab.graphics.chart.primitive.Line
            end

            % extract axis and other non-Name/Value pair arguments
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes')
                axs = varargin{1}; varargin(1) = []; 
            else, axs = gca;
            end

            % plot
            plot_args = struct2nvpair(plot_args);
            h = plot(axs, self.pos(1,:), self.pos(3,:), varargin{:}, plot_args{:});
        end        
    end

    % interfaces
    methods
        % get SIMUS medium params
        function p = getSIMUSParam(self)
            p = struct('c', self.c0);
            if ~isnan(self.alpha0), p.attenuation = self.alpha0; end
        end
    end

    % overloaded arithmetic
    methods
        function c = plus(a, b)
            % PLUS - Add two Scatterers
            %
            % out = scat_lhs + scat_rhs adds the Scatterers scat_lhs and 
            % the Scatterers scat_rhs to create a new Scatterers out. The
            % Scatterers must have matching medium properties i.e. c0 and
            % alpha0 must match.
            %
            %
            arguments
                a Scatterers
                b Scatterers
            end
            c = copy(a);
            if prop_match(a, b)
                c.pos = [a.pos, b.pos];
                c.amp = [a.amp, b.amp];
            else
                error('Cannot combine Scatterers whos properties do not match.');
            end

            
        end
    end

    % internal type checking convenience methods
    methods(Hidden)
        function tf = prop_match(self, other)
            % PROP_MATCH - Property matching utility
            %
            % tf = PROP_MATCH(self, other) returns true if the set of 
            % propertiesrequired to combine the objects self and other 
            % match or false if otherwise.
            % 
            % This means that properties that must be identical to
            % meaningfully combine two objects.
            %
            % For Scatterers, this is the c0 and alpha0 properties
            %
            %
            arguments
                self Scatterers
                other Scatterers
            end
            
            tf = true; 
            % require that c0 matches
            tf = tf && isequaln(self.c0, other.c0);

            % require the alpha0 matches
            tf = tf && isequaln(self.alpha0, other.alpha0);
        end
    end

    
    % dependent variables
    methods
        function S = get.numScat(self)
            S = unique([size(self.amp,2), size(self.pos,2)]);
            assert(isscalar(S), 'Cannot infer number of scatterers - check input sizes!');
        end
        function b = get.bounds(self), b = cat(1, self.xb, self.yb, self.zb); end
        function b = get.xb(self), b = [min(self.pos(1,:)), max(self.pos(1,:))]; end
        function b = get.yb(self), b = [min(self.pos(2,:)), max(self.pos(2,:))]; end
        function b = get.zb(self), b = [min(self.pos(3,:)), max(self.pos(3,:))]; end
    end    
end