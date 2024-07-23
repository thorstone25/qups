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
        pos (3,:) double = zeros([3,0]) % positions for scatterers
        % AMP - amplitude of the scatterers
        %
        % AMP defines the amplitude of the reflection from each scatterer.
        % It must be a vector with the same number of scatterers as are
        % defined by the Scatterers/pos property.
        %
        % See also SCATTERERS/POS SCATTERERS/C0
        amp (1,:) double = ones([1,0]) % amplitude for the scatterers
        % C0 - ambient sound speed
        %
        % C0 defines the speed of sound of the homogeneous medium on which
        % the scatterers are defined.
        %
        % See also SCATTERERS/POS SCATTERERS/AMP SCATTERERS/ALPHA0
        c0  (1,1) double = 1540; % sound speed

        % ALPHA0 - ambient attenuation
        %
        % ALPHA0 defines an attenuation constant for simulation routines 
        % that support one. If alpha0 is NaN, no attenuation is used.
        %
        % To follow convention and maintain consistent units of distance
        % and time, this value is in dB/m/Hz by default and supports
        % scaling e.g. for typical soft tissue attenuation
        %
        %                       0.5   dB / cm / MHz  is equivalent to
        % (0.5 / 1e-2 / 1e+6) = 50e-6 dB /  m /  Hz (SI units) or 
        % (0.5 / 1e+1 / 1   ) = 0.05  dB / mm / MHz (mm / us / MHz)
        % 
        % Example:
        % % FieldII required
        % if ~exist('field_init','file')
        %     error("This example requires FieldII");
        % end
        % 
        % % Create a system to simulate
        % us = UltrasoundSystem(); % default
        % us.fs = 10*us.fs; % increase for FieldII discretization
        % [us.scan.dx, us.scan.dz] = deal(us.lambda / 4); % imaging resolution
        % scat = Scatterers.Grid([11 1 11]',2e-3, 1e-3*[0 0 20]); % targets
        % scat_att = copy(scat);
        % scat_att.alpha0 = 3 * 50e-6; % 10x soft tissue attenuation
        % 
        % % Simulate
        % chd = calc_scat_multi(us, [scat, scat_att]); % simulate
        % chd = subD(hilbert(chd, 2*chd.T), 1:chd.T, chd.tdim);
        % b = DAS(us, chd);
        % 
        % % Display the B-mode images
        % figure; 
        % him    = imagesc(us.scan, b(:,:,1), nexttile()); colorbar; 
        % him(2) = imagesc(us.scan, b(:,:,2), nexttile()); colorbar;
        % colormap gray;
        % title(him(1).Parent,   "No Attenuation"); 
        % title(him(2).Parent, "With Attenuation");
        % linkaxes([him.Parent]);
        % linkprop([him.Parent], "CLim");
        % caxis([-60 0] + max(gather(mod2db(b(:)))));
        % 
        % See also SCATTERERS/C0 SCATTERERS/POS SCATTERERS/AMP
        alpha0 (1,1) double = nan; % attenuation
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
        function scat = Scatterers(kwargs)
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
                kwargs.?Scatterers
            end
            
            % Name/Value initialization
            for f = string(fieldnames(kwargs))', scat.(f) = kwargs.(f); end

            % whether data was set
            [Sa, Sp] = deal(isfield(kwargs, 'amp'), isfield(kwargs, 'pos'));

            % data sizing
            [Na, Np] = deal(size(scat.amp,2), size(scat.pos,2));

            % check/default amplitudes / positions
            if      Sp &&  Sa % both set: validate
                if Np ~= Na
                    if     Np == 1, scat.pos = repmat(scat.pos, [1,Na]);
                    elseif Na == 1, scat.amp = repmat(scat.amp, [1,Np]);
                    else,           error('Number of scatterers must match!');
                    end
                end 
            elseif  Sp && ~Sa % pos set: default amplitude
                scat.amp = ones([1,Np]);
            elseif ~Sp &&  Sa % amp set: default positions
                scat.pos = zeros([3,Na]);
            else % ~Sp && ~Sa % none set: default pos and amp
                scat.pos = [0;0;30e-3];
                scat.amp = 1;
            end
        end
        
        function scat = scale(scat, kwargs)
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
            % scat = scale(scat, 'dist', 1e3, 'time', 1e6) % mm, us
            %
            %

            arguments
                scat Scatterers
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            scat = copy(scat); % copy semantics
            cscale = 1; % speed scaling
            if isfield(kwargs, 'dist') % scale distance (e.g. m -> mm)
                scat.pos    = scat.pos    * kwargs.dist;
                cscale      = cscale      * kwargs.dist; % scale speed
            end
            if isfield(kwargs, 'time')
                cscale = cscale / kwargs.time; % scale speed
            end
            scat.c0     = scat.c0     * cscale; % scale speed
            scat.alpha0 = scat.alpha0 / cscale; % scale attenuation
        end
    end

    % conversion
    methods
        function s = obj2struct(scat)
            % OBJ2STRUCT - Convert a QUPS object into a native MATLAB struct
            %
            % scat = OBJ2STRUCT(scat) converts the Scatterers scat and all
            % of it's properties into native MATLAB structs.
            %
            % Example:
            % % Create a Scatterers
            % scat = Scatterers()
            %
            % % convert to a MATLAB struct
            % scat = obj2struct(scat)
            %
            arguments
                scat Scatterers
            end

            W = warning('off', "MATLAB:structOnObject"); % squash warnings
            s = struct(scat); % convert scan
            s.class = class(scat); % append class info
            warning(W); % restore warnings
        end
    end

    % plot methods
    methods 
        function h = plot(scat, varargin, plot_args)
            % PLOT - overload the plot function
            % 
            % PLOT(scat) plots the locations of the Scatterers scat.
            %
            % PLOT(scat, ax) uses the axes ax instead of the current axes.
            % 
            % PLOT(..., Name, Value, ...) passes name-value pair arguments
            % to the built-in plot function so that name value pairs that 
            % are valid for plot are valid here. 
            % 
            % h = PLOT(...) returns the handle to the plot.
            % 
            % Plots only the x-z slice.
            % 
            % Example:
            % scat = Scatterers('pos', (-5 : 5) .* [1 0 0]' + [0 0 30]');
            % figure;
            % plot(scat, 'r.');
            % 
            % See also MEDIUM/IMAGESC
            arguments
                scat (1,1) Scatterers
            end 
            arguments(Repeating)
                varargin
            end
            arguments
                plot_args.?matlab.graphics.chart.primitive.Line
                plot_args.DisplayName = 'Scatterers'
            end

            % extract axis and other non-Name/Value pair arguments
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes')
                axs = varargin{1}; varargin(1) = []; 
            else, axs = gca;
            end

            % default plotting style
            if isempty(varargin), varargin{1} = '.'; end

            % plot
            plot_args = struct2nvpair(plot_args);
            if any(scat.pos(2,:))
                h = plot3(axs, scat.pos(1,:), scat.pos(3,:), scat.pos(2,:), varargin{:}, plot_args{:});
            else
                h = plot( axs, scat.pos(1,:),                scat.pos(3,:), varargin{:}, plot_args{:});
            end
        end        
    end

    % interfaces
    methods
        % get SIMUS medium params
        function p = getSIMUSParam(scat)
            % GETSIMUSPARAM - Create a MUST compatible parameter struct
            %
            % p = GETSIMUSPARAM(scat) returns a structure with properties
            % for a call to simus().
            %
            % See also ULTRASOUNDSYSTEM/SIMUS 

            p = struct('c', scat.c0);
            if ~isnan(scat.alpha0), p.attenuation = (1e6*1e-2) * scat.alpha0; end % convert to dB/cm/MHz
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
            % Example:
            % left_scat  = Scatterers('pos', [-5 0 30]');
            % right_scat = Scatterers('pos', [+5 0 30]');
            % both_scats = left_scat + right_scat;
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
        function tf = prop_match(this, that)
            % PROP_MATCH - Property matching utility
            %
            % tf = PROP_MATCH(this, that) returns true if the set of 
            % properties required to combine the objects this and that 
            % match or false if otherwise.
            % 
            % For Scatterers, the c0 and alpha0 properties must match.
            %
            % Example:
            % a = Scatterers('c0', 1500);
            % b = Scatterers('c0', 1500);
            % c = Scatterers('c0', 1450);
            % d = Scatterers('c0', 1450, 'alpha0', 0);
            % e = Scatterers('c0', 1450, 'alpha0', nan);
            % 
            % scts = [a b c d e];
            % % compute `scts' == scts` with explicit broadcasting
            % repmat(scts',1,5) == repmat(scts,5,1)
            % prop_match(scts', scts)
            % 
            % See also SCATTERERS.PLUS
            arguments
                this Scatterers
                that Scatterers
            end
            
            % explicit broadcast to full size
            dims = 1:max(ndims(this), ndims(that));
            this = repmat(this, max(1,size(that,dims) ./ size(this,dims)));
            that = repmat(that, max(1,size(this,dims) ./ size(that,dims)));
            
            % check properties
            tf = true;
            for f = ["c0", "alpha0"]
                tf = tf & cellfun(@isequaln, {this.(f)}, {that.(f)});
            end
            tf = reshape(tf, size(this));
        end
    end

    methods(Static)
        function scat = Grid(sz, dp, p0, kwargs)
            % GRID - Create a grid of scatterers
            %
            % scat = Scatterers.Grid([C P R]) creates a 3D grid of point
            % scatterers with C columns (x), P pages (y), and R rows (z).
            % 
            % scat = Scatterers.Grid([C P R], [dx dy dz]) uses a spacing of
            % dx, dy and dz in the x, y, and z dimensions. The default is
            % [C P R] ./ 11.
            % 
            % scat = Scatterers.Grid([C P R], [dx dy dz], p0) uses an
            % initial point p0 as the primary point. The default is 
            % [0 0 0]'.
            % 
            % scat = Scatterers.Grid(..., Name, Value) sets other name
            % value pairs for constructing a Scatterers. The 'pos' property
            % will be overridden.
            % 
            % Example:
            % us = UltrasoundSystem();
            % scat = Scatterers.Grid([5 1 5]', 5*us.lambda, [0 0 10e-3]');
            % 
            % figure;
            % plot(us);
            % hold on;
            % plot(scat, '.');
            % 
            % See also: SCATTERERS.DIFFUSE
            arguments
                sz (3,1) = [1 1 1]'
                dp (3,1) = sz ./ 11
                p0 (3,1) = [0 0 0]'
                kwargs.?Scatterers   
            end
            
            ax = arrayfun(@(x0, dx, n) {x0 + dx * ((0 : n-1) - ((n-1)/2))}, p0, dp, sz); % get axes
            [ax{:}] = ndgrid(ax{:}); % make into a grid
            [ax{:}] = dealfun(@(x)x(:), ax{:}); % vectorize
            p = [ax{:}]'; % 3 x S
            kwargs.pos = p; % set positions
            args = namedargs2cell(kwargs); % non-position args
            scat = Scatterers(args{:}); % set
        end

        function scat = Diffuse(grid, N, dB, kwargs)
            % DIFFUSE - Construct Diffuse scatterers
            %
            % scat = Diffuse(grid, N) constructs a Scatterers scat with N
            % point scatterers within the ScanCartesian grid.
            %
            % scat = Diffuse(grid) uses a default value of 
            % N = 5^D * prod(grid.size); where D is the number of
            % non-singular dimensions e.g. D = 2 if y = 0.
            %
            % scat = Diffuse(grid, N, dB) sets the (relative) scattering
            % intensity in dB. The default is 0.
            % 
            % scat = Diffuse(..., 'c0', c0) sets the sound speed as c0.
            % 
            % scat = Diffuse(..., 'alpha0', alpha0) sets the attenuation
            % coefficient.
            % 
            % Example:
            % us = UltrasoundSystem();
            % grid = copy(us.scan); % simulation region == imaging region
            % N = 25 * prod(range([grid.xb; grid.zb],2) ./ us.lambda); % 25 scats / wavelength
            % scat = Scatterers.Diffuse(grid, N),
            % 
            % See also: SCATTERERS.GRID
            arguments
                grid (1,1) ScanCartesian
                N (1,1) {mustBePositive, mustBeInteger} = prod(5 * grid.size(grid.size > 1)); % default to 5x the grid resolution
                dB (1,:) {mustBeReal} = 0 
                kwargs.c0
                kwargs.alpha0
            end
            pb = cat(1, grid.xb, grid.yb, grid.zb); % boundaires of the grid
            [p0, dp] = deal(min(pb,[],2), range(pb, 2)); % offset / range
            kwargs.pos = p0 + dp.*rand([3 N]); % uniform random positions
            kwargs.amp = db2pow(dB) .* abs(randn([1 N])); % gaussian random amplitudes
            args = namedargs2cell(kwargs); % get optional arguments
            scat = Scatterers(args{:});
        end
    
        function scat = Verasonics(Media, kwargs)
            % Verasonics - Create a Scatterers from a Verasonics Media struct
            %
            % scat = Scatterers.Verasonics(Media) imports the point
            % scatterers in the verasonics Media struct.
            %
            % scat = Scatterers.Verasonics(..., 'c0', c0) uses sound speed
            % c0. The default is 1540.
            % 
            % scat = Scatterers.Verasonics(..., 'scale', scale) scales the
            % locations of the point targets by scale. The default is 1.
            % 
            % Example:
            % % Import the Transducer and Scatterers
            % c0   = Resource.Parameters.speedOfSound;
            % xdc  = Transducer.Verasonics(Trans, c0);
            % scat = Scatterers.Verasonics(Media,'c0',c0,'scale',c0/xdc.fc)
            % 
            % % Display
            % figure;
            % plot(xdc);
            % hold on;
            % plot(scat, 'r.');
            % set(gca, 'YDir', 'reverse')
            % legend();
            % 
            % See also UltrasoundSystem.Verasonics
            arguments
                Media struct
                kwargs.c0 (1,1) double {mustBePositive} = 1540
                kwargs.scale (1,1) double {mustBePositive} = 1
            end
            
            % short-circuit on empty
            if isempty(Media)
                scat = reshape(Scatterers.empty, size(Media)); 
                return;
            end
            
            % set default attenuation
            if ~isfield(Media, 'attenuation')
                [Media.attenuation] = deal(nan); 
            end

            % create a Scatterers for each Media
            scat = arrayfun(@(m) Scatterers(...
                "pos", m.MP(:,1:3)'.*kwargs.scale, ...
                "amp",m.MP(:,4)', ...
                "c0", kwargs.c0, ...
                "alpha0", 1e-4*m.attenuation ...
                ), Media);
        end
    end
    
    % dependent variables
    methods
        function S = get.numScat(scat)
            S = unique([size(scat.amp,2), size(scat.pos,2)]);
            assert(isscalar(S), 'Cannot infer number of scatterers - check input sizes!');
        end
        function b = get.bounds(scat), b = cat(1, scat.xb, scat.yb, scat.zb); end
        function b = get.xb(scat), b = [min(scat.pos(1,:)), max(scat.pos(1,:))]; end
        function b = get.yb(scat), b = [min(scat.pos(2,:)), max(scat.pos(2,:))]; end
        function b = get.zb(scat), b = [min(scat.pos(3,:)), max(scat.pos(3,:))]; end
    end    
end