% MEDIUM - Heterogenous medium definition class
%
% The Medium class stores definitions and provides convenience methods
% for representing a distribution of sound speed and density values 
% representing different materials.
%
% The Medium definition can be used to simulate ChannelData from an
% UltrasoundSystem using a finite difference simulation call, such as
% the UltrasoundSystem/kspaceFirstOrder method.
%
% See also SCATTERERS ULTRASOUNDSYSTEM/KSPACEFIRSTORDER CHANNELDATA

classdef Medium < matlab.mixin.Copyable
    properties
        % C0 - ambient sound speed
        c0 (1,1) double {mustBePositive} = 1540; % ambient sound speed
        % RHO0 - ambient density
        rho0 (1,1) double {mustBePositive} = 1020; % ambient density
        % BOA0 - ambient non-linearity
        BoA0 (1,1) double = NaN; % ambient non-linearity
        % ALPHA0 - ambient attenuation factor
        alpha0 (1,1) double = NaN; % ambient power law absorption factor (dB/cm/MHz)
        % ALPHAP0 - global attenuation power 
        alphap0 (1,1) double = 1.01;     % global power law absorption exponent
        % PERTREG - perturbation regions
        % 
        % PERTREG defines regions of alternate properties. A region can be 
        % an element containing a tuple {fun, [c, [rho, [BoA, [alpha]]]}
        % where fun is a filtering function that accepts a
        % multi-dimensional array of points with x/y/z in the first
        % dimension and returns a logicial array of the same size in all
        % dimensions >= 2. All points for which fun is true will have the
        % properties in the 1 x {1,2,3,4} array of properties. 
        % 
        % A region can also be an element containing a function fun that 
        % accepts a multi-dimensional array of points with x/y/z in the 
        % first dimension and returns the corresponding sound speed,
        % density, non-linearity, and attenuation values for each point as
        % a multi-dimensional array of the same size in all dimensions >= 2.
        % For invalid points or properties, it must return NaN.
        %
        % If two or more regions overlap, the last region defined
        % determines the final value at that pixel.
        %
        % Example:
        %
        % reg1 = @(p) deal(...
        % 1500 + 100 .* (sub(p,3,1) ./ 50e-3), ...
        % 1000 - 100 .* (sub(p,3,1) ./ 50e-3) ...
        % ); % make the medium change gradually with depth
        % 
        % reg2 = {@(p)sub(p,3,1) > 50e-3, [1600, 900]}; % for all points z > 50e-3
        %
        % scan = ScanCartesian(...
        % 'z', linspace(0, 60e-3, 241), ...
        % 'x',linspace(-20e-3, 20e-3, 161) ...
        % );
        %
        % med = Medium('pertreg', {reg1, reg2});
        % figure;
        % imagesc(med, scan, ["c", "rho"]);
        pertreg (1,:) cell = cell.empty([1,0]); % perturbation regions
    end
    
    methods
        % constructor
        function self = Medium(kwargs)
            % MEDIUM - MEDIUM constructor
            %
            % med = MEDIUM(...,'c0', c0) constructs a medium with an
            % ambient sound speed of c0.
            %
            % med = MEDIUM(...,'rho0', rho0) constructs a medium with an
            % ambient density of rho0.
            %
            % med = MEDIUM(...,'BoA0', BoA0) constructs a medium with an
            % ambient non-linearity parameter of BoA0. If BoA0 is NaN, it
            % is unused when calling a simulation routine.
            %
            % med = MEDIUM(...,'alpha0', alpha0) constructs a medium with
            % an ambient attenuation parameter of alpha0. If alpha0 is NaN,
            % the attenuation parameters are unused when calling a 
            % simulation routine.
            %
            % med = MEDIUM(...,'alphap0', alphap0) constructs a medium with
            % an ambient attenuation power of alphap0. If alphap0 is NaN,
            % the attenuation parameters are unused when calling a 
            % simulation routine.
            %
            % med = MEDIUM(..., 'pertreg', reg) additionally defines the
            % perturbation functions of the regions where the Medium
            % varies. A perturbation region can be used to define a
            % distribution such that it can be sampled arbitrarily.
            % 
            % See also MEDIUM.SAMPLED MEDIUM/PERTREG
            arguments
                kwargs.c0 (1,1) double
                kwargs.rho0 (1,1) double
                kwargs.BoA0 (1,1) double
                kwargs.alpha0 (1,1) double
                kwargs.alphap0 (1,1) double
                kwargs.pertreg (1,:) cell = cell.empty([1,0])
            end
            
            % Name/Value pair initialization
            for f = string(fieldnames(kwargs))'
                self.(f) = kwargs.(f);
            end
        end
        
        function [c, rho, BoA, alpha] = props(self, scan, prop) %#ok<STOUT> 
            % PROPS - Return the properties of the medium
            %
            % c = PROPS(self, scan) returns the sound speed of the Medium
            % defined on the Scan scan.
            % 
            % [c, rho] = PROPS(...) also returns the density of the Medium.
            %
            % [c, rho, BoA] = PROPS(...) also returns the non-linearity
            % paramemter B/A (B over A). If undefined, the values are all
            % NaN.
            %
            % [c, rho, BoA, alpha] = PROPS(...) also returns the
            % attenuation parameter alpha.
            % 
            % [p1, p2, ...] = PROPS(self, scan, PROP) returns the requested
            % properties of PROP only. PROP must be a string array
            % containing the names of the output property variables 
            % (e.g. ["c", "rho"]).
            % 
            % Example:
            % % Define the properties on a Scan
            % scan = ScanCartesian();
            % c = rand(scan.size);
            % rho = rand(scan.size);
            % 
            % % Construct the Medium
            % med = Medium.Sampled(scan, c, rho);
            % isequal(c  , props(med, scan, 'c'  ))
            % isequal(rho, props(med, scan, 'rho'))
            % 
            % See also TARGET ARGN GETPROPERTYMAP

            arguments
                self Medium
                scan Scan
                prop (1,:) string {mustBeMember(prop, ["c", "rho", "BoA", "alpha"])}
            end

            % get the grid points on the scan
            pts = scan.getImagingGrid(); 

            % stack points in the first dimension
            pts = cellfun(@(x) {shiftdim(x, -1)}, pts);
            pts = cat(1, pts{:});

            % get property map for all points (dimensions are shifted on
            % output)
            nms = ["c", "rho", "BoA", "alpha"]; % property names
            prps(1,:) = cellstr(nms);
            [prps{2,:}] = getPropertyMap(self, pts); % values
            prps = struct(prps{:}); %#ok<NASGU> % make a struct for easier output mapping

            % remap outputs based on prop input
            % TODO: there's gotta be a better way to get dynamic variables
            % ... or we can just use varargout
            if nargout == 1 && nargin >= 3
                for i = 1:numel(prop) % map each requested output
                    eval(nms(i) + "=prps.('" + prop(i) + "');"); 
                end
            else % map all outputs directly
                for n = nms, eval(n + "=prps.('" + n + "');"); end
            end
        end
    end

    methods(Access=private)
        % get properties map - private call that's a lot more involved
        function [c, rho, BoA, alpha] = getPropertyMap(self, points)
            
            
            assert(size(points,1) == 3) % points is 3 x N x ...
            
            % preallocate output matrix
            sz = size(points);
            sz(1) = 1; % functions collapse points in dimension 1
            c      = repmat(self.c0,       sz);
            rho    = repmat(self.rho0,     sz);
            BoA    = repmat(self.BoA0,     sz);
            alpha  = repmat(self.alpha0,   sz);
            
            if ~isempty(self.pertreg)
                % check if for any region, the properties should be
                % changed
                nfout_max = 0;
                for reg = 1:numel(self.pertreg)
                    if isa(self.pertreg{reg}, 'cell') && numel(self.pertreg{reg}) == 2 % this is a masked region
                        % get points within region
                        fun = self.pertreg{reg}{1};
                        ind = gather(fun(points));
                        
                        % modify the property
                        nfout = length(self.pertreg{reg}{2});
                        if nfout >= 1, c  (ind)    = self.pertreg{reg}{2}(1); end
                        if nfout >= 2, rho(ind)    = self.pertreg{reg}{2}(2); end
                        if nfout >= 3, BoA(ind)    = self.pertreg{reg}{2}(3); end
                        if nfout >= 4, alpha(ind)  = self.pertreg{reg}{2}(4); end

                    elseif isa(self.pertreg{reg}, 'function_handle') % this is a functional region
                        % get the values corresponding to the input points
                        fun = self.pertreg{reg};

                        % MATLAB does not promise the number of outputs,
                        % nor provide a convenient way of figuring that out
                        % from the function handle itself, so we just try 
                        % 5 or less until we get it right
                        for nfout = 4:-1:1
                            out = cell(1, nfout);
                            try
                                [out{:}] = fun(points);
                                break;
                            catch
                            end
                        end
                        
                        % expand to all 5 outputs, adding empty cells at the end
                        out_ = [out, cell(1, 4-nfout)]; 
                        [out_{nfout+1:end}] = deal(0); % fill empty cells with dummy value
                        
                        % assign to each input
                        [c_r, rho_r, BoA_r, alpha_r] = deal(out_{:});
                        
                        % set the value for valid points
                        ind = cellfun(@isnan, {c_r, rho_r, BoA_r, alpha_r}, 'UniformOutput', false);
                        if nfout >= 1, c(~ind{1})      = c_r(~ind{1});      end
                        if nfout >= 2, rho(~ind{2})    = rho_r(~ind{2});    end
                        if nfout >= 3, BoA(~ind{3})    = BoA_r(~ind{3});    end
                        if nfout >= 4, alpha(~ind{4})  = alpha_r(~ind{4});  end
                    end

                    % update the number of arguments modified
                    nfout_max = max(nfout_max, nfout);
                end
            end
            
            % TODO: sizing should be handled outside of this function
            % restore sizing
            c      = shiftdim(c     , 1);
            rho    = shiftdim(rho   , 1);
            BoA    = shiftdim(BoA   , 1);
            alpha  = shiftdim(alpha , 1);

            % restrict output
            %%% TODO: restrict output to only nfout_max modified values %%%
        end
    end

    % fullwave interface
    methods
        function maps = getFullwaveMap(self, scan)
            % GETFULLWAVEMAP - Get Fullwave compatible map structure
            %
            % maps = getFullwaveMap(self, scan) returns a map sampled on
            % the Scan scan.
            %
            % See also SCANCARTESIAN MEDIUM/PROPS

            arguments
                self Medium
                scan Scan
            end

            % sample all maps on the grid points
            [c, rho, BoA, alpha] = props(self, scan);

            % set the map properties
            eta = 1 + BoA./2;
            maps = struct('cmap', c, 'rmap', rho, 'amap', alpha, 'nmap', eta);
            maps.nmap(isnan(maps.nmap)) = 0; % set invalid non-linearity to 0
            maps.amap(isnan(maps.amap)) = 0; % set invalid attenuation to 0

            % Use 0 for invalid properties in fullwave(?)
            for f = string(fieldnames(maps))', maps.(f) = nan2zero(maps.(f)); end
        end
    end

    % k-Wave interface
    methods
        function kmedium = getMediumKWave(self, scan)
            % GETMEDIUMKWAVE - Get a kWave compatible medium struct
            %
            % kmedium = GETMEDIUMKWAVE(self, scan) creates a kWave
            % compatible struct from the Medium self and the ScanCartesian
            % scan.
            %
            %
            arguments
                self Medium
                scan Scan
            end

            % get properties in original dimensions
            [c, rho, BoA, alpha] = self.props(scan);

            % get k-Wave order
            ord = arrayfun(@(d) find(d == scan.order), 'ZXY'); % place in this order for kWave

            % move to k-Wave dimensions
            [kmedium.sound_speed, kmedium.density, kmedium.BonA, kmedium.alpha_coeff] = ...
                dealfun(@(x)permute(x, ord), c, rho, BoA, alpha);

            % remove higher order terms if the coefficients are all 0s
            if all(isnan(kmedium.alpha_coeff)), kmedium = rmfield(kmedium, "alpha_coeff"); end
            if all(isnan(kmedium.BonA)),        kmedium = rmfield(kmedium, "BonA"); end

            % set alpha power if alpha coefficient is set
            if isfield(kmedium, 'alpha_coeff'), kmedium.alpha_power = self.alphap0; end
        end
    end

    % Easy constructor
    methods(Static)
        function medium = Sampled(scan, c, rho, BoA, alpha, kwargs, med_args)
            % SAMPLED - Create a Medium from a sampled distribution
            %
            % medium = MEDIUM.SAMPLED(scan, c, rho) creates a Medium with
            % sound speed distrubtion c and density distribution rho
            % defined on the ScanCartesian scan. c and rho must be
            % multi-dimensional arrays with a size and order matching scan.
            % 
            % medium = MEDIUM.SAMPLED(scan, c, rho, BoA) additionally
            % specifies the non-linearity parameter "B/A" a.k.a. 
            % "B over A" as a multi-dimensional array. If BoA is empty, it
            % is ignored.
            %
            % medium = MEDIUM.SAMPLED(scan, c, rho, BoA, alpha)
            % additionally specifies the attenuation parameter alpha as a
            % multi-dimensional array. If alpha is empty, it is ignored.
            %
            % medium = MEDIUM.SAMPLED(..., 'interp', method) specifies the
            % interpolation method. Must be one of 
            % {"nearest"*, "linear", "spline", "cubic", "makima"}
            % 
            % medium = MEDIUM.SAMPLED(...,Name,Value) forwards following
            % arguments to the Medium constructor.
            %
            % Example:
            %
            % % Setup a scan on which to define the Medium
            % sscan = ScanCartesian(...
            % 'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            % 'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % 
            % % Define the Medium
            % [c0, rho0] = deal(1.5e3, 1e3);
            % [c, rho] = deal(c0*ones(sscan.size), rho0*ones(sscan.size));
            % [Xg, ~, Zg] = sscan.getImagingGrid();
            %
            % % Define isoimpedance layers
            % z0 = rho0 * c0; % ambient impedance
            % [c(Zg > 15e-3), rho(Zg > 15e-3)] = deal(1.4e3, z0/1.4e3); % isoimpedance
            % [c(Zg > 25e-3), rho(Zg > 25e-3)] = deal(1.6e3, z0/1.6e3); % isoimpedance
            % [c(Zg > 35e-3), rho(Zg > 35e-3)] = deal(1.4e3, z0/1.4e3); % isoimpedance
            % [c(Zg > 45e-3), rho(Zg > 45e-3)] = deal(1.5e3, z0/1.5e3); % isoimpedance
            %
            % % Define density scatterers
            % rho(Xg == 0 & Zg == 10e-3) = rho0*2;
            % rho(Xg == 0 & Zg == 20e-3) = rho0*2;
            % rho(Xg == 0 & Zg == 30e-3) = rho0*2;
            % rho(Xg == 0 & Zg == 40e-3) = rho0*2;
            % rho(Xg == 0 & Zg == 50e-3) = rho0*2;
            %
            % % Construct the Medium
            % med = Medium.Sampled(sscan, c, rho, 'c0', c0, 'rho0', rho0);
            % 
            % See also MEDIUM/MEDIUM

            arguments
                scan ScanCartesian
                c {mustBeNumeric}       = double.empty()
                rho {mustBeNumeric}     = double.empty()
                BoA {mustBeNumeric}     = double.empty()
                alpha {mustBeNumeric}   = double.empty()
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "spline", "cubic", "makima"])} = "nearest"
                med_args.?Medium % forward settable properties for the Medium class
            end

            grid = {scan.x, scan.y, scan.z};

            % we'll need to expand the grid if a dimension is singular
            sdims = cellfun(@numel, grid) == 1;
            repfun = @(x) repmat(x, 1+sdims); % 1 replicate if singular

            % we'll need to reorder from scan size to 'XYZ' and expand 
            ord = arrayfun(@(c) find(c ==  scan.order),'XYZ');
            szfun  = @(x) repfun(permute(x, ord));
            
            % expand the grid: ideally the user does this somehow
            grid(sdims) = cellfun(@(x) (x + [-1,1]), grid(sdims), 'UniformOutput', false);

            % no property function
            nullfun = @(p) nan(size(sub(p,1,1)));

            if nargin >= 2 && ~isempty(c),
                cterp = griddedInterpolant(grid, szfun(c), kwargs.interp, 'none');
                cfun = @(p) cterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                cfun = nullfun;
            end
            if nargin >= 3 && ~isempty(rho),
                rterp = griddedInterpolant(grid, szfun(rho), kwargs.interp, 'none');
                rfun = @(p) rterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                rfun = nullfun;
            end
            if nargin >= 4 && ~isempty(BoA),
                bterp = griddedInterpolant(grid, szfun(BoA), kwargs.interp, 'none');
                bfun = @(p) bterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                bfun = nullfun;
            end
            if nargin >= 5 && ~isempty(alpha),
                aterp = griddedInterpolant(grid, szfun(alpha), kwargs.interp, 'none');
                afun = @(p) aterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                afun = nullfun;
            end

            % call constructor
            med_args = struct2nvpair(med_args); 
            medium = Medium(med_args{:});

            % add perterbation
            medium.pertreg = [{@(p) dealret(p, cfun, rfun, bfun, afun)}, medium.pertreg];
        end
    end

    % visualization methods
    methods
        function h = imagesc(self, scan, axs, im_args, kwargs)
            % IMAGESC - Image the Medium (without scatterers)
            %
            % h = IMAGESC(self, scan) plots the Medium on the region
            % defined by the Scan (currently just the sound speed).
            %
            % h = IMAGESC(self, scan, ax) uses the axes ax for plotting
            % instead of the current axes
            %
            % h = IMAGESC(..., 'props', props) plots the properties in the string 
            % array or cell array of characters prop instead of just the 
            % sound speed. The values must each be one of 
            % {'c', 'rho', 'BoA', 'alpha'}.
            %
            % h = IMAGESC(..., 'linkaxs', true) links the x and y axes for 
            % all of the images.
            %
            % h = IMAGESC(..., Name, Value, ...) passes following arguments
            % to the built-in IMAGESC function
            %
            % Example:
            % % Define the properties on a Scan
            % scan = ScanCartesian();
            % c = rand(scan.size);
            % rho = rand(scan.size);
            % 
            % % Construct and image the Medium
            % med = Medium.Sampled(scan, c, rho);
            % figure;
            % imagesc(med, scan, ["c", "rho"], 'linkaxs', true);
            % 
            % See also SCAN/IMAGESC IMAGESC LINKAXES
            arguments
                self Medium
                scan Scan
                axs  (1,:) matlab.graphics.axis.Axes = gca
                im_args.?matlab.graphics.primitive.Image
                kwargs.props (1,:) string {mustBeMember(kwargs.props, ["c", "rho", "BoA", "alpha"])} = "c"
                kwargs.linkaxs (1,1) logical = false
            end
            
            % compute the properties on the grid
            [c, rho, BoA, alpha] = self.props(scan);

            % place properties into a cell array of plot objects
            x = {}; % init
            for p = kwargs.props % for each argument
                switch p
                    case 'c',       x = [x, {c}];
                    case 'rho',     x = [x, {rho}];
                    case 'BoA',     x = [x, {BoA}];
                    case 'alpha',   x = [x, {alpha}];
                end
            end

            N = numel(x); % number of plots to show
            if N > 1 % make a subplot for each property
                axs = arrayfun(@(i) {{subplot(1,N,i)}}, 1:N); % create subplots
            else
                axs = {{axs}}; % place in a nested cell array
            end

            % plot the sound speed on the scan (no scaling!)
            % TODO: allow user to toggle properties by linking data or
            % making a GUI or something
            im_args = struct2nvpair(im_args);
            h = cellfun(@(x, axs) imagesc(scan, real(x), axs{:}, im_args{:}), x, axs);

            if kwargs.linkaxs && N > 1, linkaxes(axs); end
        end
    end
end

% deal return
function varargout = dealret(x, varargin)
% DEALRET - Deal return of each function called on the data x
varargout = cellfun(@(v) v(x), varargin, 'UniformOutput', false);
end
