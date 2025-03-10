% MEDIUM - Heterogenous medium definition class
%
% The Medium class stores definitions and provides convenience methods
% for representing a distribution of sound speed and density values 
% representing different materials.
%
% The Medium definition can be used to simulate ChannelData from an
% UltrasoundSystem using a finite difference simulation call, such as
% the UltrasoundSystem.kspaceFirstOrder method.
%
% See also SCATTERERS ULTRASOUNDSYSTEM.KSPACEFIRSTORDER CHANNELDATA

classdef Medium < matlab.mixin.Copyable
    properties
        % C0 - ambient sound speed
        c0 (1,1) double {mustBePositive} = 1540; % ambient sound speed
        % RHO0 - ambient density
        rho0 (1,1) double {mustBePositive} = 1020; % ambient density
        % BOA0 - ambient non-linearity
        BoA0 (1,1) double = NaN; % ambient non-linearity
        % ALPHA0 - ambient attenuation coefficient
        %
        % ALPHA0 defines the attenuation pre-factor coefficient for
        % simulation routines that support one. If alpha0 is NaN, no
        % attenuation is used.
        %
        % To follow convention and maintain consistent units of distance
        % and time, this value is in dB/m/Hz by default and supports
        % scaling e.g. for typical soft tissue attenuation
        %
        %                       0.5   dB / cm / MHz  is equivalent to
        % (0.5 / 1e-2 / 1e+6) = 50e-6 dB /  m /  Hz (SI units) or 
        % (0.5 / 1e+1 / 1   ) = 0.05  dB / mm / MHz (mm / us / MHz)
        % 
        % See also: ALPHA_POWER BOA0
        alpha0 (1,1) double = NaN; % ambient power law absorption pre-factor
        % alpha_power - global attenuation power law (kWave-only)
        alpha_power (1,1) double {mustBeInRange(alpha_power,0,3)} = 1.01; % global power law absorption exponent
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
        % 'x', linspace(-20e-3, 20e-3, 161) ...
        % );
        %
        % med = Medium('pertreg', {reg1, reg2});
        % figure;
        % imagesc(med, scan, "props", ["c", "rho"]);
        pertreg (1,:) cell = cell.empty([1,0]); % perturbation regions
    end
    properties(Hidden)
        pscale (1,1) double = 1 % scale the distances of the points
        cscale (1,1) double = 1 % scale the sound speed in the medium
        rscale (1,1) double = 1 % scale the density of the medium
    end
    methods
        % constructor
        function med = Medium(kwargs)
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
            % med = MEDIUM(...,'alpha_power', alpha_power) constructs a medium with
            % an ambient attenuation power of alpha_power. If alpha_power is NaN,
            % the attenuation parameters are unused when calling a 
            % simulation routine.
            %
            % med = MEDIUM(..., 'pertreg', reg) additionally defines the
            % perturbation functions of the regions where the Medium
            % varies. A perturbation region can be used to define a
            % distribution such that it can be sampled arbitrarily.
            % 
            % See also MEDIUM.SAMPLED MEDIUM.PERTREG
            arguments
                kwargs.c0 (1,1) double
                kwargs.rho0 (1,1) double
                kwargs.BoA0 (1,1) double
                kwargs.alpha0 (1,1) double
                kwargs.alpha_power (1,1) double
                kwargs.pertreg (1,:) cell = cell.empty([1,0])
            end
            
            % Name/Value pair initialization
            for f = string(fieldnames(kwargs))'
                med.(f) = kwargs.(f);
            end
        end
        
        function [c, rho, BoA, alpha] = props(med, scan, prop) %#ok<STOUT> 
            % PROPS - Return the properties of the medium
            %
            % c = PROPS(med, scan) returns the sound speed of the Medium
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
            % [p1, p2, ...] = PROPS(med, scan, PROP) returns the requested
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
            % See also MEDIUM.SAMPLED

            arguments
                med Medium
                scan Scan
                prop (1,:) string {mustBeMember(prop, ["c", "rho", "BoA", "alpha"])} = ["c", "rho", "BoA", "alpha"]
            end

            % get the grid points on the scan
            pts = scan.positions(); 

            % get property map for all points (dimensions are shifted on
            % output)
            nms = ["c", "rho", "BoA", "alpha"]; % property names
            prps(1,:) = cellstr(nms);
            [prps{2,:}] = getPropertyMap(med, pts); % values, positions are scaled
            prps = struct(prps{:});  %#ok<NASGU> % make a struct for easier output mapping

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
    
        function med = scale(med, kwargs)
            % SCALE - Scale the units of the Medium
            %
            % med = SCALE(med, 'dist', pscale) scales the values of
            % distance by pscale. This includes both the input units and
            % the output units i.e. the sound speed and density are scaled
            % as well.
            %
            % med = SCALE(med, 'time', tscale) scales the values of
            % time by tscale.
            %
            % med = SCALE(med, 'mass', mscale) scales the values of
            % mass by mscale.
            %
            % Example:
            % % Define a system in meters, seconds
            % scan = ScanCartesian(...
            %   'x', linspace(-20e-3, 20e-3, 1+40*2^3), ...
            %   'z', linspace(-02e-3, 38e-3, 1+40*2^3)  ...
            % );
            % c = rand(scan.size);
            % rho = rand(scan.size);
            % med = Medium.Sampled(scan, c, rho); % defined in meters
            %
            % % Scale the values to millimiters
            % med_mm  = scale(med , 'dist', 1e3);
            % scan_mm = scale(scan, 'dist', 1e3);
            %
            % % Display the Mediums
            % figure;
            % subplot(1,2,1);
            % imagesc(med, scan);
            % subplot(1,2,2);
            % imagesc(med_mm, scan_mm);
            %
            %
            arguments
                med Medium
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
                kwargs.mass (1,1) double
            end
            med = copy(med);
            
            % compute cumulative scaling
            cscale_ = 1; % speed -> dist / time
            rscale_ = 1; % density -> mass / dist^3
            % ascale_ = 1; % attenuation -> time / dist
            
            if isfield(kwargs, 'dist')
                cscale_    = cscale_     *  kwargs.dist;
                rscale_    = rscale_     / (kwargs.dist ^ 3);
                med.pscale = med.pscale ./  kwargs.dist; % pos -> dist
            end
            if isfield(kwargs, 'time')
                cscale_ = cscale_ / kwargs.time;
            end
            if isfield(kwargs, 'mass')
                rscale_ = rscale_ * kwargs.mass;
            end

            % apply scaling
            med.cscale = med.cscale * cscale_;
            med.rscale = med.rscale * rscale_;
            med.c0     = med.c0     * cscale_;
            med.rho0   = med.rho0   * rscale_;
            med.alpha0 = med.alpha0 / cscale_;
        end
    end

    methods(Access=private)
        % get properties map - private call that's a lot more involved
        function [c, rho, BoA, alpha] = getPropertyMap(med, points)
            
            
            assert(size(points,1) == 3) % points is 3 x N x ...

            % preallocate output matrix
            sz = size(points);
            sz(1) = 1; % functions collapse points in dimension 1
            c      = repmat(med.c0,       sz);
            rho    = repmat(med.rho0,     sz);
            BoA    = repmat(med.BoA0,     sz);
            alpha  = repmat(med.alpha0,   sz);
            
            if ~isempty(med.pertreg)
                % check if for any region, the properties should be
                % changed
                nfout_max = 0;
                for reg = 1:numel(med.pertreg)
                    if isa(med.pertreg{reg}, 'cell') && numel(med.pertreg{reg}) == 2 % this is a masked region
                        % get points within region
                        fun = med.pertreg{reg}{1};
                        ind = gather(fun(med.pscale * points));
                        
                        % modify the property
                        nfout = length(med.pertreg{reg}{2});
                        if nfout >= 1, c  (  ind) = med.pertreg{reg}{2}(1) * med.cscale; end
                        if nfout >= 2, rho(  ind) = med.pertreg{reg}{2}(2) * med.rscale; end
                        if nfout >= 3, BoA(  ind) = med.pertreg{reg}{2}(3)              ; end
                        if nfout >= 4, alpha(ind) = med.pertreg{reg}{2}(4) / med.cscale; end

                    elseif isa(med.pertreg{reg}, 'function_handle') % this is a functional region
                        % get the values corresponding to the input points
                        fun = med.pertreg{reg};

                        % MATLAB does not promise the number of outputs,
                        % nor provide a convenient way of figuring that out
                        % from the function handle itself, so we just try 
                        % 5 or less until we get it right
                        for nfout = 4:-1:1
                            out = cell(1, nfout);
                            try
                                [out{:}] = fun(med.pscale * points);
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
                        if nfout >= 1, c(    ~ind{1}) = c_r(    ~ind{1}) * med.cscale; end
                        if nfout >= 2, rho(  ~ind{2}) = rho_r(  ~ind{2}) * med.rscale; end
                        if nfout >= 3, BoA(  ~ind{3}) = BoA_r(  ~ind{3})              ; end
                        if nfout >= 4, alpha(~ind{4}) = alpha_r(~ind{4}) / med.cscale; end
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
        function maps = getFullwaveMap(med, grd)
            % GETFULLWAVEMAP - Get Fullwave compatible map structure
            %
            % maps = getFullwaveMap(med, scan) returns a map sampled on
            % the Scan scan.
            %
            % Example:
            % % Define properties for two-layers
            % [c0, r0] = deal(1500, 1000); % sound speed / density
            % [c1, r1] = deal(1600, 1200); % sound speed / density
            % [a0, b0] = deal(50e-6, 6  ); % attenuation / non-linearity
            % 
            % % Define the region
            % grd = ScanCartesian();
            % z01 = median(grd.z);
            % f   = @(p) sub(p,3,1) > z01; % 2nd region predicate
            % 
            % % Define the Medium
            % med = Medium('c0',1500,'rho0',r0,'alpha0',a0,'BoA0',b0);
            % med.pertreg{1} = {f, [c1, r1]};
            % 
            % % Display the Medium
            % figure;
            % ttls = ["Sound Speed", "Density", "Attenuation", "Non-Linearity"];
            % h = imagesc(med, grd, "props", ["c", "rho", "alpha", "BoA"]);
            % arrayfun(@colorbar, [h.Parent])
            % arrayfun(@title   , [h.Parent], ttls);
            % 
            % % Get a Fullwave compatible map
            % maps = getFullwaveMap(med, grd),
            % 
            % See also MEDIUM.GETMEDIUMKWAVE MEDIUM.PROPS

            arguments
                med Medium
                grd ScanCartesian
            end

            % sample all maps on the grid points
            [c, rho, BoA, alpha] = props(med, grd);

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
        function kmedium = getMediumKWave(med, grd)
            % GETMEDIUMKWAVE - Get a kWave compatible medium struct
            %
            % kmedium = GETMEDIUMKWAVE(med, scan) creates a kWave
            % compatible struct from the Medium med and the ScanCartesian
            % scan.
            %
            % Example:
            % % Define properties for two-layers
            % [c0, r0] = deal(1500, 1000); % sound speed / density
            % [c1, r1] = deal(1600, 1200); % sound speed / density
            % [a0, b0] = deal(50e-6, 6  ); % attenuation / non-linearity
            % 
            % % Define the region
            % grd = ScanCartesian();
            % z01 = median(grd.z);
            % f   = @(p) sub(p,3,1) > z01; % 2nd region predicate
            % 
            % % Define the Medium
            % med = Medium('c0',1500,'rho0',r0,'alpha0',a0,'BoA0',b0);
            % med.pertreg{1} = {f, [c1, r1]};
            % 
            % % Display the Medium
            % figure;
            % ttls = ["Sound Speed", "Density", "Attenuation", "Non-Linearity"];
            % h = imagesc(med, grd, "props", ["c", "rho", "alpha", "BoA"]);
            % arrayfun(@colorbar, [h.Parent])
            % arrayfun(@title   , [h.Parent], ttls);
            % 
            % % Get a k-Wave compatible map
            % kmedium = getMediumKWave(med, grd),
            % 
            % See also ULTRASOUNDSYSTEM.KSPACEFIRSTORDER MEDIUM.PROPS
            arguments
                med Medium
                grd ScanCartesian
            end

            % get properties in original dimensions
            [c, rho, BoA, alpha] = med.props(grd);

            % get k-Wave order
            ord = arrayfun(@(d) find(d == grd.order), 'ZXY'); % place in this order for kWave

            % move to k-Wave dimensions
            [kmedium.sound_speed, kmedium.density, kmedium.BonA, kmedium.alpha_coeff] = ...
                dealfun(@(x)permute(x, ord), c, rho, BoA, (1e6*1e-2)*alpha); % (Hz^y * m)^-1 -> (MHz^y * cm)^-1

            % remove higher order terms if the coefficients are all 0s
            if all(isnan(kmedium.alpha_coeff)), kmedium = rmfield(kmedium, "alpha_coeff"); end
            if all(isnan(kmedium.BonA)),        kmedium = rmfield(kmedium, "BonA"); end

            % set alpha power if alpha coefficient is set
            if isfield(kmedium, 'alpha_coeff')
                kmedium.alpha_power = med.alpha_power;
                if med.alpha_power == 1
                    warning( ...
                        "QUPS:getMediumKWave:noDispersion", ...
                        "Deactivating dispersion for an exponential coefficient of 1." ...
                        );
                    kmedium.alpha_mode = 'no_dispersion';
                end
            end
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
            % See also MEDIUM.MEDIUM MEDIUM.DIFFUSE

            arguments
                scan (1,1) ScanCartesian
                c {mustBeNumeric}       = double.empty()
                rho {mustBeNumeric}     = double.empty()
                BoA {mustBeNumeric}     = double.empty()
                alpha {mustBeNumeric}   = double.empty()
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "spline", "cubic", "makima"])} = "nearest"
                med_args.?Medium % forward settable properties for the Medium class
            end

            % check that the input data is either empty, or has the correct
            % size
            for arg = {c,rho,BoA,alpha}
                assert(isempty(arg{1}) || isequal(size(arg{1},1:3), scan.size),...
                    "The size of all arguments must match the scan size (" ... 
                    + strjoin(scan.size + ",") + "). Given argument of size ("  ...
                    + strjoin(size(arg{1},1:max(3,ndims(arg{1}))) + ",") + ").")
            end

            % get the grid axes for the gridded interpolant - must be in
            % X/Y/Z order for the 3D points
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
    
        function med = Diffuse(us, grd, kwargs)
            % DIFFUSE - Generate a diffuse (density) scattering Medium
            %
            % med = DIFFUSE(us, grd) generates a diffuse scaterring Medium
            % med on the ScanCartesian grd for the UltrasoundSystem us.
            %
            % rho = DIFFUSE(..., 'scat_per_cell', n) sets the number of
            % scatterers per resolution voxel. A rule of thumb minimum is 12. The
            % default is 25. 
            % 
            % Note: The grid must be sampled finely enough to produce
            % sub-wavelength scattering characteristics. 
            %
            % rho = DIFFUSE(..., 'ampdB', amp) scales the
            % intensity of the variation by amp in dB. The default is 0.
            %
            % med = DIFFUSE(..., Name, Value) sets additionally properties
            % for the Medium med via name value pairs.
            %
            % Example:
            % % Get a default UltrasoundSystem and Medium
            % us = UltrasoundSystem();
            % med = Medium();
            %
            % % Get a simulation region
            % grd = copy(us.scan);
            % [grd.dx, grd.dy, grd.dz] = deal(us.lambda / 8); % set grid resolution
            %
            % % Create a diffuse scattering distribution
            % med = Medium.Diffuse(us, grd);
            %
            % % display
            % figure;
            % him = imagesc(med, grd, "props", ["c",     "rho"    ]);
            % arrayfun(@title,    [him.Parent], ["Speed", "Density"]);
            % arrayfun(@colormap, [him.Parent], ["jet",   "parula" ]);
            % 
            % See also Medium.Sampled
            arguments
                us UltrasoundSystem
                grd ScanCartesian
                kwargs.rho0 (1,1) {mustBeFloat} = us.seq.c0 / 1.5 % heuristic
                kwargs.c0 (1,1) {mustBeFloat} = us.seq.c0
                kwargs.?Medium
                kwargs.scat_per_cell (1,1) {mustBeInteger, mustBePositive} = 25;
                kwargs.ampdB (1,1) {mustBeReal} = 1;
            end

            % alias
            [rho0, c0] = deal(kwargs.rho0, kwargs.c0);

            % archive
            %{
            % transmit signal properties
            % ncycles = min(1/us.xdc.bw_frac, us.seq.pulse.duration * us.xdc.fc); % Number of Cycles as Function of Bandwidth
            wv = conv(conv(us.tx.impulse, us.rx.impulse), us.seq.pulse); % total temporal impulse response
            ncycles = wv.duration * us.xdc.fc;
            ncycles = min(ncycles, 1/us.xdc.bw_frac); % Number of Cycles as Function of Bandwidth

            % wavelength
            lambda = c0 / us.xdc.fc;

            % lower bound on the imaging resolution
            % TODO: compute based on either arbitrary Transducer, or frequency alone
            if isa(us.xdc, 'TransducerArray')
                D = us.xdc.aperture_size;
            else
                warning('Expected a linear array - please check this formula.');
                pn = us.xdc.positions;
                D = norm(range(pn, 2)); % heuristic upper bound on aperture size
            end
            H = us.xdc.height;

            % resolution cell size (minimum?)
            res_cell = [ ...
                (lambda / D) * (range(grd.xb) / grd.dx), ...
                (lambda / H) * (range(grd.yb) / grd.dy), ...
                (lambda / 2) * (    ncycles   / grd.dz) ...
                ];

            % total number of scatterers
            S = round(kwargs.scat_per_cell * grd.nPix ./ prod(res_cell(res_cell > 0)));
            
            %}

            % number of voxels per dim
            wvl = range([grd.xb; grd.yb; grd.zb], 2) ./ us.lambda;
            S = ceil(kwargs.scat_per_cell * prod(wvl(wvl > 0)));

            % generate uniform random positions and normally distributed amplitudes
            N   = arrayfun(@(p) grd.("n"+p), lower(grd.order)); % size in each dimension
            ind = arrayfun(@(N) {randi(N, [S,1], 'like', rho0([]))}, N); % position indices
            as  =                rand(    [S,1], 'like', rho0([]));      % amplitude values

            % assign perturbation
            rho = zeros(grd.size, 'like', rho0);
            rho(sub2ind(grd.size, ind{:})) = as;

            % add to base
            rho = rho0 + db2pow(kwargs.ampdB) * rho0 * rho;
            c = c0 * ones(grd.size);

            % sample the Medium on the grid
            args = namedargs2cell(rmfield(kwargs, ["scat_per_cell","ampdB"]));
            med = Medium.Sampled(grd, c, rho, args{:});
        end
    end

    % visualization methods
    methods
        function h = imagesc(med, scan, varargin, im_args, kwargs)
            % IMAGESC - Image the Medium (without scatterers)
            %
            % h = IMAGESC(med, scan) plots the Medium on the region
            % defined by the Scan (currently just the sound speed).
            %
            % h = IMAGESC(med, scan, ax) uses the axes ax for plotting
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
            % imagesc(med, scan, "props", ["c", "rho"], 'linkaxs', true);
            % 
            % See also SCAN/IMAGESC IMAGESC LINKAXES
            arguments
                med (1,1) Medium
                scan (1,1) Scan
            end
            arguments(Repeating)
                varargin
            end
            arguments
                im_args.?matlab.graphics.primitive.Image
                kwargs.props (1,:) string {mustBeMember(kwargs.props, ["c", "rho", "BoA", "alpha"])} = "c"
                kwargs.linkaxs (1,1) logical = false
            end

            % find axs varargs
            if numel(varargin) >= 1 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                axs = varargin{1}(:)'; varargin(1) = [];
            else % default
                axs = gca;
            end

            % compute the properties on the grid
            [c, rho, BoA, alpha] = med.props(scan);

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
            h = cellfun(@(x, axs) imagesc(scan, real(x), axs{:}, varargin{:}, im_args{:}), x, axs);

            if kwargs.linkaxs && N > 1, axs = [axs{:}]; linkaxes([axs{:}]); end
        end
    end
    properties(Hidden, Dependent)
        alphap0
    end
    methods
        function P = get.alphap0(med), warning("QUPS:Medium:syntaxDeprecated","Medium.alphap0 is deprecated. Use the .alpha_power property instead."); P = med.alpha_power; end
        function set.alphap0(med, P ), warning("QUPS:Medium:syntaxDeprecated","Medium.alphap0 is deprecated. Use the .alpha_power property instead."); med.alpha_power = P; end
    end
end

% deal return
function varargout = dealret(x, varargin)
% DEALRET - Deal return of each function called on the data x
varargout = cellfun(@(v) v(x), varargin, 'UniformOutput', false);
end
