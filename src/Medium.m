classdef Medium < handle

    properties(Access=public)
        c0 = 1540;          % reference speed of sound (m/s) - defaults to average sound speed in tissue
        rho0 = 1020;        % reference density (??/??) - defaults to some parameter thingy maybe?
        BoA0 = NaN;         % reference non-linearity - use 9 in tissue?
        alpha0 = NaN;       % reference power law absorption factor (dB/cm/MHz) - use 0.5 in tissue?
        alphap0 = 1.01;     % global power law absorption exponent
                            % regions of alternate properties given in
                            % a cell array of perturbation regions. A 
                            % region can be a {fun, [c, rho, BoA, alpha]} 
                            % tuple where fun is a filtering function that 
                            % accepts an nd-array of points with x/y/z in 
                            % the first dimension and returns a logicial 
                            % array of the same size in dimensions 2+ 
                            % or a region can be a {fun} that accepts an
                            % nd-array of points with x/y/z in the first
                            % dimension and returns the corresponding sound
                            % speed and density values etc. or NaN if the 
        pertreg = {};       % inputs are outside the perturbation region
    end
    
    methods
        % constructor
        function self = Medium(varargin)
            for i = 1:2:nargin
                switch lower(varargin{i})
                    case {'soundspeed', 'c0'}
                        self.c0 = varargin{i+1};
                    case {'density', 'rho0'}
                        self.rho0 = varargin{i+1};
                    case {'non-linearity', 'boa0'}
                        self.BoA0 = varargin{i+1};
                    case {'absorption', 'alpha0'}
                        self.alpha0 = varargin{i+1};
                    case {'absorption-power', 'alphap0'}
                        self.alphap0 = varargin{i+1};
                    case {'perturbation-regions', 'pertreg'}
                        self.pertreg = varargin{i+1};
                end
            end
        end
        
        function [c, rho, BoA, alpha, alphap] = props(self, scan, prop)
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
            % [c, rho, BoA, alpha, alphap] = PROPS(...) also returns the
            % global attenuation power law parameter alphap. This usage may
            % be deprecated
            %
            % [p1, p2, ...] = PROPS(self, scan, PROP) returns the requested
            % properties of PROP only. PROP must be a string array
            % containing the names of the output property variables 
            % (e.g. ["c", "rho"]).
            % 
            % See also TARGET ARGN GETPROPERTYMAP

            % parse inputs
            if nargin >= 3, prop = string(prop); end % enforce string 

            % get the grid points on the scan
            pts = scan.getImagingGrid(); 

            % stack points in the first dimension
            pts = cellfun(@(x) {shiftdim(x, -1)}, pts);
            pts = cat(1, pts{:});

            % get property map for all points (dimensions are shifted on
            % output)
            nms = ["c", "rho", "BoA", "alpha", "alphap"]; % property names
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

        % get properties map
        function [c, rho, BoA, alpha, alphap] = getPropertyMap(self, points)
            
            
            assert(size(points,1) == 3) % points is 3 x N x ...
            
            % preallocate output matrix
            sz = size(points);
            sz(1) = 1; % functions collapse points in dimension 1
            c      = repmat(self.c0,       sz);
            rho    = repmat(self.rho0,     sz);
            BoA    = repmat(self.BoA0,     sz);
            alpha  = repmat(self.alpha0,   sz);
            alphap = repmat(self.alphap0,  sz);
            
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
                        if nfout >= 5, alphap(ind) = self.pertreg{reg}{2}(5); end

                    elseif isa(self.pertreg{reg}, 'function_handle') % this is a functional region
                        % get the values corresponding to the input points
                        fun = self.pertreg{reg};

                        % MATLAB does not promise the number of outputs,
                        % nor provide a convenient way of figuring that out
                        % from the function handle itself, so we just try 
                        % 5 or less until we get it right
                        for nfout = 5:-1:1
                            out = cell(1, nfout);
                            try
                                [out{:}] = fun(points);
                                break;
                            catch
                            end
                        end
                        
                        % expand to all 5 outputs, adding empty cells at the end
                        out_ = [out, cell(1, 5-nfout)]; 
                        [out_{nfout+1:end}] = deal(0); % fill empty cells with dummy value
                        
                        % assign to each input
                        [c_r, rho_r, BoA_r, alpha_r, alphap_r] = deal(out_{:});
                        
                        % set the value for valid points
                        ind = cellfun(@isnan, {c_r, rho_r, BoA_r, alpha_r, alphap_r}, 'UniformOutput', false);
                        if nfout >= 1, c(~ind{1})      = c_r(~ind{1});      end
                        if nfout >= 2, rho(~ind{2})    = rho_r(~ind{2});    end
                        if nfout >= 3, BoA(~ind{3})    = BoA_r(~ind{3});    end
                        if nfout >= 4, alpha(~ind{4})  = alpha_r(~ind{4});  end
                        if nfout >= 5, alphap(~ind{5}) = alphap_r(~ind{5}); end
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
            alphap = shiftdim(alphap, 1);

            % restrict output
            %%% TODO: restrict output to only nfout_max modified values %%%
        end
    
    end

    % fullwave interface
    methods
        % get Fullwave compatible map struct
        function maps = getFullwaveMap(self, scan)
            % GETFULLWAVEMAP - Get Fullwave compatible map structure
            %
            % maps = getFullwaveMap(self, scan) returns a map sampled on
            % the Scan scan.
            %
            % See also SCANCARTESIAN MEDIUM/PROPS

            % sample all maps on the grid points
            [c, rho, BoA, alpha, ~] = props(self, scan);

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
        function medium = Sampled(scan, c, rho, BoA, alpha, alphap0, varargin)
            % SAMPLED - Create a medium from an array
            %
            % medium = Medium.SAMPLED(scan, c, rho, BoA, alpha, alphap0)
            % creates a Medium with the properties defined by the inputs.
            % They must be empty to use the ambient/default parameters.
            %
            % medium = Medium.SAMPLED(...,Name,Value) forwards following
            % arguments to the constructor.
            %
            % See also: MEDIUM/MEDIUM

            % TODO: use a Scan[Cartesian] instead of a grid

            if ~isa(scan, 'ScanCartesian')
                error('Data must be defined on a ScanCartesian');
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
                cterp = griddedInterpolant(grid, szfun(c), 'nearest', 'none');
                cfun = @(p) cterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                cfun = nullfun;
            end
            if nargin >= 3 && ~isempty(rho),
                rterp = griddedInterpolant(grid, szfun(rho), 'nearest', 'none');
                rfun = @(p) rterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                rfun = nullfun;
            end
            if nargin >= 4 && ~isempty(BoA),
                bterp = griddedInterpolant(grid, szfun(BoA), 'nearest', 'none');
                bfun = @(p) bterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                bfun = nullfun;
            end
            if nargin >= 5 && ~isempty(alpha),
                aterp = griddedInterpolant(grid, szfun(alpha), 'nearest', 'none');
                afun = @(p) aterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                afun = nullfun;
            end
            % returns a scalar: cannot have a distributed power law
            if nargin >= 6 && ~isempty(alphap0),
                apfun = @(p) alphap0; 
            else
                apfun = @(p) nullfun(0);
            end

            % call constructor
            medium = Medium(varargin{:});

            % add perterbation
            medium.pertreg{end+1} = @(p) dealret(p, cfun, rfun, bfun, afun, apfun);

        end
    end

    % visualization methods
    methods
        function h = imagesc(self, scan, varargin)
            % IMAGESC - Image the Medium (without scatterers)
            %
            % h = IMAGESC(self, scan) plots the Medium on the region
            % defined by the Scan (currently just the sound speed).
            %
            % h = IMAGESC(self, scan, ax) uses the axes ax for plotting
            % instead of the current axes
            %
            % h = IMAGESC(..., prop, ...) plots the property prop instead
            % of the sound speed. Prop must be one of {'c', 'rho', 'BoA', 
            % 'alpha', 'alphap'}.
            %
            % h = IMAGESC(..., arg1, arg2, ...) passes following arguments
            % to the built-in IMAGESC function
            %
            % See also SCAN/IMAGESC, IMAGESC
            
            % compute the properties on the grid
            [c, rho, BoA, alpha, alphap] = self.props(scan);

            x = {}; % init
            i = 1; % manual iteration
            axs = {}; % no known axes for plotting
            while i <= numel(varargin), % go through each argument
                v = varargin{i}; % extract it
                if isa(v, 'matlab.graphics.axis.Axes')
                    axs = arrayfun(@(v){{v}},v); varargin(i) = []; continue; % we found the axes to plot on
                elseif (ischar(v) || isstring(v) || iscellstr(v)) % look for a keyword
                    v_ = string(v); % convert all to a string for easier parsing
                    for j = 1:numel(v_) % for each argument
                    switch v_(j)
                        case 'c',       x = [x, {c}];
                        case 'rho',     x = [x, {rho}];
                        case 'BoA',     x = [x, {BoA}];
                        case 'alpha',   x = [x, {alpha}];
                        case 'alphap',  x = [x, {alphap}];
                        otherwise, i = i + 1; continue; % not a prop -> move along
                    end
                    end
                    % was a prop: delete it - don't pass to imagesc
                    varargin(i) = []; 
                else % not a prop - move along
                    i = i + 1; 
                end
            end
            
            % if no props found, just show the sound speed
            if isempty(x), x = {c}; end
            N = numel(x); % number of plots to show
            if isempty(axs), % make a subplot for each property
                if N == 1
                    axs = {{}}; % don't pass an argument still
                else
                for i = 1:N, axs{i} = {subplot(1,N,i)}; end % create subplots
                end
            end

            % plot the sound speed on the scan (no scaling!)
            % TODO: allow user to toggle properties by linking data or
            % making a GUI or something

            h = cellfun(@(x, axs) imagesc(scan, real(x), axs{:}, varargin{:}), x, axs);
        end
    end
end

% deal return
function varargout = dealret(x, varargin)
varargout = cellfun(@(v) v(x), varargin, 'UniformOutput', false);
end
