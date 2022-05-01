classdef Medium < handle

    properties(Access=public)
        c0 = 1540;          % reference speed of sound (m/s) - defaults to average sound speed in tissue
        rho0 = 1020;        % reference density (??/??) - defaults to some parameter thingy maybe?
        BoA0 = NaN;         % reference non-linearity - use 9 in tissue?
        alpha0 = NaN;       % reference power law absorption factor (dB/cm/MHz) - use 0.5 in tissue?
        alphap0 = 1.01;     % global power law absorption exponent
        pertreg = {};       % regions of alternate properties given in
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
                            % inputs are outside the perturbation region
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
        
        % get sound speed map
        function [c, rho, BoA, alpha, alphap] = getPropertyMap(self, points)
            % points is 3 x N x ...
            assert(size(points,1) == 3)
            
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
                        % from the function handle itself: so we just try 
                        % % up to 5 outputs until we get it right
                        for nfout = 1:5
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
                        ind = ~isnan(c_r) & ~isnan(rho_r) & ~isnan(BoA_r) & ~isnan(alpha_r) & ~isnan(alphap_r);
                        if nfout >= 1, c(ind)      = c_r(ind);      end
                        if nfout >= 2, rho(ind)    = rho_r(ind);    end
                        if nfout >= 3, BoA(ind)    = BoA_r(ind);    end
                        if nfout >= 4, alpha(ind)  = alpha_r(ind);  end
                        if nfout >= 5, alphap(ind) = alphap_r(ind); end
                    end

                    % update the number of arguments modified
                    nfout_max = max(nfout_max, nfout);
                end
            end
            
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
    methods(Static)
        function medium = Sampled(grid, c, rho, BoA, alpha, alphap0, varargin)
            % MEDIUM/SAMPLED - Create a medium from sampled data
            %
            % medium = Medium.SAMPLED(grid, c, rho, BoA, alpha, alphap0, varargin)
            %
            % create a medium from the already sampled medium
            %
            % grid is {x, y, z} grid for sampling
            %
            % See also: MEDIUM/MEDIUM

            % TODO: use a Scan[Cartesian] instead of a grid

            nullfun = @(p) nan(size(sub(p,1,1)));

            if nargin >= 2 && ~isempty(c), 
                cterp = griddedInterpolant(grid, c, 'linear', 'none');
                cfun = @(p) cterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                cfun = nullfun;
            end
            if nargin >= 3 && ~isempty(rho), 
                rterp = griddedInterpolant(grid, rho, 'linear', 'none');
                rfun = @(p) rterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                rfun = nullfun;
            end
            if nargin >= 4 && ~isempty(BoA), 
                bterp = griddedInterpolant(grid, BoA, 'linear', 'none');
                bfun = @(p) bterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                bfun = nullfun;
            end
            if nargin >= 5 && ~isempty(alpha), 
                aterp = griddedInterpolant(grid, alpha, 'linear', 'none');
                afun = @(p) aterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                afun = nullfun;
            end
            if nargin >= 6 && ~isempty(alphap0), 
                apterp = griddedInterpolant(grid, alphap0, 'linear', 'none');
                apfun = @(p) apterp(sub(p,1,1), sub(p,2,1), sub(p,3,1));
            else
                apfun = nullfun;
            end

            % call constructor
            medium = Medium(varargin{:});

            % add perterbation
            medium.pertreg{end+1} = @(p) dealret(p, cfun, rfun, bfun, afun, apfun);

        end    
    end
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
            % h = IMAGESC(..., arg1, arg2, ...) passes following arguments
            % to the imagesc function
            %
            % See also SCAN/IMAGESC, IMAGESC

            % get the imaging grid
            grd = cell(1,3); % to gather X, Y, Z
            [grd{:}] = scan.getImagingGrid();

            % shift dimensions to 1 x ...
            grd = cellfun(@(x) shiftdim(x, -1), grd, 'UniformOutput', false);
            grd = cat(1,grd{:}); % 3 x ...
            
            % compute the properties on the grid
            [c, rho, BoA, alpha, alphap] = self.getPropertyMap(grd);

            % plot the sound speed on the scan (no scaling!)
            % TODO: allow user to toggle properties by linking data or
            % making a GUI or something
            h = imagesc(scan, real(c), varargin{:}); 
        end
    end
end

% deal return
function varargout = dealret(x, varargin)
varargout = cellfun(@(v) v(x), varargin, 'UniformOutput', false);
end
