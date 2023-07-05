% ScanPolar - Imaging region with polar coordinates
%
% The ScanPolar class defines a Scan with polar image coordinates.
%
% See also SCAN SCANCARTESIAN SCANGENERIC
classdef ScanPolar < Scan
    properties
        r (1,:) {mustBeVector} = 1e-3*linspace(0,40,161)  % image range values
        a (1,:) {mustBeVector} = linspace(-45,45,91)      % image angle values (deg)
        y (1,:) {mustBeVector} = 0                        % image elevation values
        order = 'RAY';                
        origin (3,1) = [0;0;0];             % center of the coordinate system with respect to cartesian coordinates (m)
    end
    
    % dependent parameters
    properties(Dependent)
        rb                  % image bounds in range
        ab                  % image bounds in angle (deg)
        yb                  % image bounds in elevation
        dr                  % step size in range if linearly spaced
        da                  % step size in angle if linearly spaced
        dy                  % step size in elevation if linearly spaced
    end
    
    properties(Dependent, Hidden)
        nr                  % number of samples in range
        na                  % number of samples in angle
        ny                  % number of samples in elevation
        rdim                % range dimension
        adim                % angle dimension
        ydim                % elevation dimension
    end

    properties
        rlabel (1,1) string = "Range"       % plot label for the r axis
        alabel (1,1) string = "Angle ^o"    % plot label for the a axis
        ylabel (1,1) string = "Elevation"   % plot label for the y axis
    end

    
    % get/set & constructor
    methods
        % constructor
        function scan = ScanPolar(kwargs)
            % SCANPOLAR - Construct a ScanCartesian
            %
            % scan = SCANPOLAR(Name,Value,...) constructs a
            % ScanPolar using name/value pairs.
            %
            % See also SCANCARTESIAN

            % initialize with name-value pairs
            arguments, kwargs.?ScanPolar; end
            for f = string(fieldnames(kwargs))'
                scan.(f) = kwargs.(f);
            end
        end
        
        % scaling
        function self = scale(self, kwargs)
            arguments
                self ScanPolar
                kwargs.dist (1,1) double
            end
            self = copy(self);
            if isfield(kwargs, 'dist')
                w = kwargs.dist;
                % scale distance (e.g. m -> mm)
                [self.y, self.r] = deal(w*self.y, w*self.r);
                [self.origin] = w*self.origin;
            end
        end
    end

    methods
        function uscan = QUPS2USTB(scan)
            uscan = uff.sector_scan(...
                'apex', uff.point('xyz', scan.origin'), ...
                'azimuth_axis', deg2rad(scan.a(:)), ...
                'depth_axis'  , scan.r(:) ...
                );
        end
    end

    methods(Static)
        function scan = UFF(uscan)
            arguments, uscan uff.sector_scan; end
            scan = arrayfun(@(u) ScanPolar( ...
                'origin',     u.apex.xyz, ...
                'a', rad2deg( u.azimuth_axis), ...
                'r',          u.depth_axis ...
                ), uscan ...
            );
        end
    end

    % imaging computations
    methods
        function [X, Y, Z] = getImagingGrid(self, kwargs)
            arguments
                self ScanPolar
                kwargs.vector (1,1) logical = false;
            end
            [R, A, Y] = self.getImagingGridPolar();
            og = self.origin;
            [Z, X, Y] = pol2cart(deg2rad(A), R, Y);
            [X, Y, Z] = deal(X + og(1), Y + og(2), Z + og(3));
            if nargout == 1
                if kwargs.vector
                    X = cellfun(@(x) {shiftdim(x,-1)}, {X, Y, Z}); X = cat(1, X{:}); % return 3 x perm(X x Y x Z) NDarray
                else
                    X = {X, Y, Z}; % return (1 x 3) cell array 
                end
            end % pack if 1 output requested
        end
        function [R, A, Y] = getImagingGridPolar(self)
            % GETIMAGINGGRIDPOLAR - Return ND-arrays of cartesian coordinates
            %
            % [R, A, Y, sz] = GETIMAGINGGRIDPOLAR(self) returns 
            % multidimensional arrays R, A, and Y corresponding to the 
            % polar pixel coordinates of the ScanPolar self and the size of
            % the Scan sz.
            %
            % G = GETIMAGINGGRIDPOLAR(self) returns the R/A/Y arrays in a 1 x 3
            % cell (as in G = {R, A, Y}).
            %
            % Outputs:
            %   R -    r coordinate (m)
            %   A -    a coordinate (deg)
            %   Y -    y coordinate (m)
            %   sz -   size of the R, A, and Y multidimensional arrays
            %
            % See also GETIMAGINGGRID

            grid = arrayfun(@(c) {self.(c)}, lower(self.order)); % get axes
            [grid{:}] = ndgrid(grid{:}); % expand in proper order
            [~, ord] = ismember('RAY', self.order);
            [R, A, Y] = deal(grid{ord}); % send to variables
            assert(all(size(R,1:3) == self.size), 'Internal error: size mismatch.') 
            if nargout == 1, R = {R, A, Y}; end % pack if 1 output requested
        end
                
        function [b_cart, scanC] = scanConvert(self, b_pol, scanC)
            % SCANCONVERT - Move the image onto a ScanCartesian
            %
            % [b_cart, scanC] = SCANCONVERT(self, b_pol) converts the data
            % b_pol on the ScanPolar self to data b_cart on a new cartesian 
            % scan scanC. The new scan encompasses the ScanPolar self.
            % 
            % b_cart = SCANCONVERT(self, b_pol, scanC) uses a given
            % ScanCartesian scanC instead of creating one.
            % 
            % See also SCANCARTESIAN

            % create an output scan if one not given
            if nargin < 3, scanC = ScanCartesian(self); end

            % get the cartesian points for the output scan
            [X, Y, Z] = scanC.getImagingGrid();
            
            % convert to polar
            og = self.origin;
            [A, R, Y] = cart2pol(Z - og(3), X - og(1), Y - og(2));

            % sample the data at these coordinates
            % TODO: handle different data orders / dimensionality via 
            % permution / page functions
            assert(self.order == "RAY", "Data must be in order 'RAY'.");
            [RG, AG, YG] = ndgrid(self.r, deg2rad(self.a), self.y); % use interp2 for compatibility
            b_cart = cellfun(@(b) {interp2(AG, RG, b, A, R, 'linear', nan)}, num2cell(b_pol, [1,2]));
            b_cart = reshape(cat(3, b_cart{:}), [size(A,1:2), size(b_pol, 3:max(3,ndims(b_pol)))]); %#ok<CPROPLC> 
            % b_cart = interp3(RG, AG, YG, b_pol, R, A, Y, 'linear', nan); % use interp3 for compatibility

            % let nans be nans? or make zero?
            % b_cart(isnan(b_cart)) = 0;
        end

        function scan = ScanCartesian(self)
            % SCANCARTESIAN - Create a ScanCartesian object from a ScanPolar
            %
            % scan = SCANCARTESIAN(self) creates a ScanCartesian object
            % encompassing the region covered by the ScanPolar self.
            %
            % See also SCANCARTESIAN/SCANCARTESIAN

            p = self.positions();
            grd = num2cell([min(p,[],2:4), max(p,[],2:4)],2); % 3 x 2 min/max

            % create a new scan with these boundaries
            scan = ScanCartesian('xb', grd{1}, 'yb', grd{2}, 'zb', grd{3});

            % set the spatial resolution to the axial resolution
            if isfinite(self.dr), [scan.dx, scan.dy, scan.dz] = deal(self.dr); end
        end
    end

    % 
    methods(Hidden)
        function scan = scanCartesian(scan)
            warning("QUPS:ScanPolar:deprecatedSyntax","'scanCartesian' is deprecated: please use 'ScanCartesian' instead.");
            scan = ScanCartesian(scan); 
        end
    end
    
   % dependent methods
    methods
        % image sizing
        function n = get.nr(self), n = numel(self.r); end
        function n = get.na(self), n = numel(self.a); end
        function n = get.ny(self), n = numel(self.y); end
        
        % change number of points -> resample linearly, preserve endpoints
        function set.nr(self, n), self.r = linspace(min(self.r), max(self.r), n); end
        function set.na(self, n), self.a = linspace(min(self.a), max(self.a), n); end
        function set.ny(self, n), self.y = linspace(min(self.y), max(self.y), n); end

        % get boundaries
        function b = get.rb(self), b = [min(self.r), max(self.r)]; end
        function b = get.ab(self), b = [min(self.a), max(self.a)]; end
        function b = get.yb(self), b = [min(self.y), max(self.y)]; end

        % change boundaries -> resample linearly, preserve number of points
        function set.rb(self, b), self.r = linspace(min(b), max(b), self.nr); end
        function set.ab(self, b), self.a = linspace(min(b), max(b), self.na); end
        function set.yb(self, b), self.y = linspace(min(b), max(b), self.ny); end

        % get step size - Inf for scalar axes, NaN if not regularly spaced
        function d = get.dr(self), d = uniquetol(diff(self.r)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end
        function d = get.da(self), d = uniquetol(diff(self.a)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end
        function d = get.dy(self), d = uniquetol(diff(self.y)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end

        % set step size - preserve/expand the image bounds, but gaurantee
        % spacing - also gaurantee passes through zero
        function set.dr(self, dr), if isinf(dr), self.r = 0; else, self.r = dr * (floor(min(self.r) / dr) : ceil(max(self.r) / dr)); end, end
        function set.da(self, da), if isinf(da), self.a = 0; else, self.a = da * (floor(min(self.a) / da) : ceil(max(self.a) / da)); end, end
        function set.dy(self, dy), if isinf(dy), self.y = 0; else, self.y = dy * (floor(min(self.y) / dy) : ceil(max(self.y) / dy)); end, end

        % get the named dimension
        function d = get.rdim(self), d = find(self.order == 'R'); end
        function d = get.adim(self), d = find(self.order == 'A'); end
        function d = get.ydim(self), d = find(self.order == 'Y'); end
    end

    % overloads
    methods(Access=protected)
        function sc = copyElement(scan)
            sc = ScanPolar('a', scan.a, 'y', scan.y, 'r', scan.r ...
                ,'alabel', scan.alabel, 'ylabel', scan.ylabel, 'rlabel', scan.rlabel ...
                , 'order', scan.order, 'origin', scan.origin ...
                );
        end
    end
    
end

