classdef ScanPolar < Scan
    properties
        r = 1e-3*linspace(0,40,128)   % image range values
        a = linspace(-45,45,128)      % image angle values
        y = 0                         % image elevation values
        order = 'RAY';                % order of the dimensions
        origin = [0;0;0];             % center of the coordinate system
    end
    
    % dependent parameters
    properties(Dependent)
        size                % size of the final image
        rb                  % image bounds in range
        ab                  % image bounds in angle
        yb                  % image bounds in elevation
        nPix                 % number of pixels in the imaging grid
    end
    
    properties(Dependent, Hidden)
        nr                  % number of samples in range
        na                  % number of samples in angle
        ny                  % number of samples in elevation
        res                 % resolution in 3D
        resr                % resolution in range
        resa                % resolution in angle
        resy                % resolution in elevation
    end
    
    % get/set & constructor
    methods
        % constructor
        function self = ScanPolar(varargin)
            % initialize with name-value pairs
            for i = 1:2:nargin
                self.(varargin{i}) = varargin{i+1};
            end            
        end
        
        % image defs
        function setImagingGrid(self, r, a, y, origin)
            if nargin >= 5, self.origin = origin; end % set the origin first if provided
            [self.r, self.a, self.y] = deal(r, a, y);
        end
        function scan = getUSTBScan(self)
            warning('Not tested! Please edit the code here.');
            scan = uff.sector_scan(...
                'azimuth_axis', deg2rad(self.a),...
                'depth_axis'  , self.r ...
                );
        end
    end
        
    % imaging computations
    methods
        function [X, Y, Z, sz] = getImagingGrid(self)
            % returns multidimensional arrays corresponding to the
            % positions of the imaging grid and the size of the arrays
            % outputs:
            %   - X:    x coordinate (m)
            %   - Y:    y coordinate (m)
            %   - Z:    x coordinate (m)
            %   - sz:   size of the X, Y, and Z multidimensional arrays

            [R, A, Y, sz] = self.getImagingGridPolar();
            og = self.origin;
            [Z, X, Y] = pol2cart(deg2rad(A), R, Y);
            [X, Y, Z] = deal(X + og(1), Y + og(2), Z + og(3));
        end
        function [R, A, Y, sz] = getImagingGridPolar(self)
            % returns multidimensional arrays corresponding to the
            % positions of the imaging grid and the size of the arrays
            % outputs:
            %   - R:    r coordinate (m)
            %   - A:    a coordinate (m)
            %   - Y:    y coordinate (m)
            %   - sz:   size of the R, A, and Y multidimensional arrays

            ord = self.getPermuteOrder(); % get order of variables
            iord = arrayfun(@(o) find(o == ord), [1,2,3]); % inverse ordering of variables
            grid = {self.r, self.a, self.y}; % get axis
            grid = grid(ord); % reorder
            [grid{:}] = ndgrid(grid{:}); % expand in proper order
            grid = grid(iord); % undo reorder
            [R, A, Y] = deal(grid{:}); % send to variables
            sz = self.size; % output image size
            assert(all(size(R,1:3) == sz), 'Internal error: size mismatch.') %#ok<CPROP> 
        end
                
        function setImageGridOnTarget(self, target, margin)
            % sets the imaging grid around the boundary of target leaves
            % unchanged the number of points on the grid, so the resolution
            % may change whenever this function is called
            % inputs:
            %   target:     Target object
            %   margin:     a 3 x 2 matrix of r/a/y  min/max bounds for the
            %               imaging grid

            if(nargin < 3 || isempty(margin))
                % margin expansion (m)
                margin = [0 3e-3; ...
                          0 0;...
                         -2 2;];
            end

            % convert the target boundaries to their polar equivalent
            og = self.origin;
            [Xb, Yb, Zb] = ndgrid(target.xb, target.yb, target.zb);
            [a_, r_, y_] = cart2pol(Xb(:) - og(1), Yb(:)- og(2), Zb(:) - og(3));

            % copy the boundaries - gaurd against unreasonable values
            self.rb = [min(min(r_), 0   ),     max(r_)      ] + margin(1,:);
            self.ab = [min(min(a_), -180), max(max(a_), 180)] + margin(2,:);
            self.yb = [    min(y_),            max(y_)      ] + margin(3,:);
        end
    
        function ord = getPermuteOrder(self)
            ord = arrayfun(@(c) find(c == 'RAY'), self.order);
        end

        function [b_cart, scanC] = scanConvert(self, b_pol, scanC)
            % SCANCONVERT - Move scan onto a cartesian grid
            %
            % [b_cart, scanC] = SCANCONVERT(self, b_pol) converts the data
            % on the polar scan b_pol to data b_cart on a new cartesian 
            % scan scanC. The new scan encompasses the ScanPolar self
            % 
            % b_cart = SCANCONVERT(self, b_pol, scanC) uses a given
            % ScanCartesian instead of creating one.
            % 

            % create an output scan if one not given
            if nargin < 3, scanC = scanCartesian(self); end

            % get the cartesian points for the output scan
            [X, Y, Z] = scanC.getImagingGrid();
            
            % convert to polar
            og = self.origin;
            [A, R, Y] = cart2pol(Z - og(3), X - og(1), Y - og(2));

            % sample the data at these coordinates
            % TODO: handle different data orders? check data order?
            terp = griddedInterpolant({self.r, self.a}, gather(b_pol), 'cubic', 'none');
            b_cart = terp(R, rad2deg(A));

            % let nans be nans? or make zero?
            % b_cart(isnan(b_cart)) = 0;
        end

        function scan = scanCartesian(self)
            grd = cell(1,3);
            [grd{:}] = self.getImagingGrid(); % get x,y,z values
            grd = cellfun(@(g) {[min(g(:)), max(g(:))]}, grd(:)); % 3 x 2 min/max

            % create a new scan with these boundaries (leave the rest as
            % default)
            scan = ScanCartesian('xb', grd{1}, 'yb', grd{2}, 'zb', grd{3});
        end
    end
    
   % dependent methods
    methods
        % image sizing
        function n = get.nr(self), n = numel(self.r); end
        function n = get.na(self), n = numel(self.a); end
        function n = get.ny(self), n = numel(self.y); end
        function n = get.nPix(self), n = self.na * self.ny * self.nr; end
        function sz = get.size(self),
            sz = [self.nr, self.na, self.ny];
            sz = sz(self.getPermuteOrder());            
        end
        
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
        function rlim(self, rb), self.rb = rb; end
        function alim(self, ab), self.ab = ab; end
        function ylim(self, yb), self.yb = yb; end


        % get resolution
        function r = get.resr(self), r = diff(self.rb) / (self.nr - 1); end
        function r = get.resa(self), r = diff(self.ab) / (self.na - 1); end
        function r = get.resy(self), r = diff(self.yb) / (self.ny - 1); end
        function r = get.res(self), r = cat(1, self.resr, self.resa, self.resy); end

        % set resolution -> change number of points
        function set.resr(self,r)
            if ~isnan(r), self.nr = ceil(diff(self.rb) / r) + 1; end
        end
        function set.resa(self,r)
            if ~isnan(r), self.na = ceil(diff(self.ab) / r) + 1; end
        end
        function set.resy(self,r)
            if ~isnan(r), self.ny = ceil(diff(self.yb) / r) + 1; end
        end
        function set.res(self,r), [self.resr, self.resa, self.resy] = deal(r(1), r(2), r(3)); end
    end

    % overloads
    methods(Access=protected)
        function sc = copyElement(self)
            sc = ScanPolar('a', self.a, 'y', self.y, 'r', self.r, 'order', self.order, 'origin', self.origin);
        end
    end
    
end

