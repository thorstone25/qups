classdef ScanCartesian < Scan
    properties
        x = 1e-3*linspace(-20,20,128) % image x values
        y = 0                         % image y values
        z = 1e-3*linspace(0,40,128)   % image z values
        order = 'ZXY';                % order of the dimensions for display
    end
    
    % dependent parameters
    properties(Dependent)
        size                % size of the final image
        xb                  % image bounds in x
        yb                  % image bounds in y
        zb                  % image bounds in z
        nPix                 % number of pixels in the imaging grid
    end
    
    properties(Dependent, Hidden)
        nx                  % number of samples in x
        ny                  % number of samples in y
        nz                  % number of samples in z
        res                 % resolution in 3D
        resx                % resolution in x
        resy                % resolution in y
        resz                % resolution in z
    end
    
    % get/set & constructor
    methods
        % constructor
        function self = ScanCartesian(varargin)
            % initialize with name-value pairs
            for i = 1:2:nargin
                self.(varargin{i}) = varargin{i+1};
            end            
        end
        
        % image defs
        function setImagingGrid(self, x, y, z)
            [self.x, self.y, self.z] = deal(x, y, z);
        end
        function setImagingBounds(self, x, y, z)
            [self.xb, self.yb, self.zb] = deal(x, y, z);
        end
        function scan = getUSTBScan(self)
            scan = uff.linear_scan(...
                'x_axis', self.x, ...
                'z_axis', self.z ...
                );
        end
    end
        
    % imaging computations
    methods
        function [X, Y, Z, sz] = getImagingGrid(self)
            ord = self.getPermuteOrder(); % get order of variables
            iord = arrayfun(@(o) find(o == ord), [1,2,3]); % inverse ordering of variables
            grid = {self.x, self.y, self.z}; % get axis
            grid = grid(ord); % reorder
            [grid{:}] = ndgrid(grid{:}); % expand in proper order
            grid = grid(iord); % undo reorder
            [X, Y, Z] = deal(grid{:}); % send to variables
            sz = self.size; % output image size
            assert(all(size(X,1:3) == sz), 'Internal error: size mismatch.') %#ok<CPROP> 
        end
                
        function setImageGridOnTarget(self, target, margin)
            % sets the imaging grid around the boundary of target leaves
            % unchanged the number of points on the grid, so the resolution
            % may change whenever this function is called
            % inputs:
            %   target:     Target object
            %   margin:     a 3 x 2 matrix of x/y/z  min/max bounds for the
            %               imaging grid

            if(nargin < 3 || isempty(margin))
                % margin expansion (m)
                margin = [-3e-3 3e-3; ...
                    0    0;...
                    -2e-3 2e-3;];
            end

            % copy the boundaries
            self.xb = target.xb + margin(1,:);
            self.yb = target.yb + margin(2,:);
            self.zb = target.zb + margin(3,:);

        end
    
        function ord = getPermuteOrder(self)
            ord = arrayfun(@(c) find(c == 'XYZ'), self.order);
        end
    
    end
    
   % dependent methods
    methods
        % image sizing
        function n = get.nx(self), n = numel(self.x); end
        function n = get.ny(self), n = numel(self.y); end
        function n = get.nz(self), n = numel(self.z); end
        function n = get.nPix(self), n = self.nx * self.ny * self.nz; end
        function sz = get.size(self),
            sz = [self.nx, self.ny, self.nz];
            sz = sz(self.getPermuteOrder());            
        end
        
        % change number of points -> resample linearly, preserve endpoints
        function set.nx(self, n), self.x = linspace(min(self.x), max(self.x), n); end
        function set.ny(self, n), self.y = linspace(min(self.y), max(self.y), n); end
        function set.nz(self, n), self.z = linspace(min(self.z), max(self.z), n); end

        % get boundaries
        function b = get.xb(self), b = [min(self.x), max(self.x)]; end
        function b = get.yb(self), b = [min(self.y), max(self.y)]; end
        function b = get.zb(self), b = [min(self.z), max(self.z)]; end

        % change boundaries -> resample linearly, preserve number of points
        function set.xb(self, b), self.x = linspace(min(b), max(b), self.nx); end
        function set.yb(self, b), self.y = linspace(min(b), max(b), self.ny); end
        function set.zb(self, b), self.z = linspace(min(b), max(b), self.nz); end
        function xlim(self, xb), self.xb = xb; end
        function ylim(self, yb), self.yb = yb; end
        function zlim(self, zb), self.zb = zb; end


        % get resolution
        function r = get.resx(self), r = diff(self.xb) / (self.nx - 1); end
        function r = get.resy(self), r = diff(self.yb) / (self.ny - 1); end
        function r = get.resz(self), r = diff(self.zb) / (self.nz - 1); end
        function r = get.res(self), r = cat(1, self.resx, self.resy, self.resz); end

        % set resolution -> change number of points
        function set.resx(self,r)
            if ~isnan(r), self.nx = ceil(diff(self.xb) / r) + 1; end
        end
        function set.resy(self,r)
            if ~isnan(r), self.ny = ceil(diff(self.yb) / r) + 1; end
        end
        function set.resz(self,r)
            if ~isnan(r), self.nz = ceil(diff(self.zb) / r) + 1; end
        end
        function set.res(self,r), [self.resx, self.resy, self.resz] = deal(r(1), r(2), r(3)); end
    end
    
end

