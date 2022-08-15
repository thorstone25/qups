% SCANCARTESIAN - Defines an imaging region with Cartesian coordinates
%
% The ScanCartesian class defines a Scan for Cartesian image coordinates.
%
% The Scan class stores definitions for the imaging region. A Scan
% provides a method to return the image pixel coordinates as an ND-array of
% up to 3 dimensions as well as the row vectors for each individual
% dimension. It also provides convenience methods for defining apodization
% array defined for the Scan.
%
% See also SCAN SCANPOLAR
classdef ScanCartesian < Scan
    properties
        x = 1e-3*linspace(-20,20,128) % image x values (m)
        y = 0                         % image y values (m)
        z = 1e-3*linspace(0,40,128)   % image z values (m)
        order = 'ZXY';                % order of the dimensions for display
    end
    
    % dependent parameters
    properties(Dependent)
        size                % size of the final image
        xb                  % image bounds in x (m)
        yb                  % image bounds in y (m)
        zb                  % image bounds in z (m)
        nPix                % number of pixels in the imaging grid
        dx                  % step size in x if linearly spaced (m)
        dy                  % step size in y if linearly spaced (m)
        dz                  % step size in z if linearly spaced (m)
    end
    
    properties(Dependent, Hidden)
        nx                  % number of samples in x
        ny                  % number of samples in y
        nz                  % number of samples in z
        res                 % resolution in all coordinates
        resx                % resolution in x (m)
        resy                % resolution in y (m)
        resz                % resolution in z (m)
    end
    
    % get/set & constructor
    methods
        % constructor
        function self = ScanCartesian(varargin)
            % SCANCARTESIAN - Construct a ScanCartesian
            %
            % scan = SCANCARTESIAN(Name,Value,...) constructs a
            % ScanCartesian using name/value pairs.
            %
            % See also SCANPOLAR

            % initialize with name-value pairs
            for i = 1:2:nargin
                self.(varargin{i}) = varargin{i+1};
            end            
        end
        
        % image defs
        function setImagingGrid(self, x, y, z)
            % SETIMAGINGGRID - Set image axes directly
            %
            % SETIMAGINGGRID(self, x, y, z) sets the image grid row vectors
            % for the ScanCartesian self in all coordinates.
            %
            % See also SETIMAGINGBOUNDS
            [self.x, self.y, self.z] = deal(x, y, z);
        end
        function setImagingBounds(self, x, y, z)
            % SETIMAGINGBOUNDS - Set image axes boundaries
            %
            % SETIMAGINGBOUNDS(self, x, y, z) sets the image grid bounds 
            % for the ScanCartesian self in all coordinates.
            %
            % See also SETIMAGINGGRID
            [self.xb, self.yb, self.zb] = deal(x, y, z);
        end
        % scaling
        function self = scale(self, kwargs)
            arguments
                self ScanCartesian
                kwargs.dist (1,1) double
            end
            self = copy(self);
            if isfield(kwargs, 'dist')
                w = kwargs.dist;
                % scale distance (e.g. m -> mm)
                [self.x, self.y, self.z] = deal(w*self.x, w*self.y, w*self.z);
            end
        end
    end

    % USTB interface methods
    methods
        function scan = getUSTBScan(self)
            scan = uff.linear_scan(...
                'x_axis', self.x, ...
                'z_axis', self.z ...
                );
        end
    end

    % k-Wave interface methods
    methods
        function [kgrid, Npml] = getkWaveGrid(scan, varargin)
            % GETKWAVEGRID - Create a kWaveGrid object
            %
            % [kgrid, offset] = GETKWAVEGRID(self, medium, varargin) 
            % creates a kWaveGrid kgrid and the offset between the grid and
            % the original coordinates of the Scan
            %
            %

            % defaults
            kwargs.PML = [4 48]; % 1x1 or 1x2 vector of PML size or range bounds
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % choose PML size
            PML_buf = @(n) (argmin(arrayfun(@(a)max(factor(n+2*a)), kwargs.PML(1):kwargs.PML(end))) + kwargs.PML(1) - 1);
            Npml = arrayfun(PML_buf, scan.size); % PML size in each dim

            % translate into additional outputs in kwave format
            ind_map = [3,1,2]; % UltrasoundSystem to k-wave coordinate mapping ('ZXY')
            [~, o] = ismember(scan.order, 'ZXY'); % should be 1st, 2nd, 3rd
            assert(all(o == (1:numel(o))), 'The Scan must have order ''ZXY'''); % QUPS -> k-Wave mapping
            
            % get the grid sizing with k-Wave mapping
            kgrid_args(1,:) = num2cell(scan.size); % size args
            kgrid_args(2,:) = num2cell([scan.dz, scan.dx, scan.dy]); % step args
            dims = find(scan.size ~= 1, 1, 'last'); % dimension given by last non-singleton dimension
            
            % create a kWaveGrid
            kgrid_args = kgrid_args(:,1:dims);
            kgrid = kWaveGrid(kgrid_args{:});
        end
    end



    % imaging computations
    methods
        function [X, Y, Z, sz] = getImagingGrid(self, kwargs)
            arguments
                self ScanCartesian
                kwargs.vector (1,1) logical = false; 
            end
            ord = self.getPermuteOrder(); % get order of variables
            iord = arrayfun(@(o) find(o == ord), [1,2,3]); % inverse ordering of variables
            grid = {self.x, self.y, self.z}; % get axis
            grid = grid(ord); % reorder
            [grid{:}] = ndgrid(grid{:}); % expand in proper order
            grid = grid(iord); % undo reorder
            [X, Y, Z] = deal(grid{:}); % send to variables
            sz = self.size; % output image size
            assert(all(size(X,1:3) == sz), 'Internal error: size mismatch.') %#ok<CPROPLC,CPROP> 
            if nargout == 1
                if kwargs.vector
                    X = cellfun(@(x) {shiftdim(x,-1)}, {X, Y, Z}); X = cat(1, X{:}); % return 3 x perm(X x Y x Z) NDarray
                else
                    X = {X, Y, Z}; % return (1 x 3) cell array 
                end
            end % pack if 1 output requested
        end
                
        function setImageGridOnTarget(self, target, margin)
            % SETIMAGEGRIDONTARGET - Set imaging grid to surround a Target
            %
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
            % GETPERMUTEORDER - Get the permutation of 'XYZ'
            %
            % ord = GETPERMUTEORDER(self) returns the permutation of the
            % scan.
            % 
            % This functon is likely to be deprecated.
            %
            %

            ord = arrayfun(@(c) find(c == 'XYZ'), self.order);
        end

        function setImageGridOnSequence(self, seq)
            % setImageGridOnSequence Align Scan to a Sequence
            %
            % setImageGridOnSequence(self, seq) modifies the Scan so that
            % it aligns with the Sequence seq. 
            %

            % soft validate the transmit sequence type: it should be focused
            if seq.type ~= "VS", warning(...
                    "Expected sequence type to be VS but instead got " + seq.type + ". This may produce unexpected results."...
                    );
            end

            self.x = sub(seq.focus,1,1); % use these x-values
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
        function set.size(self, sz)
            iord = arrayfun(@(c) find(c == self.order), 'XYZ');
            [self.nx, self.ny, self.nz] = deal(sz(iord(1)), sz(iord(2)), sz(iord(3)));
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
        % set boundaries in x
        function ylim(self, yb), self.yb = yb; end
        % set boundaries in y
        function zlim(self, zb), self.zb = zb; end
        % set boundaries in z

        % get step size - Inf for scalar axes, NaN if not regularly spaced
        function d = get.dx(self), d = uniquetol(diff(self.x)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end
        function d = get.dy(self), d = uniquetol(diff(self.y)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end
        function d = get.dz(self), d = uniquetol(diff(self.z)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end

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

    % overloads
    methods(Access=protected)
        function sc = copyElement(self)
            sc = ScanCartesian('x', self.x, 'y', self.y, 'z', self.z, 'order', self.order);
        end
    end
    
end

