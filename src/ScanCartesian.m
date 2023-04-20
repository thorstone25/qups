% SCANCARTESIAN - Imaging region with Cartesian coordinates
%
% The ScanCartesian class defines a Scan with Cartesian image coordinates.
%
% See also SCAN SCANPOLAR SCANGENERIC
classdef ScanCartesian < Scan
    properties
        x (1,:) {mustBeVector} = 1e-3*linspace(-20,20,128) % image x values
        y (1,:) {mustBeVector} = 0                         % image y values
        z (1,:) {mustBeVector} = 1e-3*linspace(0,40,128)   % image z values
        order = 'ZXY'; 
    end
    
    % dependent parameters
    properties(Dependent)
        xb                  % image bounds in x
        yb                  % image bounds in y
        zb                  % image bounds in z
        dx                  % step size in x if linearly spaced
        dy                  % step size in y if linearly spaced
        dz                  % step size in z if linearly spaced
    end
    
    properties(Dependent, Hidden)
        nx                  % number of samples in x
        ny                  % number of samples in y
        nz                  % number of samples in z
        xdim                % x-dimension
        ydim                % y-dimension
        zdim                % z-dimension
    end

    properties
        xlabel (1,1) string = "Lateral"     % plot label for the x axis
        ylabel (1,1) string = "Elevation"   % plot label for the y axis
        zlabel (1,1) string = "Axial"       % plot label for the z axis
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
            grid = arrayfun(@(c) {self.(c)}, lower(self.order)); % get axes
            [grid{:}] = ndgrid(grid{:}); % expand in proper order
            [~, ord] = ismember('XYZ', self.order);
            [X, Y, Z] = deal(grid{ord}); % send to variables
            assert(all(size(X,1:3) == self.size), 'Internal error: size mismatch.')
            if nargout == 1
                if kwargs.vector
                    X = cellfun(@(x) {shiftdim(x,-1)}, {X, Y, Z}); X = cat(1, X{:}); % return 3 x perm(X x Y x Z) NDarray
                else
                    X = {X, Y, Z}; % return (1 x 3) cell array 
                end
            end % pack if 1 output requested
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
            % it aligns with the virtual source Sequence seq. 
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

        % get step size - Inf for scalar axes, NaN if not regularly spaced
        function d = get.dx(self), d = uniquetol(diff(self.x)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end
        function d = get.dy(self), d = uniquetol(diff(self.y)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end
        function d = get.dz(self), d = uniquetol(diff(self.z)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end

        % set step size - preserve/expand the image bounds, but gaurantee
        % spacing - also gaurantee passes through zero
        function set.dx(self, dx), if isinf(dx), self.x = 0; else, self.x = dx * (floor(min(self.x) / dx) : ceil(max(self.x) / dx)); end, end
        function set.dy(self, dy), if isinf(dy), self.y = 0; else, self.y = dy * (floor(min(self.y) / dy) : ceil(max(self.y) / dy)); end, end
        function set.dz(self, dz), if isinf(dz), self.z = 0; else, self.z = dz * (floor(min(self.z) / dz) : ceil(max(self.z) / dz)); end, end
    
        % get the named dimension
        function d = get.xdim(self), d = find(self.order == 'X'); end
        function d = get.ydim(self), d = find(self.order == 'Y'); end
        function d = get.zdim(self), d = find(self.order == 'Z'); end
    end

    % overloads
    methods(Access=protected)
        function sc = copyElement(scan)
            sc = ScanCartesian('x', scan.x, 'y', scan.y, 'z', scan.z, 'order', scan.order ...
                ,'xlabel', scan.xlabel, 'ylabel', scan.ylabel, 'zlabel', scan.zlabel ...
                );
        end
    end
    
end

