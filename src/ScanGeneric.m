% SCANGENERIC - Defines a generic imaging region by each pixel's Cartesian coordinates
%
% The ScanGeneric class defines a Scan for any image coordinates.
%
% The Scan class stores definitions for the imaging region. A Scan
% provides a method to return the image pixel coordinates as an ND-array of
% up to 3 dimensions as well as the row vectors for each individual
% dimension.
%
% See also SCAN SCANPOLAR SCANCARTESIAN
classdef ScanGeneric < Scan
    properties
        u (1,:) {mustBeVector} = 1e-3*linspace(  0,50,160)   % image u values
        v (1,:) {mustBeVector} = 1e-3*linspace(-20,20,128)   % image v values
        w (1,:) {mustBeVector} = 0                           % image w values
        trans function_handle = @(u,v,w) cat(1, u, v, w);    % (u,v,w) -> [x;y;z] transform
        pos (3,:,:,:) = []; % image coordinates (3 x U x V x W)
        order = 'UVW'; % order of the dimensions
    end
    
    % dependent parameters
    properties(Dependent)
        size                % size of the final image
        nPix                % number of pixels in the imaging grid
        ub                  % image bounds in u
        vb                  % image bounds in v
        wb                  % image bounds in w
        du                  % step size in u if linearly spaced
        dv                  % step size in v if linearly spaced
        dw                  % step size in w if linearly spaced
    end
    
    properties(Dependent, Hidden)
        nu                  % number of samples in u
        nv                  % number of samples in v
        nw                  % number of samples in w
        udim                % u-dimension
        vdim                % v-dimension
        wdim                % w-dimension
    end

    properties
        ulabel (1,1) string = "U"
        vlabel (1,1) string = "V"
        wlabel (1,1) string = "W"
    end
    
    % get/set & constructor
    methods
        % constructor
        function scan = ScanGeneric(varargin)
            % SCANCARTESIAN - Construct a ScanCartesian
            %
            % scan = SCANCARTESIAN(Name,Value,...) constructs a
            % ScanCartesian using name/value pairs.
            %
            % See also SCANCARTESIAN SCANPOLAR

            % initialize with name-value pairs
            for i = 1:2:nargin
                scan.(varargin{i}) = varargin{i+1};
            end

            % initialize positions if not given explicitly
            if ~any(cellfun(@(c) (ischar(c) || isstring(c)) && (c == "pos"), varargin))
                scan.pos = pixelCube(scan);
            end
        end

        function pos = pixelCube(scan)
            % PIXELCUBE - Create the tensor of pixel positions from the
            % axes
            %
            % pos = PIXELCUBE(scan) returns a 4D array of pixels positions
            % corresponding to the scan axes
            %
            % See also SCAN.POS
            arguments
                scan ScanGeneric
            end
            % get the grid axes order w.r.t. {'u','v','w'}
            iord([scan.udim, scan.vdim, scan.wdim]) = 1:3; % inverse ordering of variables

            % get data in order
            axs = sub({scan.u, scan.v, scan.w}, iord, 2);

            % make a full grid
            [axs{:}] = ndgrid(axs{:});

            % undo order and shift dimension up 1
            axs = cellfun(@(x){reshape(x, [1,size(x)])}, axs(iord)); %#ok<CPROP> 

            % apply (u,v,w) -> [x,y,z]' transform
            pos = scan.trans(axs{:});
        end
        
        % image defs
        function setImagingGrid(scan, u, v, w)
            % SETIMAGINGGRID - Set image axes directly
            %
            % SETIMAGINGGRID(scan, u, v, w) sets the image grid row vectors
            % for the ScanCartesian scan in all coordinates.
            %
            % See also SETIMAGINGBOUNDS
            [scan.u, scan.v, scan.w] = deal(u, v, w);
        end
        function setImagingBounds(scan, u, v, w)
            % SETIMAGINGBOUNDS - Set image axes boundaries
            %
            % SETIMAGINGBOUNDS(scan, u, v, w) sets the image grid bounds 
            % for the ScanCartesian scan in all coordinates.
            %
            % See also SETIMAGINGGRID
            [scan.ub, scan.vb, scan.wb] = deal(u, v, w);
        end
        % scaling
        function scan = scale(scan, kwargs)
            arguments
                scan ScanGeneric
                kwargs.dist (1,1) double
            end
            scan = copy(scan);
            if isfield(kwargs, 'dist')
                d = kwargs.dist;
                % scale distance (e.g. m -> mm)
                [scan.u, scan.v, scan.w, scan.pos] = deal(d*scan.u, d*scan.v, d*scan.w, d*scan.pos);
            end
        end
    end

    % USTB interface methods
    methods
        function scan = getUSTBScan(scan)
            error('Not implemented.');
            scan = uff.scan(...
                'u_axis', scan.u, ...
                'w_axis', scan.w ...
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
            error('Not implemented');

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
            kgrid_args(2,:) = num2cell([scan.dw, scan.du, scan.dv]); % step args
            dims = find(scan.size ~= 1, 1, 'last'); % dimension given by last non-singleton dimension
            
            % create a kWaveGrid
            kgrid_args = kgrid_args(:,1:dims);
            kgrid = kWaveGrid(kgrid_args{:});
        end
    end

    % imaging computations
    methods
        function [X, Y, Z, sz] = getImagingGrid(scan, kwargs)
            arguments
                scan ScanGeneric
                kwargs.vector (1,1) logical = false; 
            end

            sz = scan.size;
            X = reshape(sub(scan.pos, 1, 1), sz);
            Y = reshape(sub(scan.pos, 2, 1), sz);
            Z = reshape(sub(scan.pos, 3, 1), sz);
            
            if nargout == 1
                if kwargs.vector
                    X = scan.pos; % return 3 x perm(U x V x W) NDarray
                else
                    X = {X, Y, Z}; % return (1 x 3) cell array
                end
            end
        end                
    end

   % dependent methods
    methods
        % self validation
        function tf = isconsistent(scan)
            % ISCONSISTENT - Check if a ScanGeneric's properties are consistent
            %
            % tf = ISCONSISTENT(scan) returns true if the dimensions of the
            % position property scan.pos matches the size of the axes
            % properties scan.u, scan.v, and scan.w in the corresponding
            % dimensions and returns false otherwise.
            %
            % Example:
            %
            % scan = ScanGeneric('u', 1:4, 'v', 1:5, 'w', 1:6, 'order','UVW');
            % scan.pos = rand([3 4 5 6]); 
            % isconsistent(scan) % returns true
            % 
            % scan.pos = rand([2 4 5 6]); 
            % isconsistent(scan) % returns false
            % 
            % scan.pos = rand([3 5 5 6]); 
            % isconsistent(scan) % returns false
            % 
            % 
            % 
            tf = true ...
                && (~scan.nu || size(scan.pos, 1+scan.udim) == scan.nu) ...
                && (~scan.nv || size(scan.pos, 1+scan.vdim) == scan.nv) ...
                && (~scan.nw || size(scan.pos, 1+scan.wdim) == scan.nw) ...
            ; %#ok<CPROP> 
        end

        % image sizing
        function n = get.nu(scan), n = numel(scan.u); end
        function n = get.nv(scan), n = numel(scan.v); end
        function n = get.nw(scan), n = numel(scan.w); end
        function n = get.nPix(scan), n = scan.nu * scan.nv * scan.nw; end
        function sz = get.size(scan),
            sz = [scan.nu, scan.nv, scan.nw];
            ord([scan.udim, scan.vdim, scan.wdim]) = 1:3;
            sz = sz(ord);
        end
        function set.size(scan, sz)
            sz(numel(sz)+1:3) = 1; % send omitted dimensions to size 1
            iord = arrayfun(@(c) find(c == scan.order), 'UVW');
            [scan.nu, scan.nv, scan.nw] = deal(sz(iord(1)), sz(iord(2)), sz(iord(3)));
        end
        
        % change number of points -> resample linearly, preserve endpoints
        function set.nu(scan, n), scan.u = linspace(min(scan.u), max(scan.u), n); end
        function set.nv(scan, n), scan.v = linspace(min(scan.v), max(scan.v), n); end
        function set.nw(scan, n), scan.w = linspace(min(scan.w), max(scan.w), n); end

        % get boundaries
        function b = get.ub(scan), b = [min(scan.u), max(scan.u)]; end
        function b = get.vb(scan), b = [min(scan.v), max(scan.v)]; end
        function b = get.wb(scan), b = [min(scan.w), max(scan.w)]; end

        % change boundaries -> resample linearly, preserve number of points
        function set.ub(scan, b), scan.u = linspace(min(b), max(b), scan.nu); end
        function set.vb(scan, b), scan.v = linspace(min(b), max(b), scan.nv); end
        function set.wb(scan, b), scan.w = linspace(min(b), max(b), scan.nw); end

        % get step size - Inf for scalar axes, NaN if not regularly spaced
        function d = get.du(scan), d = uniquetol(diff(scan.u)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end
        function d = get.dv(scan), d = uniquetol(diff(scan.v)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end
        function d = get.dw(scan), d = uniquetol(diff(scan.w)); if isempty(d), d = Inf; elseif ~isscalar(d), d = NaN; end, end

        % set step size - preserve/expand the image bounds, but gaurantee
        % spacing - also gaurantee passes through zero
        function set.du(scan, du), if isinf(du), scan.u = 0; else, scan.u = du * (floor(min(scan.u) / du) : ceil(max(scan.u) / du)); end, end
        function set.dv(scan, dv), if isinf(dv), scan.v = 0; else, scan.v = dv * (floor(min(scan.v) / dv) : ceil(max(scan.v) / dv)); end, end
        function set.dw(scan, dw), if isinf(dw), scan.w = 0; else, scan.w = dw * (floor(min(scan.w) / dw) : ceil(max(scan.w) / dw)); end, end
    
        % get the named dimension
        function d = get.udim(scan), d = find(scan.order == 'U'); end
        function d = get.vdim(scan), d = find(scan.order == 'V'); end
        function d = get.wdim(scan), d = find(scan.order == 'W'); end
    end

    % overloads
    methods(Access=protected)
        function sc = copyElement(scan)
            sc = ScanGeneric('u', scan.u, 'v', scan.v, 'w', scan.w, 'order', scan.order, 'trans', scan.trans);
        end
    end
    
end

