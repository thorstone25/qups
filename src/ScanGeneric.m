% SCANGENERIC - Imaging region with arbitrary coordinates
%
% The ScanGeneric class defines a Scan with arbitrary image coordinates.
%
% The ScanGeneric class contains 3 arbitrary axes u, v, and w and
% uses the ScanGeneric.trans property to define the coordinate transform
% (u,v,w) -> (x,y,z).
% 
% See also SCAN SCANPOLAR SCANCARTESIAN
classdef ScanGeneric < Scan
    properties
        u (1,:) {mustBeVector} = 1e-3*linspace(  0,50,160)   % image u values
        v (1,:) {mustBeVector} = 1e-3*linspace(-20,20,128)   % image v values
        w (1,:) {mustBeVector} = 0                           % image w values
        % TRANS - Coordinate transform
        % 
        % SCANGENERIC.TRANS is a function handle that accepts 3 ND-arrays
        % u, v, and w each with size 1 x [scan.size] and returns an 
        % ND-array p of size 3 x [scan.size]. The first dimension of the 
        % output p should contain the x/y/z coordinates.
        %
        % Example:
        % % make a new ScanGeneric and figure
        % scan = ScanGeneric('u',1e-3*(0:20),'v',1e-3*(-5:5));
        % 
        % % 1) the trivial transform: (x = v, y = w, z = u)
        % scan.trans = @(u,v,w) [v; w; u];
        % 
        % figure;
        % plot(scan, nexttile(), '.'); view(2); axis equal;
        % 
        % % 2) define a rotation matrix and translation vector 
        % theta = -30;
        % R = [-sind(theta) 0 cosd(theta); 0 1 0; cosd(theta), 0 sind(theta)];
        % p0 = 1e-3*[-2; 0; 5]; % offset in x/y/z coordinates
        % 
        % % rotation _then_ translation transform
        % scan.trans = @(u,v,w) p0 + pagemtimes(R, [v; w; u]);
        % 
        % plot(scan, nexttile(), '.'); view(2); axis equal;
        %
        % % 3) linearly scale the x/z components across the scaling dimension w.
        % % (x = v * (c0 / c), y = c, z = u * (c0 / c))
        % scan.w = 1500 ./ ([1450 : 10 : 1650]); % distance scaling
        % scan.wlabel = "Scaling";
        % scan.trans = @(u,v,w) [v; w; u];
        %   
        % plot(scan, nexttile(), '.'); view(2); axis equal;
        % 
        % See also SCANGENERIC.POS FUNCTION_HANDLE
        trans function_handle {mustBeScalarOrEmpty} = @(u,v,w) cat(1, v, w, u);  % (u,v,w) -> [x;y;z] transform
        order = 'UVW';
    end
    
    % dependent parameters
    properties(Dependent)
        % POS - Coordinate positions
        % 
        % SCANGENERIC.POS is a depedent property that returns the pixel
        % positions in x/y/z. 
        % 
        % scan.pos = p, stores the pixel position ND-array p as a 
        % 3 x [scan.size]) array. It sets the hidden property scan.pos_, 
        % which can optionally be set directly.
        %
        % p = scan.pos, returns the result of applying the transform
        % scan.trans to the u/v/w axes if scan.pos_ is empty. Otherwise it
        % returns the positions stored at scan.pos_.
        %
        % Example:
        % scan = ScanGeneric('u',1e-3*(0:20),'v',1e-3*(-5:5), 'w', 1e-3*(-2:2), 'order', 'UVW');
        % [~, u, v, w] = ndgrid(0, scan.u, scan.v, scan.w); % 1 x U x V x W
        % 
        % % 1) the trivial transform: (x = v, y = w, z = u)
        % scan.pos = [v; w; u]; % 3 x U x V x W
        % 
        % figure;
        % plot(scan, nexttile(), '.'); view(2); axis equal;
        % 
        % % 2) define a rotation matrix and translation vector 
        % theta = -30;
        % R = [-sind(theta) 0 cosd(theta); 0 1 0; cosd(theta), 0 sind(theta)];
        % p0 = 1e-3*[-2; 0; 5]; % offset in x/y/z coordinates
        % 
        % % rotation _then_ translation transform
        % scan.pos = p0 + pagemtimes(R, [v; w; u]);
        % 
        % plot(scan, nexttile(), '.'); view(2); axis equal;
        % 
        % % 3) linearly scale the x/z components across the scaling dimension w.
        % % (x = v * (c0 / c), y = c, z = u * (c0 / c))
        % scan.w = 1500 ./ ([1450 : 10 : 1650]); % distance scaling
        % scan.wlabel = "Scaling";
        % 
        % [~, u, v, w] = ndgrid(0, scan.u, scan.v, scan.w); % 1 x U x V x W
        % scan.pos = [v; w; u];
        % 
        % plot(scan, nexttile(), '.'); view(2); axis equal;
        % 
        % See also SCANGENERIC.TRANS FUNCTION_HANDLE
        pos (3,:,:,:)       % pixel positions (3 x [scan.size])
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

    properties(Hidden)
        pos_  (3,:,:,:)     % pixel positions (3 x [scan.size])
    end

    properties
        ulabel (1,1) string = "U" % plot label for the u axis
        vlabel (1,1) string = "V" % plot label for the v axis
        wlabel (1,1) string = "W" % plot label for the w axis
    end
    
    % get/set & constructor
    methods
        % constructor
        function scan = ScanGeneric(kwargs)
            % SCANGENERIC - Construct a ScanGeneric
            %
            % scan = SCANGENERIC(Name,Value,...) constructs a
            % ScanCartesian using name/value pairs.
            %
            % See also SCANCARTESIAN SCANPOLAR
            arguments
                kwargs = ?ScanGeneric
            end

            % initialize with name-value pairs
            for f = string(fieldnames(kwargs))'
                scan.(f) = kwargs.(f);
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
            axs = arrayfun(@(c) {scan.(c)}, lower(scan.order)); % e.g. {scan.u, scan.v, scan.w}
            
            % make a full grid
            [axs{:}] = ndgrid(axs{:});

            % undo order and shift dimension up 1
            axs = cellfun(@(x){reshape(x, [1,size(x)])}, axs(iord)); 

            % apply (u,v,w) -> [x,y,z]' transform
            if isempty(scan.trans)
                pos = []; 
            else 
                pos = scan.trans(axs{:}); 
            end
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

            if ~isconsistent(scan)
                warning("The ScanGeneric is inconsistent: unexpected errors may occur.");
            end

            sz = scan.size;
            P = scan.pos();
            X = reshape(sub(P, 1, 1), sz);
            Y = reshape(sub(P, 2, 1), sz);
            Z = reshape(sub(P, 3, 1), sz);
            
            if nargout == 1
                if kwargs.vector
                    X = P; % return 3 x perm(U x V x W) NDarray
                else
                    X = {X, Y, Z}; % return (1 x 3) cell array
                end
            end
        end                
    end

   % dependent methods
    methods
        % load/store or compute positions 
        function p = get.pos(scan)
            if isempty(scan.pos_), p = pixelCube(scan); else; p = scan.pos_; end
        end
        function set.pos(scan, p), scan.pos_ = p; end

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
            sc = ScanGeneric('u', scan.u, 'v', scan.v, 'w', scan.w ...
                , 'order', scan.order, 'trans', scan.trans ...
                ,'ulabel', scan.ulabel, 'vlabel', scan.vlabel, 'wlabel', scan.wlabel ...
                );
        end
    end
    
end

