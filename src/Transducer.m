% TRANSDUCER - Abstract Transducer class
%
% Superclass for a diagnostic ultrasound transducer. Any TRANSDUCER has
% methods for the positions and orientations of it's elements as well as
% the characteristics of the elements (width/height, impulse response,
% central frequency, etc.). This class offers definitions for common
% transducers and conversion functions between real transducers (mainly
% from Verasonics) and simulation programs (k-Wave, FieldII, MUST).
%
% See also TRANSDUCERARRAY TRANSDUCERCONVEX TRANSDUCERMATRIX

classdef (Abstract) Transducer < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    properties
        fc (1,1) double = 5e6        % center frequency
        bw (1,2) double = [3.5e6 6.5e6] % bandwidth
        width (1,1) double = 1.5e-4 % width of an element
        height (1,1) double = 6e-3   % height of an element
        numel (1,1) double = 128     % number of elements
        offset (3,1) double = [0;0;0] % the offset from the origin
        rot (1,2) double = [0 0] % [azimuth, elevation] rotation after translation (deg)
        impulse Waveform {mustBeScalarOrEmpty} = Waveform.empty() % the impulse response function of the element
    end

    properties
        el_focus (1,1) double = Inf % elevation focal depth
    end

    properties(GetAccess=public, SetAccess=protected, Dependent)
        area            % area of an element
    end

    properties(Dependent)
        bw_frac         % fractional bandwidth
    end
    properties(Dependent, Hidden)
        origin
    end

    % constructor
    methods
        function xdc = Transducer(kwargs)
            % TRANSDUCER - Transducer constructor
            %
            % xdc = TRANSDUCER(Name, Value, ...) constructs a
            % Transducer using name/value pair arguments.
            %
            % A Transducer is an abstract class and cannot be instantiated
            % directly.
            % 
            % See also TRANSDUCERARRAY TRANSDUCERCONVEX

            arguments
                kwargs.?Transducer
            end
            % name-value pair initialization
            % if width set but not height, default to 20x
            if isfield(kwargs,'width') && ~isfield(kwargs, 'height')
                kwargs.height = 20*kwargs.width;
            end

            % set all properties except bw_frac
            for f = setdiff(string(fieldnames(kwargs)), {'bw_frac'})', xdc.(f) = kwargs.(f); end

            % if frequency but not bandwidth or fractional bandwidth is
            % set, default to 60%
            if ~isfield(kwargs, 'bw') && ~isfield(kwargs, 'bw_frac')
                kwargs.bw_frac = 0.6;
            end

            % set fractional bandwidth last
            if isfield(kwargs, 'bw_frac'), xdc.bw_frac = kwargs.bw_frac; end


            % regardless of input, if impulse is empty, initialize it
            if isempty(xdc.impulse)
                xdc.impulse = xdc.xdcImpulse();
            end
        end
    
        function s = obj2struct(xdc)
            % OBJ2STRUCT - Convert a QUPS object into a native MATLAB struct
            %
            % xdc = OBJ2STRUCT(xdc) converts the Transducer xdc and all of 
            % it's properties into native MATLAB structs.
            %
            % Example:
            %
            % % Create a Transducer
            % xdc = TransducerArray()
            %
            % % convert to a MATLAB struct
            % xdc = obj2struct(xdc)
            %
            arguments, xdc Transducer {mustBeScalarOrEmpty}; end
            
            W = warning('off', "MATLAB:structOnObject"); % squash warnings
            s = struct(xdc); % convert xdc
            if ~isempty(s), s.impulse = obj2struct(s.impulse); end % convert impulse
            s.class = class(xdc); % append class info
            warning(W); % restore warnings
        end
    end

    % manipulation
    methods
        % scaling
        function xdc = scale(xdc, kwargs)
            % SCALE - Scale units
            %
            % xdc = SCALE(xdc, 'dist', factor) scales the distance of the
            % properties by factor. This can be used to convert from meters
            % to millimeters for example.
            %
            % xdc = SCALE(xdc, 'time', factor) scales the temporal
            % properties by factor. This can be used to convert from
            % seconds to microseconds and hertz to megahertz.
            %
            % Example:
            %
            % % Create a Transducer
            % xdc = TransducerArray() % m, s, Hz
            %
            % % convert from meters to millimeters, hertz to megahertz
            % xdc = scale(xdc, 'dist', 1e3, 'time', 1e6) % mm, us, MHz
            %

            arguments
                xdc Transducer
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            xdc = copy(xdc);
            if isfield(kwargs, 'dist')
                w = kwargs.dist;
                % scale distance (e.g. m -> mm)
                [xdc.width, xdc.height, xdc.offset, xdc.el_focus] = deal(w*xdc.width, w*xdc.height, w*xdc.offset, w*xdc.el_focus);
            end
            if isfield(kwargs, 'time')
                w = kwargs.time;
                % scale time (e.g. s -> us / Hz -> MHz)
                [xdc.fc, xdc.bw] = deal(xdc.fc/w, xdc.bw/w);
                xdc.impulse = scale(xdc.impulse, 'time', w);
            end
        end
    end

    % transducer specific methods
    methods (Abstract)
        % POSITIONS - Positions of the elements
        %
        % p = POSITIONS(xdc) returns a 3 x N array representing the 
        % positions of the center of each element.
        % 
        % See also ORIENTATIONS TRANSDUCER.PLOT
        p = positions(xdc); 

        % ORIENTATIONS - Orientation of the elements
        %
        % [theta, phi] = ORIENTATIONS(xdc) returns the azimuth and
        % elevation angles theta and phi of the elements of the transducer.
        % The angles are in degrees.
        %
        % [theta, phi, normal, width, height] = ORIENTATIONS(xdc) also
        % returns 3D vectors representing the normal vector, a vector with
        % the direction and magnitude of the width of the element, and a 
        % vector with the direction and mangitude of the height of the
        % element.
        %
        % See also POSITIONS TRANSDUCER.PLOT
        [theta, phi, normal, width, height] = orientations(xdc); % compute the orientations
    end        

    % define methods derived from position and orientation properties
    methods
        function pch = patches(xdc, sub_div)
            % PATCHES - Compute a matrix of sub-elements
            %
            % pch = PATCHES(xdc,sub_div) computes a cell matrix where each
            % element contains {X,Y,Z,C} tuples of 2x2 matrices specifying four
            % corners of each patch, where each patch is a subdivision of the
            % elements into sub-elements specified by the 1 x 2 array sub_div.
            % pch is a [Nel x Ndiv] matrix where Nel is the number of elements
            % and Ndiv is the number of sub-elements given by the
            % prod(sub_div).
            %
            % See also GETBARYCENTERS POSITIONS ORIENTATIONS
            arguments
                xdc Transducer
                sub_div (1,2) {mustBeInteger, mustBePositive} = [1,1]
            end
            % get the matrix of sub element patches
            % get unweighted difference in width and height from element center
            [dx, dy] = ndgrid(...
                ((0 : sub_div(1) - 1) / sub_div(1)) - 0.5, ...
                ((0 : sub_div(2) - 1) / sub_div(2)) - 0.5 ...
                ); % E1 x E2
            [dxs, dys] = deal(1 / sub_div(1), 1 / sub_div(2)); % sub difference

            % reshape to E x 4
            dx = dx(:); dy = dy(:);
            dx = [dx, dx + dxs, dx      , dx + dxs];
            dy = [dy, dy      , dy + dys, dy + dys];

            % move up to 1 x 1 x E x 4
            [dx, dy] = deal(shiftdim(dx, -2), shiftdim(dy, -2));

            % get the difference matrices (3 x N)
            [~, ~, ~, w, h] = xdc.orientations; % difference vectors

            % get the positions
            pn = xdc.positions;

            % get the element dimensions
            wel = xdc.width;
            hel = xdc.height;

            % get patch locations
            pc = pn + wel .* w .* dx + hel .* h .* dy; % 3 x N x E x 4
            pc(4,:) = 1; % x/y/z/c x N x E x 4 - set apodization to 1 
            % TODO: for soft baffle, this should be some distribution
            % instead
            
            % convert to patch structure {{2 x 2} x x/y/z/c} x 1 x N x E
            pc = permute(pc, [5,4,1,2,3]);  % 1 x 4 x x/y/z/c x N x E
            pc = reshape(pc, [2,2,size(pc,3:5)]); % 2 x 2 x x/y/z/c x N x E
            pch = cellfun(@(pc) shiftdim(num2cell(pc, 1:2),1), num2cell(pc, 1:3), 'UniformOutput',false);
            pch = shiftdim(pch, 3);
                        
        end
    
        function pb = bounds(xdc)
            % BOUNDS - Compute the boundaries of the Transducer
            %
            % pb = BOUNDS(xdc) returns a 3 x 2 array of the minimum and maximum
            % cartesian coordinate in x/y/z for the Transducer xdc.
            %
            % See also POSITIONS ORIENTATIONS PATCHES
    
            arguments, xdc Transducer, end

            % transducer patches of {x,y,z,c} bound tuples
            pch = xdc.patches([1,1]);

            % get min/max bounds of the tx by iterating over each patch
            pb = [inf(3,1), -inf(3,1)];
            for i = 1:xdc.numel
                pchi = pch{i}(1:3);
                pb(:,1) = min(pb(:,1), cellfun(@(pch) min(pch, [], 'all'), pchi(:)));
                pb(:,2) = max(pb(:,2), cellfun(@(pch) max(pch, [], 'all'), pchi(:)));
            end
        end
    
        function p = getBaryCenters(xdc, el_sub_div)
            % GETBARYCENTERS - computes barycenters of the sub-elements
            % 
            % p = GETBARYCENTERS(xdc) computes the positions of the 
            % transducer sub-elements by computing the barycenters 
            % (geometric mean) of the corners of the sub-elements computed 
            % xdc.patches and returns them as a 3 x N matrix of positions
            % p.
            %
            % p = GETBARYCENTERS(xdc, sub_div) uses the 2-element
            % array sub_div to divide each element into 
            % sub_div(1) x sub_div(2) sub-elements and returns the 
            % positions as a 3 x N x E array of positions p, where 
            % E = prod(el_sub_div). The default is [1,1].
            %
            % See also PATCHES POSITIONS ORIENTATIONS

            arguments
                xdc Transducer
                el_sub_div (1,2) {mustBeInteger, mustBePositive} = [1,1]
            end

            el_patches = xdc.patches(el_sub_div);

            p = cellfun(@(patch) cellfun(@(comp) mean(comp(:)),...
                patch(1:3)'), el_patches, 'UniformOutput', false);
            p = reshape(cat(1, p{:}), [3, size(p)]);
        end

        function p = transPos(xdc, p, method)
            % transPos - Transform positions
            % 
            % p = transPos(xdc, p) rotates and translates the matrix of
            % positions p according to the rotation xdc.rot and translation
            % xdc.offset. The positions must be a 3 x N array in cartesian
            % coordinates.
            %
            % p = transPos(xdc, p, "quat") uses a quaternionic expression
            % to performs the rotation.
            % 
            % p = transPos(xdc, p, "mat") uses a matrix expression to
            % performs the rotation. The default is "mat". 
            %
            % This method is for internal use, and may be moved or changed
            % in a future version.
            % 
            % See also TRANSDUCER.ROT TRANSDUCER.OFFSET

            arguments
                xdc Transducer
                p (3,:)
                method (1,1) string {mustBeMember(method, ["quat", "mat"])} = "mat" % rotate with quaternions
            end
            switch method
                case "quat"
                q = prod(quaternion([-xdc.rot(2),0,0;0,xdc.rot(1),0], 'rotvecd'));
                p = rotatepoint(q, p')' + xdc.offset;
                case "mat" % rotate with matrices
                [az, el] = deal(xdc.rot(1), xdc.rot(2)); % azimuth, elevation rotation
                Raz = [cosd(az), 0, sind(az); 0 1 0; -sind(az), 0, cosd(az)]; % azimuth   rotation (about y)
                Rel = [1 0 0; 0, cosd(el), sind(el); 0, -sind(el), cosd(el)]; % elevation rotation (about x)
                p = (Rel * Raz) * p + xdc.offset; % rotation and translation
            end
        end
        
        function [pf, nf] = focActive(xdc, apd, r)
            % focActive - Create foci for the active apertures
            %
            % pf = focActive(xdc, apd, r) creates an array of foci pf at a
            % focal depth r from the Transducer xdc with the apodization
            % matrix apd. 
            % 
            % The array apd must be a (N x S) array of weights where 
            % N == xdc.numel and S is the number of transmit pulses. The
            % median of the non-zero elements of each transmit in apd are
            % used generate the beam origins, and the foci are placed at a
            % range r normal to the surface of the transducer from the beam
            % origins. The range r can be a scalar or a (1 x S) array of
            % ranges per transmit pulse.
            % 
            % [pf, nf] = focActive(...) additionaly returns the normal
            % vector from center of the active aperture to the focal
            % points.
            % 
            % Note: a negative value of r will define a diverging wave.
            % 
            % Example:
            % % Create a Transducer
            % xdc = TransducerConvex.C5_2v();
            % 
            % % Create a walking aperture of 64 elements each
            % apd = Sequence.apWalking(xdc.numel, 64);
            % 
            % % Create a focused Sequence at a range of 50mm
            % pf = xdc.focActive(apd, 50e-3);
            % seq = SequenceRadial( ...
            %     'type','FC', 'focus',pf, 'apd',apd, 'apex',xdc.center ...
            % ); 
            % 
            % % Create and plot the system
            % us = UltrasoundSystem('xdc', xdc, 'seq', seq);
            % plot(us);
            % 
            % See also SEQUENCE.APWALKING
            arguments
                xdc (1,1) Transducer
                apd (:,:) {mustBeNumericOrLogical} % apodization (N x S)
                r (1,:) {mustBeReal, mustBeFinite} = 0
            end

            % central element of active apertures
            ic = cellfun(@(a) median(find(a)), num2cell(apd,1)); % 1 x S

            % compute beams based on transducer geometry
            if any(arrayfun(@(s)isa(xdc,s),["TransducerArray", "TransducerMatrix"])) % linear interpolation
                pn = xdc.positions(); % element position (3 x N)
                [~,~,nn] = xdc.orientations(); % element normals (3xN)
                pnc = (pn(:,floor(ic)) + pn(:,ceil(ic))) ./ 2; % mean position
                nf = (nn(:,floor(ic)) + nn(:,ceil(ic))) ./ 2; % mean normal
                pf = pnc + r .* nf; % create focal positions
            elseif isa(xdc,"TransducerConvex") % angular interpolation
                th = xdc.orientations(); % element angles (1xN)
                thc = (th(:,floor(ic)) + th(:,ceil(ic))) ./ 2; % mean angle (1 x S)
                nf = [sind(thc); 0*thc; cosd(thc)]; % normal vectors
                pf = (xdc.radius + r) * nf + xdc.center; % extend radius from beam origins
            else
                error("Sequence:focActive:unsupportedTransducer", ...
                    "A "+class(xdc)+" is not supported.");
            end
        end
    end

    % toolbox conversion functions
    methods
        function aperture = getFieldIIAperture(xdc, sub_div, focus)
            % GETFIELDIIAPERTURE - Create a FieldII aperture object
            %
            % ap = GETFIELDIIAPERTURE(xdc) creates a FieldII aperture ap
            % from the Transducer xdc. If FieldII is not initialized before
            % calling this function, it will be initiated.
            %
            % ap = GETFIELDIIAPERTURE(xdc, sub_div) divides the width and
            % height of each element by into a sub_div(1) x sub_div(2)
            % sub-elements. The default is [1 1].
            %
            % ap = GETFIELDIIAPERTURE(xdc, sub_div, focus) creates an aperture
            % with an elevational focus at the point focus. A focal depth of
            % Inf represents no focus. The default is [0 0 Inf].
            %
            % When creating the aperture, an infinite focus in instead set to
            % realmax('single').
            %
            % See also ULTRASOUNDSYSTEM.CALC_SCAT_ALL, FIELD_INIT
            arguments
                xdc Transducer
                sub_div (1,2) double = [1,1]
                focus (1,3) double = [0 0 realmax('single')]
            end

            focus(isinf(focus)) = realmax('single') .* sign(focus(isinf(focus))); % make focus finite
            sdiv = sub_div; % alias
            aperture = arrayfun(@make_fieldII_aperture, xdc);
            
            function aperture = make_fieldII_aperture(xdc)
                pch = xdc.patches(sdiv); % [Nel x Ndv] array with  {X / Y / Z / C} tuples
                r = zeros([size(pch'),19]); % Ndv x Nel x 19
                for i = 1 : size(pch,1) % each element
                    for j = 1 : size(pch,2) % each subelement
                        pchij = pch{i,j}; % get tuple
                        p = reshape(permute(cat(3, pchij{1:3}), [3,1,2]), 3, 4); % get as 3 x 4 array
                        p = p(:,[1,2,4,3]); % swap 4th<->3rd for clockwise ordering
                        r(j,i,:) = [i, p(:)', 1, [xdc.width, xdc.height] ./ sdiv, mean(p,2)']; % get rectangle
                    end
                end

                % reshape arguments
                r = reshape(r, [numel(pch) 19]); %#ok<CPROPLC> % rectangles: [sdiv x element] x 19
                c = double(gather(xdc.positions()')); % element centers

                % Generate aperture for emission
                try evalc('field_info'); catch, field_init(-1); end
                aperture = xdc_rectangles(r, c, focus);
            end
        end

        function probe = QUPS2USTB(xdc)
            % QUPS2USTB - Create a USTB compatible uff.probe object
            %
            % probe = QUPS2USTB(xdc) creates a uff.probe probe from the
            % Transducer xdc. USTB must be on the path.
            %
            % Example:
            %
            % probe = QUPS2USTB(TransducerArray.L12_3v());
            %
            % See also UFF.PROBE
            arguments, xdc Transducer, end
            probe = arrayfun(@(xdc) uff.probe(...
                'geometry', [ ...
                xdc.positions(); ...
                deg2rad(argn(1, @orientations, xdc)); ...
                deg2rad(argn(1, @orientations, xdc)); ...
                repmat([ ...
                xdc.width; ...
                xdc.height ...
                ], [1, xdc.numel]) ...
                ]', ...
                'origin', uff.point('xyz', xdc.offset(:)') ...
                ), xdc);
        end
    end

    % Verasonics conversion functions
    methods (Static)
        function xdc = Verasonics(Trans, c0)
            % VERASONICS - Construct a Transducer from a Verasonics struct
            %
            % xdc = Transducer.VERASONICS(Trans) constructs a Transducer 
            % from the properties defined in the Verasonics 'Trans' struct.
            %
            % xdc = Transducer.VERASONICS(Trans, c0) uses c0 in m/s as the 
            % sound speed when converting from wavelengths to meters. This
            % is typicaly set by the Verasonics property
            % 'Resource.Parameters.speedOfSound'. Be sure to explicitly set
            % this if other than 1540. 
            %
            % Example:
            %
            % Trans = struct('name', 'L12-3v', 'units', 'mm');
            % Trans = computeTrans(Trans);
            % xdc = Transducer.Verasonics(Trans);
            %
            % See also SCAN.VERASONICS TRANSDUCER.UFF TRANSDUCER.QUPS2USTB
            arguments
                Trans (1,1) struct
                c0 (1,1) double = 1540
            end

            % get distance property scaling
            switch Trans.units
                case 'wavelengths', scale = c0 / Trans.frequency * 1e-6; % lambda -> m
                case 'mm', scale = 1e-3; % mm -> m
            end

            % set the element length: search for an indication of the element height
            if isfield(Trans, 'elementLength') % best - leave as is
            elseif isfield(Trans, 'elevationApertureMm') % backup
                Trans.elementLength = Trans.elevationApertureMm * (1e-3 / scale); 
            else % fallback - make it square
                Trans.elementLength = Trans.elementWidth; 
            end

            % set required arguments
            switch Trans.type % 0=Lin(y=z=0),1=CrvdLn(y=0),2=2D,3=ann,4=R/C
                case 0,     xdc = TransducerArray.Verasonics(Trans  , c0);
                case 1,     xdc = TransducerConvex.Verasonics(Trans , c0);
                case {2,4}, xdc = TransducerMatrix.Verasonics(Trans , c0);
                otherwise,  xdc = TransducerGeneric.Verasonics(Trans, c0);
            end
            
            % set optional arguments
            % parse the impulse response
            if isfield(Trans, 'IR1wy')
                h = Trans.IR1wy; % impulse response
                t0 = - (argmax(hilbert(h))-1) / 250e6; % offset to peak time
                xdc.impulse = Waveform('t', t0 + (0:numel(h)-1) / 250e6, 'samples',h); % impulse response
            else
                xdc.impulse = Waveform.Delta();
            end

            % set the elevation focus
            if isfield(Trans, 'elevationFocusMm'), xdc.el_focus = 1e-3 * Trans.elevationFocusMm; end

            % set the bandwidth
            if isfield(Trans, 'Bandwidth'), xdc.bw = 1e6*Trans.Bandwidth([1 end]); end
        end
    end

    % SIMUS functions
    methods
        function getSIMUSParam(xdc)
        % GETSIMUSPARAM - Create a MUST compatible parameter struct
        %
        % p = GETSIMUSPARAM(xdc) returns a structure with properties
        % for a call to simus().
        %
        % Example:
        % 
        % xdc = TransducerArray();
        % p = xdc.getSIMUSParam();
        % 
        % See also ULTRASOUNDSYSTEM/SIMUS
            error( ...
                "QUPS:Transducer:unsupportedTransducer", ...
                "SIMUS does not support a " + class(xdc) + "." ...
                );
        end
    end

    % UFF constructor
    methods(Static)
        function xdc = UFF(probe)
            % UFF - Construct a Transducer from a UFF probe
            %
            % xdc = TRANSDUCER.UFF(probe) converts the uff.probe probe to a
            % Transducer xdc.
            %
            % See also TRANSDUCER.QUPS2USTB
            arguments, probe uff.probe; end
            switch class(probe) % dispatch
                case 'uff.linear_array',     xdc = TransducerArray.UFF(probe);
                case 'uff.curvilinear_array',xdc = TransducerConvex.UFF(probe);
                case 'uff.matrix_array',     xdc = TransducerMatrix.UFF(probe);
                case 'uff.probe',            xdc = TransducerGeneric.UFF(probe);
                otherwise,                   xdc = TransducerGeneric.UFF(probe);
            end
        end
    end

    % FDTD functions
    methods
        function [mask, el_weight, el_dist, el_ind] = elem2grid(xdc, grid, el_sub_div)
            % ELEM2GRID - Transducer element to grid mapping
            %
            % [mask, el_weight, el_dist, el_ind] = ELEM2GRID(xdc, grid)
            % returns a binary mask and a set of indices, weights, and
            % distances defined for each element give a Transducer xdc and
            % a ScanCartesian grid. The element indicies are defined on the
            % vectorized non-zero indices of the mask i.e. on the indices
            % of `find(mask)`.
            %
            % [...] = ELEM2GRID(xdc, grid, el_sub_div) subdivides the
            % elements in width and height by el_sub_div. The default is
            % [1,1].
            % 
            % See also ULTRASOUNDSYSTEM/KSPACEFIRSTORDER

            arguments
                xdc (1,1) Transducer
                grid (1,1) ScanCartesian
                el_sub_div (1,2) {mustBeInteger, mustBePositive} = [1,1]
            end

            % get the sensor and source masks. This is the hard part: how 
            % do I do this for a convex probe on a grid surface?

            vec = @(x)x(:); % define locally for compatibility

            % cast grid points to single type for efficiency and base grid on the origin
            [gxv, gyv, gzv] = deal(single(grid.x), single(grid.y), single(grid.z)); % grid {dim} vector
            [gx0, gy0, gz0] = deal(grid.x(1), grid.y(1), grid.z(1)); % grid {dim} first point
            [gxd, gyd, gzd] = deal(single(grid.dx), single(grid.dy), single(grid.dz)); % grid {dim} delta
            % pg = scan.positions(); % 3 x Z x X x Y
            pdims = arrayfun(@(c){find(c == grid.order)}, 'XYZ'); % dimensions mapping
            [xdim, ydim, zdim] = deal(pdims{:});
            iord = arrayfun(@(d)find([pdims{:}] == d), [1,2,3]); % inverse mapping

            % get array sizing - size '1' and step 'inf' for sliced dimensions
            [Nx, Ny, Nz] = deal(grid.nx, grid.ny, grid.nz); %#ok<ASGLU>

            % get local element size reference
            [width_, height_, sz_] = deal(xdc.width, xdc.height, grid.size); %#ok<ASGLU>

            % get regions for each receive transducer element
            el_cen = xdc.positions(); % center of each element
            [theta, phi, el_dir, el_wid, el_ht] = xdc.orientations(); %#ok<ASGLU> % orientations vectors for each element

            % reorder to map between kwave and conventional ultrasound coordinates
            % [el_dir, el_cen, el_wid, el_ht] = dealfun(@(v)v([3 1 2],:), el_dir, el_cen, el_wid, el_ht);

            % get edges in x/y/z (i.e. patches) for each sub element
            patches = xdc.patches(el_sub_div)'; % (E x N)
            nSub = size(patches, 1); % number of subelements

            % convert to cell array of points with x/y/z in 1st dimension
            p_patches = cellfun(@(pch) vec(cellfun(@(p)mean(p(:)), pch(1:3))), patches, 'UniformOutput', false);

            % send grid data to workers if we have a process pool
            if(~isempty(gcp('nocreate')) && isa(gcp('nocreate'), 'parallel.ProcessPool'))
                sendDataToWorkers = @parallel.pool.Constant;
            else
                sendDataToWorkers = @(x) struct('Value', x);
            end
            [gxv, gyv, gzv] = dealfun(sendDataToWorkers, gxv, gyv, gzv);

            % get sensing map for all patches for all elements
            parfor el = 1:xdc.numel
                % set variables for the element
                [el_dir_, el_cen_, el_wid_, el_ht_, p_patches_] = deal(el_dir(:,el), el_cen(:,el), el_wid(:,el), el_ht(:,el), p_patches(:,el)); %#ok<ASGLU>

                % initialize
                % mask{el} = false(sz_);

                for i = nSub:-1:1
                    % get the center of the subelement
                    pcen = p_patches_{i}; % (3 x 1)

                    % get zero crossing in indices
                    xind = 1 + (pcen(1) - gx0) / gxd;
                    yind = 1 + (pcen(2) - gy0) / gyd;
                    zind = 1 + (pcen(3) - gz0) / gzd;

                    % get integer index on both sides of zero crossing
                    [xind, yind, zind] = dealfun(@(n) vec((floor(n) + [0; 1])), xind, yind, zind);

                    % set index of sliced dimensions
                    if isinf(gxd), xind(:) = 1; end
                    if isinf(gyd), yind(:) = 1; end
                    if isinf(gzd), zind(:) = 1; end

                    % shift to appropriate dimensions for the scan
                    pind = {xind, yind, zind};
                    pind = pind(iord); % shift to proper order given by the scan
                    [pind{:}] = ndgrid(pind{:}); % outer product expansion
                    [pind{:}] = dealfun(vec, pind{:}); % vectorize
                    ind_msk = sub2ind(sz_, pind{:}); % get linear indices of the scan (J x 1)
                    pind_ = [pind{:}]; % combine
                    [xind, yind, zind] = deal(pind_(:,xdim), pind_(:,ydim), pind_(:,zdim)); % split

                    % get vector from element to the grid pixels
                    pgrd  = [gxv.Value(xind'); gyv.Value(yind'); gzv.Value(zind')]; %#ok<PFBNS> all data passed is necessary
                    vec_ = gather(pgrd) - pcen; % 3 x J

                    % get plane wave phase shift distance as the inner product
                    % sign is whether in front or behind
                    d = (el_dir_' * vec_).'; % J x 1

                    % get subelement apodization accounting for cosine
                    % distribution along the transducer surface
                    % I don't ... actually know how to do this ...

                    % a = cosd(90 * 2 * el_wid_' * (pcen - el_cen_) / width_ ) ...
                    %   * cosd(90 * 2 * el_ht_'  * (pcen - el_cen_) / height_);
                    % a = a + false(size(d));

                    % get subelement apodization accounting for grid
                    % interpolation
                    a = prod(1 - (abs(pgrd - pcen) ./ [gxd; gyd; gzd]),1)'; % J x 1
                    % a = ones(size(d)); %%% DEBUG %%%

                    % save outputs
                    mask_ind{i,el} = ind_msk; % indices of the scan
                    weights {i,el} = a; % amplitudes
                    dists   {i,el} = d; % distances for phase shifting response
                end
            end

            % check that each (sub-)element has some an associated grid index
            assert(~any(cellfun(@isempty, mask_ind), 'all'), 'Unable to assign all (sub-)elements to the grid.');

            % reduce to make full recording sensor mask
            mask = false(grid.size);
            mask(unique(cat(1, mask_ind{:}))) = true;

            % get the translation map for grid mask index to receiver element
            find_mask = find(mask);
            parfor m = 1:numel(mask_ind) %#ok<CPROPLC> 
                ind_el{m} = gather(argn(2, @ismember, mask_ind{m}, find_mask));
            end
            ind_el = reshape(ind_el, size(mask_ind));
            
            % make all into arrays, defined by mask indices
            parfor el = 1:xdc.numel
                el_weight{el} = cat(1, weights{:,el});
                el_dist  {el} = cat(1, dists  {:,el});
                el_ind   {el} = cat(1, ind_el {:,el});
            end

            % make numeric arrays (J x N)
            el_weight = cat(2, el_weight{:});
            el_dist   = cat(2, el_dist  {:});
            el_ind    = cat(2, el_ind   {:});
        end
    end

    % kWaveArray constructor
    methods
        function karray = kWaveArray(xdc, dim, og, karray_args)
            % KWAVEARRAY - Create a kWaveArray object from the Transducer
            %
            % karray = KWAVEARRAY(xdc, dim, og) creates a kWaveArray from the
            % Transducer xdc defined on a kWaveGrid with dim dimensions and
            % origin og.
            %
            % karray = KWAVEARRAY(..., Name, Value, ...) forwards Name/Value
            % pair arguments to the kWaveArrray constructor.
            %
            % Inputs:
            %   Axisymmetric
            %   BLITolerance
            %   BLIType
            %   UpsamplingRate
            %
            % See also KWAVEARRAY
            arguments
                xdc Transducer
                dim {mustBeInteger, mustBeInRange(dim, 2,3)}
                og (:,1) {mustBeNumeric}
                karray_args.Axisymmetric (1,1) logical
                karray_args.UpsamplingRate (1,1) double
                karray_args.BLITolerance (1,1) double {mustBeInRange(karray_args.BLITolerance, 0, 1)}
                karray_args.BLIType (1,1) string {mustBeMember(karray_args.BLIType, ["sinc", "exact"])}
            end

            % TODO: set axisymmetric if it makes sense
            % TODO: allow for elevation

            % initialize
            if isfield(karray_args,'BLIType')
                karray_args.BLIType = char(karray_args.BLIType);
            end
            karray_args = struct2nvpair(karray_args);
            karray = kWaveArray(karray_args{:});

            % get positions and orientations
            p = xdc.positions; % TODO: translate the entire array?
            [th, phi] = xdc.orientations;
            n = xdc.numel;
            if any(phi), warning('kWaveArray initialization with elevation angle is untested.'); end

            % build the rotation vectors
            % v = [cosd(phi); cosd(phi); sind(phi)] .* [cosd(th); sind(th); 1];
            elrot = [phi; th; th] .* [1;1;0];

            % convert to kWaveGrid coorindates: axial x lateral x elevation
            switch dim
                case 2 % we only rotate in x-z plane
                    [p, elrot] = deal(p([3,1  ],:), elrot(2,:));
                    [arr_off, arr_rot] = deal(-og(1:2), 0);
                case 3
                    [p, elrot] = deal(p([3,1,2],:), elrot([3,1,2],:)); % I think this is right?
                    [arr_off, arr_rot] = deal(-og(1:3), [0;0;0]);
            end

            % add each element to the array
            for i = 1:n
                karray.addRectElement(p(:,i), xdc.width, xdc.height, elrot(:,i));
            end

            karray.setArrayPosition(arr_off, arr_rot);
        end
    end

    % fullwave functions
    methods(Hidden)
        % GETFULLWAVETRANSDUCER - define a Fullwave transducer structure
        %
        % xdc_fw = GETFULLWAVETRANSDUCER(xdc, scan) creates a fullwave
        % compatible structure xdc_fw defined on the ScanCartesian scan.
        %
        % xdc_fw is a Fullwave compatible struct with (at least) the
        % following fields:
        %       - npx           number of elements
        %       - inmap         mask of the input pixels
        %       - nInPx         number of input pixels
        %       - nOutPx        number of output pixels
        %       - incoords      (x,y,1,el) coordinate pairs of the input pixels
        %       - outcoords     (x,y,1,el) coordinate pairs of the output pixels
        %
        % Note: fullwave support is incomplete, and is subject to change or
        % removal.
        % 
        % See also ULTRASOUNDSYSTEM/FULLWAVESIM
        function xdc_fw = getFullwaveTransducer(xdc, scan) %#ok<INUSD,STOUT>
            error( ...
                "QUPS:Transducer:notImplemented", ...
                "Fullwave support is not implemented for a " + class(xdc) + "." ...
                );
        end
    end

    % Field II calls - rely on being able to get a FieldII aperture
    methods
        function p = getFieldIIPositions(xdc)
            % GETFIELDIIPOSITIONS - get an array of element positions
            % 
            % p = getFieldIIPositions(xdc) returns a 3 x N vector of the 
            % positions of the N elements as represented in FieldII.
            %
            % See also GETFIELDIIAPERTURE FIELD_INFO

            ap = xdc.getFieldIIAperture();
            data = xdc_get(ap, 'rect');
            p = data([24 25 26], :) + xdc.offset;
            xdc_free(ap);
        end

        function el_patches = getFieldIIPatches(xdc, el_sub_div)
            % GETFIELDIIPATCHES - get cell array of patches
            %
            % el_patches = GETFIELDIIPATCHES(xdc) returns the
            % corners of the mathematical elements represented in FieldII.
            % Each element is represented by a {x,y,z,c} tuples (a.k.a. a
            % "patch") which holds [x1, x2; x3 x4], which are the corners
            % of the elements.
            %
            % el_patches = GETFIELDIIPATCHES(xdc, el_sub_div) uses the 
            % 2-element array el_sub_div to determine the number of 
            % mathematical elements. The result is a (N x E) array where N
            % is the number elements and E is the number of subelements.
            %
            % See also GETFIELDIIBARYCENTERS
            arguments
                xdc Transducer
                el_sub_div (1,2) double {mustBeInteger, mustBePositive} = [1,1]
            end

            % get the aperture
            ap = xdc.getFieldIIAperture(el_sub_div);

            % extract corners
            data = xdc_get(ap, 'rect');
            ind = 1:size(data,2);

            % release the aperture
            xdc_free(ap);

            % extract the corners into {x,y,z,c} tuples (a.k.a. a patch)
            % where x/y/z/c are each 2 x 2 arrays of points
            corner_map = @(i){...
                xdc.offset(1) + ...
                [   data(11,i), data(20,i); ...
                data(14,i), data(17,i)], ...
                xdc.offset(2) + ...
                [   data(12,i), data(21,i); ...
                data(15,i), data(18,i)], ...
                xdc.offset(3) + ...
                [   data(13,i), data(22,i);...
                data(16,i), data(19,i)], ...
                data(5, i) * ones(2,2),...
                };
            el_patches = arrayfun(corner_map, ind, 'UniformOutput', false); % array of {x,y,z,c} tuples

            % reshape by element (N x E), E is subelements
            el_patches = reshape(el_patches, prod(el_sub_div), [])';
        end
    end

    % internal subroutines
    methods
        function impulse = xdcImpulse(xdc)
            % XDCIMPULSE - create an impulse response Waveform
            %
            % impulse = XDCIMPULSE(xdc) creates a gaussian
            % pulse Waveform with the bandwidth and fractional bandwidth of
            % the Transducer xdc.
            %
            % See also WAVEFORM.DELTA()

            % defaults
            arguments
                xdc Transducer
            end

            % get energy strength cutoff time
            bwr = -6;
            tpr = -80;

            tc = gauspuls('cutoff',xdc.fc, xdc.bw_frac, bwr, tpr);
            xfc = xdc.fc;
            bwf = xdc.bw_frac;

            % get the function
            impulse_fun = Transducer.cgauspulsfun(xfc, bwf, bwr);

            % make a Waveform object
            impulse = Waveform('fun', impulse_fun, 't0', -tc, 'tend', tc);
        end
    
        function sub_div = getLambdaSubDiv(xdc, c0, p)
            % GETLAMBDASUBDIV - Get subelement divisions w.r.t. wavelength
            %
            % sub_div = GETLAMBDASUBDIV(xdc, c0, p) returns the element
            % subdivision sizes corresponding to a proportion p of the
            % wavelength given sound speed c0 (m/s).
            %
            % sub_div = GETLAMBDASUBDIV(xdc, c0, [pw ph]) uses a proportion
            % pw in the width dimension and ph in the height dimension.
            %
            % sub_div = GETLAMBDASUBDIV(xdc, c0) uses a default value of 
            % [pw ph] == p == 0.1 (10%). 
            %
            % Example:
            % % Get a system
            % us = UltrasoundSystem(); % a system
            % scat = Scatterers(); % a scaterrer
            % 
            % % Get subdivisions of rectangles of <= [lambda / 10, lambda]
            % sub_div = us.xdc.getLambdaSubDiv(scat.c0, [1/10 1]),
            % 
            % % Simulate
            % us.fs = single(us.fs); % reduce workload
            % chd = greens(us, scat, sub_div);
            % 
            % See also
            arguments
                xdc (1,1) Transducer
                c0 (1,1) {mustBePositive}
                p (1,2) {mustBePositive} = 0.1
            end

            % make odd
            plus2Odd = @(x) x + 1 - rem(x, 2); % odd -> odd, even -> odd

            % get divisions
            sub_div = plus2Odd(ceil([xdc.width, xdc.height] ./ (p .* c0 ./ xdc.fc)));
        end
    end

    % dependent methods
    methods
        function a = get.area(xdc), a = xdc.width * xdc.height; end

        function b = get.bw_frac(xdc), b = range(xdc.bw) ./ xdc.fc; end

        function set.bw_frac(xdc, bf), xdc.bw = mean(xdc.fc) + mean(xdc.fc) * bf * ([-1 1] ./ 2); end

        function o = get.origin(xdc), o = - xdc.offset; end

        function set.origin(xdc, o), xdc.offset = -o; end
    end

    % heterogeneous support
    methods (Static,Sealed,Access = protected)
        function xdc = getDefaultScalarElement()
            xdc = TransducerGeneric(); % default heterogeneous instance
        end
    end

    % plot functions
    methods
        function hp = patch(xdc, varargin, patch_args, xdc_args)
            % PATCH - Overload of patch for a transducer
            %
            % PATCH(xdc) plots the aperture from its patch representation.
            % It overloads the patch function so that arguments valid for
            % patch are valid here.
            %
            % PATCH(xdc, ax) plots on the axes ax.
            %
            % PATCH(..., 'sub_div', sub_div) specifices the number of
            % element subdivisions in width and height. The default is 
            % [1 1]. 
            %
            % PATCH(..., Name, Value, ...) passes Name/Value pair arguments
            % to patch.
            %
            % h = PATCH(...) returns a patch object handle
            %
            % Example:
            % xdc = TransducerMatrix.PO1921();
            % xdc.rot = [30, 15]
            % patch(xdc)
            % 
            % See also PATCH

            arguments
                xdc (1,1) Transducer
            end
            arguments(Repeating)
                varargin
            end
            arguments
                patch_args.?matlab.graphics.primitive.Patch
                patch_args.DisplayName = 'Elements'
                patch_args.FaceColor   = 'flat'
                patch_args.EdgeColor   = 'black'
                xdc_args.el_sub_div (1,2) double {mustBeInteger, mustBePositive} = [1 1];
            end

            % extract axes
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes') %#ok<CPROPLC> use built-in numel
                axs = varargin{1}; varargin(1) = [];
            else, axs = gca;
            end

            % get the aperture as a cell array of {x,y,z,c} pairs
            el_patches = xdc.patches(xdc_args.el_sub_div);

            % map point ordering to generate a square with patch; the
            % corners must be specified in sequential order
            % make column vector, drop apodization,
            ind_map = [1 2 4 3]';
            el_patches = cellfun(@(patch) cell2mat(cellfun(@(component) ...
                reshape(component(ind_map),[],1), ... .* 1e3, ...
                patch, 'UniformOutput', false)), ...
                el_patches, 'UniformOutput', false);

            % convert to get vertex x index x X/Y/Z array
            el_patches = permute(cat(3, el_patches{:}), [1,3,2]);

            % get X/Y/Z/C components
            Xp = el_patches(:,:,1);
            Yp = el_patches(:,:,2);
            Zp = el_patches(:,:,3);
            % Cp = el_patches(:,:,4); % apodization - not used
            % TODO: add option for showing apodization weighting

            % plot a patch: use depth as the color
            args = struct2nvpair(patch_args);
            hp = patch(axs, 'XData', Xp, 'YData', Zp, 'ZData', Yp, 'CData', Zp, varargin{:}, args{:});

            % set default axis arguments
            xlabel(axs, 'x');
            ylabel(axs, 'z');
            zlabel(axs, 'y');
            grid  (axs, 'on');
            % zlim  (axs, [-1e-3, 1e-3] + [min(Zp(:)), max(Zp(:))])
            axis  (axs, 'equal')
            view  (axs, 3);
            
        end

        function varargout = plot(xdc, varargin, plot_args)
            % PLOT - overload the plot function
            %
            % PLOT(xdc) plots the locations of the Transducer xdc.
            %
            % PLOT(xdc, ax) uses the axes ax instead of the current axes.
            %
            % PLOT(..., Name, Value, ...) passes name-value pair arguments
            % to the built-in plot function so that name value pairs that
            % are valid for plot are valid here.
            %
            % h = PLOT(...) returns the handle to the plot.
            %
            % Plots only the x-z slice.
            %
            % See also TRANSDUCER.PATCH 

            arguments
                xdc (1,1) Transducer
            end
            arguments(Repeating)
                varargin
            end
            arguments
                plot_args.?matlab.graphics.chart.primitive.Line
                plot_args.DisplayName = 'Elements'
                plot_args.LineStyle = 'none'
                plot_args.Marker = 'square'
            end
            
            % extract axes
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes') %#ok<CPROPLC> use built-in numel 
                axs = varargin{1}; varargin(1) = []; 
            else, axs = gca;
            end

            % get positions
            p = xdc.positions();

            % plot
            plot_args = struct2nvpair(plot_args);
            hp = plot(axs, p(1,:), p(3,:), varargin{:}, plot_args{:});

            % return
            if nargout > 0, varargout{1} = hp; end
        end
    end

    % to be deprecated functions
    methods(Static, Access=private)
        function f = cgauspulsfun(fc, bw, bwr)
            isig = ((4*pi*pi*(-bw*bw*fc*fc / (8*log(10.^(bwr/20))))) ./ 2);
            f = @(t) exp(-t .* t .* isig) .* (cospi(2*fc*t) + 1i .* sinpi(2*fc*t));
        end
    end
    methods(Hidden)
        function varargout = ultrasoundTransducerImpulse(varargin)
            warning("QUPS:Transducer:DeprecatedMethod", "Transducer.ultrasoundTransducerImpulse is deprecated - use xdcImpulse instead.");
            [varargout{1:nargout}] = xdcImpulse(varargin{:});
        end
    end
end