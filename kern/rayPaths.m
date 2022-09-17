function [w, d] = rayPaths(xg, zg, pj, pa, kwargs)
% RAYPATHS - ray path integral weights
%
% w = RAYPATHS(xg, zg, pj, pa) computes a set of weights w defined on the 
% grid {xg, zg} for all pairs of line segments from the each point {pj} to
% each point {pa}. w is a sparse matrix of size [Z x X x J] x A where Z =
% numel(zg), X = numel(xg), J = numel(pj)/2 and A = numel(pa)/2.
%
% [w, d] = RAYPATHS(...) also returns a full matrix d of the length of each
% ray as a J x A array.
% 
% [...] = RAYPATHS(..., 'method', method) chooses the interpolation weight
% method. Must be either "xiaolin" or "bilerp". The "xiaolin" method uses
% Xiaolin Wu's line drawing algorithm. The "bilerp" method gives bilinear
% interpolation weights. The default is "bilerp".
%
% [...] = RAYPATHS(..., 'gpu', tf) chooses whether to use a default. The
% default is true if one is available.
%
% [...] = RAYPATHS(..., 'ord', "XZ") swaps the order of the output array so
% that the x-grid is placed before the z-grid i.e. so that w is of size [X
% x Z x J] x A.
%
% [...] = RAYPATHS(..., 'bsize', B) computes in blocks of at most B points 
% at a time. A larger value of B offers better compute performance, but
% uses more memory. A lower value prevents out-of-memory (OOM) errors.
%
% [...] = RAYPATHS(..., 'verbose', true) prints updates about the progress
% of the computation.
%
% Example:
%
% % Create a grid
% x = -5:5;
% z = -5:5;
% [xa, ya] = deal(-4, +1);
% [xb, yb] = deal(+3, -2);
% 
% % Create a series of rays
% w = rayPaths(x, z, [[xa; ya], [xb; yb]], [0;0], 'ord', "ZX");
% w = reshape(full(w), numel(z), numel(x), []);
% 
% % Display it
% figure; 
% subplot(1,2,1); pcolor(x, z, w(:,:,1));
% title("Weights from (0,0) to ("+xa+","+ya+")."); 
% subplot(1,2,2); pcolor(x, z, w(:,:,2));
% title("Weights from (0,0) to ("+xb+","+yb+").");
% 
% See also WBILERP WBILERPG XIAOLINWU_K_SCALED

arguments
    xg (1,:) {mustBeFinite, mustBeReal} % 1 x X
    zg (1,:) {mustBeFinite, mustBeReal} % 1 x Z
    pj (2,:,:) {mustBeReal} % 2 x J x M x ...
    pa (2,:,:) {mustBeReal} % 2 x J x M x ...
    kwargs.bsize (1,1) {mustBePositive, mustBeInteger} = max(1,floor(2^32 /... 
         (max(size(pj,3),size(pa,3)) * ((1 + numel(xg) + numel(zg)) * 4) * 4 * 8) ... output size per batch in btyes
         ))
    kwargs.ord (1,1) string {mustBeMember(kwargs.ord, ["XZ", "ZX"])} = "ZX"
    kwargs.gpu (1,1) logical = logical(gpuDeviceCount())
    kwargs.method (1,1) string {mustBeMember(kwargs.method, ["bilerp", "xiaolin"])} = "bilerp"
    kwargs.verbose (1,1) logical = false;
    kwargs.prototype = zeros([0, 0], 'like', [pj(end), pa(end)]) % prototype
end

% output: I x J x M, I = Z x X
dproto = kwargs.prototype;

[X, Z] = dealfun(@numel, xg, zg);
I = X * Z;
szjm = max(size(pj,1:3), size(pa,1:3));
[J, M] = deal(szjm(2), szjm(3)); 

% buffer size: 4*(Nx + Nz + 1)
switch kwargs.method
    case "bilerp",  Kmax = 4*(X + Z + 1);
    case "xiaolin", Kmax = 2*(X + Z + 1); 
end

% output data sizing
switch kwargs.ord
    case "ZX", [xstride, zstride] = deal(Z, 1);
    case "XZ", [xstride, zstride] = deal(1, X);
end

% get input data sizing
osz = [Kmax, J, M];
[ix, iz] = deal(zeros(osz, 'int32'));
iw = zeros(osz, 'like', dproto);
d = zeros([osz(2:end)], 'like', dproto);

% split up data into blocks
js = num2cell(uint64((1:kwargs.bsize)' + (0:kwargs.bsize:J-1)), 1);
js{end} = js{end}(js{end} <= J);

% allocate output data
for j = numel(js):-1:1
    if kwargs.verbose, fprintf('\n'); end

    % get input data
    if size(pj,2) == 1,
        pj_ = pj;
    else
        pj_ = sub(pj,js{j},2);
    end
    if size(pa,2) == 1,
        pa_ = pa;
    else
        pa_ = sub(pa,js{j},2);
    end

    Jp = numel(js{j}); % == J'

    % convert to pixel coordinates and brodcast to all angles
    % (1 x J' x M)
    D = max(ndims(pa_), ndims(pj_));
    psz = [1, max(size(pa_,2:D),size(pj_,2:D))]; % broadcast size
    ZERO = zeros(psz, 'like', dproto); % broadcasting identity
    pax = sub(pa_,1,1) + ZERO;
    paz = sub(pa_,2,1) + ZERO;
    pbx = sub(pj_,1,1) + ZERO;
    pbz = sub(pj_,2,1) + ZERO;

    % get output sizing and data
    josz = [osz(1), Jp, osz(3)];

    % apply function for all start/end points
    tt = tic;
    switch kwargs.method
        case "bilerp"
            % preallocate output
            [ixo, izo] = deal(zeros(josz, 'like', ix));
            iwo = zeros(josz, 'like', iw);

            if kwargs.gpu && exist('wbilerp.cu', 'file') % call the GPU version
                [iwo(:), ixo(:), izo(:)] = wbilerpg(xg, zg, pax, paz, pbx, pbz);
            else % call the CPU version in parallel with defaults
                parfor (jm = (1:numel(pax)))
                    [iwo(:,jm), ixo(:,jm), izo(:,jm)] = wbilerp(xg, zg, pax(jm), paz(jm), pbx(jm), pbz(jm));
                end
            end
        
        case "xiaolin"
            % get the coordinate transform to 0-based coordinates
            [dx, dz] = deal(mean(diff(xg)), mean(diff(zg)));
            [x0, z0] = deal(xg(1), deal(zg(1)));

            % preallocate output
            josz(1) = Kmax/2;
            [ix1, iz1, ix2, iz2] = deal(zeros(josz, 'like', ix));
            [ c1,  c2] = deal(zeros(josz, 'like', iw));

            % use no parpool on GPU, default parpool on CPU, 
            if kwargs.gpu, pul = {0}; else, pul = {}; end
            
            % Get the interpolation weights for each part of the line
            parfor (k = 1 : Kmax, pul{:})
            % for k = 1:Kmax/2
                [c1k, ix1k, iz1k, c2k, ix2k, iz2k] ...
                    = arrayfun(@xiaolinwu_k_scaled, ...
                    round((pax-x0)/dx), round((paz-z0)/dz), ...
                    round((pbx-x0)/dx), round((pbz-z0)/dz), ...
                    k+ZERO, dx+ZERO, dz+ZERO);
                [ix1(k,:), iz1(k,:), c1(k,:), ix2(k,:), iz2(k,:), c2(k,:)] = ...
                deal(ix1k(:), iz1k(:), c1k(:), ix2k(:), iz2k(:), c2k(:));
            end

            iwo = [c1; c2]; % combine output
            [~, ixo] = ismember([ix1; ix2], round((xg-x0)/dx)); % find x-indices of the grid
            [~, izo] = ismember([iz1; iz2], round((zg-z0)/dz)); % find y-indices of the grid

    end
    
    % check that all values positive or nan
    assert(all(iwo >= 0, 'all'), ...
        "Found negative values - this is likely due to numerical precision issues. " ...
        + "Consider changing the grid axes and ensure that the grid contains all endpoints. " ...
        );
    assert(all(~isnan(iwo), 'all'), 'Detected nan values - this is a bug :''(');
    

    % move results to MATLAB indexing % (K x J' x M)
    ix = ixo - 1;
    iz = izo - 1;
    iw = iwo;

    % get ray lengths
    d(js{j},:) = gather(vecnorm(pa - pj_, 2, 1)); % (1 x J' x M)
    
    % filter indices outside range / nan values
    % (Is x J' x M)
    val = 0 <= ix & ix < X ...
        & 0 <= iz & iz < Z ...
        & ~isnan(iw);

    % delete invalid indices
    ix(~val) = [];
    iz(~val) = [];
    iw(~val) = [];

    % get linear indexing
    val = int32(find(val));

    % compute sparse tensor as matrix indices (& add
    % identical indices implicitly)
    [~,jv,mv] = ind2sub([Kmax,Jp,M], val);
    ijv = sub2ind([I,Jp],(1 + zstride .* iz(:) + xstride .* ix(:)), jv);
    % [indij, indm] = ind2sub([int32(Kmax)*spsz(2), spsz(3)], find(val));
    w{j} = sparse(ijv, mv, double(iw(:)), double(I*Jp), double(M)); % ([I x J'] x M)

    if kwargs.verbose, 
        fprintf('Completed batch %i of %i in %0.5f seconds.\n', ...
        j, numel(js), toc(tt)); 
    end
end

% concatenate all sparse arrays to make a final sparse array
% ([I x J'] x M) -> % ([I x J] x M)
w = cat(1, w{:}); 

if kwargs.verbose, fprintf('done!\n'); end

end
