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
    pj (2,:,1,:,:) {mustBeReal} % 2 x J x 1 x ...
    pa (2,:,1,:,:) {mustBeReal} % 2 x A x 1 x ...
    kwargs.bsize (1,1) {mustBePositive, mustBeInteger} = max(1,floor(2^33 /... 
         ((numel(pa) / 2) * ((1 + numel(xg) + numel(zg)) * 4) * 4 * 8) ... output size per batch in btyes
         ))
    kwargs.ord (1,1) string {mustBeMember(kwargs.ord, ["XZ", "ZX"])} = "ZX"
    kwargs.gpu (1,1) logical = logical(gpuDeviceCount())
    kwargs.method (1,1) string {mustBeMember(kwargs.method, ["bilerp", "xiaolin"])} = "bilerp"
    kwargs.verbose (1,1) logical = false;
    kwargs.prototype = zeros([0, 0], 'like', [pj(end), pa(end)]) % prototype
end

% output: I x J x A, I = Z x X
dproto = kwargs.prototype;

% buffer size: 4*(Nx + Nz + 1)
[X, Z, J, M] = dealfun(@numel, xg, zg, pj, pa);
[I, J, M] = deal(X*Z, J/2, M/2);

% move A to dim 3 - TODO: make this dim ndims(pa)+1
pa = permute(pa, [1 3 2 4:ndims(pa)]);

% output size from grid size
switch kwargs.method
    case "bilerp",  Kmax = 4*(X + Z + 1);
    case "xiaolin", Kmax = 2*(X + Z + 1); 
end

% get data sizing
osz = [Kmax, J, M];
[ix, iz] = deal(zeros(osz, 'int32'));
iw = zeros(osz, 'like', dproto);
d = zeros([osz(2:end)], 'like', dproto);

% split up data into blocks
js = num2cell(uint64((1:kwargs.bsize)' + (0:kwargs.bsize:J)), 1);
js{end} = js{end}(js{end} <= J);

% allocate output data
for j = 1:numel(js)
    if kwargs.verbose, fprintf('\n'); end

    % get input data
    pj_ = pj(:,js{j});

    % convert to pixel coordinates and brodcast to all angles
    % (1 x J' x M)
    D = max(ndims(pa), ndims(pj_));
    psz = [1, max(size(pa,2:D),size(pj_,2:D))]; % broadcast size
    ZERO = zeros(psz, 'like', dproto); % broadcasting identity
    pax = sub(pa,1,1)  + ZERO;
    paz = sub(pa,2,1)  + ZERO;
    pbx = sub(pj_,1,1) + ZERO;
    pbz = sub(pj_,2,1) + ZERO;

    % get output sizing and data
    josz = osz; 
    josz(2) = numel(js{j});
    % apply function for all start/end points
    tt = tic;
    switch kwargs.method
        case "bilerp"
            % preallocate output
            [ixo, izo] = deal(zeros(josz, 'like', ix));
            iwo = zeros(josz, 'like', iw);

            if kwargs.gpu && exist('wbilerp.cu', 'file') % call the GPU version
                [iwo(:), ixo(:), izo(:)] = wbilerpg(xg, zg, pax, paz, pbx, pbz);
            else % call the CPU version in parallel
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
            % parfor (k = 1 : Kmax, pul{:})
            for k = 1:Kmax/2
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
    if kwargs.verbose, 
        fprintf('Completed batch %i of %i in %0.5f seconds.\n', ...
        j, numel(js), toc(tt)); 
    end

    % check that all values positive or nan
    assert(all(iwo >= 0 | isnan(iwo), 'all'), 'Negative or nan values - this is a bug :''(');

    % move results to the CPU % (K x J' x M)
    ix(:,js{j},:) = gather(ixo-1);
    iz(:,js{j},:) = gather(izo-1);
    iw(:,js{j},:) = gather(iwo);

    % get ray lengths
    d(js{j},:) = gather(vecnorm(pa - pj_, 2, 1)); % (1 x J' x M)
end

% collect results in cells as sparse matrices in I
if kwargs.verbose, fprintf('\nStoring results ... '); end

% data sizing
switch kwargs.ord
    case "ZX", [xstride, zstride] = deal(Z, 1);
    case "XZ", [xstride, zstride] = deal(1, X);
end

% filter indices outside range / nan values
% (Is x J x M)
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
[~,jv,mv] = ind2sub([Kmax,J,M], val);
ijv = sub2ind([I,J],(1 + zstride .* iz(:) + xstride .* ix(:)), jv);
% [indij, indm] = ind2sub([int32(Kmax)*spsz(2), spsz(3)], find(val));
w = sparse(ijv, mv, double(iw(:)), double(I*J), double(M)); % ([I x J] x M)

if kwargs.verbose, fprintf('done!\n'); end

end