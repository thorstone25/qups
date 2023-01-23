function [c1, x1, y1, c2, x2, y2] = xiaolinwu_k_scaled(xa, ya, xb, yb, k, dx, dy)
% XIAOLINWU_K_SCALED - Kth line segment weights of Xiaolin Wu's algorithm
%
% [c1, x1, y1, c2, x2, y2] = XIAOLINWU_K_SCALED(xa, ya, xb, yb, k)
% computes the weights for the kth set of pixels described by the line
% segment from (xa,ya) to (xb,yb) and outputs them in the triples
% (x1,y1,c1) and (x2,y2,c2), where 0 <= ci <= 1 if ci is valid. For invalid
% points, ci is NaN.
% 
% All points are in pixel coordinates i.e. the x-grid is implicitly defined
% as 
% 
% x = floor(min(xa,xb)):ceil(max(xa,xb))
% 
% and is defined similarly for the y-grid. The calling function should 
% determine the total number of points necessary to traverse this grid.
%
% [...] = XIAOLINWU_K_SCALED(xa, ya, xb, yb, k, dxy) scales the size of
% the grid in x and y so that each pixel step is of size dxy in x and dxy 
% in y. The outputs ci are then scaled by the length of the line through
% the pixel. 
%
% [...] = XIAOLINWU_K_SCALED(xa, ya, xb, yb, k, dx, dy) scales the size of
% the grid in x and y so that each pixel step is of size dx in x and dy in
% y. The outputs ci are then scaled by the length of the line through the
% pixel.
%
% Example:
% % Create a grid
% [x, y] = deal(-5:5, -5:5); % 11 x 11 grid
% [xa, ya] = deal(-4, +1);
% [xb, yb] = deal(+3, -2);
% 
% % Get the interpolation weights
% K = max(numel(x), numel(y)) + 1; % total number
% [c1, ix1, iy1, c2, ix2, iy2] = arrayfun(...
%     @(k) xiaolinwu_k_scaled(xa, ya, xb, yb, k),...
%     (1:K)); % arrayfun computes this in parallel if on a gpuArray
% cxy = [c1, c2]; % combine output
% [~, ixo] = ismember([ix1 ix2], x); % find x-indices of the grid
% [~, iyo] = ismember([iy1 iy2], y); % find y-indices of the grid
% 
% % Create a sparse matrix, implicitly summing weights from neighboring
% % line segments
% val = (ixo ~= 0) & (iyo ~= 0) & ~isnan(cxy); % filter out invalid indices
% s = sparse(ixo(val), iyo(val), cxy(val), numel(x), numel(y));            
%
% figure;
% pcolor(x, y, full(s)')
% title("Weights from ("+xa+","+ya+"), to ("+xb+","+yb+").");
% 
% See also WBILERP SPARSE

% all argument validation and type checks disabled in order to run in
% parallel on GPU
%{
% arguments
%     xa (1,1) {mustBeReal, mustBeFinite, mustBeNumeric, mustBeFloat}
%     ya (1,1) {mustBeReal, mustBeFinite, mustBeNumeric, mustBeFloat}
%     xb (1,1) {mustBeReal, mustBeFinite, mustBeNumeric, mustBeFloat}
%     yb (1,1) {mustBeReal, mustBeFinite, mustBeNumeric, mustBeFloat}
%     k  (1,1) {mustBeInteger}
%     dx (1,1) {mustBeReal, mustBeFinite, mustBeNumeric, mustBeFloat} = ones(1,'like', xa);
%     dy (1,1) {mustBeReal, mustBeFinite, mustBeNumeric, mustBeFloat} = dx;
% end

% defaults
% if nargin < 6, dx = ones('like', xa); end
% if nargin < 7, dy = ones('like', dx); end

% Type checks and casting
% assert(all(cellfun(@isscalar, {xa,ya,xb,yb,dx,dy})), 'All coordinates and distances must be scalar.')
% assert(all(cellfun(@isfloat , {xa,ya,xb,yb,dx,dy})), 'All coordinates and distances must be a floating point type.');
%}

% convert to 0-based indexing
k = k - 1;

% integer zero type
izero = 0*k;

% whether to swap x/y axis
steep = abs(yb - ya) > abs(xb - xa);

% whether to swap a/b order
if steep, reverse = ya > yb; else, reverse = xa > xb; end

% sort coordinates and endpoints
if      steep &&  reverse
    [ux, uy, vx, vy, dx, dy] = deal6(yb, xb, ya, xa, dy, dx);
elseif  steep && ~reverse
    [ux, uy, vx, vy, dx, dy] = deal6(ya, xa, yb, xb, dy, dx);
elseif ~steep &&  reverse
    [ux, uy, vx, vy, dx, dy] = deal6(xb, yb, xa, ya, dx, dy);
else % ~steep && ~reverse
    [ux, uy, vx, vy, dx, dy] = deal6(xa, ya, xb, yb, dx, dy);
end
    
% now, we have x0 <= x1, dy < dx, 0 < dy/dx <= 1
% get the gradient
if vx == ux 
    g = 1; 
else
    g = double((vy - uy) / (vx - ux));
end
l = hypot(dx, g*dy); % length through the pixel
sx = double(floor(ux+1/2)) + izero; % first (index) endpoint
ex = double(floor(vx+1/2)) + izero; % last (index) endpoint

% get fraction values for the pair of pixels
ix = sx + k; % x value
ixf = double(ix) + 0*xa; % x value
yf = uy + g * (ixf - ux); % y-fractional
iy = double(floor(yf)) + izero; % y value
cf = yf - floor(yf); % c-fractional
% x-fractional
if(k == 0 && k == (ex - sx)) % special case: less than one pixel
    xgap = vx - ux;
elseif(k == 0) % first pixel in line
    xgap = (1 - ((ux+1/2) - ixf));
elseif(0 < k && k < (ex - sx)) % mid line
    xgap = 1 + 0*xa;
elseif(k == (ex - sx)) % last pixel in line
    xgap = (vx+1/2) - ixf;
else % beyond the line
    xgap = nan + 0*xa;
end

% write outputs
[x1, y1, c1, ~] = deal4(ix, iy, (1-cf) * l * xgap, 0);
[x2, y2, c2, ~] = deal4(ix, iy+1,  cf  * l * xgap, 0);

% set correct output for x and y in original coordinates
if steep, [x1, y1, x2, y2] = deal4(y1, x1, y2, x2); end

function [a,b,c,d] = deal4(a,b,c,d), end
function [a,b,c,d,e,f] = deal6(a,b,c,d,e,f), end

end

