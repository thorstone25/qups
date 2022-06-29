function z = slsc(x, dim, L, method, varargin)
% SLSC - Short-lag Spatial Coherence 
%
% INPUTS:
%       x       - input data with the dimensions [tx, rx, ??]
%       dim     - Aperture dimension (default 2)
%       L       - max number of lags to use (default 10)
%       method  - method for normalization of input data (default ensamble)
%   
% OPTIONAL INPUTS:
%   
%
% OUTPUTS:
%       z       - SLSC image
%
% EXAMPLE USE:
%
%       z = SLSC(x)                 computes the short-lag spatial coherence (SLSC) 
%                                   across the data x. x can be any ND-array.
%
%       z = SLSC(x, dim)            operate in dimension dim. The default is 2.
% 
%       z = SLSC(x, dim, L)         uses up to L lags
% 
%       z = SLSC(x, dim, L, method) uses the specified method. Must be one of
%                                   {"ensemble"* | "average"}.
% 
% References:
% D. Hyun, A. C. Crowley and J. J. Dahl, 
% "Efficient Strategies for Estimating the Spatial Coherence of Backscatter,"
% in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
% vol. 64, no. 3, pp. 500-513, March 2017.  
% doi: 10.1109/TUFFC.2016.2634004
%
% M. A. Lediju, G. E. Trahey, B. C. Byram and J. J. Dahl, 
% "Short-lag spatial coherence of backscattered echoes: imaging characteristics," 
% in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
% vol. 58, no. 7, pp. 1377-1388, July 2011.  
% doi: 10.1109/TUFFC.2011.1957
%

% defaults
if nargin < 2 || isempty(dim   ), dim    = 2;          end
if nargin < 3 || isempty(L     ), L      = 10;         end
if nargin < 4 || isempty(method), method = "ensemble"; end

% default parpool: this is data heavy, so only use a threadpool
if isa(gcp('nocreate'), 'parallel.threadPool'), clu = gcp(); else, clu = 0; end
kwargs.parcluster = clu;

% assign optional inputs
for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

% parse inputs
clu = kwargs.parcluster;

% get weighting filter across receiver pairs
K = gather(size(x,dim)); % the full size of the aperture
[M, N] = ndgrid(1:K, 1:K);
H = abs(M - N); % lags for each receiver/cross-receiver pair
S = (0 < H & H <= L); % valid lags for adding (i.e. the short lags)

% choose average or ensemble
switch method
    case "average"
        % normalize magnitude per sample (averaging)
        x = x ./ abs(x);
        x(isnan(x)) = 0; % 0/0 -> 0
        % TODO: test if abs(x) close to zero instead of assuming all nans
        % are from computing 0/0 

        % get weighting filter across receiver pairs
        W = S ./ (K - H); % final weights per pair (debiased)

        % place cross-receiver across cells
        % xc = num2cell(x, setdiff(1:ndims(x), dim));

        % place weights as cross receiver over cells, receiver in dim
        W = num2cell(shiftdim(W, -(dim-2)), dim);

        % correlation sum across receivers per cross-receiver kernel
        %vn = @(w,xc) sum(w .* conj(xc) .* x, dim);

        % compute product and sum per cross-receiver
        z = 0;
        parfor (i = 1:numel(W), clu)
            % z = z + vn(W{i}, xc{i});
            z = z + sum(W{i} .* conj(sub(x,i,dim)) .* x, dim);
        end

    case "ensemble"
        % place cross-receiver across cells
        % xc = num2cell(x, setdiff(1:ndims(x), dim));
        for i = gather(size(x,dim)):-1:1, xc{i,1} = sub(x,i,dim); end
        xc = shiftdim(xc, 1-dim);

        % place weights as cross receiver over cells, receiver in "dim"
        W = num2cell(shiftdim(S, -(dim-2)), dim);

        % correlation across receivers per cross-receiver kernel
        vn = @(w,xc) deal(...
            sum(w .* conj(xc) .* x , dim),...
            sum(w .* conj(x ) .* x , dim),...
            sum(w .* conj(xc) .* xc, dim) ...
            );

        % compute and sum for each cross-receiver
        [z, a, b] = deal(0);
        parfor (i = 1:numel(W), clu)
            [zr, ar, br] = vn(W{i}, xc{i});
            z = z + zr;
            a = a + ar;
            b = b + br;
        end

        % get final image
        z = z ./ sqrt(a .* b);
    otherwise, error("Unrecognized method '" + method + "'");
end
end

