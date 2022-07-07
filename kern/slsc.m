function z = slsc(x, dim, L, method, kdim)
% SLSC - Short-lag Spatial Coherence 
%
% z = SLSC(x) computes the short-lag spatial coherence (SLSC) across the
% data x. x can be any ND-array.
%
% z = SLSC(x, dim) operates in dimension dim. The default is the first 
% non-singleton dimension.
% 
% z = SLSC(x, dim, L) uses up to L lags. If L is an array of indices, the
% indices matching L are used. Select a single lag by making an array
% including 0 e.g. L = [0 5]. The default is 
% max(1, floor(1/4 * size(x,dim)));
% 
% z = SLSC(x, dim, L, method) uses the specified method. Must be one of
% {"ensemble"* | "average"}.
% 
% z = SLSC(x, dim, L, method, kdim) designates a temporal sample dimension
% used in calculating the average SLSC. The default is ndims(x) + 1.
% 
% Example:
% 
% % generate a sine wave with random phase error
% T = 2^4; % number of samples
% fc = 2; % normalized frequency
% t = (0:T-1)'/T; % time vector
% phi = pi/2*rand([1,64]); % random phase differences
% x = exp(1j*(2*pi*fc*t + phi)); % sine wave
%
% % get the SLSC across the phases
% z = slsc(x, 2, [0 3], 'average')
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

% defaults
arguments
    x
    dim (1,1) {mustBePositive, mustBeInteger} = prod(find(size(x,1:ndims(x)) ~= 1, 1, 'first')); % lag dimension
    L (1,:) {mustBeNonnegative, mustBeInteger} = max(1, floor(size(x,dim)/4)); % lags or maximum lag
    method (1,1) string {mustBeMember(method, ["average", "ensemble"])} = "average"; % average or ensemble estimator
    kdim (1,1) {mustBePositive, mustBeInteger} = ndims(x) + 1; % time sample dimension
end

% default parpool: this is data heavy, so only use a threadpool
% don't use any pool for gpu data
if isa(gcp('nocreate'), 'parallel.threadPool') ...
        && ~isa(x, 'gpuArray') && ~isa(x, 'tall'), 
    clu = gcp(); else, clu = 0; 
end

% get weighting filter across receiver pairs
A = gather(size(x,dim)); % the full size of the aperture
[M, N] = ndgrid(1:A, 1:A);
H = abs(M - N); % lags for each receiver/cross-receiver pair
if isscalar(L), lags = 1:L; else, lags = L; end % chosen lags for adding (i.e. the short lags)
S = ismember(H,lags); % selection mask for the lags
L = length(lags); % number of lags - this normalizes the sum of the average estimates
K = gather(size(x,kdim)); % number of samples in time

% choose average or ensemble
switch method
    case "average"
        % normalize magnitude per sample (averaging)
        x = x ./ vecnorm(x,2,kdim) / K;
        x(isnan(x)) = 0; % 0/0 -> 0
        % TODO: test if norm(x) close to zero instead of assuming all nans
        % are from computing 0/0 

        % get weighting / filter across receiver pairs
        W = S ./ (A - H) / 2 / L; % weights per pair (debiased, pos & neg averaged, multi-correlation-averaged)

        % place weights as cross receiver over cells, receiver in dim
        W = num2cell(shiftdim(W, -(dim-2)), dim);

        % compute product and sum per cross-receiver
        z = 0;
        parfor (i = 1:A, clu)
            z = z + sum(W{i} .* conj(sub(x,i,dim)) .* x, [dim, kdim], 'omitnan');
        end

    case "ensemble"
        % scale the problem to maintain numerical accuracy
        % we do this as precisely as possible by shifting by a power of 2.
        x = x .* 2.^nextpow2(1/mean(vecnorm(x,2,dim),"all",'omitnan'));

        % place weights as cross receiver over cells, receiver in "dim"
        W = num2cell(shiftdim(S, -(dim-2)), dim);

        % compute and sum for each cross-receiver
        [z, a, b] = deal(0);
        %parfor (i = 1:numel(W), clu)
        for i = 1:A
            w = W{i}; % weight for this set of indices
            xc = sub(x,i,dim); % cross-receiver
            z = z + sum(w .* conj(xc) .* x , [dim, kdim],'omitnan'); % inner-product
            a = a + sum(w .* conj(x ) .* x , [dim, kdim],'omitnan'); % a-norm
            b = b + sum(w .* conj(xc) .* xc, [dim, kdim],'omitnan'); % b-norm
        end

        % normalize by the norm of the vector
        z = z .* rsqrt(a) .* rsqrt(b);
end
end

