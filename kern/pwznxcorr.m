function y = pwznxcorr(x, lags, W, U, kwargs)
% PWZNXCORR - Pair-wise Windowed Zero-Normalized Cross-Correlation 
% 
% y = PWZNXCORR(x, lags) correlates pairs of channels in the ND-array x
% within a moving time-domain window and with the given temporal sample
% offsets lags. The windowed vectors are optionally debiased and normalized
% prior to correlation.
% 
% This calculation is roughly equivalent to the following pseudo-code: 
% ```
% [T, N, F] = size(x); % time x channels x frames
% [L, W, S] = deal(numel(lags), numel(w), 1); % lags, window, stride
% for l = 1:L, for f = 1:F, for n = 1:N-S, for t = 1:T, % lags, frames, channels, time
%     u = x(t + (1:W)          , n  , f); % reference  signal (W x 1)
%     v = x(t + (1:W) + lags(l), n+S, f); % comparison signal (W x 1)
%     ...                                 % debiasing and normalization (optional)
%     y(t,n,f,l) = v'*diag(w)*u;          % weighted inner-product (1 x 1)
% end,         end,         end,           end,
% ```
% 
% y = PWZNXCORR(x, L) where L is a scalar specifies lags -L:L. 
% 
% y = PWZNXCORR(x, L, W) or y = PWZNXCORR(x, lags, W) where W is a scalar
% uses a window of size W in time to localize the correlation. W specifies
% the full length of the window. The default is ceil(max(abs(lags)) / 2).
%
% y = PWZNXCORR(x, L, w) or y = PWZNXCORR(x, lags, w) where w is a vector
% or NDarray uses the filter weights w to localize the correlation. If w is
% non-scalar across the channel dimension, adjacent channels are used up to
% the size of w in that dimension. The weights may be non-scalar only in
% the time and channel dimensions.
% 
% y = PWZNXCORR(x, L, w, U) or y = PWZNXCORR(x, lags, w, U) upsamples the
% comparison signal `v` by a factor U using 'pchip' interpolation prior to
% computing the discrete lag. The lags are then scaled by 1/U. The default
% is 1.
% 
% y = PWZNXCORR(..., 'zero', true) debiases the windowed vectors u and v
% such that sum(u) == sum(v) == 0. The default is true.
% 
% y = PWZNXCORR(..., 'norm', true) normalizes the windowed vectors u and v
% such that |u||_2 == ||v||_2 == 1. The default is true.
% 
% y = PWZNXCORR(..., 'tdim', tdim) uses dimension tdim as the time
% dimension. The default is 1.
%
% y = PWZNXCORR(..., 'ndim', ndim) uses dimension ndim as the channel
% dimension. The default is 2.
%
% y = PWZNXCORR(..., 'ldim', ldim) uses dimension ldim as the (output) lag
% dimension. The default is ndims(x) + 1.
%
% y = PWZNXCORR(..., 'lvec', true) vectorizes the lag dimension rather than
% computing each lag iteratively. The default is false.
%
% y = PWZNXCORR(..., 'ref', 'neighbor') compares each channel n to the
% neighboring channel n + 1. This is the default.
%
% y = PWZNXCORR(..., 'ref', 'neighbor', 'stride', S) compares each channel
% n to the channel n + S. The default is 1.
%
% y = PWZNXCORR(..., 'ref', 'center') compares each channel n to the median
% channel n0 = N/2. If the number of channels is odd, the mean of the
% center two channels [floor(N/2), ceil(N/2)] is used. This is roughly
% equivalent to modifying the psuedo-code above to contain the following:
% ```
% for l = 1:L, for f = 1:F, for n = 1:N, for t = 1:T, % lags, frames, channels, time
%     n0 = unique([floor(N/2), ceil(N/2)]);      % reference channel(s) (W x [1|2])
%     u =      x(t + (1:W)          , n , f)   ; % reference  signal (W x 1)
%     v = mean(x(t + (1:W) + lags(l), n0, f),2); % comparison signal (W x 1)
%     y(t,n,l,f) = v'*diag(w)*u;                 % weighted inner-product (1 x 1)
% end,         end,         end,         end,
% ```
%
% y = PWZNXCORR(..., 'ref', 'center', 'multi', true) compares each channel n 
% to the central M == size(w, ndim) channels rather than the median channel. 
% The default is false.
% 
% y = PWZNXCORR(..., 'ref', 'x0', 'x0', signal) compares each channel to
% the given signal x0. The dimensions of x0 must be compatible with x. This
% is roughly equivalent to modifying the psuedo-code above to contain the
% following:
% ```    
% for l = 1:L, for f = 1:F, for n = 1:N-1, for t = 1:T, % lags, frames, channels, time
%     u = x( t + (1:W)          , n, f); % reference  signal (W x 1)
%     v = x0(t + (1:W) + lags(l), :, f); % comparison signal (W x 1)
%     y(t,n,l,f) = v'*diag(w)*u;         % weighted inner-product (1 x 1)
% end,         end,         end,           end,
% ```
%
% y = PWZNXCORR(..., 'pad', false) disables 0-padding prior to computing
% the cross-correlation. If `imfilter` (Image Processing Toolbox) is not
% available, this introduces wrap-around artefacts, but prevents the need
% to copy x, reducing the memory footprint slightly. The default is true.
%
% y = PWZNXCORR(..., 'iflt', true) utilizes `imfilter` to apply convolution
% rather than native convolution.  The Image Processing Toolbox is required
% to use this option. The default is true if `imfilter` exists and false
% otherwise. 
% 
% Example:
% % Create correlated channels with gaussian phase noise
% [T, N, P] = deal(512, 64, 10); % time, channels, samples-per-period
% [fc, sig] = deal(2, 2); % frequency, phase noise
% fs = P*fc;
% tn = randi([-1 1]*floor(0.25*P), [1 N]); % random shift - differnce < 1/2 period
% t  = ((0 : T-1)' + tn) ./ fs;
% theta = fc * t;
% eta = 1e-2*sig*randn(size(theta)); % std of 5%
% x = sinpi(2*(theta + eta)); 
% 
% % Get correlation
% L = floor(0.5*P); % max temporal lag - don't exceed more than 1/2 period
% W = floor(5*P); % temporal window size (in samples)
% y = pwznxcorr(x, L, W);
% 
% % The peak lag across channels should be identical to the random shift
% [v, i] = max(y, [], 3);
% lags = -L : L;
% assert(isequal(lags(mode(i,1)), -diff(tn)));
% 
% % Show the peak lag
% figure;
% imagesc(v);
% title('Peak Correlation');
% xlabel('Channels (#)');
% ylabel('Time Sample (#)');
% 
% See also CONVN XCORR

arguments
    x {mustBeNumeric, mustBeFloat} % data
    lags (1,:) double {mustBeFinite} % lags
    W {mustBeNumeric} = max([ceil(max(abs(lags)) / 2), 1]) % window size or weighting vector
    U (1,1) double {mustBeInteger, mustBePositive} = 1 % upsampling ratio
    kwargs.pad (1,1) logical  = true; % whether to zero-pad the data prior to correlation
    kwargs.zero (1,1) logical = true; % whether to 0-mean the data
    kwargs.norm (1,1) logical = true; % whether to normalize the data
    kwargs.ref (1,1) string {mustBeMember(kwargs.ref, ["neighbor", "center", "x0"])} = "neighbor" % reference channel
    kwargs.stride (1,1) double {mustBeInteger, mustBePositive} = 1
    kwargs.x0 {mustBeNumeric, mustBeFloat} = [] % reference data
    kwargs.tdim (1,1) int64 {mustBeNumeric, mustBeInteger, mustBePositive} = 1
    kwargs.ndim (1,1) int64 {mustBeNumeric, mustBeInteger, mustBePositive} = 2
    kwargs.ldim (1,1) int64 {mustBeNumeric, mustBeInteger, mustBePositive} = ndims(x) + 1;
    kwargs.multi (1,1) logical = true; % whether to reference multiple channels
    kwargs.lvec (1,1) logical = true; % whether to vectorize lags
    kwargs.iflt (1,1) logical = exist('imfilter', 'file'); % requires image processing toolbox
end
% alias
x0 = kwargs.x0; 

% verify that w/W operates in the correct dimensions
wsz = size(W,1:max(kwargs.tdim, kwargs.ndim));
wsz_chk = wsz; wsz_chk([kwargs.tdim, kwargs.ndim]) = []; % size of W except in tdim,ndim
if any(wsz_chk ~= 1)
    error("QUPS:pwznxcorr:incompatibleWeightSize", ...
        "The filter weights w must be scalar in all dimensions except time (" ...
        +kwargs.tdim+") and channel ("+kwargs.ndim+")." ...
        );
end

% parse lags
isint = all(lags == floor(lags),'all');
if isscalar(lags) && isint, lags = -lags:lags; end

% parse window
if isscalar(W)
    wsz = ones(1,ndims(x));
    wsz(kwargs.tdim) = W; % window size is 1 x ... x W 
    w = ones(wsz, 'like', x); % averaging window (no need to scale)
else % W is the window itself
    [w, W] = deal(W, numel(W)); %#ok<ASGLU> (for debug)
end

% window sizing
[C, Wn] = deal(1, wsz(kwargs.ndim)); % number of channels: to correlate | normalization
if kwargs.multi, [C, Wn] = deal(Wn, C); end 

% vector accumulation function
if kwargs.iflt % whether to use imfilter
    if kwargs.pad, iopt = 0; else, iopt = 'circular'; end % 0-padding
    kernfun = @(x) imfilter(x, double(w), 'same', iopt); 
else
    kernfun = @(x) cast(convn(x, w, 'same'), 'like', x);
end

% lag-sampling function
if isint
    lsampfun = @(x, l, dim) arrayfun(@(l){circshift(x, -l, dim)},l);
else
    lsampfun = @(x, l, dim) {interpd(x, l+swapdim(0:size(x,dim)-1,2,dim), dim)};
end

% manual temporal 0-padding (if imfilter not available)
if kwargs.pad && ~kwargs.iflt
    P = ceil(max(abs(lags))); % padding length
    D = max([ndims(x), kwargs.tdim, kwargs.ndim, kwargs.ldim]); % max dimension
    psz = size(x,1:D); % data dimensions
    psz(kwargs.tdim) = P; % pad in this dimension
    x = cat(kwargs.tdim, x, zeros(psz, 'like', x)); % 0-pad at end
    psz = size(x0,1:D); % data dimensions
    psz(kwargs.tdim) = P; % pad in this dimension
    x0 = cat(kwargs.tdim, x0, zeros(psz, 'like', x)); % 0-pad at end
end

% chose which method selects the reference channel for correlation
N = size(x, kwargs.ndim); % number of channels to correlate
switch kwargs.ref
    case "neighbor"
        % specify left/right channels
        S = kwargs.stride;
        xl = sub(x, 1:N-S, kwargs.ndim); % (T x N-S x ...)
        xr = sub(x, 1+S:N, kwargs.ndim); % (T x N-S x ...)

    case "center"
        % find the middle channel
        mid = (N + 1 - C + 1)/2 + (0:C-1);
        n = unique([floor(mid), ceil(mid)]); % mean signal approach
        % n = floor(mid); % pick one signal, although spatially biased
        xl = x; % (T x N x ...)
        xr = sub(x, n, kwargs.ndim);
        if ~kwargs.multi, xr = mean(xr, kwargs.ndim); end % (T x 1 x ...)
        
    case "x0"
        % TODO: validate data size / cast type etc.
        % use a reference channel to correlate against
        xl = x; % (T x [1|N] x ...)
        xr = x0; % ([1|T] x N x ...)
end

% splice 
CEN  = kwargs.zero; % 0-mean
NORM = kwargs.norm; % normalization
tdim = kwargs.tdim; % operation dimension

% upsample
if U > 1, [xr, tr] = resample(xr, 1:size(xr,tdim), U, 'spline', 'Dimension', tdim); end
T = size(xr, tdim);  % time samples
assert(numel(1:U:T) == size(xl, tdim), "Sizing mismatch - this is a bug.");

% normalization (with centering) for the left channels
if CEN,  xlz = xl - kernfun(xl); else,  xlz = xl; end 
if NORM, xln = kernfun(real(xlz .* conj(xlz))); end % normalization power

% windowed inner product for each lag
% function y = windowed_inner_lag(lag)
for i = numel(lags) : -1 : 1
    % delay right channel by lag
    if kwargs.lvec
        xr_l = lsampfun(xr, swapdim(lags,2,kwargs.ldim), tdim);
        xr_l = conj(cat(kwargs.ldim, xr_l{:}));
    else
        xr_l = conj(subsref(lsampfun(xr, lags(i), tdim),substruct('{}',{1})));
    end

    % select portion
    if U > 1, xr_l = sub(xr_l, 1:U:T, tdim); end

    % mean centering for Z(N)CC
    if CEN, xrz_l = xr_l - kernfun(xr_l); %   0-centered
    else,   xrz_l = xr_l;                 % non-centered
    end

    % non-normalized correlation
    if kwargs.multi
        y = kernfun(convd(xlz, flip(xrz_l,2), 2, "same")); % (T x [N|N-S] x ...)
    else
        y = kernfun(xlz .* xrz_l); % (T x [N|N-S] x ...)
    end

    % normalization for (Z)NCC
    if NORM
        % normalization for the lagged right channels
        xrn_l = kernfun(real(xrz_l .* conj(xrz_l)));

        % normalization denominator - complex for (garbage) negative amplitudes
        if kwargs.multi
            r = sqrt(complex(convd(xln, flip(xrn_l,2), 2, "same"))) .* sqrt(Wn); % denom
        else
            r = sqrt(complex(xln)) .* sqrt(complex(xrn_l)) .* sqrt(Wn); % denom
        end

        % normalize
        y = y ./ r;
    end

    % store
    yi{i} = y;

    % short-circuit
    if kwargs.lvec, break; end
end

% compute for all lags and unpack
% y = arrayfun(@windowed_inner_lag, lags, 'UniformOutput', false);
if ~kwargs.lvec, y = cast(cat(kwargs.ldim, yi{:}), 'like', x); end % (T x [N|N-S] x ... x lags)

% undo manual 0-padding
if kwargs.pad && ~kwargs.iflt, y = sub(y, 1:size(y,kwargs.tdim)-P, kwargs.tdim); end


