
I = 64;
T = 32;
N = 4;
M = 3;
F = 2;
t = (0:T-1)'/T;
n = (0:N-1);
m = shiftdim((0:M-1),-2);
f = shiftdim((0:F-1),-1);

wm = 1/2 + (m/(M-1));
wf = 1/2 + (f/(F-1));
x0 = exp(2j*pi*(1/2+f/2.*n/4).*t);
tau = 2 + (T-4) * randn([T,N,F,1]);
tau = 2 + (T-4) * ((1:I)'./I .* (1+n)./N .* (1+m)./M);

w = 1;
% w = w .* wm;
w = w .* wf;

clear y ;
for dsum = {[], [2], [3,4], 6}
for wvec = {1, wf, wm .* wf} % different outer-product weights
for type = ["double", "single", "halfT"] % different precisions
for slicedim = [8]; % [8, 4, 3, 2]
for i = 1:2 % gpu/cpu
    switch i,case 1, x = gather(x0); case 2, x = gpuArray(x0); end
    x = sub(x, 1, slicedim);
    w = sub(wvec{1}, 1, slicedim);
    

    [x, tau, w] = dealfun(str2func(type), x, tau, w);
    y{i} = wsinterpd(x, tau, 1, w, dsum{1}, 'nearest', 0);

    
%     figure;
%     subplot(2,1,1); plot(real(gather(y(:,:))), '.-');
%     subplot(2,1,2); plot(imag(gather(y(:,:))), '.-');
end
assert(isalmostn(y{2}, y{1}));
end
end
end
end