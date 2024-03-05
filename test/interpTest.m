classdef(TestTags = ["Github", "full"]) interpTest < matlab.unittest.TestCase
    % KERNTEST - Test the compute kernels within QUPS
    %
    % This class test the various compute kernels within QUPS
    properties
        x0
        tau
        t1
        t2
        inds
    end

    properties(ClassSetupParameter)
        dsize = {[16 32 4 3 2]} % data test set size
        type = {'double', 'single'}%, 'halfT'} % different precisions
        dev = {'CPU', 'GPU'} % gpu/cpu
    end
    methods(TestClassSetup, ParameterCombination = 'exhaustive')
        function setupData(test, dsize,type, dev)

            if type == "halfT"
                test.assumeTrue(logical(exist('halfT', 'class')));
            end
    	    if dev == "GPU"
        		test.assumeTrue(gpuDeviceCount > 0);
    	    end

            % get sizing
            tmp = num2cell(dsize);
            [I,T,N,M,F] = deal(tmp{:});
            
            % data vectors
            i = (0:I-1)';
            t = (0:T-1)'/T;
            n = (0:N-1);
            m = shiftdim((0:M-1),-1);
            f = shiftdim((0:F-1),-2);
            t1_ = 4 + (T-8) * ((1+i)./I .*   1./N    .* (1+m)./M); % I x 1 x M x 1
            t2_ = 4 + (T-8) *   1./I    .* (1+n)./N .*      1    ; % 1 x N x 1 x 1
            tau_ = t1_ + t2_; % I x N x M x 1

            % create data and sampling matrices
            x0_ = exp(2j*pi*(1/2+f/2.*n/4).*t); % T x N x 1 x F

            cfun = str2func(type);
            x0_  = cfun(x0_ );
            tau_ = cfun(tau_);

            % gpu/cpu
            switch dev, 
                case "CPU", x0_ = gather(  x0_); 
                case "GPU", x0_ = gpuArray(x0_);
            end

            [test.x0, test.tau, test.t1, test.t2, test.inds] = deal( ...
                x0_, tau_, t1_, t2_,{i,t,n,m,f});

            % cast data type
        end
    end
    methods(TestClassTeardown)
    end

    methods(TestMethodSetup)
        % Setup for each test
    end
    properties(TestParameter)
        % all permutations to test
        % ord = {[1,2,3,4], [1,3,2,4], [1,4,2,3]} % permute the data?
        terp = {'cubic', 'nearest', 'linear'}; % overlapping methods
        dsum = {[], [2], [3,4], 6} % summation dimensions
        wvecd = {[], 3, [3,4]} % different outer-product weights dimensions
    end
    methods(Test, ParameterCombination = 'exhaustive')
        function wsinterpdTest(test,terp,dsum,wvecd)
            [x0_, tau_, t1_, t2_] = deal(test.x0, test.tau, test.t1, test.t2);
            [i,t,n,m,f] = test.inds{:};
            type_ = class(gather(x0_([])));
            cfun = str2func(type_);

            % expand w as requested
            wsz = ones([1 max(wvecd)]);
            if ~isempty(wvecd), wsz(wvecd) = max(size(x0_,wvecd), size(tau_, wvecd)); end
            w = cfun(rand(wsz));

            % compute split and joint
            y1 = wsinterpd( x0_,   tau_,   1, w, dsum, terp, 0);
            y2 = wsinterpd2(x0_, t1_, t2_, 1, w, dsum, terp, 0);

            % use single precision for comparing to half values
            if type_ == "halfT", [x0_,tau_] = dealfun(@(x)single(gather(x)), x0_,tau_); end

            % use a (slow) for loop
            y0 = zeros(size(i+n+m+f)) .* y1(1); %#ok<SZARLOG> 
            for fi = 1+flip(f(:)'), for mi = 1+flip(m(:)), for ni = 1+flip(n(:))'
                y0(:,ni,mi,fi) = interp1(x0_(:,ni,1,fi), 1+tau_(:,ni,mi,1), terp, 0);
            end, end, end
            if isempty(dsum), dsum = ndims(y0) + 1; end % avoid summing empty dimensions
            y0 = sum(w .* y0, dsum);

            % compare answers
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.IsEqualTo;
            tol = 1e2; 
            if terp == "cubic"
            switch type_ % needs much more room because the algs aren't the same
                case "double", tol = tol * 1e16; 
                case "single", tol = tol * 1e8;
                case "halfT",  tol = tol * 1;
            end
            end

            switch type_
                case "halfT"
                    test.assertThat(double(y1), IsEqualTo(double(y0), 'Within', RelativeTolerance(tol*double(max(eps(y0(:)))))))
                    test.assertThat(double(y2), IsEqualTo(double(y0), 'Within', RelativeTolerance(tol*double(max(eps(y0(:)))))))
                otherwise
                    test.assertThat(y1, IsEqualTo(y0, 'Within', RelativeTolerance(tol*(max(eps(y0(:)))))))
                    test.assertThat(y2, IsEqualTo(y0, 'Within', RelativeTolerance(tol*(max(eps(y0(:)))))))
            end
        end
    end
end
