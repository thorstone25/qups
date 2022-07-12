classdef(TestTags = ["Github", "full"]) interpTest < matlab.unittest.TestCase
    % KERNTEST - Test the compte kernels within QUPS
    %
    % This class test the various compute kernels within QUPS

    properties(ClassSetupParameter)
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
        % Shared setup for the entire test class
        function setup(test), 
            addpath(fullfile(interpTest.base_dir(), 'kern'));
            addpath(fullfile(interpTest.base_dir(), 'src'));
            addpath(fullfile(interpTest.base_dir(), 'utils'));
        end
    end
    methods(TestClassTeardown)
        function teardown(test), 
            rmpath(fullfile(interpTest.base_dir(), 'kern'));
            rmpath(fullfile(interpTest.base_dir(), 'src'));
            rmpath(fullfile(interpTest.base_dir(), 'utils'));
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end
    properties(TestParameter)
        % all permutations to test
        % ord = {[1,2,3,4], [1,3,2,4], [1,4,2,3]} % permute the data?
        dsize = {[16 32 4 3 2]} % data test set size
        terp = {'cubic', 'nearest', 'linear'}; % overlapping methods
        dsum = {[], [2], [3,4], 6} % summation dimensions
        wvecd = {[], 3, [3,4]} % different outer-product weights dimensions
        type = {"double", "single"}%, "halfT"} % different precisions
        dev = {'CPU'}%, 'GPU'} % gpu/cpu
    end
    methods(Test, ParameterCombination = 'exhaustive')
        function wsinterpdTest(test,dsize,terp,dsum,wvecd,type,dev)

            % get sizing
            tmp = num2cell(dsize);
            [I,T,N,M,F] = deal(tmp{:});
            
            % data vectors
            i = (0:I-1)';
            t = (0:T-1)'/T;
            n = (0:N-1);
            m = shiftdim((0:M-1),-1);
            f = shiftdim((0:F-1),-2);

            % create data and sampling matrices
            x0 = exp(2j*pi*(1/2+f/2.*n/4).*t); % T x N x 1 x F
            tau = 4 + (T-8) * ((1+i)./I .* (1+n)./N .* (1+m)./M); % I x N x M x 1

            % expand w as requested
            w = 1; for d = wvecd, w = w + shiftdim(rand([max(size(x0,d),size(t,d)),1]),1-d); end

            % gpu/cpu
            switch dev, 
                case "CPU", x = gather(x0); 
                case "GPU", x = gpuArray(x0);
            end

            % cast data type
            [x, tau, w] = dealfun(str2func(type), x, tau, w);

            % compute using optimized routine
            y = wsinterpd(x, tau, 1, w, dsum, terp, 0);

            % use single precision for comparing to half values
            if type == "halfT", [x,tau] = dealfun(@(x)single(gather(x)), x,tau); end

            % use a (slow) for loop
            y0 = zeros(size(i+n+m+f)) .* y(1);
            for fi = 1+flip(f(:)'), for mi = 1+flip(m(:)), for ni = 1+flip(n(:))'
                y0(:,ni,mi,fi) = interp1(x(:,ni,1,fi), 1+tau(:,ni,mi,1), terp, 0);
            end, end, end
            if isempty(dsum), dsum = ndims(y0) + 1; end % avoid summing empty dimensions
            y0 = sum(w .* y0, dsum);


            % compare answers
            import matlab.unittest.constraints.RelativeTolerance;
            import matlab.unittest.constraints.IsEqualTo;
            tol = 1e2; 
            if terp == "cubic"
            switch type % needs much more room because the algs aren't the same
                case "double", tol = tol * 1e10; 
                case "single", tol = tol * 1e2;
                case "halfT",  tol = tol * 1;
            end
            end
            switch type
                case "halfT"
                    test.assertThat(double(y), IsEqualTo(double(y0), 'Within', RelativeTolerance(tol*double(max(eps(y0(:)))))))
                otherwise
                    test.assertThat(y, IsEqualTo(y0, 'Within', RelativeTolerance(tol*(max(eps(y0(:)))))))
            end
            

        end
    end
    methods(Static)
        function d = base_dir(), [d] = fullfile(fileparts(mfilename('fullpath')), '..'); end
    end
end