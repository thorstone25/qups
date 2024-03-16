classdef KernTest < matlab.unittest.TestCase
    % KernTest - Function tests class
    %
    % This class test that all kernel functions operate properly

     properties(TestParameter)
        prec = {'double', 'single'}%, 'halfT'}
        complexity = struct('real', 'real', 'complex', 'complex');
        dev = getdevs()
     end

     properties
         dargs % data arguments
     end

    methods(TestMethodSetup)
        % Setup for each test
        function set_data(test)
            % create random values for interpd
            [I, T, N, M, F] = deal(2^5, 2^8, 2, 3, 4);
            fc = shiftdim((0 : F - 1), -2); % 1 x 1 x 1 x F
            t = (0 : T - 1)' ./ T; % T x 1 x 1 x 1
            x = (cospi(2*fc.*t) + 1i*sinpi(2*fc.*t)) + 0.01*(complex(rand([1,N]), rand([1,N])) - (0.5 + 0.5i)); % T x N x 1 x F
            tau = (4 + (T - 8)) .* rand([I, N, M, 1]); % I x N x M x 1
            test.dargs{1} = {I, T, N, M, F, fc, t, x, tau};

            % random values for convd
            [N,M] = deal(2^10, 1);
            A = complex(randn([N,M]), randn([N,M]));
            B = complex(randn([N,M]), randn([N,M]));
            test.dargs{2} = {A, B, N, M};

            % random values for slsc / pcf / dmas
            us = UltrasoundSystem();
            [us.xdc.numel, us.seq.numPulse] = deal(8);
            [us.scan.nx, us.scan.nz] = deal(32);
            sz = [us.scan.size us.xdc.numel];
            b = complex(randn(sz), randn(sz));
            test.dargs{3} = {us, b};
        end
    end

    % dispatch
    % Github test routine
    methods(Test, ParameterCombination = 'exhaustive', TestTags={'Github'})
        function github_runconvd(test, prec, complexity)
            runconvd(test, prec, complexity, "none")% forward remaing
        end
        function github_runinterpd(test, prec, complexity)
            runinterpd(test, prec, complexity, "none")
        end
        function github_apred(test, dev)
            aperture_reduction(test, dev, "single");
        end
    end

    % Full test routine
    methods(Test, ParameterCombination = 'exhaustive', TestTags={'full', 'build'})
        function full_convd(test, prec, complexity, dev), runconvd(test, prec, complexity, dev); end % forward all
        function full_interpd(test, prec, complexity, dev), runinterpd(test, prec, complexity, dev); end % forward all
        function full_apred(test, prec, dev), aperture_reduction(test, dev, prec); end % forward all
    end


    % test implementations
    methods
        function runconvd(test, prec, complexity, dev)

            % import matlab.unittest.TestCase;
            import matlab.unittest.constraints.NumericComparator;
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            if prec == "halfT", test.assumeTrue(logical(exist('halfT', 'class'))); end

            % random values
            [A, B, N, M] = deal(test.dargs{2}{:});
            [A,B] = dealfun(str2func(prec), A, B); % cast to type
            switch prec
                case "double", tol = 1e14;
                case "single", tol = 1e4; 
                case "halfT" , tol = 1e2;
            end
            switch complexity, case "real", [A,B] = deal(real(A), real(B)); end

            % check each type of argument
            % TODO: move this to setup params
            shapes = ["full", "same", "valid"];
            dims = [1,2,3]; 

            % convolve and check
            for s = shapes, for d = dims %#ok<ALIGN> 
                z1 = cell2mat(cellfun(@(A,B) {conv(A, B, char(s))}, num2cell(A,1), num2cell(B,1)));
                z2 = gather(convd(shiftdim(A,1-d),shiftdim(B,1-d),d,s,'gpu',dev == "gpu", 'ocl', dev == "ocl"));
                test.assertThat(...
                    double(shiftdim(z2,d-1)), ...
                    IsEqualTo(double(z1), 'Using', NumericComparator('Within', RelativeTolerance(tol*double(eps(max(z1(:))))))) ...
                    );
            end,end
        end
        function runinterpd(test, prec, complexity, dev)
            % test the interpd fuction behaves like interp1 but better
            import matlab.unittest.constraints.NumericComparator;
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;
            if prec == "halfT", test.assumeTrue(logical(exist('halfT', 'class'))); end
            if prec == "halfT", test.assumeFalse(dev == "ocl"); end % not supported

            [I, T, N, M, F, fc, t, x, tau] = deal(test.dargs{1}{:});
            if dev == "gpu", [x, t, tau] = dealfun(@gpuArray, x, t, tau); end

            % all permutations
            terp = ["nearest", "linear"];%, "cubic"];
            ords = [1,2,3,4; 1,2,4,3; 2,1,3,4; 3,4,1,2; 4,3,2,1];
            switch prec
                case "double", tol = 1e12;
                case "single", tol = 1e5;
                case "halfT" , tol = 1e4;
            end
            switch complexity, case "real", x = real(x); end

            % cast data
            for terp_ = terp
                [x,t,tau] = dealfun(str2func(prec), x,t,tau);
                % get matching data via interp1
                clear z0;
                for f = F:-1:1, for m = M:-1:1 %#ok<ALIGN> 
                        xf = num2cell(x  (:,:,:,f),1);
                        tf = num2cell(tau(:,:,m,:),1);
                        z0(1,:,m,f) = cellfun(@(x,t) {interp1(x,1+t,terp_,0)}, xf, tf);
                end, end
            z0 = reshape(cat(2,z0{:}), [I,N,M,F]); % I x N x M x F

            for ord = ords'
                [xp, tp] = dealfun(@(x)permute(x, ord), x, tau);
                z1 = wsinterpd(xp, tp, find(ord == 1), 1, [], terp_, 0);
                z1 = ipermute(z1, ord);
                test.assertThat(...
                    double(gather(z1)), IsEqualTo(double(gather(z0)), 'Using', ...
                    NumericComparator('Within', RelativeTolerance(tol*gather(double(eps(cast(1,'like',z0)))))) ...
                    ));
            end
            end
        end
        function aperture_reduction(test, dev, prec)
            if prec == "halfT", return; end % not supported
            [us, b] = deal(test.dargs{3}{:});
            f = str2func(prec); b = f(b);
            if dev == "gpu", b = gpuArray(b); end
            z = slsc(b, 'ocl', dev == "ocl");
            w = pcf(b);
            m = dmas(b);
            r = cohfac(b);
            for d = [2 4]
                z = slsc(b, d, [0 2], "ensemble", 'ocl', dev == "ocl");
                w = pcf(b, d);
                m = dmas(b, d);
                r = cohfac(b, d);
            end
            
            if dev ~= "ocl"
                z = slsc(repmat(b,[ones(1,6),2]), 4, 2, "ensemble", 7);
            end
        end
    end

    methods(Static)
        % PROJ_FOLDER - Identifies the base folder for the project
        function f = proj_folder(), f = fullfile(fileparts(mfilename('fullpath')), '..'); end
    end
end

% test GPU devices if we can
function s = getdevs()
s = "none";
if gpuDeviceCount, s(end+1) = "gpu"; end
if exist('oclDeviceCount', 'file') && oclDeviceCount()
    devs = oclDeviceTable();
    devs = sortrows(devs, ["Type", "MaxComputeUnits"]);
    oclDevice(devs{end,"Index"}); % select best device
    s(end+1) = "ocl";
end
s = cellstr(s);
end
