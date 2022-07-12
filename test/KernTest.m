classdef KernTest < matlab.unittest.TestCase
    % KernTest - Function tests class
    %
    % This class test that all functions operate properly

    properties(ClassSetupParameter)
        gpu   = getgpu()
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
        % Shared setup for the entire test class
        function setupQUPS(test, gpu)
            cd(KernTest.proj_folder); % setup relative to here
            setup(gpu{:}); % setup with/without GPU copmilation support: each option combo should work
        end
    end
    methods(TestClassTeardown)
        function teardownQUPS(test)
            cd(KernTest.proj_folder); 
            teardown; % basic teardown should run
        end
    end
    properties(TestParameter)
        prec = {'double', 'single', 'halfT'}
        complexity = struct('real', 'real', 'complex', 'complex');
        gdev = getdevs() %struct('no_dev', 0); %, 'dev', -1');
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    % dispatch
    % Github test routine
    methods(Test, ParameterCombination = 'exhaustive', TestTags={'Github'})
        function github_runconvd(test, prec, complexity, gdev)
            switch prec, case {'halfT'}, return; end % not supported
            runconvd(test, prec, complexity, gdev)% forward remaing
        end
        function github_runinterpd(test, prec, complexity, gdev)
            switch prec, case {'halfT'}, return; end % not supported
            runinterpd(test, prec, complexity, gdev)
        end
    end

    % Full test routine
    methods(Test, ParameterCombination = 'exhaustive', TestTags={'full'})
        function full_convd(test, prec, complexity, gdev), runconvd(test, prec, complexity, gdev); end % forward all
        function full_interpd(test, prec, complexity, gdev), runinterpd(test, prec, complexity, gdev); end % forward all
    end


    % test implementations
    methods
        function runconvd(test, prec, complexity, gdev)

            % import matlab.unittest.TestCase;
            import matlab.unittest.constraints.NumericComparator;
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;

            % create random values
            [N,M] = deal(2^10, 1);
            A = complex(randn([N,M]), randn([N,M]));
            B = complex(randn([N,M]), randn([N,M]));
            [A,B] = dealfun(str2func(prec), A, B); % cast to type
            switch prec
                case "double", tol = 1e10;
                case "single", tol = 1e2; 
                case "halfT" , tol = 1e2;
            end
            switch complexity, case "real", [A,B] = deal(real(A), real(B)); end

            % check each type of argument
            % TODO: move this to setup params
            shapes = ["full", "same", "valid"];
            dims = [1,2,3]; 

            % convolve and check
            for s = shapes, for d = dims, %#ok<ALIGN> 
                z1 = cell2mat(cellfun(@(A,B) {conv(A, B, char(s))}, num2cell(A,1), num2cell(B,1)));
                z2 = convd(shiftdim(A,1-d),shiftdim(B,1-d),d,s,'device',gdev);
                test.assertThat(...
                    double(shiftdim(z2,d-1)), ...
                    IsEqualTo(double(z1), 'Using', NumericComparator('Within', RelativeTolerance(tol*double(eps(max(z1(:))))))) ...
                    );
            end,end
        end
        function runinterpd(test, prec, complexity, gdev)
            import matlab.unittest.constraints.NumericComparator;
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.RelativeTolerance;

            % test the interpd fuction behaves like interp1 but better
            % create random values
            [I, T, N, M, F] = deal(2^5, 2^8, 2, 3, 4);
            fc = shiftdim((0 : F - 1), -2); % 1 x 1 x 1 x F
            t = (0 : T - 1)' ./ T; % T x 1 x 1 x 1
            x = (cospi(2*fc.*t) + 1i*sinpi(2*fc.*t)) + 0.01*(complex(rand([1,N]), rand([1,N])) - (0.5 + 0.5i)); % T x N x 1 x F
            tau = (4 + (T - 8)) .* rand([I, N, M, 1]); % I x N x M x 1
            if gdev, [x, t, tau] = deal(gpuArray(x), gpuArray(t), gpuArray(tau)); end

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
                z1 = interpd(xp, tp, find(ord == 1), terp_, 0);
                z1 = ipermute(z1, ord);
                test.assertThat(...
                    double(gather(z1)), IsEqualTo(double(gather(z0)), 'Using', ...
                    NumericComparator('Within', RelativeTolerance(tol*gather(double(eps(cast(1,'like',z0)))))) ...
                    ));
            end
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
s.nodev = 0;
if gpuDeviceCount, s.dev = -1; end
end

function s = getgpu()
s.no_ptx = {};
if gpuDeviceCount, s.ptx = {'CUDA'}; end
end