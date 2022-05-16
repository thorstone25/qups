classdef InitTest < matlab.unittest.TestCase
    % INITTEST - Initialization tests class
    %
    % This class test that all objects initialize properly

    properties(ClassSetupParameter)
        par   = struct('no_parpool', {{}}, 'def_parpool', {{'parallel'}});
        cache = struct('no_cache', {{}}, 'cache', {{'cache'}});
        gpu   = struct('no_ptx', {{}}, 'ptx', {{'CUDA'}});
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
        % Shared setup for the entire test class
        function setupQUPS(test, par, cache, gpu)
            cd(InitTest.proj_folder); % setup relative to here
            setup(par{:}, gpu{:}, cache{:});
        end
    end
    methods(TestClassTeardown)
        function teardownQUPS(test)
            cd(InitTest.proj_folder); 
            teardown; % basic teardown should run
            try delete(gcp('nocreate')); end % undo parpool if exists
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end
    methods(Test)
        function initxdc(test)
            % INITXDC - Assert that Transducer constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(TransducerArray() , IsInstanceOf('Transducer'));
            test.assertThat(TransducerConvex(), IsInstanceOf('Transducer'));
        end
        function initchd(test)
            % INITCHD - Assert that a ChannelData constructor initializes
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(ChannelData(), IsInstanceOf('ChannelData'));
        end
        function initseq(test)
            % INITSEQ - Assert that Sequence constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(Sequence(), IsInstanceOf('Sequence'));
            test.assertThat(SequenceRadial(), IsInstanceOf('Sequence'));
        end
        function initscan(test)
            % INITSCAN - Assert that Scan constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(ScanCartesian(), IsInstanceOf('Scan'));
            test.assertThat(ScanPolar(), IsInstanceOf('Scan'));
        end
        function initus(test)
            % INITUS - Assert that an UltrasoundSystem constructor
            % initializes without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(UltrasoundSystem(), IsInstanceOf('UltrasoundSystem'));
        end
    end

    methods(Static)
        % PROJ_FOLDER - Identifies the base folder for the project
        function f = proj_folder(), f = fullfile(fileparts(mfilename('fullpath')), '..'); end
    end
end