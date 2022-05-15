classdef InitTest < matlab.unittest.TestCase
    % INITTEST - Initialization tests class
    %
    % This class test that all objects initialize properly
    methods(TestClassSetup)
        % Shared setup for the entire test class
        function setupQUPS(test)
            cd(InitTest.proj_folder); % setup relative to here
            setup; % basic setup should be run
        end
    end
    methods(TestClassTeardown)
        function teardownQUPS(test)
            cd(InitTest.proj_folder); 
            teardown; % basic teardown should run
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    properties(TestParameter)
        par   = struct('none', {{}}, 'some', {{'parallel'}});
        cache = struct('none', {{}}, 'some', {{'cache'}});
        gpu   = struct('none', {{}}, 'some', {{'CUDA'}});
    end

    methods(Test, ParameterCombination = 'exhaustive')
        % Initialization test methods

        % TODO: write an actual test script for SETUP that
        % * accounts for when no parpool exists prior
        %   whether in teardown or repeated call
        % * tests that the proper folders are added to the path
        function initQUPS(test, par, gpu, cache)
            % INITQUPS - Assert that the QUPS setup scripts run
            cd(InitTest.proj_folder); % move to project folder
            setup(par{:}, gpu{:}, cache{:});
        end
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