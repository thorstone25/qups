classdef (TestTags = ["Github", "full"])InitTest < matlab.unittest.TestCase
    % INITTEST - Initialization tests class
    %
    % This class test that all objects initialize properly

    properties(TestParameter)
        par   = struct('no_parpool', {{}}, 'def_parpool', {{'parallel'}});
        cache = struct('no_cache', {{}}, 'cache', {{'cache'}});
        gpu   = struct('no_ptx', {{}}, 'ptx', {{'CUDA'}});
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
        % Shared setup for the entire test class
        function setupQUPS(test)
            cd(InitTest.proj_folder); % setup relative to here
            setup; % setup alone to add paths
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
        function initQUPS(test, par, cache, gpu)
            cd(InitTest.proj_folder); % setup relative to here
            setup(par{:}, gpu{:}, cache{:}); % setup with each option combo should work
        end

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
        function inittarg(test)
            % INITTARG - Assert that a Target and it's subclasses
            % initialize without arguments
            
            % TODO: this will be reorganized so the Medium/Scatterers are
            % subclasses of Target, at which point the test must be changed
            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(Scatterers(), IsInstanceOf('Scatterers'));
            test.assertThat(Medium(), IsInstanceOf('Medium'));
            test.assertThat(Target(), IsInstanceOf('Target'));
        end
    end

    methods(Static)
        % PROJ_FOLDER - Identifies the base folder for the project
        function f = proj_folder(), f = fullfile(fileparts(mfilename('fullpath')), '..'); end
    end
end