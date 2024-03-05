classdef (TestTags = ["Github", "full"])InitTest < matlab.unittest.TestCase
    % INITTEST - Initialization tests class
    %
    % This class test that all objects initialize properly

    properties(TestParameter)
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
    end
    methods(TestClassTeardown)
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
            test.assertThat(TransducerMatrix(), IsInstanceOf('Transducer'));
            test.assertThat(TransducerGeneric(), IsInstanceOf('Transducer'));
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
            test.assertThat(SequenceGeneric(), IsInstanceOf('Sequence'));
        end
        function initscan(test)
            % INITSCAN - Assert that Scan constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(ScanCartesian(), IsInstanceOf('Scan'));
            test.assertThat(ScanPolar(), IsInstanceOf('Scan'));
            test.assertThat(ScanSpherical(), IsInstanceOf('Scan'));
            test.assertThat(ScanGeneric(), IsInstanceOf('Scan'));
        end
        function initus(test)
            % INITUS - Assert that an UltrasoundSystem constructor
            % initializes without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(UltrasoundSystem(), IsInstanceOf('UltrasoundSystem'));
        end
        function initmedscat(test)
            % INITMEDSCAT - Assert that a Scatterers and Mediums
            % initialize without arguments

            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(Scatterers(), IsInstanceOf('Scatterers'));
            test.assertThat(Medium(), IsInstanceOf('Medium'));
        end
        function staticTest(test)
            % static constructors
            cls = [
                "Transducer" + ["Array","Convex", "Matrix", "Generic"], ...
                "Sequence" + ["", "Radial", "Generic"], ...
                "Scan" + ["Cartesian", "Polar", "Spherical", "Generic"], ...
                "Scatterers", "Medium", "Waveform", ...
                "UltrasoundSystem", "ChannelData", ...
                ];

            % all class escriptions
            mc = arrayfun(@meta.class.fromName, cls);
            for i = 1:numel(mc)
                mthd = mc(i).MethodList([mc(i).MethodList.Static]);
                n = arrayfun(@(m)length(m.InputNames), mthd);
                mthd = mthd(n == 0 & {mthd.Access}' == "public"); % 0 inputs only TODO: handle varargin?
                com = "@()"+cls(i) + "." + {mthd.Name}' + "()";
                arrayfun(@(s) test.assertWarningFree(str2func(s)), com);
            end 
        end
    end
end