classdef (TestTags = ["Github", "full", "build"]) InitTest < matlab.unittest.TestCase
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
            xdcs = [TransducerArray, TransducerConvex, TransducerMatrix, TransducerGeneric]; % init / hetero
            arrayfun(@(xdc) test.assertThat(xdc, IsInstanceOf('Transducer')), xdcs); % class maintained
            arrayfun(@(xdc) scale(xdc, 'dist', 1e3, 'time', 1e6), xdcs, "UniformOutput",false); % can scale
            xdcs(end+2) = xdcs(1); % implicit empty value            
       end
        function initseq(test)
            % INITSEQ - Assert that Sequence constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            seqs = [Sequence(), SequenceRadial(), SequenceGeneric()];
            arrayfun(@(seq) test.assertThat(seq, IsInstanceOf('Sequence')), seqs);
            arrayfun(@(scn) scale(scn, 'dist', 1e3), seqs, "UniformOutput",false); % can scale
            seqs(end+2) = seqs(1); % implicit empty value 
        end
        function initscan(test)
            % INITSCAN - Assert that Scan constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            scns = [ScanCartesian(), ScanPolar(), ScanSpherical(),ScanGeneric()]; % init / hetero
            arrayfun(@(scn) test.assertThat(scn, IsInstanceOf('Scan')), scns);
            arrayfun(@(scn) scale(scn, 'dist', 1e3), scns, "UniformOutput",false); % can scale
            scns(end+2) = scns(1); % implicit empty value            
        end
        function initchd(test)
            % INITCHD - Assert that a ChannelData constructor initializes
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            chd = ChannelData();
            test.assertThat(chd, IsInstanceOf('ChannelData'));
            scale(chd, 'time', 1e6);
        end
        function initus(test)
            % INITUS - Assert that an UltrasoundSystem constructor
            % initializes without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            us = UltrasoundSystem();
            test.assertThat(us, IsInstanceOf('UltrasoundSystem'));
            us = scale(us, "dist",1e3, "time",1e6);
        end
        function initmedscat(test)
            % INITMEDSCAT - Assert that a Scatterers and Mediums
            % initialize without arguments

            import matlab.unittest.constraints.IsInstanceOf;
            test.assertThat(Scatterers(), IsInstanceOf('Scatterers'));
            test.assertThat(Medium(),     IsInstanceOf('Medium'));
            scale(Medium()    , "dist", 1e3, "time", 1e6, "mass", 1e6);
            scale(Scatterers(), "dist", 1e3, "time", 1e6);

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