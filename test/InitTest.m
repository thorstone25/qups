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

            % deprecation
            seq = Sequence();
            test.assertWarning(@()setfield(seq,'delays_',     zeros(4)), "QUPS:Sequence:DeprecatedProperty")
            test.assertWarning(@()setfield(seq,'apodization_',zeros(4)), "QUPS:Sequence:DeprecatedProperty")
            test.assertWarning(@()seq.delays_                          , "QUPS:Sequence:DeprecatedProperty")
            test.assertWarning(@()seq.apodization_                     , "QUPS:Sequence:DeprecatedProperty")
            
            % on construction only?
            test.assertWarning(@()Sequence('type','VS')                , "QUPS:Sequence:DeprecatedValue"   )
            
        end
        function initscan(test)
            % INITSCAN - Assert that Scan constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            scns = [ScanCartesian(), ScanPolar(), ScanSpherical(), ScanGeneric()]; % init / hetero
            arrayfun(@(scn) test.assertThat(scn, IsInstanceOf('Scan')), scns);
            arrayfun(@(scn) scale(scn, 'dist', 1e3), scns, "UniformOutput",false); % can scale
            scns(end+2) = scns(1); % implicit empty value

            % test assigning / retrieving sizing
            dvars = 'XYZRUVW'; % dist
            avars = 'AE'; % ang
            for scn = scns
                ax = lower(scn.order);
                for a = ax, scn.("n"+a) =   4; end % dist or ang
                test.assertTrue(all(scn.size == 4), "Setting the scan size failed for a " + class(scn) + "!");
                for a = ax
                    if a == 'r',                      bd = [  0 40]; % range
                    elseif ismember(upper(a), dvars), bd = [-20 20]; % dist
                    elseif ismember(upper(a), avars), bd = [ -5  5]; % ang
                    else,                             bd = [ -5  5]; 
                    end
                    scn.(a+"b") = bd;
                    test.assertTrue(isalmostn(scn.(a+"b"), bd), "Setting the scan boundaries failed for a " + class(scn) + "!");
                end
                for a = ax, scn.("d"+a) = 1/2; end % dist or ang
                test.assertTrue(all(arrayfun(@(c) scn.("d"+c),ax) == 1/2), "Setting the scan resolution failed for a " + class(scn) + "!");
            end



        end
        function initchd(test)
            % INITCHD - Assert that a ChannelData constructor initializes
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            chd = ChannelData();
            test.assertThat(chd, IsInstanceOf('ChannelData'));
            scale(chd, 'time', 1e6);
            [chd.tdim, chd.ndim, chd.mdim, chd.T, chd.N, chd.M]; % get

            % deprecated properties
            test.assertWarning(@()setfield(chd,'ord','TMN'), "QUPS:ChannelData:syntaxDeprecated")
            test.assertWarning(@()chd.ord                  , "QUPS:ChannelData:syntaxDeprecated")
        end
        function initus(test)
            % INITUS - Assert that an UltrasoundSystem constructor
            % initializes without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            us = UltrasoundSystem();
            test.assertThat(us, IsInstanceOf('UltrasoundSystem'));
            us = scale(us, "dist",1e3, "time",1e6);
            fld = us.tmp_folder;
            clear us; % destroy
            test.assertFalse(logical(exist(fld,"dir")));
            us = UltrasoundSystem('copybin',true);
            us = UltrasoundSystem('recompile',true);

            % deprecation
            test.assertWarning(@()setfield(us,'sequence',Sequence()),"QUPS:UltrasoundSystem:syntaxDeprecated")
            test.assertWarning(@()us.sequence                       ,"QUPS:UltrasoundSystem:syntaxDeprecated")

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