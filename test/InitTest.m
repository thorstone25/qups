classdef (TestTags = ["Github", "full", "build", "syntax"]) InitTest < matlab.unittest.TestCase
    % INITTEST - Initialization tests class
    %
    % This class test that all objects initialize properly

    %#ok<*NASGU> unused variables

    properties(TestParameter)
        % all classes to test
        clss = cellstr([
            "Transducer" + ["Array","Convex", "Matrix", "Generic"], ...
            "Sequence" + ["", "Radial", "Generic"], ...
            "Scan" + ["Cartesian", "Polar", "Spherical", "Generic"], ...
            "Scatterers", "Medium", "Waveform", ...
            "UltrasoundSystem", "ChannelData", ...
            ]);
    end

    properties
    end

    methods(TestClassSetup, ParameterCombination = 'exhaustive')
    end
    methods(TestClassTeardown)
    end

    methods(TestMethodSetup)
        % Setup for each test
        function fig(~), figure; end
    end
    methods(TestMethodTeardown)
        function cls(~), close; end
    end
    methods(Test)
        function initxdc(test)
            % INITXDC - Assert that Transducer constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            xdcs = [TransducerArray, TransducerConvex, TransducerMatrix, TransducerGeneric]; % init / hetero
            arrayfun(@(xdc) test.assertThat(xdc, IsInstanceOf('Transducer')), xdcs); % class maintained
            arrayfun(@(xdc) scale(xdc, 'dist', 1e3, 'time', 1e6), xdcs, "UniformOutput",false); % can scale
            arrayfun(@obj2struct, xdcs, 'UniformOutput', false); % supports specialized struct conversion
            arrayfun(@plot, xdcs); % supports plotting
            arrayfun(@patch, xdcs); % supports patch
            arrayfun(@(x)patch(x, nexttile(), 'el_sub_div', [2 2]), xdcs); % supports args
            xdcs(end+2) = xdcs(1), %#ok<NOPRT> implicit empty value, display
            arrayfun(@disp, xdcs) % display scalar
            arrayfun(@(x)disp(x([])), xdcs) % display empty
            x = copy(xdcs); x.delete(); x, arrayfun(@disp, x); % display deleted
            test.assertWarning(@()xdcs(1).ultrasoundTransducerImpulse(), "QUPS:Transducer:DeprecatedMethod");
            [xdcs.origin] = deal([0 0 -10e-3]); % offset (to be deprecated?)
            [xdcs.rot] = deal([20 -10]); % offset (to be deprecated?)

        end
        function initseq(test)
            % INITSEQ - Assert that Sequence constructors initialize
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            typs = ["FSA","PW","FC","DV","VS"];
            seqs = [arrayfun(@(t) Sequence(      "type",t), typs(1:end)), ...
                arrayfun(@(t) SequenceRadial("type",t), typs(2:end)), ...
                SequenceGeneric(), ...
                SequenceGeneric('apd', ones(4), 'del', zeros(4)), ...
                ];
            [seqs([seqs.type] == "FSA").numPulse] = deal(1);
            arrayfun(@(seq) test.assertThat(seq, IsInstanceOf('Sequence')), seqs);
            arrayfun(@(scn) scale(scn, 'dist', 1e3), seqs, "UniformOutput",false); % can scale
            seqs(end+2) = seqs(1), %#ok<NOPRT> implicit empty value, display
            seqs(end-1).numPulse = 1; % fix implicit seq
            s = arrayfun(@obj2struct, seqs, 'UniformOutput', false); % supports specialized struct conversion
            arrayfun(@plot, seqs, 'UniformOutput',false); % supports plotting
            cellfun(@(s) test.assertThat(s.pulse, IsInstanceOf('struct')), s); % recursive check
            arrayfun(@disp, seqs) % display scalar
            arrayfun(@(x)disp(x([])), seqs) % display empty
            arrayfun(@splice, seqs, 'UniformOutput', false); % can splice
            x = copy(seqs); x.delete(); x, arrayfun(@disp, x); % display deleted

            % polar: manipulate range/angles
            seq = SequenceRadial("focus",[0 0 1]'); % fine
            seq = SequenceRadial("focus",[-0.5 0 0.5; 0 0 sqrt(2)/2; 0.5 0 0.5]'); % fine
            seq = SequenceRadial("angles",-1:1); % fine
            seq = SequenceRadial("ranges",[1 1 1]); % fine
            seq = SequenceRadial("angles",-1:1, "ranges",[1 1 1]); % fine
            seq = SequenceRadial("angles",   0, "ranges",[1 1 1]); % fine
            seq = SequenceRadial("angles",-1:1, "ranges",[1]); % fine
            seq.angles = [-5 0 5]; % fine
            seq.ranges = [2 3 2]; % fine
            test.assertError(@()SequenceRadial("focus",[0 0 1]',"angles",0),""); % bad
            test.assertError(@()SequenceRadial("focus",[0 0 1]',"ranges",1),""); % bad
            seq = SequenceRadial("angles",-1:1, "ranges",2*[1 1 1]); % fine
            seq.moveApex([0 0 -1]); % moving apex should not move range/angle
            test.assertTrue(isalmostn(seq.angles, -1:1))
            test.assertTrue(isalmostn(seq.ranges, 2*[1 1 1]))
            test.assertTrue(isalmostn(seq.apex, [0 0 -1]'))

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
            arrayfun(@obj2struct, scns, 'UniformOutput', false); % supports specialized struct conversion
            arrayfun(@plot, scns); % supports plotting
            arrayfun(@(s) imagesc(s,randn(s.size)), scns); % supports imagesc
            scns(end+2) = scns(1),%#ok<NOPRT> implicit empty value, display
            arrayfun(@disp, scns) % display scalar
            arrayfun(@(x)disp(x([])), scns) % display empty
            x = copy(scns); x.delete(); x, arrayfun(@disp, x); % display deleted

            % supports req'd overload methods
            [~,~,~] = arrayfun(@getImagingGrid, scns, "UniformOutput",false);
            [~]     = arrayfun(@positions     , scns, "UniformOutput",false);

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

            % convert polar to cart
            ScanCartesian(ScanPolar()),
            ScanCartesian(ScanSpherical()),

            % deprecations
            test.assertWarning(@()scanCartesian(ScanPolar()    ), "QUPS:ScanPolar:syntaxDeprecated"    );
            test.assertWarning(@()scanCartesian(ScanSpherical()), "QUPS:ScanSpherical:syntaxDeprecated");

        end
        function initwv(test)
            % INITWV - Assert that the Waveform constructor initializes
            import matlab.unittest.constraints.IsInstanceOf;
            wv = Waveform(); % init
            test.assertThat(wv, IsInstanceOf('Waveform'));

            % multiple ways to init
            dt = 1/8;
            fc = 1/2;
            t = 0 : dt : 1;
            f = @(t)sinpi(2*fc*t);
            x = f(t);
            wv = Waveform();
            wv = Waveform("samples", x, "t",  t);
            wv = Waveform("samples", x, "t0", t(1), "dt", dt  , "tend",t(end));
            wv = Waveform("samples", x, "t0", t(1), "fs", 1/dt, "tend",t(end));
            wv = Waveform("fun"    , f, "t",  t);
            wv = Waveform("fun"    , f, "t0", t(1), "dt", dt  , "tend",t(end));
            wv = Waveform("fun"    , f, "t0", t(1), "fs", 1/dt, "tend",t(end));

            scale(wv, 'time', 1e6); % supports scaling
            obj2struct(wv); % supports specialized struct conversion
            plot(wv); % supports plot
            plot(wv, "freqs", true); % supports frequency plot
        end
        function initchd(test)
            % INITCHD - Assert that a ChannelData constructor initializes
            % without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            chd = ChannelData(), %#ok<NOPRT> implicit display
            test.assertThat(chd, IsInstanceOf('ChannelData'));
            scale(chd, 'time', 1e6);
            arrayfun(@obj2struct, chd, 'UniformOutput', false); % supports specialized struct conversion
            arrayfun(@imagesc, chd); % supports imagesc

            % deprecated properties
            test.assertWarning(@()setfield(chd,'ord','TMN'), "QUPS:ChannelData:syntaxDeprecated")
            test.assertWarning(@()chd.ord                  , "QUPS:ChannelData:syntaxDeprecated")
        end
        function initus(test)
            % INITUS - Assert that an UltrasoundSystem constructor
            % initializes without arguments
            import matlab.unittest.constraints.IsInstanceOf;
            us = UltrasoundSystem(), %#ok<NOPRT> implicit display
            test.assertThat(us, IsInstanceOf('UltrasoundSystem'));

            fld = us.tmp_folder;
            test.assertTrue(logical(exist(fld,"dir")));
            UltrasoundSystem('copybin',  true);
            UltrasoundSystem('recompile',true);
            us2 = copy(us); fld2 = us2.tmp_folder; % new instance
            test.assertFalse(fld == fld2); % different folders
            clear us2; test.assertFalse(logical(exist(fld2,'dir'))); % deleted on clear
            us2 = copy(us); fld2 = us2.tmp_folder; % new instance
            copyfile(which("Qups.prj"), fld2); % dirty
            test.addTeardown(@()delete(fullfile(fld2,"Qups.prj"))); % cleanup
            try clear us2 % should error
                [~, lid] = lastwarn;
                test.assertTrue(lid == "MATLAB:class:DestructorError", ...
                    "Attempt to clear " + class(us) + " did not produce a warning." ...
                    );
                % test.assertFail("Attempt to clear " + class(us) + " did not produce an error.");
            catch ME
               test.assertTrue(ME.identifier == "QUPS:UltrasoundSystem:dirtyDirectory", ...
                   "Attempt to clear " + class(us) + " did not produce an error." ...
                   );
            end

            us.xdc; % should be same
            usb = copy(us); usb.rx = copy(usb.tx); % switch to bistatic aperture
            for us = [us usb]
                scale(us, "dist",1e3, "time",1e6);
                plot(us); % supports plotting
                s = obj2struct(us);
                flds = intersect(["tx","rx","xdc","seq","scan"], fieldnames(s)); % class properties
                arrayfun(@(p) test.assertThat(s.(p), IsInstanceOf('struct')), flds); % sub-class conversion worked
                test.assertSize(us.fc, [1 1], "The 'fc' property should be scalar for identical frequency transducers.")
            end
            test.assertError(@() usb.xdc, "QUPS:UltrasoundSystem:ambiguousTransducer")

            uss = copy([us us us]); % array
            uss(end+2) = copy(us); % implicit creation
            arrayfun(@disp, uss) % display scalar
            arrayfun(@(x)disp(x([])), uss) % display empty

            % deprecation
            test.assertWarning(@()setfield(us,'sequence',Sequence()),"QUPS:UltrasoundSystem:syntaxDeprecated")
            test.assertWarning(@()us.sequence                       ,"QUPS:UltrasoundSystem:syntaxDeprecated")

            % destroy
            clear us;
        end
        function initmedscat(test)
            % INITMEDSCAT - Assert that a Scatterers and Mediums
            % initialize without arguments

            import matlab.unittest.constraints.IsInstanceOf;
            sct = Scatterers(), %#ok<NOPRT> implicit display
            med = Medium(),     %#ok<NOPRT> implicit display
            test.assertThat(sct, IsInstanceOf('Scatterers'));
            test.assertThat(med, IsInstanceOf('Medium'    ));
            scale(sct, "dist", 1e3, "time", 1e6);
            scale(med, "dist", 1e3, "time", 1e6, "mass", 1e6);
            struct(med); % no-specialized struct conversion yet ...
            struct(sct); % specialized struct conversion
            imagesc(med, ScanCartesian()); % support image on cartesian grid
            plot(sct); % supports plot

            % special constructors
            grd = ScanCartesian();
            Medium.Sampled(grd)
            Scatterers.Diffuse(grd)
            Scatterers.Grid()

            % deprecated properties
            test.assertWarning(@()setfield(med,'alphap0',1.2), "QUPS:Medium:syntaxDeprecated")
            test.assertWarning(@()med.alphap0                , "QUPS:Medium:syntaxDeprecated")

        end
        function staticTest(test, clss)
            % staticTest - Test static functions with no inputs (constructors)
            mc = meta.class.fromName(clss); % all class descriptions
            mthd = mc.MethodList([mc.MethodList.Static]);
            n = arrayfun(@(m)length(m.InputNames), mthd);
            mthd = mthd(n == 0 & {mthd.Access}' == "public"); % 0 inputs only TODO: handle varargin?
            com = "@()"+clss + "." + {mthd.Name}' + "()";
            arrayfun(@(s) test.assertWarningFree(str2func(s)), com);
        end
        function propTest(test, clss)
            % Property set/get works as expected

            % all class descriptions
            mc = meta.class.fromName(clss);
            if(~mc.Abstract), return; end % ensure non-abstract classes only
            
            % instantiate the class
            x = eval(clss);

            % get all getable properties
            prp = mc.PropertyList({mc.PropertyList.GetAccess} == "public");

            % request all properties
            cellfun(@(f) {x.(f)}, {prp.Name});

            % get all setable properties
            prp = mc.PropertyList({mc.PropertyList.SetAccess} == "public");

            % set each property, somehow, in order to check that
            % (re)assignment works as expected
            for j = 1 : numel(prp)
                % property field name
                f = string(prp(j).Name);

                % assign via some strategy
                if prp(j).HasDefault
                    % if it has a default, just reassign the default
                    x.(f) = prp(j).DefaultValue;

                elseif ~isempty(prp(j).Validation) && ~isempty(prp(j).Validation.Class)
                    % if the type is known, assign based on that type
                    switch prp(j).Validation.Class.Name
                        case "Scan",        x.(f) = ScanCartesian();
                        case "Sequence",    x.(f) = Sequence();
                        case "Transducer",  x.(f) = TransducerArray();
                        case "Waveform",    x.(f) = Waveform.Delta();
                        case "double",      x.(f) = 5;
                        otherwise
                            test.log(3, "No known test value for `"+mc.Name+"."+prp(j).Name+"`.")
                    end
                else  % no known type - manually validate?
                            test.log(3, "No known test value for `"+mc.Name+"."+prp(j).Name+"`.")
                end
            end
            % arrayfun(@(s) test.assertWarningFree(str2func(s)), com);

        end
    end
end