classdef(TestTags = ["full","Github","build","syntax"]) InteropTest < matlab.unittest.TestCase
    % INTEROPTEST - Test interoperability for extension packages (FieldII, USTB, Verasonics)

    %#ok<*ASGLU,*NASGU> ignore unused outputs
    properties(TestParameter)
        xdcs = struct( ...
            "Linr", TransducerArray("numel", 3), ...
            "Conv", TransducerConvex("numel", 3), ...
            "Mtrx", TransducerMatrix("numd", [3 3]), ...
            "Gnrc", TransducerGeneric("az", [-25 25],"el",[0 0],"pos", 1e-3*[-5 0 0; 5 0 0]') ...
            );
        seqs = struct( ...
            "FSA", SequenceGeneric(), ...
            "PW", SequenceRadial("angles",-10:10:10), ...
            "FC", Sequence('type','FC','focus',1e-3*([0 0  50]' + [-10 0 10] * [1 0 0]')), ...
            "DV", Sequence('type','DV','focus',1e-3*([0 0 -10]' + [-10 0 10] * [1 0 0]')),  ...
            "VS", Sequence('type','VS','focus',1e-3*([0 0  50]' + [-10 0 10] * [1 0 0]')) ...
            );
        scns = struct( ...
            "Cart", scale(ScanCartesian('x',-5:5:5,'z',0:5:10           ),'dist', 1e-3), ...
            "Polr", scale(ScanPolar(    'a',-5:5:5,'r',0:5:10           ),'dist', 1e-3), ...
            "Sphr", scale(ScanSpherical('e',-5:5:5,'r',0:5:10,'a',-5:5:5),'dist', 1e-3) ...
        );
    end
    properties(TestParameter)
        vsx_files % VSX - Verasonics data structs
    end

    methods(Static, TestParameterDefinition)
        function vsx_files = loadVSXData()
            % load files in data/VSX-*.mat into the class
            prj = currentProject();
            vfls = dir(fullfile(prj.RootFolder, "data","VSX-test","**","*.mat")); % VSX test files
            val = false(size(vfls));
            % vs = {};
            for i = 1 : numel(vfls)
                fl = vfls(i);

                % lazy load
                dat = matfile(fullfile(fl.folder, fl.name), 'Writable',false);

                % verify
                flds = string(fieldnames(dat));
                val(i) = all(ismember(["Trans", "TW", "TX", "Receive", "PData", "Resource"], flds)); % req'd for testing
                val(i) = val(i) && any(ismember(flds, ["RcvData", "RData"])); % need one or the other

                % check for specific post-processed properties
                val(i) = val(i) && any(isfield(dat.TW, "TriLvlWvfm"+["","_Sim"]));
                val(i) = val(i) &&     isfield(dat.Receive, "demodFrequency");

                % import
                if val(i)
                    % move RData to RcvData
                    if isfield(dat, "RData")
                        if ~isfield(dat, "RcvData")
                            dat.RcvData = {dat.RData};
                        else
                            warning( ...
                                "QUPS:ExampleTest:duplicateData", ...
                                "Found 'RData' and 'RcvData' in file " + fl.name +...
                                " - RData will be ignored." ...
                                );
                        end
                    end

                    % only save specified properties
                    % dat = rmfield(dat, setdiff(flds, ["Trans", "TW", "TX", "Receive", "PData", "Resource", "RcvData"]));

                    % save
                    % vs{i} = dat;
                end
            end
            vsx_files = cellstr(string(fullfile({vfls(val).folder},{vfls(val).name})));
            if isempty(vsx_files), vsx_files = {''}; end
        end
    end
    methods(TestClassSetup)
        function silenceAcceptableWarnings(tst)
            lids = [...
                "QUPS:Sequence:DeprecatedValue"    ... in Sequence.UFF
                "QUPS:QUPS2USTB:ambiguousSequence" ... in Sequence.QUPS2USTB
                ];
            W = warning(); % get current state
            arrayfun(@(l) warning('off', l), lids); % silence
            tst.addTeardown(@() warning(W)); % restore on exit
            if ~isempty(gcp('nocreate')) % any pool - execute on each worker
                ws = parfevalOnAll(@warning, 1); wait(ws);% current state
                ws = fetchOutputs(ws); % retrieve
                [~, i] = unique(string({ws.identifier}), 'stable');
                ws = ws(i); % select first unique set (assumer identical)
                wait(arrayfun(@(l) parfevalOnAll(@warning, 0, 'off', l), lids)); % silence
                tst.addTeardown(@() parfevalOnAll(@warning, 0, ws)); % restore on exit
            end
        end
    end
    methods(Test)
        function fieldII_xdc(tst, xdcs)
            tst.assumeTrue(logical(exist('field_init','file'))) % need FieldII for this
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;

            xdcs.patches(); % assumes defaults
            pchq = xdcs.patches([4 2]);
            pchf = xdcs.getFieldIIPatches([4 2]);
            [pchq, pchf] = dealfun(@(x) cell2mat(swapdim(cellfun(@(x) {cat(3, x{:})}, x),1:2,4:5)), pchq, pchf);

            pnf = xdcs.getFieldIIPositions();
            tst.assertThat(xdcs.positions(), IsEqualTo(pnf , "Within", AbsoluteTolerance(1e-6)));
            % tst.assertThat(pchq            , IsEqualTo(pchf, "Within", AbsoluteTolerance(1e-6))); % extra transposition

        end
        function MUST_xdc(tst, xdcs)
            if any(arrayfun(@(t)isa(xdcs,t), "Transducer" + ["Array","Convex","Matrix"])) % valid
                p = xdcs.getSIMUSParam(); % valid - should run
            else
                tst.assertError(@()xdcs.getSIMUSParam(), ...
                    "QUPS:Transducer:unsupportedTransducer") % should error
            end
        end
        function ustb_xdc(tst, xdcs)
            tst.assumeTrue(logical(exist('uff','class'))) % need USTB for this
            uxdc = QUPS2USTB(xdcs);
            qxdc = Transducer.UFF(uxdc);
            tst.assertEqual(xdcs.positions(), qxdc.positions());
        end
        function ustb_scan(tst, scns)
            tst.assumeTrue(logical(exist('uff','class'))) % need USTB for this
            uscn = QUPS2USTB(scns);
            qscn = Scan.UFF(uscn);
            tst.assertEqual(reshape(scns.positions(),3,[]), reshape(qscn.positions(),3,[]));
        end
        function ustb_seq(tst, seqs, xdcs)
            tst.assumeTrue(logical(exist('uff','class'))) % need USTB for this
            if seqs.type == "FSA", seqs.numPulse = xdcs.numel; end
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;
            useq = QUPS2USTB(seqs, xdcs);
            qseq = Sequence.UFF(useq);
            tst.assertThat(seqs.focus   , IsEqualTo(qseq.focus   , "Within", AbsoluteTolerance(1e-9)));
            tst.assertThat(seqs.numPulse, IsEqualTo(qseq.numPulse, "Within", AbsoluteTolerance(1e-9)));
        end        
        function ustb_sys(tst, xdcs, seqs, scns)
            tst.assumeTrue(logical(exist('uff','class'))) % need USTB for this
            us = UltrasoundSystem('scan', scns, 'seq', seqs, 'xdc', xdcs);
            if seqs.type == "FSA", seqs.numPulse = xdcs.numel; end
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;
            [T,N,M,F] = deal(16, xdcs.numel, seqs.numPulse, 2);

            chd = ChannelData('data',zeros([T,N,M,F],'single'), 'order','TNMF', 't0',randn([1,1,M,1]));
            chds = [chd, permuteD(chd,[1,3,2,4])]; % swap tx/rx
            for chd = chds 
                [uchd, uscn] = QUPS2USTB(us, chd, 0);
                [~, chd1] = UltrasoundSystem.UFF(uchd, uscn);
                chd2 = ChannelData.UFF(uchd);
                for chdo = [chd1, chd2]
                chdo.order(4) = 'F';
                [~, ord] = ismember(chd.order, chdo.order);
                chd3 = permuteD(chdo, ord);
                tst.assertThat(chd3.t0, IsEqualTo(chd.t0, "Within", AbsoluteTolerance(1e-3)));
                tst.assertEqual(size(chd3.data), size(chd.data));
                end
            end
        end
        function ustb_ext(tst) % test objects in USTB that don't have a direct QUPS analogy
            tst.assumeTrue(logical(exist('uff','class'))) % need USTB for this

            % test scans
            uscn = cell(1,1);

            % 3D plane scan (deprecated)
            %{
            sca = uff.linear_3D_scan();
            sca.radial_axis=linspace(-20e-3,20e-3,256)';
            sca.axial_axis=linspace(0e-3,40e-3,256)';
            sca.roll=0;
            uscn{1} = sca;
            %}

            sca = uff.linear_scan();
            sca.x_axis=linspace(-20e-3, 20e-3, 256)';
            sca.z_axis=linspace(  0e-3, 40e-3, 256)';
            uscn{1} = sca;

            % rotated scan? (deprecated)
            %{
            sca = uff.linear_scan_rotated();
            sca.x_axis=linspace(-20e-3,20e-3,256)';
            sca.z_axis=linspace(0e-3,40e-3,256)';
            sca.rotation_angle = 20*pi/180;
            sca.center_of_rotation = [0 0 0]';
            uscn{2} = sca;
            %}

            
            for i = 1:numel(uscn)
                sca = uscn{i};
                scn = Scan.UFF(sca);
                tst.assertEqual(reshape(scn.positions(),3,[])', sca.xyz, "Conversion from " + class(sca) + " to " + class(scn) + " failed.");
                tst.assertEqual(sca.xyz, scn.QUPS2USTB().xyz, "Conversion from " + class(scn) + " to " + class(sca) + " failed.");
            end

            % curvilinear matrix array
            uxdc = uff.curvilinear_matrix_array();
            uxdc.radius_x = 40e-3;
            uxdc.element_width  = 0.3e-3;
            uxdc.element_height = 0.3e-3;
            uxdc.N_x = 5;
            uxdc.N_y = 5;
            uxdc.pitch_x = 0.5e-3;
            uxdc.pitch_y = 0.5e-3;

            xdc = Transducer.UFF(uxdc);
            tst.assertEqual(xdc.positions()', uxdc.geometry(:,1:3), "Conversion from " + class(xdc) + " to " + class(uxdc) + " failed.");
            tst.assertEqual(uxdc.geometry(:,1:3), xdc.QUPS2USTB().geometry(:,1:3), "Conversion from " + class(uxdc) + " to " + class(xdc) + " failed.");

        end
        function vsx_ext(~, vsx_files)
            % check if valid
            if isempty(vsx_files), return; end % nothing to test

            % select and extract data
            v = load(vsx_files);
            c0 = v.Resource.Parameters.speedOfSound;

            % enforce RcvData property from RData
            if ~isfield(v, "RcvData"), v.RcvData = {v.RData}; end

            % TODO: identify multiplexing

            % US loader
            opts = reshape({'c0', c0, 'PData', v.PData, 'Receive', v.Receive, 'RcvData', v.RcvData},2,[]);
            N = size(opts, 2); % number of options
            for i = uint64(0 : 2 ^ N - 1) % for each option permutation
                j = logical(bitget(i, 1:N)); % use bit mask to select options permutation
                [us, chd] = UltrasoundSystem.Verasonics(v.Trans, v.TX,       opts{:,j}); % should run
                [us, chd] = UltrasoundSystem.Verasonics(v.Trans, v.TX, v.TW, opts{:,j}); % should run
            end

            % Chd loader
            opts = reshape({'frames', 2, 'insert0s', false},2,[]);
            N = size(opts, 2); % number of options
            for i = uint64(0 : 2 ^ N - 1) % for each option permutation
                j = logical(bitget(i, 1:N)); % use bit mask to select options permutation
                [chd, fmod, smode] = ChannelData.Verasonics(v.RcvData, v.Receive,          opts{:,j});
                [chd, fmod, smode] = ChannelData.Verasonics(v.RcvData, v.Receive, v.Trans, opts{:,j});
            end

            % Sequence loader
            opts = reshape({'c0', c0, 'xdc', us.xdc},2,[]);
            N = size(opts, 2); % number of options
            for i = uint64(0 : 2 ^ N - 1) % for each option permutation
                j = logical(bitget(i, 1:N)); % use bit mask to select options permutation
                [seq, t0] = Sequence.Verasonics(v.TX, v.Trans,       opts{:,j});
                [seq, t0] = Sequence.Verasonics(v.TX, v.Trans, v.TW, opts{:,j});
            end

            % Transudcer
            xdc = Transducer.Verasonics(v.Trans);
            xdc = Transducer.Verasonics(v.Trans, c0);

            % Scan
            scan = Scan.Verasonics(v.PData, c0 / xdc.fc);

            % Waveform
            [wvtri, wvm1wy, wvm2wy] = Waveform.Verasonics(v.TW);

            % Scatterers
            if isfield(v, 'Media')
                sct = Scatterers.Verasonics(v.Media, "c0",c0, "scale", c0 / xdc.fc);
                Scatterers.Verasonics(v.Media);
                if ~isfield(v.Media, 'attenuation')
                    v.Media.attenuation = 0;
                    Scatterers.Verasonics(v.Media);
                    v.Media.attenuation = 0.5;
                    Scatterers.Verasonics(v.Media);
                end
            end

        end
    end
end