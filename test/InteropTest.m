classdef(TestTags = ["full","Github","build"]) InteropTest < matlab.unittest.TestCase
    properties(TestParameter)
        xdcs = struct( ...
            "Linr", TransducerArray("numel", 3), ...
            "Conv", TransducerConvex("numel", 3), ...
            "Mtrx", TransducerMatrix("numd", [3 3]), ...
            "Gnrc", TransducerGeneric("az", [-25 25],"el",[0 0],"pos", 1e-3*[-5 0 0; 5 0 0]') ...
            );
        seqs = struct( ...
            "FSA", Sequence(), ...
            "PW", SequenceRadial("angles",-10:10:10), ...
            "FC", Sequence('type','FC','focus',1e-3*([0 0  50]' + [-10 0 10] * [1 0 0]')), ...
            "DV", Sequence('type','DV','focus',1e-3*([0 0 -10]' + [-10 0 10] * [1 0 0]'))  ...
            );
        scns = struct( ...
            "Cart", scale(ScanCartesian('x',-5:5:5,'z',0:5:10           ),'dist', 1e-3), ...
            "Polr", scale(ScanPolar(    'a',-5:5:5,'r',0:5:10           ),'dist', 1e-3), ...
            "Sphr", scale(ScanSpherical('e',-5:5:5,'r',0:5:10,'a',-5:5:5),'dist', 1e-3) ...
        );
    end

    methods(TestClassSetup)
    end
    methods(Test)
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
            tst.assertEqual(scns.positions(), qscn.positions());
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
                [uso, chdo] = UltrasoundSystem.UFF(uchd, uscn);
                chdo.order(4) = 'F';
                [~, ord] = ismember(chd.order, chdo.order);
                chdo = permuteD(chdo, ord);
                tst.assertThat(chdo.t0, IsEqualTo(chd.t0, "Within", AbsoluteTolerance(1e-3)));
                tst.assertEqual(size(chdo.data), size(chd.data));
            end
        end
        function ustb_ext(tst) % test objects in USTB that don't have a direct QUPS analogy
            % test scans
            uscn = cell(1,2);

            % 3D plane scan
            sca = uff.linear_3D_scan();
            sca.radial_axis=linspace(-20e-3,20e-3,256)';
            sca.axial_axis=linspace(0e-3,40e-3,256)';
            sca.roll=0;
            uscn{1} = sca;

            % rotated scan?
            sca = uff.linear_scan_rotated();
            sca.x_axis=linspace(-20e-3,20e-3,256)';
            sca.z_axis=linspace(0e-3,40e-3,256)';
            sca.rotation_angle = 20*pi/180;
            sca.center_of_rotation = [0 0 0]';
            uscn{2} = sca;
            
            for i = 1:numel(uscn)
                sca = uscn{i};
                scn = Scan.UFF(sca);
                tst.assertEqual(scn.positions()', sca.xyz, "Conversion from " + class(sca) + " to " + class(scn) + " failed.");
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
    end
end