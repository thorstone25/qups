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
        end
        function ustb_scan(tst, scns)
            tst.assumeTrue(logical(exist('uff','class'))) % need USTB for this
            uscn = QUPS2USTB(scns);
            qscn = Scan.UFF(uscn);
        end
        function ustb_seq(tst, seqs, xdcs)
            tst.assumeTrue(logical(exist('uff','class'))) % need USTB for this
            if seqs.type == "FSA", seqs.numPulse = xdcs.numel; end
            useq = QUPS2USTB(seqs, xdcs);
            qseq = Sequence.UFF(useq);
        end        
        function ustb_sys(tst, xdcs, seqs, scns)
            tst.assumeTrue(logical(exist('uff','class'))) % need USTB for this
            us = UltrasoundSystem('scan', scns, 'seq', seqs, 'xdc', xdcs);
            if seqs.type == "FSA", seqs.numPulse = xdcs.numel; end
            [T,N,M,F] = deal(16, xdcs.numel, seqs.numPulse, 2);

            chd = ChannelData('data',zeros([T,N,M,F],'single'), 'order','TNMF');
            chds = [chd, permuteD(chd,[1,3,2,4])]; % swap tx/rx
            for chd = chds 
                [uchd, uscn] = QUPS2USTB(us, chd, 0);
                [uso, chdo] = UltrasoundSystem.UFF(uchd, uscn);
            end
        end
    end
end