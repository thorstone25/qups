classdef USTest < matlab.unittest.TestCase
    % USTEST - Test the beamforming and simulation methods complete as expected

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    % sim
    properties(TestParameter)
        seq
        tx = struct('Array',TransducerArray('numel',5),'Convex',TransducerConvex('numel',5), 'Matrix', TransducerMatrix('numd',[3 5]), 'Generic', TransducerGeneric("pos", 0.5e-3*randn([3,5])));
        rx = struct('Array',TransducerArray('numel',9),'Convex',TransducerConvex('numel',9), 'Matrix', TransducerMatrix('numd',[3 3]), 'Generic', TransducerGeneric("pos", 0.5e-3*randn([3,9])));
        scan = struct("Cart", ScanCartesian, "Pol", ScanPolar, "Sphr", ScanSpherical, "Gen", ScanGeneric);
        simulator
        apod
    end

    % beamform
    properties(TestParameter)
        beamformer = struct("bfDAS", @bfDAS, "DAS", @DAS, "LUT", @bfDASLUT, "Eik", @bfEikonal);
        apodizer = {"none","angle","fnumb","linem","lines",}
        scan_size = struct("I2D", [15 15 1], "I3D", [15, 15, 5], "I1D", [15 1 1]);
        scan_order = {[1 2 3],[2 1 3],[1 3 2]}; % num2cell(perms([1 2 3]),2)'
        % data_order = {[1 2 3],[2 1 3],[1 3 2]}; % iterate in method
    end

    % generate params
    methods(Static, TestParameterDefinition)
        function simulator = getSimulators()
            simulator.greens = @greens;
            if exist("calc_scat_multi.m",'file')
                simulator.FieldIIm = @calc_scat_multi;
                simulator.FieldIIa = @calc_scat_all;
            end
            if exist("kspaceFirstOrder2D.m", 'file')
                simulator.kWave = @kspaceFirstOrder;
            end
            if exist("simus3.m", "file")
                simulator.MUST = @simus;
            end
        end
        function seq = getSequences()
            typs = ["FSA","PW","FC","DV"]'; % for all types
            sqs = [Sequence, SequenceRadial, SequenceGeneric]; % for all classes
            for i = numel(typs):-1:1
                for j = numel(sqs):-1:1
                    seq{j,i} = setfield(copy(sqs(j)),'type',typs(i));
                end
            end
            nms = typs+"_"+ extractAfter(arrayfun(@(s)string(class(s)), sqs), "Sequence");
            seq = cell2struct(seq(:), nms(:), 1);
        end
        function apod = getScanApodFcns()
            f = string(methods("UltrasoundSystem")); % all methods
            f = f(startsWith(f, "ap"));
            apod = cellstr(f);
        end
    end

    % github settings
    methods(Test, ParameterCombination='pairwise', TestTags = ["Github"])
        function simgeneric_github(testCase, tx, rx, seq, simulator)
            simgeneric(       testCase, tx, rx, seq, simulator);
        end
        function bfusgeneric_github(testCase, tx, rx, seq, scan, beamformer)
                 bfusgeneric(       testCase, tx, rx, seq, scan, beamformer);
        end
        function bfordgeneric_github(testCase, rx, scan, beamformer, scan_size, scan_order)
                 bfordgeneric(       testCase, rx, scan, beamformer, scan_size, scan_order)
        end
        % Not ready ...
        % function apgen_github(testCase, rx, scan, seq, scan_size, scan_order, apod)
        %          apgen(       testCase, rx, scan, seq, scan_size, scan_order, apod)
        % end
    end

    % full settings
    methods(Test, ParameterCombination='exhaustive', TestTags = ["full", "build"])
        function simgeneric_full(testCase, tx, rx, seq, simulator)
                 simgeneric(     testCase, tx, rx, seq, simulator);
        end
        function bfusgeneric_full(testCase, tx, rx, seq, scan, beamformer)
                 bfusgeneric(     testCase, tx, rx, seq, scan, beamformer);
        end
        function bfordgeneric_full(testCase, rx, scan, beamformer, scan_size, scan_order)
                 bfordgeneric(     testCase, rx, scan, beamformer, scan_size, scan_order);
        end
        function apgen_full(testCase, rx, scan, seq, scan_size, scan_order, apod)
                 apgen(     testCase, rx, scan, seq, scan_size, scan_order, apod)
        end
    end

    % Test methods
    methods
        function simgeneric(testCase, tx, rx, seq, simulator)
            % set numpulse for fsa
            if seq.type == "FSA", seq.numPulse = tx.numel; end
            % rx = TransducerArray; % const to reduce combinatorics

            % shift transducers
            [tx.offset, rx.offset] = deal([-2 0 0]*1e-3, [2 0 0]*1e-3);

            % filter incompatible transducer
            for ap = [tx, rx] % do for tx, rx
                if isa(ap, 'TransducerGeneric') && func2str(simulator) == "simus", return; end % incompatible
            end
            
            % filter for identical transducer where required
            if func2str(simulator) == "simus"
                if string(class(tx)) ~= string(class(rx))
                    return; % incompatible
                else
                    tx = rx; % enforce identical
                end
            end

            % construct each tx, rx, seq, scan
            us = UltrasoundSystem('tx', tx, 'rx', rx, 'seq', seq,'recompile',false,'fs',single(4*rx.fc));
            [lb, ub] = bounds([tx.bounds, rx.bounds], 2);
            scat = Scatterers("pos", (lb+ub)./2 + [0 0 50*us.lambda]');
            switch func2str(simulator)
                case {"greens", "calc_scat_multi","calc_scat_all","simus"} % point target
                    chd = simulator(us, scat);
                case {"kwave", "fullwave"} % medium
                    pb = [rx.bounds, tx.bounds]; pb = [min(pb,[],2), max(pb,[],2)];
                    grid = ScanCartesian('xb', pb(1,:),'zb',pb(3,:),'yb',pb(2,:));
                    [grid.dx, grid.dy, grid.dz] = deal(min(us.lambda)/4);
                    med = Medium("c0",scat.c0);
                    med.pertreg{1} = {@(p)argmin(vecnorm(p - scat.pos,2,1)), [scat.c0, 2*med.rho0]}; % point target
                    chd = simulator(us, med, grid);
            end
        end

        function bfusgeneric(testCase, tx, rx, seq, scan, beamformer)
            scan = copy(scan);
            seq = copy(seq);

            % enforce sequence length
            if seq.type == "FSA", seq.numPulse = tx.numel; end

            % shrink image size
            scan.size = [15 15 1];

            % system and apodization
            us = UltrasoundSystem('tx', tx, 'rx', rx, 'seq', seq, 'scan', scan,'recompile',false,'fs',single(4*rx.fc));

            % get dummy data
            T = 32;
            sz = [T, rx.numel, seq.numPulse];
            chd0 = ChannelData('order','TNMF','data',zeros([sz 2],'single'));
            ords = perms(1:3);

            % dummy apodization
            a = ones([scan.size chd0.N, 1], 'single');

            % simulation grid
            pb = [rx.bounds, tx.bounds]; pb = [min(pb,[],2), max(pb,[],2)];
            pb = pb + [-1 1] .* 3*max(us.lambda);
            grid = ScanCartesian('x', pb(1,:),'z',pb(3,:),'y',pb(2,:));
            [grid.dx, grid.dy, grid.dz] = deal(min(us.lambda)/4);

            % beamform
            bscan = us.scan; % default
            for i = 1:size(ords,1)
                chd = permuteD(chd0, [ords(i,:), 4]);
                switch func2str(beamformer)
                    case "bfMigration" % "lin-pw" only
                        if us.seq.type ~= "PW" || ~any(arrayfun(@(s)isa(us.seq.xdc,s), "Transducer"+["Array","Matrix"]))
                            return; % compatibility
                        end
                        [b, bscan] = beamformer(us,chd,'apod',a);
                    case "bfEikonal" % fsa only, grid
                        if us.seq.type ~= "FSA", return; end % compatibility
                        if nnz(grid.size > 1) == 3, return; end % segfaults in 3D
                        b = beamformer(us,chd,Medium(),grid,'apod',a);
                    case "bfDASLUT"
                        [~, tn, tm] = bfDAS(us, chd, 'delay_only',true);
                        b = beamformer(us,chd,tn,tm,'apod',a);
                    otherwise
                        b = beamformer(us,chd,'apod',a);
                end
            end

            % imagesc(us.scan, b);

            % TODO: add verification
            testCase.assertEqual(size(b,1:3), double(bscan.size));
            % testCase.verifyFail("Unimplemented test");
        end
    
        function bfordgeneric(testCase, rx, scan, beamformer, scan_size, scan_order)
            scan = copy(scan);
            seq = SequenceRadial('type', 'PW', 'angles', [-10 0 10]);
            if isa(rx, "TransducerConvex"), seq.apex = rx.center; end

            % shrink image size
            scan.size = scan_size;
            scan.order = scan.order(scan_order);

            % system and apodization
            tx = TransducerArray('numel', 9);
            us = UltrasoundSystem('tx', tx, 'rx', rx, 'seq', seq, 'scan', scan,'recompile',false,'fs',single(4*rx.fc));

            % get dummy data
            T = 32;
            sz = [T, rx.numel, seq.numPulse];
            ord = 'TNM';
            chd = ChannelData('order',[ord,'F'],'data',zeros([sz 2],'single'));
            assert(chd.T == T);
            assert(chd.N == us.rx.numel);
            assert(chd.M == us.seq.numPulse);

            % simulation grid
            pb = [rx.bounds, tx.bounds]; pb = [min(pb,[],2), max(pb,[],2)];
            pb = pb + [-1 1] .* 3*max(us.lambda);
            grid = ScanCartesian('x', pb(1,:),'z',pb(3,:),'y',pb(2,:));
            [grid.dx, grid.dy, grid.dz] = deal(min(us.lambda)/4);

            % get unique list of apodization sizes
            bitctr = @(b,n)logical(mod(idivide(-1+int64(b^0:b^n)', int64(b.^(0:n-1)), "floor"),b)); % options base b choices n
            scl = bitctr(2, 5); % true/false, 5 dims
            asz = repmat([us.scan.size, us.rx.numel, us.seq.numPulse], [size(scl,1), 1]);
            for i = 1:size(asz,1), asz(i,scl(i,:)) = 1; end % scalarize
            asz = unique(asz, 'rows');

            % beamform
            bscan = us.scan; % default
            switch func2str(beamformer)
                case "bfMigration" % "lin-pw" only
                    if us.seq.type ~= "PW" || ~any(arrayfun(@(s)isa(us.seq.xdc,s), "Transducer"+["Array","Matrix"]))
                        return; % compatibility
                    end
                    [b, bscan] = beamformer(us,chd); % no apodization
                case "bfEikonal" % fsa only, grid
                    if us.seq.type ~= "FSA", return; end % compatibility
                    if nnz(grid.size > 1) == 3, return; end % segfaults in 3D
                    for j = 1:size(asz,1)
                        a = randn(asz(j,:),'single');
                        b = beamformer(us,chd,Medium(),grid,'apod',a);
                    end
                case "bfDASLUT"
                    [~, tn, tm] = bfDAS(us, chd, 'delay_only',true);
                    for j = 1:size(asz,1)
                        a = randn(asz(j,:),'single');
                        b = beamformer(us,chd,tn,tm,'apod',a);
                    end
                case "bfAdjoint"
                    s = asz == 1; % apd scalar
                    i = s(:,4) | s(:,5) | all(s(:,1:3),2); % valid apodization
                    asz(~i,:) = []; % delete invalid
                    for j = 1:size(asz,1)
                        a = randn(asz(j,:),'single');
                        b = beamformer(us,chd,'apod',a);
                    end
                otherwise
                    for j = 1:size(asz,1)
                        a = randn(asz(j,:),'single');
                        b = beamformer(us,chd,'apod',a);
                    end
            end

        end

        function apgen(testCase, rx, scan, seq, scan_size, scan_order, apod)
            
            scan = copy(scan);
            seq  = copy(seq);

            % shrink image size
            scan.size = scan_size;
            scan.order = scan.order(scan_order);

            % make system
            if seq.type == "FSA", seq.numPulse = rx.numel; end
            us = UltrasoundSystem('xdc', rx, 'seq', seq, 'scan', scan, 'fs',single(4*rx.fc));

            % try to make apodization
            f = str2func(apod);
            try ap = f(us); % success!
            catch ME % failure - QUPS should tell the user
                warning(ME.identifier, '%s', ME.message);
                testCase.assertTrue(ME.identifier == "QUPS:UltrasoundSystem:UnsupportedScan", "Caught " + ME.identifier + ":" + newline + ME.message);
                return;
            end

            testCase.assertTrue(all(any(size(ap,1:5) == [[us.scan.size, us.xdc.numel, us.seq.numPulse]; ones(1,5)],1),2));

            % TODO: apply


        end
    end

end