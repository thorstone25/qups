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
    end
    % data for iteration
    properties
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
            if exist("fullwave2_executable", "file") ...
            && exist("mapToCoords", "file")
                simulator.fullwave = @fullwaveSim;
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
    methods(Test, ParameterCombination="pairwise", TestTags = ["Github", "build"])
        function simgeneric_github(testCase, tx, rx, seq, simulator)
            simgeneric(       testCase, tx, rx, seq, simulator);
        end
        function bfusgeneric_github(testCase, tx, rx, seq, scan, beamformer)
                 bfusgeneric(       testCase, tx, rx, seq, scan, beamformer);
        end
        function bfordgeneric_github(testCase, rx, scan, beamformer)
                 bfordgeneric(       testCase, rx, scan, beamformer)
        end
        function apgen_github(testCase, rx, scan, seq, apod)
                 apgen(       testCase, rx, scan, seq, apod)
        end
    end

    % full settings
    methods(Test, ParameterCombination="exhaustive", TestTags = ["full"])
        function simgeneric_full(testCase, tx, rx, seq, simulator)
                 simgeneric(     testCase, tx, rx, seq, simulator);
        end
        function bfusgeneric_full(testCase, tx, rx, seq, scan, beamformer)
                 bfusgeneric(     testCase, tx, rx, seq, scan, beamformer);
        end
        function bfordgeneric_full(testCase, rx, scan, beamformer)
                 bfordgeneric(     testCase, rx, scan, beamformer);
        end
        function apgen_full(testCase, rx, scan, seq, apod)
                 apgen(     testCase, rx, scan, seq, apod)
        end
    end

    % Test methods
    methods
        function simgeneric(testCase, tx, rx, seq, simulator)
            import matlab.unittest.fixtures.TemporaryFolderFixture;
            import matlab.unittest.fixtures.CurrentFolderFixture;
            simname = func2str(simulator); % alias

            [tx, rx, seq] = dealfun(@copy, tx, rx, seq); % copy semantics

            % set numpulse for fsa
            if seq.type == "FSA", seq.numPulse = tx.numel; end
            % rx = TransducerArray; % const to reduce combinatorics

            % shift transducers
            [tx.offset, rx.offset] = deal([-1 0 0]*1e-3, [1 0 0]*1e-3);

            % filter for identical transducer where required
            if simname == "simus"
                if string(class(tx)) ~= string(class(rx))
                    return; % incompatible
                else
                    tx = rx; % enforce identical
                end
               if isa(tx, 'TransducerGeneric'), return; end % incompatible
            end

            % construct each tx, rx, seq, scan
            us = UltrasoundSystem('tx', tx, 'rx', rx, 'seq', seq,'recompile',false,'fs',single(4*rx.fc));
            [lb, ub] = bounds([tx.bounds, rx.bounds], 2);
            scat = Scatterers("pos", (lb+ub)./2 + [0 0 50*us.lambda]');
            switch simname
                case {"greens", "calc_scat_multi","calc_scat_all","simus"} % point target
                    chd = simulator(us, scat);
                    if simname == "greens"
                        chd = simulator(us, scat, 'device', 0, 'tall', 0);
                        % chd = simulator(us, scat, 'device', 0, 'tall',1); % strange behaviour as a test?
                    end
                case {"kspaceFirstOrder", "fullwaveSim"} % medium
                    pb = [rx.bounds, tx.bounds]; pb = [min(pb,[],2), max(pb,[],2)]; % transducer boundaries
                    pb = pb + [-1 1] .* us.lambda; % tolerance in x/y/z
                    pb(3,2) = max(pb(3,2),2*us.lambda); % add depth
                    if ~(isa(us.tx, "TransducerMatrix") || isa(us.rx, "TransducerMatrix")), pb(2,:) = 0; end % unless actually 3D ...
                    zero2nan = @(x) x.*(x./x); % HACK
                    dp = min(arrayfun(@(ap) min(zero2nan(vecnorm(ap.positions - swapdim(ap.positions,2,3),2,1)),[],'all','omitnan'), [us.tx, us.rx])); % min element spacing
                    grid = ScanCartesian('x', pb(1,:),'z',pb(3,:),'y',pb(2,:));
                    [grid.dx, grid.dy, grid.dz] = deal(min([us.lambda, dp])/4);
                    us.scan = grid; % same
                    med = Medium("c0",scat.c0);
                    med.pertreg{1} = {@(p)argmin(vecnorm(p - scat.pos,2,1)), [scat.c0, 2*med.rho0]}; % point target
                    switch simname
                        case "kspaceFirstOrder" 
                            % should always work
                            us.fs = 2*us.fc; % unstable, but faster
                            chd = simulator(us, med, "CFL_max", 10);

                            % further test element mapping
                            for elm = ["nearest", "linear", "karray-direct", "karray-depend"]
                                kspaceFirstOrder(us, med, "CFL_max", 10, "ElemMapMethod", elm, "parenv", 0);
                            end
                            
                        case "fullwaveSim"
                            % only test convex/array, single xdc, 2D
                            typ = "Transducer"+["Convex","Array",];
                            % testCase.assumeTrue(us.tx == us.rx, "The transmit and receive must be identical for " + simname);
                            testCase.assumeTrue(nnz(grid.size > 1) == 2, "Only 2D sims are supported by " + simname + "(requested "+nnz(grid.size > 1)+"D).");
                            testCase.assumeTrue( ...
                                any(arrayfun(@(T)isa(us.tx,T), typ)), ...
                                "A "+class(us.tx)+" transmitter is not supported for " + simname ...
                                );
                            testCase.assumeTrue( ...
                                any(arrayfun(@(T)isa(us.rx,T), typ)), ...
                                "A "+class(us.rx)+" receiver is not supported for " + simname ...
                                );
                            
                            chd = simulator(us, med, "CFL_max", 10); % supported
                    end
            end
        end

        function bfusgeneric(test, tx, rx, seq, scan, beamformer)
            scan = copy(scan);
            seq = copy(seq);

            % enforce sequence length
            if seq.type == "FSA"
                seq.numPulse = tx.numel; 
            elseif seq.type == "PW"
                th = [-25 0 25];
                seq.focus = [sind(th); 0*th; cosd(th)];
            elseif seq.type ~= "DV"
                th = [-25 0 25];
                seq.focus = 50e-3 * [sind(th); 0*th; cosd(th)];
            else
                th = [-25 0 25];
                seq.focus = -5e-3 * [sind(th); 0*th; cosd(th)];                
            end

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
                            return; % incompatible - should draw an error
                        else
                            [b, bscan] = beamformer(us,chd); % works
                            b = test.assertWarning(@() beamformer(us,chd), ...
                                "QUPS:bfMigration:artefacts"); % draw a warning for no second output
                        end
                    case "bfEikonal" % fsa only, grid
                        if us.seq.type ~= "FSA", return; end % compatibility
                        if nnz(grid.size > 1) == 3, return; end % segfaults in 3D
                        b = beamformer(us,chd,Medium(),grid,a);
                        [~, tn, tm] = beamformer(us); % should work too
                        test.assertWarning(@()beamformer(us), ... why would you even ...
                            "QUPS:bfEikonal:"+["unexpectedOutput","emptyChannelData"]);

                    case "bfDASLUT"
                        [~, tn, tm] = bfDAS(us, 'delay_only', true); % should draw empty chd default
                        b = beamformer(us, chd, tn, tm, a); % works
                        if us.seq.type == "FSA" && chd.N == chd.mag2db
                            b = beamformer(us,chd,tn ,a); % shoudl work too
                        end

                        % fail if any dim is sub-indexed wrong
                        dn = find(size(tn) > 2);
                        dm = find(size(tm) > 2);
                        for d = dn(:)'
                            test.assertError(@()beamformer(us,chd,sub(tn,1:2,d),tm,a), ...
                                "QUPS:UltrasoundSystem:bfDASLUT:incompatibleReceiveDelayTable" ...
                                ); % bad sizing on rx
                        end
                        for d = dm(:)'
                            test.assertError(@()beamformer(us,chd,tn, sub(tm,1:2,d),a), ...
                                "QUPS:UltrasoundSystem:bfDASLUT:incompatibleTransmitDelayTable" ...
                                ); % bad sizing on tx
                        end
                        % fail if chd changes size
                        test.assertError(@()beamformer(us,[chd, subD(chd,1:2,chd.ndim)],tn,tm,a), ...
                            "QUPS:UltrasoundSystem:bfDASLUT:nonUniqueReceiverSize" ...
                            ); % bad sizing on rx
                        test.assertError(@()beamformer(us,[chd, subD(chd,1:2,chd.mdim)],tn,tm,a), ...
                            "QUPS:UltrasoundSystem:bfDASLUT:nonUniqueTransmitSize" ...
                            ); % bad sizing on tx
                    case "DAS"
                        if canUseGPU(), d = [0 gpuDevice().Index]; else, d = 0; end
                        for d = d
                            b = beamformer(us,chd,a, 'device',d);
                            b = beamformer(us,chd,a, 'device',d,'keep_tx',1);
                            b = beamformer(us,chd,a, 'device',d,'keep_rx',1);
                            b = beamformer(us,chd,a, 'device',d,'keep_tx',1,'keep_rx',1);
                        end
                    otherwise
                        b = beamformer(us,chd,a);
                end
            end

            % imagesc(us.scan, b);

            % TODO: add verification
            test.assertEqual(size(b,1:3), double(bscan.size));
            % testCase.verifyFail("Unimplemented test");
        end
    
        function bfordgeneric(testCase, rx, scan, beamformer)
            [scan0, scan] = dealfun(@copy, scan);
            seq = SequenceRadial('type', 'PW', 'angles', [-10 0 10]);
            if isa(rx, "TransducerConvex"), seq.apex = rx.center; end

            % Iterate through sizes
            for f = string(fieldnames(testCase.scan_size))'
            for j = 1:numel(testCase.scan_order)
            
            % shrink image size
            scan.size = testCase.scan_size.(f);
            scan.order = scan.order(testCase.scan_order{j});

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
                        b = beamformer(us,chd,Medium(),grid,a);
                    end
                case "bfDASLUT"
                    [~, tn, tm] = bfDAS(us, chd, 'delay_only',true);
                    for j = 1:size(asz,1)
                        a = randn(asz(j,:),'single');
                        b = beamformer(us,chd,tn,tm,a);
                    end
                case "bfAdjoint"
                    s = asz == 1; % apd scalar
                    i = s(:,4) | s(:,5) | all(s(:,1:3),2); % valid apodization
                    asz(~i,:) = []; % delete invalid
                    for j = 1:size(asz,1)
                        a = randn(asz(j,:),'single');
                        b = beamformer(us,chd,a);
                    end
                otherwise
                    for j = 1:size(asz,1)
                        a = randn(asz(j,:),'single');
                        b = beamformer(us,chd,a);
                    end
            end
            end
            end
        end

        function apgen(testCase, rx, scan, seq, apod)
            
            scan = copy(scan);
            seq  = copy(seq);

            % Iterate through sizes
            for i = string(fieldnames(testCase.scan_size))'
            for j = 1:numel(testCase.scan_order)
                
            % shrink image size
            scan.size = testCase.scan_size.(i);
            scan.order = scan.order(testCase.scan_order{j});

            % make system
            if seq.type == "FSA", seq.numPulse = rx.numel; end
            us = UltrasoundSystem('xdc', rx, 'seq', seq, 'scan', scan, 'fs',single(4*rx.fc));

            % try to make apodization
            f = str2func(apod);
            try ap = f(us); % success?!
            catch ME % failure - QUPS should tell the user
                warning(ME.identifier, '%s', ME.message);
                testCase.assertTrue(any(ME.identifier == [ ... any of the following errors are fine
                    "QUPS:UltrasoundSystem:UnsupportedScan" ...
                    "MATLAB:noSuchMethodOrField" ... we warn before hand if we aren't sure it's there 
                    ... "MATLAB:getReshapeDims:notSameNumel" ... probably a bug ...
                    ]), "Caught " + ME.identifier + ":" + newline + ME.message ...
                    );
                return;
            end

            dsz = [us.scan.size, us.xdc.numel, us.seq.numPulse]; % data size
            testCase.assertTrue(all(size(ap,1:5) == dsz | size(ap,1:5) == 1), ...
                "The output apodization of size ("   + join(string(size(ap,1:5)), " x ") +") "+... 
                " is not compatible with the data (" +join(string(dsz), " x ") +").");

            % TODO: apply
            end
            end

        end
    end

end