classdef(TestTags = ["full", "Github", "build", "syntax"]) ChdTest < matlab.unittest.TestCase
    % CHDTEST - This class test that ChannelData methods function properly 
    
    %#ok<*NASGU> The outputs will not be televised

    methods(Test)
        function typeCheck(tst)
            chd = ChannelData('data', complex(randn([16 8 4 2], 'single'),randn([16 8 4 2], 'single')));

            % all precision types
            chdo = doubleT(chd);
            chdo = singleT(chd);
            chdo = uint64T(chd);
            chdo = uint32T(chd);
            chdo = uint16T(chd);
            chdo =  uint8T(chd);
            chdo =  int64T(chd);
            chdo =  int32T(chd);
            chdo =  int16T(chd);
            chdo =   int8T(chd);

            % data types
            z = double(chd);
            z = single(chd);
            z = uint64(chd);
            z = uint32(chd);
            z = uint16(chd);
            z =  uint8(chd);
            z =  int64(chd);
            z =  int32(chd);
            z =  int16(chd);
            z =   int8(chd);

            % tall, sparse
            chdo = tall(chd);
            sparse(ChannelData('data', randi([0 1], [16 32])))

            % unit conversion check
            chdo = angle(complex(chd));
            chdo = deg2rad(rad2deg(chd));
            chdo = mag2db(abs(chd));
            chdo = mod2db(    chd );
            chdo = real(imag(complex(chd)));

            % apply fun to data same as apply fun directly
            tst.assertEqual(conj(chd.data)  , conj(chd).data);
            tst.assertEqual(mod2db(chd.data), mod2db(chd).data);

            % qualifiers
            isreal(chd); 
            istall(chd);
            classUnderlying(chd);
            underlyingType(chd);


        end

        function freqDomainCheck(tst)
            [T, N, M, F] = deal(128, 8, 4, 2);
            chd = ChannelData('data', randn([T M N F], 'single'), 't0', randn([1 1 N F]), 'fs', 1, 'order', 'TMNF');

            % fft
             fftshift( fft(chd));
            ifftshift(ifft(chd));
            
            % Filtering
            D1 = chd.getLowpassFilter(0.2, 5);
            D2 = chd.getPassbandFilter([0.1 0.4], 5);
            D3 = chd.getLowpassFilter(0.2);
            D4 = chd.getPassbandFilter([0.1 0.4]);
            for D = [D1 D2]
                filter(  chd, D)
                filtfilt(chd, D)
                fftfilt( chd, D)
            end

            % sampling
            downsample(downmix(resample(chd,5), 1/4), 2);
            

        end

        function arithmetic(tst)
            [T,M,N,F] = deal(16,8,4,2);
            a = ChannelData('data', rand([T M N F]), 't0', randn([1 M 1 F]), 'order', 'TMNF');
            b = ChannelData('data', rand([T M N F]), 't0', a.t0            , 'order', 'TMNF');
            c = copy(a); c.t0 = c.t0 + rand(size(c.t0)); % incompatible
            d = copy(b); d.fs = d.fs * 2; % incompatible
            v = randn([T 1 1 1]);
            u = randn([T 1 N 1]);

            for f = {@plus, @minus, @times, @rdivide, @ldivide}
                % should work
                f{1}(5, a);
                f{1}(a, 5);
                f{1}(a, b);

                % should fail
                tst.assertError(@() f{1}(a, c), "");
                tst.assertError(@() f{1}(a, d), "");
            end

            % should work
            -a;  +b; %#ok<VUNUS>
           

            % should fail
        end

        function transforming(tst)
            [T,M,N,F] = deal(16,8,4,2);
            chd = ChannelData('data', rand([T M N F]), 't0', randn([1 M 1 F]), 'order', 'TMNF', 'fs', 25);

            % time alignment
            chd.rectifyt0();
            chd.alignInt();
            chd.time = randn(size(chd.t0)) + swapdim((0 : chd.T - 1) ./ chd.fs, 2, chd.tdim);

            % time reversal and convolution
            wv = Waveform('t0', -1/32, 'tend', 1/8, 'fs', chd.fs, 'fun', @(t) cospi(2j*pi*100*t));
            wvr = reverse(wv); 
            tst.assertEqual(-wvr.tend, wv.t0);
            tst.assertEqual(-wv.tend, wvr.t0);
            tst.assertEqual(wv.samples, reverse(reverse(wv)).samples);
            for s = ["full", "same", "valid"], convt(chd, wv, s), end
            
            % join/splice
            ChannelData.empty().join(4)
            join(chd.splice(),4)
            chd.splice(4)
            chd.splice(1,4)

            % data order
            ord = [1,4,2,3];
            chdp = chd.permuteD(ord);
            chdi = chdp.ipermuteD(ord);
            tst.assertEqual(chd.order, chdi.order);
            tst.assertEqual(chd.data , chdi.data );
            tst.assertEqual(chd.t0   , chdi.t0   );

            % sampling
            tst.assertEqual(downsample(chd,2).fs, chd.fs/2);
            tst.assertEqual(subD(chd, 1:2:chd.T, chd.tdim).fs, chd.fs/2);
            tst.assertEqual(subD(chd, {1:4:chd.T, [2,3,4]}, [chd.tdim, chd.ndim]).fs, chd.fs/4);

            tst.assertError(@() subD(chd, {[1 2 4 7]}, chd.tdim), 'QUPS:ChannelData:nonuniformTemporalIndexing'); % bad spacing

        end
    end
end