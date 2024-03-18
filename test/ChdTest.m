classdef(TestTags = ["full", "Github", "build", "syntax"]) ChdTest < matlab.unittest.TestCase
    % CHDTEST - This class test that ChannelData methods function properly 
    
    %#ok<*NASGU> The outputs will not be televised

    methods(Test)
        function typeCheck(tst)
            chd = ChannelData('data', randn([16 8 4 2], 'single'));

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
            sparse(ChannelData('data', randn([16 32])))

            % unit conversion check
            chdo = angle(complex(chd));
            chdo = deg2rad(rad2deg(chd));
            chdo = mag2db(abs(chd));
            chdo = real(imag(complex(chd)));

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
            for D = [D1 D2]
                filter(chd, D)
                filtfilt(chd, D)
                fftfilt(chd, D)
            end

            % sampling
            downsample(downmix(resample(chd,5), 1/4), 2);
        end

        function arithmetic(tst)
            [T,M,N,F] = deal(16,8,4,2);
            a = ChannelData('data', rand([T M N F]), 't0', randn([1 M 1 F]), 'order', 'TMNF');
            b = ChannelData('data', rand([T M N F]), 't0', a.t0            , 'order', 'TMNF');

            5 + a, 
            a + 5, 
            a + b,
            5 .* a
            a .* 5
            a .* b
        end
    end
end