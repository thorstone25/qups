// # include "helper_math.h" // vector math

// # include "sizes.cu" // size defines

# include "interpolators.cl" // samplers using constant sizing

// # include "half2_math.h" // vector math for half types only 

/*
* Delay and sum the given data at the given pixels
*
* Given a set of pixels, (virtual or plane wave) transmitter locations, 
* receiver locations, as well as a datacube equipped with a time, 
* transmitter and receiver axis, perform simple delay-and-sum beamforming. 
* The data is linearly interpolated at the sample time. An image is 
* generated for each receiver element. Summation across the transmitters 
* and receivers is implicit.
*
* All positions are in vector coordinates. 
* 
* If the virtual transmitter 
* normal has a fourth component that is 0, this indicates that the 
* transmission should be treated as a plane wave transmission instead of a 
* virtual source transmission. 
* 
* The value of t = 0 must be the time when the peak of the wavefront 
* reaches the virtual source location. Because this time must be the same 
* for all transmits, the datacube must be stitched together in such a way 
* that for all transmits, the same time axis is used
*
* Inputs:
*  y:           complex pixel values per channel (I)
*  Pi:          pixel positions (3 x I)
*  Pr:          receiver positions (3 x N)
*  Pv:          (virtual) transmitter positions (3 x M)
*  Nv:          (virtual) transmitter normal (3 x M)
*  x:           datacube of complex sample values (T x M x N)
*  t0:          initial time for the data
*  fs:          sampling frequency of the data
*  cinv:        inverse of the speed of sound used for beamforming
* 
* I -> pixels, M -> transmitters, N -> receivers, T -> time samples
*
*/

#   if QUPS_INTERPD_PRECISION == 32
typedef float4 U4;
# elif QUPS_INTERPD_PRECISION == 64 
typedef double4 U4;
# elif QUPS_INTERPD_PRECISION == 16 
typedef float4 U4;
# endif

kernel void DAS(volatile global T2 * y, 
    const global U4 * pi, const global U4 * pr,
    const global U4 * pv, const global U4 * nv, 
    const global T2 * a,  const global V  * cinv, const global ulong * acstride, 
    const global T2 * x, const int flag, const global V t0fsfc[2]) {
    
    // unpack inputs
    #ifndef QUPS_BF_FLAG
    const int QUPS_BF_FLAG = flag; // set only if not fixed at compile time
    #endif
    // const V t0   = t0fsfc[0]; // start time
    const V fs   = t0fsfc[0]; // sampling frequency
    const V fc   = t0fsfc[1]; // modulation frequency
    const global ulong * astride = acstride;
    const global ulong * cstride = acstride + 5;

    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, T = QUPS_T, I = QUPS_I;
    const size_t I1 = QUPS_I1, I2 = QUPS_I2, I3 = QUPS_I3;
        
    // get starting index of this pixel
    const size_t tid = get_global_id(0); // pixel index
    const size_t kI  = get_global_size(0); // number of pixels per call

    // temp vars
    const T2 zero_v = (T2)(0, 0);
    T2 w = (T2)(1, 0);
    V dv, dr, tau;
    T2 val, pix;
    U4 rv = (U4)(0);

    // if valid pixel, for each tx/rx
    for(size_t i = tid; i < I; i += kI){        
        // get image coordinates
        const size_t i1 = (i             % QUPS_I1); // index in I1
        const size_t i2 = (i /  QUPS_I1) % QUPS_I2 ; // index in I2
        const size_t i3 = (i / (I1 * I2) % QUPS_I3); // index in I3
        const size_t abase = i1 * astride[0] + i2 * astride[1] + i3 * astride[2]; // base index for this pixel
        const size_t cbase = i1 * cstride[0] + i2 * cstride[1] + i3 * cstride[2]; // base index for this pixel

        // reset accumulator
        pix = zero_v;

        for(size_t m = 0; m < M; ++m){
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                // const U3 pvm = {pv[m].x,pv[m].y,pv[m].z}; // declared for MSVC (2019)
                rv.xyz = pi[i].xyz - pv[m].xyz; // (virtual) transmit to pixel vector

                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), (QUPS_DV ? 1.f : dot(rv.xyz, nv[m].xyz))) // virtual source
                    : dot(rv.xyz, nv[m].xyz); // plane wave

                dr = length(pi[i].xyz - pr[n].xyz); // rx path length

                // data/time index number
                const V ci = cinv[cbase + n * cstride[3] + m * cstride[4]];
                tau = (ci * (dv + dr) - pv[m].w);

                // apply demodulation if non zero
                if (fc) {w.x = cospi(2*fc*tau); w.y = sinpi(2*fc*tau);}

                // sample the trace
                val = sample(&x[(n + m * N) * T], tau * fs, flag & 7, zero_v); // out of bounds: extrap 0
                // const int t = (int) (tau * fs);
                // val = (0 <= t & t < T) ? x[(size_t) (t + (n + m * N) * T)] : zero_v;

                // apply apodization (requires complex multiplication)
                val = cmul(val, cmul(w, a[abase + n * astride[3] + m * astride[4]]));

                // choose the accumulation
                const int sflag = ((int)QUPS_BF_FLAG) & 24; // extract bits 5,4
                if(sflag == 8)
                    y[i + n*I] += val; // sum over tx, store over rx
                else if (sflag == 16)
                    y[i + m*I] += val; // sum over rx, store over tx
                else if (sflag == 24)
                    y[i + n*I + m*N*I] = val; // store over tx/rx
                else
                    pix += val; // sum over all
            }
        }
        if (!(((int)QUPS_BF_FLAG) & 24)) y[i] = pix; // output value here if accumulating over all
    }
}
