// # include "helper_math.h" // vector math

// # include "sizes.cu" // size defines

# include "interpolators.cl" // samplers using constant sizing

#   if QUPS_PRECISION == 32
typedef float4 U4;
# elif QUPS_PRECISION == 64 
typedef double4 U4;
# elif QUPS_PRECISION == 16 
typedef float4 U4;
# endif

kernel void greens(global T2 * __restrict__ y, 
    global const U4 * __restrict__ pi, global const T2 * __restrict__ a, 
    global const U4 * __restrict__ pr, global const U4 * __restrict__ pv, 
    global const T2 * __restrict__ x,  global const V * __restrict__ sb, global const ulong * iblock,
	global const V s0t0fscinv[5],
    global const int * E, const int iflag
    ) {

    // extract time parameters
    const V s0   = s0t0fscinv[0];
    const V t0   = s0t0fscinv[1];
    const V fs   = s0t0fscinv[2];
    const V fsr  = s0t0fscinv[3];
    const V cinv = s0t0fscinv[4];
    
    // get starting index of this scatterer
    const ulong s = get_global_id(0); // time 
    const ulong n = get_global_id(1); // rx
    const ulong m = get_global_id(2); // tx

    // rename for readability
    const ulong N = QUPS_N, S = QUPS_S, M = QUPS_M, I = QUPS_I; //, T = QUPS_T;
    // rxs, num scat, output time size, txs, kernel time size, 
    // S is size of output, T is size of input kernel, I the number of scats
    
    // temp vars
    // const float ts = s + s0*fs; // compute the time index for this thread
    const T2 zero_v = (T2)(0, 0), fpow = (T2)(fsr, fsr); // OOB value, power scaling
    V r, tau; // length, time (tmp values)
    T2 val = zero_v; // accumulator

    // if valid scat, for each tx/rx
    if(s < S){
        for(ulong i = iblock[2*get_group_id(0)+0]; i <= iblock[2*get_group_id(0)+1] && i < I; ++i){ // for each scatterer
            if(s >= sb[2*i+0]){ // if within sampling window
                # pragma unroll 
                for(ulong me = 0; me < E[1]; ++me){ // for each tx sub-aperture
                    # pragma unroll 
                    for(ulong ne = 0; ne < E[0]; ++ne){ // for each rx sub-aperture
    
                        // 2-way path distance
                        r = (length(pi[i] - pr[n + ne*N]) + length(pi[i] - pv[m + me*M])); // (virtual) transmit to pixel vector
                        
                        // get kernel delay for the scatterer
                        tau = (V)s - (cinv * r + t0 - s0)*fs;
                        
                        // sample the kernel and add to the signal at this time
                        // fsr applies a 'stretch' operation to the sample time, because the 
                        // input data x is sampled at sampling frequency fsr * fs
                        val += cmul(a[i], sample(x, fsr*tau, iflag, zero_v)); // out of bounds: extrap 0
                    }
                }
            }
        }
        
        // output signal when all scatterers and sub-apertures are sampled
        // normalize by the discrete length of the signal
        y[s + n*S + m*N*S] = val / fpow;
    }
}

