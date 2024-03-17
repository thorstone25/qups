// # include "helper_math.h" // vector math

// # include "sizes.cu" // size defines

# include "interpd.cu" // samplers using constant sizing

// # include "half2_math.h" // vector math for half types only 
template<typename T2, typename U, typename U3>
__device__ void greens_temp(T2 * __restrict__ y, 
    const U * __restrict__ Pi, const T2 * __restrict__ a, 
    const U * __restrict__ Pr, const U * __restrict__ Pv, 
    const T2 * __restrict__ x, const U * __restrict__ sb, const size_t * iblock,
	const U s0t0fscinv[6],
    const int * E, const int iflag
    ) {

    // extract time parameters
    const U s0   = s0t0fscinv[0];
    const U t0   = s0t0fscinv[1];
    const U fs   = s0t0fscinv[2];
    const U fsr  = s0t0fscinv[3];
    const U cinv = s0t0fscinv[4];
    const U R0   = s0t0fscinv[5]; // min distance for div by 0
    
    // get starting index of this scatterer
    const size_t s = threadIdx.x + blockIdx.x * blockDim.x; // time 
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y; // rx
    const size_t m = threadIdx.z + blockIdx.z * blockDim.z; // tx

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const U3 * pi = reinterpret_cast<const U3*>(Pi); // 3 x I
    const U3 * pr = reinterpret_cast<const U3*>(Pr); // 3 x N x E
    const U3 * pv = reinterpret_cast<const U3*>(Pv); // 3 x M x E

    // rename for readability
    const size_t N = QUPS_N, S = QUPS_S, M = QUPS_M, I = QUPS_I;//, T = QUPS_T;
    // rxs, num scat, output time size, txs, kernel time size, 
    // S is size of output, T is size of input kernel, I the number of scats
    
    // temp vars
    // const float ts = s + s0*fs; // compute the time index for this thread
    const T2 zero_v = {0, 0}; // OOB value, power scaling
    // const U A0 = (R0) ? 1 / (R0 * R0) : 1; // max power scaling
    U r1, r2, tau; // length, time (tmp values)
    T2 val = zero_v; // accumulator

    // if valid scat, for each tx/rx
    if(s < S){
        for(size_t i = iblock[2*blockIdx.x+0]; i <= iblock[2*blockIdx.x+1] && i < I; ++i){ // for each scatterer
            if(s >= sb[2*i+0]){ // if within sampling window
                # pragma unroll 
                for(size_t me = 0; me < E[1]; ++me){ // for each tx sub-aperture
                    # pragma unroll 
                    for(size_t ne = 0; ne < E[0]; ++ne){ // for each rx sub-aperture
    
                        // 2-way path distance
                        r1 = length(pi[i] - pr[n + ne*N]); // rx propagation
                        r2 = length(pi[i] - pv[m + me*M]); // tx propagation
                        
                        // get kernel delay for the scatterer
                        tau = (U)s - (cinv * (r1 + r2) + t0 - s0)*fs;

                        // limit to minimum dist of R0
                        if(R0){
                            r1 /= R0; r1 = (r1 < 1) ? 1 : r1;
                            r2 /= R0; r2 = (r2 < 1) ? 1 : r2;
                        } else {
                            r1 = 1;
                            r2 = 1;
                        }
                        
                        // sample the kernel and add to the signal at this time
                        // fsr applies a 'stretch' operation to the sample time, because the 
                        // input data x is sampled at sampling frequency fsr * fs
                        val += a[i] * sample(x, fsr*tau, iflag, zero_v) / (r1 * r2); // out of bounds: extrap 0
                    }
                }
            }
        }
        
        // output signal when all scatterers and sub-apertures are sampled
        // normalize by the discrete length of the signal, and the upsampling factor
        y[s + n*S + m*N*S] =  val / (R0*R0*fsr);
    }
}

__global__ void greensf(float2 * __restrict__ y, 
    const float * __restrict__ Pi, const float2 * __restrict__ a, 
    const float * __restrict__ Pr, const float * __restrict__ Pv, 
    const float2 * __restrict__ x, const float * __restrict__ sb, const size_t * iblock,
	const float s0t0fscinv[6],
    const int * E, const int iflag
    ) {
    greens_temp<float2, float, float3>(y, Pi, a, Pr, Pv, x, sb, iblock, s0t0fscinv, E, iflag);
}

__global__ void greens(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double2 * __restrict__ a, 
    const double * __restrict__ Pr, const double * __restrict__ Pv, 
    const double2 * __restrict__ x, const double * __restrict__ sb, const size_t * iblock,
	const double s0t0fscinv[6],
    const int * E, const int iflag
    ) {
    greens_temp<double2, double, double3>(y, Pi, a, Pr, Pv, x, sb, iblock, s0t0fscinv, E, iflag);
}

#if (__CUDA_ARCH__ >= 530)
inline __host__ __device__ half2 operator/(half2 a, float b){
    return make_half2((float)a.x / b, (float)a.y / b); // scale amplitude in FP32 precision
}

__global__ void greensh(ushort2 * __restrict__ y, 
    const float * __restrict__ Pi, const short2 * __restrict__ a, 
    const float * __restrict__ Pr, const float * __restrict__ Pv, 
    const ushort2 * __restrict__ x, const float * __restrict__ sb,const size_t * iblock,
	const float s0t0fscinv[6],
    const int * E, const int iflag
    ) {
    greens_temp<half2, float, float3>((half2 *)y, Pi, (const half2 *)a, Pr, Pv, (const half2 *)x, sb, iblock, s0t0fscinv, E, iflag);
}
#endif
