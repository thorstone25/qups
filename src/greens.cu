// # include "helper_math.h" // vector math

// # include "sizes.cu" // size defines

# include "interpd.cu" // samplers using constant sizing

// # include "half2_math.h" // vector math for half types only 
template<typename T2, typename U, typename U3>
__device__ void greens_temp(T2 * __restrict__ y, 
    const U * __restrict__ Pi, const T2 * __restrict__ a, 
    const U * __restrict__ Pr, const U * __restrict__ Pv, 
    const T2 * __restrict__ x, const U s0,
	const U t0, const U fs, const U cinv,
    const int * E, const int iflag
    ) {

    // get starting index of this scatterer
    const size_t s = threadIdx.x + blockIdx.x * blockDim.x; // time 
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y; // rx
    const size_t m = threadIdx.z + blockIdx.z * blockDim.z; // tx

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const U3 * pi = reinterpret_cast<const U3*>(Pi); // 3 x I
    const U3 * pr = reinterpret_cast<const U3*>(Pr); // 3 x N x E
    const U3 * pv = reinterpret_cast<const U3*>(Pv); // 3 x M x E

    // rename for readability
    const size_t N = QUPS_N, I = QUPS_I, S = QUPS_S, M = QUPS_M;//, T = QUPS_T;
    // rxs, num scat, output time size, txs, kernel time size, 
    // S is size of output, T is size of input kernel, I the number of scats
    
    // temp vars
    // const float ts = s + s0*fs; // compute the time index for this thread
    const T2 zero_v = {0, 0}; // OOB value
    U r, tau;
    T2 val = zero_v;

    // if valid scat, for each tx/rx
    if(s < S){
        #pragma unroll
        for(size_t i = 0; i < I; ++i){ // for each scatterer
            # pragma unroll 
            for(size_t me = 0; me < E[1]; ++me){ // for each tx sub-aperture
                # pragma unroll 
                for(size_t ne = 0; ne < E[0]; ++ne){ // for each rx sub-aperture

                    // 2-way path time
                    r = cinv * (length(pi[i] - pr[n + ne*N]) + length(pi[i] - pv[m + me*M])); // (virtual) transmit to pixel vector
                    
                    // get kernel delay for the scatterer
                    tau = (float)s + (s0 - r - t0)*fs;
                    
                    // sample the kernel and add to the signal at this time
                    val += a[i] * sample(x, tau, iflag, zero_v); // out of bounds: extrap 0            
                }
            }
        }

        // output signal when all scatterers and sub-apertures are sampled
         y[s + n*S + m*N*S] = val;
    }
}

__global__ void greensf(float2 * __restrict__ y, 
    const float * __restrict__ Pi, const float2 * __restrict__ a, 
    const float * __restrict__ Pr, const float * __restrict__ Pv, 
    const float2 * __restrict__ x, const float s0,
	const float t0, const float fs, const float cinv,
    const int * E, const int iflag
    ) {
    greens_temp<float2, float, float3>(y, Pi, a, Pr, Pv, x, s0, t0, fs, cinv, E, iflag);
}

__global__ void greens(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double2 * __restrict__ a, 
    const double * __restrict__ Pr, const double * __restrict__ Pv, 
    const double2 * __restrict__ x, const double s0,
	const double t0, const double fs, const double cinv,
    const int * E, const int iflag
    ) {
    greens_temp<double2, double, double3>(y, Pi, a, Pr, Pv, x, s0, t0, fs, cinv, E, iflag);
}

#if (__CUDA_ARCH__ >= 530)

__global__ void greensh(ushort2 * __restrict__ y, 
    const float * __restrict__ Pi, const short2 * __restrict__ a, 
    const float * __restrict__ Pr, const float * __restrict__ Pv, 
    const ushort2 * __restrict__ x, const float s0,
	const float t0, const float fs, const float cinv,
    const int * E, const int iflag
    ) {
    greens_temp<half2, float, float3>((half2 *)y, Pi, (const half2 *)a, Pr, Pv, (const half2 *)x, s0, t0, fs, cinv, E, iflag);
}
#endif

/*

__global__ void greensf(float2 * __restrict__ y, 
    const float * __restrict__ Pi, const float * __restrict__ a, 
    const float * __restrict__ Pr, const float * __restrict__ Pv, 
    const float2 * __restrict__ x, const float s0,
	const float t0, const float fs, const float cinv,
    const int * E, const int iflag
    ) {

    // get starting index of this scatterer
    const size_t s = threadIdx.x + blockIdx.x * blockDim.x; // time 
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y; // rx
    const size_t m = threadIdx.z + blockIdx.z * blockDim.z; // tx

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 3 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N x E
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M x E

    // rename for readability
    const size_t N = QUPS_N, I = QUPS_I, S = QUPS_S, M = QUPS_M;//, T = QUPS_T;
    // rxs, num scat, output time size, txs, kernel time size, 
    // S is size of output, T is size of input kernel, I the number of scats
    
    // temp vars
    // const float ts = s + s0*fs; // compute the time index for this thread
    const float2 zero_v = make_float2(0.0f); // OOB value
    float r, tau;
    float2 val = zero_v;

    // if valid scat, for each tx/rx
    if(s < S){
        #pragma unroll
        for(size_t i = 0; i < I; ++i){ // for each scatterer
            # pragma unroll 
            for(size_t me = 0; me < E[1]; ++me){ // for each tx sub-aperture
                # pragma unroll 
                for(size_t ne = 0; ne < E[0]; ++ne){ // for each rx sub-aperture

                    // 2-way path time
                    r = cinv * (length(pi[i] - pr[n + ne*N]) + length(pi[i] - pv[m + me*M])); // (virtual) transmit to pixel vector
                    
                    // get kernel delay for the scatterer
                    tau = (float)s + (s0 - r - t0)*fs;
                    
                    // sample the kernel and add to the signal at this time
                    val += a[i] * sample(x, tau, iflag, zero_v); // out of bounds: extrap 0            
                }
            }
        }

        // output signal when all scatterers and sub-apertures are sampled
         y[s + n*S + m*N*S] = val;
    }
}

__global__ void greens(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double * __restrict__ a, 
    const double * __restrict__ Pr, const double * __restrict__ Pv, 
    const double2 * __restrict__ x, const double s0,
	const double t0, const double fs, const double cinv,
    const int * E, const int iflag
    ) {

    // get starting index of this scatterer
    const size_t s = threadIdx.x + blockIdx.x * blockDim.x; // time 
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y; // rx
    const size_t m = threadIdx.z + blockIdx.z * blockDim.z; // tx

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 3 x I
    const double3 * pr = reinterpret_cast<const double3*>(Pr); // 3 x N x E
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M x E

    // rename for readability
    const size_t N = QUPS_N, I = QUPS_I, S = QUPS_S, M = QUPS_M;//, T = QUPS_T;
    // rxs, num scat, output time size, txs, kernel time size, 
    // S is size of output, T is size of input kernel, I the number of scats
    
    // temp vars
    // const double ts = s + s0*fs; // compute the time index for this thread
    const double2 zero_v = make_double2(0.0f); // OOB value
    double r, tau;
    double2 val = zero_v;

    // if valid scat, for each tx/rx
    if(s < S){
        #pragma unroll
        for(size_t i = 0; i < I; ++i){ // for each scatterer
            # pragma unroll 
            for(size_t me = 0; me < E[1]; ++me){ // for each tx sub-aperture
                # pragma unroll 
                for(size_t ne = 0; ne < E[0]; ++ne){ // for each rx sub-aperture

                    // 2-way path time
                    r = cinv * (length(pi[i] - pr[n + ne*N]) + length(pi[i] - pv[m + me*M])); // (virtual) transmit to pixel vector
                    
                    // get kernel delay for the scatterer
                    tau = (double)s + (s0 - r - t0)*fs;
                    
                    // sample the kernel and add to the signal at this time
                    val += a[i] * sample(x, tau, iflag, zero_v); // out of bounds: extrap 0            
                }
            }
        }

        // output signal when all scatterers and sub-apertures are sampled
         y[s + n*S + m*N*S] = val;
    }
}

#if (__CUDA_ARCH__ >= 530)
__global__ void greensh(ushort2 * __restrict__ y, 
    const float * __restrict__ Pi, const unsigned short * __restrict__ a, 
    const float * __restrict__ Pr, const float * __restrict__ Pv, 
    const ushort2 * __restrict__ x, const float s0,
	const float t0, const float fs, const float cinv,
    const int * E, const int iflag
    ) {

    // get starting index of this scatterer
    const size_t s = threadIdx.x + blockIdx.x * blockDim.x; // time 
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y; // rx
    const size_t m = threadIdx.z + blockIdx.z * blockDim.z; // tx

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 3 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N x E
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M x E

    // rename for readability
    const size_t N = QUPS_N, I = QUPS_I, S = QUPS_S, M = QUPS_M;//, T = QUPS_T;
    // rxs, num scat, output time size, txs, kernel time size, 
    // S is size of output, T is size of input kernel, I the number of scats
    
    // temp vars
    // const double ts = s + s0*fs; // compute the time index for this thread
    const half2 zero_v = make_half2(0.0f, 0.0f); // OOB value
    float r, tau;
    half2 val = zero_v;

    // if valid scat, for each tx/rx
    if(s < S){
        #pragma unroll
        for(size_t i = 0; i < I; ++i){ // for each scatterer
            # pragma unroll 
            for(size_t me = 0; me < E[1]; ++me){ // for each tx sub-aperture
                # pragma unroll 
                for(size_t ne = 0; ne < E[0]; ++ne){ // for each rx sub-aperture

                    // 2-way path time
                    r = cinv * (length(pi[i] - pr[n + ne*N]) + length(pi[i] - pv[m + me*M])); // (virtual) transmit to pixel vector
                    
                    // get kernel delay for the scatterer
                    tau = (float)s + (s0 - r - t0)*fs;
                    
                    // sample the kernel and add to the signal at this time
                    val += u2h(a[i]) * sample((half2 *)x, (half)tau, iflag, zero_v); // out of bounds: extrap 0            
                }
            }
        }

        // output signal when all scatterers and sub-apertures are sampled
         y[s + n*S + m*N*S] = h2u(val);
    }
}
#endif
*/
