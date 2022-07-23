/**
 @file gpuBF/interp1d.cuh
 @author Dongwoon Hyun (dongwoon.hyun@stanford.edu)
 @date 2021-03-08

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

# include "helper_math.h" // vector math
        
# include "sizes.cu" // size defines

// Modified for use as a stand-alone ptx
// data is (T x N x F)
// sample times are (I x N x M)

#if (__CUDA_ARCH__ >= 530)
/// @brief Device function to reinterpret ushort values as half
inline __device__ half2 u2h(const ushort2 a){
    return __halves2half2(__ushort_as_half(a.x), __ushort_as_half(a.y));
}

inline __device__ half u2h(const ushort a){
    return __ushort_as_half(a);
}

/// @brief Device function to reinterpret half values as ushort
inline __device__ ushort2 h2u(const half2 a){
    return make_ushort2(__half_as_ushort(a.x), __half_as_ushort(a.y));
}

inline __device__ ushort u2h(const half a){
    return __half_as_ushort(a);
}

#endif

int __device__ modfi(const double x, double* r){
    double i; // temp var
    (*r) = modf(x, &i);  // x is the fractional part, i the integer part
    return (int) i;
}

int __device__ modfi(const float x, float* r){
    float i; // temp var
    (*r) = modff(x, &i);  // x is the fractional part, i the integer part
    return (int) i;
}

#if (__CUDA_ARCH__ >= 530)
int __device__ modfi(const half x, half* r){
    float i; // temp var
    (*r) = modff((float) x, &i);  // x is the fractional part, i the integer part
    return (int) i;
}
#endif

/// @brief Device function for nearest-neighbor interpolation
template<typename U, typename T2>
__device__ T2 nearest(const T2 * x, U tau, T2 no_v) {
    const int ti = (int) roundf(tau); // round to nearest integer
    return (0 <= ti && ti < QUPS_T) ? x[ti] : no_v;
}

/// @brief Device function for linear interpolation
template<typename U, typename T2>
__device__ T2 linear(const T2 * x, U tau, T2 no_v) {
    U tf;
    
    // fractional and integer part
    const int ti = modfi(tau, &tf); 
                
    // if in bounds, linearly interpolate by ratio tau at time-index ti[+1]
    return (0 <= ti && ti + 1 < QUPS_T) ? lerp(x[ti], x[ti + 1], tf) : no_v;
}

/// @brief Device function for cubic Hermite interpolation
template<typename U, typename T2>
__device__ T2 cubic(const T2 * x, U tau, T2 no_v) {
  U u;
  int ti = modfi(tau, &u); // u is the fractional part, ti the integer part
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;

  T2 s0 = x[ti - 1];
  T2 s1 = x[ti + 0];
  T2 s2 = x[ti + 1];
  T2 s3 = x[ti + 2];

  // Cubic Hermite interpolation (increased precision using fused multiply-adds)
  U a0 = 0 + u * (-1 + u * (+2 * u - 1));
  U a1 = 2 + u * (+0 + u * (-5 * u + 3));
  U a2 = 0 + u * (+1 + u * (+4 * u - 3));
  U a3 = 0 + u * (+0 + u * (-1 * u + 1));
  // // Cubic Hermite interpolation (naive, less precise implementation)
  // float a0 = -1 * u * u * u + 2 * u * u - 1 * u + 0;
  // float a1 = +3 * u * u * u - 5 * u * u + 0 * u + 2;
  // float a2 = -3 * u * u * u + 4 * u * u + 1 * u + 0;
  // float a3 = +1 * u * u * u - 1 * u * u + 0 * u + 0;
  return (s0 * a0 + s1 * a1 + s2 * a2 + s3 * a3) * 0.5f;
}

/// @brief Inline helper code for Lanczos 3-lobe interpolation
template<typename T>
__device__ T lanczos_helper(T u, int a) {
  const T PI = (atanf(1) * 4); // TODO: isn't there a built-in PI value?
  T ret;
  if (u == (T) 0.)
      ret = (T) 1.;
  else {
      const T spu  = sinpif(u    );
      const T spua = sinpif(u / a);
      ret = 2. * spu * spua / (PI * PI * u * u);
  }
  return ret;
  // return (u == (half)0.f) ? (half) 1.f : (half) 2.f * sinpif(__half2float(u)) * sinpif(__half2float(u) / a) / (PI * PI * u * u);
}

/// @brief Device function for Lanczos 3-lobe interpolation
template<typename U, typename T2>
__device__ T2 lanczos3(const T2 * x, U tau, T2 no_v) {
  constexpr int a = 2;  // a=2 for 3-lobe Lanczos resampling
  U u;
  const int ti = modfi(tau, &u); // u is the fractional part, ti the integer part
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;
  
  T2 s0 = x[ti - 1];
  T2 s1 = x[ti + 0];
  T2 s2 = x[ti + 1];
  T2 s3 = x[ti + 2];
  U a0 = lanczos_helper(u + 1, a);
  U a1 = lanczos_helper(u + 0, a);
  U a2 = lanczos_helper(u - 1, a);
  U a3 = lanczos_helper(u - 2, a);
  return s0 * a0 + s1 * a1 + s2 * a2 + s3 * a3;
}

template <typename U, typename T2>
__device__ T2 sample(const T2 * x, U tau, int flag, const T2 no_v){
    // sample according to the flag
         if (flag == 0)
        return nearest  (x, tau, no_v);
    else if (flag == 1)
        return linear   (x, tau, no_v);
    else if (flag == 2)
        return cubic    (x, tau, no_v);
    else if (flag == 3)
        return lanczos3 (x, tau, no_v);
    else if (flag == 4)
        return linear   (x, tau, no_v);
    else 
        return no_v;     
}

template<typename T2, typename U, typename V> // channel data type, time data type, time-sampling type
__device__ void interpd_temp(T2 * __restrict__ y, 
    const T2 * __restrict__ x, const U * __restrict__ tau, const int flag, const T2 no_v) {
    
    // get sampling index
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    // const size_t m = threadIdx.z + blockIdx.z * blockDim.z;
    
    // rename for readability
    const size_t I = QUPS_I, M = QUPS_M, N = QUPS_N, T = QUPS_T, F = QUPS_F;

    // remap indices
    const size_t i = tid % I;
    const size_t m = tid / I;

    // if valid sample, for each tx/rx
    if(i < I && n < N && m < M){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // per transmit
            y[i + n*I + m*N*I + f*M*N*I] = sample(&x[n*T + f*N*T], (V)(tau[i + n*I + m*I*N]), flag, no_v);
        }
    }
}

template<typename T2, typename U, typename V> // channel data type, time data type, time-sampling type
__device__ void wsinterpd_temp(T2 * __restrict__ y, 
    const T2 * __restrict__ w, const T2 * __restrict__ x, 
    const U * __restrict__ tau, const size_t * sizes, 
    const size_t * wstride, const size_t * ystride, const int flag, 
    const half2 no_v) {

    // get sampling index
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    // const size_t m = threadIdx.z + blockIdx.z * blockDim.z;
    
    // rename for readability
    const size_t I = QUPS_I, M = QUPS_M, N = QUPS_N, T = QUPS_T, F = QUPS_F, S = QUPS_S;
    size_t k,l,sz; // weighting / output indexing

    // remap indices
    const size_t i = tid % I;
    const size_t m = tid / I;

    // if valid sample, per each i,n,f
    if(i < I && n < N && m < M){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // for m
            // global index
            const size_t j = (i + n*I + m*N*I + f*N*I*M);

            // get weight vector and output indices
            k = 0, l = 0; 
            # pragma unroll
            for(size_t s = 0; s < S; ++s){ // for each dimension s
                // calculate the indexing stride for dimension s i.e. 
                // size of all prior dimensions
                sz = 1;
                for(size_t sp = 0; sp < s; ++sp)
                    sz *= sizes[sp]; 

                const size_t js = ((j / sz) % sizes[s]); // index for this dimension
                k += js * wstride[s]; // add pitched index for this dim (weights)
                l += js * ystride[s]; // add pitched index for this dim (outputs)
            }

            const T2 val = w[k] * sample(&x[n*T + f*N*T], (V)tau[i + n*I + m*I*N], flag, no_v); // weighted sample
            atomicAdd(&y[l], val); // store
        }
    }
}


#if (__CUDA_ARCH__ >= 530)

__global__ void interpdh(ushort2 * __restrict__ y, 
    const ushort2 * __restrict__ x, const unsigned short * __restrict__ tau, const int flag) {
    interpd_temp<half2, half, float>(
        (half2 *)y, (const half2 *)x, (const half *)tau, flag, (const half2) make_half2(0,0)
    );
}

__global__ void wsinterpd_temp(ushort2 * __restrict__ y, 
    const ushort2 * __restrict__ w, const ushort2 * __restrict__ x, 
    const unsigned short * __restrict__ tau, const size_t * sizes, 
    const size_t * wstride, const size_t * ystride, const int flag, 
    const ushort2 no_v) {
    wsinterpd_temp<half2, half, float>((half2 *)y, (const half2 *)w, (const half2 *)x,
             (half *)tau, sizes, wstride, ystride, flag, (const half2) make_half2(0,0));
}
#endif
/*
    // get sampling index
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    // const size_t m = threadIdx.z + blockIdx.z * blockDim.z;
    const half2 no_v = make_half2(0.0f, 0.0f);

    // rename for readability
    const size_t I = QUPS_I, M = QUPS_M, N = QUPS_N, T = QUPS_T, F = QUPS_F;

    // remap indices
    const size_t i = tid % I;
    const size_t m = tid / I;

    // if valid sample, for each tx/rx
    if(i < I && n < N && m < M){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // per transmit
            y[i + n*I + m*N*I + f*M*N*I] = h2u(sample((const half2 *)&x[n*T + f*N*T], (float)u2h(tau[i + n*I + m*I*N]), flag, no_v));
        }
    }
}

__global__ void wsinterpdh(ushort2 * __restrict__ y, 
    const ushort2 * __restrict__ w, const ushort2 * __restrict__ x, 
    const unsigned short * __restrict__ tau, const size_t * sizes, 
    const size_t * wstride, const size_t * ystride, const int flag
    ) {

    // get sampling index
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    // const size_t m = threadIdx.z + blockIdx.z * blockDim.z;
    const half2 no_v = make_half2(0.0f, 0.0f);

    // rename for readability
    const size_t I = QUPS_I, M = QUPS_M, N = QUPS_N, T = QUPS_T, F = QUPS_F, S = QUPS_S;
    size_t k,l,sz; // weighting / output indexing

    // remap indices
    const size_t i = tid % I;
    const size_t m = tid / I;

    // cast MATLAB alias to CUDA half type
    half2 * yh = reinterpret_cast<half2 *>(y);

    // if valid sample, per each i,n,f
    if(i < I && n < N && m < M){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // for m
            // global index
            const size_t j = (i + n*I + m*N*I + f*N*I*M);

            // get weight vector and output indices
            k = 0, l = 0; 
            # pragma unroll
            for(size_t s = 0; s < S; ++s){ // for each dimension s
                // calculate the indexing stride for dimension s i.e. 
                // size of all prior dimensions
                sz = 1;
                for(size_t sp = 0; sp < s; ++sp)
                    sz *= sizes[sp]; 

                const size_t js = ((j / sz) % sizes[s]); // index for this dimension
                k += js * wstride[s]; // add pitched index for this dim (weights)
                l += js * ystride[s]; // add pitched index for this dim (outputs)
            }

            const half2 val = u2h(w[k]) * sample(&x[n*T + f*N*T], u2h(tau[i + n*I + m*I*N]), flag, no_v); // weighted sample
            atomicAdd(&yh[l], val); // store
        }
    }
}

*/
__global__ void interpdf(float2 * __restrict__ y, 
    const float2 * __restrict__ x, const float * __restrict__ tau, const int flag) {

    // get sampling index
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    // const size_t m = threadIdx.z + blockIdx.z * blockDim.z;
    float2 no_v = make_float2(0.0f, 0.0f);

    // rename for readability
    const size_t I = QUPS_I, M = QUPS_M, N = QUPS_N, T = QUPS_T, F = QUPS_F;

    // remap indices
    const size_t i = tid % I;
    const size_t m = tid / I;
    
    // if valid sample, for each time index
    if(i < I && n < N && m < M){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // per frame of data
            y[i + n*I + m*N*I + f*M*N*I] = sample(&x[n*T + f*N*T], tau[i + n*I + m*I*N], flag, no_v);
        }
    }
}

__global__ void wsinterpdf(float2 * __restrict__ y, 
    const float2 * __restrict__ w, const float2 * __restrict__ x, 
    const float * __restrict__ tau, const size_t * sizes, 
    const size_t * wstride, const size_t * ystride, const int flag
    ) {

    // get sampling index
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    // const size_t m = threadIdx.z + blockIdx.z * blockDim.z;
    float2 no_v = make_float2(0.0f, 0.0f);

    // rename for readability
    const size_t I = QUPS_I, M = QUPS_M, N = QUPS_N, T = QUPS_T, F = QUPS_F, S = QUPS_S;
    size_t k,l,sz; // weighting / output indexing

    // remap indices
    const size_t i = tid % I;
    const size_t m = tid / I;
    
    // if valid sample, per each i,n,f
    if(i < I && n < N && m < M){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // for m
            // global index
            const size_t j = (i + n*I + m*N*I + f*N*I*M);

            // get weight vector and output indices
            k = 0, l = 0; 
            # pragma unroll
            for(size_t s = 0; s < S; ++s){ // for each dimension s
                // calculate the indexing stride for dimension s i.e. 
                // size of all prior dimensions
                sz = 1;
                for(size_t sp = 0; sp < s; ++sp)
                    sz *= sizes[sp]; 

                const size_t js = ((j / sz) % sizes[s]); // index for this dimension
                k += js * wstride[s]; // add pitched index for this dim (weights)
                l += js * ystride[s]; // add pitched index for this dim (outputs)
            }

            const float2 val = w[k] * sample(&x[n*T + f*N*T], tau[i + n*I + m*I*N], flag, no_v); // weighted sample
            atomicAdd(&y[l].x, val.x); // store
            atomicAdd(&y[l].y, val.y); // store
        }
    }
}

__global__ void interpd(double2 * __restrict__ y, 
    const double2 * __restrict__ x, const double * __restrict__ tau, const int flag) {

    // get sampling index
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    // const size_t m = threadIdx.z + blockIdx.z * blockDim.z;
    double2 no_v = make_double2(0.0, 0.0);

    // rename for readability
    const size_t I = QUPS_I, M = QUPS_M, N = QUPS_N, T = QUPS_T, F = QUPS_F;

    // remap indices
    const size_t i = tid % I;
    const size_t m = tid / I;

    // if valid sample, for each sample index
    if(i < I && n < N && m < M){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // loop over data
            y[i + n*I + m*N*I + f*M*N*I] = sample(&x[n*T + f*N*T], tau[i + n*I + m*I*N], flag, no_v);
        }
    }
}

__global__ void wsinterpd(double2 * __restrict__ y, 
    const double2 * __restrict__ w, const double2 * __restrict__ x, 
    const double * __restrict__ tau, const size_t * sizes, 
    const size_t * wstride, const size_t * ystride, const int flag
    ) {

    // get sampling index
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    // const size_t f = threadIdx.z + blockIdx.z * blockDim.z;
    double2 no_v = make_double2(0.0, 0.0);
    
    // rename for readability
    const size_t I = QUPS_I, M = QUPS_M, N = QUPS_N, T = QUPS_T, F = QUPS_F, S = QUPS_S;
    size_t k,l,sz; // weighting / output indexing

    // remap indices
    const size_t i = tid % I;
    const size_t m = tid / I;
    
    // if valid sample, per each i,n,m (indexing the time delays)
    if(i < I && n < N && m < M){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // loop over the frames of data
            // global index
            const size_t j = (i + n*I + m*N*I + f*N*I*M);

            // get weight vector and output indices
            k = 0, l = 0; 
            # pragma unroll
            for(size_t s = 0; s < S; ++s){ // for each dimension s
                // calculate the indexing stride for dimension s i.e. 
                // size of all prior dimensions
                sz = 1;
                for(size_t sp = 0; sp < s; ++sp)
                    sz *= sizes[sp]; 

                const size_t js = ((j / sz) % sizes[s]); // index for this dimension
                k += js * wstride[s]; // add pitched index for this dim (weights)
                l += js * ystride[s]; // add pitched index for this dim (outputs)
            }

            const double2 val = w[k] * sample(&x[n*T + f*N*T], tau[i + n*I + m*I*N], flag, no_v); // weighted sample
            atomicAdd(&y[l].x, val.x); // store
            atomicAdd(&y[l].y, val.y); // store
        }
    }
}

