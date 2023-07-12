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

// atomicAdd natively supported @ CC 6.x+
#if __CUDA_ARCH__ < 600
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

// half type supported @ CC 5.3+
#if (__CUDA_ARCH__ >= 530) 

/// @brief Device function to reinterpret ushort values as half
inline __device__ half2 ui2h(const unsigned int i){
    union {
        unsigned int i;
        half2 h;
    } v;
    v.i = i;
    return __halves2half2(__ushort_as_half(v.h.x), __ushort_as_half(v.h.y));
}

inline __device__ unsigned int h2ui(const half2 a){
    union {
        unsigned int i;
        half2 h;
    } v;
    v.h = a;
    return v.i;
}

__device__ void atomicAdd(half2* address, half2 val)
{
    unsigned int * address_as_ull = (unsigned int*) address;
    unsigned int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, h2ui(val + ui2h(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    // return 	__ushort_as_half(old);
}
#endif
#endif

// half type supported @ CC 5.3+
#if (__CUDA_ARCH__ >= 530)
inline __device__ void atomicAddStore(half2 * y, const half2 val){
    atomicAdd(y, val);
}
#endif
inline __device__ void atomicAddStore(float2 * y, const float2 val){
    atomicAdd(&y[0].x, val.x);
    atomicAdd(&y[0].y, val.y);
}
inline __device__ void atomicAddStore(double2 * y, const double2 val){
    atomicAdd(&y[0].x, val.x);
    atomicAdd(&y[0].y, val.y);
}


__device__ size_t global_offset(size_t * dind, const size_t * sizes, const char * iflags){
// global index
    // init
    size_t dsz[3] = {1,1,1}; // {I,N,F} index cumulative sizes
    size_t sz = 1, j = 0; 

    // find offset
    for(size_t s = 0; s < QUPS_S; ++s){
        const char iflg = iflags[s]; // which label
        dsz[iflg] *= sizes[s]; // increase size for this label
        j += sz * (dind[iflg] %  dsz[iflg]); // add offset
                   dind[iflg] /= dsz[iflg]; // fold index
        sz *= sizes[s]; // increase indexing stride
    }

    return j;
}

template<typename T2, typename U, typename V> // channel data type, time data type, time-sampling type
__device__ void wsinterpd_temp(T2 * __restrict__ y, 
    const T2 * __restrict__ w, const T2 * __restrict__ x, 
    const U * __restrict__ tau, const size_t * sizes, 
    const char * iflags, const size_t * dstride, const int flag, 
    const T2 no_v, const U omega) {

    // get sampling index
    const size_t i = threadIdx.x + blockIdx.x * blockDim.x; // tid
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y + gridDim.y * (threadIdx.z + blockIdx.z * blockDim.z);
    // const size_t m = ;
    
    // rename for readability
    const size_t I = QUPS_I,N = QUPS_N, F = QUPS_F, S = QUPS_S, T = QUPS_T; // , M = QUPS_M
    size_t u,v,k,l,sz; // weighting / output indexing

    // remap indices
    // const size_t i = tid % I;
    // const size_t m = tid / I;

    // if valid sample, per each i,n,f
    if(i < I && n < N){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // for f
            // global index
            size_t dind[3] = {i,n,f};
            const size_t j = global_offset(dind, sizes, iflags);

            // get weight vector and output indices
            k = 0, l = 0; u = 0; v = 0;
            # pragma unroll
            for(size_t s = 0; s < S; ++s){ // for each dimension s
                // calculate the indexing stride for dimension s i.e. 
                // size of all prior dimensions
                sz = 1;
                for(size_t sp = 0; sp < s; ++sp)
                    sz *= sizes[sp]; 

                const size_t js = ((j / sz) % sizes[s]); // index for this dimension
                k += js * dstride[0 + 4*s]; // add pitched index for this dim (weights)
                l += js * dstride[1 + 4*s]; // add pitched index for this dim (outputs)
                u += js * dstride[2 + 4*s]; // add pitched index for this dim (time)
                v += js * dstride[3 + 4*s]; // add pitched index for this dim (samples)
            }
            
            const T2 a = {cosf(omega * tau[u]), sinf(omega * tau[u])}; // modulation phasor
            const T2 val = a * w[k] * sample(&x[(v)*T], (V)tau[u], flag, no_v); // weighted sample
            atomicAddStore(&y[l], val); // store
        }
    }
}

template<typename T2, typename U, typename V> // channel data type, time data type, time-sampling type
__device__ void wsinterpd2_temp(T2 * __restrict__ y, 
    const T2 * __restrict__ w, const T2 * __restrict__ x, 
    const U * __restrict__ tau1, const U * __restrict__ tau2, const size_t * sizes, 
    const char * iflags, const size_t * dstride, const int flag, 
    const T2 no_v, const U omega) {

    // get sampling index
    const size_t i = threadIdx.x + blockIdx.x * blockDim.x; // tid
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y + gridDim.y * (threadIdx.z + blockIdx.z * blockDim.z);
    // const size_t m = ;
    
    // rename for readability
    const size_t I = QUPS_I,N = QUPS_N, F = QUPS_F, S = QUPS_S, T = QUPS_T; // , M = QUPS_M
    size_t r,u,v,k,l,sz; // weighting / output indexing

    // remap indices
    // const size_t i = tid % I;
    // const size_t m = tid / I;

    // if valid sample, per each i,n,f
    if(i < I && n < N){
        # pragma unroll
        for(size_t f = 0; f < F; ++f){ // for f
            // global index
            size_t dind[3] = {i,n,f};
            const size_t j = global_offset(dind, sizes, iflags);

            // get weight vector and output indices
            k = 0, l = 0; r = 0; u = 0; v = 0;
            # pragma unroll
            for(size_t s = 0; s < S; ++s){ // for each dimension s
                // calculate the indexing stride for dimension s i.e. 
                // size of all prior dimensions
                sz = 1;
                for(size_t sp = 0; sp < s; ++sp)
                    sz *= sizes[sp]; 

                const size_t js = ((j / sz) % sizes[s]); // index for this dimension
                k += js * dstride[0 + 5*s]; // add pitched index for this dim (weights)
                l += js * dstride[1 + 5*s]; // add pitched index for this dim (outputs)
                r += js * dstride[2 + 5*s]; // add pitched index for this dim (time-1)
                u += js * dstride[3 + 5*s]; // add pitched index for this dim (time-2)
                v += js * dstride[4 + 5*s]; // add pitched index for this dim (samples)
            }
            const U  t = tau1[r] + tau2[u]; // time
            const T2 a = {cosf(omega * t), sinf(omega * t)}; // modulation phasor
            const T2 val = a * w[k] * sample(&x[(v)*T], (V)t, flag, no_v); // weighted sample
            atomicAddStore(&y[l], val); // store
        }
    }
}


/* interd kernels */
#if (__CUDA_ARCH__ >= 530)
__global__ void interpdh(ushort2 * __restrict__ y, const ushort2 * __restrict__ x, 
             const unsigned short * __restrict__ tau, const int flag) {
    interpd_temp<half2, half, float>(
        (half2 *)y, (const half2 *)x, (const half *)tau, flag, (const half2) make_half2(0,0)
    );
}
#endif

__global__ void interpdf(float2 * __restrict__ y, const float2 * __restrict__ x, 
                         const float * __restrict__ tau, const int flag) {
    interpd_temp<float2, float, float>(y, x, tau, flag, (const float2) make_float2(0,0));
}

__global__ void interpd(double2 * __restrict__ y, const double2 * __restrict__ x, 
                        const double * __restrict__ tau, const int flag) {
    interpd_temp<double2, double, double>(y, x, tau, flag, (const double2) make_double2(0,0));
}


/* wsinterpd kernels */
#if (__CUDA_ARCH__ >= 530)
__global__ void wsinterpdh(ushort2 * __restrict__ y, 
    const ushort2 * __restrict__ w, const ushort2 * __restrict__ x, 
    const unsigned short * __restrict__ tau, const size_t * sizes, const char * iflags, 
    const size_t * wstride, const int flag, const unsigned short omega) {
    wsinterpd_temp<half2, half, float>((half2 *)y, (const half2 *)w, (const half2 *)x,
             (half *)tau, sizes, iflags, wstride, flag, (const half2) make_half2(0,0), u2h(omega));
}
#endif

__global__ void wsinterpdf(float2 * __restrict__ y, 
    const float2 * __restrict__ w, const float2 * __restrict__ x, 
    const float * __restrict__ tau, const size_t * sizes, const char * iflags, 
    const size_t * wstride, const int flag, const float omega
    ) {
    wsinterpd_temp<float2, float, float>(y, w, x, tau, 
    sizes, iflags, wstride, flag, (const float2) make_float2(0,0), omega);
}

__global__ void wsinterpd(double2 * __restrict__ y, 
    const double2 * __restrict__ w, const double2 * __restrict__ x, 
    const double * __restrict__ tau, const size_t * sizes, const char * iflags, 
    const size_t * wstride,  const int flag, const double omega
    ) {
    wsinterpd_temp<double2, double, double>(y, w, x, tau, 
    sizes, iflags, wstride, flag, (const double2) make_double2(0,0), omega);
}

/* wsinterpd2 kernels */
#if (__CUDA_ARCH__ >= 530)
__global__ void wsinterpd2h(ushort2 * __restrict__ y, 
    const ushort2 * __restrict__ w, const ushort2 * __restrict__ x, 
    const unsigned short * __restrict__ tau1, const unsigned short * __restrict__ tau2, const size_t * sizes, const char * iflags, 
    const size_t * wstride, const int flag, const unsigned short omega) {
    wsinterpd2_temp<half2, half, float>((half2 *)y, (const half2 *)w, (const half2 *)x,
             (half *)tau1, (half *)tau2, sizes, iflags, wstride, flag, (const half2) make_half2(0,0), u2h(omega));
}
#endif

__global__ void wsinterpd2f(float2 * __restrict__ y, 
    const float2 * __restrict__ w, const float2 * __restrict__ x, 
    const float * __restrict__ tau1, const float * __restrict__ tau2, const size_t * sizes, const char * iflags, 
    const size_t * wstride, const int flag, const float omega
    ) {
    wsinterpd2_temp<float2, float, float>(y, w, x, tau1, tau2, 
    sizes, iflags, wstride, flag, (const float2) make_float2(0,0), omega);
}

__global__ void wsinterpd2(double2 * __restrict__ y, 
    const double2 * __restrict__ w, const double2 * __restrict__ x, 
    const double * __restrict__ tau1, const double * __restrict__ tau2, const size_t * sizes, const char * iflags, 
    const size_t * wstride,  const int flag, const double omega
    ) {
    wsinterpd2_temp<double2, double, double>(y, w, x, tau1, tau2, 
    sizes, iflags, wstride, flag, (const double2) make_double2(0,0), omega);
}
