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

/// @brief Device function for nearest-neighbor interpolation
__device__ half2 nearesth(const ushort2 * x, half tau, half2 no_v) {
    const int ti = (int) roundf(tau); // round to nearest integer
    return (0 <= ti && ti < QUPS_T) ? u2h(x[ti]) : no_v;
}

/// @brief Device function for linear interpolation
__device__ half2 linearh(const ushort2 * x, half tau, half2 no_v) {
    float tf;
    
    // fractional and integer part
    tau = modf(tau, &tf); 
    const int ti = (int) tf; 
                
    // if in bounds, linearly interpolate by ratio tau at time-index ti[+1]
    return (0 <= ti && ti + 1 < QUPS_T) ? lerp(u2h(x[ti]), u2h(x[ti + 1]), tau) : no_v;
}

/// @brief Device function for cubic Hermite interpolation
__device__ half2 cubich(const ushort2 * x, half tau, half2 no_v) {
  float tf;
  float u = modff(tau, &tf);  // u is the fractional part, xf the integer part
  const int ti = (int) tf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;

  half2 s0 = u2h(x[ti - 1]);
  half2 s1 = u2h(x[ti + 0]);
  half2 s2 = u2h(x[ti + 1]);
  half2 s3 = u2h(x[ti + 2]);

  // Cubic Hermite interpolation (increased precision using fused multiply-adds)
  half a0 = 0 + u * (-1 + u * (+2 * u - 1));
  half a1 = 2 + u * (+0 + u * (-5 * u + 3));
  half a2 = 0 + u * (+1 + u * (+4 * u - 3));
  half a3 = 0 + u * (+0 + u * (-1 * u + 1));
  // // Cubic Hermite interpolation (naive, less precise implementation)
  // float a0 = -1 * u * u * u + 2 * u * u - 1 * u + 0;
  // float a1 = +3 * u * u * u - 5 * u * u + 0 * u + 2;
  // float a2 = -3 * u * u * u + 4 * u * u + 1 * u + 0;
  // float a3 = +1 * u * u * u - 1 * u * u + 0 * u + 0;
  return (s0 * a0 + s1 * a1 + s2 * a2 + s3 * a3) * 0.5f;
}

/// @brief Inline helper code for Lanczos 3-lobe interpolation
inline __device__ half lanczos_helper(half u, int a) {
  const half PI = (atanf(1) * 4); // TODO: isn't there a built-in PI value?
  half ret;
  if (u == (half) 0.)
      ret = (half) 1.;
  else {
      const half spu  = (half) sinpif(__half2float(u)    );
      const half spua = (half) sinpif(__half2float(u) / a);
      ret = (half) 2. * spu * spua / (PI * PI * u * u);
  }
  return ret;
  // return (u == (half)0.f) ? (half) 1.f : (half) 2.f * sinpif(__half2float(u)) * sinpif(__half2float(u) / a) / (PI * PI * u * u);
}

/// @brief Device function for Lanczos 3-lobe interpolation
__device__ half2 lanczos3h(const ushort2 * x, half tau, half2 no_v) {
  constexpr int a = 2;  // a=2 for 3-lobe Lanczos resampling
  float xf;
  float u = modff(tau, &xf);  // u is the fractional part, xf the integer part
  const int ti = (int) xf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;
  
  half2 s0 = u2h(x[ti - 1]);
  half2 s1 = u2h(x[ti + 0]);
  half2 s2 = u2h(x[ti + 1]);
  half2 s3 = u2h(x[ti + 2]);
  half a0 = lanczos_helper(u + 1, a);
  half a1 = lanczos_helper(u + 0, a);
  half a2 = lanczos_helper(u - 1, a);
  half a3 = lanczos_helper(u - 2, a);
  return s0 * a0 + s1 * a1 + s2 * a2 + s3 * a3;
}

/// @brief Device function for nearest-neighbor interpolation
__device__ float2 nearestf(const float2 * x, float tau, float2 no_v) {
    const int ti = (int) roundf(tau); // round to nearest integer
    return (0 <= ti && ti < QUPS_T) ? x[ti] : no_v;
}

/// @brief Device function for linear interpolation
__device__ float2 linearf(const float2 * x, float tau, float2 no_v) {
    float tf;
    
    // fractional and integer part
    tau = modf(tau, &tf); 
    const int ti = (int) tf; 
                
    // if in bounds, linearly interpolate by ratio tau at time-index ti[+1]
    return (0 <= ti && ti + 1 < QUPS_T) ? lerp(x[ti], x[ti + 1], tau) : no_v;
}

/// @brief Device function for cubic Hermite interpolation
__device__ float2 cubicf(const float2 * x, float tau, float2 no_v) {
  float tf;
  float u = modff(tau, &tf);  // u is the fractional part, xf the integer part
  const int ti = (int) tf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;

  float2 s0 = x[ti - 1];
  float2 s1 = x[ti + 0];
  float2 s2 = x[ti + 1];
  float2 s3 = x[ti + 2];

  // Cubic Hermite interpolation (increased precision using fused multiply-adds)
  float a0 = 0 + u * (-1 + u * (+2 * u - 1));
  float a1 = 2 + u * (+0 + u * (-5 * u + 3));
  float a2 = 0 + u * (+1 + u * (+4 * u - 3));
  float a3 = 0 + u * (+0 + u * (-1 * u + 1));
  // // Cubic Hermite interpolation (naive, less precise implementation)
  // float a0 = -1 * u * u * u + 2 * u * u - 1 * u + 0;
  // float a1 = +3 * u * u * u - 5 * u * u + 0 * u + 2;
  // float a2 = -3 * u * u * u + 4 * u * u + 1 * u + 0;
  // float a3 = +1 * u * u * u - 1 * u * u + 0 * u + 0;
  return (s0 * a0 + s1 * a1 + s2 * a2 + s3 * a3) * 0.5f;
}

/// @brief Inline helper code for Lanczos 3-lobe interpolation
inline __device__ float lanczos_helper(float u, int a) {
  const float PI = atanf(1) * 4; // TODO: isn't there a built-in PI value?
  return (u == 0.f) ? 1.f : 2.f * sinpif(u) * sinpif(u / a) / (PI * PI * u * u);
}

/// @brief Device function for Lanczos 3-lobe interpolation
__device__ float2 lanczos3f(const float2 * x, float tau, float2 no_v) {
  constexpr int a = 2;  // a=2 for 3-lobe Lanczos resampling
  float xf;
  float u = modff(tau, &xf);  // u is the fractional part, xf the integer part
  const int ti = (int) xf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;
  
  float2 s0 = x[ti - 1];
  float2 s1 = x[ti + 0];
  float2 s2 = x[ti + 1];
  float2 s3 = x[ti + 2];
  float a0 = lanczos_helper(u + 1, a);
  float a1 = lanczos_helper(u + 0, a);
  float a2 = lanczos_helper(u - 1, a);
  float a3 = lanczos_helper(u - 2, a);
  return s0 * a0 + s1 * a1 + s2 * a2 + s3 * a3;
}

/// @brief Device function for nearest-neighbor interpolation
__device__ double2 nearest(const double2 * x, double tau, double2 no_v) {
    const int ti = (int) roundf(tau); // round to nearest integer
    return (0 <= ti && ti < QUPS_T) ? x[ti] : no_v;
}

/// @brief Device function for linear interpolation
__device__ double2 linear(const double2 * x, double tau, double2 no_v) {
    double tf;
    
    // fractional and integer part
    tau = modf(tau, &tf); 
    const int ti = (int) tf; 
                
    // if in bounds, linearly interpolate by ratio tau at time-index ti[+1]
    return (0 <= ti && ti + 1 < QUPS_T) ? lerp(x[ti], x[ti + 1], tau) : no_v;
}

/// @brief Device function for cubic Hermite interpolation
__device__ double2 cubic(const double2 * x, double tau, double2 no_v) {
  double tf;
  double u = modf(tau, &tf);  // u is the fractional part, xf the integer part
  const int ti = (int) tf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;

  double2 s0 = x[ti - 1];
  double2 s1 = x[ti + 0];
  double2 s2 = x[ti + 1];
  double2 s3 = x[ti + 2];

  // Cubic Hermite interpolation (increased precision using fused multiply-adds)
  double a0 = 0 + u * (-1 + u * (+2 * u - 1));
  double a1 = 2 + u * (+0 + u * (-5 * u + 3));
  double a2 = 0 + u * (+1 + u * (+4 * u - 3));
  double a3 = 0 + u * (+0 + u * (-1 * u + 1));
  // // Cubic Hermite interpolation (naive, less precise implementation)
  // double a0 = -1 * u * u * u + 2 * u * u - 1 * u + 0;
  // double a1 = +3 * u * u * u - 5 * u * u + 0 * u + 2;
  // double a2 = -3 * u * u * u + 4 * u * u + 1 * u + 0;
  // double a3 = +1 * u * u * u - 1 * u * u + 0 * u + 0;
  return (s0 * a0 + s1 * a1 + s2 * a2 + s3 * a3) * 0.5;
}

/// @brief Inline helper code for Lanczos 3-lobe interpolation
inline __device__ double lanczos_helper(double u, int a) {
  const double PI = atanf(1) * 4; // TODO: isn't there a built-in PI value?
  return (u == 0.) ? 1. : 2. * sinpif(u) * sinpif(u / a) / (PI * PI * u * u);
}

/// @brief Device function for Lanczos 3-lobe interpolation
__device__ double2 lanczos3(const double2 * x, double tau, double2 no_v) {
  constexpr int a = 2;  // a=2 for 3-lobe Lanczos resampling
  double xf;
  double u = modf(tau, &xf);  // u is the fractional part, xf the integer part
  const int ti = (int) xf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;
  
  double2 s0 = x[ti - 1];
  double2 s1 = x[ti + 0];
  double2 s2 = x[ti + 1];
  double2 s3 = x[ti + 2];
  double a0 = lanczos_helper(u + 1, a);
  double a1 = lanczos_helper(u + 0, a);
  double a2 = lanczos_helper(u - 1, a);
  double a3 = lanczos_helper(u - 2, a);
  return s0 * a0 + s1 * a1 + s2 * a2 + s3 * a3;
}

inline __device__ half2 sampleh(const ushort2 * x, half tau, int flag, const half2 no_v){
    half2 val;
// sample according to the flag
         if (flag == 0)
        val = nearesth  (x, tau, no_v);
    else if (flag == 1)
        val = linearh   (x, tau, no_v);
    else if (flag == 2)
        val = cubich    (x, tau, no_v);
    else if (flag == 3)
        val = lanczos3h (x, tau, no_v);
    else if (flag == 4)
        val = linearh   (x, tau, no_v);
    else 
        val = no_v; 
    
    return val;
}


__global__ void interpdh(ushort2 * __restrict__ y, 
    const ushort2 * __restrict__ x, const unsigned short * __restrict__ tau, const int flag) {

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
            y[i + n*I + m*N*I + f*M*N*I] = h2u(sampleh(&x[n*T + f*N*T], u2h(tau[i + n*I + m*I*N]), flag, no_v));
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

            const half2 val = u2h(w[k]) * sampleh(&x[n*T + f*N*T], u2h(tau[i + n*I + m*I*N]), flag, no_v); // weighted sample
            atomicAdd(&yh[l], val); // store
        }
    }
}




inline __device__ float2 samplef(const float2 * x, float tau, int flag, const float2 no_v){
    float2 val;
// sample according to the flag
         if (flag == 0)
        val = nearestf  (x, tau, no_v);
    else if (flag == 1)
        val = linearf   (x, tau, no_v);
    else if (flag == 2)
        val = cubicf    (x, tau, no_v);
    else if (flag == 3)
        val = lanczos3f (x, tau, no_v);
    else if (flag == 4)
        val = linearf   (x, tau, no_v);
    else 
        val = no_v; 
    
    return val;
}


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
            y[i + n*I + m*N*I + f*M*N*I] = samplef(&x[n*T + f*N*T], tau[i + n*I + m*I*N], flag, no_v);
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

            const float2 val = w[k] * samplef(&x[n*T + f*N*T], tau[i + n*I + m*I*N], flag, no_v); // weighted sample
            atomicAdd(&y[l].x, val.x); // store
            atomicAdd(&y[l].y, val.y); // store
        }
    }
}



inline __device__ double2 sample(const double2 * x, double tau, int flag, const double2 no_v){
    double2 val;
    // sample according to the flag
         if (flag == 0)
        val = nearest  (x, tau, no_v);
    else if (flag == 1)
        val = linear   (x, tau, no_v);
    else if (flag == 2)
        val = cubic    (x, tau, no_v);
    else if (flag == 3)
        val = lanczos3 (x, tau, no_v);
    else 
        val = no_v; 
    
    return val;
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

