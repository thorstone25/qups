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

// Modified for use as a stand-alone ptx
// data is (T x N)
// sample times are (I x N)
// # ifndef T
// __constant__ size_t T;
// # endif
// # ifndef I
// __constant__ size_t I;
// # endif
// # ifndef N
// __constant__ size_t N;
// # endif
// # ifndef M
// __constant__ size_t M;
// # endif


/// @brief Device function for nearest-neighbor interpolation
__device__ float2 nearestf(const float2 * x, float tau, float2 no_v) {
    const int ti = (int) roundf(tau); // round to nearest integer
    return (0 <= ti && ti < T) ? x[ti] : no_v;
}

/// @brief Device function for linear interpolation
__device__ float2 linearf(const float2 * x, float tau, float2 no_v) {
    float tf;
    
    // fractional and integer part
    tau = modf(tau, &tf); 
    const int ti = (int) tf; 
                
    // if in bounds, linearly interpolate by ratio tau at time-index ti[+1]
    return (0 <= ti && ti + 1 < T) ? lerp(x[ti + 1], x[ti], tau) : no_v;
}

/// @brief Device function for cubic Hermite interpolation
__device__ float2 cubicf(const float2 * x, float tau, float2 no_v) {
  float tf;
  float u = modff(tau, &tf);  // u is the fractional part, xf the integer part
  const int ti = (int) tf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < T))
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
  
  if (!(0 <= (ti - 1) && (ti + 2) < T))
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
    return (0 <= ti && ti < T) ? x[ti] : no_v;
}

/// @brief Device function for linear interpolation
__device__ double2 linear(const double2 * x, double tau, double2 no_v) {
    double tf;
    
    // fractional and integer part
    tau = modf(tau, &tf); 
    const int ti = (int) tf; 
                
    // if in bounds, linearly interpolate by ratio tau at time-index ti[+1]
    return (0 <= ti && ti + 1 < T) ? lerp(x[ti + 1], x[ti], tau) : no_v;
}

/// @brief Device function for cubic Hermite interpolation
__device__ double2 cubic(const double2 * x, double tau, double2 no_v) {
  double tf;
  double u = modf(tau, &tf);  // u is the fractional part, xf the integer part
  const int ti = (int) tf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < T))
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
  
  if (!(0 <= (ti - 1) && (ti + 2) < T))
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

    // get sampleing index
    const size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    float2 no_v = make_float2(0.0f, 0.0f);

    // if valid sample, for each tx/rx
    if(i < I && n < N){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){ // output
            y[i + n*I + m*I*N] = samplef(&x[n*T], tau[i + n*I + m*I*N], flag, no_v);
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

    // get sampleing index
    const size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    const size_t n = threadIdx.y + blockIdx.y * blockDim.y;
    double2 no_v = make_double2(0.0, 0.0);

    // if valid sample, for each tx/rx
    if(i < I && n < N){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){ // output
            y[i + n*I + m*I*N] = sample(&x[n*T], tau[i + n*I + m*I*N], flag, no_v);
        }
    }
}


