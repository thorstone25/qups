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

// Modified for use as a stand-alone OpenCL file (Thurston Brevett <tbrevett@stanford.edu>)
// data is (T x N x F)
// sample times are (I x N x M)


// enable integer atomics(must be supported!)
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable 

// DEBUG: constant defines for debugging / if not specified at compile time
// default precision
// # ifndef QUPS_PRECISION
// # define QUPS_PRECISION 32
// # endif

// DEBUG: to avoid linter errrors
/*
# ifndef QUPS_T
# define QUPS_T 1024
# define QUPS_N 1
# define QUPS_M 1
# define QUPS_F 1
# define QUPS_I 1
# define QUPS_S 1024
# endif
*/

// should be JIT compiles with (T2, U, V) set to the complex-data type, real-data type, and time type
#   if QUPS_PRECISION == 32
typedef float2  T2;
typedef float   U ;
typedef float   V ;
# define PI_VAL M_PI_F
# elif QUPS_PRECISION == 64 
#pragma OPENCL EXTENSION cl_khr_fp64 : enable // must enable double precision
typedef double2 T2;
typedef double  U ;
typedef double  V ;
# define PI_VAL M_PI
# elif QUPS_PRECISION == 16
#pragma OPENCL EXTENSION cl_khr_fp16 : enable // must enable half precision
typedef half2   T2;
typedef half    U ;
typedef half    V ;
# define PI_VAL M_PI_H
# endif

// complex multiplication 
T2 cmul(const T2 a, const T2 b){
    return (T2)(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

/// @brief Device function for nearest-neighbor interpolation
T2 nearest(global const T2 * x, U tau, T2 no_v) {
    const int ti = (int) round(tau); // round to nearest integer
    return (0 <= ti && ti < QUPS_T) ? x[ti] : no_v;
}

T2 lerp(const T2 a, const T2 b, const U t) {
    return a + t*(b-a); 
}

/// @brief Device function for linear interpolation
T2 linear(global const T2 * x, U tau, T2 no_v) {
    U tif; // integer part
    
    // fractional and integer part
    const U tf = modf(tau, &tif); 
    const int ti = (int) tif;
                
    // if in bounds, linearly interpolate by ratio tau at time-index ti[+1]
    return (0 <= ti && ti + 1 < QUPS_T) ? lerp(x[ti], x[ti + 1], tf) : no_v;
}

/// @brief Device function for cubic Hermite interpolation
T2 cubic(global const T2 * x, U tau, T2 no_v) {
  U tf;
  const U u = modf(tau, &tf); // u is the fractional part, tf the integer part
  const int ti = (int) tf;
  
  if (!(0 <= (ti - 1) && (ti + 2) < QUPS_T))
      return no_v;

  T2 s0 = x[ti - 1];
  T2 s1 = x[ti + 0];
  T2 s2 = x[ti + 1];
  T2 s3 = x[ti + 2];

  // Cubic Hermite interpolation (increased precision using fused multiply-adds)
  // (Catmull-Rom)
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
U lanczos_helper(U u, int a) {
  return (u == 0.f ? 1.f : (2.f * sinpi(u)*sinpi(u/a) / (PI_VAL*PI_VAL * u*u)));
}


/// @brief Device function for Lanczos 3-lobe interpolation
T2 lanczos3(global const T2 * x, U tau, T2 no_v) {
  const int a = 2;  // a=2 for 3-lobe Lanczos resampling
  U tf;
  const U u = modf(tau, &tf); // u is the fractional part, tf the integer part
  const int ti = (int) tf;
  
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

T2 sample(global const T2 * x, U tau, int flag, const T2 no_v){
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

