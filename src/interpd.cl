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

// constant defines for debugging / if not specified at compile time
// default precision
# ifndef QUPS_INTERPD_PRECISION
# define QUPS_INTERPD_PRECISION 32
# endif

// should be JIT compiles with (T2, U, V) set to the complex-data type, real-data type, and time type
#   if QUPS_INTERPD_PRECISION == 32
typedef float2  T2;
typedef float   U ;
typedef float   V ;
# define PI_VAL M_PI_F
# elif QUPS_INTERPD_PRECISION == 64 
#pragma OPENCL EXTENSION cl_khr_fp64 : enable // must enable double precision
typedef double2 T2;
typedef double  U ;
typedef double  V ;
# define PI_VAL M_PI
# elif QUPS_INTERPD_PRECISION == 16
#pragma OPENCL EXTENSION cl_khr_fp16 : enable // must enable half precision
typedef half2   T2;
typedef half    U ;
typedef half    V ;
# define PI_VAL M_PI_H
# endif

// complex number support (single/double)
// # include "clcomplex.h" // cmul

// default to avoid compiler/linter errors
# ifndef QUPS_INTERPD_FLAG
# define QUPS_INTERPD_FLAG 2
# endif
# ifndef QUPS_INTERPD_NO_V
# define QUPS_INTERPD_NO_V 0.f
# endif
# ifndef QUPS_INTERPD_OMEGA
# define QUPS_INTERPD_OMEGA 0
# endif

// DEBUG: to avoid linter errrors
// /*
# ifndef QUPS_T
# define QUPS_T 1024
# define QUPS_N 1
# define QUPS_M 1
# define QUPS_F 1
# define QUPS_I 1
# define QUPS_S 1024
# endif
// */

// complex multiplication 
T2 cmul(const T2 a, const T2 b){
    return (T2)(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

/// @brief Device function for nearest-neighbor interpolation
// template<typename U, typename T2>
T2 nearest(global const T2 * x, U tau, T2 no_v) {
    const int ti = (int) round(tau); // round to nearest integer
    return (0 <= ti && ti < QUPS_T) ? x[ti] : no_v;
}

T2 lerp(const T2 a, const T2 b, const U t) {
    return a + t*(b-a); 
}

/// @brief Device function for linear interpolation
// template<typename U, typename T2>
T2 linear(global const T2 * x, U tau, T2 no_v) {
    U tif; // integer part
    
    // fractional and integer part
    const U tf = modf(tau, &tif); 
    const int ti = (int) tif;
                
    // if in bounds, linearly interpolate by ratio tau at time-index ti[+1]
    return (0 <= ti && ti + 1 < QUPS_T) ? lerp(x[ti], x[ti + 1], tf) : no_v;
}

/// @brief Device function for cubic Hermite interpolation
// template<typename U, typename T2>
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
// template<typename T>
// float
U lanczos_helper(U u, int a) {
  return (u == 0.f ? 1.f : (2.f * sinpi(u)*sinpi(u/a) / (PI_VAL*PI_VAL * u*u)));
}


/// @brief Device function for Lanczos 3-lobe interpolation
// template<typename U, typename T2>
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

// template <typename U, typename T2>
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

// template<typename T2, typename U, typename V> // channel data type, time data type, time-sampling type
kernel void interpd(global T2 * y, global const T2 * x, global const U * tau, const int flag) {
    
    // get sampling index
    const size_t tid = get_global_id(0);
    const size_t n   = get_global_id(1);
    // const size_t m = threadIdx.z + blockIdx.z * blockDim.z;
    
    // rename for readability
    const size_t I = QUPS_I, M = QUPS_S, N = QUPS_N, T = QUPS_T, F = QUPS_F;
    const T2 no_v  = QUPS_INTERPD_NO_V;

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

// added from internet ------------------------------
/*
void __attribute__((always_inline)) atomic_add_f(volatile global float* addr, const float val) {
    union {
        uint  u32;
        float f32;
    } next, expected, current;
    current.f32 = *addr;
    do {
        next.f32 = (expected.f32=current.f32)+val; // ...*val for atomic_mul_f()
        current.u32 = atomic_cmpxchg((volatile global uint*)addr, expected.u32, next.u32);
    } while(current.u32!=expected.u32);
}

#ifdef cl_khr_int64_base_atomics
// #prgma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
void __attribute__((always_inline)) atomic_add_d(volatile global double* addr, const double val) {
    union {
        ulong  u64;
        double f64;
    } next, expected, current;
    current.f64 = *addr;
    do {
        next.f64 = (expected.f64=current.f64)+val; // ...*val for atomic_mul_d()
        current.u64 = atom_cmpxchg((volatile global ulong*)addr, expected.u64, next.u64);
    } while(current.u64!=expected.u64);
}
#endif
*/
// --------------------------------------------------------------------

// atomicAdd not natively supported for floating point types
// for double
#ifdef cl_khr_fp64
double atomicAdd(volatile global double* address, double val)
{
    // view data as a double (.f) or long int (.i)
    union {double f; long i;} old, curr, nxt; 
    old.f = *address; // previous value
    
    do {
        curr.i = old.i; // current value
        nxt.f = (curr.f + val); // desired next value
        old.i = atomic_cmpxchg((volatile global long *) address, curr.i, nxt.i); // attempt to swap

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (curr.i != old.i);

    return old.f;
}
#endif

// for float
float atomicAddf(volatile global float* address, float val)
{
    // view data as a float (.f) or int (.i)
    union {float f; int i;} old, curr, nxt; 
    old.f = *address; // previous value
    
    do {
        curr.i = old.i; // current value
        nxt.f = (curr.f + val); // desired next value
        old.i = atomic_cmpxchg((volatile global int *)address, curr.i, nxt.i); // attempt to swap

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (curr.i != old.i);

    return old.f;
}

// for half2
#ifdef cl_khr_fp16
half2 atomicAddh(volatile global half2* address, half2 val)
{
    // view data as a half2 (.f) or int (.i)
    union {half2 f; int i;} old, curr, nxt; 
    old.f = *address; // previous value
    
    do {
        curr.i = old.i; // current value
        nxt.f = (curr.f + val); // desired next value
        old.i = atomic_cmpxchg((volatile global int *) address, curr.i, nxt.i); // attempt to swap

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (curr.i != old.i);

    return old.f;
}
#endif


// add using the atomic functions for each data type
#ifdef cl_khr_fp16
inline void atomicAddStoreh(volatile global half2 * y, const half2 val){
    atomicAddh(y, val);
}
#endif
inline void atomicAddStoref(volatile global float2 * y, const float2 val){
    atomicAddf(((volatile global float*) y) + 0, val.x);
    atomicAddf(((volatile global float*) y) + 1, val.y);
}
#ifdef cl_khr_fp64
inline void atomicAddStore(volatile global double2 * y, const double2 val){
    atomicAdd(((volatile global double*) y) + 0, val.x);
    atomicAdd(((volatile global double*) y) + 1, val.y);
}
#endif


size_t global_offset(const size_t * dind, const global size_t * sizes, const global uchar * iflags){
// global index
    // init
    size_t dsz[3] = {1,1,1}; // {I,N,F} index cumulative sizes
    size_t str, j = 0; // stride, output index

    // find offset
    # pragma unroll
    for(char i = 0; i < 3; ++i){ // each label
        str = 1; // reset stride
        # pragma unroll
        for(size_t s = 0; s < QUPS_S; ++s){ // for each data dimension
            if(i == iflags[s]){ // matching label
                const size_t k = (dind[i] / dsz[i]) % sizes[s]; // get sub-index
                dsz[i] *= sizes[s]; // increase size for this label
                j += str * k; // add offset
            }
            str *= sizes[s]; // increase indexing stride
        }
    }

    return j;
}

// template<typename T2, typename U, typename V> // channel data type, time data type, time-sampling type
kernel void wsinterpd(volatile global T2 * y, 
    const global T2 * w, const global T2 * x, 
    const global U  * tau, const global ulong * sizes, 
    const global uchar * iflags, const global ulong * dstride){
    // , const int flag, const T2 no_v, const U omega){
    
    // alias constant inputs
    const int flag = QUPS_INTERPD_FLAG;
    const T2 no_v  = QUPS_INTERPD_NO_V;
    const U omega  = QUPS_INTERPD_OMEGA;

    // get sampling index
    const ulong i = get_global_id(0);
    const ulong n = get_global_id(1);
    
    // rename for readability
    const ulong I = QUPS_I, N = QUPS_N, F = QUPS_F, S = QUPS_S, T = QUPS_T; // , M = QUPS_M
    ulong u,v,k,l,sz; // weighting / output indexing

    // if valid sample, per each i,n,f
    if(i < I && n < N){
        # pragma unroll
        for(ulong f = 0; f < F; ++f){ // for f
            // global index
            const ulong dind[3] = {i,n,f};
            const ulong j = global_offset(dind, sizes, iflags);

            // get weight vector and output indices
            k = 0, l = 0; u = 0; v = 0;
            # pragma unroll
            for(ulong s = 0; s < S; ++s){ // for each dimension s
                // calculate the indexing stride for dimension s i.e. 
                // size of all prior dimensions
                sz = 1;
                # pragma unroll
                for(ulong sp = 0; sp < s; ++sp)
                    sz *= sizes[sp]; 

                const ulong js = ((j / sz) % sizes[s]); // index for this dimension
                k += js * dstride[0 + 4*s]; // add pitched index for this dim (weights)
                l += js * dstride[1 + 4*s]; // add pitched index for this dim (outputs)
                u += js * dstride[2 + 4*s]; // add pitched index for this dim (time)
                v += js * dstride[3 + 4*s]; // add pitched index for this dim (samples)
            }
            
            const T2 a = (omega != 0) ? (T2)(cos(omega * tau[u]), sin(omega * tau[u])) : (T2)(1.f, 0.f); // modulation phasor
            const T2 val = cmul(a, cmul(w[k],sample(&x[(v)*T], (V)tau[u], flag, no_v))); // weighted sample
            // y[l] += (T2)(1.0f, 1.0f);

#   if QUPS_INTERPD_PRECISION == 64
            atomicAddStore(y + l, val); // store
# elif QUPS_INTERPD_PRECISION == 32
            atomicAddStoref(y + l, val); // store
# elif QUPS_INTERPD_PRECISION == 16
            atomicAddStoreh(y + l, val); // store
# endif

        }
    }
}

// template<typename T2, typename U, typename V> // channel data type, time data type, time-sampling type
kernel void wsinterpd2(volatile global T2 * y, 
    const global T2 *  w, const global T2 * x, 
    const global U * tau1, const global U * tau2, const global ulong * sizes, 
    const global uchar * iflags, const global ulong * dstride){
    // , const int flag, const T2 no_v, const U omega) {

    // alias constant inputs
    const int flag = QUPS_INTERPD_FLAG;
    const T2 no_v  = QUPS_INTERPD_NO_V;
    const U omega  = QUPS_INTERPD_OMEGA;

    // get sampling index
    const size_t i = get_global_id(0);
    const size_t n = get_global_id(1);
    
    // rename for readability
    const size_t I = QUPS_I,N = QUPS_N, F = QUPS_F, S = QUPS_S, T = QUPS_T; // , M = QUPS_M
    size_t r,u,v,k,l,sz; // weighting / output indexing

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
            const T2 a = (omega != 0) ? (T2)(cos(omega * t), sin(omega * t)) : (T2)(1.f, 0.f); // modulation phasor
            const T2 val = cmul(a, cmul(w[k], sample(&x[(v)*T], (V)t, flag, no_v))); // weighted sample
#   if QUPS_INTERPD_PRECISION == 64
            atomicAddStore(&y[l], val); // store
# elif QUPS_INTERPD_PRECISION == 32
            atomicAddStoref(&y[l], val); // store
# elif QUPS_INTERPD_PRECISION == 16
            atomicAddStoreh(&y[l], val); // store
# endif
        }
    }
}

/*

// interd kernels
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


// wsinterpd kernels
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

// wsinterpd2 kernels
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
*/

