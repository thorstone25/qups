# include "interpolators.cl"

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
