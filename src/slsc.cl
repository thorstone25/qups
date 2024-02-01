// define precision
#if SLSC_PRECISION_SINGLE
typedef float  T ;
typedef float2 T2;

# elif SLSC_PRECISION_DOUBLE
#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
typedef double  T ;
typedef double2 T2;

# elif SLSC_PRECISION_HALF
#ifdef cl_khr_fp16
#pragma OPENCL EXTENSION cl_khr_fp16 : enable
#endif
typedef half  T ;
typedef half2 T2;

#endif

// complex maths
// conjugate
inline T2 conj(T2 a){
    return (T2)(a.x, -a.y);
}
// multiplication
inline T2 cmul(T2 a, T2 b){
    return (T2)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}
// power
inline T cpow(T2 a){
    return a.x*a.x + a.y*a.y;
}

kernel void slsc(
  global T * yi, global T * xi, global uint * lags
//  , uint L, uint I, uint N, uint J
){
  
    // indexing: input (I x N x J), output (I x 1 x J)
    uint i = get_global_id(0); // vector-index (stride is "I")
    uint j = get_global_id(1); // slice-index  (stride is "I*N")
    
    const uint istr = 1; // element stride
    const uint nstr = I; // channel stride
    const uint jstr = I*N; // batch stride

    // variable casting
    global T2 * x = (global T2 *) xi;
    global T2 * y = (global T2 *) yi;

    // accumulators
    T2 z = (T2)(0,0);
    T  a = 0, b = 0;

    // compute ensemble
    for(char l = 0; l < (char)L; ++l) // lags
        for(uint n = 0; n < (uint)N; ++n) // channels / cross channel
            for(char s = 0; s < 2; ++s) // positive / negative lags
            { 
                uint m = n + (s ? -lags[l] : lags[l]); // cross-channel
                if(m < (uint)N){
                    T2 u  = x[i*istr + n*nstr + j*jstr];
                    T2 v  = x[i*istr + m*nstr + j*jstr];
                    if(SLSC_FLAG == 0){ // ensemble
                        a += cpow(u);
                        b += cpow(v);
                        z += cmul(u, conj(v));
                    } else if(SLSC_FLAG == 1) { // average
                        const T w = (((T)N) - lags[l]) * (T)L * 2; // averaging weights
                        const T un = rsqrt(cpow(u)); // power normalization
                        const T vn = rsqrt(cpow(v)); // power normalization
                        z += cmul(u, conj(v)) * un * vn / w;
                    }
                }
            }
        

    if (SLSC_FLAG == 0) z *= rsqrt(a * b); // final power correction for ensemble only
    y[i + I*j] = z;

}

