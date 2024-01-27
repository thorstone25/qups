// DEBUG: for linter purposes
/*
# ifndef QUPS_INTERPD_PRECISION
# define QUPS_INTERPD_PRECISION 32
# endif

# ifndef QUPS_COMPLEX
# define QUPS_COMPLEX 
# endif
*/

#   if QUPS_INTERPD_PRECISION == 32
    # ifdef QUPS_COMPLEX
typedef float2 T;
    # else
typedef float T;
    # endif
# elif QUPS_INTERPD_PRECISION == 16
#pragma OPENCL EXTENSION cl_khr_fp16 : enable // must enable half precision
    # ifdef QUPS_COMPLEX
typedef half2 T;
    # else
typedef half T;
    # endif
# elif QUPS_INTERPD_PRECISION == 64 
#pragma OPENCL EXTENSION cl_khr_fp64 : enable // must enable double precision
    # ifdef QUPS_COMPLEX
typedef double2 T;
    # else
typedef double T;
    # endif
# endif



// complex multiplication
# ifdef QUPS_COMPLEX
inline T cmul(const T a, const T b){ return (T)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);}
# else
inline T cmul(const T a, const T b){ return a * b;}
# endif

/*
* Compute the cross correlation of two sets of data. The data will be
* correlated in the first dimension. M >= N must be satisfied.
*
* Inputs:
*  x:         first signal  (M x S)
*  y:         second signal (N x S)
*  sizes:     sizing info (undocumented)
*
* Outputs:
*  z:         resulting cross correlation
*
*
*/

kernel void conv(global const T* x, global const T* y, global T* z, const global int * sizes){
    // get stride and len
    const int istr  = sizes[0];
    const int ilen  = sizes[1];
    const int jstr  = sizes[2];
    const int jlen  = sizes[3];
    const int lstr  = sizes[4];
    const int llen  = sizes[5];
    const int clen  = sizes[6];
    
    const int I = ilen * clen;
    const int J = jlen * clen;
    const int L = llen * clen;
    
    const int c = get_global_id(0); // column output index 
    const int l = get_global_id(1); // strided output index 
    const int s = get_global_id(2); // batch index 
    
    // accumulator
    T za = (T) 0.f;

    // if valid lag indices, multiply and accumulate in-place
    // TODO: apply work-groups and striding for faster accumulation
    if(l < L && c < clen)
    # pragma unroll
        for(int i = 0, j = (int)(QUPS_L0) - l, jr = jlen - 1 - j; i < ilen || jr >= 0; ++i, --jr) // forward and reverse iterate over both signals
            if(0 <= i && i < ilen && 0 <= jr && jr < jlen) // if both signals in bounds
                za += cmul(x[c + i*istr + s*I] , y[c + jr*jstr + s*J]); // accumulate the cross product

    // output result:
    if(l < L && c < clen)
        z[c + l*lstr + s*L] = za;
}

