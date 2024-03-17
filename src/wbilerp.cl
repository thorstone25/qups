// linter ...
# ifndef QUPS_PRECISION
# define QUPS_PRECISION 32
# endif

#   if QUPS_PRECISION == 32
typedef float  U;
typedef float2 U2;
# elif QUPS_PRECISION == 64 
#pragma OPENCL EXTENSION cl_khr_fp64 : enable // must enable double precision
typedef double  U;
typedef double2 U2;
# elif QUPS_PRECISION == 16
#pragma OPENCL EXTENSION cl_khr_fp16 : enable // must enable half precision
typedef half  U;
typedef half2 U2;
# endif

// convert integer to U type
inline U2 i2U(const int2 a){return (U2)(a.x, a.y);}

kernel void wbilerp(
    global int2 * ixy, global U * cxy, global const U2 * pall, 
    global const U * xg, global const U * yg, 
    const ulong X, const ulong Y, const ulong I,
    global const U2 * pml, global const U2 * pmu
    ){

    const ulong i = get_global_id(0);
    const int N = X + Y + 1; // total number of points in output, 
    int ix = 0; // grid x position (iterator)
    int iy = 0; // grid y position (iterator)

    if(i < I){
        // index for this batch of pall
        const int ip = i * (N + 1); // pall is length N+1

        // compute all integrals one segment at a time
        for(int n = 0; n < N; ++n){
            // get both points (x/y by pair)
            const U2 p1   = pall[ip + n  ];
            const U2 p2   = pall[ip + n+1];
        
            // length of the line segment
            // TODO: switch to hypot for numerical accuracy
            const U l = length(p2 - p1); // hypot(p2(1) - p1(1), p2(2) - p1(2));
        
            // if points are within the support of the grid / line
            const bool val = 
                  pml[i].x <= p1.x & p1.x <= pmu[i].x
                & pml[i].y <= p1.y & p1.y <= pmu[i].y
                & pml[i].x <= p2.x & p2.x <= pmu[i].x
                & pml[i].y <= p2.y & p2.y <= pmu[i].y;
        
            // skip integral if out of grid / line support, the segment length is zero, or if co-located in x
            if ((!val) || (l == (U)0) || (p2.x == p1.x)) continue;
        
            // if any point not in x/y-bounds, move to next x/y-grid interval
            while ((p1.x >= xg[ix+1]) && (p2.x >= xg[ix+1]) && (ix+2 < X)) ++ix;
            while ((p1.y >= yg[iy+1]) && (p2.y >= yg[iy+1]) && (iy+2 < Y)) ++iy;
            
            // get the line segment region's size
            const U2 dp = (U2)(xg[ix+1] - xg[ix], yg[iy+1] - yg[iy]);
        
            // get the indices for each of the four grid points affected by this
            // line segment
            # pragma unroll
            for(int iyy = 0; iyy < 2; ++iyy){
                # pragma unroll
                for(int ixx = 0; ixx < 2; ++ixx){
                    // ixloc = [ix, ix+1, ix, ix+1]; // x index
                    // iyloc = [iy, iy, iy+1, iy+1]; // y index
                    const int ixloc = ix + ixx; // x index
                    const int iyloc = iy + iyy; // y index
        
                    // grid points
                    // q = [xg(ixloc)'; yg(iyloc)'];
                    const U2 q = (U2)(xg[ixloc], yg[iyloc]);
        
                    ////// inline the integral //////
                    // s_d == -1 if in positive quadrant of dimension d
                    const int2 s = (int2)(
                        ( (p1.x <= q.x) && (p2.x <= q.x) ) ? 1 : -1,
                        ( (p1.y <= q.y) && (p2.y <= q.y) ) ? 1 : -1
                    );
        
                    // coefficients for the integral
                    // TODO: use copysign instead of multiplying by +/- 1
                    // TODO: write out equations to avoid accidental complex overload
                    const U2 ONE = (U2)(1.f, 1.f);
                    const U2 a = ONE + i2U(s) * (p1 - q ) / dp; // (non-complex) vector division
                    const U2 b =       i2U(s) * (p1 - p2) / dp; // (non-complex) vector division
        
                    // evaluation of the integral: subject to numerical precision issues
                    const U c = a.x*a.y - (a.x*b.y + a.y*b.x) / (U)2 + (b.x*b.y) / (U)3;
        
                    // apply line integral formula w/r to quadrant
                    // TODO: redo write indexing to vectorize writes
                    ixy[i*4*N + 4*n + 2*iyy + ixx] = (int2)(ixloc, iyloc); // output indices
                    cxy[i*4*N + 4*n + 2*iyy + ixx] = l * c; // length in this grid scaled by grid size
                }
            }
        }
    }
}
