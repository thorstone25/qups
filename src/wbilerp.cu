# include "helper_math.h" // vector math

// int2 times float2 vectorized
inline __host__ __device__ float2 operator*(int2 a, float2 b)
{
    return make_float2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ double2 operator*(int2 a, double2 b)
{
    return make_double2(a.x * b.x, a.y * b.y);
}
#if (__CUDA_ARCH__ >= 530)
// /*
inline __host__ __device__ half2 operator*(int2 a, half2 b)
{
    return make_half2((half)a.x * b.x, (half)a.y * b.y);
} 
// */
#endif

template <typename U, typename U2>
__device__ void wbilerp_temp(
    int2 * ixy, U * cxy, const U2 * pall, 
    const U * xg, const U * yg, 
    const size_t X, const size_t Y, const size_t I,
    const U2 * pml, const U2 * pmu
    ){

    const size_t i = threadIdx.x + blockIdx.x * blockDim.x;
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
            const U2 dp{xg[ix+1] - xg[ix], yg[iy+1] - yg[iy]};
        
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
                    const U2 q{xg[ixloc], yg[iyloc]};
        
                    ////// inline the integral //////
                    // s_d == -1 if in positive quadrant of dimension d
                    const int2 s = make_int2(
                        ( (p1.x <= q.x) && (p2.x <= q.x) ) ? 1 : -1,
                        ( (p1.y <= q.y) && (p2.y <= q.y) ) ? 1 : -1
                    );
        
                    // coefficients for the integral
                    // TODO: use copysign instead of multiplying by +/- 1
                    // TODO: write out equations to avoid accidental complex overload
                    const U2 ONE{1.0, 1.0};
                    const U2 a = ONE + s * (p1 - q ) / dp; // (non-complex) vector division
                    const U2 b =       s * (p1 - p2) / dp; // (non-complex) vector division
        
                    // evaluation of the integral: subject to numerical precision issues
                    const U c = a.x*a.y - (a.x*b.y + a.y*b.x) / (U)2.0 + (b.x*b.y) / (U)3.0;
        
                    // apply line integral formula w/r to quadrant
                    // TODO: redo write indexing to vectorize writes
                    ixy[i*4*N + 4*n + 2*iyy + ixx] = make_int2(ixloc, iyloc); // output indices
                    cxy[i*4*N + 4*n + 2*iyy + ixx] = l * c; // length in this grid scaled by grid size
                }
            }
        }
    }
}


__global__ void wbilerpf(
    int2 * ixy, float * cxy, const float2 * pall, 
    const float * xg, const float * yg, 
    const size_t X, const size_t Y, const size_t I,
    const float2 * pml, const float2 * pmu
    ){
    wbilerp_temp<float, float2>(ixy, cxy, pall, xg, yg, X, Y, I, pml, pmu);
}

__global__ void wbilerp(
    int2 * ixy, double * cxy, const double2 * pall, 
    const double * xg, const double * yg, 
    const size_t X, const size_t Y, const size_t I,
    const double2 * pml, const double2 * pmu
    ){
    wbilerp_temp<double, double2>(ixy, cxy, pall, xg, yg, X, Y, I, pml, pmu);
}

#if (__CUDA_ARCH__ >= 530)
__global__ void wbilerph(
    int2 * ixy, unsigned short * cxy, const ushort2 * pall, 
    const unsigned short * xg, const unsigned short * yg, 
    const size_t X, const size_t Y, const size_t I,
    const ushort2 * pml, const ushort2 * pmu
    ){
    wbilerp_temp<half, half2>(
        ixy, (half *) cxy, (half2 *) pall, 
        (half *) xg, (half *) yg, X, Y, I, 
        (half2 *) pml, (half2 *) pmu
    );
}
#endif
