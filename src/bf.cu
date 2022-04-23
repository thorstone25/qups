# include "helper_math.h" // vector math


// data sizes: I -> pixels, M -> transmitters, N -> receivers, 
// T -> time samples, S -> signal
        
# ifndef T
__constant__ size_t T;
# endif
# ifndef M
__constant__ size_t M;
# endif
# ifndef N
__constant__ size_t N;
# endif
# ifndef I
__constant__ size_t I;
# endif
# ifndef S
__constant__ size_t S;
# endif

// integral params
# ifndef W
__constant__ size_t W; // temporal integration width (one-sided)
# endif
# ifndef D
__constant__ size_t D; // number of points in the spatial integration
# endif

// Beamforming transmit mode: focal point or plane-wave
# ifndef VS
__constant__ bool VS; // whither virtual source mode
# endif


/* Creates the positions for a linear array aperture given it's description
*
* All positions are in projetive coordinates.
*
* Must be run with a kernel size equal to the number of elements in the 
* array.
*
* Inputs:
*   Pn:     Vector of positions
*   Pn0:    Initial position
*   dPn:    Interelement difference
*
*
*/ 
        
__global__ void pos_step_rng_lenf(float3 * Pn, const float3 Pn0, const float3 dPn){
    const uint idx = threadIdx.x + blockIdx.x * blockDim.x;
    Pn[idx] = Pn0 + idx * dPn;
}

__global__ void pos_step_rng_len(double3 * Pn, const double3 Pn0, const double3 dPn){
    const uint idx = threadIdx.x + blockIdx.x * blockDim.x;
    Pn[idx] = Pn0 + idx * dPn;
}

/*
* Delay and sum the given data at the given pixels
*
* Given a set of pixels, (virtual or plane wave) transmitter locations, 
* receiver locations, as well as a datacube equipped with a time, 
* transmitter and receiver axis, perform simple delay-and-sum beamforming. 
* The data is linearly interpolated at the sample time. An image is 
* generated for each receiver element. Summation across the transmitters 
* and receivers is implicit.
*
* All positions are in vector coordinates. 
* 
* If the virtual transmitter 
* normal has a fourth component that is 0, this indicates that the 
* transmission should be treated as a plane wave transmission instead of a 
* virtual source transmission. 
* 
* The value of t = 0 must be the time when the peak of the wavefront 
* reaches the virtual source location. Because this time must be the same 
* for all transmits, the datacube must be stitched together in such a way 
* that for all transmits, the same time axis is used
*
* Inputs:
*  y:           complex pixel values per channel (I)
*  Pi:          pixel positions (3 x I)
*  Pr:          receiver positions (3 x N)
*  Pv:          (virtual) transmitter positions (3 x M)
*  Nv:          (virtual) transmitter normal (3 x M)
*  x:           datacube of complex sample values (T x M x N)
*  t0:          initial time for the data
*  fs:          sampling frequency of the data
*  cinv:        inverse of the speed of sound used for beamforming
* 
* I -> pixels, M -> transmitters, N -> receivers, T -> time samples
*
*/

__global__ void DASf(float2 * __restrict__ y, 
    const float * __restrict__ Pi, const float * __restrict__ Pr, 
    const float * __restrict__ Pv, const float * __restrict__ Nv, 
	const float2 * __restrict__ x, 
	const float t0, const float fs, const float cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 4 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M
    
    // temp vars
    const float2 zero_v = make_float2(0.0f);
    float rf, dv, dr;
    float2 val, pix = zero_v;
    float3 rv;
    union {float f; int i;} t;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // fractional and integer part
                rf = modf(rf, &t.f);
                t.i = (int)(t.f);
                
                // if in bounds, linearly interpolate by ratio rf at time-index t.i[+1]
                val = (0 <= t.i && t.i + 1 < T) ?            
                   lerp(x[t.i + 1 + (n + m * N) * T],
                        x[t.i +     (n + m * N) * T], rf) 
                : zero_v; // out of bounds: extrap 0

                // accumulate tx here: add to pixel value
                pix += val;
            }
        }
        y[tid] = pix; // output value 
    }
}

__global__ void DAS(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double * __restrict__ Pr, 
    const double * __restrict__ Pv, const double * __restrict__ Nv, 
	const double2 * __restrict__ x, 
	const double t0, const double fs, const double cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 4 x I
    const double3 * pr = reinterpret_cast<const double3*>(Pr); // 3 x N
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M
    const double3 * nv = reinterpret_cast<const double3*>(Nv); // 3 x M
    
    // temp vars
    const double2 zero_v = make_double2(0.0);
    double rf, dv, dr;
    double2 val, pix = zero_v;
    double3 rv;
    union {double f; int i;} t;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // fractional and integer part
                rf = modf(rf, &t.f);
                t.i = (int)(t.f);
                
                // if in bounds, linearly interpolate by ratio rf at time-index t.i[+1]
                val = (0 <= t.i && t.i + 1 < T) ?            
                   lerp(x[t.i + 1 + (n + m * N) * T],
                        x[t.i +     (n + m * N) * T], rf) 
                : zero_v; // out of bounds: extrap 0

                // accumulate tx here: add to pixel value
                pix += val;
            }            
        }
        y[tid] = pix; // output value
    }
}

/*
* Delay the given data at the given pixels and sum over transmits
*
* Given a set of pixels, (virtual or plane wave) transmitter locations, 
* receiver locations, as well as a datacube equipped with a time, 
* transmitter and receiver axis, perform simple delay-and-sum beamforming. 
* The data is linearly interpolated at the sample time. An image is 
* generated for each receiver element. 
*
* All positions are in vector coordinates. 
* 
* If the virtual transmitter 
* normal has a fourth component that is 0, this indicates that the 
* transmission should be treated as a plane wave transmission instead of a 
* virtual source transmission. 
* 
* The value of t = 0 must be the time when the peak of the wavefront 
* reaches the virtual source location. Because this time must be the same 
* for all transmits, the datacube must be stitched together in such a way 
* that for all transmits, the same time axis is used
*
* Inputs:
*  y:           complex pixel values per channel (N x I)
*  Pi:          pixel positions (3 x I)
*  Pr:          receiver positions (3 x N)
*  Pv:          (virtual) transmitter positions (3 x M)
*  Nv:          (virtual) transmitter normal (3 x M)
*  x:           datacube of complex sample values (T x M x N)
*  t0:          initial time for the data
*  fs:          sampling frequency of the data
*  cinv:        inverse of the speed of sound used for beamforming
* 
* I -> pixels, M -> transmitters, N -> receivers, T -> time samples
*
*/

__global__ void SYNf(float2 * __restrict__ y, 
    const float * __restrict__ Pi, const float * __restrict__ Pr, 
    const float * __restrict__ Pv, const float * __restrict__ Nv, 
	const float2 * __restrict__ x, 
	const float t0, const float fs, const float cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 4 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M
    
    // temp vars
    const float2 zero_v = make_float2(0.0f);
    float rf, dv, dr;
    float2 val;
    float3 rv;
    union {float f; int i;} t;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // fractional and integer part
                rf = modf(rf, &t.f);
                t.i = (int)(t.f);
                
                // if in bounds, linearly interpolate by ratio rf at time-index t.i[+1]
                val = (0 <= t.i && t.i + 1 < T) ?            
                   lerp(x[t.i + 1 + (n + m * N) * T],
                        x[t.i +     (n + m * N) * T], rf) 
                : zero_v; // out of bounds: extrap 0

                // accumulate tx here: add to pixel value
                y[tid + n*I] += val; // index as I x N
            }
        }
    }
}

__global__ void SYN(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double * __restrict__ Pr, 
    const double * __restrict__ Pv, const double * __restrict__ Nv, 
	const double2 * __restrict__ x, 
	const double t0, const double fs, const double cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 3 x I
    const double3 * pr = reinterpret_cast<const double3*>(Pr); // 3 x N
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M
    const double3 * nv = reinterpret_cast<const double3*>(Nv); // 3 x M
    
    // temp vars
    const double2 zero_v = make_double2(0.0);
    double rf, dv, dr;
    double2 val;
    double3 rv;
    union {double f; int i;} t;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // fractional and integer part
                rf = modf(rf, &t.f);
                t.i = (int)(t.f);
                
                // if in bounds, linearly interpolate by ratio rf at time-index t.i[+1]
                val = (0 <= t.i && t.i + 1 < T) ?            
                   lerp(x[t.i + 1 + (n + m * N) * T],
                        x[t.i +     (n + m * N) * T], rf) 
                : zero_v; // out of bounds: extrap 0

                // accumulate tx here: add to pixel value
                y[tid + n*I] += val; // index as I x N
            }
        }
    }
}


/*
* Beamform the data at the given pixels.
*
* Given a set of pixels, (virtual or plane wave) transmitter locations, 
* receiver locations, as well as a datacube equipped with a time, 
* transmitter and receiver axis, beamforming the data without summation. 
* The data is linearly interpolated at the sample time. 
*
* All positions are in vector coordinates. 
* 
* If the virtual transmitter normal has a fourth component that is 0, this 
* indicates that the transmission should be treated as a plane wave 
* transmission instead of a virtual source (focused) transmission. 
* 
* The value of t = 0 must be the time when the peak of the wavefront 
* reaches the virtual source location. Because this time must be the same 
* for all transmits, the datacube must be stitched together in such a way 
* that for all transmits, the same time axis is used.
*
* Inputs:
*  y:           complex pixel values per transmit/channel (M x N x I)
*  Pi:          pixel positions (3 x I)
*  Pr:          receiver positions (3 x N)
*  Pv:          (virtual) transmitter positions (3 x M)
*  Nv:          (virtual) transmitter normal (3 x M)
*  x:           datacube of complex sample values (T x M x N)
*  t0:          initial time for the data
*  fs:          sampling frequency of the data
*  cinv:        inverse of the speed of sound used for beamforming
* 
* I -> pixels, M -> transmitters, N -> receivers, T -> time samples
*
*/

__global__ void BFf(float2 * __restrict__ y, 
    const float * __restrict__ Pi, const float * __restrict__ Pr, 
    const float * __restrict__ Pv, const float * __restrict__ Nv, 
	const float2 * __restrict__ x, 
	const float t0, const float fs, const float cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 4 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M
    
    // temp vars
    const float2 zero_v = make_float2(0.0f);
    float rf, dv, dr;
    float2 val;
    float3 rv;
    union {float f; int i;} t;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // fractional and integer part
                rf = modf(rf, &t.f);
                t.i = (int)(t.f);
                
                // if in bounds, linearly interpolate by ratio rf at time-index t.i[+1]
                val = (0 <= t.i && t.i + 1 < T) ?            
                   lerp(x[t.i + 1 + (n + m * N) * T],
                        x[t.i +     (n + m * N) * T], rf) 
                : zero_v; // out of bounds: extrap 0

                // output value
                y[tid + n * I + m * I * N] = val; // index as I x N x M
            }            
        }
    }
}


__global__ void BF(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double * __restrict__ Pr, 
    const double * __restrict__ Pv, const double * __restrict__ Nv, 
	const double2 * __restrict__ x, 
	const double t0, const double fs, const double cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 4 x I
    const double3 * pr = reinterpret_cast<const double3*>(Pr); // 3 x N
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M
    const double3 * nv = reinterpret_cast<const double3*>(Nv); // 3 x M
    
    // temp vars
    const double2 zero_v = make_double2(0.0f);
    double rf, dv, dr;
    double2 val;
    double3 rv;
    union {double f; int i;} t;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // fractional and integer part
                rf = modf(rf, &t.f);
                t.i = (int)(t.f);
                
                // if in bounds, linearly interpolate by ratio rf at time-index t.i[+1]
                val = (0 <= t.i && t.i + 1 < T) ?            
                   lerp(x[t.i + 1 + (n + m * N) * T],
                        x[t.i +     (n + m * N) * T], rf) 
                : zero_v; // out of bounds: extrap 0

                // output value
                y[tid + n * I + m * I * N] = val; // index as I x N x M
            }            
        }
    }
}


/*
* Beamforming delays at the given pixels.
*
* Given a set of pixels, (virtual or plane wave) transmitter locations, 
* receiver locations, as well as a datacube equipped with a time, 
* transmitter and receiver axis, compute the sample times corresponding to 
* when the peak of the response from an ideal scatterer arrives at the 
* receiver. 
*
* All positions are in vector coordinates. 
* 
* If the virtual transmitter normal has a fourth component that is 0, this 
* indicates that the transmission should be treated as a plane wave 
* transmission instead of a virtual source (focused) transmission. 
* 
* The value of t = 0 must be the time when the peak of the wavefront 
* reaches the virtual source location. Because this time must be the same 
* for all transmits, the datacube must be stitched together in such a way 
* that for all transmits, the same time axis is used.
*
* Inputs:
*  tau:         sample time per transmit/channel/pixel (M x N x I)
*  Pi:          pixel positions (3 x I)
*  Pr:          receiver positions (3 x N)
*  Pv:          (virtual) transmitter positions (3 x M)
*  Nv:          (virtual) transmitter normal (3 x M)
*  t0:          initial time for the data
*  fs:          sampling frequency of the data
*  cinv:        inverse of the speed of sound used for beamforming
* 
* I -> pixels, M -> transmitters, N -> receivers, T -> time samples
*
*/

__global__ void delaysf(float * __restrict__ tau, 
    const float * __restrict__ Pi, const float * __restrict__ Pr, 
    const float * __restrict__ Pv, const float * __restrict__ Nv, 
	const float cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 3 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M
    
    // temp vars
    float dv, dr;
    float3 rv;
    
    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // output time
                tau[tid + n * I + m * I * N] = cinv * (dv + dr);
            }            
        }
    }
}


__global__ void delays(double * __restrict__ tau, 
    const double * __restrict__ Pi, const double * __restrict__ Pr, 
    const double * __restrict__ Pv, const double * __restrict__ Nv, 
	const double cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 3 x I
    const double3 * pr = reinterpret_cast<const double3*>(Pr); // 3 x N
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M
    const double3 * nv = reinterpret_cast<const double3*>(Nv); // 3 x M
    
    // temp vars
    double dv, dr;
    double3 rv;
    
    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // output time
                tau[tid + n * I + m * I * N] = cinv * (dv + dr);
            }
        }
    }
}


/*
* Amplitude at a given pixel assuming huygens principle on transmit.
*
* Given a set of pixels, (virtual or plane wave) transmitter locations, 
* receiver locations, as well as a datacube equipped with a time, 
* transmitter and receiver axis, compute the amplitude corresponding to 
* the sum of the actual transmit aperture wavefronts, assuming no 
* absorption. The amplitude at the time 
*
* All positions are in vector coordinates. 
* 
* If the virtual transmitter normal has a fourth component that is 0, this 
* indicates that the transmission should be treated as a plane wave 
* transmission instead of a virtual source (focused) transmission.
* 
* The value of t = 0 must be the time when the peak of the wavefront 
* reaches the virtual source location. Because this time must be the same 
* for all transmits, the datacube must be stitched together in such a way 
* that for all transmits, the same time axis is used.
*
* Inputs:
*  Ai:          amplitude as predicted by huygens at the pixel (M x I)
*  Pi:          pixel positions (3 x I)
*  Pv:          virtual transmitter positions (3 x M)
*  Nv:          virtual transmitter normals (3 x M)
*  Pt:          physical transmitter positions (3 x N)
*  At:          transmit apodization (N)
*  s:           complex transmitted signal (S)
*  t0:          initial time for the signal
*  fs:          sampling frequency of the signal
*  cinv:        inverse of the speed of sound
* 
* I -> pixels, M -> virtual transmitters, N -> physical transmitters, T -> time samples
*
*/

__global__ void huygen_ampf(float * __restrict__ Ai, 
    const float * __restrict__ Pi, const float * __restrict__ Pv,
    const float * __restrict__ Nv, const float * __restrict__ Pt, 
    const float * __restrict__ At, const float * __restrict__ s,
    const float t0, const float fs, const float cinv){

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 3 x I
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M
    const float3 * pt = reinterpret_cast<const float3*>(Pt); // 3 x N
    const float2 * sj = reinterpret_cast<const float2*>(s);  // S
    
    // temp vars
    float dv, dt, dtau;
    const float2 zero_v = make_float2(0.0f);
    float2 v = zero_v;
    float3 rv;
    union {float f; int i;} t;
    
    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            // virtual transmit to pixel ray
            rv = pi[tid] - pv[m];
            dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
            // physical transmit to pixel ray
            for(size_t n = 0; n < N; ++n){
                // ray length 
                dt = length(pi[tid] - pt[n]);
                
                // time difference in input signal 
                dtau = (cinv * (dv - dt) - t0) * fs;
                
                // fractional and integer part
                dtau = modf(dtau, &t.f);
                t.i = (int)(t.f);
                
                // add value
                v += (0 <= t.i && t.i + 1 < S) ? 
                   (At[m] * lerp(sj[t.i + 1], sj[t.i], dtau)) // in bounds: linear interp
                : (zero_v); // out of bounds: extrap 0
            }
            
            // output amplitude
            Ai[m + M * tid] = length(v);
        }
    }
}


__global__ void huygen_amp(double * __restrict__ Ai, 
    const double * __restrict__ Pi, const double * __restrict__ Pv,
    const double * __restrict__ Nv, const double * __restrict__ Pt, 
    const double * __restrict__ At, const double * __restrict__ s,
    const double t0, const double fs, const double cinv){

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 3 x I
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M
    const double3 * nv = reinterpret_cast<const double3*>(Nv); // 3 x M
    const double3 * pt = reinterpret_cast<const double3*>(Pt); // 3 x N
    const double2 * sj = reinterpret_cast<const double2*>(s);  // S
    
    // temp vars
    double dv, dt, dtau;
    const double2 zero_v = make_double2(0.0f);
    double2 v = zero_v;
    double3 rv;
    union {double f; int i;} t;
    
    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            // virtual transmit to pixel ray
            rv = pi[tid] - pv[m];
            dv = VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
            // physical transmit to pixel ray
            for(size_t n = 0; n < N; ++n){
                // ray length 
                dt = length(pi[tid] - pt[n]);
                
                // time difference in input signal 
                dtau = (cinv * (dv - dt) - t0) * fs;
                
                // fractional and integer part
                dtau = modf(dtau, &t.f);
                t.i = (int)(t.f);
                
                // add complex value
                v += (0 <= t.i && t.i + 1 < S) ? 
                   At[m] * lerp(sj[t.i + 1], sj[t.i], dtau) // in bounds: linear interp
                : zero_v; // out of bounds: extrap 0
            }
            
            // output amplitude
            Ai[m + M * tid] = length(v);
        }
    }
}



/*
* Estimate the phase abberation at the given pixels for all element pairs.
*
* Given a set of pixels, (virtual or plane wave) transmitter locations, 
* receiver locations, as well as a datacube equipped with a time, 
* transmitter and receiver axis, perform simple delay-and-sum beamforming. 
* The data is linearly interpolated at the sample time. An image is 
* generated for each receiver element. Summation across the transmitters is
* implicit.
*
* All positions are in vector coordinates. 
* 
* If the virtual transmitter 
* normal has a fourth component that is 0, this indicates that the 
* transmission should be treated as a plane wave transmission instead of a 
* virtual source transmission. 
* 
* The value of t = 0 must be the time when the peak of the wavefront 
* reaches the virtual source location. Because this time must be the same 
* for all transmits, the datacube must be stitched together in such a way 
* that for all transmits, the same time axis is used
*
* Inputs:
*  r:           complex phasor values per channel pair (N x I)
*  Pi:          pixel positions (3 x I)
*  Pr:          receiver positions (3 x N)
*  Pv:          (virtual) transmitter positions (3 x M)
*  Nv:          (virtual) transmitter normal (3 x M)
*  x:           datacube of complex sample values (T x M x N)
*  Rd:          spatial integration offsets
*  Gd:          spatial integration weights
*  t0:          initial time for the data
*  fs:          sampling frequency of the data
*  cinv:        inverse of the speed of sound used for beamforming
* 
* I -> pixels, M -> transmitters, N -> receivers, T -> time samples
*
*/
        
