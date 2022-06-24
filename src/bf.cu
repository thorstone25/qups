// # include "helper_math.h" // vector math

// # include "sizes.cu" // size defines

# include "interpd.cu" // samplers using constant sizing

// # include "half2_math.h" // vector math for half types only 


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
	const float2 * __restrict__ a, const size_t * astride, 
    const float2 * __restrict__ x, const int iflag,
	const float t0, const float fs, const float cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // get image coordinates
    const size_t I1 = QUPS_I1, I2 = QUPS_I2, I3 = QUPS_I3; // rename for readability
    const size_t i1 = (tid             % I1); // index in I1
    const size_t i2 = (tid /  I1     ) % I2 ; // index in I2
    const size_t i3 = (tid / (I1 * I2) % I3); // index in I3
    const size_t abase = i1 * astride[0] + i2 * astride[1] + i3 * astride[2]; // base index for this pixel

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 3 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M
    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, T = QUPS_T, I = QUPS_I;
            
    // temp vars
    const float2 zero_v = make_float2(0.0f);
    float2 w;
    float rf, dv, dr, tau;
    float2 val, pix = zero_v;
    float3 rv;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                tau = (cinv * (dv + dr) - t0);
                w = make_float2(1.0f, 0.0f); // TODO: enable demod: make_float2(cospi(2*fs/4*tau), -sinpi(2*fs/4*tau));
                rf =  tau * fs;

                // sample the trace
                val = samplef(&x[(n + m * N) * T], rf, iflag, zero_v); // out of bounds: extrap 0

                // accumulate tx here: add to pixel value
                pix += val * w * a[abase + n * astride[3] + m * astride[4]];
            }
        }
        y[tid] = pix; // output value 
    }
}

__global__ void DAS(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double * __restrict__ Pr, 
    const double * __restrict__ Pv, const double * __restrict__ Nv, 
    const double2 * __restrict__ a, const size_t * astride,
	const double2 * __restrict__ x, const int iflag,
	const double t0, const double fs, const double cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // get image coordinates
    const size_t I1 = QUPS_I1, I2 = QUPS_I2, I3 = QUPS_I3; // rename for readability
    const size_t i1 = (tid             % I1); // index in I1
    const size_t i2 = (tid /  I1     ) % I2 ; // index in I2
    const size_t i3 = (tid / (I1 * I2) % I3); // index in I3
    const size_t abase = i1 * astride[0] + i2 * astride[1] + i3 * astride[2]; // base index for this pixel

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 4 x I
    const double3 * pr = reinterpret_cast<const double3*>(Pr); // 3 x N
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M
    const double3 * nv = reinterpret_cast<const double3*>(Nv); // 3 x M
    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, T = QUPS_T, I = QUPS_I;
            
    // temp vars
    const double2 zero_v = make_double2(0.0);
    double2 w;
    double rf, dv, dr, tau;
    double2 val, pix = zero_v;
    double3 rv;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                tau = (cinv * (dv + dr) - t0);
                w = make_double2(1.0,0.0); // TODO: enable demod: make_double2(cospi(2*fs/4*tau), -sinpi(2*fs/4*tau));
                rf =  tau * fs;

                // sample the trace
                val = sample(&x[(n + m * N) * T], rf, iflag, zero_v); // out of bounds: extrap 0

                // accumulate tx here: add to pixel value
                pix += val * w * a[abase + n * astride[3] + m * astride[4]];
            }            
        }
        y[tid] = pix; // output value
    }
}

__global__ void DASh(ushort2 * __restrict__ y, 
    const float * __restrict__ Pi, const float * __restrict__ Pr, 
    const float * __restrict__ Pv, const float * __restrict__ Nv, 
	const ushort2 * __restrict__ a, const size_t * astride, 
    const ushort2 * __restrict__ x, const int iflag,
	const float t0, const float fs, const float cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // get image coordinates
    const size_t I1 = QUPS_I1, I2 = QUPS_I2, I3 = QUPS_I3; // rename for readability
    const size_t i1 = (tid             % I1); // index in I1
    const size_t i2 = (tid /  I1     ) % I2 ; // index in I2
    const size_t i3 = (tid / (I1 * I2) % I3); // index in I3
    const size_t abase = i1 * astride[0] + i2 * astride[1] + i3 * astride[2]; // base index for this pixel

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 3 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M

    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, T = QUPS_T, I = QUPS_I;
            
    // temp vars
    const half2 zero_v = make_half2(0.f, 0.f);
    half2 w;
    float rf, dv, dr, tau;
    half2 val, pix = zero_v;
    float3 rv;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                tau = (cinv * (dv + dr) - t0);
                w = make_half2(1.0f, 0.0f); // TODO: enable demod: make_half2(cospi(2*fs/4*tau), -sinpi(2*fs/4*tau));
                rf =  tau * fs;

                // sample the trace
                val = sampleh(&x[(n + m * N) * T], rf, iflag, zero_v); // out of bounds: extrap 0

                // accumulate tx here: add to pixel value
                pix += val * w * u2h(a[abase + n * astride[3] + m * astride[4]]);
            }
        }
        y[tid] = h2u(pix); // output value 
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
    const float2 * __restrict__ a, const size_t * astride,
    const float2 * __restrict__ x, const int iflag,
	const float t0, const float fs, const float cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // get image coordinates (Isz)
    const size_t I1 = QUPS_I1, I2 = QUPS_I2, I3 = QUPS_I3; // rename for readability
    const size_t i1 = (tid             % I1); // index in I1
    const size_t i2 = (tid /  I1     ) % I2 ; // index in I2
    const size_t i3 = (tid / (I1 * I2) % I3); // index in I3
    const size_t abase = i1 * astride[0] + i2 * astride[1] + i3 * astride[2]; // base index for this pixel

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 4 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M
    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, T = QUPS_T, I = QUPS_I;
            
    // temp vars
    const float2 zero_v = make_float2(0.0f);
    float rf, dv, dr;
    float2 val;
    float3 rv;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // sample the trace
                val = samplef(&x[(n + m * N) * T], rf, iflag, zero_v); // out of bounds: extrap 0

                // apply apodization
                val *= a[abase + n * astride[3] + m * astride[4]]; // index as I x N

                // accumulate tx here: add to pixel value
                y[tid + n*I] += val;
            }
        }
    }
}

__global__ void SYN(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double * __restrict__ Pr, 
    const double * __restrict__ Pv, const double * __restrict__ Nv, 
	const double2 * __restrict__ a, const size_t * astride,
	const double2 * __restrict__ x, const int iflag,
	const double t0, const double fs, const double cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // get image coordinates
    const size_t I1 = QUPS_I1, I2 = QUPS_I2, I3 = QUPS_I3; // rename for readability
    const size_t i1 = (tid             % I1); // index in I1
    const size_t i2 = (tid /  I1     ) % I2 ; // index in I2
    const size_t i3 = (tid / (I1 * I2) % I3); // index in I3
    const size_t abase = i1 * astride[0] + i2 * astride[1] + i3 * astride[2]; // base index for this pixel

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 3 x I
    const double3 * pr = reinterpret_cast<const double3*>(Pr); // 3 x N
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M
    const double3 * nv = reinterpret_cast<const double3*>(Nv); // 3 x M
    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, T = QUPS_T, I = QUPS_I;
            
    // temp vars
    const double2 zero_v = make_double2(0.0);
    double rf, dv, dr;
    double2 val;
    double3 rv;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // sample the trace
                val = sample(&x[(n + m * N) * T], rf, iflag, zero_v); // out of bounds: extrap 0

                // accumulate tx here: add to pixel value
                y[tid + n*I] += val * a[abase + n * astride[3] + m * astride[4]]; // index as I x N
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
	const float2 * __restrict__ a, const size_t * astride,
	const float2 * __restrict__ x, const int iflag,
	const float t0, const float fs, const float cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // get image coordinates
    const size_t I1 = QUPS_I1, I2 = QUPS_I2, I3 = QUPS_I3; // rename for readability
    const size_t i1 = (tid             % I1); // index in I1
    const size_t i2 = (tid /  I1     ) % I2 ; // index in I2
    const size_t i3 = (tid / (I1 * I2) % I3); // index in I3
    const size_t abase = i1 * astride[0] + i2 * astride[1] + i3 * astride[2]; // base index for this pixel

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const float3 * pi = reinterpret_cast<const float3*>(Pi); // 4 x I
    const float3 * pr = reinterpret_cast<const float3*>(Pr); // 3 x N
    const float3 * pv = reinterpret_cast<const float3*>(Pv); // 3 x M
    const float3 * nv = reinterpret_cast<const float3*>(Nv); // 3 x M
    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, T = QUPS_T, I = QUPS_I;
            
    // temp vars
    const float2 zero_v = make_float2(0.0f);
    float rf, dv, dr;
    float2 val;
    float3 rv;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // sample the trace
                val = samplef(&x[(n + m * N) * T], rf, iflag, zero_v); // out of bounds: extrap 0

                // output value
                y[tid + n * I + m * I * N] = val * a[abase + n * astride[3] + m * astride[4]]; // index as I x N x M
            }            
        }
    }
}


__global__ void BF(double2 * __restrict__ y, 
    const double * __restrict__ Pi, const double * __restrict__ Pr, 
    const double * __restrict__ Pv, const double * __restrict__ Nv, 
	const double2 * __restrict__ a, const size_t * astride,
	const double2 * __restrict__ x, const int iflag,
	const double t0, const double fs, const double cinv) {

    // get starting index of this pixel
    const size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // get image coordinates
    const size_t I1 = QUPS_I1, I2 = QUPS_I2, I3 = QUPS_I3; // rename for readability
    const size_t i1 = (tid             % I1); // index in I1
    const size_t i2 = (tid /  I1     ) % I2 ; // index in I2
    const size_t i3 = (tid / (I1 * I2) % I3); // index in I3
    const size_t abase = i1 * astride[0] + i2 * astride[1] + i3 * astride[2]; // base index for this pixel

    // reinterpret inputs as vector pointers (makes loading faster and indexing easier)
    const double3 * pi = reinterpret_cast<const double3*>(Pi); // 4 x I
    const double3 * pr = reinterpret_cast<const double3*>(Pr); // 3 x N
    const double3 * pv = reinterpret_cast<const double3*>(Pv); // 3 x M
    const double3 * nv = reinterpret_cast<const double3*>(Nv); // 3 x M
    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, T = QUPS_T, I = QUPS_I;
            
    // temp vars
    const double2 zero_v = make_double2(0.0f);
    double rf, dv, dr;
    double2 val;
    double3 rv;

    // if valid pixel, for each tx/rx
    if(tid < I){
        # pragma unroll
        for(size_t m = 0; m < M; ++m){
            # pragma unroll
            for(size_t n = 0; n < N; ++n){
                // 2-way virtual path distance
                rv = pi[tid] - pv[m]; // (virtual) transmit to pixel vector 
                
                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // data/time index number
                rf = (cinv * (dv + dr) - t0) * fs;

                // sample the trace
                val = sample(&x[(n + m * N) * T], rf, iflag, zero_v); // out of bounds: extrap 0

                // output value
                y[tid + n * I + m * I * N] = val * a[abase + n * astride[3] + m * astride[4]]; // index as I x N x M
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
    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, I = QUPS_I;

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
                
                dv = QUPS_VS ? // tx path length
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
    
    // rename for readability
    const size_t N = QUPS_N, M = QUPS_M, I = QUPS_I;

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
                
                dv = QUPS_VS ? // tx path length
                    copysign(length(rv), dot(rv, nv[m])) // virtual source
                    : dot(rv, nv[m]); // plane wave
                
                dr = length(pi[tid] - pr[n]); // rx path length

                // output time
                tau[tid + n * I + m * I * N] = cinv * (dv + dr);
            }
        }
    }
}
