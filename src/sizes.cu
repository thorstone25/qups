// data sizes: I -> pixels, M -> transmitters, N -> receivers, 
// T -> time samples, F -> Frames, S -> signal
// I[1-3] -> pixel dimension sizing
        

#ifndef QUPS_SIZES
#define QUPS_SIZES

// integral params
# ifndef W
__constant__ size_t W; // temporal integration width (one-sided)
# endif
# ifndef D
__constant__ size_t D; // number of points in the spatial integration
# endif

// Beamforming transmit mode: focal point or plane-wave
# ifndef QUPS_VS
__constant__ bool QUPS_VS; // whither virtual source mode
# endif
        
# ifndef QUPS_T
__constant__ size_t QUPS_T;
# endif
# ifndef QUPS_M
__constant__ size_t QUPS_M;
# endif
# ifndef QUPS_N
__constant__ size_t QUPS_N;
# endif
# ifndef QUPS_I
__constant__ size_t QUPS_I;
# endif
# ifndef QUPS_F
__constant__ size_t QUPS_F;
# endif
# ifndef QUPS_S
__constant__ size_t QUPS_S;
# endif

# ifndef QUPS_I1
__constant__ size_t QUPS_I1;
# endif
# ifndef QUPS_I2
__constant__ size_t QUPS_I2;
# endif
# ifndef QUPS_I3
__constant__ size_t QUPS_I3;
# endif
# endif // QUPS_SIZES
