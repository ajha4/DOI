#include<stdio.h>

#define CUDA_SUCCESS cudaSuccess

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define SIGN(x) ( (x > 0) ? 1 : -1 )

#ifdef __CUDACC__
#define MY_GLOBAL __global__
#define MY_DEVICE __device__
#else
#define GLOBAL 
#define MY_DEVICE _
#endif

// Macro used to check CUDA function returning value and error
#define MY_SAFE_CALL(call) \
 do{ \
   int err = call; \
   if(err != CUDA_SUCCESS) { \
     fprintf(stderr, "CUDA driver error %x in file %s in line %i.\n",err, __FILE__, __LINE__);          \
     fprintf(stderr, "CUDA Error message: %s.\n",cudaGetErrorString((cudaError_t) err));                \
     exit(EXIT_FAILURE); \
 	} \
} while(0)


// Macro to check the result of calloc(...), malloc(...), and realloc(...)
#define CHECK_ALLOC(call) \
 do{ \
 void *pointer = call; \
 if(pointer == NULL) { \
   fprintf(stderr, "ERROR: Can't allocate memory!\n"); \
   fprintf(stderr, "(File %s, line %d)\n", __FILE__, __LINE__); \
   exit(1); \
 } \
} while(0)

#define NORM_THETA(__ang) (__ang > 0 ? __ang : 2*M_PI + __ang)
#define NORM_PHI(__ang) (__ang > 0 ? __ang : 2*M_PI + __ang)


//#define THETA_ANG_IND(x) (max((int)rintf((THETA_ANG_RES-1) * (x +M_PI)/(2*M_PI)),0))
//#define PHI_ANG_IND(x)  (max((int)rintf((PHI_ANG_RES-1) * (x +M_PI)/(2*M_PI)),0))

#define THETA_ANG_IND(x) ((int)rintf((THETA_ANG_RES-1) * (x )/(M_PI)))
#define PHI_ANG_IND(x)   ((int)rintf((PHI_ANG_RES-1) * (x )/(2*M_PI)))

#define SPH_HARM_IND(__cnt,__theta,__phi) (__cnt * (THETA_ANG_RES * PHI_ANG_RES) + THETA_ANG_IND(__theta)*PHI_ANG_RES + PHI_ANG_IND(__phi))
/* Converts from the voxel position and spherical harmonic value to a linear index */
//#define VOX_TO_SPIND(__i, __j,__nvox)   (__j * __nvox + __i )
#define VOX_TO_SPIND(__i, __j,__ncoeff)   (__i * __ncoeff + __j )
#define DIAG_VOX_TO_SPIND(__i, __j, __k, __nL)  (__i* (__nL+1)*(__nL+1)*(__nL+1)*(__nL+1) + __j * (__nL+1)*(__nL+1) + __k)

#define SPH_HARM_VOX_TO_IND(__i, __j)  (__i*THETA_ANG_RES * PHI_ANG_RES + __j)
#define SPH_TO_IND(__i, __j)  (__i*ANG_RES*ANG_RES + __j)

#define SPIND_TO_VOX(__r, __cnt, __spind, __nvox) do{ \
 __r = (__spind % __nvox); \
 __cnt = (__spind / __nvox); \
 } while(0)

#define RIND_TO_PHANIND(__i, __j, __k, __r,__stridej, __stridek, __offseti,__offsetj, __offsetk) do{ \
 __k = (__r % (__stridek) - __offsetk); \
 __j = (((__r /__stridek )% __stridej) - __offsetj); \
 __i = (__r / (__stridek*__stridej) - __offseti); \
} while(0) 

/* Converts the voxel linear index, without the boundary of the phantom, to the corresponding x,y and z indices, again without the boundary */
#define TIND_TO_VOX_Z (__ind)  (__ind%(geom->nZ)) 
#define TIND_TO_VOX_Y (__ind)  (__ind/(geom->nZ)%geom->nY)
#define TIND_TO_VOX_X (__ind)  (__ind/(geom->nZ*geom->nY))

#define TIND_TO_VOX_XYZ(__ind, __i, __j, __k, __nZ, __nY, __nX) do{ \
__i = (__ind/(__nX*__nY)); \
__j = (__ind -__i * (__nX*__nY)) / __nX ; \
__k = ( __ind - __i * (__nX*__nY) - __j * __nX );\
} while (0)
 
