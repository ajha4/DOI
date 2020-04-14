#define C 1
//2.9979e10
#define INF 10000000 /* A large number */
#define SMALL_NUM 0.0000000001 /* A very small number */

#define ANG_RES 100 // The angular resolution, used while computing the self sub voxel term
#define BLK_ANG_SIZE 5
#define PHI_ANG_RES ANG_RES
#define THETA_ANG_RES ANG_RES

#ifdef __cplusplus
    #define CHECK_EXT extern "C"
#else
    #define CHECK_EXT
#endif

#include "complex_arith.h"
#include "iniparser.h"
#include "defines.h"

#define MAX_TISS_TYPE 4
#define MAX_TISS_NUM 4
#define MAX_NL 5
#define BLK_SIZE 10
#define MAX_PROP 0 // No extra voxels on boundary. 
#define BLK_SIZE_Z 5
#define BLK_SRC_SIZE (BLK_SIZE + 2*MAX_PROP) *(BLK_SIZE + 2*MAX_PROP)*(BLK_SIZE + 2*MAX_PROP)

#define BLK_SELF_SUB_VOX 10
#define MY_PREF_DEVICE 0
#define NUM_DEVICES 4


typedef double flt_doub;
typedef cudacomplex complex_double;
typedef cudacomplex doublecomplex;
typedef short SHORT;
typedef short SIGNED_SHORT;
typedef unsigned char byte;

typedef struct {
    SHORT nX, nY,nZ; /* No of voxels in the phantom in each dimension */
    int no_vox; /* Total No of voxels */
    float x_min, y_min, z_min; /* start and end of the phantom in the r space*/
    float x_max, y_max, z_max;
    flt_doub *xs, *ys, *zs; /* Center coordinates of each voxel in the r space*/
    flt_doub delX, delY, delZ; /* The dimensions of each voxel */
    float sub_thresh;
    float prop_thresh;
    SHORT sub_vox;
    SHORT self_sub_vox;
    float correction_fact;
    SHORT bounX, bounY, bounZ; /* The boundary pading that we have to do since we have a SIMD model. */
    SHORT subbounX, subbounY,subbounZ;

}Geometry;

typedef struct {
    flt_doub g;   /* The scattering anisotropy parameter */
    flt_doub n;   /* The refractive index of the medium */
    byte *tiss_type; /* The tissue type*/
	SHORT no_tiss;
	doublecomplex *mu_abs, *mu_sc;
} Phantom;

typedef struct {
    flt_doub pos_x,pos_y,pos_z; /* Position of the source in the r space */
    flt_doub diam;
    flt_doub theta;
    flt_doub phi;
    flt_doub mag;
    flt_doub modulation;  // in MHz
} Source;

typedef struct {
    SHORT nX, nY,nZ; /* No of voxels in the phantom in each dimension */
    int no_vox; /* Total No of voxels */
    float x_min, y_min, z_min; /* start and end of the phantom in the r space*/
    float x_max, y_max, z_max;
	flt_doub delX, delY, delZ; /* The dimensions of each voxel */
    float sub_thresh;
    float prop_thresh;
    SHORT sub_vox;
    SHORT self_sub_vox;
    float correction_fact;
    SHORT bounX, bounY, bounZ;
	flt_doub g;
	flt_doub n;
	SHORT no_tiss;
    float cm; // Effective speed of light in the media. Again computed and stored.
    SHORT subbounZ, subbounY, subbounX; 
	//doublecomplex mu_abs[MAX_TISS_TYPE]; // = {0.0, 0.0, 0.0, 0.0};
	//doublecomplex mu_sc[MAX_TISS_TYPE];// = {0.0, 0.0, 0.0, 0.0};
	complex_double mu_tot[MAX_TISS_NUM];
	complex_double mu_sc[MAX_TISS_NUM];
    SHORT boun_blkX, boun_blkY,boun_blkZ;
    SHORT boun_blk_sizeZ, boun_blk_sizeY, boun_blk_sizeX;
} Info_Stat;

typedef struct{
	/*doublecomplex *mu_abs;
	doublecomplex *mu_sc; */
	byte *tiss_type;
}Info_Dyn;

typedef struct{
    SHORT i;
    SHORT j;
    SHORT k;
} COR;

typedef struct {
	int device_index;
    complex_double *src_host;
	complex_double *out_host;
	SHORT num_layers;
	SHORT layer_start;
	SHORT flag; // The flag to tell whether all the spherical harmonics should be evaluated, or only nL = 0, for the g = 0 case. 
} THREAD_PARAMETERS;


void generate_source_beam(complex_double *src_dist);
complex_double * generate_trans_image(complex_double *out_dist, int flag );
double * generate_phase_image(complex_double *out_dist, int flag );
complex_double * generate_ref_image(complex_double *out_dist);

complex_double* alloc_dist();
complex_double* alloc_dist_dev();
void add_dist( complex_double *W1, complex_double *W2, complex_double *out);
void copy_dist(complex_double *W1, complex_double *W2);
void scale_dist(complex_double *W1, double fac);
complex_double * generate_trans_image(complex_double *out_dist);
void output_image(complex_double *trans_img,char *name);
void output_phase_image(double *phase_img,char *name);

Geometry *LoadGeometry(dictionary *ini);
Source *LoadSource(dictionary *ini);
Phantom *LoadPhantom(dictionary *ini, int flag);
void PropScatmu1(Geometry *geom,Phantom *phan,int nL,complex_double *src);

void AssocLegendreArray(int degree, int order,int N,flt_doub *in,flt_doub *out);
flt_doub AssocLegendre(flt_doub x,int degree, int order);

void FreeSource();
void FreeGeometry();
void FreePhantom();

/* Functions that could be invoked from both host and device */

__global__ void add_dist_dev( complex_double *W1, complex_double *W2, complex_double *out, Geometry *geom, int nL);
__global__ void copy_dist_dev (doublecomplex *W2,doublecomplex *W1);
void prop_scat_dev (doublecomplex *src_dist_dev,doublecomplex *out_dist_dev, Info_Dyn *);

void prop_scat(doublecomplex *src_dist_dev,doublecomplex *out_dist_dev);

void checkCUDAError(const char *msg);
//__global__ void copy_dist_test(complex_double *B, complex_double *A);

__device__ flt_doub AssocLegendre_dev(flt_doub x,int degree, int order);
__device__ complex_double SpherHarmonic_dev(flt_doub theta,flt_doub phi,int degree,int order);
__global__ void generate_diag_terms_dev(Geometry *geom, Phantom *phan, int nL, complex_double *diag_terms);

void SpherHarmonicArray(int degree,int order,int N,flt_doub *theta,flt_doub *phi,complex_double *out); 
void generate_diag_terms_host();
void generate_diag_terms_2();
void generate_sph_harm_terms();

void Neumann(complex_double* src_host, complex_double *out_host);
//__device__ complex_double mu_int (COR * cor_src, COR * cor_dest, COR * subvox_src, COR * subvox_dest, Info_Dyn *info_dyn);

//__global__ void prop_abs (doublecomplex *src_dist, doublecomplex *out_dist,Geometry *geom, Phantom *phan, int nL, complex_double *diag_terms);
void * prop_abs(void *arg);
complex_double ComputePhotonNumber(int num,complex_double *W);
