#include <complex.h>

/* Structure for a distribution fcn */
typedef struct {
	int nL;
	
	complex double **dist;
} Dist;

/* Helper structure that keeps track of 
non-zero elements in a distribution fcn.*/
typedef struct{
	int num;
	
	int *is;
	int *js;
	int *ks;
	int *cnts;
} NonZeros;

/* Storage structure for the r-r' terms in 
the forward propagator*/
typedef struct {
	int nL;
	// These are different from nX nY and nZ
	int NX,NY,NZ;
	
	complex double **Ylm;
	complex double *fact;
} Terms;

/* Storage structure for the r-r' terms in 
the forward propagator
New Attempt at these terms*/
typedef struct {
	int nL;
	// These are different from nX nY and nZ
	int NX,NY,NZ;
	
	complex double ***fact;
} Terms2;

/* Information about the geometry */
typedef struct {
	int nX;		// number of voxels
	int nY;
	int nZ;
	
	double *xs;	// center positions of voxels
	double *ys;
	double *zs;
	
	double xMin;	// minimum edge
	double xMax;	// maximum edge
	double yMin;
	double yMax;
	double zMin;
	double zMax;
	
	double delX;	// voxel size
	double delY;
	double delZ;
	
	int subVox;
	int self_sub_vox;
	int subR;
	int subTheta;
	int subPhi;
	double subThresh;
	double propThresh;
} Geometry;

typedef struct {
	double posX;
	double posY;
	double posZ;
	double diam;
	double theta;
	double phi;
	double mag;
	double modulation;  // in MHz	
} Source;

/* Information about the uniform phantom */
typedef struct{
	complex double mua;
	complex double mus;
	double g;
	double n;
	double c;
} Phantom;

#define THRESHOLD 0.0
#define NSTOP 0.9999 // Neumann series stopping criterion (0.0<NSTOP<1.0)


/* All the function prototypes are below */
void PropScat(Geometry *geom,Phantom *phan,int nL,Dist *src,Dist *W);
//void PropScatmu1(Geometry *geom,Phantom *phan,int nL,Dist *src);
//void PropAbs(Geometry *geom,Terms *terms,int nL,Dist *src,Dist *W);
void PropAbs2(Geometry *geom,Terms2 *terms,int nL,Dist *src,Dist *W);
//Terms *GenerateTerms(Geometry *geom,Phantom *phan,int nL);
//void FreeTerms(Terms *fac);
Terms2 *GenerateTerms2(Geometry *geom,Phantom *phan,int nL);
void FreeTerms2(Terms2 *fac);
//NonZeros *GenerateNonZeros(Geometry *geom,Dist *dist);
//void FreeNonZeros(NonZeros *zer);
int StopIteration(Geometry *geom,double *fi,int n,Dist *W,Dist *out);

Geometry *LoadGeometry(dictionary *ini);
Geometry *LoadGeometryZPlane(dictionary *ini);
void FreeGeometry(Geometry *geom);
Source *LoadSource(dictionary *ini);
void FreeSource(Source *src);
Phantom *LoadPhantom(dictionary *ini);
void FreePhantom(Phantom *phan);

void OutputImage(Geometry *geom,complex double *transImg,char *name);
void OutputPhaseImage(Geometry *geom,double *transImg,char *name);
double *ComputePhaseImage(Geometry *geom,Phantom *phan,Source *source,int nL,int num,Dist *W);
complex double *ComputeTransImage(Geometry *geom,Phantom *phan,Source *source,int nL,int num,Dist *W);
//complex double *ComputeTransImageNoSource(Geometry *geom,Phantom *phan,int nL,int num,Dist *W);
//complex double *ComputeTransImageLW(Geometry *geom,Phantom *phan,Source *source,int nL,int num,Dist *W);
complex double *ComputeReflImage(Geometry *geom,Phantom *phan,int nL,int num,Dist *W);
complex double ComputePhotonNumber(Geometry *geom,Phantom *phan,Source *source,int nL,int num,Dist *W);
void ComputeSideFaceImages(Geometry *geom,Phantom *phan,int nL, int num, Dist *W);
//void OutputDist(Geometry *geom,Dist *dist,char *name);
//void GenerateSource(Geometry *geom,Source *source,int nL,Dist *src);
//void GenerateSPECTSource(Geometry *geom,Phantom *phan,Source *source,int nL,Dist *src);
void GenerateSourceBeam(Geometry *geom,Phantom *phan,Source *source,int nL,Dist *src);
Dist *AllocDist(Geometry *geom,int nL);
void FreeDist(Dist *dst);
void AddDist(Geometry *geom,int nL,Dist *W,Dist *out);
void SubDist(Geometry *geom,int nL,Dist *W,Dist *W1,Dist *out);
void CopyDist(Geometry *geom,int nL,Dist *W,Dist *out);
void CombineOuts(Geometry *geom,int nL,Dist *out1,Dist *out2);
void ScaleDist(Geometry *geom,int nL,Dist *W,double factor);
void AdjModify(Geometry *geom, int nL,Dist *W);
void DotDist(Geometry *geom,int nL,Dist *W,Dist *W1,Dist *out);

complex double SpherHarmonic(double theta,double phi,int degree,int order);
void SpherHarmonicArray(int degree,int order,int N,double *theta,double *phi,complex double *out);
double AssocLegendre(double x,int l, int m);
void AssocLegendreArray(int degree, int order,int N,double *in,double *out);
//complex double TestSPH(int l,int m,int lp,int mp);
