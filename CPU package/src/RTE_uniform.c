#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "iniparser.h"
#include "RTE.h"

#pragma warning(disable:981)

// Questions:
//   Should I include the delx * dely * delz in the propagator?  
//

int main(int argc, char** argv )
{
	dictionary *ini;
	Geometry *geom;
	Phantom *phan;
	Source *source;
	complex double *transImg;
	complex double *refImg;
	Dist *src,*W,*W1,*out,*out2;
	int nL;   // Degree of spherical harmonics
	Terms2 *terms2;
	char *bottomName;
	char *topName;
	char *phaseName;
	char *distName;
	double *phaseImg;
	int n,nTerms;
	double fi;
	
	if (argc != 2) {
		printf("\n RTE_uniform file.ini\n");
		printf("     file.ini is the initialization.\n\n");
		exit(1);
	}
	
	// Load in the initialization file
	ini = iniparser_load(argv[1]);
		
	// Set up the geometry, phantom, etc
	printf("Loading in geometry and phantom information...\n");
	geom = LoadGeometry(ini);	
	phan = LoadPhantom(ini);
	source = LoadSource(ini);
	
	nL = iniparser_getint(ini,"Algorithm:nL",-1);
	if (nL == -1) {
		printf("Change [SphericalHarmonics] to [Algorithm] for nL term\n");
		exit(1);
	}
	nTerms = iniparser_getint(ini,"Algorithm:nTerms",1);
	bottomName = iniparser_getstring(ini, "Runtime:TransFile", NULL);
	topName = iniparser_getstring(ini, "Runtime:ReflectFile", NULL);
	phaseName = iniparser_getstring(ini,"Runtime:PhaseFile",NULL);
	distName = iniparser_getstring(ini,"Runtime:DistFile","FinalDist.out");
	
	// Fixed values
	geom->self_sub_vox = 1;
	geom->subVox = 0;
	geom->subThresh = 0.0;
			
	// The first run is for a real muabs
	out = AllocDist(geom,nL);	
	out2 = AllocDist(geom,nL);

	// Create the source
	src = AllocDist(geom,nL);
	W = AllocDist(geom,nL);
	W1 = AllocDist(geom,nL);	
	
	// The following function actually does the first propagation and scatter
	printf("Computing the time-independent solution...\n");
	printf("Creating the source...\n");
	GenerateSourceBeam(geom,phan,source,nL,src);
	AddDist(geom,nL,src,out);
	AddDist(geom,nL,src,W);
	
	// Generate the Ylm matrix to speed up the PropAbs function
	printf("Generating terms for propagator...\n");
	terms2 = GenerateTerms2(geom,phan,nL);

	// Old version of Neumann series iteration
	/*printf("Neuman series for %d terms...\n",nTerms);
	for (n=0;n<nTerms-1;n++) {
		PropAbs2(geom,terms2,nL,W,W1);
                PropScat(geom,phan,nL,W1,W);
                AddDist(geom,nL,W,out2);
	}*/

	// Automatic Neumann series iteration
	printf("Neuman series iteration\n");
	fi=1.0;
	n=1;
	while( !StopIteration(geom,&fi,n,W,out) ){
		if( phan->g ){
                        PropAbs2(geom,terms2,nL,W,W1);
                        PropScat(geom,phan,nL,W1,W);
                } else{
                        PropAbs2(geom,terms2,0,W,W1);
                        PropScat(geom,phan,0,W1,W);
                }
                AddDist(geom,nL,W,out);
		n++;
	}
	// Final x-ray transform
	PropAbs2(geom,terms2,nL,out,W);
	printf("Neumann series: %d\n",n);

	FreeTerms2(terms2);
	CopyDist(geom,nL,W,out);
	FreeDist(W);
	FreeDist(W1);
	FreeDist(src);

	if (source->modulation != 0.0) {
		src = AllocDist(geom,nL);
		W = AllocDist(geom,nL);
		W1 = AllocDist(geom,nL);	
		
		GenerateSourceBeam(geom,phan,source,nL,src);
		phan->mua = phan->mua + (source->modulation/(phan->c/phan->n))*I;
	
		AddDist(geom,nL,src,out2);
		AddDist(geom,nL,src,W);
		
		// Generate the Ylm matrix to speed up the PropAbs function
		printf("Computing the time-dependent solution...\n");
		printf("Generating terms for propagator...\n");
		terms2 = GenerateTerms2(geom,phan,nL);
		
		/*printf("Neuman series for %d terms...\n",nTerms);
		for (n=0;n<nTerms-1;n++) {
			PropAbs2(geom,terms2,nL,W,W1);
			PropScat(geom,phan,nL,W1,W);
			AddDist(geom,nL,W,out2);
		}*/

		printf("Neuman series iteration\n");
       		fi=1.0;
        	n=1;
        	while( !StopIteration(geom,&fi,n,W,out) ){
			if( phan->g ){
                        	PropAbs2(geom,terms2,nL,W,W1);
                        	PropScat(geom,phan,nL,W1,W);
                	} else{
                        	PropAbs2(geom,terms2,0,W,W1);
                        	PropScat(geom,phan,0,W1,W);
                	}
                	AddDist(geom,nL,W,out);
                	n++;
        	}
		// Final x-ray transform
		PropAbs2(geom,terms2,nL,out2,W);
        	printf("Neumann series: %d\n",n);

		CopyDist(geom,nL,W,out2);
		FreeDist(W);
		FreeDist(W1);
		FreeDist(src);
		FreeTerms2(terms2);	
	} else {
		CopyDist(geom,nL,out,out2);
	}

	// Must be done before CombineOuts
	phaseImg = ComputePhaseImage(geom,phan,source,nL,101,out2);
	OutputPhaseImage(geom,phaseImg,phaseName);
	
        if (source->modulation != 0.0)
                CombineOuts(geom,nL,out,out2);
        else
                CopyDist(geom,nL,out,out2);

//	OutputDist(geom,out2,distName);
	printf("Generating output...\n");
	transImg = ComputeTransImage(geom,phan,source,nL,101,out2);
	OutputImage(geom,transImg,bottomName);
	refImg = ComputeReflImage(geom,phan,nL,101,out2);
	OutputImage(geom,refImg,topName);		
	ComputeSideFaceImages(geom,phan,nL,101,out2);		// Output side face images
	printf("Number of photons total: %f\n",creal(ComputePhotonNumber(geom,phan,source,nL,101,out2)));
	printf("Neumann series: %d\n", n);

	// Free memory 
	free(phaseImg);
	free(transImg);
	free(refImg);
	FreeDist(out);
	FreeDist(out2);
	FreeSource(source);
	FreePhantom(phan);
	FreeGeometry(geom);
	iniparser_freedict(ini);
	return(0);
}
