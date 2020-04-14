#include <stdio.h>
#include <stdlib.h>
#include <errno.h>


#include "iniparser.h"
#include "RTE.h"

#pragma warning(disable:981)

/*
LoadPhantom

Load up the information about the uniform phantom
*/
Phantom *LoadPhantom(dictionary *ini) {
	Phantom *phan;
	
	if ((phan = (Phantom *)malloc(sizeof(Phantom)*1))== NULL) {
		printf("Error allocation space for phantom\n");
		exit(1);
	}
	
	phan->mua = iniparser_getdouble(ini,"Phantom:mua",-1) + 0*I;
	phan->mus = iniparser_getdouble(ini,"Phantom:mus",-1) + 0*I;
	phan->g = iniparser_getdouble(ini,"Phantom:g",-1);
	phan->n = iniparser_getdouble(ini,"Phantom:n",-1);
	phan->c = iniparser_getdouble(ini,"Phantom:c",2.9979e10);
	
	return(phan);
}



/*
FreePhantom

Free up the storage
*/
void FreePhantom(Phantom *phan) {
	free(phan);
}


Source *LoadSource(dictionary *ini) {
	Source *ret;
	
	if ((ret = (Source *)malloc(sizeof(Source)*1))== NULL) {
		printf("Error allocation space for phantom\n");
		exit(1);
	}
	
	ret->posX = iniparser_getdouble(ini,"Source:posX",0.0);
	ret->posY = iniparser_getdouble(ini,"Source:posY",0.0);
	ret->posZ = iniparser_getdouble(ini,"Source:posZ",0.0);
	ret->diam = iniparser_getdouble(ini,"Source:diam",0.2) + 1.0e-16;
	ret->theta = iniparser_getdouble(ini,"Source:theta",0.0);
	ret->phi = iniparser_getdouble(ini,"Source:phi",0.0);
	ret->mag = iniparser_getdouble(ini,"Source:magnitude",1.0);
	ret->modulation = iniparser_getdouble(ini,"Source:modulation",0.0);
	return(ret);
	
}


void FreeSource(Source *src) {
	free(src);
}


/*
LoadGeometry

Load in the information about the phantom geometry
*/
Geometry *LoadGeometry(dictionary *ini) {
	Geometry *geom;
	int i;
	
	if ((geom = (Geometry *)malloc(sizeof(Geometry)*1)) == NULL) {
		printf("Error allocation space for geometry.\n");
		exit(1);
	}
	
	geom->nX = iniparser_getint(ini,"Geometry:nX",-1);
	geom->nY = iniparser_getint(ini,"Geometry:nY",-1);
	geom->nZ = iniparser_getint(ini,"Geometry:nZ",-1);
	
	geom->xMin = iniparser_getdouble(ini,"Geometry:xMin",-1.0);
	geom->yMin = iniparser_getdouble(ini,"Geometry:yMin",-1.0);
	geom->zMin = iniparser_getdouble(ini,"Geometry:zMin",-1.0);
	geom->xMax = iniparser_getdouble(ini,"Geometry:xMax",-1.0);
	geom->yMax = iniparser_getdouble(ini,"Geometry:yMax",-1.0);
	geom->zMax = iniparser_getdouble(ini,"Geometry:zMax",-1.0);
		
	if ((geom->xs = (double *)malloc(sizeof(double)*geom->nX)) == NULL) {
		printf("Lack of memory!\n");
	}
	if ((geom->ys = (double *)malloc(sizeof(double)*geom->nY)) == NULL) {
		printf("Lack of memory!\n");
	}
	if ((geom->zs = (double *)malloc(sizeof(double)*geom->nZ)) == NULL) {
		printf("Lack of memory!\n");
	}
	
	geom->delX = (geom->xMax - geom->xMin) / ((double)geom->nX);
	geom->delY = (geom->yMax - geom->yMin) / ((double)geom->nY);
	geom->delZ = (geom->zMax - geom->zMin) / ((double)geom->nZ);
		
	printf("%f %f %f\n",geom->delX,geom->delY,geom->delZ);
		
	// Set up the geometry
	// We don't use a + delX in the for loop to mitigate floating-point errors
	for (i = 0;i < geom->nX; i++) {
		geom->xs[i] = (double)i * (geom->delX) + geom->xMin + geom->delX/2.0;
	}
	for (i = 0;i < geom->nY; i++) {
		geom->ys[i] = (double)i * (geom->delY) + geom->yMin + geom->delY/2.0;
	}
	for (i = 0;i < geom->nZ; i++) {
		geom->zs[i] = (double)i * (geom->delZ) + geom->zMin + geom->delZ/2.0;
	}
	
	geom->subVox = iniparser_getint(ini,"Geometry:subVox",1);
	geom->subR = iniparser_getint(ini,"Geometry:subR",1);
	geom->subTheta = iniparser_getint(ini,"Geometry:subTheta",1);
	geom->subPhi = iniparser_getint(ini,"Geometry:subPhi",1);
	geom->self_sub_vox = iniparser_getint(ini,"Geometry:selfSubVox",1);
	
	geom->subThresh = iniparser_getdouble(ini,"Geometry:subThresh",0.0);
	geom->propThresh = iniparser_getdouble(ini,"Geometry:propThresh",0.0);
	
	return(geom);
}



/*
FreeGeometry

Free up the memory used by a Geometry
*/
void FreeGeometry(Geometry *geom) {
	free(geom->xs);
	free(geom->ys);
	free(geom->zs);
	free(geom);
}
