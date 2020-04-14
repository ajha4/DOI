//#include <stdio.h>
//#include <stdlib.h>
//#include <errno.h>
//#include <math.h>
//#include <cuda.h>

#include "rte.h"
//#include "iniparser.h"

extern Geometry *geom;
extern Phantom *phan;
extern Source *beam_src;
extern SHORT nL;

//#pragma warning(disable:981)

/*
LoadPhantom

Load up the information about the phantom
*/
Phantom *LoadPhantom(dictionary *ini, int flag) { // if flag = 1, read correct tissue map   
	
	int i,j,k,r_ind,r_ind_phan;
	char tiss_abs_str[50];
	char tiss_sc_str[50];
	if ((phan = (Phantom *)malloc(sizeof(Phantom)*1))== NULL) {
		printf("Error allocation space for phantom\n");
		exit(1);
	}

	/* temporary for now */	
//	phan->mua = iniparser_getdouble(ini,"Phantom:mua",-1) + 0*I;
//	phan->mus = iniparser_getdouble(ini,"Phantom:mus",-1) + 0*I;

	phan->g = iniparser_getdouble(ini,"Phantom:g",-1);
	phan->n = iniparser_getdouble(ini,"Phantom:n",-1);


	char *filename_corr_tissmap;
  

    phan->no_tiss = iniparser_getint(ini, "Phantom:NoTiss", 1) + 1 ; // +1 for the outside boundary tissue type.
	phan->mu_abs = (flt_doub *) malloc (sizeof(flt_doub)*phan->no_tiss);
	phan->mu_sc = (flt_doub *) malloc (sizeof(flt_doub)*phan->no_tiss);


//	tiss_abs_str = (char *)(malloc (sizeof(char)*16));
//	tiss_sc_str = (char *)(malloc (sizeof(char)*16));

	printf("Reading sc and abs coeff \n");
	
	for (i=1; i < phan->no_tiss; i++){
		sprintf(tiss_abs_str,"Phantom:TissAbs%d",i);
		sprintf(tiss_sc_str,"Phantom:TissSc%d",i);
		phan->mu_abs[i] = iniparser_getdouble(ini, tiss_abs_str, 0.01);
		phan->mu_sc[i] = iniparser_getdouble(ini, tiss_sc_str, 1.0);
	}

	phan->mu_abs[0] = 0.0;
    phan->mu_sc[0] = 0.0;
  
	printf("Done Reading sc and abs coeff \n");

    byte *tiss_map;
	FILE *fp_tissmap;
/*

    phan->mu_abs[0] = 0.0;
    phan->mu_sc[0] = 0.0; 
    // Hard coding stuff here
    for (i=1; i < phan->no_tiss; i++){
		phan->mu_abs[i] = 0.03;
		phan->mu_sc[i] = 1.5;
    }

	phan->mu_abs[2] = 0.1;
	phan->mu_sc[2] = 1.5;

*/

	phan->tiss_type= (byte *)malloc(sizeof(byte)*geom->no_vox);
    memset(phan->tiss_type, 0, geom->no_vox);
 
//    fp_tissmap = fopen(filename_tissmap, "r");

//    fread(tiss_map, sizeof(int), geom->nX * geom->nY * geom->nZ, fp_tissmap);

/*
      for(i=0; i < geom->nZ + 2*geom->bounZ; i++){
      	for(j=0; j < geom->nY + 2*geom->bounY; j++){
		      for(k=0; k < geom->nX + 2*geom->bounX; k++){
				  r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
                  phan->tiss_type[r_ind] = 0;
      }}}
*/
    
	tiss_map = (byte *)malloc(sizeof(byte)*geom->nX*geom->nY*geom->nZ);
	if(flag == 1){				
		filename_corr_tissmap = iniparser_getstring(ini, "Phantom:TissMapFile",NULL);
		if(filename_corr_tissmap == NULL){
			printf("Please enter the filename which contains the correct tissue map \n");
			exit(0);
		}
		fp_tissmap = fopen(filename_corr_tissmap, "r");
        if(fp_tissmap == NULL){
           printf("Could not open file %s \n", filename_corr_tissmap);
           exit(0);
        }
		fread(tiss_map, sizeof(byte), geom->nX*geom->nY*geom->nZ, fp_tissmap); 
	 }
	
      for (i=geom->bounZ; i<geom->nZ + geom->bounZ; i++){
        for (j=geom->bounY; j<geom->nY + geom->bounY; j++){
           for (k=geom->bounX; k<geom->nX+geom->bounX; k++){
			    r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
				if(flag == 1){
					r_ind_phan = (i - geom->bounZ)* geom->nX* geom->nY + (j - geom->bounY)* geom->nX + (k-geom->bounZ);
				    phan->tiss_type[r_ind] =  tiss_map[r_ind_phan];
					if(phan->tiss_type[r_ind] != 1 && j== geom->nY/2 && k==geom->nX/2)
						printf("Tissue type %d read for r_ind as (%d %d %d) voxel with r_ind = %d \n", phan->tiss_type[r_ind], i,j,k,r_ind); 
				}
				else
				    phan->tiss_type[r_ind] = 1;
				
             }
          }
      }

       complex_double mut_sum = 0 + 0*I;
       for (i=geom->bounZ; i<geom->nZ + geom->bounZ; i++){
        for (j=geom->bounY; j<geom->nY + geom->bounY; j++){
           for (k=geom->bounX; k<geom->nX+geom->bounX; k++){
                r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
                mut_sum = mut_sum + phan->mu_abs[phan->tiss_type[r_ind]] + phan->mu_sc[phan->tiss_type[r_ind]];
            }
         }
       }

     phan->mut_mean = mut_sum/(geom->nX*geom->nY*geom->nZ);

 
                  
	  
//	free(tiss_abs_str);
//	free(tiss_sc_str);
	if(flag == 1){
		free(tiss_map);
	}

	printf("Done reading phantom information \n");
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
	
	ret->pos_x = iniparser_getdouble(ini,"Source:posX",0.0);
	ret->pos_y = iniparser_getdouble(ini,"Source:posY",0.0);
	ret->pos_z = iniparser_getdouble(ini,"Source:posZ",0.0);
	ret->diam = iniparser_getdouble(ini,"Source:diam",0.2);
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
	
	geom->x_min = iniparser_getdouble(ini,"Geometry:xMin",-1.0);
	geom->y_min = iniparser_getdouble(ini,"Geometry:yMin",-1.0);
	geom->z_min = iniparser_getdouble(ini,"Geometry:zMin",-1.0);
	geom->x_max = iniparser_getdouble(ini,"Geometry:xMax",-1.0);
	geom->y_max = iniparser_getdouble(ini,"Geometry:yMax",-1.0);
	geom->z_max = iniparser_getdouble(ini,"Geometry:zMax",-1.0);
	
	geom->delX = (geom->x_max - geom->x_min) / ((flt_doub)geom->nX);
	geom->delY = (geom->y_max - geom->y_min) / ((flt_doub)geom->nY);
	geom->delZ = (geom->z_max - geom->z_min) / ((flt_doub)geom->nZ);
		
	geom->sub_vox = iniparser_getint(ini,"Geometry:subVox",1);
	geom->sub_thresh = iniparser_getdouble(ini,"Geometry:subThresh",0.0);
	geom->prop_thresh = iniparser_getdouble(ini,"Geometry:propThresh",0.0);
	geom->correction_fact = iniparser_getdouble(ini,"Geometry:correctionFactor",1.0);
	geom->self_sub_vox = iniparser_getdouble(ini,"Geometry:selfsubVox",1.0);

	if(geom->sub_thresh > geom->prop_thresh){
		printf("Invalid subthreshhold %f, as it is greater than prop threshold. Setting it equal to prop threshold \n", geom->prop_thresh);
	    geom->sub_thresh = geom->prop_thresh;
	}

    geom->bounX = 0; //ceil(geom->prop_thresh/geom->delX);
    geom->bounY = 0;//ceil(geom->prop_thresh/geom->delY);
    geom->bounZ = 0;//ceil(geom->prop_thresh/geom->delZ);
    geom->subbounX = 0;//ceilf ((geom->sub_thresh)/(geom->delX));
    geom->subbounY = 0;//ceilf ((geom->sub_thresh)/(geom->delY));
    geom->subbounZ = 0;//ceilf ((geom->sub_thresh)/(geom->delZ));

	geom->no_vox = (geom->nX + 2*geom->bounX)*(geom->nY + 2*geom->bounY)* (geom->nZ + 2*geom->bounZ);
    	
	if ((geom->xs = (flt_doub *)malloc(sizeof(flt_doub)*(geom->nX + 2*geom->bounX))) == NULL) {
		printf("Lack of memory!\n");
	}
	if ((geom->ys = (flt_doub *)malloc(sizeof(flt_doub)*(geom->nY + 2*geom->bounY))) == NULL) {
		printf("Lack of memory!\n");
	}
	if ((geom->zs = (flt_doub *)malloc(sizeof(flt_doub)*(geom->nZ + 2*geom->bounZ))) == NULL) {
		printf("Lack of memory!\n");
	}
	
		
	//printf("%f %f %f\n",geom->delX,geom->delY,geom->delZ);
		
	// Set up the geometry
	// We don't use a + delX in the for loop to mitigate floating-point errors
	for (i = 0;i < geom->nX + 2*geom->bounX; i++) {
		geom->xs[i] = (flt_doub) (i-geom->bounX) * (geom->delX) + geom->x_min + geom->delX/2.0;
	}
	for (i = 0;i < geom->nY + 2*geom->bounY; i++) {
		geom->ys[i] = (flt_doub) (i-geom->bounY) * (geom->delY) + geom->y_min + geom->delY/2.0;
	}
	
	for (i = 0;i < geom->nZ + 2*geom->bounZ; i++) {
		geom->zs[i] = (flt_doub) (i-geom->bounZ) * (geom->delZ) + geom->z_min + geom->delZ/2.0;
	}
	
	return(geom);
}



/*
FreeGeometry

Free up the memory used by a Geometry */
void FreeGeometry(Geometry *geom) {
	free(geom->xs);
	free(geom->ys);
	free(geom->zs);
	free(geom);
}
