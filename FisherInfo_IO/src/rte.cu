//#include<cuda.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <errno.h>
//#include <math.h>
#include "rte.h"
//#include <pthread.h>
//#include <cutil.h>
//#include <cuda_runtime.h>
//#include <cutil_inline.h>



Geometry *geom;
Phantom *phan;
Source *beam_src;
complex_double *diag_terms_host;
complex_double *sph_harm;
Info_Stat *info_stat_host;
SHORT nL;
int nTerms;
__host__ Info_Stat * populate_info_dev();


__host__ Info_Stat * populate_info_dev(){

    Info_Stat *info_stat_host;
    info_stat_host = (Info_Stat *) malloc (sizeof(Info_Stat));

    info_stat_host->nX = geom->nX;
    info_stat_host->nY = geom->nY;
    info_stat_host->nZ= geom->nZ;

    info_stat_host->bounX = geom->bounX;
    info_stat_host->bounY = geom->bounY;
    info_stat_host->bounZ= geom->bounZ;

    info_stat_host->subbounX = geom->subbounX; 
    info_stat_host->subbounY =  geom->subbounY;
    info_stat_host->subbounZ =  geom->subbounZ;

    info_stat_host->delX = geom->delX;
    info_stat_host->delY = geom->delY;
    info_stat_host->delZ= geom->delZ;

    info_stat_host->x_min = geom->x_min;
    info_stat_host->y_min = geom->y_min;
    info_stat_host->z_min = geom->z_min;

    info_stat_host->x_max = geom->x_max;
    info_stat_host->y_max = geom->y_max;
    info_stat_host->z_max = geom->z_max;

    info_stat_host->sub_thresh = geom->sub_thresh;
    info_stat_host->prop_thresh = geom->prop_thresh;
    info_stat_host->sub_vox = geom->sub_vox;
    info_stat_host->self_sub_vox = geom->self_sub_vox;

    info_stat_host->g = phan->g;
    info_stat_host->n = phan->n;
    info_stat_host->no_tiss = phan->no_tiss;

	info_stat_host->cm = C/phan->n;
  
    info_stat_host->no_vox = geom->no_vox;

/*
    info_stat_host->boun_blkX = ceilf(geom->bounX/BLK_SIZE);
    info_stat_host->boun_blkY = ceilf(geom->bounY/BLK_SIZE);
    info_stat_host->boun_blkZ = ceilf(geom->bounZ/BLK_SIZE);

	info_stat_host->boun_blk_sizeZ = geom->boun_blkZ*BLK_SIZE;
	info_stat_host->boun_blk_sizeY = geom->boun_blkY*BLK_SIZE;
	info_stat_host->boun_blk_sizeX = geom->boun_blkX*BLK_SIZE;
*/
    int i;
    for(i=0; i < phan->no_tiss; i++){
        info_stat_host->mu_tot[i] = 0.0;
        info_stat_host->mu_sc[i] = 0.0;
    }
    for(i=0; i < phan->no_tiss; i++){
        info_stat_host->mu_tot[i] = phan->mu_abs[i] + phan->mu_sc[i];
        info_stat_host->mu_sc[i] = phan->mu_sc[i];
    }

    //info_stat->nL = nL;

    return info_stat_host;
}


int main(int argc, char **argv){

	time_t start_time, end_time, start_time_diag, end_time_diag;
    dictionary *ini;
    int size;
	complex_double *W, *src, *out;
    complex_double *trans_img,*ref_img;
    double *phase_img;
	char *bottomName,*topName,*phaseName;

	time(&start_time);
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
    phan= LoadPhantom(ini,1);
    beam_src = LoadSource(ini);

    if(((phan->mut_mean).real()*(geom->z_max - geom->z_min)) > 4){
    	geom->prop_thresh = -log(THRESH_PROP)/(phan->mut_mean.real());
        printf("Diffuse medium with prop_thresh = %f \n", geom->prop_thresh);
     }

    nL = iniparser_getint(ini,"Algorithm:nL",1000);
    if (nL == 1000) {
        printf("nL = 1000. Exiting\n");
        exit(1);
    }
    nTerms = iniparser_getint(ini,"Algorithm:nTerms",1);
    bottomName = iniparser_getstring(ini, "Runtime:TransFile", NULL);
    topName = iniparser_getstring(ini, "Runtime:ReflectFile", NULL);
    phaseName = iniparser_getstring(ini,"Runtime:PhaseFile",NULL);


    info_stat_host = populate_info_dev();


    size = (nL+1)*(nL+1)* geom->no_vox;

	src = alloc_dist();
    out = alloc_dist();
	W = alloc_dist();

    printf("Generating the source beam distribution \n");
    generate_source_beam(src); // Generates the KXE term as the distribution


	int i;
    copy_dist(src,out);
    copy_dist(src,W);
	
#if 0
	printf("KXE terms: \n");
	for(i=0; i<size; i++){
	    if(src[i].real() != 0.0){
			printf("%e %e %d \n", W[i].real(), src[i].real() ,i);
		}
	}
#endif
   
    printf("Generating the spherical harmonic terms \n");
    generate_sph_harm_terms();

 
	printf("Generating the diagonal terms \n");
	time(&start_time);
	generate_diag_terms_host();
    time(&end_time);
    printf("Time taken for generating the diagonal terms :%f ms \n", end_time - start_time);

#if 0
   int j,k,r_ind,cnt;
   for (i= 0;i<geom->nZ + 2*geom->bounZ;i++) {
      for (j=0;j<geom->nY + 2*geom->bounY;j++) {
        for (k=0;k<geom->nX + 2*geom->bounX;k++) {
           for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt = cnt +2) {
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
            if((out[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real() != (out[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real()){
                 printf("%e % e  %e , %d (%d %d %d) %d \n", out[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), src[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(),cnt, i-geom->bounZ,j-geom->bounY,k-geom->bounX, VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1)));
   }}}}}
#endif

	Neumann( W, out,1);

	printf("\n \n \n \n");
#if 0
   /*for (i=geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
      for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
        for (k=geom->bounX;k<geom->nX + geom->bounX;k++) { */
   for (i= geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
      for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
        for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
   		   for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt ++) {
			r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
//            if((out[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real() != 0 && i == geom->nZ/2 && j==geom->nY/2 && k == geom->nX/2){
			  if(j==geom->nY/2 && k == geom->nX/2 && cnt==0){
                 printf("%e % e  %e , %d (%d %d %d) %d \n", out[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), src[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(),cnt, i-geom->bounZ,j-geom->bounY,k-geom->bounX, VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1)));
   }}}}}
#endif



    phase_img = generate_phase_image(out,1);
    output_phase_image(phase_img,phaseName);


    printf("Generating the output image \n");
    trans_img = generate_trans_image( out,1 );
    output_image(trans_img,bottomName);
    ref_img = generate_ref_image( out);
    output_image(ref_img,topName);

#if 0	
	 for(i=0; i< geom->nY; i++){
        for(j=0; j <geom->nX; j++){
			if(trans_img[i*geom->nX+j].real()!= 0.0)
            printf("%e at (%d %d) \n", (trans_img[i*geom->nX+j].real()),i,j);
    }}
#endif
	time(&end_time);
	printf("\n*------------------------------------------*\n");
	printf("\nThe total time taken by the code = %d sec \n", end_time - start_time);
	printf("\n*------------------------------------------*\n");
	
	//printf("Number of photons total: %f\n",(ComputePhotonNumber(101,out)).real());

    free(W);
    free(out);
    free(src);
    free(beam_src);
    free(geom);
    free(phan);
	

	
    return(0);
}

