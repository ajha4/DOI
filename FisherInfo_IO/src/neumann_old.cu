//#include <cuda.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <errno.h>
//#include <math.h>
#include "rte.h"
//#include <pthread.h>
//#include <cutil.h>
//#include <cuda_runtime.h>
//#include <cutil_inline.h>

extern Geometry *geom;
extern Phantom *phan;
extern Source *beam_src;
extern complex_double *diag_terms_host;
extern complex_double *sph_harm;
extern Info_Stat *info_stat_host;
extern SHORT nL;
extern int nTerms;



void Neumann(complex_double* src_host, complex_double *out_host, int flag) // if flag = 1 ,run the entire Neumann series, else just the final absorption term

{ 

    unsigned int timer;
    timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    int i, j,k,r_ind,cnt;

	int n;
	int thread_index, tid;
	int num_layers_gpu[NUM_DEVICES];
	pthread_t thread_data[NUM_DEVICES];

	THREAD_PARAMETERS thread_parameters[NUM_DEVICES];
	int size_layer = ( geom->nX + 2*geom->bounX ) * ( geom->nY + 2*geom->bounY ) * ( nL+1) * (nL+1);

    int size = size_layer * (geom->nZ + 2*geom->bounZ);

    int num_layers_per_gpu = (int) floorf(geom->nZ / ( NUM_DEVICES));
    int rem = geom->nZ % (NUM_DEVICES);


    for (thread_index = 0; thread_index < NUM_DEVICES; thread_index++){
    	num_layers_gpu[thread_index] = num_layers_per_gpu;
		if(rem > thread_index){
			num_layers_gpu[thread_index] += 1;
		}
	}		
	complex_double *W_out_host;

    W_out_host = (complex_double *) malloc ( sizeof(complex_double)*size);

    n =2;
    double fi=1.0;

    if(flag == 1){
    while(!StopIteration(&fi,n,src_host,out_host)){
       n++;
//      for (n=0;n<nTerms-2;n++) {
//    printf("%d term of Neumann series \n",n);
		memset(W_out_host, 0, sizeof(complex_double)*size);

#if 0
    printf("Before the prop_abs function \n");
     for (i=geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
      for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
        for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
           for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt = cnt +2) {
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
         //   if((W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real() != (W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real()){
                if((j==geom->nY/2 && k == geom->nX/2 && cnt==0)|| src_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real() > 1000.0 || W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real() > 1000.0){
                 printf("%e % e  %e , %d (%d %d %d) %d \n", out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), src_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(),cnt, i-geom->bounX,j-geom->bounY,k-geom->bounZ, VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1)));
   }}}}}
#endif

    for(thread_index = 0; thread_index < NUM_DEVICES; ++thread_index) {
        thread_parameters[thread_index].device_index = thread_index + MY_START_DEVICE;
        thread_parameters[thread_index].src_host = src_host;
            thread_parameters[thread_index].num_layers = num_layers_gpu[thread_index];
            thread_parameters[thread_index].layer_start = 0 ;
            for(tid = 0; tid < thread_index; tid++){
                thread_parameters[thread_index].layer_start += num_layers_gpu[tid];
            }
        thread_parameters[thread_index].out_host = W_out_host + ( thread_parameters[thread_index].layer_start + geom->bounZ) * size_layer ;
        if(phan->g != 0.0){
        //    printf("Medium is anisotropic \n");
 			thread_parameters[thread_index].flag = 1;
		}
		else
	        thread_parameters[thread_index].flag = 0 ;
		thread_parameters[thread_index].flag_scat = 1 ;
        pthread_create(& thread_data[thread_index], NULL, prop_abs, &thread_parameters[thread_index]);
    }

		for(thread_index = 0; thread_index < NUM_DEVICES; thread_index++) {
	    	pthread_join(thread_data[thread_index], NULL);
		}

#if 0
	printf("Before the prop_scat function \n");
	 for (i=geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
      for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
        for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
           for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt = cnt +2) {
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
         //   if((W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real() != (W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real()){
                if((j==geom->nY/2 -1 && k == geom->nX/2 && cnt==0)){//|| fabs(src_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].imag()) > 0 || fabs(W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].imag()) > 0.0){
                 printf("%e % e  %e , %d (%d %d %d) %d \n", out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].imag(),cnt, i-geom->bounX,j-geom->bounY,k-geom->bounZ, VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1)));
   }}}}}
#endif
       
	//	prop_scat(W_out_host, src_host);
		 // copy_dist(W_out_host,src_host);
#if 0
		printf("After the prop_scat function \n");
   for (i=geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
      for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
        for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
           for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt = cnt +2) {
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
         //   if((W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real() != (W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real()){
				if((j==geom->nY/2 && k == geom->nX/2 && cnt==0)|| src_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real() > 1000.0 || W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real() > 1000.0){
                 printf("%e % e  %e , %d (%d %d %d) %d \n", out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), src_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(),cnt, i-geom->bounX,j-geom->bounY,k-geom->bounZ, VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1)));
   }}}}}
#endif

        //add_dist( src_host, out_host, out_host);
        add_dist( W_out_host, out_host, out_host);
		copy_dist(W_out_host,src_host);
#if 0
       printf("\n \n \n \n \n \n \n After the add_dist function \n");
   for (i=0;i<geom->nZ + 2*geom->bounZ;i++) {
      for (j=0;j<geom->nY + 2*geom->bounY;j++) {
        for (k=0;k<geom->nX + 2*geom->bounX;k++) {
           for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt = cnt +2) {
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
//            if((out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real() != (out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real()){
            if((j==geom->nY/2 && k == geom->nX/2 && cnt==0) || out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real() > 1000.0){
                 printf("%e % e  %e , %d (%d %d %d) %d \n", out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), src_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(),cnt, i-geom->bounX,j-geom->bounY,k-geom->bounZ, VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1)));
   }}}}}
#endif
    //n++;
    //printf("Neumann series for %d terms complete \n", n);
    }

#if 0
   printf("\n \n \n \n \n \n \n Before the final prop_abs function \n");
   for (i=0;i<geom->nZ + 2*geom->bounZ;i++) {
      for (j=0;j<geom->nY + 2*geom->bounY;j++) {
        for (k=0;k<geom->nX + 2*geom->bounX;k++) {
           for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt = cnt +2) {
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
//            if((out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real() != (out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real()){
            if((j==geom->nY/2 && k == geom->nX/2 && cnt==0) || out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real() > 1000.0){
                 printf("%e % e  %e , %d (%d %d %d) %d \n", out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), src_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(),cnt, i-geom->bounX,j-geom->bounY,k-geom->bounZ, VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1)));
   }}}}}
#endif
   }
    //printf("Calling the absorption kernel for the last time\n");
	memset(W_out_host, 0, sizeof(complex_double)*size);
    cutilCheckError(cutResetTimer(timer));
    cutilCheckError(cutStartTimer(timer));
 
	for(thread_index = 0; thread_index < NUM_DEVICES; ++thread_index) {
		thread_parameters[thread_index].device_index = thread_index + MY_START_DEVICE;
		thread_parameters[thread_index].src_host = out_host;
			thread_parameters[thread_index].num_layers = num_layers_gpu[thread_index];
			thread_parameters[thread_index].layer_start = 0 ;
			for(tid = 0; tid < thread_index; tid++){
				thread_parameters[thread_index].layer_start += num_layers_gpu[tid];
			}
		thread_parameters[thread_index].out_host = W_out_host + (thread_index*num_layers_per_gpu + geom->bounZ) * size_layer ;
			thread_parameters[thread_index].flag = 1 ;
			thread_parameters[thread_index].flag_scat = 0 ;
		pthread_create(& thread_data[thread_index], NULL, prop_abs, &thread_parameters[thread_index]);
	}
  
	for(thread_index = 0; thread_index < NUM_DEVICES; ++thread_index) {
    	pthread_join(thread_data[thread_index], NULL);
	}
#if 0
   printf("\n \n \n \n \n \n \n After the final prop_abs function \n");
   for (i=0;i<geom->nZ + 2*geom->bounZ;i++) {
      for (j=0;j<geom->nY + 2*geom->bounY;j++) {
        for (k=0;k<geom->nX + 2*geom->bounX;k++) {
           for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt = cnt +2) {
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
//            if((out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real() != (out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))]).real()){
            if(j==geom->nY/2 && k == geom->nX/2 && cnt==0){
                 printf("%e % e  %e , %d (%d %d %d) %d \n", out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), W_out_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(), src_host[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))].real(),cnt, i-geom->bounX,j-geom->bounY,k-geom->bounZ, VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1)));
   }}}}}
#endif
//    printf("Time taken for the absorption kernel :%f ms \n", cutGetTimerValue(timer));	
    copy_dist(W_out_host, out_host);
   
}
