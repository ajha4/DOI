#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "rte.h"
#include <pthread.h>


__host__ Info_Stat * populate_info_dev();
Geometry *geom;
Phantom *phan;
Source *beam_src;
complex_double *diag_terms_host;
complex_double *sph_harm;
Info_Stat *info_stat_host;
SHORT nL;
int nTerms;

__host__ int get_vind_phanind_host(int dep, int row, int col){

  return ((geom->bounZ + dep) * (geom->nX + 2 * geom->bounX ) * (geom->nY + 2 * geom->bounY ) /* reached the correct layer */ + ( geom->bounY + row)* (geom->nX + 2 * geom->bounX ) + (geom->bounX + col));

}


__host__ Info_Stat * populate_info_dev(){

    Info_Stat *info_stat_host;
    info_stat_host = (Info_Stat *) malloc (sizeof(Info_Stat));

    info_stat_host->nX = geom->nX;
    info_stat_host->nY = geom->nY;
    info_stat_host->nZ= geom->nZ;

    info_stat_host->bounX = geom->bounX;
    info_stat_host->bounY = geom->bounY;
    info_stat_host->bounZ= geom->bounZ;

    info_stat_host->subbounX = ceilf ((geom->sub_thresh)/(geom->delX));
    info_stat_host->subbounY = ceilf ((geom->sub_thresh)/(geom->delY));
    info_stat_host->subbounZ = ceilf ((geom->sub_thresh)/(geom->delZ));

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

    int i;
    for(i=0; i < phan->no_tiss; i++){
        info_stat_host->mu_tot[i] = phan->mu_abs[i] + phan->mu_sc[i];
        info_stat_host->mu_sc[i] = phan->mu_sc[i];
    }

    return info_stat_host;
}


int main(int argc, char** argv )
{
		
	time_t start_time, end_time;
	time(&start_time);
    dictionary *ini;

    complex_double *transImg;
    complex_double *gbar,*dgbar;
    complex_double *W, *W1, *out, *out2, *out3, *src, *tmp1, *df;
    

    double *g;
    double *grada,*grads;
    complex_double delg;

    char *bottomName;
    char *topName;
    char *phaseName;
    double *phaseImg;
    char numStr[100];
    char tmpStr[200];
    int n;
    FILE *gFile,*fid,*fid1,*resFID;
    int itt;
    int j,jj;
    int jnk;
    int tiss_idx;
    int numIT;
    int iip, jjp, kkp;
    int r_ind,i;
    int size;

        if (argc != 2) {
        printf("\n InverseRTE file.par\n");
        printf("     file.par is the parameter file.\n\n");
        exit(1);
    }

    // Load in the initialization file
    ini = iniparser_load(argv[1]);

    // Set up the geometry, phantom, etc
    printf("Loading in geometry and phantom information...\n");
    geom = LoadGeometry(ini);
    phan = LoadPhantom(ini,1);

    beam_src = LoadSource(ini);

	printf("Done reading source information \n");

    nL = iniparser_getint(ini,"Algorithm:nL",-1);
    nTerms = iniparser_getint(ini,"Algorithm:nTerms",1);
    bottomName = iniparser_getstring(ini, "Runtime:TransFile", NULL);
    topName = iniparser_getstring(ini, "Runtime:ReflectFile", NULL);
    phaseName = iniparser_getstring(ini,"Runtime:PhaseFile",NULL);

    numIT = iniparser_getint(ini,"Runtime:numIterations",1);

    double stepsizea, stepsizes;
    stepsizea = iniparser_getdouble(ini,"Runtime:stepsizea",1e3);
    stepsizes = iniparser_getdouble(ini,"Runtime:stepsizes",1e3);

	printf("Done reading integers \n");

	int jnk2;
    if ( (gFile = fopen(iniparser_getstring(ini,"Runtime:gFile",NULL),"r") ) == NULL){
		printf("Error in opening gfile. Exiting \n");
		exit(0);
	}
	printf("%s is gFile \n", gFile);
    fread(&jnk,sizeof(int),1,gFile);
    fread(&jnk2,sizeof(int),1,gFile);
	printf("Done reading gfile integers %d and %d \n", jnk, jnk2);
    g = (flt_doub *) malloc(sizeof(flt_doub) * geom->nX * geom->nY);
    fread(g,sizeof(flt_doub),geom->nX * geom->nY,gFile);
    fclose(gFile);

	printf("Done reading gfile \n");

#if 0
    for (i=0;i<geom->nY;i++) {
        for (j=0;j<geom->nX;j++) {
            if((g[i*geom->nX+j])!=0.0){
                   printf("%f (%d %d) \n", g[i*geom->nX+j], i,j);
	}}}
#endif

    info_stat_host = populate_info_dev();

	size = (nL+1)*(nL+1)* geom->no_vox;


    printf("Generating the spherical harmonic terms \n");
    generate_sph_harm_terms();

    unsigned int timer;
	int cnt,k;
	double tmp;
	double thresh_abs, thresh_sc;
	thresh_abs = 0.0000001;
	thresh_sc =  0.0000001;
     

	byte *flag_grada, *flag_grads;

	flag_grada = (byte*) malloc ( sizeof(byte)* phan->no_tiss);
	flag_grads = (byte*) malloc ( sizeof(byte)* phan->no_tiss);

	memset(flag_grada, 0, phan->no_tiss);
	memset(flag_grads, 0, phan->no_tiss);
 
	byte flag_net; 
	int r_ind_phan;
    resFID = fopen("Residual.out","w");
        for (itt = 0;itt < numIT;itt++ ) {

			flag_net = 1;

			for(tiss_idx = 1; tiss_idx < phan->no_tiss; tiss_idx++){
				if ( flag_grads[tiss_idx] == 0 || flag_grada[tiss_idx] == 0){
					printf("%d tissue type can still change \n", tiss_idx);
					flag_net = 0;
					break;
				}
			}
    
			if(flag_net == 1){
               printf("Time to terminate the iterations at iteration number %d \n", itt);
			   break;
			}
				
		    generate_diag_terms_host();

            sprintf(numStr,"%4d",itt);
            j = 0;
            while (numStr[j] != '\0') {
                if (numStr[j] == ' ') numStr[j] = '0';
                j++;
            }

            out2 = alloc_dist();

            src = alloc_dist();
            W = alloc_dist();

            generate_source_beam(src);

            copy_dist(src,out2);
            copy_dist(src,W);

            printf("Computing gbar\n");
            Neumann(W,out2);

            free(W);
            free(src);


            // THE GRADIENT of mua
            //___________________________________   
            gbar = generate_trans_image(out2,1);
//          OutputImage(geom,gbar,"Hi.dat");
  

#if 0 
			tmp = 0; 
			for (i=0;i<geom->nY;i++) {
		        for (j=0;j<geom->nX;j++) {
        		    if((fabs(g[i*geom->nX+j] - gbar[i*geom->nX+j].real())) > 0.000001 ){
                  		//printf("%e %e (%d %d) \n", g[i*geom->nX+j], gbar[i*geom->nX+j].real(), i,j);
					}
					tmp = tmp + g[i*geom->nX+j] - gbar[i*geom->nX+j].real();
					}}
			printf("tmp is %f \n", tmp);
#endif

            out3 = alloc_dist();

            grada = (double *)malloc(sizeof(double)*phan->no_tiss);
            printf("Compute the gradient with respect to mua");
            for (tiss_idx = 1; tiss_idx < phan->no_tiss; tiss_idx++) {

				if(flag_grada[tiss_idx] == 0){

	               copy_dist(out2,out3);

				   for (iip=0; iip<geom->nZ; iip++) {
				      for (jjp=0; jjp<geom->nY; jjp++) {
				        for (kkp=0; kkp<geom->nX; kkp++) {

	            	    r_ind = (iip + geom->bounZ)* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + (jjp + geom->bounY)* (geom->nX + 2*geom->bounX) + (kkp + geom->bounX);
						r_ind_phan =iip* geom->nX * geom->nY + jjp * geom->nX + kkp;

                    	if ( phan->tiss_type[r_ind_phan] != tiss_idx) {
						
                        	for (n=0;n<(nL+1)*(nL+1);n++) {
                            	out3[VOX_TO_SPIND(r_ind, n,(nL+1)*(nL+1))] = 0+0*I;
	                        }
    	                }
        	        }}}

            	    scale_dist(out3,-1.0*C);
                	src = alloc_dist();
	                copy_dist(out3,src);
    	            Neumann(src,out3);
        	        free(src);
            	    dgbar = generate_trans_image(out3,0);
#if 0
            		for (i=0;i<geom->nY;i++) {
		              for (j=0;j<geom->nX;j++) {
        	           if((fabs(dgbar[i*geom->nX+j].real())) > 0.00000001 && i == 4 && j == 4 && tiss_idx%(geom->nX*geom->nY) == 44){
                        //printf("g = %f, gbar =  %f, dgbar =  %f (%d %d) \n",g[j], gbar[j].real(), dgbar[i*geom->nX+j].real(), i,j);
	                    }
                    }}
#endif
                	grada[tiss_idx] = 0.0;
	                for (j = 0;j<geom->nX*geom->nY;j++) {
        	            grada[tiss_idx] = grada[tiss_idx] - 2.0*(g[j]-gbar[j].real())*dgbar[j].real();
						if ( j == 4*geom->nX + 4 && tiss_idx%(geom->nX*geom->nY) == 4*geom->nX + 4 ){
						//	printf("g = %f, gbar =  %f, dgbar =  %f \n", g[j], gbar[j].real(), dgbar[j].real()); 
						}
                	}
					//if(fabs(grada[tiss_idx]) > 0.0){
		                printf("\n mu_a tiss type:  %d of %d.  Value of grada = %e                   \n",tiss_idx,phan->no_tiss,grada[tiss_idx]);
					//}
    	            free(dgbar);
        	    }
			}
/*
            sprintf(tmpStr,"GradA%s.out",numStr);
            fid = fopen(tmpStr,"w");
            fwrite(grada,sizeof(double),phan->no_tiss,fid);
            fclose(fid);
*/
           	printf("\n");
	        free(out3);

            printf("Compute the gradient with respect to mus \n");
            grads = (double *)malloc(sizeof(double)* phan->no_tiss);
#if 1
            // THE GRADIENT of mus
            //___________________________________   


            out3 = alloc_dist();
            tmp1 = alloc_dist();

            for (tiss_idx = 1; tiss_idx < phan->no_tiss ; tiss_idx++) {


			   if(flag_grads[tiss_idx] == 0){
	                copy_dist(out2,out3);
	
    	            scale_dist(out3,-1.0*C);

        	        copy_dist(out2,tmp1);
            	    PropScatmu1(geom,phan,nL,tmp1);

    	            add_dist(tmp1,out3,out3);
				   
					for (iip=0; iip<geom->nZ; iip++) {
				      for (jjp=0; jjp<geom->nY; jjp++) {
				        for (kkp=0; kkp<geom->nX; kkp++) {

						r_ind_phan =iip* geom->nX * geom->nY + jjp * geom->nX + kkp;

                    	if ( phan->tiss_type[r_ind_phan] != tiss_idx) {
        		                for (n=0;n<(nL+1)*(nL+1);n++) {
                		            out3[VOX_TO_SPIND(r_ind, n, (nL+1)*(nL+1))] = 0+0*I;
                        		}

                    		}
	                }}}


        	        src = alloc_dist();
            	    copy_dist(out3,src);
                	Neumann(src,out3);
	                free(src);
    	            dgbar = generate_trans_image(out3,0);

        	        grads[tiss_idx] = 0.0;
            	    for (j = 0;j<geom->nX*geom->nY;j++) {
						grads[tiss_idx] = grads[tiss_idx] - 2.0*(g[j]-gbar[j].real())*dgbar[j].real();
                	}

					//if(fabs(grads[tiss_idx]) > 0.0){
		                printf("mu_s tiss_type:  %d of %d.  Value of grads= %e \n",tiss_idx,phan->no_tiss,grads[tiss_idx]);
					//}
    	            free(dgbar);
        	    }
			}

            printf("\n");
            free(out3);
            free(tmp1);
#endif

            free(out2);
		    sprintf(tmpStr,"tissmapabs%s.out",numStr);
            fid = fopen(tmpStr,"w");
		    sprintf(tmpStr,"tissmapsc%s.out",numStr);
            fid1 = fopen(tmpStr,"w");
			printf("Iteration No %d \n", itt);
            for (tiss_idx = 0;tiss_idx< phan->no_tiss;tiss_idx++) {

                if (fabs(grada[tiss_idx]) < thresh_abs)
					 flag_grada[tiss_idx] = 1;
				if ( fabs(grads[tiss_idx]) < thresh_sc)
					flag_grads[tiss_idx] = 1;

				if ( flag_grada[tiss_idx] == 0 ||  flag_grads[tiss_idx] == 0){
				printf("Old value of mua = %e mus = %e for tisstype = %d\n", phan->mu_abs[tiss_idx].real(), phan->mu_sc[tiss_idx].real(),tiss_idx);
					
				if ( flag_grada[tiss_idx] == 0){
	        	   phan->mu_abs[tiss_idx] = phan->mu_abs[tiss_idx] + stepsizea*grada[tiss_idx];
				   if(phan->mu_abs[tiss_idx].real() < 0)
						phan->mu_abs[tiss_idx] = phan->mu_abs[tiss_idx] - 2*stepsizea*grada[tiss_idx];
				}
				if ( flag_grads[tiss_idx] == 0){
    	           phan->mu_sc[tiss_idx] = phan->mu_sc[tiss_idx] + stepsizes*grads[tiss_idx];
				   if(phan->mu_sc[tiss_idx].real() < 0)
						phan->mu_sc[tiss_idx] = phan->mu_sc[tiss_idx] - 2*stepsizes*grads[tiss_idx];
				}
            
				
				printf("Value of mua grad = %e musgrad = %e for tissue type %d \n", grada[tiss_idx], grads[tiss_idx], tiss_idx);
				printf("New value of mua = %e mus = %e for tissue type = %d\n", phan->mu_abs[tiss_idx].real(), phan->mu_sc[tiss_idx].real(),tiss_idx);
				}

				fwrite(&(phan->tiss_type[j]),sizeof(byte),1,fid);
	            fwrite(&(phan->tiss_type[j]),sizeof(byte),1,fid1);

            }
    
			/*updating the absorption and scatter coeff. values in the device code */

            populate_info_dev(); 

   
            fclose(fid);
            fclose(fid1);


            free(grads);
            free(grada);

            delg = 0.0;
            for (i=0;i<geom->nX*geom->nY;i++) {
                delg = delg + (gbar[i]-g[i])*(gbar[i]-g[i]);
            }
            printf("Residual = %e\n",(delg).real());
            fwrite(&delg,sizeof(double),1,resFID);
            free(gbar);




        }

        fclose(resFID);

        free(g);
        free(beam_src);
        free(phan);
        free(geom);
        iniparser_freedict(ini);

		time(&end_time);
	
		printf("\n*------------------------------------------*\n");
		printf("\nThe total time taken by the code = %d sec \n", end_time - start_time);
		printf("\n*------------------------------------------*\n");

        return(0);
}

#if 0

void Neumann(complex_double* src_host, complex_double *out_host){
    

   unsigned int timer;
    timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
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
    memset(W_out_host, 0, sizeof(complex_double)*size);

    for (n=0;n<nTerms-1;n++) {

        for(thread_index = 0; thread_index < NUM_DEVICES; thread_index++) {
            thread_parameters[thread_index].device_index =thread_index;
            thread_parameters[thread_index].src_host = src_host;
            thread_parameters[thread_index].num_layers = num_layers_gpu[thread_index];
            thread_parameters[thread_index].layer_start = 0 ;
            for(tid = 0; tid < thread_index; tid++){
                thread_parameters[thread_index].layer_start += num_layers_gpu[tid];
            }
            thread_parameters[thread_index].out_host = W_out_host + (thread_parameters[thread_index].layer_start + geom->bounZ) * size_layer ;
            pthread_create(& thread_data[thread_index], NULL, prop_abs, &thread_parameters[thread_index]);
        }


        for(thread_index = 0; thread_index < NUM_DEVICES; thread_index++) {
            pthread_join(thread_data[thread_index], NULL);
        }

        prop_scat(W_out_host, src_host);
        add_dist( src_host, out_host, out_host);

    }

    printf("Calling the absorption kernel \n");
    cutilCheckError(cutResetTimer(timer));
    cutilCheckError(cutStartTimer(timer));

    for(thread_index = 0; thread_index < NUM_DEVICES; ++thread_index) {
        thread_parameters[thread_index].device_index = thread_index;
        thread_parameters[thread_index].src_host = out_host;
            thread_parameters[thread_index].num_layers = num_layers_gpu[thread_index];
            thread_parameters[thread_index].layer_start = 0 ;
            for(tid = 0; tid < thread_index; tid++){
                thread_parameters[thread_index].layer_start += num_layers_gpu[tid];
            }
        thread_parameters[thread_index].out_host = W_out_host + (thread_index*num_layers_per_gpu + geom->bounZ) * size_layer ;
        pthread_create(& thread_data[thread_index], NULL, prop_abs, &thread_parameters[thread_index]);
    }

    for(thread_index = 0; thread_index < NUM_DEVICES; ++thread_index) {
        pthread_join(thread_data[thread_index], NULL);
    }

   printf("Time taken for the absorption kernel :%f ms \n", cutGetTimerValue(timer));

    copy_dist(W_out_host, out_host);


	int cnt,i,j,k,r_ind;
#if 0
		   for (cnt=0;cnt< (nL+1)*(nL+1) ;cnt = cnt +2) {
			  for (i=geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
			      for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
			         for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
			            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
            			if((out[VOX_TO_SPIND(r_ind,cnt,geom->no_vox)]).real()!=0.0){
                        printf("% e , %e i, %d (%d %d %d) \n", out[VOX_TO_SPIND(r_ind,cnt,geom->no_vox)].real(), W[VOX_TO_SPIND(r_ind,cnt,geom->no_vox)].real(),cnt, i-geom->bounX,j-geom->bounY,k-geom->bounZ);
		    }}}}}   
#endif     
}
#endif

