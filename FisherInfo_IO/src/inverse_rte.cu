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
    complex_double *gbar,*dgbara, *dgbars;
    complex_double *W, *W1, *out, *out2, *out3, *src, *tmp1, *df, *src2, *out4;
    

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
    float av_mua; 
    float av_mus;

    av_mua = 0.2;
    av_mus = 1;

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

	int jnk2;
    if ((gFile = fopen(iniparser_getstring(ini,"Runtime:gFile",NULL),"r")) == NULL){
		printf("Error in opening gfile. Exiting \n");
		exit(0);
	}
	//printf("%s is gFile \n", gFile);
    fread(&jnk,sizeof(int),1,gFile);
    fread(&jnk2,sizeof(int),1,gFile);
	printf("Done reading gfile integers %d and %d \n", jnk, jnk2);
    g = (flt_doub *) malloc(sizeof(flt_doub) * geom->nX * geom->nY);
    fread(g,sizeof(flt_doub),geom->nX * geom->nY,gFile);
    fclose(gFile);


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

    double fisher11, fisher12, fisher21, fisher22;
     
    grada = (double *)malloc(sizeof(double)*phan->no_tiss);
    grads = (double *)malloc(sizeof(double)* phan->no_tiss);

	byte *flag_grada, *flag_grads;

	flag_grada = (byte*) malloc ( sizeof(byte)* phan->no_tiss);
	flag_grads = (byte*) malloc ( sizeof(byte)* phan->no_tiss);

	memset(flag_grada, 0, phan->no_tiss);
	memset(flag_grads, 0, phan->no_tiss);
 
	byte flag_net; 
	int r_ind_phan;
    W = alloc_dist();
    out2 = alloc_dist();
    out3 = alloc_dist();
    out4 = alloc_dist();
    src = alloc_dist();
    src2 = alloc_dist();
    tmp1 = alloc_dist();
    resFID = fopen("Residual.out","w");
    float mus_hat, mua_hat;

    for (itt = 0;itt < numIT;itt++ ) {

	   flag_net = 1;
	   for(tiss_idx = 1; tiss_idx < phan->no_tiss; tiss_idx++){
	 	 if ( flag_grads[tiss_idx] == 0 || flag_grada[tiss_idx] == 0){
			flag_net = 0;
			break;
		 }
	   }
    
   	   if(flag_net == 1){
			   break;
		}
				
	   generate_diag_terms_host();

       sprintf(numStr,"%4d",itt);
       j = 0;
       while (numStr[j] != '\0') {
          if (numStr[j] == ' ') numStr[j] = '0';
             j++;
       }

       generate_source_beam(src);
       copy_dist(src,out2);
       copy_dist(src,W);
       Neumann(W,out2,1);

       // THE GRADIENT of mua
       gbar = generate_trans_image(out2,1);
  
       //printf("Compute the gradient with respect to mua");
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
	             copy_dist(out3,src);
                 Neumann(src,out3,1);
            	 dgbara = generate_trans_image(out3,0);
                 grada[tiss_idx] = 0.0;
                 for (j = 0;j<geom->nX*geom->nY;j++) {
                    if(g[j] ||  gbar[j].real()){
        	            grada[tiss_idx] = grada[tiss_idx] - 2*phan->mu_abs[tiss_idx]*(log(g[j])-log(gbar[j].real()))*dgbara[j].real()/(gbar[j].real());
                      }
                	}
        	    }
			}

            /*sprintf(tmpStr,"GradA%s.out",numStr);
            fid = fopen(tmpStr,"w");
            fwrite(grada,sizeof(double),phan->no_tiss,fid);
            fclose(fid); */


            //printf("Compute the gradient with respect to mus \n");
            // THE GRADIENT of mus
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
                	Neumann(src,out3,1);
    	            dgbars = generate_trans_image(out3,0);

        	        grads[tiss_idx] = 0.0;
            	    for (j = 0;j<geom->nX*geom->nY;j++) {
                        if(g[j] || gbar[j].real()){
						grads[tiss_idx] = grads[tiss_idx] - 2*phan->mu_sc[tiss_idx]*(log(g[j])-log(gbar[j].real()))*dgbars[j].real()/(gbar[j].real());
                	  } 
                    }
        	    }
			}

			fisher11 = 0;
			fisher12 = 0;
			fisher21 = 0;
			fisher22 = 0;
            
            for(j = 0;j<geom->nX*geom->nY;j++){
               if(gbar[j].real()){
                fisher11 = fisher11 + (1/gbar[j].real()) * dgbars[j].real()*dgbars[j].real();
                fisher12 = fisher12 + (1/gbar[j].real()) * dgbars[j].real()*dgbara[j].real();
                fisher21 = fisher21 + (1/gbar[j].real()) * dgbars[j].real()*dgbara[j].real();
                fisher22 = fisher22 + (1/gbar[j].real()) * dgbara[j].real()*dgbara[j].real(); 
                
               }
            }
 
            printf("%f %f %f %f \n", fisher11, fisher12, fisher21, fisher22);
         
		    sprintf(tmpStr,"tissmapabs%s.out",numStr);
            fid = fopen(tmpStr,"w");
		    sprintf(tmpStr,"tissmapsc%s.out",numStr);
            fid1 = fopen(tmpStr,"w");
            for (tiss_idx = 1;tiss_idx< phan->no_tiss;tiss_idx++) {
                if (fabs(grada[tiss_idx]) < thresh_abs)
					 flag_grada[tiss_idx] = 1;
				if ( fabs(grads[tiss_idx]) < thresh_sc)
					flag_grads[tiss_idx] = 1;

				if ( flag_grada[tiss_idx] == 0 ||  flag_grads[tiss_idx] == 0){
				printf("Old value of mua = %e mus = %e for tisstype = %d\n", phan->mu_abs[tiss_idx], phan->mu_sc[tiss_idx],tiss_idx);
					
				if ( flag_grada[tiss_idx] == 0){
                   mua_hat = log(phan->mu_abs[tiss_idx]/av_mua) - stepsizea*grada[tiss_idx];
	        	   phan->mu_abs[tiss_idx] = av_mua*exp(mua_hat);
				}
				if ( flag_grads[tiss_idx] == 0){
                   mus_hat = log(phan->mu_sc[tiss_idx]/av_mus) - stepsizes*grads[tiss_idx];
	        	   phan->mu_sc[tiss_idx] = av_mus*exp(mus_hat);
    	           /*phan->mu_sc[tiss_idx] = phan->mu_sc[tiss_idx] - av_mus*stepsizes*grads[tiss_idx];
				   if(phan->mu_sc[tiss_idx] < 0)
						phan->mu_sc[tiss_idx] = phan->mu_sc[tiss_idx] + 2*av_mus*stepsizes*grads[tiss_idx]; */
				}
            
				printf("Value of mua grad = %e musgrad = %e for tissue type %d \n", grada[tiss_idx], grads[tiss_idx], tiss_idx);
				printf("New value of mua = %e mus = %e for tissue type = %d\n", phan->mu_abs[tiss_idx], phan->mu_sc[tiss_idx],tiss_idx);
				}

				fwrite(&(phan->tiss_type[j]),sizeof(byte),1,fid);
	            fwrite(&(phan->tiss_type[j]),sizeof(byte),1,fid1);

            }
    
			/*updating the absorption and scatter coeff. values in the device code */

            populate_info_dev(); 

   
            fclose(fid);
            fclose(fid1);


            delg = 0.0;
            for (i=0;i<geom->nX*geom->nY;i++) {
                delg = delg + (log(gbar[i].real())-log(g[i]))*(log(gbar[i].real())-log(g[i]));
            }
            printf("Residual = %e\n",(delg).real());
            fwrite(&delg,sizeof(double),1,resFID);

        }

        fclose(resFID);

        free(g);
        free(beam_src);
        free(phan);
        free(geom);
        free(gbar);
        free(dgbara);
        free(dgbars);
        free(src);
        free(src2);
        free(out);
        free(out2);
        free(out3);
        free(out4);
        iniparser_freedict(ini);

		time(&end_time);
	
		printf("\n*------------------------------------------*\n");
		printf("\nThe total time taken by the code = %d sec \n", end_time - start_time);
		printf("\n*------------------------------------------*\n");

        return(0);
}


