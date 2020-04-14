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

    complex_double *gbar,*dgbar;
    complex_double *W, *W1, *out, *out2, *out3, *src, *tmp1, *df, *src2, *out4;
    

    flt_doub *g;
    flt_doub *grada,*grads;
    flt_doub *grada2,*grads2;
    flt_doub delg;

    int n;
    int j,jj;
    int jnk;
    int tiss_idx;
    int iip, jjp, kkp;
    int r_ind,i;
    int size;
    float av_mua; 
    float av_mus;

    av_mua = 0.01;
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


    double stepsizea, stepsizes;
    stepsizea = iniparser_getdouble(ini,"Runtime:stepsizea",1e3);
    stepsizes = iniparser_getdouble(ini,"Runtime:stepsizes",1e3);

	int jnk2;
    FILE *gFile;
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
	flt_doub tmp;
     
    grada = (flt_doub *)malloc(sizeof(flt_doub)*phan->no_tiss);
    grads = (flt_doub *)malloc(sizeof(flt_doub)* phan->no_tiss);

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
    float mus_hat, mua_hat;

    int abs_ind, sc_ind;
    int max_ind = 10;
    flt_doub *mua, *mus;
    mua = (flt_doub *) malloc(sizeof(flt_doub)*max_ind);  
    mus = (flt_doub *) malloc(sizeof(flt_doub)*max_ind); 
 
    FILE *abs_fid, *sc_fid, *grada_fid, *grads_fid, *res_fid1, *res_fid2,*grada_fid2,*grads_fid2;
    abs_fid = fopen("abs_co.dat","w");
    sc_fid = fopen("sc_co.dat","w");
    grada_fid = fopen("grada_terms_high.dat","w");
    grads_fid = fopen("grads_terms_high.dat","w");
    grada_fid2 = fopen("grada_terms2_high.dat","w");
    grads_fid2 = fopen("grads_terms2_high.dat","w");
    res_fid1 = fopen("res_terms1_high.dat","w");
    res_fid2 = fopen("res_terms2_high.dat","w");

    for(sc_ind = 0; sc_ind < max_ind; sc_ind++){
        mus[sc_ind] = sc_ind*0.25 + 0.25;
    }
    for(abs_ind=0; abs_ind< max_ind; abs_ind++){
        mua[abs_ind] = abs_ind*0.0025 + 0.0025;
    }
                
    fwrite(mua, sizeof(flt_doub),max_ind, abs_fid);
    fwrite(mus, sizeof(flt_doub),max_ind, sc_fid);
   
    for(sc_ind = 0; sc_ind < max_ind; sc_ind++){
      for(abs_ind=0; abs_ind< max_ind; abs_ind++){
	   for(tiss_idx = 1; tiss_idx < phan->no_tiss; tiss_idx++){
        phan->mu_abs[tiss_idx] = mua[abs_ind];
        phan->mu_sc[tiss_idx] = mus[sc_ind];
       }
            
       populate_info_dev(); 
	   generate_diag_terms_host();
       generate_source_beam(src);
       copy_dist(src,out2);
       copy_dist(src,W);
       Neumann(W,out2,1);

       // THE GRADIENT of mua
       gbar = generate_trans_image(out2,1);

#if 1  
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
            	 dgbar = generate_trans_image(out3,0);
                 grada[tiss_idx] = 0.0;
                 grada2[tiss_idx] = 0.0;
                 for (j = 0;j<geom->nX*geom->nY;j++) {
                    if(g[j] ||  gbar[j].real()){
        	            grada[tiss_idx] = grada[tiss_idx] - 2*phan->mu_abs[tiss_idx]*(log(g[j])-log(gbar[j].real()))*dgbar[j].real()/(gbar[j].real());
        	            grada2[tiss_idx] = grada[tiss_idx] - 2*((g[j])-(gbar[j].real()))*dgbar[j].real();
                      }
                	}
        	    }
			}

            //printf("Compute the gradient with respect to mus \n");
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
    	            dgbar = generate_trans_image(out3,0);

        	        grads[tiss_idx] = 0.0;
            	    for (j = 0;j<geom->nX*geom->nY;j++) {
                        if(g[j] || gbar[j].real()){
						grads[tiss_idx] = grads[tiss_idx] - 2*phan->mu_sc[tiss_idx]*(log(g[j])-log(gbar[j].real()))*dgbar[j].real()/(gbar[j].real());
        	            grads2[tiss_idx] = grads[tiss_idx] - 2*((g[j])-(gbar[j].real()))*dgbar[j].real();
                	  } 
                    }
        	    }
			}

            for (tiss_idx = 1;tiss_idx< phan->no_tiss;tiss_idx++) {
				printf("Old value of mua = %e mus = %e for tissue type = %d\n", phan->mu_abs[tiss_idx], phan->mu_sc[tiss_idx],tiss_idx);
                   mua_hat = log(phan->mu_abs[tiss_idx]/av_mua) - stepsizea*grada[tiss_idx];
	        	   phan->mu_abs[tiss_idx] = av_mua*exp(mua_hat);
                   mus_hat = log(phan->mu_sc[tiss_idx]/av_mus) - stepsizes*grads[tiss_idx];
	        	   phan->mu_sc[tiss_idx] = av_mus*exp(mus_hat);
				printf("Value of mua grad = %e musgrad = %e for tissue type %d \n", grada[tiss_idx], grads[tiss_idx], tiss_idx);
				printf("New value of mua = %e mus = %e for tissue type = %d\n", phan->mu_abs[tiss_idx], phan->mu_sc[tiss_idx],tiss_idx);
		     }

#endif
    
            delg = 0.0;
            delg2 = 0.0
            for (i=0;i<geom->nX*geom->nY;i++) {
                if(g[j] || gbar[j].real()){
                   delg = delg + (log(gbar[i].real())-log(g[i]))*(log(gbar[i].real())-log(g[i]));
                   delg2 = delg2 + ((gbar[i].real())-(g[i]))*((gbar[i].real())-(g[i]));
                }
            }
            tiss_idx = 1;
            printf("Residual1 = %e mua grad = %e musgrad = %e \n",(delg),  grada[tiss_idx], grads[tiss_idx]);
            fwrite(grada, sizeof(flt_doub),1, grada_fid);
            fwrite(grads, sizeof(flt_doub),1, grads_fid);
            fwrite(&delg,sizeof(flt_doub),1,res_fid1);
            fwrite(&delg2,sizeof(flt_doub),1,res_fid2);
           

         }
        }

        fclose(abs_fid); 
        fclose(sc_fid); 
        fclose(grada_fid); 
        fclose(grads_fid); 
        fclose(res_fid1); 
        fclose(res_fid2); 

        free(g);
        free(beam_src);
        free(phan);
        free(geom);
        free(gbar);
        free(dgbar);
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


