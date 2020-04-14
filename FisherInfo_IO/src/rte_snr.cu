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

    complex_double *gbar1,*gbar0, *gbar, *dgbara, *dgbars;
    complex_double *W, *W1, *out, *out2, *out3, *src, *tmp1, *df, *src2, *out4;
    

    flt_doub *g;
    flt_doub *grada,*grads;
    flt_doub delg;

    int n;
    int j,jj;
    int jnk;
    int tiss_idx;
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
 
	int r_ind_phan;
    W = alloc_dist();
    out2 = alloc_dist();
    out3 = alloc_dist();
    out4 = alloc_dist();
    src = alloc_dist();
    src2 = alloc_dist();
    tmp1 = alloc_dist();

    int abs_ind, sc_ind;
    int max_ind = 10;
    flt_doub *mua, *mus;
    mua = (flt_doub *) malloc(sizeof(flt_doub)*max_ind);  
    mus = (flt_doub *) malloc(sizeof(flt_doub)*max_ind);

    flt_doub *fisher_mat;
    fisher_mat = (flt_doub*) malloc(sizeof(flt_doub)*4);

    flt_doub *x,*y,*z;
    x = (flt_doub *)malloc(sizeof(flt_doub)*geom->nX);
    y = (flt_doub *)malloc(sizeof(flt_doub)*geom->nY);
    z = (flt_doub *)malloc(sizeof(flt_doub)*geom->nZ);

    flt_doub obj_fov_size;
    obj_fov_size =  geom->x_max - geom->x_min;
    for(i=0; i<geom->nX; i++){
      x[i] = i*obj_fov_size/(geom->nX) - obj_fov_size/2 + obj_fov_size/(2*geom->nX);
      y[i] = i*obj_fov_size/(geom->nY) - obj_fov_size/2 + obj_fov_size/(2*geom->nY); // Change both of these to obj->nX-1 and obj->nY-1, respectively?
      z[i] = i*obj_fov_size/(geom->nZ) + obj_fov_size/(2*geom->nZ); // Change both of these to obj->nX-1 and obj->nY-1, respectively?
    }
 
    FILE *abs_fid, *sc_fid, *grada_fid, *grads_fid, *res_fid, *fim_fid, *snr_fid;
    abs_fid = fopen("abs_co_fim.dat","w");
    sc_fid = fopen("sc_co_fim.dat","w");
    grada_fid = fopen("grada_terms_fim.dat","w");
    grads_fid = fopen("grads_terms_fim.dat","w");
    res_fid = fopen("res_terms_fim.dat","w");
    fim_fid = fopen("fim_terms.dat", "w");
    snr_fid = fopen("snr_terms.dat", "w");

    for(sc_ind = 0; sc_ind < max_ind; sc_ind++){
        mus[sc_ind] = sc_ind*0.25 + 0.25;
    }
    for(abs_ind=0; abs_ind< max_ind; abs_ind++){
        mua[abs_ind] = abs_ind*0.0025 + 0.0025;
    }
                
    fwrite(mua, sizeof(flt_doub),max_ind, abs_fid);
    fwrite(mus, sizeof(flt_doub),max_ind, sc_fid);

    int rad_ind, max_rad_ind;
    float* rad_sig;
    max_rad_ind = 5;
    rad_sig = (float *) malloc(sizeof(float)*max_rad_ind);
    for (i=0; i< max_rad_ind; i++){
      rad_sig[i] = (i+1)*obj_fov_size/(2*geom->nX);
    }
    rad_ind = 2;
    int loc_indx, loc_indy, loc_indz;
    int hyp_no;
    flt_doub small_sig = 0.1;
    flt_doub *snr_sq;
    int initz = 2;
    snr_sq = (flt_doub* ) malloc(sizeof(flt_doub)*(geom->nZ-initz));
    int ind;
   
    for(loc_indz = initz; loc_indz < 16; loc_indz = loc_indz++){
//     for(loc_indy = 0; loc_indy < geom->nY; loc_indy++){
//      for(loc_indx = 0; loc_indx < geom->nX; loc_indx++){
         loc_indy = geom->nY/2;
         loc_indx = geom->nX/2;
         for (i=0; i<geom->nZ; i++){
          for (j=0; j<geom->nY; j++){
           for (k=0; k<geom->nX; k++){
             r_ind = i*geom->nX*geom->nY + j*geom->nX + k;
             if(((x[k] - x[loc_indx])*(x[k]-x[loc_indx]) + (y[j]-y[loc_indy])*(y[j]-y[loc_indy]) + (z[i] - z[loc_indz])*(z[i] - z[loc_indz])) < rad_sig[rad_ind]*rad_sig[rad_ind])
                phan->tiss_type[r_ind] = 2;
             else
                phan->tiss_type[r_ind] = 1;
         }}}
         for(hyp_no=1; hyp_no >= 0; hyp_no--){
	       for(tiss_idx = 1; tiss_idx < phan->no_tiss; tiss_idx++){
               if(tiss_idx == 2){
                phan->mu_abs[tiss_idx] = 0.25 + small_sig*hyp_no; // mua[abs_ind];
                phan->mu_sc[tiss_idx] = 1.0; 
               }
              }
       	      populate_info_dev(); 
 		      generate_diag_terms_host();
   	          generate_source_beam(src);
     	      copy_dist(src,out2);
 	          copy_dist(src,W);
    	      Neumann(W,out2,1);

          if(hyp_no == 1){
             printf("Generating signal present image\n");
       	     gbar1 = generate_ref_image(out2);
          }
          else{
             printf("Generating signal absent image \n");
             gbar0 = generate_ref_image(out2);
           }
         }
        
       snr_sq[loc_indz-initz] = 0.0;
       printf("Computing the SNR \n"); 
       for (ind=0; ind<geom->nX*geom->nY; ind++){
            snr_sq[loc_indz-initz] += ((gbar1[ind].real() - gbar0[ind].real())*(gbar1[ind].real() - gbar0[ind].real()))/(gbar0[ind].real());
    //        if(fabs(gbar1[ind].real() - gbar0[ind].real())>1e-15)
      //          printf("%lf  %lf \n ", gbar1[ind].real(), gbar0[ind].real());
           if(ind == geom->nX*(geom->nY/2) + geom->nX/2 )
              printf("%f is gbar1 and %f is gbar0 for voxel %d \n", gbar1[ind].real(), gbar0[ind].real(), ind); 
       }
       printf("%e is snr_sq for %d as z index\n", snr_sq[loc_indz-initz], loc_indz);
#if 0
       //printf("Compute the gradient with respect to mua");
       for (tiss_idx = 1; tiss_idx < phan->no_tiss; tiss_idx++) {
	            copy_dist(out2,out3);
				for (iip=0; iip<geom->nZ; iip++) {
			       for (jjp=0; jjp<geom->nY; jjp++) {
				      for (kkp=0; kkp<geom->nX; kkp++) {
	            	    r_ind = (iip + geom->bounZ)* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + (jjp + geom->bounY)* (geom->nX + 2*geom->bounX) + (kkp + geom->bounX);
                        if(!(iip == loc_indz && jjp == loc_indy && kkp == loc_indx)){
                        	for (n=0;n<(nL+1)*(nL+1);n++) {
                            	out3[VOX_TO_SPIND(r_ind, n,(nL+1)*(nL+1))] = 0+0*I;
	                        }
    	                }
                        else
                           printf("%d %d %d \n", iip, jjp, kkp);
        	     }}}
            	 scale_dist(out3,-1.0*C);
	             copy_dist(out3,src);
                 Neumann(src,out3,1);
            	 dgbara = generate_trans_image(out3,0);
                 grada[tiss_idx] = 0.0;
                 for (j = 0;j<geom->nX*geom->nY;j++) {
        	            grada[tiss_idx] = grada[tiss_idx] - 2*(g[j]-gbar[j].real())*dgbara[j].real();
                	}
			}

            //printf("Compute the gradient with respect to mus \n");
            for (tiss_idx = 1; tiss_idx < phan->no_tiss ; tiss_idx++) {
	                copy_dist(out2,out3);
    	            scale_dist(out3,-1.0*C);
        	        copy_dist(out2,tmp1);
            	    PropScatmu1(geom,phan,nL,tmp1);
    	            add_dist(tmp1,out3,out3);
					for (iip=0; iip<geom->nZ; iip++) {
				      for (jjp=0; jjp<geom->nY; jjp++) {
				        for (kkp=0; kkp<geom->nX; kkp++) {
                        if(!(iip == loc_indz && jjp == loc_indy && kkp == loc_indx)){
	            	            r_ind = (iip + geom->bounZ)* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + (jjp + geom->bounY)* (geom->nX + 2*geom->bounX) + (kkp + geom->bounX);
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
                        grads[tiss_idx] = grads[tiss_idx] - (g[j]-gbar[j].real())*dgbars[j].real();
                    }
            }

            memset(fisher_mat, 0, sizeof(flt_doub)*4);

            for(j = 0;j<geom->nX*geom->nY;j++){
               if(gbar[j].real()){
                fisher_mat[0] += (1/gbar[j].real()) * dgbars[j].real()*dgbars[j].real();
                fisher_mat[1] += (1/gbar[j].real()) * dgbars[j].real()*dgbara[j].real();
                fisher_mat[2] += (1/gbar[j].real()) * dgbars[j].real()*dgbara[j].real();
                fisher_mat[3] += (1/gbar[j].real()) * dgbara[j].real()*dgbara[j].real();

               }
            }

             
            fwrite(fisher_mat, sizeof(flt_doub), 4, fim_fid);
            fwrite(grads, sizeof(flt_doub), 1, grads_fid);
            fwrite(grada, sizeof(flt_doub), 1, grada_fid);

    
            delg = 0.0;
            for (i=0;i<geom->nX*geom->nY;i++) {
                delg = delg + (gbar[i].real()-g[i])*(gbar[i].real()-g[i]);
            }
            printf("Residual = %e\n",(delg));
            fwrite(&delg,sizeof(flt_doub),1,res_fid);
          }
         }
#endif
        }
        fwrite(snr_sq, sizeof(flt_doub),geom->nZ-initz, snr_fid);

        fclose(snr_fid); 

        fclose(abs_fid); 
        fclose(sc_fid); 
        fclose(grada_fid); 
        fclose(grads_fid); 
        fclose(res_fid); 
        fclose(fim_fid); 

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
        free(fisher_mat);
        iniparser_freedict(ini);

		time(&end_time);
	
		printf("\n*------------------------------------------*\n");
		printf("\nThe total time taken by the code = %d sec \n", end_time - start_time);
		printf("\n*------------------------------------------*\n");

        return(0);
}


