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

    complex_double *gbar_trans, *gbar_ref,*dgbara_trans,*dgbara_ref,*dgbars_ref,*dgbars_trans;
//    complex_double *gbar_trans_orig, *gbar_ref_orig,*dgbara_trans_orig,*dgbara_ref_orig,*dgbars_ref_orig,*dgbars_trans_orig;
    complex_double *W, *W1, *out, *out2, *out3, *src, *tmp1, *df, *src2, *out4;
    

    flt_doub *g;
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
    phan = LoadPhantom(ini,0);

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

 
    FILE *abs_fid, *sc_fid, *grada_fid, *grads_fid, *res_fid, *fim_fid;
    abs_fid = fopen("abs_co_fim.dat","w");
    sc_fid = fopen("sc_co_fim.dat","w");
    res_fid = fopen("res_terms_fim.dat","w");
    fim_fid = fopen("fim_terms_size.dat", "w");

    for(sc_ind = 0; sc_ind < max_ind; sc_ind++){
        mus[sc_ind] = sc_ind*0.25 + 0.25;
    }
    for(abs_ind=0; abs_ind< max_ind; abs_ind++){
        mua[abs_ind] = abs_ind*0.1 + 0.1;
    }
                
    fwrite(mua, sizeof(flt_doub),max_ind, abs_fid);
    fwrite(mus, sizeof(flt_doub),max_ind, sc_fid);

    int loc_indx, loc_indy, loc_indz;
    int initz = 1;
  
/* 
    int ds_ind, max_ds_ind;
    max_ds_ind = 1;
    int ds_fac[max_ds_ind];
    for(i = 0; i < max_ds_ind; i++){
      ds_fac[i] = pow(2,i);
    //  printf("%d is ds_fac \n", ds_fac[i]);
    }
*/
    int rad_ind, max_rad_ind;
    float* rad_sig;
    max_rad_ind = 5;
    rad_sig = (float *) malloc(sizeof(float)*max_rad_ind);
    for (i=0; i< max_rad_ind; i++){
      rad_sig[i] = (i+1)*obj_fov_size/(2*geom->nX);
    }
    
//     for(sc_ind = 0; sc_ind < max_ind; sc_ind++){
//     for(loc_indy = 0; loc_indy < geom->nY; loc_indy++){
//      for(loc_indx = 0; loc_indx < geom->nX; loc_indx++){
  //       loc_indz = 10;
     //    sc_ind = 5;
         rad_ind = 2;
         loc_indy = geom->nY/2;
         loc_indx = geom->nX/2;
	   for(tiss_idx = 1; tiss_idx < phan->no_tiss; tiss_idx++){
        phan->mu_abs[tiss_idx] = 0.01; // mua[abs_ind];
        phan->mu_sc[tiss_idx] = 1.0;//0.25*sc_ind + 0.25;
       }
            
       populate_info_dev(); 
	   generate_diag_terms_host();
       memset(src, 0, sizeof(complex_double ) * geom->no_vox  * pow(nL+1,2));
       generate_source_beam(src);
       copy_dist(src,out2);
       copy_dist(src,W);
       Neumann(W,out2,1);

       gbar_ref = generate_ref_image(out2);
       //gbar_ref = (complex_double*) malloc(sizeof(complex_double)*geom->nX*geom->nY);
       //memcpy(gbar_ref, gbar_ref_orig,sizeof(complex_double)* geom->nX*geom->nY);
       //gbar_ref = downsample_image(gbar_ref_orig,ds_fac[ds_ind]);
      
/*       for(j=0; j<geom->nX*geom->nY; j++){
          if (j == geom->nX*(geom->nY/2-1) + geom->nX/2-1)
          printf("%f is gbar_ref[%d] \n", gbar_ref[j],j); 
       }*/
       //gbar_trans = generate_trans_image(out2,1);
       
      // THE GRADIENT of mua
     
    //  for(loc_indz = initz; loc_indz < geom->nZ-2; loc_indz = loc_indz+1){
    //   for(rad_ind = 0; rad_ind < max_rad_ind; rad_ind++){
         loc_indz = initz;
         rad_ind = 0;
//       printf("Compute the gradient with respect to mua \n");
       for (tiss_idx = 1; tiss_idx < phan->no_tiss; tiss_idx++) {
	            copy_dist(out2,out3);
				for (iip=0; iip<geom->nZ; iip++) {
			       for (jjp=0; jjp<geom->nY; jjp++) {
				      for (kkp=0; kkp<geom->nX; kkp++) {
	            	    r_ind = (iip + geom->bounZ)* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + (jjp + geom->bounY)* (geom->nX + 2*geom->bounX) + (kkp + geom->bounX);
             //           if(!(iip == loc_indz && jjp == loc_indy && kkp == loc_indx)){
             //     if( !(((iip == loc_indz) || (iip == loc_indz - 1) || (iip == loc_indz +1)) && (jjp == loc_indy || jjp == loc_indy - 1 || jjp == loc_indy + 1)  && (kkp == loc_indx || kkp == loc_indx + 1 || kkp == loc_indx-1))){
                          if( ((x[kkp]-x[loc_indx])*(x[kkp]-x[loc_indx]) + (y[jjp]-y[loc_indy])*(y[jjp]-y[loc_indy]) + (z[iip] - z[loc_indz])*(z[iip] - z[loc_indz])) > rad_sig[rad_ind]*rad_sig[rad_ind]){
                        	for (n=0;n<(nL+1)*(nL+1);n++) {
                            	out3[VOX_TO_SPIND(r_ind, n,(nL+1)*(nL+1))] = 0+0*I;
	                        }
    	                }
                      //}}
                        //else
                          //printf("%d %d %d \n", iip, jjp, kkp);
        	     }}} 
            	 scale_dist(out3,-1.0*C);
	             copy_dist(out3,src);
                 Neumann(src,out3,1);
            	 dgbara_ref = generate_ref_image(out3);
       	/*		  for(j=0; j<geom->nX*geom->nY; j++){
                   if (j == geom->nX*(geom->nY/2-1) + geom->nX/2-1)
                     printf("%f is dgbara_ref[%d] \n", dgbara_ref[j],j); 
                 }*/
            	 //dgbara_ref = downsample_image(dgbara_ref_orig,ds_fac[ds_ind]);
            	 //dgbara_trans = generate_trans_image(out3,1);
			}
#if 0
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
                      //  if(!(iip == loc_indz && jjp == loc_indy && kkp == loc_indx)){
                      //    if( !(((iip == loc_indz) || (iip == loc_indz - 1) || (iip == loc_indz +1)) && (jjp == loc_indy || jjp == loc_indy - 1 || jjp == loc_indy + 1)  && (kkp == loc_indx || kkp == loc_indx + 1 || kkp == loc_indx-1))){
                          if( ((x[kkp]-x[loc_indx])*(x[kkp]-x[loc_indx]) + (y[jjp]-y[loc_indy])*(y[jjp]-y[loc_indy]) + (z[iip] - z[loc_indz])*(z[iip] - z[loc_indz])) > rad_sig[rad_ind]*rad_sig[rad_ind]){
	            	            r_ind = (iip + geom->bounZ)* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + (jjp + geom->bounY)* (geom->nX + 2*geom->bounX) + (kkp + geom->bounX);
        		                for (n=0;n<(nL+1)*(nL+1);n++) {
                		            out3[VOX_TO_SPIND(r_ind, n, (nL+1)*(nL+1))] = 0+0*I;
                        		}
                    		}
                      //}}
	                }}}

        	        src = alloc_dist();
            	    copy_dist(out3,src);
                	Neumann(src,out3,1);
    	            dgbars_ref = generate_ref_image(out3);
            	    //dgbars_ref = downsample_image(dgbars_ref_orig,ds_fac[ds_ind]);
    	            //dgbars_trans = generate_trans_image(out3,1);
            }
#endif

            memset(fisher_mat, 0, sizeof(flt_doub)*4);

            for(j = 0;j<geom->nX*geom->nY;j++){
               if(gbar_ref[j].real()){
/*                fisher_mat[0] += (1/gbar_ref[j].real()) * dgbars_ref[j].real()*dgbars_ref[j].real();
                fisher_mat[1] += (1/gbar_ref[j].real()) * dgbars_ref[j].real()*dgbara_ref[j].real();
                fisher_mat[2] += (1/gbar_ref[j].real()) * dgbars_ref[j].real()*dgbara_ref[j].real(); */ 
                fisher_mat[3] += (1/gbar_ref[j].real()) * dgbara_ref[j].real()*dgbara_ref[j].real(); 
               }
             /*  if(gbar_trans[j].real()){
                fisher_mat[0] += (1/gbar_trans[j].real()) * dgbars_trans[j].real()*dgbars_trans[j].real();
                fisher_mat[1] += (1/gbar_trans[j].real()) * dgbars_trans[j].real()*dgbara_trans[j].real();
                fisher_mat[2] += (1/gbar_trans[j].real()) * dgbars_trans[j].real()*dgbara_trans[j].real();
                fisher_mat[3] += (1/gbar_trans[j].real()) * dgbara_trans[j].real()*dgbara_trans[j].real();
               } */
            }
            printf("%e is fisher abs for %d as z index and %f as radius \n", fisher_mat[3], loc_indz, rad_sig[rad_ind]);
            //printf("%e is fisher abs for %f as depth and %f as sc co \n", fisher_mat[3], z[loc_indz], phan->mu_sc[1]);

             
            //fwrite(&fisher_mat[3], sizeof(flt_doub), 1, fim_fid);

         //}}}
        // }
        //}
 

        fclose(abs_fid); 
        fclose(sc_fid); 
        fclose(res_fid); 
        fclose(fim_fid); 

        free(g);
        free(beam_src);
        free(phan);
        free(geom);
        free(gbar_ref);
//        free(gbar_ref_orig);
        free(gbar_trans);
        free(dgbara_ref);
        free(dgbars_ref);
//        free(dgbara_ref_orig);
//        free(dgbars_ref_orig);
        free(dgbara_trans);
        free(dgbars_trans);
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


