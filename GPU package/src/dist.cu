#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include<cuda.h>
#include <omp.h>

#include "iniparser.h"
#include "rte.h"

#pragma warning(disable:981)


extern SHORT nL;
extern Geometry *geom;
extern Phantom *phan;
extern Source *beam_src;
extern complex_double *diag_terms_host;
extern complex_double *sph_harm;

complex_double* alloc_dist(){
	complex_double *dist;
	if ( (dist = (complex_double *) malloc(sizeof(complex_double ) * geom->no_vox  * pow(nL+1,2))) == NULL){
		printf("Error in memory allocation for distribution \n");
		exit(0);
	}
    memset(dist, 0, sizeof(complex_double ) * geom->no_vox  * pow(nL+1,2));
	return dist;
}

complex_double* alloc_dist_dev(){
    complex_double *dist;
    MY_SAFE_CALL(cudaMalloc(&dist,sizeof(complex_double ) * geom->no_vox * pow(nL+1,2)));
    return dist;
}

void add_dist(complex_double *W1, complex_double *W2, complex_double *out){
	int i, cnt;

	for ( i = 0; i < geom->no_vox; i++){
		for ( cnt = 0; cnt < pow(nL+1,2); cnt ++){
			out[VOX_TO_SPIND(i,cnt,(nL+1)*(nL+1))] = W1[VOX_TO_SPIND(i,cnt,(nL+1)*(nL+1))] + W2[VOX_TO_SPIND(i,cnt,(nL+1)*(nL+1))];
		}
	}
}

void copy_dist(complex_double *W1, complex_double *W2){
	int r_ind, cnt;
	for ( r_ind = 0; r_ind < geom->no_vox; r_ind++){
		for (cnt = 0; cnt < pow(nL+1,2); cnt++){
			W2[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))] = W1[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))];
		}
	}	
}

void scale_dist(complex_double *W1, double fac){
	int r_ind, cnt;
	for ( r_ind = 0; r_ind < geom->no_vox; r_ind++){
		for (cnt = 0; cnt < pow(nL+1,2); cnt++){
			W1[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))] = fac*W1[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))];
		}
	}	
}


void generate_source_beam(complex_double *src_dist){
	int i,j,k,cnt, l,m, r_ind;
	int r_indp, ip;
	doublecomplex tmp;
	byte tiss_type;

	for (i=geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
      for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
	    for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
			if ( sqrt(pow(beam_src->pos_x - geom->xs[k],2)  + pow(beam_src->pos_y - geom->ys[j],2)) < beam_src->diam/2.0){
//			if( k == geom->bounX + geom->nX/2.0 && j == geom->bounY + geom->nY/2.0){
			cnt = 0;
			for (l = 0; l <= nL; l++){
				for ( m = -l; m <= l; m++){
					if ( m == 0){ 
						r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
						tmp = 1.0; 
                        for(ip = geom->bounZ; ip < i; ip++){
							r_indp = ip* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k; 
							tiss_type = phan->tiss_type[r_indp];
							tmp = tmp* (cexp(-(phan->mu_sc[tiss_type] + phan->mu_abs[tiss_type]) * geom->delZ));
                            //tmp = tmp * ((-(phan->mu_sc[r_indp] +phan->mu_abs[r_indp])  * geom->delZ ).cexp());
                        }
						tiss_type = phan->tiss_type[r_ind];
					//	printf("tmp is %e for (%d %d %d) voxel\n", tmp.real(),i,j,k);
						src_dist[VOX_TO_SPIND(r_ind,cnt,(nL+1)*(nL+1))] = pow(phan->g,l) * phan->mu_sc[tiss_type] / (phan->mu_sc[tiss_type] + phan->mu_abs[tiss_type]) * sqrt(((2.0*(double)l)+1.0)/(4.0*M_PI))  * tmp * (1 - cexp(-(phan->mu_sc[tiss_type]+phan->mu_abs[tiss_type])*geom->delZ)) / (geom->delX * geom->delY * geom->delZ);
//					printf("%e + %e is the source dist generated after the KXE term in the (%d %d %d) voxel whose linear ind is %d for cnt = %d \n",  src_dist[VOX_TO_SPIND(r_ind,cnt,geom->no_vox)].real(), (src_dist[VOX_TO_SPIND(r_ind,cnt,geom->no_vox)] ).imag(), i-geom->bounZ,j-geom->bounY,k-geom->bounX, VOX_TO_SPIND(r_ind,cnt,geom->no_vox),cnt);
					}
				cnt++;
				}
			}
			}
			}
		}
	}
}

complex_double * generate_trans_image(complex_double *out_dist, int flag ){ // flag =1 => compute XE term else dont. 
	int i,j,k;
	int num_xs;
    int sp, l;
	num_xs = 100;
	complex_double val;

	complex_double *trans_image;
	flt_doub *xs, *out;

	xs = (flt_doub *)malloc(num_xs*sizeof(flt_doub));
    out =  (flt_doub *)malloc(num_xs*sizeof(flt_doub));

    for (i=0;i<num_xs;i++) {
        xs[i] = (flt_doub)i/((flt_doub)num_xs - 1.0);
    }

	trans_image = (complex_double *)malloc (sizeof(complex_double)*geom->nX * geom->nY);

	int r_ind, r_indp, ip;
	complex_double tmp;
	i = geom->bounZ + geom->nZ -1;
    for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
       for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
			val = 0.0 + 0.0*I ;
			r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;	
			if (flag == 1){			
			if ( sqrt(pow(beam_src->pos_x - geom->xs[k],2)  + pow(beam_src->pos_y - geom->ys[j],2)) <= beam_src->diam/2.0){
//			if( k == geom->bounX + geom->nX/2.0 && j == geom->bounY + geom->nY/2.0){
				tmp = beam_src->mag*1.0;
				for(ip=geom->bounZ ; ip < geom->nZ + geom->bounZ; ip++){
	  				 r_indp = ip* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
				     tmp = tmp * cexp(-(phan->mu_sc[phan->tiss_type[r_indp]] +phan->mu_abs[phan->tiss_type[r_indp]]) * geom->delZ);
				}

				printf("%e is the value of the 0 term for voxel (%d %d %d) with r_ind as %d \n", tmp.real(),i,j,k,r_ind);
                val = val+ tmp;   // mu_int(cor_src,cor_dest,subvox_src, subvox_dest);
			}
			}
			for (l=0;l<= nL;l++) {
/*                if(out_dist[VOX_TO_SPIND(r_ind,(l*l+l),geom->no_vox)].real() > 0.0){
					printf("out_dist is %e for voxel (%d %d %d) with r_ind as %d and cnt = %d, l = %d \n", out_dist[VOX_TO_SPIND(r_ind,(l*l+l),geom->no_vox)], i,j,k,r_ind,l*l+l,l);
				} */
                AssocLegendreArray(l,0,num_xs,xs,out);
                for (sp = 0;sp<num_xs;sp++) {
                    val = val+ 2.0*M_PI*geom->delX * geom->delY / sqrt(4.0*M_PI) * sqrt(2.0*(flt_doub)l +1.0) *
                        out_dist[VOX_TO_SPIND(r_ind,(l*l+l),(nL+1)*(nL+1))] *
                        out[sp]*xs[sp] / ((flt_doub)num_xs - 1.0);
                }
            }
			trans_image[(j-geom->bounY)*geom->nX + (k-geom->bounX)] = val;
		}
	}
	free(out);
	free(xs);
	return trans_image;
}

double * generate_phase_image(complex_double *out_dist, int flag ){ // flag =1 => compute XE term else dont. 
	int i,j,k;
	int num_xs;
    int sp, l;
	num_xs = 100;
	complex_double val;

	complex_double *ret;
	double *phase_image;
	flt_doub *xs, *out;

	xs = (flt_doub *)malloc(num_xs*sizeof(flt_doub));
    out =  (flt_doub *)malloc(num_xs*sizeof(flt_doub));

    for (i=0;i<num_xs;i++) {
        xs[i] = (flt_doub)i/((flt_doub)num_xs - 1.0);
    }

	phase_image = (double *)malloc (sizeof(double)*geom->nX * geom->nY);
	ret = (complex_double *)malloc (sizeof(complex_double)*geom->nX * geom->nY);

	int r_ind, r_indp, ip;
	complex_double tmp;
	i = geom->bounZ + geom->nZ -1;
    for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
       for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
			val = 0.0 + 0.0*I ;
			r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;	
			if (flag == 1){			
			if ( sqrt(pow(beam_src->pos_x - geom->xs[k],2)  + pow(beam_src->pos_y - geom->ys[j],2)) < beam_src->diam/2.0){
				tmp = beam_src->mag*1.0;
				for(ip=geom->bounZ ; ip < geom->nZ + geom->bounZ; ip++){
	  				 r_indp = ip* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
				     tmp = tmp * cexp(-(phan->mu_sc[phan->tiss_type[r_indp]] +phan->mu_abs[phan->tiss_type[r_indp]]) * geom->delZ);
				}

				printf("%e is the value of the 0 term for voxel (%d %d %d) with r_ind as %d \n", tmp.real(),i,j,k,r_ind);
                val = val+ tmp;   // mu_int(cor_src,cor_dest,subvox_src, subvox_dest);
			}
			}
			for (l=0;l<= nL;l++) {
/*                if(out_dist[VOX_TO_SPIND(r_ind,(l*l+l),geom->no_vox)].real() > 0.0){
					printf("out_dist is %e for voxel (%d %d %d) with r_ind as %d and cnt = %d, l = %d \n", out_dist[VOX_TO_SPIND(r_ind,(l*l+l),geom->no_vox)], i,j,k,r_ind,l*l+l,l);
				} */
                AssocLegendreArray(l,0,num_xs,xs,out);
                for (sp = 0;sp<num_xs;sp++) {
                    val = val+ 2.0*M_PI*geom->delX * geom->delY / sqrt(4.0*M_PI) * sqrt(2.0*(flt_doub)l +1.0) *
                        out_dist[VOX_TO_SPIND(r_ind,(l*l+l),(nL+1)*(nL+1))] *
                        out[sp]*xs[sp] / ((flt_doub)num_xs - 1.0);
                }
            }
			ret[(j-geom->bounY)*geom->nX + (k-geom->bounX)] = val;
		}
	}

    for(i=0; i< geom->nY; i++){
		for(j=0; j<geom->nX; j++){
			phase_image[i*geom->nX+j] = carg(ret[i*geom->nX+j]);
		}
	}
    
	free(out);
	free(xs);
	return phase_image;
}

complex_double * generate_ref_image(complex_double *out_dist){  
	int i,j,k;
	int num_xs;
    int sp, l;
	num_xs = 100;
	complex_double val;

	complex_double *ref_image;
	flt_doub *xs, *out;

	xs = (flt_doub *)malloc(num_xs*sizeof(flt_doub));
    out =  (flt_doub *)malloc(num_xs*sizeof(flt_doub));

    for (i=0;i<num_xs;i++) {
        xs[i] = (flt_doub)i/((flt_doub)num_xs - 1.0);
    }

	ref_image = (complex_double *)malloc (sizeof(complex_double)*geom->nX * geom->nY);

	int r_ind, r_indp, ip;
	complex_double tmp;
	i = geom->bounZ;
    for (j=geom->bounY;j<geom->nY + geom->bounY;j++) {
       for (k=geom->bounX;k<geom->nX + geom->bounX;k++) {
			val = 0.0 + 0.0*I ;
			r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;	
			for (l=0;l<= nL;l++) {
                AssocLegendreArray(l,0,num_xs,xs,out);
                for (sp = 0;sp<num_xs;sp++) {
                    val = val+ 2.0*M_PI*geom->delX * geom->delY / sqrt(4.0*M_PI) * sqrt(2.0*(flt_doub)l +1.0) *
                        out_dist[VOX_TO_SPIND(r_ind,(l*l+l),(nL+1)*(nL+1))] *
                        out[sp]*xs[sp] / ((flt_doub)num_xs - 1.0);
                }
            }
			ref_image[(j-geom->bounY)*geom->nX + (k-geom->bounX)] = val;
		}
	}
	free(out);
	free(xs);
	return ref_image;
}


void output_phase_image(double *transImg,char *name) {
    FILE *fle;
    int i,j;
    double val;

    fle = fopen(name,"w");

   int wid, ht;

    wid = (int) geom->nY;
    ht = (int) geom->nX;

    fwrite(&wid,sizeof(int),1,fle);
    fwrite(&ht,sizeof(int),1,fle);

    for (i=0;i<geom->nY;i++) {
        for (j=0;j<geom->nX;j++) {
            val = transImg[i*geom->nX+j];
            fwrite(&val,sizeof(double),1,fle);
        }
    }

    fclose(fle);
}

void output_image(complex_double *trans_img,char *name) {
    FILE *fle;
    int i,j;
    flt_doub val;

    fle = fopen(name,"w");
    int wid, ht;

    wid = (int) geom->nY;
	ht = (int) geom->nX;

    fwrite(&wid,sizeof(int),1,fle);
    fwrite(&ht,sizeof(int),1,fle);

    for (i=0;i<geom->nY;i++) {
		for (j=0;j<geom->nX;j++) {
            val = (trans_img[i*geom->nX+j]).real();
            fwrite(&val,sizeof(flt_doub),1,fle);
//			printf("%f (%d %d) \n", val, i,j);
        }
    }
    for (i=0;i<geom->nY;i++) {
		for (j=0;j<geom->nX;j++) {
            val = (trans_img[i*geom->nX+j]).imag();
            fwrite(&val,sizeof(flt_doub),1,fle);
        }
    }

    fclose(fle);

}

void generate_diag_terms_host_invoke(){

/*
	THREAD_PAR_DIAG thread_parameters[NUM_CPU_THREADS];

	int i;

    for(i = 0; i < (int) ceilf(phan->no_tiss/NUM_CPU_THREADS); i++){
		
	}
*/


}
void generate_diag_terms_host() {

    int ang_res = ANG_RES;
    int theta_count, phi_count, omega_count;
    flt_doub *theta_self, *phi_self;
    float cube_x, cube_y, cube_z;
    int face_calc;
    flt_doub *sph_x, *sph_y, *sph_z, face_x, face_y, face_z;
    flt_doub dist_self;
    flt_doub del_phi, del_theta;
    int r_ind_self;
    int ang_ind ;
    complex_double **fact_self_vox;
    complex_double *rt_self, *rtp_self;
    complex_double sub_v_sum_self;
    float ii_self, jj_self, kk_self;
    int face = 1; // Defines the face of the cube to be considered later.
    int l,m,lp,mp,cnt,cntp;

    int i;


    diag_terms_host = (complex_double *)malloc(sizeof(complex_double )* MAX_TISS_NUM*MAX_NL*MAX_NL*MAX_NL*MAX_NL);
    
	fact_self_vox = (complex_double **)malloc ( sizeof(complex_double *) * pow(geom->self_sub_vox,3));
    for ( i = 0; i < pow(geom->self_sub_vox,3); i++){
        fact_self_vox[i] = (complex_double *)malloc(sizeof(complex_double) * pow(ang_res,2));
        if (fact_self_vox[i] == 0L) {
            printf("Memory allocation error \n");
        }
    }


    sph_x = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));
    sph_y = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));
    sph_z = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));

    theta_self = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));
    phi_self = (flt_doub *)malloc(sizeof(flt_doub) * pow(ang_res,2));
    rt_self = (complex_double *)malloc(sizeof(complex_double) * pow(ang_res,2));
    rtp_self = (complex_double *)malloc(sizeof(complex_double) * pow(ang_res,2));

    omega_count = 0;
//	#pragma omp parallel for 
    for ( theta_count = 0; theta_count < ang_res; theta_count++){
        for ( phi_count = 0; phi_count  < ang_res; phi_count++){
			omega_count = theta_count * ang_res + phi_count;
            theta_self[omega_count] = theta_count * M_PI / ang_res ;
            phi_self[omega_count] = phi_count * 2.0*M_PI / ang_res ;
            sph_x[omega_count] = cos(phi_self[omega_count]) * sin(theta_self[omega_count]);
            sph_y[omega_count] = sin(phi_self[omega_count]) * sin(theta_self[omega_count]) ;
            sph_z[omega_count] = cos(theta_self[omega_count]) ;
        }
    }


    del_theta = theta_self[ang_res] - theta_self[0];
    del_phi = phi_self[1] - phi_self[0];
    int tiss_num;

    for (tiss_num = 0; tiss_num < phan->no_tiss; tiss_num++){
		for(cnt=0; cnt<(nL+1)*(nL+1); cnt++){
			for(cntp=0; cntp<(nL+1)*(nL+1); cntp++){
				 diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num] = 0;

	}}}

	if(geom->self_sub_vox){
    for (tiss_num = 1; tiss_num < phan->no_tiss; tiss_num++){
              r_ind_self = 0;
#if 1
             for (ii_self = -geom->self_sub_vox/2.0 +0.5; ii_self <= geom->self_sub_vox/2.0 -0.5; ii_self++) {
                for (jj_self = -geom->self_sub_vox/2.0 +0.5; jj_self <= geom->self_sub_vox/2.0-0.5; jj_self++) {
                   for (kk_self = -geom->self_sub_vox/2.0 +0.5; kk_self <= geom->self_sub_vox/2.0 -0.5; kk_self++) {
			 		  // #pragma omp parallel for 
                       for ( omega_count = 0; omega_count  < ang_res*ang_res; omega_count++){
                       /* Now find the intersection with the different faces of lines from this point. We have to determine where any line at an angle (theta, phi) from the point in inside the cube intersects the face, given that we already know that one of the intersections coordinates. */
                        if (sph_x[omega_count] != 0.0){
                                for ( face_calc = 0; face_calc <2; face_calc++){
                                    face_x = face_calc ==0 ? face:-face;
                                    face_y = (face_x - ii_self*2.0/geom->self_sub_vox) * sph_y[omega_count]/ sph_x[omega_count] + jj_self*2.0/geom->self_sub_vox;
                                    face_z = (face_x - ii_self*2.0/geom->self_sub_vox) * sph_z[omega_count]/ sph_x[omega_count] + kk_self*2.0/geom->self_sub_vox;
                                    if (face_x <= face && face_x >= -face &&  face_y <= face && face_y >= -face  && face_z <= face && face_z >= -face  && sph_x[omega_count] * face_x >=0){
                                        cube_x = face_x;
                                        cube_y = face_y;
                                        cube_z = face_z;
                                    }
                                }
                                }


                                if(sph_y[omega_count] != 0.0){
                                for ( face_calc = 0; face_calc <2; face_calc++){
                                    face_y = face_calc ==0 ? face:-face;
                                    face_z = (face_y - jj_self*2.0/geom->self_sub_vox) * sph_z[omega_count]/ sph_y[omega_count] + kk_self*2.0/geom->self_sub_vox;
                                    face_x = (face_y - jj_self*2.0/geom->self_sub_vox) * sph_x[omega_count]/ sph_y[omega_count] + ii_self*2.0/geom->self_sub_vox;
 if (face_x <= face && face_x >= -face &&  face_y <= face && face_y >= -face  && face_z <= face && face_z >= -face && sph_y[omega_count] * face_y >= 0){
                                        cube_x = face_x;
                                        cube_y = face_y;
                                        cube_z = face_z;
                                    }
                                }
                                }

                                if(sph_z[omega_count] !=0.0){
                                for ( face_calc = 0; face_calc <2; face_calc++){
                                    face_z = face_calc ==0 ? face:-face;
                                    face_x = (face_z - kk_self*2.0/geom->self_sub_vox) * sph_x[omega_count]/ sph_z[omega_count] + ii_self*2.0/geom->self_sub_vox;
                                    face_y = (face_z - kk_self*2.0/geom->self_sub_vox) * sph_y[omega_count]/ sph_z[omega_count] + jj_self*2.0/geom->self_sub_vox;
                                    if (face_x <= face && face_x >= -face &&  face_y <= face && face_y >= -face  && face_z <= face && face_z >= -face && sph_z[omega_count] * face_z >=0){
                                        cube_x = face_x;
                                        cube_y = face_y;
                                        cube_z = face_z;
                                    }
                                }
                                }
						    dist_self = sqrt(pow(ii_self*2.0/geom->self_sub_vox - cube_x,2.0) + pow(jj_self*2.0/geom->self_sub_vox - cube_y,2.0) + pow(kk_self*2.0/geom->self_sub_vox - cube_z,2.0)) * geom->delX/2.0;
                            fact_self_vox[r_ind_self][omega_count] = ( 1 - cexp( -(phan->mu_abs[tiss_num] + phan->mu_sc[tiss_num]) * dist_self)) * sin(theta_self[omega_count]);
                            }
                             r_ind_self++;
                            }
                        }
                    }
#endif
			        cnt = 0;
                    for (l = 0; l <= nL; l++) {
                        for (m = -l; m <= l; m++) {
                            cntp = 0;
                            SpherHarmonicArray(l, m, powf(ang_res,2), theta_self, phi_self, rt_self);
                            for (lp = 0; lp <= nL; lp++) {
                                for (mp = -lp; mp <= lp; mp++) {
                                    sub_v_sum_self = 0.0 + 0.0*I;
                                    SpherHarmonicArray(lp, mp, pow(ang_res,2), theta_self, phi_self, rtp_self);
                                    for ( r_ind_self = 0; r_ind_self < pow(geom->self_sub_vox,3); r_ind_self++){
                                        for ( ang_ind = 0; ang_ind < ang_res * ang_res; ang_ind++){
                                            sub_v_sum_self = sub_v_sum_self +  ~(rt_self[ang_ind]) * rtp_self[ang_ind] *  fact_self_vox[r_ind_self][ang_ind];
                                        }
                                    }
                                    diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num] = sub_v_sum_self *del_theta * del_phi / (pow((double)geom->self_sub_vox,3) * (phan->mu_abs[tiss_num] + phan->mu_sc[tiss_num]) * geom->delX * geom->delY * geom->delZ);
								if(cnt == cntp){
									printf("The diagonal term is %e +%e i for tiss = %d, cnt = %d and cntp = %d \n", diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num].real() * geom->delX * geom->delY * geom->delZ, diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num].imag(), tiss_num, cnt, cntp);
								}
                                    cntp++;
                                }
                        }
                        cnt++;
                      }
                    }
#if 0            
        cnt = 0;
                    for (l = 0; l <= nL; l++) {
                        for (m = -l; m <= l; m++) {
                            cntp = 0;
                            SpherHarmonicArray(l, m, powf(ang_res,2), theta_self, phi_self, rt_self);
                            for (lp = 0; lp <= nL; lp++) {
                                for (mp = -lp; mp <= lp; mp++) {
								SpherHarmonicArray(lp, mp, pow(ang_res,2), theta_self, phi_self, rtp_self);
								sub_v_sum_self = 0 + 0*I;
	
						    for ( omega_count = 0; omega_count  < pow(ang_res,2); omega_count++){
                        		   sub_v_sum_self = sub_v_sum_self +  ~(rt_self[omega_count]) * rtp_self[omega_count] * sin(theta_self[omega_count]) * del_theta * del_phi;
                    }
                    diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num] = sub_v_sum_self ;
  printf("The diagonal term is %e +%e i for tiss = %d, cnt = %d and cntp = %d \n", diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num].real(), diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num].imag(), tiss_num, cnt, cntp);
			cntp++;
		}}
 	cnt++;
	}}
#endif
        }
	}
}

void generate_sph_harm_terms() /* Currently this code is for equal values of THETA_ANG_RES and PHI_ANG_RES only */
{
    int l,m, i,j;
    double *theta, *phi;
    complex_double *out;
    int cnt;

    theta = (double *) malloc ( sizeof(double)*THETA_ANG_RES);
    phi = (double *) malloc ( sizeof(double)*PHI_ANG_RES);
    out = (complex_double *) malloc(sizeof(complex_double) * PHI_ANG_RES);

    sph_harm = (complex_double *) malloc (sizeof(complex_double )*(nL+1)*(nL+1) * THETA_ANG_RES * PHI_ANG_RES);

    cnt = 0;
    for ( l = 0; l <= nL; l++){
        for ( m = -l; m <=l ; m++){
            for ( i = 0; i < THETA_ANG_RES; i++){
                for ( j = 0; j < PHI_ANG_RES; j++){
                    theta[j]= (i * M_PI )/ (THETA_ANG_RES-1);
                    phi[j] = (j * 2.0* M_PI) / (PHI_ANG_RES-1);
                }
                SpherHarmonicArray(l, m, PHI_ANG_RES, theta, phi, out);
                memcpy(&sph_harm[cnt*THETA_ANG_RES*PHI_ANG_RES + i*PHI_ANG_RES],out,PHI_ANG_RES*sizeof(complex_double));
            }
            cnt++;
        }
    }
}


#if 0
void generate_diag_terms_2() {

    int ang_res = ANG_RES;
    int theta_count, phi_count, omega_count;
    flt_doub *theta_self, *phi_self;
    float cube_x, cube_y, cube_z;
    int face_calc;
    flt_doub *sph_x, *sph_y, *sph_z, face_x, face_y, face_z;
    flt_doub dist_self;
    flt_doub del_phi, del_theta;
    int r_ind_self;
    int ang_ind ;
    complex_double fact_self_vox;
    complex_double sub_v_sum_self;
    float ii_self, jj_self, kk_self;
    int face = 1; // Defines the face of the cube to be considered later.
    int l,m,lp,mp,cnt,cntp;
    flt_doub cm;

    int i;
    cm = C / phan->n;


    diag_terms_host = (complex_double *)malloc(sizeof(complex_double )* MAX_TISS_NUM*MAX_NL*MAX_NL*MAX_NL*MAX_NL);
    


    sph_x = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));
    sph_y = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));
    sph_z = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));

    theta_self = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));
    phi_self = (flt_doub *)malloc(sizeof(flt_doub) * pow(ang_res,2));

    omega_count = 0;
    for ( theta_count = 0; theta_count < ang_res; theta_count++){
        for ( phi_count = 0; phi_count  < ang_res; phi_count++){
            theta_self[omega_count] = theta_count * M_PI /(ang_res-1);
            phi_self[omega_count] = phi_count * 2.0*M_PI / (ang_res-1);
//            sph_x[omega_count] = cos(phi_self[omega_count]) * sin(theta_self[omega_count]);
//            sph_y[omega_count] = sin(phi_self[omega_count]) * sin(theta_self[omega_count]) ;
//            sph_z[omega_count] = cos(theta_self[omega_count]) ;
            omega_count++;
        }
    }


    del_theta = theta_self[ang_res] - theta_self[0];
    del_phi = phi_self[1] - phi_self[0];
    int tiss_num;

    for (tiss_num = 1; tiss_num < phan->no_tiss; tiss_num++){
		     cnt = 0;
             for (l = 0; l <= nL; l++) {
               for (m = -l; m <= l; m++) {
                               cntp = 0;
							    for (lp = 0; lp <= nL; lp++) {
	                             for (mp = -lp; mp <= lp; mp++) {
                r_ind_self = 0;
        	    sub_v_sum_self = 0.0 + 0.0*I;
/*
                for (ii_self = -geom->self_sub_vox/2.0 +0.5; ii_self <= geom->self_sub_vox/2.0 -0.5; ii_self++) {
                  for (jj_self = -geom->self_sub_vox/2.0 +0.5; jj_self <= geom->self_sub_vox/2.0-0.5; jj_self++) {
                   for (kk_self = -geom->self_sub_vox/2.0 +0.5; kk_self <= geom->self_sub_vox/2.0 -0.5; kk_self++) {
                       for ( omega_count = 0; omega_count  < pow(ang_res,2); omega_count++){
                        if (sph_x[omega_count] != 0.0){
                                for ( face_calc = 0; face_calc <2; face_calc++){
                                    face_x = face_calc ==0 ? face:-face;
                                    face_y = (face_x - ii_self*2.0/geom->self_sub_vox) * sph_y[omega_count]/ sph_x[omega_count] + jj_self*2.0/geom->self_sub_vox;
                                    face_z = (face_x - ii_self*2.0/geom->self_sub_vox) * sph_z[omega_count]/ sph_x[omega_count] + kk_self*2.0/geom->self_sub_vox;
                                    if (face_x <= face && face_x >= -face &&  face_y <= face && face_y >= -face  && face_z <= face && face_z >= -face  && sph_x[omega_count] * face_x >=0){
                                        cube_x = face_x;
                                        cube_y = face_y;
                                        cube_z = face_z;
                                    }
                                }
                                }


                                if(sph_y[omega_count] != 0.0){
                                for ( face_calc = 0; face_calc <2; face_calc++){
                                    face_y = face_calc ==0 ? face:-face;
                                    face_z = (face_y - jj_self*2.0/geom->self_sub_vox) * sph_z[omega_count]/ sph_y[omega_count] + kk_self*2.0/geom->self_sub_vox;
                                    face_x = (face_y - jj_self*2.0/geom->self_sub_vox) * sph_x[omega_count]/ sph_y[omega_count] + ii_self*2.0/geom->self_sub_vox;
 if (face_x <= face && face_x >= -face &&  face_y <= face && face_y >= -face  && face_z <= face && face_z >= -face && sph_y[omega_count] * face_y >= 0){
                                        cube_x = face_x;
                                        cube_y = face_y;
                                        cube_z = face_z;
                                    }
                                }
                                }

                                if(sph_z[omega_count] !=0.0){
                                for ( face_calc = 0; face_calc <2; face_calc++){
                                    face_z = face_calc ==0 ? face:-face;
                                    face_x = (face_z - kk_self*2.0/geom->self_sub_vox) * sph_x[omega_count]/ sph_z[omega_count] + ii_self*2.0/geom->self_sub_vox;
                                    face_y = (face_z - kk_self*2.0/geom->self_sub_vox) * sph_y[omega_count]/ sph_z[omega_count] + jj_self*2.0/geom->self_sub_vox;
                                    if (face_x <= face && face_x >= -face &&  face_y <= face && face_y >= -face  && face_z <= face && face_z >= -face && sph_z[omega_count] * face_z >=0){
                                        cube_x = face_x;
                                        cube_y = face_y;
                                        cube_z = face_z;
                                    }
                                }
                                }
						    dist_self = sqrt(pow(ii_self*2.0/geom->self_sub_vox - cube_x,2.0) + pow(jj_self*2.0/geom->self_sub_vox - cube_y,2.0) + pow(kk_self*2.0/geom->self_sub_vox - cube_z,2.0)) * geom->delX/2.0;
                            fact_self_vox =  ( 1 - cexp( -(phan->mu_abs[tiss_num] + phan->mu_sc[tiss_num]) * dist_self))  *  sin(theta_self[omega_count]) ;
                            sub_v_sum_self = sub_v_sum_self +  ~(sph_harm[SPH_HARM_IND(cnt,theta_self[omega_count],phi_self[omega_count])]) * sph_harm[SPH_HARM_IND(cntp,theta_self[omega_count], phi_self[omega_count])]  * fact_self_vox ;
                            }
                            r_ind_self++;
                            }

                        }
                    }
                       diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num] = sub_v_sum_self *del_theta * del_phi / (cm * pow((double)geom->self_sub_vox,3) * (phan->mu_abs[tiss_num] + phan->mu_sc[tiss_num]) * geom->delX * geom->delY * geom->delZ);
*/
                    for ( omega_count = 0; omega_count  < pow(ang_res,2); omega_count++){
                           sub_v_sum_self = sub_v_sum_self + ~(sph_harm[SPH_HARM_IND(cnt,theta_self[omega_count],phi_self[omega_count])]) * sph_harm[SPH_HARM_IND(cntp,theta_self[omega_count], phi_self[omega_count])] * sin(theta_self[omega_count]) * del_theta * del_phi;
					}
                    diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num] = sub_v_sum_self ;
					printf("The orthonormality using the 2nd method is %e +%e i for tiss = %d, cnt = %d and cntp = %d \n", diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num].real(), diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num].imag(), tiss_num, cnt, cntp);
					cntp++;
					}
					}
					cnt++;
					}
					}
        }
}
#endif

void PropScatmu1(Geometry *geom,Phantom *phan,int nL,complex_double *src) {
    int i;
    int cnt,l;

     for ( i = 0; i < geom->no_vox; i++){
     	for ( cnt = 0; cnt < pow(nL+1,2); cnt ++){
			l = (int) sqrtf(cnt);
            src[VOX_TO_SPIND(i,cnt,(nL+1)*(nL+1))] = pow(phan->g,(double)l) * 1.0 * src[VOX_TO_SPIND(i,cnt,(nL+1)*(nL+1))];
        }
    }

}


void prop_scat(complex_double *src, complex_double *out) {
    int i,l;
    int cnt;


    for ( i = 0; i < geom->no_vox; i++){
	  for ( cnt = 0; cnt < pow(nL+1,2); cnt++){
	       l = (int) sqrtf(cnt);
           out[VOX_TO_SPIND(i,cnt,(nL+1)*(nL+1))] =pow(phan->g,(double)l) * phan->mu_sc[phan->tiss_type[i]]*src[VOX_TO_SPIND(i,cnt,(nL+1)*(nL+1))];
        }
    }

}

complex_double ComputePhotonNumber(int num,complex_double *W) {
    complex_double ret;
    complex_double **tmp;
    double *theta;
    double *phi;
    int l,m,i,j,k;
    int ip,jp;
    int cnt;
    int nCoef;
    double posX,posY,diam,mag,rad;
    complex_double val,val_tmp;
    double delphi,deltheta;

    ret = 0.0 + 0*I;

    nCoef = (nL+1)*(nL+1);
    tmp = (complex_double **) malloc(sizeof(complex_double*)*nCoef);
    for (m=0;m<nCoef;m++) {
        tmp[m] = (complex_double *) malloc(sizeof(complex_double)*num*num);
    }

    theta = (double *)malloc(sizeof(double) * num*num);
    phi = (double *) malloc(sizeof(double)*num*num);

  // Bottom Face
    deltheta = (M_PI/2.0)/(double)(num);
    delphi = (2.0*M_PI)/(double)(num);
    for (i=0;i<num;i++) {
        for (j=0;j<num;j++) {
            theta[i+j*num] = (double)i*(M_PI/2.0)/(double)(num);
            phi[i+j*num] = (double)j*(M_PI*2.0)/(double)(num);
        }
    }
    complex_double tmp_var;
    int r_ind,r_indp;
	int lp,mp;


    cnt = 0;
    for (l=0;l<=nL;l++) {
        for (m=-l;m<=l;m++) {
            SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
            cnt++;
		}
	 }
    i = geom->bounZ + geom->nZ - 1;

    for (j=geom->bounY;j<geom->nY+geom->bounY;j++) {
        for (k=geom->bounX;k<geom->nX+geom->bounX;k++) {
		    val = 0.0 + 0.0*I ;
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;

            rad = sqrt((geom->xs[k] - beam_src->pos_x) * (geom->xs[k] - beam_src->pos_x) +
                (geom->ys[j] - beam_src->pos_y) * (geom->ys[j] - beam_src->pos_y));
          if (rad <= beam_src->diam/2.0) {
//            if (k == geom->nX/2 + geom->bounX && j == geom->nY/2 + geom->bounY ) {
			    tmp_var = beam_src->mag*1.0 + 0.0*I;
                for(ip=geom->bounZ ; ip < geom->nZ + geom->bounZ; ip++){
                     r_indp = ip* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
                     tmp_var = tmp_var * cexp(-(phan->mu_sc[phan->tiss_type[r_indp]] +phan->mu_abs[phan->tiss_type[r_indp]]) * geom->delZ);
                }

                printf("%e is the value of the 0 term for voxel (%d %d %d) with r_ind as %d \n", tmp_var.real(),i,j,k,r_ind);
                val = val+ tmp_var;  	
            }

            for (m=0;m<nCoef;m++) {
			   val_tmp = 0 + 0*I;
	           for (lp = 0;lp < num;lp++) {
                for (mp=0;mp<num;mp++) {
                        val_tmp = val_tmp + W[VOX_TO_SPIND(r_ind,m,nCoef)]*tmp[m][lp+mp*num] * sin(theta[lp+mp*num])*deltheta *delphi*cos(theta[lp+mp*num])*geom->delX * geom->delY * C/phan->n;
                    }
                }
	            val = val_tmp + val;
            }
            ret = ret + val;
        }
    }

	complex_double ret_tmp;
    ret_tmp = ret;
	printf("Contribution from bottom face is %e \n", ret_tmp.real());

    // Top Face
    deltheta = (M_PI/2.0)/(double)(num);
    delphi = (2.0*M_PI)/(double)(num);
    for (i=0;i<num;i++) {
        for (j=0;j<num;j++) {
            theta[i+j*num] = (double)i*(M_PI/2.0)/(double)(num-1) + (M_PI/2.0);
            phi[i+j*num] = (double)j*(M_PI*2.0)/(double)(num-1);
        }
    }

   cnt = 0;
    for (l=0;l<=nL;l++) {
        for (m=-l;m<=l;m++) {
            SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
            cnt++;
        }
    }

    i = geom->bounZ;
    for (j=geom->bounY;j<geom->nY+geom->bounY;j++) {
        for (k=geom->bounX;k<geom->nX+geom->bounX;k++) {
            val = 0 + 0*I;
			r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
            for (m=0;m<nCoef;m++) {
				val_tmp = 0 + 0*I;
            	for (ip = 0;ip < num;ip++) {
                	for (jp=0;jp<num;jp++) {
                        val_tmp = val_tmp + W[VOX_TO_SPIND(r_ind,m,nCoef)] * tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi * (-1) * cos(theta[ip+jp*num]) * geom->delX * geom->delY  * C/phan->n;
                    }
                }
                val = val_tmp + val;
            }
            ret = ret + val;
        }
    }
    
	printf("Contribution from top face is %e \n", ret - ret_tmp );
	ret_tmp = ret;

  // Side Face 1
    deltheta = (M_PI)/(double)(num);
    delphi = (M_PI)/(double)(num);
    for (i=0;i<num;i++) {
        for (j=0;j<num;j++) {
            theta[i+j*num] = (double)i*(M_PI)/(double)(num);
            phi[i+j*num] = (double)j*(M_PI)/(double)(num)+ (M_PI/2);
        }
    }

    cnt = 0;
    for (l=0;l<=nL;l++) {
        for (m=-l;m<=l;m++) {
            SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
            cnt++;
        }
    }

    k = geom->bounX;
    for (i=geom->bounZ;i<geom->nZ+geom->bounZ;i++) {
   		for (j=geom->bounY;j<geom->nY+geom->bounY;j++) {
            val = 0 + 0*I;
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
            for (ip = 0;ip < num;ip++) {
                for (jp=0;jp<num;jp++) {
                    for (m=0;m<nCoef;m++) {
				       val = val + W[VOX_TO_SPIND(r_ind,m,nCoef)] * tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi * sin(theta[ip+jp*num]) * cos(phi[ip+jp*num]) * (-1) * geom->delX * geom->delY  * C/phan->n;
                    }
                }
            }
            ret = ret + val;
        }
    }

	printf("Contribution from side face is %e \n", ret - ret_tmp );
	ret_tmp = ret;

    // Side Face 2
    deltheta = (M_PI)/(double)(num);
    delphi = (M_PI)/(double)(num);
    for (i=0;i<num;i++) {
        for (j=0;j<num;j++) {
            theta[i+j*num] = (double)i*(M_PI)/(double)(num);
            phi[i+j*num] = (double)j*(M_PI)/(double)(num)- (M_PI/2);
        }
    }

    cnt = 0;
    for (l=0;l<=nL;l++) {
        for (m=-l;m<=l;m++) {
            SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
            cnt++;
        }
    }

    k = geom->nX-1 + geom->bounX;
    for (i=geom->bounZ;i<geom->nZ+geom->bounZ;i++) {
        for (j=geom->bounY;j<geom->nY+geom->bounY;j++) {
            val = 0 + 0*I;
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;

            for (ip = 0;ip < num;ip++) {
                for (jp=0;jp<num;jp++) {
                    for (m=0;m<nCoef;m++) {
                        val = val + W[VOX_TO_SPIND(r_ind,m,nCoef)]*tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi* sin(theta[ip+jp*num]) * cos(phi[ip+jp*num]) * geom->delX * geom->delY  * C/phan->n ;
                    }
                }
            }
            ret = ret + val;
        }
    }

	printf("Contribution from side face 2 is %e \n", ret - ret_tmp );
	ret_tmp = ret;

    // Side Face 3
    deltheta = (M_PI)/(double)(num);
    delphi = (M_PI)/(double)(num);
    for (i=0;i<num;i++) {
        for (j=0;j<num;j++) {
            theta[i+j*num] = (double)i*(M_PI)/(double)(num);
            phi[i+j*num] = (double)j*(M_PI)/(double)(num)+ (M_PI);
        }
    }

    cnt = 0;
    for (l=0;l<=nL;l++) {
        for (m=-l;m<=l;m++) {
            SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
            cnt++;
        }
    }

    j = geom->bounY;
    for (i=geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
        for (k=geom->bounX;k<geom->nX+geom->bounX;k++) {
            val = 0 + 0*I;
			r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;
            for (ip = 0;ip < num;ip++) {
                for (jp=0;jp<num;jp++) {
                    for (m=0;m<nCoef;m++) {
                        val = val + W[VOX_TO_SPIND(r_ind,m,nCoef)] *
                            tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi * sin(theta[ip+jp*num]) * (-1) * sin(phi[ip+jp*num]) * geom->delX * geom->delY  * C/phan->n;
                    }
                }
            }
            ret = ret + val; 
        }
    }
	printf("Contribution from side face 3 is %e \n", ret - ret_tmp );
	ret_tmp = ret;

   // Side Face 4
    deltheta = (M_PI)/(double)(num);
    delphi = (M_PI)/(double)(num);
    for (i=0;i<num;i++) {
        for (j=0;j<num;j++) {
            theta[i+j*num] = (double)i*(M_PI)/(double)(num);
            phi[i+j*num] = (double)j*(M_PI)/(double)(num);
        }
    }

    cnt = 0;
    for (l=0;l<=nL;l++) {
        for (m=-l;m<=l;m++) {
            SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
            cnt++;
        }
    }

    j = geom->nY - 1 + geom->bounY;
    for (i=geom->bounZ;i<geom->nZ + geom->bounZ;i++) {
        for (k=geom->bounX;k<geom->nX+geom->bounX;k++) {
            val = 0 + 0*I;
            r_ind = i* (geom->nX + 2*geom->bounX )* (geom->nY + 2*geom->bounY) + j* (geom->nX + 2*geom->bounX) + k;

            val = 0 + 0*I;
            for (ip = 0;ip < num;ip++) {
                for (jp=0;jp<num;jp++) {
                    for (m=0;m<nCoef;m++) {
					     val = val + W[VOX_TO_SPIND(r_ind,m,nCoef)] *
                            tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi *  sin(theta[ip+jp*num]) * sin(phi[ip+jp*num]) * geom->delX * geom->delY  * C/phan->n ;
                    }
                }
            }
            ret = ret + val;
        }
    }
	
	printf("Contribution from side face 4 is %e \n", ret - ret_tmp );
	ret_tmp = ret;

    // Clean up
    for (m=0;m<nCoef;m++){
        free(tmp[m]);
    }
    free(tmp);
    free(theta);
    free(phi);


    return(ret);
}

