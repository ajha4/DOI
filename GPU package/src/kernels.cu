#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "rte.h"
#include <pthread.h>

#include "complex_arith.cu"

extern Geometry *geom;
extern Phantom *phan;
extern Source *beam_src;
extern complex_double *diag_terms_host;
extern complex_double *sph_harm;
extern Info_Stat *info_stat_host;
extern SHORT nL;
extern int nTerms;
__constant__ Info_Stat info_stat;
__constant__ SHORT nL_dev;
__constant__ complex_double diag_terms_dev[MAX_TISS_NUM*MAX_NL*MAX_NL*MAX_NL*MAX_NL];


/* Just think of the vox indices as (blk_dep*BLK_SIZE + dep + info_stat.bounZ, blk_row*BLK_SIZE + row + info_stat.bounY, blk_col * BLK_SIZE + col + info_stat.nX) */
__device__ int get_voxind_spind(int blk_dep, int blk_row,int blk_col, int dep, int row, int col, int cnt){

     return VOX_TO_SPIND( ((blk_dep + dep + info_stat.bounZ)*(2*info_stat.bounX + info_stat.nX) * (2*info_stat.bounY + info_stat.nY) + (blk_row*BLK_SIZE + row + info_stat.bounY)* (2*info_stat.bounX + info_stat.nX) + (blk_col * BLK_SIZE + col + info_stat.bounX)), cnt, (nL_dev+1)*(nL_dev+1));

}

__device__ int get_vind_phanind(int dep, int row, int col){

  return ((info_stat.bounZ + dep) * (info_stat.nX + 2 * info_stat.bounX ) * (info_stat.nY + 2 * info_stat.bounY ) /* reached the correct layer */ + ( info_stat.bounY + row)* (info_stat.nX + 2 * info_stat.bounX ) + (info_stat.bounX + col));

}

__device__ int get_voxind_phanind(int blk_dep, int blk_row,int blk_col,int blk_offset_i, int blk_offset_j, int blk_offset_k, int dep, int row, int col){
      
	return (info_stat.bounZ  + blk_dep + blk_offset_i * BLK_SIZE+ dep) * (info_stat.nX + 2 * info_stat.bounX ) * (info_stat.nY + 2 * info_stat.bounY ) /* reached the correct layer */ + (blk_row * BLK_SIZE + info_stat.bounY + blk_offset_j * BLK_SIZE + row)* (info_stat.nX + 2 * info_stat.bounX ) + (info_stat.bounX + blk_col * BLK_SIZE + blk_offset_k * BLK_SIZE + col);
}

__device__ int get_ind_phanind(int dep, int row, int col){

   return dep* (2*info_stat.boun_blk_sizeX + BLK_SIZE)* (BLK_SIZE + 2*info_stat.boun_blk_sizeY) + row* (BLK_SIZE + 2*info_stat.boun_blk_sizeX) + col;

}


__device__ complex_double
mu_int (COR cor_src, COR cor_dest, float r, Info_Dyn info_dyn)
{
  complex_double mu_tot_path;
  float alpha_x_curr, alpha_y_curr, alpha_z_curr;
  float alpha_xinc, alpha_yinc, alpha_zinc;
  SHORT i_curr, j_curr, k_curr;
  float alpha_prev;
  byte flag = 1;

  alpha_x_curr = 0;
  alpha_y_curr = 0;
  alpha_z_curr = 0;

  i_curr = cor_src.i;
  j_curr = cor_src.j;
  k_curr = cor_src.k;

  alpha_prev = 0;
  mu_tot_path = 0.0 + 0.0*I;
  SHORT x_inc_end,y_inc_end,z_inc_end;



  if (cor_dest.i != cor_src.i)
    alpha_zinc = fabs (0.5 / (cor_dest.i - cor_src.i)); // The initial increment along the z axis. 
  else
    alpha_zinc = INF;
  if (cor_dest.j != cor_src.j)
    alpha_yinc = fabs (0.5 / (cor_dest.j - cor_src.j));
  else
    alpha_yinc = INF;
  if (cor_dest.k != cor_src.k)
    alpha_xinc = fabs (0.5 / (cor_dest.k - cor_src.k));
  else
    alpha_xinc = INF;

  while (flag == 1 && alpha_prev < 1) // Hack for now to avoid infinite loops
    {

    if (alpha_z_curr + alpha_zinc  <= alpha_x_curr + alpha_xinc  && alpha_z_curr + alpha_zinc <= alpha_y_curr + alpha_yinc)
    {

      alpha_z_curr += alpha_zinc;
      if ( i_curr == cor_src.i)
            alpha_zinc *=2; // We have taken the first step along the z axis, which was half a voxel. Now every step will be one voxel. 

      mu_tot_path = mu_tot_path + info_stat.mu_tot[info_dyn.tiss_type[get_vind_phanind(i_curr,j_curr,k_curr)]] * r * (alpha_z_curr -  alpha_prev);
      i_curr = (cor_src.i < cor_dest.i) ? i_curr + 1 : i_curr - 1;

      alpha_prev = alpha_z_curr;
	 if (i_curr == cor_dest.i && j_curr == cor_dest.j && k_curr == cor_dest.k)
        {
          mu_tot_path = mu_tot_path + info_stat.mu_tot[info_dyn.tiss_type[get_vind_phanind(i_curr,j_curr,k_curr)]] * r * alpha_zinc/2.0;
          flag = 0;
          return mu_tot_path;
        }
    }

    else if (alpha_y_curr + alpha_yinc  < alpha_z_curr + alpha_zinc  && alpha_y_curr + alpha_yinc <= alpha_x_curr + alpha_xinc )
    {
      alpha_y_curr += alpha_yinc;

      if ( j_curr == cor_src.j)
            alpha_yinc *=2; 

      mu_tot_path = mu_tot_path + info_stat.mu_tot[info_dyn.tiss_type[get_vind_phanind(i_curr,j_curr,k_curr)]] * r * (alpha_y_curr -  alpha_prev);
      j_curr = (cor_src.j < cor_dest.j) ? j_curr + 1 : j_curr - 1;

      alpha_prev = alpha_y_curr;
      if (i_curr == cor_dest.i && j_curr == cor_dest.j && k_curr == cor_dest.k)
        {
          mu_tot_path = mu_tot_path + info_stat.mu_tot[info_dyn.tiss_type[get_vind_phanind(i_curr,j_curr,k_curr)]] * r * alpha_yinc/2.0;
          flag = 0;
          return mu_tot_path;
        }
    }
    else if (alpha_x_curr + alpha_xinc < alpha_y_curr + alpha_yinc && alpha_x_curr + alpha_xinc < alpha_z_curr + alpha_zinc )
    {
      alpha_x_curr += alpha_xinc;
      if ( k_curr == cor_src.k)
           alpha_xinc *=2; 

      mu_tot_path = mu_tot_path +  (info_stat.mu_tot[info_dyn.tiss_type[get_vind_phanind(i_curr,j_curr,k_curr)]]) * r * (alpha_x_curr -  alpha_prev);
      k_curr = (cor_src.k < cor_dest.k) ? k_curr + 1 : k_curr - 1;
      alpha_prev = alpha_x_curr;
      if (i_curr == cor_dest.i && j_curr == cor_dest.j && k_curr == cor_dest.k)
        {
          mu_tot_path = mu_tot_path + info_stat.mu_tot[info_dyn.tiss_type[get_vind_phanind(i_curr,j_curr,k_curr)]] * r * alpha_xinc/2.0;
          flag = 0;
          return mu_tot_path;
        }
    }

    }
  return mu_tot_path;
}


__global__ void compute_diagonal_abs (doublecomplex *src_dist,doublecomplex *out_dist,Info_Dyn info_dyn, int cnt,SHORT block_dep, SHORT layer_start, SHORT flag)
{
  int block_row = blockIdx.y;
  int block_col = blockIdx.x;
  int i = threadIdx.z;
  int j = threadIdx.y;
  int k = threadIdx.x;
  int sp_ind, sp_ind_src, r_ind;
  SHORT cntp;

  int nL_tmp;
  flt_doub cm_tmp = 1.0;

  nL_tmp = (flag == 0) ? 0 : nL_dev;

//  cm_tmp = (flag == 0) ? 1 : info_stat.cm

  r_ind = get_vind_phanind(block_dep + layer_start + i, block_row * BLK_SIZE + j, block_col * BLK_SIZE+ k);

   for(cnt=0; cnt < (nL_tmp+1)*(nL_tmp+1); cnt++){
  	sp_ind = get_voxind_spind(block_dep - info_stat.bounZ, block_row, block_col ,i,j,k,cnt);
    out_dist[sp_ind] = 0 + 0*I;
   	for(cntp=0; cntp<(nL_tmp+1)*(nL_tmp+1) ; cntp++){
  			sp_ind_src = get_voxind_spind(block_dep + layer_start, block_row, block_col,i,j,k,cntp);
            out_dist[sp_ind] =  out_dist[sp_ind] + (1.0/cm_tmp)*diag_terms_dev[cnt * (MAX_NL) * MAX_NL * MAX_TISS_NUM + cntp * MAX_TISS_NUM + info_dyn.tiss_type[r_ind]] * src_dist[sp_ind_src];
	}
   }
}

#if 0
//extern __shared__ char array_tisstype[];
__global__ void compute_subvox_2(complex_double* src_dist_dev,complex_double* out_dist_dev, Info_Dyn info_dyn, complex_double* sph_harm_dev, SHORT cnt, COR subvox_src,COR subvox_dest, flt_doub dz_sub, flt_doub dy_sub, flt_doub dx_sub, SHORT blk_dep, SHORT start_layer){


  SHORT block_row = blockIdx.y;
  SHORT block_col = blockIdx.x;
  SHORT i = threadIdx.z;
  SHORT j = threadIdx.y;
  SHORT k = threadIdx.x;

  SHORT cntp;
  complex_double out_tmp, tmp;
  COR cor_src, cor_dest;

  SHORT ip, jp, kp;
  int sp_ind;
  flt_doub theta,phi,dist,dx,dy,dz;


  int sp_ind_src;

#if 0
  cor_src.i = i  + info_stat.boun_blk_sizeZ;
  cor_src.j = j  + info_stat.boun_blk_sizeY;
  cor_src.k = k  + info_stat.boun_blk_sizeX;

#else
  cor_src.i = i + blk_dep + start_layer;
  cor_src.j = j  + block_row*BLK_SIZE;
  cor_src.k = k  + block_col*BLK_SIZE;
#endif

//  __shared__ complex_double src_tmp[BLK_SRC_SIZE];


#if 0
  int blk_offset_i, blk_offset_j,blk_offset_k;
   byte *tisstype_tmp = (byte *) array_tisstype;

    for(blk_offset_i=0; blk_offset_i< 1 + 2*info_stat.boun_blkZ; blk_offset_i++){
        for(blk_offset_j=0; blk_offset_j< 1 + 2*info_stat.boun_blkY; blk_offset_j++){
            for(blk_offset_k=0; blk_offset_k< 1 + 2*info_stat.boun_blkX; blk_offset_k++){

                    tisstype_tmp[get_ind_phanind(blk_offset_i*BLK_SIZE + i, blk_offset_j*BLK_SIZE + j, blk_offset_k*BLK_SIZE + k)] = info_dyn.tiss_type[get_voxind_phanind(blk_dep, block_row, block_col,blk_offset_i - info_stat.boun_blkZ,blk_offset_j - info_stat.boun_blkY,blk_offset_k - info_stat.boun_blkX,i,j,k)];

//                  src_tmp[get_ind_phanind(blk_offset_i*BLK_SIZE + i, blk_offset_j*BLK_SIZE + j, blk_offset_k*BLK_SIZE + k)] = src_dist_dev[VOX_TO_SPIND(get_voxind_phanind(blk_dep, block_row, block_col,(blk_offset_i - info_stat.bounZ/BLK_SIZE),(blk_offset_j - info_stat.bounY/BLK_SIZE),(blk_offset_k - info_stat.bounX/BLK_SIZE),i,j,k),cntp,info_stat.no_vox)];

            }
        }
    }
    __syncthreads();
#endif

  out_tmp = 0 + 0*I;
  flt_doub sub_dist;
 
 
  for(ip= i - info_stat.subbounZ; ip <= i + info_stat.subbounZ; ip++){
	dz = -(ip-i)*info_stat.delZ + dz_sub;
  	for(jp= j - info_stat.subbounY; jp <= j + info_stat.subbounY; jp++){
		dy = -(jp-j)*info_stat.delY + dy_sub;
		  for(kp= k - info_stat.subbounX; kp <= k + info_stat.subbounX; kp++){
			dx = -(kp-k)*info_stat.delX + dx_sub;

            dist = sqrt((i-ip)*(i-ip)*info_stat.delZ*info_stat.delZ + (j-jp)*(j-jp)*info_stat.delY*info_stat.delY + (k-kp)*(k-kp)*info_stat.delX*info_stat.delX);

			if( dist <= info_stat.sub_thresh && ( i != ip || j != jp || k != kp)){
            
			  sub_dist =  sqrt(dx*dx + dy*dy + dz*dz);

#if 0
              cor_dest.i = ip + info_stat.boun_blk_sizeZ;
              cor_dest.j = jp + info_stat.boun_blk_sizeY;
              cor_dest.k = kp + info_stat.boun_blk_sizeX;
#else
              cor_dest.i = ip +blk_dep + start_layer;
              cor_dest.j = jp + block_row*BLK_SIZE;
              cor_dest.k = kp +  block_col*BLK_SIZE;
#endif


              theta = atan(sqrt(dx*dx + dy*dy)/dz );
              phi = atan2(dy,dx);

   			  for(cnt=0; cnt < (nL_dev+1)*(nL_dev+1); cnt++){
  			    sp_ind = get_voxind_spind(blk_dep - info_stat.bounZ, block_row, block_col,i,j,k,cnt);

                tmp = 0 + 0*I;
              	for(cntp=0; cntp< (nL_dev+1)*(nL_dev+1); cntp++){
                	sp_ind_src = get_voxind_spind(blk_dep + start_layer, block_row, block_col,ip,jp,kp,cntp);
	                tmp = tmp + src_dist_dev[sp_ind_src] * sph_harm_dev[SPH_HARM_IND(cntp,theta,phi)];
    	        }
				tmp = tmp * ~sph_harm_dev[SPH_HARM_IND(cnt,theta,phi)];
        	    tmp = tmp* cexp_dev(-mu_int(cor_src, cor_dest, subvox_src, subvox_dest,sub_dist, info_dyn )); ;//cexp_dev(-(1.01+0.0*I)*sub_dist);// cexp_dev(-mu_int(cor_src, cor_dest, subvox_src, subvox_dest,sub_dist, tisstype_tmp,info_dyn )); //cexp_dev(-(1.01+0.0*I)*dist)
            	tmp = tmp * (1.0/(  info_stat.cm * sub_dist*sub_dist *  __powf(__int2float_ru(info_stat.sub_vox),6.0)));
    			out_dist_dev[sp_ind] =  out_dist_dev[sp_ind] + tmp;
	              //out_tmp = out_tmp + tmp;
    	      }
			 }
            }
        }
      }
    __syncthreads();

}

#endif              
			  
__global__ void compute_propabs(complex_double *src_dist_dev, complex_double *out_dist_dev, Info_Dyn info_dyn,complex_double *sph_harm_dev,SHORT cnt, SHORT blk_dep, SHORT start_layer, SHORT flag){ // If flag =1, only then evaluate all the other spherical harmonics, else nl = 0.

  SHORT block_row = blockIdx.y;
  SHORT block_col = blockIdx.x;
  SHORT i = threadIdx.z;
  SHORT j = threadIdx.y;
  SHORT k = threadIdx.x;
  int sp_ind;
  complex_double out_tmp,tmp;
  COR cor_src, cor_dest;
  SHORT cntp;



  flt_doub theta,phi,dist,dx,dy,dz;
  int ip,jp,kp;



  int sp_ind_src; 
  int flag_error=0;
			  
  int nL_tmp;
  flt_doub cm_tmp = 1.0;

  nL_tmp = (flag == 0) ? 0 : nL_dev;

  cor_src.i = i + blk_dep + start_layer;
  cor_src.j = j + block_row*BLK_SIZE; 
  cor_src.k = k + block_col*BLK_SIZE; 

  int sp_ind2;
  for (ip = 0; ip < info_stat.nZ ; ip++){
    dz = -(ip-cor_src.i)*info_stat.delZ;
	  for (jp = 0; jp < info_stat.nY ; jp++){
       dy = -(jp-cor_src.j)*info_stat.delY;
       for (kp = 0; kp < info_stat.nX; kp++){
          dx = -(kp-cor_src.k)*info_stat.delX;

          dist =  sqrt(dx*dx + dy*dy + dz*dz);

          if((ip != cor_src.i || jp != cor_src.j || kp != cor_src.k) && dist < info_stat.prop_thresh){
              theta = acos(sqrt(dx*dx + dy*dy)/dz );

              phi = atan2(dy,dx);
			  if(phi < 0)
					phi = phi + 2*M_PI;

			  cor_dest.i = ip;
			  cor_dest.j = jp;
			  cor_dest.k = kp;


			  for(cnt = 0; cnt < (nL_tmp+1)*(nL_tmp+1); cnt++){

    			  sp_ind = get_voxind_spind(blk_dep - info_stat.bounZ, block_row, block_col,i,j,k,cnt);
				  tmp = 0 + 0*I;
				  for(cntp=0; cntp< (nL_tmp+1)*(nL_tmp+1); cntp++){
				  	  sp_ind_src = get_voxind_spind(0, 0, 0,ip,jp,kp,cntp);
					  if(theta > M_PI || phi > 2*M_PI )
						flag_error = 1;
	        	      tmp = tmp + src_dist_dev[sp_ind_src]* sph_harm_dev[SPH_HARM_IND(cntp,theta,phi)];
				  } 
    	          tmp = tmp * ~sph_harm_dev[SPH_HARM_IND(cnt,theta,phi)];

        	      tmp = tmp * cexp_dev(-mu_int(cor_src, cor_dest, dist, info_dyn )); //cexp_dev(-1*(1.01 + 0*I)*dist);
   		    	  tmp = tmp * (1.0/(cm_tmp*dist*dist));
    
				  out_dist_dev[sp_ind] = out_dist_dev[sp_ind] + tmp;

	              out_tmp = out_tmp + tmp; 
			   }
             }
           }
        }
      }

    __syncthreads();
}

__global__ void scale_dist_dev (doublecomplex *W,double scale_fac,SHORT cnt,SHORT block_dep )
{

  int sp_ind;

  int block_row = blockIdx.y;
  int block_col = blockIdx.x;
  int i = threadIdx.z;
  int j = threadIdx.y;
  int k = threadIdx.x;

	for(cnt=0; cnt < (nL_dev+1)*(nL_dev+1); cnt++){
        sp_ind = get_voxind_spind(block_dep - info_stat.bounZ, block_row , block_col ,i,j,k,cnt);
        W[sp_ind] = W[sp_ind]*scale_fac;
	}
}

__global__ void write_dist_dev (doublecomplex *W,doublecomplex val,SHORT cnt,SHORT block_dep)
{

  int sp_ind;

  int block_row = blockIdx.y;
  int block_col = blockIdx.x;
  int i = threadIdx.z;
  int j = threadIdx.y;
  int k = threadIdx.x;

	for(cnt=0; cnt < (nL_dev+1)*(nL_dev+1); cnt++){
        sp_ind = get_voxind_spind(block_dep - info_stat.bounZ, block_row, block_col,i,j,k,cnt);
        W[sp_ind] = val;
	}
}



__global__ void copy_dist_dev (doublecomplex *W1,doublecomplex *W2)
{

  int sp_ind;
  int cnt=0;
  int block_dep=0;

  int block_row = blockIdx.y;
  int block_col = blockIdx.x;
  int i = threadIdx.z;
  int j = threadIdx.y;
  int k = threadIdx.x;


  for(block_dep = 0; block_dep < info_stat.nZ; block_dep = block_dep + BLK_SIZE){
  	for (cnt = 0; cnt < (nL_dev+1)*(nL_dev+1); cnt++){
		sp_ind = get_voxind_spind(block_dep - info_stat.bounZ, block_row, block_col,i,j,k,cnt);
		 W2[sp_ind] = W1[sp_ind];
		 __syncthreads();
	 }
  }
}


__global__ void add_dist_dev (doublecomplex *W1,doublecomplex *W2, doublecomplex *out ) 
{
  int sp_ind,cnt,block_dep;

  int block_row = blockIdx.y;
  int block_col = blockIdx.x;
  int i = threadIdx.z; 
  int j = threadIdx.y; 
  int k = threadIdx.x; 


  for(block_dep = 0; block_dep < info_stat.nZ; block_dep = block_dep + BLK_SIZE){
 	 for (cnt = 0; cnt < (nL_dev+1)*(nL_dev+1); cnt++){
	  	 sp_ind = get_voxind_spind(block_dep -info_stat.bounZ, block_row, block_col,i,j,k,cnt);
    	 out[sp_ind] = W1[sp_ind] + W2[sp_ind];
	 }
  }
}

__global__ void prop_scat_dev (doublecomplex *src_dist,doublecomplex *out_dist, Info_Dyn info_dyn)
{

  int sp_ind,cnt,block_dep,l,r_ind;

  int block_row = blockIdx.y;
  int block_col = blockIdx.x;
  int i = threadIdx.z;
  int j = threadIdx.y;
  int k = threadIdx.x;


  for(block_dep = 0; block_dep < info_stat.nZ; block_dep = block_dep + BLK_SIZE){
  	for (cnt = 0; cnt < (nL_dev+1)*(nL_dev+1); cnt++){
         sp_ind = get_voxind_spind(block_dep -info_stat.bounZ, block_row, block_col,i,j,k,cnt);
		 r_ind = get_vind_phanind(block_dep + i, block_row * BLK_SIZE + j, block_col * BLK_SIZE +k);
		 l = (int) sqrtf(cnt); 
         out_dist[sp_ind] = pow (__int2double_rn(info_stat.g),__int2double_rn(l))   * src_dist[sp_ind] * info_stat.mu_sc[info_dyn.tiss_type[r_ind]] ;
     }
  }
}

__global__ void compute_sph_coord(flt_doub* theta_self, flt_doub* phi_self, flt_doub* sph_x, flt_doub* sph_y, flt_doub* sph_z, int theta_blk, int phi_blk)
{
    int theta_count = threadIdx.x + theta_blk;
    int phi_count = threadIdx.y + phi_blk;

    int omega_count;

    omega_count = theta_count*ANG_RES + phi_count;

    theta_self[omega_count] = theta_count * M_PI / ANG_RES ;
    phi_self[omega_count] = phi_count * 2.0*M_PI / ANG_RES ;
    sph_x[omega_count] = cos(phi_count * 2.0*M_PI / ANG_RES) * sin(theta_count * M_PI / ANG_RES);
    sph_y[omega_count] = sin(phi_count * 2.0*M_PI / ANG_RES) * sin(theta_count * M_PI / ANG_RES) ;
    sph_z[omega_count] = cos(theta_count * M_PI / ANG_RES) ;
}

__global__ void compute_diag_selfsub(complex_double *fact_self_vox, flt_doub *sph_x, flt_doub *sph_y, flt_doub *sph_z, flt_doub *theta_self, int omega_count, SHORT tiss_num, int blk_dep)

{

    int blk_row = blockIdx.y;
    int blk_col = blockIdx.x; 
    int z_ind = threadIdx.z + blk_dep;
    int y_ind = threadIdx.y +blk_row* BLK_SELF_SUB_VOX;
    int x_ind = threadIdx.x + blk_col* BLK_SELF_SUB_VOX;
    int face_calc;
    int face = 1;
    flt_doub face_x, face_y, face_z, cube_x, cube_y, cube_z, dist_self ;


    //int r_ind_self = (threadIdx.z + blk_dep) * (info_stat.self_sub_vox)*(info_stat.self_sub_vox) + (threadIdx.y +blk_row* BLK_SELF_SUB_VOX)  * info_stat.self_sub_vox + (threadIdx.x + blk_col* BLK_SELF_SUB_VOX);
    int r_ind_self = (z_ind) * (info_stat.self_sub_vox)*(info_stat.self_sub_vox) + (y_ind)  * info_stat.self_sub_vox + (x_ind);

    flt_doub ii_self = -info_stat.self_sub_vox/2.0 +0.5 + z_ind;
    flt_doub jj_self = -info_stat.self_sub_vox/2.0 +0.5 + y_ind;
    flt_doub kk_self = -info_stat.self_sub_vox/2.0 +0.5 + x_ind;

    face_x = 0;
    face_y = 0;
    face_z = 0;
    cube_x = 0;
    cube_y = 0;
    cube_z = 0;

    if (sph_x[omega_count] != 0.0){
       for ( face_calc = 0; face_calc <2; face_calc++){
          face_x = face_calc ==0 ? face:-face;
          face_y = (face_x - ii_self*2.0/info_stat.self_sub_vox) * sph_y[omega_count]/ sph_x[omega_count] + jj_self*2.0/info_stat.self_sub_vox;
          face_z = (face_x - ii_self*2.0/info_stat.self_sub_vox) * sph_z[omega_count]/ sph_x[omega_count] + kk_self*2.0/info_stat.self_sub_vox;
          if (face_x <= face && face_x >= -face &&  face_y <=face  && face_y >= -face  && face_z <= face && face_z >= -face  && sph_x[omega_count] * face_x >=0){
              cube_x = face_x;
              cube_y = face_y;
              cube_z = face_z;
          }
       }
     }

#if 1
     if(sph_y[omega_count] != 0.0){
        for ( face_calc = 0; face_calc <2; face_calc++){
           face_y = face_calc ==0 ? face:-face;
           face_z = (face_y - jj_self*2.0/info_stat.self_sub_vox) * sph_z[omega_count]/ sph_y[omega_count] + kk_self*2.0/info_stat.self_sub_vox;
           face_x = (face_y - jj_self*2.0/info_stat.self_sub_vox) * sph_x[omega_count]/ sph_y[omega_count] + ii_self*2.0/info_stat.self_sub_vox;
           if (face_x <= face && face_x >= -face &&  face_y <= face && face_y >= -face  && face_z <= face && face_z >= -face && sph_y[omega_count] * face_y >= 0){
               cube_x = face_x;
               cube_y = face_y;
               cube_z = face_z;
           }
        }
     }

     if(sph_z[omega_count] != 0.0){
         for ( face_calc = 0; face_calc <2; face_calc++){
             face_z = face_calc ==0 ? face:-face;
             face_x = (face_z - kk_self*2.0/info_stat.self_sub_vox) * sph_x[omega_count]/ sph_z[omega_count] + ii_self*2.0/info_stat.self_sub_vox;
             face_y = (face_z - kk_self*2.0/info_stat.self_sub_vox) * sph_y[omega_count]/ sph_z[omega_count] + jj_self*2.0/info_stat.self_sub_vox;
             if (face_x <= face && face_x >= -face &&  face_y <= face && face_y >= -face  && face_z <= face && face_z >= -face && sph_z[omega_count] * face_z >=0){
                 cube_x = face_x;
                 cube_y = face_y;
                 cube_z = face_z;
             }
         }
     }
#endif
	 dist_self = sqrt( (ii_self*2.0/info_stat.self_sub_vox - cube_x)*(ii_self*2.0/info_stat.self_sub_vox - cube_x) + (jj_self*2.0/info_stat.self_sub_vox- cube_y)*(jj_self*2.0/info_stat.self_sub_vox- cube_y) + (kk_self*2.0/info_stat.self_sub_vox - cube_z)*(kk_self*2.0/info_stat.self_sub_vox - cube_z)) * info_stat.delX/2.0; //square voxel approx.

     fact_self_vox[omega_count * info_stat.self_sub_vox * info_stat.self_sub_vox * info_stat.self_sub_vox + r_ind_self ] =  ( 1 - cexp( -(info_stat.mu_tot[tiss_num]) * dist_self)) *   sin(theta_self[omega_count]);

}

#if 0

void generate_diag_terms_dev() {

    int ang_res = ANG_RES;
    int omega_count,ang_ind;
    int r_ind_self;
    complex_double *rt_self, *rtp_self;
    int l,m,lp,mp,cnt,cntp;
    flt_doub *theta_self, *phi_self, *sph_x_self, *sph_y_self, *sph_z_self;
    flt_doub cm, del_theta, del_phi;
    complex_double sub_v_sum_self;

    int i;
    cm = C / phan->n;


    diag_terms_host = (complex_double *)malloc(sizeof(complex_double )* MAX_TISS_NUM*MAX_NL*MAX_NL*MAX_NL*MAX_NL);

    theta_self = (flt_doub *)malloc(sizeof(flt_doub) * pow( ang_res,2));
    phi_self = (flt_doub *)malloc(sizeof(flt_doub) * pow(ang_res,2));
    sph_x_self = (flt_doub *)malloc(sizeof(flt_doub) * pow(ang_res,2));
    sph_y_self = (flt_doub *)malloc(sizeof(flt_doub) * pow(ang_res,2));
    sph_z_self = (flt_doub *)malloc(sizeof(flt_doub) * pow(ang_res,2));
    rt_self = (complex_double *)malloc(sizeof(complex_double) * pow(ang_res,2));
    rtp_self = (complex_double *)malloc(sizeof(complex_double) * pow(ang_res,2));

   flt_doub *theta_self_dev, *phi_self_dev,*sph_x_dev,*sph_y_dev,*sph_z_dev;
   complex_double *fact_self_vox_dev, *fact_self_vox_host;

   cudaMalloc(&theta_self_dev, sizeof(flt_doub)*pow( ang_res,2));
   MY_SAFE_CALL(cudaMalloc(&phi_self_dev, sizeof(flt_doub)*pow( ang_res,2)));
   MY_SAFE_CALL(cudaMalloc(&sph_x_dev,sizeof(flt_doub)*pow( ang_res,2)));
   MY_SAFE_CALL(cudaMalloc(&sph_y_dev,sizeof(flt_doub)*pow( ang_res,2)));
   MY_SAFE_CALL(cudaMalloc(&sph_z_dev,sizeof(flt_doub)*pow( ang_res,2)));
   MY_SAFE_CALL(cudaMalloc(&fact_self_vox_dev, sizeof(complex_double ) * pow(geom->self_sub_vox,3) * pow( ang_res,2)));
   fact_self_vox_host = (complex_double *) malloc (sizeof(complex_double ) * pow(geom->self_sub_vox,3) * pow( ang_res,2));

   if(fact_self_vox_host == NULL){
       printf("error in memory allocation \n");
       exit(0);
   }

   dim3 dim_block_1(BLK_ANG_SIZE,BLK_ANG_SIZE,1);
   dim3 dim_grid_1(1,1);

   dim3 dim_block_2(BLK_SELF_SUB_VOX,BLK_SELF_SUB_VOX,BLK_SELF_SUB_VOX);
   dim3 dim_grid_2(geom->self_sub_vox/BLK_SELF_SUB_VOX,geom->self_sub_vox/BLK_SELF_SUB_VOX);

   int theta_count, phi_count;
   for(theta_count = 0; theta_count < ANG_RES; theta_count = theta_count + BLK_ANG_SIZE){
      for( phi_count=0; phi_count < ANG_RES; phi_count = phi_count + BLK_ANG_SIZE){
	      compute_sph_coord<<<dim_grid_1, dim_block_1>>>(theta_self_dev, phi_self_dev, sph_x_dev, sph_y_dev, sph_z_dev, theta_count, phi_count);
		  checkCUDAError("Kernel invocation in compute_sph_coord\n");
      }
  }

   cudaMemcpy(theta_self, theta_self_dev, sizeof(flt_doub)*pow( ang_res,2), cudaMemcpyDeviceToHost);
   cudaMemcpy(phi_self, phi_self_dev, sizeof(flt_doub)*pow( ang_res,2), cudaMemcpyDeviceToHost);

   omega_count = 0;
/*
   for(theta_count = 0; theta_count < ANG_RES; theta_count = theta_count + BLK_ANG_SIZE){
      for( phi_count=0; phi_count < ANG_RES; phi_count = phi_count + BLK_ANG_SIZE){
		  omega_count = theta_count * ANG_RES + phi_count;
//		  printf("%f %f %f  \n", sph_x_self[omega_count], sph_y_self[omega_count],sph_z_self[omega_count], omega_count);
	  }
   }
*/
    del_theta =  M_PI / ANG_RES;
    del_phi =  2*M_PI / ANG_RES;

    int tiss_num;
    int blk_dep;

    omega_count = 0;
    for (tiss_num = 1; tiss_num < phan->no_tiss; tiss_num++){
           for ( omega_count = 0; omega_count  < pow(ang_res,2); omega_count++){
              for(blk_dep=0; blk_dep < geom->self_sub_vox; blk_dep = blk_dep + BLK_SELF_SUB_VOX){
                  compute_diag_selfsub<<<dim_grid_2, dim_block_2>>>(fact_self_vox_dev, sph_x_dev, sph_y_dev,sph_z_dev, theta_self_dev, omega_count, tiss_num,blk_dep);
   	     		  checkCUDAError("Kernel invocation in compute_diag_selfsub\n");
         	  }
     	   }
           cudaMemcpy(fact_self_vox_host, fact_self_vox_dev, sizeof(complex_double) * pow(geom->self_sub_vox,3) * pow( ang_res,2), cudaMemcpyDeviceToHost);
                    cnt = 0;
                    for (l = 0; l <= nL; l++) {
                        for (m = -l; m <= l; m++) {
                            cntp = 0;
                            SpherHarmonicArray(l, m, powf(ang_res,2), theta_self, phi_self, rt_self);
                            for (lp = 0; lp <= nL; lp++) {
                                for (mp = -lp; mp <= lp; mp++) {
                                    sub_v_sum_self = 0.0 + 0.0*I;
                                    SpherHarmonicArray(lp, mp, pow(ang_res,2), theta_self, phi_self, rtp_self);
                                    for ( omega_count = 0; omega_count < ang_res * ang_res; omega_count++){
                                    	for ( r_ind_self = 0; r_ind_self < pow(geom->self_sub_vox,3); r_ind_self++){
                                            sub_v_sum_self = sub_v_sum_self +  ~(rt_self[omega_count]) * rtp_self[omega_count] *  fact_self_vox_host[omega_count * (int)pow(geom->self_sub_vox,3) + r_ind_self];
                                        }
                                    }
                                    diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num] = sub_v_sum_self  * del_theta * del_phi / (cm * pow((double)geom->self_sub_vox,3) * (phan->mu_abs[tiss_num] + phan->mu_sc[tiss_num]) * geom->delX * geom->delY * geom->delZ) ;
                        if(cnt == cntp){
                                    printf("The diagonal term is %e +%e i for tiss = %d, cnt = %d and cntp = %d \n", diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num].real(), diag_terms_host[cnt*MAX_TISS_NUM*(MAX_NL)*(MAX_NL) + cntp * MAX_TISS_NUM + tiss_num].imag(), tiss_num, cnt, cntp);
						}
                                    cntp++;
                                }
                        }
                        cnt++;
                      }
                    }
   }
   cudaFree(sph_x_dev);
   cudaFree(sph_y_dev);
   cudaFree(sph_z_dev);
   cudaFree(theta_self_dev);
   cudaFree(phi_self_dev);
   exit(0);
}
#endif

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg,
                                  cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }
}


void * prop_abs(void *arg){
#if 1
    SHORT cnt;
    SHORT block_dep;
    COR subvox_src, subvox_dest;
    SIGNED_SHORT i,j,k,ip,jp,kp;
	flt_doub dx_sub, dy_sub, dz_sub;

    dim3 dim_block(BLK_SIZE, BLK_SIZE, BLK_SIZE_Z);
    dim3 dim_grid(geom->nX/dim_block.x,geom->nY/dim_block.y);

//    printf("% d and %d are no of blocks per grid \n", dim_grid.y, dim_grid.x);
//    printf("% d %d and %d  are no of threads per block \n", dim_block.x, dim_block.y, dim_block.z);


	const THREAD_PARAMETERS parameters = *((THREAD_PARAMETERS *) arg);
	
	int size_layer = ( geom->nX + 2*geom->bounX ) * ( geom->nY + 2*geom->bounY ) * ( nL+1) * (nL+1) ;
	int size = ( geom->nX + 2*geom->bounX ) * ( geom->nY + 2*geom->bounY ) * (geom->nZ + 2*geom->bounZ) * ( nL+1) * (nL+1);
    
	complex_double *src_dev, *out_dev;

	if(cudaSetDevice(parameters.device_index) != cudaSuccess){
		printf("Error in setting up device %d \n", parameters.device_index);
		exit(0);
	}

	
    MY_SAFE_CALL(cudaMalloc(&src_dev, sizeof(complex_double)*size)); 
    MY_SAFE_CALL(cudaMalloc(&out_dev, sizeof(complex_double)*size_layer*(parameters.num_layers))); 
    MY_SAFE_CALL(cudaMemset(out_dev, 0, sizeof(complex_double)*size_layer*(parameters.num_layers)));

	MY_SAFE_CALL(cudaMemcpy(src_dev, parameters.src_host,sizeof(complex_double)*size, cudaMemcpyHostToDevice));

    Info_Dyn info_dyn_dev;
    MY_SAFE_CALL(cudaMalloc(&(info_dyn_dev.tiss_type), sizeof(byte)*geom->no_vox));
    MY_SAFE_CALL(cudaMemcpy(info_dyn_dev.tiss_type, phan->tiss_type,sizeof(byte)*geom->no_vox, cudaMemcpyHostToDevice));
	
	MY_SAFE_CALL(cudaMemcpyToSymbol(info_stat,info_stat_host,sizeof(Info_Stat) ));
    MY_SAFE_CALL(cudaMemcpyToSymbol(nL_dev,&nL,sizeof(SHORT) ));
    MY_SAFE_CALL(cudaMemcpyToSymbol(diag_terms_dev,diag_terms_host, sizeof(complex_double)*MAX_TISS_NUM*MAX_NL*MAX_NL*MAX_NL*MAX_NL));

	complex_double *sph_harm_dev;
    MY_SAFE_CALL(cudaMalloc(&sph_harm_dev, sizeof(complex_double)*(nL+1)*(nL+1) * THETA_ANG_RES * PHI_ANG_RES));
    MY_SAFE_CALL(cudaMemcpy(sph_harm_dev,sph_harm, sizeof(complex_double)*(nL+1)*(nL+1) * THETA_ANG_RES * PHI_ANG_RES,cudaMemcpyHostToDevice));


    //for(cnt=0; cnt < (nL+1)*(nL+1); cnt++){

//        printf("Invoking compute_diagonal_abs with cnt = %d  \n", cnt);
//        for(block_dep = 0; block_dep < parameters.num_layers; block_dep = block_dep + BLK_SIZE_Z){
//				write_dist_dev<<<dim_grid, dim_block>>>(src_dev,1+0.0*I,cnt,block_dep);
//		}
#if 1
        for(block_dep = 0; block_dep < parameters.num_layers; block_dep = block_dep + BLK_SIZE_Z){
                compute_diagonal_abs<<<dim_grid, dim_block>>>(src_dev,out_dev,info_dyn_dev,cnt,block_dep, parameters.layer_start, parameters.flag);
                checkCUDAError("Kernel invocation in compute_diagonal_abs\n");
        }
#endif


    /* The prop_thresh condition. Again run thread for all the voxels */

#if 1
//        printf("Invoking compute_propabs with cnt = %d  \n", cnt);
        for(block_dep = 0; block_dep < parameters.num_layers; block_dep = block_dep + BLK_SIZE_Z){
                    compute_propabs<<<dim_grid, dim_block>>>(src_dev, out_dev, info_dyn_dev, sph_harm_dev, cnt, block_dep,parameters.layer_start,parameters.flag);
				    //printf("%d operation complete \n", block_dep/parameters.num_layers);
                    checkCUDAError("Kernel invocation in compute_propabs\n");
        }

#endif


        for(block_dep = 0; block_dep < parameters.num_layers; block_dep = block_dep + BLK_SIZE_Z){
                scale_dist_dev<<<dim_grid, dim_block>>>(out_dev, geom->delX * geom->delY * geom->delZ,cnt,block_dep);
                checkCUDAError("Kernel invocation in scale_dist_dev\n");
                cudaThreadSynchronize();
        }
        cudaThreadSynchronize();
	//}
    
	MY_SAFE_CALL(cudaMemcpy(parameters.out_host, out_dev, sizeof(complex_double)*size_layer*(parameters.num_layers), cudaMemcpyDeviceToHost)); 
    MY_SAFE_CALL(cudaMemset(out_dev, 0, sizeof(complex_double)*size_layer*(parameters.num_layers)));
    MY_SAFE_CALL(cudaMemset(src_dev, 0, sizeof(complex_double)*size));
    MY_SAFE_CALL(cudaMemset(sph_harm_dev, 0, sizeof(complex_double)*(nL+1)*(nL+1) * THETA_ANG_RES * PHI_ANG_RES));
	cudaFree(src_dev);
    cudaFree(out_dev);
	cudaFree(sph_harm_dev);
	cudaThreadExit();
    pthread_exit(NULL);

#endif 
}

