#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <complex.h>

#include "iniparser.h"
#include "RTE.h"

#pragma warning(disable:981)

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define THETA_ANG_RES 100
#define PHI_ANG_RES 200


/*
PropScat

Propaage the src using the scatter transform.  This the the scatter propagator.
In spherical harmonics, this is pretty easy.  Input is src.  Output
is W.  This is for a uniform phantom.
*/
void PropScat(Geometry *geom,Phantom *phan,int nL,Dist *src,Dist *W) {
	int nCoef;
	int i,j,k,l,m;
	int cnt;
	
	nCoef = (nL+1)*(nL+1);
	cnt = 0;
	for (l=0;l<=nL;l++) {
		for (m=-l;m<=l;m++) {
			for (i=0;i<geom->nX;i++) {
				for (j=0;j<geom->nY;j++) {
					for (k=0;k<geom->nZ;k++) {
						W->dist[cnt][i+ j*geom->nX+k*geom->nX*geom->nY] = (phan->c/phan->n) * pow(phan->g,(double)l) * phan->mus * 
							src->dist[cnt][i+ j*geom->nX+k*geom->nX*geom->nY];
					}
				}
			}
			cnt++;
		}
	}
	
}


/*
PropScat

Propaage the src using the scatter transform.  This the the scatter propagator.
In spherical harmonics, this is pretty easy.  Input is src.  Output
is W.  This is for a uniform phantom.
*/
void PropScatmu1(Geometry *geom,Phantom *phan,int nL,Dist *src) {
	int nCoef;
	int i,j,k,l,m;
	int cnt;
	
	nCoef = (nL+1)*(nL+1);
	cnt = 0;
	for (l=0;l<=nL;l++) {
		for (m=-l;m<=l;m++) {
			for (i=0;i<geom->nX;i++) {
				for (j=0;j<geom->nY;j++) {
					for (k=0;k<geom->nZ;k++) {
						src->dist[cnt][i+ j*geom->nX+k*geom->nX*geom->nY] = (phan->c/phan->n) * pow(phan->g,(double)l) * 1.0 * 
							src->dist[cnt][i+ j*geom->nX+k*geom->nX*geom->nY];
					}
				}
			}
			cnt++;
		}
	}
	
}



/*
PropAbs

Propaage the src using the x-ray transform.  This the the absorption propagator.
In spherical harmonics, this takes the bulk of the time.  Input is src.  Output
is W.  This is for a uniform phantom.

Note: We don't need the phantom here because it was used to generate all the values
in Terms *terms.
*/
void PropAbs(Geometry *geom,Terms *terms,int nL,Dist *src,Dist *W) {
	int i,j,k;
	int ip,jp,kp; // p stands for prime
	int l,m;
	int n;
	int cnt,cntp;
	int pcntage;
	complex double tmp;
	complex double b;
	
	// for speed
	int fact1,fact2;
	int FACT1,FACT2;
	NonZeros *zer;
	
	zer = GenerateNonZeros(geom,src);
	
	fact1 = geom->nX;
	fact2 = geom->nX * geom->nY;
	FACT1 = terms->NX;
	FACT2 = terms->NX*terms->NY;
	// Start the loops.

	pcntage = 0;

#pragma omp parallel private(tmp,b,cnt,cntp,n,j,k,l,m,ip,jp,kp)	
{
#pragma omp for	
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			for (k=0;k<geom->nZ;k++) {
				cnt = 0;
				for (l=0;l<=nL;l++) {
					for (m=-l;m<=l;m++) {
						tmp = 0 + 0 * I;
						// Only analyze non-zero terms--to speed up
						for (n=0;n<zer->num;n++) {
							ip = zer->is[n];
							jp = zer->js[n];
							kp = zer->ks[n];
							cntp = zer->cnts[n];
							if ((ip==i) && (jp==j) && (kp==k) && (cntp != cnt)) {
								b = 0 + 0*I;
							} else {
								b = terms->fact[geom->nX-i+ip-1 + (geom->nY-j+jp-1)*FACT1 + (geom->nZ-k+kp-1)*FACT2]
									* conj(terms->Ylm[cnt][geom->nX-i+ip-1 + (geom->nY-j+jp-1)*FACT1 + (geom->nZ-k+kp-1)*FACT2]);
							}
							tmp = tmp + b * 
								terms->Ylm[cntp][geom->nX-i+ip-1 + (geom->nY-j+jp-1)*FACT1 + (geom->nZ-k+kp-1)*FACT2] *
								src->dist[cntp][ip+jp*fact1+kp*fact2];
						}
						W->dist[cnt][i+ j*fact1+k*fact2] = tmp * geom->delX*geom->delY*geom->delZ; //* // NEW!
//										geom->delX * geom->delY * geom->delZ;
						cnt++;
					}
				}
			}
		}
		pcntage++;
		//printf("   %2.0f percent complete       \r",(double)pcntage/(double)(geom->nX)*100.0);
		fflush(stdout);
	}
}
//	printf("   100 percent complete    \n");
	FreeNonZeros(zer);
}



/*
PropAbs2

Propaage the src using the x-ray transform.  This the the absorption propagator.
In spherical harmonics, this takes the bulk of the time.  Input is src.  Output
is W.  This is for a uniform phantom.

Note: We don't need the phantom here because it was used to generate all the values
in Terms *terms.
*/
void PropAbs2(Geometry *geom,Terms2 *terms,int nL,Dist *src,Dist *W) {
	int i,j,k;
	int ip,jp,kp; // p stands for prime
	int l,m;
	int lp,mp;
	int n;
	int cnt,cntp;
	int pcntage;
	complex double tmp;
	complex double b;
	int nearX,nearY,nearZ;
	int ipStart,ipEnd;
	int jpStart,jpEnd;
	int kpStart,kpEnd;
	
	// for speed
	int fact1,fact2;
	int FACT1,FACT2;
	NonZeros *zer;
	
//	zer = GenerateNonZeros(geom,src);
	
	fact1 = geom->nX;
	fact2 = geom->nX * geom->nY;
	FACT1 = terms->NX;
	FACT2 = terms->NX*terms->NY;
	// Start the loops.

	nearX = ceil(geom->propThresh / (double)geom->delX);
	nearY = ceil(geom->propThresh / (double)geom->delY);
	nearZ = ceil(geom->propThresh / (double)geom->delZ);

	pcntage = 0;

#pragma omp parallel private(tmp,b,cnt,cntp,n,j,k,l,m,lp,mp,ip,jp,kp,ipStart,ipEnd,jpStart,jpEnd,kpStart,kpEnd)
{
#pragma omp for	
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			for (k=0;k<geom->nZ;k++) {
				cnt = 0;
				ipStart = max(0,i-nearX);
				ipEnd = min(geom->nX-1,i+nearX);
				jpStart = max(0,j-nearY);
				jpEnd = min(geom->nY-1,j+nearY);
				kpStart = max(0,k-nearZ);
				kpEnd = min(geom->nZ-1,k+nearZ);
				for (l=0;l<=nL;l++) {
					for (m=-l;m<=l;m++) {
						tmp = 0 + 0 * I;
						for (ip = ipStart;ip <=ipEnd;ip++) {
							for (jp = jpStart;jp <=jpEnd;jp++) {
								for (kp=kpStart;kp<=kpEnd;kp++) {
									cntp = 0;
									for (lp = 0;lp <=nL;lp++) {
										for (mp=-lp;mp<=lp;mp++) {
											tmp = tmp + terms->fact[cnt][cntp][geom->nX-i+ip-1 + (geom->nY-j+jp-1)*FACT1 + (geom->nZ-k+kp-1)*FACT2] *
												src->dist[cntp][ip+jp*fact1+kp*fact2];
											cntp++;
										}
									}
								}
							}
						}
						// Only analyze non-zero terms--to speed up
/*						for (n=0;n<zer->num;n++) {
							ip = zer->is[n];
							jp = zer->js[n];
							kp = zer->ks[n];
							cntp = zer->cnts[n];
							tmp = tmp + terms->fact[cnt][cntp][geom->nX-i+ip-1 + (geom->nY-j+jp-1)*FACT1 + (geom->nZ-k+kp-1)*FACT2] *
								src->dist[cntp][ip+jp*fact1+kp*fact2];
						} */
						W->dist[cnt][i+ j*fact1+k*fact2] = tmp; // * geom->delX*geom->delY*geom->delZ;
						cnt++;
					}
				}
			}
		}
		pcntage++;
//		printf("   %2.0f percent complete       \r",(double)pcntage/(double)(geom->nX)*100.0);
		fflush(stdout);
	}
}
//	printf("   100 percent complete    \n");
//	FreeNonZeros(zer);
}




/*
GenerateTerms

Generate the terms needed for this the Propagator term.  These terms are a function
of r - r' and are stored as such.
*/
Terms *GenerateTerms(Geometry *geom,Phantom *phan,int nL) {
	Terms *ret;
	double *theta,*phi;
	int NX,NY,NZ;
	int l,m;
	int nCoef;
	int i,j,k;
	int cnt;
	double *delx,*dely,*delz;
	double xmin,xmax,ymin,ymax,zmax;
	double r;
	double c;
	complex double mutotal;
	// Variables for sub-voxelizations
	int ii,jj,kk;
	int iip,jjp,kkp;
	double subVoxSpacing;
	double dx,dy,dz;
	double subVSum;
	double subVR;
	
	
	
	nCoef = (nL+1)*(nL+1);
	NX = geom->nX + geom->nX - 1;
	NY = geom->nY + geom->nY - 1;
	NZ = geom->nZ + geom->nZ - 1;
	
	c = phan->c / phan->n;
	mutotal = (phan->mua + phan->mus);
	
	ret = malloc(sizeof(Terms)*1);
	ret->nL = nL;
	ret->NX = NX;
	ret->NY = NY;
	ret->NZ = NZ;
	ret->Ylm = malloc(sizeof(complex double*)*nCoef);
	for (i=0;i<nCoef;i++) {
		ret->Ylm[i] = malloc(sizeof(complex double)*NX*NY*NZ);
	}
	ret->fact = malloc(sizeof(complex double)*NX*NY*NZ);
	
	
	delx = malloc(sizeof(double)*NX);
	dely = malloc(sizeof(double)*NY);
	delz = malloc(sizeof(double)*NZ);
	theta = malloc(sizeof(double)*NX*NY*NZ);
	phi = malloc(sizeof(double)*NX*NY*NZ);
	xmin = geom->xs[0];
	xmax = geom->xs[geom->nX - 1];
	ymin = geom->ys[0];
	ymax = geom->ys[geom->nY - 1];
	zmax = geom->zs[geom->nZ - 1];
	
	for (i=0;i<NX;i++) {
		delx[i] = (2.0*xmin - 2.0*xmax)/((double)NX - 1.0) * (double)i + 2.0*xmax;
	}
	for (j=0;j<NY;j++) {
		dely[j] = (2.0*ymin - 2.0*ymax)/((double)NY - 1.0) * (double)j + 2.0*ymax;
	}
	for (k=0;k<NZ;k++) {
		delz[k] = (-2.0*zmax)/((double)NZ - 1.0) * (double)k + zmax;
	}
		
	for (i=0;i<NX;i++) {
		for (j=0;j<NY;j++) {
			for (k=0;k<NZ;k++) {
				theta[i + NX*j + NX*NY*k] = atan2(sqrt(delx[i]*delx[i] + dely[j]*dely[j]),delz[k]);
				phi[i+NX*j + NX*NY*k] = atan2(dely[j],delx[i]);
				
				r = sqrt(delx[i]*delx[i] + dely[j]*dely[j] + delz[k]*delz[k]);
				if (r == 0.0) {
					ret->fact[i+NX*j+NX*NY*k] = 
						//(1.0/(geom->delX*geom->delY*geom->delZ*mutotal*c)) * (1.0-cexp(-mutotal*geom->delX/2.0));    // EC/MK 5/19
						(.07398629)/(c*geom->delX*geom->delY*geom->delZ);
				} else {
					if (r < geom->subThresh) {
						subVSum = 0.0;
						for (ii = 0;ii < geom->subVox;ii++) {
							for (jj=0;jj<geom->subVox;jj++) {
								for (kk=0;kk<geom->subVox;kk++) {
									for (iip=0;iip<geom->subVox;iip++) {
										dx = delx[i] - geom->delX + (geom->delX/(double)geom->subVox)* (double)(1+ii+iip);
										for (jjp=0;jjp<geom->subVox;jjp++) {
											dy = dely[j] - geom->delY + (geom->delY/(double)geom->subVox)* (double)(1+jj+jjp);
											for (kkp=0;kkp<geom->subVox;kkp++) {
												dz = delz[k] - geom->delZ + (geom->delZ/(double)geom->subVox)* (double)(1+kk+kkp);
												subVR = sqrt(dx*dx+dy*dy+dz*dz);
												subVSum += (1.0/(subVR*subVR))*exp(-mutotal*subVR);												
											}
										}
									}
								}
							}
						}
						subVSum = subVSum/(pow((double)(geom->subVox),6.0)*c);
						ret->fact[i+NX*j+NX*NY*k] = subVSum;
						//printf("%f %e %e %e\n",r,subVSum,creal(ret->fact[i+NX*j+NX*NY*k]),subVSum/creal(ret->fact[i+NX*j+NX*NY*k]));
					} else {
						ret->fact[i+NX*j+NX*NY*k] = (1.0/c) * (1.0/(r*r)) * cexp(-mutotal*r);
					}
				}
			}
		}
	}
	cnt = 0;
	for (l=0;l<=nL;l++) {
		for (m=-l;m<=l;m++) {
			SpherHarmonicArray(l,m,NX*NY*NZ,theta,phi,ret->Ylm[cnt]);
			cnt++;
		}
	}

	for (i=0;i<NX;i++) {
		for (j=0;j<NY;j++) {
			for (k=0;k<NZ;k++) {
				if ((delx[i]==0) && (dely[j]==0) && (delz[k]==0)) {
					cnt = 0;
					for (l=0;l<=nL;l++) {
						for (m=-l;m<=l;m++) {
							ret->Ylm[cnt][i + j*NX + k*NX*NY] = 1.0 + 0*I;
							cnt++;
						}
					}
				}
			}
		}
	}
		
	free(delx);
	free(dely);
	free(delz);
	free(theta);
	free(phi);
	
	return(ret);
}




/*
FreeTerms

Free up the memory used by the r-r' factors
*/
void FreeTerms(Terms *fac) {
	int i;
	int nCoef;
	
	nCoef = (fac->nL+1)*(fac->nL+1);
	
	for (i=0;i<nCoef;i++) {
		free(fac->Ylm[i]);
	}
	free(fac->Ylm);
	free(fac->fact);
	free(fac);
}


/*
GenerateTerms2

Generate the terms needed for this the Propagator term.  These terms are a function
of r - r' and are stored as such.
*/
Terms2 *GenerateTerms2(Geometry *geom,Phantom *phan,int nL) {
	Terms2 *ret;
	int NX,NY,NZ;
	int l,m,lp,mp;
	int nCoef;
	int i,j,k;
	int cnt,cntp;
	double *delx,*dely,*delz;
	double xmin,xmax,ymin,ymax,zmax,zmin;
	double r;
	double c;
	complex double mutotal;
	// Variables for sub-voxelizations
	int ii,jj,kk;
	int iip,jjp,kkp;
	double subVoxSpacing;
	double dx,dy,dz;
	complex double subVSum;
	double subVR;
	complex double *subVVal;
	double theta,phi;
	complex double *rt,*rtp;
	double thetan,phin;
	int subCounter;
	int totalNumber;

  //Declaring variables for self-sub-vox
    int totalNumber_self;
    int ang_res = 100;
    int theta_count, phi_count, omega_count;
    double *theta_self, *phi_self, *flag_self;
    float cube_x, cube_y, cube_z;
    complex double *subVVal_self, *rt_self, *rtp_self;
    int face_calc;
    float *sph_x, *sph_y, *sph_z, face_x, face_y, face_z;
    double dist;
    complex double **fact_self_vox;
    double del_phi, del_theta;
    int r_ind;
    int ang_ind, total_r;
    complex double sub_v_sum;
    float ii_self, jj_self, kk_self;
    int face = 1; // Defines the face of the cube to be considered later.
    complex double tmp_self;
    //Declaration ends here

    
	nCoef = (nL+1)*(nL+1);
	NX = geom->nX + geom->nX - 1;
	NY = geom->nY + geom->nY - 1;
	NZ = geom->nZ + geom->nZ - 1;
	
	c = phan->c / phan->n;
	mutotal = (phan->mua + phan->mus);
	
	ret = malloc(sizeof(Terms2)*1);
	ret->nL = nL;
	ret->NX = NX;
	ret->NY = NY;
	ret->NZ = NZ;
	ret->fact = malloc(sizeof(complex double*)*nCoef);
	for (i=0;i<nCoef;i++) {
		ret->fact[i] = malloc(sizeof(complex double*)*nCoef);
		for (j=0;j<nCoef;j++) {
			ret->fact[i][j] = malloc(sizeof(complex double)*NX*NY*NZ);
			if (ret->fact[i][j] == 0L) {
				printf("Memory allocation error\n");
			}
		}
	}
	
	
	delx = malloc(sizeof(double)*NX);
	dely = malloc(sizeof(double)*NY);
	delz = malloc(sizeof(double)*NZ);
	xmin = geom->xs[0];
	xmax = geom->xs[geom->nX - 1];
	ymin = geom->ys[0];
	ymax = geom->ys[geom->nY - 1];
	zmin = geom->zs[0];
	zmax = geom->zs[geom->nZ - 1];
	
	totalNumber = geom->subVox*geom->subVox* geom->subVox*geom->subVox* geom->subVox*geom->subVox;
	subVVal = malloc(sizeof(complex double) * totalNumber);
	rt = malloc(sizeof(complex double)*totalNumber);
	rtp = malloc(sizeof(complex double)*totalNumber);
	
	for (i=0;i<NX;i++) {
		delx[i] = (2.0*xmin - 2.0*xmax)/((double)NX - 1.0) * (double)i + xmax - xmin;
	}
	for (j=0;j<NY;j++) {
		dely[j] = (2.0*ymin - 2.0*ymax)/((double)NY - 1.0) * (double)j + ymax - ymin;
	}
	for (k=0;k<NZ;k++) {
		delz[k] = zmax-zmin - (double)k * geom->delZ;
		// Used to be zmin - zmax + (double)k*geom->delZ
		//  AJ :  10/11/2010
	}

	int ip, jp, kp, nVox;
    float delTheta, delPhi, delR;
    delTheta =  M_PI/(geom->subTheta);
    delPhi = 2*M_PI/(geom->subPhi);
    delR = geom->delX*sqrt(3)/(geom->subR);
    float Rmin;
    Rmin = geom->delX < geom->delY ? geom->delX : geom->delY;
    Rmin = Rmin < geom->delZ ? Rmin : geom->delZ;

    //Memory allocation for self sub voxelization
    fact_self_vox = malloc ( sizeof(complex double *) * pow(2*geom->self_sub_vox+1,3));
    for ( i = 0; i < pow(2*geom->self_sub_vox+1,3); i++){
        fact_self_vox[i] = malloc(sizeof(complex double) * pow(ang_res,2));
        if (fact_self_vox[i] == 0L) {
            printf("Memory allocation error \n");
        }
    }


    sph_x = malloc(sizeof(double) * pow( ang_res,2));
    sph_y = malloc(sizeof(double) * pow( ang_res,2));
    sph_z = malloc(sizeof(double) * pow( ang_res,2));

    theta_self = malloc(sizeof(double) * pow( ang_res,2));
    phi_self = malloc(sizeof(double) * pow(ang_res,2));
    rt_self = malloc(sizeof(complex double) * pow(ang_res,2));
    rtp_self = malloc(sizeof(complex double) * pow(ang_res,2));


    // Memory allocation ends here

    //for(r = Rmin/(2.0 * geom->subVox); r < geom->propThresh; r += delR){
    for(r = 0; r < geom->propThresh; r += delR){
		for(theta = 0; theta < M_PI; theta += delTheta){
			for(phi = 0; phi < 2*M_PI; phi += delPhi){
    			  for (ii=0;ii < geom->subVox;ii++) {
				        for (jj=0;jj < geom->subVox;jj++) {
	    					  for (kk=0;kk < geom->subVox;kk++) {
             						ip = (int) round((ii + 0.5 - geom->subVox/2.0)/geom->subVox + r*sin(theta)*cos(phi)/geom->delX) + (NX-1)/2; // Adding the offset of (NX-1)/2.0 since array indices should not be negative. 
									jp = (int) round((jj + 0.5 - geom->subVox/2.0)/geom->subVox + r*sin(theta)*sin(phi)/geom->delY) + (NY-1)/2;
									kp = (int) round((kk + 0.5 - geom->subVox/2.0)/geom->subVox + r*cos(theta)/geom->delZ) + (NZ-1)/2;
							    	if(ip >=0 && ip < NX && jp >= 0 && jp < NY && kp >= 0 && kp < NZ){ 
    												cnt = 0;
												    for (l=0;l<=nL;l++) {
												       for (m=-l;m<=l;m++) {
												           cntp = 0;
												           for (lp=0;lp<=nL;lp++) {
												             for (mp=-lp;mp<=lp;mp++) {
																ret->fact[cnt][cntp][ip + NX * jp + NX * NY * kp] += cexp(-mutotal*r)* conj(SpherHarmonic(theta,phi,l,m)) * SpherHarmonic(theta,phi,lp,mp) * sin(theta); 
															}
														}
													}
									  			}
									}
						    }
     				     }
					}
			}
		}
	}

    for(cnt = 0; cnt < (nL+1)*(nL+1); cnt++){
    	for(cntp = 0; cntp < (nL+1)*(nL+1); cntp++){
			for(ip = 0; ip < NX*NY*NZ; ip++){
						ret->fact[cnt][cntp][ip ] *= delR * delTheta * delPhi/(c * pow(geom->subVox,3));
	}}}


    for(cnt = 0; cnt < (nL+1)*(nL+1); cnt++){
        for(cntp = 0; cntp < (nL+1)*(nL+1); cntp++){
            for(ip = 0; ip < NX; ip++){
                for(jp=0; jp< NY; jp++){
                    for(kp=0; kp<NZ; kp++){
                        if(ip <= (NX+1)/2 && ip >= (NX-3)/2 &&  jp >= (NY-3)/2 && jp <= (NY+1)/2 && kp >= (NZ-3)/2 && kp <= (NZ+1)/2 ){
                                printf(" (%d %d %d) : %e + %e i  \n",ip -(NX-1)/2,jp - (NY-1)/2,kp - (NZ-1)/2, creal(ret->fact[cnt][cntp][ip + NX * jp + NX * NY * kp]), cimag(ret->fact[cnt][cntp][ip + NX * jp + NX * NY * kp]));
    }}}}}}
   
	omega_count = 0;
    for ( theta_count = 0; theta_count < ang_res; theta_count++){
        for ( phi_count = 0; phi_count  < ang_res; phi_count++){
            theta_self[omega_count] = theta_count * M_PI / ang_res;
            phi_self[omega_count] = phi_count * 2.0*M_PI / ang_res;
            sph_x[omega_count] = cos(phi_self[omega_count]) * sin(theta_self[omega_count]);
            sph_y[omega_count] = sin(phi_self[omega_count]) * sin(theta_self[omega_count]) ;
            sph_z[omega_count] = cos(theta_self[omega_count]) ;
            omega_count++;
        }
    }


    del_theta = theta_self[ang_res] - theta_self[0];
    del_phi = phi_self[1] - phi_self[0];



                    sub_v_sum = 0.0 + 0 * I;
                    subCounter = 0;
                    r_ind = 0;
                    //printf("%f %f is ii range \n", -geom->self_sub_vox/2.0 +0.5, -geom->self_sub_vox/2.0 +0.5);
					#pragma omp parallel private(jp,kp,ii_self,jj_self, kk_self, face_x, face_y, face_z, cube_x, cube_y, cube_z,face_calc, r_ind, omega_count) 
					{
					#pragma omp for
					for(ip = 0; ip < geom->self_sub_vox; ip++){
						for(jp = 0; jp < geom->self_sub_vox; jp++){ 
							for(kp = 0;kp < geom->self_sub_vox; kp++){
/*                    for (ii_self = -geom->self_sub_vox/2.0 +0.5; ii_self <= geom->self_sub_vox/2.0 -0.5; ii_self++) {
                        for (jj_self = -geom->self_sub_vox/2.0 +0.5; jj_self <= geom->self_sub_vox/2.0-0.5; jj_self++) {
                            for (kk_self = -geom->self_sub_vox/2.0 +0.5; kk_self <= geom->self_sub_vox/2.0 -0.5; kk_self++) { // First subvoxelise the cube. Assume here that cube is of dimension [-1 1]. */
                                ii_self = ip - geom->self_sub_vox/2.0 +0.5;
								jj_self = jp - geom->self_sub_vox/2.0 +0.5;
								kk_self = kp - geom->self_sub_vox/2.0 +0.5;
								r_ind = ip* (geom->self_sub_vox*geom->self_sub_vox) + jp* geom->self_sub_vox + kp;
								//r_ind = (int) (ii_self - 0.5 + geom->self_sub_vox/2.0)* (geom->self_sub_vox*geom->self_sub_vox) + (jj_self - 0.5 + geom->self_sub_vox/2.0) * geom->self_sub_vox + kk_self - 0.5 + geom->self_sub_vox/2.0;
                            //  printf("Enters the self sub vox loop %f %f %f \n", ii_self, jj_self, kk_self);
                                for ( omega_count = 0; omega_count  < pow(ang_res,2); omega_count++){
                                // Now find the intersection with the different faces of lines from this point. We have to determine where any line at an angle (theta, phi) from the point in inside the cube intersects the face, given that we already know that one of the intersections coordinates.
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

                                dist = sqrt(pow(ii_self*2.0/geom->self_sub_vox - cube_x,2.0) + pow(jj_self*2.0/geom->self_sub_vox - cube_y,2.0) + pow(kk_self*2.0/geom->self_sub_vox - cube_z,2.0)) * geom->delX/2.0;
                                fact_self_vox[r_ind][omega_count] = ( 1 - cexp( -mutotal* dist)) *  sin(theta_self[omega_count]);
                            }
                             //r_ind++;
                            }
                        }
                    }
					}
                    cnt = 0;
                    for (l = 0; l <= nL; l++) {
                        for (m = -l; m <= l; m++) {
                            cntp = 0;
                            SpherHarmonicArray(l, m, pow(ang_res,2), theta_self, phi_self, rt_self);
                            for (lp = 0; lp <= nL; lp++) {
                                for (mp = -lp; mp <= lp; mp++) {
                                    sub_v_sum = 0 + 0*I;
                                    SpherHarmonicArray(lp, mp, pow(ang_res,2), theta_self, phi_self, rtp_self);
                                    for ( r_ind = 0; r_ind < pow(geom->self_sub_vox,3); r_ind++){
                                        for ( ang_ind = 0; ang_ind < ang_res * ang_res; ang_ind++){
                                            sub_v_sum +=  conj(rt_self[ang_ind]) * rtp_self[ang_ind] * fact_self_vox[r_ind][ang_ind];
                                        }
                                    }
                                    //  printf("subVsum in the %d loop = %e \n'", r_ind, creal(sub_v_sum));
                                ret->fact[cnt][cntp][(NX-1)/2 + NX * (NY-1)/2 + NX * NY * (NZ-1)/2] = sub_v_sum *del_theta * del_phi / (c * pow((double)geom->self_sub_vox,3) * mutotal) ;
 #if 1
                                 if(cnt == cntp)
                                 printf("subVSum value for the self voxel is %e + %e i for (l,m) = (%d, %d) and (lp,mp) = (%d, %d) \n", creal(ret->fact[cnt][cntp][(NX-1)/2 + NX * (NY-1)/2 + NX * NY * (NZ-1)/2]), cimag(ret->fact[cnt][cntp][i + NX * j + NX * NY * k]), l,m,lp,mp);
 #endif
                                cntp++;
                                }
                        }
                        cnt++;
                      }
                    }
// // Self-sub-voxelization ends
// #if 1
// // This is how it happened previously.
//                     cnt = 0;
//                     for (l=0;l<=nL;l++) {
//                        for (m=-l;m<=l;m++) {
//                             cntp = 0;
//                             for (lp=0;lp<=nL;lp++) {
//                                 for (mp=-lp;mp<=lp;mp++) {
//                                     if (cnt == cntp) {
//                                         tmp_self = /*geom->correctionFact * */ 1*
//                                             (1.0/(geom->delX*geom->delY*geom->delZ*c)) * (1.0-cexp(-mutotal*geom->delX/2.0))/ mutotal;
//                                         printf("Approximate value for the self voxel is %e + %e i\n",creal(tmp_self), cimag(tmp_self));
//                                     cntp++;
//                                 }
//                             cnt++;
//                         }
//                     }
//                     }
//                   }
// #endif
// Diagonal term done

#if 0	
				} else {
					if (r < geom->subThresh) {
						subVSum = 0.0 + 0* I;
						subCounter = 0;
						for (ii = 0;ii < geom->subVox;ii++) {
							for (jj=0;jj<geom->subVox;jj++) {
								for (kk=0;kk<geom->subVox;kk++) {
									for (iip=0;iip<geom->subVox;iip++) {
										dx = delx[i] - geom->delX + (geom->delX/(double)geom->subVox)* (double)(1+ii+iip);
										for (jjp=0;jjp<geom->subVox;jjp++) {
											dy = dely[j] - geom->delY + (geom->delY/(double)geom->subVox)* (double)(1+jj+jjp);
											for (kkp=0;kkp<geom->subVox;kkp++) {
												dz = delz[k] - geom->delZ + (geom->delZ/(double)geom->subVox)* (double)(1+kk+kkp);
												phi[subCounter] = atan2(dy,dx);
												theta[subCounter] = atan2(sqrt(dx*dx + dy*dy),dz);
												subVR = sqrt(dx*dx+dy*dy+dz*dz);
												subVVal[subCounter] = (1.0/(subVR*subVR)) * cexp(-mutotal*subVR);
												subCounter++;
												//subVSum += (1.0/(subVR*subVR))*cexp(-mutotal*subVR)*conj(SpherHarmonic(theta,phi,l,m)) * SpherHarmonic(theta,phi,lp,mp);
											}
										}
									}
								}
							}
						}
						cnt = 0;
						for (l=0;l<=nL;l++) {
							for (m=-l;m<=l;m++) {
								cntp = 0;
								SpherHarmonicArray(l,m,totalNumber,theta,phi,rt);
								for (lp=0;lp<=nL;lp++) {
									for (mp=-lp;mp<=lp;mp++) {
										SpherHarmonicArray(lp,mp,totalNumber,theta,phi,rtp);
										subVSum = 0.0 + 0.0*I;	
										for (subCounter=0;subCounter<totalNumber;subCounter++) {
											subVSum += subVVal[subCounter] * conj(rt[subCounter]) * rtp[subCounter];
										}
										subVSum = subVSum/(pow((double)(geom->subVox),6.0)*c);
										ret->fact[cnt][cntp][i+NX*j+NX*NY*k] = subVSum;
										
										cntp++;
									}
								}
								cnt++;
							}
						}
					} else {
						thetan = atan2(sqrt(delx[i]*delx[i] + dely[j]*dely[j]),delz[k]);
						phin = atan2(dely[j],delx[i]);
						cnt = 0;
						for (l=0;l<=nL;l++) {
							for (m=-l;m<=l;m++) {
								cntp = 0;
								for (lp=0;lp<=nL;lp++) {
									for (mp=-lp;mp<=lp;mp++) {
										ret->fact[cnt][cntp][i + NX*j + NX*NY*k] = (1.0)/(c*r*r)*cexp(-mutotal*r)*conj(SpherHarmonic(thetan,phin,l,m)) * SpherHarmonic(thetan,phin,lp,mp);
										cntp++;
									}
								}
								cnt++;
							}
						}
					}
				}
			}
		}
	}
#endif
	
	
	free(delx);
	free(dely);
	free(delz);
	free(subVVal);
	free(rt);
	free(rtp);
	
	return(ret);
}




/*
FreeTerms2

Free up the memory used by the r-r' factors
*/
void FreeTerms2(Terms2 *fac) {
	int i,j;
	int nCoef;
	
	nCoef = (fac->nL+1)*(fac->nL+1);
	
	for (i=0;i<nCoef;i++) {
		for (j=0;j<nCoef;j++) {
			free(fac->fact[i][j]);
		}
		free(fac->fact[i]);
	}
	free(fac->fact);
	free(fac);
}




/*
GenerateNonZeros

Find the non-zero regions in a distribution function.  This is
used to speed up the propagator.  Especially duing the first 
iteration.
*/
NonZeros *GenerateNonZeros(Geometry *geom,Dist *dist) {
	int nCoef;
	int i,j,k,n;
	NonZeros *ret;
	int total,cnt;
	int bigtotal;
	
	nCoef = (dist->nL+1) * (dist->nL+1);
	ret = malloc(sizeof(NonZeros)*1);
	// One time through to count the non-zeros
	
	total = 0;
	bigtotal=0;
	for (n=0;n<nCoef;n++) {
		for (i=0;i<geom->nX;i++) {
			for (j=0;j<geom->nY;j++) {
				for (k=0;k<geom->nZ;k++) {
					bigtotal++;
					if (cabs(dist->dist[n][i+j*geom->nX + k *geom->nX*geom->nY]) > THRESHOLD) {
						total++;
					}
				}
			}
		}
	}
	//printf("fraction non zero = %f\n",(double)total/(double)bigtotal);
	ret->num = total;
	ret->is = malloc(sizeof(int)*total);
	ret->js = malloc(sizeof(int)*total);
	ret->ks = malloc(sizeof(int)*total);
	ret->cnts = malloc(sizeof(int)*total);
	
	cnt = 0;
	for (n=0;n<nCoef;n++) {
		for (i=0;i<geom->nX;i++) {
			for (j=0;j<geom->nY;j++) {
				for (k=0;k<geom->nZ;k++) {
					if (cabs(dist->dist[n][i+j*geom->nX + k *geom->nX*geom->nY]) > THRESHOLD) {
						ret->is[cnt] = i;
						ret->js[cnt] = j;
						ret->ks[cnt] = k;
						ret->cnts[cnt] = n;
						cnt++;
					} 
				}
			}
		}
	}

	
	return(ret);
}

/*
FreeNonZeros

Free up the memory used by the non-zero structure.
*/
void FreeNonZeros(NonZeros *zer) {
	
	free(zer->is);
	free(zer->js);
	free(zer->ks);
	free(zer->cnts);
	free(zer);
}

