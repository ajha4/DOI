#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <complex.h>

#include "iniparser.h"
#include "RTE.h"

#pragma warning(disable:981)

/* 
Output Intensity Image

*/
void OutputImage(Geometry *geom,complex double *transImg,char *name) {
	FILE *fle;
	int i,j;
	double val;

	fle = fopen(name,"w");
	
	fwrite(&geom->nX,sizeof(int),1,fle);
	fwrite(&geom->nY,sizeof(int),1,fle);
	
	for (j=0;j<geom->nY;j++) {
		for (i=0;i<geom->nX;i++) {
			val = creal(transImg[i+j*geom->nX]);
			fwrite(&val,sizeof(double),1,fle);
		}
	}
	for (j=0;j<geom->nY;j++) {
		for (i=0;i<geom->nX;i++) {
			val = cimag(transImg[i+j*geom->nX]);
			fwrite(&val,sizeof(double),1,fle);
		}
	}

	fclose(fle);
}


/* 
Output Phase-delay Image

*/
void OutputPhaseImage(Geometry *geom,double *transImg,char *name) {
	FILE *fle;
	int i,j;
	double val;

	fle = fopen(name,"w");
	
	fwrite(&geom->nX,sizeof(int),1,fle);
	fwrite(&geom->nY,sizeof(int),1,fle);
	
	for (j=0;j<geom->nY;j++) {
		for (i=0;i<geom->nX;i++) {
			val = transImg[i+j*geom->nX];
			fwrite(&val,sizeof(double),1,fle);
		}
	}

	fclose(fle);
}

/*
Compute Phase-delay Image

*/
double *ComputePhaseImage(Geometry *geom,Phantom *phan,Source *source,int nL,int num,Dist *W) {
	double *phse;
	complex double *ret;
	int l,i,j,k;
	int ip;
	double rad;
	complex double val;
	int NumXs = 100;
	double *xs;
	double *out;
	
	ret = malloc(sizeof(complex double) * geom->nX * geom->nY);
	phse = malloc(sizeof(double) * geom->nX * geom->nY);
	xs = malloc(NumXs*sizeof(double));
	out = malloc(NumXs*sizeof(double));
	for (i=0;i<NumXs;i++) {
		xs[i] = (double)i/((double)NumXs);
	}
	
	k = geom->nZ - 1;
	
	/* We have to account for the very first term which has 
	the source propagated.  We do this analytically.
	*/
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			val = 0 + 0*I;
			rad = sqrt((geom->xs[i] - source->posX) * (geom->xs[i] - source->posX) + 
				(geom->ys[j] - source->posY) * (geom->ys[j] - source->posY));
			if (rad <= source->diam/2.0) {
				val = val + source->mag * cexp(-(phan->mus+phan->mua) * (geom->zs[k]+geom->delZ/2.0+0*I));
			}
			
			for (l=0;l<= nL;l++) {
				AssocLegendreArray(l,0,NumXs,xs,out);
				for (ip = 0;ip<NumXs;ip++) {
					val = val + 2.0*M_PI*geom->delX * geom->delY * (phan->c/phan->n) / sqrt(4.0*M_PI) * sqrt(2.0*(double)l +1.0) *
					 	W->dist[l*l+l][i+j*geom->nX + k*geom->nX*geom->nY] *
						out[ip]*xs[ip] / ((double)NumXs - 1.0);
				}
			}
						
			ret[i+j*geom->nX] = val;
		}
	}
	
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			phse[i+j*geom->nX] = carg(ret[i+j*geom->nX]) * 180.0/M_PI;
		}
	}
	
	free(ret);
	free(xs);
	free(out);
	return(phse);
}

/*
ComputeReflImage

Compute the reflected image
*/
complex double *ComputeReflImage(Geometry *geom,Phantom *phan,int nL,int num,Dist *W) {
	complex double *ret;
	int l,i,j,k;
	int ip;
	complex double val;
	int NumXs = 100;
	double *xs;
	double *out;
	
	ret = calloc(geom->nX * geom->nY,sizeof(complex double));
	xs = malloc(NumXs*sizeof(double));
	out = malloc(NumXs*sizeof(double));
	for (i=0;i<NumXs;i++) {
		xs[i] = (double)i/((double)NumXs) - 1.0;
	}
	
	k = 0;
	
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			val = 0 + 0*I;			
			for (l=0;l<= nL;l++) {
				AssocLegendreArray(l,0,NumXs,xs,out);
				for (ip = 0;ip<NumXs;ip++) {
					val = val - 2.0*M_PI*geom->delX * geom->delY * (phan->c/phan->n) / sqrt(4.0*M_PI) * sqrt(2.0*(double)l +1.0) *
					 	W->dist[l*l+l][i+j*geom->nX + k*geom->nX*geom->nY] *
						out[ip]*xs[ip] / ((double)NumXs);
				}
			}
			
			ret[i+j*geom->nX] = val;
		}
	}
	
	free(xs);
	free(out);
	return(ret);
}

/*
ComputeTransImage

Compute the transmitted image
*/
complex double *ComputeTransImage(Geometry *geom,Phantom *phan,Source *source,int nL,int num,Dist *W) {
	complex double *ret;
	int l,i,j,k;
	int ip;
	double rad;
	complex double val;
	int NumXs = 100;
	double *xs;
	double *out;
	
	ret = calloc(geom->nX * geom->nY,sizeof(complex double));
	xs = malloc(NumXs*sizeof(double));
	out = malloc(NumXs*sizeof(double));
	for (i=0;i<NumXs;i++) {
		xs[i] = (double)i/((double)NumXs);
	}
	
	k = geom->nZ - 1;
	
	/* We have to account for the very first term which has 
	the source propagated.  We do this analytically.
	*/
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			val = 0 + 0*I;
			rad = sqrt((geom->xs[i] - source->posX) * (geom->xs[i] - source->posX) + 
				(geom->ys[j] - source->posY) * (geom->ys[j] - source->posY));
			//if (i == geom->nX/2 && j == geom->nY/2 ) {

			if (rad <= source->diam/2.0) {
				val = val + source->mag * cexp(-(phan->mus+phan->mua) * (geom->zs[k]+geom->delZ/2.0+0*I));
			}

			for (l=0;l<= nL;l++) {
				AssocLegendreArray(l,0,NumXs,xs,out);
				for (ip = 0;ip<NumXs;ip++) {
					val = val + 2.0*M_PI*geom->delX * geom->delY * (phan->c/phan->n) / sqrt(4.0*M_PI) * sqrt(2.0*(double)l +1.0) *
					 	W->dist[l*l+l][i+j*geom->nX + k*geom->nX*geom->nY] *
						out[ip]*xs[ip] / ((double)NumXs);
				}
			}
			//val = val + 2.0*M_PI*geom->delX * geom->delY * (phan->c/phan->n) / sqrt(4.0*M_PI) * ((0.5*W->dist[0][i+j*geom->nX + k*geom->nX*geom->nY]) + 
			//	((sqrt(3.0)/3.0) * W->dist[2][i+j*geom->nX + k*geom->nX*geom->nY]));
			
			ret[i+j*geom->nX] = val;
		}
	}

	free(xs);
	free(out);
	return(ret);
}


/*
ComputeTransImageNoSource

Compute the transmitted image

complex double *ComputeTransImageNoSource(Geometry *geom,Phantom *phan,int nL,int num,Dist *W) {
	complex double *ret;
	int l,i,j,k;
	int ip;
	complex double val;
	int NumXs = 100;
	double *xs;
	double *out;
	
	ret = calloc(geom->nX * geom->nY,sizeof(complex double));
	xs = malloc(NumXs*sizeof(double));
	out = malloc(NumXs*sizeof(double));
	for (i=0;i<NumXs;i++) {
		xs[i] = (double)i/((double)NumXs - 1.0);
	}

	k = geom->nZ - 1;
	
	// We have to account for the very first term which has 
	// the source propagated.  We do this analytically.
	
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			val = 0 + 0*I;			
			for (l=0;l<= nL;l++) {
				AssocLegendreArray(l,0,NumXs,xs,out);
				for (ip = 0;ip<NumXs;ip++) {
					val = val + 2.0*M_PI*geom->delX * geom->delY * (phan->c/phan->n) / sqrt(4.0*M_PI) * sqrt(2.0*(double)l +1.0) *
					 	W->dist[l*l+l][i+j*geom->nX + k*geom->nX*geom->nY] *
						out[ip]*xs[ip] / ((double)NumXs - 1.0);
				}
			}			
			ret[i+j*geom->nX] = val;
		}
	}
	
	free(xs);
	free(out);
	return(ret);
}
*/

/*
ComputeTransImage

Compute the transmitted image

complex double *ComputeTransImageLW(Geometry *geom,Phantom *phan,Source *source,int nL,int num,Dist *W) {
	complex double *ret;
	complex double **tmp;
	double *theta;
	double *phi;
	int l,m,i,j,k;
	int ip,jp;
	int cnt;
	int nCoef;
	complex double val;
	double delphi,deltheta;
	
	ret = calloc(geom->nX * geom->nY,sizeof(complex double));
	
	k = geom->nZ - 1;
	nCoef = (nL+1)*(nL+1);
	tmp = malloc(sizeof(complex double*)*nCoef);
	for (m=0;m<nCoef;m++) {
		tmp[m] = malloc(sizeof(complex double)*num*num);
	}
	
	theta = malloc(sizeof(double) * num*num);
	phi = malloc(sizeof(double)*num*num);
	deltheta = (M_PI/2.0)/(double)(num);
	delphi = (2.0*M_PI)/(double)(num);
	for (i=0;i<num;i++) {
		for (j=0;j<num;j++) {
			theta[i+j*num] = (double)i*(M_PI/2.0)/(double)(num);
			phi[i+j*num] = (double)j*(M_PI*2.0)/(double)(num);
		}
	}
	
	cnt = 0;
	for (l=0;l<=nL;l++) {
		for (m=-l;m<=l;m++) {
			SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
			cnt++;
		}
	}
	
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			val = 0 + 0*I;
			for (ip = 0;ip < num;ip++) {
				for (jp=0;jp<num;jp++) {
					for (m=0;m<nCoef;m++) {
						val = val + W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] *
							tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi;							
					}
				}
			}
			ret[i+j*geom->nX] = val * geom->delX * geom->delY * phan->c;
		}
	}
		
	for (m=0;m<nCoef;m++){
		free(tmp[m]);
	}
	free(tmp);
	free(theta);
	free(phi);
	return(ret);
}
*/

/*
Compute and output the side wall images

*/
void ComputeSideFaceImages(Geometry *geom,Phantom *phan,int nL, int num, Dist *W){
	complex double *sideImg;
	complex double **tmp;
	double *theta;
	double *phi;
	int l,m,i,j,k;
	int ip,jp;
	int cnt;
	int nCoef;
	complex double val,val_tmp;
	double delphi,deltheta;

	nCoef = (nL+1)*(nL+1);
	tmp = malloc(sizeof(complex double*)*nCoef);
	for (m=0;m<nCoef;m++) {
		tmp[m] = malloc(sizeof(complex double)*num*num);
	}
	theta = malloc(sizeof(double) * num*num);
	phi = malloc(sizeof(double)*num*num);
	sideImg = calloc(geom->nZ * geom->nY,sizeof(complex double));

	// Side Face 1 (-X)
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

	i = 0;
	for (k=0;k<geom->nZ;k++) {
		for (j=0;j<geom->nY;j++) {
          		val = 0 + 0*I;
           		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp +  geom->delZ * geom->delY * phan->c/phan->n*
                            			tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi * sin(theta[ip+jp*num]) * cos(phi[ip+jp*num]) * (-1);	// Along the negative x axis
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
			}
           		sideImg[k+j*geom->nZ] = val;
        	}
    	}
	// Output image
	OutputImage(geom, sideImg, "SideFace1(-X).out");

	// Side Face 2 (+X)
        deltheta = (M_PI)/(double)(num);
        delphi = (M_PI)/(double)(num);
        for (i=0;i<num;i++) {
                for (j=0;j<num;j++) {
                        theta[i+j*num] = (double)i*(M_PI)/(double)(num);
                        phi[i+j*num] = (double)j*(M_PI)/(double)(num) - (M_PI/2);
                }
        }
        cnt = 0;
        for (l=0;l<=nL;l++) {
                for (m=-l;m<=l;m++) {
                        SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
                        cnt++;
                }
        }

        i = geom->nX-1;
        for (k=0;k<geom->nZ;k++) {
                for (j=0;j<geom->nY;j++) {
                        val = 0 + 0*I;
                        for (m=0;m<nCoef;m++) {
                                val_tmp = 0 + 0*I;
                                for (ip = 0;ip < num;ip++) {
                                        for (jp=0;jp<num;jp++) {
                                                val_tmp = val_tmp +  geom->delZ * geom->delY * phan->c/phan->n*
                                                tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi * sin(theta[ip+jp*num]) * cos(phi[ip+jp*num]);	// Along the positive x axis
                                        }
                                }
                                val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
                        }
                        sideImg[k+j*geom->nZ] = val;
                }
        }
        // Output image
        OutputImage(geom, sideImg, "SideFace2(+X).out");
	free(sideImg);

	// Side Face 3 (-Y)
        sideImg = calloc(geom->nZ * geom->nX,sizeof(complex double));
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

        j = 0;
        for (i=0;i<geom->nX;i++) {
                for (k=0;k<geom->nZ;k++) {
                        val = 0 + 0*I;
                        for (m=0;m<nCoef;m++) {
                                val_tmp = 0 + 0*I;
                                for (ip = 0;ip < num;ip++) {
                                        for (jp=0;jp<num;jp++) {
						val_tmp = val_tmp +  geom->delZ * geom->delX * phan->c/phan->n*
                            			tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi *sin(theta[ip+jp*num]) * sin(phi[ip+jp*num]) * (-1);		// Along the negative y axis
                                        }
                                }
                                val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
                        }
                        sideImg[k+i*geom->nZ] = val;
                }
        }
        // Output image
        OutputImage(geom, sideImg, "SideFace3(-Y).out");

	// Side Face 4 (+Y)
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

        j = geom->nY - 1;
        for (i=0;i<geom->nX;i++) {
                for (k=0;k<geom->nZ;k++) {
                        val = 0 + 0*I;
                        for (m=0;m<nCoef;m++) {
                                val_tmp = 0 + 0*I;
                                for (ip = 0;ip < num;ip++) {
                                        for (jp=0;jp<num;jp++) {
                                                val_tmp = val_tmp +  geom->delZ * geom->delX * phan->c/phan->n*
                                                tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi *sin(theta[ip+jp*num]) * sin(phi[ip+jp*num]); 	// Along the positive y axis
                                        }
                                }
                                val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
                        }
                        sideImg[k+i*geom->nZ] = val;
                }
        }
        // Output image
        OutputImage(geom, sideImg, "SideFace4(+Y).out");

	// Clean up	
	for (m=0;m<nCoef;m++){
		free(tmp[m]);
	}
	free(tmp);
	free(theta);
	free(phi);
	free(sideImg);
}

/*
ComputePhotonNumber

*/
complex double ComputePhotonNumber(Geometry *geom,Phantom *phan,Source *source,int nL,int num,Dist *W) {
	complex double ret,ret_tmp,ret_in, ret_tmp_in;
	complex double **tmp;
	double *theta;
	double *phi;
	int l,m,i,j,k;
	int ip,jp;
	int cnt;
	int nCoef;
	double rad;
	complex double val;
	double delphi,deltheta;
	
	ret = 0 + 0*I;
    	ret_in = 0 + 0*I;
	ret_tmp_in = 0 + 0*I;
	
	nCoef = (nL+1)*(nL+1);
	tmp = malloc(sizeof(complex double*)*nCoef);
	for (m=0;m<nCoef;m++) {
		tmp[m] = malloc(sizeof(complex double)*num*num);
	}
	
	theta = malloc(sizeof(double) * num*num);
	phi = malloc(sizeof(double)*num*num);

   	complex double val_tmp;
	// Bottom Face
	k = geom->nZ - 1;
	deltheta = (M_PI/2.0)/(double)(num);
	delphi = (2.0*M_PI)/(double)(num);
	for (i=0;i<num;i++) {
		for (j=0;j<num;j++) {
			theta[i+j*num] = (double)i*(M_PI/2.0)/(double)(num);
			phi[i+j*num] = (double)j*(M_PI*2.0)/(double)(num);
		}
	}
	cnt = 0;
	for (l=0;l<=nL;l++) {
		for (m=-l;m<=l;m++) {
			SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
			cnt++;
		}
	}

   	complex double trans = 0 + 0*I;
	printf("Bottom face outgoing flux \n");
//    	FILE *fp;
//    	fp = fopen("bottomface.dat","w");

	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			val = 0 + 0*I;
			rad = sqrt((geom->xs[i] - source->posX) * (geom->xs[i] - source->posX) + 
				(geom->ys[j] - source->posY) * (geom->ys[j] - source->posY));
			if (rad <= source->diam/2.0) {
			//if (i == geom->nX/2 && j == geom->nY/2 ) {
				val = val + source->mag * cexp(-(phan->mus+phan->mua) * (geom->zs[k]+geom->delZ/2.0+0*I));
				//printf("%e is XE term \n", creal(val));
				trans = val;
			}
			for (m=0;m<nCoef;m++) {
				val_tmp = 0 + 0*I;
				for (ip = 0;ip < num; ip++) {
					for (jp=0;jp < num; jp++) {
						val_tmp = val_tmp + 
							tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi *cos(theta[ip+jp*num]) * geom->delX * geom->delY * phan->c /phan->n;							
					}
				}
				val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
				if(i == geom->nX/2 && j == geom->nY/2){// && (m == 0 || m ==2 || m== 6 || m == 12)){
		    		printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
				}
//				if(m == 0 || m==2){
//					fwrite(&creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), sizeof(double), 1, fp);
//				}
//				if(creal(val_tmp) < 0.0)
//    					printf("Contribution from bottom face is %e for %d spherical harmonic for (%d %d) voxel \n", creal(val_tmp), m,i,j );
			}
			ret += val; 
		}
	}
//    	fclose(fp);
	printf("Contribution from bottom face is %e \n", creal(ret));

	ret_tmp = ret;
    	printf("Top face incoming flux: \n");
	k = 0;
    	//fp = fopen("topface.dat","w");

    	for (i=0;i<geom->nX;i++) {
        	for (j=0;j<geom->nY;j++) {
            		val = 0 + 0*I;
		    	for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num; ip++) {
                    			for (jp=0;jp < num; jp++) {
                        			val_tmp = val_tmp + tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi *cos(theta[ip+jp*num]) * geom->delX * geom->delY * phan->c /phan->n;   
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
                		if(i == geom->nX/2 && j == geom->nY/2){// && (m == 0 || m ==2)){
	                		printf(" A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
 //               			fwrite(&creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), sizeof(double),1,fp);
				}
            		}
			ret_in +=val;
        	}
    	}
    	//fclose(fp);
	printf("Contribution from top face incoming is %e \n", creal(ret_in - ret_tmp_in));

    	ret_tmp_in = ret_in;

	// Top Face
	k = 0;
    	printf("Top face outgoing flux: \n");

	deltheta = (M_PI/2.0)/(double)(num);
	delphi = (2.0*M_PI)/(double)(num);
	for (i=0;i<num;i++) {
		for (j=0;j<num;j++) {
			theta[i+j*num] = (double)i*(M_PI/2.0)/(double)(num) + (M_PI/2.0);
			phi[i+j*num] = (double)j*(M_PI*2.0 )/(double)(num) ;
		}
	}
	
	cnt = 0;
	for (l=0;l<=nL;l++) {
		for (m=-l;m<=l;m++) {
			SpherHarmonicArray(l,m,num*num,theta,phi,tmp[cnt]);
			cnt++;
		}
	}

/*
	printf("\n \n ********************* \n \n");

	m = 1;
	ip = 0;
    	fp = fopen("sph_harm.dat", "w");
			val_tmp = 0 + 0*I;
			for (jp=0;jp<num;jp++) {
				val_tmp = val_tmp +  tmp[m][ip+jp*num]; // * sin(theta[ip+jp*num])*deltheta *delphi ; // *cos(theta[ip+jp*num]);
				//printf("sph_harm is %e + %e i for theta = %f phi = %f \n", creal(tmp[m][ip+jp*num]), cimag(tmp[m][ip+jp*num]),theta[ip+jp*num],phi[ip+jp*num]);
		    	fwrite(&(creal(tmp[m][ip+jp*num])), sizeof(double), 1,fp);
			    fwrite(&(cimag(tmp[m][ip+jp*num])), sizeof(double), 1,fp);
			}
		}
	}
	printf("\n \n ********************* \n \n");
   	fclose(fp);
*/
	for (i=0;i<geom->nX;i++) {
		for (j=0;j<geom->nY;j++) {
			val = 0 + 0*I;
			for (m=0;m<nCoef;m++) {
				val_tmp = 0 + 0*I;
				for (ip = 0;ip < num;ip++) {
					for (jp=0;jp<num;jp++) {
						val_tmp = val_tmp + geom->delX * geom->delY * phan->c/phan->n * tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta *delphi *cos(theta[ip+jp*num])*(-1);//-1 because nhat.shat = -cos(theta)						
					}
				}
				val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
				if(i == geom->nX/2 && j == geom->nY/2)// && ( m == 0 || m ==2))
	    				printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
			}
			ret += val;  // divide by n??
		}
	}
    	printf("Contribution from top face outgoing is %e \n", creal(ret - ret_tmp));
    	ret_tmp = ret;
	
	k = geom->nZ - 1;

	printf("Bottom face incoming flux: \n");
    	for (i=0;i<geom->nX;i++) {
        	for (j=0;j<geom->nY;j++) {
            		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num; ip++) {
                    			for (jp=0;jp < num; jp++) {
                        			val_tmp = val_tmp + tmp[m][ip+jp*num] * sin(theta[ip+jp*num])*deltheta *delphi * (-1) *cos(theta[ip+jp*num]) * geom->delX * geom->delY * phan->c /phan->n;   
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
                		if(i == geom->nX/2 && j == geom->nY/2){// && (m == 0 || m ==2 || m== 6 || m == 12)){
           //     			if(m==0 || m ==2 ){
	    				printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
				}
            		}
            		ret_in +=val;
        	}
    	}
	printf("Contribution from bottom face incoming is %e \n", creal(ret_in - ret_tmp_in));
    	ret_tmp_in = ret_in;

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
	
	i = 0;
	for (k=0;k<geom->nZ;k++) {
		for (j=0;j<geom->nY;j++) {
          		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp+ geom->delZ* geom->delY* phan->c/phan->n*tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta *delphi* sin(theta[ip+jp*num])* cos(phi[ip+jp*num])*(-1);// Along the negative x axis
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
            			if(k == geom->nZ/2 && j == geom->nY/2)// && ( m == 0 || m ==2))
                			printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
            		}
            		ret += val;  // divide by n??
        	}
    	}
    	printf("Outward flux normal to side face 1 along negative x axis  is %e \n", creal(ret - ret_tmp));
    	ret_tmp = ret;

	i = geom->nX-1;
    	for (k=0;k<geom->nZ;k++) {
        	for (j=0;j<geom->nY;j++) {
          		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp+ geom->delZ* geom->delY* phan->c/phan->n*tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta *delphi *sin(theta[ip+jp*num])* cos(phi[ip+jp*num])*(-1);// Along the negative x axis;                     
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
            			if(k == geom->nZ/2 && j == geom->nY/2 )//&& ( m == 0 || m ==2))
                			printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
            		}
            		ret_in += val;  // divide by n??
        	}
    	}
    	printf("Inward flux normal to side face 2 along negative x axis  is %e \n", creal(ret_in - ret_tmp_in));
    	ret_tmp_in = ret_in;

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

	i = geom->nX-1;
	for (k=0;k<geom->nZ;k++) {
		for (j=0;j<geom->nY;j++) {
          		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp+ geom->delZ* geom->delY* phan->c/phan->n*tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta *delphi*sin(theta[ip+jp*num])* cos(phi[ip+jp*num]);// Along the negative x axis                    
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
            			if(k == geom->nZ/2 && j == geom->nY/2)// && ( m == 0 || m ==2))
                			printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
            		}
            		ret += val;  // divide by n??
        	}
    	}
    	printf("Outward flux normal to side face 2 along positive x axis is %e \n", creal(ret - ret_tmp));
    	ret_tmp = ret;

	i = 0;
    	for (k=0;k<geom->nZ;k++) {
        	for (j=0;j<geom->nY;j++) {
          		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp+ geom->delZ* geom->delY* phan->c/phan->n*tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta*delphi* sin(theta[ip+jp*num])* cos(phi[ip+jp*num]);// Along the negative x axis                    
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
            			if(k == geom->nZ/2 && j == geom->nY/2)// && ( m == 0 || m ==2))
                			printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
            		}
            		ret_in += val;  // divide by n??
        	}
    	}
    	printf("Inward flux normal to side face 1 along positive x axis is %e \n", creal(ret_in - ret_tmp_in));
    	ret_tmp_in = ret_in;

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

	j = 0;
	for (i=0;i<geom->nX;i++) {
		for (k=0;k<geom->nZ;k++) {
          		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp+ geom->delX* geom->delZ* phan->c/phan->n*tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta*delphi*sin(theta[ip+jp*num])* sin(phi[ip+jp*num])*(-1);// Along the negative y axis                    
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
            			if(i == geom->nX/2 && k == geom->nZ/2)// && ( m == 0 || m ==2))
                			printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
            		}
            		ret += val;  // divide by n??
        	}
    	}
    	printf("Outward flux normal to side face 3 along negative y axis  is %e \n", creal(ret - ret_tmp));
    	ret_tmp = ret;

   	j = geom->nY-1;
   	for (i=0;i<geom->nX;i++) {
        	for (k=0;k<geom->nZ;k++) {
          		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp+ geom->delX* geom->delZ* phan->c/phan->n*tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta*delphi *sin(theta[ip+jp*num])* sin(phi[ip+jp*num])*(-1);// Along the negative y axis          
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
            			if(i == geom->nX/2 && k == geom->nZ/2)// && ( m == 0 || m ==2))
                			printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
            		}
            		ret_in += val;  // divide by n??
        	}
    	}
    	printf("Inward flux normal to side face 4 along negative y axis  is %e \n", creal(ret_in - ret_tmp_in));
    	ret_tmp_in = ret_in;

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
	
	j = geom->nY - 1;
	for (i=0;i<geom->nX;i++) {
		for (k=0;k<geom->nZ;k++) {
          		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp+ geom->delX* geom->delZ* phan->c/phan->n*tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta*delphi*sin(theta[ip+jp*num])* sin(phi[ip+jp*num]);// Along the positive y axis         
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
            			if(i == geom->nX/2 && k == geom->nZ/2)// && ( m == 0 || m ==2))
                			printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
            		}
            		ret += val;  // divide by n??
        	}
    	}
    	printf("Outward flux normal to side face 4 along positive y axis is %e \n", creal(ret - ret_tmp));
    	ret_tmp = ret;

	j = 0;
    	for (i=0;i<geom->nX;i++) {
        	for (k=0;k<geom->nZ;k++) {
          		val = 0 + 0*I;
            		for (m=0;m<nCoef;m++) {
                		val_tmp = 0 + 0*I;
                		for (ip = 0;ip < num;ip++) {
                    			for (jp=0;jp<num;jp++) {
                        			val_tmp = val_tmp+ geom->delX* geom->delZ* phan->c/phan->n*tmp[m][ip+jp*num]* sin(theta[ip+jp*num])*deltheta*delphi*sin(theta[ip+jp*num])* sin(phi[ip+jp*num]);// Along the positive y axis        
                    			}
                		}
                		val = val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY] + val;
            			if(i == geom->nX/2 && k == geom->nZ/2)
                			printf("A(%d) = %e,w(%d)= %e. F(%d) = %e \n", m,creal(val_tmp),m,creal(W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]), m,creal(val_tmp*W->dist[m][i+j*geom->nX + k*geom->nX*geom->nY]) );
            		}
            		ret_in += val;  // divide by n??
        	}
    	}
    	printf("Inward flux normal to side face 3 along positive y axis is %e \n", creal(ret_in - ret_tmp_in));
    	ret_tmp_in = ret_in;

	// Clean up	
	for (m=0;m<nCoef;m++){
		free(tmp[m]);
	}
	free(tmp);
	free(theta);
	free(phi);
	
    	printf("Transmitted flux is %e \n",creal(trans) );
    	printf("Total outward flux is %e \n", creal(ret));
    	printf("Total inward flux is %e \n", creal(ret_in));
	
	return(ret + ret_in);
}


/*
OutputDist

Output the distribution function to a file for viewing in Matlab
*/
void OutputDist(Geometry *geom,Dist *dist,char *name) {
	FILE *fle;
	int nCoef;
	int i,j,k,n;
	double val;
	
	nCoef = (dist->nL+1) * (dist->nL+1);
	
	fle = fopen(name,"w");
	
	fwrite(&geom->nX,sizeof(int),1,fle);
	fwrite(&geom->nY,sizeof(int),1,fle);
	fwrite(&geom->nZ,sizeof(int),1,fle);
	fwrite(&dist->nL,sizeof(int),1,fle);
	
	for (n=0;n<nCoef;n++) {
		for (k=0;k<geom->nZ;k++) {
			for (j=0;j<geom->nY;j++) {
				for (i=0;i<geom->nX;i++) {
					val = creal(dist->dist[n][i+j*geom->nX + k*geom->nY*geom->nX]);
					fwrite(&val,sizeof(double),1,fle);
				}
			}
		}
	}
	
	for (n=0;n<nCoef;n++) {
		for (k=0;k<geom->nZ;k++) {
			for (j=0;j<geom->nY;j++) {
				for (i=0;i<geom->nX;i++) {
					val = cimag(dist->dist[n][i+j*geom->nX + k*geom->nY*geom->nX]);
					fwrite(&val,sizeof(double),1,fle);
				}
			}
		}
	}
	
	fclose(fle);
}

/* 
GenerateSPECTSource

Read in the paramater file and generate the source distribution function.

void GenerateSPECTSource(Geometry *geom,Phantom *phan,Source *source,int nL,Dist *src) {
	int i,j,k;
	int l,m;
	int cnt;
	double rad;
	int total;
					
	total = 0;
	for (k=0;k<geom->nZ;k++) {
		for (i=0;i<geom->nX;i++) {
			for (j=0;j<geom->nY;j++) {
				rad = sqrt((geom->xs[i] - source->posX) * (geom->xs[i] - source->posX) + 
					(geom->ys[j] - source->posY) * (geom->ys[j] - source->posY) + 
					(geom->zs[k] - source->posZ) * (geom->zs[k] - source->posZ));
				if (rad <= source->diam/2.0) {
					total++;
					cnt = 0;
					for (l=0;l<=nL;l++) {
						for (m=-l;m<=l;m++) {
							src->dist[cnt][i + j*geom->nX + k*geom->nY*geom->nX] =  source->mag * pow(phan->g,(double)l);
							cnt++;
						}
					}
				}
			}
		}
	}
}
*/

/* 
   GenerateSourceBeam

Read in ithe paramater file and generate the source distribution function.
*/
void GenerateSourceBeam(Geometry *geom,Phantom *phan,Source *source,int nL,Dist *src) {
	int i,j,k;
	int l,m;
	int cnt;
	double rad;
	double fact;
	int total;
	double tmpz;

	total = 0;
	for (k=0;k<geom->nZ;k++) {
		for (i=0;i<geom->nX;i++) {
			for (j=0;j<geom->nY;j++) {
				rad = sqrt((geom->xs[i] - source->posX) * (geom->xs[i] - source->posX) + 
					(geom->ys[j] - source->posY) * (geom->ys[j] - source->posY));
				if (rad <= source->diam/2.0) {
					total++;
					cnt = 0;
					for (l=0;l<=nL;l++) {
						for (m=-l;m<=l;m++) {
							tmpz = geom->zs[k]-geom->delZ/2.0;
							if (m==0) {
								src->dist[cnt][i + j*geom->nX + k*geom->nY*geom->nX] =
									pow(phan->g,(double)l) * (phan->mus/(phan->mus+phan->mua)) *  sqrt(((2.0*(double)l)+1.0)/(4.0*M_PI))
									 * (cexp(-(phan->mus+phan->mua)*tmpz) - cexp(-(phan->mus+phan->mua)*(tmpz+geom->delZ)))
									/ (geom->delX * geom->delY * geom->delZ); 
							}
							cnt++;
						}
					}
				}
			}
		}
	}

	// Normalize the total number of incident photons
	fact = source->mag / ((double)total/(double)geom->nZ);
	source->mag = fact;
	for (k=0;k<geom->nZ;k++) {
		for (i=0;i<geom->nX;i++) {
			for (j=0;j<geom->nY;j++) {
				cnt = 0;
				for (l=0;l<=nL;l++) {
					for (m=-l;m<=l;m++) {
						src->dist[cnt][i+j*geom->nX + k*geom->nY*geom->nX] *= fact;
						cnt++;
					}
				}
			}
		}
	}
}


/* 
GenerateSource

Read in the paramater file and generate the source distribution function.

void GenerateSource(Geometry *geom,Source *source,int nL,Dist *src) {
	int i,j,k;
	int l,m;
	int cnt;
	double rad;
	double theta,phi;
	double fact;
	int total;
					
	total = 0;
	for (k=0;k<geom->nZ;k++) {
		if (geom->zs[k] - source->posZ <= 1e-20) {
			for (i=0;i<geom->nX;i++) {
				for (j=0;j<geom->nY;j++) {
					rad = sqrt((geom->xs[i] - source->posX) * (geom->xs[i] - source->posX) + 
						(geom->ys[j] - source->posY) * (geom->ys[j] - source->posY));
					if (rad <= source->diam/2.0) {
						total++;
						cnt = 0;
						for (l=0;l<=nL;l++) {
							for (m=-l;m<=l;m++) {
								src->dist[cnt][i+j*geom->nX + k*geom->nY*geom->nX] = conj(SpherHarmonic(theta,phi,l,m));								
								cnt++;
							}
						}
					}
				}
			}
		}
	}
	
	// We must ensure that the total number of incident photons is normalized for each study 
	// Otherwise, a change in geometry will result in a change in results 
	fact = source->mag / ((double)total * geom->delX * geom->delY * geom->delZ);
	for (k=0;k<geom->nZ;k++) {
		if (geom->zs[k] - source->posZ <= 1e-20) {
			for (i=0;i<geom->nX;i++) {
				for (j=0;j<geom->nY;j++) {
					cnt = 0;
					for (l=0;l<=nL;l++) {
						for (m=-l;m<=l;m++) {
							src->dist[cnt][i+j*geom->nX + k*geom->nY*geom->nX] *= fact;
							cnt++;
						}
					}
				}
			}
		}
	}
}
*/
	

/*
AllocDist

Allocate space for a distribution function.  The representation is in spherical
harmonics.  The returned distribution function has all zeros.
*/
Dist *AllocDist(Geometry *geom,int nL) {
	int i,nCoef;
	Dist *dst;
	
	nCoef = (nL+1) * (nL+1);
	
	if ((dst = malloc(sizeof(Dist)*1)) == NULL) {
		printf("Error allocating memory!\n");
		exit(1);
	}
	dst->nL = nL;
	
	if ((dst->dist = malloc(sizeof(complex double *)*nCoef)) == NULL) {
		printf("Error allocating memory!\n");
		exit(1);
	}
	for (i=0;i<nCoef;i++) {
		if ((dst->dist[i] = calloc((geom->nX * geom->nY * geom->nZ),sizeof(complex double))) == NULL) {
			printf("Error allocating memory!\n");
			exit(1);
		}
	}
	return dst;
}


/*
FreeDist

Free up the space allocated by a particular distribution function.
*/
void FreeDist(Dist *dst) {
	int nCoef,i;
	
	nCoef = (dst->nL+1) * (dst->nL+1);
	
	for (i=0;i<nCoef;i++) {
		free(dst->dist[i]);
	}
	free(dst->dist);
	free(dst);
	
}

void AddDist(Geometry *geom,int nL,Dist *W,Dist *out) {
	int i,n,nCoef;
	
	nCoef = (nL+1)*(nL+1);
	for (n=0;n<nCoef;n++) {
		for (i=0;i< geom->nX*geom->nY*geom->nZ;i++) {
			out->dist[n][i] += W->dist[n][i];
		}			
	}
}

void ScaleDist(Geometry *geom,int nL,Dist *W,double factor) {
	int i,n,nCoef;
	
	nCoef = (nL+1)*(nL+1);
	for (n=0;n<nCoef;n++) {
		for (i=0;i< geom->nX*geom->nY*geom->nZ;i++) {
			W->dist[n][i] *= factor;
		}			
	}
}

void DotDist(Geometry *geom,int nL,Dist *W,Dist *W1,Dist *out) {
	int i,n,nCoef;
	
	nCoef = (nL+1)*(nL+1);
	for (n=0;n<nCoef;n++) {
		for (i=0;i< geom->nX*geom->nY*geom->nZ;i++) {
			out->dist[n][i] = conj(W->dist[n][i]) * W1->dist[n][i];
		}			
	}
}

void SubDist(Geometry *geom,int nL,Dist *W,Dist *W1,Dist *out) {
	int i,n,nCoef;
	
	nCoef = (nL+1)*(nL+1);
	for (n=0;n<nCoef;n++) {
		for (i=0;i< geom->nX*geom->nY*geom->nZ;i++) {
			out->dist[n][i] = W->dist[n][i] - W1->dist[n][i];
		}			
	}
}

void CopyDist(Geometry *geom,int nL,Dist *W,Dist *out) {
	int i,n,nCoef;
	
	nCoef = (nL+1)*(nL+1);
	for (n=0;n<nCoef;n++) {
		for (i=0;i< geom->nX*geom->nY*geom->nZ;i++) {
			out->dist[n][i] = W->dist[n][i];
		}			
	}
}

void CombineOuts(Geometry *geom,int nL,Dist *out1,Dist *out2) {
	int i,n,nCoef;
	
	nCoef = (nL+1)*(nL+1);
	for (n=0;n<nCoef;n++) {
		for (i=0;i< geom->nX*geom->nY*geom->nZ;i++) {
			out2->dist[n][i] =  0.5*out1->dist[n][i] + 0.25*out2->dist[n][i] + 0.25 * conj(out2->dist[n][i]);
		}			
	}
}

void AdjModify(Geometry *geom, int nL,Dist *W) {
	int l,m,n,i;
	
	double fact;
	n = 0;
	for (l=0;l<nL;l++) {
		for (m=-l;m<=l;m++) {
			fact = pow(-1.0,(double)m);
			if (fact != 1.0) {
				for (i=0;i< geom->nX*geom->nY*geom->nZ;i++) {
					W->dist[n][i] = fact * W->dist[n][i];
				}
			}
			n++;
		}
	}
}

