#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iniparser.h"
#include "rte.h"

#pragma warning(disable:981)


/*
SpherHarmonicArray

Computes the spherical harmonic at N values of (theta,phi) and with degree
and order

Should be orthonormal 
*/
void SpherHarmonicArray(int degree,int order,int N,double *theta,double *phi,complex_double *out) {
	int tmporder;
	double nrm,fact;
	int i,j;
	complex_double ret;
	double *costheta,*Ymn;

	if ((costheta = (double *)malloc(sizeof(double)*N)) == NULL) {
		printf("Error allocating memory!\n");
		exit(1);
	}
	if ((Ymn = (double *) malloc(sizeof(double)*N)) == NULL) {
		printf("Error allocating memory!\n");
		exit(1);
	}
	
	if (order < 0) {
		tmporder = -order;
	} else {
		tmporder = order;
	}
	
	fact = 1.0;
	for (i=0;i < 2 * tmporder; i++) {
		fact = fact * (double)(tmporder + degree -i);
	}
	nrm = sqrt(((2.0*(double)degree + 1.0)/(4.0*M_PI * fact)));
	
	for (j=0;j<N;j++) {
		costheta[j] = cos(theta[j]);
	}
	AssocLegendreArray(degree,tmporder,N,costheta,Ymn);
	
	for (j=0;j<N;j++) {
		ret = nrm * Ymn[j] * cos((double)tmporder * phi[j]) + nrm *Ymn[j] * sin((double)tmporder * phi[j]) * I;
	
		if (order < 0) {
			ret = ~(ret);
			ret = ret / pow(-1.0,(double)tmporder);
		}
		out[j] = ret;
	}
	
	free(costheta);
	free(Ymn);
}

/*
SpherHarmonic

Computes the spherical harmonic at one value (theta,phi) and with degree
and order

Should be orthonormal 
*/
complex_double SpherHarmonic(double theta,double phi,int degree,int order) {
	int tmporder;
	complex_double ret;
	double Ymn,nrm,fact;
	int i;
	
	if (order < 0) {
		tmporder = -order;
	} else {
		tmporder = order;
	}
	
	fact = 1.0;
	for (i=0;i < 2 * tmporder; i++) {
		fact = fact * (double)(tmporder + degree -i);
	}
	nrm = sqrt(((2.0*(double)degree + 1.0)/(4.0*M_PI * fact)));
	
	Ymn = AssocLegendre(cos(theta),degree,tmporder);
	
	ret = nrm * Ymn * cos((double)tmporder * phi) + nrm *Ymn * sin((double)tmporder * phi) * I;
	
	if (order < 0) {
		ret = ~(ret);
		ret = ret / pow(-1.0,(double)tmporder);
	}
	
	return(ret);
}


/*
AssocLegendre

Compute the associated Legendre functions of degree and order at value x (between 0 and 1).
This uses recursion relationships found in "Abramowitz and Segun"
 
*/
double AssocLegendre(double x,int degree, int order) {
	int tmporder;
	int len;
	double old,new_var,tmp;
	double fact;
	int i;
	int l;
	
	if (order < 0) {
		tmporder = -order;
	} else {
		tmporder = order;
	}
	
	len = degree - tmporder + 2;
	old = 0.0;
	
	fact = 1.0;
	if (degree == 0) {
		new_var = 1.0;
	} else {
		for (i=(2*tmporder-1);i>0;i-=2) {
			fact *= (double)i;
		}
		new_var = fact*pow(-1.0,(double)tmporder);
		new_var = new_var*pow(1.0-pow(x,2.0),((double)tmporder/2.0));
	}
	
	l = tmporder + 1;
	for (i = 2;i < len;i++) {		
		tmp = x * (2.0 * (double)l - 1.0) * new_var/((double)(l-tmporder)) - (double)(l+tmporder-1)/((double)(l-tmporder))*old;
		old = new_var;
		new_var = tmp;
		l++;
	}
	
	if (order < 0) {
		fact = 1.0;
		for (i=0;i < 2 * tmporder; i++) {
			fact = fact * (double)(tmporder + degree -i);
		}
		tmp = pow(-1.0,(double)tmporder);
		tmp = tmp / fact * new_var;
		new_var = tmp;
	}
	
	return(new_var);
}

/*
AssocLegendreArray

Compute the associated Legendre functions of degree and order at all values in *in.  The
output is placed in *out.
N is the length of in and out.

This uses recursion relationships found in "Abramowitz and Segun"
 
*/

void AssocLegendreArray(int degree, int order,int N,double *in,double *out) {
	int tmporder;
	int len;
	double *old,*tmp;
	double fact;
	int i,j;
	int l;
	double fact1;
	
	if (order < 0) {
		tmporder = -order;
	} else {
		tmporder = order;
	}
	
	len = degree - tmporder + 2;
	
	if ((old = (double *) calloc(N,sizeof(double))) == NULL) {
		printf("Error allocating memory!\n");
		exit(1);
	}
	if ((tmp = (double *) malloc(sizeof(double)*N)) == NULL) {
		printf("Error allocating memory!\n");
		exit(1);
	}
		
	fact = 1.0;
	if (degree == 0) {
		for (j=0;j<N;j++) {
			out[j] = 1.0;
		}
	} else {
		for (i=(2*tmporder-1);i>0;i-=2) {
			fact *= (double)i;
		}
		for (j=0;j<N;j++) {
			out[j] = fact*pow(-1.0,(double)tmporder);
			out[j] = out[j]*pow(1.0-pow(in[j],2.0),((double)tmporder/2.0));
		}
	}
	
	l = tmporder + 1;
	for (i = 2;i < len;i++) {
		for (j=0;j<N;j++) {
			tmp[j] = in[j] * (2.0 * (double)l - 1.0) * out[j]/((double)(l-tmporder)) - (double)(l+tmporder-1)/((double)(l-tmporder))*old[j];
			old[j] = out[j];
			out[j] = tmp[j];
		}
		l++;
	}
	
	if (order < 0) {
		fact = 1.0;
		for (i=0;i < 2 * tmporder; i++) {
			fact = fact * (double)(tmporder + degree -i);
		}
		fact1 = pow(-1.0,(double)(tmporder));
		for (j=0;j<N;j++) {
			tmp[j] = (fact1 * tmp[j]) / fact * out[j];
			out[j] = tmp[j];
		}
	}
	
	free(old);
	free(tmp);
}

