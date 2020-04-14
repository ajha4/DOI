/*This code processes the .his file generated using tmcimg to produce the image. 
The code was originally provided by Matthew Kupinski, PhD, University of Arizona, and was modified by Abhinav Jha, PhD. 
Please cite the following references when you use this code. 
1. A. K. Jha, M. A. Kupinski, H. H. Barrett, E. Clarkson, J. H. Hartman, “A three-dimensional Neumann-series approach to model light transport in non-uniform media”, J. Opt. Soc. Amer. A, 29(8), 1885-99, 2012 PMCID: PMC3963433 (link)
2. A. K. Jha, M. A. Kupinski, T. Masumura, E. Clarkson, A. A. Maslov, H. H. Barrett, “Simulating photon transport in uniform media using the radiative transport equation: A study using the Neumann-series approach”, J. Opt. Soc. Amer. A, 29(8), 1741-1757, 2012. (Top 10 downloaded articles in Aug. 2012 edition.) PMCID: PMC3985394 (link)
3. A. K. Jha, E. Clarkson, and M. A. Kupinski, “An ideal-observer framework to investigate signal detectability in diffuse optical imaging”, Biomed. Opt. Express, 4(10), 2107-23, 2013 PMCID: PMC3799670 (link)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Will make this a variable later
#define DIM 51
#define BUFSIZE 1000000

int main(int argc,char **argv) {
	int n;
	int i,j;
	int ii;
	FILE *fle;
	FILE *outfle;
	double *mus;
	double *img;
	float *buffer;
	int nmus;
	int packet_size;
	size_t numread;
	double x;
	double y;
	double z;
	double len;
	double delta;
	
	// Check the arguments
	if (argc < 4) {
		printf("ProcessHisFil filename [mua1 mua2 ...] outputname\n");
		exit(1);
	}
	
	// Open the file
	fle = fopen(argv[1],"r");
	if (fle == NULL) {
		printf("Cannot open file %s\n",argv[1]);
		exit(1);
	}
	
	nmus = argc - 3;  // 3 because of program name, and two file arguments
	mus = malloc(sizeof(double)*nmus);
	
	
	// Read the mu_a coefficient from the command line
	for (i = 2;i<2+nmus;i++) {
		mus[i] = atof(argv[i]);
		printf("%f ",mus[i]);
	}
	printf("\n");
	
	packet_size = nmus+4; // 4 includes x,y,z,t
	
	delta = 100.0/((double)DIM);
	
	img = malloc(sizeof(double)*DIM*DIM);
	for (n = 0;n<DIM*DIM;n++){ img[n] = 0.0;}
	buffer = malloc(sizeof(float)*BUFSIZE*packet_size);
	
	do {
		numread = fread(buffer,sizeof(float),BUFSIZE*packet_size,fle);
		for (n = 0;n<(numread/packet_size);n++) {
			z = (double)buffer[n*packet_size + 2];
			if (z == 42.0) {
				x = (double)buffer[n*packet_size + 0];
				y = (double)buffer[n*packet_size + 1];
				i = (int)floor((x-2.0)/delta);
				j = (int)floor((y-2.0)/delta);
				if (i<0) i=0;
				if (i>(DIM-1)) i=DIM-1;
				if (j<0) j=0;
				if (j>(DIM-1)) j=DIM-1;
				for (ii = 0;ii<nmus;ii++) {
					len = (double)buffer[n*packet_size + ii+4];
					img[j*DIM+i] += exp(-mus[ii]*len);
				}
			}
		}
	} while (numread == BUFSIZE*packet_size);
	
	outfle = fopen(argv[argc - 1],"w");
	fwrite(img,sizeof(double),DIM*DIM,outfle);
	fclose(outfle);
	
	free(mus);
	free(buffer);
	free(img);
	fclose(fle);
}
