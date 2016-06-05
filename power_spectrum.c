#include        <stdio.h>
#include       <stdlib.h>
#include         <math.h>
#include    <libconfig.h>
#include        <fftw3.h>
#include       <string.h>
#include <gsl/gsl_sort.h>

#include "constansts_and_structures.h"
#include                 "functions.h"

/*
  FFTW_MEASURE: find the optimal plan by actually computing several FFTs 
  FFTW_ESTIMATE: do not run any FFT and provide a "reasonable" plan
  FFTW_OUT_OF_PLACE: a plan assumes that the in and out are distinct 
  FFTW_IN_PLACE: a plan assumes that the in and out are same
*/

int main(int argc, char *argv[]){
  int i, j, k, l; // Counter in the X, Y, and Z axis.
  long int id_cell;
  int n[3]; // Number of grids in each dimension.
  double kx, ky, kz; // Wavenumber componentst
  int Nbins;
  double *denk2 = NULL;
  int *kn = NULL;
  double Pk, PkError, kdist;
  size_t *order=NULL;
  double *kMag=NULL; /* Array with the magnitude of the 
			wavenumber vector */
  double *kpos; /* Array of positions values according 
		   to FFTW k-position convention */ 
  double (*W2_k)(double) = NULL; /* Memory addres to the function for 
				    correction of the Fourier mode for
				    the effect of mass assignment scheme
				    according to Montesano et al. 2012 */
  FILE *fout = NULL;


  //////////////////////////////
  //* READING PARAMETER FILE *//
  //////////////////////////////
  if(argc<2){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("%s: You must specify the name of the parameter file\n",argv[0]);
    printf("For example: %s pathOfFile/parameterFile.txt\n",argv[0]);
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  // Reading parameter file and verifying there is no error.
  switch( read_parameters(argv[1]) ){
  case -1 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad path to the parameter file.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  case -2 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad settings in the parameter file.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ////////////////////////////////
  //* READING CELL BINARY FILE *//
  ////////////////////////////////
  switch( readBinaryFile() ){
  case -1 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: The parameter file could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  case -2 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: FFTW arrays could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ///////////////////////////////////////////////////////
  //* DEFINING CORRECTION TERM FOR THE POWER SPECTRUM *//
  ///////////////////////////////////////////////////////
  if(      strcmp(GV.SCHEME, "NGP") == 0 ){
    W2_k = Sum_W2_NGP; // NGP
  }
  else if( strcmp(GV.SCHEME, "CIC") == 0 ){
    W2_k = Sum_W2_CIC; // CIC
  }
  else if( strcmp(GV.SCHEME, "TSC") == 0 ){
    W2_k = Sum_W2_TSC; // TSC
  }
  else if( strcmp(GV.SCHEME, "D20") == 0 ){
    W2_k = Sum_W2_D20; // D20
  }
  
  
  
  ////////////////////////////////
  //* EXECUTING FFTW3 ROUTINES *//
  ////////////////////////////////
  
  /* Number of grids in each axis */
  n[X] = n[Y] = n[Z] = GV.NGRID;
  
  /* FFTW plan for a 3D Fourier transform */
  forwardPlan = fftw_plan_dft(3, n, 
			      denConX, denConK,
			      FFTW_FORWARD, FFTW_ESTIMATE);

  /* Do forward FFT */
  fftw_execute(forwardPlan);
  printf("\n-----------------------------------------------\n");
  printf("Fourier Transform succes.\n");
  
  /* Magnitud of the k vector */
  kMag = (double *) calloc(GV.NGRID3, sizeof(double));
  if(kMag == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Wavenumber magnitude array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  /* Position array for storing in the densityContrast */
  kpos = (double *) calloc(GV.NGRID, sizeof(double));
  if(kpos == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: kpos array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  /* Setting index and positions according to FFTW convention  */
  /* REMEMBER kF is (2.0*PI)/L */
  for( i=0; i<GV.NGRID; i++ )
    kpos[i] = (i<GV.NGRID/2) ? GV.KF*i : GV.KF*(i-GV.NGRID);
  
  
  /* REMEMBER kF is (2.0*PI)/L; */
  for(i=0; i<GV.NGRID; i++){
    // Momentum coordinate in the X-Axis
    //kx = (i<=GV.NGRID/2) ? GV.KF*i : GV.KF*(i-GV.NGRID);
    kx = kpos[i];
    
    for(j=0; j<GV.NGRID; j++){
      // Momentum coordinate in the Y-Axis
      //ky = (j<=GV.NGRID/2) ? GV.KF*j : GV.KF*(j-GV.NGRID);
      ky = kpos[j];
      
      for(k=0; k<GV.NGRID; k++){
	// Momentum coordinate in the Z-Axis
	//kz = (k<=GV.NGRID/2) ? GV.KF*k : GV.KF*(k-GV.NGRID);
	kz = kpos[k];

	// Distance from kX=0, kY=0, kZ=0.
	kMag[INDEX(i,j,k)] = VECTORMAG(kx,ky,kz);

      }// for k
    }// for j
  }// for i 

  printf("\n-----------------------------------------------\n");
  printf("Fourier space values assigned.\n");
  
  
  
  ////////////////////////////////////////////////////////
  //* SORTING THE ARRAY ACCORDING TO THE KMAG DISTANCE *//
  ////////////////////////////////////////////////////////
  
  /* Array whith the order in encreasing form */
  order = malloc(GV.NGRID3 * sizeof(size_t));
  if(order == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Index array for sorting could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  /* GSL routine for index sorting */
  gsl_sort_index(order, kMag, 1, GV.NGRID3);
  
  printf("\n-----------------------------------------------\n");
  printf("Sort succes.\n");
  
  
  
  //////////////////////////////////////////////////
  //* BINNING IN ORDER TO GET THE POWER SPECTRUM *//
  //////////////////////////////////////////////////

  /* Binning until the last value of Kmag */
  Nbins = (int) floor( kMag[ order[GV.NGRID3-1L] ] / GV.DELTA_K );

  printf("\n-----------------------------------------------\n");
  printf("KMag_max = %lf\n", kMag[ order[GV.NGRID3-1L] ]);
  printf("Nbins    = %d\n",  Nbins);
  
  denk2 = (double *) calloc(Nbins, sizeof(double));
  if(denk2 == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Array of contrast density to square for binning could not be allocated.\n");
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  kn = (int *) calloc(Nbins, sizeof(int));
  if(kn == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Array of bin positions could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  /* Lineal binning  */
  id_cell = 0L;
  for(l=0; l<Nbins; l++){
    while( kMag[ order[id_cell] ] <= (GV.DELTA_K * (l+1)) ){
      /* Aliasing correction */
      k = order[id_cell]%GV.NGRID;
      j = ((order[id_cell]-k)/GV.NGRID)%GV.NGRID;
      i = (((order[id_cell]-k)/GV.NGRID)-j)/GV.NGRID;

      denk2[l] += COMPLEXMAG(denConK,order[id_cell])/( W2_k(kpos[i]) *
						       W2_k(kpos[j]) *
						       W2_k(kpos[k]) );
      kn[l] += 1;
      id_cell++;
    }
  }//for l

  
  /* Saving data in the outfile */
  printf("\n-----------------------------------------------\n");
  printf("Saving data in %s\n", GV.OUTPUT);

  /* Opening file */
  fout = fopen(GV.OUTPUT, "w");
  if(fout == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Outfile could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  /* Writing header */
  fprintf(fout,"# NGRID          = %d\n",  GV.NGRID);
  fprintf(fout,"# GADGET VERSION = %d\n",  GV.GADGET_VERSION);
  fprintf(fout,"# L              = %lf\n", GV.L);
  fprintf(fout,"# SIM VOL        = %lf\n", GV.SIM_VOL);
  fprintf(fout,"# NP TOT         = %ld\n", GV.NP_TOT);
  fprintf(fout,"# TOTAL MASS     = %lf\n", GV.TOTAL_MASS);
  fprintf(fout,"# RHO MEAN       = %lf\n", GV.RHO_MEAN);
  fprintf(fout,"# VOL_CELL       = %lf\n", GV.VOL_CELL);
  fprintf(fout,"# H              = %lf\n", GV.H);
  fprintf(fout,"# DELTA k        = %lf\n", GV.DELTA_K);
  fprintf(fout,"# kF             = %lf\n", GV.KF);
  fprintf(fout,"# kN             = %lf\n", GV.KN);
  fprintf(fout,"# Shot Noise     = %lf\n", GV.SHOT_NOISE);
  fprintf(fout,"# SCHEME         = %s\n",  GV.SCHEME);
  fprintf(fout,"# OMEGA_M0       = %lf\n", GV.OMEGA_M0);
  fprintf(fout,"# OMEGA_L0       = %lf\n", GV.OMEGA_L0);
  fprintf(fout,"# ZRS            = %lf\n", GV.ZRS);
  fprintf(fout,"# HUBBLEPARAM    = %lf\n", GV.HUBBLEPARAM);
  fprintf(fout,"\n");

  /* Saving power spectrum values */
  fprintf(fout,"#%19s %20s %20s\n","k", "Pk(k)", "Pk_Error(k)");

  for(j=0; j<Nbins; j++){
    
    kdist = GV.DELTA_K * (j+0.5);
    
    if( kn[j] > 0){
      
      /* Mean denk2 at k */
      Pk  = denk2[j]/(1.0*kn[j]);
      
      /* Normaliztion of the power spectrum */
      Pk *= (GV.H * GV.H * GV.H) / GV.NGRID3;
      
      /* Correction of the Fourier mode for the effect 
	 of the Mass Assignment Scheme according to 
	 Montesano et al. 2010 */
      //Pk /= C1(kdist);
      
      /* substracting shot noise term */
      Pk -= GV.SHOT_NOISE;

      /* Error of the estimation */
      PkError = (Pk*(GV.DELTA_K/kdist))/sqrt(2.0*M_PI);
      
    }else{

      /* If counts are zero, set Pk to zero */
      Pk = 0.0;
      PkError = 0.0;
      
    }

    fprintf(fout,"%20lf %20e %20e\n",kdist, Pk, PkError);
    
  }
  
  
  
  ///////////////////
  //* FREE MEMORY *//
  ///////////////////
  free(kMag);  
  fftw_free(denConX);
  fftw_free(denConK);
  fftw_destroy_plan(forwardPlan);
  
  return 0;
}
