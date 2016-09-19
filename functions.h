
double h[] = { 2.667005790055555358661744877130858277192498290851289932779975e-02,
	       1.881768000776914890208929736790939942702546758640393484348595e-01,
	       5.272011889317255864817448279595081924981402680840223445318549e-01,
	       6.884590394536035657418717825492358539771364042407339537279681e-01,
	       2.811723436605774607487269984455892876243888859026150413831543e-01,
	       -2.498464243273153794161018979207791000564669737132073715013121e-01,
	       -1.959462743773770435042992543190981318766776476382778474396781e-01,
	       1.273693403357932600826772332014009770786177480422245995563097e-01,
	       9.305736460357235116035228983545273226942917998946925868063974e-02,
	       -7.139414716639708714533609307605064767292611983702150917523756e-02,
	       -2.945753682187581285828323760141839199388200516064948779769654e-02,
	       3.321267405934100173976365318215912897978337413267096043323351e-02,
	       3.606553566956169655423291417133403299517350518618994762730612e-03,
	       -1.073317548333057504431811410651364448111548781143923213370333e-02,
	       1.395351747052901165789318447957707567660542855688552426721117e-03,
	       1.992405295185056117158742242640643211762555365514105280067936e-03,
	       -6.858566949597116265613709819265714196625043336786920516211903e-04,
	       -1.164668551292854509514809710258991891527461854347597362819235e-04,
	       9.358867032006959133405013034222854399688456215297276443521873e-05,
	       -1.326420289452124481243667531226683305749240960605829756400674e-05 };

/*
 * Function:  read_parameters
 * --------------------
 * Reads the parameter file in which are the main parameters 
 * necessary to run the code.
 *
 * The information loaded are:
 * FILE_NAME:      File name path of the GADGET binary file.
 * OUTPUT:         Path of the output file.
 * DELTA_K:        Width of space sampled for the calculation of
 *                 the power spectrum. The value is given in terms
 *                 of the fundamental frequency kF. The value
 *                 should be bigger or equal than 1.
 * 
 * The parameter file is read with the help of the library libconfig,
 * for more information of the library use go to the manual:
 * http://www.hyperrealm.com/libconfig/libconfig_manual.html
 *
 *
 *  param_file_name: String with the name of the parameter file.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the parameter file.
 *           -2 --> There is an error whith the settings of the 
 *                  parameter file.
 */
int read_parameters(char param_file_name[]){
  
  FILE *cfg=NULL;  // Stream to the parameter (config) file
  int len = 200;   // Len of the read parameter
  char *buf =NULL; // buf variables to be used to read strings variables
  char *buf1=NULL;
  char *buf2=NULL; 
  char *dumb=NULL;
  
  
  if( (cfg=fopen(param_file_name,"r"))==NULL ){
    printf("%s not found.\n", param_file_name);
    // Value -1 means there is an error loading the param file
    return -1;
  }
  
  buf  = (char *) malloc( len*sizeof(char) );
  buf1 = (char *) malloc( len*sizeof(char) );
  buf2 = (char *) malloc( len*sizeof(char) );

  /* Reading FILE_NAME parameter */  
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    printf("No 'FILE_NAME' setting in configuration file.\n");
    return -2;
  }
  else{
    GV.FILE_NAME = strdup(buf2);
    printf("Reading from File: %s\n", GV.FILE_NAME);
  }

  /* Reading OUTPUT parameter */
  //dumb=fgets(buf,len,cfg);
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    printf("No 'OUTPUT' setting in configuration file.\n");
    return -2;
  }
  else{
    GV.OUTPUT = strdup(buf2);
    printf("Output File: %s\n", GV.OUTPUT);
  }
  
  /* Reading S_KF parameter */
  //dumb=fgets(buf,len,cfg);
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    printf("No 'S_KF' setting in configuration file.\n");
    return -2;
  }
  else{
    GV.S_KF=atof(buf2);
    if(GV.S_KF >= 1.0){
      printf("Binning width in terms of the fundamental frequency kF: %lf\n", GV.S_KF);
    }
    else{
      printf("Invalid 'S_KF' setting in configuration file.\n");
      return -2;
    }
  }
  
  if(dumb==NULL){}
  
  fclose(cfg);
  free(buf);
  free(buf1);
  free(buf2);
    
  return 0;
}



/*
 * Function:  readBinaryFile
 -------------------- 
 * Reads a binary file with the information of the cell and store 
 * the information in the data structure variable *part* it also 
 * returns the total number of particles.
 *
 *  There are no arguments in the routiene.
 *                                                                               
 *  returns: Integer value.
 *            0 --> There is no error.
 *           -1 --> There is an error loading file
 *           -2 --> Structure cell could not be allocated.
 */
int readBinaryFile(){
  FILE *fdata;
  long int id_cell;
  double density_Contrast;
  size_t err;

  printf("\n-----------------------------------------------\n");
  printf("Reading file:   %s\n", GV.FILE_NAME);
  
  fdata = fopen(GV.FILE_NAME,"rb");
  if(fdata == NULL){
    printf("File %s cannot be open\n", GV.FILE_NAME);
    return -1;
  }


  /* Getting cosmological parameters of the simulation */
  err = fread(&GV.OMEGA_M0,    sizeof(double), 1, fdata);
  err = fread(&GV.OMEGA_L0,    sizeof(double), 1, fdata);
  err = fread(&GV.ZRS,         sizeof(double), 1, fdata);
  err = fread(&GV.HUBBLEPARAM, sizeof(double), 1, fdata);

  /* Getting simulation parameters */
  err = fread(&GV.NGRID,          sizeof(int),      1, fdata);
  err = fread(&GV.GADGET_VERSION, sizeof(int),      1, fdata);
  err = fread(&GV.L,              sizeof(double),   1, fdata);
  err = fread(&GV.NP_TOT,         sizeof(long int), 1, fdata);
  err = fread(&GV.TOTAL_MASS,     sizeof(double),   1, fdata);
  err = fread(&GV.RHO_MEAN,       sizeof(double),   1, fdata);
  err = fread(&GV.VOL_CELL,       sizeof(double),   1, fdata);
  err = fread(&GV.H,              sizeof(double),   1, fdata);
  err = fread(&(GV.SCHEME[0]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[1]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[2]),    sizeof(char),     1, fdata);
  GV.SCHEME[3] = '\0';

  GV.NGRID3     = (1L*GV.NGRID) * (1L*GV.NGRID) * (1L*GV.NGRID);
  GV.SIM_VOL    = GV.L * GV.L * GV.L;
  GV.KF         = (2.0*M_PI) / GV.L;
  GV.DELTA_K    = GV.S_KF * GV.KF;
  GV.SHOT_NOISE = GV.VOL_CELL / GV.NP_TOT;
  GV.KN         = M_PI / GV.H;

  printf("\n-----------------------------------------------\n");
  printf("The original snapshot has a total of %ld particles\n", GV.NP_TOT);
  printf("----------------------------------------\n");
  printf(" * Redshift...     %16.8lf\n", GV.ZRS);
  printf(" * Omega0...       %16.8lf\n", GV.OMEGA_M0);
  printf(" * OmageLa...      %16.8lf\n", GV.OMEGA_L0);
  printf(" * Hubbleparam...  %16.8lf\n", GV.HUBBLEPARAM);
  printf("----------------------------------------\n");
  printf(" * Boxsize...      %16.8lf\n", GV.L);
  printf(" * Ngrid...        %16d\n",    GV.NGRID);
  printf(" * SimMass...      %16.8e\n", GV.TOTAL_MASS);
  printf(" * Scheme...       %16s\n",    GV.SCHEME);
  printf("----------------------------------------\n");
  printf(" * kF...           %16.8lf\n", GV.KF);
  printf(" * kN...           %16.8lf\n", GV.KN);
  printf(" * DELTA_k...      %16.8lf\n", GV.DELTA_K);
  printf(" * PSshotNoise...  %16.8e\n", GV.SHOT_NOISE);

  /* Memory allocation for the input and output FFTW arrays,
     for the space, initialize input array to (1.,0.) */
  denConX = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  denConK = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  if(denConX == NULL || denConK == NULL){
    printf("FFTW arrays could not be allocated\n");
    return -2;
  }//if 


  /****************************/
  /* Getting density contrast */
  /****************************/
  for(id_cell=0L; id_cell<GV.NGRID3; id_cell++){
    err = fread(&density_Contrast, sizeof(double), 1, fdata);
    denConX[id_cell][0] = density_Contrast;
    denConX[id_cell][1] = 0.0;
  }//for id_cell

  if(err){};

  return 0;
}

double pow2(double x){
  return x*x;
}

double pow3(double x){
  return x*x*x;
}

double pow4(double x){
  return x*x*x*x;
}

//Sum square window function
double Sum_W2_NGP(double k){
  return 1.0;
}

double Sum_W2_CIC(double k){
  double sine = sin( (M_PI*0.5*k)/GV.KN );
  return 1.0 - 0.666666666666667*pow2( sine );
}

double Sum_W2_TSC(double k){
  double sine = sin( (M_PI*0.5*k)/GV.KN );
  return 1.0 - pow2( sine ) + 0.133333333333333*pow4( sine );
}

void Hfun(double w, double *real, double *imag)
{
  int k;
  
  *real=0.0;
  *imag=0.0;
  
  for(k=0; k<20; k++)
    {
      *real +=  h[k]*cos(w*(k*1.0));
      *imag += -h[k]*sin(w*(k*1.0));
    }
  
  *real *= M_SQRT1_2;
  *imag *= M_SQRT1_2;
}

double Sum_W2_D20(double k)
{
  double arg = (M_PI*k)/GV.KN;
  double prod0[] = {1.0, 0.0};
  double Heval[2], prod1[2];

  while(1)
    {
      arg *= 0.5;
      Hfun( arg, &(Heval[0]), &(Heval[1]) );
      prod1[0] = prod0[0]*Heval[0]-prod0[1]*Heval[1];
      prod1[1] = prod0[0]*Heval[1]+prod0[1]*Heval[0];

      if( fabs(prod0[0]-prod1[0])<1E-10 && fabs(prod0[1]-prod1[1])<1E-10 )
	break;

      prod0[0] = prod1[0];
      prod0[1] = prod1[1];
    }

  return ( (prod1[0]*prod1[0]) + (prod1[1]*prod1[1]) );
}


//double Sum_W2_D20(double k){
//return 1.0;
//}
