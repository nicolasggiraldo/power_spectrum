////////////////////////////////////////////////////
// HEADER WITH ALL DATA SRUCTURES FOR THE PROGRAM //
////////////////////////////////////////////////////



/////////////////////////////
// PREPROCESSOR DIRECTIVES //
/////////////////////////////

/* Index preprocessor for the C-Order */
#define INDEX(i,j,k)     (1L*k) + GV.NGRID * ( (1L*j) + GV.NGRID * (1L*i) )
#define COMPLEXMAG(A,i)  ( (A[i][0] * A[i][0]) + (A[i][1] * A[i][1]) )
#define VECTORMAG(x,y,z) sqrt( (x * x) + (y * y) + (z * z) )

#define X      0
#define Y      1
#define Z      2



////////////////////////////////////////////
// GLOBAL VARIABLES FOR THE FFTW ROUTINES //
////////////////////////////////////////////

fftw_complex *denConX = NULL; // Density contrast in X-space
fftw_complex *denConK = NULL; // Density contrast in K-space
fftw_plan forwardPlan;



///////////////////////////////
// STRUCTURES OF THE PROGRAM //
///////////////////////////////

/* Global variales */
struct globalVariables{
  int      NGRID;          // Number of cell in each axis.
  long int NGRID3;         // Total number of cells (NGRID3 = NGRID^3)
  int      GADGET_VERSION; // GADGET version of the snapshot
  double   L;              // Lenght of the simulation in Mpc
  double   SIM_VOL;        // Volume of the simulation
  long int NP_TOT;         // Total number of particles in the simulation
  double   TOTAL_MASS;     // Total mass of all particles in the simulation
  double   RHO_MEAN;       // Mean density of ALL the simulation
  double   VOL_CELL;       // Volume of each cell
  double   H;              // Size of the cell
  double   S_KF;           // Constant to define the value of DELTA_K = S_KF * KF
  double   DELTA_K;        // Delta k for the binning
  double   KF;             // Fundamental frequency kF
  double   SHOT_NOISE;     // shotNoise of the Power Spectrum
  double   KN;             // Nyquist frequency
  char     *FILE_NAME;     // Path of the GADGET binary
  char     *OUTPUT;        // Name of the outputfile
  char     SCHEME[4];      // Scheme used for grid assignation

  /* COSMOLOGICAL PARAMETERS OF THE SIMULATION */
  double OMEGA_M0;         //Omega matter at present time
  double OMEGA_L0;         //Omega Lambda at present time
  double ZRS;              //Redshift of the simulation
  double HUBBLEPARAM;      //Hubble parameter of the simulation
  
}GV;
