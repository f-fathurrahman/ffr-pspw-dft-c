#include <string>
#include <vector>
#include <complex.h>

#ifdef _USE_CUBLAS2
  #include "cublas_v2.h"
#endif

using namespace std;

//
// Several shortcut variables
//
int I_ONE = 1;
int I_TWO = 2;
double complex Z_ZERO = 0.0 + I*0.0;
double complex Z_ONE = 1.0 + I*0.0;
double complex MZ_ONE = -1.0 + I*0.0;

//
// Global variables
//
int ISCF; // 0=SCF, 1=ChebySCF

int CHEBYPOL_DEGREE;

int IDIAG; // 0=Davidson 1=LOBPCG
int MAX_DIAG_ITER;
bool USER_MAX_DIAG_ITER=false;
bool DIAG_VERBOSE=false;

int NSP;
int *NA=NULL;
int NAX;
double A1[3], A2[3], A3[3];
double ALAT;
double OMEGA;
double B1[3], B2[3], B3[3];
vector<string> ATOMTYP;
double *ZV=NULL;
double *RGAUSS=NULL;
int NLMAX;
int *L_MAX=NULL;
int *L_LOC=NULL;
double *TAU=NULL;

double ECUT;
int NEMPTY;
int NX;
int NSTATX;
double NEL;
int NR1, NR2, NR3, NNR;
bool USER_FFT_GRID=false;

int NKPT;
double *XK=NULL;
double *WKPT=NULL;

double TPIBA;
double TPIBA2;
double GCUT;
double GCUTW;
int NGWX;
int NGX;
int NG;
double *G=NULL;
double *GG=NULL;
double *GGK=NULL;
int *N1=NULL, *N2=NULL, *N3=NULL;
int *IN1=NULL, *IN2=NULL, *IN3=NULL;
int *NGW=NULL;
int *IGK=NULL;
double *XKG=NULL;
unsigned int *N123=NULL;

bool TMETAL;
double EKT;
double *FOCC=NULL;

int *MMAX=NULL;
double *CLOG=NULL;
double *R=NULL;
double *PSI=NULL;
double *VION=NULL;

double complex *C0=NULL;
double complex *EI1=NULL, *EI2=NULL, *EI3=NULL;
double complex *EIGR=NULL;
double complex *SFAC=NULL;

double *VPS=NULL;
double *RHOPS=NULL;
double *PKG=NULL;
double *PKG_A=NULL;
double *WNL=NULL;

double ESELF;
double ESR;
double ETOT;
double EKIN;
double EHT;
double EPSEU;
double ENL;
double EXC;
double ERHO;

double *EIG=NULL;

double complex* FNL=NULL;

int MAX_SCF_ITER;
double SCF_CONV_CRIT;
double RHO_MIX_FAC;
int IMIX;
double BETAMAX;
double KERKER_Q0;
int MIXDIM;

double CPUMEM=0.0;
double GPUMEM=0.0;

//
// Arrays allocated on GPU
//
bool USE_CUDA=false;
unsigned int *d_N123;
double *d_XKG;
double complex *d_EIGR;
double *d_PKG_A;
double *d_rho;
double complex *d_c_fft;
double complex *d_eigr_pkg;

#ifdef _USE_CUBLAS2
cublasHandle_t CUBLAS_HANDLE;
#endif
