// eFeFeR (20910015), October 2011

#ifndef _COMMON_PSPW_H_
#define _COMMON_PSPW_H_

// Standard headers
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string.h> // Work around for GCC 4.4.0
#include <vector>
#include <math.h>
#include <complex.h>

using namespace std;

#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))
#define IDX3F(i,j,k,DIM1,DIM2) (((k)-1)*(DIM1)*(DIM2) + ((j)-1)*(DIM1) + ((i)-1))
#define IDX4F(i,j,k,l,DIM1,DIM2,DIM3) (((l)-1)*(DIM1)*(DIM2)*(DIM3) + ((k)-1)*(DIM1)*(DIM2) + ((j)-1)*(DIM1) + ((i)-1))
#define IDX5F(i,j,k,l,m,DIM1,DIM2,DIM3,DIM4) (((m)-1)*(DIM1)*(DIM2)*(DIM3)*(DIM4) + ((l)-1)*(DIM1)*(DIM2)*(DIM3) + ((k)-1)*(DIM1)*(DIM2) + ((j)-1)*(DIM1) + ((i)-1))

#define SQUARE(x) ((x)*(x))
#define POW3(x) ((x)*(x)*(x))

//
// Several parameters
//
extern int I_ONE;
extern int I_TWO;
extern double complex Z_ZERO;
extern double complex Z_ONE;
extern double complex MZ_ONE;

// Maximum angular momentum component (d)
#define MAXAM 3
// Maximum radial grid number
#define MAXR 570


//
// Global variables
//
extern int ISCF;
extern int CHEBYPOL_DEGREE;

extern int IDIAG;
extern bool USER_MAX_DIAG_ITER;
extern int MAX_DIAG_ITER;
extern bool DIAG_VERBOSE;

extern int NSP;
extern int *NA;
extern int NAX;
extern double A1[3], A2[3], A3[3];
extern double ALAT;
extern double OMEGA;
extern double B1[3], B2[3], B3[3];
extern vector<string> ATOMTYP;
extern double *ZV;
extern double *RGAUSS;
extern int NLMAX;
extern int *L_MAX;
extern int *L_LOC;
extern double *TAU;

extern double ECUT;
extern int NEMPTY;
extern int NX;
extern int NSTATX;
extern double NEL;
extern int NR1, NR2, NR3, NNR;
extern bool USER_FFT_GRID;

extern int NKPT;
extern double *XK;
extern double *WKPT;

extern double TPIBA;
extern double TPIBA2;
extern double GCUT;
extern double GCUTW;
extern int NGWX;
extern int NGX;
extern int NG;
extern double *G;
extern double *GG;
extern double *GGK;
extern int *N1, *N2, *N3;
extern int *IN1, *IN2, *IN3;
extern int *NGW;
extern int *IGK;
extern double *XKG;
extern unsigned int *N123;

extern bool TMETAL;
extern double EKT;
extern double *FOCC;

extern int *MMAX;
extern double *CLOG;
extern double *R;
extern double *PSI;
extern double *VION;

extern double complex *C0;
extern double complex *EI1, *EI2, *EI3;
extern double complex *EIGR;
extern double complex *SFAC;

extern double *VPS;
extern double *RHOPS;
extern double *PKG;
extern double *PKG_A;
extern double *WNL;

extern double ESELF;
extern double ESR;
extern double ETOT;
extern double EKIN;
extern double EHT;
extern double EPSEU;
extern double ENL;
extern double EXC;
extern double ERHO;

extern double *EIG;

extern double complex *FNL;

extern int MAX_SCF_ITER;
extern double SCF_CONV_CRIT;
extern double RHO_MIX_FAC;
extern int IMIX;
extern double BETAMAX;
extern double KERKER_Q0;
extern int MIXDIM;

extern double CPUMEM;
extern double GPUMEM;

//
// Global arrays allocated on GPU
//
extern unsigned int *d_N123;
extern double *d_XKG;
extern double complex *d_EIGR;
extern double *d_PKG_A;
extern double *d_rho;

//
// PSPW_DFT functions
//
void ReadInputFile();
void Setup();
void CalcRecipLatt();
void GVectorsGen0();
void ReadPseudopotentials();
void AllocateArrays();
void InitOcc();
void GVectorsGen();
void SortGVectors(int N, double *Gvec, int *indexG);
void driver_orderf(int N, double *Gvec, int *indexG);
void GkVectorsGen();
void PrepareSCF();
void PhaseFactor();
void StructureFactor();
void FormFactor();
void NLFormFactor();
void SCFLoop();
void FormFactorAtomic(double complex* c_fft);
void NormalizeRho(double *rho, int nnr, double omega, double nel);
void EwaldEnergy();
void EvalEnergy(double *rhoe);
void xc_PZ(double *rho, double *vxc, double &Exc, int NNR);
void PrintSystemInfo();
void ApplyHam(double complex *vin, double complex *vout, int NGW_ik, int ik, double *rho, int nnr);
void ApplyHam_block(double complex *vin, double complex *vout, int NGW_ik, int ik, double *rho, int nnr,
    int nstates);
void EvalRhoPsi(int ik, double *rhoe, double complex *eigvec, int NGW_ik);
void EvalRhoPsi(int ik, double *rhoe, double complex *eigvec, int ldEigvec, int NGW_ik);
void EvalENL();
void GenPrec(int ik, double *prec, int nbasis, double tau);

void ChebySCFLoop();
void lanczos(int ik, double complex *v0, int k, double *T, int ldT, double complex *f,
    double *rhoe, int nnr);
void chebyfilt(int ik, double complex *X, int degree, double lb, double ub, double *rhoe, int nnr);

void SimpleRhoMix(double beta, double *rhoe, double *rho_old, int nnr, double &d);
void AndersonRhoMix(int iter, double beta, double *rhoe, double *rho_old,
    double *f, double *g, int NNR, double &d);
void BroydenRhoMix(int iter, double alpha, double w0, double *rho, double *rho_old,
    double *f, double *df, double *u, double *a, int mixdim, int nnr, double &d);
void PulayRhoMix(int iter, double *rhoe, double *rho_old, double *f, int mixdim, int nnr, double &d);
void KerkerRhoMix(double beta, double q0, double *rhoe, double *rho_old, int nnr, double &d);
void KerkerPrec(double *diff_rhoe, int idx1, int nnr, double beta, double q0);

void lobpcg_cpu(int ik, double *lambda, double complex *X, int nbasis, int nstates,
    double *rhoe, int nnr, double *prec);
void Diagon_lobpcg_cpu(int ik, double *eigval, double complex *eigvec, int NGW_ik,
    double *rhoe, int nnr, double tau);

void davidson_cpu(int ik, double *evals, double complex *X, int nbasis, int nstates,
    double *rhoe, int nnr, double *prec);
void Diagon_davidson_cpu(int ik, double *eigval, double complex *eigvec, int NGW_ik,
    double *rhoe, int nnr, double tau);

void lobpcg_cuda(int ik, double *lambda, double complex *X, int nbasis, int nstates,
    double *prec);
void Diagon_davidson_cuda(int ik, double *eigval, double complex *eigvec, int NGW_ik, double tau);
void Diagon_lobpcg_cuda(int ik, double *eigval, double complex *eigvec, int NGW_ik, double tau);

void test_ApplyHam(double *rhoe);
void test_ApplyHam_cuda(double *rhoe);

void cufft3d_z2z(double complex *d_data,int DIM1, int DIM2, int DIM3, bool inverse);
void eig_zheevd(double complex *A, int ldA, double *lambda, int N);
void eig_zhegv(double complex *A, int ldA, double complex *B, int ldB,
    double complex *evec, int ldE, int N);
void eig_zhegv_eval(double complex *A, int ldA, double complex *B, int ldB,
    double *eval, double complex *evec, int ldE, int N);
void ortho_qr(double complex *X, int nbasis, int nstates);

//
// MAGMA
//
typedef int magma_int_t;
extern "C"
magma_int_t magma_zheevd_gpu( char jobz, char uplo,
                              magma_int_t n,
                              double complex *da, magma_int_t ldda,
                              double *w,
                              double complex *wa,  magma_int_t ldwa,
                              double complex *work, magma_int_t lwork,
                              double *rwork, magma_int_t lrwork,
                              magma_int_t *iwork, magma_int_t liwork,
                              magma_int_t *info);
extern "C" magma_int_t magma_zgeqrf_gpu( magma_int_t m, magma_int_t n,
                              double complex *dA,  magma_int_t ldda,
                              double complex *tau, double complex *dT,
                              magma_int_t *info);
extern "C" magma_int_t magma_zungqr_gpu( magma_int_t m, magma_int_t n, magma_int_t k,
                              double complex *da, magma_int_t ldda,
                              double complex *tau, double complex *dwork,
                              magma_int_t nb, magma_int_t *info );
extern "C" int magma_get_zhetrd_nb(int m);
extern "C" int magma_get_zgeqrf_nb(int m);
extern "C" magma_int_t magma_zpotrf_gpu( char uplo,  magma_int_t n,
                              double complex *dA, magma_int_t ldda, magma_int_t *info);


void SetupCUDA();
void ShutdownCUDA();
int ApplyHam_block_cuda(double complex *d_vin, double complex *d_vout, int NGW_ik, int ik, int nstates);
int SetupArrays_ApplyHam_cuda(unsigned int *N123, double *XKG, double complex *EIGR, double *PKG_A);
int SetupRho_ApplyHam_cuda();
int CopyRho_ApplyHam_cuda(double *rho);
void FreeArrays_ApplyHam_cuda();
void FreeRho_ApplyHam_cuda();
void eig_zhegv_eval_magma(double complex *A, int ldA, double complex *B, int ldB,
    double *eval, double complex *evec, int ldE, int N);
void ortho_qr_magma(double complex *dX, int nbasis, int nstates);
void eig_zheevd_magma(double complex *dA, int lddA, double *lambda, int N);
void eig_zhegv_magma(double complex *dA, int lddA, double complex *dB, int lddB,
    double complex *devec, int lddE, int N);


//
// BLAS Level 1
//
extern "C" double ddot_(int *n, double *x, int *incx, double *y, int *incy);
extern "C" double dnrm2_(int *n, double *x, int *incx);
extern "C" void zdotc_(double complex *pres, int *n, double complex *x, int *incx, double complex *y, int *incy);
extern "C" void zdotu_(double complex *pres, int *n, double complex *x, int *incx, double complex *y, int *incy);
extern "C" void zscal_(int *n, double complex *a, double complex *x, int *incx);
extern "C" void zdscal_(int *n, double *a, double complex *x, int *incx);
extern "C" void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
extern "C" void zcopy_(int *n, double complex *x, int *incx, double complex *y, int *incy);
extern "C" void zaxpy_(int *n, double complex *a, double complex *x, int *incx, double complex *y, int *incy);


//
// BLAS Level 2
//
extern "C" void zgemv_(char *trans, int *m, int *n,
    double complex *alpha, double complex *a, int *lda,
    double complex *x, int *incx,
    double complex *beta, double complex *y, int *incy);
extern "C" void dgemv_(char *trans, int *m, int *n,
    double *alpha, double *a, int *lda,
    double *x, int *incx,
    double *beta, double *y, int *incy);


//
// BLAS Level 3
//
extern "C" void zgemm_(char *transa, char *transb,
    int *m, int *n, int *k, double complex *alpha, double complex *A, int *lda,
    double complex *b, int *ldb, double complex *beta, double complex *C, int *ldc);
extern "C" void ztrsm_(char *side, char *uplo, char *transa, char *diag,
    int *m, int *n, double complex *alpha, double complex *A, int *ldA, double complex *B, int *ldB);

//
// LAPACK
//
extern "C" void zheevd_(char *jobz, char *uplo, int *n, double complex *A, int *ldA, double *w,
    double complex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);
extern "C" void zgeqrf_(int *m, int *n, double complex *X, int *lda, double complex *tau,
    double complex *work, int *lwork,int *info);
extern "C" void zungqr_(int *m, int *n, int *k, double complex *A, int *ldA, double complex *tau,
    double complex *work, int *lwork, int *info);
extern "C" void zpotrf_(char *uplo, int *n, double complex *A, int *ldA, int *info);
extern "C" void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern "C" void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
extern "C" void dsysv_(char *uplo, int *n, int *nrhs, double *a, int *lda,
    int *ipiv, double *b, int *ldb, double *work, int *lwork, int *info);
extern "C"void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda,
    double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr,
    double *work, int *lwork, int *info);
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
    double *work, int *lwork, int *info);


//
// Utility functions
//
void PrintVector(double complex *V, int nelem);
void PrintVector(double *V, int nelem);
void PrintMatrix(double complex *A, int NROW, int NCOL);
void PrintMatrix(double *A, int NROW, int NCOL);
void Print3DArray(double *A, int dim1, int dim2, int dim3, string filename);

double determinant3x3(double A1[3], double A2[3], double A3[3]);

double cpu_time();
char *timestring();
void timestamp();

void sort1(int N, double *arrin, int *index);
void SumArray(int *a, int N, int &sum);
void SumArray(unsigned int *a, unsigned int N, unsigned int &sum);
void SumArray(double *a, int N, double &sum);
void SumArray(double complex *a, int N, double complex &sum);
void matinv3x3(double *HM, double *HI);
double simpson(int N, double *F, double h);
double dsum(int N, double *x, int incx);
void fft_fftw3(double complex *data, int NR1, int NR2, int NR3, bool inverse);

#endif 

