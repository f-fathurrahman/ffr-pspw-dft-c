// eFeFeR (20910015)

#include "common_pspw_cuda.h"

void eig_zheevd(double complex *A, int ldA, double *lambda, int N)
{
  // Local
  int lwork, lrwork, liwork;
  double complex *work=NULL;
  double *rwork=NULL;
  int *iwork=NULL;
  int info;
  char jobz, uplo;

  lwork = N*N + 2*N;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;
  work = (double complex*)malloc(lwork*sizeof(double complex));
  rwork = (double*)malloc(lrwork*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));

  jobz = 'V';
  uplo = 'U';
  zheevd_(&jobz, &uplo, &N, A, &ldA, lambda, 
      work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  if(info != 0) {
    printf("ERROR calling zheevd_ : info = %d\n", info);
    return;
  }

  free(work); work=NULL;
  free(rwork); rwork=NULL;
  free(iwork); iwork=NULL;
}

// Only diagonalize A
void eig_zheevd(double complex *A, int ldA, int N)
{
  // Local
  int lwork, lrwork, liwork;
  double complex *work=NULL;
  double *rwork=NULL;
  int *iwork=NULL;
  int info;
  char jobz, uplo;
  double *lambda=NULL;

  lwork = N*N + 2*N;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;
  work = (double complex*)malloc(lwork*sizeof(double complex));
  rwork = (double*)malloc(lrwork*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));
  lambda = (double*)malloc(N*sizeof(double));

  jobz = 'V';
  uplo = 'U';
  zheevd_(&jobz, &uplo, &N, A, &ldA, lambda, 
      work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  if(info != 0) {
    printf("ERROR calling zheev_ : info = %d\n", info);
    return;
  }

  free(lambda); lambda=NULL;
  free(work); work=NULL;
  free(rwork); rwork=NULL;
  free(iwork); iwork=NULL;
}

