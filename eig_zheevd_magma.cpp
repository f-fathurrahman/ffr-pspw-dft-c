// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void eig_zheevd_magma(double complex *dA, int lddA, double *lambda, int N)
{
  // Local
  int lwork, lrwork, liwork;
  int nb;
  double complex *work=NULL, *hA=NULL;
  double *rwork=NULL;
  int *iwork=NULL;
  int info, ldhA;
  char jobz, uplo;

  nb = magma_get_zhetrd_nb(N);
  lwork = N*N + 2*N*nb;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;
  hA = (double complex*)malloc(N*N*sizeof(double complex));
  work = (double complex*)malloc(lwork*sizeof(double complex));
  rwork = (double*)malloc(lrwork*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));
  ldhA = lddA;

  jobz = 'V';
  uplo = 'U';
  magma_zheevd_gpu(jobz, uplo, N, dA, lddA, lambda, hA, ldhA,
      work, lwork, rwork, lrwork, iwork, liwork, &info);
  if(info != 0) {
    printf("ERROR calling magma_zheevd_gpu : info = %d\n", info);
    return;
  }

  free(hA); hA=NULL;
  free(work); work=NULL;
  free(rwork); rwork=NULL;
  free(iwork); iwork=NULL;
}

