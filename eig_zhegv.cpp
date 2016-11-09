// eFeFeR (20910015), November 2011

#include "common_pspw_cuda.h"

void eig_zhegv(double complex *A, int ldA, double complex *B, int ldB,
    double complex *evec, int ldE, int N)
{
  const double SMALL = 12.220446049250313e-16;
  int I_ONE=1;
  int lwork, lrwork, liwork, info, i;
  char jobz, uplo, transa, transb;
  double complex *work=NULL;
  double *rwork=NULL;
  int *iwork=NULL;
  double *eval=NULL;
  double scale;
  int nn;
  double complex Z_ONE = 1.0 + I*0.0;
  double complex Z_ZERO = 0.0 + I*0.0;

  lwork = N*N + 2*N;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;
  work = (double complex*)malloc(lwork*sizeof(double complex));
  rwork = (double*)malloc(lrwork*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));
  eval = (double*)malloc(N*sizeof(double));

  jobz = 'V';
  uplo = 'U';
  zheevd_(&jobz, &uplo, &N, B, &ldB, eval,
      work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  if(info != 0) {
    printf("ERROR calling zheevd in eig_zhegv : info = %d\n", info);
    return;
  }

  nn = 0;
  for(i=1; i<=N; i++) {
    if(eval[i-1] > SMALL) {
      nn = nn + 1;
      scale = 1.0/sqrt(eval[i-1]);
      zdscal_(&N, &scale, &B[IDX2F(1,i,ldB)], &I_ONE);
    } else {
      printf("Small eigenvalue detected = %f\n", eval[i-1]);
    }
  }
  if(nn < N) {
    printf("Warning: Number of linearly independent vectors = %d\n", nn);
    printf("         while size of the problem: %d\n", N);
  }

  transa = 'N';
  transb = 'N';
  zgemm_(&transa, &transb, &N, &N, &N, &Z_ONE, A, &ldA, B, &ldB, &Z_ZERO, evec, &ldE);
  //
  transa = 'C';
  transb = 'N';
  zgemm_(&transa, &transb, &N, &N, &N, &Z_ONE, B, &ldB, evec, &ldE, &Z_ZERO, A, &ldA);

  // Diagonalize transformed A
  jobz = 'V';
  uplo = 'U';
  zheevd_(&jobz, &uplo, &N, A, &ldA, eval,
      work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  if(info != 0) {
    printf("ERROR calling zheevd in eig_zhegv : info = %d\n", info);
    return;
  }
  // Backtransform eigenvectors
  transa = 'N';
  transb = 'N';
  zgemm_(&transa, &transb, &N, &N, &N, &Z_ONE, B, &ldB, A, &ldA, &Z_ZERO, evec, &ldE);

  //for(i=1; i<=N; i++) {
  //  printf("Eigenvalues = %f\n", eval[i-1]);
  //}

  free(eval); eval=NULL;
  free(work); work=NULL;
  free(rwork); rwork=NULL;
  free(iwork); iwork=NULL;
}


void eig_zhegv_eval(double complex *A, int ldA, double complex *B, int ldB,
    double *eval, double complex *evec, int ldE, int N)
{
  const double SMALL = 12.220446049250313e-16;
  int I_ONE=1;
  int lwork, lrwork, liwork, info, i;
  char jobz, uplo, transa, transb;
  double complex *work=NULL;
  double *rwork=NULL;
  int *iwork=NULL;
  double scale;
  int nn;
  double complex Z_ONE = 1.0 + I*0.0;
  double complex Z_ZERO = 0.0 + I*0.0;

  lwork = N*N + 2*N;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;
  work = (double complex*)malloc(lwork*sizeof(double complex));
  rwork = (double*)malloc(lrwork*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));

  jobz = 'V';
  uplo = 'U';
  zheevd_(&jobz, &uplo, &N, B, &ldB, eval,
      work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  if(info != 0) {
    printf("ERROR calling zheevd in eig_zhegv : info = %d\n", info);
    return;
  }

  nn = 0;
  for(i=1; i<=N; i++) {
    if(eval[i-1] > SMALL) {
      nn = nn + 1;
      scale = 1.0/sqrt(eval[i-1]);
      zdscal_(&N, &scale, &B[IDX2F(1,i,ldB)], &I_ONE);
    } else {
      printf("Small eigenvalue detected = %f\n", eval[i-1]);
    }
  }
  if(nn < N) {
    printf("Warning: Number of linearly independent vectors = %d\n", nn);
    printf("         while size of the problem: %d\n", N);
  }

  transa = 'N';
  transb = 'N';
  zgemm_(&transa, &transb, &N, &N, &N, &Z_ONE, A, &ldA, B, &ldB, &Z_ZERO, evec, &ldE);
  //
  transa = 'C';
  transb = 'N';
  zgemm_(&transa, &transb, &N, &N, &N, &Z_ONE, B, &ldB, evec, &ldE, &Z_ZERO, A, &ldA);

  // Diagonalize transformed A
  jobz = 'V';
  uplo = 'U';
  zheevd_(&jobz, &uplo, &N, A, &ldA, eval,
      work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  if(info != 0) {
    printf("ERROR calling zheevd in eig_zhegv : info = %d\n", info);
    return;
  }
  // Backtransform eigenvectors
  transa = 'N';
  transb = 'N';
  zgemm_(&transa, &transb, &N, &N, &N, &Z_ONE, B, &ldB, A, &ldA, &Z_ZERO, evec, &ldE);

  //for(i=1; i<=N; i++) {
  //  printf("Eigenvalues = %f\n", eval[i-1]);
  //}

  free(work); work=NULL;
  free(rwork); rwork=NULL;
  free(iwork); iwork=NULL;
}
