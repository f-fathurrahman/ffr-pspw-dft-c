// eFeFeR (20910015), December 2011
//
// CUBLAS should have been initialized before elsewhere

#include "common_pspw.h"

void eig_zhegv_magma(double complex *dA, int lddA, double complex *dB, int lddB, 
    double complex *devec, int lddE, int N)
{
  const double SMALL = 12.220446049250313e-16;
  int I_ONE=1;
  int lwork, lrwork, liwork, info, i, nb, ldhA;
  char jobz, uplo, transa, transb;
  double complex *work=NULL, *hA=NULL;
  double *rwork=NULL;
  int *iwork=NULL;
  double *eval=NULL;
  double scale;
  int nn;
  double complex Z_ONE = make_double complex(1.0,0.0);
  double complex Z_ZERO = make_double complex(0.0,0.0);

  nb = magma_get_zhetrd_nb(N);
  lwork = N*N + 2*N*nb;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;
  hA = (double complex*)malloc(N*N*sizeof(double complex));
  work = (double complex*)malloc(lwork*sizeof(double complex));
  rwork = (double*)malloc(lrwork*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));
  eval = (double*)malloc(N*sizeof(double));
  ldhA = N;

  // Diagonalize matrix dB
  jobz = 'V';
  uplo = 'U';
  magma_zheevd_gpu(jobz, uplo, N, dB, lddB, eval, hA, ldhA,
      work, lwork, rwork, lrwork, iwork, liwork, &info);
  if(info != 0) {
    printf("ERROR calling magma_zheevd_gpu in eig_zhegv_gpu: info = %d\n", info);
    return;
  }

  nn = 0;
  for(i=1; i<=N; i++) {
    if(fabs(eval[i-1]) > SMALL) {
      nn = nn + 1;
      scale = 1.0/sqrt(eval[i-1]);
#ifdef _USE_CUBLAS2
      CALL( cublasZdscal(CUBLAS_HANDLE, N, &scale, &dB[IDX2F(1,i,lddB)], I_ONE) );
			cudaThreadSynchronize();
#else
      cublasZdscal(N, scale, &dB[IDX2F(1,i,lddB)], I_ONE);
#endif
    }
  }
  if(nn < N) {
    printf("ERROR: Number of linearly independent vectors = %d\n", nn);
    printf("       while size of the problem: %d\n", N);
    return;
  }

  transa = 'N';
  transb = 'N';
#ifdef _USE_CUBLAS2
  CALL( cublasZgemm(CUBLAS_HANDLE,
        CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &Z_ONE, dA, lddA, dB, lddB, &Z_ZERO, devec, lddE) );
	cudaThreadSynchronize();
#else
  cublasZgemm(transa, transb, N, N, N, Z_ONE, dA, lddA, dB, lddB, Z_ZERO, devec, lddE);
#endif
  //
  transa = 'C';
  transb = 'N';
#ifdef _USE_CUBLAS2
  CALL( cublasZgemm(CUBLAS_HANDLE,
        CUBLAS_OP_C, CUBLAS_OP_N, N, N, N, &Z_ONE, dB, lddB, devec, lddE, &Z_ZERO, dA, lddA) );
	cudaThreadSynchronize();
#else
  cublasZgemm(transa, transb, N, N, N, Z_ONE, dB, lddB, devec, lddE, Z_ZERO, dA, lddA);
#endif

  // Diagonalize transformed A
  jobz = 'V';
  uplo = 'U';
  magma_zheevd_gpu(jobz, uplo, N, dA, lddA, eval, hA, ldhA,
      work, lwork, rwork, lrwork, iwork, liwork, &info);
  if(info != 0) {
    printf("ERROR calling magma_zheevd_gpu in eig_zhegv : info = %d\n", info);
    return;
  }
  // Backtransform eigenvectors
  transa = 'N';
  transb = 'N';
#ifdef _USE_CUBLAS2
  CALL( cublasZgemm(CUBLAS_HANDLE,
        CUBLAS_OP_N,CUBLAS_OP_N, N,N,N, &Z_ONE,dB,lddB, dA,lddA, &Z_ZERO, devec,lddE) );
	cudaThreadSynchronize();
#else
  cublasZgemm(transa,transb, N,N,N, Z_ONE,dB,lddB, dA,lddA, Z_ZERO, devec,lddE);
#endif

  //for(i=1; i<=N; i++) {
  //  printf("Eigenvalues = %f\n", eval[i-1]);
  //}

  free(eval); eval=NULL;
  free(hA); hA=NULL;
  free(work); work=NULL;
  free(rwork); rwork=NULL;
  free(iwork); iwork=NULL;
}


//
void eig_zhegv_eval_magma(double complex *dA, int lddA, double complex *dB, int lddB, 
    double *eval, double complex *devec, int lddE, int N)
{
  const double SMALL = 12.220446049250313e-16;
  int I_ONE=1;
  int lwork, lrwork, liwork, info, i, nb, ldhA;
  char jobz, uplo, transa, transb;
  double complex *work=NULL, *hA=NULL;
  double *rwork=NULL;
  int *iwork=NULL;
  double scale;
  int nn;
  double complex Z_ONE = make_double complex(1.0,0.0);
  double complex Z_ZERO = make_double complex(0.0,0.0);

  nb = magma_get_zhetrd_nb(N);
  lwork = N*N + 2*N*nb;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;
  hA = (double complex*)malloc(N*N*sizeof(double complex));
  work = (double complex*)malloc(lwork*sizeof(double complex));
  rwork = (double*)malloc(lrwork*sizeof(double));
  iwork = (int*)malloc(liwork*sizeof(int));
  ldhA = N;

  // Diagonalize matrix dB
  jobz = 'V';
  uplo = 'U';
  magma_zheevd_gpu(jobz, uplo, N, dB, lddB, eval, hA, ldhA,
      work, lwork, rwork, lrwork, iwork, liwork, &info);
  if(info != 0) {
    printf("ERROR calling magma_zheevd_gpu in eig_zhegv_gpu: info = %d\n", info);
    return;
  }

  nn = 0;
  for(i=1; i<=N; i++) {
    if(fabs(eval[i-1]) > SMALL) {
      nn = nn + 1;
      scale = 1.0/sqrt(eval[i-1]);
#ifdef _USE_CUBLAS2
      CALL( cublasZdscal(CUBLAS_HANDLE, N, &scale, &dB[IDX2F(1,i,lddB)], I_ONE) );
#else
      cublasZdscal(N, scale, &dB[IDX2F(1,i,lddB)], I_ONE);
#endif
    }
  }
  if(nn < N) {
    printf("ERROR: Number of linearly independent vectors = %d\n", nn);
    printf("       while size of the problem: %d\n", N);
    return;
  }

  transa = 'N';
  transb = 'N';
#ifdef _USE_CUBLAS2
  CALL( cublasZgemm(CUBLAS_HANDLE,
        CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &Z_ONE, dA, lddA, dB, lddB, &Z_ZERO, devec, lddE) );
#else
  cublasZgemm(transa, transb, N, N, N, Z_ONE, dA, lddA, dB, lddB, Z_ZERO, devec, lddE);
#endif
  
  //
  transa = 'C';
  transb = 'N';
#ifdef _USE_CUBLAS2
  CALL( cublasZgemm(CUBLAS_HANDLE,
        CUBLAS_OP_C, CUBLAS_OP_N, N, N, N, &Z_ONE, dB, lddB, devec, lddE, &Z_ZERO, dA, lddA) );
#else
  cublasZgemm(transa, transb, N, N, N, Z_ONE, dB, lddB, devec, lddE, Z_ZERO, dA, lddA);
#endif

  // Diagonalize transformed A
  jobz = 'V';
  uplo = 'U';
  magma_zheevd_gpu(jobz, uplo, N, dA, lddA, eval, hA, ldhA,
      work, lwork, rwork, lrwork, iwork, liwork, &info);
  if(info != 0) {
    printf("ERROR calling magma_zheevd_gpu in eig_zhegv : info = %d\n", info);
    return;
  }
  // Backtransform eigenvectors
  transa = 'N';
  transb = 'N';
#ifdef _USE_CUBLAS2
  CALL( cublasZgemm(CUBLAS_HANDLE,
        CUBLAS_OP_N,CUBLAS_OP_N, N,N,N, &Z_ONE,dB,lddB, dA,lddA, &Z_ZERO, devec,lddE) );
#else
  cublasZgemm(transa,transb, N,N,N, Z_ONE,dB,lddB, dA,lddA, Z_ZERO, devec,lddE);
#endif

  //for(i=1; i<=N; i++) {
  //  printf("Eigenvalues = %f\n", eval[i-1]);
  //}

  free(hA); hA=NULL;
  free(work); work=NULL;
  free(rwork); rwork=NULL;
  free(iwork); iwork=NULL;
}

