#include <stdio.h>
#include <stdlib.h>
#include <cuda/cuComplex.h>

#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))

typedef cuDoubleComplex complex16;

extern "C" void zheevd_(char *jobz, char *uplo, int *N, 
  complex16 *A, int *lda, double *w,
  complex16 *work, int *lwork,
  double *rwork, int *lrwork,
  int *iwork, int *liwork, int *info);

int main(int argc, char **argv)
{
  int N=4;
  complex16 *A=NULL, *work=NULL;
  double *w=NULL;
  double *rwork=NULL;
  int *iwork=NULL;
  int lwork, lrwork, liwork;
  char jobz='V', uplo='U';
  int info;
  int i,j;

  lwork = N*N + 2*N;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;

//
// Allocate memory
//
  A = (complex16*)malloc( N*N*sizeof(complex16) );
  w = (double*)malloc( N*sizeof(double) );
  work = (complex16*)malloc( lwork*sizeof(complex16) );
  rwork = (double*)malloc( lrwork*sizeof(double) );
  iwork = (int*)malloc( liwork*sizeof(int) );

  // Initialize matrix A
  A[0] = make_cuDoubleComplex(1.0,0.0);
  A[1] = make_cuDoubleComplex(2.0,0.0);
  A[2] = make_cuDoubleComplex(3.0,0.0);
  A[3] = make_cuDoubleComplex(4.0,0.0);
  //
  A[4] = make_cuDoubleComplex(2.0,0.0);
  A[5] = make_cuDoubleComplex(2.0,0.0);
  A[6] = make_cuDoubleComplex(3.0,0.0);
  A[7] = make_cuDoubleComplex(4.0,0.0);
  //
  A[8] = make_cuDoubleComplex(3.0,0.0);
  A[9] = make_cuDoubleComplex(3.0,0.0);
  A[10] = make_cuDoubleComplex(3.0,0.0);
  A[11] = make_cuDoubleComplex(4.0,0.0);
  //
  A[12] = make_cuDoubleComplex(4.0,0.0);
  A[13] = make_cuDoubleComplex(4.0,0.0);
  A[14] = make_cuDoubleComplex(4.0,0.0);
  A[15] = make_cuDoubleComplex(4.0,0.0);

  printf("Matrix A:\n");
  for(i=1; i<=N; i++) {
    for(j=1; j<=N; j++) {
      printf("(%f,%f) ", A[IDX2F(i,j,N)].x, A[IDX2F(i,j,N)].y);
    }
    printf("\n");
  }

  zheevd_(&jobz, &uplo, &N, A, &N, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  if(info != 0) {
    printf("Error calling zheevd\n");
    abort();
  }
  //
  printf("Eigenvalues:\n");
  for(i=1; i<=N; i++) {
    printf("%f\n", w[i-1]);
  }
  printf("Eigenvectors:\n");
  for(i=1; i<=N; i++) {
    for(j=1; j<=N; j++) {
      printf("(%f,%f) ", A[IDX2F(i,j,N)].x, A[IDX2F(i,j,N)].y);
    }
    printf("\n");
  }

  free(A); A=NULL;
  free(w); w=NULL;
  free(work); work=NULL;
  free(rwork); rwork=NULL;
  free(iwork); iwork=NULL;

  return 0;
}  
