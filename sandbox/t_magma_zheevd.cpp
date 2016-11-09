#include <stdio.h>
#include <stdlib.h>
#include <cuda/cuComplex.h>
#include "../magma/include/magma.h"
#include "../magma/include/magma_lapack.h"

#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))

typedef cuDoubleComplex complex16;

static void HandleError( cudaError_t err,const char *file,int line )
{
  if (err != cudaSuccess) {
    printf("%s in %s at line %d\n", cudaGetErrorString(err),file,line);
    exit(EXIT_FAILURE);
  }
}

#define CALL( err ) (HandleError( err, __FILE__, __LINE__ ))

#define HANDLE_NULL( a ) {if (a == NULL) { \
                            printf( "Host memory failed in %s at line %d\n", \
                                    __FILE__, __LINE__ ); \
                            exit( EXIT_FAILURE );}}

int main(int argc, char **argv)
{
  int N=4;
  complex16 *A=NULL, *d_A;
  double *w=NULL, *d_w;
  complex16 *d_work=NULL;
  double *d_rwork=NULL;
  int *d_iwork=NULL;
  int lwork, lrwork, liwork;
  char jobz='V', uplo='U';
  int info;
  int i,j;

  magma_int_t = magma_get_zhetrd_nb(N);
  lwork = 2*N*nb + N*N;
  lrwork = 2*N*N + 5*N + 1;
  liwork = 5*N + 3;

//
// Allocate memory
//
  // CPU
  A = (complex16*)malloc( N*N*sizeof(complex16) );
  w = (double*)malloc( N*sizeof(double) );
  // GPU
  CALL( cudaMalloc((void**)&d_A, N*N*sizeof(complex16)) );
  CALL( cudaMalloc((void**)&d_w, N*sizeof(double)) );
  CALL( cudaMalloc((void**)&d_work, lwork*sizeof(complex16)) );
  CALL( cudaMalloc((void**)&d_rwork, lrwork*sizeof(double)) );
  CALL( cudaMalloc((void**)&d_iwork, liwork*sizeof(int)) );

  // Initialize matrix A
  A[0] = make_cuDoubleComplex(1.0,0.0);
  A[1] = make_cuDoubleComplex(2.0,-1.0);
  A[2] = make_cuDoubleComplex(3.0,-1.0);
  A[3] = make_cuDoubleComplex(4.0,-1.0);
  //
  A[4] = make_cuDoubleComplex(2.0,1.0);
  A[5] = make_cuDoubleComplex(2.0,0.0);
  A[6] = make_cuDoubleComplex(3.0,-2.0);
  A[7] = make_cuDoubleComplex(4.0,-2.0);
  //
  A[8] = make_cuDoubleComplex(3.0,1.0);
  A[9] = make_cuDoubleComplex(3.0,2.0);
  A[10] = make_cuDoubleComplex(3.0,0.0);
  A[11] = make_cuDoubleComplex(4.0,-3.0);
  //
  A[12] = make_cuDoubleComplex(4.0,1.0);
  A[13] = make_cuDoubleComplex(4.0,2.0);
  A[14] = make_cuDoubleComplex(4.0,3.0);
  A[15] = make_cuDoubleComplex(4.0,0.0);

  printf("Matrix A:\n");
  for(i=1; i<=N; i++) {
    for(j=1; j<=N; j++) {
      printf("(%f,%f) ", A[IDX2F(i,j,N)].x, A[IDX2F(i,j,N)].y);
    }
    printf("\n");
  }

  CALL( cudaMemcpy(d_A, A, N*N*sizeof(complex16), cudaMemcpyHostToDevice) );

  magma_zheevd(jobz, uplo, N, A, N, w, work, lwork, rwork, lrwork, iwork, liwork, &info);
  if(info != 0) {
    printf("Error calling magma_zheevd\n");
    abort();
  }

  CALL( cudaMemcpy(w, d_w, N*sizeof(double), cudaMemcpyDeviceToHost) );
  CALL( cudaMemcpy(A, d_A, N*N*sizeof(complex16), cudaMemcpyDeviceToHost) );
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

  CALL( cudaFree(d_A) );
  CALL( cudaFree(d_w) );
  CALL( cudaFree(d_work) );
  CALL( cudaFree(d_rwork) );
  CALL( cudaFree(d_iwork) );

  free(A); A=NULL;
  free(w); w=NULL;

  return 0;
}  
