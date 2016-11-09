#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda/cublas.h>

#define M 6
#define N 5
#define IDX2F(i,j,ld) ((((j)-1)*(ld)) + ((i)-1))

static __inline__ void modify(
    double *m, int ldm, int n, int p, int q, double alpha, double beta)
{
  printf("Calling cublasDscal\n");
  printf("n-p+1 = %d\n",n-p+1);
  //cublasDscal(n-p, alpha, &m[IDX2F(p,q,ldm)], ldm);
  cublasDscal(ldm-p+1, beta, &m[IDX2F(p,q,ldm)], 1);
  printf("End calling cublasDscal\n");
}


int main(int argc, char **argv)
{
  int i, j;
  cublasStatus stat;
  double *d_A;
  double *A = 0;

  A = (double*)malloc(M*N*sizeof(*A));
  if(!A) {
    printf("Host memory allocation failed\n");
    return EXIT_FAILURE;
  }
  for(i=1; i<=M; i++) {
    for(j=1; j<=N; j++) {
      A[IDX2F(i,j,M)] = (double)((i-1)*M+j);
    }
  }

  printf("Matrix A before modify:\n");
  for(i=1; i<=M; i++) {
    for(j=1; j<=N; j++) {
      printf("%7.0f", A[IDX2F(i,j,M)]);
    }
    printf("\n");
  }

  cublasInit();
  printf("Allocating memory in GPU.\n");
  stat = cublasAlloc(M*N, sizeof(*A), (void**)&d_A);
  if(stat != CUBLAS_STATUS_SUCCESS) {
    printf("Device memory allocation failed\n");
    cublasShutdown();
    return EXIT_FAILURE;
  }
  printf("Copy matrix to GPU.\n");
  stat = cublasSetMatrix(M,N,sizeof(*A),A,M,d_A,M);
  if(stat != CUBLAS_STATUS_SUCCESS) {
    printf("Data download failed.\n");
    cublasFree(d_A);
    cublasShutdown();
    return EXIT_FAILURE;
  }

  printf("Modifying matrix A\n");
  modify(d_A, M, N, 2, 3, 16.0, 12.0);
  printf("Copy matrix to CPU\n");
  stat = cublasGetMatrix(M,N,sizeof(*A), d_A, M, A, M);
  if(stat != CUBLAS_STATUS_SUCCESS) {
    printf("Data upload failed\n");
    cublasFree(d_A);
    cublasShutdown();
    return EXIT_FAILURE;
  }

  cublasFree(d_A);
  cublasShutdown();

  printf("Matrix A after modify:\n");
  for(i=1; i<=M; i++) {
    for(j=1; j<=N; j++) {
      printf("%7.0f", A[IDX2F(i,j,M)]);
    }
    printf("\n");
  }

  free(A);

  printf("Program ended normally\n");
  return EXIT_SUCCESS;
}

