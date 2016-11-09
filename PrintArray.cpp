// eFeFeR (20910015), October 2011

#include "common_pspw.h"

// Several utilities for printing arrays
//
// For this time, I don't use template function because it
// will be difficult to implement. I am using printf, not cout :)

//
// Print full double complex array
//
void PrintVector(double complex *V, int nelem)
{
  printf("nelem=%d\n",nelem);
  for(int i=0; i<nelem; i++) {
    printf("(%f,%f)\n", creal(V[i]), cimag(V[i]));
  }
}

//
// Print full double array
//
void PrintVector(double *V, int nelem)
{
  printf("nelem=%d\n",nelem);
  for(int i=0; i<nelem; i++) {
    printf("%f\n", V[i]);
  }
}

//
// Print full double complex matrix
//
void PrintMatrix(double complex *A, int NROW, int NCOL)
{
  printf("NROW=%d, NCOL=%d\n", NROW, NCOL);
  for(int i=1; i<=NROW; i++) {
    for(int j=1; j<=NCOL; j++) {
      printf("(%f,%f) ", creal( A[IDX2F(i,j,NROW)]), cimag(A[IDX2F(i,j,NROW)]) );
    }
    printf("\n");
  }
}

//
// Print full double matrix
//
void PrintMatrix(double *A, int NROW, int NCOL)
{
  printf("NROW=%d, NCOL=%d\n", NROW, NCOL);
  for(int i=1; i<=NROW; i++) {
    for(int j=1; j<=NCOL; j++) {
      printf("%f ", A[IDX2F(i,j,NROW)]);
    }
    printf("\n");
  }
}
