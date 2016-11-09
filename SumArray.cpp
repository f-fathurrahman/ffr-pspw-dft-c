// eFeFeR (20910015), November 2011

#include "common_pspw_cuda.h"

void SumArray(int *a, int N, int &sum)
{
  int i;
  sum = 0;
  for(i=0; i<N; i++) {
    sum = sum + a[i];
  }
}

void SumArray(unsigned int *a, unsigned int N, unsigned int &sum)
{
  int i;
  sum = 0;
  for(i=0; i<N; i++) {
    sum = sum + a[i];
  }
}

void SumArray(double *a, int N, double &sum)
{
  int i;
  sum = 0.0;
  for(i=0; i<N; i++) {
    sum = sum + a[i];
  }
}

void SumArray(double complex *a, int N, double complex &sum)
{
  int i;
  sum = Z_ZERO;
  for(i=0; i<N; i++) {
    sum = sum + a[i];
  }
}

