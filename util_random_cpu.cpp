// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

int myrandom_cpu(long &a)
{
  int lim = 100;
  
  a = (((a * 214013L + 2531011L) >> 16) & 32767);
  
  return ((a % lim) + 1);
}

void random_col(double complex *A, int nelem, long a)
{
  long seed=a;
  double A_re, A_im;
  for(int i=0; i<nelem; i++) {
    A_re = myrandom_cpu(seed)/(double)seed;
    A_im = myrandom_cpu(seed)/(double)seed;
    A[i] = A_re + I*A_im;
  }
}


