// eFeFeR (20910015), November 2011

#include "common_pspw.h"

double dsum(int N, double *x, int incx)
{
  int i, Nincx, Lincx;
  double res;

  if(N <= 0) {
    printf("ERROR: N <= 0 in dsum\n");
    abort();
  }

  Nincx = N*incx;
  Lincx = incx + 1;
  res = x[0];
  for(i=Lincx; i<=Nincx; i=i+incx) {
    res = res + x[i-1];
  }
  return res;
}

