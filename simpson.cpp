// eFeFeR (20910015), November 2011

#include "common_pspw.h"

double simpson(int N, double *F, double h)
{
  int nm12 = (N-1)/2;
  if(nm12*2 != N-1) {
    printf("ERROR: N must be odd in simpson\n");
    abort();
  }

  double simps = 4.0*dsum(nm12,&F[1],2) + 2.0*dsum(nm12-1,&F[2],2) + F[0] + F[N-1];
  return simps*h/3.0;
}
