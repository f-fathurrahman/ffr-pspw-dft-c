// eFeFeR (20910015), November 2011

#include <cstdio>
#include <cstdlib>
#include <cmath>

void NormalizeRho(double *rho, int nnr, double omega, double nel)
{
  const double SMALL=1.0e-15;
  double sum = 0.0;
  double truesum = nnr*nel/omega;
  int i;

  for(i=1; i<=nnr; i++) {
    if(rho[i-1] < 0.0) rho[i-1] = 0.0;
    sum = sum + rho[i-1];
  }

  if(sum < SMALL) {
    printf("\nERROR in NormalizeRho: sum(rho)=%f < %f\n", sum, SMALL);
    abort();
  }

  if(fabs(truesum-sum)*omega/nnr > 1.0e-5) {
    printf("\nWARNING in NormalizeRho: difference between correct and actual\n");
    printf("number of electrons: %f\n", (truesum-sum)*omega/nnr);
  }

  for(i=1; i<=nnr; i++) {
    rho[i-1] = rho[i-1]*truesum/sum;
  }

}

