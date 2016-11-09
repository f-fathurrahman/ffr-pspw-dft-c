// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void SimpleRhoMix(double beta, double *rhoe, double *rho_old, int nnr, double &d)
{
  // Local
  int ir;

  printf("Simple charge density mixing with mixing parameter %15.10f\n", beta);
  d = 0.0;
  for(ir=1; ir<=nnr; ir++) {
    d = d + (rhoe[ir-1] - rho_old[ir-1])*(rhoe[ir-1] - rho_old[ir-1]);
    rhoe[ir-1] = beta*rhoe[ir-1] + (1.0 - beta)*rho_old[ir-1];
    rho_old[ir-1] = rhoe[ir-1];
  }
  d = sqrt(d/(double)nnr);

  printf("In SimpleRhoMix: d = %f\n", d);
}

