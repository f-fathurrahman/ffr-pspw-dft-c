// eFeFeR (20910015), December 2011

#include "common_pspw.h"

// real(8) rho_old(NNR,3)
// real(8) f(NNR,3)
// real(8) g(NNR,2)
void AndersonRhoMix(int iter, double beta, double *rhoe, double *rho_old,
    double *f, double *g, int NNR, double &d)
{
  // Local
  int j0, j1, j2, k;
  double d11, r1, d22, d12, r2, t1, t2;

  printf("Anderson mixing with beta = %f\n", beta);

  if(iter==1) {
    dcopy_(&NNR, rhoe,&I_ONE, &rho_old[IDX2F(1,1,NNR)],&I_ONE);
    for(k=1; k<=NNR; k++) f[IDX2F(k,1,NNR)] = 0.0;
    d = 1.0;
    return;
  }

  j0 = (iter-1)%3 + 1;
  j1 = (iter-2)%3 + 1;
  j2 = (iter-3)%3 + 1;
  d = 0.0;

  for(k=1; k<=NNR; k++) {
    f[IDX2F(k,j0,NNR)] = rhoe[k-1] - rho_old[IDX2F(k,j1,NNR)];
    d = d + SQUARE(f[IDX2F(k,j0,NNR)]);
  }

  if(iter==2) {
    for(k=1; k<=NNR; k++) {
      rhoe[k-1] = beta*rhoe[k-1] + (1.0 - beta)*rho_old[IDX2F(k,j1,NNR)];
      rho_old[IDX2F(k,j0,NNR)] = rhoe[k-1];
    }
    return;
  }

  for(k=1; k<=NNR; k++) {
    g[IDX2F(k,1,NNR)] = f[IDX2F(k,j1,NNR)] - f[IDX2F(k,j0,NNR)];
  }
  d11 = ddot_(&NNR, &g[IDX2F(1,1,NNR)],&I_ONE, &g[IDX2F(1,1,NNR)],&I_ONE);
  r1 = ddot_(&NNR, &f[IDX2F(1,j0,NNR)],&I_ONE, &g[IDX2F(1,1,NNR)],&I_ONE); 

  if(iter==3) {
    t1 = -r1/d11;
    for(k=1; k<=NNR; k++) {
      rhoe[k-1] = (1.0-t1)*rho_old[IDX2F(k,j1,NNR)] + t1*rho_old[IDX2F(k,j2,NNR)] + 
        beta*( f[IDX2F(k,j0,NNR)] + t1*g[IDX2F(k,1,NNR)] );
    }
  }
  else {
    for(k=1; k<=NNR; k++) {
      g[IDX2F(k,2,NNR)] = f[IDX2F(k,j2,NNR)] - f[IDX2F(k,j0,NNR)];
    }
    d22 = ddot_(&NNR, &g[IDX2F(1,2,NNR)],&I_ONE, &g[IDX2F(1,2,NNR)],&I_ONE);
    d12 = ddot_(&NNR, &g[IDX2F(1,1,NNR)],&I_ONE, &g[IDX2F(1,2,NNR)],&I_ONE);
    r2 = ddot_(&NNR, &f[IDX2F(1,j0,NNR)],&I_ONE, &g[IDX2F(1,2,NNR)],&I_ONE);
    t1 = d11*d22;
    t2 = t1 - d12*d12;
    if(fabs(t2) > fabs(t1)*1.0e-8) {
      t1 = (-r1*d22 + r2*d12)/t2;
      t2 = (r1*d12 - r2*d11)/t2;
    }
    else {
      t1 = 0.0;
      t2 = 0.0;
    }
    // Mix
    for(k=1; k<=NNR; k++) {
      rhoe[k-1] = (1.0-t1-t2)*rho_old[IDX2F(k,j1,NNR)] + t1*rho_old[IDX2F(k,j2,NNR)] +
        t2*rho_old[IDX2F(k,j0,NNR)] + beta*( f[IDX2F(k,j0,NNR)] +
            t1*g[IDX2F(k,1,NNR)] + t2*g[IDX2F(k,1,NNR)] );
    }
  }

  // rho_old(:,j0) = rhoe(:)
  dcopy_(&NNR, rhoe,&I_ONE, &rho_old[IDX2F(1,j0,NNR)],&I_ONE);

  printf("In AndersonRhoMix: d = %18.10f\n", d);

}

