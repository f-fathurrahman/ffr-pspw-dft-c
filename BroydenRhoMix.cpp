// eFeFeR (20910015), January 2012

#include "common_pspw_cuda.h"

// real(8) rhoe_old(nnr,2)
// real(8) f(nnr,2)
// real(8) df(nnr,mixdim)
// real(8) u(nnr,mixdim)
// real(8) a(mixdim,mixdim)
void BroydenRhoMix(int iter, double alpha, double w0, double *rho, double *rho_old,
    double *f, double *df, double *u, double *a, int mixdim, int nnr, double &d)
{
  // Local
  double t1, sum;
  int i,j,k,l,m,info,jc,kp,kc;
  int *ipiv;
  double *c, *beta, *gamma, *work;

  c = (double*)malloc(mixdim*sizeof(double));
  beta = (double*)malloc(mixdim*mixdim*sizeof(double));
  gamma = (double*)malloc(mixdim*sizeof(double));

  ipiv = (int*)malloc(mixdim*sizeof(int));
  work = (double*)malloc(mixdim*sizeof(double));

  printf("Broyden mixing with MIXDIM=%5d\n",mixdim);
  printf("alpha = %8.5f w0 = %8.5f\n", alpha, w0);

  if(mixdim < 2) {
    printf("ERROR in BroydenRhoMix: mixdim < 2: %d\n", mixdim);
    exit(1);
  }

  m = min(iter+1,mixdim); // current subspace dimension
  jc = iter%m + 1; // current index modulo m
  kp = (iter-1)%2 + 1; // previous index modulo 2
  kc = iter%2 + 1; // current index modulo 2

  for(i=1; i<=nnr; i++) f[IDX2F(i,kc,nnr)] = rho[i-1] - rho_old[IDX2F(i,kp,nnr)];
  d = 0.0;
  for(i=1; i<=nnr; i++) {
    d = d + f[IDX2F(i,kc,nnr)]*f[IDX2F(i,kc,nnr)];
  }
  d = sqrt(d/(double)nnr);

  for(i=1; i<=nnr; i++) df[IDX2F(i,jc,nnr)] = f[IDX2F(i,kc,nnr)] - f[IDX2F(i,kp,nnr)];
  t1 = dnrm2_(&nnr, &df[IDX2F(1,jc,nnr)],&I_ONE);
  if(t1 > 1.e-8) t1 = 1.0/t1;

  for(i=1; i<=nnr; i++) {
    df[IDX2F(i,jc,nnr)] = t1*df[IDX2F(i,jc,nnr)];
    u[IDX2F(i,jc,nnr)] = alpha*df[IDX2F(i,jc,nnr)] + t1*(rho_old[IDX2F(i,kp,nnr)] - rho_old[IDX2F(i,kc,nnr)]);
  }

  for(k=1; k<=m; k++) {
    c[k-1] = ddot_(&NNR, &df[IDX2F(1,k,nnr)],&I_ONE, &f[IDX2F(1,kc,nnr)],&I_ONE);
  }
  
  for(k=1; k<=m; k++) {
    a[IDX2F(k,jc,mixdim)] = ddot_(&NNR, &df[IDX2F(1,jc,nnr)],&I_ONE, &df[IDX2F(1,k,nnr)],&I_ONE);
    a[IDX2F(jc,k,mixdim)] = a[IDX2F(k,jc,mixdim)];
  }

  for(i=0; i<mixdim*mixdim; i++) beta[i] = a[i];
  for(k=1; k<=m; k++) {
    beta[IDX2F(k,k,mixdim)] = beta[IDX2F(k,k,mixdim)] + w0*w0;
  }

  // Invert beta
  dgetrf_(&m,&m,beta,&mixdim,ipiv,&info);
  if(info != 0) {
    printf("ERROR calling dgetrf_ in BroydenRhoMix, info=%d\n",info);
    exit(1);
  }
  dgetri_(&m,beta,&mixdim,ipiv,work,&m,&info);
  if(info != 0) {
    printf("ERROR calling dgetri_ in BroydenRhoMix, info=%d\n",info);
    exit(1);
  }

  for(l=1; l<=m; l++) {
    gamma[l-1] = 0.0;
    for(k=1; k<=m; k++) {
      gamma[l-1] = gamma[l-1] + c[k-1]*beta[IDX2F(k,l,mixdim)];
    }
  }

  for(i=1; i<=nnr; i++) rho[i-1] = rho_old[IDX2F(i,kp,nnr)] + alpha*f[IDX2F(i,kc,nnr)];
  for(l=1; l<=m; l++) {
    t1 = -gamma[l-1];
    for(i=1; i<=nnr; i++) {
      rho[i-1] = rho[i-1] + t1*u[IDX2F(i,l,nnr)];
    }
  }

  for(i=1; i<=nnr; i++) rho_old[IDX2F(i,kc,nnr)] = rho[i-1];

  printf("In BroydenMix: d = %f\n", d);

  free(c); c=NULL;
  free(beta); beta=NULL;
  free(gamma); gamma=NULL;
  free(ipiv); ipiv=NULL;
  free(work); work=NULL;
  
  return;
}

