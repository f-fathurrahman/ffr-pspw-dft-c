// eFeFeR (20910015), January 2012

#include "common_pspw.h"

void PulayRhoMix(int iter, double *rhoe, double *rho_old, double *f, int mixdim, int nnr,
    double &d)
{
  // Local
  int jc,jn,i,j,m,info;
  double beta=0.25;
  int *ipiv=NULL;
  double *alpha=NULL, *a=NULL, *work=NULL;

  printf("Pulay mixing with MIXDIM = %d\n", mixdim);

  if(mixdim < 2) {
    printf("ERROR in PulayRhoMix: mixdim < 2 : %d\n", mixdim);
    exit(1);
  }

  d = 0.0;
  if(iter < mixdim) {
    for(i=1; i<=nnr; i++) {
      rhoe[i-1] = beta*rhoe[i-1] + (1.0-beta)*rho_old[IDX2F(i,iter,nnr)];
      f[IDX2F(i,iter+1,nnr)] = rhoe[i-1] - rho_old[IDX2F(i,iter,nnr)];
    }
    KerkerPrec(f, iter+1, nnr, 1.0, KERKER_Q0); // TODO: PRECONDITION?
    for(i=1; i<=nnr; i++) {
      rho_old[IDX2F(i,iter+1,nnr)] = rhoe[i-1];
      d = d + f[IDX2F(i,iter+1,nnr)]*f[IDX2F(i,iter+1,nnr)];
    }
    d = sqrt(d/(double)nnr);
    printf("In PulayRhoMix: d = %f\n",d);
    return;
  }

  jc = iter%mixdim + 1; // current index
  jn = (iter+1)%mixdim + 1; // next index
  m = min(iter+1,mixdim) + 1; // matrix size

  ipiv = (int*)malloc(m*sizeof(int));
  alpha = (double*)malloc(m*sizeof(double));
  a = (double*)malloc(m*m*sizeof(double));
  work = (double*)malloc(m*sizeof(double));

  // Compute f and RMS difference
  d = 0.0;
  for(i=1; i<=nnr; i++) {
    f[IDX2F(i,jc,nnr)] = rhoe[i-1] - rho_old[IDX2F(i,jc,nnr)];
  }
  KerkerPrec(f, jc, nnr, 1.0, KERKER_Q0); // TODO: PRECONDITION?
  for(i=1; i<=nnr; i++) {
    d = d + f[IDX2F(i,jc,nnr)]*f[IDX2F(i,jc,nnr)];
  }
  d = sqrt(d/(double)nnr);

  // Solve the linear system
  for(i=0; i<m*m; i++) a[i] = 0.0;
  for(i=1; i<=m-1; i++) {
    for(j=i; j<=m-1; j++) {
      a[IDX2F(i,j,m)] = a[IDX2F(i,j,m)] +
        ddot_(&NNR,&f[IDX2F(1,i,nnr)],&I_ONE,&f[IDX2F(1,j,nnr)],&I_ONE);
    }
    a[IDX2F(i,m,m)] = 1.0;
  }
  for(i=0; i<m; i++) alpha[i] = 0.0;
  alpha[m-1] = 1.0;

  char U='U';
  dsysv_(&U,&m,&I_ONE, a,&m, ipiv, alpha,&m, work,&m, &info);
  if(info != 0) {
    printf("ERROR calling dsysv in PulayRhoMix: info = %d\n",info);
    exit(1);
  }

  for(i=1; i<=nnr; i++) rhoe[i-1] = 0.0;
  for(j=1; j<=m-1; j++) {
    for(i=1; i<=nnr; i++) {
      rhoe[i-1] = rhoe[i-1] + alpha[j-1]*(rho_old[IDX2F(i,j,nnr)] + f[IDX2F(i,j,nnr)]);
    }
  }
  NormalizeRho(rhoe, NNR, OMEGA, NEL);
  for(i=1; i<=nnr; i++) rho_old[IDX2F(i,jn,nnr)] = rhoe[i-1];

  printf("In PulayRhoMix: d = %f\n", d);

  free(ipiv); ipiv=NULL;
  free(alpha); alpha=NULL;
  free(a); a=NULL;
  free(work); work=NULL;
}

