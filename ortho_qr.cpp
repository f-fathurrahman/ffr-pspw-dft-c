#include "common_pspw.h"

void ortho_qr(double complex *X, int nbasis, int nstates)
{
  // Local
  int m, n, k, lwork, info;
  double complex *tau=NULL, *work=NULL;

  m = nbasis;
  n = nstates;
  k = nstates;
  lwork = n;

  // Allocate memory for arrays
  tau = (double complex*)malloc(k*sizeof(double complex));
  work = (double complex*)malloc(lwork*sizeof(double complex));

  // Compute QR factorization
  zgeqrf_(&m,&n,X,&m,tau,work,&lwork,&info);
  if(info != 0) {
    printf("ERROR calling zgeqrf_ : info = %d\n",info);
    return;
  }

  // Construct Q explicitly
  zungqr_(&m,&k,&k,X,&m,tau,work,&lwork,&info);
  if(info != 0) {
    printf("ERROR calling zungqr_ : info = %d\n", info);
    return;
  }

  free(tau); tau=NULL;
  free(work); work=NULL;
}

