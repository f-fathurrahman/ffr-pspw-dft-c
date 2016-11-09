// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void ortho_qr_magma(double complex *dX, int nbasis, int nstates)
{
  // Local
  int m, n, k, lwork, info, dimT, nb;
  double complex *tau=NULL;
  double complex *dT;

  m = nbasis;
  n = nstates;
  k = nstates;
  lwork = n;

  nb = magma_get_zgeqrf_nb(m);
  dimT = (2*k + (n+31)/32*32)*nb;
  // Allocate memory for arrays
  tau = (double complex*)malloc(k*sizeof(double complex));
  // TODO: Guard this statement
  cudaMalloc((void**)&dT,dimT*sizeof(double complex));

  // Compute QR factorization
  magma_zgeqrf_gpu(m,n,dX,m,tau,dT,&info);
  if(info != 0) {
    printf("ERROR calling magma_zgeqrf_gpu : info = %d\n",info);
    return;
  }

  // Construct Q explicitly
  magma_zungqr_gpu(m,n,k,dX,m,tau,dT,nb,&info);
  if(info != 0) {
    printf("ERROR calling magma_zungqr_gpu : info = %d\n", info);
    return;
  }

  free(tau); tau=NULL;
  cudaFree(dT);
}

