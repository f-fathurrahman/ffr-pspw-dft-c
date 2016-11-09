// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void Diagon_davidson_cuda(int ik, double *eigval, double complex *eigvec, int NGW_ik, double tau)
{
  // Local variables
  double complex *eigvec_tmp=NULL;
  double *prec=NULL, *d_prec=NULL;
  int ist, ig;

  prec = (double*)malloc(NGW_ik*sizeof(double));
  CALL( cudaMalloc((void**)&eigvec_tmp, NGW_ik*NX*sizeof(double complex)) );
  CALL( cudaMalloc((void**)&d_prec, NGW_ik*sizeof(double)) );

  GenPrec(ik,prec,NGW_ik,tau);
  CALL( cudaMemcpy(d_prec, prec, NGW_ik*sizeof(double), cudaMemcpyHostToDevice) );

  davidson_cuda(ik, eigval, eigvec_tmp, NGW_ik, NX, d_prec);


  // Assign eigenvalues and eigenvectors
  CALL( cublasGetMatrix(NGW_ik,NX,sizeof(double complex),eigvec_tmp,NGW_ik,eigvec,NGWX) );
  for(ist=1; ist<=NX; ist++) {
    for(ig=NGW_ik+1; ig<=NGWX; ig++) eigvec[IDX2F(ig,ist,NGWX)] = Z_ZERO;
  }

  free(prec); prec=NULL;
  CALL( cudaFree(eigvec_tmp) );
  CALL( cudaFree(d_prec) );
}

void Diagon_lobpcg_cuda(int ik, double *eigval, double complex *eigvec, int NGW_ik, double tau)
{
  // Local variables
  double complex *eigvec_tmp=NULL;
  double *prec=NULL, *d_prec=NULL;
  int ist, ig;

  prec = (double*)malloc(NGW_ik*sizeof(double));
  CALL( cudaMalloc((void**)&eigvec_tmp, NGW_ik*NX*sizeof(double complex)) );
  CALL( cudaMalloc((void**)&d_prec, NGW_ik*sizeof(double)) );

  GenPrec(ik,prec,NGW_ik,tau);
  CALL( cudaMemcpy(d_prec, prec, NGW_ik*sizeof(double), cudaMemcpyHostToDevice) );

  lobpcg_cuda(ik, eigval, eigvec_tmp, NGW_ik, NX, d_prec);


  // Assign eigenvalues and eigenvectors
  CALL( cublasGetMatrix(NGW_ik,NX,sizeof(double complex),eigvec_tmp,NGW_ik,eigvec,NGWX) );
  for(ist=1; ist<=NX; ist++) {
    for(ig=NGW_ik+1; ig<=NGWX; ig++) eigvec[IDX2F(ig,ist,NGWX)] = Z_ZERO;
  }

  free(prec); prec=NULL;
  CALL( cudaFree(eigvec_tmp) );
  CALL( cudaFree(d_prec) );
}

