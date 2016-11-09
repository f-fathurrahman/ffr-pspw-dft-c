// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void Diagon_lobpcg_cpu(int ik, double *eigval, double complex *eigvec, int NGW_ik,
    double *rhoe, int nnr, double tau)
{
  // Local variables
  double complex *eigvec_tmp=NULL;
  double *prec=NULL;
  int ist, ig;

  eigvec_tmp = (double complex*)malloc(NGW_ik*NX*sizeof(double complex));
  prec = (double*)malloc(NGW_ik*sizeof(double));

  //for(ist=1; ist<=NGW_ik; ist++) prec[ist-1] = 1.0;
  GenPrec(ik,prec,NGW_ik,tau);

  lobpcg_cpu(ik, eigval, eigvec_tmp, NGW_ik, NX, rhoe, nnr, prec);

  // Assign eigenvalues and eigenvectors
  for(ist=1; ist<=NX; ist++) {
    zcopy_(&NGW_ik, &eigvec_tmp[IDX2F(1,ist,NGW_ik)],&I_ONE,&eigvec[IDX2F(1,ist,NGWX)],&I_ONE);
    for(ig=NGW_ik+1; ig<=NGWX; ig++) eigvec[IDX2F(ig,ist,NGWX)] = Z_ZERO;
  }

  free(eigvec_tmp); eigvec_tmp = NULL;
  free(prec); prec=NULL;

}

void Diagon_davidson_cpu(int ik, double *eigval, double complex *eigvec, int NGW_ik,
    double *rhoe, int nnr, double tau)
{
  // Local variables
  double complex *eigvec_tmp=NULL;
  double *prec=NULL;
  int ist, ig;

  eigvec_tmp = (double complex*)malloc(NGW_ik*NX*sizeof(double complex));
  prec = (double*)malloc(NGW_ik*sizeof(double));

  //for(ist=1; ist<=NGW_ik; ist++) prec[ist-1] = 1.0;
  GenPrec(ik,prec,NGW_ik,tau);

  davidson_cpu(ik, eigval, eigvec_tmp, NGW_ik, NX, rhoe, nnr, prec);

  // Assign eigenvalues and eigenvectors
  for(ist=1; ist<=NX; ist++) {
    zcopy_(&NGW_ik, &eigvec_tmp[IDX2F(1,ist,NGW_ik)],&I_ONE,&eigvec[IDX2F(1,ist,NGWX)],&I_ONE);
    for(ig=NGW_ik+1; ig<=NGWX; ig++) eigvec[IDX2F(ig,ist,NGWX)] = Z_ZERO;
  }

  free(eigvec_tmp); eigvec_tmp = NULL;
  free(prec); prec=NULL;

}

