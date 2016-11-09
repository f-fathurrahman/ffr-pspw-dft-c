// eFeFeR (20910015), January 2012

#include "common_pspw.h"

// X(1:NGW(ik),1:NX+2)
void chebyfilt(int ik, double complex *X, int degree, double lb, double ub, double *rhoe, int nnr)
{
  double e = (ub-lb)/2.0;
  double c = (ub+lb)/2.0;
  double sigma = e/(lb-ub);
  double sigma1 = sigma;
  double sigma2;
  double scale;
  int nelem;
  int i, id;
  double complex za = Z_ZERO; 

  double complex *Y=NULL, *Y1=NULL;

  Y = (double complex*)malloc(NGW[ik-1]*(NX+2)*sizeof(double complex));
  Y1 = (double complex*)malloc(NGW[ik-1]*(NX+2)*sizeof(double complex));

  nelem = NGW[ik-1]*(NX+2);

  // Y = H*X
  ApplyHam_block(X, Y, NGW[ik-1], ik, rhoe, nnr, NX+2);
  // Y <-- Y - c*X
  /*for(i=1; i<=nelem; i++) {
    Y[i-1] = Y[i-1] - c*X[i-1];
  }*/
  //
  // It seems that using BLAS Level 1 only gives
  // slight increase in performance (as expected)
  // USE BLAS Level1
  za = -c + I*0.0;
  zaxpy_(&nelem, &za, X,&I_ONE, Y,&I_ONE);
  // Y = Y*sigma1/e
  scale = sigma1/e;
  zdscal_(&nelem,&scale,Y,&I_ONE);

  for(id=2; id<=degree; id++) {
    sigma2 = 1.0/(2.0/sigma1 - sigma);
    // Y1 <-- H*Y
    ApplyHam_block(Y,Y1,NGW[ik-1], ik, rhoe, nnr, NX+2);
    // Y1 <-- (H*Y - c*Y)*2.0*sigma2/e - sigma*sigma2*X
    /*for(i=1; i<=nelem; i++) {
      Y1[i-1] = (Y1[i-1] - c*Y[i-1])*2.0*sigma2*(1.0/e) - sigma*sigma2*X[i-1];
    }*/
    //
    // Use BLAS Level 1
    //
    scale = 2.0*sigma2/e;
    zdscal_(&nelem,&scale,Y1,&I_ONE);
    //
    za = -c*2.0*sigma2/e + I*0.0;
    zaxpy_(&nelem,&za,Y,&I_ONE,Y1,&I_ONE);
    //
    za = -sigma*sigma2 + I*0.0;
    zaxpy_(&nelem,&za,X,&I_ONE,Y1,&I_ONE);

    
    // X <-- Y
    zcopy_(&nelem,Y,&I_ONE,X,&I_ONE);
    // Y <-- Y1
    zcopy_(&nelem,Y1,&I_ONE,Y,&I_ONE);
    sigma = sigma2;
  }

  // X <-- Y
  zcopy_(&nelem,Y,&I_ONE,X,&I_ONE);

  free(Y); Y=NULL;
  free(Y1); Y1=NULL;
}

