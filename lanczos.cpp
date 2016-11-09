// eFeFeR (20910015), January 2012

#include "common_pspw_cuda.h"

// T(1:k,1:k)
// f(1:NGW_ik)
// v0(1:NGW_ik)
void lanczos(int ik, double complex *v0, int k, double *T, int ldT, double complex *f,
    double *rhoe, int nnr)
{
  // Local
  double complex *V=NULL;
  double complex *Hv=NULL;
  double complex *h=NULL, *s=NULL;
  double complex tz, za;
  double beta;
  int i, j, ig;
  int NGW_ik = NGW[ik-1];

  V = (double complex*)malloc(NGW_ik*k*sizeof(double complex));
  Hv = (double complex*)malloc(NGW_ik*sizeof(double complex));
  h = (double complex*)malloc(k*sizeof(double complex));
  s = (double complex*)malloc(k*sizeof(double complex));

  // Zero out T
  for(i=0; i<k*k; i++) T[i] = 0.0;
  
  zdotc_(&tz,&NGW_ik,v0,&I_ONE,v0,&I_ONE);
  beta = 1.0/sqrt(creal(tz));
  zdscal_(&NGW_ik,&beta,v0,&I_ONE);

  // Copy initial vector to history
  zcopy_(&NGW_ik, v0,&I_ONE, &V[IDX2F(1,1,NGW_ik)],&I_ONE);

  // Apply Hamiltonian
  ApplyHam(&V[IDX2F(1,1,NGW_ik)], Hv,NGW_ik,ik,rhoe,nnr);

  zdotc_(&tz,&NGW_ik,&V[IDX2F(1,1,NGW_ik)],&I_ONE,Hv,&I_ONE);
  h[0] = tz;
  T[IDX2F(1,1,ldT)] = creal(tz);

  for(ig=1; ig<=NGW_ik; ig++) {
    f[ig-1] = Hv[ig-1] - h[0]*V[IDX2F(ig,1,NGW_ik)];
  }

  // One step of reorthogonalization
  zdotc_(&s[0],&NGW_ik,&V[IDX2F(1,1,NGW_ik)],&I_ONE,f,&I_ONE);
  h[0] = h[0] + s[0];
  /*for(ig=1; ig<=NGW_ik; ig++) {
    f[ig-1] = f[ig-1] - s[0]*V[IDX2F(ig,1,NGW_ik)];
  }*/
  //za.x = -s[0].x; za.y = -s[0].y;
  za = -s[0];
  zaxpy_(&NGW_ik,&za, &V[IDX2F(1,1,NGW_ik)],&I_ONE,f,&I_ONE);

  for(j=2; j<=k; j++) {
    // T(j,j-1) = norm(f)
    zdotc_(&tz,&NGW_ik,f,&I_ONE,f,&I_ONE);
    T[IDX2F(j,j-1,ldT)] = sqrt(creal(tz));
    // V(:,j) = v = f/norm(f)
    beta = 1.0/sqrt(creal(tz));
    zcopy_(&NGW_ik,f,&I_ONE,&V[IDX2F(1,j,NGW_ik)],&I_ONE);
    zdscal_(&NGW_ik,&beta,&V[IDX2F(1,j,NGW_ik)],&I_ONE);
    // Hv = H*v = H*V(:,j)
    ApplyHam(&V[IDX2F(1,j,NGW_ik)], Hv,NGW_ik,ik,rhoe,nnr);
    //
    zcopy_(&NGW_ik,Hv,&I_ONE,f,&I_ONE);
    for(i=1; i<=j; i++) {
      zdotc_(&h[i-1], &NGW_ik,&V[IDX2F(1,i,NGW_ik)],&I_ONE, Hv,&I_ONE);
      /*for(ig=1; ig<=NGW_ik; ig++) {
        f[ig-1] = f[ig-1] - h[i-1]*V[IDX2F(ig,i,NGW_ik)];
      }*/
      //za.x = -h[i-1].x; za.y = -h[i-1].y;
      za = -h[i-1];
      zaxpy_(&NGW_ik,&za, &V[IDX2F(1,i,NGW_ik)],&I_ONE,f,&I_ONE);
    }
    // One step of reorthogonalization
    for(i=1; i<=j; i++) {
      zdotc_(&s[i-1], &NGW_ik,&V[IDX2F(1,i,NGW_ik)],&I_ONE, f,&I_ONE);
    }
    for(i=1; i<=j; i++) {
      h[i-1] = h[i-1] + s[i-1];
    }
    for(i=1; i<=j; i++) {
      /*for(ig=1; ig<=NGW_ik; ig++) {
        f[ig-1] = f[ig-1] - s[i-1]*V[IDX2F(ig,i,NGW_ik)];
      }*/
      //za.x = -s[i-1].x; za.y = -s[i-1].y;
      za = -s[i-1];
      zaxpy_(&NGW_ik,&za, &V[IDX2F(1,i,NGW_ik)],&I_ONE,f,&I_ONE);
    }
    for(i=1; i<=j; i++) {
      T[IDX2F(i,j,ldT)] = creal(h[i-1]);
    }
  }

  free(V); V=NULL;
  free(Hv); Hv=NULL;
  free(h); h=NULL;
  free(s); s=NULL;
}

