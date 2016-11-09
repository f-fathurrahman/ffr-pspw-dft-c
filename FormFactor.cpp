// eFeFeR (20910015), November 2011

#include "common_pspw_cuda.h"

void FormFactor()
{
  // Several constant
  const double FPIBOM = 4.0*M_PI/OMEGA;
  const double SQ2PIM1 = 1.0/sqrt(2.0*M_PI);
  // Local
  int is, ir, ig;
  double *vloc=NULL, *fun=NULL;
  double vv, g_old, g_new, arg, facc, q1;

  vloc = (double*)malloc(MAXR*sizeof(double));
  fun = (double*)malloc(MAXR*sizeof(double));

  // Calculate self-energy of the ionic pseudocharges
  ESELF = 0.0;
  for(is=1; is<=NSP; is++) {
    ESELF = ESELF + ZV[is-1]*ZV[is-1]/RGAUSS[is-1]*NA[is-1];
  }
  ESELF = ESELF*SQ2PIM1;

  // Form factors of the ionic pseudocharge
  for(is=1; is<=NSP; is++) {
    facc = 0.25*TPIBA2*RGAUSS[is-1]*RGAUSS[is-1];
    for(ig=2; ig<=NG; ig++) {
      RHOPS[IDX2F(is,ig,NSP)] = -ZV[is-1]/OMEGA * exp(-facc*G[ig-1]);
    }
  }

  // Calculate form factors of pseudopotential
  // i.e. part of pseudopotential plus field of Gaussian pseudocharges
  int ll;
  for(is=1; is<=NSP; is++) {
    ll = L_LOC[is-1];
    if(ll < 1) {
      printf("ERROR in formf: ll < 1: %d\n", ll);
      abort();
    }
    // Local pseudopotential
    for(ir=1; ir<=MMAX[is-1]; ir++) {
      vv = -ZV[is-1]/R[IDX2F(ir,is,MAXR)]*erf( R[IDX2F(ir,is,MAXR)]/RGAUSS[is-1] );
      vloc[ir-1] = VION[IDX3F(ir,is,ll,MAXR,NSP)] - vv;
    }
    // G=0 component
    for(ir=1; ir<=MMAX[is-1]; ir++) {
      fun[ir-1] = vloc[ir-1]*POW3( R[IDX2F(ir,is,MAXR)] );
    }
    q1  = simpson(MMAX[is-1], fun, CLOG[is-1]);
    VPS[IDX2F(is,1,NSP)] = q1*FPIBOM;

    g_old = -1.0;
    for(ig=2; ig<=NG; ig++) {
      g_new = G[ig-1];
      if(g_new != g_old) {
        g_old = g_new;
        for(ir=1; ir<=MMAX[is-1]; ir++) {
          arg = R[IDX2F(ir,is,MAXR)]*sqrt(g_new)*TPIBA;
          fun[ir-1] = vloc[ir-1]*sin(arg)/arg*POW3( R[IDX2F(ir,is,MAXR)] );
        } // ir
        q1 = simpson(MMAX[is-1], fun, CLOG[is-1]);
      }
      VPS[IDX2F(is,ig,NSP)] = q1*FPIBOM;
    } // ig
  } // is

  // Free memory
  free(vloc); vloc=NULL;
  free(fun); fun=NULL;

}
