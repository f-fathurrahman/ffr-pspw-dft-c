// eFeFeR (20910015), December 2011

#include "common_pspw.h"

void GenPrec(int ik, double *prec, int nbasis, double tau)
{
  // Local
  int ig;
  double x,y;

  if(tau <= 0.0) {
    for(ig=1; ig<=nbasis; ig++) {
      x = XKG[IDX2F(ig,ik,NGWX)]*SQUARE(2.0*M_PI/ALAT);
      y = 27.0 + x*(18.0 + x*(12.0 +8.0*x));
      prec[ig-1] = y/(y + 16.0*x*x*x*x);
    }
  }
  else {
    for(ig=1; ig<=nbasis; ig++) {
      prec[ig-1] = tau/(0.5*(XKG[IDX2F(ig,ik,NGWX)]+1.0)*TPIBA2);
    }
  }

}

