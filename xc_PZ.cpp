// eFeFeR, December 2011

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void xc_PZ(double *rho, double *vxc, double &Exc, int NNR)
{
  const double A = -0.1423;
  const double B = 1.0529;
  const double C = 0.3334;
  const double a = 0.0311;
  const double b = -0.0480;
  const double c = 0.0020;
  const double d = -0.0116;

  const double SMALL=1.e-13;
  double ax, r, ux, uc, ex, ec, rs, rrs, bpr, rsl;
  int i;

  Exc = 0.0;
  ax = -(3.0/4.0)*pow(3.0/2.0/M_PI,2.0/3.0);

  for(i=0; i<NNR; i++) {
    r = rho[i];
    //
    if(r < -0.1) { // Is this safe?
      printf("***WARNING: Negative electron density: %18.10f\n", r);
      //abort();
    }
    //
    else if(r <= SMALL) {
      ux = 0.0;
      uc = 0.0;
      ec = 0.0;
      ex = 0.0;
    }
    //
    else {
      rs = pow( 3.0/(4.0*M_PI*r), 1.0/3.0);
      rrs = sqrt(rs);
      // Exchange
      ex = ax/rs;
      ux = 4.0/3.0*ex;
      // Correlation for rs >= 1
      if(rs >= 1.0) {
        bpr = 1.0 + B*rrs + C*rs;
        ec = A/bpr;
        uc = (1.0 + 7.0/6.0*B*rrs + 4.0/3.0*C*rs)*ec/bpr;
      }
      // Correlation for rs < 1
      else {
        rsl = log(rs);
        ec = a*rsl + b + c*rs*rsl + d*rs;
        uc = (2.0*d - c)/3.0*rs + b-a/3.0 + a*rsl + 2.0*c/3.0*rs*rsl;
      }
    }

    vxc[i] = ux + uc;
    Exc = Exc + r*(ex + ec);
  }

}

