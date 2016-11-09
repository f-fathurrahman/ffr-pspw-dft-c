// eFeFeR (20910015), November 2011

#include "common_pspw.h"

void EwaldEnergy()
{
  int i, j, k, ir, is, ia, is1, ia1, ir2, is2, ia2;
  double gkl, erre2, arg, fac, addesr, addpre;
  double repand, rcut2, r2, sumh;
  double tau1[3], tau2[3], tau3[3], t[3];
  double *rlat=NULL, *tau_ws=NULL;
  int nrvec, n_cell;
  const int NRX=40000;
  const double RCUT=30.0;
  double xmat[9], xmati[9];
  bool t_nearest_neighbor_cell;

  rlat = (double*)malloc(NRX*3*sizeof(double));
  tau_ws = (double*)malloc(3*NAX*NSP*sizeof(double));

  for(i=0; i<NRX*3; i++) rlat[i] = 0.0;
  for(i=0; i<3*NAX*NSP; i++) tau_ws[i] = 0.0;

  // Determine reciprocal lattice vectors
  for(i=1; i<=3; i++) {
    xmat[IDX2F(i,1,3)] = A1[i-1];
    xmat[IDX2F(i,2,3)] = A2[i-1];
    xmat[IDX2F(i,3,3)] = A3[i-1];
  }
  matinv3x3(xmat, xmati);

  // Determine atomic positions in lattice coordinates
  for(is=1; is<=NSP; is++) {
  for(ia=1; ia<=NA[is-1]; ia++) {
    for(i=1; i<=3; i++) {
      sumh =        xmati[IDX2F(i,1,3)]*TAU[IDX3F(1,ia,is,3,NAX)];
      sumh = sumh + xmati[IDX2F(i,2,3)]*TAU[IDX3F(2,ia,is,3,NAX)];
      sumh = sumh + xmati[IDX2F(i,3,3)]*TAU[IDX3F(3,ia,is,3,NAX)];
      tau1[i-1] = sumh - floor(sumh); // TODO: Check this
    }
    for(i=1; i<=3; i++) {
      tau_ws[IDX3F(i,ia,is,3,NAX)] = A1[i-1]*tau1[0] + A2[i-1]*tau1[1] + A3[i-1]*tau1[2];
    }
  }
  }

  // Calculate lattice vectors with length less than \verb|RCUT|
  rcut2 = RCUT*RCUT;
  nrvec = 0;
  n_cell = max( (int)round(5.0*RCUT/pow(OMEGA,1.0/3.0)), 2);

  for(i=-n_cell; i<=n_cell; i++) {
  for(j=-n_cell; j<=n_cell; j++) {
  for(k=-n_cell; k<=n_cell; k++) {
    r2 = 0.0;
    for(ir=1; ir<=3; ir++) {
      t[ir-1] = (double)i*A1[ir-1] + (double)j*A2[ir-1] + (double)k*A3[ir-1];
      r2 = r2 + t[ir-1]*t[ir-1];
    }
    t_nearest_neighbor_cell = (abs(i) <= 1) && (abs(j) <= 1) && (abs(k) <= 1);
    if(r2 <= rcut2 || t_nearest_neighbor_cell) {
      nrvec = nrvec + 1;
      if(nrvec > NRX) {
        printf("ERROR in EwaldEnergy: nrvec > nrx = %d > %d\n", nrvec,NRX);
        abort();
      }
      rlat[IDX2F(nrvec,1,NRX)] = t[0];
      rlat[IDX2F(nrvec,2,NRX)] = t[1];
      rlat[IDX2F(nrvec,3,NRX)] = t[2];
    }
  }
  }
  }

  // Ewald summation with setting = 0 of forces
  ESR = 0.0;

  for(is1=1; is1<=NSP; is1++) {
    for(ia1=1; ia1<=NA[is1-1]; ia1++) {
      for(i=1; i<=3; i++) tau1[i-1] = tau_ws[IDX3F(i,ia1,is1,3,NAX)];
      for(ir2=1; ir2<=nrvec; ir2++) {
        for(is2=1; is2<=NSP; is2++) {
          gkl = 1.0/sqrt(SQUARE(RGAUSS[is1-1]) + SQUARE(RGAUSS[is2-1]));
          for(ia2=1; ia2<=NA[is2-1]; ia2++) {
            for(i=1; i<=3; i++) tau2[i-1] = rlat[IDX2F(ir2,i,NRX)] + tau_ws[IDX3F(i,ia2,is2,3,NAX)];
            if(tau1[0]==tau2[0] && tau1[1]==tau2[1] && tau1[2]==tau2[2]) continue;
            for(i=0; i<3; i++) tau3[i] = tau1[i] - tau2[i];
            erre2 = 0.0;
            for(i=0; i<3; i++) {
              erre2 = erre2 + SQUARE(tau3[i]);
            }
            arg = sqrt(erre2);
            fac = ZV[is1-1]*ZV[is2-1]/arg*0.5;
            arg = arg*gkl;
            addesr = erf(arg);
            addesr = (1.0 - addesr)*fac;
            ESR = ESR + addesr;
						//printf("ESR = %18.10f\n", ESR);
            addpre = exp(-arg*arg)*gkl;
            addpre = 1.0/sqrt(M_PI)*ZV[is1-1]*ZV[is2-1]*addpre;
            repand = (addesr + addpre)/erre2;
            } // ia2
        } // is2
      } // ir2
    } // ia1
  } //is1

  free(rlat); rlat=NULL;
  free(tau_ws); tau_ws=NULL;
}

