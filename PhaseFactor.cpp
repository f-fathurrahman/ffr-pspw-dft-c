// eFeFeR (20910015), November 2011

#include "common_pspw.h"

void PhaseFactor()
{
  int is, ia, i, j, k, ik, ig, igp;
  double taup[3];
  double sum1;
  double xmat[9];
  double xmati[9];
  int idxzero1, idxzero2, idxzero3;
  double s,c;
  double ar1, ar2, ar3;

  idxzero1 = (NR1+1)/2 + 1;
  idxzero2 = (NR2+1)/2 + 1;
  idxzero3 = (NR3+1)/2 + 1;
  int NR11 = NR1 + 1;
  int NR21 = NR2 + 1;
  int NR31 = NR3 + 1;

  for(i=1; i<=3; i++) {
    xmat[IDX2F(i,1,3)] = A1[i-1];
    xmat[IDX2F(i,2,3)] = A2[i-1];
    xmat[IDX2F(i,3,3)] = A3[i-1];
  }
  matinv3x3(xmat, xmati);

  for(is=1; is<=NSP; is++) {
    for(ia=1; ia<=NA[is-1]; ia++) {
      for(i=1; i<=3; i++) {
        sum1 = xmati[IDX2F(i,1,3)]*TAU[IDX3F(1,ia,is,3,NAX)];
        sum1 = sum1 + xmati[IDX2F(i,2,3)]*TAU[IDX3F(2,ia,is,3,NAX)];
        sum1 = sum1 + xmati[IDX2F(i,3,3)]*TAU[IDX3F(3,ia,is,3,NAX)];
        taup[i-1] = sum1;
      }
      EI1[IDX3F(idxzero1,ia,is,NR11,NAX)] = Z_ONE;
      EI2[IDX3F(idxzero2,ia,is,NR21,NAX)] = Z_ONE;
      EI3[IDX3F(idxzero3,ia,is,NR31,NAX)] = Z_ONE;
      ar1 = -2.0*M_PI*taup[0];
      ar2 = -2.0*M_PI*taup[1];
      ar3 = -2.0*M_PI*taup[2];
      for(i=1; i<=(NR1+1)/2; i++) {
        c = cos(i*ar1);
        s = sin(i*ar1);
        EI1[IDX3F(idxzero1+i,ia,is,NR11,NAX)] = c + I*s;
        EI1[IDX3F(idxzero1-i,ia,is,NR11,NAX)] = c - I*s;
      }
      for(j=1; j<=(NR2+1)/2; j++) {
        c = cos(j*ar2);
        s = sin(j*ar2);
        EI2[IDX3F(idxzero2+j,ia,is,NR21,NAX)] = c + I*s;
        EI2[IDX3F(idxzero2-j,ia,is,NR21,NAX)] = c - I*s;
      }
      for(k=1; k<=(NR3+1)/2; k++) {
        c = cos(k*ar3);
        s = sin(k*ar3);
        EI3[IDX3F(idxzero3+k,ia,is,NR31,NAX)] = c + I*s;
        EI3[IDX3F(idxzero3-k,ia,is,NR31,NAX)] = c - I*s;
      }
    }
  }

  // Calculate EIGR
  for(ik=1; ik<=NKPT; ik++) {
    for(is=1; is<=NSP; is++) {
      for(ia=1; ia<=NAX; ia++) {
        for(ig=1; ig<=NGW[ik-1]; ig++) {
          igp = IGK[IDX2F(ig,ik,NGWX)];
          EIGR[IDX4F(ig,ia,is,ik,NGWX,NAX,NSP)] =
            EI1[IDX3F(idxzero1 + IN1[igp-1],ia,is,NR11,NAX)]*
            EI2[IDX3F(idxzero2 + IN2[igp-1],ia,is,NR21,NAX)]*
            EI3[IDX3F(idxzero3 + IN3[igp-1],ia,is,NR31,NAX)];
        }
      }
    }
  }

}
