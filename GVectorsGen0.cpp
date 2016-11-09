// eFeFeR (209100150), October 2011

#include "common_pspw.h"

void GVectorsGen0()
{
  // Local variables
  int ik;
  int i, j, k;
  int jg;
  int n1_start, n2_start;
  double g1tmp, g2tmp, g3tmp, gtmp;
  double xtmp;

  for(ik=1; ik<=NKPT; ik++) {
    NGW[ik-1] = 1;
  }

  jg = 2;
  //
  for(k=0; k<=NR3/2; k++) {
    if(k==0) {
      n2_start = 0;
    } else {
      n2_start = -NR2/2;
    }
    //
    for(j=n2_start; j<=NR2/2; j++) {
      if(j==0 && k==0) {
        n1_start = 1;
      } else {
        n1_start = -NR1/2;
      }
      //
      for(i=n1_start; i<=NR1/2; i++) {
        g1tmp = i*B1[0] + j*B2[0] + k*B3[0];
        g2tmp = i*B1[1] + j*B2[1] + k*B3[1];
        g3tmp = i*B1[2] + j*B2[2] + k*B3[2];
        gtmp = g1tmp*g1tmp + g2tmp*g2tmp + g3tmp*g3tmp;
        //printf("%f %f %f %f\n", g1tmp, g2tmp, g3tmp, gtmp);
        for(ik=1; ik<=NKPT; ik++) {
          xtmp = SQUARE(g1tmp - XK[IDX2F(1,ik,3)]) + SQUARE(g2tmp - XK[IDX2F(2,ik,3)]) + SQUARE(g3tmp - XK[IDX2F(3,ik,3)]);
          //printf("ik = %d xtmp = %f\n", ik, xtmp);
          if(xtmp <= GCUTW) NGW[ik-1] = NGW[ik-1] + 1;
          xtmp = SQUARE(g1tmp + XK[IDX2F(1,ik,3)]) + SQUARE(g2tmp + XK[IDX2F(2,ik,3)]) + SQUARE(g3tmp + XK[IDX2F(3,ik,3)]);
          //printf("ik = %d xtmp = %f\n", ik, xtmp);
          if(xtmp <= GCUTW) NGW[ik-1] = NGW[ik-1] + 1;
        }
        if(gtmp < GCUT) {
          jg = jg + 1;
        }
      } // i
    } // j
  }// k

  NGX = 2*(jg-1) + 1;
  //printf("NGX = %d\n", NGX);

  for(ik=1; ik<=NKPT; ik++) {
    //printf("NGW(%d) = %d\n",ik,NGW[ik-1]);
    if(NGX < 8*NGW[ik-1]) NGX = 8*NGW[ik-1];
  }

}

