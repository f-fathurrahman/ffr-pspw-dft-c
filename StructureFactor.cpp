// eFeFeR (20910015), November 2011

#include "common_pspw_cuda.h"

void StructureFactor()
{
  int is, ig, ia;
  double complex sf;
  int idxzero1 = (NR1+1)/2 + 1;
  int idxzero2 = (NR2+1)/2 + 1;
  int idxzero3 = (NR3+1)/2 + 1;
  int NR11 = NR1 + 1;
  int NR21 = NR2 + 1;
  int NR31 = NR3 + 1;

  for(is=1; is<=NSP; is++) {
    for(ig=1; ig<=NG; ig++) {
      sf = Z_ZERO;
      for(ia=1; ia<=NA[is-1]; ia++) {
        sf = sf + EI1[IDX3F(idxzero1 + IN1[ig-1],ia,is,NR11,NAX)]*
            EI2[IDX3F(idxzero2 + IN2[ig-1],ia,is,NR21,NAX)]*
            EI3[IDX3F(idxzero3 + IN3[ig-1],ia,is,NR31,NAX)];
      } // ia
      SFAC[IDX2F(is,ig,NSP)] = sf;
    } // ig
  } // is

  //printf("EXIT: %s\n", __FILE__);
}

