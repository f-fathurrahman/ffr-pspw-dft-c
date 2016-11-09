// eFeFeR (20910015), December 2011

#include "common_pspw.h"

void test_ApplyHam(double *rhoe)
{
  double complex *vin=NULL;
  double complex *vout=NULL;
  int ik=1;
  int i;
  double sumd;
  double complex sumz;

  vin = (double complex*)malloc(NGW[ik-1]*sizeof(double complex));
  vout = (double complex*)malloc(NGW[ik-1]*sizeof(double complex));

  SumArray(rhoe,NNR,sumd);
  printf("sum(rhoe) = %18.10f\n", sumd);

  for(i=0; i<NGW[ik-1]; i++) vin[i] = make_double complex(1.0,2.0);

  ApplyHam_block(vin,vout,NGW[ik-1],ik,rhoe,NNR,1);

  free(vin); vin=NULL;
  free(vout); vout=NULL;
}

