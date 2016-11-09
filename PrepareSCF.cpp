// eFeFeR (20910015), November 2011

#include "common_pspw.h"

void PrepareSCF()
{
  double complex zs;
  double ds;

  printf("\n");
  printf("Calculating phases factor ....");
  PhaseFactor();
  printf("....Done\n");
  /*SumArray(EI1,(NR1+1)*NAX*NSP,zs); printf("sum(EI1) = (%f,%f)\n", zs.x, zs.y);
  SumArray(EI2,(NR2+1)*NAX*NSP,zs); printf("sum(EI2) = (%f,%f)\n", zs.x, zs.y);
  SumArray(EI3,(NR3+1)*NAX*NSP,zs); printf("sum(EI3) = (%f,%f)\n", zs.x, zs.y);
  SumArray(EIGR,NGWX*NAX*NSP*NKPT,zs); printf("sum(EIGR) = (%f,%f)\n", zs.x, zs.y);*/

  printf("\n");
  printf("Calculating structure factor ....");
  StructureFactor();
  printf("....Done\n");
  /*SumArray(SFAC,NSP*NG,zs); printf("sum(SFAC)= (%f,%f)\n", zs.x, zs.y);*/

  printf("\n");
  printf("Calculating form factor of local pseudopotentials ....");
  FormFactor();
  printf("....Done\n");
  /*printf("ESELF = %f\n", ESELF);
  SumArray(RHOPS,NSP*NGX,ds); printf("SUM(RHOPS) = %f\n", ds);
  SumArray(VPS,NSP*NGX,ds); printf("SUM(VPS) = %f\n", ds);*/

  printf("\n");
  printf("Calculating form factor of non-local pseudopotentials ....");
  NLFormFactor();
  printf("....Done\n");
  /*SumArray(WNL,NSP*NLMAX,ds); printf("sum(WNL) = %f\n", ds);
  SumArray(PKG,NGWX*NSP*NKPT*NLMAX,ds); printf("sum(PKG) = %f\n", ds);
  SumArray(PKG_A,NGWX*NSP*NKPT*NLMAX,ds); printf("sum(PKG_A) = %f\n", ds);*/

}


