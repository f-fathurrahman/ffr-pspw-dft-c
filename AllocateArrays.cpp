// eFeFeR (20910015), October 2011

#include "common_pspw_cuda.h"

void AllocateArrays()
{
  // Local
  int i;

  // Occupation numbers
  FOCC = (double*)malloc(NX*NKPT*sizeof(double)); for(i=0; i<NX*NKPT; i++) FOCC[i] = 0.0;

  // G vectors
  G = (double*)malloc(NGX*sizeof(double)); for(i=0; i<NGX; i++) G[i] = 0.0;
  GG = (double*)malloc((NGX+1)*3*sizeof(double)); for(i=0; i<(NGX+1)*3; i++) GG[i] = 0.0;
  N1 = (int*)malloc(NGX*sizeof(int)); for(i=0; i<NGX; i++) N1[i] = 0;
  N2 = (int*)malloc(NGX*sizeof(int)); for(i=0; i<NGX; i++) N2[i] = 0;
  N3 = (int*)malloc(NGX*sizeof(int)); for(i=0; i<NGX; i++) N3[i] = 0;
  IN1 = (int*)malloc(NGX*sizeof(int)); for(i=0; i<NGX; i++) IN1[i] = 0;
  IN2 = (int*)malloc(NGX*sizeof(int)); for(i=0; i<NGX; i++) IN2[i] = 0;
  IN3 = (int*)malloc(NGX*sizeof(int)); for(i=0; i<NGX; i++) IN3[i] = 0;

  // (G+k) vectors
  IGK = (int*)malloc(NGWX*NKPT*sizeof(int)); for(i=0; i<NGWX*NKPT; i++) IGK[i] = 0;
  XKG = (double*)malloc(NGWX*NKPT*sizeof(double)); for(i=0; i<NGWX*NKPT; i++) XKG[i] = 0.0;
  N123 = (unsigned int*)malloc(NGWX*NKPT*sizeof(unsigned int)); for(i=0; i<NGWX*NKPT; i++) N123[i] = 0;
  GGK = (double*)malloc(3*NGWX*NKPT*sizeof(double)); for(i=0; i<3*NGWX*NKPT; i++) GGK[i] = 0.0;

  EI1 = (double complex*)malloc( (NR1+1)*NAX*NSP*sizeof(double complex) );
  for(i=0; i<(NR1+1)*NAX*NSP; i++) EI1[i] = Z_ZERO;
  EI2 = (double complex*)malloc( (NR2+1)*NAX*NSP*sizeof(double complex) );
  for(i=0; i<(NR2+1)*NAX*NSP; i++) EI2[i] = Z_ZERO;
  EI3 = (double complex*)malloc( (NR3+1)*NAX*NSP*sizeof(double complex) );
  for(i=0; i<(NR3+1)*NAX*NSP; i++) EI3[i] = Z_ZERO;
  EIGR = (double complex*)malloc( NGWX*NAX*NSP*NKPT*sizeof(double complex) );
  for(i=0; i<NGWX*NAX*NSP*NKPT; i++) EIGR[i] = Z_ZERO;

  // Structure factor
  SFAC = (double complex*)malloc( NSP*NGX*sizeof(double complex) );
  for(i=0; i<NSP*NGX; i++) SFAC[i] = Z_ZERO;

  RHOPS = (double*)malloc( NSP*NGX*sizeof(double) );
  for(i=0; i<NSP*NGX; i++) RHOPS[i] = 0.0;
  VPS = (double*)malloc( NSP*NGX*sizeof(double) );
  for(i=0; i<NSP*NGX; i++) VPS[i] = 0.0;

  //
  PKG = (double*)malloc( NGWX*NSP*NKPT*NLMAX*sizeof(double) );
  for(i=0; i<NGWX*NSP*NKPT*NLMAX; i++) PKG[i] = 0.0;
  PKG_A = (double*)malloc( NLMAX*NGWX*NSP*NKPT*sizeof(double) );
  for(i=0; i<NLMAX*NGWX*NSP*NKPT; i++) PKG_A[i] = 0.0;
  WNL = (double*)malloc( NSP*NLMAX*sizeof(double) );
  for(i=0; i<NSP*NLMAX; i++) WNL[i] = 0.0;

  //
  EIG = (double*)malloc(NX*NKPT*sizeof(double)); for(i=0; i<NX*NKPT; i++) EIG[i] = 0.0; 

  //
  FNL = (double complex*)malloc(NX*NKPT*NSP*NAX*NLMAX*sizeof(double complex));
  for(i=0; i<NX*NKPT*NSP*NAX*NLMAX; i++) FNL[i] = Z_ZERO;

  //
  C0 = (double complex*)malloc(NGWX*NX*NKPT*sizeof(double complex));
  for(i=0; i<NGWX*NX*NKPT; i++) C0[i] = Z_ZERO;

  CPUMEM += NX*NKPT*sizeof(double);
  CPUMEM += NGX*sizeof(double) + (NGX+1)*3*sizeof(double) + 6*NGX*sizeof(int);
  CPUMEM += 2*NGWX*NKPT*sizeof(int) + 4*NGWX*NKPT*sizeof(double);
  CPUMEM += (NR1 + NR2 + NR3 + 3)*NAX*NSP*sizeof(double complex);
  CPUMEM += NSP*NGX*sizeof(double complex);
  CPUMEM += 2*NSP*NGX*sizeof(double);
  CPUMEM += NGWX*NAX*NSP*NKPT*sizeof(double complex);
  CPUMEM += (2*NGWX*NSP*NKPT*NLMAX + NSP*NLMAX)*sizeof(double);
  CPUMEM += NX*NKPT*sizeof(double);
  CPUMEM += NGWX*NX*NKPT*sizeof(double complex);

  //printf("\n*****************************************\n");
  //printf("Allocated CPU memory for arrays: %f MB\n", CPUMEM/1024./1024.);
  //printf("*****************************************\n");

  //printf("EXIT: %s\n", __FILE__);
}

