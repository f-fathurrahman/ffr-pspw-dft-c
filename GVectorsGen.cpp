// eFeFeR (20910015), October 2011

#include "common_pspw.h"

void GVectorsGen()
{
  // Local
  int i, j, k;
  int *index_tmp=NULL, *idum=NULL;
  double *dum;
  double gtmp, g1tmp, g2tmp, g3tmp;
  int NGX21=NGX/2+1;
  int ng_stop, jg, n2_start, n1_start;

  index_tmp = (int*)malloc(NGX21*sizeof(int)); for(i=0; i<NGX21; i++) index_tmp[i] = 0;
  idum = (int*)malloc(NGX21*3*sizeof(int)); for(i=0; i<NGX21*3; i++) idum[i] = 0;
  dum = (double*)malloc(NGX21*sizeof(double)); for(i=0; i<NGX21; i++) dum[i] = 0.0;

// Generate G vectors of the first half
  jg = 2;

  for(k=0; k<=NR3/2; k++) {
    if(k==0) {
      n2_start = 0;
    } else {
      n2_start = -NR2/2;
    }
    for(j=n2_start; j<=NR2/2; j++) {
      if(j==0 && k==0) {
        n1_start = 1;
      } else {
        n1_start = -NR1/2;
      }
      for(i=n1_start; i<=NR1/2; i++) {
        g1tmp = i*B1[0] + j*B2[0] + k*B3[0];
        g2tmp = i*B1[1] + j*B2[1] + k*B3[1];
        g3tmp = i*B1[2] + j*B2[2] + k*B3[2];
        gtmp = g1tmp*g1tmp + g2tmp*g2tmp + g3tmp*g3tmp;
        if(gtmp <= GCUT) {
          if(jg > NGX21) { printf("ERROR: jg > NGX/2+1\n"); abort(); }
          dum[jg-1] = gtmp;
          idum[IDX2F(jg,1,NGX21)] = i;
          idum[IDX2F(jg,2,NGX21)] = j;
          idum[IDX2F(jg,3,NGX21)] = k;
          jg = jg + 1;
        }
      } // i
    } // j
  } // k

  ng_stop = jg - 1;
  
  // Sort G vectors in order of increasing magnitude
  sort1(ng_stop,dum,index_tmp);
  //SortGVectors(ng_stop,dum,index_tmp);
  //driver_orderf(ng_stop,dum,index_tmp);

  // The first G vector must be zero vector
  if(dum[index_tmp[0]-1] != 0) {
    printf("ERROR in GVectorsGen: first G != 0\n");
    abort();
  }
  // Index of the zero vector must be 1
  if(index_tmp[0] != 1) {
    printf("ERROR in GVectorsGen: first index != 1\n");
    abort();
  }
  // Fill the actual G vectors array
  // The first vector is zero vector
  G[0] = 0.0;
  IN1[0] = 0; IN2[0] = 0; IN3[0] = 0;
  GG[IDX2F(1,1,NGX+1)] = 0.0;
  GG[IDX2F(1,2,NGX+1)] = 0.0;
  GG[IDX2F(1,3,NGX+1)] = 0.0;
  N1[0] = 1; N2[0] = 1; N3[0] = 1;
  
  // Extend G vector to other half
  NG = 2; // FIXME: Due to difference between sort
  for(jg=2; jg<=ng_stop; jg++) {
    // Check wheter current number of G vectors is allowed
    if(NG+1 > NGX) {
      printf("ERROR in GVectorsGen: NG+1=%d > NGX=%d\n", NG+1, NGX);
      abort();
    }
    i = idum[IDX2F(index_tmp[jg-1],1,NGX21)];
    j = idum[IDX2F(index_tmp[jg-1],2,NGX21)];
    k = idum[IDX2F(index_tmp[jg-1],3,NGX21)];
    gtmp = dum[index_tmp[jg-1]-1];
    // Indices of G vectors: first half and other half
    IN1[NG-1] = i;
    IN2[NG-1] = j;
    IN3[NG-1] = k;
    IN1[NG] = -i;
    IN2[NG] = -j;
    IN3[NG] = -k;
    // The G vectors themself: first half and other half
    GG[IDX2F(NG,1,NGX+1)] = i*B1[0] + j*B2[0] + k*B3[0];
    GG[IDX2F(NG,2,NGX+1)] = i*B1[1] + j*B2[1] + k*B3[1];
    GG[IDX2F(NG,3,NGX+1)] = i*B1[2] + j*B2[2] + k*B3[2];
    GG[IDX2F(NG+1,1,NGX+1)] = -GG[IDX2F(NG,1,NGX+1)];
    GG[IDX2F(NG+1,2,NGX+1)] = -GG[IDX2F(NG,2,NGX+1)];
    GG[IDX2F(NG+1,3,NGX+1)] = -GG[IDX2F(NG,3,NGX+1)];
    // The magnitude of G vectors: first half and other half
    G[NG-1] = gtmp;
    G[NG] = gtmp;
    // Indices: first half and other half
    N1[NG-1] = i + 1;
    if(i < 0) N1[NG-1] = N1[NG-1] + NR1;
    N2[NG-1] = j + 1;
    if(j < 0) N2[NG-1] = N2[NG-1] + NR2;
    N3[NG-1] = k + 1;
    if(k < 0) N3[NG-1] = N3[NG-1] + NR3;
    i = -i;
    j = -j;
    k = -k;
    N1[NG] = i + 1;
    if(i < 0) N1[NG] = N1[NG] + NR1;
    N2[NG] = j + 1;
    if(j < 0) N2[NG] = N2[NG] + NR2;
    N3[NG] = k + 1;
    if(k < 0) N3[NG] = N3[NG] + NR3;
    NG = NG + 2;
  }

  NG = NG - 1 ;

  free(index_tmp); index_tmp=NULL;
  free(idum); idum=NULL;
  free(dum); dum=NULL;

  /*printf("Several generated G vectors\n");
  for(i=NG; i>=NG-20; i--) {
    printf("%d (%d,%d,%d) (%f,%f,%f) %f\n", i, IN1[i-1],IN2[i-1],IN3[i-1],
        GG[IDX2F(i,1,NGX+1)],GG[IDX2F(i,2,NGX+1)],GG[IDX2F(i,3,NGX+1)],G[i-1]);
  }*/

}
