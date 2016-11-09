// eFeFeR (20910015), October 2011

#include "common_pspw.h"

void Setup()
{
  // Local variables
  int is, ia, ik;
  int lm;

  // Calculate NLMAX
  NLMAX = 0;
  for(is=1; is<=NSP; is++) {
    lm = L_MAX[is-1]*L_MAX[is-1] - 2*L_LOC[is-1] + 1;
    if(lm > NLMAX) NLMAX = lm;
  }
  if(NLMAX < 4) NLMAX=4;

  // Calculate Fourier mesh
  if(!USER_FFT_GRID) {
    NR1 = 2.0*sqrt(ECUT*(A1[0]*A1[0] + A1[1]*A1[1] + A1[2]*A1[2]))/M_PI;
    NR2 = 2.0*sqrt(ECUT*(A2[0]*A2[0] + A2[1]*A2[1] + A2[2]*A2[2]))/M_PI;
    NR3 = 2.0*sqrt(ECUT*(A3[0]*A3[0] + A3[1]*A3[1] + A3[2]*A3[2]))/M_PI;
    printf("Minimum Fourier grid: (%d,%d,%d)\n", NR1,NR2,NR3);
    // TODO
    NR1 = NR1 + 2;
    NR2 = NR2 + 2;
    NR3 = NR3 + 2;
    // Force the grid size to be even ?
    if(NR1%2 == 1) NR1 = NR1 + 1;
    if(NR2%2 == 1) NR2 = NR2 + 1;
    if(NR3%2 == 1) NR3 = NR3 + 1;
    NNR = NR1*NR2*NR3;
  }
  if(NR1%2 == 1) printf("WARNING: FFT grid is not even NR1=%d\n",NR1);
  if(NR2%2 == 1) printf("WARNING: FFT grid is not even NR2=%d\n",NR2);
  if(NR3%2 == 1) printf("WARNING: FFT grid is not even NR3=%d\n",NR3);

  ALAT = A1[0];

  // Calculate unit cell volume
  OMEGA = determinant3x3(A1,A2,A3);
  //printf("Unit cell volume = %f\n", OMEGA);

  // Calculate reciprocal lattice vectors
  CalcRecipLatt();
  //printf("Reciprocal lattice vectors:\n");
  //PrintVector(B1, 3);
  //PrintVector(B2, 3);
  //PrintVector(B3, 3);

  // Calculate NGX and NGWX
  NGW = (int*)malloc(NKPT*sizeof(int));
  TPIBA = 2.*M_PI/ALAT;
  TPIBA2 = TPIBA*TPIBA;
  GCUT = 4.*ECUT/TPIBA/TPIBA;
  GCUTW = GCUT/4.0;
  //printf("TPIBA = %f\n", TPIBA);
  //printf("TPIBA2 = %f\n", TPIBA2);
  //printf("GCUT = %f\n", GCUT);
  //printf("GCUTW = %f\n", GCUTW);

  GVectorsGen0();
  NGWX = NGX/8 + 8;
  NGX = NGWX*8 + 8;
  NGWX = NGWX + 1;
  //printf("Setup of variables NGWX and NGX:\n");
  //printf("NGWX = %d\n", NGWX);
  //printf("NGX = %d\n", NGX);

  // Read pseudopotentials
  ReadPseudopotentials();

  // Calculate number of electrons
  NEL = 0.0;
  for(is=1; is<=NSP; is++) {
    NEL = NEL + NA[is-1]*ZV[is-1];
  }
  //printf("Number of electrons = %f\n", NEL);

  // Calculate number of states per k-point
  NX = NEMPTY + (int)((NEL + 1.0)/2.0);
  //printf("Number of states = %d\n", NX);

  // Total number of states
  NSTATX = NX*NKPT;
  //printf("Number of states for all k-points = %d\n", NSTATX);

//
// Check input
//
  
  // Number of empty states
  if(TMETAL && (NX/NEMPTY > 50)) {
    printf("Give more empty states for TMETAL:\n");
    printf("Current number of empty states: %d\n", NEMPTY);
    abort();
  }
  
  // Check weight of k-points
  double weight_k_ges = 0.0;
  for(ik=1; ik<=NKPT; ik++) {
    weight_k_ges = weight_k_ges + WKPT[ik-1];
  }
  if(fabs(weight_k_ges-1.0) > 1.e-3) {
    printf("----------------------------------------------------------");
    printf("WARNING: Weight of all k-points : %f\n", weight_k_ges);
    printf("Program will continue because symmetry is not used anyway\n");
    printf("----------------------------------------------------------");
  }
  if(fabs(weight_k_ges-1.0) > 1.e-8) {
    //printf("Rescale k-point weights\n");
    for(ik=1; ik<=NKPT; ik++) WKPT[ik-1] = WKPT[ik-1]/weight_k_ges;
  }

  // Angular momentum component
  for(is=1; is<=NSP; is++) {
    if(L_LOC[is-1] > L_MAX[is-1]) {
      printf("L_LOC=%d must not be larger than L_MAX=%d for is = %d\n", is, L_LOC[is-1], L_MAX[is-1]);
      abort();
    }
  }

  AllocateArrays();

  InitOcc();

  GVectorsGen();
  /*int sumInt; double sumDouble;
  SumArray(G,NGX*3,sumDouble); printf("sum(G) = %f\n", sumDouble);
  SumArray(IN1,NGX,sumInt); printf("sum(IN1) = %d\n", sumInt);
  SumArray(IN2,NGX,sumInt); printf("sum(IN2) = %d\n", sumInt);
  SumArray(IN3,NGX,sumInt); printf("sum(IN3) = %d\n", sumInt);
  SumArray(GG,NGX,sumDouble); printf("sum(GG) = %f\n", sumDouble);
  SumArray(N1,NGX,sumInt); printf("sum(N1) = %d\n", sumInt);
  SumArray(N2,NGX,sumInt); printf("sum(N2) = %d\n", sumInt);
  SumArray(N3,NGX,sumInt); printf("sum(N3) = %d\n", sumInt);*/

  GkVectorsGen();
  /*SumArray(IGK,NGWX*NKPT,sumInt); printf("sum(IGK) = %d\n", sumInt);
  SumArray(XKG,NGWX*NKPT,sumDouble); printf("sum(XKG) = %f\n", sumDouble);
  SumArray(N123,NGWX*NKPT,sumInt); printf("sum(N123) = %d\n", sumInt);
  SumArray(GGK,3*NGWX*NKPT,sumDouble); printf("sum(GGK) = %f\n", sumDouble);*/

}


