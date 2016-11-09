// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void PrintSystemInfo()
{
  int i;
  int ia, is, ik;
  int sumi;

  printf("\n");
  printf("------------------------------------------------\n");
  printf("           System Description:\n");
  printf("------------------------------------------------\n");

  printf("\n");
  printf("Number of species:     %5d\n", NSP);
  for(is=1; is<=NSP; is++) {
    printf("Species name: %s\n", &(ATOMTYP[is-1][0]));
    printf("ZV(%d)      = %14.7f\n", is,ZV[is-1]);
    printf("RGAUSS(%d)  = %14.7f\n", is,RGAUSS[is-1]);
    printf("L_MAX(%d)   = %5d\n", is,L_MAX[is-1]);
    printf("L_LOC(%d)   = %5d\n", is,L_LOC[is-1]);
  }
  
  printf("\n");
  printf("NLMAX:                 %5d\n", NLMAX);
  printf("NAX:                   %5d\n", NAX);

  sumi = 0;
  for(is=1; is<=NSP; is++) sumi = sumi + NA[is-1];
  printf("\nTotal number of atoms: %5d\n", sumi);
  printf("Atomic coordinates (in bohr):\n");
  for(is=1; is<=NSP; is++) {
    for(ia=1; ia<=NA[is-1]; ia++) {
      printf("%s %18.10f %18.10f %18.10f\n", &(ATOMTYP[is-1][0]),TAU[IDX3F(1,ia,is,3,NAX)],
          TAU[IDX3F(2,ia,is,3,NAX)],TAU[IDX3F(3,ia,is,3,NAX)]);
    }
  }

  // Print Lattice Information
  printf("\nLattice vectors (in bohr)\n");
  printf("%18.10f %18.10f %18.10f\n", A1[0], A2[0], A3[0]);
  printf("%18.10f %18.10f %18.10f\n", A1[1], A2[1], A3[1]);
  printf("%18.10f %18.10f %18.10f\n", A1[2], A2[2], A3[2]);

  printf("\nUnit cell volume: %18.10f\n", OMEGA);

  printf("\nReciprocal lattice vectors (in 2*pi/ALAT = %18.10f)\n", 2.0*M_PI/ALAT);
  printf("%18.10f %18.10f %18.10f\n", B1[0], B2[0], B3[0]);
  printf("%18.10f %18.10f %18.10f\n", B1[1], B2[1], B3[1]);
  printf("%18.10f %18.10f %18.10f\n", B1[2], B2[2], B3[2]);

  printf("\nNumber of k-point: %d\n", NKPT);
  for(ik=1; ik<=NKPT; ik++) {
    printf("(%14.10f,%14.10f,%14.10f) %14.10f\n",
        XK[IDX2F(1,ik,3)], XK[IDX2F(2,ik,3)], XK[IDX2F(3,ik,3)], WKPT[ik-1]);
  }

  printf("\n");
  printf("Number of electrons (NEL):         %6d\n", (int)NEL);
  printf("Number of occupied states:         %6d\n", (int)(NEL+1)/2);
  printf("Number of empty states  (NEMPTY):  %6d\n", NEMPTY);
  printf("Number of states per k-point (NX): %6d\n", NX);
  printf("Total number of states (NSTATX):   %6d\n", NSTATX);

  printf("\nCut off energy (ECUT): %18.10f Ry\n", ECUT);
  printf("Fourier grid: (%8d,%8d,%8d)\n", NR1,NR2,NR3);
  printf("NGWX     = %8d\n", NGWX);
  printf("NGX      = %8d\n", NGX);
  printf("NG       = %8d\n", NG);
  printf("GCUT     = %18.10f\n", GCUT);
  printf("GCUTW    = %18.10f\n", GCUTW);
  printf("TPIBA    = %18.10f\n", TPIBA);
  printf("TPIBA2   = %18.10f\n", TPIBA2);
  for(ik=1; ik<=NKPT; ik++) {
    printf("NGW(%3d) = %8d\n", ik, NGW[ik-1]);
  }
}

