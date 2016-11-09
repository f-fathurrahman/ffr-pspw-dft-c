// eFeFeR (20910015), October 2011

#include "common_pspw.h"

void ReadPseudopotentials()
{
  // Local variables
  int is, i, il, ir;
  int nang;
  string line;

  // Allocate memory
  ZV = (double*)malloc(NSP*sizeof(double));
  MMAX = (int*)malloc(NSP*sizeof(int));
  CLOG = (double*)malloc(NSP*sizeof(double));
  R = (double*)malloc(MAXR*NSP*sizeof(double));
  PSI = (double*)malloc(MAXR*NSP*MAXAM*sizeof(double));
  VION = (double*)malloc(MAXR*NSP*MAXAM*sizeof(double));

  for(is=1; is<=NSP; is++) {
    //printf("Reading pseudopotential for species: %s\n", &ATOMTYP[is-1][0]);
    ifstream pspFile(&ATOMTYP[is-1][0]);
    if(!pspFile.is_open()) {
      printf("ERROR: Cannot open pseudopotential for species: %s\n", &ATOMTYP[is-1][0]);
      abort();
    }
    // Read ion charge and number of angular momentum
    getline(pspFile,line);
    sscanf(&line[0], "%lf %d", &ZV[is-1], &nang);
    //printf("ZV(%d) = %f nang = %d\n", is, ZV[is-1], nang);

    if(nang < L_MAX[is-1]) {
      printf("ERROR: Pseudopotential does not contain enough angular momentum\n");
      printf("for L_MAX(%d) = %d\n", is, L_MAX[is-1]);
      printf("Provided in pseudopotential only up to %d\n", nang);
      abort();
    }

    // Skip reading Gaussian parameters for pseudopotentials
    for(i=1; i<=10; i++) {
      getline(pspFile,line);
    }

    if(nang > L_MAX[is-1]) nang = L_MAX[is-1];
    // Read atomic wavefunctions and pseudopotential on radial mesh
    for(il=1; il<=nang; il++) {
      //printf("il = %d\n", il);
      getline(pspFile,line);
      sscanf(&line[0], "%d %lf", &MMAX[is-1], &CLOG[is-1]);
      //printf("%d %f\n", MMAX[is-1], CLOG[is-1]);
      if(MMAX[is-1] > MAXR) {
        printf("ERROR: MMAX(%d)=%d > MAXR=%d\n", is, MMAX[is-1], MAXR);
        printf("Please increase value of MAXR in common_pspw.h\n");
        abort();
      }
      CLOG[is-1] = log(CLOG[is-1]);
      for(ir=1; ir<=MMAX[is-1]; ir++) {
        getline(pspFile,line);
        sscanf(&line[0], "%d %lf %lf %lf", &i, &R[IDX2F(ir,is,MAXR)], &PSI[IDX3F(ir,is,il,MAXR,NSP)], &VION[IDX3F(ir,is,il,MAXR,NSP)]);
        //printf("%d %f %f %f\n", i, R[IDX2F(ir,is,MAXR)], PSI[IDX3F(ir,is,il,MAXR,NSP)], VION[IDX3F(ir,is,il,MAXR,NSP)]);
      }
    } // il
    // TODO: Handle nonlinear core-correction
    pspFile.close();
  }
}

