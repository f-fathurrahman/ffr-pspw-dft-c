// eFeFeR (20910015), November 2011

#include "common_pspw_cuda.h"

void FormFactorAtomic(double complex* c_fft)
{

  int i, is, ir, ig;
  double complex* qat=NULL;
  double *rscal=NULL, *fun=NULL, *fun1=NULL;
  double fs, fp, fd, scalat, scaler, q1, arg;

//
// Allocate memory
//
  qat = (double complex*)malloc(NGX*sizeof(double complex));
  for(i=0; i<NGX; i++) qat[i] = Z_ZERO;
  //
  rscal = (double*)malloc(NSP*sizeof(double));
  for(i=0; i<NSP; i++) rscal[i] = 0.0;
  //
  fun = (double*)malloc(MAXR*sizeof(double));
  for(i=0; i<MAXR; i++) fun[i] = 0.0;
  //
  fun1 = (double*)malloc(MAXR*sizeof(double));
  for(i=0; i<MAXR; i++) fun1[i] = 0.0;

  for(is=1; is<=NSP; is++) {
    fs = 2.0;
    fp = 0.0;
    fd = 0.0;
    scalat = 0.0;
    rscal[is-1] = 0.0;
    //
    if(fabs(ZV[is-1]-1.0) < 1.e-4) {
      fs = 1.0;
    }
    else if(fabs(ZV[is-1]-0.75) < 1.e-4) {
      fs = 0.75;
    }
    else if(fabs(ZV[is-1]-1.25) < 1.e-4 ) {
      fs = 1.25;
    }
    else if(fabs(ZV[is-1]-2.0) < 1.e-4) {
      fs = 2.0;
      fp = 0.0;
    }
    else if(fabs(ZV[is-1]-3.0) < 1.e-4) {
      // hybridisation for gallium
      fs = 2.0;
      fp = 1.0;
    }
    else if(fabs(ZV[is-1]-4.0) < 1.e-4) {
      fs = 2.0;
      fp = 2.0;
    }
    else if(fabs(ZV[is-1]-5.0) < 1.e-4) {
      // hybridisation for arsenic
      fs = 2.0;
      fp = 3.0;
    }
    else if(fabs(ZV[is-1]-6.0) < 1.e-4) {
      fs = 2.0;
      fp = 4.0;
    }
    else if(fabs(ZV[is-1]-7.0) < 1.e-4) {
      fp = 5.0;
    }
    else if(fabs(ZV[is-1]-8.0) < 1.e-4) {
      fp = 6.0;
    }
    else if(fabs(ZV[is-1]-9.0) < 1.e-4) {
      fs = 1.0;
      fd = 8.0;
    }
    else if(fabs(ZV[is-1]-10.0) < 1.e-4) {
      fs = 2.0;
      fd = 8.0;
      printf("** occupation 4d8 5s2  is assumed **\n");
    }
    else if(fabs(ZV[is-1]-11.0) < 1.e-4) {
      fs = 1.0;
      fd = 10.0;
    }
    else if(fabs(ZV[is-1]-12.0) < 1.e-4) {
      fs = 2.0;
      fd = 10.0;
    }
    else if(fabs(ZV[is-1]-13.0) < 1.e-4) {
      fs = 2.0;
      fp = 1.0;
      fd = 10.0;
    }
    else if(fabs(ZV[is-1]-14.0) < 1.e-4) {
      fs = 2.0;
      fp = 2.0;
      fd = 10.0;
    }
    else if(fabs(ZV[is-1]-15.0) < 1.e-4) {
      fs = 2.0;
      fp = 3.0;
      fd = 10.0;
    }
    else if(fabs(ZV[is-1]-16.0) < 1.e-4) {
      fs = 2.0;
      fp = 4.0;
      fd = 10.0;
    }
    else if(fabs(ZV[is-1]-17.0) < 1.e-4) {
      fs = 2.0;
      fp = 5.0;
      fd = 10.0;
    }
    else if(fabs(ZV[is-1]-18.0) < 1.e-4) {
      fs = 2.0;
      fp = 6.0;
      fd = 10.0;
    }
    else {
      printf("Error in FormFactorAtomic not implemented for ZV > 18: %f\n", ZV[is-1]);
      abort();
    }

    // Scale the shape of atomic charge densities
    for(ir=1; ir<=MMAX[is-1]; ir++) {
      if(L_MAX[is-1]==3) {
        fun[ir-1] = ( fs*SQUARE(PSI[IDX3F(ir,is,1,MAXR,NSP)]) + fp*SQUARE(PSI[IDX3F(ir,is,2,MAXR,NSP)]) + fd*SQUARE(PSI[IDX3F(ir,is,3,MAXR,NSP)]) )*R[IDX2F(ir,is,MAXR)]/OMEGA;
      }
      else if(L_MAX[is-1]==2) {
        fun[ir-1] = ( fs*SQUARE(PSI[IDX3F(ir,is,1,MAXR,NSP)]) + fp*SQUARE(PSI[IDX3F(ir,is,2,MAXR,NSP)]) )*R[IDX2F(ir,is,MAXR)]/OMEGA;
      }
      else if(L_MAX[is-1]==1) {
        fun[ir-1] = fs*SQUARE(PSI[IDX3F(ir,is,1,MAXR,NSP)])*R[IDX2F(ir,is,MAXR)]/OMEGA;
      }
      //
      if(scalat > 0.0) {
        scaler = 1.0/(1.0 + exp(scalat*(R[IDX2F(ir,is,MAXR)]-rscal[is-1])));
        fun[ir-1] = fun[ir-1]*scaler;
      }
    } // ir

    // For ig=1
    q1 = simpson(MMAX[is-1], fun, CLOG[is-1]);
    qat[0] = qat[0] + q1*SFAC[IDX2F(is,1,NSP)];
    printf("FormFactorAtomic: rho of atom %d is %f\n", is, q1*OMEGA);

    for(ig=2; ig<=NG; ig++) {
      for(ir=1; ir<=MMAX[is-1]; ir++) {
        arg = R[IDX2F(ir,is,MAXR)]*sqrt(G[ig-1])*TPIBA;
        fun1[ir-1] = fun[ir-1]*sin(arg)/arg;
      }
      q1 = simpson(MMAX[is-1],fun1, CLOG[is-1]);
      qat[ig-1] = qat[ig-1] + q1*SFAC[IDX2F(is,ig,NSP)];
    } // ig
  } // is

  // Copy qat to c_fft (map to 3D grid)
  for(ig=1; ig<=NG; ig++) {
    c_fft[IDX3F(N1[ig-1],N2[ig-1],N3[ig-1],NR1,NR2)] = qat[ig-1];
  }

  free(qat); qat=NULL;
  free(rscal); rscal=NULL;
  free(fun); fun=NULL;
  free(fun1); fun1=NULL;

}
