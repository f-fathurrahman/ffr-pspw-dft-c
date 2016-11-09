// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void ApplyHam_block(double complex *vin, double complex *vout, int NGW_ik, int ik, double *rho, int nnr,
    int nstates)
{
  // Local
  double complex *c_fft=NULL;
  double complex *t_eigr_pkg=NULL;
  double rr1, rr2, rr3, rr4, sumd;
  double complex tt1, tt2, tt3, tt4;
  double rr[NLMAX];
  int i, ig, is, ia, il, lm_end, ist;
  double complex sumz;
  double complex tz;

  if(nnr != NNR) {
    printf("nnr != NNR : %d %d\n", nnr,NNR);
    abort();
  }

  // Allocate memory
  c_fft = (double complex*)malloc(NNR*sizeof(double complex));
  t_eigr_pkg = (double complex*)malloc(NGWX*NLMAX*sizeof(double complex));

  // Begin loop over states
  for(ist=1; ist<=nstates; ist++) {

//
// Local part
//
  for(i=0; i<NNR; i++) c_fft[i] = Z_ZERO;
  for(ig=1; ig<=NGW_ik; ig++) {
    c_fft[N123[IDX2F(ig,ik,NGWX)]-1] = vin[IDX2F(ig,ist,NGW_ik)];
  }
  /*zdotc_(&tz, &NNR, c_fft,&I_ONE, c_fft,&I_ONE);
  printf("After set c_fft: tz = (%18.10f,%18.10f)\n", tz.x, tz.y);*/

  // Inverse FFT
  fft_fftw3(c_fft, NR1,NR2,NR3, true);
  /*zdotc_(&tz, &NNR, c_fft,&I_ONE, c_fft,&I_ONE);
  printf("After inverse FFT: tz = (%18.10f,%18.10f)\n", tz.x, tz.y);*/

  // Multiply local potential and wavefunction in real space
  for(i=1; i<=NNR; i++) {
    c_fft[i-1] = rho[i-1]*c_fft[i-1];
  }
  /*zdotc_(&tz, &NNR, c_fft,&I_ONE, c_fft,&I_ONE);
  printf("After multiply local potential: tz = (%18.10f,%18.10f)\n", tz.x, tz.y);*/

  // Transform back to momentum space
  fft_fftw3(c_fft, NR1,NR2,NR3, false);
  /*zdotc_(&tz, &NNR, c_fft,&I_ONE, c_fft,&I_ONE);
  printf("After forward FFT: tz = (%18.10f,%18.10f)\n", tz.x, tz.y);*/

  // Add kinetic contribution
  for(ig=1; ig<=NGW_ik; ig++) {
    vout[IDX2F(ig,ist,NGW_ik)] = 0.5*XKG[IDX2F(ig,ik,NGWX)]*TPIBA2*vin[IDX2F(ig,ist,NGW_ik)] +
      c_fft[N123[IDX2F(ig,ik,NGWX)]-1];
  }
  /*zdotc_(&tz, &NGW_ik, &vout[IDX2F(1,ist,NGW_ik)],&I_ONE, &vout[IDX2F(1,ist,NGW_ik)],&I_ONE);
  printf("After add kinetic contribution: tz = (%18.10f,%18.10f)\n", tz.x, tz.y);*/

//
// Apply non-local part
//
  for(is=1; is<=NSP; is++) {
    lm_end = L_MAX[is-1]*L_MAX[is-1] - 2*L_LOC[is-1] + 1;
    if(lm_end > 4) {
      printf("ERROR: lm_end > 4 is not implemented\n");
      abort();
    }
    if(lm_end == 4) {
      rr1 = WNL[IDX2F(is,1,NSP)];
      rr2 = WNL[IDX2F(is,2,NSP)];
      rr3 = WNL[IDX2F(is,3,NSP)];
      rr4 = WNL[IDX2F(is,4,NSP)];
    } else {
      for(il=1; il<=lm_end; il++) {
        rr[il-1] = WNL[IDX2F(is,il,NSP)];
      }
    }

    //
    for(ia=1; ia<=NA[is-1]; ia++) {
      //
      for(il=1; il<=lm_end; il++) {
        for(ig=1; ig<=NGW_ik; ig++) {
          t_eigr_pkg[IDX2F(ig,il,NGWX)] = creal(EIGR[IDX4F(ig,ia,is,ik,NGWX,NAX,NSP)]) *
            PKG_A[IDX4F(il,ig,is,ik,NLMAX,NGWX,NSP)] +
            I*cimag(EIGR[IDX4F(ig,ia,is,ik,NGWX,NAX,NSP)]) *
            -PKG_A[IDX4F(il,ig,is,ik,NLMAX,NGWX,NSP)];
        }
      }
      //
      if(lm_end==4) {
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,1,NGWX)],&I_ONE);
        tt1 = rr1*tz;
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,2,NGWX)],&I_ONE);
        tt2 = rr2*tz;
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,3,NGWX)],&I_ONE);
        tt3 = rr3*tz;
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,4,NGWX)],&I_ONE);
        tt4 = rr4*tz;
        for(ig=1; ig<=NGW_ik; ig++) {
          vout[IDX2F(ig,ist,NGW_ik)] = vout[IDX2F(ig,ist,NGW_ik)] +
            EIGR[IDX4F(ig,ia,is,ik,NGWX,NAX,NSP)]*(
              PKG_A[IDX4F(1,ig,is,ik,NLMAX,NGWX,NSP)]*tt1 +
              PKG_A[IDX4F(2,ig,is,ik,NLMAX,NGWX,NSP)]*tt2 +
              PKG_A[IDX4F(3,ig,is,ik,NLMAX,NGWX,NSP)]*tt3 +
              PKG_A[IDX4F(4,ig,is,ik,NLMAX,NGWX,NSP)]*tt4 );
        }
      }
      //
      else if(lm_end==3) {
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,1,NGWX)],&I_ONE);
        tt1 = rr[0]*tz;
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,2,NGWX)],&I_ONE);
        tt2 = rr[1]*tz;
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,3,NGWX)],&I_ONE);
        tt3 = rr[2]*tz;
        for(ig=1; ig<=NGW_ik; ig++) {
          vout[IDX2F(ig,ist,NGW_ik)] = vout[IDX2F(ig,ist,NGW_ik)] +
            EIGR[IDX4F(ig,ia,is,ik,NGWX,NAX,NSP)]*(
              PKG_A[IDX4F(1,ig,is,ik,NLMAX,NGWX,NSP)]*tt1 +
              PKG_A[IDX4F(2,ig,is,ik,NLMAX,NGWX,NSP)]*tt2 +
              PKG_A[IDX4F(3,ig,is,ik,NLMAX,NGWX,NSP)]*tt3 );
        }
      }
      //
      else if(lm_end==2) {
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,1,NGWX)],&I_ONE);
        tt1 = rr[0]*tz;
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,2,NGWX)],&I_ONE);
        tt2 = rr[1]*tz;
        for(ig=1; ig<=NGW_ik; ig++) {
          vout[IDX2F(ig,ist,NGW_ik)] = vout[IDX2F(ig,ist,NGW_ik)] +
            EIGR[IDX4F(ig,ia,is,ik,NGWX,NAX,NSP)]*(
              PKG_A[IDX4F(1,ig,is,ik,NLMAX,NGWX,NSP)]*tt1 +
              PKG_A[IDX4F(2,ig,is,ik,NLMAX,NGWX,NSP)]*tt2 );
        }
      }
      //
      else if(lm_end==1) {
        zdotu_(&tz, &NGW_ik, &vin[IDX2F(1,ist,NGW_ik)],&I_ONE,
            &t_eigr_pkg[IDX2F(1,1,NGWX)],&I_ONE);
        tt1 = rr[0]*tz;
        for(ig=1; ig<=NGW_ik; ig++) {
          vout[IDX2F(ig,ist,NGW_ik)] = vout[IDX2F(ig,ist,NGW_ik)] +
            EIGR[IDX4F(ig,ia,is,ik,NGWX,NAX,NSP)]*(
              PKG_A[IDX4F(1,ig,is,ik,NLMAX,NGWX,NSP)]*tt1 );
        }
      }
    } // ia
  } // is

  } // end of loop over states

  // Free memory
  free(c_fft); c_fft=NULL;
  free(t_eigr_pkg); t_eigr_pkg=NULL;
}
