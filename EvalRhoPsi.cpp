// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void EvalRhoPsi(int ik, double *rhoe, double complex *eigvec, int NGW_ik)
{
  int ig, ist, i1,i2,i3, igp;
  double xkin, wdotf, wdotfbo, wdot2fb;
  double ss, sk, scg, tt;
  double complex c_t, cc;
  double complex *c_fft=NULL;

  c_fft = (double complex*)malloc(NNR*sizeof(double complex));

  xkin = 0.0;

  // Calculate new wavevectors C0 in planewave basis
  // and also kinetic energy and charge density
  for(ist=1; ist<=NX; ist++) {
    wdotf = WKPT[ik-1]*FOCC[IDX2F(ist,ik,NX)];
    wdotfbo = wdotf/OMEGA;
    wdot2fb = wdotf*0.5;
    ss = 0.0;
    sk = 0.0;

    // Zero out FFT grid values
    for(int i=0; i<NNR; i++) c_fft[i] = Z_ZERO;

    // Zero out C0(:,ist,ik)
    for(ig=1; ig<=NGWX; ig++) {
      C0[IDX3F(ig,ist,ik,NGWX,NX)] = Z_ZERO;
    }

    // Assign new eigenvectors
    for(ig=1; ig<=NGWX; ig++) {
      C0[IDX3F(ig,ist,ik,NGWX,NX)] = eigvec[IDX2F(ig,ist,NGWX)];
    }

    for(ig=1; ig<NGW_ik; ig++) {
      c_t = C0[IDX3F(ig,ist,ik,NGWX,NX)];
      igp = IGK[IDX2F(ig,ik,NGWX)];
      // TODO: Any other ideas?
      c_fft[IDX3F(N1[igp-1],N2[igp-1],N3[igp-1],NR1,NR2)] = c_t;
      //scg = SQUARE(c_t.x) + SQUARE(c_t.y);
      scg = SQUARE( cabs( c_t ) );
      ss = ss + scg;
      sk = sk + XKG[IDX2F(ig,ik,NGWX)]*scg;
    }

    // Check charge sum
    if(fabs(ss-1.0) > 6.e-3) {
      printf("Wrong charge density in G-space: %18.10f\n", ss);
      printf("for k-point %d and states %d\n", ik, ist);
      abort();
    }

    // Add kinetic energy of this band to total kinetic energy
    xkin = xkin + sk*wdot2fb;

    // Inverse FFT (transform to real space)
    fft_fftw3(c_fft,NR1,NR2,NR3,true);

    for(i3=1; i3<=NR3; i3++) {
      for(i2=1; i2<=NR2; i2++) {
        for(i1=1; i1<=NR1; i1++) {
          cc = c_fft[IDX3F(i1,i2,i3,NR1,NR2)];
          //tt = ( SQUARE(cc.x) + SQUARE(cc.y) )*wdotfbo;
          tt = SQUARE(cabs(cc)) * wdotfbo;
          rhoe[IDX3F(i1,i2,i3,NR1,NR2)] = rhoe[IDX3F(i1,i2,i3,NR1,NR2)] + tt;
        }
      }
    }
  } // end of loop over states

  // Calculate kinetic energy
  EKIN = EKIN + xkin*TPIBA2;

  free(c_fft); c_fft=NULL;
}


//
void EvalRhoPsi(int ik, double *rhoe, double complex *eigvec, int ldEigvec, int NGW_ik)
{
  int ig, ist, i1,i2,i3, igp;
  double xkin, wdotf, wdotfbo, wdot2fb;
  double ss, sk, scg, tt;
  double complex c_t, cc;
  double complex *c_fft=NULL;

  c_fft = (double complex*)malloc(NNR*sizeof(double complex));

  xkin = 0.0;

  // Calculate new wavevectors C0 in planewave basis
  // and also kinetic energy and charge density
  for(ist=1; ist<=NX; ist++) {
    wdotf = WKPT[ik-1]*FOCC[IDX2F(ist,ik,NX)];
    wdotfbo = wdotf/OMEGA;
    wdot2fb = wdotf*0.5;
    ss = 0.0;
    sk = 0.0;

    // Zero out FFT grid values
    for(int i=0; i<NNR; i++) c_fft[i] = Z_ZERO;

    // Zero out C0(:,ist,ik)
    for(ig=1; ig<=NGWX; ig++) {
      C0[IDX3F(ig,ist,ik,NGWX,NX)] = Z_ZERO;
    }

    // Assign new eigenvectors
    for(ig=1; ig<=NGWX; ig++) {
      C0[IDX3F(ig,ist,ik,NGWX,NX)] = eigvec[IDX2F(ig,ist,ldEigvec)];
    }

    for(ig=1; ig<NGW_ik; ig++) {
      c_t = C0[IDX3F(ig,ist,ik,NGWX,NX)];
      igp = IGK[IDX2F(ig,ik,NGWX)];
      // TODO: Any other ideas?
      c_fft[IDX3F(N1[igp-1],N2[igp-1],N3[igp-1],NR1,NR2)] = c_t;
      //scg = SQUARE(c_t.x) + SQUARE(c_t.y);
      scg = SQUARE(cabs(c_t));
      ss = ss + scg;
      sk = sk + XKG[IDX2F(ig,ik,NGWX)]*scg;
    }

    // Check charge sum
    if(fabs(ss-1.0) > 6.e-3) {
      printf("Wrong charge density in G-space: %18.10f\n", ss);
      printf("for k-point %d and states %d\n", ik, ist);
      abort();
    }

    // Add kinetic energy of this band to total kinetic energy
    xkin = xkin + sk*wdot2fb;

    // Inverse FFT (transform to real space)
    fft_fftw3(c_fft,NR1,NR2,NR3,true);

    for(i3=1; i3<=NR3; i3++) {
      for(i2=1; i2<=NR2; i2++) {
        for(i1=1; i1<=NR1; i1++) {
          cc = c_fft[IDX3F(i1,i2,i3,NR1,NR2)];
          //tt = ( SQUARE(cc.x) + SQUARE(cc.y) )*wdotfbo;
          tt = SQUARE(cabs(cc)) * wdotfbo;
          rhoe[IDX3F(i1,i2,i3,NR1,NR2)] = rhoe[IDX3F(i1,i2,i3,NR1,NR2)] + tt;
        }
      }
    }
  } // end of loop over states

  // Calculate kinetic energy
  EKIN = EKIN + xkin*TPIBA2;

  free(c_fft); c_fft=NULL;
}


