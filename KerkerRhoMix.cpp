// eFeFeR (20910015), January 2012

#include "common_pspw.h"

void KerkerRhoMix(double beta, double q0, double *rhoe, double *rho_old, int nnr, double &d)
{
  int ik=1; // TODO: consider only the first k-point in list?
  double complex *c_fft=NULL;
  int i,ig;
  double gg;

  printf("Kerker rho mixing with beta = %f\n", beta);

  c_fft = (double complex*)malloc(nnr*sizeof(double complex));

  for(i=0; i<nnr; i++) {
    c_fft[i] = rhoe[i] - rho_old[i] + I*0.0;
  }

  // Transform to reciprocal space (forward FFT)
  fft_fftw3(c_fft,NR1,NR2,NR3,false);

  // Filter
  for(ig=1; ig<=NGW[ik-1]; ig++) {
    gg = XKG[IDX2F(ig,ik,NGWX)];
    /*c_fft[N123[IDX2F(ig,ik,NGWX)]-1] = (0.8*gg/(gg+0.5) - 0.8)*c_fft[N123[IDX2F(ig,ik,NGWX)]-1];
    if(gg<=1.e-10) { // small G
      c_fft[N123[IDX2F(ig,ik,NGWX)]-1] = -0.5*c_fft[N123[IDX2F(ig,ik,NGWX)]-1];
    }*/
    c_fft[N123[IDX2F(ig,ik,NGWX)]-1] = (beta*gg/(gg + q0))*c_fft[N123[IDX2F(ig,ik,NGWX)]-1];
  }

  // Back to real space
  fft_fftw3(c_fft,NR1,NR2,NR3,true);

  // Set new electron density
  d = 0.0;
  for(i=0; i<nnr; i++) {
    d = d + (rhoe[i] - rho_old[i])*(rhoe[i] - rho_old[i]);
    //rhoe[i] = rho_old[i] + c_fft[i].x + beta*(rhoe[i] - rho_old[i]);
    //rhoe[i] = rho_old[i] + beta*c_fft[i].x;
    rhoe[i] = rho_old[i] + creal(c_fft[i]);
    rho_old[i] = rhoe[i];
  }
  d = sqrt(d/(double)nnr);
  printf("In KerkerRhoMix: d = %f\n", d);

  // Free memory
  free(c_fft); c_fft=NULL;

}

void KerkerPrec(double *diff_rhoe, int idx1, int nnr, double beta, double q0)
{
  double complex *c_fft=NULL;
  int ig,i;
  int ik=1; // TODO: consider only the first k-point in list?
  double gg;

  c_fft = (double complex*)malloc(nnr*sizeof(double complex));

  for(i=1; i<=nnr; i++) {
    c_fft[i-1] = diff_rhoe[IDX2F(i,idx1,nnr)] + I*0.0;
  }
  // Transform to reciprocal space (forward FFT)
  fft_fftw3(c_fft,NR1,NR2,NR3,false);
  // Filter
  for(ig=1; ig<=NGW[ik-1]; ig++) {
    gg = XKG[IDX2F(ig,ik,NGWX)];
    c_fft[N123[IDX2F(ig,ik,NGWX)]-1] = (beta*gg/(gg + q0))*c_fft[N123[IDX2F(ig,ik,NGWX)]-1];
  }
  // Back to real space
  fft_fftw3(c_fft,NR1,NR2,NR3,true);
  // Copy to diff_rhoe
  for(i=1; i<=nnr; i++) {
    diff_rhoe[IDX2F(i,idx1,nnr)] = creal(c_fft[i-1]);
  }

  // Free memory
  free(c_fft); c_fft=NULL;
}

