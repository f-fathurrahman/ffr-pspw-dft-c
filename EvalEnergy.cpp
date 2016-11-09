// eFeFeR, December 2011

#include "common_pspw_cuda.h"

void EvalEnergy(double *rhoe)
{
  printf("Evaluating Kohn-Sham energy terms:\n");

  // Local
  int i, is, ig;
  double complex *osum=NULL, *vtemp=NULL, *v_rhets=NULL, *v_vcgs=NULL;
  double complex *c_fft=NULL;
  double *Vxc_vec=NULL;
  double Vxc;
  double complex eh, eh2, eg2, vp, eps, rp, rps, rhet, rhets, rhog, rhogs, vcg;
  double complex sumz;
  double sumd;
  double vh0, fpibg;

// Allocate arrays
  osum = (double complex*)malloc(3*NAX*NSP*sizeof(double complex));
  c_fft = (double complex*)malloc(NNR*sizeof(double complex));
  vtemp = (double complex*)malloc(NGX*sizeof(double complex));
  v_rhets = (double complex*)malloc(NGX*sizeof(double complex));
  v_vcgs = (double complex*)malloc(NGX*sizeof(double complex));
  Vxc_vec = (double*)malloc(NNR*sizeof(double));

  for(i=0; i<3*NAX*NSP; i++) osum[i] = Z_ZERO;
  for(i=0; i<NNR; i++) c_fft[i] = Z_ZERO;
  for(i=0; i<NGX; i++) vtemp[i] = Z_ZERO;
  for(i=0; i<NGX; i++) v_rhets[i] = Z_ZERO;
  for(i=0; i<NGX; i++) v_vcgs[i] = Z_ZERO;
  for(i=0; i<NNR; i++) Vxc_vec[i] = 0.0;

  for(i=0; i<NNR; i++) c_fft[i] = rhoe[i] + I*0.0;
  /*SumArray(rhoe, NNR, sumd);
  printf("sum(rhoe) = %f\n", sumd);
  SumArray(c_fft, NNR, sumz);
  printf("Before: sum_c_fft = (%f,%f)\n", sumz.x, sumz.y);*/

  fft_fftw3(c_fft, NR1,NR2,NR3, false); // forward FFT

  /*SumArray(c_fft, NNR, sumz);
  printf("sum_c_fft = (%f,%f)\n", sumz.x, sumz.y);*/

  eh = Z_ZERO;
  eh2 = Z_ZERO;
  eg2 = Z_ZERO;
  vp = Z_ZERO;

  for(is=1; is<=NSP; is++) {
    vp = vp + SFAC[IDX2F(is,1,NSP)]*VPS[IDX2F(is,1,NSP)];
  }
  //printf("vp = (%f,%f)\n", vp.x, vp.y);

  // eps is G=0 Coulomb energy of the charge density and the local
  // pseudopotential field of Gaussian pseudocharges
  vh0 = 0.0;
  eps = vp*conj(c_fft[0]);
  vtemp[0] = vp + vh0; // need this?

  for(ig=2; ig<=NG; ig++) {
    vp = Z_ZERO;
    rp = Z_ZERO;

    for(is=1; is<=NSP; is++) {
      vp = vp + SFAC[IDX2F(is,ig,NSP)]*VPS[IDX2F(is,ig,NSP)];
      rp = rp + SFAC[IDX2F(is,ig,NSP)]*RHOPS[IDX2F(is,ig,NSP)];
    }

    // rhet = charge density in Fourier space
    rhet = c_fft[IDX3F(N1[ig-1],N2[ig-1],N3[ig-1],NR1,NR2)];

    // rhog = sum of electronic and pseudo charge
    rhog = rhet + rp;
    rhets = conj(rhet);
    v_rhets[ig-1] = rhets;
    rhogs = conj(rhog);
    rps = conj(rp);

    fpibg = 4.0*M_PI/(TPIBA2*G[ig-1]);

    // vcg = Coulomb potential at G
    vcg = fpibg*rhog;
    v_vcgs[ig-1] = conj(vcg);

    // Sum up Hartree (electron + pseudoatom charge density)
    eh = eh + 0.5*vcg*rhogs;

    // Sum up Hartee (only electrons) energy
    eh2 = eh2 + 0.5*fpibg*rhet*rhets;

    // Sum up Hartree (only pseudoatom charge) energy
    eg2 = eg2 + 0.5*fpibg*rp*rps;

    // Sum up pseudopotential energy of electronic charge
    eps = eps + rhets*vp;

    // vtemp = Hartree potential + local part of pseudopotential
    vtemp[ig-1] = vcg + vp;
  } // ig

  /*SumArray(vtemp,NGX,sumz);
  printf("sum(vtemp) = (%f,%f)\n", sumz.x, sumz.y);
  SumArray(v_vcgs,NGX,sumz);
  printf("sum(v_vcgs) = (%f,%f)\n", sumz.x, sumz.y);*/

  // TODO
  // Compute forces here?

  // Transform electrostatic potential (Hartree + local part of pseudopotential)
  for(ig=1; ig<=NG; ig++) {
    c_fft[IDX3F(N1[ig-1],N2[ig-1],N3[ig-1],NR1,NR2)] = vtemp[ig-1];
  }
  fft_fftw3(c_fft,NR1,NR2,NR3,true); // inverse transform

  //SumArray(c_fft,NNR,sumz);
  //printf("sum(c_fft) = (%f,%f)\n", sumz.x, sumz.y);

  // TODO
  // Calculate correction to surface calculation

  // Evaluate exchange-correlation energy and potential
  xc_PZ(rhoe,Vxc_vec,EXC,NNR);

  // rhoe contains total effective potential now
  Vxc = 0.0;
  for(i=0; i<NNR; i++) {
    Vxc = Vxc + rhoe[i]*Vxc_vec[i];
    rhoe[i] = creal(c_fft[i]) + Vxc_vec[i];
  }

  // TODO: Write out effective and electrostatic potential to file

  // TODO: Compute contribution to local forces to total forces

//
// Calculation of contribution of local potential to energy
//
  Vxc = Vxc*OMEGA/NNR;
  EXC = EXC*OMEGA/NNR;
  EHT = creal(eh)*OMEGA + ESR - ESELF;
  EPSEU = creal(eps)*OMEGA;
  ETOT = EKIN + EHT + EPSEU + ENL + EXC;


  printf("\n");
  printf("ETOT   = %18.10f\n", ETOT);
  printf("EKIN   = %18.10f\n", EKIN);
  printf("EHT    = %18.10f\n", EHT);
  printf("EPSEU  = %18.10f\n", EPSEU);
  printf("ENL    = %18.10f\n", ENL);
  printf("EXC    = %18.10f\n", EXC);

  // Free memory
  free(osum); osum=NULL;
  free(c_fft); c_fft=NULL;
  free(vtemp); vtemp=NULL;
  free(v_rhets); v_rhets=NULL;
  free(v_vcgs); v_vcgs=NULL;
  free(Vxc_vec); Vxc_vec=NULL;

}
