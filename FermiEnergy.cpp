// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void FermiEnergy(int istepe)
{
  // Local
  double *eigv=NULL, *eigvb=NULL, *fik=NULL;
  double *dfocc=NULL;
  int *kpoint=NULL;
  const double AU = 2.0*13.6058;
  const int niter=200;
  double efermi, all1;
  int ik, ist, ikn;
  double epsn = 1.0e-7;
  int ic;

  // Allocate memory
  eigv = (double*)malloc(NSTATX*sizeof(double));
  eigvb = (double*)malloc(NSTATX*sizeof(double));
  fik = (double*)malloc(NSTATX*sizeof(double));
  //
  dfocc = (double*)malloc(NX*NKPT*sizeof(double));
  //
  kpoint = (int*)malloc(NSTATX*sizeof(int));

  // For hydrogen atom
  if(NX==1 && NKPT==1) {
    FOCC[IDX2F(1,1,NX)] = NEL;
    efermi = EIG[IDX2F(1,1,NX)]*AU;
  }

  if(!TMETAL) return;

  for(ikn=1; ikn<=NSTATX; ikn++) {
    eigv[ikn-1] = 100000.0; // arbitrary large number
  }

  for(ik=1; ik<=NKPT; ik++) {
    for(ist=1; ist<=NX; ist++) {
      ikn = (ik-1)*NX + ist;
      eigv[ikn-1] = EIG[IDX2F(ist,ik,NX)]*AU;
    }
  }

  // if(t_first)
  if(istepe==2) {
    sort2(NSTATX,eigv,eigvb,kpoint);
    all1 = 0.0;
    for(ikn=1; ikn<=NSTATX; ikn++) {
      d0 = NEL - all1;
      all1 = all1 + WKPT[(kpoint[ikn-1]-1)/NX];
      d1 = NEL - all1;
      if(fabs(d1) > fabs(d0)) {
        nfermi = ikn - 1;
        if(NEL/2.0 - (int)(NEL/2.0) < 0.1) {
          efermi = (eigvb[nfermi-1] + eigvb[nfermi])/2.0;
        } else {
          efermi = eigvb[nfermi-1];
        }
        goto LABEL2000;
      }
    }
  }

LABEL2000:
  
  ic = 0;
  def = 0.0;
  max_efermi =  40000.0;
  min_efermi = -40000.0;

LABEL1000:
  for(ikn=1; ikn<=NSTATX; ikn++) {
    x = (eigv[ikn-1] - efermi)/EKT;
    if(x <= -30.0) {
      fik[ikn-1] = 2.0;
    }
    else if(x < 30.0) {
      fik[ikn-1] = 2.0/(1.0 + exp(x));
    }
    else {
      fik[ikn-1] = 0.0;
    }
  }

  fsum = 0.0;
  dsums = 0.0;
  for(ikn=1; ikn<=NSTATX; ikn++) {
    k = (ikn-1)/NX + 1;
    fsum = fsum + fik[ikn-1]*WKPT[k-1];
    // How many electrons?
    dsums = dsums + fik[ikn-1]*(1.0 - 0.5*fik[ikn-1])*WKPT[k-1];
  }

  // Test if Fermi energy is already converged
  if(fabs(fsum-NEL) >= epsn) {
    ic = ic + 1;
    if(ic > niter) {
      printf("ERROR: Fermi energy is not found after maximum number of iterations: %5d\n");
      abort();
    }
    else {
      if(fsum-NEL > 0.0 && efermi < max_efermi) max_efermi = efermi;
      if(fsum-NEL < 0.0 && efermi > min_efermi) min_efermi = efermi;
      efermi = (max_efermi + min_efermi)/2.0;
      goto LABEL1000;
    }
  }

  seq = 0.0;
  dpmax = 0.0;
  for(ikn=1; ikn<=NSTATX; ikn++) {
    ik = (ikn-1)/NX + 1;
    in1 = ikn - (ik-1)*NX;
    pp = fik[ikn-1];
    dp = pp - FOCC[IDX2F(in1,ik,NX)];
    if(fabs(dp) > dpmax) dpmax = fabs(dp);
    dfocc[IDX2F(in1,ik,NX)] = dp;
  }

  // Free memory
  free(eigv);
  free(eigvb);
  free(fik);
  free(dfocc);
  free(kpoint);
}


