// eFeFeR (20910015), November 2011

#include "common_pspw_cuda.h"

void SCFLoop()
{
  int i;
  double complex *c_fft=NULL;
  double *rhoe2=NULL, *rho_old=NULL, *rhoe=NULL;
  double complex zs;
  double ds;
  int ngi, iterSCF, ik;
  double complex *eigvec=NULL;
  double *new_eig=NULL;
  double etot_old, detot;
  double d_rms_rho;
  // Arrays for Anderson mixing (IMIX=2)
  double *f_anderson=NULL, *g_anderson=NULL, *hist_anderson=NULL;
  // Arrays for Pulay mixing (IMIX=3)
  double *hist_p=NULL, *f_p=NULL;
  // Arrays for Broyden mixing (IMIX=4)
  double *f_b=NULL, *df_b=NULL, *u_b=NULL, *a_b=NULL, *hist_b=NULL;


//
// Allocate memory
//
  c_fft = (double complex*)malloc(NNR*sizeof(double complex));
  for(i=0; i<NNR; i++) c_fft[i] = Z_ZERO;
  //
  rhoe = (double*)malloc(NNR*sizeof(double));
  for(i=0; i<NNR; i++) rhoe[i] = 0.0;
  //
  rhoe2 = (double*)malloc(NNR*sizeof(double));
  for(i=0; i<NNR; i++) rhoe2[i] = 0.0;
  //
  rho_old = (double*)malloc(NNR*sizeof(double));
  for(i=0; i<NNR; i++) rho_old[i] = 0.0;
  //
  new_eig = (double*)malloc(NX*sizeof(double));
  for(i=0; i<NX; i++) new_eig[i] = 0.0;
  //
  eigvec = (double complex*)malloc(NGWX*NX*sizeof(double complex));
  for(i=0; i<NGWX*NX; i++) eigvec[i] = Z_ZERO;

//
// Allocate memory for rho mixing
//
  // Anderson mixing
  if(IMIX==2) {
    f_anderson = (double*)malloc(NNR*3*sizeof(double));
    for(i=0; i<NNR*3; i++) f_anderson[i] = 0.0;
    g_anderson = (double*)malloc(NNR*2*sizeof(double));
    for(i=0; i<NNR*2; i++) g_anderson[i] = 0.0;
    hist_anderson = (double*)malloc(NNR*3*sizeof(double));
    for(i=0; i<NNR*3; i++) hist_anderson[i] = 0.0;
  }
  if(IMIX==3) {
    f_p = (double*)malloc(NNR*MIXDIM*sizeof(double));
    for(i=0; i<NNR*MIXDIM; i++) f_p[i] = 0.0;
    hist_p = (double*)malloc(NNR*MIXDIM*sizeof(double));
    for(i=0; i<NNR*MIXDIM; i++) hist_p[i] = 0.0;
  }
  // Broyden mixing
  if(IMIX==4) {
    f_b = (double*)malloc(NNR*2*sizeof(double));
    for(i=0; i<NNR*2; i++) f_b[i] = 0.0;
    df_b = (double*)malloc(NNR*MIXDIM*sizeof(double));
    for(i=0; i<NNR*MIXDIM; i++) df_b[i] = 0.0;
    u_b = (double*)malloc(NNR*MIXDIM*sizeof(double));
    for(i=0; i<NNR*MIXDIM; i++) u_b[i] = 0.0;
    a_b = (double*)malloc(MIXDIM*MIXDIM*sizeof(double));
    for(i=0; i<MIXDIM*MIXDIM; i++) a_b[i] = 0.0;
    hist_b = (double*)malloc(NNR*2*sizeof(double));
    for(i=0; i<NNR*2; i++) hist_b[i] = 0.0;
  }

  // Form starting density from atomic charges
  printf("\n");
  printf("Starting density from superposition of atomic charges\n");
  FormFactorAtomic(c_fft);

  // Inverse FFT
  fft_fftw3(c_fft, NR1, NR2, NR3, true);
  for(i=0; i<NNR; i++) rhoe2[i] = creal(c_fft[i]);

  // Normalize charge density
  printf("Normalize charge density ...");
  NormalizeRho(rhoe2,NNR,OMEGA,NEL);
  printf("...Done\n");

  Print3DArray(rhoe2,NR1,NR2,NR3,string("STARTING_RHO.dat"));

  // Save initial rho to rho_old
  dcopy_(&NNR, rhoe2, &I_ONE, rho_old, &I_ONE);

  // Determine number of planewaves
  // NEED THIS?
  ngi = 0;
  for(i=1; i<=NG; i++) {
    if(G[i-1] < GCUT) {
      ngi = ngi + 1;
      if(ngi > NGX) {
        printf("ERROR: ngi > NGX: %d %d\n", ngi, NGX);
      }
    }
  }

  EwaldEnergy();

  EvalEnergy(rhoe2);

  etot_old = ETOT;

  //test_ApplyHam(rhoe2);
  //abort();

#ifdef _USE_CUDA
  /*if(USE_CUDA) {
    CopyRho_ApplyHam_cuda(rhoe2);
    test_ApplyHam_cuda(rhoe2);
  }
  ShutdownCUDA();
  printf("After testing Apply_Ham_cuda\n");
  exit(1);*/
#endif

//
// Begin SCF iteration
//
  for(iterSCF=1; iterSCF<=MAX_SCF_ITER; iterSCF++) {

    printf("\n");
    printf("=====SCF iteration: %d=====\n", iterSCF);

    // Reset kinetic energy component
    EKIN = 0.0;

    // Reset output charge density
    for(i=0; i<NNR; i++) rhoe[i] = 0.0;

    // Loop over k-point
    for(ik=1; ik<=NKPT; ik++) {
      printf("\n");
      printf("----------- k-point: %5d -----------\n", ik);
      printf("\n");

      double t1,t2;
      t1 = cpu_time();
#if _USE_CUDA
      if(USE_CUDA) {
        CopyRho_ApplyHam_cuda(rhoe2);
        if(IDIAG==1) {
          Diagon_lobpcg_cuda(ik, new_eig, eigvec, NGW[ik-1], -10.0);
        } else {
          Diagon_davidson_cuda(ik, new_eig, eigvec, NGW[ik-1], -10.0);
        }
      } else {
        if(IDIAG==1) {
          Diagon_lobpcg_cpu(ik, new_eig, eigvec, NGW[ik-1], rhoe2, NNR, -10.0);
        } else {
          Diagon_davidson_cpu(ik, new_eig, eigvec, NGW[ik-1], rhoe2, NNR, -10.0);
        }
      }
#else
      if(IDIAG==1) {
        Diagon_lobpcg_cpu(ik, new_eig, eigvec, NGW[ik-1], rhoe2, NNR, -10.0);
      } else {
        Diagon_davidson_cpu(ik, new_eig, eigvec, NGW[ik-1], rhoe2, NNR, -10.0);
      }
#endif
      t2 = cpu_time();

      printf("\nTime for diagonalization: %f\n", t2-t1);

      for(i=1; i<=NX; i++) {
        EIG[IDX2F(i,ik,NX)] = new_eig[i-1];
      }

      // Evaluate new wavefunctions and electron density
      EvalRhoPsi(ik, rhoe, eigvec, NGW[ik-1]);

      // Evaluate non local pseudopotential contribution to energy
      EvalENL();

    } // end of loop over k-point

    // TODO: Call FermiOcc if tmetal

    // Copy new potential
    dcopy_(&NNR, rhoe,&I_ONE, rhoe2,&I_ONE);

    // Mix potential
    if(IMIX==1) {
      KerkerRhoMix(RHO_MIX_FAC,KERKER_Q0,rhoe2,rho_old,NNR,d_rms_rho);
    }
    else if(IMIX==2) {
      AndersonRhoMix(iterSCF,RHO_MIX_FAC,rhoe2,hist_anderson,f_anderson,g_anderson,NNR,d_rms_rho);
    }
    else if(IMIX==3) {
      PulayRhoMix(iterSCF,rhoe2,hist_p,f_p,MIXDIM,NNR,d_rms_rho);
    }
    else if(IMIX==4) {
      BroydenRhoMix(iterSCF,0.25,0.01,rhoe2,hist_b,f_b,df_b,u_b,a_b,MIXDIM,NNR,d_rms_rho);
    }
    else {
      SimpleRhoMix(RHO_MIX_FAC,rhoe2,rho_old,NNR,d_rms_rho);
    }
    dcopy_(&NNR, rhoe2,&I_ONE, rhoe,&I_ONE);

    // Evaluate local energy contributions
    EvalEnergy(rhoe2);

    detot = fabs(ETOT - etot_old);
    printf("\n");
    printf("! %5d %18.10f %18.10f %18.10f\n", iterSCF, ETOT, detot, d_rms_rho);
    printf("\n");
    etot_old = ETOT;

    if(detot <= SCF_CONV_CRIT) {
      printf("***Convergence achieved in %d iterations.\n", iterSCF);
      break;
    }

  } // end of SCF iteration

  // Print electron density to file
  // print before rhoe2 becomes total effective potential
  Print3DArray(rhoe,NR1,NR2,NR3,string("FINAL_RHO.dat"));

  free(c_fft); c_fft=NULL;
  free(rhoe); rhoe=NULL;
  free(rhoe2); rhoe2=NULL;
  free(rho_old); rho_old=NULL;
  free(new_eig); new_eig=NULL;
  free(eigvec); eigvec=NULL;

  if(IMIX==2) {
    free(f_anderson); f_anderson=NULL;
    free(g_anderson); g_anderson=NULL;
    free(hist_anderson); hist_anderson=NULL;
  }

  if(IMIX==3) {
    free(f_p); f_p=NULL;
    free(hist_p); hist_p=NULL;
  }

  if(IMIX==4) {
    free(f_b); f_b=NULL;
    free(df_b); df_b=NULL;
    free(u_b); u_b=NULL;
    free(a_b); a_b[i]=NULL;
    free(hist_b); hist_b=NULL;
  }

}
