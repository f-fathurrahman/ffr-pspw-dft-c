// eFeFeR (20910015), January 2012

#include "common_pspw.h"

extern void symmetrize(double complex *v, int ldv, int nrow, int ncol);

double induced_norm_mat(double *A, int ldA, int N)
{
  if(ldA < N) {
    printf("ERROR: ldA < N = %d %d\n", ldA, N);
    exit(1);
  }
  // Create a unit-norm vector
  double d = 1.0/sqrt((double)N);
  double *v1=NULL;
  v1 = (double*)malloc(N*sizeof(double));
  for(int i=0; i<N; i++) v1[i] = d;
  // Result vector
  double *v=NULL;
  v = (double*)malloc(N*sizeof(double));

  // Calculate matrix-vector multiplication
  double D_ONE = 1.0;
  double D_ZERO = 0.0;
  char trans='N';
  dgemv_(&trans,&N,&N, &D_ONE,A,&ldA, v1,&I_ONE, &D_ZERO,v,&I_ONE);
  double res = ddot_(&N, v,&I_ONE, v,&I_ONE);

  free(v); v=NULL;
  free(v1); v=NULL;

  return res;

}

void ChebySCFLoop()
{
  int i;
  double complex *c_fft=NULL;
  double *rhoe2=NULL, *rho_old=NULL, *rhoe=NULL;
  double complex tz;
  double ds;
  int ngi, iterSCF, ik;
  //double complex *eigvec=NULL;
  double *new_eig=NULL;
  double etot_old, detot, d_rms_rho;
  double complex *v0;
  double t1, t2;
  // Arrays for Anderson mixing (IMIX=2)
  double *f_anderson=NULL, *g_anderson=NULL, *hist_anderson=NULL;
  // Arrays for Pulay mixing (IMIX=3)
  double *hist_p=NULL, *f_p=NULL;
  // Arrays for Broyden mixing (IMIX=4)
  double *f_b=NULL, *df_b=NULL, *u_b=NULL, *a_b=NULL, *hist_b=NULL;

  int nelem;
  int k;
  double complex *f=NULL;
  double *rT=NULL;
  double *work=NULL, *eval_rT=NULL;
  int lwork;
  char jobz; char uplo;
  int info;
  double lb, ub, norm_rT;
  double complex *X=NULL;
  double complex *HX=NULL;
  double complex *RES=NULL;
  int ig,is;
  double complex *T=NULL, *G=NULL;
  double *eval_T=NULL;
  char tA; char tB;
  int NX2 = NX + 2;

  if(NKPT > 1) {
    printf("Sorry: ChebySCFLoop only works for 1 k-point only for this moment ...\n");
    exit(1);
  }
  ik=1;

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
  /*eigvec = (double complex*)malloc(NGWX*NX*sizeof(double complex));
  for(i=0; i<NGWX*NX; i++) eigvec[i] = Z_ZERO;*/

  double r_re, r_im;
  v0 = (double complex*)malloc(NGWX*sizeof(double complex));
  for(i=0; i<NGWX; i++) {
    r_re = rand()/(double)RAND_MAX;
    r_im = rand()/(double)RAND_MAX;
    v0[i] = r_re + I*r_im;
  }

  X = (double complex*)malloc(NGW[ik-1]*NX2*sizeof(double complex));
  HX = (double complex*)malloc(NGW[ik-1]*NX2*sizeof(double complex));
  G = (double complex*)malloc(NX*NX*sizeof(double complex));
  T = (double complex*)malloc(NX2*NX2*sizeof(double complex));
  eval_T = (double*)malloc(NX2*sizeof(double));

  k = max(2*NX,10);
  rT = (double*)malloc(k*k*sizeof(double));
  f = (double complex*)malloc(NGWX*sizeof(double complex));

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
  
  EwaldEnergy();
  EvalEnergy(rhoe2);
  etot_old = ETOT;

  // Lanczos iteration
  t1 = cpu_time(); 
  lanczos(ik, v0, k, rT, k, f, rhoe2, NNR);
  t2 = cpu_time();
  printf("\nTime for 1st Lanczos iteration: %f\n", t2-t1);
  norm_rT = induced_norm_mat(rT,k,k);
  printf("norm_rT = %f\n", norm_rT);
  zdotc_(&tz, &NGW[ik-1], f,&I_ONE,f,&I_ONE);
  printf("norm(f) = %f\n", sqrt(creal(tz)));
  //
  lwork=3*k-1;
  work = (double*)malloc(lwork*sizeof(double));
  eval_rT = (double*)malloc(k*sizeof(double));
  //
  jobz='N'; uplo='U';
  dsyev_(&jobz, &uplo, &k, rT, &k, eval_rT, work, &lwork, &info);
  if(info != 0) {
    printf("ERROR calling dsyev_: info = %d\n", info);
    exit(1);
  }
  //
  lb = eval_rT[NX+1]; printf("lb = %f\n", lb);
  ub = eval_rT[2*NX-1]; printf("ub = %f\n", ub);

  // Initial wavefunctions
  // Random
  if(IDIAG < 0) {
    for(is=1; is<=NX+2; is++) {
      for(ig=1; ig<=NGW[ik-1]; ig++) {
        r_re = rand()/(double)RAND_MAX;
        r_im = rand()/(double)RAND_MAX;
        X[IDX2F(ig,is,NGW[ik-1])] = r_re + I*r_im;
      }
    }
    ortho_qr(X,NGW[ik-1],NX+2);
  }
  else {
    // First SCF iteration, using Davidson method
    iterSCF = 1;
    printf("\n");
    printf("First SCF iteration using Davidson diagonalization\n");
    printf("\n");
    EKIN = 0.0;
    for(i=0; i<NNR; i++) rhoe[i] = 0.0;
    double *prec=NULL;
    prec = (double*)malloc(NGW[ik-1]*sizeof(double));
    GenPrec(ik,prec,NGW[ik-1],-10.0);
    davidson_cpu(ik,eval_T,X,NGW[ik-1],NX2,rhoe2,NNR,prec);
    free(prec); prec=NULL;
    // Copy eigenvalues
    for(i=1; i<=NX; i++) {
      EIG[IDX2F(i,ik,NX)] = eval_T[i-1];
    }
    // Copy eigenvectors
    /*for(is=1; is<=NX; is++) {
      zcopy_(&NGW[ik-1],&X[IDX2F(1,is,NGW[ik-1])],&I_ONE,&eigvec[IDX2F(1,is,NGWX)],&I_ONE);
    }
    // Evaluate new wavefunctions and electron density
    EvalRhoPsi(ik, rhoe, eigvec, NGW[ik-1]);*/
    EvalRhoPsi(ik, rhoe, X, NGW[ik-1], NGW[ik-1]);
    // Evaluate non local pseudopotential contribution to energy
    EvalENL();
    // Copy new potential
    dcopy_(&NNR, rhoe,&I_ONE, rhoe2,&I_ONE);
    // Mix charge density
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
  }

  int iterStart;
  if(IDIAG < 0) {
    iterStart = 1;
  } else {
    iterStart = 2;
  }

  for(iterSCF=iterStart; iterSCF<=MAX_SCF_ITER; iterSCF++) {
    
    printf("\n");
    printf("=====ChebySCF iteration: %d=====\n", iterSCF);

    // Reset kinetic energy component
    EKIN = 0.0;

    // Reset output charge density
    for(i=0; i<NNR; i++) rhoe[i] = 0.0;
    
    printf("\n");
    printf("Calling chebyfilt with bounds: (%f,%f) and degree %d\n",lb,ub,CHEBYPOL_DEGREE);
    t1 = cpu_time();
    chebyfilt(ik, X, CHEBYPOL_DEGREE, lb, ub, rhoe2, NNR);
    t2 = cpu_time();
    printf("Time for chebyfilt: %f\n",t2-t1);

    // Orthonormalize the filtered wave functions
    ortho_qr(X,NGW[ik-1],NX2);

    // Compute Rayleigh quotient
    ApplyHam_block(X,HX,NGW[ik-1], ik, rhoe2, NNR, NX2);
    //
    tA = 'C'; tB='N';
    zgemm_(&tA,&tB,&NX2,&NX2,&NGW[ik-1], &Z_ONE,X,&NGW[ik-1], HX,&NGW[ik-1], &Z_ZERO,T,&NX2);
    symmetrize(T,NX2,NX2,NX2);
    eig_zheevd(T,NX2,eval_T,NX2);
    /*printf("Eigenvalue of T\n");
    for(is=1; is<=NX; is++) {
      printf("%5d %f\n", is, eval_T[is-1]);
    }*/
    tA = 'N'; tB='N';
    zgemm_(&tA,&tB,&NGW[ik-1],&NX2,&NX2, &Z_ONE,X,&NGW[ik-1], T,&NX2, &Z_ZERO,HX,&NGW[ik-1]);
    nelem = NGW[ik-1]*NX2;
    zcopy_(&nelem,HX,&I_ONE,X,&I_ONE);

    // Copy eigenvalues
    for(is=1; is<=NX; is++) {
      EIG[IDX2F(is,ik,NX)] = eval_T[is-1];
    }
    // Copy eigenvectors
    /*for(is=1; is<=NX; is++) {
      zcopy_(&NGW[ik-1],&X[IDX2F(1,is,NGW[ik-1])],&I_ONE,&eigvec[IDX2F(1,is,NGWX)],&I_ONE);
    }
    EvalRhoPsi(ik, rhoe, eigvec, NGW[ik-1]);*/
    EvalRhoPsi(ik, rhoe, X, NGW[ik-1], NGW[ik-1]);
    EvalENL();

    // Copy new electron density
    dcopy_(&NNR, rhoe,&I_ONE, rhoe2,&I_ONE);

    // Mix charge density
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

    // Update lb and ub
    lb = eval_T[NX+1];
    if(k > 10) {
      lanczos(ik, v0, 10, rT, k, f, rhoe2, NNR); // why it is 10?
    } else {
      lanczos(ik, v0, k, rT, k, f, rhoe2, NNR);
    }
    //lanczos(ik, v0, k, rT, k, f, rhoe2, NNR);
    zdotc_(&tz,&NGW[ik-1],f,&I_ONE,f,&I_ONE);
    if(k > 10) {
      ub = induced_norm_mat(rT,k,10) + sqrt(creal(tz));
    } else {
      ub = induced_norm_mat(rT,k,k) + sqrt(creal(tz));
    }
    //ub = induced_norm_mat(rT,k,k) + sqrt(tz.x);
  }

  printf("End of ChebySCFLoop\n");
  
  Print3DArray(rhoe,NR1,NR2,NR3,string("FINAL_RHO.dat"));

  free(v0); v0=NULL;
  free(T); T=NULL;
  free(eval_T); eval_T=NULL;
  //
  free(work); work=NULL;
  free(eval_rT); eval_rT=NULL;
  //
  free(rT); rT=NULL;
  free(f); f=NULL;
  // Free memory
  free(c_fft); c_fft=NULL;
  free(rhoe); rhoe=NULL;
  free(rho_old); rho_old=NULL;
  free(new_eig); new_eig=NULL;
  //free(eigvec); eigvec=NULL;

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


  //
  // Check residual error
  //
  printf("\n");
  printf("ChebySCF: Calculating residual error:\n");
  printf("\n");
  // Allocate memory
  RES = (double complex*)malloc(NGW[ik-1]*NX*sizeof(double complex));
  ApplyHam_block(X,HX,NGW[ik-1], ik, rhoe2, NNR, NX);
  tA = 'C'; tB='N';
  zgemm_(&tA,&tB,&NX,&NX,&NGW[ik-1], &Z_ONE,X,&NGW[ik-1], HX,&NGW[ik-1], &Z_ZERO,G,&NX);
  symmetrize(G,NX,NX,NX);
  tA = 'N'; tB='N';
  // RES <-- X*G
  zgemm_(&tA,&tB,&NGW[ik-1],&NX,&NX, &Z_ONE,X,&NGW[ik-1], G,&NX, &Z_ZERO,RES,&NGW[ik-1]);
  for(is=1; is<=NX; is++) {
    for(ig=1; ig<=NGW[ik-1]; ig++) {
      RES[IDX2F(ig,is,NGW[ik-1])] = HX[IDX2F(ig,is,NGW[ik-1])] - RES[IDX2F(ig,is,NGW[ik-1])];
    }
    zdotc_(&tz, &NGW[ik-1], &RES[IDX2F(1,is,NGW[ik-1])],&I_ONE, &RES[IDX2F(1,is,NGW[ik-1])],&I_ONE);
    printf("Residual: %5d %18.10f\n", is, sqrt(creal(tz)));
  }

  free(rhoe2); rhoe2=NULL;
  free(X); X=NULL;
  free(HX); HX=NULL;
  free(RES); RES=NULL;
  free(G); G=NULL;
}

