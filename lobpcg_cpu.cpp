// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

int random3(int lim)
{
  static long a = 3;
  
  a = (((a * 214013L + 2531011L) >> 16) & 32767);
  
  return ((a % lim) + 1);
}


void symmetrize(double complex *v, int ldv, int nrow, int ncol)
{
  int i, j;
  double complex tmp;
  // Loop over upper tridiagonal part
  for(j=1; j<=ncol; j++) {
    for(i=1; i<=j; i++) {
      tmp = 0.5*(v[IDX2F(i,j,ldv)] + conj(v[IDX2F(j,i,ldv)]));
      v[IDX2F(i,j,ldv)] = tmp;
      v[IDX2F(j,i,ldv)] = conj(tmp);
    }
  }
}

void random_col(double complex *A, int nelem, long a);

void lobpcg_cpu(int ik, double *lambda, double complex *X, int nbasis, int nstates,
    double *rhoe, int nnr, double *prec)
{
  // Constants
  const double TOL = 1.0e-6;
  const double TFUDGE = 1.e10;
  int maxiter;
  if(!USER_MAX_DIAG_ITER) {
    maxiter = 30;
    if(maxiter < nstates/2) maxiter = nstates/2;
  }
  else {
    maxiter = MAX_DIAG_ITER;
  }

  // Local
  int i, j, ist;
  int nelem;
  char transA, transB, uplo, side, diag;
  int iter, info;
  double complex *Q=NULL, *HQ=NULL;
  double complex *temp1=NULL, *T=NULL, *G=NULL, *tempX=NULL, *U=NULL;
  double *resnrm=NULL;
  double rnorm;
  double complex sumz, tz;
  double sumd;
  int nconv, nlock;
  int nstates2 = nstates*2;
  int nstates3 = nstates*3;

  Q = (double complex*)malloc(nbasis*nstates3*sizeof(double complex)); //
  HQ = (double complex*)malloc(nbasis*nstates3*sizeof(double complex)); //
  temp1 = (double complex*)malloc(nstates*nstates*sizeof(double complex)); //
  T = (double complex*)malloc(nstates3*nstates3*sizeof(double complex)); //
  G = (double complex*)malloc(nstates3*nstates3*sizeof(double complex)); //
  tempX = (double complex*)malloc(nbasis*nstates*sizeof(double complex)); //
  U = (double complex*)malloc(nstates3*nstates3*sizeof(double complex)); //
  resnrm = (double*)malloc(nstates*sizeof(double));

  double mem = (7.0*nbasis*nstates3 + nstates*nstates + 3.0*nstates3*nstates3)*16.0;
  mem = mem + nstates*8.0;
  printf("Allocated dynamic memory in lobpcg_cpu: %15.10f MB\n", mem/1024./1024.);

//
// Create initial guess wavefunction
//
  /*double re, im;
  for(i=0; i<nbasis*nstates; i++) {
    Q[i].x = random3(111)/(double)111;
    Q[i].y = random3(111)/(double)111;
  }*/
  long seed = 100;
  for(ist=1; ist<=nstates; ist++) {
    random_col(&Q[IDX2F(1,ist,nbasis)],nbasis,seed);
    seed = seed + 1;
  }
  //for(i=0; i<nbasis*nstates; i++) Q[i] = Z_ONE;
  //for(i=1; i<=nstates; i++) Q[IDX2F(i,i,nbasis)] = 10*Z_ONE;
  // Orthogonalize the wavefunctions
  ortho_qr(Q,nbasis,nstates);
   
//
// Apply Hamiltonian
//
  ApplyHam_block(&Q[IDX2F(1,1,nbasis)], &HQ[IDX2F(1,1,nbasis)],nbasis,ik,rhoe,nnr,nstates);

//
// First iteration, pulled out of the loop
//

  // Calculate XHX = X* HX
  transA = 'C'; transB = 'N';
  zgemm_(&transA,&transB,&nstates,&nstates,&nbasis,&Z_ONE,Q,&nbasis,HQ,&nbasis,&Z_ZERO,
      temp1,&nstates);
  // Calculate residual vectors
  nelem = nbasis*nstates;
  zcopy_(&nelem, HQ,&I_ONE, &Q[IDX2F(1,nstates+1,nbasis)],&I_ONE); // W = HX
  //
  transA = 'N'; transB = 'N';
  zgemm_(&transA,&transB,&nbasis,&nstates,&nstates, &MZ_ONE,Q,&nbasis, temp1,&nstates,
      &Z_ONE,&Q[IDX2F(1,nstates+1,nbasis)],&nbasis);
  // Diagonalize temp1
  eig_zheevd(temp1, nstates, lambda, nstates);

//
// Check convergence
//
  nconv = 0;
  nlock = 0;
  rnorm = 0.0;
  for(i=1; i<=nstates; i++) {
    zdotc_(&tz, &nbasis, &Q[IDX2F(1,nstates+i,nbasis)],&I_ONE, &Q[IDX2F(1,nstates+i,nbasis)],&I_ONE);
    resnrm[i-1] = sqrt(creal(tz));
    if(resnrm[i-1] < TOL) nconv = nconv + 1;
    if(resnrm[i-1] < TOL/TFUDGE) nlock = nlock + 1;
    rnorm = rnorm + resnrm[i-1];
  }
  rnorm = rnorm/(double)nstates;

  // TODO: nconv >= number of occupied states
  if(nconv >= nstates-NEMPTY) {
    printf("Convergence achieved for %d eigenvalues in LOBPCG at first iteration.\n", nconv);
    goto LABEL10;
  }

  if(rnorm < TOL) {
    printf("Convergence achieved for %d eigenvalues in LOBPCG (rnorm).\n", nconv);
    goto LABEL10;
  }

  // Apply preconditioner
  for(i=1; i<=nstates; i++) {
    for(int ib=1; ib<=nbasis; ib++) {
      Q[IDX2F(ib,nstates+i,nbasis)] = prec[ib-1]*Q[IDX2F(ib,nstates+i,nbasis)];
    }
  }

  // HW = H*W
  ApplyHam_block(&Q[IDX2F(1,nstates+1,nbasis)], &HQ[IDX2F(1,nstates+1,nbasis)],
      nbasis,ik,rhoe,nnr,nstates);

  // C = W* W
  transA = 'C'; transB = 'N';
  zgemm_(&transA,&transB,&nstates,&nstates,&nbasis,
      &Z_ONE,&Q[IDX2F(1,nstates+1,nbasis)],&nbasis, &Q[IDX2F(1,nstates+1,nbasis)],&nbasis,
      &Z_ZERO,temp1,&nstates);
  // C = (C + C*)/2 TODO: Need this to avoid spurious imaginary part?
  symmetrize(temp1, nstates, nstates, nstates);

  // Cholesky decomposition
  uplo = 'U';
  zpotrf_(&uplo,&nstates,temp1,&nstates,&info);
  if(info != 0) {
    printf("ERROR calculating Cholesky decomposition\n");
    abort();
  }

//
// Solve system of linear equations
//
  side = 'R'; uplo = 'U'; transA = 'N'; diag = 'N';

  // W = W/C
  ztrsm_(&side,&uplo,&transA,&diag, &nbasis,&nstates, &Z_ONE, temp1,&nstates,
      &Q[IDX2F(1,nstates+1,nbasis)],&nbasis);

  // HW = HW/C
  ztrsm_(&side,&uplo,&transA,&diag, &nbasis,&nstates, &Z_ONE, temp1,&nstates,
      &HQ[IDX2F(1,nstates+1,nbasis)],&nbasis);

  // T = Q* HQ
  transA = 'C'; transB = 'N';
  zgemm_(&transA,&transB,&nstates2,&nstates2,&nbasis, &Z_ONE,Q,&nbasis, HQ,&nbasis,
      &Z_ZERO,T,&nstates3);
  symmetrize(T,nstates3, nstates2,nstates2);

  // G = Q* Q
  transA = 'C'; transB = 'N';
  zgemm_(&transA,&transB,&nstates2,&nstates2,&nbasis, &Z_ONE,Q,&nbasis, Q,&nbasis,
      &Z_ZERO,G,&nstates3);
  symmetrize(G,nstates3, nstates2,nstates2);

//
// Diagonalize T wrt G
//
  eig_zhegv(T,nstates3, G,nstates3, U,nstates3, nstates2);
  sumz = Z_ZERO; 

  nelem = nbasis*nstates;
  // X = Q U
  transA = 'N'; transB = 'N';
  zgemm_(&transA,&transB,&nbasis,&nstates,&nstates2, &Z_ONE,Q,&nbasis,
      U,&nstates3, &Z_ZERO,tempX,&nbasis);
  zcopy_(&nelem, tempX,&I_ONE, Q,&I_ONE);
  // HX = HQ U
  transA = 'N'; transB = 'N';
  zgemm_(&transA,&transB,&nbasis,&nstates,&nstates2, &Z_ONE,HQ,&nbasis,
      U,&nstates3, &Z_ZERO,tempX,&nbasis);
  zcopy_(&nelem, tempX,&I_ONE, HQ,&I_ONE);
  // P = W
  zcopy_(&nelem, &Q[IDX2F(1,nstates+1,nbasis)],&I_ONE, &Q[IDX2F(1,nstates2+1,nbasis)],&I_ONE);
  // HP = HW
  zcopy_(&nelem, &HQ[IDX2F(1,nstates+1,nbasis)],&I_ONE, &HQ[IDX2F(1,nstates2+1,nbasis)],&I_ONE);

//
// Main LOBPCG iterations
//
  for(iter=2; iter<=maxiter; iter++) {
    //
    if(DIAG_VERBOSE) printf("iter = %5d %5d %18.10f\n", iter, nconv, rnorm);
    //
    // Calculate XHX = X* HX
    transA = 'C'; transB = 'N';
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis,&Z_ONE,Q,&nbasis,HQ,&nbasis,&Z_ZERO,
        temp1,&nstates);
    // Calculate residual vectors
    nelem = nbasis*nstates;
    zcopy_(&nelem, HQ,&I_ONE, &Q[IDX2F(1,nstates+1,nbasis)],&I_ONE); // W = HX
    //
    transA = 'N'; transB = 'N';
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates, &MZ_ONE,Q,&nbasis, temp1,&nstates,
        &Z_ONE,&Q[IDX2F(1,nstates+1,nbasis)],&nbasis);
    // Diagonalize temp1
    eig_zheevd(temp1, nstates, lambda, nstates);

    //
    // Check convergence
    //
    nconv = 0;
    nlock = 0;
    rnorm = 0.0;
    for(i=1; i<=nstates; i++) {
      zdotc_(&tz, &nbasis, &Q[IDX2F(1,nstates+i,nbasis)],&I_ONE,
          &Q[IDX2F(1,nstates+i,nbasis)],&I_ONE);
      resnrm[i-1] = sqrt(creal(tz));
      if(resnrm[i-1] < TOL) nconv = nconv + 1;
      if(resnrm[i-1] < TOL/TFUDGE) nlock = nlock + 1;
      rnorm = rnorm + resnrm[i-1];
    }
    rnorm = rnorm/(double)nstates;

    // TODO: nconv >= number of occupied states
    if(nconv >= nstates-NEMPTY) {
      printf("Convergence achieved for %d eigenvalues in LOBPCG at iter = %d\n", nconv, iter);
      goto LABEL10;
    }

    if(rnorm < TOL) {
      printf("Convergence achieved for %d eigenvalues in LOBPCG (rnorm) at iter = %d\n",
          nconv, iter);
      goto LABEL10;
    }

    // Apply preconditioner
    for(i=1; i<=nstates; i++) {
      for(int ib=1; ib<=nbasis; ib++) {
        Q[IDX2F(ib,nstates+i,nbasis)] = prec[ib-1]*Q[IDX2F(ib,nstates+i,nbasis)];
      }
    }

    // HW = H*W
    ApplyHam_block(&Q[IDX2F(1,nstates+1,nbasis)], &HQ[IDX2F(1,nstates+1,nbasis)],
      nbasis,ik,rhoe,nnr,nstates);

    // C = W* W
    transA = 'C'; transB = 'N';
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis,
        &Z_ONE,&Q[IDX2F(1,nstates+1,nbasis)],&nbasis, &Q[IDX2F(1,nstates+1,nbasis)],&nbasis,
        &Z_ZERO,temp1,&nstates);
    // C = (C + C*)/2 TODO: Need this to avoid spurious imaginary part?
    symmetrize(temp1,nstates,nstates,nstates);

    // Cholesky decomposition
    uplo = 'U';
    zpotrf_(&uplo,&nstates,temp1,&nstates,&info);
    if(info != 0) {
      printf("ERROR calculating Cholesky decomposition\n");
      abort();
    }
    // Solve system of linear equations
    side = 'R'; uplo = 'U'; transA = 'N'; diag = 'N';
    // W = W/C
    ztrsm_(&side,&uplo,&transA,&diag, &nbasis,&nstates, &Z_ONE, temp1,&nstates,
        &Q[IDX2F(1,nstates+1,nbasis)],&nbasis);
    // HW = HW/C
    ztrsm_(&side,&uplo,&transA,&diag, &nbasis,&nstates, &Z_ONE, temp1,&nstates,
        &HQ[IDX2F(1,nstates+1,nbasis)],&nbasis);

    // T = Q* HQ
    transA = 'C'; transB = 'N';
    zgemm_(&transA,&transB,&nstates3,&nstates3,&nbasis, &Z_ONE,Q,&nbasis, HQ,&nbasis,
        &Z_ZERO,T,&nstates3);
    symmetrize(T,nstates3,nstates3,nstates3);

    // G = Q* Q
    transA = 'C'; transB = 'N';
    zgemm_(&transA,&transB,&nstates3,&nstates3,&nbasis, &Z_ONE,Q,&nbasis, Q,&nbasis,
        &Z_ZERO,G,&nstates3);
    symmetrize(G,nstates3,nstates3,nstates3);

    //
    // Diagonalize T wrt G
    //
    eig_zhegv(T,nstates3, G,nstates3, U,nstates3, nstates3);

    nelem = nbasis*nstates;
    // X = Q U
    transA = 'N'; transB = 'N';
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates3, &Z_ONE,Q,&nbasis,
        U,&nstates3, &Z_ZERO,tempX,&nbasis);
    zcopy_(&nelem, tempX,&I_ONE, Q,&I_ONE);
    // HX = HQ U
    transA = 'N'; transB = 'N';
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates3, &Z_ONE,HQ,&nbasis,
        U,&nstates3, &Z_ZERO,tempX,&nbasis);
    zcopy_(&nelem, tempX,&I_ONE, HQ,&I_ONE);
    //P = matmul(W,U((nstates+1):(2*nstates),:)) + matmul(P,U((2*nstates+1):(3*nstates),:))
    transA = 'N'; transB = 'N';
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates,&Z_ONE,&Q[IDX2F(1,nstates2+1,nbasis)],&nbasis,
        &U[IDX2F(nstates2+1,1,nstates3)],&nstates3,
        &Z_ZERO,tempX,&nbasis);
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates,&Z_ONE,&Q[IDX2F(1,nstates+1,nbasis)],&nbasis,
        &U[IDX2F(nstates+1,1,nbasis)], &nstates3,
        &Z_ZERO,&Q[IDX2F(1,nstates2+1,nbasis)],&nbasis);
    zaxpy_(&nelem, &Z_ONE, tempX,&I_ONE, &Q[IDX2F(1,nstates2+1,nbasis)],&I_ONE);
    //
    // HP = matmul(HW,U((nstates+1):(2*nstates),:)) + matmul(HP,U((2*nstates+1):(3*nstates),:))
    transA = 'N'; transB = 'N';
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates,&Z_ONE,&HQ[IDX2F(1,nstates2+1,nbasis)],&nbasis,
        &U[IDX2F(nstates2+1,1,nstates3)],&nstates3,
        &Z_ZERO,tempX,&nbasis);
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates,&Z_ONE,&HQ[IDX2F(1,nstates+1,nbasis)],&nbasis,
        &U[IDX2F(nstates+1,1,nbasis)], &nstates3,
        &Z_ZERO,&HQ[IDX2F(1,nstates2+1,nbasis)],&nbasis);
    zaxpy_(&nelem, &Z_ONE, tempX,&I_ONE, &HQ[IDX2F(1,nstates2+1,nbasis)],&I_ONE);
    //
    // C = P* P
    transA = 'C'; transB = 'N';
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis,
        &Z_ONE,&Q[IDX2F(1,nstates2+1,nbasis)],&nbasis,&Q[IDX2F(1,nstates2+1,nbasis)],&nbasis,
        &Z_ZERO,temp1,&nstates);
    symmetrize(temp1,nstates,nstates,nstates);

    //
    // Cholesky decomposition
    //
    uplo = 'U';
    zpotrf_(&uplo,&nstates,temp1,&nstates,&info);
    if(info != 0) {
      printf("ERROR calculating Cholesky decomposition\n");
      abort();
    }
    // Solve system of linear equations
    side = 'R'; uplo = 'U'; transA = 'N'; diag = 'N';
    // P = P/C
    ztrsm_(&side,&uplo,&transA,&diag, &nbasis,&nstates, &Z_ONE, temp1,&nstates,
        &Q[IDX2F(1,nstates2+1,nbasis)],&nbasis);
    // HP = HP/C
    ztrsm_(&side,&uplo,&transA,&diag, &nbasis,&nstates, &Z_ONE, temp1,&nstates,
        &HQ[IDX2F(1,nstates2+1,nbasis)],&nbasis);
  }

LABEL10:
//
// Copy wavefunction to output argument
//
  nelem = nbasis*nstates;
  zcopy_(&nelem, Q,&I_ONE, X,&I_ONE);

  printf("Final LOBPCG result\n");
  for(i=1; i<=nstates; i++) {
    printf("%d %18.10f %18.10f\n", i, lambda[i-1], resnrm[i-1]);
  }

  free(Q); Q=NULL;
  free(HQ); HQ=NULL;
  free(temp1); temp1=NULL;
  free(T); T=NULL;
  free(G); G=NULL;
  free(U); U=NULL;
  free(tempX); tempX=NULL;
  free(resnrm); resnrm=NULL;

  return;
}


