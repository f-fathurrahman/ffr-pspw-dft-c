// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void random_col(double complex *A, int nelem, long a);
// Implemented in LOBPCG
//void symmetrize(double complex *v, int ldv, int nrow, int ncol);

void davidson_cpu(int ik, double *evals, double complex *X, int nbasis, int nstates,
    double *rhoe, int nnr, double *prec)
{
  // Local variables
  double rnorm;
  int ist, ig, istep, max_dir, nconv, nelem;
  double machine_zero, TOLERANCE;
  double *res_tol=NULL, *res_norm=NULL, *evals_red=NULL;
  double complex *h_mat=NULL, *o_mat=NULL, *evec=NULL;
  double complex *HX=NULL, *R=NULL, *HR=NULL, *Xtemp=NULL;
  double complex tz;
  char transA, transB;
  int nstates2 = 2*nstates;

  res_tol = (double*)malloc(nstates*sizeof(double));
  res_norm = (double*)malloc(nstates*sizeof(double));
  evals_red = (double*)malloc(nstates2*sizeof(double));
  //
  h_mat = (double complex*)malloc(nstates2*nstates2*sizeof(double complex));
  o_mat = (double complex*)malloc(nstates2*nstates2*sizeof(double complex));
  evec = (double complex*)malloc(nstates2*nstates2*sizeof(double complex));
  //
  HX = (double complex*)malloc(nbasis*nstates*sizeof(double complex));
  R = (double complex*)malloc(nbasis*nstates*sizeof(double complex));
  HR = (double complex*)malloc(nbasis*nstates*sizeof(double complex));
  Xtemp = (double complex*)malloc(nbasis*nstates*sizeof(double complex));

  double mem = 4.0*nbasis*nstates*sizeof(double complex);
  mem = mem + 3.0*nstates2*nstates2*sizeof(double complex);
  mem = mem + 4.0*nstates*sizeof(double);

  printf("Allocated dynamic memory in davidson_cpu: %15.10f MB\n", mem/1024./1024.);

  // Initial guess of eigenvectors: NOT GOOD?
  /*for(int i=0; i<nbasis*nstates; i++) X[i] = Z_ZERO;
  for(ist=1; ist<=nstates; ist++) {
    X[IDX2F(ist,ist,nbasis)] = Z_ONE;
  }*/
  long seed = 100;
  for(ist=1; ist<=nstates; ist++) {
    random_col(&X[IDX2F(1,ist,nbasis)],nbasis,seed);
    seed = seed + 1;
  }
  // Orthogonalize the wavefunctions
  ortho_qr(X,nbasis,nstates);

  // Apply Hamiltonian to wavefunction
  ApplyHam_block(&X[IDX2F(1,1,nbasis)], &HX[IDX2F(1,1,nbasis)],nbasis,ik,rhoe,nnr,nstates);
  
  // Calculate Rayleigh quotient
  for(ist=1; ist<=nstates; ist++) {
    printf("State: %5d\n", ist);
    zdotc_(&tz, &nbasis, &X[IDX2F(1,ist,nbasis)],&I_ONE, &HX[IDX2F(1,ist,nbasis)],&I_ONE);
    evals[ist-1] = creal(tz);
  }

  // Calculate matrix of residual vector
  for(ist=1; ist<=nstates; ist++) {
    for(ig=1; ig<=nbasis; ig++) {
      R[IDX2F(ig,ist,nbasis)] = evals[ist-1]*X[IDX2F(ig,ist,nbasis)] -
        HX[IDX2F(ig,ist,nbasis)];
    }
    zdotc_(&tz, &nbasis, &R[IDX2F(1,ist,nbasis)],&I_ONE, &R[IDX2F(1,ist,nbasis)],&I_ONE);
    res_tol[ist-1] = sqrt( creal(tz) );
  }

  if(!USER_MAX_DIAG_ITER) {
    max_dir = 30;
    if(max_dir < nstates/2) max_dir = nstates/2;
  }
  else {
    max_dir = MAX_DIAG_ITER;
  }
  machine_zero = 2.220446049250313e-16;
  TOLERANCE = 1.0e-6;
  rnorm = 1.0;

  nconv = 0;
  for(istep=1; istep<=max_dir; istep++) {
    if(DIAG_VERBOSE) printf("davidson_cpu: %5d %5d %18.10f\n", istep, nconv, rnorm);

    for(ist=1; ist<=nstates; ist++) res_norm[ist-1] = 1.0;

    for(ist=1; ist<=nstates; ist++) {
      if(res_tol[ist-1] > machine_zero) res_norm[ist-1] = 1.0/res_tol[ist-1];
    }

    // Scale the residual vectors
    for(ist=1; ist<=nstates; ist++) {
      zdscal_(&nbasis, &res_norm[ist-1], &R[IDX2F(1,ist,nbasis)],&I_ONE);
    }

    // Apply preconditioner
    for(ist=1; ist<=nstates; ist++) {
      for(ig=1; ig<=nbasis; ig++) {
        R[IDX2F(ig,ist,nbasis)] = prec[ig-1]*R[IDX2F(ig,ist,nbasis)];
      }
    }


// Construct the reduced hamiltonian. The reduced hamiltonian has dimensions
//  2nb x 2nb and is constructed by filling in four nb x nb blocks one at a time:
// __ 
//|  |
//| <x|H|x>   <x|H|r>  |
//    h_mat = |  |
//| *******   <r|H|r>  | 
//|__|

    // Apply Hamiltonian to residual matrix
    ApplyHam_block(&R[IDX2F(1,1,nbasis)], &HR[IDX2F(1,1,nbasis)],nbasis,ik,rhoe,nnr,nstates);

    if(istep == 1) {
      transA = 'C'; transB = 'N';
      //h_mat(1:nstates,1:nstates) = cmat
      zgemm_(&transA,&transB,&nstates,&nstates,&nbasis, &Z_ONE,X,&nbasis, HX,&nbasis,
          &Z_ZERO,&h_mat[IDX2F(1,1,nstates2)],&nstates2);
    } else {
      //h_mat(1:nstates,1:nstates) = diag(evals)
      for(int i=1; i<=nstates; i++) {
        for(int j=1; j<=nstates; j++) {
          if(i==j) h_mat[IDX2F(j,i,nstates2)] = evals[i-1]*Z_ONE;
          else h_mat[IDX2F(j,i,nstates2)] = Z_ZERO;
        }
      }
    }
    
    // <x|H|r> --> cmat
    transA = 'C'; transB = 'N';
    //h_mat(1:nstates,nstates+1:2*nstates) = cmat
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis, &Z_ONE,X,&nbasis, HR,&nbasis,
        &Z_ZERO,&h_mat[IDX2F(1,nstates+1,nstates2)],&nstates2);
    //h_mat(nstates+1:2*nstates,1:nstates) = cmat'
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis, &Z_ONE,HR,&nbasis, X,&nbasis,
        &Z_ZERO,&h_mat[IDX2F(nstates+1,1,nstates2)],&nstates2);
    
    // <r|H|r> --> cmat
    // h_mat(nstates+1:2*nstates,nstates+1:2*nstates) = cmat
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis, &Z_ONE,R,&nbasis, HR,&nbasis,
        &Z_ZERO,&h_mat[IDX2F(nstates+1,nstates+1,nstates2)],&nstates2);
    

// Construct the reduced overlap matrix which has dimenstions 2nb x 2nb
//   and is constructed by filling in four nb x nb blocks one at a time:
// _   _ 
//|     |
//|  <v|v>   <v|r>  |
//    o_mat = |     |
//|  *****   <r|r>  | 
//|_   _|

    // o_mat(1:nstates,1:nstates) = IDENTITY
    for(int i=1; i<=nstates; i++) {
      for(int j=1; j<=nstates; j++) {
        if(i==j) o_mat[IDX2F(j,i,nstates2)] = Z_ONE;
        else o_mat[IDX2F(j,i,nstates2)] = Z_ZERO;
      }
    }

    // <X|R> --> cmat
    // o_mat(1:nstates,nstates+1:2*nstates) = cmat
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis, &Z_ONE,X,&nbasis, R,&nbasis,
        &Z_ZERO,&o_mat[IDX2F(1,nstates+1,nstates2)],&nstates2);
    // o_mat(nstates+1:2*nstates,1:nstates) = cmat'
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis, &Z_ONE,R,&nbasis, X,&nbasis,
        &Z_ZERO,&o_mat[IDX2F(nstates+1,1,nstates2)],&nstates2);
    
    // <r|r> --> cmat
    // o_mat(nstates+1:2*nstates,nstates+1:2*nstates) = cmat
    zgemm_(&transA,&transB,&nstates,&nstates,&nbasis, &Z_ONE,R,&nbasis, R,&nbasis,
        &Z_ZERO,&o_mat[IDX2F(nstates+1,nstates+1,nstates2)],&nstates2);

    // Diagonalize reduced problems
    eig_zhegv_eval(h_mat,nstates2, o_mat,nstates2, evals_red, evec,nstates2, nstates2);

    // evals = evals_red(1:nstates)
    dcopy_(&nstates, evals_red,&I_ONE, evals,&I_ONE);

    // cmat = evec(1:nstates,1:nstates)
    // v*cmat --> v
    transA = 'N'; transB = 'N';
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates, &Z_ONE,X,&nbasis, evec,&nstates2,
        &Z_ZERO,Xtemp,&nbasis);
    // v = Xtemp
    nelem = nbasis*nstates;
    zcopy_(&nelem, Xtemp,&I_ONE, X,&I_ONE);
    
    // HX = HX*cmat
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates, &Z_ONE,HX,&nbasis, evec,&nstates2,
        &Z_ZERO,Xtemp,&nbasis);
    // HX = Xtemp;
    zcopy_(&nelem, Xtemp,&I_ONE, HX,&I_ONE);

    // cmat = evec(nstates+1:2*nstates,1:nstates)
    // X = X + R*cmat
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates, &Z_ONE,R,&nbasis,
        &evec[IDX2F(nstates+1,1,nstates2)],&nstates2, &Z_ONE,X,&nbasis);
    // HX = HX + HR*cmat
    zgemm_(&transA,&transB,&nbasis,&nstates,&nstates, &Z_ONE,HR,&nbasis,
        &evec[IDX2F(nstates+1,1,nstates2)],&nstates2, &Z_ONE,HX,&nbasis);

    // Calculate matrix of residual vector
    for(ist=1; ist<=nstates; ist++) {
      for(ig=1; ig<=nbasis; ig++) {
        R[IDX2F(ig,ist,nbasis)] = evals[ist-1]*X[IDX2F(ig,ist,nbasis)] -
          HX[IDX2F(ig,ist,nbasis)];
      }
      zdotc_(&tz, &nbasis, &R[IDX2F(1,ist,nbasis)],&I_ONE, &R[IDX2F(1,ist,nbasis)],&I_ONE);
      res_tol[ist-1] = sqrt( creal(tz) );
    }
    
    nconv = 0;
    for(ist=1; ist<=nstates; ist++) {
      if(res_tol[ist-1] < TOLERANCE) nconv = nconv + 1;
    }

    if(nconv >= nstates-NEMPTY) {
      printf("In iter = %5d: nconv = %5d\n", istep, nconv);
      break;
    }

    rnorm = 0.0;
    for(ist=1; ist<=nstates; ist++) {
      rnorm = rnorm + res_tol[ist-1];
    }
    rnorm = rnorm/(double)nstates;

    if(rnorm < TOLERANCE/10.0) {
      printf("davidson_cpu: Convergence achieved in iter = %5d: nconv = %5d\n", istep, nconv);
      break;
    }

  }
  
  //!rnorm = sum(res_tol)/real(nstates,8)

  rnorm = 0.0;
  printf("\n");
  printf("Final Davidson result nconv=%d from %d states:\n", nconv, nstates);
  for(ist=1; ist<=nstates; ist++) {
    printf("%5d %18.10f %18.10f\n", ist, evals[ist-1], res_tol[ist-1]);
    rnorm = rnorm + res_tol[ist-1];
  }
  rnorm = rnorm/(double)nstates;
  printf("End of Davidson iteration: rnorm = %18.10f\n", rnorm);
  
  // Free memory
  free(res_tol); res_tol=NULL;
  free(res_norm); res_norm=NULL;
  free(h_mat); h_mat=NULL;
  free(o_mat); o_mat=NULL;
  free(evec); evec=NULL;
  free(evals_red); evals_red=NULL;
  //
  free(HX); HX=NULL;
  free(R); R=NULL;
  free(HR); HR=NULL;
  free(Xtemp); Xtemp=NULL;
}


