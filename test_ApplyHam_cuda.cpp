// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void test_ApplyHam_cuda(double *rhoe)
{
  double complex *vin=NULL, *d_vin;
  double complex *vout=NULL, *vout_gpu=NULL, *d_vout;
  int ik=1;
  int i;
  double sumd;
  double complex tz;
  double t1, t2;
  int ntry=1;
  unsigned int sumi;

  SumArray(N123,(unsigned int)NGWX*NKPT,sumi);
  printf("sum(N123) = %d\n", sumi);

  printf("\n");
  printf("***Calling test_ApplyHam_cuda\n");

  if(ik > NKPT) {
    printf("ik > NKPT: %d %d\n",ik,NKPT);
    abort();
  }

  // Allocate memory on CPU
  vin = (double complex*)malloc(NGW[ik-1]*sizeof(double complex));
  vout = (double complex*)malloc(NGW[ik-1]*sizeof(double complex));
  vout_gpu = (double complex*)malloc(NGW[ik-1]*sizeof(double complex));
  // Allocate memory on GPU
  CALL( cudaMalloc((void**)&d_vin, NGW[ik-1]*sizeof(double complex)) );
  CALL( cudaMalloc((void**)&d_vout, NGW[ik-1]*sizeof(double complex)) );

  // Initialize vin
  for(i=0; i<NGW[ik-1]; i++) vin[i] = make_double complex(1.0,2.0);

  // Copy vin to GPU
  CALL( cudaMemcpy(d_vin, vin, NGW[ik-1]*sizeof(double complex), cudaMemcpyHostToDevice) );
  

  int G_NT = 512; // TODO: Tune these parameters....
  int G_NB = NGW[ik-1]/G_NT + ( (NGW[ik-1]%G_NT) ? 1:0 );
  int R_NT = 512;
  int R_NB = NNR/R_NT + ( (NNR%R_NT) ? 1:0 );

  printf("ntry = %d\n", ntry);
  printf("NGW  = %d\n", NGW[ik-1]);
  printf("(G_NB,G_NT) = (%d,%d) %d\n", G_NB,G_NT, G_NB*G_NT);
  printf("NNR = %d\n", NNR);
  printf("(R_NB,R_NT) = (%d,%d) %d\n", R_NB,R_NT, R_NB*R_NT);

  // Test ApplyHam
  t1 = cpu_time();
  for(i=1; i<=ntry; i++) {
    ApplyHam_block(vin,vout,NGW[ik-1],ik,rhoe,NNR,1);
  }
  t2 = cpu_time();
  printf("\n");
  printf("Time elapsed for ApplyHam_block: %f\n", t2-t1);
  zdotc_(&tz, &NGW[ik-1], vout,&I_ONE, vout,&I_ONE);
  printf("***CHECK: dot(vout,vout) = (%18.10f,%18.10f)\n", tz.x, tz.y);

  t1 = cpu_time();
  // Test ApplyHam_block_cuda
  for(i=1; i<=ntry; i++) {
    ApplyHam_block_cuda(d_vin, d_vout, NGW[ik-1],ik,1);
  }
  t2 = cpu_time();
  // Copy result back to CPU
  CALL( cudaMemcpy(vout_gpu, d_vout, NGW[ik-1]*sizeof(double complex), cudaMemcpyDeviceToHost) );
  printf("\n");
  printf("Time elapsed for ApplyHam_block_cuda: %f\n", t2-t1);
  zdotc_(&tz, &NGW[ik-1], vout_gpu,&I_ONE, vout_gpu,&I_ONE);
  printf("***CHECK: dot(vout_gpu,vout_gpu) = (%18.10f,%18.10f)\n", tz.x, tz.y);

  // Free memory
  free(vin); vin=NULL;
  free(vout); vout=NULL;
  free(vout_gpu); vout_gpu=NULL;
  //
  CALL( cudaFree(d_vin) );
  CALL( cudaFree(d_vout) );
}


