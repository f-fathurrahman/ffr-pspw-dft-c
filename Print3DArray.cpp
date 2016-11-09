// eFeFeR, January 2012

#include "common_pspw_cuda.h"

void Print3DArray(double *A, int dim1, int dim2, int dim3, string filename)
{
  int N = dim1*dim2*dim3;

  FILE *outfile = fopen(&filename[0], "w");
  if(outfile==NULL) {
    printf("Error opening file\n");
    exit(1);
  }

  for(int i=0; i<N; i++) {
    fprintf(outfile,"%18.10f\n",A[i]);
  }

  fclose(outfile);

  printf("Data written successfully\n");
}

