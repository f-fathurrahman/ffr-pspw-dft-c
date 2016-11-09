#include <stdio.h>
#include <stdlib.h>

#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))

void tqli(double *d, double *e, int n, double *z);
void eigsrt(double *d, double *v, int n);

int main(int argc, char **argv)
{
  int N = 4;
  double *d, *sd;
  double *evec;
  int i, j;

  d = (double*)malloc(N*sizeof(double));
  sd = (double*)malloc(N*sizeof(double));
  evec = (double*)malloc(N*N*sizeof(double));

  d[0] = 1.0;
  d[1] = 4.0;
  d[2] = 9.0;
  d[3] = 16.0;
  sd[0] = 0.0;
  sd[1] = 1.0;
  sd[2] = 2.0;
  sd[3] = 3.0;

  for(j=1; j<=N*N; j++) evec[j-1] = 0.0;
  for(j=1; j<=N; j++) evec[IDX2F(j,j,N)] = 1.0;

  for(i=1; i<=N; i++) {
    for(j=1; j<=N; j++) {
      printf("%f ",evec[IDX2F(i,j,N)]);
    }
    printf("\n");
  }

  printf("Calling tqli ...");
  tqli(d, sd, N, evec);
  printf("....Done\n");

  eigsrt(d, evec, N);

  printf("Eigenvalues:\n");
  for(i=1; i<=N; i++) {
    printf("%f\n",d[i-1]);
  }

  printf("Eigenvectors:\n");
  for(i=1; i<=N; i++) {
    for(j=1; j<=N; j++) {
      printf("%f ",evec[IDX2F(i,j,N)]);
    }
    printf("\n");
  }

  return 0;
}
