#include <stdio.h>
#include <stdlib.h>

extern "C" double ddot_(int *n, double *x, int *incx, double *y, int *incy);

int main(int argc, char **argv)
{
  double s;
  int N=4;
  double *x=NULL, *y=NULL;
  int I_ONE=1;

  x = (double*)malloc(N*sizeof(double));
  y = (double*)malloc(N*sizeof(double));

  for(int i=0; i<N; i++) {
    x[i] = 1.0;
    printf("%f\n", x[i]);
  }
  for(int i=0; i<N; i++) {
    y[i] = 2.0;
    printf("%f\n", y[i]);
  }

  s = ddot_(&N, x,&I_ONE, y,&I_ONE);
  printf("s = %f\n", s);

  free(x); x=NULL;
  free(y); y=NULL;

  return 0;
}
