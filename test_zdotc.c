#include <complex.h>
#include <stdio.h>
#define N 5

int main()
{
  int n = N, inca = 1, incb = 1, i;
  double complex a[N], b[N], c;

  for( i = 0; i < n; i++ ){
    a[i] = (double)i + i*3.3*I;
    b[i] = (double)(n - i) + i * 2.0*I;
    printf("a = %f %f\n", creal(a[i]), cimag(a[i]));
    printf("b = %f %f\n", creal(b[i]), cimag(b[i]));
  }

  zdotc_( &c, &n, &a[0], &inca, &b[0], &incb );

  printf( "The complex dot product is: ( %6.2f, %6.2f)\n", creal(c), cimag(c) );

  return 0;

}

