#include "mkl.h"
#include <stdio.h>

typedef complex16;

extern "C" void zdotc(complex16*, int *, complex16 *, int *, complex16 *, int *);

#define N 5

int main(int argc, char argv)
{
  int n, inca = 1, incb = 1, i;
  complex16 a[N], b[N], c;
  
  n = N;
  for( i = 0; i < n; i++ ){
    a[i].x = (double)i; a[i].y = (double)i * 2.0;
    b[i].x = (double)(n - i); b[i].y = (double)i * 2.0;
  }
  
  zdotc(&c, &n, a, &inca, b, &incb );
  
  printf( "The complex dot product is: (%6.2f, %6.2f)\n", c.x, c.y);

  return 0;
}

