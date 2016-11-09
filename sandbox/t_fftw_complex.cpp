#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cuda/cuComplex.h>
#include <fftw3.h>

typedef cuDoubleComplex complex16;

int main(int argc, char **argv)
{
  int i;
  int N=10;

  printf("sizeof(complex16) = %d\n", sizeof(complex16));
  printf("sizeof(fftw_complex) = %d\n", sizeof(fftw_complex));

  complex16 *data1=NULL;
  data1 = (complex16*)malloc(N*sizeof(complex16));
  for(i=0; i<10; i++) {
    data1[i] = make_cuDoubleComplex(i+1.0,i+2.0);
    printf("data1[%d] = (%f,%f) (%d,%d)\n", i, data1[i].x, data1[i].y, &data1[i].x, &data1[i].y);
  }

  fftw_complex *data2;
  data2 = (fftw_complex*)fftw_malloc(N*sizeof(fftw_complex));
  memcpy(data2, data1, N*sizeof(cuDoubleComplex));
  for(i=0; i<10; i++) {
    /*data2[i][0] = data1[i].x;
    data2[i][1] = data1[i].y;*/
    printf("data2[%d] = (%f,%f) (%d,%d)\n", i, data2[i][0], data2[i][1], &data2[i][0], &data2[i][1]);
  }

  free(data1); data1=NULL;
  fftw_free(data2);

  return 0;
}
