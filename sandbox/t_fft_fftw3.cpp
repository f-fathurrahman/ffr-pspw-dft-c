// eFeFeR (20910015), November 2011

#include <cstdio>
#include <cstdlib>
#include <cuda/cuComplex.h>
#include <complex>
#include <fftw3.h>

using namespace std;

typedef cuDoubleComplex complex16;
#define IDX3F(i,j,k,DIM1,DIM2) (((k)-1)*(DIM1)*(DIM2) + ((j)-1)*(DIM1) + ((i)-1))

void fft_fftw3(complex16 *data, int NR1, int NR2, int NR3, bool inverse);

int main(int argc, char **argv)
{
  int NR1=3, NR2=3, NR3=2;
  int i, j, k;
  complex16 *dat=NULL;

  // Initialize data
  dat = (complex16*)malloc(NR1*NR2*NR3*sizeof(complex16));
  for(i=1; i<=NR1*NR2*NR3; i++) {
    dat[i-1] = make_cuDoubleComplex(1.0,2.0);
  }
  //dat[IDX3F(1,2,1,NR1,NR2)] = make_cuDoubleComplex(2.1,3.3);
  //dat[IDX3F(2,2,1,NR1,NR2)] = make_cuDoubleComplex(1.1,4.3);

  // Call FFT
  fft_fftw3(dat, NR1, NR2, NR3, true);

  printf("Result\n");
  for(k=1; k<=NR3; k++) {
    for(j=1; j<=NR2; j++) {
      for(i=1; i<=NR1; i++) {
        //printf("(%d,%d,%d): (%f,%f)\n", i,j,k, dat[IDX3F(i,j,k,NR1,NR2)][0], dat[IDX3F(i,j,k,NR1,NR2)][1]);
        printf("(%d,%d,%d): (%f,%f)\n", i,j,k, dat[IDX3F(i,j,k,NR1,NR2)].x, dat[IDX3F(i,j,k,NR1,NR2)].y);
	//printf("(%d,%d,%d): (%f,%f)\n", i,j,k, dat[IDX3F(i,j,k,NR1,NR2)].real(), dat[IDX3F(i,j,k,NR1,NR2)].imag());
      }
    }
  }

  free(dat); dat=NULL;

  printf("Program ended normally\n");
  return 0;
}


