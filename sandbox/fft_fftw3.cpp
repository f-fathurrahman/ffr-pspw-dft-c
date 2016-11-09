// Fadjar Fathurrahman (20910015), November 2011

#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cuda/cuComplex.h>
#include <fftw3.h>

using namespace std;

typedef cuDoubleComplex complex16;

extern "C" void zdscal_(int *n, double *a, fftw_complex *x, int *incx);

void fft_fftw3(complex16 *data, int NR1, int NR2, int NR3, bool inverse)
{
  int N = NR1*NR2*NR3;
  double scale = 1.0/double(N);
  int incx=1;
	fftw_complex *dat=NULL;
	dat = (fftw_complex*)fftw_malloc(N*sizeof(fftw_complex));
	memcpy(dat,data,N*sizeof(fftw_complex));

  if(!inverse) {
    fftw_plan plan_forward;
    plan_forward = fftw_plan_dft_3d(NR3, NR2, NR1, dat, dat, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);
    zdscal_(&N, &scale, dat, &incx);
    fftw_destroy_plan(plan_forward);
  } else {
    fftw_plan plan_backward;
    plan_backward = fftw_plan_dft_3d(NR3, NR2, NR1, dat, dat, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backward);
    fftw_destroy_plan(plan_backward);
  }

	memcpy(data,dat,N*sizeof(fftw_complex));
	fftw_free(dat); dat=NULL;

}

