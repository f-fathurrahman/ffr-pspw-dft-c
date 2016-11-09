// eFeFeR (20910015), November 2011

#include <iostream>
#include <mkl/mkl_dfti.h>
#include <cuda/cuComplex.h>

typedef cuDoubleComplex complex16;
extern "C" void zdscal_(int *n, double *a, complex16 *x, int *incx);

void fft_mkl(complex16 *data, int NR1, int NR2, int NR3, bool inverse)
{
  DFTI_DESCRIPTOR_HANDLE fft_handle;
  MKL_LONG status;
  int dim[3];
  dim[0] = NR1; dim[1] = NR2; dim[2] = NR3;
  int N = NR1*NR2*NR3;
  double scale = 1.0/(double)N;
  int incx=1;

  // Create FFTI descriptor
  status = DftiCreateDescriptor(&fft_handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, dim);
  if(status != 0) {
    printf("ERROR creating FFTI descriptor.\n");
    abort();
  }
  // Committing FFT descriptor
  status = DftiCommitDescriptor(fft_handle);
  if(status != 0) {
    printf("ERROR committing FFTI descriptor.\n");
    abort();
  }
  // Compute FFT
  if(!inverse) { // forward FFT
    status = DftiComputeForward(fft_handle,data);
    zdscal_(&N, &scale, data, &incx);
    if(status != 0) {
      printf("ERROR computing forward FFT.\n");
      abort();
    }
  } else { // Backward FFT
    status = DftiComputeBackward(fft_handle,data);
    if(status != 0) {
      printf("ERROR computing backward FFT.\n");
      abort();
    }
  } 
  status = DftiFreeDescriptor(&fft_handle);
  if(status != 0) {
    printf("ERROR freeing FFT descriptor.\n");
    abort();
  }

}

