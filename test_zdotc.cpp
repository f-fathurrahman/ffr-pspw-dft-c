#include <complex>
#include <iostream>
#define N 5

extern "C" void zdotc_(std::complex<double> *pres, int *n, std::complex<double> *x,
    int *incx, std::complex<double> *y, int *incy);

int main()
{
  int n, inca = 1, incb = 1, i;
  std::complex<double> a[N], b[N], c;
  n = N;
       
  for( i = 0; i < n; i++ ){
    a[i] = std::complex<double>(i,i*2.3);
    b[i] = std::complex<double>(i,i*2.0);
  }
  zdotc_(&c, &n, a, &inca, b, &incb );
  std::cout << "The complex dot product is: " << c << std::endl;
  return 0;
}


