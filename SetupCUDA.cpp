// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void SetupCUDA()
{
#ifdef _USE_CUBLAS2
  CALL( cublasCreate(&CUBLAS_HANDLE) );
#else
  cublasInit();
#endif 
  SetupArrays_ApplyHam_cuda(N123, XKG,EIGR, PKG_A);
  SetupRho_ApplyHam_cuda();
}

