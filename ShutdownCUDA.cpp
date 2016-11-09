// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void ShutdownCUDA()
{
  printf("Shutdown CUDA ....\n");
#ifdef _USE_CUBLAS2
  CALL( cublasDestroy(CUBLAS_HANDLE) );
#else
  cublasShutdown();
#endif
  FreeArrays_ApplyHam_cuda();
  FreeRho_ApplyHam_cuda();
  printf("....done...\n");
}

