// eFeFeR (20910015), October 2011

#include "common_pspw_cuda.h"

int main(int argc, char **argv)
{
  // Local variables
  double t1,t2;
  double time_total;

  // TODO: use strcmp?
  if(argc==2) {
    string st(&argv[1][0]);
    if(st == "--use-cuda") {
      USE_CUDA = true;
    } else {
      printf("ERROR: Unrecognized option: %s\n", argv[1]);
      return 0;
    }
  }

  if(argc > 2) {
    printf("!!!!! WARNING: Arguments ignored, only one argument needed.\n");
  }

  printf("Program started at: %s\n", timestring());

  t1 = cpu_time();

  // Read input file
  ReadInputFile();

  // Setup several global variables and allocate arrays
  Setup();
  
  PrintSystemInfo();

  PrepareSCF();

#ifdef _USE_CUDA
  if(USE_CUDA) {
    SetupCUDA();
  }
#endif

#ifdef _USE_CUDA
  if(USE_CUDA) {
    ShutdownCUDA();
  }
#endif

  printf("\n");
  printf("Program ended normally at: %s\n", timestring());
  t2 = cpu_time();
  
  time_total = t2-t1;

  printf("Total elapsed time for main program: %18.5f s\n", time_total);

  return 0;
}


