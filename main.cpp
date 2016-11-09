// eFeFeR (20910015), October 2011

#include "common_pspw.h"

int main(int argc, char **argv)
{
  // Local variables
  double t1,t2, t3,t4;
  double time_total, time_scf;

  if(argc > 2) {
    printf("!!!!! WARNING: Arguments ignored, only one argument needed.\n");
  }

  printf("\n");
  printf("Program started at: %s\n", timestring());
  printf("\n");

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

  t3 = cpu_time();
  if(ISCF==1) {
    ChebySCFLoop();
  } else {
    SCFLoop();
  }
  t4 = cpu_time();
  time_scf = t4-t3;

#ifdef _USE_CUDA
  if(USE_CUDA) {
    ShutdownCUDA();
  }
#endif

  printf("\n");
  printf("Program ended normally at: %s\n", timestring());
  t2 = cpu_time();
  
  time_total = t2-t1;

  printf("Total time for SCF loop:             %18.5f s\n", time_scf);
  printf("Total elapsed time for main program: %18.5f s\n", time_total);

  return 0;
}

