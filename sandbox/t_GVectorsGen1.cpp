#include "../common_pspw_cuda.h"

void GVectorsGen1();

int main(int argc, char **argv)
{
  NR1 = 20; NR2 = 20; NR3 = 20;
  GCUT = 88.855355;
	B1[0] = 1.0; B1[1] = 0.0; B1[2] = 0.0;
	B2[0] = 0.0; B2[1] = 1.0; B2[2] = 0.0;
	B3[0] = 0.0; B3[1] = 0.0; B3[2] = 1.0;

  GVectorsGen1();

  return 0;
}

