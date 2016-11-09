#include "common_pspw.h"

void PrintGlobals()
{
	int ik;

	printf("k point list:\n");
  for(ik=1; ik<=NKPT; ik++) { 
    printf("%f %f %f\n", XK[IDX2F(1,ik,3)],XK[IDX2F(2,ik,3)], XK[IDX2F(3,ik,3)]);
  }
}
