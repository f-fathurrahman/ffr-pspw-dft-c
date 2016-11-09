// eFeFeR (20910015), October 2011

#include "common_pspw_cuda.h"

void CalcRecipLatt()
{
  double scale;
  
  scale = ALAT/OMEGA;
  
  B1[0] = ( A2[1]*A3[2] - A2[2]*A3[1] )*scale;
  B1[1] = ( A1[2]*A3[1] - A1[1]*A3[2] )*scale;
  B1[2] = ( A1[1]*A2[2] - A1[2]*A2[1] )*scale;
  
  B2[0] = ( A2[2]*A3[0] - A2[0]*A3[2] )*scale;
  B2[1] = ( A1[0]*A3[2] - A1[2]*A3[0] )*scale;
  B2[2] = ( A1[2]*A2[0] - A1[0]*A2[2] )*scale;

  B3[0] = ( A2[0]*A3[1] - A2[1]*A3[0] )*scale;
  B3[1] = ( A1[1]*A3[0] - A1[0]*A3[1] )*scale;
  B3[2] = ( A1[0]*A2[1] - A1[1]*A2[0] )*scale;

}
