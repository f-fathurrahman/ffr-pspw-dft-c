// eFeFeR (20910015), October 2011

#include "common_pspw.h"

void InitOcc()
{
  int ik, i;
  int n;

  //printf("Setting up initial occupation number\n");

  n = ceil(NEL/2.0);
  //printf("Highest occupied band = %d\n", n);

  // Fill the highest occupied band
  if((int)NEL%2 != 0) {
    for(ik=1; ik<=NKPT; ik++) FOCC[IDX2F(n,ik,NX)] = 1.0;
  } else {
    for(ik=1; ik<=NKPT; ik++) FOCC[IDX2F(n,ik,NX)] = 2.0;
  }
  // Fill remaining occupied bands
  for(i=1; i<=n-1; i++) {
    for(ik=1; ik<=NKPT; ik++) FOCC[IDX2F(i,ik,NX)] = 2.0;
  }
  // The empty bands
  for(i=n+1; i<=NX; i++) {
    for(ik=1; ik<=NKPT; ik++) FOCC[IDX2F(i,ik,NX)] = 0.0;
  }
  
  // Check initial occupation numbers
  double sum_el = 0.0;
  for(ik=1; ik<=NKPT; ik++) {
    for(i=1; i<=NX; i++) {
      sum_el = sum_el + FOCC[IDX2F(i,ik,NX)]*WKPT[ik-1];
    }
  }
  if(fabs(NEL-sum_el) > 1.e-5) {
    printf("ERROR: wrong initial occupation numbers\n");
    printf("Sum of occupation number = %f\n", sum_el);
    printf("is not equal to number of electrons = %f\n", NEL);
    printf("Check your k-point weights.\n");
    abort();
  }

  /*for(i=1; i<=NX; i++) {
    printf("Band number %d\n", i);
    for(ik=1; ik<=NKPT; ik++) {
      printf("%f ", FOCC[IDX2F(i,ik,NX)]);
    }
    printf("\n");
  }*/

  //printf("EXIT: %s\n", __FILE__);
}

