// eFeFeR (20910015), October 2011
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void sort1(int N, double *A, int *index);
void driver_orderf(int N, double *A, int *index);
void SortGVectors(int N, double *A, int *index);

int main(int argc, char **argv)
{
  int N=5;
  double A[N];
  int index[N];

  for(int i=0; i<N; i++) {
    A[i] = (double)rand()/(double)RAND_MAX;
    index[i] = 0;
    printf("A[%d] = %f index[%d] = %d\n", i, A[i], i, index[i]);
  }

  //sort1(N, A, index);
	//driver_orderf(N,A,index);
	SortGVectors(N,A,index);

	//printf("After sort:\n");
	//printf("After driver_orderf:\n");
	//printf("After SortGVectors:\n");
  for(int i=0; i<N; i++) {
    printf("A[%d] = %f index[%d] = %d\n", i, A[i], i, index[i]);
  }

  return 0;
}
