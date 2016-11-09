#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))
#define IDX3F(i,j,k,DIM1,DIM2) (((k)-1)*(DIM1)*(DIM2) + ((j)-1)*(DIM1) + ((i)-1))
#define IDX4F(i,j,k,l,DIM1,DIM2,DIM3) (((l)-1)*(DIM1)*(DIM2)*(DIM3) + ((k)-1)*(DIM1)*(DIM2) + ((j)-1)*(DIM1) + ((i)-1))

int main(int argc, char **argv)
{
  double ZVAL;
  int nang;

  printf("argc = %d\n", argc);
  if(argc != 2) {
    printf("One argument (input file) must be given\n");
    return -1;
  }

  ifstream inFile(argv[1]);
  if(!inFile.is_open()) {
    printf("ERROR reading file: %s\n", argv[1]);
    return -1;
  }

  string line;

  getline(inFile,line);
  sscanf(&line[0],"%lf %d", &ZVAL, &nang);
  printf("ZVAL = %f\n", ZVAL);
  printf("nang = %d\n", nang);

  // Skip ten lines
  // We are not using Gaussian parameters for pseudopotentials
  for(int i=1; i<=10; i++) {
    getline(inFile,line);
  }

  int is = 1;
  int mmax;
  double clog;
  // Read atomic wavefunction and pseudopotential on radial mesh
  for(il=1; il<=nang) {
    getline(inFile,line);
    sscanf(&line[0],"%d %lf",&mmax, &clog);
    printf("mmax = %d\n", mmax);
    printf("clog = %f\n", clog);
    for(ir=1; ir<=mmax; ir++) {
      getline(inFile,line);
      //sscanf(&line[0], "%d %lf %lf %lf", idummy, R[IDX2F(ir,1,)]);
    }
  }
  inFile.close();

  printf("\nProgram ended normally\n");
  return 0;
}
