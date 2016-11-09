#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))
#define IDX3F(i,j,k,DIM1,DIM2) (((k)-1)*(DIM1)*(DIM2) + ((j)-1)*(DIM1) + ((i)-1))
#define IDX4F(i,j,k,l,DIM1,DIM2,DIM3) (((l)-1)*(DIM1)*(DIM2)*(DIM3) + ((k)-1)*(DIM1)*(DIM2) + ((j)-1)*(DIM1) + ((i)-1))

using namespace std;

int main(int argc, char **argv)
{
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
  vector<string> SPECIES;
  string speciesName;
  int NSP, NAX=0;
  int *NA;

  double ecut;
  double a1[3], a2[3], a3[3];
  int natom;
  string coord_type;
  int *LMAX, *LLOC;
  double *TAU;

  while(!inFile.eof()) {
    getline(inFile,line);
    if(line == "ecut") {
      inFile >> ecut;
      printf("ecut = %f\n", ecut);
    }
    else if(line == "latvec") {
      getline(inFile,line);
      sscanf(&line[0],"%lf %lf %lf", &a1[0], &a1[1], &a1[2]);
      getline(inFile,line);
      sscanf(&line[0],"%lf %lf %lf", &a2[0], &a2[1], &a2[2]);
      getline(inFile,line);
      sscanf(&line[0],"%lf %lf %lf", &a3[0], &a3[1], &a3[2]);
      printf("a1 = %f %f %f\n", a1[0], a1[1], a1[2]);
      printf("a2 = %f %f %f\n", a2[0], a2[1], a2[2]);
      printf("a3 = %f %f %f\n", a3[0], a3[1], a3[2]);
    }
    else if(line == "atoms") {
      getline(inFile,line);
      sscanf(&line[0],"%d", &NSP);
      printf("nsp = %d\n", NSP);
      NA = (int*)malloc(NSP*sizeof(int));
      LMAX = (int*)malloc(NSP*sizeof(int));
      LLOC = (int*)malloc(NSP*sizeof(int));
      NAX = 0;
      // Read number of atoms per species
      for(int is=0; is<NSP; is++) {
        inFile >> NA[is];
        if(NAX < NA[is]) NAX = NA[is];
      }
      printf("NAX = %d\n", NAX);
      // Allocate array for atomic coordinates
      TAU = (double*)malloc(3*NAX*NSP*sizeof(double));
      for(int is=0; is<NSP; is++) {
        inFile >> speciesName;
        SPECIES.push_back(speciesName);
        printf("SPECIES[%d] = %s\n", is, &(SPECIES[is][0]));
        // Read LMAX and LLOC
        getline(inFile,line);
        getline(inFile,line);
        printf("line = %s\n", &line[0]);
        sscanf(&line[0],"%d %d", &LMAX[is], &LLOC[is]);
        printf("LMAX[%d]=%d, LLOC[%d]=%d\n", is, LMAX[is], is, LLOC[is]);
        // Read coord type
        inFile >> coord_type;
        printf("Coord_type = %s\n", &coord_type[0]);
        getline(inFile,line);
        for(int ia=0; ia<NA[is]; ia++) {
          getline(inFile,line);
          sscanf(&line[0],"%lf %lf %lf", &TAU[IDX3F(1,ia,is,3,NAX)],
              &TAU[IDX3F(2,ia,is,3,NAX)], &TAU[IDX3F(3,ia,is,3,NAX)]);
          printf("%f %f %f\n",TAU[IDX3F(1,ia,is,3,NAX)],
              TAU[IDX3F(2,ia,is,3,NAX)], TAU[IDX3F(3,ia,is,3,NAX)]);
        }
      }
    } else if(line == "") {
      // Pass
    }
    else {
      printf("Line = %s\n", &line[0]);
    }
  }

  inFile.close();
  printf("Program ended normally\n");
  return 0;
}

