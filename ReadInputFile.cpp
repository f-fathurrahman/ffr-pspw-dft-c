// eFeFeR (20910015), October 2011

#include "common_pspw_cuda.h"

void ReadInputFile()
{
  printf("Reading input file ....");

  string inputFileName = "INPUT";
  string line;
  string speciesName;
  string coord_type;
  int ia, is, ik;

  ifstream inFile(&inputFileName[0]);
  if(!inFile.is_open()) {
    printf("ERROR opening file: %s\n", &inputFileName[0]);
    abort();
  }

// Set up default values
  NEMPTY = 5;
  TMETAL = false;
  EKT = 0.004;
  MAX_SCF_ITER = 30;
  SCF_CONV_CRIT = 1.e-6;
  RHO_MIX_FAC = 0.3;
  KERKER_Q0 = 0.5;
  BETAMAX = 0.5;
  MIXDIM = 4;
  IDIAG = 0;
  MAX_DIAG_ITER = 30;
  ISCF = 0;
  CHEBYPOL_DEGREE = 10;

  while(!inFile.eof()) {
    getline(inFile,line);
    //
    if(line == "ecut") {
      inFile >> ECUT;
    }
    //
    else if(line == "latvec") {
      getline(inFile,line);
      sscanf(&line[0],"%lf %lf %lf", &A1[0], &A1[1], &A1[2]);
      getline(inFile,line);
      sscanf(&line[0],"%lf %lf %lf", &A2[0], &A2[1], &A2[2]);
      getline(inFile,line);
      sscanf(&line[0],"%lf %lf %lf", &A3[0], &A3[1], &A3[2]);
    }
    //
    else if(line == "atoms") {
      getline(inFile,line);
      sscanf(&line[0],"%d", &NSP);
      //printf("NSP = %d\n", NSP);
      // Allocate memory for several variables that depend on NSP
      NA = (int*)malloc(NSP*sizeof(int));
      L_MAX = (int*)malloc(NSP*sizeof(int));
      L_LOC = (int*)malloc(NSP*sizeof(int));
      RGAUSS = (double*)malloc(NSP*sizeof(double));
      for(is=1; is<=NSP; is++) RGAUSS[is-1] = 1.0; // default value
      NAX = 0;
      // Read number of atoms per species
      for(is=1; is<=NSP; is++) {
        inFile >> NA[is-1];
        if(NAX < NA[is-1]) NAX = NA[is-1];
      }
      // Allocate array for atomic coordinates
      TAU = (double*)malloc(3*NAX*NSP*sizeof(double));
      for(int i=0; i<3*NAX*NSP; i++) TAU[i] = 0.0;
      for(is=1; is<=NSP; is++) {
        inFile >> speciesName;
        ATOMTYP.push_back(speciesName);
        // Read LMAX and LLOC
        getline(inFile,line);
        getline(inFile,line);
        sscanf(&line[0],"%d %d", &L_MAX[is-1], &L_LOC[is-1]);
        // Read coord type
        inFile >> coord_type;
        getline(inFile,line);
        for(ia=1; ia<=NA[is-1]; ia++) {
          getline(inFile,line);
          sscanf(&line[0],"%lf %lf %lf", &TAU[IDX3F(1,ia,is,3,NAX)],
              &TAU[IDX3F(2,ia,is,3,NAX)], &TAU[IDX3F(3,ia,is,3,NAX)]);
        }
      }
    }
    //
    else if(line == "k_points") {
      getline(inFile,line);
      sscanf(&line[0],"%d",&NKPT);
      WKPT = (double*)malloc(NKPT*sizeof(double));
      XK = (double*)malloc(3*NKPT*sizeof(double));
      for(ik=1; ik<=NKPT; ik++) {
        getline(inFile,line);
        sscanf(&line[0],"%lf %lf %lf %lf\n",&XK[IDX2F(1,ik,3)],
            &XK[IDX2F(2,ik,3)], &XK[IDX2F(3,ik,3)], &WKPT[ik-1]);
      }
    }
    //
    else if(line == "rho_mix_fac") {
      getline(inFile,line);
      sscanf(&line[0], "%lf", &RHO_MIX_FAC);
    }
    //
    else if(line == "kerker_q0") {
      getline(inFile,line);
      sscanf(&line[0], "%lf", &KERKER_Q0);
    }
    //
    else if(line == "nempty") {
      getline(inFile,line);
      sscanf(&line[0], "%d", &NEMPTY);
    }
    //
    else if(line == "tmetal") {
      getline(inFile,line);
      if(line == " .true.") { // TODO: trim the string
        TMETAL = true;
      }
    }
    //
    else if(line == "ekt") {
      getline(inFile,line);
      sscanf(&line[0], "%lf", &EKT);
    }
    //
    else if(line == "imix") {
      getline(inFile,line);
      sscanf(&line[0], "%d", &IMIX);
    }
    else if(line == "mixdim") {
      getline(inFile,line);
      sscanf(&line[0], "%d", &MIXDIM);
    }
    //
    else if(line == "max_scf_iter") {
      getline(inFile,line);
      sscanf(&line[0], "%d", &MAX_SCF_ITER);
    }
    //
    else if(line == "scf_conv_crit") {
      getline(inFile,line);
      sscanf(&line[0], "%lf", &SCF_CONV_CRIT);
    }
    //
    else if(line == "use_cuda") {
      getline(inFile,line);
      if(line == " .true.") { // TODO: trim the string
        USE_CUDA = true;
      }
    }
    //
    else if(line == "fft_grid") {
      getline(inFile,line);
      USER_FFT_GRID = true;
      sscanf(&line[0], "%d %d %d", &NR1,&NR2,&NR3);
      NNR = NR1*NR2*NR3;
    }
    else if(line == "idiag") {
      getline(inFile,line);
      sscanf(&line[0], "%d", &IDIAG);
    }
    else if(line == "max_diag_iter") {
      getline(inFile,line);
      sscanf(&line[0], "%d", &MAX_DIAG_ITER);
      USER_MAX_DIAG_ITER = true;
    }
    else if(line == "iscf") {
      getline(inFile,line);
      sscanf(&line[0], "%d", &ISCF);
    }
    else if(line == "chebypol_degree") {
      getline(inFile,line);
      sscanf(&line[0], "%d", &CHEBYPOL_DEGREE);
    }
    else if(line == "diag_verbose") {
      getline(inFile,line);
      if(line == " .true.") { // TODO: trim the string
        DIAG_VERBOSE = true;
      }
    }
    // Default
    else if(line == "") {
      // Pass
    }
    else {
      printf("***Unknown line = %s\n", &line[0]);
      abort();
    }
  }

  inFile.close();

  if(coord_type=="angstrom") {
    for(is=1; is<=NSP; is++) {
      for(ia=1; ia<=NA[is-1]; ia++) { // convert to bohr
        TAU[IDX3F(1,ia,is,3,NAX)] = TAU[IDX3F(1,ia,is,3,NAX)]/0.529;
        TAU[IDX3F(2,ia,is,3,NAX)] = TAU[IDX3F(2,ia,is,3,NAX)]/0.529;
        TAU[IDX3F(3,ia,is,3,NAX)] = TAU[IDX3F(3,ia,is,3,NAX)]/0.529;
      }
    }
  }
  else if(coord_type=="crystal") { // only works for orthogonal lattice?
    for(is=1; is<=NSP; is++) {
      for(ia=1; ia<=NA[is-1]; ia++) { // convert to bohr
        TAU[IDX3F(1,ia,is,3,NAX)] = TAU[IDX3F(1,ia,is,3,NAX)]*A1[0] +
          TAU[IDX3F(2,ia,is,3,NAX)]*A2[0] + TAU[IDX3F(3,ia,is,3,NAX)]*A3[0];
        TAU[IDX3F(2,ia,is,3,NAX)] = TAU[IDX3F(1,ia,is,3,NAX)]*A1[1] + 
          TAU[IDX3F(2,ia,is,3,NAX)]*A2[1] + TAU[IDX3F(3,ia,is,3,NAX)]*A3[1];
        TAU[IDX3F(3,ia,is,3,NAX)] = TAU[IDX3F(1,ia,is,3,NAX)]*A1[2] + 
          TAU[IDX3F(2,ia,is,3,NAX)]*A2[2] + TAU[IDX3F(3,ia,is,3,NAX)]*A3[2];
        /*TAU[IDX3F(1,ia,is,3,NAX)] = TAU[IDX3F(1,ia,is,3,NAX)]*A1[0] +
          TAU[IDX3F(2,ia,is,3,NAX)]*A1[1] + TAU[IDX3F(3,ia,is,3,NAX)]*A1[2];
        TAU[IDX3F(2,ia,is,3,NAX)] = TAU[IDX3F(1,ia,is,3,NAX)]*A2[0] +
          TAU[IDX3F(2,ia,is,3,NAX)]*A2[1] + TAU[IDX3F(3,ia,is,3,NAX)]*A2[2];
        TAU[IDX3F(3,ia,is,3,NAX)] = TAU[IDX3F(1,ia,is,3,NAX)]*A3[0] +
          TAU[IDX3F(2,ia,is,3,NAX)]*A3[1] + TAU[IDX3F(3,ia,is,3,NAX)]*A3[2];*/
      }
    }
  }

  printf("...Done reading input file\n");
}

