include platform/make.inc.x86_64_intel

SRC = \
  PrintArray.cpp  global_variables.cpp timing.cpp \
  ReadInputFile.cpp Setup.cpp dsum.cpp \
  PrintGlobals.cpp SumArray.cpp Print3DArray.cpp \
  determinant3x3.cpp CalcRecipLatt.cpp \
  GVectorsGen0.cpp ReadPseudopotentials.cpp \
  AllocateArrays.cpp InitOcc.cpp GVectorsGen.cpp GkVectorsGen.cpp \
  sort1.cpp PrepareSCF.cpp PhaseFactor.cpp matinv3x3.cpp SortGVectors.cpp \
  StructureFactor.cpp simpson.cpp FormFactor.cpp \
  NLFormFactor.cpp SCFLoop.cpp FormFactorAtomic.cpp fft_fftw3.cpp \
  NormalizeRho.cpp EvalEnergy.cpp EwaldEnergy.cpp xc_PZ.cpp \
  PrintSystemInfo.cpp ApplyHam.cpp ApplyHam_block.cpp \
  eig_zheevd.cpp eig_zhegv.cpp ortho_qr.cpp \
  EvalRhoPsi.cpp EvalENL.cpp GenPrec.cpp SimpleRhoMix.cpp \
  AndersonRhoMix.cpp util_random_cpu.cpp \
  BroydenRhoMix.cpp PulayRhoMix.cpp KerkerRhoMix.cpp \
  lobpcg_cpu.cpp Diagon_cpu.cpp davidson_cpu.cpp \
  lanczos.cpp ChebySCFLoop.cpp chebyfilt.cpp

OBJ = $(SRC:.cpp=.o) $(SRC:.f=.o)

# Suffix rules
.SUFFIXES: .o .cpp
.cpp.o:
	$(CXX) $(CXX_OPTS) -c $<

.SUFFIXES: .o .cu
.cu.o:
	$(NVCC) $(NVCC_OPTS) -c $<

.SUFFIXES: .o .f
.f.o:
	ifort -c $<

# Targets
lib: $(OBJ)
	ar rcs libmain.a *.o

main:
	$(CXX) $(CXX_OPTS) $(CXX_INC) main.cpp libmain.a $(LIBS) -o ffr_pspw_dft.x

clean:
	rm -rf *.o libmain.a

