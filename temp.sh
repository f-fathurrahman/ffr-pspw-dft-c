g++ -D_USE_CUDA -I /home/marisa1709efefer/mysoftwares/cuda/include/ main.cpp test_ApplyHam_cuda.cpp libmain.a -L/opt/intel/Compiler/11.1/038/mkl/lib/32/ -lmkl_intel -lmkl_core -lmkl_sequential -lfftw3 -L/home/marisa1709efefer/mysoftwares/cuda/lib -lcudartemu -lcufftemu -lcublasemu

