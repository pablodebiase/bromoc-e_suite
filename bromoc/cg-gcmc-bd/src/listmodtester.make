gfortran -O0 -Wall -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fbacktrace -g -fdump-core -fcheck=all -fdefault-real-8 -o listmodtester listmod.f90 random.f90 listmodtester.f90
