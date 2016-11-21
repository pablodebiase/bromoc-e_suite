#! /usr/bin/make
#
# Main program and default target

OBJ    = imc-macro
PAR    = $(OBJ)-par
PAR2   = $(OBJ)-par-2
OBJ2   = imc-e
all : $(OBJ) $(PAR) $(PAR2) $(OBJ2)

# System dependent compilation commands

FCOMPL = gfortran #gnu
#FFLAGC = -fdefault-real-8 # gnu real*8 
FFLAGC =  # gnu real*4 (default)
CCOMPL = gcc # gnu
CFLAGC =
OPT = -march=native -ffast-math -fno-protect-parens -ffree-line-length-none -O3 -fopenmp #-Wall
SYSDEP =
UTILS  = 
LDFLAG = 

# Implicit rules

.SUFFIXES: .vec .f .c .s .o .a .lib .f90

.f.o:
	$(FCOMPL) -c $(FFLAGC) $(OPT) -o $@ $<

.f90.o:
	$(FCOMPL) -c $(FFLAGC) $(OPT) -o $@ $<

.s.o:
	as -o $@ $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) $(OPT) -o $@ $<

.lib.a:
	make -f $(<:.lib=.Makefile) $@
	make -f $(<:.lib=.Makefile) clean

# Object files

OBJECTS = xerbla.o lsame.o ieeeck.o iparmq.o sger.o slaswp.o sswap.o ilaenv.o isamax.o sgemm.o slamch.o sscal.o strsm.o sgetf2.o sgetrs.o sgetrf.o sgesv.o imc-macro.o

OBJECTS1 = xerbla.o lsame.o ieeeck.o iparmq.o sger.o slaswp.o sswap.o ilaenv.o isamax.o sgemm.o slamch.o sscal.o strsm.o sgetf2.o sgetrs.o sgetrf.o sgesv.o imc-macro-par.o

OBJECTS1P = xerbla.o lsame.o ieeeck.o iparmq.o sger.o slaswp.o sswap.o ilaenv.o isamax.o sgemm.o slamch.o sscal.o strsm.o sgetf2.o sgetrs.o sgetrf.o sgesv.o imc-macro-par-2.o

OBJECTS2P = xerbla.o lsame.o ieeeck.o iparmq.o sger.o slaswp.o sswap.o ilaenv.o isamax.o sgemm.o slamch.o sscal.o strsm.o sgetf2.o sgetrs.o sgetrf.o sgesv.o mathmod.o imc-e.o

#IFILES = 

#IMC2XYZ

OBJECTS2 = imc2xyz.o

#CONVERTRST

OBJECTS3 = convertrst.o

# Clean OBJects

clean:
	@rm -f $(OBJECTS) $(OBJECTS1) $(OBJECTS1P) $(OBJECTS2P) $(OBJECTS2) $(OBJECTS3) $(OBJ) $(PAR) $(PAR2) $(OBJ2) *.mod imc2xyz convertrst

# Main program target

debuge: OPT = -Wall -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fbacktrace -g -fdump-core -fcheck=all -fopenmp
debuge: $(OBJ2)

$(OBJ): $(SYSDEP) $(VECTOR) $(OBJECTS) $(FFILES) $(UTILS)
	$(FCOMPL) $(FFLAGC) $(LDFLAG) $(OPT) -o $(OBJ) $(OBJECTS) $(SYSDEP) $(UTILS)
	@rm -f invmc.mod math.mod

$(PAR): $(SYSDEP) $(VECTOR) $(OBJECTS1) $(FFILES) $(UTILS)
	$(FCOMPL) $(FFLAGC) $(LDFLAG) $(OPT) -o $(PAR) $(OBJECTS1) $(SYSDEP) $(UTILS)
	@rm -f invmc.mod math.mod

$(PAR2):$(SYSDEP) $(VECTOR) $(OBJECTS1P) $(FFILES) $(UTILS)
	$(FCOMPL) $(FFLAGC) $(LDFLAG) $(OPT) -o $(PAR2) $(OBJECTS1P) $(SYSDEP) $(UTILS)
	@rm -f invmc.mod math.mod

$(OBJ2):$(SYSDEP) $(VECTOR) $(OBJECTS2P) $(FFILES) $(UTILS)
	$(FCOMPL) $(FFLAGC) $(LDFLAG) $(OPT) -o $(OBJ2) $(OBJECTS2P) $(SYSDEP) $(UTILS)
	@rm -f invmc.mod math.mod mathmod.mod

imc2xyz: $(SYSDEP) $(VECTOR) $(OBJECTS2) $(FFILES) $(UTILS)
	 $(FCOMPL) $(FFLAGC) $(LDFLAG) $(OPT) -o imc2xyz $(OBJECTS2) $(SYSDEP) $(UTILS)

convertrst: $(SYSDEP) $(VECTOR) $(OBJECTS3) $(FFILES) $(UTILS)
	    $(FCOMPL) $(FFLAGC) $(LDFLAG) $(OPT) -o convertrst $(OBJECTS3) $(SYSDEP) $(UTILS)
 
install:
	@ if [ -f $(OBJ) ]; then mv $(OBJ) ../bin/ ; fi
	@ if [ -f $(PAR) ]; then mv $(PAR) ../bin/ ; fi
	@ if [ -f $(PAR2) ]; then mv $(PAR2) ../bin/ ; fi
	@ if [ -f $(OBJ2) ]; then mv $(OBJ2) ../bin/ ; fi
	@ if [ ! -e ../../bin/$(OBJ) ]; then ln -s ../imc-macro/bin/$(OBJ) ../../bin/$(OBJ) ; fi
	@ if [ ! -e ../../bin/$(PAR) ]; then ln -s ../imc-macro/bin/$(PAR) ../../bin/$(PAR) ; fi
	@ if [ ! -e ../../bin/$(PAR2) ]; then ln -s ../imc-macro/bin/$(PAR2) ../../bin/$(PAR2) ; fi
	@ if [ ! -e ../../bin/$(OBJ2) ]; then ln -s ../imc-macro/bin/$(OBJ2) ../../bin/$(OBJ2) ; fi
