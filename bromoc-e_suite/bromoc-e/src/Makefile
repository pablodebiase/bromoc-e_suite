#! /usr/bin/make
#
# Main program and default target

OBJ    = bromoc
OBJI   = bromoc-intel
all : $(OBJ)

# System dependent compilation commands

FCOMPL = gfortran
#FFLAGC = -C -fno-align-commons
#Double Precision Real
FFLAGC = -fno-align-commons -fdefault-real-8 
OPT    = -O3 -march=native -ffast-math -fno-protect-parens -ffree-line-length-none

# Implicit rules

.SUFFIXES: .vec .f .c .s .o .a .lib .f90

.f90.o: 
	$(FCOMPL) -c $(FFLAGC) $(OPT) -o $@ $<

.f.o:
	$(FCOMPL) -c $(FFLAGC) $(OPT) -o $@ $<

.s.o:
	as -o $@ $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) $(OPT) -o $@ $<

.lib.a:
	make -f $(<:.lib=.Makefile) $@
	make -f $(<:.lib=.Makefile) clean

# Vectorization of selected files where available
# Object files

OBJECTSMOD = 	ioxmod.o charfuncmod.o listmod.o constamod.o strtoolsmod.o errormod.o nucleotmod.o\
		stdiomod.o apfmod.o gsbpmod.o efpmod.o grandmod.o splinemod.o sevalmod.o charmmmod.o

OBJECTS = 	charman.o cputime.o error.o discretize.o splinepot.o cosmphi2.o sinmphi2.o getopt.o\
		ioinit.o lualloc.o lunalloc.o m3.o dm3.o rpowerl2.o cross_product.o factori.o find.o\
		angles.o bbonds.o count.o proxdiff.o dihedral.o countions.o repwalls.o nucenergy.o\
		fixcoor.o go_qq.o grand.o header.o mapsplot.o help.o insert.o membrane.o metropolis.o\
		move.o native_structure.o rfpar.o feedseed.o chrden.o readphi.o rgauss.o energy.o\
		staefield.o shell_simul.o simul1.o simul.o spline.o staticf1.o switch1.o switch2.o\
		switch3.o vdwgd1spln.o vdwgd1trln.o wrttraj.o angbond.o angdih.o polynom.o readpsf.o\
		readcharmm.o fatnam.o parinteract.o dynamics0nuc.o dynamics0.o dynamics1.o dynamics2.o\
		dynamics3.o charmmenergy.o

ORANDOM =	random.o
ORANDOMI =	random-i.o

# Clean objects

clean:
	@rm -f core $(OBJ) *.mod *.o

cleano:
	@rm -f core *.o *.mod

# GFORTRAN Single Precision
single: FFLAGC = -fno-align-commons
single: $(OBJ)

# GFORTRAN Debug
debug: OPT = -Wall -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fbacktrace -g -fdump-core -fcheck=all
debug: $(OBJ)	

# INTEL Single Precision
intelsp: FCOMPL = ifort
intelsp: FFLAGC = -real-size 32 -integer-size 32 -heap-arrays 0
intelsp: OPT = -ip -O3 -xHost
intelsp: ORANDOM = random-i.o
intelsp: $(OBJI)

# INTEL Debug
inteldbg: FCOMPL = ifort
inteldbg: FFLAGC = -real-size 64 -integer-size 32 -heap-arrays 0
inteldbg: OPT = -warn all -O0 -g -traceback -fpe:0 -check all -fp-stack-check -debug all
inteldbg: ORANDOM = random-i.o
inteldbg: $(OBJI)

# INTEL Optimized double precision
intel: FCOMPL = ifort
intel: FFLAGC = -real-size 64 -integer-size 32 -heap-arrays 0
intel: OPT = -ip -O3 -xHost
intel: ORANDOM = random-i.o
intel: $(OBJI)

# Main program target

$(OBJ):	$(OBJECTSMOD) $(OBJECTS) $(ORANDOM)
	$(FCOMPL) $(FFLAGC) $(OPT) -o $(OBJ) $(OBJECTSMOD) $(OBJECTS) $(ORANDOM)

$(OBJI): $(OBJECTSMOD) $(OBJECTS) $(ORANDOMI)
	$(FCOMPL) $(FFLAGC) $(OPT) -o $(OBJ) $(OBJECTSMOD) $(OBJECTS) $(ORANDOMI)

install:
	@mv -f $(OBJ $(OBJI)) ../bin/
	@if [ ! -e ../../bin/$(OBJ) ] ; then ln -s ../cg-gcmc-bd/bin/$(OBJ) ../../bin/$(OBJ); fi
	@if [ ! -e ../../bin/$(OBJI) ] ; then ln -s ../cg-gcmc-bd/bin/$(OBJI) ../../bin/$(OBJI); fi
 
