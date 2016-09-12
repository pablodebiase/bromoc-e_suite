#! /usr/bin/make
#
# Main program and default target

OBJ    = bromoc
all : $(OBJ)

# System dependent compilation commands

FCOMPL = gfortran
#FFLAGC = -C -fno-align-commons
#Double Precision Real
FFLAGC = -fno-align-commons -fdefault-real-8 
CCOMPL = gcc
CFLAGC =
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
		move.o native_structure.o rfpar.o feedseed.o random.o chrden.o readphi.o rgauss.o\
		staefield.o shell_simul.o simul1.o simul.o spline.o staticf1.o switch1.o switch2.o\
		switch3.o vdwgd1spln.o vdwgd1trln.o wrttraj.o angbond.o angdih.o polynom.o readpsf.o\
		readcharmm.o fatnam.o parinteract.o dynamics0nuc.o dynamics0.o dynamics1.o dynamics2.o\
		dynamics3.o charmmenergy.o energy.o

# Clean objects

clean:
	@rm -f core $(OBJECTSMOD) $(OBJECTS) $(OBJ) *.mod

cleano:
	@rm -f core *.o *.mod

debug: OPT = -Wall -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fbacktrace -g -fdump-core -fcheck=all
debug: $(OBJ)	

# Main program target

$(OBJ):	$(OBJECTSMOD) $(OBJECTS)
	$(FCOMPL) $(FFLAGC) $(OPT) -o $(OBJ) $(OBJECTSMOD) $(OBJECTS)

install:
	@mv $(OBJ) ../bin/
	@if [ ! -e ../../bin/$(OBJ) ] ; then ln -s ../cg-gcmc-bd/bin/$(OBJ) ../../bin/$(OBJ); fi
 
