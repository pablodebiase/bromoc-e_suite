#! /usr/bin/make
#
# Main program and default target

OBJ    = bromoc-e
INTEL  = _intel
OMP    = _omp
all: $(OBJ)

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
		angles.o bbonds.o count.o proxdiff.o dihedral.o countions.o repwalls.o fixcoor.o\
		go_qq.o grand.o header.o mapsplot.o help.o insert.o membrane.o metropolis.o move.o\
		native_structure.o rfpar.o feedseed.o chrden.o readphi.o rgauss.o staefield.o fatnam.o\
		shell_simul.o simul1.o simul.o spline.o switch1.o switch2.o switch3.o wrttraj.o\
		angbond.o angdih.o polynom.o readpsf.o readcharmm.o parinteract.o charmmenergy.o
OBJSER = 	dynamics0nuc.o dynamics0.o dynamics1.o dynamics2.o dynamics3.o energy.o staticf1.o\
		vdwgd1trln.o vdwgd1spln.f90 rfparion.f90 nucenergy.f90
PARALLEL =	dynamics0nuc-omp.f90 dynamics0-omp.f90 dynamics1-omp.f90 dynamics2-omp.f90 dynamics3-omp.f90\
		energy-omp.f90 staticf1-omp.f90 vdwgd1trln-omp.f90 vdwgd1spln-omp.f90 rfparion-omp.f90 nucenergy-omp.f90
ORANDOM =	random.o
ORANDOMI =	random-i.o

# Clean objects

clean:
	@rm -f core $(OBJ) $(OBJ)$(INTEL) $(OBJ)$(OMP) $(OBJ)$(INTEL)$(OMP) *.mod *.o

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
intelsp: $(OBJ)$(INTEL)

# INTEL Debug
inteldbg: FCOMPL = ifort
inteldbg: FFLAGC = -real-size 64 -integer-size 32 -heap-arrays 0
inteldbg: OPT = -warn all -O0 -g -traceback -fpe:0 -check all -fp-stack-check -debug all
inteldbg: ORANDOM = random-i.o
inteldbg: $(OBJ)$(INTEL)

# INTEL Optimized double precision
intel: FCOMPL = ifort
intel: FFLAGC = -real-size 64 -integer-size 32 -heap-arrays 0
intel: OPT = -ip -O3 -xHost
intel: ORANDOM = random-i.o
intel: $(OBJ)$(INTEL)

# Main program target

$(OBJ):	$(OBJECTSMOD) $(OBJECTS) $(ORANDOM) $(OBJSER)
	$(FCOMPL) $(FFLAGC) $(OPT) -o $(OBJ) $(OBJECTSMOD) $(OBJECTS) $(ORANDOM) $(OBJSER)
	$(FCOMPL) $(FFLAGC) $(OPT) -fopenmp -o $(OBJ)$(OMP) $(OBJECTSMOD) $(OBJECTS) $(ORANDOM) $(PARALLEL)

$(OBJ)$(INTEL): $(OBJECTSMOD) $(OBJECTS) $(ORANDOMI) $(OBJSER)
	$(FCOMPL) $(FFLAGC) $(OPT) -o $(OBJ)$(INTEL) $(OBJECTSMOD) $(OBJECTS) $(ORANDOMI) $(OBJSER)
	$(FCOMPL) $(FFLAGC) $(OPT) -openmp -o $(OBJ)$(INTEL)$(OMP) $(OBJECTSMOD) $(OBJECTS) $(ORANDOMI) $(PARALLEL)

install:
	@
	if [ ! -d ../bin ]; then mkdir ../bin; fi
	if [ ! -d ../../bin ]; then mkdir ../../bin; fi
	if [ -e $(OBJ) ]; then mv -f $(OBJ) ../bin/; fi
	if [ -e $(OBJ)$(OMP) ]; then mv -f $(OBJ)$(OMP) ../bin/; fi
	if [ -e $(OBJ)$(INTEL) ]; then mv -f $(OBJ)$(INTEL) ../bin/; fi
	if [ -e $(OBJ)$(INTEL)$(OMP) ]; then mv -f $(OBJ)$(INTEL)$(OMP) ../bin/; fi
	if [ ! -e ../../bin/$(OBJ) -a -e ../bromoc-e/bin/$(OBJ) ] ; then ln -s ../bromoc-e/bin/$(OBJ) ../../bin/$(OBJ); fi
	if [ ! -e ../../bin/$(OBJ)$(OMP) -a -e ../bromoc-e/bin/$(OBJ)$(OMP) ] ; then ln -s ../bromoc-e/bin/$(OBJ)$(OMP) ../../bin/$(OBJ)$(OMP); fi
	if [ ! -e ../../bin/$(OBJ)$(INTEL) -a -e ../bromoc-e/bin/$(OBJ)$(INTEL) ] ; then ln -s ../bromoc-e/bin/$(OBJ)$(INTEL) ../../bin/$(OBJ)$(INTEL); fi
	if [ ! -e ../../bin/$(OBJ)$(INTEL)$(OMP) -a -e ../bromoc-e/bin/$(OBJ)$(INTEL)$(OMP) ] ; then ln -s ../bromoc-e/bin/$(OBJ)$(INTEL)$(OMP) ../../bin/$(OBJ)$(INTEL)$(OMP); fi

uninstall:
	@
	rm -f ../bin/$(OBJ)
	rm -f ../bin/$(OBJ)$(OMP)
	rm -f ../bin/$(OBJ)$(INTEL)
	rm -f ../bin/$(OBJ)$(INTEL)$(OMP)
	rm -f ../../bin/$(OBJ)
	rm -f ../../bin/$(OBJ)$(OMP)
	rm -f ../../bin/$(OBJ)$(INTEL)
	rm -f ../../bin/$(OBJ)$(INTEL)$(OMP)
