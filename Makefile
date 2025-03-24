## Created on 2024-05-15 at 15:21:09 CEST by David Gaspard <david.gaspard@espci.fr>
## Makefile of the Fortran program RecurGreen aiming at solving the wave equation in a 2D disordered waveguide.
##########################
## FILE PARAMETERS
##########################
SRCDIR   = src
BINDIR   = bin

## Sources in dependency 'use' order (from low to high level):
MODULES  = base_utils.f90 \
           random_normal.f90 \
           waveguide.f90

MODLIST = $(MODULES:%.f90=$(SRCDIR)/%.f90)
BINLIST = $(MODULES:%.f90=$(BINDIR)/%.o)

##########################
## COMPILATION OPTIONS
##########################
FORT    = gfortran
STD     = -std=f2008
OMP     = -fopenmp
DEBUG   = -g -O0 -fbacktrace -fcheck=all
## -fsanitize=address
FFLAGS  = -J$(BINDIR) $(OMP) -Wall -Wno-tabs
## -O3 -march=native
LIBS    = -lblas -llapack

all: recurgreen plothisto

recurgreen: $(BINLIST) bin/main.o
	$(FORT) $(FFLAGS) $^ $(LIBS) -o $@

plothisto: bin/base_utils.o bin/plot_histogram.o
	$(FORT) $(FFLAGS) $^ $(LIBS) -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.f90
	$(FORT) $(FFLAGS) -c $< -o $@

##########################
## TEST BUILD FUNCTIONS
##########################
TESTSRCLIST = $(shell find $(SRCDIR) -name "*.test.f90")
TESTEXELIST = $(TESTSRCLIST:$(SRCDIR)/%.f90=%)

test: $(BINLIST) $(TESTEXELIST)

%.test: $(BINLIST) $(BINDIR)/%.test.o
	$(FORT) $(FFLAGS) $^ $(LIBS) -o $@

##################################
## CLEAN ALL BUILDS AND TESTS
##################################
clean:
	rm -rfv recurgreen plothisto $(BINDIR)/*.o $(BINDIR)/*.mod $(TESTEXELIST)

###### END OF FILE ######
