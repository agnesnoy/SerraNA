# Set up the flags you want in FFLAGS and compiler in FC

# flags for GNU 
FFLAGS = -O0 #-g -Wall -fcheck=all
FC = gfortran

# Suffixes rules to control how .o and .mod files are treated
.SUFFIXES:
.SUFFIXES: .f90 .o .mod 
.f90.o:
	$(FC) $(FFLAGS)  -c $< 
.f90.mod:
	$(FC) $(FFLAGS)  -c $< 

# Objects (OBJ) to be used

OBJ = SerraNA.f90 parms.o io_mod.o functions_mod.o
OBJ2 = Analysis.f90 parms.o io_mod.o functions_mod.o
OBJ3 = Extract.f90 parms.o io_mod.o functions_mod.o

#Compiles SerraLINE
all:	SerraNA Analysis Extract

SerraNA:	$(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@

Analysis:	$(OBJ2)
	$(FC) $(FFLAGS) $(OBJ2) -o $@

Extract:	$(OBJ3)
	$(FC) $(FFLAGS) $(OBJ3) -o $@

parms.o parms.mod: parms.f90
	$(FC) $(FFLAGS) -c parms.f90

clean:
	rm *.o *.mod SerraNA Analysis Extract

