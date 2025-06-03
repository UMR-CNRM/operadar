# This is the Makefile for the T-matrix program

# Compiler settings are here
#include ./makefile_include

# Objects to compile
OBJECTS = ./Tmatrix.f

# Belenos :
#LIBS = -L ../LIB/LAPACK  -llapack_LINUX -lblas_LINUX

# Local (ubuntu)
# Should be something like : LIBS = -L ./my_lapack_path -llapack -lrefblas
LIBS = -L ./libFortran/lapack -llapack -lrefblas

CF = gfortran

CFLAGS = -pipe -Wall 


all: Tmat

Tmat: $(OBJECTS) 
	$(CF) $(CFLAGS) $(OBJECTS) $(LIBS) -w -o ./Tmat 

clean:
	rm -f *.o Tmat
