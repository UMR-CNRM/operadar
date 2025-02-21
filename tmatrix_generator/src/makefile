# This is the Makefile for the T-matrix program

# Compiler settings are here
#include ./makefile_include

# Objects to compile
OBJECTS = ./Tmatrix.f

# belenos :
#LIBS = -L ../LIB/LAPACK  -llapack_LINUX -lblas_LINUX

# tests pxprecip9 (ubuntu)
#LIBS = -L ./LIB/lapack-3.10.0  -lliblapack.a -llibrefblas.a
LIBS = -L ./lib/lapack-3.11  -llapack -lrefblas
#LIBS = -L /usr/lib/x86_64-linux-gnu/lapack -lliblapack

CF = gfortran

CFLAGS = -pipe -Wall 


all: Tmat

Tmat: $(OBJECTS) 
	$(CF) $(CFLAGS) $(OBJECTS) $(LIBS) -w -o ./Tmat 

clean:
	rm -f *.o Tmat
