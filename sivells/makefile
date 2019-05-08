#
#  makefile to compile Contur's code by Sivells AEDC-TR-78-63(1978)
#
FC=ifort
FFLAGS= -ansi-alias -diag-enable=all -I. -f77rtl -fast -fno-common \
-noextend-source -no-fast-transcendentals -fixed -fp-model precise 

modules=\
AXIAL.o \
BOUND.o \
CONIC.o \
CUBIC.o \
FMV.o \
FVDGE.o \
HEAT.o \
NEO.o \
OFELD.o \
OREZ.o \
PERFC.o \
PLATE.o \
SCOND.o \
SORCE.o \
SPLIND.o \
TRANS.o \
TORIC.o \
TWIXT.o \
XYZ.o

default: ${modules} MAIN.o
	$(FC) $(FFLAGS) ${modules} MAIN.o -o MAIN.exe

.SUFFIXES: .f .o

clean:
	/bin/rm MAIN.o MAIN.exe ${modules}

