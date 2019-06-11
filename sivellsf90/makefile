#
#  makefile to compile Contur's code by Sivells AEDC-TR-78-63(1978)
#
FC=ifort
FFLAGS= -ansi-alias -diag-enable=all -I. -f77rtl -fast -fno-common \
-noextend-source -no-fast-transcendentals -fixed -fp-model precise \
-warn all -check all

modules=\
kinddefine.o\
mod_cline.o\
mod_coef.o\
mod_contr.o\
mod_fg.o\
mod_gg.o\
mod_httr.o\
mod_coord.o\
mod_corr.o\
mod_prop.o\
mod_param.o\
mod_jack.o\
mod_troat.o\
mod_work.o\
axial.o \
bound.o \
conic.o \
cubic.o \
fmv.o \
fvdge.o \
heat.o \
neo.o \
ofeld.o \
perfc.o \
plate.o \
scond.o \
sorce.o \
splind.o \
toric.o \
trans.o \
twixt.o \
xyz.o

default: ${modules} main.o
	$(FC) $(FFLAGS) ${modules} main.o -o main.exe

.SUFFIXES: .f .o

clean:
	/bin/rm main.o main.exe ${modules} *.mod

