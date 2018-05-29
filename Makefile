# $HeadURL$
# $Id$
FFLAGS= -g -O3 -fdefault-real-8 -I.
#programs=weno3.x weno5.x
programs=weno5.x
common=prms.o conv.o rhs.o io.o fluxes.o ic.o matrices.o maths.o viscous.o

all: $(programs)

main3.o: FFLAGS += -DWENOORDER=3
main3.o: main.F90
	$(FC) $(FFLAGS) -c -o $@ $<

weno3.x: main3.o reconstruct3.o $(common)
	$(LD) -o $@ $^

main5.o: FFLAGS += -DWENOORDER=5
main5.o: main.F90
	$(FC) $(FFLAGS) -c -o $@ $<

weno5.x: main5.o reconstruct5.o $(common)
	$(LD) -o $@ $^

# Module dependencies
prms.F90:		
io.F90:               prms.mod matrices.mod
rhs.F90:              prms.mod
maths.F90:            prms.mod
fluxes.F90:           prms.mod
mv.F90:               prms.mod
main.F90:             prms.mod conv.mod rhs.mod matrices.mod maths.mod io.mod fluxes.mod
ic.F90:               prms.mod io.mod conv.mod
reconstruct3.F90:     prms.mod
reconstruct5.F90:     prms.mod
viscous.F90:          prms.mod
matrices.F90:         prms.mod conv.mod maths.mod

clean:
	@rm -fv  *.mod *.o *.x *__genmod.f90 *__genmod.mod

FC=gfortran
LD=${FC} 
RANLIB=touch
AR=ar r

.PHONY: clean docs

.SUFFIXES:
.SUFFIXES: .f .o
.SUFFIXES: .F90 .o
.SUFFIXES: .F90 .mod

.f.o:
	$(FC) $(FFLAGS) -c $<

.F90.o:
	$(FC) $(FFLAGS) -c $<

.F90.mod:
	$(FC) $(FFLAGS) -c $<
