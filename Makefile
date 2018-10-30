.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h
MAKE     = make
F90=gfortran 
#FLAGS= -c -O3 -ffree-line-length-none
FLAGS= -c -ffree-line-length-none -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -fbounds-check -ffpe-trap=zero
#FLAGS=-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

MUMPS    = ../mumps/lib
MUMPS2   = ../mumps/libseq
PORD     = ../mumps/PORD

LIBS = \
-L$(MUMPS) -ldmumps -lmumps_common\
-L$(MUMPS2) -lmpiseq \
-L$(PORD) -lpord \
-lpthread -lblas 


OBJECTS2D =\
blas_routines.o\
linpack_d.o\
output_for_paraview.o\
solve_uzawa1.o\
solve_uzawa2.o\
solve_uzawa3.o\
solve_linpack.o\
int_to_char.o\
inverse_icon.o\
elemental_to_nodal.o\
solve_uzawa_MUMPS.o\
equation_of_state.o\
simplefem_thermal_streamlined.o 

.f.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f
.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

code:	$(OBJECTS2D)
	$(F90) $(OPTIONS) $(OBJECTS2D) $(LIBS) -o simplefem

clean: 
	rm -f *.o
	rm -f *.dat
	rm OUT/*.dat
	rm OUT/*.vtu


