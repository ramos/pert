F90 = f90 -CB -check all
F90LD = f90 -CB -check all

all:
	@make geometry.o bases.o pert.out

bases.o: bases.f90
	$(F90) -c -I/home/alberto/fortran/lib bases.f90

geometry.o: geometry.f90
	$(F90) -c -I/home/alberto/fortran/lib geometry.f90

pert.o: pert.f90
	$(F90) -c -I/home/alberto/fortran/lib pert.f90

pert.out: bases.o geometry.o pert.o
	$(F90LD) -o pert.out bases.o geometry.o pert.o \
-L/home/alberto/fortran/lib -lf90

cleanall:
	rm *.o *.mod *.out
