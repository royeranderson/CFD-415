FC	= gfortran 

#FFLAGS= -m64 -g -ffpe-trap=invalid,zero -fbacktrace -fdump-core
FFLAGS= -m64 

#LDFLAGS = -m64 -g 
LDFLAGS = -m64 
TARGET_ARCH =

EXE   = JST

.SUFFIXES:
.SUFFIXES: .f90 .o

SRC  =					\
	eulersol.f90 matrices.f90 eulerBCs.f90 RK.f90 cellarea.f90 residuals.f90 dissapation.f90 initial.f90 ghostfill.f90 ghostcell.f90 		\

OBJ = $(SRC:.f90=.o)

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LIBS) -o $(EXE)

%.o  : %.f90 
	$(FC) $(FFLAGS) -c $<

clean: 
	rm -f *.mod *~ core
	rm -f *.o