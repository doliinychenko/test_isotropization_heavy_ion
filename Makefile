# The compiler
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbacktrace -Wall
#FCFLAGS = -O3

PROJDIR = .
MDIR = ${PROJDIR}/modules
OBJ = ${PROJDIR}/obj

# List of executables to be built within the package
PROGRAMS = ${MDIR}/eispack.o ${MDIR}/kinds.o ${MDIR}/cla.o ${MDIR}/Land_Eck.o thermalization

all: $(PROGRAMS)

thermalization: thermalization.o eispack.o Land_Eck.o kinds.o cla.o
	$(FC) $(FCFLAGS) -o $@ $^ -I${OBJ}

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ -I${OBJ}

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -J${OBJ} -I${OBJ}

clean:
	rm -f *.o *.mod *.MOD



