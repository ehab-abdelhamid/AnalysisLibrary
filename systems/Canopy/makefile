GXX=mpic++
FLAGS=-std=c++11

all: canopy_mpi

canopy_mpi: canopy_mpi.cpp
	${GXX} ${FLAGS} canopy_mpi.cpp -o canopy_mpi

run: all
	mpiexec -n 2 ./canopy_mpi -n 10

clean:
	rm canopy_mpi

.PHONY: all clean run
