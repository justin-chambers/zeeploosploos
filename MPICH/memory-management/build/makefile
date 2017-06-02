#General config environment variables...
CC=mpiCC
CXXFLAGS=-Wall
LIBS=-lpmi
FILE=../src/mm.cpp
#INCLUDES=-I../include
OUTPUT=../bin/mm

# MPI-Specific environment variables...
MPI_HOSTS=slurmless/mpi_hosts

all: $(OUTPUT)

$(OUTPUT): $(FILE)
	$(CC) $(CXXFLAGS) -o $(OUTPUT) $(FILE) $(LIBS)

run: $(OUTPUT)
	mpirun -n 4 $(OUTPUT)

clean:
	-@if rm -rf $(OUTPUT) || true; then echo "Cleaned"; else echo "Nothing to clean"; fi
