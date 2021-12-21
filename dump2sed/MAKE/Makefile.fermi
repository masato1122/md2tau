#make

#--- For intel compiler
CC=mpiicpc

CFLAGS = -O3 -xHOST -std=c++11 \
	 -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64 \
         -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mpi/lib64 \
         -L/opt/intel/compilers_and_libraries/linux/lib/intel64 \
         -L/opt/intel/compilers_and_libraries/linux/mpi/lib64 \
         -L/opt/intel/mkl/lib/intel64 \
         -L/home/ohnishi/usr/local/lib

LFLAGS = -I/opt/intel/compilers_and_libraries_2018.3.222/linux/mpi/include64 \
         -L/home/ohnishi/usr/local/lib \
         -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_lapack95_lp64 -lmkl_core -lpthread \
         -liomp5 -lmpi -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

IFFTW = -I/home/ohnishi/usr/local/include

# --------------------------
# general rules
# --------------------------
.cpp.o:
	${CC} ${CFLAGS} ${LFLAGS} -c $<

PROG = anphon

cxxsrc = dump2sed.cpp myfftw3d.cpp file_handling.cpp 
objs=${cxxsrc:.cpp=.o}

default: dump2sed

all: ${PROG}

dump2sed: ${objs}
	${CC} ${CFLAGS} ${LFLAGS} -o $@ ${objs}

clean:
	rm -f ${objs}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend

