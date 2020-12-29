#Common
ROOT           = ~/

#Compiler options
PREPROC        = -cpp
#OPT            = -O3  -ffree-line-length-none -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread -I/usr/lib/x86_64-linux-gnu/openmpi/lib

#OPT           = -g  -fcheck='all'   -ffree-line-length-none
#OPT           = -g -O0 -traceback -CB
OPT           = -O3 

#FFTW related
INC_FFTW       = -I/home/andreas/Local/fftw/fftw-3.3.9/build/include/
LIB_FFTW       = /home/andreas/Local/fftw/fftw-3.3.9/build/lib/libfftw3.a
INC_MPI        = -I/home/andreas/Local/openmpi/openmpi-4.0.5/build/include
LIB_MPI        = /home/andreas/Local/openmpi/openmpi-4.0.5/build/include/lib
FC             = /home/andreas/Local/openmpi/openmpi-4.0.5/build/bin/mpif90
#Petsc related

#Definition of compiler and linker
F90	 =  /home/andreas/Local/openmpi/openmpi-4.0.5/build/bin/mpif90 $(PREPROC) $(OPT) $(INC_FFTW) $(INC_MPI) -L$(LIB_MPI)

ff  	 = /home/andreas/Local/openmpi/openmpi-4.0.5/build/bin/mpif90
LNK	 = /home/andreas/Local/openmpi/openmpi-4.0.5/build/bin/mpif90  $(FLINKER) $(OPT) $(INC_FFTW)  -L$(LIB_FFTW) -L$(LIB_MPI) -lmpi -I/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi

