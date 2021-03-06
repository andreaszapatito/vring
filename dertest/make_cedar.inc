#Common
ROOT           = ~/

#Compiler options
PREPROC        = -cpp
#OPT            = -O3  -ffree-line-length-none -I/home/andreas/Local/openmpi/openmpi-4.0.5/build/include -pthread -I/home/andreas/Local/oopenmpi/openmpi-4.0.5/build/lib/

#OPT           = -g -CB -fcheck='all'   -ffree-line-length-none
OPT           = -O3 

#FFTW related
INC_FFTW       = -I/home/andreas/Local/fftw/fftw-3.3.9/build/include/
LIB_FFTW       = /home/andreas/Local/fftw/fftw-3.3.9/build/lib/libfftw3.a
INC_MPI        = -I/Local/openmpi/openmpi-4.0.5/build/include
LIB_MPI        = /home/andreas/Local/openmpi/openmpi-4.0.5/build/lib/
FC             = /home/andreas/Local/openmpi/openmpi-4.0.5/build/bin/mpif90
#Petsc related

#Definition of compiler and linker
F90	 =    $(FC) $(PREPROC) $(OPT) $(INC_FFTW) $(INC_MPI) -L$(LIB_MPI)

ff  	 = $(FC)
LNK	 = $(FC) $(FLINKER) $(OPT) $(INC_FFTW)  -L$(LIB_FFTW) -L$(LIB_MPI) -lmpi $(INC_MPI) -L/usr/lib -L/home/andreas/Local/openmpi/openmpi-4.0.5/build/lib/  -lmpi_mpifh -lmpi -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi

