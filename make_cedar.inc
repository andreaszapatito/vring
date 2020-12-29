#Common
ROOT           = ~/

#Compiler options
PREPROC        = -cpp
#OPT            = -O3  -ffree-line-length-none -I/home/andreas/Local/openmpi/build/include -pthread -I/home/andreas/Local/openmpi/build/lib/

#OPT           = -g  -fcheck='all'   -ffree-line-length-none
OPT           = -O3 

#FFTW related
INC_FFTW       = -I/home/andreas/Local/fftw/fftw-3.3.8/install/include/
LIB_FFTW       = /home/andreas/Local/fftw/fftw-3.3.8/install/lib/libfftw3.a
INC_MPI        = -I/Local/openmpi/build/include
LIB_MPI        = /home/andreas/Local/openmpi/build/lib/
FC             = /home/andreas/Local/openmpi/build/bin/mpif90
#Petsc related

#Definition of compiler and linker
F90	 =    $(FC) $(PREPROC) $(OPT) $(INC_FFTW) $(INC_MPI) -L$(LIB_MPI)

ff  	 = $(FC)
LNK	 = $(FC) $(FLINKER) $(OPT) $(INC_FFTW)  -L$(LIB_FFTW) -L$(LIB_MPI) -lmpi $(INC_MPI) -L/usr/lib -L/home/andreas/Local/openmpi/build/lib/  -lmpi_mpifh -lmpi -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi

>>>>>>> cd992c7187a781c0ef3a4e67e3567bd0507a1d5e

