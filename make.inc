#Common
ROOT           = ~/

#Compiler options
PREPROC        = -cpp
#OPT            = -O3  -ffree-line-length-none -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread -I/usr/lib/x86_64-linux-gnu/openmpi/lib

#OPT           = -g  -fcheck='all'   -ffree-line-length-none
OPT           = -O3  -ffree-line-length-none -pthread

#FFTW related
INC_FFTW       = -I/home/andreas/Local/fftw/fftw-3.3.8/install/include/
LIB_FFTW       = /home/andreas/Local/fftw/fftw-3.3.8/install/lib/libfftw3.a
INC_MPI        = -I/Local/openmpiintel/install/include
LIB_MPI        = /home/andreas/Local/openmpiintel/install/lib/
FC             = mpif90
#Petsc related

#Definition of compiler and linker
F90	 =  mpif90 $(PREPROC) $(OPT) $(INC_FFTW) $(INC_MPI) -L$(LIB_MPI)

ff  	 = mpif90
LNK	 = mpif90 $(FLINKER) $(OPT) $(INC_FFTW)  -L$(LIB_FFTW) -L$(LIB_MPI) -lmpi -pthread -I/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi

