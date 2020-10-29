include make.inc
LOCALOBJ= \
	gen_tools.o \
	cons_tools.o \
	comm_tools.o \
	para_tools.o \
	mesh_tools.o \
	solp_tools.o \
	solu_toolsII.o \
	part_tools.o \
	run.o \
        fastfouriertransform.o \
        triseq.o \
        tripar.o \
        btrithom.o \
        poisson.o \
        poisson_solve.o \
        momentum_coef.o \
        momentum_terms.o \
        step.o \
        communicate.o \
        divergence.o \
        bessel.o \
        loop.o 

a.exe: main.o
	mkdir -p VR2D
	$(LNK) -o  VR_00.EXE $(OPT)  main.o $(LOCALOBJ) $(PETSC_LIB) $(LIB_FFTW) 
	mv VR_00.EXE VR2D

# Here are the compile steps

main.o:main.F90 $(LOCALOBJ) 

%.o: %.F90
	$(F90) -c  $(OPT) $<

%.o: %.f
	$(ff) -c $(OPT) $<


# This entry allows you to type " make cl " to get rid of
# all object and module files
cl:
	 rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo  V*.f *.dbg album F.err *~

nuovtk:
	\rm -rf EXE/VTKs/
