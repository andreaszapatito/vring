ifort -c -g -r8 fla_postpro.F90 
ifort -c -g -r8 powel.F90 
ifort -c -g -r8 svd.for 
ifort -o fla.x svd.o fla_postpro.o  powel.o

