ifort -c -r8 fla_postpro.F90 
ifort -c -r8 svd.for 
ifort -o fla.x svd.o fla_postpro.o 

