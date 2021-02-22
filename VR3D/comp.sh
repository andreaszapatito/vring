ifort -c -O3 -r8 fla_postpro_III.F90 
ifort -c -O3 -r8 powel.F90 
ifort -c -O3 -r8 svd.for 
ifort -O3 -o fla.x svd.o fla_postpro_III.o  powel.o

