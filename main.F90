program main 
!#include <petsc/finclude/petscdef.h>
!  use petsc

  use run_tools

  use cons_tools
  use para_tools
  use comm_tools
  use mesh_tools
  use solu_tools
  use solp_tools

  use poisson
  use poisson_solve
  use loop
  use part_tools
  implicit none

  integer               :: code,c,b

  type(run)             :: var

  integer               :: ierr


! initialise 
  call mpi_init(ierr)
 
  call mpi_comm_size   (MPI_COMM_WORLD, var%com%np, code)
  call mpi_comm_rank   (MPI_COMM_WORLD, var%com%ip, code)
  var%com%comm_0=MPI_COMM_WORLD
  var%par%verbose=1



  call verbose(var,"Initialization called")
! create constants / input
  call create_cons    (var%con)
  call verbose(var,"Create constants")
! create communicate / input 
  call create_comm    (var%com,var%con)
  call verbose(var,"Create communication")
! create parameters / input
  call create_para    (var%par,var%com,var%con)
  call verbose(var,"Create parameters")

! create mesh /input
  call create_mesh_a  (var%msh,var%par,var%com,var%con)
  call verbose(var,"Create mesh")

  call verbose(var,"Objects create called")

! output case setup
  call echo_var(var)


! allocate pressure and velocity fields
  call alct_solp      (var%sp,var%msh)                  ;  call verbose(var,"alct_solp called")
  call alct_solu      (var%su,var%msh,var%com)          ;  call verbose(var,"alct_solu called")

  call init_solp      (var%sp,var%msh)                  ;  call verbose(var,"init_solp called")
  call init_solu      (var%su,var%msh,var%par,var%com)  ;  call verbose(var,"init_solu called")

  if (var%par%restrt.eq.1) call restart_load (var);     ;  call verbose(var,"restart   called")

  call prep_poisson_solver(var)                         ;  call verbose(var,"prep_poisson_solver called")
  var%su%dq1=0.0
  var%su%dq2=0.0
  var%su%dq3=0.0

  call init_part(var%prt,var%msh)
!!  call ffemexport(var,0)
  call iterate(var)                                     ;  call verbose(var,"iterate      ")
 
  var%con%timming0= MPI_Wtime()
  var%con%timming1= MPI_Wtime()
  var%con%timming2= MPI_Wtime()
 
  var%con%timmingPO=var%con%timming2-var%con%timming1

  call MPI_Finalize(ierr)

  contains
 
end program main 

