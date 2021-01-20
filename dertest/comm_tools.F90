module comm_type
  implicit none
#include <mpif.h>


  type :: comm
    integer                             :: comm_0         ! global comm
    integer                             :: np             ! number of tasks
    integer                             :: ip             ! global task id

    integer                             :: dd             ! number of partitioning dimensions
    integer,  dimension(:), pointer     :: comm_a         ! communicator per partitioning dimension
    integer,  dimension(:), allocatable :: ip_a           ! ranks per comm_A
    integer,  dimension(:), allocatable :: np_a           ! tasks per comm_A
  end type comm
end module comm_type


module comm_tools
  use comm_type
  use cons_tools
  use gen_tools
  implicit none

  contains

  subroutine create_comm(com,con)
    implicit none
    type(cons)         :: con
    type(comm)         :: com

    integer :: i, ierr
    logical :: flag
    character :: a
    
    open(unit=22,file=trim(adjustl(con%inp_file)), form='formatted', status='unknown')

    call read_until(22,'===Comm===',flag)
!    if (.not.flag) call err_n_stop('failure in input_comm',(/com%ip/))


!  Read number of splitting directions (3)
    read(22,*) a,com%dd
    allocate(com%np_a(com%dd))
    allocate(com%ip_a(com%dd))
!  Read number of tasks per direction
    read(22,*) a,com%np_a
    close(22)

    if (com%np.ne.com%np_a(1)*com%np_a(2)*com%np_a(3)) then
      write (*,*) 'error in partitioning'
      stop
    endif
!  Create partitioning per direction 
    call create_par_comm(com%comm_0,com%comm_a,com%np_a)
    do i = 1, com%dd
      call mpi_comm_rank(com%comm_a (i),com%ip_a(i),ierr)
    end do
    write (*,*) " RANKS ",com%ip,com%ip_a
  end subroutine create_comm

  subroutine create_par_comm(comm_in, comm_out, ndim,opt_reorder)
    implicit none
    integer,                intent(in)    :: comm_in
    integer,  dimension(:), pointer       :: comm_out
    integer,  dimension(:), intent(inout) :: ndim
    integer,  dimension(:), optional      :: opt_reorder

    integer   :: nprocs, rank, code, p, i
    integer   :: dd

    integer   :: comm_cart
    integer,  dimension(:), pointer    :: coord_cart
    logical,  dimension(size(ndim))    :: period, remains
    logical                            :: reorder
    integer,  dimension(size(ndim))    :: act_ndim



    integer*8 c,b,one
    call mpi_comm_size(comm_in, nprocs, code)
    call mpi_comm_rank(comm_in, rank  , code)


    dd = size(ndim)
    p = 1
    do i  = 1, size(ndim)
       p = p*ndim(i)
    end do

    if (nprocs/=p) then
       call warning('changing procs distribution',(/rank/))
       p = p/ndim(dd)
       if ( (p.le.nprocs) .and. (mod(nprocs,p)==0) ) then
          ndim(dd) = nprocs/p
       else
          ndim = 1
!!$       ndim(dd-1) = nprocs
          ndim(dd) = nprocs
       end if
    end if

    allocate(comm_out(dd), coord_cart(dd))
    period  = .false.
    reorder = .false.

    if (present(opt_reorder)) then
       act_ndim = ndim(opt_reorder)
    else
       act_ndim = ndim
    end if


    call mpi_cart_create(comm_in, dd, act_ndim, period, reorder, comm_cart, code)
    call mpi_comm_rank(comm_cart, rank, code)
    call mpi_cart_coords(comm_cart, rank, dd, coord_cart, code)

    do i = 1, dd
       remains=.false.
       if (present(opt_reorder)) then
          remains(opt_reorder(i)) = .true.
       else
          remains(i) = .true.
       end if
       call mpi_cart_sub(comm_cart, remains, comm_out(i), code)

       

    end do
!    call MPI_ALLGATHER(c, one, MPI_INTEGER, b, one, MPI_INTEGER, petsc_comm_world, code)
  end subroutine create_par_comm

end module comm_tools
