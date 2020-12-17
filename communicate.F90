module communicate
!#include <petsc/finclude/petscdef.h>
!  use petsc

  use run_tools
  use cons_tools
  use para_tools
  use comm_tools
  use mesh_tools
  use solu_tools

  implicit none

  contains

  subroutine communicateeq(var)

    implicit none
    type(run)         :: var

!    call commq(var,var%su%ru1,2)
!    call commq(var,var%su%ru2,2)
!    call commq(var,var%su%ru3,2)

    call commq(var,var%su%ru1,3)
    call commq(var,var%su%ru2,3)
    call commq(var,var%su%ru3,3)
    call commq(var,var%su%q1,2)
    call commq(var,var%su%q1,3)

    call commq(var,var%su%q2,2)
    call commq(var,var%su%q2,3)

    call commq(var,var%su%q3,2)
    call commq(var,var%su%q3,3)

    call commq(var,var%su%dq1,2)
    call commq(var,var%su%dq1,3)

    call commq(var,var%su%dq2,2)
    call commq(var,var%su%dq2,3)

    call commq(var,var%su%dq3,2)
    call commq(var,var%su%dq3,3)

    call commq(var,var%su%qcap,2)
    call commq(var,var%su%qcap,3)

    call commq(var,var%sp%dph,2)
    call commq(var,var%sp%dph,3)


  end subroutine communicateeq


  subroutine commq(var,fld,idir)
    implicit none
    type(run)                           :: var
    real(kind=8)              :: fld (:,0:,0:)

    integer                   :: idir
    integer                   :: iibuff,ifbuff,jibuff,jfbuff,nbuff,mbuff
    integer                   :: i,j
    integer                   :: tag,msgsize
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    real(kind=8), allocatable :: send_bf(:,:),recv_bf(:,:)

    integer                   :: ierr

    nbuff=var%msh%nthetac
    if (idir.eq.2) then
      mbuff=var%msh%nzc
    elseif (idir.eq.3) then
      mbuff=var%msh%nrc
    endif

    allocate (send_bf(nbuff,mbuff))
    allocate (recv_bf(nbuff,mbuff))

    call mpi_barrier(var%com%comm_0, ierr)
!   communicate east nort

    msgsize=nbuff*mbuff
    if (var%com%ip_a(idir).lt.var%com%np_a(idir)-1) then
      do i=1,nbuff
        do j=1,mbuff
          if (idir.eq.2) then
            send_bf(i,j)=fld(i,var%msh%nrc,j)
          elseif (idir.eq.3) then
            send_bf(i,j)=fld(i,j,var%msh%nzc)
          endif
        enddo
      enddo
      tag = 1
      call mpi_send(send_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip_a(idir)+1,tag,var%com%comm_a(idir),ierr)
    endif

    if (var%com%ip_a(idir).gt.              0) then
      tag = 1
      call mpi_recv(recv_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip_a(idir)-1,tag,var%com%comm_a(idir),status,ierr)

      do i=1,nbuff
        do j=1,mbuff
          if (idir.eq.2) then
             fld(i,0,j)=recv_bf(i,j)
          elseif (idir.eq.3) then
             fld(i,j,0)=recv_bf(i,j)
          endif
        enddo
      enddo
    endif

!   communicate west sout
    call mpi_barrier(var%com%comm_0, ierr)

    if (var%com%ip_a(idir).gt.              0) then
      do i=1,nbuff
        do j=1,mbuff
          if (idir.eq.2) then
            send_bf(i,j)=fld(i,1,j)
          elseif (idir.eq.3) then
            send_bf(i,j)=fld(i,j,1)
          endif
        enddo
      enddo
      tag = 2
      call mpi_send(send_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip_a(idir)-1,tag,var%com%comm_a(idir),ierr)
    endif

    if (var%com%ip_a(idir).lt.var%com%np_a(idir)-1) then
      tag = 2
      call mpi_recv(recv_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip_a(idir)+1,tag,var%com%comm_a(idir),status,ierr)

      do i=1,nbuff
        do j=1,mbuff
          if (idir.eq.2) then
            fld(i,var%msh%nrc+1,j)=recv_bf(i,j)
          elseif (idir.eq.3) then
            fld(i,j,var%msh%nzc+1)=recv_bf(i,j)
          endif
        enddo
      enddo
    endif

    call mpi_barrier(var%com%comm_0, ierr)
    msgsize=nbuff

    if ((var%com%ip_a(2).gt.0).and.(var%com%ip_a(3).lt.var%com%np_a(3)-1)) then
      do i=1,nbuff
        send_bf(i,1)=fld(i,1,var%msh%nzc)
      enddo
      tag = 2
      call mpi_send(send_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip-var%com%np_a(3)+1,tag,var%com%comm_0,ierr)
    endif


    if ((var%com%ip_a(3).gt.0).and.(var%com%ip_a(2).lt.var%com%np_a(2)-1)) then
      tag = 2
      call mpi_recv(recv_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip+var%com%np_a(3)-1,tag,var%com%comm_0,status,ierr)
      do i=1,nbuff
        fld(i,var%msh%nrc+1,0)=recv_bf(i,1)
      enddo
    endif


    msgsize=nbuff

    if ((var%com%ip_a(3).gt.0).and.(var%com%ip_a(2).lt.var%com%np_a(2)-1)) then
      do i=1,nbuff
        send_bf(i,1)=fld(i,var%msh%nrc,1)
      enddo
      tag = 2
      call mpi_send(send_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip+var%com%np_a(3)-1,tag,var%com%comm_0,ierr)
    endif


    if ((var%com%ip_a(2).gt.0).and.(var%com%ip_a(3).lt.var%com%np_a(3)-1)) then
      tag = 2
      call mpi_recv(recv_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip-var%com%np_a(3)+1,tag,var%com%comm_0,status,ierr)
      do i=1,nbuff
        fld(i,0,var%msh%nzc+1)=recv_bf(i,1)
      enddo
    endif

    msgsize=nbuff

    if ((var%com%ip_a(3).gt.0).and.(var%com%ip_a(2).gt.0)) then
      do i=1,nbuff
        send_bf(i,1)=fld(i,1,1)
      enddo
      tag = 2
      call mpi_send(send_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip-var%com%np_a(3)-1,tag,var%com%comm_0,ierr)
    endif


    if ((var%com%ip_a(2).lt.var%com%np_a(2)-1).and.(var%com%ip_a(3).lt.var%com%np_a(3)-1)) then
      tag = 2
      call mpi_recv(recv_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip+var%com%np_a(3)+1,tag,var%com%comm_0,status,ierr)
      do i=1,nbuff
        fld(i,var%msh%nrc+1,var%msh%nzc+1)=recv_bf(i,1)
      enddo
    endif

    msgsize=nbuff

    if ((var%com%ip_a(2).lt.var%com%np_a(2)-1).and.(var%com%ip_a(3).lt.var%com%np_a(3)-1)) then
      do i=1,nbuff
        send_bf(i,1)=fld(i,var%msh%nrc,var%msh%nzc)
      enddo
      tag = 2
      call mpi_send(send_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip+var%com%np_a(3)+1,tag,var%com%comm_0,ierr)
    endif


    if ((var%com%ip_a(2).gt.0).and.(var%com%ip_a(3).gt.0)) then
      tag = 2
      call mpi_recv(recv_bf,msgsize,MPI_DOUBLE_PRECISION,var%com%ip-var%com%np_a(3)-1,tag,var%com%comm_0,status,ierr)
      do i=1,nbuff
        fld(i,0,0)=recv_bf(i,1)
      enddo
    endif

   


    deallocate (send_bf)
    deallocate (recv_bf)
    call mpi_barrier(var%com%comm_0, ierr)



  end subroutine commq
end module communicate
