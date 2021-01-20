module para_type
!#include <petsc/finclude/petscdef.h>
  implicit none

  type :: para
    real(kind=8)              :: u1,u2,u3,p,Re,Pe
    real(kind=8)              :: uring,uswrl,ufluc,kfluc
    character(len=4)          :: poisson_type
    integer                   :: verbose
    integer                   :: totstp
    integer                   :: expstp 
    integer                   :: anlstp 
    integer                   :: restrt


    integer                   :: istep  
    integer                   :: nstep
    integer                   :: nstepi  
    integer                   :: nstepf  
    integer                   :: nrkstep
    integer                   :: irkstep
    real(kind=8)              :: crkalm(3)
    real(kind=8)              :: crkgam(3)
    real(kind=8)              :: crkrom(3)

    real(kind=8)              :: dt
    real(kind=8)              :: unom1
    real(kind=8)              :: unom0
    real(kind=8)              :: itime
    real(kind=8)              :: ntime
    real(kind=8)              :: ntimei
    real(kind=8)              :: ntimef
    real(kind=8)              :: r0
    real(kind=8)              :: rw
    real(kind=8)              :: t0
    real(kind=8),dimension(20):: stat

    real(kind=8)              :: slip
  end type para

end module para_type


module para_tools

  use para_type

  use comm_tools
  use cons_tools

  use gen_tools
  implicit none
  contains

  subroutine create_para(par,com,con)
    implicit none
   
    type(para)         :: par
    type(comm)         :: com
    type(cons)         :: con
    
    logical            :: flag
    character          :: a
    integer            :: i
    
    open(unit=22,file=trim(adjustl(con%inp_file)), form='formatted', status='unknown')

    call read_until(22,'===Para===',flag)
!    if (.not.flag) call err_n_stop('failure in create_para',(/com%ip/))

    read(22,*) a,par%dt
    read(22,*) a,par%totstp
    read(22,*) a,par%expstp
    read(22,*) a,par%anlstp
    read(22,*) a,par%restrt
    read(22,*) a,par%verbose
    read(22,*) a,par%u1
    read(22,*) a,par%u2
    read(22,*) a,par%u3
    read(22,*) a,par%uring,par%uswrl,par%ufluc,par%kfluc
    read(22,*) a,par%rw,par%r0,par%t0
    read(22,*) a,par%p
    read(22,*) a,par%Re
    read(22,*) a,par%Pe
    read(22,*) a,par%nrkstep
    read(22,*) a,par%slip
    close(22)

    par%Re=par%Re

    open(unit=22,file=trim(adjustl(con%inp_file)), form='formatted', status='unknown')

    call read_until(22,'===Flag===',flag)
!    if (.not.flag) call err_n_stop('failure in create_para',(/com%ip/))
    read(22,*) a,par%poisson_type
    close(22)

    par%nstepi=0   
    par%nstepf=par%nstepi+par%totstp
    par%ntimei=0.0
    par%ntimef=par%ntimei+par%totstp*par%dt

    if (par%nrkstep.eq.3) then
       par%crkgam(1)=8.d0/15.d0
       par%crkgam(2)=5.d0/12.d0
       par%crkgam(3)=3.d0/4.d0
       par%crkrom(1)=0.d0
       par%crkrom(2)=-17.d0/60.d0
       par%crkrom(3)=-5.d0/12.d0

       do i=1,3
         par%crkalm(i)=par%crkgam(i)+par%crkrom(i)
       enddo


    endif
    if (par%nrkstep.eq.1) then
       par%crkgam(1)=1.5d0
       par%crkgam(2)=0.d0
       par%crkgam(3)=0.d0
       par%crkrom(1)=-0.5d0
       par%crkrom(2)=0.d0
       par%crkrom(3)=0.d0

       do i=1,3
         par%crkalm(i)=par%crkgam(i)+par%crkrom(i)
       enddo
    endif

    
  end subroutine create_para

end module para_tools
