module solu_type
!#include <petsc/finclude/petscdef.h>
!  use petsc
  use gen_tools
  use comm_type
  implicit none
  save
! Poisson solver object

  type :: solu

! Data
    real(kind=8), allocatable :: u1(:,:,:), u2(:,:,:), u3(:,:,:) 
    real(kind=8), allocatable :: q1(:,:,:), q2(:,:,:), q3(:,:,:) 
    real(kind=8), allocatable :: dq1(:,:,:),dq2(:,:,:),dq3(:,:,:)

    real(kind=8), allocatable :: q1s(:,:),  q2s(:,:),  q3s(:,:) ,dq3no(:,:)
    real(kind=8), allocatable :: qcap(:,:,:),pres(:,:,:)

    real(kind=8), allocatable :: qrhs(:,:,:),ru1(:,:,:),ru2(:,:,:),ru3(:,:,:),prg(:,:,:,:)
    real(kind=8), allocatable :: qg(:,:)
    real(kind=8), allocatable :: ql(:,:)




! Parameters
    integer                                  :: nsys

! Utility arrays
    real(kind=8), dimension(:), allocatable  :: bb
    real(kind=8)                             :: r1
    real(kind=8)                             :: r2

    real(kind=8), dimension(:), allocatable  :: atri
    real(kind=8), dimension(:), allocatable  :: btri
    real(kind=8), dimension(:), allocatable  :: ctri
    real(kind=8), dimension(:), allocatable  :: rtri
    real(kind=8), dimension(:), allocatable  :: stri


    real(kind=8), dimension(:), allocatable  :: h1d11, h1d22d, h1d22n, h1dare
    real(kind=8), dimension(:), allocatable  :: h2d11, h2d22m, h2d22p, h2dare

    real(kind=8), dimension(:), allocatable  :: s1t11, s1t2n, s1t2s, s1t2c
    real(kind=8), dimension(:), allocatable  :: s2t11, s2t2p, s2t2m, s2t2c

    real(kind=8), dimension(:), allocatable  :: s3t2n, s3t2s
    real(kind=8), dimension(:), allocatable  :: amj1, acj1, apj1
    real(kind=8), dimension(:), allocatable  :: amj2, acj2, apj2
    real(kind=8), dimension(:), allocatable  :: amj3, acj3, apj3
    real(kind=8), dimension(:), allocatable  :: amjs, acjs, apjs
    real(kind=8), dimension(:), allocatable  :: amk1, ack1, apk1
    real(kind=8), dimension(:), allocatable  :: amk2, ack2, apk2
    real(kind=8), dimension(:), allocatable  :: amk3, ack3, apk3
    real(kind=8), dimension(:), allocatable  :: amks, acks, apks

    real(kind=8), dimension(:,:), allocatable:: corr1n,corr2n,corr3n,csort,tconv
    real(kind=8), dimension(:,:), allocatable:: corr1s,corr2s,corr3s

    real(kind=8), dimension(:,:), allocatable  :: ami, aci, api, fi, gam2


 


  end type solu

end module solu_type


module solu_tools
  use mesh_type
  use cons_type
  use comm_type
  use para_type
  use solu_type
  implicit none

  contains

  subroutine alct_solu(su,msh,com)
    implicit none
    type(solu)         :: su
    type(comm)         :: com
    type(mesh_a)       :: msh

    allocate(su%qcap(msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%pres(msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))

!   Allocate q1 q2 q3 and allow for communication buffer
    if (com%ip.eq.0)  allocate (su%qg  (msh%nrg,msh%nzg))
    allocate(su%ql  (msh%nr,msh%nz))

    allocate(su%u1  (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%q1  (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%dq1 (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%ru1 (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))

    allocate(su%u2  (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%q2  (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%dq2 (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%ru2 (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))

    allocate(su%u3  (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%q3  (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%dq3 (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))
    allocate(su%dq3no (msh%nthetac,0:msh%nrc+1))
    allocate(su%ru3 (msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))

    allocate(su%qrhs(msh%nthetac,0:msh%nrc+1,0:msh%nzc+1))

    allocate(su%prg(msh%ntheta,msh%nr,msh%nz,3))


    allocate(su%q1s(msh%nthetac,0:msh%nrc))
    allocate(su%q2s(msh%nthetac,0:msh%nrc))
    allocate(su%q3s(msh%nthetac,0:msh%nrc))



    allocate(su%h1d11    (msh%iirc1:msh%ifrc1))
    allocate(su%h1d22d   (msh%iirc1:msh%ifrc1))
    allocate(su%h1d22n   (msh%iirc1:msh%ifrc1))
    allocate(su%h1dare   (msh%iirc1:msh%ifrc1))

    allocate(su%h2d11    (msh%iirc1:msh%ifrc1))
    allocate(su%h2d22m   (msh%iirc1:msh%ifrc1))
    allocate(su%h2d22p   (msh%iirc1:msh%ifrc1))
    allocate(su%h2dare   (msh%iirc1:msh%ifrc1))

    allocate(su%s1t11    (msh%iirc1:msh%ifrc1))
    allocate(su%s1t2n    (msh%iirc1:msh%ifrc1))
    allocate(su%s1t2s    (msh%iirc1:msh%ifrc1))
    allocate(su%s1t2c    (msh%iirc1:msh%ifrc1))

    allocate(su%s2t11    (msh%iirc2:msh%ifrc2))
    allocate(su%s2t2p    (msh%iirc2:msh%ifrc2))
    allocate(su%s2t2m    (msh%iirc2:msh%ifrc2))
    allocate(su%s2t2c    (msh%iirc2:msh%ifrc2))

    allocate(su%s3t2n    (msh%iirc1:msh%ifrc1))
    allocate(su%s3t2s    (msh%iirc1:msh%ifrc1))

    allocate(su%amj1     (0:msh%nr+1))
    allocate(su%acj1     (0:msh%nr+1))
    allocate(su%apj1     (0:msh%nr+1))

    allocate(su%amj2     (0:msh%nr+1))
    allocate(su%acj2     (0:msh%nr+1))
    allocate(su%apj2     (0:msh%nr+1))

    allocate(su%amj3     (0:msh%nr+1))
    allocate(su%acj3     (0:msh%nr+1))
    allocate(su%apj3     (0:msh%nr+1))

    allocate(su%amjs     (0:msh%nr+1))
    allocate(su%acjs     (0:msh%nr+1))
    allocate(su%apjs     (0:msh%nr+1))

    allocate(su%amk1     (0:msh%nzc))
    allocate(su%ack1     (0:msh%nzc))
    allocate(su%apk1     (0:msh%nzc))

    allocate(su%amk2     (0:msh%nzc))
    allocate(su%ack2     (0:msh%nzc))
    allocate(su%apk2     (0:msh%nzc))

    allocate(su%amk3     (0:msh%nzc))
    allocate(su%ack3     (0:msh%nzc))
    allocate(su%apk3     (0:msh%nzc))

    allocate(su%amks     (0:msh%nzc))
    allocate(su%acks     (0:msh%nzc))
    allocate(su%apks     (0:msh%nzc))

    allocate(su%atri     (1:max(msh%nzc,msh%nrc)))
    allocate(su%btri     (1:max(msh%nzc,msh%nrc)))
    allocate(su%ctri     (1:max(msh%nzc,msh%nrc)))
    allocate(su%rtri     (1:max(msh%nzc,msh%nrc)))
    allocate(su%stri     (1:max(msh%nzc,msh%nrc)))

    allocate(su%ami      (0:msh%nrc,0:msh%nthetac))
    allocate(su%aci      (0:msh%nrc,0:msh%nthetac))
    allocate(su%api      (0:msh%nrc,0:msh%nthetac))
    allocate(su%fi       (0:msh%nrc,0:msh%nthetac))
    allocate(su%gam2     (0:msh%nrc,0:msh%nthetac))



    allocate(su%corr1s   (msh%nthetac,msh%nrc))
    allocate(su%corr2s   (msh%nthetac,msh%nrc))
    allocate(su%corr3s   (msh%nthetac,msh%nrc))

    allocate(su%corr1n   (msh%nthetac,msh%nrc))
    allocate(su%corr2n   (msh%nthetac,msh%nrc))
    allocate(su%corr3n   (msh%nthetac,msh%nrc))
    allocate(su%csort    (msh%nthetac,msh%nrc))
    allocate(su%tconv    (msh%nthetac,msh%nrc))




  end subroutine alct_solu

  subroutine init_solu(su,msh,par,com)
    implicit none
    type (solu)  :: su
    type (mesh_a):: msh
    type (para)  :: par
    type (comm)  :: com

    real(kind=8)      :: rr, zz, th, po2,phi,rc,zc,radius,ur,uv,vv,wv,t1
    integer           :: i,j,k,ierr
!    par%uring=1.0
!    par%uswrl=0.2
!    par%ufluc=0.05
!    par%kfluc=8.0

!    if (com%ip.eq.0) then
!    
!    do k = 1,100
!          zc=k*0.01
!          rc=1.0
!          t1=0.1
!          write (*,*) zc,uvring(rc/par%r0,zc/par%r0,t1,par%r0/par%rw)
!    enddo
!    endif
!    call mpi_barrier(com%comm_0, ierr)
!    stop -1
    do k = msh%iizc1, msh%ifzc1
       do j = msh%iirc1, msh%ifrc1
          rc=msh%rm(j)
          zc=msh%zm(k)-msh%zmaxg/8.0
          if (par%uswrl.gt.0.1) zc=msh%zm(k)-msh%zmaxg/2.0

          !zc=msh%zm(k)-4.0
          wv=wvring(rc/par%r0,zc/par%r0,par%t0,par%r0/par%rw)!*par%Re
          do i = msh%iithetac1, msh%ifthetac1
             su%u1(i,j,k) = par%u1
             su%q1(i,j,k) = par%u1

             su%u1(i,j,k) = (wv*(par%uswrl+par%ufluc*cos(par%kfluc*msh%thc(i))))
             su%q1(i,j,k) = (wv*(par%uswrl+par%ufluc*cos(par%kfluc*msh%thc(i))))
!!             su%u1(i,j,k) = 0.0
!!             su%q1(i,j,k) = 0.0

         end do
!!       write (1000+com%ip,*) "wvring",msh%rm(j),msh%zm(k),wvring(rc,zc)
       end do
!!      write (1000+com%ip,*) 
    end do

    do k = msh%iizc2, msh%ifzc2
       do j = msh%iirc2, msh%ifrc2
          rc=msh%rc(j)
          zc=msh%zm(k)-msh%zmaxg/8.0
          if (par%uswrl.gt.0.1) zc=msh%zm(k)-msh%zmaxg/2.0
          !zc=msh%zm(k)-4.0
          !vv=1.0*(0.742/0.549)*vvring(rc/par%r0,zc/par%r0,par%t0,par%r0/par%rw)!par%Re
          vv=1.0*(0.8/0.5)*vvring(rc/par%r0,zc/par%r0,par%t0,par%r0/par%rw)!par%Re
          do i = msh%iithetac2, msh%ifthetac2

             su%u2(i,j,k) =-vv*(par%uring+par%ufluc*sin(par%kfluc*msh%thc(i)))
             su%q2(i,j,k) =-vv*(par%uring+par%ufluc*sin(par%kfluc*msh%thc(i)))*msh%rc(j)
             


          end do
!!       write (2000+com%ip,*) "vvring",msh%rm(j),msh%zm(k),vvring(rc,zc)
       end do
!!     write (2000+com%ip,*) 
    end do

!!    do k = msh%iizc3, msh%ifzc3
!!       do j = msh%iirc3, msh%ifrc3
    do k = 1, msh%ifzc3
       do j = 1, msh%ifrc3
          rc=msh%rm(j)
          zc=msh%zc(k)-msh%zmaxg/8.0
          if (par%uswrl.gt.0.1) zc=msh%zc(k)-msh%zmaxg/2.0
          !zc=msh%zm(k)-4.0
!          uv=1.0*(0.742/0.549)*uvring(rc/par%r0,zc/par%r0,par%t0,par%r0/par%rw)!*par%Re
          uv=1.0*(0.8/0.5)*uvring(rc/par%r0,zc/par%r0,par%t0,par%r0/par%rw)!*par%Re
          do i = msh%iithetac3, msh%ifthetac3
             su%u3(i,j,k) = par%u3
             su%q3(i,j,k) = par%u3





             su%u3(i,j,k) =-par%uring*uv
             su%q3(i,j,k) =-par%uring*uv


          end do
!!          if (k.gt.0.and.k.lt.msh%nzc.and.j.gt.0.and.j.lt.msh%ifrc3) write (3000+com%ip,*) "uvring",msh%rm(j),msh%zm(k),uvring(rc,zc)
       end do
!!      write (3000+com%ip,*) 
    end do



  end subroutine init_solu

end module solu_tools
