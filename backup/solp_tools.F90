module solp_type
!#include <petsc/finclude/petscdef.h>
!  use petsc
  implicit none
  save
! Poisson solver object

  type :: solp

! Data
    real(kind=8), allocatable :: rhs(:,:,:), sol(:,:,:), ref(:,:,:), dph(:,:,:), qcap(:,:,:)
!   real(kind=8), allocatable :: act_rhs(:), act_sol(:)

! Parameters
    integer                                  :: nsys
    real(kind=8)                             :: btri_mult
    character(len=5)                         :: meth_tripar

! Petsc objects
!    KSP,     dimension(:), allocatable       :: ksp
!    logical, dimension(:), allocatable       :: used
!    Vec,     dimension(:), allocatable       :: vx
!    Vec,     dimension(:), allocatable       :: vb

! Utility arrays
    real(kind=8), dimension(:), allocatable  :: bb_base, aa_base, cc_base

    real(kind=8)                             :: r1
    real(kind=8)                             :: r2
  end type solp

end module solp_type


module solp_tools
  use solp_type
  use mesh_type
  use cons_type
  use comm_type
  use para_type
  implicit none

  contains

  subroutine alct_solp(sp,msh)
    implicit none
    type(solp)         :: sp 
    type(mesh_a)       :: msh
    type(cons)         :: con
    type(comm)         :: com


    allocate(sp%rhs(msh%ntheta,msh%nr,msh%nz))
!    allocate(sp%qcap(msh%ntheta,msh%nr,msh%nz))
    allocate(sp%sol(msh%ntheta,msh%nr,msh%nz))
    allocate(sp%ref(msh%ntheta,msh%nr,msh%nz))
    allocate(sp%dph(msh%ntheta,0:msh%nr+1,0:msh%nz+1))
!   allocate(pct_rhs(msh%nr*msh%nz))
!   allocate(act_sol(msh%nr*msh%nz))


  end subroutine alct_solp

  subroutine init_solp(sp,msh)
    implicit none
    type (solp)  :: sp
    type (mesh_a):: msh

    real(kind=8)      :: rr, zz, th
    integer           :: i,j,k

    sp%sol = 0.d0
    DO k = 1, msh%nz
       DO j = 1, msh%nr
          DO i = 1, msh%ntheta
             rr = msh%rm(j)
             zz = msh%zm(k)
             th = msh%thm(i)
  
             sp%ref(i,j,k) = sol_ex(th,rr,zz)
             sp%rhs(i,j,k) = rhs_ex(th,rr,zz)

             sp%rhs(i,j,k) = 0.d0
             
             

          END DO
       END DO
    END DO



  end subroutine init_solp

  function sol_ex(th,rr,zz) result(ff)
    implicit none
    real(kind=8)      :: rr, zz, th
    real(kind=8)      :: ff,pi
 
    pi = 4.0 * atan (1.0)


    ff = cos(2.d0*pi*zz)/(4.d0*pi**2)*(rr**3-1.5d0*rr**2)*&
        (1.d0-0.5d0*cos(5.d0*th))


!    ff = cos(2.d0*pi*zz)/(4.d0*pi**2)*(rr**3-1.5d0*rr**2)


  end function sol_ex

  function rhs_ex(th,rr,zz) result(ff)
    implicit none
    real(kind=8)      :: rr, zz, th
    real(kind=8)      :: ff, pi

    pi = 4.0 * atan (1.0)

   ff = cos(2.d0*pi*zz)*(rr**3-1.5d0*rr**2)*(1.d0-0.5d0*cos(5.d0*th))  + &
       cos(2.d0*pi*zz)*(6.d0-9.d0*rr)/(4.d0*pi**2)*(1.d0-0.5d0*cos(5.d0*th)) + &
        cos(2.d0*pi*zz)/(4.d0*pi**2)*(rr**3-1.5d0*rr**2)*(-0.5d0*25.d0*cos(5.d0*th))/rr**2

!    ff = cos(2.d0*pi*zz)*(rr**3-1.5d0*rr**2) + &
!         cos(2.d0*pi*zz)*(6.d0-9.d0*rr)/(4.d0*pi**2)

  end function rhs_ex



end module solp_tools
