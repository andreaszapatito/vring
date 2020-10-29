module mesh_type
  implicit none
!#include <petsc/finclude/petscdef.h>
!#include <mpi.h>

  type :: mesh_a

     integer                   :: nr, ntheta, nz
     integer                   :: nrc, nthetac, nzc
     integer                   :: nrg, nthetag, nzg
     integer                   :: iizc1,ifzc1,iirc1,ifrc1,iithetac1,ifthetac1
     integer                   :: iizc2,ifzc2,iirc2,ifrc2,iithetac2,ifthetac2
     integer                   :: iizc3,ifzc3,iirc3,ifrc3,iithetac3,ifthetac3
     real(kind=8)              :: rmin, rmax, zmin, zmax, thmin, thmax
     real(kind=8)              :: rmaxg, zmaxg
     real(kind=8)              :: axsym
     integer                   :: if_grad_null
     real(kind=8)              :: dr, dtheta, dz, dx1, dx2, dx3
     real(kind=8)              :: dx1q, dx2q, dx3q
     real(kind=8), allocatable :: rc(:), thc(:), zc(:), da(:),qc(:)
     real(kind=8), allocatable :: rm(:), thm(:), zm(:)
     real(kind=8), allocatable :: dmr(:)
     real(kind=8), allocatable :: dcr(:)

     integer,      allocatable :: ntrz(:,:)
     integer,      allocatable :: rz_ig_start(:)
     integer                   :: rz_rk_neigh(4), rz_ig_start_sc


     integer,      allocatable :: imv(:),ipv(:),kmv(:),kpv(:),jmv(:),jpv(:),jpc(:),jmc(:)
     integer,      allocatable :: isym(:)
  end type mesh_a          

  type :: mesh_b

     integer                   :: nx,  ny,  nz
     integer                   :: nxg, nyg, nzg
     real(kind=8)              :: xmin, xmax
     real(kind=8)              :: ymin, ymax
     real(kind=8)              :: zmin, zmax
     real(kind=8)              :: dx, dy, dz
     real(kind=8), allocatable :: xc(:), yc(:), zc(:)


  end type mesh_b

end module mesh_type

module mesh_tools
!  use petsc
  use mesh_type

  use para_tools
  use comm_tools
  use cons_tools

  use gen_tools

  
  implicit none
!INCLUDE "mpif.h"

contains

  subroutine create_mesh_a(msh,par,com,con)
    implicit none
    type(mesh_a)               :: msh
    type(cons)                 :: con
    type(comm)                 :: com
    type(para)                 :: par
    logical                    :: flag

    integer   :: i, j, k
    integer     :: c,b
    integer     :: one
    integer   :: ic, kc, jc
    integer   :: ig, jg, kg
    integer   :: code
    integer   :: nd
    integer, allocatable :: tmp_ntrz(:,:)
    character :: a

    open(unit=22, file=trim(adjustl(con%inp_file)), form='formatted', &
         status='unknown')
    call read_until(22, '===Mesh===', flag)

!   Read number of global nodes in theta, r and z

    if (flag) then
       read(22,*) a,msh%nthetag,msh%nrg,msh%nzg
       read(22,*) a,msh%rmax
       read(22,*) a,msh%zmax




       msh%axsym=1.0
       msh%if_grad_null=0   !original
       if (msh%nthetag.eq.2) msh%axsym=0.0

!   Number of global nodes should divide by the number of processors, 
! if not number of global nodes is changed

       if (mod(msh%nthetag,com%np_a(1))/=0) then
          call warning('ntheta_glob has been changed to fit nprocs, see param.out',com%ip_a)
          msh%nthetag = msh%nthetag-mod(msh%nthetag,com%np_a(1))+com%np_a(1)

       end if
       if (mod(msh%nrg,com%np_a(2))/=0) then
          call warning('nr_glob has been changed to fit nprocs, see param.out',com%ip_a)
          msh%nrg = msh%nrg-mod(msh%nrg,com%np_a(2))+com%np_a(2)

       end if
       if (mod(msh%nzg,com%np_a(3))/=0) then
          call warning('nz_glob has been changed to fit nprocs, see param.out',com%ip_a)
          msh%nzg = msh%nzg-mod(msh%nzg,com%np_a(3))+com%np_a(3)

       end if


!   Calculate number of local nodes ie ntheta, nr and nz

       msh%ntheta = msh%nthetag/com%np_a(1)
       msh%nr     = msh%nrg/com%np_a(2)
       msh%nz     = msh%nzg/com%np_a(3)


       msh%nthetac=msh%ntheta
       if (com%ip_a(2).eq.com%np_a(2)-1) then
         msh%nrc    =msh%nr+1
       else
         msh%nrc    =msh%nr
       endif

       if (com%ip_a(3).eq.com%np_a(3)-1) then
         msh%nzc    =msh%nz+1
       else
         msh%nzc    =msh%nz
       endif

!      Q1 dimensions

       msh%iithetac1=1
       msh%ifthetac1=msh%ntheta

       msh%iirc1=1
       msh%ifrc1=msh%nr
!       if (com%ip_a(2).eq.0) then
!         msh%iirc1=0
!         msh%ifrc1=msh%nr+1
!       else
!         msh%iirc1=1
!         msh%ifrc1=msh%nr
!       endif

       msh%iizc1=1 
       msh%ifzc1=msh%nz
!       if (com%ip_a(3).eq.0) then
!         msh%iizc1=1
!         msh%ifzc1=msh%nz+1
!       else
!         msh%iizc1=1
!         msh%ifzc1=msh%nz
!       endif

!      Q2 dimensions

       msh%iithetac2=1
       msh%ifthetac2=msh%ntheta

       if (com%ip_a(2).eq.com%np_a(2)-1) then
         msh%iirc2=1
         msh%ifrc2=msh%nr
       else
         msh%iirc2=1
         msh%ifrc2=msh%nr+1
       endif

       msh%iizc2=1 
       msh%ifzc2=msh%nz
!       if (com%ip_a(3).eq.0) then
!         msh%iizc2=1
!         msh%ifzc2=msh%nz+1
!       else
!         msh%iizc2=1
!         msh%ifzc2=msh%nz+1
!       endif

!      Q3 dimensions

       msh%iithetac3=1
       msh%ifthetac3=msh%ntheta

       msh%iirc3=1
       msh%ifrc3=msh%nr
!       if (com%ip_a(2).eq.0) then
!         msh%iirc3=0
!         msh%ifrc3=msh%nr+1
!       else
!         msh%iirc3=1
!         msh%ifrc3=msh%nr
!       endif

       if (com%ip_a(3).eq.0) then
         msh%iizc3=0
         msh%ifzc3=msh%nz
       else
         msh%iizc3=1
         msh%ifzc3=msh%nz
       endif


       msh%rmaxg = msh%rmax
       msh%zmaxg = msh%zmax

       msh%dtheta= 2.d0*acos(-1.d0)/msh%nthetag
       msh%dr    = msh%rmax/msh%nrg
       msh%dz    = msh%zmax/msh%nzg
       msh%thmin = com%ip_a(1)*2.d0*acos(-1.d0)/com%np_a(1)
       msh%thmax = (com%ip_a(1)+1)*2.d0*acos(-1.d0)/com%np_a(1)


!   Not sure for the following 

       msh%rmin  = com%ip_a(2)*msh%rmax/com%np_a(2)
       msh%rmax  = (com%ip_a(2)+1)*msh%rmax/com%np_a(2)
       msh%zmin  = com%ip_a(3)*msh%zmax/com%np_a(3)
       msh%zmax  = (com%ip_a(3)+1)*msh%zmax/com%np_a(3)





!   Indices des mailles
!   Direction theta (periodique: maille n1=maille 1
!                              : maille  0=maille n1m)
       allocate(msh%thm(1:msh%ntheta))
       allocate(msh%rm (0:msh%nr+1))      ! radius of the collocation point, account for neighboring cells
       allocate(msh%zm (0:msh%nz+1))

       allocate(msh%thc(msh%iithetac1:msh%ifthetac1))
       allocate(msh%rc (msh%iirc2-1:msh%ifrc2+1)) ! radius of face point, account for ghost cell
       allocate(msh%zc (msh%iizc3-1:msh%ifzc3+1))

       allocate(msh%qc(msh%iirc1-1:msh%ifrc1+1))
       allocate(msh%da(msh%iirc2-1:msh%ifrc2+1))


       allocate(msh%dmr(msh%iirc1-1:msh%ifrc1+1))
       allocate(msh%dcr(msh%iirc2-1:msh%ifrc2+1))

       allocate(msh%imv(msh%ntheta),msh%ipv(msh%ntheta),msh%kmv(msh%nz),msh%kpv(msh%nz),msh%jmv(msh%nr),msh%jpv(msh%nr),msh%jpc(msh%nr),msh%jmc(msh%nr))
       allocate(msh%isym(msh%ntheta))


      do  ic=2,msh%ntheta
        msh%imv(ic)=ic-1
      enddo
      do  ic=1,msh%ntheta-1
        msh%ipv(ic)=ic+1
      enddo

      msh%imv(1)=msh%ntheta
      msh%ipv(msh%ntheta)=1

!   Direction     z (pour gradients nuls a la frontiere)

      do  kc=1,msh%nz
        msh%kmv(kc)=kc-1
      enddo
      do  kc=1,msh%nz
        msh%kpv(kc)=kc+1
      enddo
      
      if (com%ip_a(3).eq.0)             msh%kmv(1)=1
      if (com%ip_a(3).eq.com%np_a(3)-1) msh%kmv(msh%nz)=msh%nz

!   Direction     r (condition de glissement=gradients nuls )

      do  jc=1,msh%nr
        msh%jmv(jc)=jc-1
      enddo
      do  jc=1,msh%nr
       msh%jpv(jc)=jc+1
      enddo

      if (com%ip_a(2).eq.0)             msh%jmv(1)=1
      if (com%ip_a(2).eq.com%np_a(2)-1) msh%jpv(msh%nr)=msh%nr


!   Direction     r (condition de glissement )

      do  jc=1,msh%nr
       msh%jpc(jc)=msh%jpv(jc)-jc
       msh%jmc(jc)=jc-msh%jmv(jc)
      enddo

!   Symetrie par rapport a l'axe r=0

      nd=msh%ntheta/2
      do ic=1,nd
       msh%isym(ic)=ic+nd
      enddo
      do ic=nd+1,msh%ntheta
       msh%isym(ic)=ic-nd
      enddo

!   Allocate central coordinates


       msh%dx1 = 1.d0/msh%dtheta
       msh%dx2 = 1.d0/msh%dr
       msh%dx3 = 1.d0/msh%dz
       msh%dx1q = msh%dx1**2
       msh%dx2q = msh%dx2**2
       msh%dx3q = msh%dx3**2
       do k = 0, msh%nz+1
          msh%zm(k) = msh%zmin+(k-0.5d0)*msh%dz
       end do

       do j = msh%iirc1-1,msh%ifrc1+1
          msh%dmr(j) = 1.d0
       end do
       do j = msh%iirc2-1,msh%ifrc2+1
          msh%dcr(j) = 1.d0
       end do

       do j = 0, msh%nr+1
          msh%rm(j) = msh%rmin+(j-0.5d0)*msh%dr
       end do

       do i = 1,msh%ntheta
          msh%thm(i) = msh%thmin+(i-0.5d0)*msh%dtheta
       end do

       do k = 0,msh%nz+1
          msh%zc(k) = msh%zmin+(k-1.d0)*msh%dz
       end do

       do j = 0,msh%nr+1
          msh%rc(j) = msh%rmin+(j-1.d0)*msh%dr
       end do

       do i = msh%iithetac1,msh%ifthetac1
          msh%thc(i) = msh%thmin+(i-1.d0)*msh%dtheta
       end do


!      rzplane is .TRUE.

       allocate(msh%rz_ig_start(com%np_a(2)*com%np_a(3)), msh%ntrz(3,com%np_a(1)*com%np_a(2)*com%np_a(3)))
       msh%rz_rk_neigh = -1
       if (com%ip_a(2).gt.0) then
          msh%rz_rk_neigh(1) = (com%ip_a(2)-1)*com%np_a(3) + com%ip_a(3)
       end if
       if (com%ip_a(2).lt.com%np_a(2)-1) then
          msh%rz_rk_neigh(2) = (com%ip_a(2)+1)*com%np_a(3) + com%ip_a(3)
       end if
       if (com%ip_a(3).gt.0) then
          msh%rz_rk_neigh(3) = com%ip_a(2)*com%np_a(3) + com%ip_a(3)-1
       end if
       if (com%ip_a(3).lt.com%np_a(3)-1) then
          msh%rz_rk_neigh(4) = com%ip_a(2)*com%np_a(3) + com%ip_a(3)+1
       end if

       call MPI_Allgather((/msh%ntheta,msh%nr,msh%nz/), 3, MPI_INTEGER, msh%ntrz, 3, MPI_INTEGER, com%comm_0, code)
!       call MPI_Allgather((/msh%ntheta,msh%nr,msh%nz/), 3, MPI_INTEGER, msh%ntrz, 3, MPI_INTEGER, com%comm_0, code)
!       call MPI_Allgather(c, 1, MPI_INTEGER, b, 1, MPI_INTEGER, com%comm_0, code)
!       one=1
!       call MPI_ALLGATHER(c, 1, MPI_INTEGER, b, 1, MPI_INTEGER, com%comm_0, code)

!       MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, COMM, IERROR)

       msh%rz_ig_start=0
       do k = 1, com%np_a(2)*com%np_a(3)-1
          msh%rz_ig_start(k+1) = msh%rz_ig_start(k) + msh%ntrz(2,k)*msh%ntrz(3,k)
       end do
       msh%rz_ig_start_sc = msh%rz_ig_start(com%ip_a(2)*com%np_a(3)+com%ip_a(3)+1)

!      rzplane was .TRUE.



    else
       print*, 'failure in create_mesh'
    end if

    close(22)


  end subroutine create_mesh_a


end module mesh_tools
