module poisson
!#include <petsc/finclude/petscdef.h>
!  use petsc
  use cons_type  
  use para_type  
  use comm_type  
  use mesh_type 
  use solp_type
  use gen_tools
  use run_type
  use run_tools
!  use petsc_mat
!  use fftw3

  implicit none
!#include <mpif.h>
contains

  subroutine prep_poisson_solver(var)
    implicit none
    type(run)        :: var


    integer      :: k, nhalf, i

    select case(var%par%poisson_type)
!    case('fft1')
!       nhalf = var%msh%ntheta/2
!       var%sp%nsys=nhalf+1
!       call init_solvers(var%sp)
!       call verbose(var,"init_solvers called")
!
!       do k = 0, nhalf
!          call init_specific_solver(var%com,var%sp,1+k,0)
!          call assign_matrix_3df_cyl(k+1, var%com, var%msh, var%sp, 'MatLapCylNModif', k)
!!!$          call assign_matrix(k+1, comm, com, msh, 'matlapcylnmodif', k)
!       end do
!    case('fft2')
!       nhalf = var%msh%ntheta/2
!       var%sp%nsys=var%msh%nz*(nhalf+1)
!       call init_solvers(var%sp)
!       
!       do i = 0, nhalf
!          do k = 1, var%msh%nz
!             call init_specific_solver(var%com,var%sp,i*var%msh%nz+k,2)
!             call assign_matrix_3dff_cyl(i*var%msh%nz+k,var%com,var%msh,var%sp, (/i, &
!                  var%com%ip_a(3)*var%msh%nz+k-1/))
!          end do
!       end do
    case('fftp')
       call init_tripar_solvers(var%msh,var%sp)
    case('btri')
       call init_btri_solver(var%msh,var%sp)
    end select
       

  end subroutine prep_poisson_solver

!  subroutine init_solvers(sp)
!    implicit none
!    type(solp) sp
!    allocate(sp%ksp  (sp%nsys))
!    allocate(sp%used (sp%nsys))
!    allocate(sp%vx   (sp%nsys))
!    allocate(sp%vb   (sp%nsys))
!    sp%used=.false.
!
!  end subroutine init_solvers
!
!
!  subroutine init_specific_solver(com,sp,nn,icom)
!    implicit none
!    integer              :: nn,icom
!    integer              :: com1
!    integer              :: rank, code
!    integer              :: ierr
!    type (comm)          :: com
!    type (solp)          :: sp
!
!    if (icom.eq.0) com1=com%comm_0
!    if (icom.gt.0) com1=com%comm_a(icom)
!
!    if ( (nn.gt.sp%nsys) .or. (nn.lt.1) ) then
!       call warning('solver not initialized, wrong index', (/com%ip/))
!       return
!    end if
!
!    if (sp%used(nn)) then
!       print*, 'warning, replacing solver', nn
!       call kspdestroy(sp%ksp(nn),ierr)
!
!       call vecdestroy(sp%vx(nn), ierr)
!       call vecdestroy(sp%vb(nn), ierr)
!
!    end if
!
!    call kspcreate(com1,sp%ksp(nn),ierr)
!    call veccreate(com1,sp%vx(nn), ierr)
!    call veccreate(com1,sp%vb(nn), ierr)
!    sp%used(nn)=.true.
!
!
!  end subroutine init_specific_solver


!  subroutine assign_matrix_3d_cart(nn,com1,sp,com,msh,typemat,options)
!    implicit none
!    integer              :: nn
!    MPI_Comm             :: com1
!    type(comm)           :: com
!    type(mesh_a)         :: msh
!    type(solp)           :: sp
!    character(len=*)     :: typemat
!    integer, optional    :: options(:)
!
!    print*, 'oups, not yet available'
!
!  end subroutine assign_matrix_3d_cart

  subroutine assign_matrix_3dff_cyl(nn,com,msh,sp,modes)
    implicit none
    integer              :: nn
    integer              :: com1
    type(comm)           :: com
    type(mesh_a)         :: msh
    type(solp)           :: sp
    integer, optional    :: modes(2)

 !   Mat                  :: matrix
 !   PC                   :: prec
 !   integer              :: nloc, nglob
 !   integer              :: idxn(3)
 !   real(kind=8)         :: valm(3), valx(3)
 !   PetscErrorCode       :: ierr
 !   real(kind=8)         :: mth, mz, mult,  pi
 !   real(kind=8)         :: rr
 !   integer              :: j, k, jg, kg, idxg, ll, jc, kc, rkc
 !   integer              :: rk2d

 !   com1=com%comm_a(2)

 !   call matcreate(com1, matrix, ierr)
 !   nloc  = msh%nr
 !   nglob = msh%nrg
 !   call matsetsizes(matrix, nloc, nloc, nglob, nglob, ierr)
 !   call matsettype(matrix, MATMPIAIJ,ierr)
 !   call matmpiaijsetpreallocation(matrix, 3, petsc_null_integer, 1, &
 !        petsc_null_integer, ierr)
 !   call matzeroentries(matrix, ierr)


 !   mth = modes(1)**2
 !   mz  = (2.d0*acos(-1.d0)*modes(2)/2.d0*msh%zmaxg)**2

 !   do j = 1, msh%nr
 !      jg = com%ip_a(2)*msh%nr+j-1
 !      idxn = (/ jg, jg-1, jg+1 /)
 !      if (jg+1.ge.msh%nrg) then
 !         idxn(3) = -1
 !      end if
 !      rr = msh%rm(j)
 !      valx = (/ 2.d0/msh%dr**2 + mz + mth/rr**2, &
 !           -1.d0/msh%dr**2 + 1.d0/(2.d0*rr*msh%dr), &
 !           -1.d0/msh%dr**2 - 1.d0/(2.d0*rr*msh%dr) /)

 !      if (jg==0) then
 !         valx(1) = valx(1)+valx(2)
 !      end if
 !      if (jg==msh%nrg-1) then
 !         if ( (modes(2)==0) .and. (modes(1)==0) ) then
 !            valx(1) = valx(1) - valx(3)
 !         else
 !            valx(1) = valx(1) + valx(3)
 !         end if
 !      end if

 !      call matsetvalues(matrix, 1, jg, 3, idxn, valx, insert_values, ierr)
 !   end do
 !   call matassemblybegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
 !   call matassemblyend(matrix, MAT_FINAL_ASSEMBLY, ierr)

 !   call kspsetoperators(sp%ksp(nn), matrix, matrix, ierr)
 !   call kspgetpc(sp%ksp(nn), prec, ierr)
 !   call pcsettype(prec, PCLU, ierr)
 !   call kspsettype(sp%ksp(nn), KSPPREONLY, ierr)
 !   call pcfactorsetmatsolverpackage(prec, MATSOLVERMUMPS, ierr)
!!    call pcfactorsetmatsolverpackage(prec, superlu, ierr)

 !   call kspsetfromoptions(sp%ksp(nn),ierr)

 !   call vecsetsizes(sp%vx(nn), nloc, nglob, ierr)
 !   call vecsetsizes(sp%vb(nn), nloc, nglob, ierr)
 !   call vecsetfromoptions(sp%vx(nn), ierr)
 !   call vecsetfromoptions(sp%vb(nn), ierr)

 !   call matdestroy(matrix,ierr)


  end subroutine assign_matrix_3dff_cyl

  subroutine assign_matrix_3df_cyl(nn,com,msh,sp,typemat,mode)
    implicit none
    integer              :: nn
    integer              :: com1
    type(comm)           :: com
    type(solp)           :: sp
    type(mesh_a)         :: msh
    character(len=*)     :: typemat
    integer, optional    :: mode

!    Mat                  :: matrix
!    PC                   :: prec
!    integer              :: nloc, nglob
!    integer              :: idxn(5)
!    real(kind=8)         :: valm(5), valx(5)
!    PetscErrorCode       :: ierr
!    real(kind=8)         :: mult, pi
!    real(kind=8)         :: rr
!    integer              :: j, k, jg, kg, idxg, ll, jc, kc, rkc
!    integer              :: rk2d
!
!    com1=com%comm_0
!    call matcreate(com1, matrix, ierr)
!    nloc  = msh%nr*msh%nz
!    nglob = msh%nrg*msh%nzg
!    call matsetsizes(matrix, nloc, nloc, nglob, nglob, ierr)
!    call matsettype(matrix, MATMPIAIJ,ierr)
!    call matmpiaijsetpreallocation(matrix, 5, petsc_null_integer, 2, &
!         petsc_null_integer, ierr)
!    call matzeroentries(matrix, ierr)
!
!
!    pi = acos(-1.d0)
!    if (trim(adjustl(typemat))=='MatLapCylNModif') then
!       mult = 2.d0/msh%dtheta**2*(1.d0-cos(2.d0*pi*real(mode)/msh%nthetag))
!    else
!      mult = mode**2
!    end if
!
!
!    do k = 1, msh%nz
!       do j = 1, msh%nr
!          rr = msh%rm(j)
!          valm = (/2.d0/msh%dz**2+2.d0/msh%dr**2+mult/rr**2, &
!               -1.d0/msh%dr**2+1.d0/(2.d0*rr*msh%dr), -1.d0/msh%dr**2-1.d0/(2.d0*rr*msh%dr),&
!               -1.d0/msh%dz**2, -1.d0/msh%dz**2/)
!
!          valm = (/2.d0/msh%dz**2+2.d0/msh%dr**2+mult/rr**2, &
!               -1.d0/msh%dr**2+1.d0/(2.d0*rr*msh%dr), -1.d0/msh%dr**2-1.d0/(2.d0*rr*msh%dr),&
!               -1.d0/msh%dz**2, -1.d0/msh%dz**2/)
!          idxn=-1
!          idxg = msh%rz_ig_start_sc + (k-1)*msh%nr + j-1
!          idxn(1) = idxg
!          rk2d = com%ip_a(2)*com%np_a(3) + com%ip_a(3)
!
!
!          jg = com%ip_a(2)*msh%nr + j
!          kg = com%ip_a(3)*msh%nz + k
!
!          !left
!          jc = j-1
!          kc = k
!          rkc = rk2d
!          if (jc==0) then
!             rkc = msh%rz_rk_neigh(1)
!             if (rkc.ge.0) then
!                jc = msh%ntrz(2,rkc+1)
!             end if
!          end if
!          if (rkc.ge.0) then
!             idxn(2) = msh%rz_ig_start(rkc+1) + (kc-1)*msh%ntrz(2,rkc+1) + jc-1
!          else
!             idxn(2) = -1
!          end if
!
!          !right
!          jc = j+1
!          kc = k
!          rkc = rk2d
!          if (jc==msh%nr+1) then
!             rkc = msh%rz_rk_neigh(2)
!             jc = 1
!          end if
!          if (rkc.ge.0) then
!             idxn(3) = msh%rz_ig_start(rkc+1) + (kc-1)*msh%ntrz(2,rkc+1) + jc-1
!          else
!             idxn(3) = -1
!          end if
!
!          !down
!          jc = j
!          kc = k-1
!          rkc=rk2d
!          if (kc==0) then
!             rkc = msh%rz_rk_neigh(3)
!             if (rkc.ge.0) then
!                kc = msh%ntrz(3,rkc+1)
!             end if
!          end if
!          if (rkc.ge.0) then
!             idxn(4) = msh%rz_ig_start(rkc+1) + (kc-1)*msh%ntrz(2,rkc+1) + jc-1
!          else
!             idxn(4) = -1
!          end if
!
!          !up
!          jc = j
!          kc = k+1
!          rkc=rk2d
!          if (kc==msh%nz+1) then
!             rkc = msh%rz_rk_neigh(4)
!             kc = 1
!          end if
!          if (rkc.ge.0) then
!             idxn(5) = msh%rz_ig_start(rkc+1) + (kc-1)*msh%ntrz(2,rkc+1) + jc-1
!          else
!             idxn(5) = -1
!          end if
!
!          if (jg==1) then
!             valm(1) = valm(1)+valm(2)
!          end if
!          if (jg==msh%nrg) then
!             valm(1) = valm(1)+valm(3)
!          end if
!          if (kg==1) then
!             valm(1) = valm(1)+valm(4)
!          end if
!          if (kg==msh%nzg) then
!             valm(1) = valm(1)+valm(5)
!          end if
!
!          if ( (jg==1).and.(kg==1).and.(mode==0) ) then
!             valm = (/2.d0/msh%dz**2+2.d0/msh%dr**2           ,  -1.d0/msh%dr**2+1.d0/(2.d0*rr*msh%dr), -1.d0/msh%dr**2-1.d0/(2.d0*rr*msh%dr), -1.d0/msh%dz**2, -1.d0/msh%dz**2/)
!              valm(1) = valm(1) - valm(4)
!!            valm(1) =         + valm(3) + valm(5)
!!             valm(1) = 0.0                         
!          end if
!
!          call matsetvalues(matrix,1,idxg,5,idxn,valm,insert_values,ierr)
!       end do
!    end do
!
!
!    call matassemblybegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
!    call matassemblyend(matrix, MAT_FINAL_ASSEMBLY, ierr)
!    call KSPSetOperators(sp%ksp(nn), matrix, matrix, ierr)
!!!$ call kspgetpc(sp%ksp(nn), prec, ierr)
!!!$ call pcsettype(prec, pclu, ierr)
!!!$ call kspsettype(sp%ksp(nn), ksppreonly, ierr)
!!!$ call pcfactorsetmatsolverpackage(prec, matsolverpastix, ierr)
!
!!!$ call pcsettype(prec, pcjacobi, ierr)
!    call kspsettype(sp%ksp(nn),KSPBCGS, ierr)
!    call kspsettolerances(sp%ksp(nn), 1.d-16, 1.d-16, PETSC_DEFAULT_REAL, 10000, ierr)
!    call kspsetfromoptions(sp%ksp(nn),ierr)
!    call vecsetsizes(sp%vx(nn), nloc, nglob, ierr)
!    call vecsetsizes(sp%vb(nn), nloc, nglob, ierr)
!    call vecsetfromoptions(sp%vx(nn), ierr)
!    call vecsetfromoptions(sp%vb(nn), ierr)
!
!    call matdestroy(matrix,ierr)


  end subroutine assign_matrix_3df_cyl


  subroutine init_tripar_solvers(msh,sp)
    implicit none
    type(mesh_a)  :: msh
    type(solp)           :: sp

    if (allocated(sp%bb_base)) then
       deallocate(sp%bb_base,sp%aa_base,sp%cc_base)
    end if

    allocate(sp%bb_base(msh%nr), sp%aa_base(msh%nr), sp%cc_base(msh%nr))

    sp%bb_base(:) = 2.d0/msh%dr**2
    sp%aa_base(:) = -1.d0/msh%dr**2 + 1.d0/(2.d0*msh%rm(1:msh%nr)*msh%dr)
    sp%cc_base(:) = -1.d0/msh%dr**2 - 1.d0/(2.d0*msh%rm(1:msh%nr)*msh%dr)
!    write (*,*) 'initialising base', sp%cc_base(1),-1.d0/(msh%dr)**2.0 - 1.d0/(2.d0*msh%rm*msh%dr)



    sp%meth_tripar = 'wvnex'

  end subroutine init_tripar_solvers

  subroutine init_btri_solver(msh,sp)
    implicit none
    type(mesh_a)           :: msh
    type(solp)           :: sp
    integer                :: k

    if (allocated(sp%bb_base)) then
       deallocate(sp%bb_base,sp%aa_base,sp%cc_base)
    end if

    allocate(sp%bb_base(msh%nr*msh%nz), sp%aa_base(msh%nr*msh%nz), &
         sp%cc_base(msh%nr*msh%nz))

    sp%bb_base = 2.d0/msh%dr**2 + 2.d0/msh%dz**2
    do k = 1, msh%nz
       sp%aa_base((k-1)*msh%nr+1:k*msh%nr) = -1.d0/msh%dr**2 + 1.d0/(2.d0*msh%rm*msh%dr)
       sp%cc_base((k-1)*msh%nr+1:k*msh%nr) = -1.d0/msh%dr**2 - 1.d0/(2.d0*msh%rm*msh%dr)
    end do

    sp%meth_tripar = 'wvnex'
    sp%btri_mult   = -1.d0*msh%dz**2

  end subroutine init_btri_solver



end module poisson
