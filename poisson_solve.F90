module poisson_solve
!#include <petsc/finclude/petscdef.h>
  use poisson
  use comm_type
  use mesh_type
  use fftw3
!  use user_mod
  implicit none
!#include <mpif.h>

contains

  subroutine laplace_solve(sp,com,msh,par)
    implicit none
    type(solp)                     :: sp
    type(comm)                     :: com
    type(mesh_a)                   :: msh
    type(para)                     :: par

    if      (par%poisson_type=='fft1') then
!       call solve_lap_fft1(sp,com,msh)
    else if (par%poisson_type=='fft2') then
!       call solve_lap_fft2(sp,com,msh)
    else if (par%poisson_type=='fftp') then
       call solve_lap_tripar(sp,com,msh)
    else
!       call solve_lap_btri(sp,com,msh)
    end if
    sp%dph(1:msh%ntheta,1:msh%nr,1:msh%nz)=sp%sol(1:msh%ntheta,1:msh%nr,1:msh%nz)

    
  end subroutine laplace_solve

!  subroutine solve_lap_fft2(sp,com,msh)
!    implicit none
!    type(comm)                     :: com
!    type(mesh_a)                   :: msh
!    type(solp)                     :: sp
!
!    real(kind=8), allocatable      :: act_rhs(:,:,:), act_sol(:,:)
!    integer                        :: j, k, m, nhalf, nrl, nrt
!    real(kind=8)                   :: t0, t1
!    PetscErrorCode        :: ierr
!
!    allocate(act_rhs(msh%ntheta,msh%nr,msh%nz))
!    act_rhs(1:msh%ntheta,1:msh%nr,1:msh%nz) = sp%rhs(1:msh%ntheta,1:msh%nr,1:msh%nz)
!    nhalf = msh%ntheta/2
!    call fft_ptf_multi(act_rhs, act_rhs)
!
!
!    do m = 1, msh%ntheta
!       call fft_cos_multi_2ways(act_rhs(m,:,:),act_rhs(m,:,:),com%comm_a(3),1)
!    end do
!
!
!
!    do k = 1, msh%nz
!       call solve_ksp_rk(k,com%ip_a(2), act_rhs(1,:,k),sp%sol(1,:,k),sp)
!    end do
!
!    do m = 1, nhalf
!       do k = 1, msh%nz
!          call mpi_barrier(com%comm_0, ierr)
!
!          call solve_ksp_rk(m*msh%nz+k, com%ip_a(2), act_rhs(2*m,:,k),sp%sol(2*m,:,k),sp)
!          call solve_ksp_rk(m*msh%nz+k, com%ip_a(2), act_rhs(2*m+1,:,k), sp%sol(2*m+1,:,k),sp)
!       end do
!    end do
!
!    do m = 1, msh%ntheta
!       call fft_cos_multi_2ways(sp%sol(m,:,:),sp%sol(m,:,:),com%comm_a(3),-1)
!    end do
!    
!    call fft_ftp_multi(sp%sol,sp%sol)
!
!    deallocate(act_rhs)
!
!  end subroutine solve_lap_fft2
!
!  subroutine solve_lap_fft1(sp,com,msh)
!    implicit none
!    type(comm)                     :: com
!    type(mesh_a)                   :: msh
!    type(solp)                     :: sp
!
!    real(kind=8), allocatable      :: act_rhs(:,:,:), act_sol(:,:,:)
!    integer                        :: j, k, m, nhalf
!    real(kind=8)                   :: t0, t1
!    allocate(act_rhs(size(sp%rhs,1),size(sp%rhs,2),size(sp%rhs,3)))
!    act_rhs = sp%rhs
!
!    nhalf = (msh%ntheta)/2
!    call fft_ptf_multi(act_rhs,act_rhs)
!
!    !!mode 0
!    call solve_ksp_2d(1,com,msh,sp,act_rhs(1,:,:),sp%sol(1,:,:))
!
!    !!autres modes
!    do m = 1, nhalf
!       call solve_ksp_2d(m+1,com,msh,sp,act_rhs(2*m,:,:)  ,sp%sol(2*m,:,:))
!       call solve_ksp_2d(m+1,com,msh,sp,act_rhs(2*m+1,:,:),sp%sol(2*m+1,:,:))
!    end do
!
!    call fft_ftp_multi(sp%sol,sp%sol)
!    deallocate(act_rhs)
!
!  end subroutine solve_lap_fft1

!  subroutine solve_lap_btri(sp,com,msh)
!    implicit none
!    type(comm)                     :: com
!    type(mesh_a)                   :: msh
!    type(solp)                     :: sp
!
!    real(kind=8), allocatable      :: act_rhs(:,:,:), act_sol(:,:,:)
!    integer                        :: j, k, m, nhalf
!    real(kind=8)                   :: t0, t1
!
!    allocate(act_rhs(size(sp%rhs,1),size(sp%rhs,2),size(sp%rhs,3)))
!    act_rhs = sp%rhs
!
!    nhalf = msh%ntheta/2
!     call fft_ptf_multi(act_rhs, act_rhs)
!
!    !!mode 0
!    call solve_btri(0, com, msh, sp, act_rhs(1,:,:), sp%sol(1,:,:))
!
!    !!autres modes
!    do m = 1, nhalf
!       call solve_btri(m, com, msh, sp, act_rhs(2*m,:,:), sp%sol(2*m,:,:))
!       call solve_btri(m, com, msh, sp, act_rhs(2*m+1,:,:), sp%sol(2*m+1,:,:))
!    end do
!
!    call fft_ftp_multi(sp%sol,sp%sol)
!    deallocate(act_rhs)
!
!  end subroutine solve_lap_btri


  subroutine solve_lap_tripar(sp,com,msh)
    implicit none
    type(comm)                     :: com
    type(mesh_a)                   :: msh
    type(solp)                     :: sp
    real(kind=8), allocatable      :: act_rhs(:,:,:)

    integer                        :: j, k, m, nhalf, nrl, nrt, mz
    real(kind=8)                   :: t0, t1

    allocate(act_rhs(size(sp%rhs,1),size(sp%rhs,2),size(sp%rhs,3)))
    act_rhs = sp%rhs

    nhalf = msh%ntheta/2
    call fft_ptf_multi(act_rhs, act_rhs)


    do m = 1, msh%ntheta
       call fft_cos_multi_2ways(act_rhs(m,:,:),act_rhs(m,:,:),com%comm_a(3),1)
    end do


    do k = 1, msh%nz
     mz = com%ip_a(3)*msh%nz+k-1
     call solve_radial( (/0,mz/), com, msh, sp, act_rhs(1,:,k), sp%sol(1,:,k), (/1,-1/))
      do m = 1, nhalf
       call solve_radial( (/m,mz/), com, msh, sp, act_rhs(2*m,:,k), sp%sol(2*m,:,k), (/1,1/))
       call solve_radial( (/m,mz/), com, msh, sp, act_rhs(2*m+1,:,k), sp%sol(2*m+1,:,k), (/1,1/))
      end do
    end do

    do m = 1, msh%ntheta
       call fft_cos_multi_2ways(sp%sol(m,:,:),sp%sol(m,:,:),com%comm_a(3),-1)
    end do

    call fft_ftp_multi(sp%sol,sp%sol)

    deallocate(act_rhs)

  end subroutine solve_lap_tripar

!  subroutine solve_ksp_rk(nn, rank, rhs, xx,sp)
!    implicit none
!    integer               :: nn
!    integer               :: rank
!    real(kind=8)          :: rhs(:), xx(:)
!    type(solp)            :: sp
!    integer               :: k, nloc, flag
!    integer               :: nit
!    real(kind=8)          :: t0, t1, resid, toto
!    PetscErrorCode        :: ierr
!
!    nloc = size(rhs)
!    call vecsetvalues(sp%vb(nn), nloc, (/ (rank*nloc+k, k=0,nloc-1) /), rhs, insert_values, ierr)
!    call vecassemblybegin(sp%vb(nn), ierr)
!    call vecassemblyend(sp%vb(nn), ierr)
!    call vecnorm(sp%vb(nn), norm_infinity, toto, ierr)
!
!    t0 = mpi_wtime()
!    call kspsolve(sp%ksp(nn), sp%vb(nn), sp%vx(nn), ierr)
!    t1 = mpi_wtime()
!
!!$    call kspgetconvergedreason(my_ksp(nn), flag, ierr)
!!$    call kspgetiterationnumber(my_ksp(nn), nit, ierr)
!!$    call kspgetresidualnorm(my_ksp(nn), resid, ierr)
!!$    call par_print('===='//num2str(nn)//'====',(/rank/))
!!$    call par_print('norm rhs  '//rke2str(toto),(/rank/))
!!$    call par_print('ksp resid '//rke2str(resid),(/rank/))
!!$    call par_print('ksp nit   '//num2str(nit),(/rank/))
!!$    call par_print('ksp solve '//rke2str(t1-t0),(/rank/))
!!$    call par_print('ksp reason'//num2str(flag),(/rank/))
!!$    call par_print('===='//num2str(nn)//'====', (/rank/))


!    call vecgetvalues(sp%vx(nn), nloc, (/(rank*nloc+k, k=0,nloc-1) /), xx, ierr)


!  end subroutine solve_ksp_rk

!  subroutine solve_ksp_2d(nn,com,msh,sp,rhs,xx)
!    implicit none
!    integer               :: nn
!    type(comm)            :: com
!    type(mesh_a)          :: msh
!    type(solp)            :: sp
!    real(kind=8)          :: rhs(:,:), xx(:,:)
!
!    integer               :: k, nloc, flag, j, nit
!    integer               :: idxm(msh%nr)
!    PetscErrorCode        :: ierr
!    Vec :: vt
!    real(kind=8) :: toto, t0, t1, resid
!    nloc = size(rhs,1)*size(rhs,2)
!    do k = 1, msh%nz
!       idxm = (/ (msh%rz_ig_start_sc + (k-1)*msh%nr +j, j=0,msh%nr-1) /)
!
!       call vecsetvalues(sp%vb(nn), msh%nr,idxm, rhs(:,k),&
!            insert_values, ierr)
!    end do
!
!    call vecassemblybegin(sp%vb(nn), ierr)
!    call vecassemblyend(sp%vb(nn), ierr)
!
!    call vecnorm(sp%vb(nn), norm_infinity, toto, ierr)
!
!    t0 = mpi_wtime()
!    call kspsolve(sp%ksp(nn), sp%vb(nn), sp%vx(nn), ierr)
!    t1 = mpi_wtime()
!
!    call kspgetconvergedreason(sp%ksp(nn), flag, ierr)
!    call kspgetiterationnumber(sp%ksp(nn), nit, ierr)
!    call kspgetresidualnorm(sp%ksp(nn), resid, ierr)
!!!$    call par_print('====', pcomm%ranks)
!!!$    call par_print('norm rhs  '//rke2str(toto),pcomm%ranks)
!!!$    call par_print('ksp resid '//rke2str(resid),pcomm%ranks)
!!!$    call par_print('ksp nit   '//num2str(nit),pcomm%ranks)
!!!$    call par_print('ksp solve '//rke2str(t1-t0),pcomm%ranks)
!!!$    call par_print('ksp reason'//num2str(flag),pcomm%ranks)
!!!$    call par_print('====', pcomm%ranks)
!
!!!$    print*, '!!!!!flag!!!!!', flag
!    do k = 1, msh%nz
!       idxm = (/ (msh%rz_ig_start_sc + (k-1)*msh%nr +j, j=0,msh%nr-1) /)
!      call vecgetvalues(sp%vx(nn), msh%nr, idxm, xx(:,k), ierr)
!    end do
!
!
!
!  end subroutine solve_ksp_2d

  subroutine solve_radial(wavenumbers, com, msh, sp, rhs, sol, dn)
    use tripar
    implicit none
    real(kind=8), dimension(:)        :: rhs, sol
    integer,      dimension(2)        :: wavenumbers
    integer,      dimension(2)        :: dn
    type(mesh_a)                      :: msh
    type(comm)                        :: com
    type(solp)                        :: sp

    real(kind=8), dimension(size(rhs)):: aa, bb, cc
    real(kind=8)                      :: wth, wz, pi
    integer                           :: bcl, bcr,ii

    pi = acos(-1.d0)
!    write (*,*) sp%meth_tripar
    if (trim(adjustl(sp%meth_tripar))=='wvnex') then
       wth = real(wavenumbers(1))**2
       wz  = (pi*real(wavenumbers(2))/msh%zmaxg)**2
       wz  = 2.d0/msh%dz**2*(1.d0-cos(wavenumbers(2)*pi/msh%nzg)) !Francy 1907
    else
       wth = 2.d0/msh%dtheta**2*(1.d0-cos(pi*real(wavenumbers(1))/real(msh%nthetag)))
       wz  = 2.d0/msh%dz**2*(1.d0-cos(wavenumbers(2)*pi/msh%nzg))
    end if

!       wth = real(wavenumbers(1))**2
!       wz  = (pi*real(wavenumbers(2))/msh%zmaxg)**2
!       wz  = 2.d0/msh%dz**2*(1.d0-cos(wavenumbers(2)*pi/msh%nzg)) !Francy 1907
    aa = sp%aa_base
    bb = sp%bb_base + wth/msh%rm(1:msh%nr)**2 + wz
    cc = sp%cc_base
!    write (*,*) 'station 1', (aa(ii),ii=1,3)
!    write (*,*) 'station 2', (bb(ii),ii=1,3)
!    write (*,*) 'station 3', (cc(ii),ii=1,3)
    bcl = 1
    if (sum(wavenumbers)/=0) then
       bcr = 1
    else
       bcr = dn(2)
    end if


    if (com%ip_a(2)==0) then
       bb(1) = bb(1) + bcl*aa(1)
    end if
    if (com%ip_a(2)==com%np_a(2)-1) then
       bb(msh%nr) = bb(msh%nr) + bcr*cc(msh%nr)
    end if

    call solve_lap_par(rhs, sol,aa,bb,cc,msh%nrg, msh%nr, com%comm_a(2))



  end subroutine solve_radial
!
!  subroutine solve_btri(wavenumber, com, msh, sp, rhs, sol)
!    use btrithom
!    implicit none
!    real(kind=8), dimension(:,:)      :: rhs, sol
!    integer                           :: wavenumber
!    type(mesh_a)                      :: msh
!    type(comm)                        :: com
!    type(solp)                        :: sp
!
!    real(kind=8), dimension(msh%nr*msh%nz) :: b1, b2, b3, r, s
!    real(kind=8)                               :: wth
!    integer                                    :: k
!
!    b1 = sp%btri_mult*sp%aa_base
!    b3 = sp%btri_mult*sp%cc_base
!
!    do k = 1, msh%nz
!       b2((k-1)*msh%nr+1:k*msh%nr) = sp%btri_mult*(sp%bb_base((k-1)*msh%nr+1:k*msh%nr) + &
!            (1.d0*wavenumber)**2/msh%rm**2)
!
!       r((k-1)*msh%nr+1:k*msh%nr) = sp%btri_mult*rhs(:,k)
!
!    end do
!
!    !neumann in z
!    b2(1:msh%nr) = b2(1:msh%nr)+1.d0
!    b2((msh%nz-1)*msh%nr+1:msh%nr*msh%nz) = b2((msh%nz-1)*msh%nr+1:msh%nr*msh%nz) + 1.d0
!    !neumann in r
!    do k = 1, msh%nz
!       b2((k-1)*msh%nr+1) = b2((k-1)*msh%nr+1) + b1((k-1)*msh%nr+1)
!       b2(k*msh%nr)       = b2(k*msh%nr) + b3(k*msh%nr)
!    end do
!
!    if (wavenumber==0) then
!       b2(1) = b2(1)-2.d0
!    end if
!
!
!    call solve_lap_bthom(r,s,b1,b2,b3,msh%nr,msh%nz)
!
!
!    do k = 1, msh%nz
!       sol(:,k) = s((k-1)*msh%nr+1:k*msh%nr)
!    end do
!
!  end subroutine solve_btri



end module poisson_solve 
