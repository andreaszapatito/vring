MODULE fftw3
!#include <petsc/finclude/petscdef.h>
  use comm_tools
  USE, INTRINSIC :: iso_c_binding
  implicit none
!#include <mpif.h>
#include <fftw3.f03>

CONTAINS
  SUBROUTINE fft_ptf_multi(tab_in, tab_out)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:) :: tab_in, tab_out
    COMPLEX(KIND=8), ALLOCATABLE :: four(:,:,:)
    INTEGER          :: N, nhalf, Nfft, Nd2, Nd3
    INTEGER          :: k
    INTEGER          :: istride, ostride, idist, odist, inembed, onembed
    INTEGER(KIND=8)  :: plan
    REAL(KIND=8)     :: t0, t1

    N = SIZE(tab_in,1)
    nhalf = N/2+1
    Nd2=SIZE(tab_in,2)
    Nd3=SIZE(tab_in,3)
    ALLOCATE(four(nhalf,Nd2,Nd3))
    Nfft=Nd2*Nd3
    idist = N
    odist = nhalf
    istride = 1
    ostride = 1
    inembed = N
    onembed = nhalf


    CALL dfftw_plan_many_dft_r2c(plan,1, (/N/), Nfft, tab_in, &
         inembed, istride, idist, four, onembed, ostride, odist, FFTW_ESTIMATE)

    CALL dfftw_execute(plan)

    tab_out(1,:,:) = REAL(four(1,:,:))/N
    DO k = 1, nhalf-1
       tab_out(2*k,:,:)   = 2.d0*REAL (four(k+1,:,:),KIND=8)/N
       tab_out(2*k+1,:,:) = -2.d0*AIMAG(four(k+1,:,:))/N
    END DO

    DEALLOCATE(four)
    CALL dfftw_destroy_plan(plan)

  END SUBROUTINE fft_ptf_multi


  SUBROUTINE fft_ftp_multi(tab_in, tab_out)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:,:) :: tab_in, tab_out
    COMPLEX(KIND=8), ALLOCATABLE :: four(:,:,:)
    INTEGER          :: N, nhalf, Nfft, Nd2, Nd3
    INTEGER          :: k
    INTEGER          :: istride, ostride, idist, odist, inembed, onembed
    INTEGER(KIND=8)  :: plan
    REAL(KIND=8)     :: t0, t1

    N = SIZE(tab_in,1)
    nhalf = N/2+1
    Nd2=SIZE(tab_in,2)
    Nd3=SIZE(tab_in,3)
    ALLOCATE(four(nhalf,Nd2,Nd3))
    Nfft=Nd2*Nd3
    idist = N
    odist = nhalf
    istride = 1
    ostride = 1
    inembed = N
    onembed = nhalf

    four(1,:,:) = tab_in(1,:,:)
    DO k = 1, nhalf-1
       four(k+1,:,:) = CMPLX(tab_in(2*k,:,:),-1.d0*tab_in(2*k+1,:,:),KIND=8)/2.d0
    END DO

    CALL dfftw_plan_many_dft_c2r(plan,1, (/N/), Nfft, four, onembed, ostride, odist,&
         tab_out, inembed, istride, idist, FFTW_ESTIMATE)

    CALL dfftw_execute(plan)

    DEALLOCATE(four)
    CALL dfftw_destroy_plan(plan)

  END SUBROUTINE fft_ftp_multi

  SUBROUTINE fft_cos_multi_forward(tab_in, tab_out, comm1d)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: tab_in, tab_out
    integer                      :: comm1d

    INTEGER                      :: nr, nz
    REAL(KIND=8), ALLOCATABLE    :: tmp(:,:), tmpf(:,:), trsfr(:)
    INTEGER(KIND=8)  :: plan

    INTEGER          :: nprocs, nrloc, rk, code
    INTEGER          :: k, j, i, ig

    nr = SIZE(tab_in, 1)
    nz = SIZE(tab_in, 2)

    CALL MPI_COMM_SIZE(comm1d, nprocs, code)
    CALL MPI_COMM_RANK(comm1d, rk,     code)

    IF (MOD(nr,nprocs)/=0) THEN
       nrloc = nr/nprocs+1
    ELSE
       nrloc = nr/nprocs
    END IF
    ALLOCATE(trsfr(nprocs*nz*nrloc))
    trsfr = 0.d0
  
    DO k = 0, nprocs-1
       DO j = 1, nz
          DO i = 1, nrloc
             ig = k*nrloc + i
             IF (ig.LE.nr) THEN
                trsfr( k*nrloc*nz + (j-1)*nrloc + i ) = tab_in(ig,j)
             END IF
          END DO
       END DO
    END DO

    CALL MPI_ALLTOALL(trsfr, nrloc*nz, MPI_DOUBLE_PRECISION, trsfr, nrloc*nz, &
         MPI_DOUBLE_PRECISION, comm1d, code)
    
    ALLOCATE(tmp(nz*nprocs,nrloc), tmpf(nz*nprocs,nrloc))
    DO j = 1, nz*nprocs
       DO i = 1, nrloc
          tmp(j,i) = trsfr( (j-1)*nrloc+i )
       END DO
    END DO


    CALL dfftw_plan_many_r2r(plan, 1, (/nz*nprocs/), nrloc, tmp, nz*nprocs, 1, nz*nprocs, &
         tmpf, nz*nprocs, 1, nz*nprocs, FFTW_REDFT10, FFTW_ESTIMATE)
    CALL dfftw_execute(plan)
    

    DO j = 1, nz*nprocs
       DO i = 1, nrloc
          trsfr( (j-1)*nrloc+i ) = tmpf(j,i)
       END DO
    END DO

    CALL MPI_ALLTOALL(trsfr, nrloc*nz, MPI_DOUBLE_PRECISION, trsfr, nrloc*nz, &
         MPI_DOUBLE_PRECISION, comm1d, code)


    DO k = 0, nprocs-1
       DO j = 1, nz
          DO i = 1, nrloc
             ig = k*nrloc + i
             IF (ig.LE.nr) THEN
                tab_out(ig,j) =  trsfr( k*nrloc*nz + (j-1)*nrloc + i )
             END IF
          END DO
       END DO
    END DO

    DEALLOCATE(tmp, tmpf, trsfr)
    CALL dfftw_destroy_plan(plan)   
    tab_out = tab_out/(2.d0*nz*nprocs)


  END SUBROUTINE fft_cos_multi_forward

  SUBROUTINE fft_cos_multi_backward(tab_in, tab_out, comm1d)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: tab_in, tab_out
    integer                      :: comm1d
    INTEGER                      :: nr, nz
    REAL(KIND=8), ALLOCATABLE    :: tmp(:,:), tmpf(:,:), trsfr(:)
    INTEGER(KIND=8)  :: plan
    INTEGER          :: nrloc, nprocs, rk, ig, i, j, k, code

    nr = SIZE(tab_in, 1)
    nz = SIZE(tab_in, 2)


    CALL MPI_COMM_SIZE(comm1d, nprocs, code)
    CALL MPI_COMM_RANK(comm1d, rk,     code)

    IF (MOD(nr,nprocs)/=0) THEN
       nrloc = nr/nprocs+1
    ELSE
       nrloc = nr/nprocs
    END IF
    ALLOCATE(trsfr(nprocs*nz*nrloc))
    trsfr = 0.d0
  
    DO k = 0, nprocs-1
       DO j = 1, nz
          DO i = 1, nrloc
             ig = k*nrloc + i
             IF (ig.LE.nr) THEN
                trsfr( k*nrloc*nz + (j-1)*nrloc + i ) = tab_in(ig,j)
             END IF
          END DO
       END DO
    END DO

    CALL MPI_ALLTOALL(trsfr, nrloc*nz, MPI_DOUBLE_PRECISION, trsfr, nrloc*nz, &
         MPI_DOUBLE_PRECISION, comm1d, code)
    
    ALLOCATE(tmp(nz*nprocs,nrloc), tmpf(nz*nprocs,nrloc))
    DO j = 1, nz*nprocs
       DO i = 1, nrloc
          tmp(j,i) = trsfr((j-1)*nrloc+i )
       END DO
    END DO


    CALL dfftw_plan_many_r2r(plan, 1, (/nz*nprocs/), nrloc, tmp, nz*nprocs, 1, nz*nprocs,&
         tmpf, nz*nprocs, 1, nz*nprocs, FFTW_REDFT01, FFTW_ESTIMATE)
    CALL dfftw_execute(plan)


    DO j = 1, nz*nprocs
       DO i = 1, nrloc
          trsfr( (j-1)*nrloc+i ) = tmpf(j,i)
       END DO
    END DO

    CALL MPI_ALLTOALL(trsfr, nrloc*nz, MPI_DOUBLE_PRECISION, trsfr, nrloc*nz, &
         MPI_DOUBLE_PRECISION, comm1d, code)


    DO k = 0, nprocs-1
       DO j = 1, nz
          DO i = 1, nrloc
             ig = k*nrloc + i
             IF (ig.LE.nr) THEN
                tab_out(ig,j) =  trsfr( k*nrloc*nz + (j-1)*nrloc + i )
             END IF
          END DO
       END DO
    END DO

    
    DEALLOCATE(tmp, tmpf, trsfr)
    CALL dfftw_destroy_plan(plan)  


  END SUBROUTINE fft_cos_multi_backward


  SUBROUTINE fft_cos_multi_2ways(tab_in, tab_out, comm1d, fft_dir)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: tab_in, tab_out
    integer                      :: comm1d
    INTEGER                      :: fft_dir
    INTEGER                      :: nr, nz
    REAL(KIND=8), ALLOCATABLE    :: tmp(:,:), tmpf(:,:), trsfr(:), trsfr_out(:)
    INTEGER(KIND=8)  :: plan
    INTEGER          :: nrloc, nprocs, rk, ig, i, j, k, code

    nr = SIZE(tab_in, 1)
    nz = SIZE(tab_in, 2)


    CALL MPI_COMM_SIZE(comm1d, nprocs, code)
    CALL MPI_COMM_RANK(comm1d, rk,     code)

    IF (MOD(nr,nprocs)/=0) THEN
       nrloc = nr/nprocs+1
    ELSE
       nrloc = nr/nprocs
    END IF
    ALLOCATE(trsfr(nprocs*nz*nrloc), trsfr_out(nprocs*nz*nrloc))
    trsfr = 0.d0
  

    DO k = 0, nprocs-1
       DO j = 1, nz
          DO i = 1, nrloc
             ig = k*nrloc + i
             IF (ig.LE.nr) THEN
                trsfr( k*nrloc*nz + (j-1)*nrloc + i ) = tab_in(ig,j)
             END IF
          END DO
       END DO
    END DO
    ALLOCATE(tmp(nz*nprocs,nrloc), tmpf(nz*nprocs,nrloc))

    CALL MPI_ALLTOALL(trsfr, nrloc*nz, MPI_DOUBLE_PRECISION, trsfr_out, nrloc*nz, &
         MPI_DOUBLE_PRECISION, comm1d, code)
    
    DO j = 1, nz*nprocs
       DO i = 1, nrloc
          tmp(j,i) = trsfr_out((j-1)*nrloc+i )
       END DO
    END DO


    IF (fft_dir==1) THEN
       CALL dfftw_plan_many_r2r(plan, 1, (/nz*nprocs/), nrloc, tmp, nz*nprocs, 1, nz*nprocs,&
            tmpf, nz*nprocs, 1, nz*nprocs, FFTW_REDFT10, FFTW_ESTIMATE)
    ELSE
       CALL dfftw_plan_many_r2r(plan, 1, (/nz*nprocs/), nrloc, tmp, nz*nprocs, 1, nz*nprocs,&
            tmpf, nz*nprocs, 1, nz*nprocs, FFTW_REDFT01, FFTW_ESTIMATE)
    END IF
    CALL dfftw_execute(plan)
!    tmpf=tmp
    
    trsfr_out = 0.d0
    DO j = 1, nz*nprocs
       DO i = 1, nrloc
          trsfr_out( (j-1)*nrloc+i ) = tmpf(j,i)
       END DO
    END DO

    CALL MPI_ALLTOALL(trsfr_out, nrloc*nz, MPI_DOUBLE_PRECISION, trsfr, nrloc*nz, &
         MPI_DOUBLE_PRECISION, comm1d, code)


    DO k = 0, nprocs-1
       DO j = 1, nz
          DO i = 1, nrloc
             ig = k*nrloc + i
             IF (ig.LE.nr) THEN
                tab_out(ig,j) =  trsfr( k*nrloc*nz + (j-1)*nrloc + i )
             END IF
          END DO
       END DO
    END DO

    IF (fft_dir==1) THEN
       tab_out = tab_out/(2.d0*nz*nprocs)
    END IF

    
    DEALLOCATE(tmp, tmpf, trsfr,trsfr_out)
    CALL dfftw_destroy_plan(plan)  


  END SUBROUTINE fft_cos_multi_2ways



END MODULE fftw3
