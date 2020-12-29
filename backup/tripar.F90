MODULE tripar
  use comm_tools
  IMPLICIT NONE
!#include <petsc/finclude/petscdef.h>
!  USE petsc
!  USE triseq
!  USE trithom
!#include <mpif.h>
CONTAINS

  SUBROUTINE solve_lap_par(rhs,sol,aa,bb,cc,nmat, nloc, comm1d)
    IMPLICIT NONE
    INTEGER                         :: nmat,nloc,i,ierr
    REAL(KIND=8), DIMENSION(nloc)   :: rhs, sol, bb, aa, cc
    integer                         :: comm1d

    REAL(KIND=8), DIMENSION(nloc)   :: rhstmp
    REAL(KIND=8)                    :: al1, be1, ga1
    REAL(KIND=8)                    :: al2, be2, ga2

    REAL(KIND=8), ALLOCATABLE       :: aag(:), bbg(:), ccg(:), rhsg(:), solg(:)
    REAL(KIND=8)                    :: aloc(2), bloc(2), cloc(2), rhsloc(2)

    INTEGER       :: rank, nprocs, code, k
    REAL(KIND=8)  :: rr

    CALL MPI_COMM_RANK(comm1d, rank,   code)
    CALL MPI_COMM_SIZE(comm1d, nprocs, code)

    IF (rank==0) THEN
       aa(1) = 0.d0
    END IF
    IF (rank==nprocs-1) THEN
       cc(nloc) = 0.d0
    END IF
    
    rhstmp = rhs


    CALL solve_lap_thom(rhstmp,sol,aa,bb,cc,nloc)
    do i=1,nloc
      !write (*,*) 'solver',i,aa(i),bb(i),cc(i),rhs(i),sol(i)
    enddo
    al1 = sol(1)
    al2 = sol(nloc)
    rhstmp(1) = rhs(1)-aa(1)
    CALL solve_lap_thom(rhstmp,sol,aa,bb,cc,nloc)
    be1 = sol(1)-al1
    be2 = sol(nloc)-al2
    rhstmp(1) = rhs(1)
    rhstmp(nloc)=rhs(nloc)-cc(nloc)
    CALL solve_lap_thom(rhstmp,sol,aa,bb,cc,nloc)
    ga1 = sol(1)-al1
    ga2 = sol(nloc)-al2

    ALLOCATE(aag(2*nprocs), bbg(2*nprocs), ccg(2*nprocs), rhsg(2*nprocs), solg(2*nprocs))

    IF (rank==0) THEN
       aloc(1) = 0.d0
       bloc(1) = 1.d0
       cloc(1) = 0.d0
       rhsloc(1) = 0.d0
       aloc(2) = 0.d0
       bloc(2) = 1.d0
       cloc(2) = -1.d0*ga2
       rhsloc(2) = al2
    ELSE IF (rank==nprocs-1) THEN
       aloc(1) = -1.d0*be1
       bloc(1) = 1.d0
       cloc(1) = 0.d0
       rhsloc(1) = al1
       aloc(2) = 0.d0
       bloc(2) = 1.d0
       cloc(2) = 0.d0
       rhsloc(2) = 0.d0
    ELSE
       aloc(1) = ga1*be2-ga2*be1
       bloc(1) = ga2
       cloc(1) = -1.d0*ga1
       rhsloc(1) = ga2*al1-ga1*al2
       aloc(2) = -1.d0*be2
       bloc(2) = be1
       cloc(2) = be2*ga1-be1*ga2
       rhsloc(2) = be1*al2-be2*al1
    END IF

    CALL MPI_ALLGATHER(aloc, 2, MPI_DOUBLE_PRECISION, aag, 2, MPI_DOUBLE_PRECISION, &
         comm1d, code)
    CALL MPI_ALLGATHER(bloc, 2, MPI_DOUBLE_PRECISION, bbg, 2, MPI_DOUBLE_PRECISION, &
         comm1d, code)
    CALL MPI_ALLGATHER(cloc, 2, MPI_DOUBLE_PRECISION, ccg, 2, MPI_DOUBLE_PRECISION, &
         comm1d, code)
    CALL MPI_ALLGATHER(rhsloc, 2, MPI_DOUBLE_PRECISION, rhsg, 2, MPI_DOUBLE_PRECISION, &
         comm1d, code)

    bbg(1) = 1.d0
    ccg(1) = 0.d0
    rhsg(1) = 0.d0
    bbg(2*nprocs) = 1.d0
    aag(2*nprocs) = 0.d0
    rhsg(2*nprocs) = 0.d0

    CALL solve_lap_thom(rhsg,solg,aag,bbg,ccg,2*nprocs)

    rhstmp = rhs
    IF (rank/=0) THEN
       rhstmp(1) = rhs(1) - aa(1)*solg(2*rank)
    END IF
    IF (rank/=nprocs-1) THEN
       rhstmp(nloc) = rhs(nloc) - cc(nloc)*solg(2*rank+3)
    END IF

    CALL solve_lap_thom(rhstmp,sol,aa,bb,cc,nloc)


    DEALLOCATE(aag,bbg,ccg,rhsg,solg)

  END SUBROUTINE solve_lap_par

  SUBROUTINE solve_lap_thom(rhs,sol,aa,bb,cc,nmat)
    IMPLICIT NONE
    INTEGER                         :: nmat,i
    REAL(KIND=8), DIMENSION(nmat)   :: rhs, sol, bb, aa, cc

    REAL(KIND=8)      :: beta(nmat), gamm(nmat)
    INTEGER           :: k

    beta(1) = bb(1)
    gamm(1) = rhs(1)/bb(1)
    DO k = 2, nmat
       beta(k) = bb(k) - cc(k-1)*aa(k)/beta(k-1)
       gamm(k) = (rhs(k)-aa(k)*gamm(k-1))/beta(k)
    END DO

    sol(nmat) = gamm(nmat)
    DO k = 1, nmat-1
       sol(nmat-k) = gamm(nmat-k) - cc(nmat-k)*sol(nmat+1-k)/beta(nmat-k)
    END DO

    do i=1,nmat
!      write (*,*) 'solver tri',i,aa(i),bb(i),cc(i),rhs(i),sol(i)
    enddo


  END SUBROUTINE solve_lap_thom



END MODULE tripar
