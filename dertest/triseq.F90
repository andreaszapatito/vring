MODULE triseq
!#include <petsc/finclude/petscdef.h>
  IMPLICIT NONE
#include "mpif.h"
CONTAINS

  SUBROUTINE solve_lap_tri(rhs,sol,aa,bb,cc,nmat)
    IMPLICIT NONE
    INTEGER                         :: nmat
    REAL(KIND=8), DIMENSION(nmat)   :: rhs, sol, bb, aa, cc
    

    REAL(KIND=8), SAVE, ALLOCATABLE  :: as(:), bs(:), cs(:)
    REAL(KIND=8), SAVE, ALLOCATABLE  :: yy(:)
    INTEGER             :: k
    REAL(KIND=8)        :: t0, t1, t2, tini, tfin
    LOGICAL, SAVE       :: once=.TRUE.

    tini = MPI_Wtime()
    t0 = MPI_Wtime()
    IF (.NOT.once) THEN
       IF (SIZE(yy)/=nmat) THEN
          DEALLOCATE(yy,bs,as,cs)
          once = .TRUE.
       END IF
    END IF
    
    IF (once) THEN
       once=.FALSE.
       ALLOCATE(yy(nmat), bs(nmat), as(nmat), cs(nmat))
       yy = 0.d0
       as = aa
       bs(1) = bb(1)
       cs(1) = cc(1)/bb(1)
       DO k = 2, nmat-1
          bs(k) = bb(k) - cs(k-1)*aa(k)
          cs(k) = cc(k)/bs(k)
       END DO
       bs(nmat) = bb(nmat)-cs(nmat-1)*as(nmat)
    END IF

    t1 = MPI_Wtime()
!    PRINT*, 'solve_lap_tri, factorisation', t1-t0

    yy(1) = rhs(1)/bs(1)
    DO k = 2, nmat
       yy(k) = (rhs(k)-as(k)*yy(k-1))/bs(k)
    END DO

    sol(nmat) = yy(nmat)
    DO k = 1, nmat-1
       sol(nmat-k) = yy(nmat-k) - cs(nmat-k)*sol(nmat-k+1)
    END DO
    t0 = MPI_Wtime()
!    PRINT*, 'solve_lap_tri, resolution   ', t0-t1
    tfin = MPI_Wtime()
!    PRINT*, '===TRI===', tfin-tini



  END SUBROUTINE solve_lap_tri

END MODULE triseq
