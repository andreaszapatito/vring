MODULE trithom
  IMPLICIT NONE

CONTAINS

  SUBROUTINE solve_lap_thom(rhs,sol,aa,bb,cc,nmat)
    IMPLICIT NONE
    INTEGER                         :: nmat
    REAL(KIND=8), DIMENSION(nmat)   :: rhs, sol, bb, aa, cc

    REAL(KIND=8)      :: beta(nmat), gamm(nmat)
    INTEGER           :: k
    REAL(KIND=8)      :: t0, t1, tini, tfin

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



  END SUBROUTINE solve_lap_thom

END MODULE trithom
