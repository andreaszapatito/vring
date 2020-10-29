MODULE btrithom
  IMPLICIT NONE
CONTAINS

  SUBROUTINE solve_lap_bthom(rhs,sol,b1,b2,b3,m,n)

!   m block size
!   n block number
!   b1,b2,b3 diagonal block bands

    IMPLICIT NONE
    INTEGER                        :: m,n
    REAL(KIND=8), DIMENSION(n*m)   :: rhs, sol, b1, b2, b3

    REAL(KIND=8)      :: delta(m,m,n),beta(m,m)
    REAL(KIND=8)      :: gamm(m,n)
    INTEGER           :: k,l
    REAL(KIND=8)      :: t0, t1, tini, tfin

    
    delta=0.0
    beta=0.0
    do k = 1,m
      beta(k,k)=b2(k)
    end do
    do k = 2,m
      beta(k,k-1)=b1(k)
    end do
    do k = 1,m-1
      beta(k,k+1)=b3(k)
    end do

    gamm(1:m,1) = rhs(1:m)
    call inverse(beta,delta(:,:,1),m)

    do k =2,n
      beta=0.0
      do l =1,m
        beta(l,l)   = b2((k-1)*m+l)
      end do
      do l =2,m
        beta(l,l-1) = b1((k-1)*m+l)
      end do
      do l =1,m-1
        beta(l,l+1) = b3((k-1)*m+l)
      end do

      beta(:,:) = beta(:,:)-1.0*delta(:,:,k-1)
!      gamm(1:m,k) = rhs((k-1)*m+1:k*m) - gamm(1:m,k-1)
      gamm(1:m,k) = rhs((k-1)*m+1:k*m) -matmul(delta(:,:,k-1),gamm(:,k-1))

      call inverse(beta,delta(:,:,k),m)

      do l =1,m
!        write (*,*) "gamm :",gamm(l,k),k,rhs((k-1)*m+l)
      end do
    end do

    sol((n-1)*m+1:n*m) = -1.0*matmul(delta(:,:,n),-(gamm(:,n)))
!    sol((n-1)*m+1:n*m) =                           gamm(:,n)  

    do k =(n-1)*m+1,n*m
    !  write (*,*) "solution :",sol(k)
    end do

    do l = 1,n-1
      k=n-l
      sol((k-1)*m+1:k*m) = -1.0*matmul(delta(:,:,k),sol(k*m+1:(k+1)*m)-gamm(:,k))
    end do

  END SUBROUTINE solve_lap_bthom

!!$  SUBROUTINE solve_lap_thom(rhs,sol,aa,bb,cc,nmat)
!!$    IMPLICIT NONE
!!$    INTEGER                         :: nmat
!!$    REAL(KIND=8), DIMENSION(nmat)   :: rhs, sol, bb, aa, cc
!!$
!!$    REAL(KIND=8)      :: beta(nmat), gamm(nmat)
!!$    INTEGER           ::  k
!!$    REAL(KIND=8)      :: t0, t1, tini, tfin
!!$
!!$    beta(1) = bb(1)
!!$    gamm(1) = rhs(1)/bb(1)
!!$    DO k = 2, nmat
!!$       beta(k) = bb(k) - cc(k-1)*aa(k)/beta(k-1)
!!$       gamm(k) = (rhs(k)-aa(k)*gamm(k-1))/beta(k)
!!$    END DO
!!$
!!$    sol(nmat) = gamm(nmat)
!!$    DO k = 1, nmat-1
!!$       sol(nmat-k) = gamm(nmat-k) - cc(nmat-k)*sol(nmat+1-k)/beta(nmat-k)
!!$    END DO
!!$
!!$
!!$
!!$  END SUBROUTINE solve_lap_thom


  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
 REAL(KIND=8)  ::a(n,n), c(n,n)
 REAL(KIND=8)  ::L(n,n), U(n,n), b(n), d(n), x(n)
 REAL(KIND=8)  ::coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

END MODULE btrithom
