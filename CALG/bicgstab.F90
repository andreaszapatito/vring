 function BiCGStab(A,b) result(x)
        implicit none

!--------------------------PARAMETER AND VARIABLE-------------------------------!
        real    (kind=8), intent(in )                   :: A (:,:)
        real    (kind=8), intent(in )                   :: b ( : )
        real    (kind=8), dimension(1:size(b, dim=1))   :: x

        real    (kind=8), dimension(1:size(b, dim=1))   :: r, rs, v, p, s, t
        real    (kind=8), parameter                     :: e = 1d-33
        real    (kind=8)                                :: rho      , rho_prev
        real    (kind=8)                                :: alpha    , omega   , beta
        real    (kind=8)                                :: norm_r   , norm_b       
        real    (kind=8)                                :: summesion, temp

        integer                                         :: it=0,err
!------------------------END PARAMETER AND VARIABLE-----------------------------!  

        if(size(A, dim=1) /= size(A, dim=2)) stop &
        "Error: Improper dimension of matrix A in BiCGStab."





!-------------------------------------------------------!
        x  = 0.0d0                                     !-------> INITIAL GUESS
!-------------------------------------------------------!
        r  = b - matmul(A,x)                            !-------> LINE 1
        rs = r                                          !
!-------------------------------------------------------!
        rho   = 1.0d0; alpha = 1.0d0; omega = 1.0d0  !-------> LINE 2
!-------------------------------------------------------!
        v  = 0.0d0; p  = 0.0d0                        !-------> LINE 3
!                                                       !
        norm_r = sqrt(dot_product(r,r))                 !
        norm_b = sqrt(dot_product(b,b))                 !
!-------------------------------------------------------!





        do while(norm_r .GT. e*norm_b)                          !-------> START OF LOOP

        !-------------------------------------------------------!
            rho_prev = rho                                      !-------> LINE 5
            rho      = dot_product(rs,r)                        !
        !-------------------------------------------------------!
            beta     = (rho/rho_prev) * (alpha/omega)           !-------> LINE 6
        !-------------------------------------------------------!
            p        = r + beta * (p - omega*v)                 !-------> LINE 7
        !-------------------------------------------------------!
            v        = matmul(A,p)                              !-------> LINE 8
        !-------------------------------------------------------!
            alpha    = rho/dot_product(rs,v)                    !-------> LINE 9
        !-------------------------------------------------------!
            s        = r - alpha*v                              !-------> LINE 10
        !-------------------------------------------------------!
            t        = matmul(A,s)                              !-------> LINE 11
        !-------------------------------------------------------!
            omega    = dot_product(t,s)/dot_product(t,t)        !-------> LINE 12
        !-------------------------------------------------------!
            x        = x + alpha*p + omega*s                    !-------> LINE 13
        !-------------------------------------------------------!
            r        = s - omega*t                              !-------> LINE 17
        !-------------------------------------------------------!
            norm_r   = sqrt(dot_product(r,r))                   !
            norm_b   = sqrt(dot_product(b,b))                   !
        !-------------------------------------------------------!
            it = it + 1                                         !
        !-------------------------------------------------------!

        end do                                                      !-------> END OF LOOP

        print*, "Iteration Required :", it

return
end function BiCGStab     

