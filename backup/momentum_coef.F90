! module momentum_coef

! calculation of local coefficients for the assembly of the linear 
! systems for the factorisation steps AF1 AF2 and AF3 for all the 
! three momentum variables q1 q2 q3



module momentum_coef

!#include <petsc/finclude/petscdef.h>

!  use petsc
  use run_tools
  use cons_tools
  use para_tools
  use comm_tools
  use mesh_tools
  use solu_tools


  implicit none
  contains
    subroutine ucoeff(var)
    implicit none
    type (run)   :: var

      integer jc,ji
      do jc=1,var%msh%nr

!   Coeff necessaires en HDNL1             
        var%su%h1d11(jc) =0.25d0*var%msh%dx1/(var%msh%rm(jc)*var%msh%axsym)
        var%su%h1d22d(jc)=0.25d0*var%msh%dx2/(var%msh%rm(jc)*var%msh%dmr(jc))
        var%su%h1d22n(jc)=0.125d0/(var%msh%rm(jc)**2)
        var%su%h1dare(jc)=var%msh%dx1/(var%par%Re*var%msh%rm(jc)**3)

!  Coeff necessaires en INVTR1            
        var%su%s1t11(jc)=var%msh%dx1q/((var%msh%rm(jc)**2)*var%msh%axsym)
        var%su%s1t2n(jc)=var%msh%rc(jc+1)**3/(var%msh%dcr(jc+1))
        var%su%s1t2s(jc)=var%msh%rc(jc)**3/(var%msh%dcr(jc))
        var%su%s1t2c(jc)=var%msh%dx2q/(var%msh%dmr(jc)*var%msh%rm(jc)**2)
!  Coeff necessaires en INVTR3


        var%su%s3t2n(jc)=var%msh%dx2q*var%msh%rc(jc+1)/(var%msh%rm(jc)*var%msh%dcr(jc+1)*var%msh%dmr(jc))
        var%su%s3t2s(jc)=var%msh%dx2q*var%msh%rc(jc)/(var%msh%rm(jc)*var%msh%dcr(jc)*var%msh%dmr(jc))


      enddo
      

      ji=1; if (var%com%ip_a(2).eq.0) ji=2
      do jc=ji,var%msh%nr
!  Coeff necessaires en HDNL2
        var%su%h2d11(jc) =0.25d0*var%msh%dx1/var%msh%rc(jc)*var%msh%axsym
        var%su%h2d22m(jc)=0.25d0*var%msh%dx2/(var%msh%rm(jc-1)*var%msh%dcr(jc))
        var%su%h2d22p(jc)=0.25d0*var%msh%dx2/(var%msh%rm(jc)*var%msh%dcr(jc))
        var%su%h2dare(jc)=var%msh%dx1/(var%par%Re*var%msh%rc(jc)*var%msh%axsym)
!  Coeff necessaires en INVTR2

        var%su%s2t11(jc)=var%msh%axsym*var%msh%dx1q/(var%msh%rc(jc)**2)
        var%su%s2t2p(jc)=var%msh%dx2q*var%msh%rc(jc)/var%msh%rm(jc)/(var%msh%dcr(jc)*var%msh%dmr(jc))
        var%su%s2t2m(jc)=var%msh%dx2q*var%msh%rc(jc)/var%msh%rm(jc-1)/(var%msh%dcr(jc)*var%msh%dmr(jc-1))
        var%su%s2t2c(jc)=var%su%s2t2p(jc)+var%su%s2t2m(jc)


      enddo      

    end subroutine ucoeff

    subroutine coesys(var)
    implicit none
    type (run)    :: var
    integer jc,kc
    real (kind=8) :: betadx,a22p,a22m
    real (kind=8) :: apj10,apj20,apj30
    real (kind=8) :: apj1n,apj2n,apj3n
    real (kind=8) :: apk10,apk20,apk30
    real (kind=8) :: apk1n,apk2n,apk3n
    real (kind=8) :: s1t2n0
    real (kind=8) :: s1t2c0
    real (kind=8) :: s2t2p0
    real (kind=8) :: s3t2n0
    integer                   :: ip,np
    integer                   :: tag
!    integer                   :: status
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
!    PetscErrorCode            :: ierr
    integer                   :: ierr



    betadx=0.5d0*var%par%crkalm(var%par%irkstep)*var%par%dt/var%par%Re


!   Inversion de l'eq de q1 sur J

!   Coeff du systeme
!   Pres de l'axe r=0 (jc=1) (dq_1/dr=0)


    if (var%com%ip_a(2).eq.0) then

      a22p=betadx*var%su%s1t2c(1)*var%su%s1t2n(1)    
      var%su%amj1(1)=0.d0
      var%su%acj1(1)=1.d0+a22p/var%msh%rm(1)
      var%su%apj1(1)=-a22p/var%msh%rm(2)
    else 
       a22p=betadx*var%su%s1t2c(1)*var%su%s1t2n(1)
       a22m=betadx*var%su%s1t2c(1)*var%su%s1t2s(1)

! apj1 for inerrmost cell for outer partitions
       s1t2n0=var%msh%rc(1)**3/(var%msh%dcr(1))
       s1t2c0=var%msh%dx2q/(var%msh%dmr(0)*var%msh%rm(0)**2)
       apj10=-betadx*s1t2c0*s1t2n0/var%msh%rm(1)

       var%su%apj1(1)= -a22p/var%msh%rm(2)
       var%su%acj1(1)=1.d0+(a22p+a22m)/var%msh%rm(1)
       var%su%amj1(1)= -a22m/var%msh%rm(0)
    endif

    if (var%com%ip_a(2).eq.var%com%np_a(2)-1)  then
!   A la frontiere exterieure (dq_1/dr=0)
      a22m=betadx*var%su%s1t2c(var%msh%nr)*var%su%s1t2s(var%msh%nr)
      var%su%apj1(var%msh%nr)=0.d0
      var%su%acj1(var%msh%nr)=1.d0+a22m/var%msh%rm(var%msh%nr)-betadx*var%su%s1t2c(var%msh%nr)*var%msh%rc(var%msh%nrc)*var%msh%rc(var%msh%nrc)*(var%msh%rc(var%msh%nrc)/var%msh%rm(var%msh%nr)*var%par%slip-var%msh%rc(var%msh%nrc)/var%msh%rm(var%msh%nrc))/var%msh%dcr(var%msh%nrc) ! maille fictive
      var%su%amj1(var%msh%nr)=-a22m/var%msh%rm(var%msh%nr-1)
    else
       a22p=betadx*var%su%s1t2c(var%msh%nr)*var%su%s1t2n(var%msh%nr)
       a22m=betadx*var%su%s1t2c(var%msh%nr)*var%su%s1t2s(var%msh%nr)
       var%su%apj1(var%msh%nr)= -a22p/var%msh%rm(var%msh%nr+1)
       var%su%acj1(var%msh%nr)=1.d0+(a22p+a22m)/var%msh%rm(var%msh%nr)
       var%su%amj1(var%msh%nr)= -a22m/var%msh%rm(var%msh%nr-1)
    endif




!   Dans le domaine calculable
      do jc=2,var%msh%nr-1
        a22p=betadx*var%su%s1t2c(jc)*var%su%s1t2n(jc)
        a22m=betadx*var%su%s1t2c(jc)*var%su%s1t2s(jc)
        var%su%apj1(jc)= -a22p/var%msh%rm(jc+1)
        var%su%acj1(jc)=1.d0+(a22p+a22m)/var%msh%rm(jc)
        var%su%amj1(jc)= -a22m/var%msh%rm(jc-1)
      enddo    

!   Factorisation L*U de la matrice

!      do ip=0,var%com%np_a(2)-1
!        if (var%com%ip_a(2).eq.ip) then
!
!          if (var%com%ip_a(2).eq.0) then
!            var%su%apj1(1)=var%su%apj1(1)/var%su%acj1(1)
!          else
!            var%su%acj1(1)=var%su%acj1(1)-var%su%amj1(1)*apj10
!            var%su%apj1(1)=var%su%apj1(1)/var%su%acj1(1)
!          endif
!
!          do jc=2,var%msh%nr
!            var%su%acj1(jc)=var%su%acj1(jc)-var%su%amj1(jc)*var%su%apj1(jc-1)
!            var%su%apj1(jc)=var%su%apj1(jc)/var%su%acj1(jc)
!          enddo
!          apj1n=var%su%apj1(var%msh%nr)
!        endif
!        call mpi_barrier(var%com%comm_a(2), ierr)
!
!
!        if (ip.lt.var%com%np_a(2)-1) then
!          tag = 10
!          if (var%com%ip_a(2).eq.ip  )                 call mpi_send(apj1n,1,MPI_REAL8,ip+1,tag,var%com%comm_a(2),ierr)
!          if (var%com%ip_a(2).eq.ip+1)                 call mpi_recv(apj10,1,MPI_DOUBLE_PRECISION,ip,tag,var%com%comm_a(2),status,ierr)
!        endif
!
!      enddo




!  Inversion de l'eq de q2 sur J
!  Coeff du systeme

    s2t2p0=var%msh%dx2q*var%msh%rc(0)/var%msh%rm(0)/(var%msh%dcr(0)*var%msh%dmr(0))
    apj20=-betadx*s2t2p0
    if (var%com%ip_a(2).eq.0) then
      var%su%acj2(1)=1.0d0
      var%su%amj2(1)=0.0d0
      var%su%apj2(1)=0.0d0
    else
      var%su%apj2(1)=-betadx*var%su%s2t2p(1)
      var%su%amj2(1)=-betadx*var%su%s2t2m(1)   
      var%su%acj2(1)=1.0d0+betadx*var%su%s2t2c(1)
    endif
    if (var%com%ip_a(2).eq.var%com%np_a(2)-1)  then
      var%su%acj2(var%msh%nrc)=1.0d0
      var%su%amj2(var%msh%nrc)=0.0d0
      var%su%apj2(var%msh%nrc)=0.0d0
    else
      var%su%apj2(var%msh%nrc)=-betadx*var%su%s2t2p(var%msh%nrc)
      var%su%amj2(var%msh%nrc)=-betadx*var%su%s2t2m(var%msh%nrc)   
      var%su%acj2(var%msh%nrc)=1.0d0+betadx*var%su%s2t2c(var%msh%nrc)
    endif

    do jc=1,var%msh%nr
      var%su%apj2(jc)=-betadx*var%su%s2t2p(jc)
      var%su%amj2(jc)=-betadx*var%su%s2t2m(jc)   
      var%su%acj2(jc)=1.0d0+betadx*var%su%s2t2c(jc)         
    enddo


!  Factorisation L*U de la matrice
!    do ip=0,var%com%np_a(2)-1
!      if (var%com%ip_a(2).eq.ip) then
!        if (var%com%ip_a(2).eq.0) var%su%apj2(1)=var%su%apj2(1)/var%su%acj2(1)
!        if (var%com%ip_a(2).ne.0) var%su%acj2(1)=var%su%acj2(1)-var%su%amj2(1)*apj20
!        if (var%com%ip_a(2).ne.0) var%su%apj2(1)=var%su%apj2(1)/var%su%acj2(1)
!      
!        if (var%com%ip_a(2).eq.0) then
!          do jc=2,var%msh%nrc
!            var%su%acj2(jc)=var%su%acj2(jc)-var%su%amj2(jc)*var%su%apj2(jc-1)
!            var%su%apj2(jc)=var%su%apj2(jc)/var%su%acj2(jc)
!          enddo
!        else 
!          do jc=2,var%msh%nrc
!            var%su%acj2(jc)=var%su%acj2(jc)-var%su%amj2(jc)*var%su%apj2(jc-1)
!            var%su%apj2(jc)=var%su%apj2(jc)/var%su%acj2(jc)
!          enddo
!        endif
!        apj2n=var%su%apj2(var%msh%nrc)
!      endif
!
!      call mpi_barrier(var%com%comm_a(2), ierr)
!
!
!      if (ip.lt.var%com%np_a(2)-1) then
!        tag = 10
!        if (var%com%ip_a(2).eq.ip  )                 call mpi_send(apj2n,1,MPI_DOUBLE_PRECISION,ip+1,tag,var%com%comm_a(2),ierr)
!        if (var%com%ip_a(2).eq.ip+1)                 call mpi_recv(apj20,1,MPI_DOUBLE_PRECISION,ip,tag,var%com%comm_a(2),status,ierr)
!      endif
!
!    enddo




!  Inversion de l'eq de q3 sur J
!  Coeff du systeme
      s3t2n0=var%msh%dx2q*var%msh%rc(1)/(var%msh%rm(0)*var%msh%dcr(1)*var%msh%dmr(0))

      apj30=-betadx*s3t2n0
      do jc=1,var%msh%nr
       var%su%apj3(jc)=-betadx*var%su%s3t2n(jc)
       var%su%amj3(jc)=-betadx*var%su%s3t2s(jc)*float(var%msh%jmc(jc))
       var%su%acj3(jc)=1.0d0-(var%su%apj3(jc)+var%su%amj3(jc))
      enddo

      if (var%com%ip_a(2).eq.var%com%np_a(2)-1) then 
        var%su%apj3(var%msh%nr)= 0.d0                      
        var%su%amj3(var%msh%nr)=-betadx*var%su%s3t2s(var%msh%nr)         
        var%su%acj3(var%msh%nr)=1.0d0-(betadx*var%su%s3t2n(var%msh%nr)*(var%par%slip-1.d0)+var%su%amj3(var%msh%nr))  ! maille fictive
      endif

!   Factorisation L*U de la matrice


!      do ip=0,var%com%np_a(2)-1
!        if (var%com%ip_a(2).eq.ip) then
!          if (var%com%ip_a(2).eq.0)  var%su%apj3(1)=var%su%apj3(1)/var%su%acj3(1)
!          if (var%com%ip_a(2).ne.0)  var%su%acj3(1)=var%su%acj3(1)-var%su%amj3(1)*apj30
!          if (var%com%ip_a(2).ne.0)  var%su%apj3(1)=var%su%apj3(1)/var%su%acj3(1)
!
!          do jc=2,var%msh%nr
!            var%su%acj3(jc)=var%su%acj3(jc)-var%su%amj3(jc)*var%su%apj3(jc-1)
!            var%su%apj3(jc)=var%su%apj3(jc)/var%su%acj3(jc)
!          enddo
!          apj3n=var%su%apj3(var%msh%nr)
!        endif
!
!        call mpi_barrier(var%com%comm_a(2), ierr)
!
!
!        if (ip.lt.var%com%np_a(2)-1) then
!          tag = 10
!          if (var%com%ip_a(2).eq.ip  )                 call mpi_send(apj3n,1,MPI_DOUBLE_PRECISION,ip+1,tag,var%com%comm_a(2),ierr)
!          if (var%com%ip_a(2).eq.ip+1)                 call mpi_recv(apj30,1,MPI_DOUBLE_PRECISION,ip,tag,var%com%comm_a(2),status,ierr)
!        endif
!
!      enddo




      do jc=1,var%msh%nr
!        write (500+var%com%ip,"(100(E32.16))") var%msh%rc(jc), var%su%amj1(jc),var%su%acj1(jc),var%su%apj1(jc),var%su%amj2(jc),var%su%acj2(jc),var%su%apj2(jc),var%su%amj3(jc),var%su%acj3(jc),var%su%apj3(jc)
      enddo



!   Inversion du scalaire   sur J
      betadx=betadx*var%par%Re/var%par%Pe
!   Coeff du systeme
      do jc=1,var%msh%nr
       var%su%apjs(jc)=-betadx*var%su%s3t2n(jc)*float(var%msh%jpc(jc))
       var%su%amjs(jc)=-betadx*var%su%s3t2s(jc)*float(var%msh%jmc(jc))
       var%su%acjs(jc)=1.0d0-(var%su%apjs(jc)+var%su%amjs(jc))
      enddo
!   Factorisation L*U de la matrice
      var%su%apjs(1)=var%su%apjs(1)/var%su%acjs(1)
      do jc=2,var%msh%nr
        var%su%acjs(jc)=var%su%acjs(jc)-var%su%amjs(jc)*var%su%apjs(jc-1)
        var%su%apjs(jc)=var%su%apjs(jc)/var%su%acjs(jc)
      enddo

      betadx=0.5d0*var%msh%dx3q*var%par%crkalm(var%par%irkstep)*var%par%dt/var%par%Re

!   Dans le domaine calculable
      do kc=1,var%msh%nz
        var%su%amk1(kc)=-betadx
        var%su%ack1(kc)=1.d0+2.d0*betadx
        var%su%apk1(kc)=-betadx
      enddo

!    A la frontiere nord
     if (var%com%ip_a(3).eq.(var%com%np_a(3)-1)) then
         if (var%msh%if_grad_null.eq.1) then
           var%su%amk1(var%msh%nzc)=1.0d0
           var%su%ack1(var%msh%nzc)=-1.0d0
           var%su%apk1(var%msh%nzc)=0.0d0
         else
           var%su%amk1(var%msh%nzc)=1.0d0
           var%su%ack1(var%msh%nzc)=1.0d0
           var%su%apk1(var%msh%nzc)=0.0d0
         endif
      endif
!    A la frontiere sud
      if (var%com%ip_a(3).eq.0) then
      var%su%amk1(1)=-8.d0/3.d0*betadx !! quand on impose Vtheta (voir coesys.f)
      var%su%ack1(1)= 1.d0+4.d0*betadx
      var%su%apk1(1)=-4.d0/3.d0*betadx
      endif
!    FACTORISATION L*U

!      do ip=0,var%com%np_a(3)-1
!        if (var%com%ip_a(3).eq.ip) then
!         if (var%com%ip_a(3).eq.0) var%su%apk1(1)=var%su%apk1(1)/var%su%ack1(1)
!         if (var%com%ip_a(3).ne.0) var%su%ack1(1)=var%su%ack1(1)-var%su%amk1(1)*apk10
!         if (var%com%ip_a(3).ne.0) var%su%apk1(1)=var%su%apk1(1)/var%su%ack1(1)
!
!         do kc=2,var%msh%nzc
!           var%su%ack1(kc)=var%su%ack1(kc)-var%su%amk1(kc)*var%su%apk1(kc-1)
!           var%su%apk1(kc)=var%su%apk1(kc)/var%su%ack1(kc)
!         enddo
!         apk1n=var%su%apk1(var%msh%nzc)
!        endif
!
!        call mpi_barrier(var%com%comm_a(3), ierr)
!
!
!        if (ip.lt.var%com%np_a(3)-1) then
!          tag = 10
!          if (var%com%ip_a(3).eq.ip  )                 call mpi_send(apk1n,1,MPI_DOUBLE_PRECISION,ip+1,tag,var%com%comm_a(3),ierr)
!          if (var%com%ip_a(3).eq.ip+1)                 call mpi_recv(apk10,1,MPI_DOUBLE_PRECISION,ip,tag,var%com%comm_a(3),status,ierr)
!        endif
!
!      enddo










!    Inversion de l'eq de q2 sur K
!    Dans le domaine calculable

         do kc=1,var%msh%nz
           var%su%amk2(kc)=-betadx
           var%su%ack2(kc)=1.d0+2.d0*betadx
           var%su%apk2(kc)=-betadx
         enddo

!     A la frontiere nord
       if (var%com%ip_a(3).eq.(var%com%np_a(3)-1)) then

         if (var%msh%if_grad_null.eq.1) then
           var%su%amk2(var%msh%nzc)=1.0d0
           var%su%ack2(var%msh%nzc)=-1.0d0
           var%su%apk2(var%msh%nzc)=0.0d0
         else
           var%su%amk2(var%msh%nzc)=1.0d0
           var%su%ack2(var%msh%nzc)=1.0d0
           var%su%apk2(var%msh%nzc)=0.0d0
         endif
       endif
!     A la frontiere sud

       if (var%com%ip_a(3).eq.0) then
         var%su%amk2(1)= -8.d0/3.d0*betadx !! quand on impose q2s (voir step2.f) 
         var%su%ack2(1)=1.d0+4.d0*betadx
         var%su%apk2(1)=-4.d0/3.d0*betadx
       endif
!     FACTORISATION L*U

!      do ip=0,var%com%np_a(3)-1
!        if (var%com%ip_a(3).eq.ip) then
!         if (var%com%ip_a(3).eq.0) var%su%apk2(1)=var%su%apk2(1)/var%su%ack2(1)       
!         if (var%com%ip_a(3).ne.0) var%su%ack2(1)=var%su%ack2(1)-var%su%amk2(1)*apk20
!         if (var%com%ip_a(3).ne.0) var%su%apk2(1)=var%su%apk2(1)/var%su%ack2(1)
!
!         do kc=2,var%msh%nzc
!           var%su%ack2(kc)=var%su%ack2(kc)-var%su%amk2(kc)*var%su%apk2(kc-1)
!           var%su%apk2(kc)=var%su%apk2(kc)/var%su%ack2(kc)
!         enddo
!         apk2n=var%su%apk2(var%msh%nzc)
!        endif
!
!        call mpi_barrier(var%com%comm_a(3), ierr)
!
!
!        if (ip.lt.var%com%np_a(3)-1) then
!          tag = 10
!          if (var%com%ip_a(3).eq.ip  )                 call mpi_send(apk2n,1,MPI_DOUBLE_PRECISION,ip+1,tag,var%com%comm_a(3),ierr)
!          if (var%com%ip_a(3).eq.ip+1)                 call mpi_recv(apk20,1,MPI_DOUBLE_PRECISION,ip,tag,var%com%comm_a(3),status,ierr)
!        endif
!
!      enddo








!     Inversion de l'eq de q3 sur K
!     Dans le domaine calculable

         do kc=1,var%msh%nz
           var%su%amk3(kc)=-betadx
           var%su%ack3(kc)=1.d0+2.d0*betadx
           var%su%apk3(kc)=-betadx
         enddo


!     A la fronaiere nord
         if (var%com%ip_a(3).eq.(var%com%np_a(3)-1)) then
         if (var%msh%if_grad_null.eq.1) then
           var%su%amk3(var%msh%nzc)=1.0d0
           var%su%ack3(var%msh%nzc)=-1.0d0
           var%su%apk3(var%msh%nzc)=0.0d0
         else 
           var%su%amk3(var%msh%nzc)=0.0d0
           var%su%ack3(var%msh%nzc)=1.0d0
           var%su%apk3(var%msh%nzc)=0.0d0
         endif
         endif
!     A la frontiere sud
         if (var%com%ip_a(3).eq.0) then
         var%su%amk3(1)=0.0d0
         var%su%ack3(1)=1.0d0
         var%su%apk3(1)=0.0d0
         endif
!     FACTORISATION L*U


!      do ip=0,var%com%np_a(3)-1
!        if (var%com%ip_a(3).eq.ip) then
!         if (var%com%ip_a(3).eq.0) var%su%apk3(1)=var%su%apk3(1)/var%su%ack3(1)       
!         if (var%com%ip_a(3).ne.0) var%su%ack3(1)=var%su%ack3(1)-var%su%amk3(1)*apk30
!         if (var%com%ip_a(3).ne.0) var%su%apk3(1)=var%su%apk3(1)/var%su%ack3(1)
!
!         do kc=2,var%msh%nzc
!           var%su%ack3(kc)=var%su%ack3(kc)-var%su%amk3(kc)*var%su%apk3(kc-1)
!           var%su%apk3(kc)=var%su%apk3(kc)/var%su%ack3(kc)
!         enddo
!         apk3n=var%su%apk3(var%msh%nzc)
!        endif
!
!        call mpi_barrier(var%com%comm_a(3), ierr)
!
!
!        if (ip.lt.var%com%np_a(3)-1) then
!          tag = 10
!          if (var%com%ip_a(3).eq.ip  )                 call mpi_send(apk3n,1,MPI_DOUBLE_PRECISION,ip+1,tag,var%com%comm_a(3),ierr)
!          if (var%com%ip_a(3).eq.ip+1)                 call mpi_recv(apk30,1,MPI_DOUBLE_PRECISION,ip,tag,var%com%comm_a(3),status,ierr)
!        endif
!
!      enddo


          


      do kc=1,var%msh%nz
!        write (600+var%com%ip,"(100(E32.16))") var%msh%zc(kc), var%su%amk1(kc),var%su%ack1(kc),var%su%apk1(kc),var%su%amk2(kc),var%su%ack2(kc),var%su%apk2(kc),var%su%amk3(kc),var%su%ack3(kc),var%su%apk3(kc)
      enddo

!     Inversion du scalaire   sur K           c
      betadx=0.5d0*var%msh%dx3q*var%par%crkalm(var%par%irkstep)*var%par%dt/var%par%Re
!     Dans le domaine calculable
      do  kc=2,var%msh%nz
        var%su%amks(kc)=-betadx
        var%su%acks(kc)=1.d0+2.d0*betadx
        var%su%apks(kc)=-betadx
      enddo
!     A la frontiere nord

         if (var%msh%if_grad_null.eq.1) then
           var%su%amks(var%msh%nzc)=1.0d0
           var%su%acks(var%msh%nzc)=-1.0d0
           var%su%apks(var%msh%nzc)=0.0d0
         else
           var%su%amks(var%msh%nzc)=1.0d0
           var%su%acks(var%msh%nzc)=1.0d0
           var%su%apks(var%msh%nzc)=0.0d0
         endif
           var%su%amks(var%msh%nzc)=0.0d0
           var%su%acks(var%msh%nzc)=1.0d0
           var%su%apks(var%msh%nzc)=0.0d0
!     A la frontiere sud

         var%su%amks(1)=-8.d0/3.d0*betadx
         var%su%acks(1)=1.d0+4.d0*betadx
         var%su%apks(1)=-4.d0/3.d0*betadx

!     FACTORISATION L*U
        var%su%apks(1)=var%su%apks(1)/var%su%acks(1)       
       do  kc=2,var%msh%nzc
        var%su%acks(kc)=var%su%acks(kc)-var%su%amks(kc)*var%su%apks(kc-1)
        var%su%apks(kc)=var%su%apks(kc)/var%su%acks(kc)
       enddo



    end subroutine coesys




end module momentum_coef
