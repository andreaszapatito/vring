! module convective

module momentum_terms

  use run_tools
  use cons_tools
  use para_tools
  use comm_tools
  use mesh_tools
  use solu_tools


  implicit none
  contains
  subroutine hdnl1(var)

!   Calcul des termes explicites pour
!   l'equation de q_1
!   point de calcul (i,j+1/2,k+1/2)

    implicit none
    type (run)   :: var

    integer ic,kc,jc,ip,im,ik

    REAL(KIND=8)  ::qnn,qss
    REAL(KIND=8)  ::h11,h12d,h12n,h13,d1q2

    do ic=1,var%msh%ntheta
      ip=var%msh%ipv(ic)
      im=var%msh%imv(ic)
      do kc=1,var%msh%nz
        do jc=1,var%msh%nr
!  Terme 1/r*d/dth(q1*q1) a i,j+1/2,k+1/2 
          qnn=var%su%q1(ip,jc,kc)+var%su%q1(ic,jc,kc)
          qss=var%su%q1(im,jc,kc)+var%su%q1(ic,jc,kc)
          h11=(qnn*qnn-qss*qss)*var%su%h1d11(jc)

!  Terme 1/r*d/dr(q1*q2) a i,j+1/2,k+1/2

          qnn=(var%su%q2(ic,jc+1,kc)+var%su%q2(im,jc+1,kc))*(var%su%q1(ic,jc+1,kc)+var%su%q1(ic,jc,kc))
          qss=(var%su%q1(ic,jc,kc)+var%su%q1(ic,var%msh%jmv(jc),kc))*(var%su%q2(ic,jc,kc)+var%su%q2(im,jc,kc))
          h12d=(qnn-qss)*var%su%h1d22d(jc)
!  Terme 1/r^2*(q1*q2) a i,j+1/2,k+1/2

          qnn=(var%su%q1(ip,jc,kc)+var%su%q1(ic,jc,kc))*(var%su%q2(ic,jc+1,kc)+var%su%q2(ic,jc,kc))
          qss=(var%su%q1(ic,jc,kc)+var%su%q1(im,jc,kc))*(var%su%q2(im,jc+1,kc)+var%su%q2(im,jc,kc))
          h12n=(qnn+qss)*var%su%h1d22n(jc)

!  Terme d/dz(q1*q3) a i,j+1/2,k+1/2 

          qnn=0.25d0*(var%su%q3(ic,jc,kc+1)+var%su%q3(im,jc,kc+1))*(var%su%q1(ic,jc,kc+1)+var%su%q1(ic,jc,kc))

          if(var%com%ip_a(3).eq.0.and.kc.eq.1) then
            qss=0.5d0*(var%su%q3s(ic,jc)+var%su%q3s(im,jc))*var%su%q1s(ic,jc)
          else
            qss=0.25d0*(var%su%q3(ic,jc,kc)+var%su%q3(im,jc,kc))*(var%su%q1(ic,jc,kc)+var%su%q1(ic,jc,kc-1))
          endif
          h13=(qnn-qss)*var%msh%dx3
!  Terme total convectif C1
!  Terme visqueux 2/Re*1/r^3*d/dth(q2)
 
          qnn=(var%su%q2(ic,jc+1,kc)+var%su%q2(ic,jc,kc))
          qss=(var%su%q2(im,jc+1,kc)+var%su%q2(im,jc,kc))
          d1q2=(qnn-qss)*var%su%h1dare(jc)
!  Terme explicite total H1
 
          var%su%dq1(ic,jc,kc)=d1q2-(h11+(h12d+h12n)+h13)



       enddo
     enddo
   enddo


  end subroutine hdnl1
  subroutine hdnl2(var)

!   Calcul des termes explicites pour
!   l'equation de q_1
!   point de calcul (i+1/2,j,k+1/2)     c

    implicit none
    type (run)   :: var
    integer ic,kc,jc,ip,im,ij
     REAL(KIND=8)  ::qnn,qss,h21,h22,h23,h11n,d1q1

!   ik=1; if(var%com%ip_a(3).eq.var%com%np_a(3)-1) ik=2
    ij=1; if(var%com%ip_a(2).eq.0) ij=2



    do ic=1,var%msh%ntheta
      ip=var%msh%ipv(ic)
      im=var%msh%imv(ic)
      do kc=1,var%msh%nz
        do jc=ij,var%msh%nr
!  Terme 1/r*d/dth(q1*q2) a i+1/2,j,k+1/2
          qnn=(var%su%q1(ip,jc,kc)+var%su%q1(ip,jc-1,kc))*(var%su%q2(ip,jc,kc)+var%su%q2(ic,jc,kc))
          qss=(var%su%q1(ic,jc,kc)+var%su%q1(ic,jc-1,kc))*(var%su%q2(ic,jc,kc)+var%su%q2(im,jc,kc))
          h21=(qnn-qss)*var%su%h2d11(jc)

!  Terme d/dr(q2^2/r) a i+1/2,j,k+1/2
          qnn=var%su%q2(ic,jc+1,kc)+var%su%q2(ic,jc,kc)
          qss=var%su%q2(ic,jc,kc)+var%su%q2(ic,jc-1,kc)
          h22=qnn*qnn*var%su%h2d22p(jc)-qss*qss*var%su%h2d22m(jc)

!  Terme q1*q1 a i+1/2,j,k+1/2   
          qnn=0.5d0*(var%su%q1(ic,jc,kc)+var%su%q1(ip,jc,kc))
          qss=0.5d0*(var%su%q1(ic,jc-1,kc)+var%su%q1(ip,jc-1,kc))
          h11n=0.25d0*(qnn+qss)*(qnn+qss)

!  Terme d/dz(q2*q3) a i+1/2,j,k+1/2  
          qnn=0.25d0*(var%su%q3(ic,jc,kc+1)+var%su%q3(ic,jc-1,kc+1))*(var%su%q2(ic,jc,kc+1)+var%su%q2(ic,jc,kc))
          if(var%com%ip_a(3).eq.0.and.kc.eq.1) then
            qss=0.5d0*(var%su%q3s(ic,jc)+var%su%q3s(ic,jc-1))*var%su%q2s(ic,jc)
          else
            qss=0.25d0*(var%su%q3(ic,jc,kc)+var%su%q3(ic,jc-1,kc))*(var%su%q2(ic,jc,kc)+var%su%q2(ic,jc,kc-1))
          endif
          h23=(qnn-qss)*var%msh%dx3
!  Terme total convectif  C2          
!  Terme diffusif -2/(r*Re)*d/dth(q1)
          qnn=(var%su%q1(ip,jc,kc)+var%su%q1(ip,jc-1,kc))
          qss=(var%su%q1(ic,jc,kc)+var%su%q1(ic,jc-1,kc))
          d1q1=-(qnn-qss)*var%su%h2dare(jc)
!  Terme explicite total a i+1/2,j,k+1/2
          var%su%dq2(ic,jc,kc)=d1q1-(h21+h22+h23-h11n)



       enddo
     enddo
   enddo
  end subroutine hdnl2

  subroutine hdnl3(var)

!   Calcul des termes explicites pour
!   l'equation de q_1
!   point de calcul (i+1/2,j+1/2,k)

    implicit none
    type (run)   :: var
    integer ic,kc,jc,ip,im,ik
    REAL(KIND=8)  ::qnn,qss,h31,h32,h33

!   ik=1; if(var%com%ip_a(3).eq.var%com%np_a(3)-1) ik=2
    ik=1; if(var%com%ip_a(3).eq.0) ik=2
    do ic=1,var%msh%ntheta
      ip=var%msh%ipv(ic)
      im=var%msh%imv(ic)
      do kc=ik,var%msh%nz
        do jc=1,var%msh%nr

!   Terme 1/r*d/dr(q2*q3) a i+1/2,j+1/2,k
          qnn=(var%su%q3(ic,jc+1,kc)+var%su%q3(ic,jc,kc))*(var%su%q2(ic,jc+1,kc)+var%su%q2(ic,jc+1,kc-1))

          qss=(var%su%q3(ic,var%msh%jmv(jc),kc)+var%su%q3(ic,jc,kc))*(var%su%q2(ic,jc,kc)+var%su%q2(ic,jc,kc-1))
          h32=(qnn-qss)*var%su%h1d22d(jc)
!   Terme 1/r*d/dth(q3*q1) a i+1/2,j+1/2,k
          qnn=(var%su%q3(ic,jc,kc)+var%su%q3(ip,jc,kc))*(var%su%q1(ip,jc,kc)+var%su%q1(ip,jc,kc-1))
          qss=(var%su%q3(im,jc,kc)+var%su%q3(ic,jc,kc))*(var%su%q1(ic,jc,kc)+var%su%q1(ic,jc,kc-1))
          h31=(qnn-qss)*var%su%h1d11(jc)
!   Terme d/dz(q3*q3) a i+1/2,j+1/2,k
          qnn=0.25d0*(var%su%q3(ic,jc,kc)+var%su%q3(ic,jc,kc+1))*(var%su%q3(ic,jc,kc)+var%su%q3(ic,jc,kc+1))
          qss=0.25d0*(var%su%q3(ic,jc,kc)+var%su%q3(ic,jc,kc-1))*(var%su%q3(ic,jc,kc)+var%su%q3(ic,jc,kc-1))
          h33=(qnn-qss)*var%msh%dx3
!   Terme total convectif C3 a i+1/2,j+1/2,k
          var%su%dq3(ic,jc,kc)=-(h31+h32+h33)
       enddo
     enddo
   enddo
  end subroutine hdnl3


  subroutine hdnbc(var)

    implicit none
    type (run)   :: var
    integer ic,kc,jc,ip,im,ierr
    REAL(KIND=8)  ::qnn,qss,h31,h32,h33
    REAL(KIND=8)  ::red,qout,qinf,cor,dt,t1,t0,ts,r,th
    REAL(KIND=8)  ::redall,qoutall,qinfall


    dt=var%par%dt
    ts=0.5
    t0=var%par%ntime+var%par%crkrom(var%par%irkstep)*dt
    t1=var%par%ntime+(var%par%crkgam(var%par%irkstep)+var%par%crkrom(var%par%irkstep))*dt

!   no inlet
!    t1=t0

    r=0.1
    th=1.0
    var%par%unom0 = profile(t0,r,th)
    var%par%unom1 = profile(t1,r,th) 
   

!    var%su%corr3n(:,:)= 0.d0
    var%su%corr3s(:,:)= 0.d0
    do ic=1,var%msh%ntheta
     do jc=1,var%msh%nr
       r=var%msh%rc(jc)
       th=var%msh%thc(ic)
!      if (var%msh%rc(jc).le.0.1999) var%su%corr3s(ic,jc)= (var%par%unom1-var%par%unom0)*(var%par%crkgam(var%par%irkstep)+var%par%crkrom(var%par%irkstep))
      if (var%msh%rc(jc).le.0.4999) var%su%corr1s(ic,jc)= (profilet(t1,r,th)-profilet(t0,r,th))
      if (var%msh%rc(jc).le.0.4999) var%su%corr2s(ic,jc)= (profiler(t1,r,th)-profiler(t0,r,th))
      if (var%msh%rc(jc).le.0.4999) var%su%corr3s(ic,jc)= (var%par%unom1-var%par%unom0)
     enddo
    enddo
 
    do ic=1,var%msh%ntheta
     do jc=1,var%msh%nr
      var%su%q1s(ic,jc)=0.0
      var%su%q2s(ic,jc)=0.0*var%msh%rc(jc)
      if (var%com%ip_a(3).eq.0) var%su%q3s(ic,jc)=0.0
      if (var%com%ip_a(3).ne.0) var%su%q3s(ic,jc)=0.0
     enddo
    enddo

    do  ic=1,var%msh%ntheta
     do  jc=1,var%msh%nr
!       var%su%csort(ic,jc)=0.0
      if (var%su%q3(ic,jc,var%msh%nzc).gt.0.0) var%su%csort(ic,jc)=0.6*var%su%q3(ic,jc,var%msh%nzc)
!     var%su%csort(ic,jc)=0.0
     enddo
    enddo


    do ic=1,var%msh%ntheta
     do jc=1,var%msh%nrc
       var%su%tconv(ic,jc)=-var%su%csort(ic,jc)*var%msh%dx3*(var%su%q3(ic,jc,var%msh%nzc)-var%su%q3(ic,jc,var%msh%nz))
       var%su%corr3n(ic,jc)=(var%par%crkgam(var%par%irkstep)*var%su%tconv(ic,jc)+var%par%crkrom(var%par%irkstep)*var%su%dq3no(ic,jc))*var%par%dt
       var%su%dq3no(ic,jc)=var%su%tconv(ic,jc)
!       var%su%corr3n(ic,jc)=0.d0
     enddo
    enddo
!
! FRONTIERE SUD  --> voir fluxin pour corr3s(ic,jc)
! Correction pour conserver la masse    
! Calcul de la loi de redistribution
!
      do jc=1,var%msh%nr
       var%msh%da(jc)=var%msh%rm(jc)*var%msh%dmr(jc)/(var%msh%dx1*var%msh%dx2)
!!       var%msh%qc(jc)=1.d0-((var%msh%rm(jc)-var%msh%rc(1))/(var%msh%rc(var%msh%nrc)-var%msh%rc(1)))**2
       var%msh%qc(jc)=1.d0-(var%msh%rm(jc)/var%msh%rmaxg)**2
      enddo

      red=0.d0
      do jc=1,var%msh%ntheta
       do  ic=1,var%msh%nr
       enddo
      enddo
! Calcul du flux de sortie
      qout=0.d0
      qinf=0.d0
      do jc=1,var%msh%nr
       do ic=1,var%msh%ntheta
        red=red+var%msh%qc(jc)*var%msh%qc(jc)*var%msh%da(jc)
        if(var%com%ip_a(3).eq.0)                 qinf=qinf+var%su%corr3s(ic,jc)*var%msh%da(jc)
        if(var%com%ip_a(3).eq.var%com%np_a(3)-1) qout=qout+var%su%corr3n(ic,jc)*var%msh%da(jc)
       enddo
      enddo

      call mpi_allreduce(red,  redall,1,mpi_double_precision,mpi_sum,var%com%comm_0,ierr)

      red =redall/float(var%com%np_a(3))

      call mpi_allreduce(qout,qoutall,1,mpi_double_precision,mpi_sum,var%com%comm_0,ierr)
      call mpi_allreduce(qinf,qinfall,1,mpi_double_precision,mpi_sum,var%com%comm_0,ierr)

 
!      call mpi_bcast(qoutall, 1, mpi_double_precision,var%com%np_a(3)-1,var%com%comm_a(3),ierr)     
!      call mpi_bcast(qinfall, 1, mpi_double_precision,0                ,var%com%comm_a(3),ierr)     

      qinf=qinfall
      qout=qoutall
      cor=(qinfall-qoutall)/red
      cor=(qinfall-0.0   )/(redall+1.e-15)
    do  jc=1,var%msh%nrc
!      write (  *,"(1I5,100(F20.8))") jc,var%su%corr3n(1,jc),var%su%corr3s(1,jc),var%msh%qc(jc),t0,t1,cor,var%msh%qc(jc),red,qinf,var%su%tconv(1,jc),cor
       do  ic=1,var%msh%ntheta
        var%su%corr3n(ic,jc)=var%su%corr3n(ic,jc)+var%msh%qc(jc)*cor
       enddo
!      write (20+var%par%irkstep+3*var%com%ip,"(1I5,100(F20.8))") jc,var%su%corr3n(1,jc),var%su%corr3s(1,jc),var%msh%qc(jc),redall
!      if (jc.eq.1) write (  *,"(1I5,100(F20.8))") jc,qinfall,qoutall,var%su%corr3n(1,jc),var%su%corr3s(1,jc),var%msh%qc(jc),t0,t1,cor,var%msh%qc(jc),red,qinf,var%su%tconv(1,jc),cor,var%su%q3(1,jc,var%msh%nzc),var%su%q3(1,jc,var%msh%nz),var%su%dq3(1,jc,var%msh%nzc)
     enddo

!


   end subroutine hdnbc

end module momentum_terms
