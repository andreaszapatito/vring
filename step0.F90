! module convective

module inverse

  use run_tools
  use cons_tools
  use para_tools
  use comm_tools
  use mesh_tools
  use solu_tools
  use gen_tools
  use tripar

  implicit none
  contains
  subroutine step(var)

    implicit none
    integer      :: ieq
    type (run)   :: var

    integer ic,jc,kc,ji,ki
    integer ip,im,jm


    REAL(KIND=8)  :: qnn,qss,r,po2
    REAL(KIND=8)  :: d11q1,d33q1,d22q1
    REAL(KIND=8)  :: d11q2,d33q2,d22q2
    REAL(KIND=8)  :: d11q3,d33q3,d22q3
      !write (*,*) 'in step'
    do ic=1,var%msh%ntheta
      ip=var%msh%ipv(ic)
      im=var%msh%imv(ic)
      do jc=1,var%msh%nr
        jm=var%msh%jmv(jc)
        do kc=1,var%msh%nz

          d11q1=(var%su%q1(ip,jc,kc)-2.0d0*var%su%q1(ic,jc,kc)+var%su%q1(im,jc,kc))*var%su%s1t11(jc)
          if(var%com%ip_a(3).eq.0.and.kc.eq.1) then
            d33q1=4.d0/3.d0*(var%su%q1(ic,jc,kc+1)-3.d0*var%su%q1(ic,jc,kc)+2.d0*var%su%q1s(ic,jc))*var%msh%dx3q
          else
            d33q1=(var%su%q1(ic,jc,kc+1)-2.0d0*var%su%q1(ic,jc,kc)+var%su%q1(ic,jc,kc-1))*var%msh%dx3q
          endif
          qnn=var%su%q1(ic,jc+1,kc)/var%msh%rm(jc+1)-var%su%q1(ic,jc,kc)/var%msh%rm(jc)
          qnn=qnn*var%su%s1t2n(jc)

          qss=var%su%q1(ic,jc,kc)/var%msh%rm(jc)-var%su%q1(ic,jm,kc)/var%msh%rm(jm)
          qss=qss*var%su%s1t2s(jc)

          d22q1=(qnn-qss)*var%su%s1t2c(jc)

!         RHS of First pass of the factorisation AF1 for q1
    
          var%su%qrhs(ic,jc,kc)=(var%par%crkgam(var%par%irkstep)*var%su%dq1(ic,jc,kc)+var%par%crkrom(var%par%irkstep)*var%su%ru1(ic,jc,kc)-var%par%crkalm(var%par%irkstep)*var%su%prg(ic,jc,kc,1)+var%par%crkalm(var%par%irkstep)/var%par%Re*(d11q1+d33q1+d22q1))*var%par%dt
!          write (var%com%ip+300,"(10(E16.8))") var%su%qrhs(ic,jc,kc),var%su%dq1(ic,jc,kc),d11q1,d33q1,d22q1,var%su%q1(ic,jc,kc+1),var%su%q1(ic,jc,kc),var%su%q1(ic,jc,kc-1),var%msh%dx3q,var%su%q1s(ic,jc)

 
          var%su%ru1(ic,jc,kc)=var%su%dq1(ic,jc,kc)

        enddo
      enddo
    enddo


    call solqi(1,var)
    call solqj(1,var)
    call solqk(1,var)

    do ic=1,var%msh%ntheta
      do kc=1,var%msh%nzc
        if(var%com%ip_a(2).eq.var%com%np_a(2)-1) var%su%dq1(ic,var%msh%nrc,kc)=var%par%slip*var%su%dq1(ic,var%msh%nr,kc)
      enddo
    enddo


    ki=1; if(var%com%ip_a(3).eq.0) ki=2

    do ic=1,var%msh%ntheta
      ip=var%msh%ipv(ic)
      im=var%msh%imv(ic)
      do kc=ki,var%msh%nz
        do jc=1,var%msh%nr

          d11q3=(var%su%q3(ip,jc,kc)-2.0d0*var%su%q3(ic,jc,kc)+var%su%q3(im,jc,kc))*var%su%s1t11(jc)
          d33q3=(var%su%q3(ic,jc,kc+1)-2.0d0*var%su%q3(ic,jc,kc)+var%su%q3(ic,jc,kc-1))*var%msh%dx3q

          qnn=(var%su%q3(ic,jc+1,kc)-var%su%q3(ic,jc,kc))*var%su%s3t2n(jc)
          qss=(var%su%q3(ic,jc,kc)-var%su%q3(ic,var%msh%jmv(jc),kc))*var%su%s3t2s(jc)

          r=var%msh%rm(jc)
          po2=2.0*atan(1.0)
          d22q3=qnn-qss !-(po2*sin(po2*r)/r+po2**2.d0*cos(po2*r))

!         RHS of First pass of the factorisation AF1 for q3

          var%su%qrhs(ic,jc,kc)=(var%par%crkgam(var%par%irkstep)*var%su%dq3(ic,jc,kc)+var%par%crkrom(var%par%irkstep)*var%su%ru3(ic,jc,kc)-var%par%crkalm(var%par%irkstep)*var%su%prg(ic,jc,kc,3)+(var%par%crkalm(var%par%irkstep)/var%par%Re)*(d11q3+d33q3+d22q3))*var%par%dt
          var%su%ru3(ic,jc,kc)=var%su%dq3(ic,jc,kc)

        enddo
      enddo
    enddo



    call solqi(3,var)
    call solqj(3,var)
    call solqk(3,var)
   
    do ic=1,var%msh%ntheta
      do kc=1,var%msh%nzc
        if(var%com%ip_a(2).eq.var%com%np_a(2)-1) var%su%dq3(ic,var%msh%nrc,kc)=var%par%slip*var%su%dq3(ic,var%msh%nr,kc)
      enddo
    enddo



    ji=1; if(var%com%ip_a(2).eq.0) ji=2

    do ic=1,var%msh%ntheta
      ip=var%msh%ipv(ic)
      im=var%msh%imv(ic)
      do kc=1,var%msh%nz
        do jc=ji,var%msh%nr

          d11q2=(var%su%q2(ip,jc,kc)-2.0d0*var%su%q2(ic,jc,kc)+var%su%q2(im,jc,kc))*var%su%s2t11(jc)
          if(var%com%ip_a(3).eq.0.and.kc.eq.1) then
            d33q2=4.d0/3.d0*(var%su%q2(ic,jc,kc+1)-3.d0*var%su%q2(ic,jc,kc)+2.d0*var%su%q2s(ic,jc))*var%msh%dx3q
          else
            d33q2=(var%su%q2(ic,jc,kc+1)-2.0d0*var%su%q2(ic,jc,kc)+var%su%q2(ic,jc,kc-1))*var%msh%dx3q
          endif
          d22q2=var%su%q2(ic,jc+1,kc)*var%su%s2t2p(jc)-var%su%q2(ic,jc,kc)*var%su%s2t2c(jc)+var%su%q2(ic,jc-1,kc)*var%su%s2t2m(jc)

!         RHS of First pass of the factorisation AF1 for q2
          var%su%qrhs(ic,jc,kc)=(var%par%crkgam(var%par%irkstep)*var%su%dq2(ic,jc,kc)+var%par%crkrom(var%par%irkstep)*var%su%ru2(ic,jc,kc)-var%par%crkalm(var%par%irkstep)*var%su%prg(ic,jc,kc,2)+var%par%crkalm(var%par%irkstep)/var%par%Re*(d11q2+d33q2+d22q2))*var%par%dt

          var%su%ru2(ic,jc,kc)=var%su%dq2(ic,jc,kc)

        enddo
      enddo
    enddo

    call solqi(2,var)
    call solqk(2,var)
    call solqj(2,var)
  
    do ic=1,var%msh%ntheta
      do kc=1,var%msh%nzc
        if(var%com%ip_a(2).eq.var%com%np_a(2)-1) var%su%dq2(ic,var%msh%nrc,kc)=var%par%slip*var%su%dq2(ic,var%msh%nr,kc)
      enddo
    enddo

!
  do jc=1,var%msh%nr
    do kc=1,var%msh%nz
!      write (500+3*var%par%istep+var%par%irkstep+var%com%ip,"(100(F20.8))") var%msh%thc(2),var%msh%rc(jc),var%msh%zc(kc) ,var%su%dq1(2,jc,kc),var%su%dq2(2,jc,kc),var%su%dq3(2,jc,kc),var%su%q1(2,jc,kc),var%su%q2(2,jc,kc),var%su%q3(2,jc,kc),var%su%qrhs(2,jc,kc),d11q2+d22q2+d33q2,var%su%ru2(2,jc,kc),var%su%prg(2,jc,kc,2)
    enddo
!    write (500+3*var%par%istep+var%par%irkstep+var%com%ip,"(100(E16.8))")
  enddo


!      do jc=1,var%msh%nr
!        do kc=1,var%msh%nz
!          write (500+var%com%ip,"(100(F16.8))") var%msh%thc(2),var%msh%rc(jc),var%msh%zc(kc),var%su%q1(2,jc,kc),var%su%q2(2,jc,kc),var%su%q3(2,jc,kc),var%su%pres(2,jc,kc),var%sp%dph(2,jc,kc),var%su%dq3(2,jc,kc),var%su%stri(kc),var%su%rtri(kc),var%su%atri(kc),var%su%btri(kc),var%su%ctri(kc),var%su%qrhs(2,jc,kc),var%su%prg(2,jc,kc,3)
!        enddo
!        write (500+var%com%ip,"(100(F16.8))")
!        write (*,*) "hi there"
!      enddo





  end subroutine step

  subroutine solqi(ieq,var)
    implicit none
    integer      :: ieq
    type (run)   :: var

    integer i,j,k
    integer ip,im,jm
    integer ik,ij,itrip

    REAL(KIND=8)  :: al,betadx,acc
    betadx=0.5d0*var%par%crkalm(var%par%irkstep)*var%par%dt/var%par%Re

!   First Pass of the Factorisation AF1
    itrip=1
    if (ieq.eq.1) then 
      ik=1 
      ij=1
    elseif (ieq.eq.2) then 
      ik=1 
      ij=1; if(var%com%ip_a(2).eq.0) ij=2
      itrip=1 !caution
    elseif (ieq.eq.3) then 
      ik=1; if(var%com%ip_a(3).eq.0) ik=2
      ij=1
    endif
       
 
    do k=ik,var%msh%nz
      do j=ij,var%msh%nr
        do i=1,var%msh%nthetac
!          write (var%com%ip+400,"(10(E16.8))") var%su%qrhs(i,j,k)
        enddo
        if (ieq.eq.1) acc=betadx*var%su%s1t11(j)
        if (ieq.eq.2) acc=betadx*var%su%s2t11(j)
        if (ieq.eq.3) acc=betadx*var%su%s1t11(j)
        do i=1,var%msh%nthetac
          var%su%ami(j,i)=-acc
          var%su%api(j,i)=-acc
          var%su%aci(j,i)=1.d0+2.d0*acc
!          write (*,*) jc,var%msh%nrc,var%com%ip,var%com%ip_a(2)
          var%su%fi(j,i)=var%su%qrhs(i,j,k)
        enddo
      enddo
      call tripij(1,var%msh%ntheta,itrip,var%msh%nr,var)
      do j=ij,var%msh%nr
        do i=1,var%msh%ntheta
          var%su%qrhs(i,j,k)=var%su%fi(j,i)
!         write (*,"(A,I5,10(F16.8))") "step 1st",ieq, var%su%qrhs(i,j,k)
        enddo
      enddo
    enddo
  end subroutine solqi


  subroutine solqj(ieq,var)

  implicit none
  integer      :: ieq,ik,fk,ij
  type (run)   :: var

  integer ic,jc,kc,ierr,nrloc

  REAL(KIND=8)  :: al


    if     (ieq.eq.1) then
      ik=1
      ij=1
      fk=var%msh%nz
      nrloc=var%msh%nr
    elseif (ieq.eq.2) then
      ik=1
      ij=1; if(var%com%ip_a(2).eq.0) ij=1
      fk=var%msh%nzc
      nrloc=var%msh%nrc ! important
    elseif (ieq.eq.3) then
      ik=1; if(var%com%ip_a(3).eq.0) ik=1
      ij=1;
      fk=var%msh%nz
      nrloc=var%msh%nr ! caution
!      nrloc=var%msh%nr
    endif


  do ic=1,var%msh%ntheta
    do kc=ik,fk
      do jc=ij,var%msh%nrc
       var%su%rtri(jc)=var%su%qrhs(ic,jc,kc)
       if (ieq.eq.1) then
         var%su%atri(jc)=var%su%amj1(jc)
         var%su%btri(jc)=var%su%acj1(jc)
         var%su%ctri(jc)=var%su%apj1(jc)
       elseif (ieq.eq.2) then
         var%su%atri(jc)=var%su%amj2(jc)
         var%su%btri(jc)=var%su%acj2(jc)
         var%su%ctri(jc)=var%su%apj2(jc)
       elseif (ieq.eq.3) then
         var%su%atri(jc)=var%su%amj3(jc)
         var%su%btri(jc)=var%su%acj3(jc)
         var%su%ctri(jc)=var%su%apj3(jc)
       endif 
      enddo

      if (ieq.eq.2) then
        if (var%com%ip_a(2).eq.0)  then 
          var%su%rtri(1)=0.0
          var%su%atri(1)=0.0
          var%su%btri(1)=1.0
          var%su%ctri(1)=0.0
        endif
        if (var%com%ip_a(2).eq.var%com%np_a(2)-1)  then
          var%su%rtri(var%msh%nrc)=0.0
          var%su%atri(var%msh%nrc)=0.0
          var%su%btri(var%msh%nrc)=1.0
          var%su%ctri(var%msh%nrc)=0.0
        endif
      endif 
      call solve_lap_par(var%su%rtri,var%su%stri,var%su%atri,var%su%btri,var%su%ctri,nrloc,nrloc,var%com%comm_a(2))
      if (ieq.eq.2) then
        do jc=1,var%msh%nr
          var%su%dq2(ic,jc,kc)=var%su%stri(jc)+var%su%q2(ic,jc,kc)
        enddo
      else
        do jc=1,var%msh%nr ! question
          var%su%qrhs(ic,jc,kc)=var%su%stri(jc)
!          if (ieq.eq.3.and.ic.eq.2) write (500+var%com%ip+3*var%par%istep+var%par%irkstep,"(100(E20.8))")  var%msh%thc(2),var%msh%rc(jc),var%msh%zc(kc),var%su%qrhs(ic,jc,kc),var%su%amj2(jc),var%su%acj2(jc),var%su%apj2(jc)
        enddo
!		if (ieq.eq.3.and.ic.eq.2) write (500+var%com%ip+3*var%par%istep+var%par%irkstep,"(100(E20.8))")
      endif

       if (ieq.eq.3) then
       if (ic.eq.1) then
       if (kc.eq.40) then
        do jc=1,var%msh%nr
!          write (500+var%com%ip*10+var%par%irkstep,"(100(F16.8))") var%msh%thc(1),var%msh%rc(jc),var%msh%zc(jc),var%su%q1(1,jc,kc),var%su%q2(1,jc,kc),var%su%q3(1,jc,kc),var%su%dq1(1,jc,kc),var%su%dq2(1,jc,kc),var%su%dq3(1,jc,kc),var%su%stri(jc),var%su%rtri(jc),var%su%atri(jc),var%su%btri(jc),var%su%ctri(jc),var%su%qrhs(1,jc,jc)
         enddo
!       write (500+var%com%ip,"(100(F16.8))")
!       write (*,*) "hi there"
       endif
       endif
       endif


    enddo

  enddo
  end subroutine solqj

  subroutine solqk (ieq,var)
  implicit none

  integer      :: ieq,ierr,nzloc
  type (run)   :: var

  integer ic,jc,kc,ij,ik,fk
  real po2
!  integer ip,im,jm
  REAL(KIND=8)  :: al,uin,uot,dt,ts,pi,cc
  po2=2.0*atan(1.0)
  pi=4.0*atan(1.0)
  
 ! write (*,*) 'In SOLQK',ieq,var%msh%ntheta,var%msh%nrc
!   Third Pass of the Factorisation AF3
      
    if     (ieq.eq.1) then
      ij=1
      dt=var%par%dt
      uin=(var%par%crkgam(var%par%irkstep)+var%par%crkrom(var%par%irkstep))*var%msh%rc(jc)*(profile(var%par%ntime)-profile(var%par%ntime-dt))
!      uin=0.d0
      uot=0.d0
      nzloc=var%msh%nzc
      ik=1; if(var%com%ip_a(3).eq.0) ik=1

      fk=var%msh%nzc
    elseif (ieq.eq.2) then
      ij=1; if(var%com%ip_a(2).eq.0) ij=1
      uin=0.d0
      uot=0.d0
      nzloc=var%msh%nzc ! francky visit
      ik=1
      fk=var%msh%nz
    elseif (ieq.eq.3) then
      ij=1
      dt=var%par%dt
      ts=0.5
!      uin= 1.0*((3.0*(var%par%ntime/ts)**2.0 - 2.0*(var%par%ntime/ts)**3.0)-(3.0*((var%par%ntime-dt)/ts)**2.0 - 2.0*((var%par%ntime-dt)/ts)**3.0))
       uin=(var%par%crkgam(var%par%irkstep)+var%par%crkrom(var%par%irkstep))*(profile(var%par%ntime)-profile(var%par%ntime-dt))
!      if (var%par%ntime.gt.ts) uin=0.0
       uot=uin*(1.0-var%msh%rc(jc)**2/1.0**2)*2.0*(0.2)**2
! +0.5*(sin(2.0*pi*10.0*var%par%ntime/ts)-sin(10.0*2.0*pi*(var%par%ntime-dt)/ts))
!      uin= (24.0*(var%par%ntime)**1.0 - 3.0*32.0*(var%par%ntime)**2.0)*var%par%dt
!      uin=0.01
!      if (uin.le.0.0) uin=0.0000
      nzloc=var%msh%nzc !caution  !changed this to avoid u component at R0
      ik=1; if(var%com%ip_a(3).eq.0) ik=1
      fk=var%msh%nzc
    endif

!  ik=1
!  write (*,*) 'In SOLQK II',ieq,ik,fk
  uin=0
  uot=0
  do ic=1,var%msh%ntheta
    do jc=ij,var%msh%nrc
      var%su%rtri(:)=0.0

      if (ieq.eq.1) then
      dt=var%par%dt
      uin=1.0*(var%par%crkgam(var%par%irkstep)+var%par%crkrom(var%par%irkstep))*var%msh%rc(jc)*(profile(var%par%ntime)-profile(var%par%ntime-dt))
      uot=0.d0
      endif  
  uin=0
  uot=0

      do kc=ik,fk
       var%su%rtri(kc)=var%su%qrhs(ic,jc,kc)
       if (ieq.eq.1) then
         var%su%atri(kc)=var%su%amk1(kc)
         var%su%btri(kc)=var%su%ack1(kc)
         var%su%ctri(kc)=var%su%apk1(kc)
       elseif (ieq.eq.2) then
         var%su%atri(kc)=var%su%amk2(kc)
         var%su%btri(kc)=var%su%ack2(kc)
         var%su%ctri(kc)=var%su%apk2(kc)
       elseif (ieq.eq.3) then
         var%su%atri(kc)=var%su%amk3(kc)
         var%su%btri(kc)=var%su%ack3(kc)
         var%su%ctri(kc)=var%su%apk3(kc)
       endif
      enddo
      cc=0.0
      if (var%par%istep.eq.1) cc=0.0
!      if (ieq.eq.3) then
      if(var%msh%rc(jc).le.0.1999.and.var%com%ip_a(3).eq.0) var%su%rtri(1)=0.0*uin                          ! quand Vtheta_sud impose (voir coesys.f)
!      if(var%msh%rc(jc).gt.0.1999.and.var%com%ip_a(3).eq.0) var%su%rtri(1)=0                    ! quand Vtheta_sud impose (voir coesys.f)
 !     if(ieq.eq.3.and.var%com%ip_a(3).eq.0) var%su%rtri(1)=cc*(1.0+cos(po2*var%msh%rc(jc)))    ! quand Vtheta_sud impose (voir coesys.f)
!      write (100+var%par%itime,*) jc,uin,uot
!      write (*,*) 'In SOLQK III',jc,uin,uot
      if(var%com%ip_a(3).eq.var%com%np_a(3)-1) var%su%rtri(var%msh%nzc)=+1.*uot                    ! quand Vtheta_sud impose (voir coesys.f)
!      if(ieq.eq.3.and.var%com%ip_a(3).eq.var%com%np_a(3)-1) var%su%rtri(var%msh%nzc)= cc*(1.0+ cos(po2*var%msh%rc(jc)))     ! quand Vtheta_sud impose (voir coesys.f)
!      if(ieq.eq.3.and.var%com%ip_a(3).eq.var%com%np_a(3)-1) var%su%rtri(var%msh%nzc)=+1.*var%su%corr3n(ic,jc)                    ! quand Vtheta_sud impose (voir coesys.f)
!      if(ieq.eq.3.and.var%com%ip_a(3).eq.0) var%su%rtri(var%msh%nzc)=+1.*var%su%corr3s(ic,jc)                    ! quand Vtheta_sud impose (voir coesys.f)
!      if(var%msh%rc(jc).lt.2.0.and.var%com%ip_a(3).eq.var%com%np_a(3)-1) var%su%rtri(var%msh%nz)=+1.*uot                    ! quand Vtheta_sud impose (voir coesys.f)
!      endif

      call solve_lap_par(var%su%rtri,var%su%stri,var%su%atri,var%su%btri,var%su%ctri,nzloc,nzloc,var%com%comm_a(3))

       if (ieq.eq.3) then
           do kc=1,var%msh%nzc
!             write (2001,*) jc,kc,var%su%stri(kc),var%su%rtri(kc)
           enddo
!           write (2001,*) 
       endif

      do kc=1,var%msh%nzc
       if(ieq.eq.1) var%su%dq1(ic,jc,kc)=var%su%q1(ic,jc,kc)+var%su%stri(kc)
       if(ieq.eq.2) var%su%qrhs(ic,jc,kc)=var%su%stri(kc)
       if(ieq.eq.3) var%su%dq3(ic,jc,kc)=var%su%q3(ic,jc,kc)+var%su%stri(kc)
     enddo
    enddo
  enddo

  end subroutine solqk

  subroutine tripij(i1,i2,mi,mf,var)

    implicit none
    integer      :: ieq
    type (run)   :: var

    integer i1,i2,mi,mf
    integer ib,ie,ibw,im,ip,j,i
    REAL(KIND=8)  :: za1,za2

    do j=mi,mf
      var%su%aci(j,i1)=1.d0/var%su%aci(j,i1)
      var%su%gam2(j,i1)=-var%su%ami(j,i1)*var%su%aci(j,i1)
      var%su%ami(j,i1)=var%su%fi(j,i1)*var%su%aci(j,i1)
    enddo
    ib=i1+1
    ie=i2-1
    do j=mi,mf
      do i=ib,ie
        im=i-1
        var%su%api(j,im)=var%su%api(j,im)*var%su%aci(j,im)
        var%su%aci(j,i)=1.d0/(var%su%aci(j,i)-var%su%ami(j,i)*var%su%api(j,im))
        var%su%gam2(j,i)=-var%su%ami(j,i)*var%su%gam2(j,im)*var%su%aci(j,i)
        var%su%ami(j,i)=(var%su%fi(j,i)-var%su%ami(j,i)*var%su%ami(j,im))*var%su%aci(j,i)
      enddo
      var%su%gam2(j,ie)=var%su%gam2(j,ie)-var%su%api(j,ie)*var%su%aci(j,ie)
    enddo
  
    ibw=i1+i2-1
    do j=mi,mf
      var%su%fi(j,ie)=var%su%ami(j,ie)
      var%su%aci(j,ie)=var%su%gam2(j,ie)
      do i=ib,ie
        ip=ibw-i
        var%su%fi(j,ip)=var%su%ami(j,ip)-var%su%api(j,ip)*var%su%fi(j,ip+1)
        var%su%aci(j,ip)=var%su%gam2(j,ip)-var%su%api(j,ip)*var%su%aci(j,ip+1)

      enddo
    enddo
    do j=mi,mf
      za1=var%su%fi(j,i2)-var%su%api(j,i2)*var%su%fi(j,i1)-var%su%ami(j,i2)*var%su%fi(j,ie)
      za2=var%su%aci(j,i2)+var%su%api(j,i2)*var%su%aci(j,i1)+var%su%ami(j,i2)*var%su%aci(j,ie)
      var%su%fi(j,i2)=za1/za2
      do i=i1,ie
        var%su%fi(j,i)=var%su%fi(j,i)+var%su%aci(j,i)*var%su%fi(j,i2)
      enddo
    enddo
  end subroutine tripij

end module inverse



