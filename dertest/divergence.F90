module divergence

!#include <petsc/finclude/petscdef.h>
  use run_tools
  use cons_tools
  use para_tools
  use comm_tools
  use mesh_tools
  use solu_tools
  use tripar
  use communicate


  implicit none
  contains
  subroutine divq(var)

    implicit none
    integer      :: ieq
    type (run)   :: var

    integer ic,jc,kc
    integer ip,im,jm


    REAL(KIND=8)  :: qnn,qss,dqcap1,dqcap2,dqcap3
    REAL(KIND=8)  :: d11q1,d33q1,d22q1
    REAL(KIND=8)  :: d11q2,d33q2,d22q2
    REAL(KIND=8)  :: d11q3,d33q3,d22q3

    do ic=1,var%msh%ntheta
      ip=var%msh%ipv(ic)
      do jc=1,var%msh%nr
        do kc=1,var%msh%nz
!          var%su%dq2(ic,jc+1,kc)=0.0
!          var%su%dq2(ic,jc,kc)=0.0
          var%sp%rhs(ic,jc,kc)=-(                                                                                                      &
                                 (var%su%dq1(ip,jc,kc)-var%su%dq1(ic,jc,kc))*var%msh%dx1*var%msh%axsym/var%msh%rm(jc)                 &
                                +(var%su%dq2(ic,jc+1,kc)-var%su%dq2(ic,jc,kc))*var%msh%dx2/(var%msh%rm(jc)*var%msh%dmr(jc))           &
                                +(var%su%dq3(ic,jc,kc+1)-var%su%dq3(ic,jc,kc))*var%msh%dx3                                            &
                               )/(var%par%crkalm(var%par%irkstep)*var%par%dt)
          
!          var%sp%rhs(ic,jc,kc)=0.0
!          if (kc.eq.1) var%sp%rhs(ic,jc,kc)=-1.0
!          if (kc.eq.var%msh%nz) var%sp%rhs(ic,jc,kc)=1.0
!
!          if (var%com%ip.eq.0) var%sp%rhs(ic,jc,1)=-1.0

!                                 dqcap1= (var%su%dq1(ip,jc,kc)-var%su%dq1(ic,jc,kc))*var%msh%dx1*var%msh%axsym/var%msh%rm(jc)    
!                                 dqcap2= (var%su%dq2(ic,jc+1,kc)-var%su%dq2(ic,jc,kc))*var%msh%dx2/(var%msh%rm(jc)*var%msh%dmr(jc))           
!                                 dqcap3= (var%su%dq3(ic,jc,kc+1)-var%su%dq3(ic,jc,kc))*var%msh%dx3                                            
!                                 if (kc.eq.1) var%sp%rhs(ic,jc,kc)=0.d0 
!                                 if (kc.eq.var%msh%nz) var%sp%rhs(ic,jc,kc)=0.d0 
                                                          
!          var%sp%rhs(ic,jc,kc)=cos(0.01*float(kc))
!           if (ic.eq.1) write (550+3*var%par%istep+var%par%irkstep,"(100(F30.16))") var%msh%rm(jc),var%msh%zm(kc),var%su%q1(1,jc,kc),var%su%q2(1,jc,kc),var%su%q3(1,jc,kc),var%su%dq1(1,jc,kc),var%su%dq2(1,jc,kc),var%su%dq3(1,jc,kc),var%sp%rhs(ic,jc,kc),dqcap1,dqcap2,dqcap3,var%su%dq3(ic,jc,kc+1),var%su%dq3(ic,jc,kc+1)-var%su%dq3(ic,jc,kc)
        enddo
!        if (ic.eq.1) write (550+3*var%par%istep+var%par%irkstep,"(100(F16.8))")
      enddo
    enddo
!    stop -1
 
!  call commq(var,var%sp%rhs,3)
!  call commq(var,var%sp%rhs,2)
    call communicateeq(var)

  end subroutine divq

  subroutine upsol(var)

    implicit none
    integer      :: ieq
    type (run)   :: var

    integer ic,jc,kc,ii
    integer im


    REAL(KIND=8)  :: ag,uin

      do  kc=1,var%msh%nzc
      do  jc=1,var%msh%nr
      do  ic=1,var%msh%ntheta
!        read (1000+var%par%irkstep+3*var%par%istep,"(3I5,100(E30.16))") ii,ii,ii,var%sp%dph(ic,jc,kc)
!        write (2000+var%par%irkstep,"(3I5,100(E20.8))") ic,jc,kc,var%sp%dph(ic,jc,kc)
      enddo
!      read (1000+var%par%irkstep+3*var%par%istep,"(3I5,100(F20.8))")
      enddo
      enddo



    if (var%com%ip_a(3).eq.var%com%np_a(3)-1) then
    do ic=1,var%msh%ntheta
      do jc=1,var%msh%nr
          var%sp%dph(ic,jc,var%msh%nzc)= var%sp%dph(ic,jc,var%msh%nz)
      enddo
    enddo
    endif

    call commq(var,var%sp%dph,3)
!    call communicateeq(var)


    do ic=1,var%msh%ntheta
      im=var%msh%imv(ic)
      do jc=1,var%msh%nr
        do kc=1,var%msh%nzc
           var%su%q1(ic,jc,kc)=var%su%dq1(ic,jc,kc)-(var%par%crkalm(var%par%irkstep)*var%par%dt*var%msh%dx1/var%msh%rm(jc)) *(var%sp%dph(ic,jc,kc)-var%sp%dph(im,jc,kc))
        enddo
      enddo
    enddo

    call communicateeq(var)
    if (var%com%ip_a(2).eq.var%com%np_a(2)-1) then    
    do ic=1,var%msh%ntheta
      do kc=1,var%msh%nzc
        var%su%q1(ic,var%msh%nrc,kc)=var%par%slip*var%su%q1(ic,var%msh%nr,kc)
      enddo 
    enddo 
    endif 
    
    call commq(var,var%su%q1,3)
    call commq(var,var%su%q1,2)

!    call communicateeq(var)

    if (var%com%ip_a(2).eq.0) then
      do ic=1,var%msh%ntheta
        do jc=2,var%msh%nr
          do kc=1,var%msh%nzc
            var%su%q2(ic,jc,kc)=var%su%dq2(ic,jc,kc)-(var%par%crkalm(var%par%irkstep)*var%par%dt)*var%msh%rc(jc)*var%msh%dx2*(var%sp%dph(ic,jc,kc)-var%sp%dph(ic,var%msh%jmv(jc),kc))/var%msh%dcr(jc)
         !   var%su%q2(ic,jc,kc)=0.d0
          enddo
        enddo
      enddo
    else
      do ic=1,var%msh%ntheta
        do jc=1,var%msh%nr
          do kc=1,var%msh%nzc
            var%su%q2(ic,jc,kc)=var%su%dq2(ic,jc,kc)-(var%par%crkalm(var%par%irkstep)*var%par%dt)*var%msh%rc(jc)*var%msh%dx2*(var%sp%dph(ic,jc,kc)-var%sp%dph(ic,var%msh%jmv(jc),kc))/var%msh%dcr(jc)
         !   var%su%q2(ic,jc,kc)=0.d0
          enddo
        enddo
      enddo
    endif

    do kc=1,var%msh%nzc
      do ic=1,var%msh%ntheta
       if (var%com%ip_a(2).eq.0) var%su%q2(ic,1,kc)=0.d0
       if (var%com%ip_a(2).eq.var%com%np_a(2)-1) var%su%q2(ic,var%msh%nrc,kc)=0.d0
      enddo
    enddo

    call communicateeq(var)
    do ic=1,var%msh%ntheta
      do jc=0,var%msh%nrc 
        var%su%q2s(ic,jc)=var%su%q2s(ic,jc)!var%su%corr2s(ic,jc)
      enddo
    enddo

    call commq(var,var%su%q2,3)
    call commq(var,var%su%q2,2)
!    call communicateeq(var)
 

    if (var%com%ip_a(3).eq.0) then
      do kc=2,var%msh%nzc
        do jc=1,var%msh%nr
          do ic=1,var%msh%ntheta
            var%su%q3(ic,jc,kc)= var%su%dq3(ic,jc,kc)-(var%par%crkalm(var%par%irkstep)*var%par%dt) * var%msh%dx3*( var%sp%dph(ic,jc,kc)- var%sp%dph(ic,jc,kc-1))
          enddo
        enddo
      enddo
    else
      do kc=1,var%msh%nzc
        do jc=1,var%msh%nr
          do ic=1,var%msh%ntheta
            var%su%q3(ic,jc,kc)= var%su%dq3(ic,jc,kc)-(var%par%crkalm(var%par%irkstep)*var%par%dt) * var%msh%dx3*( var%sp%dph(ic,jc,kc)- var%sp%dph(ic,jc,kc-1))
          enddo
        enddo
      enddo
    endif
!   pressure field
!    call communicateeq(var)

    call commq(var,var%su%q3,3)
    call commq(var,var%su%q3,2)


    uin= 0.01*(3*(var%par%ntime/0.5)**2.0 - 2*(var%par%ntime/0.5)**3.0)
    if (var%par%ntime.gt.0.5) uin=0.0100
    do jc=0,var%msh%nrc
      do ic=1,var%msh%ntheta
        if (var%com%ip_a(3).eq.0)  then
           var%su%q3(ic,jc,1)=var%su%dq3(ic,jc,1)   ! vitesse imposee, deja prise en compte
!           if (var%msh%rc(jc).lt.0.49) var%su%q3s(ic,jc)=var%su%q3s(ic,jc)+uin                  ! on rajoute la contrib de l'injection
        endif
      enddo
    enddo


    if (var%com%ip_a(3).eq.0) then
      do ic=1,var%msh%ntheta
        do kc=2,var%msh%nzc
         if (var%com%ip_a(2).eq.var%com%np_a(2)-1) var%su%q3(ic,var%msh%nrc,kc)=var%par%slip*var%su%q3(ic,var%msh%nr,kc)
        enddo
      enddo 
    else
      do ic=1,var%msh%ntheta
        do kc=1,var%msh%nzc
          if (var%com%ip_a(2).eq.var%com%np_a(2)-1) var%su%q3(ic,var%msh%nrc,kc)=var%par%slip*var%su%q3(ic,var%msh%nr,kc)
        enddo
      enddo 
    endif

    do kc=1,var%msh%nzc
      do jc=1,var%msh%nr
        do ic=1,var%msh%ntheta
!      if (ic.eq.2) write (650+3*var%par%istep+var%par%irkstep,"(3I5,100(F20.8))") 2,jc,kc,var%su%q1(2,jc,kc),var%su%q2(2,jc,kc),var%su%q3(2,jc,kc),var%su%dq1(2,jc,kc),var%su%dq2(2,jc,kc),var%su%dq3(2,jc,kc),var%sp%dph(ic,jc,kc)
    enddo 
!      if (ic.eq.2) write (650+3*var%par%istep+var%par%irkstep,"(3I5,100(F20.8))") 
    enddo 
    enddo 
!--> Reactualisation du scalaire passif                
!
!      do 50 ic=1,n1m
!        do 50 jc=1,n2m
!          scasud(ic,jc)=scasud(ic,jc)+corrss(ic,jc)
!        enddo
!      enddo
    call commq(var,var%su%q3,3)
   call commq(var,var%su%q3,2)
!    call communicateeq(var)

  end subroutine upsol


  subroutine uppres(var)

    implicit none
    integer      :: ieq
    type (run)   :: var

    integer ic,jc,kc
    integer ip,jp,kp
    integer im,jm,km,ij,ik


    REAL(KIND=8)  :: d11q,d33q,d22q
    REAL(KIND=8)  :: qnn,qss,be


    do ic=1,var%msh%ntheta
      do jc=1,var%msh%nr
          var%sp%dph(ic,jc,var%msh%nzc)= var%sp%dph(ic,jc,var%msh%nz)
      enddo
    enddo


!  CALCUL DU GRADIENT DE PRESSION         
!  Gp^l+1 - Gp^l = GPhi-alf*dt/2*A(GPhi)   
   be=0.5d0*(var%par%crkalm(var%par%irkstep)*var%par%dt)/var%par%Re


!  pressure gradient theta


!  Composante du gradient suivant theta
   if(var%msh%ntheta.gt.1) then
!  Gradient du Phi suivant theta GPhi= 1/r*d/dth(Phi) en i,j+1/2,k+1/2
      do ic=1,var%msh%ntheta
        im=var%msh%imv(ic)
        do jc=1,var%msh%nr
          do kc=1,var%msh%nzc
            var%su%qcap(ic,jc,kc)=(var%sp%dph(ic,jc,kc)-var%sp%dph(im,jc,kc))*var%msh%dx1/var%msh%rm(jc)
          enddo
        enddo
      enddo

    call commq(var,var%su%qcap,3)
    call commq(var,var%su%qcap,2)
!    call communicateeq(var)

!   Calcul du A(GPhi)
      do kc=1,var%msh%nz ! decide
        do jc=1,var%msh%nr
          do ic=1,var%msh%ntheta
            km=var%msh%kmv(kc)
            kp=kc+1

            jm=var%msh%jmv(jc)
            jp=var%msh%jpv(jc)

            ip=var%msh%ipv(ic)
            im=var%msh%imv(ic)
!   Calcul du A1th(GPhi)=1/r^2*d^2/dth^2(GPhi) en i,j+1/2,k+1/2
            d11q=(var%su%qcap(ip,jc,kc)-2.d0*var%su%qcap(ic,jc,kc)+var%su%qcap(im,jc,kc))*var%su%s1t11(jc)

!   Calcul du A1z(GPhi)=d^2/dz^2(GPhi) en i,j+1/2,k+1/2
            d33q=(var%su%qcap(ic,jc,kp)-2.d0*var%su%qcap(ic,jc,kc)+var%su%qcap(ic,jc,km))*var%msh%dx3q 
!   Calcul du A1r(GPhi)=1/r^2d/dr(r^3*d/dr(GPhi/r))
            qnn=(var%su%qcap(ic,jp,kc)/var%msh%rm(jc+1) -var%su%qcap(ic,jc,kc)/var%msh%rm(jc))
            qss=(var%su%qcap(ic,jc,kc)/var%msh%rm(jc) -var%su%qcap(ic,jm,kc)/var%msh%rm(jm))
            d22q=(qnn*var%su%s1t2n(jc)-qss*var%su%s1t2s(jc))*var%su%s1t2c(jc)
!   Nouveau gradient de pression, composante theta
            var%su%prg(ic,jc,kc,1)=var%su%prg(ic,jc,kc,1)+var%su%qcap(ic,jc,kc)-be*(d11q+d33q+d22q)
          enddo
        enddo
      enddo
    endif

!   Composante du gradient suivant r     

!   Gradient du Phi suivant r GPhi= r*d/dr(Phi) en i+1/2,j,k+1/2
    do ic=1,var%msh%ntheta
      do kc=1,var%msh%nz
        if (var%com%ip_a(2).eq.                0) var%su%qcap(ic,1,kc)=0.d0
        if (var%com%ip_a(2).eq.var%com%np_a(2)-1) var%su%qcap(ic,var%msh%nrc,kc)=0.d0
      enddo
    enddo

    if (var%com%ip_a(2).eq.0) then
      do kc=1,var%msh%nzc
        do jc=2,var%msh%nr
          do ic=1,var%msh%ntheta
            var%su%qcap(ic,jc,kc)=(var%sp%dph(ic,jc,kc)-var%sp%dph(ic,jc-1,kc))* var%msh%dx2*var%msh%rc(jc)/var%msh%dcr(jc)    
          enddo
        enddo
      enddo
    else
      do kc=1,var%msh%nzc
        do jc=1,var%msh%nr
          do ic=1,var%msh%ntheta
            var%su%qcap(ic,jc,kc)=(var%sp%dph(ic,jc,kc)-var%sp%dph(ic,jc-1,kc))* var%msh%dx2*var%msh%rc(jc)/var%msh%dcr(jc)    
          enddo
        enddo
      enddo
    endif

    call commq(var,var%su%qcap,3)
    call commq(var,var%su%qcap,2)
!    call communicateeq(var)
!    call commq(var,var%sp%prg,3)
!    call commq(var,var%sp%prg,2)
!    call communicateeq(var)
!   pressure gradient r

!   Calcul du A(GPhi)
      ij=1;  if (var%com%ip_a(2).eq.0) ij=2
      do kc=1,var%msh%nz ! decide
        do jc=ij,var%msh%nr  !decide
          do ic=1,var%msh%ntheta
            km=var%msh%kmv(kc)
            kp=kc+1
            ip=var%msh%ipv(ic)
            im=var%msh%imv(ic)
!   Calcul du A2th(GPhi)=1/r^2*d^2/dth^2(GPhi) en i+1/2,j,k+1/2
            d11q=(var%su%qcap(ip,jc,kc)-2.d0*var%su%qcap(ic,jc,kc)+var%su%qcap(im,jc,kc))*var%su%s2t11(jc)
!   Calcul du A2z(GPhi)=d^2/dz^2(GPhi) en i+1/2,j,k+1/2
            d33q=(var%su%qcap(ic,jc,kp)-2.d0*var%su%qcap(ic,jc,kc)+var%su%qcap(ic,jc,km))*var%msh%dx3q
!   Calcul du A2r(GPhi)=r*d/dr(1/r*d/dr(GPhi)) en i+1/2,j,k+1/2
            qnn=var%su%qcap(ic,jc+1,kc)-var%su%qcap(ic,jc,kc)
            qss=var%su%qcap(ic,jc,kc)-var%su%qcap(ic,jc-1,kc)
            d22q=qnn*var%su%s2t2p(jc)-qss*var%su%s2t2m(jc)
!   Nouveau gradient de pression suivant r
            var%su%prg(ic,jc,kc,2)=var%su%prg(ic,jc,kc,2)+var%su%qcap(ic,jc,kc)-be*(d11q+d22q+d33q)
!        if (ic.eq.1) write (650+10*var%com%ip+var%par%irkstep,"(1I5,100(F20.8))") 2,var%msh%rc(jc),var%msh%zc(kc),var%su%q1(1,jc,kc),var%su%q2(1,jc,kc),var%su%q3(1,jc,kc),var%su%dq1(1,jc,kc),var%su%dq2(1,jc,kc),var%su%dq3(1,jc,kc),qnn,qss,var%su%qcap(ic,jc,kp)

        enddo
      enddo
!     if (ic.eq.1) write (650+10*var%com%ip+var%par%irkstep,"(3I5,100(F20.8))") 
    enddo

!   pressure gradient z


!   Composante du gradient suivant z     
!   Gradient du Phi suivant z GPhi= d/dz(Phi) en i+1/2,j+1/2,k

    ij=1;  if (var%com%ip_a(2).eq.0) ij=2

    do jc=ij,var%msh%nr
      do ic=1,var%msh%ntheta
        if (var%com%ip_a(3).eq.                0) var%su%qcap(ic,jc,1)=0.d0
      enddo
    enddo

 
    ik=1;  if (var%com%ip_a(3).eq.0) ik=2

    do kc=ik,var%msh%nzc !decide
      do jc=1,var%msh%nr
        do ic=1,var%msh%ntheta
          var%su%qcap(ic,jc,kc)=(var%sp%dph(ic,jc,kc)-var%sp%dph(ic,jc,kc-1))* var%msh%dx3
          if (isnan(var%su%qcap(ic,jc,kc))) then
                  write (*,*) "Flow field crash"
                  stop
          endif
        enddo
      enddo
    enddo

    call commq(var,var%su%qcap,3)
    call commq(var,var%su%qcap,2)
    !call communicateeq(var)
    !call communicateeq(var)

    ik=1;  if (var%com%ip_a(3).eq.0) ik=2
    do kc=ik,var%msh%nz ! decide
      do jc=1,var%msh%nr
        do ic=1,var%msh%ntheta
          km=var%msh%kmv(kc)
          kp=kc+1

            jm=var%msh%jmv(jc)
            jp=var%msh%jpv(jc)


          ip=var%msh%ipv(ic)
          im=var%msh%imv(ic)

!   Calcul du A3th(GPhi)=1/r^2*d^2/dth^2(GPhi) en i+1/2,j+1/2,k
          d11q=(var%su%qcap(ip,jc,kc)-2.d0*var%su%qcap(ic,jc,kc)+var%su%qcap(im,jc,kc))*var%su%s1t11(jc)
!   Calcul du A3z(GPhi)=d^2/dz^2(GPhi)  en i+1/2,j+1/2,k
          d33q=(var%su%qcap(ic,jc,kc+1)-2.d0*var%su%qcap(ic,jc,kc)+var%su%qcap(ic,jc,kc-1))*var%msh%dx3q
!   Calcul du A3r(GPhi)=1/r(d/dr(r*d/dr(GPhi))) en i+1/2,j+1/2,k
          qnn=(var%su%qcap(ic,jp,kc)-var%su%qcap(ic,jc,kc))
          qss=(var%su%qcap(ic,jc,kc)-var%su%qcap(ic,jm,kc))

          d22q=qnn*var%su%s3t2n(jc)-qss*var%su%s3t2s(jc)
!   Nouveau gradient sur z
          var%su%prg(ic,jc,kc,3)=var%su%prg(ic,jc,kc,3)+var%su%qcap(ic,jc,kc)-be*(d11q+d22q+d33q)

          

        enddo
      enddo
    enddo

!   CALCUL DE LA PRESSION COMME SCALAIRE   p^l+1 - p^l = Phi-alf*dt/2*L(Phi)
!   Calcul du L(Phi)
!    call communicateeq(var)

    do kc=1,var%msh%nz  !decide
      do jc=1,var%msh%nr  !decide
        do ic=1,var%msh%ntheta
          km=var%msh%kmv(kc)
          kp=kc+1

          jm=var%msh%jmv(jc)
          jp=var%msh%jpv(jc)

          ip=var%msh%ipv(ic)
          im=var%msh%imv(ic)

!   Calcul de 1/r^2*d^2/dth^2(Phi)
          d11q=(var%sp%dph(ip,jc,kc)+var%sp%dph(im,jc,kc)-2.d0*var%sp%dph(ic,jc,kc))*var%su%s1t11(jc)
!   Calcul de 1/r*d/dr(r*d/dr(Phi))
          qnn=var%sp%dph(ic,jp,kc)-var%sp%dph(ic,jc,kc)
          qss=var%sp%dph(ic,jc,kc)-var%sp%dph(ic,jm,kc)
          d22q=qnn*var%su%s3t2n(jc)-qss*var%su%s3t2s(jc)
!   Calcul de d^2/dz^2(Phi)
          d33q=(var%sp%dph(ic,jc,kp)+var%sp%dph(ic,jc,km)-2.d0*var%sp%dph(ic,jc,kc))*var%msh%dx3q
!   Calcul de la pression
          var%su%pres(ic,jc,kc)=var%su%pres(ic,jc,kc)+var%sp%dph(ic,jc,kc)-be*(d11q+d22q+d33q)
!              integer iproc=var%com%ip
!              if (iproc.eq.1.and.iz.eq.1.and.ir.eq.32) write (*,*) 'press 1',var%su%pres(ith,ir,iz),var%sp%dph(ith,ir,iz)
!              if (iproc.eq.1.and.iz.eq.2.and.ir.eq.32) write (*,*) 'press 2',var%su%pres(ith,ir,iz),var%sp%dph(ith,ir,iz)
!              if (iproc.eq.1.and.iz.eq.3.and.ir.eq.32) write (*,*) 'press 3',var%su%pres(ith,ir,iz),var%sp%dph(ith,ir,iz)
!              if (var%com%ip.eq.1.and.kc.eq.1.and.jc.eq.32) write (*,*) 'press 1',var%su%pres(1,jc,kc),var%sp%dph(1,jc,kc),var%sp%dph(1,jc,km),km
!              if (var%com%ip.eq.1.and.kc.eq.2.and.jc.eq.32) write (*,*) 'press 2',var%su%pres(1,jc,kc),var%sp%dph(1,jc,kc),var%sp%dph(1,jc,km),km
!              if (var%com%ip.eq.1.and.kc.eq.3.and.jc.eq.32) write (*,*) 'press 3',var%su%pres(1,jc,kc),var%sp%dph(1,jc,kc),var%sp%dph(1,jc,km),km
        enddo     
      enddo     
    enddo     

    call commq(var,var%su%pres,2)
    call commq(var,var%su%pres,3)
  end subroutine uppres
end module divergence
