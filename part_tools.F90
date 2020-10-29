module part_type
  implicit none
  save

  type :: part

! Data
    real(kind=8), allocatable      :: u(:,:), ddt(:,:)
    real(kind=8), dimension(3)     :: vel
    real(kind=8), dimension(3,3)   :: du
    real(kind=8), dimension(3,3,3) :: ddu
    real(kind=8), dimension(3,3)   :: jac, jrt
    real(kind=8), dimension(3,3,3) :: hes, hrt
    real(kind=8), dimension(27,3)  :: xint
    real(kind=8), dimension(27,3)  :: yint
    real(kind=8), dimension(27)    :: rint
    real(kind=8), dimension(9,27)  :: mint
    real(kind=8), dimension(9)     :: acoef
    real(kind=8), dimension(9)     :: mp


! Parameters
    integer      :: npart
    integer      :: nvar
    integer      :: ipos, ivel, ijac,ijrt,ihes,ihrt
    integer      :: nint
    real(kind=8) :: st
  end type part

end module part_type


module part_tools
  use part_type
  use mesh_type
  use solu_type
  use comm_type
  use para_type
  use para_tools
  use cons_tools

  implicit none

  contains


  subroutine alct_part(prt)
    implicit none
    type(part)        :: prt

    allocate(prt%u(prt%npart,prt%nvar))
    allocate(prt%ddt(prt%npart,prt%nvar))

  end subroutine alct_part

  subroutine init_part(prt,msh)
    implicit none
    type(part)       :: prt
    type(mesh_a)     :: msh

    real(kind=8)     :: x,y,z,inpos,dx,dy,dz
    integer          :: ip,ip1

    prt%ipos=0
    prt%ivel=3
    prt%ijac=6
    prt%ijrt=15
    prt%ihes=24
    prt%ihrt=51
    
    prt%nvar=78
    prt%npart=10000
    prt%st=0.1d0

    call alct_part(prt)
    prt%u(:,:)=0.d0
    
    ip=0  
    ip1=0
    do  while (ip.le.prt%npart)
      ip1=ip1+1
      if (ip1.gt.1000000000) write (*,*) "error while count exceeds one billion"
    
           call random_number(inpos)
       !   rnd=2.0*rand(52358)-1.0
           x=(2.0*inpos-1.0)*msh%rmax
       !   rnd=2.0*rand(66535)-1.0
           call random_number(inpos)
       !    y=(2*inpos-1)*msh%rmax
           y=inpos*msh%rmax
           y=0.d0
       !   rnd=2.0*rand(21524)-1.0
           call random_number(inpos)
           z=(2.0*inpos-1.0)*msh%zmax
          
     
        !   if (x**2+y**2.lt.msh%rmax**2) then
           if (x**2+(z-1.0)**2.lt.0.64.and.x.lt.0.0)  then
       
                   ip=ip+1
          
                   prt%u(ip,prt%ipos+1)=x
                   prt%u(ip,prt%ipos+2)=y
                   prt%u(ip,prt%ipos+3)=z
           endif  
    enddo 
    dx=msh%dtheta
    dy=msh%dr
    dz=msh%dz
    prt%acoef(:)=0.d0

    prt%xint=reshape(              &
             (/  -dx,  -dy,  -dz , &
                0.d0,  -dy,  -dz , &
                  dx,  -dy,  -dz , &
                 -dx, 0.d0,  -dz , &
                0.d0, 0.d0,  -dz , &
                  dx, 0.d0,  -dz , &
                 -dx,   dy,  -dz , &
                0.d0,   dy,  -dz , &
                  dx,   dy,  -dz , &
                 -dx,  -dy, 0.d0 , &
                0.d0,  -dy, 0.d0 , &
                  dx,  -dy, 0.d0 , &
                 -dx, 0.d0, 0.d0 , &
                0.d0, 0.d0, 0.d0 , &
                  dx, 0.d0, 0.d0 , &
                 -dx,   dy, 0.d0 , &
                0.d0,   dy, 0.d0 , &
                  dx,   dy, 0.d0 , &
                 -dx,  -dy,   dz , &
                0.d0,  -dy,   dz , &
                  dx,  -dy,   dz , &
                 -dx, 0.d0,   dz , &
                0.d0, 0.d0,   dz , &
                  dx, 0.d0,   dz , &
                 -dx,   dy,   dz , &
                0.d0,   dy,   dz , &
                  dx,   dy,   dz   &
             /),(/27,3/))

    prt%mint(1,:)= (/ &
    -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx), -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx), -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx), -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx), -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx), -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx), -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx), -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx), -(1.d0/(18.d0*dx)), 0.d0, 1.d0/(18.d0*dx)  & 
                  /)

    prt%mint(2,:)=(/ &
            -(1.d0/(18.d0*dy)), -(1.d0/(18.d0*dy)), -(1.d0/(18.d0*dy)), 0.d0, 0.d0, 0.d0, 1.d0/(18.d0*dy), 1.d0/(18.d0*dy), 1.d0/(18.d0*dy), -(1.d0/(18.d0*dy)), -(1.d0/(18.d0*dy)), -(1.d0/(18.d0*dy)), 0.d0, 0.d0, 0.d0, 1.d0/(18.d0*dy), 1.d0/(18.d0*dy), 1.d0/(18.d0*dy), -(1.d0/(18.d0*dy)), -(1.d0/(18.d0*dy)), -(1.d0/(18.d0*dy)), 0.d0, 0.d0, 0.d0, 1.d0/(18.d0*dy), 1.d0/(18.d0*dy), 1.d0/(18.d0*dy) &
                  /)

    prt%mint(3,:)=(/ &
     -(1.d0/(18.d0*dz)), -(1.d0/(18.d0*dz)), -(1.d0/(18.d0*dz)), -(1.d0/(18.0*dz)), -(1.d0/(18.d0*dz)), -(1.d0/(18.d0*dz)), -(1.d0/(18.d0*dz)), -(1.d0/(18.d0*dz)), -(1.d0/(18.d0*dz)), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/(18.d0*dz), 1.d0/(18.d0*dz), 1.d0/(18.d0*dz), 1.d0/(18.d0*dz), 1.d0/(18.d0*dz), 1.d0/(18.d0*dz), 1.d0/(18.d0*dz), 1.d0/(18.d0*dz), 1.d0/(18.d0*dz) &
                  /)

    prt%mint(4,:)=(/ &
      1.d0/(42.d0*dx**2), -(2.d0/(21.d0*dx**2)), 1.d0/(42.d0*dx**2), 1.d0/(14.d0*dx**2), -(1.0/(21.d0*dx**2)), 1.d0/(14.d0*dx**2), 1.d0/(42.d0*dx**2),-(2.d0/(21.d0*dx**2)), 1.d0/(42.d0*dx**2), 1.d0/(14.d0*dx**2), -(1.d0/(21.d0*dx**2)), 1.d0/(14.d0*dx**2), 5.d0/(42.d0*dx**2), 0.d0, 5.0/(42.d0*dx**2), 1.d0/(14.d0*dx**2), -(1.d0/(21.d0*dx**2)), 1.d0/(14.d0*dx**2), 1.d0/(42.d0*dx**2), -(2.d0/(21.d0*dx**2)), 1.d0/(42.d0*dx**2), 1.d0/(14.d0*dx**2), -(1.d0/(21.0*dx**2)), 1.d0/(14.d0*dx**2), 1.d0/(42.d0*dx**2), -(2.d0/(21.d0*dx**2)), 1.d0/(42.d0*dx**2) &
                  /)

    prt%mint(5,:)=(/ &
      1.d0/(12.d0*dx*dy), 0.d0, -(1.d0/(12.d0*dx*dy)), 0.d0, 0.d0, 0.d0, -(1.d0/(12.d0*dx*dy)), 0.d0, 1.d0/(12.d0*dx*dy), 1.d0/(12.d0*dx*dy), 0.d0, -(1.d0/(12.d0*dx*dy)), 0.d0, 0.d0, 0.d0, -(1.d0/(12.d0*dx*dy)), 0.d0, 1.d0/(12.d0*dx*dy), 1.d0/(12.d0*dx*dy), 0.d0, -(1.d0/(12.d0*dx*dy)), 0.d0, 0.d0, 0.d0, -(1.d0/(12.d0*dx*dy)), 0.d0, 1.d0/(12.d0*dx*dy) &
                  /)

    prt%mint(6,:)=(/ &
      1.d0/(42.d0*dy**2), 1.d0/(14.d0*dy**2), 1.d0/(42.d0*dy**2), -(2.d0/(21.d0*dy**2)), -(1.d0/(21.d0*dy**2)), -(2.d0/(21.0*dy**2)), 1.d0/(42.d0*dy**2), 1.d0/(14.d0*dy**2), 1.d0/(42.d0*dy**2), 1.d0/(14.d0*dy**2), 5.d0/(42.d0*dy**2), 1.d0/(14.d0*dy**2), -(1.d0/(21.d0*dy**2)), 0.d0, -(1.d0/(21.d0*dy**2)), 1.d0/(14.d0*dy**2), 5.d0/(42.d0*dy**2), 1.d0/(14.d0*dy**2), 1.d0/(42.d0*dy**2), 1.d0/(14.d0*dy**2), 1.d0/(42.d0*dy**2), -(2.d0/(21.d0*dy**2)), -(1.d0/(21.d0*dy**2)), -(2.d0/(21.d0*dy**2)), 1.d0/(42.d0*dy**2), 1.d0/(14.d0*dy**2), 1.d0/(42.d0*dy**2)  &
                  /)

    prt%mint(7,:)=(/ &
      1.d0/(12.d0*dy*dz), 1.d0/(12.d0*dy*dz), 1.d0/(12.d0*dy*dz), 0.d0, 0.d0, 0.d0, -(1.d0/(12.d0*dy*dz)), -(1.d0/(12.d0*dy*dz)), -(1.d0/(12.d0*dy*dz)), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(1.d0/(12.d0*dy*dz)), -(1.d0/(12.d0*dy*dz)), -(1.d0/(12.d0*dy*dz)), 0.d0, 0.d0, 0.d0, 1.d0/(12.d0*dy*dz), 1.d0/(12.d0*dy*dz), 1.d0/(12.d0*dy*dz) &
                  /)

    prt%mint(8,:)=(/ &
      1.d0/(42.d0*dz**2), 1.d0/(14.d0*dz**2), 1.d0/(42.d0*dz**2), 1.d0/(14.d0*dz**2), 5.d0/(42.d0*dz**2), 1.d0/(14.d0*dz**2), 1.d0/(42.d0*dz**2), 1.d0/(14.d0*dz**2), 1.d0/(42.d0*dz**2), -(2.d0/(21.d0*dz**2)), -(1.d0/(21.d0*dz**2)), -(2.d0/(21.d0*dz**2)), -(1.d0/(21.d0*dz**2)), 0.d0, -(1.d0/(21.d0*dz**2)), -(2.d0/(21.d0*dz**2)), -(1.d0/(21.d0*dz**2)), -(2.d0/(21.d0*dz**2)), 1.d0/(42.d0*dz**2), 1.d0/(14.d0*dz**2), 1.d0/(42.d0*dz**2), 1.d0/(14.d0*dz**2), 5.d0/(42.d0*dz**2), 1.d0/(14.d0*dz**2), 1.d0/(42.d0*dz**2), 1.d0/(14.d0*dz**2), 1.d0/(42.d0*dz**2) &
                  /)
    prt%mint(9,:)=(/ &
      1.d0/(12.d0*dx*dz), 0.d0, -(1.d0/(12.d0*dx*dz)), 1.d0/(12.d0*dx*dz), 0.d0, -(1.d0/(12.d0*dx*dz)), 1.d0/(12.d0*dx*dz), 0.d0, -(1.d0/(12.d0*dx*dz)), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(1.d0/(12.d0*dx*dz)), 0.d0, 1.d0/(12.d0*dx*dz), -(1.d0/(12.d0*dx*dz)), 0.d0, 1.d0/(12.d0*dx*dz), -(1.d0/(12.d0*dx*dz)), 0.d0, 1.d0/(12.d0*dx*dz) &
                  /)
    prt%nint=100
  end subroutine init_part
  subroutine interpolation (xint,acoef,rint)
         implicit none
         real(kind=8) :: xint(27,3)
         real(kind=8) :: yint(27)
         real(kind=8) :: rint(27)
         real(kind=8) :: acoef(9)

         rint=xint(:,1)*acoef(1)+xint(:,2)*acoef(2)+xint(:,3)*acoef(3) &
             +xint(:,1)*xint(:,1)*acoef(4)+xint(:,1)*xint(:,2)*acoef(5)+xint(:,2)*xint(:,2)*acoef(6) &
             +xint(:,2)*xint(:,3)*acoef(7)+xint(:,3)*xint(:,3)*acoef(8)+xint(:,3)*xint(:,1)*acoef(9)

  end subroutine interpolation
         
         
 subroutine iter_part(prt,msh,com,su,par)
    implicit none
    type(part)      :: prt
    type(mesh_a)    :: msh
    type(comm)      :: com
    type(solu)      :: su
    type(para)      :: par


    integer :: it, ir, iz
    integer :: ip,i,j,k
    integer :: ifst1,ifst2,ifst3
    integer :: iint
    integer :: irc,izc
    integer :: irof,itof,izof
    integer :: iv,ivar

    real(kind=8) :: dt
    real(kind=8) :: x,y,z,r,theta,ut,ur,uz
    real(kind=8) :: tp
    real(kind=8) :: sumj,sumh
    real(kind=8) :: drp,dtp,dzp
    real(kind=8) :: fx,fy,f
    real(kind=8) :: ft,ftt,ftr,ftz
    real(kind=8) :: fr,frt,frr,frz
    real(kind=8) :: fz,fzt,fzr,fzz

    character(len=8)          :: fmt5 ! format descriptor
    character(5)              :: timechar

            fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
          write (timechar,fmt5) par%nstep
          
          if (mod((par%nstep),par%expstp).eq.0) then
           open (unit=22,file="part"//trim(timechar)//".dat", form='formatted', position='rewind')
            write (22,*) 'title = "particles"'
            write (22,*) 'variables = "xp", "yp", "zp", "up", "vp", "wp" '
          endif
        
        do ip=1,prt%npart
          x=prt%u(ip,prt%ipos+1)
          y=prt%u(ip,prt%ipos+2)
          z=prt%u(ip,prt%ipos+3)
          r=dsqrt(x**2+y**2)
          if (x.ge.0.d0.and.y.ge.0.d0) then 
             theta=dacos(x/r)
          elseif (x.lt.0.d0.and.y.ge.0.d0) then  
             theta=dacos(x/r)
          elseif (x.lt.0.d0.and.y.lt.0.d0) then  
             theta=dacos(-x/r)+4.0*datan(1.d0)
          elseif (x.ge.0.d0.and.y.lt.0.d0) then  
             theta=dacos(-x/r)+4.0*datan(1.d0)
          endif  
           
          it=nint(theta/msh%dtheta)+1-com%ip_a(1)*msh%ntheta
          ir=nint(r/msh%dr)+1-com%ip_a(2)*msh%nr
          iz=nint(z/msh%dz)+1-com%ip_a(3)*msh%nz
          if (it.ge.1.and.it.le.msh%ntheta.and.ir.ge.1.and.ir.le.msh%nr.and.iz.ge.1.and.iz.le.msh%nz) then 
          ut=su%q1(it,ir,iz)
          ur=su%q2(it,ir,iz)/msh%rm(ir)
          uz=su%q3(it,ir,iz)
          irc=0
          izc=0

          if (ir.ge.msh%nr) irc=-1
          if (ir.le.1) irc=1
          if (iz.ge.msh%nz) izc=-1
          if (iz.le.1) izc=1

          
          do ifst1=-1,1,1
            do ifst2=-1,1,1
              do ifst3=-1,1,1
                 iv=(ifst1+1)+3*(ifst2+1)+9*(ifst3+1)
                 itof=it+ifst1
                 if (itof.gt.msh%ntheta) itof=1
                 if (itof.lt.1) itof=msh%ntheta
                 irof=ir+ifst2+irc
                 izof=iz+ifst3+izc
                 
                 prt%yint(iv,1)=su%q1(itof,irof,izof)
                 prt%yint(iv,2)=su%q2(itof,irof,izof)/msh%rc(irof)
                 prt%yint(iv,3)=su%q3(itof,irof,izof)
              enddo
            enddo
          enddo
          do ivar=1,3
            do iint=1,prt%nint
              call interpolation (prt%xint,prt%acoef,prt%rint)
              prt%mp=-matmul(prt%mint,prt%rint-prt%yint(:,ivar))
              prt%acoef=prt%acoef-prt%mp
              if (norm2(prt%mp).lt.0.000001) exit
            enddo
          
            if (ivar.eq.1) then
              dtp=theta-msh%thc(it)
              drp=r-msh%rm(it)
              dzp=z-msh%zm(it)
            endif
              
            if (ivar.eq.2) then
              dtp=theta-msh%thm(it)
              drp=r-msh%rc(it)
              dzp=z-msh%zm(it)
            endif
              
            if (ivar.eq.3) then
              dtp=theta-msh%thm(it)
              drp=r-msh%rm(it)
              dzp=z-msh%zc(it)
            endif
            if (dtp.gt.8.d0*datan(1.d0)) dtp=dtp-8.d0*datan(1.d0)
            if (dtp.lt.-8.d0*datan(1.d0)) dtp=dtp+8.d0*datan(1.d0)
             f=dtp*prt%acoef(1)+drp*prt%acoef(2)+dzp*prt%acoef(3)+dtp*dtp*prt%acoef(4)+dtp*drp*prt%acoef(5)+drp*drp*prt%acoef(6)+drp*dzp*prt%acoef(7)+dzp*dzp*prt%acoef(8)+dzp*dtp*prt%acoef(9)
             ft=prt%acoef(1)+2.0*dtp*prt%acoef(4)+drp*prt%acoef(5)+dzp*prt%acoef(9)
             ftt=2.0*prt%acoef(4)
             ftr=drp*prt%acoef(5)
             ftz=dzp*prt%acoef(9)
             fr=prt%acoef(2)+dtp*prt%acoef(5)+2.0*drp*prt%acoef(6)+dzp*prt%acoef(7)
             frt=prt%acoef(5)
             frr=2.0*prt%acoef(6)
             frz=prt%acoef(7)
             fz=prt%acoef(3)+drp*prt%acoef(7)+2.0*dzp*prt%acoef(8)+dtp*prt%acoef(9)
             fzt=prt%acoef(9)
             fzr=prt%acoef(7)
             fzz=2.0*prt%acoef(8)

             fx=dcos(theta)*fr-sin(theta)*ft/msh%rm(ir)
             fy=dsin(theta)*fr+cos(theta)*ft/msh%rc(ir)
          enddo

          prt%vel(1)=ur*dcos(theta)-ut*dsin(theta)
          prt%vel(2)=ur*dsin(theta)+ut*dcos(theta)
          prt%vel(3)=uz
          ! 2D problem
     !    prt%vel(1)=ur
     !    prt%vel(2)=0.d0
     !    prt%vel(3)=uz
  else
          prt%vel(1)=0.d0
          prt%vel(2)=0.d0
          prt%vel(3)=0.d0
  endif
        !  prt%du(1,1)=(su%q1(it+1,ir,iz)-su%q1(it-1,ir,iz))/(msh%dtheta*msh%rm(ir)*2.0)
          prt%du(:,:)=0.d0
          prt%ddu(:,:,:)=0.d0

          do i=1,3
                prt%ddt(ip,prt%ipos+i)=prt%u(ip,prt%ivel+i)
                prt%ddt(ip,prt%ivel+i)=(1.0/prt%st)*(prt%vel(i)-prt%u(ip,prt%ivel+i))
          enddo
          do i=1,3
            do j=1,3
              prt%ddt(ip,prt%ijac+i+3*j-3)=prt%u(ip,prt%ijrt+i+3*j-3)
              sumj=prt%jac(1,j)*prt%du(i,1)+prt%jac(2,j)*prt%du(i,2)+prt%jac(3,j)*prt%du(i,3)
              prt%ddt(ip,prt%ijrt+i+3*j-3)=(1.0/prt%st)*(sumj-prt%u(ip,prt%ijac+i+3*j-3))
          enddo
        enddo

        do i=1,3
          do j=1,3
            do k=1,3
              prt%ddt(ip,prt%ihes+i+3*j-3+9*k-9)=prt%u(ip,prt%ihrt+i+3*j-3+9*k-9)
              sumh=prt%hes(1,j,k)*prt%du(i,1)+prt%hes(2,j,k)*prt%du(i,2)+prt%hes(3,j,k)*prt%du(i,3)+prt%jac(1,j)*prt%jac(1,k)*prt%ddu(i,1,1)+prt%jac(2,j)*prt%jac(1,k)*prt%ddu(i,2,1)+prt%jac(3,j)*prt%jac(1,k)*prt%ddu(i,3,1)+prt%jac(1,j)*prt%jac(2,k)*prt%ddu(i,1,2)+prt%jac(1,j)*prt%jac(3,k)*prt%ddu(i,1,3)+prt%jac(2,j)*prt%jac(2,k)*prt%ddu(i,2,2)+prt%jac(2,j)*prt%jac(3,k)*prt%ddu(i,2,3)+prt%jac(3,j)*prt%jac(2,k)*prt%ddu(i,3,2)+prt%jac(3,j)*prt%jac(3,k)*prt%ddu(i,3,3)
       
              prt%ddt(ip,prt%ihrt+i+3*j-3+9*k-9)=(1.0/prt%st)*(sumh-prt%u(ip,prt%ihrt+i+3*j-3+9*k-9))
            enddo
          enddo
        enddo
        do i=1,6
                prt%u(ip,i)=prt%u(ip,i)+prt%ddt(ip,i)*par%dt
        enddo
                if (mod((par%nstep),par%expstp).eq.0) then
          write (22,*) prt%u(ip,1),prt%u(ip,2),prt%u(ip,3),prt%u(ip,4),prt%u(ip,5),prt%u(ip,6) 
          endif
       enddo

          if (mod((par%nstep),par%expstp).eq.0) close(22) 

  end subroutine iter_part



end module part_tools
