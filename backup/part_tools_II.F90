module part_type
  implicit none
  save
  type :: part

! Data
    real(kind=8), allocatable      :: u(:,:), ddt(:,:)
    real(kind=8), dimension(3)     :: vel
    real(kind=8), dimension(3,3)   :: dvel
    real(kind=8), dimension(3,3,3) :: ddvel
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
  save
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
    prt%st=1.0d0

    !call alct_part(prt)
    allocate(prt%u(prt%npart,prt%nvar))
    allocate(prt%ddt(prt%npart,prt%nvar))
    prt%u(:,:)=0.d0
    prt%ddt(:,:)=0.d0
    
    ip=0  
    ip1=0
    do  while (ip.lt.prt%npart)
      ip1=ip1+1
      if (ip1.gt.1000000000) write (*,*) "error while count exceeds one billion"
    
       !    call random_number(inpos)
       !   rnd=2.0*rand(52358)-1.0
           inpos=rand(0)
           x=(2.0*inpos-1.0)*msh%rmax
       !   rnd=2.0*rand(66535)-1.0
       !    call random_number(inpos)
           inpos=rand(0)
       !    y=(2*inpos-1)*msh%rmax
           y=inpos*msh%rmax
           y=0.d0
       !   rnd=2.0*rand(21524)-1.0
       !    call random_number(inpos)
           inpos=rand(0)
           z=(2.0*inpos-1.0)*msh%zmax
          
     
        !   if (x**2+y**2.lt.msh%rmax**2) then
           if (x**2+(z-1.0)**2.lt.0.25.and.x.lt.0.0)  then
       
                   ip=ip+1
          
                   prt%u(ip,prt%ipos+1)=x
                   prt%u(ip,prt%ipos+2)=y
                   prt%u(ip,prt%ipos+3)=z
                   prt%u(ip,prt%ijac+1)=1.d0
                   prt%u(ip,prt%ijac+5)=1.d0
                   prt%u(ip,prt%ijac+9)=1.d0
          if (ip.eq.12) write (*,*) 'jacobian ', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
           endif  
    enddo 
    dx=msh%dtheta
    dy=msh%dr
    dz=msh%dz
    prt%acoef(:)=0.d0
 
    if (ip.eq.12) write (*,*) 'part 12 ', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
    prt%xint( 1,:)= (/  -dx,  -dy,  -dz /)
    prt%xint( 2,:)= (/ 0.d0,  -dy,  -dz /)
    prt%xint( 3,:)= (/   dx,  -dy,  -dz /)
    prt%xint( 4,:)= (/  -dx, 0.d0,  -dz /)
    prt%xint( 5,:)= (/ 0.d0, 0.d0,  -dz /)
    prt%xint( 6,:)= (/   dx, 0.d0,  -dz /)
    prt%xint( 7,:)= (/  -dx,   dy,  -dz /)
    prt%xint( 8,:)= (/ 0.d0,   dy,  -dz /)
    prt%xint( 9,:)= (/   dx,   dy,  -dz /)
    prt%xint(10,:)= (/  -dx,  -dy, 0.d0 /)
    prt%xint(11,:)= (/ 0.d0,  -dy, 0.d0 /)
    prt%xint(12,:)= (/   dx,  -dy, 0.d0 /)
    prt%xint(13,:)= (/  -dx, 0.d0, 0.d0 /)
    prt%xint(14,:)= (/ 0.d0, 0.d0, 0.d0 /)
    prt%xint(15,:)= (/   dx, 0.d0, 0.d0 /)
    prt%xint(16,:)= (/  -dx,   dy, 0.d0 /)
    prt%xint(17,:)= (/ 0.d0,   dy, 0.d0 /)
    prt%xint(18,:)= (/   dx,   dy, 0.d0 /)
    prt%xint(19,:)= (/  -dx,  -dy,   dz /)
    prt%xint(20,:)= (/ 0.d0,  -dy,   dz /)
    prt%xint(21,:)= (/   dx,  -dy,   dz /)
    prt%xint(22,:)= (/  -dx, 0.d0,   dz /)
    prt%xint(23,:)= (/ 0.d0, 0.d0,   dz /)
    prt%xint(24,:)= (/   dx, 0.d0,   dz /)
    prt%xint(25,:)= (/  -dx,   dy,   dz /)
    prt%xint(26,:)= (/ 0.d0,   dy,   dz /)
    prt%xint(27,:)= (/   dx,   dy,   dz /)
  
    if (ip.eq.12) write (*,*) 'part 12 A', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
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
    if (ip.eq.12) write (*,*) 'part 12 D', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
  end subroutine init_part
  subroutine interpolation (xint,acoef,rint)
         implicit none
         integer i
         real(kind=8) :: xint(27,3)
         real(kind=8) :: yint(27)
         real(kind=8) :: rint(27)
         real(kind=8) :: acoef(9)
         do i=1,27
         rint(i)=xint(i,1)*acoef(1)+xint(i,2)*acoef(2)+xint(i,3)*acoef(3) &
             +xint(i,1)*xint(i,1)*acoef(4)+xint(i,1)*xint(i,2)*acoef(5)+xint(i,2)*xint(i,2)*acoef(6) &
             +xint(i,2)*xint(i,3)*acoef(7)+xint(i,3)*xint(i,3)*acoef(8)+xint(i,3)*xint(i,1)*acoef(9)
         enddo
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
    integer :: id,id1,id2,im

    real(kind=8) :: dt
    real(kind=8) :: Jdet
    real(kind=8) :: y0(3)
    real(kind=8) :: x,y,z,r,theta,ut,ur,uz
    real(kind=8) :: tp
    real(kind=8) :: sumj,sumh
    real(kind=8) :: drp,dtp,dzp
    real(kind=8) :: fx,fy,f
    real(kind=8) :: sn,sn2,cs,cs2
    real(kind=8) :: ff(3),dff(3,3),ddff(3,3,3)
    real(kind=8) :: ft,ftt,ftr,ftz
    real(kind=8) :: fr,frt,frr,frz
    real(kind=8) :: fz,fzt,fzr,fzz
    real(kind=8) :: M(12,12),FI(12)
    
    character(len=8)          :: fmt5 ! format descriptor
    character(5)              :: timechar
    ip=12
    if (ip.eq.12) write (*,*) 'part 12 E', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)

        fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
        write (timechar,fmt5) par%istep
         
          if (ip.eq.12) write (*,*) 'particle 12 ', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
        if (mod((par%istep),par%expstp).eq.0) then
          open (unit=22,file="part"//trim(timechar)//".dat", form='formatted', position='rewind')
          write (22,*) 'title = "particles"'
          write (22,*) 'variables = "xp", "yp", "zp", "up", "vp", "wp" '
        endif
        do ip=1,prt%npart
        if (ip.eq.12) write (*,*) 'part 12 F', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
          x=prt%u(ip,prt%ipos+1)
          y=prt%u(ip,prt%ipos+2)
          z=prt%u(ip,prt%ipos+3)
          if (ip.eq.1) write (*,*) 'jacobian ', prt%u(1,prt%ijac+1)
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
           
          do ivar=1,3
            if (ivar==1) then 
              it=nint(theta/msh%dtheta)+1-com%ip_a(1)*msh%ntheta
              if (msh%ntheta.eq.1) it=1
              ir=nint((r-0.5*msh%dr)/msh%dr)+1-com%ip_a(2)*msh%nr
              iz=nint((z-0.5*msh%dz)/msh%dz)+1-com%ip_a(3)*msh%nz
            endif

            if (ivar==2) then 
              it=nint((theta-0.5*msh%dtheta)/msh%dtheta)+1-com%ip_a(1)*msh%ntheta
              if (msh%ntheta.eq.1) it=1
              ir=nint(r/msh%dr)+1-com%ip_a(2)*msh%nr
              iz=nint((z-0.5*msh%dz)/msh%dz)+1-com%ip_a(3)*msh%nz
            endif

            if (ivar==3) then 
              it=nint((theta-0.5*msh%dtheta)/msh%dtheta)+1-com%ip_a(1)*msh%ntheta
              if (msh%ntheta.eq.1) it=1
              ir=nint((r-0.5*msh%dr)/msh%dr)+1-com%ip_a(2)*msh%nr
              iz=nint(z/msh%dz)+1-com%ip_a(3)*msh%nz
            endif

          if (ip.eq.12) write (*,*) 'particle A 12 ', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
            if (it.ge.1.and.it.le.msh%ntheta.and.ir.ge.1.and.ir.le.msh%nr.and.iz.ge.1.and.iz.le.msh%nz) then 
              irc=0
              izc=0
              if (ip.eq.12) write (*,*) 'particle 1', ut,ur,uz
              if (ir.ge.msh%nr) irc=-1
              if (ir.le.1) irc=1
              if (iz.ge.msh%nz) izc=-1
              if (iz.le.1) izc=1
  
            
              do ifst3=-1,1,1
                do ifst2=-1,1,1
                  do ifst1=-1,1,1
                    iv=1+(ifst1+1)+3*(ifst2+1)+9*(ifst3+1)
                    itof=it+ifst1
                    if (itof.gt.msh%ntheta) itof=1
                    if (itof.lt.1) itof=msh%ntheta
                    irof=ir+ifst2+irc
                    izof=iz+ifst3+izc
                   
                    if (ivar==1) prt%yint(iv,1)=su%q1(itof,irof,izof)
                    if (ivar==2) prt%yint(iv,2)=su%q2(itof,irof,izof)/msh%rc(irof)
                    if (ivar==3) prt%yint(iv,3)=su%q3(itof,irof,izof)
                    if (ivar==1) prt%yint(iv,1)=0.d0
!                    if (ivar==2) prt%yint(iv,2)=0.1d0
!                    if (ivar==3) prt%yint(iv,3)=0.1d0
                    if (ivar==2.and.irof==1) prt%yint(iv,2)=0.d0
            if (ip.eq.12) write (*,*) 'particle interp 1',msh%rc(irof),irof,ir,irc,ifst2
                    if (ip.eq.12) write (*,*) 'particle interp 1',iv,ut,ur,uz,itof,irof,izof,it,ir,iz,prt%yint(iv,1),prt%yint(iv,2),prt%yint(iv,3)
                  enddo
                enddo
              enddo
              y0(ivar)=prt%yint(14,ivar)
              do iv=1,27
                prt%yint(iv,ivar)=prt%yint(iv,ivar)-y0(ivar)
              enddo
              prt%acoef(:)=0.d0
             if (ip.eq.12) write (*,*) 'particle B 12 ', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
              do iint=1,prt%nint
                call interpolation (prt%xint,prt%acoef,prt%rint)
                do i=1,9
                  prt%mp(i)=0.0
                  do k=1,27
                    prt%mp(i)=prt%mp(i)-prt%mint(i,k)*(prt%rint(k)-prt%yint(k,ivar))
                  enddo
                enddo
                do i=1,9
                  prt%acoef(i)=prt%acoef(i)+prt%mp(i)
                enddo
                if (ip.eq.12) write (*,*) 'particle interp 2a  a',iint,prt%acoef
                if (ip.eq.12) write (*,*) 'particle interp 2c  r',iint,prt%rint
                if (ip.eq.12) write (*,*) 'particle interp 2d  y',iint,prt%yint
                do i=1,27
                  if (ip.eq.12) write (*,*) 'particle interp 2b  x',iint,(prt%xint(i,j),j=1,3)
                enddo
                do i=1,9
                  if (ip.eq.12) write (*,*) 'particle interp 2e  m',i,(prt%mint(i,j),j=1,27)
                enddo
                if (ip.eq.12) write (*,*) 'particle interp acoef',iint,prt%acoef
                if (ip.eq.12) write (*,*) 'particle interp mp',prt%mp
                if (norm2(prt%mp).lt.0.000001) exit
              enddo
            
              if (ivar.eq.1) then
                dtp=theta-msh%thc(it)
                drp=r-msh%rm(ir)
                dzp=z-msh%zm(iz)
              endif
              
              if (ivar.eq.2) then
                dtp=theta-msh%thm(it)
                drp=r-msh%rc(ir)
                dzp=z-msh%zm(iz)
              endif
                
              if (ivar.eq.3) then
                dtp=theta-msh%thm(it)
                drp=r-msh%rm(ir)
                dzp=z-msh%zc(iz)
              endif
              if (ip.eq.12) write (*,*) 'particle transf 1',dtp,drp,dzp,theta-msh%thc(it)
              if (dtp.gt.msh%dtheta) dtp=dtp-8.d0*datan(1.d0)
              if (dtp.lt.-msh%dtheta) dtp=dtp+8.d0*datan(1.d0)
              if (ip.eq.12) write (*,*) 'particle transf 2',dtp,drp,dzp
              if (ip.eq.12) write (*,*) 'particle acoef 2',prt%acoef(:)
              f=dtp*prt%acoef(1)+drp*prt%acoef(2)+dzp*prt%acoef(3)+dtp*dtp*prt%acoef(4)+dtp*drp*prt%acoef(5)+drp*drp*prt%acoef(6)+drp*dzp*prt%acoef(7)+dzp*dzp*prt%acoef(8)+dzp*dtp*prt%acoef(9)+y0(ivar)
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
 
              sn=dsin(theta)
              sn2=dsin(2.0*theta)
              cs=dcos(theta)
              cs2=dcos(2.0*theta)
              FI=(/ ft, fr, fz, ftt, ftr, ftz, frt, frr, frz, fzt, fzr, fzz /)
              if (ip.eq.1) write (*,*) 'particle transf c ff',ivar,f,y0(ivar)
              M(1,1:12)= (/    -sn/r,      cs,0.d0,       0.d0,    0.d0,0.d0,    0.d0,  0.d0,0.d0, 0.d0,0.d0,0.d0 /)
              M(2,1:12)= (/     cs/r,      sn,0.d0,       0.d0,    0.d0,0.d0,    0.d0,  0.d0,0.d0, 0.d0,0.d0,0.d0 /)
              M(3,1:12)= (/     0.d0,    0.d0,1.d0,       0.d0,    0.d0,0.d0,    0.d0,  0.d0,0.d0, 0.d0,0.d0,0.d0 /)
              M(4,1:12)= (/ sn2/r**2, sn*sn/r,0.d0, sn*sn/r**2,-cs*sn/r,0.d0,-cs*sn/r, cs*cs,0.d0, 0.d0,0.d0,0.d0 /)
              M(5,1:12)= (/-cs2/r**2,-cs*sn/r,0.d0,-cs*sn/r**2,-sn*sn/r,0.d0, cs*cs/r, cs*sn,0.d0, 0.d0,0.d0,0.d0 /)
              M(6,1:12)= (/     0.d0,    0.d0,0.d0,       0.d0,    0.d0,-sn/r,   0.d0,  0.d0,  cs, 0.d0,0.d0,0.d0 /)
              M(7,1:12)= (/-cs2/r**2,-cs*sn/r,0.d0,-cs*sn/r**2, cs*cs/r,0.d0,-sn*sn/r, cs*sn,0.d0, 0.d0,0.d0,0.d0 /)
              M(8,1:12)= (/-sn2/r**2, cs*cs/r,0.d0, cs*cs/r**2, cs*sn/r,0.d0, cs*sn/r, sn*sn,0.d0, 0.d0,0.d0,0.d0 /)
              M(9,1:12)= (/     0.d0,    0.d0,0.d0,       0.d0,    0.d0, cs/r,   0.d0,  0.d0,  sn, 0.d0,0.d0,0.d0 /)
              M(10,1:12)=(/     0.d0,    0.d0,0.d0,       0.d0,    0.d0,0.d0,    0.d0,  0.d0,0.d0,-sn/r,  cs,0.d0 /)
              M(11,1:12)=(/     0.d0,    0.d0,0.d0,       0.d0,    0.d0,0.d0,    0.d0,  0.d0,0.d0, cs/r,  sn,0.d0 /)
              M(12,1:12)=(/     0.d0,    0.d0,0.d0,       0.d0,    0.d0,0.d0,    0.d0,  0.d0,0.d0, 0.d0,0.d0,1.d0 /)
              if (ip.eq.12) write (*,*) 'particle transf 3',M
              if (ip.eq.12) write (*,*) 'particle 12 ff',f,r,theta
              ff(ivar)=f
              do id=1,3
                dff(ivar,id)=0.d0
                do im=1,12
                  dff(ivar,id)=dff(ivar,id)+M(id,im)*FI(im)
                enddo
              enddo
              do id1=1,3
                do id2=1,3
                  ddff(ivar,id1,id2)=0.d0
                  do im=1,12
                    ddff(ivar,id1,id2)=ddff(ivar,id1,id2)+M(3+id1+3*(id2-1),im)*FI(im)
                  enddo
                enddo
              enddo
              prt%vel(1)=ff(2)*dcos(theta)-ff(1)*dsin(theta)
              prt%vel(2)=ff(2)*dsin(theta)+ff(1)*dcos(theta)
              prt%vel(3)=ff(3)
            
              if (ip.eq.1) write (*,*) 'particle u,v,w', prt%vel(1),prt%vel(2),prt%vel(3)
              do id=1,3
                prt%dvel(1,id)=dff(2,id)*dcos(theta)-dff(1,id)*dsin(theta)
                prt%dvel(2,id)=dff(2,id)*dsin(theta)+dff(1,id)*dcos(theta)
                prt%dvel(3,id)=dff(3,id)
              enddo
              if (ip.eq.1) write (*,*) 'particle du,dv,dw', prt%dvel(1,id)
              if (ip.eq.1) write (*,*) 'particle du,dv,dw', prt%dvel(2,id)
              if (ip.eq.1) write (*,*) 'particle du,dv,dw', prt%dvel(3,id)
              do id1=1,3
                do id2=1,3
                  prt%ddvel(1,id1,id2)=ddff(2,id1,id2)*dcos(theta)-ddff(1,id1,id2)*dsin(theta)
                  prt%ddvel(2,id1,id2)=ddff(2,id1,id2)*dsin(theta)+ddff(1,id1,id2)*dcos(theta)
                  prt%ddvel(3,id1,id2)=ddff(3,id1,id2)
                enddo
              enddo
             if (ip.eq.12) write (*,*) 'particle DD12 ',ff(1),ff(2),ff(3),theta

            else
              prt%vel(:)=0.d0
              prt%dvel(:,:)=0.d0
              prt%ddvel(:,:,:)=0.d0
            endif
          enddo
             if (ip.eq.12) write (*,*) 'particle C 12 ', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2),&
             prt%u(ip,prt%ivel+3),prt%vel(1),prt%u(ip,prt%ivel+1)
          do i=1,3
            prt%ddt(ip,prt%ipos+i)=prt%u(ip,prt%ivel+i)
            prt%ddt(ip,prt%ivel+i)=(1.0/prt%st)*(prt%vel(i)-prt%u(ip,prt%ivel+i))
          enddo
          do i=1,3
            do j=1,3
            prt%jac(i,j)=prt%u(ip,prt%ijac+i+3*j-3)
            prt%jrt(i,j)=prt%u(ip,prt%ijrt+i+3*j-3)
            do k=1,3
              prt%hes(i,j,k)=prt%u(ip,prt%ihes+i+3*j-3+9*k-9)
              prt%hrt(i,j,k)=prt%u(ip,prt%ihrt+i+3*j-3+9*k-9)
            enddo
            enddo
          enddo

          do i=1,3
            do j=1,3
              prt%ddt(ip,prt%ijac+i+3*j-3)=prt%u(ip,prt%ijrt+i+3*j-3)
              sumj=prt%jac(1,j)*prt%dvel(i,1)+prt%jac(2,j)*prt%dvel(i,2)+prt%jac(3,j)*prt%dvel(i,3)
              prt%ddt(ip,prt%ijrt+i+3*j-3)=(1.0/prt%st)*(sumj-prt%u(ip,prt%ijrt+i+3*j-3))
            enddo
          enddo
!          write (*,*) 'jacobian b' &
!          ,ip,prt%u(ip,prt%ijac+9),Jdet,prt%ddt(ip,prt%ijac+9),prt%ddt(ip,prt%ijrt+9),prt%u(ip,prt%ijrt+9),sumj,prt%dvel(3,3)
!          write (*,*) 'sumj', sumj,prt%dvel(i,3),prt%dvel(3,3)
          do i=1,3
            do j=1,3
              do k=1,3
                prt%ddt(ip,prt%ihes+i+3*j-3+9*k-9)=prt%u(ip,prt%ihrt+i+3*j-3+9*k-9)
                sumh=prt%hes(1,j,k)*prt%dvel(i,1)+prt%hes(2,j,k)*prt%dvel(i,2)+prt%hes(3,j,k)*prt%dvel(i,3)+prt%jac(1,j)*prt%jac(1,k)*prt%ddvel(i,1,1)+prt%jac(2,j)*prt%jac(1,k)*prt%ddvel(i,2,1)+prt%jac(3,j)*prt%jac(1,k)*prt%ddvel(i,3,1)+prt%jac(1,j)*prt%jac(2,k)*prt%ddvel(i,1,2)+prt%jac(1,j)*prt%jac(3,k)*prt%ddvel(i,1,3)+prt%jac(2,j)*prt%jac(2,k)*prt%ddvel(i,2,2)+prt%jac(2,j)*prt%jac(3,k)*prt%ddvel(i,2,3)+prt%jac(3,j)*prt%jac(2,k)*prt%ddvel(i,3,2)+prt%jac(3,j)*prt%jac(3,k)*prt%ddvel(i,3,3)
       
                prt%ddt(ip,prt%ihrt+i+3*j-3+9*k-9)=(1.0/prt%st)*(sumh-prt%u(ip,prt%ihrt+i+3*j-3+9*k-9))
              enddo
            enddo
          enddo
         if (ip.eq.12) write (*,*) 'particle F 12 ', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3),prt%ddt(ip,prt%ivel+2)
!          write (*,*) 'jacobian c' &
!          ,ip,prt%u(ip,prt%ijac+9),Jdet,prt%ddt(ip,prt%ijac+9),prt%ddt(ip,prt%ijrt+9),prt%u(ip,prt%ijrt+9),sumj,prt%dvel(3,3)
          do i=1,prt%nvar
            prt%u(ip,i)=prt%u(ip,i)+prt%ddt(ip,i)*par%dt
          enddo
          Jdet=prt%u(ip,prt%ijac+1)*prt%u(ip,prt%ijac+5)*prt%u(ip,prt%ijac+9) &
              +prt%u(ip,prt%ijac+2)*prt%u(ip,prt%ijac+6)*prt%u(ip,prt%ijac+7) &
              +prt%u(ip,prt%ijac+4)*prt%u(ip,prt%ijac+8)*prt%u(ip,prt%ijac+3) &
              -prt%u(ip,prt%ijac+7)*prt%u(ip,prt%ijac+5)*prt%u(ip,prt%ijac+3) &
              -prt%u(ip,prt%ijac+8)*prt%u(ip,prt%ijac+6)*prt%u(ip,prt%ijac+1) &
              -prt%u(ip,prt%ijac+4)*prt%u(ip,prt%ijac+2)*prt%u(ip,prt%ijac+9) 

!          write (*,*) 'jacobian d' &
!          ,ip,prt%u(ip,prt%ijac+1),Jdet,prt%ddt(ip,prt%ijac+9),prt%ddt(ip,prt%ijrt+9),prt%u(ip,prt%ijrt+9),sumj,prt%dvel(3,3)
                 
          if (mod((par%istep),par%expstp).eq.0) then
            write (22,*) prt%u(ip,1),prt%u(ip,2),prt%u(ip,3),prt%u(ip,4),prt%u(ip,5),prt%u(ip,6),Jdet
          endif
         if (ip.eq.12) write (*,*) 'particle G 12 ', prt%u(ip,prt%ivel+1), prt%u(ip,prt%ivel+2), prt%u(ip,prt%ivel+3)
        enddo

      if (mod((par%istep),par%expstp).eq.0) close(22) 

  end subroutine iter_part



end module part_tools
