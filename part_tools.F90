module part_type
  implicit none
  save

  type :: part

! Data
    real(kind=8), allocatable      :: u(:,:), ddt(:,:,:), v(:,:),inflow(:,:),d(:)
    real(kind=8), allocatable      :: send_buff(:),recv_buff(:)
    real(kind=8), dimension(3)     :: vel
    real(kind=8), dimension(3,3)   :: dvel
    real(kind=8), dimension(3,3,3) :: ddvel
    
    real(kind=8), dimension(3,3)   :: du
    real(kind=8), dimension(3,3,3) :: ddu
    real(kind=8), dimension(3,3)   :: jac, jrt
    real(kind=8), dimension(3,3,3) :: hes, hrt
    real(kind=8), dimension(27,3)  :: xint
    real(kind=8), dimension(27)    :: yint
    real(kind=8), dimension(27)    :: rint
    real(kind=8), dimension(9,27)  :: mint
    real(kind=8), dimension(9)     :: acoef
    real(kind=8), dimension(9)     :: mp


! Parameters
    integer      :: id
    integer      :: npart
    integer      :: tpart
    integer      :: nfree
    integer      :: naloc
    integer      :: nlast
    integer      :: nvar
    integer      :: nder
    integer      :: ipos, ivel, ijac,ijrt,ihes,ihrt
    integer      :: jtme, jcon, jid
    integer      :: nitr
    real(kind=8) :: nprcl
    integer      :: inj
    integer      :: expstp
    integer      :: nee,nww,nss,nnn,nne,nnw,nse,nsw,ndl

    integer, allocatable   :: isfre(:),ifree(:),iee(:),iww(:),iss(:),inn(:),ine(:),inw(:),ise(:),isw(:),idl(:)
    real(kind=8) :: c
    real(kind=8) :: st
    real(kind=8) :: r
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


  subroutine alct_part(prt,msh)
    implicit none
    type(part)       :: prt
    type(mesh_a)     :: msh

    allocate(prt%u(prt%naloc,prt%nvar))
    allocate(prt%d(prt%nder))
    allocate(prt%ddt(3,prt%naloc,prt%nvar))
    allocate(prt%v(prt%naloc,3))
    allocate(prt%ifree(prt%naloc))
    allocate(prt%isfre(prt%naloc))
    allocate(prt%inflow(msh%ntheta,msh%nr))

    allocate(prt%iee(prt%naloc))
    allocate(prt%iww(prt%naloc))
    allocate(prt%iss(prt%naloc))
    allocate(prt%inn(prt%naloc))
    allocate(prt%ine(prt%naloc))
    allocate(prt%inw(prt%naloc))
    allocate(prt%ise(prt%naloc))
    allocate(prt%isw(prt%naloc))
    allocate(prt%send_buff(prt%nvar))
    allocate(prt%recv_buff(prt%nvar))

    allocate(prt%idl(prt%naloc))
  end subroutine alct_part

  subroutine inlt_part(par,prt,msh,su,com)
    implicit none
    type(para)       :: par
    type(part)       :: prt
    type(mesh_a)     :: msh
    type(solu)       :: su
    type(comm)       :: com

    real(kind=8)     :: t,r,z,x,y,inpos,u,v,w,a,ut,ur,uz

    integer          :: n,ir,it,ip,ipnew
    do ir=1,msh%nr
      if (msh%rc(ir)<par%r0.and.com%ip_a(3)==0) then
        a=msh%dtheta*(msh%rc(ir+1)**2-msh%rc(ir)**2)
        do it=1,msh%ntheta
          ut=su%q1(it,ir,1)
          ur=su%q2(it,ir,1)/msh%rm(ir)
          ut=0.d0
          ur=0.d0
          uz=su%q3(it,ir,1)
  
          prt%inflow(it,ir)=prt%inflow(it,ir)+(prt%c*a*su%q3(it,ir,1)*par%dt*float(prt%inj)/prt%nprcl)
          n=floor(prt%inflow(it,ir))
 !         n=1
        !  if (ir==8.and.it==1) write (*,*) "injection",n,it,ir,prt%inflow(it,ir),prt%c,a,su%q3(it,ir,2),prt%inj,prt%nprcl
         ! write (*,*) "injection",n,it,ir,prt%inflow(it,ir),prt%c,a,su%q3(it,ir,2),prt%inj,prt%nprcl
          prt%inflow(it,ir)=prt%inflow(it,ir)-float(n)
          
          do ip=1,n
            call random_number(inpos)
            t=msh%thc(it)+msh%dtheta*inpos
            call random_number(inpos)
            r=msh%rc(ir)+msh%dr*inpos
            call random_number(inpos)
            z=msh%zc(1)+msh%dz*inpos
            x=r*dcos(t)
            y=r*dsin(t)
            u=ur*dcos(t)-ut*dsin(t)
            v=ur*dsin(t)+ut*dcos(t)
            w=uz
!            if (prt%id==1) write (*,*) 'injecting',prt%id,prt%npart,prt%inj,n,x,y,z
            call injc_part(x,y,z,prt%c,u,v,w,prt,ipnew,1)

          enddo
        enddo
      endif
    enddo
  end subroutine inlt_part

  subroutine injc_part(x,y,z,c,u,v,w,prt,ip,icall)
    implicit none
    type(part)       :: prt
    type(mesh_a)     :: msh

    real(kind=8)     :: x,y,z,c,u,v,w
    integer          :: ip,it,icall
    it=0
    prt%npart=prt%npart+1
    if (it==1) prt%tpart=prt%npart+1
    if (prt%nfree>0) then
            ip=prt%ifree(prt%nfree)
            prt%isfre(ip)=0
            prt%ifree(prt%nfree)=0
            prt%nfree=prt%nfree-1
    else
            prt%nlast=prt%nlast+1
            ip=prt%nlast
            if (ip>prt%naloc) write (*,*) "icall", icall,x,y,z
            prt%isfre(ip)=0
    endif


!    if (ip==0) write (*,*) 'ip is zero in injc',ip,prt%nfree,prt%nlast,prt%npart
    prt%u(ip,:)=0.0
    prt%v(ip,:)=0.0

    prt%ddt(:,ip,:)=0.0

    prt%u(ip,1)=x
    prt%u(ip,2)=y
    prt%u(ip,3)=z

    prt%u(ip,4)=u
    prt%u(ip,5)=v
    prt%u(ip,6)=w

    prt%u(ip,prt%ijac+1)=1.0
    prt%u(ip,prt%ijac+5)=1.0
    prt%u(ip,prt%ijac+9)=1.0

    prt%v(ip,prt%jid)=float(prt%npart)
    prt%v(ip,prt%jtme)=0.0
    prt%v(ip,prt%jcon)=c
  end subroutine injc_part

  subroutine remv_part(ip,prt)
    implicit none
    type(part)       :: prt
    type(mesh_a)     :: msh

    real(kind=8)     :: x,y,z,inpos,dx,dy,dz
    integer          :: ip,ip1
    prt%npart=prt%npart-1
    prt%tpart=prt%tpart-1
    if (ip==prt%nlast) then
            if (prt%nfree>0) prt%ifree(prt%nfree+1)=ip
            prt%isfre(ip)=0
            prt%nlast=prt%nlast-1
    else
            prt%nfree=prt%nfree+1
            prt%ifree(prt%nfree)=ip
            prt%isfre(ip)=1
            if (ip==0) write (*,*) 'ip is zero',ip
    endif

    prt%u(ip,:)=0.0
    prt%v(ip,:)=0.0
    prt%ddt(:,ip,:)=0.0

  end subroutine remv_part


  subroutine init_part(id,prt,msh,par)
    implicit none
    save
    type(part)       :: prt
    type(para)       :: par
    type(mesh_a)     :: msh

    real(kind=8)     :: x,y,z,inpos,dx,dy,dz
    integer          :: ip,ip1,id

    prt%id=id

    prt%ipos=0
    prt%ivel=3
    prt%ijac=6
    prt%ijrt=15
    prt%ihes=24
    prt%ihrt=51
    prt%nvar=78
    prt%nder=10

    prt%npart=0
    prt%tpart=0
    prt%naloc=100000

    prt%nfree=0
    prt%nlast=0

    prt%r=0.1*(10.0**(0.5*float(id-1)))*(10.0**(-6.0))/0.02
    prt%st=(1.0/18.0)*(998.0/1.2)*(prt%r**2)*par%Re

    prt%c=1000.0*1000.0*0.02**3
    prt%nprcl=0.001
    prt%inj=10
    prt%expstp=100
    prt%jtme=1
    prt%jcon=2
    prt%jid=3
    call alct_part(prt,msh)
    prt%u(:,:)=0.d0
    prt%v(:,:)=0.d0
    prt%ddt(:,:,:)=0.d0
    prt%ifree(:)=0
    prt%isfre(:)=0
    prt%inflow(:,:)=0.d0
    
    ip=0  
    ip1=0
    dx=msh%dtheta
    dy=msh%dr
    dz=msh%dz
    prt%acoef(:)=0.d0
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
    prt%nitr=100
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
         
         
 subroutine iter_part(irk,prt,msh,com,su,par)
    implicit none
    type(part)      :: prt
    type(mesh_a)    :: msh
    type(comm)      :: com
    type(solu)      :: su
    type(para)      :: par

    logical         :: conda1,conda2,conda3,conda4
    logical         :: condb1,condb2,condb3,condb4
    logical         :: overshoot
    integer :: it, ir, iz
    integer :: ip,i,j,k,iploc
    integer :: ifst1,ifst2,ifst3
    integer :: iint
    integer :: irk
    integer :: irc,izc
    integer :: irof,itof,izof
    integer :: iv,ivar,ii
    integer :: id,id1,id2,im
    real(kind=8) :: partdt
    integer :: isubdt,nsubdt
    integer :: idir,isde,isd1,isd2,ierr,tag,nspart,nrpart,ipnew,itask
    integer, dimension(MPI_STATUS_SIZE) :: status
    real(kind=8) :: dt
    real(kind=8) :: y0(3)

    real(kind=8) :: x,y,z,r,theta,ut,ur,uz
    real(kind=8) :: tp
    real(kind=8) :: sumj,sumh
    real(kind=8) :: drp,dtp,dzp
    real(kind=8) :: fx,fy,f
    real(kind=8) :: ft,ftt,ftr,ftz
    real(kind=8) :: fr,frt,frr,frz
    real(kind=8) :: fz,fzt,fzr,fzz
    real(kind=8) :: sn,sn2,cs,cs2
    real(kind=8) :: ff(3),dff(3,3),ddff(3,3,3)
    real(kind=8) :: M(12,12),FI(12)


    character(len=8)          :: fmt2 ! format descriptor
    character(len=8)          :: fmt5 ! format descriptor
    character(2)              :: batchchar
    character(5)              :: timechar
    character(5)              :: partchar

!    if (mod((par%nstep),par%expstp).eq.0) then
!      fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
!      write (timechar,fmt5) par%nstep
!      write (partchar,fmt5) com%ip
!      do itask=0,com%np-1
!        if (com%ip==itask) then
!          open (unit=22+com%ip,file="part"//trim(timechar)//'p'//trim(partchar)//".dat", form='formatted', position='rewind')
!!          if (itask.eq.0) then 
!!            open (unit=22+com%ip,file="part"//trim(timechar)//p//trim(partchar)//".dat", form='formatted', position='rewind')
!!            write (22,*) 'title = "particles"'
!!            write (22,*) 'variables = "xp", "yp", "zp", "up", "vp", "wp" '
!!          endif
!!          if (itask.ne.0) open (unit=22,file="part"//trim(timechar)//p//trim(partchar)".dat", form='formatted', position='rewind')
!          do ip=1,prt%nlast
!            if (prt%isfre(ip)==0) then
!              write (22+com%ip,*) prt%v(ip,1),prt%v(ip,2),prt%v(ip,3),prt%u(ip,1),prt%u(ip,2),prt%u(ip,3),prt%u(ip,4),prt%u(ip,5),prt%u(ip,6) 
!            endif
!          enddo
!          close (22+com%ip)
!        endif
!        call MPI_BARRIER(com%comm_0,ierr)
!      enddo
!    endif

    if (mod((par%nstep),prt%expstp).eq.0) then
      fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
      fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left
      write (timechar,fmt5) par%nstep
      write (batchchar,fmt2) prt%id
      do itask=0,com%np-1
        if (com%ip==itask) then
          if (itask.eq.0) then 
                          open (unit=prt%id+22,file="part"//trim(timechar)//"size"//trim(batchchar)//".dat", form='formatted', position='rewind')
            write (prt%id+22,"(A32,E20.8,A12,E20.8,A2)") 'title = "particles diameter: ',prt%r,' stokes: ',prt%st,'"'
            write (prt%id+22,"(A32,100(A4,I2.2))") 'variables = "id", "time", "c",',(('u',ii),ii=1,prt%nvar),(('d',ii),ii=1,prt%nder)
          endif
          if (itask.ne.0) open (unit=prt%id+22,file="part"//trim(timechar)//"size"//trim(batchchar)//".dat", form='formatted', position='append')
          do ip=1,prt%nlast
            if (prt%isfre(ip)==0) then

              prt%d(1)=prt%u(ip,prt%ijac+1)*prt%u(ip,prt%ijac+5)*prt%u(ip,prt%ijac+9) &
                      +prt%u(ip,prt%ijac+2)*prt%u(ip,prt%ijac+6)*prt%u(ip,prt%ijac+7) &
                      +prt%u(ip,prt%ijac+3)*prt%u(ip,prt%ijac+4)*prt%u(ip,prt%ijac+8) &
                      -prt%u(ip,prt%ijac+7)*prt%u(ip,prt%ijac+5)*prt%u(ip,prt%ijac+3) &
                      -prt%u(ip,prt%ijac+4)*prt%u(ip,prt%ijac+2)*prt%u(ip,prt%ijac+9) &
                      -prt%u(ip,prt%ijac+8)*prt%u(ip,prt%ijac+6)*prt%u(ip,prt%ijac+1) 

              write (prt%id+22,"(100(E20.8))") (prt%v(ip,ii),ii=1,3),(prt%u(ip,ii),ii=1,prt%nvar),(prt%d(ii),ii=1,prt%nder)
            endif
          enddo
          close (prt%id+22)
        endif
        call MPI_BARRIER(com%comm_0,ierr)
      enddo
    endif


!      do itask=0,com%np-1
!        if (com%ip==itask) then
!          open (unit=23,file="trajectory.dat", form='formatted', position='append')
!          do ip=1,prt%npart
!            if (prt%isfre(ip)==0) then
!              write (23,*) com%ip,prt%v(ip,1),prt%v(ip,2),prt%v(ip,3),prt%u(ip,1),prt%u(ip,2),prt%u(ip,3),prt%u(ip,4),prt%u(ip,5),prt%u(ip,6)
!            endif
!          enddo
!          close (23)
!        endif
!        call MPI_BARRIER(com%comm_0,ierr)
!      enddo

if (prt%id==1) nsubdt=100
if (prt%id==2) nsubdt=10
if (prt%id==3) nsubdt=5
if (prt%id==4) nsubdt=2
if (prt%id==5) nsubdt=1
if (prt%id==6) nsubdt=1
if (prt%id==7) nsubdt=1
if (prt%id==8) nsubdt=1

  do isubdt=1,nsubdt
    partdt=par%dt/float(nsubdt)

    prt%nee=0
    prt%nww=0
    prt%nss=0
    prt%nnn=0
    prt%nne=0
    prt%nnw=0
    prt%nse=0
    prt%nsw=0
    prt%ndl=0
    do ip=1,prt%nlast
      if (prt%isfre(ip)==0) then
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

!          if (it.ge.1.and.it.le.msh%nthetag.and.ir.ge.1.and.ir.le.msh%nrg.and.iz.ge.1.and.iz.le.msh%nzg) then 
            irc=0
            izc=0
    
            if (ir.ge.msh%nr) irc=-1
            if (ir.le.1) irc=1
!            if (iz.ge.msh%nz) izc=-1
!            if (iz.le.1) izc=1
    
            do ifst1=-1,1,1
              do ifst2=-1,1,1
                do ifst3=-1,1,1
                   iv=1+(ifst1+1)+3*(ifst2+1)+9*(ifst3+1)
                   itof=it+ifst1
                   if (itof.gt.msh%ntheta) itof=1
                   if (itof.lt.1) itof=msh%ntheta
                   irof=ir+ifst2+irc
                   izof=iz+ifst3+izc
                   
                   if (ivar==1) prt%yint(iv)=su%q1(itof,irof,izof)
                   if (ivar==2.and.irof.gt.1) prt%yint(iv)=su%q2(itof,irof,izof)/msh%rc(irof)
                   if (ivar==2.and.irof.eq.1) prt%yint(iv)=0.0
                   if (ivar==3) prt%yint(iv)=su%q3(itof,irof,izof)
!                   prt%yint(iv,1)=float(izof)+float(irof)
!                   prt%yint(iv,2)=float(irof)+float(izof)
!                   prt%yint(iv,3)=1.0
!                   prt%yint(iv,1)=0.0
!                   prt%yint(iv,2)=0.1
!                   prt%yint(iv,3)=1.0
        if (isnan(prt%yint(iv))) then
                write (*,*) 'yint ivar is nan',iv,ivar,itof,irof,izof,ifst1,ifst2,ifst3
           stop 
        endif
                enddo
              enddo
            enddo

            y0(ivar)=prt%yint(14)
            do iv=1,27
              prt%yint(iv)=prt%yint(iv)-y0(ivar)
            enddo
            prt%acoef(:)=0.d0
            do iint=1,prt%nitr
              call interpolation (prt%xint,prt%acoef,prt%rint)
              do i=1,9
                prt%mp(i)=0.0
                do k=1,27
                  prt%mp(i)=prt%mp(i)-prt%mint(i,k)*(prt%rint(k)-prt%yint(k))
                enddo
              enddo
              do i=1,9
                prt%acoef(i)=prt%acoef(i)+prt%mp(i)
              enddo
!                if (ip.eq.1) write (*,*) 'particle acoef interp step',iint,prt%acoef(:),norm2(prt%mp)
              if (norm2(prt%mp).lt.0.000001) exit
            enddo



            if (dtp.gt.8.d0*datan(1.d0)) dtp=dtp-8.d0*datan(1.d0)
            if (dtp.lt.-8.d0*datan(1.d0)) dtp=dtp+8.d0*datan(1.d0)


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
            if (dtp.gt.msh%dtheta) dtp=dtp-8.d0*datan(1.d0)
            if (dtp.lt.-msh%dtheta) dtp=dtp+8.d0*datan(1.d0)
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
!            if (ip.eq.1) write (*,*) 'particle transf c ff',ivar,f,y0(ivar)
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
!            if (ip.eq.12) write (*,*) 'particle transf 3',M
!            if (ip.eq.12) write (*,*) 'particle 12 ff',f,r,theta
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
!          else
!            ff(ivar)=0.d0
!          endif
        enddo

        prt%vel(1)=ff(2)*dcos(theta)-ff(1)*dsin(theta)
        prt%vel(2)=ff(2)*dsin(theta)+ff(1)*dcos(theta)
!        prt%vel(1)=0.d0
!        prt%vel(2)=0.d0
        prt%vel(3)=ff(3)
        do id=1,3
          prt%dvel(1,id)=dff(2,id)*dcos(theta)-dff(1,id)*dsin(theta)
          prt%dvel(2,id)=dff(2,id)*dsin(theta)+dff(1,id)*dcos(theta)
          prt%dvel(3,id)=dff(3,id)
        enddo
        if (isnan(prt%vel(1))) then
           write (*,*) 'particle du,dv,dw',ip
           stop 
        endif
!        if (ip.eq.1) write (*,*) 'particle du,dv,dw', prt%dvel(1,id)
!        if (ip.eq.1) write (*,*) 'particle du,dv,dw', prt%dvel(2,id)
!        if (ip.eq.1) write (*,*) 'particle du,dv,dw', prt%dvel(3,id)
        do id1=1,3
          do id2=1,3
            prt%ddvel(1,id1,id2)=ddff(2,id1,id2)*dcos(theta)-ddff(1,id1,id2)*dsin(theta)
            prt%ddvel(2,id1,id2)=ddff(2,id1,id2)*dsin(theta)+ddff(1,id1,id2)*dcos(theta)
            prt%ddvel(3,id1,id2)=ddff(3,id1,id2)
          enddo
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
          prt%ddt(1,ip,prt%ipos+i)=prt%u(ip,prt%ivel+i)
          prt%ddt(1,ip,prt%ivel+i)=(1.0/prt%st)*(prt%vel(i)-prt%u(ip,prt%ivel+i))
        enddo
        do i=1,3
          do j=1,3
            prt%ddt(1,ip,prt%ijac+i+3*j-3)=prt%u(ip,prt%ijrt+i+3*j-3)
            sumj=prt%jac(1,j)*prt%dvel(i,1)+prt%jac(2,j)*prt%dvel(i,2)+prt%jac(3,j)*prt%dvel(i,3)
            prt%ddt(1,ip,prt%ijrt+i+3*j-3)=(1.0/prt%st)*(sumj-prt%u(ip,prt%ijrt+i+3*j-3))
          enddo
        enddo
        do i=1,3
          do j=1,3
            do k=1,3
              prt%ddt(1,ip,prt%ihes+i+3*j-3+9*k-9)=prt%u(ip,prt%ihrt+i+3*j-3+9*k-9)
              sumh=prt%hes(1,j,k)*prt%dvel(i,1)+prt%hes(2,j,k)*prt%dvel(i,2)+prt%hes(3,j,k)*prt%dvel(i,3)+prt%jac(1,j)*prt%jac(1,k)*prt%ddvel(i,1,1)+prt%jac(2,j)*prt%jac(1,k)*prt%ddvel(i,2,1)+prt%jac(3,j)*prt%jac(1,k)*prt%ddvel(i,3,1)+prt%jac(1,j)*prt%jac(2,k)*prt%ddvel(i,1,2)+prt%jac(1,j)*prt%jac(3,k)*prt%ddvel(i,1,3)+prt%jac(2,j)*prt%jac(2,k)*prt%ddvel(i,2,2)+prt%jac(2,j)*prt%jac(3,k)*prt%ddvel(i,2,3)+prt%jac(3,j)*prt%jac(2,k)*prt%ddvel(i,3,2)+prt%jac(3,j)*prt%jac(3,k)*prt%ddvel(i,3,3)

              prt%ddt(1,ip,prt%ihrt+i+3*j-3+9*k-9)=(1.0/prt%st)*(sumh-prt%u(ip,prt%ihrt+i+3*j-3+9*k-9))
            enddo
          enddo
        enddo



!        do i=1,3
!          do j=1,3
!            prt%jac(i,j)=prt%u(ip,prt%ijac+i+3*j-3)
!            prt%jrt(i,j)=prt%u(ip,prt%ijrt+i+3*j-3)
!            do k=1,3
!              prt%hes(i,j,k)=prt%u(ip,prt%ihes+i+3*j-3+9*k-9)
!              prt%hrt(i,j,k)=prt%u(ip,prt%ihrt+i+3*j-3+9*k-9)
!            enddo
!          enddo
!        enddo
!        do i=1,3
!          prt%ddt(1,ip,prt%ipos+i)=prt%u(ip,prt%ivel+i)
!          prt%ddt(1,ip,prt%ivel+i)=(1.0/prt%st)*(prt%vel(i)-prt%u(ip,prt%ivel+i))
!        enddo
!        do i=1,3
!          do j=1,3
!            prt%ddt(1,ip,prt%ijac+i+3*j-3)=prt%jrt(i,j)
!            sumj=prt%jac(1,j)*prt%du(i,1)+prt%jac(2,j)*prt%du(i,2)+prt%jac(3,j)*prt%du(i,3)
!            prt%ddt(1,ip,prt%ijrt+i+3*j-3)=(1.0/prt%st)*(sumj-prt%u(ip,prt%ijac+i+3*j-3))
!          enddo
!        enddo
!    
!        do i=1,3
!          do j=1,3
!            do k=1,3
!              prt%ddt(1,ip,prt%ihes+i+3*j-3+9*k-9)=prt%u(ip,prt%ihrt+i+3*j-3+9*k-9)
!              sumh=prt%hes(1,j,k)*prt%du(i,1)+prt%hes(2,j,k)*prt%du(i,2)+prt%hes(3,j,k)*prt%du(i,3)+prt%jac(1,j)*prt%jac(1,k)*prt%ddu(i,1,1)+prt%jac(2,j)*prt%jac(1,k)*prt%ddu(i,2,1)+prt%jac(3,j)*prt%jac(1,k)*prt%ddu(i,3,1)+prt%jac(1,j)*prt%jac(2,k)*prt%ddu(i,1,2)+prt%jac(1,j)*prt%jac(3,k)*prt%ddu(i,1,3)+prt%jac(2,j)*prt%jac(2,k)*prt%ddu(i,2,2)+prt%jac(2,j)*prt%jac(3,k)*prt%ddu(i,2,3)+prt%jac(3,j)*prt%jac(2,k)*prt%ddu(i,3,2)+prt%jac(3,j)*prt%jac(3,k)*prt%ddu(i,3,3)
!              prt%ddt(1,ip,prt%ihrt+i+3*j-3+9*k-9)=(1.0/prt%st)*(sumh-prt%u(ip,prt%ihrt+i+3*j-3+9*k-9))
!            enddo
!          enddo
!        enddo
        overshoot=.false.
!        if (sqrt(prt%ddt(1,ip,prt%ivel+1)**2+prt%ddt(1,ip,prt%ivel+2)**2+prt%ddt(1,ip,prt%ivel+3)**2)>0.1*msh%dr/partdt**2) then
!        if (sqrt(prt%ddt(1,ip,prt%ivel+1)**2+prt%ddt(1,ip,prt%ivel+2)**2+prt%ddt(1,ip,prt%ivel+3)**2)*partdt+sqrt(prt%u(ip,prt%ivel+1)**2+prt%u(ip,prt%ivel+2)**2+prt%u(ip,prt%ivel+3)**2)>0.5*msh%dr/partdt) then
!          do i=1,3
!            prt%ddt(1,ip,prt%ipos+i)=prt%vel(i)
!          enddo
!          overshoot=.true.
!        endif
!        if (overshoot) then
!          do i=1,3
!            prt%u(ip,prt%ivel+i)=prt%vel(i)
!          enddo
!        else
          do i=1,prt%nvar
            prt%u(ip,i)=prt%u(ip,i)+par%crkgam(irk)*prt%ddt(1,ip,i)*partdt+par%crkrom(irk)*prt%ddt(2,ip,i)*partdt
            prt%ddt(2,ip,i)=prt%ddt(1,ip,i)
          enddo
!        endif
        prt%v(ip,prt%jtme)=prt%v(ip,prt%jtme)+partdt
      endif
    enddo
    if (irk==3) then
    do ip=1,prt%nlast
      if (prt%isfre(ip)==0) then
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

        it=floor(theta/msh%dtheta)+1-com%ip_a(1)*msh%ntheta
        ir=floor(r/msh%dr)+1-com%ip_a(2)*msh%nr
        iz=floor(z/msh%dz)+1-com%ip_a(3)*msh%nz

        if ((z>msh%dzmax).or.(z<0.0).or.(r>msh%drmax).or.(r<0.0)) then
!          write (*,*) "Station out of bounds",ip,r,z
          prt%ndl=prt%ndl+1
          prt%idl(prt%ndl)=ip
        else
          if (iz.le.msh%nz.and.iz.ge.1) then
            if (ir.lt.1) then
              prt%nee=prt%nee+1
              prt%iee(prt%nee)=ip
            elseif (ir.gt.msh%nr) then
              prt%nww=prt%nww+1
              prt%iww(prt%nww)=ip
            endif
          elseif (ir.le.msh%nr.and.ir.ge.1) then
            if (iz.lt.1) then
              prt%nss=prt%nss+1
              prt%iss(prt%nss)=ip
            elseif (iz.gt.msh%nz) then
              prt%nnn=prt%nnn+1
              prt%inn(prt%nnn)=ip
            endif
          elseif (iz.gt.msh%nz.and.ir.gt.msh%nr) then
              prt%nnw=prt%nnw+1
              prt%inw(prt%nnw)=ip
          elseif (iz.gt.msh%nz.and.ir.lt.1) then
              prt%nne=prt%nne+1
              prt%ine(prt%nne)=ip
          elseif (iz.lt.1.and.ir.gt.msh%nr) then
              prt%nsw=prt%nsw+1
              prt%isw(prt%nsw)=ip
          elseif (iz.lt.1.and.ir.lt.1) then
              prt%nse=prt%nse+1
              prt%ise(prt%nse)=ip
          endif
        endif
      endif
    enddo
    do ip=1,prt%ndl
      call remv_part(prt%idl(ip),prt)
    enddo

!   communicate n w e s
 
    do isde=-1,1,2
      do idir=2,3
        if (idir==2.and.isde==-1) nspart=prt%nee
        if (idir==3.and.isde==-1) nspart=prt%nss
        if (idir==2.and.isde==1) nspart=prt%nww
        if (idir==3.and.isde==1) nspart=prt%nnn
        
        if ((com%ip_a(idir).lt.com%np_a(idir)-1.and.isde==1).or.(com%ip_a(idir).gt.0.and.isde==-1)) then
          tag = 1
          call mpi_send(nspart,1,MPI_INTEGER,com%ip_a(idir)+isde,tag,com%comm_a(idir),ierr)
          do ip=1,nspart
           if (idir==2.and.isde==-1) iploc=prt%iee(ip)
           if (idir==3.and.isde==-1) iploc=prt%iss(ip)
           if (idir==2.and.isde==1) iploc=prt%iww(ip)
           if (idir==3.and.isde==1) iploc=prt%inn(ip)
            do iv=1,prt%nvar
              prt%send_buff(iv)=prt%u(iploc,iv)
            enddo
            tag = 1+ip
            call mpi_send(prt%send_buff,prt%nvar,MPI_DOUBLE_PRECISION,com%ip_a(idir)+isde,tag,com%comm_a(idir),ierr)
            do iv=1,3
              prt%send_buff(iv)=prt%v(iploc,iv)
            enddo
            tag = 1+nspart+ip
            call mpi_send(prt%send_buff(1),3,MPI_DOUBLE_PRECISION,com%ip_a(idir)+isde,tag,com%comm_a(idir),ierr)
            call remv_part(iploc,prt)
          enddo
        endif
    
        if ((com%ip_a(idir).lt.com%np_a(idir)-1.and.isde==-1).or.(com%ip_a(idir).gt.0.and.isde==1)) then
          tag = 1
          call mpi_recv(nrpart,1,MPI_INTEGER,com%ip_a(idir)-isde,tag,com%comm_a(idir),status,ierr)
          do ip=1,nrpart
            tag = 1+ip
            call mpi_recv(prt%recv_buff,prt%nvar,MPI_DOUBLE_PRECISION,com%ip_a(idir)-isde,tag,com%comm_a(idir),status,ierr)
            call injc_part(prt%recv_buff(1),prt%recv_buff(2),prt%recv_buff(3),prt%c,prt%recv_buff(4),prt%recv_buff(5),prt%recv_buff(6),prt,ipnew,2)
            do iv=1,prt%nvar
              prt%u(ipnew,iv)=prt%recv_buff(iv)
            enddo
            tag = 1+nrpart+ip
            call mpi_recv(prt%recv_buff(1),3,MPI_DOUBLE_PRECISION,com%ip_a(idir)-isde,tag,com%comm_a(idir),status,ierr)
            do iv=1,3
              prt%v(ipnew,iv)=prt%recv_buff(iv)
            enddo
          enddo
        endif
        call MPI_BARRIER(com%comm_0,ierr)
      enddo
    enddo
        call MPI_BARRIER(com%comm_0,ierr)
    write (*,*) "A npart:",prt%npart," nlast:",prt%nlast," nfree",prt%nfree," ip ",com%ip
        call MPI_BARRIER(com%comm_0,ierr)

    do isd2=-1,1,2
      do isd1=-1,1,2
        conda1=((isd1==-1).and.(isd2==-1))
        conda2=((isd1==-1).and.(isd2==1))
        conda3=((isd1==1).and.(isd2==-1))
        conda4=((isd1==1).and.(isd2==1))

        condb1=((com%ip_a(2).lt.com%np_a(2)-1).and.(com%ip_a(3).lt.com%np_a(3)-1))
        condb2=((com%ip_a(2).lt.com%np_a(2)-1).and.(com%ip_a(3).gt.0))
        condb3=((com%ip_a(2).gt.0).and.(com%ip_a(3).lt.com%np_a(3)-1))
        condb4=((com%ip_a(2).gt.0).and.(com%ip_a(3).gt.0))
        if (isd1==-1.and.isd2==-1) nspart=prt%nse
        if (isd1==-1.and.isd2==1) nspart=prt%nsw
        if (isd1==1.and.isd2==-1) nspart=prt%nne
        if (isd1==1.and.isd2==1) nspart=prt%nnw
        if  ((conda1.and.condb4).or.(conda2.and.condb2).or.(conda3.and.condb3).or.(conda4.and.condb1)) then
          tag = 1
          call mpi_send(nspart,1,MPI_INTEGER,com%ip+isd2*com%np_a(3)+isd1,tag,com%comm_0,ierr)
          do ip=1,nspart
            if (isd1==-1.and.isd2==-1) iploc=prt%ise(ip)
            if (isd1==-1.and.isd2==1) iploc=prt%isw(ip)
            if (isd1==1.and.isd2==-1) iploc=prt%ine(ip)
            if (isd1==1.and.isd2==1) iploc=prt%inw(ip)
            prt%send_buff=prt%u(iploc,:)
            tag = 1+ip
            call mpi_send(prt%send_buff,prt%nvar,MPI_DOUBLE_PRECISION,com%ip+isd2*com%np_a(3)+isd1,tag,com%comm_0,ierr)
            prt%send_buff(1:3)=prt%v(iploc,1:3)
            tag = 1+nspart+ip
            call mpi_send(prt%send_buff,3,MPI_DOUBLE_PRECISION,com%ip+isd2*com%np_a(3)+isd1,tag,com%comm_0,ierr)

            call remv_part(iploc,prt)
          enddo
        endif

        if  ((conda1.and.condb1).or.(conda2.and.condb3).or.(conda3.and.condb2).or.(conda4.and.condb4)) then
          tag = 1
          call mpi_recv(nrpart,1,MPI_INTEGER,com%ip-isd2*com%np_a(3)-isd1,tag,com%comm_0,status,ierr)
          do ip=1,nrpart
            tag = 1+ip
            call mpi_recv(prt%recv_buff,prt%nvar,MPI_DOUBLE_PRECISION,com%ip-isd2*com%np_a(3)-isd1,tag,com%comm_0,status,ierr)
            call injc_part(0.d0,0.d0,0.d0,prt%c,0.d0,0.d0,0.d0,prt,ipnew,3)
            prt%u(ipnew,:)=prt%recv_buff
            tag = 1+ip+nrpart
            call mpi_recv(prt%recv_buff,3,MPI_DOUBLE_PRECISION,com%ip-isd2*com%np_a(3)-isd1,tag,com%comm_0,status,ierr)
            prt%v(ipnew,:)=prt%recv_buff(1:3)
          enddo
        endif
        call MPI_BARRIER(com%comm_0,ierr)
      enddo
    enddo  
  endif
  enddo
end subroutine iter_part

end module part_tools
