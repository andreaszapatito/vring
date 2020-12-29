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

! Parameters
    integer      :: npart
    integer      :: nvar
    integer      :: ipos, ivel, ijac,ijrt,ihes,ihrt
    real(kind=8) :: st
  end type part

end module part_type


module part_tools
  use part_type
  use mesh_type
  use solu_type
  use comm_type
  use para_type
!  use run_tools
!  use run_type
  use para_tools

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

    real(kind=8)     :: x,y,z,rnd
    integer          :: ip,ip1

    prt%ipos=0
    prt%ivel=3
    prt%ijac=6
    prt%ijrt=15
    prt%ihes=24
    prt%ihrt=51
    
    prt%nvar=78
    prt%npart=10000
    prt%st=0.5d0

    call alct_part(prt)
    prt%u(:,:)=0.d0
    
    ip=0  
    ip1=0
    do  while (ip.le.prt%npart)
      ip1=ip1+1
      if (ip1.gt.1000000000) write (*,*) "error while count exceeds one billion"
      
      rnd=2.0*rand(52358)-1.0
      x=rnd*msh%rmax
      rnd=2.0*rand(52358)-1.0
      y=rnd*msh%rmax
      rnd=2.0*rand(52358)-1.0
      z=rnd*msh%zmax
      if (x**2+y**2.lt.msh%rmax**2) then
        ip=ip+1
        prt%u(ip,prt%ipos+1)=x
        prt%u(ip,prt%ipos+2)=y
        prt%u(ip,prt%ipos+3)=z
        
      endif
    enddo 

  end subroutine init_part

 subroutine iter_part(prt,msh,com,su)
    implicit none
    type(part)      :: prt
    type(mesh_a)    :: msh
    type(comm)      :: com
    type(solu)      :: su


    integer :: it, ir, iz
    integer :: ip,i,j,k

    real(kind=8) :: dt
    real(kind=8) :: x,y,z,r,theta,ut,ur,uz
    real(kind=8) :: tp
    real(kind=8) :: sumj,sumh
   
        
        do ip=1,prt%npart
          x=prt%u(ip,prt%ipos+1)
          y=prt%u(ip,prt%ipos+2)
          z=prt%u(ip,prt%ipos+3)
          r=dsqrt(x**2+y**2)
          if (x.ge.0.d0.and.y.ge.0.d0) then 
             theta=dacos(y/r)
          elseif (x.lt.0.d0.and.y.ge.0.d0) then  
             theta=dacos(y/r)
          elseif (x.lt.0.d0.and.y.lt.0.d0) then  
             theta=dacos(-y/r)+4.0*datan(1.d0)
          elseif (x.ge.0.d0.and.y.lt.0.d0) then  
             theta=dacos(-y/r)+4.0*datan(1.d0)
          endif  
          it=floor(theta/msh%dtheta)+1-com%ip_a(1)*msh%ntheta
          ir=floor(r/msh%dr)+1-com%ip_a(2)*msh%nr
          iz=floor(z/msh%dz)+1-com%ip_a(3)*msh%nz
          ut=su%q1(it,ir,iz)
          ur=su%q2(it,ir,iz)/msh%rm(ir)
          uz=su%q3(it,ir,iz)
          prt%vel(1)=ur*dcos(theta)-ut*dsin(theta)
          prt%vel(2)=ur*dsin(theta)+ut*dcos(theta)
          prt%vel(3)=uz
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
      enddo

  end subroutine iter_part



end module part_tools
