module run_type
!#include <petsc/finclude/petscdef.h>
!  use petsc
  use comm_type 
  use cons_type
  use para_type
  use mesh_type
  use solu_type
  use solp_type
  use part_type

  implicit none
  save
  type :: run
    type (comm)  :: com
    type (cons)  :: con
    type (para)  :: par
    type (mesh_a):: msh
    type (solu)  :: su
    type (solp)  :: sp
    type (part)  :: prt
  end type run
end module run_type

module run_tools
  use run_type
!#include <mpif.h>  
  contains
    subroutine echo_var(var)
    implicit none
    type (run) :: var


    if (var%com%ip==0) then

      open(unit=22, file=trim(adjustl(var%con%out_file)), form='formatted', status='replace')
  
      write(22,"(A,3(I16.8))") 'Domain decomposition__(comm)_'
      write(22,"(A,3(I16.8))") 'n. processors (th,r,z):', var%com%np_a

      write(22,"(A,3(I16.8))") 'Mesh characteristics__(mesh)_'
      write(22,"(A,3(I16.8))") 'nth,nr,nz   local     :', var%msh%ntheta , var%msh%nr , var%msh%nz
      write(22,"(A,3(I16.8))") 'nth,nr,nz   global    :', var%msh%nthetag, var%msh%nrg, var%msh%nzg
      write(22,"(A,3(F16.8))") 'dth,dr,dz             :', var%msh%dtheta, var%msh%dr, var%msh%dz
      write(22,"(A,3(F16.8))") 'rmax                  :', var%msh%rmax, var%msh%rm(1), var%msh%rm(var%msh%nr)
      write(22,"(A,3(F16.8))") 'zmax                  :', var%msh%zmax, var%msh%zm(1), var%msh%zm(var%msh%nz)
      write(22,"(A,3(I16.8))") 'Problem constants_____(cons)_'
      write(22,"(A,3(I16.8))") 'solver                :'

      write(22,"(A,3(F16.8))") 'Problem parameters____(para)_'
      write(22,"(A,3(I16.8))") 'total steps           :', var%par%totstp
      write(22,"(A,3(I16.8))") 'export every          :', var%par%expstp
      write(22,"(A,3(I16.8))") 'analys every          :', var%par%anlstp
      write(22,"(A,3(I16.8))") 'restart               :', var%par%restrt
      write(22,"(A,3(F16.8))") 'u1,v2,w3    nominal   :', var%par%u1,var%par%u2,var%par%u3
      write(22,"(A,4(F16.8))") 'ur,wr,uf,kf nominal   :', var%par%uring,var%par%uswrl,var%par%ufluc,var%par%kfluc
      write(22,"(A,4(F16.8))") 'rw,r0,t0              :', var%par%rw,var%par%r0,var%par%t0
      write(22,"(A,3(F16.8))") 'p           nominal   :', var%par%p                       
      write(22,"(A,3(F16.8))") 'Re                    :', var%par%Re
      write(22,"(A,3(F16.8))") 'Pe                    :', var%par%Pe
      write(22,"(A,3(F16.8))") 'dt                    :', var%par%dt
      write(22,"(A,3(I16.8))") 'rkstep                :', var%par%nrkstep
      write(22,"(A,3(F16.8))") 'slip                  :', var%par%slip
 

 
      close(22)
    endif
    end subroutine echo_var

 subroutine ffemexport (var,id)
  implicit none

  type(run)                 :: var
  integer                   :: i,ilocal
  integer                   :: j,jlocal
  integer                   :: k,klocal
  integer                   :: ip,jp,id
  integer                   :: nvert,nedge,ntrig
  integer                   :: ipart,jpart,tpart

  real(kind=8)              :: qan,dwdr,dvdz

  character(len=8)          :: fmt5,fmt2 ! format descriptor
  character(5)              :: timechar
  character(2)              :: rkchar
  character(2)              :: prchar
  integer                   :: tag,ierr
  integer, dimension(MPI_STATUS_SIZE) :: status
  fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
  fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left


  write (timechar,fmt5) var%par%nstep
  write (rkchar,fmt2) id
  write (prchar,fmt2) var%com%ip

  if (var%com%ip.eq.0.and.id.eq.0) then
    open(unit=21,file="FFEM"//trim(adjustl(var%con%fld_file))//trim(timechar)//"ID"//trim(rkchar)//".msh", form='formatted', position='rewind')

    nvert=(var%msh%nrg+1)*(var%msh%nzg+1)
    nedge=2*var%msh%nrg+2*var%msh%nzg
    ntrig=2*var%msh%nrg*var%msh%nzg


    write (21,"(3(I10))") nvert, ntrig, nedge
    ! verticies
    do j=1,var%msh%nzg+1
      do i=1,var%msh%nrg+1
        write (21,"(2(E16.8),I5)") (i-1)*var%msh%dr,(j-1)*var%msh%dz, 0
      enddo
    enddo
    ! triangles
    do j=1,var%msh%nzg
      do i=1,var%msh%nrg
        write (21,"(4(I10))") (var%msh%nrg+1)*(j-1)+i,(var%msh%nrg+1)*(j-1)+i+1,(var%msh%nrg+1)*j+i+1, 10
        write (21,"(4(I10))") (var%msh%nrg+1)*(j-1)+i,(var%msh%nrg+1)*(j)+i+1,  (var%msh%nrg+1)*j +i, 10
      enddo
    enddo
    !edges
    do i=1,var%msh%nrg
      write (21,"(4(I10))") i,i+1,1
    enddo

    do j=1,var%msh%nzg
      write (21,"(4(I10))") (var%msh%nrg+1)*j,(var%msh%nrg+1)*(j+1),2
    enddo

    do i=var%msh%nrg+1,2,-1
      write (21,"(4(I10))") (var%msh%nrg+1)*var%msh%nzg+i,(var%msh%nrg+1)*var%msh%nzg+i-1,3
    enddo

    do j=var%msh%nzg,1,-1
      write (21,"(4(I10))") (var%msh%nrg+1)*j+1,(var%msh%nrg+1)*(j-1)+1,4
    enddo



    close (21)
  endif



  if (id.ne.0) then
    if (var%com%ip.eq.0) then
      open(unit=21,file="FFEM"//trim(adjustl(var%con%fld_file))//trim(timechar)//"ID"//trim(rkchar)//".dat", form='formatted', position='rewind')
      write (21,"(3(I10))") (var%msh%nrg+1)*(var%msh%nzg+1)
    endif

    ! verticies
      
    do klocal=1,var%msh%nz
      do jlocal=1,var%msh%nr
        do ilocal=1,var%msh%ntheta
           qan=0.0
           i=ilocal
           j=jlocal+var%com%ip_a(2)*var%msh%nr
           k=klocal+var%com%ip_a(3)*var%msh%nz
           if (id.eq.3) then
               if (j.eq.1) qan=var%su%q3(ilocal,jlocal,klocal)
               if (j.gt.1) qan=0.5*(var%su%q3(ilocal,jlocal,klocal)+var%su%q3(ilocal,jlocal-1,klocal))
           elseif (id.eq.2) then
               if (k.eq.1.and.j.gt.1) qan=var%su%q2(ilocal,jlocal,klocal)/var%msh%rc(jlocal)
               if (k.gt.1.and.j.gt.1) qan=0.5*(var%su%q2(ilocal,jlocal,klocal)/var%msh%rc(jlocal)+var%su%q2(ilocal,jlocal,klocal-1)/var%msh%rc(jlocal))
               if (j.eq.1) qan=0.0
           elseif (id.eq.1) then
               if (k.eq.1.and.j.eq.1) qan=var%su%q3(ilocal,jlocal,klocal)
               if (k.gt.1.and.j.eq.1) qan=0.5*(var%su%q1(ilocal,jlocal,klocal)+var%su%q1(ilocal,jlocal,klocal-1))
               if (k.eq.1.and.j.gt.1) qan=0.5*(var%su%q1(ilocal,jlocal,klocal)+var%su%q1(ilocal,jlocal-1,klocal))
               if (k.gt.1.and.j.gt.1) qan=0.25*(var%su%q1(ilocal,jlocal,klocal)+var%su%q1(ilocal,jlocal-1,klocal)+var%su%q1(ilocal,jlocal,klocal-1)+var%su%q1(ilocal,jlocal-1,klocal-1))
           elseif (id.eq.4) then
               if (k.eq.1.and.j.eq.1) qan=var%su%pres(ilocal,jlocal,klocal)
               if (k.gt.1.and.j.eq.1) qan=0.5*(var%su%pres(ilocal,jlocal,klocal)+var%su%pres(ilocal,jlocal,klocal-1))
               if (k.eq.1.and.j.gt.1) qan=0.5*(var%su%pres(ilocal,jlocal,klocal)+var%su%pres(ilocal,jlocal-1,klocal))
               if (k.gt.1.and.j.gt.1) qan=0.25*(var%su%pres(ilocal,jlocal,klocal)+var%su%pres(ilocal,jlocal-1,klocal)+var%su%pres(ilocal,jlocal,klocal-1)+var%su%pres(ilocal,jlocal-1,klocal-1))
           elseif (id.eq.5) then
             if (j.eq.1.or.j.eq.var%msh%nrg+1.or.k.eq.1.or.k.eq.var%msh%nzg+1) then
               dwdr=0.0
               dvdz=0.0
             else
               dwdr=(var%su%q3(i,jlocal,klocal)-var%su%q3(i,jlocal-1,klocal))/(var%msh%rm(jlocal)-var%msh%rm(jlocal-1))
               dvdz=(var%su%q2(i,jlocal,klocal)/var%msh%rc(jlocal)-var%su%q2(i,jlocal,klocal-1)/var%msh%rc(jlocal))/(var%msh%zm(klocal)-var%msh%zm(klocal-1))
             endif
             qan=-dwdr+dvdz
           endif
           var%su%ql(jlocal,klocal)=qan
        enddo
      enddo
    enddo

    do ip=0,var%com%np_a(2)-1
      do jp=0,var%com%np_a(3)-1
        ipart=jp+ip*var%com%np_a(2)
        tag=ipart
!        write(*,*) var%com%ip_a(2), var%com%ip_a(3), var%com%ip
        call mpi_barrier(var%com%comm_0,ierr)
        
        if(ip.ne.0.or.jp.ne.0) then
          if (var%com%ip_a(2).eq.ip.and.var%com%ip_a(3).eq.jp) call mpi_send(var%su%ql,var%msh%nr*var%msh%nz,MPI_DOUBLE_PRECISION,0,    tag,var%com%comm_0,ierr)
          if (var%com%ip_a(2).eq.0.and.var%com%ip_a(3).eq.0) then
                                                               call mpi_recv(var%su%ql,var%msh%nr*var%msh%nz,MPI_DOUBLE_PRECISION,ipart,tag,var%com%comm_0,status,ierr)
            var%su%qg(ip*var%msh%nr+1:(ip+1)*var%msh%nr,jp*var%msh%nz+1:(jp+1)*var%msh%nz)=var%su%ql
          endif
        else
          if (var%com%ip_a(2).eq.0.and.var%com%ip_a(3).eq.0) then 
           var%su%qg(ip*var%msh%nr+1:(ip+1)*var%msh%nr,jp*var%msh%nz+1:(jp+1)*var%msh%nz)=var%su%ql
          endif
        endif
       call mpi_barrier(var%com%comm_0, ierr)
      enddo
    enddo

    if (var%com%ip.eq.0) then
    do j=1,var%msh%nzg+1
      do i=1,var%msh%nrg+1
        if (i.eq.var%msh%nrg+1.or.j.eq.var%msh%nzg+1) then
          write (21,"(2(E16.8),I5)") 0.0
        else
          write (21,"(2(E16.8),I5)") var%su%qg(i,j)
        endif
      enddo
    enddo
    
      close (21)
    endif
  endif
!    a single line giving NV, NE, NT, the number of vertices, boundary edges, and triangles
!    NV lines, each containing the (X,Y) coordinates and integer label for a node.
!    NT lines, each containing three vertex indices and integer label for a triangle.
!    NE lines, each containing two vertex indices and integer label for a boundary edge.





  end subroutine ffemexport

  subroutine analysis (var)
  implicit none

  type(run)                 :: var
  integer                   :: i
  integer                   :: iz
  integer                   :: ir
  integer                   :: ipass
  integer                   :: ipr,ipz
  integer                   :: izl,irl
  integer                   :: ith
  integer                   :: ip
  real(kind=8)              :: xc,yc,ux,uy,l2(3),l2all(3),q3sol,q2sol,q1sol,q1m,q2m,q3m
  real(kind=8)              :: vortmax,swrlmax,dwdr,dvdz,vort,swrl,xvort,yvort,xswrl,yswrl,vortmaxall,swrlmaxall,vvort,vswrl,uvort,uswrl,xvortold,xswrlold,xvortnew,xswrlnew
  real(kind=8),dimension(20):: statall
  integer                   :: ierr

  character(len=8)          :: fmt5,fmt2 ! format descriptor
  character(5)              :: timechar
  character(2)              :: rkchar
  character(2)              :: prchar
  fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
  fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left

  
  write (timechar,fmt5) var%par%nstep
  write (rkchar,fmt2) var%par%irkstep
  write (prchar,fmt2) var%com%ip

  do ipass=1,2
    do ip=0,var%com%np-1
      call mpi_barrier(var%com%comm_0, ierr)
      call mpi_barrier(var%com%comm_0, ierr)
      if (ip.ne.-1.and.ip.eq.var%com%ip) then
        if (mod((var%par%nstep),var%par%expstp).eq.0) then
          if (ip.eq.0) then
            open(unit=21,file="2D"//trim(adjustl(var%con%fld_file))//trim(timechar)//"RK"//trim(rkchar)//".dat", form='formatted', position='rewind')
            open(unit=22,file="3D"//trim(adjustl(var%con%fld_file))//trim(timechar)//"RK"//trim(rkchar)//"PR"//trim(prchar)//".dat", form='formatted', position='rewind')
          else
            open(unit=21,file="2D"//trim(adjustl(var%con%fld_file))//trim(timechar)//"RK"//trim(rkchar)//".dat", form='formatted', position='append')
            open(unit=22,file="3D"//trim(adjustl(var%con%fld_file))//trim(timechar)//"RK"//trim(rkchar)//"PR"//trim(prchar)//".dat", form='formatted', position='rewind')
          endif
    
!          write (21,"(a)") 'title = "sample mesh"'
!          write (21,"(a)") 'variables = "x", "y", "z", "u", "v", "w", "p", "ut", "ur", "vr", "sr0", "vr0" '
!          write (21,"(a,I5,a,I5,a,I5,a)") 'zone i=',1             ,', j=',var%msh%nr,',k=',var%msh%nz,' f=point'
    
          write (22,"(a)") 'title = "sample mesh"'
          write (22,"(a)") 'variables = "x", "y", "z", "u", "v", "w", "p", "ut", "ur", "vr", "sr0", "vr0" '
          write (22,"(a,I5,a,I5,a,I5,a)") 'zone i=',var%msh%ntheta,', j=',var%msh%nr,',k=',var%msh%nz,' f=point'
        endif
        
        vortmax=-10000.0
        swrlmax=-10000.0
        xvort=0.0
        yvort=0.0
        vvort=0.0
        uvort=0.0

        xswrl=0.0
        yswrl=0.0
        vswrl=0.0
        uswrl=0.0

        do iz=1,var%msh%nz
          do ir=1,var%msh%nr
            do ith=1,var%msh%ntheta
              if (ir.eq.1.and.var%com%ip_a(2).eq.0) var%su%q3(ith,0,iz)=var%su%q3(ith,1,iz)
              if (ir.eq.1.and.var%com%ip_a(2).eq.0) var%su%q2(ith,0,iz)=var%su%q2(ith,1,iz)
              if (ir.eq.1.and.var%com%ip_a(2).eq.0) var%su%q1(ith,0,iz)=var%su%q1(ith,1,iz)
!!            if (ir.eq.1) var%su%q2(ith,0,iz)=var%su%q2(ith,1,iz)
!!            if (ir.eq.1) var%su%q1(ith,0,iz)=var%su%q1(ith,1,iz)
!!            if (ir.eq.1.and.iz.eq.var%msh%nz) var%su%q3(ith,0,iz+1)=var%su%q3(ith,1,iz+1)
              xc=var%msh%rm(ir)*cos(var%msh%thm(ith))
              yc=var%msh%rm(ir)*sin(var%msh%thm(ith))
              if (var%msh%ntheta.gt.1) q1m=0.5*(var%su%q1(ith,ir,iz)+var%su%q1(ith+1,ir,iz))
              if (var%msh%ntheta.eq.1) q1m=1.0*(var%su%q1(ith,ir,iz))
              q2m=0.5*(var%su%q2(ith,ir,iz)+var%su%q2(ith,ir+1,iz))
              q3m=0.5*(var%su%q3(ith,ir,iz)+var%su%q3(ith,ir,iz+1))
              
              dwdr=0.5*((var%su%q3(ith,ir+1,iz)+var%su%q3(ith,ir+1,iz+1))-(var%su%q3(ith,ir-1,iz)+var%su%q3(ith,ir-1,iz+1)))/(var%msh%rm(ir+1)-var%msh%rm(ir-1))
              dvdz=0.5*((var%su%q2(ith,ir,iz+1)/var%msh%rm(ir)+var%su%q2(ith,ir+1,iz+1)/var%msh%rm(ir+1))-(var%su%q2(ith,ir+1,iz-1)/var%msh%rm(ir+1)+var%su%q2(ith,ir,iz-1)/var%msh%rm(ir)))/(var%msh%zm(iz+1)-var%msh%zm(iz-1))

              vort=+(dvdz-dwdr)
              swrl=q1m

              ux=q2m*cos(var%msh%thm(ith))/var%msh%rm(ir)-q1m*sin(var%msh%thm(ith))
              uy=q2m*sin(var%msh%thm(ith))/var%msh%rm(ir)+q1m*cos(var%msh%thm(ith))
              if ((var%msh%zm(iz).lt.0.90*var%msh%zmaxg).and.(var%msh%zm(iz).gt.0.10*var%msh%zmaxg).and.(var%msh%rm(ir).lt.0.80*var%msh%rmaxg).and.(var%msh%rm(ir).gt.var%msh%rmaxg*0.05).and.(ir.ne.1.or.var%com%ip_a(2).ne.0)) then
              if (ipass.eq.1) vortmax=max(vort,vortmax)             
              if (ipass.eq.1) swrlmax=max(swrl,swrlmax)             
              if (ipass.eq.2) then
                if (vort.gt.0.00*vortmaxall) then 
                  xvort=xvort+var%msh%zm(iz)*vort**2*(var%msh%rm(ir)**1)*var%msh%dx1*var%msh%dx2*var%msh%dx3
                  vvort=vvort+               vort**2*(var%msh%rm(ir)**1)*var%msh%dx1*var%msh%dx2*var%msh%dx3
                  yvort=yvort+               vort**2*(var%msh%rm(ir)**1)*var%msh%dx1*var%msh%dx2*var%msh%dx3
                  uvort=uvort+               vort**2                    *var%msh%dx1*var%msh%dx2*var%msh%dx3
                endif

                if (swrl.gt.0.0*swrlmaxall) then 
                  xswrl=xswrl+var%msh%zm(iz)*swrl**2*(var%msh%rm(ir)**1)*var%msh%dx1*var%msh%dx2*var%msh%dx3
                  vswrl=vswrl+               swrl**2*(var%msh%rm(ir)**1)*var%msh%dx1*var%msh%dx2*var%msh%dx3
                  yswrl=yswrl+               swrl**2*(var%msh%rm(ir)**1)*var%msh%dx1*var%msh%dx2*var%msh%dx3
                  uswrl=uswrl+               swrl**2                    *var%msh%dx1*var%msh%dx2*var%msh%dx3
                endif
              endif
              endif
  
              !if (mod((var%par%nstep-1),var%par%expstp).eq.0.and.ipass.eq.2) write (22,"(20(E16.8))") xc,yc,var%msh%zm(iz),ux,uy,var%su%q3(ith,ir,iz),var%su%pres(ith,ir,iz),q1m,q2m/var%msh%rm(ir),vort,q1m/(swrlmaxall+0.000001),vort/(vortmaxall+0.000001)
              if (mod((var%par%nstep),var%par%expstp).eq.0.and.ipass.eq.2) write (22,"(20(E16.8))") xc,yc,var%msh%zm(iz),ux,uy,var%su%q3(ith,ir,iz),var%su%pres(ith,ir,iz),q1m,q2m/var%msh%rm(ir),vort,q1m/(swrlmaxall+0.000001),vort/(vortmaxall+0.000001)
              if (mod((var%par%nstep),var%par%expstp).eq.0.and.ipass.eq.2.and.ith.eq.1) write (21,"(20(E16.8))") xc,yc,var%msh%zm(iz),ux,uy,var%su%q3(ith,ir,iz),var%su%pres(ith,ir,iz),q1m,q2m/var%msh%rm(ir),vort,q1m/(swrlmaxall+0.000001),vort/(vortmaxall+0.000001)
            enddo
 
          enddo
          if (mod((var%par%nstep),var%par%expstp).eq.0.and.ipass.eq.2) write (21,"(a)")                           
        enddo
        if (mod((var%par%nstep),var%par%expstp).eq.0.and.ipass.eq.2) write (21,"(a)")                           
          
!        if (ipass.eq.2) write(*,*) vortmaxall,swrlmaxall,vort,uvort,vvort,xvort,uvort
      if (mod((var%par%nstep),var%par%expstp).eq.0)  close(22)
      if (mod((var%par%nstep),var%par%expstp).eq.0)  close(21)
         
    endif
    call mpi_barrier(var%com%comm_0, ierr)
    call mpi_barrier(var%com%comm_0, ierr)
  enddo
  if (ipass.eq.1) call mpi_allreduce(vortmax,vortmaxall,1,mpi_double_precision,mpi_max,var%com%comm_0,ierr)
  if (ipass.eq.1) call mpi_allreduce(swrlmax,swrlmaxall,1,mpi_double_precision,mpi_max,var%com%comm_0,ierr)
  if (ipass.eq.2) then
    xvortold=var%par%stat(1)/var%par%stat(3)
    xswrlold=var%par%stat(5)/var%par%stat(7)
    var%par%stat(1:8)=(/xvort,yvort,vvort,uvort,xswrl,yswrl,vswrl,uswrl/)
    call mpi_allreduce(var%par%stat,statall,20,mpi_double_precision,mpi_sum,var%com%comm_0,ierr)
    var%par%stat=statall
    xvortnew=var%par%stat(1)/var%par%stat(3)
    xswrlnew=var%par%stat(5)/var%par%stat(7)
    var%par%stat(9)=(xvortnew-xvortold)/(var%par%dt*var%par%anlstp)
    var%par%stat(10)=(xswrlnew-xswrlold)/(var%par%dt*var%par%anlstp)
  endif
 
  enddo

  if (var%com%ip.eq.0) then
    if (var%par%nstep.eq.1) open(unit=23,file="analysis.dat", form='formatted', position='rewind')
    if (var%par%nstep.ne.1) open(unit=23,file="analysis.dat", form='formatted', position='append')
    write (23,"(I5,20(E16.8))") var%par%nstep,var%par%ntime,vortmaxall,swrlmaxall,(var%par%stat(i)**(1.0/float(i))/var%par%stat(i+2)**(1.0/float(i)),i=1,2),(var%par%stat(i+2),i=1,2),((var%par%stat(i+4)**(1.0/float(i)))/var%par%stat(i+6)**(1.0/float(i)),i=1,2),(var%par%stat(i+6)**(1.0/float(i)),i=1,2),var%par%stat(9),var%par%stat(10)
    close(23)
  endif
  end subroutine analysis

  subroutine restart_save (var)
  implicit none

  type(run)                :: var
  integer                  :: iz
  integer                  :: ir
  integer                  :: ipr,ipz
  integer                  :: izl,irl
  integer                  :: ith
  integer                  :: ip
  real(kind=8)             :: xc,yc,ux,uy
  integer                  :: ierr

  character(len=8)         :: fmt5,fmt2 ! format descriptor
  character(5)             :: timechar
  character(2)             :: rkchar
  fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
  fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left

  

  write (timechar,fmt5) var%par%nstep
  write (rkchar,fmt2) var%par%irkstep

!  do ip=0,var%com%np-1
!    if (ip.eq.var%com%ip) then



  if (var%com%ip.eq.0) then
    open(unit=22,file="restart"//trim(timechar)//"RK"//trim(rkchar)//".dat", form='unformatted', position='rewind')
!    write (22,"(I5,E25.16,I5,I5,I5)") var%par%nstep,var%par%ntime,var%msh%nthetag,var%msh%nrg,var%msh%nzg
    write (22) var%par%nstep,var%par%ntime,var%msh%nthetag,var%msh%nrg,var%msh%nzg
    close(22)
  endif


  do iz=1,var%msh%nzg
    do ir=1,var%msh%nrg
      ipz=floor(float(iz-1)/float(var%msh%nz))
      ipr=floor(float(ir-1)/float(var%msh%nr))
!      izl=mod(iz-1,var%msh%nz)+1
!      irl=mod(ir-1,var%msh%nr)+1
      irl=ir-ipr*var%msh%nr
      izl=iz-ipz*var%msh%nz
!      if (iz.eq.var%msh%nzg+1) then
!        ipz=var%com%ip_a(3)-1
!        izl=var%msh%nz+1
!      endif
!      if (ir.eq.var%msh%nrg+1) then
!        ipr=var%com%ip_a(2)-1
!        irl=var%msh%nr+1
!      endif

      if ((var%com%ip_a(2).eq.ipr).and.(var%com%ip_a(3).eq.ipz)) then
        if (irl.eq.1) then
          open(unit=22,file="restart"//trim(timechar)//"RK"//trim(rkchar)//".dat", form='unformatted', position='append')
        endif 

        do ith=1,var%msh%ntheta
!          write (22,"(10(E25.16))") var%su%q1(ith,ir,iz),var%su%q2(ith,ir,iz),var%su%q3(ith,ir,iz)
          write (22) var%su%q1(ith,irl,izl),var%su%q2(ith,irl,izl),var%su%q3(ith,irl,izl)
        enddo
        if (irl.eq.var%msh%nr) then
          close(22)
        endif 
      endif 
      call mpi_barrier(var%com%comm_0, ierr)
    enddo
  enddo
!      endif 
!      call mpi_barrier(var%com%comm_0, ierr)
!      enddo



  end subroutine restart_save


  subroutine restart_load (var)
  implicit none

  type(run)                :: var
  integer                  :: iz
  integer                  :: ir
  integer                  :: ipr,ipz
  integer                  :: izl,irl
  integer                  :: ith
  integer                  :: ip
  real(kind=8)             :: xc,yc,ux,uy,buf
  integer                  :: ierr

  character(len=8)         :: fmt5,fmt2 ! format descriptor
  character(5)             :: timechar
  character(2)             :: rkchar
  fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
  fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left

  

 write (timechar,fmt5) var%par%nstep
 write (rkchar,fmt2) var%par%irkstep
 do ip=0,var%com%np-1
   if (ip.eq.var%com%ip) then
     open(unit=22,file="RESTART" , form='unformatted', position='rewind')
!        read (22,"(I5,E25.16,I5,I5,I5)") var%par%nstepi,var%par%ntime,var%msh%nthetag,var%msh%nrg,var%msh%nzg
        read (22) var%par%nstepi,var%par%ntime,var%msh%nthetag,var%msh%nrg,var%msh%nzg
         !  write (*,*) "open",var%com%ip_a(2),var%com%ip_a(3),ir,iz


      do iz=1,var%msh%nzg
        do ir=1,var%msh%nrg
          ipz=floor(float(iz-1)/float(var%msh%nz))
          ipr=floor(float(ir-1)/float(var%msh%nr))
!          izl=mod(iz-1,var%msh%nz)+1
!          irl=mod(ir-1,var%msh%nr)+1
          irl=ir-ipr*var%msh%nr
          izl=iz-ipz*var%msh%nz

!          if (iz.eq.var%msh%nzg+1) then
!            ipz=var%com%ip_a(3)-1
!            izl=var%msh%nz+1
!          endif
!          if (ir.eq.var%msh%nrg+1) then
!            ipr=var%com%ip_a(2)-1
!            irl=var%msh%nr+1
!          endif

          if ((var%com%ip_a(2).eq.ipr).and.(var%com%ip_a(3).eq.ipz)) then
            do ith=1,var%msh%ntheta
!              read  (22,"(10(E25.16))") var%su%q1(ith,ir,iz),var%su%q2(ith,ir,iz),var%su%q3(ith,ir,iz)
              read  (22) var%su%q1(ith,irl,izl),var%su%q2(ith,irl,izl),var%su%q3(ith,irl,izl)
            enddo
!            write (*,*) "part",var%com%ip_a(2),var%com%ip_a(3),ir,iz
          else
            do ith=1,var%msh%ntheta
!            read  (22,"(10(E25.16))") buf,buf,buf
              read  (22) buf,buf,buf
            enddo
          endif
        enddo
      enddo
      close(22)
    endif
  call mpi_barrier(var%com%comm_0, ierr)
 enddo



  end subroutine restart_load



  subroutine verbose(var,message)
  type(run)                :: var
  character (len=*)        :: message

  integer                  :: ierr

    if ((var%com%ip.eq.0).and.(var%par%verbose.eq.1)) write (*,"(100a)") ":::",message


    call mpi_barrier(var%com%comm_0, ierr)
  end subroutine verbose


end module run_tools
