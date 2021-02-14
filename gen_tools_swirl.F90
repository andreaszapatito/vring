module gen_tools
!#include <petsc/finclude/petscdef.h>
!  use petsc

  implicit none

  interface
   real(kind=8) function bessj0(x)
    real(kind=8), intent(in) :: x
   end function bessj0
  end interface

  interface
   real(kind=8) function bessj1(x)
    real(kind=8), intent(in) :: x
   end function bessj1
  end interface

  interface
   real(kind=8) function bessi0(x)
    real(kind=8), intent(in) :: x
   end function bessi0
  end interface

  interface
   real(kind=8) function bessi1(x)
    real(kind=8), intent(in) :: x
   end function bessi1
  end interface

  interface
   real(kind=8) function bessk1(x)
    real(kind=8), intent(in) :: x
   end function bessk1
  end interface



          contains

          subroutine read_until(unit, string, flag)
            implicit none
            integer          :: unit
            character(len=*) :: string
            logical          :: flag

            character(len=128) :: ligne
            integer            :: iostatus

            rewind(unit)
            flag=.true.
            do while(flag)
               read(unit,'(128a)',iostat=iostatus) ligne

               if (trim(adjustl(ligne))==trim(adjustl(string))) return
               if (iostatus.gt.0) then
                  print*, 'read error'
                  flag=.false.
               end if
               if (iostatus.lt.0) then
                  print*, 'end of file'
                  flag=.false.
               end if
            end do

          end subroutine read_until


          subroutine err_n_stop(string, ranks)
            implicit none
            character(len=*)      :: string
            integer, dimension(:) :: ranks
            integer          :: ierr

            !PetscErrorCode   :: ierr

            if (sum(ranks)==0) then
               print*, 'error: '//trim(adjustl(string))
            end if
            call mpi_finalize(ierr)
            stop

          end subroutine err_n_stop

          subroutine warning(string, ranks)
            implicit none
            character(len=*)      :: string
            integer, dimension(:) :: ranks
            integer          :: ierr

            !PetscErrorCode   :: ierr

            if (sum(ranks)==0) then
               print*, 'warning: '//trim(adjustl(string))
            end if

          end subroutine warning

          subroutine par_print(string, ranks)
            implicit none
            character(len=*)      :: string
            integer, dimension(:) :: ranks

!            PetscErrorCode   :: ierr
            integer          :: ierr

            if (sum(ranks)==0) then
               print*, string
            end if

          end subroutine par_print

          function num2str(n) result(str)
            implicit none
            integer          :: n
            character(len=4) :: str

            write(str,'(i4)') n

          end function  num2str

          function rke2str(xx) result(str)
            implicit none
            real(kind=8)      :: xx
            character(len=64) :: str

            write(str,*) xx

          end function  rke2str

          function profilet(t,r,th) result(u)
            implicit none
            real(kind=8)      :: t,tt,r,th,pi
            real(kind=8)      :: u
            tt=t
            u=0.0
            pi=4.0*atan(1.0)
            if (r<0.5) then
              u=min(0.2,t*r)
              u=-0.13*4.0*r*r*sin(10.0*2.0*pi*t)*sin(10.0*2.0*pi*r)*cos(16.0*th)
            endif
           u=0.0
          end function  profilet
          function profiler(t,r,th) result(u)
            implicit none
            real(kind=8)      :: t,tt,r,th,pi
            real(kind=8)      :: u
            tt=t
            u=0.0
            pi=4.0*atan(1.0)
            if (r<0.5) then
              u=min(0.2,t*r)
              u=0.13*4.0*r*r*sin(5.0*2.0*pi*t)*cos(5.0*2.0*pi*r)*sin(16.0*th)
            endif

            u=0.0
          end function  profiler

          function profile(t,r,th) result(u)
            implicit none
            real(kind=8)      :: t,tt,r,th,pi,n,t1,t2
            real(kind=8)      :: u
            tt=t
!           tt=mod(t,2.0)
            pi=4.0*atan(1.0)
            n=2.0
            t1=0.1*n
            t2=0.9*n
            
            if (tt.gt.0.0.and.tt.lt.t1) u=3.0*(tt/t1)**2-2.0*(tt/t1)**3
            if (tt.gt.t1.and.tt.lt.n) u=1.0
            if (tt.gt.n.and.tt.lt.n+t1) u=3.0*((n+t1-tt)/t1)**2-2.0*((n+t1-tt)/t1)**3
            if (tt.gt.n+t1) u=0.0
            u=u*(1.0+0.d0*0.13*sin(5.0*2.0*pi*t))
          end function  profile

          function uvring(r1,x1,t1,epsilon) result(u)
            implicit none
            real(kind=8)      :: r,x,du,mu,th,EF,K1,I1,I0,dm,pi,ss,epsilon,x1,r1,t1
            real(kind=8)      :: u
            integer           :: nm=1000,im

            th=1.0/sqrt(2.0*t1)
!            th=1.0
            th=t1
            dm=10.0/float(nm)
            u=0.0
            pi=4.0*atan(1.0)


            do im=1,nm
              mu=(im-0.5)*dm
              ss=1.0
              EF=exp((abs(x1))*th*mu)*ERFC((mu+(abs(x1))*th)/sqrt(2.0))+exp(-(abs(x1))*th*mu)*ERFC((mu-(abs(x1))*th)/sqrt(2.0))
              du=mu*(-0.25*(th**2)*EF*BESSJ0((r1)*th*mu)*BESSJ1(th*mu)+(1.0/(2.0*pi))*(BESSK1(mu/epsilon)/(BESSI1(mu/epsilon)))*BESSI0((r1)*mu)*BESSI1(mu)*cos(mu*(x1)))
              u=u+ss*du*dm
            enddo

            if (isnan(u)) u=0.0
          end function  uvring

          function vvring(r1,x1,t1,epsilon) result(v)
            implicit none
            real(kind=8)      :: r,x,dv,mu,th,EF,K1,I1,I0,dm,pi,ss,epsilon,x1,r1,t1
            real(kind=8)      :: v
            integer           :: nm=1000,im

            th=1.0/sqrt(2.0*t1)
!            th=1.0
            th=t1
            dm=10.0/float(nm)
            v=0.0
            pi=4.0*atan(1.0)

            do im=1,nm
              mu=(im-0.5)*dm
              ss=1.0
              EF=exp(((x1))*th*mu)*ERFC((mu+((x1))*th)/sqrt(2.0))-exp(-((x1))*th*mu)*ERFC((mu-((x1))*th)/sqrt(2.0))
              dv=mu*(0.25*(th**2)*EF*BESSJ1((r1)*th*mu)*BESSJ1(th*mu)+(1.0/(2.0*pi))*(BESSK1(mu/epsilon)/(BESSI1(mu/epsilon)))*BESSI1((r1)*mu)*BESSI1(mu)*sin(mu*(x1)))
              v=v+ss*dv*dm
            enddo
            if (isnan(v)) v=0.0
          end function  vvring

         function wvring(r1,x1,t1,epsilon) result(w)
            implicit none
            real(kind=8)      :: r,x,dv,mu,th,EF,K1,I1,I0,dm,pi,ss,epsilon,x1,r1,t1
            real(kind=8)      :: w
            integer           :: nm=1000,im

            th=1.0/sqrt(2.0*t1)
!           th=1.0
            th=t1
            dm=10.0/float(nm)
            w=0.0
            pi=4.0*atan(1.0)

            w=(th**3)/sqrt(2*pi)*exp(-(((r1)**2+(x1)**2+1.0)*th**2)/2.0)* BESSI1((r1)*th**2)
            if (isnan(w)) w=0.0

          end function  wvring


end module gen_tools


