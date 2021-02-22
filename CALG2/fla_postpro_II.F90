program flapost

implicit none
interface 

   function jacnorm(j) result(hnorm)
     real(kind=8)             :: hnorm
     real(kind=8), intent(in) :: j(3,3)
   end function jacnorm

   function jachess(jxi,alfa,beta,gama) result(jhess)
     real(kind=8)             :: jhess(3,3)
     real(kind=8), intent(in) :: jxi(3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function jachess

   function jacgrad(jxi,alfa,beta,gama) result(jgrad)
     real(kind=8)             :: jgrad(3)
     real(kind=8), intent(in) :: jxi(3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function jacgrad

   function hessnorm(h) result(hnorm)
     real(kind=8)             :: hnorm
     real(kind=8), intent(in) :: h(3,3,3)
   end function hessnorm

   function hesshess(hxi,alfa,beta,gama) result(hhess)
     real(kind=8)             :: hhess(3,3)
     real(kind=8), intent(in) :: hxi(3,3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function hesshess

   function hessgrad(hxi,alfa,beta,gama) result(hgrad)
     real(kind=8)             :: hgrad(3)
     real(kind=8), intent(in) :: hxi(3,3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function hessgrad

   function hessinvr(h) result(hinv)
    real(kind=8)              :: hinv(3,3)
    real(kind=8),  intent(in) :: h(3,3)
   end function hessinvr

   function hesscg(h) result(hinv)
    real(kind=8)              :: hinv(3,3)
    real(kind=8),  intent(in) :: h(3,3)
   end function hesscg

 function BiCGStab(A,b) result(x)
        implicit none
        real    (kind=8), intent(in )                   :: A (:,:)
        real    (kind=8), intent(in )                   :: b ( : )
        real    (kind=8), dimension(1:size(b, dim=1))   :: x
 end
   function hesseta(hxi,alfa,beta,gama) result (h)
     real(kind=8)             :: h(3,3,3)
     real(kind=8), intent(in) :: hxi(3,3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function hesseta

   function jaceta(jxi,alfa,beta,gama) result (j)
     real(kind=8)             :: j(3,3)
     real(kind=8), intent(in) :: jxi(3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function jaceta


end interface 

    integer, parameter        :: npmax=50000
    integer                   :: nstep,tstep,istep,id,ii,jj,kk,ip,np,ir,nr,itest

    character(len=8)          :: fmt2 ! format descriptor
    character(len=8)          :: fmt5 ! format descriptor
    character(2)              :: batchchar
    character(5)              :: timechar
    character(5)              :: partchar
    character(500)             :: a,b,d
    real(kind=8)              :: r,st,c
    real(kind=8)              :: ra(3,3),rb(3,3),rc(3,3),vec(3,3),dir(npmax,6)
    real(kind=8)              :: v(npmax,3),x(npmax,3),u(npmax,3)
    real(kind=8)              :: j(npmax,3,3),h(npmax,3,3,3)
    real(kind=8)              :: o(npmax,3,3),p(npmax,3,3,3)
    real(kind=8)              :: detj(npmax,1),hmag(npmax,1),j1(npmax),j2(npmax)
    real(kind=8)              :: alfa,beta,gama
    real(kind=8)              :: hxi(3,3,3),heta(3,3,3)
    real(kind=8)              :: jxi(3,3),jeta(3,3)
    real(kind=8)              :: hgrad(3),jgrad(3)
    real(kind=8)              :: hhess(3,3),hhinv(3,3)
    real(kind=8)              :: jhess(3,3),jhinv(3,3)
    real(kind=8)              :: hnorm,hnorm0,hnorm1
    real(kind=8)              :: jnorm,jnorm0,jnorm1
    real(kind=8)              :: jnormmin
    real(kind=8)              :: hnormmin
    real(kind=8)              :: jnormmax
    real(kind=8)              :: hnormmax
    real(kind=8)              :: alpha1(3)
    integer                   :: itry,ikeep(3),jtry,ktry
    real(kind=8)              :: pi4,cmag 
    real(kind=8)              :: error,prod
    real(kind=8)              :: c1,c01,c001,c0001
    logical                   :: reading,ch
    real(kind=8)              :: ftol,fret,par(3,3,3),ppowel(3),xi(3,3)
    integer                   :: dof,iter,nskip
ch=.true.
tstep=20
nr=100
fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left
pi4=datan(1.d0)
nskip=1
do istep=10,tstep
  do id=1,1
    nstep=istep*50
    write (timechar,fmt5) nstep
    write (batchchar,fmt2) id
    open (unit=100,file="part"//trim(timechar)//"size"//trim(batchchar)//".dat", form='formatted', position='rewind')
    read (100,"(A200,E20.8,A200,E20.8,A200)") a,r,b,st,c
    write (*,*) "station 1"
    read (100,"(A500,100(A4,I2.2))") d
    write (*,*) "station 2"
    reading=.TRUE.
    np=0
    do ip=1,npmax
      read (100,"(200(E20.8))",END=100) (v(ip,ii),ii=1,3),(x(ip,ii),ii=1,3),(u(ip,ii),ii=1,3),((j(ip,ii,jj),jj=1,3),ii=1,3),((o(ip,ii,jj),jj=1,3),ii=1,3),(((h(ip,ii,jj,kk),ii=1,3),jj=1,3),kk=1,3),(((p(ip,ii,jj,kk),jj=1,3),ii=1,3),kk=1,3)
      np=np+1
    enddo
    100 continue  
    close (100)
!    h(:,:,:,:)=0.d0
!    h(1,1,1,1)=sqrt(3.d0)/2.d0
!    h(1,1,2,1)=sqrt(1.d0)/2.d0
!    h(1,2,1,1)=-sqrt(1.d0)/2.d0
!    h(1,2,2,1)=sqrt(3.d0)/2.d0
!
    do ip=1,np,nskip
      detj(ip,1)=(j(ip,1,1)*j(ip,2,2)*j(ip,3,3)+j(ip,1,2)*j(ip,2,3)*j(ip,3,1)+j(ip,1,3)*j(ip,2,1)*j(ip,3,2)   &
                 -j(ip,3,1)*j(ip,2,2)*j(ip,1,3)-j(ip,3,2)*j(ip,2,3)*j(ip,1,1)-j(ip,3,3)*j(ip,2,1)*j(ip,1,2)   )
!      do itest=1,100
!      h(1,3,3,3)=1.0
      do ii=1,3
        do jj=1,3
          jxi(ii,jj)=j(ip,ii,jj)
          do kk=1,3
            hxi(ii,jj,kk)=h(ip,ii,jj,kk)
          enddo
        enddo
      enddo
      alfa=0.d0;
      beta=0.d0;
      gama=0.d0;



  if (.true.) then
      hnormmin=10000.0
      hnormmax=-10000.0
      jnormmin=10000.0
      jnormmax=-10000.0
      do itry=0,2
        do jtry=0,2
          do ktry=0,2
      
      alfa=float(itry)*1.0d0
      beta=float(jtry)*1.0d0
      gama=float(ktry)*2.0d0


      do ir=1,nr
!        write (*,"(A10,100(E20.8))") 'hxi',hxi(1,:,:)
!        write (*,"(A10,100(E20.8))") 'hxi',hxi(2,:,:)
!        write (*,"(A10,100(E20.8))") 'hxi',hxi(3,:,:)
        heta=hesseta(hxi,alfa,beta,gama)
        jeta=jaceta(jxi,alfa,beta,gama)
!        write (*,"(A10,100(E20.8))") 'het',heta(1,:,:)
!        write (*,"(A10,100(E20.8))") 'het',heta(2,:,:)
!        write (*,"(A10,100(E20.8))") 'het',heta(3,:,:)
        !heta=hesseta(hxi,alfa,beta,gama)
        hnorm=hessnorm(heta)
        jnorm=jacnorm(jeta)
!        write (*,"(A10,100(E20.8))") 'hnr',hnorm
        hgrad=hessgrad(hxi,alfa,beta,gama)
        jgrad=jacgrad(jxi,alfa,beta,gama)
!        write (*,"(A10,100(E20.8))") 'hgr',hgrad
!        write (*,"(A10,100(E20.8))") 'alfa',alfa,beta,gama
        hhess=hesshess(hxi,alfa,beta,gama)
        jhess=jachess(jxi,alfa,beta,gama)
!       hhess(1,1)=1.d0
!       hhess(1,2)=0.50d0
!       hhess(1,3)=0.00
!
!       hhess(2,1)=0.0d0
!       hhess(2,2)=1.0d0
!       hhess(2,3)=0.0d0
!
!       hhess(3,1)=0.0d0
!       hhess(3,2)=0.0d0
!       hhess(3,3)=0.0d0
!
!        write (*,"(A10,100(E20.8))") 'hhs',hhess(1,:)
!        write (*,"(A10,100(E20.8))") 'hhs',hhess(2,:)
!        write (*,"(A10,100(E20.8))") 'hhs',hhess(3,:)
if(ch) then
       error=hgrad(1)**2+hgrad(2)**2+hgrad(3)**2
else
       error=jgrad(1)**2+jgrad(2)**2+jgrad(3)**2
endif

!       error=hgrad(1)**2
!        write (*,*) 'error 1',error,"alfa,beta,gama",alfa,beta,gama
!        hhinv=hesscg(hhess)
if(ch) then
        hhinv=hesscg(hhess)
        prod=0.d0
        do ii=1,3
          prod=prod+hhinv(1,ii)*hgrad(ii)
        enddo
        alfa=alfa-prod
        alfa=dmod(alfa,pi4*4.d0)
        prod=0.d0
        do ii=1,3
          prod=prod+hhinv(2,ii)*hgrad(ii)
        enddo
        beta=beta-prod
       beta=dmod(beta,pi4*4.d0)

        prod=0.d0
        do ii=1,3
          prod=prod+hhinv(3,ii)*hgrad(ii)
        enddo
        gama=gama-prod
        gama=dmod(gama,pi4*8.d0)
else
        jhinv=hesscg(jhess)
        prod=0.d0
        do ii=1,3
          prod=prod+jhinv(1,ii)*jgrad(ii)
        enddo
        alfa=alfa-prod
        alfa=dmod(alfa,pi4*4.d0)
        prod=0.d0
        do ii=1,3
          prod=prod+jhinv(2,ii)*jgrad(ii)
        enddo
        beta=beta-prod
       beta=dmod(beta,pi4*4.d0)

        prod=0.d0
        do ii=1,3
          prod=prod+jhinv(3,ii)*jgrad(ii)
        enddo
        gama=gama-prod
        gama=dmod(gama,pi4*8.d0)
endif
!        write (*,"(A10,100(E20.8))") 'hhi',hhinv(1,:)
!        write (*,"(A10,100(E20.8))") 'hhi',hhinv(2,:)
!        write (*,"(A10,100(E20.8))") 'hhi',hhinv(3,:)
       ! j1(ip)=abs(jeta(1,1))
       ! j1(ip)=abs(jeta(1,1))
       ! j1(ip)=abs(detj(ip,1))
       ! j2(ip)=abs(jeta(2,2)*jeta(3,3)-jeta(3,2)*jeta(2,3))
!        j2(ip)=abs(detj(ip,1)/j1(ip))
       ! j2(ip)=1.d0
!        write (*,*) "rotated",j1(ip),j2(ip),jeta
        if (error.lt.1.0e-7) exit

      enddo
      if (ch) then 
        if (hnorm>hnormmax) then
          ikeep(1)=itry
          ikeep(2)=jtry
          ikeep(3)=ktry
          hnormmax=hnorm
          alpha1(1)=alfa
          alpha1(2)=beta
          alpha1(3)=gama
        endif
      else

       if (jnorm>jnormmax) then
          ikeep(1)=itry
          ikeep(2)=jtry
          ikeep(3)=ktry
          jnormmax=jnorm
          alpha1(1)=alfa
          alpha1(2)=beta
          alpha1(3)=gama
        endif
     endif

!        write (*,*) 'part',ip,"nits",ir,"norm",hnorm,jnorm,"err",error,"alfa,beta,gama",alfa*45.d0/pi4,beta*45.d0/pi4,gama*45.d0/pi4
      enddo
      enddo
      enddo
      endif
      alfa=alpha1(1)
      beta=alpha1(2)
      gama=alpha1(3)
      heta=hesseta(hxi,alfa,beta,gama)
      jeta=jaceta(jxi,alfa,beta,gama)
      hnorm=hessnorm(heta)
      jnorm=jacnorm(jeta)

!      j1(ip)=abs(jeta(1,1))
      j1(ip)=abs(detj(ip,1))
     
!     j2(ip)=abs(jeta(2,2)*jeta(3,3)-jeta(3,2)*jeta(2,3))
      j2(ip)=abs(detj(ip,1)/j1(ip))
!      j2(ip)=1.0


      write (*,*) 'istep',istep,'id',id,'part',ip,"nkeep",ikeep,"norm",hnorm,"jnorms j1 j2 jdet",j1(ip),j2(ip),detj(ip,1),"err",error,"alfa,beta,gama",alfa*45.d0/pi4,beta*45.d0/pi4,gama*45.d0/pi4
      hmag(ip,1)=abs(heta(1,1,1))
      hmag(ip,1)=abs(heta(1,1,1)*heta(1,2,2)*heta(1,3,3)+heta(1,1,2)*heta(1,2,3)*heta(1,3,1)+heta(1,1,3)*heta(1,2,1)*heta(1,3,2)   &
                 -heta(1,3,1)*heta(1,2,2)*heta(1,1,3)-heta(1,3,2)*heta(1,2,3)*heta(1,1,1)-heta(1,3,3)*heta(1,2,1)*heta(1,1,2)   )
      call rabc(vec,alfa,beta,gama)

      dir(ip,1)=vec(1,1)
      dir(ip,2)=vec(2,1)
      dir(ip,3)=vec(3,1)
      if (dir(ip,3).gt.1.0) write (*,*) "what?",alfa,beta,gama,vec
    enddo
    open (unit=200,file="post"//trim(timechar)//"size"//trim(batchchar)//".dat", form='formatted', position='rewind')
    do ip=1,np,nskip


      r=0.01d0
      if (abs(j1(ip))**2>2.d0*hmag(ip,1)*r) then
        c1=(1.0/j2(ip))*2.0/(sqrt(abs(j1(ip))**2+2.d0*hmag(ip,1)*r)+sqrt(abs(j1(ip))**2-2.d0*hmag(ip,1)*r)) 
      else
        c1=(sqrt(abs(j1(ip))**2+2.d0*hmag(ip,1)*r)+sqrt(-abs(j1(ip))**2+2.d0*hmag(ip,1)*r))/(2.d0*hmag(ip,1)*r*j2(ip))
      endif

!      if (c1.gt.1.0/abs(j1(ip))) then
!              write (*,*) "Why?",c1,abs(j1(ip)),hmag(ip,1),r,abs(j1(ip))**2>2.d0*hmag(ip,1)*r
!      endif
      r=0.001d0
      if (abs(j1(ip))**2>2.d0*hmag(ip,1)*r) then
        c01=(1.0/j2(ip))*2.0/(sqrt(abs(j1(ip))**2+2.d0*hmag(ip,1)*r)+sqrt(abs(j1(ip))**2-2.d0*hmag(ip,1)*r)) 
      else
        c01=(sqrt(abs(j1(ip))**2+2.d0*hmag(ip,1)*r)+sqrt(-abs(j1(ip))**2+2.d0*hmag(ip,1)*r))/(2.d0*hmag(ip,1)*r*j2(ip))
      endif


      r=0.0001d0
      if (abs(j1(ip))**2>2.d0*hmag(ip,1)*r) then
        c001=(1.0/j2(ip))*2.0/(sqrt(abs(j1(ip))**2+2.d0*hmag(ip,1)*r)+sqrt(abs(j1(ip))**2-2.d0*hmag(ip,1)*r)) 
      else
        c001=(sqrt(abs(j1(ip))**2+2.d0*hmag(ip,1)*r)+sqrt(-abs(j1(ip))**2+2.d0*hmag(ip,1)*r))/(2.d0*hmag(ip,1)*r*j2(ip))
      endif


      r=0.00001d0
      if (abs(j1(ip))**2>2.d0*hmag(ip,1)*r) then
        c0001=(1.0/j2(ip))*2.0/(sqrt(abs(j1(ip))**2+2.d0*hmag(ip,1)*r)+sqrt(abs(j1(ip))**2-2.d0*hmag(ip,1)*r)) 
      else
        c0001=(sqrt(abs(j1(ip))**2+2.d0*hmag(ip,1)*r)+sqrt(-abs(j1(ip))**2+2.d0*hmag(ip,1)*r))/(2.d0*hmag(ip,1)*r*j2(ip))
      endif

      !cmag=1.d0/abs(j1(ip))
      write (200,"(100(E20.8))") (x(ip,ii),ii=1,3),j1(ip),1.d0/abs(detj(ip,1)),hmag(ip,1),c1,c01,c001,c0001,dir(ip,1),dir(ip,2),dir(ip,3),hmag(ip,1)*dir(ip,:)
    enddo
    close (200)

  enddo
enddo


end program flapost

function jacgrad(jxi,alfa,beta,gama) result (jgrad)
implicit none

interface
   function jacnorm(jac) result (jnorm)
     real(kind=8)             :: jnorm
     real(kind=8), intent(in) :: jac(3,3)
   end function jacnorm

   function jaceta(jxi,alfa,beta,gama) result (jac)
     real(kind=8)             :: jac(3,3)
     real(kind=8), intent(in) :: jxi(3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function jaceta
end interface
    real(kind=8)             :: jgrad(3)
    real(kind=8), intent(in) :: jxi(3,3)
    real(kind=8), intent(in) :: alfa,beta,gama

    real(kind=8)              :: da
    real(kind=8)  :: jnorm(3,3,3)
    real(kind=8)  :: jac(3,3)
    integer                   :: i,j,k
    integer                   :: l,m,n

    da=0.01d0
    do i=1,3
      do j=1,3
        do k=1,3
          jac(:,:)=jaceta(jxi,alfa+dfloat(i-2)*da,beta+dfloat(j-2)*da,gama+dfloat(k-2)*da)
!          jnorm(i,j,k)=(jac(1,1))**2
          jnorm(i,j,k)=1.d0/(jac(1,1))
        enddo
      enddo
    enddo

    jgrad(1)=(jnorm(3,2,2)-jnorm(1,2,2))/(2.d0*da)
    jgrad(2)=(jnorm(2,3,2)-jnorm(2,1,2))/(2.d0*da)
    jgrad(3)=(jnorm(2,2,3)-jnorm(2,2,1))/(2.d0*da)

end function jacgrad

function hessgrad(hxi,alfa,beta,gama) result (hgrad)
implicit none

interface
   function hessnorm(h) result (hnorm)
     real(kind=8)             :: hnorm
     real(kind=8), intent(in) :: h(3,3,3)
   end function hessnorm

   function hesseta(hxi,alfa,beta,gama) result (h)
     real(kind=8)             :: h(3,3,3)
     real(kind=8), intent(in) :: hxi(3,3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function hesseta
end interface
    real(kind=8)             :: hgrad(3)
    real(kind=8), intent(in) :: hxi(3,3,3)
    real(kind=8), intent(in) :: alfa,beta,gama

    real(kind=8)              :: da
    real(kind=8)  :: hnorm(3,3,3)
    real(kind=8)  :: h(3,3,3)
    integer                   :: i,j,k
    integer                   :: l,m,n

    da=0.01d0
!    write(*,*) 'angles',alfa,beta,gama
    do i=1,3
      do j=1,3
        do k=1,3
          h(:,:,:)=hesseta(hxi,alfa+dfloat(i-2)*da,beta+dfloat(j-2)*da,gama+dfloat(k-2)*da)
!          write (*,*) "alfas",alfa+dfloat(i-2)*da,beta+dfloat(j-2)*da,gama+dfloat(k-2)*da,h(1,1,1)
          hnorm(i,j,k)=h(1,1,1)
        enddo
      enddo
    enddo

    hgrad(1)=(hnorm(3,2,2)-hnorm(1,2,2))/(2.d0*da)
    hgrad(2)=(hnorm(2,3,2)-hnorm(2,1,2))/(2.d0*da)
    hgrad(3)=(hnorm(2,2,3)-hnorm(2,2,1))/(2.d0*da)
end function hessgrad

function jachess(jxi,alfa,beta,gama) result(jhess)
implicit none

interface
   function jacnorm(jac) result (jnorm)
     real(kind=8)             :: jnorm
     real(kind=8), intent(in) :: jac(3,3)
   end function jacnorm

   function jaceta(jxi,alfa,beta,gama) result (jac)
     real(kind=8)             :: jac(3,3)
     real(kind=8), intent(in) :: jxi(3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function jaceta


end interface

    real(kind=8)             :: jhess(3,3)
    real(kind=8), intent(in) :: jxi(3,3)
    real(kind=8), intent(in) :: alfa,beta,gama

    real(kind=8)  :: jnorm(3,3,3)
    real(kind=8)  :: jac(3,3)
    real(kind=8)              :: da
    integer                   :: i,j,k
    integer                   :: l,m,n
    da=0.01d0
    do i=1,3
      do j=1,3
        do k=1,3
          jac=jaceta(jxi,dfloat(i-2)*da+alfa,dfloat(j-2)*da+beta,dfloat(k-2)*da+gama)
!          jnorm(i,j,k)=(jac(1,1))**2
          jnorm(i,j,k)=1.0/jac(1,1)
        enddo
      enddo
    enddo

    jhess(1,1)=(jnorm(3,2,2)-2.0*jnorm(2,2,2)+jnorm(1,2,2))/(da**2.d0)
    jhess(1,2)=(jnorm(3,3,2)+jnorm(1,1,2)-jnorm(3,1,2)-jnorm(1,3,2))/(4.d0*da**2.d0)
    jhess(1,3)=(jnorm(3,2,3)+jnorm(1,2,1)-jnorm(3,2,1)-jnorm(1,2,3))/(4.d0*da**2.d0)

    jhess(2,1)=(jnorm(3,3,2)+jnorm(1,1,2)-jnorm(3,1,2)-jnorm(1,3,2))/(4.d0*da**2.d0)
    jhess(2,2)=(jnorm(2,3,2)-2.0*jnorm(2,2,2)+jnorm(2,1,2))/(da**2.d0)
    jhess(2,3)=(jnorm(2,3,3)+jnorm(2,1,1)-jnorm(2,3,1)-jnorm(2,1,3))/(4.d0*da**2.d0)

    jhess(3,1)=(jnorm(3,2,3)+jnorm(1,2,1)-jnorm(3,2,1)-jnorm(1,2,3))/(4.d0*da**2.d0)
    jhess(3,2)=(jnorm(2,3,3)+jnorm(2,1,1)-jnorm(2,3,1)-jnorm(2,1,3))/(4.d0*da**2.d0)
    jhess(3,3)=(jnorm(2,2,3)-2.0*jnorm(2,2,2)+jnorm(2,2,1))/(da**2.d0)

end function jachess


function hesshess(hxi,alfa,beta,gama) result(hhess)
implicit none

interface
   function hessnorm(h) result (hnorm)
     real(kind=8)             :: hnorm
     real(kind=8), intent(in) :: h(3,3,3)
   end function hessnorm

   function hesseta(hxi,alfa,beta,gama) result (h)
     real(kind=8)             :: h(3,3,3)
     real(kind=8), intent(in) :: hxi(3,3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function hesseta


end interface

    real(kind=8)             :: hhess(3,3)
    real(kind=8), intent(in) :: hxi(3,3,3)
    real(kind=8), intent(in) :: alfa,beta,gama

    real(kind=8)  :: hnorm(3,3,3)
    real(kind=8)  :: h(3,3,3)
    real(kind=8)              :: da
    integer                   :: i,j,k
    integer                   :: l,m,n
    da=0.01d0
    do i=1,3
      do j=1,3
        do k=1,3
          h=hesseta(hxi,dfloat(i-2)*da+alfa,dfloat(j-2)*da+beta,dfloat(k-2)*da+gama)
          hnorm(i,j,k)=h(1,1,1)
        enddo
      enddo
    enddo
    
    hhess(1,1)=(hnorm(3,2,2)-2.0*hnorm(2,2,2)+hnorm(1,2,2))/(da**2.d0)
    hhess(1,2)=(hnorm(3,3,2)+hnorm(1,1,2)-hnorm(3,1,2)-hnorm(1,3,2))/(4.d0*da**2.d0)
    hhess(1,3)=(hnorm(3,2,3)+hnorm(1,2,1)-hnorm(3,2,1)-hnorm(1,2,3))/(4.d0*da**2.d0)

    hhess(2,1)=(hnorm(3,3,2)+hnorm(1,1,2)-hnorm(3,1,2)-hnorm(1,3,2))/(4.d0*da**2.d0)
    hhess(2,2)=(hnorm(2,3,2)-2.0*hnorm(2,2,2)+hnorm(2,1,2))/(da**2.d0)
    hhess(2,3)=(hnorm(2,3,3)+hnorm(2,1,1)-hnorm(2,3,1)-hnorm(2,1,3))/(4.d0*da**2.d0)

    hhess(3,1)=(hnorm(3,2,3)+hnorm(1,2,1)-hnorm(3,2,1)-hnorm(1,2,3))/(4.d0*da**2.d0)
    hhess(3,2)=(hnorm(2,3,3)+hnorm(2,1,1)-hnorm(2,3,1)-hnorm(2,1,3))/(4.d0*da**2.d0)
    hhess(3,3)=(hnorm(2,2,3)-2.0*hnorm(2,2,2)+hnorm(2,2,1))/(da**2.d0)

end function hesshess

subroutine rabc(vec,alfa,beta,gama) 
implicit none
    real(kind=8), intent(in)  :: alfa,beta,gama
    real(kind=8)              :: ra(3,3),rb(3,3),rc(3,3),re(3,3)
    real(kind=8)              :: rab(3,3),rba(3,3),vec(3,3)
    real(kind=8)              :: nx(3,3),xn(3,3)
    integer                   :: i,j,k
    integer                   :: l,m,n

    real(kind=8)              :: c1,c2,c3
    real(kind=8)              :: s1,s2,s3

    c1=cos(alfa)
    c2=cos(beta)
    c3=cos(gama)

    s1=sin(alfa)
    s2=sin(beta)
    s3=sin(gama)

    

    re(1,1)=c2*c3-c1*s2*s3;  re(1,2)=s2*c3+c1*c2*s3; re(1,3)=s1*s3;
    re(2,1)=-c2*s3-c1*s2*c3; re(2,2)=-s2*s3+c1*c2*c3;re(2,3)=s1*c3;
    re(3,1)=s1*s2;           re(3,2)=-s1*c2;         re(3,3)=c1;




    do i=1,3
      do j=1,3
        vec(i,j)=re(j,i)
      enddo
    enddo

    
end subroutine rabc

function jaceta(jxi,alfa,beta,gama) result (jac)
implicit none
    real(kind=8)              :: jac(3,3)
    real(kind=8), intent(in)  :: jxi(3,3)
    real(kind=8), intent(in)  :: alfa,beta,gama
    real(kind=8)              :: ra(3,3),rb(3,3),rc(3,3),re(3,3)
    real(kind=8)              :: rab(3,3),rabc(3,3),rba(3,3),rcba(3,3)
    real(kind=8)              :: nx(3,3),xn(3,3)
    integer                   :: i,j,k
    integer                   :: l,m,n
    real(kind=8)              :: c1,c2,c3
    real(kind=8)              :: s1,s2,s3

    c1=cos(alfa)
    c2=cos(beta)
    c3=cos(gama)

    s1=sin(alfa)
    s2=sin(beta)
    s3=sin(gama)

    re(1,1)=c2*c3-c1*s2*s3;  re(1,2)=s2*c3+c1*c2*s3; re(1,3)=s1*s3;
    re(2,1)=-c2*s3-c1*s2*c3; re(2,2)=-s2*s3+c1*c2*c3;re(2,3)=s1*c3;
    re(3,1)=s1*s2;           re(3,2)=-s1*c2;         re(3,3)=c1;


    do i=1,3
      do j=1,3
        nx(i,j)=re(i,j)
        xn(i,j)=re(j,i)
      enddo
    enddo



    do i=1,3
      do j=1,3
          jac(i,j)=0.d0
          do l=1,3
            do m=1,3
                jac(i,j)=jac(i,j)+nx(i,l)*jxi(l,m)*xn(m,j)
            enddo
          enddo
      enddo
    enddo
 !   write (*,*) 'xn',xn,'nx',nx
end function jaceta


function hesseta(hxi,alfa,beta,gama) result (h)
implicit none
    real(kind=8)              :: h(3,3,3)
    real(kind=8), intent(in)  :: hxi(3,3,3)
    real(kind=8), intent(in)  :: alfa,beta,gama
    real(kind=8)              :: ra(3,3),rb(3,3),rc(3,3),re(3,3)
    real(kind=8)              :: rab(3,3),rabc(3,3),rba(3,3),rcba(3,3)
    real(kind=8)              :: nx(3,3),xn(3,3)
    integer                   :: i,j,k
    integer                   :: l,m,n
    real(kind=8)              :: c1,c2,c3
    real(kind=8)              :: s1,s2,s3

    c1=cos(alfa)
    c2=cos(beta)
    c3=cos(gama)

    s1=sin(alfa)
    s2=sin(beta)
    s3=sin(gama)

    re(1,1)=c2*c3-c1*s2*s3;  re(1,2)=s2*c3+c1*c2*s3; re(1,3)=s1*s3;
    re(2,1)=-c2*s3-c1*s2*c3; re(2,2)=-s2*s3+c1*c2*c3;re(2,3)=s1*c3;
    re(3,1)=s1*s2;           re(3,2)=-s1*c2;         re(3,3)=c1;





    do i=1,3
      do j=1,3
        nx(i,j)=re(i,j)
        xn(i,j)=re(j,i)
      enddo
    enddo
    do i=1,3
      do j=1,3
        do k=1,3
          h(i,j,k)=0.d0
          do l=1,3
            do m=1,3
              do n=1,3
                h(i,j,k)=h(i,j,k)+nx(i,l)*hxi(l,m,n)*xn(n,k)*xn(m,j)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
end function hesseta

function jacnorm(j) result(jnorm)
implicit none
    real(kind=8)             :: jnorm
    real(kind=8), intent(in) :: j(3,3)

!    jnorm=(j(1,1))**2
    jnorm=1.d0/(j(1,1))

end function jacnorm


function hessnorm(h) result(hnorm)
implicit none
    real(kind=8)             :: hnorm
    real(kind=8), intent(in) :: h(3,3,3)

    hnorm=h(1,1,1)    
end function hessnorm


function hessinvr(h) result(hinv)
implicit none
    integer :: stype
    real(kind=8)              :: hinv(3,3)
    real(kind=8),  intent(in) :: h(3,3)
    real(kind=8)              :: c(3,3)

    real(kind=8)              :: d,d1,d2,d3



      d= h(1,1)*h(2,2)*h(3,3)  &
        -h(1,1)*h(2,3)*h(3,2)  &
        -h(1,2)*h(2,1)*h(3,3)  &
        +h(1,2)*h(2,3)*h(3,1)  &
        +h(1,3)*h(2,1)*h(3,2)  &
        -h(1,3)*h(2,2)*h(3,1)

      if (d>1.0e-8) then

        c(1,1)=+(h(2,2)*h(3,3)-h(2,3)*h(3,2))
        c(1,2)=-(h(2,1)*h(3,3)-h(2,3)*h(3,1))
        c(1,3)=+(h(2,1)*h(3,2)-h(2,2)*h(3,1))
        c(2,1)=-(h(1,2)*h(3,3)-h(1,3)*h(3,2))
        c(2,2)=+(h(1,1)*h(3,3)-h(1,3)*h(3,1))
        c(2,3)=-(h(1,1)*h(3,2)-h(1,2)*h(3,1))
        c(3,1)=+(h(1,2)*h(2,3)-h(1,3)*h(2,2))
        c(3,2)=-(h(1,1)*h(2,3)-h(1,3)*h(2,1))
        c(3,3)=+(h(1,1)*h(2,2)-h(1,2)*h(2,1))
    
        hinv=transpose(c)/d
        stype=1
      else
        d1= h(1,1)*h(2,2)  &
           -h(1,2)*h(2,1)  
        d2= h(1,1)*h(3,3)  &
           -h(1,3)*h(3,1)  
        d3= h(2,2)*h(3,3)  &
           -h(3,2)*h(2,3)  

        
        if (dabs(d1)>dabs(d2).and.dabs(d1)>dabs(d3).and.dabs(d1)>1.e-8) then
                stype=2
        elseif  (dabs(d2)>dabs(d1).and.dabs(d2)>dabs(d3).and.dabs(d2)>1.e-8) then
                stype=3
        elseif  (dabs(d3)>dabs(d1).and.dabs(d3)>dabs(d2).and.dabs(d3)>1.e-8) then
                stype=4
        elseif  (1.e-8>dabs(d1).and.1.e-8>dabs(d2).and.1.e-8>dabs(d3)) then
                stype=5
        endif


        if (stype.eq.2) then
          c(1,1)=+(h(2,2))
          c(1,2)=-(h(2,1))
          c(1,3)=0.d0
          c(2,1)=-(h(1,2))
          c(2,2)=+(h(1,1))
          c(2,3)=0.d0
          c(3,1)=0.d0
          c(3,2)=0.d0
          c(3,3)=0.d0
  
          hinv=transpose(c)/d1
       elseif (stype.eq.3) then
         c(1,1)=+(h(3,3))
         c(1,2)=0.d0
         c(1,3)=+(-h(3,1))
         c(2,1)=0.d0
         c(2,2)=0.d0
         c(2,3)=0.d0
         c(3,1)=+(-h(1,3))
         c(3,2)=0.d0
         c(3,3)=+(h(1,1))

         hinv=transpose(c)/d2
       elseif (stype.eq.4) then
         c(1,1)=0.d0
         c(1,2)=0.d0
         c(1,3)=0.d0
         c(2,1)=0.d0
         c(2,2)=+(h(3,3))
         c(2,3)=-(h(3,2))
         c(3,1)=0.d0
         c(3,2)=-(h(2,3))
         c(3,3)=+(h(2,2))

         hinv=transpose(c)/d3
       elseif (stype.eq.5) then
         hinv(:,:)=0.d0
       endif
     endif
!     write (*,*) "System inverse type",stype
end function hessinvr



function hesscg(h) result(hinv)
implicit none
!interface
!function BiCGStab(A,b) result(x)
!        implicit none
!        real    (kind=8), intent(in )                   :: A (:,:)
!        real    (kind=8), intent(in )                   :: b ( : )
!        real    (kind=8), dimension(1:size(b, dim=1))   :: x
! end
! end interface

    integer                   :: i,j,k
    real(kind=8),  intent(in) :: h(3,3)
    real(kind=8)              :: hinv(3,3)
    real(kind=8)              :: s1uc(3)
    real(kind=8)              :: a(3,3),u(3,3),v(3,3),w(3),c(3),wu(3,3)
 !   real(kind=8)              :: c(3),ax(3),r(3),p(3),Ap(3),x(3),x1(3)
 !   real(kind=8)              :: r1(3),rr1
 !   real(kind=8)              :: rr,pAp,a,b
 !   integer                   :: ic,ik,i,j

!    do ic=1,3
!      c(1)=0.d0
!      c(2)=0.d0
!      c(3)=0.d0
!      x(1)=0.d0
!      x(2)=0.d0
!      x(3)=0.d0
!      c(ic)=1.d0
!      do i=1,3
!        ax(i)=0.d0
!        do j=1,3
!          ax(i)=ax(i)+h(i,j)*x(j)
!        enddo
!        r(i)=c(i)-ax(i)
!        r1(i)=r(i)
!        p(i)=r(i)
!      enddo
!      rr1=0.d0
!      do i=1,3
!        rr1=rr1+r1(i)**2
!      enddo
!      ik=0
!      do while (rr1.gt.1.0e-8)
!        ik=ik+1
!        rr=0.d0
!        do i=1,3
!          rr=rr+r(i)*r(i)
!        enddo
!
!        do i=1,3
!          Ap(i)=0.d0
!          do j=1,3
!            Ap(i)=Ap(i)+h(i,j)*p(j)
!          enddo
!        enddo
!
!
!        pAp=0.d0
!        do i=1,3
!          pAp=pAp+p(i)*Ap(i)
!        enddo
!        write (*,*),ic,ik,"rr",rr,"pAp",pAp,"p",p,"Ap",Ap,"mat",h
!        a=rr/pAp
!
!        do i=1,3
!          x1(i)=x(i)+a*p(i)
!          r1(i)=r(i)-a*Ap(i)
!        enddo
!        rr1=0.d0
!        do i=1,3
!          rr1=rr1+r1(i)*r1(i)
!        enddo
!        if (rr1.lt.1.0e-2) exit
!        if (ik.gt.100) exit
!
!        b=rr1/rr
!
!        do i=1,3
!          p(i)=r1(i)+b*p(i)
!        enddo
!        do i=1,3
!          r(i)=r1(i)
!          x(i)=x1(i)
!        enddo
!      enddo
!      
!      call r8ge_cg(3,a,c,x)
!      x=BiCGStab(h,c)
!      do i=1,3
!        hinv(i,ic)=x(i)
!      enddo
!    write (*,*) c,x
!    enddo
!c    Given a matrix (1:m,1:n) with physical dimensions mp by np,
!c    this routine computes its singular value decomposition,
!c    A = U W VT.  The matrix U replaces a on output.  The diagonal
!c    matrix of singular values W is output as a vector w(1:n)
!c    The matrix V (not the transpose VT) is the output as v(1:n,1:n) 
!ca

      do i=1,3
        do j=1,3
          u(i,j)=h(i,j)
        enddo
      enddo

!      write (*,*) "a",u
      call svdcmp(u,3,3,3,3,w,v)
!      write (*,*) "u",u
!      write (*,*) "v",v
!      write (*,*) "w",w


      hinv=0.d0
      do i=1,3
        if (w(i).gt.1.0e-20) then
          do j=1,3
            wu(i,j)=u(j,i)/w(i)
          enddo
        else  
          do j=1,3
            wu(i,j)=0.d0
          enddo
        endif
      enddo

      do i=1,3
        do j=1,3
          do k=1,3
            hinv(i,j)=hinv(i,j)+v(i,k)*wu(k,j)
          enddo
        enddo
      enddo

!      write (*,*) "w",wu,hinv
!      do ic=1,3
!        do i=1,3
!          do j=1,3
!           u(i,j)=h(i,j)
!          enddo
!        enddo
!
!        c=0.d0
!        c(ic)=1.d0
!
!        write (*,*) "a",u
!        write (*,*) "c",c
!        call svdcmp(u,3,3,3,3,w,v)
!        write (*,*) "u",u
!        write (*,*) "v",v
!        write (*,*) "w",w
!
!        do i=1,3
!          s1uc(i)=0.d0
!          if (abs(w(i))>1.0e-12) then
!            do j=1,3
!              s1uc(i)=s1uc(i)+u(j,i)*c(j)
!            enddo
!            if (abs(w(i))>1.0e-12) then
!
!            s1uc(i)=s1uc(i)/w(i)
!    endif
!          endif
!        enddo
!        do i=1,3
!          hinv(i,ic)=0.d0
!          do j=1,3
!            hinv(i,ic)=hinv(i,ic)+v(i,j)*s1uc(j)
!          enddo
!        enddo
!      enddo

end function hesscg


 function BiCGStab(A,b) result(x)
        implicit none

!--------------------------PARAMETER AND VARIABLE-------------------------------!
        real    (kind=8), intent(in )                   :: A (:,:)
        real    (kind=8), intent(in )                   :: b ( : )
        real    (kind=8), dimension(1:size(b, dim=1))   :: x

        real    (kind=8), dimension(1:size(b, dim=1))   :: r, rs, v, p, s, t
        real    (kind=8), parameter                     :: e = 1d-8
        real    (kind=8)                                :: rho      , rho_prev
        real    (kind=8)                                :: alpha    , omega   , beta
        real    (kind=8)                                :: norm_r   , norm_b       
        real    (kind=8)                                :: summesion, temp

        integer                                         :: it=0,err
!------------------------END PARAMETER AND VARIABLE-----------------------------!  

        if(size(A, dim=1) /= size(A, dim=2)) stop &
        "Error: Improper dimension of matrix A in BiCGStab."

!-------------------------------------------------------!
        x  = 0.0d0                                     !-------> INITIAL GUESS
!-------------------------------------------------------!
        r  = b - matmul(A,x)                            !-------> LINE 1
        rs = r                                          !
!-------------------------------------------------------!
        rho   = 1.0d0; alpha = 1.0d0; omega = 1.0d0  !-------> LINE 2
!-------------------------------------------------------!
        v  = 0.0d0; p  = 0.0d0                        !-------> LINE 3
!                                                       !
        norm_r = sqrt(dot_product(r,r))                 !
        norm_b = sqrt(dot_product(b,b))                 !
!-------------------------------------------------------!


        write (*,*) "in bcgstab x",x
        write (*,*) "in bcgstab b",b
        write (*,*) "in bcgstab a",A
        write (*,*) "in bcgstab r",r
        write (*,*) "in bcgstab v",v

       
        do while(norm_r .GT. e*norm_b.and.norm_b.gt.e)          !-------> START OF LOOP

        !-------------------------------------------------------!
            rho_prev = rho                                      !-------> LINE 5
            rho      = dot_product(rs,r)                        !
        !-------------------------------------------------------!
            beta     = (rho/rho_prev) * (alpha/omega)           !-------> LINE 6
        !-------------------------------------------------------!
            p        = r + beta * (p - omega*v)                 !-------> LINE 7
        !-------------------------------------------------------!
            v        = matmul(A,p)                              !-------> LINE 8
        !-------------------------------------------------------!
            alpha    = rho/dot_product(rs,v)                    !-------> LINE 9
        !-------------------------------------------------------!
            s        = r - alpha*v                              !-------> LINE 10
        !-------------------------------------------------------!
            t        = matmul(A,s)                              !-------> LINE 11
        !-------------------------------------------------------!
            omega    = dot_product(t,s)/dot_product(t,t)        !-------> LINE 12
        !-------------------------------------------------------!
            x        = x + alpha*p + omega*s                    !-------> LINE 13
        !-------------------------------------------------------!
            r        = s - omega*t                              !-------> LINE 17
        !-------------------------------------------------------!
            norm_r   = sqrt(dot_product(r,r))                   !
            norm_b   = sqrt(dot_product(b,b))                   !
        !-------------------------------------------------------!
            it = it + 1                                         !
        !-------------------------------------------------------!
        write (*,*) "in bcgstab x it",it,x,norm_r,norm_b

        end do                                                      !-------> END OF LOOP

        print*, "Iteration Required :", it

return
end function BiCGStab     


REAL*8 FUNCTION FUNC(P,PAR)
        implicit none
        interface
   function hesseta(hxi,alfa,beta,gama) result (h)
     real(kind=8)             :: h(3,3,3)
     real(kind=8), intent(in) :: hxi(3,3,3)
     real(kind=8), intent(in) :: alfa,beta,gama
   end function hesseta
   end interface
REAL*8 P(3),PAR(3,3,3),h(3,3,3)



h=hesseta(PAR,P(1),P(2),P(3))
FUNC=-h(1,1,1)
RETURN
END



