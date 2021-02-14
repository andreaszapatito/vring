program flacond
implicit none

    integer, parameter        :: npmax=50000
    integer, parameter        :: nbin=100
    integer                   :: nstep,tstep,istep,id,ip,np,nb,nskip,ii,ic,ibin

    character(len=8)          :: fmt2 ! format descriptor
    character(len=8)          :: fmt5 ! format descriptor
    character(2)              :: batchchar
    character(5)              :: timechar
    character(5)              :: partchar
    real(kind=8)              :: c(npmax,10),cc(10,nbin),dc,x(npmax,3),jac(npmax),hmag(npmax,3),cbin
tstep=40
fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left
nskip=1
do istep=40,tstep
  do id=1,8
    nstep=istep*100
    write (timechar,fmt5) nstep
    write (batchchar,fmt2) id
    open (unit=200,file="post"//trim(timechar)//"size"//trim(batchchar)//".dat", form='formatted', position='rewind')
    do ip=1,npmax
      read (200,"(100(E20.8))",END=100) (x(ip,ii),ii=1,3),jac(ip),c(ip,1),hmag(ip,1),c(ip,2),c(ip,3),c(ip,4),c(ip,5)
      np=np+1
    enddo
    100 continue  
    close (200)
    dc=0.5
    nb=nbin
    do ip=1,np
      do ic=1,5
        ibin=min(floor(c(ip,ic)/dc)+1,nbin)
 !       write (*,*) "read",ic,ibin

        cc(ic,ibin)=cc(ic,ibin)+1.0/float(np)
      enddo
    enddo
    open (unit=300,file="cond"//trim(timechar)//"size"//trim(batchchar)//".dat", form='formatted', position='rewind')
    do ibin=1,nbin
      cbin=(float(ibin)-0.5)*dc
      write (300,"(100(E20.8))") ibin,cbin-0.5*dc,cbin,cbin+0.5*dc,cc(1,ibin),cc(2,ibin),cc(3,ibin),cc(4,ibin),cc(5,ibin)
    enddo
    close (300)

  enddo
enddo


end program flacond

