program flapost
implicit none

integer nstep,tstep,istep

nstep=100
tstep=80
fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left

do istep=1,tstep

  write (timechar,fmt5) par%nstep
  write (batchchar,fmt2) prt%id
  open (unit=22,file="part"//trim(timechar)//"size"//trim(batchchar)//".dat", form='formatted', position='rewind')
  read (22,"(A32,E20.8,A12,E20.8,A2)") a,r,b,st,c
  read (22,"(A32,100(A4,I2.2))") a
  do ip=1,np
    read (prt%id+22,"(100(E20.8))") (prt%v(ip,ii),ii=1,3),(prt%u(ip,ii),ii=1,prt%nvar),(prt%d(ii),ii=1,prt%nder)
  enddo
  close (22)


enddo


end program flapost
