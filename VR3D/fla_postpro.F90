program flapost
implicit none

    integer                   :: nstep,tstep,istep

    character(len=8)          :: fmt2 ! format descriptor
    character(len=8)          :: fmt5 ! format descriptor
    character(2)              :: batchchar
    character(5)              :: timechar
    character(5)              :: partchar
    character(5)              :: a
    real(kind=8)              :: r,st,c
    real(kind=8)              :: v(50000,3),u(50000,100),d(50000,10)


nstep=100
tstep=80
fmt5 = '(I5.5)' ! an integer of width 5 with zeros at the left
fmt2 = '(I2.2)' ! an integer of width 5 with zeros at the left

do istep=1,tstep
  do id=1,8
  write (timechar,fmt5) nstep
  write (batchchar,fmt2) id
  open (unit=22,file="part"//trim(timechar)//"size"//trim(batchchar)//".dat", form='formatted', position='rewind')
  read (22,"(A32,E20.8,A12,E20.8,A2)") a,r,b,st,c
  read (22,"(A32,100(A4,I2.2))") a
  do ip=1,np
    read (22,"(100(E20.8))") (v(ip,ii),ii=1,3),(u(ip,ii),ii=1,78),(prt%d(ii),ii=1,10)
  enddo
  close (22)


enddo


end program flapost
