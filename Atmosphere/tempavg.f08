program tempavg
!*******************************************************************************
!* Created by Stephen J. Houston 4.11.18
!*******************************************************************************
!* Averages two temperature profiles gathered by Juno's infrared camera.
!*******************************************************************************
implicit real*8(a-h,o-z)

parameter(nlinesmax=1000)

real*8,allocatable,dimension(:,:) :: alt,press,temp

open(unit=100,file="./Tp_april2016_70N120W.dat",status='old')
open(unit=101,file="./Tp_april2016_70N180W.dat",status='old')
open(unit=102,file="./IRTF-TEXES_april2016_temp.dat",status='unknown')

read(100,*) !Skip the header
read(101,*)

do i=1,nlinesmax
  read(100,*,end=1000)
  nlines=nlines+1
end do
1000 continue
rewind(100)
allocate(alt(2,nlines))
allocate(press(2,nlines))
allocate(temp(2,nlines))
read(100,*) !Skip the header again
do i=1,nlines
  read(100,*) alt(1,i),press(1,i),temp(1,i)
  read(101,*) alt(2,i),press(2,i),temp(2,i)
end do
close(100)
close(101)
write(102,*) '# Alt. above 1bar Pressure[mbar] Temp[K]'
do i=1,nlines
  write(102,10000) alt(1,i),press(1,i),(temp(1,i)+temp(2,i))/2
end do
close(102)
10000 format(4x,F10.3,6x,ES11.5E2,F10.3)
end program
