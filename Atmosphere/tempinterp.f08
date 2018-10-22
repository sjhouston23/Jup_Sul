program tempinterp
!*******************************************************************************
!* Created by Stephen J. Houston 4.13.18
!*******************************************************************************
!* Reads in temperature, pressure, and density profiles and creates an
!* interpolated atmosphere with a delta_z of 2 km. Also calculates the
!* corresponding scale height.
!*******************************************************************************

implicit real*8(a-h,o-z)

real*8 Kb
parameter(nmaxlines=100000,Kb=1.38064852E-19) !Boltzmann constant in cm^2

real*8 CD_H2, CD_He, CD_CH4, CD_H,TotCD
real*8,allocatable,dimension(:) :: AltJ,PressJ,TempJ,AltS,PressS,TempS,TempSI
real*8,allocatable,dimension(:) :: AltSP,Temp2,AltI,PressI,TempI,PressI2,TempI2
real*8,dimension(1544) :: AltIS,PressIS,TempIS,nH2I,nHeI,nCH4I,nHI,nTotI,H
real*8,dimension(1400) :: nH2,nHe,nCH4,nH,nTot,mrH2,mrHe,mrCH4,mrH,mrTot

!Input Juno's reading of the lower atmosphere temperatures
open(unit=100,file='./IRTF-TEXES_april2016_temp.dat',status='old')
open(unit=101,file='./Temp.dat',status='old') !Input upper atmosphere temp
open(unit=102,file='./TempI.dat',status='unknown') !Output interp. temp
open(unit=103,file='./AltDensPress.dat',status='old') !Input
!Input current atmosphere density profile
open(unit=104,file='./Input/atmosphere2km.dat',status='old')
!The rest are outputs
open(unit=105,file='./Input/JunoAtmosphere_2km.dat',status='unknown')
open(unit=106,file='./Input/JunoColumnDensity_2km.dat',status='unknown')

read(100,*) !Skip header
read(103,*)
nlinesJ=0;nlinesS=0;nlinesSP=0
do i=1,nmaxlines !Count the number of lines in this file
  read(100,*,end=1000)
  nlinesJ=nlinesJ+1
end do
write(*,*) 'nmaxlines exceeded 1'
1000 continue
do i=1,nmaxlines !Count the number of lines in this file
  read(101,*,end=1001) dum
  nlinesS=nlinesS+1
end do
write(*,*) 'nmaxlines exceeded 2'
1001 continue
do i=1,nmaxlines !Count the number of lines in this file
  read(103,*,end=1003) dum
  nlinesSP=nlinesSP+1
end do
write(*,*) 'nmaxlines exceeded 3'
1003 continue
rewind(100) !Go back to the beginning of the files
rewind(101)
rewind(103)
read(100,*) !Skip header again
read(103,*)
!Allocate the arrays based on file length
allocate(AltJ(nlinesJ))  !Juno atmosphere
allocate(PressJ(nlinesJ)) !Juno pressure
allocate(TempJ(nlinesJ)) !Juno temperature
allocate(AltS(nlinesS)) !Altitude that goes with my temp
allocate(PressS(nlinesSP)) !My previous pressures (these will stay the same)
allocate(TempS(nlinesS)) !Discrete temperature that I had (from Galileo)
allocate(TempSI(nlinesSP)) !Temp to be put into same alt bins as "PressS"
allocate(AltSP(nlinesSP)) !Altitude bins from "PressS"
allocate(Temp2(nlinesS)) !2nd derivatives for temp interpolation

do i=1,nlinesJ !Read in the variables
  read(100,*) AltJ(i),PressJ(i),TempJ(i)
end do

do i=1,nlinesS !Read in the variables
  read(101,*) AltS(i),TempS(i) !"AltS" for Altitude-Stephen
end do

do i=1,nlinesSP !Read in the variables
  read(103,*) AltSP(i),dum,PressS(i) !"AltSP" for Altitude-Stephen Pressure
end do

close(100) !Close all the files
close(101)
close(103)

call spline(AltS,TempS,nlinesS,Temp2) !2nd derivative for spline interpolation

do i=1,nlinesSP !Interpolate the temperature of the upper atmosphere to match
  call splineinterp(AltSP(i),AltS,TempS,nlinesS,Temp2,TempSI(i)) !press alt bins
end do

allocate(AltI(nlinesSP+61)) !Allocate the arrays to combine both data
allocate(PressI(nlinesSP+61))
allocate(TempI(nlinesSP+61))
allocate(TempI2(nlinesSP+61)) !2nd derivatives
allocate(PressI2(nlinesSP+61)) !2nd derivatives
TempI2=0.0 !Initialize
PressI2=0.0
do i=1,nlinesSP+61 !Create alt, temp, and press arrays of upper and lower atm.
  if(i.lt.62)then
    AltI(i)=AltJ(i)
    PressI(i)=PressJ(i)
    TempI(i)=TempJ(i)
  elseif(i.ge.62)then
    AltI(i)=AltSP(i-61)
    PressI(i)=PressS(i-61)
    TempI(i)=TempSI(i-61)
  end if
end do
do i=1,1544
  AltIS(i)=-90+2*i !Create the altitude bins to be interpolated
end do
call spline(AltI,TempI,1461,TempI2) !2nd derivative calcuations
call spline(AltI,PressI,1461,PressI2)
do i=1,1544 !Interpolate the temperature and pressure
  call splineinterp(AltIS(i),AltI,TempI,nlinesSP+61,TempI2,TempIS(i))
  call splineinterp(AltIS(i),AltI,PressI,nlinesSP+61,PressI2,PressIS(i))
end do
do i=1,1544 !Write out press and temp interpolation for a check, this should be
  write(102,*) AltIS(i),PressIS(i),TempIS(i) !unneeded if everything works.
end do
close(102) !Close that file
do i=1,1400 !Read in previous atmosphere density profile
  read(104,*) dum,nH2(i),nHe(i),nCH4(i),nH(i),nTot(i)
  mrH2(i)=nH2(i)/nTot(i) !Calculate the mixing ratios
  mrHe(i)=nHe(i)/nTot(i)
  mrCH4(i)=nCH4(i)/nTot(i)
  mrH(i)=nH(i)/nTot(i)
  mrTot(i)=mrH2(i)+mrHe(i)+mrCH4(i)+mrH(i)
end do
close(104) !Close that file
nTotI=PressIS/(TempIS*Kb) !Calculate total density array of previous atmosphere
do i=1,1544
  if(i.le.144)then
    nH2I(i)=mrH2(10)*nTotI(i) !Distribute the constituents through the lower
    nHeI(i)=mrHe(10)*nTotI(i) !atmosphere based on a mixing ratio at the
    nCH4I(i)=mrCH4(10)*nTotI(i) !bottom of the atmosphere, which is assumed to
    nHI(i)=mrH(10)*nTotI(i) !be in the homopause.
  else
    nH2I(i)=nH2(i-144) !Leave the upper atmosphere as how I previously had it
    nHeI(i)=nHe(i-144) !Juno was unable to measure the upper atmosphere temp.
    nCH4I(i)=nCH4(i-144) !at the creation of this program.
    nHI(i)=nH(i-144)
  end if
end do
PressIS=PressIS(1544:1:-1) !Flip the array for scale height calculation
do i=1,1544 !Calculate scale height
  if(i.lt.1544)H(i)=PressIS(i)/(-(PressIS(i)-PressIS(i+1))/(2.0e5))
  if(i.eq.1544)H(i)=PressIS(i)/(-(PressIS(i)-10306.7)/(2.0e5))
end do
PressIS=PressIS(1544:1:-1) !Flip the arrays back
H=H(1544:1:-1) !Flip the arrays back
write(105,*) '# Alt [km]  H2 [cm^-3]   He [cm^-3]  CH4 [cm^-3]   H [cm^-3]   &
              Tot [cm^-3]  Temp [K]  Press [mBar] Scale height [km]'
do i=1,1544
  write(105,10005) AltIS(i),nH2I(i),nHeI(i),nCH4I(i),nHI(i),nTotI(i),&
                   TempIS(i),PressIS(i),H(i)
end do
close(105) !Close that file
!Flip arrays again for the column density calculation
nH2I  = nH2I(1544:1:-1)
nHeI  = nHeI(1544:1:-1)
nCH4I = nCH4I(1544:1:-1)
nHI   = nHI(1544:1:-1)
write(106,*) '# Alt [km]  H2 [cm^-2]   He [cm^-2]  CH4 [cm^-2]   H [cm^-2]   &
              Tot [cm^-2]'
do i = 1544,1,-1 !From 2998 km to -88 km
  !For a given altitude (alt) compute the column density for each species using
  !COLDENS subroutine
  call COLDENS2(AltIS(i),nH2I,nHeI,nCH4I,nHI,CD_H2, CD_He, CD_CH4, CD_H,TotCD)
  !write results to output file
  write(106,10006) AltIS(i), CD_H2, CD_He, CD_CH4, CD_H, TotCD
end do
close(106) !Close that file
10000 format(F8.2,2(2x,ES11.3E2),2x,ES11.3E2)
10005 format(1x,F8.2,2X,5(ES11.3E2,2X),F8.2,2x,ES11.3E2,5x,ES11.3E2)
10006 format(1x,F8.2,2x,5(ES11.3E2,2X))

end program
