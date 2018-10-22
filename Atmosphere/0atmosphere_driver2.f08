subroutine atmos_driver!(TotalCD,Altitude,Delta)
!*******************************************************************************
!* Created by Stephen J. Houston 2.20.18 from previous work
!*******************************************************************************
!************************************************************************
!	  The purpose of this driver is to find the density of the atmosphere
!   in Jupiter as a function of altitude and to calculate the column density
!***********************************************************************
! Returns:
!	  species 1
!     nH2 --> H2 density; Type: real*8; Units: cm^-3
!     CD_H2 --> H2 column density; Type: real*8; Units: cm^-2

!	  species 2
!     nHe --> He density; Type: real*8; Units: cm^-3
!	  CD_He --> He column density; Type: real*8; Units: cm^-2

!	  species 3
!     nCH4 --> CH4 density; Type: real*8; Units: cm^-3
!     CD_CH4 --> CH4 column density; Type: real*8; Units: cm^-2

!	  species 4
!     nH --> H density; Type: real*8; Units: cm^-3
!     CD_H --> H column density; Type: real*8; Units: cm^-2

!	  TotCD--> Total column density; Type: real*8; Units: cm^-2

!************************************************************************

implicit none

!*****************INTERFACE BLOCK FOR ALLOCATABLE ARRAYS*****************
!You must add an interface block for any allocatable arrays that you want to
!pass from one routine/subroutine to another.

interface
  subroutine atmos_data(species,neutdens,alt,arrsize,y2)
    integer i, iostatus, N !N is one plus the size of the array. The array is N-1
    integer, intent(in) :: species
    integer, intent(out) :: arrsize
    real*8 yp1, ypn
    real*8 x, y
    real*8, allocatable,intent(out) :: neutdens(:), alt(:), y2(:)
  end subroutine
end interface

!*************************DECLARATION OF VARIABLES***********************


!real*8 NH2, NHE, NCH4, NH, N
integer i, species, lH2, lHe, lCH4, lH, atmoslen, j
real*8 CD_H2, CD_He, CD_CH4, CD_H,TotCD, alt,dz,altinit
parameter (dz = 0.1, atmoslen=2800*1/dz+1) !dz = size of altidude bins; atmoslen = number of altitude bins to cover the atmosphere
real*8, allocatable :: nH2dat(:),nHedat(:),nCH4dat(:),nHdat(:),altH2dat(:),&
altHedat(:), altCH4dat(:), altHdat(:),nH2y2(:), nHey2(:), nCH4y2(:), nHy2(:)
real*8 Z(atmoslen),nH2(atmoslen),nHE(atmoslen),nCH4(atmoslen),nH(atmoslen),&
TotalCD(atmoslen), altitude(atmoslen), delta(atmoslen), TotDens(atmoslen)
real*8 dum,P(1401),H(1400)

external atmosphere, coldens2

!*********************************MAIN PROGRAM***************************

open(unit=9, file='./Atmosphere/Input/atmosphere2km.dat', status='UNKNOWN') !output file with neutral densities as a function of altitude
open(unit=10, file='./Atmosphere/Input/columnDens2km.dat', status='UNKNOWN') !output file with column density as a function of altitude
open(unit=11, file='./AltDensPress.dat',status='old')
read(11,*)
P(1)=0.142
do i=2,1401
  read(11,*) dum,dum,P(i)
  !write(*,*) P(i)
end do
P=P(1401:1:-1)
do i=1,1400
  H(i)=P(i)/(-(P(i)-P(i+1))/(2.0e5)) !Scale height
!  write(*,*) P(i), H(i)
end do
H=H(1400:1:-1)
!*********************************ATMOS_DATA*****************************
! first calculate the neutral densities:
do species = 1,4 !for each different species
  !write(*,*) 'species driver=', species
  if (species .eq. 1) then !for H2
    call atmos_data(species,nH2dat,altH2dat,lH2,nH2y2)
	else if (species .eq. 2) then !for H2
    call atmos_data(species,nHedat,altHedat,lHe,nHey2)
	else if (species .eq. 3) then !for CH4
    call atmos_data(species,nCH4dat,altCH4dat,lCH4,nCH4y2)
	else if (species .eq. 4) then
    call atmos_data(species,nHdat,altHdat,lH,nHy2)
	end if
end do
!*********************************ATMOSPHERE*****************************

altinit = 200 !lowest altitude in the atmosphere = 200 km
do i=1,atmoslen !200, 3000
  Z(i) = ((i-1)*dz)+altinit  !find the altitude where you are calculating things
  call atmosphere(Z(i),altH2dat,nH2dat,lH2,nH2y2,nH2(i))
  call atmosphere(Z(i),altHedat,nHedat,lHe,nHey2,nHe(i))
  call atmosphere(Z(i),altCH4dat,nCH4dat,lCH4,nCH4y2,nCH4(i))
  call atmosphere(Z(i),altHdat,nHdat,lH,nHy2,nH(i))
  TotDens(i)=nH2(i)+nHE(i)+nCH4(i)+nH(i)
!  write(*,200) Z(i), nH2(i)
!  call atmosphere(Z(i),nH2dat,altH2dat,lH2,nHedat,altHedat,lHe,&
!  nCH4dat,altCH4dat,lCH4,nHdat,altHdat,lH,nH2(i), nHE(i),&
!  nCH4(i), nH(i))  !for the given altitude point Z(i) get the density for all the species
	if(i.eq.1)then
    j=1
		write(9, 101) Z(i), nH2(i), nHE(i), nCH4(i), nH(i), TotDens(i),H(j)
	else if(mod(i-1,20).eq.0)then !write every 10 km, but you can change this
    j=j+1
		write(9, 101) Z(i), nH2(i), nHE(i), nCH4(i), nH(i), TotDens(i),H(j)
  end if
end do

!*******************************COLUMN DENSITY***************************

!FLIP ARRAYS
nH2  = nH2(atmoslen:1:-1)
nHe  = nHe(atmoslen:1:-1)
nCH4 = nCH4(atmoslen:1:-1)
nH   = nH(atmoslen:1:-1)
! now compute the column density (cm^{-2})
j=1
do i = 2998,200,-2 !From 2990 km (to match 2-stream code) to 500 km in increments of 10 km
  alt = i*1.0
  call COLDENS2(alt,nH2,nHe,nCH4,nH,CD_H2, CD_He, CD_CH4, CD_H,TotCD) !for a given altitude (alt) compute the column density for each species using COLDENS subroutine
  write(10, 102) alt, CD_H2, CD_He, CD_CH4, CD_H, TotCD !write results to output file
  TotalCD(j)=TotCD
  altitude(j)=alt
  delta(j)=2
!  write(*,*) TotalCD (j), altitude(j), delta(j), size(delta)
  j=j+1
end do
goto 100
j=j-1 !This is here so the 500 km altitude isn't counted twice
do i = 500,400,-2 !From 500 km to 400 km in increments of 5 km
  alt = i*1.0
  call COLDENS2(alt,nH2,nHe,nCH4,nH,CD_H2, CD_He, CD_CH4, CD_H,TotCD) !for a given altitude (alt) compute the column density for each species using COLDENS subroutine
!  write(10, 102) alt, CD_H2, CD_He, CD_CH4, CD_H, TotCD !write results to output file
  TotalCD(j)=TotCD
  altitude(j)=alt
  delta(j)=2
  j=j+1
end do
j=j-1 !This is here so the 400 km altitude isn't counted twice
do i = 400,200,-2 !From 400 km to 200 km in increments of 2 km
  alt = i*1.0
  call COLDENS2(alt,nH2,nHe,nCH4,nH,CD_H2, CD_He, CD_CH4, CD_H,TotCD) !for a given altitude (alt) compute the column density for each species using COLDENS subroutine
!  write(10, 102) alt, CD_H2, CD_He, CD_CH4, CD_H, TotCD !write results to output file
  TotalCD(j)=TotCD
  altitude(j)=alt
  delta(j)=2
  j=j+1
end do
j=j-1
!write(*,*) '3'
!write(*,*) j
!do i =1,j
!  write(*,*) TotalCD(i), altitude(i), delta(i)
!end do
101  FORMAT(1x,7(E11.4,1X))
102  FORMAT(1x,6(E11.4,1X))
200  FORMAT(3e12.4)
100 continue

end subroutine
