program ProbDist
!*******************************************************************************
!* Created by Stephen J. Houston 10.23.18
!*******************************************************************************
!* Interpolates the singly differential ejected electron energy/angle
!* cross-section data with a loglog linear interpolator.
!* Number of interpolation points is determined by the number of points you
!* want between each data point, using "nInterp"
!* Uses data from Schultz et al., 2019 to determine the probability
!* distribution functions for energy and angle (0-180).
!*
!* Electron Ejecting Collision Types:
!*    1) Single Ionization
!*    2) Double Ionization
!*    3) Transfer Ionization
!*    4) Double-Capture Auto Ionization
!*    5) Single Stripping
!*    6) Double stripping
!*
!* Possible Ion Energies:
!*    1, 10, 50, 75, 100, 200, 500, 1000, 2000
!*
!*******************************************************************************

implicit none!real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer nInterp
parameter(nInterp=10) !Number of interpolation points between each data point

integer nEnergies,nProc,neEnergies,nInterpEng,neAngles,nInterpAng,nChS,iEng
integer iAng,ionEng,i,j
integer Proc,ChS,Eng,Ang
integer SI,DI,TI,DCAI,SS,DS !Electron producing processes
real*8 pi

parameter(nEnergies=9) !Number of inital ion energies
parameter(nProc=6) !Number of electron ejection processes
parameter(neEnergies=85) !Number of electron ejection energies (0.005 - 6525 eV)
parameter(nInterpEng=nInterp*(neEnergies-1)+1) !Number of interpolated energies
parameter(neAngles=46) !Number of electron ejection angles (0.05° - 179.5°)
parameter(nInterpAng=nInterp*(neAngles-1)+1) !Number of interpolated angles
parameter(nChS=17) !Number of charge states from 0-16
parameter(SI=1,DI=2,TI=3,DCAI=4,SS=5,DS=6) !Electron producing processes
parameter(pi=4.0*atan(1.0d0)) !pi to convert degrees to radians

real*8 integralSDXSe,integralSDXSa
!* Schultz et al., 2019 data
real*8,dimension(neEnergies) :: eEnergy,dE,tmpe !The electron ejection energies
real*8,dimension(neAngles) :: eAngle,dA,tmpa,eAngleRad !The electron angles
real*8,dimension(nProc,nChS,nEnergies,neEnergies) :: SDXSe !SDXS energy
real*8,dimension(nProc,nChS,nEnergies,neAngles) :: SDXSa !SDXS angle
real*8,dimension(nProc,nChS,nEnergies) :: NSIMxs !NSIM cross-sections
!* Interpolated data
real*8,dimension(nInterpEng) :: InterpEnergy
real*8,dimension(nInterpAng) :: InterpAngle,InterpAngleRad
real*8,dimension(nProc,nChS,nEnergies,nInterpEng) :: InterpSDXSe
real*8,dimension(nProc,nChS,nEnergies,nInterpAng) :: InterpSDXSa

character(len=4),dimension(nProc) :: Processes
!****************************** Data Declaration *******************************
data Processes/'  SI','  DI','  TI','DCAI','  SS','  DS'/
!******************* Get singly differential cross-sections ********************
open(unit=101,file='./SDXSeall.dat') !If 0.0 values, need to change to 1.000E-30
open(unit=102,file='./SDXSaall.dat') !If 0.0 values, need to change to 1.000E-30
read(101,1000) eEnergy !Electron energies
read(102,1000) eAngle !Electron angles
read(101,1001) SDXSe !Singly differential cross-section data for electron energy
read(102,1001) SDXSa !Singly differential cross-section data for electron angle
close(101)
close(102)
1000 format(10(F8.3,1x))
1001 format(10(ES9.3E2,1x))
!******************************* Read NSIM data ********************************
open(unit=104,file='NSIM.txt',status='old') !All NSIM data
do ChS=1,nChS !Loop through every charge state
  do i=1,20
    read(104,*) !Skip a bunch of blank lines
  end do
  do Eng=1,nEnergies !Loop through every energy
    read(104,*) !Skip a line
    do Proc=1,nProc
      read(104,*) NSIMxs(Proc,ChS,Eng)
    end do
    do i=1,3
      read(104,*)
    end do
  end do
end do
close(104)
!**************************** Initialize Variables *****************************
InterpSDXSe=0.0;InterpSDXSa=0.0;InterpEnergy=0.0;InterpAngle=0.0
InterpAngleRad=0.0;eAngleRad=0.0
!***** Calculate the change in energies and angles for Simpson's rule
do Eng=1,neEnergies
  if(Eng.eq.1)then
    dE(Eng)=eEnergy(Eng)-0.0
  else
    dE(Eng)=eEnergy(Eng)-eEnergy(Eng-1)
  end if
end do
eAngleRad=eAngle*pi/180.0 !Convert angles to radians
do Ang=1,neAngles
  if(Eng.eq.1)then
    dA(Ang)=eAngleRad(Ang)-0.0 !Want dA in radians
  else
    dA(Ang)=eAngleRad(Ang)-eAngleRad(Ang-1)
  end if
end do
!******************* Create interpolated energies and angles *******************
j=0
do i=1,nInterpEng !Calculate interpolated energy bins
  if(mod(i-1,nInterp).eq.0)j=j+1
  InterpEnergy(i)=eEnergy(j)+mod(i-1,nInterp)*(eEnergy(j+1)-eEnergy(j))/nInterp
end do
j=0
do i=1,nInterpAng !Calculate interpolated angle bins
  if(mod(i-1,nInterp).eq.0)j=j+1
  InterpAngle(i)=eAngle(j)+mod(i-1,nInterp)*(eAngle(j+1)-eAngle(j))/nInterp
end do
!******************************** Main Program *********************************
do Proc=1,nProc !Loop through every process
  write(*,*) '---------- NEW PROCESS ----------'
  do ChS=1,nChS !Loop through every charge state
    write(*,*) '---------- NEW CHARGE STATE ----------'
    write(*,*) " Proc   ChS Energy Simpson's      NSIM  Ratio"
    do ionEng=1,nEnergies !Loop through every initial ion energy
      Eng=1 !Reset index
      do iEng=1,nInterpEng !Loop through all interpolated energies
        if(InterpEnergy(iEng).gt.eEnergy(Eng+1)) Eng=Eng+1 !Go to next Energy
        if(iEng.eq.nInterpEng) Eng=neEnergies-1 !Don't want to go out of bounds
        InterpSDXSe(Proc,ChS,ionEng,iEng)=log(SDXSe(Proc,ChS,ionEng,Eng))+&
          (log(InterpEnergy(iEng))-log(eEnergy(Eng)))*&
          (log(SDXSe(Proc,ChS,ionEng,Eng+1))-log(SDXSe(Proc,ChS,ionEng,Eng)))/&
          (log(eEnergy(Eng+1))-log(eEnergy(Eng)))
        ! write(*,"(F8.3,1x,ES9.3E2,2(1x,F8.3),2(1x,ES9.3E2))") InterpEnergy(iEng),&
        ! exp(InterpSDXSe(Proc,ChS,ionEng,iEng)),&
        ! eEnergy(Eng+1),eEnergy(Eng),SDXSe(Proc,ChS,ionEng,Eng+1),&
        ! SDXSe(Proc,ChS,ionEng,Eng)
      end do !End of energy interpolation do-loop
      Ang=1 !Reset index
      do iAng=1,nInterpAng !Loop through all interpolated angles
        if(InterpAngle(iAng).gt.eAngle(Ang+1)) Ang=Ang+1 !Go to next Angle
        if(iAng.eq.nInterpAng) Ang=neAngles-1 !Don't want Ang to go out of bounds
        InterpSDXSa(Proc,ChS,ionEng,iAng)=log(SDXSa(Proc,ChS,ionEng,Ang))+&
          (log(InterpAngle(iAng))-log(eAngle(Ang)))*&
          (log(SDXSa(Proc,ChS,ionEng,Ang+1))-log(SDXSa(Proc,ChS,ionEng,Ang)))/&
          (log(eAngle(Ang+1))-log(eAngle(Ang)))
        ! write(*,"(F8.3,1x,ES9.3E2,2(1x,F8.3),2(1x,ES9.3E2))") InterpAngle(iAng),&
        ! exp(InterpSDXSa(Proc,ChS,ionEng,iAng)),&
        ! eAngle(Ang+1),eAngle(Ang),SDXSa(Proc,ChS,ionEng,Ang+1),&
        ! SDXSa(Proc,ChS,ionEng,Ang)
      end do
      do Eng=1,neEnergies
        tmpe(Eng)=SDXSe(Proc,ChS,ionEng,Eng) !Create a single array for Simp
        ! write(*,'(F10.3,2x,ES10.2E2,2x,F10.4,2x,ES10.2E2)')&
        ! eEnergy(Eng),tmpe(Eng),dE(Eng),dE(Eng)*tmpe(Eng)
      end do
      do Ang=1,neAngles
        tmpa(Ang)=SDXSa(Proc,ChS,ionEng,Ang)*sin(eAngleRad(Ang))!*2*pi
      end do
      call simp(neEnergies,dE,tmpe,integralSDXSe)
      call simp(neAngles,dA,tmpa,integralSDXSa)
      ! write(*,'(2x,A4,1x,I5,1x,1x,I5,2x,ES8.2E2,2x,ES8.2,1x,F6.2)') &
      ! Processes(Proc),ChS-1,ionEng,integralSDXSe,NSIMxs(Proc,ChS,ionEng),&
      ! integralSDXSe/NSIMxs(Proc,ChS,ionEng)
      write(*,'(2x,A4,1x,I5,1x,1x,I5,2x,ES8.2E2,2x,ES8.2,1x,F6.2)') &
      Processes(Proc),ChS-1,ionEng,integralSDXSa*2*pi,NSIMxs(Proc,ChS,ionEng),&
      integralSDXSa*2*pi/NSIMxs(Proc,ChS,ionEng)
    end do !End of ion energy do-loop
  end do !End of charge state do-loop
end do !End of processes do-loop









end program

!******************************* SIMPSON'S RULE ********************************
subroutine simp (N,H,F,INTEG)
!     integration via simpson's rule
!     N: number of intervals for the integration
!     H: integration interval dx (can't be variable)
!     F: vector with y values to be integrated
!     INTEG: resultant integral
!     This subroutine uses Simpson's rule to integrate f(x)

implicit none

integer N,i
real*8 F(N), F0, F2(N/2), H(N), INTEG, SF1,SF2
real*8,allocatable,dimension(:) :: F1

if(mod(N,2).eq.0)then
  allocate(F1(N/2-1))
else
  allocate(F1(N/2))
end if
F0 = 0.0
F1 = 0.0 !SUM OF ODD F(I)
F2 = 0.0 !SUM OF EVEN F(I)
SF1 = 0.0
SF2 = 0.0
!!!write(*,*) F
F1=F(3:N:2)
F2=F(2:N:2)
do i=1,size(F1)
  SF1=SF1+F1(i)*H((i*2)+1)
end do
do i=1,size(F2)
  SF2=SF2+F2(i)*H(i*2)
end do
!SF1=sum(F1) !Summation of all odd elements
!SF2=sum(F2) !Summation of all even elements
F0 = H(1)*F(1) + H(N)*F(N)

INTEG = (F0 + 2*SF2 + 4*SF1)/3
deallocate(F1)

return
end
