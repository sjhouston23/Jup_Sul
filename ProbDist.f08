program ProbDist
!*******************************************************************************
!* Created by Stephen J. Houston 10.23.18
!*******************************************************************************
!* Interpolates the singly differential ejected electron energy/angle
!* cross-section data with a loglog linear interpolator.
!* Number of interpolation points is determined by the number of points you
!* want between each data point, using "nInterp"
!* Uses data from Schultz et al., 2019 to determine the probability
!* distribution functions for energy (2 stream bins) and angle (0-180).
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

integer nEnergies,nProc,neEnergies,nInterpEng,neAngles,nInterpAng,nChS,iEng,ionEng,i,j
integer Proc,ChS,Eng,Ang
integer SI,DI,TI,DA,SS,DS !Electron producing processes
real*8 pi

parameter(nEnergies=9) !Number of inital ion energies
parameter(nProc=6) !Number of electron ejection processes
parameter(neEnergies=85) !Number of electron ejection energies (0.005 - 6525 eV)
parameter(nInterpEng=nInterp*(neEnergies-1)+1) !Number of interpolated energies
parameter(neAngles=46) !Number of electron ejection angles (0.05° - 179.5°)
parameter(nInterpAng=nInterp*(neAngles-1)+1) !Number of interpolated angles
parameter(nChS=17) !Number of charge states from 0-16
parameter(SI=1,DI=2,TI=3,DA=4,SS=5,DS=6) !Electron producing processes
parameter(pi=4.0*atan(1.0d0)) !pi to convert degrees to radians
!* Schultz et al., 2019 data
real*8,dimension(neEnergies) :: eEnergy !The electron ejection energies
real*8,dimension(neAngles) :: eAngle !The electron ejection agles
real*8,dimension(nProc,nChS,nEnergies,neEnergies) :: SDXSe !SDXS energy
real*8,dimension(nProc,nChS,nEnergies,neAngles) :: SDXSa !SDXS angle
!* Interpolated data
real*8,dimension(nInterpEng) :: InterpEnergy
real*8,dimension(nInterpAng) :: InterpAngle,InterpAngleRad
real*8,dimension(nProc,nChS,nEnergies,nInterpEng) :: InterpSDXSe
real*8,dimension(nProc,nChS,nEnergies,nInterpAng) :: InterpSDXSa
!******************* Get singly differential cross-sections ********************
open(unit=101,file='./SDXSeall.dat')
open(unit=102,file='./SDXSaall.dat')
read(101,1000) eEnergy !Electron energies
read(102,1000) eAngle !Electron angles
read(101,1001) SDXSe !Singly differential cross-section data for electron energy
read(102,1001) SDXSa !Singly differential cross-section data for electron angle
close(101)
close(102)
1000 format(10(F8.3,1x))
1001 format(10(ES9.3E2,1x))
!**************************** Initialize Variables *****************************
InterpSDXSe=0.0;InterpSDXSa=0.0;InterpEnergy=0.0;InterpAngle=0.0
InterpAngleRad=0.0
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
do Proc=1,1!nProc !Loop through every process
  do ChS=1,1!nChS !Loop through every charge state
    do ionEng=1,1!nEnergies !Loop through every initial ion energy
      Eng=1
      do iEng=1,31!nInterpEng !Loop through all interpolated energies
        if(InterpEnergy(iEng).ge.eEnergy(Eng+1)) Eng=Eng+1 !Go to next Energy when appropriate
        if(iEng.eq.nInterpEng) Eng=neEnergies!Don't want Eng to go out of bounds
        InterpSDXSe(Proc,ChS,ionEng,iEng)=log(SDXSe(Proc,ChS,ionEng,Eng))+&
        (log(InterpEnergy(iEng))-log(eEnergy(Eng+1)))*&
        (log(SDXSe(Proc,ChS,ionEng,Eng+1))-log(SDXSe(Proc,ChS,ionEng,Eng)))/&
        (log(eEnergy(Eng+1))-log(eEnergy(Eng)))
        write(*,"(F8.3,1x,ES9.3E2,2(1x,F8.3),2(1x,ES9.3E2))") InterpEnergy(iEng),&
        exp(InterpSDXSe(Proc,ChS,ionEng,iEng)),&
        eEnergy(Eng+1),eEnergy(Eng),SDXSe(Proc,ChS,ionEng,Eng+1),&
        SDXSe(Proc,ChS,ionEng,Eng)
        ! write(*,"(2(ES12.3E2,1x,F8.3,1x))") log(SDXSe(Proc,ChS,ionEng,Eng)),&
        ! (log(InterpEnergy(iEng))-log(eEnergy(Eng+1))),&
        ! (log(SDXSe(Proc,ChS,ionEng,Eng+1))-log(SDXSe(Proc,ChS,ionEng,Eng))),&
        ! (log(eEnergy(Eng+1))-log(eEnergy(Eng)))
      end do !End of interpolation do-loop
    end do !End of ion energy do-loop
  end do !End of charge state do-loop
end do !End of processes do-loop
! do i=1,neEnergies
!   write(*,"(F8.3,1x,ES9.3E2)") eEnergy(i),SDXSe(1,1,1,i)
! end do
! do i=1,31!nInterpEng
!   write(*,"(F8.3,1x,ES9.3E2)") InterpEnergy(i),exp(InterpSDXSe(1,1,1,i))
! end do








end program
