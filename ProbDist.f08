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

integer,dimension(nEnergies) :: Energy !Each initial energy

real*8 integralSDXSe,integralSDXSa,integralSDXSei,integralSDXSai,EsXS,AsXS
!* Schultz et al., 2019 data
real*8,dimension(neEnergies) :: eEnergy,dE,tmpe !The electron ejection energies
real*8,dimension(neAngles) :: eAngle,dA,tmpa,eAngleRad !The electron angles
real*8,dimension(nProc,nChS,nEnergies,neEnergies) :: SDXSe !SDXS energy
real*8,dimension(nProc,nChS,nEnergies,neAngles) :: SDXSa !SDXS angle
real*8,dimension(nProc,nChS,nEnergies) :: NSIMxs !NSIM cross-sections
!* Interpolated data
real*8,dimension(nInterpEng) :: InterpEnergy,dEi,tmpei
real*8,dimension(nInterpAng) :: InterpAngle,InterpAngleRad,dAi,tmpai
real*8,dimension(nProc,nChS,nEnergies,nInterpEng) :: InterpSDXSe,eProbFunc
real*8,dimension(nProc,nChS,nEnergies,nInterpAng) :: InterpSDXSa,aProbFunc

real*8,allocatable,dimension(:) :: dEtmp,dAtmp,temp

character(len=4),dimension(nProc) :: Processes
!****************************** Data Declaration *******************************
data Energy/1,10,50,75,100,200,500,1000,2000/
data Processes/'  SI','  DI','  TI','DCAI','  SS','  DS'/
!******************* Get singly differential cross-sections ********************
! open(unit=101,file='./Electron_Dist/SDXSeall.dat') !0.0 values, change 1.000E-30
! open(unit=102,file='./Electron_Dist/SDXSaall.dat') !0.0 values, change 1.000E-30
open(unit=101,file='../SulfurXS/SDXSeall.dat')
open(unit=102,file='../SulfurXS/SDXSaall.dat')
read(101,1000) eEnergy !Electron energies
read(102,1000) eAngle !Electron angles
read(101,1001) SDXSe !Singly differential cross-section data for electron energy
read(102,1001) SDXSa !Singly differential cross-section data for electron angle
close(101)
close(102)
1000 format(10(F8.3,1x))
1001 format(10(ES9.3E2,1x))
!******************************* Read NSIM data ********************************
open(unit=104,file='NSIMall.txt',status='old') !All NSIM data
do proc=1,nProc
  read(104,*)
  read(104,*)
  do Eng=1,nEnergies
    read(104,1002) (NSIMxs(Proc,ChS,Eng),ChS=1,nChS)
  end do
end do
1002 format(7x,17(2x,ES8.2E2))
!**************************** Initialize Variables *****************************
InterpSDXSe=0.0;InterpSDXSa=0.0;InterpEnergy=0.0;InterpAngle=0.0
InterpAngleRad=0.0;eAngleRad=0.0;eProbFunc=0.0;aProbFunc=0.0
integralSDXSe=0.0;integralSDXSa=0.0;integralSDXSei=0.0;integralSDXSai=0.0
!***** Calculate the change in energies and angles for Simpson's rule
do Eng=1,neEnergies
  if(Eng.lt.neEnergies)then
    dE(Eng)=eEnergy(Eng+1)-eEnergy(Eng)
  else
    dE(Eng)=6725.0-eEnergy(Eng)
  end if
end do
eAngleRad=eAngle*pi/180.0 !Convert angles to radians
do Ang=1,neAngles
  if(Ang.lt.neAngles)then
    dA(Ang)=eAngleRad(Ang+1)-eAngleRad(Ang) !Want dA in radians
  else
    dA(Ang)=pi-eAngleRad(Ang)
  end if
end do
!******************* Create interpolated energies and angles *******************
j=0
do i=1,nInterpEng !Calculate interpolated energy bins
  if(mod(i-1,nInterp).eq.0)j=j+1
  if(j.eq.neEnergies)then
    InterpEnergy(i)=eEnergy(j)+mod(i-1,nInterp)*(eEnergy(j)-eEnergy(j))/nInterp
    goto 10
  end if
  InterpEnergy(i)=eEnergy(j)+mod(i-1,nInterp)*(eEnergy(j+1)-eEnergy(j))/nInterp
  10 continue
end do
j=0
do i=1,nInterpAng !Calculate interpolated angle bins
  if(mod(i-1,nInterp).eq.0)j=j+1
  if(j.eq.neAngles)then
    InterpAngle(i)=eAngle(j)+mod(i-1,nInterp)*(eAngle(j)-eAngle(j))/nInterp
    goto 20
  end if
  InterpAngle(i)=eAngle(j)+mod(i-1,nInterp)*(eAngle(j+1)-eAngle(j))/nInterp
  20 continue
end do
!***** Calculate the change in energies and angles for Simpson's rule
do Eng=1,nInterpEng
  if(Eng.lt.nInterpEng)then
    dEi(Eng)=InterpEnergy(Eng+1)-InterpEnergy(Eng)
  else
    dEi(Eng)=dEi(Eng-1)
  end if
end do
InterpAngleRad=InterpAngle*pi/180.0 !Convert angles to radians
do Ang=1,nInterpAng
  if(Ang.lt.nInterpAng)then
    dAi(Ang)=InterpAngleRad(Ang+1)-InterpAngleRad(Ang) !Want dA in radians
  else
    dAi(Ang)=pi-InterpAngleRad(Ang)
  end if
end do
!******************************** Main Program *********************************
do Proc=1,nProc !Loop through every process
  write(*,*) '---------- NEW PROCESS ----------'
  do ChS=1,nChS !Loop through every charge state
    write(*,*) '---------- NEW CHARGE STATE ----------'
    write(*,*) " Proc   ChS Energy Simpson's    Interp      NSIM  Ratio Ratioi"
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
      end do
      do Ang=1,neAngles
        tmpa(Ang)=SDXSa(Proc,ChS,ionEng,Ang)*sin(eAngleRad(Ang))*2*pi
      end do
      do Eng=1,nInterpEng
        tmpei(Eng)=exp(InterpSDXSe(Proc,ChS,ionEng,Eng)) !Create array for Simp
      end do
      do Ang=1,nInterpAng
        tmpai(Ang)=exp(InterpSDXSa(Proc,ChS,ionEng,Ang))*&
                   sin(InterpAngleRad(Ang))*2*pi
      end do
      call simp(neEnergies,dE,tmpe,integralSDXSe)
      call simp(neAngles,dA,tmpa,integralSDXSa)
      call simp(nInterpEng,dEi,tmpei,integralSDXSei)
      call simp(nInterpAng,dAi,tmpai,integralSDXSai)
      if(integralSDXSei/NSIMxs(Proc,ChS,ionEng).lt.0.8.or.&
      integralSDXSei/NSIMxs(Proc,ChS,ionEng).gt.1.20)&
        write(*,'(2x,A4,"-E",I4,2x,I5,3(2x,ES8.2E2),2(1x,F6.2))') &
        Processes(Proc),ChS-1,Energy(ionEng),integralSDXSe,integralSDXSei,&
        NSIMxs(Proc,ChS,ionEng),integralSDXSe/NSIMxs(Proc,ChS,ionEng),&
        integralSDXSei/NSIMxs(Proc,ChS,ionEng)
      if(integralSDXSai/NSIMxs(Proc,ChS,ionEng).lt.0.8.or.&
      integralSDXSai/NSIMxs(Proc,ChS,ionEng).gt.1.20)&
        write(*,'(2x,A4,"-A",I4,2x,I5,3(2x,ES8.2E2),2(1x,F6.2))')&
        ! write(*,'(2x,A4,"-A",I4,2x,I5,2x,ES8.2E2,2x,ES8.2E2,11x,F6.2,1x,F6.2)')&
        Processes(Proc),ChS-1,Energy(ionEng),integralSDXSa,integralSDXSai,&
        NSIMxs(Proc,ChS,ionEng),integralSDXSa/NSIMxs(Proc,ChS,ionEng),&
        integralSDXSai/NSIMxs(Proc,ChS,ionEng)
!************************** Probability Distribution ***************************
      do iEng=1,nInterpEng !Loop through all interpolated energies
        allocate(dEtmp(iEng),temp(iEng))
        dEtmp=0.0;temp=0.0
        do i=1,iEng
          dEtmp(i)=dEi(i) !For Simpson's rule
          temp(i)=tmpei(i) !For Simpson's rule
        end do
        call simp(iEng,dEtmp,temp,EsXS) !Integrate from 1 to Es
        eProbFunc(Proc,ChS,ionEng,iEng)=1-(EsXS/integralSDXSei) !Prob dist func
        if(eProbFunc(Proc,ChS,ionEng,iEng).lt.0)then !Issue with going negative
          eProbFunc(Proc,ChS,ionEng,iEng)=-eProbFunc(Proc,ChS,ionEng,iEng)
        endif
        deallocate(dEtmp,temp)
      end do
      do iAng=1,nInterpAng !Loop through all interpolated energies
        allocate(dAtmp(iAng),temp(iAng))
        dAtmp=0.0;temp=0.0
        do i=1,iAng
          dAtmp(i)=dAi(i) !For Simpson's rule
          temp(i)=tmpai(i) !For Simpson's rule
        end do
        call simp(iAng,dAtmp,temp,AsXS) !Integrate from 1 to Es
        aProbFunc(Proc,ChS,ionEng,iAng)=1-(AsXS/integralSDXSai) !Prob dist func
        if(aProbFunc(Proc,ChS,ionEng,iAng).lt.0)then !Issue with going negative
          aProbFunc(Proc,ChS,ionEng,iAng)=-aProbFunc(Proc,ChS,ionEng,iAng)
        endif
        deallocate(dAtmp,temp)
      end do
    end do !End of ion energy do-loop
  end do !End of charge state do-loop
end do !End of processes do-loop

open(unit=200,file='./Electron_Dist/eProbFunc.dat') !Electron energy prob dist
open(unit=201,file='./Electron_Dist/aProbFunc.dat') !Electron angle prob dist
write(200,*) nInterpEng !Write out number of energy interpolations
write(201,*) nInterpAng !Write out number of angle interpolations
write(200,1000) InterpEnergy
write(201,1000) InterpAngle
write(200,1001) eProbFunc
write(201,1001) aProbFunc
close(200)
close(201)

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
