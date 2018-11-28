subroutine EjectedElectron(E,Proc,ChS,electron_energy,electron_angle,eBin)
!*******************************************************************************
!* Created by Stephen J. Houston 10.31.18
!*******************************************************************************
!* This subroutine will calculate the energy and angle of an ejected
!* electron for a given ion energy, ion charge state, & collision type.
!* Taking the probability distribution calculated (see ProbDist.f08)
!* using the singly differential cross sections as a function of the
!* ejected electron energy or angle, a monte carlo method similar to
!* that used in CollisionSim.f08 is used to calculate the ejected
!* electron's energy and angle.
!*******************************************************************************
!*
!* 	Input:
!*		E --> Energy of ions
!*		Type: Real
!*		Units: keV/u
!*
!*		Proc --> Process
!*		Type: Integer
!*		Units: None
!*
!*		ChS --> Charge state
!*		Type: Integer
!*		Units: None
!*
!*  Returns:
!*    electron_energy --> energy of ejected electron
!*    Type: Real
!*    Units: eV
!*
!*    electron_angle --> angle of ejected electron
!*    Type: Real
!*    Units: Degrees (0-180)
!*
!*    eBin --> 2-Stream electron energy bin
!*    Type: integer
!*    Units: None
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer,intent(in) :: Proc,ChS !Process and charge state
integer,intent(out) :: eBin !Electron energy bin in 2-stream format
real*8,intent(in) :: E !Ion energy
real*8,intent(out) :: electron_energy,electron_angle

parameter(neProc=6) !Number of processes that eject electrons
parameter(nChS=17) !Number of charge states from 0-16
parameter(nEnergies=9) !Number of inital energies
parameter(nE2strBins=260) !Number of 2 stream bins

integer Eng,eEng,eAng
integer iBin !IonEnergy bin
save :: neEnergies,neAngles
real,dimension(2) :: ranVecC
real*8,dimension(nEnergies) :: IonEnergy
real*8,allocatable,dimension(:),save :: eEnergy,eAngle
real*8,allocatable,dimension(:,:,:,:),save :: eProbFunc,aProbFunc
real*8,save :: es,del(nE2strBins),E2str(nE2strBins) !2-Stream energy bins
!****************************** Data Declaration *******************************
!* Initial ion energy input:
data IonEnergy/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0/
!dE for each 2-Stream energy bin. Must match two stream code binning
data del/20*0.5,70*1.0,10*2.0,20*5.0,10*10.0,20*10.0,10*50.0,10*100.0,40*200.0,&
         10*400,10*1000,10*2000,10*5000,10*10000.0/
!**************************** 2-Stream Bin Creation ****************************
if(E2str(nE2strBins).gt.10000.0)goto 1001
!2-Stream energy bins:
do i=1,nE2strBins
  es=es+del(i)
  E2str(i)=Es
end do
1001 continue
!*********************** Ejected Electron Probabilities ************************
if(neEnergies.gt.100)goto 1000 !Only want to read the files once
!* Created by ReadElectDist.f08 and ProbDist.f08
open(unit=100,file='./Electron_Dist/eProbFunc.dat',status='old')
open(unit=101,file='./Electron_Dist/aProbFunc.dat',status='old')
read(100,*) neEnergies !Number of ejected electron energies
read(101,*) neAngles !Number of ejected electron angles
allocate(eEnergy(neEnergies)) !Allocate electron energy array
allocate(eAngle(neAngles)) !Allocate electron angle array
allocate(eProbFunc(neProc,nChS,nEnergies,neEnergies)) !Electron energy prob func
allocate(aProbFunc(neProc,nChS,nEnergies,neAngles)) !Electron angle prob func
read(100,10000) eEnergy !Electron energy array
read(101,10000) eAngle !Electron angle array
read(100,10001) eProbFunc !Electron energy probability distribution function
read(101,10001) aProbFunc !Electron angle probability distribution function
close(100) !Close eProbFunc file
close(101) !Close aProbFunc file
10000 format(10(F8.3,1x)) !Formatting for energies and angles (ProbDist.f08)
10001 format(10(ES9.3E2,1x)) !Formatting for probabilities (ProbDist.f08)
1000 continue
!******************************** Main Program *********************************
!* Initialize:
iBin=0;k=0;f=0.0;ranVecC=0.0;eBin=0
do Eng=nEnergies,1,-1 !Loop through bins to get the correct ion energy bin
  if(E.ge.real(IonEnergy(Eng-1)+(IonEnergy(Eng)-IonEnergy(Eng-1))/2.0))then
    iBin = Eng !Get the ion energy bin number
    k = 1 !Used to get an interpolation function
    goto 2000
  elseif(E.ge.real(IonEnergy(Eng-1)))then
    iBin = Eng-1 !Get the ion energy bin number
    k = 2 !Used to get an interpolation function
    goto 2000
  endif
end do
2000 continue
if (iBin.eq.0) write(206,*) 'EjectedElectron.f08: Error! iBin=0'
!* Want to use f to somewhat interpolate the cross-section for ion energies that
!* lie between energy bins.
if (k.eq.1) f=(E-IonEnergy(iBin-1))/(IonEnergy(iBin)-IonEnergy(iBin-1))
if (k.eq.2) f=(E-IonEnergy(iBin))/(IonEnergy(iBin+1)-IonEnergy(iBin))
! write(206,*) E,IonEnergy(iBin),iBin,k,f

call ranlux(ranVecC,2) !Random number vector
!* Initialize:
electron_energy = 0.0 !Initialize ejected electron energy
electron_energytmp1 = 0.0 !Will calculate two different electron energies for
electron_energytmp2 = 0.0 ! a simple linear interpolation
electron_angle = 0.0 !Initialize ejected electron angle
electron_angletmp1 = 0.0 !Will calculate two different electron angles for
electron_angletmp2 = 0.0 ! a simple linear interpolation

do eEng=1,neEnergies !Loop through all of the electron ejection energies
  if(ranVecC(1).ge.eProbFunc(Proc,ChS,iBin,eEng))then !Random number vs. prob
    electron_energytmp1=eEnergy(eEng) !Get the corresponding electron energy
    goto 3000 !Get out of do loop
  end if !End of Monte Carlo if statement
end do !End first electron energy sample do loop
3000 continue !Continue on to find the next energy sample

do eEng=1,neEnergies !Loop through all of the electron ejection energies
  if(f.lt.0.5)then !Need to know which ion energy bin we're interpolating with
    if(ranVecC(1).ge.eProbFunc(Proc,ChS,iBin+1,eEng))then
      electron_energytmp2=eEnergy(eEng)
      electron_energy=(1-f)*electron_energytmp1+f*electron_energytmp2
      goto 4000 !Get out of do loop
    end if
  elseif(f.ge.0.5)then
    if(ranVecC(1).ge.eProbFunc(Proc,ChS,iBin-1,eEng))then
      electron_energytmp2=eEnergy(eEng)
      electron_energy=(1-f)*electron_energytmp2+f*electron_energytmp1
      goto 4000 !Get out of do loop
    end if
  end if !End of Monte Carlo if statement
end do !End second electron energy sample do loop
4000 continue
!* Loop through all of the angle probabilities to find the electron ejection
!* angle.
do eAng=1,neAngles !Loop through all of the electron ejection angles
  if(ranVecC(1).ge.aProbFunc(Proc,ChS,iBin,eAng))then !Random number vs. prob
    electron_angletmp1=eAngle(eAng) !Get the corresponding electron angle
    goto 5000 !Get out of do loop
  end if !End Monte Carlo if statement
end do !End first electron angle sample do loop
5000 continue !Continue on to find the next angle sample

do eAng=1,neAngles !Loop through all of the electron ejection angles
  if(f.lt.0.5)then !Need to know which ion energy bin we're interpolating with
    if(ranVecC(2).ge.aProbFunc(Proc,ChS,iBin+1,eAng))then
      electron_angletmp2=eAngle(eAng)
      electron_angle=(1-f)*electron_angletmp1+f*electron_angletmp2
      goto 6000 !Get out of do loop
    end if
  elseif(f.ge.0.5)then
    if(ranVecC(2).ge.aProbFunc(Proc,ChS,iBin-1,eAng))then
      electron_angletmp2=eAngle(eAng)
      electron_angle=(1-f)*electron_angletmp2+f*electron_angletmp1
      goto 6000 !Get out of do loop
    end if
  end if !End of Monte Carlo if statement
end do !End second electron angle sample do loop
6000 continue
! write(206,*) Proc,ChS,IonEnergy(iBin),electron_energy,electron_angle
do eEng=1,nE2strBins !Loop through the 2-stream energy bins
  if(electron_energy.le.E2str(eEng))then !Pick out which bin we need
    eBin=eEng !Assumes E2str(eEng) is the upper bound of the energy bin
    goto 7000
  end if
end do
7000 continue
end subroutine
