program SulfurIonPrecip
!*******************************************************************************
!* Created by Stephen J. Houston 10.19.18
!*******************************************************************************
!* This program simulates the path of an energetic sulfur ion as it penetrates
!* into the Jovian atmosphere. Different initial energies are considered.
!* I use 1-2,000 keV/u as initial energies for this exploratory work.
!* A random pitch angle is considered for the precipitating ion following a
!* cosine distribution.
!* A Monte Carlo simulation is used to determine the type of collision
!* and where the collision occurs. To calculate the collision path dN
!* I use 1-Prob = exp{-sigtot*dN}. After each collision I track the secondary
!* electrons produced, as well as the new ion energy. The results are binned by
!* column density which correspond to a particular altitude.
!******************************
!* External files required:
!* The secondary electron distributions are calculated by other codes based on
!* cross sections calculated by Dave Schultz.
!* The singly-differential cross sections for the ion+H2 processes are read in
!* at the beginning of the model. These are as function of energy.
!* The distribution functions for ejected electrons are read as eprobfunc and
!* aprobfunc. These are as function of energy and angle, respectively.
!* Ejected electron energy is calculated in this model by using the cross
!* sections and energy loss model presented in Schultz et al., 2018.
!* See input files for the files that are needed.
!******************************
!* Goals:
!* 1. Read in the distribution functions for all available collision types,
!* charge states and energies and store this info in matrices
!* 2. Set a matrix with the angular distribution to determine wether the
!* electron will be scattered forward or backward in a collision. The incident
!* ion pitch angle is added to the ejected electron angle.
!* 3. Read in all the total xs calculated by Dave.
!* 4. Create altitude bins and find the corresponding column density.
!* 5. Set the initial conditions for the ion -> charge state and initial energy,
!* incident angle (normally kept at 0) and initial pitch angle
!* 6. Follow ion as it penetrates the atmosphere and has collisions determined
!* by the MC.
!* 7. Track the charge state of the ion, number and energy of electrons
!* produced, and ion energy at each altitude bin in the atmosphere, until the
!* ion runs out of energy (E<1 keV/u)
!*******************************************************************************

use,intrinsic :: ISO_FORTRAN_ENV !Used for int64 integers
!use formatting !Formatting module to avoid cluttering the end of the program
implicit none

!**************************** Variable Declaration *****************************
integer i,j,k,l,m,n,run,ion !Do-loop variables
!* Parameter variables:
integer atmosLen !"Length" of the atmosphere (3000 - -88 km) with 2 km steps
integer nProjProc !Number of projectile processes
integer nTargProc !Number of target processes
integer nEnergies !Number of inital energies
integer nInterpEnergies !Number of interpolated energies
integer nChS !Number of sulfur charge states from 0-16
integer SS,DS,SPEX,DPEX !Projectile processes
integer SI,DI,TI,DA,SC,DC,TEX !Target processes

parameter(nProjProc=4) !Number of projectile processes
parameter(nTargProc=7) !Number of target processes
parameter(nEnergies=9) !Number of inital energies
parameter(nInterpEnergies=2000) !Number of interpolated energies
parameter(nChS=17) !Number of charge states from 0-16
parameter(SI=1,DI=2,TI=3,DA=4,SC=5,DC=6,TEX=7) !Target process numbers
parameter(SS=1,DS=2,SPEX=3,DPEX=4) !Projectile process numbers
parameter(atmosLen=1544) !Length of the atmosphere

real*8 dum

!* Computational time variables:
integer t1,t2,clock_maxTotal,clock_rateTotal !Used to calculate comp. time
integer t3,t4,clock_max,clock_rate !Used to calculate comp. time
real*8 hrs,min,sec

!* Sulfur ion variables
integer energy,trial,initChS,ChS,oldChS,dpt,excite,elect,disso,PID(2),nIons
integer numSim
real*8 incB,kappa,dE,dNTot,dZTot,E,pangle
real*8,dimension(nEnergies) :: IonEnergy
real*8,dimension(nChS,nInterpEnergies) :: SIMxs_Total
real*8,dimension(1+nProjProc,nChS,nInterpEnergies) :: SIMxs_Totaltmp
real*8,dimension(nTargProc,1+nProjProc,nChS,nInterpEnergies) :: SIMxs
!* SIMxs has an additonal projectile process which is no projectile process
!* i.e. SI, SI+SS, SI+DS, SI+SPEX, SI+DPEX (respectively)

!* Random Number Generator
integer k1,k2,lux,in
parameter(k1=0,k2=0,lux=3) !lux set to 3 for optimal randomness and timeliness
real ranVecA(10002)
real,allocatable :: angle(:)

!* Total column density, array of altitude in km, alt bin size, scale height
real*8,dimension(atmosLen) :: totalCD,altitude,altDelta,totalDens,H

!****************************** Data Declaration *******************************
!* Initial ion enegy input:
data IonEnergy/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0/
! data IonEnergy/10.625,15.017,20.225,29.783,46.653,59.770,77.522,120.647,&
!                218.125,456.250/ !Juno energy bins from JEDI
! data IonEnergy/10.625,11.619,12.656,13.786,15.017,16.177,17.427,18.774,20.225,&
!                22.280,24.543,27.036,29.783,33.319,37.276,41.702,46.653,49.634,&
!                52.806,56.180,59.770,63.785,68.070,72.642,77.522,86.586,96.710,&
!                108.018,120.647,139.90,162.223,188.108,218.125,262.319,315.467,&
!                379.384,456.250/ !JEDI energy bins interpolated
!********************************** Run Time ***********************************
!Calculate the total computational run time of the model:
call system_clock (t1,clock_rateTotal,clock_maxTotal)
!**************************** Initialize Variables *****************************
altitude=0.0;totalCD=0.0;totalDens=0.0;H=0.0;altDelta=0.0;SIMxs=0.0
SIMxs_Total=0.0;SIMxs_Totaltmp=0.0
!**************************** Create the Atmosphere ****************************
open(unit=200,file='./Atmosphere/Input/JunoColumnDensity_2km.dat',status='old')
open(unit=201,file='./Atmosphere/Input/JunoAtmosphere_2km.dat',status='old')
read(200,*);read(201,*) !Skip header lines
do i=1,atmosLen
  read(200,*)altitude(i),dum,dum,dum,dum,totalCD(i)!Read in the atmosphere
  read(201,*)dum,dum,dum,dum,dum,totalDens(atmosLen-i+1),dum,dum,H(atmosLen-i+1)
  altDelta(i)=2.0
end do
close(200) !Close column density file
close(201) !Close atmosphere file
!*************************** Get SIM cross-sections ****************************
open(unit=203,file='./SIMXSInterp/SIMXSInterpAll.dat',status='old')
read(203,20300) SIMxs !SIM cross-section data
close(203) !Close SIM cross-section data
20300 format (20(ES9.3E2,1x)) !Formatting for the file
SIMxs_Totaltmp=sum(SIMxs,dim=1) !Intermediate summing step
SIMxs_Total=sum(SIMxs_Totaltmp,dim=1) !Sum of cross-sections
!*********************** Ejected Electron Probabilities ************************
!* Created by ReadElectDist.f08 and ProbDist.f08
! open(unit=202,file='./NewElectronDist/ProbDistFunc/eprobfunc.dat',status='old')
! open(unit=203,file='./NewElectronDist/ProbDistFunc/aprobfunc.dat',status='old')
! read(202,*) eProbFunc
! read(203,*) aProbFunc
! close(202)
! close(203)

!*******************************************************************************
!******************************** MAIN PROGRAM *********************************
!*******************************************************************************
!* The following run number corresponds to the energy in keV/u
!* Juno:
!* 1=10, 2=15, 3=20, 4=30, 5=45, 6=60, 7=75, 8=120, 9=220, 10=450, 11=500,
!* 12=750, 13=1000, 14=1250, 15=1500, 16=1750, 17=2000, 18=2500, 19=3000,
!* 20=4000, 21=5000, 22=10000, 23=25000
!* Regular:
!* 1=1, 2=10, 3=50, 4=75, 5=100, 6=200, 7=500, 8=1000, 9=2000, 10=5000,
!* 11=10000, 12=25000
!*******************************************************************************
nIons=10!0 !Number of ions that are precipitating
trial=3 !The seed for the RNG
do run=4,4!nEnergies !Loop through different initial ion energies
  call system_clock(t3,clock_rate,clock_max) !Comp. time of each run
  energy=int(IonEnergy(run))
  write(*,*) "Number of ions:         ", nIons
  write(*,*) "Initial energy:         ", energy, 'keV'
  write(*,*) "Trial number (RNG Seed):", trial
  write(*,*) "***************************************************************&
              ****************************"
!*************************** Random Number Generator ***************************
  !k1=0,k2=0 Should be set to zero unless restarting at a break (See ranlux.f08)
  in=trial !RNG seed
  call rluxgo(lux,in,k1,k2) !Seed the RNG
  allocate(angle(nIons)) !Want the same number of angles as ions
  call ranlux(angle,nIons) !Calculate all the angles to be used
!********************* Reset Counters For New Ion Energies *********************
  ! totHp =0;totalElect=0;tElectFwd =0;tElectBwd =0;SPvsEng=0.0    ;nSPions=0
  ! totH2p=0;eCounts   =0;electFwdA =0;electBwdA =0;SigTotvsEng=0.0;maxDpt=0
  ! H2Ex  =0;oxygen    =0;electFwdAE=0;electBwdAE=0;dEvsEng=0.0
!************************ Ion Precipitation Begins Here ************************
  ! write(*,*) 'Starting Ion Precipitiaton: ', energy,'keV/u' !Double check energy
  do ion=1,nIons !Each ion starts here
    !*****************************
    !Initial Conditions:
    pangle=0.0         !Reset the pitch angle for every run
    incB  =0.0         !Incident B-field
    kappa =0.0         !Used to account for pitch angle
    numSim=energy*1000 !Number of simulations for a single ion. Must be great !~
                       !enough to allow the ion to lose all energy
    E=IonEnergy(run)   !Start with initial ion energy
    dE=0.0             !Energy loss
    initChS=2          !1 is an initial charge state of 0, 2 is +1
    ChS=initChS        !Set the charge state variable that will be changed
    oldChS=initChS     !Need another charge state variable for energyLoss.f08
    dNTot=0.0          !Reset the column density to the top of the atm.
    dZTot=3000.0       !Start from the top of the atmosphere
    dpt=4              !Depth of penetration for bins. (integer value)
    !Beginning scale height (H) at 4 or 5 seems to be more accurate than 1-3
    l=0                !Used as index for dN calculation (ranVecA(l))
    excite=0           !CollisionSim outputs
    PID=0              !Process identification numbers
    !*****************************
    pangle=(2.0*atan(1.0))-acos(angle(ion)) !Pitch angle calculation has a
    !cosine dist. Straight down is pitch angle of 0, random number must be 0
    ! write(*,*) 'Ion Number: ', ion, 'Pitch angle: ', pangle*90/acos(0.0)
    kappa=1.0/(cos(pangle)*cos(incB)) !Used to convert from ds to dz
    call ranlux(ranVecA,10002) !Get a random vector for collisions
    do i=1,1!numSim !This loop repeats after each collision until E < 1 keV/u
      ! !*****************************
      ! !Reset Variables:
      ! eEnergy=0.0;eEnergyTmp=0.0 !Ejected electron energy
      ! eAngle =0.0;ds        =1   !Double stripping electron angle variable
      ! addElect=0 ;processE  =0   !Ejected electron integers
      ! eAngleSS=0.0;eEnergySS=0.0
      ! eAngleDS=0.0;eEnergyDS=0.0
      ! !*****************************
      call CollisionSim(int(E),SIMxs,SIMxs_Total,ChS,excite,elect,disso,PID)
      write(*,*) E,PID,ChS,excite,elect,disso
    end do !End of i=1,numSim loop (E < 1 keV/u)
  end do !End of ion=1,nIons loop









  call system_clock(t4,clock_rate,clock_max) !Elapsed time for a single energy
  hrs=int(real(t4-t3)/clock_rate/3600.0)
  min=int(((real(t4-t3)/clock_rate)-hrs*3600)/60)
  sec=mod(real(t4-t3)/clock_rate,60.0)
  write(*,*) 'Individual run elapsed real time = ',hrs,':',min,':',sec
  deallocate(angle) !Angle variable is reallocated for each energy
end do !run=1,nEnergies


call system_clock (t2,clock_rateTotal,clock_maxTotal) !Total elapsed time
hrs=int(real(t2-t1)/clock_rateTotal/3600.0)
min=int(((real(t2-t1)/clock_rateTotal)-hrs*3600)/60)
sec=mod(real(t2-t1)/clock_rateTotal,60.0)
write(*,*) '***************************************************************&
            ****************************'
write(*,*) 'Total elapsed real time =          ',hrs,':',min,':',sec


end program
