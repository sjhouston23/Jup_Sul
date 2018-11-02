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
!* Do-loop variables:
integer i,j,k,l,m,n,run,ion

!* Computational time variables:
integer t1,t2,clock_maxTotal,clock_rateTotal !Used to calculate comp. time
integer t3,t4,clock_max,clock_rate !Used to calculate comp. time
real*8 hrs,min,sec

!* Atmosphere variables:
integer atmosLen !"Length" of the atmosphere (3000 - -88 km) with 2 km steps]
parameter(atmosLen=1544) !Length of the atmosphere
real*8,dimension(atmosLen) :: totalCD,altitude,altDelta,totalDens,H
!* Total column density, array of altitude in km, alt bin size, scale height
real*8 dum

!* Sulfur ion variables:
integer nProjProc !Number of projectile processes
parameter(nProjProc=4)

integer nTargProc !Number of target processes
parameter(nTargProc=7)

integer nChS !Number of sulfur charge states from 0-16
parameter(nChS=17)

integer nEnergies !Number of inital ion energies
integer nInterpEnergies !Number of interpolated ion energies
parameter(nEnergies=9,nInterpEnergies=2000)

integer SI,DI,TI,DA,SC,DC,TEX !Target processes
parameter(SI=1,DI=2,TI=3,DA=4,SC=5,DC=6,TEX=7)

integer SS,DS,SPEX,DPEX !Projectile processes
parameter(SS=2,DS=3,SPEX=4,DPEX=5)

integer nSulEngBins !Number of energy bins for charge state fractions
real*8 SulEngBinSize !Size of energy bins for charge state fractions
parameter(nSulEngBins=2000,SulEngBinSize=1.0)

real*8 mass !Atomic mass of sulfur (32.065)
parameter(mass=32.065)

integer energy,trial,ChS_init,ChS,ChS_old,dpt,excite,elect,disso,PID(2),nIons
integer numSim,maxDpt

integer(kind=int64) :: SulVsEng(nChS,nSulEngBins) !Charge state fractions
real*8 engBins(nSulEngBins),SulEngBins(nSulEngBins) !Sulfur energy bins

real*8 incB,kappa,dE,dN,dNTot,dZ,dZTot,E,pangle,SIMxsTotSP
real*8,dimension(nEnergies) :: IonEnergy
real*8,dimension(nChS,nInterpEnergies) :: SIMxs_Total
real*8,dimension(nTargProc,1+nProjProc) :: collisions !Collision counter
real*8,dimension(1+nProjProc,nChS,nInterpEnergies) :: SIMxs_Totaltmp
real*8,dimension(nTargProc,1+nProjProc,nChS,nInterpEnergies) :: SIMxs
!* SIMxs has an additonal projectile process which is no projectile process
!* i.e. SI, SI+SS, SI+DS, SI+SPEX, SI+DPEX (respectively)

!* Ejected electron variables:
integer neProc !Number of processes that eject electrons
integer nE2strBins !Number of 2 stream bins
integer eSI,eDI,eTI,eDA,eSS,eDS !Electron ejection processes

parameter(neProc=6) !Number of processes that eject electrons
parameter(nE2strBins=260) !Number of 2 stream bins
parameter(eSI=1,eDI=2,eTI=3,eDA=4,eSS=5,eDS=6) !Electron ejection processes

integer eProc !Electron process number
integer nElect !Number of electron counter
integer eBin !Electron energy bin in 2-stream format
integer DSelect !Double stripping electron tracker for transform
real*8 eEnergyTmp,eAngleTmp,eEnergy,eAngle
real*8 eAngleSS,eEnergySS,eAngleDS(2),eEnergyDS(2) !SS/DS elec transforms
integer(kind=int64) :: totalElect !Total Electrons
integer(kind=int64),dimension(atmosLen,nE2strBins) :: electFwd,electBwd

!* Random Number Generator:
integer k1,k2,lux,in
parameter(k1=0,k2=0,lux=3) !lux set to 3 for optimal randomness and timeliness
real ranVecA(10002)
real,allocatable :: angle(:)


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
data engBins/nSulEngBins*SulEngBinSize/ !Used for sulfur binning
!********************************** Run Time ***********************************
!Calculate the total computational run time of the model:
call system_clock (t1,clock_rateTotal,clock_maxTotal)
!**************************** Initialize Variables *****************************
altitude=0.0;totalCD=0.0;totalDens=0.0;H=0.0;altDelta=0.0;SIMxs=0.0
SIMxs_Total=0.0;SIMxs_Totaltmp=0.0;SulEngBins=0.0
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
!**************************** Various Bin Creation *****************************
! !2-Stream energy bins:
! do i=1,nE2strBins
!   es=es+del(i)
!   E2str(i)=Es
! end do
!Sulfur bins for charge state fractions:
SulEngBins(1)=SulEngBinSize
do i=2,nSulEngBins
  SulEngBins(i)=SulEngBins(i-1)+engBins(i) !1-2000 keV/u
end do
! !Stopping power bins:
! es=1.0
! do i=1,nStopPowerEBins
!   es=es+delSP(i)
!   stopPowerEBins(i)=es
! end do
! !Bins for ejected electron energy:
! es=0.0
! do i=1,790
!   es=es+delAVGe(i)
!   elecEbins(i)=es
! end do
! !Bins for ejected electron angle:
! do i=1,180
!   elecAbins(i)=i
! end do
!*******************************************************************************
!******************************** MAIN PROGRAM *********************************
!*******************************************************************************
!* The following run number corresponds to the energy in keV/u
!* Juno:
!* 1=10, 2=15, 3=20, 4=30, 5=45, 6=60, 7=75, 8=120, 9=220, 10=450, 11=500,
!* 12=750, 13=1000, 14=1250, 15=1500, 16=1750, 17=2000, 18=2500, 19=3000,
!* 20=4000, 21=5000, 22=10000, 23=25000
!* Regular:
!* 1=1, 2=10, 3=50, 4=75, 5=100, 6=200, 7=500, 8=1000, 9=2000
!*******************************************************************************
nIons=150 !Number of ions that are precipitating
trial=3 !The seed for the RNG
do run=nEnergies,nEnergies !Loop through different initial ion energies
  call system_clock(t3,clock_rate,clock_max) !Comp. time of each run
  energy=int(IonEnergy(run))
  write(*,*) "Number of ions:         ", nIons
  write(*,*) "Initial energy:         ", energy, 'keV'
  write(*,*) "Trial number (RNG Seed):", trial
  write(*,*) "***************************************************************&
             &****************************"
!*************************** Random Number Generator ***************************
  !k1=0,k2=0 Should be set to zero unless restarting at a break (See ranlux.f08)
  in=trial !RNG seed
  call rluxgo(lux,in,k1,k2) !Seed the RNG
  allocate(angle(nIons)) !Want the same number of angles as ions
  call ranlux(angle,nIons) !Calculate all the angles to be used
!********************* Reset Counters For New Ion Energies *********************
! totHp =0;totalElect=0;tElectFwd =0;tElectBwd =0;SPvsEng    =0.0;nSPions  =0
! totH2p=0;eCounts   =0;electFwdA =0;electBwdA =0;SigTotvsEng=0.0;maxDpt   =0
! H2Ex  =0;oxygen    =0;electFwdAE=0;electBwdAE=0;dEvsEng    =0.0;totalHp  =0.0
! pHp =0.0;totalH2p=0.0;collisions=0;npHp      =0;npH2p        =0;OxyVsEng =0.0
! pH2p=0.0;totO      =0;dNvsEng =0.0;oxygenCX=0.0;prode2stF  =0.0;prode2stB=0.0
! NSIM  =0;SIM       =0
SulVsEng=0
!************************ Ion Precipitation Begins Here ************************
  write(*,*) 'Starting Ion Precipitiaton: ', energy,'keV/u' !Double check energy
  do ion=1,nIons !Each ion starts here
    !*****************************
    !Initial Conditions:
    pangle=0.0         !Reset the pitch angle for every run
    incB  =0.0         !Incident B-field
    kappa =0.0         !Used to account for pitch angle
    numSim=energy*1000 !Number of simulations for a single ion. Must be great !~
                       !enough to allow the ion to lose all energy
    E=IonEnergy(run)   !Start with initial ion energy
    ChS_init=nChS         !1 is an initial charge state of 0, 2 is +1
    ChS=ChS_init       !Set the charge state variable that will be changed
    ChS_old=ChS_init   !Need another charge state variable for energyLoss.f08
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
    write(*,*) 'Ion Number: ', ion, 'Pitch angle: ', pangle*90/acos(0.0)
    kappa=1.0/(cos(pangle)*cos(incB)) !Used to convert from ds to dz
    call ranlux(ranVecA,10002) !Get a random vector for collisions
    do i=1,numSim !This loop repeats after each collision until E < 1 keV/u
      !*****************************
      !Reset Variables:
      dN=0.0;dZ=0.0;dE=0.0 !Change in column density, altitude, energy
      !*****************************
      call CollisionSim(nint(E),SIMxs,SIMxs_Total,ChS,excite,elect,disso,PID)
      collisions(PID(1),PID(2))=collisions(PID(1),PID(2))+1 !Count collisions
1000 continue
      l=l+1
      if(l.ge.10000)then
        !Filling ranVecA with a huge amount of numbers is a big time waster
        call ranlux(ranVecA,10002) !Only get ranVecA as needed
        l=1 !Reset l back to 1 (Start at 1 because ranVecA(l) is called next)
      end if
      !Calculate how far ion moves before a collision (dN)
      dN=-log(1-ranVecA(l))/SIMxs_Total(ChS_old,nint(E))
      !Sometimes ranVecA is small enough to make DN 0
      if(dN.lt.1.0)goto 1000 !Get a new dN
      SIMxsTotSP=SIMxs_Total(ChS_old,nint(E)) !Used for stopping power calc.
      dNTot=dNTot+dN !Total change in column density
      do j=1,atmosLen !Loop through all of the atmosphere
        if(dNTot.le.totalCD(j+1))then !Move to proper CD of atmosphere
          !Calculate change in z based on the movement through the column dens.
          dZ=log((cos(pangle)*dN/(totalDens(dpt)*H(dpt)))+1)*H(dpt)
          dZTot=dZTot-dZ*1e-5 !Convert to km and keep subtracting from alt.
          do k=1,atmosLen !Loop through the atmosphere again
            if(dZTot.gt.altitude(k))then
              dNTot=totalCD(k) !Check the purpose - So ion doesn't get stuck in a bin??
              dpt=k !dpt is now the bin corresponding to depth
              if(dpt.gt.maxDpt) maxDpt=dpt !Used to see how deep we go
              goto 2000 !Get out of the do-loop that finds depth of penetration
            end if !Altitude if-statemet
          end do !Altitude do-loop
          !If we get here, then the ion has went through the entire atmosphere
          write(*,*)"JupSulPrecip.f08: WARNING: Ion exited the bottom of the &
                    &atmosphere, proceeding to next ion."
          goto 4000 !Continue on to the next ion
        end if !Column density if-statement
      end do !Column density do-loop
2000 continue
!*********************** Secondary Electron Calculations ***********************
      !*****************************
      !Reset Variables:
      eEnergy =0.0;eEnergyTmp=0.0 !Ejected electron energy (eV)
      eAngle  =0.0;eAngleTmp =0.0 !Ejected electron angle (°)
      eAngleSS=0.0;eEnergySS =0.0 !Single stripping transform variables (eV)
      eAngleDS=0.0;eEnergyDS =0.0 !Double stripping transform variables (°)
      nElect  =0  ;eProc     =0   !Ejected electron integers
      DSelect =1                  !Double stripping electron counter
      !*****************************
      if(PID(1).eq.SC.or.PID(1).eq.DC.or.PID(1).eq.TEX)nElect=12 !No electrons
      do j=1,elect !Loop through all of the ejected electrons
        if(PID(1).eq.SI.and.nElect.le.10)then !Single Ionization
          eProc=eSI !1
          nElect=nElect+10 !After one time, don't want to come back in here
        elseif(PID(1).eq.DI.and.nElect.le.10)then !Double Ionization
          eProc=eDI !2
          nElect=nElect+5 !After two times, don't want to come back in here
        elseif(PID(1).eq.TI.and.nElect.le.10)then !Transfer Ionization
          eProc=eTI !3
          nElect=nElect+10 !After one time, don't want to come back in here
        elseif(PID(1).eq.DA.and.nElect.le.10)then !Double Capture Autoionization
          eProc=eDA !4
          nElect=nElect+10 !After one time, don't want to come back in here
        elseif(PID(2).eq.SS.and.nElect.ge.11)then !Single Stripping
          eProc=eSS !5
        elseif(PID(2).eq.DS.and.nElect.ge.11)then !Double Stripping
          eProc=eDS !6
        end if
        call EjectedElectron(E,eProc,ChS_old,eEnergyTmp,eAngleTmp,eBin)
        nElect=nElect+1
        if(eProc.eq.5)then
          eEnergySS=eEnergyTmp !Units of eV
          eAngleSS=eAngleTmp !Need the ejection angle for energy transformation
        end if
        if(eProc.eq.6)then
          eEnergyDS(DSelect)=eEnergyTmp !Units of eV
          eAngleDS(DSelect)=eAngleTmp
          DSelect=DSelect+1 !Double stripping electron
        end if
        totalElect=totalElect+1 !Total number of electrons produced
        eAngle=eAngleTmp+(pangle*90/acos(0.0))
        !Must add the pitch angle to ejected elect angle. pangle = [0,acos(0.0)]
        if(eAngle.le.90.0)then !Counting electrons going forward (downward)
          electFwd(dpt,eBin)=electFwd(dpt,eBin)+1 !Elect fwd vs. alt. and eng.
        elseif(eAngle.le.270.0)then !Electrons going backward (0 is down)
          electBwd(dpt,eBin)=electBwd(dpt,eBin)+1 !Elect bwd vs. alt. and eng.
        else !If the electron is ejected so far backward it's going fwd again
          write(*,*) "JupSulPrecip.f08: WARNING: Elect ejection angle &
                     &greater than 270 degrees."
        end if
        !Only want to add the electron energies for the some processes since SS
        !and DS have to be transformed into a different reference frame
        if(eProc.le.4)eEnergy=eEnergy+eEnergyTmp !Units of eV
      end do !End of electron ejection do loop (j=1,elect)
!************************** Energy Loss Calculations ***************************
      call energyloss(E,ChS_old,eEnergy,PID,eEnergySS,eAngleSS,&
                     &eEnergyDS,eAngleDS,dE)
      ! dEsp=(dE)/dN !stopping power (calc before dE is recalculated)
      ! dEold=dE
      dE=(1/mass)*(1.0e-3)*dE*kappa !Total dE function
      ! if(dN.lt.0.0)then !Change in column density should never be less than 0
      !   write(206,10001) E,dEsp,dE,dN,dEold,process,PID(1),PID(2),ChS_old
      ! end if
!********************** Sulfur Charge State Distribution ***********************
      do j=1,nSulEngBins
        if(E.le.SulEngBins(j))then
          SulVsEng(ChS,j)=SulVsEng(ChS,j)+1
          goto 3000
        end if
      end do
3000 continue
      E=E-dE
      ChS_old=ChS !Assign newly acquired charge state to old variable
      if(E.lt.1.0) goto 4000 !Stop once the energy is less than 1 keV/u
      if(i.eq.numSim)then
        write(*,*) 'JupSulPrecip.f08: ERROR: numSim not large enough.'
        write(*,*) 'JupSulPrecip.f08: Ion energy was: ',E
        goto 4000
      end if
    end do !End of i=1,numSim loop (E < 1 keV/u)
4000 continue
  end do !End of ion=1,nIons loop
  open(unit=100,file='./Output/ChargeStateDistribution.dat')
  ! totO=sum(SulVsEng,dim=1)
  ! write(104,H04) !Sul vs energy header
  do i=1,nSulEngBins !Sulfur charge state distribution
    write(100,10000) SulEngBins(i)-(SulEngBinSize/2.0), &
                    &(real(SulVsEng(j,i))/real(sum(SulVsEng(:,i))),j=1,nChS)
  end do
  close(100)
  10000 format (1x,F8.2,5x,17(2x,ES8.2))

  call system_clock(t4,clock_rate,clock_max) !Elapsed time for a single energy
  hrs=int(real(t4-t3)/clock_rate/3600.0)
  min=int(((real(t4-t3)/clock_rate)-hrs*3600)/60)
  sec=mod(real(t4-t3)/clock_rate,60.0)
  ! write(*,*) 'Individual run elapsed real time = ',hrs,':',min,':',sec
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
