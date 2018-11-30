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
use formatting !Formatting module to avoid cluttering the end of the program
implicit none

!**************************** Variable Declaration *****************************
!* Do-loop variables:
integer i,j,k,l,run,ion

!* Computational time variables:
integer t1,t2,clock_maxTotal,clock_rateTotal !Used to calculate comp. time
integer t3,t4,clock_max,clock_rate !Used to calculate comp. time
integer hrs,min
real sec

!* Atmosphere variables:
integer atmosLen !"Length" of the atmosphere (3000 - -88 km) with 2 km steps]
real*8 dN,dNTot,dZ,dZTot !Change in atmosphere
parameter(atmosLen=1544) !Length of the atmosphere
real*8,dimension(atmosLen) :: totalCD,altitude,altDelta,totalDens,H
!* Total column density, array of altitude in km, alt bin size, scale height
real*8 dum

!* Sulfur ion variables:
!PROJECTILE PROCESSES (pProc)
integer pProc,nProjProc !Number of projectile processes
integer noPP,SS,DS,SPEX,DPEX !Projectile processes + no projectile process(noPP)
parameter(nProjProc=4)
parameter(noPP=1,SS=2,DS=3,SPEX=4,DPEX=5)
!TARGET PROCESSES (tProc)
integer tProc,nTargProc !Number of target processes
integer SI,DI,TI,DA,SC,DC,TEX !Target processes
parameter(nTargProc=7)
parameter(SI=1,DI=2,TI=3,DA=4,SC=5,DC=6,TEX=7)
!CHARGE STATES (ChS)
integer ChS,ChS_init,ChS_old,nChS !Number of sulfur charge states from 0-16
parameter(nChS=17)
!ENERGY
integer Eng,energy,nEnergies !Number of inital ion energies
integer nInterpEnergies !Number of interpolated ion energies
real*8 E,dE
parameter(nEnergies=9,nInterpEnergies=2000)
real*8,dimension(nEnergies) :: IonEnergy !Initial ion energies
!CHARGE STATE DISTRIBUTION
integer nSulEngBins !Number of energy bins for charge state fractions
real*8 SulEngBinSize !Size of energy bins for charge state fractions
parameter(nSulEngBins=2000,SulEngBinSize=1.0)
integer(kind=int64) :: SulVsEng(nChS,nSulEngBins) !Charge state fractions
real*8 engBins(nSulEngBins),SulEngBins(nSulEngBins) !Sulfur energy bins
!STOPPING POWER (SP)
integer nSPBins !Number of stopping power bins vs. energy
real*8 SPBinSize !Size of bins
parameter(nSPBins=2000,SPBinSize=1.0)
real*8 dEsp,SIMxsTotSP,dEold
integer(kind=int64),dimension(nSPBins) :: nSPions
real*8,dimension(nSPBins) :: SPBins,SPvsEng,SIMxsTotvsEng,dEvsEng,dNvsEng
!MASS
real*8 mass !Atomic mass of sulfur (32.065)
parameter(mass=32.065)
!OTHER VARIABLES
integer trial,excite,elect,disso,PID(2),nIons
integer numSim,dpt,maxDpt
character(len=100) arg !For reading in trial number on the cluster

integer(kind=int64),dimension(nTargProc,1+nProjProc) :: collisions !Counter
integer(kind=int64),dimension(nTargProc,1+nProjProc,nChS,atmosLen) :: Sulfur
real*8 incB,kappa,pangle
real*8,dimension(nChS,nInterpEnergies) :: SIMxs_Total !SigTot for dN calculation
real*8,dimension(1+nProjProc,nChS,nInterpEnergies) :: SIMxs_Totaltmp
real*8,dimension(nTargProc,1+nProjProc,nChS,nInterpEnergies) :: SIMxs !All SIMxs
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

!* Output production variables:
integer(kind=int64),dimension(atmosLen) :: Hp,H2p !H+ and H2+ counts
integer(kind=int64) :: H2Ex(atmosLen) !Excited H2 counter
real*8 norm !Normalization variable to per ion per cm
real dissRan !Random number to determine dissociation probability

!* Random Number Generator:
integer k1,k2,lux,in
parameter(k1=0,k2=0,lux=3) !lux set to 3 for optimal randomness and timeliness
real ranVecA(1002)
real,allocatable :: angle(:)

!* Output Variables:
integer nOutputFiles !Number of output files
parameter(nOutputFiles=10)
character(len=100) filename,files(nOutputFiles) !Output file names
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
data files/'ChargeStateDistribution','H+_Prod','H2+_Prod','H2*_Prod',&
           'Collisions','Photons_CX','Photons_DE','Stopping_Power',&
           '2Str_Elect_Fwd','2Str_Elect_Bwd'/
!********************************** Run Time ***********************************
!Calculate the total computational run time of the model:
call system_clock (t1,clock_rateTotal,clock_maxTotal)
!**************************** Initialize Variables *****************************
altitude=0.0;totalCD=0.0;totalDens=0.0;H=0.0;altDelta=0.0;SIMxs=0.0
SIMxs_Total=0.0;SIMxs_Totaltmp=0.0;SulEngBins=0.0;SPBins=0.0
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
open(unit=203,file='./SIMXSInterp/SIMXSInterpAll.txt',status='old')
do tProc=1,nTargProc
  do pProc=1,nProjProc+1
    do i=1,2
      read(203,*)
    end do
    do Eng=1,nInterpEnergies
      read(203,20300) dum,(SIMxs(tProc,pProc,ChS,Eng),ChS=1,nChS)
    end do
  end do
end do
20300 format(F7.2,17(1x,ES9.3E2))
! read(203,20300) SIMxs !SIM cross-section data
! close(203) !Close SIM cross-section data
! 20300 format (20(ES9.3E2,1x)) !Formatting for the file
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
!Stopping power bins:
SPBins(1)=1.0
do i=2,nSPBins
  SPBins(i)=SPBins(i-1)+SPBinSize
end do
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
nIons=100 !Number of ions that are precipitating
call get_command_argument(1,arg)
read(arg,'(I100)') trial !The seed for the RNG
do run=9,9!,nEnergies !Loop through different initial ion energies
  call system_clock(t3,clock_rate,clock_max) !Comp. time of each run
  energy=int(IonEnergy(run))
  write(filename,"('../scratch/Jup_Sul/Output/',I0,'/Seeds.dat')") energy
  open(unit=205,file=filename,status='unknown',access='append',action='write')
  write(205,*) trial
  close(205)
  write(filename,"('../scratch/Jup_Sul/Output/',I0,'/Overview',I0,'.dat')")&
    energy,trial
  open(unit=206,file=filename,status='unknown')
  write(206,*) "Number of ions:         ",nIons
  write(206,*) "Initial energy:         ",energy,'keV/u'
  write(206,*) "Trial number (RNG Seed):",trial
  write(206,F3) !'**'
!*************************** Random Number Generator ***************************
  !k1=0,k2=0 Should be set to zero unless restarting at a break (See ranlux.f08)
  in=trial !RNG seed
  call rluxgo(lux,in,k1,k2) !Seed the RNG
  allocate(angle(nIons)) !Want the same number of angles as ions
  call ranlux(angle,nIons) !Calculate all the angles to be used
!********************* Reset Counters For New Ion Energies *********************
  Hp =0;totalElect=0
  H2p=0;Sulfur    =0;electFwd  =0;electBwd  =0;maxDpt   =0
  H2Ex  =0;collisions=0;SulVsEng  =0
  SPvsEng=0.0;SIMxsTotvsEng=0.0;dEvsEng=0.0;dNvsEng=0.0;nSPions=0
  ! pHp =0.0;totalH2p=0.0;collisions=0;npHp      =0;npH2p        =0
  ! pH2p=0.0;totO      =0;dNvsEng =0.0;oxygenCX=0.0;prode2stF  =0.0;prode2stB=0.0
  !SPvsEng    =0.0;nSPions  =0;totalHp =0.0;dEvsEng    =0.0;SIMxsTotvsEng=0.0
!************************ Ion Precipitation Begins Here ************************
  write(206,*) 'Starting Ion Precipitiaton: ', energy,'keV/u' !Double check energy
  flush(206)
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
    write(206,*) 'Ion Number: ',ion,' Pitch angle: ',pangle*90/acos(0.0)
    flush(206)
    kappa=1.0/(cos(pangle)*cos(incB)) !Used to convert from ds to dz
    call ranlux(ranVecA,1002) !Get a random vector for collisions
    do i=1,numSim !This loop repeats after each collision until E < 1 keV/u
      !*****************************
      !Reset Variables:
      dN=0.0;dZ=0.0;dE=0.0 !Change in column density, altitude, energy
      dEsp=0.0;SIMxsTotSP=0.0;dEold=0.0 !Stopping power variables
      !*****************************
      call CollisionSim(nint(E),SIMxs,SIMxs_Total,ChS,excite,elect,disso,PID)
      collisions(PID(1),PID(2))=collisions(PID(1),PID(2))+1 !Count collisions
1000 continue
      l=l+1
      if(l.ge.1000)then
        !Filling ranVecA with a huge amount of numbers is a big time waster
        call ranlux(ranVecA,1002) !Only get ranVecA as needed
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
          write(206,*)"JupSulPrecip.f08: WARNING: Ion exited the bottom of the &
                    &atmosphere, proceeding to next ion."
          goto 5000 !Continue on to the next ion
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
          write(206,*) "JupSulPrecip.f08: WARNING: Elect ejection angle &
                     &greater than 270 degrees."
        end if
        !Only want to add the electron energies for the some processes since SS
        !and DS have to be transformed into a different reference frame
        if(eProc.le.4)eEnergy=eEnergy+eEnergyTmp !Units of eV
      end do !End of electron ejection do loop (j=1,elect)
!************************* Counting Photon Production **************************
!* Note:
!*  A photon count at a specific altitude and charge state means that there was
!*  a photon producing collision at that specific altitude and the resultant ion
!*  was at the recorded charge state. That means, a collision that goes from
!*  S^8+ to S^7+ will be recorded as S^7+; therefore, the S^16+ bin will never
!*  produce a photon. The last bin should ALWAYS be 0. Each processes can only
!*  create one photon, if it's a photon producing collision.
!*  Direct excitation (photonsDE) producing collisions:
!*    TEX+SPEX(27) ,SI+SPEX(29), DI+SPEX(32)
!*  Charge exchange (photonsCX) producing collisions:
!*    SC+SS(19), TI(25), SC(30)
!*  This is all done in the writing of the output files. Photon productions are
!*  files 118 and 119 using the sulfur variable.
!*******************************************************************************
!********************** Counting Sulfur & H/H2 Production **********************
      Sulfur(PID(1),PID(2),ChS,dpt)=Sulfur(PID(1),PID(2),ChS,dpt)+1 !Sulfur prod
      !Sulfur variable is what shows photon production, depending on processes
      if(disso.eq.2)then
        Hp(dpt)=Hp(dpt)+2 !Number of H^+ produced
      elseif(disso.eq.1)then
        call ranlux(dissRan,1) !Random number to determine dissociation
        if(dissRan.le.0.1)then !10% chance of dissociation (H + H^+)
          Hp(dpt)=Hp(dpt)+1
        else !90% chance of no dissociation
          H2p(dpt)=H2p(dpt)+1
        end if
      elseif(disso.eq.0)then !TEX never dissociates, result is H2*
        H2Ex(dpt)=H2Ex(dpt)+1
      end if
!************************** Energy Loss Calculations ***************************
      call energyloss(E,ChS_old,eEnergy,PID,eEnergySS,eAngleSS,&
                     &eEnergyDS,eAngleDS,dE)
      dEsp=(dE)/dN !stopping power (calc before dE is recalculated)
      dEold=dE
      dE=(1/mass)*(1.0e-3)*dE*kappa !Total dE function
      if(dN.lt.0.0)then !Change in column density should never be less than 0
        write(206,10001) E,dEsp,dE,dN,dEold,PID(1),PID(2),ChS_old
      end if
!********************** Sulfur Charge State Distribution ***********************
      do j=1,nSulEngBins
        if(E.le.SulEngBins(j))then
          SulVsEng(ChS,j)=SulVsEng(ChS,j)+1
          goto 3000
        end if
      end do
3000 continue
      do j=1,nSPBins
        if(E.le.SPBins(j))then
          SPvsEng(j)=SPvsEng(j)+dEsp !Stopping power vs energy
          SIMxsTotvsEng(j)=SIMxsTotvsEng(j)+SIMxsTotSP !Total cross-section vs energy
          dEvsEng(j)=dEvsEng(j)+dEold !Change in energy vs energy
          dNvsEng(j)=dNvsEng(j)+dN !Change in column density vs energy
          ! ProcessdE(j,processC(process),tempQold)=&
          !   ProcessdE(j,processC(process),tempQold)+dEold
          nSPions(j)=nSPions(j)+1 !Number of ions in each energy bin
          ! pnSPions(j,processC(process),tempQold)=&
          !   pnSPions(j,processC(process),tempQold)+1
          goto 4000
        end if
      end do
4000 continue
      E=E-dE
      ChS_old=ChS !Assign newly acquired charge state to old variable
      if(E.lt.1.0) goto 5000 !Stop once the energy is less than 1 keV/u
      if(i.eq.numSim)then
        write(206,*) 'JupSulPrecip.f08: ERROR: numSim not large enough.'
        write(206,*) 'JupSulPrecip.f08: Ion energy was: ',E
        goto 5000
      end if
    end do !End of i=1,numSim loop (E < 1 keV/u)
5000 continue
  end do !End of ion=1,nIons loop
!******************************** Output Header ********************************
  write(206,*) '--------------------------NEW RUN---------------------------'
  write(206,*) 'Number of ions: ', nIons
  write(206,*) 'Initial energy: ', energy, 'keV'
  write(206,*) 'Trial number:   ', trial
  write(206,F3) !'**'
  !******* Check various electron counters
  write(206,*) 'Sum of total electrons foward:         ',sum(electFwd)
  write(206,*) 'Sum of total electrons backward:       ',sum(electBwd)
  write(206,*) 'Sum of total electrons foward+backward:',sum(electFwd+electBwd)
  write(206,*) 'Sum of total electrons:                ',totalElect
  write(206,F1)'Max Depth:',altitude(maxDpt)
!********** Open output data files for each set of initial energies ************
  do i=1,nOutputFiles
    write(filename,"('../scratch/Jup_Sul/Output/',I0,'/',A,'-',I0,'.dat')")&
      energy,trim(files(i)),trial
    open(unit=100+i,file=trim(filename))
  end do
!***************************** Write out to files ******************************
  norm=nIons*2e5 !Normalization condition to per ion per cm
!*** Charge state distribution
  write(101,H01) !Sulfur charge state distribution header
  do i=2,nSulEngBins !Sulfur charge state distribution
    write(101,F01) SulEngBins(i)-(SulEngBinSize/2.0),&
      (real(SulVsEng(j,i))/real(sum(SulVsEng(:,i))),j=1,nChS)
  end do
!*** H production
  write(102,H02) !H^+ Header
  write(103,H03) !H_2^+ Header
  write(104,H04) !H_2^* Header
  do i=1,atmosLen !Loop through the atmosphere
    write(102,F02) altitude(i),Hp(i)/norm !H^+ production
    write(103,F02) altitude(i),H2p(i)/norm !H_2^+ production
    write(104,F02) altitude(i),H2Ex(i)/norm !H_2^* production
  end do
!*** Collision counters
  write(105,F03) (ProjColl(i),i=1,nProjProc) !Collisions header
  do i=1,nTargProc !Total number of each type of collision
    write(105,F04) TargColl(i),(collisions(i,j),j=1,5),sum(collisions(i,:))
  end do
  write(105,F4) !'--'
  write(105,F04) 'Sum ',(sum(collisions(:,i)),i=1,nProjProc+1),sum(collisions)
  write(105,*) ''
  write(105,F03) (ProjColl(i),i=1,nProjProc) !Collisions precentage header
  do i=1,nTargProc !Total number of each type of collision
    write(105,F05) TargColl(i),&
      (real(collisions(i,j))/real(sum(collisions))*100,j=1,5),&
      real(sum(collisions(i,:)))/real(sum(collisions))*100
  end do
  write(105,F4) !'--'
  write(105,F05) 'Sum ',&
    (real(sum(collisions(:,i)))/real(sum(collisions))*100,i=1,nProjProc+1),&
    real(sum(collisions))/real(sum(collisions))*100
!*** Photon production
  write(106,N01) !CX note
  write(107,N02) !DE note
  do i=1,2 !Loop through CX and DE headers
    write(105+i,*) !Blank space
    write(105+i,H09) !Initial input header
    write(105+i,*) !Blank space
    write(105+i,H05) !Altitude integrated photon production header
    write(105+i,H06) !Charge state header
  end do
  !* Altitude integrated photon production
  write(106,F06) altDelta(1),& !CX - TI, SC, SC+SPEX
    (real(sum(sulfur(TI,noPP,ChS,:))+sum(sulfur(SC,noPP,ChS,:))+&
    sum(sulfur(SC,SPEX,ChS,:)))/norm,ChS=1,nChS)
  write(107,F06) altDelta(1),& !DE - SI+SPEX, DI+SPEX, TEX+SPEX
    (real(sum(sulfur(SI,SPEX,ChS,:))+sum(sulfur(DI,SPEX,ChS,:))+&
    sum(sulfur(TEX,SPEX,ChS,:)))/norm,ChS=1,nChS)
  do i=1,2 !Loop through CX and DE headers
    write(105+i,*) !Blank space
    write(105+i,H07) !Photon production vs. altitude header
    write(105+i,H08) !Charge state header
  end do
  do i=1,atmosLen
    write(106,F06) altitude(i),& !CX - TI, SC, SC+SPEX
     (real(sulfur(TI,noPP,ChS,i)+sulfur(SC,noPP,ChS,i)+sulfur(SC,SPEX,ChS,i))/&
     norm,ChS=1,nChS)
    write(107,F06) altitude(i),& !DE - SI+SPEX, DI+SPEX, TEX+SPEX
     (real(sulfur(SI,SPEX,ChS,i)+sulfur(DI,SPEX,ChS,i)+sulfur(TEX,SPEX,ChS,i))/&
     norm,ChS=1,nChS)
  end do
!*** Stopping power
  write(108,H10) !Stopping power header
  do i=2,nSPBins !Loop through ever stopping power bin
    write(108,F07) SPBins(i)-(SPBinSize/2.0),&
      SPvsEng(i)/real(nSPions(i)**2),&
      SIMxsTotvsEng(i)/real(nSPions(i)),&
      dEvsEng(i)/real(nSPions(i)),&
      dNvsEng(i)/real(nSPions(i)),&
      (SIMxsTotvsEng(i)*dEvsEng(i))/real(nSPions(i)**2),&
      nSPions(i)
  end do
!***************************** Secondary Electrons *****************************
do j=1,nE2strBins !2-stream electrons, forward and backward
  write(109,F2str) (real(electFwd(i,j))/norm,i=atmosLen,1,-1)
  write(110,F2str) (real(electBwd(i,j))/norm,i=atmosLen,1,-1)
end do
!******************************* Close all files *******************************
  do i=1,nOutputFiles
    close(100+i)
  end do

  call system_clock(t4,clock_rate,clock_max) !Elapsed time for a single energy
  hrs=int(real(t4-t3)/clock_rate/3600.0)
  min=int(((real(t4-t3)/clock_rate)-hrs*3600)/60)
  sec=mod(real(t4-t3)/clock_rate,60.0)
  ! write(206,F2) 'Individual run elapsed real time = ',hrs,':',min,':',sec
  deallocate(angle) !Angle variable is reallocated for each energy
end do !run=1,nEnergies


call system_clock (t2,clock_rateTotal,clock_maxTotal) !Total elapsed time
hrs=int(real(t2-t1)/clock_rateTotal/3600.0)
min=int(((real(t2-t1)/clock_rateTotal)-hrs*3600)/60)
sec=mod(real(t2-t1)/clock_rateTotal,60.0)
write(206,F3) !'**'
write(206,F2) 'Total elapsed real time = ',hrs,min,sec
close(206)
write(filename,"('../scratch/Jup_Sul/Output/',I0,'/Elapsed_Times.dat')") energy
open(unit=207,file=filename,access='append',action='write')
write(207,*) trial,hrs,':',min,':',sec
close(207)
10001 format(F8.2,3(2x,ES10.3E2),2x,F9.2,4(2x,I2))
end program
