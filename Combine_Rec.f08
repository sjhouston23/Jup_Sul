program Combine_Output_Recursive
!*******************************************************************************
!* Created by Stephen J. Houston 11.29.18
!*******************************************************************************
!* This program combines multiple files of the same information into one.
!* Used after making a large run on multiple cores (CRC) to simulate as if a
!* single run had taken place. This loops through various energies.
!*******************************************************************************

use, intrinsic :: ISO_FORTRAN_ENV
use formatting
implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer energy,atmosLen,ChS,err,start,run

integer nProjProc !Number of projectile processes
integer nTargProc !Number of target processes
parameter(nProjProc=4,nTargProc=7)

parameter(nChS=17,atmosLen=1544)
parameter(nOutputFiles=10,MaxnTrials=1000,MaxnLines=100000)
parameter(nSulEngBins=2000,nSPBins=2000)
parameter(nEnergiesNorm=23,nEnergiesJuno=36)
parameter(nE2strBins=260) !Number of 2 stream bins

integer trial(MaxnTrials),nLines(nOutputFiles) !Number of trials/lines in a file
integer(kind=int64),dimension(nTargProc,1+nProjProc)::collisions,collisionsComb
integer(kind=int64),dimension(nSPBins) :: nSPions,nSPionsComb

real*8,allocatable,dimension(:) :: IonEnergy !Initial ion energies once decided
real*8,dimension(nEnergiesJuno) :: IonEnergyJuno !Initial ion energies for Juno
real*8,dimension(nEnergiesNorm) :: IonEnergyNorm !Initial ion energies normally
real*8,dimension(nSulEngBins) :: SulEngBins
real*8,dimension(nChS,nSulEngBins) :: SulVsEng,SulVsEngComb
real*8,dimension(atmosLen) :: Hp,H2p,H2Ex,HpComb,H2pComb,H2ExComb,altitude
real*8,dimension(nChS,atmosLen) :: PhotonsCX,PhotonsDE
real*8,dimension(nChS,atmosLen) :: PhotonsCXComb,PhotonsDEComb
real*8,dimension(nSPBins) :: SPBins,SPvsEng,SIMxsTotvsEng
real*8,dimension(nSPBins) :: dEvsEng,dNvsEng,SIMxsTotxdEvsEng
real*8,dimension(nSPBins) :: SPvsEngComb,SIMxsTotvsEngComb
real*8,dimension(nSPBins) :: dEvsEngComb,dNvsEngComb,SIMxsTotxdEvsEngComb
real*8,dimension(atmosLen,nE2strBins) :: electFwd,electBwd
real*8,dimension(atmosLen,nE2strBins) :: electFwdComb,electBwdComb

character(len=100) filename,files(nOutputFiles) !Output file names
!****************************** Data Declaration *******************************
!* Initial ion enegy input:
data IonEnergyNorm/10.0,50.0,75.0,100.0,125.0,150.0,175.0,200.0,250.0,&
     300.0,350.0,400.0,450.0,500.0,600.0,700.0,800.0,900.0,1000.0,1250.0,&
     1500.0,1750.0,2000.0/
!* Initial ion enegy input from interpoalted JEDI bins:
data IonEnergyJuno/5.312,6.062,6.893,7.759,8.714,9.766,11.140,12.271,13.518,&
     14.892,16.660,18.638,20.851,23.326,24.817,26.403,28.090,29.885,31.892,&
     34.035,36.321,38.761,43.293,48.355,54.009,60.324,69.950,81.112,94.054,&
     109.062,131.160,157.734,189.692,228.125,270.78,312.5/ !Interpolated energies
data files/'ChargeStateDistribution','H+_Prod','H2+_Prod','H2*_Prod',&
     'Collisions','Photons_CX','Photons_DE','Stopping_Power',&
     '2Str_Elect_Fwd','2Str_Elect_Bwd'/
!********************************* Initialize **********************************
energy=0;nTrials=0;trial=0;nEnergies=0
!*******************************************************************************
EnergySwitch=1 !1 for normal energy bins, 2 for Juno energy bins
if(EnergySwitch.eq.1)then !Normal energy bins
  nEnergies=nEnergiesNorm
  allocate(IonEnergy(nEnergies))
  IonEnergy=IonEnergyNorm
elseif(EnergySwitch.eq.2)then !JEDI interpolated energy bins
  nEnergies=nEnergiesJuno
  allocate(IonEnergy(nEnergies))
  IonEnergy=IonEnergyJuno
end if
!******************************** Main Program *********************************
! write(*,*) "What energy [kev] for Elapsed_Times.dat file?"
! read(*,*) energy !Input the initial ion energy
! write(filename,'("../scratch/Jup_Sul/Output/",I0,"/Elapsed_Times.dat")')&
!  energy
! ! write(filename,'("../scratch/Jup_Sul/Output/Trials.dat")')
! open(unit=100,file=filename,status='old')
! do i=1,MaxnTrials
!   read(100,*,end=1000) trial(i)
!   nTrials=nTrials+1
! end do
! 1000 continue
! close(100)
do run=nEnergies,nEnergies
  nTrials=0
  energy=nint(IonEnergy(run))
  write(filename,'("../scratch/Jup_Sul/Output/",I0,"/Elapsed_Times.dat")')&
   energy
  open(unit=100,file=filename,status='old')
  do i=1,MaxnTrials
   read(100,*,end=1000) trial(i)
   nTrials=nTrials+1
  end do
  1000 continue
  close(100)
  write(*,*) 'Reading in and combining files...'
  write(*,*) 'Number of files: ',nTrials,'At an energy of: ',energy,'keV/u.'
!********************************* Initialize **********************************
  nerr=0;start=1;HpComb=0.0;H2pComb=0.0;H2ExComb=0.0;SulVsEngComb=0.0
  collisionsComb=0;PhotonsCXComb=0.0;PhotonsDEComb=0.0;SPvsEngComb=0.0
  SIMxsTotvsEngComb=0.0;dEvsEngComb=0.0;dNvsEngComb=0.0;SIMxsTotxdEvsEngComb=0.0
  nSPionsComb=0;electFwdComb=0.0;electBwdComb=0.0
!********** Open output data files for each set of initial energies ************
  1002 continue
  do n=start,nTrials
!********************************* Initialize **********************************
    nLines=0;SulEngBins=0.0;SulVsEng=0.0;altitude=0.0;Hp=0.0;H2p=0.0;H2Ex=0.0
    collisions=0;PhotonsCX=0.0;PhotonsDE=0.0;SPBins=0.0;SPvsEng=0.0;
    SIMxsTotvsEng=0.0;dEvsEng=0.0;dNvsEng=0.0;SIMxsTotxdEvsEng=0.0;nSPions=0
    electFwd=0.0;electBwd=0.0
!******************************** Read in Files ********************************
    do i=1,nOutputFiles !Open all of the files
      write(filename,'("../scratch/Jup_Sul/Output/",I0,"/",A,"-",I0,".dat")') &
            energy,trim(files(i)),trial(n) !File output name
      filename=trim(filename)
      open(unit=100+i,file=filename,status='old',iostat=err) !Open the files
      if(err.gt.0)then !If there's an error opening the file
        write(*,*) 'File:',n,'Trial:',trial(n),'ERROR! File:',filename !Error note
        start=n+1 !If there's an error, want to go to next trial
        nerr=nerr+1 !Count the number of errors accrued
        goto 1002 !Go to the next trial
      end if !End error opening file if statemen
      do j=1,MaxnLines !Do loop to determine how many lines are in each file
        read(100+i,*,end=1001) !Read the line, go to 1001 when finished
        nLines(i)=nLines(i)+1 !Count the number of lines
      end do !End line counting do loop
  1001 continue !Continue once reached the end of the file
      if(nLines(i).eq.0)then !Make sure a file isn't empty
        start=n+1 !If empty, skip it
        write(*,*) 'File:',n,'Trial:',trial(n),'nLines ERROR! File:',filename
        goto 1002 !Go on to the next trial
      end if !End nLines=0 if statement
      rewind(100+i) !Rewind all files after finding length
    end do !End do loop for all files of a specific trial
    write(*,*) 'File:',n,'Trial:',trial(n)
  !*** Charge state equilibrium fractions
    do i=1,2
      read(101,*) !Skip charge state distribution header
    end do
    do i=2,nSulEngBins
      read(101,F01) SulEngBins(i),(SulVsEng(j,i),j=1,nChS)
    end do
    SulVsEngComb=SulVsEngComb+SulVsEng !Charge state fractions combined
  !*** H production
    do i=1,3
      read(102,*) !Skip H^+ Header
      read(103,*) !Skip H_2^+ Header
      read(104,*) !Skip H_2^* Header
    end do
    read(104,*) !H_2^* has an additional line
    do i=1,atmosLen !Loop through the atmosphere
      read(102,F02) altitude(i),Hp(i) !H^+ production
      read(103,F02) altitude(i),H2p(i) !H_2^+ production
      read(104,F02) altitude(i),H2Ex(i) !H_2^* production
    end do
    HpComb=HpComb+Hp !H^+ production combined
    H2pComb=H2pComb+H2p !H_2^+ production combined
    H2ExComb=H2ExComb+H2Ex !H_2^* production combined
  !*** Collision counters
    read(105,*) !Skip first line
    do i=1,nTargProc
      read(105,F04a) (collisions(i,j),j=1,nProjProc+1)
    end do
    collisionsComb=collisionsComb+collisions
  !*** Photon production
    do i=1,13
      read(106,*) !Photon_CX header
      read(107,*) !Photon_DE header
    end do
    read(106,*) !Photon_CX has an additional line
    do i=1,atmosLen
      read(106,F06) altitude(i),(PhotonsCX(j,i),j=1,nChS) !CX-TI, SC, SC+SS
      read(107,F06) altitude(i),(PhotonsDE(j,i),j=1,nChS) !DE-(SI, DI, TEX)+SPEX
    end do
    PhotonsCXComb=PhotonsCXComb+PhotonsCX
    PhotonsDEComb=PhotonsDEComb+PhotonsDE
  !*** Stopping power
    do i=1,5
      read(108,*) !Stopping power header
    end do
    do i=2,nSPBins
      read(108,F07) SPBins(i),SPvsEng(i),SIMxsTotvsEng(i),dEvsEng(i),&
        dNvsEng(i),SIMxsTotxdEvsEng(i),nSPions(i)
    end do
    SPvsEngComb=SPvsEngComb+SPvsEng
    SIMxsTotvsEngComb=SIMxsTotvsEngComb+SIMxsTotvsEng
    dEvsEngComb=dEvsEngComb+dEvsEng
    dNvsEngComb=dNvsEngComb+dNvsEng
    SIMxsTotxdEvsEngComb=SIMxsTotxdEvsEngComb+SIMxsTotxdEvsEng
    nSPionsComb=nSPionsComb+nSPions
  !*** 2-Stream electrons
    do j=1,nE2strBins
      read(109,F2Str) (electFwd(i,j),i=atmosLen,1,-1)
      read(110,F2Str) (electBwd(i,j),i=atmosLen,1,-1)
    end do
    electFwdComb=electFwdComb+electFwd
    electBwdComb=electBwdComb+electBwd
    do i=1,nOutputFiles !Close all of the files
      close(100+i)
    end do
  end do !End number of trials loop
!********************************** Write Out **********************************
  write(*,*) 'Writing output files...'
  do i=1,nOutputFiles !Open the final combined files
    write(filename,'("./Output/",I0,"/",A,"-Comb.dat")') &
          energy,trim(files(i))
    filename=trim(filename)
    open(unit=200+i,file=filename,status='unknown')
  end do
  norm=nTrials-nerr !Normalization condition to per ion per cm
  write(*,*) norm,nTrials,nerr
  !*** Charge state distribution
  write(201,H01) !Sulfur charge state distribution header
  do i=2,nSulEngBins !Sulfur charge state distribution
    write(201,F01) SulEngBins(i),(SulVsEngComb(j,i)/norm,j=1,nChS)
  end do
  !*** H production
  write(202,H02) !H^+ Header
  write(203,H03) !H_2^+ Header
  write(204,H04) !H_2^* Header
  do i=1,atmosLen !Loop through the atmosphere
    write(202,F02) altitude(i),HpComb(i)/norm !H^+ production
    write(203,F02) altitude(i),H2pComb(i)/norm !H_2^+ production
    write(204,F02) altitude(i),H2ExComb(i)/norm !H_2^* production
  end do
  !*** Collision counters
  write(205,F03) (ProjColl(i),i=1,nProjProc) !Collisions header
  do i=1,nTargProc !Total number of each type of collision
    write(205,F04) TargColl(i),(collisionsComb(i,j),j=1,nProjProc+1),&
      sum(collisionsComb(i,:))
  end do
  write(205,F4) !'--'
  write(205,F04) 'Sum ',(sum(collisionsComb(:,i)),i=1,nProjProc+1),&
    sum(collisionsComb)
  write(205,*) ''
  write(205,F03) (ProjColl(i),i=1,nProjProc) !Collisions precentage header
  do i=1,nTargProc !Total number of each type of collision
    write(205,F05) TargColl(i),&
      (real(collisionsComb(i,j))/real(sum(collisionsComb))*100,j=1,nProjProc+1),&
      real(sum(collisionsComb(i,:)))/real(sum(collisionsComb))*100
  end do
  write(205,F4) !'--'
  write(205,F05) 'Sum ',&
    (real(sum(collisionsComb(:,i)))/real(sum(collisionsComb))*100,i=1,5),&
    real(sum(collisionsComb))/real(sum(collisionsComb))*100
  !*** Photon production
  write(206,N01) !CX note
  write(207,N02) !DE note
  do i=1,2 !Loop through CX and DE headers
    write(205+i,*) !Blank space
    write(205+i,H09) !Initial input header
    write(205+i,*) !Blank space
    write(205+i,H05) !Altitude integrated photon production header
    write(205+i,H06) !Charge state header
  end do
  !* Altitude integrated photon production
  write(206,F06) 2.0,& !CX - TI, SC, SC+SPEX
    (sum(PhotonsCXComb(ChS,:)*2e5)/norm,ChS=1,nChS)
  write(207,F06) 2.0,& !DE - SI+SPEX, DI+SPEX, TEX+SPEX
    (sum(PhotonsDEComb(ChS,:)*2e5)/norm,ChS=1,nChS)
  do i=1,2 !Loop through CX and DE headers
    write(205+i,*) !Blank space
    write(205+i,H07) !Photon production vs. altitude header
    write(205+i,H08) !Charge state header
  end do
  do i=1,atmosLen
    write(206,F06) altitude(i),& !CX - TI, SC, SC+SPEX
     (PhotonsCXComb(ChS,i)/norm,ChS=1,nChS)
    write(207,F06) altitude(i),& !DE - SI+SPEX, DI+SPEX, TEX+SPEX
     (PhotonsDEComb(ChS,i)/norm,ChS=1,nChS)
  end do
  !*** Stopping power
  write(208,H10) !Stopping power header
  do i=2,nSPBins !Loop through ever stopping power bin
    write(208,F07) SPBins(i),&
      SPvsEngComb(i)/norm,&
      SIMxsTotvsEngComb(i)/norm,&
      dEvsEngComb(i)/norm,&
      dNvsEngComb(i)/norm,&
      (SIMxsTotxdEvsEngComb(i))/norm,&
      nSPionsComb(i)
  end do
  !*** 2-Stream electrons
  write(*,*) norm,nTrials,nerr
  do j=1,nE2strBins
    write(209,F2Str) (electFwdComb(i,j)/real(norm),i=atmosLen,1,-1)
    write(210,F2Str) (electBwdComb(i,j)/real(norm),i=atmosLen,1,-1)
  end do
  do i=1,nOutputFiles !Close all of the files
    close(200+i)
  end do
end do !End of energies do loop


end program
