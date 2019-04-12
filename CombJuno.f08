program CombOxySulJuno
!*******************************************************************************
!* Created by Stephen J. Houston 9.4.18
!*******************************************************************************
!* This program combines the output of oxygen and sulfur precipitation from
!* JEDI measurements.
!*******************************************************************************

use, intrinsic :: ISO_FORTRAN_ENV !Used for kind=int64
use formatting !Used for formatting.f08
implicit real*8(a-h,o-z) !i,j,k,l,m,n are assumed to be integers

!*******************************************************************************
integer energy,atmosLen,ChS,run,Oxy,Sul

integer pProc,nProjProc !Number of projectile processes
integer tProc,nTargProc !Number of target processes
parameter(nProjProc=4,nTargProc=7)

parameter(nChSOxy=10,nChSSul=17,atmosLen=1544,nSpec=2) !nSpecies = 2, oxy and sul
parameter(Oxy=1,Sul=2)
parameter(nOutputFiles=8)!,MaxnTrials=1000,MaxnLines=100000)
parameter(nE2strBins=260) !Number of 2 stream bins
real*8,dimension(atmosLen) :: altitude
real*8,dimension(nSpec,atmosLen) :: Hp,H2p,H2Ex
real*8,dimension(nSpec,nChSOxy,atmosLen) :: PhotonsCXOxy,PhotonsDEOxy
real*8,dimension(nSpec,nChSSul,atmosLen) :: PhotonsCXSul,PhotonsDESul
real*8,dimension(nSpec,atmosLen,nE2strBins) :: electFwd,electBwd
!* JEDI variables
real*8,dimension(nSpec) :: Jflux
!*   Jflux - JEDI measured flux intensities converted to [counts/cm^2/s]

!Same as previous variables except has a leading "J"
! real*8,dimension(atmosLen) :: JHp,JH2p,JH2Ex
! real*8,dimension(nChS,atmosLen) :: JPhotonsCX,JPhotonsDE
! real*8,dimension(atmosLen,nE2strBins) :: JelectFwd,JelectBwd

character*100 filename
character*100,dimension(nOutputFiles) :: OxyFiles,SulFiles !I/O file names
character*1 random_number_file
character*4 dumChar
!****************************** Data Declaration *******************************
data OxyFiles/'H+_Prod','H2+_Prod','H2_Excite_Prod','XRay_CX','XRay_DE',&
      '2Str_Elect_Fwd','2Str_Elect_Bwd','Total_Photon_Prod'/ !Oxy input files
data SulFiles/'H+_Prod','H2+_Prod','H2*_Prod','Photons_CX','Photons_DE',&
      '2Str_Elect_Fwd','2Str_Elect_Bwd','Photons_Total'/ !Sul input files
!********************************* Initialize **********************************
altitude=0.0;Hp=0.0;H2p=0.0;H2Ex=0.0;electFwd=0.0;electBwd=0.0
PhotonsCX=0.0;PhotonsDE=0.0
!******************* Open output data files for each species *******************
!********** OXYGEN
write(*,*) 'Opening and reading oxygen preciptation files...'
do i=1,nOutputFiles-1 !Open all of the files
  write(filename,'("../Jup_Oxy/Output/Juno/PJ7_Paper/",A,".dat")') &
        trim(OxyFiles(i))
  filename=trim(filename)
  write(*,*) filename
  open(unit=100+i,file=filename,status='old')
end do
!*** H production
do i=1,3
  read(101,*) !Read the headers of the first four files
  read(102,*)
  read(103,*)
end do
read(103,*)
do i=1,atmosLen !Hydrogen Ionization/Excitation vs. altitude
  read(101,*) altitude(i),(dum,j=1,31),Hp(Oxy,i)
  read(102,*) dum,(dum,j=1,11),H2p(Oxy,i)
  read(103,*) dum,H2Ex(Oxy,i)
end do
!*** Photon production
do i=1,13
  read(104,*) !Photon_CX header
  read(105,*) !Photon_DE header
end do
read(104,*) !Photon_CX has an additional line
do i=1,atmosLen !CX - TI, SC, SC+SPEX !DE - (SI, DI, TEX)+SPEX
  read(104,*) altitude(i),(PhotonsCXOxy(Oxy,j,i),j=1,nChSOxy)
  read(105,*) altitude(i),(PhotonsDEOxy(Oxy,j,i),j=1,nChSOxy)
end do
!*** 2-Stream electrons
do j=1,nE2strBins !Electron production rate for 2-stream
  read(106,F2Str) (electFwd(Oxy,i,j),i=atmosLen,1,-1)
  read(107,F2str) (electBwd(Oxy,i,j),i=atmosLen,1,-1)
end do
do i=1,nOutputFiles-1
  close(100+i) !close files
end do
!********** SULFUR
write(*,*) 'Opening and reading sulfur preciptation files...'
do i=1,nOutputFiles-1 !Open all of the files
  write(filename,'("./Output/Juno/PJ7_Paper/",A,".dat")') &
        trim(SulFiles(i))
  filename=trim(filename)
  write(*,*) filename
  open(unit=100+i,file=filename,status='old')
end do
!*** H production
do i=1,3
  read(101,*) !Skip H^+ Header
  read(102,*) !Skip H_2^+ Header
  read(103,*) !Skip H_2^* Header
end do
read(103,*) !H_2^* has an additional line
do i=1,atmosLen !Loop through the atmosphere
  read(101,*) altitude(i),Hp(Sul,i) !H^+ production
  read(102,*) altitude(i),H2p(Sul,i) !H_2^+ production
  read(103,*) altitude(i),H2Ex(Sul,i) !H_2^* production
end do
!*** Photon production
do i=1,13
  read(104,*) !Photon_CX header
  read(105,*) !Photon_DE header
end do
read(104,*) !Photon_CX has an additional line
do i=1,atmosLen !CX - TI, SC, SC+SPEX !DE - (SI, DI, TEX)+SPEX
  read(104,*) altitude(i),(PhotonsCXSul(Sul,j,i),j=1,nChSSul)
  read(105,*) altitude(i),(PhotonsDESul(Sul,j,i),j=1,nChSSul)
end do
!*** 2-Stream electrons
do j=1,nE2strBins
  read(106,F2Str) (electFwd(Sul,i,j),i=atmosLen,1,-1)
  read(107,F2Str) (electBwd(Sul,i,j),i=atmosLen,1,-1)
end do
do i=1,nOutputFiles-1 !Close all of the files
  close(100+i)
end do
!*********************** Write Out combined Productions ************************
write(*,*) 'Writing combined oxygen and suflur preciptation files...'
do i=1,nOutputFiles-1 !Open all of the files
  write(filename,'("./Output/Juno/PJ7_Paper/OxySulComb/",A,"-Comb.dat")') &
        trim(SulFiles(i))
  filename=trim(filename)
  write(*,*) filename
  open(unit=200+i,file=filename)
end do
!*** H production
write(201,H02) !H^+ Header
write(202,H03) !H_2^+ Header
write(203,H04) !H_2^* Header
do i=1,atmosLen !Loop through the atmosphere
  write(201,F02) altitude(i),Hp(Oxy,i)+Hp(Sul,i) !H^+ production
  write(202,F02) altitude(i),H2p(Oxy,i)+H2p(Sul,i) !H_2^+ production
  write(203,F02) altitude(i),H2Ex(Oxy,i)+H2Ex(Sul,i) !H_2^* production
end do
!*** Photon production
! write(204,N01) !CX note
! write(205,N02) !DE note
! do i=1,2 !Loop through CX and DE headers
!   if(i.le.2)then
!     write(203+i,*) !Blank space
!     write(203+i,H09) !Initial input header
!     write(203+i,*) !Blank space
!   end if
!   write(203+i,H05) !Altitude integrated photon production header
!   write(203+i,H06) !Charge state header
! end do
! do i=1,atmosLen !CX - TI, SC, SC+SPEX !DE - (SI, DI, TEX)+SPEX
!   read(104,*) altitude(i),(PhotonsCXSul(Sul,j,i),j=1,nChSSul)
!   read(105,*) altitude(i),(PhotonsDESul(Sul,j,i),j=1,nChSSul)
! end do
!*** 2-Stream electrons
do j=1,nE2strBins
  write(206,F2Str) ((electFwd(Oxy,i,j)+electFwd(Sul,i,j)),i=atmosLen,1,-1)
  write(207,F2Str) ((electBwd(Oxy,i,j)+electBwd(Sul,i,j)),i=atmosLen,1,-1)
end do
do i=1,nOutputFiles-1 !Close all of the files
  close(100+i)
end do

end program
