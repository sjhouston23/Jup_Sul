program JEDISulPrecip
!*******************************************************************************
!* Created by Stephen J. Houston 9.4.18
!*******************************************************************************
!* This program reads in a JEDI spectrum (*.d2s) and normalizes it based on the
!* energy and the energy bin width. It then multiplies the new input ion flux
!* by normalized (1 input ion/cm^2/s), previously calculated results that
!* correspond to the JEDI energy bins ~(6,8,10,15,24,30,39,61,109,228 keV).
!* Some of the higher energy JEDI bins overlap - I treat them as if they don't.
!* To change the number of interpolated data points, one would need to change
!* the variables number_of_energies, the Eion data points, and then go into
!* the JEDIInterpolator subroutine and make changes.
!*******************************************************************************

use, intrinsic :: ISO_FORTRAN_ENV !Used for kind=int64
use formatting !Used for formatting.f08
implicit real*8(a-h,o-z) !i,j,k,l,m,n are assumed to be integers

!*******************************************************************************
integer energy,atmosLen,ChS,run

integer pProc,nProjProc !Number of projectile processes
integer tProc,nTargProc !Number of target processes
parameter(nProjProc=4,nTargProc=7)

parameter(nChS=17,atmosLen=1544,nEng=34)
parameter(nOutputFiles=5)!,MaxnTrials=1000,MaxnLines=100000)
real*8 Eion(nEng) !Ion energies
! parameter(nSulEngBins=2000,nSPBins=2000)

! integer trial(MaxnTrials),nLines(nOutputFiles) !Number of trials/lines in a file
! integer(kind=int64),dimension(nTargProc,1+nProjProc)::collisions
! integer(kind=int64),dimension(nSPBins) :: nSPions

! real*8,dimension(nSulEngBins) :: SulEngBins
real*8,dimension(atmosLen) :: altitude
! real*8,dimension(nSPBins) :: SPBins

! real*8,dimension(nEng,nChS,nSulEngBins) :: SulVsEng
real*8,dimension(nEng,atmosLen) :: Hp,H2p,H2Ex
real*8,dimension(nEng,nChS,atmosLen) :: PhotonsCX,PhotonsDE
! real*8,dimension(nEng,nSPBins) :: SPvsEng,SIMxsTotvsEng
! real*8,dimension(nEng,nSPBins) :: dEvsEng,dNvsEng,SIMxsTotxdEvsEng
!* JEDI variables
real*8,dimension(nEng) :: Jflux
!*   Jflux - JEDI measured flux intensities converted to [counts/cm^2/s]

!Same as previous variables except has a leading "J"
! real*8,dimension(nChS,nSulEngBins) :: JSulVsEng
real*8,dimension(atmosLen) :: JHp,JH2p,JH2Ex
real*8,dimension(nChS,atmosLen) :: JPhotonsCX,JPhotonsDE
! real*8,dimension(nSPBins) :: JSPvsEng,JSIMxsTotvsEng
! real*8,dimension(nSPBins) :: JdEvsEng,JdNvsEng,JSIMxsTotxdEvsEng

character(len=100) filename,files(nOutputFiles) !Output file names
character(len=1) random_number_file
character(len=4) dumChar
character(len=10) date,version
character(len=12) time
!****************************** Data Declaration *******************************
!* Initial ion enegy input:
!data Eion/10.625,15.017,20.225,29.783,46.653,59.770,77.522,120.647,218.125,&
!          456.250/ !Original Juno energy bins from JEDI (keV/u, u=16).
data Eion/5.312,6.062,6.893,7.759,8.714,9.766,11.140,12.271,13.518,14.892,&
     16.660,18.638,20.851,23.326,24.817,26.403,28.090,29.885,31.892,34.035,&
     36.321,38.761,43.293,48.355,54.009,60.324,69.950,81.112,94.054,109.062,&
     131.160,157.734,189.692,228.125/ !Interpolated energies in KeV/u (u=32)
data files/'H+_Prod','H2+_Prod','H2*_Prod','Photons_CX','Photons_DE'/!,&
     ! '2Str_Elect_Fwd','2Str_Elect_Bwd'/ !Filenames
!* Width of JEDI energy bins: (May eventually need to be adjusted)
!data Jebins/66.0,71.0,105.0,216.0,346.0,251.0,300.0,880.0,2280.0,5340.0/
!********************************* Initialize **********************************
pi=4.0d0*atan(1.0d0);Jflux=0.0
!*************************** Call JEDI Interpolator ****************************
write(*,*)
write(*,*) "What is the name of the JEDI ion spectrum file you wish to open?"
write(*,*) "Don't include the extension in the file name (e.g. if you want to &
 analyze PJ7-1.ds2, type PJ7-1)."
read(*,*) version
!write(version,'("v1")') !Filename of a JEDI spectrum (.d2s file)
call JEDIInterpolator(version,Jflux)
!********************************* Initialize **********************************
energy=0;altitude=0.0;Hp=0.0;H2p=0.0;H2Ex=0.0!;prode2stF=0.0;prode2stB=0.0
PhotonsCX=0.0;PhotonsDE=0.0
!********** Open output data files for each set of initial energies ************
write(*,*) 'Opening sulfur preciptation files for all energies...'
do run=1,nEng !Loop through each initial ion energy
  energy=nint(Eion(run))
  do i=1,nOutputFiles !Open all of the files
    write(filename,'("./Output/Juno/",I0,"/",A,"-Comb.dat")') &
          energy,trim(files(i))
    filename=trim(filename)
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
    write(101,F02) altitude(i),Hp(run,i) !H^+ production
    write(102,F02) altitude(i),H2p(run,i) !H_2^+ production
    write(103,F02) altitude(i),H2Ex(run,i) !H_2^* production
  end do
!*** Photon production
  do i=1,13
    read(104,*) !Photon_CX header
    read(105,*) !Photon_DE header
  end do
  read(104,*) !Photon_CX has an additional line
  do i=1,atmosLen !CX - TI, SC, SC+SPEX !DE - (SI, DI, TEX)+SPEX
    read(104,F06) altitude(i),(PhotonsCX(run,j,i),j=1,nChS)
    read(105,F06) altitude(i),(PhotonsDE(run,j,i),j=1,nChS)
  end do
  do i=1,nOutputFiles
    close(100+i)
  end do
end do
!********************************* Initialize **********************************
JHp=0.0;JH2p=0.0;JH2Ex=0.0;JPhotonsCX=0.0;JPhotonsDE=0.0
!************************* Calculate JEDI Productions **************************
write(*,*) 'Calculating JEDI production rates...'
do run=1,nEng !Loop through every energy bin
  do i=1,atmosLen !Loop through entire atmosphere
    JHp(i)=JHp(i)+Hp(run,i)*Jflux(run) !Total H+
    JH2p(i)=JH2p(i)+H2p(run,i)*Jflux(run) !Total H2+
    JH2Ex(i)=JH2Ex(i)+H2Ex(run,i)*Jflux(run) !Total H2*
    do j=1,nChS
      JPhotonsCX(j,i)=JPhotonsCX(j,i)+PhotonsCX(run,j,i)*Jflux(run)
      JPhotonsDE(j,i)=JPhotonsDE(j,i)+PhotonsDE(run,j,i)*Jflux(run)
    end do
  end do !End atmsophere loop
end do !End energy bin loop
!************************* Write Out JEDI Productions **************************
write(*,*) 'Writing output files...'
do i=1,nOutputFiles
  write(filename,'("./Output/Juno/",A,"/",A,".dat")') &
        trim(version),trim(files(i))
  filename=trim(filename)
  open(unit=200+i,file=filename,status='unknown')
end do
!*** H production
write(201,H02) !H^+ Header
write(202,H03) !H_2^+ Header
write(203,H04) !H_2^* Header
do i=1,atmosLen !Loop through the atmosphere
  write(201,F02) altitude(i),JHp(i) !H^+ production
  write(202,F02) altitude(i),JH2p(i) !H_2^+ production
  write(203,F02) altitude(i),JH2Ex(i) !H_2^* production
end do
!*** Photon production
write(204,N01) !CX note
write(205,N02) !DE note
do i=1,2 !Loop through CX and DE headers
  write(203+i,H06) !Charge state header
end do
!* Altitude integrated photon production
write(204,F06) 2.0,& !CX - TI, SC, SC+SPEX
  (sum(JPhotonsCX(ChS,:))*2.0e5/norm,ChS=1,nChS)
write(205,F06) 2.0,& !DE - SI+SPEX, DI+SPEX, TEX+SPEX
  (sum(JPhotonsDE(ChS,:))*2.0e5/norm,ChS=1,nChS)
do i=1,2 !Loop through CX and DE headers
  write(203+i,*) !Blank space
  write(203+i,H07) !Photon production vs. altitude header
  write(203+i,H08) !Charge state header
end do
do i=1,atmosLen
  write(204,F06) altitude(i),& !CX - TI, SC, SC+SPEX
   (JPhotonsCX(ChS,i)/norm,ChS=1,nChS)
  write(205,F06) altitude(i),& !DE - SI+SPEX, DI+SPEX, TEX+SPEX
   (JPhotonsDE(ChS,i)/norm,ChS=1,nChS)
end do
do i=1,nOutputFiles
  close(200+i)
end do

end program
