program IonComb
!*******************************************************************************
!* Created by Stephen J. Houston 2.25.18
!*******************************************************************************
!* This code combines the H+ and H2+ production rate profiles into one and
!* guesses at an ion production rate for He+, CH4+, CH3+, CH2+, CH+, and C+.
!*******************************************************************************

implicit real*8(a-h,o-z)

!*******************************************************************************
integer energy,atmosLen,run,eng,alt

parameter(nFiles=2,nIonEnergies=24,atmosLen=1544,nProc=36)

real*8 Eion(nIonEnergies)
real*8,dimension(atmosLen) :: Altitude,H2,H,He,CH4,CH3,CH2,CH,C !Ion species


character(len=100) filename,filenames(nFiles)

!****************************** Data Declaration *******************************
data Eion/10.0,25.0,50.0,75.0,100.0,125.0,150.0,175.0,200.0,250.0,300.0,350.0,&
     400.0,450.0,500.0,600.0,700.0,800.0,900.0,1000.0,1250.0,1500.0,1750.0,&
     2000.0/
data filenames/'H2+_Prod','H+_Prod'/
!********************************* Initialize **********************************

!******************************** Main Program *********************************
do eng=1,nIonEnergies
  !*** Initialize the ion production rates
  H2=0.0;H=0.0;He=0.0;CH4=0.0;CH3=0.0;CH2=0.0;CH=0.0;C=0.0
  energy=int(Eion(eng))
  do i=1,nFiles
    write(filename,'("./Output/",I0,"keV/",A,"_Comb.dat")') &
          energy,trim(filenames(i))
    open(100+i,file=filename,status='old')
    do j=1,3
      read(100+i,*) !Skip the header
    end do
  end do
  do alt=1,atmosLen !Read in H2+ and H+ production rates
    read(101,*) Altitude(alt),(dum,i=1,11),H2(alt)
    ! write(*,*) Altitude(alt),H2(alt)
    read(102,*) (dum,i=1,32),H(alt)
    if(H2(alt).lt.1e-20)H2(alt)=1e-20
    if(H(alt).lt.1e-20)H(alt)=1e-20
    He(alt)=1.0e-20
    CH4(alt)=1.0e-20
    CH3(alt)=1.0e-20
    CH2(alt)=1.0e-20
    CH(alt)=1.0e-20
    C(alt)=1.0e-20
  end do
  close(101)
  close(102)
  write(filename,"('./Output/',I0,'keV/primprod.dat')") energy
  open(200,file=filename)
  write(200,2000) !Header
  do alt=1,atmosLen
    write(200,2001) Altitude(alt),H2(alt),H(alt),He(alt),CH4(alt),CH3(alt),&
    CH2(alt),CH(alt),C(alt)
  end do
  close(200)
end do

2000 format(' Alt[km]',5x,'H2+','[cm^-3]','H+',8x,'He+',6x,'CH4+',6x,'CH3+',6x,'CH2+',&
            6x,'CH+',8x,'C+')
2001 format(F8.2,8(2x,ES8.2))











end program
