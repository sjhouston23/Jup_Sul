program LymanWernerSpectrum

implicit real*8 (a-h,o-z)

integer zbins, lambdabins, dz, IonEnergy
parameter (zbins=1544, lambdabins=308)
real*8 Alt(zbins), LyDE(zbins), LyC(zbins), LyI(zbins), WDE(zbins), WI(zbins), LyAlphaH2(zbins), LyAlphaH(zbins)
real*8 LWProd(zbins), Wavelength(lambdabins+1), IntensityRatio(lambdabins), LWSpectrum(lambdabins), LWIntensity(lambdabins)
real*8 WavelengthM(lambdabins)
real*8 dummy, CD_H2(zbins), CD_He(zbins), CD_CH4(zbins), CD_H(zbins), TotCD(zbins)
real*8 tau(zbins), sigCH4(lambdabins), LWProdTemp(zbins), ColRat(24)
real*8 Lyman(zbins), Werner(zbins), LymanTemp(zbins), WernerTemp(zbins)
real*8 LymanSpectrum(lambdabins), LymanIntensity(lambdabins),  WernerSpectrum(lambdabins), WernerIntensity(lambdabins)
real*8 LymanA(zbins), LWRay(lambdabins)
real*8 Angle
!real*8 CD_M1I(zbins) ,CD_M2I(zbins) ,CD_M3I(zbins)
real Eion(24), Eelec(6)
!****************************************************************
!   LyDE - Lyman Bands from Direct Excitation
!   LyC - Lyman Bands from Cascading
!   LyI - Lyman Bands from Ion Excitation
!   WDE - Werner Bands from Direct Excitation
!   WI - Werner Bands from Ion Excitation
!****************************************************************v
character(len=100) filename, filenameSpec
data Eion/10.0,25.0,50.0,75.0,100.0,125.0,150.0,175.0,200.0,250.0,&
     300.0,350.0,400.0,450.0,500.0,600.0,700.0,800.0,900.0,1000.0,1250.0,&
     1500.0,1750.0,2000.0/
data Eelec/5.0, 10.0, 20.0, 50.0, 100.0, 200.0/

! open(unit=1000,file='COLDENS_test2.DAT', status='unknown')
! open(unit=1000,file='./Atmosphere/Input/JunoColumnDensity_2km.dat',status='old')
open(unit=1000,file='./Atmosphere/Input/JupiterColumnDensity_2km_EqMixing.dat',status='old')
! open(unit=1000,file='GerardMethaneInterpCD.dat', status='unknown')
open(unit=1001,file='./2Stream/output/airglow/NormalizedSpectrum.dat',status='old')
! open(unit=1002,file='./2Stream/output/airglow/colorratio.dat',status='unknown')
! open(unit=1003,file='./Opac.dat',status='unknown')
!open(unit=1004,file='./JunoData/Output/Airglow/LWSpectrumV1.dat',status='old')
!open(unit=1004,file='./Prod5keVbeam.dat',status='unknown')

pi=4.0d0*atan(1.0d0)
read(1000,*)
do i=zbins,1,-1 !This file is from 3000->-88 km, the others are -88->2998 km
  ! read(1000,*) dummy, CD_H2(i), CD_He(i), CD_CH4(i), CD_H(i), TotCD(i)
 read(1000,*) dummy,dummy1,dummy,CD_CH4(i)
end do

do i=1,lambdabins
  read(1001,*) Wavelength(i), IntensityRatio(i), sigCH4(i)
!  if (i.lt.5.or.i.gt.304) write(*,*) Wavelength(i), IntensityRatio(i), sigCH4(i)
end do
Wavelength=Wavelength/10.0 !Convert from Angstrom to nm.
Wavelength(lambdabins+1)=1700.0d0
do i=1,lambdabins
  WavelengthM(i)=Wavelength(i)+0.5*(Wavelength(i+1)-Wavelength(i))
end do

close(1000)
close(1001)
! do i=1,24
  angle=80.0
  CR=0.0d0
  CRNum=0.0d0
  CRDenom=0.0d0
  LWSpectrum=0.0d0
  LWProd=0.0d0
  LWIntensity=0.0d0
  Lyman=0.0d0
  Werner=0.0d0
  LymanTemp=0.0d0
  WernerTemp=0.0d0
  LymanSpectrum=0.0d0
  WernerSpectrum=0.0d0
  LymanIntensity=0.0d0
  WernerIntensity=0.0d0
  TotalRayleigh=0.0d0
  ! IonEnergy = int(Eelec(i))!int(Eion(i))
  ! write(filename,'("./Output/Airglow/ME/LWSpectrum",I0,".dat")') IonEnergy
  ! write(filenameSpec,'("./Output/Airglow/ME/UVSpectrum",I0,"-",I0,"eq.dat")') IonEnergy,int(angle)
  IonEnergy = 111!int(Eion(i))
  write(filename,'("./Output/Airglow/SE/LWSpectrum",I0,"-eq.dat")') IonEnergy
  write(filenameSpec,'("./Output/Airglow/UVSpectra/spec",I0,"-",I0,"eq.dat")') IonEnergy,int(angle)
!J  IonEnergy=23 !Doesn't matter for JEDI data.
!J  write(filename,'("./JunoData/Output/Airglow/LWSpectrumv1.dat")')
!J  write(filenameSpec,'("./JunoData/Output/Airglow/Spectrum60v1.dat")')
  filename=trim(filename)
  filenameSpec=trim(filenameSpec)
  open(unit=100+i,file=filename,status='old')
  open(unit=200+i,file=filenameSpec,status='unknown')
  ! ConversionFactor=IonEnergy*1000*16*1.60218E-19*0.5*1000000*10000*1000
  ! write(*,*) ConversionFactor
  read(100+i,*)
  do dz=1,zbins
    read(100+i,*) Alt(dz), LyDE(dz), LyC(dz), LyI(zbins), WDE(dz), WI(zbins), LyAlphaH2(dz), LyAlphaH(dz)
  end do
  do dz=1,zbins
    LWProd(dz)=LyDE(dz)+LyC(dz)+LyI(dz)+WDE(dz)+WI(dz)
    Lyman(dz)=LyDE(dz)+LyC(dz)+LyI(dz)
    Werner(dz)=WDE(dz)+WI(dz)
    LymanA(dz)=LyAlphaH(dz)+LyAlphaH2(dz)
  end do
  write(*,*) IonEnergy,sum(LyAlphaH2)*2e5
!  write(*,*) ''
!  write(*,'(5ES10.3)') sum(LyDe)*2e11, sum(LyC)*2e11, sum(WDe)*2e11, sum(LyAlphaH2)*2e11, sum(LyAlphaH)*2e11
!  write(*,'(ES10.3)') (sum(LyDe)+sum(LyC)+sum(WDe)+sum(LyAlphaH2)+sum(LyAlphaH))*2e11
!  write(*,*) ''
  do j=1,lambdabins
    tau=0.0d0
    do k=1,zbins
      tau(k)=sigCH4(j)*CD_CH4(k)/cos(angle*pi/180.0d0)
      LWIntensity(j)=LWIntensity(j)+LWProd(k)*IntensityRatio(j)*exp(-tau(k))*2e5/(Wavelength(j+1)-Wavelength(j))
      LWRay(j)=LWRay(j)+(LWProd(k)+LymanA(k))*IntensityRatio(j)*exp(-tau(k))*2e5/(Wavelength(j+1)-Wavelength(j))
      LymanIntensity(j)=LymanIntensity(j)+Lyman(k)*IntensityRatio(j)*exp(-tau(k))/(Wavelength(j+1)-Wavelength(j))
      WernerIntensity(j)=WernerIntensity(j)+Werner(k)*IntensityRatio(j)*exp(-tau(k))/(Wavelength(j+1)-Wavelength(j))
      LWProdTemp(k)=LWProd(k)*exp(-tau(k))
      LymanTemp(k)=Lyman(k)*exp(-tau(k))
      WernerTemp(k)=Werner(k)*exp(-tau(k))
!      write(*,*) exp(-tau(k)), tau(k), Lyman(k), LymanTemp(k)
!      if(j.eq.1) then
!        write(1003,*) alt(k), exp(-tau(k))
!        write(1004,*) alt(k)/1e5, LWProdTemp(k), LWProd(k)
!E        write(200+i,*) alt(k)/1e5, LWProdTemp(k), LWProd(k) !Used for Monoenergetic Electron Beams
!      end if
    end do
    LWSpectrum(j)=sum(LWProd)*IntensityRatio(j)*2e5/(Wavelength(j+1)-Wavelength(j))
    LymanSpectrum(j)=sum(Lyman)*IntensityRatio(j)
    WernerSpectrum(j)=sum(Werner)*IntensityRatio(j)
    write(200+i,*) Wavelength(j), LWSpectrum(j), LWIntensity(j)
!    write(*,*) Wavelength(j), LWIntensity(j)/LWSpectrum(j), IonEnergy, sum(LWProdTemp)/sum(LWProd)
    TotalRayleigh=TotalRayleigh+LWRay(j)*(Wavelength(j+1)-Wavelength(j))
  end do
!  write(*,'(2ES10.2)') sum(LymanIntensity)*2e11/ConversionFactor, sum(WernerIntensity)*2e11/ConversionFactor
!  write(*,'(2ES10.2)') sum(Lyman)*2e11/ConversionFactor, sum(Werner)*2e11/ConversionFactor
!  write(*,'(3ES10.2)') sum(LWSpectrum), sum(LWProd)*2e5, sum(LWIntensity)
!  write(*,*) sum(LWProd)*2e5, sum(LWSpectrum), sum(LWIntensity), sum(LWIntensity)/sum(LWSpectrum)
  close(100+i)
  close(200+i)
  do j=258,283 !1549.6 - 1619.2 Angstrom
    CRNum=CRNum+LWIntensity(j)*(Wavelength(j+1)-Wavelength(j))!LWIntensity(j)
  end do
  do j=137,169 !1229.2 - 1300.8 Angstrom
    CRDenom=CRDenom+LWIntensity(j)*(Wavelength(j+1)-Wavelength(j))!LWIntensity(j)
  end do
  CR=CRNum/CRDenom
  ColRat(i)=CR
  write(*,*) filenameSpec
  write(*,*) 'Color Ratio:                         ', CR
  write(*,*) 'Number of kRayleighs [10^9 photons]: ', TotalRayleigh/1e9
!I  write(*,*) 'Number of kRayleighs [input 1 mW/m^2]: ', TotalRayleigh*0.5d0/(Eion(i)*1.60218*16), Eion(i)
  write(*,*) 'Number of kRayleighs [input 1 mW/m^2]: ', TotalRayleigh/1e9, IonEnergy!/(Eelec(i)*1.60218*0.5), Eelec(i)
!  write(*,*) 'Lyman Alpha Production:              ', sum(LymanA)*2e5
! end do
! write(filename,"('./Output/Airglow/SE/ColorRatio',I0,'eq.dat')") int(angle)
! open(unit=300,file=filename)
! do i=1,24
!  write(300,*) EIon(i), ColRat(i)
! end do
!close(300)
!write(*,*) ColRat
!write(*,*) ColRat
end
