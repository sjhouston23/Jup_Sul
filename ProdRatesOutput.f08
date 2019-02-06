program ProductionRates
!*******************************************************************************
!* Created by Stephen J. Houston 2.4.19
!*******************************************************************************
!* This program loops through all of the output files and combines the
!* production rates of various energies into output files by charge state for
!* easier plotting.
!*******************************************************************************

use formatting !Formatting module to avoid cluttering the end of the program
implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer Eng,ChS,Alt,AtmosLen

parameter(nFiles=10,nChS=17,AtmosLen=1544,nEnergies=9)

integer,dimension(nEnergies) :: Eion

real*8 Altitude(AtmosLen)
real*8,dimension(AtmosLen,nEnergies,nChS,2) :: Production

character(len=100) filename,files(2)
!****************************** Data Declaration *******************************
data Eion/10,50,75,100,200,300,500,1000,2000/
!******************************** Main Program *********************************
do Eng=1,nEnergies
  write(filename,'("./OutputNorm2/",I0,"/Photons_CX-Comb.dat")') Eion(Eng)
  open(unit=100,file=trim(filename),status='unknown') !Open X-Ray DE
  write(filename,'("./OutputNorm2/",I0,"/Photons_DE-Comb.dat")') Eion(Eng)
  open(unit=101,file=trim(filename),status='unknown') !Open X-Ray CX
  do i=1,13
    read(100,*) !Skip the header
    read(101,*) !Skip the header
  end do
  read(100,*) !Has an additional header line
  do i=1,AtmosLen
    read(100,*) Altitude(i),(Production(i,Eng,ChS,1),ChS=1,nChS) !CX production
    read(101,*) Altitude(i),(Production(i,Eng,ChS,2),ChS=1,nChS) !DE production
  end do
  close(100)
  close(101)
end do
do ChS=1,nChS
  write(filename,'("./OutputNorm2/TotalProductions/XRay_CX_",I0,".dat")') ChS-1
  open(unit=100,file=trim(filename),status='unknown') !Open X-Ray DE
  write(filename,'("./OutputNorm2/TotalProductions/XRay_DE_",I0,".dat")') ChS-1
  open(unit=101,file=trim(filename),status='unknown') !Open X-Ray CX
  write(100,1000) 'ΔAlt [km]', (Eion(Eng),Eng=1,nEnergies) !CX
  write(101,1000) 'ΔAlt [km]', (Eion(Eng),Eng=1,nEnergies) !DE
  write(100,1001) 2.0,(sum(Production(:,Eng,ChS,1))*2e5,Eng=1,nEnergies) !CX
  write(101,1001) 2.0,(sum(Production(:,Eng,ChS,2))*2e5,Eng=1,nEnergies) !DE
  write(100,*)
  write(101,*)
  write(100,1002) ' Alt [km] ', (Eion(Eng),Eng=1,nEnergies) !CX
  write(101,1002) ' Alt [km] ', (Eion(Eng),Eng=1,nEnergies) !DE
  do Alt=1,AtmosLen
    write(100,1001) Altitude(Alt),(Production(Alt,Eng,ChS,1),Eng=1,nEnergies)
    write(101,1001) Altitude(Alt),(Production(Alt,Eng,ChS,2),Eng=1,nEnergies)
  end do
  close(100)
  close(101)
end do
write(filename,'("./OutputNorm2/LaTeX/CX_Prod.dat")')
open(unit=100,file=trim(filename))
write(filename,'("./OutputNorm2/LaTeX/DE_Prod.dat")')
open(unit=101,file=trim(filename))
write(100,1003)
write(101,1003)
write(100,1004)
write(101,1004)
write(100,1005)
write(101,1005)
write(100,1006)
write(101,1006)
write(100,1004)
write(101,1004)
do ChS=1,nChS-1
  if(ChS.eq.1)then
    write(100,1007) (sum(Production(:,Eng,ChS,1))*2e5,Eng=1,5)
    write(101,1007) (sum(Production(:,Eng,ChS,2))*2e5,Eng=1,5)
  elseif(ChS.eq.2)then
    write(100,1008) (sum(Production(:,Eng,ChS,1))*2e5,Eng=1,5)
    write(101,1008) (sum(Production(:,Eng,ChS,2))*2e5,Eng=1,5)
  elseif(ChS.eq.3)then
    write(100,1009) (sum(Production(:,Eng,ChS,1))*2e5,Eng=1,5)
    write(101,1009) (sum(Production(:,Eng,ChS,2))*2e5,Eng=1,5)
  else
    write(100,1010) ChS-1,(sum(Production(:,Eng,ChS,1))*2e5,Eng=1,5)
    write(101,1010) ChS-1,(sum(Production(:,Eng,ChS,2))*2e5,Eng=1,5)
  end if
  if(ChS.eq.nChS-1)then
    write(100,1004)
    write(101,1004)
    write(100,1004)
    write(101,1004)
    write(100,1005)
    write(101,1005)
    write(100,1011)
    write(101,1011)
  end if
end do
do ChS=1,nChS-1
  if(ChS.eq.1)then
    write(100,1017) (sum(Production(:,Eng,ChS,1))*2e5,Eng=6,9)
    write(101,1017) (sum(Production(:,Eng,ChS,2))*2e5,Eng=6,9)
  elseif(ChS.eq.2)then
    write(100,1018) (sum(Production(:,Eng,ChS,1))*2e5,Eng=6,9)
    write(101,1018) (sum(Production(:,Eng,ChS,2))*2e5,Eng=6,9)
  elseif(ChS.eq.3)then
    write(100,1019) (sum(Production(:,Eng,ChS,1))*2e5,Eng=6,9)
    write(101,1019) (sum(Production(:,Eng,ChS,2))*2e5,Eng=6,9)
  else
    write(100,1020) ChS-1,(sum(Production(:,Eng,ChS,1))*2e5,Eng=6,9)
    write(101,1020) ChS-1,(sum(Production(:,Eng,ChS,2))*2e5,Eng=6,9)
  end if
  if(ChS.eq.nChS-1)then
    write(100,1004)
    write(101,1004)
    write(100,1012)
    write(101,1012)
    write(100,1013)
    write(101,1013)
    write(100,1014)
    write(101,1014)
  end if
end do

1000 format(A10,9(1x,'x',I8))
1001 format(F10.2,9(2x,ES8.2E2))
1002 format(A10,9(2x,I8))
1003 format('\begin{table}',/,TR4,'\centering',/,TR4,&
'\begin{tabular}{l|c|c|c|c|c}')
1004 format(TR4,'\hline')
1005 format(TR4,'& Energy & & & & \\')
1006 format(TR4,'Ion Charge State & 10 keV/u & 50 & 75 & 100 & 200 \\')
1007 format(TR4,'S        ',5(' & ',ES8.2E2),' \\')
1008 format(TR4,'S$^+$    ',5(' & ',ES8.2E2),' \\')
1009 format(TR4,'S$^{ ++}$',5(' & ',ES8.2E2),' \\')
1010 format(TR4,'S$^{',I2,'+}$',5(' & ',ES8.2E2),' \\')
1011 format(TR4,'Ion Charge State & 300 keV/u & 500 & 1000 & 2000 & \\')
1012 format(TR4,'\end{tabular}')
1013 format(TR4,'\caption{Altitude integrated X-ray production [cm$^{-2}$&
& s$^{-1}$] from direct excitation collisions (i.e. SI+SPEX, DI+SPEX, TEX+SPEX)&
& for sulfur with incident ion energies between 10 and 25000 keV/u with no&
& opacity effects considered. Everything has been normalized to a single&
& incident ion/cm$^2$/s.}')
1014 format(TR4,'\label{tab:SulDEProd}',/,'\end{table}')
1017 format(TR4,'S        ',4(' & ',ES8.2E2),' & \\')
1018 format(TR4,'S$^+$    ',4(' & ',ES8.2E2),' & \\')
1019 format(TR4,'S$^{ ++}$',4(' & ',ES8.2E2),' & \\')
1020 format(TR4,'S$^{',I2,'+}$',4(' & ',ES8.2E2),' & \\')



end program
