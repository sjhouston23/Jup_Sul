program OtherProcesses
!*******************************************************************************
!* Created by Stephen J Houston 12.04.18
!*******************************************************************************
!* This program is designed to read in sulfur integral cross-section data
!* that neither CX's or removes electrons and add them together to look at the
!* total cross-section from those processes.
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer ChS,Eng,Proc

parameter(nChS=17,nProc=13,nEng=2000)

real*8,dimension(nChS,nEng) :: TotalXS
real*8,dimension(nProc,nChS,nEng) :: ProcessXS

character(len=100) :: files(nProc),filename

!****************************** Data Declaration *******************************

data files/'SI','SI+SPEX','SI+DPEX','TI+SS','DI','DI+SPEX','DI+DPEX','DCAI+SS',&
           'SC+SS','DC+DS','TEX','TEX+SPEX','TEX+DPEX'/

!******************************** Main Program *********************************
ProcessXS=0.0;TotalXS=0.0
do Proc=1,nProc
  write(filename,'("./SIMXSInterp/",A,".dat")') trim(files(Proc))
  open(unit=100+Proc,file=trim(filename),status='old')
  read(100+Proc,*) !Skip the header
  do Eng=1,nEng
    read(100+Proc,1000) dum,(ProcessXS(Proc,ChS,Eng),ChS=1,nChS)
  end do
end do
do Proc=1,nProc
  close(100+Proc)
end do
TotalXS=sum(ProcessXS,dim=1)
open(unit=200,file='./SIMXSInterp/OtherProcesses.dat')
do Eng=1,nEng
  write(200,1000) real(Eng),(TotalXS(ChS,Eng),ChS=1,nChS)
end do
close(200)

1000 format(F7.2,17(1x,ES9.3E2))


end program
