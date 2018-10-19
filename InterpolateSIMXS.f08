program InterpolateSIMXS
!*******************************************************************************
!* Created by Stephen J. Houston 10.18.18
!*******************************************************************************
!* This program reads in SIM cross-sections for sulfur from .txt files and
!* interpolates them using both loglog linear and loglog spline interpolation.
!* All the cross-sections are from Schultz et al. (2019).
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer ChS,Eng,E !Charge state and energy
integer pProc,tProc !Projectile and target process
integer SS,DS,SPEX,DPEX !Projectile processes
integer SI,DI,TI,DA,SC,DC,TEX !Target processes

parameter(nProjProc=4) !Number of projectile processes
parameter(nTargProc=7) !Number of target processes
parameter(nEnergies=9) !Number of inital energies
parameter(nInterpEnergies=2000) !Number of interpolated energies
parameter(nChS=17) !Number of charge states from 0-16
parameter(SI=1,DI=2,TI=3,DA=4,SC=5,DC=6,TEX=7) !Target process numbers
parameter(SS=1,DS=2,SPEX=3,DPEX=4) !Projectile process numbers

real*8,dimension(nEnergies) :: Energy !Each initial energy
real*8,dimension(nChS,nEnergies,nTargProc,1+nProjProc) :: SIMxs !SIMxs points
real*8,dimension(nTargProc,1+nProjProc,nChS,nInterpEnergies) ::SIMxsInterp
!* SIMxs has an additonal projectile process which is no projectile process
!* i.e. SI, SI+SS, SI+DS, SI+SPEX, SI+DPEX (respectively)

character(len=100) :: filename
character(len=4),dimension(nTargProc) :: TargProcCap !Target processes
character(len=5),dimension(nProjProc+1) :: ProjProcCap !Projectile processes

!****************************** Data Declaration *******************************
data Energy/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0/
data TargProcCap/'SI  ','DI  ','TI  ','DCAI','SC  ','DC  ','TEX '/
data ProjProcCap/'     ','+SS  ','+DS  ','+SPEX','+DPEX'/
!**************************** Initialize Variables *****************************
SIMxs=0.0;SIMxsInterp=0.0
!******************************** Main Program *********************************
open(unit=100,file='./SIMXS/SIMXSall.dat',status='old') !All SIMxs data in one file
read(100,*) !Skip first line
read(100,1000) SIMxs !Read in all the data
close(100) !Close the file

do tProc=1,nTargProc !Loop through every target process
  do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
    do ChS=1,nChS !Loop through every charge state (0-16)
      Eng=1 !Set Eng variable back to 1
      do E=1,nInterpEnergies !Interpolation loop (1-2000 keV/u)
        if(E.ge.Energy(Eng+1)) Eng=Eng+1 !Go to next Energy when appropriate
        SIMxsInterp(tProc,pProc,ChS,E)=log(SIMxs(ChS,Eng,tProc,pProc))+&
        (log(real(E))-log(Energy(Eng)))*&
        (log(SIMxs(ChS,Eng+1,tProc,pProc))-log(SIMxs(ChS,Eng,tProc,pProc)))/&
        (log(Energy(Eng+1))-log(Energy(Eng)))
        ! write(*,*) E,exp(SIMxsInterp(tProc,pProc,ChS,E)),Energy(Eng)
      end do !End interpolation loop (1-2000 keV/u)
    end do !End loop through every charge state (0-16)
  end do !End loop through every projectile process plus 1
end do !End loop through every target process

do tProc=1,nTargProc !Loop through every target process
  do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
    write(filename,"('./SIMXSInterp/',A,A,'.dat')")&
    trim(TargProcCap(tProc)),trim(ProjProcCap(pProc))
    open(unit=101,file=trim(filename))
    write(101,1002) (i-1,i=1,17)
    do E=1,nInterpEnergies !Interpolation loop (1-2000 keV/u)
      write(101,1001) real(E),(exp(SIMxsInterp(tProc,pProc,ChS,E)),ChS=1,nChS)
    end do !End interpolation loop (1-2000 keV/u)
    close(101)
    write(filename,"('./SIMXS/',A,A,'p.dat')")&
    trim(TargProcCap(tProc)),trim(ProjProcCap(pProc))
    open(unit=102,file=trim(filename))
    write(102,1002) (i-1,i=1,17)
    do Eng=1,nEnergies
      write(102,1001) Energy(Eng),(SIMxs(ChS,Eng,tProc,pProc),ChS=1,nChS)
    end do
    close(102)
  end do !End loop through every projectile process plus 1
end do !End loop through every target process

1000 format(17(ES9.3E2,1x))
1001 format(F7.2,17(1x,ES9.3E2))
1002 format(' Energy',4x,9('S^',I0,'+',6x),8('S^',I0,'+',5x))

end program
