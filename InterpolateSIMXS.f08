program InterpolateSIMXS
!*******************************************************************************
!* Created by Stephen J. Houston 10.18.18
!*******************************************************************************
!* This program reads in SIM cross-sections for sulfur from .txt files and
!* interpolates them using both loglog linear and loglog spline interpolation.
!* All the cross-sections are from Schultz et al. (2019).
!* The initial cross-sections are stored in SIMXSall.dat, which is generated
!* with SIMXSCalculation.f08 in the SulfurXS directory.
!*******************************************************************************

implicit none!real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer i,j,k,l,m,n
integer nProjProc,nTargProc,nEnergies,nInterpEnergies,nChS,ndum

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
real*8,dimension(nEnergies) :: SIMxs_tmp,SIMxs_tmp2
real*8 SIMxs_tmpI,SpE !Spline variables
real*8,dimension(nEnergies,nChS) :: SigTot
real*8,dimension(nInterpEnergies,nChS) :: SigTotInterp
real*8,dimension(nChS,nEnergies,nTargProc,1+nProjProc) :: SIMxs !SIMxs points
real*8,dimension(nTargProc,1+nProjProc,nChS,nInterpEnergies) :: SIMxsInterp
!* SIMxs has an additonal projectile process which is no projectile process
!* i.e. SI, SI+SS, SI+DS, SI+SPEX, SI+DPEX (respectively)

character(len=100) :: filename
character(len=4),dimension(nTargProc) :: TargProcCap !Target processes
character(len=4),dimension(nTargProc) :: TargProcCap2 !Target processes
character(len=5),dimension(nProjProc+1) :: ProjProcCap !Projectile processes

!****************************** Data Declaration *******************************
data Energy/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0/
data TargProcCap/'SI  ','DI  ','TI  ','DCAI','SC  ','DC  ','TEX '/
data TargProcCap2/'  SI','  DI','  TI','DCAI','  SC','  DC',' TEX'/
data ProjProcCap/'     ','+SS  ','+DS  ','+SPEX','+DPEX'/
!**************************** Initialize Variables *****************************
SIMxs=0.0;SIMxsInterp=0.0;SigTot=0.0;SigTotInterp=0.0;SIMxs_tmp=0.0
SIMxs_tmp2=0.0;SIMxs_tmpI=0.0
!******************************** Main Program *********************************
! open(unit=100,file='./SIMXS/SIMXSall.dat',status='old')!All SIMxs data in a file
! read(100,*) !Skip first line
! read(100,1000) SIMxs !Read in all the data
! close(100) !Close the file

open(unit=100,file='./SIMXS/SIMXSfinal1.txt',status='old')
do tProc=1,nTargProc
  do pProc=1,nProjProc+1
    do i=1,2
      read(100,*)
    end do
    do Eng=1,nEnergies
      read(100,*) ndum,(SIMxs(ChS,Eng,tProc,pProc),ChS=1,nChS)
    end do
  end do
end do
close(100)

open(unit=204,file='./SIMXS_TotalOG.dat')
write(204,*) 'Energy     S        S^+       S^++      S^3+      S^4+      &
&S^5+      S^6+      S^7+      S^8+      S^9+     S^10+     S^11+     S^12+&
&     S^13+     S^14+     S^15+     S^16+'
do i=1,nEnergies
  write(204,20400) int(Energy(i)),(sum(SIMxs(ChS,i,:,:)),ChS=1,nChS)
  do ChS=1,nChS
    SigTot(i,ChS)=sum(SIMxs(ChS,i,:,:))
  end do
end do
20400 format(I7,17(2x,ES8.2E2))
close(204)
open(unit=206,file='./XS.dat')
do tProc=1,nTargProc !Loop through every target process
  do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
    do ChS=1,nChS !Loop through every charge state (0-16)
      Eng=1 !Set Eng variable back to 1
      SpE=0.0
      do E=1,nInterpEnergies !Interpolation loop (1-2000 keV/u)
        SpE=SpE+1.0
        if(E.ge.Energy(Eng+1)) Eng=Eng+1 !Go to next Energy when appropriate
        if(E.eq.2000) Eng=8 !Don't want Eng to go out of bounds
        SIMxsInterp(tProc,pProc,ChS,E)=log(SIMxs(ChS,Eng,tProc,pProc))+&
        (log(real(E))-log(Energy(Eng)))*&
        (log(SIMxs(ChS,Eng+1,tProc,pProc))-log(SIMxs(ChS,Eng,tProc,pProc)))/&
        (log(Energy(Eng+1))-log(Energy(Eng)))

        SIMxsInterp(tProc,pProc,ChS,E)=exp(SIMxsInterp(tProc,pProc,ChS,E))

        if(E.ge.10.and.E.le.1000)then
          if(E.eq.10)then
            do i=1,nEnergies
              SIMxs_tmp(i)=SIMxs(ChS,i,tProc,pProc) !Create a vector for spline
            end do
          end if
          call spline(log(Energy),log(SIMxs_tmp),nEnergies,SIMxs_tmp2)
          call splineinterp(log(SpE),log(Energy),log(SIMxs_tmp),nEnergies,&
                            SIMxs_tmp2,SIMxs_tmpI)
          SIMxsInterp(tProc,pProc,ChS,E)=exp(SIMxs_tmpI)
          ! write(*,*) E,SpE,Energy(Eng),SIMxs_tmp(Eng),nEnergies,SIMxs_tmp2(Eng),exp(SIMxs_tmpI),SIMxs_tmpI
        end if
        ! write(*,*) E,SIMxsInterp(tProc,pProc,ChS,E)
        ! if(tProc.eq.1.and.pProc.eq.1)then
        !   SigTotInterp(E,ChS)=log(SigTot(Eng,ChS))+&
        !   (log(real(E))-log(Energy(Eng)))*&
        !   (log(SigTot(Eng+1,ChS))-log(SigTot(Eng,ChS)))/&
        !   (log(Energy(Eng+1))-log(Energy(Eng)))
        !   SigTotInterp(E,ChS)=exp(SigTotInterp(E,ChS))
        !   ! write(*,*) E,SigTotInterp(E,ChS)!,Energy(Eng)
        ! end if
      end do !End interpolation loop (1-2000 keV/u)
    end do !End loop through every charge state (0-16)
  end do !End loop through every projectile process plus 1
end do !End loop through every target process
! SIMxsInterp=exp(SIMxsInterp) !Revert back from log
! do i=1,nInterpEnergies
!   write(206,20400) i,(sum(SIMxsInterp(:,:,ChS,i)),ChS=10,10)
! end do

do tProc=1,nTargProc !Loop through every target process
  do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
    write(filename,"('./SIMXSInterp/',A,A,'.dat')")&
    trim(TargProcCap(tProc)),trim(ProjProcCap(pProc))
    open(unit=101,file=trim(filename))
    write(101,1002) (i-1,i=1,17)
    do E=1,nInterpEnergies !Energy interpolation loop (1-2000 keV/u)
      do ChS=1,nChS !Loop through every charge state
        if(isnan(SIMxsInterp(tProc,pProc,ChS,E)))& !Get rid of any NaN values
        SIMxsInterp(tProc,pProc,ChS,E)=0.0
      end do !End loop through every charge state
      write(101,1001) real(E),(SIMxsInterp(tProc,pProc,ChS,E),ChS=1,nChS)
    end do !End energy interpolation loop (1-2000 keV/u)
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

open(unit=103,file='./SIMXSInterp/SIMXSInterpAll.txt')
do tProc=1,nTargProc !Loop through every target process
  do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
    write(103,1111) TargProcCap2(tProc),ProjProcCap(pProc)
    write(103,*) 'Energy     S        S^+       S^++      S^3+      S^4+      &
    &S^5+      S^6+      S^7+      S^8+      S^9+     S^10+     S^11+     S^12+&
    &     S^13+     S^14+     S^15+     S^16+'
    do E=1,nInterpEnergies !Energy interpolation loop (1-2000 keV/u)
      do ChS=1,nChS !Loop through every charge state
        if(isnan(SIMxsInterp(tProc,pProc,ChS,E)))& !Get rid of any NaN values
        SIMxsInterp(tProc,pProc,ChS,E)=0.0
      end do !End loop through every charge state
      write(103,1001) real(E),&
        (SIMxsInterp(tProc,pProc,ChS,E),ChS=1,nChS)
    end do !End energy interpolation loop (1-2000 keV/u)
  end do !End loop through every projectile process plus 1
end do !End loop through every target process
close(103)
! open(unit=103,file='./SIMXSInterp/SIMXSInterpAll.dat')
! write(103,1100) exp(SIMxsInterp) !Write out every cross-section to a single file
! close(103)
do eng=1,nInterpEnergies
  do ChS=1,nChS
    SigTotInterp(Eng,ChS)=sum(SIMxsInterp(:,:,ChS,Eng))
  end do
end do
open(unit=205,file='./SIMXSInterp_TotalOG.dat')
write(205,*) 'Energy     S        S^+       S^++      S^3+      S^4+      &
&S^5+      S^6+      S^7+      S^8+      S^9+     S^10+     S^11+     S^12+&
&     S^13+     S^14+     S^15+     S^16+'
do i=1,nInterpEnergies
  write(205,20400) i,(SigTotInterp(i,ChS),ChS=1,nChS)
end do
close(205)

1000 format(17(ES9.3E2,1x))
1001 format(F7.2,17(1x,ES9.3E2))
1002 format(' Energy',4x,9('S^',I0,'+',6x),8('S^',I0,'+',5x))
1100 format(20(ES9.3E2,1x))
1111 format(A4,A5,'----------')
1112 format(I7,17(1x,ES9.3E2))

end program
