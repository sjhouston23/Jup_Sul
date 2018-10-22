subroutine coldens2(alt,nH2,nHe,nCH4,nH,ColDensH2,ColDensHe,ColDensCH4,&
  ColDensH,TotColDens)
!*******************************************************************************
!* Created by Stephen J. Houston 2.20.18 from previous work
!*******************************************************************************
!************************************************************************
!     This soubroutine is designed to calculate the column denstity in
!     the Jovian atmosphere for H2,He,H and CH4 as a function of height
!     alt. First you must call the atmosphere subroutine to get the neutral densities
!     and then pass those values into this subroutine
!************************************************************************

! Input:
!     alt--> Height
!       Type: real*8
!       Units: km

!	  species 1
!     nH2--> H2 density
!     Type: real*8
!     Units: cm^-3

!	  species 2
!     nHe--> He density
!     Type: real*8
!     Units: cm^-3

!	  species 3
!     nCH4--> CH4 density
!     Type: real*8
!     Units: cm^-3

!	  species 4
!     nH--> H density
!     Type: real*8
!     Units: cm^-3

! Returns:
!		ColDensH2--> H2 Column density
!		Type: real*8
!		Units: cm^-2

!		ColDensHE--> HE Column density
!		Type: real*8
!		Units: cm^-2

!		ColDensCH4--> CH4 Column density
!		Type: real*8
!		Units: cm^-2

!   ColDensH--> H Column density
! 	Type: real*8
!   Units: cm^-2

!*************************DECLARATION OF VARIABLES***********************
implicit none

real*8 DZ,H2,HE,CH4
real*8 alt, ColDensH2, ColDensHe, ColDensCH4, ColDensH,TotColDens, H
integer steps, atmoslen
parameter (atmoslen = 1544)
real*8 nH2(atmoslen),nHe(atmoslen),nCH4(atmoslen),nH(atmoslen)


!************************************************************************
! Initialize sums to zero
DZ = 2.0 ! make sure this interval matches the dz (altitude bin size) from the neutral density arrays
steps=int((3000.0-alt)/dz)
H=DZ*1e5
!write(*,*) alt,steps,dz
!if(steps.eq.0)then
!  ColDensH2=nH2(1)*H
!  ColDensHe=nHe(1)*H
!  ColDensCH4=nCH4(1)*H
!  ColDensH=nH(1)*H
!  goto 100
!end if
!write(*,*) '1'
if(mod(steps,2).eq.0)call simp(steps,H,nH2,ColDensH2) !for H2
if(mod(steps,2).eq.0)call simp(steps,H,nHe,ColDensHe) !for He
if(mod(steps,2).eq.0)call simp(steps,H,nCH4,ColDensCH4)!for CH4
if(mod(steps,2).eq.0)call simp(steps,H,nH,ColDensH) !for H
if(mod(steps,2).eq.1)call simp2(steps,H,nH2,ColDensH2) !for H2
if(mod(steps,2).eq.1)call simp2(steps,H,nHe,ColDensHe) !for He
if(mod(steps,2).eq.1)call simp2(steps,H,nCH4,ColDensCH4)!for CH4
if(mod(steps,2).eq.1)call simp2(steps,H,nH,ColDensH) !for H
!100 continue
!write(*,*) '2'
TotColDens = ColDensH2 + ColDensHE + ColDensCH4 + ColDensH !total column density
return
end


!**********************************************************************
subroutine simp (N,H,F,INTEG)
!     integration via simpson's rule
!     N: number of intervals for the integration
!     H: integration interval dx
!     F: vector with y values to be integrated
!     INTEG: resultant integral
!     This subroutine uses Simpson's rule to integrate f(x)

implicit none

integer N
real*8 F(N), F0, F1(N/2-1), F2(N/2), H, INTEG, SF1,SF2
integer K
F0 = 0.0
F1 = 0.0 !SUM OF ODD F(I)
F2 = 0.0 !SUM OF EVEN F(I)

F1=F(3:N:2)
F2=F(2:N:2)
SF1=sum(F1) !Summation of all odd elements
SF2=sum(F2) !Summation of all even elements
F0 = F(1) + F(N)

INTEG = H * (F0 + 2*SF2 + 4*SF1)/3

return
end

!**********************************************************************
subroutine simp2 (N,H,F,INTEG)
!     integration via simpson's rule
!     N: number of intervals for the integration
!     H: integration interval dx (can't be variable)
!     F: vector with y values to be integrated
!     INTEG: resultant integral
!     This subroutine uses Simpson's rule to integrate f(x)

implicit none

integer N
real*8 F(N), F0, F1(N/2), F2(N/2), H, INTEG, SF1,SF2
F0 = 0.0
F1 = 0.0 !SUM OF ODD F(I)
F2 = 0.0 !SUM OF EVEN F(I)
!!!write(*,*) F
F1=F(3:N:2)
F2=F(2:N:2)
SF1=sum(F1) !Summation of all odd elements
SF2=sum(F2) !Summation of all even elements
F0 = F(1) + F(N)

INTEG = H * (F0 + 2*SF2 + 4*SF1)/3

return
end
