subroutine JEDIInterpolator(version,Iflux)
!*******************************************************************************
!* Created by Stephen J. Houston 9.5.18
!*******************************************************************************
!* This program interpolates JEDI ion spectra into a finer energy grid. It
!* outputs the intensity, energy bin width, and then a renormalized ion flux
!* which is too be input into the oxygen ion precipitation model. Alternatively,
!* the flux can instead be multiplied by a normalized (1 ion/cm^2/s) output from
!* JupOxyPrecip.
!*******************************************************************************

implicit real*8(a-h,o-z)

!*******************************************************************************
character(len=10),intent(in) :: version
real*8,dimension(37),intent(out) :: Iflux

parameter(nJebins=10,nSteps=3) !Number of JEDI energy bins
!nSteps is the number of energy data points you want in between the JEDI points

!* JEDI variables
real*8,dimension(nJebins) :: Jenergy,Jintensity,Jebins,Jflux
!*   Jenergy - JEDI energy bin [keV]
!*   Jintensity - JEDI ion flux [c/s/ster/cm^2/keV]
!*   Jebins - Size of JEDI energy bins [keV]
!*   Jflux - Jintensity values converted to [counts/cm^2/s]

character(len=10) date
character(len=12) time
character(len=100) filename
!* Interpolated variables
real*8,allocatable,dimension(:) ::Ienergy,Iintensity,Iebins

!****************************** Data Declaration *******************************
!* Width of JEDI energy bins: (May eventually need to be adjusted)
data Jebins/66.0,71.0,105.0,216.0,346.0,251.0,300.0,880.0,2280.0,5340.0/
!********************************* Initialize **********************************
!nInterp includes the JEDI energy data points as well
nInterp=nJebins+(nJebins-1)*nSteps
allocate(Ienergy(nInterp),Iintensity(nInterp),Iebins(nInterp))
pi=4.0d0*atan(1.0d0)
Jenergy=0.0;Jintensity=0.0;Jbins=0.0;Jflux=0.0
Ienergy=0.0;Iintensity=0.0;Ibins=0.0;Iflux=0.0
!*************************** Open JEDI Ion Spectrum ****************************
!write(version,'("v1")') !Filename of a JEDI spectrum (.d2s file)
write(filename,'("./JunoData/Spectra/",A,".d2s")') trim(version)
open(unit=100,file=trim(filename),status='old')
write(*,*)
do i=1,25 !Reading in the data measured by JEDI
  if(i.le.2.or.i.ge.4.and.i.le.15)read(100,*)
  if(i.eq.3)read(100,1001) date, time !Read the date and time of the flyby
  if(i.eq.3)write(*,1000) trim(version),date,time !Write to screen
  if(i.ge.16)read(100,1002) Jenergy(i-15),Jintensity(i-15)
end do
write(*,*)
write(*,1003)'Energy Bin:','JEDI Intensity:','Energy Bin Width:',&
             'Normalized Flux:' !Write out general information
do i=1,nJebins !Convert to [counts/cm^2/s]
!* The first 3 energy bins include both sulfur and oxygen. I'm assuming a 2:1
!* ratio of oxygen:sulfur (from SO_2)
  if(i.le.3)Jflux(i)=Jintensity(i)*2*pi*Jebins(i)*2/3
  if(i.ge.4)Jflux(i)=Jintensity(i)*2*pi*Jebins(i)
  write(*,1004)Jenergy(i),Jintensity(i),Jebins(i),Jflux(i)
end do
close(100) !Close JEDI measurement file
!**************************** Interpolate the data *****************************
j=0
do i=1,nInterp !Create new energy points that include the original points
  if(mod(i-1,nSteps+1).eq.0)j=j+1
  Ienergy(i)=log(Jenergy(j))+& !Want a log interpolation
            (log(Jenergy(j+1))-log(Jenergy(j)))*mod(i-1,nSteps+1)/(nSteps+1)
end do
Ienergy=exp(Ienergy) !Put back into a normal number
j=0
do i=1,nInterp !Loglog linearly interpolate the intensities
  if(mod(i-1,nSteps+1).eq.0)j=j+1
  Iintensity(i)=log(Jintensity(j))+(log(Ienergy(i))-log(Jenergy(j)))*&
                (log(Jintensity(j+1))-log(Jintensity(j)))/&
                (log(Jenergy(j+1))-log(Jenergy(j)))
  Iintensity(i)=exp(Iintensity(i)) !Put back into a normal number
end do
!************************* Find new energy bin widths **************************
lower_bound=145.0 !Lowest energy detectable by JEDI [keV]
upper_bound=10000.0 !Highest energy detectable by JEDI [keV]
do i=1,nInterp !Assume the new energy points are the midpoint of the bins
  if(i.eq.1)then !Need to account for lowest energy value
    Iebins(i)=(Ienergy(i+1)-Ienergy(i))/2.0+(Ienergy(i)-lower_bound)
  elseif(i.eq.nInterp)then !Need to account for highest energy value
    Iebins(i)=(upper_bound-Ienergy(i))+(Ienergy(i)-Ienergy(i-1))/2.0
  else
    Iebins(i)=(Ienergy(i+1)-Ienergy(i))/2.0+(Ienergy(i)-Ienergy(i-1))/2.0
  end if
end do
!************************ Calculate interpolated fluxes ************************
write(*,*)
write(*,*) 'Interpolated Values'
write(*,1003)'Energy Bin:','JEDI Intensity:','Energy Bin Width:',&
             'Normalized Flux:' !Write out general information
do i=1,nInterp !Convert to [counts/cm^2/s]
!* The first 3 JEDI energy bins include both sulfur and oxygen. I'm assuming
!* a 2:1 ratio of oxygen:sulfur (from SO_2). The energy 387 is what I
!* received in an email from Dennis Haggerty; however, Table 19 from the SSR
!* publication concering JEDI would suggest it's more like 322.30.
  if(Iebins(i).le.387)then
    Iflux(i)=Iintensity(i)*2*pi*Iebins(i)*2/3
  else
    Iflux(i)=Iintensity(i)*2*pi*Iebins(i)
  end if
  write(*,1004)Ienergy(i),Iintensity(i),Iebins(i),Iflux(i)
end do

1000 format('FILE: ',A5,'.d2s',/,'DATE: ',A10,/,'TIME: ',A12)
1001 format(39X,A10,1X,A12)
1002 format(5X,ES12.9,1X,ES13.10)
1003 format(3x,A11,2x,A15,2x,A17,2x,A16)
1004 format(F11.3,1x,F15.4,2x,F17.2,2x,F16.3)

end subroutine
