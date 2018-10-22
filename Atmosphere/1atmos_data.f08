subroutine atmos_data(species,neutdens,alt,arrsize,y2)
!*******************************************************************************
!* Created by Stephen J. Houston 2.20.18 from previous work
!*******************************************************************************
!************************************************************************
! This subroutine is designed to read the data files that contain the
!	neutral density for H_2, H, He and CH4 and return all the data
!************************************************************************

! Input: species; integer
!	  species 1
!     nH2 --> H2 density
!     Type: real*8 array
!     Units: cm^-3

!	  species 2
!     nHe --> He density
!     Type: real*8 array
!     Units: cm^-3

!	  species 3
!     nCH4 --> CH4 density
!     Type: real*8 array
!     Units: cm^-3

!	  species 4
!     nH --> H density
!     Type: real*8 array
!     Units: cm^-3

! Returns:
!	  neutdens --> neutral density data for the given species
!     alt --> altitude data for the corresponding densities
!     arrsize --> number of altitude/density points (size of alt, and neutdens)
!************************************************************************
!
implicit none
!
!*****************************DECLARATION OF VARIABLES*******************

integer i, iostatus, N !N is one plus the size of the array. The array is N-1
integer, intent(in) :: species
integer, intent(out) :: arrsize
real*8 x, y
real*8, allocatable, intent(out) :: neutdens(:), alt(:), y2(:)

! open the neutral density file for the given species
if (species .eq. 1) then
  open(unit=species+25, file='./Atmosphere/Input/H2_M.txt',&
  status='old') !H2	!This has the data used in Ozak et al. model which comes from Maurellis et al. papers
!  open(unit=species+25, file='H2_JG.txt', status='old') !This data comes from Denis Grodent

elseif (species .eq. 2) then
  open(unit=species+25, file='./Atmosphere/Input/HE_M.txt',&
  status='old') !Helium density
!  open(unit=species+25, file='HE_JG.txt', status='old')

elseif (species .eq. 3) then
	 open(unit=species+25, file='./Atmosphere/Input/CH4_M.txt',&
   status='old') !CH4 density
!  open(unit=species+25, file='CH4_JG.txt', status='old')

elseif (species .eq. 4) then
  open(unit=species+25, file='./Atmosphere/Input/H_M.txt',&
  status='old') !H density
!  open(unit=species+25, file='H_JG.txt', status='old')

endif
! determine the size of the arrays
N=0
iostatus = 0
do while (iostatus.eq.0)
  read (species+25,*,iostat=iostatus) x,y
  N=N+1
end do
arrsize=N-1

! allocate the neutral density and altitude arrays to the correct size
allocate (neutdens(arrsize),alt(arrsize),y2(arrsize))
neutdens = 0.0
alt = 0.0
y2 = 0.0
rewind (species+25)

! read in the neutral density and altitude data
do i=1,arrsize
  read (species+25,*) neutdens(i), alt(i)
end do
call spline(alt, neutdens, arrsize, y2)
!write(*,*) y2
close(species+25)
return
stop
end
