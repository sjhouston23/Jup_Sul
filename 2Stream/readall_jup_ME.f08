PROGRAM READALL

IMPLICIT NONE

INTEGER ns, nz, ne, nn
PARAMETER(ns=1544, nz=1544, ne=260, nn=10)
REAL*8  PHIUP(NZ,NE),PHIDWN(NZ,NE),ENERGY(NE),PRODH2(2,NZ),PRODH(1,NZ),PRODHE(1,NZ)
REAL*8  PRODCH4(6,NZ),z(nz),phiu(NE,nz),phid(NE,nz),ev(NE),A,chi,pii,RAD(NZ),S(NZ)
REAL*8  AGLW(NZ,NN, NN), AGLWH(NZ,NN, NN), totalescapeflux, dele(ne), totalelectrons
REAL*8  LymanH2(nz), WernerH2(nz), LymanAlphaH2(nz), LymanCascH2(nz), LymanH(nz)
REAL*8  TProdH2(nz), TProdH(nz), TProdHe(nz), TProdCH4(nz)
REAL*8  dummy1(nz), dummy2(nz), dummy3(nz), HExcite(nz), H2Excite(nz), TIonE(nz)
REAL*8  dummy11(nz), dummy21(nz), dummy31(nz)
real*8  CH4opac(nz), alt(nz), attenLyman(nz), attenWerner(nz)
INTEGER I, J, K, NSPEC, NEXC, IQQ, NSPECH, NEXCH
character*3 SZA
CHARACTER*100 flnm2, flnm3, flnm4, flnm5
DATA dele/20*0.5,70*1.0,10*2.0,20*5.0,10*10.0,20*10.0,10*50.0,10*100.0,40*200.0,&
        10*400,10*1000,10*2000,10*5000,10*10000.0/

dele=dele(ne:1:-1)
pii=4.0*atan(1.0)
open(50,FILE='input/elect/jelect.par',STATUS='UNKNOWN')
do i = 1,3
	read(50,*) !skip some lines
end do
!c Check which sza is being worked on to output to the correct folder
read(50,*) chi
if(chi .eq. 0.0)then
	SZA = '00/'
elseif(chi .eq. 60.0)then
	SZA = '60/'
elseif(chi .eq. 80.0)then
	SZA = '80/'
elseif(chi .eq. 90.0)then
	SZA = '90/'
elseif(chi .eq. 100.0)then
	SZA = 'ME/'
elseif(chi .eq. 200.0)then
	SZA = 'SE/'
end if
print*, 'reading file with chi= ', chi
read(50,*)
read(50,9999) flnm3 !output/SE/flux3d_JUP[Energy]keVJG.dat
read(50,9999) flnm2 !output/SE/airglow_JUP[Energy]keVJG.dat
do i = 1,2
	read(50,*) !skip some lines
end do
read(50,9999) flnm4 !../Output/[Energy]keV/H2_Excite_Prod_Comb.dat
write(*,*) flnm4
open(1, FILE='output/'//SZA//'phi_JG.chk', STATUS='OLD')
open(11,FILE='output/jelect/'//SZA//'PHIZ_2000keV_test.DAT',status='unknown')
!c      open(2,FILE='output/'//SZA//'airglow_JUP20keV.dat',STATUS='OLD')
open(2,FILE=flnm2,STATUS='OLD')
!c      open(2,FILE='common/output/'//SZA//'airglow1MeV.dat',STATUS='OLD')
open(23,FILE='common/output/SE/airglow5000.dat',status='unknown')
open(22,FILE='common/output/SE/secprd5000.dat',status='unknown')
open(24,file='common/output/SE/escapeflux5000.dat',status='unknown')
open(25,file='common/output/SE/escapecurrent5000.dat')
!c	  open(3,file='output/'//sza//'flux3d_JUP.dat',status='old')
open(3,file=flnm3,status='old')
open(33,file='output/jelect/'//sza//'ALTvsFLUX5000.dat')
!c	  open(4,file='OUTPUT/heatrate.chk',status='old')
!c	  open(44,file='my output/'//folder//'ElecHeatRate.dat')
open(4,file=flnm4,status='old') !H2 Excitation
open(1002,file='../Output/Airglow/SE/LymanWerner5000.dat',status='unknown')
open(1000,file='../Output/Airglow/SE/LWSpectrum5000.dat',status='unknown')
!open(1000,file='output/airglow/Electron_Beams/20keVbeam.dat',status='unknown')
!open(1001,file='CH4opacity.dat',status='unknown')
!do i=1,nz
!  read(1001,*) alt(i), CH4opac(i)
!end do

!C Read electron flux and then output flux up/dwn for altitudes selected
!do I = 1,3
!  read(1,*)
!end do
!do I = 1, nz
!  do J = 1, ne
!    read(1,*) ENERGY(J), PHIUP(I,J)
!  end do
!  read(1,*)
!  do J = 1, ne
!    read(1,*) ENERGY(J), PHIDWN(I,J)
!  end do
!  read(1,*)
!end do
!write(11,*)'Phi up'
!write(11,*) 'Energy 250km 300km 350km 400km 500km 750km 1000km 1500km 2000km 2500km 3000km'
!do J = 2,ne
!  write(11,*)ENERGY(J),PHIUP(5,J), PHIUP(10,J),PHIUP(15,J),PHIUP(20,J),PHIUP(30,J),&
!  PHIUP(55,J),PHIUP(80,J),PHIUP(130,J),PHIUP(180,J),PHIUP(230,J),PHIUP(280,J)
!end do
!write(11,*)
!write(11,*)'Phi down'
!write(11,*) 'Energy 250km 300km 350km 400km 500km 750km 1000km 1500km 2000km 2500km 3000km'
!
!do J = 2,ne
!  write(11,*)ENERGY(J),PHIDWN(5,J),PHIDWN(10,J),PHIDWN(15,J),PHIDWN(20,J),PHIDWN(30,J),&
!  PHIDWN(55,J),PHIDWN(80,J),PHIDWN(130,J),PHIDWN(180,J),PHIDWN(230,J),PHIDWN(280,J)
!end do
!C  production rate totals

!C first read in all the airglow emissions
do I = 1,4
  read(2,*)
end do
read(2,102) NSPEC, NEXC
print*, NSPEC, NEXC
do I =1,5
	read(2,*)
end do
do I = 1, NZ
	read(2,679)RAD(I), S(I), LymanH2(i), WernerH2(i), LymanAlphaH2(i), LymanCascH2(i)!(AGLW(I,NSPEC,IQQ),IQQ=1,NEXC)
!  write(*,*)rad(i)*1e-5, LymanH2(i), WenerH2(i), LymanAlphaH2(i), LymanCascH2(i)!(aglw(i,nspec,iqq),iqq=1,nexc)
end do
do I= 1,8
	read(2,*)
end do
read(2,102) NSPECH, NEXCH
print*, NSPEC, NEXC,NSPECH, NEXCH
do I =1,5
	read(2,*)
end do
!c     write(23,*)

do I = 1, NZ
	read(2,679)RAD(I),S(I),LymanH(i)!(AGLWH(I,NSPECH,IQQ),IQQ=1,NEXCH)
!	write(*,*)rad(i)*1e-5, LymanH(i)!(aglw(i,nspec,iqq),iqq=1,nexc),(aglwh(i,nspech,iqq),iqq=1,nexch)
end do
do i=1,10
	read(2,*)
end do
write(*,'(A,2ES10.3)') 'Total Production of Lyman from H2:', sum(LymanH2)*2e5*10**6, sum(LymanH2)*2e5*10**6

!C      NSPEC=1
!C      NEXC=4
!C      NSPECH=3
!C      NEXCH = 1
!C      Do i = 1,nz
!C       read(2,*) rad(i), (aglw(i,nspec,iqq),iqq=1,nexc)
!C       print*, i, rad(i),(aglw(i,nspec,iqq),iqq=1,nexc)
!C      enddo
!C      print*, 'space'
!C      read(2,*)
!C      do i = 1,nz
!C        read(2,*)  rad(i),(aglwh(i,nspech,iqq),iqq=1,nexch)
!C        print*, i, rad(i),(aglwh(i,nspech,iqq),iqq=1,nexch)
!C      enddo
!C	  do I = 1, NZ
!C      	read(2,*)RAD(I),S(I),(AGLWH(I,NSPECH,IQQ),IQQ=1,NEXCH)
!C      	write(23,*)rad(i)*1e-5, (aglw(i,nspec,iqq),iqq=1,nexc),
!C     + 	(aglwh(i,nspech,iqq),iqq=1,nexch)
!C      end do
!C
 679 format(5X,1P10E11.3)
!* Then read all the ion production rates:
do I = 1, NZ
  read(2,100) PRODH2(1,I), PRODH2(2,I)
!  write(*,*) i, prodh2(1,I)
end do
do I = 1,4
  read(2,*)
end do
do I = 1, NZ
  read(2,100) PRODHE(1,I), PRODH(1,I)
!  write(*,*) i, prodhe(1,I)
end do
do I = 1,4
  read(2,*)
end do
do I = 1, NZ
  read(2,101) (PRODCH4(J,I),J=1,6)
end do

write(22,*)'    Z','         H2+','        H+','      HE+','       CH4+','      CH3+',&
'     CH2+','       CH+','       C+'
do I = 1,NZ
write(22,103)RAD(I),PRODH2(1,I),PRODH2(2,I)+PRODH(1,I)+PRODCH4(4,I),PRODHE(1,I),&
PRODCH4(1,I),PRODCH4(2,I),PRODCH4(3,I),PRODCH4(5,I),PRODCH4(6,I)
TProdH2(i)=PRODH2(1,I)
TProdH(i)=PRODH2(2,I)+PRODH(1,I)+PRODCH4(4,I)
TProdHe(i)=PRODHE(1,I)
TProdCH4(i)=PRODCH4(1,I)
end do
do i=1,4
  read(4,*)
end do
H2Excite=0.0
do i=1,nz
  ! read(4,*) dummy1(i), H2Excite(i)
end do
H2Excite=H2Excite(nz:1:-1)
write(1000,*) '# Z[km]  Ly Dir.Ex. Ly Casc.   Ly - Ions  W Dir. Ex. ',&
              'W Ions     Ly AlphaH2 Ly Alpha H'
!TIonE=0.0d0
do i=1,nz !This is in units of photons/cm^3/s
  write(1000,104) RAD(I)/1e5,LymanH2(i),LymanCascH2(i),H2Excite(i)/2,WernerH2(i),H2Excite(i)/2,LymanAlphaH2(i),LymanH(i)
end do
do i=1,nz
  attenLyman(i)=(LymanH2(i)+LymanCascH2(i)+LymanAlphaH2(i)+LymanH(i)+(H2Excite(i)/2))*CH4opac(i)
  attenWerner(i)=(WernerH2(i)+(H2Excite(i)/2))*CH4opac(i)
end do
write(*,*) ''
write(*,'(2ES10.3)') sum(attenLyman)*2e5*10**6, sum(attenWerner)*2e5*10**6
write(1002,*) 'Alt [km]  Lyman-H2    Werner-H2   Ly Alpha H2 Ly Casc H2  ',&
'Ly Alpha H  H2 Excite   Total [photons/cm^-3/s]'
do i=1,nz
  write(1002,105) rad(i)/1e5,LymanH2(i),WernerH2(i),LymanAlphaH2(i),LymanCascH2(i),LymanH(i),H2Excite(i),&
  LymanH2(i)+WernerH2(i)+LymanAlphaH2(i)+LymanCascH2(i)+LymanH(i)+H2Excite(i)
end do
write(*,*) ''
write(*,*) 'Ion and airglow production: '
write(*,'(11ES10.3)') sum(TProdH2)*2e5*10**6, sum(TProdH)*2e5*10**6, sum(TProdHe)*2e5*10**6, &
sum(TProdCH4)*2e5*10**6, sum(LymanH2)*2e5*10**6, sum(LymanCascH2)*2e5*10**6, sum(H2Excite)*2e5*10**6/2, &
sum(WernerH2)*2e5*10**6, sum(H2Excite)*2e5*10**6/2, sum(LymanAlphaH2)*2e5*10**6, sum(LymanH)*2e5*10**6
write(*,*) ''
write(*,*) ''
write(*,*) 'Ion and airglow production: '
write(*,'(11ES10.3)') sum(TProdH2)*2e5, sum(TProdH)*2e5, sum(TProdHe)*2e5, &
sum(TProdCH4)*2e5, sum(LymanH2)*2e5, sum(LymanCascH2)*2e5, sum(H2Excite)*2e5/2, &
sum(WernerH2)*2e5, sum(H2Excite)*2e5/2, sum(LymanAlphaH2)*2e5, sum(LymanH)*2e5
write(*,*) ''
write(*,*) '******************************************************'
write(*,'(A,ES10.3)') 'Ion Lyman/Werner band emission:', sum(H2Excite)*2e5*10**6/2
write(*,*) '******************************************************'
write(*,*) ''

 100 format(39X,ES9.3,2X,ES9.3)
 101 format(39X,6(ES9.3,2X))
 102 format(61X,I4,I4)
 103 format(9ES10.3)
 104 format(F8.2,7(2x,ES10.3))
 105 format(F8.2,7(2x,ES10.3))
9999 format(1x,A50)
!C  altitude versus flux data
!* Read in data from flux3d.dat
!c	    read(3,* )ev(i)
do i=1,ne
  read(3,200)ev(i)
  read(3,*)
  do j=1,nz
    read(3,391)z(j),phiu(i,j),phid(i,j)
  enddo
enddo
!* Write altitude and energy labels
write(33,*)'#Alt',(ev(i),i=1,ne) !use # to prevent python from reading this data
do j=1,nz
  write(33,*)rad(j)*1e-5,(phiu(i,j),i=1,ne) !Write the altitude and flux up for all E at that altitude
enddo
write(33,*) !Leave a blank space between up and down
!* Write altitude and energy labels
write(33,*)'#Alt',(ev(i),i=1,ne) !use # to prevent python from reading this data
do j=1,nz
  write(33,*)'#',rad(j)*1e-5,(phid(i,j),i=1,ne) !Write the altitude and flux down for all E at that altitude
	!use # to prevent python from reading this data
enddo
write(24,*) '# Escape (Altitude=3000km) flux of electrons across the entire energy spectrum.'
totalescapeflux=0.0
!write(*,*) nz, rad(nz)
do i=1,ne
	j=nz
	write(24,*) ev(i), phiu(i,j), dele(i) !Hold J constant at the top of the atmosphere
	totalescapeflux=totalescapeflux+phiu(i,j)*dele(i)*ev(i)*pii*2
  totalelectrons=totalelectrons+phiu(i,j)*dele(i)*pii*2
end do
write(*,*) 'Total escape flux of electrons: ', totalescapeflux, 'eV/cm^2/s'
write(*,*) 'Total electrons: ', totalelectrons, 'electrons/cm^2/s'
write(*,*) 'Total current density: ', totalelectrons*1.60217662e-19*10**6, 'A/cm^2'
write(25,*) 'Escape flux [eV/cm^2/s]  Electrons [electrons/cm^2/s]  Current [A/cm^2]'
write(25,250) totalescapeflux, totalelectrons, totalelectrons*1.60217662e-19

close(1)
close(2)
close(3)
close(4)
close(11)
close(22)
close(23)
close(24)
close(25)
close(33)
close(50)
close(1000)
close(1002)

 200 format(61X,F10.2)
 250 format(8x,ES10.4,18x,ES10.4,13x,ES10.4)
 391 format(3e12.3)

!C  electron heating rate


print*,'read done'
stop
end
