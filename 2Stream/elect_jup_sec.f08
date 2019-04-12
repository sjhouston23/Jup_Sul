PROGRAM JELECT

!C     ===============  JELECT  ========================================
!C     This program is the electron code used originally for
!C     Titan and comets, and now it is being implemented for Jupiter.
!C     Nataly Ozak, January 10, 2011.
!C 	   This is a new variable of the code made to read in secondary electron
!C     prod rates calculated by the ion precipitation MC code
!C     There is an electron production that goes upward and another one downward
!C     which are not the same -> the production is not isotropic unlike the photoelectrons
!C     April 5, 2012. We will ingnore the photoelectrons for now.
!******************
!C   This code was made suitable for comets
!C   09/29/04   Ina P. Robertson

!C   this version ( >;46) reads in an electron temperature profile and then
!C   extrapolates it onto a specific field line of the calculation.
!C 			Mar 1, 1991

!C  telect.for;32 is the first working version of titan's two-stream code.
!C  This version has the following features:
!C   1. it is not for optically-thin only.
!C   2. it has the 3-d output as featured by velect.for
!C   3. it has the airglow calculations.
!C   4. it has a spherical symmetric atmosphere. The
!C      neutral densities follows Ned Keller's model. The neutrals considered
!C      are N2 and CH4.
!C      For the Jupiter case, the atmosphere is adapted
!C	   from Maurellis et al. 2000. It has H2, CH4, H and He. Use files
!C      H2_M.txt, He_M.txt, CH4_M.txt, H_M.txt, ATMDATA_M.f, ATM_M.f
!C   5. Ne = Ni is taken from N.Keller's photo- chemical equilibrium model.
!C      Te is taken to be 1000 K, and to be changed later when the
!C      electron energy code is running.
!C   6. 1st assumption will take magnetospheric electrons as a maxwellian
!C      with  Nme = 0.1 cm-3, and Tme = 500 eV.
!C   8. This version assumes a spherically summetric photoelectron production
!C      rate. The vertical profile taken is from ephot output. A 60sza run is
!C      used.  The improved
!C      model might have simple fit photoelectron production rate of several
!C      sza runs, and to incorporate the 20 deg. angle between the sun and the
!C      corotating magnetic field/plasma of Saturn.
!C                                 Sept. 5, 1990
!C   9. The program is split to :  telect.for + telect.sub
!C      after version 6.			Sept. 11, 1990
!C  10. This version assumes the center of titan is located at the center
!C 	of the tangential circle that intersects the parabola at the sub-
!C	solar pt. For this geometry, the radius at the flanks is twice the
!C	distance of the subsolar radius.
!C        The atmosphere parameter outside of r=5100km is from exponential
!C      extrapolation of the lower atmosphere.
!C					Sept. 19, 1990.
!C   11. e impact ion  prod rates computed in sub eimpir. Oct.3,1990
!C   12. Modified to allow parameter file to update solar flux file.
!C                                       J. Clark  Feb 24, 2004
!C *********************************************************************!
integer IonEnergy
PARAMETER(ns=1544, nz=1544, ne=260, nn=10) !Make sure that these values are updated if there is a change in the alt. bins.
!c Make sure to change in EIMPIR sub. and Impit sub.
!c 	  PARAMETER(ns=1544, nz=1544, ne=260, nn=10)
!c	ns = number of bins along the magnetic field line
!c	nz = number of vertical altitude bins
!c 	ne = number of energy bins
!c	nn = number of max neutral species?

real EE(ne),DELE(ne),PHIINF(ne),PE(nn,ne),PI(nn,ne),SIGA(nn,ne),SIGS(nn,ne),SEC(nn,ne)
real PHIUP(ne,ns),PHIDWN(ne,ns),PRODUP(ns,ne),PRODWN(ns,ne),UFX(INT(5E3),ne)
! arrays of variable E(nergy)
real S(ns),RAD(ns),AREA(ns),TE(ns),EHEAT(ns),ALPHA(ns),BETA(ns),GAMMA(ns),T1(ns),T2(ns),TPROD(ns)
real TSIGNE(ns),SECION(ns),PHINET(ns),SUMUP(ns),SUMDWN(ns),SUMNET(ns),SUMDIF(ns),ZSPEC(nn,ns)
real SION(nn,ns),PHI(ns),DUMMY(ns),AGLW(ns,4,10),DENSTH(ns),EBARTH(ns),SNE(ns),H2,HE,H,CH4,R
real SIGHFAST,tprodf(ns),tprodb(ns), sprodf(nz,ne), sprodb(nz,ne),dummy1
!........ arrays of variable S (along a B line)
real rkm(nz),ZNE(nz),PROD(nz),PEPROD(nz,ne),PEPROD1(nz),DUMM(nz),dneu(nz,nn),dneutemp(nz,nn)
real xnn(nz),ter(nz),rte(nz),eprodf(nz,ne),eprodb(nz,ne),seprod1(nz),seprod2(nz),prodf(nz)
real prodb(nz), eprodfR2(nz,ne), eprodbR2(nz,ne),eprodbR3(nz,ne), eprodfR3(nz,ne)
!........ arrays of variable Z (alt. from ephot)
real TION(nn),TSA(nn),SECP(nn),ALOSS(nn)
real Q(50),PHHD(50), NES, delz, CHI
REAL del(ne)
dimension NMM(20),DDDE(20),ID(nn)
INTEGER*2 SPEC(nn)
INTEGER PARA,choice,dummy2,COUNT, j, V
INTEGER ENERGY, jbins, nbins, jen
character*9 aglwlab(nn,20)
COMMON/TWO/ALPHA,BETA,GAMMA,IMAX,DELS,FAC,PHI,PHIINF,J
COMMON/RABBIT/NAPP(4),CRAB(60,10,4),ECRAB(60,10,4),NCRAB(10,4)
common/atmos/s,rad,ramagl(ns),ee,dele,zspec, delz  	! for sub. eimpir
common /angle/theta,rp,ramdir			! for sub. para2
CHARACTER*100 FLNM4,FLNM11,FLNM1,FLNM15,FlNM21,FLNM28,flnm14
CHARACTER*100 FLNM19,flnm20,flnm22,flnm2323
CHARACTER*50 tag1,tag2
CHARACTER*2 SZA
data pii/3.1415926/
DATA DEL/20*0.5,70*1.0,10*2.0,20*5.0,10*10.0,20*10.0,10*50.0,10*100.0,40*200.0,&
        10*400,10*1000,10*2000,10*5000,10*10000.0/

!C     COMET/NEW VARIABLES
REAL DELTA, TAU, QZ, DHEL, NEUTDENS, T, FIONZ, SKC(INT(1E5))
REAL RKC(INT(1E5)),XKC(INT(1E5)), YKC(INT(1E5)),Z, alt
INTEGER COUNT2

EXTERNAL ATM, ATMDATA
COMMON/KC/SKC, RKC, XKC, YKC

DHEL = 5.2 !Heliospheric distance to the object in AU
QZ = 5.0E27 !USED IN THE NEUTD SUBROUTINE, WHICH IS NOW DIFFERENT
QZ = QZ/DHEL/DHEL !NOT USED.
!C     ==================================================
!C                     SECTION 1: Initializaiton
!C                     -------------------------

!C  *** MODIFIED BY JC 8/20 TO give an intro and make a little
!C      user friendly
PRINT*, 'PROGRAM = JELECT'
PRINT*, 'CALCULATES THE SECONDARY PRODUCTION RATES OF'
PRINT*, 'THE MAJOR IONS'

9999    FORMAT(1x,A50)

print *
print *,' Initializing...'

OPEN(1,FILE='input/elect/jelect.par',STATUS='OLD')
!c        OPEN(1,FILE='celect.par',STATUS='OLD')
READ(1,*),nprodf  !include solar electron/secondary electron flux from photoe code/ion precip model

if(nprodf .eq. 1)then
!c        	print*, 'photoelectron input'
	print*, 'secondary electron input'
elseif(nprodf .eq. 0)then
	print*, 'no photoelectrons or sec electrons included'
endif

READ(1,*),NES  !
print*, 'nes = ', nes
READ(1,*),PARA !Parabolic (1) or radial (0) field line

if(para .eq. 1)then
	print*, 'PARABOLIC field line'
elseif(para .eq. 0)then
	print*, 'RADIAL field line'
endif
! 111 format(I4)
read(1,*) IonEnergy
write(*,*) '--------------------------------------------------------------------------------'
write(*,*) "Ion energy of:",IonEnergy,"keV/u"
write(*,*) '--------------------------------------------------------------------------------'
READ(1,*) CHI
print*, 'SZA =', chi

if(chi .eq. 0)then
	SZA = '00'
elseif(chi .eq. 60)then
	SZA = '60'
elseif(chi .eq. 80)then
	SZA = '80'
elseif(chi .eq. 90)then
	SZA = '90'
elseif(chi .eq. 100)then
	SZA = 'ME'
elseif(chi .eq. 200)then
	SZA = 'SE'
ENDIF

READ(1,9999),flnm11 !output/SE/pespect_j2MeVJG.dat
print*, flnm11, '11'
write(flnm28,"('output/',A,'/flux3d_JUP',I0,'keVJG-eq.dat')") SZA,IonEnergy
READ(1,9999)!,flnm28 !output/SE/flux3d_JUP1000keVJG.dat
print*, flnm28, '28'
write(flnm15,"('output/',A,'/airglow_JUP',I0,'keVJG-eq.dat')") SZA,IonEnergy
READ(1,9999)!,flnm15 !output/SE/airglow_JUP1000keVJG.dat
print*, flnm15, '15'
write(flnm20,"('../Output/',I0,'/2Str_Elect_Fwd-Comb.dat')") IonEnergy
read(1,9999)!, flnm20 !common/output/SE/prdelf2str1000keV.dat
print*, flnm20, '20'
write(flnm22,"('../Output/',I0,'/2Str_Elect_Bwd-Comb.dat')") IonEnergy
read(1,9999)!, flnm22 !common/output/SE/prdelb2str1000keV.dat
print*, flnm22, '22'
read(1,*) !Skip H2* Prod
write(flnm2323,"('common/output/',A,'/ForwardElecProd',I0,'-eq.dat')") SZA,IonEnergy
read(1,9999)!, flnm2323 !Forward electron production output file
print*, flnm2323, '2323'
CLOSE(1)
!* Check the sza for the file labeling
!c	OPEN(UNIT=33, FILE='PARAB32.DAT', STATUS='OLD') !This file is read if there is a parabola for B

!c	DO I = 1,2
!c	   READ(33,*)
!c	ENDDO

!c	DO I = 1, 1000
!c	   READ(33,*)XKC(I), YKC(I), RKC(I), SKC(I)
!c	ENDDO


!C
!CCCC  $1.1:  open input files.
!C  ------------------------
OPEN(3,FILE='input/elect/geopar.t',STATUS='OLD')	!geometric parameters
OPEN(7,FILE='common/input/XGRIDNEW_EXTEND.jup',STATUS='OLD')	!extended energy grid
!C		OPEN(7,FILE='common/input/xgridnew_ext.jupiter',STATUS='OLD')	!extended energy grid
!C	  OPEN(9,FILE='common/output/xsect_nobcksc.jupiter',STATUS='OLD')	!XSECT output
OPEN(9,FILE='common/output/xsect.jupiter',STATUS='OLD')	!XSECT output
! OPEN(11,FILE=FLNM11,STATUS='OLD')	!p.e.prod.rates.
write(*,*) 'NO PHOTOELECTRONS CONSIDERED'
OPEN(13,FILE='input/elect/eaglw.jup',STATUS='OLD')	!airglow cs
OPEN(17,FILE='common/input/Tempjup_JG.txt', STATUS='OLD') !Neutral Temperature profiles from Maurellis
!	  OPEN(19, FILE='common/output/'//SZA//'jpce2.out',
!     +  STATUS='OLD') !Chemical code output for electron densities !!!!!Commented out to not use Chemical Code Input
!***1MeV/u
!c      open(20,file='common/output/'//SZA//'atmosJG/prodelf2str1MeVJG.dat
!c     +' ,status='old')
!c      open(21, file='common/output/'//SZA//'prodelf2str1MeV2.dat',
!c     + status='old')
!c      open(22,file='common/output/'//SZA//'atmosJG/prodelb2str1MeVJG.dat
!c     +',status='old')
!c      open(23, file='common/output/'//SZA//'prodelb2str1MeV2.dat',
!c     + status='old')
!c      open(24, file='common/output/'//SZA//'prodelf2str1MeV3.dat',
!c     + status='old')
!c      open(25, file='common/output/'//SZA//'prodelb2str1MeV3.dat',
!c     + status='old')

!***2MeV/u
!c	  open(20, file='common/output/'//SZA//'prodelf2str2MeV1JG.dat',
!c     + status='old')
!c      open(21, file='common/output/'//SZA//'prodelf2str2MeV2.dat',
!c     + status='old')
!c      open(22, file='common/output/'//SZA//'prodelbk2str2MeV1JG.dat',
!c     + status='old')
!c      open(23, file='common/output/'//SZA//'prodelbk2str2MeV2.dat',
!c     + status='old')
!***1.5MeV/u
open(20,FILE=flnm20,STATUS='OLD')
!c   	  open(20, file='common/output/'//SZA//'atmosJG/prdelf2str15MeVJG.da
!c     +dat', status='old')
!c      open(21, file='common/output/'//SZA//'prdelf2str15MeV2.dat',
!c     + status='old')
open(22,FILE=flnm22,STATUS='OLD')
!c      open(22, file='common/output/'//SZA//'atmosJG/prdelb2str15MeVJG.dat
!c     +t', status='old')
!c      open(23, file='common/output/'//SZA//'prdelb2str15MeV2.dat',
!c     + status='old')
open(2323,file=flnm2323,status='unknown')
!C
!CCCC  $1.2: open output files.
!C   ------------------------
OPEN(34,FILE='output/'//SZA//'/heatrate_JG.chk',STATUS='UNKNOWN')
OPEN(8,FILE='output/'//SZA//'/FLUX_JG.chk',STATUS='UNKNOWN')
OPEN(12,FILE='output/'//SZA//'/ATMOS_JG.chk',STATUS='UNKNOWN')
OPEN(14,FILE='output/'//SZA//'/phi_JG.chk',STATUS='UNKNOWN')
!OPEN(15,FILE='output/SE/airglow_JUP1000keVJG.dat',STATUS='UNKNOWN')
OPEN(15,FILE=FLNM15,STATUS='UNKNOWN')  !output/airglow_JUP.dat
OPEN(16,FILE='output/'//SZA//'/CHECK_JG.chk',STATUS='UNKNOWN')
OPEN(18,FILE='output/'//SZA//'/SUTH_JG.chk', STATUS='UNKNOWN')
OPEN(28,FILE=FLNM28,STATUS='UNKNOWN') !output/flux3d_JUP.dat
!OPEN(28,FILE='output/SE/flux3d_JUP1000keVJG.dat',STATUS='UNKNOWN')
!c	  OPEN(50,FILE='R_ALPHA.DAT',STATUS='UNKNOWN')
OPEN(51,FILE='output/jelect/CELECT_DENSJG.DAT',STATUS='UNKNOWN')
OPEN(52,FILE='output/jelect/EIMPIRION_2JG.DAT',STATUS='UNKNOWN')
OPEN(53,FILE='output/jelect/DISSOC_JG.DAT', STATUS='UNKNOWN')
OPEN(54,FILE='output/jelect/'//SZA//'/TotIONprod1_JG.dat',STATUS='unknown')!, access ='append')
!C      ====================================================
!C                    SECTION 2: input section
!C		     ------------------------
!C	print *,' Start input section...'
!* From geopar.t file:
      READ(3,518) TES,IGEO !Secondary electron energy, ?
518   FORMAT(1PE10.3,I5)
      READ(3,301)IMAX
!C                       .......	IMAX is the number of s increments
      READ(3,301) IPMAX
!C     			... IPMAX is the number of alt. increments in EPHOTC.
      READ(3,303) HSTEP
      READ(3,303) HMIN
!C     			... HSTEP,HMIN are altitude step and number in EPHOTC.
      READ(3,301)JMAX,JMIN
      print*, 'jmax in elect=', jmax
!C     JMAX IS THE NUMBER OF ENERGY INCREMENTS
      READ(3,317) DIPANG,AVMU,TINF
      READ(3,303)DELS
      delz = DELS		! *** dels, z <=> dels, s
!C     DELS IS THE arc length s INCREMENT IN CM
      READ(3,303) S(1)

!C     read in the altitude of the subsolar point HSUBS in cm.
      READ(3,303)HSUBS

      READ(3,*) dummy2,JLOCAL !NPRODF -> should match the value on jelect.par
      IF(dummy2 .ne. nprodf) THEN
PRINT*, 'JELECT.PAR AND GEOPAR.T HAVE DIFFERENT PHOTOELECRTRON INPUT CRITERIA', dummy2, nprodf
     	STOP
      ENDIF

      WRITE(16,131)DELS,S(1),HSUBS,HSTEP,HMIN
 131  FORMAT(1X,'DELS,S1,HSUBS,HSTEP,HMIN',5E12.3)
!* From xgridnew_extend.jup:
      READ(7,601) NNSPEC
!C     NNSPEC IS THE NUMBER OF SPECIES BEING CONSIDERED.
  601 FORMAT(I5)
      READ(7,750)M6,NBINS
!C 		M6 is the number of different bin sizes, NBINS is the number of energy bins
      READ(7,752)(NMM(II),DDDE(II),II=1,M6) !number of bins of size NMM and corresponding size DDDE
      READ(7,602) (SPEC(N),N=1,NNSPEC)
!C       SPEC(I) is a two letter designation identifying each species.
!****** 1- H2; 2- He; 3- H; 4- CH4

  602 FORMAT(10A3)

!*Read in the neutral/electron temperature at each altitude bin
	  print*, 'check that you are using the right T profile'
      ! DO I=1,NS
      ! 	READ(17,*) Z, TE(I)
      ! ENDDO

!C      ======================================================
!C                 SECTION 3. Initial Setup
!C                 ------------------------

      DO 515 I=1,IMAX
         AREA(I)=1.0
 515  GAMMA(I)=0.
!C      $3.1: Set up energy grids.
!C            --------------------
          EEE=0.0
          IV=0
      DO 754 II=1,M6
         DD7=DDDE(II)
         ND7=NMM(II)
         DO 756 II2=1,ND7
             EEE=EEE+DD7
             IV=IV+1
             DELE(IV)=DD7
  756    EE(IV)=EEE-.5*DD7
  754 CONTINUE
!**** Un-comment this statement if you want to see the list of energy bins
!c  		do I = 1,nbins
!c         	print*, I, EE(I)
!c       enddo

      IF (IV .EQ. NBINS) GO TO 755
          WRITE(16,757)
  757 FORMAT('ERROR. # of bins does not agree with value for JMAX.')
  755 CONTINUE
      IF(IV .NE. NBINS) STOP 'iv .ne. nbins; stop'
      WRITE(34,1955)NNSPEC,M6,NBINS,IMAX,IPMAX,JMAX,JMIN,IV
 1955 FORMAT(2X,'NNSPEC,M6,NBINS,IMAX,IPMAX,JMAX,JMIN,IV',/,8I6)
      write(16,*)' EE(j), j=1,',nbins
      WRITE(16,133)(EE(IEGAN),IEGAN=1,NBINS)
 133  FORMAT(1X,10E12.3)

      DO I=2,IMAX
      	S(I)=S(I-1)+DELS
      enddo

!C         $3.2:
!C     --------- The followng loop computes the radius of all the spatial
!CCC    grids with arc length S(I) and subsolar altitude HSUBS

      DO 151 I=1,IMAX
        ARCLTH=S(I)
	  	IF (PARA .EQ. 1) THEN
           CALL PARA2(HSUBS,ARCLTH,RADIUS)
           RAD(I)=RADIUS
           ramagl(i)=ramdir
        ELSE
	   		RAD(I)=S(I)+HSUBS
!	   		print*, rad(i), i
		END IF
 151  CONTINUE

!C             $3.3: read in AIRGLOW CS.
!C----------------------------------------------------------------
!C         READ IN CROSS-SECTIONS FOR THE PRODUCTION OF VARIOUS
!C 	AIRGLOW FEATURES,THEN PUT THEM INTO THE LABELED
!C 	COMMON/RABBIT/ FOR USE BY THE FUNCTION SR SIGFE
!C..................................
      READ(13,450)(NAPP(IV),IV=1,nnspec) !Number of processes for each species
      DO 452 IB=1,nnspec
        NAPPP=NAPP(IB)
        DO 454 IBB=1,NAPPP
!c        print*, ibb
          READ(13,459)NXX,aglwlab(IB,ibb) !read the number of points and the label
!c          print*, nxx
  459	  format(i5,13x,a11)
          NCRAB(IBB,IB)=NXX
          READ(13,457)(ECRAB(IV,IBB,IB),IV=1,NXX) !read the energy points for the airglow xs
!c          print*, (ECRAB(IV,IBB,IB),IV=1,NXX)
          READ(13,458)(CRAB(IV,IBB,IB),IV=1,NXX) !read the airglow xs for each species at each energy
!c          print*, (CRAB(IV,IBB,IB),IV=1,NXX)
          DO 453 II=1,ns
  453     AGLW(II,IB,IBB)=0.0
  454   CONTINUE
  452 CONTINUE

  450 FORMAT(4I5)
  458 FORMAT(8E10.3)
  457 FORMAT(10F8.2)

!CCCCCC     $3.4: READ ELASTIC CROSS-SECTIONS from xsect output.
!C--------------------------------------------------------
      DO 701 N=1,NNSPEC
      	READ(9,304)(SIGS(N,J),J=1,NBINS)
      	READ(9,355)(PE(N,J),J=1,NBINS)
  701 READ(9,355)(PI(N,J),J=1,NBINS)


!*Add do loops to turn off the backscattering, i.e. set pe and pi eq. 0.0 if wanted for comparison NOM
!C	  DO N=1,NNSPEC
!C	  	DO J = 1,NBINS
!C	  		PE(N,J) = 1E-20
!C	  		PI(N,J) = 1E-20
!C	  	ENDDO
!C      ENDDO

!C       NOTE:  SINCE I DO NOT READ IN RKM, I WILL SET THE ARRAY HERE
!C       RKM SHOULD BE IN CM
!C       RKM = HSUBS + DELTA_Z
!C NOM   I think that this is only used for the comet code. commented out for Jupiter
        RKM(1) = HSUBS
        DO I = 2, IPMAX
           RKM(I) = RKM(I-1) + HSTEP
        ENDDO

!C       SET THE NEUTRAL AND ELECTRON DENSITIES AND ELECTRON TEMP
!C       I USED RAD(I), WHICH IS IN CM
!c nom   FIONZ = 1.0E-6/DHEL/DHEL   NOT USED IN THE JUPITER CODE, ONLY FOR COMETS
        print*, 'elect code: reading atmospheric data.'
!      CALL ATMDATA !Gets the data points to interpolate for the neutral species in Jupiter's atmos.
!	  do i =1,7
!     	read(19,*) !skip lines in the jpce output to then read the electron densities  !!!!!Commented out to not use Chemical Code Input
!      enddo
! open(unit=223,file='../Atmosphere/Input/JunoAtmosphere_2km.dat',status='old')
open(unit=223,file='../Atmosphere/Input/Jupiter_Atmosphere_EqMix.dat',status='old')
read(223,*)
! DO IP = 1, IPMAX !loop over all altitude increments (Normal Atmosphere)
do IP=IPMAX,1,-1 !Equal mixing atmosphere
	R = (RAD(IP)*1E-5)
	! read(223,*) dummy1,H2,HE,CH4,H,dummy1,TE(IP) !Normal atmosphere
	read(223,*) dummy1,H2,HE,CH4,dummy1,TE(IP) !Equal mixing atmosphere
  H=1.0E-15
	ZSPEC(1,IP) = H2
	ZSPEC(2,IP) = HE
	ZSPEC(3,IP) = H
	ZSPEC(4,IP) = CH4
!* This read statement should only be done after all the codes have been run once with the guessed density
!C		  read(19,*) alt, sne(ip)
!c		   read(19,555) sne(ip) !Read the electron densities produced by jpce
!c 555  format(5x,1E9.3)
!C 		print*, sne(ip)
!C nom            SNE(IP) = (FIONZ*NEUTDENS/DELTA)**0.5
!c				 SNE(IP) = 100 !This should be the guessed density used initially, i.e. the first time you run the code
	SNE(IP) = 1e5
ENDDO
!* output the density into a file to check it
      WRITE(51,*) 'DENSITIES VS R           H2               HE              H                CH4'
DO IP = 1, IPMAX
	WRITE(51,*) RAD(IP),ZSPEC(1,IP), ZSPEC(2,IP),ZSPEC(3,IP),ZSPEC(4,IP)
ENDDO

!c	$3.4: setup solar wind flux as boundary condition
!c 		-----------------------------------------
      FAC=0.
!c      	print*, EE(242)
!c       write(54,*) EE(20)
!* For Monoenergetic (auroral) electron beams chose the energy bin (i)
!* for desired energy and then input the flux for phiinf(i). The rest of the bins will be set to 0.0
      DO 765 I=1,JMAX
	  	PHIINF(I)=0.0 !No magnetospheric electrons.
!*****SJH Uncomment the next lines for a monogenergetic beam
!* Here I have allowed for a quick normalization to 1 mW/m^2. It takes two runs,
!* But the second run tells you want is needed for phiinf to have 1 mW/m^2.
			! IF(I .EQ. 226)THEN
			! 		print*,''
      ! 	  print*, 'energy bin for monoenergetic beam [keV] = ', EE(I)*100.0
      !    PHIINF(I) = 64015.4766 !Means: Input flux of 'something'cm-2s-1eV-1 at E bin chosen
      !    ! I = 185 and phiinf = 1273777.38 for 5keV
      !    ! I = 210 and phiinf = 630455.438 for 10keV
      !    ! I = 226 and phiinf = 64015.4766 for 20keV
      !    ! I = 242 and phiinf = 4847.77393 for 50keV
      !    ! I = 251 and phiinf = 1260.91101 for 100keV
      !    ! I = 260 and phiinf = 660.477173 for 200keV
      !    print*, 'Energy flux [eV-1 cm-2 s-1] = ',phiinf(i)
			! 		print*, 'Delta E [eV] = ', DEL(I)
			! 		print*, 'Energy flux [mW m-2] = ', phiinf(i)*DEL(I)*EE(I)*1.60217662e-12*0.5
			! 		print*, 'For 1 mW/m^2, use phiinf = ', 2/(DEL(I)*EE(I)*1.60217662e-12)
			! 		print*, ''
      ! ENDIF
!*****NOM Uncomment the next lines for a monogenergetic beam
!c       IF(I .EQ. 227)THEN
!c          PHIINF(I) = 4.85e4 !Means: Input flux of 'something'cm-2s-1eV-1 at E bin chosen
!c          print*, 'Energy flux= ',phiinf(i)
!c       ENDIF
!********NOM
!c	  IF (NES .GE. 0.0) THEN
!c      PHIINF(I)=3.3456E7*NES*EE(I)*EXP(-EE(I)/TES)/TES**1.5
!c	  PHIINF(I) = PHIINF(I)/(DHEL**2)
!c	  ELSE
!c	  	READ(19,*)PHIINF(I) !External flux input may be used
!c	  	PHIINF(I)=PHIINF(I)/2.0
!C ...........Only half of the s.w. electrons is going down.
!c	  ENDIF
 765  CONTINUE

      WRITE(16,*)' BC: phiinf(i) = 0.0'
      write(16,'(1x,1p10e11.2)')phiinf

        WRITE(12,800)
  800     FORMAT(3X,'       is      s(cm)           rad(cm)            mza               &
					n(H2)            n(He)           n(H)             n(CH4)...       n(e)             T(e)')
do i=1,imax
	WRITE(12,*)i,s(i),rad(i),ramagl(i),(ZSPEC(M,I),M=1,NNSPEC),SNE(I),TE(I)
 805      FORMAT(I4,1X,1p8E11.3)
enddo

!c      READ(3,*) dummy2,JLOCAL
!c      print*, dummy2, JLOCAL
!CCC
      SINDIP=SIN(DIPANG)
      RAVMU=1./AVMU
      RMUSIN=1./SINDIP/AVMU
!C  AVMU IS MEAN ANGLE MU
      SIGO=1.E-16
      QO=6.51E-14


!C      $3.5:  format storage
!C      ------------------------
 301  FORMAT(2I10)
 302  FORMAT(3F8.2,4E10.3)
 303  FORMAT(E10.3)
 304  FORMAT(8E10.3)
 305  FORMAT (10F8.4)
 306  FORMAT(1H ,5X,1HZ,11X,3HZC2,9X,3HZCO,9X,3HZOX,9X,3HZNE)
 307  FORMAT(1H ,5E12.3)
 310  FORMAT(1H ,5X,1HZ,11X,6HPHIDWN,5X,5HPHIUP,7X,5HPHIDD,7X,5HPHIUU,8X,5HGAMMA,9X,2HT1,9X,2HT2,6X,6HTSIGNE)
 311  FORMAT(1H ,9E12.3)
 314  FORMAT(1H1,15X,'THIS IS THE PHOTOELECTRON FLUX FOR ENERGY J=',F10.2,' EV')
 317  FORMAT (3F10.3)
 319  FORMAT(1H ,5X,19HTHIS IS FOR ENERGY=,F10.2,2HEV)
 320  FORMAT(1H ,7X,1HZ,13X,6HPRODUP,11X,6HPRODWN,12X,4HPROD)
 321  FORMAT(1H ,5E15.7)
 323  FORMAT(1H1,10X,1HZ,8X,6HSECION,6X,6HSC2ION,6X,6HSCOION,6X,6HSOXION)
 324  FORMAT(1H1,10X,1HZ,8X,6HSECION,6X,5HTPROD,61X,6HSC2ION,6X,6HSUMDIF,4X,6HSUMNET)
 325  FORMAT(1H ,5X,1P5E12.3)
 326  FORMAT(1H ,5E15.3)
 327  FORMAT(1H1,15X,6HC2SIGA,12X,6HCOSIGA,12X,6HOXSIGA)
 328  FORMAT(1H ,15X,3E15.4)
 329  FORMAT(1H ,3E15.3)
 330  FORMAT(5I10)
 331  FORMAT(1H ,E12.3)
 332  FORMAT(1H ,2E20.5)
 343  FORMAT(10E8.2)
 345  FORMAT(1X,5F12.4)
 344  FORMAT(1X,1P10E12.2)
 355  FORMAT(8F10.6)
 750  FORMAT(6I5)
 776  FORMAT(6I5)
 752  FORMAT(1X,8(I3,F7.2))
 778  FORMAT(8E10.3)
 804  FORMAT(5X,3E9.2)
 806  FORMAT(4X,1PE10.3,4X,1PE10.3,5X,1PE10.3)
 901  FORMAT(1X,1P10E11.3)

!C  	$3.6:   preset stage ........
!C-----------------------------------
      ALPHA(1)=0.
      BETA(1)=0.
      GAMMA(1)=0.
      J=JMAX
      KPROD=1
      DO 16 I=1,IMAX
        EHEAT(I)=0.0
        TPROD(I)=0.0
        PHINET(I)=0.
        SUMUP(I)=0.
        SUMDWN(I)=0.
        SUMNET(I)=0.
        DO 17 N=1,NNSPEC
   17   SION(N,I)=0.0
   16 SECION(I)=0.0

      DO 22 I=1,IMAX
        DO 21 JJ=1,JMAX
          PRODUP(I,JJ)=1.0E-20
   21     PRODWN(I,JJ)=1.0E-20
   22 CONTINUE

      NBINS=NBINS+1
			write (*,*) nbins
      IF(NPRODF .NE. 1) GO TO 23
		print*, 'p.e. loop'
!C      $2.7  read in photoelectron prod. rate which is a function of
!C        alt. (ip) and energy (j).
!C-----------------------------------------------------------------
	  ! print*, 'reading photoel. files'
	  ! read(11,*)
	  ! read(11,*)
	  ! read(11,911)IPMAX,JBINS
!C.........................jbins: energy bins used in ephot.
 ! 911  FORMAT(2I6)
 !
	!   if (jbins .lt. jmax) then
	!   	write(6,*)'jbins,jmax=',jbins,jmax
	!   	stop 'not enough energy bins in ephot! '
	!   endif
 !
 !      IF(JBINS .eq. JMAX) GO TO 125
 !
 ! 731   READ(11,912) (dumm(I),I=1,IPMAX)  ! dummy
 912  FORMAT(1X,1P10E11.3)
 !      JBINS=JBINS-1
 !      IF(JBINS .GT. JMAX) GO TO 731
!C.......................

!C**********-------------
 125    CONTINUE
!****nom commented out so we don't read these photoelectron production rates
!c      DO JL=1,JBINS
!c      	J1=JBINS-JL+1
!c        read(11,912)(peprod(I,j1),I=1,IPMAX)
!c	  enddo
!c	  open(80, file='pecheck.chk', status='unknown')
!c	  	DO JL=1,JBINS
!c     	J1=JBINS-JL+1
!c        write(80,912)(peprod(I,j1),I=1,IPMAX)
!c	  enddo
!***nom Read in the electron production from the ion precipitation
!		print*, 'reading sec elect prod files'
		eprodf = 0.0;sprodf = 0.0
		eprofb = 0.0;sprofb = 0.0
	  do JL=1,260
!c	  	J1=JBINS-JL+1
!		J1 = JL !This is ONLY done for the sec elects from Ion MC code because of output format
	   read(20,912)(sprodf(I,jl),I=1,IPMAX)
	   read(22,912)(sprodb(I,jl),I=1,IPMAX)
	  enddo
    write(*,*) sum(sprodf),sum(sprodb),sum(sprodf)+sum(sprodb)
    ! stop
		do i=1,ipmax !Looping through altitude bins
		  do jen=1,ne !Looping through energy bins
		    eprodf(i,jen) = (sprodf(i,jen)/del(jen)) !nz&ne+1 to not include labels on first row and first column
		    eprodb(i,jen) = (sprodb(i,jen)/del(jen))
		  enddo
		enddo
	  open(80, file='secelcheck.chk', status='unknown')
123 format(A11,8(11x,'x',F5.0))
		write(2323,123) ' # Energy', rad(195)*1E-5, rad(220)*1E-5, rad(245)*1E-5, rad(295)*1E-5,&
		rad(420)*1E-5, rad(545)*1E-5, rad(795)*1E-5, rad(1045)*1E-5
	  DO JL=1,jmax
      write(80,912)(sprodf(I,jl),I=1,IPMAX)
			if (sprodf(220,jl).gt.0.0) then
				write(2323,*) ee(jl), eprodf(195,jl), eprodf(220,jl), eprodf(245,jl),&
				eprodf(295,jl), eprodf(420,jl), eprodf(545,jl), eprodf(795,jl), eprodf(1045,jl)
			end if
		enddo
	  write(80,*) '*********************************'
	  DO JL=1,jmax
        write(80,912)(eprodb(I,jl),I=1,IPMAX)
	  enddo
!C	**************************************************
!C	C                                               CC
!C	C   SECTION 4:   ENERGY LOOP BEGINS HERE        CC
!C	C						CC
!C	C*************************************************

   23 CONTINUE		! J loop.

!	 print*, 'reading sec elect prod files'

!CCCCCCC ... READ ON GRID INELASTIC X-SECTIONS AT J FOR EACH K LOSS
      READ(9,201) (ID(N),N=1,nnspec),IIMIN
 201  FORMAT(6I5)

      DO 600 I=1,NNSPEC
        IIN=ID(I)
  600 READ(9,778)(SIGA(I,IV),IV=1,IIN)
!c      print*, i, (siga(i,iv),iv=1,iin)
      READ(9,776)IIMAXX
      READ(9,*)(TION(N),N=1,NNSPEC)

      DO 605 N=1,NNSPEC
  605 READ(9,778)(SEC(N,K),K=1,IIMAXX)
!C..............................................

!C.........keeping track of the energy bin of XSECT data
      NBINS=NBINS-1
      IF(NBINS .GT. JMAX) GO TO 23
!C.........clean prod for this j loop...
      DO 6 I=1,IMAX
   6  PROD(I)=1.0E-20

!C.........decide wether to include photoelectron production
      IF (NPRODF-1) 18,8,8
   8  CONTINUE

!C.........prod(i) at s(i) is extrapolated assuming spherical symmetric
!C         alt. dependent p.e. production. Linear extrapolation is used.
!c....nom	  do ip=1,ipmax
!c	   	peprod1(ip)=peprod(ip,j)
!c	  enddo
	  do ip=1,ipmax
!c	  	seprod1(ip) = (eprodf(ip,j)+eprodfR2(ip,j))/2.
!c	  	seprod1(ip) = (eprodf(ip,j)+eprodfR2(ip,j)+eprodfR3(ip,j))/3.
	  	seprod1(ip) = (eprodf(ip,j))
!c	  	seprod2(ip) = (eprodb(ip,j)+eprodbR2(ip,j))/2.
!c	  	seprod2(ip) = (eprodb(ip,j)+eprodbR2(ip,j)+eprodbR3(ip,j))/3.
	  	seprod2(ip) = (eprodb(ip,j))
	  enddo

	  ires=1
	  do I=1,IMAX
!c...nom for photoes	   call extrap(ipmax,rkm,peprod1,rad(i),prod(i),1,ires)
		call extrap(ipmax,rkm,seprod1,rad(i),prodf(i),1,ires)
		call extrap(ipmax,rkm,seprod2,rad(i),prodb(i),1,ires)
!C       EXTRAPOLATION ONLY FOR COMETS
!c nom        IF(RAD(I).GT. RKM(IPMAX)) THEN
!c nom           PROD(I) = PEPROD1(IPMAX)*(RKM(IPMAX)/RAD(I))**2
!c nom        ENDIF
!C       END EXTRAPOLATION FOR COMETS ONLY

      	ires=0
!c      	PROD(I)=AREA(I)*PROD(I)/DELE(J)*(RAVMU)/SINDIP
      	prodf(I)=AREA(I)*prodf(I)/DELE(J)*(RAVMU)/SINDIP
      	prodb(i)=AREA(I)*prodb(I)/DELE(J)*(RAVMU)/SINDIP
	  enddo

   18 KK=J
!C...........flag to write phi+,- .....
      IF(J .GE. (JMAX-5)) GO TO 855
      IF(J .LE. 4) GO TO 855
      IF(MOD(J,5) .NE. 0) GO TO 850
  855 WRITE(16,319) EE(J)
      WRITE(16,320)
      DO 519 I=1,IMAX
        IF(MOD(I,2) .EQ. 1)then
					WRITE(16,321) S(I),PRODUP(I,KK),PRODWN(I,KK),PRODf(I),prodb(I)
				endif
  519 CONTINUE

!C.....................
  850 CONTINUE
!C..............set other part of siga
      DO 781 I=1,NNSPEC
        III=ID(I)+1
        DO 782 IV=III,IIMIN
  782     SIGA(I,IV)=1.0E-30
  781 CONTINUE

!C.....
      DO 783 I=1,NNSPEC
          ALOSS(I)=0.0
  783 TSA(I)=0.0

!C..... calculate total cs at this energy and, if necessary, the loss ftn.
      L=J-1
      DO 25  K=1,L
        LJK=J-K
        IF (LJK .LT. 1) LJK = 1
        DO 784 I=1,NNSPEC
  784   TSA(I)=TSA(I)+ SIGA(I,K)*(DELE(LJK)/DELE(J))
!C........... SIGA(I,K) is the CS for (K)eV of energy loss
   25 CONTINUE

!C      DAG=EE(J)-EE(JJJ4)
	  dag=dele(j)		! this statement avoids a dip at the change
!C                                 of different bin size.

!C.........compute e - e equivalent Coulomb CS.
DO 11 I=1,IMAX
	ET=8.618E-5*TE(I)
	EET=EE(J)-ET
	IF (EET .gt. 0.0) then
		TSIGNE(I)=((3.37E-12*sne(I)**0.97)/(EE(J)**0.94))*((EET)/(EE(J)-(0.53*ET)))**2.36
		TSIGNE(I)=TSIGNE(I)*(RAVMU)/SINDIP/DAG
	else
		TSIGNE(I)=0.0
	endif

   11 CONTINUE
!C,,,,,,,,, dicide if to write to flux.dat
      IF (J .GE.(JMAX-5)) GO TO 820
      IF(J .LE. 4) GO TO 820
      IF(MOD(J,5).NE. 0) GO TO 821
  820    WRITE(16,856)
  856    FORMAT(//31HTOTAL ENERGY LOSS CROSS-SECTION)
         WRITE(16,902)
  902    FORMAT(1X,'TSA: H2, He, H, CH4')
         WRITE(16,901)(TSA(N),N=1,NNSPEC)
         WRITE(16,857)
  857    FORMAT(//27HTOTAL ELASTIC CROSS-SECTION)
         WRITE(16,903)
  903    FORMAT(1X,'SIGS: H2, He, H, CH4')
         WRITE(16,901)(SIGS(N,J),N=1,NNSPEC)
  821 CONTINUE

!C..........start solving equations..........
      DO 30 I=1,IMAX
         T1(I)=0.0
         T2(I)=0.0
         DO 31 IV=1,NNSPEC
            T1(I)=T1(I)+ZSPEC(IV,I)*SIGS(IV,J)*PE(IV,J)
   31       T2(I)=T2(I)+ZSPEC(IV,I)*(SIGS(IV,J)*PE(IV,J)+TSA(IV))
         T1(I)=T1(I)*(RAVMU)/SINDIP
         T2(I)=T2(I)*(RAVMU)/SINDIP +TSIGNE(I)
   30 CONTINUE

      IF(JLOCAL .EQ. 1) GO TO 449
      I1=IMAX-1
DO 40 I=2,I1
	PHI(I)=1.0
	ALPHA(I)=-(T1(I+1)-T1(I-1))/DELS/T1(I)/2.
	BETA(I)=T2(I)*(T1(I+1)-T1(I-1))/T1(I)/DELS/2.-(T2(I+1)-T2(I-1))/DELS/2.-T2(I)**2+T1(I)**2
!C
!C  ...........**  THIS IS FOR ZENITH ANGLES .LE. 90 DEG.  **
!C  ...........**   TO AVOID DIVIDING BY ZERO PROD OR PRODWN  **
!C
!c...nom         IF(PROD(I).LT. 1.E-30) PROD(I)=1.E-30
IF(prodf(I).LT. 1.E-30) prodf(I)=1.E-30
IF(prodb(i) .lt. 1.e-30) prodb(i)=1.e-30
!c         IF (PRODWN(I,J).LT. 1.E-30) PRODWN(I,J)=1.E-30
!C
!c ...nom        GAMMA(I)=(PROD(I)/2.0)*(-T1(I)-T2(I)-ALPHA(I)-(PROD(I+1)-
!c ...nom        PROD(I-1))/PROD(I)/DELS/2.0)+PRODWN(I,J)*(-ALPHA(I)-T2(I)-
!changed this to use the forward scattered flux calculated by ion model for sec-elecs
GAMMA(I)=(prodf(I))*(-T1(I)-T2(I)-ALPHA(I)-(prodf(I+1)-prodf(I-1))/prodf(I)/DELS)
GAMMA(I)=GAMMA(I)+PRODWN(I,J)*(-ALPHA(I)-T2(I)-(PRODWN(I+1,J)-PRODWN(I-1,J))/PRODWN(I,J)/DELS/2.0)
GAMMA(I)=GAMMA(I)+PRODUP(I,J)*(-T1(I))
!c            if (produp(i,j) .lt. 0) then
!c            	print*, produp(i,j)
!c            endif
   40 CONTINUE

!C**********------------------
!C           new boundary conditions.
!C
!C    check this section for convergence ***
!C                     --------------------------
!c      PHIDWN(J,1)=PHIINF(J)
      TTT=T2(1)-T1(1)
      IF (T2(1) .LE. T1(1)) TTT=1.0E-12
!c  DON'T UNCOMMENT THIS LINE: IF(PHIINF(J) .LT. 1.E-10) THEN
!C NOM         PHIDWN(J,1)=(PROD(1)/2.+PRODWN(1,J))/TTT
!c DON'T UNCOMMENT THIS LINE:     ENDIF
	  PHIDWN(J,1)=(prodf(1)+prodwn(1,j))/TTT	!changed this to use forward sec el flux as downward
      DUMMY(1)=PHIDWN(J,1)
      I3=0
      I3MAX=25
      Q(1)=0.0
      PHHD(1)=DUMMY(1)
      EPS=1.0E-04
41    I3=I3+1
IF (I3 .EQ. 1) GOTO 42
IF (I3 .EQ. 2 .AND. Q(1) .LT. 0.0) PHHD(2)=0.5*PHHD(1)
IF (I3 .EQ. 2 .AND. Q(1) .GE. 0.0) PHHD(2)=2.0*PHHD(1)
IF (I3 .LT. 3) GO TO 42
PHHD(I3)=PHHD(I3-2)-Q(I3-2)*((PHHD(I3-1)-PHHD(I3-2))/(Q(I3-1)-Q(I3-2)+1.0E-09))
TX=(PHHD(I3)-PHHD(I3-1))/(PHHD(I3)+1.0E-09)
TX=ABS(TX)
IF (TX .LE. EPS) GO TO 43
IF (PHHD(I3) .LT. 0.0) PHHD(I3)=PHHD(I3-1)/10.
42    CONTINUE
DUMMY(1)=PHHD(I3)
PHIDWN(J,1)=DUMMY(1)

CALL IMPIT(DUMMY)

!c... nom      R1=(T1(2)*DUMMY(2)+(PROD(2)+2.0*PRODUP(2,J))/2.)/T2(2)
	  R1=(T1(2)*DUMMY(2)+(prodb(2)+PRODUP(2,J)))/T2(2) !changed this to use bckscattered el flux for upward flux
      PHIUP(J,2)=R1+(DUMMY(1)-R1)*EXP(-T2(2)*DELS)
      Q(I3)=(PHIUP(J,2)+DUMMY(2)-2.0*PHHD(I3))/(PHHD(I3)+1.0E-09)
      EXPP=3.0E-04
      WRITE(16,919)I3,PHIDWN(J,1),PHIUP(J,2),DUMMY(2),Q(I3)
      QQ=ABS(Q(I3))
      IF (QQ .LT. EXPP) GO TO 43
      IF (I3 .GE. I3MAX) GO TO 43
      GO TO 41
  43  CONTINUE

 919  FORMAT('  convergence test :',I4,1P4E10.3)
!CCC                                   -----------------------------------
      DO 50 I=1,IMAX
   50 PHIDWN(J,I)=DUMMY(I)
      PHIUP(J,1)=PHIDWN(J,1)
      DO 60 I=2,IMAX
!C...nom      R1=(T1(I)*PHIDWN(J,I)+(PROD(I)+2.*PRODUP(I,J))/2.)/T2(I)
      R1=(T1(I)*PHIDWN(J,I)+(prodb(I)+1.*PRODUP(I,J)))/T2(I) !changed this to match the upward flux with back scattered flux
   60 PHIUP(J,I)=R1+(PHIUP(J,I-1)-R1)*EXP(-T2(I)*DELS)
      GO TO 451
  449 CONTINUE
      DO 616 I=1,IMAX
      IF (T2(I) .LE. T1(I))then
      	 T2(I)=T1(I)*1.00001
      endif
!c...nom      PHIUP(J,I)=(PROD(I)/2.0+PRODUP(I,J))/(T2(I)-T1(I))
!c ...nom      PHIDWN(J,I)=(PROD(I)/2.0+PRODWN(I,J))/(T2(I)-T1(I))
	  if (T2(i)-T1(i)  .le. 0.) then
!	  	print*, T2(i)-T1(i),i
	  endif
      PHIUP(J,I)=(prodb(I)+PRODUP(I,J))/(T2(I)-T1(I)) !Electrons leaving the atmosphere (going to higher altitudes)
      PHIDWN(J,I)=(prodf(I)+PRODWN(I,J))/(T2(I)-T1(I)) !Electrons going deeper into the atmosphere (going to lower altitudes)
  616 CONTINUE
  451 CONTINUE
      DO 118 I=1,IMAX
      PHINET(I)=(PHIUP(J,I)-PHIDWN(J,I))*AVMU
      SUMUP(I)=SUMUP(I)+PHIUP(J,I)
118   SUMDWN(I)=SUMDWN(I)+PHIDWN(J,I)

	  DO I = 1, IMAX
	    if (phidwn(j,i) .lt. 0.0) then
			phidwn(j,i) = -1.0*phidwn(j,i)
      	endif
	     IF(PARA .EQ. 1) THEN
!C         PARABOLIC:
	      UFX(2*I-1,J)=PHIUP(J,I)/2.0/pii
	      UFX(2*I,J)=PHIDWN(J,I)/2.0/pii
	     ELSE
!C       RADIAL:
  	    UFX(2*I-1,J)=PHIUP(J,I)/2.0/pii
	      UFX(2*I,J)=PHIDWN(J,I)/2.0/pii
	     ENDIF
	     if(ufx(2*i,j) .lt. 0.0)then
!	     	print*, ufx(2*I,J),i,j
	     endif
	     if(ufx(2*i-1,j) .lt. 0.0)then
!	     	print*, ufx(2*I-1,J),i,j, phiup(j,i),produp(i,j)
	     endif

	  ENDDO

!C........write output
      IF (J .GE.(JMAX-5)) GO TO 854
      IF (J .LE. 4) GO TO 854
      IF(MOD(J,5).NE.0) GO TO 852
  854  WRITE(16,314)EE(J)
       WRITE(16,310)
       DO 513 I=1,IMAX
	     if ( (i.le.3).or.(i.ge.(imax-2)).or.(mod(i,10).eq.1) ) then
         PHIDD=PHIDWN(J,I)/AREA(I)
         PHIUU=PHIUP(J,I)/AREA(I)
         WRITE(16,311) S(I),PHIDWN(J,I),PHIUP(J,I),PHIDD,PHIUU,GAMMA(I),T1(I),T2(I),TSIGNE(I)
	     endif
  513  CONTINUE
  852 CONTINUE

!C...........write to 28: flux3d
write(28,314)ee(j)
write(28,390)
 390  FORMAT(' r/s  phiup(cm-2s-1eV-1sr-1)  phidwn(cm-2s-1eV-1sr-1) ')
write(28,391)(s(i),phiup(j,i)/(2*pii),phidwn(j,i)/(2*pii),i=1,imax)
 391	format(1p3e12.3)

!C............record produp, proddwn for lower energy levels
99    DO 71 I=1,IMAX
      L=MAX(ID(1),ID(2))-1
!CCCCC we found that L = 0 for low energies, which messes up the do loop
!CCCCC to correct we added the following if loop
      IF (L .EQ. 0)THEN
      	L=1
      ENDIF

DO 70 K=1,L
	LL=J-K
	IF( LL .LT. 1) GO TO 70
		PRODDA=0.0
		PRODUA=0.0

		DO 72 N=1,NNSPEC
			PRODUA= PRODUA +ZSPEC(N,I)*(SIGA(N,K)*PI(N,J)*PHIDWN(J,I)+(1.-PI(N,J))*SIGA(N,K)*PHIUP(J,I))
			PRODDA=PRODDA+ZSPEC(N,I)*(SIGA(N,K)*PI(N,J)*PHIUP(J,I)+(1.-PI(N,J))*SIGA(N,K)*PHIDWN(J,I))
            if(produa .lt. 0)then
!print statement            	print*, 1-pi(n,j),pi(n,j),n,j
            endif
   72 	CONTINUE
      	PRODUP(I,LL)=PRODUP(I,LL)+PRODUA *(RAVMU)/SINDIP
      	PRODWN(I,LL)=PRODWN(I,LL)+PRODDA*(RAVMU)/SINDIP
   70 CONTINUE
   71 CONTINUE
      KK=J-1
      JJJ=J-1
      IF (KK .LE. 0) GO TO 73
      IF(JJJ.LE.0)JJJ=1
      DO 75 I=1,IMAX
      	PRODUP(I,KK)=PRODUP(I,KK)+TSIGNE(I)*PHIUP(J,I)*(DELE(J)/DELE(JJJ))
   75   PRODWN(I,KK)=PRODWN(I,KK)+TSIGNE(I)*PHIDWN(J,I)*(DELE(J)/DELE(JJJ))
   73 CONTINUE
      DAG=DELE(J)
      DO 76 I=1,IMAX
!c...nom      TPROD(I)=TPROD(I)+PROD(I) !I don't think this is used anywehere else in code
		 tprodf(i) = tprodf(i) + prodf(i) !Changed it to add forward prod
		 tprodb(i) = tprodb(i) + prodb(i) !Changed it to add backward prod
   76 EHEAT(I)=EHEAT(I)+TSIGNE(I)*(PHIUP(J,I)+PHIDWN(J,I))*DAG*DAG/AREA(I)
!CCCCCCCCC  CALCULATION OF AIRGLOW FEATURES
      EBMM=EE(J)
      NAPPP=NAPP(1) !For species 1 = H2
      DO 470 IBB=1,NAPPP
!* for hydrogen we have H2 dissociation that leads to h lyman emission.
!* the sigfe() function will give only the slow portion of this Hlyman emission
!* to get the total emission this must be added to the Hfast cross section
!* calculated by SIGDIS
     	if(IBB .eq. 3)then
     		CALL SIGDIS(EBMM,SIGHFAST)
        	SIGAG = (SIGFE(EBMM,IBB,1)+ SIGHFAST)*DELE(J)
        else
        	SIGAG=SIGFE(EBMM,IBB,1)*DELE(J)
        endif
      	DO 471 II=1,IMAX
  471 		AGLW(II,1,IBB)=AGLW(II,1,IBB)+(PHIUP(J,II)+PHIDWN(J,II))*SIGAG*ZSPEC(1,II)/AREA(II)
     	if(IBB .eq. 4)then
     		AGLW(II,1,IBB) = 0.95 * AGLW(II,1,IBB) !This will give the cascade emission.
				! stop
!* Might have to play with the percent that cascades to Ly alpha band
		endif
  470 	CONTINUE
			write(53,*) EBMM, (SIGFE(EBMM,ibb,1), ibb=1,nappp)
      NAPPP=NAPP(2) !For species 2 = He
      DO 475 IBB=1,NAPPP
      	SIGAG=SIGFE(EBMM,IBB,2)*DELE(J)
      DO 477 II=1,IMAX
  477 	AGLW(II,2,IBB)=AGLW(II,2,IBB)+(PHIUP(J,II)+PHIDWN(J,II))*SIGAG*ZSPEC(2,II)/AREA(II)
  475 CONTINUE
      NAPPP=NAPP(3) !For species 3 = H
      DO 481 IBB=1,NAPPP
      	SIGAG=SIGFE(EBMM,IBB,3)*DELE(J)
      DO 482 II=1,IMAX
  482 	AGLW(II,3,IBB)=AGLW(II,3,IBB)+(PHIUP(J,II)+PHIDWN(J,II))*SIGAG*ZSPEC(3,II)/AREA(II)
  481 CONTINUE
      NAPPP=NAPP(4) !For species 4 = CH4
      DO 491 IBB=1,NAPPP
      	SIGAG=SIGFE(EBMM,IBB,4)*DELE(J)
      DO 492 II=1,IMAX
  492 	AGLW(II,4,IBB)=AGLW(II,4,IBB)+(PHIUP(J,II)+PHIDWN(J,II))*SIGAG*ZSPEC(4,II)/AREA(II)
  491 CONTINUE
!CCCCCCCCCCC   READ IN IONIZATION X-SECT FOR SECONDARIES AT K BIN
!CCCCCCCCCCC          FOR ENERGY J (BIN)       ADD TO PRODUCTION
      DO 87 K=1,IIMAXX
      	DO 86 I=1,IMAX
      		DO 85 N=1,NNSPEC
      			SECP(N)=SEC(N,K)*ZSPEC(N,I)*(PHIUP(J,I)+PHIDWN(J,I))
      			SION(N,I)=SION(N,I)+SECP(N)*DELE(K)
      			SECION(I)=SECION(I)+SECP(N)*DELE(K)
      			PRODUP(I,K)=PRODUP(I,K)+(SECP(N)*.5*RMUSIN)
   85 			PRODWN(I,K)=PRODWN(I,K)+(SECP(N)*.5 *RMUSIN)
   86 	CONTINUE
   87 CONTINUE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      EQ1=EHEAT(6)*AVMU*SINDIP/AREA(6)
      EQ2=EHEAT(21)*AVMU*SINDIP/AREA(21)
      IF(J .GE.(JMAX-5)) GO TO 822
      IF(J .LE. 4) GO TO 822
      IF(MOD(J,5) .NE. 0) GO TO 823
  822 WRITE(16,901)(TION(N),N=1,NNSPEC),EQ1, EQ2
  823 CONTINUE
 535	format( ' j:',i4,'. cpu:',f9.2,'s.')
	  J=J-1
      IF (J .LT. JMIN) GO TO 80
      GO TO 23
!C	*****************************************
!C       *					*
!C	*	ENERGY LOOP ENDS HERE     	*
!C	*					*
!C	*****************************************

   80 CONTINUE

      print *,' post-loop calculation...'

!C***********************************************************************
!C***********************************************************************
!C  Following segment is commented out to reformat phi.chk data file.
!C***********************************************************************
!C***********************************************************************
!C
!C      WRITE(14,812)NPRODF,NES,TES,HSUBS
!C  812 FORMAT(1H ,5X,'FIRST PHIUP AND THEN PHIDWN AT SELECTED ALTITUDES',
!C      & /5X,'NPRODF,NES,TES,HSUBS= ',I2,1p3E11.3)
!C      DO 815 J4=1,6
!C      ZL4=0.
!C      IF(J4 .EQ. 1) ZL4=1.E-5*S(2)
!C      IF(J4 .EQ. 3) ZL4=1.E-5*S(50)
!C      IF(J4 .EQ. 5) ZL4=1.E-5*S(151)
!C      IF(ZL4 .NE. 0) GO TO 801
!C      GO TO 803
!C  801 WRITE(14,802) ZL4
!C  802 FORMAT(2X,11HALTITUDE = ,1PE9.2,3H KM)
!C  803 WRITE(14,810)(UFX(J4,II),II=JMIN,JMAX)
!C  810 FORMAT(1P10E11.3)
!C  815 CONTINUE
!C
!C***********************************************************************
!C***********************************************************************
!C         The following will reformat phi.chk data file to output
!C    Energy, UFX	for plotting.  In addition, altitude is
!C    adjusted to 725km above Titan.
!C
!C					J.Clark 10/9/02.
!C
!C***********************************************************************
!C***********************************************************************

WRITE(14,812)NPRODF,NES,TES,HSUBS
812 FORMAT(1H ,5X,'FIRST PHIUP AND THEN PHIDWN AT SELECTED ALTITUDES',/5X,'NPRODF,NES,TES,HSUBS= ',I2,1p3E11.3)
      DO 815 J4=1, 2*imax
	 	ZL4 = 0
	 	IF(MOD(J4,2).EQ.1) THEN
	    	ZL4 = RAD((J4+1)/2)
	 	ELSE
	    	ZL4 = RAD(J4/2)
	 	ENDIF

		IF(ZL4 .NE. 0) GO TO 801
      	GO TO 803
  801 	WRITE(14,802) ZL4
  802   FORMAT(2X,'RADIUS = ', 1PE13.6, 2X, 'CM')
  803	WRITE(14,810)(EE(II), UFX(J4,II), II=JMIN,JMAX)
!c  		do ii =jmin, jmax
!c  			if(UFX(J4,II) .lt. 0.0)then
!c  				print*, ii, ufx(j4,ii)
!c  			endif
!c  		enddo
  810 	FORMAT(1P2E10.3)
  815 CONTINUE
!C
!C***********************************************************************
!C  End modification.  J. Clark 10/9/02
!C
!C***********************************************************************

!CCCCCCCCCC   WRITE OUT AIRGLOW VOLUME EMISSION RATES
WRITE(15,*) 'RAD','               S', '          SION(1)          SION(2)', '          SION(3)'

      DO 677 II=1,nnspec
        NAPPP=NAPP(II)
        WRITE(15,674)II,NAPPP
	if (nappp .ne. 0) then
          WRITE(15,675)(aglwlab(ibb,iaglw),iaglw=1,nappp)
          DO 682 IM=1,IMAX
            SCOOP=SION(1,IM)/AREA(IM)
           IF(II .EQ. 2) SCOOP=SION(2,IM)/AREA(IM)
           IF(II .EQ. 3) SCOOP=SION(3,IM)/AREA(IM)
           IF(II .NE. 1) SCOOP=0.
           WRITE(15,679)RAD(IM),S(IM),(AGLW(IM,II,IQQ),IQQ=1,NAPPP)
  682     CONTINUE
	endif
  677 CONTINUE
  674 FORMAT(///,30X,26HSPECIES NUMBER and species,2I6,/)
  675 FORMAT(//25x,'AIRGLOW VOLUME EMISSION RATES (photons cm-3 s-1)',/,10X,'rad(cm) ','s(cm) ',7(2x,a11),'  ion.rate')
  679 FORMAT(5X,1P10E11.3)

      DO 95 I=1,IMAX
        EHEAT(I)=EHEAT(I)*AVMU*SINDIP
        SUMDIF(I)=(SUMUP(I)-SUMDWN(I))
   95   SUMNET(I)=(SUMUP(I)-SUMDWN(I))*AVMU

      DO 514 I=1,IMAX
        SECION(I)=SECION(I)/AREA(I)
        DO 514 N=1,NNSPEC
  514   SION(N,I)=SION(N,I)/AREA(I)

      WRITE(34,835)HSUBS
  835 FORMAT(4X,1PE10.3,10X,1HZ,6X,1HR,6X,21HELECTRON HEATING RATE)
      WRITE(34,806) (S(I),(RAD(I)*1E-5),EHEAT(I),I=1,IMAX)

!CCCCCC--------------------------------------------
!CCCCCC    Suprathermal Density Calculation
!CCCCCC--------------------------------------------
      DO 370 IDS=1,IMAX
        DENSTH(IDS)=0.0
        EBARTH(IDS)=0.0
!CCCC------- Bins no. 1 & 2 are not accurate. E starts at 1.25 eV ------
      DO 371 JD=3,JMAX
        DENST=(PHIUP(JD,IDS)+PHIDWN(JD,IDS))*DELE(JD)/EE(JD)**0.5
        DENSTH(IDS)=DENSTH(IDS)+DENST
        EBARTH(IDS)=EBARTH(IDS)+DENST*EE(JD)
 371  CONTINUE
        EBARTH(IDS)=2.0/3.0*EBARTH(IDS)/DENSTH(IDS)
!C----- EBAR is now actually the temperature corresponding to energy if Max'ed.--C
      DENSTH(IDS)=DENSTH(IDS)*1.686E-8
!CC              unit: cm-3
 370  CONTINUE
      WRITE(14,374)
 374  FORMAT(1X,'Suprathermal Density of photoelectrons N(s)')
      WRITE(14,375)(DENSTH(IDS),IDS=1,IMAX)
 375  FORMAT(1X,1P10E12.4)
      WRITE(14,376)
 376  FORMAT(1X,'Average Temperature of Suprathermal p.e.e. E,(s)')
      WRITE(14,375)(EBARTH(IDS),IDS=1,IMAX)

!C........ compute electron impact ion production rates
      print *,' start ion production calculation...'

	  call eimpir(phiup,phidwn)


!C......... close input and output files
	  close(3)
	  close(7)
	  close(9)
	  close(11)
	  close(13)
	  close(17)
	  close(21)
	  close(24)
	  close(12)
	  close(14)
	  close(15)
	  close(16)
	  close(28)

!C......... write phiup(*,ns) to phiout.dat
	  open(24,file='output/phiout.ts',status='old',access='append')
	  write(24,810)(ufx(5,je),je=jmin, jmax)
	  close(24)
      STOP
      END

!**************************END OF PROGRAM**************************
!**************************END OF PROGRAM**************************
!**************************END OF PROGRAM**************************

!C -------------------------------------------------------------
!C    This module contains subroutines for telect.for. The following
!C  subroutines are included:
!C  1. subroutine IMPIT(DEN): The two-stream equation solver.
!C  2. function SIGFE(E,NEXC,NSPEC): airglow cross section.
!C  3. subroutine PARASR(X0,S,R): s - r relation of parabola. the center
!C         (r=0) is located at the focal pt. of the parabola.
!C  3a. subroutine PARA2(X0,S,R): this subroutine computes r(s) of a parabola
!C         where the center (r=0) is located at twice the focal pt, which is
!C         co-centered with the tangential circle of the subsolar pt.
!C  4. subroutine EXTRAP(N,X,Y,XP,YP,IEXP,IRESET): extrapolation.
!C  5. function timer():  cpu timer.
!C -------------------------------------------------------------

      SUBROUTINE IMPIT(DEN)
      parameter (ns=1544)
      DIMENSION DEN(ns),ALPHA(ns),BETA(ns),GAMMA(ns)
      DIMENSION K(ns),L(ns),A(ns),B(ns),C(ns),D(ns),PHI(ns)
      DIMENSION PHIINF(260)
      COMMON/TWO/ALPHA,BETA,GAMMA,IMAX,DELZ,FAC,PHI,PHIINF,J
      REAL K,L
      J1=IMAX-1
      DO 1 JJ=1,J1
      	A(JJ)=PHI(JJ)/DELZ/DELZ+ALPHA(JJ)/2./DELZ
      	B(JJ)=-2.*PHI(JJ)/DELZ/DELZ+BETA(JJ)
      	C(JJ)=PHI(JJ)/DELZ/DELZ-ALPHA(JJ)/2./DELZ
    1 	D(JJ)=GAMMA(JJ)
      K(2)=(D(2)-C(2)*DEN(1))/B(2)
      L(2)=A(2)/B(2)
      DO 2 JJ=3,J1
      	DEM=B(JJ)-C(JJ)*L(JJ-1)
      	K(JJ)=(D(JJ)-C(JJ)*K(JJ-1))/DEM
    2 	L(JJ)=A(JJ)/DEM
      DEN(J1)=(K(J1)-L(J1)*PHIINF(J))/(1.+L(J1)*FAC)
      DEN(IMAX)=DEN(J1)
      J3=IMAX-3
      DO 3 KK=1,J3
      	JK=J1-KK
      	DEN(JK)=K(JK)-L(JK)*DEN(JK+1)
    3 CONTINUE
      RETURN
      END

!C ------------------------------------------------------------
      FUNCTION SIGFE(E,NEXC,NSPEC)
!CCC  THIS FUNCTION RETURNS THE ELECTION IMPACT CROSS-SECTIONS AT
!CCC  ENERGY E FOR SPECIES NSPEC AND AIRGLOW EXCITATION NUMBER NEXC
      COMMON/RABBIT/ NAPP(4),CRAB(60,10,4),ECRAB(60,10,4),NCRAB(10,4)
!C
      NNZ=NCRAB(NEXC,NSPEC)
      IF(E .LE. ECRAB(1,NEXC,NSPEC))GO TO 90
      IF(E .GE. ECRAB(NNZ,NEXC,NSPEC)) GO TO 95
      IJAB=NNZ
   60 CONTINUE
      IF (E.GE.ECRAB(IJAB,NEXC,NSPEC)) GO TO 80
      IJAB=IJAB-2
      IF(IJAB .LE. 1) IJAB=1
      GO TO 60
   80 CONTINUE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      I1=IJAB+2
      I0=IJAB+1
      IM=IJAB
      F1=CRAB(I1,NEXC,NSPEC)
      F0=CRAB(I0,NEXC,NSPEC)
      FM=CRAB(IM,NEXC,NSPEC)
      X1=ECRAB(I1,NEXC,NSPEC)
      X0=ECRAB(I0,NEXC,NSPEC)
      XM=ECRAB(IM,NEXC,NSPEC)
      DET=X0*(X1*X1-X1*X0)+XM*(XM*X1-X1*X1)+X0*(XM*X0-XM*XM)
      AA=X0*(X1*X1-X0*X1)*FM+XM*(XM*X1-X1*X1)*F0+X0*(X0*XM-XM*XM)*F1
      BB=FM*(X0*X0-X1*X1)+F0*(X1*X1-XM*XM)+F1*(XM*XM-X0*X0)
      CC=FM*(X1-X0)+F0*(XM-X1)+F1*(X0-XM)
      AA=AA/DET
      BB=BB/DET
      CC=CC/DET
      SIGFE=AA+BB*E+CC*E*E
      RETURN
   90 CONTINUE
      SIGFE=0.0
      RETURN
   95 CONTINUE
!CCC            PUT IN  FORBIDDEN STATES FIRST IN
!CCC             XSECT LISTS FOR EA SPECIES
      IF(NEXC .EQ. 1) GO TO 97
      CRBM=CRAB(NNZ,NEXC,NSPEC)*ECRAB(NNZ,NEXC,NSPEC)/ALOG(ECRAB(NNZ,NEXC,NSPEC))
      SIGFE=CRBM*ALOG(E)/E
      RETURN
   97 CRBM=ECRAB(NNZ,NEXC,NSPEC)*CRAB(NNZ,NEXC,NSPEC)
      SIGFE=CRBM/E
      RETURN
      END

SUBROUTINE EIMPIR(phiup,phidwn)

!C  This subroutine uses phiup and phidwn calculated in two-stream code
!C  to compute Electron Impact Ion Production Rates.

parameter(ns=1544,ne=260,nmx=10,imx=10,ned=100)
!C               nmx: max. number of neutral species.
!C               imx: max. number of ionization processes for each neutral.
!C               ned: max. number of data points of the data file.
real phiup(ne,ns),phidwn(ne,ns)
real sion(nmx,imx,ne), reii(nmx,imx,ns), ptot(nmx,imx)
real ed(ned),siond(imx,ned),siond1(ned)
real ionxs(nmx,imx), SIGHFAST
integer jneut(nmx), energy
character*5 label(nmx,imx)

!C               common/ionn/ passes GOS parameters to sub. sigion()
!C               common/atmos/passes geo- & e- grid and neutral densities from
!C                            main program.
common/ionn/ak(nmx,imx),akb(nmx,imx),aj(nmx,imx),ajb(nmx,imx),ts(nmx,imx),ta(nmx,imx)
common/ionn/tb(nmx,imx),gams(nmx,imx),gamb(nmx,imx),thi(nmx,imx),nion(nmx)
common/atmos/s(ns),r(ns),ramagl(ns),ee(ne),dele(ne),zspec(nmx,ns), delz

!C $1. read in ionization cross sections
open(23,file='input/eimpir.dat',status='old')
read(23,*)nneut

do in=1,nneut
	read(23,*)
	read(23,*)jneut(in),nion(in),npt,indicator,nskip !,csunit
! Ion spec #, number of ioniz. states NNIN, @ data pts., indicator?, number of lines to skip after, units
	read(23,201)(label(in,j),j=1,nion(in))
	write(6,201)(label(in,j),j=1,nion(in))
 201      format(10(1x,a5))

	do iskip=1,nskip
		read(23,*)
	enddo

	if ( (indicator.eq.1).or.(indicator.eq.0) ) then
		if (indicator.eq.1) then
			do je=1,npt
				read(23,*)ed(je),(siond(ii,je),ii=1,nion(in))
			enddo
		else
			read(23,*)(ed(je),je=1,npt)
			do ii=1,nion(in)
				read(23,*)(siond(ii,je),je=1,npt)
			enddo
		endif

		WRITE(52,*) 'SPECIES= ', IN
!C.........extrapolate
		do ii=1,nion(in)
			WRITE(52,*) 'IONIZATION STATE=', II
			do je=1,npt
				siond1(je)=siond(ii,je)
			enddo
			ires=1
			do ie=1,ne
				call extrap(npt,ed,siond1,ee(ie),sion1,1,ires)
				sion(in,ii,ie)=sion1*csunit
				ires=0
			enddo
		enddo
	else if (indicator.eq.2) then
		do ii=1,nion(in)
			read(23,504)thi(in,ii),ak(in,ii),akb(in,ii),aj(in,ii),ajb(in,ii),ts(in,ii),ta(in,ii),tb(in,ii),gams(in,ii),gamb(in,ii)
			write(16,504)thi(in,ii),ak(in,ii),akb(in,ii),aj(in,ii),ajb(in,ii),ts(in,ii),ta(in,ii),tb(in,ii),gams(in,ii),gamb(in,ii)
!c 504          format(2x,10f7.3)
  504 FORMAT(2X,5F7.3,F7.4,F7.1,F7.2,F7.3,F7.3)
			do ie=1,ne
				sion(in,ii,ie)=sigion(in,ii,ee(ie),0.0,0.0,-1)
			enddo
		enddo
	else
		print *,' indicator=',indicator,' has no meaning. Stop.'
	endif
enddo

      DO I=1,NNEUT
	   write(52,*) 'NGAS = ',I
      	DO ENERGY = 200000,1,-10
      		ETJ = ENERGY * 1.0
	   		DO J = 1,NION(I)
	   			IONXS(I,J)=SIGION(I,J,ETJ,0.,0.,-1)
	   		ENDDO
	   		IF(ETJ .LT. 100.0 .AND. MOD(ETJ,5.) .EQ. 0.0)THEN
	   			WRITE(52,*), ETJ, (IONXS(I,J), J=1,NION(I))
	   		ENDIF
	   		IF(ETJ .GT. 100 .AND. MOD(ETJ,50.) .EQ. 0.0) THEN
	   	   		WRITE(52,*), ETJ, (IONXS(I,J), J=1,NION(I))
      		ENDIF
      		IF (I .EQ. 1)THEN !H2 DISSOCIATION -> fast H ly alpha
      			CALL SIGDIS(ETJ,SIGHFAST)
!C      			PRINT*, ETJ, DISSOC
      			IF(ETJ .LT. 200.0 .AND. MOD(ETJ,5.) .EQ. 0.0)THEN
!c	   				WRITE(53,*), ETJ, SIGHFAST
	   			ENDIF
      		ENDIF
	   	ENDDO
      ENDDO

        close(7)
do i=1,260
!	write(*,*) dele(i)
!	write(*,*) zspec(1,i), zspec(2,i), zspec(3,i), zspec(4,i)
end do
!C $2. compute ion production rates.
        do in=1,nneut
          do ii=1,nion(in)
          	ptot(in,ii) = 0.0
            do is=1,ns
              prion=0.0
              do ie=1,ne
                prion=prion+sion(in,ii,ie)*(phiup(ie,is)+phidwn(ie,is))*dele(ie)
              enddo     ! ie
              jneutral=jneut(in)
              reii(in,ii,is)=zspec(jneutral,is)*prion
			  			ptot(in,ii) = ptot(in, ii) + reii(in,ii,is)* delz
            enddo       ! is altitude bins
          enddo         ! ii ionization states
          write(54,*) 'spec=', in, (ptot(in,ii), ii=1,nion(in))
        enddo           ! in neutral species




!C $3. write output
        write(15,211)((label(in,ii),ii=1,nion(in)),in=1,1)
 211    format(///20X,'ELECTRON IMPACT ION PRODUCTION RATES:',/,2X,'is     s          r        mza ',10(6x,a5))
        do is=1,ns
          write(15,213)is,s(is),r(is)*1e-5,ramagl(is),((reii(in,ii,is),ii=1,nion(in)),in=1,1)
 213      format(i4,1p13e11.3)
        enddo

        write(15,215)((label(in,ii),ii=1,nion(in)),in=2,3)
 215    format(//20X,'ELECTRON IMPACT ION PROD RATES (cont.):',/,2X,'is     s          r        mza ',10(6x,a5))
        do is=1,ns
          write(15,213)is,s(is),r(is)*1e-5,ramagl(is),((reii(in,ii,is),ii=1,nion(in)),in=2,3)
        enddo

       write(15,215)((label(in,ii),ii=1,nion(in)),in=4,4)
       do is=1,ns
         write(15,213)is,s(is),r(is)*1e-5,ramagl(is),((reii(in,ii,is),ii=1,nion(in)),in=4,4)
        enddo

!C   sum up processes for a ion prod. output the sub-ram (is=1) rate to
!C   output file: EIONPR.SUBRAM

        open(25,file='output/eionpr.subram',status='old',access='append')
!*** This was the Titan output
!c        DN2P=reii(1,2,1)
!c        DCH4P=reii(2,1,1)
!c        DCH3P=reii(2,2,1)
!c        DCH2P=reii(2,3,1)
!c        DCHP =reii(2,5,1)
!c        DHP  =reii(2,4,1)+reii(5,1,1)+reii(4,1,1)
!c        DH2P =reii(3,1,1)
!c        DNP  =reii(1,1,1)+reii(6,1,1)
!c        DC2HP=reii(8,2,1)+reii(10,6,1)
!c        DC2H2P=reii(8,1,1)+reii(9,3,1)+reii(10,5,1)
!c        DC2H3P=reii(9,2,1)+reii(10,4,1)
!c        DC2H4P=reii(9,1,1)+reii(10,3,1)
!c        DC2H5P=reii(10,2,1)
!c        DC2H6P=reii(10,1,1)
!c        DCP   =reii(2,6,1)
!c        DHCNP =reii(7,1,1)
!c        DCNP  =reii(7,2,1)
!*** Production rates for Jovian neutrals
		DH2P  = reii(1,1,1)  !H2 ionization
		DHP   = reii(1,2,1) + reii(3,1,1) + reii(4,4,1) !H2->H+, CH4->H+, H->H+
		DHEP  = reii(2,1,1) !He->He+
		DCH4P = reii(4,1,1) !CH4->CH4+
    DCH3P = reii(4,2,1) !CH4->CH3+
    DCH2P = reii(4,3,1) !CH4->CH2+
    DCHP  = reii(4,5,1)  !CH4->CH+
    DCP   = reii(4,6,1)   !CH4->C+
        write(25,231)r(1),DH2P, DHP, DHEP, DCH4P, DCH3P, DCH2P, DCHP, DCP
 231    format(1x,1p14e9.2,/10X,1P10E9.2)

        return
        end

!CG==============================================================
      FUNCTION SIGION(I,ML,E,E1,E2,ICAP)
!CG--------------------------------------------------------------
!CG-----  Electron impact ionization cross-sections     ---------
!CG-----
!CCCCCCCC     ICAP=0 FIRST TIME FOR ENERGY E
!CCCCCCCC      ICAP LT 0 TOTAL X-SECT FOR I,ML,E
!CCCCCCCC      ICAP GT 0 SECOND ENERGIES BETWEEN E1 AND E2

parameter(nmx=10,imx=10,ned=100)
COMMON/IONN/AK(nmx,imx),AKB(nmx,imx),AJ(nmx,imx),AJB(nmx,imx),TS(nmx,imx),TA(nmx,imx),TB(nmx,imx)
COMMON/IONN/GAMS(nmx,imx),GAMB(nmx,imx),THI(nmx,imx),NINN(nmx)
COMMON A,GG,TZ,T12
NAG=NINN(I)
IF(E.LE.THI(I,ML))GO TO 30
QQNN=6.51E-14
QQN=1.0E-16
IF(ICAP.GT.0)GO TO 8
AK1=AK(I,ML)
AKB1=AKB(I,ML)
AJ1=AJ(I,ML)
AJB1=AJB(I,ML)
TS1=TS(I,ML)
TA1=TA(I,ML)
TB1=TB(I,ML)
GAMS1=GAMS(I,ML)
GAMB1=GAMB(I,ML)
IF(NAG.LE.0)GO TO 8
EXM=-AJ1*AJB1+AJ1
IF(E.LE.EXM)GO TO 30
S=QQN*AK1*ALOG(AJB1+E/AJ1)
A=S/(E+AKB1)
TZ=TS1-TA1/(E+TB1)
GG=(GAMS1*E)/(E+GAMB1)
8 CONTINUE
NAG=NINN(I)
IF(NAG.LT.0)GO TO 55
IF(ICAP.LT.0)GO TO 13
TTL=(E-THI(I,ML))/2.0
TTL1=TTL-0.01
IF(E1.GE.TTL1)GO TO 30
IF(E2.GE.TTL)E2=TTL
ABB=(E2-TZ)/GG
ABC=(E1-TZ)/GG
AL2=GG*GG*(ABB*ABB+1.0)
AL1=GG*GG*(ABC*ABC+1.0)
ABB=ATAN(ABB)-ATAN(ABC)
T12=TZ+0.5*GG*(ALOG(AL2)-ALOG(AL1))/ABB
SIGION=A*GG*ABB
RETURN
13 CONTINUE
TM=(E-THI(I,ML))/2.0
IF(E.LE.TM)GO TO 30
SIGION=A*GG*(ATAN((TM-TZ)/GG)+ATAN(TZ/GG))
RETURN
   55 CONTINUE
TTL=(E-THI(I,ML))/2.0
TTL1=TTL-0.01
IF(E1.GE.TTL1)GO TO 30
IF(E2.GE.TTL)E2=TTL
TTHH1=THI(I,ML)+E1
TTHH2=THI(I,ML)+E2
TTHH=THI(I,ML)
AVA=(QQNN*AK1*TTHH**TS1)/(E**AKB1)
IF(AJ1.NE.1.0)GO TO 77
AA1=1.0+TS1-AKB1
AA2=1.0+TS1-AKB1-AJB1
AA1=-AA1
AA2=-AA2
IF(ICAP.LT.0) GO TO 40
ABB=+(TTHH1**AA1)/AA1-(TTHH1**AA2)/(AA2*E**AJB1)
ABBB=+(TTHH2**AA1)/AA1-(TTHH2**AA2)/(AA2*E**AJB1)
SIGION=AVA*(ABBB-ABB)
CA1=AA1+1.0
CA2=AA2+1.0
CBB=(TTHH1**CA1)/CA1-(TTHH1**CA2)/(CA2*E**AJB1)
CBBB=(TTHH2**CA1)/CA1-(TTHH2**CA2)/(CA2*E**AJB1)
T12=(CBBB-CBB)*AVA/SIGION-TTHH
IF(ICAP.GE.0)GO TO 87
   40 CONTINUE
TTL=TTL+TTHH
ABBB=+(TTL**AA1)/AA1-(TTL**AA2)/(AA2*E**AJB1)
ABB=(TTHH**AA1)/AA1-(TTHH**AA2)/(AA2*E**AJB1)
SIGION=AVA*(ABBB-ABB)
87 CONTINUE
RETURN
77 CONTINUE
TTHH1=(TTHH1+TTHH2)/2.0
AA1=2.0+TS1-AKB1
ABB=(1.0-(TTHH1/E)**AJB1)**AJ1
SIGION=(E2-E1)*AVA*ABB/(TTHH1**AA1)
RETURN
30 SIGION=0.0
RETURN
END
!C
!C*********************************************************************
	  SUBROUTINE SIGDIS(E,SIGHFAST)
!C  This function will calculate the H2 dissociation cross section for
!C the process e + H2 -> H(Ly-alpha)
!C The analytical function is taken from Ajello et al. 1995 JGR
!c I found this paper confusing. To get the xs for the fast H ly alpha
!C use 40% of the second state and 60% of the third state and add them
!C together. I could not get this to work for the slow H alpha xs, so I
!C used the data points from the plot found in the paper. Once this xs
!C is calculated, the total dissociation xs must be obtained by adding
!C fast and slow.
!C

REAL COEF(9,3), EIJ(3), SIGD(3), OMEGA(3), XE(3), SIGMA(3)
REAL E, SUM(3), UNITS, SIGHFAST
DATA EIJ/14.67,23.0,30.2/
!C     sLOW
DATA COEF/ 0.80040,0.14727,-0.00592,0.01464,0.40418,-1.1253,1.1253,0.53674,0.14454,-.067184, 0.034857,-0.099768,&
 					 0.50682,-0.92647,0.20852,-0.20852,0.0,0.25704,-.067184,0.034857,-0.099768,0.50682,-0.92647,0.20852,&
					 -0.20852, 0.0, 0.25704/


      UNITS = 8.79E-17*13.6

      DO I = 1,3
        IF(E .LT. EIJ(I))THEN
        	SIGHFAST = 1e-42
        	GOTO 50
        ENDIF
      	XE(I) = E/EIJ(I)
      	SIGMA(I) = COEF(1,I)*(1-1/XE(I))/XE(I)/XE(I)

      	SUM(I) =0.0
      	DO J =2,5
      	  SUM(I)=SUM(I) + COEF(J,I)*(XE(I)-1)*EXP(-J*COEF(9,I)*XE(I))
!C      	PRINT*, COEF(I), SUM
      	ENDDO

      	OMEGA(I) = SIGMA(I) + SUM(I) + COEF(6,I)+ COEF(7,I)/XE(I)+COEF(8,I)*ALOG(XE(I))
!C      PRINT*, 'OMEGA=', OMEGA

      	SIGD(I) = OMEGA(I)/EIJ(I)/XE(I)* UNITS
      	IF(SIGD(I) .LE. 0.0)THEN
      		SIGD(I) = 1E-40
      	ENDIF
      ENDDO

      SIGHFAST = 0.4 * SIGD(2) + 0.6 * SIGD(3)

   50 CONTINUE
      RETURN
      END


!**********************************************************************
!C
!C  3a. subroutine PARA2(X0,S,R): this subroutine computes r(s) of a parabola
!C         where the center (r=0) is located at twice the focal pt, which is
!C         co-centered with the tangential circle of the subsolar pt.

!C	subroutine para2(x0,s,r)
!C	common /angle/theta,rp,alpha		! for checking

!C	x1=x0*0.5
!C theta, rp are angle and radius with respect to the focal point of parabola.
!C	call parasr(x1,s,rp)
!C	call parasr(x0,s,r)
!C	theta=2.0*acos(sqrt(x1/rp))

!C	r=rp**2+x1**2+rp*x0*cos(theta)
!C	r=sqrt(r)

!C alpha is the sza of (s,r) point measured from the center of the planet.
!C	alpha=asin(rp/r*sin(theta))
!C	alpha=alpha*180.0/3.1415926
!C	rflank=sqrt(2.0)*x0
!C	if (r.gt.rflank) alpha=180.0-alpha
!C	write(50,*) alpha, r
!C	return
!C	end


SUBROUTINE PARA2(X0,ARCL, RADIUS)

real X0,LOW_S,HIGH_S,LOW_R,HIGH_R,X,Y,RADIUS,LOW_X,HIGH_X,LOW_Y,HIGH_Y,ARCL
real YKC(INT(1E5)),SKC(INT(1E5)),RKC(INT(1E5)),XKC(INT(1E5))
INTEGER I

   	  COMMON/KC/SKC, RKC, XKC, YKC
	  LOW_S = SKC(1)
	  HIGH_S = SKC(2)
	  LOW_R = RKC(1)
	  HIGH_R = RKC(2)
	  LOW_X = XKC(1)
	  HIGH_X = XKC(2)
	  LOW_Y = YKC(1)
	  HIGH_Y = YKC(2)

	  DO I = 2, 1000
	   IF(ARCL.GT.SKC(I)) THEN
	      LOW_S = SKC(I)
	      HIGH_S = SKC(I+1)
	      LOW_R = RKC(I)
	      HIGH_R = RKC(I+1)
	      LOW_X = XKC(I)
	      HIGH_X = XKC(I+1)
	      LOW_Y = YKC(I)
	      HIGH_Y = YKC(I+1)
	   ENDIF
	  ENDDO

!C       LET LOW_S = LS, LET HIGH_S = HS, LET ARCL = S
!C       LET LOW_R = LR, LET HIGH_R = HR, LET RADIUS = R
!C       LET X BE THE POSITION BETWEEN S AND S+1, WHERE ARCL IS LOCATED
!C       THEN
!C       S = LS +(HS-LS)/((I+1)-I)*(X-I)
!C       S = LS +(HS-LS)*(X-I)
!C       X-I = (S-LS)/(HS-LS)
!C       NOW
!C       R = LR +(HR-LR)/((I+1)-I)*(X-I)
!C       R = LR +(HR-LR)*(X-I)
!C       CONCLUSION:

	  X = LOW_X+(HIGH_X - LOW_X)*(ARCL-LOW_S)/(HIGH_S - LOW_S)
	  Y = LOW_Y+(HIGH_Y - LOW_Y)*(ARCL-LOW_S)/(HIGH_S - LOW_S)
  	  IF (X.LT.0) THEN
	     ALPHA = 180 + ATAN(Y/X)*180.0/3.1415926
	  ELSE
	     ALPHA = ATAN(Y/X)*180.0/3.1415926
	  ENDIF

	  RADIUS = LOW_R +(HIGH_R-LOW_R)*(ARCL-LOW_S)/(HIGH_S-LOW_S)

	  WRITE(50,*) X, Y, RADIUS, ALPHA

	  RETURN
 	  END
!C
!**********************************************************************
!C
!CCCC   Subroutine to calculate the radius of each spatial grid
!CCCC   calls at line ???.
      SUBROUTINE PARASR(X0,S,R)
!C         X0:   subsolar altitude of the magnetic field line.
!C         S:    arc length between the subsolar point and grid.
!C         R:    radius of grid point.
      IF (S .LE. (X0*0.3)) GOTO 404
!C     first guess.
        R=X0+S
        S1=((R-X0)*R)**0.5+X0*ALOG(((R-X0)**0.5+R**0.5)/X0**0.5)
        R1=(S-X0*ALOG(((R-X0)**0.5+R**0.5)/X0**0.5))**2/R+X0
       DO 401 I=1,50
        S2=((R1-X0)*R1)**0.5+X0*ALOG(((R1-X0)**0.5+R1**0.5)/X0**0.5)
        IF (ABS((S2-S)/S) .LT. 1.E-4) GOTO 402
      IF (ABS((S2-S1)/S) .LE. 1.E-4) GOTO 402
        R2=(R1-R)*(S-S1)/(S2-S1)+R
        IF (R2 .LE. X0) R2=X0
        R=R1
        S1=S2
        R1=R2
 401   CONTINUE
 402   CONTINUE
       R=R1
       GOTO 405
 404   R=0.2492*S**2/X0+X0
 405   CONTINUE
      RETURN
      END
!C
!**********************************************************************
	SUBROUTINE EXTRAP(N,X,Y,XP,YP,IEXP,IRESET)
!C ......   Program to interpolate data at xp. X(N),Y(N) are N dimensional
!C ......   data set. Given a x point at xp, the program returns yp for y
!C ......   value.
!C ......   IEXP=0 : linear interpolation.
!C ......   IEXP=1 : exponential-exponential linear interpolation.
!C ......   IRESET=0  : keep old J to start search.
!C ......   IRESET=1  : reset J to 1.

	REAL X(N),Y(N)
	IF (IRESET.EQ.1) J=1
	CALL FINDJ(N,X,XP,J,IFLAG)
	IF (IEXP.EQ.0) THEN
	  CALL LINEAR(N,X,Y,XP,J,IFLAG,YP)
	ELSE IF (IEXP.EQ.1) THEN
	  CALL EXPO(N,X,Y,XP,J,IFLAG,YP)
	ENDIF
	RETURN
	END

	SUBROUTINE LINEAR(N,X,Y,XP,J,IFLAG,YP)
	REAL X(N),Y(N)
	IF (IFLAG .EQ. 1) THEN
	  YP=Y(J)
	ELSE
 		IF (IFLAG .EQ. -1) THEN
		  JL=1
		  JR=2
		ELSE IF (IFLAG .EQ. -2) THEN
		  JL=N-1
		  JR=N
		ELSE IF (IFLAG .EQ. 0) THEN
		  JL=J
		  JR=J+1
		ELSE
		    PRINT *,'IFLAG not returned right'
		ENDIF
		  YP=(Y(JR)-Y(JL))*(XP-X(JL))/(X(JR)-X(JL))+Y(JL)
	ENDIF

	RETURN
	END
!C
!**********************************************************************
	SUBROUTINE EXPO(N,X,Y,XP,J,IFLAG,YP)
	REAL X(N),Y(N),XX1,XX2,YR,YL,XR,XL,YP,XPL
	INTEGER JL,JR
	IF (IFLAG .EQ. 1) THEN
	  YP=Y(J)
	ELSE
 		IF (IFLAG .EQ. -1) THEN
		  JL=1
		  JR=2
		ELSE IF (IFLAG .EQ. -2) THEN
		  JL=N-1
		  JR=N
		ELSE IF (IFLAG .EQ. 0) THEN
		  JL=J
		  JR=J+1
		ELSE
		    PRINT *,'IFLAG not returned right'
		ENDIF

		if ( (y(jr).gt.1.e-30) .and. (y(jl).gt.1.e-30) .and. (x(jr).gt.1.e-30) .and. (x(jl).gt.1.e-30) ) then
		  YR=ALOG(Y(JR))
	   	  YL=ALOG(Y(JL))
		  XR=ALOG(X(JR))
	 	  XL=ALOG(X(JL))
		  XPL=ALOG(XP)
		  YP=(YR-YL)*(XPL-XL)/(XR-XL)+YL
	 	  YP=EXP(YP)
		else
		  YP=0.0
		endif
	ENDIF

	RETURN
	END
!C
!**********************************************************************
	SUBROUTINE FINDJ(NPT,XX,X,J,IFLAG)
	PARAMETER (NMX=500)
	REAL XX(NPT)
	IF (X .LT. XX(1)) GOTO 1
	IF (X .GT. XX(NPT)) GOTO 2
!C	J=1
 13	IF (X .EQ. XX(J)) GOTO 14
	IF ((X .GT. XX(J)) .AND. (X .LT. XX(J+1))) GOTO 15
	J=J+1
	IF (J .LE. NPT) GOTO 13
	PRINT *,' ## INFINITE LOOP IN FINDJ, X=',X
	STOP
 1	IFLAG=-1
	RETURN
 2	IFLAG=-2
	RETURN
 14	IFLAG=1
	RETURN
 15	IFLAG=0
	RETURN
	END
!C
!**********************************************************************
        SUBROUTINE NEUTD(FLAG,QQ,DD,RR,DENS)
!**********************************************************************
!* This subroutine calculates the neutral density of the Enceladus plume
!* based on the model by Saur et al. 2008 equation (1). The angle \Theta
!* is assumed to be zero for this case.
!**********************************************************************
!C
!C Input:
!C		RR --> distance from the surface
!C		Units: cm
!C		Type: Real
!C
!C Output:
!C     	DENS --> neutral density in the plume
!C       Units: cm^-3
!C		Type: Real
!C
!**********************************************************************
!C        REAL UN, TAU1, TAU2, TAU3, LAMBDA1, LAMBDA2
!C        REAL LAMBDA3, DENSTOT, DENS, PI
		REAL RE, N, HD, DENSTOT, DENS
        INTEGER FLAG
        PARAMETER (RE = 2.52E7) !Enceladus Radius
	    PARAMETER (N = 2.5E9) !Density at the center of the plume for the nominal case
	    PARAMETER (HD = 9.48E7) !4X Hill radius in units of cm

!C       FLAG = 1 -> H20
!C       FLAG = 2 -> CO2
!C       FLAG = 3 -> CO

!C        PI = 2.*ASIN(1.0)
!C        UN = 1.0E+5
!C        TAU1 = 1.0E+6*DD*DD
!C        TAU2 = TAU1
!C        TAU3 = TAU1
!C        LAMBDA1 = TAU1*UN
!C        LAMBDA2 = TAU2*UN
!C        LAMBDA3 = TAU3*UN
!C        DENSTOT = QQ/(RR*RR*UN*4*PI)
         DENSTOT =  N* (RE/RR)*(RE/RR)*EXP(-(RR-RE)/HD)
        IF(FLAG.EQ.1) THEN
!C           DENS = 0.85*DENSTOT*EXP(-RR/LAMBDA1)
			DENS =  DENSTOT
	ELSEIF(FLAG.EQ.2) THEN
!C           DENS = 0.08*DENSTOT*EXP(-RR/LAMBDA2)
			DENS = 0.05/0.91 * DENSTOT
        ELSE
!C           DENS = 0.07*DENSTOT*EXP(-RR/LAMBDA3)
			DENS = 0.04/0.91 * DENSTOT
        ENDIF

        RETURN
        END
