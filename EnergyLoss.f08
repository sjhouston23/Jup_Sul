subroutine EnergyLoss(E,ChS,eEnergy,PID,eEnergySS,eAngleSS,eEnergyDS,eAngleDS,dE)
!*******************************************************************************
!* Created by Stephen J. Houston 11.1.18
!*******************************************************************************
!* This routine calculates the energy loss of a precpitating ion.
!* The energy loss obtained by energy loss models based on models
!* calculated by Schultz et. al. 2018 in: Data for secondary electron production
!* from ion precipitation at Jupiter II: Simultaneous and non-simultaneous
!* target and projectile processes in collisions of O^(q+)+H_2 (q=0-8)
!* Table A, pg. 10, Tables 1 and 2, pg. 31 and 32.
!* Also from Schultz et. al. 2016 Ionization of molecular hydrogen and stripping
!* of oxygen atoms and ions... Table 4, pg. 37.
!*******************************************************************************
!*
!* 	Input:
!*		E --> Energy of ions
!*		Type: Real
!*		Units: keV/u
!*
!*		ChS --> Charge state
!*		Type: Integer
!*		Units: None
!*
!*    eEnergy --> Energy of secondary electron(s)
!*    Type: Real*8
!*    Units: eV
!*
!*    PID --> Process ID (1-7,0-4) see CollisionSim.f08
!*    Type: 2-Component Array, Integer
!*    Units: None
!*
!*    eEnergySS --> Energy of single stripping electrons
!*    Type: Real*8
!*    Units: eV
!*
!*    eAngleSS --> Angle of single stripping electrons
!*    Type: Real*8
!*    Units: Degrees
!*
!*    eEnergyDS --> Energy of double stripping electrons
!*    Type: 2-Component array,, real*8
!*    Units: eV
!*
!*    eAngleDS --> angle of double stripping electrons
!*    Type: 2-Component array, real*8
!*    Units: Degrees
!*
!* Returns:
!*    dE --> Delta Energy
!*    Type: Real
!*    Units: eV
!*
!*******************************************************************************
!*
!* Note: For single stripping and double stripping, the electron energies need
!* to be boosted into the projectile frame rather than the target frame,
!* produced by the SDXS. This calculation can be found in the appendix of
!* Schultz et. al. 2018.
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer,intent(in) :: ChS,PID(2)
real*8,intent(in) :: E,eEnergy
real*8,intent(in) :: eEnergySS,eAngleSS,eEnergyDS(2),eAngleDS(2)
real*8,intent(out) :: dE !Change in energy

integer SI,DI,TI,DA,SC,DC,TEX !Target processes
integer SS,DS,SPEX,DPEX !Projectile processes
real*8 IP1,IP2 !First two ionization potentials for H2
real*8 pi,c,sMass !Speed of light and oxygen mass (29.8685e6 keV/c^2)
real*8 auCon !eV to a.u. conversion (1 a.u. = 27.2116 eV)
! real*8 alpha !Fine-structure constant (1/137)
real*8 ThetaP !Theta in the projectile frame

parameter(nChS=17) !Number of charge states from 0-16
parameter(nEnergies=9) !Number of inital energies
parameter(SI=1,DI=2,TI=3,DA=4,SC=5,DC=6,TEX=7) !Target process numbers
parameter(SS=2,DS=3,SPEX=4,DPEX=5) !Projectile process numbers
parameter(c=137.036) !Speed of light in a.u.
parameter(sMass=29.8685e6,auCon=27.2116) !sMass conversion to keV/c^2
parameter(pi=4.0*atan(1.0d0))
parameter(IP1=15.4254,IP2=16.4287) !Ionization potentials for hydrogen
!IP Data from J. Liu et al, Determination of the ionization and dissociation
!energy of the hydrogen molecule, J. Chem. Phys. 130, 174306 (2009).

integer Eng
integer iBin !IonEnergy bin

real*8 Vproj !Projectile velocity (precipitating ion)
real*8 Vz !Electron velocity in projectile frame
real*8 Vsquared !The square of the ejected electron's vel. in projectile frame
real*8 electEnergySS,electEnergyDS1,electEnergyDS2,electEnergyAU1,electEnergyAU2
real*8,dimension(nEnergies) :: IonEnergy !Initial ion energy
real*8,dimension(nChS-1) :: IP !Ionization potentials for sulfur (q=0-15)
real*8,dimension(nEnergies,nChS) :: dETI,dESC,dEDC,dESPEX,dEDPEX !Process dE
!****************************** Data Declaration *******************************
!* Initial ion energy input:
data IonEnergy/1.0,10.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0/
!* Ionization potentials of sulfur (eV):
data IP/10.36001,23.33788,34.86,47.222,72.5945,88.0529,280.954,328.794,379.84,&
        447.74,504.55,564.41,651.96,706.994,3223.7807,3494.1879/
!* New energy loss/gain is in Schultz et al. 2019 Tables 1-11 (eV)
!* Energy loss/gain for Transfer Ionization (Schultz et al. 2019 Table 4) (eV)
data dETI/&
0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !S^ 0+
5.55,  9.13, 40.86, 64.45, 91.42,186.88,297.97,379.29,218.02,& !S^ 1+
5.60,  9.58, 34.31, 57.17, 83.33,165.57,330.81,372.12,305.63,& !S^ 2+
3.57,  9.94, 28.90, 49.74, 72.97,159.18,259.17,361.15,273.79,& !S^ 3+
4.05, 17.99, 24.64, 43.64, 63.83,136.56,271.94,333.08,228.32,& !S^ 4+
2.73, 10.52, 21.49, 38.88, 58.15,123.68,242.42,240.37,392.93,& !S^ 5+
1.98, 10.91, 19.38, 35.61, 53.26,110.33,213.06,254.78,324.96,& !S^ 6+
1.85, 10.31, 17.72, 33.30, 50.65,102.54,200.74,242.69,262.27,& !S^ 7+
1.35, 11.54, 16.93, 31.54, 48.16,101.67,187.74,226.29,215.77,& !S^ 8+
1.50, 11.62, 15.14, 29.38, 46.71,101.81,187.91,200.71,157.65,& !S^ 9+
1.90, 11.50, 14.20, 27.30, 44.40,100.10,175.40,198.10, 95.50,& !S^10+
1.66, 12.16, 13.22, 25.55, 42.47,101.09,185.73,212.21,133.89,& !S^11+
1.70, 11.98, 12.28, 24.03, 40.95,100.36,184.81,255.40,175.67,& !S^12+
1.77, 12.53, 11.62, 22.28, 39.27, 99.77,185.37,206.61,146.47,& !S^13+
1.65, 12.25, 11.83, 20.45, 37.47,100.13,191.71,202.74,147.68,& !S^14+
2.14, 12.32, 11.06, 18.99, 35.79, 99.56,194.50,213.89,146.72,& !S^15+
1.01, 13.25, 11.13, 17.94, 35.52,100.78,193.82,232.76,189.94/  !S^16+

!* Energy loss/gain for Single Capture Charge Transfer (Schultz 2019 Table 4)
data dESC/&
  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,& !S^ 0+
 11.34,  16.24,  38.03,  51.64,  65.26, 119.72, 283.10, 555.41,   0.00,& !S^ 1+
-41.21, -36.31, -14.52,  -0.91,  12.71,  67.17, 230.56, 502.86,1047.48,& !S^ 2+
-47.90, -43.00, -21.22,  -7.60,   6.02,  60.48, 223.86, 496.17,1040.79,& !S^ 3+
-47.90, -43.00,  26.93,  40.55,   6.02,  60.48, 223.86, 496.17,1040.79,& !S^ 4+
-47.90, -43.00, -21.22,  -7.60,   6.02,  60.48, 223.86, 496.17,1040.79,& !S^ 5+
-14.64,  -9.74,  12.04,  36.68,  50.29,  93.74, 209.36, 481.67,1026.29,& !S^ 6+
-10.70,  -5.80,  15.99,  29.60,  43.22,  97.68, 189.77, 462.08,1006.70,& !S^ 7+
-18.86, -13.96,  18.47,  32.08,  45.70,  89.52, 233.31, 505.62, 984.91,& !S^ 8+
-28.11, -23.21,  12.04,  25.66,  39.27,  55.47, 218.85, 416.33, 960.95,& !S^ 9+
-21.82, -16.92,   4.86,  28.50,  42.12,  96.58, 233.31, 390.18, 934.80,& !S^10+
-29.76, -24.86,   9.06,  30.55,  44.16,  98.63, 242.00, 457.15,1001.77,& !S^11+
-24.01, -19.11,  12.04,  25.66,  39.27,  93.74, 233.31, 437.59, 982.21,& !S^12+
-30.96, -26.05,   6.73,  27.88,  41.50,  88.42, 223.86, 416.33, 960.95,& !S^13+
-25.70, -20.80,   0.99,  23.35,  36.96,  91.43, 233.31, 453.37, 997.99,& !S^14+
-31.86, -26.96,   4.86,  25.66,  39.27,  93.74, 225.26, 437.59, 982.21,& !S^15+
-27.03, -22.13,  -0.34,  27.49,  41.10,  89.52, 233.31, 420.72, 965.34/  !S^16+

!* Energy loss/gain for Double Capture Charge Transfer (Schultz 2019 Table 4)
data dEDC/&
   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,& !S^ 0+
   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,& !S^ 1+
 -81.44, -71.64, -28.07,  -0.84,  26.39, 135.32, 462.09,   0.00,   0.00,& !S^ 2+
 -94.83, -85.03, -41.46, -14.23,  13.01, 121.93, 448.70, 993.32,   0.00,& !S^ 3+
 -94.83, -85.03,  54.84,  82.07,  13.01, 121.93, 448.70,   0.00,   0.00,& !S^ 4+
 -47.90, -43.00, -21.22,  -7.60,   6.02,  60.48, 223.86, 496.17,   0.00,& !S^ 5+
-123.82,-114.02, -70.45, -43.22, -15.99,  92.93, 419.70, 964.32,   0.00,& !S^ 6+
 -20.42, -10.62,  32.95,  60.18,  87.41, 196.34, 380.51,   0.00,   0.00,& !S^ 7+
 -36.75, -26.94,  37.91,  65.14,  92.37, 180.01, 467.60,1012.21,1970.80,& !S^ 8+
 -55.25, -45.45,  25.06,  52.29,  79.52, 111.91, 438.68, 833.64,   0.00,& !S^ 9+
 -42.67, -32.87,  10.70,  57.98,  85.22, 194.14, 467.60, 781.34,   0.00,& !S^10+
 -58.55, -48.74,  19.09,  62.07,  89.30, 198.23, 484.98, 915.27,   0.00,& !S^11+
 -47.05, -37.25,  25.06,  52.29,  79.52, 188.45, 467.60, 876.16,   0.00,& !S^12+
 -60.94, -51.13,  14.43,  56.74,  83.97, 177.82, 448.70, 833.64,   0.00,& !S^13+
 -50.42, -40.62,   2.95,  47.67,  74.90, 183.83, 467.60, 907.72,   0.00,& !S^14+
 -62.75, -52.95,  10.70,  52.29,  79.52, 188.45, 451.49, 876.16,   0.00,& !S^15+
 -53.09, -43.28,   0.29,  55.95,  83.18, 180.01, 467.60, 842.41,   0.00/  !S^16+

 !* Energy loss/gain for Single Projectile Excitation (Schultz 2019 Table 4)
 data dESPEX/&
   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,& !S^ 0+
   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,& !S^ 1+
  10.52,  23.85,  62.78,  83.23, 102.06, 165.06,   0.00,   0.00,   0.00,& !S^ 2+
   8.76,  20.56,  61.00,  83.60, 100.82, 170.41,   0.00,   0.00,   0.00,& !S^ 3+
   8.95,  21.98,  66.46,  92.48, 107.74, 190.76,   0.00,   0.00,   0.00,& !S^ 4+
  12.22,  25.13,  76.86, 109.55, 139.15, 227.62, 453.67,   0.00,   0.00,& !S^ 5+
  11.54,  26.35,  89.71, 128.13, 141.17, 227.38, 427.39,   0.00,   0.00,& !S^ 6+
  14.20, 210.66, 330.19, 388.62, 434.49, 290.45, 409.37,   0.00,   0.00,& !S^ 7+
  61.99, 240.74, 359.49, 419.55, 274.19, 379.66, 827.15,   0.00,   0.00,& !S^ 8+
 115.84, 282.49, 405.43, 469.72, 350.97, 646.79, 925.54,   0.00,   0.00,& !S^ 9+
 154.60, 338.90, 468.50, 358.40, 593.20, 733.20, 919.40,   0.00,   0.00,& !S^10+
 177.18, 383.04, 517.46, 386.84, 440.00, 564.96,1031.14,   0.00,   0.00,& !S^11+
 206.61, 247.85, 585.66, 425.45, 481.87, 609.63,1285.96, 971.61,   0.00,& !S^12+
 262.02, 307.05, 428.37, 498.09, 557.95, 704.77, 977.38,1296.54,   0.00,& !S^13+
 289.68, 335.53, 460.06, 533.12, 596.35, 750.41,1030.42,1312.19,   0.00,& !S^14+
2799.27,2959.63,3291.60,3444.64,3575.69,3948.17,4667.92,5619.64,   0.00,& !S^15+
3028.94,3195.59,3539.87,3697.81,3833.95,4221.84,4979.64,5959.04,   0.00/  !S^16+

!* Energy loss for Double Projectile Excitation (Schultz 2019 Table 2)
data dEDPEX/&
  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,& !S^ 0+
 11.34,  16.24,  38.03,  51.64,  65.26, 119.72, 283.10, 555.41,   0.00,& !S^ 1+
-41.21, -36.31, -14.52,  -0.91,  12.71,  67.17, 230.56, 502.86,1047.48,& !S^ 2+
-47.90, -43.00, -21.22,  -7.60,   6.02,  60.48, 223.86, 496.17,1040.79,& !S^ 3+
-47.90, -43.00,  26.93,  40.55,   6.02,  60.48, 223.86, 496.17,1040.79,& !S^ 4+
-47.90, -43.00, -21.22,  -7.60,   6.02,  60.48, 223.86, 496.17,1040.79,& !S^ 5+
-14.64,  -9.74,  12.04,  36.68,  50.29,  93.74, 209.36, 481.67,1026.29,& !S^ 6+
-10.70,  -5.80,  15.99,  29.60,  43.22,  97.68, 189.77, 462.08,1006.70,& !S^ 7+
-18.86, -13.96,  18.47,  32.08,  45.70,  89.52, 233.31, 505.62, 984.91,& !S^ 8+
-28.11, -23.21,  12.04,  25.66,  39.27,  55.47, 218.85, 416.33, 960.95,& !S^ 9+
-21.82, -16.92,   4.86,  28.50,  42.12,  96.58, 233.31, 390.18, 934.80,& !S^10+
-29.76, -24.86,   9.06,  30.55,  44.16,  98.63, 242.00, 457.15,1001.77,& !S^11+
-24.01, -19.11,  12.04,  25.66,  39.27,  93.74, 233.31, 437.59, 982.21,& !S^12+
-30.96, -26.05,   6.73,  27.88,  41.50,  88.42, 223.86, 416.33, 960.95,& !S^13+
-25.70, -20.80,   0.99,  23.35,  36.96,  91.43, 233.31, 453.37, 997.99,& !S^14+
-31.86, -26.96,   4.86,  25.66,  39.27,  93.74, 225.26, 437.59, 982.21,& !S^15+
-27.03, -22.13,  -0.34,  27.49,  41.10,  89.52, 233.31, 420.72, 965.34/  !S^16+
!******************************** Main Program *********************************
!* Initialize:
iBin=0;k=0;f=0.0
do Eng=nEnergies,1,-1 !Loop through bins to get the correct ion energy bin
  if(E.ge.real(IonEnergy(Eng-1)+(IonEnergy(Eng)-IonEnergy(Eng-1))/2.0))then
    iBin = Eng !Get the ion energy bin number
    k = 1 !Used to get an interpolation function
    goto 1000
  elseif(E.ge.real(IonEnergy(Eng-1)))then
    iBin = Eng-1 !Get the ion energy bin number
    k = 2 !Used to get an interpolation function
    goto 1000
  endif
end do
1000 continue
if (iBin.eq.0) write(*,*) 'EnergyLoss.f08: Error! iBin=0'
!* Want to use f to somewhat interpolate the cross-section for ion energies that
!* lie between energy bins.
if (k.eq.1) f=(E-IonEnergy(iBin-1))/(IonEnergy(iBin)-IonEnergy(iBin-1))
if (k.eq.2) f=(E-IonEnergy(iBin))/(IonEnergy(iBin+1)-IonEnergy(iBin))

!* Initialize:
dE=0.0;electEnergyAU1=0.0;Vproj=0.0;Vz=0.0;Vsquared=0.0;electEnergySS=0.0
ThetaP=0.0;electEnergyDS1=0.0;electEnergyDS2=0.0;electEnergyAU2=0.0

!* Find the dE for the given process based on Schultz et al. 2018 Table 1
!* Go through the "target" processes
if(PID(1).eq.SI)then !Single Ionization
  dE=IP1+eEnergy
elseif(PID(1).eq.DI)then !Double Ionization
  dE=IP1+IP2+eEnergy !Both electron energies have been added together
elseif(PID(1).eq.TI)then !Transfer Ionization
  if(f.ge.0.5)dE=IP1+eEnergy+(f*dETI(iBin,ChS)+(1-f)*dETI(iBin-1,ChS))
  if(f.lt.0.5)dE=IP1+eEnergy+(f*dETI(iBin+1,ChS)+(1-f)*dETI(iBin,ChS))
elseif(PID(1).eq.DA)then !Double-Capture Autoionization
  dE=eEnergy
elseif(PID(1).eq.SC)then !Single Capture
  if(f.ge.0.5)dE=(f*dESC(iBin,ChS)+(1-f)*dESC(iBin-1,ChS))
  if(f.lt.0.5)dE=(f*dESC(iBin+1,ChS)+(1-f)*dESC(iBin,ChS))
elseif(PID(1).eq.DC)then !Double Capture
  if(f.ge.0.5)dE=(f*dEDC(iBin,ChS)+(1-f)*dEDC(iBin-1,ChS))
  if(f.lt.0.5)dE=(f*dEDC(iBin+1,ChS)+(1-f)*dEDC(iBin,ChS))
elseif(PID(1).eq.TEX)then !Target Excitation
  dE=7.7
end if

!* Go through the "projectile" processes
!* To review SS and DS frame of reference transformations refer to Schutlz,
!* et al. 2018, Appendix B
if(PID(2).eq.SS)then !Single Stripping
  electEnergyAU1=eEnergySS/auCon !Convert to a.u.
  Vproj=sqrt(2.0*E*16.0/sMass)*c !In a.u.
  Vz=sqrt(2.0*electEnergyAU1)*cos(eAngleSS*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU1)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electEnergySS=(0.5)*Vsquared*auCon !Convert back to eV
  ThetaP=acos(Vz/sqrt(Vsquared))*180/pi
  dE=dE+IP(ChS)+electEnergySS
elseif(PID(2).eq.DS)then !Double Stripping
  electEnergyAU1=eEnergyDS(1)/auCon !Convert to a.u.
  Vproj=sqrt(2.0*E*16.0/sMass)*c !In a.u.
  Vz=sqrt(2.0*electEnergyAU1)*cos(eAngleDS(1)*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU1)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electEnergyDS1=(0.5)*Vsquared*auCon !Convert back to eV
  electEnergyAU2=eEnergyDS(2)/auCon !Convert to a.u.
  Vproj=sqrt(2.0*E*16.0/sMass)*c !In a.u.
  Vz=sqrt(2.0*electEnergyAU2)*cos(eAngleDS(2)*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU2)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electenergyDS2=(0.5)*Vsquared*auCon !Convert back to eV
  dE=dE+IP(ChS)+IP(ChS+1)+electenergyDS1+electenergyDS2
elseif(PID(2).eq.SPEX)then !Single Projectile Excitation
  if(f.ge.0.5)dE=dE+(f*dESPEX(iBin,ChS)+(1-f)*dESPEX(iBin-1,ChS))
  if(f.lt.0.5)dE=dE+(f*dESPEX(iBin+1,ChS)+(1-f)*dESPEX(iBin,ChS))
elseif(PID(2).eq.DPEX)then !Double Projectile Excitation
  if(f.ge.0.5)dE=dE+(f*dEDPEX(iBin,ChS)+(1-f)*dEDPEX(iBin-1,ChS))
  if(f.lt.0.5)dE=dE+(f*dEDPEX(iBin+1,ChS)+(1-f)*dEDPEX(iBin,ChS))
end if


end subroutine
