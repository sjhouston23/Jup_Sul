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
real*8,dimension(nEnergies,nChS) :: dETI,dESC,dEDC,dESPEX,dEDPEX,dETEX !Process
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

 data dETEX/&
6.90,8.00,7.60,7.40,7.30,7.20,7.10,7.10,7.10,& !S^ 0+
6.20,7.70,7.50,7.40,7.30,7.10,6.90,6.90,6.80,& !S^ 1+
5.60,7.10,7.40,7.30,7.20,7.10,6.80,6.60,6.50,& !S^ 2+
5.50,6.70,7.30,7.20,7.20,7.20,6.70,6.50,6.20,& !S^ 3+
5.00,6.40,7.20,7.20,7.30,7.20,6.90,6.50,6.20,& !S^ 4+
5.00,5.80,7.10,7.20,7.20,7.30,7.00,6.70,6.30,& !S^ 5+
5.00,5.80,6.90,7.10,7.20,7.30,7.10,6.80,6.40,& !S^ 6+
4.80,5.70,6.80,7.10,7.20,7.30,7.20,6.90,6.50,& !S^ 7+
4.90,5.50,6.70,7.00,7.20,7.40,7.20,6.90,6.50,& !S^ 8+
4.90,5.40,6.60,6.80,7.10,7.40,7.30,7.00,6.70,& !S^ 9+
4.70,5.20,6.60,6.90,7.10,7.30,7.30,7.10,6.80,& !S^10+
4.70,5.30,6.50,6.80,7.10,7.30,7.40,7.10,6.80,& !S^11+
4.90,5.10,6.30,6.80,7.00,7.30,7.40,7.10,6.80,& !S^12+
4.70,5.30,6.40,6.60,7.00,7.30,7.40,7.20,6.90,& !S^13+
4.90,5.00,6.30,6.60,6.90,7.30,7.40,7.20,6.90,& !S^14+
4.60,4.90,6.10,6.70,6.90,7.30,7.40,7.20,6.90,& !S^15+
4.60,5.00,6.30,6.60,6.80,7.30,7.40,7.20,6.90/  !S^16+

 !* Energy loss/gain for Single Projectile Excitation (Schultz 2019 Table 4)
 data dESPEX/&
  9.75,   9.81,   9.79,   9.78,   9.79,   9.78,   9.78,   9.77,   9.78,& !S^ 0+
 20.45,  21.02,  20.98,  20.97,  20.95,  20.94,  20.95,  20.92,  20.93,& !S^ 1+
 27.63,  29.24,  29.47,  29.38,  29.39,  29.30,  29.23,  29.17,  29.17,& !S^ 2+
 35.23,  36.09,  37.38,  37.23,  37.16,  36.97,  36.83,  36.78,  36.75,& !S^ 3+
 53.77,  54.06,  56.84,  56.93,  56.90,  56.68,  56.27,  56.03,  55.77,& !S^ 4+
 57.39,  58.85,  63.93,  64.69,  64.75,  64.22,  63.64,  63.51,  63.15,& !S^ 5+
204.49, 191.54, 201.87, 206.15, 208.64, 215.28, 217.94, 216.39, 215.07,& !S^ 6+
224.37, 206.76, 220.74, 225.86, 230.50, 239.33, 244.42, 242.60, 241.58,& !S^ 7+
241.58, 221.46, 236.81, 244.08, 249.49, 261.81, 270.94, 268.42, 266.55,& !S^ 8+
250.30, 252.98, 265.26, 274.22, 281.05, 295.83, 310.77, 309.47, 307.45,& !S^ 9+
228.76, 258.44, 280.93, 288.79, 296.52, 314.14, 334.86, 333.95, 330.23,& !S^10+
228.00, 272.18, 292.77, 301.14, 310.34, 330.15, 356.80, 356.27, 352.88,& !S^11+
291.24, 314.02, 330.62, 337.03, 347.64, 371.02, 403.97, 406.32, 404.43,& !S^12+
299.07, 315.17, 326.68, 336.48, 347.10, 372.13, 407.16, 416.50, 409.43,& !S^13+
2438.0,2121.46,1927.67,1889.13,1935.80,2065.67,2113.66,2182.27,2262.92,& !S^14+
2630.0,2278.54,2061.61,1964.24,2057.49,2114.38,2202.98,2296.46,2365.67,& !S^15+
  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00/  !S^16+

!* Energy loss for Double Projectile Excitation (Schultz 2019 Table 2)
data dEDPEX/&
  38.27,  42.93,  44.19,  44.84,  45.66,  46.05,  45.50,  46.87,  43.55,& !S^ 0+
  58.58,  59.22,  62.05,  62.37,  62.03,  61.36,  62.94,  62.12,  62.13,& !S^ 1+
  70.71,  70.82,  71.33,  71.18,  71.10,  71.18,  71.19,  70.92,  70.56,& !S^ 2+
 102.92, 104.56, 111.87, 112.65, 113.39, 113.25, 112.71, 112.70, 112.05,& !S^ 3+
 124.25, 125.00, 128.32, 129.71, 129.68, 129.00, 128.66, 128.44, 128.30,& !S^ 4+
 321.07, 323.21, 332.73, 338.68, 352.57, 390.98, 484.74, 507.05, 499.67,& !S^ 5+
 507.28, 509.42, 518.94, 521.60, 541.29, 567.71, 619.70, 634.14, 645.05,& !S^ 6+
 539.53, 541.67, 551.19, 569.14, 590.16, 610.96, 666.58, 678.00, 673.28,& !S^ 7+
 596.26, 598.40, 607.92, 615.64, 633.18, 661.83, 689.80, 706.83, 704.81,& !S^ 8+
 630.79, 632.93, 642.45, 667.17, 676.76, 691.18, 726.36, 735.05, 727.62,& !S^ 9+
 627.48, 629.62, 639.14, 648.07, 665.55, 706.20, 727.74, 725.80, 746.74,& !S^10+
 722.54, 724.68, 734.20, 740.15, 758.01, 812.17, 853.55, 857.33, 869.33,& !S^11+
 707.76, 709.90, 719.42, 725.37, 774.72, 809.16, 850.28, 877.29, 859.97,& !S^12+
2694.52,2696.66,2706.18,2712.13,2718.08,2741.88,2813.28,2921.60,3127.02,& !S^13+
   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,& !S^14+
   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,& !S^15+
   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00/  !S^16+
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
if (iBin.eq.0) write(206,*) 'EnergyLoss.f08: Error! iBin=0'
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
  if(f.ge.0.5)dE=(f*dETEX(iBin,ChS)+(1-f)*dETEX(iBin-1,ChS))
  if(f.lt.0.5)dE=(f*dETEX(iBin+1,ChS)+(1-f)*dETEX(iBin,ChS))
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
