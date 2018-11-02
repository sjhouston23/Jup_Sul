module formatting
implicit none

!****************************** Header Variables *******************************
character(len=*),parameter :: &
&H01="(' Charge state equilibrium fractions.',/,' Energy [keV/u]',8x,'S',8x,&
 &'S+',7x,'S2+',7x,'S3+',7x,'S4+',7x,'S5+',7x,'S6+',7x,'S7+',7x,'S8+',7x,'S9+',&
 &6x,'S10+',6x,'S11+',6x,'S12+',6x,'S13+',6x,'S14+',6x,'S15+',6x,'S16+')",&
&H02="(' H^+ production rate [cm^-3 s^-1] from sulfur ion precipitation.',/,&
 &' Input of 1 ion/cm^2/s.',/,&
 &' Alt [km] Production Rate')",&
&H03="(' H_2^+ production rate [cm^-3 s^-1] from sulfur ion precipitation.',/,&
 &' Input of 1 ion/cm^2/s.',/,&
 &' Alt [km] Production Rate')",&
&H04="(' H_2 excitation rate [cm^-3 s^-1] from sulfur ion precipitation.',/,&
 &' Input of 1 ion/cm^2/s.',/,&
 &' H_2 is only excited from a single NSIM process, TEX.',/,&
 &' Alt [km] Production Rate')"
!************************** Data Formatting Variables **************************
!* Notes:
!*   1x to create a space before the data
!*   F7.2 is used for atmosphere ==> "3000.00" or " -88.00"
!*   F8.2 is used for energy ==> "25000.00" or "  150.00"
!*   ES9.2 is used for large/small data values ==> "-1.34E+04" or " 2.31E-11"
!*   Electron two-stream formatting - 912 (F2Str)
!*******************************************************************************
character(len=*),parameter :: &
&F1="(1x,A10,43x,F7.2)",&
&F2="(1x,A26,I4,':',I2,':',F5.2)",&
&F3="(1x,60('*'))",&
&F4="(1x,85('-'))",&
&F01="(1x,F8.2,5x,17(2x,ES8.2))",&
&F02="(1x,F7.2,5x,ES8.2)",&
&F03="(32x,A8,5x,A8,3x,A8,5x,A8,6x,'Sum')",&
&F04="(1x,A4,4x,6(I12,1x))",&
&F05="(1x,A4,4x,6(F12.2,1x))"
!****************************** Processes Header *******************************
character(len=8) TargColl(7),ProjColl(4)
data TargColl/'SI      ','DI      ','TI      ','DCAI    ','SC      ',&
              'DC      ','TEX     '/
data ProjColl/'SS      ','DS      ','SPEX    ','DPEX    '/

end module
