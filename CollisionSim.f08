subroutine CollisionSim(E,SIMxs,SIMxs_Total,ChS,excite,elect,disso,PID)
!*******************************************************************************
!* Created by Stephen J. Houston 10.22.18
!*******************************************************************************
!* This subroutine reads in each individual cross-section and the total cross-
!* section for all the processes for each energy and charge state for 35
!* different collision types (see processes below). Then uses a random
!* number from 0-1 to and compare it to the probability for each
!* collision process (individual cross-section divided by total cross-section).
!* This probability will be calcultated from the cross-sections for the
!* different collision processes. The probability that is larger than the random
!* number will be the selected outcome process of the collision.
!*******************************************************************************
!*
!* 	Input:
!*		E --> Energy of the ion
!*		 Type: Integer
!*		 Units: keV/u
!*
!*		SIMxs --> SIM cross-sections
!*		 Type: Real*8 Matrix
!*		 Units: cm^-2
!*
!*		SIMxs_Total --> Sum of SIM cross-sections vs. energy and charge state
!*		 Type: Real*8 Matrix
!*		 Units: cm^-2
!*
!*		ChS --> Charge state
!*		 Type: Integer
!*		 Units: None
!*
!*   Returns:
!*    ChS --> Replaced with new charge state
!*		 Type: Integer
!*		 Units: None
!*
!*    excite --> Number of ion excitations
!*     Type: Integer
!*     Units: None
!*
!*    elect --> Number of electrons ejected
!*     Type: Integer
!*     Units: None
!*
!*    disso --> Number indicating whether dissociation is possible (0-2)
!*     Type: Integer
!*     Units: None
!*
!*    PID --> Collision Type (1-7,0-4)
!*     Type: 2-Component Array, Integer
!*     Units: None
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer,intent(in) :: E!,PID(2)
integer,intent(inout) :: ChS
integer,intent(out) :: excite,elect,disso,PID(2)

integer SI,DI,TI,DCAI,SC,DC,TEX !Target processes
integer SS,DS,SPEX,DPEX !Projectile processes

parameter(nProjProc=4) !Number of projectile processes
parameter(nTargProc=7) !Number of target processes
parameter(nInterpEnergies=2000) !Number of interpolated energies
parameter(nChS=17) !Number of charge states from 0-16
parameter(SI=1,DI=2,TI=3,DCAI=4,SC=5,DC=6,TEX=7) !Target process numbers
parameter(SS=2,DS=3,SPEX=4,DPEX=5) !Projectile process numbers
parameter(k1=0,k2=0,lux=3) !lux set to 3 for optimal randomness and timeliness

real*8,dimension(nChS,nInterpEnergies),intent(in) :: SIMxs_Total
real*8,dimension(nTargProc,1+nProjProc,nChS,nInterpEnergies),intent(in) :: SIMxs

integer tProc,pProc

real*8 sumProb
real*8,dimension(nTargProc,1+nProjProc) :: Prob
!* Prob has an additonal projectile process which is no projectile process
!* i.e. SI, SI+SS, SI+DS, SI+SPEX, SI+DPEX (respectively)

real ranVecB(1) !Number from RNG
!**************************** Initialize Variables *****************************
sumProb=0.0;Prob=0.0;excite=0;elect=0;diss=0;PID=0
!******************************** Main Program *********************************
!******************** Collision-Type Probability Calculation *******************
!* Calculate the transition probabilities by taking the collision process XS and
!* dividing it by the total cross-section.
!*******************************************************************************
! goto 2000
do tProc=1,nTargProc !Loop through every target process
  do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
    Prob(tProc,pProc)=SIMxs(tProc,pProc,ChS,E)/SIMxs_Total(ChS,E)
  end do !End loop through every projectile process plus 1
end do !End loop through every target process
if(sum(Prob).ge.1.00001.or.sum(Prob).le.0.9999)then !Warning if probability is bad
  write(206,*) "CollisionSim.f08: WARNING: Normalized collision probability is &
             &not close enough to 1. The value is: ", sum(Prob)
  write(206,*) "tempQ: ", ChS, "Energy: ", E
end if
!*******************************************************************************
!************************* Collision-Type Determination ************************
!*******************************************************************************
!* Now that the collision probability is determined, a simple monte carlo
!* technique is used to determine the collision type. After collision-type is
!* determined the final charge state and collision-type are recorded as
!* tempQ and process, respectively.
!*
!* Algorithm:
!*  Step 1) Generate a random number on the interval [0,1]
!*  Step 2) Parce a sub-interval to cover the probability of a single
!*          collision type; i.e. for Single Ionization from 0->P(SI) [P(SI)<1]
!*  Step 3) Ask if the randomly generated number is in the interval [0,P(SI)],
!*          if so write out the final charge state and process
!*  Step 4) if not, parce the next sub-interval of [0,1] from [P(SI),P(DI)]
!*  Step 5) Repeat step 3.
!*  Step 6) Repreat step 4&5 through all 35 collision types until collision is
!*          determined.
!*
!* PID(1) Target Processes           Abbreviation  Charge Change  Sulfur Excite
!*   1    Single Ionization              (SI)            0              0
!*   2    Double Ionization              (DI)            0              0
!*   3    Transfer Ionization            (TI)           -1              0
!*   4    Double-Capture Autoionization  (DCAI)         -1            2->0
!*   5    Single Capture                 (SC)           -1              1
!*   6    Double Capture                 (DC)           -2              2
!*   7    Target Excitation              (TEX)           0              0
!* PID(2) Projectile Processes
!*   1    No projectile process occurred
!*   2    Single Stripping               (SS)           +1              0
!*   3    Double Stripping               (DS)           +2              0
!*   4    Single Projectile Excitation   (SPEX)          0              1
!*   5    Double Projectile Excitation   (DPEX)          0              2
!*
!* Some processes will dissociate H2 always, sometimes, or never. This is
!* followed with the "disso" variable.
!* disso = 0, never dissociates
!* disso = 1, 10% chance of dissociation (determined in the main program)
!* disso = 2, always dissociates
!*
!* PID is used as a processes identification
!*
!*******************************************************************************
call ranlux(ranVecB,1) !Only need 1 random number every time there's a collision
if(ranVecB(1).gt.0.99999)ranVecB(1)=ranVecB(1)-0.00001
do tProc=1,nTargProc !Loop through every target process
  do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
    sumProb=sumProb+Prob(tProc,pProc) !Keep adding the probability
    if(ranVecB(1).le.sumProb)goto 1000 !If it falls within range, get out
  end do !End loop through every projectile process plus 1
end do !End loop through every target process
!If we get into this portion, it means that ranVecB was .gt. sumProb
write(206,*)'CollisionSim.f08: ERROR: Random number greater than normalized &
&probability:', ranVecB(1), sumProb
STOP 'CollisionSim.f08: Stopping program...'
1000 continue
PID(1)=tProc
PID(2)=pProc !1 is nothing, 2-5 are SS, DS, SPEX, DPEX
! 2000 continue
! tProc=PID(1)
! pProc=PID(2)
if(tProc.eq.SI)then !Single Ionization
  ChS=ChS !Charge state stays the same
  excite=0 !Ion isn't excited
  elect=1 !One electron ejected
  disso=1 !10% chance of dissociation
elseif(tProc.eq.DI)then !Double Ionization
  ChS=ChS
  excite=0
  elect=2
  disso=2
elseif(tProc.eq.TI)then !Transfer Ionization
  ChS=ChS-1
  excite=0
  elect=1
  disso=2
elseif(tProc.eq.DCAI)then !Double Capture Autoionization
  ChS=ChS-1
  excite=0
  elect=1
  disso=2
elseif(tProc.eq.SC)then !Single Capture
  ChS=ChS-1
  excite=1
  elect=0
  disso=1
elseif(tProc.eq.DC)then !Double Capture
  ChS=ChS-2
  excite=2
  elect=0
  disso=2
elseif(tProc.eq.TEX)then !Target Excitation
  ChS=ChS
  excite=0
  elect=0
  disso=0
end if !End of target processes
if(pProc.eq.SS)then !Single Stripping
  ChS=ChS+1 !Strips a single electron - charge state changes
  excite=excite !No additional excitation
  elect=elect+1 !Ejects an additional electron
elseif(pProc.eq.DS)then !Double Stripping
  ChS=ChS+2
  excite=excite
  elect=elect+2
elseif(pProc.eq.SPEX)then !Single Projectile Excitation
  ChS=ChS
  excite=excite+1
  elect=elect
elseif(pProc.eq.DPEX)then !Double Projectile Excitation
  ChS=ChS
  excite=excite+2
  elect=elect
end if !End of projectile processes
end subroutine
