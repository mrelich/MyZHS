C ____________________________________________________________
      FUNCTION RPAIR(E)
c
c      Parametrization of total cross section including
c      LPM in moon.
c      Fit made with MINUIT (CERN's subroutine)
C ____________________________________________________________
C   TOTAL CROSS SECTION FOR PAIR PRODUCTION IN WATER
C   WRITTEN BY E. ZAS 31 JAN 1991 INCLUDES LPM EFFECT IN Moon
C   Rewritten by J.Alvarez-Munhiz, 17 Oct 2000 
C ____________________________________________________________
c      COMMON/BLMASS/e_MASS,e2_MAS
      COMMON / PAIR_XSECTION / PP_E(400),PP_X(400)

      RPAIR=divdif(PP_X,PP_E,400,E,2)
      return
      end

C ____________________________________________________________
      SUBROUTINE RPAI4(E,E1,AP)
C ____________________________________________________________
      COMMON/BLMASS/e_MASS,e2_MAS
c      R=s_rndm(0)*0.5            !  1/2 so that U<=1/2
c      E1=e_MASS+(E/2.-e_MASS)*R

      R=s_rndm(0)            !  1/2 so that U<=1/2
      E1=(e_MASS+AP)*(1.-R)+(E-(e_MASS+AP))*R
      RETURN
      END

C ____________________________________________________________
C ____________________________________________________________

      SUBROUTINE RPAMFC(E1,U,NO,ALF)
c      SUBROUTINE RPAMFC(Z,E1,U,NO,ALF)
C ____________________________________________________________
C ************************************************************
C **********PAIR PRODUCTION.SCREENING PARAMETER LE 2.*********
C ************************************************************
C Z IS THE ATOMIC NUMBER OF THE MEDIUM
C E1 IS THE PHOTON ENERGY
C U IS THE ENERGY OF ONE OF THE PRODUCED ELECTRONS
C NO IS THE NUMBER OF SAMPLES OF RKSI BEFORE REJECTION
C FCZ IS THE COULOMB CORRECTION TERM
C ALF IS THE SUM OF COEFFICIENS ALFAI
C ************************************************************
      COMMON/BLMASS/e_MASS,e2_MAS
c      COMMON/FCO/FCZ
c      COMMON/VLZ/ALZ
      GAMA(A_A,B,C_C)=51.1/(A_A*B*(1.-B)*C_C**.3333)
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
     
C **************GAMMA IS THE SCREENING PARAMETER**************
      FF1(B1)=20.838-2.647*B1-1.496*B1**2+2.417*B1**3-.860*B1**4
      FF2(B2)=20.170-1.176*B2-2.382*B2**2+2.097*B2**3-.486*B2**4
      F12(B3)=21.12-4.184*ALOG(B3+.952)
      ALFA1=.6666666-1./(36.*(ALOG(183./Z**.3333333)-FCZ))
      ALFA2=(1.333333+.1111111/(ALOG(183./Z**.3333333)-FCZ))/12.
      ALF=ALFA1+ALFA2
      NO=0
  100 CONTINUE
      G1=s_rndm(0)
      IF(G1-ALFA1/(ALFA1+ALFA2)) 1,2,2
    1 F1=1
C ---------------------------------------------------------------
C FIRST TERM
C ---------------------------------------------------------------
    6 G2=s_rndm(0)*0.5            !  1/2 so that U<=1/2
      NO=NO+1
      G3=s_rndm(0)
      RKSI=G2
      if (rksi.eq.0..or.rksi.eq.1.) then 
         NO=NO-1
         goto 6
      end if
      GG=1.36*GAMA(E1,RKSI,Z)
      IF(GG-1) 3,3,4
C ---------------------------------------------------------------
C REJECTION FOLLOWS
C ---------------------------------------------------------------
    3 IF(G3-(9.*FF1(GG)+3.*FF2(GG)-16.*ALZ-48.*FCZ)/
     1       (9.*FF1(0.)+3.*FF2(0.)-16.*ALZ-48.*FCZ)) 5,100,100
    4 IF(G3-(9.*F12(GG)+3.*F12(GG)-16.*ALZ-48.*FCZ)/
     1       (9.*F12(0.)+3.*F12(0.)-16.*ALZ-48.*FCZ)) 5,100,100
C ---------------------------------------------------------------
C SECOND TERM
C ---------------------------------------------------------------
    2 G2=s_rndm(0)*0.5            !  1/2 so that U<=1/2
      NO=NO+1
      G3=s_rndm(0)
      RKSI=-((0.125-G2/4.)**0.33333333)+0.5
      if (rksi.eq.0..or.rksi.eq.1.) then 
         NO=NO-1
         goto 2
      end if
      GG=1.36*GAMA(E1,RKSI,Z)
      IF(GG-1) 7,7,8
C ---------------------------------------------------------------
C REJECTION FOLLOWS
C ---------------------------------------------------------------
    7 IF(G3-(9.*FF1(GG)-3.*FF2(GG)-8.*ALZ-24.*FCZ)/(9.*FF1(0.)-3.*
     1FF2(0.)-8.*ALZ-24.*FCZ)) 5,100,100
    8 IF(G3-(9.*F12(GG)-3.*F12(GG)-8.*ALZ-24.*FCZ)/(9.*F12(0.)-3.*
     1F12(0.)-8.*ALZ-24.*FCZ)) 5,100,100
C ---------------------------------------------------------------
    5 U=RKSI
      RETURN
      END

C _____________________________________________________________________
C *********************************************************************
      SUBROUTINE RPAIME(E1,U,NO,ALF)
c       SUBROUTINE RPAIME(Z,E1,U,NO,ALF)      
C *********************************************************************
C _____________________________________________________________________
C *********************************************************************
C *************PAIR PRODUCTION.SCREENING PARAMETER LE 2.***************
C *********************************************************************
C Z IS THE ATOMIC NUMBER OF THE MEDIUM
C E1 IS THE PHOTON ENERGY
C U IS THE ENERGY OF ONE OF THE PRODUCET ELECTRONS
C NO IS THE NUMBER OF SAMPLES OF RKSI BEFORE REJECTION
C ALF IS THE SUM OF COEFFICIENS ALFAI
C *********************************************************************
c      COMMON/VLZ/ALZ
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
      GAMA(A_A,B,C_C)=51.1/(A_A*B*(1.-B)*C_C**.3333)
C ---------------------------------------------------------------------
C GAMMA IS THE SCREENING PARAMETER
C ---------------------------------------------------------------------
      FF1(B1)=20.838-2.647*B1-1.496*B1**2+2.417*B1**3-.860*B1**4
      FF2(B2)=20.170-1.176*B2-2.382*B2**2+2.097*B2**3-.486*B2**4
      F12(B3)=21.12-4.184*ALOG(B3+.952)
      ALFA1=.66666-1./(36.*ALOG(183./Z**0.3333))
      ALFA2=      (1.3333+1/(9*ALOG(183/Z** 0.3333 )))/12.
      ALF=ALFA1+ALFA2
      NO=0
  100 CONTINUE
      G1=s_rndm(0)
      IF(G1-ALFA1/(ALFA1+ALFA2)) 1,2,2
    1 F1=1
C ----------------------------------------------------------------
C FIRST TERM
C ----------------------------------------------------------------
    6 G2=s_rndm(0)*0.5            !  1/2 so that U<=1/2
      NO=NO+1
      G3=s_rndm(0)
      RKSI=G2
      if (rksi.eq.0..or.rksi.eq.1.) then 
         NO=NO-1
         goto 6
      end if
      GG=1.36*GAMA(E1,RKSI,Z)
      IF(GG-1) 3,3,4
C ----------------------------------------------------------------
C REJECTION FOLLOWS
C ----------------------------------------------------------------
    3 IF(G3-(9*FF1(GG)+3*FF2(GG)-16*ALZ)/(9*FF1(0.)+3*FF2(0.)-16*ALZ))
     #   5,100,100
    4 IF(G3-(9*F12(GG)+3*F12(GG)-16*ALZ)/(9*F12(0.)+3*F12(0.)-16*ALZ))
     #   5,100,100                                               
C -----------------------------------------------------------------
C SECOND TERM
C -----------------------------------------------------------------
    2 G2=s_rndm(0)*0.5            !  1/2 so that U<=1/2
      NO=NO+1
      G3=s_rndm(0)
      RKSI=-((0.125-G2/4.)**0.3333333)+0.5
      if (rksi.eq.0..or.rksi.eq.1.) then 
         NO=NO-1
         goto 2
      end if
      GG=1.36*GAMA(E1,RKSI,Z)
      IF(GG-1) 7,7,8
C -----------------------------------------------------------------
C REJECTION FOLLOWS
C -----------------------------------------------------------------
    7 IF(G3-(9*FF1(GG)-3*FF2(GG)-8*ALZ)/(9*FF1(0.)-3*FF2(0.)-8*ALZ))
     #  5,100,100
    8 IF(G3-(9*F12(GG)-3*F12(GG)-8*ALZ)/(9*F12(0.)-3*F12(0.)-8*ALZ))
     #  5,100,100
C ------------------------------------------------------------------
    5 U=RKSI
      RETURN
      END

C ___________________________________________________________________
      SUBROUTINE RPAIMS(E, U,NO,ALF)
c      SUBROUTINE RPAIMS(Z, E, U,NO,ALF)
C ___________________________________________________________________
C *******************************************************************    
C **************PAIR PRODUCTION.FULL SCREENING VERSION.**************
C *******************************************************************
C Z IS THE ATOMIC NUMBER OF THE MEDIUM
C U IS THE ENERGY OF ONE OF THE PRODUCET ELECTRONS
C NO IS THE NUMBER OF SAMPLES OF RKSI BEFORE REJECTION
C ALF IS THE SUM OF COEFFICIENS ALFAI
C *******************************************************************    
      COMMON/BLMASS/e_MASS,e2_MAS
c      COMMON / BLRADX / RADLEN
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI
      NO=0
      ALFA1=.6666-1./(36.*ALOG(183./Z**.3333))
      ALFA2=(1.3333+1/(9.*ALOG(183./Z**0.3333)))/12.
100   G1=s_rndm(0)
C     
C CHOISE OF THE TERM TO MODEL

      NO=NO+1

      IF(G1-ALFA1/(ALFA1+ALFA2) .LE. 0.) THEN    !  1,2,2
C     
C FIRST TERM
222     G2=s_rndm(0)*0.5            !  1/2 so that U<=1/2
        RKSI=G2
        if (rksi.eq.0..or.rksi.eq.1.) goto 222 

C     REJECTION
        DO J=1,100
        S=1.37E3 * SQRT( e_MASS/E/RKSI/(1.-RKSI)*RADLEN/XI(SOLD,Z) )
        IF (ABS((S-SOLD)/S) .LT. 1.E-5) GO TO 1
        SOLD=S
        END DO
        WRITE (*,*) 'LPM DID NOT CONVERGE MAC'
1       REJECT=(1.5*PSIFUN(S)-0.5*PHIFUN(S))*XI(S,Z)
        IF (s_rndm(0) .GT. REJECT) GO TO 100
        
      ELSE
C SECOND TERM
333     G2=s_rndm(0)*0.5            !  1/2 so that U<=1/2
        RKSI=-((0.125-G2/4.)**.33333333)+0.5
        if (rksi.eq.0..or.rksi.eq.1.) goto 333 
 
C     REJECTION
        DO J=1,100
        S=1.37E3 * SQRT( e_MASS/E/RKSI/(1.-RKSI)*RADLEN/XI(SOLD,Z) )
        IF (ABS((S-SOLD)/S) .LT. 1.E-5) GO TO 2
        SOLD=S
        END DO
        WRITE (*,*) 'LPM DID NOT CONVERGE MAC'
2       XIS=XI(S,Z)
        IF (s_rndm(0) .GT. PHIFUN(S)*XIS) GO TO 100

      END IF
      U=RKSI
      ALF=ALFA1+ALFA2
      RETURN
      END

C ____________________________________________________________
C ____________________________________________________________
      REAL FUNCTION RBRAI(EINIT,EMIN)
C Dont't need emin
C ____________________________________________________________
C *************************************************************
C **TOTAL CROSS-SECTION FOR BREMSSTRAHLUNG IN THE Moon ON 1 RL**
C ************************************************************* 
C  EMIN is minimum photon energy in MEV  1.E-3 <  EMIN  <  1.0
C  Uses 2 data matrices & logarithmic interpolation
C  works best for 1.E-3 & 1.E-2 
C  Low E and high EMIN ==> low x sections, not too accurate
C  Includes low E correction from Koch and Motz 0.1 < E < 10
C  No Elwert correction, DO NOT USE for E below ~ 100 keV
C  For EMIN not multiple of 10 it interpolates
C            Written E Zas 31 jan 1991
c Modified by Jaime Alvarez-Munhiz 17 October 2000
C ************************************************************

c To be filled in the external code:      
      COMMON/ BREMSS_XSECTION / BR_E(400),BR_X(400)
      EMIN=EMIN*1. 
      RBRAI=divdif(BR_X,BR_E,400,EINIT,2)
      return
      end
 

C _____________________________________________________________
      SUBROUTINE RBRMFC(E,U,N1,NO,ALF)
c      SUBROUTINE RBRMFC(Z,E,U,N1,NO,ALF)
C _____________________________________________________________
C *************************************************************
C BREMSSTRAHLUNG.SCREENING PARAMETER LE 2.
C Z IS THE ATOMIC NUMBER OF THE MEDIUM
C E  IS THE ELECTRON ENERGY
C U IS RELATIVE ENERGY OF THE RADIATED PHOTON
C N1 - CUT OFF.TRESHOLD ENERGY ETR = E*2**-N1
C NO IS THE NUMBER OF SAMPLES OF RKSI BEFORE REJECTION
C FCZ IS THE COULOMB CORRECTION TERM
C ALF IS THE SUM OF COEFFICIENS ALFAI
C *************************************************************
c      COMMON/VLZ/ALZ
c      COMMON/FCO/FCZ
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
     
      F12(B3)=21.12-4.184*ALOG(B3+.952)
      FF1(B1)=20.838-2.647*B1-1.496*B1**2+2.417*B1**3-.860*B1**4
      FF2(B2)=20.170-1.176*B2-2.382*B2**2+2.097*B2**3-.486*B2**4
      GAMA(A_A,B,C_C)=51.1*B/(A_A*(1.-B)*C_C**.3333)
C GAMMA IS THE SCREENING PARAMETER
      GGG(A_B)=1./(2.-A_B)
      NO=0
      ALFAI=.6923*(1.3333+.1111/(5.203-.3333*ALZ-FCZ))
  100 CONTINUE
      G1=s_rndm(0)
      I=1+INT(G1*N1)
C     
C CHOISE OF THE TERM TO MODEL
C     
      IF(G1-.5/(.5+N1*ALFAI)) 18,28,28
C FIRST TERM
   18 G1=s_rndm(0)
      NO=NO-1
      G2=s_rndm(0)
      G3=s_rndm(0)
      IF(G1-G2) 19,20,20
   19 RKSI=1.-G2+G1
      if (rksi.eq.0..or.rksi.eq.1.) then
        NO=NO-1
        goto 18
      end if
      GO TO 51
   20 RKSI=1.-G1+G2
      if (rksi.eq.0..or.rksi.eq.1.) then
        NO=NO-1
        goto 18
      end if
C     
C REJECTION FOLLOWS
   51 GG=1.36*GAMA(E,RKSI,Z)
      IF(GG-1) 21,22,22
   21 IF(G3-(3.*FF1(GG)-4.*ALZ-12.*FCZ)/(3.*FF1(0.)-4.*ALZ-12.*
     1FCZ)) 17,100,100
   22 IF(G3-(3.*F12(GG)-4.*ALZ-12.*FCZ)/(3.*F12(0.)-4.*ALZ-12.*
     1FCZ)) 17,100,100
C     
C SECOND TERM
   28 G2=s_rndm(0)
      NO=NO+1
      IF(G2-.72) 6,8,8
    6 G3=s_rndm(0)
      G5=s_rndm(0)
      RKSI=(1./2.**I)*G5
      if (rksi.eq.0..or.rksi.eq.1.) then
        NO=NO-1
        goto 28
      end if

      GG=1.36*GAMA(E,RKSI,Z)
      IF(GG-1.) 25,26,26
   25 IF(G3-(9.*FF1(GG)-3.*FF2(GG)-8.*ALZ-24.*FCZ)/(9.*FF1(0.)-3.*
     1FF2(0.)-8.*ALZ-24.*FCZ)) 17,100,100
   26 IF(G3-(9.*F12(GG)-3.*F12(GG)-8.*ALZ-24.*FCZ)/(9.*F12(0.)-3.*
     1F12(0.)-8.*ALZ-24.*FCZ)) 17,100,100
    8 G2=s_rndm(0)
      RKS=G2**.5
      BUMM=GGG(RKS)
      G3=s_rndm(0)
      IF(G3-BUMM) 10,10,8
   10 RKSI=1./2. **(I-1)-(1./2.**I)*RKS
      if (rksi.eq.0..or.rksi.eq.1.) then
        NO=NO-1
        goto 8
      end if

C     
C REJECTION FOLLOWS
c   62 GG=1.36*GAMA(E,RKSI,Z)
      GG=1.36*GAMA(E,RKSI,Z)
      G3=s_rndm(0)
      IF(GG-1) 15,16,16
   15 IF(G3-(9.*FF1(GG)-3.*FF2(GG)-8.*ALZ-24.*FCZ)/(9.*FF1(0.)-3.*
     1FF2(0.)-8.*ALZ-24.*FCZ)) 17,100,100
   16 IF(G3-(9.*F12(GG)-3.*F12(GG)-8.*ALZ-24.*FCZ)/(9.*F12(0.)-3.*
     1F12(0.)-8.*ALZ-24.*FCZ)) 17,100,100
   17 U=RKSI
      ALF=.5+N1*ALFAI
C ALF IS FUNCTION OF THE TRESHOLD ENERGY
      RETURN
      END

C _________________________________________________________________
      SUBROUTINE RBRSME(E,U,N1,NO,ALF)
c       SUBROUTINE RBRSME(Z,E,U,N1,NO,ALF)
C _________________________________________________________________
C *****************************************************************    
C BREMSSTRAHLUNG.SCREENING PARAMETER LE 2.
C Z IS THE ATOMIC NUMBER OF THE MEDIUM
C E  IS THE ELECTRON ENERGY
C U IS RELATIVE ENERGY OF THE RADIATED PHOTON
C N1 - CUT OFF.TRESHOLD ENERGY ETR = E*2**-N1
C NO IS THE NUMBER OF SAMPLES OF RKSI BEFORE REJECTION
C ALF IS THE SUM OF COEFFICIENS ALFAI
     
c      COMMON/VLZ/ALZ
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
     
      F12(B3)=21.12-4.184*ALOG(B3+.952)
      FF1(B1)=20.838-2.647*B1-1.496*B1**2+2.417*B1**3-.860*B1**4
      FF2(B2)=20.170-1.176*B2-2.382*B2**2+2.097*B2**3-.486*B2**4
      GAMA(A_A,B,C_C)=51.1*B/(A_A*(1.-B)*C_C**.3333)
C GAMMA IS THE SCREENING PARAMETER
      GGG(A_B)=1./(2.-A_B)
      NO=0
      ALFAI=.6923*(1.3333+.11111/(5.203-.3333*ALZ))
  100 CONTINUE
      G1=s_rndm(0)
      I=1+INT(G1*N1)
C     
C CHOISE OF THE TERM TO MODEL
C     
      IF(G1-.5/(.5+N1*ALFAI)) 18,28,28
C FIRST TERM
   18 G1=s_rndm(0)
      NO=NO+1
      G2=s_rndm(0)
      G3=s_rndm(0)
      IF(G1-G2) 19,20,20
   19 RKSI=1.-G2+G1
      if (rksi.eq.0..or.rksi.eq.1.) then
        NO=NO-1
        goto 18
      end if
      GO TO 51
   20 RKSI=1.-G1+G2
      if (rksi.eq.0..or.rksi.eq.1.) then
        NO=NO-1
        goto 18
      end if
C     
   51 GG=1.36*GAMA(E,RKSI,Z)
      IF(GG-1) 21,22,22
   21 IF(G3-(3.*FF1(GG)-4**ALZ)/(3.*FF1(0.)-4.*ALZ)) 17,100,100
   22 IF(G3-(3.*F12(GG)-4.*ALZ)/(3.*F12(0.)-4.*ALZ)) 17,100,100
C     
C SECOND TERM
   28 G2=s_rndm(0)
      NO=NO+1
      IF(G2-.72) 6,8,8
    6 G3=s_rndm(0)
      G5=s_rndm(0)
      RKSI=(1./2.**I)*G5
      if (rksi.eq.0..or.rksi.eq.1.) then
        NO=NO-1
        goto 28
      end if
      GG=1.36*GAMA(E,RKSI,Z)
      IF(GG-1) 25,26,26
   25 IF(G3-(9*FF1(GG)-3*FF2(GG)-8*ALZ)/(9*FF1(0.)-3*FF2(0.)-8*ALOG
     #(Z))) 17,100,100
   26 IF(G3-(9*F12(GG)-3*F12(GG)-8*ALZ)/(9*F12(0.)-3*F12(0.)-8*ALOG
     #(Z))) 17,100,100
    8 G2=s_rndm(0)
      RKS=G2**.5
      BUMM=GGG(RKS)
      G3=s_rndm(0)
      IF(G3-BUMM) 10,10,8
   10 RKSI=1./2.**(I-1)-(1./2.**I)*RKS
      if (rksi.eq.0..or.rksi.eq.1.) then
        NO=NO-1
        goto 8
      end if
C     
C REJECTION FOLLOWS
c   62 GG=1.36*GAMA(E,RKSI,Z)
      GG=1.36*GAMA(E,RKSI,Z)
      G3=s_rndm(0)
      IF(GG-1) 15,16,16
   15 IF(G3-(9*FF1(GG)-3*FF2(GG)-8*ALZ)/(9*FF1(0.)-3*FF2(0.)-8*ALOG
     #(Z))) 17,100,100
   16 IF(G3-(9*F12(GG)-3*F12(GG)-8*ALZ)/(9*F12(0.)-3*F12(0.)-8*ALOG
     #(Z))) 17,100,100
   17 U=RKSI
      ALF=.5+N1*ALFAI
C ALF IS FUNCTION OF THE TRESHOLD ENERGY
      RETURN
      END

C __________________________________________________________________
      SUBROUTINE RBRSMS(E,U,N1,NO,ALF)
c       SUBROUTINE RBRSMS(Z,E,U,N1,NO,ALF)
C __________________________________________________________________
C ******************************************************************
C BREMSSTRAHLUNG. REWRITTEN E. ZAS   13-DEC-1990
C LPM EFFECT FOR FULL SCREENING VERSION. 
C Z IS THE ATOMIC NUMBER OF THE MEDIUM
C U IS RELATIVE ENERGY OF THE RADIATED PHOTON
C NO IS THE NUMBER OF SAMPLES OF RKSI BEFORE REJECTION
C N1 - CUT OFF.TRESHOLD ENERGY ETR = E*2**-N1
C UCRIT, NCRIT - BREAKPOINTS REGIONS 2 & 3
C ALF IS THE SUM OF COEFFICIENS ALFAI
C WHEN U < UCRIT THE LPM EFFECT IS APROXIMATED AS ~ U^{-1/2}
C THIS CORRESPONDS TO SECOND TERM
C ******************************************************************
c     ALN2=log(10.)
      PARAMETER (PI = 3.141592654, ALN2 = 0.6931472, B = 8.859E-10)
C       B is the relevant factor that appears in LPM calculation (E0 in MeV) 

      COMMON/BLMASS/e_MASS,e2_MAS
c      COMMON / BLRADX / RADLEN 
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI
      GGG(A_A)=1./(2.-A_A)
      N1=N1*1
      NO=0
      UCRIT=B*E/(1.+B*E)
      NCRIT=-LOG(UCRIT)/ALN2
      ZETA=UCRIT*2.**float(NCRIT)
      ZETREF=(1.-ZETA)/ZETA
      ALFACT=1.3333333+0.11111111/(5.2523-.33333*ALOG(Z))
      SUMA3=0.5
      SUMA2=SUMA3+5.19291E4*ALFACT/SQRT(E)*
     #      ( ASIN(SQRT(UCRIT))+SQRT(UCRIT*(1.-UCRIT)) )
      SUMA1=SUMA2+(-LOG(ZETA)-(1.-ZETA))*ALFACT
      SUMALF=SUMA1+ALN2*ALFACT*NCRIT
100   G1=s_rndm(0)
      G2=s_rndm(0)
C      I=1+INT(G2*INTN1)
      NO=NO+1

C           CHOICE OF THE TERM TO MODEL

      IF (G1 .LT. SUMA3/SUMALF) THEN
C     THIRD TERM
222     G1=s_rndm(0)
        G2=s_rndm(0)
        RKSI=1.-ABS(G2-G1)
        if (rksi.eq.0..or.rksi.eq.1.) goto 222
C     REJECTION
        DO J=1,100
        S=1.37E3 * SQRT( e_MASS/E*RKSI/(1.-RKSI)*RADLEN/XI(SOLD,Z) )
        IF (ABS((S-SOLD)/S) .LT. 1.E-5) GO TO 3
        SOLD=S
        END DO
        WRITE (*,*) 'LPM DID NOT CONVERGE MAC'
3       XIS=XI(S,Z)
        IF (s_rndm(0) .GT. PSIFUN(S)*XIS) GO TO 100

      ELSE IF (G1 .LT. SUMA2/SUMALF) THEN
C     SECOND TERM (low u region of LPM)

22      G1=s_rndm(0)
        RKSI=UCRIT*G1**2
        if (rksi.eq.0..or.rksi.eq.1.) goto 22
c Error found by Enrique Marques 1./SQRT(1.-RKSI) -> SQRT(1.-RKSI)
        IF (s_rndm(0) .GT. SQRT(1.-RKSI)) GO TO 22  ! Rejection

C     REJECTION
        DO J=1,100
        S=1.37E3 * SQRT( e_MASS/E*RKSI/(1.-RKSI)*RADLEN/XI(SOLD,Z) )
        IF (ABS((S-SOLD)/S) .LT. 1.E-5) GO TO 2
        SOLD=S
        END DO
        WRITE (*,*) 'LPM DID NOT CONVERGE MAC'
2       IF (s_rndm(0) .GT. PHIFUN(S)/6./S*SQRT(XI(S,Z)/2.)) GO TO 100

C     FIRST TERM                   !  (FROM UCRIT TO 2^-NCRIT)     
      ELSE IF (G1 .LT. SUMA1/SUMALF) THEN
C                  Choose variate form (1-X)
11      G1=s_rndm(0)
        G2=s_rndm(0)
        RKSI=ABS(G1-G2)
        if (rksi.eq.0..or.rksi.eq.1.) goto 11
        IF (s_rndm(0) .GT. 1./(1.+ZETREF*RKSI)) GO TO 11  ! Rejection

C     REJECTION
        DO J=1,100
        S=1.37E3 * SQRT( e_MASS/E*RKSI/(1.-RKSI)*RADLEN/XI(SOLD,Z) )
        IF (ABS((S-SOLD)/S) .LT. 1.E-5) GO TO 1
        SOLD=S
        END DO
        WRITE (*,*) 'LPM DID NOT CONVERGE MAC'
1       IF (s_rndm(0) .GT. PHIFUN(S)*XI(S,Z)) GO TO 100

      ELSE              !  ZEROTH TERM, 
C       Total of NCRIT terms
        I=1.+s_rndm(0)*NCRIT    ! Select which of the NCRIT terms
        G2=s_rndm(0)

        IF (G2 .LT. 0.72) THEN 
444       G3=s_rndm(0)
          RKSI=(1./2.**I)*G3
          if (rksi.eq.0..or.rksi.eq.1.) goto 444 
        ELSE
16        G2=s_rndm(0)
          RKS=G2**.5
          BUMM=GGG(RKS)
          G3=s_rndm(0)
          IF(G3-BUMM) 17,17,16
17        RKSI=1./2.**(I-1.)-(1./2.**I)*RKS
          if (rksi.eq.0..or.rksi.eq.1.) goto 16 
        END IF
C     REJECTION
        IF (RKSI .LT. UCRIT) GO TO 100        
        DO J=1,100
        S=1.37E3 * SQRT( e_MASS/E*RKSI/(1.-RKSI)*RADLEN/XI(SOLD,Z) )
        IF (ABS((S-SOLD)/S) .LT. 1.E-5) GO TO 10
        SOLD=S
        END DO
        WRITE (*,*) 'LPM DID NOT CONVERGE MAC'
10      IF (s_rndm(0) .GT. PHIFUN(S)*XI(S,Z)) GO TO 100
      END IF

      U=RKSI
      ALF=SUMALF
      RETURN
      END


C ______________________________________________________________
      REAL FUNCTION RCOMP(E)
c       REAL FUNCTION RCOMP(Z,E)
C ______________________________________________________________
C       cEl: pxlHOTO KOMpTxHOBO CE~EHiE
C      BXOd: Z-ATOMEH HOMEP
C        E - Energy HA fOTOHA(MEV)

      COMMON/BLMASS/e_MASS,e2_MAS
c      COMMON /BLNORM/ FACNOR   !  =109.9667/ALZ3/(Z+VSI) 
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS

      EK1=E/e_MASS
      C1=1./EK1**2
      C2=1.-2.*(1.+EK1)*C1
      C3=(1.+2.*EK1)*C1
      E1=E/(1.+2.*E/e_MASS)
      E2=E
      EE1=E1/E
      EE2=E2/E
      R=C1*(1./EE1-1./EE2)+C2*ALOG(EE2/EE1)
     #+EE2*(C3+EE2/2.)-EE1*(C3+EE1/2.)
C     
      RCOMP=FACNOR/2./E*R          !  FACNOR=109.9667/ALZ3/(Z+VSI)
C     
      RETURN
      END

C ___________________________________________________________________
      SUBROUTINE COMPTU(E,U,V,N1,NO,ALF)
C ___________________________________________________________________
C *******************************************************************
C COMPTON EFFECT
C E  IS THE PHOTON ENERGY
C U IS THE ENERGY OF SECONDARY PHOTON
C V IS THE SECONDARY PHOTON PRODUCTION ANGLE COS
C NO IS THE NUMBER OF SAMPLES OF RKSI BEFORE REJECTION
C ALF IS THE SUM OF COEFFICIENS ALFAI
C REWRITTEN by E. Zas 8 Feb 1991 following EGS-4
C ********************************************************************
      COMMON/BLMASS/e_MASS,e2_MAS
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI
      GAMA(A,B)=1.-(A*B)/(1.+A**2)
      N1=N1*1
      NO=0
      ERATIO=E/e_MASS
      ENOUGHT=1./(1.+2.*ERATIO)
      ALFA1=-LOG(ENOUGHT)
      ALFA2=(1.-ENOUGHT**2)/2.
      ALF=ALFA1+ALFA2

    1 G1=s_rndm(0)
      NO=NO+1

C        Choice of term in probability decomposition

      IF (G1 .LE. ALFA1/ALF) THEN    !! first term
        RKSI=ENOUGHT*EXP(ALFA1*s_rndm(0))
      ELSE                           !! second term  
        G1=s_rndm(0)
        G2=s_rndm(0)
          IF (G1 .GT. E/(E+e_MASS) ) THEN  !! two more terms for second term
            EPRIME=G2
          ELSE
            EPRIME=MAX(G2,s_rndm(0))
          END IF
        RKSI=(1.-ENOUGHT)*EPRIME+ENOUGHT
      END IF

      if (rksi.eq.0.) then 
         NO=NO-1
         goto 1
      end if      

C        Rejection test follows

      T1_COS=(1.-RKSI)/RKSI/ERATIO
      SIN2=T1_COS*(2.-T1_COS)            
      IF ( s_rndm(0) .GT. GAMA(RKSI,SIN2) ) GO TO 1

      U=RKSI
      V=1.-T1_COS

      RETURN
      END

C _______________________________________________________________
      REAL FUNCTION RBHA(E)
c       REAL FUNCTION RBHA(E,AE,Z)
C _______________________________________________________________
C     
C    BHABHA - SCATTERING
C     
C     
C  THE TOTAL CROSS SECTION OF BHABHA-SCATTERING
C   -E  IS THE INCIDENT POSITRON ENERGY (MEV)
C   -AE  IS THE CUTOFF ENERGY (MEV)
C   -Z  IS ELEMENT'S ATOM NUMBER ( OF ABSORBER )
C     
      COMMON/BLMASS/e_MASS,e2_MAS
c      COMMON /BLNORM/ FACNOR   !  =109.9667/ALZ3/(Z+VSI) 
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
      IF (E.LE.CH_THR) THEN 
        RBHA=0.0
        RETURN
      END IF
      G=E/e_MASS
      BETA2=1.-1./G/G
      E1=(CH_THR-e_MASS)/(E-e_MASS)
      E2=1.0
      Y=1.0/(1.0+G)
      B4=(1.0-2.0*Y)**3
      B3=B4+(1.0-2.0*Y)*(1.0-2.0*Y)
      B2=(1.0-2.0*Y)*(3.0+Y*Y)
      B1=2.0-Y*Y
      RBHA=FACNOR*((1.0/E1-1.0/E2)/BETA2-B1*ALOG(E2/E1)+B2*(E2-E1)+
     #(E2*B4/3.0-B3/2.0)*E2*E2-(E1*B4/3.0-B3/2.0)*E1*E1)/(E-e_MASS)
      RETURN
      END

C __________________________________________________________________
      SUBROUTINE RBHAS(E,EP,EE,V3,V4)
c      SUBROUTINE RBHAS(E,AE,EP,EE,V3,V4)
C __________________________________________________________________
C ******************************************************************    
C  DIFFERENTIAL CROSS SECTION OF BHABHA-SCATTERING
C   -E  IS THE INCIDENT POSITRON ENERGY (MEV)
C   -EE  IS THE SCATTERED ELECTRON ENERGY (MEV)
C   -EP  IS THE SCATTERED POSITRON ENERGY (MEV)
C   -AE  IS THE CUTOFF ENERGY (MEV)
C   -V3  AND  V4  ARE SECONDARY PRODUCTION ANGLES (COS)
C   -V3 - ANGLE FOR SCATTERED POSITRON (COS)
C   -V4 - ANGLE FOR SCATTERED ELECTRON (COS)
C ******************************************************************    
C     
C  TPqbBA  E  dA E pO-gOlqMO OT  AE !!!
C     
      COMMON/BLMASS/e_MASS,e2_MAS
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI

      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
     
      IF (CH_THR.GE.E) THEN 
        WRITE(3,2)
        RETURN
      END IF 
      N=0
      G=E/e_MASS
      GINV=1./G
      BETA2=1.0-1.0/G/G
      E0=(CH_THR-e_MASS)/(E-e_MASS)
      Y=1.0/(1.0+G)
      B4=(1.0-2.0*Y)**3
      B3=B4+(1.0-2.0*Y)**2
      B2=(1.0-2.0*Y)*(3.0+Y*Y)
      B1=2.0-Y*Y
    1 PSI=s_rndm(0)
      ER=E0/(1.0-(1.0-E0)*PSI)
      GE=(1.0-E0)*(1.0/BETA2-ER*(B1-ER*(B2-ER*(B3-ER*B4))))
      HI=s_rndm(0)
      IF(HI.GT.GE) GO TO 1
      EE=ER*(E-e_MASS)+e_MASS
      EP=(E+e_MASS)-EE

C      P0=E*E-e2_MAS
C      P1=EP*EP-e2_MAS
C      P2=EE*EE-e2_MAS
C      V3=(P0+P1-P2)/2.0/SQRT(P0*P1)
C      V4=(P0+P2-P1)/2.0/SQRT(P0*P2)

      G1G1=(1.+GINV)/(1.-GINV)
      U1G=ER*(1.-GINV)
      V3=SQRT( G1G1*(1.-U1G-GINV)/(1.-U1G+GINV) )   ! e+ 
      V4=SQRT( G1G1*U1G/(U1G+2.*GINV) )             ! SCATTERED e-
2     FORMAT(//5X,' Incomprhensible junk RBHAS',//5X,'TPqbBA  E>CH_THR')
      RETURN
      END

C ________________________________________________________________
      REAL FUNCTION RMOL(E)
c       REAL FUNCTION RMOL(Z,AE,E        )
C ________________________________________________________________
C ****************************************************************    
C          cEl:  pxlHOTO MxOlEPOBO CE~EHiE
C          BXOd:   Z - ATOMHiq HOMEP HA CPEdATA
C                  E - pxlHATA Energy HA ElEKTPOHA (MEV)
C                  AE - pPAgOBATA Energy (MEV); pPAgxT zA
C                  diCKPETHOTO MxOlEPOBO PAzCEjBAHE E
C                  ETH = 2.*AE - e_MASS
C     
      COMMON/BLMASS/e_MASS,e2_MAS
c      COMMON /BLNORM/ FACNOR   !  =109.9667/ALZ3/(Z+VSI) 
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
   
      ETH=2.*CH_THR - e_MASS
         IF (E.LE.ETH) THEN 
           RMOL=0.
           RETURN
         END IF

      GAMA=E/e_MASS
      GINV=1./GAMA
      G2=GAMA*GAMA
      TK0=E-e_MASS
      BETA2=1.-1./G2
      C1=(1.-GINV)**2
      C2=(2.-GINV)/GAMA
      E1=(CH_THR-e_MASS)/TK0
      E2=0.5
      E11=1.-E1
      E21=1.-E2
      R=FACNOR/TK0/BETA2
C     
      RMOL=R*(C1*(E2-E1)+1./E1-1./E2+1./E21-1./E11-
     *        C2*ALOG(E2*E11/E1/E21))
C     
      RETURN
      END

C ______________________________________________________________
      SUBROUTINE RMOLS(E,ES,V3,V4)
c            SUBROUTINE RMOLS(E,ES,AE,V3,V4      )
C ______________________________________________________________
C          MOLLER SCATTERING
C          E  IS  THE INCIDENT ELECTRON ENERGY (MEV)
C          ES  IS  THE SCATTERED ELECTRON ENERGY (MEV)
C          AE  IS  THE CUTOFF ENERGY
C          V3  ANGLE COSINE FOR LEAST ENERGETIC SCATTERED ELECTRON (ES)
C          V4  ANGLE COSINE FOR MOST ENERGETIC ELECTRON (E)
C          REWRITTEN BY E. ZAS  31 - JAN - 1991
C     
      COMMON/BLMASS/e_MASS,e2_MAS
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI

      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS

      TK0=E-e_MASS
      GAMA=E/e_MASS
      GINV=1./GAMA
      GA=GAMA*GAMA
      BETA2=1.-1./GA
      TE=CH_THR-e_MASS
      ETH=2.* CH_THR     -e_MASS
C     
      E0=TE/TK0
      G1=(1.-2.*E0)/BETA2
      G2=(1.-GINV)**2
      G3=(2.-GINV)/GAMA
C     
    1 S=s_rndm(0)
      U=TE/(TK0-(E-ETH)*S)
      R=U/(1.-U)
C     
C        REJECTION FOLLOWS
C     
      GG=G1*(1.+G2*U*U+R*(R-G3))
      S=s_rndm(0)
        IF(S .GT. GG)  GO TO 1
      ES=U*TK0+e_MASS
     
C      P1=SQRT(E*E-e2_MAS)
C      P3=SQRT(ES*ES-e2_MAS)
C      V3=((E+e_MASS)*ES-E*e_MASS-e2_MAS)/P1/P3
C      E4=E-ES+e_MASS
C      P4=SQRT(E4*E4-e2_MAS)
C      V4=((E+e_MASS)*E4-E*e_MASS-e2_MAS)/P1/P4

      G1G1=(1.+GINV)/(1.-GINV)
      U1G=U*(1.-GINV)
      V3=SQRT( G1G1*U1G/(U1G+2.*GINV) )           ! SMALL EN., SCATTERED e 
      V4=SQRT( G1G1*(1.-U1G-GINV)/(1.-U1G+GINV) ) ! LARGE ENERGY
     
      RETURN
      END

C ___________________________________________________________________
c      REAL FUNCTION RANH(E,Z)
      REAL FUNCTION RANH(E)
C ___________________________________________________________________
C *******************************************************************
C  THE TOTAL CROSS SECTION OF TWO-PHOTON POSITRON-ELECTRON
C  ANNIHILATION
C    E - ENERGY OF INCIDENT POSITRON (MEV)
C    Z - ELEMENT'S  ATOM  NUMBER (OF  ABSORBER )
C *******************************************************************
      COMMON/BLMASS/e_MASS,e2_MAS
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
c      COMMON /BLNORM/ FACNOR   !  =109.9667/ALZ3/(Z+VSI) 
      G=E/e_MASS
      GS=SQRT(G*G-1.)
      ANH=(ALOG(G+GS)*(G*G+4.*G+1.)/(G*G-1.)-(G+3.)/GS)
      RANH=FACNOR*ANH/(G+1.)/2./e_MASS
      RETURN
      END

C __________________________________________________________________
      SUBROUTINE RANHS(E,EBIG,ESMAL,CBIG,CSMAL)
C __________________________________________________________________
     
C Program zA iz~iClqBAHE HA CE~EHiETO zA dBufOTOHHA AHiXilAciq
C     
C pAPAMETxPxT  E  CE zAdABA B  MEV  i TPqbBA dA E pO-gOlqM OT  1 !
C pAPAMETPiTE  EBIG  i  ESMAL  CA B  MEV.
C pAPAMETPiTE  CBIG  i  CSMAL  CA CxOTBETHO KOCiHuCiTE HA BiCOKOEHEPgiTi
C HiCKOEHEPgiTi~Hiq fOTOH ( OT pOlqPHiTE xgli ) , KATA AziMuTAlHiTE xgli
C PABHOMEPHO PAzpPEdElEHi B iHTEPBAlA  (0,2p)
     
      COMMON/BLMASS/e_MASS,e2_MAS
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI

      N=0
      GAMA=E/e_MASS
      A=GAMA+1.0
      T=GAMA-1.0
      P=SQRT(A*T)
      COSFAC=SQRT(A/T)   !   A/P
      AL=ALOG(GAMA+P)
    1 PSS=s_rndm(0)
      IF(PSS.LT.0.) WRITE(*,*) 'PSS < 0 !!!!'
      EE=EXP(PSS*AL)/(A+P)
      GE=1.0-EE+(2.0*GAMA-1.0/EE)/(A**2)
      HI=s_rndm(0)
      IF (HI.GT.GE) GO TO 1
      EA=1.0-EE
      EBIG=(E+e_MASS)*MAX(EA,EE)
      ESMAL=(E+e_MASS)-EBIG
      CBIG=COSFAC*(1.-e_MASS/EBIG)
      CSMAL=MAX(COSFAC*(1.-e_MASS/ESMAL),-1.)

      RETURN
      END

C ________________________________________________________________
      REAL FUNCTION RPHOT(E)
c       REAL FUNCTION RPHOT(E,Z)
C ________________________________________________________________
C     
C     
C        cEl:  pxlHOTO CE~EHiE zA fOTOEfEKTA
C        BXOd:  E - EnergyTA HA fOTOHA (MEV)
C               Z - ATOMHiq HOMEP HA CPEdATA
C        METOd: ApP. fOPMulA HA HUBBELL
C     
c      COMMON/VT0/T0
c      COMMON/VAKO/AKO
C     
C               T0 - dxlviHATA HA PAd.Ed. B G/SM**2
C               AKO - KOEficiEHT zA pPEXOdA B G/CM**2
C        AKO=NA/A*10**(-24)
c      COMMON/VLZ/ALZ
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS      
      
      DIMENSION AN(4),BN(4),CN(4),PN(4)
      DATA AN/1.6268E-09,1.5274E-09,1.1330E-09,-9.12E-11/,
     *     BN/-2.683E-12,-5.110E-13,-2.177E-12,0./,
     *     CN/4.173E-02,1.027E-02,2.013E-02,0./,
     *     PN/1.,2.,3.5,4./
C     
      R=1.+(0.01481-0.000788*ALZ)*ALZ**2
      SUM=0.
      DO 1 I=1,4
      SUM=SUM+(AN(I)+BN(I)*Z)/(1.+CN(I)*Z)/E**PN(I)
    1 CONTINUE
C     
      RPHOT=AKO*T0*SUM*R*Z**5
      RETURN
      END

C _______________________________________________________________
      SUBROUTINE RIONFN(I,E,T,TE1,DEI,AP)
c      SUBROUTINE RIONFN(C,A,RM,X0,X1,Z,I,AI,
c     #E,T,TE1,DEI,AP)
c EGS4 pages 66-75
C _______________________________________________________________
C ***************************************************************
C     CALCULATES THE IONIZATION LOSS BY FORD & NELSON
C        E  -  PARTICLE ENERGY
C        DEI - jOH.zAgubi(MEV) HA pxT  T (B PAd.Ed.)
C        I = 0  -  ELECTRON
C        I = 1  -  POSITRON
C        AI = I/M
C        AP = Photon threshold
C     TE1  E  KiHETic  pPAgOBA   Energy ( B  EdiHici  MC2)
C ***************************************************************
      COMMON/BLMASS/e_MASS,e2_MAS
c      COMMON /BLNORM/ FACNOR   !  =109.9667/ALZ3/(Z+VSI) 
      COMMON/BET2/BETA2,BETA2LN
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
     
      G=E/e_MASS
c      write(*,*) G
      G2=G**2
      BETA2=1.-1./G2
       IF (G .LT. 10.) THEN
        BETA2LN=LOG(BETA2)
       ELSE
        BETA2LN=-1./G2*(1.+(0.5+0.33333333/G2)/G2)
       END IF
      ETA=SQRT(BETA2)*G
      TAU=G-1.
      Y=1./(G+1.)
      TMAX1=TAU/2.
           IF(I.EQ.1) TMAX1=TAU
      DE=AMIN1(TMAX1,TE1)
      D2=DE*DE
C     
           IF(I.EQ.1) GO TO 1
C     
C     AKO  I=0  ~ACTicATA  E  ElEKTPOH
      F=ALOG((TAU-DE)*DE)+TAU/(TAU-DE)+(D2/2.+(2.*TAU+1.)*ALOG(1.-DE/
     *TAU))/(G   *G   )   -1.-BETA2
           GO TO 2
C     AKO  I=1  ~ACTicATA  E  pOziTPOH
    1 F=ALOG(TAU*DE)-BETA2*(TAU+2.*DE-1.5*D2*Y-(DE-DE*D2/3.)*Y*Y-
     *(D2/2.-TAU*DE*D2/3.+D2*D2/4.)*Y*Y*Y)/TAU
    2 CONTINUE
C     
C     
      X=ALOG10(ETA)
      IF(X.GE.X1) D=4.606*X+C
      IF((X.LT.X1).AND.(X.GE.X0)) D=4.606*X+C+A*(X1-X)**RM
      IF(X.LT.X0) D=0.
C     
C     
      AL=ALOG(2.*(TAU+2.)/(AI*AI))
C     jOHizAciOHHiTE  zAgubi  HA  pxT  T (B  PAd.Ed.)
C     
        DEI=T*FACNOR*((AL+F-D)/BETA2    ! FACNOR=109.9667/ALZ3/(Z+VSI) 
     #      + AP*((1.-0.5*AP/E)*(4.+1./LOG(183./2.))+(AP/E)**2)/3.)
c      WRITE(*,*) I,E,T,FACNOR,AL,F,D,BETA2,AP,DEI
      RETURN
      END

C ________________________________________________________________
      SUBROUTINE MULTSC(E1,X,R3)
C ________________________________________________________________
C ****************************************************************
C MULTIPLE COULOMB SCATTERING PROGRAM. DESIGN FOR AIR
C E1 IS THE ELECTRON ENERGY
C X SI THE SCATTERING ANGLE
C R3 IS THE PATH
C *****************************************************************
c      COMMON/BBB/BEE
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI
      COMMON/BET2/BETA2,BETA2LN
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
     
      FF1(A_A)=2.*EXP(-A_A**2)
      FF2(B,D)=2.*EXP(-B**2)*(B**2-1.)*(D-ALOG(B**2))-
     #(2.*(1.-2.*EXP(-B**2)))
C -----------------------------------------------------------------
C BE IS "B" FOR AIR. SEE THE LONG WRITE-UP.
C -----------------------------------------------------------------
      BE=ALOG(R3)+BEE-BETA2LN       ! EGS4 (2.14.17-18)
       IF(BE.LT.0.9) GO TO 40
C -----------------------------------------------------------------
C TOO THIN LAYER. ONLY THE MEAN SCATT.ANGLE IS CALCULATED
C -----------------------------------------------------------------

      F1=BE+ALOG(BE+ALOG(BE+ALOG(BE+ALOG(BE))))
      IF(F1-5.) 1,1,2
    1 ALFA1=0.
      GO TO 33
    2 ALFA1=1.-5./F1
   33 ALFA2=4.1509/F1
      ALFA3=2.1946/F1
C ------------------------------------------------------------------
C CHOISE OF ONE OF THE THREE TERMS
C ------------------------------------------------------------------
   30 G2=s_rndm(0)
      IF(G2-ALFA1/(ALFA1+ALFA2+ALFA3)) 4,4,13
   13 IF(G2-(ALFA1+ALFA2)/(ALFA1+ALFA2+ALFA3)) 6,6,7
C ---------------------
C FIRST TERM
C ---------------------
    4 G1=s_rndm(0)
      if (g1.eq.0.) goto 4
      RKSI = ( -ALOG(G1))**.5
      if (rksi.eq.0.) goto 4
      GO TO 11
C ---------------------
C SECOND TERM
C ---------------------
    6 RKSI=s_rndm(0)
      G3=s_rndm(0)
      if (rksi.eq.0.) goto 6
      GO TO 9
C ---------------------
C THIRD TERM
C ---------------------
    7 G1=s_rndm(0)
      G3=s_rndm(0)
      if (g1.eq.0.) goto 7
      RKSI = 1./SQRT(G1)
      GO TO 10
C ------------------------------------------------------------------
C REJECTION FOR THE SECOND AND THIRD TERM
C ------------------------------------------------------------------
    9 CALL TUPAR(RKSI,BI)
      GAMA2=(RKSI/4.1509)*(5.*FF1(RKSI)+FF2(RKSI,BI))
      IF ( G3-GAMA2 ) 11,11,30
   10 IF(RKSI.GE.8.) GO TO 30
      CALL TUPAR(RKSI,BI)
      GAMA3=(RKSI**4/4.3892)*(5.*FF1(RKSI)+FF2(RKSI,BI))
      IF ( G3-GAMA3 ) 11,11,30
c Change For Different Medium
c EGS4 manual pages 77,79
c X=capital theta=theta*B**0.5*Kappa_c and
c what it is sampled is Kappa_cc*sqrt(distance)
c Kappa_cc Sqrt(T0 in cm) = FACTOR_MS 
c FACTOR_MS is calculated in zetas_salt.f
c See EGS4 manual (2.14.24-25) and Note in between.
c Recalculated by Jaime Alvarez-Muniz 18 October 2000
c   11 FACTOR_MS=5.1249361
   11 X=RKSI*FACTOR_MS*SQRT(R3*F1)/E1/BETA2
      G4=s_rndm(0)
      IF((G4*G4*X).GE.(SIN(X))) GO TO 30
   41 RETURN
   40 X=SQRT(R3)*21./E1/BETA2
      GO TO 41
      END

C _______________________________________________________________
      SUBROUTINE TUPAR(RKSI,BETA)
C _______________________________________________________________
      N=1
      BETA=0
      Y=RKSI**2
      IF(Y-4.) 10,11,11
   10 Y2=1.
      DO 5 I=1,10
c    9 N=N*I
      N=N*I
      Y2=Y*Y2
      BETAI=Y2/N/I
      BETA=BETA+BETAI
      IF(BETAI.LT.1.E-20) GO TO 6
    5 CONTINUE
    6 BETA=BETA+ALOG(Y)+.57721
      RETURN
   11 Y2=1
      DO 15 I=1,5
      N=N*I
      Y2=Y2*Y
      BETAI=N/Y2
   15 BETA=BETA+BETAI
      BETA=((BETA+1)/Y)*EXP(Y)
      RETURN
      END

C ____________________________________________________________________
      SUBROUTINE GEOM(E1,DEI,T1,T,TIX,T2,DL,DM,DN,DNHALF,DX,DY,DZ,dCTi)
C ____________________________________________________________________
      DIMENSION A_arr(9),B(3),R9(3)
      COMMON/BLMASS/e_MASS,e2_MAS
c      COMMON/VT0/T0
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
      QQ(P,U)=SQRT(P*P+U*U)

      T2=T2*1.
C      IF(DN.GT.1.) DN=1.
C      T2=ACOS(DN)

      F=QQ(DL,DM)
      A_arr(1)=DL*DN/F
      A_arr(2)=-DM/F
      A_arr(3)=DL
      A_arr(4)=DM*DN/F
      A_arr(5)=DL/F
      A_arr(6)=DM
      A_arr(7)=-F
      A_arr(8)=0.
      A_arr(9)=DN

    1 R1=s_rndm(0)
      R2=s_rndm(0)
      RR=R1*R1+R2*R2
      IF(RR.GE.1) GO TO 1
      C_val=(R1*R1-R2*R2)/RR
      S=2*R1*R2/RR
      R=s_rndm(0)
      IF(R.LT..5) S=-S
      SINTIX=SIN(TIX)
      B(1)=C_val*SINTIX
      B(2)=S*SINTIX
      B(3)=COS(TIX)
      CALL MTXMLT(A_arr,B,R9)
      DLHALF=0.5*(DL+R9(1))
      DMHALF=0.5*(DM+R9(2))
      DNHALF=0.5*(DN+R9(3))

      DT=T*T0
      DX=DT*DLHALF
      DY=DT*DMHALF
      DZ=DT*DNHALF
      T1=T1+DT*DNHALF

      DL=R9(1)
      DM=R9(2)
      DN=R9(3)
      
      IF ((E1 .LE. 5.) .OR. (E1 .LT. 100.*DEI)) THEN
        RATIO2=(e_MASS/E1)**2
        dCTi=DT*( E1/DEI*( SQRT(1.-RATIO2)
     #                -SQRT((1.-DEI/E1)**2-RATIO2) ) -DNHALF)
      ELSE IF (E1 .LT. 100.) THEN
        RATIO2=(e_MASS/E1)**2
        dCTi=DT*(0.5*RATIO2*(1.+0.75*RATIO2+DEI/E1)+1.-DNHALF)
      ELSE
        dCTi=DT*(1.-DNHALF)
      END IF

      RETURN
      END

C ____________________________________________________________________
      SUBROUTINE GEOGAM(T1,T,DL,DM,DN,DX,DY,dCT)
C ____________________________________________________________________
c      COMMON/VT0/T0
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
     
      DT=T*T0
      DX=DT*DL
      DY=DT*DM
      DZ=DT*DN
      T1=T1+DZ
      dCT=dCT+DT*(1.-DN)       ! Possible error dCT=dCT+DT*(1.-DN)
                               ! instead of     dCT=dCT+DZ*(1.-DN)

      RETURN
      END

C ________________________________________________________________
      SUBROUTINE GEONEW(DL,DM,DN,CT3,CT4,DL3,DM3,DN3,DL4,DM4,DN4)
C ________________________________________________________________
C ****************************************************************
C  Program, PE{ABA}A geometric zAdA~A B Clu~Eq, KOgATO
C  CE PAvdAT dBE HOBi ~ACTici
C    BXOdq}i pAPAMETPi:
C      -DL,DM,DN -  COS-diPEKTOPiTE HA pxPBic ~ACTicA
C      -CT3,CT4 - KOCiHuCiTE HA pOlqPHiTE xgli HA BTOPi~HiTE ~ACTici
C    izXOdq}i pAPAMETPi:
C                      xgxl  CT3
C      -DL3,DM3,DN3 -  COS-diPEKTOPiTE HA BTOPic ~ACTicA C pOlqPEH
C                      C pOlqPEH xgxl  CT4
C      -CT3,CT4,ST3,ST4 CA CxOTBETHO  COS  i  SIN  OT pOlqPHiTE xgli
C       HA dBETE BTOPi~Hi ~ACTici ( TiTA3  i  TiTA4 )
C      -CF3,CF4,SF3,SF4  CA CxOTBETHO  COS  i  SIN  OT AziMuTAlHiTE xgli
C       HA BTOPi~HiTE ~ACTici ( fi3  i  fi4 )
C *****************************************************************
      DIMENSION A(9),B(3),R9(3)
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI
c    5 CONTINUE
      F=SQRT(DL*DL+DM*DM)
c    6 CONTINUE
      A(1)=DL*DN/F
      A(2)=-DM/F
      A(3)=DL
      A(4)=DM*DN/F
      A(5)=DL/F
      A(6)=DM
      A(7)=-F
      A(8)=0.
      A(9)=DN
c    7 CONTINUE
    1 R1=s_rndm(0)
      R2=s_rndm(0)
c    8 CONTINUE
      RR=R1*R1+R2*R2
      IF(RR.GE.1.0) GO TO 1
c    9 CONTINUE
      CF3=(R1*R1-R2*R2)/RR
      SF3=2.0*R1*R2/RR
c   10 CONTINUE
      R=s_rndm(0)
      IF(R.LT.0.5)  SF3=-SF3
c   11 CONTINUE
      CT3R2=MIN(1.,CT3**2)
      ST3=SQRT(1.0-CT3R2)
      B(1)=CF3*ST3
      B(2)=SF3*ST3
      B(3)=CT3
c   12 CONTINUE
      CALL MTXMLT(A,B,R9)
c   13 CONTINUE
      DL3=R9(1)
      DM3=R9(2)
      DN3=R9(3)
      SF4=-SF3
      CF4=-CF3
      CT4R2=MIN(1.,CT4**2)
      ST4=SQRT(1.0-CT4R2)
      B(1)=ST4*CF4
      B(2)=ST4*SF4
      B(3)=CT4
c   14 CONTINUE
      CALL MTXMLT(A,B,R9)
      DL4=R9(1)
      DM4=R9(2)
      DN4=R9(3)
c   20 CONTINUE
      RETURN
      END

C ___________________________________________________________________
      SUBROUTINE MTXMLT(A,B,R)
C ___________________________________________________________________
      DIMENSION A(9),B(3),R(3)
      M=3
      N=1
      L=3
      DO 2 I=1,M
      DO 2 J=1,N
      KR=(I-1)*N+J
      R(KR)=0.
      DO 1 I1=1,L
      KB=(I1-1)*N+J
      KA=(I-1)*L+I1
    1 R(KR)=R(KR)+A(KA)*B(KB)
    2 CONTINUE
      RETURN
      END

C ________________________________________________________________
      SUBROUTINE DIST(NO,SIGM1,T)
C ________________________________________________________________
C     COMMON IX,CP,CQ,CS,TMETR,T2,ETR,HUI
      A=1.
      N1=0
      DO 1 I=1,NO
      R=s_rndm(0)
      N1=N1+1
      S=R
      IF(A.LE.1.E-34) GO TO 2
    1 A=S*A
      T=-ALOG(A)/SIGM1
      RETURN
    2 T1=ALOG(A)
      N2=NO-N1
      DO 4 I=1,N2
C      IX=IY
    4 T1=T1+ALOG(R)
      T=-T1/SIGM1
      RETURN
      END

C ___________________________________________________________________
C      SUBROUTINE RIONIZ(Z,E,U,T)
C ___________________________________________________________________
C      COMMON/VLZ/ALZ
C      AZ=ALZ
C      ALFA=104.946/((5.2096-0.33333*AZ)*Z)
C      DE=ALFA*(20.2+3.0*(ALOG(E)+0.67139)-2.*AZ)*T
C      U=(E-DE)/E
C      RETURN
C      END

C ______________________________________________________________________
C      REAL FUNCTION TCSPHP(X)
C ______________________________________________________________________
C     
C  Program, iz~iClqBA}A pxlHOTO CE~EHiE HA BzAjMOdEjCTBiE  ( fOTOEfE
C  KOMpTxHOB EfEKT + PAvdAHE HA ElEKTPOH-pOziTPOHHi dBOjKi )  NA fOTOH
C  B OlOBO HA EdHA PAdiAciOHHA EdiHicA
C   -X- pxlHATA Energy HA fOTOHA ( MEV )
C     
C      IF(5000.0.GT.X) GO TO 1
C      TCSPHP=X*1.2E-06+0.762
C      RETURN
C    1 IF(2000.0.GT.X) GO TO 2
C      TCSPHP=X*4.33E-06+0.7463
C      RETURN
C    2 IF(300.0.GE.X) GO TO 3
C      TCSPHP=0.5994+X*(0.468556E-03-X*(0.650237E-06-X*(0.422104E-09-X*
C     *0.973374E-13)))
C      RETURN
C    3 IF(60.0.GT.X) GO TO 4
C      TCSPHP=0.351425+X*(0.44749E-02-X*(0.27283E-04-X*(0.813791E-07-X*
C     *0.919342E-10)))
C      RETURN
C    4 IF(10.0.GE.X) GO TO 5
C      TCSPHP=0.20398+X*(0.012108-X*(0.164176E-03-X*0.918497E-06))
C      RETURN
C    5 IF(5.0.GT.X) GO TO 6
C      TCSPHP=0.261978-X*(0.452396E-02-X*(0.159871E-02-X*0.666096E-04))
C     RETURN
C    6 IF(1.5.GT.X) GO TO 7
C      TCSPHP=0.720468-X*(0.480939-X*(0.190856-X*(0.0337602-X*
C     *0.224457E-02)))
C      RETURN
C    7 TCSPHP=0.0
C      WRITE(3,8)
C    8 FORMAT(//10X,'APguMEHTA E izBxH OblACTTA HA pPilOvEHiE HA TCSPHE')
C      RETURN
C      END
