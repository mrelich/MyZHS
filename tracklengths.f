      SUBROUTINE TRACKLENGTHS (E_MEV)
C--------------------------------------------------------
C
C...Compute the total tracklength of a shower 
C   generated by a gamma, e- or e+ of initial 
C   energy E (MeV) 
C
C--------------------------------------------------------
      COMMON / HYBRID_TRACK / TOT_TCK_HYBRID
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS

      COMMON/ECRIT/E_CRIT
      COMMON/BLMASS/e_MASS,e2_MAS

      DATA XLOG/5./ 

c Doesn't count track of particles below the critical energy 
      IF (E_MEV .LE. E_CRIT) RETURN
      
c Tracklengths as a function of particle threshold
c      MXJETH=MIN(INT(XLOG*LOG10((E_MEV-e_MASS)/BR_THR))+1,10)

c Total tracklength = 5.919 m/GeV
c Convert it to radiation lengths      
      TOTAL_TRACK=(5.919*100.*RHO/T0)*(E_MEV*1.e-3)

      TOT_TCK_HYBRID = TOT_TCK_HYBRID + TOTAL_TRACK 

      RETURN
      END