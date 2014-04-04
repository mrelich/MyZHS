
c      DOUBLE PRECISION FUNCTION PSRAN(DUM)
C***********************************************************************
C
C  interface function to avoid conflicting calling conventions for RNDM
C
C***********************************************************************
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      SAVE

c      REAL S_RNDM

c 10   continue
c        Xi = S_RNDM(Idummy)
c      if((Xi.le.0.D0).or.(Xi.ge.1.D0)) goto 10
c      PSRAN = Xi

c      END


c      FUNCTION RNDM(Xdum)
C***********************************************************************
C
C  interface function to avoid conflicting calling conventions for RNDM
C
C***********************************************************************
c      RNDM = S_RNDM(Idummy)
c      END


      FUNCTION S_RNDM(IDUMMY)
C...Generator  from the LUND montecarlo
C...Purpose: to generate random numbers uniformly distributed between
C...0 and 1, excluding the endpoints.
      COMMON/S_LUDATR/MRLU(6),RRLU(100)
      SAVE /S_LUDATR/
      EQUIVALENCE (MRLU1,MRLU(1)),(MRLU2,MRLU(2)),(MRLU3,MRLU(3)),
     &(MRLU4,MRLU(4)),(MRLU5,MRLU(5)),(MRLU6,MRLU(6)),
     &(RRLU98,RRLU(98)),(RRLU99,RRLU(99)),(RRLU00,RRLU(100))
 
      IDUMMY=IDUMMY*1
C...Initialize generation from given seed.
      IF(MRLU2.EQ.0) THEN
        IF (MRLU1 .EQ. 0)  MRLU1 = 19780512    ! initial seed
        write(6,'(1x,a,i12)') 'S_RNDM: initialization seed: ',MRLU1
        IJ=MOD(MRLU1/30082,31329)
        KL=MOD(MRLU1,30082)
        I=MOD(IJ/177,177)+2
        J=MOD(IJ,177)+2
        K=MOD(KL/169,178)+1
        L=MOD(KL,169)
        DO 110 II=1,97
        S=0.
        T=0.5
        DO 100 JJ=1,24
        M=MOD(MOD(I*J,179)*K,179)
        I=J
        J=K
        K=M
        L=MOD(53*L+1,169)
        IF(MOD(L*M,64).GE.32) S=S+T
        T=0.5*T
  100   CONTINUE
        RRLU(II)=S
  110   CONTINUE
        TWOM24=1.
        DO 120 I24=1,24
        TWOM24=0.5*TWOM24
  120   CONTINUE
        RRLU98=362436.*TWOM24
        RRLU99=7654321.*TWOM24
        RRLU00=16777213.*TWOM24
        MRLU2=1
        MRLU3=0
        MRLU4=97
        MRLU5=33
      ENDIF
 
C...Generate next random number.
  130 RUNI=RRLU(MRLU4)-RRLU(MRLU5)
      IF(RUNI.LT.0.) RUNI=RUNI+1.
      RRLU(MRLU4)=RUNI
      MRLU4=MRLU4-1
      IF(MRLU4.EQ.0) MRLU4=97
      MRLU5=MRLU5-1
      IF(MRLU5.EQ.0) MRLU5=97
      RRLU98=RRLU98-RRLU99
      IF(RRLU98.LT.0.) RRLU98=RRLU98+RRLU00
      RUNI=RUNI-RRLU98
      IF(RUNI.LT.0.) RUNI=RUNI+1.
      IF(RUNI.LE.0.OR.RUNI.GE.1.) GOTO 130
 
C...Update counters. Random number to output.
      MRLU3=MRLU3+1
      IF(MRLU3.EQ.1000000000) THEN
        MRLU2=MRLU2+1
        MRLU3=0
      ENDIF
      S_RNDM=RUNI
      RETURN
      END

