      FUNCTION GAUSS(F,A,B,EPS) 
      save
      LOGICAL MFLAG,RFLAG   
      EXTERNAL F    
      DIMENSION W(12),X(12) 
C   
C     ******************************************************************    
C   
C     ADAPTIVE GAUSSIAN QUADRATURE. 
C   
C     GAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF    
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER   
C     EPS.  
C   
C     ******************************************************************    
C   
      DATA W    
     */1.01228536E-01, 2.22381034E-01, 3.13706646E-01,  
     * 3.62683783E-01, 2.71524594E-02, 6.22535239E-02,  
     * 9.51585117E-02, 1.24628971E-01, 1.49595989E-01,  
     * 1.69156519E-01, 1.82603415E-01, 1.89450610E-01/  
    
      DATA X    
     */9.60289856E-01, 7.96666477E-01, 5.25532410E-01,  
     * 1.83434642E-01, 9.89400935E-01, 9.44575023E-01,  
     * 8.65631202E-01, 7.55404408E-01, 6.17876244E-01,  
     * 4.58016778E-01, 2.81603551E-01, 9.50125098E-02/  
C   
C     ******************************************************************    
C   
C  START.   
      GAUSS=0.  
      IF(B.EQ.A) RETURN 
      CONST=0.005/(B-A) 
      BB=A  
C   
C  COMPUTATIONAL LOOP.  
    1 AA=BB 
      BB=B  
    2    C1=0.5*(BB+AA) 
         C2=0.5*(BB-AA) 
         S8=0.  
         DO 3 I=1,4 
            U=C2*X(I)   
            S8=S8+W(I)*(F(C1+U)+F(C1-U))    
    3    CONTINUE   
         S8=C2*S8   
         S16=0. 
         DO 4 I=5,12    
            U=C2*X(I)   
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE   
         S16=C2*S16 
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5   
         BB=C1  
         IF( 1.+ABS(CONST*C2) .NE. 1. ) GO TO 2 
      GAUSS=0.  
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG)  
      IF(MFLAG) THEN    
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)  
         ELSE   
            WRITE(LGFILE,6) 
         ENDIF  
      ENDIF 
      IF(.NOT. RFLAG) CALL ABEND    
      RETURN    
    5 GAUSS=GAUSS+S16   
      IF(BB.NE.B) GO TO 1   
      RETURN    
C   
    6 FORMAT( 4X, 'FUNCTION GAUSS ... TOO HIGH ACCURACY REQUIRED')  
      END   
