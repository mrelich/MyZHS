C ---------------------------------------------------------------------
C Thin the same shower with different thinning levels simultaneously
C and also perform a hybrid simulation (very useful to check accuracy
C of thinning levels)
C ---------------------------------------------------------------------

C #############################################################
c   IF A -> B+C
c   IF any B,C gt thin max, always keep them (same if both lt min)
C if A gt max, then keep b or c seperately
c if A, B, and C between min and max, then keep one of B and C (one test)
c otherwise, keep them seperately
C #############################################################

C Last update-> Jaime Alvarez-Muniz Jan 17 2008
C Last thinning level equivalent to hybrid simulation
C -------------------------------------------------------------
      SUBROUTINE THINNING
     #(ITYPE2,E2_AV,T2,ITYPE3,E3_AV,T3,Wgt2,Wgt3,KEEP2,KEEP3,Wgt0)

      PARAMETER (NTHIN_MAX=10)

      COMMON / THIN / THIN_MAX(NTHIN_MAX), THIN_MIN(NTHIN_MAX)
      COMMON / HYBRID / E_HYBRID
      COMMON / CONSERVATION / E_CONS

      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
      COMMON / PHOTON / TRP

      COMMON / THIN_LEVELS / NTHIN

      DOUBLE PRECISION E_tot,InvWgt2(NTHIN_MAX),InvWgt3(NTHIN_MAX),
     #E2_AV,E3_AV
      DOUBLE PRECISION Wgt0(NTHIN_MAX)
      DOUBLE PRECISION Wgt2(NTHIN_MAX),Wgt3(NTHIN_MAX)

      INTEGER KEEP2(NTHIN_MAX+1),KEEP3(NTHIN_MAX+1)
c Keep in mind that 0='neither', 1='TWO', 2='THREE',3='both'  

c Loop over number of thinning levels
      if (NTHIN.eq.0) goto 333
      do ith=1,NTHIN

c If mother particle has weight 0 then daughters have weight 0
      if (Wgt0(ith).eq.0) then
         keep2(ith)=0
         keep3(ith)=0
         Wgt2(ith)=0.
         Wgt3(ith)=0.
         goto 10 
      end if
 
      IF (E2_AV .LE. 0.) THEN
        KEEP2(ith) = 0

        IF (E3_AV .LE. 0.) THEN
         KEEP3(ith) = 0
        ELSE
         KEEP3(ith) = 1
         INVWGT3(ith) = 1.
         WGT3(ith)=1./INVWGT3(ith)
        END IF

      ELSE IF (E3_AV .LE. 0.) THEN

        KEEP2(ith) = 1
        InvWgt2(ith) = 1.
        Wgt2(ith)=1./InvWgt2(ith)
        KEEP3(ith) = 0

      ELSE

        E_TOT = E2_AV + E3_AV

        IF (E_TOT .LT. THIN_MIN(ith)) THEN
          KEEP2(ith) = 1
          KEEP3(ith) = 1
          InvWgt2(ith) = 1.
          InvWgt3(ith) = 1.
          Wgt2(ith)=1./InvWgt2(ith)
          Wgt3(ith)=1./InvWgt3(ith)
        ELSE

          IF (E_TOT .GT. THIN_MAX(ith)) THEN

            IF (E3_AV .GT. THIN_MAX(ith)) THEN
              KEEP3(ith) = 1
              InvWgt3(ith) = 1.
              Wgt3(ith)=1./InvWgt3(ith)
            ELSE
              InvWgt3(ith) = MAX(E3_AV, THIN_MIN(ith))/THIN_MAX(ith)
              Wgt3(ith)=1./InvWgt3(ith)

              IF(s_rndm(0) .LT. InvWgt3(ith)) THEN
                KEEP3(ith) = 1
              ELSE
                KEEP3(ith) = 0
              END IF

            END IF

            IF (E2_AV .GT. THIN_MAX(ith)) THEN
              KEEP2(ith) = 1
              InvWgt2(ith) = 1.
              Wgt2(ith)=1./InvWgt2(ith)
            ELSE
              InvWgt2(ith) = MAX(E2_AV, THIN_MIN(ith))/THIN_MAX(ith)
              Wgt2(ith)=1./InvWgt2(ith)
              IF(s_rndm(0) .LT. InvWgt2(ith)) THEN
                KEEP2(ith) = 1
              ELSE
                KEEP2(ith) = 0
              END IF
            END IF

          ELSE

          IF ((E2_AV .GT. THIN_MIN(ith)) .AND. 
     #(E3_AV .GT. THIN_MIN(ith))) THEN
            InvWgt2(ith) = E2_AV/E_TOT
            Wgt2(ith)=1./InvWgt2(ith)

            IF (s_rndm(0) .LT. InvWgt2(ith)) THEN
              KEEP2(ith) = 1
              KEEP3(ith) = 0
            ELSE
              KEEP2(ith) = 0
              KEEP3(ith) = 1
              InvWgt3(ith) = E3_AV/E_TOT
              Wgt3(ith)=1./InvWgt3(ith)
            END IF

          ELSE
            InvWgt2(ith) = MAX(E2_AV, THIN_MIN(ith))/E_TOT
            InvWgt3(ith) = MAX(E3_AV,THIN_MIN(ith))/E_TOT
            Wgt2(ith)=1./InvWgt2(ith)
            Wgt3(ith)=1./InvWgt3(ith)

            IF (s_rndm(0) .LT. InvWgt2(ith)) THEN
              KEEP2(ith) = 1
            ELSE
              KEEP2(ith) = 0
            END IF

            IF (s_rndm(0) .LT. InvWgt3(ith)) THEN
              KEEP3(ith) = 1
            ELSE
              KEEP3(ith) = 0
            END IF

          END IF

          END IF
        END IF
      END IF

10    continue

      end do  ! Ends loop over thinning levels

333   continue
c #############      
c Hybrid shower
c #############      
c If mother particle not kept, do not keep daughters      
      if (Wgt0(NTHIN+1).eq.0) then
       KEEP2(NTHIN+1)=0
       Wgt2(NTHIN+1)=0.
       KEEP3(NTHIN+1)=0
       Wgt3(NTHIN+1)=0.
       goto 444
      end if

c Stop following particles below the hybrid threshold      
      if(e2_av.lt.e_hybrid) then
          KEEP2(NTHIN+1)=0
          Wgt2(NTHIN+1)=0.
          if (itype2.eq.1) then  ! e- or e+
             E2_TOT=E2_AV+CH_THR
             CALL SIZE_GREISEN(ITYPE2,E2_TOT,T2)
             CALL TRACKLENGTHS(E2_TOT)
             E_CONS=E_CONS+E2_TOT
          else  ! Photon
             E2_TOT=E2_AV+TRP
             CALL SIZE_GREISEN(ITYPE2,E2_TOT,T2)
             CALL TRACKLENGTHS(E2_TOT)
             E_CONS=E_CONS+E2_TOT
          end if
c          write(*,*) '2 ',itype2,E2_TOT,T2
      else
          KEEP2(NTHIN+1)=1
          Wgt2(NTHIN+1)=1.
      end if
 
      if(e3_av.lt.e_hybrid) then
          KEEP3(NTHIN+1)=0
          Wgt3(NTHIN+1)=0.
          if (itype3.eq.1) then  ! e- or e+
             E3_TOT=E3_AV+CH_THR
             CALL SIZE_GREISEN(ITYPE3,E3_TOT,T3)
             CALL TRACKLENGTHS(E3_TOT)
             E_CONS=E_CONS+E3_TOT
          else  ! Photon
             E3_TOT=E3_AV+TRP
             CALL SIZE_GREISEN(ITYPE3,E3_TOT,T3)
             CALL TRACKLENGTHS(E3_TOT)
             E_CONS=E_CONS+E3_TOT
          end if
c          write(*,*) '3 ',itype3,E3_TOT,T3
      else
          KEEP3(NTHIN+1)=1
          Wgt3(NTHIN+1)=1.
      end if
c #############      

444   continue
      RETURN
      END
