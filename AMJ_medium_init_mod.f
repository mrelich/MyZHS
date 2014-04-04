      subroutine Medium_Properties
c     program init_only

      parameter(pi=3.1415926)
      parameter(alpha=1./137.)
 
      parameter (Maximum_Z=100)
      parameter (iA_PER_MOLECULE=2)
      parameter (MAX_N_MOLECULES=6)
      parameter (N_ENERGIES=400)
      INTEGER Z_IN_MOLECULES(Max_N_MOLECULES*iA_PER_MOLECULE)
      DIMENSION COMPOSITION(Max_N_MOLECULES*iA_PER_MOLECULE)
      DIMENSION ABUNDANCE(Max_N_MOLECULES*iA_PER_MOLECULE)

      INTEGER MOL_ENTRY_METHOD, PHASE, N_MOLECULES, Z1,Z2
      REAL MEAN_MASS, TRP
      
      COMMON / RADIATION / RADLOWZ(4), RADPRIM(4)
      COMMON / BCUT / BR_THR   !  MINIMUM PHOTON ENERGY IN MEV
      COMMON / MEDIUM_COMPOSITION / Max_Z, FRAC_NUM(Maximum_Z),
     #ZATOM(Maximum_Z), FRAC_WEIGHT(Maximum_Z),
     #FC_Zi(Maximum_Z), Xi_Zi(Maximum_Z),
     #FZ4_50(Maximum_Z),
     #AATOM(Maximum_Z),  XI_ion(Maximum_Z), PHASE
     
c     These cross-sections are also passed in to ZHS_multimedia. _E = energy, _x = xsection.     
      COMMON/ BREMSS_XSECTION / BR_E(400), BR_X(400)
      COMMON / PAIR_XSECTION / PP_E(400), PP_X(400)
      
c     This common gets passed into 'ZHS_multimedia'. It contains all medium-specific properties, other than the cross-sections (which are contained in BREMSS_XSECTION and PAIR_XSECTION Note that 'Z' is mean Z weighted by number for most of the routine, but is set equal to Z_Wgt (ie mean Z by weight) before the routine finished.
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR, Z_Wgt, FACTOR_MS 
      
      COMMON / THIS_ENERGY / ENERGY

      COMMON/ZFUNCT/Z_B,Z_F,Z_F50,Z_AB,Z_A,Z_T,
     #Z_V,Z_G,Z_U,Z_P,Z_VP50,Z_UG50
     #,Z_VG50,Z_UP50, Z_W
 
        DATA RADLOWZ/5.310, 4.790, 4.740, 4.710/
        DATA RADPRIM/6.144, 5.621, 5.805, 5.924/
        CHARACTER*16 Medium_File
        DATA E_MS/ 21.0 /
      
c      see EGS4 manual eqn 2.7.17 (p31). Accurate to 4 digits for Z < Uranium
        FCZFUN(ALPZ2)=ALPZ2*(1./(1.+ALPZ2)+
     # 0.20206+ALPZ2*(-0.0369+ALPZ2*(0.0083-0.002*ALPZ2)))

      Max_Z=Maximum_Z

C ##################################################################################
c      MEDIUM PROPERTIES. THE FOLLOWING MUST BE MANUALLY CHANGED IN THE INPUT FILE!
C ##################################################################################

c   file containing the properties of the medium
!      PRINT *, 'Which medium file?'
      READ(*,*) Medium_file



      OPEN(33,FILE=Medium_File, STATUS='UNKNOWN')
      write(*,*) Medium_File
      OPEN(29,FILE='Medium_Init.out',STATUS='UNKNOWN')
      OPEN(30,FILE='PP_xsections.out',STATUS='UNKNOWN')
      OPEN(31,FILE='BR_xsections.out',STATUS='UNKNOWN')
      OPEN(32,FILE='Medium_Init_info.out',STATUS='UNKNOWN')
      write(29,*) Medium_file


c Refractive index of the medium
      READ(33,*) REF_N
      write(29,*) 'REF_N= ',REF_N
      
c     Density in gcm^-3      
      READ(33,*) RHO
      write(29,*) 'RHO= ',RHO
      
c     0=solid,1=liquid,2=gas. Currently, solids and liquids are equivalent.
      READ(33,*) PHASE
      write(29,*) 'PHASE= ',PHASE
      
c 1=the molecules are entered as % by number, and the program calculates % by weight and frac_number.
c 0=the molecules are entered as % by weight, and the program calculates fraction by number
c Note that fraction by number is NOT % by number!!!!!!!!
      READ(33,*) MOL_ENTRY_METHOD
      write(29,*) 'MOL_ENTRY_METHOD= ',MOL_ENTRY_METHOD
      
c The number of molecules present in the material. If this is greater than the current value
c of 'max_nmolecules' (see the parameters), then please change thus also
      READ(33,*) N_MOLECULES
      write(29,*) 'N_MOLECULES= ',N_MOLECULES
    
c #####################################
c HOW TO ENTER THE MEDIUM COMPOSITION    
c #####################################
      
c The entry method for molecule n is as follows:
c #1  Choose one element from the molecule.
c #2  Enter in the Z of element 1 (Z1) into 'Z_IN_MOLECULES(2*n-1)', ie Z_IN_MOLECULES(2*n-1)=Z1
c #3 Enter in the number of atoms of element 1 (N1) per molecule, ie COMPOSITION(2*n-1)= N1
c #4 Repeat #2 and #3 for the second element in the molecule, using index (2*n).
c #5 Enter in the % abundance of the molecule by weight / number into ABUNDANCE(i)
c MAKE SURE THE ABUNDANCE ENTERED CORRESPONDS TO THE VALUE OF 'MOL_ENTRY_METHOD' USED!


      DO I=1, N_MOLECULES
        READ(33,*) Z_IN_MOLECULES(2*I-1), Z_IN_MOLECULES(2*I),
     #COMPOSITION(2*I-1), COMPOSITION(2*i),ABUNDANCE(I)
      END DO

      write(29,*) (Z_IN_MOLECULES(J),J=1,2*N_MOLECULES)
      write(29,*) (COMPOSITION(J),J=1,2*N_MOLECULES)
      write(29,*) (ABUNDANCE(j),J=1,N_MOLECULES)

c      DO I=N_MOLECULES,MAX_N_MOLECULES
c        ABUNDANCE(I)=0
c      END DO

c ################################################################################
c
c     END OF MANUAL ENTRY. THE PROGRAM CALCULATES THE REST
c
c ################################################################################	


c max energy in MeV - 1.e15 should cover it (=1e21 eV!!!!)
      E_MAX=1.e18
c     1= 'use a log-linear scale for calculating energies at which to parameterise the x-sections'. Anything else= 'use a linear scale'.
      ILOG=1
   
c Mean molecular mass=Sum((Fraction%/100.)*Molecular mass)=66.41 g/mol
c Calculated doing:  FRAC_NUM(i)=(Fraction_i/100./A_i)*Mean molecular mass	

c Following data from slac-r-265-ch02b.pdf p10 of 23

c Atomic weights in AMUs
      AATOM(1)=1.00797
      AATOM(2)=4.00260
      AATOM(3)=6.939
      AATOM(4)=9.01220
      AATOM(5)=10.811
      AATOM(6)=12.01115
      AATOM(7)=14.00670
      AATOM(8)=15.9994
      AATOM(9)=18.9984
      AATOM(10)=20.183
      AATOM(11)=22.9898
      AATOM(12)=24.312
      AATOM(13)=26.9815
      AATOM(14)=28.088
      AATOM(15)=30.9738
      AATOM(16)=32.064
      AATOM(17)=35.453
      AATOM(18)=39.948
      AATOM(19)=39.102
      AATOM(20)=40.08
      AATOM(21)=44.956
      AATOM(22)=47.900
      AATOM(23)=50.942
      AATOM(24)=51.998
      AATOM(25)=54.938
      AATOM(26)=55.847
      AATOM(27)=58.9332
      AATOM(28)=58.71
      AATOM(29)=63.54
      AATOM(30)=65.37
      AATOM(31)=69.72
      AATOM(32)=72.59
      AATOM(33)=74.9216
      AATOM(34)=78.96
      AATOM(35)=79.808
      AATOM(36)=83.8
      AATOM(37)=85.47
      AATOM(38)=87.62
      AATOM(39)=88.905
      AATOM(40)=91.22
      AATOM(41)=92.91
      AATOM(42)=95.94
      AATOM(43)=99.00
      AATOM(44)=101.07
      AATOM(45)=102.91
      AATOM(46)=106.40
      AATOM(47)=107.87
      AATOM(48)=112.40
      AATOM(49)=114.82
      AATOM(50)=118.69
      AATOM(51)=121.75
      AATOM(52)=127.60
      AATOM(53)=126.90
      AATOM(54)=131.30
      AATOM(55)=132.91
      AATOM(56)=137.34
      AATOM(57)=138.91
      AATOM(58)=140.12
      AATOM(59)=140.91
      AATOM(60)=144.24
      AATOM(61)=147.00
      AATOM(62)=150.35
      AATOM(63)=151.98
      AATOM(64)=157.25
      AATOM(65)=158.92
      AATOM(66)=162.50
      AATOM(67)=164.93
      AATOM(68)=167.26
      AATOM(69)=168.93
      AATOM(70)=173.04
      AATOM(71)=174.97
      AATOM(72)=178.49
      AATOM(73)=180.95
      AATOM(74)=183.85
      AATOM(75)=186.20
      AATOM(76)=190.20
      AATOM(77)=192.20
      AATOM(78)=195.08
      AATOM(79)=196.99
      AATOM(80)=200.59
      AATOM(81)=204.37
      AATOM(82)=207.19
      AATOM(83)=208.98
      AATOM(84)=210.00
      AATOM(85)=210.00
      AATOM(86)=222.00
      AATOM(87)=223.00
      AATOM(88)=226.00
      AATOM(89)=227.00
      AATOM(90)=232.04
      AATOM(91)=231.00
      AATOM(92)=238.03
      AATOM(93)=237.00
      AATOM(94)=242.00
      AATOM(95)=243.00
      AATOM(96)=247.00
      AATOM(97)=247.00
      AATOM(98)=248.00
      AATOM(99)=254.00
      AATOM(100)=253.00

c Ionisation potentials in eV
      XI_ION(1)=19.2
      XI_ION(2)=41.8
      XI_ION(3)=40.0
      XI_ION(4)=63.7
      XI_ION(5)=76.0
      XI_ION(6)=78.0
      XI_ION(7)=82.0
      XI_ION(8)=95.0
      XI_ION(9)=115.0
      XI_ION(10)=137.0
      XI_ION(11)=149.0
      XI_ION(12)=156.0
      XI_ION(13)=166.0
      XI_ION(14)=173.0
      XI_ION(15)=173.0
      XI_ION(16)=180.0
      XI_ION(17)=174.0
      XI_ION(18)=188.0
      XI_ION(19)=190.0
      XI_ION(20)=191.0
      XI_ION(21)=216.0
      XI_ION(22)=233.0
      XI_ION(23)=245.0
      XI_ION(24)=257.0
      XI_ION(25)=272.0
      XI_ION(26)=286.0
      XI_ION(27)=297.0
      XI_ION(28)=311.0
      XI_ION(29)=322.0
      XI_ION(30)=330.0
      XI_ION(31)=334.0
      XI_ION(32)=350.0
      XI_ION(33)=347.0
      XI_ION(34)=348.0
      XI_ION(35)=357.0
      XI_ION(36)=352.0
      XI_ION(37)=363.0
      XI_ION(38)=366.0
      XI_ION(39)=379.0
      XI_ION(40)=393.0
      XI_ION(41)=417.0
      XI_ION(42)=424.0
      XI_ION(43)=428.0
      XI_ION(44)=441.0
      XI_ION(45)=449.0
      XI_ION(46)=470.0
      XI_ION(47)=470.0
      XI_ION(48)=469.0
      XI_ION(49)=488.0
      XI_ION(50)=488.0
      XI_ION(51)=487.0
      XI_ION(52)=485.0
      XI_ION(53)=491.0
      XI_ION(54)=482.0
      XI_ION(55)=488.0
      XI_ION(56)=491.0
      XI_ION(57)=501.0
      XI_ION(58)=523.0
      XI_ION(59)=535.0
      XI_ION(60)=546.0
      XI_ION(61)=560.0
      XI_ION(62)=574.0
      XI_ION(63)=580.0
      XI_ION(64)=591.0
      XI_ION(65)=614.0
      XI_ION(66)=628.0
      XI_ION(67)=650.0
      XI_ION(68)=658.0
      XI_ION(69)=674.0
      XI_ION(70)=684.0
      XI_ION(71)=694.0
      XI_ION(72)=705.0
      XI_ION(73)=718.0
      XI_ION(74)=727.0
      XI_ION(75)=736.0
      XI_ION(76)=746.0
      XI_ION(77)=757.0
      XI_ION(78)=790.0
      XI_ION(79)=790.0
      XI_ION(80)=800.0
      XI_ION(81)=810.0
      XI_ION(82)=823.0
      XI_ION(83)=823.0
      XI_ION(84)=830.0
      XI_ION(85)=825.0
      XI_ION(86)=794.0
      XI_ION(87)=827.0
      XI_ION(88)=826.0
      XI_ION(89)=841.0
      XI_ION(90)=847.0
      XI_ION(91)=878.0
      XI_ION(92)=890.0
      XI_ION(93)=902.0
      XI_ION(94)=921.0
      XI_ION(95)=934.0
      XI_ION(96)=939.0
      XI_ION(97)=952.0
      XI_ION(98)=966.0
      XI_ION(99)=980.0
      XI_ION(100)=994.0

c End 'data taken from...'

c Initialises variables
      do I=1,Max_Z
        FRAC_NUM(I)=0.0
        FRAC_WEIGHT(I)=0.0
      end do
      TOTAL_WEIGHT=0.0
      TOTAL_NUM=0.0
      MEAN_MASS=0.0
      T0_INV=0.0

c if entered by number - there is currently no independent check on this section!
      IF (MOL_ENTRY_METHOD .EQ. 1) THEN
        DO I=1,N_MOLECULES
          Z1=Z_IN_MOLECULES(2*I-1)
          Z2=Z_IN_MOLECULES(2*I)
          FRAC_NUM(Z1) = FRAC_NUM(Z1)+COMPOSITION(2*I-1)*ABUNDANCE(I)
          FRAC_NUM(Z2) = FRAC_NUM(Z2)+COMPOSITION(2*I)*ABUNDANCE(I)
          TOTAL_WEIGHT=TOTAL_WEIGHT+AATOM(Z1)*COMPOSITION(2*I-1)
          TOTAL_WEIGHT=TOTAL_WEIGHT+AATOM(Z2)*COMPOSITION(2*I)
        END DO
        DO I=1,Max_Z
          FRAC_WEIGHT(I)=FRAC_NUM(I)*AATOM(I)/TOTAL_WEIGHT
        END DO
      END IF

c if entered by weight
      IF (MOL_ENTRY_METHOD .EQ. 0) THEN
        DO I=1,N_MOLECULES
          Z1=Z_IN_MOLECULES(2*I-1)
          Z2=Z_IN_MOLECULES(2*I)
          WEIGHT1=COMPOSITION(2*I-1)*AATOM(Z1)
          WEIGHT2=COMPOSITION(2*I)*AATOM(Z2)
          FRAC_WEIGHT(Z1) = FRAC_WEIGHT(Z1)+
     #WEIGHT1*ABUNDANCE(I)/(WEIGHT1+WEIGHT2)
          FRAC_WEIGHT(Z2) = FRAC_WEIGHT(Z2)+
     #WEIGHT2*ABUNDANCE(I)/(WEIGHT1+WEIGHT2)
          MEAN_MASS=MEAN_MASS+(WEIGHT1+WEIGHT2)*ABUNDANCE(I)
        END DO

c   fraction by number = fraction by weight *(mean molecular weight/atomic weight)
      DO I=1,Max_Z
          FRAC_NUM(I)=FRAC_WEIGHT(I)*MEAN_MASS/AATOM(I)
        END DO
      END IF   
      
      DO I=1,MAX_Z
        write(32,*) "Frac num ", I, " = ", FRAC_NUM(I)
        write(32,*) "Frac weight ", I, " = ", FRAC_WEIGHT(I)
      END DO
      
c  Z_W is the mean Z by weight - it is not calculated anywhere else!	
        Z=0.
        Z_Wgt=0.
        Z_W=0.
        X_M=0.
        Z_S=0.
        Z_E=0.
        Z_X=0.
        Z_T=0.
        Z_A=0.
        Z_B=0.
        Z_F=0.
        Z_P=0.
        Z_G=0.
        Z_U=0.
        Z_V=0.
        Z_AB=0.
c -----------------------

        psi_MS=1.   ! EGS manual p 76. Taken as a fudge factor in mult. scatt. 

c Avogadros number
        XN_A=6.02252e23

c Classical electron radius in cm
        R0_e=2.817940e-13
 
c Electron mass in MeV
        e_MASS=0.51099906

c AKO and EKE are used in photoelectric effect subroutine
c which is not taken into account in ZHS. 
      AKO=2.907E-3
      EKE=.88E-1
      AP=.1E-1

c For eqns for Z_x, see EGS4 manual p23, table 2.6.1 (contd)
       do i=1,Max_Z
          IF (FRAC_NUM(i) .NE. 0.0) THEN
          ZI=REAL(I)
          Z=Z+ZI*FRAC_NUM(i)
          Z_W=Z_W+ZI*FRAC_WEIGHT(I)
          X_M=X_M+aatom(i)*FRAC_NUM(i)
          Z_S=Z_S+ZI*(ZI+psi_MS)*FRAC_NUM(i)
          Z_E=Z_E+ZI*(ZI+psi_MS)*FRAC_NUM(i)*
     #log(ZI**(-2./3.))
          Z_X=Z_X+ZI*(ZI+psi_MS)*FRAC_NUM(i)*
     #log(1.+3.34*((1./137.)*ZI)**2.)

c Radiation length in g cm^-2 - T0_INV is the inverse, in cm^2/g. Individually calculated according
c to Dahl - see pdg 2004, 10:26. This is about 5% too low for Z=2 (helium) however.
      T0_INV=T0_INV+
     #FRAC_WEIGHT(I)*(ZI*(ZI+1)*LOG(287./SQRT(ZI)))
     #/(716.4*AATOM(I))
     
          IF (ZI .GT. 4.) THEN
             RADL=5.2157506-LOG(ZI)/3.                  ! LOG (184.15)
             RADP=7.0850643-0.66666667*LOG(ZI)          ! LOG (1194)
          ELSE
             RADL=RADLOWZ(int(ZI))
             RADP=RADPRIM(int(ZI))
          END IF

          FC_Zi(I)=FCZFUN( (ZI/137.)**2 )
          Xi_Zi(I)=RADP/(RADL-FC_Zi(I))
          ZI3LOG=LOG(ZI)/3.
          Z2CSI=FRAC_NUM(I)*ZI*(ZI+Xi_Zi(I))

          Z_T=Z_T+Z2CSI
          Z_B=Z_B-Z2CSI*ZI3LOG
          Z_F=Z_F+Z2CSI*FC_Zi(I)
          Z_AB=Z_AB+Z2CSI*RADL
          END IF
      end do

        Z_A=Z_T*log(183.)
        Z_P=Z_B/Z_A
        Z_G=Z_B/Z_T
        Z_U=(Z_B-Z_F)/Z_A
        Z_V=(Z_B-Z_F)/Z_T
        Z_Wgt=Z_W
c
c Coulomb correction term EGS4 (2.7.17)
c
      AZ=Z_W/137.
      FCZ=AZ*AZ*(1./(1.+AZ*AZ)+
     #           .20206-.0369*AZ*AZ+.0083*AZ**4.-.002*AZ**6.)
      ALZ3=ALOG(183./Z_W**.3333333)-FCZ
      ALZ=ALOG(Z_W)
c
c EGS4 (2.7.21)
c
      Xi_Z=(7.0850642-0.66666667*ALZ)/(5.2157506-ALZ/3.-FCZ)
      
     
      XN_e=XN_A*RHO*Z/X_M
      write(32,*) 'Z',Z
      write(32,*) 'X_M',X_M
      T0=1./T0_INV
      write(32,*) 'T0',T0
      write(32,*) 'N_electron',XN_E
      write(32,*) 'Z_S',Z_S
      write(32,*) 'Z_E',Z_E
      write(32,*) 'Z_X',Z_X
      write(32,*) 'Z_W',Z_W

      RADLEN=T0/RHO
c LPM in TeV - from Stanevs paper
      E_LPM=61.5*T0/RHO
c conversion to MeV
      E_LPM = E_LPM*1.e6
      write(32,*) 'E_LPM [MeV] = ',E_LPM
c once v below c/n, no cherenkov radiation produced
      CH_THR=e_MASS/sqrt(1.-1./(REF_N**2))

c ########################################
c WARNING - WARNING - WARNING - WARNING - 
c ########################################
c Next line overrides calculation of Cherenkov threshold
c to obtain electric fields for particles below Cher. threshold
c Fixes kinetic energy threshold to 100 keV
c Note: affects calculation of continuous energy losses. 
      CH_THR=e_MASS + 0.1
      write(32,*) 'Cher_Thresh [MeV] = ',CH_THR

C Bremsstrahlung threshold: cherenkov threshold - electron mass
      BR_THR=CH_THR-e_MASS

c Critical energy [MeV]
      E_CRIT=610./(Z_W+1.24)
      write(32,*) 'Critical energy [MeV] ',E_CRIT

c Moliere radius on 100: T0*E_MS (21 MeV)/E_CRIT
      R0=T0*E_MS/E_CRIT/100.
      write(32,*) 'Moliere radius /100 (g.cm^2)',R0

      BEE=log(T0/RHO)+log(6702.33*Z_S*RHO/X_M)+((Z_E-Z_X)/Z_S)
      write(32,*) 'BEE',bee

c Factor appearing in multiscattering. See EGS4 manual (2.14.24-25) and Note in between.
c	xkappa_cc=sqrt(T0/RHO)*(22.696/(180/pi))*sqrt(RHO*Z_S/X_M)
c	write(32,*) 'KAPPA_CC*SQRT(T0/RHO)',xkappa_cc
        FACTOR_MS=sqrt(T0/RHO)*(22.696/(180/pi))*sqrt(RHO*Z_S/X_M)
        write(32,*) 'FACTOR_MS= ', FACTOR_MS

c Ionization loss parameters
        do i=1,Max_Z
         xnum=xnum+FRAC_NUM(i)*I*(log(XI_ion(i)))
        end do

        XlogI_adj=xnum/Z
        write(32,*) 'I_adj',exp(XlogI_adj)

c Ionization potential/m_electron
        AI=exp(XlogI_adj)/(e_MASS*1.e6)
        write(32,*) 'AI',AI

c Plancks constant (not reduced) in eV*s
        xh=2.*pi*6.582122*1.e-22*1.e6

c Speed of light in cm/s
        c_light=3.e10

c Plasma frequency
        xnu_p=sqrt(xn_e*r0_e*c_light*c_light/pi)

c Plasma energy
        energy_p=xh*xnu_p
        write(32,*) 'Plasma energy (eV) ',energy_p

c C parameter - see EGS4 manual, p73, 2.13.21 & 2.13.22
        C=-1.*(2.*log(exp(XlogI_adj)/energy_p)+1.)
        write(32,*) 'C',C

c The remaining parameters can be read directly in EGS4 manual, p74
c NOTE: For some media, these have been experimentally obtained - the values presented here are rough! (compare measured values, p71-72, with the prescription below).
        MS=3.
        RM=MS
        IF (phase .EQ. 1 .OR. PHASE .EQ. 0) THEN

          IF (exp(XlogI_adj) .LT. 1e-4) THEN
            X1=2.0
            IF (-C .LT. 3.681) THEN
              X0=0.2
            ELSE
              X0=-0.362*C-1 
            END IF
          ELSE
             X1=3.0
            IF (-C .LT. 5.215) THEN
              X0=0.2
            ELSE
              X0=-0.326*C-1.5
            END IF
          END IF

        ELSE

          IF (-C .LT. 10.0) THEN
            X0=1.6
            X1=4.0
          END IF
          IF (-C .GE. 10.0 .AND. -C .LT. 10.5) THEN
            X0=1.7
            X1=4.0
          END IF
          IF (-C .GE. 10.5 .AND. -C .LT. 11.0) THEN
            X0=1.8
            X1=4.0
          END IF
          IF (-C .GE. 11.0 .AND. -C .LT. 11.5) THEN
            X0=1.9
            X1=4.0
          END IF
          IF (-C .GE. 11.5 .AND. -C .LT. 12.25) THEN
            X0=2.0
            X1=4.0
          END IF
          IF (-C .GE. 12.25 .AND. -C .LT. 13.804) THEN
            X0=2.0
            X1=5.0
          END IF
          IF (-C .GE. 13.804) THEN
            X0=-0.326*C-2.5
            X1=5.0
          END IF
        END IF

        A=(-C-2*log(10.)*X0)/(X1-X0)**3.

c  write these values down!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c FACT in MeV/cm
        FACT=2.*pi*R0_e*R0_e*e_MASS*XN_e
        write(32,*) 'Factor = 2*Pi*r_e^2*m_e*n_e = ',FACT,' MeV/cm'

        write(32,*) Z_AB,Z_F

c Radiation length as defined in table 2.6.1 in EGS4 manual
        X0_EGS=X_M/(RHO*4.*alpha*(Z_AB-Z_F))
        X0_EGS=X0_EGS/XN_A/R0_e/R0_e

c FACNOR Normalization of cross sections involving electrons as target
        FACNOR=FACT*X0_EGS
        write(32,*) 'X0_EGS = ',X0_EGS*RHO,' g/cm^2'
        write(32,*) 'FACNOR = ',FACNOR,' MeV'

c     Ensures the minimum energy is the lowest of the thresholds for PP & BREMSS
c      E_MIN=min(CH_THR,2.*e_MASS)

      TRP=0.5*(BR_THR+SQRT(BR_THR*(BR_THR+2.*e_MASS)))
      E_MIN=min(CH_THR, TRP)
c ### for loop in energy begins here
      EPSILN=0.1*E_MIN
      X0L=ALOG(E_MIN-EPSILN)
      X1L=ALOG(E_MAX)
      
      BR_MAX=-999999.      
      
      DO J=0,N_ENERGIES-1
          IF (ILOG .EQ. 1) THEN
             ENERGY=X0L+(X1L-X0L)*J/(N_ENERGIES-1)
             ENERGY=EXP(ENERGY)
          ELSE 
             ENERGY=E_MIN+(E_MAX-E_MIN)*J/(N_ENERGIES-1)
          END IF
          BR_E(J+1)=ENERGY
          PP_E(J+1)=ENERGY
c For energy smaller than 50 MeV cross section is different than for
c energy >50 MeV.
c energy >50 MeV Extreme relativistic Coulomb corrected cross section.
      IF (ENERGY .GT. 50.) THEN
        DO I=1,Max_Z
         FZ4_50(I)=4.*FC_Zi(I)
        END DO
        Z_F50=Z_F
        Z_UP50=Z_U
        Z_VG50=Z_V
      ELSE
        DO I=1,Max_Z
         FZ4_50(I)=0.
        END DO
        Z_F50=0.
        Z_UP50=Z_P
        Z_VG50=Z_G
      END IF

C ### CALCULATE XSECTIONS

      IF (Energy .lt. 2.*e_MASS) THEN
        PP_X(J+1) = 0.0
      ELSE
        PP_X(J+1)=max(0.,sumint_P(ENERGY))
      ENDIF
      write(30,*) PP_E(J+1), PP_X(J+1)
      IF (ENERGY .LT. CH_THR) THEN
        BR_X(J+1)=0.0
      ELSE
        BR_X(J+1)=max(0.0001,sumint_B(BR_E(J+1))) 
      END IF  
        IF(BR_X(J+1).GT.BR_MAX) THEN
         BR_MAX=BR_X(J+1)
         E_LPM_BREMSS=BR_E(J+1)
      END IF
      write(31,*) BR_E(J+1),BR_X(J+1)
C #### for loop in energy ends here
      END DO

c Fraction by number is not used in ZHS code only fraction by weight     
      Z=Z_W
  
      write(29,*) E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_W,FACTOR_MS
     
      CLOSE(33)
      CLOSE(29)
      CLOSE(30)
      CLOSE(31)
      CLOSE(32)
     
      RETURN
      END

c *********************************************************************
c *********************************************************************
c Total cross section for bremsstrahlung in the Moon, obtained integrating
c differential cross section in epsilon variable.
c The integration goes from some minimun photon energy to maximun photon
c energy. See EGS4 manual.
c **********************************************
c Can be used for other media. Only change medium parameters as well as
c Radiation length in cm and mean Z (mean by weight)
c **********************************************
c Based on the EGS4 manual
c *********************************************************************
      FUNCTION SUMINT_B(ENERGY)  !  TOTAL X-S AS FUNCTION OF KIN. E (MEV)
      save
c      COMMON / ENERGY / ENERGY   !  INITIAL ENERGY IN MEV
      COMMON / BCUT / BR_THR   !  MINIMUM PHOTON ENERGY IN MEV
      COMMON / RATIO / RATEM
      parameter (Maximum_Z=100)
      
      COMMON / MEDIUM_COMPOSITION / Max_Z, FRAC_NUM(Maximum_Z), 
     #ZATOM(Maximum_Z), FRAC_WEIGHT(Maximum_Z),
     #FC_Zi(Maximum_Z), Xi_Zi(Maximum_Z),
     #FZ4_50(Maximum_Z),
     # AATOM(Maximum_Z),  XI_ion(Maximum_Z), PHASE
      COMMON/ZFUNCT/Z_B,Z_F,Z_F50,Z_AB,Z_A,Z_T,Z_V,
     #Z_G,Z_U,Z_P,Z_VP50,Z_UG50
     #,Z_VG50,Z_UP50, Z_W
      COMMON / RADIATION / RADLOWZ(4), RADPRIM(4)
      DATA TOLER/1.E-5/

      EXTERNAL DXDGBR1         ! DIFFERENTIAL X-S IN Ephoton/Einitial

c Coulomb correction term page 31 EGS4
      FCZFUN(ALPZ2)=ALPZ2*(1./(1.+ALPZ2)+
     # 0.20206+ALPZ2*(-0.0369+ALPZ2*(0.0083-0.002*ALPZ2)))


        RATEM=0.511/ENERGY
        IF (RATEM .GE. 1.D0 .OR. ENERGY .LE. BR_THR) THEN
          SUMINT_B=0.D0
          RETURN
        END IF
c Integration limits. Integral is performed logarythmically.
        IF (ENERGY .LE. 1.E1) THEN
         GILMAX=LOG(1.-RATEM)
        ELSE
         GILMAX=-RATEM*(1.+RATEM*(0.5+RATEM*
     #(0.33333333+0.25*RATEM))) ! Taylor expansion of log(1.-me/E)
        END IF
        GILMIN=LOG(BR_THR/ENERGY)
c Integral in the fraction of energy carried by the emitted photon.

        SUMINT_B=GAUSS(DXDGBR1,GILMIN,GILMAX,TOLER)

      RETURN
      END


      FUNCTION SUMINT_P(ENERGY)  !  TOTAL X-S AS FUNCTION OF E (MEV)
      save
      COMMON / RATIO / RATEM
      parameter (Maximum_Z=100)
      
      COMMON / MEDIUM_COMPOSITION / Max_Z, FRAC_NUM(Maximum_Z), 
     #ZATOM(Maximum_Z), FRAC_WEIGHT(Maximum_Z),
     #FC_Zi(Maximum_Z), Xi_Zi(Maximum_Z),
     #FZ4_50(Maximum_Z),
     # AATOM(Maximum_Z),  XI_ion(Maximum_Z), PHASE
      COMMON/ZFUNCT/Z_B,Z_F,Z_F50,Z_AB,Z_A,Z_T,
     #Z_V,Z_G,Z_U,Z_P,Z_VP50,Z_UG50
     #,Z_VG50,Z_UP50, Z_W
C      DIMENSION RADLOWZ(4),RADPRIM(4)
      COMMON / RADIATION / RADLOWZ(4), RADPRIM(4)
      DATA TOLER/1.E-6/
      
      EXTERNAL DXDGPA         ! DIFFERENTIAL X-S IN E(e-,e+)/ENERGYial

c Coulomb correction term page 31 EGS4
      FCZFUN(ALPZ2)=ALPZ2*(1./(1.+ALPZ2)+
     # 0.20206+ALPZ2*(-0.0369+ALPZ2*(0.0083-0.002*ALPZ2)))

        RATEM=0.511/ENERGY
        IF (RATEM .GT. 0.5) THEN
         SUMINT_P=0.
         RETURN
        END IF
c Integration limits
        GILMAX=1.-RATEM
        GILMIN=RATEM
        SUMINT_P=GAUSS(DXDGPA,GILMIN,GILMAX,TOLER)

      RETURN
      END

C **************************************************************************
c Differential bremsstrahlung cross section including corrections for
c E>50 MeV, E<50 MeV and LPM corrections.
c **************************************************************************
      FUNCTION DXDGBR1(GILN)
      
      parameter (Maximum_Z=100)
      
      COMMON / MEDIUM_COMPOSITION / Max_Z, FRAC_NUM(Maximum_Z), 
     #ZATOM(Maximum_Z), FRAC_WEIGHT(Maximum_Z),
     #FC_Zi(Maximum_Z), Xi_Zi(Maximum_Z),
     #FZ4_50(Maximum_Z),
     # AATOM(Maximum_Z),  XI_ion(Maximum_Z), PHASE
      
      COMMON/ZFUNCT/Z_B,Z_F,Z_F50,Z_AB,Z_A,Z_T,Z_V,Z_G,
     #Z_U,Z_P,Z_VP50,Z_UG50
     #,Z_VG50,Z_UP50, Z_W
      
      COMMON / BCUT / BR_THR
      COMMON / THIS_ENERGY / ENERGY   !  INITIAL TOTAL ENERGY IN MEV
      COMMON / RATIO / RATEM
      
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR, Z_Wgt, FACTOR_MS 
      
      GI=EXP(GILN)
        IF (GI .LT. 1. .AND. GI .GT. 0.) THEN
         DELTA=GI/(1.-GI)*RATEM*136.*EXP(Z_G)    ! (2.7.50-52)

C  ####################################
c Compute s parameter for LPM corrections
c Stanev et al. PRD 25,5 pp. 1291-1304.
        SOLD=1.
        DO J=1,100
        S=1.37E3 *                ! expr. 7 p.1292 Stanev et al.
     #SQRT( 0.511/ENERGY*GI/(1.-GI)*RADLEN/XI(SOLD,Z_W) )
        IF (ABS((S-SOLD)/S) .LT. 1.E-5) GO TO 1
        SOLD=S
        END DO
        WRITE (*,*) 'LPM DID NOT CONVERGE MAC'
1       PHIFAC=MIN(XI(S,Z_W)*PHIFUN(S),1.)        !  ~ 6 S
        PSIFAC=MIN(PSIFUN(S)*XI(S,Z_W),1.)        !  ~ 4 S

        ELSE
         DELTA=1.E25
         PHIFAC=1.
         PSIFAC=1.
        END IF
      DXDGBR1=((12.+1./(LOG(183.))/(1.+Z_UP50))/9.*
     #AFUN(DELTA)*(1.-GI)*
     # PHIFAC + BFUN(DELTA)*GI**2*PSIFAC ) * 
     #(Z_A+Z_B-Z_F50)/(Z_AB-Z_F)
      IF ( ENERGY .LE. 50.D0) THEN
        X=ENERGY/0.511
c Koch and Motz empirical
c correction for bremsstrahlung (depends on Z)
C Has to be corrected. It is a small correction though.
        DXDGBR1=DXDGBR1*
     #(1.D0+3.130/(EXP(SQRT(2.5/X))+EXP(SQRT(2.5*X))))
      END IF
      RETURN
      END


C **********************************
C Differential cross section in the Moon
c as a function of energy fraction carried
c by electron or positron.
c **********************************
      FUNCTION DXDGPA(GI)
      save
      parameter (Maximum_Z=100)
      COMMON / MEDIUM_COMPOSITION / Max_Z, FRAC_NUM(Maximum_Z), 
     #ZATOM(Maximum_Z), FRAC_WEIGHT(Maximum_Z),
     #FC_Zi(Maximum_Z), Xi_Zi(Maximum_Z),
     #FZ4_50(Maximum_Z),
     # AATOM(Maximum_Z),  XI_ion(Maximum_Z),PHASE
     
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR, Z_Wgt, FACTOR_MS 
     
      COMMON/ZFUNCT/Z_B,Z_F,Z_F50,Z_AB,Z_A,Z_T,Z_V,Z_G,
     #Z_U,Z_P,Z_VP50,Z_UG50
     #,Z_VG50,Z_UP50, Z_W
      COMMON / RATIO / RATEM
      COMMON / THIS_ENERGY / ENERGY   !  ENERGY IN MEV

        IF (GI .NE. 1. .AND. GI .NE. 0.) THEN
c delta (see EGS4 manual) (2.7.50-52).
         DELTA=1./GI/(1.-GI)*RATEM*136.*EXP(Z_G)
c LPM effect. See Stanev et al. Phys. Rev. D 25, 5 (1985) pp.1291-1304
        SOLD=1.
c Iterative calculus of s parameter.
        DO J=1,100
        S=1.37E3 * 
     #SQRT( 0.511/ENERGY/GI/(1.-GI)*RADLEN/XI(SOLD,Z_W) ) !Stanev et al. Formula 8
                                                       !p.1293 Physical Review D
        IF (ABS((S-SOLD)/S) .LT. 1.E-5) GO TO 1
        SOLD=S
        END DO
        WRITE (*,*) 'LPM DID NOT CONVERGE MAC'
1       PHIFAC=MIN(XI(S,Z_W)*PHIFUN(S),1.)
        PSIFAC=MIN((1.5*PSIFUN(S)-PHIFUN(S)/2.)*XI(S,Z_W),1.)

        ELSE
         DELTA=1.E25
         PHIFAC=1.
         PSIFAC=2./3.
        END IF
C      ELSE
C      DXDGPA=0.
C      RETURN
C      END IF

      DXDGPA=( (24.-1./(LOG(183.))/(1.+Z_UP50))/36.*
     #CFUN(DELTA)*PHIFAC+                        ! (2.7.98) + LPM corrections
     # PSIFAC*(12.+1./(LOG(183.))/
     #(1.+Z_UP50))/9.*AFUN(DELTA)*(GI-0.5)**2)
     #          *(Z_A+Z_B-Z_F50)/(Z_AB-Z_F)
      RETURN
      END

C ##################################
c Page 41. EGS4 manual. (2.7.60)
c **********************************
      FUNCTION AFUN(DELTA)   !  SCREENING FUNCTIONS
      save
      COMMON/ZFUNCT/Z_B,Z_F,Z_F50,Z_AB,Z_A,Z_T,Z_V,Z_G,
     #Z_U,Z_P,Z_VP50,Z_UG50
     #,Z_VG50,Z_UP50, Z_W
      AFUN=(3.*PHI1(DELTA)-PHI2(DELTA)+8.*Z_VG50)/
     #     (0.66666667+8.*(LOG(183.)+Z_VG50))
      RETURN
      END
C ##################################
c Page 42. EGS4 manual. (2.7.61)
c **********************************
      FUNCTION BFUN(DELTA)   !  SCREENING FUNCTIONS
      save
      COMMON/ZFUNCT/Z_B,Z_F,Z_F50,Z_AB,Z_A,Z_T,Z_V,Z_G,
     #Z_U,Z_P,Z_VP50,Z_UG50
     #,Z_VG50,Z_UP50, Z_W
      BFUN=0.25*(PHI1(DELTA)+Z_VG50) / (LOG(183.)+Z_VG50)
      RETURN
      END
C ##################################
c Page 47. EGS4 manual. (2.7.100)
c **********************************
      FUNCTION CFUN(DELTA)   !  SCREENING FUNCTIONS
      save
      COMMON/ZFUNCT/Z_B,Z_F,Z_F50,Z_AB,Z_A,Z_T,Z_V,Z_G,
     #Z_U,Z_P,Z_VP50,Z_UG50
     #,Z_VG50,Z_UP50, Z_W
      CFUN=(3.*PHI1(DELTA)+PHI2(DELTA)+16.*Z_VG50)/
     #     (-0.66666667+16.*(LOG(183.)+Z_VG50))
      RETURN
      END
C ##################################
      FUNCTION CGAMMA(GAMMA)   !  SCREENING FUNCTIONS
      save
      DIMENSION GMI(9),GAMFUN(9),DELGI(8)
      DATA GMI/ 2., 2.5, 3., 4., 5., 6., 8., 10., 15./
      DATA DELGI/ 0.1, 0.06, 0.04, 0.025, 0.015, 0.01, 0.005, 0.002/
      DATA GAMFUN/ 0.21, 0.16, 0.13, 0.09, 0.065, 0.05,0.03,0.02,0.01/
      DO I=2,8
      DELGAM=GAMMA-GMI(I)
      IF (DELGAM .LE. 0.)  GO TO 10
      END DO
      CGAMMA=GAMFUN(I-1)-DELGAM*DELGI(I-1)
      CGAMMA=0.8844787*CGAMMA          !  correction for continuity
      RETURN
10    CGAMMA=GAMFUN(I)-DELGAM*DELGI(I-1)
      CGAMMA=0.8844787*CGAMMA          !  correction for continuity
      RETURN
      END
C ##################################
c Page 30. EGS4 manual. (2.7.14)
c **********************************
      FUNCTION PHI1(DELTA)   !  SCREENING FUNCTIONS
      save
        IF (DELTA .LE. 1.) THEN
          PHI1=20.867+DELTA*(-3.242+DELTA*0.625)
        ELSE
          PHI1=21.12-4.184*LOG(DELTA+0.952)
        END IF
      RETURN
      END
c ##################################
c Page 30. EGS4 manual. (2.7.15)
c **********************************
      FUNCTION PHI2(DELTA)   !  SCREENING FUNCTIONS
      save
        IF (DELTA .LE. 1.) THEN
          PHI2=20.029-DELTA*(1.930+DELTA*0.086)
        ELSE
          PHI2=21.12-4.184*LOG(DELTA+0.952)
        END IF
      RETURN
      END
C ##################################
c Stanev et al expression 16
c **********************************
      FUNCTION XI(S,Z)   !  SCREENING FUNCTIONS
      save
      S1=(Z**0.333333333/191.)**2
      IF (S .GT. 1.) THEN
        XI=1.
      ELSE IF (S .GE. S1) THEN
        XI=1.+LOG(S)/LOG(S1)
      ELSE
        XI=2.
      END IF
      RETURN
      END
C ##################################
c Stanev et al expression 15
c **********************************
      FUNCTION PSIFUN(S)   ! [G+2PHI]/3, LPM-SCREENING FUNCTIONS
      save
      IF (S .LE. 0.01) THEN
        PSIFUN=4.*S
      ELSE IF (S .LT. 2.) THEN
        PSIFUN=1.-EXP(-4.*S*
     #         (1+2.*S/(1.+S*(3.936+S*(4.97-S*(0.05-7.50*S))))) )
      ELSE
        PSIFUN=1.
      END IF
      RETURN
      END
C ##################################
c Stanev et al expression 14
c **********************************
      FUNCTION PHIFUN(S)   ! PHI, LPM-SCREENING FUNCTIONS
      save
      PARAMETER (PI = 3.141592654)
      IF (S .LE. 0.01) THEN
        PHIFUN=6.*S*(1.-PI*S)
      ELSE IF (S .LT. 2.) THEN
        PHIFUN=1.-EXP(-6.*S*(1.+(3.-PI)*S)+
     #         S**3./(0.623+S*(0.796+S*0.658)) )
      ELSE
        PHIFUN=1.
      END IF
      RETURN
      END
                                



