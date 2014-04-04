C **************************************************
C **************************************************
C  ================================================*
C            PROGRAM  EMCASNEW                     *
C  ================================================*
C    ALL ELECTROMAGNETIC PROCESSES ARE INCLUDED    *
C **************************************************
c 
c ################################################# 
c ################################################# 
c WARNING:
c This version optimized and developed for beam tests 
c in a small volume (< 10 m^3) filled with air.
c The number of substeps KL has been increased by 
c a factor NSTEPS_FAC, otherwise the code does not 
c give correct results. When used in a denser medium
c and/or larger volume the code might be very slow.
c J.Alvarez-Muniz 22-Jun-2010
c 
c Uses Approximation C - very slow !
c J.Alvarez-Muniz 2-Mar-2011
c 
c Wire antenna orientation and pattern implemented
c J.Alvarez-Muniz 2-Mar-2011
c ################################################# 
c ################################################# 
c
c Computes electric field in the frequency-domain
c using 2 methods:
c ZHS - delays are calculated by hand
c Brute-force - relative phases are calculated from
c space-time positions of particles      
c
c and also in the time domain
c see notes Andres Romero and myself

c 2-level thinning implemented

c Runs the same shower with different thinning levels 
c simultaneously and also performs a hybrid simulation 
c If no thinning levels are chosen (NTHIN=0) performs
c hybrid simulation only 

c Full 3D version i.e. observer can be placed at any
c phi angle not limited to phi=0  (3 May 2010)

c Fresnel calculation possible given grid of observers
c in (x,y,z). Fraunhofer approximation still used in 
c formula for the electric field.  
C ---------------------------------------------------------------------

C Last updated: J. Alvarez-Muniz 03-May-2010

      PARAMETER (PI=3.141592653589793238)
      PARAMETER (maxnu=2000)                   ! Max. number of frequencies
      PARAMETER (MAXJ=50)                     ! Max. number of theta angles
      PARAMETER (MAXK=5)                      ! Max. number of phi angles
      PARAMETER (MAXNA=200)                     ! Max. number of antenna positions (Fresnel) 
      PARAMETER (MAXT=10000)                  ! Max. number of time bins 
      PARAMETER (MAXTHALF=MAXT/2)             
      PARAMETER (MAXNR=10000)                 ! Max. number of rad. lengths
      PARAMETER (MAXNESTACK=100000)           ! Max. number of e- + e+ in stack 
      PARAMETER (MAXNGSTACK=1000000)          ! Max. number of gammas in stack 
      PARAMETER (MAXNEM=100000000)            ! Max. number of EM particles from beam      
      PARAMETER (NTHIN_MAX=2)                ! Number of simultaneous thinning levels 

      PARAMETER (APPROX_A=1, APPROX_B1=0, APPROX_B2=-1, APPROX_C=2)
      dimension XNEXCES(MAXNR,NTHIN_MAX)
      dimension XNELEC(MAXNR,NTHIN_MAX)
      integer*8 NEXCESNOTHIN(MAXNR,NTHIN_MAX)
      integer*8 NELECNOTHIN(MAXNR,NTHIN_MAX)
      dimension XNEPSUM(MAXNR,NTHIN_MAX)
      integer*8 NEPSUMNOTHIN(MAXNR,NTHIN_MAX)
      dimension XNPOSI(MAXNR,NTHIN_MAX)
      dimension RMIDDL(20)
      integer*8 NPOSINOTHIN(MAXNR,NTHIN_MAX)

      CHARACTER*80 flnamefreq

c Frequency domain - Fraunhofer      
      COMPLEX*16 ZSUM(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM_X(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM_Y(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM_Z(MAXJ,MAXK,MAXNU,NTHIN_MAX)

      COMPLEX*16 ZSUM2(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_X(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_Y(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_Z(MAXJ,MAXK,MAXNU,NTHIN_MAX)

c Frequency domain - Fraunhofer      
      COMPLEX*16 ZSUM2_FR(MAXNA,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_X_FR(MAXNA,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_Y_FR(MAXNA,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_Z_FR(MAXNA,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_ANT_FR(MAXNA,MAXNU,NTHIN_MAX)

c Time domain - Fraunhofer      
      DOUBLE PRECISION ASUM(MAXT,MAXJ,MAXK,NTHIN_MAX)
      DOUBLE PRECISION ASUM_X(MAXT,MAXJ,MAXK,NTHIN_MAX)
      DOUBLE PRECISION ASUM_Y(MAXT,MAXJ,MAXK,NTHIN_MAX)
      DOUBLE PRECISION ASUM_Z(MAXT,MAXJ,MAXK,NTHIN_MAX)
      DOUBLE PRECISION ASUM_PERP(MAXT,MAXJ,MAXK,NTHIN_MAX)
      DOUBLE PRECISION EF(MAXK),EF_PERP(MAXK) 
      DOUBLE PRECISION EF_X(MAXK),EF_Y(MAXK),EF_Z(MAXK)

c Time domain - Fresnel      
      DOUBLE PRECISION ASUM_FR(MAXT,MAXNA,NTHIN_MAX)
      DOUBLE PRECISION ASUM_X_FR(MAXT,MAXNA,NTHIN_MAX)
      DOUBLE PRECISION ASUM_Y_FR(MAXT,MAXNA,NTHIN_MAX)
      DOUBLE PRECISION ASUM_Z_FR(MAXT,MAXNA,NTHIN_MAX)
      DOUBLE PRECISION ASUM_ANT_FR(MAXT,MAXNA,NTHIN_MAX)
      DOUBLE PRECISION EF_X_FR(MAXNA),EF_Y_FR(MAXNA),EF_Z_FR(MAXNA)
      DOUBLE PRECISION EF_FR(MAXNA)
      DOUBLE PRECISION EF_ANT_FR(MAXNA)

      DOUBLE PRECISION FACTOR_T,FAC_T,DT_BIN_NS,DT_BIN_S

      DOUBLE PRECISION FACFRQ,FRQCOM,THETA,SINTET,SINTET2
      DOUBLE PRECISION COSTET,COSTET2,SINMU,CSMU_1,COSMU
      DOUBLE PRECISION PHI,COSPHI,SINPHI
      DOUBLE PRECISION THETMX,THETDG,RADEGR
      DOUBLE PRECISION ANGLES(MAXNU)

c Frequency domain - Fraunhofer      
      DOUBLE PRECISION RADIUS(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  PHASE(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  FREQ(MAXNU),CHERAN,CHERDG,CHECOS,
     #  WIDTH(MAXNU),FACTOR
      DOUBLE PRECISION RADIUS_X(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  PHASE_X(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION RADIUS_Y(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  PHASE_Y(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION RADIUS_Z(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  PHASE_Z(MAXJ,MAXK,MAXNU,NTHIN_MAX)

      DOUBLE PRECISION ARADIUS(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  APHASE(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION ARADIUS_X(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  APHASE_X(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION ARADIUS_Y(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  APHASE_Y(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION ARADIUS_Z(MAXJ,MAXK,MAXNU,NTHIN_MAX),
     #  APHASE_Z(MAXJ,MAXK,MAXNU,NTHIN_MAX)

c Frequency domain - Fresnel      
      DOUBLE PRECISION REFIDX 
      DOUBLE PRECISION XANT,YANT,ZANT
      DOUBLE PRECISION T0_ANT(MAXNA)
      DOUBLE PRECISION ETA_XANT,ETA_YANT,ETA_ZANT
      DOUBLE PRECISION ETA_PERP_XANT,ETA_PERP_YANT,ETA_PERP_ZANT
      DOUBLE PRECISION ARADIUS_FR(MAXNA,MAXNU,NTHIN_MAX),
     #  APHASE_FR(MAXNA,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION ARADIUS_X_FR(MAXNA,MAXNU,NTHIN_MAX),
     #  APHASE_X_FR(MAXNA,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION ARADIUS_Y_FR(MAXNA,MAXNU,NTHIN_MAX),
     #  APHASE_Y_FR(MAXNA,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION ARADIUS_Z_FR(MAXNA,MAXNU,NTHIN_MAX),
     #  APHASE_Z_FR(MAXNA,MAXNU,NTHIN_MAX)
      DOUBLE PRECISION ARADIUS_ANT_FR(MAXNA,MAXNU,NTHIN_MAX),
     #  APHASE_ANT_FR(MAXNA,MAXNU,NTHIN_MAX)


      DOUBLE PRECISION GCM_NS

      double precision zobs
      

C Random number generator common
      COMMON/S_LUDATR/MRLU(6),RRLU(100)
      COMMON/TTRE/T97(20)                
      COMMON/TTR/TRR(20)
      COMMON/NLE/XLE(20,MAXNR,NTHIN_MAX),
     # XLP(20,MAXNR,NTHIN_MAX)/NLG/LG(20)
      COMMON/NLENOTHIN/XLENOTHIN(20,MAXNR,NTHIN_MAX),
     #XLPNOTHIN(20,MAXNR,NTHIN_MAX)
      COMMON/ZOB/ZOBS

      integer*8 NRE,NRPOS,NDELTE
      COMMON/NVE/NRE(0:MAXNR,NTHIN_MAX),NRPOS(0:MAXNR,NTHIN_MAX),
     #           NDELTE

      integer*8 NRENOTHIN,NRPOSNOTHIN
      COMMON/NVENOTHIN/NRENOTHIN(0:MAXNR,NTHIN_MAX),
     #NRPOSNOTHIN(0:MAXNR,NTHIN_MAX)

      COMMON/NVG/NRG
      COMMON / THIN / THIN_MAX(NTHIN_MAX), THIN_MIN(NTHIN_MAX)
      COMMON/WR0/RNUMB

      COMMON / MOMENTS / 
     # WXE(MAXNR,NTHIN_MAX),WYE(MAXNR,NTHIN_MAX),
     # WRE(MAXNR,NTHIN_MAX),
     # WX2E(MAXNR,NTHIN_MAX),WY2E(MAXNR,NTHIN_MAX),
     # WR2E(MAXNR,NTHIN_MAX),
     # WXP(MAXNR,NTHIN_MAX),WYP(MAXNR,NTHIN_MAX),
     # WRP(MAXNR,NTHIN_MAX),
     # WX2P(MAXNR,NTHIN_MAX),WY2P(MAXNR,NTHIN_MAX),
     # WR2P(MAXNR,NTHIN_MAX),
     # WZE(MAXNR,NTHIN_MAX),WZ2E(MAXNR,NTHIN_MAX),
     # WSUME(MAXNR,NTHIN_MAX), 
     # WZP(MAXNR,NTHIN_MAX),WZ2P(MAXNR,NTHIN_MAX),
     # WSUMP(MAXNR,NTHIN_MAX) 

      integer*8 NCALLS_EMPSUM_W
      integer*8 NCALLS_EMPSUM_W_2
      integer*8 NCALLS_EMPSUM_T
      COMMON / THIN_QUALITY_W / NCALLS_EMPSUM_W(NTHIN_MAX)
      COMMON / THIN_QUALITY_W_2 / NCALLS_EMPSUM_W_2(NTHIN_MAX)
      COMMON / THIN_QUALITY_T / NCALLS_EMPSUM_T(NTHIN_MAX)

      COMMON/BLMASS/e_MASS,e2_MAS
      
      COMMON /GRID/ FRQCOM(MAXNU), THETA(MAXJ), PHI(MAXK), NUMAX,
     # JEND(MAXNU),JMAX,KMAX,JCHER,J90(MAXNU),J180(MAXNU)
      COMMON / SPEEDUP / SINTET(MAXJ),COSTET(MAXJ),SINTET2(MAXJ),
     # COSPHI(MAXK),SINPHI(MAXK),
     # COSTET2(MAXJ),SINMU(MAXJ),CSMU_1(MAXJ),COSMU(MAXJ),FACFRQ(MAXNU)

      COMMON /ANT1/ NAMAX
      COMMON /ANT2/ XANT(MAXNA),YANT(MAXNA),ZANT(MAXNA)
      COMMON /ANT3/ ETA_XANT(MAXNA),ETA_YANT(MAXNA),ETA_ZANT(MAXNA)
      COMMON /ANT4/ ETA_PERP_XANT(MAXNA),ETA_PERP_YANT(MAXNA),
     #              ETA_PERP_ZANT(MAXNA)

      COMMON /CONSTANTS_T/ FACTOR_T
      COMMON /BIN_T/ DT_BIN_NS,DT_BIN_S 
      COMMON /TIMES/ TIME_MIN,TIME_MAX 

      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS

      COMMON /REFR/ REFIDX

      COMMON /UNITS_T/ GCM_NS

      COMMON/ECRIT/E_CRIT

      COMMON / PHOTON / TRP 

      COMMON/ BREMSS_XSECTION / BR_E(400), BR_X(400)
      COMMON / PAIR_XSECTION / PP_E(400), PP_X(400)

      COMMON / TRACKS_W / IWRITE_W
      COMMON / TRACKS_W_2 / IWRITE_W_2
      COMMON / TRACKS_T / IWRITE_T
      COMMON / SUBTRACKS / NZEROS,NIN,NCROSS1,NCROSS 

      COMMON / HYBRID / E_HYBRID
      COMMON / HYBRID_SIZE / XNEPSUM_HYBRID(MAXNR)
      COMMON / HYBRID_TRACK / TOT_TCK_HYBRID
      COMMON / CONSERVATION / E_CONS

      COMMON / MAX_DEPTH / NRMAX 

      COMMON / THIN_LEVELS / NTHIN

      COMMON / FREQ_ARRAY / FREQ

      INTEGER APPROX,KEEP2(NTHIN_MAX),KEEP3(NTHIN_MAX),
     #KEEP4(NTHIN_MAX),KEEPE(NTHIN_MAX)

C     These hold (energy,time, 
      DIMENSION DEN(10,MAXNESTACK),DGN(10,MAXNGSTACK)
      DIMENSION WDEN(MAXNESTACK,NTHIN_MAX),
     #          WDGN(MAXNGSTACK,NTHIN_MAX)
      DIMENSION DEN8(3,MAXNESTACK),DGN8(3,MAXNGSTACK)
      DIMENSION ETHRKE(11)
      DOUBLE PRECISION TCKSUMPRO(10,NTHIN_MAX),TCKABS(10,NTHIN_MAX)
      DOUBLE PRECISION TCKPRO(10,NTHIN_MAX),DdCT
      REAL CT_1,CT_1_N,CT_2
      DOUBLE PRECISION DELTA_CT,DELTA_CT_STEP
      DOUBLE PRECISION TCKSUMNOPRO(10,NTHIN_MAX)

      DOUBLE PRECISION Wgt(NTHIN_MAX)
      DOUBLE PRECISION Wgt2(NTHIN_MAX),Wgt3(NTHIN_MAX),
     #                 WgtE(NTHIN_MAX)
      DOUBLE PRECISION Wgt4(NTHIN_MAX),THIN_METH
      DOUBLE PRECISION E2_AV, E3_AV, E4_AV, EE_AV
      CHARACTER*4 SUFFIX
      CHARACTER*80 bunch_file 
      
!      CHARACTER*12 SENTEN(3)
!      CHARACTER*20 Input_File
C
C *************************************************************************
C ********************************* Data *********************************
C *************************************************************************
c --------------------------------------------------------------------------
      DATA e_MASS/0.51099906/,e2_MAS/0.26112004/   !  electron mass parameters
C --------------------------------------------------------------------------
      DATA TI/0./                 ! Initial thickness 
      DATA DLI/0.1E-8/,DMI/0.1E-8/,DNI/0.99999999/   ! Initial direc. cos
ccc      DATA WgtI/1.000000000/                         ! Weight of primary particle
      DATA RNUMB/5./  
      DATA RFRAC/5./    ! fraction of range allowed for ionization loss steps
      DATA XLOG/5./     ! Number of intervals per decade in TCKSUMPRO
      DATA EPSILN/1.E-7/ ! CH_THR-Epsilon = Energy of e in last step
c

C 
C ********************* multiple scattering *************************  
      DATA TMAX/0.1/       ! Maximum thickness for multiscatt. in T0 units
      DATA RG1/0.55/,RG2/0.9841/,RG3/3./   ! Parameters for range approx.
      DATA E_MS/ 21.0 /    ! Multiple scattering energy (21 MeV)
C --------------------------------------------------------------------------
     
      call Medium_Properties       ! Fills in COMMONs: MEDIUM, BREMSS_XSECTION, PAIR_XSECTION
c Critical energy of medium in MeV
      E_CRIT=T0*E_MS/R0/100.

      NRG=0
      RG=RG1/RFRAC/T0              !  Normalization of 20% range function
      RFACT=10.**(1./RNUMB)        !  Multiplicative factor for R0
      RINIT=2.*RFACT/(1.+RFACT)

      Edep_int = 0 ! Integrated energy deposition

      DO 60 I=1,20
c Units of Moliere radius/100
c i.e. multiply by Moliere radius/100 to convert to g cm^-2      
       RMIDDL(I)=RINIT*RFACT**(I-1)
       T97(I)=0.
       TRR(I)=0.

   60 LG(I)=0

      THIN_METH=1.0000000d0
      MAX_NEN=0;
      MAX_NGN=0


      PRINT *,' Initialization of random number gen. Give seed : '
      READ(*,*) MRLU(1)
      WRITE(*,*) MRLU(1)

      PRINT *, 'Run under which approximation (A=1, B1=0, B2=-1, C=2)?'
      READ(*,*) APPROX
      WRITE(*,*) APPROX

      PRINT *,'Record Particle Tracks (y=1,n=0)?'
      READ(*,*) IWRITE_W_2,IWRITE_W,IWRITE_T
      WRITE(*,*) IWRITE_W_2,IWRITE_W,IWRITE_T

      PRINT *,'Compute E-field in freq. and time domain'
      PRINT *,'using explicit tracks & a la ZHS & time domain (y=1,n=0)'
      PRINT *,'Fraunhofer'
      READ(*,*) IFIELD_W_2,IFIELD_W,IFIELD_T
      WRITE(*,*) IFIELD_W_2,IFIELD_W,IFIELD_T

      PRINT *,'Compute E-field in freq. and time domain'
      PRINT *,'using explicit tracks & a la ZHS & time domain (y=1,n=0)'
      PRINT *,'Fresnel'
      READ(*,*) IFIELD_W_2_FRESNEL,IFIELD_T_FRESNEL
      WRITE(*,*) IFIELD_W_2_FRESNEL,IFIELD_T_FRESNEL

      PRINT *,'Binning in time [ns] when computing field in time domain'
      READ(*,*) DT_BIN_NS
      WRITE(*,*) DT_BIN_NS

      TIME_MIN=-300.  ! Min. time [ns] for output 
      TIME_MAX=300.   ! Max. time [ns] for output 

      DT_BIN_S=DT_BIN_NS*1.D-9

      WRITE(*,41)
   41 FORMAT(2X,'ENTER E0 (IN MEV)')
      READ(*,*) EI
      WRITE(*,*) EI

      PRINT *,'Number of thinning levels'
      READ(*,*) NTHIN
      WRITE(*,*) NTHIN 

      IF (NTHIN.EQ.0) GOTO 444

c Loop over thinning levels
      do ith=1,NTHIN
        PRINT*, 'ENTER thin max and min threshold THIN_MAX and THIN_MIN'
        READ(*,*) THIN_MAX(ith),THIN_MIN(ith)
        write(*,*) THIN_MAX(ith),THIN_MIN(ith)
      end do

444   continue

      do ith=NTHIN+2,NTHIN_MAX
         THIN_MAX(ith)=0.
         THIN_MIN(ith)=0.
      end do

      PRINT*,'ENTER particle type'
      READ(*,*) TYPI
      WRITE(*,*) TYPI

      PRINT*,'ENTER ZOBS'
      READ(*,*) ZOBS
      write(*,*) ZOBS

      WRITE(*,42)
   42 FORMAT(2X,'ENTER THE NUMB.OF SHOWERS')
      READ(*,*) NOSHOW
      write(*,*) NOSHOW

       WRITE(*,*) 'ENTER PROCESSES TO IGNORE: COMPTON, ANNIHILATION,'
     #            ,' MOELLER, BHABHA'
       READ(*,*) ICOMPT
       READ(*,*) IANNHI
       READ(*,*) IBHABH
       READ(*,*) IMOELL
       write(*,*) ICOMPT,IANNHI,IBHABH,IMOELL

       WRITE(*,44)
   44 FORMAT(2X,'ENTER SUFFIX #### CREATES 11 FILES',/
     #   5X,'DR####.EMP',5X,'Run data and Radial Distribution',/
     #   5X,'DT####.EMP',5X,'Depth Distribution',/
     #   5X,'DTNOTHIN####.EMP',5X,'Depth Distribution kept particles',/
     #   5X,'TK####.EMP',5X,'Tracklength Information',/
     #   5X,'PS####.EMP',5X,'Pulse Angular Distribution Mod+Phas',/
     #   5X,'PF####.EMP',5X,'Pulse Frequency Spectra Mod+Phas',/
     #   5X,'ZS####.EMP',5X,'Pulse Angular Distribution Real+Im',/
     #   5X,'ZF####.EMP',5X,'Pulse Frequency Spectra Real+Im',/
     #   5X,'MM####.EMP',5X,'Moments of inertia of shower vs depth',/
     #   5X,'TQ####.EMP',5X,'Thinning quality factor (whole shower)',/
     #   5X,'TQM####.EMP',5X,'Thinning quality factor (at shower max)')

      READ(*,45) SUFFIX
      write(*,45) SUFFIX
   45 FORMAT(A4)

! Simple way of computing the time of the run: open a file
! at the beginning of the run and close it immediately and 
! compare the time at which this file was created with the
! time at which the output files (created at the end of the
! run) were created.
       open(unit=77,status='unknown',file='time'//SUFFIX//'.EMP')
       close(77)

       open(unit=88,status='unknown',
     #file='first_interaction'//SUFFIX//'.EMP')

      IF (IWRITE_W .EQ. 1) THEN
        do ith=1,NTHIN  ! One file with track info per thinning level 
          OPEN(UNIT=39+ith, FILE='TRACKS_W_'//SUFFIX//'.DAT',
     #         STATUS='UNKNOWN')
        end do
      END IF

      IF (IWRITE_W_2 .EQ. 1) THEN
        do ith=1,NTHIN  ! One file with track info per thinning level 
          OPEN(UNIT=59+ith, FILE='TRACKS_W_2_'//SUFFIX//'.DAT',
     #         STATUS='UNKNOWN')
        end do
      END IF

      IF (IWRITE_T .EQ. 1) THEN
        do ith=1,NTHIN  ! One file with track info per thinning level 
          OPEN(UNIT=69+ith, FILE='TRACKS_T_'//SUFFIX//'.DAT',
     #         STATUS='UNKNOWN')
        end do
      END IF


c Read grid in frequencies (Fraunhofer and Fresnel calculations) 
c and theta and phi angles (Fraunhofer only)
c and antenna positions (Fresnel only)

c Reads grid of frequencies from a file grid_frequency.inp
      READ(*,*) flnamefreq
      open(unit=77,status='old',file=flnamefreq)
      WRITE(*,*) ' ENTER # OF FREQ POINTS, (FREQ,WIDTH-degrees)'
      READ(77,*) NUMAX
      DO III=1,NUMAX
        READ(77,*) FREQ(III)
      END DO
      write(*,*) NUMAX,(FREQ(NU),NU=1,NUMAX)

      WRITE(*,*) ' ENTER # THETA ANGLES (offset w.r.t. Cherenkov)'
      READ(*,*) JMAX,(ANGLES(J),J=1,JMAX)
      write(*,*) JMAX,(ANGLES(J),J=1,JMAX)

      WRITE(*,*) ' ENTER # PHI ANGLES'
      READ(*,*) KMAX,(PHI(K),K=1,KMAX)
      write(*,*) KMAX,(PHI(K),K=1,KMAX)

      WRITE(*,*) ' ENTER # OF ANTENNA POSITIONS '
      READ(*,*) NAMAX
      READ(*,*) (XANT(N),YANT(N),ZANT(N),N=1,NAMAX)
      WRITE(*,*) NAMAX
      WRITE(*,*) (XANT(N),YANT(N),ZANT(N),N=1,NAMAX)

      WRITE(*,*) ' ENTER ANTENNA ORIENTATIONS '
      READ(*,*) (ETA_XANT(N),ETA_YANT(N),ETA_ZANT(N),N=1,NAMAX)
      WRITE(*,*) (ETA_XANT(N),ETA_YANT(N),ETA_ZANT(N),N=1,NAMAX)
      WRITE(*,*) ' ENTER PERPENDICULAR ORIENTATIONS TO ANTENNA'
      READ(*,*) 
     #(ETA_PERP_XANT(N),ETA_PERP_YANT(N),ETA_PERP_ZANT(N),N=1,NAMAX)

c Normalize to modulus=1 if needed
      DO IA=1,NAMAX
        ETA_MOD=SQRT(ETA_XANT(IA)*ETA_XANT(IA)+
     #               ETA_YANT(IA)*ETA_YANT(IA)+
     #               ETA_ZANT(IA)*ETA_ZANT(IA))

        ETA_XANT(IA)=ETA_XANT(IA)/ETA_MOD
        ETA_YANT(IA)=ETA_YANT(IA)/ETA_MOD
        ETA_ZANT(IA)=ETA_ZANT(IA)/ETA_MOD
        write(*,*) IA,ETA_XANT(IA),ETA_YANT(IA),ETA_ZANT(IA),ETA_MOD

        ETA_PERP_MOD=SQRT(ETA_PERP_XANT(IA)*ETA_PERP_XANT(IA)+
     #                    ETA_PERP_YANT(IA)*ETA_PERP_YANT(IA)+
     #                    ETA_PERP_ZANT(IA)*ETA_PERP_ZANT(IA))

        ETA_PERP_XANT(IA)=ETA_PERP_XANT(IA)/ETA_PERP_MOD
        ETA_PERP_YANT(IA)=ETA_PERP_YANT(IA)/ETA_PERP_MOD
        ETA_PERP_ZANT(IA)=ETA_PERP_ZANT(IA)/ETA_PERP_MOD
        write(*,*) IA,ETA_PERP_XANT(IA),ETA_PERP_YANT(IA),
     #             ETA_PERP_ZANT(IA),ETA_PERP_MOD
      END DO

C --------------------------------------------------------------------
c Compute hybrid energy: energy below which shower is hybrid
c      E_AUX=MIN(E_LPM,E_LPM_BREMSS)
      E_HYBRID=MIN(E_LPM/10.,EI/100.)
      write(*,*) 'E_LPM [TeV] ',E_LPM/1.e6
      write(*,*) 'E_Hybrid [MeV] ',E_HYBRID

      THIN_MIN(NTHIN+1)=E_HYBRID
      THIN_MAX(NTHIN+1)=E_HYBRID

C --------------------------------------------------------------------
      REFIDX=REF_N
      RADEGR=180.D0/PI          !  Conversion factor degrees <--> radians
c factor=2 pi/(rho c) used to convert electric field to S.I. units.
      FACTOR=2.D0*PI/RHO/2.99792458D4 !  Frequency must be in MHz
      FACTOR_T=1.D0/RHO/2.99792458D10 
      FACTOR_F=100.*RHO          ! Fresnel calculation
                                 ! Distance is in g cm-2 and field
                                 ! in freq. domain should be in V/MHz/m
                                 ! and in time domain should be in V/m 
      CHECOS=1.D0/REF_N
      CHERAN=ACOS(CHECOS)        ! Radians
      CHERDG=CHERAN*RADEGR       ! Degrees

      GCM_NS=2.99792458D1*RHO    ! From g/cm^2 to ns
C
C **************************************************************************
C GRID in frequency and observation angle with respect to shower axis.
C It also calculates factors that are frequently used and puts them in 
C an array in order to speed up the code.
C **************************************************************************
C
      DO NU=1,NUMAX
        FACFRQ(NU)=FACTOR*FREQ(NU)
        FRQCOM(NU)=FREQ(NU)
        THETDG=CHERDG+WIDTH(NU)  ! Degrees
        THETMX=THETDG/RADEGR     ! Radians
        JSTART=1
        JEND(NU)=JMAX
      END DO
      
      DO J=1,JMAX
        IF (abs(ANGLES(J)) .LT. 1e-4) THEN JCHEC=J
        IF (ANGLES(J) .LT. -CHERDG) THEN
          THETA(J) = 0.
        ELSE
          THETA(J)=CHERDG+ANGLES(J)
        END IF
      END DO
      
      DO J=1,JMAX
        THETA(J)=THETA(J)/RADEGR
        SINTET(J)=SIN(THETA(J))
        SINTET2(J)=SIN(THETA(J))*SIN(THETA(J))
        SINMU(J)=REF_N*SINTET(J)
        COSTET(J)=COS(THETA(J))
        COSTET2(J)=COS(THETA(J))*COS(THETA(J))
        CSMU_1(J)=REF_N*COSTET(J)-1.D0
        COSMU(J)=REF_N*COSTET(J)
      END DO

      DO K=1,KMAX
        PHI(K)=PHI(K)/RADEGR
        SINPHI(K)=SIN(PHI(K))
        COSPHI(K)=COS(PHI(K))
      END DO      

      NO0=0
      NRMAX=ZOBS/T0
      BR_THR=CH_THR-e_MASS           ! Kinetic Threshold energy in MeV
      AP=BR_THR

      DO J=1,10
        ETHRKE(J)=BR_THR*10.**(REAL(-1.+J)/XLOG)
      END DO
      ETHRKE(11)=1.E30           ! Optimizes if for portions of track to sums

      TE1=(CH_THR/e_MASS)-1.     ! Threshold energy K.E. in Me units  Tth/Me
      TRP=0.5*(BR_THR+SQRT(BR_THR*(BR_THR+2.*e_MASS))) ! Photon threshold
      CH_THR_EP=CH_THR*(1.-EPSILN)     ! e Energy at end of last step

C ******************************************************************
C Write medium info     
      WRITE(*,717) TYPI,EI,TI,DLI,DMI,DNI,  Z,FCZ,AZ,ALZ,  
     #             ALZ3,Xi_Z,FACNOR,     T0,BEE,    
     #             C,A,RM,X0,X1,AI  ,TMAX,R0  ,AKO,EKE
  717 FORMAT(1X,'DATA',/ 1X,'INITIAL:   TYPE=',F2.0,'    E=',
     #       E8.2,'    DEPTH=',F5.2,'    DIR-COS=',3F6.3  
     #       / 1X,'MEDIUM:   Z=',F6.3,'   fc(Z)=',F6.3,
     #      '   AZ=Alpha*Z=',F6.4,'   ALZ=ln(Z)=',F6.4,
     #       / 1X,'       ALZ3=ln(183/Z^.3)-fc(Z)=',F5.3,
     #      '      Xi_Z=ln(1150/Z^.7)/ALZ3=', F6.3
     #       / 1X,'       FACNOR=110/(Z+Xi_Z)=',F6.3,
     #      '    RAD LENGTH=',F5.2,'    BEE (MULT SCATT)=',F5.2,
     #       / 1X,'      IONIZ PAR C=',F6.3,'  A=',F6.3,'  Rm=',F5.3,
     #      '  X0=',F4.2,' X1=',F3.1,'  AI=',E8.3 /
     #       / 1X,'MONTECARLO  MAX STEP=TMAX=',F5.3,
     #      '   R0=LAT SPREAD UNIT (gcm-2)=',F6.3, /
     #       / 1X,'AKO=',E9.3,'  EKE=',E9.3/)
      WRITE(*,520) CH_THR,TE1,AP,TRP,ZOBS
  520 FORMAT(1X,'MeV THRESH:   ETH(e)=',F9.3,'  K.E. TE1(e)=',
     #       F9.3,' AP(gam)=',F6.4,' TRP=',F8.4 
     #       / 1X, 'OBSERVATION DEPTH=',F7.0,' gr cm^-2')


C ************************************************************************
C ************************************************************************
C ***********************========================*************************
C ***********************      SHOWER STARTS     *************************
C ***********************========================*************************
C ************************************************************************
C ************************************************************************

c  Initialization 
   50 CONTINUE
      ifirst=0  ! Flag for the primary particle
      write(*,*) '#######################################'
      write(*,*) 'Shower ',NO0+1
      E_CONS=0.           ! Check energy conservation in hybrid shower
      do ith=1,NTHIN_MAX
        DO NRR=1,NRMAX
          XNEXCES(NRR,ith)=0
          NEXCESNOTHIN(NRR,ith)=0
          XNELEC(NRR,ith)=0.
          NELECNOTHIN(NRR,ith)=0
          XNEPSUM(NRR,ith)=0
          NEPSUMNOTHIN(NRR,ith)=0
          XNPOSI(NRR,ith)=0
          NPOSINOTHIN(NRR,ith)=0
        END DO
      end do

c Hybrid shower      
      DO NRR=1,NRMAX
         XNEPSUM_HYBRID(NRR)=0.
      END DO

      do ith=1,NTHIN_MAX
        DO J=1,10
          TCKSUMNOPRO(J,ith)=0.D0
          TCKSUMPRO(J,ith)=0.D0
          TCKABS(J,ith)=0.D0
          TCKPRO(J,ith)=0.D0
        END DO
      end do

c Fields in frequency domain - Fraunhofer
      do ith=1,NTHIN_MAX
        DO NU=1,NUMAX
          DO J=1,JMAX
            DO K=1,KMAX
            ZSUM_X(J,K,NU,ith)=0.
            ZSUM_Y(J,K,NU,ith)=0.
            ZSUM_Z(J,K,NU,ith)=0.
            ZSUM(J,K,NU,ith)=0.
            RADIUS(J,K,NU,ith)=0.
            RADIUS_X(J,K,NU,ith)=0.
            RADIUS_Y(J,K,NU,ith)=0.
            RADIUS_Z(J,K,NU,ith)=0.
            ZSUM2_X(J,K,NU,ith)=0.
            ZSUM2_Y(J,K,NU,ith)=0.
            ZSUM2_Z(J,K,NU,ith)=0.
            ZSUM2(J,K,NU,ith)=0.
            ARADIUS(J,K,NU,ith)=0.
            ARADIUS_X(J,K,NU,ith)=0.
            ARADIUS_Y(J,K,NU,ith)=0.
            ARADIUS_Z(J,K,NU,ith)=0.
            END DO
          END DO
        END DO
      end do

c Fields in time domain - Fraunhofer
      do ith=1,NTHIN_MAX
        DO IT=1,MAXT
          DO J=1,JMAX
            DO K=1,KMAX
            ASUM(IT,J,K,ith)=0.
            ASUM_X(IT,J,K,ith)=0.
            ASUM_Y(IT,J,K,ith)=0.
            ASUM_Z(IT,J,K,ith)=0.
            END DO
          END DO
        END DO
      end do

c Fields in frequency domain - Fraunhofer
      do ith=1,NTHIN_MAX
        DO NU=1,NUMAX
          DO NA=1,NAMAX
            ZSUM2_X_FR(NA,NU,ith)=0.
            ZSUM2_Y_FR(NA,NU,ith)=0.
            ZSUM2_Z_FR(NA,NU,ith)=0.
            ZSUM2_FR(NA,NU,ith)=0.
            ZSUM2_ANT_FR(NA,NU,ith)=0.
            ARADIUS_FR(NA,NU,ith)=0.
            ARADIUS_X_FR(NA,NU,ith)=0.
            ARADIUS_Y_FR(NA,NU,ith)=0.
            ARADIUS_Z_FR(NA,NU,ith)=0.
            ARADIUS_ANT_FR(NA,NU,ith)=0.
          END DO
        END DO
      end do

c Fields in time domain - Fraunhofer
      do ith=1,NTHIN_MAX
        DO IT=1,MAXT
          DO NA=1,NAMAX
            ASUM_FR(IT,NA,ith)=0.
            ASUM_X_FR(IT,NA,ith)=0.
            ASUM_Y_FR(IT,NA,ith)=0.
            ASUM_Z_FR(IT,NA,ith)=0.
            ASUM_ANT_FR(IT,NA,ith)=0.
          END DO
        END DO
      end do

c ----------------------------------
      do ith=1,NTHIN_MAX
        NCALLS_EMPSUM_W(ITH)=0
        NCALLS_EMPSUM_W_2(ITH)=0
        NCALLS_EMPSUM_T(ITH)=0
        do nr=1,maxnr
          XNEXCES(NR,ITH)=0
          XNELEC(NR,ITH)=0
          NEXCESNOTHIN(NR,ITH)=0
          NELECNOTHIN(NR,ITH)=0
          XNEPSUM(NR,ITH)=0
          NEPSUMNOTHIN(NR,ITH)=0
          XNPOSI(NR,ITH)=0
      
          WXE(NR,ITH)=0.
          WYE(NR,ITH)=0.
          WRE(NR,ITH)=0.
          WX2E(NR,ITH)=0.
          WY2E(NR,ITH)=0.
          WR2E(NR,ITH)=0.
          WXP(NR,ITH)=0.
          WYP(NR,ITH)=0.
          WRP(NR,ITH)=0.
          WX2P(NR,ITH)=0.
          WY2P(NR,ITH)=0.
          WR2P(NR,ITH)=0.
          WZE(NR,ITH)=0.
          WZ2E(NR,ITH)=0.
          WSUME(NR,ITH)=0. 
          WZP(NR,ITH)=0.
          WZ2P(NR,ITH)=0.
          WSUMP(NR,ITH)=0. 
        end do
      end do

      do k=1,20
        LG(k)=0
        do nr=1,maxnr
          do ith=1,NTHIN_MAX
            XLE(k,NR,ITH)=0.
            XLP(k,NR,ITH)=0.
            XLENOTHIN(k,NR,ITH)=0.
            XLPNOTHIN(k,NR,ITH)=0.
          end do
        end do
      end do

      do nr=0,nrmax
        do ith=1,NTHIN_MAX
          NRE(NR,ITH)=0
          NRPOS(NR,ITH)=0
          NRENOTHIN(NR,ITH)=0
          NRPOSNOTHIN(NR,ITH)=0
        end do
      end do

c #############################################################
c #############################################################
c Open file containing EM particles information, fill stack and 
c pick one particle to start the shower.
c #############################################################
c #############################################################
      PRINT *, 'Which input file?'
      READ(*,*) bunch_file
      WRITE(*,*) bunch_file
      OPEN(38,FILE=bunch_file,STATUS='OLD')

      nem=0
      do iem=1,1000000000
         nem=nem+1
         if (nem.gt.MAXNEM) then 
           write(*,*) ' WARNING: Max number of EM particles reached '
         end if

         READ(38,*,end=111) WgtI,type_EM,ene_EM,xinj,yinj,time_EM  ! MeV and g/cm2 and ns

         depth_EM=0.  ! zinj in g/cm^2

         ctime_EM=time_EM*1.e-9/FACTOR_T   ! ns to g/cm2

        write(*,454) nem,WgtI,type_EM,ene_EM,xinj,yinj,depth_EM,time_EM 
454      FORMAT(1X,I10,1P,1E10.2,0P,1F4.0,1P,5E14.4)

c         write(*,*) 'Injection time converted to g cm^-2 ',ctime_EM
c REMINDER: Particle type -1=positron, 0=photon, 1=electron

c Fill particle into its corresponding stack
c Assume it travels parallel to shower axis 
c and that the time delay with respect to a particle travelling along showers axis at
c the speed of light is 0 i.e. dCT = rho*c (t - z/c) (affects calculation a la ZHS) 
c However assume the starting time of the particle is the correct one in the bunch
c (affects calculation using "explicit tracks").
          if (type_EM.eq.0) then        ! Photon

            if (ene_EM.gt.TRP) then
            NGN=NGN+1
            IF(NGN.GT.MAXNGSTACK) GO TO 9501
              DGN(1,NGN)=ene_EM
              DGN(2,NGN)=depth_EM
              DGN(4,NGN)=0.0
              DGN(6,NGN)=ctime_EM
              DGN8(1,NGN)=0.1E-8   ! All photons injected parallel to z-axis
              DGN8(2,NGN)=0.1E-8
              DGN8(3,NGN)=0.99999999
              DGN(8,NGN)=xinj        
              DGN(9,NGN)=yinj
              do ith=1,NTHIN+1
                   WDGN(NGN,ith)=WgtI
              end do
            end if

          else if (type_EM.eq.1.or.type_EM.eq.-1) then   ! Electron or positron 

            if (ene_EM.gt.CH_THR) then
            ene_EM=ene_EM+e_MASS
            NEN=NEN+1
            IF(NEN.GT.MAXNESTACK) GO TO 9501
              DEN(1,NEN)=ene_EM
              DEN(2,NEN)=depth_EM
              DEN(4,NEN)=0.0
              DEN(5,NEN)=type_EM
              DEN(6,NEN)=ctime_EM
              DEN8(1,NEN)=0.1E-8        ! All electrons and positrons injected parallel to z-axis
              DEN8(2,NEN)=0.1E-8        
              DEN8(3,NEN)=0.99999999
              DEN(8,NEN)=xinj             
              DEN(9,NEN)=yinj
              do ith=1,NTHIN+1
                 WDEN(NEN,ith)=WgtI 
              end do
            end if

          else                          ! this type of particle shouldn't be here

              write(*,*) 'Wrong input EM particle, ID =  ',type_EM

          end if   ! End loop in filling stacks
        
         EM_eqv = EM_eqv + ene_EM      ! total EM energy injected
  
      end do

111   continue


         write(*,*) '--------------------------------------------------'
         write(*,*) 'Total beam Energy [MeV] = ',EM_eqv
c Number of particles initially in stack
         write(*,*) 'Number of particles initially in stack '
         write(*,*) 'Gammas = ',NGN 
         write(*,*) 'e- & e+ = ',NEN 
         write(*,*) '--------------------------------------------------'

c Initiate shower by selecting the 
c last particle stored in the photon stack

c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c Start with a photon if photon stack is not empty
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       if (NGN.gt.0) then
          TYPI=0
          EI=DGN(1,NGN)
          TI=DGN(2,NGN)
          dCTSTAI=DGN(4,NEN)
          CT_1I=DGN(6,NEN)
          DLI=DGN8(1,NGN)
          DMI=DGN8(2,NGN)
          DNI=DGN8(3,NGN)
          X=DGN(8,NGN)
          Y=DGN(9,NGN)
          do ith=1,NTHIN+1
             Wgt(ith)=WDGN(NGN,ith)
          end do
          NGN=NGN-1
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c Start with a e- or e+ if e- + e+ stack is not empty
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       else if (nen.gt.0) then 
          EI=DEN(1,NEN)
          TI=DEN(2,NEN)
          dCTSTAI=DEN(4,NEN)
          DLI=DEN8(1,NEN)
          DMI=DEN8(2,NEN)
          DNI=DEN8(3,NEN)
          TYPEI=DEN(5,NEN)
          CT_1I=DEN(6,NEN)
          X=DEN(8,NEN)
          Y=DEN(9,NEN)
          do ith=1,NTHIN+1
             Wgt(ith)=WDEN(NEN,ith)
          end do
          NEN=NEN-1
       end if

c Initial values equal to those of first particle injected       
       DL=DLI
       DM=DMI
       DN=DNI
       E1=EI
       T1=TI
       TYPEP=TYPEI
       ZSTART=TI
       XSTART=X
       YSTART=Y
       dCTSTA=dCTSTAI
       E_A=EI
       CT_1=CT_1I
       CT_2=0.D0 
       DELTA_CT=0.0D0

c Loop to assign weights to primary particle
c in the different thinning levels
       do ith=1,NTHIN+1
         Wgt(ith)=WgtI
       end do

       do ith=NTHIN+2,NTHIN_MAX
         Wgt(ith)=0.
       end do

c Further division of number of substeps 
c in electron and positron propagation
c See WARNING at the header of this code.       
       NSTEPS_FAC=100

C -------------------------------------------------------------------
      IF(TYPEP) 400,1,450
C
C ******************************************************************
C  400 IS THE BLOCK FOR POSITRONS
C  1 - FOR PHOTONS
C  450 - FOR ELECTRONS
C ********************************************************************
C
C ********************************************************************
C ****************************===========*****************************
C **************************** ELECTRONS *****************************
C ****************************===========*****************************
C ********************************************************************
C
  
C --------------------------------------------------------------------
C
C *******************************************************************
C * Electron propagation: 
C           Finds distance to next interaction (bremss, Moeller)
C           Calculates ionization losses along that distance
C           Accounts for electron multiple scattering
C *******************************************************************
C
C    STOT=total x-section, RBRAI=Bremmstrahlung, RMOL=Moeller
C
450   STOT=RBRAI(E1,AP)+IMOELL*RMOL(E1)

      R=s_rndm(0)

C Distance to next interaction 
      T=-ALOG(R)/STOT

C Write depth of first interaction
      if (ifirst.eq.0) then
          ifirst=1
          write(88,*) T*T0 
      end if

      IF (T .NE. 0.) THEN

C                        EK1=kinetic energy
         E_A=E1
         EK1=E1-e_MASS

c Number of steps in which the distance to next interaction is divided
         KL=NSTEPS_FAC*(1+T/( MIN(RG*EK1*(1.-RG2/(1.+RG3*EK1)),TMAX)))

c         write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@' 
c         write(*,*) T,KL
         T=T/KL
c         write(*,*) T,KL

c Loop over steps
      DO I=1,KL

C Ionization losses
        CALL RIONFN(0,E1,T,TE1,DEI,AP)

        IF (E1-CH_THR .LE. DEI) THEN
          DEIPR=E1-CH_THR_EP
          T=T*DEIPR/DEI
          DEI=DEIPR
        END IF      
          E4=E1-DEI
          E8=(E1+E4)/2.
          E_B=E4


c Keep track of absolute times
          G=E8/e_MASS
          G2=G*G
          BETA=SQRT(1.-(1./G2))
  
C Multiple scattering
        CALL MULTSC(E8,TIX,T)
        TOLD=T1
        DNOLD=DN

C Geometry after scattering
        CALL GEOM(E1,DEI,T1,T,TIX,T2,DL,DM,DN,DNHALF,DX,DY,DZ,dCTi)
        DdCT=DdCT+DBLE(dCTi)
        dCT=dCTSTA+DdCT

c Keep track of absolute times
        DELTA_CT_STEP=DBLE(SQRT(DX*DX+DY*DY+DZ*DZ)/BETA)
        DELTA_CT=DELTA_CT+DELTA_CT_STEP
        CT_2=CT_1+DELTA_CT
c        write(99,*) CT_1,CT_2,DZ,T1,SQRT(DX*DX+DY*DY+DZ*DZ),DELTA_CT,
c     #SQRT(DX*DX+DY*DY+DZ*DZ)/BETA

        MXJETH=MIN(INT(XLOG*LOG10((E4-e_MASS)/BR_THR))+1,10)

C Calculate tracklength as a function of electron energy
         do ith=1,NTHIN+1
           DO JJ=1,MXJETH
              TCKSUMNOPRO(JJ,ith)=TCKSUMNOPRO(JJ,ith)-T*Wgt(ith)
              TCKSUMPRO(JJ,ith)=TCKSUMPRO(JJ,ith)-T*DNHALF*Wgt(ith)
              TCKABS(JJ,ith)=TCKABS(JJ,ith)+T*Wgt(ith)
              TCKPRO(JJ,ith)=TCKPRO(JJ,ith)+T*DNHALF*Wgt(ith)
           END DO
         end do

          EKINET=E1-e_MASS
        IF (EKINET .GT. ETHRKE(MXJETH+1)) THEN
         MAXPAR=INT(XLOG*LOG10(EKINET/ETHRKE(MXJETH)))

        do ith=1,NTHIN+1
          DO JJ=MXJETH+1,MIN(MXJETH+MAXPAR,10)
            TFRACT=(EKINET-ETHRKE(JJ))/DEI
            DNAVER=DNOLD+(DN-DNOLD)*MAX(0.,TFRACT-0.5)
            TSEGME=T*TFRACT*Wgt(ith)
            TPROJE=TSEGME*DNAVER
            TCKSUMNOPRO(JJ,ith)=TCKSUMNOPRO(JJ,ith)-TSEGME
            TCKSUMPRO(JJ,ith)=TCKSUMPRO(JJ,ith)-TPROJE
            TCKABS(JJ,ith)=TCKABS(JJ,ith)+TSEGME
            TCKPRO(JJ,ith)=TCKPRO(JJ,ith)+TPROJE
          END DO
        end do

        END IF
c #################################################################
c #################################################################
c #################################################################
c Energy losses as a function of position - for ROOT ntuple
c        write(333,*) 
c     #  (X+XOLD)/2./RHO/100.,  ! 1
c     #  (Y+YOLD)/2./RHO/100.,  ! 2
c     #  (T1+TOLD)/2./RHO/100., ! 3
c     #  (E1-E4)                ! 4
c #################################################################
c #################################################################
c #################################################################

        XOLD=X
        YOLD=Y
        X=X+DX
        Y=Y+DY
c        write(*,*) TOLD,T1
        NROLD=INT(TOLD/T0)
        NRAD=MIN(INT(T1/T0),NRMAX)
          IF (T1 .LE. TOLD) THEN                  ! Back-scattering
            do ith=1,NTHIN+1
               NRE(NROLD,ith)=NRE(NROLD,ith)+Wgt(ith)
              if(Wgt(ith).gt.0.) 
     #        NRENOTHIN(NROLD,ith)=NRENOTHIN(NROLD,ith)+1
            end do
          ELSE IF (NRAD .NE. NROLD) THEN
            CALL DISTRS(1,TOLD,T1,XOLD,YOLD,X,Y,E4,E1,NRAD,NROLD,Wgt)
          END IF
        E1=E4
         
c --------------------------------------------------------------------
c Approx C
        if (APPROX.EQ.APPROX_C) THEN

        do ith=1,NTHIN+1   ! One call to EMPSUM per thinning level

        CT_1_N=CT_2-DELTA_CT_STEP

c Fraunhofer
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),TOLD,T1-TOLD,XOLD,X-XOLD, 
     # YOLD,Y-YOLD,dCT-DdCT,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),TOLD,T1,XOLD,X,YOLD,Y, 
     # CT_1_N,CT_2,E1,E4,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,
     # Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),TOLD,T1,XOLD,X, 
     # YOLD,Y,CT_1_N,CT_2,E1,E4,ASUM_X,ASUM_Y,ASUM_Z,
     # Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),TOLD,T1,XOLD,X,YOLD,Y,
     # CT_1_N,CT_2,E1,E4,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)

c         if (ith.ne.NTHIN+1)
c     #   write(98,*) INT(TYPEP),TOLD,T1,XOLD,X,YOLD,Y,
c     #              CT_1_N,CT_2,E1,E4,Wgt(ith)

         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),TOLD,T1,XOLD,X, 
     # YOLD,Y,CT_1_N,CT_2,E1,E4,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if

        end do

        end if ! Ends approx. C

c --------------------------------------------------------------------

        IF ((E1.LE.CH_THR) .OR. (T1.GE.ZOBS) .OR. (T1.LE.0.)) THEN

c        write(*,*) 'Electron just outside volume ',T1

        if (APPROX.NE.APPROX_C) THEN

        do ith=1,NTHIN+1   ! One call to EMPSUM per thinning level

c Fraunhofer
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART, 
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     # CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


           IF (ith.eq.NTHIN+1) THEN 
               CALL SIZE_GREISEN(INT(TYPEP),E1,T1)   ! Hybrid shower
               CALL TRACKLENGTHS(E1)           
           END IF

        end do

        end if  

        GOTO 66  

        END IF

      END DO   ! Ends loop in KL

      IF (APPROX .EQ. APPROX_B2) THEN   ! In Approx_b2 EMPSUM is called even if
                                        ! the interaction is not real

        do ith=1,NTHIN+1    ! One call to EMPSUM per thinning level

c Fraunhofer
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART,
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


        end do

            XSTART=X
            YSTART=Y
            ZSTART=T1
            dCTSTA=dCT
            DdCT=0.0D0
            CT_1=CT_2     ! New track begins => update absolute time of particle
            DELTA_CT=0.0D0

      END IF


C ---------------------------------------------------------------------
C Find if interaction actually happens
C Takes into account that the energy of the electron decreases 
C as it propagates and the interaction at the new energy might
C not occur (only valid for cross sections that increase with energy)
C ---------------------------------------------------------------------
      ST1=RBRAI(E1,AP)
      ST2=IMOELL*RMOL(E1)
      ST=ST1+ST2
      R=s_rndm(0)
      IF(R.GT.(ST/STOT)) GO TO 450
C ---------------------------------------------------------------------
C  450 - THE INTERACTION IS NOT REAL
C ---------------------------------------------------------------------
      END IF
      
c      IF THE INTERACTION OCCURS, EMPSUM MUST BE CALLED IN APPROXIMATION
C      B - THE TYPE OF INTERACTION DOES NOT MATTER!!!!!!!!!      
C        APPROXIMATION B SHOULD BE CALLED HERE - WE KNOW THE FINAL DETAILS!


      IF (APPROX .EQ. APPROX_B1) THEN
       do ith=1,NTHIN+1     ! One call to EMPSUM per thinning level

c Fraunhofer
          if (ifield_w.eq.1) then
            CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART, 
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


       end do

            XSTART=X
            YSTART=Y
            ZSTART=T1
            dCTSTA=dCT
            DdCT=0.0D0
            CT_1=CT_2   
            DELTA_CT=0.D0

      END IF

      R=s_rndm(0)
      IF(R.GE.(ST1/ST)) THEN
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ------------------------ MOELLER SCATERING --------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
        CALL RMOLS(E1,E3,V3,V4)
        E2=E1-E3+e_mass

        E2_AV=MAX(E2-THIN_METH*CH_THR,0.)
        E3_AV=MAX(E3-THIN_METH*CH_THR,0.)

        CALL THINNING
     #(1,E2_AV,T1,1,E3_AV,T1,Wgt2,Wgt3,KEEP2,KEEP3,Wgt)
        CALL GEONEW(DL,DM,DN,V3,V4,DL3,DM3,DN3,DL,DM,DN)    

C Loop over thinning levels to count in how many the particle is rejected
        irejected2=0
        irejected3=0
        do ith=1,NTHIN+1

        IF ((KEEP2(ith) .EQ. 0)) irejected2=irejected2+1
        IF ((KEEP3(ith) .EQ. 0)) irejected3=irejected3+1

C If neither particle is being kept at a particular thinning level 
c sum old track in approx A 	
        IF ((KEEP2(ith) .EQ. 0).AND.(KEEP3(ith).eq.0)) THEN 

          IF (APPROX .EQ. APPROX_A) THEN    ! One call to EMPSUM per thinning level

c Fraunhofer    
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART,
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


c  No resetting of x,y,z etc starts needed as we are going to 66
          END IF

        END IF

        end do

c If the 2 particles are rejected at all thinning levels select new particle
          if (irejected2.eq.NTHIN+1.and.irejected3.eq.NTHIN+1) GOTO 66

C
C *******************************************************************
C Stack "new" electron - arbitrarily picked to be the third
C *******************************************************************
          if (irejected3.ne.NTHIN+1) then
          TYPEP=1.
          NEN=NEN+1
          nen2=nen2+1
          IF(NEN.GT.MAXNESTACK) GO TO 9501
            DEN(1,NEN)=E3
            DEN(2,NEN)=T1
            DEN(4,NEN)=dCT
            DEN(5,NEN)=TYPEP
            DEN(6,NEN)=CT_2
            DEN8(1,NEN)=DL3
            DEN8(2,NEN)=DM3
            DEN8(3,NEN)=DN3
            DEN(8,NEN)=X
            DEN(9,NEN)=Y
            do ith=1,NTHIN+1
              if (keep3(ith).eq.0) then     ! Update weights and stack particle
                 WDEN(NEN,ith)=0.           ! Assign 0 weight if particle is rejected
                                            ! in a particular thinning level
              else
                 WDEN(NEN,ith)=Wgt(ith)*Wgt3(ith)
              end if
            end do
          end if    
C
C ---------------------------------------------------------------------
C ------------------------- THE SECOND ELECTRON -----------------------
C ---------------------------------------------------------------------
C      NOTE: the one being kept is always called the 'old' one due to 
c indistinguishability. 

c If we get here it means that at least one keep2 is not 0 and we should follow the electron
           E1=E2
           do ith=1,NTHIN+1
             if (keep2(ith).eq.0) then
               Wgt(ith)=0.
             else
               Wgt(ith)=Wgt(ith)*Wgt2(ith)
             end if
           end do
           GO TO 450

      ELSE
C
C ---------------------------------------------------------------------
C ------------------------ BREMMSTRAHLUNG --------------------------
C ---------------------------------------------------------------------
C

        N1=1.442695*ALOG(E1/AP)+1
C  9200   CONTINUE
C ----------------------------------------------------------------------
C Differential bremss x-sections depending on the energy of the electron
C ----------------------------------------------------------------------
C
        IF (E1 .LE. 50.) THEN
          CALL RBRSME(E1,U,N1,NO,ALF)
        ELSE IF (E1 .LE. E_LPM_BREMSS) THEN
          CALL RBRMFC(E1,U,N1,NO,ALF)
        ELSE
          CALL RBRSMS(E1,U,N1,NO,ALF)  
        END IF
C ---------------------------------------------------------------------
C *********************************************************************
C  E2 IS THE ENERGY OF THE PHOTON
C  E3 IS THE ENERGY OF THE ELECTRON
C    V2 - COS  OF THE POLAR ANGLE OF THE PHOTON 
C    V3 - COS  OF THE POLAR ANGLE OF THE ELECTRON 
C *********************************************************************
C ---------------------------------------------------------------------
        E2=E1*U
        E3=E1-E2
        V2=COS(e_MASS/E1)
        V3=1.0

        E2_AV=MAX(E2-THIN_METH*TRP,0.)
        E3_AV=MAX(E3-THIN_METH*CH_THR,0.)

        CALL THINNING
     #(0,E2_AV,T1,1,E3_AV,T1,Wgt2,Wgt3,KEEP2,KEEP3,Wgt)
        CALL GEONEW(DL,DM,DN,V2,V3,DL2,DM2,DN2,DL,DM,DN)

C  The gamma, if kept, will always be stacked no matter what	
       istack2=0
       do ith=1,NTHIN+1
c If at least one keep2 is 1 the gamma is stacked
         IF(KEEP2(ith) .EQ. 1) istack2=1
       end do

c Stack photon if at least one of the keep2 is equal to 1
       if (istack2.eq.1) then
          NGN=NGN+1
          IF(NGN.GT.MAXNGSTACK) GO TO 9501
            DGN(1,NGN)=E2
            DGN(2,NGN)=T1
            DGN(4,NGN)=dCT
            DGN(6,NGN)=CT_2
            DGN8(1,NGN)=DL2
            DGN8(2,NGN)=DM2
            DGN8(3,NGN)=DN2
            DGN(8,NGN)=X
            DGN(9,NGN)=Y
            do ith=1,NTHIN+1
              if (keep2(ith).eq.0) then
                 WDGN(NGN,ith)=0.
              else
                 WDGN(NGN,ith)=Wgt(ith)*Wgt2(ith)
              end if 
            end do
        end if 
       
c Electron
       irejected3=0
       do ith=1,NTHIN+1
c For each thinning level in which the electron is rejected calculate field in Approx. A
       IF (KEEP3(ith) .EQ. 0) THEN 

          irejected3=irejected3+1    ! Count in how many thinning levels the electron is rejected
          IF (APPROX .EQ. APPROX_A) THEN      ! One call to EMPSUM per thinning level

c Fraunhofer    
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART, 
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


          END IF

       END IF
       end do

c If the electron is rejected in all thinning levels pick other particle 
          if (irejected3.eq.NTHIN+1) GOTO 66

!       END IF

c If we get here it means that at least one of the keep3 is not 0
C  CONTINUE WITH ELECTRON
          E1=E3
          do ith=1,NTHIN+1
            if (keep3(ith).eq.0) then
               Wgt(ith)=0.
            else 
               Wgt(ith)=Wgt(ith)*Wgt3(ith)
            end if
          end do
          GO TO 450

      END IF
C
C********************************************************************
C****************************===========*****************************
C**************************** POSITRONS *****************************
C****************************===========*****************************
C********************************************************************
C
C Same comments as in the electron block 
C
  400 CONTINUE
C -------------------------------------------------------------------
C      IF(E1.GT.1000.) GOTO 2700
C -------------------------------------------------------------------
C
C                        RBHA=Scattering Bhabha
C
      STOT=RBRAI(E1,AP)+IBHABH*RBHA(E1)+IANNHI*RANH(E1)
c      write(77,*) E1,rbrai(e1,ap),rbha(e1),ranh(e1)
      R=s_rndm(0)
      T=-ALOG(R)/STOT

C Write depth of first interaction
      if (ifirst.eq.0) then
          ifirst=1
          write(88,*) T*T0 
      end if

      IF (T .NE. 0.) THEN

         E_A=E1
         EK1=E1-e_MASS
         KL=NSTEPS_FAC*(1+T/( MIN(RG*EK1*(1.-RG2/(1.+RG3*EK1)),TMAX)))
         T=T/KL

      DO I=1,KL
        CALL RIONFN(1,E1,T,TE1,DEI,AP)
          IF (E1-CH_THR .LE. DEI) THEN
            DEIPR=E1-CH_THR_EP
            T=T*DEIPR/DEI
            DEI=DEIPR
          END IF      
c        write(97,*) e1,dei/(t*t0)
         E4=E1-DEI
         E8=(E1+E4)/2.
         E_B=E4

c Positron times
         G=E8/e_MASS
         G2=G*G
         BETA=sqrt(1.-1./G2)

         CALL MULTSC(E8,TIX,T)
         TOLD=T1
         DNOLD=DN
         CALL GEOM(E1,DEI,T1,T,TIX,T2,DL,DM,DN,DNHALF,DX,DY,DZ,dCTi)
         DdCT=DdCT+DBLE(dCTi)
         dCT=dCTSTA+DdCT
         DELTA_CT_STEP=DBLE(SQRT(DX*DX+DY*DY+DZ*DZ)/BETA)
         DELTA_CT=DELTA_CT+DELTA_CT_STEP
         CT_2=CT_1+DELTA_CT
   
         MXJETH=MIN(INT(XLOG*LOG10((E4-e_MASS)/BR_THR))+1,10)
          do ith=1,NTHIN+1
            DO JJ=1,MXJETH
              TCKSUMNOPRO(JJ,ith)=TCKSUMNOPRO(JJ,ith)+T*Wgt(ith)
              TCKSUMPRO(JJ,ith)=TCKSUMPRO(JJ,ith)+T*DNHALF*Wgt(ith)
              TCKABS(JJ,ith)=TCKABS(JJ,ith)+T*Wgt(ith)
              TCKPRO(JJ,ith)=TCKPRO(JJ,ith)+T*DNHALF*Wgt(ith)
            END DO
          end do

        EKINET=E1-e_MASS
        IF (EKINET .GT. ETHRKE(MXJETH+1)) THEN
          MAXPAR=INT(XLOG*LOG10(EKINET/ETHRKE(MXJETH)))
          do ith=1,NTHIN+1
            DO JJ=MXJETH+1,MIN(MXJETH+MAXPAR,10)
             TFRACT=(EKINET-ETHRKE(JJ))/DEI
             DNAVER=DNOLD+(DN-DNOLD)*MAX(0.,TFRACT-0.5)
             TSEGME=T*TFRACT*Wgt(ith)
             TPROJE=TSEGME*DNAVER
             TCKSUMNOPRO(JJ,ith)=TCKSUMNOPRO(JJ,ith)+TSEGME
             TCKSUMPRO(JJ,ith)=TCKSUMPRO(JJ,ith)+TPROJE
             TCKABS(JJ,ith)=TCKABS(JJ,ith)+TSEGME
             TCKPRO(JJ,ith)=TCKPRO(JJ,ith)+TPROJE
            END DO
          end do
        END IF

C ------------------------------------------------------------------
c #################################################################
c #################################################################
c #################################################################
c Energy losses as a function of position - for ROOT ntuple
c        write(333,*) 
c     #  (X+XOLD)/2./RHO/100.,  ! 1
c     #  (Y+YOLD)/2./RHO/100.,  ! 2
c     #  (T1+TOLD)/2./RHO/100., ! 3
c     #  (E1-E4)                ! 4
c #################################################################
c #################################################################
c #################################################################
C ------------------------------------------------------------------


        XOLD=X
        YOLD=Y
        X=X+DX
        Y=Y+DY
        NROLD=INT(TOLD/T0)
        NRAD=MIN(INT(T1/T0),NRMAX)
          IF (T1 .LE. TOLD) THEN ! Back-scattering
            do ith=1,NTHIN+1
              NRPOS(NROLD,ith)=NRPOS(NROLD,ith)+Wgt(ith)
              if(Wgt(ith).gt.0.) 
     #        NRPOSNOTHIN(NROLD,ith)=NRPOSNOTHIN(NROLD,ith)+1
            end do
          ELSE IF (NRAD .NE. NROLD) THEN
            CALL DISTRS(-1,TOLD,T1,XOLD,YOLD,X,Y,E4,E1,NRAD,NROLD,Wgt)
          END IF
        E1=E4

c --------------------------------------------------------------------
c Approx C
        if (APPROX.EQ.APPROX_C) THEN

        do ith=1,NTHIN+1   ! One call to EMPSUM per thinning level

        CT_1_N=CT_2-DELTA_CT_STEP

c Fraunhofer
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),TOLD,T1-TOLD,XOLD,X-XOLD, 
     # YOLD,Y-YOLD,dCT-DdCT,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),TOLD,T1,XOLD,X,YOLD,Y, 
     # CT_1_N,CT_2,E1,E4,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,
     # Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),TOLD,T1,XOLD,X, 
     # YOLD,Y,CT_1_N,CT_2,E1,E4,ASUM_X,ASUM_Y,ASUM_Z,
     # Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),TOLD,T1,XOLD,X,YOLD,Y,
     # CT_1_N,CT_2,E1,E4,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

c         if (ith.ne.NTHIN+1)
c         write(98,*) INT(TYPEP),TOLD,T1,XOLD,X,YOLD,Y,
c     #              CT_1_N,CT_2,E1,E4,Wgt(ith)

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),TOLD,T1,XOLD,X, 
     # YOLD,Y,CT_1_N,CT_2,E1,E4,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if

        end do

        end if ! Ends approx. C


C ---------------------------------------------------------------------	
        IF ((E1.LT.CH_THR) .OR. (T1.GE.ZOBS) .OR. (T1.LE.0.)) THEN

c        write(*,*) 'Positron just outside volume ',T1

        if (APPROX.NE.APPROX_C) THEN

        do ith=1,NTHIN+1       ! One call to EMPSUM per thinning level

c Fraunhofer
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART, 
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


        end do

        end if

        GOTO 66

        END IF

      END DO
C ----------------------------------------------------------------

      IF (APPROX .EQ. APPROX_B2) THEN   ! In Approx_b2 EMPSUM is called even if
                                        ! the interaction is not real

        do ith=1,NTHIN+1        ! One call to EMPSUM per thinning level

c Fraunhofer        
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART, 
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


        end do
            XSTART=X
            YSTART=Y
            ZSTART=T1
            dCTSTA=dCT
            DdCT=0.0D0
            CT_1=CT_2
            DELTA_CT=0.0D0
      END IF

      ST1=RBRAI(E1,AP)
      ST2=IBHABH*RBHA(E1)   
      ST3=IANNHI*RANH(E1)
      ST=ST1+ST2+ST3     
      R=s_rndm(0)
      IF(R.GT.(ST/STOT)) GO TO 400
C ----------------------------------------------------------------
C  400 - THE INTERACTION IS NOT REAL
C        CHOISE OF INTERACTION
C ----------------------------------------------------------------
      END IF
      
      IF (APPROX .EQ. APPROX_B1) THEN

       do ith=1,NTHIN+1       ! One call to EMPSUM per thinning level

c Fraunhofer
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART, 
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


       end do

            XSTART=X
            YSTART=Y
            ZSTART=T1
            dCTSTA=dCT
            DdCT=0.0D0
            CT_1=CT_2
            DELTA_CT=0.0D0
      END IF
      
      R=s_rndm(0)
      IF(R.LT.(ST1/ST)) THEN
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ------------------------ BREMSSTRAHLUNG ------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
        N1=1.442695*ALOG(E1/AP)+1

        IF (E1 .LE. 50.) THEN
          CALL RBRSME(E1,U,N1,NO,ALF)
        ELSE IF (E1 .LE. E_LPM_BREMSS) THEN
          CALL RBRMFC(E1,U,N1,NO,ALF)
        ELSE
          CALL RBRSMS(E1,U,N1,NO,ALF)
        END IF

C ---------------------------------------------------------------------
C *********************************************************************
C  E2 IS THE ENERGY OF THE EMITTED PHOTON
C  E3 IS THE POSITRON ENERGY
C    V2 - COS  OF THE POLAR ANGLE OF THE PHOTON
C    V3 - COS  OF THE POLAR ANGLE OF THE POSITRON
C *********************************************************************
C ---------------------------------------------------------------------
        E2=E1*U
        E3=E1-E2
        V2=COS(e_MASS/E1)
        V3=1.0
        CALL GEONEW(DL,DM,DN,V2,V3,DL2,DM2,DN2,DL,DM,DN)

        E2_AV=MAX(E2-THIN_METH*TRP,0.)
        E3_AV=MAX(E3-THIN_METH*CH_THR,0.)

        CALL THINNING
     #(0,E2_AV,T1,1,E3_AV,T1,Wgt2,Wgt3,KEEP2,KEEP3,Wgt)

C Stack the gamma if it is being kept i.e. if at least one of the keep2 is equal to 1
      istack2=0
      do ith=1,NTHIN+1
        IF (KEEP2(ith) .EQ. 1) istack2=1
      end do

      if (istack2.eq.1) then 
          NGN=NGN+1
          IF(NGN.GT.MAXNGSTACK) GO TO 9501
            DGN(1,NGN)=E2
            DGN(2,NGN)=T1
            DGN(4,NGN)=dCT
            DGN(6,NGN)=CT_2
            DGN8(1,NGN)=DL2
            DGN8(2,NGN)=DM2
            DGN8(3,NGN)=DN2
            DGN(8,NGN)=X
            DGN(9,NGN)=Y 
            do ith=1,NTHIN+1
              if (keep2(ith).eq.0) then
                WDGN(NGN,ith)=0.
              else
                WDGN(NGN,ith)=Wgt(ith)*Wgt2(ith)
              end if 
            end do
      end if

c The positron
        irejected3=0
        do ith=1,NTHIN+1
        IF (KEEP3(ith) .EQ. 0) THEN 
          irejected3=irejected3+1      ! Count in how many thinning levels the particle is rejected
          IF (APPROX .EQ. APPROX_A) THEN   ! One call to EMPSUM per thinning level

c Fraunhofer    
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART, 
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


          END IF
        END IF
        end do

c If positron is rejected for all thinning levels get another particle
        if (irejected3.eq.NTHIN+1) GOTO 66

          E1=E3
          do ith=1,NTHIN+1
             if (keep3(ith).eq.0) then
               Wgt(ith)=0.
             else
               Wgt(ith)=Wgt(ith)*Wgt3(ith)
             end if 
          end do
          GO TO 400
!jam       END IF
      
       ELSE IF (R.LT.((ST1+ST2)/ST)) THEN
C
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C ---------------------- BHABHA SCATTERING ----------------------------
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
C
        CALL RBHAS(E1,E3,EE,V3,V4)
c	        CALL RBHAS(E1,AE,E3,EE,V3,V4)
C ---------------------------------------------------------------------
C  E3  IS THE ENERGY OF THE SCATTERED POSITRON
C  EE  IS THE ENERGY OF THE EMITTED ELECTRON
C ---------------------------------------------------------------------

        CALL GEONEW(DL,DM,DN,V3,V4,DL,DM,DN,DLE,DME,DNE)
        EE=E1-E3+e_mass

        EE_AV=MAX(0.,EE-THIN_METH*CH_THR)
        E3_AV=MAX(E3-THIN_METH*CH_THR,0.)

        CALL THINNING
     #(1,EE_AV,T1,1,E3_AV,T1,WgtE,Wgt3,KEEPE,KEEP3,Wgt)


        irejectedp=0
        do ith=1,NTHIN+1
        IF (KEEP3(ith) .EQ. 0) THEN
          irejectedp=irejectedp+1       ! Count in how many thinning levels the e+ is rejected
          IF (APPROX .EQ. APPROX_A) THEN     ! One call to EMPSUM per thinning level

c Fraunhofer    
          if (ifield_w.eq.1) then
           CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART, 
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
          end if

          if (ifield_W_2.eq.1) then
           CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
          end if

          if (ifield_T.eq.1) then
           CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
          end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


          END IF
        END IF
        end do  

        irejectede=0
        do ith=1,NTHIN+1
c Count in how many thinning levels the e- is rejected
          IF (KEEPE(ith).EQ. 0) irejectede=irejectede+1     
        end do

c Both electron and positron rejected for all thinning levels get another particle
        if (irejectedp.eq.NTHIN+1.and.irejectede.eq.NTHIN+1) goto 66 

c Extract electron

c If the positron is rejected for all thinning levels propagate electron
        if (irejectedp.eq.NTHIN+1) then 

            XSTART=X
            YSTART=Y
            ZSTART=T1
            dCTSTA=dCT
            DdCT=0.0D0
            CT_1=CT_2
            DELTA_CT=0.0D0
            do ith=1,NTHIN+1
              if (keepe(ith).eq.0) then
                 Wgt(ith)=0.
              else
                 Wgt(ith)=Wgt(ith)*WgtE(ith)
              end if
            end do
            E1=EE
            TYPEP=1.
            DL=DLE
            DM=DME
            DN=DNE
            GO TO 450

        else

c store the electron if it is not rejected for all thinning levels... 
          IF (irejectede .ne. NTHIN+1) THEN
            NEN=NEN+1
            nen2=nen2+1
            IF(NEN.GT.MAXNESTACK) GO TO 9501
              DEN(1,NEN)=EE
              DEN(2,NEN)=T1
              DEN(4,NEN)=dCT
              DEN(5,NEN)=1.
              DEN(6,NEN)=CT_2
              DEN8(1,NEN)=DLE
              DEN8(2,NEN)=DME
              DEN8(3,NEN)=DNE
              DEN(8,NEN)=X
              DEN(9,NEN)=Y
              do ith=1,NTHIN+1
                if (keepe(ith).eq.0) then
                  WDEN(NEN,ith)=0.
                else
                  WDEN(NEN,ith)=Wgt(ith)*WgtE(ith)
                end if
              end do
          END IF

c ...and propagate positron
           E1=E3
           do ith=1,NTHIN+1
              if (keep3(ith).eq.0) then
                 Wgt(ith)=0.
              else 
                 Wgt(ith)=Wgt(ith)*Wgt3(ith)
              end if
           end do
           GO TO 400

        end if

      ELSE
C
C --------------------------------------------------------------------
C --------------------------------------------------------------------
C -------------------------- ANNIHILATION ----------------------------
C --------------------------------------------------------------------
C --------------------------------------------------------------------
C
        CALL RANHS(E1,E2,ESM,CBI,CSM)
c        IF(CBI.GE.1.) WRITE(*,*) E1,E2,ESM,CBI,CSM
        CALL GEONEW(DL,DM,DN,CBI,CSM,DL3,DM3,DN3,DL4,DM4,DN4)

        E2_AV=MAX(E2-THIN_METH*TRP,0.)
        E3_AV=MAX(ESM-THIN_METH*TRP,0.)

        CALL THINNING
     #(0,E2_AV,T1,0,E3_AV,T1,Wgt2,Wgt3,KEEP2, KEEP3,Wgt)

C always need to call empsum in approx a!
        do ith=1,NTHIN+1        ! One call to EMPSUM per thinning level
           IF (APPROX .EQ. APPROX_A) THEN 

c Fraunhofer   
            if (ifield_w.eq.1) then
             CALL EMPSUM_W (INT(TYPEP),ZSTART,T1-ZSTART,XSTART,X-XSTART,
     # YSTART,Y-YSTART,dCTSTA,DdCT,
     # ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt(ith),ith)
            end if

            if (ifield_W_2.eq.1) then
             CALL EMPSUM_W_2 (INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y, 
     #CT_1,CT_2,E_A,E_B,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt(ith),ith)
            end if

            if (ifield_T.eq.1) then
             CALL EMPSUM_T (INT(TYPEP),ZSTART,T1,XSTART,X, 
     #YSTART,Y,CT_1,CT_2,E_A,E_B,ASUM_X,ASUM_Y,ASUM_Z,Wgt(ith),ith)
            end if

c Fresnel
         if (ifield_W_2_fresnel.eq.1) then
         CALL EMPSUM_W_2_FRESNEL(INT(TYPEP),ZSTART,T1,XSTART,X,YSTART,Y,
     # CT_1,CT_2,E_A,E_B,
     # ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,Wgt(ith),ith)
         end if

          if (ifield_T_fresnel.eq.1) then
          CALL EMPSUM_T_FRESNEL (INT(TYPEP),ZSTART,T1,XSTART,X, 
     # YSTART,Y,CT_1,CT_2,E_A,E_B,
     # ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,Wgt(ith),ith)
          end if


           END IF
        end do


c First Photon
        istack2=0
        do ith=1,NTHIN+1
           IF (KEEP2(ith) .EQ. 1) istack2=1
        end do

        if (istack2.eq.1) then 
          NGN=NGN+1
          IF(NGN.GT.MAXNGSTACK) GO TO 9501
          DGN(1,NGN)=E2
          DGN(2,NGN)=T1
          DGN(4,NGN)=dCT
          DGN(6,NGN)=CT_2
          DGN8(1,NGN)=DL3
          DGN8(2,NGN)=DM3
          DGN8(3,NGN)=DN3
          DGN(8,NGN)=X
          DGN(9,NGN)=Y
          do ith=1,NTHIN+1
             if (keep2(ith).eq.0) then
                WDGN(NGN,ith)=0.
             else
                WDGN(NGN,ith)=Wgt(ith)*Wgt2(ith)
             end if
          end do
        end if

c Second Photon
        istack3=0
        do ith=1,NTHIN+1
           IF (KEEP3(ith).eq.1) istack3=1
        end do

        if (istack3.eq.1) then
          NGN=NGN+1
          IF(NGN.GT.MAXNGSTACK) GO TO 9501
            DGN(1,NGN)=ESM
            DGN(2,NGN)=T1
            DGN(4,NGN)=dCT
            DGN(6,NGN)=CT_2
            DGN8(1,NGN)=DL4
            DGN8(2,NGN)=DM4
            DGN8(3,NGN)=DN4
            DGN(8,NGN)=X  
            DGN(9,NGN)=Y
            do ith=1,NTHIN+1
               if (keep3(ith).eq.0) then
                  WDGN(NGN,ith)=0.
               else
                  WDGN(NGN,ith)=Wgt(ith)*Wgt3(ith)
               end if
            end do
        end if

c Pick next particle
        GOTO 66
      END IF

C
C --------------------------------------------------------------------
C --------------------------------------------------------------------
C ---------------------- ONE PHOTON ANNIHILATION ---------------------
C --------------------------------------------------------------------
C --------------------------------------------------------------------
C
C  500 E2=E1+e_MASS
C          NGN=NGN+1
C          IF(NGN.GT.10000) GO TO 9501
C            DGN(1,NGN)=E2
C            DGN(2,NGN)=T1
C            DGN(4,NGN)=dCT
C            DGN8(1,NGN)=DL3
C            DGN8(2,NGN)=DM3
C            DGN8(3,NGN)=DN3
C            DGN(8,NGN)=X
C            DGN(9,NGN)=Y
C      GO TO 6

C *** Extract photons from the stack
   24 IF(NGN.LT.1) STOP 9999
      IF(NGN .GT. MAX_NGN) MAX_NGN=NGN
        E1=DGN(1,NGN)
        T1=DGN(2,NGN)
        dCT=DGN(4,NGN)
        DL=DGN8(1,NGN)
        DM=DGN8(2,NGN)
        DN=DGN8(3,NGN)
        X=DGN(8,NGN)
        Y=DGN(9,NGN)
        CT_1=DGN(6,NGN)
        DELTA_CT=0.0D0
        do ith=1,NTHIN+1
           Wgt(ith)=WDGN(NGN,ith)
        end do
        NGN=NGN-1
C
C *********************************************************************
C *****************************=========*******************************
C ***************************** PHOTONS *******************************
C *****************************=========*******************************
C *********************************************************************
C
    1 ST1=RPAIR(E1)
      ST2=ICOMPT*RCOMP(E1)

C      ST3=0.
C ---------------------------------------------------------------------
C     ST3=RPHOT(E1)       
C NOTE THAT FOR AIR ST3=0.
C ---------------------------------------------------------------------
      ST=ST1+ST2
C ---------------------------------------------------------------------
      R=s_rndm(0)
      T=-ALOG(R)/ST

C Write depth of first interaction
      if (ifirst.eq.0) then
          ifirst=1
          write(88,*) T*T0 
      end if

      IF (T .NE. 0.) THEN
      CALL GEOGAM(T1,T,DL,DM,DN,DX,DY,dCT)
      X=X+DX
      Y=Y+DY
      CT_2=CT_1+ABS(T)*T0    ! Note BETA = 1

c      do ith=1,NTHIN
c        if (IWRITE_W_2.eq.1) then  ! Info on tracks written in separate files for each thinning level
c          write(59+ith,31) 
c     # ityp,Wgt0,E1,x-dx,y-dy,t1-(t*dn*t0),E1,x,y,t1,CT_1,CT_2
c31    format(1X,I5,1P,11E18.8)
c        end if

c        if (IWRITE_T.eq.1) then  ! Info on tracks written in separate files for each thinning level
c          write(69+ith,31) 
c     # ityp,Wgt0,E1,x-dx,y-dy,t1-(t*dn*t0),E1,x,y,t1,CT_1,CT_2
c        end if

c      end do


C      NROLD=INT(TOLD/T0)
C      NRAD=MIN(INT(T1/T0),NRMAX)
C      IF (T1 .LE. TOLD) THEN ! Back-scattering
C        NRG=NRG+1
C        GO TO 6
C      ELSE IF (NRAD .NE. NROLD) THEN
C       CALL DISTRS(0,TOLD,T1,XOLD,YOLD,X,Y,E4,E1,NRAD,NROLD)
C      END IF
      IF ((T1.GE.ZOBS) .OR. (T1.LE.0.)) GOTO 66
      END IF
      R=s_rndm(0)
C     ST1=ST-ST2-ST3
C      ST1=ST-ST2 
C -------------------------------------------------------------------
      IF(R.LT.(ST1/ST)) THEN
C
C -------------------------------------------------------------------
C -------------------------------------------------------------------
C ----------------------- PAIR PRODUCTION ---------------------------
C -------------------------------------------------------------------
C -------------------------------------------------------------------
C
c Corrected to include LPM in pair differential distribution
c
        IF(E1.GT. E_LPM) THEN
          call rpaims(e1,u,nO,alf)
          e3=e1*u
        ELSE IF (e1.gt.50.) THEN
        CALL RPAMFC(E1,U,NO,ALF)
          E3=E1*U
        ELSE IF(E1.LT.4.) THEN 
          CALL RPAI4(E1,E3,AP)
        ELSE
          CALL RPAIME(E1,U,NO,ALF)
c18        E3=E1*U
          E3=E1*U
        END IF
C -------------------------------------------------------------------
        E4=E1-E3
        V3=COS(e_MASS/E1)
        V4=V3
        CALL GEONEW(DL,DM,DN,V3,V4,DL3,DM3,DN3,DL,DM,DN)

        E4_AV=MAX(0.,E4-THIN_METH*CH_THR)
        E3_AV=MAX(E3-THIN_METH*CH_THR,0.)

        CALL THINNING
     #(1,E4_AV,T1,1,E3_AV,T1,Wgt4,Wgt3,KEEP4,KEEP3,Wgt)
        

C -------------------------------------------------------------------
        irejected3=0
        irejected4=0
        do ith=1,NTHIN+1
           IF (KEEP3(ith) .EQ. 0) irejected3=irejected3+1 
           IF (KEEP4(ith) .EQ. 0) irejected4=irejected4+1 
        end do

c Both electron and positron are rejected in all thinning levels
        if (irejected3.eq.NTHIN+1.and.irejected4.eq.NTHIN+1) goto 66

C -------------------------------------------------------------------
C  TRACKING SOMETHING - GENERATE A TYPEP FOR E3
        R=s_rndm(0)
        IF(R .GT. .5) THEN
          TYPEP = -1.
        ELSE
          TYPEP = 1.
        END IF



      if (irejected4.ne.NTHIN+1) then
C -------------------------- TRACKING E4 ------------------------------   

        if (irejected3.ne.NTHIN+1) then
C ------------------------- STACKING E3 ---------------------------	 
          NEN=NEN+1

          IF (TYPEP .EQ. 1) THEN
             nen2=nen2+1
          ELSE
             npn2=npn2+1
          END IF

          IF(NEN.GT.MAXNESTACK) GO TO 9501
            DEN(1,NEN)=E3
            DEN(2,NEN)=T1
            DEN(4,NEN)=dCT
            DEN(5,NEN)=TYPEP
            DEN(6,NEN)=CT_2
            DEN8(1,NEN)=DL3
            DEN8(2,NEN)=DM3
            DEN8(3,NEN)=DN3
            DEN(8,NEN)=X
            DEN(9,NEN)=Y
            do ith=1,NTHIN+1
               if (keep3(ith).eq.0) then
                  WDEN(NEN,ith)=0.
               else
                  WDEN(NEN,ith)=Wgt(ith)*Wgt3(ith)
               end if
            end do

        end if


c Track E4
        E1=E4
        TYPEP=-TYPEP

        IF (TYPEP .EQ. 1) THEN
           nen2=nen2+1
        ELSE
           npn2=npn2+1
        END IF

        XSTART=X
        YSTART=Y
        ZSTART=T1
        dCTSTA=dCT
        DdCT=0.D0
        CT_1=CT_2
        DELTA_CT=0.0D0
        do ith=1,NTHIN+1
           if (keep4(ith).eq.0) then
              Wgt(ith) = 0. 
           else
              Wgt(ith) = Wgt(ith)*Wgt4(ith)
           end if
        end do

        IF (TYPEP) 400,1,450

      else

C -------------------------- TRACKING E3 ------------------------------         
        E1=E3
        IF (TYPEP .EQ. 1) THEN
           nen2=nen2+1
        ELSE
           npn2=npn2+1
        END IF

        XSTART=X
        YSTART=Y
        ZSTART=T1
        DL=DL3
        DM=DM3
        DN=DN3
        dCTSTA=dCT
        DdCT=0.D0
        CT_1=CT_2
        DELTA_CT=0.0D0
        do ith=1,NTHIN+1
           if (keep3(ith).eq.0) then
              Wgt(ith) = 0. 
           else
              Wgt(ith) = Wgt(ith)*Wgt3(ith)
           end if
        end do

        IF (TYPEP) 400,1,450
       END IF

      ELSE
C
C -------------------------------------------------------------------
C -------------------------------------------------------------------
C ------------------------- COMPTON  EFFECT -------------------------
C -------------------------------------------------------------------
C -------------------------------------------------------------------
C
C        N2=1.443*ALOG(3.914*E1+1.) +2
        CALL COMPTU(E1,U,V2,N2,NC,ALF)
C -------------------------------------------------------------------
C  E2,V2 - ENERGY AND COS OF THE POLAR ANGLE OF SCATTERED PHOTON
C  E3,V3 - THE SAME FOR THE ELECTRON
C -------------------------------------------------------------------
        E2=E1*U
        E3=E1-E2+e_MASS
C        P1=E1*E1
C        P2=E2*E2
C        P3=E3*E3-e2_MAS
C        V2=COS(TIX)
        RATMAS=e_MASS/E1
        V3=(1.+RATMAS)*SQRT( (1.-U)/(1-U+2.*RATMAS) )
        CALL GEONEW(DL,DM,DN,V2,V3,DL2,DM2,DN2,DL,DM,DN)

        E2_AV=MAX(E2-THIN_METH*TRP,0.)
        E3_AV=MAX(E3-THIN_METH*CH_THR,0.)

        CALL THINNING
     #(0,E2_AV,T1,1,E3_AV,T1,Wgt2,Wgt3,KEEP2,KEEP3,Wgt)

C  keeping at least the photon
        istack2=0
        do ith=1,NTHIN+1
          IF (KEEP2(ith) .EQ. 1) istack2=1
        end do        

        if (istack2.eq.1) then
            NGN=NGN+1
            IF(NGN.GT.MAXNGSTACK) GO TO 9501
              DGN(1,NGN)=E2
              DGN(2,NGN)=T1
              DGN(4,NGN)=dCT
              DGN(6,NGN)=CT_2
              DGN8(1,NGN)=DL2
              DGN8(2,NGN)=DM2
              DGN8(3,NGN)=DN2
              DGN(8,NGN)=X
              DGN(9,NGN)=Y
              do ith=1,NTHIN+1
                 if (keep2(ith).eq.0) then
                    WDGN(NGN,ith)=0.
                 else
                    WDGN(NGN,ith)=Wgt(ith)*Wgt2(ith)
                 end if
              end do
        end if
C -------------------------------------------------------------------	
        irejected3=0
        do ith=1,NTHIN+1
           IF (KEEP3(ith) .EQ. 0) irejected3=irejected3+1
        end do 

c Electron rejected in all thinning levels
        if (irejected3.eq.NTHIN+1) then
           GOTO 66
        else
c Electron kept in at least one thinning level: keep propagating 
          E1=E3
          TYPEP=1.
          XSTART=X
          YSTART=Y
          ZSTART=T1
          do ith=1,NTHIN+1
            if (keep3(ith).eq.0) then    ! Update weights before propagation
              Wgt(ith)=0.
            else
              Wgt(ith)=Wgt(ith)*Wgt3(ith)
            end if
          end do
          nen2=nen2+1
          DdCT=0.D0
          dCTSTA=dCT
          CT_1=CT_2
          DELTA_CT=0.D0
          GO TO 450
        end if

      END IF
C
C -----------------------------------------------------------------
C -----------------------------------------------------------------
C ---------------------- PHOTOELECTRIC EFFECT ---------------------
C -----------------------------------------------------------------
C -----------------------------------------------------------------
C
C  455 E3=E1-EKE+e_MASS
C      IF(E3.LE.CH_THR) GO TO 6
C      E1=E3
C      GO TO 450
C *****************************************************************
C *****************************************************************
C ----------------------------------------------------------------
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

66    IF(NEN.NE.0) THEN
C Extract electrons from the stack 
        IF(NEN.LT.1) STOP 9999
        IF (NEN .GT. MAX_NEN) MAX_NEN = NEN
          E1=DEN(1,NEN)
          T1=DEN(2,NEN)
          dCTSTA=DEN(4,NEN)
          DL=DEN8(1,NEN)
          DM=DEN8(2,NEN)
          DN=DEN8(3,NEN)
          do ith=1,NTHIN+1
             Wgt(ith)=WDEN(NEN,ith)
          end do
          TYPEP=DEN(5,NEN)
          X=DEN(8,NEN)
          Y=DEN(9,NEN)
          CT_1=DEN(6,NEN)
          DELTA_CT=0.0D0
          NEN=NEN-1
          XSTART=X
          YSTART=Y
          ZSTART=T1
          DdCT=0.D0
        IF (TYPEP) 400,1,450
      ELSE IF (NGN.NE.0) THEN 
        GO TO 24
      END IF
C ------------------------------------------------------------------

   53 WRITE(*,56) NGN, NEN
   56 FORMAT(1X,' GAMMA BUFFER COUNTER',I5,
     #          ' ELECTRON POSITRON COUNTER',I5)
      if (ifield_t.eq.1) then 
      NSUB=NIN+NCROSS1+NCROSS
      WRITE(*,*) 'TOTAL NUMBER OF SUBTRACKS',NSUB
      WRITE(*,*) 'SUBTRACKS WITH ZERO LENGTH',NZEROS
      WRITE(*,*) '% OF TRACKS INSIDE BIN',float(NIN)*100./float(NSUB)
      WRITE(*,*) '% OF TRACKS CROSSING 1 BIN',
     #           float(NCROSS1)*100./float(NSUB)
      WRITE(*,*) '% OF TRACKS CROSSING >1 BIN',
     #           float(NCROSS)*100./float(NSUB)
      end if
C
C ***************************************************************
C *************************=============*************************
C ************************* SHOWER ENDS *************************
C *************************=============*************************
C ***************************************************************
      IF (IWRITE_W .EQ. 1) THEN 
         do ith=1,NTHIN+1
            CLOSE(39+ith)
         end do
      END IF

      IF (IWRITE_W_2 .EQ. 1) THEN 
         do ith=1,NTHIN+1
            CLOSE(59+ith)
         end do
      END IF

      IF (IWRITE_T .EQ. 1) THEN 
         do ith=1,NTHIN+1
            CLOSE(69+ith)
         end do
      END IF

      write(*,*) 'Shower Ended'
      write(*,*) 'Energy conservation hybrid shower ',E_CONS/EI
      write(*,*) 'Max number of e- / e+ on stack: ',MAX_NEN
      write(*,*) 'Max number of gammas on stack: ',MAX_NGN

      write(*,*) '#######################################'
      write(*,*) 'Computing fields and writing files.....'
      write(*,*) '#######################################'
c -----------      
c Fraunhofer
c -----------      
C Calculates Electric Field Modulus and Phase
C (one field per thinning level, frequency and observation angle)
      do ith=1,NTHIN+1
        DO NU=1,NUMAX
          DO J=1,JMAX
            DO K=1,KMAX
c 9.6064136D-12=(e mu_r)1.e6/(2 pi epsilon_0 c)
c Factor accompanying the dimensionless part of the electric field.
c The rest of this factor goes into the subroutine EMPSUM in the 
c variable FACFRQ.
            ZSUM_X(J,K,NU,ith)=ZSUM_X(J,K,NU,ith)*9.6064136D-12
            ZSUM_Y(J,K,NU,ith)=ZSUM_Y(J,K,NU,ith)*9.6064136D-12
            ZSUM_Z(J,K,NU,ith)=ZSUM_Z(J,K,NU,ith)*9.6064136D-12
            ZSUM(J,K,NU,ith)=ZSUM(J,K,NU,ith)*9.6064136D-12

            RADIUS(J,K,NU,ith)=SQRT(
     #      ZSUM_X(J,K,NU,ith)*CONJG(ZSUM_X(J,K,NU,ith))+
     #      ZSUM_Y(J,K,NU,ith)*CONJG(ZSUM_Y(J,K,NU,ith))+
     #      ZSUM_Z(J,K,NU,ith)*CONJG(ZSUM_Z(J,K,NU,ith)))

            PHASE(J,K,NU,ith)=
     #      ATAN2(DIMAG(ZSUM(J,K,NU,ith)),DREAL(ZSUM(J,K,NU,ith)))

c x,y,z components
            RADIUS_X(J,K,NU,ith)=SQRT(
     #      ZSUM_X(J,K,NU,ith)*CONJG(ZSUM_X(J,K,NU,ith)))
            PHASE_X(J,K,NU,ith)=
     #      ATAN2(DIMAG(ZSUM_X(J,K,NU,ith)),DREAL(ZSUM_X(J,K,NU,ith)))

            RADIUS_Y(J,K,NU,ith)=SQRT(
     #      ZSUM_Y(J,K,NU,ith)*CONJG(ZSUM_Y(J,K,NU,ith)))
            PHASE_Y(J,K,NU,ith)=
     #      ATAN2(DIMAG(ZSUM_Y(J,K,NU,ith)),DREAL(ZSUM_Y(J,K,NU,ith)))

            RADIUS_Z(J,K,NU,ith)=SQRT(
     #      ZSUM_Z(J,K,NU,ith)*CONJG(ZSUM_Z(J,K,NU,ith)))
            PHASE_Z(J,K,NU,ith)=
     #      ATAN2(DIMAG(ZSUM_Z(J,K,NU,ith)),DREAL(ZSUM_Z(J,K,NU,ith)))


c --------------------------------
c Second calculation in frequency 
c --------------------------------
            ZSUM2_X(J,K,NU,ith)=ZSUM2_X(J,K,NU,ith)*9.6064136D-12
            ZSUM2_Y(J,K,NU,ith)=ZSUM2_Y(J,K,NU,ith)*9.6064136D-12
            ZSUM2_Z(J,K,NU,ith)=ZSUM2_Z(J,K,NU,ith)*9.6064136D-12
            ZSUM2(J,K,NU,ith)=ZSUM2(J,K,NU,ith)*9.6064136D-12

            ARADIUS(J,K,NU,ith)=SQRT(
     #      ZSUM2_X(J,K,NU,ith)*CONJG(ZSUM2_X(J,K,NU,ith))+
     #      ZSUM2_Y(J,K,NU,ith)*CONJG(ZSUM2_Y(J,K,NU,ith))+
     #      ZSUM2_Z(J,K,NU,ith)*CONJG(ZSUM2_Z(J,K,NU,ith)))


c Minus sign so that phase agrees with ZHS PRD 92             
            APHASE(J,K,NU,ith)=
     #   -ATAN2(DIMAG(ZSUM2(J,K,NU,ith)),DREAL(ZSUM2(J,K,NU,ith)))

c x,y,z components
            ARADIUS_X(J,K,NU,ith)=SQRT(
     #      ZSUM2_X(J,K,NU,ith)*CONJG(ZSUM2_X(J,K,NU,ith)))
            APHASE_X(J,K,NU,ith)=
     #     -ATAN2(DIMAG(ZSUM2_X(J,K,NU,ith)),DREAL(ZSUM2_X(J,K,NU,ith)))

            ARADIUS_Y(J,K,NU,ith)=SQRT(
     #      ZSUM2_Y(J,K,NU,ith)*CONJG(ZSUM2_Y(J,K,NU,ith)))
            APHASE_Y(J,K,NU,ith)=
     #     -ATAN2(DIMAG(ZSUM2_Y(J,K,NU,ith)),DREAL(ZSUM2_Y(J,K,NU,ith)))

            ARADIUS_Z(J,K,NU,ith)=SQRT(
     #      ZSUM2_Z(J,K,NU,ith)*CONJG(ZSUM2_Z(J,K,NU,ith)))
            APHASE_Z(J,K,NU,ith)=
     #     -ATAN2(DIMAG(ZSUM2_Z(J,K,NU,ith)),DREAL(ZSUM2_Z(J,K,NU,ith)))


            END DO
          END DO
        END DO

c Time domain
c FAC_T = (e mu_r)/(8 pi epsilon_0 c)        
        FAC_T=2.40160339973116d-18
        DO IT=1,MAXT
          DO J=1,JMAX
            DO K=1,KMAX
            ASUM_X(IT,J,K,ith)=FAC_T*ASUM_X(IT,J,K,ith)
            ASUM_Y(IT,J,K,ith)=FAC_T*ASUM_Y(IT,J,K,ith)
            ASUM_Z(IT,J,K,ith)=FAC_T*ASUM_Z(IT,J,K,ith)
c Minus sign in modulus so that the electric field has
c the correct bipolarity. Note that modulus gives always
c positive sign. No need to put sign by hand in the x,y,z
c components of vector potential.
            ASUM(IT,J,K,ith)=-sqrt(
     #                     ASUM_X(IT,J,K,ith)*ASUM_X(IT,J,K,ith)+
     #                     ASUM_Y(IT,J,K,ith)*ASUM_Y(IT,J,K,ith)+
     #                     ASUM_Z(IT,J,K,ith)*ASUM_Z(IT,J,K,ith))
            END DO
          END DO
        END DO

      end do


c -----------      
c Fresnel
c -----------      
C Calculates Electric Field Modulus and Phase
C (one field per thinning level, frequency and observation angle)
      do ith=1,NTHIN+1
        DO NU=1,NUMAX
          DO NA=1,NAMAX
c 9.6064136D-12=(e mu_r)1.e6/(2 pi epsilon_0 c)
c Factor accompanying the dimensionless part of the electric field.
c The rest of this factor goes into the subroutine EMPSUM in the 
c variable FACFRQ.
c There is also FACTOR_F to convert V/MHz/g cm^-2 to V/MHz/m since 
c the distance to the antenna (R) is in g cm^-2          
            ZSUM2_X_FR(NA,NU,ith)=
     #      ZSUM2_X_FR(NA,NU,ith)*9.6064136D-12
            ZSUM2_Y_FR(NA,NU,ith)=
     #      ZSUM2_Y_FR(NA,NU,ith)*9.6064136D-12
            ZSUM2_Z_FR(NA,NU,ith)=
     #      ZSUM2_Z_FR(NA,NU,ith)*9.6064136D-12
            ZSUM2_ANT_FR(NA,NU,ith)=
     #      ZSUM2_ANT_FR(NA,NU,ith)*9.6064136D-12

            ZSUM2_X_FR(NA,NU,ith)=FACTOR_F*ZSUM2_X_FR(NA,NU,ith)
            ZSUM2_Y_FR(NA,NU,ith)=FACTOR_F*ZSUM2_Y_FR(NA,NU,ith)
            ZSUM2_Z_FR(NA,NU,ith)=FACTOR_F*ZSUM2_Z_FR(NA,NU,ith)
            ZSUM2_ANT_FR(NA,NU,ith)=FACTOR_F*ZSUM2_ANT_FR(NA,NU,ith)

            ARADIUS_FR(NA,NU,ith)=SQRT(
     #      ZSUM2_X_FR(NA,NU,ith)*CONJG(ZSUM2_X_FR(NA,NU,ith))+
     #      ZSUM2_Y_FR(NA,NU,ith)*CONJG(ZSUM2_Y_FR(NA,NU,ith))+
     #      ZSUM2_Z_FR(NA,NU,ith)*CONJG(ZSUM2_Z_FR(NA,NU,ith)))


c Minus sign so that phase agrees with ZHS PRD 92             
            APHASE_FR(NA,NU,ith)=
     #   -ATAN2(DIMAG(ZSUM2_FR(NA,NU,ith)),DREAL(ZSUM2_FR(NA,NU,ith)))

c x,y,z components
            ARADIUS_X_FR(NA,NU,ith)=SQRT(
     #      ZSUM2_X_FR(NA,NU,ith)*CONJG(ZSUM2_X_FR(NA,NU,ith)))
            APHASE_X_FR(NA,NU,ith)=
     #     -ATAN2(
     #      DIMAG(ZSUM2_X_FR(NA,NU,ith)),DREAL(ZSUM2_X_FR(NA,NU,ith)))

            ARADIUS_Y_FR(NA,NU,ith)=SQRT(
     #      ZSUM2_Y_FR(NA,NU,ith)*CONJG(ZSUM2_Y_FR(NA,NU,ith)))
            APHASE_Y_FR(NA,NU,ith)=
     #     -ATAN2(
     #      DIMAG(ZSUM2_Y_FR(NA,NU,ith)),DREAL(ZSUM2_Y_FR(NA,NU,ith)))

            ARADIUS_Z_FR(NA,NU,ith)=SQRT(
     #      ZSUM2_Z_FR(NA,NU,ith)*CONJG(ZSUM2_Z_FR(NA,NU,ith)))
            APHASE_Z_FR(NA,NU,ith)=
     #     -ATAN2(
     #      DIMAG(ZSUM2_Z_FR(NA,NU,ith)),DREAL(ZSUM2_Z_FR(NA,NU,ith)))

c Field at the antenna
            ARADIUS_ANT_FR(NA,NU,ith)=SQRT(
     #      ZSUM2_ANT_FR(NA,NU,ith)*CONJG(ZSUM2_ANT_FR(NA,NU,ith)))
            APHASE_ANT_FR(NA,NU,ith) =
     #     -ATAN2(DIMAG(ZSUM2_ANT_FR(NA,NU,ith)),
     #      DREAL(ZSUM2_ANT_FR(NA,NU,ith)))



          END DO
        END DO

c Time domain
c FAC_T = (e mu_r)/(8 pi epsilon_0 c)        
        FAC_T=2.40160339973116d-18
        DO IT=1,MAXT
          DO NA=1,NAMAX
            ASUM_X_FR(IT,NA,ith)=FAC_T*FACTOR_F*ASUM_X_FR(IT,NA,ith)
            ASUM_Y_FR(IT,NA,ith)=FAC_T*FACTOR_F*ASUM_Y_FR(IT,NA,ith)
            ASUM_Z_FR(IT,NA,ith)=FAC_T*FACTOR_F*ASUM_Z_FR(IT,NA,ith)
            ASUM_ANT_FR(IT,NA,ith)=FAC_T*FACTOR_F*ASUM_ANT_FR(IT,NA,ith)
c Minus sign in modulus so that the electric field has
c the correct bipolarity. Note that modulus gives always
c positive sign. No need to put sign by hand in the x,y,z
c components of vector potential.
            ASUM_FR(IT,NA,ith)=-sqrt(
     #          ASUM_X_FR(IT,NA,ith)*ASUM_X_FR(IT,NA,ith)+
     #          ASUM_Y_FR(IT,NA,ith)*ASUM_Y_FR(IT,NA,ith)+
     #          ASUM_Z_FR(IT,NA,ith)*ASUM_Z_FR(IT,NA,ith))

          END DO
        END DO

      end do

C ---------------------------------------------------------------------
      NO0=NO0+1
      IF(NO0.LE.NOSHOW.AND.NO0.NE.1) GO TO 550

C *********************************************************************
C *********************************************************************
C Write shower output 
C *********************************************************************
C *********************************************************************
c Open files only once (after the 1st shower is done)
      OPEN(21,FILE='DR'//SUFFIX//'.EMP',STATUS='UNKNOWN')
      OPEN(22,FILE='DT'//SUFFIX//'.EMP',STATUS='UNKNOWN')
      OPEN(23,FILE='TK'//SUFFIX//'.EMP',STATUS='UNKNOWN')
      OPEN(28,FILE='MM'//SUFFIX//'.EMP',STATUS='UNKNOWN')
      OPEN(29,FILE='TQ'//SUFFIX//'.EMP',STATUS='UNKNOWN')
      OPEN(30,FILE='ST'//SUFFIX//'.EMP',STATUS='UNKNOWN')
      OPEN(31,FILE='DTNOTHIN'//SUFFIX//'.EMP',STATUS='UNKNOWN')

c Fourier electric field original ZHS - Fraunhofer     
      if (ifield_w.eq.1) then
        OPEN(27,FILE='PF'//SUFFIX//'.EMP',STATUS='UNKNOWN')
        OPEN(25,FILE='ZF'//SUFFIX//'.EMP',STATUS='UNKNOWN')
        OPEN(26,FILE='ZF'//SUFFIX//'_XYZ.EMP',STATUS='UNKNOWN')
      end if

c Fourier electric field explicit tracks - Fraunhofer
      if (ifield_w_2.eq.1) then 
        OPEN(34,FILE='PF'//SUFFIX//'_W_2.EMP',STATUS='UNKNOWN')
        OPEN(35,FILE='ZF'//SUFFIX//'_W_2.EMP',STATUS='UNKNOWN')
        OPEN(36,FILE='ZF'//SUFFIX//'_W_2_XYZ.EMP',STATUS='UNKNOWN')
      end if

c Time domain - Fraunhofer     
      if (ifield_t.eq.1) then 
c Vector potential time domain      
        OPEN(45,FILE='VP'//SUFFIX//'_T_PLOT.EMP',STATUS='UNKNOWN')
        OPEN(46,FILE='VP'//SUFFIX//'_T.EMP',STATUS='UNKNOWN')
        OPEN(48,FILE='VP'//SUFFIX//'_T_XYZ.EMP',STATUS='UNKNOWN')
        OPEN(50,FILE='VP'//SUFFIX//'_T_PERP.EMP',STATUS='UNKNOWN')
c Electric field time domain      
        OPEN(47,FILE='EF'//SUFFIX//'_T.EMP',STATUS='UNKNOWN')
        OPEN(49,FILE='EF'//SUFFIX//'_T_XYZ.EMP',STATUS='UNKNOWN')
        OPEN(51,FILE='EF'//SUFFIX//'_T_PERP.EMP',STATUS='UNKNOWN')
      end if

c Fourier electric field explicit tracks - Fraunhofer
      if (ifield_w_2_fresnel.eq.1) then 
      OPEN(65,FILE='ZF'//SUFFIX//'_W_2_FRESNEL.EMP',STATUS='UNKNOWN')
      OPEN(66,
     #     FILE='ZF'//SUFFIX//'_W_2_XYZ_FRESNEL.EMP',STATUS='UNKNOWN')
      OPEN(67,
     #FILE='ZF'//SUFFIX//'_W_2_ANT_FRESNEL.EMP',STATUS='UNKNOWN')
      OPEN(68,
     #FILE='ZF'//SUFFIX//'_W_2_SIDE_FRESNEL.EMP',STATUS='UNKNOWN')
      end if

c Time domain - Fraunhofer     
      if (ifield_t_fresnel.eq.1) then 
c Vector potential time domain      
      OPEN(76,FILE='VP'//SUFFIX//'_T_FRESNEL.EMP',STATUS='UNKNOWN')
      OPEN(78,FILE='VP'//SUFFIX//'_T_XYZ_FRESNEL.EMP',STATUS='UNKNOWN')
      OPEN(79,FILE='VP'//SUFFIX//'_T_ANT_FRESNEL.EMP',STATUS='UNKNOWN')
c Electric field time domain      
      OPEN(87,FILE='EF'//SUFFIX//'_T_FRESNEL.EMP',STATUS='UNKNOWN')
      OPEN(89,FILE='EF'//SUFFIX//'_T_XYZ_FRESNEL.EMP',STATUS='UNKNOWN')
      OPEN(90,FILE='EF'//SUFFIX//'_T_ANT_FRESNEL.EMP',STATUS='UNKNOWN')
      end if

550   continue
c Converts from g/cm^2 to m.
      CONV=T0/100./RHO

      do ith=1,NTHIN+1

      if (ith.eq.NTHIN+1) then
          write(21,*) ' '
          write(21,'(A9,1X,I5)') '# Shower ',NO0
          write(21,'(A9,1X,1P,1E9.1)') '# Hybrid ',E_HYBRID
          write(22,*) ' '
          write(22,'(A9,1X,I5)') '# Shower ',NO0
          write(22,'(A9,1X,1P,1E9.1)') '# Hybrid ',E_HYBRID
          write(22,48)
          write(23,*) ' '
          write(23,'(A9,1X,I5)') '# Shower ',NO0
          write(23,'(A9,1X,1P,1E9.1)') '# Hybrid ',E_HYBRID
          write(31,*) ' '
          write(31,'(A9,1X,I5)') '# Shower ',NO0
          write(31,'(A9,1X,1P,1E9.1)') '# Hybrid ',E_HYBRID
          goto 333
      end if

      write(21,*) ' '
      write(21,'(A9,1X,I5)') '# Shower ',NO0
      write(21,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
      write(22,*) ' '
      write(22,'(A9,1X,I5)') '# Shower ',NO0
      write(22,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
      write(22,48)
      write(23,*) ' '
      write(23,'(A9,1X,I5)') '# Shower ',NO0
      write(23,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
      write(31,*) ' '
      write(31,'(A9,1X,I5)') '# Shower ',NO0
      write(31,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)

333   JLIM=0
      NRLIM=0
      DO NR=1,NRMAX

        DO J=1,20
          XNELEC(NR,ith)=XNELEC(NR,ith)+INT(XLE(J,NR,ith))
          NELECNOTHIN(NR,ith)=
     # NELECNOTHIN(NR,ith)+INT(XLENOTHIN(J,NR,ith))
          XNPOSI(NR,ith)=XNPOSI(NR,ith)+INT(XLP(J,NR,ith))
          NPOSINOTHIN(NR,ith)=
     # NPOSINOTHIN(NR,ith)+INT(XLPNOTHIN(J,NR,ith))
          IF ((J .GT. JLIM) .AND. 
     #      (XLP(J,NR,ith) .NE. 0 .OR. XLE(J,NR,ith) .NE. 0)) JLIM=J
C         NGAMM=NGAMM+LG(J)
        END DO

c Size Hybrid shower  
        IF(ith.eq.NTHIN+1) THEN  
         XNEPSUM(NR,NTHIN+1)=XNEPSUM(NR,NTHIN+1)+XNEPSUM_HYBRID(NR)
        END IF 

        IF (XNELEC(NR,ith) .GT. 0.1 .OR. XNPOSI(NR,ith) .GT. 0.1 .OR.
     #     XNEPSUM(NR,ith) .GT. 0.1) NRLIM=NR
        XNEXCES(NR,ith)=XNELEC(NR,ith)-XNPOSI(NR,ith)
        NEXCESNOTHIN(NR,ith)=NELECNOTHIN(NR,ith)-NPOSINOTHIN(NR,ith)
        XNEPSUM(NR,ith)=XNELEC(NR,ith)+XNPOSI(NR,ith)
        NEPSUMNOTHIN(NR,ith)=NELECNOTHIN(NR,ith)+NPOSINOTHIN(NR,ith)

c Size Hybrid shower ---- Do we need this line? CHECK please !!! 
        IF(ith.eq.NTHIN+1) THEN  
          XNEPSUM(NR,NTHIN+1)=XNEPSUM(NR,NTHIN+1)+XNEPSUM_HYBRID(NR)
        END IF 
  
      END DO
c Tracklength Hybrid shower  
         TCKABS(1,NTHIN+1)=TCKABS(1,NTHIN+1)+TOT_TCK_HYBRID

      NLOOPS=(NRLIM-1)/10
      WRITE(21,223) 

      DO NLO=0,NLOOPS-1
        WRITE(21,228) (NLO*10+NR,NR=1,10)
        WRITE(21,227) (RMIDDL(J),(INT(XLE(J,NR,ith)+XLP(J,NR,ith)) ,
     #NR=NLO*10+1,NLO*10+10),J=1,JLIM)
      END DO

      WRITE(21,228) (NRAD,NRAD=NLOOPS*10+1,NRLIM)
      DO J=1,JLIM
        WRITE(21,227) RMIDDL(J),
     #(INT(XLE(J,NR,ith)+XLP(J,NR,ith)),NR=NLOOPS*10+1,NRLIM)
      END DO

      WRITE(21,226) 
      DO NLO=0,NLOOPS-1
        WRITE(21,228) (NLO*10+NR,NR=1,10)
        WRITE(21,227) (RMIDDL(J),
     #(INT(XLE(J,NR,ith)-XLP(J,NR,ith)),NR=NLO*10+1,NLO*10+10),J=1,JLIM)
      END DO

      WRITE(21,228) (NRAD,NRAD=NLOOPS*10+1,NRLIM)
      DO J=1,JLIM
        WRITE(21,227) RMIDDL(J),
     #(INT(XLE(J,NR,ith)-XLP(J,NR,ith)),NR=NLOOPS*10+1,NRLIM)
      END DO

  223 FORMAT(2X,'\n# ELECTRONS + POSITRONS (LATER.DISTRIB.)')
  226 FORMAT(2X,'\n# ELECTRON EXCESS (LATER.DISTRIB.)')
  227 FORMAT(1X,1F10.5,10I18)
  228 FORMAT(1X,'\n# r/r0, Rad lengths ->',I4,9I7)

      WRITE(22,49) (NR,XNEPSUM(NR,ith),XNELEC(NR,ith),XNEXCES(NR,ith),
     #              XNEXCES(NR,ith)/MAX( 0.5,XNEPSUM(NR,ith) ), 
     #              NRE(NR,ith),NRPOS(NR,ith),NR=1,NRLIM)

      WRITE(31,51) (NR,NEPSUMNOTHIN(NR,ith),NELECNOTHIN(NR,ith),
     # NEXCESNOTHIN(NR,ith),
     # REAL(NEXCESNOTHIN(NR,ith))/MAX( 0.5,REAL(NEPSUMNOTHIN(NR,ith)) ),
     # NRENOTHIN(NR,ith),NRPOSNOTHIN(NR,ith),NR=1,NRLIM)

   48 FORMAT('# (e-)+(e+)',5X,'(e-)',8X,'e-p',4X,' excess',
     #       4X,'NRE(e-)',5X,'NRPOS(e+)',3X,'Backscattered')
   49 FORMAT((1X,I5,1P,4E15.5,2I18))
   51 FORMAT((1X,I5,I18,I18,I18,F9.5,I18,I18))

      WRITE(23,718)
      WRITE(23,719) (ETHRKE(J),ETHRKE(J)+e_MASS,
     #   TCKABS(J,ith)*CONV,TCKPRO(J,ith)*CONV,
     #   -TCKSUMNOPRO(J,ith)*CONV,-TCKSUMPRO(J,ith)*CONV,J=1,10)
     
      end do

  718 FORMAT('# KINETIC & TOTAL ENERGY (MeV)',
     # ' TRACK-SUMS: TOTAL, TOTAL PROJ, EXCESS & EXCESS PROJ (meters)')
  719 FORMAT(1X,1P,2E12.4,1P,4E15.8)

c ######################################################################## 
c Write "moments of inertia" of shower
c ######################################################################## 
      do ith=1,NTHIN+1 

      if (ith.eq.NTHIN+1) then
          write(28,*) ' '
          write(28,'(A9,1X,I5)') '# Shower ',NO0
          write(28,'(A9,1X,1P,1E9.1)') '# Hybrid ',E_HYBRID
          goto 334
      end if

        write(28,*) ' '
        write(28,'(A9,1X,I5)') '# Shower ',NO0
        write(28,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)

334     write(28,*) 
     #'# (e-,e+) Sum(w)'
        do NR=1,NRLIM
        write(28,'(1X,I5,1P,2E14.5)') NR,WSUME(NR,ith),WSUMP(NR,ith)
        end do

        write(28,*) 
     #'# (e-,e+) Sum(w*x)/Sum(w) Sum(w*y)/Sum(w) Sum(w*r)/Sum(w)'
        do NR=1,NRLIM
        if (WSUME(NR,ith).GT.0..and.WSUMP(NR,ith).GT.0.)
     #  write(28,'(1X,I5,1P,6E14.5)') 
     # NR,WXE(NR,ith)/WSUME(NR,ith),WYE(NR,ith)/WSUME(NR,ith),
     # WRE(NR,ith)/WSUME(NR,ith),WXP(NR,ith)/WSUMP(NR,ith),
     # WYP(NR,ith)/WSUMP(NR,ith),WRP(NR,ith)/WSUMP(NR,ith)
        end do

        write(28,*) 
     #'# (e-,e+) Sum(w*x^2)/Sum(w) Sum(w*y^2)/Sum(w) Sum(w*r^2)/Sum(w)'
        do NR=1,NRLIM
        if (WSUME(NR,ith).GT.0..and.WSUMP(NR,ith).GT.0.)
     #  write(28,'(1X,I5,1P,6E14.5)') 
     # NR,WX2E(NR,ith)/WSUME(NR,ith),WY2E(NR,ith)/WSUME(NR,ith),
     # WR2E(NR,ith)/WSUME(NR,ith),WX2P(NR,ith)/WSUMP(NR,ith),
     # WY2P(NR,ith)/WSUMP(NR,ith),WR2P(NR,ith)/WSUMP(NR,ith)
        end do

        write(28,*) 
     #'# (e-,e+) Sum(w*z)/Sum(w) Sum(w*z^2)/Sum(w)'
        do NR=1,NRLIM
        if (WSUME(NR,ith).GT.0..and.WSUMP(NR,ith).GT.0.)
     #  write(28,'(1X,I5,1P,4E14.5)') 
     # NR,WZE(NR,ith)/WSUME(NR,ith),WZ2E(NR,ith)/WSUME(NR,ith),
     # WZP(NR,ith)/WSUMP(NR,ith),WZ2P(NR,ith)/WSUMP(NR,ith)
        end do

      end do

c Write "moments of inertia" of shower
       write(29,*) ' ' 
       write(29,'(A9,1X,I5)') '# Shower ',NO0
       write(29,*) '# Moments of inertia whole shower'
       write(29,*) 
     #'# (e-p) Thin_min,Thin_max,Calls to EMPSUM,Excess track,
     # Sum(w*x)/Sum(w),Sum(w*y)/Sum(w),Sum(w*r)/Sum(w),
     # Sum(w*x^2)/Sum(w),Sum(w*y^2)/Sum(w),Sum(w*r^2)/Sum(w),
     # Sum(w*z)/Sum(w),Sum(w*z^2)/Sum(w)'

       write(30,*) ' ' 
       write(30,'(A9,1X,I5)') '# Shower ',NO0
       write(30,*) '# Moments of inertia at shower maximum'
       write(30,*) 
     #'# (e-p) Thin_min,Thin_max,Calls to EMPSUM,(Xmax/X0),
     #Sum(w*x)/Sum(w),Sum(w*y)/Sum(w),Sum(w*r)/Sum(w),
     #Sum(w*x^2)/Sum(w),Sum(w*y^2)/Sum(w),Sum(w*r^2)/Sum(w),
     #Sum(w*z)/Sum(w),Sum(w*z^2)/Sum(w)'

      do ith=1,NTHIN+1
         WXETOT=0.
         WYETOT=0.
         WRETOT=0.
         WX2ETOT=0.
         WY2ETOT=0.
         WR2ETOT=0.
         WZETOT=0.
         WZ2ETOT=0.
         WETOT=0.

         WXPTOT=0.
         WYPTOT=0.
         WRPTOT=0.
         WX2PTOT=0.
         WY2PTOT=0.
         WR2PTOT=0.
         WZPTOT=0.
         WZ2PTOT=0.
         WPTOT=0.


         WEPMAX=0. 
         do nr=1,nrlim

             if ((WSUME(NR,ith)-WSUMP(NR,ith)).gt.WEPMAX) then 
                 NRLXMAX=NR    ! Find bin of shower max.
                 WEPMAX=WSUME(NR,ith)-WSUMP(NR,ith)
             end if

             WETOT=WETOT+WSUME(NR,ith)

             WXETOT=WXETOT+WXE(NR,ith)
             WYETOT=WYETOT+WYE(NR,ith)
             WRETOT=WRETOT+WRE(NR,ith)

             WX2ETOT=WX2ETOT+WX2E(NR,ith)
             WY2ETOT=WY2ETOT+WY2E(NR,ith)
             WR2ETOT=WR2ETOT+WR2E(NR,ith)

             WZETOT=WZETOT+WZE(NR,ith)
             WZ2ETOT=WZ2ETOT+WZ2E(NR,ith)
c ---------------------------------------------------------
             WPTOT=WPTOT+WSUMP(NR,ith)

             WXPTOT=WXPTOT+WXP(NR,ith)
             WYPTOT=WYPTOT+WYP(NR,ith)
             WRPTOT=WRPTOT+WRP(NR,ith)

             WX2PTOT=WX2PTOT+WX2P(NR,ith)
             WY2PTOT=WY2PTOT+WY2P(NR,ith)
             WR2PTOT=WR2PTOT+WR2P(NR,ith)

             WZPTOT=WZPTOT+WZP(NR,ith)
             WZ2PTOT=WZ2PTOT+WZ2P(NR,ith)

         end do

         WEXCESS=WETOT-WPTOT

         if (WEXCESS.GT.0) THEN
           write(29,'(1X,1P,2E9.1,I18,1P,9E14.5)') 
     #THIN_MIN(ith),THIN_MAX(ith),
     #NCALLS_EMPSUM_W(ith),
     #-TCKSUMNOPRO(1,ith)*CONV,
     #(WXETOT-WXPTOT)/WEXCESS,
     #(WYETOT-WYPTOT)/WEXCESS,
     #(WRETOT-WRPTOT)/WEXCESS,
     #(WX2ETOT-WX2PTOT)/WEXCESS,
     #(WY2ETOT-WY2PTOT)/WEXCESS,
     #(WR2ETOT-WR2PTOT)/WEXCESS,
     #(WZETOT-WZPTOT)/WEXCESS,
     #(WZ2ETOT-WZ2PTOT)/WEXCESS

      WEXCESS=WSUME(NRLXMAX,ith)-WSUMP(NRLXMAX,ith)
          
           write(30,'(1X,1P,2E9.1,I7,I20,1P,8E14.5)') 
     #THIN_MIN(ith),THIN_MAX(ith),NRLXMAX,NCALLS_EMPSUM_W(ith),
     #(WXE(NRLXMAX,ith)-WXP(NRLXMAX,ith))/WEXCESS,
     #(WYE(NRLXMAX,ith)-WYP(NRLXMAX,ith))/WEXCESS,
     #(WRE(NRLXMAX,ith)-WRP(NRLXMAX,ith))/WEXCESS,
     #(WX2E(NRLXMAX,ith)-WX2P(NRLXMAX,ith))/WEXCESS,
     #(WY2E(NRLXMAX,ith)-WY2P(NRLXMAX,ith))/WEXCESS,
     #(WR2E(NRLXMAX,ith)-WR2P(NRLXMAX,ith))/WEXCESS,
     #(WZE(NRLXMAX,ith)-WZP(NRLXMAX,ith))/WEXCESS,
     #(WZ2E(NRLXMAX,ith)-WZ2P(NRLXMAX,ith))/WEXCESS

         END IF
         
      end do 
 
c ######################################################################## 
c Fourier electric field original ZHS - Fraunhofer     
c ######################################################################## 
      if (ifield_w.eq.1) then 
        DO ITH=1,NTHIN
          write(27,'(A9,1X,I5)') '# Shower ',NO0
          write(27,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO NU=1,NUMAX
            WRITE(27,'(A42,1F10.3)') 
     # '# Angular distribution at Frequency [MHz] ',FREQ(NU)
            WRITE(27,*) 
     # '# at azimuth angles [deg]',(PHI(K)*RADEGR,K=1,KMAX)
            DO J=1,JMAX  
              WRITE(27,903) THETA(J)*RADEGR,
     # (RADIUS(J,K,NU,ITH),PHASE(J,K,NU,ITH),K=1,KMAX)
            END DO
            WRITE(27,*) ''
          END DO
        END DO
      
c ------------------------------------------------------------------------ 
        DO ITH=1,NTHIN
          write(25,'(A9,1X,I5)') '# Shower ',NO0
          write(25,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO J=1,JMAX
            WRITE(25,'(A36,1F14.6)') 
     #'# Frequency Spectrum at Angle [deg] ',THETA(J)*RADEGR
            WRITE(25,*) 
     # '# at azimuth angles [deg]',(PHI(K)*RADEGR,K=1,KMAX)

            DO NU=1,NUMAX
              WRITE(25,903) FREQ(NU),
     #(RADIUS(J,K,NU,ITH),PHASE(J,K,NU,ITH),K=1,KMAX)
            END DO
            WRITE(25,*) ''
          END DO
        END DO

c ------------------------------------------------------------------------ 
        DO ITH=1,NTHIN
          write(26,'(A9,1X,I5)') '# Shower ',NO0
          write(26,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO J=1,JMAX
            WRITE(26,'(A55,1F14.6)') 
     #'# Frequency Spectrum (x,y,z) components at Angle [deg] ',
     # THETA(J)*RADEGR
            WRITE(26,'(A25,5(1F14.6))') 
     # '# at azimuth angles [deg]',(PHI(K)*RADEGR,K=1,KMAX)

            DO NU=1,NUMAX
              WRITE(26,933) FREQ(NU),
     #(RADIUS_X(J,K,NU,ITH),PHASE_X(J,K,NU,ITH),
     # RADIUS_Y(J,K,NU,ITH),PHASE_Y(J,K,NU,ITH),
     # RADIUS_Z(J,K,NU,ITH),PHASE_Z(J,K,NU,ITH),
     # K=1,KMAX)
            END DO
            WRITE(26,*) ''
          END DO
        END DO


      end if
  903 FORMAT(1X,1F11.3,1P,5(1E18.10,0P,1F14.8))
  905 FORMAT(1X,1F11.3,1P,5(1E18.10))
  933 FORMAT(1X,1F11.3,1P,5(6E18.10))
      
c ######################################################################## 
c Fourier electric field explicit tracks - Fraunhofer     
c ######################################################################## 
      if (ifield_w_2.eq.1) then 
        DO ITH=1,NTHIN
          write(34,'(A9,1X,I5)') '# Shower ',NO0
          write(34,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO NU=1,NUMAX
            WRITE(34,'(A45,1F10.3)') 
     # '# Angular distribution at Frequency [MHz] ',FREQ(NU)
            WRITE(34,*) 
     # '# at azimuth angles [deg]',(PHI(K)*RADEGR,K=1,KMAX)

            DO J=1,JMAX  
              WRITE(34,903) THETA(J)*RADEGR,
     #(ARADIUS(J,K,NU,ITH),APHASE(J,K,NU,ITH),K=1,KMAX)
            END DO
            WRITE(34,*) ''
          END DO
        END DO
      
c ------------------------------------------------------------------------ 
        DO ITH=1,NTHIN
          write(35,'(A9,1X,I5)') '# Shower ',NO0
          write(35,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO J=1,JMAX
            WRITE(35,'(A38,1F14.6)') 
     #'# Frequency Spectrum at Angle [deg] ',THETA(J)*RADEGR
            WRITE(35,*) 
     # '# at azimuth angles [deg]',(PHI(K)*RADEGR,K=1,KMAX)
            DO NU=1,NUMAX
              WRITE(35,903) FREQ(NU),
     #(ARADIUS(J,K,NU,ITH),APHASE(J,K,NU,ITH),K=1,KMAX)
            END DO
            WRITE(35,*) ''
          END DO
        END DO

c ------------------------------------------------------------------------ 
        DO ITH=1,NTHIN
          write(36,'(A9,1X,I5)') '# Shower ',NO0
          write(36,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO J=1,JMAX
            WRITE(36,'(A55,1F14.6)') 
     #'# Frequency Spectrum (x,y,z) components at Angle [deg] ',
     # THETA(J)*RADEGR




            WRITE(36,*) 
     # '# at azimuth angles [deg]',(PHI(K)*RADEGR,K=1,KMAX)

            DO NU=1,NUMAX
              WRITE(36,933) FREQ(NU),
     #(ARADIUS_X(J,K,NU,ITH),APHASE_X(J,K,NU,ITH),
     # ARADIUS_Y(J,K,NU,ITH),APHASE_Y(J,K,NU,ITH),
     # ARADIUS_Z(J,K,NU,ITH),APHASE_Z(J,K,NU,ITH),
     # K=1,KMAX)
            END DO
            WRITE(36,*) ''
          END DO
        END DO


      end if

c ######################################################################## 
c Vector potential time domain - Fraunhofer     
c ######################################################################## 
      if (ifield_t.eq.1) then 

        DO ITH=1,NTHIN
          write(45,'(A9,1X,I5)') '# Shower ',NO0
          write(46,'(A9,1X,I5)') '# Shower ',NO0
          write(48,'(A9,1X,I5)') '# Shower ',NO0
          write(50,'(A9,1X,I5)') '# Shower ',NO0
          write(45,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          write(46,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          write(48,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          write(50,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)

          DO J=1,JMAX
            WRITE(45,'(A46,1F14.6)') 
     #'# Time domain vector potential at Angle [deg] ',THETA(J)*RADEGR
            WRITE(45,*)
     #'# at azimuth angles [deg] ',(PHI(K)*RADEGR,K=1,KMAX)
            WRITE(46,'(A46,1F14.6)') 
     #'# Time domain vector potential at Angle [deg] ',THETA(J)*RADEGR
            WRITE(46,*)
     #'# at azimuth angles [deg] ',(PHI(K)*RADEGR,K=1,KMAX)
            WRITE(48,'(A40,1F14.6)') 
     #'# Vector potential x,y,z at Angle [deg] ',THETA(J)*RADEGR
            WRITE(48,*)
     #'# at azimuth angles [deg] ',(PHI(K)*RADEGR,K=1,KMAX)
            WRITE(50,'(A40,1F14.6)') 
     #'# Vector potential perp. at Angle [deg] ',THETA(J)*RADEGR
            WRITE(50,*)
     #'# at azimuth angles [deg] ',(PHI(K)*RADEGR,K=1,KMAX)

            DO IT=1,MAXT
              TIME=(IT-MAXTHALF)*DT_BIN_NS+0.5*DT_BIN_NS
              IF (TIME.GT.TIME_MAX.OR.TIME.LT.TIME_MIN) GOTO 700 
                WRITE(46,904) TIME,(ASUM(IT,J,K,ITH),K=1,KMAX)
                WRITE(48,944) TIME,
     #          (ASUM_X(IT,J,K,ITH),
     #           ASUM_Y(IT,J,K,ITH),
     #           ASUM_Z(IT,J,K,ITH),K=1,KMAX)
c k x (k x asum) 
                DO K=1,KMAX
                  UX = SINTET(J)*COSPHI(K)
                  UY = SINTET(J)*SINPHI(K)
                  UZ = COSTET(J)

        ASUM_PERP_X = -(UY*UY+UZ*UZ)*ASUM_X(IT,J,K,ITH) 
     #                + UX*UY*ASUM_Y(IT,J,K,ITH)  
     #                + UX*UZ*ASUM_Z(IT,J,K,ITH)
        ASUM_PERP_Y =    UX*UY*ASUM_X(IT,J,K,ITH) 
     #                - (UX*UX+UZ*UZ)*ASUM_Y(IT,J,K,ITH) 
     #                + UY*UZ*ASUM_Z(IT,J,K,ITH)
        ASUM_PERP_Z = UX*UZ*ASUM_X(IT,J,K,ITH) 
     #              + UY*UZ*ASUM_Y(IT,J,K,ITH) 
     #              - (UX*UX+UY*UY)*ASUM_Z(IT,J,K,ITH)

        ASUM_PERP(IT,J,K,ITH)=-SQRT(ASUM_PERP_X*ASUM_PERP_X +
     #                              ASUM_PERP_Y*ASUM_PERP_Y +
     #                              ASUM_PERP_Z*ASUM_PERP_Z)

                END DO
                WRITE(50,904) TIME,(ASUM_PERP(IT,J,K,ITH),K=1,KMAX)
c Output for plotting purposes with gnuplot        
                WRITE(45,904) (IT-MAXTHALF)*DT_BIN_NS,
     #(ASUM(IT,J,K,ITH),K=1,KMAX)
                WRITE(45,904) (IT+1-MAXTHALF)*DT_BIN_NS,
     #(ASUM(IT,J,K,ITH),K=1,KMAX)
700           CONTINUE  
            END DO

            WRITE(45,*) ''
            WRITE(46,*) ''
            WRITE(48,*) ''
            WRITE(50,*) ''
          END DO

        END DO

c ######################################################################## 
c Electric field time domain - Fraunhofer     
c Minus time derivative of vector potential to obtain electric field
c ######################################################################## 
        DO ITH=1,NTHIN
          write(47,'(A9,1X,I5)') '# Shower ',NO0
          write(49,'(A9,1X,I5)') '# Shower ',NO0
          write(51,'(A9,1X,I5)') '# Shower ',NO0
          write(47,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          write(49,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          write(51,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)

          DO J=1,JMAX
            WRITE(47,'(A44,1F14.6)') 
     #'# Time domain electric field at Angle [deg] ',THETA(J)*RADEGR
            WRITE(47,*) 
     #'# at azimuth angles [deg] ',(PHI(K)*RADEGR,K=1,KMAX)
            WRITE(49,'(A38,1F14.6)') 
     #'# Electric field x,y,z at Angle [deg] ',THETA(J)*RADEGR
            WRITE(49,*) 
     #'# at azimuth angles [deg] ',(PHI(K)*RADEGR,K=1,KMAX)
            WRITE(51,'(A38,1F14.6)') 
     #'# Electric field perp. at Angle [deg] ',THETA(J)*RADEGR
            WRITE(51,*) 
     #'# at azimuth angles [deg] ',(PHI(K)*RADEGR,K=1,KMAX)

           DO IT=1,MAXT-1
             DO K=1,KMAX
             EF_X(K)=-(ASUM_X(IT+1,J,K,ith)-ASUM_X(IT,J,K,ith))/DT_BIN_S
             EF_Y(K)=-(ASUM_Y(IT+1,J,K,ith)-ASUM_Y(IT,J,K,ith))/DT_BIN_S
             EF_Z(K)=-(ASUM_Z(IT+1,J,K,ith)-ASUM_Z(IT,J,K,ith))/DT_BIN_S
             EF_PERP(K)=
     #       -(ASUM_PERP(IT+1,J,K,ith)-ASUM_PERP(IT,J,K,ith))/DT_BIN_S
             EF(K)=-(ASUM(IT+1,J,K,ith)-ASUM(IT,J,K,ith))/DT_BIN_S
             END DO
             TIME=(IT-MAXTHALF)*DT_BIN_NS+0.5*DT_BIN_NS
             IF (TIME.GT.TIME_MAX.OR.TIME.LT.TIME_MIN) GOTO 600 
               WRITE(47,904) TIME,(EF(K),K=1,KMAX)
               WRITE(49,944) TIME,(EF_X(K),EF_Y(K),EF_Z(K),K=1,KMAX)
               WRITE(51,904) TIME,(EF_PERP(K),K=1,KMAX)
600          CONTINUE
           END DO

            WRITE(47,*) ''
            WRITE(49,*) ''
            WRITE(51,*) ''
          END DO

         END DO

      end if


904   FORMAT (1X,1F15.5,1P,5(1E18.10))
944   FORMAT (1X,1F15.5,1P,5(3E18.10))

c ######################################################################## 
c Fourier electric field explicit tracks - Fresnel     
c ######################################################################## 
      if (ifield_w_2_fresnel.eq.1) then 
c ------------------------------------------------------------------------ 
        DO ITH=1,NTHIN
          write(65,'(A9,1X,I5)') '# Shower ',NO0
          write(65,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO NA=1,NAMAX
            WRITE(65,'(A58,3F14.6)') 
     #'# Frequency Spectrum at antenna position (x,y,z) [g cm^-2]',
     #XANT(NA),YANT(NA),ZANT(NA)
            DO NU=1,NUMAX
              WRITE(65,903) FREQ(NU),
     #ARADIUS_FR(NA,NU,ITH),APHASE_FR(NA,NU,ITH)
            END DO
            WRITE(65,*) ''
          END DO
        END DO

c ------------------------------------------------------------------------ 
        DO ITH=1,NTHIN
          write(67,'(A9,1X,I5)') '# Shower ',NO0
          write(67,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO NA=1,NAMAX
            WRITE(67,'(A58,3F14.6)') 
     #'# Frequency Spectrum at antenna position (x,y,z) [g cm^-2]',
     #XANT(NA),YANT(NA),ZANT(NA)
            DO NU=1,NUMAX
              WRITE(67,905) FREQ(NU),ARADIUS_ANT_FR(NA,NU,ITH),
     # APHASE_ANT_FR(NA,NU,ITH)
            END DO
            WRITE(67,*) ''
          END DO
        END DO

c ------------------------------------------------------------------------ 
        DO ITH=1,NTHIN
          DO NU=1,NUMAX
            DO NA=1,NAMAX
              WRITE(68,'(4F13.6,1P,3E18.10)') 
     #        FREQ(NU),XANT(NA),YANT(NA),ZANT(NA),
     #        ARADIUS_X_FR(NA,NU,ITH),
     #        ARADIUS_Y_FR(NA,NU,ITH),
     #        ARADIUS_Z_FR(NA,NU,ITH)
            END DO
            WRITE(68,*) ''
          END DO
        END DO


c ------------------------------------------------------------------------ 
        DO ITH=1,NTHIN
          write(66,'(A9,1X,I5)') '# Shower ',NO0
          write(66,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          DO NA=1,NAMAX
            WRITE(66,'(A60,3F14.6)') 
     #'# Frequency Spectrum (x,y,z) components at antenna position ',
     #XANT(NA),YANT(NA),ZANT(NA)  

            DO NU=1,NUMAX
              WRITE(66,933) FREQ(NU),
     # ARADIUS_X_FR(NA,NU,ITH),APHASE_X_FR(NA,NU,ITH),
     # ARADIUS_Y_FR(NA,NU,ITH),APHASE_Y_FR(NA,NU,ITH),
     # ARADIUS_Z_FR(NA,NU,ITH),APHASE_Z_FR(NA,NU,ITH)
            END DO
            WRITE(66,*) ''
          END DO
        END DO


      end if

906    FORMAT (1X,1F20.5,1P,1(1E18.10),0P,1F20.5)
946    FORMAT (1X,1F15.5,1P,5(3E18.10),0P,1F20.5)
c ######################################################################## 
c Vector potential time domain - Fresnel      
c ######################################################################## 
      if (ifield_t_fresnel.eq.1) then 

c --------------------------------------------------------------------
c Times are referred to a signal travelling at c/n from the injection point
c of the shower to the observer. Compute that time and include in output file
c to refer to time w.r.t. injection point (t=0).
       DO NA=1,NAMAX
c Distance from injection point to antenna position
           RANT=SQRT(XANT(NA)*XANT(NA)+
     #               YANT(NA)*YANT(NA)+
     #               ZANT(NA)*ZANT(NA))

c Time between injection point (at t=0 in ZHS) and arrival of light at antenna position
           T0_ANT(NA)=REFIDX*RANT*FACTOR_T   ! seconds
           T0_ANT(NA)=T0_ANT(NA)*1.e9   ! nano seconds
        END DO


        DO ITH=1,NTHIN
          write(76,'(A9,1X,I5)') '# Shower ',NO0
          write(78,'(A9,1X,I5)') '# Shower ',NO0
          write(79,'(A9,1X,I5)') '# Shower ',NO0
          write(76,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          write(78,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          write(79,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)

          DO NA=1,NAMAX
            WRITE(76,'(A58,3F14.6)') 
     #'# Time domain vector potential at antenna position (x,y,z)',
     # XANT(NA),YANT(NA),ZANT(NA)
            WRITE(78,'(A50,3F14.6)') 
     #'# Vector potential (x,y,z) at antenna position ',
     # XANT(NA),YANT(NA),ZANT(NA)
            WRITE(78,'(A50,3F14.6)') 
     #'# Vector potential (x,y,z) at antenna position ',
     # XANT(NA),YANT(NA),ZANT(NA)

            DO IT=1,MAXT
              TIME=(IT-MAXTHALF)*DT_BIN_NS+0.5*DT_BIN_NS
              IF (TIME.GT.TIME_MAX.OR.TIME.LT.TIME_MIN) GOTO 500 
                WRITE(76,906) TIME,ASUM_FR(IT,NA,ITH),T0_ANT(NA)
                WRITE(78,946) TIME,
     #           ASUM_X_FR(IT,NA,ITH),
     #           ASUM_Y_FR(IT,NA,ITH),
     #           ASUM_Z_FR(IT,NA,ITH),T0_ANT(NA)
                WRITE(79,904) TIME,ASUM_ANT_FR(IT,NA,ITH)

500           CONTINUE  
            END DO

            WRITE(76,*) ''
            WRITE(78,*) ''
            WRITE(79,*) ''
          END DO

        END DO

c ######################################################################## 
c Electric field time domain - Fresnel     
c Minus time derivative of vector potential to obtain electric field
c ######################################################################## 
        DO ITH=1,NTHIN
          write(87,'(A9,1X,I5)') '# Shower ',NO0
          write(89,'(A9,1X,I5)') '# Shower ',NO0
          write(87,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)
          write(89,'(A10,1X,1P,2E9.1)') 
     #'# Thinning',THIN_MAX(ith),THIN_MIN(ith)

          DO NA=1,NAMAX
            WRITE(87,'(A56,3F14.6)') 
     #'# Time domain electric field at antenna position (x,y,z)',
     # XANT(NA),YANT(NA),ZANT(NA)
            WRITE(89,'(A56,3F14.6)') 
     #'# Time domain electric field at antenna position (x,y,z)',
     # XANT(NA),YANT(NA),ZANT(NA)

           DO IT=1,MAXT-1
             EF_X_FR(NA)=
     #-(ASUM_X_FR(IT+1,NA,ith)-ASUM_X_FR(IT,NA,ith))/DT_BIN_S
             EF_Y_FR(NA)=
     #-(ASUM_Y_FR(IT+1,NA,ith)-ASUM_Y_FR(IT,NA,ith))/DT_BIN_S
             EF_Z_FR(NA)=
     #-(ASUM_Z_FR(IT+1,NA,ith)-ASUM_Z_FR(IT,NA,ith))/DT_BIN_S
             EF_FR(NA)=
     #-(ASUM_FR(IT+1,NA,ith)-ASUM_FR(IT,NA,ith))/DT_BIN_S
             EF_ANT_FR(NA)=
     #-(ASUM_ANT_FR(IT+1,NA,ith)-ASUM_ANT_FR(IT,NA,ith))/DT_BIN_S
             TIME=(IT-MAXTHALF)*DT_BIN_NS+0.5*DT_BIN_NS
             IF (TIME.GT.TIME_MAX.OR.TIME.LT.TIME_MIN) GOTO 800 
               WRITE(87,906) TIME,EF_FR(NA),T0_ANT(NA)
               WRITE(89,946) TIME,
     #                       EF_X_FR(NA),EF_Y_FR(NA),EF_Z_FR(NA),
     #                       T0_ANT(NA)
               WRITE(90,904) TIME,EF_ANT_FR(NA)
800          CONTINUE
           END DO

            WRITE(87,*) ''
            WRITE(89,*) ''
          END DO

        END DO

       end if

c ######################################################################## 
c WRITING OF SHOWER NO0 OUTPUT ENDS
c ######################################################################## 
      write(*,*) 'Shower ',NO0,' written'

      IF(NO0.LT.NOSHOW) GO TO 50

c Close all files after the last shower is simulated
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(25)
      CLOSE(26)
      CLOSE(28)
      CLOSE(29)
      CLOSE(30)
      CLOSE(31)
      CLOSE(27)
      CLOSE(25)
      CLOSE(34)
      CLOSE(35)
      CLOSE(36)
      CLOSE(46)
      CLOSE(45)
      CLOSE(47)
      CLOSE(48)
      CLOSE(49)
      CLOSE(50)
      CLOSE(51)
      CLOSE(65)
      CLOSE(66)
      CLOSE(76)
      CLOSE(78)
      CLOSE(87)
      CLOSE(89)

c CPU time
      call CPU_Time(sectot)
      hours=int(sectot/3600.0)
      dummy=sectot-hours*3600.
      minutes=int(dummy/60)
      seconds=dummy-minutes*60.

!      call HostNm(host,status)
!      if (status.eq.0) then
!        write(*,*) 'Host: ',host
!      endif

       write(*,*) 'Cpu time for this run : ',hours,
     #  ' h: ',minutes,' m: ',seconds,' s'

C *********************************************************************
C Write output ends
C *********************************************************************
      CONTINUE
C ----------------------------------------------------------------
      STOP
 9501 PRINT*,'*************** BUFFER OVERFLOW ********************'
      PRINT*,'*************** BUFFER OVERFLOW ********************'
      GO TO 53
      END
C -------------------------------------------------------------------



C *******************************************************************
      SUBROUTINE DISTRS
     #       (ITYPE,TOLD,T1,XOLD,YOLD,X,Y,E4,E1,NRAD,NROLD,Wgt)
C *******************************************************************
C--------------------------------------------------------------------
C
C  THIS SUBR. PICKS UP INFORMATION ONLY FOR
C  ONE LEVEL - ZOBS
C 
C -------------------------------------------------------------------
      PARAMETER (MAXNR=10000)           ! Max. number of rad. lengths
      PARAMETER (NTHIN_MAX=2)              ! Number of simultaneous thinning levels 
      DOUBLE PRECISION Wgt(NTHIN_MAX) 
      COMMON/NLE/XLE(20,MAXNR,NTHIN_MAX),
     #           XLP(20,MAXNR,NTHIN_MAX)/NLG/LG(20)
      COMMON/NLENOTHIN/XLENOTHIN(20,MAXNR,NTHIN_MAX),
     #XLPNOTHIN(20,MAXNR,NTHIN_MAX)
      COMMON/TTRE/T97(20)
      COMMON/TTR/TRR(20)
      common / ZOB / ZOBS
      COMMON / THIN_LEVELS / NTHIN

      COMMON/NVG/NRG
      COMMON/WR0/RFACT
      COMMON / MEDIUM / E_LPM, E_LPM_BREMSS, REF_N,
     #R0,Z,T0,RHO,BEE,C,A,RM,X0,X1,AI,AKO,EKE,AZ,ALZ2,ALZ,
     #Xi_Z,FACNOR,RADLEN, ALZ3, CH_THR,Z_Wgt, FACTOR_MS
      COMMON / MOMENTS / 
     #WXE(MAXNR,NTHIN_MAX),WYE(MAXNR,NTHIN_MAX),
     #WRE(MAXNR,NTHIN_MAX),
     #WX2E(MAXNR,NTHIN_MAX),WY2E(MAXNR,NTHIN_MAX),
     #WR2E(MAXNR,NTHIN_MAX),
     #WXP(MAXNR,NTHIN_MAX),WYP(MAXNR,NTHIN_MAX),
     #WRP(MAXNR,NTHIN_MAX),
     #WX2P(MAXNR,NTHIN_MAX),WY2P(MAXNR,NTHIN_MAX),
     #WR2P(MAXNR,NTHIN_MAX),
     #WZE(MAXNR,NTHIN_MAX),WZ2E(MAXNR,NTHIN_MAX),
     #WSUME(MAXNR,NTHIN_MAX), 
     #WZP(MAXNR,NTHIN_MAX),WZ2P(MAXNR,NTHIN_MAX),
     #WSUMP(MAXNR,NTHIN_MAX) 
C -------------------------------------------------------------------
C      IF ((TOLD. LE. T1) .AND. (T1 .LE. ZOBS)) RETURN

      DO NR=NROLD+1,NRAD
      TINT=T1-TOLD
      T=T0*NR-TOLD
      Tpos=T0*NR
      IF (ITYPE.EQ.0) THEN
        E2=E1
      ELSE
        E2=E1-T*(E1-E4)/TINT
      END IF

C     K=MIN(20.,1.4428*ALOG(E2/CH_THR)+1)

      IF (E2 .LE. CH_THR) THEN
C           IF(K.LT.1) RETURN
            IF (E2 .LT. 0.) WRITE(*,*) 'E2 < 0 '
            RETURN
      END IF

        X_val=XOLD+T*(X-XOLD)/TINT
        Y_val=YOLD+T*(Y-YOLD)/TINT
        R1=SQRT(X_val*X_val+Y_val*Y_val)
        J=MAX(1.,MIN(20.,RFACT*ALOG10(R1/R0)+1))

C          IF (ITYPE.EQ.0) THEN
C            LG(J)=LG(J)+1
C          ELSE IF (ITYPE .EQ. 1) THEN

c Loop over thinning levels
        do ith=1,NTHIN+1
          IF (ITYPE.EQ.1) THEN
            XLE(J,NR,ith)=XLE(J,NR,ith)+Wgt(ith)       ! Electron count
          if (Wgt(ith).gt.0.) XLENOTHIN(J,NR,ith)=XLENOTHIN(J,NR,ith)+1.   ! Electron count no thinning
c Moments of inertia for electrons
            WXE(NR,ith)=WXE(NR,ith)+Wgt(ith)*abs(X_val)      ! Sum of weight*x 
            WYE(NR,ith)=WYE(NR,ith)+Wgt(ith)*abs(Y_val)      ! Sum of weight*y 
            WRE(NR,ith)=WRE(NR,ith)+Wgt(ith)*R1              ! Sum of weight*r 
            WX2E(NR,ith)=WX2E(NR,ith)+Wgt(ith)*X_val*X_val   ! Sum of weight*x^2 
            WY2E(NR,ith)=WY2E(NR,ith)+Wgt(ith)*Y_val*Y_val   ! Sum of weight*y^2 
            WR2E(NR,ith)=WR2E(NR,ith)+Wgt(ith)*R1*R1         ! Sum of weight*r^2 
            WZE(NR,ith)=WZE(NR,ith)+Wgt(ith)*abs(Tpos)       ! Sum of weight*z 
            WZ2E(NR,ith)=WZ2E(NR,ith)+Wgt(ith)*Tpos*Tpos     ! Sum of weight*z^2 
            WSUME(NR,ith)=WSUME(NR,ith)+Wgt(ith)             ! Sum of weights 
C            TRR(J)=TRR(J)+R1
          ELSE
            XLP(J,NR,ith)=XLP(J,NR,ith)+Wgt(ith)       ! Positron count
          if (Wgt(ith).gt.0.) XLPNOTHIN(J,NR,ith)=XLPNOTHIN(J,NR,ith)+1.   ! Positron count no thinning
c Moments of inertia for positrons 
            WXP(NR,ith)=WXP(NR,ith)+Wgt(ith)*abs(X_val)      ! Sum of weight*y 
            WYP(NR,ith)=WYP(NR,ith)+Wgt(ith)*abs(Y_val)      ! Sum of weight*y 
            WRP(NR,ith)=WRP(NR,ith)+Wgt(ith)*R1              ! Sum of weight*r 
            WX2P(NR,ith)=WX2P(NR,ith)+Wgt(ith)*X_val*X_val   ! Sum of weight*x^2 
            WY2P(NR,ith)=WY2P(NR,ith)+Wgt(ith)*Y_val*Y_val   ! Sum of weight*y^2 
            WR2P(NR,ith)=WR2P(NR,ith)+Wgt(ith)*R1*R1         ! Sum of weight*r^2 
            WZP(NR,ith)=WZP(NR,ith)+Wgt(ith)*abs(Tpos)       ! Sum of weight*z 
            WZ2P(NR,ith)=WZ2P(NR,ith)+Wgt(ith)*Tpos*Tpos     ! Sum of weight*z^2 
            WSUMP(NR,ith)=WSUMP(NR,ith)+Wgt(ith)             ! Sum of weights 
C            TRR(J)=TRR(J)+R1
          END IF
C        T97(J)=T97(J)+E2      ! Energy sum
       end do

      END DO
      RETURN
      
      END

C ******************************************************************
c Electric field from an individual track in the frequency domain      
c a la ZHS - Fraunhofer
C ******************************************************************
      SUBROUTINE EMPSUM_W(ITYP,ZSTART,DELTAZ,XSTART,DELTAX,
     # YSTART,DELTAY,dCTSTA,DdCT,ZSUM_X,ZSUM_Y,ZSUM_Z,ZSUM,Wgt0,ith)
C ******************************************************************
C ------------------------------------------------------------------
C Electric field calculation following ZHS (1992) PRD paper
C ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (PI=3.141592653589793238)
      PARAMETER (maxnu=2000)                   ! Max. number of frequencies
      PARAMETER (MAXJ=50)                     ! Max. number of theta angles
      PARAMETER (MAXK=5)                      ! Max. number of phi angles
      PARAMETER (NTHIN_MAX=2)                ! Number of simultaneous thinning levels 

      REAL ZSTART,DELTAZ,XSTART,DELTAX,dCTSTA,YSTART,DELTAY
      DOUBLE PRECISION Wgt0,CONST1
c      DOUBLE PRECISION CONST2,CONST3
      INTEGER ith 
      COMPLEX*16 ZSUM(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM_X(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM_Y(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM_Z(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 PHASE1,DPHASE

      COMMON /GRID/ FRQCOM(MAXNU), THETA(MAXJ), PHI(MAXK), NUMAX,
     # JEND(MAXNU),JMAX,KMAX,JCHER,J90(MAXNU),J180(MAXNU)
      COMMON / SPEEDUP / SINTET(MAXJ),COSTET(MAXJ),SINTET2(MAXJ),
     #  COSPHI(MAXK),SINPHI(MAXK),
     #  COSTET2(MAXJ),SINMU(MAXJ),CSMU_1(MAXJ),COSMU(MAXJ),FACFRQ(MAXNU)

      COMMON / TRACKS_W / IWRITE_W

      integer*8 NCALLS_EMPSUM_W
      COMMON / THIN_QUALITY_W / NCALLS_EMPSUM_W(NTHIN_MAX)

      if (Wgt0.le.1.e-10) return    ! Do not compute EMPSUM when Wgt0 = 0 (protected for accuracy)

      NCALLS_EMPSUM_W(ith)=NCALLS_EMPSUM_W(ith)+1  ! # Calls to EMPSUM for particles with non-zero weight 

      CONST1=-Wgt0*ITYP
      DO J=1,JMAX    ! Loop in theta
      
        DO K=1,KMAX  ! Loop in phi
c        write(*,*) '1 ',k,cosphi(k),sinphi(k)

        ADVdCT=SINMU(J)*COSPHI(K)*XSTART
     #        +SINMU(J)*SINPHI(K)*YSTART
     #        +CSMU_1(J)*ZSTART-dCTSTA
        DELdCT=SINMU(J)*COSPHI(K)*DELTAX
     #        +SINMU(J)*SINPHI(K)*DELTAY
     #        +CSMU_1(J)*DELTAZ-DdCT


c Perpendicular TCK is in g/cm^2 
c WARNING
c k x (k x delta r) is actually -tck_perp and NOT tck_perp
c which accounts for the minus sign needed in -q beta_perp for q>0
        UX = SINTET(J)*COSPHI(K)
        UY = SINTET(J)*SINPHI(K)
        UZ = COSTET(J)

        TCKX=-(UY*UY+UZ*UZ)*DELTAX + UX*UY*DELTAY + UX*UZ*DELTAZ
        TCKY=UX*UY*DELTAX - (UX*UX+UZ*UZ)*DELTAY + UY*UZ*DELTAZ
        TCKZ=UX*UZ*DELTAX + UY*UZ*DELTAY - (UX*UX+UY*UY)*DELTAZ

        TCK=SQRT(TCKX*TCKX + TCKY*TCKY + TCKZ*TCKZ)
c        IF (TCKZ.LE.0.) TCK=-TCK 

c NOTE Gives: TCK=DELTAZ*SINTET(J)-DELTAX*COSTET(J) when PHI=0 deg.


      DO NU=1,NUMAX
          PH1IMG=FACFRQ(NU)*ADVdCT
          DPHIMG=FACFRQ(NU)*DELdCT
          PHASE1=PH1IMG*(0.D0,1.D0)
          DPHASE=DPHIMG*(0.D0,1.D0)
  
c Loop over thinning levels
            IF (ABS(DPHIMG) .LT. 1.D-8) THEN

c        CONST2=CONST1*TCK
              ZSUM_X(J,K,NU,ith)=
     #       ZSUM_X(J,K,NU,ith)+CONST1*TCKX*FACFRQ(NU)*EXP(PHASE1)*
     #       ((0.5D0,0.D0)*DPHIMG-(0.D0,1.D0))

              ZSUM_Y(J,K,NU,ith)=
     #       ZSUM_Y(J,K,NU,ith)+CONST1*TCKY*FACFRQ(NU)*EXP(PHASE1)*
     #       ((0.5D0,0.D0)*DPHIMG-(0.D0,1.D0))

              ZSUM_Z(J,K,NU,ith)=
     #       ZSUM_Z(J,K,NU,ith)+CONST1*TCKZ*FACFRQ(NU)*EXP(PHASE1)*
     #       ((0.5D0,0.D0)*DPHIMG-(0.D0,1.D0))

c Note ZSUM is wrong (especially at theta=0 deg.)
c due to modulus of TCK which is always positive
c so that cancellations of components are not correct
c RADIUS and PHASE however should be ok 
              ZSUM(J,K,NU,ith)=
     #       ZSUM(J,K,NU,ith)+CONST1*TCK*FACFRQ(NU)*EXP(PHASE1)*
     #       ((0.5D0,0.D0)*DPHIMG-(0.D0,1.D0))

          ELSE
c              CONST3=CONST2/DELdCT
              ZSUM_X(J,K,NU,ith)=
     #         ZSUM_X(J,K,NU,ith)+(CONST1*TCKX/DELdCT)*
     #         EXP(PHASE1)*((1.D0,0.D0)-EXP(DPHASE))

              ZSUM_Y(J,K,NU,ith)=
     #         ZSUM_Y(J,K,NU,ith)+(CONST1*TCKY/DELdCT)*
     #         EXP(PHASE1)*((1.D0,0.D0)-EXP(DPHASE))

              ZSUM_Z(J,K,NU,ith)=
     #         ZSUM_Z(J,K,NU,ith)+(CONST1*TCKZ/DELdCT)*
     #         EXP(PHASE1)*((1.D0,0.D0)-EXP(DPHASE))

c Note ZSUM is wrong (especially at theta=0 deg.)
c due to modulus of TCK which is always positive
c so that cancellations of components are not correct
c RADIUS and PHASE however should be ok 
              ZSUM(J,K,NU,ith)=
     #         ZSUM(J,K,NU,ith)+(CONST1*TCK/DELdCT)*
     #         EXP(PHASE1)*((1.D0,0.D0)-EXP(DPHASE))

            END IF

          END DO   ! End loop in freq 
        END DO     ! End loop in phi 
      END DO       ! End loop in theta

c WARNING: If iwrite_w is set to 1 a huge file can be created (not recommended above 10 TeV primary energy)
      if (IWRITE_W.eq.1) then  ! Info on tracks written in separate files for each thinning level
         write(39+ith,21) 
     # ityp,xstart,ystart,zstart,deltax,deltay,deltaz,dCTSTA,DdCT,Wgt0
21    format(1X,I3,1P,9E12.4)
      end if

      RETURN
      END


C ******************************************************************
c Electric field from individual track in the frequency domain
c using absolute times and positions of tracks
c
c Mimics what is done in the original ZHS
c
c Fraunhofer
C ******************************************************************
      SUBROUTINE EMPSUM_W_2(ITYP,Z1,Z2,X1,X2,Y1,Y2,
     #  CT1,CT2,E1,E2,ZSUM2_X,ZSUM2_Y,ZSUM2_Z,ZSUM2,Wgt0,ith)
C ******************************************************************
C ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (PI=3.141592653589793238)
      PARAMETER (maxnu=2000)                   ! Max. number of frequencies
      PARAMETER (MAXJ=50)                     ! Max. number of theta angles
      PARAMETER (MAXK=5)                      ! Max. number of phi angles
      PARAMETER (NTHIN_MAX=2)                ! Number of simultaneous thinning levels 
      PARAMETER (e_M=0.51099906)      

      REAL Z1,Z2,X1,X2,Y1,Y2,E1,E2
      REAL CT1,CT2
      DOUBLE PRECISION Wgt0,CONST1
      INTEGER ith 
      COMPLEX*16 ZSUM2_X(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_Y(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_Z(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2(MAXJ,MAXK,MAXNU,NTHIN_MAX)
      COMPLEX*16 PHASE1,DPHASE
      COMPLEX*16 EPHASE1,EDPHASE

      COMMON /GRID/ FRQCOM(MAXNU), THETA(MAXJ), PHI(MAXK), NUMAX,
     # JEND(MAXNU),JMAX,KMAX,JCHER,J90(MAXNU),J180(MAXNU)
      COMMON / SPEEDUP / SINTET(MAXJ),COSTET(MAXJ),SINTET2(MAXJ),
     #  COSPHI(MAXK),SINPHI(MAXK),
     #  COSTET2(MAXJ),SINMU(MAXJ),CSMU_1(MAXJ),COSMU(MAXJ),FACFRQ(MAXNU)

      common / zob / zobs

      COMMON / TRACKS_W_2 / IWRITE_W_2
      COMMON / SUBTRACKS / NZEROS,NIN,NCROSS1,NCROSS 

      integer*8 NCALLS_EMPSUM_W_2
      COMMON / THIN_QUALITY_W_2 / NCALLS_EMPSUM_W_2(NTHIN_MAX)

      if (Wgt0.le.1.e-10) return    ! Do not compute EMPSUM when Wgt0 = 0 (protected for accuracy)

      NCALLS_EMPSUM_W_2(ith)=NCALLS_EMPSUM_W_2(ith)+1  ! # Calls to EMPSUM for particles with non-zero weight 

c ######      
c  SIGN 
c ######      
c According to Feynman's formula a positive charge q>0
c travelling along z>0 produces an electric field      
c proportional to -q beta_perp

c Minus sign because ITYP=1 is an electron
c                    ITYP=-1 is a positron 
      CONST1=-Wgt0*ITYP

c Average beta in track
      E_AV=(E1+E2)/2.
      G=E_AV/e_M
      G2=G*G
      BETA=SQRT(1.-1./G2)
c BETA is parallel to the vector (r2-r1)
c Compute unitary vector parallel to (r2-r1)
      R2_R1=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
      IF (R2_R1.EQ.0) THEN
c        NZEROS=NZEROS+1
        RETURN
      END IF 
      BETA_X=BETA*(X2-X1)/R2_R1
      BETA_Y=BETA*(Y2-Y1)/R2_R1
      BETA_Z=BETA*(Z2-Z1)/R2_R1

c Assume azimuth observation angle phi=0

      DO J=1,JMAX   ! Loop in theta

        DO K=1,KMAX   ! Loop in phi
c        write(*,*) '2 ',k,cosphi(k),sinphi(k)
              
        CHERFAC_1 = CT1
     #              - SINMU(J)*COSPHI(K)*X1
     #              - SINMU(J)*SINPHI(K)*Y1
     #              - COSMU(J)*Z1              ! g cm^-2 
        DENOM= 1. 
     #         - SINMU(J)*COSPHI(K)*BETA_X
     #         - SINMU(J)*SINPHI(K)*BETA_Y
     #         - COSMU(J)*BETA_Z

c k x (k x beta) 
        UX = SINTET(J)*COSPHI(K)
        UY = SINTET(J)*SINPHI(K)
        UZ = COSTET(J)

        BETA_PERP_X=-(UY*UY+UZ*UZ)*BETA_X + UX*UY*BETA_Y + UX*UZ*BETA_Z
        BETA_PERP_Y=UX*UY*BETA_X - (UX*UX+UZ*UZ)*BETA_Y + UY*UZ*BETA_Z
        BETA_PERP_Z=UX*UZ*BETA_X + UY*UZ*BETA_Y - (UX*UX+UY*UY)*BETA_Z

c WARNING
c k x (k x beta) is actually -beta_perp and NOT beta_perp 
c which accounts for the minus sign needed in -q beta_perp for q>0
        BETA_PERP=sqrt(BETA_PERP_X*BETA_PERP_X+
     #                 BETA_PERP_Y*BETA_PERP_Y+
     #                 BETA_PERP_Z*BETA_PERP_Z)


c FACFRQ = 2*pi*frequency_MHz/(c*rho)
      DO NU=1,NUMAX
          PHI1=FACFRQ(NU)*CHERFAC_1
          DPHI=FACFRQ(NU)*(CT2-CT1)*DENOM
          PHASE1=PHI1*(0.D0,1.D0)
          DPHASE=DPHI*(0.D0,1.D0)
        
          EDPHASE=EXP(DPHASE)
          EPHASE1=EXP(PHASE1)

          IF (ABS(DPHI).LT. 1.D-8) THEN
            ZSUM2_X(J,K,NU,ith) = ZSUM2_X(J,K,NU,ith) + 
     #CONST1*BETA_PERP_X*FACFRQ(NU)*(CT2-CT1)*EPHASE1*(0.D0,1.D0)

            ZSUM2_Y(J,K,NU,ith) = ZSUM2_Y(J,K,NU,ith) + 
     #CONST1*BETA_PERP_Y*FACFRQ(NU)*(CT2-CT1)*EPHASE1*(0.D0,1.D0)

            ZSUM2_Z(J,K,NU,ith) = ZSUM2_Z(J,K,NU,ith) + 
     #CONST1*BETA_PERP_Z*FACFRQ(NU)*(CT2-CT1)*EPHASE1*(0.D0,1.D0)

c Note ZSUM2 is wrong (especially at theta=0 deg.)
c due to modulus of BETA_PERP which is always positive
c so that cancellations of components are not correct
c ARADIUS and APHASE however should be ok 
            ZSUM2(J,K,NU,ith) = ZSUM2(J,K,NU,ith) + 
     #CONST1*BETA_PERP*FACFRQ(NU)*(CT2-CT1)*EPHASE1*(0.D0,1.D0)

          ELSE

            ZSUM2_X(J,K,NU,ith) = ZSUM2_X(J,K,NU,ith) + 
     #       CONST1*BETA_PERP_X*EPHASE1*(EDPHASE-(1.D0,0.D0))/DENOM

            ZSUM2_Y(J,K,NU,ith) = ZSUM2_Y(J,K,NU,ith) + 
     #       CONST1*BETA_PERP_Y*EPHASE1*(EDPHASE-(1.D0,0.D0))/DENOM

            ZSUM2_Z(J,K,NU,ith) = ZSUM2_Z(J,K,NU,ith) + 
     #       CONST1*BETA_PERP_Z*EPHASE1*(EDPHASE-(1.D0,0.D0))/DENOM

c Note ZSUM2 is wrong (especially at theta=0 deg.)
c due to modulus of BETA_PERP which is always positive
c so that cancellations of components are not correct
c ARADIUS and APHASE however should be ok 
            ZSUM2(J,K,NU,ith) = ZSUM2(J,K,NU,ith) + 
     #       CONST1*BETA_PERP*EPHASE1*(EDPHASE-(1.D0,0.D0))/DENOM
          END IF

          END DO    ! End loop in frequency
        END DO     ! End loop in theta
       END DO      ! End loop in phi

c WARNING: If iwrite_w is set to 1 a huge file can be created (not recommended above 10 TeV primary energy)
      if (IWRITE_W_2.eq.1 .and. z2 .ge. zobs) then  ! Info on tracks written in separate files for each thinning level
         write(59+ith,31) 
     # ityp,Wgt0,E1,x1,y1,z1,E2,x2,y2,z2,CT1,ct2
31    format(1X,I5,1P,11E18.8)
      end if



      RETURN
      END

C ******************************************************************
c Vector potential from individual track in the time domain
c using absolute times and positions of tracks - Fraunhofer
C ******************************************************************
      SUBROUTINE EMPSUM_T(ITYP,Z1,Z2,X1,X2,
     #   Y1,Y2,CT1,CT2,E1,E2,ASUM_X,ASUM_Y,ASUM_Z,Wgt0,ith)
C ******************************************************************
C ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (PI=3.141592653589793238)
      PARAMETER (maxnu=2000)                   ! Max. number of frequencies
      PARAMETER (MAXJ=50)                     ! Max. number of theta angles
      PARAMETER (MAXK=5)                      ! Max. number of phi angles
      PARAMETER (MAXT=10000)                  ! Max. number of time bins 
      PARAMETER (MAXTHALF=MAXT/2)              
      PARAMETER (NTHIN_MAX=2)                ! Number of simultaneous thinning levels 
      PARAMETER (e_M=0.51099906)      

      DOUBLE PRECISION ASUM_X(MAXT,MAXJ,MAXK,NTHIN_MAX)
      DOUBLE PRECISION ASUM_Y(MAXT,MAXJ,MAXK,NTHIN_MAX)
      DOUBLE PRECISION ASUM_Z(MAXT,MAXJ,MAXK,NTHIN_MAX)
      REAL Z1,Z2,X1,X2,Y1,Y2,E1,E2
      REAL CT1,CT2
c      DOUBLE PRECISION Wgt0,CONST1_T
      INTEGER ith 

      COMMON /GRID/ FRQCOM(MAXNU), THETA(MAXJ), PHI(MAXK), NUMAX,
     # JEND(MAXNU),JMAX,KMAX,JCHER,J90(MAXNU),J180(MAXNU)
      COMMON / SPEEDUP / SINTET(MAXJ),COSTET(MAXJ),SINTET2(MAXJ),
     #  COSPHI(MAXK),SINPHI(MAXK),      
     #  COSTET2(MAXJ),SINMU(MAXJ),CSMU_1(MAXJ),COSMU(MAXJ),FACFRQ(MAXNU)

      COMMON /CONSTANTS_T/ FACTOR_T
      COMMON /BIN_T/ DT_BIN_NS,DT_BIN_S 

      COMMON / TRACKS_T / IWRITE_T
      COMMON / SUBTRACKS / NZEROS,NIN,NCROSS1,NCROSS 

      integer*8 NCALLS_EMPSUM_T
      COMMON / THIN_QUALITY_T / NCALLS_EMPSUM_T(NTHIN_MAX)

      if (Wgt0.le.1.d-10) return    ! Do not compute EMPSUM when Wgt0 = 0 (protected for accuracy)

c ######      
c  SIGN 
c ######      
c According to Feynman's formula a positive charge q>0
c travelling along z>0 produces an electric field      
c proportional to -q beta_perp or equivalently
c the vector potential (computed in this routine) 
c is proportional to q beta_perp.
c For instance for observation of a single track
c with q>0 at theta = 90 deg, the z-component
c of the vector potential should be _ a positive
c step followed by a negative one _| |_

c Factor 2 comes from the sum of the 2 sign functions after 
c the transformation to the time domain
c Minus sign because ITYP=1 is an electron
c                    ITYP=-1 is a positron 
      CONST1=-2.d0*Wgt0*DBLE(ITYP)

c Average beta in track
      E_AV=(E1+E2)/2.
      G=E_AV/e_M
      G2=G*G
      BETA=SQRT(1.-1./G2)
c BETA is parallel to the vector (r2-r1)
c Compute unitary vector parallel to (r2-r1)
      R2_R1=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
      IF (R2_R1.EQ.0) THEN
        NZEROS=NZEROS+1
        RETURN
      END IF 
      BETA_X=BETA*(X2-X1)/R2_R1
      BETA_Y=BETA*(Y2-Y1)/R2_R1
      BETA_Z=BETA*(Z2-Z1)/R2_R1

c Assume azimuth observation angle phi=0

c FACTOR_T = 1/(c*rho)
      DO J=1,JMAX    ! Loop in theta

        DO K=1,KMAX  ! Loop in phi

c        write(*,*) '3 ',k,cosphi(k),sinphi(k)
        AUX1 = SINMU(J)*COSPHI(K)*X1 
     #       + SINMU(J)*SINPHI(K)*Y1 
     #       + COSMU(J)*Z1 
        AUXB = SINMU(J)*COSPHI(K)*BETA_X 
     #       + SINMU(J)*SINPHI(K)*BETA_Y 
     #       + COSMU(J)*BETA_Z
        DENOM = 1.d0 - AUXB 

c        T1=FACTOR_T*CT1
c        T2=FACTOR_T*CT2
c        DT1=T1*DENOM
c        DT2=T2*DENOM

        CHERFAC_1 = CT1 - AUX1                      ! g cm^-2 
        DELTAT_1 = FACTOR_T*CHERFAC_1               ! s
        CHERFAC_2 = CT2 - AUXB*(CT2-CT1) - AUX1     ! g cm^-2 
        DELTAT_2 = FACTOR_T*CHERFAC_2               ! s

c        DT = DELTAT_2-DELTAT_1                      ! s
c        DTT = FACTOR_T*(CT2-CT1)*DENOM              ! s
        DTT = FACTOR_T*(CT2-CT1)                     ! s

c k x (k x beta) 
        UX = SINTET(J)*COSPHI(K)
        UY = SINTET(J)*SINPHI(K)
        UZ = COSTET(J)

        BETA_PERP_X=-(UY*UY+UZ*UZ)*BETA_X + UX*UY*BETA_Y + UX*UZ*BETA_Z
        BETA_PERP_Y=UX*UY*BETA_X - (UX*UX+UZ*UZ)*BETA_Y + UY*UZ*BETA_Z
        BETA_PERP_Z=UX*UZ*BETA_X + UY*UZ*BETA_Y - (UX*UX+UY*UY)*BETA_Z

        BETA_PERP=SQRT(BETA_PERP_X*BETA_PERP_X +
     #                 BETA_PERP_Y*BETA_PERP_Y +
     #                 BETA_PERP_Z*BETA_PERP_Z)

c WARNING
c     k x (k x beta) is actually -beta_perp and NOT beta_perp => change signs
        BETA_PERP_X=-BETA_PERP_X
        BETA_PERP_Y=-BETA_PERP_Y
        BETA_PERP_Z=-BETA_PERP_Z


c Find time bins in which the contribution to vector potential is non-zero
        IF (DELTAT_1.LT.0.) THEN
          IT1=MAXTHALF + INT(DELTAT_1/DT_BIN_S)-1
        ELSE
          IT1=MAXTHALF + INT(DELTAT_1/DT_BIN_S)
        END IF

        IF (DELTAT_2.LT.0.) THEN
          IT2=MAXTHALF + INT(DELTAT_2/DT_BIN_S)-1
        ELSE
          IT2=MAXTHALF + INT(DELTAT_2/DT_BIN_S)
        END IF
c        if (it1.ne.it2) write(*,*) DELTAT_1*1.d9,DELTAT_2*1.d9,IT1,IT2

        IF (IT1.LE.IT2) THEN 
           IT_STA=IT1
           IT_END=IT2
        ELSE
           IT_STA=IT2
           IT_END=IT1
        END IF   

c Get rid of times outside time range
        IF (IT_STA.GT.MAXT.OR.IT_END.LT.1) GOTO 200  ! Do not compute field
        IF (IT_STA.LT.1) IT_STA=1
        IF (IT_END.GT.MAXT) IT_END=MAXT

        CONST2=CONST1
c Change of sign when DELTAT_2 < DELTAT_1
        IF (DELTAT_2.LT.DELTAT_1) THEN
           CONST2=-CONST1
c        write(*,*) DELTAT_2*1.d9-DELTAT_1*1.d9,IT1,IT2,ITYP,CONST2 
        END IF 

c Count types of subtracks        
        IF(IT_STA.EQ.IT_END) THEN
            NIN=NIN+1
        ELSE IF ((IT_END-IT_STA).eq.1) THEN
            NCROSS1=NCROSS1+1
        ELSE 
            NCROSS=NCROSS+1
        END IF

c Note ASUM_i is multiplied by e*mu_r/(8.*pi*epsilon_0*c)
c when writing it to the file. One factor 1/c goes into beta 

c Distribute vector potential among bins in time
c and calculate average vector potential in each bin        
c (see my notes: MC implementation)
c Since we take the derivative at the end the important 
c point is to conserve the structure from bin to bin.
        IF (DELTAT_2.LT.DELTAT_1) THEN 
            DELTAT_STA=DELTAT_2
            DELTAT_END=DELTAT_1
        ELSE
            DELTAT_STA=DELTAT_1
            DELTAT_END=DELTAT_2
        END IF


c        if (abs(denom).lt.1.d-15) then 
c          write(*,*) j,denom,auxb,beta_x,beta_z
c          write(*,*) it_sta,it_end,abs(dtt)
c        end if


        DO IT=IT_STA,IT_END,1
        
          IF(IT_STA.EQ.IT_END) THEN   ! Subtrack contained in bin

           if (it_sta.gt.maxt.or.it_end.gt.maxt) write(*,*)
     #     'IT_STA>MAXT or IT_END<MAXT ',      
     #       j,it_sta,it_end,delta_sta,delta_end

          F=(DELTAT_END-DELTAT_STA)/DT_BIN_S

            IF (ABS(DENOM).GT.1.D-15) THEN
              ASUM_X(IT,J,K,ith) = ASUM_X(IT,J,K,ith) + 
     #                         ABS(F)*CONST2*BETA_PERP_X/DENOM
              ASUM_Y(IT,J,K,ith) = ASUM_Y(IT,J,K,ith) + 
     #                         ABS(F)*CONST2*BETA_PERP_Y/DENOM
              ASUM_Z(IT,J,K,ith) = ASUM_Z(IT,J,K,ith) + 
     #                         ABS(F)*CONST2*BETA_PERP_Z/DENOM
c             write(*,*)
c     # DELTAT_1*1.d9,DELTAT_2*1.d9,IT,IT_STA,IT_END,ABS(F)

             ELSE   ! Cherenkov angle approximation
              ASUM_X(IT,J,K,ith) = ASUM_X(IT,J,K,ith) + 
     #                         ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_X
              ASUM_Y(IT,J,K,ith) = ASUM_Y(IT,J,K,ith) + 
     #                         ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_Y
              ASUM_Z(IT,J,K,ith) = ASUM_Z(IT,J,K,ith) + 
     #                         ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_Z
 
             END IF

          ELSE  ! Subtrack crossing one or more bins in time

c Correct for bin edges               

               F=((IT_STA+1-MAXTHALF)*DT_BIN_S-DELTAT_STA)/DT_BIN_S
c (Bug corrected March 2010)
c               IF (ABS(DENOM).GT.1.D-15) THEN  
c We should also account for ABS(DENOM)<1.D-15 in this case               
             IF (IT.EQ.IT_STA) THEN ! Start of subtrack 
               
                 ASUM_X(IT,J,K,ith) = ASUM_X(IT,J,K,ith) + 
     #                          ABS(F)*CONST2*BETA_PERP_X/DENOM
                 ASUM_Y(IT,J,K,ith) = ASUM_Y(IT,J,K,ith) + 
     #                          ABS(F)*CONST2*BETA_PERP_Y/DENOM
                 ASUM_Z(IT,J,K,ith) = ASUM_Z(IT,J,K,ith) + 
     #                          ABS(F)*CONST2*BETA_PERP_Z/DENOM
c             write(*,*)
c     # DELTAT_STA*1.d9,DELTAT_END*1.d9,IT,IT_STA,IT_END,ABS(F)

             ELSE IF (IT.EQ.IT_END) THEN  ! End of subtrack

               F=((IT_END+1-MAXTHALF)*DT_BIN_S-DELTAT_END)/DT_BIN_S
               ASUM_X(IT,J,K,ith) = ASUM_X(IT,J,K,ith) + 
     #                        (1.d0-ABS(F))*CONST2*BETA_PERP_X/DENOM
               ASUM_Y(IT,J,K,ith) = ASUM_Y(IT,J,K,ith) + 
     #                        (1.d0-ABS(F))*CONST2*BETA_PERP_Y/DENOM
               ASUM_Z(IT,J,K,ith) = ASUM_Z(IT,J,K,ith) + 
     #                        (1.d0-ABS(F))*CONST2*BETA_PERP_Z/DENOM
c             write(*,*)
c     # DELTAT_STA*1.d9,DELTAT_END*1.d9,IT,IT_STA,IT_END,(1.-ABS(F))

             ELSE

               F=1.d0
               ASUM_X(IT,J,K,ith) = ASUM_X(IT,J,K,ith) +
     #                          F*CONST2*BETA_PERP_X/DENOM
               ASUM_Y(IT,J,K,ith) = ASUM_Y(IT,J,K,ith) +
     #                          F*CONST2*BETA_PERP_Y/DENOM
               ASUM_Z(IT,J,K,ith) = ASUM_Z(IT,J,K,ith) +
     #                          F*CONST2*BETA_PERP_Z/DENOM

             END IF   

          END IF

        END DO    

200    CONTINUE

         END DO      ! End loop in phi 
       END DO      ! End loop in theta

c WARNING: If iwrite_t is set to 1 a huge file can be created (not recommended above 10 TeV primary energy)
      if (IWRITE_T.eq.1) then  ! Info on tracks written in separate files for each thinning level
         write(69+ith,31) 
     # ityp,Wgt0,E1,x1,y1,z1,E2,x2,y2,z2,CT1,CT2
31    format(1X,I5,1P,11E18.8)
      end if


      RETURN
      END

C ******************************************************************
c Electric field from individual track in the frequency domain
c using absolute times and positions of tracks
c
c Fresnel regime. Computes field in a grid of observers
c placed at (x,y,z) in the ZHS reference frame.
c
c FT of Washington's implementation in time-domain (see notes)
c
c Fresnel      
C ******************************************************************
      SUBROUTINE EMPSUM_W_2_FRESNEL(ITYP,Z1,Z2,X1,X2,Y1,Y2,
     #  CT1,CT2,E1,E2,ZSUM2_X_FR,ZSUM2_Y_FR,ZSUM2_Z_FR,ZSUM2_ANT_FR,
     #  Wgt0,ith)
C ******************************************************************
C ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (PI=3.141592653589793238)
      PARAMETER (maxnu=2000)                   ! Max. number of frequencies
      PARAMETER (MAXJ=50)                     ! Max. number of theta angles
      PARAMETER (MAXK=5)                      ! Max. number of phi angles
      PARAMETER (MAXNA=200)                     ! Number of antenna positions
      PARAMETER (NTHIN_MAX=2)                ! Number of simultaneous thinning levels 
      PARAMETER (e_M=0.51099906)      

      REAL Z1,Z2,X1,X2,Y1,Y2,E1,E2
      real zmin1,zmax1
      real zmin2,zmax2
      real zmin3,zmax3
      REAL CT1,CT2
      DOUBLE PRECISION Wgt0,CONST1
      INTEGER ith 
      COMPLEX*16 ZSUM2_X_FR(MAXNA,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_Y_FR(MAXNA,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_Z_FR(MAXNA,MAXNU,NTHIN_MAX)
      COMPLEX*16 ZSUM2_ANT_FR(MAXNA,MAXNU,NTHIN_MAX)
      COMPLEX*16 PHASE1,PHASE2
      COMPLEX*16 EPHASE1,EPHASE2,DEPHASE
      complex*16 tmp


      COMMON /GRID/ FRQCOM(MAXNU), THETA(MAXJ), PHI(MAXK), NUMAX,
     # JEND(MAXNU),JMAX,KMAX,JCHER,J90(MAXNU),J180(MAXNU)

      COMMON / SPEEDUP / SINTET(MAXJ),COSTET(MAXJ),SINTET2(MAXJ),
     #  COSPHI(MAXK),SINPHI(MAXK),
     #  COSTET2(MAXJ),SINMU(MAXJ),CSMU_1(MAXJ),COSMU(MAXJ),FACFRQ(MAXNU)

      COMMON /ANT1/ NAMAX
      COMMON /ANT2/ XANT(MAXNA),YANT(MAXNA),ZANT(MAXNA)
      COMMON /ANT3/ ETA_XANT(MAXNA),ETA_YANT(MAXNA),ETA_ZANT(MAXNA)
      COMMON /ANT4/ ETA_PERP_XANT(MAXNA),ETA_PERP_YANT(MAXNA),
     #              ETA_PERP_ZANT(MAXNA)

      COMMON /REFR/ REFIDX

      COMMON /CONSTANTS_T/ FACTOR_T

      COMMON /FREQ_ARRAY/ FREQ(MAXNU)

      if (Wgt0.le.1.e-10) return    ! Do not compute EMPSUM when Wgt0 = 0 (protected for accuracy)

c #################################################################
c #################################################################
c #################################################################
c Compute field for tracks contained between zmin and zmax (g/cm^2)
c Divide camera in 4 pieces from 0 - 0.468 g/cm^2
c First piece
      zmin1=0.0
      zmax1=0.117
c      if (z1.lt.zmin1.or.z1.gt.zmax1.or.z2.lt.zmin1.or.z2.gt.zmax1) 
c     #return
c Second piece
      zmin2=0.117
      zmax2=2.*0.117
c      if (z1.lt.zmin2.or.z1.gt.zmax2.or.z2.lt.zmin2.or.z2.gt.zmax2) 
c     #return
c Third piece
      zmin3=2.*0.117
      zmax3=3.*0.117
c      if (z1.lt.zmin3.or.z1.gt.zmax3.or.z2.lt.zmin3.or.z2.gt.zmax3) 
c     #return
c Fourth piece
      zmin4=3.*0.117
      zmax4=4.*0.117
c      if (z1.lt.zmin4.or.z1.gt.zmax4.or.z2.lt.zmin4.or.z2.gt.zmax4) 
c     #return

c Remove tracks crossing boundaries of camera pieces 
c even when the whole camera is accounted for
c      if (z1.lt.zmax1.and.z2.gt.zmax1) return
c      if (z1.lt.zmax2.and.z2.gt.zmax2) return
c      if (z1.lt.zmax3.and.z2.gt.zmax3) return
c      if (z1.lt.zmax4.and.z2.gt.zmax4) return
c #################################################################
c #################################################################
c #################################################################


      rho_air=1.2e-3

c ######      
c  SIGN 
c ######      
c According to Feynman's formula a positive charge q>0
c travelling along z>0 produces an electric field      
c proportional to -q beta_perp

c Minus sign because ITYP=1 is an electron
c                    ITYP=-1 is a positron 
      CONST1=-Wgt0*ITYP

c Average beta in track
      E_AV=(E1+E2)/2.
      G=E_AV/e_M
      G2=G*G
      BETA=SQRT(1.-1./G2)
c BETA is parallel to the vector (r2-r1)
c Compute unitary vector parallel to (r2-r1)
      R2_R1=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
      IF (R2_R1.EQ.0) THEN
c        NZEROS=NZEROS+1
        RETURN
      END IF 
      BETA_X=BETA*(X2-X1)/R2_R1
      BETA_Y=BETA*(Y2-Y1)/R2_R1
      BETA_Z=BETA*(Z2-Z1)/R2_R1

c ----------------------------------------------------------
      DO NA=1,NAMAX   ! Loop in antenna positions 

c Vector k from average point of track to observer 
        UX = XANT(NA)-(X1+X2)/2.     
        UY = YANT(NA)-(Y1+Y2)/2.
        UZ = ZANT(NA)-(Z1+Z2)/2.       ! g cm^-2

        UX1 = XANT(NA)-X1     
        UY1 = YANT(NA)-Y1
        UZ1 = ZANT(NA)-Z1       ! g cm^-2

        UX2 = XANT(NA)-X2     
        UY2 = YANT(NA)-Y2
        UZ2 = ZANT(NA)-Z2       ! g cm^-2

c Distance from average point of track to observer 
        R = SQRT(UX*UX+UY*UY+UZ*UZ)  ! g cm^-2
c Distance from beginning of track to observer 
        R1 = SQRT(UX1*UX1+UY1*UY1+UZ1*UZ1)  ! g cm^-2
c Distance from end of track to observer 
        R2 = SQRT(UX2*UX2+UY2*UY2+UZ2*UZ2)  ! g cm^-2

c Antenna position w.r.t. injection point of the shower        
        RANT = SQRT(XANT(NA)*XANT(NA)+
     #              YANT(NA)*YANT(NA)+
     #              ZANT(NA)*ZANT(NA))  ! g cm^-2 

c Unit vector k  
        UX=UX/R
        UY=UY/R
        UZ=UZ/R

c Unit vector k1  
        UX1=UX1/R1
        UY1=UY1/R1
        UZ1=UZ1/R1

c Unit vector k2 
        UX2=UX2/R2
        UY2=UY2/R2
        UZ2=UZ2/R2

c        write(*,'(1X,5F14.3)') 
c     #xant(na),yant(na),zant(na),R,acos(uz)*180./pi

        AUXB = REFIDX*UX*BETA_X 
     #       + REFIDX*UY*BETA_Y 
     #       + REFIDX*UZ*BETA_Z

        AUX1 = REFIDX*UX1*BETA_X 
     #       + REFIDX*UY1*BETA_Y 
     #       + REFIDX*UZ1*BETA_Z


c        DENOM= 1.d0 - AUX1 
        DENOM= 1.d0 - AUXB 

c        write(101,*) na,z1,ux,uy,uz,acos(auxb/refidx/beta)*180./pi

c k x (k x beta) 
        BETA_PERP_X=-(UY*UY+UZ*UZ)*BETA_X + UX*UY*BETA_Y + UX*UZ*BETA_Z
        BETA_PERP_Y=UX*UY*BETA_X - (UX*UX+UZ*UZ)*BETA_Y + UY*UZ*BETA_Z
        BETA_PERP_Z=UX*UZ*BETA_X + UY*UZ*BETA_Y - (UX*UX+UY*UY)*BETA_Z

c WARNING
c k x (k x beta) is actually -beta_perp and NOT beta_perp 
c which accounts for the minus sign needed in -q beta_perp for q>0
        BETA_PERP=sqrt(BETA_PERP_X*BETA_PERP_X+
     #                 BETA_PERP_Y*BETA_PERP_Y+
     #                 BETA_PERP_Z*BETA_PERP_Z)


c FACFRQ = 2*pi*frequency_MHz/(c*rho)
         DO NU=1,NUMAX
        
c iw(t1+nR/c) different for each track in Fresnel
          PHASE1=FACFRQ(NU)*(CT1+REFIDX*R1)*(0.D0,1.D0) 

c iw(t2+nR/c) different for each track in Fresnel
           PHASE2=FACFRQ(NU)*
     #(CT2+REFIDX*R1-AUXB*(CT2-CT1))*(0.D0,1.D0)  

c e^iw(t1+nR/c) different for each track in Fresnel
          EPHASE1=EXP(PHASE1)

c e^iw(t2+nR/c) different for each track in Fresnel
          EPHASE2=EXP(PHASE2)

          DEPHASE=EPHASE2-EPHASE1

c          if (nu.eq.1.and.ith.eq.1) 
c     #write(*,*) X1,X2,Y1,Y2,Z1,Z2,CT1,CT2,DEPHASE

c Set ifilter=0 if  
c Filter is NOT desired.
c Filter is linear with freq between 0 and 2 GHz and 1/freq after that.
        ifilter=1
        if (ifilter.eq.1) then   ! Filter activated

          if (freq(nu) .ge. 2000.) then

            ZSUM2_X_FR(NA,NU,ith) = ZSUM2_X_FR(NA,NU,ith) + 
     #CONST1*BETA_PERP_X*DEPHASE/DENOM/R/freq(nu)

            ZSUM2_Y_FR(NA,NU,ith) = ZSUM2_Y_FR(NA,NU,ith) + 
     #CONST1*BETA_PERP_Y*DEPHASE/DENOM/R/freq(nu)

            ZSUM2_Z_FR(NA,NU,ith) = ZSUM2_Z_FR(NA,NU,ith) + 
     #CONST1*BETA_PERP_Z*DEPHASE/DENOM/R/freq(nu)

          else

c 0.25e-6 so that filter is continuous at 2 GHz
            ZSUM2_X_FR(NA,NU,ith) = ZSUM2_X_FR(NA,NU,ith) + 
     #0.25e-6*freq(nu)*CONST1*BETA_PERP_X*DEPHASE/DENOM/R

            ZSUM2_Y_FR(NA,NU,ith) = ZSUM2_Y_FR(NA,NU,ith) + 
     #0.25e-6*freq(nu)*CONST1*BETA_PERP_Y*DEPHASE/DENOM/R

            ZSUM2_Z_FR(NA,NU,ith) = ZSUM2_Z_FR(NA,NU,ith) + 
     #0.25e-6*freq(nu)*CONST1*BETA_PERP_Z*DEPHASE/DENOM/R

           end if

         else if (ifilter.eq.0) then   ! No Filter

            ZSUM2_X_FR(NA,NU,ith) = ZSUM2_X_FR(NA,NU,ith) + 
     #CONST1*BETA_PERP_X*DEPHASE/DENOM/R

            ZSUM2_Y_FR(NA,NU,ith) = ZSUM2_Y_FR(NA,NU,ith) + 
     #CONST1*BETA_PERP_Y*DEPHASE/DENOM/R

            ZSUM2_Z_FR(NA,NU,ith) = ZSUM2_Z_FR(NA,NU,ith) + 
     #CONST1*BETA_PERP_Z*DEPHASE/DENOM/R

         end if  ! Ends filter selection

          ZSUM2_ANT_FR(NA,NU,ith)=ZSUM2_ANT_FR(NA,NU,ith)+
     #      ( (CONST1*BETA_PERP_X*DEPHASE/DENOM/R)*ETA_XANT(NA)+   
     #       (CONST1*BETA_PERP_Y*DEPHASE/DENOM/R)*ETA_YANT(NA)+   
     #       (CONST1*BETA_PERP_Z*DEPHASE/DENOM/R)*ETA_ZANT(NA) )

c           CTHETA=-UX*ETA_PERP_XANT(NA)
c     #            -UY*ETA_PERP_YANT(NA)
c     #            -UZ*ETA_PERP_ZANT(NA)

c #################################################################
c #################################################################
c #################################################################
c Write z-component of field contributed by each sub-track - for ROOT ntuple
          if (nu.eq.0.and.ith.eq.1) then  ! Uncomment this line and comment next to deactivate writing
c          if (nu.eq.1.and.ith.eq.1) then
          t_arr=((CT1+REFIDX*R1+CT2+REFIDX*R2)/2./rho_air/100./3.e8)  ! Arrival time at antenna (s)
          t_arr = t_arr*1.e9 ! in ns
c          write(*,*) ct1,refidx*r1,ct2,refidx*r2,rho_air,Tarr
          tck=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
          tmp=CONST1*BETA_PERP_Z*DEPHASE/DENOM/R 
          Ez_diff=SQRT(tmp*CONJG(tmp))/tck  ! "Differential" field - Normalized to track length
          Ez_acc=
     #    SQRT(ZSUM2_Z_FR(NA,NU,ith)*CONJG(ZSUM2_Z_FR(NA,NU,ith)))  ! "Accumulated" field 
c Assume isotropic emission (i.e. remove denom)
          tmp=CONST1*BETA_PERP_Z*DEPHASE/R  
          Ez_iso=SQRT(tmp*CONJG(tmp))/tck ! Normalized to track length
          Edep=Edep+(E2-E1)
          write(111,*) 
     #    (x1+x2)/2./rho_air/100.,  ! 1
     #    (y1+y2)/2./rho_air/100.,  ! 2
     #    (z1+z2)/2./rho_air/100.,  ! 3
     #    t_arr,                    ! 4
     #    R/rho_air/100.,           ! 5
     #    tck,                      ! 6
     #    denom,                    ! 7
     #    FREQ(NU),                 ! 8
     #    Ez_diff,                  ! 9
     #    Ez_acc,                   ! 10 
     #    Ez_iso                    ! 11
          end if  
c #################################################################
c #################################################################
c #################################################################

         END DO    ! End loop in frequency

        END DO     ! End loop in antenna positions 

      RETURN
      END


C ******************************************************************
c Vector potential from individual track in the time domain
c using absolute times and positions of tracks - Fresnel
c Washington's implementation with CT0=CT1 i.e. w.r.t. start of track
C ******************************************************************
      SUBROUTINE EMPSUM_T_FRESNEL(ITYP,Z1,Z2,X1,X2,
     #   Y1,Y2,CT1,CT2,E1,E2,ASUM_X_FR,ASUM_Y_FR,ASUM_Z_FR,ASUM_ANT_FR,
     #   Wgt0,ith)
C ******************************************************************
C ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (PI=3.141592653589793238)
      PARAMETER (maxnu=2000)
      PARAMETER (MAXNA=200)                    ! Max. number of antennas
      PARAMETER (MAXT=10000)                 ! Max. number of time bins 
      PARAMETER (MAXTHALF=MAXT/2)            
      PARAMETER (NTHIN_MAX=2)               ! Number of simultaneous thinning levels 
      PARAMETER (e_M=0.51099906)      

      DOUBLE PRECISION ASUM_X_FR(MAXT,MAXNA,NTHIN_MAX)
      DOUBLE PRECISION ASUM_Y_FR(MAXT,MAXNA,NTHIN_MAX)
      DOUBLE PRECISION ASUM_Z_FR(MAXT,MAXNA,NTHIN_MAX)
      DOUBLE PRECISION ASUM_ANT_FR(MAXT,MAXNA,NTHIN_MAX)
      REAL Z1,Z2,X1,X2,Y1,Y2,E1,E2
      REAL CT1,CT2
      REAL CT1_OBS,CT2_OBS
      INTEGER ith 

      COMMON /ANT1/ NAMAX
      COMMON /ANT2/ XANT(MAXNA),YANT(MAXNA),ZANT(MAXNA)
      COMMON /ANT3/ ETA_XANT(MAXNA),ETA_YANT(MAXNA),ETA_ZANT(MAXNA)
      COMMON /ANT4/ ETA_PERP_XANT(MAXNA),ETA_PERP_YANT(MAXNA),
     #              ETA_PERP_ZANT(MAXNA)


      COMMON /CONSTANTS_T/ FACTOR_T
      COMMON /BIN_T/ DT_BIN_NS,DT_BIN_S 

      COMMON /REFR/ REFIDX


      if (Wgt0.le.1.d-10) return    ! Do not compute EMPSUM when Wgt0 = 0 (protected for accuracy)

      rho_air=1.2e-3

c ######      
c  SIGN 
c ######      
c According to Feynman's formula a positive charge q>0
c travelling along z>0 produces an electric field      
c proportional to -q beta_perp or equivalently
c the vector potential (computed in this routine) 
c is proportional to q beta_perp.
c For instance for observation of a single track
c with q>0 at theta = 90 deg, the z-component
c of the vector potential should be _ a positive
c step followed by a negative one _| |_

c Factor 2 comes from the sum of the 2 sign functions after 
c the transformation to the time domain
c Minus sign because ITYP=1 is an electron
c                    ITYP=-1 is a positron 
      CONST1=-2.d0*Wgt0*DBLE(ITYP)

c Average beta in track
      E_AV=(E1+E2)/2.
      G=E_AV/e_M
      G2=G*G
      BETA=SQRT(1.-1./G2)
c BETA is parallel to the vector (r2-r1)
c Compute unitary vector parallel to (r2-r1)
      R2_R1=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
      IF (R2_R1.EQ.0) THEN
        NZEROS=NZEROS+1
        RETURN
      END IF 
      BETA_X=BETA*(X2-X1)/R2_R1
      BETA_Y=BETA*(Y2-Y1)/R2_R1
      BETA_Z=BETA*(Z2-Z1)/R2_R1

c Assume azimuth observation angle phi=0

c FACTOR_T = 1/(c*rho)
      DO NA=1,NAMAX    ! Loop in antenna positions 

c Vector k from starting point of track to observer 
        UX = XANT(NA)-(X1+X2)/2.     
        UY = YANT(NA)-(Y1+Y2)/2.
        UZ = ZANT(NA)-(Z1+Z2)/2.       ! g cm^-2

        UX1 = XANT(NA)-X1     
        UY1 = YANT(NA)-Y1
        UZ1 = ZANT(NA)-Z1

        UX2 = XANT(NA)-X2     
        UY2 = YANT(NA)-Y2
        UZ2 = ZANT(NA)-Z2

c Antenna position w.r.t. injection point of the shower        
        RANT = SQRT(XANT(NA)*XANT(NA)+
     #              YANT(NA)*YANT(NA)+
     #              ZANT(NA)*ZANT(NA))  ! g cm^-2 

c        write(*,*) 't ',na,namax,xant(na),yant(na),zant(na)

c Distance from medium point of track to observer 
        R = SQRT(UX*UX+UY*UY+UZ*UZ)  ! g cm^-2
        R1 = SQRT(UX1*UX1+UY1*UY1+UZ1*UZ1)  ! g cm^-2
        R2 = SQRT(UX2*UX2+UY2*UY2+UZ2*UZ2)  ! g cm^-2

c Unit vector k  
        UX=UX/R
        UY=UY/R
        UZ=UZ/R

c Unit vector k1  
        UX1=UX1/R1
        UY1=UY1/R1
        UZ1=UZ1/R1

c Unit vector k2  
        UX2=UX2/R2
        UY2=UY2/R2
        UZ2=UZ2/R2

        AUXB = REFIDX*UX*BETA_X 
     #       + REFIDX*UY*BETA_Y 
     #       + REFIDX*UZ*BETA_Z

        AUX1 = REFIDX*UX1*BETA_X 
     #       + REFIDX*UY1*BETA_Y 
     #       + REFIDX*UZ1*BETA_Z

c        DENOM = 1.d0 - AUX1
        DENOM = 1.d0 - AUXB
       

c Times at the observer:
c Arbitrary origin of times = time at which the signal from  
c the injection point in the shower reaches the corresponding antenna        
        CT1_OBS = CT1 + REFIDX*(R1-RANT)      ! at the observer [g cm^-2] 
        DELTAT_1 = FACTOR_T*CT1_OBS           ! at the observer [s]

        CT2_OBS = CT2 + REFIDX*(R1-RANT)-AUXB*(CT2-CT1)      ! at the observer [g cm^-2] 
        DELTAT_2 = FACTOR_T*CT2_OBS           ! at the observer [s]

        DTT = FACTOR_T*(CT2-CT1)              ! s

c k x (k x beta) 
        BETA_PERP_X=-(UY*UY+UZ*UZ)*BETA_X + UX*UY*BETA_Y + UX*UZ*BETA_Z
        BETA_PERP_Y=UX*UY*BETA_X - (UX*UX+UZ*UZ)*BETA_Y + UY*UZ*BETA_Z
        BETA_PERP_Z=UX*UZ*BETA_X + UY*UZ*BETA_Y - (UX*UX+UY*UY)*BETA_Z

        BETA_PERP=SQRT(BETA_PERP_X*BETA_PERP_X +
     #                 BETA_PERP_Y*BETA_PERP_Y +
     #                 BETA_PERP_Z*BETA_PERP_Z)

c WARNING
c     k x (k x beta) is actually -beta_perp and NOT beta_perp => change signs
        BETA_PERP_X=-BETA_PERP_X
        BETA_PERP_Y=-BETA_PERP_Y
        BETA_PERP_Z=-BETA_PERP_Z


c Find time bins in which the contribution to vector potential is non-zero
        IF (DELTAT_1.LT.0.) THEN
          IT1=MAXTHALF + INT(DELTAT_1/DT_BIN_S)-1
        ELSE
          IT1=MAXTHALF + INT(DELTAT_1/DT_BIN_S)
        END IF

        IF (DELTAT_2.LT.0.) THEN
          IT2=MAXTHALF + INT(DELTAT_2/DT_BIN_S)-1
        ELSE
          IT2=MAXTHALF + INT(DELTAT_2/DT_BIN_S)
        END IF
c        if (it1.ne.it2) write(*,*) DELTAT_1*1.d9,DELTAT_2*1.d9,IT1,IT2

        IF (IT1.LE.IT2) THEN 
           IT_STA=IT1
           IT_END=IT2
        ELSE
           IT_STA=IT2
           IT_END=IT1
        END IF   

c Get rid of times outside time range
        IF (IT_STA.GT.MAXT.OR.IT_END.LT.1) GOTO 200  ! Do not compute field
        IF (IT_STA.LT.1) IT_STA=1
        IF (IT_END.GT.MAXT) IT_END=MAXT

        CONST2=CONST1
c Change of sign when DELTAT_2 < DELTAT_1
        IF (DELTAT_2.LT.DELTAT_1) THEN
           CONST2=-CONST1
c        write(*,*) DELTAT_2*1.d9-DELTAT_1*1.d9,IT1,IT2,ITYP,CONST2 
        END IF 

c Note ASUM_i is multiplied by e*mu_r/(8.*pi*epsilon_0*c)
c when writing it to the file. One factor 1/c goes into beta 

c Distribute vector potential among bins in time
c and calculate average vector potential in each bin        
c (see my notes: MC implementation)
c Since we take the derivative at the end the important 
c point is to conserve the structure from bin to bin.
        IF (DELTAT_2.LT.DELTAT_1) THEN 
            DELTAT_STA=DELTAT_2
            DELTAT_END=DELTAT_1
        ELSE
            DELTAT_STA=DELTAT_1
            DELTAT_END=DELTAT_2
        END IF

           CTHETA=-UX*ETA_PERP_XANT(NA)
     #            -UY*ETA_PERP_YANT(NA)
     #            -UZ*ETA_PERP_ZANT(NA)

        DO IT=IT_STA,IT_END,1
        
          IF(IT_STA.EQ.IT_END) THEN   ! Subtrack contained in bin

           if (it_sta.gt.maxt.or.it_end.gt.maxt) write(*,*)
     #     'IT_STA>MAXT or IT_END<MAXT ',      
     #       j,it_sta,it_end,delta_sta,delta_end

          F=(DELTAT_END-DELTAT_STA)/DT_BIN_S

            IF (ABS(DENOM).GT.1.D-15) THEN
              ASUM_X_FR(IT,NA,ith) = ASUM_X_FR(IT,NA,ith) + 
     #                         ABS(F)*CONST2*BETA_PERP_X/DENOM/R
              ASUM_Y_FR(IT,NA,ith) = ASUM_Y_FR(IT,NA,ith) + 
     #                         ABS(F)*CONST2*BETA_PERP_Y/DENOM/R
              ASUM_Z_FR(IT,NA,ith) = ASUM_Z_FR(IT,NA,ith) + 
     #                         ABS(F)*CONST2*BETA_PERP_Z/DENOM/R

c -----------------------------------------------------------------
c Write z-component of field contributed by each sub-track
          if (ith.eq.-1) then  ! Uncomment this line and comment next to deactivate writing
c          if (ith.eq.1) then
          t_arr=((CT1+REFIDX*R1+CT2+REFIDX*R2)/2./rho_air/100./3.e8)  ! Arrival time at antenna (s)
          t_arr = t_arr*1.e9 ! in ns
c          write(*,*) ct1,refidx*r1,ct2,refidx*r2,rho_air,Tarr
          tck=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
          tmp=ABS(F)*CONST2*BETA_PERP_Z/DENOM/R
          Ez_diff=tmp/tck  ! "Differential" field - Normalized to track length
          Ez_acc=ASUM_Z_FR(IT,NA,ith)    ! "Accumulated" field 
c Assume isotropic emission (i.e. remove denom)
          tmp=ABS(F)*CONST2*BETA_PERP_Z/R
          Ez_iso=tmp/tck ! Normalized to track length
          write(222,*) 
     #    (x1+x2)/2./rho_air/100.,  ! 1
     #    (y1+y2)/2./rho_air/100.,  ! 2
     #    (z1+z2)/2./rho_air/100.,  ! 3
     #    t_arr,                    ! 4
     #    R/rho_air/100.,           ! 5
     #    tck,                      ! 6
     #    denom,                    ! 7
     #    Ez_diff,                  ! 8
     #    Ez_acc,                   ! 9 
     #    Ez_iso                    ! 10
          end if  
c -----------------------------------------------------------------


c Field projected onto antenna direction and accounting for antenna pattern
c Note that pattern peaks in the direction perpendicular to the antenna orientation
c i.e. we multiply by cos(theta) with theta the angle between minus the direction 
c from the source to the antenna and the normal to the antenna.
              ASUM_ANT_FR(IT,NA,ith)=ASUM_ANT_FR(IT,NA,ith)+
     #       ((ABS(F)*CONST2*BETA_PERP_X/DENOM/R)*ETA_XANT(NA)+   
     #        (ABS(F)*CONST2*BETA_PERP_Y/DENOM/R)*ETA_YANT(NA)+   
     #        (ABS(F)*CONST2*BETA_PERP_Z/DENOM/R)*ETA_ZANT(NA))*CTHETA

c             write(*,*)
c     # DELTAT_1*1.d9,DELTAT_2*1.d9,IT,IT_STA,IT_END,ABS(F)

             ELSE   ! Cherenkov angle approximation
              ASUM_X_FR(IT,NA,ith) = ASUM_X_FR(IT,NA,ith) + 
     #                         ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_X/R
              ASUM_Y_FR(IT,NA,ith) = ASUM_Y_FR(IT,NA,ith) + 
     #                         ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_Y/R
              ASUM_Z_FR(IT,NA,ith) = ASUM_Z_FR(IT,NA,ith) + 
     #                         ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_Z/R

              ASUM_ANT_FR(IT,NA,ith)=ASUM_ANT_FR(IT,NA,ith)+
     #       ((ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_X/R)*ETA_XANT(NA)+   
     #        (ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_Y/R)*ETA_YANT(NA)+   
     #        (ABS(DTT/DT_BIN_S)*CONST2*BETA_PERP_Z/R)*ETA_ZANT(NA))*
     #         CTHETA

 
             END IF

          ELSE  ! Subtrack crossing one or more bins in time

c Correct for bin edges               

               F=((IT_STA+1-MAXTHALF)*DT_BIN_S-DELTAT_STA)/DT_BIN_S
c (Bug corrected March 2010)
c               IF (ABS(DENOM).GT.1.D-15) THEN  
c We should also account for ABS(DENOM)<1.D-15 in this case               
             IF (IT.EQ.IT_STA) THEN ! Start of subtrack 
               
                 ASUM_X_FR(IT,NA,ith) = ASUM_X_FR(IT,NA,ith) + 
     #                          ABS(F)*CONST2*BETA_PERP_X/DENOM/R
                 ASUM_Y_FR(IT,NA,ith) = ASUM_Y_FR(IT,NA,ith) + 
     #                          ABS(F)*CONST2*BETA_PERP_Y/DENOM/R
                 ASUM_Z_FR(IT,NA,ith) = ASUM_Z_FR(IT,NA,ith) + 
     #                          ABS(F)*CONST2*BETA_PERP_Z/DENOM/R

              ASUM_ANT_FR(IT,NA,ith)=ASUM_ANT_FR(IT,NA,ith)+
     #       ((ABS(F)*CONST2*BETA_PERP_X/DENOM/R)*ETA_XANT(NA)+   
     #        (ABS(F)*CONST2*BETA_PERP_Y/DENOM/R)*ETA_YANT(NA)+   
     #        (ABS(F)*CONST2*BETA_PERP_Z/DENOM/R)*ETA_ZANT(NA))*CTHETA

c             write(*,*)
c     # DELTAT_STA*1.d9,DELTAT_END*1.d9,IT,IT_STA,IT_END,ABS(F)

             ELSE IF (IT.EQ.IT_END) THEN  ! End of subtrack

               F=((IT_END+1-MAXTHALF)*DT_BIN_S-DELTAT_END)/DT_BIN_S
               ASUM_X_FR(IT,NA,ith) = ASUM_X_FR(IT,NA,ith) + 
     #                        (1.d0-ABS(F))*CONST2*BETA_PERP_X/DENOM/R
               ASUM_Y_FR(IT,NA,ith) = ASUM_Y_FR(IT,NA,ith) + 
     #                        (1.d0-ABS(F))*CONST2*BETA_PERP_Y/DENOM/R
               ASUM_Z_FR(IT,NA,ith) = ASUM_Z_FR(IT,NA,ith) + 
     #                        (1.d0-ABS(F))*CONST2*BETA_PERP_Z/DENOM/R

              ASUM_ANT_FR(IT,NA,ith)=ASUM_ANT_FR(IT,NA,ith)+
     #       (((1.d0-ABS(F))*CONST2*BETA_PERP_X/DENOM/R)*ETA_XANT(NA)+  
     #        ((1.d0-ABS(F))*CONST2*BETA_PERP_Y/DENOM/R)*ETA_YANT(NA)+  
     #        ((1.d0-ABS(F))*CONST2*BETA_PERP_Z/DENOM/R)*ETA_ZANT(NA))*
     #        CTHETA

c             write(*,*)
c     # DELTAT_STA*1.d9,DELTAT_END*1.d9,IT,IT_STA,IT_END,(1.-ABS(F))

             ELSE

               F=1.d0
               ASUM_X_FR(IT,NA,ith) = ASUM_X_FR(IT,NA,ith) +
     #                          F*CONST2*BETA_PERP_X/DENOM/R
               ASUM_Y_FR(IT,NA,ith) = ASUM_Y_FR(IT,NA,ith) +
     #                          F*CONST2*BETA_PERP_Y/DENOM/R
               ASUM_Z_FR(IT,NA,ith) = ASUM_Z_FR(IT,NA,ith) +
     #                          F*CONST2*BETA_PERP_Z/DENOM/R

              ASUM_ANT_FR(IT,NA,ith)=ASUM_ANT_FR(IT,NA,ith)+
     #       ((F*CONST2*BETA_PERP_X/DENOM/R)*ETA_XANT(NA)+   
     #        (F*CONST2*BETA_PERP_Y/DENOM/R)*ETA_YANT(NA)+   
     #        (F*CONST2*BETA_PERP_Z/DENOM/R)*ETA_ZANT(NA))*CTHETA

             END IF   

          END IF

          END DO    ! End loop in time

200    CONTINUE

c CHECKING GEOMETRY      
c z position of track
c ctheta=cos Angle between direction to antenna and normal to antenna
c cwtheta=cos Angle between track and wire antenna orientation
c Projections of beta_perp onto wire normalized to beta_perp
       icheck_geom=0
       if(icheck_geom.eq.1) then
       cwtheta=ux*eta_xant(na)
     #        +uy*eta_yant(na)
     #        +uz*eta_zant(na)
       write(98,*) z1,acos(ctheta)*180./pi,acos(cwtheta)*180./pi,
     #(beta_perp_x*eta_xant(na)+
     # beta_perp_y*eta_yant(na)+
     # beta_perp_z*eta_zant(na))/beta_perp
       end if

       END DO      ! End loop in antenna positions 

      RETURN
      END
