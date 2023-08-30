      subroutine pes(x,igrad,path,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)
      COMMON /PT1CM/ Q(3),XTER,DQ(3)
      COMMON /PT2CM/ NDER, NDUM, IDEBUG, IDUM(19)
      double precision :: tx(9), r2(3), dedx(9), drdx(3,9)
      logical, save :: first_time_data=.true.

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        tx(j)=x(iatom,idir)
      enddo
      enddo
      ! input cartesian is ClHH
      Q(1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      Q(2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      Q(3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211
      
      NDER=igrad
      DQ=0.d0

      if(first_time_data) then
      call prepot(path)
      first_time_data=.false.
      endif

      call pot

      XTER=XTER*27.211386
      DQ=DQ*51.422067

      r2=q*0.529177211
      call evdrdx(tx, r2, drdx)
      dedx=0.d0
      do i=1,9
      do j=1,3
        dedx(i)=dedx(i)+DQ(j)*drdx(j,i)
      enddo
      enddo


      if (igrad==0) then
        do istate=1,nstates
          p(istate)=XTER
        enddo
      else if (igrad==1) then
        do istate=1,nstates
          p(istate)=XTER
        enddo
        do iatom=1,natoms
        do idir=1,3
          j=(iatom-1)*3+idir
          g(1,iatom,idir)=dedx(j)
        enddo
        enddo
      else if (igrad==2) then
        write (*,*) 'Only energy and gradient are available'
      endif

      endsubroutine

      subroutine EvdRdX(X,r,drdx)

      integer i,j
      double precision, intent(in) :: X(9), R(3)
      double precision, intent(out) :: dRdX(3,9)

! Initialize dRdX(3,9)
      do i=1,3
        do j=1,9
          dRdX(i,j)=0.0d0
        enddo
      enddo

      dRdX(1,1)=(x(1)-x(4))/r(1)
      dRdX(1,2)=(x(2)-x(5))/r(1)
      dRdX(1,3)=(x(3)-x(6))/r(1)
      dRdX(1,4)=-dRdX(1,1)
      dRdX(1,5)=-dRdX(1,2)
      dRdX(1,6)=-dRdX(1,3)

      dRdX(2,4)=(x(4)-x(7))/r(2)
      dRdX(2,5)=(x(5)-x(8))/r(2)
      dRdX(2,6)=(x(6)-x(9))/r(2)
      dRdX(2,7)=-dRdX(2,4)
      dRdX(2,8)=-dRdX(2,5)
      dRdX(2,9)=-dRdX(2,6)

      dRdX(3,1)=(x(1)-x(7))/r(3)
      dRdX(3,2)=(x(2)-x(8))/r(3)
      dRdX(3,3)=(x(3)-x(9))/r(3)
      dRdX(3,7)=-dRdX(3,1)
      dRdX(3,8)=-dRdX(3,2)
      dRdX(3,9)=-dRdX(3,3)

      endsubroutine

C
      SUBROUTINE PREPOT(PATH)
C
C   System:          ClH2
C   Functional form: Isaacson and Muckerman's diatomics-in-molecules version S
C   Common name:     ClH2 DIMS
C   Reference:       unpublished
C   Cross Reference: S. C. Tucker, D. G. Truhlar,
C                    B. C. Garrett, and A. D. Isaacson
C                    J. Chem. Phys. 82, 4102 (1985)
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy and the derivatives 
C               of the energy with respect to the coordinates for a given 
C               geometric configuration
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohr
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground electronic state
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Cl
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Parameters:
C      Set in the BLOCK DATA subprogram PTPARM and read in from the file
C      potclh2dims.dat that is linked to FORTRAN unit IPTRD1
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Cl-H)
C                            R(2) = R(H-H)
C                            R(3) = R(Cl-H)
C
C   Common Blocks (used between the calling program and this potential):
C      /PT1CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, ground state electronic energy, and 
C        derivatives of the ground electronic state energy with respect 
C        to the coordinates.
C      /PT2CM/ NDER, NDUM, IDEBUG, IDUM(19)
C        passes the control flags where
C        NDER   = 0 => no derivatives are computed
C               = 1 => derivatives of the energy for the ground electronic 
C                     state with respect to the coordinates are computed
C        NDUM   - not used 
C        IDEBUG - control flag for priniting extra information
C        IDEBUG = 0 => do not print extra debug information
C        IDEBUG = 1 => print intermediate values to the file linked 
C                      to FORTRAN unit IPTPRT.
C        IDUM   - not used 
C      /PT4CM/ IPTPRT, IPTRD1
C        passes the FORTRAN unit number used for potential output
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the 
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to 
C        Re(AB), the equilibrium bond length for the AB diatomic.  
C        In this potential the AB valley represents H infinitely far from
C        the ClH diatomic and R(ClH) equal to Re(ClH).  Similarly, the terms
C        EASYBC and EASYAC represent the energies in the H2 and the other ClH 
C        valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NDER             1 
C      IDEBUG           0
C      IPTPRT           6
C      IPTRD1           4
C
C*****
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character(len=1024), intent(in) :: path
      COMMON /PT4CM/ IPTPRT, IPTRD1
      COMMON /PT2CM/ NDER, NDUM, IDEBUG, IDUM(19)
      COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
      DIMENSION TH(6,6),VEC(6,6),W(6),EN(6),GVEC(6),DH(6,6),R(3),DD(3)
      DIMENSION DH11(3),DH12(3),DH13(3),DH14(3),DH15(3),DH16(3),DH22(3)
      DIMENSION DH24(3),DH25(3),DH26(3),DH33(3),DH34(3),DH44(3),DH45(3)
      DIMENSION DH46(3),DH55(3),DH56(3),DH66(3)
      DIMENSION INDX(3)
      CHARACTER*80 TITLE(2)
      COMMON /PT1CM/ Q(3),XTER,DQ(3)
      COMMON /POTCM2/ C(3)
      COMMON/POTCM3/ZERO,ONE,TWO,TRE,HALF,FORTH,SQRT3,
     $DPUX,SSGX,DSGH,SSGH
      COMMON/POTCM4/X(121),
     $SABFIT(121),S2SAB(121),S3SAB(121),DELSAB(121),
     $TABFIT(121),S2TAB(121),S3TAB(121),DELTAB(121),
     $GABFIT(121),S2GAB(121),S3GAB(121),DELGAB(121),
     $UABFIT(121),S2UAB(121),S3UAB(121),DELUAB(121),
     $S11FIT(121),S2S11(121),S3S11(121),DELS11(121),
     $S12FIT(121),S2S12(121),S3S12(121),DELS12(121),
     $S22FIT(121),S2S22(121),S3S22(121),DELS22(121),
     $TSXFIT(121),S2TSX(121),S3TSX(121),DELTSX(121),
     $SPXFIT(121),S2SPX(121),S3SPX(121),DELSPX(121),
     $TPXFIT(121),S2TPX(121),S3TPX(121),DELTPX(121),
     $DSXFIT(121),S2DSX(121),S3DSX(121),DELDSX(121),NS
      EQUIVALENCE (TH(1,1),DH(1,1))
      character(len=1024) :: file_path1

      DATA NDIM/6/,NROOT/1/
      DATA TEST/0.999999999D0/
      DATA INDX /1, 2, 3/
      DATA TITLE(1) /'    Final BNL DIM surface for Cl+H2 -- 5 and 7 kca
     *l/mole barriers'/
      DATA TITLE(2) /'    Spline fits used in atomic units'/
C
C   Echo the name of the potential to the file linked to FORTRAN unit IPTPRT
      !WRITE (IPTPRT, 600)
      !WRITE (IPTPRT, 610) TITLE

      file_path1 = trim(path)//"/ClH2/potclh2dims.dat "
C   Open the file with the potential input data
      OPEN (UNIT=IPTRD1, FILE=file_path1, FORM='FORMATTED', 
     *      STATUS='OLD', ERR=100)
C     SET UP SPLINE FITS TO H2 CURVES
      CALL SETSPL
C
C   Close the input data file
      CLOSE (UNIT = IPTRD1)
C
C   Initialize the energies in the three asymptotic valleys
      EASYAB = C(1)
      EASYBC = C(2)
      EASYAC = C(3)
C
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'ClH2 DIMS potential energy surface')
610   FORMAT(/,2X,T5,'Title cards for the potential routine:',
     *       /,80A,/,80A,//,1X,'*****')
6000  FORMAT(/,2X,T5,'Error opening input data file')
      RETURN
100   WRITE (IPTPRT, 6000)
      STOP 'PREPOT 1'
C
      ENTRY POT
C
C***********************************************************************
C
C     ENTRY POT GIVES NROOT-TH POTENTIAL ENERGY FOR THE CLH2 SYSTEM
C     ENTRY POT ALSO GIVES PARTIALS W.R.T. THREE DISTANCES
C
C     INTERNUCLEAR DISTANCES STORED IN Q(3) AND R(3)
C     NROOT-TH STATE ENERGY STORED IN XTER (OR TERM FOR DIATOMIC)
C     THREE PARTIALS STORED IN DD(3)
C
C***********************************************************************
C
C   Check the value of NDER 
         IF (NDER .GT. 1) THEN
             WRITE (IPTPRT, 900) NDER
             STOP 'POT 1'
         ENDIF
C
C     IAPP=1
C     SET FLAG FOR EIGENVECTOR OUTPUT
      IFLAG=SIGN(1,NROOT)
      NROOT=ABS(NROOT)
C
C     FOLLOWING COORDINATE DEFINITIONS PROVIDE FOR PROPER SYSTEM
      DO 15 IX=1,3
   15 R(IX)=Q(INDX(IX))
C
C     ZERO OUT UNUSED HAMILTONIAN ELEMENTS
      TH(2,3)=ZERO
      TH(3,5)=ZERO
      TH(3,6)=ZERO
C
C     COMPUTE SINE AND COSINE OF ANGLE AXB, AND POWERS
      CT=(R(1)*R(1)+R(3)*R(3)-R(2)*R(2))/(TWO*R(1)*R(3))
      IF(CT.GT.-TEST .AND. CT.LT.TEST) GO TO 10                         07AUG83
      CT = SIGN(TEST,CT)                                                07AUG83
   10 CONTINUE
      ARG=ONE-CT*CT
      IF(ARG.LT.ZERO) ARG=ZERO
      ST=SQRT(ARG)
      CT2=CT*CT
      ST2=ST*ST
      CS=CT*ST
C     COMPUTE DIATOMIC ENERGIES FOR R(3) ARRAY
      CALL SPLINE(NS,X,S11FIT,R(1),SS11A,DSS11A,S3,EPS,S2S11,S3S11,
     $DELS11)
      CALL SPLINE(NS,X,S11FIT,R(3),SS11B,DSS11B,S3,EPS,S2S11,S3S11,
     $DELS11)
      CALL SPLINE(NS,X,S12FIT,R(1),SS12A,DSS12A,S3,EPS,S2S12,S3S12,
     $DELS12)
      CALL SPLINE(NS,X,S12FIT,R(3),SS12B,DSS12B,S3,EPS,S2S12,S3S12,
     $DELS12)
      CALL SPLINE(NS,X,S22FIT,R(1),SS22A,DSS22A,S3,EPS,S2S22,S3S22,
     $DELS22)
      CALL SPLINE(NS,X,S22FIT,R(3),SS22B,DSS22B,S3,EPS,S2S22,S3S22,
     $DELS22)
      CALL SPLINE(NS,X,TSXFIT,R(1),TSXA,DTSXA,S3,EPS,S2TSX,S3TSX,DELTSX)
      CALL SPLINE(NS,X,TSXFIT,R(3),TSXB,DTSXB,S3,EPS,S2TSX,S3TSX,DELTSX)
      CALL SPLINE(NS,X,SPXFIT,R(1),SPXA,DSPXA,S3,EPS,S2SPX,S3SPX,DELSPX)
      CALL SPLINE(NS,X,SPXFIT,R(3),SPXB,DSPXB,S3,EPS,S2SPX,S3SPX,DELSPX)
      CALL SPLINE(NS,X,TPXFIT,R(1),TPXA,DTPXA,S3,EPS,S2TPX,S3TPX,DELTPX)
      CALL SPLINE(NS,X,TPXFIT,R(3),TPXB,DTPXB,S3,EPS,S2TPX,S3TPX,DELTPX)
      CALL SPLINE(NS,X,DSXFIT,R(1),DSXA,DDSXA,S3,EPS,S2DSX,S3DSX,DELDSX)
      CALL SPLINE(NS,X,DSXFIT,R(3),DSXB,DDSXB,S3,EPS,S2DSX,S3DSX,DELDSX)
      CALL SPLINE(NS,X,SABFIT,R(2),SSAB,DSSAB,S3,EPS,S2SAB,S3SAB,DELSAB)
      CALL SPLINE(NS,X,TABFIT,R(2),TSAB,DTSAB,S3,EPS,S2TAB,S3TAB,DELTAB)
      CALL SPLINE(NS,X,GABFIT,R(2),DGAB,DDGAB,S3,EPS,S2GAB,S3GAB,DELGAB)
      CALL SPLINE(NS,X,UABFIT,R(2),DUAB,DDUAB,S3,EPS,S2UAB,S3UAB,DELUAB)
C
C     CONSTRUCT DIM HAMILTONIAN IN ORTHOGONAL BASIS
C
      TH(1,1)=SS11A+FORTH*((SSAB+TRE*TSAB)+CT2*(SS11B+TRE*TSXB)
     1 +ST2*(SPXB+TRE*TPXB))-DPUX-TWO*DSGH
      TH(1,2)=FORTH*SQRT3*(SSAB-TSAB-CT2*(SS11B-TSXB)-ST2*(SPXB-TPXB))
      TH(1,3)=SS12A
      TH(1,4)=-HALF*CT*SS12B
      TH(1,5)=FORTH*CS*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)
      TH(1,6)=-FORTH*SQRT3*CS*(SS11B-TSXB-SPXB+TPXB)
      TH(2,2)=TSXA+FORTH*((TRE*SSAB+TSAB)+CT2*(TRE*SS11B+TSXB)
     1 +ST2*(TRE*SPXB+TPXB))-DPUX-TWO*DSGH
      TH(2,4)=HALF*SQRT3*CT*SS12B
      TH(2,5)=TH(1,6)
      TH(2,6)=FORTH*CS*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)
      TH(3,3)=SS22A+HALF*(DGAB+DUAB)+DSXB-SSGX-SSGH-DSGH
      TH(3,4)=HALF*(DGAB-DUAB)
      TH(4,4)=DSXA+HALF*(DGAB+DUAB)+SS22B-SSGX-SSGH-DSGH
      TH(4,5)=-HALF*ST*SS12B
      TH(4,6)=HALF*SQRT3*ST*SS12B
      TH(5,5)=SPXA+FORTH*((SSAB+TRE*TSAB)+ST2*(SS11B+TRE*TSXB)
     1 +CT2*(SPXB+TRE*TPXB))-DPUX-TWO*DSGH
      TH(5,6)=FORTH*SQRT3*(SSAB-TSAB-ST2*(SS11B-TSXB)-CT2*(SPXB-TPXB))
      TH(6,6)=TPXA+FORTH*((TRE*SSAB+TSAB)+ST2*(TRE*SS11B+TSXB)
     1 +CT2*(TRE*SPXB+TPXB))-DPUX-TWO*DSGH
C
C     FILL OUT MATRIX AND DIAGONALIZE
C
      DO 30 I=1,NDIM
      DO 30 J=1,I
      TH(I,J)=TH(J,I)
   30 CONTINUE
C
      CALL EIGN(NDIM,TH,VEC,EN,W,NDIM)
C     FOLLOWING STATEMENT ASSUMES ENERGY ZERO AT A + BC(R=RE)
      XTER=EN(NROOT)+C(2)
      IF(IFLAG.LT.ZERO) WRITE(IPTPRT,998) (VEC(I,NROOT),I=1,NDIM)
998   FORMAT(/,2X,T5,'Vector:',(T15,3E20.8))
C
C************************* END OF FPOT3N *******************************
C
      IF (NDER .NE. 1) RETURN
C     IF(IAPP.EQ.0) RETURN
C
C     SET  VECTOR FOR STATE OF INTEREST
      DO 40 I=1,NDIM
      GVEC(I)=VEC(I,NROOT)
   40 CONTINUE
C     COMPUTE PARTIALS OF U=COS(THETA)
      DTDR1=(R(1)*R(1)+R(2)*R(2)-R(3)*R(3))/(TWO*R(1)*R(1)*R(3))
      DTDR2=-R(2)/(R(1)*R(3))
      DTDR3=(R(2)*R(2)+R(3)*R(3)-R(1)*R(1))/(TWO*R(3)*R(3)*R(1))
      FACTR=-ONE/ST
      IF(ABS(ST).LT.1.D-6) WRITE(IPTPRT,905) ST
905   FORMAT(/,2X,T5,'Warning: A serious error may result because ',
     *               'sine(THETA) = ',1PE20.10)
C     COMPUTE PARTIALS OF HAMILTONIAN ELEMENTS
      XTRA=+HALF*CT*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)
      DH11(1)=DSS11A+XTRA*DTDR1
      DH11(2)=FORTH*(DSSAB+TRE*DTSAB)+XTRA*DTDR2
      DH11(3)=FORTH*(CT2*(DSS11B+TRE*DTSXB)+ST2*(DSPXB+TRE*DTPXB))+XTRA*
     1 DTDR3
      XTRA=-HALF*SQRT3*CT*(SS11B-TSXB-SPXB+TPXB)
      DH12(1)=XTRA*DTDR1
      DH12(2)=FORTH*SQRT3*(DSSAB-DTSAB)+XTRA*DTDR2
      DH12(3)=-FORTH*SQRT3*(CT2*(DSS11B-DTSXB)+ST2*(DSPXB-DTPXB))+XTRA*
     1 DTDR3
      DH13(1)=DSS12A
      DH13(2)=ZERO
      DH13(3)=ZERO
      XTRA=-HALF*SS12B
      DH14(1)=XTRA*DTDR1
      DH14(2)=XTRA*DTDR2
      DH14(3)=-HALF*CT*DSS12B+XTRA*DTDR3
      XTRA=FORTH*(CT2-ST2)*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)*FACTR
      DH15(1)=XTRA*DTDR1
      DH15(2)=XTRA*DTDR2
      DH15(3)=FORTH*CS*(DSS11B+TRE*DTSXB-DSPXB-TRE*DTPXB)+XTRA*DTDR3
      XTRA=-FORTH*SQRT3*(CT2-ST2)*(SS11B-TSXB-SPXB+TPXB)*FACTR
      DH16(1)=XTRA*DTDR1
      DH16(2)=XTRA*DTDR2
      DH16(3)=-FORTH*SQRT3*CS*(DSS11B-DTSXB-DSPXB+DTPXB)+XTRA*DTDR3
      XTRA=+HALF*CT*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)
      DH22(1)=DTSXA+XTRA*DTDR1
      DH22(2)=FORTH*(TRE*DSSAB+DTSAB)+XTRA*DTDR2
      DH22(3)=FORTH*(CT2*(TRE*DSS11B+DTSXB)+ST2*(TRE*DSPXB+DTPXB))+XTRA*
     1 DTDR3
      XTRA=+HALF*SQRT3*SS12B
      DH24(1)=XTRA*DTDR1
      DH24(2)=XTRA*DTDR2
      DH24(3)=HALF*SQRT3*CT*DSS12B+XTRA*DTDR3
      DH25(1)=DH16(1)
      DH25(2)=DH16(2)
      DH25(3)=DH16(3)
      XTRA=FORTH*(CT2-ST2)*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)*FACTR
      DH26(1)=XTRA*DTDR1
      DH26(2)=XTRA*DTDR2
      DH26(3)=FORTH*CS*(TRE*DSS11B+DTSXB-TRE*DSPXB-DTPXB)+XTRA*DTDR3
      DH33(1)=DSS22A
      DH33(2)=HALF*(DDGAB+DDUAB)
      DH33(3)=DDSXB
      DH34(1)=ZERO
      DH34(2)=HALF*(DDGAB-DDUAB)
      DH34(3)=ZERO
      DH44(1)=DDSXA
      DH44(2)=HALF*(DDGAB+DDUAB)
      DH44(3)=DSS22B
      XTRA=-HALF*CT*SS12B*FACTR
      DH45(1)=XTRA*DTDR1
      DH45(2)=XTRA*DTDR2
      DH45(3)=-HALF*ST*DSS12B+XTRA*DTDR3
      XTRA=HALF*SQRT3*CT*SS12B*FACTR
      DH46(1)=XTRA*DTDR1
      DH46(2)=XTRA*DTDR2
      DH46(3)=HALF*SQRT3*ST*DSS12B+XTRA*DTDR3
      XTRA=-HALF*CT*(SS11B+TRE*TSXB-SPXB-TRE*TPXB)
      DH55(1)=DSPXA+XTRA*DTDR1
      DH55(2)=FORTH*(DSSAB+TRE*DTSAB)+XTRA*DTDR2
      DH55(3)=FORTH*(ST2*(DSS11B+TRE*DTSXB)+CT2*(DSPXB+TRE*DTPXB))+XTRA*
     1 DTDR3
      XTRA=+HALF*SQRT3*CT*(SS11B-TSXB-SPXB+TPXB)
      DH56(1)=XTRA*DTDR1
      DH56(2)=FORTH*SQRT3*(DSSAB-DTSAB)+XTRA*DTDR2
      DH56(3)=-FORTH*SQRT3*(ST2*(DSS11B-DTSXB)+CT2*(DSPXB-DTPXB))+XTRA*
     1 DTDR3
      XTRA=-HALF*CT*(TRE*SS11B+TSXB-TRE*SPXB-TPXB)
      DH66(1)=DTPXA+XTRA*DTDR1
      DH66(2)=FORTH*(TRE*DSSAB+DTSAB)+XTRA*DTDR2
      DH66(3)=FORTH*(ST2*(TRE*DSS11B+DTSXB)+CT2*(TRE*DSPXB+DTPXB))+XTRA*
     1 DTDR3
C
C     LOOP OVER R(I) DERIVATIVES
C     FIRST ZERO OUT UNUSED ELEMENTS
      DH(2,3)=ZERO
      DH(3,5)=ZERO
      DH(3,6)=ZERO
      DO 50 I=1,3
C
C     SET UP DERIVATIVE MATRIX
C
      DH(1,1)=DH11(I)
      DH(1,2)=DH12(I)
      DH(1,3)=DH13(I)
      DH(1,4)=DH14(I)
      DH(1,5)=DH15(I)
      DH(1,6)=DH16(I)
      DH(2,2)=DH22(I)
      DH(2,4)=DH24(I)
      DH(2,5)=DH25(I)
      DH(2,6)=DH26(I)
      DH(3,3)=DH33(I)
      DH(3,4)=DH34(I)
      DH(4,4)=DH44(I)
      DH(4,5)=DH45(I)
      DH(4,6)=DH46(I)
      DH(5,5)=DH55(I)
      DH(5,6)=DH56(I)
      DH(6,6)=DH66(I)
C
C     FILL OUT DERIVATIVE MATRIX
      DO 70 J=1,NDIM
      DO 70 K=1,J
      DH(J,K)=DH(K,J)
   70 CONTINUE
C
      DNUMB=ZERO
      DO 80 J=1,NDIM
      DO 80 K=1,NDIM
      DNUMB=DNUMB+GVEC(J)*DH(J,K)*GVEC(K)
   80 CONTINUE
      DD(I)=DNUMB
C     END OF DERIVATIVE LOOP
   50 CONTINUE
C     PUT DERIVATIVES INTO CORRESPONDING SPOTS FOR CORRECT COORD CHOICE.
      DO 90 IX=1,3
   90 DQ(INDX(IX))=DD(IX)
C
C************************ END OF DFPOT *********************************
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
C
      RETURN
      END
      SUBROUTINE SETSPL
C     SETS UP SPLINE FITS FOR DIATOMIC CURVES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT4CM/ IPTPRT, IPTRD1
      COMMON/POTCM4/X(121),
     $SABFIT(121),S2SAB(121),S3SAB(121),DELSAB(121),
     $TABFIT(121),S2TAB(121),S3TAB(121),DELTAB(121),
     $GABFIT(121),S2GAB(121),S3GAB(121),DELGAB(121),
     $UABFIT(121),S2UAB(121),S3UAB(121),DELUAB(121),
     $S11FIT(121),S2S11(121),S3S11(121),DELS11(121),
     $S12FIT(121),S2S12(121),S3S12(121),DELS12(121),
     $S22FIT(121),S2S22(121),S3S22(121),DELS22(121),
     $TSXFIT(121),S2TSX(121),S3TSX(121),DELTSX(121),
     $SPXFIT(121),S2SPX(121),S3SPX(121),DELSPX(121),
     $TPXFIT(121),S2TPX(121),S3TPX(121),DELTPX(121),
     $DSXFIT(121),S2DSX(121),S3DSX(121),DELDSX(121),NS
      COMMON/POTCM3/ZERO,ONE,TWO,TRE,HALF,QUAR,RT3,
     $DPUX,SSGX,DSGH,SSGH
      DATA AUTERG/4.3598283D-11/,ANGTA0/1.88972634D0/
      NS=121
C     SET UP X ARRAY IN BOHR
      FAC=0.1D0*ANGTA0
      DO 20 I=1,NS
   20 X(I)=DBLE(I-1)*FAC
C     READ IN FIT DATA
  999 FORMAT(2D20.13)
      READ(IPTRD1,999) (SABFIT(I),TABFIT(I),I=1,NS)
      READ(IPTRD1,999) (GABFIT(I),UABFIT(I),I=1,NS)
      READ(IPTRD1,998) (S11FIT(I),S12FIT(I),S22FIT(I),TSXFIT(I),I=1,NS)
  998 FORMAT(4D20.13)
      READ(IPTRD1,997) (SPXFIT(I),TPXFIT(I),DSXFIT(I),I=1,NS)
  997 FORMAT(3D20.13)
C     CONVERT TO ATOMIC UNITS
      FAC=ONE/AUTERG
      DO 10 I=1,NS
      SABFIT(I)=SABFIT(I)*FAC
      TABFIT(I)=TABFIT(I)*FAC
      GABFIT(I)=GABFIT(I)*FAC
      UABFIT(I)=UABFIT(I)*FAC
      S11FIT(I)=S11FIT(I)*FAC
      S12FIT(I)=S12FIT(I)*FAC
      S22FIT(I)=S22FIT(I)*FAC
      TSXFIT(I)=TSXFIT(I)*FAC
      SPXFIT(I)=SPXFIT(I)*FAC
      TPXFIT(I)=TPXFIT(I)*FAC
      DSXFIT(I)=DSXFIT(I)*FAC
   10 CONTINUE
C     READ IN SECOND DERIVATIVES AT ENDS OF FIT CURVES
C     AND CONVERT TO ATOMIC UNITS
      FAC=FAC*1.D-16/(ANGTA0*ANGTA0)
      READ(IPTRD1,999) S2SAB(1),S2TAB(1),S2SAB(NS),S2TAB(NS)
      READ(IPTRD1,999) S2GAB(1),S2UAB(1),S2GAB(NS),S2UAB(NS)
      READ(IPTRD1,998) S2S11(1),S2S12(1),S2S22(1),S2TSX(1)
      READ(IPTRD1,998) S2S11(NS),S2S12(NS),S2S22(NS),S2TSX(NS)
      READ(IPTRD1,997) S2SPX(1),S2TPX(1),S2DSX(1)
      READ(IPTRD1,997) S2SPX(NS),S2TPX(NS),S2DSX(NS)
      S2SAB(1)=S2SAB(1)*FAC
      S2TAB(1)=S2TAB(1)*FAC
      S2GAB(1)=S2GAB(1)*FAC
      S2UAB(1)=S2UAB(1)*FAC
      S2S11(1)=S2S11(1)*FAC
      S2S12(1)=S2S12(1)*FAC
      S2S22(1)=S2S22(1)*FAC
      S2TSX(1)=S2TSX(1)*FAC
      S2SPX(1)=S2SPX(1)*FAC
      S2TPX(1)=S2TPX(1)*FAC
      S2DSX(1)=S2DSX(1)*FAC
C     GENERATE SPLINE FITS
      CALL SPLNIN(NS,X,SABFIT,R,R1,R2,R3,EPS,S2SAB,S3SAB,DELSAB)
      CALL SPLNIN(NS,X,TABFIT,R,R1,R2,R3,EPS,S2TAB,S3TAB,DELTAB)
      CALL SPLNIN(NS,X,GABFIT,R,R1,R2,R3,EPS,S2GAB,S3GAB,DELGAB)
      CALL SPLNIN(NS,X,UABFIT,R,R1,R2,R3,EPS,S2UAB,S3UAB,DELUAB)
      CALL SPLNIN(NS,X,S11FIT,R,R1,R2,R3,EPS,S2S11,S3S11,DELS11)
      CALL SPLNIN(NS,X,S12FIT,R,R1,R2,R3,EPS,S2S12,S3S12,DELS12)
      CALL SPLNIN(NS,X,S22FIT,R,R1,R2,R3,EPS,S2S22,S3S22,DELS22)
      CALL SPLNIN(NS,X,TSXFIT,R,R1,R2,R3,EPS,S2TSX,S3TSX,DELTSX)
      CALL SPLNIN(NS,X,SPXFIT,R,R1,R2,R3,EPS,S2SPX,S3SPX,DELSPX)
      CALL SPLNIN(NS,X,TPXFIT,R,R1,R2,R3,EPS,S2TPX,S3TPX,DELTPX)
      CALL SPLNIN(NS,X,DSXFIT,R,R1,R2,R3,EPS,S2DSX,S3DSX,DELDSX)
      RETURN
      END
      SUBROUTINE SPLINE(N,X,Y,T,SS,SS1,SS2,EPSLN,S2,S3,DELY)
C     JIM STINE'S SPLINE INTERPOLATOR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),S2(N),S3(N),DELY(N)
      I=1
      IF(T.GT.X(1)) GO TO 10
      SS = Y(1)
      SS1 = DELY(1)
      RETURN
   10 IF(T.LT.X(N)) GO TO 20
      SS = Y(N)
      SS1 = DELY(N)
      RETURN
   20 CONTINUE
   56 IF(T-X(I)) 60,17,57
   57 I=I+1
      GO TO 56
   59 I=N
   60 I=I-1
   17 HT1=T-X(I)
      HT2=T-X(I+1)
      PROD=HT1*HT2
      SS2=S2(I)+HT1*S3(I)
      DELSQS=(S2(I)+S2(I+1)+SS2)/6.D-00
      SS=Y(I)+HT1*DELY(I)+PROD*DELSQS
      SS1=DELY(I)+(HT1+HT2)*DELSQS+PROD*S3(I)/6.D-00
   61 RETURN
      END
      SUBROUTINE SPLNIN(N,X,Y,T,SS,SS1,SS2,EPSLN,S2,S3,DELY)
C     JIM STINE'S INITIALIZATION ROUTINE FOR SPLINE FITTING
C     N=NO. POINTS (X,Y); T=ARGUMENT; SS=VALUE; SS1=1ST DERIV;
C     SS2=2ND DERIV; S2,S3, AND DELY ARE COEFFICIENT ARRAYS
C     SECOND DERIV'S AT ENDS STORED IN S2(1) AND S2(N) ON INPUT
C     EPSLN IS OVERALL FIT CRITERION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT4CM/ IPTPRT, IPTRD1
      COMMON /PT2CM/ NDER, NDUM, IDEBUG, IDUM(19)
      DIMENSION C(121),H(121),H2(121),B(121)
      DIMENSION X(N),Y(N),S2(N),S3(N),DELY(N)
      DATA EX/1.1D0/,OMEGA/1.0717968D0/,HALF/0.5D0/
      EPSLN=6.4229D-13
      N1=N-1
      KNT=0
      DO 51 I=1,N1
      H(I)=X(I+1)-X(I)
   51 DELY(I)=(Y(I+1)-Y(I))/H(I)
      DO 52 I=2,N1
      H2(I)=H(I-1)+H(I)
      B(I)=HALF*H(I-1)/H2(I)
      DELSQY=(DELY(I)-DELY(I-1))/H2(I)
      S2(I)=DELSQY+DELSQY
   52 C(I)=DELSQY+S2(I)
    5 ETA=0.0D0
      DO 10 I=2,N1
      W=(C(I)-B(I)*S2(I-1)-(HALF-B(I))*S2(I+1)-S2(I))*OMEGA
      IF(ABS(W)-ETA) 10,10,9
    9 ETA=ABS(W)
   10 S2(I)=S2(I)+W
      KNT=KNT+1
      IF (KNT.GT.10) EPSLN=EPSLN*EX
      IF(ETA-EPSLN) 14,5,5
   14 DO 53 I=1,N1
   53 S3(I)=(S2(I+1)-S2(I))/H(I)
      IF (IDEBUG .EQ. 1) WRITE(IPTPRT,61) EPSLN
61    FORMAT(2X,T5,'Final EPS = ', 1PE15.6)
      RETURN
      END
      SUBROUTINE EIGN(NN,A,VEC,EIG,W,ND)
C     MATRIX DIAGONALIZATION ROUTINE FOR REAL SYMMETRIC CASE             VBDIM
C     HOUSEHOLDER METHOD                                                 VBDIM
C     RHO=UPPERLIMIT FOR OFF-DIAGONAL ELEMENT                            VBDIM
C     NN=SIZE OF MATRIX                                                  VBDIM
C     A=MATRIX (ONLY LOWER TRIANGLE IS USED ,THIS IS DESTROYED)          VBDIM
C     EIG=RETURNED EIGENVALUES IN ALGEBRAIC ASCENDING ORDER              VBDIM
C     VEC=RETURNED EIGENVECTORS IN COLUMNS                               VBDIM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EIG(ND),W(ND)                                            VBDIM
      DIMENSION GAMMA(6),BETA(6),BETASQ(6)                               VBDIM
      DIMENSION P(6),Q(6)                                                VBDIM
      EQUIVALENCE (P(1),BETA(1)),(Q(1),BETA(1))
      DIMENSION IPOSV(6),IVPOS(6),IORD(6)                                VBDIM
      DIMENSION A(ND,ND),VEC(ND,ND)                                      VBDIM
      RHO=1.0D-14                                                        ED1
      RHOSQ=RHO*RHO                                                      VBDIM
      N=NN                                                               VBDIM
      IF(N.EQ.0) GO TO 560                                               VBDIM
    1 N1=N-1                                                             VBDIM
      N2=N-2                                                             VBDIM
      GAMMA(1)=A(1,1)                                                    VBDIM
      IF(N2) 280,270,120                                                 VBDIM
  120 DO 260 NR=1,N2                                                     VBDIM
      B=A(NR+1,NR)                                                       VBDIM
      S=0.0D0                                                            VBDIM
      DO 130 I=NR,N2                                                     VBDIM
C     PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION                      VBDIM
  130 S=S+A(I+2,NR)**2                                                   VBDIM
      A(NR+1,NR)=0.0D0                                                   VBDIM
      IF(S) 250,250,140                                                  VBDIM
  140 S=S+B*B                                                            VBDIM
      SGN=+1.0D0                                                         VBDIM
      IF(B) 150,160,160                                                  VBDIM
  150 SGN=-1.0D0                                                         VBDIM
  160 SQRTS=SQRT(S)                                                      VBDIM
      D=SGN/(SQRTS+SQRTS)                                                VBDIM
      TEMP=SQRT(.5D0+B*D)                                                VBDIM
      W(NR)=TEMP                                                         VBDIM
      A(NR+1,NR)=TEMP                                                    VBDIM
      D=D/TEMP                                                           VBDIM
      B=-SGN*SQRTS                                                       VBDIM
C     D IS FACTOR OF PROPORTIONALITY  NOW COMPUTE AND SAVE W VECTOR      VBDIM
C     EXTRA SINGLY SUBSCRIPTED W VECTOR USED FOR SPEED                   VBDIM
      DO 170 I=NR,N2                                                     VBDIM
      TEMP=D*A(I+2,NR)                                                   VBDIM
      W(I+1)=TEMP                                                        VBDIM
  170 A(I+2,NR)=TEMP                                                     VBDIM
C     PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR                VBDIM
C     SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K)            VBDIM
      WTAW=0.0D0                                                         VBDIM
      DO 220 I=NR,N1                                                     VBDIM
      SUM=0.0D0                                                          VBDIM
      DO 180 J=NR,I                                                      VBDIM
  180 SUM=SUM+A(I+1,J+1)*W(J)                                            VBDIM
      I1=I+1                                                             VBDIM
      IF(N1-I1) 210,190,190                                              VBDIM
  190 DO 200 J=I1,N1                                                     VBDIM
  200 SUM=SUM+A(J+1,I+1)*W(J)                                            VBDIM
  210 P(I)=SUM                                                           VBDIM
  220 WTAW=WTAW+SUM*W(I)                                                 VBDIM
C     P VECTOR AND SCALAR K NOW STORED, NEXT COMPUTE Q VECTOR.           VBDIM
      DO 230 I=NR,N1                                                     VBDIM
C     NOW FORM PAP MATRIX, REQUIRED PART                                 VBDIM
  230 Q(I)=P(I)-WTAW*W(I)                                                VBDIM
      DO 240 J=NR,N1                                                     VBDIM
      QJ=Q(J)                                                            VBDIM
      WJ=W(J)                                                            VBDIM
      DO 240 I=J,N1                                                      VBDIM
  240 A(I+1,J+1)=A(I+1,J+1)-2.0D0*(W(I)*QJ+WJ*Q(I))                      VBDIM
  250 BETA(NR)=B                                                         VBDIM
      BETASQ(NR)=B*B                                                     VBDIM
  260 GAMMA(NR+1)=A(NR+1,NR+1)                                           VBDIM
  270 B=A(N,N-1)                                                         VBDIM
      BETA(N-1)=B                                                        VBDIM
      BETASQ(N-1)=B*B                                                    VBDIM
      GAMMA(N)=A(N,N)                                                    VBDIM
  280 BETASQ(N)=0.0D0                                                    VBDIM
C     ADJOIN AN IDENTITY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS        VBDIM
      DO 290 J=1,N                                                       VBDIM
      DO 290 I=1,N                                                       VBDIM
  290 VEC(I,J)=0.0D0                                                     VBDIM
      DO 300 I=1,N                                                       VBDIM
  300 VEC(I,I)=1.0D0                                                     VBDIM
      M=N                                                                VBDIM
      SUM=0.0D0                                                          VBDIM
      NPAS=1                                                             VBDIM
      GO TO 400                                                          VBDIM
  310 SUM=SUM+SHIFT                                                      VBDIM
      COSA=1.0D0                                                         VBDIM
      G=GAMMA(1)-SHIFT                                                   VBDIM
      PP=G                                                               VBDIM
      PPBS=PP*PP+BETASQ(1)                                               VBDIM
      PPBR=SQRT(PPBS)                                                    VBDIM
      DO 370 J=1,M                                                       VBDIM
      COSAP=COSA                                                         VBDIM
      IF(PPBS.NE.0.0D0) GO TO 320                                        VBDIM
  311 SINA=0.0D0                                                         VBDIM
      SINA2=0.0D0                                                        VBDIM
      COSA=1.0D0                                                         VBDIM
      GO TO 350                                                          VBDIM
  320 SINA=BETA(J)/PPBR                                                  VBDIM
      SINA2=BETASQ(J)/PPBS                                               VBDIM
      COSA=PP/PPBR                                                       VBDIM
C     POSTMULTIPLY IDAENTITY BY P-TRANSPOSE MATRIX                       VBDIM
      NT=J+NPAS                                                          VBDIM
      IF(NT.LT.N) GO TO 330                                              VBDIM
  321 NT=N                                                               VBDIM
  330 DO 340 I=1,NT                                                      VBDIM
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)                                 VBDIM
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)                          VBDIM
  340 VEC(I,J)=TEMP                                                      VBDIM
  350 DIA=GAMMA(J+1)-SHIFT                                               VBDIM
      U=SINA2*(G+DIA)                                                    VBDIM
      GAMMA(J)=G+U                                                       VBDIM
      G=DIA-U                                                            VBDIM
      PP=DIA*COSA-SINA*COSAP*BETA(J)                                     VBDIM
      IF(J.NE.M) GO TO 360                                               VBDIM
  351 BETA(J)=SINA*PP                                                    VBDIM
      BETASQ(J)=SINA2*PP*PP                                              VBDIM
      GO TO 380                                                          VBDIM
  360 PPBS=PP*PP+BETASQ(J+1)                                             VBDIM
      PPBR=SQRT(PPBS)                                                    VBDIM
      BETA(J)=SINA*PPBR                                                  VBDIM
  370 BETASQ(J)=SINA2*PPBS                                               VBDIM
  380 GAMMA(M+1)=G                                                       VBDIM
C     TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT                      VBDIM
      NPAS=NPAS+1                                                        VBDIM
      IF(BETASQ(M).GT.RHOSQ) GO TO 410                                   VBDIM
  390 EIG(M+1)=GAMMA(M+1)+SUM                                            VBDIM
  400 BETA(M)=0.0D0                                                      VBDIM
      BETASQ(M)=0.0D0                                                    VBDIM
      M=M-1                                                              VBDIM
      IF(M.EQ.0) GO TO 430                                               VBDIM
  401 IF(BETASQ(M).LE.RHOSQ) GO TO 390                                   VBDIM
C     TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN VALUE      VBDIM
C     AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT                         VBDIM
  410 A2=GAMMA(M+1)                                                      VBDIM
      R2=.5D0*A2                                                         VBDIM
      R1=.5D0*GAMMA(M)                                                   VBDIM
      R12=R1+R2                                                          VBDIM
      DIF=R1-R2                                                          VBDIM
      TEMP=SQRT(DIF*DIF+BETASQ(M))                                       VBDIM
      R1=R12+TEMP                                                        VBDIM
      R2=R12-TEMP                                                        VBDIM
      DIF=ABS(A2-R1)-ABS(A2-R2)                                          VBDIM
      IF(DIF.LT.0.0D0) GO TO 420                                         VBDIM
  411 SHIFT=R2                                                           VBDIM
      GO TO 310                                                          VBDIM
  420 SHIFT=R1                                                           VBDIM
      GO TO 310                                                          VBDIM
  430 EIG(1)=GAMMA(1)+SUM                                                VBDIM
C     INITIALIZE AUXILARY TABLES REQUIRED FOR REARRANGING THE VECTORS    VBDIM
      DO 440 J=1,N                                                       VBDIM
      IPOSV(J)=J                                                         VBDIM
      IVPOS(J)=J                                                         VBDIM
  440 IORD(J)=J                                                          VBDIM
C     USE A TRANSPOSITON SORT TO ORDER THE EIGENVALUES                   VBDIM
      M=N                                                                VBDIM
      GO TO 470                                                          VBDIM
  450 DO 460 J=1,M                                                       VBDIM
      IF(EIG(J).LE.EIG(J+1)) GO TO 460                                   VBDIM
  451 TEMP=EIG(J)                                                        VBDIM
      EIG(J)=EIG(J+1)                                                    VBDIM
      EIG(J+1)=TEMP                                                      VBDIM
      ITEMP=IORD(J)                                                      VBDIM
      IORD(J)=IORD(J+1)                                                  VBDIM
      IORD(J+1)=ITEMP                                                    VBDIM
  460 CONTINUE                                                           VBDIM
  470 M=M-1                                                              VBDIM
      IF(M.NE.0) GO TO 450                                               VBDIM
  471 IF(N1.EQ.0) GO TO 500                                              VBDIM
  472 DO 490 L=1,N1                                                      VBDIM
      NV=IORD(L)                                                         VBDIM
      NP=IPOSV(NV)                                                       VBDIM
      IF(NP.EQ.L) GO TO 490                                              VBDIM
  473 LV=IVPOS(L)                                                        VBDIM
      IVPOS(NP)=LV                                                       VBDIM
      IPOSV(LV)=NP                                                       VBDIM
      DO 480 I=1,N                                                       VBDIM
      TEMP=VEC(I,L)                                                      VBDIM
      VEC(I,L)=VEC(I,NP)                                                 VBDIM
  480 VEC(I,NP)=TEMP                                                     VBDIM
  490 CONTINUE                                                           VBDIM
C     BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX           VBDIM
  500 DO 550 NRR=1,N                                                     VBDIM
      K=N1                                                               VBDIM
  510 K=K-1                                                              VBDIM
      IF(K.LE.0) GO TO 550                                               VBDIM
  511 SUM=0.0D0                                                          VBDIM
      DO 520 I=K,N1                                                      VBDIM
  520 SUM=SUM+VEC(I+1,NRR)*A(I+1,K)                                      VBDIM
      SUM=SUM+SUM                                                        VBDIM
      DO 530 I=K,N1                                                      VBDIM
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(I+1,K)                             VBDIM
      GO TO 510                                                          VBDIM
  550 CONTINUE                                                           VBDIM
  560 RETURN                                                             VBDIM
C  END OF EIGN
      END                                                                VBDIM
      BLOCK DATA PTPARM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT4CM/IPTPRT, IPTRD1
      COMMON /PT2CM/ NDER, NDUM, IDEBUG, IDUM(19)
      COMMON /POTCM2/ C(3)
      COMMON/POTCM3/ZERO,ONE,TWO,TRE,HALF,QUAR,RT3,
     $DPUX,SSGX,DSGH,SSGH
C   Initialize the control flags for the potenetial
      DATA IPTPRT, IPTRD1 /6, 4/
      DATA NDUM, NDER /0, 1/
      DATA IDEBUG /0/
      DATA IDUM /19*0/ 
C   Initialize the parameters for the potential
      DATA ZERO,ONE,TWO,TRE/0.D0,1.D0,2.D0,3.D0/
      DATA HALF,QUAR/0.5D0,0.25D0/
      DATA RT3/1.7320 50807 56888D0/
      DATA DPUX,SSGX,DSGH,SSGH/0.D0,-0.132774D0,0.D0,0.5D0/
      DATA C/0.169552D0,0.1744734D0,0.169552D0/
      END

      SUBROUTINE POTINFO

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
      PARAMETER (PI = 3.141592653589793D0)
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
      write(NFLAG(18),96)
 96   format(/)
      do i =1,5
         write(NFLAG(18),97) REF(i)
      end do
 97   format(2x,a75)
      WRITE(NFLAG(18),96)
      KMAX = 0
      DO I = 1,ISURF+1
         DO J = 1,ISURF+1
            IF(NASURF(I,J).NE.0.AND.KMAX.LT.MAX(I,J)) KMAX = MAX(I,J)
         ENDDO
      ENDDO
      WRITE(NFLAG(18),101) MSURF,KMAX-1
101   FORMAT(2x,' MAX. AND ACTUAL NO. OF EXCITED SURFACES: ',I3,5x,I3)
      IF(KMAX-1.GT.MSURF) THEN
         WRITE(6,*) ' WRONG INPUT ON NUMBER OF EXCITED SURFACES'
         STOP
      ENDIF
      KSDIAG = 0
      KEDIAG = 0
      DO I = 2,ISURF+1
         IF(NASURF(I,I).NE.0) THEN
            KEDIAG = I-1
            IF(KSDIAG.EQ.0) KSDIAG = I-1
         ENDIF
      ENDDO
      KSOFFD = 0
      KEOFFD = 0
      K = 0
      DO I = 1,ISURF
         DO J = I+1,ISURF+1
            K = K+1
            IF(NASURF(I,J)+NASURF(J,I).NE.0) THEN
               KEOFFD = K
               IF(KSOFFD.EQ.0) KSOFFD = K
            ENDIF
         ENDDO
      ENDDO
      WRITE(NFLAG(18),103) MDER,NDER
103   FORMAT(2x,' MAX. AND ACTUAL ORDER OF DERIVATIVES:    ',I3,5x,I3)
      IF(NDER.GT.MDER) THEN
         WRITE(6,*) ' WRONG INPUT ON ORDER OF DERIVATIVES'
         STOP
      ENDIF
      IF(NFLAG(19).EQ.1) THEN
         write(NFLAG(18),100)
 100     format(/)
         write(NFLAG(18),120)
 120     format(2x,'Cartesian coordinates are supplied by',/,
     +          2x,'the user in the array CART.',//)
         write(NFLAG(18),125)
 125     format(2x,'Provide cartesian coordinates in the',/,
     +          2x,'following order using the array CART',//,
     +          2x,' CART(1,1)...CART(1,3)   => ATOM 1',/,
     +          2x,' CART(2,1)...CART(2,3)   => ATOM 2',/,
     +          2x,' CART(3,1)...CART(3,3)   => ATOM 3',/,
     +          2x,' CART(N,1)...CART(N,3)   => ATOM N',/,
     +          2x,'CART(25,1)...CART(25,3)  => ATOM 25',/)
         write(NFLAG(18),130)
 130     format(2x,'If the user wishes to relabel the atoms,',/,
     +          2x,'set the variable IREORDER equal to 1',/,
     +          2x,'in the PARAMETER statement.  The user',/,
     +          2x,'must also specify the new labeling',/,
     +          2x,'scheme.  This is done using the array',/,
     +          2x,'NULBL in the following manner:',//,
     +          2x,'NULBL(i) = j',/,
     +          2x,'where:  i => old index',/,
     +          2x,'        j => new index',//)
         write(NFLAG(18),150)
 150     format(2x,'Cartesian coordinates can be provided to',/,
     +          2x,'the potential routine in a variety of units.',/,
     +          2x,'The input units will be converted to Bohr',/,
     +          2x,'based on the following values of the NFLAG',/,
     +          2x,'variable:',//,
     +          2x,'NFLAG(1)  =  1  =>  CARTESIANS IN BOHR (no',/,
     +          2x,'                    conversion required)',/,
     +          2x,'NFLAG(1)  =  2  =>  CARTESIANS IN ANGSTROMS',//)
         write(NFLAG(18),160)
 160     format(2x,'The value of the energy and derivatives',/,
     +          2x,'(if computed) can be reported in a variety',/,
     +          2x,'units.  A units conversion will take place',/,
     +          2x,'as specified by the following values of the',/,
     +          2x,'NFLAG variable:',//,
     +          2x,'NFLAG(2) = 1 =>  ENERGIES REPORTED IN HARTEEE',/,
     +          2x,'NFLAG(2) = 2 =>  ENERGIES REPORTED IN mHARTREE',/,
     +          2x,'NFLAG(2) = 3 =>  ENERGIES REPORTED IN eV',/,
     +          2x,'NFLAG(2) = 4 =>  ENERGIES REPORTED IN kcal/mol',/,
     +          2x,'NFLAG(2) = 5 =>  ENERGIES REPORTED IN cm**-1',//)
         write(NFLAG(18),165)
 165     format(2x,'A units conversion will take place',/,
     +       2x,'as specified by the following values of the',/,
     +       2x,'NFLAG variable:',//,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=1 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           HARTEEE/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=2 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           mHARTREE/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=3 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           eV/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=4 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           kcal/mol/BOHR',/,
     +       2x,'NFLAG(1)=1 & NFLAG(2)=5 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           cm**-1/BOHR',//)
         write(NFLAG(18),170)
 170     format(2x,'A units conversion will take place',/,
     +       2x,'as specified by the following values of the',/,
     +       2x,'NFLAG variable:',//,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=1 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           HARTEEE/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=2 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           mHARTREE/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=3 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           eV/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=4 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           kcal/mol/ANGSTROM',/,
     +       2x,'NFLAG(1)=2 & NFLAG(2)=5 => DERIVATIVES REPORTED IN',/,
     +       2x,'                           cm**-1/ANGSTROM',//)
      ENDIF
      RETURN
      END


      SUBROUTINE ANCVRT
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*75 REF(5)
      CHARACTER*3 PERIODIC_1(7,32)
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
      CHARACTER*2 NAME1(NATOM)
      CHARACTER*2 NAME2(NATOM)
      CHARACTER*1 IBLANK
      CHARACTER*20 DISTANCE
      CHARACTER*20 UNITS
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      DIMENSION IANUM(7,32)
      DIMENSION ISAVE(NATOM),JSAVE(NATOM)
      PARAMETER(        PI = 3.141592653589793D0)
      PARAMETER(    CLIGHT = 2.99792458D08)
      PARAMETER(     CMU_0 = 4.0D0*PI*1.0D-07)
      PARAMETER(CEPSILON_0 = 1.0D0/(CMU_0*CLIGHT**2))
      PARAMETER(        CE = 1.602176462D-19)
      PARAMETER(   CPLANCK = 6.62606876D-34)
      PARAMETER(      CM_E = 9.10938188D-31)
      PARAMETER(      CANG = 1.0D-10)
      PARAMETER( CAVOGADRO = 6.02214199D23)
      PARAMETER(     CKCAL = 4.184D10)
      PARAMETER(  HTOMILLH = 1000.D0)
      PARAMETER(     HTOEV = 27.2113834D0)
      PARAMETER(   HTOKCAL = 627.509470D0)
      PARAMETER(   HTOWAVE = 219474.631D0)
      PARAMETER(     HTOKJ = 2625.49962D0)
      PARAMETER(    BOHR_A = .5291772083D0)
      DO I=1,7
         DO J=1,32
            IANUM(I,J)=0
            PERIODIC_1(I,J)=' '
         END DO
      END DO
      DISTANCE = 'BOHR                '
      UNITS    = 'HARTREE             '
      IANUM(1,1)  =  1
      IANUM(1,32) =  2
      IANUM(2,1)  =  3
      IANUM(2,2)  =  4
      IANUM(2,27) =  5
      IANUM(2,28) =  6
      IANUM(2,29) =  7
      IANUM(2,30) =  8
      IANUM(2,31) =  9
      IANUM(2,32) = 10
      IANUM(3,1)  = 11
      IANUM(3,2)  = 12
      IANUM(3,27) = 13
      IANUM(3,28) = 14
      IANUM(3,29) = 15
      IANUM(3,30) = 16
      IANUM(3,31) = 17
      IANUM(3,32) = 18
      IANUM(4,1)  = 19
      IANUM(4,2)  = 20
      IANUM(4,17) = 21
      IANUM(4,18) = 22
      IANUM(4,19) = 23
      IANUM(4,20) = 24
      IANUM(4,21) = 25
      IANUM(4,22) = 26
      IANUM(4,23) = 27
      IANUM(4,24) = 28
      IANUM(4,25) = 29
      IANUM(4,26) = 30
      IANUM(4,27) = 31
      IANUM(4,28) = 32
      IANUM(4,29) = 33
      IANUM(4,30) = 34
      IANUM(4,31) = 35
      IANUM(4,32) = 36
      IANUM(5,1)  = 37
      IANUM(5,2)  = 38
      IANUM(5,17) = 39
      IANUM(5,18) = 40
      IANUM(5,19) = 41
      IANUM(5,20) = 42
      IANUM(5,21) = 43
      IANUM(5,22) = 44
      IANUM(5,23) = 45
      IANUM(5,24) = 46
      IANUM(5,25) = 47
      IANUM(5,26) = 48
      IANUM(5,27) = 49
      IANUM(5,28) = 50
      IANUM(5,29) = 51
      IANUM(5,30) = 52
      IANUM(5,31) = 53
      IANUM(5,32) = 54
      IANUM(6,1)  = 55
      IANUM(6,2)  = 56
      IANUM(6,3)  = 57
      IANUM(6,4)  = 58
      IANUM(6,5)  = 59
      IANUM(6,6)  = 60
      IANUM(6,7)  = 61
      IANUM(6,8)  = 62
      IANUM(6,9)  = 63
      IANUM(6,10) = 64
      IANUM(6,11) = 65
      IANUM(6,12) = 66
      IANUM(6,13) = 67
      IANUM(6,14) = 68
      IANUM(6,15) = 69
      IANUM(6,16) = 70
      IANUM(6,17) = 71
      IANUM(6,18) = 72
      IANUM(6,19) = 73
      IANUM(6,20) = 74
      IANUM(6,21) = 75
      IANUM(6,22) = 76
      IANUM(6,23) = 77
      IANUM(6,24) = 78
      IANUM(6,25) = 79
      IANUM(6,26) = 80
      IANUM(6,27) = 81
      IANUM(6,28) = 82
      IANUM(6,29) = 83
      IANUM(6,30) = 84
      IANUM(6,31) = 85
      IANUM(6,32) = 86
      IANUM(7,1)  = 87
      IANUM(7,2)  = 88
      IANUM(7,3)  = 89
      IANUM(7,4)  = 90
      IANUM(7,5)  = 91
      IANUM(7,6)  = 92
      IANUM(7,7)  = 93
      IANUM(7,8)  = 94
      IANUM(7,9)  = 95
      IANUM(7,10) = 96
      IANUM(7,11) = 97
      IANUM(7,12) = 98
      IANUM(7,13) = 99
      IANUM(7,14) = 100
      IANUM(7,15) = 101
      IANUM(7,16) = 102
      IANUM(7,17) = 103
      IANUM(7,18) = 104
      IANUM(7,19) = 105
      IANUM(7,20) = 106
      IANUM(7,21) = 107
      IANUM(7,22) = 108
      IANUM(7,23) = 109
      IANUM(7,24) = 110
      IANUM(7,25) = 111
      IANUM(7,26) = 112
      IANUM(7,27) = 113
      IANUM(7,28) = 114
      IANUM(7,29) = 115
      IANUM(7,30) = 116
      IANUM(7,31) = 117
      IANUM(7,32) = 120
      PERIODIC_1(1,1)   = 'H  '
      PERIODIC_1(1,32)  = 'He '
      PERIODIC_1(2,1)   = 'Li '
      PERIODIC_1(2,2)   = 'Be '
      PERIODIC_1(2,27)  = 'B  '
      PERIODIC_1(2,28)  = 'C  '
      PERIODIC_1(2,29)  = 'N  '
      PERIODIC_1(2,30)  = 'O  '
      PERIODIC_1(2,31)  = 'F  '
      PERIODIC_1(2,32)  = 'Ne '
      PERIODIC_1(3,1)   = 'Na '
      PERIODIC_1(3,2)   = 'Mg '
      PERIODIC_1(3,27)  = 'Al '
      PERIODIC_1(3,28)  = 'Si '
      PERIODIC_1(3,29)  = 'P  '
      PERIODIC_1(3,30)  = 'S  '
      PERIODIC_1(3,31)  = 'Cl '
      PERIODIC_1(3,32)  = 'Ar '
      PERIODIC_1(4,1)   = 'K  '
      PERIODIC_1(4,2)   = 'Ca '
      PERIODIC_1(4,17)  = 'Sc '
      PERIODIC_1(4,18)  = 'Ti '
      PERIODIC_1(4,19)  = 'V  '
      PERIODIC_1(4,20)  = 'Cr '
      PERIODIC_1(4,21)  = 'Mn '
      PERIODIC_1(4,22)  = 'Fe '
      PERIODIC_1(4,23)  = 'Co '
      PERIODIC_1(4,24)  = 'Ni '
      PERIODIC_1(4,25)  = 'Cu '
      PERIODIC_1(4,26)  = 'Zn '
      PERIODIC_1(4,27)  = 'Ga '
      PERIODIC_1(4,28)  = 'Ge '
      PERIODIC_1(4,29)  = 'As '
      PERIODIC_1(4,30)  = 'Se '
      PERIODIC_1(4,31)  = 'Br '
      PERIODIC_1(4,32)  = 'Kr '
      PERIODIC_1(5,1)   = 'Rb '
      PERIODIC_1(5,2)   = 'Sr '
      PERIODIC_1(5,17)  = 'Y  '
      PERIODIC_1(5,18)  = 'Zr '
      PERIODIC_1(5,19)  = 'Nb '
      PERIODIC_1(5,20)  = 'Mo '
      PERIODIC_1(5,21)  = 'Tc '
      PERIODIC_1(5,22)  = 'Ru '
      PERIODIC_1(5,23)  = 'Rh '
      PERIODIC_1(5,24)  = 'Pd '
      PERIODIC_1(5,25)  = 'Ag '
      PERIODIC_1(5,26)  = 'Cd '
      PERIODIC_1(5,27)  = 'In '
      PERIODIC_1(5,28)  = 'Sn '
      PERIODIC_1(5,29)  = 'Sb '
      PERIODIC_1(5,30)  = 'Te '
      PERIODIC_1(5,31)  = 'I  '
      PERIODIC_1(5,32)  = 'Xe '
      PERIODIC_1(5,32)  = 'Xe '
      DO I=1,NATOMS
         ISAVE(I)=0
         JSAVE(I)=0
         NAME1(I)='  '
         NAME2(I)='  '
      END DO
      IBLANK=' '
      DO IND=1,NATOMS
         DO I=1,7
            DO J=1,32
               IF(INDEXES(IND).EQ.IANUM(I,J)) THEN
                  ISAVE(IND)=I
                  JSAVE(IND)=J
               END IF
            END DO
         END DO
      END DO
 
      DO IND=1,NATOMS
         IND2=NULBL(IND)
         IF(IND2.EQ.0) IND2=IND
      END DO
      INC1=0
      DO IND=1,IRCTNT-1
         INC1=INC1+1
         NAME1(INC1)=PERIODIC_1(ISAVE(IND),JSAVE(IND))(:2)
      END DO
      INC2=0
      DO IND=IRCTNT,NATOMS
         INC2=INC2+1
         NAME2(INC2)=PERIODIC_1(ISAVE(IND),JSAVE(IND))(:2)
      END DO
      IF(NFLAG(1).EQ.2) DISTANCE = 'ANGSTROMS           '
      IF(NFLAG(2).EQ.2) THEN
         UNITS = 'MILLIHARTREE        '
      ELSEIF(NFLAG(2).EQ.3) THEN
         UNITS = 'EV                  '
      ELSEIF(NFLAG(2).EQ.4) THEN
         UNITS = 'KCAL PER MOLE       '
      ELSEIF(NFLAG(2).EQ.5) THEN
         UNITS = 'WAVENUMBERS         '
      ELSEIF(NFLAG(2).EQ.6) THEN
         UNITS = 'KILOJOULES PER MOLE '
      ENDIF
      CNVRTD = 1.D0
      CNVRTE = 1.D0
      CNVRTDE = 1.D0
      IF(NFLAG(1).EQ.2) CNVRTD = BOHR_A
      IF(NFLAG(2).EQ.2) THEN
         CNVRTE = CNVRTE*HTOMILLH
      ELSEIF(NFLAG(2).EQ.3) THEN
         CNVRTE = CNVRTE*HTOEV
      ELSEIF(NFLAG(2).EQ.4) THEN
         CNVRTE = CNVRTE*HTOKCAL
      ELSEIF(NFLAG(2).EQ.5) THEN
         CNVRTE = CNVRTE*HTOWAVE
      ELSEIF(NFLAG(2).EQ.6) THEN
         CNVRTE = CNVRTE*HTOKJ
      ENDIF
      CNVRTDE = CNVRTE/CNVRTD
      ISUM = 0
      DO INU=1,25
         ISUM=ISUM + NULBL(INU)
      END DO
      IREORDER = 0
      IF(ISUM.NE.0) IREORDER = 1
      RETURN
      END
      SUBROUTINE CARTOU
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      IF (IREORDER.EQ.1) THEN
          DO I=1,NATOMS
             DO J=1,3
                CARTNU(NULBL(I),J)=CART(I,J)/CNVRTD
             END DO
          END DO
      ELSE
          DO I=1,NATOMS
             DO J=1,3
                CARTNU(I,J)=CART(I,J)/CNVRTD
             END DO
          END DO
      END IF
      RETURN
      END
 
      SUBROUTINE CARTTOR
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF=5)
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      IF(ICARTR.EQ.1) THEN
         DO I=1,NATOMS
            IND=3*I-2
            R(IND)   = CARTNU(I,1)
            R(IND+1) = CARTNU(I,2)
            R(IND+2) = CARTNU(I,3)
         END DO
      ELSEIF(ICARTR.EQ.2) THEN
         I = 1                                                       
         DO K=1,NATOMS-1
            DO L = K+1,NATOMS                                  
               R(I) = SQRT( (CARTNU(K,1)-CARTNU(L,1))**2 +
     +                      (CARTNU(K,2)-CARTNU(L,2))**2 +
     +                      (CARTNU(K,3)-CARTNU(L,3))**2 )
               I = I + 1                  
            END DO
         ENDDO
      ELSEIF(ICARTR.EQ.3) THEN
         R(1) = SQRT( (CARTNU(1,1)-CARTNU(2,1))**2 +
     +                (CARTNU(1,2)-CARTNU(2,2))**2 +
     +                (CARTNU(1,3)-CARTNU(2,3))**2 )
         R(2) = SQRT( (CARTNU(2,1)-CARTNU(3,1))**2 +
     +                (CARTNU(2,2)-CARTNU(3,2))**2 +
     +                (CARTNU(2,3)-CARTNU(3,3))**2 )
         R(3) = SQRT( (CARTNU(1,1)-CARTNU(3,1))**2 +
     +                (CARTNU(1,2)-CARTNU(3,2))**2 +
     +                (CARTNU(1,3)-CARTNU(3,3))**2 )
      ELSEIF(ICARTR.EQ.4) THEN
      FLM=18.99840D0
      HYM=1.007825D0
      XCM1=(HYM*CARTNU(1,1)+FLM*CARTNU(2,1))/(FLM+HYM)
      YCM1=(HYM*CARTNU(1,2)+FLM*CARTNU(2,2))/(FLM+HYM)
      ZCM1=(HYM*CARTNU(1,3)+FLM*CARTNU(2,3))/(FLM+HYM)
      XCM2=(HYM*CARTNU(3,1)+FLM*CARTNU(4,1))/(FLM+HYM)
      YCM2=(HYM*CARTNU(3,2)+FLM*CARTNU(4,2))/(FLM+HYM)
      ZCM2=(HYM*CARTNU(3,3)+FLM*CARTNU(4,3))/(FLM+HYM)
      XCM3=XCM2-XCM1
      YCM3=YCM2-YCM1
      ZCM3=ZCM2-ZCM1
      XRM1=CARTNU(1,1)-XCM1
      YRM1=CARTNU(1,2)-YCM1
      ZRM1=CARTNU(1,3)-ZCM1
      THETA1=(XRM1*XCM3+YRM1*YCM3+ZRM1*ZCM3)
      THETA1=THETA1/(SQRT(XRM1**2+YRM1**2+ZRM1**2))
      THETA1=THETA1/(SQRT(XCM3**2+YCM3**2+ZCM3**2))
      IF(THETA1.GT.1.0D0)THETA1=1.0D0
      IF(THETA1.LT.-1.0D0)THETA1=-1.0D0
      THETA1=ACOS(THETA1)
      XRM2=CARTNU(3,1)-XCM2
      YRM2=CARTNU(3,2)-YCM2
      ZRM2=CARTNU(3,3)-ZCM2
      THETA2=(XRM2*(-XCM3)+YRM2*(-YCM3)+ZRM2*(-ZCM3))
      THETA2=THETA2/(SQRT(XRM2**2+YRM2**2+ZRM2**2))
      THETA2=THETA2/(SQRT(XCM3**2+YCM3**2+ZCM3**2))
      IF(THETA2.GT.1.0D0)THETA2=1.0D0
      IF(THETA2.LT.-1.0D0)THETA2=-1.0D0
      THETA2=ACOS(THETA2)
      PI=ACOS(-1.0D0)
      THETA2=PI-THETA2
      Q1=SQRT(XRM1**2+YRM1**2+ZRM1**2)
      Q2=SQRT(XRM2**2+YRM2**2+ZRM2**2)
      CMM=(XCM3**2+YCM3**2+ZCM3**2)
      CMM=SQRT(CMM)
      HHD=(CARTNU(1,1)-CARTNU(3,1))**2 +
     +    (CARTNU(1,2)-CARTNU(3,2))**2 +
     +    (CARTNU(1,3)-CARTNU(3,3))**2
      HHD=SQRT(HHD)
      Q=CMM-Q1*COS(THETA1)+Q2*COS(THETA2)
      Q3=SQRT(ABS(HHD**2-Q**2))
      Q1=Q1*SIN(THETA1)
      Q2=Q2*SIN(THETA2)
      CPHI=(Q1**2+Q2**2-Q3**2)/(2.*Q1*Q2)
      IF(CPHI.LT.-1.0D0)CPHI=-1.0D0
      IF(CPHI.GT.1.0D0)CPHI=1.0D0
      PHI=ACOS(CPHI)
 2001 FORMAT(6F12.8)
      R(1)=SQRT(XCM3**2+YCM3**2+ZCM3**2)
      R(2)=(SQRT(XRM1**2+YRM1**2+ZRM1**2))*(FLM+HYM)/FLM
      R(3)=(SQRT(XRM2**2+YRM2**2+ZRM2**2))*(FLM+HYM)/FLM
      R(4)=THETA1
      R(5)=THETA2
      R(6)=PHI
      ELSEIF(ICARTR.NE.0) THEN
         WRITE(NFLAG(18),1000) ICARTR
 1000    FORMAT(2X,'WRONG ICARTR FOR CARTNU; ICARTR =',I5//)
         STOP
      ENDIF
      RETURN
      END
 
 
      SUBROUTINE EUNITZERO
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
 
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +               DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,
     +               CNVRTDE,IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      PENGYGS = ENGYGS * CNVRTE - ANUZERO
      IF(KSDIAG.NE.0) THEN
         DO I=KSDIAG,KEDIAG
            PENGYES(I) = ENGYES(I) * CNVRTE - ANUZERO
         END DO
      ENDIF
      IF(KSOFFD.NE.0) THEN
         DO J=KSOFFD,KEOFFD
            PENGYIJ(J) = ENGYIJ(J) * CNVRTE
         END DO
      ENDIF
      RETURN
      END
 
      SUBROUTINE RTOCART
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT3CM/  EZERO(ISURF+1)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
      COMMON /UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +                DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
     +                IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      DIMENSION YGS(N3ATOM),YES(N3ATOM,ISURF),YIJ(N3ATOM,JSURF)
      IF(ICARTR.EQ.1) THEN
         DO I = 1, NATOMS
            IND=3*I-2
            DGSCARTNU(I,1) = DEGSDR(IND)
            DGSCARTNU(I,2) = DEGSDR(IND+1)
            DGSCARTNU(I,3) = DEGSDR(IND+2)
            IF(KSDIAG.NE.0) THEN
               DO J = KSDIAG,KEDIAG
                  DESCARTNU(I,1,J) = DEESDR(IND,J)
                  DESCARTNU(I,2,J) = DEESDR(IND+1,J)
                  DESCARTNU(I,3,J) = DEESDR(IND+2,J)
               END DO
            ENDIF
            IF(KEOFFD.NE.0) THEN
               DO K = KSOFFD,KEOFFD
                  DIJCARTNU(I,1,K) = DEIJDR(IND,K)
                  DIJCARTNU(I,2,K) = DEIJDR(IND+1,K)
                  DIJCARTNU(I,3,K) = DEIJDR(IND+2,K)
               END DO
            ENDIF
         END DO
      ELSEIF(ICARTR.EQ.2) THEN
         DO I = 1, NATOMS         
            DGSCARTNU(I,1) = 0.D0
            DGSCARTNU(I,2) = 0.D0
            DGSCARTNU(I,3) = 0.D0
            IF(KSDIAG.NE.0) THEN
               DO J1=KSDIAG,KEDIAG
                  DESCARTNU(I,1,J1) = 0.D0
                  DESCARTNU(I,2,J1) = 0.D0
                  DESCARTNU(I,3,J1) = 0.D0
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO J2=KSOFFD,KEOFFD
                  DIJCARTNU(I,1,J2) = 0.D0
                  DIJCARTNU(I,2,J2) = 0.D0
                  DIJCARTNU(I,3,J2) = 0.D0
               ENDDO
            ENDIF
            DO J = 1,NATOMS
               IF(J.LT.I) THEN
                  M1 = NATOMS*(J-1) - (J*(J-1))/2 + I-J
               ELSEIF(J.GT.I) THEN
                  M1 = NATOMS*(I-1) - (I*(I-1))/2 + J-I
               ELSE
                  GO TO 20
               ENDIF
               Y = DEGSDR(M1)
               TERMX = (CARTNU(I,1)-CARTNU(J,1))/R(M1)
               TERMY = (CARTNU(I,2)-CARTNU(J,2))/R(M1)
               TERMZ = (CARTNU(I,3)-CARTNU(J,3))/R(M1)
               DGSCARTNU(I,1) = DGSCARTNU(I,1) + TERMX*Y
               DGSCARTNU(I,2) = DGSCARTNU(I,2) + TERMY*Y
               DGSCARTNU(I,3) = DGSCARTNU(I,3) + TERMZ*Y
               IF(KSDIAG.GT.0) THEN
                  Y = DEESDR(M1,J1)
                  DO J1=KSDIAG,KEDIAG
                     DESCARTNU(I,1,J1)=DESCARTNU(I,1,J1) + TERMX*Y
                     DESCARTNU(I,2,J1)=DESCARTNU(I,2,J1) + TERMY*Y
                     DESCARTNU(I,3,J1)=DESCARTNU(I,3,J1) + TERMZ*Y
                  ENDDO
               ELSEIF(KSOFFD.GT.0) THEN
                  DO J2=KSOFFD,KEOFFD
                     Y = DEIJDR(M1,J2)
                     DIJCARTNU(I,1,J2)=DIJCARTNU(I,1,J2) + TERMX*Y
                     DIJCARTNU(I,2,J2)=DIJCARTNU(I,2,J2) + TERMY*Y
                     DIJCARTNU(I,3,J2)=DIJCARTNU(I,3,J2) + TERMZ*Y
                  ENDDO
               ENDIF
20             CONTINUE
            ENDDO
         ENDDO
      ELSEIF(ICARTR.EQ.3) THEN
         DO I = 1, NATOMS
            YGS(I) = DEGSDR(I)/R(I)
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  YES(I,J) = DEESDR(I,J)/R(I)
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  YIJ(I,K) = DEIJDR(I,K)/R(I)
               ENDDO
            ENDIF
         ENDDO
         DO K = 1,3
            TERM12 = CARTNU(1,K)-CARTNU(2,K)
            TERM23 = CARTNU(2,K)-CARTNU(3,K)
            TERM13 = CARTNU(1,K)-CARTNU(3,K)
            DGSCARTNU(1,K) = TERM12*YGS(1) + TERM13*YGS(3)
            DGSCARTNU(2,K) =-TERM12*YGS(1) + TERM23*YGS(2)
            DGSCARTNU(3,K) =-TERM13*YGS(3) - TERM23*YGS(2)
            IF(KSDIAG.NE.0) THEN
               DO J1=KSDIAG,KEDIAG
                 DESCARTNU(1,K,J1) = TERM12*YES(1,J1) + TERM13*YES(3,J1)
                 DESCARTNU(2,K,J1) =-TERM12*YES(1,J1) + TERM23*YES(2,J1)
                 DESCARTNU(3,K,J1) =-TERM13*YES(3,J1) - TERM23*YES(2,J1)
               ENDDO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO J2=KSOFFD,KEOFFD
                 DIJCARTNU(1,K,J2) = TERM12*YIJ(1,J2) + TERM13*YIJ(3,J2)
                 DIJCARTNU(2,K,J2) =-TERM12*YIJ(1,J2) + TERM23*YIJ(2,J2)
                 DIJCARTNU(3,K,J2) =-TERM13*YIJ(3,J2) - TERM23*YIJ(2,J2)
               ENDDO
            ENDIF
         ENDDO
      ELSEIF(ICARTR.EQ.4) THEN
      WH=1.007825D0
      WF=18.99840D0
      SUM=WH+WF
      EPS=WF/SUM
      EPSP=WH/SUM
      U1=COS(THETA1)
      U2=COS(THETA2)
      U3=COS(PHI)
      SS1=SIN(THETA1)
      SS2=SIN(THETA2)
      SS3=SIN(PHI)
      YA=0.0D0
      YB=0.0D0
      T0=R1*U1
      ZA=-EPSP*T0
      ZB=EPS*T0
      T0=R1*SS1
      XA=-EPSP*T0
      XBB=EPS*T0
      T0=R2*SS2
      T1=T0*U3
      XC=-EPSP*T1
      XD=EPS*T1
      T1=T0*SS3
      YC=-EPSP*T1
      YD=EPS*T1
      T0=R2*U2
      ZC=-EPSP*T0+RCM
      ZD=EPS*T0+RCM
      RFF=SQRT((XA-XC)**2+YC**2+(ZA-ZC)**2)
      ELSE
          WRITE(NFLAG(18),1000) ICARTR
1000      FORMAT(2X,' WRONG ICARTR FOR DERIVATIVE; ICARTR =',I5//)
          STOP
      ENDIF
      RETURN
      END
 
      SUBROUTINE DEDCOU
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*75 REF(5)
      PARAMETER (N3ATOM=75)
      PARAMETER (NATOM=25)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
      COMMON /UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     +                DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
     +                IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      IF (IREORDER.EQ.1) THEN
         DO I = 1, NATOMS
            DGSCART(I,1) = DGSCARTNU(NULBL(I),1) * CNVRTDE
            DGSCART(I,2) = DGSCARTNU(NULBL(I),2) * CNVRTDE
            DGSCART(I,3) = DGSCARTNU(NULBL(I),3) * CNVRTDE
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  DESCART(I,1,J) = DESCARTNU(NULBL(I),1,J) * CNVRTDE
                  DESCART(I,2,J) = DESCARTNU(NULBL(I),2,J) * CNVRTDE
                  DESCART(I,3,J) = DESCARTNU(NULBL(I),3,J) * CNVRTDE
               END DO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  DIJCART(I,1,K) = DIJCARTNU(NULBL(I),1,K) * CNVRTDE
                  DIJCART(I,2,K) = DIJCARTNU(NULBL(I),2,K) * CNVRTDE
                  DIJCART(I,3,K) = DIJCARTNU(NULBL(I),3,K) * CNVRTDE
               END DO
            ENDIF
         END DO
      ELSE
         DO I = 1, NATOMS
            DGSCART(I,1) = DGSCARTNU(I,1) * CNVRTDE
            DGSCART(I,2) = DGSCARTNU(I,2) * CNVRTDE
            DGSCART(I,3) = DGSCARTNU(I,3) * CNVRTDE
            IF(KSDIAG.NE.0) THEN
               DO J=KSDIAG,KEDIAG
                  DESCART(I,1,J) = DESCARTNU(I,1,J) * CNVRTDE
                  DESCART(I,2,J) = DESCARTNU(I,2,J) * CNVRTDE
                  DESCART(I,3,J) = DESCARTNU(I,3,J) * CNVRTDE
               END DO
            ENDIF
            IF(KSOFFD.NE.0) THEN
               DO K=KSOFFD,KEOFFD
                  DIJCART(I,1,K) = DIJCARTNU(I,1,K) * CNVRTDE
                  DIJCART(I,2,K) = DIJCARTNU(I,2,K) * CNVRTDE
                  DIJCART(I,3,K) = DIJCARTNU(I,3,K) * CNVRTDE
               END DO
            ENDIF
         END DO
      ENDIF
      RETURN
      END

