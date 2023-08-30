      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=6
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      PARAMETER (NATOM=25)
      PARAMETER (ISURF=5)
      PARAMETER (JSURF=INT(ISURF*(ISURF+1)/2))

      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      logical, save :: first_time_data=.true.

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      CART=0.d0
      do iatom=1,natoms
      do idir=1,3
        CART(iatom,idir)=x(iatom,idir)/0.529177211
      enddo
      enddo

      if(first_time_data) then
      call prepot
      first_time_data=.false.
      endif

      call pot

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=PENGYGS*27.211386
        enddo
      else if (igrad==1) then
        do istate=1,nstates
          p(istate)=PENGYGS*27.211386
        enddo
        do iatom=1,natoms
        do idir=1,3
          g(1,iatom,idir)=DGSCART(iatom,idir)*51.422067
        enddo
        enddo
      else if (igrad==2) then
        write (*,*) 'Only energy and gradient are available'
      endif

      endsubroutine

C
C System:          CH5
C Functional form:
C Common name:     J1
C Reference  :     T. Joseph, R. Steckler, and D. G. Truhlar
C                  J. Chem. Phys. Vol. 87, p. 7036, 1987
C 
C Number of bodies: 6
C Number of electronic states: 1
C Number of derivatives: 0
C Interface: 6-1S
C Notes: This routine requires the utility 'utility.f' located in the 
C        'utilities' section of POTLIB-online.
C
C***********************************************************************
      SUBROUTINE POT 
C
C***********************************************************************
C     SUBROUTINE TO CALCULATE THE CH4-H SURFACE
C     ALSO COMPUTES DE/DR FOR ALL R AND STORES RESULTS IN ARRAY DER
C     DEFINITION OF CARTESIAN COORDINATES X(I), Y(I), AND Z(I)
C      I = 1, 2, 3, 4   REFERS TO H1, H2, H3, H4
C      I = 5            REFERS TO C
C      I = 6            REFERS TO H'
C
C     DEFINITION OF THE R-ARRAY  AR(I)
C     H1 IS THE ATOM THAT WILL FORM H-H WITH H5 
C      I = 1, 2, 3, 4   REFERS TO C - H(I)
C      I = 5            REFERS TO C - H'
C      I = 6, 7, 8, 9   REFERS TO H' - H(I)
C      I = 10, 11, 12   REFERS TO H1 - H(2,3,4)
C      I = 13, 14       REFERS TO H2 - H(3,4)
C      I = 15           REFERS TO H3 - H4
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      DIMENSION DFR(4),DVOOP(18),Y1(6),Y2(6),F(4),DF(4),DK(6)
C
      DIMENSION COORD(18),DXX(18)
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /POTXCM/ X(6),Y(6),Z(6),ENERGY,DRDX(15,18),DEDX(18)
      COMMON /POTRCM/ AG(6),R2(15),ACS(6),ASS(6)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
      COMMON /DERVCP/ DECCR1,DECCR5,DECCR6
      COMMON /TEMCOM/ RTEM(9),DTEMP(9),ETEM
      COMMON /ALPDCM/ DALP
      DATA EVCON / 2.721161D+1 /
C
      CALL CARTOU
      CALL CARTTOR
C
C     PUT COORDINATES IN PROPER ARRAYS
C
      DO 1 I=1,18
         COORD(I)=R(I)
   1  CONTINUE
C
      IPT = 1
      DO 5 I = 1, 6
         X(I) = COORD(IPT)
         Y(I) = COORD(IPT+1)
         Z(I) = COORD(IPT+2)
         IPT = IPT+3
   5  CONTINUE
C
C COMPUTE THE INTERPARTICLE DISTANCES R
C
      CALL DISTA
C
C***********************************************************************
C CHECK FOR THE FIFTH ATOM--I.E.--CH4---(H)
C BY DEFINITION THE FIFTH H ATOM IS THE ONE WHOSE C-H DISTANCE IS MAXIM
C***********************************************************************
C
      RHMAX = AR(1)
      KMAX = 1
      IF (AR(2).LT.RHMAX) GO TO 10
      RHMAX = AR(2)
      KMAX = 2
   10 IF (AR(3).LT.RHMAX) GO TO 20
      RHMAX = AR(3)
      KMAX = 3
   20 IF (AR(4).LT.RHMAX) GO TO 30
      RHMAX = AR(4)
      KMAX = 4
   30 IF (AR(5).LT.RHMAX) GO TO 40
      RHMAX = AR(5)
      KMAX = 5
   40 CALL SWITCH (KMAX,1, V, COORD, DXX)
C
C COMPUTE THE HCH ANGLES AND STORE RESULTS IN ARRAY AG
C
      CALL ANGLE
C
      DO 50 I = 1, 15
         DER(I) = 0.0D0
   50 CONTINUE
C
C  COMPUTE THE MORSE RANGE PARAMETERS ALPH:
C
      CALL CALPHA
C
C COMPUTE THE FOUR TRIATOMIC TERMS
C
      CALL TRI (AR(1),AR(5),AR(6),DER(1),DER(5),DER(6),1,3,2,EE,DTOT)
      DER(1) = DER(1)+DTOT
      DER2 = DTOT
      DER3 = DTOT
      DER4 = DTOT
      CALL TRI (AR(2),AR(5),AR(7),DER(2),DR2,DER(7),1,3,2,E,DTOT)
      EE = EE+E
      DER(2) = DER(2)+DER2+DTOT
      DER(1) = DER(1)+DTOT
      DER3 = DER3+DTOT
      DER4 = DER4+DTOT
      DER(5) = DER(5)+DR2
      CALL TRI (AR(3),AR(5),AR(8),DER(3),DR2,DER(8),1,3,2,E,DTOT)
      EE = EE+E
      DER(3) = DER(3)+DER3+DTOT
      DER(1) = DER(1)+DTOT
      DER(2) = DER(2)+DTOT
      DER4 = DER4+DTOT
      DER(5) = DER(5)+DR2
      CALL TRI (AR(4),AR(5),AR(9),DER(4),DR2,DER(9),1,3,2,E,DTOT)
      EE = EE+E
      DER(4) = DER(4)+DER4+DTOT
      DER(1) = DER(1)+DTOT
      DER(2) = DER(2)+DTOT
      DER(3) = DER(3)+DTOT
      DER(5) = DER(5)+DR2
C
C COMPUTE FORCE CONSTANTS
C COMPUTE ATTENUATION TERMS
C
      XA = AR(1)-RE(1)
      XB = AR(2)-RE(1)
      XC = AR(3)-RE(1)
      XD = AR(4)-RE(1)
      XA2 = XA*XA
      XB2 = XB*XB
      XC2 = XC*XC
      XD2 = XD*XD
      XXA = AR(6)-RE(2)
      XXB = AR(7)-RE(2)
      XXC = AR(8)-RE(2)
      XXD = AR(9)-RE(2)
      EXA11 = EXP(-APARM(2)*R2(6))
      EXA12 = EXP(-APARM(2)*R2(7))
      EXA13 = EXP(-APARM(2)*R2(8))
      EXA14 = EXP(-APARM(2)*R2(9))
      EXA21 = EXP(-APARM(5)*XXA*XXA)
      EXA22 = EXP(-APARM(5)*XXB*XXB)
      EXA23 = EXP(-APARM(5)*XXC*XXC)
      EXA24 = EXP(-APARM(5)*XXD*XXD)
      A11 = 1.0D0-EXA11
      A12 = 1.0D0-EXA12
      A13 = 1.0D0-EXA13
      A14 = 1.0D0-EXA14
      A21 = APARM(3)+APARM(4)*EXA21
      A22 = APARM(3)+APARM(4)*EXA22
      A23 = APARM(3)+APARM(4)*EXA23
      A24 = APARM(3)+APARM(4)*EXA24
      EXF1 = EXP(-A21*XA2)
      EXF2 = EXP(-A22*XB2)
      EXF3 = EXP(-A23*XC2)
      EXF4 = EXP(-A24*XD2)
      F(1) = A11*EXF1
      F(2) = A12*EXF2
      F(3) = A13*EXF3
      F(4) = A14*EXF4
C
C DK'S ARE THE FORCE CONSTANTS
C     FIRST COMPUTE THE K0 SINCE IT IS NOW A FUNCTION OF AR(1)
C
      CALL FORCE (DADR1)
      DK(1) = APARM(1)*F(1)*F(2)
      DK(2) = APARM(1)*F(1)*F(3)
      DK(3) = APARM(1)*F(1)*F(4)
      DK(4) = APARM(1)*F(2)*F(3)
      DK(5) = APARM(1)*F(2)*F(4)
      DK(6) = APARM(1)*F(3)*F(4)
C
C COMPUTE THE EQUILIBRIUM ANGLES--THETA AND D(THETA)/DR
C
      CALL THETA (ASF,DASF,ASB,DASB)
      DAA = 0.0D0
C
C COMPUTE CONTRIBUTION TO THE TOTAL ENERGY FROM THE ANGLE TERMS
C
      SUM = 0.0D0
      IF1 = 1
      IF2 = 2
      IF3 = 3
      IB1 = 4
      IB2 = 5
      IB3 = 6
      DO 60 I = 1, 3
         O = AG(I)-ASF
         Y1(I) = O
         Y2(I) = O*O
         SUM = SUM+DK(I)*Y2(I)
   60 CONTINUE
      DO 70 I = 4, 6
         Y1(I) = AG(I)-ASB
         Y2(I) = Y1(I)*Y1(I)
         SUM = SUM+DK(I)*Y2(I)
   70 CONTINUE
      EE = EE+0.5D0*SUM
C
C COMPUTE DERIVATIVES OF THE ANGLE ATTENUATION TERMS WITH RESPECT TO R6,
C                       AND R9
C
      DFR(1) = 2.0D0*APARM(2)*AR(6)*EXA11*EXF1+A11*EXF1*XA2*2.0D0*
     *   APARM(4)*APARM(5)*XXA*EXA21
      DFR(2) = 2.0D0*APARM(2)*AR(7)*EXA12*EXF2+A12*EXF2*XB2*2.0D0*
     *   APARM(4)*APARM(5)*XXB*EXA22
      DFR(3) = 2.0D0*APARM(2)*AR(8)*EXA13*EXF3+A13*EXF3*XC2*2.0D0*
     *   APARM(4)*APARM(5)*XXC*EXA23
      DFR(4) = 2.0D0*APARM(2)*AR(9)*EXA14*EXF4+A14*EXF4*XD2*2.0D0*
     *   APARM(4)*APARM(5)*XXD*EXA24
C
C COMPUTE DERIVATIVES OF ATTENUATION TERMS WITH RESPECT TO R1,R2,R3, AND
C
      DF(1) = -2.0D0*A21*XA*F(1)
      DF(2) = -2.0D0*A22*XB*F(2)
      DF(3) = -2.0D0*A23*XC*F(3)
      DF(4) = -2.0D0*A24*XD*F(4)
C
C ADD IN CONTRIBUTION TO DER FROM DK/DR ANGLE TERMS AND D(THETA)/DR TERM
C FIRST THREE ANGLE TERMS
C
      DO 80 I = 1, 3
         STP = 0.5D0*DK(I)*Y2(I)
         W = DK(I)*Y1(I)
         IF (ABS(F(1)).GT.1.0D-10) DER(1)=DER(1)+STP*(DF(1)/F(1)+DADR1/
     *      APARM(1))+W*(-1.0D0/(AR(I+1)*ASS(I))+ACS(I)/(AR(1)*ASS(I)))
         IF (ABS(F(I+1)).GT.1.0D-10) DER(I+1)=DER(I+1)+STP*DF(I+1)/F(I+
     *      1)+W*(-1.0D0/(AR(1)*ASS(I))+ACS(I)/(ASS(I)*AR(I+1)))
         IF (ABS(F(1)).GT.1.0D-10) DER(6)=DER(6)+STP*DFR(1)/F(1)
         IF (ABS(F(I+1)).GT.1.0D-10) DER(I+6)=DER(I+6)+STP*DFR(I+1)/F(I
     *      +1)
C
         DER(I+9) = DER(I+9)+W*AR(I+9)/(AR(1)*AR(I+1)*ASS(I))
   80 CONTINUE
C
C FOURTH AND FIFTH ANGLE TERMS
C
      DO 90 I = 4, 5
         STP = 0.5D0*DK(I)*Y2(I)
         W = DK(I)*Y1(I)
         DER(1) = DER(1)+STP*DADR1/APARM(1)
         DER(2) = DER(2)+STP*DF(2)/F(2)+W*(-1.0D0/(ASS(I)*AR(I-1))+
     *      ACS(I)/(ASS(I)*AR(2)))
         DER(7) = DER(7)+STP*DFR(2)/F(2)
         IF (ABS(F(I-1)).GT.1.0D-10) DER(I-1)=DER(I-1)+STP*DF(I-1)/F(I-
     *      1)+W*(-1.0D0/(ASS(I)*AR(2))+ACS(I)/(ASS(I)*AR(I-1)))
C
         IF (ABS(F(I-1)).GT.1.0D-10) DER(I+4)=DER(I+4)+STP*DFR(I-1)/F(I
     *      -1)
C
         DER(I+9) = DER(I+9)+W*AR(I+9)/(AR(2)*AR(I-1)*ASS(I))
   90 CONTINUE
C
C CONTRIBUTION FROM SIXTH ANGLE TERM
C
      STP = 0.5D0*DK(6)*Y2(6)
      W = DK(6)*Y1(6)
      DER(1) = DER(1)+STP*DADR1/APARM(1)
      DER(3) = DER(3)+STP*DF(3)/F(3)+W*(-1.0D0/(ASS(6)*AR(4))+ACS(6)/
     *   (ASS(6)*AR(3)))
      DER(8) = DER(8)+STP*DFR(3)/F(3)
      DER(4) = DER(4)+STP*DF(4)/F(4)+W*(-1.0D0/(ASS(6)*AR(3))+ACS(6)/
     *   (ASS(6)*AR(4)))
      DER(9) = DER(9)+STP*DFR(4)/F(4)
      DER(15) = DER(15)+W*AR(15)/(ASS(6)*AR(3)*AR(4))
C
C CONTRIBUTION FROM EQ. ANGLE VARIATION TO D(ANGLE)/DR
C
      DO 100 KLM = 1, 4
         DER(KLM) = DER(KLM)-(DK(IF1)*Y1(IF1)+DK(IF2)*Y1(IF2)+DK(IF3)*Y1
     *      (IF3))*DASF-(DK(IB1)*Y1(IB1)+DK(IB2)*Y1(IB2)+DK(IB3)*Y1(IB3)
     *      )*DASB
  100 CONTINUE
      DER(5) = DER(5)-(DK(IF1)*Y1(IF1)+DK(IF2)*Y1(IF2)+DK(IF3)*Y1(IF3))*
     *   DAA-(DK(IB1)*Y1(IB1)+DK(IB2)*Y1(IB2)+DK(IB3)*Y1(IB3))*3.0D0*   
     *   SIN(ASF)*COS(ASF)*DAA/SIN(ASB)
C
      CALL SWITCH (KMAX,2, V, COORD, DXX)
C
      ENERGY = EE/EVCON+REFV
      DO 110 I = 1, 15
         DER(I) = DER(I)/EVCON
  110 CONTINUE
C
C      CALCULATE DRDX ARRAY
C
      CALL DRDXC
C
C      CALCULATE DEDX = DEDR(DER) * DRDX
C
      CALL MULT (DER,1,15,DRDX,15,18,DEDX)
      IPT = 1
      DO 120 I = 1, 6
         DXX(IPT) = DEDX(I)
         DXX(IPT+1) = DEDX(I+6)
         DXX(IPT+2) = DEDX(I+12)
         IPT = IPT+3
  120 CONTINUE
C
C     NOW ADD IN THE OUT OF PLANE BEND TERM
C
      CALL OOPB (VOOP,DVOOP,V,COORD,DXX)
      ENERGY = ENERGY+VOOP
      DO 130 I = 1, 18
         DXX(I) = DXX(I)+DVOOP(I)
  130 CONTINUE
C
      ENGYGS = ENERGY
      CALL EUNITZERO
      IF(NDER.NE.0) THEN
         do i=1,18
            DEGSDR(i)=DXX(i)
         enddo
         CALL RTOCART
         CALL DEDCOU
      ENDIF

C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE ANGLE
C
C ROUTINE TO COMPUTE CH4 ANGLES FROM THE INTERPARTICLE DISTANCES
C***********************************************************************
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /POTXCM/ XX(6),Y(6),Z(6),ENERGY,DRDX(15,18),DEDX(18)
      COMMON /POTRCM/ AG(6),R2(15),ACS(6),ASS(6)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
      DO 10 I = 1, 3
         X = (R2(1)+R2(I+1)-R2(I+9))/(2.0D0*AR(1)*AR(I+1))
         ACS(I) = X
         ASS(I) = SQRT(1.0D0-X*X)
         AG(I) = ACOS(X)
   10 CONTINUE
      DO 20 I = 1, 2
         X = (R2(2)+R2(I+2)-R2(I+12))/(2.0D0*AR(2)*AR(I+2))
         ACS(I+3) = X
         ASS(I+3) = SQRT(1.0D0-X*X)
         AG(I+3) = ACOS(X)
   20 CONTINUE
      X = (R2(3)+R2(4)-R2(15))/(2.0D0*AR(3)*AR(4))
      ACS(6) = X
      ASS(6) = SQRT(1.0D0-X*X)
      AG(6) = ACOS(X)
      RETURN
      END
C
C
      SUBROUTINE CALPHA
C
C     COMPUTE THE C-H AND C-H' MORSE RANGE PARAMETERS FOR THE SINGLET
C     AND TRIPLET ENERGY CURVES.  THESE WILL BE A FUNCTION OF THE
C     AVERAGE C-H BOND LENGTH IN METHANE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
      COMMON /PDT3CM/ FK0,FA,CA,CB
      COMMON /ALPDCM/ DALP
C
C     FIRST COMPUTE THE AVERAGE METHANE BOND LENGTH
C
      RAVG = (AR(1)+AR(2)+AR(3)+AR(4))/4.0D0
C
      HTAN = TANH(C1*(RAVG-RE(1)))
      ALPH(1) = CA+CB*((HTAN+1.0D0)/2.0D0)
      ALPH(3) = ALPH(1)
      BETA(1) = ALPH(1)
      BETA(3) = ALPH(1)
C
C     NOW THE DERIVATIVES---DALPH WILL BE THE DERIVATIVE OF ALPHA WITH
C     RESPECT TO R1,R2,R3 OR R4
C
      DALP = (C1*CB/8.0D0)*(1.0D0-HTAN**2)
C
      RETURN
      END
C
C***********************************************************************
C***********************************************************************
C
      SUBROUTINE DISTA
C
C ROUTINE TO COMPUTE THE INTERPARTICLE DISTANCES FROM THE CARTESIAN
C     COORDINATES OF THE ATOMS
C***********************************************************************
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      DIMENSION DF(15,3)
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /POTXCM/ X(6),Y(6),Z(6),ENERGY,DRDX(15,18),DEDX(18)
      COMMON /POTRCM/ AG(6),R2(15),ACS(6),ASS(6)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
      DO 10 I = 1, 4
         DF(I,1) = X(I)-X(5)
         DF(I,2) = Y(I)-Y(5)
         DF(I,3) = Z(I)-Z(5)
         R2(I) = DF(I,1)*DF(I,1)+DF(I,2)*DF(I,2)+DF(I,3)*DF(I,3)
         AR(I) = SQRT(R2(I))
   10 CONTINUE
      DF(5,1) = X(5)-X(6)
      DF(5,2) = Y(5)-Y(6)
      DF(5,3) = Z(5)-Z(6)
      R2(5) = DF(5,1)*DF(5,1)+DF(5,2)*DF(5,2)+DF(5,3)*DF(5,3)
      AR(5) = SQRT(R2(5))
      DO 20 I = 1, 4
         J = I+5
         DF(J,1) = X(I)-X(6)
         DF(J,2) = Y(I)-Y(6)
         DF(J,3) = Z(I)-Z(6)
         R2(J) = DF(J,1)*DF(J,1)+DF(J,2)*DF(J,2)+DF(J,3)*DF(J,3)
         AR(J) = SQRT(R2(J))
   20 CONTINUE
      DO 30 I = 1, 3
         J = I+1
         K = I+9
         DF(K,1) = X(1)-X(J)
         DF(K,2) = Y(1)-Y(J)
         DF(K,3) = Z(1)-Z(J)
         R2(K) = DF(K,1)*DF(K,1)+DF(K,2)*DF(K,2)+DF(K,3)*DF(K,3)
         AR(K) = SQRT(R2(K))
   30 CONTINUE
      DO 40 I = 1, 2
         J = I+2
         K = I+12
         DF(K,1) = X(2)-X(J)
         DF(K,2) = Y(2)-Y(J)
         DF(K,3) = Z(2)-Z(J)
         R2(K) = DF(K,1)*DF(K,1)+DF(K,2)*DF(K,2)+DF(K,3)*DF(K,3)
         AR(K) = SQRT(R2(K))
   40 CONTINUE
      DF(15,1) = X(3)-X(4)
      DF(15,2) = Y(3)-Y(4)
      DF(15,3) = Z(3)-Z(4)
      R2(15) = DF(15,1)*DF(15,1)+DF(15,2)*DF(15,2)+DF(15,3)*DF(15,3)
      AR(15) = SQRT(R2(15))
      RETURN
      END
      SUBROUTINE DRDXC
C
C THIS SUBROUTINE CALCULATES DRDX
C
C     ADDED BY FBB AT THE U. OF MINN.  1/24/84
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /POTXCM/ X(6),Y(6),Z(6),ENERGY,DRDX(15,18),DEDX(18)
      COMMON /POTRCM/ AG(6),R2(15),ACS(6),ASS(6)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
C
C
C ZERO OUT DRDX
C
      DO 10 I = 1, 18
         DO 10 J = 1, 15
            DRDX(J,I) = 0.0D0
   10 CONTINUE
C
C CALCULATE DRDX FOR C - H(1,2,3,4) BONDS
C
      DO 20 I = 1, 4
         DRDX(I,I) = (X(I)-X(5))/AR(I)
         DRDX(I,5) = -DRDX(I,I)
         DRDX(I,I+6) = (Y(I)-Y(5))/AR(I)
         DRDX(I,11) = -DRDX(I,I+6)
         DRDX(I,I+12) = (Z(I)-Z(5))/AR(I)
         DRDX(I,17) = -DRDX(I,I+12)
   20 CONTINUE
C
C CALCULATE DRDX FOR C - H' BOND
C
      DRDX(5,5) = (X(5)-X(6))/AR(5)
      DRDX(5,6) = -DRDX(5,5)
      DRDX(5,11) = (Y(5)-Y(6))/AR(5)
      DRDX(5,12) = -DRDX(5,11)
      DRDX(5,17) = (Z(5)-Z(6))/AR(5)
      DRDX(5,18) = -DRDX(5,17)
C
C CALCULATE DRDX FOR H - H' BONDS
C
      DO 30 I = 1, 4
         DRDX(5+I,I) = (X(I)-X(6))/AR(5+I)
         DRDX(5+I,6) = -DRDX(5+I,I)
         DRDX(5+I,I+6) = (Y(I)-Y(6))/AR(5+I)
         DRDX(5+I,12) = -DRDX(5+I,I+6)
         DRDX(5+I,I+12) = (Z(I)-Z(6))/AR(5+I)
         DRDX(5+I,18) = -DRDX(5+I,I+12)
   30 CONTINUE
C
C CALCULATE DRDX FOR H1 - H(2,3,4) BONDS
C
      DO 40 I = 1, 3
         DRDX(9+I,I+1) = -(X(1)-X(I+1))/AR(9+I)
         DRDX(9+I,1) = -DRDX(9+I,I+1)
         DRDX(9+I,I+7) = -(Y(1)-Y(I+1))/AR(9+I)
         DRDX(9+I,7) = -DRDX(9+I,I+7)
         DRDX(9+I,I+13) = -(Z(1)-Z(I+1))/AR(9+I)
         DRDX(9+I,13) = -DRDX(9+I,I+13)
   40 CONTINUE
C
C CALCULATE DRDX FOR H2 - H(3,4)
C
      DO 50 I = 1, 2
         DRDX(12+I,I+2) = -(X(2)-X(I+2))/AR(12+I)
         DRDX(12+I,2) = -DRDX(12+I,I+2)
         DRDX(12+I,I+8) = -(Y(2)-Y(I+2))/AR(12+I)
         DRDX(12+I,8) = -DRDX(12+I,I+8)
         DRDX(12+I,I+14) = -(Z(2)-Z(I+2))/AR(12+I)
         DRDX(12+I,14) = -DRDX(12+I,I+14)
   50 CONTINUE
C
C CALCULATE DRDX FOR H3 - H4
C
      DRDX(15,4) = -(X(3)-X(4))/AR(15)
      DRDX(15,3) = -DRDX(15,4)
      DRDX(15,10) = -(Y(3)-Y(4))/AR(15)
      DRDX(15,9) = -DRDX(15,10)
      DRDX(15,16) = -(Z(3)-Z(4))/AR(15)
      DRDX(15,15) = -DRDX(15,16)
C
      RETURN
      END
C
C
      SUBROUTINE FORCE (DADR1)
C
C     THIS SUBROUTINE COMPUTES THE BENDING FORCE CONSTANT AS A FUNCTION
C     OF AR(1). (APARM(1))   IT ALSO COMPUTES THE DERIVATIVE OF THIS
C     TERM WITH RESPECT TO AR(1).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
      COMMON /PDT3CM/ FK0,FA,CA,CB
C
      FEX = EXP(-B1*(AR(1)-RE(1))**2)
      APARM(1) = FK0+FA*FEX
      DADR1 = -2.0D0*FA*B1*(AR(1)-RE(1))*FEX
      RETURN
      END
      SUBROUTINE MULT (A,NA,MA,B,NB,MB,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C THIS SUBROUTINE PERFORMS THE MATRIX MULTIPLICATION :
C
C         C = A X B
C
C
      DIMENSION A(NA,MA),B(NB,MB),C(NA,MB)
C
      DO 30 I = 1, NA
C
         DO 20 J = 1, MB
C
            SUM = 0.0D0
            DO 10 K = 1, MA
C
               SUM = SUM+A(I,K)*B(K,J)
   10       CONTINUE
C
            C(I,J) = SUM
   20    CONTINUE
C
   30 CONTINUE
      RETURN
      END
      SUBROUTINE OOPB (VOOP,DVOOP,V,X,DX)
C
C     THIS SUBROUTINE CALCULTES THE OUT OF PLANE BEND POTENTIAL THAT
C     RESULTS FROM THE METHYL.  THE FORM USED IS THE SAME AS THAT
C     GIVEN BY DUCHOVIC ET. AL. (JPC 89, 1339, (1984)).
C     THE DERIVATIVES WITH RESPECT TO THE CARTESIAN COORDINATES (X)
C     ARE ALSO COMPUTED AND STORED IN DVOOP.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      DIMENSION DVOOP(18),XDIF(2),YDIF(2),ZDIF(2),RV(3,3),DRV(3,3,18)
      DIMENSION XMAGR(3),DXMAGR(3,18),DELTA(3),DDELTA(3,18)
      DIMENSION GAMMA(3),DGAMMA(3,18),TOP(3),DTOP(3,18),DCROSS(3,18)
      DIMENSION DGCROS(18),DPHI(18),DFDX(18),DHDX(18),DRDX(18)
C
C
      DIMENSION X(18),DX(18)
C
      COMMON /TESTCM/ GAMMA,DGAMMA
      COMMON /PDT2CM/ APHI,BPHI,CPHI,PCH4,FCH3,HCH3,RNOT,A3,B3,DELT, 
     *   DIJ,A0,B0,C0,D,REX
C
C     START BY COMPUTING THE NECESSARY DIFFERENCES IN THE R1,R2 AND
C     R3 VECTORS.
C
      XDIF(1) = X(7)-X(4)
      YDIF(1) = X(8)-X(5)
      ZDIF(1) = X(9)-X(6)
      XDIF(2) = X(10)-X(4)
      YDIF(2) = X(11)-X(5)
      ZDIF(2) = X(12)-X(6)
C
C     NOW COMPUTE THE INDIVIDUAL COMPONENTS OF THE R1,R2 AND R3 VECTORS
C
      J = 1
      DO 10 I = 1, 3
         RV(I,1) = X(J+3)-X(13)
         RV(I,2) = X(J+4)-X(14)
         RV(I,3) = X(J+5)-X(15)
         J = J+3
   10 CONTINUE
C
C     RSTAR IS THE C-H4 BOND LENGTH
C
      RSTAR = SQRT((X(1)-X(13))**2+(X(2)-X(14))**2+(X(3)-X(15))**2)
C
C     COMPUTE THE CROSS PRODUCT OF THESE DIFFERENCES
C
      CROSSX = YDIF(1)*ZDIF(2)-YDIF(2)*ZDIF(1)
      CROSSY = ZDIF(1)*XDIF(2)-XDIF(1)*ZDIF(2)
      CROSSZ = XDIF(1)*YDIF(2)-XDIF(2)*YDIF(1)
      GCROSS = SQRT(CROSSX**2+CROSSY**2+CROSSZ**2)
C
C     COMPUTE THE THREE TRIPLE PRODUCTS--WITH R1,R2 AND R3 (TOP)
C     AND THEN COMPUTE THE ANGLE GAMMA USED IN DELTA.
C
      DO 20 I = 1, 3
         XMAGR(I) = SQRT(RV(I,1)**2+RV(I,2)**2+RV(I,3)**2)
         TOP(I) = CROSSX*RV(I,1)+CROSSY*RV(I,2)+CROSSZ*RV(I,3)
         GAMMA(I) = +TOP(I)/(GCROSS*XMAGR(I))
   20 CONTINUE
C
C     NOW THE SECOND TERM IN DELTA--PHINOT
C     FIRST THE SWITCHING FUNCTION SPHI THEN PHINOT
C
      SPHI = 0.0D0
      REST = BPHI*(RSTAR-CPHI)**3
      IF (REST.LT.70.0D0) THEN
         EREST = EXP(REST)
         SPHI = 1.0D0-TANH(APHI*(RSTAR-RNOT)*EREST)
      ENDIF
      PHINOT = PCH4+(PCH4-PI/2.0D0)*(SPHI-1.0D0)
C
C     NOW PUT THESE TWO TERMS TOGETHER TO GET EACH DELTA
C
      DO 30 I = 1, 3
         DELTA(I) = ACOS(GAMMA(I))-PHINOT
   30 CONTINUE
C
C     NEXT COMPUTE THE PARTIALS OF DELTA W/R TO EACH X(I)
C      IN DOING THIS COMPUTATION FIRST FOR EACH TERM CALCULATED
C      ABOVE THE PARTIAL WILL BE COMPUTED AND A STORED USING THE
C      SAME VARIABLE WITH A 'D' IN FRONT.
C
C
C     START WITH THE PARTIAL OF CROSS WITH RESPECT TO EACH X(I)
C
      DO 60 I = 1, 3
         DO 40 J = 1, 3
            DCROSS(I,J) = 0.0D0
   40    CONTINUE
         DO 50 J = 13, 18
            DCROSS(I,J) = 0.0D0
   50    CONTINUE
         IF (I.EQ.1) THEN
            DCROSS(I,4) = 0.0D0
            DCROSS(I,5) = ZDIF(1)-ZDIF(2)
            DCROSS(I,6) = YDIF(2)-YDIF(1)
            DCROSS(I,7) = 0.0D0
            DCROSS(I,8) = ZDIF(2)
            DCROSS(I,9) = -YDIF(2)
            DCROSS(I,10) = 0.0D0
            DCROSS(I,11) = -ZDIF(1)
            DCROSS(I,12) = YDIF(1)
         ELSEIF (I.EQ.2) THEN
            DCROSS(I,4) = ZDIF(2)-ZDIF(1)
            DCROSS(I,5) = 0.0D0
            DCROSS(I,6) = XDIF(1)-XDIF(2)
            DCROSS(I,7) = -ZDIF(2)
            DCROSS(I,8) = 0.0D0
            DCROSS(I,9) = XDIF(2)
            DCROSS(I,10) = ZDIF(1)
            DCROSS(I,11) = 0.0D0
            DCROSS(I,12) = -XDIF(1)
         ELSEIF (I.EQ.3) THEN
            DCROSS(I,4) = YDIF(1)-YDIF(2)
            DCROSS(I,5) = XDIF(2)-XDIF(1)
            DCROSS(I,6) = 0.0D0
            DCROSS(I,7) = YDIF(2)
            DCROSS(I,8) = -XDIF(2)
            DCROSS(I,9) = 0.0D0
            DCROSS(I,10) = -YDIF(1)
            DCROSS(I,11) = XDIF(1)
            DCROSS(I,12) = 0.0D0
         ENDIF
   60 CONTINUE
C
C     NOW COMPUTE THE PARTIALS FOR THE THREE R VECTORS
C
      DO 70 I = 1, 3
         DO 70 J = 1, 3
C            DO 70 K = 1, 3
C BUGFIX:  DRV was not fully initialized to zero causing trouble on some
C          platforms.
            DO 70 K = 1, 18
               DRV(I,J,K) = 0.0D0
   70 CONTINUE
      DO 80 K = 4, 6
         DRV(1,K-3,K) = 1.0D0
   80 CONTINUE
      DO 90 K = 7, 9
         DRV(2,K-6,K) = 1.0D0
   90 CONTINUE
      DO 100 K = 10, 12
         DRV(3,K-9,K) = 1.0D0
  100 CONTINUE
      DO 110 I = 1, 3
         DO 110 J = 1, 3
            DRV(I,J,J+12) = -1.0D0
  110 CONTINUE
C
C     NOW COMPUTE THE PARTIALS FOR THE TRIPLE PRODUCT  (TOP)
C
      DO 120 I = 1, 18
         DO 120 J = 1, 3
            DTOP(J,I) = DCROSS(1,I)*RV(J,1)+CROSSX*DRV(J,1,I)+DCROSS(2,I
     *         )*RV(J,2)+CROSSY*DRV(J,2,I)+DCROSS(3,I)*RV(J,3)+CROSSZ*  
     *         DRV(J,3,I)
  120 CONTINUE
C
C     PARTIALS FOR THE MAGNITUDE OF THE CROSS PRODUCT  (GCROSS)
C
      DO 130 I = 1, 18
         DGCROS(I) = (1.0D0/GCROSS)*(DCROSS(1,I)*CROSSX+DCROSS(2,I)*    
     *      CROSSY+DCROSS(3,I)*CROSSZ)
  130 CONTINUE
C
C     PARTIALS FOR THE MAGNITUDE OF R1, R2 AND R3
C
      DO 150 I = 1, 3
         DO 150 J = 1, 18
            DXMAGR(I,J) = 0.0D0
            DO 140 K = 1, 3
               DXMAGR(I,J) = DXMAGR(I,J)+(1.0D0/XMAGR(I))*(RV(I,K)*DRV(I
     *            ,K,J))
  140       CONTINUE
  150 CONTINUE
C
C     NOW WE'RE READY FOR THE PARTIALS OF GAMMA
C
      DO 160 I = 1, 3
         DO 160 J = 1, 18
            DGAMMA(I,J) = +GCROSS*XMAGR(I)*DTOP(I,J)-TOP(I)*(XMAGR(I)*  
     *         DGCROS(J)+GCROSS*DXMAGR(I,J))
            DGAMMA(I,J) = DGAMMA(I,J)/(GCROSS*XMAGR(I))**2
            DGAMMA(I,J) = DGAMMA(I,J)*(-1.0D0)/(SQRT(1.0D0-GAMMA(I)**2)
     *         )
  160 CONTINUE
C
C     NOW THE PARTIAL FOR THE SECOND TERM IN DELTA -- PHINOT
C       FIRST COMPUTE THE PARTIAL OF S(SWITCHING FUNCTION) WITH
C       RESPECT TO RSTAR, THE THE PARTIAL OF RSTAR W/R X(I) AND
C       FINALLY THE PARTIAL OF PHINOT W/R TO S.
C
      DSDR = 0.0D0
      REST = BPHI*(RSTAR-CPHI)**3
      IF (REST.LT.70.0D0) THEN
         TEXP = EXP(REST)
         DSDR = -((SECH(APHI*(RSTAR-RNOT)*TEXP))**2)
         DSDR = DSDR*APHI*TEXP*(1.0D0+(RSTAR-RNOT)*BPHI*3.0D0*(RSTAR-   
     *      CPHI)**2)
      ENDIF
      DO 170 I = 1, 18
         DRDX(I) = 0.0D0
         DPHI(I) = 0.0D0
  170 CONTINUE
C
      CST = 1.0D0/(2.0D0*RSTAR)
      DO 180 I = 1, 3
         DRDX(I) = CST*2.0D0*(X(I)-X(I+12))
         DRDX(I+12) = -DRDX(I)
         DPHI(I) = DRDX(I)*DSDR*(PCH4-PI/2.0D0)
         DPHI(I+12) = -DPHI(I)
  180 CONTINUE
C
C     NOW PUT BOTH TERMS TOGETHER TO GET THE PARTIALS OF DELTA
C
      DO 190 I = 1, 3
         DO 190 J = 1, 18
            DDELTA(I,J) = DGAMMA(I,J)-DPHI(J)
  190 CONTINUE
C
C     WE HAVE NOW COMPUTED ALL THE DELTA TERMS.  NEXT COMPUTE THE
C       FORCE CONSTANTS (F AND H) AND THEIR DERIVATIVES.  THE
C       OUT-OF-PLANE BEND ENERGY WILL ALSO BE COMPUTED HERE (VOOP)
C
C     FIRST THE SWITCHING FUNCTION
C
      S3 = 1.0D0-TANH(A3*(RSTAR-RNOT)*(RSTAR-B3)**2)
      F = (1.0D0-S3)*FCH3
      H = (1.0D0-S3)*HCH3
      VOOP = 0.0D0
      DO 200 I = 1, 3
         VOOP = VOOP+F*(DELTA(I)**2)+H*(DELTA(I)**4)
  200 CONTINUE
C
C     PARTIAL OF THE SWITCHING FUNCTION W/R TO RSTAR
C
      DSDR = -(SECH(A3*(RSTAR-RNOT)*(RSTAR-B3)**2)**2)*A3*((RSTAR-B3)**2
     *   +(RSTAR-RNOT)*2*(RSTAR-B3))
      DFDS = -FCH3
      DHDS = -HCH3
C
C     PARTIAL OF EACH FORCE CONSTANT
C
      DO 210 I = 1, 18
         DFDX(I) = 0.0D0
         DHDX(I) = 0.0D0
  210 CONTINUE
      DO 220 I = 1, 3
         DFDX(I) = DFDS*DSDR*DRDX(I)
         DHDX(I) = DHDS*DSDR*DRDX(I)
         DFDX(I+12) = -DFDX(I)
         DHDX(I+12) = -DHDX(I)
  220 CONTINUE
C
C    NOW WE ARE READY FOR THE FULL PARTIAL DERIVATIVES OF VOOP
C
                                                                        
      DO 250 I = 1, 18
         TERM1 = 0.0D0
         DO 230 J = 1, 3
            TERM1 = TERM1+DFDX(I)*DELTA(J)**2
  230    CONTINUE
         TERM2 = 2.0D0*DELTA(1)*DDELTA(1,I)+2.0D0*DELTA(2)*DDELTA(2,I)
     *           +2.0D0*DELTA(3)*DDELTA(3,I)
         TERM2 = TERM2*F
         TERM3 = 0.0D0
         DO 240 J = 1, 3
            TERM3 = TERM3+DHDX(I)*DELTA(J)**4
  240    CONTINUE
         TERM4 = 4.0D0*DELTA(1)**3*DDELTA(1,I)+4.0D0*DELTA(2)**3*
     *           DDELTA(2,I)+4.0D0*DELTA(3)**3*DDELTA(3,I)
         TERM4 = TERM4*H
         DVOOP(I) = TERM1+TERM2+TERM3+TERM4
  250 CONTINUE
      RETURN
      END
C
      FUNCTION SECH (X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (X.GT.70.0D0) THEN
         SECH = 0.0D0
      ELSE
         SECH = 1.0D0/COSH(X)
      ENDIF
      RETURN
      END
      SUBROUTINE SWITCH (KMAX,LMG, V, COORD, DXX)
C
C **********************************************************************
C THIS ROUTINE SWITCHES MATRICES R,R2,DER TO ALLOW FOR DIFFERENT H ATOMS
C      BEING OFF THE CH4 MOLECULE.
C LMG IS A CONTROL PARAMETER
C KMAX IS THE LABEL OF THE H ATOM WITH MAXIMUM C-H DISTANCE.
C IF LMG=1, R AND R2 MATRICES ARE SWITCHED APPROPRIATELY.
C IF LMG=2, R AND R2 MATRICES ARE SWITCHED BACK AND MATRIX DER IS SWITCH
C **********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      DIMENSION COORD(18),DXX(18)
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /POTXCM/ X(6),Y(6),Z(6),ENERGY,DRDX(15,18),DEDX(18)
      COMMON /POTRCM/ AG(6),R2(15),ACS(6),ASS(6)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
      GO TO (10,20,30,40,50), KMAX
C
C H1 IS MAXIMUM
C
   10 K = 1
      KK = 5
      L = 7
      LL = 10
      M = 8
      MM = 11
      N = 9
      NN = 12
      GO TO 60
C
C H2 IS MAXIMUM
C
   20 K = 2
      KK = 5
      L = 6
      LL = 10
      M = 8
      MM = 13
      N = 9
      NN = 14
      GO TO 60
C
C H3 IS MAXIMUM
C
   30 K = 3
      KK = 5
      L = 6
      LL = 11
      M = 7
      MM = 13
      N = 9
      NN = 15
      GO TO 60
C
C H4 IS MAXIMUM
C
   40 K = 4
      KK = 5
      L = 6
      LL = 12
      M = 7
      MM = 14
      N = 8
      NN = 15
      GO TO 60
C
C H5 IS MAXIMUM
C NO SWITCHES ARE REQUIRED
C
   50 GO TO 80
C
C SWITCH R AND R2 MATRICES
C
   60 STORE = AR(K)
      AR(K) = AR(KK)
      AR(KK) = STORE
      STORE = AR(L)
      AR(L) = AR(LL)
      AR(LL) = STORE
      STORE = AR(M)
      AR(M) = AR(MM)
      AR(MM) = STORE
      STORE = AR(N)
      AR(N) = AR(NN)
      AR(NN) = STORE
      R2(K) = AR(K)*AR(K)
      R2(KK) = AR(KK)*AR(KK)
      R2(L) = AR(L)*AR(L)
      R2(LL) = AR(LL)*AR(LL)
      R2(M) = AR(M)*AR(M)
      R2(MM) = AR(MM)*AR(MM)
      R2(N) = AR(N)*AR(N)
      R2(NN) = AR(NN)*AR(NN)
      IF (LMG-1) 80, 80, 70
C
C SWITCH DER MATRIX *******
C
   70 STORE = DER(K)
      DER(K) = DER(KK)
      DER(KK) = STORE
      STORE = DER(L)
      DER(L) = DER(LL)
      DER(LL) = STORE
      STORE = DER(M)
      DER(M) = DER(MM)
      DER(MM) = STORE
      STORE = DER(N)
      DER(N) = DER(NN)
      DER(NN) = STORE
   80 RETURN
      END
                                                                        
      SUBROUTINE THETA (ASF,DASF,ASB,DASB)
C
C     THIS SUBROUTINE COMPUTES THE EQUILIBRIUM THETA AS A FUNCTION OF
C     R (THE AVERAGE C-H BOND LENGTH IN METHANE).
C     THE DERIVATIVE OF THETA WITH RESPECT TO AR(I) IS ALSO COMPUTED.
C     NOTE HERE THAT BECAUSE OF SYMMETRY EACH OF THESE PARTIALS IS
C     THE SAME.       ADDED 11/22/85
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
      COMMON /PDT2CM/ APHI,BPHI,CPHI,PCH4,FCH3,HCH3,RNOT,A3,B3,DELTA,
     *   DIJ,A0,B0,C0,D0,REX
C
      RAVG = (AR(1)+AR(2)+AR(3)+AR(4))/4.0D0
C
C     COMPUTE THETA
C
      IF (RAVG.GT.DIJ) THEN
         ASF = PI/2.0D0
         DASF = 0.0D0
      ELSEIF (RAVG.LE.RE(1)) THEN
         ASF = TAU
         DASF = 0.0D0
      ELSEIF (RAVG.GT.RE(1).AND.RAVG.LE.RE(1)+DELTA) THEN
         DELTA2 = DELTA*DELTA
         FK1 = -4.0D0*CP/(DELTA2)
         A = 3.0D0*FK1/(DELTA2)
         B = -(5.0D0*A*DELTA2+FK1)/(2.0D0*DELTA)
         C = -(20.0D0*A*DELTA2+12.0D0*B*DELTA)/6.0D0
         F = TAU
         RDIF = RAVG-RE(1)
         RDIF2 = RDIF*RDIF
         RDIF3 = RDIF2*RDIF
         RDIF4 = RDIF3*RDIF
         RDIF5 = RDIF4*RDIF
         ASF = A*RDIF5+B*RDIF4+C*RDIF3+F
         DASF = 5.0D0*A*RDIF4+4.0D0*B*RDIF3+3.0D0*C*RDIF2
      ELSEIF (RAVG.LE.DIJ-DELTA) THEN
         ASF = TAU-4.0D0*CP*(RAVG-RE(1))
         DASF = -4.0D0*CP
      ELSEIF (RAVG.LE.DIJ.AND.RAVG.GT.DIJ-DELTA) THEN
         FK1 = (TAU-4.0D0*CP*(DIJ-DELTA-RE(1))-F)/(DELTA**3)
         FK2 = (-4.0D0*CP)/(DELTA**2)
         A = -3.0D0*(2.0D0*FK1+FK2)/(DELTA**2)
         B = (FK2+5.0D0*A*DELTA**2)/(2.0D0*DELTA)
         C = (12.0D0*B*DELTA-20.0D0*A*DELTA**2)/6.0D0
         F = PI/2.0D0
         RDIF = RAVG-DIJ
         RDIF3 = RDIF*RDIF*RDIF
         RDIF4 = RDIF3*RDIF
         RDIF5 = RDIF4*RDIF
         ASF = A*RDIF5+B*RDIF4+C*RDIF3+F
         DASF = 5.0D0*A*RDIF4+4.0D0*B*RDIF3+3.0D0*C*RDIF*RDIF
      ENDIF
C
C COMPUTE BACKSIDE EQUILIBRIUM ANGLE
C
      XYZ = COS(ASF)
      CTAU = COS(TAU)
      CTAU2 = CTAU*CTAU
      A = (CTAU+0.5D0)/CTAU2
      ASB = ACOS(A*XYZ*XYZ-0.5D0)
      DASB = 2.0D0*A*XYZ*SIN(ASF)*DASF/SIN(ASB)
C
C NOW COMPUTE THE DER W/R TO RI INSTEAD OF RAVG
C
      DASF = DASF/4.0D0
      DASB = DASB/4.0D0
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE TRI (RR1,RR2,RR3,DR1,DR2,DR3,I,J,K,E,DTOT)
C
C ROUTINE TO CALCULATE THREE-BODY ENERGY FOR ABC SYSTEM
C RR1,RR2,RR3 ARE THE THREE INTERPARTICLE DISTANCES
C DR1,DR2,DR3 ARE THE THREE DERIVATIVES WITH RESPECT TO RR1,RR2, AND RR3
C I,J, AND K GIVE THE ARRAY NUMBER FOR THE APPROPRIATE POTENTIAL PARAMET
C E IS THE FINAL THREE-BODY ENERGY
C***********************************************************************
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /POTXCM/ X(6),Y(6),Z(6),ENERGY,DRDX(15,18),DEDX(18)
      COMMON /POTRCM/ AG(6),R2(15),ACS(6),ASS(6)
      COMMON /DERVCP/ DECCR1,DECCR5,DECCR6
      COMMON /ALPDCM/ DALP
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,BP1,CP1
      A = RR1-RE(I)
      B = RR2-RE(J)
      C = RR3-RE(K)
      A1 = EXP(-ALPH(I)*A)
      B1 = EXP(-ALPH(J)*B)
      C1 = EXP(-ALPH(K)*C)
      A2 = A1*A1
      B2 = B1*B1
      C2 = C1*C1
      E1AB = D1(I)*(A2-2.0D0*A1)
      E1AC = D1(J)*(B2-2.0D0*B1)
      E1BC = D1(K)*(C2-2.0D0*C1)
      DE1AB = 2.0D0*D1(I)*(A1-A2)*ALPH(I)
      DAB1 = 2.0D0*D1(I)*(A1-A2)*(A*DALP)
      DE1AC = 2.0D0*D1(J)*(B1-B2)*ALPH(J)
      DAC1 = 2.0D0*D1(J)*(B1-B2)*(B*DALP)
      DE1BC = 2.0D0*ALPH(K)*D1(K)*(C1-C2)
      AA1 = EXP(-BETA(I)*A)
      AA2 = AA1*AA1
      E3AB = D3(I)*(AA2+2.0D0*AA1)
      DE3AB = -2.0D0*D3(I)*(AA1+AA2)*BETA(I)
      DAB3 = -2.0D0*D3(I)*(AA1+AA2)*A*DALP
      BB1 = EXP(-BETA(J)*B)
      BB2 = BB1*BB1
      E3AC = D3(J)*(BB2+2.0D0*BB1)
      DE3AC = -2.0D0*D3(J)*(BB1+BB2)*BETA(J)
      DAC3 = -2.0D0*D3(J)*(BB1+BB2)*B*DALP
      CC1 = EXP(-BETA(K)*C)
      CC2 = CC1*CC1
      E3BC = D3(K)*(CC2+2.0D0*CC1)
      DE3BC = -2.0D0*BETA(K)*D3(K)*(CC1+CC2)
      QAB = (E1AB+E3AB)/2.0D0
      QAC = (E1AC+E3AC)/2.0D0
      QBC = (E1BC+E3BC)/2.0D0
      ALAB = (E1AB-E3AB)/2.0D0
      ALAC = (E1AC-E3AC)/2.0D0
      ALBC = (E1BC-E3BC)/2.0D0
      XX = ALAB-ALBC
      YY = ALBC-ALAC
      ZZ = ALAC-ALAB
      U = (XX*XX+YY*YY+ZZ*ZZ)/2.0D0
      U = SQRT(U)
      E = QAB+QBC+QAC-U
      BB = 4.0D0*U
      DXR = (DAB1-DAB3)/2.0D0
      DYR = (DAC3-DAC1)/2.0D0
      DZR = -(DXR+DYR)
      DUX = XX/(2.0D0*U)
      DUY = YY/(2.0D0*U)
      DUZ = ZZ/(2.0D0*U)
      DUR = DUX*DXR+DUY*DYR+DUZ*DZR
      DQR = (DAB1+DAC1+DAB3+DAC3)/2.0D0
      DTOT = DQR-DUR
      DR1 = (DE1AB+DE3AB)/2.0D0-(2.0D0*ALAB-ALBC-ALAC)*(DE1AB-DE3AB)/BB
      DR2 = (DE1AC+DE3AC)/2.0D0-(2.0D0*ALAC-ALAB-ALBC)*(DE1AC-DE3AC)/BB
      DR3 = (DE1BC+DE3BC)/2.0D0-(2.0D0*ALBC-ALAC-ALAB)*(DE1BC-DE3BC)/BB
      RETURN
      END
      SUBROUTINE TRIPLE
C
C     THIS SUBROUTINE COMPUTES THE C-H' SATO PARAMETER WHICH IS A
C     FUNCTION OF AR(H'-H) AND THE C-H-H' BEND ANGLE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,CP1
      COMMON /PDT2CM/ APHI,BPHI,CPHI,PCH4,FCH3,HCH3,RNOT,A3,B3,DELT, 
     *   DIJ,A,B,C0,D,REX
      COMMON /DERVCP/ DECCR1,DECCR5,DECCR6
C
C     COMPUTE C1--IT IS A FUNCTION OF R6--THE H-H' BOND DISTANCE
C
      IF (AR(6).GT.15.0D0) THEN
         R6 = AR(6)/3.0D0
         YX = (EXP(2.0D0*PI*(R6-REX)/D))**3
      ELSE
         YX = EXP(2.0D0*PI*(AR(6)-REX)/D)
      ENDIF
      C1 = -(A*YX/(YX+1.0D0)+B/YX)-C0
      DYXR = (2.0D0*PI/D)*YX
      DCYX = B/YX-(A*YX/(1+YX))/(1+YX)
      DC1R6 = DCYX*DYXR/YX
C
C     COMPUTE THE C-H4-H5 ANGLE
C
      TEM = AR(5)**2-AR(1)**2-AR(6)**2
      TEMBOT = -2.0D0*AR(1)*AR(6)
      CHECK = TEM/TEMBOT
      IF (CHECK.LE.-1.0D0) THEN
         ALPHA = 0.0D0
         CC(3) = C0
         DECCR1 = 0.0D0
         DECCR5 = 0.0D0
         DECCR6 = 0.0D0
         RETURN
      ELSE
         ALPHA = ACOS(TEM/(-2.0D0*AR(1)*AR(6)))
         ALPHA = PI-ALPHA
      ENDIF
C
C     NOW COMPUTE THE ANGLE DERIVATIVES (DERIV. ALPHA W.R. R1,R5,R6)
C
      U = TEM/(-2.0D0*AR(1)*AR(6))
      UI = -1.0D0/(SQRT(1.0D0-U**2))
      DALR5 = -UI*AR(5)/(AR(1)*AR(6))
      TEM2 = 1.0D0/(-2.0D0*AR(1)*AR(6))
      TEM22 = TEM2*TEM2
      DALR1 = -2.0D0*UI*(AR(1)*TEM2-AR(6)*TEM*TEM22)
      DALR6 = -2.0D0*UI*(AR(6)*TEM2-AR(1)*TEM*TEM22)
C
C     COMBINE THE RADIAL AND ANGLE COMPONENT TO COMPUTE THE TRIPLET
C     PARAMETER CC(3) AND THE DERIVATIVES.  HERE TEM IS JUST THE PARTIAL
C     OF CC W.R.T. ALPHA.
C
      DEG20 = PI/9.0D0
      IF (ALPHA.LE.DEG20) THEN
         CC(3) = C0+C1*(SIN(9.0D0*ALPHA/2.0D0))**4
         TEM = 18.0D0*(SIN(9.0D0*ALPHA/2.0D0))**3*COS(9.0D0*ALPHA/    
     *      2.0D0)
         DECCR1 = C1*TEM*DALR1
         DECCR5 = C1*TEM*DALR5
         DECCR6 = DC1R6*(SIN(9.0D0*ALPHA/2.0D0))**4+C1*TEM*DALR6
      ELSE
         CC(3) = C0+C1
         DECCR1 = 0.0D0
         DECCR5 = 0.0D0
         DECCR6 = DC1R6
      ENDIF
      RETURN
      END
C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   N3TMMN = 3 * NATOMS                                                         
C   NATOMS = the number of atoms represented by this potential function         
C                                                                               
C   The variable N3TMMN is the minimum value of N3TM allowed to be              
C   passed by the calling routine for the number of cartesian                   
C   coordinates needed to represent the full system represented by this         
C   potential energy surface routine.                                           
C   N3TM must be greater than or equal to N3TMMN.                               
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER(N3ATOM = 75)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
      PARAMETER (N3TMMN = 18)                                                   
C                                                                               
      COMMON/PT3CM/ EZERO(ISURF+1)                                              
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      DIMENSION DFR(4),DVOOP(18),Y1(6),Y2(6),F(4),DF(4),DK(6)
C
      DIMENSION COORD(18),DXX(18)
C
      COMMON /ICORCM/ AR(15),DER(15)
C
      COMMON /POTXCM/ X(6),Y(6),Z(6),ENERGY,DRDX(15,18),DEDX(18)
      COMMON /POTRCM/ AG(6),R2(15),ACS(6),ASS(6)
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),
     *   APARM(5),REFV,TAU,CP,B1,C1
      COMMON /DERVCP/ DECCR1,DECCR5,DECCR6
      COMMON /TEMCOM/ RTEM(9),DTEMP(9),ETEM
      COMMON /ALPDCM/ DALP
      DATA EVCON / 2.721161D+1 /
C
C  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM         
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='T. Joseph, R. Steckler, and D. G. Truhlar'
       REF(2)='J. Chem. Phys. Vol. 87, p. 7036, 1987'
C 
      INDEXES(1) = 1                                                            
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
      INDEXES(4) = 1                                                            
      INDEXES(5) = 6                                                            
      INDEXES(6) = 1                                                            
C                                                                               
C                                                                               
C                                                                               
      IRCTNT=6                                                                  
C                                                                               
C      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C
      BLOCK DATA PTPACM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
C
      COMMON/PT3CM/ EZERO(ISURF+1)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /PDATCM/ D1(3),D3(3),ALPH(3),RE(3),BETA(3),CC(3),AA(3),    
     *   APARM(5),REFV,TAU,CP,B1,C1
      COMMON /PDT2CM/ APHI,BPHI,CPHI,PCH4,FCH3,HCH3,RNOT,A3,B3,DELT, 
     *   DIJ,A,B,C0,D,REX
      COMMON /PDT3CM/ FK0,FA,CA,CB
C
      DATA NASURF /1,35*0/
      DATA NDER /1/
       DATA NFLAG /1,1,15*0,6,0,0/
C
      DATA ANUZERO /0.0D0/
      DATA ICARTR,MSURF,MDER/1,0,1/
      DATA NULBL /25*0/
      DATA NATOMS /6/
C
      DATA D1 /  4.8668D0,4.74660D0,1.1452D0 /
      DATA D3 /  1.6890D0,1.75200D0,0.8870D0 /
      DATA B1 / 3.00D0 /
      DATA C1 / 11.00D0 /
      DATA ALPH / 0.939D0,1.029955D0,0.939D0 /
      DATA RE / 2.0673D0,1.402D0,2.0673D0 /
      DATA BETA / 0.939D0,1.029955D0,0.939D0 /
      DATA AA / 0.0D0,0.0D0,-2.3603D0 /
      DATA REFV / 0.71098329D0 /
      DATA TAU / 1.9105588D0 /
      DATA CP / 0.080768D0 /
      DATA APARM / 4.0D0,0.8999975D0,0.448035D0,0.606528D0,3.239921D0 /
C
      DATA APHI / 2.7982366D-1 /
      DATA BPHI / 5.937217D-2 /
      DATA CPHI / 3.6301534D0 /
      DATA PCH4 / 1.9106D0 /
      DATA FCH3 / 2.19D-2 /
      DATA HCH3 / 4.38D-02 /
      DATA RNOT / 2.052D0 /
      DATA A3 / 2.102956D-2 /
      DATA B3 / 5.798533D-1 /
      DATA DELT / 0.1D0 /
      DATA DIJ / 3.1190D0 /
      DATA A / -1.56944D+3 /
      DATA B / -1.12008D+3 /
      DATA C0 / 1.5786151D+3 /
      DATA D / 2.67146D0 /
      DATA CA / 0.903D0 /
      DATA CB / 0.072D0 /
      DATA REX / 1.18584D0 /
      DATA FK0 / 2.55D0 /
      DATA FA / 0.8D0 /
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
