      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      COMMON /PT1CM/ R(3),ENERGY,DER(3)
      COMMON /PT2CM/ NSURF, NDER, NDUM(21)
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
      R(1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      R(2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      R(3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211
      
      NSURF=0
      NDER=igrad
      DER=0.d0     

      if(first_time_data) then
      call prepot
      first_time_data=.false.
      endif

      call pot

      ENERGY=ENERGY*27.211386
      DER=DER*51.422067

      r2=r*0.529177211
      call evdrdx(tx, r2, drdx)
      dedx=0.d0
      do i=1,9
      do j=1,3
        dedx(i)=dedx(i)+DER(j)*drdx(j,i)
      enddo
      enddo


      if (igrad==0) then
        do istate=1,nstates
          p(istate)=ENERGY
        enddo
      else if (igrad==1) then
        do istate=1,nstates
          p(istate)=ENERGY
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
      SUBROUTINE PREPOT
C
C   System:    HO2
C   Reference: C. F. Melius and R. J. Blint
C              Chem. Phys. Lett. 64, 183 (1979).
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in DATA statements.
C   and in the block data subprogram PPARMS.
C
C   Common Blocks (used between the calling program and this potential):
C      /PT1CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, ground state electronic energy, and
C        derivatives of the ground electronic state energy with respect
C        to the coordinates.
C      /PT2CM/ NSURF, NDER, NDUM(21)
C        passes the control flags where
C        NSURF = 0    ground electronic surface
C        NDER  = 0 => no derivatives are computed
C              = 1 => derivatives of the energy for the ground electronic
C                     state with respect to the coordinates are computed
C        NDUM  - not used
C      /PT4CM/ IPTPRT, IDUM(19)
C        IPTPRT passes the FORTRAN unit number used for potential output
C        IDUM   not used
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to
C        Re(AB), the equilibrium bond length for the AB diatomic.
C        Similarly, the terms EASYBC and EASYAC represent the energies
C        in the BC and the AC valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NSURF            0
C      NDER             1
C      IPTPRT           6
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         GENERAL INFORMATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          THE POTENTIAL CONSISTS OF THREE MORSE TERMS FOR
C          THE THREE DIATOMICS, HO,OO, AND HO RESPECTIVELY
C          PLUS TWO INTERACTION TERMS.
C
C          EACH POTENTIAL TERM WILL BE CALCULATED INDIVIDUALLY
C          TO SIMPLIFY READING THE CODE.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          DEFINE VARIABLES, ARRAYS, AND COMMON BLOCK
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION NUM,DENOM,GAMMA,JUNK,EXPO,DEXPO,MEW
C
C
      COMMON /PT1CM/ R(3),ENERGY,DER(3)
      COMMON /PT2CM/ NSURF, NDER, NDUM(20)
      COMMON /PT4CM/ IPRT, IDUM(19)
      COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
C
      DIMENSION DCSTH(3),DV(3),GAMMA(3),JUNK(4),C(5,3,3),A(5),DA(5,3),  
     *   REQ(3),DE(3),EXPO(3),DEXPO(3),DMEW(3),B(5,3),DB(5,3,3),DPROD1(3
     *   ),DPROD2(3),DSEC(3),DTHIRD(3),DPROD3(3)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          Initialize the constants
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DATA C / 77.45D0, -0.4071D0, -0.508D0,
     *         0.1489D0, 1.013D0, -10.495D0,
     *         9.1484D0, 10.273D0, 0.0D0, 0.0D0,
     *         -3.050D0, -33.78D0, -22.951D0, 0.0D0,
     *         0.0D0, -1.3699D0, 0.4101D0, 0.3411D0,
     *         16.906D0, 0.0005637D0, -0.6359D0,
     *         -0.06253D0, -0.1225D0, 0.0D0, 0.0D0,
     *         -0.00906D0, 0.4435D0, 0.7439D0, 0.0D0, 0.0D0,
     *         -0.3498D0, 0.6617D0, 0.8677D0, 1.3858D0,
     *         233.36D0, -4.756D0, -1.626D0,
     *         0.6337D0, 0.0D0, 0.0D0, -476.32D0,
     *         -1.2104D0, -1.2762D0, 0.0D0, 0.0D0 /
C   The array DE contains the dissociation energies in hartree atomic units.
      DATA DE / 0.1559D0, 0.1779D0, 0.1559D0/
C   The array GAMMA contains the Morse Betas in reciprocal bohrs.
      DATA GAMMA / 1.2670D0, 1.4694D0, 1.2670D0/
C   The array REQ contains the equilibrium bond lengths in bohr.
      DATA REQ / 1.8460D0, 2.3158D0, 1.8460D0/
C   The variable ZEROAD contains the constant added to the energy in hartree 
C   atomic units.  
      DATA ZEROAD /-0.022000295721D0/
      SAVE
C
C          SET ALPHA, UNITS: RECIPROCAL BOHR
      ALPHA = 0.9172D0
C
C          SET BETA, UNITS: RECIPROCAL BOHR
      BETA = 1.4694D0
C
C   Initialize the potential values in the three asymptotic valleys.
      EASYAB = DE(1)
      EASYBC = DE(2)
      EASYAC = DE(3)
C
C   Echo potential paramters to unit IPRT.
C
C   Convert ZEROAD to eV and kcal/mol and echo to unit IPRT.
      ZEV = ZEROAD*27.21106D0
      ZKCAL = ZEROAD*627.5095D0
C
CCCCCCCCCCCC
C
      RETURN
C
CCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        MOVE ONTO THE MAIN PART OF THE ROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ENTRY POT
C
C   Check the values of NSURF and NDER
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         SET RAB,RBC,AND RAC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      RBC = R(2)
C
C           RAB IS DEFINED AS THE SMALLER OF R(1) AND R(3)
C           RAC IS DEFINED AS THE LARGER OF R(1) AND R(3)
C
      IF (R(1).GT.R(3)) GO TO 40
C
      RAB = R(1)
      RAC = R(3)
      GO TO 50
C
   40 RAB = R(3)
      RAC = R(1)
   50 CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CALCULATE USEFUL EXPONENTIALS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      EXPO(1) = EXP(-ALPHA*RAB)
      EXPO(2) = EXP(-BETA*RBC)
      EXPO(3) = EXP(-ALPHA*RAC)
C
C          AND THEIR DERIVATIVES
C
         IF (NDER .EQ. 1) THEN
             DEXPO(1) = -ALPHA*EXPO(1)
             DEXPO(2) = -BETA*EXPO(2)
             DEXPO(3) = -ALPHA*EXPO(3)
         ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CALCULATE VAB AND DV(1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          THE POTENTIAL CALCULATED HERE AND THE POTENTIALS VAC AND VBC
C          CORRESPOND TO THE VOH POTENTIAL DESCRIBED IN EQUATION 3.
C          DV(1) IS THE DERIVATIVE OF THIS TERM WITH RESPECT TO RAB.
C
      JUNK(1) = EXP(GAMMA(1)*(REQ(1)-RAB))
      VAB = DE(1)*(JUNK(1)*JUNK(1)-2.0D0*JUNK(1))
      IF (NDER .EQ. 1) 
     *    DV(1) = 2.0D0*DE(1)*(1.0D0-JUNK(1))*GAMMA(1)*JUNK(1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CALCULATE VAC AND DV(3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      JUNK(1) = EXP(GAMMA(3)*(REQ(3)-RAC))
      VAC = DE(3)*(JUNK(1)*JUNK(1)-2.0D0*JUNK(1))
      IF (NDER .EQ. 1)
     *    DV(3) = 2.0D0*DE(3)*(1.0D0-JUNK(1))*GAMMA(3)*JUNK(1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CALCULATE VBC AND DV(2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      JUNK(1) = EXP(GAMMA(2)*(REQ(2)-RBC))
      VBC = DE(2)*(JUNK(1)*JUNK(1)-2.0D0*JUNK(1))
      IF (NDER .EQ. 1) 
     *    DV(2) = -2.0D0*DE(2)*(JUNK(1)-1.0D0)*GAMMA(2)*JUNK(1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CALCULATE CSTH AND DERIVATIVES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          TH IS THE BEND ANGLE OF THE THREE ATOM SYSTEM.
C          IT IS NEVER USED USED EXPLICITELY WITHIN THE PROGRAM
C          HOWEVER COSINE(TH) IS.  COSINE(TH) IS CALCULATED
C          FROM THE LAW OF COSINES)
C
      NUM = RBC*RBC+RAB*RAB-RAC*RAC
      DENOM = 2.0D0*RAB*RBC
C
      CSTH = NUM/DENOM
C
      IF (NDER .EQ. 1) THEN
          DCSTH(1) = (2.0D0*RAB/DENOM)-(2.0D0*RBC)*CSTH/DENOM
          DCSTH(2) = (2.0D0*RBC/DENOM)-(2.0D0*RAB)*CSTH/DENOM
          DCSTH(3) = -(2.0D0*RAC/DENOM)
      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CALCULATE THE B IJ 'S--SEE EQUATION 5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO 80 I = 1, 5
         DO 70 J = 1, 3
            B(I,J) = C(I,J,1)*(1.0D0+C(I,J,2)*CSTH*(1.0D0+C(I,J,3)*
     *               CSTH))
C
C               CALCULATE THE DERIVATIVES TOO.
C
            IF (NDER .EQ. 1) THEN
                DO 60 IT = 1, 3
                      DB(I,J,IT) = C(I,J,1)*C(I,J,2)*(CSTH*C(I,J,3)*
     *                             DCSTH(IT)+DCSTH(IT)*
     *                             (1.0D0+C(I,J,3)*CSTH))
   60           CONTINUE
            ENDIF
   70    CONTINUE
   80 CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          NEXT CALCULATE MEW--SEE EQUATION FOUR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      MEW = EXPO(2)-0.02538745816335D0
      IF (NDER .EQ. 1) THEN
          DMEW(2) = -BETA*EXPO(2)
          DMEW(1) = 0.0D0
          DMEW(3) = 0.0D0
      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          NOW CALCULATE A I 'S AND DERIVATIVES--SEE EQUATION 4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO 100 IT = 1, 5
         JUNK(1) = 1.0D0+B(IT,3)*MEW
         JUNK(2) = MEW*B(IT,2)
         JUNK(3) = 1.0D0+JUNK(2)*JUNK(1)
C
         A(IT) = B(IT,1)*JUNK(3)
         IF (NDER .EQ. 1) THEN
             DO 90 IT1 = 1, 3
                   DA(IT,IT1) = DB(IT,1,IT1)*JUNK(3)+B(IT,1)*
     *                          ((DB(IT,2,IT1)*MEW+B(IT,2)*DMEW(IT1))*
     *                          JUNK(1)+JUNK(2)*(DB(IT,3,IT1)*MEW+   
     *                          B(IT,3)*DMEW(IT1)))
   90        CONTINUE
         ENDIF
  100 CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CALCULATE THE SECOND TERM IN EQUATION THREE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PROD1 = EXPO(1)*EXPO(3)*A(1)
C
      PROD2 = 1.0D0+(A(3)*RAB*RAC)/(RAB+RAC)
C
      IF (NDER .EQ. 1) THEN
          DPROD1(1) = DEXPO(1)*EXPO(3)*A(1)+EXPO(1)*EXPO(3)*DA(1,1)
          DPROD1(2) = EXPO(1)*EXPO(3)*DA(1,2)
          DPROD1(3) = DEXPO(3)*EXPO(1)*A(1)+EXPO(1)*EXPO(3)*DA(1,3)
C
          DPROD2(2) = (DA(3,2)*RAB*RAC)/(RAB+RAC)
          DPROD2(1) = ((DA(3,1)*RAB*RAC+A(3)*RAC)*(RAB+RAC)-
     *                 (A(3)*RAB*RAC))/((RAB+RAC)*(RAB+RAC))
          DPROD2(3) = ((DA(3,3)*RAB*RAC+A(3)*RAB)*(RAB+RAC)-
     *                 (A(3)*RAB*RAC))/((RAB+RAC)*(RAB+RAC))
      ENDIF
C
C          CALCULATE PROD3
C
      PROD3 = 1.0D0+A(2)*(RAB+RAC)*PROD2
C
      IF (NDER .EQ. 1) THEN
          DO 120 IT = 1, 3
                 DPROD3(IT) = DA(2,IT)*(RAB+RAC)*PROD2+A(2)*
     *                        (RAB+RAC)*DPROD2(IT)
                 IF (IT.EQ.2) GO TO 110
                 DPROD3(IT) = DPROD3(IT)+A(2)*PROD2
  110            CONTINUE
  120     CONTINUE
      ENDIF
C
      SECOND = PROD1*PROD3
C
      IF (NDER .EQ. 1) THEN
          DO 130 IT = 1, 3
                 DSEC(IT) = PROD1*DPROD3(IT)+DPROD1(IT)*PROD3
  130     CONTINUE
      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CALCULATE THE THIRD TERM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      JUNK(1) = EXPO(1)+EXPO(3)
      JUNK(2) = (1.0D0+A(5)*(RAB+RAC))
C
      THIRD = JUNK(1)*EXPO(2)*A(4)*JUNK(2)
C
      IF (NDER .EQ. 1) THEN
          DO 140 IT1 = 1, 2
                 IT = 1
                 IF (IT1.EQ.2) IT = 3
                 DTHIRD(IT) = DEXPO(IT)*EXPO(2)*A(4)*JUNK(2)+
     *                        JUNK(1)*EXPO(2)*DA(4,IT)*JUNK(2)+ 
     *                        JUNK(1)*EXPO(2)*A(4)*(DA(5,IT)*
     *                        (RAB+RAC)+A(5))
  140     CONTINUE
          DTHIRD(2) = JUNK(1)*(DEXPO(2)*A(4)*JUNK(2)+EXPO(2)*(DA(4,2)*
     *                JUNK(2)+A(4)*(DA(5,2)*(RAB+RAC))))
      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          FINALLY SUM IT ALL UP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ENERGY = VAB+VAC+VBC+SECOND+THIRD+DE(2)
      ENERGY = ENERGY+ZEROAD
C
      IF (NDER .EQ. 1) THEN
          DO 150 IT = 1, 3
                 DER(IT) = DV(IT)+DSEC(IT)+DTHIRD(IT)
  150     CONTINUE
C
          IF (R(1).LT.R(3)) GO TO 160
          JUNK(1) = DER(3)
          DER(3) = DER(1)
          DER(1) = JUNK(1)
  160     CONTINUE
      ENDIF
C
CCCCCCCCCCCCCCCC
C
      RETURN
C
CCCCCCCCCCCCCCCC
C
1000  FORMAT (/, 2X, T5, 'PREPOT has been called for H + OO',
     *        /, 2X, T5, 'Potential energy function by Melius ',
     *                   'and Blint',
     *        //, 2X, T5, 'Potential energy surface parameters:',
     *        /, 2X, T10, 'Bond', T47, 'H-O', T58, 'O-O', T69, 'O-H')
 1050 FORMAT(2X,T10,'Dissociation energies (hartrees):',
     *       T44, F10.5, T55, F10.5, T66, F10.5)
 1150 FORMAT(2X,T10,'Morse betas (reciprocal bohr):',
     *       T44, F10.5, T55, F10.5, T66, F10.5)
 1200 FORMAT(2X,T10,'Equilibrium bond lengths (bohr):',
     *       T44, F10.5, T55, F10.5, T66, F10.5)
 1300 FORMAT(/,2X,T10,'Constant added to the energy:',
     *       /,T38, 1PE20.12,' hartree atomic units',
     *       /,T38, 1PE20.12,' eV',/,T38,1PE20.12,' kcal/mol') 
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
C
      END
C
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IDUM(19)
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT /6/
         DATA NSURF /0/
         DATA NDER /1/
         DATA NDUM /21*0/
C
         END
C*****
