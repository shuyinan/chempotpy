      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
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

      NDER=igrad

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
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          FH2                                                        
C   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato)                
C                    plus E3C (a three-center term).                            
C   Common name:     no. 5A                                                     
C   Reference:       R. Steckler, D. G. Truhlar, and B. C. Garrett              
C                    J. Chem. Phys. 82, 5499-5505 (1985).                       
C                                                                               
C   Notes:           In this version of the potential the AB and AC             
C                    Sato parameters are functions of PHI, where PHI is         
C                    the angle of BC-to-A vector with the BC axis.              
C                    The BC Sato parameter is a function of BC bond             
C                    length and the angle PHI.  A three-center term is          
C                    also added to correct the H-F-H barrier height.            
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in DATA statements                    
C   and in the block data subprogram PTPACM.                                    
C   Coordinates, potential energy, and derivatives are passed                   
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the FH diatomic, with the             
C   FH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,           
C   respectively.                                                               
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(F-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-F)                                                
C   The zero of energy is defined at F "infinitely" far from the H2 diatomic.   
C                                                                               
C   The flags that indicate what calculations should be carried out in          
C   the potential routine are passed through the common block PT2CM:            
C   where:                                                                      
C        NASURF - which electronic states availalble                            
C                 (1,1) = 1 as only gs state available                          
C        NDER  = 0 => no derivatives should be calculated                       
C        NDER  = 1 => calculate first derivatives                               
C        NFLAG - these integer values can be used to flag options               
C                within the potential;                                          
C                                                                               
C                                                                               
C   Potential parameters' default settings                                      
C                  Variable            Default value                            
C                  NDER                1                                        
C                  NFLAG(18)           6                                        
C                                                                               
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                    
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
         PARAMETER (CKCAU = 627.5095D0)                                         
         PARAMETER (CANGAU = 0.52917706D0)                                      
C                                                                               
C   EXLARG IS LN(LARGEST FLOATING POINT NUMBER ALLOWED)                         
C                                                                               
         COMMON /SATOCM/ DM(3), RE(3), BETA(3)                                  
C                                                                               
C      LOGICAL LCOL, LDERIV, LZERO                                              
C                                                                               
         COMMON /PRECM/ Z(3),Q(3),XJ(3),DQ(3,3),DJ(3,3),                        
     +                  DZ(3,3),RR(3),A1,A2,A,B,C,D,F,G,                        
     +                  BET,BETP,ZFH0,ZHH0,ZHHST,DEG1,                          
     +                  DEG2,DEG3,DEG4,DEG3P,C1,C2,C3,C4,                       
     +                  C5,C6,C3C,ALF,ALFP,ALFT,P3C,Q3C,                        
     +                  ZHHDIF,OMP3C,OMQ3C,HA2,G2                               
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the FH2 ',                    
     *                'extended-LEPS plus three-center term ',                  
     *        /,2X,T5,'potential energy surface no. 5A',                        
     *       //,2X,T5,'Potential energy surface parameters:',                   
     *        /,2X,T5,'Bond', T47, 'F-H', T58, 'H-H', T69, 'H-F',               
     *        /,2X,T5,'Dissociation energies (kcal/mol):',                      
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):',                  
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):',                 
     *        T44, F10.5, T55, F10.5, T66, F10.5)                               
610   FORMAT (/,2X,T5,'HH Sato parameter fit',                                  
     *             T40,'FH Sato parameter fit',                                 
     * /,2X,T5,'A',   T10,'=',T12,1PE13.5, T40,'ZFH0',T46,'=',T48,E13.5,        
     * /,2X,T5,'A1',  T10,'=',T12,E13.5, T40,'DEG1', T46,'=',T48,1PE13.5,       
     * /,2X,T5,'A2',  T10,'=',T12,E13.5, T40,'DEG2', T46,'=',T48,E13.5,         
     * /,2X,T5,'B',   T10,'=',T12,E13.5, T40,'DEG3', T46,'=',T48,E13.5,         
     * /,2X,T5,'C',   T10,'=',T12,E13.5, T40,'DEG4', T46,'=',T48,E13.5,         
     * /,2X,T5,'D',   T10,'=',T12,E13.5, T40,'DEG3P',T46,'=',T48,E13.5,         
     * /,2X,T5,'F',   T10,'=',T12,E13.5, T40,'C1',   T46,'=',T48,E13.5,         
     * /,2X,T5,'G',   T10,'=',T12,E13.5, T40,'C2',   T46,'=',T48,E13.5,         
     * /,2X,T5,'BET', T10,'=',T12,E13.5, T40,'C3',   T46,'=',T48,E13.5,         
     * /,2X,T5,'BETP',T10,'=',T12,E13.5, T40,'C4',   T46,'=',T48,E13.5,         
     * /,2X,T40,'C5',T46,'=',T48,E13.5,                                         
     * /,2X,T40,'C6',T46,'=',T48,E13.5)                                         
620   FORMAT(/, 2X, T5, 'Parameters for the three-center term:',                
     * /,2X,T5,'C3C', T10,'=',T12,E13.5,T40,'ALF', T46,'=',T48,E13.5,           
     * /,2X,T5,'ALFP',T10,'=',T12,E13.5,T40,'ALFT',T46,'=',T48,E13.5,           
     * /,2X,T5,'P3C', T10,'=',T12,E13.5,T40,'Q3C', T46,'=',T48,E13.5)           
C                                                                               
      DO 10 I = 1,3                                                             
         DM(I) = DM(I) / CKCAU                                                  
         RE(I) = RE(I) / CANGAU                                                 
         BETA(I) = BETA(I) * CANGAU                                             
   10 CONTINUE                                                                  
C                                                                               
C   Set the energy of the potential in the each of the three asymptotic regions 
C                                                                               
         EASYAB = DM(1)                                                         
         EASYBC = DM(2)                                                         
         EASYAC = DM(3)                                                         
C                                                                               
C   USEFUL CONSTANTS                                                            
C                                                                               
      ZHHDIF = ZHHST - ZHH0                                                     
      OMP3C = 1.D0 - P3C                                                        
      OMQ3C = 1.D0 - Q3C                                                        
      HA2 = 0.5D0 * A2                                                          
      G2 = G * G                                                                
C                                                                               
C                                                                               
      EZERO(1)=DM(2)                                                            
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='R. Steckler, D. G. Truhlar, B. C. Garrett,'                      
       REF(2)='J. Chem. Phys. 82, 5499(1985)'                                   
C                                                                               
      INDEXES(1) = 9                                                            
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
C                                                                               
C                                                                               
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
c      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the FH diatomic, with the             
C   FH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(F-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-F)                                                
C   The zero of energy is defined at F "infinitely" far from the H2 diatomic.   
C                                                                               
C      ENTRY POT                                                                
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                    
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                            
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                         
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                         
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
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
         PARAMETER (CKCAU = 627.5095D0)                                         
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (TDEG = 57.29577951D0)                                       
         PARAMETER (EXLARG = 80.D0)                                             
C                                                                               
         COMMON /SATOCM/ DM(3), RE(3), BETA(3)                                  
C                                                                               
      LOGICAL LCOL, LDERIV, LZERO                                               
C                                                                               
         DATA LDERIV /.FALSE./                                                  
C                                                                               
         COMMON /PRECM/ Z(3),Q(3),XJ(3),DQ(3,3),DJ(3,3),                        
     +                  DZ(3,3),RR(3),A1,A2,A,B,C,D,F,G,                        
     +                  BET,BETP,ZFH0,ZHH0,ZHHST,DEG1,                          
     +                  DEG2,DEG3,DEG4,DEG3P,C1,C2,C3,C4,                       
     +                  C5,C6,C3C,ALF,ALFP,ALFT,P3C,Q3C,                        
     +                  ZHHDIF,OMP3C,OMQ3C,HA2,G2                               
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C   Check the values of NASURF and NDER for validity.                           
C                                                                               
C   Check the value of NDER                                                     
C                                                                               
C                                                                               
C   The logical LDERIV is used to flag the derivative calculations              
C   in the code; LDERIV is set equal to .TRUE. if NDER = 1.                     
C                                                                               
         IF (NDER .EQ. 1) LDERIV = .TRUE.                                       
C                                                                               
C   TEST R'S FOR PHYSICAL VALUES                                                
C                                                                               
      DO 20 I = 1,3                                                             
         RR(I) = R(I)                                                           
   20 CONTINUE                                                                  
      DO 30 I1 = 1,3                                                            
         I2 = MOD(I1,3) + 1                                                     
         I3 = MOD(I2,3) + 1                                                     
         T = RR(I2) + RR(I3)                                                    
         IF (RR(I1) .GT. T) RR(I1) = T                                          
   30 CONTINUE                                                                  
      R12 = RR(1) * RR(1)                                                       
      R22 = RR(2) * RR(2)                                                       
      R32 = RR(3) * RR(3)                                                       
C                                                                               
C   CAPR IS THE F TO H2 DISTANCE                                                
C                                                                               
      CAPR2 = 0.5D0 * (R12 + R32 - 0.5D0 * R22)                                 
      IF (CAPR2 .LT. 0.D0) CAPR2 = 0.0D0                                        
      CAPR = SQRT(CAPR2)                                                        
C                                                                               
C   PHI IS THE ANGLE BETWEEN CAPR AND RHH                                       
C   PHI IS RESTRICTED TO BE BETWEEN O AND 90 DEGREES                            
C                                                                               
      T = 0.5D0 * (R32 - R12)                                                   
      IF (ABS(T) .GT. 1.D-7) T = T / (RR(2) * CAPR)                             
      T2 = 1.0D0                                                                
      SGNPHI = SIGN(T2,T)                                                       
      COSPHI = SGNPHI * T                                                       
      COSPHI = MIN(COSPHI,T2)                                                   
      PHI = ACOS(COSPHI)                                                        
      SINP2 = 1.D0 - COSPHI * COSPHI                                            
      IF(SINP2 .LT. 0.D0) SINP2 = 0.D0                                          
      SINPHI = SQRT(SINP2)                                                      
C                                                                               
C   LCOL IS SET TRUE IF COLLINEAR GEOMETRY                                      
C                                                                               
      LCOL = ABS(SINPHI) .LT. 1.D-6 .OR. ABS(CAPR) .LT. 1.D-6                   
C                                                                               
C   COMPUTE COLLINEAR SATO PARAMETERS                                           
C                                                                               
      Z(1) = ZFH0                                                               
      T = A * (RR(2) - B)                                                       
      IF (.NOT.(T .LT. EXLARG)) GO TO 35                                        
         T = EXP(T)                                                             
         T1 = 1.D0 / T                                                          
         HSECH = 1.D0/(T + T1)                                                  
         ZHH = A1 + HA2 * ((T - T1) * HSECH + 1.D0)                             
         GO TO 36                                                               
   35 CONTINUE                                                                  
         HSECH = 0.D0                                                           
         ZHH = A1 + A2                                                          
   36 CONTINUE                                                                  
      T = RR(2) - F                                                             
      EXX = C * EXP(D * T * T)                                                  
      ZHH = ZHH - EXX                                                           
      Z(2) = ZHH                                                                
      IF (.NOT.LDERIV) GO TO 50                                                 
         DO 40 I = 1,3                                                          
            DO 40 J = 1,3                                                       
               DZ(I,J) = 0.D0                                                   
  40  CONTINUE                                                                  
C                                                                               
C  ONLY THE HH SATO PARAMETER HAS A NONZERO DERIVATIVE - W.R.T. R2              
C                                                                               
         DZ(2,2) = 2.D0 * A * A2 * HSECH * HSECH                                
     *      - 2.D0 * D * T * EXX                                                
   50 CONTINUE                                                                  
C                                                                               
C   COMPUTE NONCOLLINEAR CONTRIBUTIONS TO SATO PARAMETERS                       
C                                                                               
      IF (LCOL) GO TO 160                                                       
C                                                                               
C   COMPUTE Z(PHI) FOR FH SATO PARAMETER                                        
C   DZPHI IS THE DERIVATIVE OF Z(PHI) W.R.T. PHI                                
C                                                                               
         DEG = TDEG * PHI                                                       
         IF (.NOT.(DEG.LT.0.D0 .OR. DEG.GT.90.D0)) GO TO 60                     
            WRITE(NFLAG(18),6000) DEG                                           
            STOP 'POT 3'                                                        
   60    CONTINUE                                                               
         IF (.NOT.(DEG .LE. DEG1)) GO TO 70                                     
            ZPHI = 0.D0                                                         
            DZPHI = 0.D0                                                        
            GO TO 140                                                           
   70    CONTINUE                                                               
            IF (.NOT.(DEG .LE. DEG2)) GO TO 80                                  
               DEGDIF = DEG - DEG1                                              
               DD2 = DEGDIF * DEGDIF                                            
               ZPHI = DD2 * (C1 + C2 * DEGDIF)                                  
               IF (LDERIV) DZPHI = TDEG *                                       
     *            DEGDIF * (2.D0 * C1 + 3.D0 * C2 * DEGDIF)                     
               GO TO 130                                                        
   80       CONTINUE                                                            
               IF (.NOT.(DEG .LE. DEG3)) GO TO 90                               
                  DEGDIF = DEG - DEG3P                                          
                  ZPHI = C3 + C4 * DEGDIF                                       
                  IF (LDERIV) DZPHI = C4 * TDEG                                 
                  GO TO 120                                                     
   90          CONTINUE                                                         
                  IF (.NOT.(DEG .LE. DEG4)) GO TO 100                           
                     ZPHI = C5 + C6 * SINPHI * SINPHI                           
                     IF (LDERIV) DZPHI = 2.D0 * C6 * SINPHI * COSPHI            
                     GO TO 110                                                  
  100             CONTINUE                                                      
                     WRITE (NFLAG(18),6001) DEG                                 
                     STOP 'POT 4'                                               
  110             CONTINUE                                                      
  120          CONTINUE                                                         
  130       CONTINUE                                                            
  140  CONTINUE                                                                 
C                                                                               
C   COMPUTE SCALE FACTOR TO DAMP OUT CORRECTION AT LARGE DISTANCE               
C                                                                               
         SCALE1 = CAPR2 / (G2 + CAPR2)                                          
         Z(1) = Z(1) + SCALE1 * ZPHI                                            
C                                                                               
C   COMPUTE SCALE FACTOR FOR HH SATO PARAMETER                                  
C                                                                               
         SCALE2 = (ZHH - ZHH0) / ZHHDIF                                         
C                                                                               
C   COMPUTE PHI DEPENDENT PART OF HH SATO PARAMETER                             
C                                                                               
         CPHI = (BET + BETP * SINP2) * SINP2                                    
C                                                                               
C   COMPUTE HH SATO PARAMTER                                                    
C                                                                               
         Z(2) = Z(2) + SCALE1 * SCALE2 * CPHI                                   
         IF (.NOT.LDERIV) GO TO 150                                             
C                                                                               
C   COMPUTE DERIVATIVES                                                         
C   FIRST, DERIVATIVES OF PHI W.R.T. R1, R2, R3                                 
C                                                                               
            T = 1.D0 / CAPR                                                     
            T1 = 0.5D0 * COSPHI * T                                             
            T2 = SGNPHI / RR(2)                                                 
            DPHI1 = RR(1) * (T1 + T2)                                           
            DPHI3 = RR(3) * (T1 - T2)                                           
            T1 = 1.D0 / SINPHI                                                  
            T2 = T1 * T                                                         
            DPHI1 = DPHI1 * T2                                                  
            DPHI3 = DPHI3 * T2                                                  
            DPHI2 = RR(2) * COSPHI * T1 * (1.D0 / R22 - 0.25D0 / CAPR2)         
            T1 = (1.D0 - SCALE1) / CAPR2                                        
            T2 = T1 * ZPHI                                                      
C                                                                               
C   NEXT, DERIVATIVES OF FH SATO PARAMETERS                                     
C                                                                               
            DZ(1,1) = SCALE1 * (T2 * RR(1) + DZPHI * DPHI1)                     
            DZ(1,2) = SCALE1 * (-0.5D0 * T2 * RR(2) + DZPHI * DPHI2)            
            DZ(1,3) = SCALE1 * (T2 * RR(3) + DZPHI * DPHI3)                     
            T1 = T1 * CPHI                                                      
C                                                                               
C   DERIVATIVE OF THE PHI DEPENDENT PART OF HH SATO W.R.T. PHI                  
C                                                                               
            DCPHI = 2.D0 * (BET + 2.D0 * BETP * SINP2) * SINPHI * COSPHI        
            T2 = SCALE1 * SCALE2                                                
C                                                                               
C   LAST, DERIVATIVES OF HH SATO PARAMETER                                      
C                                                                               
            DZ(2,1) = T2 * (T1 * RR(1) + DCPHI * DPHI1)                         
            DZ(2,2) = DZ(2,2) * (1.D0 + SCALE1 * CPHI / ZHHDIF)                 
     *         + T2 * (-0.5D0 * T1 * RR(2) + DCPHI * DPHI2)                     
            DZ(2,3) = T2 * (T1 * RR(3) + DCPHI * DPHI3)                         
  150    CONTINUE                                                               
  160 CONTINUE                                                                  
C                                                                               
C   COMPUTE POTENTIAL FROM LEPS FORM                                            
C                                                                               
      ENGYGS = 0.D0                                                             
      Z(3) = Z(1)                                                               
      IF (.NOT.LDERIV) GO TO 170                                                
         DO 165 I = 1,3                                                         
            DZ(3,I) = DZ(1,I)                                                   
           DEGSDR(I) = 0.D0                                                     
            Q(I) = 0.D0                                                         
            XJ(I) = 0.D0                                                        
            DO 162 J = 1,3                                                      
               DQ(I,J) = 0.D0                                                   
               DJ(I,J) = 0.D0                                                   
  162       CONTINUE                                                            
  165    CONTINUE                                                               
  170 CONTINUE                                                                  
      LZERO = .TRUE.                                                            
      DO 200 I = 1,3                                                            
         EX = EXP(BETA(I) * (RE(I) - RR(I)))                                    
         LZERO = LZERO .AND. EX.EQ.0.0D0                                        
         IF (EX.EQ.0.0D0) GO TO 195                                             
            ZP1 = Z(I) + 1.D0                                                   
            ZP3 = ZP1 + 2.D0                                                    
            OP3Z = 1.D0 + 3.D0 * Z(I)                                           
            T = 0.25D0 * DM(I) * EX / ZP1                                       
            T1 = ZP3 * EX                                                       
            T2 = OP3Z * EX                                                      
C                                                                               
C   Q AND XJ ARE THE COULUMB AND EXCHANGE PARTS, INDEX I RUNS OVER THE          
C   COORDINATES OF WHICH THEY ARE EXPLICIT FUNCTIONS                            
C                                                                               
            Q(I) = T * (T1 - 2.D0  * OP3Z)                                      
            XJ(I) = T * (T2 - 2.D0  * ZP3)                                      
C                                                                               
C   PUT SUM OF Q'S IN ENERGY                                                    
C                                                                               
            ENGYGS = ENGYGS + Q(I)                                              
            IF (.NOT.LDERIV) GO TO 190                                          
C                                                                               
C   DERIVATIVE OF Q AND ZJ W.R.T. R1, R2, R3                                    
C   ELEMENT I,J IS THE DERIVATIVE OF Q(I) W.R.T. RJ                             
C                                                                               
               T3 = 2.D0 * T * (EX + 2.D0) / ZP1                                
               DO 180 J = 1,3                                                   
                  DQ(I,J) = -T3 * DZ(I,J)                                       
                  DJ(I,J) = -DQ(I,J)                                            
  180          CONTINUE                                                         
               T = 2.D0*T*BETA(I)                                               
               DQ(I,I) = DQ(I,I) - T * (T1 - OP3Z)                              
               DJ(I,I) = DJ(I,I) - T * (T2 - ZP3)                               
  190       CONTINUE                                                            
  195    CONTINUE                                                               
  200 CONTINUE                                                                  
      IF (LZERO) GO TO 270                                                      
         XJX = 0.D0                                                             
C                                                                               
C   COMPUTE THE TOTAL EXCHANGE TERM                                             
C                                                                               
         DO 230 I1 = 1,3                                                        
            I2 = MOD(I1,3) + 1                                                  
            T = XJ(I1) - XJ(I2)                                                 
            XJX = XJX + T * T                                                   
            IF (.NOT.LDERIV) GO TO 220                                          
C                                                                               
C   PUT THE DERIVATIVES OF THE TOTAL EXCHANGE TERM IN DEGSDR                    
C                                                                               
               DO 210 J = 1,3                                                   
                 DEGSDR(J) =DEGSDR(J) + T * (DJ(I1,J) - DJ(I2,J))               
  210          CONTINUE                                                         
  220       CONTINUE                                                            
  230    CONTINUE                                                               
         XJX = SQRT(0.5D0 * XJX)                                                
C                                                                               
C   COMPUTE THE LEPS PART OF THE POTENTIAL ENERGY                               
C                                                                               
         ENGYGS = ENGYGS - XJX                                                  
         IF (.NOT.LDERIV) GO TO 260                                             
C                                                                               
C   COMPUTE DERIVATIVES OF THE LEPS PART                                        
C                                                                               
            IF(ABS(XJX) .GT. 1.D-14) T = 0.5D0 / XJX                            
            DO 250 J = 1,3                                                      
C                                                                               
C   THE DERIVATIVE OF THE EXCHANGE PART MUST BE DIVIDED BY 2*XJX                
C                                                                               
              DEGSDR(J) = -T *DEGSDR(J)                                         
               DO 240  I = 1,3                                                  
C                                                                               
C   THEN ADD IN THE DERIVATIVE OF Q                                             
C                                                                               
                 DEGSDR(J) =DEGSDR(J) + DQ(I,J)                                 
  240          CONTINUE                                                         
  250       CONTINUE                                                            
  260    CONTINUE                                                               
  270 CONTINUE                                                                  
C                                                                               
C   EZERO ADDED TO PUT ZERO AT REQ FOR REACTANTS                                
C                                                                               
      ENGYGS = ENGYGS + EZERO(1)                                                
C                                                                               
C   COMPUTE THE 3-CENTER TERM                                                   
C                                                                               
      RSUM = RR(1) + RR(3)                                                      
      RDIF = RR(1) - RR(3)                                                      
      RDIF2 = RDIF * RDIF                                                       
      T1 = -ALFP * RSUM                                                         
      EX1 = C3C * EXP(T1 * RSUM)                                                
      T2 = -ALF * RDIF                                                          
      EX2 = OMP3C * EXP(T2 * RDIF)                                              
      T3 = -ALFT * RDIF2                                                        
      EX3 = P3C * EXP(T3 * RDIF2)                                               
      R1I = 1.D0 / RR(1)                                                        
      R3I = 1.D0 / RR(3)                                                        
      T4 = R1I*R3I                                                              
C                                                                               
C   COSTH IS THE ANGLE BETWEEN R1 AND R3                                        
C                                                                               
      COSTH = 0.5D0 * (R12 + R32 - R22) * T4                                    
      T5 = OMQ3C * COSTH                                                        
      T = Q3C + T5 * COSTH                                                      
C                                                                               
C   E3C CONTAINS THE THREE CENTER CORRECTION TO THE ENERGY                      
C                                                                               
      E3C = (EX2 + EX3) * T * EX1                                               
      ENGYGS = ENGYGS + E3C                                                     
      IF (.NOT.LDERIV) GO TO 280                                                
C                                                                               
C   COMPUTE DERIVATIVE OF THE 3C TERM                                           
C   FIRST, DERIVATIVE OF COSTH                                                  
C                                                                               
         DCOS1 = R3I - COSTH*R1I                                                
         DCOS2 = - RR(2)*T4                                                     
         DCOS3 = R1I - COSTH*R3I                                                
         T4 = EX2 + EX3                                                         
         T5 = T4 * T5                                                           
C                                                                               
C   NOW, DERIVATIVES OF E3C                                                     
C                                                                               
         DE3C1 = 2.D0*((T2 * EX2 + 2.D0 * T3 * RDIF * EX3 + T1 * T4)            
     *      * T + T5 * DCOS1) * EX1                                             
         DE3C2 = 2.D0 * T5 * DCOS2 * EX1                                        
         DE3C3 = 2.D0*((-T2 * EX2 - 2.D0 * T3 * RDIF * EX3 + T1 * T4)           
     *      * T + T5 * DCOS3) * EX1                                             
        DEGSDR(1) = DEGSDR(1) + DE3C1                                           
        DEGSDR(2) = DEGSDR(2) + DE3C2                                           
        DEGSDR(3) = DEGSDR(3) + DE3C3                                           
  280 CONTINUE                                                                  
600   FORMAT (/,2X,T5,'PREPOT has been called for the FH2 ',                    
     *                'extended-LEPS plus three-center term ',                  
     *        /,2X,T5,'potential energy surface no. 5A',                        
     *       //,2X,T5,'Potential energy surface parameters:',                   
     *        /,2X,T5,'Bond', T47, 'F-H', T58, 'H-H', T69, 'H-F',               
     *        /,2X,T5,'Dissociation energies (kcal/mol):',                      
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):',                  
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):',                 
     *        T44, F10.5, T55, F10.5, T66, F10.5)                               
610   FORMAT (/,2X,T5,'HH Sato parameter fit',                                  
     *             T40,'FH Sato parameter fit',                                 
     * /,2X,T5,'A',   T10,'=',T12,1PE13.5, T40,'ZFH0',T46,'=',T48,E13.5,        
     * /,2X,T5,'A1',  T10,'=',T12,E13.5, T40,'DEG1', T46,'=',T48,1PE13.5,       
     * /,2X,T5,'A2',  T10,'=',T12,E13.5, T40,'DEG2', T46,'=',T48,E13.5,         
     * /,2X,T5,'B',   T10,'=',T12,E13.5, T40,'DEG3', T46,'=',T48,E13.5,         
     * /,2X,T5,'C',   T10,'=',T12,E13.5, T40,'DEG4', T46,'=',T48,E13.5,         
     * /,2X,T5,'D',   T10,'=',T12,E13.5, T40,'DEG3P',T46,'=',T48,E13.5,         
     * /,2X,T5,'F',   T10,'=',T12,E13.5, T40,'C1',   T46,'=',T48,E13.5,         
     * /,2X,T5,'G',   T10,'=',T12,E13.5, T40,'C2',   T46,'=',T48,E13.5,         
     * /,2X,T5,'BET', T10,'=',T12,E13.5, T40,'C3',   T46,'=',T48,E13.5,         
     * /,2X,T5,'BETP',T10,'=',T12,E13.5, T40,'C4',   T46,'=',T48,E13.5,         
     * /,2X,T40,'C5',T46,'=',T48,E13.5,                                         
     * /,2X,T40,'C6',T46,'=',T48,E13.5)                                         
620   FORMAT(/, 2X, T5, 'Parameters for the three-center term:',                
     * /,2X,T5,'C3C', T10,'=',T12,E13.5,T40,'ALF', T46,'=',T48,E13.5,           
     * /,2X,T5,'ALFP',T10,'=',T12,E13.5,T40,'ALFT',T46,'=',T48,E13.5,           
     * /,2X,T5,'P3C', T10,'=',T12,E13.5,T40,'Q3C', T46,'=',T48,E13.5)           
6000  FORMAT(/,2X,T5,'In the FH2 surface no. 5A, PHI = ', F10.5,                
     *               'degrees',                                                 
     *       /,2X,T5,'but PHI must be greater than zero.')                      
6001  FORMAT(/,2X,T5,'In the FH2 surface no. 5A, PHI = ', F10.5,                
     *               'degrees',                                                 
     *       /,2X,T5,'but PHI must be less than 90 degrees.')                   
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,                       
     *       /, 2X, 'This value of NDER is not allowed in this ',               
     *              'version of the potential.')                                
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C*****                                                                          
C                                                                               
         BLOCK DATA PTPACM                                                      
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                    
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
         COMMON /SATOCM/ DM(3), RE(3), BETA(3)                                  
C                                                                               
         COMMON /PRECM/ Z(3),Q(3),XJ(3),DQ(3,3),DJ(3,3),                        
     +                  DZ(3,3),RR(3),A1,A2,A,B,C,D,F,G,                        
     +                  BET,BETP,ZFH0,ZHH0,ZHHST,DEG1,                          
     +                  DEG2,DEG3,DEG4,DEG3P,C1,C2,C3,C4,                       
     +                  C5,C6,C3C,ALF,ALFP,ALFT,P3C,Q3C,                        
     +                  ZHHDIF,OMP3C,OMQ3C,HA2,G2                               
C                                                                               
C   Initialize the flags and the I/O unit numbers for the potential             
C                                                                               
      DATA NASURF /1,35*0/                                                      
      DATA NDER /0/                                                             
         DATA NFLAG /1,1,15*0,6,0,0/                                            
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/3,0,1/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /3/                                                           
C                                                                               
C   Initialize the potential parameters; the energy                             
C   parameters are in kcal/mol, and the lengths are in Angstroms.               
C                                                                               
         DATA DM /141.196D0, 109.449D0, 141.196D0/                              
         DATA RE /0.9170D0, 0.7419D0, 0.9170D0/                                 
         DATA BETA /2.2187D0, 1.9420D0, 2.2187D0/                               
C                                                                               
      DATA A1,A2,A,B,C,D,F,G                                                    
     +     /0.0395D0,0.201D0,1.26D0,1.60D0,0.0D0,0.0D0,                         
     +      0.0D0,0.2D0/                                                        
      DATA BET,BETP,ZFH0,ZHH0,ZHHST                                             
     +     /-0.101168D0,-0.46183D0,0.170D0,                                     
     +       0.120D0,0.23107D0/                                                 
      DATA DEG1,DEG2,DEG3,DEG4,DEG3P                                            
     +     /10.0D0,20.0D0,60.0D0,90.0D0,41.0D0/                                 
      DATA C1,C2,C3,C4,C5,C6                                                    
     +     /1.15D-4,-6.0D-6,0.016D0,0.0005D0,0.00069D0,                         
     +    3.308D-2/                                                             
      DATA C3C,ALF,ALFP,ALFT,P3C,Q3C                                            
     +     /0.270063D0,0.18D0,0.078D0,2.14D0,                                   
     +      0.95D0,0.48D0/                                                      
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
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
