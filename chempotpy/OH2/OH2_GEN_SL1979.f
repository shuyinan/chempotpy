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

      NDER=igrad
      NASURF=1
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
          g(1,iatom,idir)=DGSCART(iatom,idir)*27.211386
        enddo
        enddo
      else if (igrad==2) then
        write (*,*) 'Only energy and gradient are available'
      endif

      endsubroutine



      Subroutine prepot                                                         
C                                                                               
C   System:          OH2                                                        
C   Functional form:                                                            
C   Common name:     OH2SL                                                      
C   Reference:                                                                  
C                                                                               
C   Interface:       potlib2001
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
C
C   Notes: This is a collinear potential energy surface, i.e., the potential    
C          energy is not defined for non-linear geometries.                     
C          The derivatives of the potential energy with respect to the          
C          coordinates are determined numerically using a step size STEP        
C          if R(1) and R(2) are less than 50 bohr.                              
C          The variable STEP is initialized in the block data subprogram        
C          PTPACM.                                                              
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are assigned in the Block Data Subprogram PTPACM.  
C   Coordinates, potential energy, and derivatives are passed                   
C   through the common block PT1CM:                                             
C                  COMMON /PT1CM/ R(3), ENGYGS, DEGSDR(3)                       
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(O-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-O)                                                
C                                                                               
C   The zero of energy is defined at O "infinitely" far from the H2 diatomic.   
C                                                                               
C   The flags that indicate which calculations should be carried out in         
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
         IMPLICIT REAL*8 (A-H,O-Z)                                              
C        IMPLICIT REAL*16(L)                                                    
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
         PARAMETER (BOHR = 0.529177D0)                                          
         PARAMETER (HART = 27.21161D0)                                          
C                                                                               
         COMMON /SATOCM/ DE(3), XRE(3), BET(3)                                  
         COMMON /NDERCM/ STEP                                                   
         COMMON /PRECM/  C(50),A,B,ALF,GAM,AL(3),                               
     +                   X(3),S(3),VM(3),CHI(3),V(4)                            
C                                                                               
C      DIMENSION C(50),AL(3),X(3),S(3),VM(3),R(3),CHI(3)                        
C                                                                               
C   Echo potential data to the file linked to UNIT NFLAG(18)                    
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
C      WRITE (NFLAG(18),600)                                                     
C      WRITE (NFLAG(18),601) STEP                                                
C      WRITE (NFLAG(18),602)                                                     
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the OH2 ',                    
     *                'Schinke-Lester potential energy surface')                
601   FORMAT(/,2X,T5,'Step size for the numerical derivatives = ',E10.5)        
602   FORMAT(2X,T5,'The Schinke-Lester surface has analytical ',                
     *             'derivatives if R(1) and R(2) are ',                         
     *       /,2X,T5,'less than 50 bohr, otherwise the derivatives ',           
     *               'are computed numerically.')                               
C                                                                               
      J=1                                                                       
C                                                                               
C    Set the values of the classical energy in the three asymptotic valleys     
C                                                                               
             EASYAB = DE(1)/HART                                                
             EASYBC = DE(2)/HART                                                
             EASYAC = DE(3)/HART                                                
C                                                                               
C                                                                               
      EZERO(1)=DE(2)                                                            
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='No Reference'                                                    
                                                                                
C                                                                               
      INDEXES(1) = 8                                                            
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
C                                                                               
C                                                                               
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
C      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C     CALCULATE ENERGY AND DERIVATIVE                                           
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(O-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-O)                                                
C                                                                               
C   The zero of energy is defined at O "infinitely" far from the H2 diatomic.   
C                                                                               
C      ENTRY POT                                                                
C                                                                               
         IMPLICIT REAL*8 (A-H,O-Z)                                              
C                                                                               
C        IMPLICIT REAL*16(L)                                                    
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
         PARAMETER (BOHR = 0.529177D0)                                          
         PARAMETER (HART = 27.21161D0)                                          
C                                                                               
         COMMON /SATOCM/ DE(3), XRE(3), BET(3)                                  
         COMMON /NDERCM/ STEP                                                   
         COMMON /PRECM/  C(50),A,B,ALF,GAM,AL(3),                               
     +                   X(3),S(3),VM(3),CHI(3),V(4)                            
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C   Check the values of NASURF and NDER for validity.                           
C                                                                               
C      IF (NASURF(1,1) .EQ. 0) THEN                                              
C         WRITE(NFLAG(18), 900) NASURF(1,1)                                      
C         STOP                                                                   
C      ENDIF                                                                     
C         IF (NDER .GT. 1) THEN                                                  
C             WRITE (NFLAG(18), 910) NDER                                        
C             STOP 'POT 2'                                                       
C         ENDIF                                                                  
C                                                                               
      IC=4                                                                      
      IF(R(1).GT.30.D0)IC=1                                                     
      IF(R(2).GT.30.D0)IC=1                                                     
C                                                                               
   86 GO TO (84,82,83,81),IC                                                    
C                                                                               
 81   SOH=R(1)                                                                  
      SR2=R(2)                                                                  
      GO TO 85                                                                  
C                                                                               
 82   SR2=R(2)+STEP                                                             
      SOH=R(1)                                                                  
      GO TO 85                                                                  
C                                                                               
 83   SR2=R(2)                                                                  
      SOH=R(1)+STEP                                                             
      GO TO 85                                                                  
C                                                                               
   84 SR2=R(2)                                                                  
      SOH=R(1)                                                                  
      R(3)=R(3) + STEP                                                          
C                                                                               
   85 SR2=SR2*0.529177D0                                                        
      SOH=SOH*0.529177D0                                                        
      R(3)= R(3)*0.529177D0                                                     
      R(1) = SOH                                                                
      R(2) = SR2                                                                
C                                                                               
      Q = 1.0D0                                                                 
C                                                                               
      DO 10 I=1,3                                                               
         X(I)=R(I)                                                              
         CHI(I)=(1.0D0 - EXP(-ALF*X(I)**2))                                     
         S(I)=X(I)-1.D0                                                         
         VM(I)=DE(I)*(1.0D0-EXP(-BET(I)*(X(I)-XRE(I))))**2-DE(I)                
         Q=Q*(1.0D0-TANH(AL(I)*S(I)))                                           
  10  continue                                                                  
C                                                                               
      U = VM(1)*CHI(2)*CHI(3)+VM(2)*CHI(1)*CHI(3)                               
     *	+ VM(3)*CHI(1)*CHI(2)                                                    
C                                                                               
      BP = B*(1.0D0/X(1)**4+1.0D0/X(3)**4)                                      
C                                                                               
      CHECK = (X(1)*X(3)*(X(2)-0.2D0))                                          
C                                                                               
      IF(CHECK.GT.1000.D0) AP = 0.0D0                                           
      IF(CHECK.GT.1000.D0) GO TO 112                                            
C                                                                               
      AP = A/(X(1)*X(3)*(X(2)-0.2D0))**12                                       
C                                                                               
  112 U = U + AP + BP                                                           
C                                                                               
      T = C(1) + C(2)*(S(1)+S(3)) + C(3)*S(2) +                                 
     *    C(4)*(S(1)**2 + S(3)**2) +                                            
     *	  C(5)*S(2)**2 + C(6)*S(2)*(S(1) + S(3)) +                               
     *    C(7)*S(1)*S(3) + C(8)*(S(1)**3 + S(3)**3) +                           
     *    C(9)*S(2)**3 + C(10)*S(1)*S(3)*(S(1) + S(3))                          
C                                                                               
      TT = C(11)*S(2)**2*(S(1) + S(3)) +                                        
     *     C(12)*S(2)*(S(1)**2 + S(3)**2) +                                     
     *     C(13)*S(1)*S(2)*S(3) + C(14)*(S(1)**4 + S(3)**4) +                   
     *     C(15)*S(2)**4 + C(16)*(S(1)**3*S(3) + S(1)*S(3)**3)                  
C                                                                               
      WW = C(17)*S(1)**2*S(3)**2 + C(18)*S(2)*(S(1)**3 + S(3)**3) +             
     *     C(19)*S(2)**2*(S(1)**2 + S(3)**2) +                                  
     *     C(20)*S(2)**3*(S(1) + S(3)) + C(21)*S(2)**2*S(1)*S(3) +              
     *     C(22)*S(2)*(S(1)**2*S(3) + S(1)*S(3)**2)                             
C                                                                               
      W = C(23)*(S(1)**5 + S(3)**5) +                                           
     *    C(24)*(S(1)**4*S(3) + S(1)*S(3)**4) +                                 
     *    C(25)*(S(1)**3*S(3)**2 + S(1)**2*S(3)**3) +                           
     *    C(26)*S(2)*(S(1)**4 + S(3)**4) +                                      
     *    C(27)*S(2)*(S(1)**3*S(3) + S(1)*S(3)**3)                              
C                                                                               
      GG = C(28)*S(2)*S(1)**2*S(3)**2 +                                         
     *     C(29)*S(2)**2*(S(1)**3 + S(3)**3) +                                  
     *     C(30)*S(2)**2*(S(1)**2*S(3) + S(1)*S(3)**2) +                        
     *     C(31)*S(2)**3*(S(1)**2 + S(3)**2)                                    
C                                                                               
      Y = C(32)*S(2)**3*S(1)*S(3) + C(33)*S(2)**4*(S(1) + S(3)) +               
     *    C(34)*S(2)**5 + C(35)*(S(1)**6 + S(3)**6) +                           
     *    C(36)*(S(1)**5*S(3) + S(1)*S(3)**5) +                                 
     *    C(37)*(S(1)**4*S(3)**2 + S(1)**2*S(3)**4) +                           
     *    C(38)*S(1)**3*S(3)**3                                                 
C                                                                               
      Z = C(39)*S(2)*(S(1)**5 + S(3)**5) +                                      
     *    C(40)*S(2)*(S(1)**4*S(3) + S(1)*S(3)**4) +                            
     *    C(41)*S(2)*(S(1)**3*S(3)**2 + S(1)**2*S(3)**3) +                      
     *    C(42)*S(2)**2*(S(1)**4 + S(3)**4) +                                   
     *    C(43)*S(2)**2*(S(1)**3*S(3) + S(1)*S(3)**3) +                         
     *    C(44)*S(1)**2*S(2)**2*S(3)**2                                         
C                                                                               
      ZZ = C(45)*S(2)**3*(S(1)**3 + S(3)**3) +                                  
     *     C(46)*S(2)**3*(S(1)**2*S(3) + S(1)*S(3)**2) +                        
     *     C(47)*S(2)**4*(S(1)**2 + S(3)**2) +                                  
     *     C(48)*S(2)**4*S(1)*S(3) +                                            
     *     C(49)*S(2)**5*(S(1) + S(3)) + C(50)*S(2)**6                          
C                                                                               
      P = T + TT + WW + W + GG + Y + Z + ZZ                                     
C                                                                               
      TEMP = EXP(-GAM*(X(1) + X(3))**2)                                         
C                                                                               
      U = U + P*Q*EXP(-GAM*(X(1) + X(3))**2)                                    
C                                                                               
      XTER=U + EZERO(1)                                                         
      V(IC) = XTER/HART                                                         
C                                                                               
      IF(IC .LT. 4)GO TO 97                                                     
      IF(R(1).GT.30.D0)GO TO 77                                                 
      IF(R(2).LT.30.D0)GO TO 66                                                 
C                                                                               
 77   CONTINUE                                                                  
      IF(NDER.EQ.1) THEN                                                        
         DEGSDR(1)=(V(3) - V(4))/STEP                                           
         DEGSDR(2)=(V(2) - V(4))/STEP                                           
         DEGSDR(3)=(V(1) - V(4))/STEP                                           
      ENDIF                                                                     
C                                                                               
 98   IF(R(2).LT.20)GO TO 97                                                    
C                                                                               
      V(4)=V(4)*1.0D8                                                           
      V(4)=INT(V(4))                                                            
      V(4)=V(4)*1.0D-8                                                          
C                                                                               
   97 IC=IC+1                                                                   
C                                                                               
      IF(IC.LT.5) GO TO 86                                                      
C                                                                               
      ENGYGS=V(4)                                                               
C                                                                               
      GO TO 113                                                                 
C                                                                               
 66   CONTINUE                                                                  
C                                                                               
C BEGINNING OF NDER SWITCH AROUND DERIVATIVES                                   
C                                                                               
      IF(NDER.EQ.1) THEN                                                        
C                                                                               
      DA = C(2) + 2.0D0*S(1)*C(4) + S(2)*C(6) + S(3)*C(7) +                     
     *     3.0D0*S(1)**2*C(8) + (S(3)**2 + 2.0D0*S(1)*S(3))*C(10) +             
     *     S(2)**2*C(11) + 2.0D0*S(1)*S(2)*C(12) + S(2)*S(3)*C(13) +            
     *     4.0D0*S(1)**3*C(14)                                                  
C                                                                               
      DB = (3.0D0*S(1)**2*S(3) + S(3)**3)*C(16) +                               
     *      2.0D0*S(1)*S(3)**2*C(17) + 3.0D0*S(1)**2*S(2)*C(18) +               
     *      2.0D0*S(1)*S(2)**2*C(19) + S(2)**3*C(20) +                          
     *      S(2)**2*S(3)*C(21) + (2.0D0*S(1)*S(2)*S(3) +                        
     *      S(2)*S(3)**2)*C(22)                                                 
C                                                                               
      DC = 5.0D0*S(1)**4*C(23) + (4.0D0*S(1)**3*S(3) +                          
     *     S(3)**4)*C(24) + (3.0D0*S(1)**2*S(3)**2 +                            
     *     2.0D0*S(1)*S(3)**3)*C(25) + 4.0D0*S(1)**3*S(2)*C(26) +               
     *     (3.0D0*S(2)*S(1)**2*S(3) + S(2)*S(3)**3)*C(27)                       
C                                                                               
      DD = 2.0D0*S(1)*S(2)*S(3)**2*C(28) +                                      
     *     3.0D0*S(1)**2*S(2)**2*C(29) + (2.0D0*S(1)*S(2)**2*S(3) +             
     *     S(2)**2*S(3)**2)*C(30) + 2.0D0*S(1)*S(2)**3*C(31) +                  
     *     S(2)**3*S(3)*C(32) + S(2)**4*C(33) +                                 
     *     6.0D0*S(1)**5*C(35) + (5.0D0*S(1)**4*S(3) +                          
     *     S(3)**5)*C(36)                                                       
C                                                                               
      DF = (4.0D0*S(1)**3*S(3)**2 + 2.0D0*S(1)*S(3)**4)*C(37) +                 
     *      3.0D0*S(1)**2*S(3)**3*C(38) + 5.0D0*S(1)**4*S(2)*C(39) +            
     *     (4.0D0*S(2)*S(1)**3*S(3) + S(2)*S(3)**4)*C(40) +                     
     *     (3.0D0*S(1)**2*S(2)*S(3)**2 +                                        
     *      2.0D0*S(1)*S(2)*S(3)**3)*C(41) +                                    
     *      4.0D0*S(1)**3*S(2)**2*C(42)                                         
C                                                                               
      DG = (3.0D0*S(1)**2*S(2)**2*S(3) + S(2)**2*S(3)**3)*C(43) +               
     *      2.0D0*S(1)*S(2)**2*S(3)**2*C(44) +                                  
     *      3.0D0*S(1)**2*S(2)**3*C(45) + (2.0D0*S(1)*S(2)**3*S(3) +            
     *      S(2)**3*S(3)**2)*C(46) + 2.0D0*S(1)*S(2)**4*C(47) +                 
     *      S(2)**4*S(3)*C(48) + S(2)**5*C(49)                                  
C                                                                               
      D1D = DA + DB + DC + DD + DF + DG                                         
C                                                                               
      DA2 = C(3) + 2.0D0*S(2)*C(5) + (S(1) + S(3))*C(6) +                       
     *      3.0D0*S(2)**2*C(9) + 2.0D0*S(2)*(S(1) + S(3))*C(11) +               
     *     (S(1)**2 + S(3)**2)*C(12) + S(1)*S(3)*C(13) +                        
     *      4.0D0*S(2)**3*C(15) + (S(1)**3 + S(3)**3)*C(18) +                   
     *      2.0D0*S(2)*(S(1)**2 + S(3)**2)*C(19)                                
C                                                                               
      DB2 = 3.0D0*S(2)**2*(S(1) + S(3))*C(20) +                                 
     *      2.0D0*S(1)*S(2)*S(3)*C(21) + (S(1)**2*S(3) +                        
     *      S(1)*S(3)**2)*C(22) + (S(1)**4 + S(3)**4)*C(26) +                   
     *     (S(1)**3*S(3) + S(1)*S(3)**3)*C(27) +                                
     *      S(1)**2*S(3)**2*C(28) + 2.0D0*S(2)*(S(1)**3 +                       
     *      S(3)**3)*C(29)                                                      
C                                                                               
      DC2 = 2.0D0*S(2)*(S(1)**2*S(3) + S(1)*S(3)**2)*C(30) +                    
     *      3.0D0*S(2)**2*(S(1)**2 + S(3)**2)*C(31) +                           
     *      3.0D0*S(2)**2*S(1)*S(3)*C(32) + 4.0D0*S(2)**3*(S(1) +               
     *      S(3))*C(33) + 5.0D0*S(2)**4*C(34)                                   
C                                                                               
      DD2 = (S(1)**5 + S(3)**5)*C(39) + (S(1)**4*S(3) +                         
     *       S(1)*S(3)**4)*C(40) + (S(1)**3*S(3)**2 +                           
     *       S(1)**2*S(3)**3)*C(41) + 2.0D0*S(2)*(S(1)**4 +                     
     *       S(3)**4)*C(42) + 2.0D0*S(2)*(S(1)**3*S(3) +                        
     *       S(1)*S(3)**3)*C(43) + 2.0D0*S(1)**2*S(2)*S(3)**2*C(44)             
C                                                                               
      DF2 = 3.0D0*S(2)**2*(S(1)**3 + S(3)**3)*C(45) +                           
     *      3.0D0*S(2)**2*(S(1)**2*S(3) + S(1)*S(3)**2)*C(46) +                 
     *      4.0D0*S(2)**3*(S(1)**2 + S(3)**2)*C(47) +                           
     *      4.0D0*S(2)**3*S(1)*S(3)*C(48) +                                     
     *      5.0D0*S(2)**4*(S(1) + S(3))*C(49) +                                 
     *     6.0D0*S(2)**5*C(50)                                                  
C                                                                               
      D2D = DA2 + DB2 + DC2 + DD2 + DF2                                         
C                                                                               
      DA3 = C(2) + 2.0D0*S(3)*C(4) + S(2)*C(6) + S(1)*C(7) +                    
     *      3.0D0*S(3)**2*C(8) + (2.0D0*S(1)*S(3) + S(1)**2)*C(10) +            
     *      S(2)**2*C(11) + 2.0D0*S(2)*S(3)*C(12) + S(1)*S(2)*C(13) +           
     *      4.0D0*S(3)**3*C(14)                                                 
C                                                                               
      DB3 = (S(1)**3 + 3.0D0*S(1)*S(3)**2)*C(16) +                              
     *       2.0D0*S(1)**2*S(3)*C(17) + 3.0D0*S(2)*S(3)**2*C(18) +              
     *       2.0D0*S(2)**2*S(3)*C(19) + S(2)**3*C(20) +                         
     *       S(2)**2*S(1)*C(21) + S(2)*(S(1)**2 +                               
     *       2.0D0*S(1)*S(3))*C(22) + 5.0D0*S(3)**4*C(23) + (S(1)**4 +          
     *       4.0D0*S(1)*S(3)**3)*C(24)                                          
C                                                                               
      DC3 = (2.0D0*S(1)**3*S(3) + 3.0D0*S(1)**2*S(3)**2)*C(25) +                
     *       4.0D0*S(2)*S(3)**3*C(26) + S(2)*(S(1)**3 +                         
     *       3.0D0*S(1)*S(3)**2)*C(27) +                                        
     *       2.0D0*S(2)*S(1)**2*S(3)*C(28) + 3.0D0*S(2)**2*S(3)**2*C(29)        
C                                                                               
      DD3 = S(2)**2*(S(1)**2 + 2.0D0*S(1)*S(3))*C(30) +                         
     *      2.0D0*S(2)**3*S(3)*C(31) + S(2)**3*S(1)*C(32) +                     
     *      S(2)**4*C(33) + 6.0D0*S(3)**5*C(35) + (S(1)**5 +                    
     *      5.0D0*S(1)*S(3)**4)*C(36) + (2.0D0*S(1)**4*S(3) +                   
     *      4.0D0*S(1)**2*S(3)**3)*C(37) + 3.0D0*S(1)**3*S(3)**2*C(38)          
C                                                                               
      DF3 = 5.0D0*S(2)*S(3)**4*C(39) + S(2)*(S(1)**4 +                          
     *      4.0D0*S(1)*S(3)**3)*C(40) + S(2)*(2.0D0*S(1)**3*S(3) +              
     *      3.0D0*S(1)**2*S(3)**2)*C(41) +                                      
     *      4.0D0*S(2)**2*S(3)**3*C(42) + S(2)**2*(S(1)**3 +                    
     *      3.0D0*S(1)*S(3)**2)*C(43) +                                         
     *      2.0D0*S(1)**2*S(2)**2*S(3)*C(44) +                                  
     *      3.0D0*S(2)**3*S(3)**2*C(45) + S(2)**3*(S(1)**2 +                    
     *      2.0D0*S(1)*S(3))*C(46) + 2.0D0*S(2)**4*S(3)*C(47) +                 
     *      S(2)**4*S(1)*C(48) + S(2)**5*C(49)                                  
C                                                                               
      D3D = DA3 + DB3 + DC3 + DD3 + DF3                                         
C                                                                               
      DQ1=(1.0D0 - TANH(AL(2)*S(2)))*(1.0D0-TANH(AL(3)*S(3)))                   
      DQ2=(1.0D0 - TANH(AL(1)*S(1)))*(1.0D0-TANH(AL(3)*S(3)))                   
      DQ3=(1.0D0 - TANH(AL(1)*S(1)))*(1.0D0-TANH(AL(2)*S(2)))                   
C                                                                               
      DAP1=-12.0D0*AP/R(1)                                                      
      DAP2=-12.0D0*AP/(R(2) - 0.2D0)                                            
      DAP3=-12.0D0*AP/R(3)                                                      
C                                                                               
      DB1=-4.0D0*B/R(1)**5                                                      
      DB3=-4.0D0*B/R(3)**5                                                      
C                                                                               
C      IF(DCOSH(AL(1)*S(1)).LT.1.0D-15)DCRS1=1.0D30                             
C      IF(DCOSH(AL(1)*S(1)).LT.1.0D-15)GO TO 144                                
C                                                                               
C                                                                               
       DCSR1=-AL(1)/(COSH(AL(1)*S(1))**2)                                       
C                                                                               
C144   IF(COSH(AL(2)*S(2)).LT.1.0D-15)DCRS2=1.0D30                              
C      IF(COSH(AL(2)*S(2)).LT.1.0D-15)GO TO 145                                 
C                                                                               
      DCSR2=-AL(2)/(COSH(AL(2)*S(2))**2)                                        
C                                                                               
C145   IF(COSH(AL(3)*S(3)).LT.1.0D-15)DCRS3=1.0D30                              
C      IF(COSH(AL(3)*S(3)).LT.1.0D-15)GO TO 146                                 
C                                                                               
      DCSR3=-AL(3)/(COSH(AL(3)*S(3))**2)                                        
                                                                                
 146  PSI13=2.0D0*(-GAM)*TEMP*(R(1) + R(3))                                     
C                                                                               
      DCHI1=2.0D0*ALF*R(1)*EXP(-ALF*(R(1)**2))                                  
      DCHI2=2.0D0*ALF*R(2)*EXP(-ALF*(R(2)**2))                                  
      DCHI3=2.0D0*ALF*R(3)*EXP(-ALF*(R(3)**2))                                  
C                                                                               
C      IF(VM(1).GT.1.0D20)DVR1=1.0D30                                           
C      IF(VM(1).GT.1.0D20)GO TO 153                                             
      DVR1 = 2.0D0*DE(1)*(1.0D0 - EXP(-BET(1)*(R(1)-XRE(1))))*                  
     *       BET(1)*EXP(-BET(1)*(R(1)-XRE(1)))                                  
C                                                                               
C 153  IF(VM(2).GT.1.0D20)DVR2=1.0D30                                           
C      IF(VM(2).GT.1.0D20)GO TO 154                                             
C                                                                               
      DVR2 = 2.0D0*DE(2)*(1.0D0 - EXP(-BET(2)*(R(2)-XRE(2))))*                  
     *       BET(2)*EXP(-BET(2)*(R(2)-XRE(2)))                                  
C                                                                               
C 154  IF(VM(3).GT.1.0D20)DVR3=1.0D30                                           
C      IF(VM(3).GT.1.0D20)GO TO 155                                             
C                                                                               
      DVR3 = 2.0D0*DE(3)*(1.0D0 - EXP(-BET(3)*(R(3)-XRE(3))))*                  
     *       BET(3)*EXP(-BET(3)*(R(3)-XRE(3)))                                  
C                                                                               
 155  DEGSDR(1) = VM(2)*CHI(3)*DCHI1 + CHI(2)*CHI(3)*DVR1 +                     
     *          VM(3)*CHI(2)*DCHI1 + D1D*Q*TEMP + P*DQ1*TEMP*DCSR1 +            
     *          DAP1 + DB1 + P*Q*PSI13                                          
      DEGSDR(2) = DVR2*CHI(1)*CHI(3) + VM(1)*CHI(3)*DCHI2 +                     
     *          VM(3)*CHI(1)*DCHI2 + D2D*Q*TEMP + P*DQ2*TEMP*DCSR2 +            
     *          DAP2                                                            
      DEGSDR(3) = VM(2)*CHI(1)*DCHI3 + VM(1)*CHI(2)*DCHI3 +                     
     *          DVR3*CHI(1)*CHI(2) + D3D*Q*TEMP + P*DQ3*TEMP*DCSR3 +            
     *          P*Q*PSI13 + DAP3 + DB3                                          
C                                                                               
      DO 222 I=1,3                                                              
         DEGSDR(I) = DEGSDR(I)/HART*BOHR                                        
 222  continue                                                                  
C                                                                               
C END OF NDER SWITCH AROUND DERIVATIVES                                         
C                                                                               
      ENDIF                                                                     
C                                                                               
      J=J+1                                                                     
      GO TO 98                                                                  
C                                                                               
  113 CONTINUE                                                                  
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the OH2 ',                    
     *                'Schinke-Lester potential energy surface')                
601   FORMAT(/,2X,T5,'Step size for the numerical derivatives = ',E10.5)        
602   FORMAT(2X,T5,'The Schinke-Lester surface has analytical ',                
     *             'derivatives if R(1) and R(2) are ',                         
     *       /,2X,T5,'less than 50 bohr, otherwise the derivatives ',           
     *               'are computed numerically.')                               
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,                       
     *       /, 2X, 'This value of NDER is not allowed in this ',               
     *              'version of the potential.')                                
C                                                                               
       stop                                                                     
       END                                                                      
C                                                                               
         BLOCK DATA PTPACM                                                      
         IMPLICIT REAL*8 (A-H,O-Z)                                              
C         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   
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
         PARAMETER (BOHR = 0.529177D0)                                          
         PARAMETER (HART = 27.21161D0)                                          
C                                                                               
         COMMON /SATOCM/ DE(3), XRE(3), BET(3)                                  
         COMMON /NDERCM/ STEP                                                   
         COMMON /PRECM/  C(50),A,B,ALF,GAM,AL(3),                               
     +                   X(3),S(3),VM(3),CHI(3),V(4)                            
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
C   Initialize the potential parameters; the energy parameters                  
C   are in eV and the lengths are in Angstroms.                                 
C                                                                               
         DATA DE  /4.035D0,4.1662D0,4.035D0   /                                 
         DATA XRE /0.9923D0,0.7590D0,0.9923D0 /                                 
         DATA BET /2.4113D0,2.0682D0, 2.4113D0/                                 
C                                                                               
C   Initialize the step size used for the numerical derivatives                 
C   of the energy with respect to the coordinates.                              
C                                                                               
      DATA STEP/5.0D-7/                                                         
C                                                                               
      DATA A   / 0.157002D-03/                                                  
      DATA B   /-0.162301D-01/                                                  
      DATA ALF / 0.9D0       /                                                  
      DATA GAM / 0.045D0     /                                                  
      DATA AL  / 1.0D0,0.8D0,1.0D0/                                             
C                                                                               
      DATA C/0.354195D+01, 0.302780D+01, 0.577312D+00, 0.539877D+01,            
     +       0.233027D+01, 0.157276D+02, 0.701851D+01,-0.478535D+02,            
     +      -0.106915D+01, 0.268137D+02,-0.104806D+02,-0.182660D+02,            
     +       0.833941D+01, 0.872497D+02, 0.917870D+01,-0.697766D+02,            
     +       0.434862D+02,-0.156633D+02, 0.380162D+02,-0.796542D+01,            
     +       0.602854D+02, 0.500800D+02,-0.545761D+02, 0.410210D+02,            
     +      -0.890739D+01, 0.191135D+02,-0.512653D+02,                          
     +       0.365096D+02,-0.235972D+02,-0.102322D+03, 0.837770D+00,            
     +       0.258200D+02,-0.198105D+01,-0.484189D+01, 0.105377D+02,            
     +       0.372352D+01,-0.207237D+02, 0.223342D+02, 0.244787D+01,            
     +       0.979572D+01,-0.135612D+01,-0.724496D+01, 0.234984D+02,            
     +       0.136364D+03, 0.253799D+01,-0.203821D+02, 0.699846D+01,            
     +       0.198529D+02,-0.475902D+01, 0.236328D+01/                          
C                                                                               
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
