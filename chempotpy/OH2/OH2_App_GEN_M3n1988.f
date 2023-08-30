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
C   System:          OH2   A'' surface                                          
C   Common name:     M3n, symmetrized M2 potential                              
C   Functional form: Rotated-Morse-oscillator-spline (RMOS) potential           
C                    for collinear geometeries plus anti-Morse                  
C                    bending potential. Gamma is a function of the              
C                    collinear OH distance.                                     
C   Reference:       unpublished                                                
C   Cross reference: T. Joseph, D. G. Truhlar, and B. C. Garrett                
C                    J. Chem. Phys. 88, 6982-6990 (1988).                       
C                                                                               
C   Interface:       potlib2001
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
C
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in the block data subprogram PTPACM.  
C   Coordinates, potential energy, and derivatives are passed                   
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in hartree atomic units.                                                 
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(O-H)                                           
C                       R(2) = R(H-H)                                           
C                       R(3) = R(H-O)                                           
C   The classical potential energy is set equal to zero for the O               
C   infinitely far from the equilibrium H2 diatomic.                            
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
      COMMON /AMBPCM/ D(3),BETA,REQ(3),A,GAMP(3),RBB,RHH,DE                     
C                                                                               
         PARAMETER (ONE = 1.0D0)                                                
C                                                                               
         COMMON /PRECM/  DE1(3),DE2(3),DCOSDR(3),DFDR(3),IOPSW                  
C                                                                               
C   Echo the name of the potential energy surface                               
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
C      WRITE (NFLAG(18), 600)                                                    
C                                                                               
      HPI = 2.0D0*ATAN(ONE)                                                     
C                                                                               
      IOPSW = MAX(0, IOPSW)                                                     
C                                                                               
C      WRITE (NFLAG(18), 601) IOPSW                                              
C                                                                               
C      CALL PREM2                                                               
      CALL PREPOT2                                                              
C                                                                               
C   Close the file containing the potential parameters                          
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the A'''' OH2 ',              
     *                'potential energy surface M3n',                           
     *       //,2X,T5,'Potential energy surface parameters:')                   
601   FORMAT (/,2X,T5,'Switching function parameter, n=', T50,I5,               
     *        /,2X,T5,'Switching function, Fn(chi) = ',                         
     *                'sin**2((pi/2)*Fn-1(chi))',                               
     *        /,2X,T5,'where F0(chi) = sin**2(chi/2)')                          
C                                                                               
      EZERO(1)=D(2)                                                             
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='T. Joseph, D. G. Truhlar, B. C. Garrett,'                        
       REF(2)='J. Chem. Phys. 88, 6982(1988)'                                   
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
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(O-H)                                           
C                       R(2) = R(H-H)                                           
C                       R(3) = R(H-O)                                           
C   The classical potential energy is set equal to zero for the O               
C   infinitely far from the equilibrium H2 diatomic.                            
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
C   Conversion for Angstroms to Bohr (Angstroms/Bohr)                           
C                                                                               
      PARAMETER (BOHR = .52917706D0)                                            
C                                                                               
C   Conversion for kcal to Hartree (kcal/Hartree)                               
C                                                                               
      PARAMETER (HARTRE = 627.5095D0)                                           
C                                                                               
         PARAMETER (TOL = 1.D-5)                                                
         PARAMETER (ONE = 1.0D0)                                                
C                                                                               
         COMMON /PRECM/  DE1(3),DE2(3),DCOSDR(3),DFDR(3),IOPSW                  
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
         HPI = 2.0D0*ATAN(ONE)                                                  
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
C  CHECK FOR SMALL VALUES OF R                                                  
C                                                                               
      IF (R(1) .LT. 0.001D0 .OR. R(2). LT. 0.001D0 .OR.                         
     *    R(3) .LT. 0.001D0) THEN                                               
         ENGYGS = 1.D10                                                         
         IF (NDER .NE. 1) THEN                                                  
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
             RETURN                                                             
         END IF                                                                 
         DEGSDR(1) = -ENGYGS                                                    
         DEGSDR(2) = -ENGYGS                                                    
         DEGSDR(3) = -ENGYGS                                                    
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
         RETURN                                                                 
      END IF                                                                    
C                                                                               
C    EVALUATE SWITCHING FUNCTION AND ITS DERIVATIVE                             
C                                                                               
      R12 = R(1)*R(1)                                                           
      R22 = R(2)*R(2)                                                           
      R32 = R(3)*R(3)                                                           
      BIGR2 = 0.5D0*(R12 - 0.5D0*R22 + R32)                                     
      IF (BIGR2 .LE. 1.D-8) THEN                                                
         F = 0.5D0                                                              
         DO 3 I = 1,3                                                           
    3       DFDR(I) = 0.0D0                                                     
      ELSE                                                                      
         BIGR = SQRT(BIGR2)                                                     
         COSCHI = (R12 - R32)/(2.0D0*R(2)*BIGR)                                 
C                                                                               
C    DERIVATIVES OF COS                                                         
C                                                                               
         T = 0.5D0/(R(2)*BIGR2)                                                 
         IF (NDER .EQ. 1) DCOSDR(2) = -COSCHI*(R12 - R22 + R32)*T               
         T = 0.5D0*T/BIGR                                                       
         IF (NDER .EQ. 1) THEN                                                  
             DCOSDR(1) = R(1)*(R12 - R22 + 3.0D0*R32)*T                         
             DCOSDR(3) = -R(3)*(3.0D0*R12 - R22 + R32)*T                        
         ENDIF                                                                  
C                                                                               
C    EVALUATE SWITCHING FUNCTION                                                
C                                                                               
         G = COSCHI                                                             
         DGDX = 1.0D0                                                           
         IF (IOPSW .GT. 0) THEN                                                 
            DO 20 I = 1,IOPSW                                                   
               T = HPI*G                                                        
               G = SIN(T)                                                       
               IF (NDER .EQ. 1) DGDX = HPI*COS(T)*DGDX                          
   20       CONTINUE                                                            
         END IF                                                                 
         F = 0.5D0*(1.0D0 - G)                                                  
         IF (NDER .EQ. 1) THEN                                                  
             DFDX = -0.5D0*DGDX                                                 
             DO 30 I = 1,3                                                      
                   DFDR(I) = DFDX*DCOSDR(I)                                     
30           CONTINUE                                                           
         ENDIF                                                                  
      END IF                                                                    
      OMF = 1.0D0 - F                                                           
      IF (F .GT. 1.D-10) THEN                                                   
C                                                                               
C    EVALUATE POTENTIAL FOR O + H2 DIRECTION                                    
C                                                                               
C         CALL POTM2                                                            
         CALL POT2                                                              
         V1 = ENGYGS                                                            
         IF (NDER .EQ. 1) THEN                                                  
             DO 50 I = 1,3                                                      
50                 DE1(I) = DEGSDR(I)                                           
         ENDIF                                                                  
      ELSE                                                                      
         V1 = 0.0D0                                                             
         IF (NDER .EQ. 1) THEN                                                  
             DO 55 I = 1,3                                                      
55                 DE1(I) = 0.0D0                                               
         ENDIF                                                                  
      END IF                                                                    
      IF (ABS(F-0.5D0) .LT. 1.D-10) THEN                                        
         V2 = ENGYGS                                                            
         IF (NDER .EQ. 1) THEN                                                  
             DE2(1) = DEGSDR(3)                                                 
             DE2(2) = DEGSDR(2)                                                 
             DE2(3) = DEGSDR(1)                                                 
         ENDIF                                                                  
      ELSE IF (OMF .GT. 1.D-10) THEN                                            
C                                                                               
C    INTERCHANGE R(1) AND R(3) AND RECOMPUTE POTENTIAL                          
C                                                                               
         T = R(1)                                                               
         R(1) = R(3)                                                            
         R(3) = T                                                               
C         CALL POTM2                                                            
         CALL POT2                                                              
         V2 = ENGYGS                                                            
         IF (NDER .EQ. 1) THEN                                                  
             DE2(1) = DEGSDR(3)                                                 
             DE2(2) = DEGSDR(2)                                                 
             DE2(3) = DEGSDR(1)                                                 
         ENDIF                                                                  
         T = R(1)                                                               
         R(1) = R(3)                                                            
         R(3) = T                                                               
      ELSE                                                                      
         V2 = 0.0D0                                                             
         IF (NDER .EQ. 1) THEN                                                  
             DO 60 I = 1,3                                                      
60                 DE2(I) = 0.0D0                                               
         ENDIF                                                                  
      END IF                                                                    
C                                                                               
      ENGYGS = F*V1 + OMF*V2                                                    
      IF (NDER .NE. 1) THEN                                                     
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
          RETURN                                                                
      END IF                                                                    
      DELV = V1 - V2                                                            
      DO 300 I = 1,3                                                            
         DEGSDR(I) = F*DE1(I) + OMF*DE2(I) + DELV*DFDR(I)                       
  300 CONTINUE                                                                  
C                                                                               
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/,2X,T5,'POT has been called with NDER = ',I5,                     
     *       /,2X,T5,'This value of NDER is not allowed in this ',              
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
      SUBROUTINE PREPOT2                                                        
C                                                                               
C   The subprogram PREPOT2 add the bending correction to the collinear          
C   O + H2 potential subroutine.                                                
C                                                                               
C   The potential parameters for the collinear potential are                    
C   initialized in the block data subprogram PTPACM.  The                       
C   potential parameters for the bending potential are assigned                 
C   in the subprogram PREPOT22.  The subprograms PREPOT2 and                    
C   PREPOT2 do not needed to be called by the main program.                     
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
      DOUBLE PRECISION NAB,NAB1,NBC,NBC1                                        
C                                                                               
C   Conversion for Angstroms to Bohr (Angstroms/Bohr)                           
C                                                                               
      PARAMETER (BOHR = .52917706D0)                                            
C                                                                               
C   Conversion for kcal to Hartree (kcal/Hartree)                               
C                                                                               
      PARAMETER (HARTRE = 627.5095D0)                                           
C                                                                               
C  Vaule of tolerance                                                           
C                                                                               
      PARAMETER (TOL = 1.D-5)                                                   
C                                                                               
      COMMON /AMBPCM/ D(3),BETA,REQ(3),A,GAMP(3),RBB,RHH,DE                     
      COMMON /PRECM/  DE1(3),DE2(3),DCOSDR(3),DFDR(3),IOPSW                     
      COMMON /PRE2CM/ DBCOOR(2)                                                 
C                                                                               
C   Assign parameters for the bending potential                                 
C                                                                               
C      CALL PREM22                                                              
C                                                                               
      CALL PREPOT22                                                             
C                                                                               
C   Echo the potential parameters for the bending correction.                   
C                                                                               
C      WRITE (NFLAG(18),600) REQ, D                                              
C      WRITE (NFLAG(18),601) BETA,A,RBB,GAMP                                     
C                                                                               
C       CONVERT TO ATOMIC UNITS                                                 
C                                                                               
      DO 50 IT = 1,3                                                            
         REQ(IT) = REQ(IT)/BOHR                                                 
         D(IT) = D(IT)/HARTRE                                                   
   50 CONTINUE                                                                  
      BETA = BETA * BOHR                                                        
      DE = D(1)                                                                 
      A = A/BOHR                                                                
      RHH = 0.74144D0 / BOHR                                                    
      RBB = RBB /BOHR                                                           
C                                                                               
C   Initialize the energy in the asymptotic valleys                             
C                                                                               
      EASYAB = D(1)                                                             
      EASYBC = D(2)                                                             
      EASYAC = D(3)                                                             
C                                                                               
600   FORMAT(/,2X,T5,'Parameters for the bending correction',                   
     *       /,2X,T5,'Bond', T47, 'O-H', T58, 'H-H', T69, 'H-O',                
     *       /,2X,T5,'Equilibrium bond lengths (Angstroms):',                   
     *       T44, F10.5, T55, F10.5, T66, F10.5,                                
     *       /,2X,T5,'Dissociation energy (kcal/mol):',                         
     *       T44, F10.5, T55, F10.5, T66, F10.5)                                
601   FORMAT(2X,T5,'Morse beta parameter (Angstroms**-1):', T44,1PE15.8,        
     *     /,2X,T5,'Pauling parameter (Angstroms):',T44,1PE15.8,                
     *     /,2X,T5,'RBB (Angstroms):',T44,1PE15.8,                              
     *     /,2X,T5,'Gamma for the bending correction:',T44,1PE15.8,             
     *          T60,1PE15.8,/,2X,T44,1PE15.8)                                   
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT2                                                           
C                                                                               
C      ENTRY POTM2                                                              
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
      DOUBLE PRECISION NAB,NAB1,NBC,NBC1                                        
C                                                                               
C   Conversion for Angstroms to Bohr (Angstroms/Bohr)                           
C                                                                               
      PARAMETER (BOHR = .52917706D0)                                            
C                                                                               
C   Conversion for kcal to Hartree (kcal/Hartree)                               
C                                                                               
      PARAMETER (HARTRE = 627.5095D0)                                           
C                                                                               
C  Vaule of tolerance                                                           
C                                                                               
      PARAMETER (TOL = 1.D-5)                                                   
C                                                                               
      COMMON /AMBPCM/ D(3),BETA,REQ(3),A,GAMP(3),RBB,RHH,DE                     
      COMMON /PRECM/  DE1(3),DE2(3),DCOSDR(3),DFDR(3),IOPSW                     
      COMMON /PRE2CM/ DBCOOR(2)                                                 
C                                                                               
      IF (R(1) .GT. TOL .AND. R(2) .GT. TOL .AND. R(3) .GT. TOL) THEN           
         RACCOL = R(1) + R(2)                                                   
         BCOOR = 0.0D0                                                          
         IF (NDER .EQ. 1) THEN                                                  
             DBCOOR(1) = 0.0D0                                                  
             DBCOOR(2) = 0.0D0                                                  
             DAC = 0.0D0                                                        
         ENDIF                                                                  
         ROB = 1.0D0                                                            
         IF(ABS(RHH-RBB).GT.1.D-10) THEN                                        
C                                                                               
C     FIRST CALCULATE SOME USEFUL NUMBERS                                       
C     NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS                          
C     AB AND BC RESPECTIVELY. (EQUATION 2.2)                                    
C                                                                               
            NAB = EXP((REQ(1) - R(1))/A)                                        
            IF (NAB.LT.1.D-15) NAB = 1.D-15                                     
            NBC = EXP((REQ(2) - R(2))/A)                                        
            IF (NBC.LT.1.D-15) NBC = 1.D-15                                     
C                                                                               
C       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE                     
C       SEE EQUATION 2.7A, 2.7B                                                 
C                                                                               
            C = 1.D0 - NAB/NBC                                                  
            IF (ABS(C) .LT. 1.D-14) THEN                                        
               NAB1=.5D0                                                        
               NBC1=.5D0                                                        
            ELSE                                                                
               C = C/NAB                                                        
               NAB1 = (2.D0 + C - SQRT(4.D0+C*C))/(2.D0*C)                      
               NBC1 = 1.D0 - NAB1                                               
            END IF                                                              
C                                                                               
C          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS                    
C          FIRST TRAP OUT ANY ZERO ARGUMENTS                                    
C                                                                               
            IF (NAB1*NBC1 .LE. 0.0D0) THEN                                      
               ROB = 1.D0                                                       
            ELSE                                                                
               STUFF = A * LOG(NAB1*NBC1)                                       
               T = RHH - STUFF                                                  
               ROB = 1.0D0                                                      
               IF(ABS(T).GT.1.D-15) ROB = (RBB - STUFF)/T                       
            END IF                                                              
         END IF                                                                 
C                                                                               
C        CALCULATE BENDING CORRECTION                                           
C        EVALUATE GAMMA AND DERIVATIVE                                          
C                                                                               
         EX = EXP(RACCOL*(GAMP(2) + GAMP(3)*RACCOL))                            
         GAMMA = GAMP(1) + EX                                                   
         IF (NDER .EQ. 1) DGAM = (GAMP(2) + 2.D0*GAMP(3)*RACCOL)*EX             
         GAMMA = GAMMA*DE                                                       
         IF (NDER .EQ. 1) DGAM = DGAM*DE                                        
C                                                                               
C           CALCULATE V(R(3)) TERM                                              
C                                                                               
         EX = EXP(-BETA*(R(3)-REQ(3))/ROB)                                      
         EX1 = EX*(1.D0 + 0.5D0*EX)                                             
         IF (NDER .EQ. 1) DEX1 = EX*(1.D0 + EX)                                 
C                                                                               
C           CALCULATE V(R(1)+R(2)) TERM                                         
C                                                                               
         EX = EXP(BETA*(REQ(3)-R(2)-R(1))/ROB)                                  
         EX2 = EX*(1.D0 + 0.5D0*EX)                                             
         IF (NDER .EQ. 1) DEX2 = EX*(1.D0 + EX)                                 
C                                                                               
C        HERE IS THE BENDING CORRECTION                                         
C                                                                               
         BCOOR = GAMMA*(EX1 - EX2)                                              
C                                                                               
C        WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE                        
C        OF BCOOR WITH RESPECT TO R(1)  R(2) AND R(3).                          
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))                    
C                                                                               
         GAMMA = GAMMA*BETA/ROB                                                 
         IF (NDER .EQ. 1) THEN                                                  
             DAC = -GAMMA*DEX1                                                  
C                                                                               
C          CALCULATE CORRECTIONS TO DAB, DBC                                    
C                                                                               
             DBCOOR(1) = GAMMA*DEX2 + DGAM*(EX1-EX2)                            
             DBCOOR(2)= DBCOOR(1)                                               
         ENDIF                                                                  
C                                                                               
C       NOW CALCULATE THE COLLINEAR ENERGY                                      
C       STORE R(3) IN RACTMP THEN SET R(3) TO                                   
C       THE COLLINEAR GEOMETRY AND CALL POTM22                                  
C                                                                               
         RACTMP = R(3)                                                          
         R(3) = R(2) + R(1)                                                     
C                                                                               
C         CALL POTM22                                                           
         CALL POT22                                                             
                                                                                
C                                                                               
         R(3) = RACTMP                                                          
C                                                                               
         ENGYGS = ENGYGS + BCOOR                                                
C                                                                               
         IF (NDER .EQ. 1) THEN                                                  
             DO 300 IT = 1,2                                                    
                    DEGSDR(IT) = DBCOOR(IT) + DEGSDR(IT) + DEGSDR(3)            
300          CONTINUE                                                           
             DEGSDR(3) = DAC                                                    
         ENDIF                                                                  
      ELSE                                                                      
        ENGYGS = 1.D30                                                          
        IF (NDER .NE. 1) RETURN                                                 
        DEGSDR(1) = 0.0D0                                                       
        IF (R(1) .LE. TOL) DEGSDR(1) = -ENGYGS                                  
        DEGSDR(2) = 0.0D0                                                       
        IF (R(1) .LE. TOL) DEGSDR(2) = -ENGYGS                                  
        DEGSDR(3) = 0.0D0                                                       
        IF (R(1) .LE. TOL) DEGSDR(3) = -ENGYGS                                  
      END IF                                                                    
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE PREPOT22                                                       
C                                                                               
C     COLLINEAR A + BC POTENTIAL FUNCTION                                       
C     RMO-SPLINE TYPE FIT                                                       
C     COORDINATES AND ENERGIES PASSED TO POTENTIAL IN ATOMIC UNITS              
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
      PARAMETER (DETRAD = 0.0174532D0)                                          
      PARAMETER (CKCAL = 627.5095D0)                                            
      PARAMETER (BOHR = 0.529177D0)                                             
C                                                                               
      COMMON /PRECM/  DE1(3),DE2(3),DCOSDR(3),DFDR(3),IOPSW                     
      COMMON /PRE2CM/ DBCOOR(2)                                                 
      COMMON /PRE22CM/ D1INF,R1INF,B1INF,D2INF,                                 
     +                 R2INF,B2INF,A1INF,A2INF,                                 
     +                 R1S,R2S,VSP,ALPH1,ALPH2,                                 
     +                 PHI(20),BSP(20,3),CSP(22,3),                             
     +                 D1DIF,R1DIF,B1DIF,A1DIF,                                 
     +                 D2DIF,R2DIF,B2DIF,A2DIF,                                 
     +                 DASY,IAYES,NPHI                                          
C                                                                               
      CHARACTER*20 TITLE                                                        
C                                                                               
C  TITLE                                                                        
C                                                                               
      TITLE='   O+H2 MAB 9-12-79 '                                              
C                                                                               
c      WRITE (NFLAG(18), 599)                                                    
C                                                                               
C  ASSIGN COLLINEAR PARAMETERS (UNITS OF KCAL, ANGSTROM, DEGREES)               
C  IN BLOCK DATA SUBROUTINE                                                     
C                                                                               
      A1INF = 0.0D0                                                             
      A2INF = 0.0D0                                                             
      IAYES = 0                                                                 
C                                                                               
C  ASSIGN PARAMETERS FOR THE ASYMPTOTIC REACTANT AND PRODUCT REGIONS            
C  IN BLOCK DATA SUBROUTINE                                                     
C                                                                               
      IF(A2INF .EQ.0.0D0) A2INF = D1INF-D2INF                                   
C                                                                               
C  ASSIGN THE PIVOT POINT, POTENTIAL AT THE PIVOT,                              
C  FLAG, AND ASYMPTOTIC EXPONENTIAL PARAMETERS                                  
C  IN BLOCK DATA SUBROUTINE.  ECHO.                                             
                                                                                
c      WRITE (NFLAG(18), 600) TITLE                                              
c      WRITE (NFLAG(18), 601) D1INF,R1INF,B1INF,D2INF,R2INF,B2INF,A1INF,A        
C     2INF                                                                      
c      WRITE (NFLAG(18), 602) R1S,R2S,VSP,ALPH1,ALPH2                            
C                                                                               
c      IF(IAYES.NE.0) WRITE (NFLAG(18), 603)                                     
c      IF(IAYES.EQ.0) WRITE (NFLAG(18), 604)                                     
C                                                                               
      NPHIP2 = NPHI + 2                                                         
C                                                                               
C  ECHO POTENTIAL PARAMTERS ASSIGNED BY BLOCK DATA SUBROUTINE                   
C                                                                               
c      WRITE (NFLAG(18), 605) NPHI                                               
c      WRITE (NFLAG(18), 610)                                                    
C                                                                               
c      DO 30 I = 1,NPHI                                                          
c         WRITE (NFLAG(18), 606) I,PHI(I),(BSP(I,J),CSP(I,J),J=1,3)              
c   30 CONTINUE                                                                  
C                                                                               
c      WRITE (NFLAG(18), 611)                                                    
c      DO 40 I=1,2                                                               
c         K = I + NPHI                                                           
c         WRITE (NFLAG(18), 607) K,(CSP(K,J),J=1,3)                              
c   40 CONTINUE                                                                  
C                                                                               
C  CONVERT TO ATOMIC UNITS                                                      
C                                                                               
      D1INF = D1INF/CKCAL                                                       
      A1INF = A1INF/CKCAL                                                       
      R1INF = R1INF/BOHR                                                        
      B1INF = B1INF*BOHR                                                        
      D2INF = D2INF/CKCAL                                                       
      A2INF = A2INF/CKCAL                                                       
      R2INF = R2INF/BOHR                                                        
      B2INF = B2INF*BOHR                                                        
      R1S = R1S/BOHR                                                            
      R2S = R2S/BOHR                                                            
      VSP = VSP/CKCAL                                                           
      ALPH1 = ALPH1*BOHR                                                        
      ALPH2 = ALPH2*BOHR                                                        
      DO 50 I = 1,NPHI                                                          
         BSP(I,1) = BSP(I,1)/CKCAL                                              
         BSP(I,2) = BSP(I,2)/BOHR                                               
         BSP(I,3) = BSP(I,3)*BOHR                                               
   50 CONTINUE                                                                  
      NNN = NPHI + 2                                                            
      DO 60 I = 1,NNN                                                           
         CSP(I,1) = CSP(I,1)/CKCAL                                              
         CSP(I,2) = CSP(I,2)/BOHR                                               
         CSP(I,3) = CSP(I,3)*BOHR                                               
   60 CONTINUE                                                                  
C                                                                               
C  CONVERT FROM DEGREES TO RADIANS                                              
C                                                                               
      T = 1.0D0/DETRAD                                                          
      DO 70 J = 1,3                                                             
        CSP(2,J) = CSP(2,J)*T                                                   
   70 CONTINUE                                                                  
      T = T**3                                                                  
      DO 80 II = 1,NPHI                                                         
         PHI(II) = PHI(II)*DETRAD                                               
         I = II + 2                                                             
         DO 80 J = 1,3                                                          
            CSP(I,J) = CSP(I,J)*T                                               
   80 CONTINUE                                                                  
C                                                                               
C  SET UP CONSTANTS                                                             
C                                                                               
C     REACTANT ASYMPTOTE                                                        
C                                                                               
      DIS = BSP(1,1)                                                            
      REQ2 = BSP(1,2)                                                           
      BET = BSP(1,3)                                                            
      T1 = 1.0D0 - EXP(-BET*REQ2)                                               
      AER = VSP - DIS*T1*T1                                                     
      REQ2 = R2S - REQ2                                                         
      D1DIF = DIS - D1INF                                                       
      R1DIF = REQ2 - R1INF                                                      
      B1DIF = BET - B1INF                                                       
      A1DIF = AER - A1INF                                                       
C                                                                               
C     PRODUCT ASYMPTOTE                                                         
C                                                                               
      DIS = BSP(NPHI,1)                                                         
      REQ2 = BSP(NPHI,2)                                                        
      BET = BSP(NPHI,3)                                                         
      T1 = 1.0D0 - EXP(-BET*REQ2)                                               
      AER = VSP - DIS*T1*T1                                                     
      REQ2 = R1S - REQ2                                                         
      D2DIF = DIS - D2INF                                                       
      R2DIF = REQ2 - R2INF                                                      
      B2DIF = BET - B2INF                                                       
      A2DIF = AER - A2INF                                                       
C                                                                               
C  DEFINE ASYMPTOTIC ENERGY AND ZERO OF ENERGY                                  
C                                                                               
      DASY = D1INF                                                              
C                                                                               
599   FORMAT(/,2X,T5,'Parameters for the collinear part of the ',               
     *               'potential')                                               
600   FORMAT(2X,T5,'Title in potential data file: ',A20)                        
601   FORMAT(2X,T5,'D1INF, R1INF, B1INF, D2INF, R2INF, B2INF, ',                
     *               'A1INF, A2INF (kcal/mol):',                                
     *       /,(2X,T15,F10.5,T30,F10.5,T45,F10.5,T60,F10.5))                    
602   FORMAT(2X,T5,'R1S, R2S, VSP, ALPH1, ALPH2 (ANG/KCAL):',                   
     *       /,(2X,T15,F10.5,T30,F10.5,T45,F10.5))                              
603   FORMAT(2X,T5,'Note: the pivot point is fixed at VSP')                     
604   FORMAT(2X,T5,'Note: the pivot point is not fixed; VSP is ',               
     *               'the 3-body break-up energy')                              
605   FORMAT(2X,T5,'NPHI = ',I5)                                                
610   FORMAT(2X,T13,'I',T25,'PHI',T39,'BSP',T56,'CSP')                          
606   FORMAT(2X,T11,I5,T20,F10.5,(T35,F10.5,T50,F15.5))                         
611   FORMAT(2X,T11,'K',T20,'CSP')                                              
607   FORMAT(2X,T11,I5,(T20,E15.5,T40,E15.5))                                   
  800 FORMAT (A20)                                                              
  801 FORMAT (8F9.5)                                                            
  802 FORMAT (3F10.5,I5,2F10.5)                                                 
  803 FORMAT (16I5)                                                             
  804 FORMAT (4G20.10)                                                          
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT22                                                          
C                                                                               
C      ENTRY POTM22                                                             
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
      PARAMETER (DETRAD = 0.0174532D0)                                          
      PARAMETER (CKCAL = 627.5095D0)                                            
      PARAMETER (BOHR = 0.529177D0)                                             
C                                                                               
      COMMON /PRECM/  DE1(3),DE2(3),DCOSDR(3),DFDR(3),IOPSW                     
      COMMON /PRE2CM/ DBCOOR(2)                                                 
      COMMON /PRE22CM/ D1INF,R1INF,B1INF,D2INF,                                 
     +                 R2INF,B2INF,A1INF,A2INF,                                 
     +                 R1S,R2S,VSP,ALPH1,ALPH2,                                 
     +                 PHI(20),BSP(20,3),CSP(22,3),                             
     +                 D1DIF,R1DIF,B1DIF,A1DIF,                                 
     +                 D2DIF,R2DIF,B2DIF,A2DIF,                                 
     +                 DASY,IAYES,NPHI                                          
C                                                                               
      CHARACTER*20 TITLE                                                        
C                                                                               
C      R1 = R(1)                                                                
C      R2 = R(2)                                                                
C      R3 = R(3)                                                                
      IF (R(1) .GT. R1S) THEN                                                   
         IF (R(2) .GT. R2S) THEN                                                
C                                                                               
C   3-ATOM BREAK UP REGION                                                      
C                                                                               
            ENGYGS = DASY                                                       
            IF (NDER .EQ. 1) THEN                                               
                DEGSDR(1) = 0.0D0                                               
                DEGSDR(2) = 0.0D0                                               
            ENDIF                                                               
         ELSE                                                                   
C                                                                               
C   REACTANT ASYMPTOTIC REGION                                                  
C                                                                               
            DR1 = R(1) - R1S                                                    
            IF (DR1*ALPH1 .LT. -70.0D0)                                         
     +         WRITE (NFLAG(18), 6601) R(1), R(2), R(3),                        
     *         DR1, ALPH1                                                       
            T1 = EXP(-DR1*ALPH1)                                                
            T2 = -ALPH1*T1                                                      
C                                                                               
C.....INTERPOLATE FOR MORSE PARAMETERS                                          
C                                                                               
            DIS = D1INF + D1DIF*T1                                              
            DDIS = D1DIF*T2                                                     
            REQ2 = R1INF + R1DIF*T1                                             
            DREQ = R1DIF*T2                                                     
            BET = B1INF + B1DIF*T1                                              
            DBET = B1DIF*T2                                                     
            IF(IAYES.EQ.0) THEN                                                 
               AER = VSP - DIS                                                  
               DAER = -DDIS                                                     
            ELSE                                                                
               AER = A1INF + A1DIF*T1                                           
               DAER = A1DIF*T2                                                  
            END IF                                                              
            DELR = R(2) - REQ2                                                  
            IF (DELR*BET .LT. -70.0D0)                                          
     *          WRITE (NFLAG(18), 6602) R(1), R(2), R(3), DELR, BET             
            T1 = EXP(-BET*DELR)                                                 
            T2 = 1.0D0 - T1                                                     
            ENGYGS = AER + DIS*T2*T2                                            
            IF (NDER .EQ. 1) THEN                                               
               T3 = 2.0D0*DIS*T1                                                
               DEGSDR(1) =                                                      
     +               DAER + T2*(DDIS*T2 - T3*(BET*DREQ - DBET*DELR))            
               DEGSDR(2) = T3*BET*T2                                            
            ENDIF                                                               
         END IF                                                                 
      ELSE IF (R(2) .GT. R2S) THEN                                              
C                                                                               
C   PRODUCT ASYMPTOTIC REGION                                                   
C                                                                               
         DR2 = R(2) - R2S                                                       
         IF (DR2*ALPH2 .LT. -70.D0)                                             
     *       WRITE (NFLAG(18), 6603) R(1), R(2), R(3), DR2, ALPH2               
C                                                                               
C.....INTERPOLATE FOR MORSE PARAMETERS                                          
C                                                                               
         T1 = EXP(-DR2*ALPH2)                                                   
         T2 = -ALPH2*T1                                                         
         DIS = D2INF + D2DIF*T1                                                 
         DDIS = D2DIF*T2                                                        
         REQ2 = R2INF + R2DIF*T1                                                
         DREQ = R2DIF*T2                                                        
         BET = B2INF + B2DIF*T1                                                 
         DBET = B2DIF*T2                                                        
         IF(IAYES.EQ.0) THEN                                                    
            AER = VSP - DIS                                                     
            DAER = -DDIS                                                        
         ELSE                                                                   
            AER = A2INF + A2DIF*T1                                              
            DAER = A2DIF*T2                                                     
         END IF                                                                 
         DELR = R(1) - REQ2                                                     
         IF (DELR*BET .LT. -70.0D0)                                             
     *       WRITE (NFLAG(18), 6604) R(1), R(2), R(3), DELR, BET                
         T1 = EXP(-BET*DELR)                                                    
         T2 = 1.0D0 - T1                                                        
         ENGYGS = AER + DIS*T2*T2                                               
         IF (NDER .EQ. 1) THEN                                                  
             T3 = 2.0D0*DIS*T1                                                  
             DEGSDR(2) = DAER + T2*(DDIS*T2 - T3*(BET*DREQ - DBET*DELR))        
C                                                                               
             DEGSDR(1) = T3*BET*T2                                              
         ENDIF                                                                  
      ELSE                                                                      
C                                                                               
C   SPLINE INTERPOLATION REGION                                                 
C                                                                               
         DR1 = R1S - R(1)                                                       
         DR2 = R2S - R(2)                                                       
         T1 = DR1*DR1 + DR2*DR2                                                 
         RR = SQRT(T1)                                                          
         DRRD1 = -DR1/RR                                                        
         DRRD2 = -DR2/RR                                                        
         ANGLE = ATAN(DR1/DR2)                                                  
         DPHD1 = -DR2/T1                                                        
         DPHD2 = DR1/T1                                                         
         DIS = POTSPL(3,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)                 
         DDIS = POTSPL(4,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)                
         REQ2 = POTSPL(3,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)                
         DREQ = POTSPL(4,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)                
         BET =  POTSPL(3,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)                
         DBET = POTSPL(4,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)                
         IF(IAYES.EQ.0) THEN                                                    
            AER = VSP - DIS                                                     
            DAER = -DDIS                                                        
         ELSE                                                                   
            T1 = EXP(-BET*REQ2)                                                 
            T2 = 1.0D0 - T1                                                     
            AER = VSP - DIS*T2*T2                                               
            DAER = -T2*(DDIS*T2 + 2.0D0*DIS*T1*(BET*DREQ + REQ2*DBET))          
         END IF                                                                 
         IF (BET*(RR-REQ2) .GT.  70)                                            
     *       WRITE (NFLAG(18), 6600) R(1), R(2), R(3), RR, ANGLE,               
     *                          DIS, REQ2, BET                                  
         DELR = RR - REQ2                                                       
         T1 = EXP(BET*DELR)                                                     
         T2 = 1.0D0 - T1                                                        
         ENGYGS = AER + DIS*T2*T2                                               
         IF (NDER .EQ. 1) THEN                                                  
             T3 = 2.0D0*DIS*T1                                                  
             DVDP = DAER + T2*(DDIS*T2 - T3*(DBET*DELR - BET*DREQ))             
             T3 = T3*T2*BET                                                     
             DEGSDR(1) =  DVDP*DPHD1 - T3*DRRD1                                 
             DEGSDR(2) =  DVDP*DPHD2 - T3*DRRD2                                 
         ENDIF                                                                  
      END IF                                                                    
      ENGYGS = ENGYGS                                                           
      IF (NDER .NE. 1) RETURN                                                   
      DEGSDR(1) = DEGSDR(1)                                                     
      DEGSDR(2) = DEGSDR(2)                                                     
      DEGSDR(3) = 0.0D0                                                         
C                                                                               
6600  FORMAT(/,2X,T5,'R(1), R(2), R(3), RR, ANGLE, DIS, REQ2, BET:',            
     *       /,(2X,T5,1PE20.10,T30,1PE20.10,T55,1PE20.10))                      
6601  FORMAT(/,2X,T5,'R(1), R(2), R(3), DELR, ALPH1:',                          
     +       /,2X,T10,5(1PE13.5,1X))                                            
6602  FORMAT(/,2X,T5,'R(1), R(2), R(3), DR2, BET:',                             
     +       /,2X,T10,5(1PE13.5,1X))                                            
6603  FORMAT(/,2X,T5,'R(1), R(2), R(3), DELR, ALPH2:',                          
     +       /,2X,T10,5(1PE13.5,1X))                                            
6604  FORMAT(/,2X,T5,'R(1), R(2), R(3), DR2, BET:',                             
     +       /,2X,T10,5(1PE13.5,1X))                                            
C                                                                               
      RETURN                                                                    
      END                                                                       
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
         COMMON /AMBPCM/ D(3),BETA,REQ(3),A,GAMP(3),RBB,RHH,DE                  
         COMMON /PRECM/  DE1(3),DE2(3),DCOSDR(3),DFDR(3),IOPSW                  
         COMMON /PRE2CM/ DBCOOR(2)                                              
         COMMON /PRE22CM/ D1INF,R1INF,B1INF,D2INF,                              
     +                    R2INF,B2INF,A1INF,A2INF,                              
     +                    R1S,R2S,VSP,ALPH1,ALPH2,                              
     +                    PHI(20),BSP(20,3),CSP(22,3),                          
     +                    D1DIF,R1DIF,B1DIF,A1DIF,                              
     +                    D2DIF,R2DIF,B2DIF,A2DIF,                              
     +                    DASY,IAYES,NPHI                                       
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
C   Initialize the potential parameters for the collinear part of the           
C   potential; the energy parameters are in kcal/mol, and the lengths           
C   are in Angstroms.  These are the parameters for the A'' surface.            
C                                                                               
C   Dissociation energy in kcal/mol                                             
C                                                                               
         DATA D / 106.56D0, 109.472D0, 106.56D0/                                
C                                                                               
C   Morse beta in reciprocal Angstroms                                          
C                                                                               
         DATA BETA / 2.07942D0/                                                 
C                                                                               
C   Equilibrium bond lengths in Angstroms                                       
C                                                                               
         DATA REQ / 0.96966D0, 0.74144D0, 0.96966D0/                            
C                                                                               
C   Pauling parameters in Angstroms                                             
C                                                                               
         DATA A / 0.26D0/                                                       
C                                                                               
C   Gamma parameters for the bend correction                                    
C                                                                               
         DATA GAMP / 0.210820D0, 5.42027D0, -1.44235D0/                         
C                                                                               
C   RBB in Angstroms - this number is used in                                   
C   a scaling calculation for the bending correction                            
C                                                                               
         DATA RBB / 0.74144D0/                                                  
C                                                                               
         DATA IOPSW / 1          /                                              
         DATA D1INF /  95.3d0    /                                              
         DATA R1INF /    .76d0   /                                              
         DATA B1INF /   2.0202d0 /                                              
         DATA D2INF /  92.4d0    /                                              
         DATA R2INF /    .98d0   /                                              
         DATA B2INF /   2.41105d0/                                              
         DATA A1INF /   0.0d0    /                                              
         DATA A2INF /   0.0d0    /                                              
         DATA R1S   /  2.21208d0 /                                              
         DATA R2S   /  2.21531d0 /                                              
         DATA VSP   / 95.3d0     /                                              
         DATA IAYES /  0         /                                              
         DATA ALPH1 /  1.6d0     /                                              
         DATA ALPH2 /  1.6d0     /                                              
         DATA NPHI / 11 /                                                       
         DATA (PHI(I), I=1,11)                                                  
     +        / 0.0000000d+00,  0.1000000d+02,                                  
     +          0.2000000d+02,  0.3000000d+02,                                  
     +          0.3733000d+02,  0.4000000d+02,                                  
     +          0.5000000d+02,  0.6000000d+02,                                  
     +          0.7000000d+02,  0.8000000d+02,                                  
     +          0.9000000d+02/                                                  
         DATA (BSP(I,1), I=1,11)                                                
     +        / 0.9477920d+02,  0.9316200d+02,                                  
     +          0.9035320d+02,  0.8517960d+02,                                  
     +          0.8272720d+02,  0.8304320d+02,                                  
     +          0.8855580d+02,  0.9003100d+02,                                  
     +          0.9129800d+02,  0.9184900d+02,                                  
     +          0.9209100d+02/                                                  
         DATA (BSP(I,2), I=1,11)                                                
     +        / 0.1456500d+01,  0.1475500d+01,                                  
     +          0.1530500d+01,  0.1606400d+01,                                  
     +          0.1634900d+01,  0.1632300d+01,                                  
     +          0.1534600d+01,  0.1396600d+01,                                  
     +          0.1299800d+01,  0.1243800d+01,                                  
     +          0.1225600d+01/                                                  
         DATA (BSP(I,3), I=1,11)                                                
     +        / 0.2041900d+01,  0.1987500d+01,                                  
     +          0.1905400d+01,  0.1672700d+01,                                  
     +          0.1554000d+01,  0.1539400d+01,                                  
     +          0.1749300d+01,  0.2014500d+01,                                  
     +          0.2227100d+01,  0.2348800d+01,                                  
     +          0.2383800d+01/                                                  
         DATA (CSP(I,1), I=1,13)                                                
     +        / 0.9477920d+02,  -0.1478277d+00,                                 
     +         -0.1389235d-03,  -0.3580591d-03,                                 
     +          0.1450636d-02,   0.1148717d-02,                                 
     +         -0.9458414d-03,  -0.4983968d-02,                                 
     +          0.5921762d-02,  -0.2814975d-02,                                 
     +          0.1001140d-02,  -0.2747850d-03,                                 
     +         -0.5702495d-05/                                                  
         DATA (CSP(I,2), I=1,13)                                                
     +        / 0.1456500d+01,   0.1137130d-02,                                 
     +          0.7628700d-05,  -0.9772201d-05,                                 
     +         -0.1201120d-04,   0.2800151d-05,                                 
     +         -0.4114549d-04,   0.7866376d-04,                                 
     +         -0.1122307d-04,  -0.1936703d-04,                                 
     +          0.6791209d-05,  -0.1039780d-04,                                 
     +          0.8032967d-05/                                                  
         DATA (CSP(I,3), I=1,13)                                                
     +        / 0.2041900d+01,  -0.5887558d-02,                                 
     +          0.4475579d-05,  -0.5455348d-04,                                 
     +          0.1230139d-03,  -0.6997501d-04,                                 
     +          0.2822753d-03,  -0.4019424d-03,                                 
     +          0.1215636d-03,  -0.1548133d-04,                                 
     +          0.9961687d-05,   0.1813458d-04,                                 
     +         -0.1747243d-04/                                                  
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
C                                                                               
      FUNCTION POTSPL(ISW,NN,X,Y,C,D,XPT)                                       
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
C     PARAMETER (IPRT = 6)                                                      
      DIMENSION X(NN),Y(NN),C(NN+2),D(NN+2)                                     
C                                                                               
C  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE                    
C  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT                   
C  OR ITS DERIVATIVE AT A GIVEN POINT                                           
C                                                                               
C  CALLING SEQUENCE .......                                                     
C     ISW ... CONTROL OPTION                                                    
C         ISW=1  IF A CUBIC SPLINE IS TO BE FITTED TO THE SET OF KNOTS          
C                DEFINED BY THE ARRAYS X AND Y.  THE SPLINE COEFFICIENTS        
C                ARE STORED IN THE ARRAY C.                                     
C         ISW=2  IF THE SPLINE DEFINED BY THE COEFFICIENT ARRAY 'C' IS          
C                TO BE EVALUATED (INTERPOLATED) AT THE POINT DEFINED BY         
C                THE PARAMETER 'XPT'.                                           
C         ISW=3  AS IN ISW=2, ONLY THE DERIVATIVE/3.D0 IS ALSO CALCULATE        
C         ISW=4  THE DERIVATIVE CALCULATED BY THE LAST USE OF SPLINE WIT        
C                IS RETURNED.                                                   
C                                                                               
C     NN ... THE NUMBER OF KNOTS (DATA POINTS) TO WHICH THE SPLINE IS TO        
C            BE FITTED                                                          
C                                                                               
C     X,Y ... THE ARRAYS DEFINING THE KNOTS.  THE X-VALUES MUST BE IN           
C             INCREASING ORDER.  THE ARRAYS MUST BE DIMENSIONED AT LEAST        
C             NN.                                                               
C                                                                               
C     C ... THE ARRAY THAT CONTAINS THE CUBIC SPLINE COEFFICIENTS.              
C           MUST BE DIMENSIONED AT LEAST NN+2 .                                 
C                                                                               
C     D ... A WORK SPACE.  MUST BE DIMENSIONED AT LEAST NN+2 .                  
C                                                                               
C     XPT ... THE POINT AT WHICH THE INTERPOLATION IS DESIRED (IF ISW IS        
C              SET TO 2).  THE VALUE OF POTSPL IS SET TO THE                    
C              INTERPOLATED VALUE.                                              
C                                                                               
C   Output from this subprogram are written to UNIT NFLAG(18)                   
C                                                                               
C *****  USER NOTES  *****                                                      
C                                                                               
C     INTERPOLATION INVOLVES AT LEAST TWO STEPS .......                         
C                                                                               
C       A.  CALL POTSPL WITH THE KNOTS.  THIS SETS UP THE                       
C           COEFFICIENT ARRAY C.                                                
C           EG.  DUMY=POTSPL(1,NN,X,Y,C,D,XPT)                                  
C                                                                               
C       B.  CALL POTSPL WITH THE ARRAY C WHICH WAS DEFINED BY THE               
C           PREVIOUS CALL AND WILL BE USED TO FIND THE VALUE AT THE             
C           POINT 'XPT' .                                                       
C           EG.   VALUE=POTSPL(2,NN,X,Y,C,D,XPT)                                
C                                                                               
C     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.             
C     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.                        
C                                                                               
      N = NN                                                                    
      NP1 = N+1                                                                 
      NP2 = N+2                                                                 
      Z = XPT                                                                   
      IF (ISW .EQ. 1) THEN                                                      
         C(1) = Y(1)                                                            
         D(1) = 1.0D0                                                           
         C(NP1) = 0.0D0                                                         
         D(NP1) = 0.0D0                                                         
         C(NP2) = 0.0D0                                                         
         D(NP2) = 0.0D0                                                         
         DO 41 I = 2,N                                                          
            C(I) = Y(I)-Y(1)                                                    
41          D(I) = X(I)-X(1)                                                    
         DO 410 I = 3,NP2                                                       
            IF(D(I-1).EQ.0.0D0) THEN                                            
               WRITE(NFLAG(18),1001)                                            
               STOP 'POTSPL 1'                                                  
            END IF                                                              
            PIVOT = 1.0D0/D(I-1)                                                
            IF(I.LT.NP2) THEN                                                   
               SUPD = X(I-1)-X(I-2)                                             
               IF(SUPD.LT.0.0D0) THEN                                           
                  WRITE(NFLAG(18),1000)                                         
                  STOP 'POTSPL 2'                                               
               END IF                                                           
               SUPD = SUPD*SUPD*SUPD                                            
            ELSE                                                                
               SUPD = 1.0D0                                                     
            END IF                                                              
            DFACT = SUPD*PIVOT                                                  
            CFACT = C(I-1)*PIVOT                                                
            IF(I.LE.N) THEN                                                     
               DO 47 J = I,N                                                    
                  V = X(J)-X(I-2)                                               
                  C(J) = C(J)-D(J)*CFACT                                        
47                D(J) = V*V*V-D(J)*DFACT                                       
            END IF                                                              
            IF(I.LT.NP2) THEN                                                   
               C(NP1) = C(NP1)-D(NP1)*CFACT                                     
               D(NP1) = 1.0D0-D(NP1)*DFACT                                      
            END IF                                                              
            C(NP2) = C(NP2)-D(NP2)*CFACT                                        
            D(NP2) = X(I-2)-D(NP2)*DFACT                                        
410      CONTINUE                                                               
         DO 411 I = 1,N                                                         
            J = NP2-I                                                           
            IF(J.EQ.NP1) THEN                                                   
               V = 1.0D0                                                        
            ELSE                                                                
               V = X(J)-X(J-1)                                                  
               V = V*V*V                                                        
            END IF                                                              
            IF(D(J+1).EQ.0.0D0) THEN                                            
               WRITE(NFLAG(18),1001)                                            
               STOP 'POTSPL 3'                                                  
            END IF                                                              
            C(J+1) = C(J+1)/D(J+1)                                              
            C(J) = C(J)-C(J+1)*V                                                
411      CONTINUE                                                               
         IF(D(2).EQ.0.0D0) THEN                                                 
            WRITE(NFLAG(18),1001)                                               
            STOP 'POTSPL 4'                                                     
         END IF                                                                 
         C(2) = C(2)/D(2)                                                       
      ELSE IF (ISW .EQ. 2) THEN                                                 
         POTSPL = C(1) + C(2)*(Z-X(1))                                          
         DO 51 I = 1,N                                                          
            V = Z-X(I)                                                          
            IF(V.LE.0) RETURN                                                   
            POTSPL = POTSPL + C(I+2)*V**3                                       
51       CONTINUE                                                               
      ELSE IF (ISW .EQ. 3) THEN                                                 
         POTSPL = C(1) + C(2)*(Z-X(1))                                          
         DERIV = C(2)/3.D0                                                      
         DO 53 I = 1,N                                                          
            V = Z-X(I)                                                          
            IF(V.LE.0) RETURN                                                   
            V2 = V*V                                                            
            POTSPL = POTSPL + C(I+2)*V2*V                                       
            DERIV = DERIV + C(I+2)*V2                                           
   53    CONTINUE                                                               
      ELSE IF (ISW .EQ. 4) THEN                                                 
         POTSPL = 3.D0*DERIV                                                    
      END IF                                                                    
      RETURN                                                                    
1000  FORMAT(1X,5X,'***** ERROR IN POTSPL ... UNORDERED X-VALUES                
     1*****')                                                                   
1001  FORMAT(1X,5X,' ***** ERROR IN POTSPL ... DIVIDE FAULT *****')             
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
