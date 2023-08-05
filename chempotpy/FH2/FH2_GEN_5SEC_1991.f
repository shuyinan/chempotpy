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
C   Functional form: DMBE (double-many-body expansion)                          
C   Common name:     5SEC                                                       
C   References:      G. C. Lynch, R. Steckler, D.W. Schwenke, A.J.C. Varandas,  
C                    D.G. Truhlar, and B.C. Garrett,                            
C                    J. Chem. Phys. 94, 7136-7149 (1991)                        
C                                                                               
C Number of bodies: 3
C Interface: potlib2001
C Number of electronic states: 1
C Number of derivatives: 1
C
C   Protocol:                                                                   
C      PREPOT - initializes the potential's variables                           
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
C      The classical potential energy is set equal to zero for the F            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Parameters:                                                                 
C      Set in the BLOCK DATA subprogram PTPACM                                  
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R(F-H)                                      
C                            R(2) = R(H-H)                                      
C                            R(3) = R(H-F)                                      
C                                                                               
C   Common Blocks (used between the calling program and this potential):        
C        passes the coordinates, energy, and derivatives of                     
C        the energy with respect to the coordinates.                            
C        passes the control flags where                                         
C        NDER  = 0 => no derivatives are calculated                             
C        NDER  = 1 => calculate first derivatives of the energy for the         
C                     ground electronic state with respect to the coordinates   
C        NFLAG  - Control flags                                                 
C      NFLAG(18-20)                                                             
C        passes the FORTRAN unit number used for potential output               
C      /ASYCM/ EASYAB, EASYBC, EASYAC                                           
C        passes the energy in the three asymptotic valleys for an A+BC system.  
C        The energy in the AB valley, EASYAB, is equal to the energy            
C        of the C atom "infinitely" far from the AB diatomic and R(AB) set      
C        equal to Re(AB), the equilibrium bond length for the AB diatomic.      
C        In this potential the AB valley represents H infinitely far from       
C        the FH diatomic and R(FH) equal to Re(FH).                             
C        Similarly, the terms EASYBC and EASYAC represent the energies in the   
C        H2 and the other HF asymptotic valleys, respectively.                  
C                                                                               
C   Default Parameter Values:                                                   
C      Variable      Default value                                              
C      NDER             1                                                       
C      NFLAG(18)        6                                                       
C                                                                               
C*****                                                                          
C                                                                               
C**********************************************************************         
C                                                                               
C       First:  Calculate the dynamical correlation energy e(corr)              
C       Second: Calculate the Extended-Hartree-Fock type energy given           
C               by the LEPS function calibrated from SEC energies e(ehf)        
C                                                                               
C               etot    = e(ehf) + e(corr) + deh2                               
C               e(ehf)  = e(leps) + e(3c)                                       
C               e(corr) = e(corr2) + e(corr3)                                   
C                                                                               
C       Note:   The zero of energy for this surface is F infinitely far         
C               from H2 at r(H2) = re(H2) "deh2"                                
C               To get the zero of energy r1 must be at least 50 a.u.           
C                                                                               
C**********************************************************************         
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
      COMMON /MINT/   RMT(3)                                                    
      COMMON /DIST/   RM(3),R0(3)                                               
      COMMON /REF/    RO(3)                                                     
      COMMON /TPAR/   XHH(3),XHF(3),P(4)                                        
      COMMON /ATAT/   CHH(10),CHF(10)                                           
      COMMON /COEFF/  CFHH,CHHF                                                 
      COMMON /TRIAL/  ETA1(10),ETA2(10),ETA3(10),ETAP1(10),                     
     *                ETAP2(10),ETAP3(10)                                       
      COMMON /KKP/    CK1(10),CK2(10),CK3(10),CKP1(10),CKP2(10),CKP3(10)        
      COMMON /DIAT/   D(3),A1(3),A2(3),A3(3),GAMMA(3)                           
      COMMON /SWPAR/  DEC(3),RB(3),SW(3),DSW(3)                                 
      COMMON /CSI/    FCSI,FCSI1,FCSI2                                          
      COMMON /LOCAL/  C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP                            
      COMMON /DISPC/  ALPH0,ALPH1,BET0,BET1                                     
C                                                                               
      COMMON /DAMPC/  AD1,AD2,AD3,BD1,BD2,BD3                                   
      COMMON /DPC/    DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)              
      COMMON /ECORR/  TWOBC,THBC                                                
      COMMON /DER/    DCORR2(3),DCORR3(3),DLEPS(3)                              
      COMMON /ENGYGS/ ECORR2(3),ELEPS,E3C                                       
      COMMON /PAR/    BHH(5),BHF(5)                                             
      COMMON /VTANH/  RTANH(3,10),DRTANH(3,10)                                  
C                                                                               
      PARAMETER (DEH2=0.174472692D0)                                            
      PARAMETER (DEHF=0.224989333D0)                                            
C      SAVE                                                                     
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
C                                                                               
C     CONSTANTS NEEDED FOR THE DAMPING FUNCTION                                 
C     A=ALPH0/FLOAT(N)**ALPH1        N=6,8,10                                   
C     B=BET0*EXP(-BET1*FLOAT(N))     N=6,8,10                                   
C                                                                               
      AD1= ALPH0/6.0D0**ALPH1                                                   
      AD2= ALPH0/8.0D0**ALPH1                                                   
      AD3= ALPH0/10.0D0**ALPH1                                                  
      BD1= BET0*EXP(-BET1*6.0D0)                                                
      BD2= BET0*EXP(-BET1*8.0D0)                                                
      BD3= BET0*EXP(-BET1*10.0D0)                                               
C                                                                               
C   Initialize the energy in the asymptotic valleys                             
C                                                                               
      EASYAB = DEHF                                                             
      EASYBC = DEH2                                                             
      EASYAC = DEHF                                                             
C                                                                               
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',/,              
     *       /,1X,T5,'FH2 5SEC potential energy surface',                       
     *      //,1X,T10, 'XHF:',T15,F16.6,T35,F16.6,T55,F16.6)                    
601   FORMAT(1X,T10,'XHH:',T15,F16.6,T35,F16.6,T55,F16.6)                       
602   FORMAT(1X,T10,'P:',(T15,F16.6,T35,F16.6))                                 
603   FORMAT(1X,T10,'C3C, ALF:',T35,F16.6,T55,F16.6,                            
     *     /,1X,T10,'ALFP, ALFT:',T35,F16.6,T55,F16.6,                          
     *     /,1X,T10,'P3C, Q3C:',T35,F16.6,T55,F16.6,                            
     *     /,1X,T10,'BETP:',T35,F16.6,//,1X,'*****')                            
C                                                                               
      EZERO(1)=DEH2                                                             
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='G. C. Lynch, R. Steckler, D.W. Schwenke,'                        
       REF(2)='A.J.C. Varandas, D.G. Truhlar, and B.C. Garrett,'                
       REF(3)='J. Chem. Phys. 94, 7136(1991)'                                   
C                                                                               
      INDEXES(1) = 9                                                            
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
C============================================================================== 
C                                                                               
C  PART II. POT BEGINS                                                          
C                                                                               
C============================================================================== 
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   Surfaces:                                                                   
C      ground electronic state                                                  
C                                                                               
C   Zero of energy:                                                             
C      The classical potential energy is set equal to zero for the F            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R(F-H)                                      
C                            R(2) = R(H-H)                                      
C                            R(3) = R(H-F)                                      
C                                                                               
C        The energy in the AB valley, EASYAB, is equal to the energy            
C        of the C atom "infinitely" far from the AB diatomic and R(AB) set      
C        equal to Re(AB), the equilibrium bond length for the AB diatomic.      
C        In this potential the AB valley represents H infinitely far from       
C        the FH diatomic and R(FH) equal to Re(FH).                             
C        Similarly, the terms EASYBC and EASYAC represent the energies in the   
C        H2 and the other HF asymptotic valleys, respectively.                  
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
      COMMON /MINT/   RMT(3)                                                    
      COMMON /DIST/   RM(3),R0(3)                                               
      COMMON /REF/    RO(3)                                                     
      COMMON /TPAR/   XHH(3),XHF(3),P(4)                                        
      COMMON /ATAT/   CHH(10),CHF(10)                                           
      COMMON /COEFF/  CFHH,CHHF                                                 
      COMMON /TRIAL/  ETA1(10),ETA2(10),ETA3(10),ETAP1(10),                     
     *                ETAP2(10),ETAP3(10)                                       
      COMMON /KKP/    CK1(10),CK2(10),CK3(10),CKP1(10),CKP2(10),CKP3(10)        
      COMMON /DIAT/   D(3),A1(3),A2(3),A3(3),GAMMA(3)                           
      COMMON /SWPAR/  DEC(3),RB(3),SW(3),DSW(3)                                 
      COMMON /CSI/    FCSI,FCSI1,FCSI2                                          
      COMMON /LOCAL/  C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP                            
      COMMON /DISPC/  ALPH0,ALPH1,BET0,BET1                                     
C                                                                               
      COMMON /DAMPC/  AD1,AD2,AD3,BD1,BD2,BD3                                   
      COMMON /DPC/    DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)              
      COMMON /ECORR/  TWOBC,THBC                                                
      COMMON /DER/    DCORR2(3),DCORR3(3),DLEPS(3)                              
      COMMON /ENGYGS/ ECORR2(3),ELEPS,E3C                                       
      COMMON /PAR/    BHH(5),BHF(5)                                             
      COMMON /VTANH/  RTANH(3,10),DRTANH(3,10)                                  
C                                                                               
      PARAMETER (DEH2=0.174472692D0)                                            
      PARAMETER (DEHF=0.224989333D0)                                            
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C                                                                               
C   Check the value of NDER                                                     
C                                                                               
C                                                                               
C============================================================================== 
C                                                                               
C  PART IIa.  CALCULATE THE DYNAMICAL CORRELATION ENERGY AND ITS DERIVATIVES    
C                                                                               
C============================================================================== 
C                                                                               
C     FIRST: CALL THE SUBROUTINE DAMP WHICH CALCULATES THE DAMPING FUNCTION     
C            FOR THE NTH DISPERSION COEFFICIENT AND ITS DERIVATIVES             
C            N = 6, 8, 10.                                                      
C                                                                               
      CALL DAMP(R,RM)                                                           
C                                                                               
C     SECOND: CALL THE SUBROUTINE CORR2 WHICH CALCULATES THE TWO BODY           
C             CORRELATION ENERGY AND ITS DERIVATIVES                            
C                                                                               
C                                                                               
      CALL CORR2                                                                
C                                                                               
C                                                                               
C     THIRD:  CALL XTANH TO COMPUTE HYPERBOLIC TANGENT VALUES NEEDED            
C	      IN THE SUBROUTINE CORR3, THEN                                           
C	      CALL THE SUBROUTINE CORR3 WHICH CALCULATES THE THREE BODY               
C             CORRELATION ENERGY AND ITS DERIVATIVES                            
C                                                                               
      CALL XTANH(1,ETA1)                                                        
      CALL XTANH(2,ETA2)                                                        
      CALL XTANH(3,ETA3)                                                        
C                                                                               
      CALL CORR3                                                                
C                                                                               
C     FINALLY: SUM THE TWO AND THE THREE BODY CORRELATION TERMS TO FORM E(CORR) 
C              AND DE(CORR)- THE DERIVATIVE OF E(CORR)                          
C                                                                               
C                                                                               
      EC = TWOBC + THBC                                                         
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DEC1 = DCORR2(1) + DCORR3(1)                                          
          DEC2 = DCORR2(2) + DCORR3(2)                                          
          DEC3 = DCORR2(3) + DCORR3(3)                                          
      ENDIF                                                                     
C                                                                               
C	CALCULATE THE SWITCHING FUNCTION VALUES AND                                   
C       THE CORRESPONDING DERIVATIVES                                           
C                                                                               
      CALL SWITCH                                                               
C                                                                               
      CALL VLEPS                                                                
C                                                                               
       ENGYGS = EC + ELEPS + EZERO(1)                                           
C                                                                               
      IF (NDER .NE. 1) THEN                                                     
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
          RETURN                                                                
      END IF                                                                    
      DEGSDR(1) = DEC1 + DLEPS(1)                                               
      DEGSDR(2) = DEC2 + DLEPS(2)                                               
      DEGSDR(3) = DEC3 + DLEPS(3)                                               
C                                                                               
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,             
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',             
     *                'coded in this potential')                                
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
      END                                                                       
C============================================================================== 
C                                                                               
C  PART IV.  LEAST SQUARES FIT TO SEC BENDING POTENTIAL                         
C                                                                               
C============================================================================== 
      SUBROUTINE CORR2                                                          
C============================================================================== 
C                                                                               
C     CALCULATES TWO-BODY CORRELATION ENERGY.                                   
C                                                                               
C============================================================================== 
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
      COMMON /ATAT/CHH(10),CHF(10)                                              
      COMMON /DPC/DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)                  
      COMMON /ECORR/TWOBC,THBC                                                  
      COMMON /DER/DCORR2(3),DCORR3(3),DLEPS(3)                                  
      COMMON /ENGYGS/ECORR2(3),ELEPS,E3C                                        
C                                                                               
      R1P6  = R(1)**6                                                           
      R1P8  = R(1)**8                                                           
      R1P10 = R(1)**10                                                          
      R2P6  = R(2)**6                                                           
      R2P8  = R(2)**8                                                           
      R2P10 = R(2)**10                                                          
      R3P6  = R(3)**6                                                           
      R3P8  = R(3)**8                                                           
      R3P10 = R(3)**10                                                          
C                                                                               
      ECORR2(1) = -DP1(1)*CHF(6)/R1P6                                           
     1           -DP1(2)*CHF(8)/R1P8                                            
     2           -DP1(3)*CHF(10)/R1P10                                          
C                                                                               
      ECORR2(2) = -DP2(1)*CHH(6)/R2P6                                           
     1           -DP2(2)*CHH(8)/R2P8                                            
     2           -DP2(3)*CHH(10)/R2P10                                          
C                                                                               
      ECORR2(3) = -DP3(1)*CHF(6)/R3P6                                           
     1           -DP3(2)*CHF(8)/R3P8                                            
     2           -DP3(3)*CHF(10)/R3P10                                          
C                                                                               
      TWOBC = ECORR2(1)+ECORR2(2)+ECORR2(3)                                     
C                                                                               
      IF (NDER .NE. 1) RETURN                                                   
      DCORR2(1) = -CHF(6)*(DDP1(1)-6.0D0*DP1(1)/R(1))/R1P6                      
     1            -CHF(8)*(DDP1(2)-8.0D0*DP1(2)/R(1))/R1P8                      
     2            -CHF(10)*(DDP1(3)-10.0D0*DP1(3)/R(1))/R1P10                   
C                                                                               
      DCORR2(2) = -CHH(6)*(DDP2(1)-6.0D0*DP2(1)/R(2))/R2P6                      
     1            -CHH(8)*(DDP2(2)-8.0D0*DP2(2)/R(2))/R2P8                      
     2            -CHH(10)*(DDP2(3)-10.0D0*DP2(3)/R(2))/R2P10                   
C                                                                               
      DCORR2(3) = -CHF(6)*(DDP3(1)-6.0D0*DP3(1)/R(3))/R3P6                      
     1            -CHF(8)*(DDP3(2)-8.0D0*DP3(2)/R(3))/R3P8                      
     2            -CHF(10)*(DDP3(3)-10.0D0*DP3(3)/R(3))/R3P10                   
C                                                                               
      RETURN                                                                    
      END                                                                       
C============================================================================== 
C                                                                               
      SUBROUTINE DAMP(RT1,RT2)                                                  
C============================================================================== 
C                                                                               
C     CALCULATES DAMPING FUNCTION FOR THE NTH DISPERSION COEFFICIENT            
C                                                                               
C============================================================================== 
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
      DIMENSION RT1(3),RT2(3),X(3),DX(3),POL1(3),POL2(3),POL3(3)                
      DIMENSION DPOL1(3),DPOL2(3),DPOL3(3)                                      
      COMMON /DAMPC/ AD1,AD2,AD3,BD1,BD2,BD3                                    
      COMMON /DPC/ DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)                 
      COMMON /DIST/RM(3),R0(3)                                                  
C                                                                               
            X(1) = 2.0D0*RT1(1)/(RT2(1)+2.5D0*R0(1))                            
            X(2) = 2.0D0*RT1(2)/(RT2(2)+2.5D0*R0(2))                            
            X(3) = 2.0D0*RT1(3)/(RT2(3)+2.5D0*R0(3))                            
C                                                                               
          IF (NDER .EQ. 1) THEN                                                 
            DX(1) = 2.0D0/(RT2(1)+2.5D0*R0(1))                                  
            DX(2) = 2.0D0/(RT2(2)+2.5D0*R0(2))                                  
            DX(3) = 2.0D0/(RT2(3)+2.5D0*R0(3))                                  
          ENDIF                                                                 
C                                                                               
            POL1(1) = AD1*X(1)+BD1*X(1)**2                                      
            POL2(1) = AD2*X(1)+BD2*X(1)**2                                      
            POL3(1) = AD3*X(1)+BD3*X(1)**2                                      
C                                                                               
            POL1(2) = AD1*X(2)+BD1*X(2)**2                                      
            POL2(2) = AD2*X(2)+BD2*X(2)**2                                      
            POL3(2) = AD3*X(2)+BD3*X(2)**2                                      
C                                                                               
            POL1(3) = AD1*X(3)+BD1*X(3)**2                                      
            POL2(3) = AD2*X(3)+BD2*X(3)**2                                      
            POL3(3) = AD3*X(3)+BD3*X(3)**2                                      
C                                                                               
         IF (NDER .EQ. 1) THEN                                                  
            DPOL1(1) = AD1+2.0D0*BD1*X(1)                                       
            DPOL2(1) = AD2+2.0D0*BD2*X(1)                                       
            DPOL3(1) = AD3+2.0D0*BD3*X(1)                                       
C                                                                               
            DPOL1(2) = AD1+2.0D0*BD1*X(2)                                       
            DPOL2(2) = AD2+2.0D0*BD2*X(2)                                       
            DPOL3(2) = AD3+2.0D0*BD3*X(2)                                       
C                                                                               
            DPOL1(3) = AD1+2.0D0*BD1*X(3)                                       
            DPOL2(3) = AD2+2.0D0*BD2*X(3)                                       
            DPOL3(3) = AD3+2.0D0*BD3*X(3)                                       
         ENDIF                                                                  
C                                                                               
            DP1(1) = (1.0D0-EXP(-POL1(1)))**6                                   
            DP1(2) = (1.0D0-EXP(-POL2(1)))**8                                   
            DP1(3) = (1.0D0-EXP(-POL3(1)))**10                                  
C                                                                               
            DP2(1) = (1.0D0-EXP(-POL1(2)))**6                                   
            DP2(2) = (1.0D0-EXP(-POL2(2)))**8                                   
            DP2(3) = (1.0D0-EXP(-POL3(2)))**10                                  
C                                                                               
            DP3(1) = (1.0D0-EXP(-POL1(3)))**6                                   
            DP3(2) = (1.0D0-EXP(-POL2(3)))**8                                   
            DP3(3) = (1.0D0-EXP(-POL3(3)))**10                                  
C                                                                               
         IF (NDER .NE. 1) RETURN                                                
            DDP1(1) = 6.0D0*DP1(1)*EXP(-POL1(1))*DPOL1(1)*DX(1)/                
     1                (1.0D0-EXP(-POL1(1)))                                     
            DDP1(2) = 8.0D0*DP1(2)*EXP(-POL2(1))*DPOL2(1)*DX(1)/                
     1                (1.0D0-EXP(-POL2(1)))                                     
            DDP1(3) = 10.0D0*DP1(3)*EXP(-POL3(1))*DPOL3(1)*DX(1)/               
     1                (1.0D0-EXP(-POL3(1)))                                     
C                                                                               
            DDP2(1) = 6.0D0*DP2(1)*EXP(-POL1(2))*DPOL1(2)*DX(2)/                
     1                (1.0D0-EXP(-POL1(2)))                                     
            DDP2(2) = 8.0D0*DP2(2)*EXP(-POL2(2))*DPOL2(2)*DX(2)/                
     1                (1.0D0-EXP(-POL2(2)))                                     
            DDP2(3) = 10.0D0*DP2(3)*EXP(-POL3(2))*DPOL3(2)*DX(2)/               
     1                (1.0D0-EXP(-POL3(2)))                                     
C                                                                               
            DDP3(1) = 6.0D0*DP3(1)*EXP(-POL1(3))*DPOL1(3)*DX(3)/                
     1                (1.0D0-EXP(-POL1(3)))                                     
            DDP3(2) = 8.0D0*DP3(2)*EXP(-POL2(3))*DPOL2(3)*DX(3)/                
     1                (1.0D0-EXP(-POL2(3)))                                     
            DDP3(3) = 10.0D0*DP3(3)*EXP(-POL3(3))*DPOL3(3)*DX(3)/               
     1                (1.0D0-EXP(-POL3(3)))                                     
C                                                                               
      RETURN                                                                    
      END                                                                       
C============================================================================== 
C                                                                               
      SUBROUTINE CORR3                                                          
C============================================================================== 
C                                                                               
C     CALCULATES THREE-BODY CORRELATION ENERGY.                                 
C                                                                               
C============================================================================== 
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
      DIMENSION CO3(3),DC3DR1(3),DC3DR2(3),DC3DR3(3)                            
      COMMON /KKP/CK1(10),CK2(10),CK3(10),CKP1(10),CKP2(10),CKP3(10)            
      COMMON /TRIAL/ETA1(10),ETA2(10),ETA3(10),                                 
     1              ETAP1(10),ETAP2(10),ETAP3(10)                               
      COMMON /ATAT/CHH(10),CHF(10)                                              
      COMMON /DPC/DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)                  
      COMMON /ECORR/TWOBC,THBC                                                  
      COMMON /DER/DCORR2(3),DCORR3(3),DLEPS(3)                                  
      COMMON /REF/RO(3)                                                         
      COMMON /DIST/RM(3),R0(3)                                                  
      COMMON /VTANH/ RTANH(3,10),DRTANH(3,10)                                   
C                                                                               
      JJ=0                                                                      
C                                                                               
      DO 10 I=6,10,2                                                            
C                                                                               
         JJ = JJ+1                                                              
         FI = DBLE(I)                                                           
         R1PI = R(1)**I                                                         
         R2PI = R(2)**I                                                         
         R3PI = R(3)**I                                                         
C                                                                               
         G1 = 1.0D0+CK1(I)*EXP(-CKP1(I)*(R(1)-RO(1)))                           
         G2 = 1.0D0+CK2(I)*EXP(-CKP2(I)*(R(2)-RO(2)))                           
         G3 = 1.0D0+CK3(I)*EXP(-CKP3(I)*(R(3)-RO(3)))                           
C                                                                               
       IF (NDER .EQ. 1) THEN                                                    
         DG1 = -CKP1(I)*(G1-1.0D0)                                              
         DG2 = -CKP2(I)*(G2-1.0D0)                                              
         DG3 = -CKP3(I)*(G3-1.0D0)                                              
       ENDIF                                                                    
C                                                                               
         H1 = RTANH(1,I)**ETAP1(I)                                              
         H2 = RTANH(2,I)**ETAP2(I)                                              
         H3 = RTANH(3,I)**ETAP3(I)                                              
C                                                                               
       IF (NDER .EQ. 1) THEN                                                    
         DH1 = ETAP1(I)*(RTANH(1,I)**(ETAP1(I)-1))*DRTANH(1,I)                  
         DH2 = ETAP2(I)*(RTANH(2,I)**(ETAP2(I)-1))*DRTANH(2,I)                  
         DH3 = ETAP3(I)*(RTANH(3,I)**(ETAP3(I)-1))*DRTANH(3,I)                  
       ENDIF                                                                    
C                                                                               
         T1 = CHF(I)*(1.0D0-0.5D0*(G2*H3+G3*H2))                                
     1           *DP1(JJ)/R1PI                                                  
         T2 = CHH(I)*(1.0D0-0.5D0*(G3*H1+G1*H3))                                
     1           *DP2(JJ)/R2PI                                                  
         T3 = CHF(I)*(1.0D0-0.5D0*(G1*H2+G2*H1))                                
     1           *DP3(JJ)/R3PI                                                  
C                                                                               
       IF (NDER .EQ. 1) THEN                                                    
         DT1DR1 = T1*(DDP1(JJ)/DP1(JJ)-FI/R(1))                                 
         DT1DR2 = -0.5D0*CHF(I)*DP1(JJ)*(DG2*H3+G3*DH2)/R1PI                    
         DT1DR3 = -0.5D0*CHF(I)*DP1(JJ)*(G2*DH3+DG3*H2)/R1PI                    
C                                                                               
         DT2DR1 = -0.5D0*CHH(I)*DP2(JJ)*(G3*DH1+DG1*H3)/R2PI                    
         DT2DR2 = T2*DDP2(JJ)/DP2(JJ)-T2*FI/R(2)                                
         DT2DR3 = -0.5D0*CHH(I)*DP2(JJ)*(DG3*H1+G1*DH3)/R2PI                    
C                                                                               
         DT3DR1 = -0.5D0*CHF(I)*DP3(JJ)*(DG1*H2+G2*DH1)/R3PI                    
         DT3DR2 = -0.5D0*CHF(I)*DP3(JJ)*(G1*DH2+DG2*H1)/R3PI                    
         DT3DR3 = T3*(DDP3(JJ)/DP3(JJ)-FI/R(3))                                 
       ENDIF                                                                    
C                                                                               
         CO3(JJ) = T1+T2+T3                                                     
C                                                                               
       IF (NDER .EQ. 1) THEN                                                    
         DC3DR1(JJ) = DT1DR1+DT2DR1+DT3DR1                                      
         DC3DR2(JJ) = DT1DR2+DT2DR2+DT3DR2                                      
         DC3DR3(JJ) = DT1DR3+DT2DR3+DT3DR3                                      
       ENDIF                                                                    
C                                                                               
   10 CONTINUE                                                                  
C                                                                               
       THBC = CO3(1)+CO3(2)+CO3(3)                                              
C                                                                               
       IF (NDER .NE. 1) RETURN                                                  
       DCORR3(1) = DC3DR1(1)+DC3DR1(2)+DC3DR1(3)                                
       DCORR3(2) = DC3DR2(1)+DC3DR2(2)+DC3DR2(3)                                
       DCORR3(3) = DC3DR3(1)+DC3DR3(2)+DC3DR3(3)                                
C                                                                               
      RETURN                                                                    
      END                                                                       
C============================================================================== 
C                                                                               
C============================================================================== 
C                                                                               
      SUBROUTINE VLEPS                                                          
C============================================================================== 
C                                                                               
C     CALCULATES THREE-BODY EXTENDED-HARTREE-FOCK ENERGY DEFINED BY A           
C     LEPS-TYPE FUNCTION                                                        
C                                                                               
C============================================================================== 
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
      DIMENSION D(3),A1(3),A2(3),A3(3),GAMMA(3)                                 
      DIMENSION RMT(3)                                                          
      DIMENSION DCSI(3),DFCSI1(3),DFCSI2(3),EHF(3),DEHF(3)                      
      DIMENSION EHFT(3),DEHFT(3),DSUMQ(3),TI(3),DTI(3,3)                        
      DIMENSION ETRIP(3),DETRIP(3,3),DX(3),DRI(3),X(3)                          
      DIMENSION Q(3),DQ(3,3),EX(3),DEX(3,3)                                     
C                                                                               
      COMMON /TPAR/XHH(3),XHF(3),P(4)                                           
      COMMON /MINT/RMT                                                          
      COMMON /SWPAR/DEC(3),RB(3),SW(3),DSW(3)                                   
C                                                                               
      COMMON /ATAT/CHH(10),CHF(10)                                              
      COMMON /DIST/RM(3),R0(3)                                                  
      COMMON /DIAT/D,A1,A2,A3,GAMMA                                             
      COMMON /CSI/FCSI,FCSI1,FCSI2                                              
C                                                                               
      COMMON /DAMPC/AD1,AD2,AD3,BD1,BD2,BD3                                     
      COMMON /DPC/DP1(3),DP2(3),DP3(3),DDP1(3),DDP2(3),DDP3(3)                  
      COMMON /ENGYGS/ECORR2(3),ELEPS,E3C                                        
      COMMON /DER/DCORR2(3),DCORR3(3),DLEPS(3)                                  
      COMMON /PAR/BHH(5),BHF(5)                                                 
      COMMON /LOCAL/C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP                              
C     CSI IS SIMILAR TO A COSINE FUNCTION                                       
C                                                                               
      CSI = (R(1)-R(3))/R(2)                                                    
      CSIP2 = CSI**2                                                            
      CSIP3 = CSI**3                                                            
      CSIP4 = CSI**4                                                            
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DCSI(1) = 1.0D0/R(2)                                                  
          DCSI(2) = -CSI/R(2)                                                   
          DCSI(3) = -1.0D0/R(2)                                                 
      ENDIF                                                                     
C                                                                               
      FCSI1T = 1.0D0+P(1)*(1.0D0-CSIP2)+P(2)*(1.0D0-CSIP4)                      
      FCSI2T = 1.0D0+P(3)*(1.0D0-CSIP2)+P(4)*(1.0D0-CSIP4)                      
C                                                                               
      FCSI1 = FCSI1T**2                                                         
      FCSI2 = FCSI2T**2                                                         
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          TMP1 = 2.0D0*FCSI1T*(-2.0D0*P(1)*CSI-4.0D0*P(2)*CSIP3)                
          TMP2 = 2.0D0*FCSI2T*(-2.0D0*P(3)*CSI-4.0D0*P(4)*CSIP3)                
          DO 10 I=1,3                                                           
            DFCSI1(I) = TMP1*DCSI(I)                                            
            DFCSI2(I) = TMP2*DCSI(I)                                            
10        CONTINUE                                                              
      ENDIF                                                                     
C                                                                               
      R1P2 = R(1)**2                                                            
      R2P2 = R(2)**2                                                            
      R2P3 = R(2)**3                                                            
      R2P4 = R(2)**4                                                            
      R3P2 = R(3)**2                                                            
C                                                                               
C============================================================================== 
C                                                                               
C     CALCULATE THE SWITCHING FUNCTION FOR THE DIATOMIC TRIPLET STATE           
C     CURVES AND ITS DERIVATIVES : SW(I),DSW(I)                                 
C                                                                               
C     CALCULATE THE EXTENDED-HARTREE-FOCK CURVE FOR THE GROUND SINGLET          
C     STATE OF H2 OR HF AND ITS DERIVATIVES : EHF(I),DEHF(I)                    
C                                                                               
C     I => R(J) J=1,2,3                                                         
C                                                                               
C                         A.J.C. VARANDAS & J.D. SILVA,                         
C                         J. CHEM. SOC. FARADAY TRANS. II,                      
C                         82,593 (1986)                                         
C                                                                               
C     FOR AN UPDATE FOR C8 AND C10 DISPERSION COEFFICIENTS OF HF :              
C                         A.J.C. VARANDAS                                       
C                         ADV. CHEM. PHYS., 1988 IN PRESS                       
C                                                                               
C============================================================================== 
C                                                                               
      DO 20 I=1,3                                                               
            DR = R(I)-RM(I)                                                     
            DRP2 = DR**2                                                        
            EHF(I) = -D(I)*(1.0D0+A1(I)*DR+A2(I)*DRP2+A3(I)*DRP2*DR)            
     1               *EXP(-GAMMA(I)*DR)                                         
            IF (NDER .EQ. 1) DEHF(I) = -D(I)*(A1(I)+2.0D0*A2(I)*DR+             
     *          3.0D0*A3(I)*DRP2)*EXP(-GAMMA(I)*DR)-GAMMA(I)*EHF(I)             
C                                                                               
   20 CONTINUE                                                                  
C                                                                               
C     CALCULATE THE DIATOMIC TRIPLET STATE CURVES ETRIP(I)                      
C                                                                               
C     FIRST: CALCULATE ETRIP(I) FOR I=2 => FOR H2                               
C                                                                               
C            SCALE IS THE SCALING FUNCTION FOR THE DIATOMIC                     
C            TRIPLET STATE CURVES USED FOR H2 ONLY                              
C                                                                               
C            EHFT(I) IS THE POTENTIAL FOR THE LOWEST                            
C            TRIPLET STATE CURVE                                                
C                 FOR H2 : A.J.C. VARANDAS & J. BARANDO                         
C                          MOL. PHYS. 45,857 (1982)                             
C            H2 => R(2)                                                         
C                                                                               
      SCALE = 1.0D0+XHH(1)*R(2)+XHH(2)*R2P2+XHH(3)*R2P3                         
      IF (NDER .EQ. 1) DSCALE = XHH(1)+2.0D0*XHH(2)*R(2)+                       
     *                                 3.0D0*XHH(3)*R2P2                        
C                                                                               
      ARG = BHH(2)*R(2)+BHH(3)*R2P2+BHH(4)*R2P3+BHH(5)*R2P4                     
      IF (NDER .EQ. 1) DARG = BHH(2)+2.0D0*BHH(3)*R(2)+                         
     *                        3.0D0*BHH(4)*R2P2+4.0D0*BHH(5)*R2P3               
C                                                                               
      EHFT(2) = BHH(1)*EXP(-ARG)/R(2)                                           
      IF (NDER .EQ. 1) DEHFT(2) = EHFT(2)*(-DARG-1.0D0/R(2))                    
C                                                                               
      ETRIP(2) = SCALE*EHFT(2)                                                  
      IF (NDER .EQ. 1) THEN                                                     
          DETRIP(2,2) = DSCALE*EHFT(2)+SCALE*DEHFT(2)                           
          DETRIP(2,1) = 0.0D0                                                   
          DETRIP(2,3) = 0.0D0                                                   
      ENDIF                                                                     
C                                                                               
C     NEXT:                                                                     
C            CALCULATE THE EHF POTENTIAL FOR THE LOWEST SIGMA TRIPLET           
C            STATE OF HF;                                                       
C            FIT TO ab initio DATA FROM                                         
C                                      T.H. DUNNING                             
C                                      J. CHEM. PHYS.  65,3854 (1976)           
C                                                                               
C            THEN CALCULATE THE TRIPLET STATE CURVES FOR HF                     
C                                                                               
      DO 30 I=1,3,2                                                             
            ARG = BHF(2)*R(I)+BHF(3)*R(I)**2                                    
            IF (NDER .EQ. 1) DARG = BHF(2)+2.0D0*BHF(3)*R(I)                    
            EHFT(I) = BHF(1)*EXP(-ARG)/R(I)                                     
            IF (NDER .EQ. 1) DEHFT(I) = EHFT(I)*(-1.0D0/R(I)-DARG)              
C                                                                               
            X(I) = EXP(-XHF(2)*R(I))                                            
            IF (NDER .EQ. 1) DX(I) = -XHF(2)*X(I)                               
            DRI(I) = 1.0D0                                                      
            XIP2 = X(I)**2                                                      
C                                                                               
            ETRIP(I) = XHF(1)*(FCSI1*XIP2+FCSI2*XHF(3)*X(I))/R(I)               
            IF (NDER .NE. 1) GO TO 30                                           
C                                                                               
      DO 31 J=1,3                                                               
            IF (J .NE. I) THEN                                                  
                DX(J) = 0.0D0                                                   
                DRI(J) = 0.0D0                                                  
            ENDIF                                                               
            DETRIP(I,J) = XHF(1)*(DFCSI1(J)*XIP2+2.0D0*FCSI1*X(I)*DX(J)+        
     1                  DFCSI2(J)*XHF(3)*X(I)+FCSI2*XHF(3)*DX(J))/R(I)-         
     2                  ETRIP(I)*DRI(J)/R(I)                                    
   31 CONTINUE                                                                  
   30 CONTINUE                                                                  
C                                                                               
C                                                                               
C     FINALLY:                                                                  
C             COMBINE THE THREE-BODY EXTENDED-HARTREE-FOCK ENERGY BY A          
C             LEPS TYPE FUNCTION                                                
C                                                                               
      SUMQ = 0.0D0                                                              
      IF (NDER .EQ. 1) THEN                                                     
          DSUMQ(1) = 0.0D0                                                      
          DSUMQ(2) = 0.0D0                                                      
          DSUMQ(3) = 0.0D0                                                      
      ENDIF                                                                     
C                                                                               
      DO 40 I=1,3                                                               
            TI(I) = (1.0D0-SW(I))*ETRIP(I)+SW(I)*EHFT(I)                        
            Q(I) = 0.5D0*(EHF(I)+TI(I))                                         
            EX(I) = 0.5D0*(EHF(I)-TI(I))                                        
            SUMQ = SUMQ+Q(I)                                                    
            IF (NDER .NE. 1) GO TO 40                                           
      DO 41 J=1,3                                                               
            IF(J .EQ. I)THEN                                                    
               DTI(I,J) = (1.0D0-SW(I))*DETRIP(I,J)-ETRIP(I)*DSW(J)+            
     1                    EHFT(I)*DSW(J)+SW(I)*DEHFT(J)                         
               DQ(I,J) = 0.5D0*(DEHF(J)+DTI(I,J))                               
               DEX(I,J) = 0.5D0*(DEHF(J)-DTI(I,J))                              
            ELSE                                                                
               DTI(I,J) = (1.0D0-SW(I))*DETRIP(I,J)                             
               DQ(I,J) = 0.5D0*DTI(I,J)                                         
               DEX(I,J) = -DQ(I,J)                                              
            ENDIF                                                               
   41 CONTINUE                                                                  
   40 CONTINUE                                                                  
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DSUMQ(1) = DQ(1,1)+DQ(2,1)+DQ(3,1)                                    
          DSUMQ(2) = DQ(1,2)+DQ(2,2)+DQ(3,2)                                    
          DSUMQ(3) = DQ(1,3)+DQ(2,3)+DQ(3,3)                                    
      ENDIF                                                                     
C                                                                               
      FF1 = (EX(1) - EX(2)) ** 2                                                
      FF2 = (EX(2) - EX(3)) ** 2                                                
      FF3 = (EX(3) - EX(1)) ** 2                                                
C                                                                               
      EXCH = SQRT(0.5D0*(FF1+FF2+FF3))                                          
C                                                                               
      ELEPS = SUMQ - EXCH                                                       
C                                                                               
      IF (NDER .NE. 1) GO TO 700                                                
      XTERM1 = 2.0D0*EX(1)-EX(2)-EX(3)                                          
      XTERM2 = 2.0D0*EX(2)-EX(1)-EX(3)                                          
      XTERM3 = 2.0D0*EX(3)-EX(1)-EX(2)                                          
C                                                                               
      DO 50 I=1,3                                                               
            TEXCH = EXCH                                                        
            IF (TEXCH .EQ. 0.0D0)TEXCH = 1.0D0                                  
            XTERM = 0.5D0*(XTERM1*DEX(1,I)+XTERM2*DEX(2,I)+                     
     1               XTERM3*DEX(3,I))/TEXCH                                     
            DLEPS(I) = DSUMQ(I)-XTERM                                           
   50 CONTINUE                                                                  
C                                                                               
700   CONTINUE                                                                  
      CALL THCENT(DE3C1,DE3C2,DE3C3)                                            
C                                                                               
      ELEPS = ELEPS + E3C                                                       
C                                                                               
      IF (NDER .NE. 1) RETURN                                                   
      DLEPS(1) = DLEPS(1) + DE3C1                                               
      DLEPS(2) = DLEPS(2) + DE3C2                                               
      DLEPS(3) = DLEPS(3) + DE3C3                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE THCENT(DE3C1,DE3C2,DE3C3)                                      
C===============================================================================
C                                                                               
C  LOCAL THREE-BODY EHF-TYPE TERM; ADAPTED FROM BROWN ET AL, J. CHEM. PHYS.,    
C  82(1985)188 WITH MODIFICATIONS                                               
C                                                                               
C============================================================================== 
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
      LOGICAL LDERIV                                                            
      COMMON /LOCAL/C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP                              
      COMMON /ENGYGS/ECORR2(3),ELEPS,E3C                                        
      DATA LDERIV /.FALSE./                                                     
C                                                                               
      IF (NDER .EQ. 1) LDERIV = .TRUE.                                          
      OMP3C = 1.D0 - P3C                                                        
      OMQ3C = 1.D0 - Q3C                                                        
C      HA2 = 0.5D0 * A2                                                         
C      G2 = G * G                                                               
      R12 = R(1) * R(1)                                                         
      R22 = R(2) * R(2)                                                         
      R32 = R(3) * R(3)                                                         
      RSUM = R(1) + R(3)                                                        
      RDIF = R(1) - R(3)                                                        
      RDIF2 = RDIF * RDIF                                                       
      T1 = -ALFP * RSUM                                                         
      EX1 = C3C * EXP(T1 * RSUM)                                                
      T2 = -ALF * RDIF                                                          
      EX2 = OMP3C * EXP(T2 * RDIF)                                              
      T3 = -ALFT * RDIF2                                                        
      EX3 = P3C * EXP(T3 * RDIF2)                                               
      R1I = 1.D0 / R(1)                                                         
      R3I = 1.D0 / R(3)                                                         
      T4 = R1I*R3I                                                              
C                                                                               
C   COSTH IS THE ANGLE BETWEEN R1 AND R3                                        
C                                                                               
      COSTH = 0.5D0 * (R12 + R32 - R22) * T4                                    
      T5 = OMQ3C * COSTH                                                        
C                                                                               
C   SET MODIFICATION ...                                                        
C                                                                               
C      T=1.0D0                                                                  
      T = Q3C + T5 * COSTH                                                      
C                                                                               
C   E3C CONTAINS THE THREE CENTER CORRECTION TO V                               
C                                                                               
      E3C = (EX2 + EX3) * T * EX1                                               
      F3=EXP(-BETP*(R(1)+R(3)-R(2))**2)                                         
      E3C=E3C*F3                                                                
      IF (.NOT.LDERIV) GO TO 280                                                
C                                                                               
C   COMPUTE DERIVATIVE OF THE 3C TERM                                           
C   FIRST, DERIVATIVE OF COSTH                                                  
C                                                                               
         DCOS1 = R3I - COSTH*R1I                                                
         DCOS2 = - R(2)*T4                                                      
         DCOS3 = R1I - COSTH*R3I                                                
         T4 = EX2 + EX3                                                         
         T5 = T4 * T5                                                           
C                                                                               
C   NOW, DERIVATIVES OF E3C                                                     
C                                                                               
         DE3C1 = 2.D0*((T2 * EX2 + 2.D0 * T3 * RDIF * EX3 + T1 * T4)            
     *      * T + T5 * DCOS1) * EX1*F3+                                         
     1      2.0D0*E3C*(-BETP*(R(1)+R(3)-R(2)))                                  
         DE3C2 = 2.D0 * T5 * DCOS2 * EX1*F3+                                    
     1           2.0D0*E3C*(BETP*(R(1)+R(3)-R(2)))                              
         DE3C3 = 2.D0*((-T2 * EX2 - 2.D0 * T3 * RDIF * EX3 + T1 * T4)           
     *      * T + T5 * DCOS3) * EX1*F3+                                         
     1      2.0D0*E3C*(-BETP*(R(1)+R(3)-R(2)))                                  
  280 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C                                                                               
C  SUBROUTINE XTANH : THIS IS SUBROUTINE CALCULATES THE HYPERBOLIC              
C                     TANGENT VALUE AND ITS DERIVATIVES NEEDED FOR              
C                     THE THREE-BODY DYNAMICAL CORRELATION ENERGY               
C                                                                               
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C                                                                               
        SUBROUTINE XTANH(I,A)                                                   
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
        DIMENSION A(*)                                                          
        COMMON /VTANH/ RTANH(3,10),DRTANH(3,10)                                 
C                                                                               
        DO 10 J=6,10,2                                                          
              RVAL = R(I) * A(J)                                                
              IF(RVAL .GT. 0.0D0) THEN                                          
                 EX = EXP( -2.0D0 * RVAL)                                       
                 SGN = 1.0D0                                                    
              ELSE                                                              
                 EX = EXP( 2.0D0 * RVAL)                                        
                 SGN = -1.0D0                                                   
              ENDIF                                                             
              T = 1.0D0 + EX                                                    
              TH = 1.0D0/T                                                      
              IF (NDER .EQ. 1) DRTANH(I,J) = 4.0D0*A(J)*EX*TH*TH                
              RTANH(I,J) = SGN * TH * (1.0D0 - EX)                              
   10   CONTINUE                                                                
        RETURN                                                                  
        END                                                                     
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C                                                                               
C	SUBROUTINE SWITCH: CALCULATE THE SWITCHING FUNCTION AND ITS                   
C			   DERIVATIVES NEEDED FOR THE                                               
C			   DIATOMIC TRIPLET STATE CURVES                                            
C                                                                               
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C                                                                               
      SUBROUTINE SWITCH                                                         
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
      COMMON /SWPAR/ DEC(3),RB(3),SW(3),DSW(3)                                  
C                                                                               
      DO 20 I=1,3                                                               
            RMRB = R(I)-RB(I)                                                   
            RVAL = RMRB*DEC(I)                                                  
            IF (RVAL .GT. 0.0D0)THEN                                            
                EX = EXP(-2.0D0*RVAL)                                           
                SGN = 1.0D0                                                     
            ELSE                                                                
                EX = EXP(2.0D0*RVAL)                                            
                SGN = -1.0D0                                                    
            ENDIF                                                               
            T = 1.0D0 + EX                                                      
            TH = 1.0D0/T                                                        
            IF (NDER .EQ. 1) DTH = 4.0D0*DEC(I)*EX*TH*TH                        
            TH = SGN*TH*(1.0D0-EX)                                              
            SW(I) = 0.5D0*(1.0D0+TH)                                            
            IF (NDER .EQ. 1) DSW(I) = 0.5D0*DTH                                 
C                                                                               
   20 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C============================================================================== 
      BLOCK DATA PTPACM                                                         
C============================================================================== 
C     INPUT DATA FOR FH2                                                        
C                                                                               
C============================================================================== 
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
      COMMON /ATAT/   CHH(10),CHF(10)                                           
      COMMON /DIST/   RM(3),R0(3)                                               
      COMMON /REF/    RO(3)                                                     
      COMMON /DISPC/  ALPH0,ALPH1,BET0,BET1                                     
      COMMON /TRIAL/  ETA1(10),ETA2(10),ETA3(10),                               
     +                ETAP1(10),ETAP2(10),ETAP3(10)                             
      COMMON /KKP/    CK1(10),CK2(10),CK3(10),                                  
     +                CKP1(10),CKP2(10),CKP3(10)                                
      COMMON /DIAT/   D(3),A1(3),A2(3),A3(3),GAMMA(3)                           
      COMMON /MINT/   RMT(3)                                                    
      COMMON /PAR/    BHH(5),BHF(5)                                             
      COMMON /TPAR/   XHH(3),XHF(3),P(4)                                        
      COMMON /SWPAR/  DEC(3),RB(3),SW(3),DSW(3)                                 
      COMMON /LOCAL/  C3C,ALF,ALFP,ALFT,P3C,Q3C,BETP                            
C                                                                               
C   Initialize the control flags for the potential                              
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
      DATA XHF /2.633449D0, 0.417523D0, -0.090289D0/                            
      DATA XHH /-0.393828D0, 0.167810D0, -0.017081D0/                           
      DATA P  /-0.071547D0, -0.108550D0, 0.050461D0, -0.370255D0/               
      DATA C3C, ALFP /0.296668D0, 0.078606D0/                                   
C                                                                               
      DATA ALF,ALFT/0.18D0,2.14D0/                                              
      DATA P3C,Q3C,BETP/0.95D0,1.0D0,0.18D0/                                    
      DATA DEC/2.0D0,2.0D0,2.0D0/                                               
      DATA RB/6.80D0,5.00D0,6.80D0/                                             
      DATA BHH/0.448467D0,-0.056687D0,0.246831D0,-0.018419D0,0.000598D0/        
      DATA BHF/0.9149555D0,0.3222197D0,0.1273001D0,0.0D0,0.0D0/                 
      DATA CHH(6),CHH(8),CHH(10)/6.499027D0,1.243991D2,3.2858D3/                
      DATA CHF(6),CHF(8),CHF(10)/6.68799D0,1.0289812D2,2.07391451D3/            
      DATA RM(1),RM(2),RM(3)/1.7329D0,1.4010D0,1.7329D0/                        
      DATA R0(1),R0(2),R0(3)/5.9D0,6.9283D0,5.9D0/                              
      DATA RO(1),RO(2),RO(3)/1.7329D0,1.449D0,1.7329D0/                         
      DATA ALPH0,ALPH1,BET0,BET1/2.59528D1,1.1868D0,1.57381D1,9.729D-2/         
      DATA ETA1(6),ETA1(8),ETA1(10)/0.9438421608D0,0.9774440533D0,              
     1                              0.9452240321D0/                             
      DATA ETA2(6),ETA2(8),ETA2(10)/0.1085550150D1,0.1117269260D1,              
     1                              0.1138536961D1/                             
      DATA ETA3(6),ETA3(8),ETA3(10)/0.9438421608D0,0.9774440533D0,              
     1                              0.9452240321D0/                             
      DATA ETAP1(6),ETAP1(8),ETAP1(10)/6.0D0,6.0D0,6.0D0/                       
      DATA ETAP2(6),ETAP2(8),ETAP2(10)/6.0D0,6.0D0,6.0D0/                       
      DATA ETAP3(6),ETAP3(8),ETAP3(10)/6.0D0,6.0D0,6.0D0/                       
      DATA CK1(6),CK1(8),CK1(10)/-0.2610767389D-1,-0.3428701964D-1,             
     1                           -0.4858536634D-1/                              
      DATA CK2(6),CK2(8),CK2(10)/-0.9709847270D-1,-0.1248186655D 0,             
     1                           -0.1426680807D 0/                              
      DATA CK3(6),CK3(8),CK3(10)/-0.2610767389D-1,-0.3428701964D-1,             
     1                           -0.4858536634D-1/                              
      DATA CKP1(6),CKP1(8),CKP1(10)/0.9438421608D 0,0.9774440533D 0,            
     1                            0.9452240321D 0/                              
      DATA CKP2(6),CKP2(8),CKP2(10)/0.1085550150D 1,0.1117269260D 1,            
     1                            0.1138536961D 1/                              
      DATA CKP3(6),CKP3(8),CKP3(10)/0.9438421608D 0,0.9774440533D 0,            
     1                            0.9452240321D 0/                              
      DATA D(1),A1(1),A2(1),A3(1),GAMMA(1)/0.19383609D 0,2.4822787D 0,          
     1                        1.5435337D 0,0.83093855D 0,2.3992999D 0/          
      DATA D(2),A1(2),A2(2),A3(2),GAMMA(2)/0.15796326D0,2.1977034D0,            
     1                        1.2932502D0,0.64375666D0,2.1835071D0/             
      DATA D(3),A1(3),A2(3),A3(3),GAMMA(3)/0.19383609D 0,2.4822787D 0,          
     1                        1.5435337D 0,0.83093855D 0,2.3992999D 0/          
      DATA RMT/1.0D5,1.0D5,1.0D5/                                               
      END                                                                       
C============================================================================== 
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
