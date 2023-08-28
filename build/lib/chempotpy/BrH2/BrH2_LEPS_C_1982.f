C                                                                               
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

      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          BrH2                                                       
C   Functional form: Diatomics-in-molecules plus three-center term              
C   Common name:     BrH2 DIM-3C                                                
C   Interface:       potlib2001
C   References:      D. C. Clary                                                
C                    Chem. Phys. 71, 117 (1982)                                 
C                    I. Last and M. Baer                                        
C                    in Potential Energy Surfaces and Dynamics                  
C                    edited by D. G. Truhlar (Plenum, New York, 1981)           
C                    p. 519                                                     
C
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
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
C      The classical potential energy is set equal to zero for the Br           
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Parameters:                                                                 
C      Set in data statements                                                   
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R(Br-H)                                     
C                            R(2) = R(H-H)                                      
C                            R(3) = R(Br-H)                                     
C                                                                               
C   Common Blocks (used between the calling program and this potential):        
C        passes the coordinates, ground state electronic energy, and            
C        derivatives of the ground electronic state energy with respect         
C        to the coordinates.                                                    
C      NFLAG(18-20)                                                             
C        passes the FORTRAN unit number used for output from the potential      
C      /ASYCM/ EASYAB, EASYBC, EASYAC                                           
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the     
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to    
C        Re(AB), the equilibrium bond length for the AB diatomic.               
C        In this potential the AB valley represents H infinitely far from       
C        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the terms    
C        EASYBC and EASYAC represent the energies in the H2 and the HBr         
C        valleys, respectively.                                                 
C                                                                               
C   Default Parameter Values:                                                   
C      Variable      Default value                                              
C      NFLAG(18)        6                                                       
C                                                                               
C*****                                                                          
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
         COMMON /PRECM/ H(10),DH1(10),DH2(10),DH3(10),U(4,4),                   
     +                  E(4),SCR1(4),SCR2(4),                                   
     +                  XIB,C3,DENOM,RT34,RT38                                  
         COMMON /POTCM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,         
     +                  DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P                         
C                                                                               
C      LOGICAL LCOL                                                             
C      DIMENSION H(10),DH1(10),DH2(10),DH3(10),U(4,4),E(4),SCR1(4),             
C     1          SCR2(4)                                                        
C      COMMON /POTCOM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,          
C     &                DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P                          
C                                                                               
C   Echo potential parameters to UNIT NFLAG(18)                                 
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
      EASYAB = DHX                                                              
      EASYBC = DHH                                                              
      EASYAC = DHX                                                              
C                                                                               
      XIB = 0.5D0*(XIH + XIX)                                                   
      T = XIB*RHX                                                               
      C3 = 1.D0/3.D0                                                            
      DENOM = 2.D0*RHX*(1.D0+T*(1.D0+C3*T))*EXP(-T)                             
      G = G/RHX                                                                 
      RT34 = 0.25D0*SQRT(3.D0)                                                  
      RT38 = 0.5D0*RT34                                                         
C                                                                               
 600  FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',                
     *      //,1X,T5,'BrH2 DIM-3C potential energy surface',                    
     *       //,1X,T5,'Potential energy surface parameters:',                   
     *       /,1X,T5,'HH parameters (atomic units):',                           
     *       /,1X,T10,'RHH',T16,'=',T18,1PE13.5,                                
     *            T40,'AHH',T46,'=',T48,1PE13.5,                                
     *       /,1X,T10,'DHH',T16,'=',T18,1PE13.5,                                
     *            T40,'BHH',T46,'=',T48,1PE13.5,                                
     *       /,1X,T5,'HBr parameters (atomic units):',                          
     *       /,1X,T10,'RHBr',T16,'=',T18,1PE13.5,                               
     *            T40,'AHBr',T46,'=',T48,1PE13.5,                               
     *       /,1X,T10,'DHBr',T16,'=',T18,1PE13.5,                               
     *            T40,'BHBr',T46,'=',T48,1PE13.5)                               
 610  FORMAT(/,1X,T5,'Other parameters (atomic units):',                        
     *       /,1X,T10,'ETAHH',T16,'=',T18,1PE13.5,                              
     *            T40,'ETA3S',T46,'=',T48,1PE13.5,                              
     *       /,1X,T10,'ETA1P',T16,'=',T18,1PE13.5,                              
     *            T40,'ETA3P',T46,'=',T48,1PE13.5,                              
     * /,1X,T10,'XIH',T16,'=',T18,1PE13.5,T40,'XIX',T46,'=',T48,1PE13.5,        
     * /,1X,T10,'G',T16,'=',T18,1PE13.5,T40,'ALFW',T46,'=',T48,1PE13.5,         
     * //,1X,'*****')                                                           
C                                                                               
C                                                                               
      EZERO(1)=DHH                                                              
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='D. C. Clary, Chem. Phys. 71, 117 (1982)'                         
       REF(2)='I. Last and M. Baer in'                                          
       REF(3)='Potential Energy Surfaces and Dynamics'                          
       REF(4)='D. G. Truhlar (Ed.) (Plenum, New York, 1981), p. 519'            
C                                                                               
      INDEXES(1) = 35                                                           
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
C                                                                               
C                                                                               
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
      !CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   The potential parameters must be passed through the COMMON BLOCKS           
C   PT1CM, ASYCM, PT2CM,  All information passed through the                    
C   COMMON BLOCKS PT1CM, ASYCM must be in Hartree atomic units.                 
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R(Br-H)                                     
C                            R(2) = R(H-H)                                      
C                            R(3) = R(Br-H)                                     
C                                                                               
C        The energy in the AB valley, EASYAB, is equal to the energy of the     
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to    
C        Re(AB), the equilibrium bond length for the AB diatomic.               
C                                                                               
C        In this potential the AB valley represents H infinitely far from       
C        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the terms    
C        EASYBC and EASYAC represent the energies in the H2 and the HBr         
C        valleys, respectively.                                                 
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
      LOGICAL LCOL                                                              
         COMMON /PRECM/ H(10),DH1(10),DH2(10),DH3(10),U(4,4),                   
     +                  E(4),SCR1(4),SCR2(4),                                   
     +                  XIB,C3,DENOM,RT34,RT38                                  
         COMMON /POTCM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,         
     +                  DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P                         
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C   The potential routines define R1=R(Br-H), R2=R(H-Br), and R3=R(H-H);        
C   rearrange the input coordinates to match these definitions                  
C                                                                               
      R1 = R(1)                                                                 
      R2 = R(3)                                                                 
      R3 = R(2)                                                                 
C                                                                               
      IF(R1.GT.R2+R3) R1 = R2+R3                                                
      IF(R2.GT.R1+R3) R2 = R1+R3                                                
      IF(R3.GT.R1+R2) R3 = R1+R2                                                
      R1S = R1*R1                                                               
      R2S = R2*R2                                                               
      R3S = R3*R3                                                               
      R12 = R1*R2                                                               
      R13 = R1*R3                                                               
      R23 = R2*R3                                                               
      T3 = (R1S+R2S-R3S)*0.5D0                                                  
      T2 = (R1S+R3S-R2S)*0.5D0                                                  
      T1 = (R2S+R3S-R1S)*0.5D0                                                  
      CSA = T2/R13                                                              
      IF(ABS(CSA) .GT. ONE) CSA = SIGN(ONE,CSA)                                 
      CSA2 = CSA*CSA                                                            
      SNA2 = 1.D0-CSA2                                                          
      LCOL = SNA2 .LT. EPS                                                      
      SNA = SQRT(SNA2)                                                          
      T22 = 2.D0*CSA                                                            
      T11 = 0.D0                                                                
      IF(.NOT.LCOL) T11 = 2.D0*(SNA2-CSA2)/SNA                                  
      SN2A = T22*SNA                                                            
      T = T3/(R1*R13)                                                           
      DCSA21 = T22*T                                                            
      DSN2A1 = T11*T                                                            
      T = -R2/R13                                                               
      DCSA22 = T22*T                                                            
      DSN2A2 = T11*T                                                            
      T = T1/(R13*R3)                                                           
      DCSA23 = T22*T                                                            
      DSN2A3 = T11*T                                                            
      CSB = T1/R23                                                              
      IF(ABS(CSB).GT.ONE) CSB = SIGN(ONE,CSB)                                   
      CSB2 = CSB*CSB                                                            
      SNB2 = 1.D0-CSB2                                                          
      SNB = SQRT(SNB2)                                                          
      T11 = 0.D0                                                                
      IF(.NOT.LCOL) T11 = 2.D0*(SNB2-CSB2)/SNB                                  
      T22 = 2.D0*CSB                                                            
      SN2B = T22*SNB                                                            
      T = -R1/R23                                                               
      DCSB21 = T22*T                                                            
      DSN2B1 = T11*T                                                            
      T = T3/(R2*R23)                                                           
      DCSB22 = T22*T                                                            
      DSN2B2 = T11*T                                                            
      T = T2/(R23*R3)                                                           
      DCSB23 = T22*T                                                            
      DSN2B3 = T11*T                                                            
      T = T3/R12                                                                
      CSG2 = T*T                                                                
      T = 2.D0*T                                                                
      DCSG21 = T*T2/(R1*R12)                                                    
      DCSG22 = T*T1/(R12*R2)                                                    
      DCSG23 = -T*R3/R12                                                        
C                                                                               
C  DIATOMIC CURVES                                                              
C                                                                               
C    HH                                                                         
C                                                                               
      RDIF = R3 - RHH                                                           
      EX1 = EXP(-AHH*RDIF)                                                      
      RDIF2 = RDIF*RDIF                                                         
      EX2 = EXP(-BHH*RDIF*RDIF2)                                                
      T1 = DHH*EX1*EX2                                                          
      V1HH = T1*(EX1-2.D0)                                                      
      V3HH = ETAHH*T1*(EX1+2.D0)                                                
      T1 = 2.D0*AHH*T1                                                          
      T2 = 3.D0*BHH*RDIF2                                                       
      DV1HH = T1*(1.D0-EX1) - T2*V1HH                                           
      DV3HH = -T1*ETAHH*(1.D0+EX1) - T2*V3HH                                    
C                                                                               
C    HX                                                                         
C                                                                               
      CALL VHX(R1,V1S1,V3S1,V1P1,V3P1,DV1S1,DV3S1,DV1P1,DV3P1)                  
      CALL VHX(R2,V1S2,V3S2,V1P2,V3P2,DV1S2,DV3S2,DV1P2,DV3P2)                  
C                                                                               
      S11 = V1S1 + 3.D0*V3S1                                                    
      S12 = V1S2 + 3.D0*V3S2                                                    
      S21 = 3.D0*V1S1 + V3S1                                                    
      S22 = 3.D0*V1S2 + V3S2                                                    
      S31 = V1S1 - V3S1                                                         
      S32 = V1S2 - V3S2                                                         
      P11 = V1P1 + 3.D0*V3P1                                                    
      P12 = V1P2 + 3.D0*V3P2                                                    
      P21 = 3.D0*V1P1 + V3P1                                                    
      P22 = 3.D0*V1P2 + V3P2                                                    
      P31 = V1P1 - V3P1                                                         
      P32 = V1P2 - V3P2                                                         
C                                                                               
      DS11 = DV1S1 + 3.D0*DV3S1                                                 
      DS12 = DV1S2 + 3.D0*DV3S2                                                 
      DS21 = 3.D0*DV1S1 + DV3S1                                                 
      DS22 = 3.D0*DV1S2 + DV3S2                                                 
      DS31 = DV1S1 - DV3S1                                                      
      DS32 = DV1S2 - DV3S2                                                      
      DP11 = DV1P1 + 3.D0*DV3P1                                                 
      DP12 = DV1P2 + 3.D0*DV3P2                                                 
      DP21 = 3.D0*DV1P1 + DV3P1                                                 
      DP22 = 3.D0*DV1P2 + DV3P2                                                 
      DP31 = DV1P1 - DV3P1                                                      
      DP32 = DV1P2 - DV3P2                                                      
C                                                                               
C    CONSTRUCT 4X4 HAMILTONIAN PUT IN PACKED ARRAY                              
C                                                                               
      H(1) = V1HH + 0.25D0*(S11*CSA2 + P11*SNA2 + S12*CSB2 + P12*SNB2)          
      H(3) = V1HH + 0.25D0*(S11*SNA2 + P11*CSA2 + S12*SNB2 + P12*CSB2)          
      H(6) = V3HH + 0.25D0*(S21*CSA2 + P21*SNA2 + S22*CSB2 + P22*SNB2)          
      H(10) = V3HH + 0.25D0*(S21*SNA2 + P21*CSA2 + S22*SNB2 + P22*CSB2)         
      T11 = S11 - P11                                                           
      T12 = S12 - P12                                                           
      H(2) = 0.D0                                                               
      IF(.NOT.LCOL) H(2) = 0.125D0*(T11*SN2A - T12*SN2B)                        
      H(4) = RT34*(S31*CSA2 + P31*SNA2 - S32*CSB2 - P32*SNB2)                   
      T31 = S31 - P31                                                           
      T32 = S32 - P32                                                           
      H(5) = 0.D0                                                               
      IF(.NOT.LCOL) H(5) = RT38*(T31*SN2A + T32*SN2B)                           
      H(7) = H(5)                                                               
      H(8) = RT34*(S31*SNA2 + P31*CSA2 - S32*SNB2 - P32*CSB2)                   
      T21 = S21 - P21                                                           
      T22 = S22 - P22                                                           
      H(9) = 0.D0                                                               
      IF(.NOT.LCOL) H(9) = 0.125D0*(T21*SN2A - T22*SN2B)                        
C                                                                               
C   CONSTRUCT DERIVATIVE OF HAMILTONIAN MATRIX                                  
C                                                                               
      T = T11*DCSA21 + T12*DCSB21                                               
      DH1(1) = 0.25D0*(DS11*CSA2 + DP11*SNA2 + T)                               
      DH1(2) = 0.125D0*((DS11-DP11)*SN2A + T11*DSN2A1 - T12*DSN2B1)             
      DH1(3) = 0.25D0*(DS11*SNA2 + DP11*CSA2 - T)                               
      T = T11*DCSA22 + T12*DCSB22                                               
      DH2(1) = 0.25D0*(DS12*CSB2 + DP12*SNB2 + T)                               
      DH2(2) = 0.125D0*(-(DS12-DP12)*SN2B + T11*DSN2A2 - T12*DSN2B2)            
      DH2(3) = 0.25D0*(DS12*SNB2 + DP12*CSB2 - T)                               
      T = 0.25D0*(T11*DCSA23 + T12*DCSB23)                                      
      DH3(1) = DV1HH + T                                                        
      DH3(2) = 0.125D0*(T11*DSN2A3 - T12*DSN2B3)                                
      DH3(3) = DV1HH - T                                                        
      T = T21*DCSA21 + T22*DCSB21                                               
      DH1(6) = 0.25D0*(DS21*CSA2 + DP21*SNA2 + T)                               
      DH1(9) = 0.125D0*((DS21-DP21)*SN2A + T21*DSN2A1 - T22*DSN2B1)             
      DH1(10) = 0.25D0*(DS21*SNA2 + DP21*CSA2 - T)                              
      T = T21*DCSA22 + T22*DCSB22                                               
      DH2(6) = 0.25D0*(DS22*CSB2 + DP22*SNB2 + T)                               
      DH2(9) = 0.125D0*((DP22-DS22)*SN2B + T21*DSN2A2 - T22*DSN2B2)             
      DH2(10) = 0.25D0*(DS22*SNB2 + DP22*CSB2 - T)                              
      T = 0.25D0*(T21*DCSA23 + T22*DCSB23)                                      
      DH3(6) = DV3HH + T                                                        
      DH3(9) = 0.125D0*(T21*DSN2A3 - T22*DSN2B3)                                
      DH3(10) = DV3HH - T                                                       
      T = T31*DCSA21 - T32*DCSB21                                               
      DH1(4) = RT34*(DS31*CSA2 + DP31*SNA2 + T)                                 
      DH1(5) = RT38*((DS31-DP31)*SN2A + T31*DSN2A1 + T32*DSN2B1)                
      DH1(7) = DH1(5)                                                           
      DH1(8) = RT34*(DS31*SNA2 + DP31*CSA2 - T)                                 
      T = T31*DCSA22 - T32*DCSB22                                               
      DH2(4) = RT34*(-DS32*CSB2 - DP32*SNB2 + T)                                
      DH2(5) = RT38*((DS32-DP32)*SN2B + T31*DSN2A2 + T32*DSN2B2)                
      DH2(7) = DH2(5)                                                           
      DH2(8) = RT34*(-DS32*SNB2 - DP32*CSB2 - T)                                
      T = T31*DCSA23 - T32*DCSB23                                               
      DH3(4) = RT34*T                                                           
      DH3(5) = RT38*(T31*DSN2A3 + T32*DSN2B3)                                   
      DH3(7) = DH3(5)                                                           
      DH3(8) = -DH3(4)                                                          
C                                                                               
C    DIAGONALIZE H PUT EIGENVECTORS INTO U, EIGENVALUES INTO E                  
C    USE EISPACK ROUTINE RSP                                                    
C                                                                               
      CALL RSP(4,4,10,H,E,1,U,SCR1,SCR2,IERR)                                   
      IF(IERR.NE.0) THEN                                                        
         WRITE(NFLAG(18),6000) R1,R2,R3,H                                       
         WRITE(NFLAG(18),6001) V1HH,V3HH,V1S1,V3S1,V1S2,                        
     *                  V3S2,V1P1,V3P1,V1P2,V3P2                                
         STOP 'POT 3'                                                           
      END IF                                                                    
C                                                                               
C    EVALUATE DERIVATIVES OF LOWEST ROOT                                        
C                                                                               
      D1 = 0.D0                                                                 
      D2 = 0.D0                                                                 
      D3 = 0.D0                                                                 
      DO 20 L = 1,4                                                             
      T1 = U(L,1)                                                               
      LL = (L*(L-1))/2                                                          
      DO 20 K = 1,4                                                             
      T2 = T1*U(K,1)                                                            
      KK = (K*(K-1))/2                                                          
      IF(L.GE.K) INDEX = LL + K                                                 
      IF(K.GT.L) INDEX = KK + L                                                 
      D1 = D1 + T2*DH1(INDEX)                                                   
      D2 = D2 + T2*DH2(INDEX)                                                   
      D3 = D3 + T2*DH3(INDEX)                                                   
   20 CONTINUE                                                                  
C                                                                               
C   SUPPLEMENTARY TERM                                                          
C                                                                               
      T = XIH*R3                                                                
      EX1 = EXP(-T)                                                             
      SHH = (1.D0+T*(1.D0+C3*T))*EX1                                            
      DSHH = -XIH*T*(1.D0+T)*C3*EX1                                             
      T = XIB*R1                                                                
      EX1 = EXP(-T)                                                             
      SHX1 = R1*(1.D0+T*(1.D0+C3*T))*EX1/DENOM                                  
      DSHX1 = (1.D0+T*(1.D0-C3*T*T))*EX1/DENOM                                  
      T = XIB*R2                                                                
      EX1 = EXP(-T)                                                             
      SHX2 = R2*(1.D0+T*(1.D0+C3*T))*EX1/DENOM                                  
      DSHX2 = (1.D0+T*(1.D0-C3*T*T))*EX1/DENOM                                  
      RDIF = R1 - R2                                                            
      T1 = G*EXP(-ALFW*RDIF*RDIF)                                               
      T2 = SHH*(SHX1+SHX2) + SHX1*SHX2                                          
      ENGYGS = E(1) + T2*T1*CSG2 + EZERO(1)                                     
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         T3 = 2.D0*ALFW*RDIF*CSG2                                               
         DEGSDR(1) = D1 + (DSHX1*(SHH+SHX2)*CSG2 + T2*(DCSG21-T3))*T1           
         DEGSDR(3) = D2 + (DSHX2*(SHH+SHX1)*CSG2 + T2*(DCSG22+T3))*T1           
         DEGSDR(2) = D3 + (DSHH*(SHX1+SHX2)*CSG2 + T2*DCSG23)*T1                
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
 6000 FORMAT(/,2X,'Error returned by the subprogram RSP in POT',                
     *       /,2X,'R = ',T10,3(1PE13.5,1X),                                     
     *       /,2X,'H = ',(T10,4(1PE13.5,1X)))                                   
 6001 FORMAT(2X,'V1HH,V3HH,V1S1,V3S1,V1S2,V3S2,V1P1,V3P1,V1P2,V3P2=',           
     *       /,2X,(T10,4(1PE13.5,1X)))                                          
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE VHX(RR,V1S,V3S,V1P,V3P,DV1S,DV3S,DV1P,DV3P)                    
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
      COMMON /POTCM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,            
     &                DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P                           
      RDIF = RR - RHX                                                           
      RDIF2 = RDIF*RDIF                                                         
      EX1 = EXP(-AHX*RDIF)                                                      
      EX2 = EXP(-BHX*RDIF2*RDIF)                                                
      T1 = DHX*EX1*EX2                                                          
      V1S = T1*(EX1-2.D0)                                                       
      V1 = T1*(EX1+2.D0)                                                        
      V3S = ETA3S*V1                                                            
      V1P = ETA1P*V1                                                            
      V3P = ETA3P*V1                                                            
      T1 = 2.D0*AHX*T1                                                          
      T2 = 3.D0*BHX*RDIF2                                                       
      DV1S = T1*(1.D0-EX1) - T2*V1S                                             
      DV1 = -T1*(1.D0+EX1) - T2*V1                                              
      DV3S = ETA3S*DV1                                                          
      DV1P = ETA1P*DV1                                                          
      DV3P = ETA3P*DV1                                                          
      RETURN                                                                    
      END                                                                       
C                                                                               
      BLOCK DATA PTPACM                                                         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
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
      COMMON /POTCM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,            
     &               DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P                            
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
        DATA EPS /1.D-6/, ONE /1.0D0/                                           
        DATA RHH,AHH,DHH,BHH /1.4016D0, 1.0291D0, 0.17447D0, 0.018D0/           
        DATA XIH,XIX,G,ALFW /1.D0, 1.6D0, 0.22D0, 1.0D0/                        
C                                                                               
        DATA RHX /2.673D0/                                                      
        DATA AHX /0.957D0/                                                      
        DATA DHX /0.1439D0/                                                     
        DATA BHX /0.012D0/                                                      
        DATA ETAHH /0.393764D0/                                                 
        DATA ETA3S /0.322D0/                                                    
        DATA ETA1P /0.20286D0/                                                  
        DATA ETA3P /0.1771D0/                                                   
C                                                                               
      END                                                                       
C                                                                               
      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2,IERR)                           
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
      DOUBLE PRECISION A(NV),W(N),Z(NM,N),FV1(N),FV2(N)                         
C                                                                               
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                         
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)             
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)                     
C     OF A REAL SYMMETRIC PACKED MATRIX.                                        
C                                                                               
C     ON INPUT-                                                                 
C                                                                               
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL            
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM                    
C        DIMENSION STATEMENT,                                                   
C                                                                               
C        N  IS THE ORDER OF THE MATRIX  A,                                      
C                                                                               
C        NV  IS AN INTEGER VARIABLE SET EQUAL TO THE                            
C        DIMENSION OF THE ARRAY  A  AS SPECIFIED FOR                            
C        A  IN THE CALLING PROGRAM.  NV  MUST NOT BE                            
C        LESS THAN  N*(N+1)/2,                                                  
C                                                                               
C        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC                   
C        PACKED MATRIX STORED ROW-WISE,                                         
C                                                                               
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF                      
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO                  
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.            
C                                                                               
C     ON OUTPUT-                                                                
C                                                                               
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,                        
C                                                                               
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,                      
C                                                                               
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN                    
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE                   
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,                    
C                                                                               
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                           
C                                                                               
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,                
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY                 
C                                                                               
C     ------------------------------------------------------------------        
C                                                                               
      DATA ZERO,ONE/0.0D0,1.0D0/                                                
      IF (N .LE. NM) GO TO 5                                                    
      IERR = 10 * N                                                             
      GO TO 50                                                                  
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10                                   
      IERR = 20 * N                                                             
      GO TO 50                                                                  
C                                                                               
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)                                             
      IF (MATZ .NE. 0) GO TO 20                                                 
C                                                                               
C     ********** FIND EIGENVALUES ONLY **********                               
C                                                                               
      CALL  TQLRAT(N,W,FV2,IERR)                                                
      GO TO 50                                                                  
C                                                                               
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********              
C                                                                               
   20 DO 40 I = 1, N                                                            
C                                                                               
         DO 30 J = 1, N                                                         
            Z(J,I) = ZERO                                                       
   30    CONTINUE                                                               
C                                                                               
         Z(I,I) = ONE                                                           
   40 CONTINUE                                                                  
C                                                                               
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                             
      IF (IERR .NE. 0) GO TO 50                                                 
      CALL  TRBAK3(NM,N,NV,A,N,Z)                                               
   50 RETURN                                                                    
C                                                                               
C     ********** LAST CARD OF RSP **********                                    
C                                                                               
      END                                                                       
C                                                                               
      SUBROUTINE TRED3(N,NV,A,D,E,E2)                                           
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
      DOUBLE PRECISION A(NV),D(N),E(N),E2(N)                                    
C                                                                               
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,            
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.           
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).           
C                                                                               
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS                
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX                
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.                              
C                                                                               
C     ON INPUT-                                                                 
C                                                                               
C        N IS THE ORDER OF THE MATRIX,                                          
C                                                                               
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A               
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,              
C                                                                               
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC                    
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL                   
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.                             
C                                                                               
C     ON OUTPUT-                                                                
C                                                                               
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL                            
C          TRANSFORMATIONS USED IN THE REDUCTION,                               
C                                                                               
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,            
C                                                                               
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL                 
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,              
C                                                                               
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.            
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.                
C                                                                               
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,                
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY                 
C                                                                               
C     ------------------------------------------------------------------        
C                                                                               
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********                       
C                                                                               
      DATA ZERO/0.0D0/                                                          
      DO  300 II = 1, N                                                         
         I = N + 1 - II                                                         
         L = I - 1                                                              
         IZ = (I * L) / 2                                                       
         H = ZERO                                                               
         SCALE = ZERO                                                           
         IF (L .LT. 1) GO TO 130                                                
C                                                                               
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********               
C                                                                               
         DO 120 K = 1, L                                                        
            IZ = IZ + 1                                                         
            D(K) = A(IZ)                                                        
            SCALE = SCALE + ABS(D(K))                                           
  120    CONTINUE                                                               
C                                                                               
         IF (SCALE .NE. ZERO) GO TO 140                                         
  130    E(I) = ZERO                                                            
         E2(I) = ZERO                                                           
         GO TO 290                                                              
C                                                                               
  140    DO 150 K = 1, L                                                        
            D(K) = D(K) / SCALE                                                 
            H = H + D(K) * D(K)                                                 
  150    CONTINUE                                                               
C                                                                               
         E2(I) = SCALE * SCALE * H                                              
         F = D(L)                                                               
         G = -SIGN(SQRT(H),F)                                                   
         E(I) = SCALE * G                                                       
         H = H - F * G                                                          
         D(L) = F - G                                                           
         A(IZ) = SCALE * D(L)                                                   
         IF (L .EQ. 1) GO TO 290                                                
         F = ZERO                                                               
C                                                                               
         DO 240 J = 1, L                                                        
            G = ZERO                                                            
            JK = (J * (J-1)) / 2                                                
C                                                                               
C     ********** FORM ELEMENT OF A*U **********                                 
C                                                                               
            DO 180 K = 1, L                                                     
               JK = JK + 1                                                      
               IF (K .GT. J) JK = JK + K - 2                                    
               G = G + A(JK) * D(K)                                             
 180  CONTINUE                                                                  
C                                                                               
C     ********** FORM ELEMENT OF P **********                                   
C                                                                               
            E(J) = G / H                                                        
            F = F + E(J) * D(J)                                                 
  240    CONTINUE                                                               
C                                                                               
         HH = F / (H + H)                                                       
         JK = 0                                                                 
C                                                                               
C     ********** FORM REDUCED A **********                                      
C                                                                               
         DO 260 J = 1, L                                                        
            F = D(J)                                                            
            G = E(J) - HH * F                                                   
            E(J) = G                                                            
C                                                                               
            DO 260 K = 1, J                                                     
               JK = JK + 1                                                      
               A(JK) = A(JK) - F * E(K) - G * D(K)                              
  260    CONTINUE                                                               
C                                                                               
  290    D(I) = A(IZ+1)                                                         
         A(IZ+1) = SCALE * SQRT(H)                                              
  300 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
C                                                                               
C     ********** LAST CARD OF TRED3 **********                                  
C                                                                               
      END                                                                       
C                                                                               
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                            
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
      DIMENSION D(N),E2(N)                                                      
      DOUBLE PRECISION MACHEP                                                   
C                                                                               
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,           
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                        
C                                                                               
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC                      
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                             
C                                                                               
C     ON INPUT-                                                                 
C                                                                               
C        N IS THE ORDER OF THE MATRIX,                                          
C                                                                               
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,                  
C                                                                               
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE             
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.         
C                                                                               
C      ON OUTPUT-                                                               
C                                                                               
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN                  
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND                  
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE                    
C          THE SMALLEST EIGENVALUES,                                            
C                                                                               
C        E2 HAS BEEN DESTROYED,                                                 
C                                                                               
C        IERR IS SET TO                                                         
C          ZERO       FOR NORMAL RETURN,                                        
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN                       
C                     DETERMINED AFTER 30 ITERATIONS.                           
C                                                                               
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,                
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY                 
C                                                                               
C     ------------------------------------------------------------------        
C                                                                               
C                **********                                                     
C                                                                               
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/                                      
C                                                                               
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING             
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.           
C                                                                               
      MACHEP = TWO**(-37)                                                       
C                                                                               
      IERR = 0                                                                  
      IF (N .EQ. 1) GO TO 1001                                                  
C                                                                               
      DO 100 I = 2, N                                                           
  100 E2(I-1) = E2(I)                                                           
C                                                                               
      F = ZERO                                                                  
      B = ZERO                                                                  
      E2(N) = ZERO                                                              
C                                                                               
      DO 290 L = 1, N                                                           
         J = 0                                                                  
         H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))                                 
         IF (B .GT. H) GO TO 105                                                
         B = H                                                                  
         C = B * B                                                              
C                                                                               
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********         
C                                                                               
  105    DO 110 M = L, N                                                        
            IF (E2(M) .LE. C) GO TO 120                                         
C                                                                               
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT                      
C                THROUGH THE BOTTOM OF THE LOOP **********                      
C                                                                               
  110    CONTINUE                                                               
C                                                                               
  120    IF (M .EQ. L) GO TO 210                                                
  130    IF (J .EQ. 30) GO TO 1000                                              
         J = J + 1                                                              
C                                                                               
C     ********** FORM SHIFT **********                                          
C                                                                               
         L1 = L + 1                                                             
         S = SQRT(E2(L))                                                        
         G = D(L)                                                               
         P = (D(L1) - G) / (TWO * S)                                            
         RR = SQRT(P*P+ONE)                                                     
         D(L) = S / (P + SIGN(RR,P))                                            
         H = G - D(L)                                                           
C                                                                               
         DO 140 I = L1, N                                                       
  140    D(I) = D(I) - H                                                        
C                                                                               
         F = F + H                                                              
C                                                                               
C     ********** RATIONAL QL TRANSFORMATION **********                          
C                                                                               
         G = D(M)                                                               
         IF (G .EQ. ZERO) G = B                                                 
         H = G                                                                  
         S = ZERO                                                               
         MML = M - L                                                            
C                                                                               
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********                     
C                                                                               
         DO 200 II = 1, MML                                                     
            I = M - II                                                          
            P = G * H                                                           
            RR = P + E2(I)                                                      
            E2(I+1) = S * RR                                                    
            S = E2(I) / RR                                                      
            D(I+1) = H + S * (H + D(I))                                         
            G = D(I) - E2(I) / G                                                
            IF (G .EQ. ZERO) G = B                                              
            H = G * P / RR                                                      
  200    CONTINUE                                                               
C                                                                               
         E2(L) = S * G                                                          
         D(L) = H                                                               
C                                                                               
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********         
C                                                                               
         IF (H .EQ. ZERO) GO TO 210                                             
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210                                
         E2(L) = H * E2(L)                                                      
         IF (E2(L) .NE. ZERO) GO TO 130                                         
 210  P = D(L) + F                                                              
C                                                                               
C     ********** ORDER EIGENVALUES **********                                   
C                                                                               
         IF (L .EQ. 1) GO TO 250                                                
C                                                                               
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********                       
C                                                                               
         DO 230 II = 2, L                                                       
            I = L + 2 - II                                                      
            IF (P .GE. D(I-1)) GO TO 270                                        
            D(I) = D(I-1)                                                       
  230    CONTINUE                                                               
C                                                                               
  250    I = 1                                                                  
  270    D(I) = P                                                               
  290 CONTINUE                                                                  
C                                                                               
      GO TO 1001                                                                
C                                                                               
C     ********** SET ERROR -- NO CONVERGENCE TO AN                              
C                EIGENVALUE AFTER 30 ITERATIONS **********                      
C                                                                               
 1000 IERR = L                                                                  
 1001 RETURN                                                                    
C                                                                               
C     ********** LAST CARD OF TQLRAT **********                                 
C                                                                               
      END                                                                       
C                                                                               
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                          
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
      DIMENSION D(N),E(N),Z(NM,N)                                               
      DOUBLE PRECISION MACHEP                                                   
C                                                                               
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,             
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND             
C     WILKINSON.                                                                
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).           
C                                                                               
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS                    
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.                       
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO                      
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                          
C     FULL MATRIX TO TRIDIAGONAL FORM.                                          
C                                                                               
C     ON INPUT-                                                                 
C                                                                               
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL                 
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM                  
C          DIMENSION STATEMENT,                                                 
C                                                                               
C        N IS THE ORDER OF THE MATRIX,                                          
C                                                                               
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,                  
C                                                                               
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX                
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,                       
C                                                                               
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE                   
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS              
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN                
C          THE IDENTITY MATRIX.                                                 
C                                                                               
C      ON OUTPUT-                                                               
C                                                                               
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN                  
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT                  
C          UNORDERED FOR INDICES 1,2,...,IERR-1,                                
C                                                                               
C        E HAS BEEN DESTROYED,                                                  
C                                                                               
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC                   
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,             
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED               
C          EIGENVALUES,                                                         
C                                                                               
C        IERR IS SET TO                                                         
C          ZERO       FOR NORMAL RETURN,                                        
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN                       
C                     DETERMINED AFTER 30 ITERATIONS.                           
C                                                                               
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,                
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY                 
C                                                                               
C     ------------------------------------------------------------------        
C                                                                               
C                **********                                                     
C                                                                               
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/                                      
C                                                                               
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING             
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.           
C                                                                               
      MACHEP = TWO**(-37)                                                       
C                                                                               
      IERR = 0                                                                  
      IF (N .EQ. 1) GO TO 1001                                                  
C                                                                               
      DO 100 I = 2, N                                                           
  100 E(I-1) = E(I)                                                             
C                                                                               
      F = ZERO                                                                  
      B = ZERO                                                                  
      E(N) = ZERO                                                               
C                                                                               
      DO 240 L = 1, N                                                           
         J = 0                                                                  
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))                                   
         IF (B .LT. H) B = H                                                    
C                                                                               
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********                 
C                                                                               
         DO 110 M = L, N                                                        
            IF (ABS(E(M)) .LE. B) GO TO 120                                     
C                                                                               
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT                       
C                THROUGH THE BOTTOM OF THE LOOP **********                      
C                                                                               
  110    CONTINUE                                                               
C                                                                               
  120    IF (M .EQ. L) GO TO 220                                                
  130    IF (J .EQ. 30) GO TO 1000                                              
         J = J + 1                                                              
C                                                                               
C     ********** FORM SHIFT **********                                          
C                                                                               
         L1 = L + 1                                                             
         G = D(L)                                                               
         P = (D(L1) - G) / (TWO * E(L))                                         
         RR = SQRT(P*P+ONE)                                                     
         D(L) = E(L) / (P + SIGN(RR,P))                                         
         H = G - D(L)                                                           
C                                                                               
         DO 140 I = L1, N                                                       
  140    D(I) = D(I) - H                                                        
C                                                                               
         F = F + H                                                              
C                                                                               
C     ********** QL TRANSFORMATION **********                                   
C                                                                               
         P = D(M)                                                               
         C = ONE                                                                
         S = ZERO                                                               
         MML = M - L                                                            
C                                                                               
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********                     
C                                                                               
         DO 200 II = 1, MML                                                     
            I = M - II                                                          
            G = C * E(I)                                                        
            H = C * P                                                           
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150                                
            C = E(I) / P                                                        
            RR = SQRT(C*C+ONE)                                                  
            E(I+1) = S * P * RR                                                 
            S = C / RR                                                          
            C = ONE / RR                                                        
            GO TO 160                                                           
  150       C = P / E(I)                                                        
            RR = SQRT(C*C+ONE)                                                  
            E(I+1) = S * E(I) * RR                                              
            S = ONE / RR                                                        
            C = C * S                                                           
  160       P = C * D(I) - S * G                                                
            D(I+1) = H + S * (C * G + S * D(I))                                 
C                                                                               
C     ********** FORM VECTOR **********                                         
C                                                                               
            DO 180 K = 1, N                                                     
               H = Z(K,I+1)                                                     
               Z(K,I+1) = S * Z(K,I) + C * H                                    
               Z(K,I) = C * Z(K,I) - S * H                                      
  180       CONTINUE                                                            
C                                                                               
  200    CONTINUE                                                               
C                                                                               
         E(L) = S * P                                                           
         D(L) = C * P                                                           
         IF (ABS(E(L)) .GT. B) GO TO 130                                        
  220    D(L) = D(L) + F                                                        
  240 CONTINUE                                                                  
C                                                                               
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********                  
C                                                                               
      DO 300 II = 2, N                                                          
         I = II - 1                                                             
         K = I                                                                  
         P = D(I)                                                               
C                                                                               
         DO 260 J = II, N                                                       
            IF (D(J) .GE. P) GO TO 260                                          
            K = J                                                               
            P = D(J)                                                            
  260    CONTINUE                                                               
C                                                                               
         IF (K .EQ. I) GO TO 300                                                
         D(K) = D(I)                                                            
         D(I) = P                                                               
C                                                                               
         DO 280 J = 1, N                                                        
            P = Z(J,I)                                                          
            Z(J,I) = Z(J,K)                                                     
            Z(J,K) = P                                                          
  280    CONTINUE                                                               
C                                                                               
  300 CONTINUE                                                                  
C                                                                               
      GO TO 1001                                                                
C                                                                               
C     ********** SET ERROR -- NO CONVERGENCE TO AN                              
C                EIGENVALUE AFTER 30 ITERATIONS **********                      
C                                                                               
 1000 IERR = L                                                                  
 1001 RETURN                                                                    
C                                                                               
C     ********** LAST CARD OF TQL2 **********                                   
C                                                                               
      END                                                                       
C                                                                               
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)                                          
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
      DIMENSION A(NV),Z(NM,M)                                                   
C                                                                               
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,           
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.           
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).           
C                                                                               
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC                
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING                    
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.                        
C                                                                               
C     ON INPUT-                                                                 
C                                                                               
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL                 
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM                  
C          DIMENSION STATEMENT,                                                 
C                                                                               
C        N IS THE ORDER OF THE MATRIX,                                          
C                                                                               
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A               
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,              
C                                                                               
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS            
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST                        
C          N*(N+1)/2 POSITIONS,                                                 
C                                                                               
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,                
C                                                                               
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED                     
C          IN ITS FIRST M COLUMNS.                                              
C                                                                               
C     ON OUTPUT-                                                                
C                                                                               
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS                                
C          IN ITS FIRST M COLUMNS.                                              
C                                                                               
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.                        
C                                                                               
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,                
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY                 
C                                                                               
C     ------------------------------------------------------------------        
C                                                                               
      DATA ZERO/0.0D0/                                                          
      IF (M .EQ. 0) GO TO 200                                                   
      IF (N .EQ. 1) GO TO 200                                                   
C                                                                               
      DO 140 I = 2, N                                                           
         L = I - 1                                                              
         IZ = (I * L) / 2                                                       
         IK = IZ + I                                                            
         H = A(IK)                                                              
         IF (H .EQ. ZERO) GO TO 140                                             
C                                                                               
         DO 130 J = 1, M                                                        
            S = ZERO                                                            
            IK = IZ                                                             
C                                                                               
            DO 110 K = 1, L                                                     
               IK = IK + 1                                                      
               S = S + A(IK) * Z(K,J)                                           
 110  CONTINUE                                                                  
C                                                                               
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********           
C                                                                               
            S = (S / H) / H                                                     
            IK = IZ                                                             
C                                                                               
            DO 120 K = 1, L                                                     
               IK = IK + 1                                                      
               Z(K,J) = Z(K,J) - S * A(IK)                                      
  120       CONTINUE                                                            
C                                                                               
  130    CONTINUE                                                               
C                                                                               
  140 CONTINUE                                                                  
C                                                                               
  200 RETURN                                                                    
C                                                                               
C     ********** LAST CARD OF TRBAK3 **********                                 
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
