      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      PARAMETER (N3ATOM = 75) 
      PARAMETER (NATOM=25)
      PARAMETER (ISURF=5)
      PARAMETER (JSURF=INT(ISURF*(ISURF+1)/2))

      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
      logical, save :: first_time_data=.true.
      double precision :: tx(9), r1(3), drdx(3,9)
      double precision :: hr(nstates,nstates,3), gr(nstates,3)
      double precision :: hx(nstates,nstates,9), gx(nstates,9)
      integer :: istate, idir, i, j, k, l

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

      if (igrad.eq.0) then
      NDER=0
      else if (igrad.eq.1 .or. igrad.eq.2) then
      NDER=1
      end if

      NASURF=1
      if(first_time_data) then
      call prepot
      first_time_data=.false.
      endif

      call pot

      p(1)=ENGYGS*27.211386
      do istate=2,nstates
        i=istate-1
        p(istate)=ENGYES(i)*27.211386
      enddo
     
      ! now compute gradients
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        tx(j)=x(iatom,idir)
      enddo
      enddo
      ! input cartesian is HHH
      r1(1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)
      r1(2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)
      r1(3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)

      
      call evdrdx(tx, r1, drdx)
      hr=0.d0
      gr=0.d0
      do i=1,3
        hr(1,2,i)=PENGYIJ(i)/0.529177211
        hr(2,1,i)=-PENGYIJ(i)/0.529177211
      enddo
      do i=1,3
        gr(1,i)=DEGSDR(i)*51.422067
        gr(2,i)=DEESDR(i,1)*51.422067
      enddo

      hx=0.d0
      gx=0.d0
      do i=1,9
      do j=1,3
        hx(:,:,i)=hx(:,:,i)+hr(:,:,j)*drdx(j,i)
        gx(:,i)=gx(:,i)+gr(:,j)*drdx(j,i)
      enddo
      enddo

      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        g(:,iatom,idir)=gx(:,j)
        d(:,:,iatom,idir)=hx(:,:,j)
      enddo
      enddo


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
C   System:          H3                                                         
C   Functional Form: Double many-body expansion                                 
C   Common Name:     DMBE                                                       
C   Reference:       A. J. C. Varandas, F. B. Brown, C. A. Mead,                
C                    D. G. Truhlar, and N. C. Blais                             
C                    J. Chem. Phys. 86, 6258 (1987).                            
C   Interface:       potlib2001
C Number of bodies: 3
C Number of electronic states: 1
C Number of derivatives: 1

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
C      coordinates - bohrs                                                      
C      derivatives - hartrees/bohr                                              
C                                                                               
C   Surfaces:                                                                   
C      ground electronic state                                                  
C                                                                               
C   Zero of energy:                                                             
C      The classical potential energy is set equal to zero for the H            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Parameters:                                                                 
C      Set in the BLOCK DATA subprogram PTPACM                                  
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R( first H -- second H)                     
C                            R(2) = R(second H -- third H )                     
C                            R(3) = R( third H -- first H )                     
C                                                                               
C   Common Blocks (used between the calling program and this potential):        
C        passes the coordinates, ground state electronic energy, and            
C        derivatives of the ground electronic state energy with respect         
C        to the coordinates.                                                    
C        passes the control flags where                                         
C              = 1, first excited electronic state energy and derivatives       
C              = 2, ground and first excited electronic state energies and      
C                   derivatives                                                 
C        NDER  - not used                                                       
C        NFLAG - Control flags                                                  
C                                                                               
C              derivatives for the excited electronic states                    
C                                                                               
C                                                                               
C              derivatives for the nonadiabatic coupling between ground and     
C              excited states.                                                  
C                                                                               
C   Default Parameter Values:                                                   
C      Variable      Default value                                              
C      NFLAG(18)        6                                                       
C                                                                               
C***********************************************************************        
C     Calculates double many-body expansion of the H3 potential.                
C        ENGYGS  - energy of surface 1 (Ground State).        [ARRAY]           
C        DEGSDR  - derivatives for surface 1 (Ground State).  [ARRAY]           
C        ENGYES  - energy of surface 2 (Excited State).       [ARRAY]           
C        DEESDR  - derivatives for surface 2 (Excited State). [ARRAY]           
C        ENGYIJ  - nonadiabatic coupling energy.              [ARRAY]           
C        CAPF   - nonadiabatic coupling.  [SUBROUTINE TO CALCULATE ENGYIJ]      
C***********************************************************************        
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
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
C**************************************************                             
C                                                 *                             
C  units conversion for anomalous adiabatic       *                             
C  coupling pengyij                               *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /IJCONVERT/CNVRTQ,CNVRTS,CNVRTCHI,ICOUPLE                          
C                                                                               
C     DIMENSION DLEP(3), DLEP2(3), DV2(3), DV3(3), DVA(3), DC3(3)               
C     DIMENSION R(3)                                                            
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
C NON STANDARD COUPLING BETWEEN SURFACES                                        
C                                                                               
         ICOUPLE = 1                                                            
         NASURF(1,2) = 0                                                        
         NASURF(2,1) = 0                                                        
         CNVRTD = 1.D0                                                          
         IF(NFLAG(1).EQ.2) CNVRTD = .529177249D0                                
         CNVRTE = 1.D0                                                          
         IF(NFLAG(2).EQ.2) CNVRTE = 1000.D0                                     
         IF(NFLAG(2).EQ.3) CNVRTE = 27.211395D0                                 
         IF(NFLAG(2).EQ.4) CNVRTE = 627.5095D0                                  
         IF(NFLAG(2).EQ.5) CNVRTE = 219474.6D0                                  
         IF(NFLAG(2).EQ.6) CNVRTE = 2625.500D0                                  
         CNVRTQ = CNVRTE/CNVRTD                                                 
         CNVRTS = CNVRTE/CNVRTD                                                 
         CNVRTCHI = CNVRTE                                                      
C                                                                               
      NEXP = 4                                                                  
      DO 10 ID = 1,3                                                            
         NEXP = NEXP + 2                                                        
         DAMPA(ID) = ALPH0/DBLE(NEXP)**ALPH1                                    
         DAMPB(ID) = BET0*EXP(-BET1*DBLE(NEXP))                                 
   10 CONTINUE                                                                  
C                                                                               
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='A. J. C. Varandas, F. B. Brown, C. A. Mead,'                     
       REF(2)='D. G. Truhlar, N. C. Blais,'                                     
       REF(3)='J. Chem. Phys. 86, 6258(1987)'                                   
C                                                                               
      INDEXES(1) = 1                                                            
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
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
C   Surfaces:                                                                   
C      ground electronic state                                                  
C                                                                               
C   Zero of energy:                                                             
C      The classical potential energy is set equal to zero for the H            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R( first H -- second H)                     
C                            R(2) = R(second H -- third H )                     
C                            R(3) = R( third H -- first H )                     
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (NATOM = 25)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
C                                                                               
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                        
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
C**************************************************                             
C                                                 *                             
C  units conversion for anomalous adiabatic       *                             
C  coupling pengyij                               *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /IJCONVERT/CNVRTQ,CNVRTS,CNVRTCHI,ICOUPLE                          
C                                                                               
      DIMENSION DLEP(3),DLEP2(3),                                               
     +          DV2(3),DV3(3),DVA(3),DC3(3)                                     
C     DIMENSION R(3)                                                            
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C Set up symmetry coordinates and their derivatives                             
C                                                                               
      CALL H3COOR (R(1), R(2), R(3))                                            
C                                                                               
C Set up 2-body correlation energies, dispersion damping terms, and             
C    their derivatives                                                          
C                                                                               
      CALL H3COR2                                                               
C                                                                               
C Get Leps potentials and derivatives                                           
C                                                                               
      CALL H3LEPS (VLEPS, VLEPS2, DLEP, DLEP2)                                  
C                                                                               
C Get VA term and its derivatives                                               
C                                                                               
      CALL H3VA (R(1), R(2), R(3), VA, DVA)                                     
C                                                                               
C Get VII term and its derivatives                                              
C                                                                               
      CALL H3VII (R(1), R(2), R(3), VII, DV2)                                   
C                                                                               
C Get VIII term and its derivatives                                             
C                                                                               
      CALL H3VIII (VIII, DV3)                                                   
C                                                                               
C Get 3-body correlation and its derivative                                     
C                                                                               
      CALL H3COR3 (CE3, DC3)                                                    
C                                                                               
C 2-body correlation term                                                       
C                                                                               
      CE2 = CORRS(1) + CORRS(2) + CORRS(3)                                      
C                                                                               
C ground state                                                                  
C                                                                               
      IF (NASURF(1,1) .GT. 0) THEN                                              
          VC = VA + VII + VIII + CE2 + CE3                                      
          ENGYGS = VLEPS + VC + EZERO(1)                                        
          IF(NDER.EQ.1) THEN                                                    
             DO 20 I = 1,3                                                      
                T = DVA(I) + DV2(I) + DV3(I) + DCORRS(I) + DC3(I)               
                DEGSDR(I) =( DLEP(I) + T)                                       
20           CONTINUE                                                           
          ENDIF                                                                 
      ENDIF                                                                     


C                                                                               
C first electronic state                                                        
C                                                                               
      IF (NASURF(2,2) .GT. 0) THEN                                              
          INDES = 1                                                             
          VC = VA + VII + VIII + CE2 + CE3                                      
          ENGYES(INDES) = VLEPS2 + VC + EZERO(2)                                
          IF(NDER.EQ.1) THEN                                                    
             DO 30 I = 1,3                                                      
                T = DVA(I) + DV2(I) + DV3(I) + DCORRS(I) + DC3(I)               
                DEESDR(I,INDES) =( DLEP2(I) + T)                                
30           CONTINUE                                                           
          ENDIF                                                                 
      ENDIF                                                                     

C                                                                               
C nonadiabatic coupling in an anomalous form where pengyij is treated as:       
C -  a vector with <psi|dpsi/dQ>,<psi|dpsi/ds>,<psi|dpsi/dchi>                  
C -  already transformed to the user's units                                    
C    (thus no RTOCART, DEDCOU calls via setting is nasurf(i,j))                 
C -  having no coded derivatives (no matter what nder is set to )               
C                                                                               
      IF (ICOUPLE .GT. 0) THEN                                                  
         CALL CAPF(R(1), R(2), R(3))                                            
         PENGYIJ(1) = ENGYIJ(1) * CNVRTQ                                        
         PENGYIJ(2) = ENGYIJ(2) * CNVRTS                                        
         PENGYIJ(3) = ENGYIJ(3) * CNVRTCHI                                      
      ENDIF                                                                     
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
      BLOCK DATA PTPACM                                                         
C                                                                               
C***********************************************************************        
C     Data for the DMBE H3 surface.                                             
C***********************************************************************        
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
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
      DATA EZERO(1) /0.174474112D0/                                             
      DATA EZERO(2) /0.174474112D0/                                             
                                                                                
      DATA NASURF /1,35*0/                                                      
      DATA NDER /0/                                                             
      DATA NFLAG /1,1,15*0,6,0,0/                                               
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/3,1,1/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /3/                                                           
C                                                                               
      DATA RHOL /2.4848D0/                                                      
      DATA ALPH0, ALPH1, BET0, BET1                                             
     +     /25.9528D0,1.1868D0,15.7381D0,0.09729D0/                             
      DATA ALPHA5, BETA1, BETA2, BETA3                                          
     +     /8.2433033D-03,0.53302897D0,0.39156612D-01,0.69996945D0/             
      DATA ALPH2 /4.735364D-1/                                                  
      DATA AL0, AL1, AL2, AL3                                                   
     +     /0.45024472D+01,-0.62467617D+01,0.40966542D+01,                      
     +      0.21813012D+01/                                                     
      DATA AZ2 /4.453649D-4/                                                    
      DATA CD0, CD1 /6.333404D-03,-1.726839D-03/                                
      DATA CHH /  6.499027D0,      1.243991D+02,  3285.8D0/                     
      DATA CK0 / -1.372843D-01,   -1.638459D-01, -1.973814D-01/                 
      DATA CK1 /  1.011204D0,      9.988099D-01,  9.399411D-01/                 
      DATA HFD, HFA1, HFA2, HFA3, HFGAM                                         
     +     /0.15796326D0,2.1977034D0,1.2932502D0,0.64375666D0,                  
     +      2.1835071D0/                                                        
      DATA H2RM, H2R0, H2RMT /1.401D0, 6.928203D0, 7.82D0/                      
      DATA H2TA, H2TB1, H2TB2, H2TB3, H2TB4                                     
     +     /0.448467D0,-0.056687D0,0.246831D0,-0.018419D0,                      
     +      0.000598D0/                                                         
C      DATA PI / 3.14159265358979D0/                                            
      DATA SQRT3 /1.73205080756887D0/                                           
      DATA XPAR /-0.9286579D-02,  0.2811592D-03, -0.4665659D-05,                
     +            0.2069800D-07,  0.2903613D+02, -0.2934824D+01,                
     +            0.7181886D0,   -0.3753218D0,   -0.1114538D+01,                
     +            0.2134221D+01, -0.4164343D0,    0.2022584D0,                  
     +           -0.4662687D-01, -0.4818623D+02,  0.2988468D0/                  
      END                                                                       
C                                                                               
      SUBROUTINE H3COOR (R1, R2, R3)                                            
C                                                                               
C**********************************************************************         
C   Calculates D3H symmetry coordinates and derivatives                         
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
      COMMON/COUPCM/PHI,DPHI(3),SPHI3,DSPHI3(3)                                 
C                                                                               
      PER = R1+R2+R3                                                            
      PER2 = PER*PER                                                            
C                                                                               
C   QCOORD and its derivatives                                                  
C                                                                               
      R12 = R1*R1                                                               
      R22 = R2*R2                                                               
      R32 = R3*R3                                                               
      QCOORD = R12+R22+R32                                                      
      DQ1 = 2.0D0*R1                                                            
      DQ2 = 2.0D0*R2                                                            
      DQ3 = 2.0D0*R3                                                            
C                                                                               
C   RHO and its derivatives                                                     
C                                                                               
      RHO = SQRT(QCOORD/3.0D0)                                                  
      T = 1.0D0/(6.0D0*RHO)                                                     
      DRHO1 = T*DQ1                                                             
      DRHO2 = T*DQ2                                                             
      DRHO3 = T*DQ3                                                             
C                                                                               
C   S, CPHI3 (cos(phi3)), and their derivatives                                 
C                                                                               
      GAMMA = 2.0D0*R12-R22-R32                                                 
      GAM2 = GAMMA*GAMMA                                                        
      DGM1 = 2.0D0*DQ1                                                          
      DGM2 = -DQ2                                                               
      DGM3 = -DQ3                                                               
      BETA = SQRT3*(R22-R32)                                                    
      BET2 = BETA*BETA                                                          
      DBT1 = 0.0D0                                                              
      DBT2 = SQRT3*DQ2                                                          
      DBT3 = -SQRT3*DQ3                                                         
      T12 = BET2 + GAM2                                                         
      T1 = SQRT(T12)                                                            
      S = T1/QCOORD                                                             
      S2 = S*S                                                                  
      IF (S .EQ. 0.0D0) THEN                                                    
         DS1 = 0.0D0                                                            
         DS2 = 0.0D0                                                            
         DS3 = 0.0D0                                                            
C                                                                               
C   For S=0, CPHI3 and its derivative should not be                             
C   used anywhere but set to zero anyway.                                       
C                                                                               
         CPHI3 = 0.0D0                                                          
         DCPHI1 = 0.0D0                                                         
         DCPHI2 = 0.0D0                                                         
         DCPHI3 = 0.0D0                                                         
      ELSE                                                                      
         DS1 = S*((BETA*DBT1+GAMMA*DGM1)/T12 - DQ1/QCOORD)                      
         DS2 = S*((BETA*DBT2+GAMMA*DGM2)/T12 - DQ2/QCOORD)                      
         DS3 = S*((BETA*DBT3+GAMMA*DGM3)/T12 - DQ3/QCOORD)                      
         T2 = 1.0D0/(T1*T12)                                                    
         CPHI3 = GAMMA*(3.0D0*BET2 - GAM2)*T2                                   
         T3 = 3.0D0*BETA*(3.0D0*GAM2 - BET2)*T2/T12                             
         DCPHI1 = T3*(GAMMA*DBT1 - BETA*DGM1)                                   
         DCPHI2 = T3*(GAMMA*DBT2 - BETA*DGM2)                                   
         DCPHI3 = T3*(GAMMA*DBT3 - BETA*DGM3)                                   
      END IF                                                                    
      IF (NASURF(2,2).EQ.0.AND.NFLAG(4).EQ.0) GO TO 10                          
      ARG=GAMMA/T1                                                              
      IF(ARG.LE.-1.D0)ARG=-1.D0                                                 
      IF(ARG.GE.1.D0)ARG=1.D0                                                   
      THETA=ACOS(ARG)                                                           
      PHI=PI-THETA                                                              
      PHI3=3.0D0*PHI                                                            
      SPHI3=SIN(PHI3)                                                           
      RTAN=CPHI3/SPHI3                                                          
      DSPHI3(1)=-RTAN*DCPHI1                                                    
      DSPHI3(2)=-RTAN*DCPHI2                                                    
      DSPHI3(3)=-RTAN*DCPHI3                                                    
      DPHI(1)=-DCPHI1/(SPHI3*3.D0)                                              
      DPHI(2)=-DCPHI2/(SPHI3*3.D0)                                              
      DPHI(3)=-DCPHI3/(SPHI3*3.D0)                                              
10    CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3COR2                                                         
C                                                                               
C***********************************************************************        
C     Calculates 2-body correlation energies for singlet and triplet            
C     states of H2, the damping factors for the dispersion terms, and           
C     their 1st derivatives                                                     
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C      DIMENSION RR(3)                                                          
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COR2)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
C                                                                               
C   Loop over three coordinates                                                 
C                                                                               
      DO 10 J = 1,3                                                             
         RSAVE = R(J)                                                           
         CORRS(J) = 0.0D0                                                       
         DCORRS(J) = 0.0D0                                                      
         CORRT(J) = 0.0D0                                                       
         DCORRT(J) = 0.0D0                                                      
         T1 = 1.0D0/RSAVE                                                       
         T = T1*T1                                                              
         T2 = T*T                                                               
         T1 = T2*T1                                                             
         NEXP = 4                                                               
C                                                                               
C      Loop over terms in dispersion expansion                                  
C                                                                               
         DO 1 ID = 1,3                                                          
            NEXP = NEXP + 2                                                     
            FNEXP = DBLE(NEXP)                                                  
            T2 = T2*T                                                           
            T1 = T1*T                                                           
C                                                                               
C     singlet                                                                   
C                                                                               
C        Calculate damping function for the nth dispersion coefficient          
C        and its 1st derivative                                                 
C                                                                               
            DENOM = H2RM + 2.5D0*H2R0                                           
            X = 2.0D0*RSAVE/DENOM                                               
            TT = DAMPB(ID)*X                                                    
            TT1 = DAMPA(ID) + TT                                                
            POL = X*TT1                                                         
            DPOL = TT1 + TT                                                     
            TT = EXP(-POL)                                                      
            TT1 = 1.0D0 - TT                                                    
            TT2 = TT1**(NEXP-1)                                                 
            D = TT1*TT2                                                         
            DD = 2.0D0*FNEXP*DPOL*TT*TT2/DENOM                                  
C                                                                               
C     store damping factors (including CHH coeff and 1/R**NEXP)                 
C     FOR later use in computing the 3-body correlation energy.                 
C                                                                               
            DAMP(ID,J) = CHH(ID)*D*T2                                           
            DDAMP(ID,J) = CHH(ID)*(DD*T2 - FNEXP*D*T1)                          
            CORRS(J) = CORRS(J) - DAMP(ID,J)                                    
            DCORRS(J) = DCORRS(J) - DDAMP(ID,J)                                 
C                                                                               
C     triplet                                                                   
C        Calculate damping function for the nth dispersion                      
C        coefficient and its 1st derivative                                     
C                                                                               
            DENOM = H2RMT + 2.5D0*H2R0                                          
            X = 2.0D0*RSAVE/DENOM                                               
            TT = DAMPB(ID)*X                                                    
            TT1 = DAMPA(ID) + TT                                                
            POL = X*TT1                                                         
            DPOL = TT1 + TT                                                     
            TT = EXP(-POL)                                                      
            TT1 = 1.0D0 - TT                                                    
            TT2 = TT1**(NEXP-1)                                                 
            D = TT1*TT2                                                         
            DD = 2.0D0*FNEXP*DPOL*TT*TT2/DENOM                                  
            CORRT(J) = CORRT(J) - CHH(ID)*D*T2                                  
            DCORRT(J) = DCORRT(J) - CHH(ID)*(DD*T2 - FNEXP*D*T1)                
    1    CONTINUE                                                               
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3COR3 (CE3, DC3)                                              
C                                                                               
C***********************************************************************        
C     Calculates 3-body correlation energy and its 1st derivatives              
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
      DIMENSION DC3(3), G(3), GD(3), H(3), HD(3)                                
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
      CE3 = 0.0D0                                                               
      DC3(1) = 0.0D0                                                            
      DC3(2) = 0.0D0                                                            
      DC3(3) = 0.0D0                                                            
C                                                                               
C   Loop over terms in dispersion expansion; dispersion damping                 
C   FACTORS are computed in H3COR2 and passed through COMMON /H3CRCM/.          
C                                                                               
      DO 30 ID = 1,3                                                            
C                                                                               
C   Loop over 3 coordinates                                                     
C                                                                               
         DO 10 I = 1,3                                                          
            RR = R(I)                                                           
C                                                                               
C     Calculate G function and its 1st derivative                               
C                                                                               
            T = CK0(ID)*EXP(-CK1(ID)*(RR-H2RM))                                 
            G(I) = 1.0D0 + T                                                    
            GD(I) = -CK1(ID)*T                                                  
C                                                                               
C     Calculate H function and its 1st derivative.                              
C                                                                               
            T = CK1(ID)*RR                                                      
            SGNT = 1.0D0                                                        
            IF (T .LT. 0.0D0) THEN                                              
               T = -T                                                           
               SGNT = -1.0D0                                                    
            END IF                                                              
            T = EXP(-T)                                                         
            T2 = T*T                                                            
            T1 = 1.0D0/(1.0D0+T2)                                               
            HYSEC = 2.0D0*T*T1                                                  
            HYTAN = SGNT*(1.0D0-T2)*T1                                          
            T1 = HYTAN**5                                                       
            H(I) = HYTAN*T1                                                     
            HD(I) = 6.0D0*CK1(ID)*T1*HYSEC*HYSEC                                
   10    CONTINUE                                                               
         DO 20 I = 1,3                                                          
            I2 = MOD(I,3) + 1                                                   
            I3 = MOD(I+1,3) + 1                                                 
            T = 1.0D0 - 0.5D0*(G(I2)*H(I3) + G(I3)*H(I2))                       
            CE3 = CE3 + T*DAMP(ID,I)                                            
            DC3(I) = DC3(I) + T*DDAMP(ID,I)                                     
            DC3(I2) = DC3(I2) -                                                 
     *         0.5D0*(GD(I2)*H(I3)+G(I3)*HD(I2))*DAMP(ID,I)                     
            DC3(I3) = DC3(I3) -                                                 
     *         0.5D0*(G(I2)*HD(I3)+GD(I3)*H(I2))*DAMP(ID,I)                     
   20    CONTINUE                                                               
   30 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3LEPS (VLEPS, VLEPS2, DLEP, DLEP2)                            
C                                                                               
C***********************************************************************        
C     Calculates 3-body extended-Hartree-Fock energy defined by a               
C     LEPS-type function.                                                       
C***********************************************************************        
C                                                                               
C   VLEPS  IS LEPS LOWER SURFACE.                                               
C   VLEPS2 IS LEPS UPPER SURFACE.                                               
C   DLEP(3) ARE LEPS LOWER SURFACE DERIVATIVES                                  
C   DLEP2(3) "   "   UPPER    "       "                                         
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
      DIMENSION DLEP(3), DLEP2(3)                                               
      DIMENSION DF(3),DQ(3),XJ(3),                                              
     +          DXJ(3,3)                                                        
C                                                                               
C  Calculate switching term for the Hartree-Fock component of the               
C  diatomic triplet state function.  The notation of Thompson et al.            
C  (J. Chem. Phys., 82,5597,1985; and references therein) is used               
C  throughout.                                                                  
C                                                                               
      T1 = -AZ2*QCOORD                                                          
      IF (S .EQ. 0.0D0) THEN                                                    
         F = EXP(T1*QCOORD)                                                     
         T1 = 2.0D0*T1*F                                                        
         DF(1) = T1*DQ1                                                         
         DF(2) = T1*DQ2                                                         
         DF(3) = T1*DQ3                                                         
      ELSE                                                                      
         T2 = 1.0D0 + S*S2*CPHI3                                                
         F = EXP(T1*QCOORD*T2)                                                  
         T1 = T1*F                                                              
         TDQ = 2.0D0*T2                                                         
         T2 = S2*QCOORD                                                         
         TDS = 3.0D0*CPHI3*T2                                                   
         TDC = S*T2                                                             
         DF(1) = T1*(TDQ*DQ1 + TDS*DS1 + TDC*DCPHI1)                            
         DF(2) = T1*(TDQ*DQ2 + TDS*DS2 + TDC*DCPHI2)                            
         DF(3) = T1*(TDQ*DQ3 + TDS*DS3 + TDC*DCPHI3)                            
      END IF                                                                    
      OMF = 1.0D0 - F                                                           
C                                                                               
      QSUM = 0.0D0                                                              
      DO 10 I = 1,3                                                             
         DQ(I) = 0.0D0                                                          
   10 CONTINUE                                                                  
C                                                                               
C  Loop over 3 coordinates                                                      
C                                                                               
      DO 50 I = 1,3                                                             
         RR = R(I)                                                              
         RRI = 1.0D0/RR                                                         
C                                                                               
C     Calculate extended-Hartree-Fock curve for ground singlet state of         
C     H2 [A.J.C. Varandas & J.D. Silva, J. Chem. Soc. Faraday II                
C     (submitted)] and 1st derivative of 2-body extended-Hartree-Fock           
C     curve for the ground-singlet state of H2.                                 
C                                                                               
         DR = RR - H2RM                                                         
         T = EXP(-HFGAM*DR)                                                     
         ES = -HFD*(1.0D0 + DR*(HFA1 + DR*(HFA2 + DR*HFA3)))*T                  
         DESDR = -HFGAM*ES - HFD*(HFA1 + DR*(2.0D0*HFA2 +                       
     *            3.0D0*DR*HFA3))*T                                             
C                                                                               
C     Compute the HFACE (Hartree-Fock-approximate correlation energy)           
C     potential for the lowest triplet state of H2 [A.J.C. Varandas &           
C     J. Brandao Mol. Phys.,45,1982,857] without the 2-body correlation         
C     energy.                                                                   
C                                                                               
         T = RR*(H2TB1 + RR*(H2TB2 + RR*(H2TB3 + RR*H2TB4)))                    
         AT = H2TA*EXP(-T)*RRI                                                  
         T = H2TB1 + RR*(2.0D0*H2TB2 + RR*(3.0D0*H2TB3 +                        
     *      RR*4.0D0*H2TB4))                                                    
         DAT = -AT*(T + RRI)                                                    
C                                                                               
C     Add in triplet and subtract singlet 2-body correlation terms              
C                                                                               
         AT = AT + CORRT(I) - CORRS(I)                                          
         DAT = DAT + DCORRT(I) - DCORRS(I)                                      
C                                                                               
C     Calculate effective diatomic triplet state curve.                         
C                                                                               
         T = EXP(-AL3*RR)                                                       
         WE= (AL0 + RR*(AL1 + RR*AL2))*T                                        
         DWE = (AL1 + 2.0D0*AL2*RR)*T - AL3*WE                                  
C                                                                               
C     Triplet energy                                                            
C                                                                               
         ET = F*WE + OMF*AT                                                     
C                                                                               
C     Contribution to coulomb term                                              
C                                                                               
         QSUM = QSUM + 0.5D0*(ES + ET)                                          
C                                                                               
C     Contribution to exchange term                                             
C                                                                               
         XJ(I) = 0.5D0*(ES - ET)                                                
C                                                                               
C     Derivatives of coulomb and exchange parts                                 
C                                                                               
         FTERM = WE - AT                                                        
         DO 20 J = 1,3                                                          
            T = 0.5D0*FTERM*DF(J)                                               
            DQ(J) = DQ(J) + T                                                   
            DXJ(I,J) = -T                                                       
   20    CONTINUE                                                               
         T = F*DWE + OMF*DAT                                                    
         DQ(I) = DQ(I) + 0.5D0*(DESDR + T)                                      
         DXJ(I,I) = DXJ(I,I) + 0.5D0*(DESDR - T)                                
   50 CONTINUE                                                                  
C                                                                               
      F1 = (XJ(1) - XJ(2))**2                                                   
      F2 = (XJ(2) - XJ(3))**2                                                   
      F3 = (XJ(3) - XJ(1))**2                                                   
      EXCH = SQRT(0.5D0*(F1 + F2 + F3))                                         
      VLEPS = QSUM - EXCH                                                       
      VLEPS2 = QSUM + EXCH                                                      
      XTERM1 = 2.0D0*XJ(1) - XJ(2) - XJ(3)                                      
      XTERM2 = 2.0D0*XJ(2) - XJ(1) - XJ(3)                                      
      XTERM3 = 2.0D0*XJ(3) - XJ(1) - XJ(2)                                      
      DO 80 I = 1,3                                                             
         TEXCH = EXCH                                                           
         IF(TEXCH.LE.0.0D0)TEXCH = 1.0D0                                        
         XTERM = 0.5D0/TEXCH*(XTERM1*DXJ(1,I)  + XTERM2*DXJ(2,I)                
     *      + XTERM3*DXJ(3,I))                                                  
         DLEP(I) = DQ(I) - XTERM                                                
         DLEP2(I) = DQ(I) + XTERM                                               
   80 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3VA (R1, R2, R3, VA, DVA)                                     
C                                                                               
C***********************************************************************        
C     Calculates Va correction energy and its 1st derivative.                   
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
      DIMENSION DVA(3)                                                          
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      T = ALPHA5*PER2                                                           
      EXPVA = EXP(-T*PER)                                                       
      IF (EXPVA .EQ. 0.0D0) THEN                                                
         VA = 0.0D0                                                             
         DVA(1) = 0.0D0                                                         
         DVA(2) = 0.0D0                                                         
         DVA(3) = 0.0D0                                                         
      ELSE                                                                      
         VA = 0.0D0                                                             
         DV = 0.0D0                                                             
         T3 = (R1-R2)*(R2-R3)*(R3-R1)                                           
         T1 = T3*T3                                                             
         T2 = 1.0D0                                                             
         T4 = 2.0D0                                                             
         DO 1 J=1,4                                                             
            T2 = T2*T1                                                          
            VA = VA + XPAR(J)*T2                                                
            DV = DV + T4*XPAR(J)*T3                                             
            T3 = T3*T1                                                          
            T4 = T4 + 2.0D0                                                     
    1    CONTINUE                                                               
         VA = VA*EXPVA                                                          
         T1 = DV*EXPVA                                                          
         T2 = 3.0D0*T*VA                                                        
         DVA(1) = T1*(R2-R3)*(R3+R2-2.0D0*R1) - T2                              
         DVA(2) = T1*(R3-R1)*(R1+R3-2.0D0*R2) - T2                              
         DVA(3) = T1*(R1-R2)*(R1+R2-2.0D0*R3) - T2                              
      END IF                                                                    
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3VII (R1, R2, R3, E, DV2)                                     
C                                                                               
C***********************************************************************        
C     Calculates VII correction energy and its 1st derivative.                  
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
      DIMENSION DV2(3)                                                          
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
C   Compute B1 function.                                                        
C                                                                               
      R1I = 1.0D0/R1                                                            
      R2I = 1.0D0/R2                                                            
      R3I = 1.0D0/R3                                                            
      COS1 = 0.5D0*R2I*R3I*(R12-R22-R32)                                        
      COS2 = 0.5D0*R1I*R3I*(R22-R12-R32)                                        
      COS3 = 0.5D0*R1I*R2I*(R32-R12-R22)                                        
      WB = 1.0D0+COS1+COS2+COS3                                                 
      WB2 = WB*WB                                                               
C                                                                               
C   WB derivatives                                                              
C                                                                               
      WB1P = (R1*R3I - 1.0D0)*R2I - R3I - (COS2+COS3)*R1I                       
      WB2P = (R2*R3I - 1.0D0)*R1I - R3I - (COS1+COS3)*R2I                       
      WB3P = (R3*R2I - 1.0D0)*R1I - R2I - (COS1+COS2)*R3I                       
C                                                                               
C   EB1 term                                                                    
C                                                                               
      EXP1 = EXP(-BETA1*PER)                                                    
      EXP3 = EXP(-BETA3*PER)                                                    
      EB1T = (XPAR(5) + XPAR(6)*PER)*EXP1                                       
      EB3T = (XPAR(14) + XPAR(15)*PER2)*EXP3                                    
      EB1 = WB*(EB1T + EB3T)                                                    
C                                                                               
C   EB1 derivatives                                                             
C                                                                               
      EB1PR = WB*(-BETA1*EB1T - BETA3*EB3T + XPAR(6)*EXP1                       
     *   + 2.0D0*PER*XPAR(15)*EXP3)                                             
      EB1PWB = EB1T + EB3T                                                      
      EB1P1 = EB1PWB*WB1P + EB1PR                                               
      EB1P2 = EB1PWB*WB2P + EB1PR                                               
      EB1P3 = EB1PWB*WB3P + EB1PR                                               
C                                                                               
C   EB2 term                                                                    
C                                                                               
      T1 = BETA2*PER                                                            
      EXP2 = EXP(-T1*PER)                                                       
      EB2 = WB2*(XPAR(7) + WB*(XPAR(8) + WB*XPAR(9)))*EXP2                      
C                                                                               
C   EB2 derivatives                                                             
C                                                                               
      EB2PWB = WB*(2.0D0*XPAR(7) + WB*(3.0D0*XPAR(8) + WB*4.0D0*XPAR(9))        
     V)*EXP2                                                                    
      EB2PR = -2.0D0*T1*EB2                                                     
      EB2P1 = EB2PWB*WB1P + EB2PR                                               
      EB2P2 = EB2PWB*WB2P + EB2PR                                               
      EB2P3 = EB2PWB*WB3P + EB2PR                                               
C                                                                               
C   EB4 term                                                                    
C      EB4A                                                                     
C                                                                               
      T2 = XPAR(10)*EXP1                                                        
      T3 = WB*XPAR(11)*EXP2                                                     
      EB4A = WB*(T2 + T3)                                                       
C                                                                               
C      EB4A derivatives                                                         
C                                                                               
      EB4APW = T2 + 2.0D0*T3                                                    
      EB4APR = -WB*(BETA1*T2 + 2.0D0*T1*T3)                                     
      EB4AP1 = EB4APW*WB1P + EB4APR                                             
      EB4AP2 = EB4APW*WB2P + EB4APR                                             
      EB4AP3 = EB4APW*WB3P + EB4APR                                             
C                                                                               
C      EB4B                                                                     
C                                                                               
      T2 = XPAR(12)*EXP1                                                        
      T3 = XPAR(13)*EXP2                                                        
      EB4B = WB*(T2 + T3)                                                       
C                                                                               
C      EB4B derivatives                                                         
C                                                                               
      EB4BPW = T2 + T3                                                          
      EB4BPR = -WB*(BETA1*T2 + 2.0D0*T1*T3)                                     
      EB4BP1 = EB4BPW*WB1P + EB4BPR                                             
      EB4BP2 = EB4BPW*WB2P + EB4BPR                                             
      EB4BP3 = EB4BPW*WB3P + EB4BPR                                             
      DR12 = R1 - R2                                                            
      DR23 = R2 - R3                                                            
      DR31 = R3 - R1                                                            
      EQ = DR12*DR12 + DR23*DR23 + DR31*DR31                                    
      EQP1 = 2.0D0*(DR12 - DR31)                                                
      EQP2 = 2.0D0*(-DR12 + DR23)                                               
      EQP3 = 2.0D0*(-DR23 + DR31)                                               
      RI = R1I + R2I + R3I                                                      
      EB4 = EB4A*RI + EB4B*EQ                                                   
C                                                                               
C   EB4 derivatives                                                             
C                                                                               
      EB4P1 = EB4AP1*RI - EB4A/R12 + EB4BP1*EQ + EB4B*EQP1                      
      EB4P2 = EB4AP2*RI - EB4A/R22 + EB4BP2*EQ + EB4B*EQP2                      
      EB4P3 = EB4AP3*RI - EB4A/R32 + EB4BP3*EQ + EB4B*EQP3                      
      E = EB1 + EB2 + EB4                                                       
      DV2(1) = EB1P1 + EB2P1 + EB4P1                                            
      DV2(2) = EB1P2 + EB2P2 + EB4P2                                            
      DV2(3) = EB1P3 + EB2P3 + EB4P3                                            
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3VIII (VIII, DV3)                                             
C                                                                               
C**********************************************************************         
C     Calculates VIII correction energy and it 1st derivatives                  
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
      DIMENSION DV3(3)                                                          
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      IF (S .EQ. 0.0D0) THEN                                                    
         VIII = 0.0D0                                                           
         DV3(1)=0.0D0                                                           
         DV3(2)=0.0D0                                                           
         DV3(3)=0.0D0                                                           
      ELSE                                                                      
         T1 = RHO - RHOL                                                        
         T5 = ALPH2*T1                                                          
         EXPV3 = EXP(-T5*T1)                                                    
         IF (EXPV3 .EQ. 0.0D0) THEN                                             
            VIII = 0.0D0                                                        
            DV3(1) = 0.0D0                                                      
            DV3(2) = 0.0D0                                                      
            DV3(3) = 0.0D0                                                      
         ELSE                                                                   
            T1 = S2*CPHI3                                                       
            T2 = 1.0D0 + S*T1                                                   
            T3 = S2*T2                                                          
            T4 = CD0 + CD1*RHO                                                  
            VIII = T3*T4*EXPV3                                                  
            TDS   = (2.0D0*S*T2 + 3.0D0*S2*T1)*T4                               
            TDCPH = S2*S*S2*T4                                                  
            TDRHO = T3*(CD1 - 2.0D0*T5*T4)                                      
            DV3(1) = (TDS*DS1 + TDCPH*DCPHI1 + TDRHO*DRHO1)*EXPV3               
            DV3(2) = (TDS*DS2 + TDCPH*DCPHI2 + TDRHO*DRHO2)*EXPV3               
            DV3(3) = (TDS*DS3 + TDCPH*DCPHI3 + TDRHO*DRHO3)*EXPV3               
         END IF                                                                 
      END IF                                                                    
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE CAPF(R1,R2,R3)                                                 
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      COMMON/COUPCM/PHI,DPHI(3),SPHI3,DSPHI3(3)                                 
      DIMENSION DX(3), DCAPK(4),DFZ(3),                                         
     +          DFU(3),DGZ(3),F(3,JSURF),RCAPF(3),                              
     +          DS(3), DCPHI(3)                                                 
      RCAPF(1)=R1                                                               
      RCAPF(2)=R2                                                               
      RCAPF(3)=R3                                                               
      CALL H3COOR(RCAPF(1),RCAPF(2),RCAPF(3))                                   
C                                                                               
      DS(1) = DS1                                                               
      DS(2) = DS2                                                               
      DS(3) = DS3                                                               
      DCPHI(1) = DCPHI1                                                         
      DCPHI(2) = DCPHI2                                                         
      DCPHI(3) = DCPHI3                                                         
C                                                                               
      X=SQRT(QCOORD/3.0D0)                                                      
      XI=1.0D0/X                                                                
      XI3=X/3.0D0                                                               
      X2=X*X                                                                    
      CALL KAPPA(X,CAPK,DCAPK)                                                  
      FZ=0.375D0*X*DCAPK(1)                                                     
      FU=0.03515625D0*(X*DCAPK(1)-X2*DCAPK(2)+X2*XI3*DCAPK(3))                  
      GZ=0.046875D0*(-X*DCAPK(1)+X2*DCAPK(2))                                   
      DENOM=FZ+GZ*CPHI3*S+FU*S2                                                 
      ARG=GZ*SPHI3*S/DENOM                                                      
C     FCT=PHI-ATAN(ARG)                                                         
      DO 10 I=1,3                                                               
         DX(I)=RCAPF(I)*XI                                                      
         DFZ(I)=0.125D0*DX(I)*(DCAPK(1)+X*DCAPK(2))                             
         DFU(I)=0.01171875D0*DX(I)*(DCAPK(1)-X*DCAPK(2)+X2*XI3*DCAPK(4))        
         DGZ(I)=0.015625D0*DX(I)*(-DCAPK(1)+X*DCAPK(2)+X2*DCAPK(3))             
         BRAK=S*DGZ(I)+GZ*DS(I)                                                 
         DARG=-ARG*(DFZ(I)+CPHI3*BRAK+S*(GZ*DCPHI(I)+S*DFU(I)+                  
     +         2.0D0*FU*DS(I)))                                                 
         DARG=(DARG+GZ*S*DSPHI3(I)+SPHI3*BRAK)/DENOM                            
         ENGYIJ(I)=DPHI(I)-1.0D0/(1.0D0+ARG**2)*DARG                            
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H2COR2 (RR)                                                    
C                                                                               
C*****************************************************************              
C     CALCULATES 1-4 DERIVATIVES OF 2-BODY CORRELATION TERMS FOR                
C     BOTH SINGLET AND TRIPLET STATES AND RETURNS THEIR DIFFERENCE              
C******************************************************************             
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      COMMON/DMPCM/DPOL,Y,TT,TT1,PARB                                           
      COMMON/COR2CM/CORRS,CORRT,DCOR2(4)                                        
         RSAVE = RR                                                             
         CORRS = 0.0D0                                                          
         DCORRS = 0.0D0                                                         
         CORRT = 0.0D0                                                          
         DCORRT = 0.0D0                                                         
      DO 11 J=1,4                                                               
   11 DCOR2(J)=0.0D0                                                            
         T1 = 1.0D0/RSAVE                                                       
         T = T1*T1                                                              
         T2 = T*T                                                               
         T1 = T2*T1                                                             
         NEXP = 4                                                               
C      Loop over terms in dispersion expansion                                  
         DO 1 ID = 1,3                                                          
            NEXP = NEXP + 2                                                     
            FNEXP = DBLE(NEXP)                                                  
            T2 = T2*T                                                           
            T1 = T1*T                                                           
C                                                                               
C******SINGLET                                                                  
C                                                                               
            DENOM = H2RM + 2.5D0*H2R0                                           
            Y=2.0D0/DENOM                                                       
            X=RSAVE*Y                                                           
            PARB=DAMPB(ID)                                                      
            TT = DAMPB(ID)*X                                                    
            TT1 = DAMPA(ID) + TT                                                
            POL = X*TT1                                                         
            DPOL = TT1 + TT                                                     
            TT = EXP(-POL)                                                      
            TT1 = 1.0D0 - TT                                                    
            TT2 = TT1**(NEXP-1)                                                 
            D = TT1*TT2                                                         
            DD = 2.0D0*FNEXP*DPOL*TT*TT2/DENOM                                  
            CORRS = CORRS - CHH(ID)*D*T2                                        
            DCORRS = DCORRS - CHH(ID)*(DD*T2-FNEXP*D*T1)                        
            CALL DAMP(RSAVE,NEXP,D2,D3,D4)                                      
            TEM2=D2                                                             
            TEM3=D3                                                             
            TEM4=D4                                                             
C                                                                               
C*****TRIPLET                                                                   
C                                                                               
            DENOM = H2RMT + 2.5D0*H2R0                                          
            Y=2.0D0/DENOM                                                       
            X=RSAVE*Y                                                           
            TT = DAMPB(ID)*X                                                    
            TT1 = DAMPA(ID) + TT                                                
            POL = X*TT1                                                         
            DPOL = TT1 + TT                                                     
            TT = EXP(-POL)                                                      
            TT1 = 1.0D0 - TT                                                    
            TT2 = TT1**(NEXP-1)                                                 
            D = TT1*TT2                                                         
            DD = 2.0D0*FNEXP*DPOL*TT*TT2/DENOM                                  
            CORRT = CORRT - CHH(ID)*D*T2                                        
            DCORRT = DCORRT - CHH(ID)*(DD*T2 - FNEXP*D*T1)                      
            CALL DAMP(RSAVE,NEXP,D2,D3,D4)                                      
            DCOR2(2)=DCOR2(2)+CHH(ID)*(D2-TEM2)                                 
            DCOR2(3)=DCOR2(3)+CHH(ID)*(D3-TEM3)                                 
            DCOR2(4)=DCOR2(4)+CHH(ID)*(D4-TEM4)                                 
    3       CONTINUE                                                            
    1    CONTINUE                                                               
      DCOR2(1)=DCORRT-DCORRS                                                    
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE DAMP(RMP,N,D2,D3,D4)                                           
C                                                                               
C**********************************************************************         
C      CALCULATES 2ND,3RD AND 4TH DERIVATIVES OF 2-BODY CORRELATION             
C      ENERGIES FOR SINGLET AND TRIPLET STATES FOR USE IN NON-ADIABATIC         
C      COUPLINGS OF H3 SURFACES 1 AND 2. SENT TO SUBROUTINE H3COR2              
C**********************************************************************         
C                                                                               
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
C      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                
C     +                AL0,AL1,AL2,AL3,AZ2,                                     
C     +                BETA1,BETA2,BETA3,BET0,BET1,                             
C     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                            
C     +                DAMPA(3),DAMPB(3),                                       
C     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                
C     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,            
C     +                RHOL,SQRT3,XPAR(15)                                      
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
C      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                    
C     +                QCOORD,DQ1,DQ2,DQ3,                                      
C     +                RHO,DRHO1,DRHO2,DRHO3,                                   
C     +                S,S2,DS1,DS2,DS3,                                        
C     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                               
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      COMMON/DMPCM/DPOL,Y,TEXP,B,DAMPB                                          
      DPOL2=DPOL**2                                                             
      D1B=Y*DPOL*TEXP                                                           
      D2B=Y**2*(2.0D0*DAMPB-DPOL2)*TEXP                                         
      D3B=Y**3*DPOL*(DPOL2-6.0D0*DAMPB)*TEXP                                    
      D4B=Y**4*(-12.D0*DAMPB*DAMPB+DPOL2*(12.D0*DAMPB-DPOL2))*TEXP              
      T=DBLE(N)                                                                 
      TM=DBLE(N-1)                                                              
      TM2=DBLE(N-2)                                                             
      TM3=DBLE(N-3)                                                             
      TP=DBLE(N+1)                                                              
      TP2=DBLE(N+2)                                                             
      TP3=DBLE(N+3)                                                             
      TH=3.0D0                                                                  
      TW=2.0D0                                                                  
      FR=4.0D0                                                                  
      RR=1.0D0/RMP                                                              
      RRN=RR**N                                                                 
      BM4=B**(N-4)                                                              
      BM3=B*BM4                                                                 
      BM2=B*BM3                                                                 
      D2=T*BM2*(-TP*(B*RR)**2+TW*T*B*RR*D1B-TM*D1B**2-B*D2B)*RRN                
      D3=T*BM3*(TP*TP2*(B*RR)**3-TH*T*TP*(B*RR)**2*D1B+TH*TM*T*B*D1B**2*        
     1RR-TM2*TM*D1B**3+TH*T*RR*B**2*D2B-TH*TM*B*D1B*D2B-B**2*D3B)*RRN           
      D4=T*BM4*(-TP*TP2*TP3*(B*RR)**4+FR*T*TP*TP2*(B*RR)**3*D1B-TW*TH*TM        
     1*T*TP*(B*RR*D1B)**2+FR*TM2*TM*T*B*RR*D1B**3-TW*TH*T*TP*B*(B*RR)**         
     22*D2B+TH*FR*TM*T*B*B*RR*D1B*D2B-TW*TH*TM2*TM*B*D1B*D1B*D2B+FR*T*          
     3RR*B**3*D3B-TM3*TM2*TM*D1B**4-TH*TM*(B*D2B)**2-FR*TM*B*B*D1B              
     4*D3B-B**3*D4B)*RRN                                                        
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE KAPPA(QLC,CAPK,DCAPK)                                          
C                                                                               
C*****************************************************************              
C     CALCULATES FUNCTION CAPK AND ITS FIRST ROUR DERIVATIVES                   
C*****************************************************************              
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
C                                                                               
C                                                                               
      COMMON/COR2CM/CORRS,CORRT,DCOR2(4)                                        
      DIMENSION DCAPK(4)                                                        
      CALL H2COR2(QLC)                                                          
      T1 = -AZ2*QCOORD                                                          
      FZERO=EXP(T1*QCOORD)                                                      
      RR = QLC                                                                  
      RRI = 1.0D0/RR                                                            
C                                                                               
C*****CALCULATE 1-4 DERIVATIVES OF VEHF                                         
C                                                                               
      DR = RR - H2RM                                                            
      T = EXP(-HFGAM*DR)                                                        
      ES = -HFD*(1.0D0 + DR*(HFA1 + DR*(HFA2 + DR*HFA3)))*T                     
      DESDR = -HFGAM*ES - HFD*(HFA1 + DR*(2.0D0*HFA2 +                          
     *      3.0D0*DR*HFA3))*T                                                   
      P1=HFA1+2.0D0*HFA2*DR+3.0D0*HFA3*DR*DR                                    
      P3=6.0D0*HFA3                                                             
      P2=2.0D0*HFA2 +P3*DR                                                      
      D1VHF=DESDR                                                               
      D2VHF=-HFGAM*D1VHF-HFD*T*(P2-HFGAM*P1)                                    
      D3VHF=-HFGAM*D2VHF-HFD*T*(P3-2.D0*HFGAM*P2+HFGAM*HFGAM*P1)                
      D4VHF=-HFGAM*D3VHF+HFD*T*(3.D0*HFGAM*P3-3.D0*HFGAM*HFGAM*P2+              
     1 HFGAM*HFGAM*HFGAM*P1)                                                    
      DCAPK(1)=D1VHF                                                            
      DCAPK(2)=D2VHF                                                            
      DCAPK(3)=D3VHF                                                            
      DCAPK(4)=D4VHF                                                            
C                                                                               
C*****CALCULATE 1-4 DERIVATIVES OF UEHF                                         
C                                                                               
         T = RR*(H2TB1 + RR*(H2TB2 + RR*(H2TB3 + RR*H2TB4)))                    
         AT = H2TA*EXP(-T)*RRI                                                  
         T = H2TB1 + RR*(2.0D0*H2TB2 + RR*(3.0D0*H2TB3 +                        
     *      RR*4.0D0*H2TB4))+RRI                                                
         DAT = -AT*T                                                            
         Q1=-RRI**2+2.D0*H2TB2+6.D0*H2TB3*RR+12.D0*H2TB4*RR**2                  
         Q2=2.D0*RRI**3+6.D0*H2TB3+24.D0*H2TB4*RR                               
         Q3=-6.D0*RRI**4+24.D0*H2TB4                                            
         D1UEHF=DAT                                                             
         D2UEHF=-T*D1UEHF-AT*Q1                                                 
         D3UEHF=-T*D2UEHF-2.D0*Q1*D1UEHF-AT*Q2                                  
         D4UEHF=-T*D3UEHF-3.D0*Q1*D2UEHF-3.D0*Q2*D1UEHF-AT*Q3                   
C                                                                               
C     Add in triplet and subtract singlet 2-body correlation terms              
C                                                                               
         AT = AT + CORRT - CORRS                                                
C                                                                               
C*****INCLUDE CORRELATION PART OF CAPK,(USUBC-VSUBC)*(1-FZERO)                  
C                                                                               
      OMFZ=1.0D0-FZERO                                                          
      DCAPK(1)=DCAPK(1)-(D1UEHF+DCOR2(1))*OMFZ                                  
      DCAPK(2)=DCAPK(2)-(D2UEHF+DCOR2(2))*OMFZ                                  
      DCAPK(3)=DCAPK(3)-(D3UEHF+DCOR2(3))*OMFZ                                  
      DCAPK(4)=DCAPK(4)-(D4UEHF+DCOR2(4))*OMFZ                                  
C                                                                               
C*****CALCULATE 1-4 DERIVATIVES OF W-EFFECTIVE                                  
C                                                                               
         T = EXP(-AL3*RR)                                                       
         WE= (AL0 + RR*(AL1 + RR*AL2))*T                                        
         DWE = (AL1 + 2.0D0*AL2*RR)*T - AL3*WE                                  
      O2=2.0D0*AL2                                                              
      O1=AL1+O2*RR                                                              
      D1W=DWE                                                                   
      D2W=-AL3*D1W+T*(O2-AL3*O1)                                                
      D3W=-AL3*D2W+T*(AL3*O1-2.0D0*O2)*AL3                                      
      D4W=-AL3*D3W+T*(3.0D0*O2-AL3*O1)*AL3**2                                   
      DCAPK(1)=DCAPK(1)-FZERO*D1W                                               
      DCAPK(2)=DCAPK(2)-FZERO*D2W                                               
      DCAPK(3)=DCAPK(3)-FZERO*D3W                                               
      DCAPK(4)=DCAPK(4)-FZERO*D4W                                               
      CAPK=ES-FZERO*WE-OMFZ*AT                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
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
