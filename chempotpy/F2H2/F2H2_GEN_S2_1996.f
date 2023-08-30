      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=4
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
          p(istate)=PENGYGS/8065.7112013
        enddo
      else
        write (*,*) 'Only energy is available'
      endif

      endsubroutine


C                                                                               
      subroutine prepot                                                         
C                                                                               
C   System:          F2H2                                                       
C   Common name:     S2                                                         
C   Number of derivatives: 0                                                    
C   Number of bodies: 4                                                         
C   Number of electronic surfaces: 1                                            
C                                                                               
C   Interface:       potlib2001
C                                                                               
C   Note:            This surface is based on the SQSBDE surface                
C                    by M.A. Suhm.  It was refit to give better                 
C                    tunneling splitting for the ground state.                  
C                                                                               
C   Protocol:                                                                   
C                                                                               
C      PREPOT - initializes the potential's variables                           
C               must be called once before any calls to POT                     
C      POT    - driver for the evaluation of the energy                         
C                                                                               
C   Units:                                                                      
C                                                                               
C      energies    - cm-1                                                       
C      coordinates - distances: bohr; angles: degrees                           
C                                                                               
C   Surfaces:                                                                   
C                                                                               
C      ground electronic state                                                  
C                                                                               
C M. A. SUHM  SEPTEMBER 1990                                                    
C                                                                               
C INPUT : SIX INTERNAL COORDINATES                                              
C --R1,R2,R IN ATOMIC UNITS (BOHR)                                              
C --TH1,TH2,TAU IN RADIANS                                                      
C DEFINITION AS IN QUACK+SUHM MOL.PHYS. 69 (1990) P791                          
C OUTPUT : HF DIMER (SQSBDE) POTENTIAL ENERGY IN HARTREE                        
C                                                                               
C THE PARAMETERS OF THE POTENTIAL ARE READ IN FROM UNIT 1                       
C (FILENAME GA.IN)                                                              
C                                                                               
C THE 6-D POTENTIAL IS BASED ON 1070 AB INITIO POINTS FROM                      
C KOFRANEK ET AL (CHEM PHYS 121(1988) 137) UP TO 42000 CM-1                     
C AND A PERTURBATION TREATMENT BY RIJKS AND WORMER                              
C (JCP 90(1989) 6507) AS WELL AS MONOMER We,WeXe FROM                           
C HUBER/HERZBERG. 29 FREE, 15 MONOMER OR AB INITIO                              
C BASED AND 18 CONSTRAINED PARAMETERS ARE USED.                                 
C SQSBDE IS AN EMPIRICALLY REFINED VERSION WHICH HAS THE                        
C CORRECT EXPERIMENTAL HYDROGEN BOND WAVENUMBER OF                              
C 1062CM-1 (MILLER, ACC. CHEM. RES. 23 (1990) P10) AND                          
C B ROTATIONAL CONSTANT. FOR THIS PURPOSE, THE AB INITIO                        
C RAB COORDINATES WERE SCALED BY 1/1.035.                                       
C THE RMS OF THE WEIGHTED FIT (W=1 FOR R1,2=/=1.7374 OR                         
C R1,2=1.7374 AND E<2220CM-1, W=0.1 FOR R1,2=1.7374 AND                         
C 2220<E<3210CM-1, W=0.001 FOR R1,2=1.7374 AND E>3210CM-1)                      
C IS 32.8CM-1.                                                                  
C                                                                               
C E(1.7445,1.7404,5.1440,9.0007,64.1361,180)=-1.307CM-1 IS                      
C THE ABSOLUTE MINIMUM OF THE POTENTIAL                                         
C                                                                               
C M. A. SUHM  SEPTEMBER 1989                                                    
C                                                                               
C preparatory calculations for gssqs-routine                                    
C                                                                               
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
      parameter (ma=135)                                                        
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
                                                                                
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
C   References:      W.C. Necoechea, D.G. Truhlar Chem. Phys. Lett. 248,        
C                    (1995) 182                                                 
C                    M. Quack, M. Suhm Mol. Phys. 69 (1990) P791                
C                    Kofranek et al Chem. Phys. 121 (1988) 137                  
C                    Rijks and Wormer J. Chem. Phys. 90(1989) 6507              
C                                                                               
       REF(1)='W.C. Necoechea, D.G. Truhlar,'                                   
       REF(2)='Chem. Phys. Lett. 248, 182(1995)'                                
       REF(3)='M. Quack, M. Suhm, Mol. Phys. 69, P791(1990)'                    
       REF(4)='Kofranek et al., Chem. Phys. 121, 137(1988)'                     
       REF(5)='Rijks and Wormer, J. Chem. Phys. 90, 6507(1989)'                 
C                                                                               
      INDEXES(1) = 1                                                            
      INDEXES(2) = 9                                                            
      INDEXES(3) = 1                                                            
      INDEXES(4) = 9                                                            
C                                                                               
      IRCTNT=3                                                                  
C                                                                               
C      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
C prefactors for polynomials                                                    
C                                                                               
      fac(0)=1.                                                                 
      do 101 i=1,2*4                                                            
         fac(i)=fac(i-1)*i                                                      
  101 continue                                                                  
      do 102 l=0,4                                                              
         do 103 m=0,l                                                           
            pf(l,m)=(-1)**m*sqrt((2*l+1)*fac(l-m)/(4*pi*fac(l+m)))              
  103    continue                                                               
  102 continue                                                                  
C                                                                               
C 3j-symbols for potential energy surface                                       
C (abc/d-d0)=g(a,b,c,d)                                                         
C                                                                               
C                                                                               
      fct(1)=1.0d0                                                              
      efct(1)=0.0d0                                                             
      do 104 i=1,400                                                            
         x=i                                                                    
         j=i+1                                                                  
         efct(j)=efct(i)                                                        
         fct(j)=fct(i)*x                                                        
    1    if(fct(j).ge.10.0d0) then                                              
           fct(j)=0.1d0*fct(j)                                                  
           efct(j)=efct(j)+1.0d0                                                
           goto 1                                                               
         endif                                                                  
  104 continue                                                                  
C                                                                               
C                                                                               
      m2=0                                                                      
      m3=0                                                                      
      do 108 j1=0,4*2,2                                                         
         do 107 j2=0,4*2,2                                                      
            do 106 j3=iabs(j1-j2),iabs(j1+j2),4                                 
               do 105 m1=0,min(j1,j2),2                                         
                  m2=-m1                                                        
                  i1=max0(0,(j2-j3-m1)/2,(j1-j3+m2)/2)                          
                  i2=min0((j1+j2-j3)/2,(j1-m1)/2,(j2+m2)/2)                     
                  ia=i1                                                         
                  ib=(j1+j2-j3)/2-i1                                            
                  ic=(j1-m1)/2-i1                                               
                  id=(j2+m2)/2-i1                                               
                  ie=(j3-j2+m1)/2+i1                                            
                  if=(j3-j1-m2)/2+i1                                            
                  sum=(-1)**i1/(fct(ia+1)*fct(ib+1)*fct(ic+1)*fct(id+1)*        
     *                fct(ie+1)*fct(if+1))                                      
                  esum=-efct(ia+1)-efct(ib+1)-efct(ic+1)-efct(id+1)-            
     *                 efct(ie+1)-efct(if+1)                                    
                  nmax=i2-i1                                                    
                  f=1.0d0                                                       
                  if(nmax.gt.0) then                                            
                     do 109 i=nmax,1,-1                                         
                        x=(ib-i+1)*(ic-i+1)*(id-i+1)                            
                        x=x/((ia+i)*(ie+i)*(if+i))                              
                        f=1.0d0-x*f                                             
  109                continue                                                   
                  endif                                                         
                  sum=sum*f                                                     
                  f=fct((j1+j2-j3)/2+1)*fct((j1-j2+j3)/2+1)*                    
     *              fct((j2-j1+j3)/2+1)*fct((j1+m1)/2+1)*                       
     *              fct((j1-m1)/2+1)*                                           
     *              fct((j2+m2)/2+1)*fct((j2-m2)/2+1)*fct((j3+m3)/2+1)*         
     *              fct((j3-m3)/2+1)/fct((j1+j2+j3)/2+2)                        
                  e=efct((j1+j2-j3)/2+1)+efct((j1-j2+j3)/2+1)+                  
     *              efct((j2-j1+j3)/2+1)+efct((j1+m1)/2+1)+                     
     *              efct((j1-m1)/2+1)+                                          
     *              efct((j2+m2)/2+1)+efct((j2-m2)/2+1)+                        
     *              efct((j3+m3)/2+1)+                                          
     *              efct((j3-m3)/2+1)-efct((j1+j2+j3)/2+2)                      
                  cg3j=sum*dsqrt(f)*10.0d0**(esum+0.5d0*e)*                     
     *                 (-1)**((j2-j1-m3)/2)                                     
                  g(j1/2,j2/2,j3/2,m1/2)=cg3j                                   
  105          continue                                                         
  106       continue                                                            
  107    continue                                                               
  108 continue                                                                  
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
C                                                                               
C M. A. SUHM  SEPTEMBER 1989                                                    
C                                                                               
      subroutine gsrsqs                                                         
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      parameter (ndmax=1100)                                                    
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
      parameter (ma=135)                                                        
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
C      dab1=rcm1-3.0d0                                                          
C      dab2=rcm1-4.0d0                                                          
C      dab3=rcm1-4.9d0                                                          
C                                                                               
      dab1=R(1)-3.0d0                                                           
      dab2=R(1)-4.0d0                                                           
      dab3=R(1)-4.09d0                                                          
      do 301 j1=0,4,1                                                           
         do 302 j2=j1,4,1                                                       
            do 303 j3=(j2-j1),(j1+j2),2                                         
               ar(j1,j2,j3)=ac(1,j1,j2,j3)*exp(-0.20d0*dab1**2)+                
     *         ac(2,j1,j2,j3)*exp(-0.35d0*dab2**2)+                             
     *         ac(3,j1,j2,j3)*exp(-0.80d0*dab3**2)                              
               ar(j2,j1,j3)=(-1)**(j1+j2)*ar(j1,j2,j3)                          
  303       continue                                                            
  302    continue                                                               
  301 continue                                                                  
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
C                                                                               
C M. A. SUHM  SEPTEMBER 1989                                                    
C                                                                               
      subroutine pot                                                            
C                                                                               
      implicit double precision(a-h,o-z)                                        
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
      parameter (ma=135)                                                        
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
      dimension x(4),y(4),z(4)                                                  
      dimension p1(0:4,0:4),p2(0:4,0:4),cs(4)                                   
      dimension vp(0:4,0:4,0:2*4)                                               
C                                                                               
      data amh/1837.152756d0/amf/34631.9678d0/                                  
C                                                                               
C     PUT COORDINATES IN PROPER ARRAYS                                          
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C*********************************                                              
C                                *                                              
C   Assign entries from R-array  *                                              
C   to variables used internally *                                              
C   by the POT subroutine        *                                              
C                                *                                              
C*********************************                                              
C                                                                               
      rcm=R(1)                                                                  
      r1=R(2)                                                                   
      r2=R(3)                                                                   
      th1=R(4)                                                                  
      th2=R(5)                                                                  
      tau=R(6)                                                                  
C         rcm=2.5d0                                                             
C         r1=0.9d0                                                              
C         r2=0.9d0                                                              
C         th1=1.04719755120d0                                                   
C         th2=1.04719755120d0                                                   
C         tau=0.78539816340d0                                                   
      phicon=0.01745329252d0                                                    
      pif=1.0d0/(2.0d0*sqrt(pi))                                                
      fac2=amh/(amh+amf)                                                        
      x(3)=0.0d0-r1*cos(th1)*fac2                                               
      y(3)=0.0d0-r1*sin(th1)*fac2                                               
      z(3)=0.0d0                                                                
      x(1)=r1*cos(th1)*(1.0d0-fac2)                                             
      y(1)=r1*sin(th1)*(1.0d0-fac2)                                             
      z(1)=0.0d0                                                                
      x(4)=rcm-r2*cos(th2)*fac2                                                 
      y(4)=0.0d0-r2*sin(th2)*cos(tau)*fac2                                      
      z(4)=0.0d0-r2*sin(th2)*sin(tau)*fac2                                      
      x(2)=rcm+r2*cos(th2)*(1.0d0-fac2)                                         
      y(2)=0.0d0+r2*sin(th2)*cos(tau)*(1.0d0-fac2)                              
      z(2)=0.0d0+r2*sin(th2)*sin(tau)*(1.0d0-fac2)                              
      r12=sqrt((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)                    
      r23=sqrt((x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2)                    
      r14=sqrt((x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2)                    
      r34=sqrt((x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2)                    
      do 101 j=1,4                                                              
      cs(j)=cos(j*tau)                                                          
  101 continue                                                                  
      s=sin(th1)                                                                
      c=cos(th1)                                                                
      s2=s*s                                                                    
      s3=s2*s                                                                   
      s4=s3*s                                                                   
      s5=s4*s                                                                   
      c2=c*c                                                                    
      c3=c2*c                                                                   
      c4=c3*c                                                                   
      c5=c4*c                                                                   
      p1(0,0)=pf(0,0)                                                           
      p1(1,0)=pf(1,0)*c                                                         
      p1(2,0)=pf(2,0)*0.5d0*(3.0d0*c2-1.0d0)                                    
      p1(3,0)=pf(3,0)*0.5d0*(5.0d0*c3-3.0d0*c)                                  
      p1(4,0)=pf(4,0)*0.125d0*(35.0d0*c4-30.0d0*c2+3.0d0)                       
      p1(1,1)=pf(1,1)*s                                                         
      p1(2,1)=pf(2,1)*3.0d0*s*c                                                 
      p1(3,1)=pf(3,1)*1.5d0*s*(5.0d0*c2-1.0d0)                                  
      p1(4,1)=pf(4,1)*2.5d0*s*(7.0d0*c3-3.0d0*c)                                
      p1(2,2)=pf(2,2)*3.0d0*s2                                                  
      p1(3,2)=pf(3,2)*15.0d0*s2*c                                               
      p1(4,2)=pf(4,2)*7.5d0*s2*(7.0d0*c2-1.0d0)                                 
      p1(3,3)=pf(3,3)*15.0d0*s3                                                 
      p1(4,3)=pf(4,3)*105.0d0*s3*c                                              
      p1(4,4)=pf(4,4)*105.0d0*s4                                                
      s=sin(th2)                                                                
      c=cos(th2)                                                                
      s2=s*s                                                                    
      s3=s2*s                                                                   
      s4=s3*s                                                                   
      s5=s4*s                                                                   
      c2=c*c                                                                    
      c3=c2*c                                                                   
      c4=c3*c                                                                   
      c5=c4*c                                                                   
      p2(0,0)=pf(0,0)                                                           
      p2(1,0)=pf(1,0)*c                                                         
      p2(2,0)=pf(2,0)*0.5d0*(3.0d0*c2-1.0d0)                                    
      p2(3,0)=pf(3,0)*0.5d0*(5.0d0*c3-3.0d0*c)                                  
      p2(4,0)=pf(4,0)*0.125d0*(35.0d0*c4-30.0d0*c2+3.0d0)                       
      p2(1,1)=pf(1,1)*s                                                         
      p2(2,1)=pf(2,1)*3.0d0*s*c                                                 
      p2(3,1)=pf(3,1)*1.5d0*s*(5.0d0*c2-1.0d0)                                  
      p2(4,1)=pf(4,1)*2.5d0*s*(7.0d0*c3-3.0d0*c)                                
      p2(2,2)=pf(2,2)*3.0d0*s2                                                  
      p2(3,2)=pf(3,2)*15.0d0*s2*c                                               
      p2(4,2)=pf(4,2)*7.5d0*s2*(7.0d0*c2-1.0d0)                                 
      p2(3,3)=pf(3,3)*15.0d0*s3                                                 
      p2(4,3)=pf(4,3)*105.0d0*s3*c                                              
      p2(4,4)=pf(4,4)*105.0d0*s4                                                
      icu=31                                                                    
      do 102 j1=0,4,1                                                           
         do 103 j2=j1,4,1                                                       
            do 104 j3=(j2-j1),(j1+j2),2                                         
               ac(1,j1,j2,j3)=a(icu)                                            
               ac(2,j1,j2,j3)=a(icu+35)                                         
               ac(3,j1,j2,j3)=a(icu+70)                                         
               icu=icu+1                                                        
  104       continue                                                            
  103    continue                                                               
  102 continue                                                                  
C                                                                               
      call gsrsqs                                                               
C                                                                               
      do 105 j1=0,4,1                                                           
         do 106 j2=0,4,1                                                        
            do 107 j3=abs(j2-j1),(j2+j1),2                                      
               vp(j1,j2,j3)=g(j1,j2,j3,0)*p1(j1,0)*p2(j2,0)                     
               do 108 m=1,min(j1,j2),1                                          
                  vp(j1,j2,j3)=vp(j1,j2,j3)+(-1)**m*2*g(j1,j2,j3,m)*            
     *            cs(m)*p1(j1,m)*p2(j2,m)                                       
  108          continue                                                         
               vp(j1,j2,j3)=vp(j1,j2,j3)*(-1)**(j1-j2)*(2*j3+1)*pif             
  107       continue                                                            
  106    continue                                                               
  105 continue                                                                  
C                                                                               
      e=0.0d0                                                                   
C                                                                               
      do 109 j1=0,4,1                                                           
         do 110 j2=0,4,1                                                        
            do 111 j3=abs(j2-j1),(j2+j1),2                                      
               e=e+vp(j1,j2,j3)*ar(j1,j2,j3)                                    
  111       continue                                                            
  110    continue                                                               
  109 continue                                                                  
C                                                                               
c     e=e+a(1)                                                                  
C                                                                               
      if(rcm.lt.7.2d0)then                                                      
         e=e+vp(0,0,0)*exp(-((7.2d0-rcm)/rcm)**2)*a(2)/rcm**6                   
      else                                                                      
         e=e+vp(0,0,0)*a(2)/rcm**6                                              
      endif                                                                     
      if(rcm.lt.8.4d0)then                                                      
         e=e+vp(0,0,0)*exp(-((8.4d0-rcm)/rcm)**2)*a(3)/rcm**8                   
      else                                                                      
         e=e+vp(0,0,0)*a(3)/rcm**8                                              
      endif                                                                     
      if(rcm.lt.9.6d0)then                                                      
         e=e+vp(0,0,0)*exp(-((9.6d0-rcm)/rcm)**2)*a(4)/rcm**10                  
      else                                                                      
         e=e+vp(0,0,0)*a(4)/rcm**10                                             
      endif                                                                     
      if(rcm.lt.9.5d0) then                                                     
         e=e+vp(1,1,2)*a(5)/rcm**3*exp(-((9.5d0-rcm)/rcm)**2)                   
      else                                                                      
         e=e+vp(1,1,2)*a(5)/rcm**3                                              
      endif                                                                     
      if(rcm.lt.10.5d0) then                                                    
         e=e+(vp(2,1,3)-vp(1,2,3))*a(6)/rcm**4*                                 
     *     exp(-((10.5d0-rcm)/rcm)**2)                                          
      else                                                                      
         e=e+(vp(2,1,3)-vp(1,2,3))*a(6)/rcm**4                                  
      endif                                                                     
      if(rcm.lt.11.5d0) then                                                    
         e=e+vp(2,2,4)*a(7)/rcm**5*exp(-((11.5d0-rcm)/rcm)**2)                  
         e=e+(vp(3,1,4)+vp(1,3,4))*a(8)/rcm**5*                                 
     *     exp(-((11.5d0-rcm)/rcm)**2)                                          
      else                                                                      
         e=e+vp(2,2,4)*a(7)/rcm**5                                              
         e=e+(vp(3,1,4)+vp(1,3,4))*a(8)/rcm**5                                  
      endif                                                                     
      e=e+a(9)*exp(-0.5d0*r12)                                                  
      e=e+a(10)*exp(-1.5d0*r12)                                                 
      e=e+a(13)*exp(-1.5d0*r34)                                                 
      e=e+a(17)*(exp(-10.0d0*r14)+exp(-10.0d0*r23))                             
      if(rcm.lt.4.5d0)then                                                      
         e=e+a(18)/rcm**6*((3.0d0*(cos(th1))**2+1.0d0)+                         
     *     (3.0d0*(cos(th2))**2+1.0d0))*                                        
     *     exp(-((4.5d0-rcm)/rcm)**2)                                           
      else                                                                      
         e=e+a(18)/rcm**6*((3.0d0*(cos(th1))**2+1.0d0)+                         
     *     (3.0d0*(cos(th2))**2+1.0d0))                                         
      endif                                                                     
C                                                                               
C monomer terms                                                                 
C                                                                               
      e=e+47635.0d0*(1.0d0-a(19)*exp(-0.1d0*(r23-3.0d0)**2))*                   
     *         (1.0d0-exp(-1.1953d0*                                            
     *         (1.0d0-a(26)*(1.0d0-tanh(a(27)*(r14-3.0d0)))                     
     *                                                      )*                  
     *         (r1-1.7374d0*(1.0d0                                              
     *         -a(20)*                                                          
     *         (1.0d0-tanh(a(21)*(r14-3.0d0)))                                  
     *                                        ))))**2                           
      e=e+47635.0d0*(1.0d0-a(19)*exp(-0.1d0*(r14-3.0d0)**2))*                   
     *         (1.0d0-exp(-1.1953d0*                                            
     *         (1.0d0-a(26)*(1.0d0-tanh(a(27)*(r23-3.0d0)))                     
     *                                                      )*                  
     *         (r2-1.7374d0*(1.0d0                                              
     *         -a(20)*                                                          
     *         (1.0d0-tanh(a(21)*(r23-3.0d0)))                                  
     *                                        ))))**2                           
c                                                                               
c     next line corrects for asymptotic diatomic potential which is             
c     added in separately.                                                      
c                                                                               
      e=e-(ptntl(r1-1.7374d0)+ptntl(r2-1.7374d0))*219474.63067d0                
      e=e+a(28)*(exp(-4.0d0*(r14-r1)**2)+exp(-4.0d0*(r23-r2)**2))               
      e=e+a(22)*exp(-4.0d0*r14)*exp(-4.0d0*r23)                                 
      e=e+a(23)                                                                 
      e=e+a(24)*(exp(-r14)+exp(-r23))*                                          
     *  (exp(-15.0d0*r1)+exp(-15.0d0*r2))                                       
C                                                                               
C   Units conversion:  Conversion to atomic units                               
C                                                                               
      e=e/219474.63067d0                                                        
C                                                                               
      ENGYGS=e                                                                  
C                                                                               
C                                                                               
      CALL EUNITZERO                                                            
C                                                                               
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      return                                                                    
      end                                                                       
C                                                                               
      FUNCTION PTNTL(X)                                                         
C                                                                               
C     HF DIATOMIC POTENTIAL CURVE FOR SQSBDE SURFACE, CALCULATED USING          
C     PARAMETERS GIVEN IN JCP, VOL. 95(5), P. 31.                               
C     NOTICE THAT THE PARAMETER DM IS GIVEN IN ATOMIC UNITS HERE AND            
C     IN CM-1 IN THE JCP REFERENCE.                                             
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
      parameter (ma=135)                                                        
      PARAMETER (AM=1.1953D0,DM=0.21704D0,ONE=1.0D0)                            
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
                                                                                
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
C      SAVE                                                                     
C                                                                               
      PTNTL=DM*(ONE-EXP(-AM*X))**2                                              
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
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
      parameter (ma=135)                                                        
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
       DATA NASURF /1,35*0/                                                     
       DATA NDER /0/                                                            
       DATA NFLAG /1,1,15*0,6,0,0/                                              
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/4,0,0/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /4/                                                           
C                                                                               
      DATA A(1),DA(1)  / 1558.0D0,  0.0D0/                                      
      DATA A(2),DA(2)  / -202900000.0D0,  0.0D0/                                
      DATA A(3),DA(3)  / -3703000000.0D0,  0.0D0/                               
      DATA A(4),DA(4)  / -11350000000.0D0,  0.0D0/                              
      DATA A(5),DA(5)  / -1736693.094999999D0,  0.0D0/                          
      DATA A(6),DA(6)  / -4416487.722000003D0,  0.0D0/                          
      DATA A(7),DA(7)  / 16012389.50999999D0,  0.0D0/                           
      DATA A(8),DA(8)  / -6389444.854999989D0,  0.0D0/                          
      DATA A(9),DA(9)  / -866.86737555133D0,  117.5897899275601D0/              
      DATA A(10),DA(10)  / 129789.3634701688D0,  1937.563168749279D0/           
      DATA A(11),DA(11)  / 0.0D0,  0.0D0/                                       
      DATA A(12),DA(12)  / 0.0D0,  0.0D0/                                       
      DATA A(13),DA(13)  / 2622184.756658763D0,  29726.63014604442D0/           
      DATA A(14),DA(14)  / 0.0D0,  0.0D0/                                       
      DATA A(15),DA(15)  / 0.0D0,  0.0D0/                                       
      DATA A(16),DA(16)  / 0.0D0,  0.0D0/                                       
      DATA A(17),DA(17)  / 184254701843.3105D0,  17553092934.44629D0/           
      DATA A(18),DA(18)  / -320000.0D0,  0.0D0/                                 
      DATA A(19),DA(19)  / 3.326438943860243D-2,  1.353645657747024D-3/         
      DATA A(20),DA(20)  / -4.576887206887037D-3,  1.865017706495818D-4/        
      DATA A(21),DA(21)  / 0.4097306660396534D0,  2.777020547964315D-2/         
      DATA A(22),DA(22)  / 32623551943283.75D0,  3630087056915.047D0/           
      DATA A(23),DA(23)  / 0.0D0,  0.0D0/                                       
      DATA A(24),DA(24)  / 71754497894202.0D0,  2525041881368.016D0/            
      DATA A(25),DA(25)  / 0.0D0,  0.0D0/                                       
      DATA A(26),DA(26)  / 5.032387872479127D-2,  1.860431838521345D-3/         
      DATA A(27),DA(27)  / 1.46832057915281D0,  9.71324803167759D-2/            
      DATA A(28),DA(28)  / 3384.631559286528D0,  85.55678229458908D0/           
      DATA A(29),DA(29)  / 0.0D0,  0.0D0/                                       
      DATA A(30),DA(30)  / 0.0D0,  0.0D0/                                       
      DATA A(31),DA(31)  / -160267.101272529D0,  4055.338776717705D0/           
      DATA A(32),DA(32)  / 0.0D0,  0.0D0/                                       
      DATA A(33),DA(33)  / 0.0D0,  0.0D0/                                       
      DATA A(34),DA(34)  / 0.0D0,  0.0D0/                                       
      DATA A(35),DA(35)  / 0.0D0,  0.0D0/                                       
      DATA A(36),DA(36)  / -16268.89465981041D0,  541.1396773442175D0/          
      DATA A(37),DA(37)  / 0.0D0,  0.0D0/                                       
      DATA A(38),DA(38)  / 0.0D0,  0.0D0/                                       
      DATA A(39),DA(39)  / 2802.576101972474D0,  56.80382831682482D0/           
      DATA A(40),DA(40)  / 0.0D0,  0.0D0/                                       
      DATA A(41),DA(41)  / 0.0D0,  0.0D0/                                       
      DATA A(42),DA(42)  / 0.0D0,  0.0D0/                                       
      DATA A(43),DA(43)  / 0.0D0,  0.0D0/                                       
      DATA A(44),DA(44)  / 0.0D0,  0.0D0/                                       
      DATA A(45),DA(45)  / 0.0D0,  0.0D0/                                       
      DATA A(46),DA(46)  / 2860.063943395988D0,  44.36333848084973D0/           
      DATA A(47),DA(47)  / 0.0D0,  0.0D0/                                       
      DATA A(48),DA(48)  / 0.0D0,  0.0D0/                                       
      DATA A(49),DA(49)  / 0.0D0,  0.0D0/                                       
      DATA A(50),DA(50)  / 0.0D0,  0.0D0/                                       
      DATA A(51),DA(51)  / 0.0D0,  0.0D0/                                       
      DATA A(52),DA(52)  / 509.1020221188392D0,  18.73758800462167D0/           
      DATA A(53),DA(53)  / 0.0D0,  0.0D0/                                       
      DATA A(54),DA(54)  / 0.0D0,  0.0D0/                                       
      DATA A(55),DA(55)  / 0.0D0,  0.0D0/                                       
      DATA A(56),DA(56)  / 0.0D0,  0.0D0/                                       
      DATA A(57),DA(57)  / 0.0D0,  0.0D0/                                       
      DATA A(58),DA(58)  / 0.0D0,  0.0D0/                                       
      DATA A(59),DA(59)  / 0.0D0,  0.0D0/                                       
      DATA A(60),DA(60)  / 0.0D0,  0.0D0/                                       
      DATA A(61),DA(61)  / 0.0D0,  0.0D0/                                       
      DATA A(62),DA(62)  / 0.0D0,  0.0D0/                                       
      DATA A(63),DA(63)  / 0.0D0,  0.0D0/                                       
      DATA A(64),DA(64)  / 0.0D0,  0.0D0/                                       
      DATA A(65),DA(65)  / 0.0D0,  0.0D0/                                       
      DATA A(66),DA(66)  / 0.0D0,  0.0D0/                                       
      DATA A(67),DA(67)  / 0.0D0,  0.0D0/                                       
      DATA A(68),DA(68)  / -3467.234170178781D0,  52.91056956195189D0/          
      DATA A(69),DA(69)  / 0.0D0,  0.0D0/                                       
      DATA A(70),DA(70)  / 0.0D0,  0.0D0/                                       
      DATA A(71),DA(71)  / 0.0D0,  0.0D0/                                       
      DATA A(72),DA(72)  / 0.0D0,  0.0D0/                                       
      DATA A(73),DA(73)  / 2304.632858038298D0,  44.02327215328341D0/           
      DATA A(74),DA(74)  / 0.0D0,  0.0D0/                                       
      DATA A(75),DA(75)  / 0.0D0,  0.0D0/                                       
      DATA A(76),DA(76)  / 0.0D0,  0.0D0/                                       
      DATA A(77),DA(77)  / 0.0D0,  0.0D0/                                       
      DATA A(78),DA(78)  / 0.0D0,  0.0D0/                                       
      DATA A(79),DA(79)  / 0.0D0,  0.0D0/                                       
      DATA A(80),DA(80)  / 0.0D0,  0.0D0/                                       
      DATA A(81),DA(81)  / 0.0D0,  0.0D0/                                       
      DATA A(82),DA(82)  / 0.0D0,  0.0D0/                                       
      DATA A(83),DA(83)  / 0.0D0,  0.0D0/                                       
      DATA A(84),DA(84)  / -644.8150830210325D0,  18.1151965381963D0/           
      DATA A(85),DA(85)  / 0.0D0,  0.0D0/                                       
      DATA A(86),DA(86)  / 0.0D0,  0.0D0/                                       
      DATA A(87),DA(87)  / 0.0D0,  0.0D0/                                       
      DATA A(88),DA(88)  / 0.0D0,  0.0D0/                                       
      DATA A(89),DA(89)  / 0.0D0,  0.0D0/                                       
      DATA A(90),DA(90)  / 0.0D0,  0.0D0/                                       
      DATA A(91),DA(91)  / 0.0D0,  0.0D0/                                       
      DATA A(92),DA(92)  / 0.0D0,  0.0D0/                                       
      DATA A(93),DA(93)  / 0.0D0,  0.0D0/                                       
      DATA A(94),DA(94)  / 0.0D0,  0.0D0/                                       
      DATA A(95),DA(95)  / 0.0D0,  0.0D0/                                       
      DATA A(96),DA(96)  / 0.0D0,  0.0D0/                                       
      DATA A(97),DA(97)  / 0.0D0,  0.0D0/                                       
      DATA A(98),DA(98)  / 0.0D0,  0.0D0/                                       
      DATA A(99),DA(99)  / 0.0D0,  0.0D0/                                       
      DATA A(100),DA(100)  / 0.0D0,  0.0D0/                                     
      DATA A(101),DA(101)  / 21489.15561258642D0,  989.768118754484D0/          
      DATA A(102),DA(102)  / 1895.153684423094D0,  111.7819362295832D0/         
      DATA A(103),DA(103)  / 0.0D0,  0.0D0/                                     
      DATA A(104),DA(104)  / 768.3954452456346D0,  18.43253058647917D0/         
      DATA A(105),DA(105)  / -395.2301759253078D0,  12.84348126946628D0/        
      DATA A(106),DA(106)  / 2041.09421622329D0,  246.1094667134039D0/          
      DATA A(107),DA(107)  / -2422.952146150754D0,  70.83017811025411D0/        
      DATA A(108),DA(108)  / 0.0D0,  0.0D0/                                     
      DATA A(109),DA(109)  / 0.0D0,  0.0D0/                                     
      DATA A(110),DA(110)  / -297.3729718893483D0,  19.35802742071735D0/        
      DATA A(111),DA(111)  / 0.0D0,  0.0D0/                                     
      DATA A(112),DA(112)  / 0.0D0,  0.0D0/                                     
      DATA A(113),DA(113)  / 39.46815694764132D0,  11.34523118392423D0/         
      DATA A(114),DA(114)  / 0.0D0,  0.0D0/                                     
      DATA A(115),DA(115)  / 720.5552230122485D0,  33.69084291038598D0/         
      DATA A(116),DA(116)  / 0.0D0,  0.0D0/                                     
      DATA A(117),DA(117)  / 0.0D0,  0.0D0/                                     
      DATA A(118),DA(118)  / 0.0D0,  0.0D0/                                     
      DATA A(119),DA(119)  / 0.0D0,  0.0D0/                                     
      DATA A(120),DA(120)  / 0.0D0,  0.0D0/                                     
      DATA A(121),DA(121)  / 0.0D0,  0.0D0/                                     
      DATA A(122),DA(122)  / 0.0D0,  0.0D0/                                     
      DATA A(123),DA(123)  / 0.0D0,  0.0D0/                                     
      DATA A(124),DA(124)  / 0.0D0,  0.0D0/                                     
      DATA A(125),DA(125)  / 0.0D0,  0.0D0/                                     
      DATA A(126),DA(126)  / 0.0D0,  0.0D0/                                     
      DATA A(127),DA(127)  / 0.0D0,  0.0D0/                                     
      DATA A(128),DA(128)  / 0.0D0,  0.0D0/                                     
      DATA A(129),DA(129)  / 0.0D0,  0.0D0/                                     
      DATA A(130),DA(130)  / 0.0D0,  0.0D0/                                     
      DATA A(131),DA(131)  / 0.0D0,  0.0D0/                                     
      DATA A(132),DA(132)  / 0.0D0,  0.0D0/                                     
      DATA A(133),DA(133)  / 0.0D0,  0.0D0/                                     
      DATA A(134),DA(134)  / 0.0D0,  0.0D0/                                     
      DATA A(135),DA(135)  / 0.0D0,  0.0D0/                                     
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
