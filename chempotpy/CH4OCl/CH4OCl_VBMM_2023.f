      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=7
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
        CART(iatom,idir)=x(iatom,idir)
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
      else 
        write (*,*) 'Only energy is available'
      endif

      endsubroutine

c:*************************************************************************
c     SUBROUTINE SURF(V, COORD, DX, N3TM)
      SUBROUTINE POT
C
C   System:    CH3OH + Cl --> CH2OH + HCl
c                         --> CH3O + HCl
c
C              Cipriano Rangel and Joaquin Espinosa-Garcia 2023
C
C              
C   All the information passed to and from the potential energy surface 
C   routine is in hartree atomic units.  
C
C        This potential is written such that:
C                       X(1)  - X(3)  : X, Y, Z for C
C                       X(4)  - X(6)  : X, Y, Z for O 
C                       X(7)  - X(9)  : X, Y, Z for H1
C                       X(10) - X(12) : X, Y, Z for H2
C                       X(13) - X(15) : X, Y, Z for H3
C                       X(16) - X(18) : X, Y, Z for H4 
C                       X(19) - X(21) : X, Y, Z for Cl
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION COORD(21),DX(21)
c     DIMENSION dxvstr(21),dxvtor(21),dxvop(21),dxvip(21)
C
c     include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
c Start of changes by EMN.
c The lines below do not seem to be needed because the Cartesian
c coordinates are in common qpdot.
c Besides, BOHR_A in utility.f does not match the value below.
c As a consequence, q(i) changes and energy is not conserved
cc call de venus
       CALL CARTOU
       CALL CARTTOR
C
C     PUT COORDINATES IN PROPER ARRAYS
C
C Inicialization de variables
C
      DO I=1,27
          q(I) = R(I)
          pdot(I)=0.D0
      ENDDO
c
C Changing units from angstrom to bohr
C
c      DO  I = 1, 27
cc         q(I) = COORD(I)*0.52918d0
c        q(I) = R(I)*0.52918d0
c      enddo
c
c End of corrections by EMN
c
c  calculate relative coordinates and bond lengths
c
       en=0.0d0
       call coorden
c
       call torsion(vtor) 
c
c  calculate switching functions
c
       call switchf
c
c  calculate reference angles and their derivatives
c
       call refangles
c
c  calculate stretching potential
c
       call stretch(vstr)
c
c  calculate out of plane bending potential
c
       call opbend(vop)
c
c  calculate in plane bending potential
c
       call ipbend(vip)
c
c  total potential energy is vstr+vop+vip
c
       en=vstr+vop+vip+vtor
C 
c  convert from 10(5) j/mol to au
c
       en = en*0.03812D0
       V = en 
       ENGYGS = en
c
      CALL EUNITZERO
      IF(NDER.NE.0) THEN
c
c  Initialize derivatives
c
       do n=1,21
c          DX(n)=0.0d0
           DEGSDR(i)=pdot(i)*0.0201723d0
       enddo

         CALL RTOCART
         CALL DEDCOU
      ENDIF
c
       return
       end
c******************************************************
c
c------------------------------------------------------
       subroutine coorden
c------------------------------------------------------
c
c  calculates relative coordinates and bond lengths
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
c  calculate relative coordinates
c
       do ind=1,3
         tc1b(ind)=q(nc1(ind))-q(nhb(ind))
         tc2b(ind)=q(nc2(ind))-q(nhb(ind))
         tcc(ind)=q(nc1(ind))-q(nc2(ind))
       enddo
c
       do ind=1,3
         do i=1,3
           tc1h(i,ind)=q(nc1(ind))-q(n1h(i,ind))
         enddo
       enddo
c
       do ind=1,3
           tc2h(1,ind)=q(nc2(ind))-q(n2h(1,ind))
       enddo
c
       do ind=1,3
         do i=1,3
           tb1h(i,ind)=q(nhb(ind))-q(n1h(i,ind))
         enddo
           tb2h(1,ind)=q(nhb(ind))-q(n2h(1,ind))
       enddo
c
       rc1b=sqrt(tc1b(1)*tc1b(1)+tc1b(2)*tc1b(2)+tc1b(3)*tc1b(3))
       rc2b=sqrt(tc2b(1)*tc2b(1)+tc2b(2)*tc2b(2)+tc2b(3)*tc2b(3))
       rcc=sqrt(tcc(1)*tcc(1)+tcc(2)*tcc(2)+tcc(3)*tcc(3))
c
       do i=1,3
         rc1h(i)=sqrt(tc1h(i,1)*tc1h(i,1)+tc1h(i,2)*tc1h(i,2)+
     *                tc1h(i,3)*tc1h(i,3))
       enddo
c
         rc2h(1)=sqrt(tc2h(1,1)*tc2h(1,1)+tc2h(1,2)*tc2h(1,2)+
     *                tc2h(1,3)*tc2h(1,3))
c
       do i=1,3
         rb1h(i)=sqrt(tb1h(i,1)*tb1h(i,1)+tb1h(i,2)*tb1h(i,2)+
     *                tb1h(i,3)*tb1h(i,3))
       enddo
c
         rb2h(1)=sqrt(tb2h(1,1)*tb2h(1,1)+tb2h(1,2)*tb2h(1,2)+
     *                tb2h(1,3)*tb2h(1,3))
c
c crd 2013 modificaciones de las distancias
c
c
c CRD 2021 modificacion para hacer solo un parametro que
c modifique todos los r0c en funcion de rav
c
       argmax=19.d0
       rav4=(rc1h(1)+rc1h(2)+rc1h(3)+rc2h(1))/4.0D0
       P5=1.d0
          argp5=(w9*(rav4-w10))
          if (argp5.lt.argmax) then
            t5tmp=1-tanh(argp5)
          else
            t5tmp=0.d0
          endif
          P5=P5*t5tmp
       r0cc=P5*r0ccr+(1-P5)*r0ccp
       r0c1h=P5*r0c1hr+(1-P5)*r0c1hp
       r0c2h=P5*r0c2hr+(1-P5)*r0c2hp
c      r0cc= 1.4200d0
c Fin CRD
c       argmax=19.d0
c       P3=1.d0
c       do i=1,3
c          argp3=(w3*(rc1h(i)-w4))
c          if (argp3.lt.argmax) then
c            t3tmp=1-tanh(argp3)
c          else
c            t3tmp=0.d0
c          endif
c          P3=P3*t3tmp
c       enddo
c       r0c1h=P3*r0c1hr+(1-P3)*r0c1hp
cc
cc Dejamos las distancias fijas
cc
cc      r0c1h=r0c1hr
cc
c       P4=1.d0
c          argp4=(w5*(rc2h(1)-w6))
c          if (argp4.lt.argmax) then
c            t4tmp=1-tanh(argp4)
c          else
c            t4tmp=0.d0
c          endif
c          P4=P4*t4tmp
c       r0c2h=P4*r0c2hr+(1-P4)*r0c2hp
cc
cc Dejamos las distancias fijas
cc
cc      r0c2h=r0c2hr
cc
c       P5=1.d0
c          argp5=(w9*(rav4-w10))
c          if (argp5.lt.argmax) then
c            t5tmp=1-tanh(argp5)
c          else
c            t5tmp=0.d0
c          endif
c          P5=P5*t5tmp
c       r0cc=P5*r0ccr+(1-P5)*r0ccp
cc
       return
       end
c
c------------------------------------------------------
c CRD 2011 subrutina para calcular el potencial de torsion
       subroutine torsion(vtor)
c------------------------------------------------------
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
       dimension t1(3)
c      double precision pi
c      pi=4.0d0*atan(1.0d0)
c
c CRD 2011 incluimos un apartado para calcular los angulos dihedro
c para poder calcular despues vtorsion
c
       call dihedro(3,1,2,6,dihed)
       thetahcch(1)=dihed*pi/180.0d0
       call dihedro(4,1,2,6,dihed)
       thetahcch(2)=dihed*pi/180.0d0
       call dihedro(5,1,2,6,dihed)
       thetahcch(3)=dihed*pi/180.0d0
c
c CRD 2011 calculamos los terminos switching t1 y t2
c
      do i=1,3
c       t1(i)=0.5d0*(1-tanh(w1*(rc1h(i)-w2)))
c       t1(i)=(1-tanh(w1*(rc1h(i)-w2)))
        t1(i)=(1-tanh(w1*(rc1h(i)-r0c1h)))
      enddo
c       t2(1)=0.5d0*(1-tanh(w7*(rc2h(1)-w8)))
c       t2(1)=(1-tanh(w7*(rc2h(1)-w8)))
        t2=(1-tanh(w7*(rc2h(1)-r0c2h)))
c
c     write(*,*) t1(1)*t2(1),t1(2)*t2(1),t1(3)*t2(1)
c
       vtor=0.0d0
c 
c CRD 2021 incluimos terminos para pasar de 3 a 2 minimos
c
c     arg1=(t1(1)*t1(2)*t1(3)*t2)**3.0d0
      arg1=(rb1h(1)-r0h1h)*(rb1h(2)-r0h1h)*(rb1h(3)-r0h1h)*
     *(rb2h(1)-r0h2h)
      arg2=tanh(arg1)
      if (arg2.lt.0.0) then
      A=2.0d0
      else
      A=2.0d0+1.d0*arg2
      endif
c     write(45,*) arg1, arg2, A
      A=3.00d0
c
c      vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(1)))*t1(1)*t2
c      vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(2)))*t1(2)*t2
c      vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(3)))*t1(3)*t2
       vtor=vtor+(V3/3.0d0)*(1+cos(A*thetahcch(1)))*t1(1)*t2
       vtor=vtor+(V3/3.0d0)*(1+cos(A*thetahcch(2)))*t1(2)*t2
       vtor=vtor+(V3/3.0d0)*(1+cos(A*thetahcch(3)))*t1(3)*t2
       return
       end
c
c------------------------------------------------------
c      double precision function dihed(i,j,k,l)
       subroutine dihedro(i,j,k,l, dihed)
c------------------------------------------------------
c
         implicit double precision (a-h,o-z)
c        include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
         integer i, j, k,l
         integer m
c        double precision x(300)
c        double precision dista
c        double precision d12,d13,d14,d23,d24,d34
c        double precision p,q
c        double precision p
c        double precision pi
c
         double precision anorm
         double precision dotprod
         double precision rij(3),rjk(3),rkj(3),rkl(3)
         double precision zijjk(3),zijkj(3), zkjkl(3), zjkkl(3)
         double precision tmpsign(3), tmp, arg1,arg2
c
c        pi=4.0d0*atan(1.0d0)
c
         do m = 1, 3
            rij(m)=q((i-1)*3+m)-q((j-1)*3+m)
            rjk(m)=q((j-1)*3+m)-q((k-1)*3+m)
            rkl(m)=q((k-1)*3+m)-q((l-1)*3+m)
         enddo
         call crosprod(rjk,rkl,zjkkl)
         call crosprod(rij,rjk,zijjk)
         arg2 = anorm(rjk)*dotprod(rij,zjkkl)
         arg1 = dotprod(zijjk,zjkkl)
c
c
         if (abs(arg1).lt.1.d-6) then
             if (arg2.gt.0.d0) dihed = 90.d0
             if (arg2.lt.0.d0) dihed = -90.d0
             return
         endif
         if (arg1.gt.0.d0) then
                 dihed = atan(arg2/arg1)
     *           *180.d0/pi
                 return
         endif
         if (arg2.ge.0.d0) then
                 dihed = 180.d0 + atan(arg2/arg1)
     *           *180.d0/pi
                 return
         endif
         if (arg2.lt.0.d0) then
                 dihed = -180.d0 + atan(arg2/arg1)
     *           *180.d0/pi
                 return
         endif
c
        dihed=dihed*pi/180.d0
c
         return
         end
c
c------------------------------------------------------
       subroutine crosprod(x,y,z)
c------------------------------------------------------
c
         implicit none
         double precision x(3),y(3),z(3)
c
         z(1)=x(2)*y(3)-x(3)*y(2)
         z(2)=x(3)*y(1)-x(1)*y(3)
         z(3)=x(1)*y(2)-x(2)*y(1)
c
         return
         end
c
c------------------------------------------------------
         double precision function dotprod(x,y)
c------------------------------------------------------
c
         implicit none
         integer k
         double precision x(3),y(3), tmp
c
         tmp = 0.d0
         do k = 1, 3
            tmp = tmp + x(k)*y(k)
         enddo
         dotprod = tmp
c
         return
         end
c
c------------------------------------------------------
         double precision function anorm(x)
c------------------------------------------------------
         implicit none
         integer k
         double precision x(3), tmp
c
         tmp = 0.d0
         do k = 1, 3
            tmp = tmp + x(k)*x(k)
         enddo
         anorm = dsqrt(tmp)
c
         return
         end
c
c------------------------------------------------------
       subroutine refangles
c------------------------------------------------------
c
c  subroutine calculates reference angles for the "in-plane" potential
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
c      dimension sumd2(4),sumd4(4),ddr(4,4)
       tau=acos(-1.0d0/3.0d0)
c
c CRD 2012 taucc para el angulo H-C-C
c
       tauco=acos(-0.3616d0)
       tauoc=acos(-0.4558d0)
c
c      pi=4.0d0*atan(1.0d0)
c      halfpi=0.5d0*pi
       halfpi=0.6555d0*pi
       twopi=2.0d0*pi
c jcc-2010
c valor para ajustar el angulo
c      anghch=0.649d0*pi
c      tausih = pi - asin ( sin(anghch/2.d0) / sin(pi/3.d0) )
c      tausih=acos(-0.4648d0)
c jcc-2010
c
c
c  set diagonal elements to zero
c
       do i=1,6
         theta0(i,i)=0.0d0
       enddo
c
c  calculate reference angles
c
c el atomo 4 es el c2
c
       theta0(1,2)=tau+(tau-halfpi)*(sphi(1)*sphi(2)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
       theta0(1,3)=tau+(tau-halfpi)*(sphi(1)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
       theta0(1,4)=tauco+(tau-halfpi)*(sphi(1)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
       theta0(2,3)=tau+(tau-halfpi)*(sphi(2)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
       theta0(2,4)=tauco+(tau-halfpi)*(sphi(2)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
       theta0(3,4)=tauco+(tau-halfpi)*(sphi(3)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
c
c  el atomo 6 es el c1
c
       theta0(5,6)=tauoc+(tauoc-halfpi)*(sphi(5)*sphi(6)-1.0d0)
c
c  fill in the other half of the matrix
c
        do i=1,3
          do j=i+1,4
            theta0(j,i)=theta0(i,j)
          enddo
        enddo
            theta0(6,5)=theta0(5,6)
c
c      do i=1,3
c         do j=i+1,4
c          write(44,*) 'theta',i,j,'=', theta0(i,j)*180/pi
c         enddo
c      enddo
c          write(44,*) 'theta(5,6)=', theta0(5,6)*180/pi
       return
       end
c
c******************************************************
c
c
c------------------------------------------------------
       subroutine stretch(vstr)
c------------------------------------------------------
c
c  subroutine to calculate leps-type stretching potential and its
c  derivatives
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
       dimension vqb1h(3),vjb1h(3),vqb2h(1),vjb2h(1),vq(4),vj(4),
     *           vqc1h(3), vqc2h(1), vjc1h(3), vjc2h(1)
c
c  calculate avergage bond length for the methane moiety
c
       rav1=(rc1h(1)+rc1h(2)+rc1h(3))/3.0d0
       rav2=rc2h(1)
c      write(44,*) rav1, rav2
c
c  initialise:
c
       vstr=0.0d0
c
c  ach:
c
c  nb: in double precision tanh(19.0d0)=1.0d0 and we put the if statement
c  in to avoid overflow/underflow errors
c
       arga1=c1c1h*(rav1-r0c1h)
       arga2=c1c2h*(rav2-r0c2h)
       if(arga1.lt.19.0d0)then
         ach1=a1c1h+b1c1h*(tanh(arga1)+1.0d0)*0.5d0
       else
         ach1=a1c1h+b1c1h
       endif
       if(arga2.lt.19.0d0)then
         ach2=a1c2h+b1c2h*(tanh(arga2)+1.0d0)*0.5d0
       else
         ach2=a1c2h+b1c2h
       endif
c      write(44,*) 'r0c1h',r0c1h,'rav1',rav1
c      write(44,*) 'r0c2h',r0c2h, 'rav2',rav2
c      write(44,*) 'ach1',ach1,'ach2', ach2
c
c  calculate singlet: e1, triplet: e3 energies and vq and vj
c  terms for each bond
c
c      write(44,*) d1c1b,d3c2b,ac1b,rc1b,r0c1b
c      write(44,*) d1c2b,d3c2b,ac2b,rc2b,r0c2b
       e1=d1c1b*(exp(-2.0d0*ac1b*(rc1b-r0c1b))
     *              -2.0d0*exp(-ac1b*(rc1b-r0c1b)))
       e3=d3c1b*(exp(-2.0d0*ac1b*(rc1b-r0c1b))
     *              +2.0d0*exp(-ac1b*(rc1b-r0c1b)))
       vqc1b=(e1+e3)*0.5d0
       vjc1b=(e1-e3)*0.5d0
c      write(44,*) vqc1b, vjc1b
c
       e1=d1c2b*(exp(-2.0d0*ac2b*(rc2b-r0c2b))
     *               -2.0d0*exp(-ac2b*(rc2b-r0c2b)))
       e3=d3c2b*(exp(-2.0d0*ac2b*(rc2b-r0c2b))
     *               +2.0d0*exp(-ac2b*(rc2b-r0c2b)))
       vqc2b=(e1+e3)*0.5d0
       vjc2b=(e1-e3)*0.5d0
c     write(44,*) vqc2b, vjc2b
       do i=1,3
c      write(44,*) 'ach1',ach1,'rc1h',rc1h(i),r0c1h
         e1=d1c1h*(exp(-2.0d0*ach1*(rc1h(i)-r0c1h))
     *              -2.0d0*exp(-ach1*(rc1h(i)-r0c1h)))
         e3=d3c1h*(exp(-2.0d0*ach1*(rc1h(i)-r0c1h))
     *              +2.0d0*exp(-ach1*(rc1h(i)-r0c1h)))
         vqc1h(i)=(e1+e3)*0.5d0
         vjc1h(i)=(e1-e3)*0.5d0
c      write(44,*) 'vqc1h',i, vqc1h(i), 'vjc1h',i, vjc1h(i)
       enddo
c
         e1=d1c2h*(exp(-2.0d0*ach2*(rc2h(1)-r0c2h))
     *              -2.0d0*exp(-ach2*(rc2h(1)-r0c2h)))
         e3=d3c2h*(exp(-2.0d0*ach2*(rc2h(1)-r0c2h))
     *              +2.0d0*exp(-ach2*(rc2h(1)-r0c2h)))
         vqc2h(1)=(e1+e3)*0.5d0
         vjc2h(1)=(e1-e3)*0.5d0
c      write(44,*) 'vqc2h', vqc2h(1), 'vjc2h',vjc2h(1)
c
       do i=1,3
         e1=d1h1h*(exp(-2.0d0*ah1h*(rb1h(i)-r0h1h))
     *              -2.0d0*exp(-ah1h*(rb1h(i)-r0h1h)))
         e3=d3h1h*(exp(-2.0d0*ah1h*(rb1h(i)-r0h1h))
     *              +2.0d0*exp(-ah1h*(rb1h(i)-r0h1h)))
         vqb1h(i)=(e1+e3)*0.5d0
         vjb1h(i)=(e1-e3)*0.5d0
c      write(44,*) 'vqb1h',i, vqb1h(i), 'vjb1h',i, vjb1h(i)
       enddo
         e1=d1h2h*(exp(-2.0d0*ah2h*(rb2h(1)-r0h2h))
     *              -2.0d0*exp(-ah2h*(rb2h(1)-r0h2h)))
         e3=d3h2h*(exp(-2.0d0*ah2h*(rb2h(1)-r0h2h))
     *              +2.0d0*exp(-ah2h*(rb2h(1)-r0h2h)))
         vqb2h(1)=(e1+e3)*0.5d0
         vjb2h(1)=(e1-e3)*0.5d0
c      write(44,*) 'vqb2h', vqb2h(1), 'vji2h',vjb2h(1)
c
c  calculate 3 body potential
c
       do i=1,3
c       write(44,*) i,'vqc1h',vqc1h(i),'vqc1b',vqc1b, 'vqb1h',vqb1h(i)
         vq(i)=vqc1h(i)+vqc1b+vqb1h(i)
         vj(i)=-sqrt(((vjc1h(i)-vjc1b)**2+(vjc1b-vjb1h(i))**2
     *                 +(vjb1h(i)-vjc1h(i))**2)*0.5d0)
         vstr=vstr+vq(i)+vj(i)
       enddo
c
c
         vq(4)=vqc2h(1)+vqc2b+vqb2h(1)
         vj(4)=-sqrt(((vjc2h(1)-vjc2b)**2+(vjc2b-vjb2h(1))**2
     *                 +(vjb2h(1)-vjc2h(1))**2)*0.5d0)
         vstr=vstr+vq(4)+vj(4)
c
c CRD 2012 potencial morse para el enlace C-C
c
       dt=(rcc-r0cc)
       expterm=exp(-acc*dt)
       vcc=d1cc*(1.0d0-expterm)**2.0d0
       vstr=vstr + vcc
c      write(44,*) 'rcc',rcc,'r0cc',r0cc,'acc',acc
c      write(44,*) 'vcc',vcc
c      write(44,*) 'vstr',vstr
c
       return
       end
c
c
c------------------------------------------------------
       subroutine opbend(vop)
c------------------------------------------------------
c
c  subroutine calculates symmetrized vop potential and derivatives
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
       double precision norma
       dimension sumd2(4),sumd4(4)
       dimension in(3),a(3),b(3),axb(3),c(4,3),argd(4)
c      dimension sumd2(8),sumd4(8)
c      dimension in(3),a(3),b(3),axb(3),c(8,3),argd(8)
c
c
       vop=0.0d0
c
c  calculate force constants
c
       call opforce
c
c  calculate out-of-plane angle
c
c CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C2
c
       do i=1,4
         j=i+1
         if(j.gt.4)j=j-4
         k=j+1
         if(k.gt.4)k=k-4
         l=k+1
         if(l.gt.4)l=l-4
c
         call calcdelta1(i,j,k,l,sum2,sum4)
         sumd2(i)=sum2
         sumd4(i)=sum4
         vop=vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
c      write(44,*) i,j,k,l,sumd2(i)*fdelta(i),sumd4(i)*hdelta(i),vop
       enddo
c
       return
       end
c
c------------------------------------------------------
       subroutine calcdelta1(i,j,k,l,sum2,sum4)
c------------------------------------------------------
c
c  subroutine calculates out of plane angle delta, loops
c  through delta(i,j), delta(i,k), delta(i,l)
c
c   also calculates the derivatives wrt delta
c
       implicit double precision (a-h,o-z)
       double precision norma
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
       dimension  delta(8),in(3),a(3),b(3),axb(3),c(8,3),argd(8),
     *            daxb(8,3,3),cdot(8,3,3),atemp2(3)
c
c  initialise
c
       sum2=0.0d0
       sum4=0.0d0
c
c  set j,k,l indices
c
       in(1)=j
       in(2)=k
       in(3)=l
c
c CRD 2005 incluimos bucles if para asegurar que siempre
c que una de los valores j, k o l sea 4 las coordenadas
c corresponden al atomo nuevo.
c
c  vector a is rk-rj, vector b is rl-rj
c
       if(j.eq.4)then
        do ind=1,3
          a(ind)=q(n1h(k,ind))-q(nc2(ind))
          b(ind)=q(n1h(l,ind))-q(nc2(ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
       else
       endif
c
       if(k.eq.4)then
        do ind=1,3
          a(ind)=q(nc2(ind))-q(n1h(j,ind))
          b(ind)=q(n1h(l,ind))-q(n1h(j,ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
       else
       endif
c
       if(l.eq.4)then
        do ind=1,3
          a(ind)=q(n1h(k,ind))-q(n1h(j,ind))
          b(ind)=q(nc2(ind))-q(n1h(j,ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
c
c modificacion para que hiciese vxu
c      axb(1)=b(2)*a(3)-b(3)*a(2)
c      axb(2)=a(1)*b(3)-b(1)*a(3)
c      axb(3)=b(1)*a(2)-a(1)*b(2)
c
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
       else
       endif
c
c CRD 2012 eliminamos porque seria la salida del C2
c
       if(k.lt.4.and.j.lt.4.and.l.lt.4)then
        do ind=1,3
          a(ind)=q(n1h(k,ind))-q(n1h(j,ind))
          b(ind)=q(n1h(l,ind))-q(n1h(j,ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
       else
       endif
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
       if(in(ii).eq.4) then
         do ind=1,3
           c(in(ii),ind)=-tcc(ind)/rcc
         enddo
       else
         do ind=1,3
           c(in(ii),ind)=-tc1h(in(ii),ind)/rc1h(in(ii))
         enddo
       endif
       enddo
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
         delta(in(ii))=acos(argd(in(ii)))-theta0(i,in(ii))
         sum2=sum2+delta(in(ii))**2
         sum4=sum4+delta(in(ii))**4
       enddo
c
       return
       end
c
c------------------------------------------------------
       subroutine opforce
c------------------------------------------------------
c
c  calculates the out-of-plane bending force constants
c  and their derivatives
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
       dimension switch(4)
c      dimension switch(8)
c
c  calculate switching functions:
c
       switch(1)=(1.0d0-s3(1))*s3(2)*s3(3)*s3(4)
       switch(2)=(1.0d0-s3(2))*s3(3)*s3(4)*s3(1)
       switch(3)=(1.0d0-s3(3))*s3(4)*s3(1)*s3(2)
       switch(4)=(1.0d0-s3(4))*s3(1)*s3(2)*s3(3)
c      switch(5)=(1.0d0-s3(5))*s3(6)*s3(7)*s3(8)
c      switch(6)=(1.0d0-s3(6))*s3(7)*s3(8)*s3(5)
c      switch(7)=(1.0d0-s3(7))*s3(8)*s3(5)*s3(6)
c      switch(8)=(1.0d0-s3(8))*s3(5)*s3(6)*s3(7)
c      do i=1,4
c       write(44,*) 'switch',i,'=',switch(i)
c      enddo
c       write(*,*) 'pasa'
c
c  calculate the force constants 
c
c      do i=1,8
       do i=1,4
         fdelta(i)=switch(i)*fch3
         hdelta(i)=switch(i)*hch3
       enddo
c
       return
       end
c
c
c------------------------------------------------------
       subroutine ipbend(vip)
c------------------------------------------------------
c
c  subroutine calculates symmetrised in plane bend term
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
       dimension costh(6,6),theta(6,6),dth(6,6)
c
c  initialise
c
       vip=0.0d0
c
c  calculate force constants: fk0(i,j), f1(i)
c  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
c
       call ipforce
c
c  calculate theta(i,j) and in plane bend potential
c
c c1
c
       do i=1,3
         do j=i+1,4
           if(j.lt.4.0) then
           costh(i,j)=tc1h(i,1)*tc1h(j,1)+tc1h(i,2)*tc1h(j,2)
     *                       +tc1h(i,3)*tc1h(j,3)
           costh(i,j)=costh(i,j)/rc1h(i)/rc1h(j)
           theta(i,j)=acos(costh(i,j))
           dth(i,j)=theta(i,j)-theta0(i,j)
c      write(44,*) 'f1 i',i,f1(i),'f1 j',j,f1(j)
           vip=vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
         else
           costh(i,j)=tc1h(i,1)*tcc(1)+tc1h(i,2)*tcc(2)
     *                       +tc1h(i,3)*tcc(3)
           costh(i,j)=costh(i,j)/rc1h(i)/rcc
           theta(i,j)=acos(costh(i,j))
           dth(i,j)=theta(i,j)-theta0(i,j)
           vip=vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
c     write(44,*)'fk0i,j',fk0(i,j),'f1i',f1(i),'f1j',f1(j),'dt',dth(i,j)
         endif
       enddo
       enddo
c     write(44,*) 'vip do 1', vip
c
c c2
c
           costh(5,6)=tc2h(1,1)*(-tcc(1))+tc2h(1,2)*(-tcc(2))
     *                       +tc2h(1,3)*(-tcc(3))
           costh(5,6)=costh(5,6)/rc2h(1)/rcc
           theta(5,6)=acos(costh(5,6))
           dth(5,6)=theta(5,6)-theta0(5,6)
           vip=vip+0.5d0*fk0(5,6)*f1(5)*f1(6)*dth(5,6)**2
c     write(44,*) 'theta 5,6', theta(5,6), theta0(5,6)
c     write(44,*)'fk05,6',fk0(5,6),'f15',f1(5),'f16',f1(6),'dt',dth(5,6)
c     write(44,*) 'vip do 2', vip
c
       return
       end
c
c------------------------------------------------------
       subroutine ipforce
c------------------------------------------------------
c
c  calculates the symmetrised in plane bend force constants and
c  all partial derivatives involving them
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
c  set force constant at asymptotes
c
       f0=fkinf1+ak1
       f2=fkinf1
c
       fk0(1,2)=f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
       fk0(1,3)=f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
       fk0(1,4)=f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
       fk0(2,3)=f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
       fk0(2,4)=f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
       fk0(3,4)=f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)
c
c      fk0(5,6)=fkinf2*(s1(5)*s1(6))
c      fk0(5,6)=fkinf2*s1(5)
       fk0(5,6)=f0+f0*(s1(5)*s1(6)-1.0d0)+(f0-f2)*(s2(5)*s2(6)-1.0d0)
c
c f1 para los h del c1
c
       do i=1,3
c
c  calc derivatives of fk0 wrt each of the rch(i) bonds
c  calculate the terms f1(i)
c
c       write(44,*) 'rb1h',i, rb1h(i)
         arga1=aa11*rb1h(i)*rb1h(i)
         arga2=aa41*(rb1h(i)-r0hh)*(rb1h(i)-r0hh)
         a1=1.0d0-exp(-arga1)
         a2=aa21+aa31*exp(-arga2)
c       write(44,*) 'a1',i, a1,'a2',i,a2,'aa11',aa11
         f1(i)=a1*exp(-a2*(rc1h(i)-r0c1h)**2)
       enddo
c
c f1 para el c2
c
c CRD 2012 eliminamos el calculo de a1c y a2c pq es constante
c lo incluimos en el CONST
c
         f1(4)=a1c*exp(-a2c*(rcc-r0cc)**2)
         f1(6)=f1(4)
c
c f1 para los h del c2
c
c
c  calc derivatives of fk0 wrt each of the rch(i) bonds
c  calculate the terms f1(i)
c
         arga1=aa12*rb2h(1)*rb2h(1)
         arga2=aa42*(rb2h(1)-r0hh)*(rb2h(1)-r0hh)
         a1=1.0d0-exp(-arga1)
         a2=aa22+aa32*exp(-arga2)
         f1(5)=a1*exp(-a2*(rc2h(1)-r0c2h)**2)
c
       return
       end
c
c******************************************************
c
c------------------------------------------------------
       subroutine switchf
c------------------------------------------------------
c
c  calculates switching functions: s3,sphi,stheta
c
       implicit double precision (a-h,o-z)
c      include 'ch3oh_h.inc'
CC INCLUDE DEL VENUS
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
c
c  nb remember that integration units are:
c  energy in   1.0d+05 j/mol
c  time in     1.0d-14 s
c
c original     a1s=1.5313681d-7
c      a1s=1.5313681d-8
c      b1s=-4.6696246d0
c original      a2s=1.0147402d-7
c      a2s=1.0147402d-8
c      b2s=-12.363798d0
c
c  use double precision criterion:
c
c  tanh(19.0d0)=1.0d0
c
       argmax=19.0d0
c 
c  calculate s1 
c
       do i=1,3
         args1=a1s1*(rc1h(i)-r0c1h)*(rc1h(i)-b1s1)**8
         if(args1.lt.argmax)then
           s1(i)=1.0d0-tanh(args1)
         else
           s1(i)=0.0d0
         endif
       enddo
c
         args1=a2s1*(rc2h(1)-r0c2h)*(rc2h(1)-b2s1)**8
         if(args1.lt.argmax)then
           s1(5)=1.0d0-tanh(args1)
         else
           s1(5)=0.0d0
         endif
c 
c CRD 2012 todos los S de los carbonos son 1
c
         argsc1=a1scc*(rcc-r0cc)*(rcc-b1scc)**8
         if(argsc1.lt.argmax)then
           s1(4)=1.0d0-tanh(argsc1)
           s1(6)=s1(4)
         else
           s1(4)=0.0d0
           s1(6)=s1(4)
         endif
c
c CRD correccion para dejar C-C fijo
c
           s1(4)=1.0d0
           s1(6)=1.0d0
c
c  calculate s2 and ds2
c
       do i=1,3
         args2=a1s2*(rc1h(i)-r0c1h)*(rc1h(i)-b1s2)**6
         if(args2.lt.argmax)then
           s2(i)=1.0d0-tanh(args2)
         else
           s2(i)=0.0d0
         endif
       enddo
c
         args2=a2s2*(rc2h(1)-r0c2h)*(rc2h(1)-b2s2)**6
         if(args2.lt.argmax)then
           s2(5)=1.0d0-tanh(args2)
         else
           s2(5)=0.0d0
         endif
c 
         argsc2=a2scc*(rcc-r0cc)*(rcc-b2scc)**6
         if(argsc2.lt.argmax)then
           s2(4)=1.0d0-tanh(argsc2)
           s2(6)=s2(4)
         else
           s2(4)=0.0d0
           s2(6)=s2(4)
         endif
c
c CRD correccion para dejar C-C fijo
c
           s2(4)=1.0d0
           s2(6)=1.0d0
c
c  calculate s3 and ds3
c
       do i=1,3
         args3=a3s*(rc1h(i)-r0c1h)*(rc1h(i)-b3s)**2
         if (args3.lt.argmax)then
           s3(i)=1.0d0-tanh(args3)
         else
           s3(i)=0.0d0
         endif
       enddo
c
         argsc3=a3scc*(rcc-r0cc)*(rcc-b3scc)**2
         if (argsc3.lt.argmax)then
           s3(4)=1.0d0-tanh(argsc3)
         else
           s3(4)=0.0d0
         endif
c
c CRD correccion para dejar C-C fijo
c
           s3(4)=1.0d0
c
c  calculate sphi 
c
c  condition here is on the bondlength rch(i)
c  st argsphi is lt approx 19.0d0
c
       do i=1,3
         if(rc1h(i).lt.3.8d0)then
           argsphi=aphi1*(rc1h(i)-r0c1h)*exp(bphi1*(rc1h(i)-cphi1)**3)
           sphi(i)=1.0d0-tanh(argsphi)
         else
           sphi(i)=0.0d0
         endif
       enddo
c
         if(rc2h(1).lt.3.8d0)then
           argsphi=aphi2*(rc2h(1)-r0c2h)*exp(bphi2*(rc2h(1)-cphi2)**3)
           sphi(5)=1.0d0-tanh(argsphi)
         else
           sphi(5)=0.0d0
         endif
c
         if(rcc.lt.3.8d0)then
           argsphic=aphi*(rcc-r0cc)*exp(bphi*(rcc-cphi)**3)
           sphi(4)=1.0d0-tanh(argsphic)
           sphi(6)=sphi(4)
         else
           sphi(4)=0.0d0
           sphi(6)=sphi(4)
         endif
c
c          sphi(4)=1.0d0
c          sphi(6)=1.0d0
c
c      do i=1,6
c          write(44,*) 'sphi',i,sphi(i)
c      enddo
c
c  calculate stheta and dstheta
c
       do i=1,3
         if(rc1h(i).lt.3.8d0)then
       argstheta=atheta*(rc1h(i)-r0c1h)*exp(btheta*(rc1h(i)-ctheta)**3)
           stheta(i)=1.0d0-tanh(argstheta)
         else
           stheta(i)=0.0d0
         endif
       enddo
c
c      do i=1,3
c        if(rc2h(i).lt.3.8d0)then
c      argstheta=atheta*(rc2h(i)-r0c2h)*exp(btheta*(rc2h(i)-ctheta)**3)
c          stheta(i+4)=1.0d0-tanh(argstheta)
c        else
c          stheta(i+4)=0.0d0
c        endif
c      enddo
c
c CRD 2012 todos los S de los carbonos son 1
c
         if(rcc.lt.3.8d0)then
          argsthetac=athetacc*(rcc-r0cc)*exp(bthetacc*(rcc-cthetacc)**3)
          stheta(4)=1.0d0-tanh(argsthetac)
          stheta(6)=stheta(4)
         else
          stheta(4)=0.0d0
          stheta(6)=stheta(4)
         endif
c
c          stheta(4)=1.0d0
c          stheta(6)=1.0d0
c
       return
       end
C
c
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
       implicit double precision (a-h,o-z)
CC INCLUDE DEL VENUS
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
      PARAMETER (N3TMMN = 27)
C
c     COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
c     COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
c     COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
c     COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
c    +               PENGYIJ(JSURF),
c    +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
c    +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
cc call de venus
C
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
c      REF(1)='J. Espinosa-Garcia and J. C. Corchado'
c      REF(2)='J. Chem. Phys., Vol. 112, p. 5731, 2000'
C
      INDEXES(1) = 6
      INDEXES(2) = 8
      INDEXES(3) = 1
      INDEXES(4) = 1
      INDEXES(5) = 1
      INDEXES(6) = 1
      INDEXES(7) = 35
C
c     IRCTNT=9
      IRCTNT=7
C
C      CALL POTINFO
C
      CALL ANCVRT
c
c  calculate indexes for coordinates
c
       do ind=1,3
         icount=ind-3
         nc1(ind)=3*nnc1+icount
         nc2(ind)=3*nnc2+icount
         nhb(ind)=3*nnb+icount
c
         do i=1,3
           n1h(i,ind)=3*nc1h(i)+icount
         enddo
           n2h(1,ind)=3*nc2h(1)+icount
       enddo
c
c
c  convert to appropriate units:
c
c  From:
c
c  kcal/mol
c  mdyn A-1        -> 1.0d+05 j/mol...
c  mdyn A rad-1
c
c  To
c
c  energy   in 1.0d+05 j/mol
c  distance in 1.0d-10 m
c  angles   in radians
c  mass     in amu
c
       fact1=0.041840d0
       fact2=6.022045d0
c
       d1c1h=d1c1h*fact1
       d3c1h=d3c1h*fact1
c
       d1c2h=d1c2h*fact1
       d3c2h=d3c2h*fact1
c
       d1cc=d1cc*fact1
c
       d1c1b=d1c1b*fact1
       d3c1b=d3c1b*fact1
c
       d1c2b=d1c2b*fact1
       d3c2b=d3c2b*fact1
c
       d1h1h=d1h1h*fact1
       d3h1h=d3h1h*fact1
c
       d1h2h=d1h2h*fact1
       d3h2h=d3h2h*fact1
c
       fkinf1=fkinf1*fact2
       ak1=ak1*fact2
       fkinf2=fkinf2*fact2
c
       a1s1=a1s1*1.D-8
       a2s1=a2s1*1.D-8
c
       a1s2=a1s2*1.D-8
       a2s2=a2s2*1.D-8
c
       a1scc=a1scc*1.D-8
       a2scc=a2scc*1.D-8
c
      RETURN
      END
C
      BLOCK DATA PTPACM
C
CC INCLUDE DEL VENUS
      implicit double precision (a-h,o-z)
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
c  include file for CH3OH+H PES
c
c  CONST
c  
c     COMMON /POTCM/ nnc1,nnc2,nnb,nc1h(3),nc2h(1),
      COMMON /POTCM/ nnc1,nnc2,nnb,
     +               r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h,
     +               a1c1h,b1c1h,c1c1h,
     +               r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h,
     +               a1c2h,b1c2h,c1c2h,
     +               r0ccr,r0ccp,w9,w10,d1cc,acc,
     +               a1c,a2c,
     +               r0h1h,d1h1h,d3h1h,ah1h,
     +               r0h2h,d1h2h,d3h2h,ah2h,
     +               r0c1b,d1c1b,d3c1b,ac1b,
     +               r0c2b,d1c2b,d3c2b,ac2b,
     +               aphi,bphi,cphi,
     +               aphi1,bphi1,cphi1,
     +               aphi2,bphi2,cphi2,
     +               a3s,b3s,a3scc,b3scc,  
     +               atheta,btheta,ctheta,
     +               athetacc,bthetacc,cthetacc,
     +               fch3,hch3,
     +               fkinf1,ak1,aa11,aa21,aa31,aa41,
     +               fkinf2,aa12,aa22,aa32,aa42,
     +               V3,w1,w2,w7,w8,
     +               a1s1,a2s1,b1s1,b2s1,
     +               a1s2,a2s2,b1s2,b2s2,
     +               a1scc,a2scc,b1scc,b2scc
c
c Subroutine Coorden
c
      COMMON /POTCM2/ nc1h(3),nc2h(1)
       common /qpdot/ q(21),pdot(21)
       common /coord2/ nc1(3),nc2(1),nc(3),nhb(3),n1h(3,3),n2h(1,3)
       common /coord3/ tc1b(3),tc2b(3),tcc(3)
       common /coord4/ tc1h(3,3),tc2h(1,3),tb1h(3,3),tb2h(1,3)
       common /coord5/ rc1b,rc2b,rcc,rc1h(3),rc2h(1),rb1h(3),rb2h(1)
       common /coord6/ r0c1h,r0c2h,r0cc
c
c Subroutine torsion
c
       common /torsion1/ thetahcch(3),t2
c
c Subroutine refangles
c
       common /refangles1/ theta0(6,6),sphi(6),stheta(6)
c
c Subroutine vop
c
       common /vop6/ fdelta(4),hdelta(4)
c
c Subroutine ipforce
c
       common /ipforce1/ fk0(6,6),f1(6)
c
c Subroutine Switch
c
       common /switch1/ s1(6),s2(6),s3(4)
C
      DATA NASURF /1,35*0/
      DATA NDER /1/
      DATA NFLAG /1,1,15*0,6,0,0/
C
      DATA ANUZERO /0.0D0/
      DATA ICARTR,MSURF,MDER/1,0,1/
      DATA NULBL /25*0/
      DATA NATOMS /7/
C
C CONST.10
c
      data nnc1 /1/
      data nnc2 /2/
      data nnb /7/
      data nc1h /3,4,5/
      data nc2h /6/
      data r0c1hr /1.09047d0/
      data r0c1hp /1.08177d0/
      data w3 /1.0000d0/
      data w4 /1.09047d0/
      data d1c1h /105.23000d0/
      data d3c1h /27.30400d0/
      data a1c1h /1.59300d0/
      data b1c1h /0.52500d0/
      data c1c1h /0.51404d0/
      data r0c2hr /0.95897d0/
      data r0c2hp /0.96077d0/
      data w5 /1.000d0/
      data w6 /0.96897d0/
      data d1c2h /116.23000d0/
      data d3c2h /37.00400d0/
      data a1c2h /1.55300d0/
      data b1c2h /1.42500d0/
      data c1c2h /0.51404d0/
      data r0ccr /1.42007d0/
      data r0ccp /1.37500d0/
      data w9 /1.5000d0/
      data w10 /1.08507d0/
      data d1cc /109.94900d0/
      data acc /1.934d0/
      data a1c /1.20997d0/
      data a2c /1.949d0/
      data r0h1h /1.29191d0/
      data d1h1h /111.458d0/
      data d3h1h /37.064d0/
      data ah1h /1.9157d0/
      data r0h2h /1.29191d0/
      data d1h2h /109.458d0/
      data d3h2h /26.064d0/
      data ah2h /1.8057d0/
      data r0c1b /2.05897d0/
      data d1c1b /35.009d0/
      data d3c1b /23.000d0/
      data ac1b /3.1090285d0/
      data r0c2b /1.74897d0/
      data d1c2b /68.609d0/
      data d3c2b /17.400d0/
      data ac2b /1.0170285d0/
      data aphi /0.528790d0/
      data bphi /  1.200664d0/
      data cphi /  1.420074d0/
      data aphi1 /   0.19700d0/
      data bphi1 /   1.20600d0/
      data cphi1 /   1.08200d0/
      data aphi2 /   0.19700d0/
      data bphi2 /   1.20600d0/
      data cphi2 /   0.96800d0/
      data a3s /   0.12770d0/
      data b3s /   1.08200d0/
      data a3scc /   0.297000d0/
      data b3scc /  1.420664d0/
      data atheta /   1.65781d0/
      data btheta /   0.950600d0/
      data ctheta /   1.08201d0/
      data athetacc /   1.45781d0/
      data bthetacc /   0.950600d0/
      data cthetacc /   1.42001d0/
      data fch3 /   0.19770d0/
      data hch3 /   0.09600d0/
      data fkinf1 /   0.44770d0/
      data ak1 /   0.12600d0/
      data aa11 /   0.15501d0/
      data aa21 /   1.59995d0/
      data aa31 /   2.16996d0/
      data aa41 /  11.56595d0/
      data fkinf2 /   0.45770d0/
      data aa12 /   0.05001d0/
      data aa22 /   1.30095d0/
      data aa32 /   1.49996d0/
      data aa42 /  11.46595d0/
      data V3 /   0.03000d0/
      data w1 /   0.11000d0/
      data w2 /   1.08697d0/
      data w7 /   0.11000d0/
      data w8 /   1.08367d0/
      data a1s1 /   1.53137d0/
      data b1s1 /   1.09463d0/
      data a2s1 /   1.01474d0/
      data b2s1 /   0.96898d0/
      data a1s2 /   1.53137d0/
      data b1s2 /   1.09463d0/
      data a2s2 /   1.51474d0/
      data b2s2 /   0.96898d0/
      data a1scc /   1.53137d0/
      data b1scc /   1.42003d0/
      data a2scc /   1.51474d0/
      data b2scc /   1.42008d0/
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
