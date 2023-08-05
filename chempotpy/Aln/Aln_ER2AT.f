      subroutine pes(xyz,natoms,igrad,p,g,d)
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, intent(in) :: natoms
      double precision, intent(in) :: xyz(10000,3)
      integer, intent(in) :: igrad
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: cx(natoms), cy(natoms), cz(natoms)
      double precision :: v
      integer :: istate, iatom

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
         cx(iatom)=xyz(iatom,1)/0.529177211
         cy(iatom)=xyz(iatom,2)/0.529177211
         cz(iatom)=xyz(iatom,3)/0.529177211
      enddo

      call pot(cx, cy, cz, v, natoms, 10000)

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=v*27.211386
        enddo
      else
        write (*,*) 'Only energy is available'
      endif

      endsubroutine
      subroutine pot(x,y,z,v,natom,maxatom)

C   System:                        Aln
C   Functional form:               Extended Rydberg plus Axilrod-Teller
C   Common name:                   ER2+AT
C   Number of derivatives:         0
C   Number of bodies:              variable
C   Number of electronic surfaces: 1
C   Interface:                     HO-MM-0
C
c   Notes:  Many-body aluminum potential energy function.  The parameters
c           for this PEF have been parameterized against a data set of
c           aluminum cluster energies and bulk data.  
c           This PEF has a cutoff at 6.5 angstroms.
c
c  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,                          
c              and D. G. Truhlar, "Analytic Potentials Energy Functions
c              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c
c  Units:
c       energies    - hartrees
c       coordinates - bohrs
c
c  v       --- energy of the system
c  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
c                    of the system
c  natom   --- number of atoms in the system
c  maxatom --- maximal number of atoms

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)

        b1 = 3536.39472398632142/ (autoev*autoang**9)
        b2 = 0.d0
        au23= 0.207718612603810457E-01

        de= 1.71013678553981441/autoev
        re= 5.08182706399609163
        a1= 1.24074255007327805
        a2= 0.551880801172447422
        a3= 0.129970688812896917
        au2= 0.143243771372740580
        bu2= 6.50000000000000000/autoang

c     potential, includes cutoff
      v=0.d0
      do 1 i=1,natom
      do 1 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      if (rr.le.bu2) then
      co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
      rho = rr-re
      poly=1.d0+a1*rho+a2*rho**2+a3*rho**3
      ee=dexp(-a1*rho)
      v=v-de*ee*poly*co
      endif
    1 continue


c  2. calculate energy of three-body interactions
      u3=0.d0
      if(natom.eq.2) go to 7

      do 6 i=1,natom-2
      do 4 j=i+1,natom-1
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rij=dsqrt(dx*dx+dy*dy+dz*dz)
      if (rij.ge.bu2) go to 4
      do 5 m=j+1,natom
      dx=x(m)-x(j)
      dy=y(m)-y(j)
      dz=z(m)-z(j)
      rjm=dsqrt(dx*dx+dy*dy+dz*dz)
      dx=x(m)-x(i)
      dy=y(m)-y(i)
      dz=z(m)-z(i)
      rim=dsqrt(dx*dx+dy*dy+dz*dz)
      if (rjm.ge.bu2) go to 5
      if (rim.ge.bu2) go to 5

      coij=0.d0
      cojm=0.d0
      coim=0.d0
      if (rij.lt.bu2) coij=dexp(au23*(1.d0-1.d0/(1.d0-rij/bu2)))
      if (rjm.lt.bu2) cojm=dexp(au23*(1.d0-1.d0/(1.d0-rjm/bu2)))
      if (rim.lt.bu2) coim=dexp(au23*(1.d0-1.d0/(1.d0-rim/bu2)))
      co=coij*cojm*coim

      call scos(i,j,m,x,y,z,rij,rjm,rim,natom,maxatom,co1,co2,co3)

c at term
      u3=u3+b1*(1.d0+3.d0*co1*co2*co3)/((rij*rjm*rim)**3)*co

c next term
      if ((dabs(co1)-1.d0).le.10d-12) then
        arc1=0.d0
      else
        arc1=dacos(co1)
      endif
      if ((dabs(co2)-1.d0).le.10d-12) then
        arc2=0.d0
      else
        arc2=dacos(co2)
      endif
      cot12=6.d0*(dsin(arc1)*dsin(arc2)+co1*co2)                    
      cot3=-25.d0*(4.d0*(co3)**3-3.d0*co3)                          
      cot2=3.d0+5.d0*(2.d0*(co3)**2-1.d0)                           
      bra=9.d0*co3+cot3+cot12*cot2                          
      rr1=1.d0/rij
      rr2=1.d0/rim
      rr3=1.d0/rjm
      w2=b2*bra*(rr1**3)*(rr2*rr3)**4
      u3=u3+w2

    5 continue
    4 continue
    6 continue
c  calculate the potential energy

    7 v=v+u3


      return
      end
c
c  
c
      subroutine scos(i,j,m,x,y,z,rij,rjm,rim,natom,maxatom,
     &   co1,co2,co3)
c  The subroutine calculates the cosines of the bound angles 
c  using scalar products of the position vectors of two 
c  atoms originating from the third atom.
c
c  Units:
c  coordinates - bohrs
c
c  i,j,m   --- natural indices running from 1 to natom
c  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
c              of the system
c  r       --- two-dimensional array of atomic distances
c  natom   --- number of atoms in the system
c  maxatom --- maximal number of atoms
c  co1     --- cosine of the angle between bonds ij and im
c  co2     --- cosine of the angle between bonds ji and jm
c  co3     --- cosine of the angle between bonds mi and mj
c
      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
c  calculate scalar products of the position vectors of two atoms 
c  originating from the third atom 
      scalar1=(x(j)-x(i))*(x(m)-x(i))
     *+(y(j)-y(i))*(y(m)-y(i))
     *+(z(j)-z(i))*(z(m)-z(i))
      scalar2=(x(i)-x(j))*(x(m)-x(j))
     *+(y(i)-y(j))*(y(m)-y(j))
     *+(z(i)-z(j))*(z(m)-z(j))
      scalar3=(x(i)-x(m))*(x(j)-x(m))
     *+(y(i)-y(m))*(y(j)-y(m))
     *+(z(i)-z(m))*(z(j)-z(m))
c  calculate cosine of the angle between bonds ij and im
      co1=scalar1/(rij*rim)
c  calculate cosine of the angle between bonds ji and jm
      co2=scalar2/(rij*rjm)
c  calculate cosine of the angle between bonds mi and mj
      co3=scalar3/(rim*rjm)
      return
      end

