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

c   System:                        Aln
c   Functional form:               Lennard-Jones plus Axilrod-Teller
c   Common name:                   deSPH/LJAT
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c
c   Notes:  Many-body aluminum potential energy function.  For a recent
c           discussion of this potential, see:
c           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c           and D. G. Truhlar, "Analytic Potentials Energy Functions
c           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c
c  References: P. de Sainte Claire, G. H. Peslherbe, and W. L. Hase
c              J. Phys. Chem., 99, 8147 (1995).
c
c  Units:
c       energies    - hartrees
c       coordinates - bohrs
c
c  v       --- energy of the system (output)
c  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
c              of the system (input)
c  natom   --- number of atoms in the system (input)
c  maxatom --- maximal number of atoms (input)

      implicit real*8(a-h,o-z)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)
      parameter(autokpm=627.5095d0)
      dimension x(maxatom),y(maxatom),z(maxatom)
c
c  calculate the potential energy
c in kJ/mol
c      a = 2975343.77d0 / ((autoang**12) * autokpm)
c      b = -17765.823d0 / ((autoang**6)  * autokpm)
c      c = 81286.093d0  / ((autoang**9)  * autokpm)

c in eV
      a = 129023.16d0   / ((autoang**12) * autoev)
      b = -770.399d0    / ((autoang**6)  * autoev)
      c = 3524.899741d0 / ((autoang**9)  * autoev)
      sig = 2.3478d0  /   autoang
      aa = 4.6d0      /                  autoev

      v=0.d0

c  1. calculate energy of two-body interactions L-J
      u2=0.d0
      do 2 i=1,natom-1
      do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      u2 = u2 + a/rr**12 + b/rr**6
    3 continue
    2 continue

c  2. calculate energy of three-body interactions
      u3=0.d0
      if(natom.eq.2) go to 7
      do 6 i=1,natom-2
      do 4 j=i+1,natom-1
      do 5 m=j+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rij=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      dx=x(m)-x(j)
      dy=y(m)-y(j)
      dz=z(m)-z(j)
      rjm=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      dx=x(m)-x(i)
      dy=y(m)-y(i)
      dz=z(m)-z(i)
      rim=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      call scos(i,j,m,x,y,z,rij,rjm,rim,natom,maxatom,co1,co2,co3)
      u3=u3+c*(1.d0+3.d0*co1*co2*co3)/((rij*rjm*rim)**3)
    5 continue
    4 continue
    6 continue
c  calculate the potential energy
    7 v=u2+u3
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
