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
c   Common name:                   PeaTHT
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
c  References: E. Pearson, T. Takai, T. Halicioglu, and 
c                 W.A. Tiller, J. Cryst. Growth 70, 33 (1984). 
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
      dimension x(maxatom),y(maxatom),z(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)

c  calculate the potential energy
        aa= 4.d0*1.216d0/autoev
        r0= 2.520d0/autoang
        en= 6.d0
        em= 12.d0
        sig = (en/em)**(1.d0/(em-en))*r0
        zet= 2241.d0/(autoev*autoang**9)

      v=0d0
c     calculate energy of two-body interactions
      um=0.d0
      un=0.d0
      do 2 i=1,natom-1
      do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      rel=sig/rr
      um=um+rel**em
      un=un+rel**en
    3 continue
    2 continue
      u=aa*(um-un)

c  2. calculate energy of three-body interactions
      utr=0
          if(natom.eq.2) go to 7
          do 6 i=1,natom-2
          do 4 j=i+1,natom-1
          do 5 m=j+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r1=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      dx=x(m)-x(j)
      dy=y(m)-y(j)
      dz=z(m)-z(j)
      r2=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      dx=x(m)-x(i)
      dy=y(m)-y(i)
      dz=z(m)-z(i)
      r3=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      call scos(i,j,m,x,y,z,r1,r2,r3,natom,maxatom,cote1,cote2,cote3)
      utr=utr+zet*(1+3d0*cote1*cote2*cote3)/((r1*r2*r3)**3d0)
    5 continue
    4 continue
    6 continue
c  calculate the potential energy
    7 v=u+utr
      return
      end
c
c  
c
      subroutine scos(i,j,m,x,y,z,r1,r2,r3,natom,maxatom,
     &   cote1,cote2,cote3)
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
c  cote1   --- cosine of the angle between bonds ij and im
c  cote2   --- cosine of the angle between bonds ji and jm
c  cote3   --- cosine of the angle between bonds mi and mj
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
      cote1=scalar1/(r1*r3)
c  calculate cosine of the angle between bonds ji and jm
      cote2=scalar2/(r1*r2)
c  calculate cosine of the angle between bonds mi and mj
      cote3=scalar3/(r3*r2)
      return
      end
