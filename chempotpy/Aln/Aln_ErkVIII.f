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
c   Common name:                   ErkVIII
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
c  References: S. Erkoc, Docga(A1), 9, 203 (1985).
c              S. Erkoc, Empirical Potential Energy Functions
c                 Used in the Simulations of Materials Properties,
c                 Annual Reviews of Computational Physics IX,
c                 D. Stauffer, ed., World Scientific, Singapore 2001,
c                 pp. 1-103.
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
      dimension x(maxatom),y(maxatom),z(maxatom),r(maxatom,maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)

c  calculate and tabulate the distances between all atoms
      do 1 i=1,natom
      r(i,i)=0d0
      do 1 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r(i,j)=dsqrt(dx*dx+dy*dy+dz*dz)
    1 r(j,i)=r(i,j)

c  calculate the potential energy
      r0 = 2.510d0  /  autoang
      a1 = 1.860d0  /  autoev
      a2 = 3.410d0  /  autoev
      b1 = 1643.d0  / (autoev*autoang**9)
      b2 = 1921.d0  / (autoev*autoang**11)
      em = 11.d0
      en = 6.d0

      u=0.d0
      do 2 i=1,natom-1
      do 3 j=i+1,natom
      rr=r(i,j)
      rel=r0/rr
      uem=a1*rel**em
      uen=a2*rel**en
      u=u+(uem-uen)
    3 continue
    2 continue

      v1=0.d0
      if(natom.eq.2) go to 55
      do 4 i=1,natom-2
      do 5 j=i+1,natom-1
      do 6 m=j+1,natom
      call scos(i,j,m,x,y,z,r,natom,maxatom,cote1,cote2,cote3)
      rr1=1.d0/r(i,j)
      rr2=1.d0/r(i,m)
      rr3=1.d0/r(j,m)
      cot13=1.d0+3.d0*cote1*cote2*cote3
      war=dabs(cote1)-1.d0
      if(war.le.10d-12) go to 69
      arc1=dacos(cote1)
      go to 66
  69  arc1=0d0
  66  war=dabs(cote2)-1d0
      if(war.le.10d-12) go to 68
      arc2=dacos(cote2)
      go to 67
  68  arc2=0d0
  67  continue
      cot12=6.d0*(dsin(arc1)*dsin(arc2)+cote1*cote2)
      cot3=-25.d0*(4.d0*(cote3)**3-3.d0*cote3)
      cot2=3.d0+5.d0*(2.d0*(cote3)**2-1.d0)
      bra=9.d0*cote3+cot3+cot12*cot2
      w1=b1*cot13*(rr1*rr2*rr3)**3     ! AT term
      w2=b2*bra*(rr1**3)*(rr2*rr3)**4
      v1=v1+w1+w2
    6 continue
    5 continue
    4 continue 
   55 v=u+v1
      return
      end
c
c
          subroutine scos(i,j,m,x,y,z,r,natom,maxatom,cote1,cote2,cote3)
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
      dimension x(maxatom),y(maxatom),z(maxatom),r(maxatom,maxatom)
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
      cote1=scalar1/(r(i,j)*r(i,m))
c  calculate cosine of the angle between bonds ji and jm
      cote2=scalar2/(r(i,j)*r(j,m))
c  calculate cosine of the angle between bonds mi and mj
      cote3=scalar3/(r(i,m)*r(j,m))
      return
      end
