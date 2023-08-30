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
C   Functional form:               Extended Rydberg plus a generalized three-body term
C   Common name:                   ER2+CoxJM
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
      dimension x(maxatom),y(maxatom),z(maxatom),c(15)
      dimension fu2(maxatom,maxatom),r(maxatom,maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)

c two-body
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

c three-body
        c0= 2.64020517830972157
        c1= -0.641182217879824101
        c2= 1.87493893502686859
        c3= -2.79115779189057145
        c4= 0.000000000000000000E+00
        c5= 0.000000000000000000E+00
        c6= 0.000000000000000000E+00
        dd= 1.00000000000000000/autoev
        ree= 2.18563751831949205/autoang
c        xa3= 0.000000000000000000E+00
        au23= 0.806253053248656593

      oneth=dsqrt(1.d0/3.d0)
      oneha=dsqrt(1.d0/2.d0)
      twoth=dsqrt(2.d0/3.d0)
      onesi=dsqrt(1.d0/6.d0)

      do 111 i=1,natom
      r(i,i)=0d0
      do 111 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r(i,j)=dsqrt(dx*dx+dy*dy+dz*dz)
  111 r(j,i)=r(i,j)

c  calculate the potential energy
      v3=0.d0
      do 2 i=1,natom
      do 3 j=i+1,natom

c  two body term

      rr=r(i,j)
      coij = 0.d0
      if (rr.ge.bu2) go to 3
      coij=dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
      rij=(rr-ree)/ree

      do 4 k=j+1,natom
c three body term

      rr=r(k,j)
      cojk = 0.d0
      if (rr.ge.bu2) go to 4
      cojk=dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
      rjk=(rr-ree)/ree

      rr=r(k,i)
      coki = 0.d0
      if (rr.ge.bu2) go to 4
      coki=dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
      rki=(rr-ree)/ree

      q1=oneth*(rij+rjk+rki)
      q2=oneha*(rjk-rki)
      q3=twoth*rij-onesi*(rjk+rki)
      q23=q2**2+q3**2
      q332=q3**3-3.d0*q3*q2**2
      poly=c0+c1*q1+c2*q1**2+c3*q23+c4*q1**3+c5*q1*q23+c6*q332
c      poly=c0+c1*q1+c2*q1**2+c3*q23
c      fsech=1.d0/dcosh(xa3*q1)
      fsech=1.d0
      co=coij*cojk*coki
      v3=v3+dd*poly*co*fsech
    4 continue
    3 continue
    2 continue

      v=v+v3

      return
      end
