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
c   Functional form:               Extended Rydberg plus generalized three-body
c   Common name:                   reCoxJM
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c
c   Notes:  Many-body aluminum potential energy function.  The parameters 
c           for this PEF have been parameterized against a data set of                             
c           aluminum cluster energies and bulk data.
c
c  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,                          
c              and D. G. Truhlar, "Analytic Potentials Energy Functions
c              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c              H. Cox, R. L. Johnston, and J. N. Murrell,
c              Surf. Sci., 373, 67 (1997).
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

c original parameters
c      a2=7.d0
c      a3=8.d0
c      dd=0.9073d0/autoev
c      re=2.7568d0/autoang
c      c0= 0.2525d0
c      c1=-0.4671d0
c      c2= 4.4903d0
c      c3=-1.1717d0
c      c4= 1.6498d0
c      c5=-5.3579d0
c      c6= 1.6327d0

c new parameters
        a2= 4.75085490962383972
        a3= 7.73864191499755716
        dd= 1.25012212994626282/autoev
        re= 2.61700048851978506/autoang
        c0= 0.182247191011235959
        c1= -0.319491939423546645
        c2= 9.46849047386419151
        c3= -0.502198339032730923
        c4= 1.59037616023448947
        c5= -4.89252564728871508
        c6= 1.19980459208597945

      oneth=dsqrt(1.d0/3.d0)
      oneha=dsqrt(1.d0/2.d0)
      twoth=dsqrt(2.d0/3.d0)
      onesi=dsqrt(1.d0/6.d0)

      do 1 i=1,natom
      r(i,i)=0d0
      do 1 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r(i,j)=dsqrt(dx*dx+dy*dy+dz*dz)
    1 r(j,i)=r(i,j)

c  calculate the potential energy
      v=0.d0
      v2=0.d0
      v3=0.d0
      do 2 i=1,natom
      do 2 j=i+1,natom
c  two body term
      rr=r(i,j)
      rij=(rr-re)/re
      v2=v2-dd*(1.d0+a2*rij)*dexp(-a2*rij)
      do 2 k=j+1,natom
c three body term
      rjk=(r(j,k)-re)/re
      rki=(r(k,i)-re)/re
      q1=oneth*(rij+rjk+rki)
      q2=oneha*(rjk-rki)
      q3=twoth*rij-onesi*(rjk+rki)
      q23=q2**2+q3**2
      q332=q3**3-3.d0*q3*q2**2
      fsech=1.d0/dcosh(a3*q1)
      poly=c0+c1*q1+c2*q1**2+c3*q23+c4*q1**3+c5*q1*q23+c6*q332
      v3=v3+dd*poly*fsech
    2 continue

      v=v2+v3

      return
      end
