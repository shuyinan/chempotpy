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
c   Functional form:               Embedded atom
c   Common name:                   reMisFMP
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
c              Y. Mishin, D. Farkas, M. J. Mehl, and D. A. Papaconstantopoulos
c              Mat. Res. Soc. Symp. Proc. 583, 535, (1999).
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
      dimension x(maxatom),y(maxatom),z(maxatom)

c  calculate the potential energy
        rc= 6.11700048851978462/autoang
        h = 4.78798241328773777/autoang
        e1 = 1.38861748900830473/autoev
        e2 = 0.176209086468001944E-02/autoev
        r1 = 2.34000977039570124/autoang
        r2 = 5.13873961895456777/autoang
        alp1 = 1.01670737664875421*autoang
        alp2 = 1.31289692232535415*autoang
        delta = 0.109501709819247672/autoev
        a = 0.998045920859794850E-01
        r3 = 4.41108939912066411/autoang
        r4 = -6.14313629702002917/autoang
        b1 = 0.574206155349291647*autoang**2
        b2 = 0.261455788959452817*autoang
        f0 = -3.13922813873961903/autoev
        f2 = 4.71763556424035180/autoev
        q1 = -3.62530532486565704/autoev
        q2 = -23.1485100146555922/autoev
        q3 = -20.1148021494870548/autoev

c original 
c      rc = 6.78040d0   /autoang
c      h  = 1.41584d0   /autoang
c      e1 = 2.65219d2   /autoev
c      e2 = 7.67238d-3  /autoev
c      r1 = 1.0688d0    /autoang
c      r2 = 6.45771d0   /autoang
c      alp1 = 2.09043d0 *autoang
c      alp2 = 1.03062d0 *autoang
c      delta = 0.10432d0/autoev
c      a  = 5.47491d-2
c      r3 = 2.74668d0    /autoang
c      r4 = -6.92647d0   /autoang
c      b1 = 1.92116d0    *autoang**2
c      b2 = 0.42589d0    *autoang
c      f0 = -2.81468d0   /autoev
c      f2 = 5.57686d0    /autoev
c      q1 = -6.24655d0   /autoev
c      q2 = -21.52522d0   /autoev
c      q3 = -15.30493d0   /autoev

      v=0.d0
      do 2 i=1,natom
      v2=0.d0
      rhobar=0.d0
      do 3 j=1,natom
      if (i.eq.j) go to 3
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)

c 2-body part
      psiarg = (rr-rc)/h
      if (psiarg.ge.0.d0) then
        psicut = 0.d0
      else
        psicut = psiarg**4/(1.d0+psiarg**4)
      endif
      em1 = dexp(-2.d0*alp1*(rr-r1))-2.d0*dexp(-alp1*(rr-r1))
      em2 = dexp(-2.d0*alp2*(rr-r2))-2.d0*dexp(-alp2*(rr-r2))
      v2 = v2 + (e1*em1+e2*em2+delta)*psicut

c rhobar
      rho = (a*dexp(-b1*(rr-r3)**2)+dexp(-b2*(rr-r4)))*psicut
      rhobar = rhobar + rho  ! rho is 2x because we are looping j>i

    3 continue
c embedding term
      xksum = q1*(rhobar-1.d0)**3
     &      + q2*(rhobar-1.d0)**4
     &      + q3*(rhobar-1.d0)**5

      ef = f0 + 0.5d0*f2*(rhobar-1)**2 + xksum
      v = v + 0.5d0*v2 + ef
    2 continue

      return
      end
