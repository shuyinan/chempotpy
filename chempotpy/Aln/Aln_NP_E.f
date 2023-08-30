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
c   Common name:                   NP-E
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c
c   Notes:  Many-body aluminum potential energy function.  The functional
c           form is from Ref. 1.  The parameters were re-optimized in Ref. 2
c           against a data set of energies for aluminum clusters and
c           nanoparticles and bulk data.  Reference 3 provides futher
c           background but is not a required reference for this potential.
c
c  References: (1) Y. Mishin, D. Farkas, M. J. Mehl, and D. A. Papaconstantopoulos
c              Mat. Res. Soc. Symp. Proc. 583, 535, (1999).
c              (2) A. W. Jasper, N. E. Schultz, and D. G. Truhlar, "Analytic 
c              Potential Energy Functions for Simulating Aluminum 
c              Nanoparticles," in preparation.  
c              (3) A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c              and D. G. Truhlar, "Analytic Potential Energy Functions
c              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
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

        rc= 6.21230561520981794/autoang
        h = 1.23419121131286391/autoang
        e1 = 323.306599351550119/autoev
        e2 = 0.620692841774021006E-02/autoev
        r1 = 0.849576438380276855/autoang
        r2 = 4.63688776461126384/autoang
        alp1 = 2.99843114651358267*autoang
        alp2 = 1.55910921885380915*autoang
        delta = 0.794681217776774529E-01/autoev
        a = 0.637166953301545880E-01
        r3 = 2.49022329534347442/autoang
        r4 = -6.69553673168315022/autoang
        b1 = 1.42436502858883451*autoang**2
        b2 = 0.457015505777396858*autoang
        f0 = -3.50601555923931452/autoev
        f2 = 7.33764124826237207/autoev
        q1 = -5.98181779608959907/autoev
        q2 = -19.0233990007299987/autoev
        q3 = -12.7919662283307112/autoev

c original parameters
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
