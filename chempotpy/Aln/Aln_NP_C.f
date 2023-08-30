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
c   Common name:                   NP-C
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
c  References: (1) F. H. Streitz, and J. W. Mintmire, Phys. Rev. B, 50, 11996 (1994).
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
      dimension x(maxatom),y(maxatom),z(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)

      rstar = 2.85655354324007371 / autoang
      alpha = 2.13479173244937614 * autoang
      beta  = 1.16888101913568532 * autoang
      aa    = 0.716877851705074343 / autoev
      bb    = 0.667088456913949440E-01 / autoev
      cc    = 0.333776366433971738 / autoev

      v=0.d0

c pairwise term
      phi=0.d0
      ff=0.d0
      do 5 i=1,natom

      rho=0.d0
      do 6 j=1,natom
      if (j.ne.i) then
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      rrel=rr-rstar
      rho=rho+dexp(-beta*rrel)
      endif
      if(i.gt.j) then
      p1=2.d0*bb*dexp(-0.5d0*beta*rrel)
      p2=-cc*(1.d0+alpha*rrel)*dexp(-alpha*rrel)
      phi=phi+p1+p2
      endif
    6 continue

      ff=ff-aa*dsqrt(rho)
    5 continue

      v=phi+ff

      return
      end

