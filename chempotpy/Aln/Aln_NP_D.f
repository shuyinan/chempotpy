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
c   Common name:                   NP-D
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
c  References: (1) H. Gollisch, Surface Science 166 (1986), 87.
c              (2) A. W. Jasper, N. E. Schultz, and D. G. Truhlar, "Analytic 
c              Potential Energy Functions for Simulating Aluminum 
c              Nanoparticles," in preparation.  
c              (3) A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c              and D. G. Truhlar, "Analytic Potential Energy Functions
c              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c
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

c  calculate and tabulate the distances between all atoms

        aa= 772.011042684632230/autoev
        bb= 17.6185013563156758/autoev
        a= 3.01309706914042241*autoang
        b= 1.65073816975429399*autoang
        en= 0.588540914437484020

      v=0.d0
      do i=1,natom
      u=0.d0
      ud=0.d0    
      do j=1,natom
      if(i.ne.j) then
       dx=x(i)-x(j)
       dy=y(i)-y(j)
       dz=z(i)-z(j)
       rr=dsqrt(dx*dx+dy*dy+dz*dz)
       u=u+dexp(-rr*a)
       ud=ud+dexp(-rr*b)
      endif
      enddo
      v=v+aa*u-bb*ud**en
      enddo
      return
      end
