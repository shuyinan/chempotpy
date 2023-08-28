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
c   Functional form:               Pseudopotential theory
c   Common name:                   rePetW
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
c              D.G. Pettifor and M.A. Ward, Solid State Commun.
c              49, 291 (1984).
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
c
c
c  calculate the potential energy
c

c original parameters
      a1=7.964     
      a2=1.275
      a3=0.030
      alpha1=-1.38544236
      alpha2=2.61380509
      alpha3=1.35402642     
      ek1=0.28926444
      ek2=1.19414294
      ek3=1.77638033
      ekap1=1.47042756
      ekap2=1.29427294
      ekap3=0.51733832

c new parameters
      a1 = 67.1568148510014709
      a2 = 0.164142647777234973
      a3 = 57.3815339521250607
      alpha1 = 47.5280898876404478
      alpha2 = 96.8734733756717077
      alpha3 = 75.0854909623839717
      ek1 = 1.32467024914509035
      ek2 = 1.23111871030776743
      ek3 = 0.274157303370786531
      ekap1 = 3.50244259892525633
      ekap2 = 0.513678553981436248
      ekap3 = 3.57093307278944794

      v=0d0
          do 2 i=1,natom-1
          do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      v1=(18d0/rr)*(a1*dcos(ek1*rr+alpha1)*dexp(-ekap1*rr)
     *+a2*dcos(ek2*rr+alpha2)*dexp(-ekap2*rr)+a3*dcos(ek3*rr+alpha3)
     **dexp(-ekap3*rr))
      v=v+v1
    3 continue
    2 continue
      return
      end
