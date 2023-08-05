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
C   Functional form:               Extended Rydberg plus Embedded Atom
C   Common name:                   ER2+MisFMP
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
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)
      dimension x(maxatom),y(maxatom),z(maxatom)
c
c  calculate the potential energy
c

        de= 1.71013678553981441/autoev
        re= 5.08182706399609163
        a1= 1.24074255007327805
        a2= 0.551880801172447422
        a3= 0.129970688812896917
        au2= 0.143243771372740580
        bu2= 6.50000000000000000/autoang

        a= 0.264359550561797718
        r3= 3.23400097703957012/autoang
        r4= -9.98827552515876960/autoang
        b1= 2.00000000000000000*autoang**2
        b2= 1.48949682462139710*autoang
        f0= -2.08891060087933589/autoev
        f2= 7.20566682950659487/autoev
        q1= -4.45041524181729375/autoev
        q2= -18.7274059599413789/autoev
        q3= -15.3981436248168055/autoev
        au23= 2.62012701514411361

      v=0.d0
      ef = 0.d0
      ef2 = 0.d0
      do 2 i=1,natom
      v2=0.d0
      ef2 = 0.d0
      rhobar=0.d0
      do 3 j=1,natom
      if (i.eq.j) go to 3
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      if (rr.ge.bu2) go to 3

c 2-body part
c ERCO2-6.5
      rh = rr-re
      poly=1+a1*rh+a2*rh**2+a3*rh**3
      ee=dexp(-a1*rh)
      co=0.d0
      co3=0.d0
      if (rr.lt.bu2) co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
      if (rr.lt.bu2) co3=dexp(au23*(1.d0-1.d0/(1.d0-rr/bu2)))
      v2=-de*ee*poly*co

c rhobar
      rho = (a*dexp(-b1*(rr-r3)**2)+dexp(-b2*(rr-r4)))*co3
      rhobar = rhobar + rho
c subtract two-body part
      xksum2 = q1*(rho-1.d0)**3
     &       + q2*(rho-1.d0)**4
     &       + q3*(rho-1.d0)**5
      ef2 = f0 + 0.5d0*f2*(rho-1.d0)**2 + xksum2

    3 continue
c embedding term
      xksum = q1*(rhobar-1.d0)**3
     &      + q2*(rhobar-1.d0)**4
     &      + q3*(rhobar-1.d0)**5

      ef = f0 + 0.5d0*f2*(rhobar-1.d0)**2 + xksum

      v = v + 0.5d0*v2 + ef - ef2
    2 continue

      return
      end
