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
C   Common name:                   ER2+MeiD
C   Number of derivatives:         0
C   Number of bodies:              variable
C   Number of electronic surfaces: 1
C   Interface:                     HO-MM-0
C
c   Notes:  Many-body aluminum potential energy function.  The parameters
c           for this PEF have been parameterized against a data set of
c           aluminum cluster energies and bulk data.
c           This PEF has a cutoff at 6.5 angstroms.
c           Note that eq 37 of J. Phys. Chem. B 108, 8896 (2004) has two
c           typos related to this potential. A minus sign is missing in 
c           front of phi-sub-0, and the sign in front of the variable "d"
c           should be "+".
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
      dimension cc(6),s(3)


c accurate two-body
      de= 1.71013678553981441/autoev
      re= 5.08182706399609163
      a1= 1.24074255007327805
      a2= 0.551880801172447422
      a3= 0.129970688812896917
      au2= 0.143243771372740580
      bu2= 6.50000000000000000/autoang
c     potential, includes cutoff
      vER2=0.d0
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
      vER2=vER2-de*ee*poly*co
      endif
    1 continue

c  calculate the potential energy
        ec= 2.00214948705422557/autoev
        phi0= 0.110210063507572065/autoev
        r0= 2.15915974596971205/autoang
        alpha= 7.06277479237909134
        beta= 6.58251099169516340
        gamma= 6.24914509037616028
        delta= 6.39066927210552027
        cc(1)= 0.407034684904738653
        cc(2)= -6.44787493893502628
        cc(3)= 36.1221299462628238
        cc(4)= -47.8456277479237926
        cc(5)= 25.8260869565217384
        cc(6)= -5.49975574010747437
        s(1)= 12.6893014167073765
        s(2)= 6.47234978016609652
        s(3)= 26.9702002931118692
      rn=1.75d0*r0
      rc=1.95d0*r0

      v = 0.d0
      v_pair = 0.d0
      do 2 i=1,natom
      rho = 0.d0
      do 3 j=1,natom
      if (i.eq.j) go to 3
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)

c cutoff function
      if (rr.le.rn) then
             q=1.d0
      elseif(rr.ge.rc) then
             q=0.d0
      else
             xx=(rr-rn)/(rc-rn)
             q=((1.d0-xx)**3)*(1.d0+3.d0*xx+6.d0*xx**2)
      endif

c embedding potential rho
      rho_pair = 0.d0
      do l=0,5
      ef = (cc(l+1)/12.d0)*(r0/rr)**l
      rho = rho + ef*q
      rho_pair = rho_pair + ef*q
      enddo

      sterm=0.d0
      if (rho_pair.gt.0.d0) then
      do m=1,3
      xm=dsqrt(dble(m))
      sterm=sterm+s(m)*dexp(-gamma*(xm-1.d0))
     c         *(1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho_pair))
     c         *(rho_pair)**(xm*gamma/beta)
      enddo
      endif

      ab=alpha/beta
      bigef = 0.d0
      if (rho_pair.gt.0.d0) then
      bigef = -ec*(1.d0-ab*dlog(rho_pair))*rho_pair**ab 
     &      + 0.5d0*phi0*sterm
      endif

      v_pair = v_pair + bigef

    3 continue

      sterm=0.d0
      if (rho.gt.0.d0) then
      do m=1,3
      xm=dsqrt(dble(m))
      sterm=sterm+s(m)*dexp(-gamma*(xm-1.d0))
     c                *(1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho))
     c                *(rho)**(xm*gamma/beta)
      enddo
      endif

      ab=alpha/beta
      bigef = 0.d0
      if (rho.gt.0.d0) then
      bigef = -ec*(1.d0-ab*dlog(rho))*rho**ab + 0.5d0*phi0*sterm
      endif
      v = v + bigef
    2 continue

      v = v + vER2 - v_pair

      return
      end
