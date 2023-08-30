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
c   Common name:                   MeiD
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c
c   Notes:  Many-body aluminum potential energy function.
c
c  References: J. Mei and J. W. Davenport, Phys. Rev. B, 46, 21, (1992).
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
      dimension cc(6),s(3)

c  calculate the potential energy

      ec = 3.39d0/autoev
      phi0 = 0.1318d0/autoev
      r0 = 2.8638d0/autoang
      alpha = 4.6d0
      beta  = 7.10d0
      gamma = 7.34759d0
      delta = 7.35d0
      cc(1) = 0.64085d0
      cc(2) = -6.83764d0
      cc(3) = 26.75616d0
      cc(4) = -47.16495d0
      cc(5) = 36.18925d0
      cc(6) = -8.60834d0
      s(1) = 12.d0
      s(2) = 6.d0
      s(3) = 24.d0
      rn=1.75d0*r0
      rc=1.95d0*r0

      v = 0.d0
      do 2 i=1,natom
      rho = 0.d0
      v2 = 0.d0
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

c pair potential phi
      dr = rr/r0-1.d0
      phi = -0.5d0*phi0*(1.d0+delta*dr)*dexp(-gamma*dr)
      v2 = v2 + phi*q

c embedding potential rho
      do l=0,5
      ef = (cc(l+1)/12.d0)*(r0/rr)**l
      rho = rho + ef*q
      enddo
    3 continue

      sterm=0.d0
      if (rho.ne.0.d0) then
      do m=1,3
      xm=dsqrt(dble(m))
      sterm=sterm+s(m)*dexp(-gamma*(xm-1.d0))
     c                *(1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho))
     c                *(rho)**(xm*gamma/beta)
      enddo
      endif

      ab=alpha/beta
      bigef = 0.d0
      if (rho.ne.0.d0) then
      bigef = -ec*(1.d0-ab*dlog(rho))*rho**ab + 0.5d0*phi0*sterm
      endif
      v = v + v2 + bigef
    2 continue

      return
      end
