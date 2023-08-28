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
C   Functional form:               Extended Rydberg plus Extended Screening
C   Common name:                   ER2+ES
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
      dimension x(maxatom),y(maxatom),z(maxatom),c(10)
      dimension fu2(maxatom,maxatom),r(maxatom,maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)

        de= 1.71013678553981441/autoev
        re= 5.08182706399609163
        a1= 1.24074255007327805
        a2= 0.551880801172447422
        a3= 0.129970688812896917
        au2= 0.143243771372740580
        bu2= 6.50000000000000000/autoang

        xk1= 2.02364435759648265
        xk2= 0.444699560332193433
        xk3= 2.37791890571568132
        deb= 1.77156814851001476/autoev
        reb= 4.71494870542256983
        a1b= 1.12725940400586233
        a2b= 0.496384953590620404
        a3b= 0.501460674157303332
        au2b= 0.467684416218856813
        bu2b= 6.5d0/autoang
        au23= 0.665852467024914546E-01

      do 3 i=1,natom
      r(i,i)=0.d0
      fu2(i,i)=0.d0
      do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r(i,j)=dsqrt(dx*dx+dy*dy+dz*dz)
      r(j,i)=r(i,j)
      fu2(j,i)=0.d0
      fu2(i,j)=0.d0
    3 continue

c do BA term, needs xk1,xk2,xk3
      do 10 i=1,natom                          ! alpha loop
        do 15 j=i+1,natom                      ! beta loop

          s=0.d0
          do 12 l=1,natom                      ! gamma loop
          if(l.eq.j.or.l.eq.i) go to 12        ! gamma != alpha or beta
          rjl=r(l,j)
          rij=r(i,j)
          ril=r(l,i)
          if(rjl.ge.bu2) go to 12              ! cutoff for Rbeta,delta
          if(rij.ge.bu2) go to 12              ! cutoff for Rbeta,delta
          if(ril.ge.bu2) go to 12              ! cutoff for Rbeta,delta
          ac=au23*(1d0-1d0/(1d0-rjl/bu2))       ! cutoff for Rbeta,delta
          ac=ac+au23*(1d0-1d0/(1d0-rij/bu2))    ! cutoff for Rbeta,delta
          ac=ac+au23*(1d0-1d0/(1d0-ril/bu2))    ! cutoff for Rbeta,delta
          term=xk1*dexp(-xk2*((ril+rjl)/rij)**xk3+ac)        ! cutoff times Rbeta,delta part of chialpha,beta
c          if (natom.eq.3) print *,i,j,l,rjl,rij,ril,term
          s=s+term                         ! sum
   12     continue

c        e2=dexp(-2d0*s)                       ! equivalent to 1-dtanh(s)
c        fu2(i,j)=2d0*e2/(1d0+e2)              ! equivalent to 1-dtanh(s)
        fu2(i,j)=dtanh(s)
        fu2(j,i)=fu2(i,j)

c        if (natom.eq.3) print *,i,j,fu2(i,j)

   15   continue
      fu2(i,i)=1.d0                           ! shouldn't get used anywhere
   10 continue

c potential, includes cutoff
      v=0.d0
      do 1 i=1,natom
      do 1 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      rho = rr-re
      rhob = rr-reb
      poly=1.d0+a1*rho+a2*rho**2+a3*rho**3
      polyb=1.d0+a1b*rhob+a2b*rhob**2+a3b*rhob**3
      ee=dexp(-a1*rho)
      eeb=dexp(-a1b*rhob)
      co=0.d0
      cob=0.d0
      if (rr.le.bu2) co=dexp(au2*(1.d0-1.d0/(1.d0-rr/bu2)))
      if (rr.le.bu2b) cob=dexp(au2b*(1.d0-1.d0/(1.d0-rr/bu2b)))
      v=v-de*ee*poly*co
     &   +deb*eeb*polyb*cob*fu2(i,j)
    1 continue

      return

      end
	
