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
C   Functional form:               Extended Rydberg plus Extended Coordination Number
C   Common name:                   ER2+ECN
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
      dimension x(maxatom),y(maxatom),z(maxatom),gcorr(maxatom)
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

        deb= 1.26365412799218380/autoev
        reb= 4.51734245236932086
        a1b= 0.766243282852955510
        a2b= 0.416023448949682462
        a3b= 0.260429897410845157
        au2b= 1.08568637029799708
        bu2b= 6.5d0/autoang
        g1= 1.80898876404494380
        g2= 8.09721543722520742/autoang
        du2= 2.21812408402540306
        gu2= 0.853639472398632138
        gzero= 4.62188568637029817


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

c CALCULATE CN FUNCTION
c     first calculate g_alpha
      g1e=dexp(g1)
      do 20 i=1,natom                         ! alpha loop
       gcorr(i)=0d0
       do 21 j=1,natom                        ! beta loop
        if(i.eq.j) go to 21                   ! beta != alpha
        rij=r(i,j)
        if (rij.ge.g2) go to 21               ! cutoff term counts nearby atoms (g_alphabeta)
        term=dexp(-g1/(1d0-(rij/g2)))         ! cutoff term counts nearby atoms (g_alphabeta)
        gcorr(i)=gcorr(i)+term                ! sum over beta
  21   continue
       gcorr(i)=gcorr(i)*g1e                  ! g_alpha
  20  continue

c     then calculate CN term
      do 30 i=1,natom                         ! alpha loop
      do 31 j=1,natom                         ! beta loop
       if(i.eq.j) go to 31                    ! beat != alpha
       rij=r(i,j)

       term=0.d0
       if (rij.ge.g2) go to 32                ! g_alphabeta
       term=dexp(g1*(1d0-1d0/(1d0-(rij/g2)))) ! g_alphabeta
  32   continue

       cz1=(gcorr(i)-term)/gzero              ! ratio of coordination for alpha atom
       abscz1=dabs(cz1)
       if(abscz1.le.10d-12)cz1=abscz1         ! eliminates negative values due to roundoff?
       es1=1d0/(1d0+cz1**gu2)                 ! compute alpha part of f for u2 term

       cz2=(gcorr(j)-term)/gzero              ! ratio of coordination for beta atom
       abscz2=dabs(cz2)
       if(abscz2.le.10d-12)cz2=abscz2         ! eliminates negative values due to roundoff?
       ec1=1d0/(1d0+cz2**gu2)                 ! compute beta part of f for u2 term

c       fu2(i,j)=1d0+du2*(es1*ec1-1d0)         ! put alpha and beta parts together:  1 + f for regular formalism
       fu2(i,j)=du2*(es1*ec1-1d0)         ! put alpha and beta parts together:  f for extendended formalism
c       print *,es1,ec1,fu2(i,j)
   31 continue
   30 continue


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
     &   -deb*eeb*polyb*cob*fu2(i,j)
    1 continue

      return

      end
	
