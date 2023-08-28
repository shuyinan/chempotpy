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
C   Functional form:               Extended Rydberg plus Extended Bond Angle
C   Common name:                   ER2+EBA
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

        xk1= 2.24723986321445990
        xk2= 0.151001465559355147
        xk3= 0.705227161700048844
        deb= 0.766194430874450494/autoev
        reb= 5.08099658036150448
        a1b= 0.668637029799706850
        a2b= 0.813190034196385025
        a3b= 0.357889594528578381
        au2b= 0.592427943331704898
        bu2b= 6.5d0/autoang
        au23= 1.19960918417195894

      do 3 i=1,natom
      r(i,i)=0.d0
      fu2(i,i)=1.d0
      do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      r(i,j)=dsqrt(dx*dx+dy*dy+dz*dz)
      r(j,i)=r(i,j)
      fu2(j,i)=1.d0
      fu2(i,j)=1.d0
    3 continue

c do BA term, needs xk1,xk2,xk3
      do 10 i=1,natom                          ! alpha loop
          sil=0.d0
          do 11 l=1,natom                      ! gamma loop
          if(l.eq.i) go to 11                  ! gamma != alpha
          ril=r(i,l)
          if(ril.ge.bu2) go to 11              ! cutoff for Ralpha,gamma
          ac=au23*(1.d0-1.d0/(1.d0-ril/bu2))   ! cutoff for Ralpha,gamma
          term=dexp(-xk2*ril**xk3+ac)          ! cutoff times Ralpha,gamma part of chialpha,beta
          sil=sil+term                         ! sum
   11     continue

        do 15 j=i+1,natom                      ! beta loop
          slj=0.d0
          do 12 l=1,natom                      ! delta loop
          if(l.eq.j) go to 12                  ! delta != beta
          rlj=r(l,j)
          if(rlj.ge.bu2) go to 12              ! cutoff for Rbeta,delta
          ac=au23*(1d0-1d0/(1d0-rlj/bu2))      ! cutoff for Rbeta,delta
          term=dexp(-xk2*rlj**xk3+ac)          ! cutoff times Rbeta,delta part of chialpha,beta
          slj=slj+term                         ! sum
   12     continue

          rij=r(i,j)
          sij=0.d0
          if(rij.ge.bu2) go to 13              ! cutoff for Ralpha,beta
          ac=au23*(1d0-1d0/(1d0-rij/bu2))      ! cutoff for Ralpha,beta
          sij=dexp(-xk2*rij**xk3+ac)           ! cutoff times overcounted Ralpla,beta term in chi
   13     continue

c         s=chi, need to subtract off overcounted Ralpha,beta term (sij)
c         and mutiply the two sums to get the double sum
          s=xk1*(sil-sij)*(slj-sij)

          e2=dexp(-2d0*s)                     ! equivalent to 1-dtanh(s)
          fu2(i,j)=2d0*e2/(1d0+e2)            ! equivalent to 1-dtanh(s)
          fu2(i,j)=1.d0-fu2(i,j)
          fu2(j,i)=fu2(i,j)
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
	
