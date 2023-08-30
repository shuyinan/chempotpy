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
      double precision :: dvdx(natoms), dvdy(natoms), dvdz(natoms)
      double precision :: v
      integer :: istate, iatom, i

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
c      call pot(cx, cy, cz, v, dvdx, dvdy, dvdz, natoms, 10000)

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
c   Common name:                   NP-B
c   Number of derivatives:         1
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-1
c
c   Notes:  Many-body aluminum potential energy function.  The functional 
c           form is from Ref. 1.  The parameters were re-optimized in Ref. 2
c           against a data set of energies for aluminum clusters and 
c           nanoparticles and bulk data.  Reference 3 provides futher
c           background but is not a required reference for this potential.
c           This coding contains analytic first derivatives.
c           Note that eq 37 of Reference 3 has two typos related to this 
c           potential. A minus sign is missing in front of phi-sub-0, and 
c           the sign in front of the variable "d" should be "+".
c
c  References: (1) J. Mei and J. W. Davenport, Phys. Rev. B, 46, 21, (1992).
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
c  dvdx,dvdy,dvdz   --- one-dimensional arrays of the gradients with respect
c                       to the Cartesian coordinates of the system (output)
c  natom   --- number of atoms in the system (input)
c  maxatom --- maximum number of atoms (input)

      implicit double precision(a-h,o-z)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension cc(6),s(3)
      dimension d_v2(3,maxatom),d_rho(3,maxatom),  ! Derivs
     & dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)

c Original Literature Parameters (Ref. 1)
c        ec    =   3.39d0/autoev
c        phi0  =   0.1318d0/autoev
c        r0    =   2.8638d0/autoang
c        alpha =   4.6d0
c        beta  =   7.10d0
c        gamma =   7.34759d0
c        delta =   7.35d0
c        cc(1) =   0.64085d0
c        cc(2) =  -6.83764d0
c        cc(3) =  26.75616d0
c        cc(4) = -47.16495d0
c        cc(5) =  36.18925d0
c        cc(6) =  -8.60834d0
c        s(1)  =  12.d0
c        s(2)  =   6.d0
c        s(3)  =  24.d0
c        rn    =   1.75d0*r0
c        rc    =   1.95d0*r0

c Cluster-parameterized Parameters (Ref. 3)
c        ec    =   3.648900830d0/autoev
c        phi0  =   0.395759648d0/autoev
c        r0    =   2.624670249d0/autoang
c        alpha =   5.194382022d0
c        beta  =   4.748803127d0
c        gamma =   5.754323400d0
c        delta =   7.663800684d0
c        cc(1) =   0.010571568d0
c        cc(2) =  -7.999267220d0
c        cc(3) =  31.689301416d0
c        cc(4) = -43.634587200d0
c        cc(5) =  28.121641426d0
c        cc(6) =  -6.751099169d0
c        s(1)  =   8.657547631d0
c        s(2)  =   4.169614069d0
c        s(3)  =  27.910600879d0
c        rn    =   1.75d0*r0
c        rc    =   1.95d0*r0

c Nanoparticle-parameterized Parameters (Ref. 2)
        ec    =   2.833661628d0/autoev
        phi0  =   0.209474568d0/autoev
        r0    =   2.759835989d0/autoang
        alpha =   4.953631991d0
        beta  =   5.202672172d0
        gamma =   5.824302949d0
        delta =   8.968682037d0
        cc(1) =   0.433294196d0
        cc(2) =  -7.305279256d0
        cc(3) =  29.818956621d0
        cc(4) = -54.437991632d0
        cc(5) =  48.412067298d0
        cc(6) = -15.525225110d0
        s(1)  =   6.927645227d0
        s(2)  =   3.861172975d0
        s(3)  =  15.498062621d0
        rn    =   1.75d0*r0
        rc    =   1.95d0*r0

c Initialize
      v = 0.d0
      do i=1,natom
      dvdx(i) = 0.d0
      dvdy(i) = 0.d0
      dvdz(i) = 0.d0
      enddo

c Double loop (main loop)
      do 2 i=1,natom
      rho = 0.d0
      v2 = 0.d0
      do k=1,3
      do j=1,natom
      d_rho(k,j) = 0.d0
      d_v2(k,j) = 0.d0
      enddo
      enddo
      do 3 j=1,natom
      if (i.eq.j) go to 3
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)

c Cutoff function
      if (rr.le.rn) then
             q=1.d0
             dqdr=0.d0  ! Derivs
      elseif(rr.ge.rc) then
             go to 3
             q=0.d0
             dqdr=0.d0  ! Derivs
      else
             xx=(rr-rn)/(rc-rn)
             q=((1.d0-xx)**3)*(1.d0+3.d0*xx+6.d0*xx**2)
             dxxdr=1.d0/(rc-rn)  ! Derivs
             dqdr=((1.d0-xx)**3)*(3.d0+12.d0*xx)*dxxdr   ! Derivs
     &   -(3.d0*(1.d0-xx)**2)*(1.d0+3.d0*xx+6.d0*xx**2)*dxxdr
      endif

c Pair potential phi
      dr = rr/r0-1.d0
      phi = -0.5d0*phi0*(1.d0+delta*dr)*dexp(-gamma*dr)
      v2 = v2 + phi*q
      drdr = 1.d0/r0  ! Derivs
      dphidr=-0.5d0*phi0*(1.d0+delta*dr)*(-gamma*drdr)*dexp(-gamma*dr)   ! Derivs
     &       -0.5d0*phi0*delta*drdr*dexp(-gamma*dr)
      dphiqdr = dqdr*phi+dphidr*q  ! Derivs

c Embedding potential rho
      drhodr = 0.d0  ! Derivs
      do l=0,5
      ef = (cc(l+1)/12.d0)*(r0/rr)**l
      rho = rho + ef*q
      if (l.ne.0) drhodr = drhodr  ! Derivs
     &  - (q/rr)*dble(l)*(cc(l+1)/12.d0)*(r0/rr)**l+ef*dqdr
      if (l.eq.0) drhodr = drhodr + ef*dqdr  ! Derivs
      enddo

c v2 (pairwise repulsion) derivs
      d_v2(1,i) = d_v2(1,i) + dphiqdr*dx/rr
      d_v2(1,j) = d_v2(1,j) - dphiqdr*dx/rr
      d_v2(2,i) = d_v2(2,i) + dphiqdr*dy/rr
      d_v2(2,j) = d_v2(2,j) - dphiqdr*dy/rr
      d_v2(3,i) = d_v2(3,i) + dphiqdr*dz/rr
      d_v2(3,j) = d_v2(3,j) - dphiqdr*dz/rr
c rho (density) derivs
      d_rho(1,i) = d_rho(1,i) + drhodr*dx/rr
      d_rho(1,j) = d_rho(1,j) - drhodr*dx/rr
      d_rho(2,i) = d_rho(2,i) + drhodr*dy/rr
      d_rho(2,j) = d_rho(2,j) - drhodr*dy/rr
      d_rho(3,i) = d_rho(3,i) + drhodr*dz/rr
      d_rho(3,j) = d_rho(3,j) - drhodr*dz/rr
    3 continue

      sterm=0.d0
      dstermdrho=0.d0  ! Derivs
      if (rho.gt.0.d0) then
      do m=1,3
      xm=dsqrt(dble(m))
      sterm=sterm+s(m)*dexp(-gamma*(xm-1.d0))
     c                *(1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho))
     c                *(rho)**(xm*gamma/beta)
      dstermdrho = dstermdrho + s(m)*dexp(-gamma*(xm-1.d0))*(      ! Derivs
     c  (1.d0+delta*(xm-1.d0)-xm*(delta/beta)*dlog(rho))
     c                    *(xm*gamma/beta)*(rho)**(xm*gamma/beta-1.d0)
     c  +(-xm*(delta/beta)*(1.d0/rho))*(rho)**(xm*gamma/beta) )
      enddo
      endif

      ab=alpha/beta
      bigef = 0.d0
      dbigefdrho = 0.d0
      if (rho.gt.0.d0) then
      bigef = -ec*(1.d0-ab*dlog(rho))*rho**ab + 0.5d0*phi0*sterm
      dbigefdrho = -ec*(1.d0-ab*dlog(rho))*(ab)*rho**(ab-1.d0)    ! Derivs
     &             -ec*(-ab*(1.d0/rho))*rho**ab
     &             +0.5d0*phi0*dstermdrho
      endif

c Combine derivs
      do j=1,natom
      dvdx(j) = dvdx(j) + d_v2(1,j) + dbigefdrho*d_rho(1,j)
      dvdy(j) = dvdy(j) + d_v2(2,j) + dbigefdrho*d_rho(2,j)
      dvdz(j) = dvdz(j) + d_v2(3,j) + dbigefdrho*d_rho(3,j)
      enddo

      v = v + v2 + bigef
    2 continue

      return
      end
