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
c   Functional form:               Embedded atom/Effective medium theory
c   Common name:                   Jac
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c
c   Notes:  Many-body aluminum potential energy function.  For a recent
c           discussion of this potential, see:
c           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c           and D. G. Truhlar, "Analytic Potentials Energy Functions
c           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c
c  References: K. W. Jacobsen, Comments Cond. Mat. Phys. 14,  
c		  129 (1988). 
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

c  calculate the potential energy

      eta1=0.3
      en=0.007
      s=3.0
      eta=2.0

      en     =  0.04724  *autoang**3
      s      =  1.5875   /autoang
      eta    =  3.7795   *autoang
      e0     = -3.28     /autoev
      e2     =  1.12     /autoev
      e3     = -0.35     /autoev
      alpha  =  189.7   /(autoev*autoang**3)
      eta1   =  0.5669   *autoang

      pi=dacos(-1.d0)
      beta=(2.d0**(-0.5))*((16.d0*pi/3.d0)**(1.d0/3.d0))
c      beta=1.81
      eta2=(eta+eta1)/beta
      ep=eta/(eta+eta1)
      bs=beta*s
      fra=1.d0/12.d0

      v=0d0
           do 2 i=1,natom
      u=0d0
          do 3 j=1,natom
          if(j.eq.i) go to 3          
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      u=u+dexp(-eta2*(rr-bs))   ! EAS(1) w/o exponent
    3 continue
      ud=0d0
          do 4 m=1,natom
          if(m.eq.i) go to 4
      dx=x(i)-x(m)
      dy=y(i)-y(m)
      dz=z(i)-z(m)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      ud=ud+dexp(-eta*(rr/beta-s))  ! EAS(2)
    4 continue
      uno=(fra*u)**ep   ! EAS(1)
      udno=fra*ud       ! EAS(2)
      vv=e0+e2*(uno-1d0)**2+e3*(uno-1d0)**3
     *+alpha*en*(uno-udno)
    2 v=v+vv
      return
      end
