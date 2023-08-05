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
c
             subroutine pot(x,y,z,v,natom,maxatom)

c   System:                        Aln
c   Functional form:               Pseudopotential theory
c   Common name:                   PetW
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c   Notes:  Many-body aluminum potential energy function.  For a recent
c           discussion of this potential, see:
c           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c           and D. G. Truhlar, "Analytic Potentials Energy Functions
c           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c
c  References: D.G. Pettifor and M.A. Ward, Solid State Commun. 
c		  49, 291 (1984). 
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
      parameter(autoang=0.529177249d0)
      parameter(autoev=27.2113961d0)

      v=0.d0
      zz = 3.d0
      a1=7.964d0
      a2=1.275d0
      a3=0.030d0
      pi=dacos(-1.d0)
      alpha1=-0.441*pi
      alpha2= 0.832*pi
      alpha3= 0.431*pi

c      rs=2.07d0
c      xlam=(rs/pi)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
c      xlam2kf=xlam*pi**2
c      twokf=2.d0*pi/xlam2kf
      twokf=3.504d0*autoang
      xk1=  0.156d0*twokf
      xk2=  0.644d0*twokf
      xk3=  0.958d0*twokf
      xkap1=0.793d0*twokf
      xkap2=0.698d0*twokf
      xkap3=0.279d0*twokf

      do 2 i=1,natom-1
      do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      v1=((2.d0*zz**2)/rr)*(a1*dcos(xk1*rr+alpha1)*dexp(-xkap1*rr)
     *                     +a2*dcos(xk2*rr+alpha2)*dexp(-xkap2*rr)
     *                     +a3*dcos(xk3*rr+alpha3)*dexp(-xkap3*rr))
      v=v+v1
    3 continue
    2 continue
      return
      end
