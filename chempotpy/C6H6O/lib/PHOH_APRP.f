C   System:                     C6H6O
C   Functional form:            Anchor-Points Reactive Potential (APRP) 
C                               (3x3 diabatic fit)
C   Common name:                PhOH(Coupled PESs)
C   Number of derivatives:      1
C   Number of bodies:           13
C   Number of electronic surfaces: 3
C   Interface: Section-2
C
C   Reference: Ke R. Yang, Xuefei Xu, Jingjing Zheng,
C   and Donald G. Truhlar Chem. Sci. (to be published).
C
C   Notes:
C    1. This routine calculates three coupled singlet PESs for
C       nonadiabatic photodissociation of PhOH to produce
C       PhO + H.
C    2. The final version of this potential was added to POTLIB 
C       on Aug. 3, 2014. A preliminary version in the library before
C       that is not correct.
C
C   Numbering of atoms:
C        H11      H12
C         \      /
C          C5---C6
C         /      \
C   H10--C4       C1---O7 
C         \      /      \
C          C3---C2       H13
C         /      \
C        H9       H8
C
C
C   Units:      bohr for length
C               hartree for energy
C
C   Input:  
C   igrad = 0           Energy only calculation
C         = 1           Energy + Analytic gradient
C   repflag             flag to indicate wheter to use diabatic or 
C                       adiabatic representation (used by ANT program)
C   xx(1,i)             X coordinate of atom i
C   xx(2,i)             Y coordinate of atom i
C   xx(3,i)             Z coordinate of atom i
C
C   Output: 
C   uu(j,j)             diabatic potential of diabatic state j
C   uu(j,k)             diabatic coupling between state j and k 
C   guu(1,i,j,j)        X component of the gradients of diabatical 
C                       potential of state j at atom i
C   guu(2,i,j,j)        Y component of the gradients of diabatical 
C                       potential of state j at atom i
C   guu(3,i,j,j)        Z component of the gradients of diabatical
C                       potential of state j at atom i
C   guu(1,i,j,k)        X component of the gradients of diabatical         
C                       coupling between state j and k at atom i
C   guu(2,i,j,k)        Y component of the gradients of diabatical         
C                       coupling between state j and k at atom i
C   guu(3,i,j,k)        Z component of the gradients of diabatical
C                       coupling between state j and k at atom i
C
C   vv(j)               adiabatic potential of adiabatic state j
C   gvv(1,i,j)          X component of the gradients of adiabatic 
C                       potential of state j at atom i
C   gvv(2,i,j)          Y component of the gradients of adiabatic
C                       potential of state j at atom i
C   gvv(3,i,j)          Z component of the gradients of adiabatic
C                       potential of state j at atom i
C   dvec(1,i,j,k)       X component of nonadiabatic coupling between
C                       state j and k at atom i
C   dvec(2,i,j,k)       Y component of nonadiabatic coupling between
C                       state j and k at atom i
C   dvec(3,i,j,k)       Z component of nonadiabatic coupling between
C                       state j and k at atom i
C   cc(3,3)             3*3 orghtonal matrix diagonalizing UU matrix
C                       (CC^T*UU*CC yield diagonal matrix with adiabatic
C                       energies as diagonal elements)              
C
C   Note:       LAPACK library is needed to diagonize diabatic potential
C               matrix (UU) to yield adiabatic potenitials VV(3)
C
!***********************************************************************

      subroutine pes(x,igrad,p,g,d)

      implicit none
      ! number of electronic state
      integer, parameter :: nstates=3
      integer, parameter :: natoms=13
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: xx(3,natoms), uu(nstates,nstates)
      double precision :: guu(3,natoms,nstates,nstates)
      double precision :: vv(nstates), gvv(3,natoms,nstates)
      double precision :: dvec(3,natoms,nstates,nstates)
      double precision :: cc(nstates,nstates)
      integer :: iatom, idir, j, istate, jstate
      !initialize 
      p=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
        do idir=1,3
          xx(idir,iatom)=x(iatom,idir)/0.529177211
        enddo
      enddo

      call pot(igrad,xx,uu,guu,vv,gvv,dvec,cc,0)

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=vv(istate)*27.211386
        enddo
      elseif (igrad==1) then
        do istate=1,nstates
          p(istate)=vv(istate)*27.211386
        enddo
        do istate=1,nstates
        do iatom=1,natoms
        do idir=1,3
          g(istate,iatom,idir)=gvv(idir,iatom,istate)*51.422067
        enddo
        enddo
        enddo
      elseif (igrad==2) then
        do istate=1,nstates
          p(istate)=vv(istate)*27.211386
        enddo
        do istate=1,nstates
        do iatom=1,natoms
        do idir=1,3
          g(istate,iatom,idir)=gvv(idir,iatom,istate)*51.422067
        enddo
        enddo
        enddo
        do istate=1,nstates
        do jstate=1,nstates
        do iatom=1,natoms
        do idir=1,3
          d(istate,jstate,iatom,idir)=dvec(idir,iatom,istate,jstate)/0.529177211
        enddo
        enddo
        enddo
        enddo
      endif

      endsubroutine


!***********************************************************************
      subroutine pot(igrad,xx,uu,guu,vv,gvv,dvec,cc,repflag)
************************************************************************
* 
************************************************************************
      implicit none

      integer i,j,k,l,m
      integer natm,nsurft,igrad,repflag
      parameter (natm=13)
      parameter (nsurft=3)

      double precision uu(3,3),vv(3)
      double precision guu(3,natm,3,3),uup(3,3),vvp(3),vvm(3)
      double precision guup(3,natm,3,3),gvv(3,natm,nsurft)
      double precision guum(3,natm,3,3)
      double precision xx(3,natm),x(natm),y(natm),z(natm)
      double precision dvec(3,natm,nsurft,nsurft)
!     double precision cc(nsurft,nsurft),ccp(nsurft,nsurft)
      double precision cc(3,3),ccp(3,3)
      double precision stp,tmp
      parameter (stp=1.D-4)

C calculate adiabatic and diabatic energies in hartree
!     write(6,*) 'input geometry'
      do i = 1, natm
         x(i) = xx(1,i)
         y(i) = xx(2,i)
         z(i) = xx(3,i)
!        write(6,'(3F10.6)') xx(:,i)*0.5291772192d0
      enddo

!     vv(:) = 0d0
!     uu(:,:) = 0d0
!     guu(:,:,:,:) = 0d0
!     cc(:,:) = 0d0
      call poten(1,x,y,z,uu,guu,vv,gvv,cc)
! Compute the nonadiabatic coupling vectors DVEC.
      call getdvec3(natm,uu,guu,vv,dvec,cc)

!     write(6,*) 'diabatic energy '
!     write(6,'(5F10.6,a)') uu(1,1),uu(2,2),uu(3,3),uu(1,3),uu(2,3)
!     write(6,'(a,5F10.6,a)') 'adiabatic energy ', vv(:)
!     write(6,'(5F10.6,a)') vv(:)
!     write(6,*) 'diabatic gradients'
!     do i = 1,natm
!        write(6,'(5F12.8)') guu(1,i,1,1),guu(1,i,2,2),guu(1,i,3,3),
!    $                       guu(1,i,1,3),guu(1,i,2,3)
!        write(6,'(5F12.8)') guu(2,i,1,1),guu(2,i,2,2),guu(2,i,3,3),
!    $                       guu(2,i,1,3),guu(2,i,2,3)
!        write(6,'(5F12.8)') guu(3,i,1,1),guu(3,i,2,2),guu(3,i,3,3),
!    $                       guu(3,i,1,3),guu(3,i,2,3)
!     enddo
!     write(6,*) 'adiabatic gradients'
!     do i = 1,natm 
!       write(6,'(5F10.6,a)') gvv(1,i,1),gvv(1,i,2),gvv(1,i,3)
!       write(6,'(5F10.6,a)') gvv(2,i,1),gvv(2,i,2),gvv(2,i,3)
!       write(6,'(5F10.6,a)') gvv(3,i,1),gvv(3,i,2),gvv(3,i,3)
!     enddo
!     write(6,*) 'adiabatic coupling dvec'
!     do i = 1, natm
!       write(6,'(5F10.6,a)') dvec(1,i,1,3),dvec(1,i,2,3)
!       write(6,'(5F10.6,a)') dvec(2,i,1,3),dvec(2,i,2,3)
!       write(6,'(5F10.6,a)') dvec(3,i,1,3),dvec(3,i,2,3)
!     enddo
!     stop
      return
      end
      subroutine getdvec3(nclu,pemd,gpemd,pema,dvec,cc)
        implicit none
        integer, parameter :: nsurft=3
        double precision :: pema(nsurft),pemd(nsurft,nsurft)
        double precision :: gpemd(3,nclu,nsurft,nsurft)
        double precision :: dvec(3,nclu,nsurft,nsurft),cc(nsurft,nsurft)
        integer :: i1,i2,i,j,k,l,nclu
 
!     Program to compute nonadiabatic coupling vectors from the adiabatic energies and
!     the gradients of diabatic potential matrix elements for a multisurface system
!     The formula is taken from ANT07 and is described in M.D. Hack et al., J. Phys. Chem. A 103, 6309 (1999)
!
! compute d
       dvec(:,:,:,:) = 0d0
        do i1 = 1,3
         do i2 = 1,nclu
          do i = 2,nsurft
           do j = 1, i-1
!            dvec(i1,i2,i,j) = 0.d0
             do k = 1, nsurft
              do l = 1, nsurft
               dvec(i1,i2,i,j) = dvec(i1,i2,i,j)
     &                         + cc(k,i)*cc(l,j)*gpemd(i1,i2,k,l)
              enddo
             enddo
             if ((pema(j) - pema(i)) .ne. 0.0d0) then
               dvec(i1,i2,i,j) = dvec(i1,i2,i,j)/(pema(j)-pema(i))
             else
               dvec(i1,i2,i,j) = 0.0d0
             endif
!              write(6,*) "dvec ",i1,i2,i,j,dvec(i1,i1,i,j)
             dvec(i1,i2,j,i) = -dvec(i1,i2,i,j)
           enddo
          enddo
         enddo
        enddo
!       stop
        return
        end subroutine getdvec3

      subroutine poten(igrad,X,Y,Z,uu,guu,vv,gvv,cc)
************************************************************************
* Little subroutine to call evuii and evuij
************************************************************************

      implicit double precision (a-h,o-z)

      integer natm
      parameter (natm=13)
      integer nap
      parameter (nap=4)
      double precision rap(nap)
      data rap/1.93206106d0,2.49443845d0,4.27078097d0,9.44863047d0/

      integer, intent(in) :: igrad
      double precision, intent(in) ::  X(natm),Y(natm),Z(natm)
      integer ir
      double precision roh
      double precision t(nap),dtdr(nap)
      double precision u1,u2,u3
! Add by KRY
      double precision ubm
      double precision dubmdX(3*natm)
!
      double precision du1dX(3*natm),du2dX(3*natm),du3dX(3*natm)
      double precision u13,u23
C      double precision du13dX(3*natm),du23dX(3*natm)
      double precision du12dX(3*natm),du13dX(3*natm),du23dX(3*natm)
! add by J. Zheng
      integer lwork,info
      double precision, intent(out) :: uu(3,3),vv(3),
     &                                 guu(3,natm,3,3)
      double precision :: cc(3,3)
      double precision, allocatable :: work(:)
CC add by K. Yang
      double precision agvv(3*natm,3)
      double precision gvv(3, natm, 3)
      double precision fad(3*natm,3,3)

C Determine rOH and determine ir, t, and dtdr
      roh=bndlen(7,13,X,Y,Z)
      call tent(nap,rap,roh,ir,t,dtdr)

      do i=1,3*natm
        du1dX(i)=0.0d0
        du2dX(i)=0.0d0
        du3dX(i)=0.0d0
      enddo

C Determine u13 and u13 and their Cartesian gradients
C      call evuij(nap,natm,igrad,ir,t,dtdr,X,Y,Z,u13,u23,
C     $           du13dX,du23dX)
C Determin u12, u13, and u23 and their Cartesian gradients  !KRY 20140718
      call evuij(nap,natm,igrad,ir,t,dtdr,X,Y,Z,u12,u13,u23,
     $           du12dX,du13dX,du23dX)

C Determine u11, u22, u33, and their Cartesian gradients
      call evuii(igrad,natm,X,Y,Z,u1,du1dX,11)
      call evuii(igrad,natm,X,Y,Z,u2,du2dX,22)
      call evuii(igrad,natm,X,Y,Z,u3,du3dX,33)

C KRY Born-Mayer repulsion
      call evubm(igrad,natm,X,Y,Z,ubm,dubmdX)
C
C add by J. Zheng
      guu(:,:,:,:) = 0d0
      do i = 1, natm
      do j = 1, 3
C        guu(j,i,1,1) = du1dx(3*i-3+j)
C        guu(j,i,2,2) = du2dx(3*i-3+j)
C        guu(j,i,3,3) = du3dx(3*i-3+j)
        guu(j,i,1,1) = du1dx(3*i-3+j)+dubmdX(3*i-3+j)
        guu(j,i,2,2) = du2dx(3*i-3+j)+dubmdX(3*i-3+j)
        guu(j,i,3,3) = du3dx(3*i-3+j)+dubmdX(3*i-3+j)
C add by K. R. Yang 20140718
        guu(j,i,1,2) = du12dx(3*i-3+j)
        guu(j,i,2,1) = du12dx(3*i-3+j)   ! add du12dx
        guu(j,i,1,3) = du13dx(3*i-3+j)
        guu(j,i,3,1) = du13dx(3*i-3+j)
        guu(j,i,2,3) = du23dx(3*i-3+j)
        guu(j,i,3,2) = du23dx(3*i-3+j)
      enddo
      enddo

      uu(:,:) = 0d0
C      uu(1,1) = u1
C      uu(2,2) = u2
C      uu(3,3) = u3
      uu(1,1) = u1+ubm
      uu(2,2) = u2+ubm
      uu(3,3) = u3+ubm
C add by K. R. Yang 20140718
      uu(1,2) = u12
      uu(2,1) = u12       ! add u12
      uu(1,3) = u13
      uu(3,1) = u13
      uu(2,3) = u23
      uu(3,2) = u23
C diagonalization to get adiabatic energies
! Diagonalize the coupling matrix (from LAPACK)
!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  u(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of u                               !
!output: cc(n,n) = orthonormal eigenvectors of u          !
!        vv(n) = eigenvalues of u in ascending order      !
!---------------------------------------------------------!
      cc(:,:) = uu(:,:)
      lwork = -1
      allocate (work (1) )
      call dsyev('V','L',3,cc,3,vv,work,lwork,info)
      lwork=int(work(1))
      deallocate (work)
      allocate( work(lwork) )
      call dsyev('V','L',3,cc,3,vv,work,lwork,info)
      deallocate (work)
      if(info.ne.0) then
       write(6,*) 'Error in diagolization of U matrix'
       stop
      endif
C Test
      if (igrad .eq. 1) then
!     write(*,*) "C matrix:"
!     do i=1,3
!       write(*,'(3F12.8)') (cc(i,j),j=1,3)
!     enddo
     
      gvv(:,:,:) = 0d0
C      do i1=1,natm
C      do i2=1,3
      do i1=1,3
        do i2=1,natm
        do j=1,3
          do k=1,3
            gvv(i1,i2,1)=gvv(i1,i2,1)+cc(k,1)*cc(j,1)*guu(i1,i2,j,k)
            gvv(i1,i2,2)=gvv(i1,i2,2)+cc(k,2)*cc(j,2)*guu(i1,i2,j,k)
            gvv(i1,i2,3)=gvv(i1,i2,3)+cc(k,3)*cc(j,3)*guu(i1,i2,j,k)
          enddo
        enddo
      enddo
      enddo
      endif
C END

      end
 
      subroutine evuii(igrad,natm,X,Y,Z,u,dudX,ii)
************************************************************************
* Subroutine to evaluate the diabatic potentials and gradients  
* The bohr-hartree-radian unit is used throughout the program
* qtc:  redundant internal coordinates
* bmqt: B-matrix of RICs
* ii:   flag to indicate the diabatic states
*       = 11, U11
*       = 22, U22
*       = 33, U33
************************************************************************

      implicit double precision (a-h,o-z)

      integer maxdim
      parameter (maxdim=100)
      
      integer ldbmt,ldbms,ldct
      parameter (ldbmt=maxdim,ldbms=maxdim,ldbnrt=maxdim,ldct=maxdim)

      integer nstc
      parameter (nstc=9)

      double precision u1min
      parameter (u1min=0.00905097d0)

C Global variables
      integer igrad
      integer natm
      double precision u
      double precision X(*),Y(*),Z(*)
      double precision dudX(3*natm)

C Local variables
      integer nqtc
      integer nblt,nbat,ntot,noopt
      data nblt/12/,nbat/18/ntot/24/,noopt/0/
      integer iblt(2,12),ibat(3,18),itot(4,24),ioopt(4,0)
      data iblt  /1,2,
     $            1,6,
     $            1,7,
     $            2,3,
     $            2,8,
     $            3,4,
     $            3,9,
     $            4,5,
     $            4,10,
     $            5,6,
     $            5,11,
     $            6,12/
      data ibat  /2,1,6,
     $            2,1,7,
     $            6,1,7,
     $            1,2,3,
     $            1,2,8,
     $            3,2,8,
     $            2,3,4,
     $            2,3,9,
     $            4,3,9,
     $            3,4,5,
     $            3,4,10,
     $            5,4,10,
     $            4,5,6,
     $            4,5,11,
     $            6,5,11,
     $            1,6,5,
     $            1,6,12,
     $            5,6,12/
      data itot  /6,1,2,3,
     $            6,1,2,8,
     $            7,1,2,3,
     $            7,1,2,8,
     $            2,1,6,5,
     $            2,1,6,12,
     $            7,1,6,5,
     $            7,1,6,12,
     $            1,2,3,4,
     $            1,2,3,9,
     $            8,2,3,4,
     $            8,2,3,9,
     $            2,3,4,5,
     $            2,3,4,10,
     $            9,3,4,5,
     $            9,3,4,10,
     $            3,4,5,6,
     $            3,4,5,11,
     $            10,4,5,6,
     $            10,4,5,11,
     $            4,5,6,1,
     $            4,5,6,12,
     $            11,5,6,1,
     $            11,5,6,12/

      integer nbls,nbas,ntos,noops
      data nbls/1/,nbas/1/,ntos/1/,noops/0/

      integer ibls(2,1),ibas(3,1),itos(4,1),ioops(4,0)
      data ibls /7,13/
      data ibas /1,7,13/
      data itos /2,1,7,13/

      double precision qtc(maxdim),qsc(maxdim)
      double precision dutdr,dubdr,dutodr,durdr
      double precision dtodph,dubdth
      double precision dut(maxdim)
      double precision bmt(ldbmt,3*natm),bms(ldbms,3*natm)

C Evaluate internal coordinates
      call evintcoord(nblt,iblt,nbat,ibat,ntot,itot,noopt,ioopt,
     $                nqtc,qtc,X,Y,Z)

      r=bndlen(7,13,X,Y,Z)
      theta=bndang(1,7,13,X,Y,Z)
      phi=dihedr(2,1,7,13,X,Y,Z)

C Evaluate Wilson's B-matrix if igrad = 1
      if (igrad .eq. 1) then
! Added by J. Zheng
C        bmt(:,:) = 0d0
        call bmat3(nblt,iblt,nbat,ibat,ntot,itot,noopt,ioopt,natm,
     $            X,Y,Z,bmt,ldbmt)
! Added by J. Zheng
C        bms(:,:) = 0d0
        call bmat3(nbls,ibls,nbas,ibas,ntos,itos,noops,ioops,natm,
     $            X,Y,Z,bms,ldbms)
      endif

C Evaluate uii from tertiary coordinates    
      call evuiit(igrad,r,nqtc,qtc,ut,dut,dutdr,ii)

C Evaluate uii from reactive (rOH) and secondary coordinates
C (bending C-O-H and torsion C-C-O-H)
      if (ii .eq. 11) then
        call evu1r(igrad,r,ur,durdr)
        call evu1b(igrad,r,theta,ub,dubdr,dubdth)
        call evu1to(igrad,r,phi,uto,dutodr,dtodph)
      else if (ii .eq. 22) then
        call evu2r(igrad,r,ur,durdr)
        call evu2b(igrad,r,theta,ub,dubdr,dubdth)
        call evu2to(igrad,r,phi,uto,dutodr,dtodph)
      else if (ii .eq. 33) then
        call evu3r(igrad,r,ur,durdr)
        call evu3b(igrad,r,theta,ub,dubdr,dubdth)
        call evu3to(igrad,r,phi,uto,dutodr,dtodph)
      endif

C Readujst zero of energy to have U1 min to be zero
      u=ur+ub+uto+ut-u1min
      dudr=durdr+dubdr+dutodr+dutdr
      dudth=dubdth
      dudph=dtodph

C Calculated analytical Cartesian gradient if igrad = 1
      if (igrad .eq. 1) then
        do i=1,3*natm
          dudX(i)=0.0d0
        enddo
C Contribution of dudX from dudr
        do i=1,3*natm
          dudX(i)=dudX(i)+dudr*bms(1,i)
        enddo
C Contribution of dudX from dudth
        do i=1,3*natm
          dudX(i)=dudX(i)+dudth*bms(2,i)
        enddo
C Contribution of dudX from dudph        
        do i=1,3*natm
          dudX(i)=dudX(i)+dudph*bms(3,i)
        enddo
C Contribution of dudX from dut
        do i=1,3*natm
          do j=1,nqtc
            dudX(i)=dudX(i)+dut(j)*bmt(j,i)
          enddo
        enddo
      endif

      end
!      subroutine evuiit(igrad,r,ntot,nqtc,qtc,ut,dut,dutdr,ii)
      subroutine evuiit(igrad,r,nqtc,qtc,ut,dut,dutdr,ii)
************************************************************************
* Subroutine to evaluate the diabatic potentials and gradients
* The bohr-hartree-radian unit is used throughout the program
* qtc:  redundant internal coordinates
* ii:   flag to indicate the diabatic states
*       = 11, U11
*       = 22, U22
*       = 33, U33
************************************************************************

      implicit none 

      integer nap
      parameter(nap=4)
C      integer nqa,nqb
C      parameter (nqa=30,nqb=24)
      integer nqs,nqb,nqt
C      parameter (nqs=12,nqb=18,nqt=24)
      parameter (nqs=12,nqb=18,nqt=24)

C Global variables
      integer igrad,ii
      integer ntot,nqtc
      double precision r,ut,dutdr
      double precision qtc(*),dut(*)

C Local variables
      integer i,j,k
      integer ir
      double precision u0,u01,u02,u1,u2
      double precision vs,vs1,vs2
      double precision vb,vb1,vb2
      double precision vt,vt1,vt2
      double precision vbs,vbs1,vbs2
      double precision dvsdr,dvbdr,dvtdr,dvbsdr
      double precision t(nap),dtdr(nap)
      double precision rap(nap)
      double precision v0(nap)
      double precision qs(nqs),qs1(nqs),qs2(nqs)
      double precision qb(nqb),qb1(nqb),qb2(nqb)
      double precision sqb(nqb)
      double precision qt(nqt),qt1(nqt),qt2(nqt)
      double precision sqt1(nqt),cqt1(nqt)
      double precision sqt2(nqt),cqt2(nqt)
      double precision dvs(nqs),dvs1(nqs),dvs2(nqs)
      double precision dvb(nqb),dvb1(nqb),dvb2(nqb)
      double precision dvt(nqt),dvt1(nqt),dvt2(nqt)
      double precision dvbsdb(nqb),dvbsdb1(nqb),dvbsdb2(nqb)
      double precision dvbsds(nqs),dvbsds1(nqs),dvbsds2(nqs)
      double precision qse(nqs,nap)
      double precision qbe(nqb,nap)
      double precision qte(nqt,nap)
      double precision qse1(nqs),qse2(nqs)
      double precision qbe1(nqb),qbe2(nqb)
      double precision FS1(nqs,nqs),FS2(nqs,nqs)
      double precision FB1(nqb,nqb),FB2(nqb,nqb)
      double precision FT1(nqt,nqt),FT2(nqt,nqt)
      double precision FBS1(nqb,nqs),FBS2(nqb,nqs)
      double precision HS(nqs,nqs,nap)
      double precision HB(nqb,nqb,nap)
      double precision HT(nqt,nqt,nap)
      double precision HBS(nqb,nqs,nap)

C Initialize output variables
      ut=0.0d0
      dutdr=0.0d0
      do i=1,nqtc
        dut(i)=0.0d0
      enddo

C Get HS for stretchings, HB for bendings, and HT for torsions
      call gethess(nap,HS,nqs,HB,nqb,HT,nqt,HBS,ii)
C Get equilibrium parameters 
      call getgeom(nap,qse,nqs,qbe,nqb,qte,nqt,ii)

      if (ii .eq. 11) then
        rap(1)=1.93206106d0 ! rOH=1.0224 using rOH=0.964A
        rap(2)=2.49443845d0 ! rOH=1.32 A
        rap(3)=3.77945219d0 ! rOH=2.00 A
        rap(4)=9.44863047d0 ! rOH=5.00 A
C Normalized relaxation energy 
C dE = E(DFT_opt) - E(CAS_geom)
        v0(1)=0.0d0         
        v0(2)=-0.000346019d0
        v0(3)=0.000672861d0
        v0(4)=-0.001758495d0
      else if (ii .eq. 22) then
        rap(1)=1.93206106d0 ! rOH=1.0224 using rOH=0.964A
        rap(2)=2.49443845d0 ! rOH=1.32 A
        rap(3)=4.27078097d0 ! rOH=2.26 A
        rap(4)=9.44863047d0 ! rOH=5.00 A
        v0(1)=-0.008940821d0
        v0(2)=-0.00976278d0
        v0(3)=-0.00976278d0
        v0(4)=-0.00976278d0
      else if (ii .eq. 33) then
        rap(1)=1.93206106d0 ! rOH=1.0224 using rOH=0.964A
        rap(2)=2.49443845d0 ! rOH=1.32 A
        rap(3)=4.27078097d0 ! rOH=2.26 A
        rap(4)=9.44863047d0 ! rOH=5.00 A
        v0(1)=-0.012423316d0
        v0(2)=-0.015044491d0
        v0(3)=-0.014671263d0
        v0(4)=-0.01433271d0
      endif

      call tent(nap,rap,r,ir,t,dtdr)

      do i=1,nqs
        qs(i)=qtc(i)
      enddo
      do i=1,nqb
        qb(i)=qtc(i+nqs)
        sqb(i)=sin(qb(i))
      enddo
      do i=1,nqt
        qt(i)=qtc(i+nqs+nqb)
      enddo

C Initialize hessian elements of stretchings, bendings, and torsions
      do i=1,nqs
        do j=1,i
          FS1(i,j)=0.0d0
          FS2(i,j)=0.0d0
        enddo
      enddo
      do i=1,nqb
        do j=1,i
          FB1(i,j)=0.0d0
          FB2(i,j)=0.0d0
        enddo
      enddo
      do i=1,nqt
        do j=1,i
          FT1(i,j)=0.0d0
          FT2(i,j)=0.0d0
        enddo
      enddo
      do i=1,nqb
        do j=1,nqs
          FBS1(i,j)=0.0d0
          FBS2(i,j)=0.0d0
        enddo
      enddo

C Get force constant matrices for strecthings (FS) and bending (FB)
C and torsions (FT)
      if (ir .eq. 1) then
        u0=v0(1)
        do i=1,nqs
C          qs1(i)=qs(i)-qse(i,1)
          qse1(i)=qse(i,1)
          do j=1,i
            FS1(i,j)=HS(i,j,1)*qse1(i)*qse1(j)
          enddo
        enddo
        do i=1,nqb
C          qb1(i)=qb(i)-qbe(i,1)
C cos(theta0)-cos(theta) as variables
          qbe1(i)=qbe(i,1)
          qb1(i)=cos(qbe1(i))-cos(qb(i))
          do j=1,i
            FB1(i,j)=HB(i,j,1)/(sin(qbe1(i))*sin(qbe1(j)))
          enddo
        enddo
        do i=1,nqt
          qt1(i)=qt(i)-qte(i,1)
          sqt1(i)=sin(qt1(i))
          cqt1(i)=cos(qt1(i))
          do j=1,i
            FT1(i,j)=HT(i,j,1)
          enddo
        enddo
C bend-stretch couplings
        do i=1,nqb
          do j=1,nqs
            FBS1(i,j)=HBS(i,j,1)*qse1(j)/sin(qbe1(i))
          enddo
        enddo
      else if (ir .eq. nap+1) then
        u0=v0(nap)
        do i=1,nqs
          qse1(i)=qse(i,nap)
          do j=1,i
            FS1(i,j)=HS(i,j,nap)*qse1(i)*qse1(j)
          enddo
        enddo
        do i=1,nqb
C cos(theta0)-cos(theta) as variables
          qbe1(i)=qbe(i,nap)
          qb1(i)=cos(qbe1(i))-cos(qb(i))
          do j=1,i
            FB1(i,j)=HB(i,j,nap)/(sin(qbe1(i))*sin(qbe1(j)))
          enddo
        enddo
        do i=1,nqt
          qt1(i)=qt(i)-qte(i,nap)
          sqt1(i)=sin(qt1(i))
          cqt1(i)=cos(qt1(i))
          do j=1,i
            FT1(i,j)=HT(i,j,nap)
          enddo
        enddo
C bend-stretch couplings
        do i=1,nqb
          do j=1,nqs
            FBS1(i,j)=HBS(i,j,nap)*qse1(j)/sin(qbe1(i))
          enddo
        enddo
      else if ((ir .gt. 1) .and. (ir .lt. nap+1)) then
        u01=v0(ir-1)
        u02=v0(ir)
        do i=1,nqs
          qse1(i)=qse(i,ir-1)
          qse2(i)=qse(i,ir)
          do j=1,i
            FS1(i,j)=HS(i,j,ir-1)*qse1(i)*qse1(j)
            FS2(i,j)=HS(i,j,ir)*qse2(i)*qse2(j)
          enddo
        enddo
        do i=1,nqb
          qbe1(i)=qbe(i,ir-1)
          qbe2(i)=qbe(i,ir)
          qb1(i)=cos(qbe1(i))-cos(qb(i))
          qb2(i)=cos(qbe2(i))-cos(qb(i))
          do j=1,i
            FB1(i,j)=HB(i,j,ir-1)/(sin(qbe1(i))*sin(qbe1(j)))
            FB2(i,j)=HB(i,j,ir)/(sin(qbe2(i))*sin(qbe2(j)))
          enddo
        enddo
        do i=1,nqt
          qt1(i)=qt(i)-qte(i,ir-1)
          qt2(i)=qt(i)-qte(i,ir)
          sqt1(i)=sin(qt1(i))
          sqt2(i)=sin(qt2(i))
          cqt1(i)=cos(qt1(i))
          cqt2(i)=cos(qt2(i))
          do j=1,i
            FT1(i,j)=HT(i,j,ir-1)
            FT2(i,j)=HT(i,j,ir)
          enddo
        enddo
C bend-stretch couplings
        do i=1,nqb
          do j=1,nqs
            FBS1(i,j)=HBS(i,j,ir-1)*qse1(j)/sin(qbe1(i))
            FBS2(i,j)=HBS(i,j,ir)*qse2(j)/sin(qbe2(i))
          enddo
        enddo
      endif

C Get the upper parts of force matrix by symmetry
      do i=1,nqs
        do j=i+1,nqs
          FS1(i,j)=FS1(j,i)
          FS2(i,j)=FS2(j,i)
        enddo
      enddo
      do i=1,nqb
        do j=i+1,nqb
          FB1(i,j)=FB1(j,i)
          FB2(i,j)=FB2(j,i)
        enddo
      enddo
      do i=1,nqt
        do j=i+1,nqt
          FT1(i,j)=FT1(j,i)
          FT2(i,j)=FT2(j,i)
        enddo
      enddo

C      write(*,'(A,I2)') "ir = ",ir
C Calculate potential energy
      if ((ir .eq. 1) .or. (ir .eq. nap+1)) then
        vs=0.0d0
        vb=0.0d0
        vt=0.0d0
        vbs=0.0d0
        do i=1,nqs
          do j=1,nqs
C            vs=vs+0.5d0*FS1(i,j)*(qs(i)-qse1(i))*(qs(j)-qse1(j))
            vs=vs+0.5d0*FS1(i,j)
     $           *(qs(i)-qse1(i))*(qs(j)-qse1(j))/(qs(i)*qs(j))
          enddo
        enddo
        do i=1,nqb
          do j=1,nqb
            vb=vb+0.5d0*qb1(i)*FB1(i,j)*qb1(j)
          enddo
        enddo
        do i=1,nqt
          do j=1,nqt
            vt=vt+0.5d0*FT1(i,j)*sqt1(i)*sqt1(j)
          enddo
          vt=vt-0.5d0*FT1(i,i)*sqt1(i)*sqt1(i)
     $      +2.0d0*FT1(i,i)*sin(qt1(i)/2.0d0)*sin(qt1(i)/2.0d0)
        enddo
C bend-stretch couplings
        do i=1,nqb
          do j=1,nqs
C            vbs=vbs+FBS1(i,j)*qb1(i)*(qs(j)-qse1(j))
            vbs=vbs+FBS1(i,j)*qb1(i)*(qs(j)-qse1(j))/qs(j)
          enddo
        enddo

        ut=u0+vs+vb+vt+vbs
        if (igrad .eq. 1) then
C          dvsdr=0.0d0
C          dvbdr=0.0d0
C          dvtdr=0.0d0
C          dvbsdr=0.0d0
          do i=1,nqs
            dvs(i)=0.0d0
            do j=1,nqs
C              dvs(i)=dvs(i)+FS1(i,j)*(qs(j)-qse1(j))
              dvs(i)=dvs(i)+FS1(i,j)*qse1(i)
     $                     *(qs(j)-qse1(j))/(qs(i)*qs(i)*qs(j))
            enddo
          enddo
          do i=1,nqb
            dvb(i)=0.0d0
            do j=1,nqb
              dvb(i)=dvb(i)+sqb(i)*FB1(i,j)*qb1(j)
            enddo
          enddo
          do i=1,nqt
            dvt(i)=0.0d0
            do j=1,nqt
              dvt(i)=dvt(i)+FT1(i,j)*cqt1(i)*sqt1(j)
            enddo
            dvt(i)=dvt(i)-FT1(i,i)*cqt1(i)*sqt1(i)
     $            +2.0d0*FT1(i,i)*cos(qt1(i)/2.0d0)*sin(qt1(i)/2.0d0)
          enddo
C contribution from vbs (dvs and dvb)
          do i=1,nqb
            dvbsdb(i)=0.0d0
            do j=1,nqs
C              dvbsdb(i)=dvbsdb(i)+FBS1(i,j)*sqb(i)*(qs(j)-qse1(j))
              dvbsdb(i)=dvbsdb(i)+FBS1(i,j)*sqb(i)
     $                           *(qs(j)-qse1(j))/qs(j)
            enddo
          enddo
          do i=1,nqs
            dvbsds(i)=0.0d0
            do j=1,nqb
C              dvbsds(i)=dvbsds(i)+FBS1(j,i)*qb1(j)
              dvbsds(i)=dvbsds(i)+FBS1(j,i)*qse1(i)
     $                           *qb1(j)/(qs(i)*qs(i))
            enddo
          enddo

C          dutdr=dvsdr+dvbdr+dvtdr+dvbsdr
          dutdr=0.0d0
          do i=1,nqs
            dut(i)=dvs(i)+dvbsds(i)
          enddo
          do i=1,nqb
            dut(nqs+i)=dvb(i)+dvbsdb(i)
          enddo
          do i=1,nqt
            dut(nqs+nqb+i)=dvt(i)
          enddo
        endif
      else
        vs1=0.0d0
        vs2=0.0d0
        vb1=0.0d0
        vb2=0.0d0
        vt1=0.0d0
        vt2=0.0d0
        vbs1=0.0d0
        vbs2=0.0d0
        do i=1,nqs
          do j=1,nqs
C            vs1=vs1+0.5d0*FS1(i,j)*(qs(i)-qse1(i))*(qs(j)-qse1(j))
            vs1=vs1+0.5d0*FS1(i,j)
     $             *(qs(i)-qse1(i))*(qs(j)-qse1(j))/(qs(i)*qs(j))
C            vs2=vs2+0.5d0*FS2(i,j)*(qs(i)-qse2(i))*(qs(j)-qse2(j))
            vs2=vs2+0.5d0*FS2(i,j)
     $             *(qs(i)-qse2(i))*(qs(j)-qse2(j))/(qs(i)*qs(j))
          enddo
        enddo
        do i=1,nqb
          do j=1,nqb
            vb1=vb1+0.5d0*qb1(i)*FB1(i,j)*qb1(j)
            vb2=vb2+0.5d0*qb2(i)*FB2(i,j)*qb2(j)
          enddo
        enddo
        do i=1,nqt
          do j=1,nqt
            vt1=vt1+0.5d0*sqt1(i)*FT1(i,j)*sqt1(j)
            vt2=vt2+0.5d0*sqt2(i)*FT2(i,j)*sqt2(j)
          enddo
          vt1=vt1-0.5d0*FT1(i,i)*sqt1(i)*sqt1(i)
     $           +2.0d0*FT1(i,i)*sin(qt1(i)/2.0d0)*sin(qt1(i)/2.0d0)
          vt2=vt2-0.5d0*FT2(i,i)*sqt2(i)*sqt2(i)
     $           +2.0d0*FT2(i,i)*sin(qt2(i)/2.0d0)*sin(qt2(i)/2.0d0)
        enddo
C bend-stretch couplings
        do i=1,nqb
          do j=1,nqs
C            vbs1=vbs1+qb1(i)*FBS1(i,j)*(qs(j)-qse1(j))
C            vbs2=vbs2+qb2(i)*FBS2(i,j)*(qs(j)-qse2(j))
            vbs1=vbs1+FBS1(i,j)*qb1(i)*(qs(j)-qse1(j))/qs(j)
            vbs2=vbs2+FBS2(i,j)*qb2(i)*(qs(j)-qse2(j))/qs(j)
          enddo
        enddo

        u1=u01+vs1+vb1+vt1+vbs1
        u2=u02+vs2+vb2+vt2+vbs2
        ut=u1*t(ir-1)+u2*t(ir)
        if (igrad .eq. 1) then
          dutdr=u1*dtdr(ir-1)+u2*dtdr(ir)
          do i=1,nqs
            dvs1(i)=0.0d0
            dvs2(i)=0.0d0
            do j=1,nqs
C              dvs1(i)=dvs1(i)+FS1(i,j)*(qs(j)-qse1(j))
C              dvs2(i)=dvs2(i)+FS2(i,j)*(qs(j)-qse2(j))
              dvs1(i)=dvs1(i)+FS1(i,j)*qse1(i)
     $                     *(qs(j)-qse1(j))/(qs(i)*qs(i)*qs(j))
              dvs2(i)=dvs2(i)+FS2(i,j)*qse2(i)
     $                     *(qs(j)-qse2(j))/(qs(i)*qs(i)*qs(j))
            enddo
            dvs(i)=dvs1(i)*t(ir-1)+dvs2(i)*t(ir)
          enddo
          do i=1,nqb
            dvb1(i)=0.0d0
            dvb2(i)=0.0d0
            do j=1,nqb
              dvb1(i)=dvb1(i)+sqb(i)*FB1(i,j)*qb1(j)
              dvb2(i)=dvb2(i)+sqb(i)*FB2(i,j)*qb2(j)
            enddo
            dvb(i)=dvb1(i)*t(ir-1)+dvb2(i)*t(ir)
          enddo
          do i=1,nqt
            dvt1(i)=0.0d0
            dvt2(i)=0.0d0
            do j=1,nqt
              dvt1(i)=dvt1(i)+FT1(i,j)*cqt1(i)*sqt1(j)
              dvt2(i)=dvt2(i)+FT2(i,j)*cqt2(i)*sqt2(j)
            enddo
            dvt1(i)=dvt1(i)-FT1(i,i)*cqt1(i)*sqt1(i)
     $             +2.0d0*FT1(i,i)*sin(qt1(i)/2.0d0)*cos(qt1(i)/2.0d0)
            dvt2(i)=dvt2(i)-FT2(i,i)*cqt2(i)*sqt2(i)
     $             +2.0d0*FT2(i,i)*sin(qt2(i)/2.0d0)*cos(qt2(i)/2.0d0)
            dvt(i)=dvt1(i)*t(ir-1)+dvt2(i)*t(ir)
          enddo
C contribution from vbs (dvbsdb and dvbsds)
          do i=1,nqb
            dvbsdb1(i)=0.0d0
            dvbsdb2(i)=0.0d0
            do j=1,nqs
C              dvbsdb1(i)=dvbsdb1(i)+FBS1(i,j)*sqb(i)*(qs(j)-qse1(j))
C              dvbsdb2(i)=dvbsdb2(i)+FBS2(i,j)*sqb(i)*(qs(j)-qse2(j))
              dvbsdb1(i)=dvbsdb1(i)+FBS1(i,j)*sqb(i)
     $                           *(qs(j)-qse1(j))/qs(j)
              dvbsdb2(i)=dvbsdb2(i)+FBS2(i,j)*sqb(i)
     $                           *(qs(j)-qse2(j))/qs(j)
            enddo
            dvbsdb(i)=dvbsdb1(i)*t(ir-1)+dvbsdb2(i)*t(ir)
          enddo
          do i=1,nqs
            dvbsds1(i)=0.0d0
            dvbsds2(i)=0.0d0
            do j=1,nqb
C              dvbsds1(i)=dvbsds1(i)+FBS1(j,i)*qb1(j)
C              dvbsds2(i)=dvbsds2(i)+FBS2(j,i)*qb2(j)
              dvbsds1(i)=dvbsds1(i)+FBS1(j,i)*qse1(i)
     $                           *qb1(j)/(qs(i)*qs(i))
              dvbsds2(i)=dvbsds2(i)+FBS2(j,i)*qse2(i)
     $                           *qb2(j)/(qs(i)*qs(i))
            enddo
            dvbsds(i)=dvbsds1(i)*t(ir-1)+dvbsds2(i)*t(ir)
          enddo

          do i=1,nqs
            dut(i)=dvs(i)+dvbsds(i)
          enddo
          do i=1,nqb
            dut(nqs+i)=dvb(i)+dvbsdb(i)
          enddo
          do i=1,nqt
            dut(nqs+nqb+i)=dvt(i)
          enddo
        endif
      endif

      end
      subroutine gethess(nap,HSS,nqs,HBB,nqb,HTT,nqt,HBS,ii)

C Global variables
      integer ii
      integer nqs,nqb,nqt
      double precision HSS(nqs,nqs,nap)
      double precision HBB(nqb,nqb,nap)
      double precision HTT(nqt,nqt,nap)
      double precision HBS(nqb,nqs,nap)

      if (ii .eq. 11) then
C Hessian elements for anchor point 1
        HSS( 1, 1, 1) =  0.40759D+00
        HSS( 2, 1, 1) =  0.17500D-01
        HSS( 2, 2, 1) =  0.40901D+00
        HSS( 3, 1, 1) =  0.29920D-01
        HSS( 3, 2, 1) =  0.35360D-01
        HSS( 3, 3, 1) =  0.42477D+00
        HSS( 4, 1, 1) =  0.18000D-01
        HSS( 4, 2, 1) = -0.26000D-02
        HSS( 4, 3, 1) =  0.38000D-03
        HSS( 4, 4, 1) =  0.41176D+00
        HSS( 5, 1, 1) =  0.45000D-02
        HSS( 5, 2, 1) = -0.29000D-03
        HSS( 5, 3, 1) =  0.30900D-02
        HSS( 5, 4, 1) =  0.51800D-02
        HSS( 5, 5, 1) =  0.35919D+00
        HSS( 6, 1, 1) = -0.34400D-02
        HSS( 6, 2, 1) =  0.82770D-01
        HSS( 6, 3, 1) = -0.25600D-02
        HSS( 6, 4, 1) =  0.21010D-01
        HSS( 6, 5, 1) = -0.76000D-03
        HSS( 6, 6, 1) =  0.41463D+00
        HSS( 7, 1, 1) = -0.55000D-03
        HSS( 7, 2, 1) = -0.34000D-03
        HSS( 7, 3, 1) =  0.68000D-03
        HSS( 7, 4, 1) =  0.44000D-02
        HSS( 7, 5, 1) =  0.48000D-03
        HSS( 7, 6, 1) =  0.45900D-02
        HSS( 7, 7, 1) =  0.36428D+00
        HSS( 8, 1, 1) =  0.82220D-01
        HSS( 8, 2, 1) = -0.40500D-02
        HSS( 8, 3, 1) = -0.11000D-03
        HSS( 8, 4, 1) = -0.33000D-02
        HSS( 8, 5, 1) = -0.30000D-03
        HSS( 8, 6, 1) =  0.20930D-01
        HSS( 8, 7, 1) = -0.31000D-03
        HSS( 8, 8, 1) =  0.40778D+00
        HSS( 9, 1, 1) = -0.48000D-03
        HSS( 9, 2, 1) = -0.49000D-03
        HSS( 9, 3, 1) = -0.23000D-03
        HSS( 9, 4, 1) = -0.57000D-03
        HSS( 9, 5, 1) =  0.20000D-03
        HSS( 9, 6, 1) =  0.44500D-02
        HSS( 9, 7, 1) =  0.47000D-03
        HSS( 9, 8, 1) =  0.44300D-02
        HSS( 9, 9, 1) =  0.36674D+00
        HSS(10, 1, 1) = -0.18900D-02
        HSS(10, 2, 1) =  0.18160D-01
        HSS(10, 3, 1) = -0.10000D-02
        HSS(10, 4, 1) =  0.86110D-01
        HSS(10, 5, 1) = -0.74000D-03
        HSS(10, 6, 1) = -0.29100D-02
        HSS(10, 7, 1) = -0.54000D-03
        HSS(10, 8, 1) =  0.20960D-01
        HSS(10, 9, 1) = -0.51000D-03
        HSS(10,10, 1) =  0.41852D+00
        HSS(11, 1, 1) = -0.43000D-03
        HSS(11, 2, 1) = -0.50000D-03
        HSS(11, 3, 1) =  0.78000D-03
        HSS(11, 4, 1) = -0.50000D-03
        HSS(11, 5, 1) =  0.11000D-03
        HSS(11, 6, 1) = -0.28000D-03
        HSS(11, 7, 1) =  0.20000D-03
        HSS(11, 8, 1) =  0.45700D-02
        HSS(11, 9, 1) =  0.47000D-03
        HSS(11,10, 1) =  0.45000D-02
        HSS(11,11, 1) =  0.36422D+00
        HSS(12, 1, 1) = -0.40000D-03
        HSS(12, 2, 1) =  0.35500D-02
        HSS(12, 3, 1) =  0.10000D-04
        HSS(12, 4, 1) = -0.80000D-03
        HSS(12, 5, 1) =  0.23000D-03
        HSS(12, 6, 1) = -0.14000D-03
        HSS(12, 7, 1) =  0.10000D-03
        HSS(12, 8, 1) = -0.68000D-03
        HSS(12, 9, 1) =  0.14000D-03
        HSS(12,10, 1) =  0.42300D-02
        HSS(12,11, 1) =  0.44000D-03
        HSS(12,12, 1) =  0.36727D+00
C Hessian elements for anchor point 1
        HBB( 1, 1, 1) =  0.86570D-01
        HBB( 2, 1, 1) = -0.44420D-01
        HBB( 2, 2, 1) =  0.13902D+00
        HBB( 3, 1, 1) = -0.42150D-01
        HBB( 3, 2, 1) = -0.94600D-01
        HBB( 3, 3, 1) =  0.13674D+00
        HBB( 4, 1, 1) = -0.36220D-01
        HBB( 4, 2, 1) =  0.23510D-01
        HBB( 4, 3, 1) =  0.12710D-01
        HBB( 4, 4, 1) =  0.79210D-01
        HBB( 5, 1, 1) =  0.20780D-01
        HBB( 5, 2, 1) = -0.13590D-01
        HBB( 5, 3, 1) = -0.71800D-02
        HBB( 5, 4, 1) = -0.38760D-01
        HBB( 5, 5, 1) =  0.77300D-01
        HBB( 6, 1, 1) =  0.15440D-01
        HBB( 6, 2, 1) = -0.99200D-02
        HBB( 6, 3, 1) = -0.55200D-02
        HBB( 6, 4, 1) = -0.40440D-01
        HBB( 6, 5, 1) = -0.38530D-01
        HBB( 6, 6, 1) =  0.78980D-01
        HBB( 7, 1, 1) = -0.10540D-01
        HBB( 7, 2, 1) = -0.24300D-02
        HBB( 7, 3, 1) =  0.12970D-01
        HBB( 7, 4, 1) = -0.31280D-01
        HBB( 7, 5, 1) =  0.12070D-01
        HBB( 7, 6, 1) =  0.19210D-01
        HBB( 7, 7, 1) =  0.79110D-01
        HBB( 8, 1, 1) = -0.42000D-03
        HBB( 8, 2, 1) =  0.59300D-02
        HBB( 8, 3, 1) = -0.55100D-02
        HBB( 8, 4, 1) =  0.18520D-01
        HBB( 8, 5, 1) = -0.65400D-02
        HBB( 8, 6, 1) = -0.11980D-01
        HBB( 8, 7, 1) = -0.39530D-01
        HBB( 8, 8, 1) =  0.79370D-01
        HBB( 9, 1, 1) =  0.10960D-01
        HBB( 9, 2, 1) = -0.35000D-02
        HBB( 9, 3, 1) = -0.74600D-02
        HBB( 9, 4, 1) =  0.12760D-01
        HBB( 9, 5, 1) = -0.55300D-02
        HBB( 9, 6, 1) = -0.72300D-02
        HBB( 9, 7, 1) = -0.39590D-01
        HBB( 9, 8, 1) = -0.39850D-01
        HBB( 9, 9, 1) =  0.79440D-01
        HBB(10, 1, 1) =  0.69200D-02
        HBB(10, 2, 1) = -0.37400D-02
        HBB(10, 3, 1) = -0.31900D-02
        HBB(10, 4, 1) = -0.97700D-02
        HBB(10, 5, 1) =  0.10500D-01
        HBB(10, 6, 1) = -0.72000D-03
        HBB(10, 7, 1) = -0.33920D-01
        HBB(10, 8, 1) =  0.13770D-01
        HBB(10, 9, 1) =  0.20160D-01
        HBB(10,10, 1) =  0.80580D-01
        HBB(11, 1, 1) = -0.33900D-02
        HBB(11, 2, 1) =  0.23200D-02
        HBB(11, 3, 1) =  0.10700D-02
        HBB(11, 4, 1) = -0.78000D-03
        HBB(11, 5, 1) = -0.37400D-02
        HBB(11, 6, 1) =  0.45200D-02
        HBB(11, 7, 1) =  0.19900D-01
        HBB(11, 8, 1) = -0.73600D-02
        HBB(11, 9, 1) = -0.12550D-01
        HBB(11,10, 1) = -0.40290D-01
        HBB(11,11, 1) =  0.79280D-01
        HBB(12, 1, 1) = -0.35300D-02
        HBB(12, 2, 1) =  0.14200D-02
        HBB(12, 3, 1) =  0.21200D-02
        HBB(12, 4, 1) =  0.10550D-01
        HBB(12, 5, 1) = -0.67600D-02
        HBB(12, 6, 1) = -0.37900D-02
        HBB(12, 7, 1) =  0.14020D-01
        HBB(12, 8, 1) = -0.64100D-02
        HBB(12, 9, 1) = -0.76100D-02
        HBB(12,10, 1) = -0.40290D-01
        HBB(12,11, 1) = -0.38990D-01
        HBB(12,12, 1) =  0.79280D-01
        HBB(13, 1, 1) = -0.10320D-01
        HBB(13, 2, 1) =  0.13110D-01
        HBB(13, 3, 1) = -0.27900D-02
        HBB(13, 4, 1) =  0.52200D-02
        HBB(13, 5, 1) = -0.23700D-02
        HBB(13, 6, 1) = -0.28500D-02
        HBB(13, 7, 1) = -0.86400D-02
        HBB(13, 8, 1) =  0.10140D-01
        HBB(13, 9, 1) = -0.15000D-02
        HBB(13,10, 1) = -0.34100D-01
        HBB(13,11, 1) =  0.14060D-01
        HBB(13,12, 1) =  0.20050D-01
        HBB(13,13, 1) =  0.79730D-01
        HBB(14, 1, 1) =  0.10860D-01
        HBB(14, 2, 1) = -0.76500D-02
        HBB(14, 3, 1) = -0.32000D-02
        HBB(14, 4, 1) = -0.26900D-02
        HBB(14, 5, 1) =  0.11000D-02
        HBB(14, 6, 1) =  0.16000D-02
        HBB(14, 7, 1) = -0.15200D-02
        HBB(14, 8, 1) = -0.35000D-02
        HBB(14, 9, 1) =  0.50200D-02
        HBB(14,10, 1) =  0.20210D-01
        HBB(14,11, 1) = -0.76100D-02
        HBB(14,12, 1) = -0.12600D-01
        HBB(14,13, 1) = -0.39850D-01
        HBB(14,14, 1) =  0.79680D-01
        HBB(15, 1, 1) = -0.54000D-03
        HBB(15, 2, 1) = -0.54500D-02
        HBB(15, 3, 1) =  0.59900D-02
        HBB(15, 4, 1) = -0.25300D-02
        HBB(15, 5, 1) =  0.12700D-02
        HBB(15, 6, 1) =  0.12600D-02
        HBB(15, 7, 1) =  0.10150D-01
        HBB(15, 8, 1) = -0.66300D-02
        HBB(15, 9, 1) = -0.35200D-02
        HBB(15,10, 1) =  0.13890D-01
        HBB(15,11, 1) = -0.64500D-02
        HBB(15,12, 1) = -0.74400D-02
        HBB(15,13, 1) = -0.39880D-01
        HBB(15,14, 1) = -0.39830D-01
        HBB(15,15, 1) =  0.79710D-01
        HBB(16, 1, 1) = -0.36420D-01
        HBB(16, 2, 1) =  0.13970D-01
        HBB(16, 3, 1) =  0.22450D-01
        HBB(16, 4, 1) = -0.71600D-02
        HBB(16, 5, 1) = -0.22100D-02
        HBB(16, 6, 1) =  0.93700D-02
        HBB(16, 7, 1) =  0.52700D-02
        HBB(16, 8, 1) = -0.24800D-02
        HBB(16, 9, 1) = -0.27900D-02
        HBB(16,10, 1) = -0.97100D-02
        HBB(16,11, 1) =  0.10500D-01
        HBB(16,12, 1) = -0.79000D-03
        HBB(16,13, 1) = -0.31890D-01
        HBB(16,14, 1) =  0.12990D-01
        HBB(16,15, 1) =  0.18900D-01
        HBB(16,16, 1) =  0.79900D-01
        HBB(17, 1, 1) =  0.22220D-01
        HBB(17, 2, 1) = -0.60200D-02
        HBB(17, 3, 1) = -0.16200D-01
        HBB(17, 4, 1) = -0.27800D-02
        HBB(17, 5, 1) =  0.56900D-02
        HBB(17, 6, 1) = -0.29000D-02
        HBB(17, 7, 1) = -0.25500D-02
        HBB(17, 8, 1) =  0.13300D-02
        HBB(17, 9, 1) =  0.12200D-02
        HBB(17,10, 1) =  0.10510D-01
        HBB(17,11, 1) = -0.68200D-02
        HBB(17,12, 1) = -0.37000D-02
        HBB(17,13, 1) =  0.12370D-01
        HBB(17,14, 1) = -0.57100D-02
        HBB(17,15, 1) = -0.66600D-02
        HBB(17,16, 1) = -0.39780D-01
        HBB(17,17, 1) =  0.76100D-01
        HBB(18, 1, 1) =  0.14200D-01
        HBB(18, 2, 1) = -0.79500D-02
        HBB(18, 3, 1) = -0.62400D-02
        HBB(18, 4, 1) =  0.99400D-02
        HBB(18, 5, 1) = -0.34800D-02
        HBB(18, 6, 1) = -0.64600D-02
        HBB(18, 7, 1) = -0.27200D-02
        HBB(18, 8, 1) =  0.11500D-02
        HBB(18, 9, 1) =  0.15700D-02
        HBB(18,10, 1) = -0.80000D-03
        HBB(18,11, 1) = -0.36800D-02
        HBB(18,12, 1) =  0.44900D-02
        HBB(18,13, 1) =  0.19510D-01
        HBB(18,14, 1) = -0.72700D-02
        HBB(18,15, 1) = -0.12240D-01
        HBB(18,16, 1) = -0.40120D-01
        HBB(18,17, 1) = -0.36330D-01
        HBB(18,18, 1) =  0.76450D-01
C Hessian elements for anchor point 1
        HTT( 1, 1, 1) =  0.92400D-02
        HTT( 2, 1, 1) =  0.38700D-02
        HTT( 2, 2, 1) =  0.84000D-02
        HTT( 3, 1, 1) =  0.67000D-03
        HTT( 3, 2, 1) = -0.40800D-02
        HTT( 3, 3, 1) =  0.11550D-01
        HTT( 4, 1, 1) = -0.47000D-02
        HTT( 4, 2, 1) =  0.44000D-03
        HTT( 4, 3, 1) =  0.68000D-02
        HTT( 4, 4, 1) =  0.11950D-01
        HTT( 5, 1, 1) = -0.62200D-02
        HTT( 5, 2, 1) = -0.44300D-02
        HTT( 5, 3, 1) =  0.20700D-02
        HTT( 5, 4, 1) =  0.38700D-02
        HTT( 5, 5, 1) =  0.89900D-02
        HTT( 6, 1, 1) = -0.43000D-02
        HTT( 6, 2, 1) = -0.29300D-02
        HTT( 6, 3, 1) =  0.28800D-02
        HTT( 6, 4, 1) =  0.42500D-02
        HTT( 6, 5, 1) =  0.35700D-02
        HTT( 6, 6, 1) =  0.80900D-02
        HTT( 7, 1, 1) =  0.19000D-02
        HTT( 7, 2, 1) =  0.31100D-02
        HTT( 7, 3, 1) = -0.82500D-02
        HTT( 7, 4, 1) = -0.70400D-02
        HTT( 7, 5, 1) =  0.11300D-02
        HTT( 7, 6, 1) = -0.32500D-02
        HTT( 7, 7, 1) =  0.10750D-01
        HTT( 8, 1, 1) =  0.38300D-02
        HTT( 8, 2, 1) =  0.46100D-02
        HTT( 8, 3, 1) = -0.74400D-02
        HTT( 8, 4, 1) = -0.66600D-02
        HTT( 8, 5, 1) = -0.43000D-02
        HTT( 8, 6, 1) =  0.12800D-02
        HTT( 8, 7, 1) =  0.63800D-02
        HTT( 8, 8, 1) =  0.11960D-01
        HTT( 9, 1, 1) = -0.63700D-02
        HTT( 9, 2, 1) = -0.83000D-03
        HTT( 9, 3, 1) = -0.30300D-02
        HTT( 9, 4, 1) =  0.25100D-02
        HTT( 9, 5, 1) =  0.48000D-03
        HTT( 9, 6, 1) =  0.20600D-02
        HTT( 9, 7, 1) = -0.26900D-02
        HTT( 9, 8, 1) = -0.11100D-02
        HTT( 9, 9, 1) =  0.94700D-02
        HTT(10, 1, 1) = -0.40200D-02
        HTT(10, 2, 1) =  0.76000D-03
        HTT(10, 3, 1) = -0.16700D-02
        HTT(10, 4, 1) =  0.31200D-02
        HTT(10, 5, 1) =  0.20200D-02
        HTT(10, 6, 1) =  0.15300D-02
        HTT(10, 7, 1) = -0.22000D-03
        HTT(10, 8, 1) = -0.70000D-03
        HTT(10, 9, 1) =  0.31100D-02
        HTT(10,10, 1) =  0.83700D-02
        HTT(11, 1, 1) = -0.98000D-03
        HTT(11, 2, 1) = -0.53700D-02
        HTT(11, 3, 1) =  0.17300D-02
        HTT(11, 4, 1) = -0.26500D-02
        HTT(11, 5, 1) = -0.13300D-02
        HTT(11, 6, 1) =  0.68000D-03
        HTT(11, 7, 1) = -0.39000D-02
        HTT(11, 8, 1) = -0.19000D-02
        HTT(11, 9, 1) =  0.39100D-02
        HTT(11,10, 1) = -0.16900D-02
        HTT(11,11, 1) =  0.83100D-02
        HTT(12, 1, 1) =  0.13700D-02
        HTT(12, 2, 1) = -0.37800D-02
        HTT(12, 3, 1) =  0.31000D-02
        HTT(12, 4, 1) = -0.20400D-02
        HTT(12, 5, 1) =  0.22000D-03
        HTT(12, 6, 1) =  0.16000D-03
        HTT(12, 7, 1) = -0.14300D-02
        HTT(12, 8, 1) = -0.14900D-02
        HTT(12, 9, 1) = -0.24600D-02
        HTT(12,10, 1) =  0.35600D-02
        HTT(12,11, 1) =  0.27100D-02
        HTT(12,12, 1) =  0.87300D-02
        HTT(13, 1, 1) =  0.53000D-03
        HTT(13, 2, 1) = -0.16100D-02
        HTT(13, 3, 1) =  0.26100D-02
        HTT(13, 4, 1) =  0.47000D-03
        HTT(13, 5, 1) =  0.24000D-02
        HTT(13, 6, 1) =  0.89000D-03
        HTT(13, 7, 1) =  0.43000D-03
        HTT(13, 8, 1) = -0.10800D-02
        HTT(13, 9, 1) = -0.66000D-02
        HTT(13,10, 1) = -0.21000D-03
        HTT(13,11, 1) = -0.44500D-02
        HTT(13,12, 1) =  0.19400D-02
        HTT(13,13, 1) =  0.96100D-02
        HTT(14, 1, 1) =  0.20400D-02
        HTT(14, 2, 1) =  0.52000D-03
        HTT(14, 3, 1) =  0.14800D-02
        HTT(14, 4, 1) = -0.40000D-04
        HTT(14, 5, 1) =  0.97000D-03
        HTT(14, 6, 1) = -0.11900D-02
        HTT(14, 7, 1) =  0.15000D-02
        HTT(14, 8, 1) = -0.66000D-03
        HTT(14, 9, 1) = -0.44200D-02
        HTT(14,10, 1) =  0.12200D-02
        HTT(14,11, 1) = -0.29000D-02
        HTT(14,12, 1) =  0.27500D-02
        HTT(14,13, 1) =  0.37500D-02
        HTT(14,14, 1) =  0.85800D-02
        HTT(15, 1, 1) = -0.18400D-02
        HTT(15, 2, 1) = -0.32100D-02
        HTT(15, 3, 1) =  0.12300D-02
        HTT(15, 4, 1) = -0.15000D-03
        HTT(15, 5, 1) =  0.85000D-03
        HTT(15, 6, 1) =  0.14200D-02
        HTT(15, 7, 1) = -0.20600D-02
        HTT(15, 8, 1) = -0.14900D-02
        HTT(15, 9, 1) = -0.19000D-03
        HTT(15,10, 1) = -0.55100D-02
        HTT(15,11, 1) =  0.11900D-02
        HTT(15,12, 1) = -0.41300D-02
        HTT(15,13, 1) =  0.31600D-02
        HTT(15,14, 1) = -0.19400D-02
        HTT(15,15, 1) =  0.85200D-02
        HTT(16, 1, 1) = -0.33000D-03
        HTT(16, 2, 1) = -0.10800D-02
        HTT(16, 3, 1) =  0.10000D-03
        HTT(16, 4, 1) = -0.66000D-03
        HTT(16, 5, 1) = -0.59000D-03
        HTT(16, 6, 1) = -0.66000D-03
        HTT(16, 7, 1) = -0.99000D-03
        HTT(16, 8, 1) = -0.10700D-02
        HTT(16, 9, 1) =  0.19900D-02
        HTT(16,10, 1) = -0.40700D-02
        HTT(16,11, 1) =  0.27500D-02
        HTT(16,12, 1) = -0.33200D-02
        HTT(16,13, 1) = -0.26900D-02
        HTT(16,14, 1) =  0.28900D-02
        HTT(16,15, 1) =  0.34200D-02
        HTT(16,16, 1) =  0.90000D-02
        HTT(17, 1, 1) =  0.25300D-02
        HTT(17, 2, 1) =  0.10400D-02
        HTT(17, 3, 1) =  0.19000D-03
        HTT(17, 4, 1) = -0.13000D-02
        HTT(17, 5, 1) =  0.43000D-03
        HTT(17, 6, 1) = -0.16400D-02
        HTT(17, 7, 1) =  0.26600D-02
        HTT(17, 8, 1) =  0.59000D-03
        HTT(17, 9, 1) =  0.61000D-03
        HTT(17,10, 1) = -0.18300D-02
        HTT(17,11, 1) =  0.21000D-02
        HTT(17,12, 1) = -0.33000D-03
        HTT(17,13, 1) = -0.66400D-02
        HTT(17,14, 1) = -0.70000D-03
        HTT(17,15, 1) = -0.41800D-02
        HTT(17,16, 1) =  0.17600D-02
        HTT(17,17, 1) =  0.96800D-02
        HTT(18, 1, 1) =  0.74000D-03
        HTT(18, 2, 1) =  0.13600D-02
        HTT(18, 3, 1) = -0.19200D-02
        HTT(18, 4, 1) = -0.13000D-02
        HTT(18, 5, 1) = -0.16100D-02
        HTT(18, 6, 1) = -0.30500D-02
        HTT(18, 7, 1) =  0.91000D-03
        HTT(18, 8, 1) = -0.53000D-03
        HTT(18, 9, 1) =  0.21300D-02
        HTT(18,10, 1) =  0.45000D-03
        HTT(18,11, 1) =  0.15100D-02
        HTT(18,12, 1) = -0.17000D-03
        HTT(18,13, 1) = -0.40600D-02
        HTT(18,14, 1) =  0.10400D-02
        HTT(18,15, 1) = -0.23700D-02
        HTT(18,16, 1) =  0.27300D-02
        HTT(18,17, 1) =  0.31700D-02
        HTT(18,18, 1) =  0.85000D-02
        HTT(19, 1, 1) =  0.10300D-02
        HTT(19, 2, 1) = -0.10900D-02
        HTT(19, 3, 1) =  0.13200D-02
        HTT(19, 4, 1) = -0.79000D-03
        HTT(19, 5, 1) =  0.18700D-02
        HTT(19, 6, 1) =  0.45000D-03
        HTT(19, 7, 1) =  0.15900D-02
        HTT(19, 8, 1) =  0.17000D-03
        HTT(19, 9, 1) = -0.15700D-02
        HTT(19,10, 1) = -0.32700D-02
        HTT(19,11, 1) =  0.55000D-03
        HTT(19,12, 1) = -0.11400D-02
        HTT(19,13, 1) = -0.78000D-03
        HTT(19,14, 1) = -0.55300D-02
        HTT(19,15, 1) =  0.93000D-03
        HTT(19,16, 1) = -0.38200D-02
        HTT(19,17, 1) =  0.37300D-02
        HTT(19,18, 1) = -0.19300D-02
        HTT(19,19, 1) =  0.84800D-02
        HTT(20, 1, 1) = -0.77000D-03
        HTT(20, 2, 1) = -0.77000D-03
        HTT(20, 3, 1) = -0.79000D-03
        HTT(20, 4, 1) = -0.79000D-03
        HTT(20, 5, 1) = -0.18000D-03
        HTT(20, 6, 1) = -0.97000D-03
        HTT(20, 7, 1) = -0.16000D-03
        HTT(20, 8, 1) = -0.95000D-03
        HTT(20, 9, 1) = -0.50000D-04
        HTT(20,10, 1) = -0.98000D-03
        HTT(20,11, 1) = -0.50000D-04
        HTT(20,12, 1) = -0.98000D-03
        HTT(20,13, 1) =  0.18000D-02
        HTT(20,14, 1) = -0.37900D-02
        HTT(20,15, 1) =  0.27400D-02
        HTT(20,16, 1) = -0.28500D-02
        HTT(20,17, 1) = -0.27700D-02
        HTT(20,18, 1) =  0.34000D-02
        HTT(20,19, 1) =  0.28200D-02
        HTT(20,20, 1) =  0.89900D-02
        HTT(21, 1, 1) =  0.30000D-03
        HTT(21, 2, 1) =  0.19600D-02
        HTT(21, 3, 1) = -0.25100D-02
        HTT(21, 4, 1) = -0.85000D-03
        HTT(21, 5, 1) = -0.60900D-02
        HTT(21, 6, 1) = -0.58000D-03
        HTT(21, 7, 1) = -0.34300D-02
        HTT(21, 8, 1) =  0.20800D-02
        HTT(21, 9, 1) =  0.24300D-02
        HTT(21,10, 1) =  0.93000D-03
        HTT(21,11, 1) =  0.76000D-03
        HTT(21,12, 1) = -0.74000D-03
        HTT(21,13, 1) =  0.68000D-03
        HTT(21,14, 1) = -0.16400D-02
        HTT(21,15, 1) =  0.21900D-02
        HTT(21,16, 1) = -0.13000D-03
        HTT(21,17, 1) = -0.66000D-02
        HTT(21,18, 1) = -0.36000D-03
        HTT(21,19, 1) = -0.42700D-02
        HTT(21,20, 1) =  0.19700D-02
        HTT(21,21, 1) =  0.92800D-02
        HTT(22, 1, 1) = -0.16800D-02
        HTT(22, 2, 1) =  0.42000D-03
        HTT(22, 3, 1) = -0.33400D-02
        HTT(22, 4, 1) = -0.12400D-02
        HTT(22, 5, 1) = -0.50000D-03
        HTT(22, 6, 1) = -0.52300D-02
        HTT(22, 7, 1) =  0.10700D-02
        HTT(22, 8, 1) = -0.36600D-02
        HTT(22, 9, 1) =  0.80000D-03
        HTT(22,10, 1) =  0.14300D-02
        HTT(22,11, 1) = -0.13000D-02
        HTT(22,12, 1) = -0.68000D-03
        HTT(22,13, 1) =  0.22300D-02
        HTT(22,14, 1) =  0.58000D-03
        HTT(22,15, 1) =  0.16000D-02
        HTT(22,16, 1) = -0.60000D-04
        HTT(22,17, 1) = -0.44700D-02
        HTT(22,18, 1) =  0.11200D-02
        HTT(22,19, 1) = -0.28100D-02
        HTT(22,20, 1) =  0.27800D-02
        HTT(22,21, 1) =  0.36100D-02
        HTT(22,22, 1) =  0.84800D-02
        HTT(23, 1, 1) =  0.20800D-02
        HTT(23, 2, 1) =  0.16400D-02
        HTT(23, 3, 1) = -0.41000D-03
        HTT(23, 4, 1) = -0.85000D-03
        HTT(23, 5, 1) = -0.40600D-02
        HTT(23, 6, 1) =  0.83000D-03
        HTT(23, 7, 1) = -0.16900D-02
        HTT(23, 8, 1) =  0.31900D-02
        HTT(23, 9, 1) =  0.92000D-03
        HTT(23,10, 1) = -0.13400D-02
        HTT(23,11, 1) =  0.13500D-02
        HTT(23,12, 1) = -0.90000D-03
        HTT(23,13, 1) = -0.18800D-02
        HTT(23,14, 1) = -0.33700D-02
        HTT(23,15, 1) =  0.39000D-03
        HTT(23,16, 1) = -0.11000D-02
        HTT(23,17, 1) = -0.14000D-03
        HTT(23,18, 1) = -0.56500D-02
        HTT(23,19, 1) =  0.13500D-02
        HTT(23,20, 1) = -0.41600D-02
        HTT(23,21, 1) =  0.30900D-02
        HTT(23,22, 1) = -0.19500D-02
        HTT(23,23, 1) =  0.85600D-02
        HTT(24, 1, 1) =  0.10000D-03
        HTT(24, 2, 1) =  0.11000D-03
        HTT(24, 3, 1) = -0.12500D-02
        HTT(24, 4, 1) = -0.12400D-02
        HTT(24, 5, 1) =  0.15300D-02
        HTT(24, 6, 1) = -0.38300D-02
        HTT(24, 7, 1) =  0.28100D-02
        HTT(24, 8, 1) = -0.25500D-02
        HTT(24, 9, 1) = -0.71000D-03
        HTT(24,10, 1) = -0.84000D-03
        HTT(24,11, 1) = -0.71000D-03
        HTT(24,12, 1) = -0.84000D-03
        HTT(24,13, 1) = -0.32000D-03
        HTT(24,14, 1) = -0.11500D-02
        HTT(24,15, 1) = -0.20000D-03
        HTT(24,16, 1) = -0.10200D-02
        HTT(24,17, 1) =  0.19900D-02
        HTT(24,18, 1) = -0.41700D-02
        HTT(24,19, 1) =  0.28200D-02
        HTT(24,20, 1) = -0.33500D-02
        HTT(24,21, 1) = -0.25900D-02
        HTT(24,22, 1) =  0.29200D-02
        HTT(24,23, 1) =  0.35300D-02
        HTT(24,24, 1) =  0.90500D-02
C Hessian elements for anchor point 1
        HBS( 1, 1, 1) = -0.36310D-01
        HBS( 1, 2, 1) = -0.36690D-01
        HBS( 1, 3, 1) = -0.36890D-01
        HBS( 1, 4, 1) = -0.11400D-01
        HBS( 1, 5, 1) =  0.55900D-02
        HBS( 1, 6, 1) =  0.45710D-01
        HBS( 1, 7, 1) = -0.44000D-03
        HBS( 1, 8, 1) =  0.44410D-01
        HBS( 1, 9, 1) =  0.60000D-04
        HBS( 1,10, 1) = -0.10640D-01
        HBS( 1,11, 1) = -0.46000D-03
        HBS( 1,12, 1) =  0.55900D-02
        HBS( 2, 1, 1) =  0.42410D-01
        HBS( 2, 2, 1) = -0.54500D-02
        HBS( 2, 3, 1) =  0.11600D-01
        HBS( 2, 4, 1) = -0.17300D-02
        HBS( 2, 5, 1) = -0.28000D-02
        HBS( 2, 6, 1) = -0.19690D-01
        HBS( 2, 7, 1) =  0.11300D-02
        HBS( 2, 8, 1) = -0.26090D-01
        HBS( 2, 9, 1) = -0.18000D-03
        HBS( 2,10, 1) =  0.11810D-01
        HBS( 2,11, 1) = -0.69000D-03
        HBS( 2,12, 1) = -0.11800D-02
        HBS( 3, 1, 1) = -0.61000D-02
        HBS( 3, 2, 1) =  0.42140D-01
        HBS( 3, 3, 1) =  0.25290D-01
        HBS( 3, 4, 1) =  0.13130D-01
        HBS( 3, 5, 1) = -0.27900D-02
        HBS( 3, 6, 1) = -0.26020D-01
        HBS( 3, 7, 1) = -0.68000D-03
        HBS( 3, 8, 1) = -0.18310D-01
        HBS( 3, 9, 1) =  0.12000D-03
        HBS( 3,10, 1) = -0.11700D-02
        HBS( 3,11, 1) =  0.11600D-02
        HBS( 3,12, 1) = -0.44100D-02
        HBS( 4, 1, 1) = -0.29210D-01
        HBS( 4, 2, 1) = -0.82700D-02
        HBS( 4, 3, 1) =  0.16570D-01
        HBS( 4, 4, 1) = -0.31510D-01
        HBS( 4, 5, 1) = -0.97200D-02
        HBS( 4, 6, 1) = -0.12980D-01
        HBS( 4, 7, 1) =  0.48900D-02
        HBS( 4, 8, 1) =  0.44650D-01
        HBS( 4, 9, 1) = -0.64000D-03
        HBS( 4,10, 1) =  0.43490D-01
        HBS( 4,11, 1) =  0.19000D-03
        HBS( 4,12, 1) = -0.10200D-02
        HBS( 5, 1, 1) =  0.27310D-01
        HBS( 5, 2, 1) =  0.47900D-02
        HBS( 5, 3, 1) = -0.10360D-01
        HBS( 5, 4, 1) =  0.23300D-02
        HBS( 5, 5, 1) =  0.65800D-02
        HBS( 5, 6, 1) =  0.59300D-02
        HBS( 5, 7, 1) = -0.19100D-02
        HBS( 5, 8, 1) = -0.19800D-01
        HBS( 5, 9, 1) = -0.33000D-03
        HBS( 5,10, 1) = -0.24210D-01
        HBS( 5,11, 1) = -0.10000D-03
        HBS( 5,12, 1) =  0.10100D-02
        HBS( 6, 1, 1) =  0.19000D-02
        HBS( 6, 2, 1) =  0.34900D-02
        HBS( 6, 3, 1) = -0.62000D-02
        HBS( 6, 4, 1) =  0.29170D-01
        HBS( 6, 5, 1) =  0.31400D-02
        HBS( 6, 6, 1) =  0.70400D-02
        HBS( 6, 7, 1) = -0.29800D-02
        HBS( 6, 8, 1) = -0.24840D-01
        HBS( 6, 9, 1) =  0.97000D-03
        HBS( 6,10, 1) = -0.19290D-01
        HBS( 6,11, 1) = -0.90000D-04
        HBS( 6,12, 1) =  0.10000D-04
        HBS( 7, 1, 1) = -0.12120D-01
        HBS( 7, 2, 1) =  0.42330D-01
        HBS( 7, 3, 1) =  0.37900D-02
        HBS( 7, 4, 1) = -0.32360D-01
        HBS( 7, 5, 1) =  0.55000D-02
        HBS( 7, 6, 1) = -0.33060D-01
        HBS( 7, 7, 1) = -0.89600D-02
        HBS( 7, 8, 1) = -0.12200D-01
        HBS( 7, 9, 1) =  0.50200D-02
        HBS( 7,10, 1) =  0.46020D-01
        HBS( 7,11, 1) = -0.96000D-03
        HBS( 7,12, 1) = -0.80000D-04
        HBS( 8, 1, 1) =  0.67000D-02
        HBS( 8, 2, 1) = -0.23490D-01
        HBS( 8, 3, 1) = -0.78000D-03
        HBS( 8, 4, 1) =  0.30060D-01
        HBS( 8, 5, 1) = -0.33500D-02
        HBS( 8, 6, 1) =  0.37600D-02
        HBS( 8, 7, 1) =  0.46600D-02
        HBS( 8, 8, 1) =  0.52800D-02
        HBS( 8, 9, 1) = -0.20300D-02
        HBS( 8,10, 1) = -0.20610D-01
        HBS( 8,11, 1) = -0.19000D-03
        HBS( 8,12, 1) =  0.30000D-04
        HBS( 9, 1, 1) =  0.54200D-02
        HBS( 9, 2, 1) = -0.18840D-01
        HBS( 9, 3, 1) = -0.30100D-02
        HBS( 9, 4, 1) =  0.22900D-02
        HBS( 9, 5, 1) = -0.21500D-02
        HBS( 9, 6, 1) =  0.29300D-01
        HBS( 9, 7, 1) =  0.43100D-02
        HBS( 9, 8, 1) =  0.69200D-02
        HBS( 9, 9, 1) = -0.29800D-02
        HBS( 9,10, 1) = -0.25420D-01
        HBS( 9,11, 1) =  0.11500D-02
        HBS( 9,12, 1) =  0.50000D-04
        HBS(10, 1, 1) =  0.44050D-01
        HBS(10, 2, 1) =  0.43680D-01
        HBS(10, 3, 1) = -0.27100D-02
        HBS(10, 4, 1) = -0.13440D-01
        HBS(10, 5, 1) = -0.64000D-03
        HBS(10, 6, 1) = -0.32210D-01
        HBS(10, 7, 1) =  0.53500D-02
        HBS(10, 8, 1) = -0.31660D-01
        HBS(10, 9, 1) = -0.87900D-02
        HBS(10,10, 1) = -0.13050D-01
        HBS(10,11, 1) =  0.52900D-02
        HBS(10,12, 1) = -0.50000D-03
        HBS(11, 1, 1) = -0.24480D-01
        HBS(11, 2, 1) = -0.19450D-01
        HBS(11, 3, 1) =  0.11700D-02
        HBS(11, 4, 1) =  0.71100D-02
        HBS(11, 5, 1) =  0.10200D-02
        HBS(11, 6, 1) =  0.30010D-01
        HBS(11, 7, 1) = -0.31900D-02
        HBS(11, 8, 1) =  0.20900D-02
        HBS(11, 9, 1) =  0.44300D-02
        HBS(11,10, 1) =  0.60900D-02
        HBS(11,11, 1) = -0.21400D-02
        HBS(11,12, 1) = -0.33000D-03
        HBS(12, 1, 1) = -0.19570D-01
        HBS(12, 2, 1) = -0.24220D-01
        HBS(12, 3, 1) =  0.15500D-02
        HBS(12, 4, 1) =  0.63300D-02
        HBS(12, 5, 1) = -0.38000D-03
        HBS(12, 6, 1) =  0.22000D-02
        HBS(12, 7, 1) = -0.21500D-02
        HBS(12, 8, 1) =  0.29570D-01
        HBS(12, 9, 1) =  0.43600D-02
        HBS(12,10, 1) =  0.69600D-02
        HBS(12,11, 1) = -0.31500D-02
        HBS(12,12, 1) =  0.84000D-03
        HBS(13, 1, 1) =  0.42020D-01
        HBS(13, 2, 1) = -0.11720D-01
        HBS(13, 3, 1) =  0.43200D-02
        HBS(13, 4, 1) =  0.45370D-01
        HBS(13, 5, 1) = -0.30000D-04
        HBS(13, 6, 1) = -0.12300D-01
        HBS(13, 7, 1) = -0.10000D-02
        HBS(13, 8, 1) = -0.32760D-01
        HBS(13, 9, 1) =  0.50300D-02
        HBS(13,10, 1) = -0.33540D-01
        HBS(13,11, 1) = -0.89700D-02
        HBS(13,12, 1) =  0.49300D-02
        HBS(14, 1, 1) = -0.18670D-01
        HBS(14, 2, 1) =  0.51600D-02
        HBS(14, 3, 1) = -0.32600D-02
        HBS(14, 4, 1) = -0.25130D-01
        HBS(14, 5, 1) =  0.30000D-04
        HBS(14, 6, 1) =  0.69900D-02
        HBS(14, 7, 1) =  0.11900D-02
        HBS(14, 8, 1) =  0.29240D-01
        HBS(14, 9, 1) = -0.29900D-02
        HBS(14,10, 1) =  0.30100D-02
        HBS(14,11, 1) =  0.43500D-02
        HBS(14,12, 1) = -0.19200D-02
        HBS(15, 1, 1) = -0.23350D-01
        HBS(15, 2, 1) =  0.65500D-02
        HBS(15, 3, 1) = -0.10600D-02
        HBS(15, 4, 1) = -0.20240D-01
        HBS(15, 5, 1) =  0.00000D+00
        HBS(15, 6, 1) =  0.53100D-02
        HBS(15, 7, 1) = -0.18000D-03
        HBS(15, 8, 1) =  0.35200D-02
        HBS(15, 9, 1) = -0.20400D-02
        HBS(15,10, 1) =  0.30530D-01
        HBS(15,11, 1) =  0.46200D-02
        HBS(15,12, 1) = -0.30100D-02
        HBS(16, 1, 1) = -0.84300D-02
        HBS(16, 2, 1) = -0.29330D-01
        HBS(16, 3, 1) =  0.14920D-01
        HBS(16, 4, 1) =  0.43340D-01
        HBS(16, 5, 1) = -0.70000D-03
        HBS(16, 6, 1) =  0.44840D-01
        HBS(16, 7, 1) =  0.17000D-03
        HBS(16, 8, 1) = -0.12440D-01
        HBS(16, 9, 1) = -0.67000D-03
        HBS(16,10, 1) = -0.32290D-01
        HBS(16,11, 1) =  0.49100D-02
        HBS(16,12, 1) = -0.89200D-02
        HBS(17, 1, 1) =  0.44500D-02
        HBS(17, 2, 1) =  0.25980D-01
        HBS(17, 3, 1) = -0.89400D-02
        HBS(17, 4, 1) = -0.24010D-01
        HBS(17, 5, 1) =  0.10600D-02
        HBS(17, 6, 1) = -0.19750D-01
        HBS(17, 7, 1) = -0.90000D-04
        HBS(17, 8, 1) =  0.59800D-02
        HBS(17, 9, 1) = -0.29000D-03
        HBS(17,10, 1) =  0.31400D-02
        HBS(17,11, 1) = -0.18900D-02
        HBS(17,12, 1) =  0.71100D-02
        HBS(18, 1, 1) =  0.39800D-02
        HBS(18, 2, 1) =  0.33600D-02
        HBS(18, 3, 1) = -0.59800D-02
        HBS(18, 4, 1) = -0.19330D-01
        HBS(18, 5, 1) = -0.36000D-03
        HBS(18, 6, 1) = -0.25090D-01
        HBS(18, 7, 1) = -0.80000D-04
        HBS(18, 8, 1) =  0.64600D-02
        HBS(18, 9, 1) =  0.96000D-03
        HBS(18,10, 1) =  0.29140D-01
        HBS(18,11, 1) = -0.30200D-02
        HBS(18,12, 1) =  0.18100D-02
C Hessian elements for anchor point 2
        HSS( 1, 1, 2) =  0.40931D+00
        HSS( 2, 1, 2) =  0.17750D-01
        HSS( 2, 2, 2) =  0.39831D+00
        HSS( 3, 1, 2) =  0.30090D-01
        HSS( 3, 2, 2) =  0.39670D-01
        HSS( 3, 3, 2) =  0.41039D+00
        HSS( 4, 1, 2) =  0.17380D-01
        HSS( 4, 2, 2) = -0.24100D-02
        HSS( 4, 3, 2) = -0.47000D-03
        HSS( 4, 4, 2) =  0.41129D+00
        HSS( 5, 1, 2) =  0.42800D-02
        HSS( 5, 2, 2) = -0.46000D-03
        HSS( 5, 3, 2) =  0.49400D-02
        HSS( 5, 4, 2) =  0.55000D-02
        HSS( 5, 5, 2) =  0.35488D+00
        HSS( 6, 1, 2) = -0.36200D-02
        HSS( 6, 2, 2) =  0.81920D-01
        HSS( 6, 3, 2) = -0.25500D-02
        HSS( 6, 4, 2) =  0.21100D-01
        HSS( 6, 5, 2) = -0.76000D-03
        HSS( 6, 6, 2) =  0.41554D+00
        HSS( 7, 1, 2) = -0.59000D-03
        HSS( 7, 2, 2) = -0.36000D-03
        HSS( 7, 3, 2) =  0.73000D-03
        HSS( 7, 4, 2) =  0.43500D-02
        HSS( 7, 5, 2) =  0.47000D-03
        HSS( 7, 6, 2) =  0.46000D-02
        HSS( 7, 7, 2) =  0.36429D+00
        HSS( 8, 1, 2) =  0.81370D-01
        HSS( 8, 2, 2) = -0.43100D-02
        HSS( 8, 3, 2) =  0.46000D-03
        HSS( 8, 4, 2) = -0.31800D-02
        HSS( 8, 5, 2) = -0.28000D-03
        HSS( 8, 6, 2) =  0.21080D-01
        HSS( 8, 7, 2) = -0.32000D-03
        HSS( 8, 8, 2) =  0.40703D+00
        HSS( 9, 1, 2) = -0.47000D-03
        HSS( 9, 2, 2) = -0.50000D-03
        HSS( 9, 3, 2) = -0.34000D-03
        HSS( 9, 4, 2) = -0.61000D-03
        HSS( 9, 5, 2) =  0.23000D-03
        HSS( 9, 6, 2) =  0.44800D-02
        HSS( 9, 7, 2) =  0.47000D-03
        HSS( 9, 8, 2) =  0.44100D-02
        HSS( 9, 9, 2) =  0.36691D+00
        HSS(10, 1, 2) = -0.11900D-02
        HSS(10, 2, 2) =  0.18110D-01
        HSS(10, 3, 2) = -0.13500D-02
        HSS(10, 4, 2) =  0.85990D-01
        HSS(10, 5, 2) = -0.79000D-03
        HSS(10, 6, 2) = -0.28600D-02
        HSS(10, 7, 2) = -0.55000D-03
        HSS(10, 8, 2) =  0.21070D-01
        HSS(10, 9, 2) = -0.50000D-03
        HSS(10,10, 2) =  0.41886D+00
        HSS(11, 1, 2) = -0.47000D-03
        HSS(11, 2, 2) = -0.52000D-03
        HSS(11, 3, 2) =  0.68000D-03
        HSS(11, 4, 2) = -0.53000D-03
        HSS(11, 5, 2) =  0.12000D-03
        HSS(11, 6, 2) = -0.24000D-03
        HSS(11, 7, 2) =  0.21000D-03
        HSS(11, 8, 2) =  0.45400D-02
        HSS(11, 9, 2) =  0.47000D-03
        HSS(11,10, 2) =  0.45200D-02
        HSS(11,11, 2) =  0.36397D+00
        HSS(12, 1, 2) = -0.53000D-03
        HSS(12, 2, 2) =  0.41000D-02
        HSS(12, 3, 2) = -0.53000D-03
        HSS(12, 4, 2) = -0.83000D-03
        HSS(12, 5, 2) =  0.34000D-03
        HSS(12, 6, 2) = -0.14000D-03
        HSS(12, 7, 2) =  0.10000D-03
        HSS(12, 8, 2) = -0.69000D-03
        HSS(12, 9, 2) =  0.14000D-03
        HSS(12,10, 2) =  0.42400D-02
        HSS(12,11, 2) =  0.46000D-03
        HSS(12,12, 2) =  0.36700D+00
C Hessian elements for anchor point 2
        HBB( 1, 1, 2) =  0.85740D-01
        HBB( 2, 1, 2) = -0.47120D-01
        HBB( 2, 2, 2) =  0.14272D+00
        HBB( 3, 1, 2) = -0.38620D-01
        HBB( 3, 2, 2) = -0.95590D-01
        HBB( 3, 3, 2) =  0.13421D+00
        HBB( 4, 1, 2) = -0.35340D-01
        HBB( 4, 2, 2) =  0.24370D-01
        HBB( 4, 3, 2) =  0.10970D-01
        HBB( 4, 4, 2) =  0.78390D-01
        HBB( 5, 1, 2) =  0.19660D-01
        HBB( 5, 2, 2) = -0.13050D-01
        HBB( 5, 3, 2) = -0.66100D-02
        HBB( 5, 4, 2) = -0.38100D-01
        HBB( 5, 5, 2) =  0.76190D-01
        HBB( 6, 1, 2) =  0.15680D-01
        HBB( 6, 2, 2) = -0.11320D-01
        HBB( 6, 3, 2) = -0.43700D-02
        HBB( 6, 4, 2) = -0.40290D-01
        HBB( 6, 5, 2) = -0.38090D-01
        HBB( 6, 6, 2) =  0.78390D-01
        HBB( 7, 1, 2) = -0.10880D-01
        HBB( 7, 2, 2) = -0.19300D-02
        HBB( 7, 3, 2) =  0.12810D-01
        HBB( 7, 4, 2) = -0.30710D-01
        HBB( 7, 5, 2) =  0.11850D-01
        HBB( 7, 6, 2) =  0.18860D-01
        HBB( 7, 7, 2) =  0.78510D-01
        HBB( 8, 1, 2) = -0.26000D-03
        HBB( 8, 2, 2) =  0.56200D-02
        HBB( 8, 3, 2) = -0.53700D-02
        HBB( 8, 4, 2) =  0.18260D-01
        HBB( 8, 5, 2) = -0.64700D-02
        HBB( 8, 6, 2) = -0.11790D-01
        HBB( 8, 7, 2) = -0.39180D-01
        HBB( 8, 8, 2) =  0.78940D-01
        HBB( 9, 1, 2) =  0.11130D-01
        HBB( 9, 2, 2) = -0.36900D-02
        HBB( 9, 3, 2) = -0.74400D-02
        HBB( 9, 4, 2) =  0.12450D-01
        HBB( 9, 5, 2) = -0.53800D-02
        HBB( 9, 6, 2) = -0.70800D-02
        HBB( 9, 7, 2) = -0.39340D-01
        HBB( 9, 8, 2) = -0.39760D-01
        HBB( 9, 9, 2) =  0.79100D-01
        HBB(10, 1, 2) =  0.73300D-02
        HBB(10, 2, 2) = -0.43700D-02
        HBB(10, 3, 2) = -0.29600D-02
        HBB(10, 4, 2) = -0.10030D-01
        HBB(10, 5, 2) =  0.10560D-01
        HBB(10, 6, 2) = -0.53000D-03
        HBB(10, 7, 2) = -0.33760D-01
        HBB(10, 8, 2) =  0.13640D-01
        HBB(10, 9, 2) =  0.20120D-01
        HBB(10,10, 2) =  0.80240D-01
        HBB(11, 1, 2) = -0.36200D-02
        HBB(11, 2, 2) =  0.25500D-02
        HBB(11, 3, 2) =  0.10700D-02
        HBB(11, 4, 2) = -0.61000D-03
        HBB(11, 5, 2) = -0.38000D-02
        HBB(11, 6, 2) =  0.44100D-02
        HBB(11, 7, 2) =  0.19800D-01
        HBB(11, 8, 2) = -0.72600D-02
        HBB(11, 9, 2) = -0.12540D-01
        HBB(11,10, 2) = -0.40120D-01
        HBB(11,11, 2) =  0.79050D-01
        HBB(12, 1, 2) = -0.37100D-02
        HBB(12, 2, 2) =  0.18200D-02
        HBB(12, 3, 2) =  0.18900D-02
        HBB(12, 4, 2) =  0.10640D-01
        HBB(12, 5, 2) = -0.67600D-02
        HBB(12, 6, 2) = -0.38800D-02
        HBB(12, 7, 2) =  0.13960D-01
        HBB(12, 8, 2) = -0.63800D-02
        HBB(12, 9, 2) = -0.75800D-02
        HBB(12,10, 2) = -0.40110D-01
        HBB(12,11, 2) = -0.38920D-01
        HBB(12,12, 2) =  0.79040D-01
        HBB(13, 1, 2) = -0.10540D-01
        HBB(13, 2, 2) =  0.13360D-01
        HBB(13, 3, 2) = -0.28100D-02
        HBB(13, 4, 2) =  0.52100D-02
        HBB(13, 5, 2) = -0.22200D-02
        HBB(13, 6, 2) = -0.29900D-02
        HBB(13, 7, 2) = -0.87200D-02
        HBB(13, 8, 2) =  0.10150D-01
        HBB(13, 9, 2) = -0.14300D-02
        HBB(13,10, 2) = -0.34040D-01
        HBB(13,11, 2) =  0.14040D-01
        HBB(13,12, 2) =  0.19990D-01
        HBB(13,13, 2) =  0.79490D-01
        HBB(14, 1, 2) =  0.10840D-01
        HBB(14, 2, 2) = -0.78900D-02
        HBB(14, 3, 2) = -0.29500D-02
        HBB(14, 4, 2) = -0.26600D-02
        HBB(14, 5, 2) =  0.10300D-02
        HBB(14, 6, 2) =  0.16200D-02
        HBB(14, 7, 2) = -0.14800D-02
        HBB(14, 8, 2) = -0.35100D-02
        HBB(14, 9, 2) =  0.49900D-02
        HBB(14,10, 2) =  0.20150D-01
        HBB(14,11, 2) = -0.75900D-02
        HBB(14,12, 2) = -0.12560D-01
        HBB(14,13, 2) = -0.39740D-01
        HBB(14,14, 2) =  0.79580D-01
        HBB(15, 1, 2) = -0.30000D-03
        HBB(15, 2, 2) = -0.54600D-02
        HBB(15, 3, 2) =  0.57600D-02
        HBB(15, 4, 2) = -0.25500D-02
        HBB(15, 5, 2) =  0.11900D-02
        HBB(15, 6, 2) =  0.13700D-02
        HBB(15, 7, 2) =  0.10200D-01
        HBB(15, 8, 2) = -0.66400D-02
        HBB(15, 9, 2) = -0.35600D-02
        HBB(15,10, 2) =  0.13890D-01
        HBB(15,11, 2) = -0.64600D-02
        HBB(15,12, 2) = -0.74300D-02
        HBB(15,13, 2) = -0.39750D-01
        HBB(15,14, 2) = -0.39830D-01
        HBB(15,15, 2) =  0.79580D-01
        HBB(16, 1, 2) = -0.36310D-01
        HBB(16, 2, 2) =  0.15710D-01
        HBB(16, 3, 2) =  0.20610D-01
        HBB(16, 4, 2) = -0.75200D-02
        HBB(16, 5, 2) = -0.17400D-02
        HBB(16, 6, 2) =  0.92600D-02
        HBB(16, 7, 2) =  0.55500D-02
        HBB(16, 8, 2) = -0.26200D-02
        HBB(16, 9, 2) = -0.29400D-02
        HBB(16,10, 2) = -0.97400D-02
        HBB(16,11, 2) =  0.10510D-01
        HBB(16,12, 2) = -0.77000D-03
        HBB(16,13, 2) = -0.31400D-01
        HBB(16,14, 2) =  0.12890D-01
        HBB(16,15, 2) =  0.18510D-01
        HBB(16,16, 2) =  0.79420D-01
        HBB(17, 1, 2) =  0.22000D-01
        HBB(17, 2, 2) = -0.74700D-02
        HBB(17, 3, 2) = -0.14540D-01
        HBB(17, 4, 2) = -0.24600D-02
        HBB(17, 5, 2) =  0.53500D-02
        HBB(17, 6, 2) = -0.28900D-02
        HBB(17, 7, 2) = -0.26800D-02
        HBB(17, 8, 2) =  0.14000D-02
        HBB(17, 9, 2) =  0.12800D-02
        HBB(17,10, 2) =  0.10490D-01
        HBB(17,11, 2) = -0.68100D-02
        HBB(17,12, 2) = -0.36800D-02
        HBB(17,13, 2) =  0.12230D-01
        HBB(17,14, 2) = -0.56800D-02
        HBB(17,15, 2) = -0.65500D-02
        HBB(17,16, 2) = -0.39570D-01
        HBB(17,17, 2) =  0.76420D-01
        HBB(18, 1, 2) =  0.14310D-01
        HBB(18, 2, 2) = -0.82400D-02
        HBB(18, 3, 2) = -0.60700D-02
        HBB(18, 4, 2) =  0.99800D-02
        HBB(18, 5, 2) = -0.36100D-02
        HBB(18, 6, 2) = -0.63700D-02
        HBB(18, 7, 2) = -0.28700D-02
        HBB(18, 8, 2) =  0.12200D-02
        HBB(18, 9, 2) =  0.16500D-02
        HBB(18,10, 2) = -0.75000D-03
        HBB(18,11, 2) = -0.37000D-02
        HBB(18,12, 2) =  0.44500D-02
        HBB(18,13, 2) =  0.19170D-01
        HBB(18,14, 2) = -0.72100D-02
        HBB(18,15, 2) = -0.11960D-01
        HBB(18,16, 2) = -0.39850D-01
        HBB(18,17, 2) = -0.36850D-01
        HBB(18,18, 2) =  0.76700D-01
C Hessian elements for anchor point 2
        HTT( 1, 1, 2) =  0.90500D-02
        HTT( 2, 1, 2) =  0.38100D-02
        HTT( 2, 2, 2) =  0.81400D-02
        HTT( 3, 1, 2) =  0.25000D-03
        HTT( 3, 2, 2) = -0.43200D-02
        HTT( 3, 3, 2) =  0.11860D-01
        HTT( 4, 1, 2) = -0.49900D-02
        HTT( 4, 2, 2) =  0.10000D-04
        HTT( 4, 3, 2) =  0.72800D-02
        HTT( 4, 4, 2) =  0.12280D-01
        HTT( 5, 1, 2) = -0.60200D-02
        HTT( 5, 2, 2) = -0.42900D-02
        HTT( 5, 3, 2) =  0.24900D-02
        HTT( 5, 4, 2) =  0.42200D-02
        HTT( 5, 5, 2) =  0.87300D-02
        HTT( 6, 1, 2) = -0.42000D-02
        HTT( 6, 2, 2) = -0.28900D-02
        HTT( 6, 3, 2) =  0.31600D-02
        HTT( 6, 4, 2) =  0.44800D-02
        HTT( 6, 5, 2) =  0.34600D-02
        HTT( 6, 6, 2) =  0.79900D-02
        HTT( 7, 1, 2) =  0.21800D-02
        HTT( 7, 2, 2) =  0.32900D-02
        HTT( 7, 3, 2) = -0.83100D-02
        HTT( 7, 4, 2) = -0.72000D-02
        HTT( 7, 5, 2) =  0.81000D-03
        HTT( 7, 6, 2) = -0.34000D-02
        HTT( 7, 7, 2) =  0.10570D-01
        HTT( 8, 1, 2) =  0.39900D-02
        HTT( 8, 2, 2) =  0.46900D-02
        HTT( 8, 3, 2) = -0.76400D-02
        HTT( 8, 4, 2) = -0.69400D-02
        HTT( 8, 5, 2) = -0.44700D-02
        HTT( 8, 6, 2) =  0.11300D-02
        HTT( 8, 7, 2) =  0.63600D-02
        HTT( 8, 8, 2) =  0.11960D-01
        HTT( 9, 1, 2) = -0.62900D-02
        HTT( 9, 2, 2) = -0.85000D-03
        HTT( 9, 3, 2) = -0.28600D-02
        HTT( 9, 4, 2) =  0.25800D-02
        HTT( 9, 5, 2) =  0.40000D-03
        HTT( 9, 6, 2) =  0.20200D-02
        HTT( 9, 7, 2) = -0.27900D-02
        HTT( 9, 8, 2) = -0.11700D-02
        HTT( 9, 9, 2) =  0.94300D-02
        HTT(10, 1, 2) = -0.39800D-02
        HTT(10, 2, 2) =  0.72000D-03
        HTT(10, 3, 2) = -0.15100D-02
        HTT(10, 4, 2) =  0.31900D-02
        HTT(10, 5, 2) =  0.20000D-02
        HTT(10, 6, 2) =  0.15200D-02
        HTT(10, 7, 2) = -0.30000D-03
        HTT(10, 8, 2) = -0.78000D-03
        HTT(10, 9, 2) =  0.30800D-02
        HTT(10,10, 2) =  0.83100D-02
        HTT(11, 1, 2) = -0.99000D-03
        HTT(11, 2, 2) = -0.52300D-02
        HTT(11, 3, 2) =  0.17700D-02
        HTT(11, 4, 2) = -0.24700D-02
        HTT(11, 5, 2) = -0.13500D-02
        HTT(11, 6, 2) =  0.69000D-03
        HTT(11, 7, 2) = -0.39100D-02
        HTT(11, 8, 2) = -0.18800D-02
        HTT(11, 9, 2) =  0.39300D-02
        HTT(11,10, 2) = -0.16800D-02
        HTT(11,11, 2) =  0.82200D-02
        HTT(12, 1, 2) =  0.13200D-02
        HTT(12, 2, 2) = -0.36600D-02
        HTT(12, 3, 2) =  0.31200D-02
        HTT(12, 4, 2) = -0.18600D-02
        HTT(12, 5, 2) =  0.24000D-03
        HTT(12, 6, 2) =  0.19000D-03
        HTT(12, 7, 2) = -0.14300D-02
        HTT(12, 8, 2) = -0.14900D-02
        HTT(12, 9, 2) = -0.24200D-02
        HTT(12,10, 2) =  0.35600D-02
        HTT(12,11, 2) =  0.26100D-02
        HTT(12,12, 2) =  0.85900D-02
        HTT(13, 1, 2) =  0.48000D-03
        HTT(13, 2, 2) = -0.16000D-02
        HTT(13, 3, 2) =  0.26700D-02
        HTT(13, 4, 2) =  0.58000D-03
        HTT(13, 5, 2) =  0.24300D-02
        HTT(13, 6, 2) =  0.89000D-03
        HTT(13, 7, 2) =  0.40000D-03
        HTT(13, 8, 2) = -0.11500D-02
        HTT(13, 9, 2) = -0.65700D-02
        HTT(13,10, 2) = -0.18000D-03
        HTT(13,11, 2) = -0.44600D-02
        HTT(13,12, 2) =  0.19300D-02
        HTT(13,13, 2) =  0.95900D-02
        HTT(14, 1, 2) =  0.19800D-02
        HTT(14, 2, 2) =  0.50000D-03
        HTT(14, 3, 2) =  0.14700D-02
        HTT(14, 4, 2) = -0.10000D-04
        HTT(14, 5, 2) =  0.10300D-02
        HTT(14, 6, 2) = -0.11500D-02
        HTT(14, 7, 2) =  0.15100D-02
        HTT(14, 8, 2) = -0.68000D-03
        HTT(14, 9, 2) = -0.44200D-02
        HTT(14,10, 2) =  0.12200D-02
        HTT(14,11, 2) = -0.29300D-02
        HTT(14,12, 2) =  0.27100D-02
        HTT(14,13, 2) =  0.37900D-02
        HTT(14,14, 2) =  0.85500D-02
        HTT(15, 1, 2) = -0.18400D-02
        HTT(15, 2, 2) = -0.31900D-02
        HTT(15, 3, 2) =  0.13100D-02
        HTT(15, 4, 2) = -0.40000D-04
        HTT(15, 5, 2) =  0.83000D-03
        HTT(15, 6, 2) =  0.14000D-02
        HTT(15, 7, 2) = -0.21100D-02
        HTT(15, 8, 2) = -0.15400D-02
        HTT(15, 9, 2) = -0.16000D-03
        HTT(15,10, 2) = -0.54600D-02
        HTT(15,11, 2) =  0.12000D-02
        HTT(15,12, 2) = -0.41000D-02
        HTT(15,13, 2) =  0.31500D-02
        HTT(15,14, 2) = -0.19100D-02
        HTT(15,15, 2) =  0.85000D-02
        HTT(16, 1, 2) = -0.35000D-03
        HTT(16, 2, 2) = -0.10800D-02
        HTT(16, 3, 2) =  0.11000D-03
        HTT(16, 4, 2) = -0.62000D-03
        HTT(16, 5, 2) = -0.57000D-03
        HTT(16, 6, 2) = -0.64000D-03
        HTT(16, 7, 2) = -0.10000D-02
        HTT(16, 8, 2) = -0.10700D-02
        HTT(16, 9, 2) =  0.19800D-02
        HTT(16,10, 2) = -0.40600D-02
        HTT(16,11, 2) =  0.27200D-02
        HTT(16,12, 2) = -0.33200D-02
        HTT(16,13, 2) = -0.26500D-02
        HTT(16,14, 2) =  0.28600D-02
        HTT(16,15, 2) =  0.34400D-02
        HTT(16,16, 2) =  0.89500D-02
        HTT(17, 1, 2) =  0.25600D-02
        HTT(17, 2, 2) =  0.11000D-02
        HTT(17, 3, 2) =  0.12000D-03
        HTT(17, 4, 2) = -0.13400D-02
        HTT(17, 5, 2) =  0.36000D-03
        HTT(17, 6, 2) = -0.16300D-02
        HTT(17, 7, 2) =  0.26400D-02
        HTT(17, 8, 2) =  0.64000D-03
        HTT(17, 9, 2) =  0.63000D-03
        HTT(17,10, 2) = -0.18200D-02
        HTT(17,11, 2) =  0.21000D-02
        HTT(17,12, 2) = -0.35000D-03
        HTT(17,13, 2) = -0.66700D-02
        HTT(17,14, 2) = -0.74000D-03
        HTT(17,15, 2) = -0.42000D-02
        HTT(17,16, 2) =  0.17300D-02
        HTT(17,17, 2) =  0.96700D-02
        HTT(18, 1, 2) =  0.74000D-03
        HTT(18, 2, 2) =  0.13700D-02
        HTT(18, 3, 2) = -0.19900D-02
        HTT(18, 4, 2) = -0.13600D-02
        HTT(18, 5, 2) = -0.16000D-02
        HTT(18, 6, 2) = -0.30000D-02
        HTT(18, 7, 2) =  0.94000D-03
        HTT(18, 8, 2) = -0.46000D-03
        HTT(18, 9, 2) =  0.21100D-02
        HTT(18,10, 2) =  0.44000D-03
        HTT(18,11, 2) =  0.14700D-02
        HTT(18,12, 2) = -0.19000D-03
        HTT(18,13, 2) = -0.40300D-02
        HTT(18,14, 2) =  0.10300D-02
        HTT(18,15, 2) = -0.23600D-02
        HTT(18,16, 2) =  0.27100D-02
        HTT(18,17, 2) =  0.31700D-02
        HTT(18,18, 2) =  0.84700D-02
        HTT(19, 1, 2) =  0.10700D-02
        HTT(19, 2, 2) = -0.10000D-02
        HTT(19, 3, 2) =  0.13200D-02
        HTT(19, 4, 2) = -0.75000D-03
        HTT(19, 5, 2) =  0.17600D-02
        HTT(19, 6, 2) =  0.41000D-03
        HTT(19, 7, 2) =  0.15300D-02
        HTT(19, 8, 2) =  0.17000D-03
        HTT(19, 9, 2) = -0.15200D-02
        HTT(19,10, 2) = -0.32300D-02
        HTT(19,11, 2) =  0.58000D-03
        HTT(19,12, 2) = -0.11300D-02
        HTT(19,13, 2) = -0.86000D-03
        HTT(19,14, 2) = -0.55100D-02
        HTT(19,15, 2) =  0.86000D-03
        HTT(19,16, 2) = -0.37900D-02
        HTT(19,17, 2) =  0.37400D-02
        HTT(19,18, 2) = -0.18900D-02
        HTT(19,19, 2) =  0.84000D-02
        HTT(20, 1, 2) = -0.75000D-03
        HTT(20, 2, 2) = -0.74000D-03
        HTT(20, 3, 2) = -0.79000D-03
        HTT(20, 4, 2) = -0.78000D-03
        HTT(20, 5, 2) = -0.20000D-03
        HTT(20, 6, 2) = -0.96000D-03
        HTT(20, 7, 2) = -0.17000D-03
        HTT(20, 8, 2) = -0.93000D-03
        HTT(20, 9, 2) = -0.40000D-04
        HTT(20,10, 2) = -0.97000D-03
        HTT(20,11, 2) = -0.50000D-04
        HTT(20,12, 2) = -0.98000D-03
        HTT(20,13, 2) =  0.17700D-02
        HTT(20,14, 2) = -0.37400D-02
        HTT(20,15, 2) =  0.27000D-02
        HTT(20,16, 2) = -0.28100D-02
        HTT(20,17, 2) = -0.27500D-02
        HTT(20,18, 2) =  0.34100D-02
        HTT(20,19, 2) =  0.27600D-02
        HTT(20,20, 2) =  0.89200D-02
        HTT(21, 1, 2) =  0.21000D-03
        HTT(21, 2, 2) =  0.18300D-02
        HTT(21, 3, 2) = -0.26900D-02
        HTT(21, 4, 2) = -0.10600D-02
        HTT(21, 5, 2) = -0.59300D-02
        HTT(21, 6, 2) = -0.53000D-03
        HTT(21, 7, 2) = -0.32400D-02
        HTT(21, 8, 2) =  0.21600D-02
        HTT(21, 9, 2) =  0.24300D-02
        HTT(21,10, 2) =  0.91000D-03
        HTT(21,11, 2) =  0.79000D-03
        HTT(21,12, 2) = -0.73000D-03
        HTT(21,13, 2) =  0.69000D-03
        HTT(21,14, 2) = -0.16600D-02
        HTT(21,15, 2) =  0.22200D-02
        HTT(21,16, 2) = -0.12000D-03
        HTT(21,17, 2) = -0.65400D-02
        HTT(21,18, 2) = -0.37000D-03
        HTT(21,19, 2) = -0.42000D-02
        HTT(21,20, 2) =  0.19700D-02
        HTT(21,21, 2) =  0.91600D-02
        HTT(22, 1, 2) = -0.16600D-02
        HTT(22, 2, 2) =  0.40000D-03
        HTT(22, 3, 2) = -0.33700D-02
        HTT(22, 4, 2) = -0.13200D-02
        HTT(22, 5, 2) = -0.52000D-03
        HTT(22, 6, 2) = -0.51900D-02
        HTT(22, 7, 2) =  0.10800D-02
        HTT(22, 8, 2) = -0.35900D-02
        HTT(22, 9, 2) =  0.77000D-03
        HTT(22,10, 2) =  0.14000D-02
        HTT(22,11, 2) = -0.13100D-02
        HTT(22,12, 2) = -0.68000D-03
        HTT(22,13, 2) =  0.22700D-02
        HTT(22,14, 2) =  0.59000D-03
        HTT(22,15, 2) =  0.16300D-02
        HTT(22,16, 2) = -0.50000D-04
        HTT(22,17, 2) = -0.44900D-02
        HTT(22,18, 2) =  0.10700D-02
        HTT(22,19, 2) = -0.28100D-02
        HTT(22,20, 2) =  0.27500D-02
        HTT(22,21, 2) =  0.36200D-02
        HTT(22,22, 2) =  0.84200D-02
        HTT(23, 1, 2) =  0.20100D-02
        HTT(23, 2, 2) =  0.15700D-02
        HTT(23, 3, 2) = -0.60000D-03
        HTT(23, 4, 2) = -0.10400D-02
        HTT(23, 5, 2) = -0.39900D-02
        HTT(23, 6, 2) =  0.82000D-03
        HTT(23, 7, 2) = -0.15600D-02
        HTT(23, 8, 2) =  0.32500D-02
        HTT(23, 9, 2) =  0.97000D-03
        HTT(23,10, 2) = -0.13300D-02
        HTT(23,11, 2) =  0.14100D-02
        HTT(23,12, 2) = -0.89000D-03
        HTT(23,13, 2) = -0.19300D-02
        HTT(23,14, 2) = -0.34100D-02
        HTT(23,15, 2) =  0.39000D-03
        HTT(23,16, 2) = -0.11000D-02
        HTT(23,17, 2) = -0.90000D-04
        HTT(23,18, 2) = -0.56300D-02
        HTT(23,19, 2) =  0.13900D-02
        HTT(23,20, 2) = -0.41400D-02
        HTT(23,21, 2) =  0.30400D-02
        HTT(23,22, 2) = -0.18900D-02
        HTT(23,23, 2) =  0.85300D-02
        HTT(24, 1, 2) =  0.15000D-03
        HTT(24, 2, 2) =  0.14000D-03
        HTT(24, 3, 2) = -0.12900D-02
        HTT(24, 4, 2) = -0.13000D-02
        HTT(24, 5, 2) =  0.14300D-02
        HTT(24, 6, 2) = -0.38400D-02
        HTT(24, 7, 2) =  0.27700D-02
        HTT(24, 8, 2) = -0.25000D-02
        HTT(24, 9, 2) = -0.70000D-03
        HTT(24,10, 2) = -0.84000D-03
        HTT(24,11, 2) = -0.68000D-03
        HTT(24,12, 2) = -0.83000D-03
        HTT(24,13, 2) = -0.34000D-03
        HTT(24,14, 2) = -0.11700D-02
        HTT(24,15, 2) = -0.19000D-03
        HTT(24,16, 2) = -0.10200D-02
        HTT(24,17, 2) =  0.19500D-02
        HTT(24,18, 2) = -0.41900D-02
        HTT(24,19, 2) =  0.27800D-02
        HTT(24,20, 2) = -0.33600D-02
        HTT(24,21, 2) = -0.25000D-02
        HTT(24,22, 2) =  0.29000D-02
        HTT(24,23, 2) =  0.36000D-02
        HTT(24,24, 2) =  0.90000D-02
C Hessian elements for anchor point 2
        HBS( 1, 1, 2) = -0.37150D-01
        HBS( 1, 2, 2) = -0.35080D-01
        HBS( 1, 3, 2) = -0.37820D-01
        HBS( 1, 4, 2) = -0.11490D-01
        HBS( 1, 5, 2) =  0.62300D-02
        HBS( 1, 6, 2) =  0.46070D-01
        HBS( 1, 7, 2) = -0.43000D-03
        HBS( 1, 8, 2) =  0.44740D-01
        HBS( 1, 9, 2) =  0.80000D-04
        HBS( 1,10, 2) = -0.10140D-01
        HBS( 1,11, 2) = -0.44000D-03
        HBS( 1,12, 2) =  0.53800D-02
        HBS( 2, 1, 2) =  0.45070D-01
        HBS( 2, 2, 2) = -0.50600D-02
        HBS( 2, 3, 2) =  0.16200D-01
        HBS( 2, 4, 2) = -0.18400D-02
        HBS( 2, 5, 2) = -0.27200D-02
        HBS( 2, 6, 2) = -0.20110D-01
        HBS( 2, 7, 2) =  0.11400D-02
        HBS( 2, 8, 2) = -0.26860D-01
        HBS( 2, 9, 2) = -0.17000D-03
        HBS( 2,10, 2) =  0.11580D-01
        HBS( 2,11, 2) = -0.64000D-03
        HBS( 2,12, 2) = -0.16600D-02
        HBS( 3, 1, 2) = -0.79200D-02
        HBS( 3, 2, 2) =  0.40140D-01
        HBS( 3, 3, 2) =  0.21620D-01
        HBS( 3, 4, 2) =  0.13330D-01
        HBS( 3, 5, 2) = -0.35100D-02
        HBS( 3, 6, 2) = -0.25960D-01
        HBS( 3, 7, 2) = -0.71000D-03
        HBS( 3, 8, 2) = -0.17880D-01
        HBS( 3, 9, 2) =  0.10000D-03
        HBS( 3,10, 2) = -0.14400D-02
        HBS( 3,11, 2) =  0.10800D-02
        HBS( 3,12, 2) = -0.37200D-02
        HBS( 4, 1, 2) = -0.29410D-01
        HBS( 4, 2, 2) = -0.87500D-02
        HBS( 4, 3, 2) =  0.18410D-01
        HBS( 4, 4, 2) = -0.31300D-01
        HBS( 4, 5, 2) = -0.10370D-01
        HBS( 4, 6, 2) = -0.13280D-01
        HBS( 4, 7, 2) =  0.49300D-02
        HBS( 4, 8, 2) =  0.44510D-01
        HBS( 4, 9, 2) = -0.65000D-03
        HBS( 4,10, 2) =  0.43310D-01
        HBS( 4,11, 2) =  0.18000D-03
        HBS( 4,12, 2) = -0.83000D-03
        HBS( 5, 1, 2) =  0.28150D-01
        HBS( 5, 2, 2) =  0.54900D-02
        HBS( 5, 3, 2) = -0.13780D-01
        HBS( 5, 4, 2) =  0.21700D-02
        HBS( 5, 5, 2) =  0.70200D-02
        HBS( 5, 6, 2) =  0.60300D-02
        HBS( 5, 7, 2) = -0.19500D-02
        HBS( 5, 8, 2) = -0.19880D-01
        HBS( 5, 9, 2) = -0.33000D-03
        HBS( 5,10, 2) = -0.24160D-01
        HBS( 5,11, 2) = -0.12000D-03
        HBS( 5,12, 2) =  0.77000D-03
        HBS( 6, 1, 2) =  0.12600D-02
        HBS( 6, 2, 2) =  0.32500D-02
        HBS( 6, 3, 2) = -0.46300D-02
        HBS( 6, 4, 2) =  0.29130D-01
        HBS( 6, 5, 2) =  0.33500D-02
        HBS( 6, 6, 2) =  0.72400D-02
        HBS( 6, 7, 2) = -0.29800D-02
        HBS( 6, 8, 2) = -0.24630D-01
        HBS( 6, 9, 2) =  0.99000D-03
        HBS( 6,10, 2) = -0.19140D-01
        HBS( 6,11, 2) = -0.70000D-04
        HBS( 6,12, 2) =  0.60000D-04
        HBS( 7, 1, 2) = -0.11820D-01
        HBS( 7, 2, 2) =  0.41350D-01
        HBS( 7, 3, 2) =  0.31600D-02
        HBS( 7, 4, 2) = -0.32360D-01
        HBS( 7, 5, 2) =  0.57200D-02
        HBS( 7, 6, 2) = -0.33140D-01
        HBS( 7, 7, 2) = -0.90100D-02
        HBS( 7, 8, 2) = -0.12310D-01
        HBS( 7, 9, 2) =  0.50300D-02
        HBS( 7,10, 2) =  0.46100D-01
        HBS( 7,11, 2) = -0.96000D-03
        HBS( 7,12, 2) = -0.70000D-04
        HBS( 8, 1, 2) =  0.65000D-02
        HBS( 8, 2, 2) = -0.22990D-01
        HBS( 8, 3, 2) = -0.46000D-03
        HBS( 8, 4, 2) =  0.29900D-01
        HBS( 8, 5, 2) = -0.34100D-02
        HBS( 8, 6, 2) =  0.38600D-02
        HBS( 8, 7, 2) =  0.47600D-02
        HBS( 8, 8, 2) =  0.53400D-02
        HBS( 8, 9, 2) = -0.20400D-02
        HBS( 8,10, 2) = -0.20620D-01
        HBS( 8,11, 2) = -0.19000D-03
        HBS( 8,12, 2) =  0.30000D-04
        HBS( 9, 1, 2) =  0.53200D-02
        HBS( 9, 2, 2) = -0.18360D-01
        HBS( 9, 3, 2) = -0.27000D-02
        HBS( 9, 4, 2) =  0.24600D-02
        HBS( 9, 5, 2) = -0.23100D-02
        HBS( 9, 6, 2) =  0.29280D-01
        HBS( 9, 7, 2) =  0.42500D-02
        HBS( 9, 8, 2) =  0.69700D-02
        HBS( 9, 9, 2) = -0.29900D-02
        HBS( 9,10, 2) = -0.25480D-01
        HBS( 9,11, 2) =  0.11500D-02
        HBS( 9,12, 2) =  0.40000D-04
        HBS(10, 1, 2) =  0.43830D-01
        HBS(10, 2, 2) =  0.43000D-01
        HBS(10, 3, 2) = -0.29200D-02
        HBS(10, 4, 2) = -0.13560D-01
        HBS(10, 5, 2) = -0.60000D-03
        HBS(10, 6, 2) = -0.32140D-01
        HBS(10, 7, 2) =  0.53900D-02
        HBS(10, 8, 2) = -0.31440D-01
        HBS(10, 9, 2) = -0.88000D-02
        HBS(10,10, 2) = -0.13100D-01
        HBS(10,11, 2) =  0.53000D-02
        HBS(10,12, 2) = -0.54000D-03
        HBS(11, 1, 2) = -0.24340D-01
        HBS(11, 2, 2) = -0.19110D-01
        HBS(11, 3, 2) =  0.11900D-02
        HBS(11, 4, 2) =  0.71000D-02
        HBS(11, 5, 2) =  0.10500D-02
        HBS(11, 6, 2) =  0.30020D-01
        HBS(11, 7, 2) = -0.32200D-02
        HBS(11, 8, 2) =  0.19900D-02
        HBS(11, 9, 2) =  0.44100D-02
        HBS(11,10, 2) =  0.61200D-02
        HBS(11,11, 2) = -0.21600D-02
        HBS(11,12, 2) = -0.33000D-03
        HBS(12, 1, 2) = -0.19490D-01
        HBS(12, 2, 2) = -0.23890D-01
        HBS(12, 3, 2) =  0.17300D-02
        HBS(12, 4, 2) =  0.64600D-02
        HBS(12, 5, 2) = -0.45000D-03
        HBS(12, 6, 2) =  0.21300D-02
        HBS(12, 7, 2) = -0.21700D-02
        HBS(12, 8, 2) =  0.29450D-01
        HBS(12, 9, 2) =  0.43900D-02
        HBS(12,10, 2) =  0.69800D-02
        HBS(12,11, 2) = -0.31500D-02
        HBS(12,12, 2) =  0.87000D-03
        HBS(13, 1, 2) =  0.42040D-01
        HBS(13, 2, 2) = -0.12440D-01
        HBS(13, 3, 2) =  0.58600D-02
        HBS(13, 4, 2) =  0.45620D-01
        HBS(13, 5, 2) = -0.90000D-04
        HBS(13, 6, 2) = -0.12560D-01
        HBS(13, 7, 2) = -0.10300D-02
        HBS(13, 8, 2) = -0.32810D-01
        HBS(13, 9, 2) =  0.50400D-02
        HBS(13,10, 2) = -0.33390D-01
        HBS(13,11, 2) = -0.90100D-02
        HBS(13,12, 2) =  0.50300D-02
        HBS(14, 1, 2) = -0.18640D-01
        HBS(14, 2, 2) =  0.53900D-02
        HBS(14, 3, 2) = -0.35400D-02
        HBS(14, 4, 2) = -0.25230D-01
        HBS(14, 5, 2) =  0.70000D-04
        HBS(14, 6, 2) =  0.71300D-02
        HBS(14, 7, 2) =  0.12100D-02
        HBS(14, 8, 2) =  0.29250D-01
        HBS(14, 9, 2) = -0.29800D-02
        HBS(14,10, 2) =  0.28500D-02
        HBS(14,11, 2) =  0.44100D-02
        HBS(14,12, 2) = -0.19500D-02
        HBS(15, 1, 2) = -0.23400D-01
        HBS(15, 2, 2) =  0.70500D-02
        HBS(15, 3, 2) = -0.23200D-02
        HBS(15, 4, 2) = -0.20390D-01
        HBS(15, 5, 2) =  0.20000D-04
        HBS(15, 6, 2) =  0.54300D-02
        HBS(15, 7, 2) = -0.17000D-03
        HBS(15, 8, 2) =  0.35600D-02
        HBS(15, 9, 2) = -0.20600D-02
        HBS(15,10, 2) =  0.30550D-01
        HBS(15,11, 2) =  0.46000D-02
        HBS(15,12, 2) = -0.30800D-02
        HBS(16, 1, 2) = -0.75000D-02
        HBS(16, 2, 2) = -0.28090D-01
        HBS(16, 3, 2) =  0.13300D-01
        HBS(16, 4, 2) =  0.43100D-01
        HBS(16, 5, 2) = -0.88000D-03
        HBS(16, 6, 2) =  0.45050D-01
        HBS(16, 7, 2) =  0.16000D-03
        HBS(16, 8, 2) = -0.12690D-01
        HBS(16, 9, 2) = -0.70000D-03
        HBS(16,10, 2) = -0.32780D-01
        HBS(16,11, 2) =  0.49200D-02
        HBS(16,12, 2) = -0.89700D-02
        HBS(17, 1, 2) =  0.38600D-02
        HBS(17, 2, 2) =  0.25790D-01
        HBS(17, 3, 2) = -0.82600D-02
        HBS(17, 4, 2) = -0.23900D-01
        HBS(17, 5, 2) =  0.12400D-02
        HBS(17, 6, 2) = -0.19910D-01
        HBS(17, 7, 2) = -0.90000D-04
        HBS(17, 8, 2) =  0.61400D-02
        HBS(17, 9, 2) = -0.27000D-03
        HBS(17,10, 2) =  0.32800D-02
        HBS(17,11, 2) = -0.18800D-02
        HBS(17,12, 2) =  0.70400D-02
        HBS(18, 1, 2) =  0.36300D-02
        HBS(18, 2, 2) =  0.22900D-02
        HBS(18, 3, 2) = -0.50400D-02
        HBS(18, 4, 2) = -0.19200D-01
        HBS(18, 5, 2) = -0.37000D-03
        HBS(18, 6, 2) = -0.25140D-01
        HBS(18, 7, 2) = -0.80000D-04
        HBS(18, 8, 2) =  0.65600D-02
        HBS(18, 9, 2) =  0.97000D-03
        HBS(18,10, 2) =  0.29500D-01
        HBS(18,11, 2) = -0.30300D-02
        HBS(18,12, 2) =  0.19300D-02
C Hessian elements for anchor point 3
        HSS( 1, 1, 3) =  0.39742D+00
        HSS( 2, 1, 3) =  0.15810D-01
        HSS( 2, 2, 3) =  0.37782D+00
        HSS( 3, 1, 3) =  0.44160D-01
        HSS( 3, 2, 3) =  0.54830D-01
        HSS( 3, 3, 3) =  0.38872D+00
        HSS( 4, 1, 3) =  0.19410D-01
        HSS( 4, 2, 3) = -0.16000D-02
        HSS( 4, 3, 3) = -0.34300D-02
        HSS( 4, 4, 3) =  0.41358D+00
        HSS( 5, 1, 3) =  0.23700D-02
        HSS( 5, 2, 3) = -0.70000D-03
        HSS( 5, 3, 3) =  0.14900D-02
        HSS( 5, 4, 3) =  0.48900D-02
        HSS( 5, 5, 3) =  0.36213D+00
        HSS( 6, 1, 3) = -0.43800D-02
        HSS( 6, 2, 3) =  0.80920D-01
        HSS( 6, 3, 3) = -0.20700D-02
        HSS( 6, 4, 3) =  0.21020D-01
        HSS( 6, 5, 3) = -0.70000D-03
        HSS( 6, 6, 3) =  0.41580D+00
        HSS( 7, 1, 3) = -0.56000D-03
        HSS( 7, 2, 3) = -0.41000D-03
        HSS( 7, 3, 3) =  0.64000D-03
        HSS( 7, 4, 3) =  0.44200D-02
        HSS( 7, 5, 3) =  0.38000D-03
        HSS( 7, 6, 3) =  0.45800D-02
        HSS( 7, 7, 3) =  0.36425D+00
        HSS( 8, 1, 3) =  0.80680D-01
        HSS( 8, 2, 3) = -0.46300D-02
        HSS( 8, 3, 3) = -0.10400D-02
        HSS( 8, 4, 3) = -0.39500D-02
        HSS( 8, 5, 3) =  0.11000D-03
        HSS( 8, 6, 3) =  0.21280D-01
        HSS( 8, 7, 3) = -0.32000D-03
        HSS( 8, 8, 3) =  0.41002D+00
        HSS( 9, 1, 3) = -0.40000D-03
        HSS( 9, 2, 3) = -0.40000D-03
        HSS( 9, 3, 3) = -0.55000D-03
        HSS( 9, 4, 3) = -0.66000D-03
        HSS( 9, 5, 3) =  0.17000D-03
        HSS( 9, 6, 3) =  0.44100D-02
        HSS( 9, 7, 3) =  0.44000D-03
        HSS( 9, 8, 3) =  0.43500D-02
        HSS( 9, 9, 3) =  0.36777D+00
        HSS(10, 1, 3) = -0.11100D-02
        HSS(10, 2, 3) =  0.19640D-01
        HSS(10, 3, 3) = -0.29200D-02
        HSS(10, 4, 3) =  0.86160D-01
        HSS(10, 5, 3) = -0.75000D-03
        HSS(10, 6, 3) = -0.35500D-02
        HSS(10, 7, 3) = -0.56000D-03
        HSS(10, 8, 3) =  0.21260D-01
        HSS(10, 9, 3) = -0.57000D-03
        HSS(10,10, 3) =  0.41917D+00
        HSS(11, 1, 3) = -0.46000D-03
        HSS(11, 2, 3) = -0.58000D-03
        HSS(11, 3, 3) =  0.53000D-03
        HSS(11, 4, 3) = -0.55000D-03
        HSS(11, 5, 3) =  0.90000D-04
        HSS(11, 6, 3) = -0.25000D-03
        HSS(11, 7, 3) =  0.21000D-03
        HSS(11, 8, 3) =  0.45200D-02
        HSS(11, 9, 3) =  0.43000D-03
        HSS(11,10, 3) =  0.45100D-02
        HSS(11,11, 3) =  0.36393D+00
        HSS(12, 1, 3) =  0.20000D-04
        HSS(12, 2, 3) =  0.49500D-02
        HSS(12, 3, 3) = -0.25500D-02
        HSS(12, 4, 3) = -0.86000D-03
        HSS(12, 5, 3) =  0.25000D-03
        HSS(12, 6, 3) = -0.25000D-03
        HSS(12, 7, 3) =  0.90000D-04
        HSS(12, 8, 3) = -0.82000D-03
        HSS(12, 9, 3) =  0.13000D-03
        HSS(12,10, 3) =  0.39900D-02
        HSS(12,11, 3) =  0.47000D-03
        HSS(12,12, 3) =  0.36874D+00
C Hessian elements for anchor point 3
        HBB( 1, 1, 3) =  0.82270D-01
        HBB( 2, 1, 3) = -0.48250D-01
        HBB( 2, 2, 3) =  0.13087D+00
        HBB( 3, 1, 3) = -0.34020D-01
        HBB( 3, 2, 3) = -0.82620D-01
        HBB( 3, 3, 3) =  0.11664D+00
        HBB( 4, 1, 3) = -0.33450D-01
        HBB( 4, 2, 3) =  0.20400D-01
        HBB( 4, 3, 3) =  0.13040D-01
        HBB( 4, 4, 3) =  0.77460D-01
        HBB( 5, 1, 3) =  0.17380D-01
        HBB( 5, 2, 3) = -0.84000D-02
        HBB( 5, 3, 3) = -0.89800D-02
        HBB( 5, 4, 3) = -0.37880D-01
        HBB( 5, 5, 3) =  0.75680D-01
        HBB( 6, 1, 3) =  0.16060D-01
        HBB( 6, 2, 3) = -0.12000D-01
        HBB( 6, 3, 3) = -0.40600D-02
        HBB( 6, 4, 3) = -0.39580D-01
        HBB( 6, 5, 3) = -0.37810D-01
        HBB( 6, 6, 3) =  0.77390D-01
        HBB( 7, 1, 3) = -0.10450D-01
        HBB( 7, 2, 3) =  0.14900D-02
        HBB( 7, 3, 3) =  0.89500D-02
        HBB( 7, 4, 3) = -0.31240D-01
        HBB( 7, 5, 3) =  0.12780D-01
        HBB( 7, 6, 3) =  0.18460D-01
        HBB( 7, 7, 3) =  0.78900D-01
        HBB( 8, 1, 3) = -0.47000D-03
        HBB( 8, 2, 3) =  0.30400D-02
        HBB( 8, 3, 3) = -0.25700D-02
        HBB( 8, 4, 3) =  0.18370D-01
        HBB( 8, 5, 3) = -0.69100D-02
        HBB( 8, 6, 3) = -0.11460D-01
        HBB( 8, 7, 3) = -0.39290D-01
        HBB( 8, 8, 3) =  0.78920D-01
        HBB( 9, 1, 3) =  0.10920D-01
        HBB( 9, 2, 3) = -0.45400D-02
        HBB( 9, 3, 3) = -0.63800D-02
        HBB( 9, 4, 3) =  0.12870D-01
        HBB( 9, 5, 3) = -0.58800D-02
        HBB( 9, 6, 3) = -0.70000D-02
        HBB( 9, 7, 3) = -0.39610D-01
        HBB( 9, 8, 3) = -0.39630D-01
        HBB( 9, 9, 3) =  0.79240D-01
        HBB(10, 1, 3) =  0.65500D-02
        HBB(10, 2, 3) = -0.46300D-02
        HBB(10, 3, 3) = -0.19200D-02
        HBB(10, 4, 3) = -0.92300D-02
        HBB(10, 5, 3) =  0.10020D-01
        HBB(10, 6, 3) = -0.79000D-03
        HBB(10, 7, 3) = -0.34090D-01
        HBB(10, 8, 3) =  0.13830D-01
        HBB(10, 9, 3) =  0.20260D-01
        HBB(10,10, 3) =  0.80030D-01
        HBB(11, 1, 3) = -0.33400D-02
        HBB(11, 2, 3) =  0.27100D-02
        HBB(11, 3, 3) =  0.63000D-03
        HBB(11, 4, 3) = -0.98000D-03
        HBB(11, 5, 3) = -0.35000D-02
        HBB(11, 6, 3) =  0.44800D-02
        HBB(11, 7, 3) =  0.19980D-01
        HBB(11, 8, 3) = -0.73600D-02
        HBB(11, 9, 3) = -0.12620D-01
        HBB(11,10, 3) = -0.40030D-01
        HBB(11,11, 3) =  0.78930D-01
        HBB(12, 1, 3) = -0.32100D-02
        HBB(12, 2, 3) =  0.19200D-02
        HBB(12, 3, 3) =  0.13000D-02
        HBB(12, 4, 3) =  0.10210D-01
        HBB(12, 5, 3) = -0.65300D-02
        HBB(12, 6, 3) = -0.36800D-02
        HBB(12, 7, 3) =  0.14110D-01
        HBB(12, 8, 3) = -0.64600D-02
        HBB(12, 9, 3) = -0.76400D-02
        HBB(12,10, 3) = -0.40000D-01
        HBB(12,11, 3) = -0.38900D-01
        HBB(12,12, 3) =  0.78900D-01
        HBB(13, 1, 3) = -0.96900D-02
        HBB(13, 2, 3) =  0.10920D-01
        HBB(13, 3, 3) = -0.12400D-02
        HBB(13, 4, 3) =  0.49700D-02
        HBB(13, 5, 3) = -0.18200D-02
        HBB(13, 6, 3) = -0.31500D-02
        HBB(13, 7, 3) = -0.88100D-02
        HBB(13, 8, 3) =  0.10200D-01
        HBB(13, 9, 3) = -0.13900D-02
        HBB(13,10, 3) = -0.34270D-01
        HBB(13,11, 3) =  0.14180D-01
        HBB(13,12, 3) =  0.20100D-01
        HBB(13,13, 3) =  0.79280D-01
        HBB(14, 1, 3) =  0.10380D-01
        HBB(14, 2, 3) = -0.74600D-02
        HBB(14, 3, 3) = -0.29200D-02
        HBB(14, 4, 3) = -0.26400D-02
        HBB(14, 5, 3) =  0.92000D-03
        HBB(14, 6, 3) =  0.17200D-02
        HBB(14, 7, 3) = -0.14500D-02
        HBB(14, 8, 3) = -0.35600D-02
        HBB(14, 9, 3) =  0.50100D-02
        HBB(14,10, 3) =  0.20210D-01
        HBB(14,11, 3) = -0.76400D-02
        HBB(14,12, 3) = -0.12580D-01
        HBB(14,13, 3) = -0.39700D-01
        HBB(14,14, 3) =  0.79390D-01
        HBB(15, 1, 3) = -0.69000D-03
        HBB(15, 2, 3) = -0.34600D-02
        HBB(15, 3, 3) =  0.41500D-02
        HBB(15, 4, 3) = -0.23300D-02
        HBB(15, 5, 3) =  0.91000D-03
        HBB(15, 6, 3) =  0.14300D-02
        HBB(15, 7, 3) =  0.10260D-01
        HBB(15, 8, 3) = -0.66500D-02
        HBB(15, 9, 3) = -0.36100D-02
        HBB(15,10, 3) =  0.14060D-01
        HBB(15,11, 3) = -0.65400D-02
        HBB(15,12, 3) = -0.75200D-02
        HBB(15,13, 3) = -0.39580D-01
        HBB(15,14, 3) = -0.39690D-01
        HBB(15,15, 3) =  0.79270D-01
        HBB(16, 1, 3) = -0.35240D-01
        HBB(16, 2, 3) =  0.20060D-01
        HBB(16, 3, 3) =  0.15180D-01
        HBB(16, 4, 3) = -0.85100D-02
        HBB(16, 5, 3) = -0.49000D-03
        HBB(16, 6, 3) =  0.90000D-02
        HBB(16, 7, 3) =  0.56800D-02
        HBB(16, 8, 3) = -0.26300D-02
        HBB(16, 9, 3) = -0.30500D-02
        HBB(16,10, 3) = -0.90000D-02
        HBB(16,11, 3) =  0.10190D-01
        HBB(16,12, 3) = -0.12000D-02
        HBB(16,13, 3) = -0.31490D-01
        HBB(16,14, 3) =  0.13200D-01
        HBB(16,15, 3) =  0.18290D-01
        HBB(16,16, 3) =  0.78560D-01
        HBB(17, 1, 3) =  0.20310D-01
        HBB(17, 2, 3) = -0.11130D-01
        HBB(17, 3, 3) = -0.91900D-02
        HBB(17, 4, 3) = -0.11900D-02
        HBB(17, 5, 3) =  0.41100D-02
        HBB(17, 6, 3) = -0.29200D-02
        HBB(17, 7, 3) = -0.26100D-02
        HBB(17, 8, 3) =  0.13200D-02
        HBB(17, 9, 3) =  0.12900D-02
        HBB(17,10, 3) =  0.10050D-01
        HBB(17,11, 3) = -0.66300D-02
        HBB(17,12, 3) = -0.34200D-02
        HBB(17,13, 3) =  0.12780D-01
        HBB(17,14, 3) = -0.59900D-02
        HBB(17,15, 3) = -0.67900D-02
        HBB(17,16, 3) = -0.39350D-01
        HBB(17,17, 3) =  0.77640D-01
        HBB(18, 1, 3) =  0.14930D-01
        HBB(18, 2, 3) = -0.89300D-02
        HBB(18, 3, 3) = -0.60000D-02
        HBB(18, 4, 3) =  0.97000D-02
        HBB(18, 5, 3) = -0.36200D-02
        HBB(18, 6, 3) = -0.60800D-02
        HBB(18, 7, 3) = -0.30800D-02
        HBB(18, 8, 3) =  0.13100D-02
        HBB(18, 9, 3) =  0.17700D-02
        HBB(18,10, 3) = -0.10500D-02
        HBB(18,11, 3) = -0.35600D-02
        HBB(18,12, 3) =  0.46100D-02
        HBB(18,13, 3) =  0.18710D-01
        HBB(18,14, 3) = -0.72100D-02
        HBB(18,15, 3) = -0.11500D-01
        HBB(18,16, 3) = -0.39210D-01
        HBB(18,17, 3) = -0.38290D-01
        HBB(18,18, 3) =  0.77500D-01
C Hessian elements for anchor point 3
        HTT( 1, 1, 3) =  0.10290D-01
        HTT( 2, 1, 3) =  0.63900D-02
        HTT( 2, 2, 3) =  0.13000D-01
        HTT( 3, 1, 3) =  0.28000D-03
        HTT( 3, 2, 3) = -0.46200D-02
        HTT( 3, 3, 3) =  0.11260D-01
        HTT( 4, 1, 3) = -0.36100D-02
        HTT( 4, 2, 3) =  0.19900D-02
        HTT( 4, 3, 3) =  0.63700D-02
        HTT( 4, 4, 3) =  0.11970D-01
        HTT( 5, 1, 3) = -0.93800D-02
        HTT( 5, 2, 3) = -0.10580D-01
        HTT( 5, 3, 3) =  0.29200D-02
        HTT( 5, 4, 3) =  0.17200D-02
        HTT( 5, 5, 3) =  0.16730D-01
        HTT( 6, 1, 3) = -0.54800D-02
        HTT( 6, 2, 3) = -0.54200D-02
        HTT( 6, 3, 3) =  0.33200D-02
        HTT( 6, 4, 3) =  0.33700D-02
        HTT( 6, 5, 3) =  0.66800D-02
        HTT( 6, 6, 3) =  0.91600D-02
        HTT( 7, 1, 3) = -0.40000D-04
        HTT( 7, 2, 3) = -0.31000D-03
        HTT( 7, 3, 3) = -0.73300D-02
        HTT( 7, 4, 3) = -0.75900D-02
        HTT( 7, 5, 3) =  0.52600D-02
        HTT( 7, 6, 3) = -0.15300D-02
        HTT( 7, 7, 3) =  0.12050D-01
        HTT( 8, 1, 3) =  0.38600D-02
        HTT( 8, 2, 3) =  0.48500D-02
        HTT( 8, 3, 3) = -0.69300D-02
        HTT( 8, 4, 3) = -0.59400D-02
        HTT( 8, 5, 3) = -0.48000D-02
        HTT( 8, 6, 3) =  0.95000D-03
        HTT( 8, 7, 3) =  0.52700D-02
        HTT( 8, 8, 3) =  0.11020D-01
        HTT( 9, 1, 3) = -0.52600D-02
        HTT( 9, 2, 3) =  0.85000D-03
        HTT( 9, 3, 3) = -0.32600D-02
        HTT( 9, 4, 3) =  0.28400D-02
        HTT( 9, 5, 3) = -0.16300D-02
        HTT( 9, 6, 3) =  0.11100D-02
        HTT( 9, 7, 3) = -0.35000D-02
        HTT( 9, 8, 3) = -0.76000D-03
        HTT( 9, 9, 3) =  0.98000D-02
        HTT(10, 1, 3) = -0.32000D-02
        HTT(10, 2, 3) =  0.20800D-02
        HTT(10, 3, 3) = -0.19100D-02
        HTT(10, 4, 3) =  0.33700D-02
        HTT(10, 5, 3) =  0.37000D-03
        HTT(10, 6, 3) =  0.79000D-03
        HTT(10, 7, 3) = -0.83000D-03
        HTT(10, 8, 3) = -0.41000D-03
        HTT(10, 9, 3) =  0.34100D-02
        HTT(10,10, 3) =  0.86300D-02
        HTT(11, 1, 3) = -0.12800D-02
        HTT(11, 2, 3) = -0.59200D-02
        HTT(11, 3, 3) =  0.17500D-02
        HTT(11, 4, 3) = -0.28900D-02
        HTT(11, 5, 3) = -0.41000D-03
        HTT(11, 6, 3) =  0.10500D-02
        HTT(11, 7, 3) = -0.32300D-02
        HTT(11, 8, 3) = -0.17700D-02
        HTT(11, 9, 3) =  0.35500D-02
        HTT(11,10, 3) = -0.20000D-02
        HTT(11,11, 3) =  0.83000D-02
        HTT(12, 1, 3) =  0.78000D-03
        HTT(12, 2, 3) = -0.46800D-02
        HTT(12, 3, 3) =  0.31000D-02
        HTT(12, 4, 3) = -0.23600D-02
        HTT(12, 5, 3) =  0.16000D-02
        HTT(12, 6, 3) =  0.73000D-03
        HTT(12, 7, 3) = -0.56000D-03
        HTT(12, 8, 3) = -0.14300D-02
        HTT(12, 9, 3) = -0.28400D-02
        HTT(12,10, 3) =  0.32200D-02
        HTT(12,11, 3) =  0.27500D-02
        HTT(12,12, 3) =  0.88100D-02
        HTT(13, 1, 3) = -0.66000D-03
        HTT(13, 2, 3) = -0.37900D-02
        HTT(13, 3, 3) =  0.29600D-02
        HTT(13, 4, 3) = -0.16000D-03
        HTT(13, 5, 3) =  0.51200D-02
        HTT(13, 6, 3) =  0.20100D-02
        HTT(13, 7, 3) =  0.17400D-02
        HTT(13, 8, 3) = -0.13700D-02
        HTT(13, 9, 3) = -0.72600D-02
        HTT(13,10, 3) = -0.74000D-03
        HTT(13,11, 3) = -0.40600D-02
        HTT(13,12, 3) =  0.24600D-02
        HTT(13,13, 3) =  0.10500D-01
        HTT(14, 1, 3) =  0.16400D-02
        HTT(14, 2, 3) = -0.60000D-04
        HTT(14, 3, 3) =  0.14800D-02
        HTT(14, 4, 3) = -0.22000D-03
        HTT(14, 5, 3) =  0.17400D-02
        HTT(14, 6, 3) = -0.85000D-03
        HTT(14, 7, 3) =  0.18900D-02
        HTT(14, 8, 3) = -0.70000D-03
        HTT(14, 9, 3) = -0.45800D-02
        HTT(14,10, 3) =  0.11000D-02
        HTT(14,11, 3) = -0.28400D-02
        HTT(14,12, 3) =  0.28400D-02
        HTT(14,13, 3) =  0.40400D-02
        HTT(14,14, 3) =  0.85900D-02
        HTT(15, 1, 3) = -0.27400D-02
        HTT(15, 2, 3) = -0.50400D-02
        HTT(15, 3, 3) =  0.16000D-02
        HTT(15, 4, 3) = -0.70000D-03
        HTT(15, 5, 3) =  0.30800D-02
        HTT(15, 6, 3) =  0.23300D-02
        HTT(15, 7, 3) = -0.97000D-03
        HTT(15, 8, 3) = -0.17200D-02
        HTT(15, 9, 3) = -0.78000D-03
        HTT(15,10, 3) = -0.60300D-02
        HTT(15,11, 3) =  0.15700D-02
        HTT(15,12, 3) = -0.36800D-02
        HTT(15,13, 3) =  0.39000D-02
        HTT(15,14, 3) = -0.17100D-02
        HTT(15,15, 3) =  0.92100D-02
        HTT(16, 1, 3) = -0.44000D-03
        HTT(16, 2, 3) = -0.13100D-02
        HTT(16, 3, 3) =  0.11000D-03
        HTT(16, 4, 3) = -0.76000D-03
        HTT(16, 5, 3) = -0.29000D-03
        HTT(16, 6, 3) = -0.53000D-03
        HTT(16, 7, 3) = -0.81000D-03
        HTT(16, 8, 3) = -0.10500D-02
        HTT(16, 9, 3) =  0.18900D-02
        HTT(16,10, 3) = -0.41800D-02
        HTT(16,11, 3) =  0.27800D-02
        HTT(16,12, 3) = -0.33000D-02
        HTT(16,13, 3) = -0.25500D-02
        HTT(16,14, 3) =  0.28300D-02
        HTT(16,15, 3) =  0.36000D-02
        HTT(16,16, 3) =  0.89900D-02
        HTT(17, 1, 3) =  0.15600D-02
        HTT(17, 2, 3) = -0.56000D-03
        HTT(17, 3, 3) =  0.34000D-03
        HTT(17, 4, 3) = -0.17800D-02
        HTT(17, 5, 3) =  0.25100D-02
        HTT(17, 6, 3) = -0.75000D-03
        HTT(17, 7, 3) =  0.36600D-02
        HTT(17, 8, 3) =  0.39000D-03
        HTT(17, 9, 3) =  0.16000D-03
        HTT(17,10, 3) = -0.21700D-02
        HTT(17,11, 3) =  0.23300D-02
        HTT(17,12, 3) =  0.00000D+00
        HTT(17,13, 3) = -0.59100D-02
        HTT(17,14, 3) = -0.56000D-03
        HTT(17,15, 3) = -0.35500D-02
        HTT(17,16, 3) =  0.18000D-02
        HTT(17,17, 3) =  0.10140D-01
        HTT(18, 1, 3) =  0.13300D-02
        HTT(18, 2, 3) =  0.25800D-02
        HTT(18, 3, 3) = -0.22700D-02
        HTT(18, 4, 3) = -0.10200D-02
        HTT(18, 5, 3) = -0.30300D-02
        HTT(18, 6, 3) = -0.35800D-02
        HTT(18, 7, 3) =  0.33000D-03
        HTT(18, 8, 3) = -0.23000D-03
        HTT(18, 9, 3) =  0.24700D-02
        HTT(18,10, 3) =  0.75000D-03
        HTT(18,11, 3) =  0.12000D-02
        HTT(18,12, 3) = -0.52000D-03
        HTT(18,13, 3) = -0.44900D-02
        HTT(18,14, 3) =  0.90000D-03
        HTT(18,15, 3) = -0.27500D-02
        HTT(18,16, 3) =  0.26400D-02
        HTT(18,17, 3) =  0.27400D-02
        HTT(18,18, 3) =  0.87100D-02
        HTT(19, 1, 3) = -0.74000D-03
        HTT(19, 2, 3) = -0.42800D-02
        HTT(19, 3, 3) =  0.18200D-02
        HTT(19, 4, 3) = -0.17200D-02
        HTT(19, 5, 3) =  0.58900D-02
        HTT(19, 6, 3) =  0.21100D-02
        HTT(19, 7, 3) =  0.35000D-02
        HTT(19, 8, 3) = -0.28000D-03
        HTT(19, 9, 3) = -0.25100D-02
        HTT(19,10, 3) = -0.40200D-02
        HTT(19,11, 3) =  0.11200D-02
        HTT(19,12, 3) = -0.38000D-03
        HTT(19,13, 3) =  0.54000D-03
        HTT(19,14, 3) = -0.51100D-02
        HTT(19,15, 3) =  0.20600D-02
        HTT(19,16, 3) = -0.35900D-02
        HTT(19,17, 3) =  0.47900D-02
        HTT(19,18, 3) = -0.26500D-02
        HTT(19,19, 3) =  0.10440D-01
        HTT(20, 1, 3) = -0.97000D-03
        HTT(20, 2, 3) = -0.11500D-02
        HTT(20, 3, 3) = -0.78000D-03
        HTT(20, 4, 3) = -0.96000D-03
        HTT(20, 5, 3) =  0.35000D-03
        HTT(20, 6, 3) = -0.72000D-03
        HTT(20, 7, 3) =  0.18000D-03
        HTT(20, 8, 3) = -0.89000D-03
        HTT(20, 9, 3) = -0.20000D-03
        HTT(20,10, 3) = -0.10900D-02
        HTT(20,11, 3) = -0.20000D-04
        HTT(20,12, 3) = -0.90000D-03
        HTT(20,13, 3) =  0.19600D-02
        HTT(20,14, 3) = -0.36500D-02
        HTT(20,15, 3) =  0.28600D-02
        HTT(20,16, 3) = -0.27500D-02
        HTT(20,17, 3) = -0.26100D-02
        HTT(20,18, 3) =  0.33200D-02
        HTT(20,19, 3) =  0.30100D-02
        HTT(20,20, 3) =  0.89300D-02
        HTT(21, 1, 3) =  0.34400D-02
        HTT(21, 2, 3) =  0.76900D-02
        HTT(21, 3, 3) = -0.32500D-02
        HTT(21, 4, 3) =  0.10000D-02
        HTT(21, 5, 3) = -0.13360D-01
        HTT(21, 6, 3) = -0.35700D-02
        HTT(21, 7, 3) = -0.71200D-02
        HTT(21, 8, 3) =  0.26800D-02
        HTT(21, 9, 3) =  0.42100D-02
        HTT(21,10, 3) =  0.23400D-02
        HTT(21,11, 3) = -0.13000D-03
        HTT(21,12, 3) = -0.20000D-02
        HTT(21,13, 3) = -0.18100D-02
        HTT(21,14, 3) = -0.22900D-02
        HTT(21,15, 3) =  0.90000D-04
        HTT(21,16, 3) = -0.39000D-03
        HTT(21,17, 3) = -0.84600D-02
        HTT(21,18, 3) =  0.98000D-03
        HTT(21,19, 3) = -0.79800D-02
        HTT(21,20, 3) =  0.14600D-02
        HTT(21,21, 3) =  0.15990D-01
        HTT(22, 1, 3) = -0.55000D-03
        HTT(22, 2, 3) =  0.24000D-02
        HTT(22, 3, 3) = -0.36500D-02
        HTT(22, 4, 3) = -0.70000D-03
        HTT(22, 5, 3) = -0.30600D-02
        HTT(22, 6, 3) = -0.61100D-02
        HTT(22, 7, 3) = -0.17000D-03
        HTT(22, 8, 3) = -0.32100D-02
        HTT(22, 9, 3) =  0.14100D-02
        HTT(22,10, 3) =  0.19200D-02
        HTT(22,11, 3) = -0.16200D-02
        HTT(22,12, 3) = -0.11100D-02
        HTT(22,13, 3) =  0.13700D-02
        HTT(22,14, 3) =  0.37000D-03
        HTT(22,15, 3) =  0.86000D-03
        HTT(22,16, 3) = -0.15000D-03
        HTT(22,17, 3) = -0.51200D-02
        HTT(22,18, 3) =  0.15500D-02
        HTT(22,19, 3) = -0.41100D-02
        HTT(22,20, 3) =  0.25600D-02
        HTT(22,21, 3) =  0.59500D-02
        HTT(22,22, 3) =  0.90700D-02
        HTT(23, 1, 3) =  0.36700D-02
        HTT(23, 2, 3) =  0.46000D-02
        HTT(23, 3, 3) = -0.68000D-03
        HTT(23, 4, 3) =  0.25000D-03
        HTT(23, 5, 3) = -0.79000D-02
        HTT(23, 6, 3) = -0.77000D-03
        HTT(23, 7, 3) = -0.38400D-02
        HTT(23, 8, 3) =  0.32800D-02
        HTT(23, 9, 3) =  0.19400D-02
        HTT(23,10, 3) = -0.54000D-03
        HTT(23,11, 3) =  0.99000D-03
        HTT(23,12, 3) = -0.14900D-02
        HTT(23,13, 3) = -0.32100D-02
        HTT(23,14, 3) = -0.37300D-02
        HTT(23,15, 3) = -0.70000D-03
        HTT(23,16, 3) = -0.12200D-02
        HTT(23,17, 3) = -0.11700D-02
        HTT(23,18, 3) = -0.49000D-02
        HTT(23,19, 3) = -0.65000D-03
        HTT(23,20, 3) = -0.43800D-02
        HTT(23,21, 3) =  0.66800D-02
        HTT(23,22, 3) = -0.62000D-03
        HTT(23,23, 3) =  0.10360D-01
        HTT(24, 1, 3) = -0.33000D-03
        HTT(24, 2, 3) = -0.69000D-03
        HTT(24, 3, 3) = -0.10800D-02
        HTT(24, 4, 3) = -0.14500D-02
        HTT(24, 5, 3) =  0.24000D-02
        HTT(24, 6, 3) = -0.33100D-02
        HTT(24, 7, 3) =  0.31100D-02
        HTT(24, 8, 3) = -0.26100D-02
        HTT(24, 9, 3) = -0.87000D-03
        HTT(24,10, 3) = -0.97000D-03
        HTT(24,11, 3) = -0.50000D-03
        HTT(24,12, 3) = -0.60000D-03
        HTT(24,13, 3) = -0.30000D-04
        HTT(24,14, 3) = -0.10800D-02
        HTT(24,15, 3) =  0.70000D-04
        HTT(24,16, 3) = -0.98000D-03
        HTT(24,17, 3) =  0.21800D-02
        HTT(24,18, 3) = -0.43300D-02
        HTT(24,19, 3) =  0.32200D-02
        HTT(24,20, 3) = -0.32900D-02
        HTT(24,21, 3) = -0.33600D-02
        HTT(24,22, 3) =  0.25000D-02
        HTT(24,23, 3) =  0.30700D-02
        HTT(24,24, 3) =  0.89200D-02
C Hessian elements for anchor point 3
        HBS( 1, 1, 3) = -0.35620D-01
        HBS( 1, 2, 3) = -0.31870D-01
        HBS( 1, 3, 3) = -0.41940D-01
        HBS( 1, 4, 3) = -0.11410D-01
        HBS( 1, 5, 3) =  0.68300D-02
        HBS( 1, 6, 3) =  0.45490D-01
        HBS( 1, 7, 3) = -0.31000D-03
        HBS( 1, 8, 3) =  0.44930D-01
        HBS( 1, 9, 3) =  0.00000D+00
        HBS( 1,10, 3) = -0.10110D-01
        HBS( 1,11, 3) = -0.33000D-03
        HBS( 1,12, 3) =  0.43300D-02
        HBS( 2, 1, 3) =  0.49480D-01
        HBS( 2, 2, 3) = -0.69900D-02
        HBS( 2, 3, 3) =  0.24890D-01
        HBS( 2, 4, 3) = -0.87000D-03
        HBS( 2, 5, 3) = -0.66200D-02
        HBS( 2, 6, 3) = -0.20260D-01
        HBS( 2, 7, 3) =  0.89000D-03
        HBS( 2, 8, 3) = -0.28090D-01
        HBS( 2, 9, 3) = -0.10000D-04
        HBS( 2,10, 3) =  0.10450D-01
        HBS( 2,11, 3) = -0.46000D-03
        HBS( 2,12, 3) = -0.19500D-02
        HBS( 3, 1, 3) = -0.13860D-01
        HBS( 3, 2, 3) =  0.38860D-01
        HBS( 3, 3, 3) =  0.17050D-01
        HBS( 3, 4, 3) =  0.12280D-01
        HBS( 3, 5, 3) = -0.21000D-03
        HBS( 3, 6, 3) = -0.25230D-01
        HBS( 3, 7, 3) = -0.58000D-03
        HBS( 3, 8, 3) = -0.16840D-01
        HBS( 3, 9, 3) =  0.20000D-04
        HBS( 3,10, 3) = -0.34000D-03
        HBS( 3,11, 3) =  0.79000D-03
        HBS( 3,12, 3) = -0.23700D-02
        HBS( 4, 1, 3) = -0.28260D-01
        HBS( 4, 2, 3) = -0.99500D-02
        HBS( 4, 3, 3) =  0.16700D-01
        HBS( 4, 4, 3) = -0.32020D-01
        HBS( 4, 5, 3) = -0.96100D-02
        HBS( 4, 6, 3) = -0.12920D-01
        HBS( 4, 7, 3) =  0.48000D-02
        HBS( 4, 8, 3) =  0.45530D-01
        HBS( 4, 9, 3) = -0.62000D-03
        HBS( 4,10, 3) =  0.44050D-01
        HBS( 4,11, 3) =  0.10000D-03
        HBS( 4,12, 3) = -0.40000D-03
        HBS( 5, 1, 3) =  0.30480D-01
        HBS( 5, 2, 3) =  0.69800D-02
        HBS( 5, 3, 3) = -0.14860D-01
        HBS( 5, 4, 3) =  0.26800D-02
        HBS( 5, 5, 3) =  0.62000D-02
        HBS( 5, 6, 3) =  0.57700D-02
        HBS( 5, 7, 3) = -0.18500D-02
        HBS( 5, 8, 3) = -0.20690D-01
        HBS( 5, 9, 3) = -0.33000D-03
        HBS( 5,10, 3) = -0.24600D-01
        HBS( 5,11, 3) = -0.60000D-04
        HBS( 5,12, 3) =  0.44000D-03
        HBS( 6, 1, 3) = -0.22300D-02
        HBS( 6, 2, 3) =  0.29700D-02
        HBS( 6, 3, 3) = -0.18400D-02
        HBS( 6, 4, 3) =  0.29340D-01
        HBS( 6, 5, 3) =  0.34100D-02
        HBS( 6, 6, 3) =  0.71600D-02
        HBS( 6, 7, 3) = -0.29500D-02
        HBS( 6, 8, 3) = -0.24830D-01
        HBS( 6, 9, 3) =  0.94000D-03
        HBS( 6,10, 3) = -0.19440D-01
        HBS( 6,11, 3) = -0.40000D-04
        HBS( 6,12, 3) = -0.30000D-04
        HBS( 7, 1, 3) = -0.12340D-01
        HBS( 7, 2, 3) =  0.39640D-01
        HBS( 7, 3, 3) =  0.73200D-02
        HBS( 7, 4, 3) = -0.32120D-01
        HBS( 7, 5, 3) =  0.49300D-02
        HBS( 7, 6, 3) = -0.33370D-01
        HBS( 7, 7, 3) = -0.89700D-02
        HBS( 7, 8, 3) = -0.12760D-01
        HBS( 7, 9, 3) =  0.49400D-02
        HBS( 7,10, 3) =  0.45710D-01
        HBS( 7,11, 3) = -0.89000D-03
        HBS( 7,12, 3) =  0.50000D-04
        HBS( 8, 1, 3) =  0.73400D-02
        HBS( 8, 2, 3) = -0.22220D-01
        HBS( 8, 3, 3) = -0.35200D-02
        HBS( 8, 4, 3) =  0.29850D-01
        HBS( 8, 5, 3) = -0.30000D-02
        HBS( 8, 6, 3) =  0.38900D-02
        HBS( 8, 7, 3) =  0.48800D-02
        HBS( 8, 8, 3) =  0.56300D-02
        HBS( 8, 9, 3) = -0.19800D-02
        HBS( 8,10, 3) = -0.20400D-01
        HBS( 8,11, 3) = -0.22000D-03
        HBS( 8,12, 3) = -0.50000D-04
        HBS( 9, 1, 3) =  0.50100D-02
        HBS( 9, 2, 3) = -0.17420D-01
        HBS( 9, 3, 3) = -0.38000D-02
        HBS( 9, 4, 3) =  0.22700D-02
        HBS( 9, 5, 3) = -0.19200D-02
        HBS( 9, 6, 3) =  0.29480D-01
        HBS( 9, 7, 3) =  0.40900D-02
        HBS( 9, 8, 3) =  0.71400D-02
        HBS( 9, 9, 3) = -0.29600D-02
        HBS( 9,10, 3) = -0.25300D-01
        HBS( 9,11, 3) =  0.11100D-02
        HBS( 9,12, 3) = -0.10000D-04
        HBS(10, 1, 3) =  0.42740D-01
        HBS(10, 2, 3) =  0.42010D-01
        HBS(10, 3, 3) = -0.32900D-02
        HBS(10, 4, 3) = -0.13800D-01
        HBS(10, 5, 3) = -0.37000D-03
        HBS(10, 6, 3) = -0.32450D-01
        HBS(10, 7, 3) =  0.53400D-02
        HBS(10, 8, 3) = -0.31850D-01
        HBS(10, 9, 3) = -0.86200D-02
        HBS(10,10, 3) = -0.13500D-01
        HBS(10,11, 3) =  0.52600D-02
        HBS(10,12, 3) = -0.66000D-03
        HBS(11, 1, 3) = -0.23640D-01
        HBS(11, 2, 3) = -0.18640D-01
        HBS(11, 3, 3) =  0.15000D-02
        HBS(11, 4, 3) =  0.71800D-02
        HBS(11, 5, 3) =  0.84000D-03
        HBS(11, 6, 3) =  0.30210D-01
        HBS(11, 7, 3) = -0.32100D-02
        HBS(11, 8, 3) =  0.20600D-02
        HBS(11, 9, 3) =  0.42700D-02
        HBS(11,10, 3) =  0.63000D-02
        HBS(11,11, 3) = -0.21200D-02
        HBS(11,12, 3) = -0.27000D-03
        HBS(12, 1, 3) = -0.19100D-01
        HBS(12, 2, 3) = -0.23380D-01
        HBS(12, 3, 3) =  0.17900D-02
        HBS(12, 4, 3) =  0.66100D-02
        HBS(12, 5, 3) = -0.47000D-03
        HBS(12, 6, 3) =  0.22400D-02
        HBS(12, 7, 3) = -0.21200D-02
        HBS(12, 8, 3) =  0.29790D-01
        HBS(12, 9, 3) =  0.43500D-02
        HBS(12,10, 3) =  0.72000D-02
        HBS(12,11, 3) = -0.31400D-02
        HBS(12,12, 3) =  0.93000D-03
        HBS(13, 1, 3) =  0.40850D-01
        HBS(13, 2, 3) = -0.14500D-01
        HBS(13, 3, 3) =  0.10130D-01
        HBS(13, 4, 3) =  0.45630D-01
        HBS(13, 5, 3) = -0.20000D-03
        HBS(13, 6, 3) = -0.12730D-01
        HBS(13, 7, 3) = -0.96000D-03
        HBS(13, 8, 3) = -0.32840D-01
        HBS(13, 9, 3) =  0.49400D-02
        HBS(13,10, 3) = -0.32520D-01
        HBS(13,11, 3) = -0.89800D-02
        HBS(13,12, 3) =  0.51000D-02
        HBS(14, 1, 3) = -0.17990D-01
        HBS(14, 2, 3) =  0.57300D-02
        HBS(14, 3, 3) = -0.40700D-02
        HBS(14, 4, 3) = -0.25200D-01
        HBS(14, 5, 3) =  0.12000D-03
        HBS(14, 6, 3) =  0.71400D-02
        HBS(14, 7, 3) =  0.11600D-02
        HBS(14, 8, 3) =  0.29210D-01
        HBS(14, 9, 3) = -0.29300D-02
        HBS(14,10, 3) =  0.21300D-02
        HBS(14,11, 3) =  0.42700D-02
        HBS(14,12, 3) = -0.19500D-02
        HBS(15, 1, 3) = -0.22870D-01
        HBS(15, 2, 3) =  0.87800D-02
        HBS(15, 3, 3) = -0.60600D-02
        HBS(15, 4, 3) = -0.20420D-01
        HBS(15, 5, 3) =  0.70000D-04
        HBS(15, 6, 3) =  0.55900D-02
        HBS(15, 7, 3) = -0.20000D-03
        HBS(15, 8, 3) =  0.36300D-02
        HBS(15, 9, 3) = -0.20100D-02
        HBS(15,10, 3) =  0.30390D-01
        HBS(15,11, 3) =  0.47100D-02
        HBS(15,12, 3) = -0.31500D-02
        HBS(16, 1, 3) = -0.73800D-02
        HBS(16, 2, 3) = -0.25320D-01
        HBS(16, 3, 3) =  0.11080D-01
        HBS(16, 4, 3) =  0.43720D-01
        HBS(16, 5, 3) = -0.15800D-02
        HBS(16, 6, 3) =  0.45980D-01
        HBS(16, 7, 3) =  0.11000D-03
        HBS(16, 8, 3) = -0.13000D-01
        HBS(16, 9, 3) = -0.64000D-03
        HBS(16,10, 3) = -0.33620D-01
        HBS(16,11, 3) =  0.48400D-02
        HBS(16,12, 3) = -0.84200D-02
        HBS(17, 1, 3) =  0.40600D-02
        HBS(17, 2, 3) =  0.27060D-01
        HBS(17, 3, 3) = -0.92100D-02
        HBS(17, 4, 3) = -0.24320D-01
        HBS(17, 5, 3) =  0.15600D-02
        HBS(17, 6, 3) = -0.20640D-01
        HBS(17, 7, 3) = -0.50000D-04
        HBS(17, 8, 3) =  0.61500D-02
        HBS(17, 9, 3) = -0.30000D-03
        HBS(17,10, 3) =  0.36700D-02
        HBS(17,11, 3) = -0.18400D-02
        HBS(17,12, 3) =  0.60700D-02
        HBS(18, 1, 3) =  0.33200D-02
        HBS(18, 2, 3) = -0.17400D-02
        HBS(18, 3, 3) = -0.18800D-02
        HBS(18, 4, 3) = -0.19400D-01
        HBS(18, 5, 3) =  0.20000D-04
        HBS(18, 6, 3) = -0.25340D-01
        HBS(18, 7, 3) = -0.50000D-04
        HBS(18, 8, 3) =  0.68600D-02
        HBS(18, 9, 3) =  0.93000D-03
        HBS(18,10, 3) =  0.29960D-01
        HBS(18,11, 3) = -0.30000D-02
        HBS(18,12, 3) =  0.23500D-02
C Hessian elements for anchor point 4
        HSS( 1, 1, 4) =  0.38571D+00
        HSS( 2, 1, 4) =  0.12420D-01
        HSS( 2, 2, 4) =  0.38571D+00
        HSS( 3, 1, 4) =  0.51260D-01
        HSS( 3, 2, 4) =  0.51260D-01
        HSS( 3, 3, 4) =  0.38504D+00
        HSS( 4, 1, 4) =  0.20570D-01
        HSS( 4, 2, 4) = -0.17300D-02
        HSS( 4, 3, 4) = -0.32200D-02
        HSS( 4, 4, 4) =  0.41668D+00
        HSS( 5, 1, 4) =  0.50700D-02
        HSS( 5, 2, 4) =  0.48000D-03
        HSS( 5, 3, 4) = -0.24900D-02
        HSS( 5, 4, 4) =  0.39900D-02
        HSS( 5, 5, 4) =  0.36905D+00
        HSS( 6, 1, 4) = -0.36300D-02
        HSS( 6, 2, 4) =  0.82190D-01
        HSS( 6, 3, 4) = -0.13300D-02
        HSS( 6, 4, 4) =  0.21260D-01
        HSS( 6, 5, 4) = -0.91000D-03
        HSS( 6, 6, 4) =  0.41243D+00
        HSS( 7, 1, 4) = -0.51000D-03
        HSS( 7, 2, 4) = -0.39000D-03
        HSS( 7, 3, 4) =  0.41000D-03
        HSS( 7, 4, 4) =  0.44700D-02
        HSS( 7, 5, 4) =  0.43000D-03
        HSS( 7, 6, 4) =  0.45000D-02
        HSS( 7, 7, 4) =  0.36452D+00
        HSS( 8, 1, 4) =  0.82190D-01
        HSS( 8, 2, 4) = -0.36300D-02
        HSS( 8, 3, 4) = -0.13300D-02
        HSS( 8, 4, 4) = -0.41800D-02
        HSS( 8, 5, 4) = -0.36000D-03
        HSS( 8, 6, 4) =  0.20820D-01
        HSS( 8, 7, 4) = -0.28000D-03
        HSS( 8, 8, 4) =  0.41243D+00
        HSS( 9, 1, 4) = -0.39000D-03
        HSS( 9, 2, 4) = -0.39000D-03
        HSS( 9, 3, 4) = -0.57000D-03
        HSS( 9, 4, 4) = -0.62000D-03
        HSS( 9, 5, 4) =  0.12000D-03
        HSS( 9, 6, 4) =  0.43600D-02
        HSS( 9, 7, 4) =  0.43000D-03
        HSS( 9, 8, 4) =  0.43600D-02
        HSS( 9, 9, 4) =  0.36801D+00
        HSS(10, 1, 4) = -0.17300D-02
        HSS(10, 2, 4) =  0.20570D-01
        HSS(10, 3, 4) = -0.32200D-02
        HSS(10, 4, 4) =  0.86500D-01
        HSS(10, 5, 4) = -0.80000D-03
        HSS(10, 6, 4) = -0.41800D-02
        HSS(10, 7, 4) = -0.54000D-03
        HSS(10, 8, 4) =  0.21260D-01
        HSS(10, 9, 4) = -0.62000D-03
        HSS(10,10, 4) =  0.41668D+00
        HSS(11, 1, 4) = -0.39000D-03
        HSS(11, 2, 4) = -0.51000D-03
        HSS(11, 3, 4) =  0.41000D-03
        HSS(11, 4, 4) = -0.54000D-03
        HSS(11, 5, 4) =  0.90000D-04
        HSS(11, 6, 4) = -0.28000D-03
        HSS(11, 7, 4) =  0.20000D-03
        HSS(11, 8, 4) =  0.45000D-02
        HSS(11, 9, 4) =  0.43000D-03
        HSS(11,10, 4) =  0.44700D-02
        HSS(11,11, 4) =  0.36452D+00
        HSS(12, 1, 4) =  0.48000D-03
        HSS(12, 2, 4) =  0.50700D-02
        HSS(12, 3, 4) = -0.24900D-02
        HSS(12, 4, 4) = -0.80000D-03
        HSS(12, 5, 4) = -0.30000D-04
        HSS(12, 6, 4) = -0.36000D-03
        HSS(12, 7, 4) =  0.90000D-04
        HSS(12, 8, 4) = -0.91000D-03
        HSS(12, 9, 4) =  0.12000D-03
        HSS(12,10, 4) =  0.39900D-02
        HSS(12,11, 4) =  0.43000D-03
        HSS(12,12, 4) =  0.36905D+00
C Hessian elements for anchor point 4
        HBB( 1, 1, 4) =  0.80030D-01
        HBB( 2, 1, 4) = -0.40020D-01
        HBB( 2, 2, 4) =  0.92960D-01
        HBB( 3, 1, 4) = -0.40020D-01
        HBB( 3, 2, 4) = -0.52950D-01
        HBB( 3, 3, 4) =  0.92960D-01
        HBB( 4, 1, 4) = -0.33740D-01
        HBB( 4, 2, 4) =  0.17860D-01
        HBB( 4, 3, 4) =  0.15880D-01
        HBB( 4, 4, 4) =  0.77910D-01
        HBB( 5, 1, 4) =  0.19460D-01
        HBB( 5, 2, 4) = -0.10850D-01
        HBB( 5, 3, 4) = -0.86100D-02
        HBB( 5, 4, 4) = -0.38850D-01
        HBB( 5, 5, 4) =  0.77000D-01
        HBB( 6, 1, 4) =  0.14280D-01
        HBB( 6, 2, 4) = -0.70200D-02
        HBB( 6, 3, 4) = -0.72700D-02
        HBB( 6, 4, 4) = -0.39060D-01
        HBB( 6, 5, 4) = -0.38150D-01
        HBB( 6, 6, 4) =  0.77210D-01
        HBB( 7, 1, 4) = -0.89600D-02
        HBB( 7, 2, 4) = -0.20000D-04
        HBB( 7, 3, 4) =  0.89900D-02
        HBB( 7, 4, 4) = -0.32200D-01
        HBB( 7, 5, 4) =  0.13030D-01
        HBB( 7, 6, 4) =  0.19170D-01
        HBB( 7, 7, 4) =  0.79490D-01
        HBB( 8, 1, 4) = -0.12600D-02
        HBB( 8, 2, 4) =  0.34000D-02
        HBB( 8, 3, 4) = -0.21400D-02
        HBB( 8, 4, 4) =  0.18730D-01
        HBB( 8, 5, 4) = -0.69300D-02
        HBB( 8, 6, 4) = -0.11800D-01
        HBB( 8, 7, 4) = -0.39640D-01
        HBB( 8, 8, 4) =  0.79360D-01
        HBB( 9, 1, 4) =  0.10220D-01
        HBB( 9, 2, 4) = -0.33800D-02
        HBB( 9, 3, 4) = -0.68500D-02
        HBB( 9, 4, 4) =  0.13460D-01
        HBB( 9, 5, 4) = -0.60900D-02
        HBB( 9, 6, 4) = -0.73700D-02
        HBB( 9, 7, 4) = -0.39850D-01
        HBB( 9, 8, 4) = -0.39720D-01
        HBB( 9, 9, 4) =  0.79570D-01
        HBB(10, 1, 4) =  0.53800D-02
        HBB(10, 2, 4) = -0.26900D-02
        HBB(10, 3, 4) = -0.26900D-02
        HBB(10, 4, 4) = -0.85500D-02
        HBB(10, 5, 4) =  0.98800D-02
        HBB(10, 6, 4) = -0.13400D-02
        HBB(10, 7, 4) = -0.34400D-01
        HBB(10, 8, 4) =  0.14090D-01
        HBB(10, 9, 4) =  0.20310D-01
        HBB(10,10, 4) =  0.80510D-01
        HBB(11, 1, 4) = -0.26900D-02
        HBB(11, 2, 4) =  0.16300D-02
        HBB(11, 3, 4) =  0.10600D-02
        HBB(11, 4, 4) = -0.14300D-02
        HBB(11, 5, 4) = -0.33500D-02
        HBB(11, 6, 4) =  0.47900D-02
        HBB(11, 7, 4) =  0.20180D-01
        HBB(11, 8, 4) = -0.75700D-02
        HBB(11, 9, 4) = -0.12610D-01
        HBB(11,10, 4) = -0.40260D-01
        HBB(11,11, 4) =  0.79330D-01
        HBB(12, 1, 4) = -0.26900D-02
        HBB(12, 2, 4) =  0.10600D-02
        HBB(12, 3, 4) =  0.16300D-02
        HBB(12, 4, 4) =  0.99800D-02
        HBB(12, 5, 4) = -0.65300D-02
        HBB(12, 6, 4) = -0.34500D-02
        HBB(12, 7, 4) =  0.14220D-01
        HBB(12, 8, 4) = -0.65200D-02
        HBB(12, 9, 4) = -0.77100D-02
        HBB(12,10, 4) = -0.40260D-01
        HBB(12,11, 4) = -0.39080D-01
        HBB(12,12, 4) =  0.79330D-01
        HBB(13, 1, 4) = -0.89600D-02
        HBB(13, 2, 4) =  0.89900D-02
        HBB(13, 3, 4) = -0.20000D-04
        HBB(13, 4, 4) =  0.49400D-02
        HBB(13, 5, 4) = -0.20900D-02
        HBB(13, 6, 4) = -0.28500D-02
        HBB(13, 7, 4) = -0.88700D-02
        HBB(13, 8, 4) =  0.10300D-01
        HBB(13, 9, 4) = -0.14300D-02
        HBB(13,10, 4) = -0.34400D-01
        HBB(13,11, 4) =  0.14220D-01
        HBB(13,12, 4) =  0.20180D-01
        HBB(13,13, 4) =  0.79490D-01
        HBB(14, 1, 4) =  0.10220D-01
        HBB(14, 2, 4) = -0.68500D-02
        HBB(14, 3, 4) = -0.33800D-02
        HBB(14, 4, 4) = -0.27200D-02
        HBB(14, 5, 4) =  0.10700D-02
        HBB(14, 6, 4) =  0.16500D-02
        HBB(14, 7, 4) = -0.14300D-02
        HBB(14, 8, 4) = -0.36200D-02
        HBB(14, 9, 4) =  0.50500D-02
        HBB(14,10, 4) =  0.20310D-01
        HBB(14,11, 4) = -0.77100D-02
        HBB(14,12, 4) = -0.12610D-01
        HBB(14,13, 4) = -0.39850D-01
        HBB(14,14, 4) =  0.79570D-01
        HBB(15, 1, 4) = -0.12600D-02
        HBB(15, 2, 4) = -0.21400D-02
        HBB(15, 3, 4) =  0.34000D-02
        HBB(15, 4, 4) = -0.22200D-02
        HBB(15, 5, 4) =  0.10200D-02
        HBB(15, 6, 4) =  0.12000D-02
        HBB(15, 7, 4) =  0.10300D-01
        HBB(15, 8, 4) = -0.66700D-02
        HBB(15, 9, 4) = -0.36200D-02
        HBB(15,10, 4) =  0.14090D-01
        HBB(15,11, 4) = -0.65200D-02
        HBB(15,12, 4) = -0.75700D-02
        HBB(15,13, 4) = -0.39640D-01
        HBB(15,14, 4) = -0.39720D-01
        HBB(15,15, 4) =  0.79360D-01
        HBB(16, 1, 4) = -0.33740D-01
        HBB(16, 2, 4) =  0.15880D-01
        HBB(16, 3, 4) =  0.17860D-01
        HBB(16, 4, 4) = -0.83600D-02
        HBB(16, 5, 4) = -0.14300D-02
        HBB(16, 6, 4) =  0.97900D-02
        HBB(16, 7, 4) =  0.49400D-02
        HBB(16, 8, 4) = -0.22200D-02
        HBB(16, 9, 4) = -0.27200D-02
        HBB(16,10, 4) = -0.85500D-02
        HBB(16,11, 4) =  0.99800D-02
        HBB(16,12, 4) = -0.14300D-02
        HBB(16,13, 4) = -0.32200D-01
        HBB(16,14, 4) =  0.13460D-01
        HBB(16,15, 4) =  0.18730D-01
        HBB(16,16, 4) =  0.77910D-01
        HBB(17, 1, 4) =  0.19460D-01
        HBB(17, 2, 4) = -0.86100D-02
        HBB(17, 3, 4) = -0.10850D-01
        HBB(17, 4, 4) = -0.14300D-02
        HBB(17, 5, 4) =  0.48200D-02
        HBB(17, 6, 4) = -0.33900D-02
        HBB(17, 7, 4) = -0.20900D-02
        HBB(17, 8, 4) =  0.10200D-02
        HBB(17, 9, 4) =  0.10700D-02
        HBB(17,10, 4) =  0.98800D-02
        HBB(17,11, 4) = -0.65300D-02
        HBB(17,12, 4) = -0.33500D-02
        HBB(17,13, 4) =  0.13030D-01
        HBB(17,14, 4) = -0.60900D-02
        HBB(17,15, 4) = -0.69300D-02
        HBB(17,16, 4) = -0.38850D-01
        HBB(17,17, 4) =  0.77000D-01
        HBB(18, 1, 4) =  0.14280D-01
        HBB(18, 2, 4) = -0.72700D-02
        HBB(18, 3, 4) = -0.70200D-02
        HBB(18, 4, 4) =  0.97900D-02
        HBB(18, 5, 4) = -0.33900D-02
        HBB(18, 6, 4) = -0.63900D-02
        HBB(18, 7, 4) = -0.28500D-02
        HBB(18, 8, 4) =  0.12000D-02
        HBB(18, 9, 4) =  0.16500D-02
        HBB(18,10, 4) = -0.13400D-02
        HBB(18,11, 4) = -0.34500D-02
        HBB(18,12, 4) =  0.47900D-02
        HBB(18,13, 4) =  0.19170D-01
        HBB(18,14, 4) = -0.73700D-02
        HBB(18,15, 4) = -0.11800D-01
        HBB(18,16, 4) = -0.39060D-01
        HBB(18,17, 4) = -0.38150D-01
        HBB(18,18, 4) =  0.77210D-01
C Hessian elements for anchor point 4
        HTT( 1, 1, 4) =  0.84400D-02
        HTT( 2, 1, 4) =  0.31800D-02
        HTT( 2, 2, 4) =  0.76400D-02
        HTT( 3, 1, 4) =  0.17500D-02
        HTT( 3, 2, 4) = -0.29500D-02
        HTT( 3, 3, 4) =  0.10070D-01
        HTT( 4, 1, 4) = -0.35100D-02
        HTT( 4, 2, 4) =  0.15200D-02
        HTT( 4, 3, 4) =  0.53700D-02
        HTT( 4, 4, 4) =  0.10400D-01
        HTT( 5, 1, 4) = -0.54900D-02
        HTT( 5, 2, 4) = -0.37900D-02
        HTT( 5, 3, 4) =  0.12100D-02
        HTT( 5, 4, 4) =  0.29000D-02
        HTT( 5, 5, 4) =  0.84400D-02
        HTT( 6, 1, 4) = -0.37900D-02
        HTT( 6, 2, 4) = -0.24900D-02
        HTT( 6, 3, 4) =  0.23300D-02
        HTT( 6, 4, 4) =  0.36400D-02
        HTT( 6, 5, 4) =  0.31800D-02
        HTT( 6, 6, 4) =  0.76400D-02
        HTT( 7, 1, 4) =  0.12100D-02
        HTT( 7, 2, 4) =  0.23300D-02
        HTT( 7, 3, 4) = -0.71100D-02
        HTT( 7, 4, 4) = -0.59900D-02
        HTT( 7, 5, 4) =  0.17500D-02
        HTT( 7, 6, 4) = -0.29500D-02
        HTT( 7, 7, 4) =  0.10070D-01
        HTT( 8, 1, 4) =  0.29000D-02
        HTT( 8, 2, 4) =  0.36400D-02
        HTT( 8, 3, 4) = -0.59900D-02
        HTT( 8, 4, 4) = -0.52500D-02
        HTT( 8, 5, 4) = -0.35100D-02
        HTT( 8, 6, 4) =  0.15200D-02
        HTT( 8, 7, 4) =  0.53700D-02
        HTT( 8, 8, 4) =  0.10400D-01
        HTT( 9, 1, 4) = -0.60400D-02
        HTT( 9, 2, 4) = -0.50000D-03
        HTT( 9, 3, 4) = -0.37000D-02
        HTT( 9, 4, 4) =  0.18400D-02
        HTT( 9, 5, 4) =  0.10000D-04
        HTT( 9, 6, 4) =  0.17500D-02
        HTT( 9, 7, 4) = -0.23300D-02
        HTT( 9, 8, 4) = -0.59000D-03
        HTT( 9, 9, 4) =  0.94900D-02
        HTT(10, 1, 4) = -0.38600D-02
        HTT(10, 2, 4) =  0.97000D-03
        HTT(10, 3, 4) = -0.21900D-02
        HTT(10, 4, 4) =  0.26400D-02
        HTT(10, 5, 4) =  0.17900D-02
        HTT(10, 6, 4) =  0.13600D-02
        HTT(10, 7, 4) =  0.12000D-03
        HTT(10, 8, 4) = -0.31000D-03
        HTT(10, 9, 4) =  0.31000D-02
        HTT(10,10, 4) =  0.84100D-02
        HTT(11, 1, 4) = -0.65000D-03
        HTT(11, 2, 4) = -0.50600D-02
        HTT(11, 3, 4) =  0.11100D-02
        HTT(11, 4, 4) = -0.33000D-02
        HTT(11, 5, 4) = -0.17200D-02
        HTT(11, 6, 4) =  0.41000D-03
        HTT(11, 7, 4) = -0.34800D-02
        HTT(11, 8, 4) = -0.13500D-02
        HTT(11, 9, 4) =  0.38300D-02
        HTT(11,10, 4) = -0.18300D-02
        HTT(11,11, 4) =  0.83400D-02
        HTT(12, 1, 4) =  0.15300D-02
        HTT(12, 2, 4) = -0.35900D-02
        HTT(12, 3, 4) =  0.26200D-02
        HTT(12, 4, 4) = -0.25000D-02
        HTT(12, 5, 4) =  0.60000D-04
        HTT(12, 6, 4) =  0.20000D-04
        HTT(12, 7, 4) = -0.10300D-02
        HTT(12, 8, 4) = -0.10700D-02
        HTT(12, 9, 4) = -0.25600D-02
        HTT(12,10, 4) =  0.34800D-02
        HTT(12,11, 4) =  0.26800D-02
        HTT(12,12, 4) =  0.87100D-02
        HTT(13, 1, 4) =  0.74000D-03
        HTT(13, 2, 4) = -0.14900D-02
        HTT(13, 3, 4) =  0.26300D-02
        HTT(13, 4, 4) =  0.40000D-03
        HTT(13, 5, 4) =  0.23400D-02
        HTT(13, 6, 4) =  0.85000D-03
        HTT(13, 7, 4) =  0.45000D-03
        HTT(13, 8, 4) = -0.10400D-02
        HTT(13, 9, 4) = -0.68300D-02
        HTT(13,10, 4) = -0.31000D-03
        HTT(13,11, 4) = -0.45500D-02
        HTT(13,12, 4) =  0.19800D-02
        HTT(13,13, 4) =  0.96800D-02
        HTT(14, 1, 4) =  0.19800D-02
        HTT(14, 2, 4) =  0.49000D-03
        HTT(14, 3, 4) =  0.15200D-02
        HTT(14, 4, 4) =  0.30000D-04
        HTT(14, 5, 4) =  0.11000D-02
        HTT(14, 6, 4) = -0.11300D-02
        HTT(14, 7, 4) =  0.15600D-02
        HTT(14, 8, 4) = -0.67000D-03
        HTT(14, 9, 4) = -0.45100D-02
        HTT(14,10, 4) =  0.12300D-02
        HTT(14,11, 4) = -0.29800D-02
        HTT(14,12, 4) =  0.27600D-02
        HTT(14,13, 4) =  0.38800D-02
        HTT(14,14, 4) =  0.86200D-02
        HTT(15, 1, 4) = -0.14700D-02
        HTT(15, 2, 4) = -0.29800D-02
        HTT(15, 3, 4) =  0.11000D-02
        HTT(15, 4, 4) = -0.41000D-03
        HTT(15, 5, 4) =  0.54000D-03
        HTT(15, 6, 4) =  0.12400D-02
        HTT(15, 7, 4) = -0.20300D-02
        HTT(15, 8, 4) = -0.13200D-02
        HTT(15, 9, 4) = -0.36000D-03
        HTT(15,10, 4) = -0.56900D-02
        HTT(15,11, 4) =  0.11900D-02
        HTT(15,12, 4) = -0.41400D-02
        HTT(15,13, 4) =  0.30600D-02
        HTT(15,14, 4) = -0.19400D-02
        HTT(15,15, 4) =  0.84700D-02
        HTT(16, 1, 4) = -0.22000D-03
        HTT(16, 2, 4) = -0.10000D-02
        HTT(16, 3, 4) = -0.10000D-04
        HTT(16, 4, 4) = -0.78000D-03
        HTT(16, 5, 4) = -0.71000D-03
        HTT(16, 6, 4) = -0.74000D-03
        HTT(16, 7, 4) = -0.92000D-03
        HTT(16, 8, 4) = -0.95000D-03
        HTT(16, 9, 4) =  0.19700D-02
        HTT(16,10, 4) = -0.41500D-02
        HTT(16,11, 4) =  0.27600D-02
        HTT(16,12, 4) = -0.33600D-02
        HTT(16,13, 4) = -0.27400D-02
        HTT(16,14, 4) =  0.28000D-02
        HTT(16,15, 4) =  0.34700D-02
        HTT(16,16, 4) =  0.90000D-02
        HTT(17, 1, 4) =  0.23400D-02
        HTT(17, 2, 4) =  0.85000D-03
        HTT(17, 3, 4) =  0.45000D-03
        HTT(17, 4, 4) = -0.10400D-02
        HTT(17, 5, 4) =  0.74000D-03
        HTT(17, 6, 4) = -0.14900D-02
        HTT(17, 7, 4) =  0.26300D-02
        HTT(17, 8, 4) =  0.40000D-03
        HTT(17, 9, 4) =  0.56000D-03
        HTT(17,10, 4) = -0.18400D-02
        HTT(17,11, 4) =  0.20800D-02
        HTT(17,12, 4) = -0.32000D-03
        HTT(17,13, 4) = -0.64700D-02
        HTT(17,14, 4) = -0.67000D-03
        HTT(17,15, 4) = -0.40300D-02
        HTT(17,16, 4) =  0.17700D-02
        HTT(17,17, 4) =  0.96800D-02
        HTT(18, 1, 4) =  0.54000D-03
        HTT(18, 2, 4) =  0.12400D-02
        HTT(18, 3, 4) = -0.20300D-02
        HTT(18, 4, 4) = -0.13200D-02
        HTT(18, 5, 4) = -0.14700D-02
        HTT(18, 6, 4) = -0.29800D-02
        HTT(18, 7, 4) =  0.11000D-02
        HTT(18, 8, 4) = -0.41000D-03
        HTT(18, 9, 4) =  0.22500D-02
        HTT(18,10, 4) =  0.51000D-03
        HTT(18,11, 4) =  0.15300D-02
        HTT(18,12, 4) = -0.22000D-03
        HTT(18,13, 4) = -0.40300D-02
        HTT(18,14, 4) =  0.97000D-03
        HTT(18,15, 4) = -0.22600D-02
        HTT(18,16, 4) =  0.27400D-02
        HTT(18,17, 4) =  0.30600D-02
        HTT(18,18, 4) =  0.84700D-02
        HTT(19, 1, 4) =  0.11000D-02
        HTT(19, 2, 4) = -0.11300D-02
        HTT(19, 3, 4) =  0.15600D-02
        HTT(19, 4, 4) = -0.67000D-03
        HTT(19, 5, 4) =  0.19800D-02
        HTT(19, 6, 4) =  0.49000D-03
        HTT(19, 7, 4) =  0.15200D-02
        HTT(19, 8, 4) =  0.30000D-04
        HTT(19, 9, 4) = -0.17700D-02
        HTT(19,10, 4) = -0.33800D-02
        HTT(19,11, 4) =  0.51000D-03
        HTT(19,12, 4) = -0.11100D-02
        HTT(19,13, 4) = -0.67000D-03
        HTT(19,14, 4) = -0.54000D-02
        HTT(19,15, 4) =  0.97000D-03
        HTT(19,16, 4) = -0.37700D-02
        HTT(19,17, 4) =  0.38800D-02
        HTT(19,18, 4) = -0.19400D-02
        HTT(19,19, 4) =  0.86200D-02
        HTT(20, 1, 4) = -0.71000D-03
        HTT(20, 2, 4) = -0.74000D-03
        HTT(20, 3, 4) = -0.92000D-03
        HTT(20, 4, 4) = -0.95000D-03
        HTT(20, 5, 4) = -0.22000D-03
        HTT(20, 6, 4) = -0.10000D-02
        HTT(20, 7, 4) = -0.10000D-04
        HTT(20, 8, 4) = -0.78000D-03
        HTT(20, 9, 4) = -0.70000D-04
        HTT(20,10, 4) = -0.10300D-02
        HTT(20,11, 4) = -0.40000D-04
        HTT(20,12, 4) = -0.10000D-02
        HTT(20,13, 4) =  0.17700D-02
        HTT(20,14, 4) = -0.37700D-02
        HTT(20,15, 4) =  0.27400D-02
        HTT(20,16, 4) = -0.28000D-02
        HTT(20,17, 4) = -0.27400D-02
        HTT(20,18, 4) =  0.34700D-02
        HTT(20,19, 4) =  0.28000D-02
        HTT(20,20, 4) =  0.90000D-02
        HTT(21, 1, 4) =  0.10000D-04
        HTT(21, 2, 4) =  0.17500D-02
        HTT(21, 3, 4) = -0.23300D-02
        HTT(21, 4, 4) = -0.59000D-03
        HTT(21, 5, 4) = -0.60400D-02
        HTT(21, 6, 4) = -0.50000D-03
        HTT(21, 7, 4) = -0.37000D-02
        HTT(21, 8, 4) =  0.18400D-02
        HTT(21, 9, 4) =  0.27700D-02
        HTT(21,10, 4) =  0.11000D-02
        HTT(21,11, 4) =  0.10000D-02
        HTT(21,12, 4) = -0.68000D-03
        HTT(21,13, 4) =  0.56000D-03
        HTT(21,14, 4) = -0.17700D-02
        HTT(21,15, 4) =  0.22500D-02
        HTT(21,16, 4) = -0.70000D-04
        HTT(21,17, 4) = -0.68300D-02
        HTT(21,18, 4) = -0.36000D-03
        HTT(21,19, 4) = -0.45100D-02
        HTT(21,20, 4) =  0.19700D-02
        HTT(21,21, 4) =  0.94900D-02
        HTT(22, 1, 4) = -0.17200D-02
        HTT(22, 2, 4) =  0.41000D-03
        HTT(22, 3, 4) = -0.34800D-02
        HTT(22, 4, 4) = -0.13500D-02
        HTT(22, 5, 4) = -0.65000D-03
        HTT(22, 6, 4) = -0.50600D-02
        HTT(22, 7, 4) =  0.11100D-02
        HTT(22, 8, 4) = -0.33000D-02
        HTT(22, 9, 4) =  0.10000D-02
        HTT(22,10, 4) =  0.15400D-02
        HTT(22,11, 4) = -0.11800D-02
        HTT(22,12, 4) = -0.64000D-03
        HTT(22,13, 4) =  0.20800D-02
        HTT(22,14, 4) =  0.51000D-03
        HTT(22,15, 4) =  0.15300D-02
        HTT(22,16, 4) = -0.40000D-04
        HTT(22,17, 4) = -0.45500D-02
        HTT(22,18, 4) =  0.11900D-02
        HTT(22,19, 4) = -0.29800D-02
        HTT(22,20, 4) =  0.27600D-02
        HTT(22,21, 4) =  0.38300D-02
        HTT(22,22, 4) =  0.83400D-02
        HTT(23, 1, 4) =  0.17900D-02
        HTT(23, 2, 4) =  0.13600D-02
        HTT(23, 3, 4) =  0.12000D-03
        HTT(23, 4, 4) = -0.31000D-03
        HTT(23, 5, 4) = -0.38600D-02
        HTT(23, 6, 4) =  0.97000D-03
        HTT(23, 7, 4) = -0.21900D-02
        HTT(23, 8, 4) =  0.26400D-02
        HTT(23, 9, 4) =  0.11000D-02
        HTT(23,10, 4) = -0.12200D-02
        HTT(23,11, 4) =  0.15400D-02
        HTT(23,12, 4) = -0.78000D-03
        HTT(23,13, 4) = -0.18400D-02
        HTT(23,14, 4) = -0.33800D-02
        HTT(23,15, 4) =  0.51000D-03
        HTT(23,16, 4) = -0.10300D-02
        HTT(23,17, 4) = -0.31000D-03
        HTT(23,18, 4) = -0.56900D-02
        HTT(23,19, 4) =  0.12300D-02
        HTT(23,20, 4) = -0.41500D-02
        HTT(23,21, 4) =  0.31000D-02
        HTT(23,22, 4) = -0.18300D-02
        HTT(23,23, 4) =  0.84100D-02
        HTT(24, 1, 4) =  0.60000D-04
        HTT(24, 2, 4) =  0.20000D-04
        HTT(24, 3, 4) = -0.10300D-02
        HTT(24, 4, 4) = -0.10700D-02
        HTT(24, 5, 4) =  0.15300D-02
        HTT(24, 6, 4) = -0.35900D-02
        HTT(24, 7, 4) =  0.26200D-02
        HTT(24, 8, 4) = -0.25000D-02
        HTT(24, 9, 4) = -0.68000D-03
        HTT(24,10, 4) = -0.78000D-03
        HTT(24,11, 4) = -0.64000D-03
        HTT(24,12, 4) = -0.74000D-03
        HTT(24,13, 4) = -0.32000D-03
        HTT(24,14, 4) = -0.11100D-02
        HTT(24,15, 4) = -0.22000D-03
        HTT(24,16, 4) = -0.10000D-02
        HTT(24,17, 4) =  0.19800D-02
        HTT(24,18, 4) = -0.41400D-02
        HTT(24,19, 4) =  0.27600D-02
        HTT(24,20, 4) = -0.33600D-02
        HTT(24,21, 4) = -0.25600D-02
        HTT(24,22, 4) =  0.26800D-02
        HTT(24,23, 4) =  0.34800D-02
        HTT(24,24, 4) =  0.87100D-02
C Hessian elements for anchor point 4
        HBS( 1, 1, 4) = -0.31760D-01
        HBS( 1, 2, 4) = -0.31760D-01
        HBS( 1, 3, 4) = -0.37700D-01
        HBS( 1, 4, 4) = -0.11100D-01
        HBS( 1, 5, 4) =  0.40600D-02
        HBS( 1, 6, 4) =  0.43670D-01
        HBS( 1, 7, 4) = -0.35000D-03
        HBS( 1, 8, 4) =  0.43670D-01
        HBS( 1, 9, 4) = -0.10000D-04
        HBS( 1,10, 4) = -0.11100D-01
        HBS( 1,11, 4) = -0.35000D-03
        HBS( 1,12, 4) =  0.40600D-02
        HBS( 2, 1, 4) =  0.39800D-01
        HBS( 2, 2, 4) = -0.80400D-02
        HBS( 2, 3, 4) =  0.18850D-01
        HBS( 2, 4, 4) = -0.12000D-03
        HBS( 2, 5, 4) = -0.21800D-02
        HBS( 2, 6, 4) = -0.18210D-01
        HBS( 2, 7, 4) =  0.75000D-03
        HBS( 2, 8, 4) = -0.25460D-01
        HBS( 2, 9, 4) =  0.00000D+00
        HBS( 2,10, 4) =  0.11220D-01
        HBS( 2,11, 4) = -0.39000D-03
        HBS( 2,12, 4) = -0.18800D-02
        HBS( 3, 1, 4) = -0.80400D-02
        HBS( 3, 2, 4) =  0.39800D-01
        HBS( 3, 3, 4) =  0.18850D-01
        HBS( 3, 4, 4) =  0.11220D-01
        HBS( 3, 5, 4) = -0.18800D-02
        HBS( 3, 6, 4) = -0.25460D-01
        HBS( 3, 7, 4) = -0.39000D-03
        HBS( 3, 8, 4) = -0.18210D-01
        HBS( 3, 9, 4) =  0.00000D+00
        HBS( 3,10, 4) = -0.12000D-03
        HBS( 3,11, 4) =  0.75000D-03
        HBS( 3,12, 4) = -0.21800D-02
        HBS( 4, 1, 4) = -0.26690D-01
        HBS( 4, 2, 4) = -0.99100D-02
        HBS( 4, 3, 4) =  0.10200D-01
        HBS( 4, 4, 4) = -0.32600D-01
        HBS( 4, 5, 4) = -0.82100D-02
        HBS( 4, 6, 4) = -0.12130D-01
        HBS( 4, 7, 4) =  0.47800D-02
        HBS( 4, 8, 4) =  0.46360D-01
        HBS( 4, 9, 4) = -0.61000D-03
        HBS( 4,10, 4) =  0.44260D-01
        HBS( 4,11, 4) =  0.10000D-03
        HBS( 4,12, 4) = -0.43000D-03
        HBS( 5, 1, 4) =  0.27720D-01
        HBS( 5, 2, 4) =  0.58200D-02
        HBS( 5, 3, 4) = -0.92100D-02
        HBS( 5, 4, 4) =  0.33000D-02
        HBS( 5, 5, 4) =  0.59500D-02
        HBS( 5, 6, 4) =  0.56400D-02
        HBS( 5, 7, 4) = -0.18300D-02
        HBS( 5, 8, 4) = -0.20790D-01
        HBS( 5, 9, 4) = -0.31000D-03
        HBS( 5,10, 4) = -0.24570D-01
        HBS( 5,11, 4) = -0.60000D-04
        HBS( 5,12, 4) =  0.64000D-03
        HBS( 6, 1, 4) = -0.10300D-02
        HBS( 6, 2, 4) =  0.40900D-02
        HBS( 6, 3, 4) = -0.99000D-03
        HBS( 6, 4, 4) =  0.29300D-01
        HBS( 6, 5, 4) =  0.22600D-02
        HBS( 6, 6, 4) =  0.65000D-02
        HBS( 6, 7, 4) = -0.29500D-02
        HBS( 6, 8, 4) = -0.25570D-01
        HBS( 6, 9, 4) =  0.92000D-03
        HBS( 6,10, 4) = -0.19690D-01
        HBS( 6,11, 4) = -0.40000D-04
        HBS( 6,12, 4) = -0.21000D-03
        HBS( 7, 1, 4) = -0.14540D-01
        HBS( 7, 2, 4) =  0.40000D-01
        HBS( 7, 3, 4) =  0.98300D-02
        HBS( 7, 4, 4) = -0.32290D-01
        HBS( 7, 5, 4) =  0.50700D-02
        HBS( 7, 6, 4) = -0.32910D-01
        HBS( 7, 7, 4) = -0.88600D-02
        HBS( 7, 8, 4) = -0.12410D-01
        HBS( 7, 9, 4) =  0.48800D-02
        HBS( 7,10, 4) =  0.45460D-01
        HBS( 7,11, 4) = -0.87000D-03
        HBS( 7,12, 4) =  0.14000D-03
        HBS( 8, 1, 4) =  0.89100D-02
        HBS( 8, 2, 4) = -0.22340D-01
        HBS( 8, 3, 4) = -0.58300D-02
        HBS( 8, 4, 4) =  0.30270D-01
        HBS( 8, 5, 4) = -0.31400D-02
        HBS( 8, 6, 4) =  0.36400D-02
        HBS( 8, 7, 4) =  0.47300D-02
        HBS( 8, 8, 4) =  0.54700D-02
        HBS( 8, 9, 4) = -0.19600D-02
        HBS( 8,10, 4) = -0.20340D-01
        HBS( 8,11, 4) = -0.23000D-03
        HBS( 8,12, 4) = -0.12000D-03
        HBS( 9, 1, 4) =  0.56200D-02
        HBS( 9, 2, 4) = -0.17660D-01
        HBS( 9, 3, 4) = -0.39900D-02
        HBS( 9, 4, 4) =  0.20200D-02
        HBS( 9, 5, 4) = -0.19200D-02
        HBS( 9, 6, 4) =  0.29280D-01
        HBS( 9, 7, 4) =  0.41300D-02
        HBS( 9, 8, 4) =  0.69400D-02
        HBS( 9, 9, 4) = -0.29200D-02
        HBS( 9,10, 4) = -0.25120D-01
        HBS( 9,11, 4) =  0.11000D-02
        HBS( 9,12, 4) = -0.30000D-04
        HBS(10, 1, 4) =  0.42900D-01
        HBS(10, 2, 4) =  0.42900D-01
        HBS(10, 3, 4) = -0.23400D-02
        HBS(10, 4, 4) = -0.13730D-01
        HBS(10, 5, 4) = -0.64000D-03
        HBS(10, 6, 4) = -0.32580D-01
        HBS(10, 7, 4) =  0.52100D-02
        HBS(10, 8, 4) = -0.32580D-01
        HBS(10, 9, 4) = -0.85400D-02
        HBS(10,10, 4) = -0.13730D-01
        HBS(10,11, 4) =  0.52100D-02
        HBS(10,12, 4) = -0.64000D-03
        HBS(11, 1, 4) = -0.23760D-01
        HBS(11, 2, 4) = -0.19130D-01
        HBS(11, 3, 4) =  0.11700D-02
        HBS(11, 4, 4) =  0.73400D-02
        HBS(11, 5, 4) =  0.91000D-03
        HBS(11, 6, 4) =  0.30250D-01
        HBS(11, 7, 4) = -0.31400D-02
        HBS(11, 8, 4) =  0.23200D-02
        HBS(11, 9, 4) =  0.42700D-02
        HBS(11,10, 4) =  0.64000D-02
        HBS(11,11, 4) = -0.20700D-02
        HBS(11,12, 4) = -0.27000D-03
        HBS(12, 1, 4) = -0.19130D-01
        HBS(12, 2, 4) = -0.23760D-01
        HBS(12, 3, 4) =  0.11700D-02
        HBS(12, 4, 4) =  0.64000D-02
        HBS(12, 5, 4) = -0.27000D-03
        HBS(12, 6, 4) =  0.23200D-02
        HBS(12, 7, 4) = -0.20700D-02
        HBS(12, 8, 4) =  0.30250D-01
        HBS(12, 9, 4) =  0.42700D-02
        HBS(12,10, 4) =  0.73400D-02
        HBS(12,11, 4) = -0.31400D-02
        HBS(12,12, 4) =  0.91000D-03
        HBS(13, 1, 4) =  0.40000D-01
        HBS(13, 2, 4) = -0.14540D-01
        HBS(13, 3, 4) =  0.98300D-02
        HBS(13, 4, 4) =  0.45460D-01
        HBS(13, 5, 4) =  0.14000D-03
        HBS(13, 6, 4) = -0.12410D-01
        HBS(13, 7, 4) = -0.87000D-03
        HBS(13, 8, 4) = -0.32910D-01
        HBS(13, 9, 4) =  0.48800D-02
        HBS(13,10, 4) = -0.32290D-01
        HBS(13,11, 4) = -0.88600D-02
        HBS(13,12, 4) =  0.50700D-02
        HBS(14, 1, 4) = -0.17660D-01
        HBS(14, 2, 4) =  0.56200D-02
        HBS(14, 3, 4) = -0.39900D-02
        HBS(14, 4, 4) = -0.25120D-01
        HBS(14, 5, 4) = -0.30000D-04
        HBS(14, 6, 4) =  0.69400D-02
        HBS(14, 7, 4) =  0.11000D-02
        HBS(14, 8, 4) =  0.29280D-01
        HBS(14, 9, 4) = -0.29200D-02
        HBS(14,10, 4) =  0.20200D-02
        HBS(14,11, 4) =  0.41300D-02
        HBS(14,12, 4) = -0.19200D-02
        HBS(15, 1, 4) = -0.22340D-01
        HBS(15, 2, 4) =  0.89100D-02
        HBS(15, 3, 4) = -0.58300D-02
        HBS(15, 4, 4) = -0.20340D-01
        HBS(15, 5, 4) = -0.12000D-03
        HBS(15, 6, 4) =  0.54700D-02
        HBS(15, 7, 4) = -0.23000D-03
        HBS(15, 8, 4) =  0.36400D-02
        HBS(15, 9, 4) = -0.19600D-02
        HBS(15,10, 4) =  0.30270D-01
        HBS(15,11, 4) =  0.47300D-02
        HBS(15,12, 4) = -0.31400D-02
        HBS(16, 1, 4) = -0.99100D-02
        HBS(16, 2, 4) = -0.26690D-01
        HBS(16, 3, 4) =  0.10200D-01
        HBS(16, 4, 4) =  0.44260D-01
        HBS(16, 5, 4) = -0.43000D-03
        HBS(16, 6, 4) =  0.46360D-01
        HBS(16, 7, 4) =  0.10000D-03
        HBS(16, 8, 4) = -0.12130D-01
        HBS(16, 9, 4) = -0.61000D-03
        HBS(16,10, 4) = -0.32600D-01
        HBS(16,11, 4) =  0.47800D-02
        HBS(16,12, 4) = -0.82100D-02
        HBS(17, 1, 4) =  0.58200D-02
        HBS(17, 2, 4) =  0.27720D-01
        HBS(17, 3, 4) = -0.92100D-02
        HBS(17, 4, 4) = -0.24570D-01
        HBS(17, 5, 4) =  0.64000D-03
        HBS(17, 6, 4) = -0.20790D-01
        HBS(17, 7, 4) = -0.60000D-04
        HBS(17, 8, 4) =  0.56400D-02
        HBS(17, 9, 4) = -0.31000D-03
        HBS(17,10, 4) =  0.33000D-02
        HBS(17,11, 4) = -0.18300D-02
        HBS(17,12, 4) =  0.59500D-02
        HBS(18, 1, 4) =  0.40900D-02
        HBS(18, 2, 4) = -0.10300D-02
        HBS(18, 3, 4) = -0.99000D-03
        HBS(18, 4, 4) = -0.19690D-01
        HBS(18, 5, 4) = -0.21000D-03
        HBS(18, 6, 4) = -0.25570D-01
        HBS(18, 7, 4) = -0.40000D-04
        HBS(18, 8, 4) =  0.65000D-02
        HBS(18, 9, 4) =  0.92000D-03
        HBS(18,10, 4) =  0.29300D-01
        HBS(18,11, 4) = -0.29500D-02
        HBS(18,12, 4) =  0.22600D-02
      else if(ii .eq. 22) then
C Hessian elements for anchor point 1
        HSS( 1, 1, 1) =  0.35006D+00
        HSS( 2, 1, 1) =  0.13620D-01
        HSS( 2, 2, 1) =  0.36874D+00
        HSS( 3, 1, 1) =  0.29130D-01
        HSS( 3, 2, 1) =  0.28430D-01
        HSS( 3, 3, 1) =  0.46705D+00
        HSS( 4, 1, 1) = -0.29900D-02
        HSS( 4, 2, 1) =  0.16080D-01
        HSS( 4, 3, 1) = -0.77700D-02
        HSS( 4, 4, 1) =  0.38651D+00
        HSS( 5, 1, 1) =  0.40700D-02
        HSS( 5, 2, 1) = -0.10000D-02
        HSS( 5, 3, 1) =  0.38700D-02
        HSS( 5, 4, 1) =  0.32400D-02
        HSS( 5, 5, 1) =  0.36515D+00
        HSS( 6, 1, 1) =  0.33920D-01
        HSS( 6, 2, 1) =  0.46890D-01
        HSS( 6, 3, 1) = -0.39100D-02
        HSS( 6, 4, 1) = -0.52200D-02
        HSS( 6, 5, 1) =  0.40000D-03
        HSS( 6, 6, 1) =  0.36166D+00
        HSS( 7, 1, 1) = -0.40000D-03
        HSS( 7, 2, 1) = -0.70000D-03
        HSS( 7, 3, 1) =  0.45000D-03
        HSS( 7, 4, 1) =  0.39500D-02
        HSS( 7, 5, 1) =  0.42000D-03
        HSS( 7, 6, 1) =  0.39100D-02
        HSS( 7, 7, 1) =  0.37071D+00
        HSS( 8, 1, 1) =  0.44240D-01
        HSS( 8, 2, 1) =  0.34840D-01
        HSS( 8, 3, 1) = -0.58500D-02
        HSS( 8, 4, 1) =  0.24280D-01
        HSS( 8, 5, 1) = -0.11400D-02
        HSS( 8, 6, 1) =  0.14590D-01
        HSS( 8, 7, 1) = -0.98000D-03
        HSS( 8, 8, 1) =  0.36092D+00
        HSS( 9, 1, 1) = -0.57000D-03
        HSS( 9, 2, 1) = -0.54000D-03
        HSS( 9, 3, 1) =  0.41000D-03
        HSS( 9, 4, 1) = -0.31000D-03
        HSS( 9, 5, 1) =  0.17000D-03
        HSS( 9, 6, 1) =  0.39200D-02
        HSS( 9, 7, 1) =  0.50000D-03
        HSS( 9, 8, 1) =  0.39400D-02
        HSS( 9, 9, 1) =  0.36389D+00
        HSS(10, 1, 1) =  0.19730D-01
        HSS(10, 2, 1) = -0.83300D-02
        HSS(10, 3, 1) = -0.46800D-02
        HSS(10, 4, 1) =  0.63810D-01
        HSS(10, 5, 1) = -0.33000D-03
        HSS(10, 6, 1) =  0.22020D-01
        HSS(10, 7, 1) = -0.10900D-02
        HSS(10, 8, 1) = -0.21000D-02
        HSS(10, 9, 1) = -0.19000D-03
        HSS(10,10, 1) =  0.38332D+00
        HSS(11, 1, 1) = -0.70000D-03
        HSS(11, 2, 1) = -0.43000D-03
        HSS(11, 3, 1) =  0.55000D-03
        HSS(11, 4, 1) = -0.11400D-02
        HSS(11, 5, 1) =  0.40000D-04
        HSS(11, 6, 1) = -0.97000D-03
        HSS(11, 7, 1) =  0.22000D-03
        HSS(11, 8, 1) =  0.38800D-02
        HSS(11, 9, 1) =  0.48000D-03
        HSS(11,10, 1) =  0.39000D-02
        HSS(11,11, 1) =  0.37178D+00
        HSS(12, 1, 1) = -0.14000D-02
        HSS(12, 2, 1) =  0.34100D-02
        HSS(12, 3, 1) =  0.37000D-03
        HSS(12, 4, 1) = -0.13000D-03
        HSS(12, 5, 1) =  0.24000D-03
        HSS(12, 6, 1) = -0.10500D-02
        HSS(12, 7, 1) =  0.20000D-04
        HSS(12, 8, 1) =  0.55000D-03
        HSS(12, 9, 1) =  0.12000D-03
        HSS(12,10, 1) =  0.24400D-02
        HSS(12,11, 1) =  0.39000D-03
        HSS(12,12, 1) =  0.37251D+00
C Hessian elements for anchor point 1
        HBB( 1, 1, 1) =  0.82220D-01
        HBB( 2, 1, 1) = -0.42040D-01
        HBB( 2, 2, 1) =  0.14128D+00
        HBB( 3, 1, 1) = -0.40180D-01
        HBB( 3, 2, 1) = -0.99240D-01
        HBB( 3, 3, 1) =  0.13942D+00
        HBB( 4, 1, 1) = -0.36210D-01
        HBB( 4, 2, 1) =  0.22250D-01
        HBB( 4, 3, 1) =  0.13960D-01
        HBB( 4, 4, 1) =  0.66190D-01
        HBB( 5, 1, 1) =  0.21220D-01
        HBB( 5, 2, 1) = -0.13320D-01
        HBB( 5, 3, 1) = -0.79000D-02
        HBB( 5, 4, 1) = -0.33260D-01
        HBB( 5, 5, 1) =  0.72740D-01
        HBB( 6, 1, 1) =  0.14990D-01
        HBB( 6, 2, 1) = -0.89300D-02
        HBB( 6, 3, 1) = -0.60600D-02
        HBB( 6, 4, 1) = -0.32930D-01
        HBB( 6, 5, 1) = -0.39480D-01
        HBB( 6, 6, 1) =  0.72410D-01
        HBB( 7, 1, 1) = -0.67600D-02
        HBB( 7, 2, 1) = -0.46000D-02
        HBB( 7, 3, 1) =  0.11370D-01
        HBB( 7, 4, 1) = -0.26080D-01
        HBB( 7, 5, 1) =  0.93900D-02
        HBB( 7, 6, 1) =  0.16690D-01
        HBB( 7, 7, 1) =  0.64820D-01
        HBB( 8, 1, 1) = -0.24700D-02
        HBB( 8, 2, 1) =  0.77600D-02
        HBB( 8, 3, 1) = -0.52900D-02
        HBB( 8, 4, 1) =  0.15900D-01
        HBB( 8, 5, 1) = -0.49900D-02
        HBB( 8, 6, 1) = -0.10910D-01
        HBB( 8, 7, 1) = -0.31790D-01
        HBB( 8, 8, 1) =  0.73440D-01
        HBB( 9, 1, 1) =  0.92300D-02
        HBB( 9, 2, 1) = -0.31500D-02
        HBB( 9, 3, 1) = -0.60800D-02
        HBB( 9, 4, 1) =  0.10190D-01
        HBB( 9, 5, 1) = -0.44000D-02
        HBB( 9, 6, 1) = -0.57900D-02
        HBB( 9, 7, 1) = -0.33030D-01
        HBB( 9, 8, 1) = -0.41650D-01
        HBB( 9, 9, 1) =  0.74680D-01
        HBB(10, 1, 1) =  0.38400D-02
        HBB(10, 2, 1) = -0.14400D-02
        HBB(10, 3, 1) = -0.24000D-02
        HBB(10, 4, 1) = -0.57300D-02
        HBB(10, 5, 1) =  0.88600D-02
        HBB(10, 6, 1) = -0.31300D-02
        HBB(10, 7, 1) = -0.33500D-01
        HBB(10, 8, 1) =  0.13490D-01
        HBB(10, 9, 1) =  0.20010D-01
        HBB(10,10, 1) =  0.75070D-01
        HBB(11, 1, 1) = -0.18400D-02
        HBB(11, 2, 1) =  0.11400D-02
        HBB(11, 3, 1) =  0.70000D-03
        HBB(11, 4, 1) = -0.27200D-02
        HBB(11, 5, 1) = -0.30000D-02
        HBB(11, 6, 1) =  0.57200D-02
        HBB(11, 7, 1) =  0.20110D-01
        HBB(11, 8, 1) = -0.74900D-02
        HBB(11, 9, 1) = -0.12620D-01
        HBB(11,10, 1) = -0.37620D-01
        HBB(11,11, 1) =  0.78370D-01
        HBB(12, 1, 1) = -0.20000D-02
        HBB(12, 2, 1) =  0.30000D-03
        HBB(12, 3, 1) =  0.17000D-02
        HBB(12, 4, 1) =  0.84500D-02
        HBB(12, 5, 1) = -0.58700D-02
        HBB(12, 6, 1) = -0.25800D-02
        HBB(12, 7, 1) =  0.13390D-01
        HBB(12, 8, 1) = -0.60000D-02
        HBB(12, 9, 1) = -0.73900D-02
        HBB(12,10, 1) = -0.37460D-01
        HBB(12,11, 1) = -0.40760D-01
        HBB(12,12, 1) =  0.78220D-01
        HBB(13, 1, 1) = -0.61900D-02
        HBB(13, 2, 1) =  0.10630D-01
        HBB(13, 3, 1) = -0.44400D-02
        HBB(13, 4, 1) = -0.14000D-02
        HBB(13, 5, 1) =  0.11900D-02
        HBB(13, 6, 1) =  0.21000D-03
        HBB(13, 7, 1) =  0.29200D-02
        HBB(13, 8, 1) =  0.41800D-02
        HBB(13, 9, 1) = -0.71000D-02
        HBB(13,10, 1) = -0.33570D-01
        HBB(13,11, 1) =  0.13390D-01
        HBB(13,12, 1) =  0.20180D-01
        HBB(13,13, 1) =  0.65200D-01
        HBB(14, 1, 1) =  0.89900D-02
        HBB(14, 2, 1) = -0.58800D-02
        HBB(14, 3, 1) = -0.31100D-02
        HBB(14, 4, 1) =  0.77000D-03
        HBB(14, 5, 1) = -0.83000D-03
        HBB(14, 6, 1) =  0.70000D-04
        HBB(14, 7, 1) = -0.71700D-02
        HBB(14, 8, 1) = -0.51000D-03
        HBB(14, 9, 1) =  0.76800D-02
        HBB(14,10, 1) =  0.20000D-01
        HBB(14,11, 1) = -0.73300D-02
        HBB(14,12, 1) = -0.12670D-01
        HBB(14,13, 1) = -0.33080D-01
        HBB(14,14, 1) =  0.74340D-01
        HBB(15, 1, 1) = -0.28000D-02
        HBB(15, 2, 1) = -0.47500D-02
        HBB(15, 3, 1) =  0.75500D-02
        HBB(15, 4, 1) =  0.64000D-03
        HBB(15, 5, 1) = -0.36000D-03
        HBB(15, 6, 1) = -0.28000D-03
        HBB(15, 7, 1) =  0.42500D-02
        HBB(15, 8, 1) = -0.36700D-02
        HBB(15, 9, 1) = -0.58000D-03
        HBB(15,10, 1) =  0.13570D-01
        HBB(15,11, 1) = -0.60500D-02
        HBB(15,12, 1) = -0.75200D-02
        HBB(15,13, 1) = -0.32120D-01
        HBB(15,14, 1) = -0.41260D-01
        HBB(15,15, 1) =  0.73390D-01
        HBB(16, 1, 1) = -0.36890D-01
        HBB(16, 2, 1) =  0.15200D-01
        HBB(16, 3, 1) =  0.21690D-01
        HBB(16, 4, 1) =  0.32300D-02
        HBB(16, 5, 1) = -0.74000D-02
        HBB(16, 6, 1) =  0.41700D-02
        HBB(16, 7, 1) = -0.13900D-02
        HBB(16, 8, 1) =  0.69000D-03
        HBB(16, 9, 1) =  0.71000D-03
        HBB(16,10, 1) = -0.61100D-02
        HBB(16,11, 1) =  0.86800D-02
        HBB(16,12, 1) = -0.25700D-02
        HBB(16,13, 1) = -0.26960D-01
        HBB(16,14, 1) =  0.10490D-01
        HBB(16,15, 1) =  0.16470D-01
        HBB(16,16, 1) =  0.68110D-01
        HBB(17, 1, 1) =  0.22720D-01
        HBB(17, 2, 1) = -0.69100D-02
        HBB(17, 3, 1) = -0.15810D-01
        HBB(17, 4, 1) = -0.78500D-02
        HBB(17, 5, 1) =  0.81700D-02
        HBB(17, 6, 1) = -0.32000D-03
        HBB(17, 7, 1) =  0.10300D-02
        HBB(17, 8, 1) = -0.32000D-03
        HBB(17, 9, 1) = -0.71000D-03
        HBB(17,10, 1) =  0.90700D-02
        HBB(17,11, 1) = -0.60300D-02
        HBB(17,12, 1) = -0.30400D-02
        HBB(17,13, 1) =  0.97900D-02
        HBB(17,14, 1) = -0.45700D-02
        HBB(17,15, 1) = -0.52200D-02
        HBB(17,16, 1) = -0.34760D-01
        HBB(17,17, 1) =  0.72170D-01
        HBB(18, 1, 1) =  0.14170D-01
        HBB(18, 2, 1) = -0.82800D-02
        HBB(18, 3, 1) = -0.58800D-02
        HBB(18, 4, 1) =  0.46200D-02
        HBB(18, 5, 1) = -0.77000D-03
        HBB(18, 6, 1) = -0.38500D-02
        HBB(18, 7, 1) =  0.36000D-03
        HBB(18, 8, 1) = -0.37000D-03
        HBB(18, 9, 1) =  0.10000D-04
        HBB(18,10, 1) = -0.29600D-02
        HBB(18,11, 1) = -0.26500D-02
        HBB(18,12, 1) =  0.56000D-02
        HBB(18,13, 1) =  0.17170D-01
        HBB(18,14, 1) = -0.59200D-02
        HBB(18,15, 1) = -0.11250D-01
        HBB(18,16, 1) = -0.33360D-01
        HBB(18,17, 1) = -0.37420D-01
        HBB(18,18, 1) =  0.70770D-01
C Hessian elements for anchor point 1
        HTT( 1, 1, 1) =  0.59700D-02
        HTT( 2, 1, 1) =  0.40300D-02
        HTT( 2, 2, 1) =  0.53700D-02
        HTT( 3, 1, 1) = -0.14500D-02
        HTT( 3, 2, 1) = -0.35400D-02
        HTT( 3, 3, 1) =  0.68100D-02
        HTT( 4, 1, 1) = -0.33900D-02
        HTT( 4, 2, 1) = -0.21900D-02
        HTT( 4, 3, 1) =  0.47200D-02
        HTT( 4, 4, 1) =  0.59200D-02
        HTT( 5, 1, 1) = -0.57900D-02
        HTT( 5, 2, 1) = -0.45300D-02
        HTT( 5, 3, 1) =  0.17600D-02
        HTT( 5, 4, 1) =  0.30200D-02
        HTT( 5, 5, 1) =  0.64900D-02
        HTT( 6, 1, 1) = -0.38500D-02
        HTT( 6, 2, 1) = -0.41000D-02
        HTT( 6, 3, 1) =  0.21800D-02
        HTT( 6, 4, 1) =  0.19300D-02
        HTT( 6, 5, 1) =  0.34500D-02
        HTT( 6, 6, 1) =  0.52100D-02
        HTT( 7, 1, 1) =  0.13600D-02
        HTT( 7, 2, 1) =  0.27600D-02
        HTT( 7, 3, 1) = -0.62000D-02
        HTT( 7, 4, 1) = -0.48000D-02
        HTT( 7, 5, 1) = -0.78000D-03
        HTT( 7, 6, 1) = -0.23700D-02
        HTT( 7, 7, 1) =  0.65100D-02
        HTT( 8, 1, 1) =  0.32900D-02
        HTT( 8, 2, 1) =  0.31900D-02
        HTT( 8, 3, 1) = -0.57800D-02
        HTT( 8, 4, 1) = -0.58800D-02
        HTT( 8, 5, 1) = -0.38200D-02
        HTT( 8, 6, 1) = -0.61000D-03
        HTT( 8, 7, 1) =  0.49200D-02
        HTT( 8, 8, 1) =  0.81300D-02
        HTT( 9, 1, 1) = -0.23500D-02
        HTT( 9, 2, 1) = -0.73000D-03
        HTT( 9, 3, 1) =  0.17000D-03
        HTT( 9, 4, 1) =  0.17900D-02
        HTT( 9, 5, 1) =  0.15900D-02
        HTT( 9, 6, 1) =  0.12400D-02
        HTT( 9, 7, 1) = -0.84000D-03
        HTT( 9, 8, 1) = -0.11900D-02
        HTT( 9, 9, 1) =  0.25000D-02
        HTT(10, 1, 1) = -0.13700D-02
        HTT(10, 2, 1) =  0.15000D-03
        HTT(10, 3, 1) = -0.60000D-04
        HTT(10, 4, 1) =  0.14600D-02
        HTT(10, 5, 1) =  0.15100D-02
        HTT(10, 6, 1) =  0.46000D-03
        HTT(10, 7, 1) =  0.25000D-03
        HTT(10, 8, 1) = -0.80000D-03
        HTT(10, 9, 1) =  0.46000D-03
        HTT(10,10, 1) =  0.37200D-02
        HTT(11, 1, 1) = -0.36000D-03
        HTT(11, 2, 1) = -0.21100D-02
        HTT(11, 3, 1) =  0.23100D-02
        HTT(11, 4, 1) =  0.56000D-03
        HTT(11, 5, 1) =  0.29000D-03
        HTT(11, 6, 1) =  0.15000D-02
        HTT(11, 7, 1) = -0.22800D-02
        HTT(11, 8, 1) = -0.10800D-02
        HTT(11, 9, 1) =  0.83000D-03
        HTT(11,10, 1) = -0.11000D-02
        HTT(11,11, 1) =  0.26300D-02
        HTT(12, 1, 1) =  0.62000D-03
        HTT(12, 2, 1) = -0.12400D-02
        HTT(12, 3, 1) =  0.20800D-02
        HTT(12, 4, 1) =  0.22000D-03
        HTT(12, 5, 1) =  0.21000D-03
        HTT(12, 6, 1) =  0.72000D-03
        HTT(12, 7, 1) = -0.11900D-02
        HTT(12, 8, 1) = -0.69000D-03
        HTT(12, 9, 1) = -0.12000D-02
        HTT(12,10, 1) =  0.21600D-02
        HTT(12,11, 1) =  0.70000D-03
        HTT(12,12, 1) =  0.40700D-02
        HTT(13, 1, 1) = -0.11400D-02
        HTT(13, 2, 1) = -0.19000D-02
        HTT(13, 3, 1) =  0.76000D-03
        HTT(13, 4, 1) = -0.10000D-04
        HTT(13, 5, 1) =  0.16400D-02
        HTT(13, 6, 1) =  0.16200D-02
        HTT(13, 7, 1) = -0.19000D-03
        HTT(13, 8, 1) = -0.20000D-03
        HTT(13, 9, 1) = -0.21000D-02
        HTT(13,10, 1) =  0.24000D-03
        HTT(13,11, 1) = -0.13100D-02
        HTT(13,12, 1) =  0.10300D-02
        HTT(13,13, 1) =  0.49700D-02
        HTT(14, 1, 1) =  0.70000D-03
        HTT(14, 2, 1) =  0.47000D-03
        HTT(14, 3, 1) = -0.32000D-03
        HTT(14, 4, 1) = -0.55000D-03
        HTT(14, 5, 1) = -0.18000D-03
        HTT(14, 6, 1) = -0.82000D-03
        HTT(14, 7, 1) =  0.80000D-03
        HTT(14, 8, 1) =  0.16000D-03
        HTT(14, 9, 1) = -0.22000D-03
        HTT(14,10, 1) =  0.13600D-02
        HTT(14,11, 1) =  0.10000D-04
        HTT(14,12, 1) =  0.15900D-02
        HTT(14,13, 1) = -0.78000D-03
        HTT(14,14, 1) =  0.38600D-02
        HTT(15, 1, 1) = -0.21100D-02
        HTT(15, 2, 1) = -0.27800D-02
        HTT(15, 3, 1) =  0.98000D-03
        HTT(15, 4, 1) =  0.32000D-03
        HTT(15, 5, 1) =  0.17100D-02
        HTT(15, 6, 1) =  0.23900D-02
        HTT(15, 7, 1) = -0.12700D-02
        HTT(15, 8, 1) = -0.59000D-03
        HTT(15, 9, 1) = -0.80000D-04
        HTT(15,10, 1) = -0.29900D-02
        HTT(15,11, 1) =  0.61000D-03
        HTT(15,12, 1) = -0.23100D-02
        HTT(15,13, 1) =  0.26500D-02
        HTT(15,14, 1) = -0.23400D-02
        HTT(15,15, 1) =  0.55400D-02
        HTT(16, 1, 1) = -0.27000D-03
        HTT(16, 2, 1) = -0.40000D-03
        HTT(16, 3, 1) = -0.90000D-04
        HTT(16, 4, 1) = -0.22000D-03
        HTT(16, 5, 1) = -0.10000D-03
        HTT(16, 6, 1) = -0.50000D-04
        HTT(16, 7, 1) = -0.28000D-03
        HTT(16, 8, 1) = -0.22000D-03
        HTT(16, 9, 1) =  0.18000D-02
        HTT(16,10, 1) = -0.18700D-02
        HTT(16,11, 1) =  0.19300D-02
        HTT(16,12, 1) = -0.17400D-02
        HTT(16,13, 1) = -0.30900D-02
        HTT(16,14, 1) =  0.23000D-02
        HTT(16,15, 1) =  0.55000D-03
        HTT(16,16, 1) =  0.59400D-02
        HTT(17, 1, 1) =  0.13400D-02
        HTT(17, 2, 1) =  0.14100D-02
        HTT(17, 3, 1) = -0.45000D-03
        HTT(17, 4, 1) = -0.38000D-03
        HTT(17, 5, 1) = -0.94000D-03
        HTT(17, 6, 1) = -0.20500D-02
        HTT(17, 7, 1) =  0.78000D-03
        HTT(17, 8, 1) = -0.32000D-03
        HTT(17, 9, 1) =  0.13000D-02
        HTT(17,10, 1) = -0.10000D-03
        HTT(17,11, 1) =  0.12300D-02
        HTT(17,12, 1) = -0.17000D-03
        HTT(17,13, 1) = -0.44500D-02
        HTT(17,14, 1) =  0.13100D-02
        HTT(17,15, 1) = -0.30600D-02
        HTT(17,16, 1) =  0.27000D-02
        HTT(17,17, 1) =  0.48400D-02
        HTT(18, 1, 1) =  0.14900D-02
        HTT(18, 2, 1) =  0.21600D-02
        HTT(18, 3, 1) = -0.11900D-02
        HTT(18, 4, 1) = -0.51000D-03
        HTT(18, 5, 1) = -0.18100D-02
        HTT(18, 6, 1) = -0.28100D-02
        HTT(18, 7, 1) =  0.76000D-03
        HTT(18, 8, 1) = -0.23000D-03
        HTT(18, 9, 1) =  0.98000D-03
        HTT(18,10, 1) =  0.75000D-03
        HTT(18,11, 1) =  0.29000D-03
        HTT(18,12, 1) =  0.60000D-04
        HTT(18,13, 1) = -0.31600D-02
        HTT(18,14, 1) =  0.18600D-02
        HTT(18,15, 1) = -0.29300D-02
        HTT(18,16, 1) =  0.20900D-02
        HTT(18,17, 1) =  0.28200D-02
        HTT(18,18, 1) =  0.55200D-02
        HTT(19, 1, 1) = -0.50000D-03
        HTT(19, 2, 1) = -0.97000D-03
        HTT(19, 3, 1) =  0.63000D-03
        HTT(19, 4, 1) =  0.15000D-03
        HTT(19, 5, 1) =  0.88000D-03
        HTT(19, 6, 1) =  0.40000D-03
        HTT(19, 7, 1) = -0.21000D-03
        HTT(19, 8, 1) = -0.69000D-03
        HTT(19, 9, 1) = -0.59000D-03
        HTT(19,10, 1) = -0.12200D-02
        HTT(19,11, 1) = -0.10000D-03
        HTT(19,12, 1) = -0.74000D-03
        HTT(19,13, 1) =  0.13200D-02
        HTT(19,14, 1) = -0.33400D-02
        HTT(19,15, 1) =  0.19500D-02
        HTT(19,16, 1) = -0.27100D-02
        HTT(19,17, 1) = -0.93000D-03
        HTT(19,18, 1) = -0.22100D-02
        HTT(19,19, 1) =  0.37400D-02
        HTT(20, 1, 1) = -0.35000D-03
        HTT(20, 2, 1) = -0.22000D-03
        HTT(20, 3, 1) = -0.11000D-03
        HTT(20, 4, 1) =  0.20000D-04
        HTT(20, 5, 1) =  0.10000D-04
        HTT(20, 6, 1) = -0.36000D-03
        HTT(20, 7, 1) = -0.23000D-03
        HTT(20, 8, 1) = -0.60000D-03
        HTT(20, 9, 1) = -0.91000D-03
        HTT(20,10, 1) = -0.38000D-03
        HTT(20,11, 1) = -0.10400D-02
        HTT(20,12, 1) = -0.51000D-03
        HTT(20,13, 1) =  0.26000D-02
        HTT(20,14, 1) = -0.27900D-02
        HTT(20,15, 1) =  0.20800D-02
        HTT(20,16, 1) = -0.33200D-02
        HTT(20,17, 1) = -0.29500D-02
        HTT(20,18, 1) =  0.49000D-03
        HTT(20,19, 1) =  0.24600D-02
        HTT(20,20, 1) =  0.59000D-02
        HTT(21, 1, 1) =  0.19800D-02
        HTT(21, 2, 1) =  0.17000D-02
        HTT(21, 3, 1) = -0.76000D-03
        HTT(21, 4, 1) = -0.10500D-02
        HTT(21, 5, 1) = -0.29500D-02
        HTT(21, 6, 1) = -0.43000D-03
        HTT(21, 7, 1) = -0.30000D-03
        HTT(21, 8, 1) =  0.22100D-02
        HTT(21, 9, 1) = -0.98000D-03
        HTT(21,10, 1) = -0.73000D-03
        HTT(21,11, 1) = -0.69000D-03
        HTT(21,12, 1) = -0.44000D-03
        HTT(21,13, 1) =  0.11100D-02
        HTT(21,14, 1) = -0.81000D-03
        HTT(21,15, 1) =  0.87000D-03
        HTT(21,16, 1) = -0.10600D-02
        HTT(21,17, 1) = -0.20800D-02
        HTT(21,18, 1) = -0.34000D-03
        HTT(21,19, 1) = -0.15000D-03
        HTT(21,20, 1) =  0.15900D-02
        HTT(21,21, 1) =  0.28700D-02
        HTT(22, 1, 1) = -0.50000D-04
        HTT(22, 2, 1) =  0.12500D-02
        HTT(22, 3, 1) = -0.12100D-02
        HTT(22, 4, 1) =  0.10000D-03
        HTT(22, 5, 1) =  0.25000D-03
        HTT(22, 6, 1) = -0.22800D-02
        HTT(22, 7, 1) =  0.13700D-02
        HTT(22, 8, 1) = -0.11700D-02
        HTT(22, 9, 1) = -0.62000D-03
        HTT(22,10, 1) =  0.37000D-03
        HTT(22,11, 1) = -0.19600D-02
        HTT(22,12, 1) = -0.97000D-03
        HTT(22,13, 1) =  0.11300D-02
        HTT(22,14, 1) = -0.13000D-03
        HTT(22,15, 1) =  0.15000D-03
        HTT(22,16, 1) = -0.11100D-02
        HTT(22,17, 1) = -0.92000D-03
        HTT(22,18, 1) =  0.70000D-03
        HTT(22,19, 1) =  0.35000D-03
        HTT(22,20, 1) =  0.19700D-02
        HTT(22,21, 1) =  0.22000D-03
        HTT(22,22, 1) =  0.28800D-02
        HTT(23, 1, 1) =  0.18300D-02
        HTT(23, 2, 1) =  0.94000D-03
        HTT(23, 3, 1) = -0.20000D-04
        HTT(23, 4, 1) = -0.91000D-03
        HTT(23, 5, 1) = -0.20700D-02
        HTT(23, 6, 1) =  0.34000D-03
        HTT(23, 7, 1) = -0.28000D-03
        HTT(23, 8, 1) =  0.21200D-02
        HTT(23, 9, 1) = -0.66000D-03
        HTT(23,10, 1) = -0.15900D-02
        HTT(23,11, 1) =  0.26000D-03
        HTT(23,12, 1) = -0.67000D-03
        HTT(23,13, 1) = -0.18000D-03
        HTT(23,14, 1) = -0.13600D-02
        HTT(23,15, 1) =  0.74000D-03
        HTT(23,16, 1) = -0.44000D-03
        HTT(23,17, 1) = -0.50000D-04
        HTT(23,18, 1) = -0.30600D-02
        HTT(23,19, 1) =  0.11300D-02
        HTT(23,20, 1) = -0.18700D-02
        HTT(23,21, 1) =  0.11100D-02
        HTT(23,22, 1) = -0.14200D-02
        HTT(23,23, 1) =  0.41400D-02
        HTT(24, 1, 1) = -0.20000D-03
        HTT(24, 2, 1) =  0.49000D-03
        HTT(24, 3, 1) = -0.47000D-03
        HTT(24, 4, 1) =  0.23000D-03
        HTT(24, 5, 1) =  0.11300D-02
        HTT(24, 6, 1) = -0.15100D-02
        HTT(24, 7, 1) =  0.13900D-02
        HTT(24, 8, 1) = -0.12600D-02
        HTT(24, 9, 1) = -0.29000D-03
        HTT(24,10, 1) = -0.48000D-03
        HTT(24,11, 1) = -0.10100D-02
        HTT(24,12, 1) = -0.12000D-02
        HTT(24,13, 1) = -0.17000D-03
        HTT(24,14, 1) = -0.69000D-03
        HTT(24,15, 1) =  0.20000D-04
        HTT(24,16, 1) = -0.50000D-03
        HTT(24,17, 1) =  0.11100D-02
        HTT(24,18, 1) = -0.20100D-02
        HTT(24,19, 1) =  0.16300D-02
        HTT(24,20, 1) = -0.14900D-02
        HTT(24,21, 1) = -0.15300D-02
        HTT(24,22, 1) =  0.12400D-02
        HTT(24,23, 1) =  0.16100D-02
        HTT(24,24, 1) =  0.43900D-02
C Hessian elements for anchor point 1
        HBS( 1, 1, 1) = -0.30890D-01
        HBS( 1, 2, 1) = -0.30720D-01
        HBS( 1, 3, 1) = -0.30470D-01
        HBS( 1, 4, 1) = -0.89900D-02
        HBS( 1, 5, 1) =  0.47600D-02
        HBS( 1, 6, 1) =  0.39030D-01
        HBS( 1, 7, 1) = -0.44000D-03
        HBS( 1, 8, 1) =  0.39840D-01
        HBS( 1, 9, 1) = -0.25000D-03
        HBS( 1,10, 1) = -0.97500D-02
        HBS( 1,11, 1) = -0.42000D-03
        HBS( 1,12, 1) =  0.49400D-02
        HBS( 2, 1, 1) =  0.43860D-01
        HBS( 2, 2, 1) = -0.11950D-01
        HBS( 2, 3, 1) =  0.67900D-02
        HBS( 2, 4, 1) =  0.23800D-02
        HBS( 2, 5, 1) = -0.22800D-02
        HBS( 2, 6, 1) = -0.26060D-01
        HBS( 2, 7, 1) =  0.11000D-02
        HBS( 2, 8, 1) = -0.14350D-01
        HBS( 2, 9, 1) =  0.40000D-04
        HBS( 2,10, 1) =  0.41400D-02
        HBS( 2,11, 1) = -0.73000D-03
        HBS( 2,12, 1) = -0.48000D-03
        HBS( 3, 1, 1) = -0.12970D-01
        HBS( 3, 2, 1) =  0.42680D-01
        HBS( 3, 3, 1) =  0.23690D-01
        HBS( 3, 4, 1) =  0.66100D-02
        HBS( 3, 5, 1) = -0.24900D-02
        HBS( 3, 6, 1) = -0.12970D-01
        HBS( 3, 7, 1) = -0.66000D-03
        HBS( 3, 8, 1) = -0.25490D-01
        HBS( 3, 9, 1) =  0.21000D-03
        HBS( 3,10, 1) =  0.56200D-02
        HBS( 3,11, 1) =  0.11500D-02
        HBS( 3,12, 1) = -0.44600D-02
        HBS( 4, 1, 1) = -0.21190D-01
        HBS( 4, 2, 1) = -0.10770D-01
        HBS( 4, 3, 1) =  0.12660D-01
        HBS( 4, 4, 1) = -0.29880D-01
        HBS( 4, 5, 1) = -0.90800D-02
        HBS( 4, 6, 1) = -0.51700D-02
        HBS( 4, 7, 1) =  0.49300D-02
        HBS( 4, 8, 1) =  0.30850D-01
        HBS( 4, 9, 1) = -0.40000D-03
        HBS( 4,10, 1) =  0.40050D-01
        HBS( 4,11, 1) = -0.62000D-03
        HBS( 4,12, 1) =  0.12000D-03
        HBS( 5, 1, 1) =  0.25060D-01
        HBS( 5, 2, 1) =  0.46900D-02
        HBS( 5, 3, 1) = -0.97000D-02
        HBS( 5, 4, 1) =  0.24100D-02
        HBS( 5, 5, 1) =  0.62100D-02
        HBS( 5, 6, 1) =  0.30800D-02
        HBS( 5, 7, 1) = -0.21300D-02
        HBS( 5, 8, 1) = -0.14880D-01
        HBS( 5, 9, 1) = -0.39000D-03
        HBS( 5,10, 1) = -0.21940D-01
        HBS( 5,11, 1) =  0.31000D-03
        HBS( 5,12, 1) =  0.35000D-03
        HBS( 6, 1, 1) = -0.38700D-02
        HBS( 6, 2, 1) =  0.60800D-02
        HBS( 6, 3, 1) = -0.29500D-02
        HBS( 6, 4, 1) =  0.27480D-01
        HBS( 6, 5, 1) =  0.28700D-02
        HBS( 6, 6, 1) =  0.21000D-02
        HBS( 6, 7, 1) = -0.28000D-02
        HBS( 6, 8, 1) = -0.15960D-01
        HBS( 6, 9, 1) =  0.80000D-03
        HBS( 6,10, 1) = -0.18110D-01
        HBS( 6,11, 1) =  0.31000D-03
        HBS( 6,12, 1) = -0.47000D-03
        HBS( 7, 1, 1) = -0.88700D-02
        HBS( 7, 2, 1) =  0.35920D-01
        HBS( 7, 3, 1) =  0.96000D-03
        HBS( 7, 4, 1) = -0.23260D-01
        HBS( 7, 5, 1) =  0.57800D-02
        HBS( 7, 6, 1) = -0.29090D-01
        HBS( 7, 7, 1) = -0.87400D-02
        HBS( 7, 8, 1) = -0.88100D-02
        HBS( 7, 9, 1) =  0.49600D-02
        HBS( 7,10, 1) =  0.32900D-01
        HBS( 7,11, 1) =  0.60000D-04
        HBS( 7,12, 1) = -0.10800D-02
        HBS( 8, 1, 1) =  0.36600D-02
        HBS( 8, 2, 1) = -0.18680D-01
        HBS( 8, 3, 1) =  0.82000D-03
        HBS( 8, 4, 1) =  0.25710D-01
        HBS( 8, 5, 1) = -0.33200D-02
        HBS( 8, 6, 1) = -0.40000D-03
        HBS( 8, 7, 1) =  0.44400D-02
        HBS( 8, 8, 1) =  0.48100D-02
        HBS( 8, 9, 1) = -0.20200D-02
        HBS( 8,10, 1) = -0.14750D-01
        HBS( 8,11, 1) = -0.62000D-03
        HBS( 8,12, 1) =  0.56000D-03
        HBS( 9, 1, 1) =  0.52100D-02
        HBS( 9, 2, 1) = -0.17240D-01
        HBS( 9, 3, 1) = -0.17800D-02
        HBS( 9, 4, 1) = -0.24500D-02
        HBS( 9, 5, 1) = -0.24600D-02
        HBS( 9, 6, 1) =  0.29490D-01
        HBS( 9, 7, 1) =  0.43000D-02
        HBS( 9, 8, 1) =  0.40000D-02
        HBS( 9, 9, 1) = -0.29400D-02
        HBS( 9,10, 1) = -0.18150D-01
        HBS( 9,11, 1) =  0.56000D-03
        HBS( 9,12, 1) =  0.52000D-03
        HBS(10, 1, 1) =  0.36080D-01
        HBS(10, 2, 1) =  0.38610D-01
        HBS(10, 3, 1) =  0.39500D-02
        HBS(10, 4, 1) = -0.11320D-01
        HBS(10, 5, 1) = -0.89000D-03
        HBS(10, 6, 1) = -0.26740D-01
        HBS(10, 7, 1) =  0.48000D-02
        HBS(10, 8, 1) = -0.26910D-01
        HBS(10, 9, 1) = -0.88100D-02
        HBS(10,10, 1) = -0.11650D-01
        HBS(10,11, 1) =  0.47600D-02
        HBS(10,12, 1) = -0.70000D-03
        HBS(11, 1, 1) = -0.20240D-01
        HBS(11, 2, 1) = -0.17250D-01
        HBS(11, 3, 1) = -0.21300D-02
        HBS(11, 4, 1) =  0.55100D-02
        HBS(11, 5, 1) =  0.11200D-02
        HBS(11, 6, 1) =  0.27050D-01
        HBS(11, 7, 1) = -0.28200D-02
        HBS(11, 8, 1) =  0.25000D-03
        HBS(11, 9, 1) =  0.44500D-02
        HBS(11,10, 1) =  0.58600D-02
        HBS(11,11, 1) = -0.19800D-02
        HBS(11,12, 1) = -0.20000D-03
        HBS(12, 1, 1) = -0.15840D-01
        HBS(12, 2, 1) = -0.21360D-01
        HBS(12, 3, 1) = -0.18200D-02
        HBS(12, 4, 1) =  0.58100D-02
        HBS(12, 5, 1) = -0.23000D-03
        HBS(12, 6, 1) = -0.31000D-03
        HBS(12, 7, 1) = -0.19800D-02
        HBS(12, 8, 1) =  0.26660D-01
        HBS(12, 9, 1) =  0.43600D-02
        HBS(12,10, 1) =  0.57900D-02
        HBS(12,11, 1) = -0.27800D-02
        HBS(12,12, 1) =  0.90000D-03
        HBS(13, 1, 1) =  0.35090D-01
        HBS(13, 2, 1) = -0.85500D-02
        HBS(13, 3, 1) =  0.10200D-02
        HBS(13, 4, 1) =  0.32180D-01
        HBS(13, 5, 1) = -0.97000D-03
        HBS(13, 6, 1) = -0.94600D-02
        HBS(13, 7, 1) =  0.10000D-03
        HBS(13, 8, 1) = -0.29520D-01
        HBS(13, 9, 1) =  0.50200D-02
        HBS(13,10, 1) = -0.23910D-01
        HBS(13,11, 1) = -0.86300D-02
        HBS(13,12, 1) =  0.52600D-02
        HBS(14, 1, 1) = -0.16950D-01
        HBS(14, 2, 1) =  0.51900D-02
        HBS(14, 3, 1) = -0.19500D-02
        HBS(14, 4, 1) = -0.17590D-01
        HBS(14, 5, 1) =  0.48000D-03
        HBS(14, 6, 1) =  0.43300D-02
        HBS(14, 7, 1) =  0.54000D-03
        HBS(14, 8, 1) =  0.29620D-01
        HBS(14, 9, 1) = -0.29700D-02
        HBS(14,10, 1) = -0.21100D-02
        HBS(14,11, 1) =  0.42400D-02
        HBS(14,12, 1) = -0.22400D-02
        HBS(15, 1, 1) = -0.18140D-01
        HBS(15, 2, 1) =  0.33600D-02
        HBS(15, 3, 1) =  0.94000D-03
        HBS(15, 4, 1) = -0.14580D-01
        HBS(15, 5, 1) =  0.48000D-03
        HBS(15, 6, 1) =  0.51200D-02
        HBS(15, 7, 1) = -0.64000D-03
        HBS(15, 8, 1) = -0.10000D-03
        HBS(15, 9, 1) = -0.20500D-02
        HBS(15,10, 1) =  0.26020D-01
        HBS(15,11, 1) =  0.43900D-02
        HBS(15,12, 1) = -0.30200D-02
        HBS(16, 1, 1) = -0.10210D-01
        HBS(16, 2, 1) = -0.24490D-01
        HBS(16, 3, 1) =  0.11890D-01
        HBS(16, 4, 1) =  0.41280D-01
        HBS(16, 5, 1) =  0.40000D-03
        HBS(16, 6, 1) =  0.31430D-01
        HBS(16, 7, 1) = -0.65000D-03
        HBS(16, 8, 1) = -0.54400D-02
        HBS(16, 9, 1) = -0.51000D-03
        HBS(16,10, 1) = -0.27630D-01
        HBS(16,11, 1) =  0.48400D-02
        HBS(16,12, 1) = -0.85500D-02
        HBS(17, 1, 1) =  0.43200D-02
        HBS(17, 2, 1) =  0.25240D-01
        HBS(17, 3, 1) = -0.81200D-02
        HBS(17, 4, 1) = -0.22630D-01
        HBS(17, 5, 1) =  0.41000D-03
        HBS(17, 6, 1) = -0.14630D-01
        HBS(17, 7, 1) =  0.33000D-03
        HBS(17, 8, 1) =  0.32400D-02
        HBS(17, 9, 1) = -0.33000D-03
        HBS(17,10, 1) =  0.18800D-02
        HBS(17,11, 1) = -0.21000D-02
        HBS(17,12, 1) =  0.69700D-02
        HBS(18, 1, 1) =  0.58800D-02
        HBS(18, 2, 1) = -0.75000D-03
        HBS(18, 3, 1) = -0.37700D-02
        HBS(18, 4, 1) = -0.18650D-01
        HBS(18, 5, 1) = -0.80000D-03
        HBS(18, 6, 1) = -0.16800D-01
        HBS(18, 7, 1) =  0.32000D-03
        HBS(18, 8, 1) =  0.22000D-02
        HBS(18, 9, 1) =  0.85000D-03
        HBS(18,10, 1) =  0.25750D-01
        HBS(18,11, 1) = -0.27400D-02
        HBS(18,12, 1) =  0.15700D-02
C Hessian elements for anchor point 2
        HSS( 1, 1, 2) =  0.34759D+00
        HSS( 2, 1, 2) =  0.17810D-01
        HSS( 2, 2, 2) =  0.35623D+00
        HSS( 3, 1, 2) =  0.29010D-01
        HSS( 3, 2, 2) =  0.32710D-01
        HSS( 3, 3, 2) =  0.46124D+00
        HSS( 4, 1, 2) = -0.19600D-02
        HSS( 4, 2, 2) =  0.14870D-01
        HSS( 4, 3, 2) = -0.89100D-02
        HSS( 4, 4, 2) =  0.38581D+00
        HSS( 5, 1, 2) =  0.38800D-02
        HSS( 5, 2, 2) = -0.77000D-03
        HSS( 5, 3, 2) =  0.71700D-02
        HSS( 5, 4, 2) =  0.36500D-02
        HSS( 5, 5, 2) =  0.35807D+00
        HSS( 6, 1, 2) =  0.33890D-01
        HSS( 6, 2, 2) =  0.44570D-01
        HSS( 6, 3, 2) = -0.43100D-02
        HSS( 6, 4, 2) = -0.36100D-02
        HSS( 6, 5, 2) =  0.55000D-03
        HSS( 6, 6, 2) =  0.35931D+00
        HSS( 7, 1, 2) = -0.54000D-03
        HSS( 7, 2, 2) = -0.63000D-03
        HSS( 7, 3, 2) =  0.57000D-03
        HSS( 7, 4, 2) =  0.39700D-02
        HSS( 7, 5, 2) =  0.39000D-03
        HSS( 7, 6, 2) =  0.39800D-02
        HSS( 7, 7, 2) =  0.37091D+00
        HSS( 8, 1, 2) =  0.42030D-01
        HSS( 8, 2, 2) =  0.35060D-01
        HSS( 8, 3, 2) = -0.52800D-02
        HSS( 8, 4, 2) =  0.24280D-01
        HSS( 8, 5, 2) = -0.11100D-02
        HSS( 8, 6, 2) =  0.17880D-01
        HSS( 8, 7, 2) = -0.10400D-02
        HSS( 8, 8, 2) =  0.35888D+00
        HSS( 9, 1, 2) = -0.53000D-03
        HSS( 9, 2, 2) = -0.57000D-03
        HSS( 9, 3, 2) =  0.26000D-03
        HSS( 9, 4, 2) = -0.39000D-03
        HSS( 9, 5, 2) =  0.19000D-03
        HSS( 9, 6, 2) =  0.40100D-02
        HSS( 9, 7, 2) =  0.49000D-03
        HSS( 9, 8, 2) =  0.39600D-02
        HSS( 9, 9, 2) =  0.36356D+00
        HSS(10, 1, 2) =  0.19520D-01
        HSS(10, 2, 2) = -0.69300D-02
        HSS(10, 3, 2) = -0.48800D-02
        HSS(10, 4, 2) =  0.64010D-01
        HSS(10, 5, 2) = -0.43000D-03
        HSS(10, 6, 2) =  0.21230D-01
        HSS(10, 7, 2) = -0.10900D-02
        HSS(10, 8, 2) = -0.77000D-03
        HSS(10, 9, 2) = -0.19000D-03
        HSS(10,10, 2) =  0.38116D+00
        HSS(11, 1, 2) = -0.69000D-03
        HSS(11, 2, 2) = -0.50000D-03
        HSS(11, 3, 2) =  0.52000D-03
        HSS(11, 4, 2) = -0.11600D-02
        HSS(11, 5, 2) =  0.80000D-04
        HSS(11, 6, 2) = -0.10000D-02
        HSS(11, 7, 2) =  0.22000D-03
        HSS(11, 8, 2) =  0.39700D-02
        HSS(11, 9, 2) =  0.48000D-03
        HSS(11,10, 2) =  0.39800D-02
        HSS(11,11, 2) =  0.37152D+00
        HSS(12, 1, 2) = -0.16500D-02
        HSS(12, 2, 2) =  0.39400D-02
        HSS(12, 3, 2) = -0.19000D-03
        HSS(12, 4, 2) = -0.10000D-03
        HSS(12, 5, 2) =  0.35000D-03
        HSS(12, 6, 2) = -0.10400D-02
        HSS(12, 7, 2) =  0.30000D-04
        HSS(12, 8, 2) =  0.57000D-03
        HSS(12, 9, 2) =  0.12000D-03
        HSS(12,10, 2) =  0.23600D-02
        HSS(12,11, 2) =  0.40000D-03
        HSS(12,12, 2) =  0.37254D+00
C Hessian elements for anchor point 2
        HBB( 1, 1, 2) =  0.82450D-01
        HBB( 2, 1, 2) = -0.46090D-01
        HBB( 2, 2, 2) =  0.14487D+00
        HBB( 3, 1, 2) = -0.36360D-01
        HBB( 3, 2, 2) = -0.98780D-01
        HBB( 3, 3, 2) =  0.13515D+00
        HBB( 4, 1, 2) = -0.36000D-01
        HBB( 4, 2, 2) =  0.23310D-01
        HBB( 4, 3, 2) =  0.12690D-01
        HBB( 4, 4, 2) =  0.66170D-01
        HBB( 5, 1, 2) =  0.20240D-01
        HBB( 5, 2, 2) = -0.12480D-01
        HBB( 5, 3, 2) = -0.77700D-02
        HBB( 5, 4, 2) = -0.33030D-01
        HBB( 5, 5, 2) =  0.71950D-01
        HBB( 6, 1, 2) =  0.15760D-01
        HBB( 6, 2, 2) = -0.10840D-01
        HBB( 6, 3, 2) = -0.49200D-02
        HBB( 6, 4, 2) = -0.33140D-01
        HBB( 6, 5, 2) = -0.38920D-01
        HBB( 6, 6, 2) =  0.72060D-01
        HBB( 7, 1, 2) = -0.71500D-02
        HBB( 7, 2, 2) = -0.39000D-02
        HBB( 7, 3, 2) =  0.11050D-01
        HBB( 7, 4, 2) = -0.25740D-01
        HBB( 7, 5, 2) =  0.93100D-02
        HBB( 7, 6, 2) =  0.16430D-01
        HBB( 7, 7, 2) =  0.64420D-01
        HBB( 8, 1, 2) = -0.22800D-02
        HBB( 8, 2, 2) =  0.74300D-02
        HBB( 8, 3, 2) = -0.51500D-02
        HBB( 8, 4, 2) =  0.15690D-01
        HBB( 8, 5, 2) = -0.49500D-02
        HBB( 8, 6, 2) = -0.10740D-01
        HBB( 8, 7, 2) = -0.31490D-01
        HBB( 8, 8, 2) =  0.72950D-01
        HBB( 9, 1, 2) =  0.94300D-02
        HBB( 9, 2, 2) = -0.35300D-02
        HBB( 9, 3, 2) = -0.59100D-02
        HBB( 9, 4, 2) =  0.10050D-01
        HBB( 9, 5, 2) = -0.43500D-02
        HBB( 9, 6, 2) = -0.56900D-02
        HBB( 9, 7, 2) = -0.32930D-01
        HBB( 9, 8, 2) = -0.41460D-01
        HBB( 9, 9, 2) =  0.74390D-01
        HBB(10, 1, 2) =  0.45500D-02
        HBB(10, 2, 2) = -0.20200D-02
        HBB(10, 3, 2) = -0.25300D-02
        HBB(10, 4, 2) = -0.64100D-02
        HBB(10, 5, 2) =  0.91900D-02
        HBB(10, 6, 2) = -0.27700D-02
        HBB(10, 7, 2) = -0.33410D-01
        HBB(10, 8, 2) =  0.13410D-01
        HBB(10, 9, 2) =  0.20000D-01
        HBB(10,10, 2) =  0.75440D-01
        HBB(11, 1, 2) = -0.22100D-02
        HBB(11, 2, 2) =  0.12300D-02
        HBB(11, 3, 2) =  0.98000D-03
        HBB(11, 4, 2) = -0.23100D-02
        HBB(11, 5, 2) = -0.32200D-02
        HBB(11, 6, 2) =  0.55300D-02
        HBB(11, 7, 2) =  0.20020D-01
        HBB(11, 8, 2) = -0.74100D-02
        HBB(11, 9, 2) = -0.12610D-01
        HBB(11,10, 2) = -0.37790D-01
        HBB(11,11, 2) =  0.78360D-01
        HBB(12, 1, 2) = -0.23400D-02
        HBB(12, 2, 2) =  0.79000D-03
        HBB(12, 3, 2) =  0.15500D-02
        HBB(12, 4, 2) =  0.87200D-02
        HBB(12, 5, 2) = -0.59700D-02
        HBB(12, 6, 2) = -0.27500D-02
        HBB(12, 7, 2) =  0.13390D-01
        HBB(12, 8, 2) = -0.60000D-02
        HBB(12, 9, 2) = -0.73900D-02
        HBB(12,10, 2) = -0.37640D-01
        HBB(12,11, 2) = -0.40570D-01
        HBB(12,12, 2) =  0.78220D-01
        HBB(13, 1, 2) = -0.67600D-02
        HBB(13, 2, 2) =  0.10710D-01
        HBB(13, 3, 2) = -0.39500D-02
        HBB(13, 4, 2) = -0.11500D-02
        HBB(13, 5, 2) =  0.12900D-02
        HBB(13, 6, 2) = -0.14000D-03
        HBB(13, 7, 2) =  0.29100D-02
        HBB(13, 8, 2) =  0.41300D-02
        HBB(13, 9, 2) = -0.70300D-02
        HBB(13,10, 2) = -0.33760D-01
        HBB(13,11, 2) =  0.13510D-01
        HBB(13,12, 2) =  0.20240D-01
        HBB(13,13, 2) =  0.65300D-01
        HBB(14, 1, 2) =  0.91200D-02
        HBB(14, 2, 2) = -0.60500D-02
        HBB(14, 3, 2) = -0.30700D-02
        HBB(14, 4, 2) =  0.66000D-03
        HBB(14, 5, 2) = -0.85000D-03
        HBB(14, 6, 2) =  0.19000D-03
        HBB(14, 7, 2) = -0.71100D-02
        HBB(14, 8, 2) = -0.51000D-03
        HBB(14, 9, 2) =  0.76200D-02
        HBB(14,10, 2) =  0.20030D-01
        HBB(14,11, 2) = -0.73800D-02
        HBB(14,12, 2) = -0.12650D-01
        HBB(14,13, 2) = -0.33180D-01
        HBB(14,14, 2) =  0.74350D-01
        HBB(15, 1, 2) = -0.23600D-02
        HBB(15, 2, 2) = -0.46500D-02
        HBB(15, 3, 2) =  0.70200D-02
        HBB(15, 4, 2) =  0.49000D-03
        HBB(15, 5, 2) = -0.44000D-03
        HBB(15, 6, 2) = -0.50000D-04
        HBB(15, 7, 2) =  0.42000D-02
        HBB(15, 8, 2) = -0.36200D-02
        HBB(15, 9, 2) = -0.59000D-03
        HBB(15,10, 2) =  0.13730D-01
        HBB(15,11, 2) = -0.61300D-02
        HBB(15,12, 2) = -0.76000D-02
        HBB(15,13, 2) = -0.32120D-01
        HBB(15,14, 2) = -0.41170D-01
        HBB(15,15, 2) =  0.73290D-01
        HBB(16, 1, 2) = -0.37090D-01
        HBB(16, 2, 2) =  0.17990D-01
        HBB(16, 3, 2) =  0.19100D-01
        HBB(16, 4, 2) =  0.31400D-02
        HBB(16, 5, 2) = -0.70000D-02
        HBB(16, 6, 2) =  0.38600D-02
        HBB(16, 7, 2) = -0.10300D-02
        HBB(16, 8, 2) =  0.54000D-03
        HBB(16, 9, 2) =  0.48000D-03
        HBB(16,10, 2) = -0.64000D-02
        HBB(16,11, 2) =  0.87800D-02
        HBB(16,12, 2) = -0.23700D-02
        HBB(16,13, 2) = -0.26530D-01
        HBB(16,14, 2) =  0.10480D-01
        HBB(16,15, 2) =  0.16060D-01
        HBB(16,16, 2) =  0.67910D-01
        HBB(17, 1, 2) =  0.22610D-01
        HBB(17, 2, 2) = -0.90200D-02
        HBB(17, 3, 2) = -0.13590D-01
        HBB(17, 4, 2) = -0.76000D-02
        HBB(17, 5, 2) =  0.77900D-02
        HBB(17, 6, 2) = -0.19000D-03
        HBB(17, 7, 2) =  0.87000D-03
        HBB(17, 8, 2) = -0.25000D-03
        HBB(17, 9, 2) = -0.62000D-03
        HBB(17,10, 2) =  0.91600D-02
        HBB(17,11, 2) = -0.60700D-02
        HBB(17,12, 2) = -0.31000D-02
        HBB(17,13, 2) =  0.96800D-02
        HBB(17,14, 2) = -0.46000D-02
        HBB(17,15, 2) = -0.50800D-02
        HBB(17,16, 2) = -0.34720D-01
        HBB(17,17, 2) =  0.72390D-01
        HBB(18, 1, 2) =  0.14470D-01
        HBB(18, 2, 2) = -0.89700D-02
        HBB(18, 3, 2) = -0.55000D-02
        HBB(18, 4, 2) =  0.44600D-02
        HBB(18, 5, 2) = -0.79000D-03
        HBB(18, 6, 2) = -0.36700D-02
        HBB(18, 7, 2) =  0.16000D-03
        HBB(18, 8, 2) = -0.29000D-03
        HBB(18, 9, 2) =  0.13000D-03
        HBB(18,10, 2) = -0.27600D-02
        HBB(18,11, 2) = -0.27100D-02
        HBB(18,12, 2) =  0.54700D-02
        HBB(18,13, 2) =  0.16860D-01
        HBB(18,14, 2) = -0.58800D-02
        HBB(18,15, 2) = -0.10980D-01
        HBB(18,16, 2) = -0.33190D-01
        HBB(18,17, 2) = -0.37670D-01
        HBB(18,18, 2) =  0.70860D-01
C Hessian elements for anchor point 2
        HTT( 1, 1, 2) =  0.60400D-02
        HTT( 2, 1, 2) =  0.38400D-02
        HTT( 2, 2, 2) =  0.50900D-02
        HTT( 3, 1, 2) = -0.18500D-02
        HTT( 3, 2, 2) = -0.33100D-02
        HTT( 3, 3, 2) =  0.75100D-02
        HTT( 4, 1, 2) = -0.40500D-02
        HTT( 4, 2, 2) = -0.20700D-02
        HTT( 4, 3, 2) =  0.60500D-02
        HTT( 4, 4, 2) =  0.80200D-02
        HTT( 5, 1, 2) = -0.59200D-02
        HTT( 5, 2, 2) = -0.42600D-02
        HTT( 5, 3, 2) =  0.21200D-02
        HTT( 5, 4, 2) =  0.37800D-02
        HTT( 5, 5, 2) =  0.65600D-02
        HTT( 6, 1, 2) = -0.41200D-02
        HTT( 6, 2, 2) = -0.39000D-02
        HTT( 6, 3, 2) =  0.26900D-02
        HTT( 6, 4, 2) =  0.29100D-02
        HTT( 6, 5, 2) =  0.37400D-02
        HTT( 6, 6, 2) =  0.55900D-02
        HTT( 7, 1, 2) =  0.15700D-02
        HTT( 7, 2, 2) =  0.25400D-02
        HTT( 7, 3, 2) = -0.67700D-02
        HTT( 7, 4, 2) = -0.58000D-02
        HTT( 7, 5, 2) = -0.10700D-02
        HTT( 7, 6, 2) = -0.27200D-02
        HTT( 7, 7, 2) =  0.68600D-02
        HTT( 8, 1, 2) =  0.33800D-02
        HTT( 8, 2, 2) =  0.29000D-02
        HTT( 8, 3, 2) = -0.62000D-02
        HTT( 8, 4, 2) = -0.66800D-02
        HTT( 8, 5, 2) = -0.38900D-02
        HTT( 8, 6, 2) = -0.88000D-03
        HTT( 8, 7, 2) =  0.52100D-02
        HTT( 8, 8, 2) =  0.82200D-02
        HTT( 9, 1, 2) = -0.23900D-02
        HTT( 9, 2, 2) = -0.68000D-03
        HTT( 9, 3, 2) =  0.45000D-03
        HTT( 9, 4, 2) =  0.21600D-02
        HTT( 9, 5, 2) =  0.17400D-02
        HTT( 9, 6, 2) =  0.14100D-02
        HTT( 9, 7, 2) = -0.96000D-03
        HTT( 9, 8, 2) = -0.12900D-02
        HTT( 9, 9, 2) =  0.24600D-02
        HTT(10, 1, 2) = -0.15600D-02
        HTT(10, 2, 2) =  0.23000D-03
        HTT(10, 3, 2) =  0.34000D-03
        HTT(10, 4, 2) =  0.21300D-02
        HTT(10, 5, 2) =  0.17600D-02
        HTT(10, 6, 2) =  0.79000D-03
        HTT(10, 7, 2) = -0.40000D-04
        HTT(10, 8, 2) = -0.10100D-02
        HTT(10, 9, 2) =  0.57000D-03
        HTT(10,10, 2) =  0.39100D-02
        HTT(11, 1, 2) = -0.10000D-03
        HTT(11, 2, 2) = -0.19700D-02
        HTT(11, 3, 2) =  0.19800D-02
        HTT(11, 4, 2) =  0.11000D-03
        HTT(11, 5, 2) =  0.10000D-04
        HTT(11, 6, 2) =  0.11900D-02
        HTT(11, 7, 2) = -0.19700D-02
        HTT(11, 8, 2) = -0.79000D-03
        HTT(11, 9, 2) =  0.68000D-03
        HTT(11,10, 2) = -0.12900D-02
        HTT(11,11, 2) =  0.26200D-02
        HTT(12, 1, 2) =  0.72000D-03
        HTT(12, 2, 2) = -0.10500D-02
        HTT(12, 3, 2) =  0.18600D-02
        HTT(12, 4, 2) =  0.90000D-04
        HTT(12, 5, 2) =  0.30000D-04
        HTT(12, 6, 2) =  0.57000D-03
        HTT(12, 7, 2) = -0.10500D-02
        HTT(12, 8, 2) = -0.51000D-03
        HTT(12, 9, 2) = -0.12000D-02
        HTT(12,10, 2) =  0.20500D-02
        HTT(12,11, 2) =  0.64000D-03
        HTT(12,12, 2) =  0.39000D-02
        HTT(13, 1, 2) = -0.10900D-02
        HTT(13, 2, 2) = -0.19400D-02
        HTT(13, 3, 2) =  0.59000D-03
        HTT(13, 4, 2) = -0.26000D-03
        HTT(13, 5, 2) =  0.15200D-02
        HTT(13, 6, 2) =  0.15100D-02
        HTT(13, 7, 2) = -0.70000D-04
        HTT(13, 8, 2) = -0.90000D-04
        HTT(13, 9, 2) = -0.21100D-02
        HTT(13,10, 2) =  0.12000D-03
        HTT(13,11, 2) = -0.12300D-02
        HTT(13,12, 2) =  0.10000D-02
        HTT(13,13, 2) =  0.50300D-02
        HTT(14, 1, 2) =  0.70000D-03
        HTT(14, 2, 2) =  0.53000D-03
        HTT(14, 3, 2) = -0.31000D-03
        HTT(14, 4, 2) = -0.47000D-03
        HTT(14, 5, 2) = -0.22000D-03
        HTT(14, 6, 2) = -0.80000D-03
        HTT(14, 7, 2) =  0.74000D-03
        HTT(14, 8, 2) =  0.16000D-03
        HTT(14, 9, 2) = -0.19000D-03
        HTT(14,10, 2) =  0.13500D-02
        HTT(14,11, 2) = -0.10000D-04
        HTT(14,12, 2) =  0.15300D-02
        HTT(14,13, 2) = -0.81000D-03
        HTT(14,14, 2) =  0.39200D-02
        HTT(15, 1, 2) = -0.19100D-02
        HTT(15, 2, 2) = -0.28400D-02
        HTT(15, 3, 2) =  0.70000D-03
        HTT(15, 4, 2) = -0.23000D-03
        HTT(15, 5, 2) =  0.15000D-02
        HTT(15, 6, 2) =  0.21100D-02
        HTT(15, 7, 2) = -0.98000D-03
        HTT(15, 8, 2) = -0.37000D-03
        HTT(15, 9, 2) = -0.24000D-03
        HTT(15,10, 2) = -0.31800D-02
        HTT(15,11, 2) =  0.73000D-03
        HTT(15,12, 2) = -0.22200D-02
        HTT(15,13, 2) =  0.28200D-02
        HTT(15,14, 2) = -0.23300D-02
        HTT(15,15, 2) =  0.57300D-02
        HTT(16, 1, 2) = -0.12000D-03
        HTT(16, 2, 2) = -0.37000D-03
        HTT(16, 3, 2) = -0.19000D-03
        HTT(16, 4, 2) = -0.45000D-03
        HTT(16, 5, 2) = -0.25000D-03
        HTT(16, 6, 2) = -0.19000D-03
        HTT(16, 7, 2) = -0.17000D-03
        HTT(16, 8, 2) = -0.12000D-03
        HTT(16, 9, 2) =  0.16800D-02
        HTT(16,10, 2) = -0.19500D-02
        HTT(16,11, 2) =  0.19400D-02
        HTT(16,12, 2) = -0.16900D-02
        HTT(16,13, 2) = -0.30200D-02
        HTT(16,14, 2) =  0.24000D-02
        HTT(16,15, 2) =  0.58000D-03
        HTT(16,16, 2) =  0.60000D-02
        HTT(17, 1, 2) =  0.12100D-02
        HTT(17, 2, 2) =  0.15200D-02
        HTT(17, 3, 2) = -0.31000D-03
        HTT(17, 4, 2) =  0.00000D+00
        HTT(17, 5, 2) = -0.86000D-03
        HTT(17, 6, 2) = -0.18800D-02
        HTT(17, 7, 2) =  0.58000D-03
        HTT(17, 8, 2) = -0.44000D-03
        HTT(17, 9, 2) =  0.14500D-02
        HTT(17,10, 2) =  0.90000D-04
        HTT(17,11, 2) =  0.11300D-02
        HTT(17,12, 2) = -0.23000D-03
        HTT(17,13, 2) = -0.46000D-02
        HTT(17,14, 2) =  0.13000D-02
        HTT(17,15, 2) = -0.32400D-02
        HTT(17,16, 2) =  0.26500D-02
        HTT(17,17, 2) =  0.49600D-02
        HTT(18, 1, 2) =  0.14500D-02
        HTT(18, 2, 2) =  0.22200D-02
        HTT(18, 3, 2) = -0.10600D-02
        HTT(18, 4, 2) = -0.29000D-03
        HTT(18, 5, 2) = -0.17100D-02
        HTT(18, 6, 2) = -0.27100D-02
        HTT(18, 7, 2) =  0.67000D-03
        HTT(18, 8, 2) = -0.32000D-03
        HTT(18, 9, 2) =  0.10000D-02
        HTT(18,10, 2) =  0.82000D-03
        HTT(18,11, 2) =  0.19000D-03
        HTT(18,12, 2) =  0.10000D-04
        HTT(18,13, 2) = -0.32400D-02
        HTT(18,14, 2) =  0.19000D-02
        HTT(18,15, 2) = -0.30600D-02
        HTT(18,16, 2) =  0.20800D-02
        HTT(18,17, 2) =  0.29700D-02
        HTT(18,18, 2) =  0.56300D-02
        HTT(19, 1, 2) = -0.59000D-03
        HTT(19, 2, 2) = -0.96000D-03
        HTT(19, 3, 2) =  0.59000D-03
        HTT(19, 4, 2) =  0.21000D-03
        HTT(19, 5, 2) =  0.89000D-03
        HTT(19, 6, 2) =  0.43000D-03
        HTT(19, 7, 2) = -0.23000D-03
        HTT(19, 8, 2) = -0.69000D-03
        HTT(19, 9, 2) = -0.48000D-03
        HTT(19,10, 2) = -0.11500D-02
        HTT(19,11, 2) = -0.90000D-04
        HTT(19,12, 2) = -0.76000D-03
        HTT(19,13, 2) =  0.12600D-02
        HTT(19,14, 2) = -0.34500D-02
        HTT(19,15, 2) =  0.19200D-02
        HTT(19,16, 2) = -0.27800D-02
        HTT(19,17, 2) = -0.95000D-03
        HTT(19,18, 2) = -0.21800D-02
        HTT(19,19, 2) =  0.37700D-02
        HTT(20, 1, 2) = -0.35000D-03
        HTT(20, 2, 2) = -0.26000D-03
        HTT(20, 3, 2) = -0.16000D-03
        HTT(20, 4, 2) = -0.80000D-04
        HTT(20, 5, 2) =  0.40000D-04
        HTT(20, 6, 2) = -0.40000D-03
        HTT(20, 7, 2) = -0.14000D-03
        HTT(20, 8, 2) = -0.57000D-03
        HTT(20, 9, 2) = -0.93000D-03
        HTT(20,10, 2) = -0.42000D-03
        HTT(20,11, 2) = -0.10200D-02
        HTT(20,12, 2) = -0.51000D-03
        HTT(20,13, 2) =  0.26100D-02
        HTT(20,14, 2) = -0.28400D-02
        HTT(20,15, 2) =  0.21100D-02
        HTT(20,16, 2) = -0.33500D-02
        HTT(20,17, 2) = -0.29400D-02
        HTT(20,18, 2) =  0.47000D-03
        HTT(20,19, 2) =  0.25400D-02
        HTT(20,20, 2) =  0.59400D-02
        HTT(21, 1, 2) =  0.21700D-02
        HTT(21, 2, 2) =  0.15000D-02
        HTT(21, 3, 2) = -0.99000D-03
        HTT(21, 4, 2) = -0.16600D-02
        HTT(21, 5, 2) = -0.30200D-02
        HTT(21, 6, 2) = -0.69000D-03
        HTT(21, 7, 2) = -0.20000D-04
        HTT(21, 8, 2) =  0.23100D-02
        HTT(21, 9, 2) = -0.11900D-02
        HTT(21,10, 2) = -0.98000D-03
        HTT(21,11, 2) = -0.50000D-03
        HTT(21,12, 2) = -0.28000D-03
        HTT(21,13, 2) =  0.12600D-02
        HTT(21,14, 2) = -0.76000D-03
        HTT(21,15, 2) =  0.10500D-02
        HTT(21,16, 2) = -0.97000D-03
        HTT(21,17, 2) = -0.21400D-02
        HTT(21,18, 2) = -0.48000D-03
        HTT(21,19, 2) = -0.11000D-03
        HTT(21,20, 2) =  0.15500D-02
        HTT(21,21, 2) =  0.28900D-02
        HTT(22, 1, 2) =  0.27000D-03
        HTT(22, 2, 2) =  0.11300D-02
        HTT(22, 3, 2) = -0.15900D-02
        HTT(22, 4, 2) = -0.74000D-03
        HTT(22, 5, 2) = -0.60000D-04
        HTT(22, 6, 2) = -0.26200D-02
        HTT(22, 7, 2) =  0.17100D-02
        HTT(22, 8, 2) = -0.85000D-03
        HTT(22, 9, 2) = -0.85000D-03
        HTT(22,10, 2) =  0.40000D-04
        HTT(22,11, 2) = -0.17300D-02
        HTT(22,12, 2) = -0.85000D-03
        HTT(22,13, 2) =  0.12800D-02
        HTT(22,14, 2) = -0.15000D-03
        HTT(22,15, 2) =  0.40000D-03
        HTT(22,16, 2) = -0.10300D-02
        HTT(22,17, 2) = -0.10700D-02
        HTT(22,18, 2) =  0.56000D-03
        HTT(22,19, 2) =  0.37000D-03
        HTT(22,20, 2) =  0.20000D-02
        HTT(22,21, 2) =  0.44000D-03
        HTT(22,22, 2) =  0.31200D-02
        HTT(23, 1, 2) =  0.19200D-02
        HTT(23, 2, 2) =  0.79000D-03
        HTT(23, 3, 2) = -0.23000D-03
        HTT(23, 4, 2) = -0.13700D-02
        HTT(23, 5, 2) = -0.21700D-02
        HTT(23, 6, 2) =  0.14000D-03
        HTT(23, 7, 2) = -0.12000D-03
        HTT(23, 8, 2) =  0.21900D-02
        HTT(23, 9, 2) = -0.74000D-03
        HTT(23,10, 2) = -0.17100D-02
        HTT(23,11, 2) =  0.45000D-03
        HTT(23,12, 2) = -0.53000D-03
        HTT(23,13, 2) = -0.10000D-03
        HTT(23,14, 2) = -0.13700D-02
        HTT(23,15, 2) =  0.87000D-03
        HTT(23,16, 2) = -0.40000D-03
        HTT(23,17, 2) = -0.15000D-03
        HTT(23,18, 2) = -0.31500D-02
        HTT(23,19, 2) =  0.11200D-02
        HTT(23,20, 2) = -0.18800D-02
        HTT(23,21, 2) =  0.12200D-02
        HTT(23,22, 2) = -0.12000D-02
        HTT(23,23, 2) =  0.42400D-02
        HTT(24, 1, 2) =  0.30000D-04
        HTT(24, 2, 2) =  0.42000D-03
        HTT(24, 3, 2) = -0.84000D-03
        HTT(24, 4, 2) = -0.45000D-03
        HTT(24, 5, 2) =  0.79000D-03
        HTT(24, 6, 2) = -0.18000D-02
        HTT(24, 7, 2) =  0.16200D-02
        HTT(24, 8, 2) = -0.97000D-03
        HTT(24, 9, 2) = -0.39000D-03
        HTT(24,10, 2) = -0.70000D-03
        HTT(24,11, 2) = -0.79000D-03
        HTT(24,12, 2) = -0.11000D-02
        HTT(24,13, 2) = -0.80000D-04
        HTT(24,14, 2) = -0.76000D-03
        HTT(24,15, 2) =  0.22000D-03
        HTT(24,16, 2) = -0.45000D-03
        HTT(24,17, 2) =  0.92000D-03
        HTT(24,18, 2) = -0.21000D-02
        HTT(24,19, 2) =  0.16000D-02
        HTT(24,20, 2) = -0.14300D-02
        HTT(24,21, 2) = -0.12300D-02
        HTT(24,22, 2) =  0.14800D-02
        HTT(24,23, 2) =  0.18100D-02
        HTT(24,24, 2) =  0.45300D-02
C Hessian elements for anchor point 2
        HBS( 1, 1, 2) = -0.31760D-01
        HBS( 1, 2, 2) = -0.29360D-01
        HBS( 1, 3, 2) = -0.31820D-01
        HBS( 1, 4, 2) = -0.95900D-02
        HBS( 1, 5, 2) =  0.52300D-02
        HBS( 1, 6, 2) =  0.39930D-01
        HBS( 1, 7, 2) = -0.45000D-03
        HBS( 1, 8, 2) =  0.40450D-01
        HBS( 1, 9, 2) = -0.18000D-03
        HBS( 1,10, 2) = -0.96900D-02
        HBS( 1,11, 2) = -0.42000D-03
        HBS( 1,12, 2) =  0.47400D-02
        HBS( 2, 1, 2) =  0.49930D-01
        HBS( 2, 2, 2) = -0.13930D-01
        HBS( 2, 3, 2) =  0.99600D-02
        HBS( 2, 4, 2) =  0.18000D-02
        HBS( 2, 5, 2) = -0.12200D-02
        HBS( 2, 6, 2) = -0.27810D-01
        HBS( 2, 7, 2) =  0.11800D-02
        HBS( 2, 8, 2) = -0.14650D-01
        HBS( 2, 9, 2) =  0.00000D+00
        HBS( 2,10, 2) =  0.41300D-02
        HBS( 2,11, 2) = -0.69000D-03
        HBS( 2,12, 2) = -0.89000D-03
        HBS( 3, 1, 2) = -0.18180D-01
        HBS( 3, 2, 2) =  0.43280D-01
        HBS( 3, 3, 2) =  0.21860D-01
        HBS( 3, 4, 2) =  0.77900D-02
        HBS( 3, 5, 2) = -0.40100D-02
        HBS( 3, 6, 2) = -0.12120D-01
        HBS( 3, 7, 2) = -0.73000D-03
        HBS( 3, 8, 2) = -0.25800D-01
        HBS( 3, 9, 2) =  0.18000D-03
        HBS( 3,10, 2) =  0.55600D-02
        HBS( 3,11, 2) =  0.11100D-02
        HBS( 3,12, 2) = -0.38500D-02
        HBS( 4, 1, 2) = -0.20820D-01
        HBS( 4, 2, 2) = -0.12100D-01
        HBS( 4, 3, 2) =  0.15180D-01
        HBS( 4, 4, 2) = -0.29840D-01
        HBS( 4, 5, 2) = -0.97300D-02
        HBS( 4, 6, 2) = -0.57000D-02
        HBS( 4, 7, 2) =  0.49500D-02
        HBS( 4, 8, 2) =  0.31150D-01
        HBS( 4, 9, 2) = -0.44000D-03
        HBS( 4,10, 2) =  0.40270D-01
        HBS( 4,11, 2) = -0.59000D-03
        HBS( 4,12, 2) =  0.29000D-03
        HBS( 5, 1, 2) =  0.26090D-01
        HBS( 5, 2, 2) =  0.58900D-02
        HBS( 5, 3, 2) = -0.14520D-01
        HBS( 5, 4, 2) =  0.25400D-02
        HBS( 5, 5, 2) =  0.69400D-02
        HBS( 5, 6, 2) =  0.30900D-02
        HBS( 5, 7, 2) = -0.21700D-02
        HBS( 5, 8, 2) = -0.15190D-01
        HBS( 5, 9, 2) = -0.39000D-03
        HBS( 5,10, 2) = -0.22140D-01
        HBS( 5,11, 2) =  0.28000D-03
        HBS( 5,12, 2) =  0.90000D-04
        HBS( 6, 1, 2) = -0.52700D-02
        HBS( 6, 2, 2) =  0.62100D-02
        HBS( 6, 3, 2) = -0.66000D-03
        HBS( 6, 4, 2) =  0.27300D-01
        HBS( 6, 5, 2) =  0.27900D-02
        HBS( 6, 6, 2) =  0.26100D-02
        HBS( 6, 7, 2) = -0.27800D-02
        HBS( 6, 8, 2) = -0.15960D-01
        HBS( 6, 9, 2) =  0.83000D-03
        HBS( 6,10, 2) = -0.18130D-01
        HBS( 6,11, 2) =  0.31000D-03
        HBS( 6,12, 2) = -0.38000D-03
        HBS( 7, 1, 2) = -0.90800D-02
        HBS( 7, 2, 2) =  0.35870D-01
        HBS( 7, 3, 2) = -0.10000D-04
        HBS( 7, 4, 2) = -0.22920D-01
        HBS( 7, 5, 2) =  0.61900D-02
        HBS( 7, 6, 2) = -0.28970D-01
        HBS( 7, 7, 2) = -0.87500D-02
        HBS( 7, 8, 2) = -0.92900D-02
        HBS( 7, 9, 2) =  0.49700D-02
        HBS( 7,10, 2) =  0.32780D-01
        HBS( 7,11, 2) =  0.60000D-04
        HBS( 7,12, 2) = -0.10700D-02
        HBS( 8, 1, 2) =  0.36900D-02
        HBS( 8, 2, 2) = -0.18590D-01
        HBS( 8, 3, 2) =  0.13800D-02
        HBS( 8, 4, 2) =  0.25350D-01
        HBS( 8, 5, 2) = -0.35000D-02
        HBS( 8, 6, 2) = -0.42000D-03
        HBS( 8, 7, 2) =  0.45000D-02
        HBS( 8, 8, 2) =  0.50400D-02
        HBS( 8, 9, 2) = -0.20100D-02
        HBS( 8,10, 2) = -0.14680D-01
        HBS( 8,11, 2) = -0.62000D-03
        HBS( 8,12, 2) =  0.55000D-03
        HBS( 9, 1, 2) =  0.53900D-02
        HBS( 9, 2, 2) = -0.17270D-01
        HBS( 9, 3, 2) = -0.13700D-02
        HBS( 9, 4, 2) = -0.24400D-02
        HBS( 9, 5, 2) = -0.26800D-02
        HBS( 9, 6, 2) =  0.29390D-01
        HBS( 9, 7, 2) =  0.42400D-02
        HBS( 9, 8, 2) =  0.42500D-02
        HBS( 9, 9, 2) = -0.29600D-02
        HBS( 9,10, 2) = -0.18100D-01
        HBS( 9,11, 2) =  0.55000D-03
        HBS( 9,12, 2) =  0.51000D-03
        HBS(10, 1, 2) =  0.36130D-01
        HBS(10, 2, 2) =  0.38360D-01
        HBS(10, 3, 2) =  0.36300D-02
        HBS(10, 4, 2) = -0.11890D-01
        HBS(10, 5, 2) = -0.10500D-02
        HBS(10, 6, 2) = -0.26840D-01
        HBS(10, 7, 2) =  0.48100D-02
        HBS(10, 8, 2) = -0.27100D-01
        HBS(10, 9, 2) = -0.88200D-02
        HBS(10,10, 2) = -0.12140D-01
        HBS(10,11, 2) =  0.47600D-02
        HBS(10,12, 2) = -0.77000D-03
        HBS(11, 1, 2) = -0.20190D-01
        HBS(11, 2, 2) = -0.17180D-01
        HBS(11, 3, 2) = -0.21300D-02
        HBS(11, 4, 2) =  0.57400D-02
        HBS(11, 5, 2) =  0.12900D-02
        HBS(11, 6, 2) =  0.27000D-01
        HBS(11, 7, 2) = -0.28400D-02
        HBS(11, 8, 2) =  0.47000D-03
        HBS(11, 9, 2) =  0.44300D-02
        HBS(11,10, 2) =  0.60700D-02
        HBS(11,11, 2) = -0.19800D-02
        HBS(11,12, 2) = -0.18000D-03
        HBS(12, 1, 2) = -0.15940D-01
        HBS(12, 2, 2) = -0.21180D-01
        HBS(12, 3, 2) = -0.14900D-02
        HBS(12, 4, 2) =  0.61500D-02
        HBS(12, 5, 2) = -0.24000D-03
        HBS(12, 6, 2) = -0.16000D-03
        HBS(12, 7, 2) = -0.19600D-02
        HBS(12, 8, 2) =  0.26630D-01
        HBS(12, 9, 2) =  0.43900D-02
        HBS(12,10, 2) =  0.60700D-02
        HBS(12,11, 2) = -0.27800D-02
        HBS(12,12, 2) =  0.94000D-03
        HBS(13, 1, 2) =  0.35720D-01
        HBS(13, 2, 2) = -0.96900D-02
        HBS(13, 3, 2) =  0.29800D-02
        HBS(13, 4, 2) =  0.32310D-01
        HBS(13, 5, 2) = -0.94000D-03
        HBS(13, 6, 2) = -0.10140D-01
        HBS(13, 7, 2) =  0.70000D-04
        HBS(13, 8, 2) = -0.29530D-01
        HBS(13, 9, 2) =  0.50200D-02
        HBS(13,10, 2) = -0.23420D-01
        HBS(13,11, 2) = -0.86600D-02
        HBS(13,12, 2) =  0.53400D-02
        HBS(14, 1, 2) = -0.17210D-01
        HBS(14, 2, 2) =  0.56100D-02
        HBS(14, 3, 2) = -0.23300D-02
        HBS(14, 4, 2) = -0.17680D-01
        HBS(14, 5, 2) =  0.48000D-03
        HBS(14, 6, 2) =  0.47700D-02
        HBS(14, 7, 2) =  0.56000D-03
        HBS(14, 8, 2) =  0.29630D-01
        HBS(14, 9, 2) = -0.29700D-02
        HBS(14,10, 2) = -0.24000D-02
        HBS(14,11, 2) =  0.43000D-02
        HBS(14,12, 2) = -0.22600D-02
        HBS(15, 1, 2) = -0.18510D-01
        HBS(15, 2, 2) =  0.40800D-02
        HBS(15, 3, 2) = -0.65000D-03
        HBS(15, 4, 2) = -0.14630D-01
        HBS(15, 5, 2) =  0.46000D-03
        HBS(15, 6, 2) =  0.53700D-02
        HBS(15, 7, 2) = -0.63000D-03
        HBS(15, 8, 2) = -0.10000D-03
        HBS(15, 9, 2) = -0.20500D-02
        HBS(15,10, 2) =  0.25820D-01
        HBS(15,11, 2) =  0.43700D-02
        HBS(15,12, 2) = -0.30900D-02
        HBS(16, 1, 2) = -0.10200D-01
        HBS(16, 2, 2) = -0.23080D-01
        HBS(16, 3, 2) =  0.10040D-01
        HBS(16, 4, 2) =  0.41930D-01
        HBS(16, 5, 2) =  0.31000D-03
        HBS(16, 6, 2) =  0.31710D-01
        HBS(16, 7, 2) = -0.63000D-03
        HBS(16, 8, 2) = -0.56800D-02
        HBS(16, 9, 2) = -0.55000D-03
        HBS(16,10, 2) = -0.27810D-01
        HBS(16,11, 2) =  0.48500D-02
        HBS(16,12, 2) = -0.85400D-02
        HBS(17, 1, 2) =  0.41100D-02
        HBS(17, 2, 2) =  0.25000D-01
        HBS(17, 3, 2) = -0.74400D-02
        HBS(17, 4, 2) = -0.22980D-01
        HBS(17, 5, 2) =  0.60000D-03
        HBS(17, 6, 2) = -0.14860D-01
        HBS(17, 7, 2) =  0.32000D-03
        HBS(17, 8, 2) =  0.33600D-02
        HBS(17, 9, 2) = -0.32000D-03
        HBS(17,10, 2) =  0.20400D-02
        HBS(17,11, 2) = -0.21000D-02
        HBS(17,12, 2) =  0.68800D-02
        HBS(18, 1, 2) =  0.61000D-02
        HBS(18, 2, 2) = -0.19200D-02
        HBS(18, 3, 2) = -0.26000D-02
        HBS(18, 4, 2) = -0.18950D-01
        HBS(18, 5, 2) = -0.92000D-03
        HBS(18, 6, 2) = -0.16850D-01
        HBS(18, 7, 2) =  0.31000D-03
        HBS(18, 8, 2) =  0.23200D-02
        HBS(18, 9, 2) =  0.87000D-03
        HBS(18,10, 2) =  0.25770D-01
        HBS(18,11, 2) = -0.27500D-02
        HBS(18,12, 2) =  0.16600D-02
C Hessian elements for anchor point 3
        HSS( 1, 1, 3) =  0.34759D+00
        HSS( 2, 1, 3) =  0.17810D-01
        HSS( 2, 2, 3) =  0.35623D+00
        HSS( 3, 1, 3) =  0.29010D-01
        HSS( 3, 2, 3) =  0.32710D-01
        HSS( 3, 3, 3) =  0.46124D+00
        HSS( 4, 1, 3) = -0.19600D-02
        HSS( 4, 2, 3) =  0.14870D-01
        HSS( 4, 3, 3) = -0.89100D-02
        HSS( 4, 4, 3) =  0.38581D+00
        HSS( 5, 1, 3) =  0.38800D-02
        HSS( 5, 2, 3) = -0.77000D-03
        HSS( 5, 3, 3) =  0.71700D-02
        HSS( 5, 4, 3) =  0.36500D-02
        HSS( 5, 5, 3) =  0.35807D+00
        HSS( 6, 1, 3) =  0.33890D-01
        HSS( 6, 2, 3) =  0.44570D-01
        HSS( 6, 3, 3) = -0.43100D-02
        HSS( 6, 4, 3) = -0.36100D-02
        HSS( 6, 5, 3) =  0.55000D-03
        HSS( 6, 6, 3) =  0.35931D+00
        HSS( 7, 1, 3) = -0.54000D-03
        HSS( 7, 2, 3) = -0.63000D-03
        HSS( 7, 3, 3) =  0.57000D-03
        HSS( 7, 4, 3) =  0.39700D-02
        HSS( 7, 5, 3) =  0.39000D-03
        HSS( 7, 6, 3) =  0.39800D-02
        HSS( 7, 7, 3) =  0.37091D+00
        HSS( 8, 1, 3) =  0.42030D-01
        HSS( 8, 2, 3) =  0.35060D-01
        HSS( 8, 3, 3) = -0.52800D-02
        HSS( 8, 4, 3) =  0.24280D-01
        HSS( 8, 5, 3) = -0.11100D-02
        HSS( 8, 6, 3) =  0.17880D-01
        HSS( 8, 7, 3) = -0.10400D-02
        HSS( 8, 8, 3) =  0.35888D+00
        HSS( 9, 1, 3) = -0.53000D-03
        HSS( 9, 2, 3) = -0.57000D-03
        HSS( 9, 3, 3) =  0.26000D-03
        HSS( 9, 4, 3) = -0.39000D-03
        HSS( 9, 5, 3) =  0.19000D-03
        HSS( 9, 6, 3) =  0.40100D-02
        HSS( 9, 7, 3) =  0.49000D-03
        HSS( 9, 8, 3) =  0.39600D-02
        HSS( 9, 9, 3) =  0.36356D+00
        HSS(10, 1, 3) =  0.19520D-01
        HSS(10, 2, 3) = -0.69300D-02
        HSS(10, 3, 3) = -0.48800D-02
        HSS(10, 4, 3) =  0.64010D-01
        HSS(10, 5, 3) = -0.43000D-03
        HSS(10, 6, 3) =  0.21230D-01
        HSS(10, 7, 3) = -0.10900D-02
        HSS(10, 8, 3) = -0.77000D-03
        HSS(10, 9, 3) = -0.19000D-03
        HSS(10,10, 3) =  0.38116D+00
        HSS(11, 1, 3) = -0.69000D-03
        HSS(11, 2, 3) = -0.50000D-03
        HSS(11, 3, 3) =  0.52000D-03
        HSS(11, 4, 3) = -0.11600D-02
        HSS(11, 5, 3) =  0.80000D-04
        HSS(11, 6, 3) = -0.10000D-02
        HSS(11, 7, 3) =  0.22000D-03
        HSS(11, 8, 3) =  0.39700D-02
        HSS(11, 9, 3) =  0.48000D-03
        HSS(11,10, 3) =  0.39800D-02
        HSS(11,11, 3) =  0.37152D+00
        HSS(12, 1, 3) = -0.16500D-02
        HSS(12, 2, 3) =  0.39400D-02
        HSS(12, 3, 3) = -0.19000D-03
        HSS(12, 4, 3) = -0.10000D-03
        HSS(12, 5, 3) =  0.35000D-03
        HSS(12, 6, 3) = -0.10400D-02
        HSS(12, 7, 3) =  0.30000D-04
        HSS(12, 8, 3) =  0.57000D-03
        HSS(12, 9, 3) =  0.12000D-03
        HSS(12,10, 3) =  0.23600D-02
        HSS(12,11, 3) =  0.40000D-03
        HSS(12,12, 3) =  0.37254D+00
C Hessian elements for anchor point 3
        HBB( 1, 1, 3) =  0.82450D-01
        HBB( 2, 1, 3) = -0.46090D-01
        HBB( 2, 2, 3) =  0.14487D+00
        HBB( 3, 1, 3) = -0.36360D-01
        HBB( 3, 2, 3) = -0.98780D-01
        HBB( 3, 3, 3) =  0.13515D+00
        HBB( 4, 1, 3) = -0.36000D-01
        HBB( 4, 2, 3) =  0.23310D-01
        HBB( 4, 3, 3) =  0.12690D-01
        HBB( 4, 4, 3) =  0.66170D-01
        HBB( 5, 1, 3) =  0.20240D-01
        HBB( 5, 2, 3) = -0.12480D-01
        HBB( 5, 3, 3) = -0.77700D-02
        HBB( 5, 4, 3) = -0.33030D-01
        HBB( 5, 5, 3) =  0.71950D-01
        HBB( 6, 1, 3) =  0.15760D-01
        HBB( 6, 2, 3) = -0.10840D-01
        HBB( 6, 3, 3) = -0.49200D-02
        HBB( 6, 4, 3) = -0.33140D-01
        HBB( 6, 5, 3) = -0.38920D-01
        HBB( 6, 6, 3) =  0.72060D-01
        HBB( 7, 1, 3) = -0.71500D-02
        HBB( 7, 2, 3) = -0.39000D-02
        HBB( 7, 3, 3) =  0.11050D-01
        HBB( 7, 4, 3) = -0.25740D-01
        HBB( 7, 5, 3) =  0.93100D-02
        HBB( 7, 6, 3) =  0.16430D-01
        HBB( 7, 7, 3) =  0.64420D-01
        HBB( 8, 1, 3) = -0.22800D-02
        HBB( 8, 2, 3) =  0.74300D-02
        HBB( 8, 3, 3) = -0.51500D-02
        HBB( 8, 4, 3) =  0.15690D-01
        HBB( 8, 5, 3) = -0.49500D-02
        HBB( 8, 6, 3) = -0.10740D-01
        HBB( 8, 7, 3) = -0.31490D-01
        HBB( 8, 8, 3) =  0.72950D-01
        HBB( 9, 1, 3) =  0.94300D-02
        HBB( 9, 2, 3) = -0.35300D-02
        HBB( 9, 3, 3) = -0.59100D-02
        HBB( 9, 4, 3) =  0.10050D-01
        HBB( 9, 5, 3) = -0.43500D-02
        HBB( 9, 6, 3) = -0.56900D-02
        HBB( 9, 7, 3) = -0.32930D-01
        HBB( 9, 8, 3) = -0.41460D-01
        HBB( 9, 9, 3) =  0.74390D-01
        HBB(10, 1, 3) =  0.45500D-02
        HBB(10, 2, 3) = -0.20200D-02
        HBB(10, 3, 3) = -0.25300D-02
        HBB(10, 4, 3) = -0.64100D-02
        HBB(10, 5, 3) =  0.91900D-02
        HBB(10, 6, 3) = -0.27700D-02
        HBB(10, 7, 3) = -0.33410D-01
        HBB(10, 8, 3) =  0.13410D-01
        HBB(10, 9, 3) =  0.20000D-01
        HBB(10,10, 3) =  0.75440D-01
        HBB(11, 1, 3) = -0.22100D-02
        HBB(11, 2, 3) =  0.12300D-02
        HBB(11, 3, 3) =  0.98000D-03
        HBB(11, 4, 3) = -0.23100D-02
        HBB(11, 5, 3) = -0.32200D-02
        HBB(11, 6, 3) =  0.55300D-02
        HBB(11, 7, 3) =  0.20020D-01
        HBB(11, 8, 3) = -0.74100D-02
        HBB(11, 9, 3) = -0.12610D-01
        HBB(11,10, 3) = -0.37790D-01
        HBB(11,11, 3) =  0.78360D-01
        HBB(12, 1, 3) = -0.23400D-02
        HBB(12, 2, 3) =  0.79000D-03
        HBB(12, 3, 3) =  0.15500D-02
        HBB(12, 4, 3) =  0.87200D-02
        HBB(12, 5, 3) = -0.59700D-02
        HBB(12, 6, 3) = -0.27500D-02
        HBB(12, 7, 3) =  0.13390D-01
        HBB(12, 8, 3) = -0.60000D-02
        HBB(12, 9, 3) = -0.73900D-02
        HBB(12,10, 3) = -0.37640D-01
        HBB(12,11, 3) = -0.40570D-01
        HBB(12,12, 3) =  0.78220D-01
        HBB(13, 1, 3) = -0.67600D-02
        HBB(13, 2, 3) =  0.10710D-01
        HBB(13, 3, 3) = -0.39500D-02
        HBB(13, 4, 3) = -0.11500D-02
        HBB(13, 5, 3) =  0.12900D-02
        HBB(13, 6, 3) = -0.14000D-03
        HBB(13, 7, 3) =  0.29100D-02
        HBB(13, 8, 3) =  0.41300D-02
        HBB(13, 9, 3) = -0.70300D-02
        HBB(13,10, 3) = -0.33760D-01
        HBB(13,11, 3) =  0.13510D-01
        HBB(13,12, 3) =  0.20240D-01
        HBB(13,13, 3) =  0.65300D-01
        HBB(14, 1, 3) =  0.91200D-02
        HBB(14, 2, 3) = -0.60500D-02
        HBB(14, 3, 3) = -0.30700D-02
        HBB(14, 4, 3) =  0.66000D-03
        HBB(14, 5, 3) = -0.85000D-03
        HBB(14, 6, 3) =  0.19000D-03
        HBB(14, 7, 3) = -0.71100D-02
        HBB(14, 8, 3) = -0.51000D-03
        HBB(14, 9, 3) =  0.76200D-02
        HBB(14,10, 3) =  0.20030D-01
        HBB(14,11, 3) = -0.73800D-02
        HBB(14,12, 3) = -0.12650D-01
        HBB(14,13, 3) = -0.33180D-01
        HBB(14,14, 3) =  0.74350D-01
        HBB(15, 1, 3) = -0.23600D-02
        HBB(15, 2, 3) = -0.46500D-02
        HBB(15, 3, 3) =  0.70200D-02
        HBB(15, 4, 3) =  0.49000D-03
        HBB(15, 5, 3) = -0.44000D-03
        HBB(15, 6, 3) = -0.50000D-04
        HBB(15, 7, 3) =  0.42000D-02
        HBB(15, 8, 3) = -0.36200D-02
        HBB(15, 9, 3) = -0.59000D-03
        HBB(15,10, 3) =  0.13730D-01
        HBB(15,11, 3) = -0.61300D-02
        HBB(15,12, 3) = -0.76000D-02
        HBB(15,13, 3) = -0.32120D-01
        HBB(15,14, 3) = -0.41170D-01
        HBB(15,15, 3) =  0.73290D-01
        HBB(16, 1, 3) = -0.37090D-01
        HBB(16, 2, 3) =  0.17990D-01
        HBB(16, 3, 3) =  0.19100D-01
        HBB(16, 4, 3) =  0.31400D-02
        HBB(16, 5, 3) = -0.70000D-02
        HBB(16, 6, 3) =  0.38600D-02
        HBB(16, 7, 3) = -0.10300D-02
        HBB(16, 8, 3) =  0.54000D-03
        HBB(16, 9, 3) =  0.48000D-03
        HBB(16,10, 3) = -0.64000D-02
        HBB(16,11, 3) =  0.87800D-02
        HBB(16,12, 3) = -0.23700D-02
        HBB(16,13, 3) = -0.26530D-01
        HBB(16,14, 3) =  0.10480D-01
        HBB(16,15, 3) =  0.16060D-01
        HBB(16,16, 3) =  0.67910D-01
        HBB(17, 1, 3) =  0.22610D-01
        HBB(17, 2, 3) = -0.90200D-02
        HBB(17, 3, 3) = -0.13590D-01
        HBB(17, 4, 3) = -0.76000D-02
        HBB(17, 5, 3) =  0.77900D-02
        HBB(17, 6, 3) = -0.19000D-03
        HBB(17, 7, 3) =  0.87000D-03
        HBB(17, 8, 3) = -0.25000D-03
        HBB(17, 9, 3) = -0.62000D-03
        HBB(17,10, 3) =  0.91600D-02
        HBB(17,11, 3) = -0.60700D-02
        HBB(17,12, 3) = -0.31000D-02
        HBB(17,13, 3) =  0.96800D-02
        HBB(17,14, 3) = -0.46000D-02
        HBB(17,15, 3) = -0.50800D-02
        HBB(17,16, 3) = -0.34720D-01
        HBB(17,17, 3) =  0.72390D-01
        HBB(18, 1, 3) =  0.14470D-01
        HBB(18, 2, 3) = -0.89700D-02
        HBB(18, 3, 3) = -0.55000D-02
        HBB(18, 4, 3) =  0.44600D-02
        HBB(18, 5, 3) = -0.79000D-03
        HBB(18, 6, 3) = -0.36700D-02
        HBB(18, 7, 3) =  0.16000D-03
        HBB(18, 8, 3) = -0.29000D-03
        HBB(18, 9, 3) =  0.13000D-03
        HBB(18,10, 3) = -0.27600D-02
        HBB(18,11, 3) = -0.27100D-02
        HBB(18,12, 3) =  0.54700D-02
        HBB(18,13, 3) =  0.16860D-01
        HBB(18,14, 3) = -0.58800D-02
        HBB(18,15, 3) = -0.10980D-01
        HBB(18,16, 3) = -0.33190D-01
        HBB(18,17, 3) = -0.37670D-01
        HBB(18,18, 3) =  0.70860D-01
C Hessian elements for anchor point 3
        HTT( 1, 1, 3) =  0.60400D-02
        HTT( 2, 1, 3) =  0.38400D-02
        HTT( 2, 2, 3) =  0.50900D-02
        HTT( 3, 1, 3) = -0.18500D-02
        HTT( 3, 2, 3) = -0.33100D-02
        HTT( 3, 3, 3) =  0.75100D-02
        HTT( 4, 1, 3) = -0.40500D-02
        HTT( 4, 2, 3) = -0.20700D-02
        HTT( 4, 3, 3) =  0.60500D-02
        HTT( 4, 4, 3) =  0.80200D-02
        HTT( 5, 1, 3) = -0.59200D-02
        HTT( 5, 2, 3) = -0.42600D-02
        HTT( 5, 3, 3) =  0.21200D-02
        HTT( 5, 4, 3) =  0.37800D-02
        HTT( 5, 5, 3) =  0.65600D-02
        HTT( 6, 1, 3) = -0.41200D-02
        HTT( 6, 2, 3) = -0.39000D-02
        HTT( 6, 3, 3) =  0.26900D-02
        HTT( 6, 4, 3) =  0.29100D-02
        HTT( 6, 5, 3) =  0.37400D-02
        HTT( 6, 6, 3) =  0.55900D-02
        HTT( 7, 1, 3) =  0.15700D-02
        HTT( 7, 2, 3) =  0.25400D-02
        HTT( 7, 3, 3) = -0.67700D-02
        HTT( 7, 4, 3) = -0.58000D-02
        HTT( 7, 5, 3) = -0.10700D-02
        HTT( 7, 6, 3) = -0.27200D-02
        HTT( 7, 7, 3) =  0.68600D-02
        HTT( 8, 1, 3) =  0.33800D-02
        HTT( 8, 2, 3) =  0.29000D-02
        HTT( 8, 3, 3) = -0.62000D-02
        HTT( 8, 4, 3) = -0.66800D-02
        HTT( 8, 5, 3) = -0.38900D-02
        HTT( 8, 6, 3) = -0.88000D-03
        HTT( 8, 7, 3) =  0.52100D-02
        HTT( 8, 8, 3) =  0.82200D-02
        HTT( 9, 1, 3) = -0.23900D-02
        HTT( 9, 2, 3) = -0.68000D-03
        HTT( 9, 3, 3) =  0.45000D-03
        HTT( 9, 4, 3) =  0.21600D-02
        HTT( 9, 5, 3) =  0.17400D-02
        HTT( 9, 6, 3) =  0.14100D-02
        HTT( 9, 7, 3) = -0.96000D-03
        HTT( 9, 8, 3) = -0.12900D-02
        HTT( 9, 9, 3) =  0.24600D-02
        HTT(10, 1, 3) = -0.15600D-02
        HTT(10, 2, 3) =  0.23000D-03
        HTT(10, 3, 3) =  0.34000D-03
        HTT(10, 4, 3) =  0.21300D-02
        HTT(10, 5, 3) =  0.17600D-02
        HTT(10, 6, 3) =  0.79000D-03
        HTT(10, 7, 3) = -0.40000D-04
        HTT(10, 8, 3) = -0.10100D-02
        HTT(10, 9, 3) =  0.57000D-03
        HTT(10,10, 3) =  0.39100D-02
        HTT(11, 1, 3) = -0.10000D-03
        HTT(11, 2, 3) = -0.19700D-02
        HTT(11, 3, 3) =  0.19800D-02
        HTT(11, 4, 3) =  0.11000D-03
        HTT(11, 5, 3) =  0.10000D-04
        HTT(11, 6, 3) =  0.11900D-02
        HTT(11, 7, 3) = -0.19700D-02
        HTT(11, 8, 3) = -0.79000D-03
        HTT(11, 9, 3) =  0.68000D-03
        HTT(11,10, 3) = -0.12900D-02
        HTT(11,11, 3) =  0.26200D-02
        HTT(12, 1, 3) =  0.72000D-03
        HTT(12, 2, 3) = -0.10500D-02
        HTT(12, 3, 3) =  0.18600D-02
        HTT(12, 4, 3) =  0.90000D-04
        HTT(12, 5, 3) =  0.30000D-04
        HTT(12, 6, 3) =  0.57000D-03
        HTT(12, 7, 3) = -0.10500D-02
        HTT(12, 8, 3) = -0.51000D-03
        HTT(12, 9, 3) = -0.12000D-02
        HTT(12,10, 3) =  0.20500D-02
        HTT(12,11, 3) =  0.64000D-03
        HTT(12,12, 3) =  0.39000D-02
        HTT(13, 1, 3) = -0.10900D-02
        HTT(13, 2, 3) = -0.19400D-02
        HTT(13, 3, 3) =  0.59000D-03
        HTT(13, 4, 3) = -0.26000D-03
        HTT(13, 5, 3) =  0.15200D-02
        HTT(13, 6, 3) =  0.15100D-02
        HTT(13, 7, 3) = -0.70000D-04
        HTT(13, 8, 3) = -0.90000D-04
        HTT(13, 9, 3) = -0.21100D-02
        HTT(13,10, 3) =  0.12000D-03
        HTT(13,11, 3) = -0.12300D-02
        HTT(13,12, 3) =  0.10000D-02
        HTT(13,13, 3) =  0.50300D-02
        HTT(14, 1, 3) =  0.70000D-03
        HTT(14, 2, 3) =  0.53000D-03
        HTT(14, 3, 3) = -0.31000D-03
        HTT(14, 4, 3) = -0.47000D-03
        HTT(14, 5, 3) = -0.22000D-03
        HTT(14, 6, 3) = -0.80000D-03
        HTT(14, 7, 3) =  0.74000D-03
        HTT(14, 8, 3) =  0.16000D-03
        HTT(14, 9, 3) = -0.19000D-03
        HTT(14,10, 3) =  0.13500D-02
        HTT(14,11, 3) = -0.10000D-04
        HTT(14,12, 3) =  0.15300D-02
        HTT(14,13, 3) = -0.81000D-03
        HTT(14,14, 3) =  0.39200D-02
        HTT(15, 1, 3) = -0.19100D-02
        HTT(15, 2, 3) = -0.28400D-02
        HTT(15, 3, 3) =  0.70000D-03
        HTT(15, 4, 3) = -0.23000D-03
        HTT(15, 5, 3) =  0.15000D-02
        HTT(15, 6, 3) =  0.21100D-02
        HTT(15, 7, 3) = -0.98000D-03
        HTT(15, 8, 3) = -0.37000D-03
        HTT(15, 9, 3) = -0.24000D-03
        HTT(15,10, 3) = -0.31800D-02
        HTT(15,11, 3) =  0.73000D-03
        HTT(15,12, 3) = -0.22200D-02
        HTT(15,13, 3) =  0.28200D-02
        HTT(15,14, 3) = -0.23300D-02
        HTT(15,15, 3) =  0.57300D-02
        HTT(16, 1, 3) = -0.12000D-03
        HTT(16, 2, 3) = -0.37000D-03
        HTT(16, 3, 3) = -0.19000D-03
        HTT(16, 4, 3) = -0.45000D-03
        HTT(16, 5, 3) = -0.25000D-03
        HTT(16, 6, 3) = -0.19000D-03
        HTT(16, 7, 3) = -0.17000D-03
        HTT(16, 8, 3) = -0.12000D-03
        HTT(16, 9, 3) =  0.16800D-02
        HTT(16,10, 3) = -0.19500D-02
        HTT(16,11, 3) =  0.19400D-02
        HTT(16,12, 3) = -0.16900D-02
        HTT(16,13, 3) = -0.30200D-02
        HTT(16,14, 3) =  0.24000D-02
        HTT(16,15, 3) =  0.58000D-03
        HTT(16,16, 3) =  0.60000D-02
        HTT(17, 1, 3) =  0.12100D-02
        HTT(17, 2, 3) =  0.15200D-02
        HTT(17, 3, 3) = -0.31000D-03
        HTT(17, 4, 3) =  0.00000D+00
        HTT(17, 5, 3) = -0.86000D-03
        HTT(17, 6, 3) = -0.18800D-02
        HTT(17, 7, 3) =  0.58000D-03
        HTT(17, 8, 3) = -0.44000D-03
        HTT(17, 9, 3) =  0.14500D-02
        HTT(17,10, 3) =  0.90000D-04
        HTT(17,11, 3) =  0.11300D-02
        HTT(17,12, 3) = -0.23000D-03
        HTT(17,13, 3) = -0.46000D-02
        HTT(17,14, 3) =  0.13000D-02
        HTT(17,15, 3) = -0.32400D-02
        HTT(17,16, 3) =  0.26500D-02
        HTT(17,17, 3) =  0.49600D-02
        HTT(18, 1, 3) =  0.14500D-02
        HTT(18, 2, 3) =  0.22200D-02
        HTT(18, 3, 3) = -0.10600D-02
        HTT(18, 4, 3) = -0.29000D-03
        HTT(18, 5, 3) = -0.17100D-02
        HTT(18, 6, 3) = -0.27100D-02
        HTT(18, 7, 3) =  0.67000D-03
        HTT(18, 8, 3) = -0.32000D-03
        HTT(18, 9, 3) =  0.10000D-02
        HTT(18,10, 3) =  0.82000D-03
        HTT(18,11, 3) =  0.19000D-03
        HTT(18,12, 3) =  0.10000D-04
        HTT(18,13, 3) = -0.32400D-02
        HTT(18,14, 3) =  0.19000D-02
        HTT(18,15, 3) = -0.30600D-02
        HTT(18,16, 3) =  0.20800D-02
        HTT(18,17, 3) =  0.29700D-02
        HTT(18,18, 3) =  0.56300D-02
        HTT(19, 1, 3) = -0.59000D-03
        HTT(19, 2, 3) = -0.96000D-03
        HTT(19, 3, 3) =  0.59000D-03
        HTT(19, 4, 3) =  0.21000D-03
        HTT(19, 5, 3) =  0.89000D-03
        HTT(19, 6, 3) =  0.43000D-03
        HTT(19, 7, 3) = -0.23000D-03
        HTT(19, 8, 3) = -0.69000D-03
        HTT(19, 9, 3) = -0.48000D-03
        HTT(19,10, 3) = -0.11500D-02
        HTT(19,11, 3) = -0.90000D-04
        HTT(19,12, 3) = -0.76000D-03
        HTT(19,13, 3) =  0.12600D-02
        HTT(19,14, 3) = -0.34500D-02
        HTT(19,15, 3) =  0.19200D-02
        HTT(19,16, 3) = -0.27800D-02
        HTT(19,17, 3) = -0.95000D-03
        HTT(19,18, 3) = -0.21800D-02
        HTT(19,19, 3) =  0.37700D-02
        HTT(20, 1, 3) = -0.35000D-03
        HTT(20, 2, 3) = -0.26000D-03
        HTT(20, 3, 3) = -0.16000D-03
        HTT(20, 4, 3) = -0.80000D-04
        HTT(20, 5, 3) =  0.40000D-04
        HTT(20, 6, 3) = -0.40000D-03
        HTT(20, 7, 3) = -0.14000D-03
        HTT(20, 8, 3) = -0.57000D-03
        HTT(20, 9, 3) = -0.93000D-03
        HTT(20,10, 3) = -0.42000D-03
        HTT(20,11, 3) = -0.10200D-02
        HTT(20,12, 3) = -0.51000D-03
        HTT(20,13, 3) =  0.26100D-02
        HTT(20,14, 3) = -0.28400D-02
        HTT(20,15, 3) =  0.21100D-02
        HTT(20,16, 3) = -0.33500D-02
        HTT(20,17, 3) = -0.29400D-02
        HTT(20,18, 3) =  0.47000D-03
        HTT(20,19, 3) =  0.25400D-02
        HTT(20,20, 3) =  0.59400D-02
        HTT(21, 1, 3) =  0.21700D-02
        HTT(21, 2, 3) =  0.15000D-02
        HTT(21, 3, 3) = -0.99000D-03
        HTT(21, 4, 3) = -0.16600D-02
        HTT(21, 5, 3) = -0.30200D-02
        HTT(21, 6, 3) = -0.69000D-03
        HTT(21, 7, 3) = -0.20000D-04
        HTT(21, 8, 3) =  0.23100D-02
        HTT(21, 9, 3) = -0.11900D-02
        HTT(21,10, 3) = -0.98000D-03
        HTT(21,11, 3) = -0.50000D-03
        HTT(21,12, 3) = -0.28000D-03
        HTT(21,13, 3) =  0.12600D-02
        HTT(21,14, 3) = -0.76000D-03
        HTT(21,15, 3) =  0.10500D-02
        HTT(21,16, 3) = -0.97000D-03
        HTT(21,17, 3) = -0.21400D-02
        HTT(21,18, 3) = -0.48000D-03
        HTT(21,19, 3) = -0.11000D-03
        HTT(21,20, 3) =  0.15500D-02
        HTT(21,21, 3) =  0.28900D-02
        HTT(22, 1, 3) =  0.27000D-03
        HTT(22, 2, 3) =  0.11300D-02
        HTT(22, 3, 3) = -0.15900D-02
        HTT(22, 4, 3) = -0.74000D-03
        HTT(22, 5, 3) = -0.60000D-04
        HTT(22, 6, 3) = -0.26200D-02
        HTT(22, 7, 3) =  0.17100D-02
        HTT(22, 8, 3) = -0.85000D-03
        HTT(22, 9, 3) = -0.85000D-03
        HTT(22,10, 3) =  0.40000D-04
        HTT(22,11, 3) = -0.17300D-02
        HTT(22,12, 3) = -0.85000D-03
        HTT(22,13, 3) =  0.12800D-02
        HTT(22,14, 3) = -0.15000D-03
        HTT(22,15, 3) =  0.40000D-03
        HTT(22,16, 3) = -0.10300D-02
        HTT(22,17, 3) = -0.10700D-02
        HTT(22,18, 3) =  0.56000D-03
        HTT(22,19, 3) =  0.37000D-03
        HTT(22,20, 3) =  0.20000D-02
        HTT(22,21, 3) =  0.44000D-03
        HTT(22,22, 3) =  0.31200D-02
        HTT(23, 1, 3) =  0.19200D-02
        HTT(23, 2, 3) =  0.79000D-03
        HTT(23, 3, 3) = -0.23000D-03
        HTT(23, 4, 3) = -0.13700D-02
        HTT(23, 5, 3) = -0.21700D-02
        HTT(23, 6, 3) =  0.14000D-03
        HTT(23, 7, 3) = -0.12000D-03
        HTT(23, 8, 3) =  0.21900D-02
        HTT(23, 9, 3) = -0.74000D-03
        HTT(23,10, 3) = -0.17100D-02
        HTT(23,11, 3) =  0.45000D-03
        HTT(23,12, 3) = -0.53000D-03
        HTT(23,13, 3) = -0.10000D-03
        HTT(23,14, 3) = -0.13700D-02
        HTT(23,15, 3) =  0.87000D-03
        HTT(23,16, 3) = -0.40000D-03
        HTT(23,17, 3) = -0.15000D-03
        HTT(23,18, 3) = -0.31500D-02
        HTT(23,19, 3) =  0.11200D-02
        HTT(23,20, 3) = -0.18800D-02
        HTT(23,21, 3) =  0.12200D-02
        HTT(23,22, 3) = -0.12000D-02
        HTT(23,23, 3) =  0.42400D-02
        HTT(24, 1, 3) =  0.30000D-04
        HTT(24, 2, 3) =  0.42000D-03
        HTT(24, 3, 3) = -0.84000D-03
        HTT(24, 4, 3) = -0.45000D-03
        HTT(24, 5, 3) =  0.79000D-03
        HTT(24, 6, 3) = -0.18000D-02
        HTT(24, 7, 3) =  0.16200D-02
        HTT(24, 8, 3) = -0.97000D-03
        HTT(24, 9, 3) = -0.39000D-03
        HTT(24,10, 3) = -0.70000D-03
        HTT(24,11, 3) = -0.79000D-03
        HTT(24,12, 3) = -0.11000D-02
        HTT(24,13, 3) = -0.80000D-04
        HTT(24,14, 3) = -0.76000D-03
        HTT(24,15, 3) =  0.22000D-03
        HTT(24,16, 3) = -0.45000D-03
        HTT(24,17, 3) =  0.92000D-03
        HTT(24,18, 3) = -0.21000D-02
        HTT(24,19, 3) =  0.16000D-02
        HTT(24,20, 3) = -0.14300D-02
        HTT(24,21, 3) = -0.12300D-02
        HTT(24,22, 3) =  0.14800D-02
        HTT(24,23, 3) =  0.18100D-02
        HTT(24,24, 3) =  0.45300D-02
C Hessian elements for anchor point 3
        HBS( 1, 1, 3) = -0.31760D-01
        HBS( 1, 2, 3) = -0.29360D-01
        HBS( 1, 3, 3) = -0.31820D-01
        HBS( 1, 4, 3) = -0.95900D-02
        HBS( 1, 5, 3) =  0.52300D-02
        HBS( 1, 6, 3) =  0.39930D-01
        HBS( 1, 7, 3) = -0.45000D-03
        HBS( 1, 8, 3) =  0.40450D-01
        HBS( 1, 9, 3) = -0.18000D-03
        HBS( 1,10, 3) = -0.96900D-02
        HBS( 1,11, 3) = -0.42000D-03
        HBS( 1,12, 3) =  0.47400D-02
        HBS( 2, 1, 3) =  0.49930D-01
        HBS( 2, 2, 3) = -0.13930D-01
        HBS( 2, 3, 3) =  0.99600D-02
        HBS( 2, 4, 3) =  0.18000D-02
        HBS( 2, 5, 3) = -0.12200D-02
        HBS( 2, 6, 3) = -0.27810D-01
        HBS( 2, 7, 3) =  0.11800D-02
        HBS( 2, 8, 3) = -0.14650D-01
        HBS( 2, 9, 3) =  0.00000D+00
        HBS( 2,10, 3) =  0.41300D-02
        HBS( 2,11, 3) = -0.69000D-03
        HBS( 2,12, 3) = -0.89000D-03
        HBS( 3, 1, 3) = -0.18180D-01
        HBS( 3, 2, 3) =  0.43280D-01
        HBS( 3, 3, 3) =  0.21860D-01
        HBS( 3, 4, 3) =  0.77900D-02
        HBS( 3, 5, 3) = -0.40100D-02
        HBS( 3, 6, 3) = -0.12120D-01
        HBS( 3, 7, 3) = -0.73000D-03
        HBS( 3, 8, 3) = -0.25800D-01
        HBS( 3, 9, 3) =  0.18000D-03
        HBS( 3,10, 3) =  0.55600D-02
        HBS( 3,11, 3) =  0.11100D-02
        HBS( 3,12, 3) = -0.38500D-02
        HBS( 4, 1, 3) = -0.20820D-01
        HBS( 4, 2, 3) = -0.12100D-01
        HBS( 4, 3, 3) =  0.15180D-01
        HBS( 4, 4, 3) = -0.29840D-01
        HBS( 4, 5, 3) = -0.97300D-02
        HBS( 4, 6, 3) = -0.57000D-02
        HBS( 4, 7, 3) =  0.49500D-02
        HBS( 4, 8, 3) =  0.31150D-01
        HBS( 4, 9, 3) = -0.44000D-03
        HBS( 4,10, 3) =  0.40270D-01
        HBS( 4,11, 3) = -0.59000D-03
        HBS( 4,12, 3) =  0.29000D-03
        HBS( 5, 1, 3) =  0.26090D-01
        HBS( 5, 2, 3) =  0.58900D-02
        HBS( 5, 3, 3) = -0.14520D-01
        HBS( 5, 4, 3) =  0.25400D-02
        HBS( 5, 5, 3) =  0.69400D-02
        HBS( 5, 6, 3) =  0.30900D-02
        HBS( 5, 7, 3) = -0.21700D-02
        HBS( 5, 8, 3) = -0.15190D-01
        HBS( 5, 9, 3) = -0.39000D-03
        HBS( 5,10, 3) = -0.22140D-01
        HBS( 5,11, 3) =  0.28000D-03
        HBS( 5,12, 3) =  0.90000D-04
        HBS( 6, 1, 3) = -0.52700D-02
        HBS( 6, 2, 3) =  0.62100D-02
        HBS( 6, 3, 3) = -0.66000D-03
        HBS( 6, 4, 3) =  0.27300D-01
        HBS( 6, 5, 3) =  0.27900D-02
        HBS( 6, 6, 3) =  0.26100D-02
        HBS( 6, 7, 3) = -0.27800D-02
        HBS( 6, 8, 3) = -0.15960D-01
        HBS( 6, 9, 3) =  0.83000D-03
        HBS( 6,10, 3) = -0.18130D-01
        HBS( 6,11, 3) =  0.31000D-03
        HBS( 6,12, 3) = -0.38000D-03
        HBS( 7, 1, 3) = -0.90800D-02
        HBS( 7, 2, 3) =  0.35870D-01
        HBS( 7, 3, 3) = -0.10000D-04
        HBS( 7, 4, 3) = -0.22920D-01
        HBS( 7, 5, 3) =  0.61900D-02
        HBS( 7, 6, 3) = -0.28970D-01
        HBS( 7, 7, 3) = -0.87500D-02
        HBS( 7, 8, 3) = -0.92900D-02
        HBS( 7, 9, 3) =  0.49700D-02
        HBS( 7,10, 3) =  0.32780D-01
        HBS( 7,11, 3) =  0.60000D-04
        HBS( 7,12, 3) = -0.10700D-02
        HBS( 8, 1, 3) =  0.36900D-02
        HBS( 8, 2, 3) = -0.18590D-01
        HBS( 8, 3, 3) =  0.13800D-02
        HBS( 8, 4, 3) =  0.25350D-01
        HBS( 8, 5, 3) = -0.35000D-02
        HBS( 8, 6, 3) = -0.42000D-03
        HBS( 8, 7, 3) =  0.45000D-02
        HBS( 8, 8, 3) =  0.50400D-02
        HBS( 8, 9, 3) = -0.20100D-02
        HBS( 8,10, 3) = -0.14680D-01
        HBS( 8,11, 3) = -0.62000D-03
        HBS( 8,12, 3) =  0.55000D-03
        HBS( 9, 1, 3) =  0.53900D-02
        HBS( 9, 2, 3) = -0.17270D-01
        HBS( 9, 3, 3) = -0.13700D-02
        HBS( 9, 4, 3) = -0.24400D-02
        HBS( 9, 5, 3) = -0.26800D-02
        HBS( 9, 6, 3) =  0.29390D-01
        HBS( 9, 7, 3) =  0.42400D-02
        HBS( 9, 8, 3) =  0.42500D-02
        HBS( 9, 9, 3) = -0.29600D-02
        HBS( 9,10, 3) = -0.18100D-01
        HBS( 9,11, 3) =  0.55000D-03
        HBS( 9,12, 3) =  0.51000D-03
        HBS(10, 1, 3) =  0.36130D-01
        HBS(10, 2, 3) =  0.38360D-01
        HBS(10, 3, 3) =  0.36300D-02
        HBS(10, 4, 3) = -0.11890D-01
        HBS(10, 5, 3) = -0.10500D-02
        HBS(10, 6, 3) = -0.26840D-01
        HBS(10, 7, 3) =  0.48100D-02
        HBS(10, 8, 3) = -0.27100D-01
        HBS(10, 9, 3) = -0.88200D-02
        HBS(10,10, 3) = -0.12140D-01
        HBS(10,11, 3) =  0.47600D-02
        HBS(10,12, 3) = -0.77000D-03
        HBS(11, 1, 3) = -0.20190D-01
        HBS(11, 2, 3) = -0.17180D-01
        HBS(11, 3, 3) = -0.21300D-02
        HBS(11, 4, 3) =  0.57400D-02
        HBS(11, 5, 3) =  0.12900D-02
        HBS(11, 6, 3) =  0.27000D-01
        HBS(11, 7, 3) = -0.28400D-02
        HBS(11, 8, 3) =  0.47000D-03
        HBS(11, 9, 3) =  0.44300D-02
        HBS(11,10, 3) =  0.60700D-02
        HBS(11,11, 3) = -0.19800D-02
        HBS(11,12, 3) = -0.18000D-03
        HBS(12, 1, 3) = -0.15940D-01
        HBS(12, 2, 3) = -0.21180D-01
        HBS(12, 3, 3) = -0.14900D-02
        HBS(12, 4, 3) =  0.61500D-02
        HBS(12, 5, 3) = -0.24000D-03
        HBS(12, 6, 3) = -0.16000D-03
        HBS(12, 7, 3) = -0.19600D-02
        HBS(12, 8, 3) =  0.26630D-01
        HBS(12, 9, 3) =  0.43900D-02
        HBS(12,10, 3) =  0.60700D-02
        HBS(12,11, 3) = -0.27800D-02
        HBS(12,12, 3) =  0.94000D-03
        HBS(13, 1, 3) =  0.35720D-01
        HBS(13, 2, 3) = -0.96900D-02
        HBS(13, 3, 3) =  0.29800D-02
        HBS(13, 4, 3) =  0.32310D-01
        HBS(13, 5, 3) = -0.94000D-03
        HBS(13, 6, 3) = -0.10140D-01
        HBS(13, 7, 3) =  0.70000D-04
        HBS(13, 8, 3) = -0.29530D-01
        HBS(13, 9, 3) =  0.50200D-02
        HBS(13,10, 3) = -0.23420D-01
        HBS(13,11, 3) = -0.86600D-02
        HBS(13,12, 3) =  0.53400D-02
        HBS(14, 1, 3) = -0.17210D-01
        HBS(14, 2, 3) =  0.56100D-02
        HBS(14, 3, 3) = -0.23300D-02
        HBS(14, 4, 3) = -0.17680D-01
        HBS(14, 5, 3) =  0.48000D-03
        HBS(14, 6, 3) =  0.47700D-02
        HBS(14, 7, 3) =  0.56000D-03
        HBS(14, 8, 3) =  0.29630D-01
        HBS(14, 9, 3) = -0.29700D-02
        HBS(14,10, 3) = -0.24000D-02
        HBS(14,11, 3) =  0.43000D-02
        HBS(14,12, 3) = -0.22600D-02
        HBS(15, 1, 3) = -0.18510D-01
        HBS(15, 2, 3) =  0.40800D-02
        HBS(15, 3, 3) = -0.65000D-03
        HBS(15, 4, 3) = -0.14630D-01
        HBS(15, 5, 3) =  0.46000D-03
        HBS(15, 6, 3) =  0.53700D-02
        HBS(15, 7, 3) = -0.63000D-03
        HBS(15, 8, 3) = -0.10000D-03
        HBS(15, 9, 3) = -0.20500D-02
        HBS(15,10, 3) =  0.25820D-01
        HBS(15,11, 3) =  0.43700D-02
        HBS(15,12, 3) = -0.30900D-02
        HBS(16, 1, 3) = -0.10200D-01
        HBS(16, 2, 3) = -0.23080D-01
        HBS(16, 3, 3) =  0.10040D-01
        HBS(16, 4, 3) =  0.41930D-01
        HBS(16, 5, 3) =  0.31000D-03
        HBS(16, 6, 3) =  0.31710D-01
        HBS(16, 7, 3) = -0.63000D-03
        HBS(16, 8, 3) = -0.56800D-02
        HBS(16, 9, 3) = -0.55000D-03
        HBS(16,10, 3) = -0.27810D-01
        HBS(16,11, 3) =  0.48500D-02
        HBS(16,12, 3) = -0.85400D-02
        HBS(17, 1, 3) =  0.41100D-02
        HBS(17, 2, 3) =  0.25000D-01
        HBS(17, 3, 3) = -0.74400D-02
        HBS(17, 4, 3) = -0.22980D-01
        HBS(17, 5, 3) =  0.60000D-03
        HBS(17, 6, 3) = -0.14860D-01
        HBS(17, 7, 3) =  0.32000D-03
        HBS(17, 8, 3) =  0.33600D-02
        HBS(17, 9, 3) = -0.32000D-03
        HBS(17,10, 3) =  0.20400D-02
        HBS(17,11, 3) = -0.21000D-02
        HBS(17,12, 3) =  0.68800D-02
        HBS(18, 1, 3) =  0.61000D-02
        HBS(18, 2, 3) = -0.19200D-02
        HBS(18, 3, 3) = -0.26000D-02
        HBS(18, 4, 3) = -0.18950D-01
        HBS(18, 5, 3) = -0.92000D-03
        HBS(18, 6, 3) = -0.16850D-01
        HBS(18, 7, 3) =  0.31000D-03
        HBS(18, 8, 3) =  0.23200D-02
        HBS(18, 9, 3) =  0.87000D-03
        HBS(18,10, 3) =  0.25770D-01
        HBS(18,11, 3) = -0.27500D-02
        HBS(18,12, 3) =  0.16600D-02
C Hessian elements for anchor point 4
        HSS( 1, 1, 4) =  0.34759D+00
        HSS( 2, 1, 4) =  0.17810D-01
        HSS( 2, 2, 4) =  0.35623D+00
        HSS( 3, 1, 4) =  0.29010D-01
        HSS( 3, 2, 4) =  0.32710D-01
        HSS( 3, 3, 4) =  0.46124D+00
        HSS( 4, 1, 4) = -0.19600D-02
        HSS( 4, 2, 4) =  0.14870D-01
        HSS( 4, 3, 4) = -0.89100D-02
        HSS( 4, 4, 4) =  0.38581D+00
        HSS( 5, 1, 4) =  0.38800D-02
        HSS( 5, 2, 4) = -0.77000D-03
        HSS( 5, 3, 4) =  0.71700D-02
        HSS( 5, 4, 4) =  0.36500D-02
        HSS( 5, 5, 4) =  0.35807D+00
        HSS( 6, 1, 4) =  0.33890D-01
        HSS( 6, 2, 4) =  0.44570D-01
        HSS( 6, 3, 4) = -0.43100D-02
        HSS( 6, 4, 4) = -0.36100D-02
        HSS( 6, 5, 4) =  0.55000D-03
        HSS( 6, 6, 4) =  0.35931D+00
        HSS( 7, 1, 4) = -0.54000D-03
        HSS( 7, 2, 4) = -0.63000D-03
        HSS( 7, 3, 4) =  0.57000D-03
        HSS( 7, 4, 4) =  0.39700D-02
        HSS( 7, 5, 4) =  0.39000D-03
        HSS( 7, 6, 4) =  0.39800D-02
        HSS( 7, 7, 4) =  0.37091D+00
        HSS( 8, 1, 4) =  0.42030D-01
        HSS( 8, 2, 4) =  0.35060D-01
        HSS( 8, 3, 4) = -0.52800D-02
        HSS( 8, 4, 4) =  0.24280D-01
        HSS( 8, 5, 4) = -0.11100D-02
        HSS( 8, 6, 4) =  0.17880D-01
        HSS( 8, 7, 4) = -0.10400D-02
        HSS( 8, 8, 4) =  0.35888D+00
        HSS( 9, 1, 4) = -0.53000D-03
        HSS( 9, 2, 4) = -0.57000D-03
        HSS( 9, 3, 4) =  0.26000D-03
        HSS( 9, 4, 4) = -0.39000D-03
        HSS( 9, 5, 4) =  0.19000D-03
        HSS( 9, 6, 4) =  0.40100D-02
        HSS( 9, 7, 4) =  0.49000D-03
        HSS( 9, 8, 4) =  0.39600D-02
        HSS( 9, 9, 4) =  0.36356D+00
        HSS(10, 1, 4) =  0.19520D-01
        HSS(10, 2, 4) = -0.69300D-02
        HSS(10, 3, 4) = -0.48800D-02
        HSS(10, 4, 4) =  0.64010D-01
        HSS(10, 5, 4) = -0.43000D-03
        HSS(10, 6, 4) =  0.21230D-01
        HSS(10, 7, 4) = -0.10900D-02
        HSS(10, 8, 4) = -0.77000D-03
        HSS(10, 9, 4) = -0.19000D-03
        HSS(10,10, 4) =  0.38116D+00
        HSS(11, 1, 4) = -0.69000D-03
        HSS(11, 2, 4) = -0.50000D-03
        HSS(11, 3, 4) =  0.52000D-03
        HSS(11, 4, 4) = -0.11600D-02
        HSS(11, 5, 4) =  0.80000D-04
        HSS(11, 6, 4) = -0.10000D-02
        HSS(11, 7, 4) =  0.22000D-03
        HSS(11, 8, 4) =  0.39700D-02
        HSS(11, 9, 4) =  0.48000D-03
        HSS(11,10, 4) =  0.39800D-02
        HSS(11,11, 4) =  0.37152D+00
        HSS(12, 1, 4) = -0.16500D-02
        HSS(12, 2, 4) =  0.39400D-02
        HSS(12, 3, 4) = -0.19000D-03
        HSS(12, 4, 4) = -0.10000D-03
        HSS(12, 5, 4) =  0.35000D-03
        HSS(12, 6, 4) = -0.10400D-02
        HSS(12, 7, 4) =  0.30000D-04
        HSS(12, 8, 4) =  0.57000D-03
        HSS(12, 9, 4) =  0.12000D-03
        HSS(12,10, 4) =  0.23600D-02
        HSS(12,11, 4) =  0.40000D-03
        HSS(12,12, 4) =  0.37254D+00
C Hessian elements for anchor point 4
        HBB( 1, 1, 4) =  0.82450D-01
        HBB( 2, 1, 4) = -0.46090D-01
        HBB( 2, 2, 4) =  0.14487D+00
        HBB( 3, 1, 4) = -0.36360D-01
        HBB( 3, 2, 4) = -0.98780D-01
        HBB( 3, 3, 4) =  0.13515D+00
        HBB( 4, 1, 4) = -0.36000D-01
        HBB( 4, 2, 4) =  0.23310D-01
        HBB( 4, 3, 4) =  0.12690D-01
        HBB( 4, 4, 4) =  0.66170D-01
        HBB( 5, 1, 4) =  0.20240D-01
        HBB( 5, 2, 4) = -0.12480D-01
        HBB( 5, 3, 4) = -0.77700D-02
        HBB( 5, 4, 4) = -0.33030D-01
        HBB( 5, 5, 4) =  0.71950D-01
        HBB( 6, 1, 4) =  0.15760D-01
        HBB( 6, 2, 4) = -0.10840D-01
        HBB( 6, 3, 4) = -0.49200D-02
        HBB( 6, 4, 4) = -0.33140D-01
        HBB( 6, 5, 4) = -0.38920D-01
        HBB( 6, 6, 4) =  0.72060D-01
        HBB( 7, 1, 4) = -0.71500D-02
        HBB( 7, 2, 4) = -0.39000D-02
        HBB( 7, 3, 4) =  0.11050D-01
        HBB( 7, 4, 4) = -0.25740D-01
        HBB( 7, 5, 4) =  0.93100D-02
        HBB( 7, 6, 4) =  0.16430D-01
        HBB( 7, 7, 4) =  0.64420D-01
        HBB( 8, 1, 4) = -0.22800D-02
        HBB( 8, 2, 4) =  0.74300D-02
        HBB( 8, 3, 4) = -0.51500D-02
        HBB( 8, 4, 4) =  0.15690D-01
        HBB( 8, 5, 4) = -0.49500D-02
        HBB( 8, 6, 4) = -0.10740D-01
        HBB( 8, 7, 4) = -0.31490D-01
        HBB( 8, 8, 4) =  0.72950D-01
        HBB( 9, 1, 4) =  0.94300D-02
        HBB( 9, 2, 4) = -0.35300D-02
        HBB( 9, 3, 4) = -0.59100D-02
        HBB( 9, 4, 4) =  0.10050D-01
        HBB( 9, 5, 4) = -0.43500D-02
        HBB( 9, 6, 4) = -0.56900D-02
        HBB( 9, 7, 4) = -0.32930D-01
        HBB( 9, 8, 4) = -0.41460D-01
        HBB( 9, 9, 4) =  0.74390D-01
        HBB(10, 1, 4) =  0.45500D-02
        HBB(10, 2, 4) = -0.20200D-02
        HBB(10, 3, 4) = -0.25300D-02
        HBB(10, 4, 4) = -0.64100D-02
        HBB(10, 5, 4) =  0.91900D-02
        HBB(10, 6, 4) = -0.27700D-02
        HBB(10, 7, 4) = -0.33410D-01
        HBB(10, 8, 4) =  0.13410D-01
        HBB(10, 9, 4) =  0.20000D-01
        HBB(10,10, 4) =  0.75440D-01
        HBB(11, 1, 4) = -0.22100D-02
        HBB(11, 2, 4) =  0.12300D-02
        HBB(11, 3, 4) =  0.98000D-03
        HBB(11, 4, 4) = -0.23100D-02
        HBB(11, 5, 4) = -0.32200D-02
        HBB(11, 6, 4) =  0.55300D-02
        HBB(11, 7, 4) =  0.20020D-01
        HBB(11, 8, 4) = -0.74100D-02
        HBB(11, 9, 4) = -0.12610D-01
        HBB(11,10, 4) = -0.37790D-01
        HBB(11,11, 4) =  0.78360D-01
        HBB(12, 1, 4) = -0.23400D-02
        HBB(12, 2, 4) =  0.79000D-03
        HBB(12, 3, 4) =  0.15500D-02
        HBB(12, 4, 4) =  0.87200D-02
        HBB(12, 5, 4) = -0.59700D-02
        HBB(12, 6, 4) = -0.27500D-02
        HBB(12, 7, 4) =  0.13390D-01
        HBB(12, 8, 4) = -0.60000D-02
        HBB(12, 9, 4) = -0.73900D-02
        HBB(12,10, 4) = -0.37640D-01
        HBB(12,11, 4) = -0.40570D-01
        HBB(12,12, 4) =  0.78220D-01
        HBB(13, 1, 4) = -0.67600D-02
        HBB(13, 2, 4) =  0.10710D-01
        HBB(13, 3, 4) = -0.39500D-02
        HBB(13, 4, 4) = -0.11500D-02
        HBB(13, 5, 4) =  0.12900D-02
        HBB(13, 6, 4) = -0.14000D-03
        HBB(13, 7, 4) =  0.29100D-02
        HBB(13, 8, 4) =  0.41300D-02
        HBB(13, 9, 4) = -0.70300D-02
        HBB(13,10, 4) = -0.33760D-01
        HBB(13,11, 4) =  0.13510D-01
        HBB(13,12, 4) =  0.20240D-01
        HBB(13,13, 4) =  0.65300D-01
        HBB(14, 1, 4) =  0.91200D-02
        HBB(14, 2, 4) = -0.60500D-02
        HBB(14, 3, 4) = -0.30700D-02
        HBB(14, 4, 4) =  0.66000D-03
        HBB(14, 5, 4) = -0.85000D-03
        HBB(14, 6, 4) =  0.19000D-03
        HBB(14, 7, 4) = -0.71100D-02
        HBB(14, 8, 4) = -0.51000D-03
        HBB(14, 9, 4) =  0.76200D-02
        HBB(14,10, 4) =  0.20030D-01
        HBB(14,11, 4) = -0.73800D-02
        HBB(14,12, 4) = -0.12650D-01
        HBB(14,13, 4) = -0.33180D-01
        HBB(14,14, 4) =  0.74350D-01
        HBB(15, 1, 4) = -0.23600D-02
        HBB(15, 2, 4) = -0.46500D-02
        HBB(15, 3, 4) =  0.70200D-02
        HBB(15, 4, 4) =  0.49000D-03
        HBB(15, 5, 4) = -0.44000D-03
        HBB(15, 6, 4) = -0.50000D-04
        HBB(15, 7, 4) =  0.42000D-02
        HBB(15, 8, 4) = -0.36200D-02
        HBB(15, 9, 4) = -0.59000D-03
        HBB(15,10, 4) =  0.13730D-01
        HBB(15,11, 4) = -0.61300D-02
        HBB(15,12, 4) = -0.76000D-02
        HBB(15,13, 4) = -0.32120D-01
        HBB(15,14, 4) = -0.41170D-01
        HBB(15,15, 4) =  0.73290D-01
        HBB(16, 1, 4) = -0.37090D-01
        HBB(16, 2, 4) =  0.17990D-01
        HBB(16, 3, 4) =  0.19100D-01
        HBB(16, 4, 4) =  0.31400D-02
        HBB(16, 5, 4) = -0.70000D-02
        HBB(16, 6, 4) =  0.38600D-02
        HBB(16, 7, 4) = -0.10300D-02
        HBB(16, 8, 4) =  0.54000D-03
        HBB(16, 9, 4) =  0.48000D-03
        HBB(16,10, 4) = -0.64000D-02
        HBB(16,11, 4) =  0.87800D-02
        HBB(16,12, 4) = -0.23700D-02
        HBB(16,13, 4) = -0.26530D-01
        HBB(16,14, 4) =  0.10480D-01
        HBB(16,15, 4) =  0.16060D-01
        HBB(16,16, 4) =  0.67910D-01
        HBB(17, 1, 4) =  0.22610D-01
        HBB(17, 2, 4) = -0.90200D-02
        HBB(17, 3, 4) = -0.13590D-01
        HBB(17, 4, 4) = -0.76000D-02
        HBB(17, 5, 4) =  0.77900D-02
        HBB(17, 6, 4) = -0.19000D-03
        HBB(17, 7, 4) =  0.87000D-03
        HBB(17, 8, 4) = -0.25000D-03
        HBB(17, 9, 4) = -0.62000D-03
        HBB(17,10, 4) =  0.91600D-02
        HBB(17,11, 4) = -0.60700D-02
        HBB(17,12, 4) = -0.31000D-02
        HBB(17,13, 4) =  0.96800D-02
        HBB(17,14, 4) = -0.46000D-02
        HBB(17,15, 4) = -0.50800D-02
        HBB(17,16, 4) = -0.34720D-01
        HBB(17,17, 4) =  0.72390D-01
        HBB(18, 1, 4) =  0.14470D-01
        HBB(18, 2, 4) = -0.89700D-02
        HBB(18, 3, 4) = -0.55000D-02
        HBB(18, 4, 4) =  0.44600D-02
        HBB(18, 5, 4) = -0.79000D-03
        HBB(18, 6, 4) = -0.36700D-02
        HBB(18, 7, 4) =  0.16000D-03
        HBB(18, 8, 4) = -0.29000D-03
        HBB(18, 9, 4) =  0.13000D-03
        HBB(18,10, 4) = -0.27600D-02
        HBB(18,11, 4) = -0.27100D-02
        HBB(18,12, 4) =  0.54700D-02
        HBB(18,13, 4) =  0.16860D-01
        HBB(18,14, 4) = -0.58800D-02
        HBB(18,15, 4) = -0.10980D-01
        HBB(18,16, 4) = -0.33190D-01
        HBB(18,17, 4) = -0.37670D-01
        HBB(18,18, 4) =  0.70860D-01
C Hessian elements for anchor point 4
        HTT( 1, 1, 4) =  0.60400D-02
        HTT( 2, 1, 4) =  0.38400D-02
        HTT( 2, 2, 4) =  0.50900D-02
        HTT( 3, 1, 4) = -0.18500D-02
        HTT( 3, 2, 4) = -0.33100D-02
        HTT( 3, 3, 4) =  0.75100D-02
        HTT( 4, 1, 4) = -0.40500D-02
        HTT( 4, 2, 4) = -0.20700D-02
        HTT( 4, 3, 4) =  0.60500D-02
        HTT( 4, 4, 4) =  0.80200D-02
        HTT( 5, 1, 4) = -0.59200D-02
        HTT( 5, 2, 4) = -0.42600D-02
        HTT( 5, 3, 4) =  0.21200D-02
        HTT( 5, 4, 4) =  0.37800D-02
        HTT( 5, 5, 4) =  0.65600D-02
        HTT( 6, 1, 4) = -0.41200D-02
        HTT( 6, 2, 4) = -0.39000D-02
        HTT( 6, 3, 4) =  0.26900D-02
        HTT( 6, 4, 4) =  0.29100D-02
        HTT( 6, 5, 4) =  0.37400D-02
        HTT( 6, 6, 4) =  0.55900D-02
        HTT( 7, 1, 4) =  0.15700D-02
        HTT( 7, 2, 4) =  0.25400D-02
        HTT( 7, 3, 4) = -0.67700D-02
        HTT( 7, 4, 4) = -0.58000D-02
        HTT( 7, 5, 4) = -0.10700D-02
        HTT( 7, 6, 4) = -0.27200D-02
        HTT( 7, 7, 4) =  0.68600D-02
        HTT( 8, 1, 4) =  0.33800D-02
        HTT( 8, 2, 4) =  0.29000D-02
        HTT( 8, 3, 4) = -0.62000D-02
        HTT( 8, 4, 4) = -0.66800D-02
        HTT( 8, 5, 4) = -0.38900D-02
        HTT( 8, 6, 4) = -0.88000D-03
        HTT( 8, 7, 4) =  0.52100D-02
        HTT( 8, 8, 4) =  0.82200D-02
        HTT( 9, 1, 4) = -0.23900D-02
        HTT( 9, 2, 4) = -0.68000D-03
        HTT( 9, 3, 4) =  0.45000D-03
        HTT( 9, 4, 4) =  0.21600D-02
        HTT( 9, 5, 4) =  0.17400D-02
        HTT( 9, 6, 4) =  0.14100D-02
        HTT( 9, 7, 4) = -0.96000D-03
        HTT( 9, 8, 4) = -0.12900D-02
        HTT( 9, 9, 4) =  0.24600D-02
        HTT(10, 1, 4) = -0.15600D-02
        HTT(10, 2, 4) =  0.23000D-03
        HTT(10, 3, 4) =  0.34000D-03
        HTT(10, 4, 4) =  0.21300D-02
        HTT(10, 5, 4) =  0.17600D-02
        HTT(10, 6, 4) =  0.79000D-03
        HTT(10, 7, 4) = -0.40000D-04
        HTT(10, 8, 4) = -0.10100D-02
        HTT(10, 9, 4) =  0.57000D-03
        HTT(10,10, 4) =  0.39100D-02
        HTT(11, 1, 4) = -0.10000D-03
        HTT(11, 2, 4) = -0.19700D-02
        HTT(11, 3, 4) =  0.19800D-02
        HTT(11, 4, 4) =  0.11000D-03
        HTT(11, 5, 4) =  0.10000D-04
        HTT(11, 6, 4) =  0.11900D-02
        HTT(11, 7, 4) = -0.19700D-02
        HTT(11, 8, 4) = -0.79000D-03
        HTT(11, 9, 4) =  0.68000D-03
        HTT(11,10, 4) = -0.12900D-02
        HTT(11,11, 4) =  0.26200D-02
        HTT(12, 1, 4) =  0.72000D-03
        HTT(12, 2, 4) = -0.10500D-02
        HTT(12, 3, 4) =  0.18600D-02
        HTT(12, 4, 4) =  0.90000D-04
        HTT(12, 5, 4) =  0.30000D-04
        HTT(12, 6, 4) =  0.57000D-03
        HTT(12, 7, 4) = -0.10500D-02
        HTT(12, 8, 4) = -0.51000D-03
        HTT(12, 9, 4) = -0.12000D-02
        HTT(12,10, 4) =  0.20500D-02
        HTT(12,11, 4) =  0.64000D-03
        HTT(12,12, 4) =  0.39000D-02
        HTT(13, 1, 4) = -0.10900D-02
        HTT(13, 2, 4) = -0.19400D-02
        HTT(13, 3, 4) =  0.59000D-03
        HTT(13, 4, 4) = -0.26000D-03
        HTT(13, 5, 4) =  0.15200D-02
        HTT(13, 6, 4) =  0.15100D-02
        HTT(13, 7, 4) = -0.70000D-04
        HTT(13, 8, 4) = -0.90000D-04
        HTT(13, 9, 4) = -0.21100D-02
        HTT(13,10, 4) =  0.12000D-03
        HTT(13,11, 4) = -0.12300D-02
        HTT(13,12, 4) =  0.10000D-02
        HTT(13,13, 4) =  0.50300D-02
        HTT(14, 1, 4) =  0.70000D-03
        HTT(14, 2, 4) =  0.53000D-03
        HTT(14, 3, 4) = -0.31000D-03
        HTT(14, 4, 4) = -0.47000D-03
        HTT(14, 5, 4) = -0.22000D-03
        HTT(14, 6, 4) = -0.80000D-03
        HTT(14, 7, 4) =  0.74000D-03
        HTT(14, 8, 4) =  0.16000D-03
        HTT(14, 9, 4) = -0.19000D-03
        HTT(14,10, 4) =  0.13500D-02
        HTT(14,11, 4) = -0.10000D-04
        HTT(14,12, 4) =  0.15300D-02
        HTT(14,13, 4) = -0.81000D-03
        HTT(14,14, 4) =  0.39200D-02
        HTT(15, 1, 4) = -0.19100D-02
        HTT(15, 2, 4) = -0.28400D-02
        HTT(15, 3, 4) =  0.70000D-03
        HTT(15, 4, 4) = -0.23000D-03
        HTT(15, 5, 4) =  0.15000D-02
        HTT(15, 6, 4) =  0.21100D-02
        HTT(15, 7, 4) = -0.98000D-03
        HTT(15, 8, 4) = -0.37000D-03
        HTT(15, 9, 4) = -0.24000D-03
        HTT(15,10, 4) = -0.31800D-02
        HTT(15,11, 4) =  0.73000D-03
        HTT(15,12, 4) = -0.22200D-02
        HTT(15,13, 4) =  0.28200D-02
        HTT(15,14, 4) = -0.23300D-02
        HTT(15,15, 4) =  0.57300D-02
        HTT(16, 1, 4) = -0.12000D-03
        HTT(16, 2, 4) = -0.37000D-03
        HTT(16, 3, 4) = -0.19000D-03
        HTT(16, 4, 4) = -0.45000D-03
        HTT(16, 5, 4) = -0.25000D-03
        HTT(16, 6, 4) = -0.19000D-03
        HTT(16, 7, 4) = -0.17000D-03
        HTT(16, 8, 4) = -0.12000D-03
        HTT(16, 9, 4) =  0.16800D-02
        HTT(16,10, 4) = -0.19500D-02
        HTT(16,11, 4) =  0.19400D-02
        HTT(16,12, 4) = -0.16900D-02
        HTT(16,13, 4) = -0.30200D-02
        HTT(16,14, 4) =  0.24000D-02
        HTT(16,15, 4) =  0.58000D-03
        HTT(16,16, 4) =  0.60000D-02
        HTT(17, 1, 4) =  0.12100D-02
        HTT(17, 2, 4) =  0.15200D-02
        HTT(17, 3, 4) = -0.31000D-03
        HTT(17, 4, 4) =  0.00000D+00
        HTT(17, 5, 4) = -0.86000D-03
        HTT(17, 6, 4) = -0.18800D-02
        HTT(17, 7, 4) =  0.58000D-03
        HTT(17, 8, 4) = -0.44000D-03
        HTT(17, 9, 4) =  0.14500D-02
        HTT(17,10, 4) =  0.90000D-04
        HTT(17,11, 4) =  0.11300D-02
        HTT(17,12, 4) = -0.23000D-03
        HTT(17,13, 4) = -0.46000D-02
        HTT(17,14, 4) =  0.13000D-02
        HTT(17,15, 4) = -0.32400D-02
        HTT(17,16, 4) =  0.26500D-02
        HTT(17,17, 4) =  0.49600D-02
        HTT(18, 1, 4) =  0.14500D-02
        HTT(18, 2, 4) =  0.22200D-02
        HTT(18, 3, 4) = -0.10600D-02
        HTT(18, 4, 4) = -0.29000D-03
        HTT(18, 5, 4) = -0.17100D-02
        HTT(18, 6, 4) = -0.27100D-02
        HTT(18, 7, 4) =  0.67000D-03
        HTT(18, 8, 4) = -0.32000D-03
        HTT(18, 9, 4) =  0.10000D-02
        HTT(18,10, 4) =  0.82000D-03
        HTT(18,11, 4) =  0.19000D-03
        HTT(18,12, 4) =  0.10000D-04
        HTT(18,13, 4) = -0.32400D-02
        HTT(18,14, 4) =  0.19000D-02
        HTT(18,15, 4) = -0.30600D-02
        HTT(18,16, 4) =  0.20800D-02
        HTT(18,17, 4) =  0.29700D-02
        HTT(18,18, 4) =  0.56300D-02
        HTT(19, 1, 4) = -0.59000D-03
        HTT(19, 2, 4) = -0.96000D-03
        HTT(19, 3, 4) =  0.59000D-03
        HTT(19, 4, 4) =  0.21000D-03
        HTT(19, 5, 4) =  0.89000D-03
        HTT(19, 6, 4) =  0.43000D-03
        HTT(19, 7, 4) = -0.23000D-03
        HTT(19, 8, 4) = -0.69000D-03
        HTT(19, 9, 4) = -0.48000D-03
        HTT(19,10, 4) = -0.11500D-02
        HTT(19,11, 4) = -0.90000D-04
        HTT(19,12, 4) = -0.76000D-03
        HTT(19,13, 4) =  0.12600D-02
        HTT(19,14, 4) = -0.34500D-02
        HTT(19,15, 4) =  0.19200D-02
        HTT(19,16, 4) = -0.27800D-02
        HTT(19,17, 4) = -0.95000D-03
        HTT(19,18, 4) = -0.21800D-02
        HTT(19,19, 4) =  0.37700D-02
        HTT(20, 1, 4) = -0.35000D-03
        HTT(20, 2, 4) = -0.26000D-03
        HTT(20, 3, 4) = -0.16000D-03
        HTT(20, 4, 4) = -0.80000D-04
        HTT(20, 5, 4) =  0.40000D-04
        HTT(20, 6, 4) = -0.40000D-03
        HTT(20, 7, 4) = -0.14000D-03
        HTT(20, 8, 4) = -0.57000D-03
        HTT(20, 9, 4) = -0.93000D-03
        HTT(20,10, 4) = -0.42000D-03
        HTT(20,11, 4) = -0.10200D-02
        HTT(20,12, 4) = -0.51000D-03
        HTT(20,13, 4) =  0.26100D-02
        HTT(20,14, 4) = -0.28400D-02
        HTT(20,15, 4) =  0.21100D-02
        HTT(20,16, 4) = -0.33500D-02
        HTT(20,17, 4) = -0.29400D-02
        HTT(20,18, 4) =  0.47000D-03
        HTT(20,19, 4) =  0.25400D-02
        HTT(20,20, 4) =  0.59400D-02
        HTT(21, 1, 4) =  0.21700D-02
        HTT(21, 2, 4) =  0.15000D-02
        HTT(21, 3, 4) = -0.99000D-03
        HTT(21, 4, 4) = -0.16600D-02
        HTT(21, 5, 4) = -0.30200D-02
        HTT(21, 6, 4) = -0.69000D-03
        HTT(21, 7, 4) = -0.20000D-04
        HTT(21, 8, 4) =  0.23100D-02
        HTT(21, 9, 4) = -0.11900D-02
        HTT(21,10, 4) = -0.98000D-03
        HTT(21,11, 4) = -0.50000D-03
        HTT(21,12, 4) = -0.28000D-03
        HTT(21,13, 4) =  0.12600D-02
        HTT(21,14, 4) = -0.76000D-03
        HTT(21,15, 4) =  0.10500D-02
        HTT(21,16, 4) = -0.97000D-03
        HTT(21,17, 4) = -0.21400D-02
        HTT(21,18, 4) = -0.48000D-03
        HTT(21,19, 4) = -0.11000D-03
        HTT(21,20, 4) =  0.15500D-02
        HTT(21,21, 4) =  0.28900D-02
        HTT(22, 1, 4) =  0.27000D-03
        HTT(22, 2, 4) =  0.11300D-02
        HTT(22, 3, 4) = -0.15900D-02
        HTT(22, 4, 4) = -0.74000D-03
        HTT(22, 5, 4) = -0.60000D-04
        HTT(22, 6, 4) = -0.26200D-02
        HTT(22, 7, 4) =  0.17100D-02
        HTT(22, 8, 4) = -0.85000D-03
        HTT(22, 9, 4) = -0.85000D-03
        HTT(22,10, 4) =  0.40000D-04
        HTT(22,11, 4) = -0.17300D-02
        HTT(22,12, 4) = -0.85000D-03
        HTT(22,13, 4) =  0.12800D-02
        HTT(22,14, 4) = -0.15000D-03
        HTT(22,15, 4) =  0.40000D-03
        HTT(22,16, 4) = -0.10300D-02
        HTT(22,17, 4) = -0.10700D-02
        HTT(22,18, 4) =  0.56000D-03
        HTT(22,19, 4) =  0.37000D-03
        HTT(22,20, 4) =  0.20000D-02
        HTT(22,21, 4) =  0.44000D-03
        HTT(22,22, 4) =  0.31200D-02
        HTT(23, 1, 4) =  0.19200D-02
        HTT(23, 2, 4) =  0.79000D-03
        HTT(23, 3, 4) = -0.23000D-03
        HTT(23, 4, 4) = -0.13700D-02
        HTT(23, 5, 4) = -0.21700D-02
        HTT(23, 6, 4) =  0.14000D-03
        HTT(23, 7, 4) = -0.12000D-03
        HTT(23, 8, 4) =  0.21900D-02
        HTT(23, 9, 4) = -0.74000D-03
        HTT(23,10, 4) = -0.17100D-02
        HTT(23,11, 4) =  0.45000D-03
        HTT(23,12, 4) = -0.53000D-03
        HTT(23,13, 4) = -0.10000D-03
        HTT(23,14, 4) = -0.13700D-02
        HTT(23,15, 4) =  0.87000D-03
        HTT(23,16, 4) = -0.40000D-03
        HTT(23,17, 4) = -0.15000D-03
        HTT(23,18, 4) = -0.31500D-02
        HTT(23,19, 4) =  0.11200D-02
        HTT(23,20, 4) = -0.18800D-02
        HTT(23,21, 4) =  0.12200D-02
        HTT(23,22, 4) = -0.12000D-02
        HTT(23,23, 4) =  0.42400D-02
        HTT(24, 1, 4) =  0.30000D-04
        HTT(24, 2, 4) =  0.42000D-03
        HTT(24, 3, 4) = -0.84000D-03
        HTT(24, 4, 4) = -0.45000D-03
        HTT(24, 5, 4) =  0.79000D-03
        HTT(24, 6, 4) = -0.18000D-02
        HTT(24, 7, 4) =  0.16200D-02
        HTT(24, 8, 4) = -0.97000D-03
        HTT(24, 9, 4) = -0.39000D-03
        HTT(24,10, 4) = -0.70000D-03
        HTT(24,11, 4) = -0.79000D-03
        HTT(24,12, 4) = -0.11000D-02
        HTT(24,13, 4) = -0.80000D-04
        HTT(24,14, 4) = -0.76000D-03
        HTT(24,15, 4) =  0.22000D-03
        HTT(24,16, 4) = -0.45000D-03
        HTT(24,17, 4) =  0.92000D-03
        HTT(24,18, 4) = -0.21000D-02
        HTT(24,19, 4) =  0.16000D-02
        HTT(24,20, 4) = -0.14300D-02
        HTT(24,21, 4) = -0.12300D-02
        HTT(24,22, 4) =  0.14800D-02
        HTT(24,23, 4) =  0.18100D-02
        HTT(24,24, 4) =  0.45300D-02
C Hessian elements for anchor point 4
        HBS( 1, 1, 4) = -0.31760D-01
        HBS( 1, 2, 4) = -0.29360D-01
        HBS( 1, 3, 4) = -0.31820D-01
        HBS( 1, 4, 4) = -0.95900D-02
        HBS( 1, 5, 4) =  0.52300D-02
        HBS( 1, 6, 4) =  0.39930D-01
        HBS( 1, 7, 4) = -0.45000D-03
        HBS( 1, 8, 4) =  0.40450D-01
        HBS( 1, 9, 4) = -0.18000D-03
        HBS( 1,10, 4) = -0.96900D-02
        HBS( 1,11, 4) = -0.42000D-03
        HBS( 1,12, 4) =  0.47400D-02
        HBS( 2, 1, 4) =  0.49930D-01
        HBS( 2, 2, 4) = -0.13930D-01
        HBS( 2, 3, 4) =  0.99600D-02
        HBS( 2, 4, 4) =  0.18000D-02
        HBS( 2, 5, 4) = -0.12200D-02
        HBS( 2, 6, 4) = -0.27810D-01
        HBS( 2, 7, 4) =  0.11800D-02
        HBS( 2, 8, 4) = -0.14650D-01
        HBS( 2, 9, 4) =  0.00000D+00
        HBS( 2,10, 4) =  0.41300D-02
        HBS( 2,11, 4) = -0.69000D-03
        HBS( 2,12, 4) = -0.89000D-03
        HBS( 3, 1, 4) = -0.18180D-01
        HBS( 3, 2, 4) =  0.43280D-01
        HBS( 3, 3, 4) =  0.21860D-01
        HBS( 3, 4, 4) =  0.77900D-02
        HBS( 3, 5, 4) = -0.40100D-02
        HBS( 3, 6, 4) = -0.12120D-01
        HBS( 3, 7, 4) = -0.73000D-03
        HBS( 3, 8, 4) = -0.25800D-01
        HBS( 3, 9, 4) =  0.18000D-03
        HBS( 3,10, 4) =  0.55600D-02
        HBS( 3,11, 4) =  0.11100D-02
        HBS( 3,12, 4) = -0.38500D-02
        HBS( 4, 1, 4) = -0.20820D-01
        HBS( 4, 2, 4) = -0.12100D-01
        HBS( 4, 3, 4) =  0.15180D-01
        HBS( 4, 4, 4) = -0.29840D-01
        HBS( 4, 5, 4) = -0.97300D-02
        HBS( 4, 6, 4) = -0.57000D-02
        HBS( 4, 7, 4) =  0.49500D-02
        HBS( 4, 8, 4) =  0.31150D-01
        HBS( 4, 9, 4) = -0.44000D-03
        HBS( 4,10, 4) =  0.40270D-01
        HBS( 4,11, 4) = -0.59000D-03
        HBS( 4,12, 4) =  0.29000D-03
        HBS( 5, 1, 4) =  0.26090D-01
        HBS( 5, 2, 4) =  0.58900D-02
        HBS( 5, 3, 4) = -0.14520D-01
        HBS( 5, 4, 4) =  0.25400D-02
        HBS( 5, 5, 4) =  0.69400D-02
        HBS( 5, 6, 4) =  0.30900D-02
        HBS( 5, 7, 4) = -0.21700D-02
        HBS( 5, 8, 4) = -0.15190D-01
        HBS( 5, 9, 4) = -0.39000D-03
        HBS( 5,10, 4) = -0.22140D-01
        HBS( 5,11, 4) =  0.28000D-03
        HBS( 5,12, 4) =  0.90000D-04
        HBS( 6, 1, 4) = -0.52700D-02
        HBS( 6, 2, 4) =  0.62100D-02
        HBS( 6, 3, 4) = -0.66000D-03
        HBS( 6, 4, 4) =  0.27300D-01
        HBS( 6, 5, 4) =  0.27900D-02
        HBS( 6, 6, 4) =  0.26100D-02
        HBS( 6, 7, 4) = -0.27800D-02
        HBS( 6, 8, 4) = -0.15960D-01
        HBS( 6, 9, 4) =  0.83000D-03
        HBS( 6,10, 4) = -0.18130D-01
        HBS( 6,11, 4) =  0.31000D-03
        HBS( 6,12, 4) = -0.38000D-03
        HBS( 7, 1, 4) = -0.90800D-02
        HBS( 7, 2, 4) =  0.35870D-01
        HBS( 7, 3, 4) = -0.10000D-04
        HBS( 7, 4, 4) = -0.22920D-01
        HBS( 7, 5, 4) =  0.61900D-02
        HBS( 7, 6, 4) = -0.28970D-01
        HBS( 7, 7, 4) = -0.87500D-02
        HBS( 7, 8, 4) = -0.92900D-02
        HBS( 7, 9, 4) =  0.49700D-02
        HBS( 7,10, 4) =  0.32780D-01
        HBS( 7,11, 4) =  0.60000D-04
        HBS( 7,12, 4) = -0.10700D-02
        HBS( 8, 1, 4) =  0.36900D-02
        HBS( 8, 2, 4) = -0.18590D-01
        HBS( 8, 3, 4) =  0.13800D-02
        HBS( 8, 4, 4) =  0.25350D-01
        HBS( 8, 5, 4) = -0.35000D-02
        HBS( 8, 6, 4) = -0.42000D-03
        HBS( 8, 7, 4) =  0.45000D-02
        HBS( 8, 8, 4) =  0.50400D-02
        HBS( 8, 9, 4) = -0.20100D-02
        HBS( 8,10, 4) = -0.14680D-01
        HBS( 8,11, 4) = -0.62000D-03
        HBS( 8,12, 4) =  0.55000D-03
        HBS( 9, 1, 4) =  0.53900D-02
        HBS( 9, 2, 4) = -0.17270D-01
        HBS( 9, 3, 4) = -0.13700D-02
        HBS( 9, 4, 4) = -0.24400D-02
        HBS( 9, 5, 4) = -0.26800D-02
        HBS( 9, 6, 4) =  0.29390D-01
        HBS( 9, 7, 4) =  0.42400D-02
        HBS( 9, 8, 4) =  0.42500D-02
        HBS( 9, 9, 4) = -0.29600D-02
        HBS( 9,10, 4) = -0.18100D-01
        HBS( 9,11, 4) =  0.55000D-03
        HBS( 9,12, 4) =  0.51000D-03
        HBS(10, 1, 4) =  0.36130D-01
        HBS(10, 2, 4) =  0.38360D-01
        HBS(10, 3, 4) =  0.36300D-02
        HBS(10, 4, 4) = -0.11890D-01
        HBS(10, 5, 4) = -0.10500D-02
        HBS(10, 6, 4) = -0.26840D-01
        HBS(10, 7, 4) =  0.48100D-02
        HBS(10, 8, 4) = -0.27100D-01
        HBS(10, 9, 4) = -0.88200D-02
        HBS(10,10, 4) = -0.12140D-01
        HBS(10,11, 4) =  0.47600D-02
        HBS(10,12, 4) = -0.77000D-03
        HBS(11, 1, 4) = -0.20190D-01
        HBS(11, 2, 4) = -0.17180D-01
        HBS(11, 3, 4) = -0.21300D-02
        HBS(11, 4, 4) =  0.57400D-02
        HBS(11, 5, 4) =  0.12900D-02
        HBS(11, 6, 4) =  0.27000D-01
        HBS(11, 7, 4) = -0.28400D-02
        HBS(11, 8, 4) =  0.47000D-03
        HBS(11, 9, 4) =  0.44300D-02
        HBS(11,10, 4) =  0.60700D-02
        HBS(11,11, 4) = -0.19800D-02
        HBS(11,12, 4) = -0.18000D-03
        HBS(12, 1, 4) = -0.15940D-01
        HBS(12, 2, 4) = -0.21180D-01
        HBS(12, 3, 4) = -0.14900D-02
        HBS(12, 4, 4) =  0.61500D-02
        HBS(12, 5, 4) = -0.24000D-03
        HBS(12, 6, 4) = -0.16000D-03
        HBS(12, 7, 4) = -0.19600D-02
        HBS(12, 8, 4) =  0.26630D-01
        HBS(12, 9, 4) =  0.43900D-02
        HBS(12,10, 4) =  0.60700D-02
        HBS(12,11, 4) = -0.27800D-02
        HBS(12,12, 4) =  0.94000D-03
        HBS(13, 1, 4) =  0.35720D-01
        HBS(13, 2, 4) = -0.96900D-02
        HBS(13, 3, 4) =  0.29800D-02
        HBS(13, 4, 4) =  0.32310D-01
        HBS(13, 5, 4) = -0.94000D-03
        HBS(13, 6, 4) = -0.10140D-01
        HBS(13, 7, 4) =  0.70000D-04
        HBS(13, 8, 4) = -0.29530D-01
        HBS(13, 9, 4) =  0.50200D-02
        HBS(13,10, 4) = -0.23420D-01
        HBS(13,11, 4) = -0.86600D-02
        HBS(13,12, 4) =  0.53400D-02
        HBS(14, 1, 4) = -0.17210D-01
        HBS(14, 2, 4) =  0.56100D-02
        HBS(14, 3, 4) = -0.23300D-02
        HBS(14, 4, 4) = -0.17680D-01
        HBS(14, 5, 4) =  0.48000D-03
        HBS(14, 6, 4) =  0.47700D-02
        HBS(14, 7, 4) =  0.56000D-03
        HBS(14, 8, 4) =  0.29630D-01
        HBS(14, 9, 4) = -0.29700D-02
        HBS(14,10, 4) = -0.24000D-02
        HBS(14,11, 4) =  0.43000D-02
        HBS(14,12, 4) = -0.22600D-02
        HBS(15, 1, 4) = -0.18510D-01
        HBS(15, 2, 4) =  0.40800D-02
        HBS(15, 3, 4) = -0.65000D-03
        HBS(15, 4, 4) = -0.14630D-01
        HBS(15, 5, 4) =  0.46000D-03
        HBS(15, 6, 4) =  0.53700D-02
        HBS(15, 7, 4) = -0.63000D-03
        HBS(15, 8, 4) = -0.10000D-03
        HBS(15, 9, 4) = -0.20500D-02
        HBS(15,10, 4) =  0.25820D-01
        HBS(15,11, 4) =  0.43700D-02
        HBS(15,12, 4) = -0.30900D-02
        HBS(16, 1, 4) = -0.10200D-01
        HBS(16, 2, 4) = -0.23080D-01
        HBS(16, 3, 4) =  0.10040D-01
        HBS(16, 4, 4) =  0.41930D-01
        HBS(16, 5, 4) =  0.31000D-03
        HBS(16, 6, 4) =  0.31710D-01
        HBS(16, 7, 4) = -0.63000D-03
        HBS(16, 8, 4) = -0.56800D-02
        HBS(16, 9, 4) = -0.55000D-03
        HBS(16,10, 4) = -0.27810D-01
        HBS(16,11, 4) =  0.48500D-02
        HBS(16,12, 4) = -0.85400D-02
        HBS(17, 1, 4) =  0.41100D-02
        HBS(17, 2, 4) =  0.25000D-01
        HBS(17, 3, 4) = -0.74400D-02
        HBS(17, 4, 4) = -0.22980D-01
        HBS(17, 5, 4) =  0.60000D-03
        HBS(17, 6, 4) = -0.14860D-01
        HBS(17, 7, 4) =  0.32000D-03
        HBS(17, 8, 4) =  0.33600D-02
        HBS(17, 9, 4) = -0.32000D-03
        HBS(17,10, 4) =  0.20400D-02
        HBS(17,11, 4) = -0.21000D-02
        HBS(17,12, 4) =  0.68800D-02
        HBS(18, 1, 4) =  0.61000D-02
        HBS(18, 2, 4) = -0.19200D-02
        HBS(18, 3, 4) = -0.26000D-02
        HBS(18, 4, 4) = -0.18950D-01
        HBS(18, 5, 4) = -0.92000D-03
        HBS(18, 6, 4) = -0.16850D-01
        HBS(18, 7, 4) =  0.31000D-03
        HBS(18, 8, 4) =  0.23200D-02
        HBS(18, 9, 4) =  0.87000D-03
        HBS(18,10, 4) =  0.25770D-01
        HBS(18,11, 4) = -0.27500D-02
        HBS(18,12, 4) =  0.16600D-02
      else if(ii .eq. 33) then
C Hessian elements for anchor point 1
        HSS( 1, 1, 1) =  0.32935D+00
        HSS( 2, 1, 1) =  0.10630D-01
        HSS( 2, 2, 1) =  0.33545D+00
        HSS( 3, 1, 1) =  0.41480D-01
        HSS( 3, 2, 1) =  0.47860D-01
        HSS( 3, 3, 1) =  0.60075D+00
        HSS( 4, 1, 1) = -0.28200D-02
        HSS( 4, 2, 1) =  0.17340D-01
        HSS( 4, 3, 1) = -0.51100D-02
        HSS( 4, 4, 1) =  0.44790D+00
        HSS( 5, 1, 1) =  0.39700D-02
        HSS( 5, 2, 1) = -0.41000D-03
        HSS( 5, 3, 1) =  0.20000D-03
        HSS( 5, 4, 1) =  0.13100D-02
        HSS( 5, 5, 1) =  0.36210D+00
        HSS( 6, 1, 1) =  0.23160D-01
        HSS( 6, 2, 1) =  0.46150D-01
        HSS( 6, 3, 1) =  0.31100D-02
        HSS( 6, 4, 1) =  0.42160D-01
        HSS( 6, 5, 1) = -0.19000D-03
        HSS( 6, 6, 1) =  0.33957D+00
        HSS( 7, 1, 1) =  0.14000D-03
        HSS( 7, 2, 1) = -0.96000D-03
        HSS( 7, 3, 1) =  0.11700D-02
        HSS( 7, 4, 1) =  0.55900D-02
        HSS( 7, 5, 1) = -0.68000D-03
        HSS( 7, 6, 1) =  0.28700D-02
        HSS( 7, 7, 1) =  0.36718D+00
        HSS( 8, 1, 1) =  0.48190D-01
        HSS( 8, 2, 1) =  0.21300D-01
        HSS( 8, 3, 1) =  0.10200D-02
        HSS( 8, 4, 1) = -0.55600D-02
        HSS( 8, 5, 1) = -0.21100D-02
        HSS( 8, 6, 1) =  0.34440D-01
        HSS( 8, 7, 1) = -0.39000D-03
        HSS( 8, 8, 1) =  0.33723D+00
        HSS( 9, 1, 1) = -0.13000D-03
        HSS( 9, 2, 1) =  0.50000D-04
        HSS( 9, 3, 1) = -0.53000D-03
        HSS( 9, 4, 1) = -0.54000D-03
        HSS( 9, 5, 1) = -0.64000D-03
        HSS( 9, 6, 1) =  0.37800D-02
        HSS( 9, 7, 1) =  0.12000D-03
        HSS( 9, 8, 1) =  0.36400D-02
        HSS( 9, 9, 1) =  0.36496D+00
        HSS(10, 1, 1) =  0.14800D-01
        HSS(10, 2, 1) =  0.14000D-03
        HSS(10, 3, 1) = -0.42300D-02
        HSS(10, 4, 1) =  0.81590D-01
        HSS(10, 5, 1) = -0.48000D-03
        HSS(10, 6, 1) = -0.32700D-02
        HSS(10, 7, 1) = -0.71000D-03
        HSS(10, 8, 1) =  0.36040D-01
        HSS(10, 9, 1) = -0.87000D-03
        HSS(10,10, 1) =  0.46416D+00
        HSS(11, 1, 1) = -0.11800D-02
        HSS(11, 2, 1) = -0.64000D-03
        HSS(11, 3, 1) =  0.90000D-03
        HSS(11, 4, 1) = -0.78000D-03
        HSS(11, 5, 1) =  0.18000D-03
        HSS(11, 6, 1) = -0.11000D-03
        HSS(11, 7, 1) =  0.12000D-03
        HSS(11, 8, 1) =  0.28600D-02
        HSS(11, 9, 1) =  0.42000D-03
        HSS(11,10, 1) =  0.45200D-02
        HSS(11,11, 1) =  0.36758D+00
        HSS(12, 1, 1) =  0.23000D-03
        HSS(12, 2, 1) =  0.14400D-02
        HSS(12, 3, 1) =  0.13000D-03
        HSS(12, 4, 1) = -0.34000D-03
        HSS(12, 5, 1) = -0.37000D-03
        HSS(12, 6, 1) = -0.80000D-03
        HSS(12, 7, 1) =  0.40000D-04
        HSS(12, 8, 1) =  0.22000D-03
        HSS(12, 9, 1) = -0.10000D-03
        HSS(12,10, 1) =  0.23100D-02
        HSS(12,11, 1) =  0.17000D-03
        HSS(12,12, 1) =  0.37062D+00
C Hessian elements for anchor point 1
        HBB( 1, 1, 1) =  0.82090D-01
        HBB( 2, 1, 1) = -0.40800D-01
        HBB( 2, 2, 1) =  0.13448D+00
        HBB( 3, 1, 1) = -0.41290D-01
        HBB( 3, 2, 1) = -0.93680D-01
        HBB( 3, 3, 1) =  0.13497D+00
        HBB( 4, 1, 1) = -0.37560D-01
        HBB( 4, 2, 1) =  0.25500D-01
        HBB( 4, 3, 1) =  0.12070D-01
        HBB( 4, 4, 1) =  0.76550D-01
        HBB( 5, 1, 1) =  0.22450D-01
        HBB( 5, 2, 1) = -0.15660D-01
        HBB( 5, 3, 1) = -0.67900D-02
        HBB( 5, 4, 1) = -0.39660D-01
        HBB( 5, 5, 1) =  0.77760D-01
        HBB( 6, 1, 1) =  0.15120D-01
        HBB( 6, 2, 1) = -0.98400D-02
        HBB( 6, 3, 1) = -0.52700D-02
        HBB( 6, 4, 1) = -0.36880D-01
        HBB( 6, 5, 1) = -0.38100D-01
        HBB( 6, 6, 1) =  0.74990D-01
        HBB( 7, 1, 1) = -0.83500D-02
        HBB( 7, 2, 1) = -0.60600D-02
        HBB( 7, 3, 1) =  0.14410D-01
        HBB( 7, 4, 1) = -0.28860D-01
        HBB( 7, 5, 1) =  0.12030D-01
        HBB( 7, 6, 1) =  0.16840D-01
        HBB( 7, 7, 1) =  0.74530D-01
        HBB( 8, 1, 1) = -0.15000D-02
        HBB( 8, 2, 1) =  0.75600D-02
        HBB( 8, 3, 1) = -0.60500D-02
        HBB( 8, 4, 1) =  0.17370D-01
        HBB( 8, 5, 1) = -0.66700D-02
        HBB( 8, 6, 1) = -0.10710D-01
        HBB( 8, 7, 1) = -0.36970D-01
        HBB( 8, 8, 1) =  0.78530D-01
        HBB( 9, 1, 1) =  0.98600D-02
        HBB( 9, 2, 1) = -0.15000D-02
        HBB( 9, 3, 1) = -0.83500D-02
        HBB( 9, 4, 1) =  0.11490D-01
        HBB( 9, 5, 1) = -0.53600D-02
        HBB( 9, 6, 1) = -0.61300D-02
        HBB( 9, 7, 1) = -0.37550D-01
        HBB( 9, 8, 1) = -0.41560D-01
        HBB( 9, 9, 1) =  0.79110D-01
        HBB(10, 1, 1) =  0.88000D-02
        HBB(10, 2, 1) = -0.45300D-02
        HBB(10, 3, 1) = -0.42700D-02
        HBB(10, 4, 1) = -0.86000D-02
        HBB(10, 5, 1) =  0.10210D-01
        HBB(10, 6, 1) = -0.16100D-02
        HBB(10, 7, 1) = -0.33520D-01
        HBB(10, 8, 1) =  0.13620D-01
        HBB(10, 9, 1) =  0.19900D-01
        HBB(10,10, 1) =  0.75620D-01
        HBB(11, 1, 1) = -0.44200D-02
        HBB(11, 2, 1) =  0.26700D-02
        HBB(11, 3, 1) =  0.17500D-02
        HBB(11, 4, 1) = -0.13500D-02
        HBB(11, 5, 1) = -0.34500D-02
        HBB(11, 6, 1) =  0.48000D-02
        HBB(11, 7, 1) =  0.19740D-01
        HBB(11, 8, 1) = -0.73600D-02
        HBB(11, 9, 1) = -0.12380D-01
        HBB(11,10, 1) = -0.37740D-01
        HBB(11,11, 1) =  0.80040D-01
        HBB(12, 1, 1) = -0.43800D-02
        HBB(12, 2, 1) =  0.18600D-02
        HBB(12, 3, 1) =  0.25200D-02
        HBB(12, 4, 1) =  0.99500D-02
        HBB(12, 5, 1) = -0.67600D-02
        HBB(12, 6, 1) = -0.31900D-02
        HBB(12, 7, 1) =  0.13790D-01
        HBB(12, 8, 1) = -0.62600D-02
        HBB(12, 9, 1) = -0.75300D-02
        HBB(12,10, 1) = -0.37880D-01
        HBB(12,11, 1) = -0.42300D-01
        HBB(12,12, 1) =  0.80170D-01
        HBB(13, 1, 1) = -0.76900D-02
        HBB(13, 2, 1) =  0.14150D-01
        HBB(13, 3, 1) = -0.64500D-02
        HBB(13, 4, 1) =  0.32300D-02
        HBB(13, 5, 1) = -0.14400D-02
        HBB(13, 6, 1) = -0.17900D-02
        HBB(13, 7, 1) = -0.67100D-02
        HBB(13, 8, 1) =  0.90400D-02
        HBB(13, 9, 1) = -0.23400D-02
        HBB(13,10, 1) = -0.33330D-01
        HBB(13,11, 1) =  0.13660D-01
        HBB(13,12, 1) =  0.19670D-01
        HBB(13,13, 1) =  0.73530D-01
        HBB(14, 1, 1) =  0.96900D-02
        HBB(14, 2, 1) = -0.81500D-02
        HBB(14, 3, 1) = -0.15500D-02
        HBB(14, 4, 1) = -0.17500D-02
        HBB(14, 5, 1) =  0.64000D-03
        HBB(14, 6, 1) =  0.11000D-02
        HBB(14, 7, 1) = -0.22600D-02
        HBB(14, 8, 1) = -0.29800D-02
        HBB(14, 9, 1) =  0.52500D-02
        HBB(14,10, 1) =  0.19960D-01
        HBB(14,11, 1) = -0.75600D-02
        HBB(14,12, 1) = -0.12400D-01
        HBB(14,13, 1) = -0.37280D-01
        HBB(14,14, 1) =  0.78950D-01
        HBB(15, 1, 1) = -0.20000D-02
        HBB(15, 2, 1) = -0.60000D-02
        HBB(15, 3, 1) =  0.80000D-02
        HBB(15, 4, 1) = -0.14800D-02
        HBB(15, 5, 1) =  0.80000D-03
        HBB(15, 6, 1) =  0.69000D-03
        HBB(15, 7, 1) =  0.89700D-02
        HBB(15, 8, 1) = -0.60600D-02
        HBB(15, 9, 1) = -0.29100D-02
        HBB(15,10, 1) =  0.13370D-01
        HBB(15,11, 1) = -0.61000D-02
        HBB(15,12, 1) = -0.72700D-02
        HBB(15,13, 1) = -0.36250D-01
        HBB(15,14, 1) = -0.41670D-01
        HBB(15,15, 1) =  0.77920D-01
        HBB(16, 1, 1) = -0.37280D-01
        HBB(16, 2, 1) =  0.11740D-01
        HBB(16, 3, 1) =  0.25540D-01
        HBB(16, 4, 1) = -0.47500D-02
        HBB(16, 5, 1) = -0.35800D-02
        HBB(16, 6, 1) =  0.83300D-02
        HBB(16, 7, 1) =  0.29200D-02
        HBB(16, 8, 1) = -0.15600D-02
        HBB(16, 9, 1) = -0.13600D-02
        HBB(16,10, 1) = -0.89600D-02
        HBB(16,11, 1) =  0.10110D-01
        HBB(16,12, 1) = -0.11500D-02
        HBB(16,13, 1) = -0.29030D-01
        HBB(16,14, 1) =  0.11640D-01
        HBB(16,15, 1) =  0.17390D-01
        HBB(16,16, 1) =  0.77110D-01
        HBB(17, 1, 1) =  0.22780D-01
        HBB(17, 2, 1) = -0.48100D-02
        HBB(17, 3, 1) = -0.17970D-01
        HBB(17, 4, 1) = -0.38000D-02
        HBB(17, 5, 1) =  0.60900D-02
        HBB(17, 6, 1) = -0.22900D-02
        HBB(17, 7, 1) = -0.14200D-02
        HBB(17, 8, 1) =  0.93000D-03
        HBB(17, 9, 1) =  0.49000D-03
        HBB(17,10, 1) =  0.10430D-01
        HBB(17,11, 1) = -0.67600D-02
        HBB(17,12, 1) = -0.36800D-02
        HBB(17,13, 1) =  0.11310D-01
        HBB(17,14, 1) = -0.50800D-02
        HBB(17,15, 1) = -0.62300D-02
        HBB(17,16, 1) = -0.39310D-01
        HBB(17,17, 1) =  0.76690D-01
        HBB(18, 1, 1) =  0.14500D-01
        HBB(18, 2, 1) = -0.69300D-02
        HBB(18, 3, 1) = -0.75700D-02
        HBB(18, 4, 1) =  0.85400D-02
        HBB(18, 5, 1) = -0.25100D-02
        HBB(18, 6, 1) = -0.60300D-02
        HBB(18, 7, 1) = -0.15000D-02
        HBB(18, 8, 1) =  0.63000D-03
        HBB(18, 9, 1) =  0.87000D-03
        HBB(18,10, 1) = -0.14700D-02
        HBB(18,11, 1) = -0.33600D-02
        HBB(18,12, 1) =  0.48300D-02
        HBB(18,13, 1) =  0.17720D-01
        HBB(18,14, 1) = -0.65600D-02
        HBB(18,15, 1) = -0.11160D-01
        HBB(18,16, 1) = -0.37800D-01
        HBB(18,17, 1) = -0.37380D-01
        HBB(18,18, 1) =  0.75170D-01
C Hessian elements for anchor point 1
        HTT( 1, 1, 1) =  0.12410D-01
        HTT( 2, 1, 1) =  0.66300D-02
        HTT( 2, 2, 1) =  0.11350D-01
        HTT( 3, 1, 1) = -0.60000D-02
        HTT( 3, 2, 1) = -0.92200D-02
        HTT( 3, 3, 1) =  0.16800D-01
        HTT( 4, 1, 1) = -0.11780D-01
        HTT( 4, 2, 1) = -0.45100D-02
        HTT( 4, 3, 1) =  0.13580D-01
        HTT( 4, 4, 1) =  0.20850D-01
        HTT( 5, 1, 1) = -0.10450D-01
        HTT( 5, 2, 1) = -0.82200D-02
        HTT( 5, 3, 1) =  0.95200D-02
        HTT( 5, 4, 1) =  0.11750D-01
        HTT( 5, 5, 1) =  0.13600D-01
        HTT( 6, 1, 1) = -0.73200D-02
        HTT( 6, 2, 1) = -0.56100D-02
        HTT( 6, 3, 1) =  0.79400D-02
        HTT( 6, 4, 1) =  0.96500D-02
        HTT( 6, 5, 1) =  0.70500D-02
        HTT( 6, 6, 1) =  0.10410D-01
        HTT( 7, 1, 1) =  0.68700D-02
        HTT( 7, 2, 1) =  0.66900D-02
        HTT( 7, 3, 1) = -0.11920D-01
        HTT( 7, 4, 1) = -0.12100D-01
        HTT( 7, 5, 1) = -0.51800D-02
        HTT( 7, 6, 1) = -0.73000D-02
        HTT( 7, 7, 1) =  0.12490D-01
        HTT( 8, 1, 1) =  0.10000D-01
        HTT( 8, 2, 1) =  0.93000D-02
        HTT( 8, 3, 1) = -0.13500D-01
        HTT( 8, 4, 1) = -0.14200D-01
        HTT( 8, 5, 1) = -0.11730D-01
        HTT( 8, 6, 1) = -0.39400D-02
        HTT( 8, 7, 1) =  0.10370D-01
        HTT( 8, 8, 1) =  0.18160D-01
        HTT( 9, 1, 1) = -0.62000D-02
        HTT( 9, 2, 1) = -0.56000D-03
        HTT( 9, 3, 1) = -0.18500D-02
        HTT( 9, 4, 1) =  0.37900D-02
        HTT( 9, 5, 1) =  0.12100D-02
        HTT( 9, 6, 1) =  0.24300D-02
        HTT( 9, 7, 1) = -0.28800D-02
        HTT( 9, 8, 1) = -0.16600D-02
        HTT( 9, 9, 1) =  0.82000D-02
        HTT(10, 1, 1) = -0.49900D-02
        HTT(10, 2, 1) = -0.30000D-04
        HTT(10, 3, 1) =  0.60000D-03
        HTT(10, 4, 1) =  0.55600D-02
        HTT(10, 5, 1) =  0.33700D-02
        HTT(10, 6, 1) =  0.24600D-02
        HTT(10, 7, 1) = -0.18800D-02
        HTT(10, 8, 1) = -0.27900D-02
        HTT(10, 9, 1) =  0.31400D-02
        HTT(10,10, 1) =  0.89300D-02
        HTT(11, 1, 1) = -0.90000D-04
        HTT(11, 2, 1) = -0.55400D-02
        HTT(11, 3, 1) =  0.15600D-02
        HTT(11, 4, 1) = -0.39000D-02
        HTT(11, 5, 1) = -0.11400D-02
        HTT(11, 6, 1) =  0.63000D-03
        HTT(11, 7, 1) = -0.26900D-02
        HTT(11, 8, 1) = -0.91000D-03
        HTT(11, 9, 1) =  0.22400D-02
        HTT(11,10, 1) = -0.21000D-02
        HTT(11,11, 1) =  0.80000D-02
        HTT(12, 1, 1) =  0.11300D-02
        HTT(12, 2, 1) = -0.50100D-02
        HTT(12, 3, 1) =  0.40000D-02
        HTT(12, 4, 1) = -0.21300D-02
        HTT(12, 5, 1) =  0.10100D-02
        HTT(12, 6, 1) =  0.67000D-03
        HTT(12, 7, 1) = -0.17000D-02
        HTT(12, 8, 1) = -0.20400D-02
        HTT(12, 9, 1) = -0.28200D-02
        HTT(12,10, 1) =  0.36900D-02
        HTT(12,11, 1) =  0.36700D-02
        HTT(12,12, 1) =  0.10180D-01
        HTT(13, 1, 1) = -0.19600D-02
        HTT(13, 2, 1) = -0.40500D-02
        HTT(13, 3, 1) =  0.64000D-02
        HTT(13, 4, 1) =  0.43100D-02
        HTT(13, 5, 1) =  0.50600D-02
        HTT(13, 6, 1) =  0.27500D-02
        HTT(13, 7, 1) = -0.28000D-02
        HTT(13, 8, 1) = -0.51000D-02
        HTT(13, 9, 1) = -0.54600D-02
        HTT(13,10, 1) =  0.27000D-03
        HTT(13,11, 1) = -0.32400D-02
        HTT(13,12, 1) =  0.24800D-02
        HTT(13,13, 1) =  0.10020D-01
        HTT(14, 1, 1) =  0.73000D-03
        HTT(14, 2, 1) = -0.74000D-03
        HTT(14, 3, 1) =  0.28700D-02
        HTT(14, 4, 1) =  0.14100D-02
        HTT(14, 5, 1) =  0.19900D-02
        HTT(14, 6, 1) = -0.24000D-03
        HTT(14, 7, 1) = -0.30000D-04
        HTT(14, 8, 1) = -0.22600D-02
        HTT(14, 9, 1) = -0.28600D-02
        HTT(14,10, 1) =  0.17000D-02
        HTT(14,11, 1) = -0.13100D-02
        HTT(14,12, 1) =  0.32500D-02
        HTT(14,13, 1) =  0.23900D-02
        HTT(14,14, 1) =  0.74100D-02
        HTT(15, 1, 1) = -0.31700D-02
        HTT(15, 2, 1) = -0.45800D-02
        HTT(15, 3, 1) =  0.39600D-02
        HTT(15, 4, 1) =  0.25500D-02
        HTT(15, 5, 1) =  0.29100D-02
        HTT(15, 6, 1) =  0.27200D-02
        HTT(15, 7, 1) = -0.37900D-02
        HTT(15, 8, 1) = -0.39800D-02
        HTT(15, 9, 1) = -0.41000D-03
        HTT(15,10, 1) = -0.55000D-02
        HTT(15,11, 1) =  0.10800D-02
        HTT(15,12, 1) = -0.40100D-02
        HTT(15,13, 1) =  0.43100D-02
        HTT(15,14, 1) = -0.21600D-02
        HTT(15,15, 1) =  0.93900D-02
        HTT(16, 1, 1) = -0.48000D-03
        HTT(16, 2, 1) = -0.12700D-02
        HTT(16, 3, 1) =  0.43000D-03
        HTT(16, 4, 1) = -0.35000D-03
        HTT(16, 5, 1) = -0.16000D-03
        HTT(16, 6, 1) = -0.28000D-03
        HTT(16, 7, 1) = -0.10200D-02
        HTT(16, 8, 1) = -0.11300D-02
        HTT(16, 9, 1) =  0.21900D-02
        HTT(16,10, 1) = -0.40700D-02
        HTT(16,11, 1) =  0.30200D-02
        HTT(16,12, 1) = -0.32400D-02
        HTT(16,13, 1) = -0.33300D-02
        HTT(16,14, 1) =  0.28700D-02
        HTT(16,15, 1) =  0.29200D-02
        HTT(16,16, 1) =  0.91100D-02
        HTT(17, 1, 1) =  0.40200D-02
        HTT(17, 2, 1) =  0.24400D-02
        HTT(17, 3, 1) = -0.28100D-02
        HTT(17, 4, 1) = -0.43900D-02
        HTT(17, 5, 1) = -0.18700D-02
        HTT(17, 6, 1) = -0.30500D-02
        HTT(17, 7, 1) =  0.45500D-02
        HTT(17, 8, 1) =  0.33700D-02
        HTT(17, 9, 1) =  0.30000D-03
        HTT(17,10, 1) = -0.19500D-02
        HTT(17,11, 1) =  0.19700D-02
        HTT(17,12, 1) = -0.28000D-03
        HTT(17,13, 1) = -0.68200D-02
        HTT(17,14, 1) =  0.40000D-03
        HTT(17,15, 1) = -0.45700D-02
        HTT(17,16, 1) =  0.26500D-02
        HTT(17,17, 1) =  0.90000D-02
        HTT(18, 1, 1) =  0.19000D-02
        HTT(18, 2, 1) =  0.26000D-02
        HTT(18, 3, 1) = -0.39700D-02
        HTT(18, 4, 1) = -0.32700D-02
        HTT(18, 5, 1) = -0.30200D-02
        HTT(18, 6, 1) = -0.37000D-02
        HTT(18, 7, 1) =  0.25000D-02
        HTT(18, 8, 1) =  0.18200D-02
        HTT(18, 9, 1) =  0.19300D-02
        HTT(18,10, 1) =  0.43000D-03
        HTT(18,11, 1) =  0.11900D-02
        HTT(18,12, 1) = -0.31000D-03
        HTT(18,13, 1) = -0.47600D-02
        HTT(18,14, 1) =  0.15300D-02
        HTT(18,15, 1) = -0.32600D-02
        HTT(18,16, 1) =  0.30200D-02
        HTT(18,17, 1) =  0.36100D-02
        HTT(18,18, 1) =  0.88300D-02
        HTT(19, 1, 1) =  0.13400D-02
        HTT(19, 2, 1) = -0.86000D-03
        HTT(19, 3, 1) =  0.71000D-03
        HTT(19, 4, 1) = -0.14900D-02
        HTT(19, 5, 1) =  0.11900D-02
        HTT(19, 6, 1) = -0.60000D-04
        HTT(19, 7, 1) =  0.17900D-02
        HTT(19, 8, 1) =  0.54000D-03
        HTT(19, 9, 1) = -0.23000D-02
        HTT(19,10, 1) = -0.33800D-02
        HTT(19,11, 1) =  0.40000D-04
        HTT(19,12, 1) = -0.10500D-02
        HTT(19,13, 1) =  0.81000D-03
        HTT(19,14, 1) = -0.46200D-02
        HTT(19,15, 1) =  0.18900D-02
        HTT(19,16, 1) = -0.35300D-02
        HTT(19,17, 1) =  0.17900D-02
        HTT(19,18, 1) = -0.26700D-02
        HTT(19,19, 1) =  0.72000D-02
        HTT(20, 1, 1) = -0.78000D-03
        HTT(20, 2, 1) = -0.70000D-03
        HTT(20, 3, 1) = -0.45000D-03
        HTT(20, 4, 1) = -0.37000D-03
        HTT(20, 5, 1) =  0.40000D-04
        HTT(20, 6, 1) = -0.71000D-03
        HTT(20, 7, 1) = -0.27000D-03
        HTT(20, 8, 1) = -0.10200D-02
        HTT(20, 9, 1) = -0.66000D-03
        HTT(20,10, 1) = -0.10000D-02
        HTT(20,11, 1) = -0.74000D-03
        HTT(20,12, 1) = -0.10700D-02
        HTT(20,13, 1) =  0.28600D-02
        HTT(20,14, 1) = -0.34900D-02
        HTT(20,15, 1) =  0.32000D-02
        HTT(20,16, 1) = -0.31600D-02
        HTT(20,17, 1) = -0.36000D-02
        HTT(20,18, 1) =  0.25500D-02
        HTT(20,19, 1) =  0.27400D-02
        HTT(20,20, 1) =  0.89000D-02
        HTT(21, 1, 1) =  0.22500D-02
        HTT(21, 2, 1) =  0.36500D-02
        HTT(21, 3, 1) = -0.50600D-02
        HTT(21, 4, 1) = -0.36600D-02
        HTT(21, 5, 1) = -0.73600D-02
        HTT(21, 6, 1) = -0.18700D-02
        HTT(21, 7, 1) = -0.49000D-03
        HTT(21, 8, 1) =  0.50000D-02
        HTT(21, 9, 1) =  0.16800D-02
        HTT(21,10, 1) =  0.80000D-04
        HTT(21,11, 1) =  0.20000D-03
        HTT(21,12, 1) = -0.14000D-02
        HTT(21,13, 1) = -0.68000D-03
        HTT(21,14, 1) = -0.25100D-02
        HTT(21,15, 1) =  0.92000D-03
        HTT(21,16, 1) = -0.91000D-03
        HTT(21,17, 1) = -0.45300D-02
        HTT(21,18, 1) =  0.28000D-03
        HTT(21,19, 1) = -0.27000D-02
        HTT(21,20, 1) =  0.21000D-02
        HTT(21,21, 1) =  0.83700D-02
        HTT(22, 1, 1) = -0.10700D-02
        HTT(22, 2, 1) =  0.89000D-03
        HTT(22, 3, 1) = -0.33900D-02
        HTT(22, 4, 1) = -0.14300D-02
        HTT(22, 5, 1) = -0.43000D-03
        HTT(22, 6, 1) = -0.54300D-02
        HTT(22, 7, 1) =  0.17600D-02
        HTT(22, 8, 1) = -0.32400D-02
        HTT(22, 9, 1) =  0.39000D-03
        HTT(22,10, 1) =  0.10300D-02
        HTT(22,11, 1) = -0.16800D-02
        HTT(22,12, 1) = -0.10400D-02
        HTT(22,13, 1) =  0.17600D-02
        HTT(22,14, 1) = -0.15000D-03
        HTT(22,15, 1) =  0.11200D-02
        HTT(22,16, 1) = -0.78000D-03
        HTT(22,17, 1) = -0.32800D-02
        HTT(22,18, 1) =  0.99000D-03
        HTT(22,19, 1) = -0.13800D-02
        HTT(22,20, 1) =  0.29000D-02
        HTT(22,21, 1) =  0.25500D-02
        HTT(22,22, 1) =  0.78500D-02
        HTT(23, 1, 1) =  0.43900D-02
        HTT(23, 2, 1) =  0.34900D-02
        HTT(23, 3, 1) = -0.38900D-02
        HTT(23, 4, 1) = -0.47900D-02
        HTT(23, 5, 1) = -0.62100D-02
        HTT(23, 6, 1) = -0.12200D-02
        HTT(23, 7, 1) =  0.15800D-02
        HTT(23, 8, 1) =  0.65700D-02
        HTT(23, 9, 1) =  0.40000D-04
        HTT(23,10, 1) = -0.23200D-02
        HTT(23,11, 1) =  0.98000D-03
        HTT(23,12, 1) = -0.13800D-02
        HTT(23,13, 1) = -0.27500D-02
        HTT(23,14, 1) = -0.36400D-02
        HTT(23,15, 1) = -0.39000D-03
        HTT(23,16, 1) = -0.12800D-02
        HTT(23,17, 1) =  0.91000D-03
        HTT(23,18, 1) = -0.49900D-02
        HTT(23,19, 1) =  0.18000D-02
        HTT(23,20, 1) = -0.41000D-02
        HTT(23,21, 1) =  0.35200D-02
        HTT(23,22, 1) = -0.17600D-02
        HTT(23,23, 1) =  0.94700D-02
        HTT(24, 1, 1) =  0.10700D-02
        HTT(24, 2, 1) =  0.73000D-03
        HTT(24, 3, 1) = -0.22200D-02
        HTT(24, 4, 1) = -0.25600D-02
        HTT(24, 5, 1) =  0.73000D-03
        HTT(24, 6, 1) = -0.47800D-02
        HTT(24, 7, 1) =  0.38300D-02
        HTT(24, 8, 1) = -0.16800D-02
        HTT(24, 9, 1) = -0.12500D-02
        HTT(24,10, 1) = -0.13700D-02
        HTT(24,11, 1) = -0.90000D-03
        HTT(24,12, 1) = -0.10100D-02
        HTT(24,13, 1) = -0.31000D-03
        HTT(24,14, 1) = -0.12800D-02
        HTT(24,15, 1) = -0.20000D-03
        HTT(24,16, 1) = -0.11600D-02
        HTT(24,17, 1) =  0.21600D-02
        HTT(24,18, 1) = -0.42700D-02
        HTT(24,19, 1) =  0.31200D-02
        HTT(24,20, 1) = -0.33100D-02
        HTT(24,21, 1) = -0.23000D-02
        HTT(24,22, 1) =  0.35300D-02
        HTT(24,23, 1) =  0.41900D-02
        HTT(24,24, 1) =  0.10010D-01
C Hessian elements for anchor point 1
        HBS( 1, 1, 1) = -0.26190D-01
        HBS( 1, 2, 1) = -0.29060D-01
        HBS( 1, 3, 1) = -0.39330D-01
        HBS( 1, 4, 1) = -0.78300D-02
        HBS( 1, 5, 1) =  0.21800D-02
        HBS( 1, 6, 1) =  0.40640D-01
        HBS( 1, 7, 1) = -0.41000D-03
        HBS( 1, 8, 1) =  0.40740D-01
        HBS( 1, 9, 1) = -0.20000D-04
        HBS( 1,10, 1) = -0.87600D-02
        HBS( 1,11, 1) = -0.35000D-03
        HBS( 1,12, 1) =  0.45000D-02
        HBS( 2, 1, 1) =  0.36740D-01
        HBS( 2, 2, 1) = -0.10920D-01
        HBS( 2, 3, 1) =  0.13650D-01
        HBS( 2, 4, 1) =  0.34200D-02
        HBS( 2, 5, 1) =  0.62000D-03
        HBS( 2, 6, 1) = -0.29580D-01
        HBS( 2, 7, 1) =  0.13300D-02
        HBS( 2, 8, 1) = -0.12060D-01
        HBS( 2, 9, 1) =  0.19000D-03
        HBS( 2,10, 1) =  0.40800D-02
        HBS( 2,11, 1) = -0.39000D-03
        HBS( 2,12, 1) = -0.22000D-03
        HBS( 3, 1, 1) = -0.10550D-01
        HBS( 3, 2, 1) =  0.39980D-01
        HBS( 3, 3, 1) =  0.25680D-01
        HBS( 3, 4, 1) =  0.44100D-02
        HBS( 3, 5, 1) = -0.28000D-02
        HBS( 3, 6, 1) = -0.11060D-01
        HBS( 3, 7, 1) = -0.91000D-03
        HBS( 3, 8, 1) = -0.28680D-01
        HBS( 3, 9, 1) = -0.17000D-03
        HBS( 3,10, 1) =  0.46800D-02
        HBS( 3,11, 1) =  0.74000D-03
        HBS( 3,12, 1) = -0.42800D-02
        HBS( 4, 1, 1) = -0.26180D-01
        HBS( 4, 2, 1) = -0.40500D-02
        HBS( 4, 3, 1) =  0.14760D-01
        HBS( 4, 4, 1) = -0.35070D-01
        HBS( 4, 5, 1) = -0.56800D-02
        HBS( 4, 6, 1) = -0.74200D-02
        HBS( 4, 7, 1) =  0.38500D-02
        HBS( 4, 8, 1) =  0.33760D-01
        HBS( 4, 9, 1) = -0.72000D-03
        HBS( 4,10, 1) =  0.51570D-01
        HBS( 4,11, 1) =  0.20000D-04
        HBS( 4,12, 1) = -0.81000D-03
        HBS( 5, 1, 1) =  0.28180D-01
        HBS( 5, 2, 1) =  0.24300D-02
        HBS( 5, 3, 1) = -0.87400D-02
        HBS( 5, 4, 1) =  0.50900D-02
        HBS( 5, 5, 1) =  0.11940D-01
        HBS( 5, 6, 1) =  0.32000D-02
        HBS( 5, 7, 1) = -0.88000D-03
        HBS( 5, 8, 1) = -0.14770D-01
        HBS( 5, 9, 1) =  0.30000D-04
        HBS( 5,10, 1) = -0.27450D-01
        HBS( 5,11, 1) =  0.00000D+00
        HBS( 5,12, 1) =  0.11400D-02
        HBS( 6, 1, 1) = -0.20100D-02
        HBS( 6, 2, 1) =  0.16200D-02
        HBS( 6, 3, 1) = -0.60300D-02
        HBS( 6, 4, 1) =  0.29980D-01
        HBS( 6, 5, 1) = -0.62600D-02
        HBS( 6, 6, 1) =  0.42200D-02
        HBS( 6, 7, 1) = -0.29700D-02
        HBS( 6, 8, 1) = -0.18990D-01
        HBS( 6, 9, 1) =  0.69000D-03
        HBS( 6,10, 1) = -0.24120D-01
        HBS( 6,11, 1) = -0.10000D-04
        HBS( 6,12, 1) = -0.33000D-03
        HBS( 7, 1, 1) = -0.59200D-02
        HBS( 7, 2, 1) =  0.29370D-01
        HBS( 7, 3, 1) =  0.59700D-02
        HBS( 7, 4, 1) = -0.34310D-01
        HBS( 7, 5, 1) =  0.35700D-02
        HBS( 7, 6, 1) = -0.38110D-01
        HBS( 7, 7, 1) = -0.77800D-02
        HBS( 7, 8, 1) = -0.29600D-02
        HBS( 7, 9, 1) =  0.44000D-02
        HBS( 7,10, 1) =  0.44330D-01
        HBS( 7,11, 1) = -0.52000D-03
        HBS( 7,12, 1) = -0.41000D-03
        HBS( 8, 1, 1) =  0.27500D-02
        HBS( 8, 2, 1) = -0.16290D-01
        HBS( 8, 3, 1) = -0.20500D-02
        HBS( 8, 4, 1) =  0.31190D-01
        HBS( 8, 5, 1) = -0.23800D-02
        HBS( 8, 6, 1) =  0.45900D-02
        HBS( 8, 7, 1) =  0.44200D-02
        HBS( 8, 8, 1) =  0.12900D-02
        HBS( 8, 9, 1) = -0.16900D-02
        HBS( 8,10, 1) = -0.19840D-01
        HBS( 8,11, 1) = -0.28000D-03
        HBS( 8,12, 1) =  0.25000D-03
        HBS( 9, 1, 1) =  0.31600D-02
        HBS( 9, 2, 1) = -0.13080D-01
        HBS( 9, 3, 1) = -0.39200D-02
        HBS( 9, 4, 1) =  0.31100D-02
        HBS( 9, 5, 1) = -0.11900D-02
        HBS( 9, 6, 1) =  0.33520D-01
        HBS( 9, 7, 1) =  0.33700D-02
        HBS( 9, 8, 1) =  0.16700D-02
        HBS( 9, 9, 1) = -0.27000D-02
        HBS( 9,10, 1) = -0.24490D-01
        HBS( 9,11, 1) =  0.80000D-03
        HBS( 9,12, 1) =  0.17000D-03
        HBS(10, 1, 1) =  0.35610D-01
        HBS(10, 2, 1) =  0.35190D-01
        HBS(10, 3, 1) = -0.80000D-03
        HBS(10, 4, 1) = -0.15230D-01
        HBS(10, 5, 1) =  0.40000D-04
        HBS(10, 6, 1) = -0.28120D-01
        HBS(10, 7, 1) =  0.45800D-02
        HBS(10, 8, 1) = -0.27160D-01
        HBS(10, 9, 1) = -0.74500D-02
        HBS(10,10, 1) = -0.14900D-01
        HBS(10,11, 1) =  0.44000D-02
        HBS(10,12, 1) = -0.51000D-03
        HBS(11, 1, 1) = -0.19830D-01
        HBS(11, 2, 1) = -0.15670D-01
        HBS(11, 3, 1) =  0.46000D-03
        HBS(11, 4, 1) =  0.87500D-02
        HBS(11, 5, 1) =  0.12000D-03
        HBS(11, 6, 1) =  0.26580D-01
        HBS(11, 7, 1) = -0.28700D-02
        HBS(11, 8, 1) =  0.69000D-03
        HBS(11, 9, 1) =  0.35900D-02
        HBS(11,10, 1) =  0.63600D-02
        HBS(11,11, 1) = -0.16100D-02
        HBS(11,12, 1) = -0.22000D-03
        HBS(12, 1, 1) = -0.15780D-01
        HBS(12, 2, 1) = -0.19520D-01
        HBS(12, 3, 1) =  0.34000D-03
        HBS(12, 4, 1) =  0.64800D-02
        HBS(12, 5, 1) = -0.15000D-03
        HBS(12, 6, 1) =  0.15400D-02
        HBS(12, 7, 1) = -0.17100D-02
        HBS(12, 8, 1) =  0.26480D-01
        HBS(12, 9, 1) =  0.38600D-02
        HBS(12,10, 1) =  0.85400D-02
        HBS(12,11, 1) = -0.27900D-02
        HBS(12,12, 1) =  0.73000D-03
        HBS(13, 1, 1) =  0.28450D-01
        HBS(13, 2, 1) = -0.45800D-02
        HBS(13, 3, 1) =  0.32200D-02
        HBS(13, 4, 1) =  0.42370D-01
        HBS(13, 5, 1) = -0.53000D-03
        HBS(13, 6, 1) = -0.11000D-02
        HBS(13, 7, 1) = -0.47000D-03
        HBS(13, 8, 1) = -0.36660D-01
        HBS(13, 9, 1) =  0.43000D-02
        HBS(13,10, 1) = -0.34790D-01
        HBS(13,11, 1) = -0.78500D-02
        HBS(13,12, 1) =  0.43000D-02
        HBS(14, 1, 1) = -0.12720D-01
        HBS(14, 2, 1) =  0.26100D-02
        HBS(14, 3, 1) = -0.28000D-02
        HBS(14, 4, 1) = -0.23620D-01
        HBS(14, 5, 1) =  0.19000D-03
        HBS(14, 6, 1) =  0.10100D-02
        HBS(14, 7, 1) =  0.75000D-03
        HBS(14, 8, 1) =  0.32490D-01
        HBS(14, 9, 1) = -0.26900D-02
        HBS(14,10, 1) =  0.34600D-02
        HBS(14,11, 1) =  0.37900D-02
        HBS(14,12, 1) = -0.16300D-02
        HBS(15, 1, 1) = -0.15730D-01
        HBS(15, 2, 1) =  0.19700D-02
        HBS(15, 3, 1) = -0.42000D-03
        HBS(15, 4, 1) = -0.18750D-01
        HBS(15, 5, 1) =  0.33000D-03
        HBS(15, 6, 1) =  0.90000D-04
        HBS(15, 7, 1) = -0.28000D-03
        HBS(15, 8, 1) =  0.41800D-02
        HBS(15, 9, 1) = -0.16100D-02
        HBS(15,10, 1) =  0.31330D-01
        HBS(15,11, 1) =  0.40600D-02
        HBS(15,12, 1) = -0.26700D-02
        HBS(16, 1, 1) = -0.57800D-02
        HBS(16, 2, 1) = -0.26870D-01
        HBS(16, 3, 1) =  0.16180D-01
        HBS(16, 4, 1) =  0.50070D-01
        HBS(16, 5, 1) =  0.42000D-03
        HBS(16, 6, 1) =  0.34120D-01
        HBS(16, 7, 1) =  0.25000D-03
        HBS(16, 8, 1) = -0.77100D-02
        HBS(16, 9, 1) = -0.51000D-03
        HBS(16,10, 1) = -0.37450D-01
        HBS(16,11, 1) =  0.43000D-02
        HBS(16,12, 1) = -0.70600D-02
        HBS(17, 1, 1) =  0.20300D-02
        HBS(17, 2, 1) =  0.26580D-01
        HBS(17, 3, 1) = -0.98100D-02
        HBS(17, 4, 1) = -0.26570D-01
        HBS(17, 5, 1) =  0.45000D-03
        HBS(17, 6, 1) = -0.15140D-01
        HBS(17, 7, 1) = -0.15000D-03
        HBS(17, 8, 1) =  0.35400D-02
        HBS(17, 9, 1) = -0.25000D-03
        HBS(17,10, 1) =  0.60500D-02
        HBS(17,11, 1) = -0.15900D-02
        HBS(17,12, 1) =  0.61300D-02
        HBS(18, 1, 1) =  0.37500D-02
        HBS(18, 2, 1) =  0.29000D-03
        HBS(18, 3, 1) = -0.63700D-02
        HBS(18, 4, 1) = -0.23500D-01
        HBS(18, 5, 1) = -0.87000D-03
        HBS(18, 6, 1) = -0.18980D-01
        HBS(18, 7, 1) = -0.10000D-03
        HBS(18, 8, 1) =  0.41700D-02
        HBS(18, 9, 1) =  0.76000D-03
        HBS(18,10, 1) =  0.31400D-01
        HBS(18,11, 1) = -0.27100D-02
        HBS(18,12, 1) =  0.93000D-03
C Hessian elements for anchor point 2
        HSS( 1, 1, 2) =  0.30793D+00
        HSS( 2, 1, 2) =  0.25300D-02
        HSS( 2, 2, 2) =  0.31227D+00
        HSS( 3, 1, 2) =  0.72500D-01
        HSS( 3, 2, 2) =  0.70950D-01
        HSS( 3, 3, 2) =  0.63382D+00
        HSS( 4, 1, 2) =  0.56300D-02
        HSS( 4, 2, 2) =  0.20240D-01
        HSS( 4, 3, 2) = -0.11340D-01
        HSS( 4, 4, 2) =  0.43881D+00
        HSS( 5, 1, 2) =  0.31200D-02
        HSS( 5, 2, 2) = -0.58000D-03
        HSS( 5, 3, 2) = -0.14100D-02
        HSS( 5, 4, 2) =  0.19700D-02
        HSS( 5, 5, 2) =  0.36903D+00
        HSS( 6, 1, 2) =  0.93200D-02
        HSS( 6, 2, 2) =  0.48920D-01
        HSS( 6, 3, 2) =  0.11310D-01
        HSS( 6, 4, 2) =  0.43540D-01
        HSS( 6, 5, 2) = -0.31000D-03
        HSS( 6, 6, 2) =  0.34629D+00
        HSS( 7, 1, 2) = -0.79000D-03
        HSS( 7, 2, 2) = -0.53000D-03
        HSS( 7, 3, 2) =  0.11700D-02
        HSS( 7, 4, 2) =  0.55200D-02
        HSS( 7, 5, 2) =  0.15000D-03
        HSS( 7, 6, 2) =  0.35500D-02
        HSS( 7, 7, 2) =  0.36483D+00
        HSS( 8, 1, 2) =  0.48320D-01
        HSS( 8, 2, 2) =  0.88800D-02
        HSS( 8, 3, 2) =  0.10960D-01
        HSS( 8, 4, 2) = -0.62700D-02
        HSS( 8, 5, 2) = -0.47000D-03
        HSS( 8, 6, 2) =  0.30770D-01
        HSS( 8, 7, 2) = -0.33000D-03
        HSS( 8, 8, 2) =  0.34600D+00
        HSS( 9, 1, 2) = -0.27000D-03
        HSS( 9, 2, 2) = -0.90000D-04
        HSS( 9, 3, 2) = -0.37000D-03
        HSS( 9, 4, 2) = -0.49000D-03
        HSS( 9, 5, 2) = -0.90000D-04
        HSS( 9, 6, 2) =  0.39700D-02
        HSS( 9, 7, 2) =  0.32000D-03
        HSS( 9, 8, 2) =  0.40100D-02
        HSS( 9, 9, 2) =  0.36524D+00
        HSS(10, 1, 2) =  0.17870D-01
        HSS(10, 2, 2) =  0.45600D-02
        HSS(10, 3, 2) = -0.11880D-01
        HSS(10, 4, 2) =  0.81470D-01
        HSS(10, 5, 2) = -0.59000D-03
        HSS(10, 6, 2) = -0.42200D-02
        HSS(10, 7, 2) = -0.61000D-03
        HSS(10, 8, 2) =  0.41780D-01
        HSS(10, 9, 2) = -0.69000D-03
        HSS(10,10, 2) =  0.44702D+00
        HSS(11, 1, 2) = -0.59000D-03
        HSS(11, 2, 2) = -0.11700D-02
        HSS(11, 3, 2) =  0.11000D-02
        HSS(11, 4, 2) = -0.48000D-03
        HSS(11, 5, 2) =  0.30000D-04
        HSS(11, 6, 2) = -0.44000D-03
        HSS(11, 7, 2) =  0.16000D-03
        HSS(11, 8, 2) =  0.35200D-02
        HSS(11, 9, 2) =  0.40000D-03
        HSS(11,10, 2) =  0.50000D-02
        HSS(11,11, 2) =  0.36451D+00
        HSS(12, 1, 2) =  0.14000D-03
        HSS(12, 2, 2) =  0.21900D-02
        HSS(12, 3, 2) =  0.30000D-04
        HSS(12, 4, 2) = -0.45000D-03
        HSS(12, 5, 2) =  0.10000D-04
        HSS(12, 6, 2) = -0.43000D-03
        HSS(12, 7, 2) =  0.60000D-04
        HSS(12, 8, 2) =  0.80000D-04
        HSS(12, 9, 2) =  0.30000D-04
        HSS(12,10, 2) =  0.27300D-02
        HSS(12,11, 2) =  0.33000D-03
        HSS(12,12, 2) =  0.36751D+00
C Hessian elements for anchor point 2
        HBB( 1, 1, 2) =  0.81080D-01
        HBB( 2, 1, 2) = -0.44810D-01
        HBB( 2, 2, 2) =  0.14731D+00
        HBB( 3, 1, 2) = -0.36270D-01
        HBB( 3, 2, 2) = -0.10250D+00
        HBB( 3, 3, 2) =  0.13877D+00
        HBB( 4, 1, 2) = -0.36710D-01
        HBB( 4, 2, 2) =  0.25500D-01
        HBB( 4, 3, 2) =  0.11210D-01
        HBB( 4, 4, 2) =  0.75960D-01
        HBB( 5, 1, 2) =  0.22750D-01
        HBB( 5, 2, 2) = -0.18060D-01
        HBB( 5, 3, 2) = -0.46900D-02
        HBB( 5, 4, 2) = -0.38310D-01
        HBB( 5, 5, 2) =  0.72880D-01
        HBB( 6, 1, 2) =  0.13960D-01
        HBB( 6, 2, 2) = -0.74300D-02
        HBB( 6, 3, 2) = -0.65200D-02
        HBB( 6, 4, 2) = -0.37640D-01
        HBB( 6, 5, 2) = -0.34570D-01
        HBB( 6, 6, 2) =  0.72210D-01
        HBB( 7, 1, 2) = -0.88100D-02
        HBB( 7, 2, 2) = -0.49400D-02
        HBB( 7, 3, 2) =  0.13750D-01
        HBB( 7, 4, 2) = -0.28850D-01
        HBB( 7, 5, 2) =  0.11260D-01
        HBB( 7, 6, 2) =  0.17590D-01
        HBB( 7, 7, 2) =  0.74530D-01
        HBB( 8, 1, 2) = -0.15100D-02
        HBB( 8, 2, 2) =  0.76700D-02
        HBB( 8, 3, 2) = -0.61600D-02
        HBB( 8, 4, 2) =  0.17480D-01
        HBB( 8, 5, 2) = -0.61800D-02
        HBB( 8, 6, 2) = -0.11300D-01
        HBB( 8, 7, 2) = -0.36990D-01
        HBB( 8, 8, 2) =  0.78710D-01
        HBB( 9, 1, 2) =  0.10320D-01
        HBB( 9, 2, 2) = -0.27300D-02
        HBB( 9, 3, 2) = -0.75900D-02
        HBB( 9, 4, 2) =  0.11370D-01
        HBB( 9, 5, 2) = -0.50800D-02
        HBB( 9, 6, 2) = -0.62900D-02
        HBB( 9, 7, 2) = -0.37540D-01
        HBB( 9, 8, 2) = -0.41720D-01
        HBB( 9, 9, 2) =  0.79260D-01
        HBB(10, 1, 2) =  0.10400D-01
        HBB(10, 2, 2) = -0.60200D-02
        HBB(10, 3, 2) = -0.43800D-02
        HBB(10, 4, 2) = -0.90700D-02
        HBB(10, 5, 2) =  0.10230D-01
        HBB(10, 6, 2) = -0.11600D-02
        HBB(10, 7, 2) = -0.34090D-01
        HBB(10, 8, 2) =  0.13820D-01
        HBB(10, 9, 2) =  0.20260D-01
        HBB(10,10, 2) =  0.76010D-01
        HBB(11, 1, 2) = -0.53400D-02
        HBB(11, 2, 2) =  0.37500D-02
        HBB(11, 3, 2) =  0.15900D-02
        HBB(11, 4, 2) = -0.11000D-02
        HBB(11, 5, 2) = -0.35000D-02
        HBB(11, 6, 2) =  0.46000D-02
        HBB(11, 7, 2) =  0.20040D-01
        HBB(11, 8, 2) = -0.74900D-02
        HBB(11, 9, 2) = -0.12560D-01
        HBB(11,10, 2) = -0.37950D-01
        HBB(11,11, 2) =  0.79720D-01
        HBB(12, 1, 2) = -0.50600D-02
        HBB(12, 2, 2) =  0.22700D-02
        HBB(12, 3, 2) =  0.27900D-02
        HBB(12, 4, 2) =  0.10170D-01
        HBB(12, 5, 2) = -0.67300D-02
        HBB(12, 6, 2) = -0.34400D-02
        HBB(12, 7, 2) =  0.14040D-01
        HBB(12, 8, 2) = -0.63400D-02
        HBB(12, 9, 2) = -0.77100D-02
        HBB(12,10, 2) = -0.38060D-01
        HBB(12,11, 2) = -0.41770D-01
        HBB(12,12, 2) =  0.79830D-01
        HBB(13, 1, 2) = -0.81300D-02
        HBB(13, 2, 2) =  0.15080D-01
        HBB(13, 3, 2) = -0.69500D-02
        HBB(13, 4, 2) =  0.40600D-02
        HBB(13, 5, 2) = -0.19200D-02
        HBB(13, 6, 2) = -0.21400D-02
        HBB(13, 7, 2) = -0.67600D-02
        HBB(13, 8, 2) =  0.90200D-02
        HBB(13, 9, 2) = -0.22600D-02
        HBB(13,10, 2) = -0.33920D-01
        HBB(13,11, 2) =  0.13940D-01
        HBB(13,12, 2) =  0.19980D-01
        HBB(13,13, 2) =  0.73910D-01
        HBB(14, 1, 2) =  0.96700D-02
        HBB(14, 2, 2) = -0.84100D-02
        HBB(14, 3, 2) = -0.12600D-02
        HBB(14, 4, 2) = -0.20700D-02
        HBB(14, 5, 2) =  0.84000D-03
        HBB(14, 6, 2) =  0.12300D-02
        HBB(14, 7, 2) = -0.21800D-02
        HBB(14, 8, 2) = -0.30500D-02
        HBB(14, 9, 2) =  0.52300D-02
        HBB(14,10, 2) =  0.20350D-01
        HBB(14,11, 2) = -0.77300D-02
        HBB(14,12, 2) = -0.12620D-01
        HBB(14,13, 2) = -0.37510D-01
        HBB(14,14, 2) =  0.78920D-01
        HBB(15, 1, 2) = -0.15300D-02
        HBB(15, 2, 2) = -0.66700D-02
        HBB(15, 3, 2) =  0.82100D-02
        HBB(15, 4, 2) = -0.20000D-02
        HBB(15, 5, 2) =  0.10900D-02
        HBB(15, 6, 2) =  0.91000D-03
        HBB(15, 7, 2) =  0.89400D-02
        HBB(15, 8, 2) = -0.59700D-02
        HBB(15, 9, 2) = -0.29700D-02
        HBB(15,10, 2) =  0.13570D-01
        HBB(15,11, 2) = -0.62100D-02
        HBB(15,12, 2) = -0.73600D-02
        HBB(15,13, 2) = -0.36400D-01
        HBB(15,14, 2) = -0.41400D-01
        HBB(15,15, 2) =  0.77800D-01
        HBB(16, 1, 2) = -0.37830D-01
        HBB(16, 2, 2) =  0.15190D-01
        HBB(16, 3, 2) =  0.22640D-01
        HBB(16, 4, 2) = -0.53900D-02
        HBB(16, 5, 2) = -0.40000D-02
        HBB(16, 6, 2) =  0.93900D-02
        HBB(16, 7, 2) =  0.39700D-02
        HBB(16, 8, 2) = -0.18200D-02
        HBB(16, 9, 2) = -0.21500D-02
        HBB(16,10, 2) = -0.93400D-02
        HBB(16,11, 2) =  0.10420D-01
        HBB(16,12, 2) = -0.10800D-02
        HBB(16,13, 2) = -0.29170D-01
        HBB(16,14, 2) =  0.11750D-01
        HBB(16,15, 2) =  0.17420D-01
        HBB(16,16, 2) =  0.77760D-01
        HBB(17, 1, 2) =  0.23450D-01
        HBB(17, 2, 2) = -0.61200D-02
        HBB(17, 3, 2) = -0.17330D-01
        HBB(17, 4, 2) = -0.38800D-02
        HBB(17, 5, 2) =  0.65900D-02
        HBB(17, 6, 2) = -0.27000D-02
        HBB(17, 7, 2) = -0.18400D-02
        HBB(17, 8, 2) =  0.10400D-02
        HBB(17, 9, 2) =  0.79000D-03
        HBB(17,10, 2) =  0.10460D-01
        HBB(17,11, 2) = -0.67900D-02
        HBB(17,12, 2) = -0.36700D-02
        HBB(17,13, 2) =  0.11210D-01
        HBB(17,14, 2) = -0.50900D-02
        HBB(17,15, 2) = -0.61300D-02
        HBB(17,16, 2) = -0.39400D-01
        HBB(17,17, 2) =  0.75830D-01
        HBB(18, 1, 2) =  0.14380D-01
        HBB(18, 2, 2) = -0.90700D-02
        HBB(18, 3, 2) = -0.53200D-02
        HBB(18, 4, 2) =  0.92800D-02
        HBB(18, 5, 2) = -0.25900D-02
        HBB(18, 6, 2) = -0.66900D-02
        HBB(18, 7, 2) = -0.21300D-02
        HBB(18, 8, 2) =  0.77000D-03
        HBB(18, 9, 2) =  0.13500D-02
        HBB(18,10, 2) = -0.11200D-02
        HBB(18,11, 2) = -0.36200D-02
        HBB(18,12, 2) =  0.47500D-02
        HBB(18,13, 2) =  0.17950D-01
        HBB(18,14, 2) = -0.66600D-02
        HBB(18,15, 2) = -0.11290D-01
        HBB(18,16, 2) = -0.38360D-01
        HBB(18,17, 2) = -0.36420D-01
        HBB(18,18, 2) =  0.74790D-01
C Hessian elements for anchor point 2
        HTT( 1, 1, 2) =  0.87500D-02
        HTT( 2, 1, 2) =  0.39300D-02
        HTT( 2, 2, 2) =  0.86800D-02
        HTT( 3, 1, 2) = -0.25000D-03
        HTT( 3, 2, 2) = -0.40300D-02
        HTT( 3, 3, 2) =  0.69500D-02
        HTT( 4, 1, 2) = -0.50700D-02
        HTT( 4, 2, 2) =  0.71000D-03
        HTT( 4, 3, 2) =  0.31700D-02
        HTT( 4, 4, 2) =  0.89500D-02
        HTT( 5, 1, 2) = -0.63800D-02
        HTT( 5, 2, 2) = -0.48700D-02
        HTT( 5, 3, 2) =  0.26100D-02
        HTT( 5, 4, 2) =  0.41100D-02
        HTT( 5, 5, 2) =  0.87900D-02
        HTT( 6, 1, 2) = -0.41700D-02
        HTT( 6, 2, 2) = -0.31900D-02
        HTT( 6, 3, 2) =  0.27300D-02
        HTT( 6, 4, 2) =  0.37100D-02
        HTT( 6, 5, 2) =  0.34900D-02
        HTT( 6, 6, 2) =  0.76600D-02
        HTT( 7, 1, 2) =  0.23700D-02
        HTT( 7, 2, 2) =  0.28700D-02
        HTT( 7, 3, 2) = -0.43900D-02
        HTT( 7, 4, 2) = -0.39000D-02
        HTT( 7, 5, 2) =  0.60000D-04
        HTT( 7, 6, 2) = -0.32200D-02
        HTT( 7, 7, 2) =  0.66300D-02
        HTT( 8, 1, 2) =  0.45800D-02
        HTT( 8, 2, 2) =  0.45500D-02
        HTT( 8, 3, 2) = -0.42700D-02
        HTT( 8, 4, 2) = -0.43000D-02
        HTT( 8, 5, 2) = -0.52500D-02
        HTT( 8, 6, 2) =  0.95000D-03
        HTT( 8, 7, 2) =  0.33500D-02
        HTT( 8, 8, 2) =  0.95500D-02
        HTT( 9, 1, 2) = -0.59200D-02
        HTT( 9, 2, 2) = -0.90000D-03
        HTT( 9, 3, 2) = -0.18700D-02
        HTT( 9, 4, 2) =  0.31600D-02
        HTT( 9, 5, 2) =  0.12000D-02
        HTT( 9, 6, 2) =  0.22300D-02
        HTT( 9, 7, 2) = -0.27400D-02
        HTT( 9, 8, 2) = -0.17100D-02
        HTT( 9, 9, 2) =  0.79900D-02
        HTT(10, 1, 2) = -0.37700D-02
        HTT(10, 2, 2) =  0.64000D-03
        HTT(10, 3, 2) = -0.10100D-02
        HTT(10, 4, 2) =  0.33900D-02
        HTT(10, 5, 2) =  0.20800D-02
        HTT(10, 6, 2) =  0.14400D-02
        HTT(10, 7, 2) = -0.59000D-03
        HTT(10, 8, 2) = -0.12400D-02
        HTT(10, 9, 2) =  0.28400D-02
        HTT(10,10, 2) =  0.84900D-02
        HTT(11, 1, 2) = -0.84000D-03
        HTT(11, 2, 2) = -0.59000D-02
        HTT(11, 3, 2) =  0.21200D-02
        HTT(11, 4, 2) = -0.29400D-02
        HTT(11, 5, 2) = -0.38000D-03
        HTT(11, 6, 2) =  0.11900D-02
        HTT(11, 7, 2) = -0.32600D-02
        HTT(11, 8, 2) = -0.16800D-02
        HTT(11, 9, 2) =  0.26900D-02
        HTT(11,10, 2) = -0.18100D-02
        HTT(11,11, 2) =  0.80400D-02
        HTT(12, 1, 2) =  0.13200D-02
        HTT(12, 2, 2) = -0.43700D-02
        HTT(12, 3, 2) =  0.29800D-02
        HTT(12, 4, 2) = -0.27100D-02
        HTT(12, 5, 2) =  0.50000D-03
        HTT(12, 6, 2) =  0.41000D-03
        HTT(12, 7, 2) = -0.11200D-02
        HTT(12, 8, 2) = -0.12100D-02
        HTT(12, 9, 2) = -0.24600D-02
        HTT(12,10, 2) =  0.38400D-02
        HTT(12,11, 2) =  0.35400D-02
        HTT(12,12, 2) =  0.98300D-02
        HTT(13, 1, 2) =  0.53000D-03
        HTT(13, 2, 2) = -0.13800D-02
        HTT(13, 3, 2) =  0.17500D-02
        HTT(13, 4, 2) = -0.16000D-03
        HTT(13, 5, 2) =  0.18900D-02
        HTT(13, 6, 2) =  0.51000D-03
        HTT(13, 7, 2) =  0.70000D-03
        HTT(13, 8, 2) = -0.68000D-03
        HTT(13, 9, 2) = -0.54000D-02
        HTT(13,10, 2) = -0.17000D-03
        HTT(13,11, 2) = -0.33800D-02
        HTT(13,12, 2) =  0.18400D-02
        HTT(13,13, 2) =  0.80400D-02
        HTT(14, 1, 2) =  0.19700D-02
        HTT(14, 2, 2) =  0.54000D-03
        HTT(14, 3, 2) =  0.79000D-03
        HTT(14, 4, 2) = -0.65000D-03
        HTT(14, 5, 2) =  0.36000D-03
        HTT(14, 6, 2) = -0.13400D-02
        HTT(14, 7, 2) =  0.15100D-02
        HTT(14, 8, 2) = -0.19000D-03
        HTT(14, 9, 2) = -0.30400D-02
        HTT(14,10, 2) =  0.14100D-02
        HTT(14,11, 2) = -0.15300D-02
        HTT(14,12, 2) =  0.29300D-02
        HTT(14,13, 2) =  0.18300D-02
        HTT(14,14, 2) =  0.71800D-02
        HTT(15, 1, 2) = -0.16100D-02
        HTT(15, 2, 2) = -0.29100D-02
        HTT(15, 3, 2) =  0.90000D-03
        HTT(15, 4, 2) = -0.39000D-03
        HTT(15, 5, 2) =  0.10100D-02
        HTT(15, 6, 2) =  0.13000D-02
        HTT(15, 7, 2) = -0.14300D-02
        HTT(15, 8, 2) = -0.11500D-02
        HTT(15, 9, 2) = -0.28000D-03
        HTT(15,10, 2) = -0.57800D-02
        HTT(15,11, 2) =  0.10900D-02
        HTT(15,12, 2) = -0.44100D-02
        HTT(15,13, 2) =  0.28500D-02
        HTT(15,14, 2) = -0.26000D-02
        HTT(15,15, 2) =  0.83100D-02
        HTT(16, 1, 2) = -0.16000D-03
        HTT(16, 2, 2) = -0.99000D-03
        HTT(16, 3, 2) = -0.50000D-04
        HTT(16, 4, 2) = -0.88000D-03
        HTT(16, 5, 2) = -0.52000D-03
        HTT(16, 6, 2) = -0.55000D-03
        HTT(16, 7, 2) = -0.62000D-03
        HTT(16, 8, 2) = -0.66000D-03
        HTT(16, 9, 2) =  0.20700D-02
        HTT(16,10, 2) = -0.42000D-02
        HTT(16,11, 2) =  0.29500D-02
        HTT(16,12, 2) = -0.33300D-02
        HTT(16,13, 2) = -0.33700D-02
        HTT(16,14, 2) =  0.27500D-02
        HTT(16,15, 2) =  0.28600D-02
        HTT(16,16, 2) =  0.89800D-02
        HTT(17, 1, 2) =  0.19500D-02
        HTT(17, 2, 2) =  0.45000D-03
        HTT(17, 3, 2) =  0.64000D-03
        HTT(17, 4, 2) = -0.85000D-03
        HTT(17, 5, 2) =  0.51000D-03
        HTT(17, 6, 2) = -0.12300D-02
        HTT(17, 7, 2) =  0.17800D-02
        HTT(17, 8, 2) =  0.40000D-04
        HTT(17, 9, 2) =  0.55000D-03
        HTT(17,10, 2) = -0.15600D-02
        HTT(17,11, 2) =  0.21200D-02
        HTT(17,12, 2) =  0.10000D-04
        HTT(17,13, 2) = -0.55500D-02
        HTT(17,14, 2) =  0.56000D-03
        HTT(17,15, 2) = -0.34500D-02
        HTT(17,16, 2) =  0.26500D-02
        HTT(17,17, 2) =  0.80100D-02
        HTT(18, 1, 2) =  0.68000D-03
        HTT(18, 2, 2) =  0.11500D-02
        HTT(18, 3, 2) = -0.14100D-02
        HTT(18, 4, 2) = -0.93000D-03
        HTT(18, 5, 2) = -0.14500D-02
        HTT(18, 6, 2) = -0.25900D-02
        HTT(18, 7, 2) =  0.58000D-03
        HTT(18, 8, 2) = -0.57000D-03
        HTT(18, 9, 2) =  0.17600D-02
        HTT(18,10, 2) =  0.54000D-03
        HTT(18,11, 2) =  0.12600D-02
        HTT(18,12, 2) =  0.40000D-04
        HTT(18,13, 2) = -0.35300D-02
        HTT(18,14, 2) =  0.18700D-02
        HTT(18,15, 2) = -0.23200D-02
        HTT(18,16, 2) =  0.30800D-02
        HTT(18,17, 2) =  0.27300D-02
        HTT(18,18, 2) =  0.81100D-02
        HTT(19, 1, 2) =  0.50000D-03
        HTT(19, 2, 2) = -0.14600D-02
        HTT(19, 3, 2) =  0.16000D-02
        HTT(19, 4, 2) = -0.37000D-03
        HTT(19, 5, 2) =  0.20400D-02
        HTT(19, 6, 2) =  0.62000D-03
        HTT(19, 7, 2) =  0.97000D-03
        HTT(19, 8, 2) = -0.45000D-03
        HTT(19, 9, 2) = -0.18100D-02
        HTT(19,10, 2) = -0.31500D-02
        HTT(19,11, 2) =  0.26000D-03
        HTT(19,12, 2) = -0.10800D-02
        HTT(19,13, 2) =  0.67000D-03
        HTT(19,14, 2) = -0.48000D-02
        HTT(19,15, 2) =  0.20000D-02
        HTT(19,16, 2) = -0.34700D-02
        HTT(19,17, 2) =  0.19000D-02
        HTT(19,18, 2) = -0.26700D-02
        HTT(19,19, 2) =  0.73700D-02
        HTT(20, 1, 2) = -0.77000D-03
        HTT(20, 2, 2) = -0.76000D-03
        HTT(20, 3, 2) = -0.45000D-03
        HTT(20, 4, 2) = -0.44000D-03
        HTT(20, 5, 2) =  0.80000D-04
        HTT(20, 6, 2) = -0.74000D-03
        HTT(20, 7, 2) = -0.23000D-03
        HTT(20, 8, 2) = -0.10600D-02
        HTT(20, 9, 2) = -0.59000D-03
        HTT(20,10, 2) = -0.10400D-02
        HTT(20,11, 2) = -0.60000D-03
        HTT(20,12, 2) = -0.10500D-02
        HTT(20,13, 2) =  0.26900D-02
        HTT(20,14, 2) = -0.34800D-02
        HTT(20,15, 2) =  0.31300D-02
        HTT(20,16, 2) = -0.30400D-02
        HTT(20,17, 2) = -0.33700D-02
        HTT(20,18, 2) =  0.27100D-02
        HTT(20,19, 2) =  0.28000D-02
        HTT(20,20, 2) =  0.88800D-02
        HTT(21, 1, 2) =  0.11400D-02
        HTT(21, 2, 2) =  0.27200D-02
        HTT(21, 3, 2) = -0.27900D-02
        HTT(21, 4, 2) = -0.12100D-02
        HTT(21, 5, 2) = -0.59100D-02
        HTT(21, 6, 2) = -0.84000D-03
        HTT(21, 7, 2) = -0.20900D-02
        HTT(21, 8, 2) =  0.29800D-02
        HTT(21, 9, 2) =  0.14100D-02
        HTT(21,10, 2) =  0.53000D-03
        HTT(21,11, 2) = -0.26000D-03
        HTT(21,12, 2) = -0.11400D-02
        HTT(21,13, 2) =  0.60000D-03
        HTT(21,14, 2) = -0.15900D-02
        HTT(21,15, 2) =  0.14800D-02
        HTT(21,16, 2) = -0.71000D-03
        HTT(21,17, 2) = -0.54000D-02
        HTT(21,18, 2) = -0.23000D-03
        HTT(21,19, 2) = -0.32100D-02
        HTT(21,20, 2) =  0.19600D-02
        HTT(21,21, 2) =  0.79800D-02
        HTT(22, 1, 2) = -0.11900D-02
        HTT(22, 2, 2) =  0.94000D-03
        HTT(22, 3, 2) = -0.29200D-02
        HTT(22, 4, 2) = -0.79000D-03
        HTT(22, 5, 2) = -0.30000D-03
        HTT(22, 6, 2) = -0.52400D-02
        HTT(22, 7, 2) =  0.13800D-02
        HTT(22, 8, 2) = -0.35600D-02
        HTT(22, 9, 2) =  0.33000D-03
        HTT(22,10, 2) =  0.12100D-02
        HTT(22,11, 2) = -0.19200D-02
        HTT(22,12, 2) = -0.10400D-02
        HTT(22,13, 2) =  0.20500D-02
        HTT(22,14, 2) =  0.20000D-03
        HTT(22,15, 2) =  0.11700D-02
        HTT(22,16, 2) = -0.67000D-03
        HTT(22,17, 2) = -0.35500D-02
        HTT(22,18, 2) =  0.98000D-03
        HTT(22,19, 2) = -0.17100D-02
        HTT(22,20, 2) =  0.28200D-02
        HTT(22,21, 2) =  0.26300D-02
        HTT(22,22, 2) =  0.78500D-02
        HTT(23, 1, 2) =  0.24200D-02
        HTT(23, 2, 2) =  0.20100D-02
        HTT(23, 3, 2) = -0.73000D-03
        HTT(23, 4, 2) = -0.11400D-02
        HTT(23, 5, 2) = -0.39300D-02
        HTT(23, 6, 2) =  0.53000D-03
        HTT(23, 7, 2) = -0.87000D-03
        HTT(23, 8, 2) =  0.35900D-02
        HTT(23, 9, 2) =  0.18000D-03
        HTT(23,10, 2) = -0.16000D-02
        HTT(23,11, 2) =  0.61000D-03
        HTT(23,12, 2) = -0.11700D-02
        HTT(23,13, 2) = -0.14400D-02
        HTT(23,14, 2) = -0.29100D-02
        HTT(23,15, 2) =  0.33000D-03
        HTT(23,16, 2) = -0.11400D-02
        HTT(23,17, 2) = -0.80000D-04
        HTT(23,18, 2) = -0.56500D-02
        HTT(23,19, 2) =  0.14000D-02
        HTT(23,20, 2) = -0.41800D-02
        HTT(23,21, 2) =  0.27800D-02
        HTT(23,22, 2) = -0.19400D-02
        HTT(23,23, 2) =  0.84000D-02
        HTT(24, 1, 2) =  0.90000D-04
        HTT(24, 2, 2) =  0.23000D-03
        HTT(24, 3, 2) = -0.86000D-03
        HTT(24, 4, 2) = -0.71000D-03
        HTT(24, 5, 2) =  0.16700D-02
        HTT(24, 6, 2) = -0.38700D-02
        HTT(24, 7, 2) =  0.25900D-02
        HTT(24, 8, 2) = -0.29500D-02
        HTT(24, 9, 2) = -0.90000D-03
        HTT(24,10, 2) = -0.92000D-03
        HTT(24,11, 2) = -0.10500D-02
        HTT(24,12, 2) = -0.10700D-02
        HTT(24,13, 2) =  0.10000D-04
        HTT(24,14, 2) = -0.11200D-02
        HTT(24,15, 2) =  0.30000D-04
        HTT(24,16, 2) = -0.11000D-02
        HTT(24,17, 2) =  0.17600D-02
        HTT(24,18, 2) = -0.44400D-02
        HTT(24,19, 2) =  0.29000D-02
        HTT(24,20, 2) = -0.33100D-02
        HTT(24,21, 2) = -0.25700D-02
        HTT(24,22, 2) =  0.32800D-02
        HTT(24,23, 2) =  0.36800D-02
        HTT(24,24, 2) =  0.95300D-02
C Hessian elements for anchor point 2
        HBS( 1, 1, 2) = -0.25450D-01
        HBS( 1, 2, 2) = -0.24360D-01
        HBS( 1, 3, 2) = -0.45340D-01
        HBS( 1, 4, 2) = -0.86200D-02
        HBS( 1, 5, 2) =  0.45000D-02
        HBS( 1, 6, 2) =  0.42960D-01
        HBS( 1, 7, 2) = -0.49000D-03
        HBS( 1, 8, 2) =  0.43220D-01
        HBS( 1, 9, 2) = -0.40000D-04
        HBS( 1,10, 2) = -0.84400D-02
        HBS( 1,11, 2) = -0.54000D-03
        HBS( 1,12, 2) =  0.48600D-02
        HBS( 2, 1, 2) =  0.39900D-01
        HBS( 2, 2, 2) = -0.84900D-02
        HBS( 2, 3, 2) =  0.28370D-01
        HBS( 2, 4, 2) =  0.69300D-02
        HBS( 2, 5, 2) = -0.38100D-02
        HBS( 2, 6, 2) = -0.29340D-01
        HBS( 2, 7, 2) =  0.12900D-02
        HBS( 2, 8, 2) = -0.15540D-01
        HBS( 2, 9, 2) =  0.18000D-03
        HBS( 2,10, 2) =  0.28900D-02
        HBS( 2,11, 2) = -0.42000D-03
        HBS( 2,12, 2) = -0.40000D-03
        HBS( 3, 1, 2) = -0.14450D-01
        HBS( 3, 2, 2) =  0.32850D-01
        HBS( 3, 3, 2) =  0.16980D-01
        HBS( 3, 4, 2) =  0.16900D-02
        HBS( 3, 5, 2) = -0.70000D-03
        HBS( 3, 6, 2) = -0.13620D-01
        HBS( 3, 7, 2) = -0.79000D-03
        HBS( 3, 8, 2) = -0.27680D-01
        HBS( 3, 9, 2) = -0.14000D-03
        HBS( 3,10, 2) =  0.55400D-02
        HBS( 3,11, 2) =  0.95000D-03
        HBS( 3,12, 2) = -0.44600D-02
        HBS( 4, 1, 2) = -0.22690D-01
        HBS( 4, 2, 2) = -0.61800D-02
        HBS( 4, 3, 2) =  0.14530D-01
        HBS( 4, 4, 2) = -0.35450D-01
        HBS( 4, 5, 2) = -0.71600D-02
        HBS( 4, 6, 2) = -0.92200D-02
        HBS( 4, 7, 2) =  0.43300D-02
        HBS( 4, 8, 2) =  0.36360D-01
        HBS( 4, 9, 2) = -0.65000D-03
        HBS( 4,10, 2) =  0.49680D-01
        HBS( 4,11, 2) =  0.19000D-03
        HBS( 4,12, 2) = -0.97000D-03
        HBS( 5, 1, 2) =  0.21990D-01
        HBS( 5, 2, 2) =  0.37900D-02
        HBS( 5, 3, 2) = -0.71900D-02
        HBS( 5, 4, 2) =  0.47500D-02
        HBS( 5, 5, 2) =  0.95400D-02
        HBS( 5, 6, 2) =  0.43800D-02
        HBS( 5, 7, 2) = -0.13300D-02
        HBS( 5, 8, 2) = -0.15760D-01
        HBS( 5, 9, 2) = -0.15000D-03
        HBS( 5,10, 2) = -0.26690D-01
        HBS( 5,11, 2) = -0.13000D-03
        HBS( 5,12, 2) =  0.10900D-02
        HBS( 6, 1, 2) =  0.70000D-03
        HBS( 6, 2, 2) =  0.24000D-02
        HBS( 6, 3, 2) = -0.73400D-02
        HBS( 6, 4, 2) =  0.30700D-01
        HBS( 6, 5, 2) = -0.23800D-02
        HBS( 6, 6, 2) =  0.48500D-02
        HBS( 6, 7, 2) = -0.29900D-02
        HBS( 6, 8, 2) = -0.20600D-01
        HBS( 6, 9, 2) =  0.80000D-03
        HBS( 6,10, 2) = -0.22990D-01
        HBS( 6,11, 2) = -0.60000D-04
        HBS( 6,12, 2) = -0.12000D-03
        HBS( 7, 1, 2) = -0.73400D-02
        HBS( 7, 2, 2) =  0.28830D-01
        HBS( 7, 3, 2) =  0.66600D-02
        HBS( 7, 4, 2) = -0.33160D-01
        HBS( 7, 5, 2) =  0.40800D-02
        HBS( 7, 6, 2) = -0.36980D-01
        HBS( 7, 7, 2) = -0.81400D-02
        HBS( 7, 8, 2) = -0.48000D-02
        HBS( 7, 9, 2) =  0.47000D-02
        HBS( 7,10, 2) =  0.43210D-01
        HBS( 7,11, 2) = -0.64000D-03
        HBS( 7,12, 2) = -0.30000D-03
        HBS( 8, 1, 2) =  0.34500D-02
        HBS( 8, 2, 2) = -0.15860D-01
        HBS( 8, 3, 2) = -0.19800D-02
        HBS( 8, 4, 2) =  0.30660D-01
        HBS( 8, 5, 2) = -0.26700D-02
        HBS( 8, 6, 2) =  0.43400D-02
        HBS( 8, 7, 2) =  0.39600D-02
        HBS( 8, 8, 2) =  0.18400D-02
        HBS( 8, 9, 2) = -0.18900D-02
        HBS( 8,10, 2) = -0.19260D-01
        HBS( 8,11, 2) = -0.26000D-03
        HBS( 8,12, 2) =  0.16000D-03
        HBS( 9, 1, 2) =  0.38900D-02
        HBS( 9, 2, 2) = -0.12980D-01
        HBS( 9, 3, 2) = -0.46900D-02
        HBS( 9, 4, 2) =  0.25000D-02
        HBS( 9, 5, 2) = -0.14100D-02
        HBS( 9, 6, 2) =  0.32640D-01
        HBS( 9, 7, 2) =  0.41800D-02
        HBS( 9, 8, 2) =  0.29500D-02
        HBS( 9, 9, 2) = -0.28100D-02
        HBS( 9,10, 2) = -0.23950D-01
        HBS( 9,11, 2) =  0.90000D-03
        HBS( 9,12, 2) =  0.14000D-03
        HBS(10, 1, 2) =  0.32210D-01
        HBS(10, 2, 2) =  0.32440D-01
        HBS(10, 3, 2) =  0.24000D-02
        HBS(10, 4, 2) = -0.14180D-01
        HBS(10, 5, 2) = -0.19000D-03
        HBS(10, 6, 2) = -0.29000D-01
        HBS(10, 7, 2) =  0.46200D-02
        HBS(10, 8, 2) = -0.28530D-01
        HBS(10, 9, 2) = -0.80800D-02
        HBS(10,10, 2) = -0.14370D-01
        HBS(10,11, 2) =  0.45400D-02
        HBS(10,12, 2) = -0.45000D-03
        HBS(11, 1, 2) = -0.18020D-01
        HBS(11, 2, 2) = -0.14150D-01
        HBS(11, 3, 2) = -0.11000D-02
        HBS(11, 4, 2) =  0.83100D-02
        HBS(11, 5, 2) =  0.45000D-03
        HBS(11, 6, 2) =  0.27070D-01
        HBS(11, 7, 2) = -0.28400D-02
        HBS(11, 8, 2) =  0.14300D-02
        HBS(11, 9, 2) =  0.40300D-02
        HBS(11,10, 2) =  0.60300D-02
        HBS(11,11, 2) = -0.17200D-02
        HBS(11,12, 2) = -0.30000D-03
        HBS(12, 1, 2) = -0.14190D-01
        HBS(12, 2, 2) = -0.18290D-01
        HBS(12, 3, 2) = -0.13000D-02
        HBS(12, 4, 2) =  0.58700D-02
        HBS(12, 5, 2) = -0.26000D-03
        HBS(12, 6, 2) =  0.19300D-02
        HBS(12, 7, 2) = -0.17800D-02
        HBS(12, 8, 2) =  0.27100D-01
        HBS(12, 9, 2) =  0.40500D-02
        HBS(12,10, 2) =  0.83400D-02
        HBS(12,11, 2) = -0.28100D-02
        HBS(12,12, 2) =  0.75000D-03
        HBS(13, 1, 2) =  0.28370D-01
        HBS(13, 2, 2) = -0.63900D-02
        HBS(13, 3, 2) =  0.52100D-02
        HBS(13, 4, 2) =  0.42290D-01
        HBS(13, 5, 2) = -0.27000D-03
        HBS(13, 6, 2) = -0.35800D-02
        HBS(13, 7, 2) = -0.64000D-03
        HBS(13, 8, 2) = -0.36310D-01
        HBS(13, 9, 2) =  0.46600D-02
        HBS(13,10, 2) = -0.33260D-01
        HBS(13,11, 2) = -0.81200D-02
        HBS(13,12, 2) =  0.45000D-02
        HBS(14, 1, 2) = -0.12750D-01
        HBS(14, 2, 2) =  0.38700D-02
        HBS(14, 3, 2) = -0.40000D-02
        HBS(14, 4, 2) = -0.23630D-01
        HBS(14, 5, 2) =  0.12000D-03
        HBS(14, 6, 2) =  0.23100D-02
        HBS(14, 7, 2) =  0.89000D-03
        HBS(14, 8, 2) =  0.31980D-01
        HBS(14, 9, 2) = -0.28200D-02
        HBS(14,10, 2) =  0.30000D-02
        HBS(14,11, 2) =  0.42200D-02
        HBS(14,12, 2) = -0.17600D-02
        HBS(15, 1, 2) = -0.15620D-01
        HBS(15, 2, 2) =  0.25200D-02
        HBS(15, 3, 2) = -0.12100D-02
        HBS(15, 4, 2) = -0.18660D-01
        HBS(15, 5, 2) =  0.14000D-03
        HBS(15, 6, 2) =  0.12700D-02
        HBS(15, 7, 2) = -0.25000D-03
        HBS(15, 8, 2) =  0.43300D-02
        HBS(15, 9, 2) = -0.18400D-02
        HBS(15,10, 2) =  0.30260D-01
        HBS(15,11, 2) =  0.39000D-02
        HBS(15,12, 2) = -0.27400D-02
        HBS(16, 1, 2) = -0.51100D-02
        HBS(16, 2, 2) = -0.24340D-01
        HBS(16, 3, 2) =  0.16540D-01
        HBS(16, 4, 2) =  0.49110D-01
        HBS(16, 5, 2) = -0.98000D-03
        HBS(16, 6, 2) =  0.35820D-01
        HBS(16, 7, 2) =  0.33000D-03
        HBS(16, 8, 2) = -0.99500D-02
        HBS(16, 9, 2) = -0.59000D-03
        HBS(16,10, 2) = -0.36820D-01
        HBS(16,11, 2) =  0.45700D-02
        HBS(16,12, 2) = -0.76400D-02
        HBS(17, 1, 2) =  0.19500D-02
        HBS(17, 2, 2) =  0.23750D-01
        HBS(17, 3, 2) = -0.99800D-02
        HBS(17, 4, 2) = -0.26130D-01
        HBS(17, 5, 2) =  0.10100D-02
        HBS(17, 6, 2) = -0.15770D-01
        HBS(17, 7, 2) = -0.17000D-03
        HBS(17, 8, 2) =  0.47200D-02
        HBS(17, 9, 2) = -0.26000D-03
        HBS(17,10, 2) =  0.55900D-02
        HBS(17,11, 2) = -0.16700D-02
        HBS(17,12, 2) =  0.59100D-02
        HBS(18, 1, 2) =  0.31600D-02
        HBS(18, 2, 2) =  0.59000D-03
        HBS(18, 3, 2) = -0.65700D-02
        HBS(18, 4, 2) = -0.22990D-01
        HBS(18, 5, 2) = -0.30000D-04
        HBS(18, 6, 2) = -0.20050D-01
        HBS(18, 7, 2) = -0.16000D-03
        HBS(18, 8, 2) =  0.52300D-02
        HBS(18, 9, 2) =  0.85000D-03
        HBS(18,10, 2) =  0.31230D-01
        HBS(18,11, 2) = -0.29000D-02
        HBS(18,12, 2) =  0.17300D-02
C Hessian elements for anchor point 3
        HSS( 1, 1, 3) =  0.30450D+00
        HSS( 2, 1, 3) = -0.38000D-02
        HSS( 2, 2, 3) =  0.30517D+00
        HSS( 3, 1, 3) =  0.88850D-01
        HSS( 3, 2, 3) =  0.88210D-01
        HSS( 3, 3, 3) =  0.63344D+00
        HSS( 4, 1, 3) =  0.83900D-02
        HSS( 4, 2, 3) =  0.19510D-01
        HSS( 4, 3, 3) = -0.15960D-01
        HSS( 4, 4, 3) =  0.43559D+00
        HSS( 5, 1, 3) =  0.15900D-02
        HSS( 5, 2, 3) =  0.10000D-04
        HSS( 5, 3, 3) = -0.79000D-03
        HSS( 5, 4, 3) =  0.28500D-02
        HSS( 5, 5, 3) =  0.36540D+00
        HSS( 6, 1, 3) =  0.33900D-02
        HSS( 6, 2, 3) =  0.49580D-01
        HSS( 6, 3, 3) =  0.17260D-01
        HSS( 6, 4, 3) =  0.43530D-01
        HSS( 6, 5, 3) = -0.70000D-04
        HSS( 6, 6, 3) =  0.35132D+00
        HSS( 7, 1, 3) = -0.12200D-02
        HSS( 7, 2, 3) = -0.41000D-03
        HSS( 7, 3, 3) =  0.12400D-02
        HSS( 7, 4, 3) =  0.52300D-02
        HSS( 7, 5, 3) =  0.38000D-03
        HSS( 7, 6, 3) =  0.38200D-02
        HSS( 7, 7, 3) =  0.36408D+00
        HSS( 8, 1, 3) =  0.49330D-01
        HSS( 8, 2, 3) =  0.31600D-02
        HSS( 8, 3, 3) =  0.17100D-01
        HSS( 8, 4, 3) = -0.51700D-02
        HSS( 8, 5, 3) =  0.13000D-03
        HSS( 8, 6, 3) =  0.29020D-01
        HSS( 8, 7, 3) = -0.29000D-03
        HSS( 8, 8, 3) =  0.35081D+00
        HSS( 9, 1, 3) = -0.28000D-03
        HSS( 9, 2, 3) = -0.28000D-03
        HSS( 9, 3, 3) = -0.38000D-03
        HSS( 9, 4, 3) = -0.48000D-03
        HSS( 9, 5, 3) =  0.60000D-04
        HSS( 9, 6, 3) =  0.41500D-02
        HSS( 9, 7, 3) =  0.40000D-03
        HSS( 9, 8, 3) =  0.41400D-02
        HSS( 9, 9, 3) =  0.36450D+00
        HSS(10, 1, 3) =  0.19190D-01
        HSS(10, 2, 3) =  0.79700D-02
        HSS(10, 3, 3) = -0.16090D-01
        HSS(10, 4, 3) =  0.80800D-01
        HSS(10, 5, 3) = -0.65000D-03
        HSS(10, 6, 3) = -0.49200D-02
        HSS(10, 7, 3) = -0.53000D-03
        HSS(10, 8, 3) =  0.43630D-01
        HSS(10, 9, 3) = -0.50000D-03
        HSS(10,10, 3) =  0.43670D+00
        HSS(11, 1, 3) = -0.38000D-03
        HSS(11, 2, 3) = -0.12000D-02
        HSS(11, 3, 3) =  0.12500D-02
        HSS(11, 4, 3) = -0.56000D-03
        HSS(11, 5, 3) =  0.60000D-04
        HSS(11, 6, 3) = -0.36000D-03
        HSS(11, 7, 3) =  0.19000D-03
        HSS(11, 8, 3) =  0.38300D-02
        HSS(11, 9, 3) =  0.40000D-03
        HSS(11,10, 3) =  0.51900D-02
        HSS(11,11, 3) =  0.36404D+00
        HSS(12, 1, 3) =  0.24000D-03
        HSS(12, 2, 3) =  0.21300D-02
        HSS(12, 3, 3) =  0.70000D-04
        HSS(12, 4, 3) = -0.68000D-03
        HSS(12, 5, 3) =  0.14000D-03
        HSS(12, 6, 3) = -0.18000D-03
        HSS(12, 7, 3) =  0.70000D-04
        HSS(12, 8, 3) = -0.40000D-04
        HSS(12, 9, 3) =  0.70000D-04
        HSS(12,10, 3) =  0.28100D-02
        HSS(12,11, 3) =  0.41000D-03
        HSS(12,12, 3) =  0.36629D+00
C Hessian elements for anchor point 3
        HBB( 1, 1, 3) =  0.81840D-01
        HBB( 2, 1, 3) = -0.43040D-01
        HBB( 2, 2, 3) =  0.14305D+00
        HBB( 3, 1, 3) = -0.38800D-01
        HBB( 3, 2, 3) = -0.10001D+00
        HBB( 3, 3, 3) =  0.13881D+00
        HBB( 4, 1, 3) = -0.37990D-01
        HBB( 4, 2, 3) =  0.26100D-01
        HBB( 4, 3, 3) =  0.11900D-01
        HBB( 4, 4, 3) =  0.77580D-01
        HBB( 5, 1, 3) =  0.23090D-01
        HBB( 5, 2, 3) = -0.18080D-01
        HBB( 5, 3, 3) = -0.50100D-02
        HBB( 5, 4, 3) = -0.38630D-01
        HBB( 5, 5, 3) =  0.74170D-01
        HBB( 6, 1, 3) =  0.14900D-01
        HBB( 6, 2, 3) = -0.80200D-02
        HBB( 6, 3, 3) = -0.68800D-02
        HBB( 6, 4, 3) = -0.38950D-01
        HBB( 6, 5, 3) = -0.35540D-01
        HBB( 6, 6, 3) =  0.74490D-01
        HBB( 7, 1, 3) = -0.82000D-02
        HBB( 7, 2, 3) = -0.62200D-02
        HBB( 7, 3, 3) =  0.14420D-01
        HBB( 7, 4, 3) = -0.29360D-01
        HBB( 7, 5, 3) =  0.11220D-01
        HBB( 7, 6, 3) =  0.18140D-01
        HBB( 7, 7, 3) =  0.74660D-01
        HBB( 8, 1, 3) = -0.15800D-02
        HBB( 8, 2, 3) =  0.80800D-02
        HBB( 8, 3, 3) = -0.65000D-02
        HBB( 8, 4, 3) =  0.17600D-01
        HBB( 8, 5, 3) = -0.61400D-02
        HBB( 8, 6, 3) = -0.11460D-01
        HBB( 8, 7, 3) = -0.36880D-01
        HBB( 8, 8, 3) =  0.78170D-01
        HBB( 9, 1, 3) =  0.97800D-02
        HBB( 9, 2, 3) = -0.18600D-02
        HBB( 9, 3, 3) = -0.79200D-02
        HBB( 9, 4, 3) =  0.11760D-01
        HBB( 9, 5, 3) = -0.50800D-02
        HBB( 9, 6, 3) = -0.66800D-02
        HBB( 9, 7, 3) = -0.37780D-01
        HBB( 9, 8, 3) = -0.41290D-01
        HBB( 9, 9, 3) =  0.79060D-01
        HBB(10, 1, 3) =  0.10690D-01
        HBB(10, 2, 3) = -0.58300D-02
        HBB(10, 3, 3) = -0.48700D-02
        HBB(10, 4, 3) = -0.91900D-02
        HBB(10, 5, 3) =  0.10240D-01
        HBB(10, 6, 3) = -0.10500D-02
        HBB(10, 7, 3) = -0.34370D-01
        HBB(10, 8, 3) =  0.13850D-01
        HBB(10, 9, 3) =  0.20520D-01
        HBB(10,10, 3) =  0.76450D-01
        HBB(11, 1, 3) = -0.53700D-02
        HBB(11, 2, 3) =  0.34700D-02
        HBB(11, 3, 3) =  0.19000D-02
        HBB(11, 4, 3) = -0.11400D-02
        HBB(11, 5, 3) = -0.35400D-02
        HBB(11, 6, 3) =  0.46900D-02
        HBB(11, 7, 3) =  0.20260D-01
        HBB(11, 8, 3) = -0.74900D-02
        HBB(11, 9, 3) = -0.12770D-01
        HBB(11,10, 3) = -0.38210D-01
        HBB(11,11, 3) =  0.79260D-01
        HBB(12, 1, 3) = -0.53200D-02
        HBB(12, 2, 3) =  0.23500D-02
        HBB(12, 3, 3) =  0.29700D-02
        HBB(12, 4, 3) =  0.10330D-01
        HBB(12, 5, 3) = -0.67000D-02
        HBB(12, 6, 3) = -0.36300D-02
        HBB(12, 7, 3) =  0.14110D-01
        HBB(12, 8, 3) = -0.63600D-02
        HBB(12, 9, 3) = -0.77500D-02
        HBB(12,10, 3) = -0.38240D-01
        HBB(12,11, 3) = -0.41050D-01
        HBB(12,12, 3) =  0.79290D-01
        HBB(13, 1, 3) = -0.77600D-02
        HBB(13, 2, 3) =  0.14880D-01
        HBB(13, 3, 3) = -0.71200D-02
        HBB(13, 4, 3) =  0.36800D-02
        HBB(13, 5, 3) = -0.17700D-02
        HBB(13, 6, 3) = -0.19200D-02
        HBB(13, 7, 3) = -0.65800D-02
        HBB(13, 8, 3) =  0.88500D-02
        HBB(13, 9, 3) = -0.22700D-02
        HBB(13,10, 3) = -0.34350D-01
        HBB(13,11, 3) =  0.14110D-01
        HBB(13,12, 3) =  0.20240D-01
        HBB(13,13, 3) =  0.74690D-01
        HBB(14, 1, 3) =  0.95100D-02
        HBB(14, 2, 3) = -0.81200D-02
        HBB(14, 3, 3) = -0.13900D-02
        HBB(14, 4, 3) = -0.19200D-02
        HBB(14, 5, 3) =  0.75000D-03
        HBB(14, 6, 3) =  0.11700D-02
        HBB(14, 7, 3) = -0.22400D-02
        HBB(14, 8, 3) = -0.29900D-02
        HBB(14, 9, 3) =  0.52300D-02
        HBB(14,10, 3) =  0.20540D-01
        HBB(14,11, 3) = -0.77700D-02
        HBB(14,12, 3) = -0.12770D-01
        HBB(14,13, 3) = -0.37840D-01
        HBB(14,14, 3) =  0.79020D-01
        HBB(15, 1, 3) = -0.17500D-02
        HBB(15, 2, 3) = -0.67600D-02
        HBB(15, 3, 3) =  0.85100D-02
        HBB(15, 4, 3) = -0.17600D-02
        HBB(15, 5, 3) =  0.10100D-02
        HBB(15, 6, 3) =  0.75000D-03
        HBB(15, 7, 3) =  0.88200D-02
        HBB(15, 8, 3) = -0.58500D-02
        HBB(15, 9, 3) = -0.29600D-02
        HBB(15,10, 3) =  0.13810D-01
        HBB(15,11, 3) = -0.63400D-02
        HBB(15,12, 3) = -0.74700D-02
        HBB(15,13, 3) = -0.36860D-01
        HBB(15,14, 3) = -0.41180D-01
        HBB(15,15, 3) =  0.78040D-01
        HBB(16, 1, 3) = -0.38580D-01
        HBB(16, 2, 3) =  0.14110D-01
        HBB(16, 3, 3) =  0.24470D-01
        HBB(16, 4, 3) = -0.47200D-02
        HBB(16, 5, 3) = -0.41600D-02
        HBB(16, 6, 3) =  0.88700D-02
        HBB(16, 7, 3) =  0.38500D-02
        HBB(16, 8, 3) = -0.18400D-02
        HBB(16, 9, 3) = -0.20100D-02
        HBB(16,10, 3) = -0.92300D-02
        HBB(16,11, 3) =  0.10360D-01
        HBB(16,12, 3) = -0.11200D-02
        HBB(16,13, 3) = -0.29700D-01
        HBB(16,14, 3) =  0.11950D-01
        HBB(16,15, 3) =  0.17750D-01
        HBB(16,16, 3) =  0.78380D-01
        HBB(17, 1, 3) =  0.23970D-01
        HBB(17, 2, 3) = -0.47800D-02
        HBB(17, 3, 3) = -0.19190D-01
        HBB(17, 4, 3) = -0.44800D-02
        HBB(17, 5, 3) =  0.69900D-02
        HBB(17, 6, 3) = -0.25200D-02
        HBB(17, 7, 3) = -0.18900D-02
        HBB(17, 8, 3) =  0.10700D-02
        HBB(17, 9, 3) =  0.82000D-03
        HBB(17,10, 3) =  0.10360D-01
        HBB(17,11, 3) = -0.67400D-02
        HBB(17,12, 3) = -0.36100D-02
        HBB(17,13, 3) =  0.11220D-01
        HBB(17,14, 3) = -0.51000D-02
        HBB(17,15, 3) = -0.61200D-02
        HBB(17,16, 3) = -0.39180D-01
        HBB(17,17, 3) =  0.74630D-01
        HBB(18, 1, 3) =  0.14610D-01
        HBB(18, 2, 3) = -0.93300D-02
        HBB(18, 3, 3) = -0.52800D-02
        HBB(18, 4, 3) =  0.91900D-02
        HBB(18, 5, 3) = -0.28400D-02
        HBB(18, 6, 3) = -0.63600D-02
        HBB(18, 7, 3) = -0.19600D-02
        HBB(18, 8, 3) =  0.77000D-03
        HBB(18, 9, 3) =  0.11900D-02
        HBB(18,10, 3) = -0.11200D-02
        HBB(18,11, 3) = -0.36100D-02
        HBB(18,12, 3) =  0.47400D-02
        HBB(18,13, 3) =  0.18470D-01
        HBB(18,14, 3) = -0.68500D-02
        HBB(18,15, 3) = -0.11630D-01
        HBB(18,16, 3) = -0.39190D-01
        HBB(18,17, 3) = -0.35450D-01
        HBB(18,18, 3) =  0.74640D-01
C Hessian elements for anchor point 3
        HTT( 1, 1, 3) =  0.86000D-02
        HTT( 2, 1, 3) =  0.34800D-02
        HTT( 2, 2, 3) =  0.76600D-02
        HTT( 3, 1, 3) = -0.57000D-03
        HTT( 3, 2, 3) = -0.37600D-02
        HTT( 3, 3, 3) =  0.75500D-02
        HTT( 4, 1, 3) = -0.56900D-02
        HTT( 4, 2, 3) =  0.42000D-03
        HTT( 4, 3, 3) =  0.43600D-02
        HTT( 4, 4, 3) =  0.10470D-01
        HTT( 5, 1, 3) = -0.61300D-02
        HTT( 5, 2, 3) = -0.40700D-02
        HTT( 5, 3, 3) =  0.29900D-02
        HTT( 5, 4, 3) =  0.50400D-02
        HTT( 5, 5, 3) =  0.84500D-02
        HTT( 6, 1, 3) = -0.42600D-02
        HTT( 6, 2, 3) = -0.28700D-02
        HTT( 6, 3, 3) =  0.30200D-02
        HTT( 6, 4, 3) =  0.44100D-02
        HTT( 6, 5, 3) =  0.35900D-02
        HTT( 6, 6, 3) =  0.78600D-02
        HTT( 7, 1, 3) =  0.30000D-02
        HTT( 7, 2, 3) =  0.31400D-02
        HTT( 7, 3, 3) = -0.51000D-02
        HTT( 7, 4, 3) = -0.49600D-02
        HTT( 7, 5, 3) = -0.62000D-03
        HTT( 7, 6, 3) = -0.36600D-02
        HTT( 7, 7, 3) =  0.74400D-02
        HTT( 8, 1, 3) =  0.48700D-02
        HTT( 8, 2, 3) =  0.43400D-02
        HTT( 8, 3, 3) = -0.50700D-02
        HTT( 8, 4, 3) = -0.56000D-02
        HTT( 8, 5, 3) = -0.54900D-02
        HTT( 8, 6, 3) =  0.61000D-03
        HTT( 8, 7, 3) =  0.44100D-02
        HTT( 8, 8, 3) =  0.10510D-01
        HTT( 9, 1, 3) = -0.58900D-02
        HTT( 9, 2, 3) = -0.93000D-03
        HTT( 9, 3, 3) = -0.18100D-02
        HTT( 9, 4, 3) =  0.31600D-02
        HTT( 9, 5, 3) =  0.11100D-02
        HTT( 9, 6, 3) =  0.22400D-02
        HTT( 9, 7, 3) = -0.29500D-02
        HTT( 9, 8, 3) = -0.18200D-02
        HTT( 9, 9, 3) =  0.80600D-02
        HTT(10, 1, 3) = -0.39400D-02
        HTT(10, 2, 3) =  0.44000D-03
        HTT(10, 3, 3) = -0.58000D-03
        HTT(10, 4, 3) =  0.37900D-02
        HTT(10, 5, 3) =  0.23700D-02
        HTT(10, 6, 3) =  0.17100D-02
        HTT(10, 7, 3) = -0.97000D-03
        HTT(10, 8, 3) = -0.16300D-02
        HTT(10, 9, 3) =  0.28200D-02
        HTT(10,10, 3) =  0.84500D-02
        HTT(11, 1, 3) = -0.51000D-03
        HTT(11, 2, 3) = -0.53200D-02
        HTT(11, 3, 3) =  0.15400D-02
        HTT(11, 4, 3) = -0.32600D-02
        HTT(11, 5, 3) = -0.10500D-02
        HTT(11, 6, 3) =  0.78000D-03
        HTT(11, 7, 3) = -0.31000D-02
        HTT(11, 8, 3) = -0.12700D-02
        HTT(11, 9, 3) =  0.28400D-02
        HTT(11,10, 3) = -0.17800D-02
        HTT(11,11, 3) =  0.78900D-02
        HTT(12, 1, 3) =  0.14400D-02
        HTT(12, 2, 3) = -0.39600D-02
        HTT(12, 3, 3) =  0.27700D-02
        HTT(12, 4, 3) = -0.26300D-02
        HTT(12, 5, 3) =  0.21000D-03
        HTT(12, 6, 3) =  0.25000D-03
        HTT(12, 7, 3) = -0.11200D-02
        HTT(12, 8, 3) = -0.10700D-02
        HTT(12, 9, 3) = -0.24000D-02
        HTT(12,10, 3) =  0.38400D-02
        HTT(12,11, 3) =  0.32800D-02
        HTT(12,12, 3) =  0.95300D-02
        HTT(13, 1, 3) =  0.46000D-03
        HTT(13, 2, 3) = -0.12300D-02
        HTT(13, 3, 3) =  0.19200D-02
        HTT(13, 4, 3) =  0.23000D-03
        HTT(13, 5, 3) =  0.19500D-02
        HTT(13, 6, 3) =  0.60000D-03
        HTT(13, 7, 3) =  0.50000D-03
        HTT(13, 8, 3) = -0.85000D-03
        HTT(13, 9, 3) = -0.54600D-02
        HTT(13,10, 3) = -0.30000D-04
        HTT(13,11, 3) = -0.36800D-02
        HTT(13,12, 3) =  0.17400D-02
        HTT(13,13, 3) =  0.82300D-02
        HTT(14, 1, 3) =  0.18200D-02
        HTT(14, 2, 3) =  0.48000D-03
        HTT(14, 3, 3) =  0.11000D-02
        HTT(14, 4, 3) = -0.24000D-03
        HTT(14, 5, 3) =  0.60000D-03
        HTT(14, 6, 3) = -0.11700D-02
        HTT(14, 7, 3) =  0.13200D-02
        HTT(14, 8, 3) = -0.45000D-03
        HTT(14, 9, 3) = -0.31800D-02
        HTT(14,10, 3) =  0.14500D-02
        HTT(14,11, 3) = -0.17700D-02
        HTT(14,12, 3) =  0.28600D-02
        HTT(14,13, 3) =  0.21600D-02
        HTT(14,14, 3) =  0.73500D-02
        HTT(15, 1, 3) = -0.14800D-02
        HTT(15, 2, 3) = -0.25800D-02
        HTT(15, 3, 3) =  0.70000D-03
        HTT(15, 4, 3) = -0.40000D-03
        HTT(15, 5, 3) =  0.70000D-03
        HTT(15, 6, 3) =  0.11300D-02
        HTT(15, 7, 3) = -0.14600D-02
        HTT(15, 8, 3) = -0.10400D-02
        HTT(15, 9, 3) = -0.26000D-03
        HTT(15,10, 3) = -0.56100D-02
        HTT(15,11, 3) =  0.90000D-03
        HTT(15,12, 3) = -0.44500D-02
        HTT(15,13, 3) =  0.28500D-02
        HTT(15,14, 3) = -0.24300D-02
        HTT(15,15, 3) =  0.81500D-02
        HTT(16, 1, 3) = -0.11000D-03
        HTT(16, 2, 3) = -0.87000D-03
        HTT(16, 3, 3) = -0.11000D-03
        HTT(16, 4, 3) = -0.87000D-03
        HTT(16, 5, 3) = -0.64000D-03
        HTT(16, 6, 3) = -0.64000D-03
        HTT(16, 7, 3) = -0.64000D-03
        HTT(16, 8, 3) = -0.64000D-03
        HTT(16, 9, 3) =  0.20100D-02
        HTT(16,10, 3) = -0.41300D-02
        HTT(16,11, 3) =  0.28100D-02
        HTT(16,12, 3) = -0.33300D-02
        HTT(16,13, 3) = -0.32200D-02
        HTT(16,14, 3) =  0.27600D-02
        HTT(16,15, 3) =  0.28700D-02
        HTT(16,16, 3) =  0.88500D-02
        HTT(17, 1, 3) =  0.20300D-02
        HTT(17, 2, 3) =  0.64000D-03
        HTT(17, 3, 3) =  0.51000D-03
        HTT(17, 4, 3) = -0.88000D-03
        HTT(17, 5, 3) =  0.38000D-03
        HTT(17, 6, 3) = -0.12800D-02
        HTT(17, 7, 3) =  0.19000D-02
        HTT(17, 8, 3) =  0.24000D-03
        HTT(17, 9, 3) =  0.64000D-03
        HTT(17,10, 3) = -0.15500D-02
        HTT(17,11, 3) =  0.21000D-02
        HTT(17,12, 3) = -0.80000D-04
        HTT(17,13, 3) = -0.58000D-02
        HTT(17,14, 3) =  0.29000D-03
        HTT(17,15, 3) = -0.36300D-02
        HTT(17,16, 3) =  0.24600D-02
        HTT(17,17, 3) =  0.82200D-02
        HTT(18, 1, 3) =  0.77000D-03
        HTT(18, 2, 3) =  0.11200D-02
        HTT(18, 3, 3) = -0.14100D-02
        HTT(18, 4, 3) = -0.10600D-02
        HTT(18, 5, 3) = -0.15200D-02
        HTT(18, 6, 3) = -0.26700D-02
        HTT(18, 7, 3) =  0.66000D-03
        HTT(18, 8, 3) = -0.50000D-03
        HTT(18, 9, 3) =  0.17500D-02
        HTT(18,10, 3) =  0.41000D-03
        HTT(18,11, 3) =  0.13800D-02
        HTT(18,12, 3) =  0.40000D-04
        HTT(18,13, 3) = -0.36300D-02
        HTT(18,14, 3) =  0.16800D-02
        HTT(18,15, 3) = -0.23000D-02
        HTT(18,16, 3) =  0.30100D-02
        HTT(18,17, 3) =  0.28700D-02
        HTT(18,18, 3) =  0.81700D-02
        HTT(19, 1, 3) =  0.67000D-03
        HTT(19, 2, 3) = -0.10700D-02
        HTT(19, 3, 3) =  0.13300D-02
        HTT(19, 4, 3) = -0.41000D-03
        HTT(19, 5, 3) =  0.17300D-02
        HTT(19, 6, 3) =  0.49000D-03
        HTT(19, 7, 3) =  0.10800D-02
        HTT(19, 8, 3) = -0.16000D-03
        HTT(19, 9, 3) = -0.16400D-02
        HTT(19,10, 3) = -0.30300D-02
        HTT(19,11, 3) =  0.19000D-03
        HTT(19,12, 3) = -0.12000D-02
        HTT(19,13, 3) =  0.27000D-03
        HTT(19,14, 3) = -0.49100D-02
        HTT(19,15, 3) =  0.16600D-02
        HTT(19,16, 3) = -0.35200D-02
        HTT(19,17, 3) =  0.21400D-02
        HTT(19,18, 3) = -0.24400D-02
        HTT(19,19, 3) =  0.73200D-02
        HTT(20, 1, 3) = -0.59000D-03
        HTT(20, 2, 3) = -0.59000D-03
        HTT(20, 3, 3) = -0.60000D-03
        HTT(20, 4, 3) = -0.59000D-03
        HTT(20, 5, 3) = -0.17000D-03
        HTT(20, 6, 3) = -0.90000D-03
        HTT(20, 7, 3) = -0.16000D-03
        HTT(20, 8, 3) = -0.90000D-03
        HTT(20, 9, 3) = -0.53000D-03
        HTT(20,10, 3) = -0.10700D-02
        HTT(20,11, 3) = -0.54000D-03
        HTT(20,12, 3) = -0.10700D-02
        HTT(20,13, 3) =  0.24500D-02
        HTT(20,14, 3) = -0.35100D-02
        HTT(20,15, 3) =  0.29800D-02
        HTT(20,16, 3) = -0.29800D-02
        HTT(20,17, 3) = -0.32100D-02
        HTT(20,18, 3) =  0.28600D-02
        HTT(20,19, 3) =  0.27400D-02
        HTT(20,20, 3) =  0.88200D-02
        HTT(21, 1, 3) =  0.99000D-03
        HTT(21, 2, 3) =  0.20900D-02
        HTT(21, 3, 3) = -0.29700D-02
        HTT(21, 4, 3) = -0.18700D-02
        HTT(21, 5, 3) = -0.57000D-02
        HTT(21, 6, 3) = -0.91000D-03
        HTT(21, 7, 3) = -0.17600D-02
        HTT(21, 8, 3) =  0.30300D-02
        HTT(21, 9, 3) =  0.14100D-02
        HTT(21,10, 3) =  0.29000D-03
        HTT(21,11, 3) =  0.26000D-03
        HTT(21,12, 3) = -0.86000D-03
        HTT(21,13, 3) =  0.68000D-03
        HTT(21,14, 3) = -0.16300D-02
        HTT(21,15, 3) =  0.18000D-02
        HTT(21,16, 3) = -0.51000D-03
        HTT(21,17, 3) = -0.54200D-02
        HTT(21,18, 3) = -0.27000D-03
        HTT(21,19, 3) = -0.31100D-02
        HTT(21,20, 3) =  0.20400D-02
        HTT(21,21, 3) =  0.79000D-02
        HTT(22, 1, 3) = -0.98000D-03
        HTT(22, 2, 3) =  0.81000D-03
        HTT(22, 3, 3) = -0.30000D-02
        HTT(22, 4, 3) = -0.12000D-02
        HTT(22, 5, 3) = -0.56000D-03
        HTT(22, 6, 3) = -0.54200D-02
        HTT(22, 7, 3) =  0.14400D-02
        HTT(22, 8, 3) = -0.34100D-02
        HTT(22, 9, 3) =  0.22000D-03
        HTT(22,10, 3) =  0.98000D-03
        HTT(22,11, 3) = -0.16700D-02
        HTT(22,12, 3) = -0.91000D-03
        HTT(22,13, 3) =  0.21100D-02
        HTT(22,14, 3) =  0.24000D-03
        HTT(22,15, 3) =  0.13500D-02
        HTT(22,16, 3) = -0.51000D-03
        HTT(22,17, 3) = -0.36600D-02
        HTT(22,18, 3) =  0.95000D-03
        HTT(22,19, 3) = -0.18000D-02
        HTT(22,20, 3) =  0.28200D-02
        HTT(22,21, 3) =  0.28500D-02
        HTT(22,22, 3) =  0.79700D-02
        HTT(23, 1, 3) =  0.22700D-02
        HTT(23, 2, 3) =  0.16000D-02
        HTT(23, 3, 3) = -0.10300D-02
        HTT(23, 4, 3) = -0.16900D-02
        HTT(23, 5, 3) = -0.37900D-02
        HTT(23, 6, 3) =  0.50000D-03
        HTT(23, 7, 3) = -0.51000D-03
        HTT(23, 8, 3) =  0.37700D-02
        HTT(23, 9, 3) =  0.30000D-03
        HTT(23,10, 3) = -0.16900D-02
        HTT(23,11, 3) =  0.10000D-02
        HTT(23,12, 3) = -0.99000D-03
        HTT(23,13, 3) = -0.15100D-02
        HTT(23,14, 3) = -0.30400D-02
        HTT(23,15, 3) =  0.46000D-03
        HTT(23,16, 3) = -0.10700D-02
        HTT(23,17, 3) = -0.20000D-04
        HTT(23,18, 3) = -0.56100D-02
        HTT(23,19, 3) =  0.15100D-02
        HTT(23,20, 3) = -0.40900D-02
        HTT(23,21, 3) =  0.27100D-02
        HTT(23,22, 3) = -0.18100D-02
        HTT(23,23, 3) =  0.83500D-02
        HTT(24, 1, 3) =  0.29000D-03
        HTT(24, 2, 3) =  0.33000D-03
        HTT(24, 3, 3) = -0.10600D-02
        HTT(24, 4, 3) = -0.10200D-02
        HTT(24, 5, 3) =  0.13500D-02
        HTT(24, 6, 3) = -0.40100D-02
        HTT(24, 7, 3) =  0.27000D-02
        HTT(24, 8, 3) = -0.26700D-02
        HTT(24, 9, 3) = -0.90000D-03
        HTT(24,10, 3) = -0.10000D-02
        HTT(24,11, 3) = -0.94000D-03
        HTT(24,12, 3) = -0.10400D-02
        HTT(24,13, 3) = -0.90000D-04
        HTT(24,14, 3) = -0.11700D-02
        HTT(24,15, 3) =  0.10000D-04
        HTT(24,16, 3) = -0.10700D-02
        HTT(24,17, 3) =  0.17400D-02
        HTT(24,18, 3) = -0.43900D-02
        HTT(24,19, 3) =  0.28200D-02
        HTT(24,20, 3) = -0.33100D-02
        HTT(24,21, 3) = -0.23500D-02
        HTT(24,22, 3) =  0.33100D-02
        HTT(24,23, 3) =  0.38300D-02
        HTT(24,24, 3) =  0.94900D-02
C Hessian elements for anchor point 3
        HBS( 1, 1, 3) = -0.23260D-01
        HBS( 1, 2, 3) = -0.22950D-01
        HBS( 1, 3, 3) = -0.47950D-01
        HBS( 1, 4, 3) = -0.88000D-02
        HBS( 1, 5, 3) =  0.57100D-02
        HBS( 1, 6, 3) =  0.44410D-01
        HBS( 1, 7, 3) = -0.55000D-03
        HBS( 1, 8, 3) =  0.44440D-01
        HBS( 1, 9, 3) = -0.20000D-04
        HBS( 1,10, 3) = -0.86200D-02
        HBS( 1,11, 3) = -0.54000D-03
        HBS( 1,12, 3) =  0.53300D-02
        HBS( 2, 1, 3) =  0.35140D-01
        HBS( 2, 2, 3) = -0.98400D-02
        HBS( 2, 3, 3) =  0.27820D-01
        HBS( 2, 4, 3) =  0.60300D-02
        HBS( 2, 5, 3) = -0.65700D-02
        HBS( 2, 6, 3) = -0.28710D-01
        HBS( 2, 7, 3) =  0.12300D-02
        HBS( 2, 8, 3) = -0.16910D-01
        HBS( 2, 9, 3) =  0.50000D-04
        HBS( 2,10, 3) =  0.28400D-02
        HBS( 2,11, 3) = -0.71000D-03
        HBS( 2,12, 3) = -0.49000D-03
        HBS( 3, 1, 3) = -0.11880D-01
        HBS( 3, 2, 3) =  0.32790D-01
        HBS( 3, 3, 3) =  0.20130D-01
        HBS( 3, 4, 3) =  0.27700D-02
        HBS( 3, 5, 3) =  0.86000D-03
        HBS( 3, 6, 3) = -0.15690D-01
        HBS( 3, 7, 3) = -0.68000D-03
        HBS( 3, 8, 3) = -0.27530D-01
        HBS( 3, 9, 3) = -0.30000D-04
        HBS( 3,10, 3) =  0.57800D-02
        HBS( 3,11, 3) =  0.12500D-02
        HBS( 3,12, 3) = -0.48400D-02
        HBS( 4, 1, 3) = -0.23570D-01
        HBS( 4, 2, 3) = -0.60500D-02
        HBS( 4, 3, 3) =  0.14830D-01
        HBS( 4, 4, 3) = -0.35570D-01
        HBS( 4, 5, 3) = -0.79300D-02
        HBS( 4, 6, 3) = -0.10760D-01
        HBS( 4, 7, 3) =  0.47300D-02
        HBS( 4, 8, 3) =  0.37200D-01
        HBS( 4, 9, 3) = -0.64000D-03
        HBS( 4,10, 3) =  0.48250D-01
        HBS( 4,11, 3) =  0.25000D-03
        HBS( 4,12, 3) = -0.12100D-02
        HBS( 5, 1, 3) =  0.21940D-01
        HBS( 5, 2, 3) =  0.33300D-02
        HBS( 5, 3, 3) = -0.76700D-02
        HBS( 5, 4, 3) =  0.47000D-02
        HBS( 5, 5, 3) =  0.57400D-02
        HBS( 5, 6, 3) =  0.50100D-02
        HBS( 5, 7, 3) = -0.16700D-02
        HBS( 5, 8, 3) = -0.16330D-01
        HBS( 5, 9, 3) = -0.26000D-03
        HBS( 5,10, 3) = -0.25900D-01
        HBS( 5,11, 3) = -0.12000D-03
        HBS( 5,12, 3) =  0.12000D-02
        HBS( 6, 1, 3) =  0.16300D-02
        HBS( 6, 2, 3) =  0.27200D-02
        HBS( 6, 3, 3) = -0.71600D-02
        HBS( 6, 4, 3) =  0.30860D-01
        HBS( 6, 5, 3) =  0.21900D-02
        HBS( 6, 6, 3) =  0.57500D-02
        HBS( 6, 7, 3) = -0.30600D-02
        HBS( 6, 8, 3) = -0.20870D-01
        HBS( 6, 9, 3) =  0.90000D-03
        HBS( 6,10, 3) = -0.22350D-01
        HBS( 6,11, 3) = -0.13000D-03
        HBS( 6,12, 3) =  0.10000D-04
        HBS( 7, 1, 3) = -0.76900D-02
        HBS( 7, 2, 3) =  0.28950D-01
        HBS( 7, 3, 3) =  0.62100D-02
        HBS( 7, 4, 3) = -0.33090D-01
        HBS( 7, 5, 3) =  0.43500D-02
        HBS( 7, 6, 3) = -0.35980D-01
        HBS( 7, 7, 3) = -0.84400D-02
        HBS( 7, 8, 3) = -0.54300D-02
        HBS( 7, 9, 3) =  0.48600D-02
        HBS( 7,10, 3) =  0.43110D-01
        HBS( 7,11, 3) = -0.75000D-03
        HBS( 7,12, 3) = -0.33000D-03
        HBS( 8, 1, 3) =  0.32000D-02
        HBS( 8, 2, 3) = -0.15810D-01
        HBS( 8, 3, 3) = -0.15400D-02
        HBS( 8, 4, 3) =  0.30210D-01
        HBS( 8, 5, 3) = -0.26700D-02
        HBS( 8, 6, 3) =  0.42000D-02
        HBS( 8, 7, 3) =  0.39400D-02
        HBS( 8, 8, 3) =  0.20700D-02
        HBS( 8, 9, 3) = -0.19800D-02
        HBS( 8,10, 3) = -0.19200D-01
        HBS( 8,11, 3) = -0.22000D-03
        HBS( 8,12, 3) =  0.19000D-03
        HBS( 9, 1, 3) =  0.44900D-02
        HBS( 9, 2, 3) = -0.13140D-01
        HBS( 9, 3, 3) = -0.46700D-02
        HBS( 9, 4, 3) =  0.28800D-02
        HBS( 9, 5, 3) = -0.16800D-02
        HBS( 9, 6, 3) =  0.31780D-01
        HBS( 9, 7, 3) =  0.45000D-02
        HBS( 9, 8, 3) =  0.33600D-02
        HBS( 9, 9, 3) = -0.28900D-02
        HBS( 9,10, 3) = -0.23910D-01
        HBS( 9,11, 3) =  0.97000D-03
        HBS( 9,12, 3) =  0.14000D-03
        HBS(10, 1, 3) =  0.31310D-01
        HBS(10, 2, 3) =  0.31460D-01
        HBS(10, 3, 3) =  0.56600D-02
        HBS(10, 4, 3) = -0.13750D-01
        HBS(10, 5, 3) = -0.28000D-03
        HBS(10, 6, 3) = -0.29490D-01
        HBS(10, 7, 3) =  0.47400D-02
        HBS(10, 8, 3) = -0.29500D-01
        HBS(10, 9, 3) = -0.84400D-02
        HBS(10,10, 3) = -0.13800D-01
        HBS(10,11, 3) =  0.47800D-02
        HBS(10,12, 3) = -0.37000D-03
        HBS(11, 1, 3) = -0.17710D-01
        HBS(11, 2, 3) = -0.13600D-01
        HBS(11, 3, 3) = -0.28200D-02
        HBS(11, 4, 3) =  0.80200D-02
        HBS(11, 5, 3) =  0.65000D-03
        HBS(11, 6, 3) =  0.27440D-01
        HBS(11, 7, 3) = -0.28800D-02
        HBS(11, 8, 3) =  0.19800D-02
        HBS(11, 9, 3) =  0.42300D-02
        HBS(11,10, 3) =  0.57200D-02
        HBS(11,11, 3) = -0.18700D-02
        HBS(11,12, 3) = -0.35000D-03
        HBS(12, 1, 3) = -0.13600D-01
        HBS(12, 2, 3) = -0.17860D-01
        HBS(12, 3, 3) = -0.28400D-02
        HBS(12, 4, 3) =  0.57400D-02
        HBS(12, 5, 3) = -0.37000D-03
        HBS(12, 6, 3) =  0.20500D-02
        HBS(12, 7, 3) = -0.18500D-02
        HBS(12, 8, 3) =  0.27520D-01
        HBS(12, 9, 3) =  0.42000D-02
        HBS(12,10, 3) =  0.80800D-02
        HBS(12,11, 3) = -0.29100D-02
        HBS(12,12, 3) =  0.72000D-03
        HBS(13, 1, 3) =  0.28850D-01
        HBS(13, 2, 3) = -0.76600D-02
        HBS(13, 3, 3) =  0.57200D-02
        HBS(13, 4, 3) =  0.43010D-01
        HBS(13, 5, 3) = -0.30000D-03
        HBS(13, 6, 3) = -0.53400D-02
        HBS(13, 7, 3) = -0.74000D-03
        HBS(13, 8, 3) = -0.35790D-01
        HBS(13, 9, 3) =  0.48800D-02
        HBS(13,10, 3) = -0.33150D-01
        HBS(13,11, 3) = -0.84400D-02
        HBS(13,12, 3) =  0.46600D-02
        HBS(14, 1, 3) = -0.13040D-01
        HBS(14, 2, 3) =  0.45600D-02
        HBS(14, 3, 3) = -0.44600D-02
        HBS(14, 4, 3) = -0.23920D-01
        HBS(14, 5, 3) =  0.13000D-03
        HBS(14, 6, 3) =  0.32800D-02
        HBS(14, 7, 3) =  0.98000D-03
        HBS(14, 8, 3) =  0.31630D-01
        HBS(14, 9, 3) = -0.29000D-02
        HBS(14,10, 3) =  0.29700D-02
        HBS(14,11, 3) =  0.43300D-02
        HBS(14,12, 3) = -0.18500D-02
        HBS(15, 1, 3) = -0.15810D-01
        HBS(15, 2, 3) =  0.31000D-02
        HBS(15, 3, 3) = -0.12600D-02
        HBS(15, 4, 3) = -0.19090D-01
        HBS(15, 5, 3) =  0.17000D-03
        HBS(15, 6, 3) =  0.20600D-02
        HBS(15, 7, 3) = -0.24000D-03
        HBS(15, 8, 3) =  0.41600D-02
        HBS(15, 9, 3) = -0.19700D-02
        HBS(15,10, 3) =  0.30180D-01
        HBS(15,11, 3) =  0.41100D-02
        HBS(15,12, 3) = -0.28100D-02
        HBS(16, 1, 3) = -0.56300D-02
        HBS(16, 2, 3) = -0.23750D-01
        HBS(16, 3, 3) =  0.15530D-01
        HBS(16, 4, 3) =  0.48210D-01
        HBS(16, 5, 3) = -0.15500D-02
        HBS(16, 6, 3) =  0.37160D-01
        HBS(16, 7, 3) =  0.26000D-03
        HBS(16, 8, 3) = -0.10920D-01
        HBS(16, 9, 3) = -0.64000D-03
        HBS(16,10, 3) = -0.35790D-01
        HBS(16,11, 3) =  0.47000D-02
        HBS(16,12, 3) = -0.80900D-02
        HBS(17, 1, 3) =  0.27600D-02
        HBS(17, 2, 3) =  0.21560D-01
        HBS(17, 3, 3) = -0.91700D-02
        HBS(17, 4, 3) = -0.25780D-01
        HBS(17, 5, 3) =  0.13700D-02
        HBS(17, 6, 3) = -0.16120D-01
        HBS(17, 7, 3) = -0.14000D-03
        HBS(17, 8, 3) =  0.52400D-02
        HBS(17, 9, 3) = -0.26000D-03
        HBS(17,10, 3) =  0.48900D-02
        HBS(17,11, 3) = -0.16800D-02
        HBS(17,12, 3) =  0.62500D-02
        HBS(18, 1, 3) =  0.28800D-02
        HBS(18, 2, 3) =  0.21900D-02
        HBS(18, 3, 3) = -0.63500D-02
        HBS(18, 4, 3) = -0.22430D-01
        HBS(18, 5, 3) =  0.18000D-03
        HBS(18, 6, 3) = -0.21030D-01
        HBS(18, 7, 3) = -0.12000D-03
        HBS(18, 8, 3) =  0.56800D-02
        HBS(18, 9, 3) =  0.91000D-03
        HBS(18,10, 3) =  0.30890D-01
        HBS(18,11, 3) = -0.30200D-02
        HBS(18,12, 3) =  0.18400D-02
C Hessian elements for anchor point 4
        HSS( 1, 1, 4) =  0.30630D+00
        HSS( 2, 1, 4) = -0.48000D-02
        HSS( 2, 2, 4) =  0.30630D+00
        HSS( 3, 1, 4) =  0.89520D-01
        HSS( 3, 2, 4) =  0.89520D-01
        HSS( 3, 3, 4) =  0.62529D+00
        HSS( 4, 1, 4) =  0.88300D-02
        HSS( 4, 2, 4) =  0.19200D-01
        HSS( 4, 3, 4) = -0.16480D-01
        HSS( 4, 4, 4) =  0.43562D+00
        HSS( 5, 1, 4) =  0.21600D-02
        HSS( 5, 2, 4) =  0.18000D-03
        HSS( 5, 3, 4) =  0.10000D-03
        HSS( 5, 4, 4) =  0.28600D-02
        HSS( 5, 5, 4) =  0.36547D+00
        HSS( 6, 1, 4) =  0.27000D-02
        HSS( 6, 2, 4) =  0.49850D-01
        HSS( 6, 3, 4) =  0.18220D-01
        HSS( 6, 4, 4) =  0.43220D-01
        HSS( 6, 5, 4) = -0.10000D-04
        HSS( 6, 6, 4) =  0.35171D+00
        HSS( 7, 1, 4) = -0.12200D-02
        HSS( 7, 2, 4) = -0.39000D-03
        HSS( 7, 3, 4) =  0.12600D-02
        HSS( 7, 4, 4) =  0.51900D-02
        HSS( 7, 5, 4) =  0.41000D-03
        HSS( 7, 6, 4) =  0.38600D-02
        HSS( 7, 7, 4) =  0.36379D+00
        HSS( 8, 1, 4) =  0.49850D-01
        HSS( 8, 2, 4) =  0.27000D-02
        HSS( 8, 3, 4) =  0.18220D-01
        HSS( 8, 4, 4) = -0.52300D-02
        HSS( 8, 5, 4) = -0.11000D-03
        HSS( 8, 6, 4) =  0.28880D-01
        HSS( 8, 7, 4) = -0.34000D-03
        HSS( 8, 8, 4) =  0.35171D+00
        HSS( 9, 1, 4) = -0.29000D-03
        HSS( 9, 2, 4) = -0.29000D-03
        HSS( 9, 3, 4) = -0.38000D-03
        HSS( 9, 4, 4) = -0.48000D-03
        HSS( 9, 5, 4) =  0.80000D-04
        HSS( 9, 6, 4) =  0.41700D-02
        HSS( 9, 7, 4) =  0.40000D-03
        HSS( 9, 8, 4) =  0.41700D-02
        HSS( 9, 9, 4) =  0.36487D+00
        HSS(10, 1, 4) =  0.19200D-01
        HSS(10, 2, 4) =  0.88300D-02
        HSS(10, 3, 4) = -0.16480D-01
        HSS(10, 4, 4) =  0.80730D-01
        HSS(10, 5, 4) = -0.72000D-03
        HSS(10, 6, 4) = -0.52300D-02
        HSS(10, 7, 4) = -0.55000D-03
        HSS(10, 8, 4) =  0.43220D-01
        HSS(10, 9, 4) = -0.48000D-03
        HSS(10,10, 4) =  0.43562D+00
        HSS(11, 1, 4) = -0.39000D-03
        HSS(11, 2, 4) = -0.12200D-02
        HSS(11, 3, 4) =  0.12600D-02
        HSS(11, 4, 4) = -0.55000D-03
        HSS(11, 5, 4) =  0.70000D-04
        HSS(11, 6, 4) = -0.34000D-03
        HSS(11, 7, 4) =  0.20000D-03
        HSS(11, 8, 4) =  0.38600D-02
        HSS(11, 9, 4) =  0.40000D-03
        HSS(11,10, 4) =  0.51900D-02
        HSS(11,11, 4) =  0.36379D+00
        HSS(12, 1, 4) =  0.18000D-03
        HSS(12, 2, 4) =  0.21600D-02
        HSS(12, 3, 4) =  0.10000D-03
        HSS(12, 4, 4) = -0.72000D-03
        HSS(12, 5, 4) =  0.17000D-03
        HSS(12, 6, 4) = -0.11000D-03
        HSS(12, 7, 4) =  0.70000D-04
        HSS(12, 8, 4) = -0.10000D-04
        HSS(12, 9, 4) =  0.80000D-04
        HSS(12,10, 4) =  0.28600D-02
        HSS(12,11, 4) =  0.41000D-03
        HSS(12,12, 4) =  0.36547D+00
C Hessian elements for anchor point 4
        HBB( 1, 1, 4) =  0.81570D-01
        HBB( 2, 1, 4) = -0.40780D-01
        HBB( 2, 2, 4) =  0.13357D+00
        HBB( 3, 1, 4) = -0.40780D-01
        HBB( 3, 2, 4) = -0.92780D-01
        HBB( 3, 3, 4) =  0.13357D+00
        HBB( 4, 1, 4) = -0.38240D-01
        HBB( 4, 2, 4) =  0.25780D-01
        HBB( 4, 3, 4) =  0.12470D-01
        HBB( 4, 4, 4) =  0.78290D-01
        HBB( 5, 1, 4) =  0.23810D-01
        HBB( 5, 2, 4) = -0.19980D-01
        HBB( 5, 3, 4) = -0.38300D-02
        HBB( 5, 4, 4) = -0.39070D-01
        HBB( 5, 5, 4) =  0.74370D-01
        HBB( 6, 1, 4) =  0.14430D-01
        HBB( 6, 2, 4) = -0.58000D-02
        HBB( 6, 3, 4) = -0.86300D-02
        HBB( 6, 4, 4) = -0.39220D-01
        HBB( 6, 5, 4) = -0.35300D-01
        HBB( 6, 6, 4) =  0.74520D-01
        HBB( 7, 1, 4) = -0.78000D-02
        HBB( 7, 2, 4) = -0.70900D-02
        HBB( 7, 3, 4) =  0.14890D-01
        HBB( 7, 4, 4) = -0.29810D-01
        HBB( 7, 5, 4) =  0.11250D-01
        HBB( 7, 6, 4) =  0.18570D-01
        HBB( 7, 7, 4) =  0.74890D-01
        HBB( 8, 1, 4) = -0.17500D-02
        HBB( 8, 2, 4) =  0.84900D-02
        HBB( 8, 3, 4) = -0.67400D-02
        HBB( 8, 4, 4) =  0.17840D-01
        HBB( 8, 5, 4) = -0.61200D-02
        HBB( 8, 6, 4) = -0.11720D-01
        HBB( 8, 7, 4) = -0.36950D-01
        HBB( 8, 8, 4) =  0.78020D-01
        HBB( 9, 1, 4) =  0.95600D-02
        HBB( 9, 2, 4) = -0.14000D-02
        HBB( 9, 3, 4) = -0.81500D-02
        HBB( 9, 4, 4) =  0.11970D-01
        HBB( 9, 5, 4) = -0.51200D-02
        HBB( 9, 6, 4) = -0.68500D-02
        HBB( 9, 7, 4) = -0.37940D-01
        HBB( 9, 8, 4) = -0.41060D-01
        HBB( 9, 9, 4) =  0.79000D-01
        HBB(10, 1, 4) =  0.10520D-01
        HBB(10, 2, 4) = -0.52600D-02
        HBB(10, 3, 4) = -0.52600D-02
        HBB(10, 4, 4) = -0.91700D-02
        HBB(10, 5, 4) =  0.10330D-01
        HBB(10, 6, 4) = -0.11700D-02
        HBB(10, 7, 4) = -0.34390D-01
        HBB(10, 8, 4) =  0.13830D-01
        HBB(10, 9, 4) =  0.20560D-01
        HBB(10,10, 4) =  0.76580D-01
        HBB(11, 1, 4) = -0.52600D-02
        HBB(11, 2, 4) =  0.31800D-02
        HBB(11, 3, 4) =  0.20800D-02
        HBB(11, 4, 4) = -0.11400D-02
        HBB(11, 5, 4) = -0.36300D-02
        HBB(11, 6, 4) =  0.47700D-02
        HBB(11, 7, 4) =  0.20190D-01
        HBB(11, 8, 4) = -0.74500D-02
        HBB(11, 9, 4) = -0.12740D-01
        HBB(11,10, 4) = -0.38290D-01
        HBB(11,11, 4) =  0.79780D-01
        HBB(12, 1, 4) = -0.52600D-02
        HBB(12, 2, 4) =  0.20800D-02
        HBB(12, 3, 4) =  0.31800D-02
        HBB(12, 4, 4) =  0.10300D-01
        HBB(12, 5, 4) = -0.67000D-02
        HBB(12, 6, 4) = -0.36000D-02
        HBB(12, 7, 4) =  0.14200D-01
        HBB(12, 8, 4) = -0.63700D-02
        HBB(12, 9, 4) = -0.78200D-02
        HBB(12,10, 4) = -0.38290D-01
        HBB(12,11, 4) = -0.41490D-01
        HBB(12,12, 4) =  0.79780D-01
        HBB(13, 1, 4) = -0.78000D-02
        HBB(13, 2, 4) =  0.14890D-01
        HBB(13, 3, 4) = -0.70900D-02
        HBB(13, 4, 4) =  0.37200D-02
        HBB(13, 5, 4) = -0.18400D-02
        HBB(13, 6, 4) = -0.18800D-02
        HBB(13, 7, 4) = -0.66100D-02
        HBB(13, 8, 4) =  0.88400D-02
        HBB(13, 9, 4) = -0.22300D-02
        HBB(13,10, 4) = -0.34390D-01
        HBB(13,11, 4) =  0.14200D-01
        HBB(13,12, 4) =  0.20190D-01
        HBB(13,13, 4) =  0.74890D-01
        HBB(14, 1, 4) =  0.95600D-02
        HBB(14, 2, 4) = -0.81500D-02
        HBB(14, 3, 4) = -0.14000D-02
        HBB(14, 4, 4) = -0.19200D-02
        HBB(14, 5, 4) =  0.79000D-03
        HBB(14, 6, 4) =  0.11300D-02
        HBB(14, 7, 4) = -0.22300D-02
        HBB(14, 8, 4) = -0.29700D-02
        HBB(14, 9, 4) =  0.52000D-02
        HBB(14,10, 4) =  0.20560D-01
        HBB(14,11, 4) = -0.78200D-02
        HBB(14,12, 4) = -0.12740D-01
        HBB(14,13, 4) = -0.37940D-01
        HBB(14,14, 4) =  0.79000D-01
        HBB(15, 1, 4) = -0.17500D-02
        HBB(15, 2, 4) = -0.67400D-02
        HBB(15, 3, 4) =  0.84900D-02
        HBB(15, 4, 4) = -0.18000D-02
        HBB(15, 5, 4) =  0.10500D-02
        HBB(15, 6, 4) =  0.75000D-03
        HBB(15, 7, 4) =  0.88400D-02
        HBB(15, 8, 4) = -0.58700D-02
        HBB(15, 9, 4) = -0.29700D-02
        HBB(15,10, 4) =  0.13830D-01
        HBB(15,11, 4) = -0.63700D-02
        HBB(15,12, 4) = -0.74500D-02
        HBB(15,13, 4) = -0.36950D-01
        HBB(15,14, 4) = -0.41060D-01
        HBB(15,15, 4) =  0.78020D-01
        HBB(16, 1, 4) = -0.38240D-01
        HBB(16, 2, 4) =  0.12470D-01
        HBB(16, 3, 4) =  0.25780D-01
        HBB(16, 4, 4) = -0.47900D-02
        HBB(16, 5, 4) = -0.44800D-02
        HBB(16, 6, 4) =  0.92700D-02
        HBB(16, 7, 4) =  0.37200D-02
        HBB(16, 8, 4) = -0.18000D-02
        HBB(16, 9, 4) = -0.19200D-02
        HBB(16,10, 4) = -0.91700D-02
        HBB(16,11, 4) =  0.10300D-01
        HBB(16,12, 4) = -0.11400D-02
        HBB(16,13, 4) = -0.29810D-01
        HBB(16,14, 4) =  0.11970D-01
        HBB(16,15, 4) =  0.17840D-01
        HBB(16,16, 4) =  0.78290D-01
        HBB(17, 1, 4) =  0.23810D-01
        HBB(17, 2, 4) = -0.38300D-02
        HBB(17, 3, 4) = -0.19980D-01
        HBB(17, 4, 4) = -0.44800D-02
        HBB(17, 5, 4) =  0.71700D-02
        HBB(17, 6, 4) = -0.26900D-02
        HBB(17, 7, 4) = -0.18400D-02
        HBB(17, 8, 4) =  0.10500D-02
        HBB(17, 9, 4) =  0.79000D-03
        HBB(17,10, 4) =  0.10330D-01
        HBB(17,11, 4) = -0.67000D-02
        HBB(17,12, 4) = -0.36300D-02
        HBB(17,13, 4) =  0.11250D-01
        HBB(17,14, 4) = -0.51200D-02
        HBB(17,15, 4) = -0.61200D-02
        HBB(17,16, 4) = -0.39070D-01
        HBB(17,17, 4) =  0.74370D-01
        HBB(18, 1, 4) =  0.14430D-01
        HBB(18, 2, 4) = -0.86300D-02
        HBB(18, 3, 4) = -0.58000D-02
        HBB(18, 4, 4) =  0.92700D-02
        HBB(18, 5, 4) = -0.26900D-02
        HBB(18, 6, 4) = -0.65700D-02
        HBB(18, 7, 4) = -0.18800D-02
        HBB(18, 8, 4) =  0.75000D-03
        HBB(18, 9, 4) =  0.11300D-02
        HBB(18,10, 4) = -0.11700D-02
        HBB(18,11, 4) = -0.36000D-02
        HBB(18,12, 4) =  0.47700D-02
        HBB(18,13, 4) =  0.18570D-01
        HBB(18,14, 4) = -0.68500D-02
        HBB(18,15, 4) = -0.11720D-01
        HBB(18,16, 4) = -0.39220D-01
        HBB(18,17, 4) = -0.35300D-01
        HBB(18,18, 4) =  0.74520D-01
C Hessian elements for anchor point 4
        HTT( 1, 1, 4) =  0.86600D-02
        HTT( 2, 1, 4) =  0.36400D-02
        HTT( 2, 2, 4) =  0.78600D-02
        HTT( 3, 1, 4) = -0.53000D-03
        HTT( 3, 2, 4) = -0.36100D-02
        HTT( 3, 3, 4) =  0.74400D-02
        HTT( 4, 1, 4) = -0.55500D-02
        HTT( 4, 2, 4) =  0.61000D-03
        HTT( 4, 3, 4) =  0.43600D-02
        HTT( 4, 4, 4) =  0.10520D-01
        HTT( 5, 1, 4) = -0.62400D-02
        HTT( 5, 2, 4) = -0.42900D-02
        HTT( 5, 3, 4) =  0.29500D-02
        HTT( 5, 4, 4) =  0.49100D-02
        HTT( 5, 5, 4) =  0.86600D-02
        HTT( 6, 1, 4) = -0.42900D-02
        HTT( 6, 2, 4) = -0.29000D-02
        HTT( 6, 3, 4) =  0.29600D-02
        HTT( 6, 4, 4) =  0.43500D-02
        HTT( 6, 5, 4) =  0.36400D-02
        HTT( 6, 6, 4) =  0.78600D-02
        HTT( 7, 1, 4) =  0.29500D-02
        HTT( 7, 2, 4) =  0.29600D-02
        HTT( 7, 3, 4) = -0.50200D-02
        HTT( 7, 4, 4) = -0.50000D-02
        HTT( 7, 5, 4) = -0.53000D-03
        HTT( 7, 6, 4) = -0.36100D-02
        HTT( 7, 7, 4) =  0.74400D-02
        HTT( 8, 1, 4) =  0.49100D-02
        HTT( 8, 2, 4) =  0.43500D-02
        HTT( 8, 3, 4) = -0.50000D-02
        HTT( 8, 4, 4) = -0.55600D-02
        HTT( 8, 5, 4) = -0.55500D-02
        HTT( 8, 6, 4) =  0.61000D-03
        HTT( 8, 7, 4) =  0.43600D-02
        HTT( 8, 8, 4) =  0.10520D-01
        HTT( 9, 1, 4) = -0.58800D-02
        HTT( 9, 2, 4) = -0.96000D-03
        HTT( 9, 3, 4) = -0.18400D-02
        HTT( 9, 4, 4) =  0.30900D-02
        HTT( 9, 5, 4) =  0.10900D-02
        HTT( 9, 6, 4) =  0.22400D-02
        HTT( 9, 7, 4) = -0.29500D-02
        HTT( 9, 8, 4) = -0.18000D-02
        HTT( 9, 9, 4) =  0.80500D-02
        HTT(10, 1, 4) = -0.39100D-02
        HTT(10, 2, 4) =  0.46000D-03
        HTT(10, 3, 4) = -0.55000D-03
        HTT(10, 4, 4) =  0.38200D-02
        HTT(10, 5, 4) =  0.23300D-02
        HTT(10, 6, 4) =  0.17200D-02
        HTT(10, 7, 4) = -0.10300D-02
        HTT(10, 8, 4) = -0.16300D-02
        HTT(10, 9, 4) =  0.28100D-02
        HTT(10,10, 4) =  0.84200D-02
        HTT(11, 1, 4) = -0.58000D-03
        HTT(11, 2, 4) = -0.54100D-02
        HTT(11, 3, 4) =  0.14200D-02
        HTT(11, 4, 4) = -0.34100D-02
        HTT(11, 5, 4) = -0.97000D-03
        HTT(11, 6, 4) =  0.77000D-03
        HTT(11, 7, 4) = -0.29700D-02
        HTT(11, 8, 4) = -0.12200D-02
        HTT(11, 9, 4) =  0.28600D-02
        HTT(11,10, 4) = -0.18000D-02
        HTT(11,11, 4) =  0.79500D-02
        HTT(12, 1, 4) =  0.13900D-02
        HTT(12, 2, 4) = -0.39900D-02
        HTT(12, 3, 4) =  0.27100D-02
        HTT(12, 4, 4) = -0.26800D-02
        HTT(12, 5, 4) =  0.27000D-03
        HTT(12, 6, 4) =  0.26000D-03
        HTT(12, 7, 4) = -0.10500D-02
        HTT(12, 8, 4) = -0.10500D-02
        HTT(12, 9, 4) = -0.23900D-02
        HTT(12,10, 4) =  0.38100D-02
        HTT(12,11, 4) =  0.32900D-02
        HTT(12,12, 4) =  0.94800D-02
        HTT(13, 1, 4) =  0.42000D-03
        HTT(13, 2, 4) = -0.12700D-02
        HTT(13, 3, 4) =  0.19200D-02
        HTT(13, 4, 4) =  0.23000D-03
        HTT(13, 5, 4) =  0.20100D-02
        HTT(13, 6, 4) =  0.61000D-03
        HTT(13, 7, 4) =  0.51000D-03
        HTT(13, 8, 4) = -0.88000D-03
        HTT(13, 9, 4) = -0.54500D-02
        HTT(13,10, 4) = -0.30000D-04
        HTT(13,11, 4) = -0.36700D-02
        HTT(13,12, 4) =  0.17500D-02
        HTT(13,13, 4) =  0.82400D-02
        HTT(14, 1, 4) =  0.18200D-02
        HTT(14, 2, 4) =  0.50000D-03
        HTT(14, 3, 4) =  0.11400D-02
        HTT(14, 4, 4) = -0.18000D-03
        HTT(14, 5, 4) =  0.61000D-03
        HTT(14, 6, 4) = -0.11500D-02
        HTT(14, 7, 4) =  0.12900D-02
        HTT(14, 8, 4) = -0.47000D-03
        HTT(14, 9, 4) = -0.31900D-02
        HTT(14,10, 4) =  0.14500D-02
        HTT(14,11, 4) = -0.18000D-02
        HTT(14,12, 4) =  0.28400D-02
        HTT(14,13, 4) =  0.21800D-02
        HTT(14,14, 4) =  0.73800D-02
        HTT(15, 1, 4) = -0.15300D-02
        HTT(15, 2, 4) = -0.26700D-02
        HTT(15, 3, 4) =  0.64000D-03
        HTT(15, 4, 4) = -0.50000D-03
        HTT(15, 5, 4) =  0.78000D-03
        HTT(15, 6, 4) =  0.11300D-02
        HTT(15, 7, 4) = -0.13900D-02
        HTT(15, 8, 4) = -0.10500D-02
        HTT(15, 9, 4) = -0.26000D-03
        HTT(15,10, 4) = -0.55900D-02
        HTT(15,11, 4) =  0.94000D-03
        HTT(15,12, 4) = -0.43900D-02
        HTT(15,13, 4) =  0.28700D-02
        HTT(15,14, 4) = -0.24200D-02
        HTT(15,15, 4) =  0.81600D-02
        HTT(16, 1, 4) = -0.14000D-03
        HTT(16, 2, 4) = -0.90000D-03
        HTT(16, 3, 4) = -0.13000D-03
        HTT(16, 4, 4) = -0.90000D-03
        HTT(16, 5, 4) = -0.62000D-03
        HTT(16, 6, 4) = -0.64000D-03
        HTT(16, 7, 4) = -0.62000D-03
        HTT(16, 8, 4) = -0.64000D-03
        HTT(16, 9, 4) =  0.20000D-02
        HTT(16,10, 4) = -0.41100D-02
        HTT(16,11, 4) =  0.28100D-02
        HTT(16,12, 4) = -0.32900D-02
        HTT(16,13, 4) = -0.31800D-02
        HTT(16,14, 4) =  0.27800D-02
        HTT(16,15, 4) =  0.28700D-02
        HTT(16,16, 4) =  0.88300D-02
        HTT(17, 1, 4) =  0.20100D-02
        HTT(17, 2, 4) =  0.61000D-03
        HTT(17, 3, 4) =  0.51000D-03
        HTT(17, 4, 4) = -0.88000D-03
        HTT(17, 5, 4) =  0.42000D-03
        HTT(17, 6, 4) = -0.12700D-02
        HTT(17, 7, 4) =  0.19200D-02
        HTT(17, 8, 4) =  0.23000D-03
        HTT(17, 9, 4) =  0.64000D-03
        HTT(17,10, 4) = -0.15500D-02
        HTT(17,11, 4) =  0.21100D-02
        HTT(17,12, 4) = -0.80000D-04
        HTT(17,13, 4) = -0.58000D-02
        HTT(17,14, 4) =  0.26000D-03
        HTT(17,15, 4) = -0.36300D-02
        HTT(17,16, 4) =  0.24300D-02
        HTT(17,17, 4) =  0.82400D-02
        HTT(18, 1, 4) =  0.78000D-03
        HTT(18, 2, 4) =  0.11300D-02
        HTT(18, 3, 4) = -0.13900D-02
        HTT(18, 4, 4) = -0.10500D-02
        HTT(18, 5, 4) = -0.15300D-02
        HTT(18, 6, 4) = -0.26700D-02
        HTT(18, 7, 4) =  0.64000D-03
        HTT(18, 8, 4) = -0.50000D-03
        HTT(18, 9, 4) =  0.17400D-02
        HTT(18,10, 4) =  0.40000D-03
        HTT(18,11, 4) =  0.13800D-02
        HTT(18,12, 4) =  0.40000D-04
        HTT(18,13, 4) = -0.36300D-02
        HTT(18,14, 4) =  0.16600D-02
        HTT(18,15, 4) = -0.23000D-02
        HTT(18,16, 4) =  0.29900D-02
        HTT(18,17, 4) =  0.28700D-02
        HTT(18,18, 4) =  0.81600D-02
        HTT(19, 1, 4) =  0.61000D-03
        HTT(19, 2, 4) = -0.11500D-02
        HTT(19, 3, 4) =  0.12900D-02
        HTT(19, 4, 4) = -0.47000D-03
        HTT(19, 5, 4) =  0.18200D-02
        HTT(19, 6, 4) =  0.50000D-03
        HTT(19, 7, 4) =  0.11400D-02
        HTT(19, 8, 4) = -0.18000D-03
        HTT(19, 9, 4) = -0.16200D-02
        HTT(19,10, 4) = -0.30300D-02
        HTT(19,11, 4) =  0.24000D-03
        HTT(19,12, 4) = -0.11700D-02
        HTT(19,13, 4) =  0.26000D-03
        HTT(19,14, 4) = -0.49400D-02
        HTT(19,15, 4) =  0.16600D-02
        HTT(19,16, 4) = -0.35300D-02
        HTT(19,17, 4) =  0.21800D-02
        HTT(19,18, 4) = -0.24200D-02
        HTT(19,19, 4) =  0.73800D-02
        HTT(20, 1, 4) = -0.62000D-03
        HTT(20, 2, 4) = -0.64000D-03
        HTT(20, 3, 4) = -0.62000D-03
        HTT(20, 4, 4) = -0.64000D-03
        HTT(20, 5, 4) = -0.14000D-03
        HTT(20, 6, 4) = -0.90000D-03
        HTT(20, 7, 4) = -0.13000D-03
        HTT(20, 8, 4) = -0.90000D-03
        HTT(20, 9, 4) = -0.52000D-03
        HTT(20,10, 4) = -0.10800D-02
        HTT(20,11, 4) = -0.49000D-03
        HTT(20,12, 4) = -0.10500D-02
        HTT(20,13, 4) =  0.24300D-02
        HTT(20,14, 4) = -0.35300D-02
        HTT(20,15, 4) =  0.29900D-02
        HTT(20,16, 4) = -0.29700D-02
        HTT(20,17, 4) = -0.31800D-02
        HTT(20,18, 4) =  0.28700D-02
        HTT(20,19, 4) =  0.27800D-02
        HTT(20,20, 4) =  0.88300D-02
        HTT(21, 1, 4) =  0.10900D-02
        HTT(21, 2, 4) =  0.22400D-02
        HTT(21, 3, 4) = -0.29500D-02
        HTT(21, 4, 4) = -0.18000D-02
        HTT(21, 5, 4) = -0.58800D-02
        HTT(21, 6, 4) = -0.96000D-03
        HTT(21, 7, 4) = -0.18400D-02
        HTT(21, 8, 4) =  0.30900D-02
        HTT(21, 9, 4) =  0.14200D-02
        HTT(21,10, 4) =  0.31000D-03
        HTT(21,11, 4) =  0.21000D-03
        HTT(21,12, 4) = -0.90000D-03
        HTT(21,13, 4) =  0.64000D-03
        HTT(21,14, 4) = -0.16200D-02
        HTT(21,15, 4) =  0.17400D-02
        HTT(21,16, 4) = -0.52000D-03
        HTT(21,17, 4) = -0.54500D-02
        HTT(21,18, 4) = -0.26000D-03
        HTT(21,19, 4) = -0.31900D-02
        HTT(21,20, 4) =  0.20000D-02
        HTT(21,21, 4) =  0.80500D-02
        HTT(22, 1, 4) = -0.97000D-03
        HTT(22, 2, 4) =  0.77000D-03
        HTT(22, 3, 4) = -0.29700D-02
        HTT(22, 4, 4) = -0.12200D-02
        HTT(22, 5, 4) = -0.58000D-03
        HTT(22, 6, 4) = -0.54100D-02
        HTT(22, 7, 4) =  0.14200D-02
        HTT(22, 8, 4) = -0.34100D-02
        HTT(22, 9, 4) =  0.21000D-03
        HTT(22,10, 4) =  0.95000D-03
        HTT(22,11, 4) = -0.16300D-02
        HTT(22,12, 4) = -0.89000D-03
        HTT(22,13, 4) =  0.21100D-02
        HTT(22,14, 4) =  0.24000D-03
        HTT(22,15, 4) =  0.13800D-02
        HTT(22,16, 4) = -0.49000D-03
        HTT(22,17, 4) = -0.36700D-02
        HTT(22,18, 4) =  0.94000D-03
        HTT(22,19, 4) = -0.18000D-02
        HTT(22,20, 4) =  0.28100D-02
        HTT(22,21, 4) =  0.28600D-02
        HTT(22,22, 4) =  0.79500D-02
        HTT(23, 1, 4) =  0.23300D-02
        HTT(23, 2, 4) =  0.17200D-02
        HTT(23, 3, 4) = -0.10300D-02
        HTT(23, 4, 4) = -0.16300D-02
        HTT(23, 5, 4) = -0.39100D-02
        HTT(23, 6, 4) =  0.46000D-03
        HTT(23, 7, 4) = -0.55000D-03
        HTT(23, 8, 4) =  0.38200D-02
        HTT(23, 9, 4) =  0.31000D-03
        HTT(23,10, 4) = -0.16600D-02
        HTT(23,11, 4) =  0.95000D-03
        HTT(23,12, 4) = -0.10200D-02
        HTT(23,13, 4) = -0.15500D-02
        HTT(23,14, 4) = -0.30300D-02
        HTT(23,15, 4) =  0.40000D-03
        HTT(23,16, 4) = -0.10800D-02
        HTT(23,17, 4) = -0.30000D-04
        HTT(23,18, 4) = -0.55900D-02
        HTT(23,19, 4) =  0.14500D-02
        HTT(23,20, 4) = -0.41100D-02
        HTT(23,21, 4) =  0.28100D-02
        HTT(23,22, 4) = -0.18000D-02
        HTT(23,23, 4) =  0.84200D-02
        HTT(24, 1, 4) =  0.27000D-03
        HTT(24, 2, 4) =  0.26000D-03
        HTT(24, 3, 4) = -0.10500D-02
        HTT(24, 4, 4) = -0.10500D-02
        HTT(24, 5, 4) =  0.13900D-02
        HTT(24, 6, 4) = -0.39900D-02
        HTT(24, 7, 4) =  0.27100D-02
        HTT(24, 8, 4) = -0.26800D-02
        HTT(24, 9, 4) = -0.90000D-03
        HTT(24,10, 4) = -0.10200D-02
        HTT(24,11, 4) = -0.89000D-03
        HTT(24,12, 4) = -0.10100D-02
        HTT(24,13, 4) = -0.80000D-04
        HTT(24,14, 4) = -0.11700D-02
        HTT(24,15, 4) =  0.40000D-04
        HTT(24,16, 4) = -0.10500D-02
        HTT(24,17, 4) =  0.17500D-02
        HTT(24,18, 4) = -0.43900D-02
        HTT(24,19, 4) =  0.28400D-02
        HTT(24,20, 4) = -0.32900D-02
        HTT(24,21, 4) = -0.23900D-02
        HTT(24,22, 4) =  0.32900D-02
        HTT(24,23, 4) =  0.38100D-02
        HTT(24,24, 4) =  0.94800D-02
C Hessian elements for anchor point 4
        HBS( 1, 1, 4) = -0.22700D-01
        HBS( 1, 2, 4) = -0.22700D-01
        HBS( 1, 3, 4) = -0.46990D-01
        HBS( 1, 4, 4) = -0.88700D-02
        HBS( 1, 5, 4) =  0.54300D-02
        HBS( 1, 6, 4) =  0.44400D-01
        HBS( 1, 7, 4) = -0.55000D-03
        HBS( 1, 8, 4) =  0.44400D-01
        HBS( 1, 9, 4) = -0.10000D-04
        HBS( 1,10, 4) = -0.88700D-02
        HBS( 1,11, 4) = -0.55000D-03
        HBS( 1,12, 4) =  0.54300D-02
        HBS( 2, 1, 4) =  0.32950D-01
        HBS( 2, 2, 4) = -0.10250D-01
        HBS( 2, 3, 4) =  0.23500D-01
        HBS( 2, 4, 4) =  0.60000D-02
        HBS( 2, 5, 4) = -0.49200D-02
        HBS( 2, 6, 4) = -0.28130D-01
        HBS( 2, 7, 4) =  0.12500D-02
        HBS( 2, 8, 4) = -0.16270D-01
        HBS( 2, 9, 4) =  0.00000D+00
        HBS( 2,10, 4) =  0.28700D-02
        HBS( 2,11, 4) = -0.70000D-03
        HBS( 2,12, 4) = -0.51000D-03
        HBS( 3, 1, 4) = -0.10250D-01
        HBS( 3, 2, 4) =  0.32950D-01
        HBS( 3, 3, 4) =  0.23500D-01
        HBS( 3, 4, 4) =  0.28700D-02
        HBS( 3, 5, 4) = -0.51000D-03
        HBS( 3, 6, 4) = -0.16270D-01
        HBS( 3, 7, 4) = -0.70000D-03
        HBS( 3, 8, 4) = -0.28130D-01
        HBS( 3, 9, 4) =  0.00000D+00
        HBS( 3,10, 4) =  0.60000D-02
        HBS( 3,11, 4) =  0.12500D-02
        HBS( 3,12, 4) = -0.49200D-02
        HBS( 4, 1, 4) = -0.24180D-01
        HBS( 4, 2, 4) = -0.60700D-02
        HBS( 4, 3, 4) =  0.14620D-01
        HBS( 4, 4, 4) = -0.35760D-01
        HBS( 4, 5, 4) = -0.81300D-02
        HBS( 4, 6, 4) = -0.10830D-01
        HBS( 4, 7, 4) =  0.47400D-02
        HBS( 4, 8, 4) =  0.37380D-01
        HBS( 4, 9, 4) = -0.65000D-03
        HBS( 4,10, 4) =  0.48130D-01
        HBS( 4,11, 4) =  0.25000D-03
        HBS( 4,12, 4) = -0.12600D-02
        HBS( 5, 1, 4) =  0.21630D-01
        HBS( 5, 2, 4) =  0.30300D-02
        HBS( 5, 3, 4) = -0.86900D-02
        HBS( 5, 4, 4) =  0.48300D-02
        HBS( 5, 5, 4) =  0.60700D-02
        HBS( 5, 6, 4) =  0.52000D-02
        HBS( 5, 7, 4) = -0.17000D-02
        HBS( 5, 8, 4) = -0.16220D-01
        HBS( 5, 9, 4) = -0.27000D-03
        HBS( 5,10, 4) = -0.25720D-01
        HBS( 5,11, 4) = -0.13000D-03
        HBS( 5,12, 4) =  0.12100D-02
        HBS( 6, 1, 4) =  0.25500D-02
        HBS( 6, 2, 4) =  0.30300D-02
        HBS( 6, 3, 4) = -0.59400D-02
        HBS( 6, 4, 4) =  0.30930D-01
        HBS( 6, 5, 4) =  0.20600D-02
        HBS( 6, 6, 4) =  0.56300D-02
        HBS( 6, 7, 4) = -0.30400D-02
        HBS( 6, 8, 4) = -0.21160D-01
        HBS( 6, 9, 4) =  0.91000D-03
        HBS( 6,10, 4) = -0.22410D-01
        HBS( 6,11, 4) = -0.12000D-03
        HBS( 6,12, 4) =  0.50000D-04
        HBS( 7, 1, 4) = -0.78200D-02
        HBS( 7, 2, 4) =  0.29170D-01
        HBS( 7, 3, 4) =  0.57100D-02
        HBS( 7, 4, 4) = -0.32990D-01
        HBS( 7, 5, 4) =  0.46200D-02
        HBS( 7, 6, 4) = -0.35760D-01
        HBS( 7, 7, 4) = -0.84700D-02
        HBS( 7, 8, 4) = -0.55800D-02
        HBS( 7, 9, 4) =  0.48900D-02
        HBS( 7,10, 4) =  0.43210D-01
        HBS( 7,11, 4) = -0.75000D-03
        HBS( 7,12, 4) = -0.32000D-03
        HBS( 8, 1, 4) =  0.31800D-02
        HBS( 8, 2, 4) = -0.15950D-01
        HBS( 8, 3, 4) = -0.12300D-02
        HBS( 8, 4, 4) =  0.30010D-01
        HBS( 8, 5, 4) = -0.27700D-02
        HBS( 8, 6, 4) =  0.41300D-02
        HBS( 8, 7, 4) =  0.40300D-02
        HBS( 8, 8, 4) =  0.21600D-02
        HBS( 8, 9, 4) = -0.19900D-02
        HBS( 8,10, 4) = -0.19160D-01
        HBS( 8,11, 4) = -0.23000D-03
        HBS( 8,12, 4) =  0.17000D-03
        HBS( 9, 1, 4) =  0.46400D-02
        HBS( 9, 2, 4) = -0.13220D-01
        HBS( 9, 3, 4) = -0.44800D-02
        HBS( 9, 4, 4) =  0.29800D-02
        HBS( 9, 5, 4) = -0.18500D-02
        HBS( 9, 6, 4) =  0.31630D-01
        HBS( 9, 7, 4) =  0.44300D-02
        HBS( 9, 8, 4) =  0.34300D-02
        HBS( 9, 9, 4) = -0.29100D-02
        HBS( 9,10, 4) = -0.24050D-01
        HBS( 9,11, 4) =  0.98000D-03
        HBS( 9,12, 4) =  0.15000D-03
        HBS(10, 1, 4) =  0.31600D-01
        HBS(10, 2, 4) =  0.31600D-01
        HBS(10, 3, 4) =  0.63400D-02
        HBS(10, 4, 4) = -0.13710D-01
        HBS(10, 5, 4) = -0.34000D-03
        HBS(10, 6, 4) = -0.29610D-01
        HBS(10, 7, 4) =  0.47800D-02
        HBS(10, 8, 4) = -0.29610D-01
        HBS(10, 9, 4) = -0.84900D-02
        HBS(10,10, 4) = -0.13710D-01
        HBS(10,11, 4) =  0.47800D-02
        HBS(10,12, 4) = -0.34000D-03
        HBS(11, 1, 4) = -0.17850D-01
        HBS(11, 2, 4) = -0.13750D-01
        HBS(11, 3, 4) = -0.31700D-02
        HBS(11, 4, 4) =  0.80700D-02
        HBS(11, 5, 4) =  0.73000D-03
        HBS(11, 6, 4) =  0.27610D-01
        HBS(11, 7, 4) = -0.29000D-02
        HBS(11, 8, 4) =  0.20000D-02
        HBS(11, 9, 4) =  0.42400D-02
        HBS(11,10, 4) =  0.56400D-02
        HBS(11,11, 4) = -0.18700D-02
        HBS(11,12, 4) = -0.38000D-03
        HBS(12, 1, 4) = -0.13750D-01
        HBS(12, 2, 4) = -0.17850D-01
        HBS(12, 3, 4) = -0.31700D-02
        HBS(12, 4, 4) =  0.56400D-02
        HBS(12, 5, 4) = -0.38000D-03
        HBS(12, 6, 4) =  0.20000D-02
        HBS(12, 7, 4) = -0.18700D-02
        HBS(12, 8, 4) =  0.27610D-01
        HBS(12, 9, 4) =  0.42400D-02
        HBS(12,10, 4) =  0.80700D-02
        HBS(12,11, 4) = -0.29000D-02
        HBS(12,12, 4) =  0.73000D-03
        HBS(13, 1, 4) =  0.29170D-01
        HBS(13, 2, 4) = -0.78200D-02
        HBS(13, 3, 4) =  0.57100D-02
        HBS(13, 4, 4) =  0.43210D-01
        HBS(13, 5, 4) = -0.32000D-03
        HBS(13, 6, 4) = -0.55800D-02
        HBS(13, 7, 4) = -0.75000D-03
        HBS(13, 8, 4) = -0.35760D-01
        HBS(13, 9, 4) =  0.48900D-02
        HBS(13,10, 4) = -0.32990D-01
        HBS(13,11, 4) = -0.84700D-02
        HBS(13,12, 4) =  0.46200D-02
        HBS(14, 1, 4) = -0.13220D-01
        HBS(14, 2, 4) =  0.46400D-02
        HBS(14, 3, 4) = -0.44800D-02
        HBS(14, 4, 4) = -0.24050D-01
        HBS(14, 5, 4) =  0.15000D-03
        HBS(14, 6, 4) =  0.34300D-02
        HBS(14, 7, 4) =  0.98000D-03
        HBS(14, 8, 4) =  0.31630D-01
        HBS(14, 9, 4) = -0.29100D-02
        HBS(14,10, 4) =  0.29800D-02
        HBS(14,11, 4) =  0.44300D-02
        HBS(14,12, 4) = -0.18500D-02
        HBS(15, 1, 4) = -0.15950D-01
        HBS(15, 2, 4) =  0.31800D-02
        HBS(15, 3, 4) = -0.12300D-02
        HBS(15, 4, 4) = -0.19160D-01
        HBS(15, 5, 4) =  0.17000D-03
        HBS(15, 6, 4) =  0.21600D-02
        HBS(15, 7, 4) = -0.23000D-03
        HBS(15, 8, 4) =  0.41300D-02
        HBS(15, 9, 4) = -0.19900D-02
        HBS(15,10, 4) =  0.30010D-01
        HBS(15,11, 4) =  0.40300D-02
        HBS(15,12, 4) = -0.27700D-02
        HBS(16, 1, 4) = -0.60700D-02
        HBS(16, 2, 4) = -0.24180D-01
        HBS(16, 3, 4) =  0.14620D-01
        HBS(16, 4, 4) =  0.48130D-01
        HBS(16, 5, 4) = -0.12600D-02
        HBS(16, 6, 4) =  0.37380D-01
        HBS(16, 7, 4) =  0.25000D-03
        HBS(16, 8, 4) = -0.10830D-01
        HBS(16, 9, 4) = -0.65000D-03
        HBS(16,10, 4) = -0.35760D-01
        HBS(16,11, 4) =  0.47400D-02
        HBS(16,12, 4) = -0.81300D-02
        HBS(17, 1, 4) =  0.30300D-02
        HBS(17, 2, 4) =  0.21630D-01
        HBS(17, 3, 4) = -0.86900D-02
        HBS(17, 4, 4) = -0.25720D-01
        HBS(17, 5, 4) =  0.12100D-02
        HBS(17, 6, 4) = -0.16220D-01
        HBS(17, 7, 4) = -0.13000D-03
        HBS(17, 8, 4) =  0.52000D-02
        HBS(17, 9, 4) = -0.27000D-03
        HBS(17,10, 4) =  0.48300D-02
        HBS(17,11, 4) = -0.17000D-02
        HBS(17,12, 4) =  0.60700D-02
        HBS(18, 1, 4) =  0.30300D-02
        HBS(18, 2, 4) =  0.25500D-02
        HBS(18, 3, 4) = -0.59400D-02
        HBS(18, 4, 4) = -0.22410D-01
        HBS(18, 5, 4) =  0.50000D-04
        HBS(18, 6, 4) = -0.21160D-01
        HBS(18, 7, 4) = -0.12000D-03
        HBS(18, 8, 4) =  0.56300D-02
        HBS(18, 9, 4) =  0.91000D-03
        HBS(18,10, 4) =  0.30930D-01
        HBS(18,11, 4) = -0.30400D-02
        HBS(18,12, 4) =  0.20600D-02
      endif

      end
      subroutine getgeom(nap,qse,nqs,qbe,nqb,qte,nqt,ii)

      double precision zero, pi
      parameter(zero=0.0d0,pi=3.1415926535898d0)

C Global variables
      integer ii
      integer nap
      integer nqs,nqb,nqt
      double precision qse(nqs,nap),qbe(nqb,nap),qte(nqt,nap)

C Local variables
      integer i,j

      qte(1,1)=zero
      qte(2,1)=pi
      qte(3,1)=pi
      qte(4,1)=zero
      qte(5,1)=zero
      qte(6,1)=pi
      qte(7,1)=pi
      qte(8,1)=zero
      qte(9,1)=zero
      qte(10,1)=pi
      qte(11,1)=pi
      qte(12,1)=zero
      qte(13,1)=zero
      qte(14,1)=pi
      qte(15,1)=pi
      qte(16,1)=zero
      qte(17,1)=zero
      qte(18,1)=pi
      qte(19,1)=pi
      qte(20,1)=zero
      qte(21,1)=zero
      qte(22,1)=pi
      qte(23,1)=pi
      qte(24,1)=zero

      do i=2,nap
        do j=1,nqt
          qte(j,i)=qte(j,1)
        enddo
      enddo

C Now get equilibrium geometrical parameters for qae
      if (ii .eq. 11) then
C equilibrium geometry at anchor point 1
        qse( 1, 1) =  0.263789D+01
        qse( 2, 1) =  0.263870D+01
        qse( 3, 1) =  0.258038D+01
        qse( 4, 1) =  0.263573D+01
        qse( 5, 1) =  0.206310D+01
        qse( 6, 1) =  0.263336D+01
        qse( 7, 1) =  0.205963D+01
        qse( 8, 1) =  0.263920D+01
        qse( 9, 1) =  0.205742D+01
        qse(10, 1) =  0.262993D+01
        qse(11, 1) =  0.205975D+01
        qse(12, 1) =  0.205790D+01
C equilibrium geometry at anchor point 1
        qbe( 1, 1) =  0.210140D+01
        qbe( 2, 1) =  0.213715D+01
        qbe( 3, 1) =  0.204464D+01
        qbe( 4, 1) =  0.208658D+01
        qbe( 5, 1) =  0.209471D+01
        qbe( 6, 1) =  0.210190D+01
        qbe( 7, 1) =  0.210459D+01
        qbe( 8, 1) =  0.208211D+01
        qbe( 9, 1) =  0.209649D+01
        qbe(10, 1) =  0.208155D+01
        qbe(11, 1) =  0.210036D+01
        qbe(12, 1) =  0.210127D+01
        qbe(13, 1) =  0.210834D+01
        qbe(14, 1) =  0.209320D+01
        qbe(15, 1) =  0.208164D+01
        qbe(16, 1) =  0.208392D+01
        qbe(17, 1) =  0.207495D+01
        qbe(18, 1) =  0.212431D+01
C equilibrium geometry at anchor point 2
        qse( 1, 2) =  0.263724D+01
        qse( 2, 2) =  0.264694D+01
        qse( 3, 2) =  0.257138D+01
        qse( 4, 2) =  0.263516D+01
        qse( 5, 2) =  0.206589D+01
        qse( 6, 2) =  0.263264D+01
        qse( 7, 2) =  0.205963D+01
        qse( 8, 2) =  0.263972D+01
        qse( 9, 2) =  0.205726D+01
        qse(10, 2) =  0.262912D+01
        qse(11, 2) =  0.205993D+01
        qse(12, 2) =  0.205801D+01
C equilibrium geometry at anchor point 2
        qbe( 1, 2) =  0.209269D+01
        qbe( 2, 2) =  0.215700D+01
        qbe( 3, 2) =  0.203350D+01
        qbe( 4, 2) =  0.209468D+01
        qbe( 5, 2) =  0.208443D+01
        qbe( 6, 2) =  0.210408D+01
        qbe( 7, 2) =  0.210277D+01
        qbe( 8, 2) =  0.208276D+01
        qbe( 9, 2) =  0.209766D+01
        qbe(10, 2) =  0.207979D+01
        qbe(11, 2) =  0.210131D+01
        qbe(12, 2) =  0.210209D+01
        qbe(13, 2) =  0.211103D+01
        qbe(14, 2) =  0.209320D+01
        qbe(15, 2) =  0.207895D+01
        qbe(16, 2) =  0.208541D+01
        qbe(17, 2) =  0.207654D+01
        qbe(18, 2) =  0.212123D+01
C equilibrium geometry at anchor point 3
        qse( 1, 3) =  0.264670D+01
        qse( 2, 3) =  0.266084D+01
        qse( 3, 3) =  0.252665D+01
        qse( 4, 3) =  0.263224D+01
        qse( 5, 3) =  0.206010D+01
        qse( 6, 3) =  0.263214D+01
        qse( 7, 3) =  0.205982D+01
        qse( 8, 3) =  0.263685D+01
        qse( 9, 3) =  0.205662D+01
        qse(10, 3) =  0.262787D+01
        qse(11, 3) =  0.206007D+01
        qse(12, 3) =  0.205634D+01
C equilibrium geometry at anchor point 3
        qbe( 1, 3) =  0.209445D+01
        qbe( 2, 3) =  0.215439D+01
        qbe( 3, 3) =  0.203434D+01
        qbe( 4, 3) =  0.208659D+01
        qbe( 5, 3) =  0.207835D+01
        qbe( 6, 3) =  0.211824D+01
        qbe( 7, 3) =  0.211238D+01
        qbe( 8, 3) =  0.207426D+01
        qbe( 9, 3) =  0.209655D+01
        qbe(10, 3) =  0.207613D+01
        qbe(11, 3) =  0.210353D+01
        qbe(12, 3) =  0.210352D+01
        qbe(13, 3) =  0.211662D+01
        qbe(14, 3) =  0.209569D+01
        qbe(15, 3) =  0.207088D+01
        qbe(16, 3) =  0.208020D+01
        qbe(17, 3) =  0.208110D+01
        qbe(18, 3) =  0.212189D+01
C equilibrium geometry at anchor point 4
        qse( 1, 4) =  0.265159D+01
        qse( 2, 4) =  0.265159D+01
        qse( 3, 4) =  0.251838D+01
        qse( 4, 4) =  0.263062D+01
        qse( 5, 4) =  0.205598D+01
        qse( 6, 4) =  0.263458D+01
        qse( 7, 4) =  0.205967D+01
        qse( 8, 4) =  0.263458D+01
        qse( 9, 4) =  0.205648D+01
        qse(10, 4) =  0.263062D+01
        qse(11, 4) =  0.205967D+01
        qse(12, 4) =  0.205598D+01
C equilibrium geometry at anchor point 4
        qbe( 1, 4) =  0.211210D+01
        qbe( 2, 4) =  0.208554D+01
        qbe( 3, 4) =  0.208554D+01
        qbe( 4, 4) =  0.207166D+01
        qbe( 5, 4) =  0.208711D+01
        qbe( 6, 4) =  0.212441D+01
        qbe( 7, 4) =  0.211593D+01
        qbe( 8, 4) =  0.207134D+01
        qbe( 9, 4) =  0.209592D+01
        qbe(10, 4) =  0.207909D+01
        qbe(11, 4) =  0.210205D+01
        qbe(12, 4) =  0.210205D+01
        qbe(13, 4) =  0.211593D+01
        qbe(14, 4) =  0.209592D+01
        qbe(15, 4) =  0.207134D+01
        qbe(16, 4) =  0.207166D+01
        qbe(17, 4) =  0.208711D+01
        qbe(18, 4) =  0.212441D+01
      else if (ii .eq. 22) then
C equilibrium geometry at anchor point 1
        qse( 1, 1) =  0.269941D+01
        qse( 2, 1) =  0.268363D+01
        qse( 3, 1) =  0.252916D+01
        qse( 4, 1) =  0.268093D+01
        qse( 5, 1) =  0.205887D+01
        qse( 6, 1) =  0.268592D+01
        qse( 7, 1) =  0.205454D+01
        qse( 8, 1) =  0.268227D+01
        qse( 9, 1) =  0.206068D+01
        qse(10, 1) =  0.268469D+01
        qse(11, 1) =  0.205375D+01
        qse(12, 1) =  0.205424D+01
C equilibrium geometry at anchor point 1
        qbe( 1, 1) =  0.216375D+01
        qbe( 2, 1) =  0.209466D+01
        qbe( 3, 1) =  0.202477D+01
        qbe( 4, 1) =  0.204773D+01
        qbe( 5, 1) =  0.209602D+01
        qbe( 6, 1) =  0.213944D+01
        qbe( 7, 1) =  0.208156D+01
        qbe( 8, 1) =  0.210850D+01
        qbe( 9, 1) =  0.209313D+01
        qbe(10, 1) =  0.214184D+01
        qbe(11, 1) =  0.206797D+01
        qbe(12, 1) =  0.207337D+01
        qbe(13, 1) =  0.207509D+01
        qbe(14, 1) =  0.209807D+01
        qbe(15, 1) =  0.211003D+01
        qbe(16, 1) =  0.205641D+01
        qbe(17, 1) =  0.207231D+01
        qbe(18, 1) =  0.215446D+01
C equilibrium geometry at anchor point 2
        qse( 1, 2) =  0.269988D+01
        qse( 2, 2) =  0.269065D+01
        qse( 3, 2) =  0.251103D+01
        qse( 4, 2) =  0.267864D+01
        qse( 5, 2) =  0.206339D+01
        qse( 6, 2) =  0.268396D+01
        qse( 7, 2) =  0.205439D+01
        qse( 8, 2) =  0.268024D+01
        qse( 9, 2) =  0.206089D+01
        qse(10, 2) =  0.268483D+01
        qse(11, 2) =  0.205390D+01
        qse(12, 2) =  0.205412D+01
C equilibrium geometry at anchor point 2
        qbe( 1, 2) =  0.215707D+01
        qbe( 2, 2) =  0.211094D+01
        qbe( 3, 2) =  0.201518D+01
        qbe( 4, 2) =  0.205661D+01
        qbe( 5, 2) =  0.208168D+01
        qbe( 6, 2) =  0.214490D+01
        qbe( 7, 2) =  0.207623D+01
        qbe( 8, 2) =  0.211213D+01
        qbe( 9, 2) =  0.209483D+01
        qbe(10, 2) =  0.214363D+01
        qbe(11, 2) =  0.206707D+01
        qbe(12, 2) =  0.207249D+01
        qbe(13, 2) =  0.207837D+01
        qbe(14, 2) =  0.209780D+01
        qbe(15, 2) =  0.210701D+01
        qbe(16, 2) =  0.205447D+01
        qbe(17, 2) =  0.207468D+01
        qbe(18, 2) =  0.215403D+01
C equilibrium geometry at anchor point 3
        qse( 1, 3) =  0.269988D+01
        qse( 2, 3) =  0.269065D+01
        qse( 3, 3) =  0.251103D+01
        qse( 4, 3) =  0.267864D+01
        qse( 5, 3) =  0.206339D+01
        qse( 6, 3) =  0.268396D+01
        qse( 7, 3) =  0.205439D+01
        qse( 8, 3) =  0.268024D+01
        qse( 9, 3) =  0.206089D+01
        qse(10, 3) =  0.268483D+01
        qse(11, 3) =  0.205390D+01
        qse(12, 3) =  0.205412D+01
C equilibrium geometry at anchor point 3
        qbe( 1, 3) =  0.215707D+01
        qbe( 2, 3) =  0.211094D+01
        qbe( 3, 3) =  0.201518D+01
        qbe( 4, 3) =  0.205661D+01
        qbe( 5, 3) =  0.208168D+01
        qbe( 6, 3) =  0.214490D+01
        qbe( 7, 3) =  0.207623D+01
        qbe( 8, 3) =  0.211213D+01
        qbe( 9, 3) =  0.209483D+01
        qbe(10, 3) =  0.214363D+01
        qbe(11, 3) =  0.206707D+01
        qbe(12, 3) =  0.207249D+01
        qbe(13, 3) =  0.207837D+01
        qbe(14, 3) =  0.209780D+01
        qbe(15, 3) =  0.210701D+01
        qbe(16, 3) =  0.205447D+01
        qbe(17, 3) =  0.207468D+01
        qbe(18, 3) =  0.215403D+01
C equilibrium geometry at anchor point 4
        qse( 1, 4) =  0.269988D+01
        qse( 2, 4) =  0.269065D+01
        qse( 3, 4) =  0.251103D+01
        qse( 4, 4) =  0.267864D+01
        qse( 5, 4) =  0.206339D+01
        qse( 6, 4) =  0.268396D+01
        qse( 7, 4) =  0.205439D+01
        qse( 8, 4) =  0.268024D+01
        qse( 9, 4) =  0.206089D+01
        qse(10, 4) =  0.268483D+01
        qse(11, 4) =  0.205390D+01
        qse(12, 4) =  0.205412D+01
C equilibrium geometry at anchor point 4
        qbe( 1, 4) =  0.215707D+01
        qbe( 2, 4) =  0.211094D+01
        qbe( 3, 4) =  0.201518D+01
        qbe( 4, 4) =  0.205661D+01
        qbe( 5, 4) =  0.208168D+01
        qbe( 6, 4) =  0.214490D+01
        qbe( 7, 4) =  0.207623D+01
        qbe( 8, 4) =  0.211213D+01
        qbe( 9, 4) =  0.209483D+01
        qbe(10, 4) =  0.214363D+01
        qbe(11, 4) =  0.206707D+01
        qbe(12, 4) =  0.207249D+01
        qbe(13, 4) =  0.207837D+01
        qbe(14, 4) =  0.209780D+01
        qbe(15, 4) =  0.210701D+01
        qbe(16, 4) =  0.205447D+01
        qbe(17, 4) =  0.207468D+01
        qbe(18, 4) =  0.215403D+01
      else if(ii .eq. 33) then
C equilibrium geometry at anchor point 1
        qse( 1, 1) =  0.272760D+01
        qse( 2, 1) =  0.272219D+01
        qse( 3, 1) =  0.243093D+01
        qse( 4, 1) =  0.259142D+01
        qse( 5, 1) =  0.205558D+01
        qse( 6, 1) =  0.266713D+01
        qse( 7, 1) =  0.205707D+01
        qse( 8, 1) =  0.267731D+01
        qse( 9, 1) =  0.205961D+01
        qse(10, 1) =  0.258566D+01
        qse(11, 1) =  0.205666D+01
        qse(12, 1) =  0.205567D+01
C equilibrium geometry at anchor point 1
        qbe( 1, 1) =  0.209913D+01
        qbe( 2, 1) =  0.214547D+01
        qbe( 3, 1) =  0.203859D+01
        qbe( 4, 1) =  0.207342D+01
        qbe( 5, 1) =  0.205792D+01
        qbe( 6, 1) =  0.215185D+01
        qbe( 7, 1) =  0.210210D+01
        qbe( 8, 1) =  0.209263D+01
        qbe( 9, 1) =  0.208845D+01
        qbe(10, 1) =  0.211486D+01
        qbe(11, 1) =  0.208558D+01
        qbe(12, 1) =  0.208274D+01
        qbe(13, 1) =  0.209348D+01
        qbe(14, 1) =  0.208777D+01
        qbe(15, 1) =  0.210193D+01
        qbe(16, 1) =  0.208337D+01
        qbe(17, 1) =  0.205127D+01
        qbe(18, 1) =  0.214854D+01
C equilibrium geometry at anchor point 2
        qse( 1, 2) =  0.274917D+01
        qse( 2, 2) =  0.274193D+01
        qse( 3, 2) =  0.237981D+01
        qse( 4, 2) =  0.259572D+01
        qse( 5, 2) =  0.205541D+01
        qse( 6, 2) =  0.266481D+01
        qse( 7, 2) =  0.205872D+01
        qse( 8, 2) =  0.266882D+01
        qse( 9, 2) =  0.205898D+01
        qse(10, 2) =  0.259216D+01
        qse(11, 2) =  0.205873D+01
        qse(12, 2) =  0.205800D+01
C equilibrium geometry at anchor point 2
        qbe( 1, 2) =  0.206394D+01
        qbe( 2, 2) =  0.213351D+01
        qbe( 3, 2) =  0.208574D+01
        qbe( 4, 2) =  0.208869D+01
        qbe( 5, 2) =  0.205086D+01
        qbe( 6, 2) =  0.214363D+01
        qbe( 7, 2) =  0.210572D+01
        qbe( 8, 2) =  0.209443D+01
        qbe( 9, 2) =  0.208304D+01
        qbe(10, 2) =  0.211155D+01
        qbe(11, 2) =  0.208615D+01
        qbe(12, 2) =  0.208549D+01
        qbe(13, 2) =  0.209340D+01
        qbe(14, 2) =  0.208795D+01
        qbe(15, 2) =  0.210183D+01
        qbe(16, 2) =  0.210307D+01
        qbe(17, 2) =  0.204251D+01
        qbe(18, 2) =  0.213761D+01
C equilibrium geometry at anchor point 3
        qse( 1, 3) =  0.274751D+01
        qse( 2, 3) =  0.274706D+01
        qse( 3, 3) =  0.235956D+01
        qse( 4, 3) =  0.259829D+01
        qse( 5, 3) =  0.205917D+01
        qse( 6, 3) =  0.266445D+01
        qse( 7, 3) =  0.205922D+01
        qse( 8, 3) =  0.266458D+01
        qse( 9, 3) =  0.205897D+01
        qse(10, 3) =  0.259790D+01
        qse(11, 3) =  0.205917D+01
        qse(12, 3) =  0.205900D+01
C equilibrium geometry at anchor point 3
        qbe( 1, 3) =  0.204756D+01
        qbe( 2, 3) =  0.212154D+01
        qbe( 3, 3) =  0.211409D+01
        qbe( 4, 3) =  0.210600D+01
        qbe( 5, 3) =  0.204389D+01
        qbe( 6, 3) =  0.213330D+01
        qbe( 7, 3) =  0.209805D+01
        qbe( 8, 3) =  0.210013D+01
        qbe( 9, 3) =  0.208501D+01
        qbe(10, 3) =  0.211048D+01
        qbe(11, 3) =  0.208627D+01
        qbe(12, 3) =  0.208644D+01
        qbe(13, 3) =  0.209637D+01
        qbe(14, 3) =  0.208569D+01
        qbe(15, 3) =  0.210112D+01
        qbe(16, 3) =  0.210790D+01
        qbe(17, 3) =  0.204012D+01
        qbe(18, 3) =  0.213517D+01
C equilibrium geometry at anchor point 4
        qse( 1, 4) =  0.274535D+01
        qse( 2, 4) =  0.274535D+01
        qse( 3, 4) =  0.236001D+01
        qse( 4, 4) =  0.259897D+01
        qse( 5, 4) =  0.205919D+01
        qse( 6, 4) =  0.266411D+01
        qse( 7, 4) =  0.205914D+01
        qse( 8, 4) =  0.266411D+01
        qse( 9, 4) =  0.205909D+01
        qse(10, 4) =  0.259897D+01
        qse(11, 4) =  0.205914D+01
        qse(12, 4) =  0.205919D+01
C equilibrium geometry at anchor point 4
        qbe( 1, 4) =  0.204757D+01
        qbe( 2, 4) =  0.211781D+01
        qbe( 3, 4) =  0.211781D+01
        qbe( 4, 4) =  0.210758D+01
        qbe( 5, 4) =  0.204104D+01
        qbe( 6, 4) =  0.213456D+01
        qbe( 7, 4) =  0.209631D+01
        qbe( 8, 4) =  0.210152D+01
        qbe( 9, 4) =  0.208536D+01
        qbe(10, 4) =  0.211101D+01
        qbe(11, 4) =  0.208609D+01
        qbe(12, 4) =  0.208609D+01
        qbe(13, 4) =  0.209631D+01
        qbe(14, 4) =  0.208536D+01
        qbe(15, 4) =  0.210152D+01
        qbe(16, 4) =  0.210758D+01
        qbe(17, 4) =  0.204104D+01
        qbe(18, 4) =  0.213456D+01
      endif

      end
      subroutine tent(nap,rap,r,ir,t,dtdr)
************************************************************************
* Subroutine to evaluate tent functions
************************************************************************

      implicit double precision (a-h,o-z)

      integer ir
      integer nap
      double precision r
      double precision rap(*),t(*),dtdr(*)

C Local variables
      integer i
      double precision u,v

      do i=1,nap
        t(i)=0.0d0
        dtdr(i)=0.0d0
      enddo

C Check whether the point are end points
      if (r .lt. rap(1)) then
        ir=1
        t(1)=1.0d0
        dtdr(1)=0.0d0
      else if (r .gt. rap(nap)) then
        ir=nap+1
        t(nap)=1.0d0
        dtdr(nap)=0.0d0
      endif

C Determine i for which r(i-1) <= r < r(i)
      do i=2,nap
        if ((rap(i-1) .le. r) .and. (r .lt. rap(i))) then
          ir=i
          u=(r-rap(i-1))*(r-rap(i-1))*(r-rap(i-1))*(r-rap(i-1))
          v=(r-rap(i))*(r-rap(i))*(r-rap(i))*(r-rap(i))
          dudr=4*(r-rap(i-1))*(r-rap(i-1))*(r-rap(i-1))
          dvdr=4*(r-rap(i))*(r-rap(i))*(r-rap(i))
          t(i-1)=v/(u+v)
          t(i)=u/(u+v)
          dtdr(i-1)=(u*dvdr-v*dudr)/((u+v)*(u+v))
          dtdr(i)=(v*dudr-u*dvdr)/((u+v)*(u+v))
          exit
        endif
      enddo 
      
      end 
      subroutine evu1r(igrad,r,v,dvdr)
***********************************************************************
* Subroutine to evaluate Varshini potential and derivative for given r
***********************************************************************

      implicit double precision (a-h,o-z)

      double precision de,re,beta
      parameter(de=0.17596055d0,re=1.93206106d0,
     $          beta=0.63822703d0)

      integer igrad
      double precision r,v,dvdr
      double precision x,y,dvdy,dydx,dxdr

      x=r/re
      y=1.0d0-exp(-beta*(x*x-1.0d0))/x
      v=de*y*y
      
      if (igrad .eq. 1) then
        dvdy=2.0d0*de*y
        dydx=(2.0d0*beta*x*x+1.0d0)*exp(-beta*(x*x-1.0d0))/(x*x)
        dxdr=1.0d0/re
        dvdr=dvdy*dydx*dxdr
      endif

      end
      subroutine evu1b(igrad,r,theta,vbnd,dvdr,dvdth)
***********************************************************************
* Subroutine to evaluate potential energy from H-O-C bending 
* and derivatives for given rOH and aHOC
***********************************************************************

      implicit double precision (a-h,o-z)

      integer ngk2
      parameter (ngk2=1)

      double precision a2(ngk2)
      data a2/0.10449356d0/

      double precision alpha2(ngk2)
      data alpha2/0.27051185d0/

      double precision r2(ngk2)
      data r2/2.28431658d0/

      double precision cthta1,cthta2,a,r0
      data cthta1/-0.39964587d0/,cthta2/-0.24769126d0/
      data a/4.01097890d0/,r0/1.83449325d0/

      integer i,j,k
      double precision r
      double precision k2

      cthta=cos(theta)
      sthta=sin(theta)

      x=a*(r-r0)
C      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))      
      tx=(1.0d0-exp(-2.0d0*x))/(1.0d0+exp(-2.0d0*x))
      cthta0=cthta1+0.5d0*(1+tx)*(cthta2-cthta1)

      dxdr=a
      dtxdx=4.0d0/(2.0d0+exp(2.0d0*x)+exp(-2.0d0*x))
      dtxdr=dtxdx*dxdr
      dct0dr=0.5d0*(cthta2-cthta1)*dtxdr

      k2=0.0d0
      dk2dr=0.0d0
      do i=1,ngk2
        gr=a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
        k2=k2+gr
        if (igrad .eq. 1) then
          dk2dr=dk2dr-2.0d0*alpha2(i)*(r-r2(i))*gr
        endif
      enddo

      dthta=cthta-cthta0

      vbnd=k2*dthta*dthta

      if (igrad .eq. 1) then
        dvdr=dk2dr*dthta*dthta - 2.0d0*k2*dthta*dct0dr 
        dvdth=-2.0d0*k2*dthta*sthta
      endif

      end

      subroutine evu1to(igrad,r,phi,v,dvdr,dvdph)
***********************************************************************
* Subroutine to evaluate torsion potentials and derivatives 
* for given r and phi
***********************************************************************

      implicit double precision (a-h,o-z)

      double precision a1,alpha1,r1
      double precision a2,alpha2,r2
      parameter(a1=0.02371114d0,alpha1=0.96558083d0,r1=3.25505778d0)
      parameter(a2=0.00349299d0,alpha2=1.092945551d1,r2=2.10516007d0)

      integer igrad
      double precision phi,r,v,dvdph,dvdr
      double precision x,y,dvdy,dydx,dxdr

      x1=a1*exp(-alpha1*(r-r1)*(r-r1))
      x2=a2*exp(-alpha2*(r-r2)*(r-r2))
      y=1.0d0-cos(2.0d0*phi)
      v=x1*y+x2*y*y

      if (igrad .eq. 1) then
        dx1dr=-2.0d0*alpha1*(r-r1)*x1
        dx2dr=-2.0d0*alpha2*(r-r2)*x2
        dydph=2.0d0*sin(2.0d0*phi)
        dvdr=dx1dr*y+dx2dr*y*y
        dvdph=x1*dydph+2.0d0*x2*y*dydph
      endif

      end
      subroutine evu2r(igrad,r,v,dvdr)
***********************************************************************
* Subroutine to evaluate Morse potential and derivative for given r
***********************************************************************

      implicit double precision (a-h,o-z)

      double precision alpha,de,re,v0
      parameter(alpha=0.67269344d0,de=0.41866243d0,
     $          re=1.95347670d0,v0=0.17280158d0)

      integer igrad
      double precision r,v,dvdr
      double precision x,y,dvdy,dydx,dxdr

      x=r-re
      y=1.0d0-exp(-alpha*x)
      v=de*y*y+v0
      
      if (igrad .eq. 1) then
        dvdy=2.0d0*de*y
        dydx=alpha*exp(-alpha*x)
        dxdr=1.0d0
        dvdr=dvdy*dydx*dxdr
      endif

      end
      subroutine evu2b(igrad,r,theta,vbnd,dvdr,dvdth)
***********************************************************************
* Subroutine to evaluate potential energy from H-O-C bending 
* and derivatives for given rOH and aHOC
***********************************************************************

      implicit double precision (a-h,o-z)

      integer ngk2
      parameter (ngk2=1)

      double precision a2(ngk2)
      data a2/0.10090668d0/

      double precision alpha2(ngk2)
      data alpha2/0.21154184d0/

      double precision r2(ngk2)
      data r2/1.80972781d0/

      double precision cthta1,cthta2,a,r0
      data cthta1/-0.31835353d0/,cthta2/-0.22230161d0/
      data a/4.73544047d0/,r0/2.01933829d0/

      integer i,j,k
      double precision r
      double precision k2

      cthta=cos(theta)
      sthta=sin(theta)

      x=a*(r-r0)
C      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))      
      tx=(1.0d0-exp(-2.0d0*x))/(1.0d0+exp(-2.0d0*x))
      cthta0=cthta1+0.5d0*(1+tx)*(cthta2-cthta1)

      dxdr=a
      dtxdx=4.0d0/(2.0d0+exp(2.0d0*x)+exp(-2.0d0*x))
      dtxdr=dtxdx*dxdr
      dct0dr=0.5d0*(cthta2-cthta1)*dtxdr

      k2=0.0d0
      dk2dr=0.0d0
      do i=1,ngk2
        gr=a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
        k2=k2+gr
        if (igrad .eq. 1) then
          dk2dr=dk2dr-2.0d0*alpha2(i)*(r-r2(i))*gr
        endif
      enddo

      dthta=cthta-cthta0

      vbnd=k2*dthta*dthta

      if (igrad .eq. 1) then
        dvdr=dk2dr*dthta*dthta - 2.0d0*k2*dthta*dct0dr 
        dvdth=-2.0d0*k2*dthta*sthta
      endif

      end

      subroutine evu2to(igrad,r,phi,v,dvdr,dvdph)
***********************************************************************
* Subroutine to evaluate torsion potentials and derivatives 
* for given r and phi
***********************************************************************

      implicit double precision (a-h,o-z)

      double precision a1,alpha1,r1
      parameter(a1=0.01925782d0,alpha1=5.00873578d0,r1=2.22678700d0)

      integer igrad
      double precision phi,r,v,dvdph,dvdr
      double precision x,y,dvdy,dydx,dxdr

      x=a1*exp(-alpha1*(r-r1)*(r-r1))
      y=1.0d0-cos(2.0d0*phi)
      v=x*y
      
      if (igrad .eq. 1) then
        dxdr=-2.0d0*alpha1*(r-r1)*x
        dydph=2.0d0*sin(2.0d0*phi)
        dvdr=y*dxdr
        dvdph=x*dydph 
      endif

      end
      subroutine evu3r(igrad,r,v,dvdr)
***********************************************************************
* Subroutine to evaluate Morse potential and derivative for given r
***********************************************************************

      implicit double precision (a-h,o-z)

      double precision v0
      double precision a1,alpha1,r1
      double precision a2,alpha2,r2
      double precision a3,alpha3,r3
      parameter(v0=0.15976218d0)
      parameter(a1=0.06496177d0,alpha1=1.81347726d0,r1=2.63168269d0)
      parameter(a2=-0.00077880d0,alpha2=3.24769899d0,r2=3.71412291d0)
      parameter(a3=9.3151142407d2,alpha3=4.73016534d0,r3=-0.04117486d0)

      integer igrad
      double precision r,v,dvdr
      double precision x,y,dvdy,dydx,dxdr

      x1=r-r1
      y1=a1*exp(-alpha1*x1)
      x2=r-r2
      y2=a2*exp(-alpha2*x2)
      x3=r-r3
      y3=a3*exp(-alpha3*x3)
      v=v0+y1+y2+y3
      
      if (igrad .eq. 1) then
        dx1dr=1.0d0
        dx2dr=1.0d0
        dx3dr=1.0d0
        dy1dx1=-alpha1*y1
        dy2dx2=-alpha2*y2
        dy3dx3=-alpha3*y3
        dvdr=dy1dx1*dx1dr+dy2dx2*dx2dr+dy3dx3*dx3dr
      endif

      end
      subroutine evu3b(igrad,r,theta,vbnd,dvdr,dvdth)
***********************************************************************
* Subroutine to evaluate potential energy from H-O-C bending 
* and derivatives for given rOH and aHOC
***********************************************************************

      implicit double precision (a-h,o-z)

      integer ngk2,ngk3
      parameter (ngk2=1,ngk3=1)

      double precision a2(ngk2),a3(ngk3)
      data a2/3.00968695d0/
      data a3/0.06538039d0/

      double precision alpha2(ngk2),alpha3(ngk3)
      data alpha2/0.04954325d0/
      data alpha3/14.94793100d0/

      double precision r2(ngk2),r3(ngk3)
      data r2/-6.33258333d0/
      data r3/2.09315608d0/

      double precision cthta1,cthta2,a,r0
      data cthta1/-0.30362526d0/,cthta2/-0.65443171d0/
      data a/0.90134535d0/,r0/2.95896455d0/

      integer i,j,k
      double precision r
      double precision k2,k3

      cthta=cos(theta)
      sthta=sin(theta)

      x=a*(r-r0)
C      tx=(exp(x)-exp(-x))/(exp(x)+exp(-x))       
      tx=(1.0d0-exp(-2.0d0*x))/(1.0d0+exp(-2.0d0*x))
      cthta0=cthta1+0.5d0*(1+tx)*(cthta2-cthta1)

      dxdr=a
      dtxdx=4.0d0/(2.0d0+exp(2.0d0*x)+exp(-2.0d0*x))
      dtxdr=dtxdx*dxdr
      dct0dr=0.5d0*(cthta2-cthta1)*dtxdr

      k2=0.0d0
      dk2dr=0.0d0
      do i=1,ngk2
        gr=a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
        k2=k2+gr
        if (igrad .eq. 1) then
          dk2dr=dk2dr-2.0d0*alpha2(i)*(r-r2(i))*gr
        endif
      enddo
 
      k3=0.0d0
      dk3dr=0.0d0
      do i=1,ngk3
        gr=a3(i)*exp(-alpha3(i)*(r-r3(i))*(r-r3(i)))
        k3=k3+gr
        if (igrad .eq. 1) then
          dk3dr=dk3dr-2.0d0*alpha3(i)*(r-r3(i))*gr
        endif
      enddo

      dthta=cthta-cthta0
      vbnd=k2*dthta*dthta + k3*dthta*dthta*dthta
      if (igrad .eq. 1) then
        dvdr=dk2dr*dthta*dthta + dk3dr*dthta*dthta*dthta
     $      - 2.0d0*k2*dthta*dct0dr - 3.0d0*k3*dthta*dthta*dct0dr
        dvdth=-2.0d0*k2*dthta*sthta - 3.0d0*k3*dthta*dthta*sthta
      endif

      end

      subroutine evu3to(igrad,r,phi,v,dvdr,dvdph)
***********************************************************************
* Subroutine to evaluate torsion potentials and derivatives 
* for given r and phi
***********************************************************************

      implicit double precision (a-h,o-z)

      double precision a1,alpha1,r1
      double precision a2,alpha2,r2
      parameter(a1=-0.02298368d0,alpha1=1.06740430d0,r1=3.28664483d0)
      parameter(a2=0.01138152d0,alpha2=0.65339172d0,r2=2.58341105d0)

      integer igrad
      double precision phi,r,v,dvdph,dvdr
      double precision x,y,dvdy,dydx,dxdr

      x1=a1*exp(-alpha1*(r-r1)*(r-r1))
      x2=a2*exp(-alpha2*(r-r2)*(r-r2))
      y=1.0d0-cos(2.0d0*phi)
      v=x1*y+x2*y*y
      
      if (igrad .eq. 1) then
        dx1dr=-2.0d0*alpha1*(r-r1)*x1
        dx2dr=-2.0d0*alpha2*(r-r2)*x2
        dydph=2.0d0*sin(2.0d0*phi)
        dvdr=dx1dr*y+dx2dr*y*y
        dvdph=x1*dydph+2.0d0*x2*y*dydph 
      endif

      end
      subroutine evubm(igrad,natm,X,Y,Z,ubm,dubmdX)
***************************************************************************
* Subroutine to evulate the Born-Mayer repulsions between non-bonded
* C atoms
***************************************************************************

      implicit none

      double precision B,alpha
      parameter (B=6.693126d1,alpha=1.894454d0)

      integer i,j
      integer igrad,natm
      double precision ubm
      double precision X(natm),Y(natm),Z(natm)
      double precision dubmdX(3*natm)
      double precision bndlen
      external bndlen
C Local variables
      integer nbl,nba,nto,noop
      data nbl/3/,nba/0/,nto/0/,noop/0/
      integer ibl(2,3),iba(3,0),ito(4,0),ioop(4,0)
      data ibl /1,4, 2,5, 3,6/
      double precision r14,r25,r36
      double precision dudr14,dudr25,dudr36
      double precision bm(3,3*natm)

      r14=bndlen(1,4,X,Y,Z)
      r25=bndlen(2,5,X,Y,Z)
      r36=bndlen(3,6,X,Y,Z)

      ubm=B*(exp(-alpha*r14)+exp(-alpha*r25)+exp(-alpha*r36))

C 
C      write(*,*) "ubm =",ubm
C
      if (igrad .eq. 1) then
        dudr14=-alpha*B*exp(-alpha*r14)
        dudr25=-alpha*B*exp(-alpha*r25)
        dudr36=-alpha*B*exp(-alpha*r36)
        
        call bmat3(nbl,ibl,nba,iba,nto,ito,noop,ioop,natm,X,Y,Z,bm,3)
        do i=1,3*natm
          dubmdX(i)=0.d0
        enddo

        do i=1,3*natm
          dubmdX(i)=dubmdX(i)+dudr14*bm(1,i)+dudr25*bm(2,i)
     $                       +dudr36*bm(3,i)
        enddo

C 
C        do i=1,3*natm
C          write(*,*) dubmdX(i)
C        enddo
C
      endif

      end
C      subroutine evuij(nap,natm,igrad,ir,t,dtdr,X,Y,Z,u13,u23,
C     $                 du13dX,du23dX)
C Add u12
      subroutine evuij(nap,natm,igrad,ir,t,dtdr,X,Y,Z,u12,u13,u23,
     $                 du12dX,du13dX,du23dX)
************************************************************************
* Subroutine to evaluate the diabatic couplings and gradients  
* The bohr-hartree-radian unit is used throughout the program
* qtc:  redundant internal coordinates
* stc:  non-redundant internal coordinates
* Ct:   transform matrix between RICs and NRICs
*       stc=Ct*qtc
* bmqt: B-matrix of RICs
* bmst: B-matrix of NRICs
************************************************************************

      implicit double precision (a-h,o-z)

      integer maxdim
      parameter (maxdim=100)
      
      integer ldbmt,ldbms,ldct
      parameter (ldbmt=maxdim,ldbms=maxdim,ldbnrt=maxdim,ldct=maxdim)

      integer nstc
      parameter (nstc=9)

      double precision pi
      parameter (pi=3.1415926535898d0)

C Global variables
      integer igrad,ir
      integer natm
C      double precision u13,u23
      double precision u12,u13,u23
      double precision t(*),dtdr(*)
      double precision X(*),Y(*),Z(*)
C      double precision du13dX(3*natm),du23dX(3*natm)
      double precision du12dX(3*natm),du13dX(3*natm),du23dX(3*natm)

C Local variables
      integer nqtc
      integer nblt,nbat,ntot,noopt
C     data(nblt=0,nbat=0,ntot=6,noopt=6)
      data nblt/0/,nbat/0/,ntot/6/,noopt/6/
      data nblt/0/,nbat/0/,ntot/6/,noopt/6/
      integer iblt(2,0),ibat(3,0),itot(4,6),ioopt(4,6)
      data itot  /2,1,6,5,
     $            1,6,5,4,
     $            6,5,4,3,
     $            5,4,3,2,
     $            4,3,2,1,
     $            3,2,1,6/
      data ioopt /7,1,2,6,
     $            12,6,1,5,
     $            11,5,6,4,
     $            10,4,5,3,
     $            9,3,4,2,
     $            8,2,3,1/

      integer nbls,nbas,ntos,noops
C      data(nbls=1,nbas=0,ntos=1,noops=0)
      data nbls/1/,nbas/0/,ntos/1/,noops/0/

      integer ibls(2,1),ibas(3,0),itos(4,1),ioops(4,0)
      data ibls /7,13/
      data itos /2,1,7,13/

      double precision half,oossix,oostlv
      parameter(half=0.5d0,oossix=0.408248290463863d0,
     $          oostlv=0.288675134594813d0)

      double precision qtc(maxdim),qsc(maxdim)
      double precision stc(maxdim)
C KRY 20140718
      double precision d12tdr,d12sdr
      double precision d13tdr,d13sdr,du13dr,du13dph
      double precision d23tdr,d23sdr,du23dr,du23dph
C      double precision du13t(maxdim),du23t(maxdim)
      double precision du12t(maxdim),du13t(maxdim),du23t(maxdim)
      double precision bmt(ldbmt,3*natm),bms(ldbms,3*natm)
      double precision bnrt(ldbnrt,3*natm)
      double precision ct(ldct,maxdim)

C Initialize Ct matrix
      do i=1,9
        do j=1,12
          Ct(i,j)=0.0d0
        enddo
      enddo
C Initialize Qtc vectors
      do i=1,nblt+nbat+ntot+noopt
        qtc(i)=0.0d0
      enddo
C Initialize Stc vectors
      do i=1,nstc
        stc(i)=0.0d0
      enddo

      Ct(1,1)=oossix
      Ct(1,2)=-oossix
      Ct(1,3)=oossix
      Ct(1,4)=-oossix
      Ct(1,5)=oossix
      Ct(1,6)=-oossix
      Ct(2,1)=-oostlv
      Ct(2,2)=2.0d0*oostlv
      Ct(2,3)=-oostlv
      Ct(2,4)=-oostlv
      Ct(2,5)=2.0d0*oostlv
      Ct(2,6)=-oostlv
      Ct(3,1)=-half
      Ct(3,3)=half
      Ct(3,4)=-half
      Ct(3,6)=half
      Ct(4,7)=1.0d0
      Ct(5,8)=1.0d0
      Ct(6,9)=1.0d0
      Ct(7,10)=1.0d0
      Ct(8,11)=1.0d0
      Ct(9,12)=1.0d0

C Evaluate internal coordinates
!Add by J. Zheng
C      qtc(:) = 0d0   
      call evintcoord(nblt,iblt,nbat,ibat,ntot,itot,noopt,ioopt,
     $                nqtc,qtc,X,Y,Z)
C Test
C      do i=1,nqtc
C        write(*,'(A,I2,A,F8.4)') 'qtc(',i,') = ',qtc(i)*180.d0/pi
C      enddo
C End
       

!Add by J. Zheng
C      stc(:) = 0d0
      call matvect(stc,nstc,Ct,ldct,qtc,nqtc)

      r=bndlen(7,13,X,Y,Z)
      phi=dihedr(2,1,7,13,X,Y,Z)

C Evaluate Wilson's B-matrix if igrad = 1
      if (igrad .eq. 1) then
!Added by J. Zheng
C        bmt(:,:) = 0d0
        call bmat3(nblt,iblt,nbat,ibat,ntot,itot,noopt,ioopt,natm,
     $            X,Y,Z,bmt,ldbmt)
C Evalate the Wilson's B-matrix in nonredundant internal coordinates
        call matmult(bnrt,ldbnrt,Ct,ldct,nstc,bmt,ldbmt,nqtc,3*natm)
!Added by J. Zheng
C        bms(:,:) = 0d0
        call bmat3(nbls,ibls,nbas,ibas,ntos,itos,noops,ioops,natm,
     $            X,Y,Z,bms,ldbms)
      endif

C Evaluate u12 KRY 20140718
      call evuijt(nap,igrad,ir,t,dtdr,nstc,stc,u12t,du12t,d12tdr,12)
C Evaluate u13 and u23 from tertiary out-of-plane coordinates    
      call evuijt(nap,igrad,ir,t,dtdr,nstc,stc,u13t,du13t,d13tdr,13)
      call evuijt(nap,igrad,ir,t,dtdr,nstc,stc,u23t,du23t,d23tdr,23)

C Evaluate u13 and u23 from secondary out-of-plane coordinates
      call evu13to(r,phi,u13s,d13dph,d13sdr)
      call evu23to(r,phi,u23s,d23dph,d23sdr)

      u12=u12t                  ! KRY 20140718
      u13=u13t+u13s
      u23=u23t+u23s

C Calculated analytical Cartesian gradient if igrad = 1
      if (igrad .eq. 1) then
        do i=1,3*natm
          du12dX(i)=0.0d0
          du13dX(i)=0.0d0
          du23dX(i)=0.0d0
        enddo
        du12dr=d12tdr
        du13dr=d13tdr+d13sdr
        du23dr=d23tdr+d23sdr
C Contribution of duijdX from duijdr
        do i=1,3*natm
          du12dX(i)=du12dX(i)+du12dr*bms(1,i)
          du13dX(i)=du13dX(i)+du13dr*bms(1,i)
          du23dX(i)=du23dX(i)+du23dr*bms(1,i)
        enddo
C Contribution of duijdX from dijdph        
        do i=1,3*natm
          du13dX(i)=du13dX(i)+d13dph*bms(2,i)
          du23dX(i)=du23dX(i)+d23dph*bms(2,i)
        enddo
C Contribution of duijdX from duijt
        do i=1,3*natm
          do j=1,nstc
            du12dX(i)=du12dX(i)+du12t(j)*bnrt(j,i)
            du13dX(i)=du13dX(i)+du13t(j)*bnrt(j,i)
            du23dX(i)=du23dX(i)+du23t(j)*bnrt(j,i)
          enddo
        enddo
      endif

      end
      subroutine evuijt(nap,igrad,ir,t,dtdr,ntc,qtc,uij,
     $                  duij,duijdr,ij)

      implicit double precision (a-h,o-z)

      integer ij 
      integer igrad,ir
      integer ntc
      double precision t(*),dtdr(*)
      double precision qtc(*)
      double precision uij
      double precision duijdr
      double precision duij(*)

C Local variables
      integer i,j
      double precision u(nap)
      double precision a(nap,ntc)
      double precision uij0(nap)
      double precision uij1,uij2
      double precision duij1(ntc),duij2(ntc)

C Initialize uij,uij1,uij2
      uij=0.0d0
      uij1=0.0d0
      uij2=0.0d0
C Initialize duij,duij1,duij2
      do i=1,ntc
        duij(i)=0.0d0
        duij1(i)=0.0d0
        duij2(i)=0.0d0
      enddo
C Initialize duijdr
      duijdr=0.0d0

C Initialize numerical gradients of diabatic couplings
      do i=1,nap
        uij0(i)=0.0d0
        do j=1,9
          a(i,j)=0.0d0
        enddo
      enddo
C Assign u120 for anchor point 1 and 2
      if (ij .eq. 12) then
        uij0(1)=-0.000616170d0
        uij0(2)=-0.001064286d0
      endif

C Assign parameters based ij 
      if (ij .eq. 12) then
        a(3,1)=0.0d0
        a(3,2)=0.0d0
        a(3,3)=0.0d0
        a(3,4)=0.0d0
        a(3,5)=0.0d0
        a(3,6)=0.0d0
        a(3,7)=0.0d0
        a(3,8)=0.0d0
        a(3,9)=0.0d0
      else if (ij .eq. 13) then
        a(3,1)=-0.00070d0
        a(3,2)=-0.00966d0
        a(3,3)=-0.00053d0
        a(3,4)=0.00193d0
        a(3,5)=0.00056d0
        a(3,6)=-0.00362d0
        a(3,7)=0.00030d0
        a(3,8)=0.00377d0
        a(3,9)=0.00101d0
      else if (ij .eq. 23) then
        a(2,1)=0.00452d0
        a(2,2)=-0.01212d0
        a(2,3)=-0.00088d0
        a(2,4)=0.00532d0
        a(2,5)=0.00757d0
        a(2,6)=-0.00348d0
        a(2,7)=0.00023d0
        a(2,8)=-0.00209d0
        a(2,9)=-0.01501d0
      else
        write(*,*) "Only parameters for U12, U13 and U23 are allowed!"
        stop
      endif

C Check for consistence
      if (ntc .ne. 9) then
        write(*,*) "Inconsistence in the number of out-of-plane tertiary
     $ coordinates"
        stop
      endif
      if (nap .ne. 4) then
        write(*,*) "Inconsistence in the number of anchor structures"
        stop
      endif

C Now calculate coupling contributions from tertiary coordinates
      if (ir .eq. 1) then
        uij=uij0(1)
        do i=1,ntc
          uij=uij+a(1,i)*qtc(i)
        enddo
      else if (ir .eq. nap+1) then
        uij=uij0(nap)
        do i=1,ntc
          uij=uij+a(nap,i)*qtc(i)
        enddo
      else if ((ir .gt. 1) .and. (ir .lt. nap+1)) then
        uij1=uij0(ir-1)
        uij2=uij0(ir)
        do i=1,ntc
          uij1=uij1+a(ir-1,i)*qtc(i)
          uij2=uij2+a(ir,i)*qtc(i)
        enddo
        uij=uij1*t(ir-1)+uij2*t(ir)
      endif

C Now calculate the derivatives of couplings if igrad =1
      if (igrad .eq. 1) then
        if (ir .eq. 1) then
          do i=1,ntc
            duij(i)=a(1,i)
          enddo
          duijdr=uij*dtdr(1)
        else if (ir .eq. nap+1) then
          do i=1,ntc
            duij(i)=a(nap,i)
          enddo
          duijdr=uij*dtdr(nap)
        else if ((ir .gt. 1) .and. (ir .lt. nap+1)) then
          do i=1,ntc
            duij1(i)=a(ir-1,i)
            duij2(i)=a(ir,i)
            duij(i)=duij1(i)*t(ir-1)+duij2(i)*t(ir)
          enddo
          duijdr=uij1*dtdr(ir-1)+uij2*dtdr(ir) 
        endif
      endif   

      end
      subroutine evu13to(r,phi,u13to,d13dph,d13dr)
***********************************************************************
* Little subroutine to evaluate U13 from C-C-O-H torsion
* U13(r,phi)=f1(r)*sin(phi)+f3(r)*sin(phi)**3+f5(r)*sin(phi)**5
* where fi(r) are expanded with a linear combination of Gaussians
***********************************************************************

      implicit double precision (a-h,o-z)

      integer ngf
      parameter (ngf=3)

      double precision a1(ngf),a3(ngf),a5(ngf)
      data a1/0.068885d0,0.034221d0,0.001544d0/
      data a3/-0.082676d0,-0.021770d0,0.0d0/
      data a5/0.053915d0,-0.008539d0,.0d0/          


      double precision alpha1(ngf),alpha3(ngf),alpha5(ngf)
      data alpha1/0.2956970d1,0.4695508d1,0.1167016d1/
      data alpha3/0.5029657d1,0.8744002d1,0.0d0/
      data alpha5/0.7263172d1,0.27942459d2,0.0d0/

      double precision r1(ngf),r3(ngf),r5(ngf)
      data r1/0.3162075d1,0.2255351d1,0.4480427d1/
      data r3/0.3286539d1,0.2060405d1,0.0d0/
      data r5/0.3275937d1,0.2498052d1,0.0d0/

      double precision r,phi,u13to,d13dph,d13dr

C Local variables
      integer i
      double precision sphi,cphi
      double precision sphi1,sphi2,sphi3,sphi4,sphi5
      double precision f1r,f3r,f5r
      double precision df1dr,df3dr,df5dr
      double precision g1r,g3r,g5r 

      u13to=0.0d0
      d13dph=0.0d0
      d13dr=0.0d0

      f1r=0.0d0
      f3r=0.0d0
      f5r=0.0d0
      df1rdr=0.0d0
      df3rdr=0.0d0
      df5rdr=0.0d0
      do i=1,ngf
        g1r=a1(i)*exp(-alpha1(i)*(r-r1(i))*(r-r1(i)))
        f1r=f1r+g1r
        df1rdr=df1rdr-2.0d0*alpha1(i)*(r-r1(i))*g1r
        g3r=a3(i)*exp(-alpha3(i)*(r-r3(i))*(r-r3(i)))
        f3r=f3r+g3r
        df3rdr=df3rdr-2.0d0*alpha3(i)*(r-r3(i))*g3r
        g5r=a5(i)*exp(-alpha5(i)*(r-r5(i))*(r-r5(i)))
        f5r=f5r+g5r
        df5rdr=df5rdr-2.0d0*alpha5(i)*(r-r5(i))*g5r
      enddo

      sphi=dsin(phi)   
      cphi=dcos(phi)
      sphi1=sphi
      sphi2=sphi*sphi1
      sphi3=sphi*sphi2
      sphi4=sphi*sphi3
      sphi5=sphi*sphi4   

      u13to=f1r*sphi1+f3r*sphi3+f5r*sphi5
      d13dph=f1r*cphi+3.0d0*f3r*sphi2*cphi+5.0d0*f5r*sphi4*cphi
      d13dr=df1rdr*sphi1+df3rdr*sphi3+df5rdr*sphi5

      end
      subroutine evu23to(r,phi,u23to,d23dph,d23dr)
***********************************************************************
* Little subroutine to evaluate U23 from C-C-O-H torsion
* U23(r,phi)=f2(r)*sin(2*phi)+f4(r)*sin(4*phi)+f6(r)*sin(6*phi)
* where fi(r) are expanded with a linear combination of Gaussians
***********************************************************************

      implicit double precision (a-h,o-z)

      integer ngf
      parameter (ngf=1)

      double precision scale
      parameter (scale=1.0d0)

      double precision a2(ngf),a4(ngf),a6(ngf)
      data a2/0.016924d0/
      data a4/0.007145d0/
      data a6/0.003113d0/          


      double precision alpha2(ngf),alpha4(ngf),alpha6(ngf)
      data alpha2/0.10841076d2/
      data alpha4/0.12574937d2/
      data alpha6/0.25329325d2/

      double precision r2(ngf),r4(ngf),r6(ngf)
      data r2/0.2151852d1/
      data r4/0.2197436d1/
      data r6/0.2285492d1/

      double precision r,phi,u23to,d23dph,d23dr

C Local variables
      integer i
      double precision sphi,cphi
      double precision sphi1,sphi2,sphi3,sphi4,sphi5
      double precision f2r,f4r,f6r
      double precision df2dr,df4dr,df6dr
      double precision g2r,g4r,g6r 

      u23to=0.0d0
      d23dph=0.0d0
      d23dr=0.0d0

      f2r=0.0d0
      f4r=0.0d0
      f6r=0.0d0
      df2rdr=0.0d0
      df4rdr=0.0d0
      df6rdr=0.0d0
      do i=1,ngf
        g2r=a2(i)*exp(-alpha2(i)*(r-r2(i))*(r-r2(i)))
        f2r=f2r+g2r
        df2rdr=df2rdr-2.0d0*alpha2(i)*(r-r2(i))*g2r
        g4r=a4(i)*exp(-alpha4(i)*(r-r4(i))*(r-r4(i)))
        f4r=f4r+g4r
        df4rdr=df4rdr-2.0d0*alpha4(i)*(r-r4(i))*g4r
        g6r=a6(i)*exp(-alpha6(i)*(r-r6(i))*(r-r6(i)))
        f6r=f6r+g6r
        df6rdr=df6rdr-2.0d0*alpha6(i)*(r-r6(i))*g6r
      enddo

      s2phi=dsin(2.0d0*phi)   
      c2phi=dcos(2.0d0*phi)
      s4phi=dsin(4.0d0*phi)
      c4phi=dcos(4.0d0*phi)
      s6phi=dsin(6.0d0*phi)
      c6phi=dcos(6.0d0*phi)

C      u23to=f2r*s2phi+f4r*s4phi+f6r*s6phi
C      d23dph=2.0d0*f2r*c2phi+4.0d0*f4r*c4phi+6.0d0*f6r*c6phi
C      d23dr=df2rdr*s2phi+df4rdr*s4phi+df6rdr*s6phi
      u23to=scale*(f2r*s2phi+f4r*s4phi+f6r*s6phi)
      d23dph=scale*(2.0d0*f2r*c2phi+4.0d0*f4r*c4phi+6.0d0*f6r*c6phi)
      d23dr=scale*(df2rdr*s2phi+df4rdr*s4phi+df6rdr*s6phi)

      end
      subroutine evintcoord(nbl,ibl,nba,iba,nto,ito,noop,ioop,
     $                      nint,RIC,X,Y,Z)
************************************************************************
* Evaluate the the interal coordinates for given connections
* and Cartesian coordinates
* nbl:          # of bond lengths
* nba:          # of bond angles
* nto:          # of torsions
* noop:         # of out-of-plane bendings
* nint:         # of internal coordinates
* ibl:          index of bond stretchings
* iba:          index of bond angles
* ito:          index of torsions
* ioop:         index of out-of-plane bendings
* RIC:          redundant internal coordinates
* X,Y,Z:        Cartesian coordinates
************************************************************************

      implicit double precision(a-h,o-z)


      integer maxint
      parameter (maxint=1000)
      double precision btoa
      parameter(btoa=0.5291772192d0)
      double precision pi
      parameter(pi=3.1415926535898d0)

      integer i
      integer nbl,nba,noop,nto,nint
      integer ibl(2,*),iba(3,*),ito(4,*),ioop(4,*)

      double precision X(*),Y(*),Z(*)
      double precision rab(maxint),theta(maxint),tors(maxint),
     $                 oopb(maxint)
      double precision RIC(*)

C Initialize RIC
      do i=1,nbl+nba+nto+noop
        RIC(i)=0.0d0
      enddo 
C Initialize rab,theta,tors, and oopb
      do i=1,maxint
        rab(i)=0.0d0
        theta(i)=0.0d0
        tors(i)=0.0d0
        oopb(i)=0.0d0
      enddo

C Evaluate bond lengths
      do i=1,nbl
        rab(i)=bndlen(ibl(1,i),ibl(2,i),X,Y,Z)
      enddo
      
C Evaluate bend angles
      do i=1,nba
        theta(i)=bndang(iba(1,i),iba(2,i),iba(3,i),X,Y,Z)
      enddo

C Evaluate torsions
      do i=1,nto
        tors(i)=dihedr(ito(1,i),ito(2,i),ito(3,i),ito(4,i),X,Y,Z)
      enddo

C Evaluate out-of-plane bending angles
      do i=1,noop
        oopb(i)=oopba1(ioop(1,i),ioop(2,i),ioop(3,i),ioop(4,i),X,Y,Z)
      enddo 

      nint=0
      do i=1,nbl
        nint=nint+1
        RIC(nint)=rab(i)
      enddo
      do i=1,nba
        nint=nint+1
        RIC(nint)=theta(i)
      enddo
      do i=1,nto
        nint=nint+1
        RIC(nint)=tors(i)
      enddo
      do i=1,noop
        nint=nint+1
        RIC(nint)=oopb(i)
      enddo
        
      end
      subroutine matvect(Y,m,A,lda,X,n)
***********************************************************************
* Subroutien to calculate matrix A multiply vector X
***********************************************************************

      integer lda,m,n
      double precision X(*),Y(*)
      double precision A(lda,*)

C Initialize Y
      do i=1,m
        y(i)=0.0d0
      enddo

      do i=1,m
        do j=1,n
          y(i)=y(i)+A(i,j)*x(j)
        enddo
      enddo

      end
      subroutine bmat3(nbl,ibl,nba,iba,nto,ito,noop,ioop,natm,
     $                X,Y,Z,bm,ldbm)
************************************************************************
* Subroutine to evaluate the Wilson B-matrix
************************************************************************

      implicit double precision(a-h,o-z)      

C Global variables
      integer ldbm
      integer nbl,nba,noop,nto,natm
      integer ibl(2,nbl),iba(3,nba),ito(4,nto),ioop(4,noop)
      double precision X(*),Y(*),Z(*)
      double precision bm(ldbm,*)

C Local variables
      integer i,j,k,l,ij,ijk,ijkl,irow
      double precision eij(3),ejk(3),ekl(3),eji(3),ejl(3),u(3),v(3)
      double precision w(3)

C Initialize B matrix
      do i=1,ldbm
        do j=1,3*natm
          BM(i,j)=0.0d0
        enddo
      enddo
      
      irow=0
C Evaluate the B matrix contributed from bond stretchings
      do ij=1,nbl
        irow=irow+1
        i=ibl(1,ij)
        j=ibl(2,ij)

        rij=bndlen(i,j,X,Y,Z)

        eij(1)=(X(j)-X(i))/rij
        eij(2)=(Y(j)-Y(i))/rij
        eij(3)=(Z(j)-Z(i))/rij

        BM(irow,3*i-2)=-eij(1)
        BM(irow,3*i-1)=-eij(2)
        BM(irow,3*i)  =-eij(3)
        BM(irow,3*j-2)= eij(1)
        BM(irow,3*j-1)= eij(2)
        BM(irow,3*j)  = eij(3)
      enddo
        
      
C Evaluate the B matrix contributed from bendings
      do ijk=1,nba
        irow=irow+1
        i=iba(1,ijk)
        j=iba(2,ijk)
        k=iba(3,ijk)

        rji=bndlen(j,i,X,Y,Z)
        rjk=bndlen(j,k,X,Y,Z)
        aijk=bndang(i,j,k,X,Y,Z)
        cijk=cos(aijk)
        sijk=sin(aijk)

        eji(1)=(X(i)-X(j))/rji
        eji(2)=(Y(i)-Y(j))/rji
        eji(3)=(Z(i)-Z(j))/rji
        ejk(1)=(X(k)-X(j))/rjk
        ejk(2)=(Y(k)-Y(j))/rjk
        ejk(3)=(Z(k)-Z(j))/rjk
        
        BM(irow,3*i-2)=(cijk*eji(1)-ejk(1))/(rji*sijk)
        BM(irow,3*i-1)=(cijk*eji(2)-ejk(2))/(rji*sijk)
        BM(irow,3*i)  =(cijk*eji(3)-ejk(3))/(rji*sijk)
        BM(irow,3*k-2)=(cijk*ejk(1)-eji(1))/(rjk*sijk)
        BM(irow,3*k-1)=(cijk*ejk(2)-eji(2))/(rjk*sijk)
        BM(irow,3*k)  =(cijk*ejk(3)-eji(3))/(rjk*sijk)
        BM(irow,3*j-2)=-BM(irow,3*i-2)-BM(irow,3*k-2)
        BM(irow,3*j-1)=-BM(irow,3*i-1)-BM(irow,3*k-1)
        BM(irow,3*j)  =-BM(irow,3*i)  -BM(irow,3*k)
      enddo

C
C Evaluate the B matrix contributed from torsions
C The calculations follow I. H. Williams, J. Mol. Spectros. 66, 288 (1977)
C 
      do ijkl=1,nto
        irow=irow+1
        i=ito(1,ijkl)
        j=ito(2,ijkl)
        k=ito(3,ijkl)
        l=ito(4,ijkl)

        rij=bndlen(i,j,X,Y,Z)
        rjk=bndlen(j,k,X,Y,Z)
        rkl=bndlen(k,l,X,Y,Z)
        aijk=bndang(i,j,k,X,Y,Z)
        ajkl=bndang(j,k,l,X,Y,Z)
        cijk=cos(aijk)
        sijk=sin(aijk)
        cjkl=cos(ajkl)
        sjkl=sin(ajkl)


        eij(1)=(X(j)-X(i))/rij
        eij(2)=(Y(j)-Y(i))/rij
        eij(3)=(Z(j)-Z(i))/rij
        ejk(1)=(X(k)-X(j))/rjk
        ejk(2)=(Y(k)-Y(j))/rjk
        ejk(3)=(Z(k)-Z(j))/rjk
        ekl(1)=(X(l)-X(k))/rkl
        ekl(2)=(Y(l)-Y(k))/rkl
        ekl(3)=(Z(l)-Z(k))/rkl
C
C Caculate u=eij*ejk, v=ejk*ekl
C
        call xprod1(eij,ejk,u)
        call xprod1(ejk,ekl,v)
C 
        BM(irow,3*i-2)=-u(1)/(rij*sijk*sijk)
        BM(irow,3*i-1)=-u(2)/(rij*sijk*sijk)
        BM(irow,3*i)  =-u(3)/(rij*sijk*sijk)
        BM(irow,3*j-2)= u(1)*(rjk-rij*cijk)/(rij*rjk*sijk*sijk)
     $                 -v(1)*cjkl/(rjk*sjkl*sjkl)
        BM(irow,3*j-1)= u(2)*(rjk-rij*cijk)/(rij*rjk*sijk*sijk)
     $                 -v(2)*cjkl/(rjk*sjkl*sjkl)
        BM(irow,3*j)  = u(3)*(rjk-rij*cijk)/(rij*rjk*sijk*sijk)
     $                 -v(3)*cjkl/(rjk*sjkl*sjkl)
        BM(irow,3*k-2)=-v(1)*(rjk-rkl*cjkl)/(rjk*rkl*sjkl*sjkl)
     $                 +u(1)*cijk/(rjk*sijk*sijk)
        BM(irow,3*k-1)=-v(2)*(rjk-rkl*cjkl)/(rjk*rkl*sjkl*sjkl)
     $                 +u(2)*cijk/(rjk*sijk*sijk)
        BM(irow,3*k)  =-v(3)*(rjk-rkl*cjkl)/(rjk*rkl*sjkl*sjkl)
     $                 +u(3)*cijk/(rjk*sijk*sijk)
        BM(irow,3*l-2)= v(1)/(rkl*sjkl*sjkl)
        BM(irow,3*l-1)= v(2)/(rkl*sjkl*sjkl)
        BM(irow,3*l)  = v(3)/(rkl*sjkl*sjkl)
      enddo

C
C Evaluate the B matrix contributed from out-of-plane bending 
C The calculation follows D. F. Mclntosh and K. H. Michaelian,
C Can. J. Spectros. 24, 35-40 (1979).
C
      do ijkl=1,noop
        irow=irow+1
        i=ioop(1,ijkl)
        j=ioop(2,ijkl)
        k=ioop(3,ijkl)
        l=ioop(4,ijkl)

        rji=bndlen(j,i,X,Y,Z)
        rjk=bndlen(j,k,X,Y,Z)
        rjl=bndlen(j,l,X,Y,Z)
        ai=bndang(k,j,l,X,Y,Z)
        ak=bndang(i,j,l,X,Y,Z)
        al=bndang(i,j,k,X,Y,Z)
        theta=oopba1(i,j,k,l,X,Y,Z)
        ci=cos(ai)
        si=sin(ai)
        ck=cos(ak)
        sk=sin(ak)
        cl=cos(al)
        sl=sin(al)
        ctheta=cos(theta)
        ttheta=tan(theta)

        eji(1)=(X(i)-X(j))/rji
        eji(2)=(Y(i)-Y(j))/rji
        eji(3)=(Z(i)-Z(j))/rji
        ejk(1)=(X(k)-X(j))/rjk
        ejk(2)=(Y(k)-Y(j))/rjk
        ejk(3)=(Z(k)-Z(j))/rjk
        ejl(1)=(X(l)-X(j))/rjl
        ejl(2)=(Y(l)-Y(j))/rjl
        ejl(3)=(Z(l)-Z(j))/rjl

C
C Caculate u=ejk*ejl
C
        call xprod1(ejk,ejl,u)
C
        BM(irow,3*i-2)=(1.0d0/rji)*(u(1)/(ctheta*si)-ttheta*eji(1))
        BM(irow,3*i-1)=(1.0d0/rji)*(u(2)/(ctheta*si)-ttheta*eji(2))
        BM(irow,3*i)  =(1.0d0/rji)*(u(3)/(ctheta*si)-ttheta*eji(3))
        BM(irow,3*k-2)=(1.0d0/rjk)*u(1)*(ci*ck-cl)/(ctheta*si*si*si)
        BM(irow,3*k-1)=(1.0d0/rjk)*u(2)*(ci*ck-cl)/(ctheta*si*si*si)
        BM(irow,3*k)  =(1.0d0/rjk)*u(3)*(ci*ck-cl)/(ctheta*si*si*si)
        BM(irow,3*l-2)=(1.0d0/rjl)*u(1)*(ci*cl-ck)/(ctheta*si*si*si)
        BM(irow,3*l-1)=(1.0d0/rjl)*u(2)*(ci*cl-ck)/(ctheta*si*si*si)
        BM(irow,3*l)  =(1.0d0/rjl)*u(3)*(ci*cl-ck)/(ctheta*si*si*si)
        BM(irow,3*j-2)=-BM(irow,3*i-2)-BM(irow,3*k-2)-BM(irow,3*l-2)
        BM(irow,3*j-1)=-BM(irow,3*i-1)-BM(irow,3*k-1)-BM(irow,3*l-1)
        BM(irow,3*j)  =-BM(irow,3*i)  -BM(irow,3*k)  -BM(irow,3*l)
      enddo

      end
      subroutine matmult(C,ldc,A,lda,m,B,ldb,n,l)

C Global variables
      integer lda,ldb,ldc
      integer m,n,l
      double precision A(lda,*),B(ldb,*),C(ldc,*)

C Local variables
      integer i,j,k

      do i=1,m
        do j=1,l
          C(i,j)=0.0d0
        enddo
      enddo
 
      do i=1,m
        do j=1,l
          do k=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          enddo
        enddo
      enddo

      end
      double precision function bndlen(i,j,X,Y,Z)
************************************************************************
* Function to evaluate the bond length between atom i and j
* Input: i,j (indices of 2 atoms)
************************************************************************

      implicit double precision(a-h,o-z)      

      double precision tiny
      parameter(tiny=1.0d-13)     

      double precision X(*),Y(*),Z(*)

      bndlen=sqrt( (x(i)-x(j))*(x(i)-x(j))
     $           + (y(i)-y(j))*(y(i)-y(j))
     $           + (z(i)-z(j))*(z(i)-z(j)))

      if (bndlen .le. tiny) then
        write(*,9999) i,j
CDB  KRY 20140718
        write(*,'(I3,3F18.12)') i,x(i),y(i),z(i)
        write(*,'(I3,3F18.12)') j,x(j),y(j),z(j)
        stop
      endif

 9999 Format(' Atom',I3,' and',I3,' are too close to each other!')

      end
      double precision function bndang(i,j,k,X,Y,Z)
************************************************************************
* Evaluate the bond angle between atoms i, j, and k
* Input: i,j,k (indices of 3 atoms, i-j and j-k bonds)
************************************************************************

      implicit double precision(a-h,o-z)      

      double precision X(*),Y(*),Z(*)
      double precision eji(3),ejk(3),rji,rjk

      rji=bndlen(i,j,X,Y,Z)
      rjk=bndlen(j,k,X,Y,Z)

      eji(1)=(X(i)-X(j))/rji
      eji(2)=(Y(i)-Y(j))/rji
      eji(3)=(Z(i)-Z(j))/rji

      ejk(1)=(X(k)-X(j))/rjk
      ejk(2)=(Y(k)-Y(j))/rjk
      ejk(3)=(Z(k)-Z(j))/rjk

      call sprod(3,eji,ejk,a)
      bndang=acos(a)

      end
      double precision function dihedr(i,j,k,l,X,Y,Z)
************************************************************************
* Evaluate the dihedral angle between atoms i, j, k, and l
* Input: i,j,k,l (indices of 4 atoms)
*        i--j--k--l
************************************************************************

      implicit double precision(a-h,o-z)      

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk

C   eij
C i-----j
C   aijk \ ejk
C         \ 
C          \ ajkl
C           k----l
C             kl   
      rij=bndlen(i,j,X,Y,Z)
      rjk=bndlen(j,k,X,Y,Z)
      rkl=bndlen(k,l,X,Y,Z)
      aijk=bndang(i,j,k,X,Y,Z)
      ajkl=bndang(j,k,l,X,Y,Z) 

C Check whether i-j-k or j-k-l are collinear
      if ((abs(aijk) .le. tiny) .or. (abs(aijk-pi) .lt. tiny)) then 
        write(*,9999) i,j,k
      endif
      if ((abs(ajkl) .le. tiny) .or. (abs(ajkl-pi) .lt. tiny)) then
        write(*,9999) j,k,l
      endif
 9999 Format('Torsion angle not defined since atoms',3I3,' collinear')

C Calculate eij
      eij(1)=(X(j)-X(i))/rij
      eij(2)=(Y(j)-Y(i))/rij
      eij(3)=(Z(j)-Z(i))/rij
C Calculate ejk
      ejk(1)=(X(k)-X(j))/rjk
      ejk(2)=(Y(k)-Y(j))/rjk
      ejk(3)=(Z(k)-Z(j))/rjk
C Calculate ekl
      ekl(1)=(X(l)-X(k))/rkl
      ekl(2)=(Y(l)-Y(k))/rkl
      ekl(3)=(Z(l)-Z(k))/rkl

C
C u = eij*ejk normal vector in plane i-j-k
C v = ejk*ekl Normal vector in plane j-k-l
C w = u*v used to determine the sign of dihedral angle
C    
      call xprod1(eij,ejk,u)
      call xprod1(ejk,ekl,v)
      call xprod1(u,v,w)

      call sprod(3,u,v,a)

      a=a/(dsin(aijk)*dsin(ajkl))
      a=max(a,-1.0d0)
      a=min(a,1.0d0)
      dihedr=dacos(a)
      
C Determine the sign of dihedr based on u,v,ejk
C dihedr >= 0, u,v,ejk right-hand (w*ejk >=0)
C dihedr < 0,  u,v,ejk left-hand (w*ejk < 0)

      wejk=w(1)*ejk(1)+w(2)*ejk(2)+w(3)*ejk(3)
      if (wejk .lt. 0.0d0) then
        dihedr=-dihedr
      endif

C Setting the range of dihedral angle between -90.0 degree to 270 degree
C to avoid discontinuity around 180 degrees
C      if (dihedr .lt. -0.5d0*pi) then
C        dihedr=2.0d0*pi+dihedr
C      endif

      end
      double precision function oopba1(i,j,k,l,X,Y,Z)
************************************************************************
* Evaluate the out-of-plane bending angle between atoms i, j, k, and l
* Input: i,j,k,l (indices of 4 atoms)
*             l
*            /
*     i-----j
*            \
*             k
************************************************************************

      implicit double precision(a-h,o-z)      

      double precision tiny
      data tiny/1.0d-30/

      double precision X(*),Y(*),Z(*)
      double precision eji(3),ejk(3),ejl(3),u(3)
      double precision v(3),w(3)
      double precision a1,a2
      double precision a,ai,si,rji,rjk,rjl
      
C               l
C              /
C          ak /
C     i------j ai != 180 deg
C          al \
C              \
C               k
      rji=bndlen(j,i,X,Y,Z)
      eji(1)=(X(i)-X(j))/rji
      eji(2)=(Y(i)-Y(j))/rji
      eji(3)=(Z(i)-Z(j))/rji

      rjk=bndlen(j,k,X,Y,Z)
      ejk(1)=(X(k)-X(j))/rjk
      ejk(2)=(Y(k)-Y(j))/rjk
      ejk(3)=(Z(k)-Z(j))/rjk

      rjl=bndlen(j,l,X,Y,Z)
      ejl(1)=(X(l)-X(j))/rjl
      ejl(2)=(Y(l)-Y(j))/rjl
      ejl(3)=(Z(l)-Z(j))/rjl

      ai=bndang(k,j,l,X,Y,Z)

      si=sin(ai)

      call xprod1(ejk,ejl,u)
C Check if whether ai = 180
      if (u(1)**2+u(2)**2+u(3)**2 .le. tiny) then
        write(*,*) "Looks like k-j-l collinear."
        write(*,*) "Out-of-plane bending angle is not defined."
        stop
      endif
      
      call sprod(3,u,eji,a)
      a=a/si
      a=max(a,-1.0d0)
      a=min(a,1.0d0)
      oopba1=asin(a)
      
C Check whether needs to go beyond -pi/2 to pi/2
C      write(*,*) "before correction, oopba1 = ",oopba1

C      do n=1,3
C        v(n)=ejk(n)+ejl(n)
C      enddo
C      call sprod(3,eji,v,b)

C      if ((oopba1 .lt. 0.0d0) .and. (b .gt. 0.0d0)) then
C        oopba1=2.0d0*asin(-1.0d0)-oopba1
C      endif
C      if ((oopba1 .ge. 0.0d0) .and. (b .gt. 0.0d0)) then
C        oopba1=2.0d0*asin(1.0d0)-oopba1
C      endif

      end
      subroutine xprod1(x1,x2,x3)
************************************************************************
* Calculate cross-product x1 and x2 and return x3
* Input: X1,X2
* Output: X3
************************************************************************

       implicit double precision (a-h,o-z)

       double precision x1(3),x2(3),x3(3)

       x3(1)=x1(2)*x2(3)-x1(3)*x2(2)
       x3(2)=x1(3)*x2(1)-x1(1)*x2(3)
       x3(3)=x1(1)*x2(2)-x1(2)*x2(1)

       end
      subroutine sprod(n,a,b,c)
************************************************************************
* Calculate scalar-product a and b and return c
* Input: a(n),b(n)
* Output: c
************************************************************************

       implicit double precision (a-h,o-z)

       integer n
       double precision a(*),b(*),c

       c=0.0d0
       do i=1,n
         c=c+a(i)*b(i)
       enddo

       end
