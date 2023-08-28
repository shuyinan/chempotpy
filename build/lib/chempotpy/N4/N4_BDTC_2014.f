      subroutine pes(x,igrad,path,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=4
      character(len=1024), intent(in) :: path
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: x1(4), y1(4), z1(4)
      double precision :: dx1(4), dy1(4), dz1(4)
      double precision :: v
      logical, save :: first_time_data=.true.
      integer :: i, j, k, l
      integer :: istate, iatom, idir

      !initialize 
      v=0.d0
      p=0.d0
      g=0.d0
      d=0.d0


      x1(:)=x(:,1)/0.529177211
      y1(:)=x(:,2)/0.529177211
      z1(:)=x(:,3)/0.529177211

      if(first_time_data) then
      call prepot(path)
      first_time_data=.false.
      endif

      call pot(x1,y1,z1,v,dx1,dy1,dz1)

      if (igrad==0) then 
        do istate=1,nstates
          p(istate)=v*27.211386
        enddo
      else if (igrad==1) then
        do istate=1,nstates
          p(istate)=v*27.211386
        enddo
        do istate=1,nstates
          g(istate,:,1)=dx1(:)*51.422067
          g(istate,:,2)=dy1(:)*51.422067
          g(istate,:,3)=dz1(:)*51.422067
        enddo
      else
        write(*,*) "Only energy and gradients are available"
      endif

      endsubroutine

C   System:                     N4
C   Functional form:            L-IMLS-G2 local method
C   Common name:                N4(adiabatic ground state)
C   Number of derivatives:      1
C   Number of bodies:           4
C   Number of electronic surfaces: 1
C   Interface: Section-2
C
C   References: Jason D. Bender, Sriram Doraiswamy, 
C   Donald G. Truhlar, and Graham V. Candler, J. Chem. Phys.,
C    forthcoming
C
C   Notes:    PES for N4, with special emphasis for 
C             N2 + N2 --> N2 + N + N.  The surface is constructed 
C             with an L-IMLS-based local method that incorporates
C             separation of pairwise interactions, permutational
C             invariance in the basis functions and weight 
C             functions, and a statistically-correlated cutoff 
C             radius. Note that the subroutine 'prepot' must be 
C             called once before using the subroutine 'pot'.  
C             'prepot' reads important parameters from a file
C             'pes.dat', which is also provided.
C
C   Input: x_(4), y_(4), z_4,           in unit of Bohr
C   Output: e_,                         in unit of Hartree
C   Output: dedx_, dedy_, dedz_,        in unit of Hartree/Bohr
C

c     Downloaded from the POTLIB library at
c       http://comp.chem.umn.edu/potlib/

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     !!! NOTE:
c     !!! The following module code is a modified version of code
c     !!! from the full FALCONS program.
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     FALCONS: Fitting Algorithm for the Local Construction Of 
c     eNergy Surfaces
c     Jason Bender, University of Minnesota, bende194@umn.edu
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Modules of common variables.
      module falcons_vars_eval
      implicit none
      real*8, dimension(:,:), allocatable :: qi,aji
      end module falcons_vars_eval

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module falcons_vars_main
      implicit none
      integer :: i,j,k,l,nn,ind,n,m,eta,np,moltype,pesform,
     & nma,nplmax,nall,nmincrc,mcrc,qps,pwonly,cgmeth,
     & ntarmin,ntarmax,ntar,crcmult,nitcrc
      integer, dimension(6) :: sd
      integer, dimension(4,3) :: cmap
      integer, dimension(24,6) :: qp
      real*8 :: pnu,epsnu,ps,pu,epsu,e,epw,systo,r,rwfa,rwfb
      real*8, dimension(6) :: t,qr,qs,q,xr,x,aa,qeq
      real*8, dimension(4,3) :: y,dedy
      real*8, dimension(:), allocatable :: crc
      end module falcons_vars_main

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Preliminary subroutine.  Call 'prepot' once before any call 
c     of 'pot'.
      subroutine prepot(path)
      use falcons_vars_eval
      use falcons_vars_main
      implicit none
      character(len=1024), intent(in) :: path
      call readpes(path)
      end subroutine prepot

      subroutine readpes(path)
      use falcons_vars_eval
      use falcons_vars_main
      implicit none
      character(len=1024), intent(in) :: path
      character(len=1024) :: file_path1
      logical :: ex
      integer :: oridval,maskval
      integer, parameter :: fid=200
      cmap(1,:)= (/1,2,3/)
      cmap(2,:)= (/1,4,5/)
      cmap(3,:)= (/2,4,6/)
      cmap(4,:)= (/3,5,6/)
c     A4 permutational symmetry.
      qps= 24
      qp(1 ,:)= (/1,2,3,4,5,6/)
      qp(2 ,:)= (/1,3,2,5,4,6/)
      qp(3 ,:)= (/2,1,3,4,6,5/)
      qp(4 ,:)= (/2,3,1,6,4,5/)
      qp(5 ,:)= (/3,1,2,5,6,4/)
      qp(6 ,:)= (/3,2,1,6,5,4/)
      qp(7 ,:)= (/1,4,5,2,3,6/)
      qp(8 ,:)= (/1,5,4,3,2,6/)
      qp(9 ,:)= (/4,1,5,2,6,3/)
      qp(10,:)= (/4,5,1,6,2,3/)
      qp(11,:)= (/5,1,4,3,6,2/)
      qp(12,:)= (/5,4,1,6,3,2/)
      qp(13,:)= (/2,4,6,1,3,5/)
      qp(14,:)= (/2,6,4,3,1,5/)
      qp(15,:)= (/4,2,6,1,5,3/)
      qp(16,:)= (/4,6,2,5,1,3/)
      qp(17,:)= (/6,2,4,3,5,1/)
      qp(18,:)= (/6,4,2,5,3,1/)
      qp(19,:)= (/3,5,6,1,2,4/)
      qp(20,:)= (/3,6,5,2,1,4/)
      qp(21,:)= (/5,3,6,1,4,2/)
      qp(22,:)= (/5,6,3,4,1,2/)
      qp(23,:)= (/6,3,5,2,4,1/)
      qp(24,:)= (/6,5,3,4,2,1/)
      ps= 4d0
c     Next read the pes.dat file.
      file_path1 = trim(path)//"/N4/pes.dat"
      inquire(file=file_path1,exist=ex)
      if (.not. ex) then
       write(*,'(a)') 'ERROR: No pes.dat file found.'
       stop
      end if
c     Read the parameter information from the first several lines of
c     the file.
      open(fid,file=file_path1,action='read')
      read(fid,*)
      read(fid,*)
      read(fid,*)
      read(fid,*)
      read(fid,*)
c     Skip over the header lines
      do i=1,20
       read(fid,*)
      end do
c     End of header lines
      read(fid,*)
      read(fid,*)
      read(fid,*)
      read(fid,*) moltype,pesform
      read(fid,*)
      read(fid,*) n,m,eta
      read(fid,*)
      read(fid,*) pnu,epsnu
      read(fid,*)
      read(fid,*) (aa(l),l=1,6)
      read(fid,*)
      read(fid,*) (qeq(l),l=1,6)
      read(fid,*)
      read(fid,*)
      read(fid,*)
      read(fid,*) nma,nplmax
      read(fid,*)
      call readcrcparameters(fid)
      read(fid,*)
c     Read the data points.
      nall= n+nma
      allocate(qi(n,6),aji(n,m))
      read(fid,*)
      read(fid,*)
      j= 0
      if (pesform== 1) then
c      Use default coordinate system format
       do i=1,nall
        read(fid,*) oridval,e,(t(l),l=1,6),maskval
        if (maskval== 0) then
         j= j+1
         call default2cartesian(
     &    t(1),t(2),t(3),t(4),t(5),t(6),y)
         call cartesian2rawdist(y,qr)
         call rawdist2rawmorse(qr,qi(j,:),qeq,aa)
        end if
       end do
      else
       write(*,'(a)') 
     &  'ERROR: pesform /= 1 is not supported in this code.'
       stop
      end if
      if (j /= n) then
       write(*,'(a)') 'ERROR: Bad pes.dat file.'//
     & ' Inconsistent number of unmasked points.'
       stop
      end if
c     Read the polynomial coefficients aji.
      read(fid,*)
      do i= 1,n
       do j= 1,m
        read(fid,*) aji(i,j)
       end do
      end do
      close(fid)
      end subroutine readpes

      subroutine readcrcparameters(fid)
      use falcons_vars_eval
      use falcons_vars_main
      implicit none
      integer :: fid
      read(fid,*)
      read(fid,*)
      read(fid,*) nmincrc,mcrc
      if (allocated(crc)) then
c      Clear the crc vector, if it was already defined earlier. This
c      was added as a safety check, when conducting multiple tests
c      in succession.
       deallocate(crc)
      end if
      allocate(crc(mcrc))
      read(fid,*)
      read(fid,*) (crc(l),l=1,mcrc)
      read(fid,*)
      read(fid,*) ntarmin,ntarmax,ntar,crcmult,nitcrc,rwfa,rwfb
      end subroutine readcrcparameters

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Main PES subroutine.
      subroutine pot(x_,y_,z_,e_,dedx_,dedy_,dedz_)
      use falcons_vars_eval
      use falcons_vars_main
      implicit none
      real*8, dimension(4), intent(in) :: x_,y_,z_
      real*8 :: cconv,econv,gconv,eref
      real*8, dimension(4) :: x__,y__,z__
      real*8, intent(out) :: e_
      real*8, dimension(4), intent(out) :: dedx_,dedy_,dedz_

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Set core fitting parameters.
      pu= 3d0
      epsu= 1d-1
      pwonly= 0
      cgmeth= 1 
      !cgmeth= 0: Do not compute gradient, which is instead set to 0.
      !cgmeth= 1: Compute gradient analytically.
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      cconv= 0.52917721092d0
      econv= 0.159360144d-2
      gconv= 0.843297564d-3
c     Reference energy of infinitely separated N2+N2 in Hartree, 
c     (from Yuliya Paukku).
      eref= -218.4080132d0
c     Convert x_, y_, and z_ from Bohr to A
      x__= x_*cconv
      y__= y_*cconv
      z__= z_*cconv
c     Assemble the FALCONS array y
      do l=1,4
       y(l,1)= x__(l)
       y(l,2)= y__(l)
       y(l,3)= z__(l)
      end do
c     Evaluate the potential energy surface fitting function
      call calccoordinates(y,qr,xr,q,x,sd,aa,qeq,moltype,qp,qps)
      call evaluatesurface(cgmeth,y,qr,xr,q,x,moltype,pwonly,
     & cmap,qp,qps,mcrc,crc,n,m,eta,ps,pu,epsu,aa,r,np,e,dedy) 
c     Store the results from the FALCONS variables into e_, dedx_,
c     dedy_, and dedz_. 
      e_= e
      dedx_= dedy(:,1)
      dedy_= dedy(:,2)
      dedz_= dedy(:,3)
c     Convert e_ from kcal/mol to Hartrees
      e_= e_*econv+ eref
c     Convert dedx_, dedy_, and dedz_ from [kcal mol^-1 A^-1] to
c     Hartree/Bohr
      dedx_= dedx_*gconv
      dedy_= dedy_*gconv
      dedz_= dedz_*gconv
      end subroutine pot

      subroutine evaluatesurface(cg,y,qr,xr,q,x,moltype,pwonly,
     & cmap,qp,qps,mcrc,crc,n,m,eta,ps,pu,epsu,aa,r,np,e,dedy)
      use falcons_vars_eval
      implicit none
c     Input variables
      integer, intent(in) :: cg,moltype,pwonly,qps,mcrc,n,m,eta
      integer, dimension(4,3), intent(in) :: cmap
      integer, dimension(qps,6), intent(in) :: qp
      real*8, intent(in) :: ps,pu,epsu
      real*8, dimension(6), intent(in) :: qr,xr,q,x,aa
      real*8, dimension(mcrc), intent(in) :: crc
      real*8, dimension(4,3), intent(in) :: y
c     Output variables
      integer, intent(out) :: np
      real*8, intent(out) :: r,e
      real*8, dimension(4,3), intent(out) :: dedy
c     Internal variables
      integer :: i,j,k,l,a,b,c
      real*8 :: rs,drdqch,dis,wsum,wi,pi,dwidz,tmp,epw
      real*8, dimension(6) :: qii,dwidxr,dwdxrsum,dpidx,dpidxr,
     & dedxrpart,dedxr,dedqr,ddisdq,depwdqr
      real*8, dimension(m) :: aj,bj 
      real*8, dimension(6,6) :: dxdxr
      real*8, dimension(m,6) :: dbjdx
      real*8, dimension(6,4,3) :: dqrdy
c     Initialize values. NOTE: It is essential to set every one of
c     the following quantities to 0; each is incremented an
c     arbitrary number of times.
      np=   0 
      e=    0d0
      wsum= 0d0 !sum of the weight functions
      do a=1,6
       dedxrpart(a)= 0d0 
       dwdxrsum(a)=  0d0 !sums of the weight function derivatives
      end do
c     Compute the cutoff radius. 
      call calcqch(q,tmp)
      call calccr(tmp,mcrc,crc,r,rs,drdqch)
c     Calculate polynomial basis functions and derivatives
      call calcbj(bj,dbjdx,x,eta,m,1)
c     Calculate the Jacobian matrices dxdxr and dqrdy, if analytic
c     gradients are to be computed.
      if (cg /= 0) then
       call calcinvariantsjacobian(xr,x,dxdxr,moltype)
       call calcdqrdyjacobian(qr,y,dqrdy)
      end if

c     Calculate the pairwise (i.e., two-body) interaction energy
      call calcsumpairwise(qr,epw,depwdqr,moltype,cg)

      if (pwonly == 0) then
c     Begin loop over all data points. 
      do i=1,n

c      Compute the distance between the specified point q and each
c      of the data points. A data point will only be used if its 
c      distance from q is less than or equal to the cutoff radius.
       do a=1,6
        qii(a)= qi(i,a)
       end do
       call distance6ds_perm(q,qii,qp,qps,0,dis,ddisdq)

       if (dis <= rs) then
c       The data point i lies within the cutoff radius.  Increment
c       the count np.
        np= np+1  
c       Re-compute the distance, if analytic gradients are to be
c       calculated.
        if (cg /= 0) then
         call distance6ds_perm(q,qii,qp,qps,cg,dis,ddisdq)
        end if
        do j=1,m
         aj(j)= aji(i,j)
        end do
c       Compute the weight function and its derivative.
        call weightw(wi,dwidz,dis,rs,ps,pu,epsu)
c       Compute the local fit pi
        pi= 0d0
        do j=1,m
         pi= pi+ aj(j)*bj(j)
        end do
c       Increment the running sum of the energy e and the running
c       sum of the weight functions wsum.
        wsum= wsum+ wi
        e= e+ wi*pi
c ----------- Additional computations, for analytic gradients
        if (cg /= 0) then
c        Compute derivatives of the local fits dpidx
         do a=1,6
          dpidx(a)= 0d0
         end do
         do j=1,m
          do a=1,6
           dpidx(a)= dpidx(a)+ aj(j)*dbjdx(j,a)
          end do
         end do
c        Compute the derivatives of the local fits dpidxr, 
c        by matrix multiplication with the Jacobian matrix dxdxr.
         do a=1,6
          tmp= 0d0
          do k=1,6
           tmp= tmp+ dpidx(k)*dxdxr(k,a)
          end do
          dpidxr(a)= tmp
         end do
c        Compute the derivatives of the weight functions dwidxr,
c        and increment the running sums of the derivatives of the
c        weight functions.
         tmp= dis*4d0*drdqch/r
         do a=1,6
          dwidxr(a)= dwidz*(ddisdq(a)- tmp*q(a))/rs
          dwdxrsum(a)= dwdxrsum(a)+ dwidxr(a)
         end do
c        Increment the running sums for the energy e and its
c        derivatives with respect to xr(a).
c        The derivatives are computed in two terms, which will be
c        appropriately scaled later.
         do a=1,6
          dedxrpart(a)= dedxrpart(a)+ wi*dpidxr(a)+ pi*dwidxr(a) 
         end do
        end if
c ----------- End additional computations, for analytic gradients
       end if
      end do  !end loop over data points within the cutoff radius
      end if  !end if 'pwonly==0' computation

c     If wsum==0d0, then no local fits were used for the
c     evaluation, (either because no points were found within the
c     cutoff radius, or because pwonly was set to a non-zero value,
c     indicating that only pairwise interactions should be 
c     computed).  Set wsum= 1d0 for safety.  
      if (wsum== 0d0) then
       wsum= 1d0
      end if

c     Normalize the energy e.
      e= e/wsum

c --------- Additional computations, for analytic gradients
      if (cg /= 0) then
c      Normalize the derivatives of the energy dedxr
       do a=1,6
        dedxr(a)= (dedxrpart(a)- e*dwdxrsum(a))/wsum
       end do
c      Compute the derivatives of the energy dedqr.  Include the
c      component of the derivative of the pairwise interaction energy
c      depwdqr.
       do a=1,6
        dedqr(a)= -dedxr(a)*xr(a)/aa(a)+ depwdqr(a)
       end do
c      Compute the derivatives dedqr, using the Jacobian matrix dqrdy.
       do b=1,4 
        do c=1,3
         tmp= 0d0
         do k=1,3
c         Extract those indices for which dqrdy /= 0
          a= cmap(b,k)
          tmp= tmp+ dedqr(a)*dqrdy(a,b,c)
         end do
         dedy(b,c)= tmp
        end do
       end do
c --------- End additional computations, for analytic gradients
      else
       do b=1,4 
        do c=1,3
         dedy(b,c)= 0d0
        end do
       end do
      end if

c     At this point, 'e' is equal to the many-body interaction 
c     energy. Add the pairwise interaction energy 'epw' to obtain
c     the total energy.
      e= e+ epw

c     If no points were found within the cutoff radius, report a
c     warning.
c     WHEN THIS HAPPENS, THE FIRST REMEDY TO TRY IS TO INCREASE
c     THE CUTOFF RADIUS.  
      if (np==0) then
       write(*,'(a)') 
     &  ' WARNING: No points found within cutoff radius for'//
     &  ' evaluation point:'
       write(*,'(a,6(1x,f15.5))') ' qr= ',(qr(l),l=1,6)
       write(*,'(a,6(1x,f15.5))') ' q = ',(q(l),l=1,6)
       write(*,'(a,6(1x,f15.5))') ' x = ',(x(l),l=1,6)
      end if
      end subroutine evaluatesurface
  
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Compute the characteristic value 'qch', used in the cutoff 
c     radius correlation, given a length-6 coordinate vector 'q'. 
      subroutine calcqch(q,qch)
      implicit none
      integer :: i
      real*8 :: val
      real*8, dimension(6), intent(in) :: q
      real*8, intent(out) :: qch
      val= 0d0
      do i=1,6
       val= val+ q(i)**2d0
      end do
      qch= val
      end subroutine calcqch
 
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Compute the cutoff radius 'cr', the square of the cutoff radius
c     'crs', and the derivative of the cutoff radius with respect
c     to the characteristic coordinate 'dcrdqch'. The subroutine
c     requires the cutoff radius correlation stored in the length-
c     (mcrc) vector 'crc'.
      subroutine calccr(qch,mcrc,crc,cr,crs,dcrdqch) 
      implicit none
      integer :: icrc
      integer, intent(in) :: mcrc
      real*8, intent(in) :: qch
      real*8 :: pqch
      real*8, intent(out) :: cr,crs,dcrdqch
      real*8, dimension(mcrc), intent(in) :: crc
      pqch= 1d0
      cr= crc(1)
      dcrdqch= 0d0
      do icrc= 2,mcrc
c      Here, pqch equals qch**(icrc-2)
       dcrdqch= dcrdqch+ (crc(icrc)*dble(icrc-1)*pqch)
       pqch= pqch*qch 
c      Now, pqch equals qch**(icrc-1)
       cr= cr+ (crc(icrc)* pqch)
      end do
      crs= cr*cr
      end subroutine calccr

      subroutine calcbj(b,dbdx,x,eta,m,dflag)
      implicit none
      integer :: ind,i,j,k,l,ll,a
      integer, intent(in) :: eta,m,dflag
      real*8 :: tmp,tmp2,tmp3
      real*8, dimension(6), intent(in) :: x
      real*8, dimension(m), intent(out) :: b
      real*8, dimension(m,6), intent(out) :: dbdx
      real*8, dimension(6) :: xx,x2,xxx
      real*8, dimension(m) :: dbdxa
      real*8, dimension(6,6) :: xy,xy2
      real*8, dimension(6,6,6) :: xyz
      b   = 0d0
      dbdx= 0d0
      xx  = 0d0 !xx(j)     = x(j)*x(j)
      x2  = 0d0 !x2(j)     = x(j)*2d0
      xxx = 0d0 !xxx(j)    = x(j)*x(j)*x(j)
      xy  = 0d0 !xy(j,k)   = x(j)*x(k)
      xy2 = 0d0 !xy2(j)    = x(j)*x(k)*2d0
      xyz = 0d0 !xyz(j,k,l)= x(j)*x(k)*x(l)
      tmp = 0d0
      tmp2= 0d0
      tmp3= 0d0
c     ---------- COMPUTATION OF B
      ind= 0
      do j=1,6
       ind= ind+1
       b(ind)= x(j) !linear terms
      end do
c     [ind= 6 at this point]
      if (eta .ge. 2) then
       do j= 1,6
        ind= ind+1
        tmp= x(j)*x(j)
        xx(j) = tmp
        b(ind)= tmp !quadratic terms such as x1^2
        do k= j+1,6
         ind= ind+1
         tmp= x(j)*x(k)
         xy(j,k)= tmp
         b(ind) = tmp
c        quadratic cross terms such as as x1*x2
        end do
       end do
c      [ind= 27 at this point]
      end if
      if (eta .ge. 3) then
       do j= 1,6
        ind= ind+1
        tmp= xx(j)
        tmp2= tmp*x(j)
        b(ind)= tmp2 !cubic terms such as x1^3
        xxx(j)= tmp2
        do k= j+1,6
         ind= ind+1
         b(ind)= tmp*x(k)
         ind= ind+1
         b(ind)= x(j)*xx(k)
c        cubic cross terms with a quadratic factor and a linear
c        factor such as (x1^2)*(x2)
        end do
       end do
c      [ind= 63 at this point]
       do j= 1,6
        do k= j+1,6
         tmp= xy(j,k)
         do l= k+1,6
          ind= ind+1
          tmp2= tmp*x(l)
          b(ind)= tmp2
          xyz(j,k,l)= tmp2
c         cubic cross terms with three linear factors, such as
c         x1*x2*x3
         end do
        end do
       end do
c      [ind= 83 at this point]
      end if
      if (eta .ge. 4) then
       do j= 1,6
        ind= ind+1
        tmp= xxx(j)
        b(ind)= tmp*x(j) !quartic terms such as x1^4
        do k= j+1,6
         ind= ind+1
         b(ind)= tmp*x(k)
         ind= ind+1
         b(ind)= x(j)*xxx(k)
c        quartic cross terms with a cubic factor and a linear factor
c        such as (x1^3)*(x2)
        end do
       end do
c      [ind= 119 at this point]
       do j= 1,6
        tmp= xx(j)
        do k= j+1,6
         tmp2= xx(k)
         ind= ind+1
         b(ind)= tmp*tmp2
c        quartic cross terms with two quadratic factors, such as
c        (x1^2)*(x2^2)
         tmp3= xy(j,k)
         do l= k+1,6
c         --- Term 1 of 3
          ind= ind+1
          b(ind)= tmp*xy(k,l)
c         --- Term 2 of 3
          ind= ind+1
          b(ind)= tmp2*xy(j,l)
c         --- Term 3 of 3
          ind= ind+1
          b(ind)= xx(l)*tmp3
c         quartic cross terms with a quadratic factor and two linear
c         factors, such as (x1^2)*x2*x3
         end do
        end do
       end do
c      [ind= 194 at this point]
       do j=1,6
        do k= j+1,6
         tmp= xy(j,k)
         do l= k+1,6
          do ll= l+1,6
           ind= ind+1
           b(ind)= tmp*xy(l,ll)
c          quartic cross terms with four linear factors, such as
c          x1*x2*x3*x4
          end do
         end do
        end do
       end do
c      [ind= 209 at this point]
      end if
c     ---------- COMPUTATION OF DBDX (optional)
c     For each a, a= 1, ..., 6, compute the derivative dbdxa.
      if (dflag /= 0) then
      do j=1,6
       x2(j)= 2d0*x(j)
       do k=j+1,6
        xy2(j,k)= 2d0*xy(j,k)
       end do
      end do
      do a=1,6
       do j=1,m
        dbdxa(j)= 0d0
       end do
       dbdxa(a)= 1d0 !derivative of linear term 
       ind= 6
c      [ind= 6 at this point]
       if (eta .ge. 2) then
        do j= 1,6
         ind= ind+1
         if (j == a) then
          dbdxa(ind)= x2(j)
c         derivative of quadratic terms such as x1^2
         end if
         do k= j+1,6
          ind= ind+1
          if (j == a) then
           dbdxa(ind)= x(k)
          else if (k == a) then
           dbdxa(ind)= x(j)
          end if
c         derivative of quadratic cross terms with two linear 
c         factors
         end do
        end do
c       [ind= 27 at this point]
       end if
       if (eta .ge. 3) then
        do j= 1,6
         ind= ind+1
         if (j == a) then
          dbdxa(ind)= 3d0*xx(a)
c         derivative of cubic terms such as x1^3
          do k= j+1,6
           ind= ind+1
           dbdxa(ind)= xy2(j,k)
           ind= ind+1
           dbdxa(ind)= xx(k)
c          derivatives of cubic cross terms with a quadratic factor
c          and a linear factor, such as (x1^2)*(x2)
          end do
         else 
          do k= j+1,6
           if (k == a) then
            ind= ind+1
            dbdxa(ind)= xx(j)
            ind= ind+1
            dbdxa(ind)= xy2(j,k)
c           derivatives of cubic cross terms with a quadratic factor
c           and a linear factor, such as (x1^2)*(x2)
           else
            ind= ind+2
           end if
          end do
         end if
        end do
c       [ind= 63 at this point]
        do j= 1,6
         do k= j+1,6
          do l= k+1,6
           ind= ind+1
           if (j == a) then
            dbdxa(ind)= xy(k,l)
           else if (k == a) then
            dbdxa(ind)= xy(j,l)
           else if (l == a) then
            dbdxa(ind)= xy(j,k)
           end if
c          derivatives of cubic cross terms with three linear
c          factors, such as x1*x2*x3
          end do
         end do
        end do
c       [ind= 83 at this point]
       end if
       if (eta .ge. 4) then
        do j= 1,6
         ind= ind+1
         if (j == a) then
          dbdxa(ind)= 4d0*xxx(a)
c         derivative of quartic terms such as x1^4
          do k= j+1,6
           ind= ind+1
           dbdxa(ind)= 3d0*xx(j)*x(k)
           ind= ind+1
           dbdxa(ind)= xxx(k)
c          derivatives of quartic cross terms with a cubic factor
c          and a linear factor, such as (x1^3)*(x2)
          end do
         else 
          do k= j+1,6
           if (k == a) then
            ind= ind+1
            dbdxa(ind)= xxx(j)
            ind= ind+1
            dbdxa(ind)= 3d0*x(j)*xx(k)
c           derivatives of quartic cross terms with a cubic factor
c           and a linear factor, such as (x1^3)*(x2)
           else
            ind= ind+2
           end if
          end do
         end if
        end do
c       [ind= 119 at this point]
        do j= 1,6
         tmp= xx(j)
         do k= j+1,6
          tmp2= xx(k)
          ind= ind+1
          if (j == a) then
           dbdxa(ind)= x2(j)*tmp2
          else if (k == a) then
           dbdxa(ind)= tmp*x2(k)
          end if
c         derivatives of quartic cross terms with two quadratic 
c         factors, such as (x1^2)*(x2^2)
          do l= k+1,6
c          --- Term 1 of 3
           ind= ind+1
           if (j == a) then
            dbdxa(ind)= x2(j)*xy(k,l)
           else if (k == a) then
            dbdxa(ind)= tmp*x(l)
           else if (l == a) then
            dbdxa(ind)= tmp*x(k)
           end if
c          derivatives of quartic cross terms with a quadratic factor
c          and two linear factors, such as (x1^2)*x2*x3
c          --- Term 2 of 3
           ind= ind+1
           if (j == a) then
            dbdxa(ind)= tmp2*x(l)
           else if (k == a) then
            dbdxa(ind)= x2(k)*xy(j,l)
           else if (l == a) then
            dbdxa(ind)= tmp2*x(j)
           end if
c          --- Term 3 of 3
           ind= ind+1
           if (j == a) then
            dbdxa(ind)= xx(l)*x(k)
           else if (k == a) then
            dbdxa(ind)= xx(l)*x(j)
           else if (l == a) then
            dbdxa(ind)= x2(l)*xy(j,k)
           end if
          end do
         end do
        end do
c       [ind= 194 at this point]
        do j= 1,6
         do k= j+1,6
          do l= k+1,6
           do ll= l+1,6
            ind= ind+1
            if (j == a) then
             dbdxa(ind)= xyz(k,l,ll)
            else if (k == a) then
             dbdxa(ind)= xyz(j,l,ll)
            else if (l == a) then
             dbdxa(ind)= xyz(j,k,ll)
            else if (ll == a) then
             dbdxa(ind)= xyz(j,k,l)
            end if
c           derivatives of quartic cross terms with four linear 
c           factors, such as x1*x2*x3*x4
           end do
          end do
         end do
        end do
c      [ind= 209 at this point]
       end if
       do j=1,m
        dbdx(j,a)= dbdxa(j)
       end do
      end do !end 'do a=1,6'
      end if !end 'if (dflag /= 0) then'
      end subroutine calcbj
    
      subroutine distance6ds_perm(q,qi,qp,qps,cg,ds,ddsdq)
      implicit none 
      integer :: i,a,iszero
      integer, intent(in) :: cg,qps
      integer, dimension(6) :: sd
      integer, dimension(qps,6), intent(in) :: qp
      real*8 :: dsval,dssum,mult,dsvali,dsvalis,dssumi
      real*8, intent(out) :: ds
      real*8, dimension(6) :: qdif
      real*8, dimension(6), intent(in) :: q,qi
      real*8, dimension(6), intent(out) :: ddsdq
c     Initialize variables
      iszero= 0 
      dssum= 0d0
      do a=1,6
       ddsdq(a)= 0d0
      end do
c     Loop through all permutations
      do i=1,qps
       do a=1,6
        sd(a)= qp(i,a)
       end do
       do a=1,6
        qdif(a)= q(a)- qi(sd(a))
       end do
       dsval= 0d0
       do a=1,6
        dsval= dsval+ qdif(a)**2d0
       end do
c      Increment running sums 'dssum' and 'ddsdq'
       if (dsval > 0d0) then
        dsvali = 1d0/dsval
        dsvalis= dsvali*dsvali !dsvalis= dsval^(-2)
        dssum= dssum+ dsvalis
        if (cg /= 0) then
         mult= 2d0*dsvali*dsvalis !mult= 2*dsval^(-3)
         do a=1,6
          ddsdq(a)= ddsdq(a)+ mult*qdif(a)
         end do
        end if
       else
c       A zero distance was found. Update the flag variable 'iszero'.
        iszero= 1
       end if
      end do
c     The loop over all permutations is complete. Calculate final
c     values 'ds' and 'ddsdq'.
      if (iszero /= 1) then
       dssumi = 1d0/dssum
       ds     = sqrt(dssumi) !ds= dssum^(-0.5)
       if (cg /= 0) then
        mult= dssumi*ds  !mult= dssum^(-1.5)
        do a=1,6
         ddsdq(a)= mult*ddsdq(a)
        end do
       end if
      else
c      In this case, 'ds' is zero, and all derivatives 'ddsdq' are
c      also zero.
       ds= 0d0
       if (cg /= 0) then
        do a=1,6
         ddsdq(a)= 0d0
        end do
       end if
      end if
      end subroutine distance6ds_perm

      subroutine weightw(w,dwdz,ds,rs,ps,pu,epsu)
      implicit none
      real*8 :: z,uval,sval
      real*8, intent(in) :: ds,ps,pu,epsu,rs
      real*8, intent(out) :: w,dwdz
      z= ds/rs
      uval= 1d0/((z**pu)+ (epsu**pu))
      sval= (1d0- z**ps)**4d0
      w= uval*sval
      dwdz= uval*4d0*((1d0-z**ps)**3d0)*(-ps*(z**(ps-1d0)))+
     & sval*(uval*uval)*(-(pu*(z**(pu-1d0))))
      end subroutine weightw

      subroutine calccoordinates(
     & y,qr,xr,q,x,sd,aa,qeq,moltype,qp,qps)
      implicit none
      integer, intent(in) :: moltype,qps
      integer, dimension(qps,6) :: qp
      integer, dimension(6), intent(out) :: sd
      real*8, dimension(6), intent(in) :: aa,qeq
      real*8, dimension(6), intent(out) :: qr,xr,q,x
      real*8, dimension(4,3) :: y
      call cartesian2rawdist(y,qr)
      call rawdist2rawmorse(qr,xr,qeq,aa)
      q= xr
      call calcinvariants(xr,x,moltype)
      end subroutine calccoordinates

      subroutine default2cartesian(ra,rb,d,thetaa_d,thetab_d,
     &  phi_d,y)
      implicit none
      real*8, intent(in) :: ra,rb,d,thetaa_d,thetab_d,phi_d
      real*8 :: thetaa,thetab,phi,pi,y11,y31,y32
      real*8, dimension(4,3), intent(out) :: y
      pi= 3.1415926536d0
      thetaa= thetaa_d*pi/180d0
      thetab= thetab_d*pi/180d0
      phi= phi_d*pi/180d0
      y11= ra*sin(thetaa)/2d0
      y(1,1)= y11
      y(1,2)= 0d0
      y(1,3)= (d- ra*cos(thetaa))/2d0
      y(2,1)= -y11
      y(2,2)= 0d0
      y(2,3)= (d+ ra*cos(thetaa))/2d0
      y31= rb*sin(thetab)*cos(phi)/2d0
      y(3,1)= y31
      y32= rb*sin(thetab)*sin(phi)/2d0
      y(3,2)= y32
      y(3,3)= (-d+ rb*cos(thetab))/2d0
      y(4,1)= -y31
      y(4,2)= -y32
      y(4,3)= -(d+ rb*cos(thetab))/2d0
      end subroutine default2cartesian

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     'Cartesian' to 'raw internuclear distance' coordinate 
c     transformation.
c     Converts from the Cartesian coordinates y to the raw 
c     internuclear distance coordinates,
c       qr= (r12, r13, r14, r23, r24, r34).
c     We assume that all distances are in [A].
      subroutine cartesian2rawdist(y,qr)
      implicit none
      integer :: i
      real*8, dimension(6), intent(out) :: qr
      real*8, dimension(4,3), intent(in) :: y
      call distance3d(y(1,:),y(2,:),qr(1))
      call distance3d(y(1,:),y(3,:),qr(2))
      call distance3d(y(1,:),y(4,:),qr(3))
      call distance3d(y(2,:),y(3,:),qr(4))
      call distance3d(y(2,:),y(4,:),qr(5))
      call distance3d(y(3,:),y(4,:),qr(6))
      end subroutine cartesian2rawdist

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Calculates the distance between two points in three-dimensional
c     space.
      subroutine distance3d(p1,p2,d)
      implicit none
      real*8, intent(out) :: d
      real*8, dimension(3), intent(in) :: p1,p2
      d=
     & sqrt((p1(1)-p2(1))**2d0+(p1(2)-p2(2))**2d0+(p1(3)-p2(3))**2d0)
      end subroutine distance3d

      subroutine rawdist2rawmorse(qr,xr,qeq,aa)
      implicit none
      integer :: i
      real*8, dimension(6), intent(in) :: qr,qeq,aa
      real*8, dimension(6), intent(out) :: xr
      do i=1,6
       xr(i)= exp(-(qr(i)-qeq(i))/aa(i))
      end do
      end subroutine rawdist2rawmorse

      subroutine calcinvariants(vr,v,moltype)
      implicit none
      integer, intent(in) :: moltype
      integer :: i
      real*8 :: e0,e1,f0,f1,vrtmp
      real*8, dimension(6), intent(in) :: vr
      real*8, dimension(6), intent(out) :: v
      real*8, dimension(4) :: vv
      v=0d0
      select case (moltype)
      case (1)
c      A4 permutational symmetry.
c      In the following scheme, v(1), v(2), and v(3) are
c      three-body interaction terms. v(4), v(5), and v(6)
c      are four-body interaction terms.
c      1st permutationally invaraint polynomial
       v(1)= (vr(1)*vr(2)+vr(1)*vr(3)+vr(1)*vr(4)+
     &        vr(1)*vr(5)+vr(2)*vr(3)+vr(2)*vr(4)+
     &        vr(2)*vr(6)+vr(3)*vr(5)+vr(3)*vr(6)+
     &        vr(4)*vr(5)+vr(4)*vr(6)+vr(5)*vr(6))/12d0
       v(1)= v(1)       !LEGACY: **0.5d0
c      2nd permutationally invariant polynomial
       v(2)= (vr(1)*vr(2)*vr(4)+ vr(1)*vr(3)*vr(5)+
     &        vr(2)*vr(3)*vr(6)+ vr(4)*vr(5)*vr(6))/4d0
       v(2)= v(2)       !LEGACY: **(1d0/3d0)
c      3nd permutationally invariant polynomial
       v(3)= (vr(1)*vr(2)*(vr(1)+vr(2))+
     &        vr(1)*vr(3)*(vr(1)+vr(3))+
     &        vr(1)*vr(4)*(vr(1)+vr(4))+
     &        vr(1)*vr(5)*(vr(1)+vr(5))+
     &        vr(2)*vr(3)*(vr(2)+vr(3))+
     &        vr(2)*vr(4)*(vr(2)+vr(4))+
     &        vr(2)*vr(6)*(vr(2)+vr(6))+
     &        vr(3)*vr(5)*(vr(3)+vr(5))+
     &        vr(3)*vr(6)*(vr(3)+vr(6))+
     &        vr(4)*vr(5)*(vr(4)+vr(5))+
     &        vr(4)*vr(6)*(vr(4)+vr(6))+
     &        vr(5)*vr(6)*(vr(5)+vr(6)))/24d0
       v(3)= v(3)       !LEGACY: **(1d0/3d0)
c      4th permutationally invariant polynomial
       v(4)= (vr(1)*vr(2)*vr(3)+ vr(1)*vr(4)*vr(5)+
     &        vr(2)*vr(4)*vr(6)+ vr(3)*vr(5)*vr(6))/4d0
       v(4)= v(4)       !LEGACY: **(1d0/3d0)
c      5th permutationally invariant polynomial
       v(5)= (vr(1)*vr(3)*vr(4)+ vr(2)*vr(3)*vr(4)+
     &        vr(1)*vr(2)*vr(5)+ vr(2)*vr(3)*vr(5)+
     &        vr(2)*vr(4)*vr(5)+ vr(3)*vr(4)*vr(5)+
     &        vr(1)*vr(2)*vr(6)+ vr(1)*vr(3)*vr(6)+
     &        vr(1)*vr(4)*vr(6)+ vr(3)*vr(4)*vr(6)+
     &        vr(1)*vr(5)*vr(6)+ vr(2)*vr(5)*vr(6))/12d0
       v(5)= v(5)       !LEGACY: **(1d0/3d0)
c      6th permutationally invariant polynomial
       v(6)= (vr(1)*vr(2)*vr(5)*vr(6)+
     &        vr(1)*vr(3)*vr(4)*vr(6)+
     &        vr(2)*vr(3)*vr(4)*vr(5))/3d0
       v(6)= v(6)       !LEGACY: **0.25d0

      case (2)
c      A2B2 permutational symmetry.
       v(1)= (vr(2)*vr(4)+ vr(3)*vr(5))/2d0
       v(2)= (vr(2)*vr(3)+ vr(4)*vr(5))/2d0
       v(3)= (vr(1)*(vr(2)+vr(3)+vr(4)+vr(5)))/4d0
       v(4)= (vr(6)*(vr(2)+vr(3)+vr(4)+vr(5)))/4d0
       v(5)= (vr(1)*(vr(2)*vr(3)+vr(4)*vr(5)))/2d0
       v(6)= (vr(6)*(vr(2)*vr(4)+vr(3)*vr(5)))/2d0
      end select
      end subroutine calcinvariants

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Given a raw coordinate vector 'vr' and the corresponding 
c     permutationally invariant coordinate vector 'v', computes the
c     (6x6) Jacobian matrix 'dvdvr', defined by:
c        dvdvr(i,j)= d( v(i) )/ d( vr(j) )
c     This is used for computing analytic gradients in the
c     'evaluatesurface' subroutine.
      subroutine calcinvariantsjacobian(vr,v,dvdvr,moltype)
      implicit none
      integer, intent(in) :: moltype
      integer :: i
      real*8 :: mult
      real*8, dimension(6) :: tmp,tmp2
      real*8, dimension(6), intent(in) :: vr,v
      real*8, dimension(6,6), intent(out) :: dvdvr
      select case (moltype)
      case (1)
c     A4 permutational symmetry.
c     derivatives of the 4th permutationally invaraint polynomial
      mult= 1d0/4d0     !LEGACY: *3d0*(v(4)**2d0))
      tmp(1)= vr(2)*vr(3)+ vr(4)*vr(5)
      tmp(2)= vr(1)*vr(3)+ vr(4)*vr(6)
      tmp(3)= vr(1)*vr(2)+ vr(5)*vr(6)
      tmp(4)= vr(1)*vr(5)+ vr(2)*vr(6)
      tmp(5)= vr(1)*vr(4)+ vr(3)*vr(6)
      tmp(6)= vr(2)*vr(4)+ vr(3)*vr(5)
      do i=1,6
       dvdvr(4,i)= tmp(i)*mult
      end do
c     derivatives of the 2nd permutationally invaraint polynomial
      mult= 1d0/4d0     !LEGACY: *3d0*(v(2)**2d0))
      do i=1,6
       dvdvr(2,i)= tmp(7-i)*mult
      end do
c     derivatives of the 1st permutationally invaraint polynomial
      mult= 1d0/12d0    !LEGACY: *2d0*v(1))
      tmp(1)= vr(2)+vr(3)+vr(4)+vr(5)
      tmp(2)= vr(1)+vr(3)+vr(4)+vr(6)
      tmp(3)= vr(1)+vr(2)+vr(5)+vr(6)
      tmp(4)= tmp(3)
      tmp(5)= tmp(2)
      tmp(6)= tmp(1)
      do i=1,6
       dvdvr(1,i)= tmp(i)*mult
      end do
c     derivatives of the 5th permutationally invaraint polynomial
      mult= 1d0/12d0    !LEGACY: *3d0*(v(5)**2d0))
      tmp2(1)= vr(3)*vr(4)+ vr(2)*vr(5)
      tmp2(2)= vr(3)*vr(4)+ vr(1)*vr(6)
      tmp2(3)= vr(2)*vr(5)+ vr(1)*vr(6)
      tmp2(4)= tmp2(3)
      tmp2(5)= tmp2(2)
      tmp2(6)= tmp2(1)
      do i=1,6
       dvdvr(5,i)= (tmp(i)*vr(7-i)+ tmp2(i))*mult
      end do
c     derivatives of the 6th permutationally invaraint polynomial
      mult= 1d0/3d0     !LEGACY: *4d0*(v(6)**3d0))
      do i=1,6
       dvdvr(6,i)= vr(7-i)*tmp2(i)*mult
      end do
c     derivatives of the 3rd permutationally invariant polynomial
      mult= 1d0/24d0    !LEGACY: *3d0*(v(3)**2d0))
      tmp2(1) = vr(2)**2d0+vr(3)**2d0+vr(4)**2d0+vr(5)**2d0
      tmp2(2) = vr(1)**2d0+vr(3)**2d0+vr(4)**2d0+vr(6)**2d0
      tmp2(3) = vr(1)**2d0+vr(2)**2d0+vr(5)**2d0+vr(6)**2d0
      tmp2(4) = tmp2(3)
      tmp2(5) = tmp2(2)
      tmp2(6) = tmp2(1)
      do i=1,6
       dvdvr(3,i)= (tmp2(i)+ 2d0*vr(i)*tmp(i))*mult
      end do 

      case (2)
c     A2B2 permutational symmetry.
c     derivatives of the 1st permutationally invaraint polynomial
      mult= 1d0/2d0
      tmp(1)= 0d0
      tmp(2)= vr(4)
      tmp(3)= vr(5)
      tmp(4)= vr(2)
      tmp(5)= vr(3)
      tmp(6)= 0d0
      do i=1,6
       dvdvr(1,i)= tmp(i)*mult
      end do
c     derivatives of the 2nd permutationally invaraint polynomial
      tmp(1)= 0d0
      tmp(2)= vr(3)
      tmp(3)= vr(2)
      tmp(4)= vr(5)
      tmp(5)= vr(4)
      tmp(6)= 0d0
      do i=1,6
       dvdvr(2,i)= tmp(i)*mult
      end do
c     derivatives of the 5th permutationally invaraint polynomial
      tmp(1)= 2d0*v(2)
      tmp(2)= vr(1)*vr(3)
      tmp(3)= vr(1)*vr(2)
      tmp(4)= vr(1)*vr(5)
      tmp(5)= vr(1)*vr(4)
      tmp(6)= 0d0
      do i=1,6
       dvdvr(5,i)= tmp(i)*mult
      end do
c     derivatives of the 6th permutationally invaraint polynomial
      tmp(1)= 0d0
      tmp(2)= vr(6)*vr(4)
      tmp(3)= vr(6)*vr(5)
      tmp(4)= vr(6)*vr(2)
      tmp(5)= vr(6)*vr(3)
      tmp(6)= 2d0*v(1)
      do i=1,6
       dvdvr(6,i)= tmp(i)*mult
      end do
c     derivatives of the 3rd permutationally invaraint polynomial
      mult= 1d0/4d0
      tmp(1)= vr(2)+vr(3)+vr(4)+vr(5)
      tmp(2)= vr(1)
      tmp(3)= vr(1)
      tmp(4)= vr(1)
      tmp(5)= vr(1)
      tmp(6)= 0d0
      do i=1,6
       dvdvr(3,i)= tmp(i)*mult
      end do
c     derivatives of the 4th permutationally invaraint polynomial
      tmp(1)= 0d0
      tmp(2)= vr(6)
      tmp(3)= vr(6)
      tmp(4)= vr(6)
      tmp(5)= vr(6)
      tmp(6)= vr(2)+vr(3)+vr(4)+vr(5)
      do i=1,6
       dvdvr(4,i)= tmp(i)*mult
      end do
      end select
      end subroutine calcinvariantsjacobian

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Given a vector of 'raw internuclear distance' coordinates 'qr'
c     and an array of 'Cartesian' coordinates 'y', calculates the 
c     (6x12) Jacobian matrix dqrdy, defined by:
c        dqrdy(i,b,c)= d( qr(i) )/d( y(b,c) )
c     This is used for computing analytic gradients in the 
c     'evaluatesurface' subroutine.
      subroutine calcdqrdyjacobian(qr,y,dqrdy)
      implicit none
      integer :: i
      real*8 :: val
      real*8, dimension(6), intent(in) :: qr
      real*8, dimension(4,3), intent(in) :: y
      real*8, dimension(6,4,3), intent(out) :: dqrdy
      dqrdy= 0d0
      do i=1,3
c      First row: qr1= r12
       val= y(1,i)-y(2,i)
       dqrdy(1,1,i)= val
       dqrdy(1,2,i)= -val
c      Second row: qr2= r13
       val= y(1,i)-y(3,i)
       dqrdy(2,1,i)= val
       dqrdy(2,3,i)= -val
c      Third row: qr3= r14
       val= y(1,i)-y(4,i)
       dqrdy(3,1,i)= val
       dqrdy(3,4,i)= -val
c      Fourth row: qr4= r23
       val= y(2,i)-y(3,i)
       dqrdy(4,2,i)= val
       dqrdy(4,3,i)= -val
c      Fifth row: qr5= r24
       val= y(2,i)-y(4,i)
       dqrdy(5,2,i)= val
       dqrdy(5,4,i)= -val
c      Sixth row: qr6= r34
       val= y(3,i)-y(4,i)
       dqrdy(6,3,i)= val
       dqrdy(6,4,i)= -val
      end do
c     Each row is now scaled by the corresponding distance qr
      do i=1,6
       dqrdy(i,:,:)= dqrdy(i,:,:)/qr(i)
      end do
      end subroutine calcdqrdyjacobian

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Subroutine to evaluate the 2-body potential energy and gradient
c     for given r with generalized Morse potential with Scheme 2.
c      V(r) = De*(1-exp(-f(y)*(r-re)))^2
c      f(y) = c0 + c1*y + c2*y^2 + c3*y^3 + c4*y^4 + c5*y^5 + c6*y^6
c      y = (r^4 - re^4)/(r^4 + re^4)
c      
c      imol .eq. 1   Paramemters for N2 dissociation without SEC 
c                    (for Yuliya's N4 data)
c      imol .eq. 2   Paramemters for N2 dissociation with SEC 
c                    (for Zoltan's N2O2 data)
c      imol .eq. 2   Paramemters for NO dissociation with SEC 
c                    (for Zoltan's N2O2 data)
c      imol .eq. 2   Paramemters for O2 dissociation with SEC 
c                    (for Zoltan's N2O2 data and future O4 data)
c     Units: r is in A, v2 is in kcal/mol, and dv2dr is in 
c     kcal*(mol^-1)*(A^-1).
c     For all of the potentials, the zero energy is defined as the
c     energy of the fully dissociated diatom.
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcpairwise(r,v2,dv2dr,imol,igrad) 
      implicit none 
      integer, intent(in) :: igrad,imol
      real*8, dimension(0:10) :: c
      real*8, intent(in) :: r
      real*8, intent(out) :: v2,dv2dr
      real*8 :: re,de,fy,dfdy,dfdr,u,y,dydr,r2,r4,re2,re4,y2,y4,v

      select case (imol)
      case (1)
c     Parameter for N2 dissociation without SEC (for Yuliya's N4 data)
        re=1.098d0
        de=228.7d0
        c(0) = 2.70963254293d0
        c(1) = 1.32620177271d-1
        c(2) = 2.96757048793d-1
        c(3) = 1.97112432229d-1
        c(4) =-5.02002309588d-1
        c(5) = 3.80734244606d-1
        c(6) = 1.21001628750d0
      case (2)
c     Parameter for N2 dissociation with SEC (for Zoltan's N2O2)
        re=1.098d0
        de=228.4d0
        c(0) = 2.71405774451d0
        c(1) = 1.32757649829d-1
        c(2) = 2.66756890408d-1
        c(3) = 1.95350725241d-1
        c(4) =-4.08663480982d-1
        c(5) = 3.92451705557d-1
        c(6) = 1.13006674877d0
      case (3)
c     Parameter for NO dissociation with SEC (for Zoltan's N2O2)
        re=1.1508d0
        de=149.9d0
        c(0) = 2.81134495569d0
        c(1) = 1.43241169611d-1
        c(2) = 1.35760038863d-2
        c(3) = 3.92892178507d-1
        c(4) = 9.29495534058d-1
        c(5) = 2.66966672332d-1
        c(6) =-3.68118714223d-1
      case (4)
c     Parameter for O2 dissociation with SEC (for Zoltan's N2O2 and KRY's O4)
        re=1.208d0
        de=120.243d0
        c(0) = 2.69132890094d0
        c(1) = 3.39550045614d-1
        c(2) = 3.46401777195d-1
        c(3) =-7.76983671636d-1
        c(4) =-3.29632972405d-1
        c(5) = 2.40883331247d0
        c(6) = 2.09264029009d0
      end select

      r2 = r*r
      r4 = r2*r2
      re2= re*re
      re4= re2*re2
      v  = r4+ re4
      y  = (r4 - re4)/v
      y2 = y*y
      y4 = y2*y2

      fy = c(0) + c(1)*y + c(2)*y2 + c(3)*y2*y + c(4)*y4
     &   + c(5)*y4*y + c(6)*y4*y2

      u  = exp(-fy*(r-re))

      v2 = de*(1.0d0-u)*(1.0d0-u)- de

      if (igrad .eq. 1) then
        dfdy = c(1)+ 2d0*c(2)*y+ 3d0*c(3)*y2+ 4d0*c(4)*y2*y+
     &   5d0*c(5)*y4 + 6d0*c(6)*y4*y
        dydr = 8.0d0*r2*r*re4/(v*v)
        dfdr = dfdy*dydr
        dv2dr = 2.0d0*de*(1d0-u)*u*(dfdr*(r-re)+fy)
      end if
      end subroutine calcpairwise 

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Calculates the sum of the six pairwise interaction terms,
c     one for each of the six internuclear distances.  This 
c     subroutine calls the 'calcpairwise' subroutine six times.
c     - 'rvec' is the array of internuclear distances, in A.
c     - 'v2' is the sum of the six pairwise interaction energies,
c       in kcal/mol
c     - 'dv2dr' is the derivative of 'v2' with respect to each of
c       the six internuclear distances, in kcal*(mol^-1)*(A^-1)
c     - 'systype' is a switch corresponding to the system. Use
c       systype==1 for N4, and systype==2 for N2O2.
      subroutine calcsumpairwise(rvec,v2,dv2dr,systype,igrad)
      implicit none
      integer, intent(in) :: systype,igrad
      integer :: a
      real*8 :: v2,tmp,de
      real*8, dimension(6), intent(in) :: rvec
      real*8, dimension(6), intent(out) :: dv2dr
      v2= 0d0
      do a=1,6
       dv2dr(a)= 0d0
      end do
      select case (systype)
      case (1)
c      N4 system
       do a=1,6
        call calcpairwise(rvec(a),tmp,dv2dr(a),1,igrad)
        v2= v2+ tmp
       end do
c      correct the zero of energy
       de=228.7d0
       v2= v2+ 2d0*de
      case (2)
c      N2O2 system
c      rvec(1) corresponds to the O-O bond distance.
       call calcpairwise(rvec(1),tmp,dv2dr(1),4,igrad)
       v2= v2+ tmp
c      rvec(2:5) correspond to the four N-O bond distance.
       do a=2,5
        call calcpairwise(rvec(a),tmp,dv2dr(a),3,igrad)
        v2= v2+ tmp        
       end do
c      rvec(6) corresponds to the N-N bond distance.
       call calcpairwise(rvec(6),tmp,dv2dr(6),2,igrad)
       v2= v2+ tmp
c      correct the zero of energy; add the dissociation energy for
c      N2 and the dissociation energy for O2
       v2= v2+ 228.4d0+ 120.243d0
      end select
      end subroutine calcsumpairwise
