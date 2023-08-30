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

      module falcons_vars_eval
      implicit none
      real*8, dimension(:,:), allocatable :: qi,aji
      end module falcons_vars_eval

      module falcons_vars_main
      implicit none
      integer :: i,j,k,l,nn,ind,n,m,eta,np,moltype,pesform,&
     & nma,nplmax,nall,nmincrc,mcrc,qps,pwonly,cgmeth,&
     & ntarmin,ntarmax,ntar,crcmult,nitcrc
      integer, dimension(6) :: sd
      integer, dimension(4,3) :: cmap
      integer, dimension(24,6) :: qp
      real*8 :: pnu,epsnu,ps,pu,epsu,e,epw,systo,r,rwfa,rwfb
      real*8, dimension(6) :: t,qr,qs,q,xr,x,aa,qeq
      real*8, dimension(4,3) :: y,dedy
      real*8, dimension(:), allocatable :: crc
      end module falcons_vars_main

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
      file_path1 = trim(path)//"/N4/pes.dat"
      inquire(file=file_path1,exist=ex)
      if (.not. ex) then
       write(*,'(a)') 'ERROR: No pes.dat file found.'
       stop
      end if
      open(fid,file=file_path1,action='read')
      read(fid,*)
      read(fid,*)
      read(fid,*)
      read(fid,*)
      read(fid,*)
      do i=1,20
       read(fid,*)
      end do
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
      nall= n+nma
      allocate(qi(n,6),aji(n,m))
      read(fid,*)
      read(fid,*)
      j= 0
      if (pesform== 1) then
       do i=1,nall
        read(fid,*) oridval,e,(t(l),l=1,6),maskval
        if (maskval== 0) then
         j= j+1
         call default2cartesian(&
     &    t(1),t(2),t(3),t(4),t(5),t(6),y)
         call cartesian2rawdist(y,qr)
         call rawdist2rawmorse(qr,qi(j,:),qeq,aa)
        end if
       end do
      else
       write(*,'(a)') &
     &  'ERROR: pesform /= 1 is not supported in this code.'
       stop
      end if
      if (j /= n) then
       write(*,'(a)') 'ERROR: Bad pes.dat file.'//&
     & ' Inconsistent number of unmasked points.'
       stop
      end if
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
       deallocate(crc)
      end if
      allocate(crc(mcrc))
      read(fid,*)
      read(fid,*) (crc(l),l=1,mcrc)
      read(fid,*)
      read(fid,*) ntarmin,ntarmax,ntar,crcmult,nitcrc,rwfa,rwfb
      end subroutine readcrcparameters

      subroutine pot(x_,y_,z_,e_,dedx_,dedy_,dedz_)
      use falcons_vars_eval
      use falcons_vars_main
      implicit none
      real*8, dimension(4), intent(in) :: x_,y_,z_
      real*8 :: cconv,econv,gconv,eref
      real*8, dimension(4) :: x__,y__,z__
      real*8, intent(out) :: e_
      real*8, dimension(4), intent(out) :: dedx_,dedy_,dedz_

      pu= 3d0
      epsu= 1d-1
      pwonly= 0
      cgmeth= 1 

      cconv= 0.52917721092d0
      econv= 0.159360144d-2
      gconv= 0.843297564d-3
      eref= -218.4080132d0
      x__= x_*cconv
      y__= y_*cconv
      z__= z_*cconv
      do l=1,4
       y(l,1)= x__(l)
       y(l,2)= y__(l)
       y(l,3)= z__(l)
      end do
      call calccoordinates(y,qr,xr,q,x,sd,aa,qeq,moltype,qp,qps)
      call evaluatesurface(cgmeth,y,qr,xr,q,x,moltype,pwonly,&
     & cmap,qp,qps,mcrc,crc,n,m,eta,ps,pu,epsu,aa,r,np,e,dedy) 
      e_= e
      dedx_= dedy(:,1)
      dedy_= dedy(:,2)
      dedz_= dedy(:,3)
      e_= e_*econv+ eref
      dedx_= dedx_*gconv
      dedy_= dedy_*gconv
      dedz_= dedz_*gconv
      end subroutine pot

      subroutine evaluatesurface(cg,y,qr,xr,q,x,moltype,pwonly,&
     & cmap,qp,qps,mcrc,crc,n,m,eta,ps,pu,epsu,aa,r,np,e,dedy)
      use falcons_vars_eval
      implicit none
      integer, intent(in) :: cg,moltype,pwonly,qps,mcrc,n,m,eta
      integer, dimension(4,3), intent(in) :: cmap
      integer, dimension(qps,6), intent(in) :: qp
      real*8, intent(in) :: ps,pu,epsu
      real*8, dimension(6), intent(in) :: qr,xr,q,x,aa
      real*8, dimension(mcrc), intent(in) :: crc
      real*8, dimension(4,3), intent(in) :: y
      integer, intent(out) :: np
      real*8, intent(out) :: r,e
      real*8, dimension(4,3), intent(out) :: dedy
      integer :: i,j,k,l,a,b,c
      real*8 :: rs,drdqch,dis,wsum,wi,pi,dwidz,tmp,epw
      real*8, dimension(6) :: qii,dwidxr,dwdxrsum,dpidx,dpidxr,&
     & dedxrpart,dedxr,dedqr,ddisdq,depwdqr
      real*8, dimension(m) :: aj,bj 
      real*8, dimension(6,6) :: dxdxr
      real*8, dimension(m,6) :: dbjdx
      real*8, dimension(6,4,3) :: dqrdy
      np=   0 
      e=    0d0
      wsum= 0d0 !sum of the weight functions
      do a=1,6
       dedxrpart(a)= 0d0 
       dwdxrsum(a)=  0d0 !sums of the weight function derivatives
      end do
      call calcqch(q,tmp)
      call calccr(tmp,mcrc,crc,r,rs,drdqch)
      call calcbj(bj,dbjdx,x,eta,m,1)
      if (cg /= 0) then
       call calcinvariantsjacobian(xr,x,dxdxr,moltype)
       call calcdqrdyjacobian(qr,y,dqrdy)
      end if

      call calcsumpairwise(qr,epw,depwdqr,moltype,cg)

      if (pwonly == 0) then
      do i=1,n

       do a=1,6
        qii(a)= qi(i,a)
       end do
       call distance6ds_perm(q,qii,qp,qps,0,dis,ddisdq)

       if (dis <= rs) then
        np= np+1  
        if (cg /= 0) then
         call distance6ds_perm(q,qii,qp,qps,cg,dis,ddisdq)
        end if
        do j=1,m
         aj(j)= aji(i,j)
        end do
        call weightw(wi,dwidz,dis,rs,ps,pu,epsu)
        pi= 0d0
        do j=1,m
         pi= pi+ aj(j)*bj(j)
        end do
        wsum= wsum+ wi
        e= e+ wi*pi
        if (cg /= 0) then
         do a=1,6
          dpidx(a)= 0d0
         end do
         do j=1,m
          do a=1,6
           dpidx(a)= dpidx(a)+ aj(j)*dbjdx(j,a)
          end do
         end do
         do a=1,6
          tmp= 0d0
          do k=1,6
           tmp= tmp+ dpidx(k)*dxdxr(k,a)
          end do
          dpidxr(a)= tmp
         end do
         tmp= dis*4d0*drdqch/r
         do a=1,6
          dwidxr(a)= dwidz*(ddisdq(a)- tmp*q(a))/rs
          dwdxrsum(a)= dwdxrsum(a)+ dwidxr(a)
         end do
         do a=1,6
          dedxrpart(a)= dedxrpart(a)+ wi*dpidxr(a)+ pi*dwidxr(a) 
         end do
        end if
       end if
      end do  !end loop over data points within the cutoff radius
      end if  !end if 'pwonly==0' computation

      if (wsum== 0d0) then
       wsum= 1d0
      end if

      e= e/wsum

      if (cg /= 0) then
       do a=1,6
        dedxr(a)= (dedxrpart(a)- e*dwdxrsum(a))/wsum
       end do
       do a=1,6
        dedqr(a)= -dedxr(a)*xr(a)/aa(a)+ depwdqr(a)
       end do
       do b=1,4 
        do c=1,3
         tmp= 0d0
         do k=1,3
          a= cmap(b,k)
          tmp= tmp+ dedqr(a)*dqrdy(a,b,c)
         end do
         dedy(b,c)= tmp
        end do
       end do
      else
       do b=1,4 
        do c=1,3
         dedy(b,c)= 0d0
        end do
       end do
      end if

      e= e+ epw

      if (np==0) then
       write(*,'(a)') &
     &  ' WARNING: No points found within cutoff radius for'//&
     &  ' evaluation point:'
       write(*,'(a,6(1x,f15.5))') ' qr= ',(qr(l),l=1,6)
       write(*,'(a,6(1x,f15.5))') ' q = ',(q(l),l=1,6)
       write(*,'(a,6(1x,f15.5))') ' x = ',(x(l),l=1,6)
      end if
      end subroutine evaluatesurface
  
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
       dcrdqch= dcrdqch+ (crc(icrc)*dble(icrc-1)*pqch)
       pqch= pqch*qch 
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
      ind= 0
      do j=1,6
       ind= ind+1
       b(ind)= x(j) !linear terms
      end do
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
        end do
       end do
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
        end do
       end do
       do j= 1,6
        do k= j+1,6
         tmp= xy(j,k)
         do l= k+1,6
          ind= ind+1
          tmp2= tmp*x(l)
          b(ind)= tmp2
          xyz(j,k,l)= tmp2
         end do
        end do
       end do
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
        end do
       end do
       do j= 1,6
        tmp= xx(j)
        do k= j+1,6
         tmp2= xx(k)
         ind= ind+1
         b(ind)= tmp*tmp2
         tmp3= xy(j,k)
         do l= k+1,6
          ind= ind+1
          b(ind)= tmp*xy(k,l)
          ind= ind+1
          b(ind)= tmp2*xy(j,l)
          ind= ind+1
          b(ind)= xx(l)*tmp3
         end do
        end do
       end do
       do j=1,6
        do k= j+1,6
         tmp= xy(j,k)
         do l= k+1,6
          do ll= l+1,6
           ind= ind+1
           b(ind)= tmp*xy(l,ll)
          end do
         end do
        end do
       end do
      end if
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
       if (eta .ge. 2) then
        do j= 1,6
         ind= ind+1
         if (j == a) then
          dbdxa(ind)= x2(j)
         end if
         do k= j+1,6
          ind= ind+1
          if (j == a) then
           dbdxa(ind)= x(k)
          else if (k == a) then
           dbdxa(ind)= x(j)
          end if
         end do
        end do
       end if
       if (eta .ge. 3) then
        do j= 1,6
         ind= ind+1
         if (j == a) then
          dbdxa(ind)= 3d0*xx(a)
          do k= j+1,6
           ind= ind+1
           dbdxa(ind)= xy2(j,k)
           ind= ind+1
           dbdxa(ind)= xx(k)
          end do
         else 
          do k= j+1,6
           if (k == a) then
            ind= ind+1
            dbdxa(ind)= xx(j)
            ind= ind+1
            dbdxa(ind)= xy2(j,k)
           else
            ind= ind+2
           end if
          end do
         end if
        end do
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
          end do
         end do
        end do
       end if
       if (eta .ge. 4) then
        do j= 1,6
         ind= ind+1
         if (j == a) then
          dbdxa(ind)= 4d0*xxx(a)
          do k= j+1,6
           ind= ind+1
           dbdxa(ind)= 3d0*xx(j)*x(k)
           ind= ind+1
           dbdxa(ind)= xxx(k)
          end do
         else 
          do k= j+1,6
           if (k == a) then
            ind= ind+1
            dbdxa(ind)= xxx(j)
            ind= ind+1
            dbdxa(ind)= 3d0*x(j)*xx(k)
           else
            ind= ind+2
           end if
          end do
         end if
        end do
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
          do l= k+1,6
           ind= ind+1
           if (j == a) then
            dbdxa(ind)= x2(j)*xy(k,l)
           else if (k == a) then
            dbdxa(ind)= tmp*x(l)
           else if (l == a) then
            dbdxa(ind)= tmp*x(k)
           end if
           ind= ind+1
           if (j == a) then
            dbdxa(ind)= tmp2*x(l)
           else if (k == a) then
            dbdxa(ind)= x2(k)*xy(j,l)
           else if (l == a) then
            dbdxa(ind)= tmp2*x(j)
           end if
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
           end do
          end do
         end do
        end do
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
      iszero= 0 
      dssum= 0d0
      do a=1,6
       ddsdq(a)= 0d0
      end do
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
        iszero= 1
       end if
      end do
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
      dwdz= uval*4d0*((1d0-z**ps)**3d0)*(-ps*(z**(ps-1d0)))+&
     & sval*(uval*uval)*(-(pu*(z**(pu-1d0))))
      end subroutine weightw

      subroutine calccoordinates(&
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

      subroutine default2cartesian(ra,rb,d,thetaa_d,thetab_d,&
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

      subroutine distance3d(p1,p2,d)
      implicit none
      real*8, intent(out) :: d
      real*8, dimension(3), intent(in) :: p1,p2
      d=&
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
       v(1)= (vr(1)*vr(2)+vr(1)*vr(3)+vr(1)*vr(4)+&
     &        vr(1)*vr(5)+vr(2)*vr(3)+vr(2)*vr(4)+&
     &        vr(2)*vr(6)+vr(3)*vr(5)+vr(3)*vr(6)+&
     &        vr(4)*vr(5)+vr(4)*vr(6)+vr(5)*vr(6))/12d0
       v(1)= v(1)       !LEGACY: **0.5d0
       v(2)= (vr(1)*vr(2)*vr(4)+ vr(1)*vr(3)*vr(5)+&
     &        vr(2)*vr(3)*vr(6)+ vr(4)*vr(5)*vr(6))/4d0
       v(2)= v(2)       !LEGACY: **(1d0/3d0)
       v(3)= (vr(1)*vr(2)*(vr(1)+vr(2))+&
     &        vr(1)*vr(3)*(vr(1)+vr(3))+&
     &        vr(1)*vr(4)*(vr(1)+vr(4))+&
     &        vr(1)*vr(5)*(vr(1)+vr(5))+&
     &        vr(2)*vr(3)*(vr(2)+vr(3))+&
     &        vr(2)*vr(4)*(vr(2)+vr(4))+&
     &        vr(2)*vr(6)*(vr(2)+vr(6))+&
     &        vr(3)*vr(5)*(vr(3)+vr(5))+&
     &        vr(3)*vr(6)*(vr(3)+vr(6))+&
     &        vr(4)*vr(5)*(vr(4)+vr(5))+&
     &        vr(4)*vr(6)*(vr(4)+vr(6))+&
     &        vr(5)*vr(6)*(vr(5)+vr(6)))/24d0
       v(3)= v(3)       !LEGACY: **(1d0/3d0)
       v(4)= (vr(1)*vr(2)*vr(3)+ vr(1)*vr(4)*vr(5)+&
     &        vr(2)*vr(4)*vr(6)+ vr(3)*vr(5)*vr(6))/4d0
       v(4)= v(4)       !LEGACY: **(1d0/3d0)
       v(5)= (vr(1)*vr(3)*vr(4)+ vr(2)*vr(3)*vr(4)+&
     &        vr(1)*vr(2)*vr(5)+ vr(2)*vr(3)*vr(5)+&
     &        vr(2)*vr(4)*vr(5)+ vr(3)*vr(4)*vr(5)+&
     &        vr(1)*vr(2)*vr(6)+ vr(1)*vr(3)*vr(6)+&
     &        vr(1)*vr(4)*vr(6)+ vr(3)*vr(4)*vr(6)+&
     &        vr(1)*vr(5)*vr(6)+ vr(2)*vr(5)*vr(6))/12d0
       v(5)= v(5)       !LEGACY: **(1d0/3d0)
       v(6)= (vr(1)*vr(2)*vr(5)*vr(6)+&
     &        vr(1)*vr(3)*vr(4)*vr(6)+&
     &        vr(2)*vr(3)*vr(4)*vr(5))/3d0
       v(6)= v(6)       !LEGACY: **0.25d0

      case (2)
       v(1)= (vr(2)*vr(4)+ vr(3)*vr(5))/2d0
       v(2)= (vr(2)*vr(3)+ vr(4)*vr(5))/2d0
       v(3)= (vr(1)*(vr(2)+vr(3)+vr(4)+vr(5)))/4d0
       v(4)= (vr(6)*(vr(2)+vr(3)+vr(4)+vr(5)))/4d0
       v(5)= (vr(1)*(vr(2)*vr(3)+vr(4)*vr(5)))/2d0
       v(6)= (vr(6)*(vr(2)*vr(4)+vr(3)*vr(5)))/2d0
      end select
      end subroutine calcinvariants

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
      mult= 1d0/4d0     !LEGACY: *3d0*(v(2)**2d0))
      do i=1,6
       dvdvr(2,i)= tmp(7-i)*mult
      end do
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
      mult= 1d0/3d0     !LEGACY: *4d0*(v(6)**3d0))
      do i=1,6
       dvdvr(6,i)= vr(7-i)*tmp2(i)*mult
      end do
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
      tmp(1)= 0d0
      tmp(2)= vr(3)
      tmp(3)= vr(2)
      tmp(4)= vr(5)
      tmp(5)= vr(4)
      tmp(6)= 0d0
      do i=1,6
       dvdvr(2,i)= tmp(i)*mult
      end do
      tmp(1)= 2d0*v(2)
      tmp(2)= vr(1)*vr(3)
      tmp(3)= vr(1)*vr(2)
      tmp(4)= vr(1)*vr(5)
      tmp(5)= vr(1)*vr(4)
      tmp(6)= 0d0
      do i=1,6
       dvdvr(5,i)= tmp(i)*mult
      end do
      tmp(1)= 0d0
      tmp(2)= vr(6)*vr(4)
      tmp(3)= vr(6)*vr(5)
      tmp(4)= vr(6)*vr(2)
      tmp(5)= vr(6)*vr(3)
      tmp(6)= 2d0*v(1)
      do i=1,6
       dvdvr(6,i)= tmp(i)*mult
      end do
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

      subroutine calcdqrdyjacobian(qr,y,dqrdy)
      implicit none
      integer :: i
      real*8 :: val
      real*8, dimension(6), intent(in) :: qr
      real*8, dimension(4,3), intent(in) :: y
      real*8, dimension(6,4,3), intent(out) :: dqrdy
      dqrdy= 0d0
      do i=1,3
       val= y(1,i)-y(2,i)
       dqrdy(1,1,i)= val
       dqrdy(1,2,i)= -val
       val= y(1,i)-y(3,i)
       dqrdy(2,1,i)= val
       dqrdy(2,3,i)= -val
       val= y(1,i)-y(4,i)
       dqrdy(3,1,i)= val
       dqrdy(3,4,i)= -val
       val= y(2,i)-y(3,i)
       dqrdy(4,2,i)= val
       dqrdy(4,3,i)= -val
       val= y(2,i)-y(4,i)
       dqrdy(5,2,i)= val
       dqrdy(5,4,i)= -val
       val= y(3,i)-y(4,i)
       dqrdy(6,3,i)= val
       dqrdy(6,4,i)= -val
      end do
      do i=1,6
       dqrdy(i,:,:)= dqrdy(i,:,:)/qr(i)
      end do
      end subroutine calcdqrdyjacobian

      subroutine calcpairwise(r,v2,dv2dr,imol,igrad) 
      implicit none 
      integer, intent(in) :: igrad,imol
      real*8, dimension(0:10) :: c
      real*8, intent(in) :: r
      real*8, intent(out) :: v2,dv2dr
      real*8 :: re,de,fy,dfdy,dfdr,u,y,dydr,r2,r4,re2,re4,y2,y4,v

      select case (imol)
      case (1)
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

      fy = c(0) + c(1)*y + c(2)*y2 + c(3)*y2*y + c(4)*y4&
     &   + c(5)*y4*y + c(6)*y4*y2

      u  = exp(-fy*(r-re))

      v2 = de*(1.0d0-u)*(1.0d0-u)- de

      if (igrad .eq. 1) then
        dfdy = c(1)+ 2d0*c(2)*y+ 3d0*c(3)*y2+ 4d0*c(4)*y2*y+&
     &   5d0*c(5)*y4 + 6d0*c(6)*y4*y
        dydr = 8.0d0*r2*r*re4/(v*v)
        dfdr = dfdy*dydr
        dv2dr = 2.0d0*de*(1d0-u)*u*(dfdr*(r-re)+fy)
      end if
      end subroutine calcpairwise 

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
       do a=1,6
        call calcpairwise(rvec(a),tmp,dv2dr(a),1,igrad)
        v2= v2+ tmp
       end do
       de=228.7d0
       v2= v2+ 2d0*de
      case (2)
       call calcpairwise(rvec(1),tmp,dv2dr(1),4,igrad)
       v2= v2+ tmp
       do a=2,5
        call calcpairwise(rvec(a),tmp,dv2dr(a),3,igrad)
        v2= v2+ tmp        
       end do
       call calcpairwise(rvec(6),tmp,dv2dr(6),2,igrad)
       v2= v2+ tmp
       v2= v2+ 228.4d0+ 120.243d0
      end select
      end subroutine calcsumpairwise
