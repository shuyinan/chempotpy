!=======================================================================
!     The ^2A doublet prime PES for C2N system at the MC-PDFT/AVTZ level.
!=======================================================================
      subroutine pes(x,igrad,path,p,g,d)
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: v
      double precision :: rb(3)
      integer :: iatom, idir, j, istate
      !initialize 
      v=0.d0
      p=0.d0
      g=0.d0
      d=0.d0

      ! order of r_CC,r_NC,r_NC
      rb=0.d0
      rb(1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2+&
               &(x(1,3)-x(2,3))**2)
      rb(2)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2+&
               &(x(1,3)-x(3,3))**2)
      rb(3)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2+&
               &(x(2,3)-x(3,3))**2)

      if (igrad==0) then
        call c2n_2a_dprime(rb, path, v)
      else
        write (*,*) 'Only energy is available'
      endif

      do istate=1,nstates
        p(istate)=v
      enddo

      g=0.d0
      d=0.d0

      endsubroutine


!=======================================================================
!     PIP-NN PES of total energy for the C2N(^2A") system.
!     where,
!       rNC,rCC: angstrom
!       theta: degree
!       Vpes: total energy, eV
!-----------------------------------------------------------------------
      subroutine c2n_2a_dprime(rb,path,Vpes)
      implicit none
      integer,parameter :: npip=12
      real*8,parameter :: alpha=1.5d0 ! angstrom
      real*8,parameter :: Vpescut=5.d5 ! eV
      real*8,parameter :: pi=dacos(-1.d0)
      integer :: i
      character(len=1024), intent(in) :: path
      real*8 :: Vpesa(1),Vpesb(1),Vpesc(1),Vpes
      real*8 :: rb(3),xbond(3),basis(0:npip),txinput(npip)

      if (minval(rb).lt.0.79999999d0) then
         Vpes=Vpescut
         Vpesa=Vpescut
         Vpesb=Vpescut
         Vpesc=Vpescut
      else
         xbond(:)=dexp(-rb(:)/alpha)

         basis=0.d0
         call bemsav(xbond,basis)

         do i=1,npip
            txinput(i)=basis(i)
         enddo

         call nnfit_a(npip,1,txinput,path,Vpesa)
         call nnfit_b(npip,1,txinput,path,Vpesb)
         call nnfit_c(npip,1,txinput,path,Vpesc)
         Vpes=(Vpesa(1)+Vpesb(1)+Vpesc(1))/3.0d0

         if (Vpes.gt.Vpescut) then
            Vpes=Vpescut
            Vpesa=Vpescut
            Vpesb=Vpescut
            Vpesc=Vpescut
         endif
      endif

      return
      end subroutine c2n_2a_dprime
!=======================================================================


!=======================================================================
      subroutine nnfit_a(ndim,ntot,r0,path,vx)
      implicit none
      integer,intent(in) :: ndim,ntot
      real*8,intent(in) :: r0(ndim,ntot)
      real*8,intent(out) :: vx(ntot)

      integer,parameter :: s0=12,s1=30,s2=30
      real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3
      real*8,save :: rg(2,s0),vg(2)
      character(len=1024), intent(in) :: path
      character(len=1024) :: file_path1

      real*8 :: dv(ndim,ntot)
      integer :: i
      integer,save :: init=0
      integer,save :: fid

      if (ndim.ne.s0) stop "ndim .ne. s0"

      file_path1 = trim(path)//"/C2N/c2n_2a_dprime_para_a.txt"

      if (init.eq.0) then
         init=1
         fid=123456
         open(fid,file=file_path1,status='old',action='read')
         read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
         close(fid)
      endif

      if (ntot.ge.24) then
         call nsimx(ntot,r0,vx,0,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
      else
         do i=1,ntot
            call nsim(r0(1,i),vx(i),0,dv(:,i),s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
         enddo
      endif

      return
      end subroutine
!=======================================================================


!=======================================================================
      subroutine nnfit_b(ndim,ntot,r0,path,vx)
      implicit none
      integer,intent(in) :: ndim,ntot
      real*8,intent(in) :: r0(ndim,ntot)
      real*8,intent(out) :: vx(ntot)

      integer,parameter :: s0=12,s1=30,s2=30         
      real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3
      real*8,save :: rg(2,s0),vg(2)
      character(len=1024), intent(in) :: path
      character(len=1024) :: file_path1

      real*8 :: dv(ndim,ntot)
      integer :: i 
      integer,save :: init=0
      integer,save :: fid

      if (ndim.ne.s0) stop "ndim .ne. s0"

      file_path1 = trim(path)//"/C2N/c2n_2a_dprime_para_b.txt"

      if (init.eq.0) then
         init=1
         fid=123456
         open(fid,file=file_path1,status='old',action='read')
         read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
         close(fid)
      endif

      if (ntot.ge.24) then
         call nsimx(ntot,r0,vx,0,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
      else
         do i=1,ntot
            call nsim(r0(1,i),vx(i),0,dv(:,i),s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
         enddo
      endif
      
      return
      end subroutine
!=======================================================================


!=======================================================================
      subroutine nnfit_c(ndim,ntot,r0,path,vx)
      implicit none
      integer,intent(in) :: ndim,ntot
      real*8,intent(in) :: r0(ndim,ntot)
      real*8,intent(out) :: vx(ntot)

      integer,parameter :: s0=12,s1=30,s2=30         
      real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3
      real*8,save :: rg(2,s0),vg(2)
      character(len=1024), intent(in) :: path
      character(len=1024) :: file_path1

      real*8 :: dv(ndim,ntot)
      integer :: i 
      integer,save :: init=0
      integer,save :: fid

      if (ndim.ne.s0) stop "ndim .ne. s0"

      file_path1 = trim(path)//"/C2N/c2n_2a_dprime_para_c.txt"

      if (init.eq.0) then
         init=1
         fid=123456
         open(fid,file=file_path1,status='old',action='read')
         read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
         close(fid)
      endif

      if (ntot.ge.24) then
         call nsimx(ntot,r0,vx,0,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
      else
         do i=1,ntot
            call nsim(r0(1,i),vx(i),0,dv(:,i),s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
         enddo
      endif
      
      return
      end subroutine
!=======================================================================


!=======================================================================
! Simulate a neural network with two hidden layers: n0-n1-n2-1
! Blas routines used in this subroutine: dgemv, ddot
!-----------------------------------------------------------------------
      subroutine nsim(r0,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
      implicit none
      integer,intent(in) :: n0,n1,n2,idv
      real*8,intent(in) :: r0(n0),rg(2,n0),vg(2)
      real*8,intent(in) :: w1(n0,n1),b1(n1)
      real*8,intent(in) :: w2(n1,n2),b2(n2)
      real*8,intent(in) :: w3(n2),b3
      real*8,intent(out) :: v,dv(n0)

      integer :: i,j,k
      real*8 :: r(n0),rgg(n0),vgg,ax(n1),bx(n2)
      real*8 :: dvtm,rtmp,rt1(n1),rt2(n2)
      real*8,external :: ddot

      v=0.d0
      r=r0
      vgg=vg(2)-vg(1)
      ! mapminmax [-1,1]
      do i=1,n0
        rgg(i)=rg(2,i)-rg(1,i)
        r(i)=2.d0*(r(i)-rg(1,i))/rgg(i)-1.d0
      end do
      ! 1st layer
      rt1=b1
      call dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1)
      ax=dtanh(rt1)
      ! 2nd layer
      rt2=b2
      call dgemv('t',n1,n2,1.d0,w2,n1,ax,1,1.d0,rt2,1)
      bx=dtanh(rt2)
      ! output layer
      v=b3+ddot(n2,w3,1,bx,1)
      !reverse map
      v=vgg*(v+1.d0)/2+vg(1)
      if(idv.ne.1) return

      ! calculate first derivatives, dv(i)=dv/dr(i)
      dv=0.d0
      do i=1,n0
        do k=1,n2
          dvtm=0.d0
          do j=1,n1
            dvtm=dvtm+w2(j,k)*w1(i,j)*(1-ax(j)**2)
          enddo
          dv(i)=dv(i)+w3(k)*dvtm*(1-bx(k)**2)
        enddo
        dv(i)=dv(i)*vgg/rgg(i)
      enddo

      return
      end subroutine nsim
!=======================================================================


!=======================================================================
! Simulate a neural network with two hidden layers:
!    n0-n1(tansig)-n2(tansig)-1(pureline)
! Blas routines used in this subroutine: dgemm, dgemv
! The syntax keeps the same with nsim, except input nt at first.
! Chen Jun, 2013/11/26
!-----------------------------------------------------------------------
! Description of parameters:
!    nt  : total number of points in r,v and dv
!    r   : input of nn, r(n0,nt), r(1:n0,i) is the i-th input
!    v   : output of nn, v(nt), v(i) is the i-th input
!    idv : if idv==0, calc. v only, else calc. v and derivatives dv
!    dv  : derivatives of dv/dr, dv(n0,nt), dv(j,i)=dv(i)/dr(j,i)
!    n0  : neurons of the input layer, dimension of r
!    n1  : neurons of the first hidden layer
!    n2  : neurons of the second hidden layer
!    rg  : input ranges of the training set, rg(2,n0)
!    w1  : weights connect input layer with 1-st hidden layer
!    b1  : biases of 1-st hidden layer
!    w2  : weights connect 1-st hidden layer with 2-nd hidden layer
!    b2  : biases of 2-nd hidden layer
!    w3  : weights connect 2-nd hidden layer with output layer
!    b3  : biase of output layer
!    vg  : output range of the training set, vg(2)
!-----------------------------------------------------------------------
      subroutine nsimx(nt,r,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
      implicit none
      integer,intent(in) :: nt,idv,n0,n1,n2
      real*8,intent(in) :: r(n0,nt)
      real*8,intent(out) :: v(nt),dv(n0,nt)
      real*8,intent(in) :: w1(n0,n1),b1(n1),w2(n1,n2),b2(n2),w3(n2),b3
      real*8,intent(in) :: rg(2,n0),vg(2)

      real*8 :: x0(n0,nt),x1(nt,n1),x2(nt,n2),tmp,rgg(n0),vgg
      integer :: i,j,k,n

      v=0.d0
      do i=1,n0
        rgg(i)=rg(2,i)-rg(1,i)
      enddo
      vgg=vg(2)-vg(1)
      
      x0=r
      
      do i=1,n0
        x0(i,1:nt)=(x0(i,1:nt)-rg(1,i))/rgg(i)*2.d0-1.d0
      enddo
      
      do i=1,n1
       x1(1:nt,i)=b1(i)
      enddo
      
      ! ---- x1=trans(x0(n0,nt))*w1(n0,n1)+x1(nt,n1)
      !call gpu_dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
      call dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
      
      x1=dtanh(x1)
      
      do i=1,n2
       x2(1:nt,i)=b2(i)
      enddo
      
      ! ---- x2=x1(nt,n1)*w2(n1,n2)+x2(nt,n2)
      !call gpu_dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
      call dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
      
      x2=dtanh(x2)
      
      v(1:nt)=b3
      
      ! ---- v=x2(nt,n2)*w3(n2)+v(nt)
      call dgemv('n',nt,n2,1.d0,x2,nt,w3,1,1.d0,v,1)
      
      v=(v+1.d0)/2.d0*vgg+vg(1)
      
      if (idv.ne.1) return
      ! ---- calculate first derivatives
      dv=0.d0
      do n=1,nt
        do i=1,n0
          do k=1,n2
            tmp=0.d0
            do j=1,n1
              tmp=tmp+w2(j,k)*w1(i,j)*(1-x1(n,j)**2)
            enddo
            dv(i,n)=dv(i,n)+w3(k)*tmp*(1-x2(n,k)**2)
          enddo
          dv(i,n)=dv(i,n)*vgg/(rgg(i))
        enddo
      enddo

      return
      end subroutine nsimx
!=======================================================================


!=======================================================================
!------->This is the subroutine for PIPs up to the order of 3 of A2B1
!------->type molecule.
!-----------------------------------------------------------------------
!        double precision :: x(1:3)
!        double precision :: p(0:12)
!        call bemsav(x,p)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
        subroutine bemsav(x,p)
        implicit none
        double precision,intent(in) :: x(1:3)
        double precision :: m(0:4)
        double precision,intent(out) :: p(0:12)

        call evmono(x,m)
        call evpoly(m,p)

        return
        end subroutine bemsav
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
        subroutine evmono(x,m)
        implicit none
        double precision,intent(in) :: x(1:3)
        double precision,intent(out) :: m(0:4)

        m(0) = 1.0D0
        m(1) = x(3)
        m(2) = x(2)
        m(3) = x(1)
        m(4) = m(1)*m(2)

        return
        end subroutine evmono
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
        subroutine evpoly(m,p)
        implicit none
        double precision,intent(in) :: m(0:4)
        double precision,intent(out) :: p(0:12)

        p(0) = m(0)
        p(1) = m(1)+m(2)
        p(2) = m(3)
        p(3) = m(4)
        p(4) = p(2)*p(1)
        p(5) = p(1)*p(1)-p(3)-p(3)
        p(6) = p(2)*p(2)
        p(7) = p(2)*p(3)
        p(8) = p(3)*p(1)
        p(9) = p(2)*p(5)
        p(10) = p(2)*p(4)
        p(11) = p(1)*p(5)-p(8)
        p(12) = p(2)*p(6)

        return
        end subroutine evpoly
!=======================================================================
