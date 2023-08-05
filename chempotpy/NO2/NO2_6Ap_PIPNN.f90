!********************************************************************** 
!   System:                     NO2
!   Functional form:            permutationally invariant polynomials
!   Common name:                NO2(adiabatic sextet ground state)
!   Number of derivatives:      1
!   Number of bodies:           3
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References:: Z. Varga, Y. Liu, J. Li, Y. Paukku, H. Guo, and
!                D. G. Truhlar, in preparation
!
!   Notes:    -Sextet A' (6Ap) surface 
!             -Neural network fit
!             -the code requires files weights.txt and biases.txt
!
!     O1--O2
!
!     N3
!
!********************************************************************** 

!***************************************************************************
!                         PIP-NN PES for NO2(6A')
!***************************************************************************
!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
      module nnparam
      implicit none
!***************************************************************************
!natom     ==>Number of atoms
!npes      ==>Number of PESs
!***************************************************************************
      real*8,parameter::alpha=1.0d0,PI=3.141592653589793238d0,
     %radian=PI/180.0d0,bohr=0.5291772d0
      integer,parameter::natom=3,npes=3
      integer,parameter::nbond=natom*(natom-1)/2
      integer ninput,noutput,nscale
      integer nhid3,nlayer3,ifunc3
      integer nwe3,nodemax3
      integer,allocatable::nodes3a(:)
      real*8,allocatable::weight3a(:,:,:,:),bias3a(:,:,:),
     %pdel3a(:,:),pavg3a(:,:)
      end module nnparam

      subroutine pes(x,igrad,path,p,g,d)
      use nnparam

      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: v, vpes(3)
      double precision :: tx(3,natoms), tg(3*natoms), tgn(3*natoms,3)
      integer :: iatom, idir, j, istate
      !initialize 
      p=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
        do idir=1,3
          tx(idir, iatom)=x(iatom, idir)
        enddo
      enddo
      tg=0.d0

      call pes_init(path)
      call evvdvdx(tx,v,vpes,tg,tgn,igrad)

      deallocate(nodes3a)
      deallocate(weight3a)
      deallocate(bias3a)
      deallocate(pdel3a)
      deallocate(pavg3a)

      do istate=1,nstates
        p(istate)=v
      enddo

      do istate=1,nstates
        do iatom=1,natoms
          do idir=1,3
            j=3*(iatom-1)+idir
            g(istate,iatom,idir)=tg(j)
          enddo
        enddo
      enddo

      d=0.d0

      endsubroutine

!***************************************************************************
      subroutine evvdvdx(xcart,v,vpes,dvdxa,dvdx,ndriv)
!***************************************************************************
!Subroutine to calculate the average potential energy v and gradient dvdxa
!Xcart     ==>Cartesian coordinates in angstrom
!             O1: xcart(1,1), xcart(2,1), xcart(3,1)
!             O2: xcart(1,2), xcart(2,2), xcart(3,2)
!             N3: xcart(1,3), xcart(2,3), xcart(3,3)
!v         ==>Average potential energy(in eV), v=sum(vpes)/npes
!vpes      ==>Potential energy(in eV) for each PES
!dvdxa     ==>Average gradient(in eV/ang), dvdxa=sum(dvdx)/npes
!dvdx      ==>Gradient(in eV/ang) for each PES
!ndriv     ==>ndriv=1 to compute the gradient, otherwise ndriv=0
!PES1: totrmse=0.0455, errmax=0.6321, trainrmse=0.0445, validrmse=0.0503, testrmse=0.0565eV 
!PES2: totrmse=0.0521, errmax=0.6235, trainrmse=0.0503, validrmse=0.0666, testrmse=0.0658eV 
!PES3: totrmse=0.0483, errmax=0.0657, trainrmse=0.0465, validrmse=0.0556, testrmse=0.0678eV 
!Average potential energy surface: totormse=0.0426eV
!***************************************************************************
      use nnparam
      implicit none
      integer i,j,k,ndriv
      real*8 v,vpes(npes),c(0:ninput),dvdx(1:natom*3,npes)
      real*8 xcart(1:3,1:natom),m(0:4),xmorse(1:nbond),r(1:nbond)
      real*8 dvdr(1:nbond),x(1:natom*3),drdx(1:nbond,1:natom*3)
      real*8 txinput(1:ninput),dvdg(1:ninput,npes),p(0:ninput)
      real*8 dpdr(1:nbond,0:ninput),dvdxa(1:natom*3)
      ! ::::::::::::::::::::
      do i=1,3
       x((3*(i-1)+1):(3*(i-1)+3))=xcart(:,i)
      enddo

      k=0
      do i=1,2
       do j=i+1,3
        k=k+1
       r(k)=dot_product(xcart(:,i)-xcart(:,j),
     $                  xcart(:,i)-xcart(:,j))
        r(k)=dsqrt(r(k))
        xmorse(k)=dexp(-r(k)/alpha)
       enddo
      enddo
      if(k.ne.3)stop "error"

      call bemsav(xmorse,m,p)

      do j=1,12
      txinput(j)=p(j)
      enddo

      call pot3a(txinput,vpes,dvdg,ndriv)

      v=0d0 
      do i=1,npes
       v=v+vpes(i)
      enddo      

      v=v/npes                   
 
      if(ndriv.eq.1)then
       call evdvdr(r,c,m,p,dpdr)
       call evdrdx(xcart,r,x,drdx)
       call evdvdx(dvdg,dpdr,drdx,dvdxa,dvdx)
      endif

      if(ndriv.eq.2)then
        write (*,*) 'Only energy and gradient are available'
      endif
      
      return
      end subroutine evvdvdx
!****************************************************************!
!-->  read NN weights and biases from matlab output
      subroutine pes_init(path)
      use nnparam
      implicit none
      character(len=1024), intent(in) :: path
      integer i,j,ihid,iwe,inode1,inode2,ilay1,ilay2
      integer ibasis,npd,iterm,ib,nfile1,nfile2
      character f1*80
      character(len=1024) :: file_path1, file_path2

      nfile1=4
      nfile2=7
      file_path1 = trim(path)//"/NO2/weights_NO2_6Ap_PIPNN.txt"
      file_path2 = trim(path)//"/NO2/biases_NO2_6Ap_PIPNN.txt"

      open(nfile1,file=file_path1,status='old')
      open(nfile2,file=file_path2,status='old')

!      allocate(weight3a(100,100,2:4,npes),
!     &bias3a(100,2:4,npes),nodes3a(10,npes),
!     &pdel3a(50000,npes),pavg3a(50000,npes))

       read(nfile1,*)
       read(nfile2,*) 
       read(nfile1,*)ninput,nhid3,noutput
       nscale=ninput+noutput
       nlayer3=nhid3+2 !additional one for input layer and one for output 
       allocate(nodes3a(nlayer3),pdel3a(nscale,npes),
     &pavg3a(nscale,npes))
       nodes3a(1)=ninput
       nodes3a(nlayer3)=noutput
       read(nfile1,*)(nodes3a(ihid),ihid=2,nhid3+1)
       nodemax3=0
       do i=1,nlayer3
        nodemax3=max(nodemax3,nodes3a(i))
       enddo
       allocate(weight3a(nodemax3,nodemax3,2:nlayer3,npes),
     &bias3a(nodemax3,2:nlayer3,npes))
       read(nfile1,*)ifunc3,nwe3

      do j=1,npes
       if(j.gt.1) then
        read(nfile1,*)
        read(nfile1,*)
        read(nfile1,*)
        read(nfile1,*)
        read(nfile2,*)
       endif
       read(nfile1,*)(pdel3a(i,j),i=1,nscale)
       read(nfile1,*)(pavg3a(i,j),i=1,nscale)
       iwe=0
       do ilay1=2,nlayer3
       ilay2=ilay1-1
       do inode1=1,nodes3a(ilay1)
       do inode2=1,nodes3a(ilay2) !
       read(nfile1,*)weight3a(inode2,inode1,ilay1,j)
       iwe=iwe+1
       enddo
       read(nfile2,*)bias3a(inode1,ilay1,j)
       iwe=iwe+1
       enddo
       enddo
       if (iwe.ne.nwe3) then
         write(*,*)'provided number of parameters ',nwe3
         write(*,*)'actual number of parameters ',iwe
         write(*,*)'nwe not equal to iwe, check input files or code'
         stop
       endif
      enddo

      close(nfile1)
      close(nfile2)

      end subroutine pes_init

!*************************************************************************
      subroutine pot3a(x,vpot3,dvdg,ndriv) !-->3-HOCO,a
      use nnparam
      implicit none
      integer i,ndriv,inode1,inode2,ilay1,ilay2
      integer j,k,m,neu1,neu2
      real*8 x(ninput),y(nodemax3,nlayer3,npes),vpot3(npes)
      real*8,allocatable::nxw3(:,:),nxw2(:,:,:),nxw1(:,:,:)
      real*8,allocatable::ax(:,:),bx(:,:)
      real*8 vrange,girange(ninput)
      real*8 dvdg(1:ninput,npes),dvtmp(npes)
      real*8, external::tranfun

      dvdg=0d0
      dvtmp=0d0

      do m=1,npes
       do i=1,ninput
        y(i,1,m)=(x(i)-pavg3a(i,m))/pdel3a(i,m)
       enddo

!-->.....evaluate the hidden layer
       do ilay1=2,nlayer3-1
        ilay2=ilay1-1
        do inode1=1,nodes3a(ilay1)
         y(inode1,ilay1,m)=bias3a(inode1,ilay1,m)
         do inode2=1,nodes3a(ilay2)
          y(inode1,ilay1,m)=y(inode1,ilay1,m)+y(inode2,ilay2,m)
     &*weight3a(inode2,inode1,ilay1,m)
         enddo
         y(inode1,ilay1,m)=tranfun(y(inode1,ilay1,m),ifunc3)
        enddo
       enddo

!-->.....now evaluate the output
       ilay1=nlayer3
       ilay2=ilay1-1
       do inode1=1,nodes3a(ilay1)
        y(inode1,ilay1,m)=bias3a(inode1,ilay1,m)
        do inode2=1,nodes3a(ilay2)
         y(inode1,ilay1,m)=y(inode1,ilay1,m)+y(inode2,ilay2,m)
     &*weight3a(inode2,inode1,ilay1,m)
        enddo
       enddo

!-->.....the value of output layer is the fitted potntial 
       vpot3(m)=y(nodes3a(nlayer3),nlayer3,m)
     &*pdel3a(nscale,m)+pavg3a(nscale,m)
      enddo

      if(ndriv.eq.1)then
       neu1=nodes3a(2);neu2=nodes3a(3)
       allocate(nxw1(1:ninput,1:neu1,npes),nxw2(1:neu1,1:neu2,npes),
     $ax(1:neu1,npes),bx(1:neu2,npes),nxw3(1:neu2,npes))
       do m=1,npes
        do i=1,ninput
         do j=1,neu1
          nxw1(i,j,m)=weight3a(i,j,2,m)
         enddo
        enddo
        do j=1,neu1
         do k=1,neu2
          nxw2(j,k,m)=weight3a(j,k,3,m)
         enddo
        enddo
        do j=1,neu1
         ax(j,m)=y(j,2,m)
        enddo
        do k=1,neu2
         bx(k,m)=y(k,3,m)
         nxw3(k,m)=weight3a(k,1,4,m)
        enddo
!        stop

        do i=1,ninput
         do k=1,neu2
          dvtmp=0d0
          do j=1,neu1
           dvtmp(m)=dvtmp(m)+nxw2(j,k,m)*nxw1(i,j,m)*(1d0-ax(j,m)**2)
          enddo
          dvdg(i,m)=dvdg(i,m)+nxw3(k,m)*dvtmp(m)*(1d0-bx(k,m)**2)
         enddo
         dvdg(i,m)=dvdg(i,m)*pdel3a(nscale,m)/pdel3a(i,m)
        enddo

       enddo
      endif
            
      return
      end subroutine pot3a
!**************************************************************************
        function tranfun(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun,x
c    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
c    ifunc=2, transfer function is log sigmoid function, 'logsig'
c    ifunc=3, transfer function is pure linear function, 'purelin'. It is imposed to the output layer by default
        if (ifunc.eq.1) then
        tranfun=dtanh(x)
        else if (ifunc.eq.2) then
        tranfun=1d0/(1d0+exp(-x))
        else if (ifunc.eq.3) then
        tranfun=x
        endif
        return
        end function tranfun
!**************************************************************************
        FUNCTION LOGSIG(X)
        REAL*8 X,LOGSIG
        LOGSIG=1.d0/(1.d0+DEXP(-X))
        RETURN
        END FUNCTION LOGSIG
!!**************************************************************************
!      subroutine evvdvdx(xcart,v,vpes,dvdxa,dvdx,ndriv)
!      use nnparam
!      implicit none
!      integer i,j,k,ndriv
!      real*8 v,vpes(npes),c(0:ninput),dvdx(1:natom*3,npes)
!      real*8 xcart(1:3,1:natom),m(0:4),xmorse(1:nbond),r(1:nbond)
!      real*8 dvdr(1:nbond),x(1:natom*3),drdx(1:nbond,1:natom*3)
!      real*8 txinput(1:ninput),dvdg(1:ninput,npes),p(0:ninput)
!      real*8 dpdr(1:nbond,0:ninput),dvdxa(1:natom*3)
!      ! ::::::::::::::::::::
!      do i=1,3
!       x((3*(i-1)+1):(3*(i-1)+3))=xcart(:,i)
!      enddo
!
!      k=0
!      do i=1,2
!       do j=i+1,3
!        k=k+1
!       r(k)=dot_product(xcart(:,i)-xcart(:,j),
!     $                  xcart(:,i)-xcart(:,j))
!        r(k)=dsqrt(r(k))
!        xmorse(k)=dexp(-r(k)/alpha)
!       enddo
!      enddo
!      if(k.ne.3)stop "error"
!
!      call bemsav(xmorse,m,p)
!
!      do j=1,12
!      txinput(j)=p(j)
!      enddo
!
!      call pot3a(txinput,vpes,dvdg,ndriv)
!
!      v=0d0 
!      do i=1,npes
!       v=v+vpes(i)
!      enddo      
!
!      v=v/npes                   
! 
!      if(ndriv.eq.1)then
!       call evdvdr(r,c,m,p,dpdr)
!       call evdrdx(xcart,r,x,drdx)
!       call evdvdx(dvdg,dpdr,drdx,dvdxa,dvdx)
!      endif
!      
!      return
!      end subroutine evvdvdx
!****************************************************************!
      subroutine bemsav(xmorse,m,p)
      implicit none
      real*8 xmorse(1:3)
      real*8 p(0:12)
      real*8 m(0:4)

       call evmono(xmorse,m)
       call evpoly(m,p)

       return
      end subroutine bemsav

      subroutine evmono(x,m)
      implicit none
      real*8 x(1:3)
      real*8 m(0:4)

      m(0) = 1.0D0
      m(1) = x(3)
      m(2) = x(2)
      m(3) = x(1)
      m(4) = m(1)*m(2)

      return
      end subroutine evmono

      subroutine evpoly(m,p)
      implicit none
      real*8 p(0:12)
      real*8 m(0:4)

      p(0) = m(0)
      p(1) = m(1) + m(2)
      p(2) = m(3)
      p(3) = m(4)
      p(4) = p(2)*p(1)
      p(5) = p(1)*p(1) - p(3) - p(3)
      p(6) = p(2)*p(2)
      p(7) = p(2)*p(3)
      p(8) = p(3)*p(1)
      p(9) = p(2)*p(5)
      p(10) = p(2)*p(4)
      p(11) = p(1)*p(5) - p(8)
      p(12) = p(2)*p(6)

      return
      end subroutine evpoly
!****************************************************************
      subroutine evdvdx(dvdg,dpdr,drdx,dvdxa,dvdx)
      use nnparam
      implicit none
      integer i,j,k
      real*8 dvdx(1:natom*3,npes),xcart(1:3,1:natom)
      real*8 drdx(1:nbond,1:natom*3),x(1:ninput)
      real*8 dvdr(1:nbond),dpdr(1:nbond,0:ninput),dvdp(0:ninput)
      real*8 dgdx(1:natom*3,1:ninput),dgdr(1:nbond,1:ninput)
      real*8 dvdxa(1:natom*3),dvdg(1:ninput,npes),rtmp
 
      dgdx=0d0
      dgdr(1:3,1:12)=dpdr(1:3,1:12)
      do i=1,9
       do j=1,12
        do k=1,3
         dgdx(i,j)=dgdx(i,j) + dgdr(k,j)*drdx(k,i)
        enddo
       enddo
      enddo

      dvdx=0d0;dvdxa=0
      
      do k=1,npes
       do i=1,9
        do j=1,12
         dvdx(i,k)=dvdx(i,k) + dvdg(j,k)*dgdx(i,j)
        enddo
       enddo
      enddo

      do i=1,9
       do k=1,npes
        dvdxa(i)=dvdxa(i)+dvdx(i,k)
       enddo
      enddo

      dvdxa=dvdxa/dble(npes)

      return
      end subroutine evdvdx
!*********************************************************************** 
      subroutine evdvdr(r,c,m,p,dpdr)
      use nnparam
      implicit none
      integer i,j
      real*8 r(1:nbond),dmsdr(1:nbond,1:nbond),dmdr(1:nbond,0:4)
      real*8 dvdr(1:nbond),dpdr(1:nbond,0:ninput),c(0:ninput)
      real*8 m(0:4),p(0:ninput)

      dvdr(:)=0d0
      call evdmsdr(r,dmsdr)
      call evdmdr(m,dmsdr,dmdr)
      call evdpdr(p,dmdr,dpdr)

      return
      end subroutine evdvdr
!***********************************************************************
      subroutine evdrdx(xcart,r,x,drdx)
      use nnparam
      implicit none
      integer i,j,k
      real*8 r(1:nbond),xcart(1:3,1:natom),drdx(1:nbond,1:natom*3)
      real*8 x(1:natom*3)
      drdx(:,:)=0d0
!dr1dx
      drdx(1,1)=(x(1)-x(4))/r(1)
      drdx(1,2)=(x(2)-x(5))/r(1)
      drdx(1,3)=(x(3)-x(6))/r(1)
      drdx(1,4)=-drdx(1,1)
      drdx(1,5)=-drdx(1,2)
      drdx(1,6)=-drdx(1,3)
!dr2dx
      drdx(2,1)=(x(1)-x(7))/r(2)
      drdx(2,2)=(x(2)-x(8))/r(2)
      drdx(2,3)=(x(3)-x(9))/r(2)
      drdx(2,7)=-drdx(2,1)
      drdx(2,8)=-drdx(2,2)
      drdx(2,9)=-drdx(2,3)
!dr3dx
      drdx(3,4)=(x(4)-x(7))/r(3)
      drdx(3,5)=(x(5)-x(8))/r(3)
      drdx(3,6)=(x(6)-x(9))/r(3)
      drdx(3,7)=-drdx(3,4)
      drdx(3,8)=-drdx(3,5)
      drdx(3,9)=-drdx(3,6)
 
      return
      end subroutine evdrdx
!*************************************************************************
      subroutine evdmsdr(r,dmsdr)
      use nnparam
      implicit none
      integer i,j
      real*8 dmsdr(1:nbond,1:nbond),r(1:nbond)
! ::::::::::::::::::::
      dmsdr(:,:)=0d0
 
!Morse term ms = exp(-r/alpha)
! dmsdr(i,j)=0  i!=j
! dmsdr(i,i)= -(1/alpha)*Exp(-r(i)/alpha)
 
      do i=1,nbond
       dmsdr(i,i)=-(1d0/alpha)*dexp(-r(i)/alpha)
      enddo
 
      return
      end subroutine evdmsdr
!**************************************************************************
      subroutine evdmdr(m,dmsdr,dmdr)
      use nnparam
      implicit none
      integer i,j
      real*8 dmsdr(1:nbond,1:nbond),dmdr(1:nbond,0:4),m(0:4)
 
      do i=1,nbond
       dmdr(i,0)=0.D0
       dmdr(i,1)=dmsdr(i,3)
       dmdr(i,2)=dmsdr(i,2)
       dmdr(i,3)=dmsdr(i,1)
       dmdr(i,4)=dmdr(i,1)*m(2)+m(1)*dmdr(i,2)
      enddo

      return
      end subroutine evdmdr
!*************************************************************************
      subroutine evdpdr(p,dmdr,dpdr)
      use nnparam
      implicit none
      integer i,j
      real*8 dmdr(1:nbond,0:4),dpdr(1:nbond,0:ninput)
      real*8 p(0:ninput)
 
      do i=1,nbond
       dpdr(i,0)=dmdr(i,0)
       dpdr(i,1)=dmdr(i,1)+dmdr(i,2)
       dpdr(i,2)=dmdr(i,3)
       dpdr(i,3)=dmdr(i,4)
       dpdr(i,4)=dpdr(i,2)*P(1)+P(2)*dpdr(i,1)
       dpdr(i,5)=dpdr(i,1)*p(1)+p(1)*dpdr(i,1)-dpdr(i,3)-dpdr(i,3)
       dpdr(i,6)=dpdr(i,2)*P(2)+P(2)*dpdr(i,2)
       dpdr(i,7)=dpdr(i,2)*P(3)+P(2)*dpdr(i,3)
       dpdr(i,8)=dpdr(i,3)*P(1)+P(3)*dpdr(i,1)
       dpdr(i,9)=dpdr(i,2)*P(5)+P(2)*dpdr(i,5)
       dpdr(i,10)=dpdr(i,2)*P(4)+P(2)*dpdr(i,4)
       dpdr(i,11)=dpdr(i,1)*P(5)+P(1)*dpdr(i,5)-dpdr(i,8)
       dpdr(i,12)=dpdr(i,2)*P(6)+P(2)*dpdr(i,6)
      enddo

      return
      end subroutine evdpdr
