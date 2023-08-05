       module nnparam
       implicit none
       real*8,parameter::alpha=1.0d0,PI=3.1415926d0,radian=PI/180.0d0
       real*8,parameter::PARA2=27.2116d0
       integer,parameter::nbasis=23,ndim=6,natom=4,ndata=33435
       real*8,parameter::vpesmin=-399.45863953d0
       integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
       integer nscale
       integer, allocatable::nodes(:)
       real*8, allocatable::weighta(:,:,:),biasa(:,:)
       real*8, allocatable::pdela(:),pavga(:)

       end module nnparam

      subroutine pes(x,igrad,path,p,g,d)

      use nnparam
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=4
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: v, va, vb, vc
      double precision :: tx(3,natoms), tg(3,natoms)
      integer :: iatom, idir, j, istate
      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
        do idir=1,3
          tx(idir,iatom)=x(iatom, idir)
        enddo
      enddo

        call pes_init(path)
        call evvdvdx(tx, v, tg, igrad)
        deallocate(nodes)
        deallocate(weighta)
        deallocate(biasa)
        deallocate(pdela)
        deallocate(pavga)

      do istate=1,nstates
        p(istate)=v/23.0609d0
      enddo

      do iatom=1,natoms
        do idir=1,3
          g(1,iatom,idir)=tg(idir,iatom)/23.0609d0
        enddo
      enddo

      !g=0.d0
      d=0.d0

      endsubroutine

       subroutine pipNN(ct,vpes,vpesa,vpesb,vpesc)
       use nnparam
       implicit none
       integer i,j,k,lpes
       real*8 rb(ndim),xbond(ndim),basis(nbasis),tmp1,txinput(nbasis-1)
       real*8 ct(3,natom),xvec(3,ndim),cx(3,natom),cy(3,natom)
       real*8 m(0:nbasis-1)
       real*8 vpes,vpesa,vpesb,vpesc,vpesH,vnn
       basis=0.d0

       cx=ct

       k=0
       do i=1,natom-1
        do j=i+1,natom
         k=k+1
!        write(*,*)k,i,j
         xvec(:,k)=cx(:,i)-cx(:,j)
        enddo
       enddo
       if(k.ne.ndim)stop 'error in bond dimension'
       do i=1,ndim
        rb(i)=dsqrt(dot_product(xvec(:,i),xvec(:,i)))
       enddo

       xbond(:)=dexp(-rb(:)/alpha)
       call bemsav(xbond,m,basis)

       do j=1,nbasis-1
        txinput(j)=basis(j+1)
       enddo

       call getpota(txinput,vpesa)
        vnn=vpesc
       vpes=vnn
!-->      write(*,*)vpes
       return

       end subroutine pipNN

!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weights.txt'
!-->  biases saved in 'biases.txt'
!-->  one has to call this subroutine one and only one before calling the getpot() subroutine
      subroutine pes_init(path)
      use nnparam
      implicit none
      character(len=1024), intent(in) :: path
      integer i,ihid,iwe,inode1,inode2,ilay1,ilay2
      integer ibasis,npd,iterm,ib,nfile
      character f1*80,line*80
      character(len=1024) :: file_path1, file_path2


      file_path1 = trim(path)//"/H3S/weights.txt"
      file_path2 = trim(path)//"/H3S/biase.txt"

      open(4,file=file_path2,status='old')
      read(4,'(a80)')line

        nfile=7
        open(nfile,file=file_path1)
        read(nfile,'(a80)')line
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(nodes(nlayer),pdela(nscale),pavga(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weighta(nodemax,nodemax,2:nlayer),
     %   biasa(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdela(i),i=1,nscale)
        read(nfile,*)(pavga(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weighta(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasa(inode1,ilay1)
        iwe=iwe+1
        enddo
        enddo
        
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif

        close(nfile)
        close(4)

        return
        end subroutine pes_init 

        subroutine getpota(x,vpot,dvdg,ndriv)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        integer j,k,neu1,neu2,ndriv
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8,allocatable::nxw3(:),nxw2(:,:),nxw1(:,:),ax(:),bx(:)
        real*8 dvdg(ninput),dvtmp
        real*8, external :: tranfun
!-->....set up the normalized input layer
c       write(*,*)ninput
        do i=1,ninput
          y(i,1)=(x(i)-pavga(i))/pdela(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weighta(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &*weighta(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdela(nscale)+pavga(nscale)
        if(ndriv.eq.1)then
         neu1=nodes(2);neu2=nodes(3)
         allocate(nxw1(1:ninput,1:neu1),nxw2(1:neu1,1:neu2),
     $ax(1:neu1),bx(1:neu2),nxw3(1:neu2))
         do i=1,ninput
          do j=1,neu1
           nxw1(i,j)=weighta(i,j,2)
          enddo
         enddo
         do j=1,neu1
          do k=1,neu2
           nxw2(j,k)=weighta(j,k,3)
          enddo
         enddo
         do j=1,neu1
          ax(j)=y(j,2)
         enddo
         do k=1,neu2
          bx(k)=y(k,3)
          nxw3(k)=weighta(k,1,4)
         enddo
!        stop

         dvdg=0d0
         dvtmp=0d0
         do i=1,ninput
          do k=1,neu2
           dvtmp=0d0
           do j=1,neu1
            dvtmp=dvtmp+nxw2(j,k)*nxw1(i,j)*(1d0-ax(j)**2)
           enddo
           dvdg(i)=dvdg(i)+nxw3(k)*dvtmp*(1d0-bx(k)**2)
           enddo
          dvdg(i)=dvdg(i)*pdela(nscale)/pdela(i)
         enddo
         deallocate(nxw1,nxw2,nxw3,ax,bx)
        endif


        return
        end subroutine getpota


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
        end



      subroutine evvdvdx(xcart,v,dvdx,ndriv)
      implicit none
      real*8,parameter::alpha=1d0
      integer i,j,k,ndriv
      real*8 v,dvdx(1:12),p(0:22),r(1:6)
      real*8 xcart(1:3,1:4),m(0:29),xmorse(1:6)
      real*8 dvdr(1:6),x(1:12),drdx(1:6,1:12)
      real*8 txinput(1:22),dvdg(1:22),dpdr(1:6,0:22)
      ! ::::::::::::::::::::
      do i=1,4
       x((3*(i-1)+1):(3*(i-1)+3))=xcart(:,i)
      enddo

      k=0
      do i=1,3
       do j=i+1,4
        k=k+1
       r(k)=dot_product(xcart(:,i)-xcart(:,j),
     $                  xcart(:,i)-xcart(:,j))
        r(k)=dsqrt(r(k))
        xmorse(k)=dexp(-r(k)/alpha)
       enddo
      enddo
      if(k.ne.6)stop "error"
      call bemsav(xmorse,m,p)
      do j=1,22
       txinput(j)=p(j)
      enddo
!      write(*,*)txinput(:)
!      stop
      call getpota(txinput,v,dvdg,ndriv)
    
      if(ndriv.eq.1)then
       call evdvdr(r,m,p,dpdr)
       call evdrdx(r,x,drdx)
       call evdvdx(dvdg,dpdr,drdx,dvdx)
      endif

      if(ndriv.eq.2)then  
        write(*,*) 'Only energy and gradient are available'
      endif

      return
      end subroutine evvdvdx
 
      subroutine bemsav(xmorse,m,p)
      implicit none
      real*8,dimension(1:6),intent(in)::xmorse
      real*8,dimension(0:22),intent(out)::p
      ! ::::::::::::::::::::
      real*8,dimension(0:29)::m
      call evmono(xmorse,m)
      call evpoly(m,p)
      return
      end subroutine bemsav
 
      subroutine evmono(x,m)
      implicit none
      real*8,dimension(1:6),intent(in)::x
      real*8,dimension(0:29),intent(out)::m
 
      m(0)=1.D0
      m(1)=x(6)
      m(2)=x(5)
      m(3)=x(3)
      m(4)=x(4)
      m(5)=x(2)
      m(6)=x(1)
      m(7)=m(1)*m(2)
      m(8)=m(1)*m(3)
      m(9)=m(2)*m(3)
      m(10)=m(3)*m(4)
      m(11)=m(2)*m(5)
      m(12)=m(1)*m(6)
      m(13)=m(4)*m(5)
      m(14)=m(4)*m(6)
      m(15)=m(5)*m(6)
      m(16)=m(1)*m(9)
      m(17)=m(1)*m(10)
      m(18)=m(2)*m(10)
      m(19)=m(1)*m(11)
      m(20)=m(3)*m(11)
      m(21)=m(2)*m(12)
      m(22)=m(3)*m(12)
      m(23)=m(2)*m(13)
      m(24)=m(3)*m(13)
      m(25)=m(1)*m(14)
      m(26)=m(3)*m(14)
      m(27)=m(1)*m(15)
      m(28)=m(2)*m(15)
      m(29)=m(4)*m(15)
 
      return
      end subroutine evmono
 
      subroutine evpoly(m,p)
      implicit none
      real*8,dimension(0:29),intent(in)::m
      real*8,dimension(0:22),intent(out)::p
 
       p(0)=m(0)
       p(1)=m(1)+m(2)+m(3)
       p(2)=m(4)+m(5)+m(6)
       p(3)=m(7)+m(8)+m(9)
       p(4)=m(10)+m(11)+m(12)
       p(5)=p(1)*p(2)-p(4)
       p(6)=m(13)+m(14)+m(15)
       p(7)=p(1)*p(1)-p(3)-p(3)
       p(8)=p(2)*p(2)-p(6)-p(6)
       p(9)=m(16)
       p(10)=m(17)+m(18)+m(19)+m(20)+m(21)+m(22)
       p(11)=p(2)*p(3)-p(10)
       p(12)=m(23)+m(24)+m(25)+m(26)+m(27)+m(28)
       p(13)=m(29)
       p(14)=p(1)*p(6)-p(12)
       p(15)=p(1)*p(3)-p(9)-p(9)-p(9)
       p(16)=p(1)*p(4)-p(10)
       p(17)=p(2)*p(7)-p(16)
       p(18)=p(2)*p(4)-p(12)
       p(19)=p(1)*p(8)-p(18)
       p(20)=p(2)*p(6)-p(13)-p(13)-p(13)
       p(21)=p(1)*p(7)-p(15)
       p(22)=p(2)*p(8)-p(20)
 
      return
      end subroutine evpoly



      subroutine evdvdx(dvdg,dpdr,drdx,dvdx)
      implicit none
      integer i,j,k
      real*8 dvdx(1:12),xcart(1:3,1:4),drdx(1:6,1:12),x(1:12)
      real*8 dvdr(1:6),dpdr(1:6,0:22),rtmp
      real*8 dgdx(1:12,1:22),dgdr(1:6,1:22),dvdg(1:22)
 
      dgdx=0d0
      dgdr(1:6,1:22)=dpdr(1:6,1:22)
      do i=1,12
       do j=1,22
        do k=1,6
         dgdx(i,j)=dgdx(i,j)+dgdr(k,j)*drdx(k,i)
        enddo
       enddo
      enddo
 
      dvdx=0d0
      do i=1,12
       do j=1,22
         dvdx(i)=dvdx(i)+dvdg(j)*dgdx(i,j)
       enddo
      enddo
 
      return
      end subroutine evdvdx
 
      subroutine evdvdr(r,m,p,dpdr)
      implicit none
      integer i,j
      real*8 r(1:6),dmsdr(1:6,1:6),dmdr(1:6,0:29)
      real*8 dvdr(1:6),dpdr(1:6,0:22)
      real*8 m(0:29),p(0:22)

      dvdr(:)=0d0
      call evdmsdr(r,dmsdr)
      call evdmdr(m,dmsdr,dmdr)
      call evdpdr(p,dmdr,dpdr)


      return
      end subroutine evdvdr
 
      subroutine evdrdx(r,x,drdx)
      implicit none
      integer i,j,k
      real*8 r(1:6),drdx(1:6,1:12),x(1:12)
      drdx(:,:)=0d0
      k=0
      do i=1,3
       do j=i+1,4
        k=k+1
!drdx
       drdx(k,3*(i-1)+1)=(x(3*(i-1)+1)-x(3*(j-1)+1))/r(k)
       drdx(k,3*(i-1)+2)=(x(3*(i-1)+2)-x(3*(j-1)+2))/r(k)
       drdx(k,3*(i-1)+3)=(x(3*(i-1)+3)-x(3*(j-1)+3))/r(k)
       drdx(k,3*(j-1)+1)=-drdx(k,3*(i-1)+1)
       drdx(k,3*(j-1)+2)=-drdx(k,3*(i-1)+2)
       drdx(k,3*(j-1)+3)=-drdx(k,3*(i-1)+3)
       enddo
      enddo
      if(k.ne.6)stop "error"
 
      return
      end subroutine evdrdx
 
      subroutine evdmsdr(r,dmsdr)
      implicit none
      real*8,parameter::alpha=1d0
      integer i,j
      real*8 dmsdr(1:6,1:6),r(1:6)
! ::::::::::::::::::::
      dmsdr(:,:)=0d0
 
!Morse term ms = exp(-r/alpha)
! dmsdr(i,j)=0  i!=j
! dmsdr(i,i)= -(1/alpha)*Exp(-r(i)/alpha)
 
      do i=1,6
       dmsdr(i,i)=-(1d0/alpha)*dexp(-r(i)/alpha)
      enddo
 
      return
      end subroutine evdmsdr
 
      subroutine evdmdr(m,dmsdr,dmdr)
      implicit none
      integer i,j
      real*8 dmsdr(1:6,1:6),dmdr(1:6,0:29),m(0:29)
 
      do i=1,6
       dmdr(i,0)=0.D0
       dmdr(i,1)=dmsdr(i,6)
       dmdr(i,2)=dmsdr(i,5)
       dmdr(i,3)=dmsdr(i,3)
       dmdr(i,4)=dmsdr(i,4)
       dmdr(i,5)=dmsdr(i,2)
       dmdr(i,6)=dmsdr(i,1)
       dmdr(i,7)=dmdr(i,1)*m(2)+m(1)*dmdr(i,2)
       dmdr(i,8)=dmdr(i,1)*m(3)+m(1)*dmdr(i,3)
       dmdr(i,9)=dmdr(i,2)*m(3)+m(2)*dmdr(i,3)
       dmdr(i,10)=dmdr(i,3)*m(4)+m(3)*dmdr(i,4)
       dmdr(i,11)=dmdr(i,2)*m(5)+m(2)*dmdr(i,5)
       dmdr(i,12)=dmdr(i,1)*m(6)+m(1)*dmdr(i,6)
       dmdr(i,13)=dmdr(i,4)*m(5)+m(4)*dmdr(i,5)
       dmdr(i,14)=dmdr(i,4)*m(6)+m(4)*dmdr(i,6)
       dmdr(i,15)=dmdr(i,5)*m(6)+m(5)*dmdr(i,6)
       dmdr(i,16)=dmdr(i,1)*m(9)+m(1)*dmdr(i,9)
       dmdr(i,17)=dmdr(i,1)*m(10)+m(1)*dmdr(i,10)
       dmdr(i,18)=dmdr(i,2)*m(10)+m(2)*dmdr(i,10)
       dmdr(i,19)=dmdr(i,1)*m(11)+m(1)*dmdr(i,11)
       dmdr(i,20)=dmdr(i,3)*m(11)+m(3)*dmdr(i,11)
       dmdr(i,21)=dmdr(i,2)*m(12)+m(2)*dmdr(i,12)
       dmdr(i,22)=dmdr(i,3)*m(12)+m(3)*dmdr(i,12)
       dmdr(i,23)=dmdr(i,2)*m(13)+m(2)*dmdr(i,13)
       dmdr(i,24)=dmdr(i,3)*m(13)+m(3)*dmdr(i,13)
       dmdr(i,25)=dmdr(i,1)*m(14)+m(1)*dmdr(i,14)
       dmdr(i,26)=dmdr(i,3)*m(14)+m(3)*dmdr(i,14)
       dmdr(i,27)=dmdr(i,1)*m(15)+m(1)*dmdr(i,15)
       dmdr(i,28)=dmdr(i,2)*m(15)+m(2)*dmdr(i,15)
       dmdr(i,29)=dmdr(i,4)*m(15)+m(4)*dmdr(i,15)
      enddo
 
      return
      end subroutine evdmdr
 
      subroutine evdpdr(p,dmdr,dpdr)
      implicit none
      integer i,j
      real*8 dmdr(1:6,0:29),dpdr(1:6,0:22)
      real*8 p(0:22)
 
      do i=1,6
       dpdr(i,0)=dmdr(i,0)
       dpdr(i,1)=dmdr(i,1)+dmdr(i,2)+dmdr(i,3)
       dpdr(i,2)=dmdr(i,4)+dmdr(i,5)+dmdr(i,6)
       dpdr(i,3)=dmdr(i,7)+dmdr(i,8)+dmdr(i,9)
       dpdr(i,4)=dmdr(i,10)+dmdr(i,11)+dmdr(i,12)
       dpdr(i,5)=dpdr(i,1)*p(2)+p(1)*dpdr(i,2)-dpdr(i,4)
       dpdr(i,6)=dmdr(i,13)+dmdr(i,14)+dmdr(i,15)
       dpdr(i,7)=dpdr(i,1)*p(1)+p(1)*dpdr(i,1)-dpdr(i,3)-dpdr(i,3)
       dpdr(i,8)=dpdr(i,2)*p(2)+p(2)*dpdr(i,2)-dpdr(i,6)-dpdr(i,6)
       dpdr(i,9)=dmdr(i,16)
       dpdr(i,10)=dmdr(i,17)+dmdr(i,18)+dmdr(i,19)+dmdr(i,20)+dmdr(i,21)
     &+dmdr(i,22)
       dpdr(i,11)=dpdr(i,2)*p(3)+p(2)*dpdr(i,3)-dpdr(i,10)
       dpdr(i,12)=dmdr(i,23)+dmdr(i,24)+dmdr(i,25)+dmdr(i,26)+dmdr(i,27)
     &+dmdr(i,28)
       dpdr(i,13)=dmdr(i,29)
       dpdr(i,14)=dpdr(i,1)*p(6)+p(1)*dpdr(i,6)-dpdr(i,12)
       dpdr(i,15)=dpdr(i,1)*p(3)+p(1)*dpdr(i,3)-dpdr(i,9)-dpdr(i,9)-dpdr
     &(i,9)
       dpdr(i,16)=dpdr(i,1)*p(4)+p(1)*dpdr(i,4)-dpdr(i,10)
       dpdr(i,17)=dpdr(i,2)*p(7)+p(2)*dpdr(i,7)-dpdr(i,16)
       dpdr(i,18)=dpdr(i,2)*p(4)+p(2)*dpdr(i,4)-dpdr(i,12)
       dpdr(i,19)=dpdr(i,1)*p(8)+p(1)*dpdr(i,8)-dpdr(i,18)
       dpdr(i,20)=dpdr(i,2)*p(6)+p(2)*dpdr(i,6)-dpdr(i,13)-dpdr(i,13)-dp
     &dr(i,13)
       dpdr(i,21)=dpdr(i,1)*p(7)+p(1)*dpdr(i,7)-dpdr(i,15)
       dpdr(i,22)=dpdr(i,2)*p(8)+p(2)*dpdr(i,8)-dpdr(i,20)
      enddo
 
      return
      end subroutine evdpdr
