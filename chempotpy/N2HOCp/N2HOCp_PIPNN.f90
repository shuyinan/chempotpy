!     Subroutine to generate values of the 9D PES
!     N2 + HOC+ reactions
!     (C)2019 Guo Group, UNM
!     2019.5.12 Qian Yao updated
       module nnparam
       implicit none
       real*8,parameter::alpha=1.5d0,PI=3.1415926d0,radian=PI/180.0d0
       integer,parameter::nbasis=167,ndim=10,natom=5
       integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
       integer nscale
       integer, allocatable::nodes(:)

       real*8, allocatable::weighta(:,:,:),biasa(:,:)
       real*8, allocatable::pdela(:),pavga(:)
       real*8, allocatable::weightb(:,:,:),biasb(:,:)
       real*8, allocatable::pdelb(:),pavgb(:)
       real*8, allocatable::weightc(:,:,:),biasc(:,:)
       real*8, allocatable::pdelc(:),pavgc(:)

       end module nnparam

      subroutine pes(x,igrad,path,p,g,d)

      use nnparam
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=5
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: v, va, vb, vc
      double precision :: coord(3,natoms), r(10)
      integer :: iatom, jatom, idir, j, istate
      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
        do idir=1,3
          coord(idir,iatom)=x(iatom, idir)
        enddo
      enddo

      j=0
      do iatom=1,natoms
        do jatom=iatom+1, natoms
          j=j+1
          r(j)=sqrt((coord(1,iatom)-coord(1,jatom))**2 +&
                + (coord(2,iatom)-coord(2,jatom))**2 +&
                + (coord(3,iatom)-coord(3,jatom))**2)
        enddo
      enddo


      if (igrad==0) then
        call pes_init(path)
        call N2CHOpipNN(r,v,va,vb,vc)
        deallocate(nodes)
        deallocate(weighta)
        deallocate(biasa)
        deallocate(pdela)
        deallocate(pavga)
        deallocate(weightb)
        deallocate(biasb)
        deallocate(pdelb)
        deallocate(pavgb)
        deallocate(weightc)
        deallocate(biasc)
        deallocate(pdelc)
        deallocate(pavgc)
      else
        write (*,*) 'Only energy is available'
      endif
      do istate=1,nstates
        p(istate)=v
      enddo

      g=0.d0
      d=0.d0

      endsubroutine

        double precision function distance(a,b)
        double precision :: a(3), b(3)

        distance=0.d0
        do i=1,3
        distance=distance+(a(i)-b(i))**2
        enddo
        distance=dsqrt(distance)

        return
        end

       subroutine N2CHOpipNN(rb,vpes,vpesa,vpesb,vpesc) 
       use nnparam
       implicit none
       integer i,j,k,lpes
       real*8 rb(ndim),xbond(ndim),basis(0:nbasis),tmp1,txinput(nbasis)
       real*8 vpes,vpesa,vpesb,vpesc,vpesH,vnn

       basis=0.d0

       xbond(:)=dexp(-rb(:)/alpha)

      call bemsav21113(xbond,basis)
       do j=1,nbasis
        txinput(j)=basis(j)
       enddo
      
      call getpota(txinput,vpesa)
      call getpotb(txinput,vpesb)
      call getpotc(txinput,vpesc)
        vpes=(vpesa+vpesb+vpesc)/3d0 

       return

       end subroutine N2CHOpipNN

!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weights.txt'
!-->  biases saved in 'biases.txt'
!-->  one has to call this subroutine one and only one before calling
!the getpot() subroutine
      subroutine pes_init(path)
      use nnparam
      implicit none
      character(len=1024), intent(in) :: path
      integer i,ihid,iwe,inode1,inode2,ilay1,ilay2
      integer ibasis,npd,iterm,ib,nfile
      character f1*80,line*80
      character(len=1024) :: file_path1, file_path2
      character(len=1024) :: file_path3, file_path4
      character(len=1024) :: file_path5, file_path6

      file_path1 = trim(path)//"/N2HOCp/weights.txt-1"
      file_path2 = trim(path)//"/N2HOCp/biases.txt-1"

      open(4,file=file_path2,status='old')
        nfile=7
        open(nfile,file=file_path1)
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
        allocate(weighta(nodemax,nodemax,2:nlayer),biasa(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for
!hidden layers
!-->....At this time, only an equivalent transfer function can be used
!for all hidden layers
!-->....and the pure linear function is always applid to the output
!layer.
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

      close(4);close(nfile)
 
      file_path3 = trim(path)//"/N2HOCp/weights.txt-2"
      file_path4 = trim(path)//"/N2HOCp/biases.txt-2"
      open(4,file=file_path4,status='old')
        nfile=8
        open(nfile,file=file_path3)
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(pdelb(nscale),pavgb(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weightb(nodemax,nodemax,2:nlayer),biasb(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for
!hidden layers
!-->....At this time, only an equivalent transfer function can be used
!for all hidden layers
!-->....and the pure linear function is always applid to the output
!layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdelb(i),i=1,nscale)
        read(nfile,*)(pavgb(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weightb(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasb(inode1,ilay1)
        iwe=iwe+1
        enddo
        enddo
        
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif

      close(4);close(nfile)

      file_path5 = trim(path)//"/N2HOCp/weights.txt-3"
      file_path6 = trim(path)//"/N2HOCp/biases.txt-3"
      open(4,file=file_path6,status='old')
        nfile=9
        open(nfile,file=file_path5)
        read(nfile,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 !additional one for input layer and one for output 
        allocate(pdelc(nscale),pavgc(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(nfile,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
         nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weightc(nodemax,nodemax,2:nlayer),biasc(nodemax,2:nlayer))
        read(nfile,*)ifunc,nwe
!-->....ifunc hence controls the type of transfer function used for
!hidden layers
!-->....At this time, only an equivalent transfer function can be used
!for all hidden layers
!-->....and the pure linear function is always applid to the output
!layer.
!-->....see function tranfun() for details
        read(nfile,*)(pdelc(i),i=1,nscale)
        read(nfile,*)(pavgc(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        do inode2=1,nodes(ilay2) !
        read(nfile,*)weightc(inode2,inode1,ilay1)
        iwe=iwe+1
        enddo
        read(4,*)biasc(inode1,ilay1)
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
        end 

        function tranfun(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun,x
        if (ifunc.eq.1) then
        tranfun=dtanh(x)
        else if (ifunc.eq.2) then
        tranfun=1d0/(1d0+dexp(-x))
        else if (ifunc.eq.3) then
        tranfun=x
        endif
        return
        end
!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->subroutine to accelerate~~~~~~
        subroutine getpot(x,vpot)
        use nnparam
        implicit none
        integer i,ndriv,inode1,inode2,ilay1,ilay2,j,k,neu1,neu2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer

        do i=1,ninput
          y(i,1)=(x(i)-pavga(i))/pdela(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*weighta(inode2,inode1,ilay1)
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
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*weighta(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdela(nscale)+pavga(nscale)

        return
        end subroutine getpot

      function emsav21113(x,c) result(v)
      implicit none
      real*8,dimension(1:10)::x
      real*8,dimension(0:167)::c
      real*8::v
      ! ::::::::::::::::::::
      real*8,dimension(0:167)::p
      call bemsav21113(x,p)
      v = dot_product(p,c)
      return
      end function emsav21113
 
      subroutine bemsav21113(x,p)
      implicit none
      real*8,dimension(1:10),intent(in)::x
      real*8,dimension(0:167),intent(out)::p
      ! ::::::::::::::::::::
      real*8,dimension(0:21)::m
      call evmono21113(x,m)
      call evpoly21113(m,p)
      return
      end subroutine bemsav21113
 
      subroutine evmono21113(x,m)
      implicit none
      real*8,dimension(1:10),intent(in)::x
      real*8,dimension(0:21),intent(out)::m
 
      m(0)=1.D0
      m(1)=x(10)
      m(2)=x(9)
      m(3)=x(8)
      m(4)=x(7)
      m(5)=x(4)
      m(6)=x(6)
      m(7)=x(3)
      m(8)=x(5)
      m(9)=x(2)
      m(10)=x(1)
      m(11)=m(4)*m(5)
      m(12)=m(5)*m(6)
      m(13)=m(4)*m(7)
      m(14)=m(6)*m(7)
      m(15)=m(5)*m(8)
      m(16)=m(4)*m(9)
      m(17)=m(7)*m(8)
      m(18)=m(6)*m(9)
      m(19)=m(8)*m(9)
      m(20)=m(5)*m(17)
      m(21)=m(4)*m(18)
 
      return
      end subroutine evmono21113
 
      subroutine evpoly21113(m,p)
      implicit none
      real*8,dimension(0:21),intent(in)::m
      real*8,dimension(0:167),intent(out)::p
 
      p(0)=m(0)
      p(1)=m(1)
      p(2)=m(2)
      p(3)=m(3)
      p(4)=m(4)+m(5)
      p(5)=m(6)+m(7)
      p(6)=m(8)+m(9)
      p(7)=m(10)
      p(8)=p(1)*p(2)
      p(9)=p(1)*p(3)
      p(10)=p(2)*p(3)
      p(11)=p(1)*p(4)
      p(12)=p(2)*p(4)
      p(13)=p(3)*p(4)
      p(14)=m(11)
      p(15)=p(1)*p(5)
      p(16)=p(2)*p(5)
      p(17)=p(3)*p(5)
      p(18)=m(12)+m(13)
      p(19)=m(14)
      p(20)=p(4)*p(5)-p(18)
      p(21)=p(1)*p(6)
      p(22)=p(2)*p(6)
      p(23)=p(3)*p(6)
      p(24)=m(15)+m(16)
      p(25)=m(17)+m(18)
      p(26)=m(19)
      p(27)=p(4)*p(6)-p(24)
      p(28)=p(5)*p(6)-p(25)
      p(29)=p(1)*p(7)
      p(30)=p(2)*p(7)
      p(31)=p(3)*p(7)
      p(32)=p(7)*p(4)
      p(33)=p(7)*p(5)
      p(34)=p(7)*p(6)
      p(35)=p(1)*p(1)
      p(36)=p(2)*p(2)
      p(37)=p(3)*p(3)
      p(38)=p(4)*p(4)-p(14)-p(14)
      p(39)=p(5)*p(5)-p(19)-p(19)
      p(40)=p(6)*p(6)-p(26)-p(26)
      p(41)=p(7)*p(7)
      p(42)=p(1)*p(10)
      p(43)=p(1)*p(12)
      p(44)=p(1)*p(13)
      p(45)=p(2)*p(13)
      p(46)=p(1)*p(14)
      p(47)=p(2)*p(14)
      p(48)=p(3)*p(14)
      p(49)=p(1)*p(16)
      p(50)=p(1)*p(17)
      p(51)=p(2)*p(17)
      p(52)=p(1)*p(18)
      p(53)=p(2)*p(18)
      p(54)=p(3)*p(18)
      p(55)=p(1)*p(19)
      p(56)=p(2)*p(19)
      p(57)=p(3)*p(19)
      p(58)=p(1)*p(20)
      p(59)=p(2)*p(20)
      p(60)=p(3)*p(20)
      p(61)=p(14)*p(5)
      p(62)=p(19)*p(4)
      p(63)=p(1)*p(22)
      p(64)=p(1)*p(23)
      p(65)=p(2)*p(23)
      p(66)=p(1)*p(24)
      p(67)=p(2)*p(24)
      p(68)=p(3)*p(24)
      p(69)=p(1)*p(25)
      p(70)=p(2)*p(25)
      p(71)=p(3)*p(25)
      p(72)=m(20)+m(21)
      p(73)=p(1)*p(26)
      p(74)=p(2)*p(26)
      p(75)=p(3)*p(26)
      p(76)=p(1)*p(27)
      p(77)=p(2)*p(27)
      p(78)=p(3)*p(27)
      p(79)=p(14)*p(6)
      p(80)=p(4)*p(25)-p(72)
      p(81)=p(26)*p(4)
      p(82)=p(1)*p(28)
      p(83)=p(2)*p(28)
      p(84)=p(3)*p(28)
      p(85)=p(5)*p(24)-p(72)
      p(86)=p(19)*p(6)
      p(87)=p(26)*p(5)
      p(88)=p(4)*p(28)-p(85)
      p(89)=p(1)*p(30)
      p(90)=p(1)*p(31)
      p(91)=p(2)*p(31)
      p(92)=p(1)*p(32)
      p(93)=p(2)*p(32)
      p(94)=p(3)*p(32)
      p(95)=p(7)*p(14)
      p(96)=p(1)*p(33)
      p(97)=p(2)*p(33)
      p(98)=p(3)*p(33)
      p(99)=p(7)*p(18)
      p(100)=p(7)*p(19)
      p(101)=p(7)*p(20)
      p(102)=p(1)*p(34)
      p(103)=p(2)*p(34)
      p(104)=p(3)*p(34)
      p(105)=p(7)*p(24)
      p(106)=p(7)*p(25)
      p(107)=p(7)*p(26)
      p(108)=p(7)*p(27)
      p(109)=p(7)*p(28)
      p(110)=p(1)*p(8)
      p(111)=p(1)*p(36)
      p(112)=p(1)*p(9)
      p(113)=p(2)*p(10)
      p(114)=p(1)*p(37)
      p(115)=p(2)*p(37)
      p(116)=p(1)*p(11)
      p(117)=p(2)*p(12)
      p(118)=p(3)*p(13)
      p(119)=p(1)*p(38)
      p(120)=p(2)*p(38)
      p(121)=p(3)*p(38)
      p(122)=p(14)*p(4)
      p(123)=p(1)*p(15)
      p(124)=p(2)*p(16)
      p(125)=p(3)*p(17)
      p(126)=p(4)*p(18)-p(61)
      p(127)=p(4)*p(20)-p(61)
      p(128)=p(1)*p(39)
      p(129)=p(2)*p(39)
      p(130)=p(3)*p(39)
      p(131)=p(5)*p(18)-p(62)
      p(132)=p(19)*p(5)
      p(133)=p(4)*p(39)-p(131)
      p(134)=p(1)*p(21)
      p(135)=p(2)*p(22)
      p(136)=p(3)*p(23)
      p(137)=p(4)*p(24)-p(79)
      p(138)=p(5)*p(25)-p(86)
      p(139)=p(4)*p(27)-p(79)
      p(140)=p(5)*p(28)-p(86)
      p(141)=p(1)*p(40)
      p(142)=p(2)*p(40)
      p(143)=p(3)*p(40)
      p(144)=p(6)*p(24)-p(81)
      p(145)=p(6)*p(25)-p(87)
      p(146)=p(26)*p(6)
      p(147)=p(4)*p(40)-p(144)
      p(148)=p(5)*p(40)-p(145)
      p(149)=p(1)*p(29)
      p(150)=p(2)*p(30)
      p(151)=p(3)*p(31)
      p(152)=p(7)*p(38)
      p(153)=p(7)*p(39)
      p(154)=p(7)*p(40)
      p(155)=p(1)*p(41)
      p(156)=p(2)*p(41)
      p(157)=p(3)*p(41)
      p(158)=p(7)*p(32)
      p(159)=p(7)*p(33)
      p(160)=p(7)*p(34)
      p(161)=p(1)*p(35)
      p(162)=p(2)*p(36)
      p(163)=p(3)*p(37)
      p(164)=p(4)*p(38)-p(122)
      p(165)=p(5)*p(39)-p(132)
      p(166)=p(6)*p(40)-p(146)
      p(167)=p(7)*p(41)
 
      return
      end subroutine evpoly21113

        subroutine getpota(x,vpot)
        use nnparam
        implicit none
        integer i,ndriv,inode1,inode2,ilay1,ilay2,j,k,neu1,neu2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer

        do i=1,ninput
          y(i,1)=(x(i)-pavga(i))/pdela(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasa(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*weighta(inode2,inode1,ilay1)
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
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*weighta(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdela(nscale)+pavga(nscale)

        return
        end subroutine getpota

        subroutine getpotb(x,vpot)
        use nnparam
        implicit none
        integer i,ndriv,inode1,inode2,ilay1,ilay2,j,k,neu1,neu2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer

        do i=1,ninput
          y(i,1)=(x(i)-pavgb(i))/pdelb(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasb(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*weightb(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasb(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*weightb(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelb(nscale)+pavgb(nscale)

        return
        end subroutine getpotb

        subroutine getpotc(x,vpot)
        use nnparam
        implicit none
        integer i,ndriv,inode1,inode2,ilay1,ilay2,j,k,neu1,neu2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer

        do i=1,ninput
          y(i,1)=(x(i)-pavgc(i))/pdelc(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer-1
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasc(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*weightc(inode2,inode1,ilay1)
        enddo
        y(inode1,ilay1)=tranfun(y(inode1,ilay1),ifunc)
        enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer
        ilay2=ilay1-1
        do inode1=1,nodes(ilay1)
        y(inode1,ilay1)=biasc(inode1,ilay1)
        do inode2=1,nodes(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)*weightc(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelc(nscale)+pavgc(nscale)

        return
        end subroutine getpotc

