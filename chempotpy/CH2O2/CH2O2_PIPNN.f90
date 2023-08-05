
!******************************************************************************
!******************************************************************************
!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
      module nnparam1
      implicit none
      real*8,parameter::alpha=1.8d0,PI=3.141592653589793238d0,
     $radian=PI/180.0d0,bohr=0.5291772d0
      real*8,parameter::PARA2=27.2114d0
      integer,parameter::ndata=6136
      real*8,parameter::vpesmin=-560.0228149965d0
      integer,parameter::nbasis=102,natom=5,nbond=10
      integer ninput1,noutput1,nscale1
      integer nhid3,nlayer31,ifunc31,nwe3,nodemax31
      integer,allocatable::nodes3a1(:)
      real*8,allocatable::weight3a1(:,:,:),bias3a1(:,:),
     $pdel3a1(:),pavg3a1(:)
      end module nnparam1

      subroutine pes(x,igrad,path,p,g,d)

      use nnparam1
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=5
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)


      double precision :: v, tx(3,natoms)
      integer :: iatom, idir, j, istate
      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
        do idir=1,3
          tx(idir, iatom)=x(iatom, idir)
        enddo
      enddo

      if (igrad==0) then
        call pes_init(path)
        call evvdvdx(tx,v)
        deallocate(nodes3a1)
        deallocate(weight3a1)
        deallocate(bias3a1)
        deallocate(pdel3a1)
        deallocate(pavg3a1)
      else
        write (*,*) 'Only energy is available'
      endif

      do istate=1,nstates
        p(istate)=v
      enddo

      endsubroutine


      subroutine pes_init(path)
      use nnparam1
      implicit none
      character(len=1024), intent(in) :: path
      integer i,ihid,iwe,inode1,inode2,ilay1,ilay2
      integer ibasis,npd,iterm,ib,nfile1,nfile2
      character f1*80
      character(len=1024) :: file_path1, file_path2

      nfile1=4
      nfile2=7
 
      file_path1 = trim(path)//"/CH2O2/weights.txt"
      file_path2 = trim(path)//"/CH2O2/biases.txt"

      open(nfile1,file=file_path1,status='old')
      open(nfile2,file=file_path2,status='old')
      read(nfile1,*)ninput1,nhid3,noutput1
      nscale1=ninput1+noutput1
      nlayer31=nhid3+2 !additional one for input layer and one for output 
      allocate(nodes3a1(nlayer31),pdel3a1(nscale1),pavg3a1(nscale1))
      nodes3a1(1)=ninput1
      nodes3a1(nlayer31)=noutput1
      read(nfile1,*)(nodes3a1(ihid),ihid=2,nhid3+1)
      nodemax31=0
      do i=1,nlayer31
       nodemax31=max(nodemax31,nodes3a1(i))
      enddo
      allocate(weight3a1(nodemax31,nodemax31,2:nlayer31),
     $bias3a1(nodemax31,2:nlayer31))
      read(nfile1,*)ifunc31,nwe3
      read(nfile1,*)(pdel3a1(i),i=1,nscale1)
      read(nfile1,*)(pavg3a1(i),i=1,nscale1)
      iwe=0
      do ilay1=2,nlayer31
      ilay2=ilay1-1
      do inode1=1,nodes3a1(ilay1)
      do inode2=1,nodes3a1(ilay2) !
      read(nfile1,*)weight3a1(inode2,inode1,ilay1)
      iwe=iwe+1
      enddo
      read(nfile2,*)bias3a1(inode1,ilay1)
      iwe=iwe+1
      enddo
      enddo
      if (iwe.ne.nwe3) then
        write(*,*)'provided number of parameters ',nwe3
        write(*,*)'actual number of parameters ',iwe
        write(*,*)'nwe not equal to iwe, check input files or code'
        stop
      endif

      close(nfile1)
      close(nfile2)

      end subroutine pes_init

!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->
        subroutine pot3a1(x,vpot3) ! fhcl,a
        use nnparam1
        implicit none
        integer i,ndriv,inode1,inode2,ilay1,ilay2
        integer j,k,neu1,neu2
        real*8 x(ninput1),y(nodemax31,nlayer31),vpot3
        real*8,allocatable::nxw3(:),nxw2(:,:),nxw1(:,:),ax(:),bx(:)
        real*8 dvdg(ninput1),dvtmp
        real*8, external :: tranfun1

        dvdg=0d0
        dvtmp=0d0          
        do i=1,ninput1
          y(i,1)=(x(i)-pavg3a1(i))/pdel3a1(i)
        enddo

!-->....evaluate the hidden layer
        do ilay1=2,nlayer31-1
         ilay2=ilay1-1
         do inode1=1,nodes3a1(ilay1)
          y(inode1,ilay1)=bias3a1(inode1,ilay1)
          do inode2=1,nodes3a1(ilay2)
           y(inode1,ilay1)=y(inode1,ilay1)+
     $y(inode2,ilay2)*weight3a1(inode2,inode1,ilay1)
          enddo
          y(inode1,ilay1)=tranfun1(y(inode1,ilay1),ifunc31)
         enddo
        enddo

!-->....now evaluate the output
        ilay1=nlayer31
        ilay2=ilay1-1
        do inode1=1,nodes3a1(ilay1)
        y(inode1,ilay1)=bias3a1(inode1,ilay1)
        do inode2=1,nodes3a1(ilay2)
        y(inode1,ilay1)=y(inode1,ilay1)+
     $y(inode2,ilay2)*weight3a1(inode2,inode1,ilay1)
        enddo
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot3=y(nodes3a1(nlayer31),nlayer31)*pdel3a1(nscale1)+
     &  pavg3a1(nscale1)

!        write(*,'(f16.8)')vpot3
!        stop
        return
        end subroutine pot3a1
!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->!-->

        function tranfun1(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun1,x
!    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
!    ifunc=2, transfer function is log sigmoid function, 'logsig'
!    ifunc=3, transfer function is pure linear function, 'purelin'. It is imposed to the output layer by default
        if (ifunc.eq.1) then
        tranfun1=dtanh(x)
        else if (ifunc.eq.2) then
        tranfun1=1d0/(1d0+dexp(-x))
        else if (ifunc.eq.3) then
        tranfun1=x
        endif
        return
        end function tranfun1

        FUNCTION LOGSIG1(X)
        REAL*8 X,LOGSIG1
        LOGSIG1=1.d0/(1.d0+DEXP(-X))
        RETURN
        END FUNCTION LOGSIG1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evvdvdx(xcart,v)
      implicit none
      real*8,parameter::alpha=1.8d0,PARA=0.5291772d0
      integer,parameter::natom=5
      integer i,j,k,ndriv
      real*8 v,p(0:101),r(1:10)
      real*8 xcart(1:3,1:5),m(0:38),xmorse(1:10)
      real*8 x(1:15)
      real*8 txinput(1:101)

      do i=1,natom
       x((3*(i-1)+1):(3*(i-1)+3))=xcart(:,i)
      enddo

      k=0
      do i=1,natom-1
       do j=i+1,natom
        k=k+1
       r(k)=dsqrt(dot_product(xcart(:,i)-xcart(:,j),
     $                  xcart(:,i)-xcart(:,j)))
!        r(k)=dsqrt(r(k))
        xmorse(k)=dexp(-r(k)/(PARA*alpha))
       enddo
      enddo
      if(k.ne.(natom*(natom-1)/2))stop "error"
      call bemsav1(xmorse,p)

      do j=1,101
       txinput(j)=p(j)
      enddo

!     write(66,'(1000E16.8)') txinput(1:101)
!      write(66,*) p(0:193)
!      stop

      call pot3a1(txinput,v)

      return
      end subroutine evvdvdx
!********************************************************
      subroutine bemsav1(x,p)
      implicit none
      real*8 x(1:10)
      real*8 p(0:101)
      real*8 m(0:38)
      
      call evmono1(x,m)
      call evpoly1(m,p)
      
      return
      end subroutine bemsav1
    
      subroutine evmono1(x,m)
      implicit none
      real*8 x(1:10)
      real*8 m(0:38)
      
      m(0) = 1.0D0
      m(1) = x(10)
      m(2) = x(9)
      m(3) = x(8)
      m(4) = x(7)
      m(5) = x(4)
      m(6) = x(6)
      m(7) = x(5)
      m(8) = x(3)
      m(9) = x(2)
      m(10) = x(1)
      m(11) = m(1)*m(2)
      m(12) = m(4)*m(5)
      m(13) = m(2)*m(6)
      m(14) = m(1)*m(7)
      m(15) = m(2)*m(8)
      m(16) = m(1)*m(9)
      m(17) = m(5)*m(6)
      m(18) = m(5)*m(7)
      m(19) = m(4)*m(8)
      m(20) = m(4)*m(9)
      m(21) = m(7)*m(8)
      m(22) = m(6)*m(9)
      m(23) = m(6)*m(8)
      m(24) = m(7)*m(9)
      m(25) = m(6)*m(7)
      m(26) = m(8)*m(9)
      m(27) = m(2)*m(17)
      m(28) = m(1)*m(18)
      m(29) = m(2)*m(19)
      m(30) = m(1)*m(20)
      m(31) = m(2)*m(23)
      m(32) = m(1)*m(24)
      m(33) = m(5)*m(25)
      m(34) = m(4)*m(26)
      m(35) = m(6)*m(21)
      m(36) = m(6)*m(24)
      m(37) = m(6)*m(26)
      m(38) = m(7)*m(26)
  
      return
      end subroutine evmono1
  
      subroutine evpoly1(m,p)
      implicit none
      real*8 x(1:10)
      real*8 p(0:101)
      real*8 m(0:38)
      
      p(0) = m(0)
      p(1) = m(1) + m(2)
      p(2) = m(3)
      p(3) = m(4) + m(5)
      p(4) = m(6) + m(7) + m(8) + m(9)
      p(5) = m(10)
      p(6) = m(11)
      p(7) = p(2)*p(1)
      p(8) = p(1)*p(3)
      p(9) = p(2)*p(3)
      p(10) = m(12)
      p(11) = m(13) + m(14) + m(15) + m(16)
      p(12) = p(1)*p(4) - p(11)
      p(13) = p(2)*p(4)
      p(14) = m(17) + m(18) + m(19) + m(20)
      p(15) = m(21) + m(22)
      p(16) = m(23) + m(24)
      p(17) = p(3)*p(4) - p(14)
      p(18) = m(25) + m(26)
      p(19) = p(5)*p(1)
      p(20) = p(2)*p(5)
      p(21) = p(5)*p(3)
      p(22) = p(5)*p(4)
      p(23) = p(1)*p(1) - p(6) - p(6)
      p(24) = p(2)*p(2)
      p(25) = p(3)*p(3) - p(10) - p(10)
      p(26) = p(4)*p(4) - p(18) - p(16) - p(15) - p(18) - p(16) - p(15)
      p(27) = p(5)*p(5)
      p(28) = p(2)*p(6)
      p(29) = p(6)*p(3)
      p(30) = p(2)*p(8)
      p(31) = p(10)*p(1)
      p(32) = p(2)*p(10)
      p(33) = p(6)*p(4)
      p(34) = p(2)*p(11)
      p(35) = p(2)*p(12)
      p(36) = m(27) + m(28) + m(29) + m(30)
      p(37) = p(1)*p(14) - p(36)
      p(38) = p(2)*p(14)
      p(39) = p(1)*p(15)
      p(40) = p(2)*p(15)
      p(41) = m(31) + m(32)
      p(42) = p(1)*p(16) - p(41)
      p(43) = p(2)*p(16)
      p(44) = p(3)*p(11) - p(36)
      p(45) = p(1)*p(17) - p(44)
      p(46) = p(2)*p(17)
      p(47) = p(10)*p(4)
      p(48) = p(3)*p(15)
      p(49) = p(3)*p(16)
      p(50) = p(1)*p(18)
      p(51) = p(2)*p(18)
      p(52) = m(33) + m(34)
      p(53) = m(35) + m(36) + m(37) + m(38)
      p(54) = p(3)*p(18) - p(52)
      p(55) = p(5)*p(6)
      p(56) = p(2)*p(19)
      p(57) = p(5)*p(8)
      p(58) = p(2)*p(21)
      p(59) = p(5)*p(10)
      p(60) = p(5)*p(11)
      p(61) = p(5)*p(12)
      p(62) = p(2)*p(22)
      p(63) = p(5)*p(14)
      p(64) = p(5)*p(15)
      p(65) = p(5)*p(16)
      p(66) = p(5)*p(17)
      p(67) = p(5)*p(18)
      p(68) = p(6)*p(1)
      p(69) = p(2)*p(23)
      p(70) = p(2)*p(7)
      p(71) = p(3)*p(23)
      p(72) = p(2)*p(9)
      p(73) = p(1)*p(25)
      p(74) = p(2)*p(25)
      p(75) = p(10)*p(3)
      p(76) = p(1)*p(11) - p(33)
      p(77) = p(1)*p(12) - p(33)
      p(78) = p(2)*p(13)
      p(79) = p(3)*p(14) - p(47)
      p(80) = p(3)*p(17) - p(47)
      p(81) = p(4)*p(11) - p(50) - p(41) - p(39) - p(41)
      p(82) = p(1)*p(26) - p(81)
      p(83) = p(2)*p(26)
      p(84) = p(4)*p(14) - p(52) - p(49) - p(48) - p(52)
      p(85) = p(4)*p(15) - p(53)
      p(86) = p(4)*p(16) - p(53)
      p(87) = p(3)*p(26) - p(84)
      p(88) = p(4)*p(18) - p(53)
      p(89) = p(5)*p(23)
      p(90) = p(2)*p(20)
      p(91) = p(5)*p(25)
      p(92) = p(5)*p(26)
      p(93) = p(5)*p(19)
      p(94) = p(2)*p(27)
      p(95) = p(5)*p(21)
      p(96) = p(5)*p(22)
      p(97) = p(1)*p(23) - p(68)
      p(98) = p(2)*p(24)
      p(99) = p(3)*p(25) - p(75)
      p(100) = p(4)*p(26) - p(88) - p(86) - p(85)
      p(101) = p(5)*p(27)
  
      return
      end subroutine evpoly1
  
