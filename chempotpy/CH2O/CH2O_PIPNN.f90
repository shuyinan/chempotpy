!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
       module nnparam
       implicit none
       real*8,parameter::alpha=1.0d0,PI=3.1415926d0,radian=PI/180.0d0
       integer,parameter::nbasis=18,ndim=6,natom=4
       real*8,parameter::vpesmin=-114.158373618668d0
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
        call ch2opipNN(tx,v)
        deallocate(nodes)
        deallocate(weighta)
        deallocate(biasa)
        deallocate(pdela)
        deallocate(pavga)
      else
        write (*,*) 'Only energy is available'
      endif

      do istate=1,nstates
        p(istate)=v
      enddo

      endsubroutine


       subroutine ch2opipNN(ct,vpes)
       use nnparam
       implicit none
       integer i,j,k,lpes
       real*8 rb(ndim),xbond(ndim),basis(nbasis),tmp1,txinput(nbasis-1)
       real*8 ct(3,natom),xvec(3,ndim),cx(3,natom),cy(3,natom)
       real*8 vpes,vpesa,vpesb,vpesc,vpesH,vnn,rfx(6)
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
!       if(rb(i).lt.0.65d0)then
!         vpes=5d0
!         return
!       endif
       enddo

       xbond(:)=dexp(-rb(:)/alpha)
       call bemsav(xbond,basis)

       do j=1,nbasis-1
        txinput(j)=basis(j+1)
       enddo

       call getpota(txinput,vpes)
!      vnn=min(vpes,5d0)
       vnn=vpes
       vpes=vnn

       return

       end subroutine ch2opipNN

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


      file_path1 = trim(path)//"/CH2O/weights.txt"
      file_path2 = trim(path)//"/CH2O/biases.txt"

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

        subroutine getpota(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        integer j,k,neu1,neu2,ndriv
        real*8 x(ninput),y(nodemax,nlayer),vpot
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
        tranfun=1d0/(1d0+dexp(-x))
        else if (ifunc.eq.3) then
        tranfun=x
        endif
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function emsav(x,c) result(v)
        implicit none
        real*8,dimension(1:6)::x
        real*8,dimension(0:17)::c
        real*8::v
        ! ::::::::::::::::::::
        real*8,dimension(0:17)::p
        
        call bemsav(x,p)
        v = dot_product(p,c)
        
        return
      end function emsav
      
      subroutine bemsav(x,p)
        implicit none
        real*8,dimension(1:6),intent(in)::x
        real*8,dimension(0:17),intent(out)::p
        ! ::::::::::::::::::::
        real*8,dimension(0:10)::m
        
        call evmono(x,m)
        call evpoly(m,p)
        
        return
      end subroutine bemsav
      
      subroutine evmono(x,m)
        implicit none
        real*8,dimension(1:6),intent(in)::x
        real*8,dimension(0:10),intent(out)::m
        !::::::::::::::::::::
        
        m(0) = 1.0D0
        m(1) = x(6)
        m(2) = x(5)
        m(3) = x(3)
        m(4) = x(4)
        m(5) = x(2)
        m(6) = x(1)
        m(7) = m(2)*m(3)
        m(8) = m(3)*m(4)
        m(9) = m(2)*m(5)
        m(10) = m(4)*m(5)

        return
      end subroutine evmono

      subroutine evpoly(m,p)
        implicit none
        real*8,dimension(0:10),intent(in)::m
        real*8,dimension(0:17),intent(out)::p
        !::::::::::::::::::::
        
        p(0) = m(0)
        p(1) = m(1)
        p(2) = m(2) + m(3)
        p(3) = m(4) + m(5)
        p(4) = m(6)
        p(5) = p(1)*p(2)
        p(6) = m(7)
        p(7) = p(1)*p(3)
        p(8) = m(8) + m(9)
        p(9) = m(10)
        p(10) = p(2)*p(3) - p(8)
        p(11) = p(1)*p(4)
        p(12) = p(4)*p(2)
        p(13) = p(4)*p(3)
        p(14) = p(1)*p(1)
        p(15) = p(2)*p(2) - p(6) - p(6)
        p(16) = p(3)*p(3) - p(9) - p(9)
        p(17) = p(4)*p(4)

        return
      end subroutine evpoly

