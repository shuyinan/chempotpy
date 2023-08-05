!******************************************************************
!-->  program to get potential energy for a given geometry after NN
!     fitting
!-->  global variables are declared in this module
        module nnparam
        implicit none
        integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
        integer nscale
        integer, allocatable::nodes(:)
        real*8, allocatable::weighta(:,:,:),biasa(:,:)
        real*8, allocatable::pdela(:),pavga(:)
        real*8, allocatable::weightb(:,:,:),biasb(:,:)
        real*8, allocatable::pdelb(:),pavgb(:)
        real*8, allocatable::weightc(:,:,:),biasc(:,:)
        real*8, allocatable::pdelc(:),pavgc(:)
        double precision,parameter:: vpescut=6.d0 ! eV
        end module nnparam
!******************************************************************

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
        call ho3_pes_interface(tx, v)
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

      endsubroutine

!******************************************************************
        subroutine ho3_pes_interface(q,vpes)
        use nnparam,only: vpescut
        implicit none
        integer :: lxc
        double precision :: q(3,4) ! q in HOOO order and angstrom.
        double precision :: vpes,vpesa,vpesb,vpesc ! in eV

        vpes=0.d0
        vpesa=0.d0
        vpesb=0.d0
        vpesc=0.d0
        lxc=0

        call before_vpes(q,lxc)
        if (lxc.eq.-1) then
           vpesa=vpescut
           vpesb=vpescut
           vpesc=vpescut
           vpes=vpescut
           return
        endif

        if (lxc.eq.0) then
           call ho3NN(q,vpes,vpesa,vpesb,vpesc) 
           if (vpes.lt.-1.5 .or. vpes.gt.vpescut) then
              vpes=vpescut
              vpesa=vpescut
              vpesb=vpescut
              vpesc=vpescut
           endif
           return
        endif

        return
        end 
!******************************************************************

!******************************************************************
        subroutine ho3NN(ct,vpes,vpesa,vpesb,vpesc) 
        use nnparam
        implicit none
        integer i
        real*8 rb(6),xbond(6),basis(0:22)
        real*8 txinput(22)
        real*8 ct(3,4),qvec(3,6),xct(3,4)
        real*8 vpesa,vpesb,vpesc,vpes
        double precision,parameter :: alpha=2.5d0 ! bohr
        double precision, parameter :: BtoA=0.529177249d0 
        double precision,parameter :: PI=dacos(-1.d0)
 
        basis=0.d0
        xct=ct

        qvec(:,1)=xct(:,3)-xct(:,4) !  O3->O4
        qvec(:,2)=xct(:,2)-xct(:,4) !  O2->O4
        qvec(:,3)=xct(:,1)-xct(:,4) !  H1->O4
        qvec(:,4)=xct(:,2)-xct(:,3) !  O2->O3
        qvec(:,5)=xct(:,1)-xct(:,3) !  H1->O3
        qvec(:,6)=xct(:,1)-xct(:,2) !  H1->O2 

        rb(1)=dsqrt(dot_product(qvec(:,1),qvec(:,1)))
        rb(2)=dsqrt(dot_product(qvec(:,2),qvec(:,2)))
        rb(3)=dsqrt(dot_product(qvec(:,3),qvec(:,3)))
        rb(4)=dsqrt(dot_product(qvec(:,4),qvec(:,4)))
        rb(5)=dsqrt(dot_product(qvec(:,5),qvec(:,5)))
        rb(6)=dsqrt(dot_product(qvec(:,6),qvec(:,6)))
        
        xbond(:)=dexp(-rb(:)/(alpha*BtoA))
 
        call bemsav(xbond,basis)
 
        do i=1,22
           txinput(i)=basis(i)
        enddo
 
        call getpota(txinput,vpesa)
        call getpotb(txinput,vpesb)
        call getpotc(txinput,vpesc)
        vpes=(vpesa+vpesb+vpesc)/3.0d0
 
        return
        end subroutine ho3NN
!******************************************************************

!******************************************************************
!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weightsa.txt','weightsb.txt','weightsc.txt'
!-->  biases saved in 'biasesa.txt','biasesb.txt','biasesc.txt'
!-->  one has to call this subroutine once and only once in the 
!     main code
        subroutine pes_init(path)
        use nnparam
        implicit none
        character(len=1024), intent(in) :: path
        integer ihid,iwe,inode1,inode2,ilay1,ilay2
        integer i,wfilea,bfilea,wfileb,bfileb,wfilec,bfilec
        character(len=1024) :: file_path1, file_path2        
        character(len=1024) :: file_path3, file_path4
        character(len=1024) :: file_path5, file_path6
 
        wfilea=551
        bfilea=661 
        wfileb=552
        bfileb=662
        wfilec=553
        bfilec=663 

        file_path1 = trim(path)//"/HO3/weightsa.txt"
        file_path2 = trim(path)//"/HO3/biasesa.txt"
        open(wfilea,file=file_path1,status='old')
        open(bfilea,file=file_path2,status='old')

        rewind(wfilea)
        rewind(bfilea)

        read(wfilea,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2 
        allocate(nodes(nlayer),pdela(nscale),pavga(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(wfilea,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
           nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weighta(nodemax,nodemax,2:nlayer),
     %  biasa(nodemax,2:nlayer))
        read(wfilea,*)ifunc,nwe
        read(wfilea,*)(pdela(i),i=1,nscale)
        read(wfilea,*)(pavga(i),i=1,nscale)
        iwe=0
        do ilay1=2,nlayer
           ilay2=ilay1-1
           do inode1=1,nodes(ilay1)
               do inode2=1,nodes(ilay2) 
                   read(wfilea,*)weighta(inode2,inode1,ilay1)
                   iwe=iwe+1
               enddo
               read(bfilea,*)biasa(inode1,ilay1)
               iwe=iwe+1
           enddo
        enddo
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif
        close(wfilea)
        close(bfilea)
        
        file_path3 = trim(path)//"/HO3/weightsb.txt"
        file_path4 = trim(path)//"/HO3/biasesb.txt"
        open(wfileb,file=file_path3,status='old')
        open(bfileb,file=file_path4,status='old')

        rewind(wfileb)
        rewind(bfileb)

        read(wfileb,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2
        allocate(pdelb(nscale),pavgb(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(wfileb,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
           nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weightb(nodemax,nodemax,2:nlayer),
     %  biasb(nodemax,2:nlayer))
        read(wfileb,*)ifunc,nwe
        read(wfileb,*)(pdelb(i),i=1,nscale)
        read(wfileb,*)(pavgb(i),i=1,nscale)

        iwe=0
        do ilay1=2,nlayer
           ilay2=ilay1-1
           do inode1=1,nodes(ilay1)
               do inode2=1,nodes(ilay2)
                   read(wfileb,*)weightb(inode2,inode1,ilay1)
                   iwe=iwe+1
               enddo
               read(bfileb,*)biasb(inode1,ilay1)
               iwe=iwe+1
           enddo
        enddo
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif
        close(wfileb)
        close(bfileb)


        file_path5 = trim(path)//"/HO3/weightsc.txt"
        file_path6 = trim(path)//"/HO3/biasesc.txt"
        open(wfilec,file=file_path5,status='old')
        open(bfilec,file=file_path6,status='old')

        rewind(wfilec)
        rewind(bfilec)

        read(wfilec,*)ninput,nhid,noutput
        nscale=ninput+noutput
        nlayer=nhid+2
        allocate(pdelc(nscale),pavgc(nscale))
        nodes(1)=ninput
        nodes(nlayer)=noutput
        read(wfilec,*)(nodes(ihid),ihid=2,nhid+1)
        nodemax=0
        do i=1,nlayer
           nodemax=max(nodemax,nodes(i))
        enddo
        allocate(weightc(nodemax,nodemax,2:nlayer),
     %  biasc(nodemax,2:nlayer))
        read(wfilec,*)ifunc,nwe
        read(wfilec,*)(pdelc(i),i=1,nscale)
        read(wfilec,*)(pavgc(i),i=1,nscale)

        iwe=0
        do ilay1=2,nlayer
           ilay2=ilay1-1
           do inode1=1,nodes(ilay1)
               do inode2=1,nodes(ilay2)
                   read(wfilec,*)weightc(inode2,inode1,ilay1)
                   iwe=iwe+1
               enddo
               read(bfilec,*)biasc(inode1,ilay1)
               iwe=iwe+1
           enddo
        enddo
        if (iwe.ne.nwe) then
           write(*,*)'provided number of parameters ',nwe
           write(*,*)'actual number of parameters ',iwe
           write(*,*)'nwe not equal to iwe, check input files or code'
           stop
        endif
        close(wfilec)
        close(bfilec)

        return

        end subroutine pes_init
!**************************************************************************

!**************************************************************************
        subroutine getpota(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
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
                 y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &                           *weighta(inode2,inode1,ilay1)
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
     &                        *weighta(inode2,inode1,ilay1)
           enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdela(nscale)+pavga(nscale)
        return
        end subroutine getpota
!**************************************************************************

!**************************************************************************
        subroutine getpotb(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
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
                 y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &                           *weightb(inode2,inode1,ilay1)
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
              y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &                        *weightb(inode2,inode1,ilay1)
           enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelb(nscale)+pavgb(nscale)
        return
        end subroutine getpotb
!**************************************************************************

!**************************************************************************
        subroutine getpotc(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
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
                 y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &                           *weightc(inode2,inode1,ilay1)
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
              y(inode1,ilay1)=y(inode1,ilay1)+y(inode2,ilay2)
     &                        *weightc(inode2,inode1,ilay1)
           enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelc(nscale)+pavgc(nscale)
        return
        end subroutine getpotc
!**************************************************************************

!**************************************************************************
        function tranfun(x,ifunc)
        implicit none
        integer ifunc
        real*8 tranfun,x
c    ifunc=1, transfer function is hyperbolic tangent function, 'tansig'
c    ifunc=2, transfer function is log sigmoid function, 'logsig'
c    ifunc=3, transfer function is pure linear function, 'purelin'. It
c             is imposed to the output layer by default
        if (ifunc.eq.1) then
           tranfun=dtanh(x)
        else if (ifunc.eq.2) then
           tranfun=1d0/(1d0+exp(-x))
        else if (ifunc.eq.3) then
           tranfun=x
        endif
        return
        end
!**************************************************************************

!**************************************************************************
        subroutine before_vpes(xct,lxc)
        implicit none
        integer :: lxc
        double precision :: xct(3,4),xc(3,4),xvec(3,6),rb(6)
        double precision :: rO3O4,rO2O4,rH1O4,rO2O3,rH1O3,rH1O2 
! lxc=-1,  unphysical point
! lxc=0,   normal

        lxc=0
        xc=xct

        xvec(:,1)=xc(:,3)-xc(:,4) !  O3->O4
        xvec(:,2)=xc(:,2)-xc(:,4) !  O2->O4
        xvec(:,3)=xc(:,1)-xc(:,4) !  H1->O4
        xvec(:,4)=xc(:,2)-xc(:,3) !  O2->O3
        xvec(:,5)=xc(:,1)-xc(:,3) !  H1->O3
        xvec(:,6)=xc(:,1)-xc(:,2) !  H1->O2 
        
        rO3O4=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
        rO2O4=dsqrt(dot_product(xvec(:,2),xvec(:,2)))
        rH1O4=dsqrt(dot_product(xvec(:,3),xvec(:,3)))
        rO2O3=dsqrt(dot_product(xvec(:,4),xvec(:,4)))
        rH1O3=dsqrt(dot_product(xvec(:,5),xvec(:,5)))
        rH1O2=dsqrt(dot_product(xvec(:,6),xvec(:,6)))

        rb(:)=(/rO3O4,rO2O4,rH1O4,rO2O3,rH1O3,rH1O2/)

! remove some unphysical points
        if(min(rb(3),rb(5),rb(6)).le.0.65d0.or.
     &     min(rb(1),rb(2),rb(4)).le.0.65d0) then
           lxc=-1
           return
        endif

        end subroutine before_vpes
!**************************************************************************

!------>
!This program used to generate PIPs for molecule of A3B1 type
!up to 3rd order.
!------>

!********************************************************
        subroutine bemsav(x,p)
        implicit none
        double precision,intent(in) :: x(1:6)
        double precision,intent(out) :: p(0:22)
        double precision :: m(0:29)

        call evmono(x,m)
        call evpoly(m,p)

        return
        end subroutine bemsav
!********************************************************

!********************************************************
        subroutine evmono(x,m)
        implicit none
        double precision,intent(in) :: x(1:6)
        double precision,intent(out) :: m(0:29)

        m(0) = 1.0D0
        m(1) = x(6)
        m(2) = x(5)
        m(3) = x(3)
        m(4) = x(4)
        m(5) = x(2)
        m(6) = x(1)
        m(7) = m(1)*m(2)
        m(8) = m(1)*m(3)
        m(9) = m(2)*m(3)
        m(10) = m(3)*m(4)
        m(11) = m(2)*m(5)
        m(12) = m(1)*m(6)
        m(13) = m(4)*m(5)
        m(14) = m(4)*m(6)
        m(15) = m(5)*m(6)
        m(16) = m(1)*m(9)
        m(17) = m(1)*m(10)
        m(18) = m(2)*m(10)
        m(19) = m(1)*m(11)
        m(20) = m(3)*m(11)
        m(21) = m(2)*m(12)
        m(22) = m(3)*m(12)
        m(23) = m(2)*m(13)
        m(24) = m(3)*m(13)
        m(25) = m(1)*m(14)
        m(26) = m(3)*m(14)
        m(27) = m(1)*m(15)
        m(28) = m(2)*m(15)
        m(29) = m(4)*m(15)

        return
        end subroutine evmono
!********************************************************

!********************************************************
        subroutine evpoly(m,p)
        implicit none
        double precision,intent(in) :: m(0:29)
        double precision,intent(out) :: p(0:22)

        p(0) = m(0)                                          
        p(1) = m(1) + m(2) + m(3)
        p(2) = m(4) + m(5) + m(6)
        p(3) = m(7) + m(8) + m(9)
        p(4) = m(10) + m(11) + m(12)
        p(5) = p(1)*p(2) - p(4)
        p(6) = m(13) + m(14) + m(15)
        p(7) = p(1)*p(1) - p(3) - p(3)
        p(8) = p(2)*p(2) - p(6) - p(6)
        p(9) = m(16)
        p(10) = m(17) + m(18) + m(19) + m(20) + m(21) + m(22)
        p(11) = p(2)*p(3) - p(10)
        p(12) = m(23) + m(24) + m(25) + m(26) + m(27) + m(28)
        p(13) = m(29)
        p(14) = p(1)*p(6) - p(12)
        p(15) = p(1)*p(3) - p(9) - p(9) - p(9)
        p(16) = p(1)*p(4) - p(10)
        p(17) = p(2)*p(7) - p(16)
        p(18) = p(2)*p(4) - p(12)
        p(19) = p(1)*p(8) - p(18)
        p(20) = p(2)*p(6) - p(13) - p(13) - p(13)
        p(21) = p(1)*p(7) - p(15)
        p(22) = p(2)*p(8) - p(20)                            

        return
        end subroutine evpoly 
!********************************************************
