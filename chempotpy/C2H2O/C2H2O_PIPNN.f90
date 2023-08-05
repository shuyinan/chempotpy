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
        double precision,parameter:: vpescut=5.0d0 ! eV
        end module nnparam
!******************************************************************

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
        call c2h2o_pes_interface(tx,v)
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
!       q in OCCHH order and angstrom
!       vpes in eV
        subroutine c2h2o_pes_interface(q,vpes)
        use nnparam,only: vpescut
        implicit none
        integer :: lxc
        double precision :: q(3,5)
        double precision :: vpes,vpesa,vpesb,vpesc

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
           call c2h2oNN(q,vpes,vpesa,vpesb,vpesc) 
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
        subroutine c2h2oNN(ct,vpes,vpesa,vpesb,vpesc) 
        use nnparam
        implicit none
        integer i
        real*8 rb(10),xbond(10),basis(0:101)
        real*8 txinput(101)
! ct in OCCHH order and angstrom
        real*8 ct(3,5),qvec(3,10),xct(3,5)
        real*8 vpesa,vpesb,vpesc,vpes
        double precision,parameter :: alpha=2.d0 ! bohr
! bohr to angstrom
        double precision, parameter :: BtoA=0.529177249d0 
        double precision,parameter :: PI=dacos(-1.d0)
 
        basis=0.d0
        xct=ct

        qvec(:,1)=xct(:,3)-xct(:,2) !  C3->C2
        qvec(:,2)=xct(:,4)-xct(:,2) !  H4->C2
        qvec(:,3)=xct(:,5)-xct(:,2) !  H5->C2
        qvec(:,4)=xct(:,1)-xct(:,2) !  O1->C2
        qvec(:,5)=xct(:,4)-xct(:,3) !  H4->C3
        qvec(:,6)=xct(:,5)-xct(:,3) !  H5->C3 
        qvec(:,7)=xct(:,1)-xct(:,3) !  O1->C3 
        qvec(:,8)=xct(:,5)-xct(:,4) !  H5->H4 
        qvec(:,9)=xct(:,1)-xct(:,4) !  O1->H4 
        qvec(:,10)=xct(:,1)-xct(:,5) !  O1->H5 
        
        rb(1)=dsqrt(dot_product(qvec(:,1),qvec(:,1)))
        rb(2)=dsqrt(dot_product(qvec(:,2),qvec(:,2)))
        rb(3)=dsqrt(dot_product(qvec(:,3),qvec(:,3)))
        rb(4)=dsqrt(dot_product(qvec(:,4),qvec(:,4)))
        rb(5)=dsqrt(dot_product(qvec(:,5),qvec(:,5)))
        rb(6)=dsqrt(dot_product(qvec(:,6),qvec(:,6)))
        rb(7)=dsqrt(dot_product(qvec(:,7),qvec(:,7)))
        rb(8)=dsqrt(dot_product(qvec(:,8),qvec(:,8)))
        rb(9)=dsqrt(dot_product(qvec(:,9),qvec(:,9)))
        rb(10)=dsqrt(dot_product(qvec(:,10),qvec(:,10)))
        
        xbond(:)=dexp(-rb(:)/(alpha*BtoA))
 
        call bemsav(xbond,basis)
 
        do i=1,101
           txinput(i)=basis(i)
        enddo
 
        call getpota(txinput,vpesa)
        call getpotb(txinput,vpesb)
        call getpotc(txinput,vpesc)
        vpes=(vpesa+vpesb+vpesc)/3.0d0
 
        return
        end subroutine c2h2oNN
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

        file_path1 = trim(path)//"/C2H2O/weightsa.txt"
        file_path2 = trim(path)//"/C2H2O/biasesa.txt"
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
        
        file_path3 = trim(path)//"/C2H2O/weightsb.txt"
        file_path4 = trim(path)//"/C2H2O/biasesb.txt"
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

        file_path5 = trim(path)//"/C2H2O/weightsc.txt"
        file_path6 = trim(path)//"/C2H2O/biasesc.txt"
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
        integer :: i,lxc
        double precision :: xct(3,5),xc(3,5),xvec(3,10),rb(10)
        double precision :: rC2C3,rC2H4,rC2H5,rC2O1,rC3H4
        double precision :: rC3H5,rC3O1,rH4H5,rH4O1,rH5O1 
        double precision :: min_rCH,xx,yy,zz,tt
! lxc=-1,  unphysical point
! lxc=0,   normal

        lxc=0
        xc=xct

        xvec(:,1)=xc(:,3)-xc(:,2)  !  C3->C2
        xvec(:,2)=xc(:,4)-xc(:,2)  !  H4->C2
        xvec(:,3)=xc(:,5)-xc(:,2)  !  H5->C2
        xvec(:,4)=xc(:,1)-xc(:,2)  !  O1->C2
        xvec(:,5)=xc(:,4)-xc(:,3)  !  H4->C3
        xvec(:,6)=xc(:,5)-xc(:,3)  !  H5->C3 
        xvec(:,7)=xc(:,1)-xc(:,3)  !  O1->C3 
        xvec(:,8)=xc(:,5)-xc(:,4)  !  H5->H4 
        xvec(:,9)=xc(:,1)-xc(:,4)  !  O1->H4 
        xvec(:,10)=xc(:,1)-xc(:,5) !  O1->H5 
        
        rC2C3=dsqrt(dot_product(xvec(:,1),xvec(:,1)))    !  C3->C2
        rC2H4=dsqrt(dot_product(xvec(:,2),xvec(:,2)))    !  H4->C2
        rC2H5=dsqrt(dot_product(xvec(:,3),xvec(:,3)))    !  H5->C2
        rC2O1=dsqrt(dot_product(xvec(:,4),xvec(:,4)))    !  O1->C2
        rC3H4=dsqrt(dot_product(xvec(:,5),xvec(:,5)))    !  H4->C3
        rC3H5=dsqrt(dot_product(xvec(:,6),xvec(:,6)))    !  H5->C3
        rC3O1=dsqrt(dot_product(xvec(:,7),xvec(:,7)))    !  O1->C3
        rH4H5=dsqrt(dot_product(xvec(:,8),xvec(:,8)))    !  H5->H4
        rH4O1=dsqrt(dot_product(xvec(:,9),xvec(:,9)))    !  O1->H4
        rH5O1=dsqrt(dot_product(xvec(:,10),xvec(:,10)))  !  O1->H5

        rb(:)=(/rC2C3,rC2H4,rC2H5,rC2O1,rC3H4,
     $         rC3H5,rC3O1,rH4H5,rH4O1,rH5O1/)

        min_rCH=min(rC2H4,rC2H5,rC3H4,rC3H5)

! remove some unphysical points
        do i=1,10
           if (rb(i).le.0.7d0) then
              lxc=-1
!              print*,"lxc.eq.-1, rbond too small"
              return
           endif
        enddo

! CH single bond equilibrium distance,1.1199 angstrom,
!    1.1199*150% = 1.6798 angstrom
        xx=1.68d0
! CC single bond equilibrium distance,1.5360 angstrom,
!    1.5360*150% = 2.3040 angstrom
        yy=2.30d0
! CO single bond equilibrium distance,1.4270 angstrom,
!    1.4270*150% = 2.1405 angstrom
        zz=2.14d0 
! HO single bond equilibrium distance,0.9560 angstrom,
!    0.9560*150% = 1.4340 angstrom
        tt=1.43d0

        if (min_rCH.gt.xx) then
           lxc=-1
!           print*,"lxc.eq.-1, min_rCH too large"
           return
        else
           if(rC2C3.gt.yy) then
             if((rC2H4.le.xx .and. rC2H5.le.xx .and. rC3O1.le.zz).or.
     $          (rC3H4.le.xx .and. rC3H5.le.xx .and. rC2O1.le.zz))then
                lxc=0
!                print*,"lxc.eq.0, CH2+CO channel"
                return
             else
                lxc=-1
!                print*,"lxc.eq.-1, CC bond break but not CH2+CO channel"
                return
             endif
           else
             if((rC2H4.le.xx .and.rC3H5.le.xx) .or. 
     $          (rC2H5.le.xx .and.rC3H4.le.xx))then
                lxc=0
!                print*,"lxc.eq.0, HCCH+O channel"
                return
             elseif(((rC2H4.le.xx.or.rC2H5.le.xx).and.rC3O1.le.yy ).or.
     $          ((rC3H4.le.xx.or.rC3H5.le.xx) .and. rC2O1.le.yy ))then
                lxc=0
!                print*,"lxc.eq.0, HCCO+H channel"
                return
             elseif(((rC2H4.le.xx.or.rC3H4.le.xx).and.rH5O1.le.tt ).or.
     $          ((rC2H5.le.xx.or.rC3H5.le.xx) .and. rH4O1.le.tt ))then
                lxc=0
!                print*,"lxc.eq.0, HCC+OH channel"
                return
             elseif((rC2H4.le.xx.and.rC2H5.le.xx.and.rC3O1.le.yy ).or.
     $          (rC3H4.le.xx.and.rC3H5.le.xx.and.rC2O1.le.yy ))then
                lxc=0
!                print*,"lxc.eq.0, H2CCO channel"
                return
             else
                lxc=-1
!                print*,"lxc.eq.-1, CC bond not break, but not HCCH+O/
!     $                  HCCO+H/HCC+OH/H2CCO channel"
                return
             endif
           endif
        endif

        return
        end subroutine before_vpes


!**************************************************************************
!------->
!------->This is the subroutine for PIPs up to 3rd order
!------->of A2B2C type molecule.
!------->
!********************************************************
!        double precision :: x(1:10)
!        double precision :: p(0:101)
!        call bemsav(x,p)
!********************************************************

!********************************************************
        subroutine bemsav(x,p)
        implicit none
        double precision,intent(in) :: x(1:10)
        double precision :: m(0:38)
        double precision,intent(out) :: p(0:101)

        call evmono(x,m)
        call evpoly(m,p)

        return
        end subroutine bemsav
!********************************************************

!********************************************************
        subroutine evmono(x,m)
        implicit none
        double precision,intent(in) :: x(1:10)
        double precision,intent(out) :: m(0:38)

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
        end subroutine evmono
!********************************************************

!********************************************************
        subroutine evpoly(m,p)
        implicit none
        double precision,intent(in) :: m(0:38)
        double precision,intent(out) :: p(0:101)

        p(0)  = m(0) 
        p(1)  = m(1)+m(2)
        p(2)  = m(3)
        p(3)  = m(4)+m(5)
        p(4)  = m(6)+m(7)+m(8)+m(9)
        p(5)  = m(10)
        p(6)  = m(11)
        p(7)  = p(2)*p(1)
        p(8)  = p(1)*p(3)
        p(9)  = p(2)*p(3)
        p(10) = m(12)
        p(11) = m(13)+m(14)+m(15)+m(16)
        p(12) = p(1)*p(4)-p(11)
        p(13) = p(2)*p(4)
        p(14) = m(17)+m(18)+m(19)+m(20)
        p(15) = m(21)+m(22)
        p(16) = m(23)+m(24)
        p(17) = p(3)*p(4)-p(14)
        p(18) = m(25)+m(26)
        p(19) = p(5)*p(1)
        p(20) = p(2)*p(5)
        p(21) = p(5)*p(3)
        p(22) = p(5)*p(4)
        p(23) = p(1)*p(1)-p(6)-p(6)
        p(24) = p(2)*p(2)
        p(25) = p(3)*p(3)-p(10)-p(10)   
        p(26) = p(4)*p(4)-p(18)-p(16)-p(15)-p(18)-p(16)-p(15)
        p(27) = p(5)*p(5)
        p(28) = p(2)*p(6)
        p(29) = p(6)*p(3)
        p(30) = p(2)*p(8)
        p(31) = p(10)*p(1)
        p(32) = p(2)*p(10)
        p(33) = p(6)*p(4)
        p(34) = p(2)*p(11)
        p(35) = p(2)*p(12)
        p(36) = m(27)+m(28)+m(29)+m(30)
        p(37) = p(1)*p(14)-p(36)
        p(38) = p(2)*p(14)
        p(39) = p(1)*p(15)
        p(40) = p(2)*p(15)
        p(41) = m(31)+m(32)
        p(42) = p(1)*p(16)-p(41)
        p(43) = p(2)*p(16)
        p(44) = p(3)*p(11)-p(36)
        p(45) = p(1)*p(17)-p(44)
        p(46) = p(2)*p(17)
        p(47) = p(10)*p(4)
        p(48) = p(3)*p(15)
        p(49) = p(3)*p(16)
        p(50) = p(1)*p(18)
        p(51) = p(2)*p(18)
        p(52) = m(33)+m(34)
        p(53) = m(35)+m(36)+m(37)+m(38)
        p(54) = p(3)*p(18)-p(52)
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
        p(76) = p(1)*p(11)-p(33)
        p(77) = p(1)*p(12)-p(33)
        p(78) = p(2)*p(13)
        p(79) = p(3)*p(14)-p(47)
        p(80) = p(3)*p(17)-p(47)
        p(81) = p(4)*p(11)-p(50)-p(41)-p(39)-p(41)
        p(82) = p(1)*p(26)-p(81)
        p(83) = p(2)*p(26)
        p(84) = p(4)*p(14)-p(52)-p(49)-p(48)-p(52)
        p(85) = p(4)*p(15)-p(53)
        p(86) = p(4)*p(16)-p(53)
        p(87) = p(3)*p(26)-p(84)
        p(88) = p(4)*p(18)-p(53)
        p(89) = p(5)*p(23)
        p(90) = p(2)*p(20)
        p(91) = p(5)*p(25)
        p(92) = p(5)*p(26)
        p(93) = p(5)*p(19)
        p(94) = p(2)*p(27)
        p(95) = p(5)*p(21)
        p(96) = p(5)*p(22)
        p(97) = p(1)*p(23)-p(68)
        p(98) = p(2)*p(24)
        p(99) = p(3)*p(25)-p(75)
        p(100)= p(4)*p(26)-p(88)-p(86)-p(85)
        p(101)= p(5)*p(27)                                 

        return
        end subroutine evpoly 
!********************************************************
