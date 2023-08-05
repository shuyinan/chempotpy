!-->  program to get potential energy for a given geometry after NN
!fitting
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
        double precision,parameter:: vpescut=5.0d0
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
        call hcloh_pes_interface(tx, v)
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

!**************************************************************************
        subroutine hcloh_pes_interface(ct,vpes)
        use nnparam,only: vpescut
        implicit none
        integer :: lxc
        double precision :: xct(3,4),ct(3,4),xvec(3,6),xct2(3,4)
        double precision :: vpes,vpesa,vpesb,vpesc
        double precision :: r13,r23

        vpes=0.d0
        vpesa=0.d0
        vpesb=0.d0
        vpesc=0.d0
        xct=ct
        
        xvec(:,2)=xct(:,3)-xct(:,1) !  Cl-H1
        xvec(:,4)=xct(:,3)-xct(:,2) !  H2-Cl
        r13=dsqrt(dot_product(xvec(:,2),xvec(:,2)))
        r23=dsqrt(dot_product(xvec(:,4),xvec(:,4)))
        
! check the linked atom order: make sure that
! 1st H linked or near Cl 
        if (r13 .gt. r23) then
           xct2(:,1)=xct(:,2)
           xct2(:,2)=xct(:,1)
           xct2(:,3)=xct(:,3)
           xct2(:,4)=xct(:,4)
        else
           xct2=xct
        endif        

        call before_vpes(xct2,lxc)

        if (lxc .eq. -1) then
           vpes=vpescut
           vpesa=vpescut
           vpesb=vpescut
           vpesc=vpescut
           return
        endif

        if (lxc .eq. 0) then
           call hclohNN(xct2,vpes,vpesa,vpesb,vpesc) 
           vpes=min(vpes,vpescut)
           return
        endif

        return
        end 
!**************************************************************************


!**************************************************************************
        subroutine hclohNN(ct,vpes,vpesa,vpesb,vpesc) 
! ct in HHClO order and angstrom
        use nnparam
        implicit none
        integer i
        real*8 rb(6),xbond(6),basis(0:17)
        real*8 txinput(17)
        real*8 ct(3,4),xvec(3,6),xct(3,4)
        real*8 vpesa,vpesb,vpesc,vpes
        double precision,parameter :: alpha=2.d0*0.5291772d0
 
        basis=0.d0
        xct=ct

        xvec(:,1)=xct(:,2)-xct(:,1) !  H2->H1
        xvec(:,2)=xct(:,3)-xct(:,1) !  Cl->H1
        xvec(:,3)=xct(:,4)-xct(:,1) !  O ->H1
        xvec(:,4)=xct(:,3)-xct(:,2) !  Cl->H2
        xvec(:,5)=xct(:,4)-xct(:,2) !  O ->H2
        xvec(:,6)=xct(:,4)-xct(:,3) !  O ->Cl

        rb(1)=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
        rb(2)=dsqrt(dot_product(xvec(:,2),xvec(:,2)))
        rb(3)=dsqrt(dot_product(xvec(:,3),xvec(:,3)))
        rb(4)=dsqrt(dot_product(xvec(:,4),xvec(:,4)))
        rb(5)=dsqrt(dot_product(xvec(:,5),xvec(:,5)))
        rb(6)=dsqrt(dot_product(xvec(:,6),xvec(:,6)))
 
        xbond(:)=dexp(-rb(:)/alpha)
 
        call bemsav(xbond,basis)
 
        do i=1,17
           txinput(i)=basis(i)
        enddo
 
        call getpota(txinput,vpesa)
        call getpotb(txinput,vpesb)
        call getpotc(txinput,vpesc)
        vpes=(vpesa+vpesb+vpesc)/3.0d0
 
        if(vpes.lt.-1.5d0)vpes=5.0d0
        vpes=min(vpes,5.0d0)
 
        return
 
        end subroutine hclohNN
!**************************************************************************


!**************************************************************************
!-->  read NN weights and biases from matlab output
!-->  weights saved in 'weightsa.txt','weightsb.txt','weightsc.txt'
!-->  biases saved in 'biasesa.txt','biasesb.txt','biasesc.txt'
!-->  one has to call this subroutine once and only once in the main
!code
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

        file_path1 = trim(path)//"/HClOH/weightsa.txt"
        file_path2 = trim(path)//"/HClOH/biasesa.txt"
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
       
        file_path3 = trim(path)//"/HClOH/weightsb.txt"
        file_path4 = trim(path)//"/HClOH/biasesb.txt"
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

        file_path5 = trim(path)//"/HClOH/weightsc.txt"
        file_path6 = trim(path)//"/HClOH/biasesc.txt"
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

        subroutine bemsav(x,p)
        implicit none
        double precision :: x(6),m(0:10),p(0:17)

        call evmono(x,m)
        call evpoly(m,p)

        return
        end subroutine bemsav
!**************************************************************************


!**************************************************************************
        subroutine evmono(x,m)
        implicit none
        double precision :: x(6),m(0:10)

        m(0)=1.d0
        m(1)=x(6)
        m(2)=x(5)
        m(3)=x(3)
        m(4)=x(4)
        m(5)=x(2)
        m(6)=x(1)
        m(7)=m(2)*m(3)
        m(8)=m(3)*m(4)
        m(9)=m(2)*m(5)
        m(10)=m(4)*m(5)

        return
        end subroutine evmono
!**************************************************************************


!**************************************************************************
        subroutine evpoly(m,p)
        implicit none
        double precision :: m(0:10),p(0:17)

        p(0)=m(0)
        p(1)=m(1)
        p(2)=m(2)+m(3)
        p(3)=m(4)+m(5)
        p(4)=m(6)
        p(5)=p(1)*p(2)
        p(6)=m(7)
        p(7)=p(1)*p(3)
        p(8)=m(8)+m(9)
        p(9)=m(10)
        p(10)=p(2)*p(3)-p(8)
        p(11)=p(1)*p(4)
        p(12)=p(4)*p(2)
        p(13)=p(4)*p(3)
        p(14)=p(1)*p(1)
        p(15)=p(2)*p(2)-p(6)-p(6)
        p(16)=p(3)*p(3)-p(9)-p(9)
        p(17)=p(4)*p(4)

        return
        end subroutine evpoly
!**************************************************************************


!**************************************************************************
        subroutine before_vpes(xct,lxc)
        implicit none
        integer :: lxc,i
        double precision :: xct(3,4),xct2(3,4),xvec(3,6),bt(6)
        double precision :: r12,r13,r14,r23,r24,r34

! lxc=-1, unphysical point
! lxc=0,  normal
! lxc=1,  HO+HCl, using switch function
! lxc=2,  HO+HCl, using dipole interaction correction

        lxc=0

        xct2=xct

        xvec(:,1)=xct2(:,2)-xct2(:,1) ! H2-H1
        xvec(:,2)=xct2(:,3)-xct2(:,1) ! Cl-H1
        xvec(:,3)=xct2(:,4)-xct2(:,1) ! O-H1
        xvec(:,4)=xct2(:,3)-xct2(:,2) ! H2-Cl
        xvec(:,5)=xct2(:,4)-xct2(:,2) ! H2-O
        xvec(:,6)=xct2(:,4)-xct2(:,3) ! Cl-O
        
        r12=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
        r13=dsqrt(dot_product(xvec(:,2),xvec(:,2)))
        r14=dsqrt(dot_product(xvec(:,3),xvec(:,3)))
        r23=dsqrt(dot_product(xvec(:,4),xvec(:,4)))
        r24=dsqrt(dot_product(xvec(:,5),xvec(:,5)))
        r34=dsqrt(dot_product(xvec(:,6),xvec(:,6)))

        bt(:)=(/r12,r13,r14,r23,r24,r34/)

! remove some unphysical points      
        do i=1,6
           if(bt(i).le.0.60d0)then
              lxc=-1
              return
           endif
        enddo

!--> ClO distance, do not consider  H2 + ClO channel      
!--> ClO equilibrium distance,1.596gstrom, x 130% = 2.07 angstrom
        if (r34 .le. 1.30d0 .or. r12 .le. 0.60d0) then
           lxc=-1
           return
        endif
        
        if (minval(bt) .ge. 2.00d0) then 
           lxc=-1
           return
        endif
        
        if (r13 .le. 0.75d0 .or. r23 .le. 0.75d0) then
           lxc=-1
           return
        endif

! H+O+HCl or HO + H +Cl or H+H+Cl+O      
        if (r14 .lt. r24 .and. r24 .gt. 2.0d0) then
           if (r14 .ge. 2.0d0 .or. r23 .ge. 2.0d0) then
              lxc=-1
              return
           endif
        endif
        
        if (r24 .lt. r14 .and. r14 .gt. 2.0d0) then
           if (r24 .ge. 2.0d0 .or. r13 .ge. 2.0d0) then
              lxc=-1
              return
           endif
        endif

! OH1 and OH2 are both large
        if (r14 .ge. 2.0d0 .and. r24 .ge. 2.0d0) then
           lxc=-1
           return
        endif

! ClH1 and ClH2 are both small
        if (r13 .le. 1.1d0 .and. r23 .le. 1.1d0) then
           lxc=-1
           return
        endif

        if (r23 .lt. r13 .and. r13 .gt. 1.9d0 .and. r24 .gt. 1.8d0) then
           if (r14 .ge. 1.8d0 .or. r23 .ge. 1.9d0) then   
              lxc=-1
              return
           endif
        endif
        
        if (r13 .lt. r23 .and. r23 .gt. 1.9d0 .and. r14 .gt. 1.8d0) then
           if (r13 .ge. 1.8d0 .or. r24 .ge. 1.9d0) then  
              lxc=-1
              return
           endif
        endif

        lxc=0
        return

        end subroutine before_vpes
!**************************************************************************
