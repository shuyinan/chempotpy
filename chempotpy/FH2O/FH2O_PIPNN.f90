!-->  program to get potential energy for a given geometry after NN
!fitting
!-->  global variables are declared in this module
       module nnparam
       implicit none
       real*8,parameter::alpha=1.0d0,vpesmin=-176.0260035639990d0,
     % PI=3.141592653589793238d0,radian=PI/180.0d0,bohr=0.5291772d0
       integer,parameter::nbasis=18,nrbasissoc=147,npbasissoc=288
       integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
       integer nterm(1:nbasis),nindex(1:nbasis,1:100,1:6)
       integer nscale
       integer, allocatable::nodes(:)
       real*8 r1soccoef(1:nrbasissoc),r2soccoef(1:nrbasissoc)
       real*8 p2soccoef(1:npbasissoc)

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
      integer, parameter :: natoms=4
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: v, va, vb, vc, socg, vgsoc
      double precision :: tx(3,natoms)
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

      if (igrad==0) then
        call pes_init(path)
        call fh2oNN(tx, v, va, vb, vc)
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
        write (*,*) 'Only energy and gradient are available'
      endif

      v=v  !/23.0609

      do istate=1,nstates
        p(istate)=v
      enddo

      g=0.d0
      d=0.d0

      endsubroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      subroutine fh2o_soc_nn(xt,vg,socg,vgsoc)
      use nnparam
      implicit none
      
      real*8 xt(3,4),xt2(3,4),xn(3,4),xvec(3,6),rint(6),rt(6)
      real*8 tr1basis(nrbasissoc),tr2basis(nrbasissoc) 
      real*8 r1function,r2function,rcom,the1,the2,phi
      real*8 p2function,tp2basis(npbasissoc)
      integer n,k,l,m,j,nrcom
      real*8 vg,vgsoc,socg,socs,vpesa,vpesb,vpesc
      real*8 socrg,socrs,socpg,socps
      real*8 r1switch,r2switch,p2switch,switch

      call fh2oNN(xt,vg,vpesa,vpesb,vpesc) !xt in HHFO and angstrom
      vg=vg*23.0605d0 
      vg=min(vg,100.0d0)
      if(vg.ge.100.0d0)then
        socg=0.0d0
        vgsoc=vg+socg/349.75d0
        return
      endif
      
      xt2=xt

      call orderxnt(xt2,xn) !xn 1st H closer to F

      xvec(:,1)=xn(:,2)-xn(:,1) !  H2->H1
      xvec(:,2)=xn(:,3)-xn(:,1) !  F->H1
      xvec(:,3)=xn(:,4)-xn(:,1) !  O->H1
      xvec(:,4)=xn(:,3)-xn(:,2) !  H2->F
      xvec(:,5)=xn(:,4)-xn(:,2) !  H2->O
      xvec(:,6)=xn(:,4)-xn(:,3) !  F-O 
     
      rt(1)=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
      rt(2)=dsqrt(dot_product(xvec(:,2),xvec(:,2)))
      rt(3)=dsqrt(dot_product(xvec(:,3),xvec(:,3))) !the breaking bond
      rt(4)=dsqrt(dot_product(xvec(:,4),xvec(:,4)))
      rt(5)=dsqrt(dot_product(xvec(:,5),xvec(:,5)))
      rt(6)=dsqrt(dot_product(xvec(:,6),xvec(:,6)))

cccccccc  reactant 
       xt2(:,1)=xt(:,1)
       xt2(:,2)=xt(:,4)
       xt2(:,3)=xt(:,2)
       xt2(:,4)=xt(:,3)

       call xyzTOjcb3(xt2,rint) !-->xt2 in HOHF order,angstrom, rt in angstrom and radians
       rcom=rint(3)
       the2=rint(5)
       phi=rint(6) 

       call soc_react_ground(rcom,the2,phi,socrg)

cccccccc  product       
       xt2(:,1)=xn(:,2)
       xt2(:,2)=xn(:,4)
       xt2(:,3)=xn(:,1)
       xt2(:,4)=xn(:,3)

       call cart_diJacobi(xt2,rint) !xt2 in HOHF order,angstrom, rt in angstrom and radians
       rcom=rint(3)
       the1=rint(4)
       the2=rint(5)
       phi=rint(6) 

       call soc_product_ground(rcom,the1,the2,phi,socpg)

       switch=(1.d0-dtanh(20.0d0*(rt(3)-1.6d0)))/2.0d0

       socg=socrg*switch+socpg*(1.0d0-switch)
       
      vgsoc=vg+socg/349.75d0

      vgsoc=min(100.0d0,vgsoc)
      
      return
      end subroutine fh2o_soc_nn
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine orderxnt(x1,x2)
      implicit none
      real*8 x1(3,4),x2(3,4),xvec(3,2),r1,r2,xt(3,4)
      xt=x1
      xvec(:,1)=xt(:,1)-xt(:,3)
      xvec(:,2)=xt(:,2)-xt(:,3)
      r1=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
      r2=dsqrt(dot_product(xvec(:,2),xvec(:,2)))

      if(r1.gt.r2)then
        x2(:,1)=xt(:,2)  !-->1st H closer to F
        x2(:,2)=xt(:,1)
        x2(:,3)=xt(:,3)
        x2(:,4)=xt(:,4)
      else
        x2=xt
      endif

      return
      end subroutine orderxnt

      real*8 function r1function(r,n)
      implicit none
      real*8, parameter:: r0=1.8d0,c=0.45d0
      real*8 r
      integer n

      if(n.eq.0)then
        r1function=1.0d0
        return
      else
        r1function=dtanh(dble(n)*c*(r-r0))
        return
      endif

      end function r1function
     
     
      real*8 function r2function(r,n)
      implicit none
      real*8,parameter:: r0=2.0d0,c=0.42d0
      real*8 r
      integer n

      if(n.eq.0)then
        r2function=1.0d0
        return
      else
        r2function=dtanh(dble(n)*c*(r-r0))
        return
      endif

      end function r2function
     
      real*8 function p2function(r,n)
      implicit none
      real*8,parameter:: r0=2.3d0,c=0.40d0
      real*8 r
      integer n

      if(n.eq.0)then
        p2function=1.0d0
        return
      else
        p2function=dtanh(dble(n)*c*(r-r0))
        return
      endif

      end function p2function
     

*********************SOC**************************
      subroutine soc_react_ground(rcom,t2,phi,socrg)
      use nnparam
      implicit none
      real*8 rcom,t2,phi,socrg,r1function
      real*8 tr1basis(nrbasissoc)
      integer j,n,k,l,nrcom

      tr1basis=0.0d0
      j=0
      do n=0,6
       do k=0,2
        do l=0,6
         
         j=j+1
         nrcom=n
         if(l.le.3)then
          tr1basis(j)=r1function(rcom,nrcom)*dcos(dble(k)*phi)*
     $                dcos(dble(l)*t2)
         else
          tr1basis(j)=r1function(rcom,nrcom)*dcos(dble(k)*phi)*
     $                 dsin(dble(l-3)*t2)
         endif
        enddo
       enddo
      enddo
      if(j.ne.nrbasissoc)then
       write(*,*)'error for nrbasissoc'
       stop
      endif
      socrg=dot_product(tr1basis(:),r1soccoef(:))

      return
      end subroutine soc_react_ground

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine soc_product_ground(rcom,t1,t2,phi,socpg)    
      use nnparam
      implicit none
      real*8 rcom,t1,t2,phi,socpg,socps,p2function
      real*8 tp2basis(npbasissoc)
      integer j,n,k,l,m,nrcom

      tp2basis=0.0d0
      j=0
      do n=0,5
       do k=0,2
        do l=0,3
         do m=0,3 
          j=j+1
          nrcom=n
          tp2basis(j)=p2function(rcom,nrcom)*dcos(dble(k)*phi)*
     $ dcos(dble(l)*t1)*dcos(dble(m)*t2)
         enddo
        enddo
       enddo
      enddo
      if(j.ne.npbasissoc)then
        write(*,*)'error for nrbasissoc'
        stop
      endif
       
      socps=dot_product(tp2basis(:),p2soccoef(:))
      socpg=-1.0d0*socps

      return
      end subroutine soc_product_ground

       subroutine fh2oNN(ct,vpes,vpesa,vpesb,vpesc) !ct in HHFO order and angstrom
       use nnparam
       implicit none
       integer,parameter::ndim=6
       integer j,k
       real*8 rb(ndim),xbond(ndim),basis(0:nbasis-1),tmp1
       real*8 txinput(1:nbasis-1)
       real*8 ct(3,4),xvec(3,ndim),xct(3,4),vpes
       real*8 vpesa,vpesb,vpesc,r12,r13,r14,r23,r24,r34

       basis=0.d0
       xct=ct
       
       xvec(:,1)=xct(:,2)-xct(:,1) !  H2->H1
       xvec(:,2)=xct(:,3)-xct(:,1) !  F->H1
       xvec(:,3)=xct(:,4)-xct(:,1) !  O->H1
       xvec(:,4)=xct(:,3)-xct(:,2) !  H2->F
       xvec(:,5)=xct(:,4)-xct(:,2) !  H2->O
       xvec(:,6)=xct(:,4)-xct(:,3) !  F-O 
       rb(1)=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
       rb(2)=dsqrt(dot_product(xvec(:,2),xvec(:,2)))
       rb(3)=dsqrt(dot_product(xvec(:,3),xvec(:,3)))
       rb(4)=dsqrt(dot_product(xvec(:,4),xvec(:,4)))
       rb(5)=dsqrt(dot_product(xvec(:,5),xvec(:,5)))
       rb(6)=dsqrt(dot_product(xvec(:,6),xvec(:,6))) 

       r12=rb(1);r13=rb(2);r14=rb(3);r23=rb(4);r24=rb(5);r34=rb(6)

       if(minval(rb(:)).le.0.67d0.or.minval(rb(:)).gt.1.45d0)then
         vpes=5.0d0
         return
       endif

       if(r34.le.1.5d0.or.r12.le.0.8d0.or.minval(rb).eq.r12)then
         vpes=5.0d0
         return
       endif
      
      if(r14.lt.r24.and.r24.gt.1.9d0) then
        if(r14.ge.1.8d0.or.r23.ge.1.8d0)then   
         vpes=5.0d0
         return
        endif
      endif

      if(r24.lt.r14.and.r14.gt.1.9d0) then
        if(r24.ge.1.8d0.or.r13.ge.1.8d0)then  
         vpes=5.0d0
         return
        endif
      endif

      if(r23.lt.r13.and.r13.gt.1.9d0.and.r24.gt.1.7d0) then
        if(r14.ge.1.8d0.or.r23.ge.1.8d0)then   
         vpes=5.0d0
         return
        endif
      endif

      if(r13.lt.r23.and.r23.gt.1.9d0.and.r14.gt.1.7d0) then
        if(r13.ge.1.8d0.or.r24.ge.1.8d0)then  
         vpes=5.0d0
         return
        endif
      endif

      xbond(:)=dexp(-rb(:)/alpha)
       
!      do j=1,nbasis
!       tmp1=0.0d0
!       do k=1,nterm(j)
!        tmp1=tmp1+
!    $  (xbond(1)**nindex(j,k,1))*(xbond(2)**nindex(j,k,2))
!    & *(xbond(3)**nindex(j,k,3))*(xbond(4)**nindex(j,k,4))
!    & *(xbond(5)**nindex(j,k,5))*(xbond(6)**nindex(j,k,6))
!       enddo
!       basis(j)=basis(j)+tmp1
!      enddo

       call bemsav(xbond,basis)
      
       do j=1,nbasis-1
        txinput(j)=basis(j)
       enddo

       call getpota(txinput,vpesa)
       call getpotb(txinput,vpesb)
       call getpotc(txinput,vpesc)
       vpes=(vpesa+vpesb+vpesc)/3.0d0

       if(vpes.lt.-1.5d0)vpes=5.0d0
       vpes=min(vpes,5.0d0)
       
       return
        
       end subroutine fh2oNN

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
        character f1*80,line
        character(len=1024) :: file_path1, file_path2

      file_path1 = trim(path)//"/FH2O/weights.txt"
      file_path2 = trim(path)//"/FH2O/biases.txt"

      open(4,file=file_path2,status='old')
      rewind(4)
      read(4,'(a100)')line
109   read(4,'(i4,3i3,3x,6i2)',end=121)ibasis,npd,nterm(ibasis+1),
     $iterm,(nindex(ibasis+1,iterm,ib), ib=1,6)
      if(nindex(ibasis+1,iterm,1).eq.2)then
        goto 121
      else
        goto 109
      endif
121   continue

      read(4,'(a100)')line
      do i=1,nrbasissoc
       read(4,*)r1soccoef(i)
      enddo

      read(4,'(a100)')line
      do i=1,npbasissoc
       read(4,*)p2soccoef(i)
      enddo
      read(4,*)line

        nfile=7
        open(nfile,file=file_path1)
        read(nfile,*)line
!       open(nfile,file='weights.txt-13')
!       open(4,file='biases.txt-13')
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

        read(4,*)line
        read(nfile,*)line
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
        allocate(weightb(nodemax,nodemax,2:nlayer),
     %   biasb(nodemax,2:nlayer))
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

        read(4,*)line
        read(nfile,*)line
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
        allocate(weightc(nodemax,nodemax,2:nlayer),
     % biasc(nodemax,2:nlayer))
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

        end subroutine pes_init

        subroutine getpota(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
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

        subroutine getpotb(x,vpot)
        use nnparam
        implicit none
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
c       write(*,*)ninput
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
     &*weightb(inode2,inode1,ilay1)
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
     &*weightb(inode2,inode1,ilay1)
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
        integer i,inode1,inode2,ilay1,ilay2
        real*8 x(ninput),y(nodemax,nlayer),vpot
        real*8, external :: tranfun
!-->....set up the normalized input layer
c       write(*,*)ninput
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
     &*weightc(inode2,inode1,ilay1)
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
     &*weightc(inode2,inode1,ilay1)
        enddo
!-->....the transfer function is linear y=x for output layer
!-->....so no operation is needed here
        enddo

!-->....the value of output layer is the fitted potntial 
        vpot=y(nodes(nlayer),nlayer)*pdelc(nscale)+pavgc(nscale)
        return
        end subroutine getpotc

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
 
      m(0)=1.D0
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
 
      subroutine evpoly(m,p)
      implicit none
      real*8,dimension(0:10),intent(in)::m
      real*8,dimension(0:17),intent(out)::p
 
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

      subroutine xyzTOjcb3(ct,rt)
      implicit none  
                     !   D  F
!           /|\         /
!     A  /\  |         /
!     H    \ |        /
!           \|       /        
!  ----------\------------------/ C  H
!             \ 
!              \  
!               \
!             O  B
!     
      real*8 ct(3,4),xAB(3),xABC(3),rt(6) !ct in HOH + F= HO + HF order
      real*8 xvec(3,6),crop1(3),crop2(3)
      real*8 mass(4)

      mass(1)=1.00783d0
      mass(2)=15.99491d0
      mass(3)=1.00783d0
      mass(4)=18.99840d0
      
      xAB(:)=(mass(1)*ct(:,1)+mass(2)*ct(:,2))/(mass(1)+mass(2))
      xABC(:)=(mass(1)*ct(:,1)+mass(2)*ct(:,2)+mass(3)*ct(:,3))/
     $ (mass(1)+mass(2)+mass(3))

      xvec(:,1)=ct(:,1)-ct(:,2)
      xvec(:,2)=ct(:,3)-xAB(:)
      xvec(:,3)=ct(:,4)-xABC(:)
      
      rt(1)=dsqrt(dot_product(xvec(:,1),xvec(:,1))) 
      rt(2)=dsqrt(dot_product(xvec(:,2),xvec(:,2))) 
      rt(3)=dsqrt(dot_product(xvec(:,3),xvec(:,3))) 

      rt(4)=dot_product(xvec(:,1),xvec(:,2))/rt(1)/rt(2)
      rt(5)=dot_product(xvec(:,2),xvec(:,3))/rt(2)/rt(3)

      rt(4)=dacos(max(-1.0d0,min(1.0d0,rt(4))))
      rt(5)=dacos(max(-1.0d0,min(1.0d0,rt(5))))
      
      call crossp(xvec(:,1),xvec(:,2),crop1) 
      call crossp(xvec(:,2),xvec(:,3),crop2) 
      
      rt(6)=dot_product(crop1,crop2)/dsqrt(dot_product(crop1,crop1))/
     $dsqrt(dot_product(crop2,crop2))
      
      rt(6)=dacos(max(-1.0d0,min(1.0d0,rt(6))))

      rt(6)=dacos(-1.d0)-rt(6) 

      return
      end subroutine xyzTOjcb3


!******************************************************************************
c     subroutine crossp(v1,v2,v3)
c     real(kind=8)::v1(1:3),v2(1:3),v3(1:3)
c     ! cross product of two 3-dimensional vectors
c     v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
c     v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
c     v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
c     return
c     end subroutine crossp

!******************************************************************************
      subroutine jcb3tocart(rint,cart)
      implicit none
      real*8 rint(6),cart(3,4)
      real*8 r1,r2,r3,th1,th2,phia,sth1,sth2,sth3,cth1,cth2,cth3
      real*8 wa,wb,wc,wd,abm,cdm,ama,amb,amc,amd
!-->-->rint(1-6)=r1,r2,r3,th1,th2,phi
!-->-->r1 is vector A->B, r2 is a vector pointing from the COM of AB to
!C
!-->-->r3 is a vector pointing from the COM of ABC to D
!-->-->th1 is the angle for r1 and r2
!-->-->th2 is the angle for r2 and r3
!-->-->phia is the diheral angle for r1*r2 plane and r2*r3 plane.
!-->-->cart(3,1-4)=cartesian coordinates of four atoms A,B,C,D.
      real*8 temp1,temp2,temp3

      ama=1.00783d0
      amb=15.99491d0
      amc=1.00783d0
      amd=18.99840d0

      r1=rint(1)
      r2=rint(2)
      r3=rint(3)
      th1=rint(4)
      th2=rint(5)
      phia=rint(6)
!-->-->ABC+D Jacobi coordinates 
!-->-->            Schematic diagram
c          X     
c         /|\                 - D  F
c          |                -   /
c          |              -    /
c     H A  |            -r3   /
c        \ |          -      /
c         \th1      -th2    /
c----------e---------------C-------------->Z 
c        r1|\      r2    /   H
c          | \         /
c          |  \      /
c                B O             
!c          |   
c               
!-->-->origin is the center of mass of AB.
!-->-->r2 lies on Z positive axis, AB is in the XZ positive plane
!-->-->Xa is always positive and Xb is always negative.
!-->-->Xc is always positive
      abm=r1/(ama+amb)
      wa=amb*abm
      wb=ama*abm
      sth1=dsin(th1)
      cth1=dcos(th1)
      sth2=dsin(th2)
      cth2=dcos(th2)
      cth3=dcos(phia)
      sth3=dsin(phia)
      cart(1,1)=wa*sth1
      cart(2,1)=0d0
      cart(3,1)=wa*cth1
      cart(1,2)=-wb*sth1
      cart(2,2)=0d0
      cart(3,2)=-wb*cth1
      cart(1,3)=0d0
      cart(2,3)=0d0
      cart(3,3)=r2
      cart(1,4)=r3*sth2*cth3
      cart(2,4)=r3*sth2*sth3
      cart(3,4)=r3*cth2+amc*r2/(ama+amb+amc)
      end subroutine jcb3tocart 

!******************************************************************************

      subroutine cart_diJacobi(c1,r)
! convert the diatom-diatom Jacobi coordinates to XYZ coordinates
!  c in angstrom and ABCD, r in angstrom and radians
      implicit none
      real*8 r(1:6),c(1:3,1:4),Mass(1:4),c1(1:3,1:4)
      real*8 crop1(1:3),crop2(1:3)
      real*8 vec(1:3,1:9),dotpp(1:9),bohr
      real*8 dotp,costh1,costh2,AMass,BMass,CMass,DMass,costh3

c     c(:,1)=c1(:,4)
c     c(:,2)=c1(:,1)
c     c(:,3)=c1(:,3)
c     c(:,4)=c1(:,2)
      c=c1  ! c in HOHF order angstrom
      bohr=1.0d0
      AMass=2.0141d0
      BMass=15.99491d0
      CMass=AMASS
      DMass=18.9984d0
      vec(:,1)=c(:,1)-c(:,2) !  vector B->A OH
      r(1)    =dsqrt( sum(vec(:,1)**2) )/bohr
      vec(:,2)=c(:,2)+vec(:,1)*AMass/(AMass+BMass) !vector BX (X COM of BA)
      vec(:,3)=c(:,3)-c(:,4) ! vector D->C OC
      r(2)    =dsqrt( sum(vec(:,3)**2) )/bohr
      vec(:,4)=c(:,4)+vec(:,3)*CMass/(CMass+DMass) !vector Y(Y COM of
CD)
      vec(:,5)=vec(:,4)-vec(:,2) !vecor R X->Y
      r(3)    =dsqrt( sum(vec(:,5)**2) )/bohr
!-
      vec(:,6)=c(:,3)-c(:,2) ! vector  O->C
      vec(:,7)=c(:,2)-c(:,3) ! vector  C->O
      vec(:,8)=c(:,4)-c(:,3) ! vector  C->O'
!     dompp(5)=dsqrt(dot_product(vec(:,5),vec(:,5)))
!     dotpp(6)=dsqrt(dot_product(vec(:,6),vec(:,6)))
!     dotpp(7)=dsqrt(dot_product(vec(:,7),vec(:,7)))

      dotpp(1)=dsqrt(dot_product(vec(:,1),vec(:,1)))
      dotpp(3)=dsqrt(dot_product(vec(:,3),vec(:,3)))
      dotpp(5)=dsqrt(dot_product(vec(:,5),vec(:,5)))

      if (dotpp(1)*dotpp(5).lt. 1.0d-16)then
        costh1=0.0d0
      else
        costh1=dot_product(vec(:,1),vec(:,5))/(dotpp(1)*dotpp(5))
      endif
      if (dotpp(3)*dotpp(5).lt. 1.0d-16)then
        costh2=0.0d0
      else
        costh2=dot_product(vec(:,3),vec(:,5))/(dotpp(3)*dotpp(5))
      endif
!  note the round off errors
      r(4) = dacos( max(-1.0d0,min(1.0d0,costh1))) ! in radians
      r(5) = dacos( max(-1.0d0,min(1.0d0,costh2))) ! in radians

      if (dabs(dcos(r(4))*dcos(r(5))).lt.1.0d-16)then
c       write(*,*) '   ******'
c       write(*,'(3f12.6)')c
c      write(*,*)'some 3-atom collinear, dihedral angle cant be defined'
c       write(*,*) '   ******'
!       stop
!       return
      else
! compute the Jacobi dihedral angle, based on the norm of the plane
       call crossp(vec(:,1), vec(:,5),crop1) ! crossproduct for the norm
       call crossp(vec(:,3), vec(:,5),crop2) !
       if(dsqrt( sum(crop1(:)**2))*dsqrt(sum(crop2(:)**2)).lt.1.0d-10)
     &  then
        r(6) = 0.0d0
       else
       costh3=dot_product(crop1(:),crop2(:))/(dsqrt(sum(crop1(:)**2))
     $ *dsqrt(sum(crop2(:)**2)))
!    write(99999,'(f20.10)') costh3
       r(6) = dacos( max(-1.0d0,min(1.0d0,costh3)))
       endif
      endif
      return
      end subroutine cart_diJacobi 
!******************************************************************************

      subroutine diJacobi_cart(r,c)
! convert the diatom-diatom Jacobi coordinates to XYZ coordinates
! r in angstrom and radians,  c in bohr and ABCD
      implicit none
      real*8 r(1:6),c(1:3,1:4),Mass(1:4)
      real*8 cm(3,2),vec(1:3,1:6), PI
      real*8 XA,XB,YC,YD,Cxsin,Dxsin,dotp
      real*8 RX1,RX2,RX3,RX4,RX5,RX6
      real*8 AMass,BMass,CMass,DMass

      PI=dacos(-1.0d0)
      
      AMass=1.00783d0
      BMass=15.99491d0
      CMass=AMASS
      DMass=18.99834d0

      Mass(1)=AMass
      Mass(2)=BMass
      Mass(3)=CMass
      Mass(4)=DMass
!     diatomic-diatomic AB+CD Jacobi is
!     r1, r2, R, theta1, theta2, phi
!     r1, distance AB
!     r2, distance CD
!     R,  distance between the centers of mass AB(X) and CD(Y)
!     theta1, angle from vector B->A to vector X->Y
!     theta2, angle from vector D->C to vector X->Y
!     phi, dihedral angle between the plane AXY and plane XYD    
!
!   A/- r1                D\
!     \                     \ r2
!      \ theta1    R         \
!      X\-------------------->\ Y--->  
!        \                     \ th2
!       B \                    -/ C

!   put ABY on xy plane, X origin point, X->Y (R)  on the x axis
!   and also put A in the positive y direction
!             ^y axis
!        A    |
!        /\-  |
!        r1\  |
!           \ |
!            \|     R
!    ---------X--------->Y-->x axis
!              \
!               \
!                B
      c=0.0d0
      XA=Bmass/(Amass+Bmass)*r(1)
      XB=Amass/(Amass+Bmass)*r(1)
!    A
      c(1,1)= XA*dcos(r(4))
      c(2,1)= XA*dsin(r(4))
      c(3,1)= 0.0d0
!    B
      c(1,2)= XB*dcos(PI-r(4))
      c(2,2)=-XB*dsin(PI-r(4))
      c(3,2)=0.d0

!     put ABY on xy plane, X origin point
!               ^y
!          A    |
!          /\-  |            
!          r1\  |               \D
!             \ |                \ |
!              \|     R           \|
!      ---------X----------------->Y-->x
!                \                 |\
!                 \                |-\/C
!                  B
!      C in +z and -z phase is identical, so put it in +z
      YC=Dmass/(Cmass+Dmass)*r(2)
      YD=Cmass/(Cmass+Dmass)*r(2)
      Cxsin = YC*dsin(r(5))
      Dxsin =-YD*dsin(r(5))
!     C
      c(1,3)=r(3)+YC*dcos(r(5))
      c(2,3)= Cxsin*dcos(r(6))
      c(3,3)= Cxsin*dsin(r(6))
!     D
      c(1,4)=r(3)-YD*dcos(r(5))
      c(2,4)= Dxsin*dcos(r(6))
      c(3,4)= Dxsin*dsin(r(6))

      return
      end subroutine diJacobi_cart 

!******************************************************************************
!******************************************************************************
      subroutine crossp(v1,v2,v3)
      real*8 v1(1:3),v2(1:3),v3(1:3)
      ! cross product of two 3-dimensional vectors
      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      return
      end subroutine crossp
