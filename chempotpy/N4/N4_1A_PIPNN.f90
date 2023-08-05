!**********************************************************************
!   System:                     N4
!   Functional form:            permutation-invariant polynomials
!   Common name:                N4(adiabatic ground state)
!   Number of derivatives:      1
!   Number of bodies:           4
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References:: J. Li, Z. Varga, D. G. Truhlar, and H. Guo
!                J. Chem. Theory Comput. 16, 4822 (2020).
!
!   Notes:    PES of N4 
!            - Neural network fit based on extended N4 
!              dataset with 21,406 points
!            - the N2 pairwise potentials are replaced compared
!              to PES_N4_singlet_umn_v3
!            - the gradient calculation is a built-in numerical one
!            - the code requires files weights.txt-Noconnected and 
!              biases.txt-Noconnected
!            - in the reference article this PES is called MB-PIP-NN
!
!     N1--N2
!
!     N3--N4
!   
!   Use:
!   subroutine n4pipNN(x,y,z,vpes,dedx,dedy,dedz,igrad)
!
!   Input: X(4),Y(4),Z(4)               in units of bohr
!   Input: igrad                        1 - gradient is calculated
!                                       0 - gradient is not calculated 
!   Output: Vpes                        in units of hartree
!   Output: dEdX(4),dEdY(4),dEdZ(4)     hartree/bohr
!**********************************************************************

!-->  program to get potential energy for a given geometry after NN fitting
!-->  global variables are declared in this module
       module nnparam
       implicit none
       real*8,parameter::alpha=1.0d0,PI=3.1415926d0,radian=PI/180.0d0
       integer,parameter::nbasis=306,ndim=6,natom=4,nnew=276
       integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
       integer nterm(1:nbasis),nindex(1:nbasis,1:1000,1:ndim)
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
      integer, parameter :: natoms=4
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: v
      double precision :: tx(natoms), ty(natoms), tz(natoms)
      double precision :: txg(natoms), tyg(natoms), tzg(natoms)
      integer :: iatom, idir, j, istate
      !initialize 
      p=0.d0
      g=0.d0
      d=0.d0

      ! input coordinates of tx, ty, tz are in bohr
      tx(:)=x(:,1)/0.529177211
      ty(:)=x(:,2)/0.529177211
      tz(:)=x(:,3)/0.529177211

      call pes_init(path)

      call n4pipNN(tx,ty,tz,v,txg,tyg,tzg,igrad)

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

      do istate=1,nstates
        p(istate)=v*27.211386
      enddo

      do istate=1,nstates
        g(istate,:,1)=txg(:)*51.422067
        g(istate,:,2)=tyg(:)*51.422067
        g(istate,:,3)=tzg(:)*51.422067
      enddo

      d=0.d0

      endsubroutine


!      subroutine n4pipNN(ct,vpes,vpesa,vpesb,vpesc,vn2) !ct in N4 order and angstrom
       subroutine n4pipNN(x,y,z,vpes,dedx,dedy,dedz,igrad)

       use nnparam
       implicit none
       integer i,j,k
       real*8 rb(ndim),xbond(ndim),basis(1:nbasis),tmp1,
     %  txinput(nbasis-1),bnew(1:nnew)
       real*8 ct(3,natom),xvec(3,ndim),cx(3,natom)
       real*8 vpes,vpesa,vpesb,vpesc
       real*8 vn2,vn2tmp,grad
       integer imol,igrad
! ZV
       real*8 x(4),y(4),z(4),dedx(4),dedy(4),dedz(4)
       real*8 dr, vp(6),vm(6),dedr(6), dVdx(12),dRdX(6,12)
       double precision Cconv
       double precision Econv
       double precision Gconv
       parameter(Cconv=0.52917721092d0)
       parameter(Econv=0.159360144d-2)
       parameter(Gconv=0.843297564d-3)

       basis=0.d0

! Cartesians in Ang
       do i=1,4
       cx(1,i)=x(i)*Cconv
       cx(2,i)=y(i)*Cconv
       cx(3,i)=z(i)*Cconv
!      write(1000,*) "cart ",cx(1,i),cx(2,i),cx(3,i)
       enddo

        dEdX(:)=0.0d0
        dEdY(:)=0.0d0
        dEdZ(:)=0.0d0

       k=0
       do i=1,natom-1
        do j=i+1,natom
         k=k+1
!        write(*,*)k,i,j
         xvec(:,k)=cx(:,i)-cx(:,j)
        enddo
       enddo
       if(k.ne.ndim)stop 'error in bond dimension'

       vn2=0d0
       do i=1,ndim
        rb(i)=dsqrt(dot_product(xvec(:,i),xvec(:,i)))
!       if(rb(i).le.0.6d0)then 
!         vpes=1000.0d0/23.0605d0
!         return
!       endif
!       write(1000,*) "rb(i) ",rb(i),i

        call ev2gm2(rb(i),vn2tmp,grad,1,0)
!       write(1000,*) "vn2tmp ",vn2tmp,i
        vn2=vn2+vn2tmp
       enddo

       vpes=0.0d0
       xbond(:)=dexp(-rb(:)/alpha)
       call bemsav49(xbond,basis)
       call bmx2b109(basis,bnew,nbasis,nnew)

       call getpota(bnew,vpesa)
       call getpotb(bnew,vpesb)
       call getpotc(bnew,vpesc)
       vpes=(vpesa+vpesb+vpesc)/3.0d0  !in eV
!      write(1000,*) "vmb ",vpes
! vpes in kcal/mol
       vpes=vpes*23.0605d0+vn2+457.4d0
!      write(1000,*) "vpes ",vpes
! vpes in au
       vpes=vpes*Econv

! numerical gradient w.r.t. dr, in kcal/Ang
       dr=1.0d-5 ! in Ang
       if (igrad.eq.1) then
! calculate evdrdx
       Call evdrdx(cx,rb,dRdX)

       do i=1,6
       vn2=0.0d0
       rb(i)=rb(i)+dr
        do j=1,6
        call ev2gm2(rb(j),vn2tmp,grad,1,0)
        vn2=vn2+vn2tmp
        enddo
!      write(*,*) "vn2p ",vn2,i
       xbond=0.0d0
       xbond(:)=dexp(-rb(:)/alpha)
       call bemsav49(xbond,basis)
       call bmx2b109(basis,bnew,nbasis,nnew)
       call getpota(bnew,vpesa)
       call getpotb(bnew,vpesb)
       call getpotc(bnew,vpesc)
       vp(i)=(vpesa+vpesb+vpesc)/3.0d0
! vp in kcal/mol
       vp(i)=vp(i)*23.0605d0+vn2+457.4d0
!       write(*,*) "vp(i) ",vp(i),i
       vn2=0.0d0
       rb(i)=rb(i)-2*dr
        do j=1,6
        call ev2gm2(rb(j),vn2tmp,grad,1,0)
        vn2=vn2+vn2tmp
        enddo
!      write(*,*) "vn2m ",vn2,i
       xbond=0.0d0
       xbond(:)=dexp(-rb(:)/alpha)
       call bemsav49(xbond,basis)
       call bmx2b109(basis,bnew,nbasis,nnew)
       call getpota(bnew,vpesa)
       call getpotb(bnew,vpesb)
       call getpotc(bnew,vpesc)
       vm(i)=(vpesa+vpesb+vpesc)/3.0d0
! vm in kcal/mol
       vm(i)=vm(i)*23.0605d0+vn2+457.4d0
!       write(*,*) "vm(i) ",vm(i),i
! dedr in kcal/(mol Ang)
       dedr(i)=(vp(i) - vm(i))/(2*dr)
!       write(*,*) "dedr(i) ",dedr(i),i
       rb(i)=rb(i)+dr
       enddo

      dVdX(:)=0.0d0

      do i=1,12
        do j=1,6
!         write(*,*) "dedr(j) ",dedr(j),j
!          write(*,*) "dRdX(j,i) ",dRdX(j,i),j,i
          dVdX(i)=dVdX(i) + dedr(j)*dRdX(j,i)
        enddo
!         write(*,*) "dVdX(i) ",dVdX(i)
      enddo

! dEdX, dEdY, and dEdZ in au/bohr
      do i=1,4
        dEdX(i)=dVdX(3*i-2)*Gconv
        dEdY(i)=dVdX(3*i-1)*Gconv
        dEdZ(i)=dVdX(3*i)*Gconv
!      write(*,*) "grad ",dEdX(i),dEdY(i),dEdZ(i)
      enddo
        endif

      if (igrad.eq.2) then 
        write (*,*) 'Only energy and gradient are available'
      endif

       return

       end subroutine n4pipNN

      subroutine EvdRdX(cx,r,dRdX)
!**********************************************************************
! Subroutine to evaluate dRdX for given R and X
! R:            R(6), 6 bond lengths
! X:            X(12), 12 Cartesian coordinates
!**********************************************************************

      integer i,j
      double precision cx(3,4),r(6),dRdX(6,12)

! Initialize dRdX(6,12)
      do i=1,6
        do j=1,12
          dRdX(i,j)=0.0d0
        enddo
      enddo

! Start to calculate the non-zero dRdX
! dr1dx
      dRdX(1,1)=(cx(1,1)-cx(1,2))/r(1)
      dRdX(1,2)=(cx(2,1)-cx(2,2))/r(1)
      dRdX(1,3)=(cx(3,1)-cx(3,2))/r(1)
      dRdX(1,4)=-dRdX(1,1)
      dRdX(1,5)=-dRdX(1,2)
      dRdX(1,6)=-dRdX(1,3)

! dr2dx
      dRdX(2,1)=(cx(1,1)-cx(1,3))/r(2)
      dRdX(2,2)=(cx(2,1)-cx(2,3))/r(2)
      dRdX(2,3)=(cx(3,1)-cx(3,3))/r(2)
      dRdX(2,7)=-dRdX(2,1)
      dRdX(2,8)=-dRdX(2,2)
      dRdX(2,9)=-dRdX(2,3)

! dr3dx
      dRdX(3,1)=(cx(1,1)-cx(1,4))/r(3)
      dRdX(3,2)=(cx(2,1)-cx(2,4))/r(3)
      dRdX(3,3)=(cx(3,1)-cx(3,4))/r(3)
      dRdX(3,10)=-dRdX(3,1)
      dRdX(3,11)=-dRdX(3,2)
      dRdX(3,12)=-dRdX(3,3)

! dr4dx
      dRdX(4,4)=(cx(1,2)-cx(1,3))/r(4)
      dRdX(4,5)=(cx(2,2)-cx(2,3))/r(4)
      dRdX(4,6)=(cx(3,2)-cx(3,3))/r(4)
      dRdX(4,7)=-dRdX(4,4)
      dRdX(4,8)=-dRdX(4,5)
      dRdX(4,9)=-dRdX(4,6)

! dr5dx
      dRdX(5,4)=(cx(1,2)-cx(1,4))/r(5)
      dRdX(5,5)=(cx(2,2)-cx(2,4))/r(5)
      dRdX(5,6)=(cx(3,2)-cx(3,4))/r(5)
      dRdX(5,10)=-dRdX(5,4)
      dRdX(5,11)=-dRdX(5,5)
      dRdX(5,12)=-dRdX(5,6)

! dr6dx
      dRdX(6,7)=(cx(1,3)-cx(1,4))/r(6)
      dRdX(6,8)=(cx(2,3)-cx(2,4))/r(6)
      dRdX(6,9)=(cx(3,3)-cx(3,4))/r(6)
      dRdX(6,10)=-dRdX(6,7)
      dRdX(6,11)=-dRdX(6,8)
      dRdX(6,12)=-dRdX(6,9)
! Finish the calculation of non-zero dRdX

!     do i=1,6
!       do j=1,12
!         write(*,*) "dRdX(i,j) ",dRdX(i,j),i,j
!       enddo
!     enddo

      return

      end subroutine EvdRdX
!**********************************************************************
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

      file_path1 = trim(path)//"/N4/weights.txt-Noconnected"
      file_path2 = trim(path)//"/N4/biases.txt-Noconnected"

      open(4,file=file_path2,status='old')
      rewind(4)
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

        read(4,'(a80)')line
        read(nfile,'(a80)')line
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
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
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

        read(4,'(a80)')line
        read(nfile,'(a80)')line
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
!-->....ifunc hence controls the type of transfer function used for hidden layers
!-->....At this time, only an equivalent transfer function can be used for all hidden layers
!-->....and the pure linear function is always applid to the output layer.
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


!     function emsav49(x,c) result(v)
!     implicit none
!     real*8,dimension(1:6)::x
!     real*8,dimension(0:305)::c
!     real*8::v
!     ! ::::::::::::::::::::
!     real*8,dimension(0:305)::p
!     call bemsav49(x,p)
!     v = dot_product(p,c)
!     return
!     end function emsav49

      subroutine bemsav49(x,p)
      implicit none
      real*8,dimension(1:6),intent(in)::x
      real*8,dimension(0:305),intent(out)::p
      ! ::::::::::::::::::::
      real*8,dimension(0:111)::m
      call evmono49(x,m)
      call evpoly49(m,p)
      return
      end subroutine bemsav49

      subroutine evmono49(x,m)
      implicit none
      real*8,dimension(1:6),intent(in)::x
      real*8,dimension(0:111),intent(out)::m

      m(0)=1.D0
      m(1)=x(6)
      m(2)=x(5)
      m(3)=x(4)
      m(4)=x(3)
      m(5)=x(2)
      m(6)=x(1)
      m(7)=m(3)*m(4)
      m(8)=m(2)*m(5)
      m(9)=m(1)*m(6)
      m(10)=m(1)*m(2)
      m(11)=m(1)*m(3)
      m(12)=m(2)*m(3)
      m(13)=m(1)*m(4)
      m(14)=m(2)*m(4)
      m(15)=m(1)*m(5)
      m(16)=m(3)*m(5)
      m(17)=m(4)*m(5)
      m(18)=m(2)*m(6)
      m(19)=m(3)*m(6)
      m(20)=m(4)*m(6)
      m(21)=m(5)*m(6)
      m(22)=m(1)*m(7)
      m(23)=m(2)*m(7)
      m(24)=m(1)*m(8)
      m(25)=m(2)*m(16)
      m(26)=m(2)*m(17)
      m(27)=m(3)*m(17)
      m(28)=m(1)*m(18)
      m(29)=m(1)*m(19)
      m(30)=m(1)*m(20)
      m(31)=m(3)*m(20)
      m(32)=m(1)*m(21)
      m(33)=m(2)*m(21)
      m(34)=m(1)*m(12)
      m(35)=m(1)*m(17)
      m(36)=m(2)*m(20)
      m(37)=m(3)*m(21)
      m(38)=m(1)*m(14)
      m(39)=m(1)*m(16)
      m(40)=m(2)*m(19)
      m(41)=m(4)*m(21)
      m(42)=m(2)*m(27)
      m(43)=m(1)*m(31)
      m(44)=m(1)*m(33)
      m(45)=m(1)*m(23)
      m(46)=m(1)*m(25)
      m(47)=m(1)*m(26)
      m(48)=m(1)*m(27)
      m(49)=m(1)*m(40)
      m(50)=m(1)*m(36)
      m(51)=m(2)*m(31)
      m(52)=m(1)*m(37)
      m(53)=m(2)*m(37)
      m(54)=m(1)*m(41)
      m(55)=m(2)*m(41)
      m(56)=m(3)*m(41)
      m(57)=m(1)*m(42)
      m(58)=m(1)*m(51)
      m(59)=m(1)*m(53)
      m(60)=m(1)*m(55)
      m(61)=m(1)*m(56)
      m(62)=m(2)*m(56)
      m(63)=m(1)*m(62)
      m(64)=m(2)*m(57)
      m(65)=m(3)*m(57)
      m(66)=m(4)*m(57)
      m(67)=m(5)*m(57)
      m(68)=m(1)*m(58)
      m(69)=m(3)*m(58)
      m(70)=m(4)*m(58)
      m(71)=m(1)*m(59)
      m(72)=m(2)*m(59)
      m(73)=m(1)*m(60)
      m(74)=m(2)*m(60)
      m(75)=m(1)*m(61)
      m(76)=m(2)*m(62)
      m(77)=m(3)*m(61)
      m(78)=m(3)*m(62)
      m(79)=m(4)*m(61)
      m(80)=m(4)*m(62)
      m(81)=m(5)*m(59)
      m(82)=m(5)*m(60)
      m(83)=m(5)*m(62)
      m(84)=m(6)*m(58)
      m(85)=m(6)*m(59)
      m(86)=m(6)*m(60)
      m(87)=m(6)*m(61)
      m(88)=m(2)*m(64)
      m(89)=m(3)*m(65)
      m(90)=m(4)*m(66)
      m(91)=m(5)*m(67)
      m(92)=m(1)*m(68)
      m(93)=m(3)*m(69)
      m(94)=m(4)*m(70)
      m(95)=m(1)*m(71)
      m(96)=m(2)*m(72)
      m(97)=m(1)*m(73)
      m(98)=m(2)*m(74)
      m(99)=m(1)*m(75)
      m(100)=m(2)*m(76)
      m(101)=m(3)*m(77)
      m(102)=m(3)*m(78)
      m(103)=m(4)*m(79)
      m(104)=m(4)*m(80)
      m(105)=m(5)*m(81)
      m(106)=m(5)*m(82)
      m(107)=m(5)*m(83)
      m(108)=m(6)*m(84)
      m(109)=m(6)*m(85)
      m(110)=m(6)*m(86)
      m(111)=m(6)*m(87)

      return
      end subroutine evmono49

      subroutine evpoly49(m,p)
      implicit none
      real*8,dimension(0:111),intent(in)::m
      real*8,dimension(0:305),intent(out)::p

      p(0)=m(0)
      p(1)=m(1)+m(2)+m(3)+m(4)+m(5)+m(6)
      p(2)=m(7)+m(8)+m(9)
      p(3)=m(10)+m(11)+m(12)+m(13)+m(14)+m(15)+m(16)+m(17)+m(18)+m(19)+m
     &(20)+m(21)
      p(4)=p(1)*p(1)-p(3)-p(2)-p(3)-p(2)
      p(5)=m(22)+m(23)+m(24)+m(25)+m(26)+m(27)+m(28)+m(29)+m(30)+m(31)+m
     &(32)+m(33)
      p(6)=m(34)+m(35)+m(36)+m(37)
      p(7)=m(38)+m(39)+m(40)+m(41)
      p(8)=p(1)*p(2)-p(5)
      p(9)=p(1)*p(3)-p(6)-p(7)-p(5)-p(6)-p(7)-p(5)-p(6)-p(7)
      p(10)=p(1)*p(4)-p(9)-p(8)
      p(11)=m(42)+m(43)+m(44)
      p(12)=m(45)+m(46)+m(47)+m(48)+m(49)+m(50)+m(51)+m(52)+m(53)+m(54)+
     &m(55)+m(56)
      p(13)=p(2)*p(3)-p(12)
      p(14)=p(1)*p(5)-p(12)-p(11)-p(13)-p(12)-p(11)-p(11)-p(11)
      p(15)=p(1)*p(6)-p(12)
      p(16)=p(1)*p(7)-p(12)
      p(17)=p(2)*p(2)-p(11)-p(11)
      p(18)=p(3)*p(3)-p(12)-p(11)-p(15)-p(16)-p(14)-p(12)-p(11)-p(15)-p(
     &16)-p(14)-p(12)-p(11)-p(12)-p(11)
      p(19)=p(2)*p(4)-p(14)
      p(20)=p(3)*p(4)-p(15)-p(16)-p(13)
      p(21)=p(1)*p(10)-p(20)-p(19)
      p(22)=m(57)+m(58)+m(59)+m(60)+m(61)+m(62)
      p(23)=p(1)*p(11)-p(22)
      p(24)=p(2)*p(6)
      p(25)=p(2)*p(7)
      p(26)=p(1)*p(12)-p(22)-p(24)-p(25)-p(22)-p(22)-p(22)
      p(27)=p(2)*p(5)-p(22)-p(23)-p(22)
      p(28)=p(3)*p(5)-p(22)-p(26)-p(24)-p(25)-p(23)-p(22)-p(24)-p(25)-p(
     &23)-p(22)-p(22)
      p(29)=p(3)*p(6)-p(22)-p(26)-p(22)
      p(30)=p(3)*p(7)-p(22)-p(26)-p(22)
      p(31)=p(2)*p(9)-p(26)-p(28)
      p(32)=p(1)*p(14)-p(26)-p(23)-p(28)
      p(33)=p(4)*p(6)-p(25)
      p(34)=p(4)*p(7)-p(24)
      p(35)=p(1)*p(17)-p(27)
      p(36)=p(1)*p(18)-p(29)-p(30)-p(28)
      p(37)=p(2)*p(10)-p(32)
      p(38)=p(3)*p(10)-p(33)-p(34)-p(31)
      p(39)=p(1)*p(21)-p(38)-p(37)
      p(40)=m(63)
      p(41)=m(64)+m(65)+m(66)+m(67)+m(68)+m(69)+m(70)+m(71)+m(72)+m(73)+
     &m(74)+m(75)+m(76)+m(77)+m(78)+m(79)+m(80)+m(81)+m(82)+m(83)+m(84)+
     &m(85)+m(86)+m(87)
      p(42)=p(1)*p(22)-p(40)-p(41)-p(40)-p(40)-p(40)-p(40)-p(40)
      p(43)=p(2)*p(11)-p(40)-p(40)-p(40)
      p(44)=p(2)*p(12)-p(41)
      p(45)=p(3)*p(11)-p(41)
      p(46)=p(5)*p(6)-p(41)
      p(47)=p(5)*p(7)-p(41)
      p(48)=p(6)*p(7)-p(40)-p(40)-p(40)-p(40)
      p(49)=p(4)*p(11)-p(42)
      p(50)=p(2)*p(15)-p(46)
      p(51)=p(2)*p(16)-p(47)
      p(52)=p(4)*p(12)-p(41)-p(50)-p(51)
      p(53)=p(2)*p(14)-p(42)-p(49)-p(42)
      p(54)=p(6)*p(6)-p(42)-p(42)
      p(55)=p(7)*p(7)-p(42)-p(42)
      p(56)=p(3)*p(17)-p(44)
      p(57)=p(2)*p(18)-p(48)
      p(58)=p(3)*p(14)-p(41)-p(52)-p(46)-p(47)-p(45)-p(45)
      p(59)=p(6)*p(9)-p(41)-p(52)-p(47)
      p(60)=p(7)*p(9)-p(41)-p(52)-p(46)
      p(61)=p(2)*p(20)-p(52)-p(58)
      p(62)=p(1)*p(32)-p(52)-p(49)-p(58)
      p(63)=p(6)*p(10)-p(51)
      p(64)=p(7)*p(10)-p(50)
      p(65)=p(2)*p(17)-p(43)
      p(66)=p(3)*p(18)-p(46)-p(47)-p(45)-p(59)-p(60)-p(58)
      p(67)=p(2)*p(19)-p(49)
      p(68)=p(1)*p(36)-p(59)-p(60)-p(58)-p(57)-p(66)-p(66)
      p(69)=p(2)*p(21)-p(62)
      p(70)=p(3)*p(21)-p(63)-p(64)-p(61)
      p(71)=p(1)*p(39)-p(70)-p(69)
      p(72)=p(40)*p(1)
      p(73)=p(2)*p(22)-p(72)
      p(74)=p(6)*p(11)
      p(75)=p(7)*p(11)
      p(76)=p(3)*p(22)-p(72)-p(74)-p(75)-p(72)-p(72)-p(72)
      p(77)=m(88)+m(89)+m(90)+m(91)+m(92)+m(93)+m(94)+m(95)+m(96)+m(97)+
     &m(98)+m(99)+m(100)+m(101)+m(102)+m(103)+m(104)+m(105)+m(106)+m(107
     &)+m(108)+m(109)+m(110)+m(111)
      p(78)=p(1)*p(42)-p(72)-p(76)
      p(79)=p(5)*p(11)-p(72)-p(73)-p(72)
      p(80)=p(2)*p(26)-p(76)-p(77)
      p(81)=p(6)*p(12)-p(72)-p(76)-p(72)
      p(82)=p(7)*p(12)-p(72)-p(76)-p(72)
      p(83)=p(8)*p(11)-p(72)
      p(84)=p(6)*p(17)
      p(85)=p(7)*p(17)
      p(86)=p(9)*p(11)-p(76)-p(77)
      p(87)=p(2)*p(29)-p(81)
      p(88)=p(6)*p(14)-p(75)-p(75)
      p(89)=p(2)*p(30)-p(82)
      p(90)=p(7)*p(14)-p(74)-p(74)
      p(91)=p(1)*p(48)-p(76)-p(81)-p(82)
      p(92)=p(10)*p(11)-p(78)
      p(93)=p(2)*p(33)-p(88)
      p(94)=p(2)*p(34)-p(90)
      p(95)=p(10)*p(12)-p(77)-p(93)-p(94)
      p(96)=p(2)*p(27)-p(73)-p(79)
      p(97)=p(2)*p(28)-p(76)-p(86)
      p(98)=p(1)*p(53)-p(80)-p(79)-p(97)
      p(99)=p(1)*p(54)-p(81)
      p(100)=p(1)*p(55)-p(82)
      p(101)=p(5)*p(18)-p(76)-p(91)-p(87)-p(89)-p(86)
      p(102)=p(6)*p(18)-p(74)-p(90)
      p(103)=p(7)*p(18)-p(75)-p(88)
      p(104)=p(2)*p(31)-p(77)-p(86)
      p(105)=p(2)*p(36)-p(91)-p(101)
      p(106)=p(3)*p(32)-p(77)-p(95)-p(88)-p(90)-p(86)
      p(107)=p(4)*p(29)-p(82)-p(80)-p(99)
      p(108)=p(4)*p(30)-p(81)-p(80)-p(100)
      p(109)=p(2)*p(38)-p(95)-p(106)
      p(110)=p(1)*p(62)-p(95)-p(92)-p(106)
      p(111)=p(6)*p(21)-p(94)
      p(112)=p(7)*p(21)-p(93)
      p(113)=p(1)*p(65)-p(96)
      p(114)=p(1)*p(66)-p(102)-p(103)-p(101)
      p(115)=p(2)*p(37)-p(92)
      p(116)=p(10)*p(18)-p(99)-p(100)-p(97)
      p(117)=p(2)*p(39)-p(110)
      p(118)=p(3)*p(39)-p(111)-p(112)-p(109)
      p(119)=p(1)*p(71)-p(118)-p(117)
      p(120)=p(40)*p(2)
      p(121)=p(40)*p(3)
      p(122)=p(40)*p(4)
      p(123)=p(11)*p(12)-p(121)
      p(124)=p(2)*p(42)-p(122)
      p(125)=p(6)*p(22)-p(121)
      p(126)=p(7)*p(22)-p(121)
      p(127)=p(2)*p(41)-p(121)-p(123)-p(121)
      p(128)=p(6)*p(23)-p(123)
      p(129)=p(7)*p(23)-p(123)
      p(130)=p(3)*p(41)-p(122)-p(121)-p(120)-p(125)-p(128)-p(126)-p(129)
     &-p(124)-p(123)-p(122)-p(121)-p(120)-p(125)-p(126)-p(124)-p(123)-p(
     &122)-p(121)-p(120)-p(122)-p(121)-p(120)-p(120)-p(120)-p(120)-p(120
     &)
      p(131)=p(3)*p(42)-p(121)-p(125)-p(126)-p(121)
      p(132)=p(4)*p(41)-p(121)-p(131)-p(128)-p(129)-p(127)-p(121)
      p(133)=p(1)*p(78)-p(122)-p(131)
      p(134)=p(11)*p(11)-p(120)-p(120)
      p(135)=p(2)*p(48)-p(130)
      p(136)=p(11)*p(17)-p(120)
      p(137)=p(2)*p(44)-p(123)
      p(138)=p(2)*p(45)-p(121)
      p(139)=p(11)*p(14)-p(122)-p(124)-p(122)
      p(140)=p(6)*p(27)-p(127)
      p(141)=p(2)*p(54)
      p(142)=p(7)*p(27)-p(127)
      p(143)=p(2)*p(52)-p(131)-p(132)
      p(144)=p(1)*p(81)-p(125)-p(141)-p(135)-p(125)
      p(145)=p(2)*p(55)
      p(146)=p(1)*p(82)-p(126)-p(135)-p(145)-p(126)
      p(147)=p(11)*p(18)-p(130)
      p(148)=p(6)*p(28)-p(129)-p(123)-p(143)
      p(149)=p(7)*p(28)-p(128)-p(123)-p(143)
      p(150)=p(6)*p(30)-p(121)-p(146)
      p(151)=p(11)*p(19)-p(122)
      p(152)=p(2)*p(50)-p(128)
      p(153)=p(2)*p(51)-p(129)
      p(154)=p(11)*p(20)-p(131)-p(132)
      p(155)=p(2)*p(59)-p(144)-p(148)
      p(156)=p(6)*p(32)-p(129)
      p(157)=p(2)*p(60)-p(146)-p(149)
      p(158)=p(7)*p(32)-p(128)
      p(159)=p(6)*p(34)-p(122)-p(145)-p(122)
      p(160)=p(11)*p(21)-p(133)
      p(161)=p(2)*p(63)-p(156)
      p(162)=p(2)*p(64)-p(158)
      p(163)=p(12)*p(21)-p(132)-p(161)-p(162)
      p(164)=p(2)*p(53)-p(124)-p(139)
      p(165)=p(2)*p(58)-p(131)-p(154)
      p(166)=p(3)*p(54)-p(125)-p(144)
      p(167)=p(3)*p(55)-p(126)-p(146)
      p(168)=p(3)*p(65)-p(137)
      p(169)=p(17)*p(18)-p(135)
      p(170)=p(1)*p(98)-p(143)-p(139)-p(165)
      p(171)=p(4)*p(54)-p(135)
      p(172)=p(4)*p(55)-p(135)
      p(173)=p(2)*p(66)-p(150)
      p(174)=p(1)*p(101)-p(148)-p(149)-p(147)-p(173)-p(165)-p(147)
      p(175)=p(1)*p(102)-p(150)-p(148)-p(166)
      p(176)=p(1)*p(103)-p(150)-p(149)-p(167)
      p(177)=p(2)*p(61)-p(132)-p(154)
      p(178)=p(2)*p(68)-p(159)-p(174)
      p(179)=p(3)*p(62)-p(132)-p(163)-p(156)-p(158)-p(154)
      p(180)=p(6)*p(38)-p(132)-p(163)-p(157)
      p(181)=p(7)*p(38)-p(132)-p(163)-p(155)
      p(182)=p(2)*p(70)-p(163)-p(179)
      p(183)=p(1)*p(110)-p(163)-p(160)-p(179)
      p(184)=p(6)*p(39)-p(162)
      p(185)=p(7)*p(39)-p(161)
      p(186)=p(2)*p(65)-p(136)
      p(187)=p(3)*p(66)-p(148)-p(149)-p(147)-p(175)-p(176)-p(174)
      p(188)=p(2)*p(67)-p(151)
      p(189)=p(4)*p(66)-p(166)-p(167)-p(165)
      p(190)=p(2)*p(69)-p(160)
      p(191)=p(18)*p(21)-p(171)-p(172)-p(169)
      p(192)=p(2)*p(71)-p(183)
      p(193)=p(3)*p(71)-p(184)-p(185)-p(182)
      p(194)=p(1)*p(119)-p(193)-p(192)
      p(195)=p(40)*p(5)
      p(196)=p(40)*p(6)
      p(197)=p(40)*p(7)
      p(198)=p(40)*p(8)
      p(199)=p(40)*p(9)
      p(200)=p(40)*p(10)
      p(201)=p(11)*p(22)-p(195)
      p(202)=p(12)*p(22)-p(196)-p(197)-p(195)-p(196)-p(197)-p(195)-p(196
     &)-p(197)
      p(203)=p(17)*p(22)-p(198)
      p(204)=p(6)*p(43)
      p(205)=p(7)*p(43)
      p(206)=p(11)*p(26)-p(199)-p(202)
      p(207)=p(2)*p(76)-p(199)-p(202)
      p(208)=p(2)*p(78)-p(200)
      p(209)=p(6)*p(41)-p(199)-p(195)-p(202)-p(195)
      p(210)=p(6)*p(42)-p(197)-p(197)-p(197)
      p(211)=p(7)*p(41)-p(199)-p(195)-p(202)-p(195)
      p(212)=p(7)*p(42)-p(196)-p(196)-p(196)
      p(213)=p(11)*p(29)-p(209)
      p(214)=p(11)*p(30)-p(211)
      p(215)=p(18)*p(22)-p(199)-p(213)-p(214)
      p(216)=p(2)*p(77)-p(199)-p(206)
      p(217)=p(6)*p(49)-p(205)
      p(218)=p(7)*p(49)-p(204)
      p(219)=p(3)*p(77)-p(200)-p(199)-p(198)-p(209)-p(217)-p(211)-p(218)
     &-p(207)-p(204)-p(205)-p(200)-p(199)-p(198)-p(200)-p(198)-p(200)-p(
     &198)
      p(220)=p(3)*p(78)-p(199)-p(210)-p(212)
      p(221)=p(10)*p(41)-p(199)-p(220)-p(217)-p(218)-p(216)
      p(222)=p(1)*p(133)-p(200)-p(220)
      p(223)=p(11)*p(27)-p(195)-p(203)
      p(224)=p(1)*p(134)-p(201)
      p(225)=p(2)*p(81)-p(209)
      p(226)=p(2)*p(80)-p(202)-p(206)
      p(227)=p(2)*p(82)-p(211)
      p(228)=p(2)*p(91)-p(215)-p(219)
      p(229)=p(11)*p(28)-p(199)-p(207)
      p(230)=p(6)*p(53)-p(205)
      p(231)=p(5)*p(54)-p(209)
      p(232)=p(7)*p(53)-p(204)
      p(233)=p(7)*p(54)-p(196)
      p(234)=p(5)*p(55)-p(211)
      p(235)=p(6)*p(55)-p(197)
      p(236)=p(11)*p(35)-p(198)
      p(237)=p(6)*p(65)
      p(238)=p(7)*p(65)
      p(239)=p(2)*p(86)-p(199)-p(229)
      p(240)=p(1)*p(139)-p(206)-p(229)-p(224)
      p(241)=p(17)*p(29)-p(225)
      p(242)=p(2)*p(99)-p(231)
      p(243)=p(17)*p(30)-p(227)
      p(244)=p(2)*p(95)-p(220)-p(221)
      p(245)=p(4)*p(81)-p(202)-p(242)-p(227)
      p(246)=p(2)*p(100)-p(234)
      p(247)=p(4)*p(82)-p(202)-p(225)-p(246)
      p(248)=p(11)*p(36)-p(215)-p(219)
      p(249)=p(2)*p(102)-p(233)
      p(250)=p(6)*p(58)-p(206)-p(214)-p(244)-p(214)
      p(251)=p(2)*p(103)-p(235)
      p(252)=p(7)*p(58)-p(213)-p(206)-p(244)-p(213)
      p(253)=p(1)*p(150)-p(215)-p(233)-p(235)
      p(254)=p(11)*p(37)-p(200)
      p(255)=p(2)*p(93)-p(217)
      p(256)=p(2)*p(94)-p(218)
      p(257)=p(11)*p(38)-p(220)-p(221)
      p(258)=p(2)*p(107)-p(245)-p(250)
      p(259)=p(6)*p(62)-p(218)
      p(260)=p(2)*p(108)-p(247)-p(252)
      p(261)=p(7)*p(62)-p(217)
      p(262)=p(6)*p(64)-p(200)-p(246)-p(200)
      p(263)=p(11)*p(39)-p(222)
      p(264)=p(2)*p(111)-p(259)
      p(265)=p(2)*p(112)-p(261)
      p(266)=p(12)*p(39)-p(221)-p(264)-p(265)
      p(267)=p(2)*p(98)-p(208)-p(240)
      p(268)=p(6)*p(54)-p(210)
      p(269)=p(7)*p(55)-p(212)
      p(270)=p(2)*p(96)-p(203)-p(223)
      p(271)=p(2)*p(97)-p(207)-p(229)
      p(272)=p(2)*p(101)-p(215)-p(248)
      p(273)=p(2)*p(106)-p(220)-p(257)
      p(274)=p(6)*p(59)-p(220)-p(215)-p(209)
      p(275)=p(7)*p(60)-p(220)-p(215)-p(211)
      p(276)=p(5)*p(66)-p(215)-p(253)-p(249)-p(251)-p(248)
      p(277)=p(6)*p(66)-p(213)-p(252)
      p(278)=p(7)*p(66)-p(214)-p(250)
      p(279)=p(2)*p(104)-p(216)-p(239)
      p(280)=p(2)*p(105)-p(219)-p(248)
      p(281)=p(1)*p(170)-p(244)-p(240)-p(273)
      p(282)=p(10)*p(54)-p(227)
      p(283)=p(10)*p(55)-p(225)
      p(284)=p(2)*p(114)-p(253)-p(276)
      p(285)=p(1)*p(174)-p(250)-p(252)-p(248)-p(276)-p(273)
      p(286)=p(6)*p(68)-p(217)-p(261)-p(251)
      p(287)=p(7)*p(68)-p(218)-p(259)-p(249)
      p(288)=p(2)*p(109)-p(221)-p(257)
      p(289)=p(2)*p(116)-p(262)-p(285)
      p(290)=p(3)*p(110)-p(221)-p(266)-p(259)-p(261)-p(257)
      p(291)=p(6)*p(70)-p(221)-p(266)-p(260)
      p(292)=p(7)*p(70)-p(221)-p(266)-p(258)
      p(293)=p(2)*p(118)-p(266)-p(290)
      p(294)=p(1)*p(183)-p(266)-p(263)-p(290)
      p(295)=p(6)*p(71)-p(265)
      p(296)=p(7)*p(71)-p(264)
      p(297)=p(1)*p(186)-p(270)
      p(298)=p(1)*p(187)-p(277)-p(278)-p(276)
      p(299)=p(2)*p(115)-p(254)
      p(300)=p(1)*p(189)-p(286)-p(287)-p(285)-p(284)-p(298)
      p(301)=p(2)*p(117)-p(263)
      p(302)=p(18)*p(39)-p(282)-p(283)-p(280)
      p(303)=p(2)*p(119)-p(294)
      p(304)=p(3)*p(119)-p(295)-p(296)-p(293)
      p(305)=p(1)*p(194)-p(304)-p(303)

      return
      end subroutine evpoly49



      subroutine ev2gm2(r,v,grad,imol,igrad)
!**********************************************************************
!
! Compute the diatomic potential of ground-state N2 (singlet)
!
! Input:  r      interatomic distance in Angstrom
!         imol   1 - stands for N2 for N4 PES
!         igrad  1 - gradient is calculated
!                0 - gradient is not calculated
! Output: V      potential in kcal/mol
!         grad   gradient (kcal/mol)/Angstrom
!
!**********************************************************************
      implicit none

      integer,intent(in) :: imol, igrad
      double precision :: r
      double precision :: v, grad
      double precision :: re, de, as(0:6)
      double precision :: y, fy, u
      double precision :: dfdy,rr,dydr,dfdr
      integer :: k
! Dispersion variables
      double precision :: dist(1),disp,dispdr(1)

! Original parameters
!       re=1.098d0
!       de=228.7d0
!       as(0) =  2.70963254293d0
!       as(1) =  1.32620177271d-1
!       as(2) =  2.96757048793d-1
!       as(3) =  1.97112432229d-1
!       as(4) = -5.02002309588d-1
!       as(5) =  3.80734244606d-1
!       as(6) =  1.21001628750d0

! Modified parameters for D3(BJ)
        re=1.098d0
        de=225.213d0
        as(0) =  2.7475450369759d0
        as(1) =  0.218868498415108d0
        as(2) =  0.248885765371433d0
        as(3) = -0.229295466336412d0
        as(4) = -0.653389048592838d0
        as(5) =  1.03611964035396d0
        as(6) =  1.71287482791961d0

      v=0.d0

      y=(r*r*r*r - re*re*re*re)/(r*r*r*r + re*re*re*re)

      fy = as(0) + as(1)*y + as(2)*y*y + as(3)*y*y*y + as(4)*y*y*y*y 
     &   + as(5)*y*y*y*y*y + as(6)*y*y*y*y*y*y

      u=exp(-fy*(r-re))

      v=de*(1.0d0-u)*(1.0d0-u)-de


! Add D3 dispersion correction
        dist(1)=r
        call d3disp(dist,disp,dispdr,0,imol)
        v=v+disp

! Compute the gradient if needed
      if (igrad.eq.1) then
       grad=0.d0

       dfdy = as(1) + 2.0d0*as(2)*y + 3.0d0*as(3)*y*y 
     &      + 4.0d0*as(4)*y*y*y + 5.0d0*as(5)*y*y*y*y 
     &      + 6.0d0*as(6)*y*y*y*y*y
        rr = r*r*r*r + re*re*re*re
        dydr = 8.0d0*r*r*r*re*re*re*re/(rr*rr)
        dfdr = dfdy*dydr
        grad = 2.0d0*de*(1-u)*u*(dfdr*(r-re)+fy)

! Add analytical gradient of D3 dispersion correction
        call d3disp(dist,disp,dispdr,1,imol)
        grad= grad + dispdr(1)

      endif
      return
      end subroutine ev2gm2

      subroutine d3disp(dist,disp,dispdr,igrad,imol)
!**********************************************************************
! Dispersion correction based on Grimme's D3(BJ) calculation for
! diatomic pairs
!
! Several subroutines of DFTD3 V3.1 Rev 1 by Grimme were merged into
! subroutine edisp and they have been heavily modified to calculate
! only dispersion energy corrections that are needed.
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
! and
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011),
! 1456-1465
!
! The C6 values are fixed.
!
!**********************************************************************

      double precision cn(2),s6,s8,rs6,rs8
      double precision dist(1), e6(1), e8(1), disp, dispdr(1), c6(2)
      double precision e6dr(1),e8dr(1)
      integer iz(2), mxc(94), i, j, igrad
      double precision c6ab(94,94,5,5,3)
      double precision r2r4(94)
      double precision autoang,autokcal
      integer imol

      autoang =0.52917726d0
      autokcal=627.509541d0

! Generalized parameters for BJ damping from P. Verma, B. Wang,
! L. E. Fernandez, and D. G. Truhlar, J. Phys. Chem. A 121, 2855 (2017)
      s6= 1.0d0
      s8= 2.0d0
      rs6= 0.5299d0
      rs8= 2.20d0

      do i=1,1
      dist(i)=dist(i)/autoang
      enddo

      if (imol.eq.1) then
! iz for N2 system
      iz(1)=7
      iz(2)=7

! C6 for N2 system
      c6(1)=19.7d0
      c6(2)=19.7d0
      else
! currently imol = 1 is used
      stop
      endif

! Calculate dispersion correction
      call edisp(94,5,2,dist,iz,mxc, 
     &     rs6,rs8,e6,e8,e6dr,e8dr,c6,0)

      disp = 0.0d0

      do i=1,1
      disp =disp + (-s6*e6(i)-s8*e8(i))*autokcal
      enddo

      if (igrad .eq. 1) then
      call edisp(94,5,2,dist,iz,mxc, 
     &     rs6,rs8,e6,e8,e6dr,e8dr,c6,1)

      dispdr(:) = 0.0d0

      do i=1,1
      dispdr(i) =dispdr(i) + (-s6*e6dr(i)-s8*e8dr(i))*autokcal/autoang
      enddo
      endif

      do i=1,1
      dist(i)=dist(i)*autoang
      enddo

      end subroutine d3disp

!**********************************************************************
! compute energy
!**********************************************************************
      subroutine edisp(max_elem,maxc,n,dist,iz,mxc, 
     &           rs6,rs8,e6,e8,e6dr,e8dr,c6a,igrad)

      integer n,iz(2),max_elem,maxc,mxc(max_elem)
      double precision dist(1),r2r4(max_elem),r0ab(max_elem,max_elem)
      double precision rs6,rs8,rcov(max_elem)
      double precision c6ab(max_elem,max_elem,maxc,maxc,3)
      double precision e6(1), e8(1), c6a(2), e6dr(1), e8dr(1)

      integer iat,jat,igrad
      double precision r,tmp,c6,c8,a1,a2
      double precision damp6,damp8
      double precision cn(n)
      double precision r2ab(n*n),cc6ab(n*n),dmp(n*n)
      integer step

      e6(:) =0.0d0
      e8(:) =0.0d0

      e6dr(:) =0.0d0
      e8dr(:) =0.0d0

      a1=rs6
      a2=rs8

!  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
!  the large number of digits is just to keep the results consistent
!  with older versions. They should not imply any higher accuracy than
!  the old values
      r2r4(1:94)=(/ 
     &2.00734898d0,  1.56637132d0,  5.01986934d0,  3.85379032d0,  
     &3.64446594d0,  3.10492822d0,  2.71175247d0,  2.59361680d0,  
     &2.38825250d0,  2.21522516d0,  6.58585536d0,  5.46295967d0,  
     &5.65216669d0,  4.88284902d0,  4.29727576d0,  4.04108902d0,  
     &3.72932356d0,  3.44677275d0,  7.97762753d0,  7.07623947d0,  
     &6.60844053d0,  6.28791364d0,  6.07728703d0,  5.54643096d0,  
     &5.80491167d0,  5.58415602d0,  5.41374528d0,  5.28497229d0,  
     &5.22592821d0,  5.09817141d0,  6.12149689d0,  5.54083734d0,  
     &5.06696878d0,  4.87005108d0,  4.59089647d0,  4.31176304d0,  
     &9.55461698d0,  8.67396077d0,  7.97210197d0,  7.43439917d0,  
     &6.58711862d0,  6.19536215d0,  6.01517290d0,  5.81623410d0,  
     &5.65710424d0,  5.52640661d0,  5.44263305d0,  5.58285373d0,  
     &7.02081898d0,  6.46815523d0,  5.98089120d0,  5.81686657d0,  
     &5.53321815d0,  5.25477007d0, 11.02204549d0,  0.15679528d0,  
     &9.35167836d0,  9.06926079d0,  8.97241155d0,  8.90092807d0,  
     &8.85984840d0,  8.81736827d0,  8.79317710d0,  7.89969626d0,  
     &8.80588454d0,  8.42439218d0,  8.54289262d0,  8.47583370d0,  
     &8.45090888d0,  8.47339339d0,  7.83525634d0,  8.20702843d0,  
     &7.70559063d0,  7.32755997d0,  7.03887381d0,  6.68978720d0,  
     &6.05450052d0,  5.88752022d0,  5.70661499d0,  5.78450695d0,  
     &7.79780729d0,  7.26443867d0,  6.78151984d0,  6.67883169d0,  
     &6.39024318d0,  6.09527958d0, 11.79156076d0, 11.10997644d0,  
     &9.51377795d0,  8.67197068d0,  8.77140725d0,  8.65402716d0,  
     &8.53923501d0,  8.85024712d0 /)

! these new data are scaled with k2=4./3. and converted to a_0 via
! autoang=0.52917726d0
      rcov(1:94)=(/ 
     & 0.80628308d0, 1.15903197d0, 3.02356173d0, 2.36845659d0, 
     & 1.94011865d0, 1.88972601d0, 1.78894056d0, 1.58736983d0, 
     & 1.61256616d0, 1.68815527d0, 3.52748848d0, 3.14954334d0, 
     & 2.84718717d0, 2.62041997d0, 2.77159820d0, 2.57002732d0, 
     & 2.49443835d0, 2.41884923d0, 4.43455700d0, 3.88023730d0, 
     & 3.35111422d0, 3.07395437d0, 3.04875805d0, 2.77159820d0, 
     & 2.69600923d0, 2.62041997d0, 2.51963467d0, 2.49443835d0, 
     & 2.54483100d0, 2.74640188d0, 2.82199085d0, 2.74640188d0, 
     & 2.89757982d0, 2.77159820d0, 2.87238349d0, 2.94797246d0, 
     & 4.76210950d0, 4.20778980d0, 3.70386304d0, 3.50229216d0, 
     & 3.32591790d0, 3.12434702d0, 2.89757982d0, 2.84718717d0, 
     & 2.84718717d0, 2.72120556d0, 2.89757982d0, 3.09915070d0, 
     & 3.22513231d0, 3.17473967d0, 3.17473967d0, 3.09915070d0, 
     & 3.32591790d0, 3.30072128d0, 5.26603625d0, 4.43455700d0, 
     & 4.08180818d0, 3.70386304d0, 3.98102289d0, 3.95582657d0, 
     & 3.93062995d0, 3.90543362d0, 3.80464833d0, 3.82984466d0, 
     & 3.80464833d0, 3.77945201d0, 3.75425569d0, 3.75425569d0, 
     & 3.72905937d0, 3.85504098d0, 3.67866672d0, 3.45189952d0, 
     & 3.30072128d0, 3.09915070d0, 2.97316878d0, 2.92277614d0, 
     & 2.79679452d0, 2.82199085d0, 2.84718717d0, 3.32591790d0, 
     & 3.27552496d0, 3.27552496d0, 3.42670319d0, 3.30072128d0, 
     & 3.47709584d0, 3.57788113d0, 5.06446567d0, 4.56053862d0, 
     & 4.20778980d0, 3.98102289d0, 3.82984466d0, 3.85504098d0, 
     & 3.88023730d0, 3.90543362d0 /)

! DFT-D3
      step=0
      do iat=1,n-1
         do jat=iat+1,n
         step=step+1
         r=dist(step)
         c6=c6a(step)
! r2r4 stored in main as sqrt
         c8 =3.0d0*c6*r2r4(iz(iat))*r2r4(iz(jat))

! energy for BJ damping
          tmp=sqrt(c8/c6)
          e6(step)= c6/(r**6+(a1*tmp+a2)**6)
          e8(step)= c8/(r**8+(a1*tmp+a2)**8)
! calculate gradients
         if (igrad .eq. 1) then
! grad for BJ damping
          e6dr(step)=c6*(-6*r**5)/(r**6+(a1*tmp+a2)**6)**2
          e8dr(step)=c8*(-8*r**7)/(r**8+(a1*tmp+a2)**8)**2
         endif
         enddo
      enddo

      end subroutine edisp
      subroutine bmx2b109(bm1,bm2,nobm,nob)
!**********************************************************************
!  The subroutine eliminate the 2-body terms in Bowman's approach
!  as well as the the nonconnected terms in Bowman's approach
!
!  nobm : maximum number of basis
!  nob  : number of basis
!
!**********************************************************************

      implicit double precision (a-h,o-z)

      integer i
      double precision bm1(nobm), bm2(nob)

      bm2(1)=bm1(4)

      do i=2,4
        bm2(i)=bm1(i+4)
      enddo

      bm2(5)=bm1(10)

      do i=6,11
        bm2(i)=bm1(i+6)
      enddo

      bm2(12)=bm1(19)
      bm2(13)=bm1(21)

      do i=14,26
        bm2(i)=bm1(i+9)
      enddo

      bm2(27)=bm1(37)
      bm2(28)=bm1(39)

      do i=29,53
        bm2(i)=bm1(i+12)
      enddo

      bm2(54)=bm1(67)
      bm2(55)=bm1(69)
      bm2(56)=bm1(71)

      do i=57,97
        bm2(i)=bm1(i+16)
      enddo

      bm2(98)=bm1(115)
      bm2(99)=bm1(117)
      bm2(100)=bm1(119)

      do i=101,166
        bm2(i)=bm1(i+20)
      enddo

      bm2(167)=bm1(188)
      bm2(168)=bm1(190)
      bm2(169)=bm1(192)
      bm2(170)=bm1(194)

      do i=171,272
        bm2(i)=bm1(i+25)
      enddo

      bm2(273)=bm1(299)
      bm2(274)=bm1(301)
      bm2(275)=bm1(303)
      bm2(276)=bm1(305)

!     nob=276

      return

      end
