!! 2020/4/21
!! NN input: PIP only considering the permutation between molecules.
!! 7-45-45-1, fitting to 16000 energies
!! SCF/CVTZ-F12 + Correlation[AE-CCSD(T)-F12a]/CBS[CVDZ-F12,CVTZ-F12]
!! BSSE has not been considered
!! 2020/6/26
!! updated: fitting using exp(-r*ita), and analytical derivatives
!! input: r O1-C1-O2-C2 Bohr
!! output: v (interaction energy only) cm-1


      subroutine pes(x,igrad,path,p,g,d)

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
      double precision :: r_bond(6), r2(6), dvdr(6), drdx(6,12)
      double precision :: tx(3,natoms), txx(12), dvdx(12)
      integer :: iatom, idir, i, j, k, istate
      !initialize
      v=0.d0 
      p=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
        do idir=1,3
          tx(idir,iatom)=x(iatom, idir)/0.529177211
        enddo
      enddo
      k=0
      do i=1,3
      do j=i+1,4
         k=k+1
         r_bond(k)=dsqrt(dot_product((tx(:,i)-tx(:,j)),(tx(:,i)-tx(:,j))))
      enddo
      enddo

      call pes_ococ(r_bond, v, igrad, dvdr, path)

      dvdr=dvdr/8065.7112013/0.529177211

      ! convert dvdr to dvdx 
      do iatom=1,natoms
        do idir=1,3
          j=3*(iatom-1)+idir
          txx(j)=x(iatom, idir)
        enddo
      enddo
      r2=r_bond*0.529177211
      call EvdRdX(txx,r2,drdx)
      dvdx=0.d0
      do i=1,12
        do j=1,6
          dVdX(i)=dVdX(i) + dVdR(j)*dRdX(j,i)
        enddo
      enddo

      if (igrad==0) then 
        do istate=1,nstates
          p(istate)=v/8065.7112013
        enddo
      else if (igrad==1) then 
        do istate=1,nstates
          p(istate)=v/8065.7112013
        enddo
        do iatom=1,natoms
          do idir=1,3
            j=3*(iatom-1)+idir
            g(1,iatom,idir)=dVdX(j)
          enddo
        enddo
      endif 

      d=0.d0

      endsubroutine

  subroutine EvdRdX(X, R, dRdX)
      integer i,j
      double precision :: X(12), R(6), dRdX(6,12)
      do i=1,6
        do j=1,12
          dRdX(i,j)=0.0d0
        enddo
      enddo
! Start to calculate the non-zero dRdX
! dr1dx
      dRdX(1,1)=(x(1)-x(4))/r(1)
      dRdX(1,2)=(x(2)-x(5))/r(1)
      dRdX(1,3)=(x(3)-x(6))/r(1)
      dRdX(1,4)=-dRdX(1,1)
      dRdX(1,5)=-dRdX(1,2)
      dRdX(1,6)=-dRdX(1,3)

! dr2dx
      dRdX(2,1)=(x(1)-x(7))/r(2)
      dRdX(2,2)=(x(2)-x(8))/r(2)
      dRdX(2,3)=(x(3)-x(9))/r(2)
      dRdX(2,7)=-dRdX(2,1)
      dRdX(2,8)=-dRdX(2,2)
      dRdX(2,9)=-dRdX(2,3)

! dr3dx
      dRdX(3,1)=(x(1)-x(10))/r(3)
      dRdX(3,2)=(x(2)-x(11))/r(3)
      dRdX(3,3)=(x(3)-x(12))/r(3)
      dRdX(3,10)=-dRdX(3,1)
      dRdX(3,11)=-dRdX(3,2)
      dRdX(3,12)=-dRdX(3,3)

! dr4dx
      dRdX(4,4)=(x(4)-x(7))/r(4)
      dRdX(4,5)=(x(5)-x(8))/r(4)
      dRdX(4,6)=(x(6)-x(9))/r(4)
      dRdX(4,7)=-dRdX(4,4)
      dRdX(4,8)=-dRdX(4,5)
      dRdX(4,9)=-dRdX(4,6)

! dr5dx
      dRdX(5,4)=(x(4)-x(10))/r(5)
      dRdX(5,5)=(x(5)-x(11))/r(5)
      dRdX(5,6)=(x(6)-x(12))/r(5)
      dRdX(5,10)=-dRdX(5,4)
      dRdX(5,11)=-dRdX(5,5)
      dRdX(5,12)=-dRdX(5,6)

! dr6dx
      dRdX(6,7)=(x(7)-x(10))/r(6)
      dRdX(6,8)=(x(8)-x(11))/r(6)
      dRdX(6,9)=(x(9)-x(12))/r(6)
      dRdX(6,10)=-dRdX(6,7)
      dRdX(6,11)=-dRdX(6,8)
      dRdX(6,12)=-dRdX(6,9)
! Finish the calculation of non-zero dRdX

      return
      end subroutine EvdRdX

  subroutine pes_ococ(r,v,idv,dvdr,path)
  implicit none
  character(len=1024), intent(in) :: path
  real*8,intent(in) :: r(6)
  integer,intent(in) :: idv
  real*8,intent(out) :: v,dvdr(6)
  real*8,parameter :: pi=dacos(-1.d0)
  integer,parameter :: np=7
  real*8 :: p(np,1),vp(1),dvdp(np,1),dpdr(np,6)

  if(idv.ne.0 .and. idv.ne.1) stop "Only energy and gradient are available"

  call pip_ococ(r,p(:,1),idv,dpdr)
  call nnfit_ococ(np,1,p,vp,idv,dvdp,path)
  if(idv.eq.1) call dgemv('t',7,6,1.d0,dpdr,7,dvdp,1,0.d0,dvdr,1)
  v=vp(1)
  return

  contains

  subroutine pip_ococ(r,p,idv,dpdr)
  ! molecular based fundamental invariants
  implicit none
  real*8,intent(in) :: r(6)   ! atom order: O1 C1 O2 C2
  integer,intent(in) :: idv
  real*8,intent(out) :: p(7),dpdr(7,6)
  integer :: i

  real*8 :: g(6)
  real*8,parameter :: ita=0.3d0

  do i=1,6
     g(i)=exp(-r(i)*ita)
  enddo

  p(1)=g(1)+g(6)
  p(2)=g(3)+g(4)
  p(3)=g(1)**2+g(6)**2
  p(4)=g(3)**2+g(4)**2
  p(5)=g(1)*g(3)+g(4)*g(6)
  p(6)=g(2)
  p(7)=g(5)

  p(3)=dsqrt(p(3))
  p(4)=dsqrt(p(4))
  p(5)=dsqrt(p(5))

  if(idv.eq.0) return

  dpdr=0.d0

  dpdr(1,1)=1.d0
  dpdr(1,6)=1.d0
  dpdr(2,3)=1.d0
  dpdr(2,4)=1.d0
  dpdr(3,1)=0.5d0/p(3)*2.d0*g(1)
  dpdr(3,6)=0.5d0/p(3)*2.d0*g(6)
  dpdr(4,3)=0.5d0/p(4)*2.d0*g(3)
  dpdr(4,4)=0.5d0/p(4)*2.d0*g(4)
  dpdr(5,1)=0.5d0/p(5)*g(3)
  dpdr(5,3)=0.5d0/p(5)*g(1)
  dpdr(5,4)=0.5d0/p(5)*g(6)
  dpdr(5,6)=0.5d0/p(5)*g(4)
  dpdr(6,2)=1.d0
  dpdr(7,5)=1.d0

  do i=1,6
     dpdr(:,i)=dpdr(:,i)*g(i)*(-1.d0*ita)
  enddo

  return
  end subroutine pip_ococ

  subroutine nnfit_ococ(ndim,ntot,r0,vx,idv,dv,path)
  implicit none
  character(len=1024), intent(in) :: path
  integer,intent(in) :: ndim,ntot
  real*8,intent(in) :: r0(ndim,ntot)
  real*8,intent(out) :: vx(ntot)

  integer,parameter :: s0=7,s1=45,s2=45
  real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3,rg(2,s0),vg(2)

  integer,intent(in) :: idv
  real*8,intent(out) :: dv(ndim,ntot)
  integer :: i
  integer,save :: init=0
  integer,save :: fid
  character(len=1024) :: file_path1

  file_path1 = trim(path)//"/COCO/nn_ococ_w20.txt"

  if(ndim.ne.s0) stop "ndim .ne. s0"
  if (init.eq.0) then
    init=1
    fid=870711
    open(fid,file=file_path1,status='old',action='read')
    read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
    close(fid)
  endif
  if(ntot.ge.24) then
    call nsimx(ntot,r0,vx,idv,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
  else
    do i=1,ntot
    call nsim(r0(1,i),vx(i),idv,dv(:,i),s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
    enddo
  endif
  return
  end subroutine nnfit_ococ

  subroutine nsim(r0,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
  ! ---- simulate a neural network with two hidden layers
  ! ----      n0-n1-n2-1
  ! ---- blas routines used in this subroutine: dgemv, ddot
  implicit none
  integer,intent(in) :: n0,n1,n2,idv
  real*8,intent(in)  :: r0(n0),rg(2,n0),vg(2)
  real*8,intent(in)  :: w1(n0,n1),b1(n1)
  real*8,intent(in)  :: w2(n1,n2),b2(n2)
  real*8,intent(in)  :: w3(n2),b3
  real*8,intent(out) :: v,dv(n0)
  integer :: i,j,k
  real*8  :: r(n0),rgg(n0),vgg,ax(n1),bx(n2),maxx(n1),w2maxx(n1,n2),w3mbxx(n2)
  real*8  :: dvtm(n0,n2),rt1(n1),rt2(n2)
  real*8,external :: ddot
  v=0.d0
  r=r0
  vgg=vg(2)-vg(1)
  ! mapminmax [-1,1]
  do i=1,n0
    rgg(i)=rg(2,i)-rg(1,i)
    r(i)=2.d0*(r(i)-rg(1,i))/rgg(i)-1.d0
  end do
  ! 1st layer
  rt1=b1
  call dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1)
  ax=dtanh(rt1)
  ! 2nd layer
  rt2=b2
  call dgemv('t',n1,n2,1.d0,w2,n1,ax,1,1.d0,rt2,1)
  bx=dtanh(rt2)
  ! output layer
  v=b3+ddot(n2,w3,1,bx,1)
  !reverse map
  v=vgg*(v+1.d0)/2+vg(1)
  if(idv.ne.1) return
  ! calculate first derivatives, dv(i)=dv/dr(i)
  ! updated: 2020/6/26, speed up x5
  do j=1,n1
     maxx(j)=1-ax(j)**2
  enddo
  do k=1,n2
  do j=1,n1
     w2maxx(j,k)=maxx(j)*w2(j,k)
  enddo
  enddo
  do k=1,n2
     w3mbxx(k)=w3(k)*(1-bx(k)**2)
  enddo
  dv=0.d0
  dvtm=0.d0
  call dgemm('n','n',n0,n2,n1,1.d0,w1,n0,w2maxx,n1,0.d0,dvtm,n0)
  call dgemv('n',n0,n2,1.d0,dvtm,n0,w3mbxx,1,0.d0,dv,1)
  do i=1,n0
     dv(i)=dv(i)*vgg/rgg(i)
  enddo
  return
  end subroutine nsim

  subroutine nsimx(nt,r,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
  ! ---- simulate a neural network with two hidden layers
  ! ----    n0-n1(tansig)-n2(tansig)-1(pureline)
  ! ---- blas routines used in this subroutine: dgemm, dgemv
  ! ---- the syntax keeps the same with nsim, except input nt at first
  ! ---- Chen Jun, 2013/11/26
  !
  ! ---- description of parameters
  ! ---- nt  : total number of points in r,v and dv
  ! ---- r   : input of nn, r(n0,nt), r(1:n0,i) is the i-th input
  ! ---- v   : output of nn, v(nt), v(i) is the i-th input
  ! ---- idv : if idv==0, calc. v only, else calc. v and derivatives dv
  ! ---- dv  : derivatives of dv/dr, dv(n0,nt), dv(j,i)=d v(i) / d r(j,i)
  ! ---- n0  : neurons of the input layer, dimension of r
  ! ---- n1  : neurons of the first hidden layer
  ! ---- n2  : neurons of the second hidden layer
  ! ---- rg  : input ranges of the training set, rg(2,n0)
  ! ---- w1  : weights connect input layer with 1-st hidden layer
  ! ---- b1  : biases of 1-st hidden layer
  ! ---- w2  : weights connect 1-st hidden layer with 2-nd hidden layer
  ! ---- b2  : biases of 2-nd hidden layer
  ! ---- w3  : weights connect 2-nd hidden layer with output layer
  ! ---- b3  : biase of output layer
  ! ---- vg  : output range of the training set, vg(2)
  implicit none
  integer,intent(in) :: nt,idv,n0,n1,n2
  real*8,intent(in) :: r(n0,nt)
  real*8,intent(out) :: v(nt),dv(n0,nt)
  real*8,intent(in) :: w1(n0,n1),b1(n1),w2(n1,n2),b2(n2),w3(n2),b3
  real*8,intent(in) :: rg(2,n0),vg(2)
  real*8 :: x0(n0,nt),x1(nt,n1),x2(nt,n2),tmp,rgg(n0),vgg
  integer :: i,j,k,n

  v=0.d0
  do i=1,n0
    rgg(i)=rg(2,i)-rg(1,i)
  enddo
  vgg=vg(2)-vg(1)

  x0=r

  do i=1,n0
    x0(i,1:nt)=(x0(i,1:nt)-rg(1,i))/rgg(i)*2.d0-1.d0
  enddo

  do i=1,n1
   x1(1:nt,i)=b1(i)
  enddo

  ! ---- x1=trans(x0(n0,nt))*w1(n0,n1)+x1(nt,n1)
  !call gpu_dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
  call dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)

  x1=dtanh(x1)

  do i=1,n2
   x2(1:nt,i)=b2(i)
  enddo

  ! ---- x2=x1(nt,n1)*w2(n1,n2)+x2(nt,n2)
  !call gpu_dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
  call dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)

  x2=dtanh(x2)

  v(1:nt)=b3

  ! ---- v=x2(nt,n2)*w3(n2)+v(nt)
  call dgemv('n',nt,n2,1.d0,x2,nt,w3,1,1.d0,v,1)

  v=(v+1.d0)/2.d0*vgg+vg(1)

  if (idv.ne.1) return
  ! ---- calculate first derivatives
  dv=0.d0
  do n=1,nt
    do i=1,n0
      do k=1,n2
        tmp=0.d0
        do j=1,n1
          tmp=tmp+w2(j,k)*w1(i,j)*(1-x1(n,j)**2)
        enddo
        dv(i,n)=dv(i,n)+w3(k)*tmp*(1-x2(n,k)**2)
      enddo
      dv(i,n)=dv(i,n)*vgg/(rgg(i))
    enddo
  enddo
  return
  end subroutine nsimx

  end subroutine pes_ococ

  subroutine potco(r,v,dv)
  ! 2020/4/21  C-O AE-CCSD(T)(short range) MRCI(long range) / ACVTZ Energy in cm-1
  ! updated: 2020/6/25, dv/dr
  implicit none
  real*8,intent(in) :: r !! bohr
  real*8,intent(out) :: v,dv !! cm-1, cm-1/bohr
  real*8 :: x,y1j
  integer j
  integer,parameter :: ns=8
  real*8,parameter :: w1(ns)=[  2.5731334302226245d+00,&
                               -4.1635475296365723d+00,&
                                6.2652688635816842d+00,&
                               -8.1228824615517947d+00,&
                                2.5634824563114837d+01,&
                               -1.6643878666848999d+00,&
                                1.2863481593265147d+01,&
                               -5.1395051186435685d+00]
  real*8,parameter :: b1(ns)=[ -1.6421783363943476d+00,&
                                1.5364774081764951d+00,&
                               -1.1338455332512607d+00,&
                               -1.4448051810696709d+00,&
                                7.5217573991947644d+00,&
                               -1.4005229119446290d+00,&
                                1.1053854378210930d+01,&
                               -5.9299626180269485d+00]
  real*8,parameter :: w2(ns)=[ -8.8664820626030757d-03,&
                                7.8571245473773067d-03,&
                               -3.6411047563733342d-03,&
                               -4.0358215533209145d-03,&
                                9.6640587889626451d-04,&
                               -1.4325782866595651d+00,&
                                1.2002568907875554d-02,&
                                8.3983298757280007d+00]
  real*8,parameter :: b2=6.8970075393140338d+00
  real*8,parameter :: ra=1.4000000000000000d+00,rb=7.0000000000000000d+00
  real*8,parameter :: va=9.8087308941410573d-02,vb=1.9558422718340193d+05

  v=0.d0
  dv=0.d0
  x=2.d0*(r-ra)/(rb-ra)-1.d0
  do j=1,ns
    y1j=tanh(b1(j)+w1(j)*x)
    v=v+w2(j)*y1j
    dv=dv+w2(j)*(1.d0-y1j**2)*w1(j)
  enddo
  v=v+b2
  v=(v+1.d0)*(vb-va)/2.d0+va
  v=v+0.560096850315234d0
  dv=dv*(vb-va)/(rb-ra)
  return
  end subroutine potco
