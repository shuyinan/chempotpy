  module pes_shell
  implicit none

  real,parameter::auang=0.5291772083d0

  contains

  subroutine prepot(path)
    
    character(len=1024), intent(in) :: path
    character(len=1024) :: file_path1
    integer i1,m,mr
    double precision dc0(0:5,0:5),dw0(0:5,0:5)
    double precision coef(0:3337)
    
    common/NCOE/m,mr
    common/ch3ohcoe/dc0,dw0,coef
    
    ! for complete 6th-order fit
    m=3250 ; mr=22

    file_path1 = trim(path)//"/CH3OH/coeff.dat"
    open(20,file=file_path1,status='old')
    
    read(20,*)
    read(20,*)
    read(20,*)(coef(i1),i1=0,m+4*mr-1)
    close(20)
    
    return
  end subroutine
    
  !***********************************************************

  subroutine calcpot(v,xyz)

    double precision, dimension(3,6) :: xyz
    double precision, dimension(0:2,0:5) :: gf0
    double precision v
    
    call potshell(xyz,v)
    
    return
  end subroutine calcpot
    
  !***********************************************************

  subroutine potshell(xn,V)
    
    integer m,mr
    double precision dc0(0:5,0:5),dw0(0:5,0:5)
    double precision coef(0:3337), f0, gf0(0:2,0:5)
    
    common/NCOE/m,mr
    common/ch3ohcoe/dc0,dw0,coef
    
    double precision V, xn(0:2,0:5)
    double precision d0(0:5,0:5), dw(0:5,0:5)
    double precision vec(0:m+4*mr-1)
    
    call getvec(m,mr,xn,vec)
    f0=dot_product(coef,vec)
    
    V=f0
    
    return
  end subroutine potshell
    
  !****************************************************************
     
  subroutine getd0 (nk, r0, d0)
    integer, parameter :: wp=selected_real_kind(12,300)
    integer nk
    double precision r0(0:nk-1,0:nk-1), d0(0:nk-1,0:nk-1)
    integer i, j
    do i = 0, nk-1
       d0(i,i) = 0
       do j = i+1, nk-1
          d0(i,j) = dexp(-r0(i,j)/3)
          d0(j,i) = d0(i,j)
       enddo
    enddo

    return
  end subroutine getd0

  !****************************************************************

  subroutine getr0 (nk, xn, r0)
    integer, parameter :: wp=selected_real_kind(12,300)
    integer nk
    double precision xn(0:2,0:nk-1), r0(0:nk-1,0:nk-1)
    integer i, j
    do i = 0, nk-1
       r0(i,i) = 0
       do j = i+1, nk-1
          r0(i,j) = sqrt((xn(0,j)-xn(0,i))**2+(xn(1,j)-xn(1,i))**2+ &
                    (xn(2,j)-xn(2,i))**2)
          r0(j,i) = r0(i,j)
       enddo
    enddo

    return
  end subroutine getr0
  
  !******************************************************************
    
  subroutine getrvec (ms, r, vec)
    integer, parameter :: wp=selected_real_kind(12,300)
    ! version for X4YZ
    integer nk, ms
    parameter (nk=6)
    double precision r(0:nk-1,0:nk-1), vec(0:ms-1)
    integer i, j
    double precision x(0:3), r1(0:nk-1,0:nk-1), t0, t1
    !-----------------------------------------------------------------------
    ! Test for compatibility
    if (.not.(ms.eq.1.or.ms.eq.5)) then
       stop 'getrvec - wrong dimension'
    endif
    ! Computation
    x = 0
    do i = 0, nk-1
       do j = 0, nk-1
          if (i.eq.j) then
             r1(i,j) = 0
          else
             r1(i,j) = dexp(-r(i,j))/r(i,j)
          endif
       enddo
    enddo
    ! XX distance
    t0 = 0
    do i = 0, nk-3
       do j = i+1, nk-3
          t0 = t0+r1(i,j)
       enddo
    enddo
    x(0) = t0/((nk-2)*(nk-3)/2)
    ! XY and XZ distances
    t0 = 0 ; t1 = 0
    do i = 0, nk-3
       t0 = t0+r1(i,nk-2)
       t1 = t1+r1(i,nk-1)
    enddo
    x(1) = t0/(nk-2)
    x(2) = t1/(nk-2)
    ! YZ distance
    x(3) = r1(nk-2,nk-1)
    ! set vec
    vec(0) = 1
    if (5.le.ms) then
       vec(1:4) = x
    endif
    return
  end subroutine getrvec
    
  !*****************************************************************
    
  subroutine getvec (ms, mr, xn, vec)
    integer, parameter :: wp=selected_real_kind(12,300)
    ! version for H4CO (reordering of CH3OH).
    integer nk, ms, mr
    parameter (nk=6)
    double precision xn(0:2,0:nk-1), vec(0:ms+4*mr-1)
    integer k, l
    double precision rvec(0:4), d0(0:nk-1,0:nk-1), r0(0:nk-1,0:nk-1)
    !-----------------------------------------------------------------------
    vec = 0
    call getr0 (nk, xn, r0)
    call getd0 (nk, r0, d0)
    call getv_x4yz (ms, d0, vec(0:ms-1))
    call getrvec (5, r0, rvec)
    do l = 0, mr-1
       do k = 0, 3
          vec(ms+4*l+k) = rvec(k+1)*vec(l)
       enddo
    enddo

    return
  end subroutine getvec
    
  !*************************************************************************
    
  subroutine getv_x4yz (m, d, vec)
    integer, parameter :: wp=selected_real_kind(12,300)
    integer, parameter :: nk=6
    integer m
    double precision d(0:nk-1,0:nk-1), vec(0:m-1)
    ! version for molecule X4YZ
    ! MolienSeries(0:8): 1 4 17 65 230 736 2197 6093 15864
    ! #Primaries(1:8):   4 4 4 3 0 0 0 0
    ! #Secondaries(1:8): 0 3 13 32 62 129 221 335
    integer, parameter :: l0=1, l1=l0+4, l2=l1+17, l3=l2+65,l4=l3+230, &
     l5=l4+736, l6=l5+2197, l7=l6+6093-221
    !! We haven't incorporated the secondaries at degree 7
    integer, parameter :: np1=4, np2=4, np3=4, np4=3, np5=0, np6=0, &
                          np7=0
    integer, parameter :: ns1=0, ns2=3, ns3=13, ns4=32, ns5=62, &
                          ns6=129, ns7=221
    double precision x(0:np1-1), y(0:np2-1), z(0:np3-1), u(0:np4-1)
    double precision x2(0:np1-1), x3(0:np1-1), x4(0:np1-1), x5(0:np1-1), &
                     x6(0:np1-1), x7(0:np1-1), y2(0:np2-1), y3(0:np2-1), &
                     z2(0:np3-1)
    double precision ys(0:ns2-1), zs(0:ns3-1), us(0:ns4-1), &
                     vs(0:ns5-1), ws(0:ns6-1), w7s(0:ns7-1)
    integer mdeg, i, j, k, l, i0, i1, i2, i3, j0, k0
    double precision d2(0:nk-1,0:nk-1), d3(0:nk-1,0:nk-1), d4(0:nk-1,0:nk-1), &
                     d5(0:nk-1,0:nk-1), d6(0:nk-1,0:nk-1), d7(0:nk-1,0:nk-1)
    double precision t0
    double precision pol2, pol3, pol4, pol5, pol6, pol7
    pol2(t0) = t0**2
    pol3(t0) = t0**3
    pol4(t0) = t0**4
    pol5(t0) = t0**5
    pol6(t0) = t0**6
    pol7(t0) = t0**7
    !-----------------------------------------------------------------------
    ! Test for compatibility, set mdeg
    select case (m)
      case (l0)
        mdeg = 0
      case (l1)
        mdeg = 1
      case (l2)
        mdeg = 2
      case (l3)
        mdeg = 3
      case (l4)
        mdeg = 4
      case (l5)
        mdeg = 5
      case (l6)
        mdeg = 6
      case (l7)
        mdeg = 7
      case default
        stop 'getv - wrong dimension'
    endselect
    ! auxiliary distances
    do i = 0, nk-1
       do j = i+1, nk-1
          d2(i,j) = pol2(d(i,j))
          d2(j,i) = d2(i,j)
          d3(i,j) = pol3(d(i,j))
          d3(j,i) = d3(i,j)
          d4(i,j) = pol4(d(i,j))
          d4(j,i) = d4(i,j)
          d5(i,j) = pol5(d(i,j))
          d5(j,i) = d5(i,j)
          d6(i,j) = pol6(d(i,j))
          d6(j,i) = d6(i,j)
          d7(i,j) = pol7(d(i,j))
          d7(j,i) = d7(i,j)
       enddo
    enddo
    ! Primary Invariants
    !      x = 0 ; y = 0 ; z = 0 ; u = 0 ; v = 0 ; w = 0 ; w7 = 0
    x = 0 ; y = 0 ; z = 0 ; u = 0 
    
    do i0 = 0, 3
       t0 = 0
       do i1 = 0, 3
          if (i1.ne.i0) then
             t0 = t0+d(i0,i1)/3
             x(0) = x(0)+d(i0,i1)/12
             y(0) = y(0)+d2(i0,i1)/12
             z(0) = z(0)+d3(i0,i1)/12
          endif
       enddo
       y(1) = y(1)+pol2(t0)/4
       z(1) = z(1)+pol3(t0)/4
       u(0) = u(0)+pol4(t0)/4
    enddo

    x(1) = sum(d(0:3,4))/4
    y(2) = sum(d2(0:3,4))/4
    z(2) = sum(d3(0:3,4))/4
    u(1) = sum(d4(0:3,4))/4
    x(2) = sum(d(0:3,5))/4
    y(3) = sum(d2(0:3,5))/4
    z(3) = sum(d3(0:3,5))/4
    u(2) = sum(d4(0:3,5))/4
    x(3) = d(4,5)
    ! Required powers
    do i = 0, np1-1
       x2(i) = pol2(x(i))
       x3(i) = pol3(x(i))
       x4(i) = pol4(x(i))
       x5(i) = pol5(x(i))
       x6(i) = pol6(x(i))
       x7(i) = pol7(x(i))
    enddo
    do i = 0, np2-1
       y2(i) = pol2(y(i))
       y3(i) = pol3(y(i))
    enddo
    do i = 0, np3-1
       z2(i) = pol2(z(i))
    enddo
    ! Secondary Invariants
    !! ys(0:2), zs(0:12), us(0:31), vs(0:61), ws(0:128), w7s(0:220)
    !! reducible: us(0:5), vs(0:38), ws(0:127), w7s(0:..)
    ys = 0 ; zs = 0 ; us = 0 ; vs = 0 ; ws = 0 ; w7s = 0
    ! Irreducible secondaries
    do i0 = 0, 3
       do i1 = 0, 3
          if (i1.ne.i0) then
             us(6) = us(6)+d4(i0,i1)/12
             vs(39) = vs(39)+d5(i0,i1)/12
             do i2 = 0, 3
                if (i2.ne.i1.and.i2.ne.i0) then
                   zs(0) = zs(0)+d2(i0,i1)*d(i0,i2)/24
                endif
             enddo
          endif
       enddo
    enddo
    j0 = 4 ; k0 = 5
    do i0 = 0, 3
       ys(2) = ys(2)+d(i0,j0)*d(i0,k0)/4
       zs(8) = zs(8)+d2(i0,j0)*d(i0,k0)/4
       zs(11) = zs(11)+d(i0,j0)*d2(i0,k0)/4
       us(20) = us(20)+d3(i0,j0)*d(i0,k0)/4
       us(27) = us(27)+d2(i0,j0)*d2(i0,k0)/4
       us(30) = us(30)+d(i0,j0)*d3(i0,k0)/4
       do i1 = 0, 3
          if (i1.ne.i0) then
             ys(0) = ys(0)+d(i0,i1)*d(i0,j0)/12
             ys(1) = ys(1)+d(i0,i1)*d(i0,k0)/12
             zs(1) = zs(1)+d2(i0,i1)*d(i0,j0)/12
             zs(3) = zs(3)+d(i0,i1)*d2(i0,j0)/12
             zs(4) = zs(4)+d(i0,i1)*d(i0,j0)*d(i1,j0)/12
             zs(5) = zs(5)+d2(i0,i1)*d(i0,k0)/12
             zs(7) = zs(7)+d(i0,i1)*d(i0,j0)*d(i0,k0)/12
             zs(9) = zs(9)+d(i0,i1)*d(i1,j0)*d(i0,k0)/12
             zs(10) = zs(10)+d(i0,i1)*d2(i0,k0)/12
             zs(12) = zs(12)+d(i0,i1)*d(i0,k0)*d(i1,k0)/12
             us(7) = us(7)+d3(i0,i1)*d(i0,j0)/12
             us(10) = us(10)+d2(i0,i1)*d2(i0,j0)/12
             us(12) = us(12)+d(i0,i1)*d3(i0,j0)/12
             us(13) = us(13)+d2(i0,i1)*d(i0,j0)*d(i1,j0)/12
             us(14) = us(14)+d3(i0,i1)*d(i0,k0)/12
             us(17) = us(17)+d2(i0,i1)*d(i0,j0)*d(i0,k0)/12
             us(19) = us(19)+d(i0,i1)*d2(i0,j0)*d(i0,k0)/12
             us(21) = us(21)+d2(i0,i1)*d(i1,j0)*d(i0,k0)/12
             us(23) = us(23)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/12
             us(24) = us(24)+d2(i0,i1)*d2(i0,k0)/12
             us(26) = us(26)+d(i0,i1)*d(i0,j0)*d2(i0,k0)/12
             us(28) = us(28)+d(i0,i1)*d(i1,j0)*d2(i0,k0)/12
             us(29) = us(29)+d(i0,i1)*d3(i0,k0)/12
             us(31) = us(31)+d2(i0,i1)*d(i0,k0)*d(i1,k0)/12
             vs(40) = vs(40)+d4(i0,i1)*d(i0,j0)/12
             vs(42) = vs(42)+d3(i0,i1)*d2(i0,j0)/12
             vs(44) = vs(44)+d2(i0,i1)*d3(i0,j0)/12
             vs(45) = vs(45)+d(i0,i1)*d3(i0,j0)*d(i1,j0)/12
             vs(46) = vs(46)+d4(i0,i1)*d(i0,k0)/12
             vs(48) = vs(48)+d3(i0,i1)*d(i0,j0)*d(i0,k0)/12
             vs(50) = vs(50)+d2(i0,i1)*d2(i0,j0)*d(i0,k0)/12
             vs(52) = vs(52)+d2(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/12
             vs(53) = vs(53)+d(i0,i1)*d2(i0,j0)*d(i1,j0)*d(i0,k0)/12
             vs(54) = vs(54)+d3(i0,i1)*d2(i0,k0)/12
             vs(56) = vs(56)+d2(i0,i1)*d(i0,j0)*d2(i0,k0)/12
             vs(57) = vs(57)+d2(i0,i1)*d(i1,j0)*d2(i0,k0)/12
             vs(58) = vs(58)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d2(i0,k0)/12
             vs(59) = vs(59)+d2(i0,i1)*d3(i0,k0)/12
             vs(60) = vs(60)+d(i0,i1)*d(i1,j0)*d3(i0,k0)/12
             vs(61) = vs(61)+d(i0,i1)*d3(i0,k0)*d(i1,k0)/12
             do i2 = 0, 3
                if (i2.ne.i1.and.i2.ne.i0) then
                   zs(2) = zs(2)+d(i0,i1)*d(i0,i2)*d(i0,j0)/24
                   zs(6) = zs(6)+d(i0,i1)*d(i0,i2)*d(i0,k0)/24
                   us(8) = us(8)+d2(i0,i1)*d(i0,i2)*d(i0,j0)/24
                   us(9) = us(9)+d2(i0,i1)*d(i1,i2)*d(i0,j0)/24
                   us(11) = us(11)+d(i0,i1)*d(i0,i2)*d2(i0,j0)/24
                   us(15) = us(15)+d2(i0,i1)*d(i0,i2)*d(i0,k0)/24
                   us(16) = us(16)+d2(i0,i1)*d(i1,i2)*d(i0,k0)/24
                   us(18) = us(18)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)/24
                   us(22) = us(22)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,k0)/24
                   us(25) = us(25)+d(i0,i1)*d(i0,i2)*d2(i0,k0)/24
                   vs(41) = vs(41)+d3(i0,i1)*d(i0,i2)*d(i0,j0)/24
                   vs(43) = vs(43)+d2(i0,i1)*d(i0,i2)*d2(i0,j0)/24
                   vs(47) = vs(47)+d3(i0,i1)*d(i0,i2)*d(i0,k0)/24
                   vs(49) = vs(49)+d2(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)/24
                   vs(51) = vs(51)+d2(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,k0)/24
                   vs(55) = vs(55)+d2(i0,i1)*d(i0,i2)*d2(i0,k0)/24
                   ws(128) = ws(128)+d2(i0,i1)*d2(i0,i2)*d(i1,j0)*d(i0,k0)/24
                endif
             enddo
          endif
       enddo
    enddo
    ! Reducible secondaries
    ! us(0:5)
    k = 0
    do j = 0, 2
       do i = 0, j-1
          us(k) = ys(i)*ys(j) ; k = k+1
       enddo
       us(k) = pol2(ys(j)) ; k = k+1
    enddo
    ! vs(0:38)
    do j = 0, 12
       do i = 0, 2
         vs(3*j+i) = ys(i)*zs(j)
       enddo
    enddo
    ! ws(0:127)
    ws(0) = pol2(zs(0))
    ws(1) = zs(0)*zs(1)
    ws(2) = zs(0)*zs(2)
    ws(3) = zs(0)*zs(3)
    ws(4) = zs(1)*zs(3)
    ws(5) = zs(2)*zs(3)
    ws(6) = zs(0)*zs(4)
    ws(7) = zs(1)*zs(4)
    ws(8) = zs(3)*zs(4)
    ws(9) = pol2(zs(4))
    ws(10) = zs(0)*zs(5)
    ws(11) = zs(0)*zs(6)
    ws(12) = zs(1)*zs(6)
    ws(13) = zs(0)*zs(7)
    ws(14) = zs(1)*zs(7)
    ws(15) = zs(0)*zs(8)
    ws(16) = zs(1)*zs(8)
    ws(17) = zs(2)*zs(8)
    ws(18) = zs(4)*zs(8)
    ws(19) = zs(0)*zs(9)
    ws(20) = zs(1)*zs(9)
    ws(21) = zs(2)*zs(9)
    ws(22) = zs(3)*zs(9)
    ws(23) = zs(4)*zs(9)
    ws(24) = zs(0)*zs(10)
    ws(25) = zs(1)*zs(10)
    ws(26) = zs(4)*zs(10)
    ws(27) = zs(5)*zs(10)
    ws(28) = zs(6)*zs(10)
    ws(29) = zs(0)*zs(11)
    ws(30) = zs(1)*zs(11)
    ws(31) = zs(2)*zs(11)
    ws(32) = zs(3)*zs(11)
    ws(33) = zs(4)*zs(11)
    ws(34) = zs(5)*zs(11)
    ws(35) = zs(6)*zs(11)
    ws(36) = zs(7)*zs(11)
    ws(37) = zs(8)*zs(11)
    ws(38) = zs(0)*zs(12)
    ws(39) = zs(1)*zs(12)
    ws(40) = zs(2)*zs(12)
    ws(41) = zs(3)*zs(12)
    ws(42) = zs(4)*zs(12)
    ws(43) = zs(5)*zs(12)
    ws(44) = zs(7)*zs(12)
    ws(45) = zs(8)*zs(12)
    ws(46) = zs(9)*zs(12)
    ws(47) = zs(10)*zs(12)
    ws(48) = zs(11)*zs(12)
    ws(49) = pol2(zs(12))
    do j = 0, 25
       do i = 0, 2
          ws(50+3*j+i) = ys(i)*us(6+j)
       enddo
    enddo
    ! w7s(0:..) !! to follow
    ! Compute vec(0:*).  The code was created using these parameters:
    ! MolienSeries(0:*): 1 4 17 65 230 736 2197 6093
    ! #Primaries(1:*):   4 4 4 3 0 0 0
    ! #Secondaries(1:*): 0 3 13 32 62 129 221
    ! The MolienSeries partial sums (allowed size of vec) are:
    ! Molien Sums(0:*):  1 5 22 87 317 1053 3250 9343 9347
    ! constant term
    vec(0) = 1
    ! first degree terms
    if (1.le.mdeg) then
       vec(1) = x(0)
       vec(2) = x(1)
       vec(3) = x(2)
       vec(4) = x(3)
    endif
    ! second degree terms
    if (2.le.mdeg) then
       vec(5) = x2(0)
       vec(6:8) = x(0)*vec(2:4)
       vec(9) = x2(1)
       vec(10:11) = x(1)*vec(3:4)
       vec(12) = x2(2)
       vec(13:13) = x(2)*vec(4:4)
       vec(14) = x2(3)
       vec(15:14) = x(3)*vec(5:4)
       vec(15) = y(0)
       vec(16) = y(1)
       vec(17) = y(2)
       vec(18) = y(3)
       vec(19:21) = ys(0:2)
    endif
    ! third degree terms
    if (3.le.mdeg) then
       vec(22) = x3(0)
       vec(23:25) = x2(0)*vec(2:4)
       vec(26:38) = x(0)*vec(9:21)
       vec(39) = x3(1)
       vec(40:41) = x2(1)*vec(3:4)
       vec(42:51) = x(1)*vec(12:21)
       vec(52) = x3(2)
       vec(53:53) = x2(2)*vec(4:4)
       vec(54:61) = x(2)*vec(14:21)
       vec(62) = x3(3)
       vec(63:62) = x2(3)*vec(5:4)
       vec(63:69) = x(3)*vec(15:21)
       vec(70) = z(0)
       vec(71) = z(1)
       vec(72) = z(2)
       vec(73) = z(3)
       vec(74:86) = zs(0:12)
    endif
    ! fourth degree terms
    if (4.le.mdeg) then
       vec(87) = x4(0)
       vec(88:90) = x3(0)*vec(2:4)
       vec(91:103) = x2(0)*vec(9:21)
       vec(104:151) = x(0)*vec(39:86)
       vec(152) = x4(1)
       vec(153:154) = x3(1)*vec(3:4)
       vec(155:164) = x2(1)*vec(12:21)
       vec(165:199) = x(1)*vec(52:86)
       vec(200) = x4(2)
       vec(201:201) = x3(2)*vec(4:4)
       vec(202:209) = x2(2)*vec(14:21)
       vec(210:234) = x(2)*vec(62:86)
       vec(235) = x4(3)
       vec(236:235) = x3(3)*vec(5:4)
       vec(236:242) = x2(3)*vec(15:21)
       vec(243:259) = x(3)*vec(70:86)
       vec(260) = y2(0)
       vec(261:266) = y(0)*vec(16:21)
       vec(267) = y2(1)
       vec(268:272) = y(1)*vec(17:21)
       vec(273) = y2(2)
       vec(274:277) = y(2)*vec(18:21)
       vec(278) = y2(3)
       vec(279:281) = y(3)*vec(19:21)
       vec(282) = u(0)
       vec(283) = u(1)
       vec(284) = u(2)
       vec(285:316) = us(0:31)
    endif
    ! fifth degree terms
    if (5.le.mdeg) then
       vec(317) = x5(0)
       vec(318:320) = x4(0)*vec(2:4)
       vec(321:333) = x3(0)*vec(9:21)
       vec(334:381) = x2(0)*vec(39:86)
       vec(382:546) = x(0)*vec(152:316)
       vec(547) = x5(1)
       vec(548:549) = x4(1)*vec(3:4)
       vec(550:559) = x3(1)*vec(12:21)
       vec(560:594) = x2(1)*vec(52:86)
       vec(595:711) = x(1)*vec(200:316)
       vec(712) = x5(2)
       vec(713:713) = x4(2)*vec(4:4)
       vec(714:721) = x3(2)*vec(14:21)
       vec(722:746) = x2(2)*vec(62:86)
       vec(747:828) = x(2)*vec(235:316)
       vec(829) = x5(3)
       vec(830:829) = x4(3)*vec(5:4)
       vec(830:836) = x3(3)*vec(15:21)
       vec(837:853) = x2(3)*vec(70:86)
       vec(854:910) = x(3)*vec(260:316)
       vec(911:927) = y(0)*vec(70:86)
       vec(928:944) = y(1)*vec(70:86)
       vec(945:961) = y(2)*vec(70:86)
       vec(962:978) = y(3)*vec(70:86)
       vec(979:981) = z(0)*ys(0:2)
       vec(982:984) = z(1)*ys(0:2)
       vec(985:987) = z(2)*ys(0:2)
       vec(988:990) = z(3)*ys(0:2)
       vec(991:1052) = vs(0:61)
    endif
    ! sixth degree terms
    if (6.le.mdeg) then
       vec(1053) = x6(0)
       vec(1054:1056) = x5(0)*vec(2:4)
       vec(1057:1069) = x4(0)*vec(9:21)
       vec(1070:1117) = x3(0)*vec(39:86)
       vec(1118:1282) = x2(0)*vec(152:316)
       vec(1283:1788) = x(0)*vec(547:1052)
       vec(1789) = x6(1)
       vec(1790:1791) = x5(1)*vec(3:4)
       vec(1792:1801) = x4(1)*vec(12:21)
       vec(1802:1836) = x3(1)*vec(52:86)
       vec(1837:1953) = x2(1)*vec(200:316)
       vec(1954:2294) = x(1)*vec(712:1052)
       vec(2295) = x6(2)
       vec(2296:2296) = x5(2)*vec(4:4)
       vec(2297:2304) = x4(2)*vec(14:21)
       vec(2305:2329) = x3(2)*vec(62:86)
       vec(2330:2411) = x2(2)*vec(235:316)
       vec(2412:2635) = x(2)*vec(829:1052)
       vec(2636) = x6(3)
       vec(2637:2636) = x5(3)*vec(5:4)
       vec(2637:2643) = x4(3)*vec(15:21)
       vec(2644:2660) = x3(3)*vec(70:86)
       vec(2661:2717) = x2(3)*vec(260:316)
       vec(2718:2859) = x(3)*vec(911:1052)
       vec(2860) = y3(0)
       vec(2861:2866) = y2(0)*vec(16:21)
       vec(2867:2916) = y(0)*vec(267:316)
       vec(2917) = y3(1)
       vec(2918:2922) = y2(1)*vec(17:21)
       vec(2923:2966) = y(1)*vec(273:316)
       vec(2967) = y3(2)
       vec(2968:2971) = y2(2)*vec(18:21)
       vec(2972:3010) = y(2)*vec(278:316)
       vec(3011) = y3(3)
       vec(3012:3014) = y2(3)*vec(19:21)
       vec(3015:3049) = y(3)*vec(282:316)
       vec(3050) = z2(0)
       vec(3051:3066) = z(0)*vec(71:86)
       vec(3067) = z2(1)
       vec(3068:3082) = z(1)*vec(72:86)
       vec(3083) = z2(2)
       vec(3084:3097) = z(2)*vec(73:86)
       vec(3098) = z2(3)
       vec(3099:3111) = z(3)*vec(74:86)
       vec(3112:3114) = u(0)*ys(0:2)
       vec(3115:3117) = u(1)*ys(0:2)
       vec(3118:3120) = u(2)*ys(0:2)
       vec(3121:3249) = ws(0:128)
    endif
    ! seventh degree terms
    if (7.le.mdeg) then
       vec(3250) = x7(0)
       vec(3251:3253) = x6(0)*vec(2:4)
       vec(3254:3266) = x5(0)*vec(9:21)
       vec(3267:3314) = x4(0)*vec(39:86)
       vec(3315:3479) = x3(0)*vec(152:316)
       vec(3480:3985) = x2(0)*vec(547:1052)
       vec(3986:5446) = x(0)*vec(1789:3249)
       vec(5447) = x7(1)
       vec(5448:5449) = x6(1)*vec(3:4)
       vec(5450:5459) = x5(1)*vec(12:21)
       vec(5460:5494) = x4(1)*vec(52:86)
       vec(5495:5611) = x3(1)*vec(200:316)
       vec(5612:5952) = x2(1)*vec(712:1052)
       vec(5953:6907) = x(1)*vec(2295:3249)
       vec(6908) = x7(2)
       vec(6909:6909) = x6(2)*vec(4:4)
       vec(6910:6917) = x5(2)*vec(14:21)
       vec(6918:6942) = x4(2)*vec(62:86)
       vec(6943:7024) = x3(2)*vec(235:316)
       vec(7025:7248) = x2(2)*vec(829:1052)
       vec(7249:7862) = x(2)*vec(2636:3249)
       vec(7863) = x7(3)
       vec(7864:7863) = x6(3)*vec(5:4)
       vec(7864:7870) = x5(3)*vec(15:21)
       vec(7871:7887) = x4(3)*vec(70:86)
       vec(7888:7944) = x3(3)*vec(260:316)
       vec(7945:8086) = x2(3)*vec(911:1052)
       vec(8087:8476) = x(3)*vec(2860:3249)
       vec(8477:8493) = y2(0)*vec(70:86)
       vec(8494:8618) = y(0)*vec(928:1052)
       vec(8619:8635) = y2(1)*vec(70:86)
       vec(8636:8743) = y(1)*vec(945:1052)
       vec(8744:8760) = y2(2)*vec(70:86)
       vec(8761:8851) = y(2)*vec(962:1052)
       vec(8852:8868) = y2(3)*vec(70:86)
       vec(8869:8942) = y(3)*vec(979:1052)
       vec(8943:8977) = z(0)*vec(282:316)
       vec(8978:9012) = z(1)*vec(282:316)
       vec(9013:9047) = z(2)*vec(282:316)
       vec(9048:9082) = z(3)*vec(282:316)
       vec(9083:9095) = u(0)*zs(0:12)
       vec(9096:9108) = u(1)*zs(0:12)
       vec(9109:9121) = u(2)*zs(0:12)
    !! vec(9122:9342) = w7s(0:220)
    !! Secondaries at degree 7 yet to follow
    endif
    return
  end subroutine getv_x4yz
    
  end module pes_shell


      subroutine pes(x,igrad,path,p,g,d)

      use pes_shell
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=6
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
          tx(idir, iatom)=x(iatom, idir)/0.529177211
        enddo
      enddo

      if (igrad==0) then
        call prepot(path)
        call calcpot(v,tx)
      else
        write (*,*) 'Only energy is available'
      endif

      !employ the minimum to be -115.57493896563480 hartree.
      v=(v+115.57493896563480)*27.211386
      do istate=1,nstates
        p(istate)=v
      enddo

      endsubroutine
