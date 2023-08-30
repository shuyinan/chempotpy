      subroutine pes(x,igrad,pin,grad,dvec)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=3
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: pin(nstates)
      double precision, intent(out) :: grad(nstates,natoms,3)
      double precision, intent(out) :: dvec(nstates,nstates,natoms,3)
 
      integer :: nt
      double precision :: r(1,3), r2(3)
      double precision :: u11(1), u12(1), u13(1), u22(1), u23(1), u33(1)
      double precision :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
      double precision :: de_u13(3,1), de_u23(3,1), de_u33(3,1)
      double precision :: dpem(nstates,nstates), dudr(nstates,nstates,3)
      double precision :: trans(nstates,nstates)
      double precision :: tmpmat(nstates,nstates)
      double precision :: hr(nstates,nstates,3), gr(nstates,3)
      double precision :: tx(9), drdx(3,9)
      double precision :: hx(nstates,nstates,9), gx(nstates,9)
      logical, save :: first_time_data=.true.
      integer :: iatom,idir, i, j, k, l

      !initialize 
      u=0.d0
      pin=0.d0
      grad=0.d0
      dvec=0.d0

      nt=1
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        tx(j)=x(iatom,idir)
      enddo
      enddo

      ! input cartesian is HHH
      r(1,1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      r(1,2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      r(1,3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211

      call pot(r,u11,de_u11,1,1)
      call pot(r,u12,de_u12,1,2)
      call pot(r,u22,de_u22,1,3)
      call pot(r,u13,de_u13,1,4)
      call pot(r,u23,de_u23,1,5)
      call pot(r,u33,de_u33,1,6)

      dpem(1,1)=u11(1)*27.211386
      dpem(1,2)=u12(1)*27.211386
      dpem(2,1)=u12(1)*27.211386
      dpem(2,2)=u22(1)*27.211386
      dpem(1,3)=u13(1)*27.211386
      dpem(3,1)=u13(1)*27.211386
      dpem(2,3)=u23(1)*27.211386
      dpem(3,2)=u23(1)*27.211386
      dpem(3,3)=u33(1)*27.211386

      dudr(1,1,:)=de_u11(:,1)*51.422067
      dudr(1,2,:)=de_u12(:,1)*51.422067
      dudr(2,1,:)=de_u12(:,1)*51.422067
      dudr(2,2,:)=de_u22(:,1)*51.422067
      dudr(1,3,:)=de_u13(:,1)*51.422067
      dudr(3,1,:)=de_u13(:,1)*51.422067
      dudr(2,3,:)=de_u23(:,1)*51.422067
      dudr(3,2,:)=de_u23(:,1)*51.422067
      dudr(3,3,:)=de_u33(:,1)*51.422067

      call diagonalize(nstates,dpem,pin,trans)

      do i=1,3
        tmpmat(:,:)=dudr(:,:,i)
        tmpmat=matmul(transpose(trans), matmul(tmpmat,trans))
        do j=1,nstates
          do k=j+1,nstates
            if ((pin(k)-pin(j)) > 1.d-8) then
              hr(j,k,i)=tmpmat(j,k)/(pin(k)-pin(j))
            else
              hr(j,k,i)=tmpmat(j,k)/1.d-8
            endif
            hr(k,j,i)=-hr(j,k,i)
          enddo 
        enddo
        do j=1,nstates
          gr(j,i)=tmpmat(j,j)
        enddo
      enddo
      
      r2=r(1,:)*0.529177211
      call evdrdx(tx, r2, drdx)
      hx=0.d0
      gx=0.d0
      do i=1,9
      do j=1,3
        hx(:,:,i)=hx(:,:,i)+hr(:,:,j)*drdx(j,i)
        gx(:,i)=gx(:,i)+gr(:,j)*drdx(j,i)
      enddo
      enddo
  
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        grad(:,iatom,idir)=gx(:,j)
        dvec(:,:,iatom,idir)=hx(:,:,j)
      enddo
      enddo


      endsubroutine

      subroutine EvdRdX(X,r,drdx)

      integer i,j
      double precision, intent(in) :: X(9), R(3)
      double precision, intent(out) :: dRdX(3,9)

! Initialize dRdX(3,9)
      do i=1,3
        do j=1,9
          dRdX(i,j)=0.0d0
        enddo
      enddo

      dRdX(1,1)=(x(1)-x(4))/r(1)
      dRdX(1,2)=(x(2)-x(5))/r(1)
      dRdX(1,3)=(x(3)-x(6))/r(1)
      dRdX(1,4)=-dRdX(1,1)
      dRdX(1,5)=-dRdX(1,2)
      dRdX(1,6)=-dRdX(1,3)

      dRdX(2,4)=(x(4)-x(7))/r(2)
      dRdX(2,5)=(x(5)-x(8))/r(2)
      dRdX(2,6)=(x(6)-x(9))/r(2)
      dRdX(2,7)=-dRdX(2,4)
      dRdX(2,8)=-dRdX(2,5)
      dRdX(2,9)=-dRdX(2,6)

      dRdX(3,1)=(x(1)-x(7))/r(3)
      dRdX(3,2)=(x(2)-x(8))/r(3)
      dRdX(3,3)=(x(3)-x(9))/r(3)
      dRdX(3,7)=-dRdX(3,1)
      dRdX(3,8)=-dRdX(3,2)
      dRdX(3,9)=-dRdX(3,3)

      endsubroutine

      subroutine diagonalize(n,A_ss,EV_s,U_ss)
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: A_ss(n,n)
      real*8, intent(out) :: U_ss(n,n)
      real*8, intent(out) :: EV_s(n)
      integer :: io,i
      real*8,allocatable :: work_d(:)
      integer, save :: lwork_d
      logical, save :: first_time_diag=.true.

      U_ss=A_ss
      if (first_time_diag) then
        lwork_d= -1
        allocate( work_d(2) )
        call dsyev('V','L',n,U_ss,n,EV_s,work_d,lwork_d,io)
        lwork_d=int(work_d(1))
        deallocate ( work_d )
        first_time_diag=.false.
      endif
      allocate( work_d(lwork_d) )
      call dsyev('V','L',n,U_ss,n,EV_s,work_d,lwork_d,io)
      return
      end subroutine diagonalize

C   System:          H3plus
C   Functional form: DIM
C   Common name:     H3+ KBNN with derivatives
C   Number of derivatives: 1
C   Number of bodies: 3
C   Number of electronic surfaces: 3
C   Interface: 3-2V
C
C   References: V. G. Ushakov, K. Nobusada, and V. I. Osherov,
c               Phys. Chem. Chem. Phys., Vol. 3, 63 (2001);
c               H. Kamisaka, W. Bian, K. Nobusada, and H. Nakamura,
c               J. Chem. Phys. Vol. 116, 654 (2002).
C
c   Notes: This is a 3x3 diabatic representation of the h3+ system.
c          This potential routine was coded from the functional forms and
c          parameters as described in the publications listed above
c          (see the References section), but it may not be exactly
c          equivalent to the potential that was used in the dynamics
c          calculations of the Kamisaka paper (private communication, 2002).

c
C     Coded May 28, 2002 by A. G. Anderson, A. W. Jasper, 
C                           and D. G. Truhlar
C
C     definitions for parameters:
C     r -> array with dimensions R(nt,3)
C          R(i,1) is distance between atoms 1 and 2 (A-B)
C          R(i,2) is distance between atoms 2 and 3 (B-C)
C          R(i,3) is distance between atoms 1 and 3 (A-C)
C          where 1 <= i <= nt
C          and distances are in atomic units
C     De -> array of the derivatives of the diabatic potential
C     nt -> number of geometries to be computed. in this program, nt is 
C          always = 1
C     e ->  array dimension nt containing the diabatic potential energy(ies)
C     nsurf -> which surface we want the energy(ies) for. 
C              nsurf = 1: U11
C              nsurf = 2: U12
C              nsurf = 3: U22
C              nsurf = 4: U13
C              nsurf = 5: U23
C              nsurf = 6: U33
C
c     This potential uses the functional forms of 
c     V. G. Ushakov, K. Nobusada, and V. I. Osherov,
c     Phys. Chem. Chem. Phys., 3, 63-69 (2001).
c     and the updated parameters in
c     H. Kamisaka, W. Bian, K. Nobusada, and H. Nakamura,
c     J. Chem. Phys., 116, 654-665 (2002).
C
C                           DISCLAIMER
C
C   This file was generated on 03/21/02 by the version of
C   ADIFOR compiled on June, 1998.
C
C   ADIFOR was prepared as an account of work sponsored by an
C   agency of the United States Government, Rice University, and
C   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
C   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
C   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
C   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
C   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
C   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C******************************************
C      subroutine g_prepot(g_p_)
      subroutine pot(r,  e,  de, nt, nsurf)
C
        implicit double precision (a-h, o-z)
        
        dimension r(nt, 3), de(3, nt), e(nt)
        integer :: nt, nsurf
C
C     Table 2: fitting parameters of eqns 10 and 11
        double precision c1(3), c2(3), c3(3), c4(3), c5(3), c6(3), c7(3)
        double precision c8(3), c9(3), c10(3), c11(3), c12(3), c13(3)
        double precision c14(3)
        double precision p(3, 3)
        double precision c(3, 3)
        double precision osg(3), tsg(3), tsu(3)
        double precision h(3), g(3), u(3)
        double precision potential(3, 3)
C
C     Table 2: fitting params for eqns 10 and 11
        integer g_pmax_
        integer g_i_, g_p_, ldg_r, ldg_e
        parameter (g_pmax_ = 3)
        parameter (ldg_r = 3)
        parameter (ldg_e = 3)
        double precision d3_w, d2_w, d1_w, d4_p, d3_p, d2_p, d1_p, d2_v,
     * d3_v, d4_v
        double precision d38_v, d30_b, d7_b, d8_v, d9_v, d10_v, d41_v,
     * d44_b, d13_v, d42_v
        double precision d15_v, d16_v, d17_v, d43_v, d8_b, d20_v, d16_b,
     * d10_b, d11_b, d24_v
        double precision d12_b, d21_b, d27_v, d28_v, d29_v, d2_b, d3_b, 
     *d32_v, d4_b, d34_v
        double precision d5_b, d6_b, d37_v, g_r1(g_pmax_), g_r(ldg_r, nt
     *, 3), g_r2(g_pmax_), g_r3(g_pmax_), g_p(g_pmax_, 3, 3), g_d1_w(g_p
     *max_), g_c(g_pmax_, 3, 3)
        double precision g_osg(g_pmax_, 3), g_d2_w(g_pmax_), g_d3_w(g_pm
     *ax_), g_tsg(g_pmax_, 3), g_tsu(g_pmax_, 3), g_h(g_pmax_, 3), g_g(g
     *_pmax_, 3), g_u(g_pmax_, 3), g_potential(g_pmax_, 3, 3), g_e(ldg_e
     *, nt)
        integer g_ehfid
        intrinsic dble
        data c1  /-2.03221d0, -0.51486d0, -3.14355d0/
        data c2  / 1.46569d0,  0.08741d0,  1.39001d0/
        data c3  / 1.66412d0, -0.12057d0, -0.11632d0/
        data c4  / 0.27272d0, -0.02033d0, -0.28688d0/
        data c5  / 0.06060d0,  0.08057d0, -0.42105d0/
        data c6  /-2.25049d0, -0.37742d0,  0.44070d0/
        data c7  / 0.01945d0,  0.05521d0,  0.46790d0/
        data c8  /-0.05753d0,  0.01096d0,  0.05605d0/
        data c9  /-0.01288d0, -0.00318d0, -0.05596d0/
        data c10 / 0.57271d0,  0.00257d0,  0.10729d0/
        data c11 /-0.07831d0,  0.00774d0,  0.16146d0/
        data c12 / 0.50812d0,  0.00916d0, -0.09256d0/
        data c13 /-0.07156d0, -0.02141d0, -0.15563d0/
        data c14 / 0.98170d0,  0.74157d0,  0.79817d0/
C
C     Table 1: diatomic parameters for eqns 7-9
        data d_h /-0.17445d0/, d_g /0.10272d0/, d_u /0.62871d0/
        data b_h /2.09607d0/, b_g /0.71935d0/, b_u /0.87080d0/
        data r_h /1.40104d0/, r_g /2.00313d0/, r_u /0.34405d0/
C
        data a /1.13804d0/, b /0.52961d0/
C
        save
C
        data g_ehfid /0/
C
        g_p_ = 3
        
        do i=1,nt
        do j=1,3
        do k=1,3
         if (j .eq. k) g_r(j,i,k)=1.0d0
         if (j .ne. k) g_r(j,i,k)=0.0d0
        enddo
        enddo
        enddo
C
        do ii = 1, nt

          do g_i_ = 1, g_p_
            g_r1(g_i_) = g_r(g_i_, ii, 1)
          enddo
          r1 = r(ii, 1)
C--------
          do g_i_ = 1, g_p_
            g_r2(g_i_) = g_r(g_i_, ii, 2)
          enddo
          r2 = r(ii, 2)
C--------
          do g_i_ = 1, g_p_
            g_r3(g_i_) = g_r(g_i_, ii, 3)
          enddo
          r3 = r(ii, 3)
C--------
C
C     equation 11
          do l = 1, 3
            d4_v = c4(l) * r1
            d10_v = (r2 + r3) * (r2 + r3)
            d4_p = 2.0d0 * (r2 + r3)
            d13_v = c6(l) * r2
            d16_v = c7(l) * r1
            d17_v = r2 + r3
            d20_v = r1 ** ( 3 - 2)
            d20_v =  d20_v * r1
            d3_p =  3 *  d20_v
            d20_v =  d20_v * r1
            d24_v = (r2 + r3) ** ( 3 - 2)
            d24_v =  d24_v * (r2 + r3)
            d2_p =  3 *  d24_v
            d24_v =  d24_v * (r2 + r3)
            d27_v = c10(l) * r2
            d28_v = d27_v * r3
            d29_v = r2 + r3
            d32_v = c11(l) * r1
            d34_v = (r2 + r3) * (r2 + r3)
            d1_p = 2.0d0 * (r2 + r3)
            d37_v = c12(l) * r1
            d38_v = d37_v * r2
            d41_v = c13(l) * r1
            d42_v = d41_v * r1
            d43_v = r2 + r3
            d21_b = d32_v * d1_p
            d30_b = c9(l) * d2_p
            d44_b = c5(l) * d4_p
            d5_b = c3(l) + d42_v + r3 * d37_v + d21_b + d28_v + d29_v * 
     *r3 * c10(l) + d30_b + d16_v + r3 * c6(l) + d44_b
            d6_b = c3(l) + d42_v + d38_v + d21_b + d28_v + d29_v * d27_v
     * + d30_b + d16_v + d13_v + d44_b
            d12_b = d43_v * d41_v + d43_v * r1 * c13(l) + r3 * r2 * c12(
     *l) + d34_v * c11(l) + c8(l) * d3_p + d17_v * c7(l) + d4_v + r1 * c
     *4(l) + c2(l)
            do g_i_ = 1, g_p_
              g_p(g_i_, l, 1) = d6_b * g_r3(g_i_) + d5_b * g_r2(g_i_) + 
     *d12_b * g_r1(g_i_)
            enddo
            p(l, 1) = c1(l) + c2(l) * r1 + d4_v * r1 + c5(l) * d10_v + d
     *13_v * r3 + d16_v * d17_v + c8(l) * d20_v + c9(l) * d24_v + d28_v 
     ** d29_v + d32_v * d34_v + d38_v * r3 + d42_v * d43_v + c3(l) * (r2
     * + r3)
C--------
C
            d4_v = c4(l) * r2
            d10_v = (r3 + r1) * (r3 + r1)
            d4_p = 2.0d0 * (r3 + r1)
            d13_v = c6(l) * r3
            d16_v = c7(l) * r2
            d17_v = r3 + r1
            d20_v = r2 ** ( 3 - 2)
            d20_v =  d20_v * r2
            d3_p =  3 *  d20_v
            d20_v =  d20_v * r2
            d24_v = (r3 + r1) ** ( 3 - 2)
            d24_v =  d24_v * (r3 + r1)
            d2_p =  3 *  d24_v
            d24_v =  d24_v * (r3 + r1)
            d27_v = c10(l) * r3
            d28_v = d27_v * r1
            d29_v = r3 + r1
            d32_v = c11(l) * r2
            d34_v = (r3 + r1) * (r3 + r1)
            d1_p = 2.0d0 * (r3 + r1)
            d37_v = c12(l) * r2
            d38_v = d37_v * r3
            d41_v = c13(l) * r2
            d42_v = d41_v * r2
            d43_v = r3 + r1
            d21_b = d32_v * d1_p
            d30_b = c9(l) * d2_p
            d44_b = c5(l) * d4_p
            d5_b = c3(l) + d42_v + r1 * d37_v + d21_b + d28_v + d29_v * 
     *r1 * c10(l) + d30_b + d16_v + r1 * c6(l) + d44_b
            d6_b = c3(l) + d42_v + d38_v + d21_b + d28_v + d29_v * d27_v
     * + d30_b + d16_v + d13_v + d44_b
            d12_b = d43_v * d41_v + d43_v * r2 * c13(l) + r1 * r3 * c12(
     *l) + d34_v * c11(l) + c8(l) * d3_p + d17_v * c7(l) + d4_v + r2 * c
     *4(l) + c2(l)
            do g_i_ = 1, g_p_
              g_p(g_i_, l, 2) = d6_b * g_r1(g_i_) + d5_b * g_r3(g_i_) + 
     *d12_b * g_r2(g_i_)
            enddo
            p(l, 2) = c1(l) + c2(l) * r2 + d4_v * r2 + c5(l) * d10_v + d
     *13_v * r1 + d16_v * d17_v + c8(l) * d20_v + c9(l) * d24_v + d28_v 
     ** d29_v + d32_v * d34_v + d38_v * r1 + d42_v * d43_v + c3(l) * (r3
     * + r1)
C--------
C
            d4_v = c4(l) * r3
            d10_v = (r1 + r2) * (r1 + r2)
            d4_p = 2.0d0 * (r1 + r2)
            d13_v = c6(l) * r1
            d16_v = c7(l) * r3
            d17_v = r1 + r2
            d20_v = r3 ** ( 3 - 2)
            d20_v =  d20_v * r3
            d3_p =  3 *  d20_v
            d20_v =  d20_v * r3
            d24_v = (r1 + r2) ** ( 3 - 2)
            d24_v =  d24_v * (r1 + r2)
            d2_p =  3 *  d24_v
            d24_v =  d24_v * (r1 + r2)
            d27_v = c10(l) * r1
            d28_v = d27_v * r2
            d29_v = r1 + r2
            d32_v = c11(l) * r3
            d34_v = (r1 + r2) * (r1 + r2)
            d1_p = 2.0d0 * (r1 + r2)
            d37_v = c12(l) * r3
            d38_v = d37_v * r1
            d41_v = c13(l) * r3
            d42_v = d41_v * r3
            d43_v = r1 + r2
            d21_b = d32_v * d1_p
            d30_b = c9(l) * d2_p
            d44_b = c5(l) * d4_p
            d5_b = c3(l) + d42_v + r2 * d37_v + d21_b + d28_v + d29_v * 
     *r2 * c10(l) + d30_b + d16_v + r2 * c6(l) + d44_b
            d6_b = c3(l) + d42_v + d38_v + d21_b + d28_v + d29_v * d27_v
     * + d30_b + d16_v + d13_v + d44_b
            d12_b = d43_v * d41_v + d43_v * r3 * c13(l) + r2 * r1 * c12(
     *l) + d34_v * c11(l) + c8(l) * d3_p + d17_v * c7(l) + d4_v + r3 * c
     *4(l) + c2(l)
            do g_i_ = 1, g_p_
              g_p(g_i_, l, 3) = d6_b * g_r2(g_i_) + d5_b * g_r1(g_i_) + 
     *d12_b * g_r3(g_i_)
            enddo
            p(l, 3) = c1(l) + c2(l) * r3 + d4_v * r3 + c5(l) * d10_v + d
     *13_v * r2 + d16_v * d17_v + c8(l) * d20_v + c9(l) * d24_v + d28_v 
     ** d29_v + d32_v * d34_v + d38_v * r2 + d42_v * d43_v + c3(l) * (r1
     * + r2)
C--------
          enddo
C
C     equation 10
          do l = 1, 3
            d2_b = (-c14(l)) * c14(l)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d2_b * g_r3(g_i_) + d2_b * g_r2(g_i_)
            enddo
            d1_w = (-c14(l)) * c14(l) * (r2 + r3)
            d3_v = exp(d1_w)
            d1_p =  d3_v
            d4_b = p(l, 1) * d1_p
            do g_i_ = 1, g_p_
              g_c(g_i_, l, 1) = d4_b * g_d1_w(g_i_) + d3_v * g_p(g_i_, l
     *, 1)
            enddo
            c(l, 1) = p(l, 1) * d3_v
C--------
            d2_b = (-c14(l)) * c14(l)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d2_b * g_r1(g_i_) + d2_b * g_r3(g_i_)
            enddo
            d1_w = (-c14(l)) * c14(l) * (r3 + r1)
            d3_v = exp(d1_w)
            d1_p =  d3_v
            d4_b = p(l, 2) * d1_p
            do g_i_ = 1, g_p_
              g_c(g_i_, l, 2) = d4_b * g_d1_w(g_i_) + d3_v * g_p(g_i_, l
     *, 2)
            enddo
            c(l, 2) = p(l, 2) * d3_v
C--------
            d2_b = (-c14(l)) * c14(l)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d2_b * g_r2(g_i_) + d2_b * g_r1(g_i_)
            enddo
            d1_w = (-c14(l)) * c14(l) * (r1 + r2)
            d3_v = exp(d1_w)
            d1_p =  d3_v
            d4_b = p(l, 3) * d1_p
            do g_i_ = 1, g_p_
              g_c(g_i_, l, 3) = d4_b * g_d1_w(g_i_) + d3_v * g_p(g_i_, l
     *, 3)
            enddo
            c(l, 3) = p(l, 3) * d3_v
C--------
          enddo
C
C     equations 7,8,9
          do k = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = (-b_h) * g_r(g_i_, ii, k)
            enddo
            d1_w = (-b_h) * (r(ii, k) - r_h)
            d2_v = exp(d1_w)
            d3_p =  d2_v
            d3_v = d2_v * d_h
            d9_v = (r(ii, k) - r_h) * (r(ii, k) - r_h)
            d2_p = 2.0d0 * (r(ii, k) - r_h)
            d13_v = (r(ii, k) - r_h) ** ( 3 - 2)
            d13_v =  d13_v * (r(ii, k) - r_h)
            d1_p =  3 *  d13_v
            d13_v =  d13_v * (r(ii, k) - r_h)
            d15_v = 1.0d0 + b_h * (r(ii, k) - r_h) + a * d9_v + b * d13_
     *v
            d8_b = d3_v * b * d1_p + d3_v * a * d2_p + d3_v * b_h
            d16_b = d15_v * d_h * d3_p
            do g_i_ = 1, g_p_
              g_osg(g_i_, k) = d8_b * g_r(g_i_, ii, k) + d16_b * g_d1_w(
     *g_i_)
            enddo
            osg(k) = d3_v * d15_v
C--------
            d3_b = (-2.0d0) * b_g
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d3_b * g_r(g_i_, ii, k)
            enddo
            d1_w = (-2.0d0) * b_g * (r(ii, k) - r_g)
            d3_b = dble(-1.0) * b_g
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = d3_b * g_r(g_i_, ii, k)
            enddo
            d2_w = dble(-1.0) * b_g * (r(ii, k) - r_g)
            d3_b = (-8.0d0) * b_g
            do g_i_ = 1, g_p_
              g_d3_w(g_i_) = d3_b * g_r(g_i_, ii, k)
            enddo
            d3_w = (-8.0d0) * b_g * (r(ii, k) - r_g)
            d2_v = exp(d1_w)
            d3_p =  d2_v
            d4_v = exp(d2_w)
            d2_p =  d4_v
            d8_v = exp(d3_w)
            d1_p =  d8_v
            d6_b = d_g * 0.0015d0 * d1_p
            d10_b = (-d_g) * 2.0d0 * d2_p
            d11_b = d_g * d3_p
            do g_i_ = 1, g_p_
              g_tsg(g_i_, k) = d6_b * g_d3_w(g_i_) + d10_b * g_d2_w(g_i_
     *) + d11_b * g_d1_w(g_i_)
            enddo
            tsg(k) = d_g * (d2_v - 2.0d0 * d4_v + 0.0015d0 * d8_v)
C--------
            d3_b = (-2.0d0) * b_u
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d3_b * g_r(g_i_, ii, k)
            enddo
            d1_w = (-2.0d0) * b_u * (r(ii, k) - r_u)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = (-b_u) * g_r(g_i_, ii, k)
            enddo
            d2_w = (-b_u) * (r(ii, k) - r_u)
            d2_v = exp(d1_w)
            d2_p =  d2_v
            d4_v = exp(d2_w)
            d1_p =  d4_v
            d6_b = d_u * 2.0d0 * d1_p
            d7_b = d_u * d2_p
            do g_i_ = 1, g_p_
              g_tsu(g_i_, k) = d6_b * g_d2_w(g_i_) + d7_b * g_d1_w(g_i_)
            enddo
            tsu(k) = d_u * (d2_v + 2.0d0 * d4_v)
C--------
          enddo
C
C     eqns 4,5,6
          do i = 1, 3
            do g_i_ = 1, g_p_
              g_h(g_i_, i) = g_c(g_i_, 1, i) + g_osg(g_i_, i)
            enddo
            h(i) = osg(i) + c(1, i)
C--------
            do g_i_ = 1, g_p_
              g_g(g_i_, i) = g_c(g_i_, 2, i) + g_tsg(g_i_, i)
            enddo
            g(i) = tsg(i) + c(2, i)
C--------
            do g_i_ = 1, g_p_
              g_u(g_i_, i) = g_c(g_i_, 3, i) + g_tsu(g_i_, i)
            enddo
            u(i) = tsu(i) + c(3, i)
C--------
          enddo
C
C     this 'potential' 3X3 is the potential matrix
          do g_i_ = 1, g_p_
            g_potential(g_i_, 1, 1) = 0.5d0 * g_u(g_i_, 3) + 0.5d0 * g_g
     *(g_i_, 3) + 0.5d0 * g_u(g_i_, 2) + 0.5d0 * g_g(g_i_, 2) + g_h(g_i_
     *, 1)
          enddo
          potential(1, 1) = h(1) + 0.5d0 * (g(2) + u(2) + g(3) + u(3))
C--------
          do g_i_ = 1, g_p_
            g_potential(g_i_, 2, 2) = 0.5d0 * g_u(g_i_, 3) + 0.5d0 * g_g
     *(g_i_, 3) + 0.5d0 * g_u(g_i_, 1) + 0.5d0 * g_g(g_i_, 1) + g_h(g_i_
     *, 2)
          enddo
          potential(2, 2) = h(2) + 0.5d0 * (g(1) + u(1) + g(3) + u(3))
C--------
          do g_i_ = 1, g_p_
            g_potential(g_i_, 3, 3) = 0.5d0 * g_u(g_i_, 2) + 0.5d0 * g_g
     *(g_i_, 2) + 0.5d0 * g_u(g_i_, 1) + 0.5d0 * g_g(g_i_, 1) + g_h(g_i_
     *, 3)
          enddo
          potential(3, 3) = h(3) + 0.5d0 * (g(1) + u(1) + g(2) + u(2))
C--------
          do g_i_ = 1, g_p_
            g_potential(g_i_, 1, 2) = (-0.5d0) * g_u(g_i_, 3) + 0.5d0 * 
     *g_g(g_i_, 3)
          enddo
          potential(1, 2) = 0.5d0 * (g(3) - u(3))
C--------
          do g_i_ = 1, g_p_
            g_potential(g_i_, 1, 3) = (-0.5d0) * g_u(g_i_, 2) + 0.5d0 * 
     *g_g(g_i_, 2)
          enddo
          potential(1, 3) = 0.5d0 * (g(2) - u(2))
C--------
          do g_i_ = 1, g_p_
            g_potential(g_i_, 2, 1) = (-0.5d0) * g_u(g_i_, 3) + 0.5d0 * 
     *g_g(g_i_, 3)
          enddo
          potential(2, 1) = 0.5d0 * (g(3) - u(3))
C--------
          do g_i_ = 1, g_p_
            g_potential(g_i_, 2, 3) = (-0.5d0) * g_u(g_i_, 1) + 0.5d0 * 
     *g_g(g_i_, 1)
          enddo
          potential(2, 3) = 0.5d0 * (g(1) - u(1))
C--------
          do g_i_ = 1, g_p_
            g_potential(g_i_, 3, 1) = (-0.5d0) * g_u(g_i_, 2) + 0.5d0 * 
     *g_g(g_i_, 2)
          enddo
          potential(3, 1) = 0.5d0 * (g(2) - u(2))
C--------
          do g_i_ = 1, g_p_
            g_potential(g_i_, 3, 2) = (-0.5d0) * g_u(g_i_, 1) + 0.5d0 * 
     *g_g(g_i_, 1)
          enddo
          potential(3, 2) = 0.5d0 * (g(1) - u(1))
C--------
C
          if (nsurf .eq. 1) then
            do g_i_ = 1, g_p_
              g_e(g_i_, ii) = g_potential(g_i_, 1, 1)
            enddo
            e(ii) = potential(1, 1) - d_h
C--------
          endif
          if (nsurf .eq. 2) then
            do g_i_ = 1, g_p_
              g_e(g_i_, ii) = g_potential(g_i_, 1, 2)
            enddo
            e(ii) = potential(1, 2)
C--------
          endif
          if (nsurf .eq. 3) then
            do g_i_ = 1, g_p_
              g_e(g_i_, ii) = g_potential(g_i_, 2, 2)
            enddo
            e(ii) = potential(2, 2) - d_h
C--------
          endif
          if (nsurf .eq. 4) then
            do g_i_ = 1, g_p_
              g_e(g_i_, ii) = g_potential(g_i_, 1, 3)
            enddo
            e(ii) = potential(1, 3)
C--------
          endif
          if (nsurf .eq. 5) then
            do g_i_ = 1, g_p_
              g_e(g_i_, ii) = g_potential(g_i_, 2, 3)
            enddo
            e(ii) = potential(2, 3)
C--------
          endif
          if (nsurf .eq. 6) then
            do g_i_ = 1, g_p_
              g_e(g_i_, ii) = g_potential(g_i_, 3, 3)
            enddo
            e(ii) = potential(3, 3) - d_h
C--------
          endif
          if (nsurf .gt. 6 .or. nsurf .lt. 1) then
            do g_i_ = 1, g_p_
              g_e(g_i_, ii) = 0.0d0
            enddo
            e(ii) = 0.d0
C--------
            write (6, *) 'WARN:  NSURF = ', nsurf
          endif

C
C derivatives
          do i = 1, 3
            de(i, ii) = g_e(i,ii)
          enddo

C
C end nt loop
        enddo
C
      end
C******************************************
C
C
C
