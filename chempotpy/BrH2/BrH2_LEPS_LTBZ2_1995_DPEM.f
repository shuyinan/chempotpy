      subroutine dpem(x,igrad,u,ug)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: u(nstates,nstates)
      double precision, intent(out) :: ug(nstates,nstates,natoms,3)

      double precision :: nt, r(1,3), r2(3), v(1)
      double precision :: u11(1), u12(1), u22(1)
      double precision :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
      double precision :: dudr(nstates,nstates,3)
      double precision :: dudx(nstates,nstates,9)
      double precision :: t(nstates,nstates)
      double precision :: tmpmat(nstates,nstates)
      double precision :: tx(9), drdx(3,9)
      logical, save :: first_time_data=.true.
      integer :: i, j, k, l

      !initialize 
      u=0.d0
      p=0.d0
      g=0.d0
      d=0.d0

      nt=1
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        tx(j)=x(iatom,idir)
      enddo
      enddo
      ! input cartesian is ClHH
      r(1,1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      r(1,2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      r(1,3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211

      NDER=igrad

      if(first_time_data) then
      call prepot
      first_time_data=.false.
      endif

      call pot(r,u11,de_u11,1,1)
      call pot(r,u12,de_u12,1,2)
      call pot(r,u22,de_u22,1,3)
      u(1,1)=u11(1)*27.211386
      u(1,2)=u12(1)*27.211386
      u(2,1)=u12(1)*27.211386
      u(2,2)=u22(1)*27.211386
      dudr(1,1,:)=de_u11(:,1)*51.422067
      dudr(1,2,:)=de_u12(:,1)*51.422067
      dudr(2,1,:)=de_u12(:,1)*51.422067
      dudr(2,2,:)=de_u22(:,1)*51.422067
      
      r2=r(1,:)*0.529177211
      call evdrdx(tx, r2, drdx)
      dudx=0.d0
      do i=1,9
      do j=1,3
        dudx(:,:,i)=dudx(:,:,i)+dudr(:,:,j)*drdx(j,i)
      enddo
      enddo
  
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        ug(:,:,iatom,idir)=dudx(:,:,j)
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

C
      subroutine prepot                                               
C
C
C   System:          BrH2
C   Functional form: extended-LEPS plus three-center term
C   Common name:     BrH2SECD
C   Interface:       3-2V
C   Number of electronic surfaces: 2
C   Number of derivatives: 1
C   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g. Zhao
C                    J. Phys. Chem. 99, 207 (1995)
C
C   Note:            This is a vectorized version of the BRH2SEC surface
C                    with derivatives
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy and the derivatives 
C               of the energy with respect to the coordinates for a given 
C               geometric configuration
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohr
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground and first excited electronic states, diabatic representation
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Br
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Br-H)
C                            R(2) = R(H-H)
C                            R(3) = R(Br-H)
C
C
C
      implicit none                                                     08Y20V98
      integer g_p_                                                      08Y20V98
      parameter (g_p_=3)                                                08Y20V98
      call g_prepot11(g_p_)                                             08Y20V98 
      call prepot12
      call prepot22

      return                                                           
      end

c**********************************************************************
c**********************************************************************

      subroutine pot(r,e,de,nt,nsurf)

      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt),de(3,nt)

      integer i,j,k,g_p_,ldg_r,nsurf,nt                                 08Y20V98
      parameter (g_p_=3,ldg_e =3)                                       08Y20V98
      dimension g_e(ldg_e,nt)                                           08Y20V98

      if (nsurf .eq. 1) then
         call g_pot11(g_p_, r, e, De, ldg_e, nt)                        08Y20V98
      else if (nsurf .eq. 2) then
         call pot12(r,e,De,nt)
      else if (nsurf .eq. 3) then
         call pot22(r,e,De,nt)
      end if

      return                                                            4s19m94
      end                                                               4s19m94


      subroutine g_prepot11(g_p_)                                       08Y20V98
C   this is a vectorized version of the BrH2SEC potential with derivatives
C
C   System:          BrH2
C   Functional form: extended-LEPS plus three-center term with
C                    internal-angle-dependent Sato parameters
C   Common name:     BrH2SECD
C   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, et. al.
C                    J. Phys. Chem. 99, 207-225 (1995)
C   Analytic derivatives added, Yuri Volobuev, July-Aug 1998
C   Derivatives for pot11 are obtained using Adifor 2.1
C                    
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in DATA statements.
C   For the this potential energy surface the potential energy in the 
C   BC valley is equal to the energy of the HH diatomic. 
C   The energy in the AC valley is equal to the energy of the HBr diatomic,
C   which is also the energy value in the BA valley. 
C
C   This potential is written such that:
C                  R(1) = R(Br-H)
C                  R(2) = R(H-H)
C                  R(3) = R(H-Br)
C   The zero of energy is defined at Br "infinitely" far from the H2 diatomic.
C
      implicit double precision (a-h, o-z)
      double precision r(3)
C
      parameter (cangau = 0.52917706d0)
      parameter (charkc = 627.5095d0)
C
      common /satocm/ zhhp1, zhhp2, zhbrp1, zhbrp2, z(3), dz(3, 3)
      common /v3ccom/ v3cp1, v3cp2, v3cp3, v3cp4, v3cp5, v3cp6, v3c, 
     &                  dv3c(3)
      common /lengcm/ x(3), coul(3), exch(3), rad
      common /lsatcm/ zpo(3), op3z(3), zp3(3), tzp3(3), top3z(3), 
     &                  do4z(3), beta(3)
      common /ldercm/ dxdr(3, 3), deq(3), derad(3)
      common /decom/ de(3)
C
      dimension d(3), re(3), betaa(3), b(3), rr(n, 3), e(n)
C
      save
C
      integer g_pmax_
      parameter (g_pmax_ = 3)
      integer g_i_, g_p_, ldg_rr, ldg_e
      double precision d3_p, d2_p, d2_v, d10_v, d2_b, d3_v, d4_v, 
     & d5_v, d6_v, d7_v, d3_b, d4_b, d5_b, d6_b, d7_b, d1_p, d8_v, 
     & d8_b, d1_w, d10_b, g_r(g_pmax_, 3), g_r1p2(g_pmax_), 
     & g_r2p2(g_pmax_), g_r3p2(g_pmax_), g_zchi(g_pmax_), 
     & g_capr2(g_pmax_), g_capr(g_pmax_), g_coschi(g_pmax_), 
     & g_schip2(g_pmax_), g_z(g_pmax_, 3), g_zpo(g_pmax_, 3), 
     & g_op3z(g_pmax_, 3), g_top3z(g_pmax_, 3), g_zp3(g_pmax_, 3), 
     & g_tzp3(g_pmax_, 3), g_do4z(g_pmax_, 3), g_pe(g_pmax_), 
     & g_d1_w(g_pmax_), g_x(g_pmax_, 3), g_coul(g_pmax_, 3), 
     & g_exch(g_pmax_, 3), g_rad(g_pmax_), g_r1p3(g_pmax_), 
     & g_r1m3(g_pmax_), g_r1p3m2(g_pmax_), g_r1m3p2(g_pmax_), 
     & g_r132p2(g_pmax_), g_dist(g_pmax_), g_costh(g_pmax_),
     & g_sintp2(g_pmax_), g_sintp6(g_pmax_), g_sintp8(g_pmax_), 
     & g_fcn4(g_pmax_), g_v3c(g_pmax_), g_e(ldg_e, *), g_zhhp1(g_pmax_),
     & g_zhhp2(g_pmax_), g_zhbrp1(g_pmax_), g_zhbrp2(g_pmax_), 
     & g_dz(g_pmax_, 3, 3), g_beta(g_pmax_, 3), g_v3cp1(g_pmax_), 
     & g_v3cp2(g_pmax_), g_v3cp3(g_pmax_), g_v3cp4(g_pmax_), 
     & g_v3cp5(g_pmax_), g_v3cp6(g_pmax_), g_dv3c(g_pmax_, 3)

      common /com_para/ g_myrank,g_nprocs                               08Y24V98
      integer g_myrank,g_nprocs                                         08Y24V98


      common /g_satocm/ g_zhhp1, g_zhhp2, g_zhbrp1, g_zhbrp2, g_z, g_dz
      common /g_v3ccom/ g_v3cp1, g_v3cp2, g_v3cp3, g_v3cp4, g_v3cp5, 
     * g_v3cp6, g_v3c, g_dv3c
      common /g_lengcm/ g_x, g_coul, g_exch, g_rad
      common /g_lsatcm/ g_zpo, g_op3z, g_zp3, g_tzp3, g_top3z, g_do4z,
     * g_beta
      save g_r1m3, g_r1p3m2, g_r1m3p2, g_r132p2, g_dist, g_costh, 
     & g_sintp2, g_sintp6, g_sintp8, g_fcn4, g_r1p2, g_r2p2, g_r3p2, 
     & g_zchi, g_capr2, g_capr, g_coschi, g_schip2, g_d1_w, g_r1p3

      data d /90.4473d0, 109.559d0, 90.4473d0/

      if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
      endif

        xlamb = .4568d0 / 27.21161d0
        a1 = v3cp1
        a2 = v3cp2
        a3 = v3cp3
        a4 = v3cp4
        a5 = v3cp5
        a6 = v3cp6
        data re /1.41443d0, 0.74144d0, 1.41443d0/
        data betaa /1.81088d0, 1.94198d0, 1.81088d0/
        data r2 /1.41421356d0/
C
        zhhp1 = 2.91230499d-1
        zhhp2 = 4.99160974d0
        zhbrp1 = 3.84262775d-3
        zhbrp2 = 2.87642885d-3
        v3cp1 = 1.67457884d0
        v3cp2 = 7.09518875d-1
        v3cp3 = 1.87582679d-1
        v3cp4 = 2.43517074d-1
        v3cp5 = 2.04145206d-1
        v3cp6 = 2.57833022d0
C
C                    
C   Convert the constants to atomic units
C
        do 99999 i = 1, 3
          de(i) = d(i) / charkc
          re(i) = re(i) / cangau
          beta(i) = betaa(i) * cangau
99999   continue
C
100   FORMAT(/,2X,T5,'PREPOT has been called for the BrH2 potential',
     *               ' SEC',
     *       /,2X,T5,'Potential parameters:',
     *       /,2X,T5,'Bond:',T38,'BrH',T49,'HH',T58,'HBr',
     *       /,2X,T5,'Dissociation energies:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Equilibrium bond lengths:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Morse betas:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for ZHH: ',T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for ZHBr: ',T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for V3c: ',5(T45,F10.5,T55,F10.5,/))
C

C
        return
C
        entry g_pot11(g_p_, rr, e, g_e, ldg_e, n)                       08Y20V98
C
        do 99996 k = 1, n
          do g_i_ = 1, g_p_
            g_r(g_i_, 1) = 0.0d0
          enddo
          r(1) = rr(k, 1)
C--------
          do g_i_ = 1, g_p_
            g_r(g_i_, 2) = 0.0d0
          enddo
          r(2) = rr(k, 2)
C--------
          do g_i_ = 1, g_p_
            g_r(g_i_, 3) = 0.0d0
          enddo
          r(3) = rr(k, 3)

          g_r(1,1) = 1.0d0
          g_r(2,2) = 1.0d0
          g_r(3,3) = 1.0d0
C--------
C
C   Compute bond dependent Sato parameters
          d2_b = r(1) + r(1)
          do g_i_ = 1, g_p_
            g_r1p2(g_i_) = d2_b * g_r(g_i_, 1)
          enddo
          r1p2 = r(1) * r(1)
C--------
          d2_b = r(2) + r(2)
          do g_i_ = 1, g_p_
            g_r2p2(g_i_) = d2_b * g_r(g_i_, 2)
          enddo
          r2p2 = r(2) * r(2)
C--------
          d2_b = r(3) + r(3)
          do g_i_ = 1, g_p_
            g_r3p2(g_i_) = d2_b * g_r(g_i_, 3)
          enddo
          r3p2 = r(3) * r(3)
C--------
          r2i = 1.0d0 / r(2)
          do g_i_ = 1, g_p_
            g_zchi(g_i_) = 0.0d0
          enddo
          zchi = 0.0d0
C--------
C   CAPR IS THE Br TO H2 DISTANCE
          d5_b = (-0.5d0) * 0.5d0
          do g_i_ = 1, g_p_
            g_capr2(g_i_) = d5_b * g_r2p2(g_i_) + 0.5d0 * g_r3p2(g_i_) +
     * 0.5d0 * g_r1p2(g_i_)
          enddo
          capr2 = 0.5d0 * (r1p2 + r3p2 - 0.5d0 * r2p2)
C--------
          if (capr2 .le. 0.d0) then
            do g_i_ = 1, g_p_
              g_capr2(g_i_) = 0.0d0
            enddo
            capr2 = 0.0d0
C--------
            do g_i_ = 1, g_p_
              g_capr(g_i_) = 0.0d0
            enddo
            capr = 0.0d0
C--------
          else
            d2_v = sqrt(capr2)

            if ( capr2 .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               d1_p = 0.0d0
            endif
            do g_i_ = 1, g_p_
              g_capr(g_i_) = d1_p * g_capr2(g_i_)
            enddo
            capr = d2_v
C--------
C   Compute Chi, the angle between CAPR and RHH
            d7_v = r(2) * capr
            d8_v = 0.5d0 * (r3p2 - r1p2) / d7_v
            d3_b = (-d8_v) / d7_v
            d4_b = d3_b * capr
            d5_b = d3_b * r(2)
            d6_b = 1.0d0 / d7_v * 0.5d0
            do g_i_ = 1, g_p_
              g_coschi(g_i_) = d5_b * g_capr(g_i_) + d4_b * g_r(g_i_, 2)
     * + (-d6_b) * g_r1p2(g_i_) + d6_b * g_r3p2(g_i_)
            enddo
            coschi = d8_v
C--------
            d3_b = -coschi + (-coschi)
            do g_i_ = 1, g_p_
              g_schip2(g_i_) = d3_b * g_coschi(g_i_)
            enddo
            schip2 = 1.0d0 - coschi * coschi
C--------
            d4_v = zhbrp2 * schip2
            d6_v = zhbrp1 * schip2 + d4_v * schip2
            d7_b = capr2 * d4_v + capr2 * schip2 * zhbrp2 + capr2 * zhbr
     *p1
            do g_i_ = 1, g_p_
              g_zchi(g_i_) = d7_b * g_schip2(g_i_) + d6_v * g_capr2(g_i_
     *)
            enddo
            zchi = capr2 * d6_v
C--------
          endif
          do g_i_ = 1, g_p_
            g_z(g_i_, 1) = g_zchi(g_i_)
          enddo
          z(1) = 0.1675d0 + zchi
C--------
          do g_i_ = 1, g_p_
            g_z(g_i_, 3) = g_z(g_i_, 1)
          enddo
          z(3) = z(1)
C--------
          d3_v = zhhp2 + r2p2
          d4_v = zhhp1 * r2p2 / d3_v
          d4_b = (-d4_v) / d3_v + 1.0d0 / d3_v * zhhp1
          do g_i_ = 1, g_p_
            g_z(g_i_, 2) = d4_b * g_r2p2(g_i_)
          enddo
          z(2) = d4_v
C--------
C
C   Compute the LEPS parameters
C
          do 99998 i = 1, 3
            do g_i_ = 1, g_p_
              g_zpo(g_i_, i) = g_z(g_i_, i)
            enddo
            zpo(i) = 1.0d0 + z(i)
C--------
            do g_i_ = 1, g_p_
              g_op3z(g_i_, i) = 3.0d0 * g_z(g_i_, i)
            enddo
            op3z(i) = 1.0d0 + 3.0d0 * z(i)
C--------
            do g_i_ = 1, g_p_
              g_top3z(g_i_, i) = 2.0d0 * g_op3z(g_i_, i)
            enddo
            top3z(i) = 2.0d0 * op3z(i)
C--------
            do g_i_ = 1, g_p_
              g_zp3(g_i_, i) = g_z(g_i_, i)
            enddo
            zp3(i) = z(i) + 3.0d0
C--------
            do g_i_ = 1, g_p_
              g_tzp3(g_i_, i) = 2.0d0 * g_zp3(g_i_, i)
            enddo
            tzp3(i) = 2.0d0 * zp3(i)
C--------
            d2_v = de(i) / 4.0d0 / zpo(i)
            d2_b = (-d2_v) / zpo(i)
            do g_i_ = 1, g_p_
              g_do4z(g_i_, i) = d2_b * g_zpo(g_i_, i)
            enddo
            do4z(i) = d2_v
C--------
            b(i) = beta(i) * do4z(i) * 2.0d0
20          continue
99998     continue

          do g_i_ = 1, g_p_
            g_pe(g_i_) = 0.0d0
          enddo
          pe = 0.d0
C--------
          do 99997 i = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = (-beta(i)) * g_r(g_i_, i)
            enddo
            d1_w = (-beta(i)) * (r(i) - re(i))
            d2_v = exp(d1_w)
            d1_p =  d2_v
            do g_i_ = 1, g_p_
              g_x(g_i_, i) = d1_p * g_d1_w(g_i_)
            enddo
            x(i) = d2_v
C--------
            d6_v = zp3(i) * x(i) - top3z(i)
            d7_v = do4z(i) * d6_v
            d4_b = x(i) * d6_v
            d5_b = x(i) * do4z(i)
            d8_b = d5_b * x(i)
            d3_b = d7_v + d5_b * zp3(i)
            do g_i_ = 1, g_p_
              g_coul(g_i_, i) = (-d5_b) * g_top3z(g_i_, i) + d3_b * g_x(
     *g_i_, i) + d8_b * g_zp3(g_i_, i) + d4_b * g_do4z(g_i_, i)
            enddo
            coul(i) = d7_v * x(i)
C--------
            d6_v = op3z(i) * x(i) - tzp3(i)
            d7_v = do4z(i) * d6_v
            d4_b = x(i) * d6_v
            d5_b = x(i) * do4z(i)
            d8_b = d5_b * x(i)
            d3_b = d7_v + d5_b * op3z(i)
            do g_i_ = 1, g_p_
              g_exch(g_i_, i) = (-d5_b) * g_tzp3(g_i_, i) + d3_b * g_x(g
     *_i_, i) + d8_b * g_op3z(g_i_, i) + d4_b * g_do4z(g_i_, i)
            enddo
            exch(i) = d7_v * x(i)
C--------
            do g_i_ = 1, g_p_
              g_pe(g_i_) = g_coul(g_i_, i) + g_pe(g_i_)
            enddo
            pe = pe + coul(i)
C--------
30          continue
99997     continue
          d4_v = (exch(1) - exch(2)) * (exch(1) - exch(2))
          d3_p = 2.0d0 * (exch(1) - exch(2))
          d7_v = (exch(2) - exch(3)) * (exch(2) - exch(3))
          d2_p = 2.0d0 * (exch(2) - exch(3))
          d10_v = (exch(3) - exch(1)) * (exch(3) - exch(1))
          d1_p = 2.0d0 * (exch(3) - exch(1))
          d5_b = d1_p + (-d2_p)
          d6_b = -d1_p + d3_p
          d10_b = d2_p + (-d3_p)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_exch(g_i_, 3) + d10_b * g_exch(g_i_,
     * 2) + d6_b * g_exch(g_i_, 1)
          enddo
          d1_w = d4_v + d7_v + d10_v
          d2_v = sqrt(d1_w)

          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_rad(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          rad = d2_v
C--------
          d5_b = -(1.0d0 / r2)
          do g_i_ = 1, g_p_
            g_pe(g_i_) = d5_b * g_rad(g_i_) + g_pe(g_i_)
          enddo
          pe = pe - rad / r2 + de(2)
C--------
C   Compute the three-center term
          do g_i_ = 1, g_p_
            g_r1p3(g_i_) = g_r(g_i_, 3) + g_r(g_i_, 1)
          enddo
          r1p3 = r(1) + r(3)
C--------
          a1 = v3cp1
          a2 = v3cp2
          a3 = v3cp3
          a4 = v3cp4
          a5 = v3cp5
          a6 = v3cp6
          do g_i_ = 1, g_p_
            g_r1m3(g_i_) = -g_r(g_i_, 3) + g_r(g_i_, 1)
          enddo
          r1m3 = r(1) - r(3)
C--------
          do g_i_ = 1, g_p_
            g_r1p3m2(g_i_) = -g_r(g_i_, 2) + g_r1p3(g_i_)
          enddo
          r1p3m2 = r1p3 - r(2)
C--------
          d2_b = r1m3 + r1m3
          do g_i_ = 1, g_p_
            g_r1m3p2(g_i_) = d2_b * g_r1m3(g_i_)
          enddo
          r1m3p2 = r1m3 * r1m3
C--------
          d2_b = r1p3m2 + r1p3m2
          do g_i_ = 1, g_p_
            g_r132p2(g_i_) = d2_b * g_r1p3m2(g_i_)
          enddo
          r132p2 = r1p3m2 * r1p3m2
C--------
C   Compute the bond angle from the three internal coordinates
          d4_b = -r(2) + (-r(2))
          d7_b = r(3) + r(3)
          d8_b = r(1) + r(1)
          do g_i_ = 1, g_p_
            g_dist(g_i_) = d4_b * g_r(g_i_, 2) + d7_b * g_r(g_i_, 3) + d
     *8_b * g_r(g_i_, 1)
          enddo
          dist = r(1) * r(1) + r(3) * r(3) - r(2) * r(2)
C--------
          d3_v = 2.0d0 * r(1)
          d5_v = d3_v * r(3)
          d6_v = dist / d5_v
          d2_b = 1.0d0 / d5_v
          d3_b = (-d6_v) / d5_v
          d5_b = d3_b * d3_v
          d6_b = d3_b * r(3) * 2.0d0
          do g_i_ = 1, g_p_
            g_costh(g_i_) = d5_b * g_r(g_i_, 3) + d6_b * g_r(g_i_, 1) + 
     *d2_b * g_dist(g_i_)
          enddo
          costh = d6_v
C--------
          d3_b = -costh + (-costh)
          do g_i_ = 1, g_p_
            g_sintp2(g_i_) = d3_b * g_costh(g_i_)
          enddo
          sintp2 = 1.0d0 - costh * costh
C--------
          d2_v = sintp2 * sintp2
          d3_b = d2_v + sintp2 * sintp2 + sintp2 * sintp2
          do g_i_ = 1, g_p_
            g_sintp6(g_i_) = d3_b * g_sintp2(g_i_)
          enddo
          sintp6 = d2_v * sintp2
C--------
          do g_i_ = 1, g_p_
            g_sintp8(g_i_) = sintp6 * g_sintp2(g_i_) + sintp2 * g_sintp6
     *(g_i_)
          enddo
          sintp8 = sintp6 * sintp2
C--------
          do g_i_ = 1, g_p_
            g_fcn4(g_i_) = a6 * g_sintp8(g_i_) + a5 * g_sintp2(g_i_)
          enddo
          fcn4 = 1.d0 + a5 * sintp2 + a6 * sintp8
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-a4) * g_r132p2(g_i_) + (-a3) * g_r1m3p2(g_i
     *_) + (-a2) * g_r1p3(g_i_)
          enddo
          d1_w = -(a2 * r1p3 + a3 * r1m3p2 + a4 * r132p2)
          d2_v = a1 * fcn4
          d4_v = exp(d1_w)
          d1_p =  d4_v
          d4_b = d2_v * d1_p
          d5_b = d4_v * a1
          do g_i_ = 1, g_p_
            g_v3c(g_i_) = d4_b * g_d1_w(g_i_) + d5_b * g_fcn4(g_i_)
          enddo
          v3c = d2_v * d4_v
C--------
C   Sum all the terms that make up the energy
          do g_i_ = 1, g_p_
            g_pe(g_i_) = g_v3c(g_i_) + g_pe(g_i_)
          enddo
          pe = pe + v3c
C--------
          do g_i_ = 1, g_p_
            g_e(g_i_, k) = g_pe(g_i_)
          enddo
          e(k) = pe + xlamb / 3.d0
C--------
99996   continue
C
        return
      end

c**********************************************************************
c**********************************************************************

      SUBROUTINE PREPOT12
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         dimension rr(n,3),e(*),de(3,n)
         save
          xlamb=0.4568d0/27.21161d0
         RETURN
C
      ENTRY POT12(rr,e,de,n)
      do 7777 k=1,n
        e(k)=-dsqrt(2.d0)*xlamb/3.d0
        do j =1,3                                                       08Y20V98
          de(j,k) = 0.0d0                                               08Y20V98
        enddo                                                           08Y20V98
7777  continue

      RETURN
      END

c
c
      SUBROUTINE PREPOT22
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
         double precision r(3)
         PARAMETER(CANGAU = 0.52917706D0)
         PARAMETER(CHARKC = 627.5095D0)
C
         COMMON /SATOCM/ ZHHP1, ZHHP2, ZHBRP1, ZHBRP2, 
     *                   Z(3), DZ(3,3)
         COMMON /LENGCM/X(3), COUL(3), EXCH(3), RAD
         COMMON /LSATCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),
     *                   DO4Z(3), BETA(3)
         COMMON /DECOM/ DE(3)
         common /com_para/ g_myrank,g_nprocs                            08Y20V98
         integer g_myrank,g_nprocs                                      08Y20V98

C
         DIMENSION D(3), RE(3), BETAA(3), B(3),rr(n,3),e(*),der(3,n)
C
         SAVE
C
         DATA D/90.4473D0, 109.559D0, 90.4473D0/
         DATA RE/1.41443D0, 0.74144D0, 1.41443D0/
         DATA BETAA/1.81088D0, 1.94198D0, 1.81088D0/
         DATA R2/1.41421356D0/
C
         data ccon /0.5d0/
         xlamb=0.4568d0/27.21161d0
         ZHHP1  =  2.91230499D-1
         ZHHP2  =  4.99160974D0
         ZHBRP1 =  3.84262775D-3
         ZHBRP2 =  2.87642885D-3
C
C                    
C   Convert the constants to atomic units
C
         DO 10 I = 1, 3 
               DE(I)   = D(I) / CHARKC
               RE(I)   = RE(I) / CANGAU
               BETA(I) = BETAA(I) * CANGAU
 10      CONTINUE
C
100   FORMAT(/,2X,T5,'PREPOT22 has been called for the BrH2 potential',
     *               ' SEC',
     *       /,2X,T5,'Potential parameters:',
     *       /,2X,T5,'Bond:',T38,'BrH',T49,'HH',T58,'HBr',
     *       /,2X,T5,'Dissociation energies:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Equilibrium bond lengths:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Morse betas:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for ZHH: ',T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for ZHBr: ',T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for V3c: ',5(T45,F10.5,T55,F10.5,/))
C
         RETURN
C
      ENTRY POT22(rr,e,der,n)
      do 7777 k=1,n
       r(1)=rr(k,1)
       r(2)=rr(k,2)
       r(3)=rr(k,3)
       y2=exp(-beta(2)*(r(2)-re(2)))
       y1=exp(-beta(1)*(r(1)-re(1)))             
       y3=exp(-beta(3)*(r(3)-re(3)))             
       vmorse=de(2)*(y2**2-2.d0*y2)
       vanti1=de(1)*(y1**2+2.d0*y1)
       vanti3=de(3)*(y3**2+2.d0*y3)
       pe=de(2)+vmorse+(vanti1+vanti3)*ccon+2.d0*xlamb/3.d0
       der(1,k) = -2.0d0*beta(1)*de(1)*(y1**2+y1)*ccon                  08Y20V98
       der(2,k) = -2.0d0*beta(2)*de(2)*(y2**2-y2)                       08Y20V98
       der(3,k) = -2.0d0*beta(3)*de(3)*(y3**2+y3)*ccon                  08Y20V98
       e(k)=pe
       
7777  continue
c
         RETURN
         END
