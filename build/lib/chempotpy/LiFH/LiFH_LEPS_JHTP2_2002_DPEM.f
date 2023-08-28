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
      call diagonalize(nstates,u,p,t)

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
C   System:          LiFH
C   Functional form: Modified Extended LEPS (2x2 diabatic fit)
C   Common name:     LiFH surface fit J
C   Number of derivatives: 1
C   Number of electronic surfaces: 2
C   Interface: 3-2V
C
C   References:      A. W. Jasper, M. D. Hack, D. G. Truhlar, and P. Piecuch, 
C                    J. Chem. Phys., Vol. 116, 8353 (2002).
C
c   Notes:    This fit contains long-range forces.
C
C   Protocol:
C      PREPOT - initializes the potential's variables
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy
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
C      The classical potential energy is set equal to zero for the Li
C      infinitely far from the HF diatomic and R(HF) set equal to the
C      HF equilibrium diatomic value.
C
C   Parameters:
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Li-H)
C                            R(2) = R(H-F)
C                            R(3) = R(Li-F)
C
c
c  The potential is called by first issuing a call to PREPOT which initializes
c  all surfaces.  The potential is called by calling POT with the parameters
c  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of
c  dimension nt which is an integer representing the number of geometries
c  input.  The integer nsurf selects the desired surface.  nsurf=1 is the
c  lower diabatic surface, nsurf=2 is the coupling surface, nsurf=3 is
c  the upper diabatic surface.
c
c ********************************************************************
c ********************************************************************
c
c                   LiFH-J POTENTIAL ENERGY MATRIX
c
c  Created by A. Jasper on July 27, 2001.  
c  Based on the LiFH-H potential energy matrix created by M. Hack.
c  Based on the ab initio calculations of P. Piecuch.
c
c  Reference: A. W. Jasper, M. D. Hack, D. G. Truhlar, and P. Piecuch, 
c  J. Chem. Phys., prepared for publication.
c
c ********************************************************************
c ********************************************************************

      subroutine prepot

        call prepot11
        call prepot12
        call prepot22

      return

      end

c ********************************************************************
c ********************************************************************

c      subroutine pot(r,e,nt,nsurf)
      subroutine pot(r,e,de,nt,nsurf)

        implicit none
        double precision r, e, de, g_r
        integer i,j,k,nt,nsurf
        dimension r(nt,3),e(nt),g_r(3,nt,3),de(3,nt)

c       seed matrix
        do i=1,nt
        do j=1,3
        do k=1,3
         if (j .eq. k) g_r(j,i,k)=1.0d0
         if (j .ne. k) g_r(j,i,k)=0.0d0
        enddo
        enddo
        enddo

        if (nsurf.eq.1) call pot11(3,r,g_r,3,e,de,3,nt)
        if (nsurf.eq.2) call pot12(3,r,g_r,3,e,de,3,nt)
        if (nsurf.eq.3) call pot22(3,r,g_r,3,e,de,3,nt)
        if (nsurf.gt.3.or.nsurf.lt.1) 
     &     write(6,*)' WARNING:  nsurf = ',nsurf,
     &           ' in LiFHJ potential routine!'

      return
      end

c ********************************************************************
c ********************************************************************

C                           DISCLAIMER
C
C   This file was generated on 07/27/01 by the version of
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
C                            U11 SURFACE
C
      subroutine prepot11
        implicit none
C
C MISC VARIABLES
        integer i, nt, j
        double precision cevau
        double precision r(nt, 3), e(nt)
        double precision r1s, r2s, r3s
C
C INTERACTION REGION VARIABLES 
C     LiH
        double precision s1, relih, bes1, s1c2, s1c3, s1a, xelih, bes1b,
     * s1c, rcut2del, rcut2rho, rcut2, r1r3, r1r3al
C
C     HF
        double precision s2, behf, yhf, rehf, dehf, xehf, s2m, behfm, de
     *hfm, rehfm, xehfm, s2rho, s2alp, dehfc, dehfcc, yhf_r, behf_0, beh
     *f_1, behf_2
C
C     LiF
        double precision s3, res3, bes3, des3, res3c, bes3c, s3c, s3a, x
     *es3, des3c, rcutdel, rcutrho, rcut, r1r3al2, bes3i, bes3z, bes3r
C
C LONG RANGE FORCES VARIABLES
C     HF physical properties
        double precision hfdip_y, hfdip_a, hfdip_re, hfdip_m(11), hfdip,
     * hfpol, hfpolperp, hfpolpar, hfie, hfpolavg
C
C     Li physical properties
        double precision li2spol, li2sie
C
C     misc long-range variables
        double precision rhf0did, rhf0disp, hf_did, hf_disp, e_lr_hf
        double precision mass(3), cauang, caudeb, costheta, temp, bigr, 
     *coschi, r2lrcut_a, r2lrcut_r
C
C MISC DATA
        integer g_pmax_, g_p_
        parameter (g_pmax_ = 3)
        integer g_i_, ldg_r, ldg_e
        double precision d12_b, d12_v, d1_w, d3_p, d2_p, d1_p, d11_b, d2
     *_v, d3_v, d4_v
        double precision d5_v, d6_v, d7_v, d8_v, d9_v, d10_v, d11_v, d9_
     *b, d2_b, d3_b
        double precision d4_b, d5_b, d6_b, d7_b, d8_b, g_costheta(g_pmax
     *_), g_r(ldg_r, nt, 3), g_temp(g_pmax_), g_d1_w(g_pmax_), g_bigr(g_
     *pmax_)
        double precision g_coschi(g_pmax_), g_xelih(g_pmax_), g_s1a(g_pm
     *ax_), g_s1c(g_pmax_), g_r1r3(g_pmax_), g_rcut2(g_pmax_), g_s1(g_pm
     *ax_), g_yhf(g_pmax_), g_behf(g_pmax_), g_xehf(g_pmax_)
        double precision g_s2(g_pmax_), g_dehfm(g_pmax_), g_xehfm(g_pmax
     *_), g_s2m(g_pmax_), g_xes3(g_pmax_), g_s3(g_pmax_), g_e(ldg_e, nt)
     *, g_hfdip_y(g_pmax_), g_hfdip(g_pmax_), g_hfpol(g_pmax_)
        double precision g_hf_did(g_pmax_), g_hf_disp(g_pmax_), g_e_lr_h
     *f(g_pmax_)
        integer g_ehfid
        data cevau /0.036749309d0/
C
C INTERACTION REGION DATA
C     LiH
        data relih /1.3432062d0/, bes1 /1.8058651d0/, s1c2 /12.937438d0/
        data s1c3 /11.0913d0/, bes1b /1.2111436d0/
        data rcut2rho /1.07143d0/, rcut2del /0.6d0/, r1r3al /0.4d0/
C
C     HF
        data rehf /1.733d0/, dehf /6.122d0/, rehfm /1.6739d0/, dehfc /0.
     *344574d0/, dehfcc /0.8660801d0/, s2rho /1.63d0/, s2alp /2.9941d0/,
     * behfm /0.74633431d0/, yhf_r /2.1042d0/, behf_0 /1.1622d0/, behf_1
     * /0.025647d0/, behf_2 /0.059062d0/
C
C     LiF
        data res3 /3.6d0/, bes3 /1.49061584d0/, des3 /0.17426687d0/
        data res3c /2.48571d0/, bes3c /1.92857d0/, des3c /0.25d0/
        data rcutrho /0.4d0/, rcutdel /0.89333d0/, r1r3al2 /0.7333d0/, b
     *es3z /0.d0/, bes3r /7.4799608d0/
C
C LONG RANGE FORCES DATA
C     HF dipole moment (angstroms vs au)
        data hfdip_a /1.13d0/, hfdip_re /0.924212d0/
        data hfdip_m /0.703661d0, 0.516815d0, 0.240628d0, -0.430194d0, -
     *2.213088d0, 5.535609d0, 12.872106d0, -42.060172d0, -17.398065d0, 8
     *4.207629d0, -41.961465d0/
C
C     HF physical properties (polarization in au, ie in eV)
        data hfpolperp /4.59d0/, hfpolpar /5.10d0/, hfie /16.044d0/
C
C     Li(2s) physical properties (polarization in au, ie in eV)
        data li2spol /165.d0/, li2sie /5.392d0/
C
C     misc long-range parameters (bohr)
        data rhf0did /6.0d0/, rhf0disp /6.0d0/, r2lrcut_a /2.0d0/, r2lrc
     *ut_r /3.d0/
        data cauang /0.52917706d0/, caudeb /2.54177d0/
        data mass /7.016003d0, 1.00783d0, 18.9984d0/
C
        save
C
C
        return
C
C ---------------------------------------------------------------------- 
C
        entry pot11(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do i = 1, nt
C
C PRECOMPUTE THINGS FOR USE LATER...
C     cosine of HF bond angle (i.e., Li-F-H angle)
          d2_v = r(i, 3) * r(i, 3)
          d3_p = 2.0d0 * r(i, 3)
          d4_v = r(i, 2) * r(i, 2)
          d2_p = 2.0d0 * r(i, 2)
          d7_v = r(i, 1) * r(i, 1)
          d1_p = 2.0d0 * r(i, 1)
          d9_v = (-2.d0) * r(i, 1)
          d10_v = d9_v * r(i, 2)
          d11_v = (d2_v - d4_v - d7_v) / d10_v
          d2_b = 1.0d0 / d10_v
          d3_b = (-d11_v) / d10_v
          d6_b = d3_b * r(i, 2) * (-2.d0) + (-d2_b) * d1_p
          d5_b = d3_b * d9_v + (-d2_b) * d2_p
          d11_b = d2_b * d3_p
          do g_i_ = 1, g_p_
            g_costheta(g_i_) = d6_b * g_r(g_i_, i, 1) + d5_b * g_r(g_i_,
     * i, 2) + d11_b * g_r(g_i_, i, 3)
          enddo
          costheta = d11_v
C--------
          d2_b = mass(3) / (mass(2) + mass(3))
          do g_i_ = 1, g_p_
            g_temp(g_i_) = d2_b * g_r(g_i_, i, 2)
          enddo
          temp = mass(3) / (mass(2) + mass(3)) * r(i, 2)
C--------
C     Li-FH Jacobi distance
          d2_v = r(i, 1) * r(i, 1)
          d2_p = 2.0d0 * r(i, 1)
          d4_v = temp * temp
          d1_p = 2.0d0 * temp
          d6_v = 2.d0 * r(i, 1)
          d7_v = d6_v * temp
          d7_b = (-costheta) * d6_v + d1_p
          d8_b = (-costheta) * temp * 2.d0 + d2_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-d7_v) * g_costheta(g_i_) + d7_b * g_temp(g_
     *i_) + d8_b * g_r(g_i_, i, 1)
          enddo
          d1_w = d2_v + d4_v - d7_v * costheta
          d2_v = sqrt(d1_w)

c          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
c          else
c             call ehufDO (9,d1_w, d2_v, d1_p,
c     +g_ehfid,
c     +185)
c          endif
          do g_i_ = 1, g_p_
            g_bigr(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          bigr = d2_v
C--------
C     cosine of HF Jacobi angle (i.e., Li--HF angle)
          d2_v = r(i, 1) * r(i, 1)
          d3_p = 2.0d0 * r(i, 1)
          d4_v = bigr * bigr
          d2_p = 2.0d0 * bigr
          d7_v = temp * temp
          d1_p = 2.0d0 * temp
          d10_v = (-2.d0) * bigr
          d11_v = d10_v * temp
          d12_v = (-(d2_v - d4_v - d7_v)) / d11_v
          d3_b = (-d12_v) / d11_v
          d7_b = -(1.0d0 / d11_v)
          d5_b = d3_b * d10_v + (-d7_b) * d1_p
          d6_b = d3_b * temp * (-2.d0) + (-d7_b) * d2_p
          d12_b = d7_b * d3_p
          do g_i_ = 1, g_p_
            g_coschi(g_i_) = d5_b * g_temp(g_i_) + d6_b * g_bigr(g_i_) +
     * d12_b * g_r(g_i_, i, 1)
          enddo
          coschi = d12_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = 0.0d0
          enddo
          d1_w = -1.d0
          d3_v = max (coschi, d1_w)

c          if (coschi .gt.  d1_w) then
             d1_p = 1.0d0
             d2_p = 0.0d0
c          else if (coschi .lt.  d1_w) then
c             d1_p = 0.0d0
c             d2_p = 1.0d0
c          else
c             call ehbfDO (7,coschi, d1_w, d3_v, d1_p, d2_p,
c     +g_ehfid,
c     +228)
c             d2_p = 1.0d0 -  d1_p
c          endif
          do g_i_ = 1, g_p_
            g_coschi(g_i_) = d2_p * g_d1_w(g_i_) + d1_p * g_coschi(g_i_)
          enddo
          coschi = d3_v
C--------
          d2_v = min (coschi, 1.d0)

c          if (coschi .lt.  1.d0) then
             d1_p = 1.0d0
             d2_p = 0.0d0
c          else if (coschi .gt.  1.d0) then
c             d1_p = 0.0d0
c             d2_p = 1.0d0
c          else
c             call ehbfDO (8,coschi, 1.d0, d2_v, d1_p, d2_p,
c     +g_ehfid,
c     +247)
c             d2_p = 1.0d0 -  d1_p
c          endif
          do g_i_ = 1, g_p_
            g_coschi(g_i_) = d1_p * g_coschi(g_i_)
          enddo
          coschi = d2_v
C--------
C     squared!
          r1s = r(i, 1) ** 2
          r2s = r(i, 2) ** 2
          r3s = r(i, 3) ** 2
C
C INTERACTION REGION
C     LiH
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bes1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-bes1) * (r(i, 1) - relih)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xelih(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xelih = d2_v
C--------
          do g_i_ = 1, g_p_
            g_s1a(g_i_) = s1c2 * g_xelih(g_i_)
          enddo
          s1a = s1c2 * xelih
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bes1b) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-bes1b) * (r(i, 1) - relih)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d3_b = s1c3 * d1_p
          do g_i_ = 1, g_p_
            g_s1c(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          s1c = s1c3 * d2_v
C--------
C
          do g_i_ = 1, g_p_
            g_r1r3(g_i_) = r1r3al * g_r(g_i_, i, 2) + (-g_r(g_i_, i, 3))
     * + g_r(g_i_, i, 1)
          enddo
          r1r3 = r(i, 1) - r(i, 3) + r1r3al * r(i, 2)
C--------
          d3_b = 1.0d0 / rcut2del
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r1r3(g_i_)
          enddo
          d1_w = (r1r3 - rcut2rho) / rcut2del
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_rcut2(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          rcut2 = 0.5d0 * (1.0d0 - d2_v)
C--------
C
          d3_v = s1c - s1a
          d2_b = 1.0d0 + (-rcut2)
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d3_v * g_rcut2(g_i_) + rcut2 * g_s1c(g_i_) + d2
     *_b * g_s1a(g_i_)
          enddo
          s1 = s1a + d3_v * rcut2
C--------
C
C     HF
          do g_i_ = 1, g_p_
            g_yhf(g_i_) = g_r(g_i_, i, 2)
          enddo
          yhf = r(i, 2) - yhf_r
C--------
          d4_v = yhf * yhf
          d1_p = 2.0d0 * yhf
          d5_b = behf_2 * d1_p + (-behf_1)
          do g_i_ = 1, g_p_
            g_behf(g_i_) = d5_b * g_yhf(g_i_)
          enddo
          behf = behf_0 - behf_1 * yhf + behf_2 * d4_v
C--------
          d4_v = r(i, 2) - rehf
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-behf) * g_r(g_i_, i, 2) + (-d4_v) * g_behf(
     *g_i_)
          enddo
          d1_w = (-behf) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xehf(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xehf = d2_v
C--------
          d2_v = dehf * xehf
          d3_v = xehf - 2.0d0
          d4_b = d2_v + d3_v * dehf
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d4_b * g_xehf(g_i_)
          enddo
          s2 = d2_v * d3_v
C--------
C
          d5_b = dehfc * 0.5d0
          do g_i_ = 1, g_p_
            g_dehfm(g_i_) = d5_b * g_coschi(g_i_)
          enddo
          dehfm = dehf - (dehfcc + dehfc * 0.5d0 * (1.d0 - coschi))
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-behfm) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-behfm) * (r(i, 2) - rehfm)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xehfm(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xehfm = d2_v
C--------
          d3_v = dehfm * xehfm
          d4_v = xehfm - 2.0d0
          d5_b = d4_v * xehfm
          d4_b = d3_v + d4_v * dehfm
          do g_i_ = 1, g_p_
            g_s2m(g_i_) = d4_b * g_xehfm(g_i_) + d5_b * g_dehfm(g_i_)
          enddo
          s2m = d3_v * d4_v
C--------
          d3_b = 1.0d0 / s2alp
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 3)
          enddo
          d1_w = (r(i, 3) - s2rho) / s2alp
          d4_v = (s2 - s2m) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.d0 + d6_v
          d7_b = d4_v * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s2(g_i_) + d2_b 
     ** g_s2m(g_i_)
          enddo
          s2 = s2m + d4_v * d7_v
C--------
C
C     LiF 
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bes3) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bes3) * (r(i, 3) - res3)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xes3(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xes3 = d2_v
C--------
          d2_v = des3 * xes3
          d3_v = xes3 - 2.0d0
          d4_b = d2_v + d3_v * des3
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d4_b * g_xes3(g_i_)
          enddo
          s3 = d2_v * d3_v
C--------
C
C     Sum of diatomics
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_s3(g_i_) + g_s2(g_i_) + g_s1(g_i_)
          enddo
          e(i) = s1 + s2 + s3 + dehf
C--------
C
C EVALUATE THE LONG-RANGE FORCES
C HF DIPOLE MOMENT
          d4_b = (-hfdip_a) * cauang
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (-hfdip_a) * (r(i, 2) * cauang - hfdip_re)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_hfdip_y(g_i_) = (-d1_p) * g_d1_w(g_i_)
          enddo
          hfdip_y = 1.d0 - d2_v
C--------
          do g_i_ = 1, g_p_
            g_hfdip(g_i_) = 0.0d0
          enddo
          hfdip = 0.d0
C--------
          do j = 1, 11
            d3_v = dble(j - 1)

c            if ( hfdip_y .ne. 0.0d0 ) then
               d4_v = hfdip_y ** ( d3_v - 2.0d0)
               d4_v =  d4_v * hfdip_y
               d1_p =  d3_v *  d4_v
               d4_v =  d4_v * hfdip_y
c            else
C              (hfdip_y = 0)
c               d4_v = hfdip_y **  d3_v
c
c               if (  d3_v .lt. 1.0d0 ) then
c                  call ehbfDO (10,hfdip_y, d3_v, d4_v, d1_p, 0.0d0,
c     +g_ehfid,
c     +463)
c               else if (  d3_v .lt. 2.0d0 ) then
c                  d1_p = 0.0d0
c                  call ehbfDO (10,hfdip_y, d3_v, d4_v, d1_p, 0.0d0,
c     +g_ehfid,
c     +468)
c               else
c                  d1_p = 0.0d0
c               endif
c            endif
            d6_b = 0.0d0
            d5_b = hfdip_m(j) * d1_p
            do g_i_ = 1, g_p_
              g_hfdip(g_i_) = d5_b * g_hfdip_y(g_i_) + g_hfdip(g_i_)
            enddo
            hfdip = hfdip + hfdip_m(j) * d4_v
C--------
          enddo
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = 20.d0 * g_r(g_i_, i, 2)
          enddo
          d1_w = 20.d0 * (r(i, 2) - 1.2d0)
          d2_v = hfdip * 0.5d0
          d4_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d4_v *  d4_v)
          d5_v = 1.d0 + d4_v
          d5_b = d2_v * d1_p
          d6_b = d5_v * 0.5d0
          do g_i_ = 1, g_p_
            g_hfdip(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_hfdip(g_i_)
          enddo
          hfdip = d2_v * d5_v
C--------
C
C HF POLARIZATION
          d2_v = coschi * coschi
          d2_p = 2.0d0 * coschi
          d5_v = coschi * coschi
          d1_p = 2.0d0 * coschi
          d5_b = hfpolpar * d1_p + (-hfpolperp) * d2_p
          do g_i_ = 1, g_p_
            g_hfpol(g_i_) = d5_b * g_coschi(g_i_)
          enddo
          hfpol = hfpolperp * (1.d0 - d2_v) + hfpolpar * d5_v
C--------
          hfpolavg = 2.d0 * hfpolperp / 3.d0 + hfpolpar / 3.d0
C
C HF + Li(2s) arrangement
C     DIPOLE-INDUCED-DIPOLE
          d2_v = coschi * coschi
          d2_p = 2.0d0 * coschi
          d6_v = (-0.5d0) * (3.d0 * d2_v + 1.d0) * li2spol
          d8_v = hfdip * hfdip
          d1_p = 2.0d0 * hfdip
          d4_b = d6_v * d1_p
          d9_b = d8_v * li2spol * (-0.5d0) * 3.d0 * d2_p
          do g_i_ = 1, g_p_
            g_hf_did(g_i_) = d4_b * g_hfdip(g_i_) + d9_b * g_coschi(g_i_
     *)
          enddo
          hf_did = d6_v * d8_v
C--------
          d3_v = bigr ** ( 7 - 2)
          d3_v =  d3_v * bigr
          d1_p =  7 *  d3_v
          d3_v =  d3_v * bigr
          d4_v = d3_v + rhf0did ** 7
          d5_v = hf_did / d4_v
          d6_v = d5_v / cevau
          d4_b = bigr * (1.0d0 / cevau)
          d5_b = d4_b * (1.0d0 / d4_v)
          d3_b = d6_v + d4_b * ((-d5_v) / d4_v) * d1_p
          do g_i_ = 1, g_p_
            g_hf_did(g_i_) = d3_b * g_bigr(g_i_) + d5_b * g_hf_did(g_i_)
          enddo
          hf_did = d6_v * bigr
C--------
C     cut off angular dependence in interaction region
          d3_v = (bigr / rhf0did) ** ( 4 - 2)
          d3_v =  d3_v * (bigr / rhf0did)
          d1_p =  4 *  d3_v
          d3_v =  d3_v * (bigr / rhf0did)
          d4_b = (-d1_p) * (1.0d0 / rhf0did)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_bigr(g_i_)
          enddo
          d1_w = -d3_v
          d3_v = coschi * coschi
          d2_p = 2.0d0 * coschi
          d5_v = 3.d0 * d3_v + 1.d0
          d6_v = 2.5d0 / d5_v
          d7_v = d6_v - 1.d0
          d9_v = exp(d1_w)
          d1_p =  d9_v
          d11_v = 1.d0 + d7_v * d9_v
          d7_b = hf_did * d7_v * d1_p
          d12_b = hf_did * d9_v * ((-d6_v) / d5_v) * 3.d0 * d2_p
          do g_i_ = 1, g_p_
            g_hf_did(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_coschi(g_i_
     *) + d11_v * g_hf_did(g_i_)
          enddo
          hf_did = hf_did * d11_v
C--------
C
C     DISPERSION
          d3_b = 1.0d0 / (li2sie + hfie) * ((-1.5d0) * li2sie * hfie * l
     *i2spol)
          do g_i_ = 1, g_p_
            g_hf_disp(g_i_) = d3_b * g_hfpol(g_i_)
          enddo
          hf_disp = (-1.5d0) * li2sie * hfie * li2spol * hfpol / (li2sie
     * + hfie)
C--------
          d3_v = bigr ** ( 7 - 2)
          d3_v =  d3_v * bigr
          d1_p =  7 *  d3_v
          d3_v =  d3_v * bigr
          d4_v = d3_v + rhf0disp ** 7
          d5_v = hf_disp / d4_v
          d4_b = bigr * (1.0d0 / d4_v)
          d3_b = d5_v + bigr * ((-d5_v) / d4_v) * d1_p
          do g_i_ = 1, g_p_
            g_hf_disp(g_i_) = d3_b * g_bigr(g_i_) + d4_b * g_hf_disp(g_i
     *_)
          enddo
          hf_disp = d5_v * bigr
C--------
C     cut off angular dependence in interaction region
          d3_v = (bigr / rhf0disp) ** ( 4 - 2)
          d3_v =  d3_v * (bigr / rhf0disp)
          d1_p =  4 *  d3_v
          d3_v =  d3_v * (bigr / rhf0disp)
          d4_b = (-d1_p) * (1.0d0 / rhf0disp)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_bigr(g_i_)
          enddo
          d1_w = -d3_v
          d3_v = hfpolavg / hfpol
          d4_v = d3_v - 1.d0
          d6_v = exp(d1_w)
          d1_p =  d6_v
          d8_v = 1.d0 + d4_v * d6_v
          d7_b = hf_disp * d4_v * d1_p
          d9_b = hf_disp * d6_v * ((-d3_v) / hfpol)
          do g_i_ = 1, g_p_
            g_hf_disp(g_i_) = d7_b * g_d1_w(g_i_) + d9_b * g_hfpol(g_i_)
     * + d8_v * g_hf_disp(g_i_)
          enddo
          hf_disp = hf_disp * d8_v
C--------
C
C     SUM FORCES
          do g_i_ = 1, g_p_
            g_e_lr_hf(g_i_) = g_hf_did(g_i_) + g_hf_disp(g_i_)
          enddo
          e_lr_hf = hf_disp + hf_did
C--------
C     cut off LR forces for large r[hf] i.e., r(i,2)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = r2lrcut_a * g_r(g_i_, i, 2)
          enddo
          d1_w = r2lrcut_a * (r(i, 2) - r2lrcut_r)
          d2_v = e_lr_hf * 0.5d0
          d4_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d4_v *  d4_v)
          d5_v = 1.d0 - d4_v
          d5_b = (-d2_v) * d1_p
          d6_b = d5_v * 0.5d0
          do g_i_ = 1, g_p_
            g_e_lr_hf(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_e_lr_hf(g_i
     *_)
          enddo
          e_lr_hf = d2_v * d5_v
C--------
C
C LiF + H
C     NO LONG-RANGE FORCES
C
C LiH + F
C     NO LONG-RANGE FORCES
C
C COMBINE INTERACTION REGION ENERGY AND LONG-RANGE FORCES
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_e_lr_hf(g_i_) + g_e(g_i_, i)
          enddo
          e(i) = e(i) + e_lr_hf
C--------
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_e(g_i_, i)
          enddo
          e(i) = e(i) * cevau
C--------
C
        enddo
        return
C
      end
C

C**************************************************************************
C**************************************************************************

C                           DISCLAIMER
C
C   This file was generated on 07/27/01 by the version of
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
C                      U12 COUPLING SURFACE
C
      subroutine prepot12
C
        implicit none
C
C MISC VARIABLES
        integer i, nt
        double precision r(nt, 3), e(nt)
C
C LiH VARIABLES
        integer ng1
        double precision g1, g1r, g1a, h1, h1del, h1rho
C
C HF VARIABLES
        double precision j2p, j2, j2pdel, j2rho, j2prho, j2del
C
C LiF VARIABLES
        integer ng3
        double precision g3, g3r, g3a, h3, h3del, h3rho
C
C CUTOFF VARIABLES
        double precision gcrho, cevau, gcdel
        double precision cut2, gc
C
C MISC DATA
        integer g_pmax_, g_p_
        parameter (g_pmax_ = 3)
        integer g_i_, ldg_r, ldg_e
        double precision d13_b, d12_b, d11_v, d2_p, d9_v, d8_b, d1_w, d7
     *_b, d2_v, d3_v
        double precision d4_v, d1_p, d2_b, d3_b, d4_b, d5_v, d6_v, g_d1_
     *w(g_pmax_), g_r(ldg_r, nt, 3), g_g1(g_pmax_)
        double precision g_g3(g_pmax_), g_j2(g_pmax_), g_j2p(g_pmax_), g
     *_h1(g_pmax_), g_h3(g_pmax_), g_e(ldg_e, nt), g_gc(g_pmax_), g_cut2
     *(g_pmax_)
        integer g_ehfid
        data cevau /0.036749309d0/
C
C LiH DATA
        data g1a /1.27742d0/, g1r /3.0048876d0/, ng1 /6/
        data h1rho /4.9843597d0/, h1del /2.0899315d0/
C
C LiF DATA
        data g3a /0.48d0/, g3r /3.4975562d0/, ng3 /8/
        data h3rho /2.000978d0/, h3del /0.90615835d0/
C
C HF DATA
        data j2rho /1.15484d0/, j2del /1.75806d0/
        data j2prho /1.45161d0/, j2pdel /0.98387d0/
C
C CUTOFF DATA
        data gcrho /4.87097d0/, gcdel /2.d0/
C
        save
C
c        call ehsfid(g_ehfid, 'prepot12','g_u12.f')
C
        return
C
C ---------------------------------------------------------------------- 
C
        entry pot12(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do i = 1, nt
C
          d4_b = (-dble(ng1)) * (1.0d0 / g1r)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_r(g_i_, i, 1)
          enddo
          d1_w = (-dble(ng1)) * (r(i, 1) / g1r - 1.0d0)

          if (  ng1 .ne. 0 ) then
             d3_v = (r(i, 1) / g1r) ** ( ng1-1)
             d2_p =  ng1 *  d3_v
             d3_v =  d3_v * (r(i, 1) / g1r)
          else
C            Maybe this should be  d3_v = (r(i, 1) / g1r) **  ng1
             d3_v = 1.0d0
             d2_p = 0.0d0
          endif
          d4_v = g1a * d3_v
          d6_v = exp(d1_w)
          d1_p =  d6_v
          d4_b = d4_v * d1_p
          d7_b = d6_v * g1a * d2_p * (1.0d0 / g1r)
          do g_i_ = 1, g_p_
            g_g1(g_i_) = d4_b * g_d1_w(g_i_) + d7_b * g_r(g_i_, i, 1)
          enddo
          g1 = d4_v * d6_v
C--------
C
          d4_b = (-dble(ng3)) * (1.0d0 / g3r)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_r(g_i_, i, 3)
          enddo
          d1_w = (-dble(ng3)) * (r(i, 3) / g3r - 1.0d0)

          if (  ng3 .ne. 0 ) then
             d3_v = (r(i, 3) / g3r) ** ( ng3-1)
             d2_p =  ng3 *  d3_v
             d3_v =  d3_v * (r(i, 3) / g3r)
          else
C            Maybe this should be  d3_v = (r(i, 3) / g3r) **  ng3
             d3_v = 1.0d0
             d2_p = 0.0d0
          endif
          d4_v = g3a * d3_v
          d6_v = exp(d1_w)
          d1_p =  d6_v
          d4_b = d4_v * d1_p
          d7_b = d6_v * g3a * d2_p * (1.0d0 / g3r)
          do g_i_ = 1, g_p_
            g_g3(g_i_) = d4_b * g_d1_w(g_i_) + d7_b * g_r(g_i_, i, 3)
          enddo
          g3 = d4_v * d6_v
C--------
C
          d3_b = 1.0d0 / j2del
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (r(i, 2) - j2rho) / j2del
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = 0.5d0 * d1_p
          do g_i_ = 1, g_p_
            g_j2(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          j2 = 0.5d0 * (1.0d0 + d2_v)
C--------
          d3_b = 1.0d0 / j2pdel
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (r(i, 2) - j2prho) / j2pdel
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = 0.5d0 * d1_p
          do g_i_ = 1, g_p_
            g_j2p(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          j2p = 0.5d0 * (1.0d0 + d2_v)
C--------
C
          d3_b = 1.0d0 / h1del
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 1)
          enddo
          d1_w = (r(i, 1) - h1rho) / h1del
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_h1(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          h1 = 0.5d0 * (1.0d0 - d2_v)
C--------
          d3_b = 1.0d0 / h3del
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 3)
          enddo
          d1_w = (r(i, 3) - h3rho) / h3del
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_h3(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          h3 = 0.5d0 * (1.0d0 - d2_v)
C--------
C
          d3_v = g1 * j2
          d5_v = 1.0d0 - h3
          d9_v = g3 * j2p
          d11_v = 1.0d0 - h1
          d7_b = d11_v * j2p
          d8_b = d11_v * g3
          d12_b = d5_v * j2
          d13_b = d5_v * g1
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = (-d9_v) * g_h1(g_i_) + d8_b * g_j2p(g_i_) + d
     *7_b * g_g3(g_i_) + (-d3_v) * g_h3(g_i_) + d13_b * g_j2(g_i_) + d12
     *_b * g_g1(g_i_)
          enddo
          e(i) = d3_v * d5_v + d9_v * d11_v
C--------
C
          do g_i_ = 1, g_p_
            g_gc(g_i_) = g_r(g_i_, i, 2)
          enddo
          gc = r(i, 2) - gcrho
C--------
          d2_b = 1.0d0 / gcdel
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_b * g_gc(g_i_)
          enddo
          d1_w = gc / gcdel
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_cut2(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          cut2 = 0.5d0 - 0.5d0 * d2_v
C--------
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = e(i) * g_cut2(g_i_) + cut2 * g_e(g_i_, i)
          enddo
          e(i) = e(i) * cut2
C--------
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_e(g_i_, i)
          enddo
          e(i) = e(i) * cevau
C--------
C
        enddo
        return
      end
C

C**************************************************************************
C**************************************************************************

C                           DISCLAIMER
C
C   This file was generated on 07/27/01 by the version of
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
C                            U22 SURFACE
C
      subroutine prepot22
        implicit none
C
C MISC VARIABLES
        integer i, j, k, nt
        double precision r(nt, 3), e(nt), r1s, r2s, r3s, cs, cevau, uli2
     *p
        double precision c2a, c2b, c2c, w, cplg2
C
C INTERACTION REGION VARIABLES 
C     LiH
        double precision r1r3, xelih, s1c, bes1, rcut2del, rcut2rho, s1c
     *3
        double precision s1c1, belihc, bet1, relih, s1c2, t1c1, t1c2
        double precision r1r3al, rcut2
        double precision coul1, exch1, s1, t1, reliht
C
C     HF
        double precision behf, xehfa, dehfm, yhf, xehfa_180, t2a3, t2a2,
     * bet2
        double precision t2a1, t2a4, dehf, rehf
        double precision s2alp, rehft, bet2_180, s2rho, dehfc
        double precision coul2, exch2, t2, s2, dehfcc, s2cut, xehfm, s2m
     *, s2del, s2cc, s2cc_d, s2cc_r, s2cc_b
C
C     LiF
        integer nb
        double precision bet3, bet3_180, t3a2, t3a3, t3a4, belihc1
        double precision belifff, gam, belifc2, delifc1, delifc2, belifm
        double precision s3c, xelihf_180, xelif, s3a, delifm, gc2, cut2
        double precision del2_0, del2_180, rho2_0, rho2_180, theta2_0, g
     *rho
        double precision theta2_180, beliff, belifmod, belif, xelif_180,
     * relif
        double precision belifc1, t3a1, gswt, delif, rho2, del2, theta2,
     * gdel
        double precision coul3, exch3, s3, t3, s3cut, swit, swit_r, swit
     *_a, swit_d
C
C LONG RANGE FORCES VARIABLES
C     HF physical properties
        double precision hfdip_y, hfdip_a, hfdip_re, hfdip_m(11), hfdip,
     * hfpol, hfpolperp, hfpolpar, hfie, hfpolavg
        double precision hfqua_y, hfqua_a, hfqua_re, hfqua_m(10), hfqua
C
C     Li(2p) physical properties
        double precision li2ppol, li2pie, li2pqua, li2ppolpar, li2ppolpe
     *rp, li2ppolavg
C
C     LiF physical properties
        double precision lifpol, lifie
        double precision lifdip_y, lifdip_re, lifdip_1, lifdip_2, lifdip
     *_3, lifdip, lifdip_4, lifdip_beta, lifdip_d
C
C     H physical properties
        double precision hie, hpol
C
C     misc long-range variables
        double precision rhf0did, rhf0disp, hf_disp, e_lr_hf, e_lr_lif, 
     *rlif0did, rlif0disp, lif_disp, lif_did, hf_qq, hf_dq, rhf0qq, rhf0
     *dq, hf_did
        double precision mass(3), cauang, caudeb, costhetahf, costhetali
     *f, bigrhf, bigrlif, coschihf, coschilif, temp, r2lrcut_a, r2lrcut_
     *r, r3lrcut_a, r3lrcut_r
C
C MISC DATA
        integer g_pmax_,g_p_
        parameter (g_pmax_ = 3)
        integer g_i_, ldg_r, ldg_e
        double precision d1_w, d17_b, d12_v, d12_b, d14_b, d13_v, d16_v,
     * d2_v, d3_v, d4_v
        double precision d5_v, d6_v, d7_v, d8_v, d9_v, d10_v, d11_v, d15
     *_v, d2_b, d3_b
        double precision d4_b, d5_b, d6_b, d7_b, d8_b, d9_b, d10_b, d11_
     *b, d1_p, d2_p
        double precision d3_p, g_costhetahf(g_pmax_), g_r(ldg_r, nt, 3),
     * g_temp(g_pmax_), g_d1_w(g_pmax_), g_bigrhf(g_pmax_), g_coschihf(g
     *_pmax_), g_costhetalif(g_pmax_), g_bigrlif(g_pmax_), g_coschilif(g
     *_pmax_)
        double precision g_r1s(g_pmax_), g_r2s(g_pmax_), g_r3s(g_pmax_),
     * g_cs(g_pmax_), g_xelih(g_pmax_), g_s1(g_pmax_), g_s1c(g_pmax_), g
     *_r1r3(g_pmax_), g_rcut2(g_pmax_), g_t1(g_pmax_)
        double precision g_yhf(g_pmax_), g_behf(g_pmax_), g_xehfa(g_pmax
     *_), g_s2(g_pmax_), g_s2cut(g_pmax_), g_s2cc(g_pmax_), g_dehfm(g_pm
     *ax_), g_xehfm(g_pmax_), g_s2m(g_pmax_), g_xehfa_180(g_pmax_)
        double precision g_t2(g_pmax_), g_belif(g_pmax_), g_gswt(g_pmax_
     *), g_beliff(g_pmax_), g_belifmod(g_pmax_), g_xelif(g_pmax_), g_s3a
     *(g_pmax_), g_s3cut(g_pmax_), g_swit(g_pmax_), g_delifm(g_pmax_)
        double precision g_belifm(g_pmax_), g_s3c(g_pmax_), g_rho2(g_pma
     *x_), g_theta2(g_pmax_), g_del2(g_pmax_), g_gc2(g_pmax_), g_cut2(g_
     *pmax_), g_s3(g_pmax_), g_xelif_180(g_pmax_), g_t3(g_pmax_)
        double precision g_coul1(g_pmax_), g_coul2(g_pmax_), g_coul3(g_p
     *max_), g_exch1(g_pmax_), g_exch2(g_pmax_), g_exch3(g_pmax_), g_w(g
     *_pmax_), g_cplg2(g_pmax_), g_e(ldg_e, nt), g_hfdip_y(g_pmax_)
        double precision g_hfdip(g_pmax_), g_hfpol(g_pmax_), g_hfqua_y(g
     *_pmax_), g_hfqua(g_pmax_), g_lifdip_y(g_pmax_), g_lifdip_beta(g_pm
     *ax_), g_lifdip(g_pmax_), g_li2ppol(g_pmax_), g_hf_did(g_pmax_), g_
     *hf_disp(g_pmax_)
        double precision g_hf_qq(g_pmax_), g_hf_dq(g_pmax_), g_e_lr_hf(g
     *_pmax_), g_lif_did(g_pmax_), g_lif_disp(g_pmax_), g_e_lr_lif(g_pma
     *x_)
        integer g_ehfid
        data c2a /3.5d0/, c2b /0.27362d0/, c2c /0.15d0/
        data uli2p /1.848d0/, cevau /0.036749309d0/
C
C INTERACTION REGION DATA
C     LiH
        data s1c1 /4.32258d0/, s1c2 /7.06452d0/, t1c1 /1.1935483d0/
        data t1c2 /13.548387d0/, belihc /1.36d0/, bet1 /2.41319648d0/
        data relih /1.2d0/, r1r3al /1.0d0/, rcut2rho /0.72d0/
        data rcut2del /0.5d0/, s1c3 /14.74194d0/, bes1 /0.90667d0/
        data reliht /1.203323d0/
C
C     HF
        data rehf /1.733d0/, dehf /6.122d0/
        data t2a1 /1.4242424d0/, t2a2 /14.203323d0/, t2a3 /6.15835777d0/
        data t2a4 /4.56304985d0/, bet2 /2.4105572d0/, bet2_180 /1.404692
     *1d0/
        data rehft /1.61329423d0/, s2rho /1.00293d0/
        data s2alp /4.0899d0/, dehfc /0.0486803519d0/, dehfcc /0.7262952
     *1d0/
        data s2cc_b /1.799609d0/, s2cc_r /1.60215d0/, s2cc_d /2.3841642d
     *0/
C
C     LiF
        data belifc1 /0.819648d0/, belifc2 /0.7795699d0/, delifc1 /2.352
     *884d0/
        data delifc2 /5.5904203d0/, belifff /0.13225806d0/, nb /8/
        data gam /4.6304985d0/, grho /3.382209d0/, gdel /0.49472140d0/
        data t3a1 /-.0552298d0/, t3a2 /1.8729228d0/, t3a3 /0.94662756d0/
        data t3a4 /0.2994134d0/, bet3 /0.9519061d0/, bet3_180 /0.4531769
     *d0/
        data rho2_0 /0.2017595d0/, rho2_180 /2.d0/, del2_0 /0.51710655d0
     */
        data del2_180 /0.43695014d0/, theta2_0 /0.12463343d0/
        data theta2_180 /0.14907135d0/, delif /5.947d0/, relif /2.9553d0
     */
        data swit_r /3.340762d0/, swit_a /0.56353861d0/, swit_d /1.02439
     *98d0/
C
C LONG RANGE FORCES DATA
C     HF dipole moment
        data hfdip_a /1.13d0/, hfdip_re /0.924212d0/
        data hfdip_m /0.703661d0, 0.516815d0, 0.240628d0, -0.430194d0, -
     *2.213088d0, 5.535609d0, 12.872106d0, -42.060172d0, -17.398065d0, 8
     *4.207629d0, -41.961465d0/
C
C     HF quadrupole moment
        data hfqua_a /1.37d0/, hfqua_re /0.9168d0/
        data hfqua_m /1.6527d0, 2.072d0, 2.0521d0, 3.545d0, -9.8031d0, -
     *19.3808d0, 76.4423d0, 0.d0, -161.5896d0, 105.7460d0/
C
C     HF pol (au) and ie (eV)
        data hfpolperp /4.59d0/, hfpolpar /5.10d0/, hfie /16.044d0/
C
C     Li(2p) pol (au) and ie (eV) and quadrupole moment (au)
        data li2ppolpar /131.d0/, li2ppolperp /129.d0/, li2pie /3.544d0/
     *, li2pqua /11.1d0/
C
C     LiF dipole moment (distance and moment in au)
        data lifdip_re /10.4994d0/, lifdip_1 /0.02435d0/, lifdip_2 /9.93
     *55d-8/, lifdip_3 /0.0015999d0/, lifdip_4 /4.471959d0/, lifdip_d /9
     *.30039d0/
C
C     LiF pol (au) and ie (eV)
        data lifpol /72.9d0/, lifie /11.3d0/
C
C     H pol (au) and ie (eV)
        data hie /13.598d0/, hpol /4.4997d0/
C
C     misc long-range
        data rhf0did /6.d0/, rhf0disp /6.d0/, rhf0qq /6.d0/, rhf0dq /6.d
     *0/
C      data rlif0did/7.5d0/,rlif0disp/7.5d0/,r2LRcut_a/2.0d0/,
        data rlif0did /7.d0/, rlif0disp /7.d0/, r2lrcut_a /2.0d0/, r2lrc
     *ut_r /3.d0/, r3lrcut_a /2.d0/, r3lrcut_r /4.5d0/
        data cauang /0.52917706d0/, caudeb /2.54177d0/
        data mass /7.016003d0, 1.00783d0, 18.9984d0/
C
        save
c        data g_ehfid /0/
C
c        call ehsfid(g_ehfid, 'prepot22','g_u22.f')
C
        return
C
C ---------------------------------------------------------------------- 
C
        entry pot22(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do 99999 i = 1, nt
C
C PRECOMPUTE THINGS FOR USE LATER
C     cosine of HF bond angle (i.e., Li-F-H angle)
          d2_v = r(i, 3) * r(i, 3)
          d3_p = 2.0d0 * r(i, 3)
          d4_v = r(i, 2) * r(i, 2)
          d2_p = 2.0d0 * r(i, 2)
          d7_v = r(i, 1) * r(i, 1)
          d1_p = 2.0d0 * r(i, 1)
          d9_v = (-2.d0) * r(i, 1)
          d10_v = d9_v * r(i, 2)
          d11_v = (d2_v - d4_v - d7_v) / d10_v
          d2_b = 1.0d0 / d10_v
          d3_b = (-d11_v) / d10_v
          d6_b = d3_b * r(i, 2) * (-2.d0) + (-d2_b) * d1_p
          d5_b = d3_b * d9_v + (-d2_b) * d2_p
          d11_b = d2_b * d3_p
          do g_i_ = 1, g_p_
            g_costhetahf(g_i_) = d6_b * g_r(g_i_, i, 1) + d5_b * g_r(g_i
     *_, i, 2) + d11_b * g_r(g_i_, i, 3)
          enddo
          costhetahf = d11_v
C--------
          d2_b = mass(3) / (mass(2) + mass(3))
          do g_i_ = 1, g_p_
            g_temp(g_i_) = d2_b * g_r(g_i_, i, 2)
          enddo
          temp = mass(3) / (mass(2) + mass(3)) * r(i, 2)
C--------
C     Li-FH Jacobi distance
          d2_v = r(i, 1) * r(i, 1)
          d2_p = 2.0d0 * r(i, 1)
          d4_v = temp * temp
          d1_p = 2.0d0 * temp
          d6_v = 2.d0 * r(i, 1)
          d7_v = d6_v * temp
          d7_b = (-costhetahf) * d6_v + d1_p
          d8_b = (-costhetahf) * temp * 2.d0 + d2_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-d7_v) * g_costhetahf(g_i_) + d7_b * g_temp(
     *g_i_) + d8_b * g_r(g_i_, i, 1)
          enddo
          d1_w = d2_v + d4_v - d7_v * costhetahf
          d2_v = sqrt(d1_w)

c          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
c          else
c             call ehufDO (9,d1_w, d2_v, d1_p,
c     +g_ehfid,
c     +278)
c          endif
          do g_i_ = 1, g_p_
            g_bigrhf(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          bigrhf = d2_v
C--------
C     cosine of HF Jacobi angle (i.e., Li--HF angle)
          d2_v = r(i, 1) * r(i, 1)
          d3_p = 2.0d0 * r(i, 1)
          d4_v = bigrhf * bigrhf
          d2_p = 2.0d0 * bigrhf
          d7_v = temp * temp
          d1_p = 2.0d0 * temp
          d10_v = (-2.d0) * bigrhf
          d11_v = d10_v * temp
          d12_v = (-(d2_v - d4_v - d7_v)) / d11_v
          d3_b = (-d12_v) / d11_v
          d7_b = -(1.0d0 / d11_v)
          d5_b = d3_b * d10_v + (-d7_b) * d1_p
          d6_b = d3_b * temp * (-2.d0) + (-d7_b) * d2_p
          d12_b = d7_b * d3_p
          do g_i_ = 1, g_p_
            g_coschihf(g_i_) = d5_b * g_temp(g_i_) + d6_b * g_bigrhf(g_i
     *_) + d12_b * g_r(g_i_, i, 1)
          enddo
          coschihf = d12_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = 0.0d0
          enddo
          d1_w = -1.d0
          d3_v = max (d1_w, coschihf)

c          if (d1_w .gt.  coschihf) then
c             d1_p = 1.0d0
c             d2_p = 0.0d0
c          else if (d1_w .lt.  coschihf) then
             d1_p = 0.0d0
             d2_p = 1.0d0
c          else
c             call ehbfDO (7,d1_w, coschihf, d3_v, d1_p, d2_p,
c     +g_ehfid,
c     +321)
c             d2_p = 1.0d0 -  d1_p
c          endif
          do g_i_ = 1, g_p_
            g_coschihf(g_i_) = d2_p * g_coschihf(g_i_) + d1_p * g_d1_w(g
     *_i_)
          enddo
          coschihf = d3_v
C--------
          d2_v = min (1.d0, coschihf)

c          if (1.d0 .lt.  coschihf) then
c             d1_p = 1.0d0
c             d2_p = 0.0d0
c          else if (1.d0 .gt.  coschihf) then
             d1_p = 0.0d0
             d2_p = 1.0d0
c          else
c             call ehbfDO (8,1.d0, coschihf, d2_v, d1_p, d2_p,
c     +g_ehfid,
c     +341)
c             d2_p = 1.0d0 -  d1_p
c          endif
          do g_i_ = 1, g_p_
            g_coschihf(g_i_) = d2_p * g_coschihf(g_i_)
          enddo
          coschihf = d2_v
C--------
C
C     cosine of LiF bond angle (i.e., Li-F-H angle)
          d2_v = r(i, 1) * r(i, 1)
          d3_p = 2.0d0 * r(i, 1)
          d4_v = r(i, 2) * r(i, 2)
          d2_p = 2.0d0 * r(i, 2)
          d7_v = r(i, 3) * r(i, 3)
          d1_p = 2.0d0 * r(i, 3)
          d9_v = (-2.d0) * r(i, 3)
          d10_v = d9_v * r(i, 2)
          d11_v = (d2_v - d4_v - d7_v) / d10_v
          d2_b = 1.0d0 / d10_v
          d3_b = (-d11_v) / d10_v
          d6_b = d3_b * r(i, 2) * (-2.d0) + (-d2_b) * d1_p
          d5_b = d3_b * d9_v + (-d2_b) * d2_p
          d11_b = d2_b * d3_p
          do g_i_ = 1, g_p_
            g_costhetalif(g_i_) = d6_b * g_r(g_i_, i, 3) + d5_b * g_r(g_
     *i_, i, 2) + d11_b * g_r(g_i_, i, 1)
          enddo
          costhetalif = d11_v
C--------
          d2_b = mass(1) / (mass(1) + mass(3))
          do g_i_ = 1, g_p_
            g_temp(g_i_) = d2_b * g_r(g_i_, i, 3)
          enddo
          temp = mass(1) / (mass(1) + mass(3)) * r(i, 3)
C--------
C     H-LiF Jacobi distance
          d2_v = r(i, 2) * r(i, 2)
          d2_p = 2.0d0 * r(i, 2)
          d4_v = temp * temp
          d1_p = 2.0d0 * temp
          d6_v = 2.d0 * r(i, 2)
          d7_v = d6_v * temp
          d7_b = (-costhetalif) * d6_v + d1_p
          d8_b = (-costhetalif) * temp * 2.d0 + d2_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-d7_v) * g_costhetalif(g_i_) + d7_b * g_temp
     *(g_i_) + d8_b * g_r(g_i_, i, 2)
          enddo
          d1_w = d2_v + d4_v - d7_v * costhetalif
          d2_v = sqrt(d1_w)

c          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
c          else
c             call ehufDO (9,d1_w, d2_v, d1_p,
c     +g_ehfid,
c     +398)
c          endif
          do g_i_ = 1, g_p_
            g_bigrlif(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          bigrlif = d2_v
C--------
C     cosine of LiF Jacobi angle (i.e., H--LiF angle)
          d2_v = r(i, 2) * r(i, 2)
          d3_p = 2.0d0 * r(i, 2)
          d4_v = bigrlif * bigrlif
             d2_p = 2.0d0 * bigrlif
          d7_v = temp * temp
          d1_p = 2.0d0 * temp
          d10_v = (-2.d0) * bigrlif
          d11_v = d10_v * temp
          d12_v = (-(d2_v - d4_v - d7_v)) / d11_v
          d3_b = (-d12_v) / d11_v
          d7_b = -(1.0d0 / d11_v)
          d5_b = d3_b * d10_v + (-d7_b) * d1_p
          d6_b = d3_b * temp * (-2.d0) + (-d7_b) * d2_p
          d12_b = d7_b * d3_p
          do g_i_ = 1, g_p_
            g_coschilif(g_i_) = d5_b * g_temp(g_i_) + d6_b * g_bigrlif(g
     *_i_) + d12_b * g_r(g_i_, i, 2)
          enddo
          coschilif = d12_v
C--------
          d2_v = min (1.d0, coschilif)

c             if (1.d0 .lt.  coschilif) then
c                d1_p = 1.0d0
c                d2_p = 0.0d0
c             else if (1.d0 .gt.  coschilif) then
                d1_p = 0.0d0
                d2_p = 1.0d0
c             else
c                call ehbfDO (8,1.d0, coschilif, d2_v, d1_p, d2_p,
c     +g_ehfid,
c     +437)
c                   d2_p = 1.0d0 -  d1_p
c                endif
          do g_i_ = 1, g_p_
            g_coschilif(g_i_) = d2_p * g_coschilif(g_i_)
          enddo
          coschilif = d2_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = 0.0d0
          enddo
          d1_w = -1.d0
          d3_v = max (d1_w, coschilif)

c             if (d1_w .gt.  coschilif) then
c                d1_p = 1.0d0
c                d2_p = 0.0d0
c             else if (d1_w .lt.  coschilif) then
                d1_p = 0.0d0
                d2_p = 1.0d0
c             else
c                call ehbfDO (7,d1_w, coschilif, d3_v, d1_p, d2_p,
c     +g_ehfid,
c     +460)
c                   d2_p = 1.0d0 -  d1_p
c                endif
          do g_i_ = 1, g_p_
            g_coschilif(g_i_) = d2_p * g_coschilif(g_i_) + d1_p * g_d1_w
     *(g_i_)
          enddo
          coschilif = d3_v
C--------
C     squared!
          d2_v = r(i, 1) * r(i, 1)
          d1_p = 2.0d0 * r(i, 1)
          do g_i_ = 1, g_p_
            g_r1s(g_i_) = d1_p * g_r(g_i_, i, 1)
          enddo
          r1s = d2_v
C--------
          d2_v = r(i, 2) * r(i, 2)
          d1_p = 2.0d0 * r(i, 2)
          do g_i_ = 1, g_p_
            g_r2s(g_i_) = d1_p * g_r(g_i_, i, 2)
          enddo
          r2s = d2_v
C--------
          d2_v = r(i, 3) * r(i, 3)
          d1_p = 2.0d0 * r(i, 3)
          do g_i_ = 1, g_p_
            g_r3s(g_i_) = d1_p * g_r(g_i_, i, 3)
          enddo
          r3s = d2_v
C--------
C     Li-F-H bondangle
          do g_i_ = 1, g_p_
            g_cs(g_i_) = -g_r1s(g_i_) + g_r2s(g_i_) + g_r3s(g_i_)
          enddo
          cs = r3s + r2s - r1s
C--------
          d3_v = 2.0d0 * r(i, 2)
          d5_v = d3_v * r(i, 3)
          d6_v = cs / d5_v
          d2_b = 1.0d0 / d5_v
          d3_b = (-d6_v) / d5_v
          d5_b = d3_b * d3_v
          d6_b = d3_b * r(i, 3) * 2.0d0
          do g_i_ = 1, g_p_
            g_cs(g_i_) = d5_b * g_r(g_i_, i, 3) + d6_b * g_r(g_i_, i, 2)
     * + d2_b * g_cs(g_i_)
          enddo
          cs = d6_v
C--------
          d2_v = min (cs, 1.0d0)

c          if (cs .lt.  1.0d0) then
             d1_p = 1.0d0
             d2_p = 0.0d0
c          else if (cs .gt.  1.0d0) then
c             d1_p = 0.0d0
c             d2_p = 1.0d0
c          else
c             call ehbfDO (8,cs, 1.0d0, d2_v, d1_p, d2_p,
c     +g_ehfid,
c     +521)
c             d2_p = 1.0d0 -  d1_p
c          endif
          do g_i_ = 1, g_p_
            g_cs(g_i_) = d1_p * g_cs(g_i_)
          enddo
          cs = d2_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = 0.0d0
          enddo
          d1_w = -1.0d0
          d3_v = max (cs, d1_w)

c          if (cs .gt.  d1_w) then
             d1_p = 1.0d0
             d2_p = 0.0d0
c          else if (cs .lt.  d1_w) then
c             d1_p = 0.0d0
c             d2_p = 1.0d0
c          else
c             call ehbfDO (7,cs, d1_w, d3_v, d1_p, d2_p,
c     +g_ehfid,
c     +544)
c             d2_p = 1.0d0 -  d1_p
c          endif
          do g_i_ = 1, g_p_
            g_cs(g_i_) = d2_p * g_d1_w(g_i_) + d1_p * g_cs(g_i_)
          enddo
          cs = d3_v
C--------
C
C INTERACTION REGION
C     The LiH singlet
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-belihc) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-belihc) * (r(i, 1) - relih)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xelih(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xelih = d2_v
C--------
          d2_v = xelih * xelih
          d1_p = 2.0d0 * xelih
          d4_b = s1c2 + s1c1 * d1_p
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d4_b * g_xelih(g_i_)
          enddo
          s1 = s1c1 * d2_v + s1c2 * xelih
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bes1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-bes1) * (r(i, 1) - relih)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d3_b = s1c3 * d1_p
          do g_i_ = 1, g_p_
            g_s1c(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          s1c = s1c3 * d2_v
C--------
          do g_i_ = 1, g_p_
            g_r1r3(g_i_) = r1r3al * g_r(g_i_, i, 2) + (-g_r(g_i_, i, 3))
     * + g_r(g_i_, i, 1)
          enddo
          r1r3 = r(i, 1) - r(i, 3) + r1r3al * r(i, 2)
C--------
          d3_b = 1.0d0 / rcut2del
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r1r3(g_i_)
          enddo
          d1_w = (r1r3 - rcut2rho) / rcut2del
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_rcut2(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          rcut2 = 0.5d0 * (1.0d0 - d2_v)
C--------
          do g_i_ = 1, g_p_
            g_s1(g_i_) = rcut2 * g_s1c(g_i_) + s1c * g_rcut2(g_i_) + g_s
     *1(g_i_)
          enddo
          s1 = s1 + rcut2 * s1c
C--------
C
C     The LiH triplet
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-bet1) * (r(i, 1) - reliht)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xelih(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xelih = d2_v
C--------
          d2_v = xelih * xelih
          d1_p = 2.0d0 * xelih
          d4_b = t1c2 + t1c1 * d1_p
          do g_i_ = 1, g_p_
            g_t1(g_i_) = d4_b * g_xelih(g_i_)
          enddo
          t1 = t1c1 * d2_v + t1c2 * xelih
C--------
C
C     The HF singlet
          do g_i_ = 1, g_p_
            g_yhf(g_i_) = g_r(g_i_, i, 2)
          enddo
          yhf = r(i, 2) - 2.1042d0
C--------
          d4_v = yhf * yhf
          d1_p = 2.0d0 * yhf
          d5_b = 0.059062d0 * d1_p + (-0.025647d0)
          do g_i_ = 1, g_p_
            g_behf(g_i_) = d5_b * g_yhf(g_i_)
          enddo
          behf = 1.1622d0 - 0.025647d0 * yhf + 0.059062d0 * d4_v
C--------
          d4_v = r(i, 2) - rehf
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-behf) * g_r(g_i_, i, 2) + (-d4_v) * g_behf(
     *g_i_)
          enddo
          d1_w = (-behf) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xehfa(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xehfa = d2_v
C--------
C
          d2_v = dehf * xehfa
          d3_v = xehfa - 2.0d0
          d5_b = d2_v + d3_v * dehf
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d5_b * g_xehfa(g_i_)
          enddo
          s2 = d2_v * d3_v + uli2p
C--------
          d3_b = 1.0d0 / 0.5d0
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (r(i, 2) - 3.d0) / 0.5d0
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_s2cut(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          s2cut = 0.5d0 - 0.5d0 * d2_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s2cc_b) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-s2cc_b) * (r(i, 2) - s2cc_r)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d4_v = (1.d0 - d2_v) * (1.d0 - d2_v)
          d1_p = 2.0d0 * (1.d0 - d2_v)
          d6_b = (-(s2cc_d * d1_p)) * d2_p
          do g_i_ = 1, g_p_
            g_s2cc(g_i_) = d6_b * g_d1_w(g_i_)
          enddo
          s2cc = s2cc_d * d4_v - s2cc_d
C--------
          d3_v = s2 - s2cc
          d2_b = 1.0d0 + (-s2cut)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d3_v * g_s2cut(g_i_) + s2cut * g_s2(g_i_) + d2_
     *b * g_s2cc(g_i_)
          enddo
          s2 = s2cc + d3_v * s2cut
C--------
C
          d5_b = dehfc * 0.5d0
          do g_i_ = 1, g_p_
            g_dehfm(g_i_) = d5_b * g_cs(g_i_)
          enddo
          dehfm = dehf - (dehfcc + dehfc * 0.5d0 * (1.0d0 - cs))
C--------
          d4_v = r(i, 2) - rehf
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-behf) * g_r(g_i_, i, 2) + (-d4_v) * g_behf(
     *g_i_)
          enddo
          d1_w = (-behf) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xehfm(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xehfm = d2_v
C--------
C
          d3_v = dehfm * xehfm
          d4_v = xehfm - 2.0d0
          d6_b = d4_v * xehfm
          d5_b = d3_v + d4_v * dehfm
          do g_i_ = 1, g_p_
            g_s2m(g_i_) = d5_b * g_xehfm(g_i_) + d6_b * g_dehfm(g_i_)
          enddo
          s2m = d3_v * d4_v + uli2p
C--------
          d3_b = 1.0d0 / 0.5d0
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (r(i, 2) - 2.2d0) / 0.5d0
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_s2cut(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          s2cut = 0.5d0 - 0.5d0 * d2_v
C--------
          do g_i_ = 1, g_p_
            g_s2m(g_i_) = s2m * g_s2cut(g_i_) + s2cut * g_s2m(g_i_)
          enddo
          s2m = s2m * s2cut
C--------
C
          d3_b = 1.0d0 / s2alp
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 3)
          enddo
          d1_w = (r(i, 3) - s2rho) / s2alp
          d4_v = (s2 - s2m) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.d0 + d6_v
          d7_b = d4_v * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s2(g_i_) + d2_b 
     ** g_s2m(g_i_)
          enddo
          s2 = s2m + d4_v * d7_v
C--------
C
C     The HF triplet
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet2) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-bet2) * (r(i, 2) - rehft)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xehfa(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xehfa = d2_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet2_180) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-bet2_180) * (r(i, 2) - rehft)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xehfa_180(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xehfa_180 = d2_v
C--------
C
          d2_v = xehfa * xehfa
          d2_p = 2.0d0 * xehfa
          d6_v = (t2a1 * d2_v + t2a2 * xehfa) * 0.5d0
          d8_v = 1.0d0 + cs
          d11_v = xehfa_180 * xehfa_180
          d1_p = 2.0d0 * xehfa_180
          d15_v = (t2a3 * d11_v + t2a4 * xehfa_180) * 0.5d0
          d16_v = 1.0d0 - cs
          d7_b = d16_v * 0.5d0
          d10_b = d7_b * t2a4 + d7_b * t2a3 * d1_p
          d6_b = -d15_v + d6_v
          d14_b = d8_v * 0.5d0
          d17_b = d14_b * t2a2 + d14_b * t2a1 * d2_p
          do g_i_ = 1, g_p_
            g_t2(g_i_) = d10_b * g_xehfa_180(g_i_) + d6_b * g_cs(g_i_) +
     * d17_b * g_xehfa(g_i_)
          enddo
          t2 = d6_v * d8_v + d15_v * d16_v
C--------
C
C     The LiF singlet
          d5_v = (4.6498d0 + 7.0489d0 * r(i, 3)) * (4.6498d0 + 7.0489d0 
     +* r(i, 3))
          d1_p = 2.0d0 * (4.6498d0 + 7.0489d0 * r(i, 3))
          d6_v = r(i, 3) * 103.57d0 / d5_v
          d7_b = (-d6_v) / d5_v * d1_p * 7.0489d0 + 1.0d0 / d5_v * 103.5
     *7d0
          do g_i_ = 1, g_p_
            g_belif(g_i_) = d7_b * g_r(g_i_, i, 3)
          enddo
          belif = d6_v + 0.064076d0
C--------
C
          d3_b = 1.0d0 / gdel
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (r(i, 2) - grho) / gdel
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_gswt(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          gswt = 0.5d0 - 0.5d0 * d2_v
C--------
          d2_v = gswt * 0.5d0
          d4_v = 1.0d0 + cs
          d5_b = d4_v * 0.5d0
          do g_i_ = 1, g_p_
            g_gswt(g_i_) = d2_v * g_cs(g_i_) + d5_b * g_gswt(g_i_)
          enddo
          gswt = d2_v * d4_v
C--------
          d3_v = belifff - belif
          d2_b = 1.0d0 + (-gswt)
          do g_i_ = 1, g_p_
            g_beliff(g_i_) = d3_v * g_gswt(g_i_) + d2_b * g_belif(g_i_)
          enddo
          beliff = belif + gswt * d3_v
C--------
C

          if (  nb .ne. 0 ) then
             d5_v = (r(i, 3) / gam) ** ( nb-1)
             d2_p =  nb *  d5_v
             d5_v =  d5_v * (r(i, 3) / gam)
          else
C            Maybe this should be  d5_v = (r(i, 3) / gam) **  nb
             d5_v = 1.0d0
             d2_p = 0.0d0
          endif
          d6_v = belif + d5_v

          if (  nb .ne. 0 ) then
             d9_v = (r(i, 3) / gam) ** ( nb-1)
             d1_p =  nb *  d9_v
             d9_v =  d9_v * (r(i, 3) / gam)
          else
C            Maybe this should be  d9_v = (r(i, 3) / gam) **  nb
             d9_v = 1.0d0
             d1_p = 0.0d0
          endif
          d10_v = beliff + d9_v
          d11_v = beliff * d6_v / d10_v
          d2_b = 1.0d0 / d10_v
          d3_b = (-d11_v) / d10_v
          d4_b = d3_b + d2_b * d6_v
          d8_b = d2_b * beliff
          d7_b = d3_b * d1_p * (1.0d0 / gam) + d8_b * d2_p * (1.0d0 / ga
     *m)
          do g_i_ = 1, g_p_
            g_belifmod(g_i_) = d7_b * g_r(g_i_, i, 3) + d8_b * g_belif(g
     *_i_) + d4_b * g_beliff(g_i_)
          enddo
          belifmod = d11_v
C--------
C
          d4_v = r(i, 3) - relif
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-belifmod) * g_r(g_i_, i, 3) + (-d4_v) * g_b
     *elifmod(g_i_)
          enddo
          d1_w = (-belifmod) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xelif(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xelif = d2_v
C--------
          d2_v = (delif + uli2p) * xelif
          d3_v = xelif - 2.0d0
          d5_b = d2_v + d3_v * (delif + uli2p)
          do g_i_ = 1, g_p_
            g_s3a(g_i_) = d5_b * g_xelif(g_i_)
          enddo
          s3a = d2_v * d3_v + uli2p
C--------
          d3_b = 1.0d0 / 0.5d0
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 3)
          enddo
          d1_w = (r(i, 3) - 13.d0) / 0.5d0
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_s3cut(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          s3cut = 0.5d0 - 0.5d0 * d2_v
C--------
          do g_i_ = 1, g_p_
            g_s3a(g_i_) = s3a * g_s3cut(g_i_) + s3cut * g_s3a(g_i_)
          enddo
          s3a = s3a * s3cut
C--------
C
          d3_b = 1.0d0 / swit_a
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (r(i, 2) - swit_r) / swit_a
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d5_b = (1.d0 - swit_d) * 0.5d0 * d1_p
          do g_i_ = 1, g_p_
            g_swit(g_i_) = d5_b * g_d1_w(g_i_)
          enddo
          swit = swit_d + (1.d0 - swit_d) * 0.5d0 * (1.d0 + d2_v)
C--------
          do g_i_ = 1, g_p_
            g_s3a(g_i_) = s3a * g_swit(g_i_) + swit * g_s3a(g_i_)
          enddo
          s3a = s3a * swit
C--------
C
          d5_b = delifc2 * 0.5d0 + (-(delifc1 * 0.5d0))
          do g_i_ = 1, g_p_
            g_delifm(g_i_) = d5_b * g_cs(g_i_)
          enddo
          delifm = delifc1 * 0.5d0 * (1.0d0 - cs) + delifc2 * 0.5d0 * (1
     *.0d0 + cs)
C--------
C
          d5_b = belifc2 * 0.5d0 + (-(belifc1 * 0.5d0))
          do g_i_ = 1, g_p_
            g_belifm(g_i_) = d5_b * g_cs(g_i_)
          enddo
          belifm = belifc1 * 0.5d0 * (1.0d0 - cs) + belifc2 * 0.5d0 * (1
     *.0d0 + cs)
C--------
C
          d4_v = r(i, 3) - relif
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-belifm) * g_r(g_i_, i, 3) + (-d4_v) * g_bel
     *ifm(g_i_)
          enddo
          d1_w = (-belifm) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xelif(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xelif = d2_v
C--------
          d2_v = delifm + uli2p
          d4_v = d2_v * xelif
          d5_v = xelif - 2.0d0
          d4_b = d4_v + d5_v * d2_v
          d6_b = d5_v * xelif
          do g_i_ = 1, g_p_
            g_s3c(g_i_) = d4_b * g_xelif(g_i_) + d6_b * g_delifm(g_i_)
          enddo
          s3c = d4_v * d5_v
C--------
C
          d4_b = -((rho2_180 - rho2_0) * 0.5d0)
          do g_i_ = 1, g_p_
            g_rho2(g_i_) = d4_b * g_cs(g_i_)
          enddo
          rho2 = rho2_0 + (rho2_180 - rho2_0) * 0.5d0 * (1.0d0 - cs)
C--------
          d4_b = -((theta2_180 - theta2_0) * 0.5d0)
          do g_i_ = 1, g_p_
            g_theta2(g_i_) = d4_b * g_cs(g_i_)
          enddo
          theta2 = theta2_0 + (theta2_180 - theta2_0) * 0.5d0 * (1.0d0 -
     * cs)
C--------
          d4_b = -((del2_180 - del2_0) * 0.5d0)
          do g_i_ = 1, g_p_
            g_del2(g_i_) = d4_b * g_cs(g_i_)
          enddo
          del2 = del2_0 + (del2_180 - del2_0) * 0.5d0 * (1.0d0 - cs)
C--------
C
          d2_v = cos(theta2)
          d2_p = -sin(theta2)
          d4_v = rehf - r(i, 2)
          d6_v = sin(theta2)
          d1_p = cos(theta2)
          d9_v = rho2 - r(i, 3)
          d8_b = (-d9_v) * d1_p + d4_v * d2_p
          do g_i_ = 1, g_p_
            g_gc2(g_i_) = d6_v * g_r(g_i_, i, 3) + (-d6_v) * g_rho2(g_i_
     *) + (-d2_v) * g_r(g_i_, i, 2) + d8_b * g_theta2(g_i_)
          enddo
          gc2 = d2_v * d4_v - d6_v * d9_v
C--------
          d3_v = gc2 / del2
          d2_b = 1.0d0 / del2
          d3_b = (-d3_v) / del2
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_del2(g_i_) + d2_b * g_gc2(g_i_)
          enddo
          d1_w = d3_v
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_cut2(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          cut2 = 0.5d0 - 0.5d0 * d2_v
C--------
C
          d3_v = s3a - s3c
          d2_b = 1.0d0 + (-cut2)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d3_v * g_cut2(g_i_) + cut2 * g_s3a(g_i_) + d2_b
     * * g_s3c(g_i_)
          enddo
          s3 = s3c + d3_v * cut2
C--------
C
C     The LiF triplet
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet3) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bet3) * (r(i, 3) - relif)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xelif(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xelif = d2_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet3_180) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bet3_180) * (r(i, 3) - relif)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xelif_180(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xelif_180 = d2_v
C--------
          d2_v = xelif * xelif
             d2_p = 2.0d0 * xelif
          d6_v = (t3a1 * d2_v + t3a2 * xelif) * 0.5d0
          d8_v = 1.0d0 + cs
          d11_v = xelif_180 * xelif_180
             d1_p = 2.0d0 * xelif_180
          d15_v = (t3a3 * d11_v + t3a4 * xelif_180) * 0.5d0
          d16_v = 1.0d0 - cs
          d7_b = d16_v * 0.5d0
          d10_b = d7_b * t3a4 + d7_b * t3a3 * d1_p
          d6_b = -d15_v + d6_v
          d14_b = d8_v * 0.5d0
          d17_b = d14_b * t3a2 + d14_b * t3a1 * d2_p
          do g_i_ = 1, g_p_
            g_t3(g_i_) = d10_b * g_xelif_180(g_i_) + d6_b * g_cs(g_i_) +
     * d17_b * g_xelif(g_i_)
          enddo
          t3 = d6_v * d8_v + d15_v * d16_v
C--------
C
C     LEPS
          do g_i_ = 1, g_p_
            g_coul1(g_i_) = 0.5d0 * g_t1(g_i_) + 0.5d0 * g_s1(g_i_)
          enddo
          coul1 = 0.5d0 * (s1 + t1)
C--------
          do g_i_ = 1, g_p_
            g_coul2(g_i_) = 0.5d0 * g_t2(g_i_) + 0.5d0 * g_s2(g_i_)
          enddo
          coul2 = 0.5d0 * (s2 + t2)
C--------
          do g_i_ = 1, g_p_
            g_coul3(g_i_) = 0.5d0 * g_t3(g_i_) + 0.5d0 * g_s3(g_i_)
          enddo
          coul3 = 0.5d0 * (s3 + t3)
C--------
          do g_i_ = 1, g_p_
            g_exch1(g_i_) = (-0.5d0) * g_t1(g_i_) + 0.5d0 * g_s1(g_i_)
          enddo
          exch1 = 0.5d0 * (s1 - t1)
C--------
          do g_i_ = 1, g_p_
            g_exch2(g_i_) = (-0.5d0) * g_t2(g_i_) + 0.5d0 * g_s2(g_i_)
          enddo
          exch2 = 0.5d0 * (s2 - t2)
C--------
          do g_i_ = 1, g_p_
            g_exch3(g_i_) = (-0.5d0) * g_t3(g_i_) + 0.5d0 * g_s3(g_i_)
          enddo
          exch3 = 0.5d0 * (s3 - t3)
C--------
          d4_v = (exch1 - exch2) * (exch1 - exch2)
          d3_p = 2.0d0 * (exch1 - exch2)
          d7_v = (exch2 - exch3) * (exch2 - exch3)
          d2_p = 2.0d0 * (exch2 - exch3)
          d10_v = (exch3 - exch1) * (exch3 - exch1)
          d1_p = 2.0d0 * (exch3 - exch1)
          d5_b = d1_p + (-d2_p)
          d6_b = -d1_p + d3_p
          d10_b = d2_p + (-d3_p)
          do g_i_ = 1, g_p_
            g_w(g_i_) = d5_b * g_exch3(g_i_) + d10_b * g_exch2(g_i_) + d
     *6_b * g_exch1(g_i_)
          enddo
          w = d4_v + d7_v + d10_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-c2c) * g_r(g_i_, i, 3) + (-c2c) * g_r(g_i_,
     * i, 2) + (-c2c) * g_r(g_i_, i, 1) + (-c2b) * g_w(g_i_)
          enddo
          d1_w = (-c2b) * w - c2c * (r(i, 1) + r(i, 2) + r(i, 3))
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d3_b = c2a * d1_p
          do g_i_ = 1, g_p_
            g_cplg2(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          cplg2 = c2a * d2_v
C--------
          d4_b = cplg2 + cplg2
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_cplg2(g_i_) + g_w(g_i_)
          enddo
          d1_w = w + cplg2 * cplg2
          d7_v = sqrt(d1_w)

c          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d7_v)
c          else
c             call ehufDO (9,d1_w, d7_v, d1_p,
c     +g_ehfid,
c     +1165)
c          endif
          d6_b = (-(1.0d0 / dsqrt(2.0d0))) * d1_p
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = d6_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2
     *(g_i_) + g_coul1(g_i_)
          enddo
          e(i) = coul1 + coul2 + coul3 - d7_v / dsqrt(2.0d0) + dehf
C--------
C
C LONG-RANGE FORCES * * * * * * * * * * * * * * * * * * * * * * * * * * *
C HF DIPOLE
          d4_b = (-hfdip_a) * cauang
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (-hfdip_a) * (r(i, 2) * cauang - hfdip_re)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_hfdip_y(g_i_) = (-d1_p) * g_d1_w(g_i_)
          enddo
          hfdip_y = 1.d0 - d2_v
C--------
          do g_i_ = 1, g_p_
            g_hfdip(g_i_) = 0.0d0
          enddo
          hfdip = 0.d0
C--------
          do j = 1, 11
            d3_v = dble(j - 1)

c            if ( hfdip_y .ne. 0.0d0 ) then
               d4_v = hfdip_y ** ( d3_v - 2.0d0)
               d4_v =  d4_v * hfdip_y
               d1_p =  d3_v *  d4_v
               d4_v =  d4_v * hfdip_y
c            else
C              (hfdip_y = 0)
c               d4_v = hfdip_y **  d3_v
c
c               if (  d3_v .lt. 1.0d0 ) then
c                  call ehbfDO (10,hfdip_y, d3_v, d4_v, d1_p, 0.0d0,
c     +g_ehfid,
c     +1209)
c               else if (  d3_v .lt. 2.0d0 ) then
c                  d1_p = 0.0d0
c                  call ehbfDO (10,hfdip_y, d3_v, d4_v, d1_p, 0.0d0,
c     +g_ehfid,
c     +1214)
c               else
c                  d1_p = 0.0d0
c               endif
c            endif
            d6_b = 0.0d0
            d5_b = hfdip_m(j) * d1_p
            do g_i_ = 1, g_p_
              g_hfdip(g_i_) = d5_b * g_hfdip_y(g_i_) + g_hfdip(g_i_)
            enddo
            hfdip = hfdip + hfdip_m(j) * d4_v
C--------
          enddo
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = 20.d0 * g_r(g_i_, i, 2)
          enddo
          d1_w = 20.d0 * (r(i, 2) - 1.2d0)
          d2_v = hfdip * 0.5d0
          d4_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d4_v *  d4_v)
          d5_v = 1.d0 + d4_v
          d5_b = d2_v * d1_p
          d6_b = d5_v * 0.5d0
          do g_i_ = 1, g_p_
            g_hfdip(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_hfdip(g_i_)
          enddo
          hfdip = d2_v * d5_v
C--------
C
C HF POLARIZATION
          d2_v = coschihf * coschihf
          d2_p = 2.0d0 * coschihf
          d5_v = coschihf * coschihf
          d1_p = 2.0d0 * coschihf
          d5_b = hfpolpar * d1_p + (-hfpolperp) * d2_p
          do g_i_ = 1, g_p_
            g_hfpol(g_i_) = d5_b * g_coschihf(g_i_)
          enddo
          hfpol = hfpolperp * (1.d0 - d2_v) + hfpolpar * d5_v
C--------
          hfpolavg = 2.d0 * hfpolperp / 3.d0 + hfpolpar / 3.d0
C
C HF QUADRUPOLE
          d4_b = (-hfqua_a) * cauang
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_r(g_i_, i, 2)
          enddo
          d1_w = (-hfqua_a) * (r(i, 2) * cauang - hfqua_re)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_hfqua_y(g_i_) = (-d1_p) * g_d1_w(g_i_)
          enddo
          hfqua_y = 1.d0 - d2_v
C--------
          do g_i_ = 1, g_p_
            g_hfqua(g_i_) = 0.0d0
          enddo
          hfqua = 0.d0
C--------
          do j = 1, 10
            d3_v = dble(j - 1)

c            if ( hfqua_y .ne. 0.0d0 ) then
               d4_v = hfqua_y ** ( d3_v - 2.0d0)
               d4_v =  d4_v * hfqua_y
               d1_p =  d3_v *  d4_v
               d4_v =  d4_v * hfqua_y
c            else
C              (hfqua_y = 0)
c               d4_v = hfqua_y **  d3_v
c
c               if (  d3_v .lt. 1.0d0 ) then
c                  call ehbfDO (10,hfqua_y, d3_v, d4_v, d1_p, 0.0d0,
c     +g_ehfid,
c     +1289)
c               else if (  d3_v .lt. 2.0d0 ) then
c                  d1_p = 0.0d0
c                  call ehbfDO (10,hfqua_y, d3_v, d4_v, d1_p, 0.0d0,
c     +g_ehfid,
c     +1294)
c               else
c                  d1_p = 0.0d0
c               endif
c            endif
            d6_b = 0.0d0
            d5_b = hfqua_m(j) * d1_p
            do g_i_ = 1, g_p_
              g_hfqua(g_i_) = d5_b * g_hfqua_y(g_i_) + g_hfqua(g_i_)
            enddo
            hfqua = hfqua + hfqua_m(j) * d4_v
C--------
          enddo
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = 20.d0 * g_r(g_i_, i, 2)
          enddo
          d1_w = 20.d0 * (r(i, 2) - 1.2d0)
          d2_v = hfqua * 0.5d0
          d4_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d4_v *  d4_v)
          d5_v = 1.d0 + d4_v
          d5_b = d2_v * d1_p
          d6_b = d5_v * 0.5d0
          do g_i_ = 1, g_p_
            g_hfqua(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_hfqua(g_i_)
          enddo
          hfqua = d2_v * d5_v
C--------
C
C LiF DIPOLE
          do g_i_ = 1, g_p_
            g_lifdip_y(g_i_) = g_r(g_i_, i, 3)
          enddo
          lifdip_y = r(i, 3) - lifdip_4
C--------
          d2_v = lifdip_y ** ( 8 - 2)
             d2_v =  d2_v * lifdip_y
                d1_p =  8 *  d2_v
                d2_v =  d2_v * lifdip_y
          d4_b = lifdip_3 + lifdip_2 * d1_p
          do g_i_ = 1, g_p_
            g_lifdip_beta(g_i_) = d4_b * g_lifdip_y(g_i_)
          enddo
          lifdip_beta = lifdip_1 + lifdip_2 * d2_v + lifdip_3 * lifdip_y
C--------
          d5_v = (r(i, 3) - lifdip_re) * (r(i, 3) - lifdip_re)
             d1_p = 2.0d0 * (r(i, 3) - lifdip_re)
          d5_b = (-lifdip_beta) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_r(g_i_, i, 3) + (-d5_v) * g_lifdip_b
     *eta(g_i_)
          enddo
          d1_w = (-lifdip_beta) * d5_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d3_b = lifdip_d * d1_p
          do g_i_ = 1, g_p_
            g_lifdip(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          lifdip = lifdip_d * d2_v
C--------
C
C Li(2p) POLARIZATION
          d2_v = coschihf * coschihf
          d2_p = 2.0d0 * coschihf
          d5_v = coschihf * coschihf
          d1_p = 2.0d0 * coschihf
          d5_b = li2ppolpar * d1_p + (-li2ppolperp) * d2_p
          do g_i_ = 1, g_p_
            g_li2ppol(g_i_) = d5_b * g_coschihf(g_i_)
          enddo
          li2ppol = li2ppolperp * (1.d0 - d2_v) + li2ppolpar * d5_v
C--------
          li2ppolavg = 2.d0 * li2ppolperp / 3.d0 + li2ppolpar / 3.d0
C
C HF + Li(2p)
C     DIPOLE-INDUCED-DIPOLE
          d2_v = coschihf * coschihf
          d2_p = 2.0d0 * coschihf
          d5_v = (-0.5d0) * (3.d0 * d2_v + 1.d0)
          d7_v = d5_v * li2ppol
          d9_v = hfdip * hfdip
          d1_p = 2.0d0 * hfdip
          d4_b = d7_v * d1_p
          d6_b = d9_v * d5_v
          d10_b = d9_v * li2ppol * (-0.5d0) * 3.d0 * d2_p
          do g_i_ = 1, g_p_
            g_hf_did(g_i_) = d4_b * g_hfdip(g_i_) + d6_b * g_li2ppol(g_i
     *_) + d10_b * g_coschihf(g_i_)
          enddo
          hf_did = d7_v * d9_v
C--------
          d3_v = bigrhf ** ( 7 - 2)
          d3_v =  d3_v * bigrhf
          d1_p =  7 *  d3_v
          d3_v =  d3_v * bigrhf
          d4_v = d3_v + rhf0did ** 7
          d5_v = hf_did / d4_v
          d6_v = d5_v / cevau
          d4_b = bigrhf * (1.0d0 / cevau)
          d5_b = d4_b * (1.0d0 / d4_v)
          d3_b = d6_v + d4_b * ((-d5_v) / d4_v) * d1_p
          do g_i_ = 1, g_p_
            g_hf_did(g_i_) = d3_b * g_bigrhf(g_i_) + d5_b * g_hf_did(g_i
     *_)
          enddo
          hf_did = d6_v * bigrhf
C--------
C     cut off angular dependence in interaction region
          d3_v = (bigrhf / rhf0did) ** ( 4 - 2)
          d3_v =  d3_v * (bigrhf / rhf0did)
          d1_p =  4 *  d3_v
          d3_v =  d3_v * (bigrhf / rhf0did)
          d4_b = (-d1_p) * (1.0d0 / rhf0did)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_bigrhf(g_i_)
          enddo
          d1_w = -d3_v
          d3_v = coschihf * coschihf
          d2_p = 2.0d0 * coschihf
          d5_v = 3.d0 * d3_v + 1.d0
          d6_v = 2.5d0 * li2ppolavg / d5_v
          d8_v = d6_v / li2ppol
          d9_v = d8_v - 1.d0
          d11_v = exp(d1_w)
          d1_p =  d11_v
          d13_v = 1.d0 + d9_v * d11_v
          d7_b = hf_did * d9_v * d1_p
          d8_b = hf_did * d11_v
          d10_b = d8_b * ((-d8_v) / li2ppol)
          d14_b = d8_b * (1.0d0 / li2ppol) * ((-d6_v) / d5_v) * 3.d0 * d
     *2_p
          do g_i_ = 1, g_p_
            g_hf_did(g_i_) = d7_b * g_d1_w(g_i_) + d10_b * g_li2ppol(g_i
     *_) + d14_b * g_coschihf(g_i_) + d13_v * g_hf_did(g_i_)
          enddo
          hf_did = hf_did * d13_v
C--------
C
C     DISPERSION
          d2_v = (-1.5d0) * li2pie * hfie * li2ppol
          d2_b = 1.0d0 / (li2pie + hfie)
          d4_b = d2_b * d2_v
          d5_b = d2_b * hfpol * ((-1.5d0) * li2pie * hfie)
          do g_i_ = 1, g_p_
            g_hf_disp(g_i_) = d4_b * g_hfpol(g_i_) + d5_b * g_li2ppol(g_
     *i_)
          enddo
          hf_disp = d2_v * hfpol / (li2pie + hfie)
C--------
          d3_v = bigrhf ** ( 7 - 2)
          d3_v =  d3_v * bigrhf
          d1_p =  7 *  d3_v
          d3_v =  d3_v * bigrhf
          d4_v = d3_v + rhf0disp ** 7
          d5_v = hf_disp / d4_v
          d4_b = bigrhf * (1.0d0 / d4_v)
          d3_b = d5_v + bigrhf * ((-d5_v) / d4_v) * d1_p
          do g_i_ = 1, g_p_
            g_hf_disp(g_i_) = d3_b * g_bigrhf(g_i_) + d4_b * g_hf_disp(g
     *_i_)
          enddo
          hf_disp = d5_v * bigrhf
C--------
C     cut off angular dependence in interaction region
          d3_v = (bigrhf / rhf0disp) ** ( 4 - 2)
          d3_v =  d3_v * (bigrhf / rhf0disp)
          d1_p =  4 *  d3_v
          d3_v =  d3_v * (bigrhf / rhf0disp)
          d4_b = (-d1_p) * (1.0d0 / rhf0disp)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_bigrhf(g_i_)
          enddo
          d1_w = -d3_v
          d3_v = li2ppolavg * hfpolavg / hfpol
          d5_v = d3_v / li2ppol
          d6_v = d5_v - 1.d0
          d8_v = exp(d1_w)
          d1_p =  d8_v
          d10_v = 1.d0 + d6_v * d8_v
          d7_b = hf_disp * d6_v * d1_p
          d8_b = hf_disp * d8_v
          d10_b = d8_b * ((-d5_v) / li2ppol)
          d11_b = d8_b * (1.0d0 / li2ppol) * ((-d3_v) / hfpol)
          do g_i_ = 1, g_p_
            g_hf_disp(g_i_) = d7_b * g_d1_w(g_i_) + d10_b * g_li2ppol(g_
     *i_) + d11_b * g_hfpol(g_i_) + d10_v * g_hf_disp(g_i_)
          enddo
          hf_disp = hf_disp * d10_v
C--------
C
C     QUADRUPOLE-QUADRUPOLE
          d3_v = 3.d0 / 4.d0 * hfqua * li2pqua
          d5_v = coschihf * coschihf
          d1_p = 2.0d0 * coschihf
          d7_v = 3.d0 - 7.d0 * d5_v
          d6_b = (-d3_v) * 7.d0 * d1_p
          d8_b = d7_v * li2pqua * (3.d0 / 4.d0)
          do g_i_ = 1, g_p_
            g_hf_qq(g_i_) = d6_b * g_coschihf(g_i_) + d8_b * g_hfqua(g_i
     *_)
          enddo
          hf_qq = d3_v * d7_v
C--------
          d3_v = bigrhf ** ( 6 - 2)
          d3_v =  d3_v * bigrhf
          d1_p =  6 *  d3_v
          d3_v =  d3_v * bigrhf
          d4_v = d3_v + rhf0qq ** 6
          d5_v = hf_qq / d4_v
          d6_v = d5_v / cevau
          d4_b = bigrhf * (1.0d0 / cevau)
          d5_b = d4_b * (1.0d0 / d4_v)
          d3_b = d6_v + d4_b * ((-d5_v) / d4_v) * d1_p
          do g_i_ = 1, g_p_
            g_hf_qq(g_i_) = d3_b * g_bigrhf(g_i_) + d5_b * g_hf_qq(g_i_)
          enddo
          hf_qq = d6_v * bigrhf
C--------
C     cut off angular dependence in interaction region
          d3_v = (bigrhf / rhf0qq) ** ( 4 - 2)
          d3_v =  d3_v * (bigrhf / rhf0qq)
          d1_p =  4 *  d3_v
          d3_v =  d3_v * (bigrhf / rhf0qq)
          d4_b = (-d1_p) * (1.0d0 / rhf0qq)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_bigrhf(g_i_)
          enddo
          d1_w = -d3_v
          d3_v = coschihf * coschihf
          d2_p = 2.0d0 * coschihf
          d5_v = 3.d0 - 7.d0 * d3_v
          d6_v = (-.5d0) / d5_v
          d7_v = d6_v - 1.d0
          d9_v = exp(d1_w)
          d1_p =  d9_v
          d11_v = 1.d0 + d7_v * d9_v
          d7_b = hf_qq * d7_v * d1_p
          d12_b = (-(hf_qq * d9_v * ((-d6_v) / d5_v))) * 7.d0 * d2_p
          do g_i_ = 1, g_p_
            g_hf_qq(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_coschihf(g_i
     *_) + d11_v * g_hf_qq(g_i_)
          enddo
          hf_qq = hf_qq * d11_v
C--------
C
C     DIPOLE-QUADRUPOLE
          d3_v = 3.d0 / 4.d0 * hfdip * li2pqua
          d5_b = coschihf * li2pqua * (3.d0 / 4.d0)
          do g_i_ = 1, g_p_
            g_hf_dq(g_i_) = d3_v * g_coschihf(g_i_) + d5_b * g_hfdip(g_i
     *_)
          enddo
          hf_dq = d3_v * coschihf
C--------
          d3_v = bigrhf ** ( 5 - 2)
          d3_v =  d3_v * bigrhf
          d1_p =  5 *  d3_v
          d3_v =  d3_v * bigrhf
          d4_v = d3_v + rhf0dq ** 5
          d5_v = hf_dq / d4_v
          d6_v = d5_v / cevau
          d4_b = bigrhf * (1.0d0 / cevau)
          d5_b = d4_b * (1.0d0 / d4_v)
          d3_b = d6_v + d4_b * ((-d5_v) / d4_v) * d1_p
          do g_i_ = 1, g_p_
            g_hf_dq(g_i_) = d3_b * g_bigrhf(g_i_) + d5_b * g_hf_dq(g_i_)
          enddo
          hf_dq = d6_v * bigrhf
C--------
C     cut off angular dependence in interaction region
          d3_v = (bigrhf / rhf0dq) ** ( 4 - 2)
          d3_v =  d3_v * (bigrhf / rhf0dq)
          d1_p =  4 *  d3_v
          d3_v =  d3_v * (bigrhf / rhf0dq)
          d4_b = (-d1_p) * (1.0d0 / rhf0dq)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_bigrhf(g_i_)
          enddo
          d1_w = -d3_v
          d3_v = exp(d1_w)
          d1_p =  d3_v
          d4_b = hf_qq * d1_p
          do g_i_ = 1, g_p_
            g_hf_qq(g_i_) = d4_b * g_d1_w(g_i_) + d3_v * g_hf_qq(g_i_)
          enddo
          hf_qq = hf_qq * d3_v
C--------
C
C     SUM FORCES
          do g_i_ = 1, g_p_
            g_e_lr_hf(g_i_) = g_hf_did(g_i_) + g_hf_disp(g_i_) + g_hf_qq
     *(g_i_) + g_hf_dq(g_i_)
          enddo
          e_lr_hf = hf_dq + hf_qq + hf_disp + hf_did
C--------
C     cut off at large HF distances
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = r2lrcut_a * g_r(g_i_, i, 2)
          enddo
          d1_w = r2lrcut_a * (r(i, 2) - r2lrcut_r)
          d2_v = e_lr_hf * 0.5d0
          d4_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d4_v *  d4_v)
          d5_v = 1.d0 - d4_v
          d5_b = (-d2_v) * d1_p
          d6_b = d5_v * 0.5d0
          do g_i_ = 1, g_p_
            g_e_lr_hf(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_e_lr_hf(g_i
     *_)
          enddo
          e_lr_hf = d2_v * d5_v
C--------
C
C LiF + H
C     DIPOLE-INDUCED-DIPOLE
          d2_v = coschilif * coschilif
             d2_p = 2.0d0 * coschilif
          d6_v = (-0.5d0) * (3.d0 * d2_v + 1.d0) * hpol
          d8_v = lifdip * lifdip
             d1_p = 2.0d0 * lifdip
          d4_b = d6_v * d1_p
          d9_b = d8_v * hpol * (-0.5d0) * 3.d0 * d2_p
          do g_i_ = 1, g_p_
            g_lif_did(g_i_) = d4_b * g_lifdip(g_i_) + d9_b * g_coschilif
     *(g_i_)
          enddo
          lif_did = d6_v * d8_v
C--------
          d3_v = bigrlif ** ( 6 - 2)
             d3_v =  d3_v * bigrlif
                d1_p =  6 *  d3_v
                d3_v =  d3_v * bigrlif
          d4_v = d3_v + rlif0did ** 6
          d5_v = lif_did / d4_v
          d2_b = 1.0d0 / cevau
          d3_b = d2_b * (1.0d0 / d4_v)
          d6_b = d2_b * ((-d5_v) / d4_v) * d1_p
          do g_i_ = 1, g_p_
            g_lif_did(g_i_) = d6_b * g_bigrlif(g_i_) + d3_b * g_lif_did(
     *g_i_)
          enddo
          lif_did = d5_v / cevau
C--------
C     cut off angular dependence in interaction region
          d3_v = (bigrlif / rlif0did) ** ( 4 - 2)
             d3_v =  d3_v * (bigrlif / rlif0did)
                d1_p =  4 *  d3_v
                d3_v =  d3_v * (bigrlif / rlif0did)
          d4_b = (-d1_p) * (1.0d0 / rlif0did)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_bigrlif(g_i_)
          enddo
          d1_w = -d3_v
          d3_v = coschilif * coschilif
             d2_p = 2.0d0 * coschilif
          d5_v = 3.d0 * d3_v + 1.d0
          d6_v = 2.5d0 / d5_v
          d7_v = d6_v - 1.d0
          d9_v = exp(d1_w)
          d1_p =  d9_v
          d11_v = 1.d0 + d7_v * d9_v
          d7_b = lif_did * d7_v * d1_p
          d12_b = lif_did * d9_v * ((-d6_v) / d5_v) * 3.d0 * d2_p
          do g_i_ = 1, g_p_
            g_lif_did(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_coschilif(
     *g_i_) + d11_v * g_lif_did(g_i_)
          enddo
          lif_did = lif_did * d11_v
C--------
C
C     DISPERSION
          do g_i_ = 1, g_p_
            g_lif_disp(g_i_) = 0.0d0
          enddo
          lif_disp = (-1.5d0) * lifie * hie * lifpol * hpol / (lifie + h
     *ie)
C--------
          d3_v = bigrlif ** ( 6 - 2)
             d3_v =  d3_v * bigrlif
                d1_p =  6 *  d3_v
                d3_v =  d3_v * bigrlif
          d4_v = d3_v + rlif0disp ** 6
          d5_v = lif_disp / d4_v
          d2_b = 1.0d0 / d4_v
          d5_b = (-d5_v) / d4_v * d1_p
          do g_i_ = 1, g_p_
            g_lif_disp(g_i_) = d5_b * g_bigrlif(g_i_) + d2_b * g_lif_dis
     *p(g_i_)
          enddo
          lif_disp = d5_v
C--------
C     no angular dependence
C
C     SUM FORCES
          do g_i_ = 1, g_p_
            g_e_lr_lif(g_i_) = g_lif_did(g_i_) + g_lif_disp(g_i_)
          enddo
          e_lr_lif = lif_disp + lif_did
C--------
C     cut off at large LiF distances
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = r3lrcut_a * g_r(g_i_, i, 3)
          enddo
          d1_w = r3lrcut_a * (r(i, 3) - r3lrcut_r)
          d2_v = e_lr_lif * 0.5d0
          d4_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d4_v *  d4_v)
          d5_v = 1.d0 - d4_v
          d5_b = (-d2_v) * d1_p
          d6_b = d5_v * 0.5d0
          do g_i_ = 1, g_p_
            g_e_lr_lif(g_i_) = d5_b * g_d1_w(g_i_) + d6_b * g_e_lr_lif(g
     *_i_)
          enddo
          e_lr_lif = d2_v * d5_v
C--------
C
C LiH + F
C no long-range forces included
C
C COMBINE INTERACTION REGION ENERGY AND LONG-RANGE FORCES
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_e_lr_hf(g_i_) + g_e_lr_lif(g_i_) + g_e(g_i_
     *, i)
          enddo
          e(i) = e(i) + e_lr_lif + e_lr_hf
C--------
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_e(g_i_, i)
          enddo
          e(i) = e(i) * cevau
C--------
C
10        continue
99999   continue
      end
C
C**************************************************************************
C**************************************************************************
