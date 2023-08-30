      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: nt, r(1,3), r2(3), v(1)
      double precision :: u11(1), u12(1), u22(1)
      double precision :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
      double precision :: u(nstates,nstates), dudr(nstates,nstates,3)
      double precision :: t(nstates,nstates)
      double precision :: tmpmat(nstates,nstates)
      double precision :: hr(nstates,nstates,3), gr(nstates,3)
      double precision :: tx(9), drdx(3,9)
      double precision :: hx(nstates,nstates,9), gx(nstates,9)
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

      do i=1,3
        tmpmat(:,:)=dudr(:,:,i)
        tmpmat=matmul(transpose(t), matmul(tmpmat,t))
        do j=1,nstates
          do k=j+1,nstates
            if ((p(k)-p(j)) > 1.d-8) then
              hr(j,k,i)=tmpmat(j,k)/(p(k)-p(j))
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
        g(:,iatom,idir)=gx(:,j)
        d(:,:,iatom,idir)=hx(:,:,j)
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
C   System:          MXH
C   Common name:     MXH-SL
C   Functional form: Modified extended LEPS
C   Number of derivatives: 1
C   Number of electronic surfaces: 2
C   Interface: 3-2V

C   Reference:  Y. L. Volobuev, M. D. Hack, M. S. Topaler, and D. G. Truhlar,
C               J. Chem. Phys. 112 9716 (2000).

C   Notes:
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
C      The classical potential energy is set equal to zero for M
C      infinitely far from the HX diatomic and R(HX) set equal to the
C      HX equilibrium diatomic value.
C
C   Parameters:
C
C   Coordinates:
C      Internal, Definition: R(1) = R(M-H)
C                            R(2) = R(H-X)
C                            R(3) = R(M-X))

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

      subroutine prepot

        implicit none
        integer g_p_
        parameter(g_p_=3)
        common /com_para/ g_myrank,g_nprocs
        integer g_myrank, g_nprocs


        call g_prepot11(g_p_)
        call g_prepot12(g_p_)
        call g_prepot22(g_p_)

      return

      end
c
c ********************************************************************
c ********************************************************************
      subroutine pot(r,e,de,nt,nsurf)

        implicit double precision (a-h,o-z)
        dimension r(nt,3),e(nt),de(3,nt)

        integer i,j,k,g_p_,ldg_r,nsurf,nt
        parameter (g_p_=3,ldg_r=3,ldg_e =3)
        dimension g_r(ldg_r,nt,3),g_e(ldg_e,nt)
        dimension de1(3,1),de2(3,1),de3(3,1),ddis(3)

c initialize seed matrix here

        do i=1,nt
          do j=1,3
            do k=1,3
              if (j .eq. k) then
                g_r(j,i,k)=1.0d0
              else
               g_r(j,i,k)=0.0d0
             endif
            enddo
          enddo
        enddo

        if (nsurf .eq. 1) then
          call g_pot11(g_p_, r, g_r, ldg_r, e, De, ldg_e, nt)
        else if (nsurf .eq. 2) then
          call g_pot12(g_p_, r, g_r, ldg_r, e, De, ldg_e, nt)
        else if (nsurf .eq. 3) then
          call g_pot22(g_p_, r, g_r, ldg_r, e, De, ldg_e, nt)
        else
          write(*,*)'Error!'
          stop
        endif

      return

      end
C ********************************************************************
C ********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 07/21/99 by the version of
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
C ********************************************************************
C ********************************************************************
C  U11 surface - lower diabatic surface
C
      subroutine g_prepot11(g_p_)
        implicit none
        double precision r, e, des1, bes1, res1, rehf, dehf, des3, bes3,
     * res3
        double precision s1, s2, s3, xes1, xehf, xes3, yhf, behf
        double precision cevau, cosa
        integer i, nt
C
        dimension r(nt, 3), e(nt)
C
        save
C MH singlet
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d1_p, d5_b, d2_v, d3_v, d4_b, d4_v, d1_w, g_d1_
     *w(g_pmax_), g_r(ldg_r, nt, 3), g_xes1(g_pmax_)
        double precision g_s1(g_pmax_), g_yhf(g_pmax_), g_behf(g_pmax_),
     * g_xehf(g_pmax_), g_s2(g_pmax_), g_xes3(g_pmax_), g_s3(g_pmax_), g
     *_e(ldg_e, nt)
        integer g_ehfid
        data des1 /1.67d0/, bes1 /1.7d0/, res1 /1.5d0/
C
C HF singlet
        data rehf /1.733d0/, dehf /6.122d0/
C
C MF singlet
        data des3 /0.1d0/, bes3 /1.2d0/, res3 /4.0d0/
C
C
C  misc constants
C
        data g_ehfid /0/
c
c Common blocks for the diatomic potential 
        double precision e_diat,de_diat_dr,e_asym
        common /sur_diapot/ e_diat(3,3),de_diat_dr(3,3),e_asym(3)

C
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        cevau = 1.d0 / 27.2113961d0

C
        return
C
        entry g_pot11(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do 99999 i = 1, nt
C
c          cosa = 0.5d0 * (r(i, 2) ** 2 + r(i, 3) ** 2 - r(i, 1) ** 2) / 
c     *(r(i, 2) * r(i, 3))
C==============================================================k
CThis section contains the diatomic curves
C
C The  MH singlet
          res1 = 2.0d0
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bes1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-bes1) * (r(i, 1) - res1)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xes1(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xes1 = d2_v
C--------
          d2_v = des1 * xes1
          d3_v = xes1 + 2.0d0
          d4_b = d2_v + d3_v * des1
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d4_b * g_xes1(g_i_)
          enddo
          s1 = d2_v * d3_v

          e_diat(1,1) = s1*cevau
          de_diat_dr(1,1) = g_s1(1)*cevau

C--------
C
C The HF singlet
C
          rehf = 1.733d0
          dehf = 6.122d0

          e_asym(1) = dehf*cevau
C
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
C
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

          e_diat(2,1) = s2*cevau
          de_diat_dr(2,1) = g_s2(2)*cevau

C--------
C
C The MF singlet
C
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
          des3 = 0.060
          d2_v = des3 * xes3
          d3_v = xes3 - 2.0d0
          d4_b = d2_v + d3_v * des3
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d4_b * g_xes3(g_i_)
          enddo
          s3 = d2_v * d3_v

          e_diat(3,1) = s3*cevau
          de_diat_dr(3,1) = g_s3(3)*cevau

C--------
C
C =========================================================
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_s3(g_i_) + g_s2(g_i_) + g_s1(g_i_)
          enddo
          e(i) = s1 + s2 + s3 + dehf
C--------
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_e(g_i_, i)
          enddo
          e(i) = e(i) * cevau
C--------
C
10        continue
99999   continue
        return
C
      end
C
C
C ********************************************************************
C ********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 07/21/99 by the version of
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
C ********************************************************************
C ********************************************************************
C
C U12 coupling surface
C
      subroutine g_prepot12(g_p_)
C
        implicit none
C
        double precision ap, alp, bet, rd0, r20, theta
        double precision gc, gc2, rnag, r1s, r2s, r3s
C
        double precision cevau
C
        integer nt, i
        double precision r(nt, 3), e(nt)
C
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d2_p, d2_w, d8_b, d7_b, d6_b, d5_b, d4_b, d2_v,
     * d4_v, d5_v
        double precision d1_p, d1_w, d3_v, g_r1s(g_pmax_), g_r(ldg_r, nt
     *, 3), g_r2s(g_pmax_), g_r3s(g_pmax_), g_d1_w(g_pmax_), g_rnag(g_pm
     *ax_), g_gc(g_pmax_)
        double precision g_gc2(g_pmax_), g_d2_w(g_pmax_), g_e(ldg_e, nt)
        integer g_ehfid
        data ap /0.3017d0/, alp /3.0d0/, rd0 /3.0d0/, bet /0.5d0/, r20 /
     *-1.3d0/

        double precision e_diat,de_diat_dr,e_asym
        common /sur_diapot/ e_diat(3,3),de_diat_dr(3,3),e_asym(3)
C
        save
C
        data g_ehfid /0/
C
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        theta = 15.6d0 * dacos(-1.0d0) / 180.0d0
C
        cevau = 1.d0 / 27.2113961d0
C
        e_asym(2) = 0.0d0
        e_diat(:,2) = 0.0d0
        de_diat_dr(:,2) = 0.0d0

        return
C
        entry g_pot12(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do i = 1, nt
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
C
C
          d4_b = -(2.0d0 / 9.0d0)
          d7_b = 1.0d0 / 3.0d0
          d8_b = 2.0d0 / 3.0d0
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_r2s(g_i_) + d7_b * g_r1s(g_i_) + d8_
     *b * g_r3s(g_i_)
          enddo
          d1_w = 2.0d0 / 3.0d0 * r3s + 1.0d0 / 3.0d0 * r1s - 2.0d0 / 9.0
     *d0 * r2s
          d2_v = sqrt(d1_w)

          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_rnag(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          rnag = d2_v
C--------
C
C
          d4_b = dsin(theta)
          d5_b = dcos(theta)
          do g_i_ = 1, g_p_
            g_gc(g_i_) = d4_b * g_r(g_i_, i, 2) + d5_b * g_rnag(g_i_)
          enddo
          gc = dcos(theta) * rnag + dsin(theta) * r(i, 2)
C--------
          d4_b = -dcos(theta)
          d5_b = dsin(theta)
          do g_i_ = 1, g_p_
            g_gc2(g_i_) = d4_b * g_r(g_i_, i, 2) + d5_b * g_rnag(g_i_)
          enddo
          gc2 = dsin(theta) * rnag - dcos(theta) * r(i, 2)
C--------
C
          d4_v = ((gc - rd0) / alp) ** ( 4 - 2)
          d4_v =  d4_v * ((gc - rd0) / alp)
          d1_p =  4 *  d4_v
          d4_v =  d4_v * ((gc - rd0) / alp)
          d5_b = (-d1_p) * (1.0d0 / alp)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_gc(g_i_)
          enddo
          d1_w = -d4_v
          d4_v = ((gc2 - r20) / bet) * ((gc2 - r20) / bet)
          d1_p = 2.0d0 * ((gc2 - r20) / bet)
          d5_b = (-d1_p) * (1.0d0 / bet)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = d5_b * g_gc2(g_i_)
          enddo
          d2_w = -d4_v
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d3_v = ap * d2_v
          d5_v = exp(d2_w)
          d1_p =  d5_v
          d4_b = d3_v * d1_p
          d6_b = d5_v * ap * d2_p
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = d4_b * g_d2_w(g_i_) + d6_b * g_d1_w(g_i_)
          enddo
          e(i) = d3_v * d5_v
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
C ********************************************************************
C ********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 07/21/99 by the version of
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
C ********************************************************************
C ********************************************************************
C
C the upper diabatic surface  U22
C
      subroutine g_prepot22(g_p_)
        implicit none
        double precision r, e, denah, renah, benah, rehf, dehf, demf, re
     *mf, bemf
        double precision c2a, c2b, c2c, um2p, cevau, onethd, twothd
        double precision s1, s2, s3, t1, t2, t3, s2a, s2b, s2m
        double precision coul1, coul2, coul3, exch1, exch2, exch3, cplg2
     *, w
        double precision xenah, xehfa, xehfb, xemf, sdehf
        double precision yhf, behf
        double precision cosa, rave
        integer i, j, k, nt
        dimension r(nt, 3), e(nt)
C
        save
C
C MH singlet & triplet
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d3_p, d2_p, d10_v, d10_b, d7_v, d5_b, d1_w, d7_
     *b, d2_v, d3_v
        double precision d6_b, d2_b, d3_b, d1_p, d4_v, d4_b, g_d1_w(g_pm
     *ax_), g_r(ldg_r, nt, 3), g_xenah(g_pmax_), g_s1(g_pmax_)
        double precision g_t1(g_pmax_), g_yhf(g_pmax_), g_behf(g_pmax_),
     * g_xehfa(g_pmax_), g_s2(g_pmax_), g_t2(g_pmax_), g_xemf(g_pmax_), 
     *g_t3(g_pmax_), g_s3(g_pmax_), g_coul1(g_pmax_)
        double precision g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_p
     *max_), g_exch2(g_pmax_), g_exch3(g_pmax_), g_w(g_pmax_), g_cplg2(g
     *_pmax_), g_e(ldg_e, nt)
        integer g_ehfid
        data denah /0.30d0/, renah /3.566044d0/, benah /2.1d0/
C
C HF singlet & triplet
        data rehf /1.733d0/, dehf /6.122d0/
C
C LiF singlet & triplet
        data demf /5.452d0/, remf /2.96d0/, bemf /0.85d0/
C
C
C H2 singlets
C
C misc constants
C
        data c2a /1.38661d0/, c2b /0.27362d0/, c2c /0.15d0/, um2p /0.76d
     *0/
C
        data g_ehfid /0/
c
c
        double precision e_diat,de_diat_dr,e_asym
        common /sur_diapot/ e_diat(3,3),de_diat_dr(3,3),e_asym(3)   

C
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        cevau = 1.d0 / 27.2113961d0
        onethd = 1.0d0 / 3.0d0
        twothd = 2.0d0 / 3.0d0

        e_asym(3) = (dehf+um2p)*cevau
C
        return
C
        entry g_pot22(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do 99999 i = 1, nt
C
          cosa = 0.5d0 * (r(i, 2) ** 2 + r(i, 3) ** 2 - r(i, 1) ** 2) / 
     *(r(i, 2) * r(i, 3))
          rave = 0.5d0 * (r(i, 1) + r(i, 3))
C
C =========================================================
C The diatomic curves . . .
C
C The MH singlet & triplet
C
          benah = 2.0d0
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-benah) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-benah) * (r(i, 1) - renah)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xenah(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xenah = d2_v
C--------
          d2_v = denah * xenah
          d3_v = xenah - 2.0d0
          d4_b = d2_v + d3_v * denah
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d4_b * g_xenah(g_i_)
          enddo
          s1 = d2_v * d3_v
C--------
          d2_v = 0.05d0 * denah * xenah
          d3_v = xenah + 2.0d0
          d4_b = d2_v + d3_v * (0.05d0 * denah)
          do g_i_ = 1, g_p_
            g_t1(g_i_) = d4_b * g_xenah(g_i_)
          enddo
          t1 = d2_v * d3_v

          e_diat(1,3) = t1*cevau
          de_diat_dr(1,3) = g_t1(1)*cevau

C--------
C
C The HF singlet and triplet
C
C
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
          d2_v = dehf * xehfa
          d3_v = xehfa - 2.0d0
          d4_b = d2_v + d3_v * dehf
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d4_b * g_xehfa(g_i_)
          enddo
          s2 = d2_v * d3_v

          e_diat(2,3) = s2*cevau
          de_diat_dr(2,3) = g_s2(2)*cevau

C--------
          d2_v = xehfa * xehfa
          d1_p = 2.0d0 * xehfa
          d2_b = 0.5d0 * dehf
          d5_b = d2_b * 1.1d0 + d2_b * d1_p
          do g_i_ = 1, g_p_
            g_t2(g_i_) = d5_b * g_xehfa(g_i_)
          enddo
          t2 = 0.5d0 * dehf * (d2_v + 1.1d0 * xehfa)
C--------
C
C The Modified HF singlet 
C
C The  LiF singlet & triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bemf) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bemf) * (r(i, 3) - remf)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xemf(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xemf = d2_v
C--------
          d2_v = xemf * xemf
          d1_p = 2.0d0 * xemf
          d2_b = 0.5d0 * (demf + um2p)
          d5_b = d2_b * 1.3d0 + d2_b * d1_p
          do g_i_ = 1, g_p_
            g_t3(g_i_) = d5_b * g_xemf(g_i_)
          enddo
          t3 = 0.5d0 * (demf + um2p) * (d2_v + 1.3d0 * xemf)
C--------
          d2_v = (demf + um2p) * xemf
          d3_v = xemf - 2.0d0
          d4_b = d2_v + d3_v * (demf + um2p)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d4_b * g_xemf(g_i_)
          enddo
          s3 = d2_v * d3_v

          e_diat(3,3) = s3*cevau
          de_diat_dr(3,3) = g_s3(3)*cevau

C--------
C
C ===============================================================
C The final LEPS form
C
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

          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d7_v)
          else
             d1_p = 0.0d0
          endif
          d7_b = (-(1.0d0 / dsqrt(2.0d0))) * d1_p
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = d7_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2
     *(g_i_) + g_coul1(g_i_)
          enddo
          e(i) = coul1 + coul2 + coul3 - d7_v / dsqrt(2.0d0) + um2p + de
     *hf
C--------
C
C
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
C
C
C***************************************************************************
C**************************************************************************
