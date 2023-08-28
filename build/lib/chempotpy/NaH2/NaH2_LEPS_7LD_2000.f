      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      integer :: nt
      double precision :: r(1,3), r2(3), v(1)
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

c   System:             NaH2
C   Functional form:    Modified Extended LEPS (2x2 diabatic fit)
c   Common name:        NaH2 surface set 7ld
c   Interface:          3-2V
c   Number of electronic surfaces: 2
c   Number of derivatives: 1

c   Notes:               This version contains long range forces,
c                        and also has derivatives.  This fit is based on the
c                        surface set 6 fit.

C   References:         Unpublished, M. D. Hack and D. G. Truhlar (2000).

C   Protocol:
C      PREPOT - initializes the potential's variables
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
C      Ground diabatic electronic state, first-excited diabatic electronic state
C and the potential coupling.
C
C   Zero of energy:
C      The classical potential energy is set equal to zero for the Na
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Parameters:
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Na-H)
C                            R(2) = R(H-H')
C                            R(3) = R(Na-H')

c
c This is NaH2 potential matrix 7ld, finalized on 5/17/99 
c This version contains analytic derivatives.

c
c  The potential is called by first issuing a call to PREPOT which initializes
c  all surfaces.  The potential is called by calling POT with the parameters
c  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of
c  dimension nt which is an integer representing the number of geometries
c  input.  The integer nsurf selects the desired surface.  nsurf=1 is the
c  lower diabatic surface, nsurf=2 is the coupling surface, nsurf=3 is
c  the upper diabatic surface, nsurf = 4 selects the lower adiabatic
c  surface, and nsurf = 5 selects the upper adiabatic surface 
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
          do i=1,nt
           call 
     &      g_pot11(g_p_, r(i,1), g_r(1,i,1), ldg_r, v11, De1, ldg_e, 1)
           call 
     &      g_pot12(g_p_, r(i,1), g_r(1,i,1), ldg_r, v12, De2, ldg_e, 1)
           call 
     &      g_pot22(g_p_, r(i,1), g_r(1,i,1), ldg_r, v22, De3, ldg_e, 1)
            dis = dsqrt((v11-v22)**2+4.d0*v12**2)
            do j = 1,3
              ddis(j) = (v11-v22)*(de1(j,1)-de3(j,1))+4.0d0*v12*de2(j,1)
              ddis(j) = ddis(j)/dis
            enddo
            if (nsurf .eq. 4) then
              e(i)=0.5d0*(v11+v22-dis)
              do j = 1,3
                de(j,i) = 0.5d0*(de1(j,1)+de3(j,1)-ddis(j))
              enddo
            else
              e(i)=0.5d0*(v11+v22+dis)
              do j = 1,3
                de(j,i) = 0.5d0*(de1(j,1)+de3(j,1)+ddis(j))
              enddo
            end if
          enddo
        end if

      return

      end
C ********************************************************************
C ********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 05/17/99 by the version of
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
        implicit double precision (a-h, o-z)
C
        dimension r(nt, 3), e(nt)
C
        save
C NaH singlet
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d24_b, d13_b, d20_v, d2_p, d15_v, d16_v, d15_b,
     * d2_v, d22_v, d2_b
        double precision d1_w, d2_w, d3_v, d4_v, d5_v, d6_v, d3_b, d4_b,
     * d5_b, d6_b
        double precision d1_p, d7_v, d8_v, d9_v, d10_v, d7_b, d8_b, d9_b
     *, d10_b, d11_v
        double precision d12_v, d13_v, d19_v, d11_b, d12_b, g_r1s(g_pmax
     *_), g_r(ldg_r, nt, 3), g_r2s(g_pmax_), g_r3s(g_pmax_), g_d2_w(g_pm
     *ax_)
        double precision g_d1_w(g_pmax_), g_rnag(g_pmax_), g_rnags(g_pma
     *x_), g_rnagc(g_pmax_), g_rmid(g_pmax_), g_rave(g_pmax_), g_cos1(g_
     *pmax_), g_cos2(g_pmax_), g_cosg(g_pmax_), g_rnag1(g_pmax_)
        double precision g_rnag3(g_pmax_), g_rnag1c(g_pmax_), g_rnag3c(g
     *_pmax_), g_cstsq(g_pmax_), g_cosa(g_pmax_), g_cstsq1(g_pmax_), g_c
     *stsq3(g_pmax_), g_xes1(g_pmax_), g_s1(g_pmax_), g_beh2(g_pmax_)
        double precision g_sswt(g_pmax_), g_deh2m(g_pmax_), g_xeh2(g_pma
     *x_), g_s2(g_pmax_), g_xes3(g_pmax_), g_s3(g_pmax_), g_h2ie(g_pmax_
     *), g_al1(g_pmax_), g_al2(g_pmax_), g_rnahie1(g_pmax_)
        double precision g_rnahie3(g_pmax_), g_rnahmu1(g_pmax_), g_rnahm
     *u3(g_pmax_), g_rnahpar1(g_pmax_), g_rnahpar3(g_pmax_), g_rnahperp1
     *(g_pmax_), g_rnahperp3(g_pmax_), g_rnahpol1(g_pmax_), g_rnahpol1b(
     *g_pmax_), g_rnahpol3(g_pmax_)
        double precision g_rnahpol3b(g_pmax_), g_div1(g_pmax_), g_rmuimu
     *1(g_pmax_), g_div3(g_pmax_), g_rmuimu3(g_pmax_), g_c690(g_pmax_), 
     *g_c60(g_pmax_), g_c6(g_pmax_), g_c61(g_pmax_), g_dexp1(g_pmax_)
        double precision g_c63(g_pmax_), g_dexp3(g_pmax_), g_cind1(g_pma
     *x_), g_cind3(g_pmax_), g_c6swt1(g_pmax_), g_c6swt3(g_pmax_), g_ele
     *c(g_pmax_), g_e(ldg_e, nt), g_xect(g_pmax_), g_vmid(g_pmax_)
        integer g_ehfid
        intrinsic dble
        data des1 /0.55118d0/, bes1 /1.05118d0/, res1 /3.38583d0/
C
C H2 singlet
        data reh2 /1.401d0/, bf /1.627d0/, b0 /1.194d0/, gamh2a /2.323d0
     */, gamh2b /3.126d0/, deh2 /4.74806d0/
        data des20 /0.996d0/, s2alp /3.0d0/, s2rho /7.0d0/
C
C H-Na-H term
        data dect /2.24409d0/, bect /0.57874d0/, rect /3.46457d0/, bect2
     * /0.87143d0/
C
C data blocks for the long range forces
        data h2a1 /2.38029900d0/, h2a2 /1.47400334d0/, h2a3 /-0.22633177
     *1d0/, h2a4 /-2.17339562d0/, h2a5 /9.54597998d0/, h2a6 /-23.3366749
     *d0/, h2a7 /24.7115748d0/, h2a8 /-12.8786102d0/, h2a9 /2.75711758d0
     */
        data h2b1 /1.35164094d0/, h2b2 /-0.254404564d0/, h2b3 /-10.66270
     *55d0/, h2b4 /-14.0112406d0/, h2b5 /-7.80660398d0/, h2b6 /5.4769775
     *3d0/, h2b7 /-6.95875664d0/, h2b8 /3.64334771d0/, h2b9 /-0.49250559
     *6d0/
        data h2i1 /160.484d0/, h2i2 /-20.397d0/, h2i3 /3.913d0/, h2i4 /0
     *.777d0/, h2i5 /-7.583/, h2i6 /1.501/
        data rnaie /5.139d0/, rnapol /159.267d0/, hpol /4.5d0/, hie /13.
     *598d0/
        data hnai1 /-193.341d0/, hnai2 /16.968d0/, hnai3 /-77.120d0/, hn
     *ai4 /-1.024d0/, hnai5 /6.009d0/, hnai6 /-1.241d0/
        data hnau1 /119.880d0/, hnau2 /-0.620d0/, hnau3 /31.335d0/, hnau
     *4 /70.414d0/, hnau5 /-13.016d0/, hnau6 /0.8681d0/
        data hnaa1 /0.67d0/, hnaa2 /-92.229d0/, hnaa3 /51.01019d0/, hnaa
     *4 /-130.37826d0/, hnaa5 /6.40285d0/
        data hnab1 /0.34333d0/, hnab2 /-92.229d0/, hnab3 /-48.54867d0/, 
     *hnab4 /-17.30523d0/, hnab5 /0.99821d0/
        data balp /2.07234d0/, brho /2.0d0/, rb0 /5.1998697814229d7/, rb
     *90 /5.0d8/, f90 /1.0d0/, f0 /1.0d0/, rd /9.1538750508008d8/, dalp 
     */1.43307d0/, drho /5.71654d0/
C
C misc constants
C
        data g_ehfid /0/
C
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        rac2 = 1.d0 / dsqrt(2.d0)
        cevau = 1.d0 / 27.2113961d0
        onethd = 1.0d0 / 3.0d0
        twothd = 2.0d0 / 3.0d0
C
C
        return
C
        entry g_pot11(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do 99999 i = 1, nt
C
C
C ============================================================
C useful quantities calculated first, r1^2, R2^2, R3^2, atom to diatom
C center of mass distances, the H-Na-H bond angle and the Jacobian
C angle chi. 
C
          d2_b = r(i, 1) + r(i, 1)
          do g_i_ = 1, g_p_
            g_r1s(g_i_) = d2_b * g_r(g_i_, i, 1)
          enddo
          r1s = r(i, 1) * r(i, 1)
C--------
          d2_b = r(i, 2) + r(i, 2)
          do g_i_ = 1, g_p_
            g_r2s(g_i_) = d2_b * g_r(g_i_, i, 2)
          enddo
          r2s = r(i, 2) * r(i, 2)
C--------
          d2_b = r(i, 3) + r(i, 3)
          do g_i_ = 1, g_p_
            g_r3s(g_i_) = d2_b * g_r(g_i_, i, 3)
          enddo
          r3s = r(i, 3) * r(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-0.5d0) * g_r2s(g_i_) + g_r3s(g_i_) + g_r1s(
     *g_i_)
          enddo
          d2_w = r1s + r3s - 0.5d0 * r2s
          d2_v = abs(d2_w)

          if (d2_w .gt. 0.0d0) then
             d1_p =  1.0d0
          else if (d2_w .lt. 0.0d0) then
             d1_p = -1.0d0
          else
             d1_p = 0.0d0
          endif
          d3_b = 0.5d0 * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_d2_w(g_i_)
          enddo
          d1_w = 0.5d0 * d2_v
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
          d2_b = rnag + rnag
          do g_i_ = 1, g_p_
            g_rnags(g_i_) = d2_b * g_rnag(g_i_)
          enddo
          rnags = rnag * rnag
C--------
          do g_i_ = 1, g_p_
            g_rnagc(g_i_) = rnags * g_rnag(g_i_) + rnag * g_rnags(g_i_)
          enddo
          rnagc = rnags * rnag
C--------
C
          do g_i_ = 1, g_p_
            g_rmid(g_i_) = 0.5d0 * g_r(g_i_, i, 3) + 0.5d0 * g_r(g_i_, i
     *, 1)
          enddo
          rmid = 0.5d0 * (r(i, 1) + r(i, 3))
C--------
          do g_i_ = 1, g_p_
            g_rave(g_i_) = (-0.5d0) * g_r(g_i_, i, 3) + 0.5d0 * g_r(g_i_
     *, i, 1)
          enddo
          rave = 0.5d0 * (r(i, 1) - r(i, 3))
C--------
C
C cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
C
          d7_v = 2.0d0 * r(i, 1)
          d9_v = d7_v * r(i, 2)
          d10_v = (r1s + r2s - r3s) / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = (-d10_v) / d9_v
          d5_b = d3_b * d7_v
          d6_b = d3_b * r(i, 2) * 2.0d0
          do g_i_ = 1, g_p_
            g_cos1(g_i_) = d5_b * g_r(g_i_, i, 2) + d6_b * g_r(g_i_, i, 
     *1) + (-d2_b) * g_r3s(g_i_) + d2_b * g_r2s(g_i_) + d2_b * g_r1s(g_i
     *_)
          enddo
          cos1 = d10_v
C--------
          d7_v = 2.0d0 * r(i, 3)
          d9_v = d7_v * r(i, 2)
          d10_v = (r3s + r2s - r1s) / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = (-d10_v) / d9_v
          d5_b = d3_b * d7_v
          d6_b = d3_b * r(i, 2) * 2.0d0
          do g_i_ = 1, g_p_
            g_cos2(g_i_) = d5_b * g_r(g_i_, i, 2) + d6_b * g_r(g_i_, i, 
     *3) + (-d2_b) * g_r1s(g_i_) + d2_b * g_r2s(g_i_) + d2_b * g_r3s(g_i
     *_)
          enddo
          cos2 = d10_v
C--------
          do g_i_ = 1, g_p_
            g_cosg(g_i_) = (-0.5d0) * g_cos2(g_i_) + 0.5d0 * g_cos1(g_i_
     *)
          enddo
          cosg = 0.5d0 * (cos1 - cos2)
C--------
C
          frac = 0.958d0
C H to NaH center of masses
          d5_b = frac * frac + (-frac)
          d6_b = 1.0d0 + (-frac)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = d5_b * g_r1s(g_i_) + frac * g_r3s(g_i_) + d6_
     *b * g_r2s(g_i_)
          enddo
          d2_w = r2s + frac * (r3s - r2s - r1s) + r1s * frac * frac
          d2_v = abs(d2_w)

          if (d2_w .gt. 0.0d0) then
             d1_p =  1.0d0
          else if (d2_w .lt. 0.0d0) then
             d1_p = -1.0d0
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d1_p * g_d2_w(g_i_)
          enddo
          d1_w = d2_v
          d2_v = sqrt(d1_w)

          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_rnag1(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          rnag1 = d2_v
C--------
          d5_b = frac * frac + (-frac)
          d6_b = 1.0d0 + (-frac)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = d5_b * g_r3s(g_i_) + frac * g_r1s(g_i_) + d6_
     *b * g_r2s(g_i_)
          enddo
          d2_w = r2s + frac * (r1s - r2s - r3s) + r3s * frac * frac
          d2_v = abs(d2_w)

          if (d2_w .gt. 0.0d0) then
             d1_p =  1.0d0
          else if (d2_w .lt. 0.0d0) then
             d1_p = -1.0d0
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d1_p * g_d2_w(g_i_)
          enddo
          d1_w = d2_v
          d2_v = sqrt(d1_w)

          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_rnag3(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          rnag3 = d2_v
C--------
C
          d2_v = rnag1 * rnag1
          d3_b = d2_v + rnag1 * rnag1 + rnag1 * rnag1
          do g_i_ = 1, g_p_
            g_rnag1c(g_i_) = d3_b * g_rnag1(g_i_)
          enddo
          rnag1c = d2_v * rnag1
C--------
          d2_v = rnag3 * rnag3
          d3_b = d2_v + rnag3 * rnag3 + rnag3 * rnag3
          do g_i_ = 1, g_p_
            g_rnag3c(g_i_) = d3_b * g_rnag3(g_i_)
          enddo
          rnag3c = d2_v * rnag3
C--------
C
C the shortest of any two distances . . .
          rsm1 = r(i, 3) + (r(i, 2) - r(i, 3)) * (0.5d0 + 0.5d0 * dtanh(
     *0.5d0 * (r(i, 3) - r(i, 2))))
C
          rsm2 = r(i, 3) + (r(i, 1) - r(i, 3)) * (0.5d0 + 0.5d0 * dtanh(
     *0.5d0 * (r(i, 3) - r(i, 1))))
C
          rsm3 = r(i, 1) + (r(i, 2) - r(i, 1)) * (0.5d0 + 0.5d0 * dtanh(
     *0.5d0 * (r(i, 1) - r(i, 2))))
C
          rexp1 = dexp(-r(i, 1))
          rexp2 = dexp(-r(i, 2))
          rexp3 = dexp(-r(i, 3))
C
C checks on angles . . .
          if (rnag .eq. 0.0d0) then
            do g_i_ = 1, g_p_
              g_cstsq(g_i_) = 0.0d0
            enddo
            cstsq = 1.0d0
C--------
          else
            d7_v = r(i, 2) * rnag
            d8_v = 0.5d0 * (r3s - r1s) / d7_v
            d9_v = d8_v * d8_v
            d1_p = 2.0d0 * d8_v
            d4_b = d1_p * ((-d8_v) / d7_v)
            d5_b = d4_b * rnag
            d6_b = d4_b * r(i, 2)
            d7_b = d1_p * (1.0d0 / d7_v) * 0.5d0
            do g_i_ = 1, g_p_
              g_cstsq(g_i_) = d6_b * g_rnag(g_i_) + d5_b * g_r(g_i_, i, 
     *2) + (-d7_b) * g_r1s(g_i_) + d7_b * g_r3s(g_i_)
            enddo
            cstsq = d9_v
C--------
          endif
C
          d7_v = dble(2) * r(i, 1)
          d9_v = d7_v * r(i, 3)
          d10_v = (r1s + r3s - r2s) / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = (-d10_v) / d9_v
          d5_b = d3_b * d7_v
          d6_b = d3_b * r(i, 3) * dble(2)
          do g_i_ = 1, g_p_
            g_cosa(g_i_) = d5_b * g_r(g_i_, i, 3) + d6_b * g_r(g_i_, i, 
     *1) + (-d2_b) * g_r2s(g_i_) + d2_b * g_r3s(g_i_) + d2_b * g_r1s(g_i
     *_)
          enddo
          cosa = d10_v
C--------
C
          coscos = cstsq * (0.5d0 + 0.5d0 * cosa) * (1.0d0 - dexp((-10.0
     *d0) * rnags))
C
          if (rnag1 .eq. 0.0d0) then
            do g_i_ = 1, g_p_
              g_cstsq1(g_i_) = 0.0d0
            enddo
            cstsq1 = 1.0d0
C--------
          else
            d11_v = r(i, 1) * rnag1
            d12_v = (r3s - r2s - r1s + 2.0d0 * r3s * frac) / d11_v
            d13_v = d12_v * d12_v
            d1_p = 2.0d0 * d12_v
            d3_b = 0.25d0 * d1_p
            d4_b = d3_b * (1.0d0 / d11_v)
            d5_b = d3_b * ((-d12_v) / d11_v)
            d6_b = d5_b * rnag1
            d7_b = d5_b * r(i, 1)
            d11_b = d4_b * frac * 2.0d0 + d4_b
            do g_i_ = 1, g_p_
              g_cstsq1(g_i_) = d7_b * g_rnag1(g_i_) + d6_b * g_r(g_i_, i
     *, 1) + (-d4_b) * g_r1s(g_i_) + (-d4_b) * g_r2s(g_i_) + d11_b * g_r
     *3s(g_i_)
            enddo
            cstsq1 = 0.25d0 * d13_v
C--------
          endif
C
          if (rnag3 .eq. 0.0d0) then
            do g_i_ = 1, g_p_
              g_cstsq3(g_i_) = 0.0d0
            enddo
            cstsq3 = 1.0d0
C--------
          else
            d11_v = r(i, 3) * rnag3
            d12_v = (r1s - r2s - r3s + 2.0d0 * r1s * frac) / d11_v
            d13_v = d12_v * d12_v
            d1_p = 2.0d0 * d12_v
            d3_b = 0.25d0 * d1_p
            d4_b = d3_b * (1.0d0 / d11_v)
            d5_b = d3_b * ((-d12_v) / d11_v)
            d6_b = d5_b * rnag3
            d7_b = d5_b * r(i, 3)
            d11_b = d4_b * frac * 2.0d0 + d4_b
            do g_i_ = 1, g_p_
              g_cstsq3(g_i_) = d7_b * g_rnag3(g_i_) + d6_b * g_r(g_i_, i
     *, 3) + (-d4_b) * g_r3s(g_i_) + (-d4_b) * g_r2s(g_i_) + d11_b * g_r
     *1s(g_i_)
            enddo
            cstsq3 = 0.25d0 * d13_v
C--------
          endif
C
C==============================================================k
CThis section contains the diatomic curves
C
C The first NaH singlet
C
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
          d2_v = des1 * onethd * xes1
          d3_v = xes1 + 2.0d0
          d4_b = d2_v + d3_v * (des1 * onethd)
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d4_b * g_xes1(g_i_)
          enddo
          s1 = d2_v * d3_v
C--------
C
C The H2 singlet
C
C new singlet, designed to minic lower adiabatic behavior
          reh2 = 1.401d0
          bf = 1.627d0
          b0 = 1.194d0
          gamh2a = 2.323d0
          gamh2b = 3.126d0
C
          d5_v = (r(i, 2) / gamh2b) * (r(i, 2) / gamh2b)
          d2_p = 2.0d0 * (r(i, 2) / gamh2b)
          d11_v = (r(i, 2) / gamh2b) * (r(i, 2) / gamh2b)
          d1_p = 2.0d0 * (r(i, 2) / gamh2b)
          d12_v = bf - r(i, 2) / gamh2a + d11_v
          d13_v = bf * (b0 - r(i, 2) / gamh2a + d5_v) / d12_v
          d3_b = (-d13_v) / d12_v
          d9_b = 1.0d0 / d12_v * bf
          d7_b = d3_b * d1_p * (1.0d0 / gamh2b) + (-d3_b) * (1.0d0 / gam
     *h2a) + d9_b * d2_p * (1.0d0 / gamh2b) + (-d9_b) * (1.0d0 / gamh2a)
          do g_i_ = 1, g_p_
            g_beh2(g_i_) = d7_b * g_r(g_i_, i, 2)
          enddo
          beh2 = d13_v
C--------
C
          d3_b = 1.0d0 / s2alp
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_rmid(g_i_)
          enddo
          d1_w = (rmid - s2rho) / s2alp
          d2_v = tanh (d1_w)
          d2_p = 1.0d0 - ( d2_v *  d2_v)
          d4_v = 0.5d0 * (1.0d0 - d2_v)
          d6_v = cosg * cosg
          d1_p = 2.0d0 * cosg
          d4_b = d4_v * d1_p
          d7_b = (-(d6_v * 0.5d0)) * d2_p
          do g_i_ = 1, g_p_
            g_sswt(g_i_) = d4_b * g_cosg(g_i_) + d7_b * g_d1_w(g_i_)
          enddo
          sswt = d4_v * d6_v
C--------
          d3_b = deh2 * (des20 - 1.0d0)
          do g_i_ = 1, g_p_
            g_deh2m(g_i_) = d3_b * g_sswt(g_i_)
          enddo
          deh2m = deh2 + deh2 * (des20 - 1.0d0) * sswt
C--------
C
          d4_v = r(i, 2) - reh2
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beh2) * g_r(g_i_, i, 2) + (-d4_v) * g_beh2(
     *g_i_)
          enddo
          d1_w = (-beh2) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xeh2(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xeh2 = d2_v
C--------
          d3_v = deh2m * xeh2
          d4_v = xeh2 - 2.0d0
          d5_b = d4_v * xeh2
          d4_b = d3_v + d4_v * deh2m
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d4_b * g_xeh2(g_i_) + d5_b * g_deh2m(g_i_)
          enddo
          s2 = d3_v * d4_v
C--------
C
C The second NaH singlet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bes1) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bes1) * (r(i, 3) - res1)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xes3(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xes3 = d2_v
C--------
          d2_v = des1 * onethd * xes3
          d3_v = xes3 + 2.0d0
          d4_b = d2_v + d3_v * (des1 * onethd)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d4_b * g_xes3(g_i_)
          enddo
          s3 = d2_v * d3_v
C--------
C
C
C ===============================================================
C this long section contains data pertinent to long range forces
C
C
C  H2 ionization energy:
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-h2i3) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-h2i3) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-h2i6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-h2i6) * r(i, 2)
          d3_v = h2i1 + h2i2 * r(i, 2)
          d5_v = exp(d1_w)
          d2_p =  d5_v
          d9_v = h2i4 + h2i5 * r(i, 2)
          d10_v = r2s * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d8_b = d10_v * d1_p
          d9_b = d12_v * d9_v
          d15_b = d3_v * d2_p
          d12_b = d12_v * r2s * h2i5 + d5_v * h2i2
          do g_i_ = 1, g_p_
            g_h2ie(g_i_) = d8_b * g_d2_w(g_i_) + d9_b * g_r2s(g_i_) + d1
     *5_b * g_d1_w(g_i_) + d12_b * g_r(g_i_, i, 2)
          enddo
          h2ie = d3_v * d5_v + d10_v * d12_v - s2a + hie
C--------
C
C  H2 polarizability:
C
C   al1 is the parallel polarizability for h2
C   al2 is the perpendicaular polarizability for h2
C   alt is the anisotropic? polarizability
C   alb is the mean polarizability
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-h2a1) * g_r(g_i_, i, 2)
          enddo
          d1_w = h2a1 * (h2a2 - r(i, 2))
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = h2a6 * r(i, 2)
          d12_v = h2a7 * r2s
          d15_v = h2a8 * r2s
          d16_v = d15_v * r2s
          d19_v = h2a9 * r2s
          d20_v = d19_v * r2s
          d22_v = h2a3 + h2a4 * r(i, 2) + h2a5 * r2s + d9_v * r2s + d12_
     *v * r2s + d16_v * r(i, 2) + d20_v * r2s
          d7_b = d2_v * r2s
          d12_b = d2_v * r(i, 2)
          d8_b = d2_v * d20_v + d7_b * d19_v + d7_b * r2s * h2a9 + d12_b
     * * d15_v + d12_b * r2s * h2a8 + d2_v * d12_v + d2_v * r2s * h2a7 +
     * d2_v * d9_v + d2_v * h2a5
          d13_b = d2_v * d16_v + d2_v * r2s * h2a6 + d2_v * h2a4
          d24_b = d22_v * d1_p
          do g_i_ = 1, g_p_
            g_al1(g_i_) = d8_b * g_r2s(g_i_) + d13_b * g_r(g_i_, i, 2) +
     * d24_b * g_d1_w(g_i_)
          enddo
          al1 = d2_v * d22_v + 2.0d0 * hpol
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-h2b1) * g_r(g_i_, i, 2)
          enddo
          d1_w = h2b1 * (h2b2 - r(i, 2))
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = h2b6 * r(i, 2)
          d12_v = h2b7 * r2s
          d15_v = h2b8 * r2s
          d16_v = d15_v * r2s
          d19_v = h2b9 * r2s
          d20_v = d19_v * r2s
          d22_v = h2b3 + h2b4 * r(i, 2) + h2b5 * r2s + d9_v * r2s + d12_
     *v * r2s + d16_v * r(i, 2) + d20_v * r2s
          d7_b = d2_v * r2s
          d12_b = d2_v * r(i, 2)
          d8_b = d2_v * d20_v + d7_b * d19_v + d7_b * r2s * h2b9 + d12_b
     * * d15_v + d12_b * r2s * h2b8 + d2_v * d12_v + d2_v * r2s * h2b7 +
     * d2_v * d9_v + d2_v * h2b5
          d13_b = d2_v * d16_v + d2_v * r2s * h2b6 + d2_v * h2b4
          d24_b = d22_v * d1_p
          do g_i_ = 1, g_p_
            g_al2(g_i_) = d8_b * g_r2s(g_i_) + d13_b * g_r(g_i_, i, 2) +
     * d24_b * g_d1_w(g_i_)
          enddo
          al2 = d2_v * d22_v + 2.0d0 * hpol
C--------
C
          alt = twothd * (al1 - al2)
          alb = onethd * al1 + twothd * al2
C
          alhh = alb + alt * (1.5d0 * cstsq - 0.5d0)
C
C
C  The NaH + H dispersion and NaH + H induction terms must go on the lower diaba
Ct,
C  so that the adiabat is fitted correctly in the region of the coupling
C
C NaH properties
C
C
          d9_v = hnai6 * r(i, 1)
          d11_v = hnai3 + hnai4 * r(i, 1) + hnai5 * r1s + d9_v * r1s
          d12_v = (hnai1 + hnai2 * r(i, 1)) / d11_v
          d4_b = (-d12_v) / d11_v
          d8_b = d4_b * d9_v + d4_b * hnai5
          d9_b = d4_b * r1s * hnai6 + d4_b * hnai4 + 1.0d0 / d11_v * hna
     *i2
          do g_i_ = 1, g_p_
            g_rnahie1(g_i_) = d8_b * g_r1s(g_i_) + d9_b * g_r(g_i_, i, 1
     *)
          enddo
          rnahie1 = d12_v + rnaie
C--------
C
          d9_v = hnai6 * r(i, 3)
          d11_v = hnai3 + hnai4 * r(i, 3) + hnai5 * r3s + d9_v * r3s
          d12_v = (hnai1 + hnai2 * r(i, 3)) / d11_v
          d4_b = (-d12_v) / d11_v
          d8_b = d4_b * d9_v + d4_b * hnai5
          d9_b = d4_b * r3s * hnai6 + d4_b * hnai4 + 1.0d0 / d11_v * hna
     *i2
          do g_i_ = 1, g_p_
            g_rnahie3(g_i_) = d8_b * g_r3s(g_i_) + d9_b * g_r(g_i_, i, 3
     *)
          enddo
          rnahie3 = d12_v + rnaie
C--------
C
          d10_v = hnau6 * r(i, 1)
          d12_v = hnau3 + hnau4 * r(i, 1) + hnau5 * r1s + d10_v * r1s
          d13_v = (hnau1 * r(i, 1) + hnau2 * r1s) / d12_v
          d2_b = 1.0d0 / d12_v
          d3_b = (-d13_v) / d12_v
          d7_b = d3_b * d10_v + d3_b * hnau5 + d2_b * hnau2
          d8_b = d3_b * r1s * hnau6 + d3_b * hnau4 + d2_b * hnau1
          do g_i_ = 1, g_p_
            g_rnahmu1(g_i_) = d7_b * g_r1s(g_i_) + d8_b * g_r(g_i_, i, 1
     *)
          enddo
          rnahmu1 = d13_v
C--------
C
          d10_v = hnau6 * r(i, 3)
          d12_v = hnau3 + hnau4 * r(i, 3) + hnau5 * r3s + d10_v * r3s
          d13_v = (hnau1 * r(i, 3) + hnau2 * r3s) / d12_v
          d2_b = 1.0d0 / d12_v
          d3_b = (-d13_v) / d12_v
          d7_b = d3_b * d10_v + d3_b * hnau5 + d2_b * hnau2
          d8_b = d3_b * r3s * hnau6 + d3_b * hnau4 + d2_b * hnau1
          do g_i_ = 1, g_p_
            g_rnahmu3(g_i_) = d7_b * g_r3s(g_i_) + d8_b * g_r(g_i_, i, 3
     *)
          enddo
          rnahmu3 = d13_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-hnaa1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-hnaa1) * r(i, 1)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = hnaa5 * r1s
          d11_v = hnaa2 + hnaa3 * r(i, 1) + hnaa4 * r1s + d9_v * r(i, 1)
          d9_b = d2_v * r(i, 1) * hnaa5 + d2_v * hnaa4
          d8_b = d2_v * d9_v + d2_v * hnaa3
          d13_b = d11_v * d1_p
          do g_i_ = 1, g_p_
            g_rnahpar1(g_i_) = d9_b * g_r1s(g_i_) + d8_b * g_r(g_i_, i, 
     *1) + d13_b * g_d1_w(g_i_)
          enddo
          rnahpar1 = rnapol + hpol + d2_v * d11_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-hnaa1) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-hnaa1) * r(i, 3)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = hnaa5 * r3s
          d11_v = hnaa2 + hnaa3 * r(i, 3) + hnaa4 * r3s + d9_v * r(i, 3)
          d9_b = d2_v * r(i, 3) * hnaa5 + d2_v * hnaa4
          d8_b = d2_v * d9_v + d2_v * hnaa3
          d13_b = d11_v * d1_p
          do g_i_ = 1, g_p_
            g_rnahpar3(g_i_) = d9_b * g_r3s(g_i_) + d8_b * g_r(g_i_, i, 
     *3) + d13_b * g_d1_w(g_i_)
          enddo
          rnahpar3 = rnapol + hpol + d2_v * d11_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-hnab1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-hnab1) * r(i, 1)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = hnab5 * r1s
          d11_v = hnab2 + hnab3 * r(i, 1) + hnab4 * r1s + d9_v * r(i, 1)
          d9_b = d2_v * r(i, 1) * hnab5 + d2_v * hnab4
          d8_b = d2_v * d9_v + d2_v * hnab3
          d13_b = d11_v * d1_p
          do g_i_ = 1, g_p_
            g_rnahperp1(g_i_) = d9_b * g_r1s(g_i_) + d8_b * g_r(g_i_, i,
     * 1) + d13_b * g_d1_w(g_i_)
          enddo
          rnahperp1 = rnapol + hpol + d2_v * d11_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-hnab1) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-hnab1) * r(i, 3)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = hnab5 * r3s
          d11_v = hnab2 + hnab3 * r(i, 3) + hnab4 * r3s + d9_v * r(i, 3)
          d9_b = d2_v * r(i, 3) * hnab5 + d2_v * hnab4
          d8_b = d2_v * d9_v + d2_v * hnab3
          d13_b = d11_v * d1_p
          do g_i_ = 1, g_p_
            g_rnahperp3(g_i_) = d9_b * g_r3s(g_i_) + d8_b * g_r(g_i_, i,
     * 3) + d13_b * g_d1_w(g_i_)
          enddo
          rnahperp3 = rnapol + hpol + d2_v * d11_v
C--------
C
          d3_v = rnahpar1 - rnahperp1
          d2_b = 1.0d0 + (-cstsq1)
          do g_i_ = 1, g_p_
            g_rnahpol1(g_i_) = d3_v * g_cstsq1(g_i_) + cstsq1 * g_rnahpa
     *r1(g_i_) + d2_b * g_rnahperp1(g_i_)
          enddo
          rnahpol1 = rnahperp1 + d3_v * cstsq1
C--------
          do g_i_ = 1, g_p_
            g_rnahpol1b(g_i_) = twothd * g_rnahperp1(g_i_) + onethd * g_
     *rnahpar1(g_i_)
          enddo
          rnahpol1b = onethd * rnahpar1 + twothd * rnahperp1
C--------
C
          d3_v = rnahpar3 - rnahperp3
          d2_b = 1.0d0 + (-cstsq3)
          do g_i_ = 1, g_p_
            g_rnahpol3(g_i_) = d3_v * g_cstsq3(g_i_) + cstsq3 * g_rnahpa
     *r3(g_i_) + d2_b * g_rnahperp3(g_i_)
          enddo
          rnahpol3 = rnahperp3 + d3_v * cstsq3
C--------
          do g_i_ = 1, g_p_
            g_rnahpol3b(g_i_) = twothd * g_rnahperp3(g_i_) + onethd * g_
     *rnahpar3(g_i_)
          enddo
          rnahpol3b = onethd * rnahpar3 + twothd * rnahperp3
C--------
C
C induction forces
C dipole induced dipole
C
          d2_v = rnag1c * rnag1c
          d4_v = d2_v * rnag1c + rd
          d5_v = 1.0d0 / d4_v
          d3_b = (-d5_v) / d4_v
          d4_b = d3_b * rnag1c
          d5_b = d3_b * d2_v + d4_b * rnag1c + d4_b * rnag1c
          do g_i_ = 1, g_p_
            g_div1(g_i_) = d5_b * g_rnag1c(g_i_)
          enddo
          div1 = d5_v
C--------
          d2_v = rnahmu1 * rnahmu1
          d1_p = 2.0d0 * rnahmu1
          d4_v = (-d2_v) * hpol
          d7_v = 1.5d0 * cstsq1 + 0.5d0
          d8_v = d4_v * d7_v
          d10_v = d8_v * rnag1c
          d4_b = div1 * rnag1c
          d5_b = div1 * d8_v
          d9_b = d4_b * d4_v * 1.5d0
          d12_b = (-(d4_b * d7_v * hpol)) * d1_p
          do g_i_ = 1, g_p_
            g_rmuimu1(g_i_) = d10_v * g_div1(g_i_) + d5_b * g_rnag1c(g_i
     *_) + d9_b * g_cstsq1(g_i_) + d12_b * g_rnahmu1(g_i_)
          enddo
          rmuimu1 = d10_v * div1
C--------
C
          d2_v = rnag3c * rnag3c
          d4_v = d2_v * rnag3c + rd
          d5_v = 1.0d0 / d4_v
          d3_b = (-d5_v) / d4_v
          d4_b = d3_b * rnag3c
          d5_b = d3_b * d2_v + d4_b * rnag3c + d4_b * rnag3c
          do g_i_ = 1, g_p_
            g_div3(g_i_) = d5_b * g_rnag3c(g_i_)
          enddo
          div3 = d5_v
C--------
          d2_v = rnahmu3 * rnahmu3
          d1_p = 2.0d0 * rnahmu3
          d4_v = (-d2_v) * hpol
          d7_v = 1.5d0 * cstsq3 + 0.5d0
          d8_v = d4_v * d7_v
          d10_v = d8_v * rnag3c
          d4_b = div3 * rnag3c
          d5_b = div3 * d8_v
          d9_b = d4_b * d4_v * 1.5d0
          d12_b = (-(d4_b * d7_v * hpol)) * d1_p
          do g_i_ = 1, g_p_
            g_rmuimu3(g_i_) = d10_v * g_div3(g_i_) + d5_b * g_rnag3c(g_i
     *_) + d9_b * g_cstsq3(g_i_) + d12_b * g_rnahmu3(g_i_)
          enddo
          rmuimu3 = d10_v * div3
C--------
C
C dispersion forces
          d3_v = (-1.5d0) * rnaie * h2ie * rnapol
          d6_v = rnaie + h2ie
          d8_v = rnagc * rnagc
          d10_v = d8_v * rnagc + rb90
          d11_v = d6_v * d10_v
          d12_v = d3_v * al2 / d11_v
          d2_b = 1.0d0 / d11_v
          d3_b = (-d12_v) / d11_v
          d6_b = d3_b * d6_v
          d7_b = d6_b * rnagc
          d8_b = d6_b * d8_v + d7_b * rnagc + d7_b * rnagc
          d11_b = d2_b * d3_v
          d9_b = d3_b * d10_v + d2_b * al2 * rnapol * ((-1.5d0) * rnaie)
          do g_i_ = 1, g_p_
            g_c690(g_i_) = d8_b * g_rnagc(g_i_) + d11_b * g_al2(g_i_) + 
     *d9_b * g_h2ie(g_i_)
          enddo
          c690 = d12_v
C--------
          d2_v = c690 * f90
          d4_b = rnagc * f90
          do g_i_ = 1, g_p_
            g_c690(g_i_) = d2_v * g_rnagc(g_i_) + d4_b * g_c690(g_i_)
          enddo
          c690 = d2_v * rnagc
C--------
C
          d3_v = (-1.5d0) * rnaie * h2ie * rnapol
          d6_v = rnaie + h2ie
          d8_v = rnagc * rnagc
          d10_v = d8_v * rnagc + rb0
          d11_v = d6_v * d10_v
          d12_v = d3_v * al1 / d11_v
          d2_b = 1.0d0 / d11_v
          d3_b = (-d12_v) / d11_v
          d6_b = d3_b * d6_v
          d7_b = d6_b * rnagc
          d8_b = d6_b * d8_v + d7_b * rnagc + d7_b * rnagc
          d11_b = d2_b * d3_v
          d9_b = d3_b * d10_v + d2_b * al1 * rnapol * ((-1.5d0) * rnaie)
          do g_i_ = 1, g_p_
            g_c60(g_i_) = d8_b * g_rnagc(g_i_) + d11_b * g_al1(g_i_) + d
     *9_b * g_h2ie(g_i_)
          enddo
          c60 = d12_v
C--------
          d2_v = c60 * f0
          d4_b = rnagc * f0
          do g_i_ = 1, g_p_
            g_c60(g_i_) = d2_v * g_rnagc(g_i_) + d4_b * g_c60(g_i_)
          enddo
          c60 = d2_v * rnagc
C--------
C
          d2_v = (-39.0625d0) * rnag
          d4_b = rnagc * (-39.0625d0)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_v * g_rnagc(g_i_) + d4_b * g_rnag(g_i_)
          enddo
          d1_w = d2_v * rnagc
          d3_v = c60 - c690
          d5_v = d3_v * cstsq
          d7_v = exp(d1_w)
          d1_p =  d7_v
          d8_v = 1.0d0 - d7_v
          d7_b = (-d5_v) * d1_p
          d8_b = d8_v * cstsq
          d9_b = d8_v * d3_v
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_c6(g_i_) = d7_b * g_d1_w(g_i_) + d9_b * g_cstsq(g_i_) + d8
     *_b * g_c60(g_i_) + d2_b * g_c690(g_i_)
          enddo
          c6 = c690 + d5_v * d8_v
C--------
C
          d3_v = (-1.5d0) * rnahie1 * hie
          d6_v = rnahie1 + hie
          d7_v = d3_v * rnag1c / d6_v
          d10_v = d7_v * rnahpol1 * hpol
          d4_b = div1 * hpol
          d5_b = d4_b * rnahpol1
          d6_b = d4_b * d7_v
          d7_b = d5_b * (1.0d0 / d6_v)
          d11_b = d7_b * d3_v
          d9_b = d5_b * ((-d7_v) / d6_v) + d7_b * rnag1c * hie * (-1.5d0
     *)
          do g_i_ = 1, g_p_
            g_c61(g_i_) = d10_v * g_div1(g_i_) + d6_b * g_rnahpol1(g_i_)
     * + d11_b * g_rnag1c(g_i_) + d9_b * g_rnahie1(g_i_)
          enddo
          c61 = d10_v * div1
C--------
          d2_v = (-39.0625d0) * rnag1c
          d4_b = rnag1 * (-39.0625d0)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_v * g_rnag1(g_i_) + d4_b * g_rnag1c(g_i_)
          enddo
          d1_w = d2_v * rnag1
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_dexp1(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          dexp1 = d2_v
C--------
          d4_v = rnahpol1b / rnahpol1
          d5_v = d4_v - 1.0d0
          d8_v = 1.0d0 + d5_v * dexp1
          d6_b = c61 * d5_v
          d7_b = c61 * dexp1
          d8_b = d7_b * (1.0d0 / rnahpol1)
          d9_b = d7_b * ((-d4_v) / rnahpol1)
          do g_i_ = 1, g_p_
            g_c61(g_i_) = d6_b * g_dexp1(g_i_) + d9_b * g_rnahpol1(g_i_)
     * + d8_b * g_rnahpol1b(g_i_) + d8_v * g_c61(g_i_)
          enddo
          c61 = c61 * d8_v
C--------
C
C
          d3_v = (-1.5d0) * rnahie3 * hie
          d6_v = rnahie3 + hie
          d7_v = d3_v * rnag3c / d6_v
          d10_v = d7_v * rnahpol3 * hpol
          d4_b = div3 * hpol
          d5_b = d4_b * rnahpol3
          d6_b = d4_b * d7_v
          d7_b = d5_b * (1.0d0 / d6_v)
          d11_b = d7_b * d3_v
          d9_b = d5_b * ((-d7_v) / d6_v) + d7_b * rnag3c * hie * (-1.5d0
     *)
          do g_i_ = 1, g_p_
            g_c63(g_i_) = d10_v * g_div3(g_i_) + d6_b * g_rnahpol3(g_i_)
     * + d11_b * g_rnag3c(g_i_) + d9_b * g_rnahie3(g_i_)
          enddo
          c63 = d10_v * div3
C--------
          d2_v = (-39.0625d0) * rnag3c
          d4_b = rnag3 * (-39.0625d0)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_v * g_rnag3(g_i_) + d4_b * g_rnag3c(g_i_)
          enddo
          d1_w = d2_v * rnag3
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_dexp3(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          dexp3 = d2_v
C--------
          d4_v = rnahpol3b / rnahpol3
          d5_v = d4_v - 1.0d0
          d8_v = 1.0d0 + d5_v * dexp3
          d6_b = c63 * d5_v
          d7_b = c63 * dexp3
          d8_b = d7_b * (1.0d0 / rnahpol3)
          d9_b = d7_b * ((-d4_v) / rnahpol3)
          do g_i_ = 1, g_p_
            g_c63(g_i_) = d6_b * g_dexp3(g_i_) + d9_b * g_rnahpol3(g_i_)
     * + d8_b * g_rnahpol3b(g_i_) + d8_v * g_c63(g_i_)
          enddo
          c63 = c63 * d8_v
C--------
C
          d4_v = 1.5d0 * cstsq1 + 0.5d0
          d5_v = 1.0d0 / d4_v
          d6_v = d5_v - 1.0d0
          d9_v = 1.0d0 + d6_v * dexp1
          d6_b = rmuimu1 * d6_v
          d10_b = rmuimu1 * dexp1 * ((-d5_v) / d4_v) * 1.5d0
          do g_i_ = 1, g_p_
            g_cind1(g_i_) = d6_b * g_dexp1(g_i_) + d10_b * g_cstsq1(g_i_
     *) + d9_v * g_rmuimu1(g_i_)
          enddo
          cind1 = rmuimu1 * d9_v
C--------
C
          d4_v = 1.5d0 * cstsq3 + 0.5d0
          d5_v = 1.0d0 / d4_v
          d6_v = d5_v - 1.0d0
          d9_v = 1.0d0 + d6_v * dexp3
          d6_b = rmuimu3 * d6_v
          d10_b = rmuimu3 * dexp3 * ((-d5_v) / d4_v) * 1.5d0
          do g_i_ = 1, g_p_
            g_cind3(g_i_) = d6_b * g_dexp3(g_i_) + d10_b * g_cstsq3(g_i_
     *) + d9_v * g_rmuimu3(g_i_)
          enddo
          cind3 = rmuimu3 * d9_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = balp * g_r(g_i_, i, 2)
          enddo
          d1_w = balp * (r(i, 2) - brho)
          d3_v = c6 * cevau * 0.5d0
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d6_v = 1.0d0 - d5_v
          d5_b = (-d3_v) * d1_p
          d7_b = d6_v * 0.5d0 * cevau
          do g_i_ = 1, g_p_
            g_c6(g_i_) = d5_b * g_d1_w(g_i_) + d7_b * g_c6(g_i_)
          enddo
          c6 = d3_v * d6_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dalp * g_r(g_i_, i, 1)
          enddo
          d1_w = dalp * (r(i, 1) - drho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_c6swt1(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          c6swt1 = 0.5d0 * (1.0d0 - d2_v)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dalp * g_r(g_i_, i, 3)
          enddo
          d1_w = dalp * (r(i, 3) - drho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_c6swt3(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          c6swt3 = 0.5d0 * (1.0d0 - d2_v)
C--------
C
C
          d3_b = cevau * c6swt1
          d4_b = cevau * c61
          do g_i_ = 1, g_p_
            g_c61(g_i_) = d4_b * g_c6swt1(g_i_) + d3_b * g_c61(g_i_)
          enddo
          c61 = c61 * c6swt1 * cevau
C--------
          d3_b = cevau * c6swt3
          d4_b = cevau * c63
          do g_i_ = 1, g_p_
            g_c63(g_i_) = d4_b * g_c6swt3(g_i_) + d3_b * g_c63(g_i_)
          enddo
          c63 = c63 * c6swt3 * cevau
C--------
          do g_i_ = 1, g_p_
            g_cind1(g_i_) = cind1 * g_c6swt1(g_i_) + c6swt1 * g_cind1(g_
     *i_)
          enddo
          cind1 = cind1 * c6swt1
C--------
          do g_i_ = 1, g_p_
            g_cind3(g_i_) = cind3 * g_c6swt3(g_i_) + c6swt3 * g_cind3(g_
     *i_)
          enddo
          cind3 = cind3 * c6swt3
C--------
C
          do g_i_ = 1, g_p_
            g_elec(g_i_) = g_c6(g_i_) + g_cind3(g_i_) + g_cind1(g_i_) + 
     *g_c63(g_i_) + g_c61(g_i_)
          enddo
          elec = c61 + c63 + cind1 + cind3 + c6
C--------
C
C =========================================================
C here is where it all gets put together . . .
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_s3(g_i_) + g_s2(g_i_) + g_s1(g_i_)
          enddo
          e(i) = s1 + s2 + s3 + deh2
C--------
C
          d2_b = 1.0d0 / cevau
          do g_i_ = 1, g_p_
            g_elec(g_i_) = d2_b * g_elec(g_i_)
          enddo
          elec = elec / cevau
C--------
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_elec(g_i_) + g_e(g_i_, i)
          enddo
          e(i) = e(i) + elec
C--------
C
C here we will try a post-LEPS correction term
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bect) * g_rmid(g_i_)
          enddo
          d1_w = (-bect) * (rmid - rect)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xect(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xect = d2_v
C--------
          d2_v = rave * rave
          d1_p = 2.0d0 * rave
          d3_b = (-bect2) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_rave(g_i_)
          enddo
          d1_w = (-bect2) * d2_v
          d2_v = dect * xect
          d3_v = xect - 2.0d0
          d4_v = d2_v * d3_v
          d6_v = exp(d1_w)
          d2_p =  d6_v
          d7_v = d4_v * d6_v
          d11_v = (0.5d0 * (1.0d0 - cosa)) * (0.5d0 * (1.0d0 - cosa))
          d1_p = 2.0d0 * (0.5d0 * (1.0d0 - cosa))
          d6_b = -(d7_v * d1_p * 0.5d0)
          d7_b = d11_v * d6_v
          d9_b = d11_v * d4_v * d2_p
          d12_b = d7_b * d2_v + d7_b * d3_v * dect
          do g_i_ = 1, g_p_
            g_vmid(g_i_) = d6_b * g_cosa(g_i_) + d9_b * g_d1_w(g_i_) + d
     *12_b * g_xect(g_i_)
          enddo
          vmid = d7_v * d11_v
C--------
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_vmid(g_i_) + g_e(g_i_, i)
          enddo
          e(i) = e(i) + vmid
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
C   This file was generated on 05/17/99 by the version of
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
C U12 coupling surface . . .
C
      subroutine g_prepot12(g_p_)
C
        implicit none
C
        double precision ap, alp, bet, rnagsq, rnag, cosg, rd0, r20
        double precision r1s, r2s, r3s, cos1, cos2, theta
        double precision gc, gc2
C
        double precision cevau
C
        integer nt, i, j, k
        logical leven
        double precision r(nt, 3), e(nt)
C
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d2_p, d2_w, d1_w, d4_b, d5_b, d10_v, d9_v, d2_v
     *, d8_b, d2_b
        double precision d1_p, d3_v, d4_v, d5_v, d7_b, d7_v, d6_b, d3_b,
     * g_r1s(g_pmax_), g_r(ldg_r, nt, 3)
        double precision g_r2s(g_pmax_), g_r3s(g_pmax_), g_rnagsq(g_pmax
     *_), g_rnag(g_pmax_), g_cos1(g_pmax_), g_cos2(g_pmax_), g_cosg(g_pm
     *ax_), g_gc(g_pmax_), g_gc2(g_pmax_), g_d1_w(g_pmax_)
        double precision g_d2_w(g_pmax_), g_e(ldg_e, nt)
        integer g_ehfid
        data ap /0.835d0/, alp /4.0d0/, rd0 /4.0d0/, bet /0.5d0/, r20 /-
     *1.3d0/
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
C set leven true for an even U12 surface, and false for an odd 
C surface
        leven = .true.
        return
C
        entry g_pot12(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do i = 1, nt
C
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
C here we do the function evaluations . . . 
C
          do g_i_ = 1, g_p_
            g_rnagsq(g_i_) = (-0.25d0) * g_r2s(g_i_) + 0.5d0 * g_r3s(g_i
     *_) + 0.5d0 * g_r1s(g_i_)
          enddo
          rnagsq = 0.5d0 * r1s + 0.5d0 * r3s - 0.25d0 * r2s
C--------
          if (rnagsq .gt. 0.0d0) then
            d2_v = sqrt(rnagsq)

            if ( rnagsq .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               d1_p = 0.0d0
            endif
            do g_i_ = 1, g_p_
              g_rnag(g_i_) = d1_p * g_rnagsq(g_i_)
            enddo
            rnag = d2_v
C--------
          else
            do g_i_ = 1, g_p_
              g_rnag(g_i_) = 0.0d0
            enddo
            rnag = 0.0d0
C--------
          endif
C
C cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
C
          d7_v = 2.0d0 * r(i, 1)
          d9_v = d7_v * r(i, 2)
          d10_v = (r1s + r2s - r3s) / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = (-d10_v) / d9_v
          d5_b = d3_b * d7_v
          d6_b = d3_b * r(i, 2) * 2.0d0
          do g_i_ = 1, g_p_
            g_cos1(g_i_) = d5_b * g_r(g_i_, i, 2) + d6_b * g_r(g_i_, i, 
     *1) + (-d2_b) * g_r3s(g_i_) + d2_b * g_r2s(g_i_) + d2_b * g_r1s(g_i
     *_)
          enddo
          cos1 = d10_v
C--------
          d7_v = 2.0d0 * r(i, 3)
          d9_v = d7_v * r(i, 2)
          d10_v = (r3s + r2s - r1s) / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = (-d10_v) / d9_v
          d5_b = d3_b * d7_v
          d6_b = d3_b * r(i, 2) * 2.0d0
          do g_i_ = 1, g_p_
            g_cos2(g_i_) = d5_b * g_r(g_i_, i, 2) + d6_b * g_r(g_i_, i, 
     *3) + (-d2_b) * g_r1s(g_i_) + d2_b * g_r2s(g_i_) + d2_b * g_r3s(g_i
     *_)
          enddo
          cos2 = d10_v
C--------
          do g_i_ = 1, g_p_
            g_cosg(g_i_) = (-0.5d0) * g_cos2(g_i_) + 0.5d0 * g_cos1(g_i_
     *)
          enddo
          cosg = 0.5d0 * (cos1 - cos2)
C--------
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
          d3_v = (gc - rd0) ** ( 4 - 2)
          d3_v =  d3_v * (gc - rd0)
          d1_p =  4 *  d3_v
          d3_v =  d3_v * (gc - rd0)
          d4_b = (-((1.0d0 / alp) ** 4)) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_gc(g_i_)
          enddo
          d1_w = (-((1.0d0 / alp) ** 4)) * d3_v
          d3_v = (gc2 - r20) * (gc2 - r20)
          d1_p = 2.0d0 * (gc2 - r20)
          d4_b = (-((1.0d0 / bet) ** 2)) * d1_p
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = d4_b * g_gc2(g_i_)
          enddo
          d2_w = (-((1.0d0 / bet) ** 2)) * d3_v
          d2_v = ap * cosg
          d4_v = exp(d1_w)
          d2_p =  d4_v
          d5_v = d2_v * d4_v
          d7_v = exp(d2_w)
          d1_p =  d7_v
          d4_b = d5_v * d1_p
          d7_b = d7_v * d2_v * d2_p
          d8_b = d7_v * d4_v * ap
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = d4_b * g_d2_w(g_i_) + d7_b * g_d1_w(g_i_) + d
     *8_b * g_cosg(g_i_)
          enddo
          e(i) = d5_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_e(g_i_, i)
          enddo
          e(i) = e(i) * cevau
C--------
C
C set leven equal to true for an even surface, false for an odd one
          if (leven) then
            d2_v = e(i) ** ( 4 - 2)
            d2_v =  d2_v * e(i)
            d2_p =  4 *  d2_v
            d2_v =  d2_v * e(i)
            d3_v = e(i) * e(i)
            d1_p = 2.0d0 * e(i)
            d4_v = d3_v + 1.0d-8
            d5_v = d2_v / d4_v
            d5_b = (-d5_v) / d4_v * d1_p + 1.0d0 / d4_v * d2_p
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d5_b * g_e(g_i_, i)
            enddo
            d1_w = d5_v
            d2_v = sqrt(d1_w)

            if ( d1_w .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               d1_p = 0.0d0
            endif
            do g_i_ = 1, g_p_
              g_e(g_i_, i) = d1_p * g_d1_w(g_i_)
            enddo
            e(i) = d2_v
C--------
          endif
C
        enddo
        return
      end
C
C ********************************************************************
C ********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 05/17/99 by the version of
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
        implicit double precision (a-h, o-z)
C
        dimension r(nt, 3), e(nt)
C
        save
C NaH triplets
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d3_p, d24_b, d2_p, d15_v, d16_v, d20_v, d22_v, 
     *d2_v, d3_v, d4_v
        double precision d16_b, d2_b, d3_b, d4_b, d1_w, d2_w, d5_v, d6_v
     *, d5_b, d6_b
        double precision d1_p, d7_v, d8_v, d9_v, d7_b, d8_b, d9_b, d10_v
     *, d10_b, d11_v
        double precision d12_v, d13_v, d19_v, d11_b, d12_b, d13_b, d14_b
     *, g_rmid(g_pmax_), g_r(ldg_r, nt, 3), g_r1s(g_pmax_)
        double precision g_r2s(g_pmax_), g_r3s(g_pmax_), g_d2_w(g_pmax_)
     *, g_d1_w(g_pmax_), g_rnag(g_pmax_), g_rnags(g_pmax_), g_rnagc(g_pm
     *ax_), g_cstsq(g_pmax_), g_cos1(g_pmax_), g_cos2(g_pmax_)
        double precision g_cosg(g_pmax_), g_rnag1(g_pmax_), g_rnag3(g_pm
     *ax_), g_rnag1c(g_pmax_), g_rnag3c(g_pmax_), g_cstsq1(g_pmax_), g_c
     *stsq3(g_pmax_), g_beta(g_pmax_), g_xenah(g_pmax_), g_s1(g_pmax_)
        double precision g_g10(g_pmax_), g_t10(g_pmax_), g_xet190(g_pmax
     *_), g_t190(g_pmax_), g_xet1c(g_pmax_), g_t1c(g_pmax_), g_t1m(g_pma
     *x_), g_t1(g_pmax_), g_xet20(g_pmax_), g_t20(g_pmax_)
        double precision g_xet290(g_pmax_), g_t290(g_pmax_), g_t2(g_pmax
     *_), g_beh2(g_pmax_), g_sswt(g_pmax_), g_deh2m(g_pmax_), g_xeh2(g_p
     *max_), g_s2a(g_pmax_), g_s2(g_pmax_), g_s20b(g_pmax_)
        double precision g_s3(g_pmax_), g_xet3(g_pmax_), g_t30(g_pmax_),
     * g_g30(g_pmax_), g_xet390(g_pmax_), g_t390(g_pmax_), g_xet3c(g_pma
     *x_), g_t3c(g_pmax_), g_t3m(g_pmax_), g_t3(g_pmax_)
        double precision g_h2ie(g_pmax_), g_al1(g_pmax_), g_al2(g_pmax_)
     *, g_alt(g_pmax_), g_alb(g_pmax_), g_alhh(g_pmax_), g_h2hex(g_pmax_
     *), g_h2quad(g_pmax_), g_rnahie1(g_pmax_), g_rnahie3(g_pmax_)
        double precision g_rnahmu1(g_pmax_), g_rnahmu3(g_pmax_), g_rnahp
     *ar1(g_pmax_), g_rnahpar3(g_pmax_), g_rnahperp1(g_pmax_), g_rnahper
     *p3(g_pmax_), g_rnahpol1(g_pmax_), g_rnahpol1b(g_pmax_), g_rnahpol3
     *(g_pmax_), g_rnahpol3b(g_pmax_)
        double precision g_qqint(g_pmax_), g_qq(g_pmax_), g_dexp2(g_pmax
     *_), g_qhint(g_pmax_), g_qh(g_pmax_), g_div1(g_pmax_), g_rmuimu1(g_
     *pmax_), g_div3(g_pmax_), g_rmuimu3(g_pmax_), g_c6(g_pmax_)
        double precision g_c61(g_pmax_), g_dexp1(g_pmax_), g_c63(g_pmax_
     *), g_dexp3(g_pmax_), g_cind1(g_pmax_), g_cind3(g_pmax_), g_elec(g_
     *pmax_), g_c6swt1(g_pmax_), g_c6swt3(g_pmax_), g_coul1(g_pmax_)
        double precision g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_p
     *max_), g_exch2(g_pmax_), g_exch3(g_pmax_), g_w(g_pmax_), g_cplg2(g
     *_pmax_), g_e(ldg_e, nt)
        integer g_ehfid
        data ret10 /3.22835d0/, bet10 /0.57874d0/, bet10b /2.44094d0/, f
     *10 /9.44882d0/
        data det190 /4.32441d0/, bet190 /2.43701d0/, ret190 /1.95480d0/,
     * t1alp /0.30778d0/, t1rho /3.21850d0/
        data det1c /4.92126d0/, bet1c /0.65748d0/, ret1c /2.71654d0/
C
C H2 triplets
        data det20 /3.29921d0/, bet20 /2.64567d0/, ret20 /2.44882d0/
        data det290 /3.97638d0/, bet290 /3.73937d0/, ret290 /0.96024d0/
        data t2alp /1.0d0/, t2rho /9.0d0/
C
C H2 singlets
        data des290 /0.91433d0/, s2alp /0.66778d0/, s2rho /5.30472d0/
        data deh2 /4.74806d0/
        data s20alp /0.5d0/, s20rho /6.0d0/
C
C misc constants
C
        data c2a /1.38661d0/, c2b /0.27362d0/, c2c /0.03764d0/, una2p /2
     *.1037d0/
C
C data blocks for the long range forces
        data h2a1 /2.38029900d0/, h2a2 /1.47400334d0/, h2a3 /-0.22633177
     *1d0/, h2a4 /-2.17339562d0/, h2a5 /9.54597998d0/, h2a6 /-23.3366749
     *d0/, h2a7 /24.7115748d0/, h2a8 /-12.8786102d0/, h2a9 /2.75711758d0
     */
        data h2b1 /1.35164094d0/, h2b2 /-0.254404564d0/, h2b3 /-10.66270
     *55d0/, h2b4 /-14.0112406d0/, h2b5 /-7.80660398d0/, h2b6 /5.4769775
     *3d0/, h2b7 /-6.95875664d0/, h2b8 /3.64334771d0/, h2b9 /-0.49250559
     *6d0/
        data h2i1 /160.484d0/, h2i2 /-20.397d0/, h2i3 /3.913d0/, h2i4 /0
     *.777d0/, h2i5 /-7.583/, h2i6 /1.501/
        data rnaie /3.0353/, rnapol1 /290.4d0/, hpol /4.5d0/, hie /13.59
     *8d0/, rnaquad /-16.8992d0/
        data hnai1 /-193.341d0/, hnai2 /16.968d0/, hnai3 /-77.120d0/, hn
     *ai4 /-1.024d0/, hnai5 /6.009d0/, hnai6 /-1.241d0/
        data hnau1 /119.880d0/, hnau2 /-0.620d0/, hnau3 /31.335d0/, hnau
     *4 /70.414d0/, hnau5 /-13.016d0/, hnau6 /0.8681d0/
        data hnaa1 /0.67d0/, hnaa2 /-92.229d0/, hnaa3 /51.01019d0/, hnaa
     *4 /-130.37826d0/, hnaa5 /6.40285d0/
        data hnab1 /0.34333d0/, hnab2 /-92.229d0/, hnab3 /-48.54867d0/, 
     *hnab4 /-17.30523d0/, hnab5 /0.99821d0/
        data rqq /9.766021100d4/, rqh /2.355466062927d6/, rb /9.31193833
     *71714d5/, balp /1.93701d0/, brho /1.4d0/, rd /9.1351724748364d8/, 
     *dalp /1.43307d0/, drho /5.71654d0/
        data h2hc1 /0.677d0/, h2hc2 /0.378d0/, h2hc3 /0.281d0/, h2hc4 /-
     *1.022d0/, h2hc5 /0.665d0/
        data h2qc1 /1.532d0/, h2qc2 /-.04635d0/, h2qc3 /1.826d0/, h2qc4 
     */-2.508d0/, h2qc5 /1.856d0/
C
C
C
        data g_ehfid /0/
C
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        cevau = 1.d0 / 27.2113961d0
        onethd = 1.0d0 / 3.0d0
        twothd = 2.0d0 / 3.0d0
C
        return
C
        entry g_pot22(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do 99999 i = 1, nt
C
C =========================================================
C this section calcluates useful quantities - the A to BC center
C of mass distance, the angles between A-BC and BC, the shorest
C of any of the three diatomic distances, and the bond angle
C H-Na-H
C
          do g_i_ = 1, g_p_
            g_rmid(g_i_) = 0.5d0 * g_r(g_i_, i, 3) + 0.5d0 * g_r(g_i_, i
     *, 1)
          enddo
          rmid = 0.5d0 * (r(i, 1) + r(i, 3))
C--------
          rave = 0.5d0 * (r(i, 1) - r(i, 3))
C
          d2_b = r(i, 1) + r(i, 1)
          do g_i_ = 1, g_p_
            g_r1s(g_i_) = d2_b * g_r(g_i_, i, 1)
          enddo
          r1s = r(i, 1) * r(i, 1)
C--------
          d2_b = r(i, 2) + r(i, 2)
          do g_i_ = 1, g_p_
            g_r2s(g_i_) = d2_b * g_r(g_i_, i, 2)
          enddo
          r2s = r(i, 2) * r(i, 2)
C--------
          d2_b = r(i, 3) + r(i, 3)
          do g_i_ = 1, g_p_
            g_r3s(g_i_) = d2_b * g_r(g_i_, i, 3)
          enddo
          r3s = r(i, 3) * r(i, 3)
C--------
C
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-0.5d0) * g_r2s(g_i_) + g_r3s(g_i_) + g_r1s(
     *g_i_)
          enddo
          d2_w = r1s + r3s - 0.5d0 * r2s
          d2_v = abs(d2_w)

          if (d2_w .gt. 0.0d0) then
             d1_p =  1.0d0
          else if (d2_w .lt. 0.0d0) then
             d1_p = -1.0d0
          else
             d1_p = 0.0d0
          endif
          d3_b = 0.5d0 * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_d2_w(g_i_)
          enddo
          d1_w = 0.5d0 * d2_v
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
          d2_b = rnag + rnag
          do g_i_ = 1, g_p_
            g_rnags(g_i_) = d2_b * g_rnag(g_i_)
          enddo
          rnags = rnag * rnag
C--------
          do g_i_ = 1, g_p_
            g_rnagc(g_i_) = rnags * g_rnag(g_i_) + rnag * g_rnags(g_i_)
          enddo
          rnagc = rnags * rnag
C--------
C
          if (rnag .eq. 0.0d0) then
            do g_i_ = 1, g_p_
              g_cstsq(g_i_) = 0.0d0
            enddo
            cstsq = 0.0d0
C--------
          else
            d7_v = r(i, 2) * rnag
            d8_v = 0.5d0 * (r3s - r1s) / d7_v
            d9_v = d8_v * d8_v
            d1_p = 2.0d0 * d8_v
            d4_b = d1_p * ((-d8_v) / d7_v)
            d5_b = d4_b * rnag
            d6_b = d4_b * r(i, 2)
            d7_b = d1_p * (1.0d0 / d7_v) * 0.5d0
            do g_i_ = 1, g_p_
              g_cstsq(g_i_) = d6_b * g_rnag(g_i_) + d5_b * g_r(g_i_, i, 
     *2) + (-d7_b) * g_r1s(g_i_) + d7_b * g_r3s(g_i_)
            enddo
            cstsq = d9_v
C--------
          endif
C
          cosa = (r1s + r3s - r2s) / (2 * r(i, 1) * r(i, 3))
C
          coscos = cstsq * (0.5d0 + 0.5d0 * cosa) * (1.0d0 - dexp((-10.0
     *d0) * rnags))
C
C cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up
C
          d7_v = 2.0d0 * r(i, 1)
          d9_v = d7_v * r(i, 2)
          d10_v = (r1s + r2s - r3s) / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = (-d10_v) / d9_v
          d5_b = d3_b * d7_v
          d6_b = d3_b * r(i, 2) * 2.0d0
          do g_i_ = 1, g_p_
            g_cos1(g_i_) = d5_b * g_r(g_i_, i, 2) + d6_b * g_r(g_i_, i, 
     *1) + (-d2_b) * g_r3s(g_i_) + d2_b * g_r2s(g_i_) + d2_b * g_r1s(g_i
     *_)
          enddo
          cos1 = d10_v
C--------
          d7_v = 2.0d0 * r(i, 3)
          d9_v = d7_v * r(i, 2)
          d10_v = (r3s + r2s - r1s) / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = (-d10_v) / d9_v
          d5_b = d3_b * d7_v
          d6_b = d3_b * r(i, 2) * 2.0d0
          do g_i_ = 1, g_p_
            g_cos2(g_i_) = d5_b * g_r(g_i_, i, 2) + d6_b * g_r(g_i_, i, 
     *3) + (-d2_b) * g_r1s(g_i_) + d2_b * g_r2s(g_i_) + d2_b * g_r3s(g_i
     *_)
          enddo
          cos2 = d10_v
C--------
          do g_i_ = 1, g_p_
            g_cosg(g_i_) = (-0.5d0) * g_cos2(g_i_) + 0.5d0 * g_cos1(g_i_
     *)
          enddo
          cosg = 0.5d0 * (cos1 - cos2)
C--------
C
C
          frac = 0.958d0
          d5_b = frac * frac + (-frac)
          d6_b = 1.0d0 + (-frac)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = d5_b * g_r1s(g_i_) + frac * g_r3s(g_i_) + d6_
     *b * g_r2s(g_i_)
          enddo
          d2_w = r2s + frac * (r3s - r2s - r1s) + r1s * frac * frac
          d2_v = abs(d2_w)

          if (d2_w .gt. 0.0d0) then
             d1_p =  1.0d0
          else if (d2_w .lt. 0.0d0) then
             d1_p = -1.0d0
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d1_p * g_d2_w(g_i_)
          enddo
          d1_w = d2_v
          d2_v = sqrt(d1_w)

          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_rnag1(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          rnag1 = d2_v
C--------
          d5_b = frac * frac + (-frac)
          d6_b = 1.0d0 + (-frac)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = d5_b * g_r3s(g_i_) + frac * g_r1s(g_i_) + d6_
     *b * g_r2s(g_i_)
          enddo
          d2_w = r2s + frac * (r1s - r2s - r3s) + r3s * frac * frac
          d2_v = abs(d2_w)

          if (d2_w .gt. 0.0d0) then
             d1_p =  1.0d0
          else if (d2_w .lt. 0.0d0) then
             d1_p = -1.0d0
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d1_p * g_d2_w(g_i_)
          enddo
          d1_w = d2_v
          d2_v = sqrt(d1_w)

          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
          else
             d1_p = 0.0d0
          endif
          do g_i_ = 1, g_p_
            g_rnag3(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          rnag3 = d2_v
C--------
C
          d2_v = rnag1 * rnag1
          d3_b = d2_v + rnag1 * rnag1 + rnag1 * rnag1
          do g_i_ = 1, g_p_
            g_rnag1c(g_i_) = d3_b * g_rnag1(g_i_)
          enddo
          rnag1c = d2_v * rnag1
C--------
          d2_v = rnag3 * rnag3
          d3_b = d2_v + rnag3 * rnag3 + rnag3 * rnag3
          do g_i_ = 1, g_p_
            g_rnag3c(g_i_) = d3_b * g_rnag3(g_i_)
          enddo
          rnag3c = d2_v * rnag3
C--------
C
          if (rnag1 .eq. 0.0d0) then
            do g_i_ = 1, g_p_
              g_cstsq1(g_i_) = 0.0d0
            enddo
            cstsq1 = 1.0d0
C--------
          else
            d11_v = r(i, 1) * rnag1
            d12_v = (r3s - r2s - r1s + 2.0d0 * r3s * frac) / d11_v
            d13_v = d12_v * d12_v
            d1_p = 2.0d0 * d12_v
            d3_b = 0.25d0 * d1_p
            d4_b = d3_b * (1.0d0 / d11_v)
            d5_b = d3_b * ((-d12_v) / d11_v)
            d6_b = d5_b * rnag1
            d7_b = d5_b * r(i, 1)
            d11_b = d4_b * frac * 2.0d0 + d4_b
            do g_i_ = 1, g_p_
              g_cstsq1(g_i_) = d7_b * g_rnag1(g_i_) + d6_b * g_r(g_i_, i
     *, 1) + (-d4_b) * g_r1s(g_i_) + (-d4_b) * g_r2s(g_i_) + d11_b * g_r
     *3s(g_i_)
            enddo
            cstsq1 = 0.25d0 * d13_v
C--------
          endif
C
          if (rnag3 .eq. 0.0d0) then
            do g_i_ = 1, g_p_
              g_cstsq3(g_i_) = 0.0d0
            enddo
            cstsq3 = 1.0d0
C--------
          else
            d11_v = r(i, 3) * rnag3
            d12_v = (r1s - r2s - r3s + 2.0d0 * r1s * frac) / d11_v
            d13_v = d12_v * d12_v
            d1_p = 2.0d0 * d12_v
            d3_b = 0.25d0 * d1_p
            d4_b = d3_b * (1.0d0 / d11_v)
            d5_b = d3_b * ((-d12_v) / d11_v)
            d6_b = d5_b * rnag3
            d7_b = d5_b * r(i, 3)
            d11_b = d4_b * frac * 2.0d0 + d4_b
            do g_i_ = 1, g_p_
              g_cstsq3(g_i_) = d7_b * g_rnag3(g_i_) + d6_b * g_r(g_i_, i
     *, 3) + (-d4_b) * g_r3s(g_i_) + (-d4_b) * g_r2s(g_i_) + d11_b * g_r
     *1s(g_i_)
            enddo
            cstsq3 = 0.25d0 * d13_v
C--------
          endif
C
C =========================================================
C
C The first NaH singlet
C
C new singlet, designed to minic lower adiabatic behavior
          denah = 1.97134878d0
          renah = 3.566044d0
          bf = 0.864d0
          b0 = 0.594d0
          gamnah = 7.376d0
C
          d3_v = (r(i, 1) / gamnah) ** ( 8 - 2)
          d3_v =  d3_v * (r(i, 1) / gamnah)
          d2_p =  8 *  d3_v
          d3_v =  d3_v * (r(i, 1) / gamnah)
          d7_v = (r(i, 1) / gamnah) ** ( 8 - 2)
          d7_v =  d7_v * (r(i, 1) / gamnah)
          d1_p =  8 *  d7_v
          d7_v =  d7_v * (r(i, 1) / gamnah)
          d8_v = bf + d7_v
          d9_v = bf * (b0 + d3_v) / d8_v
          d6_b = (-d9_v) / d8_v * d1_p * (1.0d0 / gamnah) + 1.0d0 / d8_v
     * * bf * d2_p * (1.0d0 / gamnah)
          do g_i_ = 1, g_p_
            g_beta(g_i_) = d6_b * g_r(g_i_, i, 1)
          enddo
          beta = d9_v
C--------
          d4_v = r(i, 1) - renah
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta) * g_r(g_i_, i, 1) + (-d4_v) * g_beta(
     *g_i_)
          enddo
          d1_w = (-beta) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xenah(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xenah = d2_v
C--------
C
          d2_v = denah * xenah
          d3_v = xenah - 2.0d0
          d4_b = d2_v + d3_v * denah
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d4_b * g_xenah(g_i_)
          enddo
          s1 = d2_v * d3_v
C--------
C
C The first NaH triplet
C
C
          do g_i_ = 1, g_p_
            g_g10(g_i_) = g_r(g_i_, i, 1)
          enddo
          g10 = r(i, 1) - ret10
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet10) * g_g10(g_i_)
          enddo
          d1_w = (-bet10) * g10
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-bet10b) * g_g10(g_i_)
          enddo
          d2_w = (-bet10b) * g10
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d5_v = exp(d2_w)
          d1_p =  d5_v
          d6_b = f10 * d2_p
          do g_i_ = 1, g_p_
            g_t10(g_i_) = d1_p * g_d2_w(g_i_) + d6_b * g_d1_w(g_i_)
          enddo
          t10 = f10 * d2_v + d5_v
C--------
          d2_b = 1.0d0 / (1.0d0 + f10)
          do g_i_ = 1, g_p_
            g_t10(g_i_) = d2_b * g_t10(g_i_)
          enddo
          t10 = t10 / (1.0d0 + f10)
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet190) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-bet190) * (r(i, 1) - ret190)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xet190(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xet190 = d2_v
C--------
          d2_v = det190 * onethd * xet190
          d3_v = xet190 + 2.0d0
          d4_b = d2_v + d3_v * (det190 * onethd)
          do g_i_ = 1, g_p_
            g_t190(g_i_) = d4_b * g_xet190(g_i_)
          enddo
          t190 = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet1c) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-bet1c) * (r(i, 1) - ret1c)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xet1c(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xet1c = d2_v
C--------
          d2_v = det1c * onethd * xet1c
          d3_v = xet1c + 2.0d0
          d4_b = d2_v + d3_v * (det1c * onethd)
          do g_i_ = 1, g_p_
            g_t1c(g_i_) = d4_b * g_xet1c(g_i_)
          enddo
          t1c = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = t1alp * g_r(g_i_, i, 2)
          enddo
          d1_w = t1alp * (r(i, 2) - t1rho)
          d4_v = (t1c - t190) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.0d0 + d6_v
          d7_b = d4_v * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_t1m(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_t1c(g_i_) + d2_
     *b * g_t190(g_i_)
          enddo
          t1m = t190 + d4_v * d7_v
C--------
C
          d3_v = t10 - t1m
          d5_v = cosg * cosg
          d1_p = 2.0d0 * cosg
          d6_b = d3_v * d1_p
          d2_b = 1.0d0 + (-d5_v)
          do g_i_ = 1, g_p_
            g_t1(g_i_) = d6_b * g_cosg(g_i_) + d5_v * g_t10(g_i_) + d2_b
     * * g_t1m(g_i_)
          enddo
          t1 = t1m + d3_v * d5_v
C--------
C
C
C The H2 triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet20) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-bet20) * (r(i, 2) - ret20)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xet20(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xet20 = d2_v
C--------
          d2_v = det20 * onethd * xet20
          d3_v = xet20 + 2.0d0
          d4_b = d2_v + d3_v * (det20 * onethd)
          do g_i_ = 1, g_p_
            g_t20(g_i_) = d4_b * g_xet20(g_i_)
          enddo
          t20 = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet290) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-bet290) * (r(i, 2) - ret290)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xet290(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xet290 = d2_v
C--------
          d2_v = det290 * onethd * xet290
          d3_v = xet290 + 2.0d0
          d4_b = d2_v + d3_v * (det290 * onethd)
          do g_i_ = 1, g_p_
            g_t290(g_i_) = d4_b * g_xet290(g_i_)
          enddo
          t290 = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = t2alp * g_rmid(g_i_)
          enddo
          d1_w = t2alp * (rmid - t2rho)
          d3_v = t290 - t20
          d5_v = cosg * cosg
          d2_p = 2.0d0 * cosg
          d6_v = 1.0d0 - d5_v
          d8_v = d3_v * d6_v * 0.5d0
          d10_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d10_v *  d10_v)
          d11_v = 1.0d0 - d10_v
          d7_b = (-d8_v) * d1_p
          d8_b = d11_v * 0.5d0
          d9_b = d8_b * d6_v
          d12_b = (-(d8_b * d3_v)) * d2_p
          d2_b = 1.0d0 + (-d9_b)
          do g_i_ = 1, g_p_
            g_t2(g_i_) = d7_b * g_d1_w(g_i_) + d12_b * g_cosg(g_i_) + d9
     *_b * g_t290(g_i_) + d2_b * g_t20(g_i_)
          enddo
          t2 = t20 + d8_v * d11_v
C--------
C
C The H2 singlet
C new singlet, designed to minic lower adiabatic behavior
          reh2 = 1.401d0
          bf = 1.627d0
          b0 = 1.194d0
          gamh2a = 2.323d0
          gamh2b = 3.126d0
C
          d5_v = (r(i, 2) / gamh2b) * (r(i, 2) / gamh2b)
          d2_p = 2.0d0 * (r(i, 2) / gamh2b)
          d11_v = (r(i, 2) / gamh2b) * (r(i, 2) / gamh2b)
          d1_p = 2.0d0 * (r(i, 2) / gamh2b)
          d12_v = bf - r(i, 2) / gamh2a + d11_v
          d13_v = bf * (b0 - r(i, 2) / gamh2a + d5_v) / d12_v
          d3_b = (-d13_v) / d12_v
          d9_b = 1.0d0 / d12_v * bf
          d7_b = d3_b * d1_p * (1.0d0 / gamh2b) + (-d3_b) * (1.0d0 / gam
     *h2a) + d9_b * d2_p * (1.0d0 / gamh2b) + (-d9_b) * (1.0d0 / gamh2a)
          do g_i_ = 1, g_p_
            g_beh2(g_i_) = d7_b * g_r(g_i_, i, 2)
          enddo
          beh2 = d13_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s2alp * g_rmid(g_i_)
          enddo
          d1_w = s2alp * (rmid - s2rho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_sswt(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          sswt = 0.5d0 * (1.0d0 - d2_v)
C--------
          d3_b = deh2 * (des290 - 1.0d0)
          do g_i_ = 1, g_p_
            g_deh2m(g_i_) = d3_b * g_sswt(g_i_)
          enddo
          deh2m = deh2 + deh2 * (des290 - 1.0d0) * sswt
C--------
C
          d4_v = r(i, 2) - reh2
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beh2) * g_r(g_i_, i, 2) + (-d4_v) * g_beh2(
     *g_i_)
          enddo
          d1_w = (-beh2) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xeh2(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xeh2 = d2_v
C--------
          d3_v = deh2m * xeh2
          d4_v = xeh2 - 2.0d0
          d6_b = d4_v * xeh2
          d5_b = d3_v + d4_v * deh2m
          do g_i_ = 1, g_p_
            g_s2a(g_i_) = d5_b * g_xeh2(g_i_) + d6_b * g_deh2m(g_i_)
          enddo
          s2a = d3_v * d4_v + una2p
C--------
C
          d2_b = 1.0d0 / 0.1d0
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-d2_b) * g_t2(g_i_) + d2_b * g_s2a(g_i_)
          enddo
          d1_w = (s2a - t2) / 0.1d0
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d6_v = s2a - t2
          d6_b = (-0.5d0) * d5_v
          d9_b = (-0.5d0) * d6_v * d1_p
          d7_b = d6_b + 0.5d0
          d8_b = -d6_b + 0.5d0
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d9_b * g_d1_w(g_i_) + d8_b * g_t2(g_i_) + d7_b 
     ** g_s2a(g_i_)
          enddo
          s2 = 0.5d0 * (s2a + t2 - d5_v * d6_v)
C--------
C
          d2_v = (-1.34d0) * beh2
          d4_v = r(i, 2) - reh2
          d5_b = d4_v * (-1.34d0)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_v * g_r(g_i_, i, 2) + d5_b * g_beh2(g_i_)
          enddo
          d1_w = d2_v * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xeh2(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xeh2 = d2_v
C--------
          d2_v = (deh2 - una2p) * xeh2
          d3_v = xeh2 - 2.0d0
          d4_b = d2_v + d3_v * (deh2 - una2p)
          do g_i_ = 1, g_p_
            g_s20b(g_i_) = d4_b * g_xeh2(g_i_)
          enddo
          s20b = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s20alp * g_rmid(g_i_)
          enddo
          d1_w = s20alp * (rmid - s20rho)
          d3_v = s20b - s2
          d5_v = cosg * cosg
          d2_p = 2.0d0 * cosg
          d7_v = d3_v * d5_v * 0.5d0
          d9_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d9_v *  d9_v)
          d10_v = 1.0d0 - d9_v
          d7_b = (-d7_v) * d1_p
          d8_b = d10_v * 0.5d0
          d9_b = d8_b * d5_v
          d11_b = d8_b * d3_v * d2_p
          d2_b = 1.0d0 + (-d9_b)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d7_b * g_d1_w(g_i_) + d11_b * g_cosg(g_i_) + d9
     *_b * g_s20b(g_i_) + d2_b * g_s2(g_i_)
          enddo
          s2 = s2 + d7_v * d10_v
C--------
C
C
C The second NaH singlet
C new singlet, designed to minic lower adiabatic behavior
          denah = 1.97134878d0
          renah = 3.566044d0
          bf = 0.864d0
          b0 = 0.594d0
          gamnah = 7.376d0
C
          d3_v = (r(i, 3) / gamnah) ** ( 8 - 2)
          d3_v =  d3_v * (r(i, 3) / gamnah)
          d2_p =  8 *  d3_v
          d3_v =  d3_v * (r(i, 3) / gamnah)
          d7_v = (r(i, 3) / gamnah) ** ( 8 - 2)
          d7_v =  d7_v * (r(i, 3) / gamnah)
          d1_p =  8 *  d7_v
          d7_v =  d7_v * (r(i, 3) / gamnah)
          d8_v = bf + d7_v
          d9_v = bf * (b0 + d3_v) / d8_v
          d6_b = (-d9_v) / d8_v * d1_p * (1.0d0 / gamnah) + 1.0d0 / d8_v
     * * bf * d2_p * (1.0d0 / gamnah)
          do g_i_ = 1, g_p_
            g_beta(g_i_) = d6_b * g_r(g_i_, i, 3)
          enddo
          beta = d9_v
C--------
          d4_v = r(i, 3) - renah
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-beta) * g_r(g_i_, i, 3) + (-d4_v) * g_beta(
     *g_i_)
          enddo
          d1_w = (-beta) * d4_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xenah(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xenah = d2_v
C--------
C
          d2_v = denah * xenah
          d3_v = xenah - 2.0d0
          d4_b = d2_v + d3_v * denah
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d4_b * g_xenah(g_i_)
          enddo
          s3 = d2_v * d3_v
C--------
C
C The second NaH triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet1) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bet1) * (r(i, 3) - ret1)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xet3(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xet3 = d2_v
C--------
          d2_v = det1 * onethd * xet3
          d3_v = xet3 + 2.0d0
          d4_b = d2_v + d3_v * (det1 * onethd)
          do g_i_ = 1, g_p_
            g_t30(g_i_) = d4_b * g_xet3(g_i_)
          enddo
          t30 = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_g30(g_i_) = g_r(g_i_, i, 3)
          enddo
          g30 = r(i, 3) - ret10
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet10) * g_g30(g_i_)
          enddo
          d1_w = (-bet10) * g30
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-bet10b) * g_g30(g_i_)
          enddo
          d2_w = (-bet10b) * g30
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d5_v = exp(d2_w)
          d1_p =  d5_v
          d6_b = f10 * d2_p
          do g_i_ = 1, g_p_
            g_t30(g_i_) = d1_p * g_d2_w(g_i_) + d6_b * g_d1_w(g_i_)
          enddo
          t30 = f10 * d2_v + d5_v
C--------
          d2_b = 1.0d0 / (1.0d0 + f10)
          do g_i_ = 1, g_p_
            g_t30(g_i_) = d2_b * g_t30(g_i_)
          enddo
          t30 = t30 / (1.0d0 + f10)
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet190) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bet190) * (r(i, 3) - ret190)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xet390(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xet390 = d2_v
C--------
          d2_v = det190 * onethd * xet390
          d3_v = xet390 + 2.0d0
          d4_b = d2_v + d3_v * (det190 * onethd)
          do g_i_ = 1, g_p_
            g_t390(g_i_) = d4_b * g_xet390(g_i_)
          enddo
          t390 = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet1c) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bet1c) * (r(i, 3) - ret1c)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xet3c(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xet3c = d2_v
C--------
          d2_v = det1c * onethd * xet3c
          d3_v = xet3c + 2.0d0
          d4_b = d2_v + d3_v * (det1c * onethd)
          do g_i_ = 1, g_p_
            g_t3c(g_i_) = d4_b * g_xet3c(g_i_)
          enddo
          t3c = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = t1alp * g_r(g_i_, i, 2)
          enddo
          d1_w = t1alp * (r(i, 2) - t1rho)
          d4_v = (t3c - t390) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.0d0 + d6_v
          d7_b = d4_v * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_t3m(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_t3c(g_i_) + d2_
     *b * g_t390(g_i_)
          enddo
          t3m = t390 + d4_v * d7_v
C--------
C
          d3_v = t30 - t3m
          d5_v = cosg * cosg
          d1_p = 2.0d0 * cosg
          d6_b = d3_v * d1_p
          d2_b = 1.0d0 + (-d5_v)
          do g_i_ = 1, g_p_
            g_t3(g_i_) = d6_b * g_cosg(g_i_) + d5_v * g_t30(g_i_) + d2_b
     * * g_t3m(g_i_)
          enddo
          t3 = t3m + d3_v * d5_v
C--------
C
C
C =========================================================
C this long section contains data pertinent to long range forces
C
C  H2 ionization energy:
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-h2i3) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-h2i3) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-h2i6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-h2i6) * r(i, 2)
          d3_v = h2i1 + h2i2 * r(i, 2)
          d5_v = exp(d1_w)
          d2_p =  d5_v
          d9_v = h2i4 + h2i5 * r(i, 2)
          d10_v = r2s * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d9_b = d10_v * d1_p
          d10_b = d12_v * d9_v
          d16_b = d3_v * d2_p
          d13_b = d12_v * r2s * h2i5 + d5_v * h2i2
          do g_i_ = 1, g_p_
            g_h2ie(g_i_) = -g_s2a(g_i_) + d9_b * g_d2_w(g_i_) + d10_b * 
     *g_r2s(g_i_) + d16_b * g_d1_w(g_i_) + d13_b * g_r(g_i_, i, 2)
          enddo
          h2ie = d3_v * d5_v + d10_v * d12_v - s2a + hie
C--------
C
C  H2 polarizability:
C
C   al1 is the parallel polarizability for h2
C   al2 is the perpendicaular polarizability for h2
C   alt is the anisotropic? polarizability
C   alb is the mean polarizability
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-h2a1) * g_r(g_i_, i, 2)
          enddo
          d1_w = h2a1 * (h2a2 - r(i, 2))
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = h2a6 * r(i, 2)
          d12_v = h2a7 * r2s
          d15_v = h2a8 * r2s
          d16_v = d15_v * r2s
          d19_v = h2a9 * r2s
          d20_v = d19_v * r2s
          d22_v = h2a3 + h2a4 * r(i, 2) + h2a5 * r2s + d9_v * r2s + d12_
     *v * r2s + d16_v * r(i, 2) + d20_v * r2s
          d7_b = d2_v * r2s
          d12_b = d2_v * r(i, 2)
          d8_b = d2_v * d20_v + d7_b * d19_v + d7_b * r2s * h2a9 + d12_b
     * * d15_v + d12_b * r2s * h2a8 + d2_v * d12_v + d2_v * r2s * h2a7 +
     * d2_v * d9_v + d2_v * h2a5
          d13_b = d2_v * d16_v + d2_v * r2s * h2a6 + d2_v * h2a4
          d24_b = d22_v * d1_p
          do g_i_ = 1, g_p_
            g_al1(g_i_) = d8_b * g_r2s(g_i_) + d13_b * g_r(g_i_, i, 2) +
     * d24_b * g_d1_w(g_i_)
          enddo
          al1 = d2_v * d22_v + 2.0d0 * hpol
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-h2b1) * g_r(g_i_, i, 2)
          enddo
          d1_w = h2b1 * (h2b2 - r(i, 2))
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = r(i, 2) ** ( 3 - 2)
          d9_v =  d9_v * r(i, 2)
          d1_p =  3 *  d9_v
          d9_v =  d9_v * r(i, 2)
          d12_v = h2b7 * r2s
          d15_v = h2b8 * r2s
          d16_v = d15_v * r2s
          d19_v = h2b9 * r2s
          d20_v = d19_v * r2s
          d22_v = h2b3 + h2b4 * r(i, 2) + h2b5 * r2s + h2b6 * d9_v + d12
     *_v * r2s + d16_v * r(i, 2) + d20_v * r2s
          d7_b = d2_v * r2s
          d12_b = d2_v * r(i, 2)
          d8_b = d2_v * d20_v + d7_b * d19_v + d7_b * r2s * h2b9 + d12_b
     * * d15_v + d12_b * r2s * h2b8 + d2_v * d12_v + d2_v * r2s * h2b7 +
     * d2_v * h2b5
          d13_b = d2_v * d16_v + d2_v * h2b6 * d1_p + d2_v * h2b4
          d24_b = d22_v * d2_p
          do g_i_ = 1, g_p_
            g_al2(g_i_) = d8_b * g_r2s(g_i_) + d13_b * g_r(g_i_, i, 2) +
     * d24_b * g_d1_w(g_i_)
          enddo
          al2 = d2_v * d22_v + 2.0d0 * hpol
C--------
C
          do g_i_ = 1, g_p_
            g_alt(g_i_) = (-twothd) * g_al2(g_i_) + twothd * g_al1(g_i_)
          enddo
          alt = twothd * (al1 - al2)
C--------
          do g_i_ = 1, g_p_
            g_alb(g_i_) = twothd * g_al2(g_i_) + onethd * g_al1(g_i_)
          enddo
          alb = onethd * al1 + twothd * al2
C--------
          d5_v = 1.5d0 * cstsq - 0.5d0
          d7_b = alt * 1.5d0
          do g_i_ = 1, g_p_
            g_alhh(g_i_) = d7_b * g_cstsq(g_i_) + d5_v * g_alt(g_i_) + g
     *_alb(g_i_)
          enddo
          alhh = alb + alt * d5_v
C--------
C
C
C here follow properties necesary for the long range forces in
C the Na + H2 channel:
C
C H2 properties:
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-h2hc1) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-h2hc1) * r(i, 2)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d8_v = h2hc4 * r(i, 2)
          d11_v = h2hc5 * r2s
          d13_v = h2hc2 * r(i, 2) + h2hc3 * r2s + d8_v * r2s + d11_v * r
     *2s
          d7_b = d2_v * d11_v + d2_v * r2s * h2hc5 + d2_v * d8_v + d2_v 
     ** h2hc3
          d11_b = d2_v * r2s * h2hc4 + d2_v * h2hc2
          d14_b = d13_v * d1_p
          do g_i_ = 1, g_p_
            g_h2hex(g_i_) = d7_b * g_r2s(g_i_) + d11_b * g_r(g_i_, i, 2)
     * + d14_b * g_d1_w(g_i_)
          enddo
          h2hex = d2_v * d13_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-h2qc1) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-h2qc1) * r(i, 2)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d8_v = h2qc4 * r(i, 2)
          d11_v = h2qc5 * r2s
          d13_v = h2qc2 * r(i, 2) + h2qc3 * r2s + d8_v * r2s + d11_v * r
     *2s
          d7_b = d2_v * d11_v + d2_v * r2s * h2qc5 + d2_v * d8_v + d2_v 
     ** h2qc3
          d11_b = d2_v * r2s * h2qc4 + d2_v * h2qc2
          d14_b = d13_v * d1_p
          do g_i_ = 1, g_p_
            g_h2quad(g_i_) = d7_b * g_r2s(g_i_) + d11_b * g_r(g_i_, i, 2
     *) + d14_b * g_d1_w(g_i_)
          enddo
          h2quad = d2_v * d13_v
C--------
C
C NaH properties
C
          d9_v = hnai6 * r(i, 1)
          d11_v = hnai3 + hnai4 * r(i, 1) + hnai5 * r1s + d9_v * r1s
          d12_v = (hnai1 + hnai2 * r(i, 1)) / d11_v
          d5_b = (-d12_v) / d11_v
          d9_b = d5_b * d9_v + d5_b * hnai5
          d10_b = d5_b * r1s * hnai6 + d5_b * hnai4 + 1.0d0 / d11_v * hn
     *ai2
          do g_i_ = 1, g_p_
            g_rnahie1(g_i_) = d9_b * g_r1s(g_i_) + d10_b * g_r(g_i_, i, 
     *1)
          enddo
          rnahie1 = d12_v + rnaie + una2p
C--------
C
          d9_v = hnai6 * r(i, 3)
          d11_v = hnai3 + hnai4 * r(i, 3) + hnai5 * r3s + d9_v * r3s
          d12_v = (hnai1 + hnai2 * r(i, 3)) / d11_v
          d5_b = (-d12_v) / d11_v
          d9_b = d5_b * d9_v + d5_b * hnai5
          d10_b = d5_b * r3s * hnai6 + d5_b * hnai4 + 1.0d0 / d11_v * hn
     *ai2
          do g_i_ = 1, g_p_
            g_rnahie3(g_i_) = d9_b * g_r3s(g_i_) + d10_b * g_r(g_i_, i, 
     *3)
          enddo
          rnahie3 = d12_v + rnaie + una2p
C--------
C
          d10_v = hnau6 * r(i, 1)
          d12_v = hnau3 + hnau4 * r(i, 1) + hnau5 * r1s + d10_v * r1s
          d13_v = (hnau1 * r(i, 1) + hnau2 * r1s) / d12_v
          d2_b = 1.0d0 / d12_v
          d3_b = (-d13_v) / d12_v
          d7_b = d3_b * d10_v + d3_b * hnau5 + d2_b * hnau2
          d8_b = d3_b * r1s * hnau6 + d3_b * hnau4 + d2_b * hnau1
          do g_i_ = 1, g_p_
            g_rnahmu1(g_i_) = d7_b * g_r1s(g_i_) + d8_b * g_r(g_i_, i, 1
     *)
          enddo
          rnahmu1 = d13_v
C--------
C
          d10_v = hnau6 * r(i, 3)
          d12_v = hnau3 + hnau4 * r(i, 3) + hnau5 * r3s + d10_v * r3s
          d13_v = (hnau1 * r(i, 3) + hnau2 * r3s) / d12_v
          d2_b = 1.0d0 / d12_v
          d3_b = (-d13_v) / d12_v
          d7_b = d3_b * d10_v + d3_b * hnau5 + d2_b * hnau2
          d8_b = d3_b * r3s * hnau6 + d3_b * hnau4 + d2_b * hnau1
          do g_i_ = 1, g_p_
            g_rnahmu3(g_i_) = d7_b * g_r3s(g_i_) + d8_b * g_r(g_i_, i, 3
     *)
          enddo
          rnahmu3 = d13_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-hnaa1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-hnaa1) * r(i, 1)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = hnaa5 * r1s
          d11_v = hnaa2 + hnaa3 * r(i, 1) + hnaa4 * r1s + d9_v * r(i, 1)
          d9_b = d2_v * r(i, 1) * hnaa5 + d2_v * hnaa4
          d8_b = d2_v * d9_v + d2_v * hnaa3
          d13_b = d11_v * d1_p
          do g_i_ = 1, g_p_
            g_rnahpar1(g_i_) = d9_b * g_r1s(g_i_) + d8_b * g_r(g_i_, i, 
     *1) + d13_b * g_d1_w(g_i_)
          enddo
          rnahpar1 = 159.267d0 + hpol + d2_v * d11_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-hnaa1) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-hnaa1) * r(i, 3)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = hnaa5 * r3s
          d11_v = hnaa2 + hnaa3 * r(i, 3) + hnaa4 * r3s + d9_v * r(i, 3)
          d9_b = d2_v * r(i, 3) * hnaa5 + d2_v * hnaa4
          d8_b = d2_v * d9_v + d2_v * hnaa3
          d13_b = d11_v * d1_p
          do g_i_ = 1, g_p_
            g_rnahpar3(g_i_) = d9_b * g_r3s(g_i_) + d8_b * g_r(g_i_, i, 
     *3) + d13_b * g_d1_w(g_i_)
          enddo
          rnahpar3 = 159.267d0 + hpol + d2_v * d11_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-hnab1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-hnab1) * r(i, 1)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = hnab5 * r1s
          d11_v = hnab2 + hnab3 * r(i, 1) + hnab4 * r1s + d9_v * r(i, 1)
          d9_b = d2_v * r(i, 1) * hnab5 + d2_v * hnab4
          d8_b = d2_v * d9_v + d2_v * hnab3
          d13_b = d11_v * d1_p
          do g_i_ = 1, g_p_
            g_rnahperp1(g_i_) = d9_b * g_r1s(g_i_) + d8_b * g_r(g_i_, i,
     * 1) + d13_b * g_d1_w(g_i_)
          enddo
          rnahperp1 = 159.267d0 + hpol + d2_v * d11_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-hnab1) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-hnab1) * r(i, 3)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d9_v = hnab5 * r3s
          d11_v = hnab2 + hnab3 * r(i, 3) + hnab4 * r3s + d9_v * r(i, 3)
          d9_b = d2_v * r(i, 3) * hnab5 + d2_v * hnab4
          d8_b = d2_v * d9_v + d2_v * hnab3
          d13_b = d11_v * d1_p
          do g_i_ = 1, g_p_
            g_rnahperp3(g_i_) = d9_b * g_r3s(g_i_) + d8_b * g_r(g_i_, i,
     * 3) + d13_b * g_d1_w(g_i_)
          enddo
          rnahperp3 = 159.267d0 + hpol + d2_v * d11_v
C--------
C
          d3_v = rnahpar1 - rnahperp1
          d2_b = 1.0d0 + (-cstsq1)
          do g_i_ = 1, g_p_
            g_rnahpol1(g_i_) = d3_v * g_cstsq1(g_i_) + cstsq1 * g_rnahpa
     *r1(g_i_) + d2_b * g_rnahperp1(g_i_)
          enddo
          rnahpol1 = rnahperp1 + d3_v * cstsq1
C--------
          do g_i_ = 1, g_p_
            g_rnahpol1b(g_i_) = twothd * g_rnahperp1(g_i_) + onethd * g_
     *rnahpar1(g_i_)
          enddo
          rnahpol1b = onethd * rnahpar1 + twothd * rnahperp1
C--------
C
          d3_v = rnahpar3 - rnahperp3
          d2_b = 1.0d0 + (-cstsq3)
          do g_i_ = 1, g_p_
            g_rnahpol3(g_i_) = d3_v * g_cstsq3(g_i_) + cstsq3 * g_rnahpa
     *r3(g_i_) + d2_b * g_rnahperp3(g_i_)
          enddo
          rnahpol3 = rnahperp3 + d3_v * cstsq3
C--------
          do g_i_ = 1, g_p_
            g_rnahpol3b(g_i_) = twothd * g_rnahperp3(g_i_) + onethd * g_
     *rnahpar3(g_i_)
          enddo
          rnahpol3b = onethd * rnahpar3 + twothd * rnahperp3
C--------
C
C electrostatic forces
C
C quadrupole quadrupole
          do g_i_ = 1, g_p_
            g_qqint(g_i_) = (-7.0d0) * g_cstsq(g_i_)
          enddo
          qqint = 3.0d0 - 7.0d0 * cstsq
C--------
          d5_v = rnagc * rnag
          d7_v = d5_v * rnag + rqq
          d8_v = 0.75d0 * rnaquad * h2quad / d7_v
          d6_b = qqint * ((-d8_v) / d7_v)
          d7_b = d6_b * rnag
          d9_b = d7_b * rnag
          d8_b = d6_b * d5_v + d7_b * rnagc
          d10_b = qqint * (1.0d0 / d7_v) * (0.75d0 * rnaquad)
          do g_i_ = 1, g_p_
            g_qq(g_i_) = d8_v * g_qqint(g_i_) + d8_b * g_rnag(g_i_) + d9
     *_b * g_rnagc(g_i_) + d10_b * g_h2quad(g_i_)
          enddo
          qq = d8_v * qqint
C--------
          d2_v = (-39.0625d0) * rnagc
          d4_b = rnag * (-39.0625d0)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_v * g_rnag(g_i_) + d4_b * g_rnagc(g_i_)
          enddo
          d1_w = d2_v * rnag
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_dexp2(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          dexp2 = d2_v
C--------
          d3_v = 3.0d0 * qqint
          d4_v = (-4.0d0) * pi / d3_v
          d5_v = d4_v - 1.0d0
          d8_v = 1.0d0 + d5_v * dexp2
          d6_b = qq * d5_v
          d9_b = qq * dexp2 * ((-d4_v) / d3_v) * 3.0d0
          do g_i_ = 1, g_p_
            g_qq(g_i_) = d6_b * g_dexp2(g_i_) + d9_b * g_qqint(g_i_) + d
     *8_v * g_qq(g_i_)
          enddo
          qq = qq * d8_v
C--------
C
C quadrupole hexadecapole
          d2_v = 147.0d0 * cstsq
          d4_b = 162.0d0 + (-d2_v) + (-cstsq) * 147.0d0
          do g_i_ = 1, g_p_
            g_qhint(g_i_) = d4_b * g_cstsq(g_i_)
          enddo
          qhint = -45.0d0 - d2_v * cstsq + 162.0d0 * cstsq
C--------
          d4_v = rnagc * rnagc
          d7_v = d4_v * rnag + rqh
          d8_v = 0.3125d0 * rnaquad * h2hex / d7_v
          d6_b = qhint * ((-d8_v) / d7_v)
          d7_b = d6_b * rnag
          d8_b = d6_b * d4_v
          d9_b = d7_b * rnagc + d7_b * rnagc
          d10_b = qhint * (1.0d0 / d7_v) * (0.3125d0 * rnaquad)
          do g_i_ = 1, g_p_
            g_qh(g_i_) = d8_v * g_qhint(g_i_) + d8_b * g_rnag(g_i_) + d9
     *_b * g_rnagc(g_i_) + d10_b * g_h2hex(g_i_)
          enddo
          qh = d8_v * qhint
C--------
          d3_v = (-40.8d0) * pi / qhint
          d4_v = d3_v - 1.0d0
          d7_v = 1.0d0 + d4_v * dexp2
          d6_b = qh * d4_v
          d8_b = qh * dexp2 * ((-d3_v) / qhint)
          do g_i_ = 1, g_p_
            g_qh(g_i_) = d6_b * g_dexp2(g_i_) + d8_b * g_qhint(g_i_) + d
     *7_v * g_qh(g_i_)
          enddo
          qh = qh * d7_v
C--------
C
C induction forces
C dipole induced dipole
C
          d2_v = rnag1c * rnag1c
          d4_v = d2_v * rnag1c + rd
          d5_v = 1.0d0 / d4_v
          d3_b = (-d5_v) / d4_v
          d4_b = d3_b * rnag1c
          d5_b = d3_b * d2_v + d4_b * rnag1c + d4_b * rnag1c
          do g_i_ = 1, g_p_
            g_div1(g_i_) = d5_b * g_rnag1c(g_i_)
          enddo
          div1 = d5_v
C--------
          d2_v = rnahmu1 * rnahmu1
          d1_p = 2.0d0 * rnahmu1
          d4_v = (-d2_v) * hpol
          d7_v = 1.5d0 * cstsq1 + 0.5d0
          d8_v = d4_v * d7_v
          d10_v = d8_v * rnag1c
          d4_b = div1 * rnag1c
          d5_b = div1 * d8_v
          d9_b = d4_b * d4_v * 1.5d0
          d12_b = (-(d4_b * d7_v * hpol)) * d1_p
          do g_i_ = 1, g_p_
            g_rmuimu1(g_i_) = d10_v * g_div1(g_i_) + d5_b * g_rnag1c(g_i
     *_) + d9_b * g_cstsq1(g_i_) + d12_b * g_rnahmu1(g_i_)
          enddo
          rmuimu1 = d10_v * div1
C--------
C
          d2_v = rnag3c * rnag3c
          d4_v = d2_v * rnag3c + rd
          d5_v = 1.0d0 / d4_v
          d3_b = (-d5_v) / d4_v
          d4_b = d3_b * rnag3c
          d5_b = d3_b * d2_v + d4_b * rnag3c + d4_b * rnag3c
          do g_i_ = 1, g_p_
            g_div3(g_i_) = d5_b * g_rnag3c(g_i_)
          enddo
          div3 = d5_v
C--------
          d2_v = rnahmu3 * rnahmu3
          d1_p = 2.0d0 * rnahmu3
          d4_v = (-d2_v) * hpol
          d7_v = 1.5d0 * cstsq3 + 0.5d0
          d8_v = d4_v * d7_v
          d10_v = d8_v * rnag3c
          d4_b = div3 * rnag3c
          d5_b = div3 * d8_v
          d9_b = d4_b * d4_v * 1.5d0
          d12_b = (-(d4_b * d7_v * hpol)) * d1_p
          do g_i_ = 1, g_p_
            g_rmuimu3(g_i_) = d10_v * g_div3(g_i_) + d5_b * g_rnag3c(g_i
     *_) + d9_b * g_cstsq3(g_i_) + d12_b * g_rnahmu3(g_i_)
          enddo
          rmuimu3 = d10_v * div3
C--------
C
C dispersion forces
C
          d3_v = rnaie + h2ie
          d4_v = (-1.5d0) * rnaie * h2ie / d3_v
          d5_v = d4_v * rnapol1
          d10_v = rnagc * rnagc + rb
          d11_v = d5_v * alhh / d10_v
          d2_b = 1.0d0 / d10_v
          d4_b = (-d11_v) / d10_v
          d5_b = d4_b * rnagc + d4_b * rnagc
          d7_b = d2_b * d5_v
          d8_b = d2_b * alhh * rnapol1
          d11_b = d8_b * ((-d4_v) / d3_v) + d8_b * (1.0d0 / d3_v) * ((-1
     *.5d0) * rnaie)
          do g_i_ = 1, g_p_
            g_c6(g_i_) = d5_b * g_rnagc(g_i_) + d7_b * g_alhh(g_i_) + d1
     *1_b * g_h2ie(g_i_)
          enddo
          c6 = d11_v
C--------
          d4_v = alb / alhh
          d5_v = d4_v - 1.0d0
          d8_v = 1.0d0 + d5_v * dexp2
          d6_b = c6 * d5_v
          d7_b = c6 * dexp2
          d8_b = d7_b * (1.0d0 / alhh)
          d9_b = d7_b * ((-d4_v) / alhh)
          do g_i_ = 1, g_p_
            g_c6(g_i_) = d6_b * g_dexp2(g_i_) + d9_b * g_alhh(g_i_) + d8
     *_b * g_alb(g_i_) + d8_v * g_c6(g_i_)
          enddo
          c6 = c6 * d8_v
C--------
C
          d4_v = rnahie1 + hie
          d5_v = (-1.5d0) * rnahie1 * hie / d4_v
          d8_v = d5_v * rnahpol1 * hpol
          d10_v = d8_v * rnag1c
          d5_b = div1 * d8_v
          d6_b = div1 * rnag1c * hpol
          d7_b = d6_b * rnahpol1
          d8_b = d6_b * d5_v
          d11_b = d7_b * ((-d5_v) / d4_v) + d7_b * (1.0d0 / d4_v) * hie 
     ** (-1.5d0)
          do g_i_ = 1, g_p_
            g_c61(g_i_) = d10_v * g_div1(g_i_) + d5_b * g_rnag1c(g_i_) +
     * d8_b * g_rnahpol1(g_i_) + d11_b * g_rnahie1(g_i_)
          enddo
          c61 = d10_v * div1
C--------
          d2_v = (-39.0625d0) * rnag1c
          d4_b = rnag1 * (-39.0625d0)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_v * g_rnag1(g_i_) + d4_b * g_rnag1c(g_i_)
          enddo
          d1_w = d2_v * rnag1
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_dexp1(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          dexp1 = d2_v
C--------
          d4_v = rnahpol1b / rnahpol1
          d5_v = d4_v - 1.0d0
          d8_v = 1.0d0 + d5_v * dexp1
          d6_b = c61 * d5_v
          d7_b = c61 * dexp1
          d8_b = d7_b * (1.0d0 / rnahpol1)
          d9_b = d7_b * ((-d4_v) / rnahpol1)
          do g_i_ = 1, g_p_
            g_c61(g_i_) = d6_b * g_dexp1(g_i_) + d9_b * g_rnahpol1(g_i_)
     * + d8_b * g_rnahpol1b(g_i_) + d8_v * g_c61(g_i_)
          enddo
          c61 = c61 * d8_v
C--------
C
          d4_v = rnahie3 + hie
          d5_v = (-1.5d0) * rnahie3 * hie / d4_v
          d8_v = d5_v * rnahpol3 * hpol
          d10_v = d8_v * rnag3c
          d5_b = div3 * d8_v
          d6_b = div3 * rnag3c * hpol
          d7_b = d6_b * rnahpol3
          d8_b = d6_b * d5_v
          d11_b = d7_b * ((-d5_v) / d4_v) + d7_b * (1.0d0 / d4_v) * hie 
     ** (-1.5d0)
          do g_i_ = 1, g_p_
            g_c63(g_i_) = d10_v * g_div3(g_i_) + d5_b * g_rnag3c(g_i_) +
     * d8_b * g_rnahpol3(g_i_) + d11_b * g_rnahie3(g_i_)
          enddo
          c63 = d10_v * div3
C--------
          d2_v = (-39.0625d0) * rnag3c
          d4_b = rnag3 * (-39.0625d0)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d2_v * g_rnag3(g_i_) + d4_b * g_rnag3c(g_i_)
          enddo
          d1_w = d2_v * rnag3
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_dexp3(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          dexp3 = d2_v
C--------
          d4_v = rnahpol3b / rnahpol3
          d5_v = d4_v - 1.0d0
          d8_v = 1.0d0 + d5_v * dexp3
          d6_b = c63 * d5_v
          d7_b = c63 * dexp3
          d8_b = d7_b * (1.0d0 / rnahpol3)
          d9_b = d7_b * ((-d4_v) / rnahpol3)
          do g_i_ = 1, g_p_
            g_c63(g_i_) = d6_b * g_dexp3(g_i_) + d9_b * g_rnahpol3(g_i_)
     * + d8_b * g_rnahpol3b(g_i_) + d8_v * g_c63(g_i_)
          enddo
          c63 = c63 * d8_v
C--------
C
          d4_v = 1.5d0 * cstsq1 + 0.5d0
          d5_v = 1.0d0 / d4_v
          d6_v = d5_v - 1.0d0
          d9_v = 1.0d0 + d6_v * dexp1
          d6_b = rmuimu1 * d6_v
          d10_b = rmuimu1 * dexp1 * ((-d5_v) / d4_v) * 1.5d0
          do g_i_ = 1, g_p_
            g_cind1(g_i_) = d6_b * g_dexp1(g_i_) + d10_b * g_cstsq1(g_i_
     *) + d9_v * g_rmuimu1(g_i_)
          enddo
          cind1 = rmuimu1 * d9_v
C--------
C
          d4_v = 1.5d0 * cstsq3 + 0.5d0
          d5_v = 1.0d0 / d4_v
          d6_v = d5_v - 1.0d0
          d9_v = 1.0d0 + d6_v * dexp3
          d6_b = rmuimu3 * d6_v
          d10_b = rmuimu3 * dexp3 * ((-d5_v) / d4_v) * 1.5d0
          do g_i_ = 1, g_p_
            g_cind3(g_i_) = d6_b * g_dexp3(g_i_) + d10_b * g_cstsq3(g_i_
     *) + d9_v * g_rmuimu3(g_i_)
          enddo
          cind3 = rmuimu3 * d9_v
C--------
C
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = balp * g_r(g_i_, i, 2)
          enddo
          d1_w = balp * (r(i, 2) - brho)
          d7_v = (qq + qh + c6 * cevau) * 0.5d0
          d9_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d9_v *  d9_v)
          d10_v = 1.0d0 - d9_v
          d5_b = (-d7_v) * d1_p
          d6_b = d10_v * 0.5d0
          d9_b = d6_b * cevau
          do g_i_ = 1, g_p_
            g_elec(g_i_) = d5_b * g_d1_w(g_i_) + d9_b * g_c6(g_i_) + d6_
     *b * g_qh(g_i_) + d6_b * g_qq(g_i_)
          enddo
          elec = d7_v * d10_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dalp * g_r(g_i_, i, 1)
          enddo
          d1_w = dalp * (r(i, 1) - drho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_c6swt1(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          c6swt1 = 0.5d0 * (1.0d0 - d2_v)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dalp * g_r(g_i_, i, 3)
          enddo
          d1_w = dalp * (r(i, 3) - drho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_c6swt3(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          c6swt3 = 0.5d0 * (1.0d0 - d2_v)
C--------
C
C
          d3_b = cevau * c6swt1
          d4_b = cevau * c61
          do g_i_ = 1, g_p_
            g_c61(g_i_) = d4_b * g_c6swt1(g_i_) + d3_b * g_c61(g_i_)
          enddo
          c61 = c61 * c6swt1 * cevau
C--------
          d3_b = cevau * c6swt3
          d4_b = cevau * c63
          do g_i_ = 1, g_p_
            g_c63(g_i_) = d4_b * g_c6swt3(g_i_) + d3_b * g_c63(g_i_)
          enddo
          c63 = c63 * c6swt3 * cevau
C--------
          do g_i_ = 1, g_p_
            g_cind1(g_i_) = cind1 * g_c6swt1(g_i_) + c6swt1 * g_cind1(g_
     *i_)
          enddo
          cind1 = cind1 * c6swt1
C--------
          do g_i_ = 1, g_p_
            g_cind3(g_i_) = cind3 * g_c6swt3(g_i_) + c6swt3 * g_cind3(g_
     *i_)
          enddo
          cind3 = cind3 * c6swt3
C--------
C
          do g_i_ = 1, g_p_
            g_elec(g_i_) = g_cind3(g_i_) + g_cind1(g_i_) + g_c63(g_i_) +
     * g_c61(g_i_) + g_elec(g_i_)
          enddo
          elec = elec + c61 + c63 + cind1 + cind3
C--------
C ===============================================================
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
          d8_v = sqrt(d1_w)

          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d8_v)
          else
             d1_p = 0.0d0
          endif
          d5_b = (-(1.0d0 / dsqrt(2.0d0))) * d1_p
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = d5_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2
     *(g_i_) + g_coul1(g_i_)
          enddo
          e(i) = coul1 + coul2 + coul3 + deh2 - d8_v / dsqrt(2.0d0)
C--------
C
          d2_b = 1.0d0 / cevau
          do g_i_ = 1, g_p_
            g_elec(g_i_) = d2_b * g_elec(g_i_)
          enddo
          elec = elec / cevau
C--------
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = cevau * g_elec(g_i_) + cevau * g_e(g_i_, i)
          enddo
          e(i) = (e(i) + elec) * cevau
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
