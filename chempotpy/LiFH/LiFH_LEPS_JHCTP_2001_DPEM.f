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
C   Common name:     LiFH surface fit H
C   Number of derivatives: 1
C   Number of electronic surfaces: 2
C   Interface: 3-2V

C   References:      A. W. Jasper, M. D. Hack, A. Chakraborty, D. G. Truhlar, 
C                    and P. Piecuch, J. Chem. Phys. 115, 7945 (2001).

c   Notes:        This fit is qualitatively correct.
c                 It is quantitatively incorrect near the saddle point.
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
C   This file was generated on 07/11/00 by the version of
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
C
        integer i, nt
C
C misc variables
        double precision r1s, r2s, r3s
C
C LiH variables
        double precision s1, relih, bes1, s1c2, s1c3, s1a, xelih, bes1b,
     * s1c
        double precision rcut2del, rcut2rho, rcut2, r1r3, r1r3al
C
C HF variables
        double precision s2, behf, yhf, rehf, dehf, xehf
C
C LiF variables
        double precision s3, res3, bes3, des3, res3c, bes3c, s3c, s3a, x
     *es3, des3c
        double precision rcutdel, rcutrho, rcut, r1r3al2
C
        double precision cevau
C
        double precision r(nt, 3), e(nt)
C
C! values of variables:
C LiH singlet:
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d1_p, d1_w, d4_v, d5_b, d4_b, d3_b, d2_b, d2_v,
     * d3_v, g_r1r3(g_pmax_)
        double precision g_r(ldg_r, nt, 3), g_d1_w(g_pmax_), g_rcut2(g_p
     *max_), g_rcut(g_pmax_), g_xelih(g_pmax_), g_s1a(g_pmax_), g_s1c(g_
     *pmax_), g_s1(g_pmax_), g_yhf(g_pmax_), g_behf(g_pmax_)
        double precision g_xehf(g_pmax_), g_s2(g_pmax_), g_xes3(g_pmax_)
     *, g_s3a(g_pmax_), g_s3c(g_pmax_), g_s3(g_pmax_), g_e(ldg_e, nt)
        integer g_ehfid
        save g_xehf, g_s2, g_xes3, g_s3a, g_s3c, g_s3
        save g_r1r3, g_d1_w, g_rcut2, g_rcut, g_xelih, g_s1a, g_s1c, g_s
     *1, g_yhf, g_behf
        data relih /1.53333d0/, bes1 /1.73333d0/, s1c2 /12.90323d0/
        data s1c3 /7.09677d0/, bes1b /1.36d0/
        data rcut2rho /1.07143d0/, rcut2del /0.6d0/, r1r3al /0.4d0/
C
C HF singlet
        data rehf /1.733d0/, dehf /6.122d0/
C
C LiF singlet
        data res3 /3.66667d0/, bes3 /1.22667d0/, des3 /0.14286d0/
        data res3c /2.48571d0/, bes3c /1.92857d0/, des3c /0.25d0/
        data rcutrho /0.4d0/, rcutdel /0.89333d0/, r1r3al2 /0.7333d0/
C
        data cevau /0.036749309d0/
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
        return
C
        entry g_pot11(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
C
        do i = 1, nt
C
          r1s = r(i, 1) ** 2
          r2s = r(i, 2) ** 2
          r3s = r(i, 3) ** 2
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
          do g_i_ = 1, g_p_
            g_r1r3(g_i_) = r1r3al2 * g_r(g_i_, i, 2) + (-g_r(g_i_, i, 3)
     *) + g_r(g_i_, i, 1)
          enddo
          r1r3 = r(i, 1) - r(i, 3) + r1r3al2 * r(i, 2)
C--------
          d3_v = (r1r3 - rcutrho) * (r1r3 - rcutrho)
          d1_p = 2.0d0 * (r1r3 - rcutrho)
          d5_b = (-(1.0d0 / rcutdel ** 2)) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_r1r3(g_i_)
          enddo
          d1_w = (-d3_v) / rcutdel ** 2
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_rcut(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          rcut = d2_v
C--------
C
C
C==============================================================k
CThis section contains the diatomic curves
C
C The  MH singlet
C
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
          d3_v = s1c - s1a
          d2_b = 1.0d0 + (-rcut2)
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d3_v * g_rcut2(g_i_) + rcut2 * g_s1c(g_i_) + d2
     *_b * g_s1a(g_i_)
          enddo
          s1 = s1a + d3_v * rcut2
C--------
C
C The HF singlet
C
          rehf = 1.733d0
          dehf = 6.122d0
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
          d2_v = des3 * xes3
          d3_v = xes3 - 2.0d0
          d4_b = d2_v + d3_v * des3
          do g_i_ = 1, g_p_
            g_s3a(g_i_) = d4_b * g_xes3(g_i_)
          enddo
          s3a = d2_v * d3_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bes3c) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-bes3c) * (r(i, 3) - res3c)
          d2_v = exp(d1_w)
          d1_p =  d2_v
          do g_i_ = 1, g_p_
            g_xes3(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xes3 = d2_v
C--------
          d2_v = des3c * xes3
          d3_v = xes3 - 2.0d0
          d4_b = d2_v + d3_v * des3c
          do g_i_ = 1, g_p_
            g_s3c(g_i_) = d4_b * g_xes3(g_i_)
          enddo
          s3c = d2_v * d3_v
C--------
C
          d3_v = s3c - s3a
          d2_b = 1.0d0 + (-rcut)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d3_v * g_rcut(g_i_) + rcut * g_s3c(g_i_) + d2_b
     * * g_s3a(g_i_)
          enddo
          s3 = s3a + d3_v * rcut
C--------
C
C ===============================================================
C The final form
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
        enddo
        return
C
      end
C
C ********************************************************************
C ********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 07/11/00 by the version of
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
        integer i, nt
        double precision r(nt, 3), e(nt)
C
C LiH variables
        integer ng1
        double precision g1, g1r, g1a, h1, h1del, h1rho
C
C HF variables
        double precision j2p, j2, j2pdel, j2rho, j2prho, j2del
C
C LiF variables
        integer ng3
        double precision g3, g3r, g3a, h3, h3del, h3rho
C
C cutoff variables
        double precision gcrho, theta, cevau, gcdel
        double precision cut2, gc
C
C LiH data
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d13_b, d12_b, d5_b, d9_v, d7_b, d8_b, d1_w, d1_
     *p, d2_v, d3_v
        double precision d4_v, d2_p, d2_b, d3_b, d4_b, d5_v, d6_v, d11_v
     *, g_d1_w(g_pmax_), g_r(ldg_r, nt, 3)
        double precision g_g1(g_pmax_), g_g3(g_pmax_), g_j2(g_pmax_), g_
     *j2p(g_pmax_), g_h1(g_pmax_), g_h3(g_pmax_), g_e(ldg_e, nt), g_gc(g
     *_pmax_), g_cut2(g_pmax_)
        integer g_ehfid
        save g_d1_w, g_g1, g_g3, g_j2, g_j2p, g_h1, g_h3, g_gc, g_cut2
        data g1a /1.27742d0/, g1r /2.5873d0/, ng1 /6/
        data h1rho /4.27097d0/, h1del /2.5d0/
C
C LiF data
        data g3a /0.48d0/, g3r /3.47619d0/, ng3 /8/
        data h3rho /2.17742d0/, h3del /0.56452d0/
C
C HF data
        data j2rho /1.15484d0/, j2del /1.75806d0/
        data j2prho /1.45161d0/, j2pdel /0.98387d0/
C
C cutoff data
        data theta /0.0d0/, gcrho /3.87097d0/, gcdel /0.45806d0/
C
        data cevau /0.036749309d0/
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
        return
C
        entry g_pot12(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
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
          d5_b = -dsin(theta)
          d7_b = dcos(theta)
          do g_i_ = 1, g_p_
            g_gc(g_i_) = d5_b * g_r(g_i_, i, 3) + d7_b * g_r(g_i_, i, 2)
          enddo
          gc = (r(i, 2) - gcrho) * dcos(theta) - (r(i, 3) - 2.9553d0) * 
     *dsin(theta)
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
C ********************************************************************
C ********************************************************************
C                           DISCLAIMER
C
C   This file was generated on 07/11/00 by the version of
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
C
C misc variables
        integer i, j, k, nt
        double precision r(nt, 3), e(nt), r1s, r2s, r3s, cs, cevau, uli2
     *p
        double precision c2a, c2b, c2c, w, cplg2
C
C LiH variables
        double precision r1r3, xelih, s1c, bes1, rcut2del, rcut2rho, s1c
     *3
        double precision s1c1, belihc, bet1, relih, s1c2, t1c1, t1c2
        double precision r1r3al, rcut2
        double precision coul1, exch1, s1, t1
C
C HF variables
        double precision behf, xehfa, dehfm, yhf, xehfa_180, t2a3, t2a2,
     * bet2
        double precision t2a1, t2a4, dehf, rehf
        double precision s2alp, rehft, bet2_180, s2rho, dehfc
        double precision coul2, exch2, t2, s2
C
C LiF variables
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
        double precision coul3, exch3, s3, t3
C
C
C LiH singlets and triplets
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d3_p, d10_v, d17_b, d11_v, d15_v, d14_b, d10_b,
     * d2_v, d16_v, d2_b
        double precision d1_p, d3_v, d4_v, d5_v, d3_b, d4_b, d5_b, d6_v,
     * d6_b, d2_p
        double precision d1_w, d7_v, d8_v, d9_v, d7_b, d8_b, d9_b, g_r1s
     *(g_pmax_), g_r(ldg_r, nt, 3), g_r2s(g_pmax_)
        double precision g_r3s(g_pmax_), g_cs(g_pmax_), g_d1_w(g_pmax_),
     * g_r1r3(g_pmax_), g_rcut2(g_pmax_), g_xelih(g_pmax_), g_s1(g_pmax_
     *), g_t1(g_pmax_), g_s1c(g_pmax_), g_yhf(g_pmax_)
        double precision g_behf(g_pmax_), g_xehfa(g_pmax_), g_dehfm(g_pm
     *ax_), g_s2(g_pmax_), g_xehfa_180(g_pmax_), g_t2(g_pmax_), g_belif(
     *g_pmax_), g_gswt(g_pmax_), g_beliff(g_pmax_), g_belifmod(g_pmax_)
        double precision g_xelif(g_pmax_), g_s3a(g_pmax_), g_delifm(g_pm
     *ax_), g_belifm(g_pmax_), g_s3c(g_pmax_), g_xelif_180(g_pmax_), g_t
     *3(g_pmax_), g_rho2(g_pmax_), g_theta2(g_pmax_), g_del2(g_pmax_)
        double precision g_gc2(g_pmax_), g_cut2(g_pmax_), g_s3(g_pmax_),
     * g_coul1(g_pmax_), g_coul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_p
     *max_), g_exch2(g_pmax_), g_exch3(g_pmax_), g_w(g_pmax_)
        double precision g_cplg2(g_pmax_), g_e(ldg_e, nt)
        integer g_ehfid
        save g_exch3, g_w, g_cplg2
        save g_theta2, g_del2, g_gc2, g_cut2, g_s3, g_coul1, g_coul2, g_
     *coul3, g_exch1, g_exch2
        save g_beliff, g_belifmod, g_xelif, g_s3a, g_delifm, g_belifm, g
     *_s3c, g_xelif_180, g_t3, g_rho2
        save g_s1c, g_yhf, g_behf, g_xehfa, g_dehfm, g_s2, g_xehfa_180, 
     *g_t2, g_belif, g_gswt
        save g_r1s, g_r2s, g_r3s, g_cs, g_d1_w, g_r1r3, g_rcut2, g_xelih
     *, g_s1, g_t1
        intrinsic dble
        data s1c1 /4.32258d0/, s1c2 /7.06452d0/, t1c1 /1.64516d0/
        data t1c2 /10.38710d0/, belihc /1.36d0/, bet1 /2.10667d0/
        data relih /1.2d0/, r1r3al /1.0d0/, rcut2rho /0.72d0/
        data rcut2del /0.5d0/, s1c3 /14.74194d0/, bes1 /0.90667d0/
C
C HF singlet & triplet
        data rehf /1.733d0/, dehf /6.122d0/
        data t2a1 /1.26667d0/, t2a2 /16.06667d0/, t2a3 /11.61290d0/
        data t2a4 /15.51613d0/, bet2 /2.10667d0/, bet2_180 /1.73333d0/
        data rehft /1.733d0/, dehfc /0.26d0/, s2rho /2.38095d0/
        data s2alp /0.5d0/
C
C LiF singlets and triplets
        data belifc1 /0.95484d0/, belifc2 /0.91613d0/, delifc1 /2.58065d
     *0/
        data delifc2 /5.25806d0/, belifff /0.25333d0/, nb /8/
        data gam /3.87097d0/, grho /1.41935d0/, gdel /3.77419d0/
        data t3a1 /0.38710d0/, t3a2 /1.80645d0/, t3a3 /0.51613d0/
        data t3a4 /0.51613d0/, bet3 /0.89333d0/, bet3_180 /0.80d0/
        data rho2_0 /0.51613d0/, rho2_180 /1.33333d0/, del2_0 /0.49032d0
     */
        data del2_180 /0.31613d0/, theta2_0 /0.09333d0/
        data theta2_180 /0.18d0/, delif /5.909d0/, relif /2.9553d0/
C
C misc constants
C
        data c2a /3.5d0/, c2b /0.27362d0/, c2c /0.15d0/
        data uli2p /1.848d0/, cevau /0.036749309d0/
C
        save
        data g_ehfid /0/
C
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        return
C
        entry g_pot22(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)
C
        do 99999 i = 1, nt
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

          if (cs .le.  1.0d0) then
             d1_p = 1.0d0
             d2_p = 0.0d0
          else if (cs .gt.  1.0d0) then
             d1_p = 0.0d0
             d2_p = 1.0d0
          endif
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

          if (cs .ge.  d1_w) then
             d1_p = 1.0d0
             d2_p = 0.0d0
          else if (cs .lt.  d1_w) then
             d1_p = 0.0d0
             d2_p = 1.0d0
          endif
          do g_i_ = 1, g_p_
            g_cs(g_i_) = d2_p * g_d1_w(g_i_) + d1_p * g_cs(g_i_)
          enddo
          cs = d3_v
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
C =========================================================
C The diatomic curves . . .
C
C The MH singlet & triplet
C
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
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-bet1) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-bet1) * (r(i, 1) - relih)
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
            g_s1(g_i_) = rcut2 * g_s1c(g_i_) + s1c * g_rcut2(g_i_) + g_s
     *1(g_i_)
          enddo
          s1 = s1 + rcut2 * s1c
C--------
C
C The HF singlet and triplet
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
            g_xehfa(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          xehfa = d2_v
C--------
C
          d3_b = 1.0d0 / s2alp
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_r(g_i_, i, 3)
          enddo
          d1_w = (r(i, 3) - s2rho) / s2alp
          d4_v = dehfc * 0.5d0 * (1.0d0 - cs) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.0d0 - d6_v
          d6_b = d4_v * d1_p
          d9_b = -((-d7_v) * 0.5d0 * (dehfc * 0.5d0))
          do g_i_ = 1, g_p_
            g_dehfm(g_i_) = d6_b * g_d1_w(g_i_) + d9_b * g_cs(g_i_)
          enddo
          dehfm = dehf - d4_v * d7_v
C--------
C
          d3_v = dehfm * xehfa
          d4_v = xehfa - 2.0d0
          d5_b = d4_v * xehfa
          d4_b = d3_v + d4_v * dehfm
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d4_b * g_xehfa(g_i_) + d5_b * g_dehfm(g_i_)
          enddo
          s2 = d3_v * d4_v
C--------
C
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
C The Modified HF singlet 
C
C The  LiF singlet & triplet
          d5_v = (4.6498d0 + 7.0489d0 * r(i, 3)) * (4.6498d0 + 7.0489d0 
     +* r(i, 3))
          d1_p = 2.0d0 * (4.6498d0 + 7.0489d0 * r(i, 3))
          d6_v = r(i, 3) * 103.57d0 / d5_v
          d7_b = (-d6_v) / d5_v * d1_p * 7.0489d0 + 1.0d0 / d5_v * 103.5
     *7d0
          do g_i_ = 1, g_p_
            g_belif(g_i_) = d7_b * g_r(g_i_, i, 3)
          enddo
          belif = d6_v + dble(0.064076)
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
          d4_b = d2_v + d3_v * (delif + uli2p)
          do g_i_ = 1, g_p_
            g_s3a(g_i_) = d4_b * g_xelif(g_i_)
          enddo
          s3a = d2_v * d3_v
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
C
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
             d1_p = 0.0d0 ! should never occur
          endif
          d7_b = (-(1.0d0 / dsqrt(2.0d0))) * d1_p
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = d7_b * g_d1_w(g_i_) + g_coul3(g_i_) + g_coul2
     *(g_i_) + g_coul1(g_i_)
          enddo
          e(i) = coul1 + coul2 + coul3 - d7_v / dsqrt(2.0d0) + uli2p + d
     *ehf
C--------
C
C
C ===============================================================
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
