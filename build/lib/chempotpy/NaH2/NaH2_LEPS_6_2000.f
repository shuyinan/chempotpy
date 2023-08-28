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

C
C   System:          NaH2
C   Functional form: Modified Extended LEPS (2x2 diabatic fit)
C   Common name:     NaH2 surface set 6
C   Interface:       3-2V
C   Number of electronic surfaces: 2
C   Number of derivatives: 1
C   References:      M. D. Hack and D. G. Truhlar, J. Chem. Phys. 
C                    112, 9716 (2000)
C
C   Notes:   This is the most accurate NaH2 fit.
C
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
C      Ground diabatic electronic state, first-excited diabatic electronic state and the potential coupling.  
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

c  The final version of surface 6, completed on 9/24/98
c  This version contains derivatives, to be used with trajectory
c  and other codes . . .

c  It is based on version 5g of the NaH2 potential, in turn based on 5f.
c  Reference:  P. Halvick and D. G. Truhlar, J. Chem. Phys. 96,2895 (1992).  
c  The "zero" of energy occurs when the Na atom is "infinitely"
c  far from the H2 diatom and the H2 diatom is at its equilibrium geometry.
c
c  The potential is called by first issuing a call to PREPOT which initializes
c  all surfaces.  The potential is called by calling POT with the parameters
c  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of
c  dimension nt which is an integer representing the number of geometries
c  input.  The integer nsurf selects the desired surface.  nsurf=1 is the
c  lower diabatic surface, nsurf=2 is the coupling surface, and nsurf=3 is
c  the upper diabatic surface. 
c

c ********************************************************************
c ********************************************************************

      subroutine prepot
        implicit none
        integer g_p_
        parameter(g_p_=3)
        integer g_myrank,g_nprocs
        common /com_para/ g_myrank,g_nprocs

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

C

        if (nsurf .eq. 1) then
          call g_pot11(g_p_, r, g_r, ldg_r, e, De, ldg_e, nt)
        else if (nsurf .eq. 2) then
          call g_pot12(g_p_, r, g_r, ldg_r, e, De, ldg_e, nt)
        else if (nsurf .eq. 3) then
          call g_pot22(g_p_, r, g_r, ldg_r, e, De, ldg_e, nt)
        else
          do 10 i=1,nt
            call g_pot11(g_p_, r, g_r, ldg_r, v11, De, ldg_e, 1)
            call g_pot12(g_p_, r, g_r, ldg_r, v12, De, ldg_e, 1)
            call g_pot22(g_p_, r, g_r, ldg_r, v22, De, ldg_e, 1)
            if (nsurf .eq. 4) then
              e(i)=0.5d0*(v11+v22-dsqrt((v11-v22)**2+4.d0*v12**2))
            else
              e(i)=0.5d0*(v11+v22+dsqrt((v11-v22)**2+4.d0*v12**2))
            end if
   10     continue
        end if

      return

      end


C                           DISCLAIMER
C
C   This file was generated on 07/24/98 by the version of
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
        integer i,j,k,nsurf,nt

C
        save
C
C data blocks for the NaH singlets . . .
C
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d3_p, d24_b, d2_p, d15_v, d20_v, d16_v, d22_v, 
     *d2_v, d16_b, d2_b
        double precision d1_w, d2_w, d3_v, d4_v, d5_v, d6_v, d3_b, d4_b,
     * d5_b, d6_b
        double precision d1_p, d7_v, d8_v, d9_v, d10_v, d7_b, d8_b, d9_b
     *, d10_b, d11_v
        double precision d12_v, d13_v, d19_v, d11_b, d12_b, d13_b, d14_b
     *, g_r1s(g_pmax_), g_r(ldg_r, nt, 3), g_r2s(g_pmax_)
        double precision g_r3s(g_pmax_), g_d2_w(g_pmax_), g_d1_w(g_pmax_
     *), g_rnag(g_pmax_), g_rnags(g_pmax_), g_rnagc(g_pmax_), g_rnag1(g_
     *pmax_), g_rnag3(g_pmax_), g_rnag1c(g_pmax_), g_rnag3c(g_pmax_)
        double precision g_rsm1(g_pmax_), g_rsm2(g_pmax_), g_rsm3(g_pmax
     *_), g_cstsq(g_pmax_), g_cosa(g_pmax_), g_coscos(g_pmax_), g_cstsq1
     *(g_pmax_), g_cstsq3(g_pmax_), g_s1a(g_pmax_), g_s10c(g_pmax_)
        double precision g_s10(g_pmax_), g_s1ic(g_pmax_), g_s1i(g_pmax_)
     *, g_s190c(g_pmax_), g_s190(g_pmax_), g_s1(g_pmax_), g_t190(g_pmax_
     *), g_t10(g_pmax_), g_t10alt(g_pmax_), g_t1i(g_pmax_)
        double precision g_t1(g_pmax_), g_s2a(g_pmax_), g_s290c(g_pmax_)
     *, g_s290(g_pmax_), g_s20c(g_pmax_), g_s20(g_pmax_), g_s2(g_pmax_),
     * g_t2a(g_pmax_), g_t20c(g_pmax_), g_t20(g_pmax_)
        double precision g_t2(g_pmax_), g_s3a(g_pmax_), g_s30c(g_pmax_),
     * g_s30(g_pmax_), g_s3ic(g_pmax_), g_s3i(g_pmax_), g_s390c(g_pmax_)
     *, g_s390(g_pmax_), g_s3(g_pmax_), g_t390(g_pmax_)
        double precision g_t30(g_pmax_), g_t30alt(g_pmax_), g_t3(g_pmax_
     *), g_t3i(g_pmax_), g_h2ie(g_pmax_), g_al1(g_pmax_), g_al2(g_pmax_)
     *, g_rnahie1(g_pmax_), g_rnahie3(g_pmax_), g_rnahmu1(g_pmax_)
        double precision g_rnahmu3(g_pmax_), g_rnahpar1(g_pmax_), g_rnah
     *par3(g_pmax_), g_rnahperp1(g_pmax_), g_rnahperp3(g_pmax_), g_rnahp
     *ol1(g_pmax_), g_rnahpol1b(g_pmax_), g_rnahpol3(g_pmax_), g_rnahpol
     *3b(g_pmax_), g_div1(g_pmax_)
        double precision g_rmuimu1(g_pmax_), g_div3(g_pmax_), g_rmuimu3(
     *g_pmax_), g_c690(g_pmax_), g_c60(g_pmax_), g_c6(g_pmax_), g_c61(g_
     *pmax_), g_dexp1(g_pmax_), g_c63(g_pmax_), g_dexp3(g_pmax_)
        double precision g_cind1(g_pmax_), g_cind3(g_pmax_), g_c6swt1(g_
     *pmax_), g_c6swt3(g_pmax_), g_elec(g_pmax_), g_coul1(g_pmax_), g_co
     *ul2(g_pmax_), g_coul3(g_pmax_), g_exch1(g_pmax_), g_exch2(g_pmax_)
        double precision g_exch3(g_pmax_), g_w(g_pmax_), g_cplg2(g_pmax_
     *), g_e(ldg_e, nt)
        integer g_ehfid
        save g_exch1, g_exch2, g_exch3, g_w, g_cplg2
        save g_c63, g_dexp3, g_cind1, g_cind3, g_c6swt1, g_c6swt3, g_ele
     *c, g_coul1, g_coul2, g_coul3
        save g_rnahpol3b, g_div1, g_rmuimu1, g_div3, g_rmuimu3, g_c690, 
     *g_c60, g_c6, g_c61, g_dexp1
        save g_rnahie3, g_rnahmu1, g_rnahmu3, g_rnahpar1, g_rnahpar3, g_
     *rnahperp1, g_rnahperp3, g_rnahpol1, g_rnahpol1b, g_rnahpol3
        save g_s3, g_t390, g_t30, g_t30alt, g_t3, g_t3i, g_h2ie, g_al1, 
     *g_al2, g_rnahie1
        save g_t20c, g_t20, g_t2, g_s3a, g_s30c, g_s30, g_s3ic, g_s3i, g
     *_s390c, g_s390
        save g_t10alt, g_t1i, g_t1, g_s2a, g_s290c, g_s290, g_s20c, g_s2
     *0, g_s2, g_t2a
        save g_s1a, g_s10c, g_s10, g_s1ic, g_s1i, g_s190c, g_s190, g_s1,
     * g_t190, g_t10
        save g_rnag1c, g_rnag3c, g_rsm1, g_rsm2, g_rsm3, g_cstsq, g_cosa
     *, g_coscos, g_cstsq1, g_cstsq3
        save g_r1s, g_r2s, g_r3s, g_d2_w, g_d1_w, g_rnag, g_rnags, g_rna
     *gc, g_rnag1, g_rnag3
        intrinsic dble
        data s1c1 /593.88235d0/, s1c2 /1.72588d0/, s1c3 /-251.80392d0/, 
     *s1c4 /89.98431d0/, s1c5 /-8.83922d0/, s1c6 /1.383d0/
        data s10c1 /514.17323d0/, s10c2 /1.76693d0/, s10c3 /-143.89764d0
     */, s10c4 /13.76378d0/, s10c5 /-2.11024d0/, s10c6 /1.383d0/, s10alp
     * /0.79843d0/, s10rho /2.47244d0/
        data s190c1 /274.97638d0/, s190c2 /2.20787d0/, s190c3 /-28.97638
     *d0/, s190c6 /1.383d0/, s190alp /0.66693d0/, s190rho /6.95276d0/
        data s101i /479.85039d0/, s102i /2.22598d0/, s103i /-24.20472d0/
     *, s104i /30.19685d0/, s105i /-9.96063d0/, s106i /1.383d0/, s1ialp 
     */0.15d0/
C
C data blocks for the NaH triplets
        data t10c1 /335.59055d0/, t10c2 /4.02126d0/, t10c3 /78.81890d0/,
     * t10c4 /-27.51969d0/, t10c5 /4.22835d0/, t10c6 /1.383d0/
        data t10ac1 /365.90551d0/, t10ac2 /2.46220d0/, t10ac3 /49.75591d
     *0/, t10ac4 /-50.11024d0/, t10ac5 /17.32283d0/, t10ac6 /1.383d0/, t
     *10alp /0.77874d0/, t10rho /4.46457d0/
        data t190c1 /511.51181d0/, t190c2 /1.91024d0/, t190c3 /-8.53543d
     *0/, t190c4 /-1.58268d0/, t190c5 /1.53543d0/, t190c6 /1.383d0/
        data t1ic1 /430.31496d0/, t1ic2 /2.05197d0/, t1ic3 /24.40945d0/,
     * t1ic4 /12.91339d0/, t1ic5 /-6.07087d0/, t1ic6 /1.383d0/
C
C data blocks for the H2 singlets
        data s2c1 /139.7160d0/, s2c2 /-123.8978d0/, s2c3 /3.4031d0/, s2c
     *4 /-6.8725d0/, s2c5 /-23.0440d0/, s2c6 /2.032d0/
        data s290c1 /131.65354d0/, s290c2 /294.56693d0/, s290c3 /5.20472
     *d0/, s290c4 /-51.69291d0/, s290c5 /2.16535d0/, s290c6 /2.032d0/, s
     *290alp /0.09843d0/, s290rho /6.8d0/
        data s20c1 /189.76378d0/, s20c2 /-158.80315d0/, s20c3 /3.60945d0
     */, s20c4 /-16.84252d0/, s20c5 /-12.59843d0/, s20c6 /2.032d0/  s20a
     *lp /0.8063d0/, s20rho /3.37008d0/
C
C data blocks for the h2 triplets
        data t2c1 /147.214d0/, t2c2 /3.915d0/, t2c3 /41.275d0/, t2c4 /10
     *.505d0/, t2c5 /4.408d0/, t2c6 /2.032d0/
        data t20c1 /654.51969d0/, t20c2 /2.45197d0/, t20c3 /-170.03937d0
     */, t20c4 /108.26772d0/, t20c5 /-17.70079d0/, t20c6 /2.032d0/  t20a
     *lp /0.05236d0/, t20rho /52.47683993d0/
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
C misc constants
        data c2a /0.67801d0/, c2b /0.83900d0/, c2c /0.15039d0/, deh2 /4.
     *74772265651191638d0/
C
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

        return
C
        entry g_pot11(g_p_, r, g_r, ldg_r, e, g_e, ldg_e, nt)

C
        do 99999 i = 1, nt
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

          rnag = dsqrt(dabs(0.5d0*(r1s+r3s-0.5d0*r2s)))

          if (rnag .gt. 0.0d0) then
            g_rnag(1)=0.5d0*r(i,1)/rnag
            g_rnag(2)=-0.25d0*r(i,2)/rnag
            g_rnag(3)=0.5d0*r(i,3)/rnag
          else
            g_rnag(1)=0.0d0
            g_rnag(2)=0.0d0
            g_rnag(3)=0.0d0
          endif

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
          frac = 0.958d0
C H to NaH center of masses

           rnag1= dsqrt(dabs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))

           if (rnag1 .gt. 0.0d0) then
             g_rnag1(1) = r(i,1)*frac*(frac-1.0d0)/rnag1
             g_rnag1(2) = r(i,2)*(1.0d0-frac)/rnag1
             g_rnag1(3) = r(i,3)*frac/rnag1
           else
             g_rnag1(1) = 0.0d0
             g_rnag1(2) = 0.0d0
             g_rnag1(3) = 0.0d0
           endif

C--------
           rnag3= dsqrt(dabs(r2s+frac*(r1s-r2s-r3s)+r3s*frac*frac))

           if (rnag3 .gt. 0.0d0) then
             g_rnag3(1) = r(i,1)*frac/rnag3
             g_rnag3(2) = r(i,2)*(1.0d0-frac)/rnag3
             g_rnag3(3) = r(i,3)*frac*(frac-1.0d0)/rnag3
           else
             g_rnag3(1) = 0.0d0
             g_rnag3(2) = 0.0d0
             g_rnag3(3) = 0.0d0
           endif

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
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-0.5d0) * g_r(g_i_, i, 2) + 0.5d0 * g_r(g_i_
     *, i, 3)
          enddo
          d1_w = 0.5d0 * (r(i, 3) - r(i, 2))
          d3_v = r(i, 2) - r(i, 3)
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d7_v = 0.5d0 + 0.5d0 * d5_v
          d8_b = d3_v * 0.5d0 * d1_p
          d2_b = 1.0d0 + (-d7_v)
          do g_i_ = 1, g_p_
            g_rsm1(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_r(g_i_, i, 2) 
     *+ d2_b * g_r(g_i_, i, 3)
          enddo
          rsm1 = r(i, 3) + d3_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-0.5d0) * g_r(g_i_, i, 1) + 0.5d0 * g_r(g_i_
     *, i, 3)
          enddo
          d1_w = 0.5d0 * (r(i, 3) - r(i, 1))
          d3_v = r(i, 1) - r(i, 3)
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d7_v = 0.5d0 + 0.5d0 * d5_v
          d8_b = d3_v * 0.5d0 * d1_p
          d2_b = 1.0d0 + (-d7_v)
          do g_i_ = 1, g_p_
            g_rsm2(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_r(g_i_, i, 1) 
     *+ d2_b * g_r(g_i_, i, 3)
          enddo
          rsm2 = r(i, 3) + d3_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-0.5d0) * g_r(g_i_, i, 2) + 0.5d0 * g_r(g_i_
     *, i, 1)
          enddo
          d1_w = 0.5d0 * (r(i, 1) - r(i, 2))
          d3_v = r(i, 2) - r(i, 1)
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d7_v = 0.5d0 + 0.5d0 * d5_v
          d8_b = d3_v * 0.5d0 * d1_p
          d2_b = 1.0d0 + (-d7_v)
          do g_i_ = 1, g_p_
            g_rsm3(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_r(g_i_, i, 2) 
     *+ d2_b * g_r(g_i_, i, 1)
          enddo
          rsm3 = r(i, 1) + d3_v * d7_v
C--------
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
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-10.0d0) * g_rnags(g_i_)
          enddo
          d1_w = (-10.0d0) * rnags
          d4_v = 0.5d0 + 0.5d0 * cosa
          d5_v = cstsq * d4_v
          d7_v = exp(d1_w)
          d1_p =  d7_v
          d8_v = 1.0d0 - d7_v
          d5_b = (-d5_v) * d1_p
          d6_b = d8_v * d4_v
          d9_b = d8_v * cstsq * 0.5d0
          do g_i_ = 1, g_p_
            g_coscos(g_i_) = d5_b * g_d1_w(g_i_) + d9_b * g_cosa(g_i_) +
     * d6_b * g_cstsq(g_i_)
          enddo
          coscos = d5_v * d8_v
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
C==============================================================k
CThis section contains the diatomic curves
C
C The first NaH singlet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s1c2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-s1c2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s1c6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-s1c6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s1c3 + s1c4 * r(i, 1) + s1c5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * s1c5
          d7_b = d12_v * d9_v + d8_b * s1c4
          d14_b = s1c1 * d2_p
          do g_i_ = 1, g_p_
            g_s1a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7
     *_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          s1a = s1c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s10c2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-s10c2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s10c6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-s10c6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s10c3 + s10c4 * r(i, 1) + s10c5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * s10c5
          d7_b = d12_v * d9_v + d8_b * s10c4
          d14_b = s10c1 * d2_p
          do g_i_ = 1, g_p_
            g_s10c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d
     *7_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          s10c = s10c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s10alp * g_rsm1(g_i_)
          enddo
          d1_w = s10alp * (rsm1 - s10rho)
          d4_v = (s10c - s1a) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.0d0 - d6_v
          d7_b = (-d4_v) * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s10(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s10c(g_i_) + d2
     *_b * g_s1a(g_i_)
          enddo
          s10 = s1a + d4_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s102i) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-s102i) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s106i) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-s106i) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s103i + s104i * r(i, 1) + s105i * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * s105i
          d7_b = d12_v * d9_v + d8_b * s104i
          d14_b = s101i * d2_p
          do g_i_ = 1, g_p_
            g_s1ic(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d
     *7_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          s1ic = s101i * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s1ialp * g_rnags(g_i_)
          enddo
          d1_w = s1ialp * rnags
          d3_v = s1ic - s1a
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d6_v = 1.0d0 - d5_v
          d7_b = (-d3_v) * d1_p
          d2_b = 1.0d0 + (-d6_v)
          do g_i_ = 1, g_p_
            g_s1i(g_i_) = d7_b * g_d1_w(g_i_) + d6_v * g_s1ic(g_i_) + d2
     *_b * g_s1a(g_i_)
          enddo
          s1i = s1a + d3_v * d6_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s190c2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-s190c2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s190c6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-s190c6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d5_v = r(i, 1) * s190c3
          d7_v = exp(d2_w)
          d1_p =  d7_v
          d6_b = d5_v * d1_p
          d7_b = d7_v * s190c3
          d9_b = s190c1 * d2_p
          do g_i_ = 1, g_p_
            g_s190c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r(g_i_, i, 1)
     * + d9_b * g_d1_w(g_i_)
          enddo
          s190c = s190c1 * d2_v + d5_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s190alp * g_r(g_i_, i, 2)
          enddo
          d1_w = s190alp * (r(i, 2) - s190rho)
          d4_v = (s190c - s1a) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.0d0 - d6_v
          d7_b = (-d4_v) * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s190(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s190c(g_i_) + 
     *d2_b * g_s1a(g_i_)
          enddo
          s190 = s1a + d4_v * d7_v
C--------
C
          d3_v = s10 - s190
          d8_v = 0.5d0 + 0.5d0 * cosa
          d13_v = 0.5d0 - 0.5d0 * cosa
          d7_b = (-s1i) * 0.5d0 + s190 * 0.5d0
          d10_b = d8_v + (-coscos)
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d13_v * g_s1i(g_i_) + d7_b * g_cosa(g_i_) + d3_
     *v * g_coscos(g_i_) + d10_b * g_s190(g_i_) + coscos * g_s10(g_i_)
          enddo
          s1 = d3_v * coscos + s190 * d8_v + s1i * d13_v
C--------
C
C
C The first NaH triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t190c2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-t190c2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t190c6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-t190c6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t190c3 + t190c4 * r(i, 1) + t190c5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * t190c5
          d7_b = d12_v * d9_v + d8_b * t190c4
          d14_b = t190c1 * d2_p
          do g_i_ = 1, g_p_
            g_t190(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d
     *7_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          t190 = t190c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t10c2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-t10c2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t10c6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-t10c6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t10c3 + t10c4 * r(i, 1) + t10c5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * t10c5
          d7_b = d12_v * d9_v + d8_b * t10c4
          d14_b = t10c1 * d2_p
          do g_i_ = 1, g_p_
            g_t10(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7
     *_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          t10 = t10c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t10ac2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-t10ac2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t10ac6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-t10ac6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t10ac3 + t10ac4 * r(i, 1) + t10ac5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * t10ac5
          d7_b = d12_v * d9_v + d8_b * t10ac4
          d14_b = t10ac1 * d2_p
          do g_i_ = 1, g_p_
            g_t10alt(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) +
     * d7_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          t10alt = t10ac1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = t10alp * g_r(g_i_, i, 2)
          enddo
          d1_w = t10alp * (r(i, 2) - t10rho)
          d3_v = t10alt - t10
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d7_v = 0.5d0 + 0.5d0 * d5_v
          d8_b = d3_v * 0.5d0 * d1_p
          d2_b = 1.0d0 + (-d7_v)
          do g_i_ = 1, g_p_
            g_t10(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_t10alt(g_i_) + 
     *d2_b * g_t10(g_i_)
          enddo
          t10 = t10 + d3_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t1ic2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-t1ic2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t1ic6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-t1ic6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t1ic3 + t1ic4 * r(i, 1) + t1ic5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * t1ic5
          d7_b = d12_v * d9_v + d8_b * t1ic4
          d14_b = t1ic1 * d2_p
          do g_i_ = 1, g_p_
            g_t1i(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7
     *_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          t1i = t1ic1 * d2_v + d10_v * d12_v
C--------
C
          d3_v = t10 - t190
          d8_v = 0.5d0 + 0.5d0 * cosa
          d13_v = 0.5d0 - 0.5d0 * cosa
          d7_b = (-t1i) * 0.5d0 + t190 * 0.5d0
          d10_b = d8_v + (-coscos)
          do g_i_ = 1, g_p_
            g_t1(g_i_) = d13_v * g_t1i(g_i_) + d7_b * g_cosa(g_i_) + d3_
     *v * g_coscos(g_i_) + d10_b * g_t190(g_i_) + coscos * g_t10(g_i_)
          enddo
          t1 = d3_v * coscos + t190 * d8_v + t1i * d13_v
C--------
C
C
C The H2 singlet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s2c3) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-s2c3) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s2c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-s2c6) * r(i, 2)
          d3_v = s2c1 + s2c2 * r(i, 2)
          d5_v = exp(d1_w)
          d2_p =  d5_v
          d9_v = s2c4 + s2c5 * r(i, 2)
          d10_v = r2s * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d7_b = d12_v * d9_v
          d13_b = d3_v * d2_p
          d10_b = d12_v * r2s * s2c5 + d5_v * s2c2
          do g_i_ = 1, g_p_
            g_s2a(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r2s(g_i_) + d13
     *_b * g_d1_w(g_i_) + d10_b * g_r(g_i_, i, 2)
          enddo
          s2a = d3_v * d5_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s290c3) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-s290c3) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s290c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-s290c6) * r(i, 2)
          d3_v = s290c1 + s290c2 * r(i, 2)
          d5_v = exp(d1_w)
          d2_p =  d5_v
          d9_v = s290c4 + s290c5 * r(i, 2)
          d10_v = r2s * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d7_b = d12_v * d9_v
          d13_b = d3_v * d2_p
          d10_b = d12_v * r2s * s290c5 + d5_v * s290c2
          do g_i_ = 1, g_p_
            g_s290c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r2s(g_i_) + d
     *13_b * g_d1_w(g_i_) + d10_b * g_r(g_i_, i, 2)
          enddo
          s290c = d3_v * d5_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s290alp * g_rnags(g_i_)
          enddo
          d1_w = s290alp * (rnags - s290rho)
          d4_v = (s290c - s2a) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.d0 - d6_v
          d7_b = (-d4_v) * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s290(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s290c(g_i_) + 
     *d2_b * g_s2a(g_i_)
          enddo
          s290 = s2a + d4_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s20c3) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-s20c3) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s20c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-s20c6) * r(i, 2)
          d3_v = s20c1 + s20c2 * r(i, 2)
          d5_v = exp(d1_w)
          d2_p =  d5_v
          d9_v = s20c4 + s20c5 * r(i, 2)
          d10_v = r2s * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d7_b = d12_v * d9_v
          d13_b = d3_v * d2_p
          d10_b = d12_v * r2s * s20c5 + d5_v * s20c2
          do g_i_ = 1, g_p_
            g_s20c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r2s(g_i_) + d1
     *3_b * g_d1_w(g_i_) + d10_b * g_r(g_i_, i, 2)
          enddo
          s20c = d3_v * d5_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s20alp * g_rsm2(g_i_)
          enddo
          d1_w = s20alp * (rsm2 - s20rho)
          d4_v = (s20c - s2a) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.d0 - d6_v
          d7_b = (-d4_v) * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s20(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s20c(g_i_) + d2
     *_b * g_s2a(g_i_)
          enddo
          s20 = s2a + d4_v * d7_v
C--------
C
          d3_v = s20 - s290
          d2_b = 1.0d0 + (-coscos)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d3_v * g_coscos(g_i_) + coscos * g_s20(g_i_) + 
     *d2_b * g_s290(g_i_)
          enddo
          s2 = s290 + d3_v * coscos
C--------
C
C
C The H2 triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t2c2) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-t2c2) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t2c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-t2c6) * r(i, 2)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t2c3 + t2c4 * r(i, 2) + t2c5 * r2s
          d10_v = r(i, 2) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 2)
          d11_b = d8_b * t2c5
          d7_b = d12_v * d9_v + d8_b * t2c4
          d14_b = t2c1 * d2_p
          do g_i_ = 1, g_p_
            g_t2a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d7
     *_b * g_r(g_i_, i, 2) + d14_b * g_d1_w(g_i_)
          enddo
          t2a = t2c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t20c2) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-t20c2) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t20c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-t20c6) * r(i, 2)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t20c3 + t20c4 * r(i, 2) + t20c5 * r2s
          d10_v = r(i, 2) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 2)
          d11_b = d8_b * t20c5
          d7_b = d12_v * d9_v + d8_b * t20c4
          d14_b = t20c1 * d2_p
          do g_i_ = 1, g_p_
            g_t20c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d
     *7_b * g_r(g_i_, i, 2) + d14_b * g_d1_w(g_i_)
          enddo
          t20c = t20c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = t20alp * g_rnags(g_i_)
          enddo
          d1_w = t20alp * (rnags - t20rho)
          d4_v = (t20c - t2a) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.d0 - d6_v
          d7_b = (-d4_v) * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_t20(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_t20c(g_i_) + d2
     *_b * g_t2a(g_i_)
          enddo
          t20 = t2a + d4_v * d7_v
C--------
C
          d3_v = t20 - t2a
          d2_b = 1.0d0 + (-coscos)
          do g_i_ = 1, g_p_
            g_t2(g_i_) = d3_v * g_coscos(g_i_) + coscos * g_t20(g_i_) + 
     *d2_b * g_t2a(g_i_)
          enddo
          t2 = t2a + d3_v * coscos
C--------
C
C The second NaH singlet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s1c2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-s1c2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s1c6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-s1c6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s1c3 + s1c4 * r(i, 3) + s1c5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * s1c5
          d7_b = d12_v * d9_v + d8_b * s1c4
          d14_b = s1c1 * d2_p
          do g_i_ = 1, g_p_
            g_s3a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7
     *_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          s3a = s1c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s10c2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-s10c2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s10c6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-s10c6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s10c3 + s10c4 * r(i, 3) + s10c5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * s10c5
          d7_b = d12_v * d9_v + d8_b * s10c4
          d14_b = s10c1 * d2_p
          do g_i_ = 1, g_p_
            g_s30c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d
     *7_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          s30c = s10c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s10alp * g_rsm3(g_i_)
          enddo
          d1_w = s10alp * (rsm3 - s10rho)
          d4_v = (s30c - s3a) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.0d0 - d6_v
          d7_b = (-d4_v) * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s30(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s30c(g_i_) + d2
     *_b * g_s3a(g_i_)
          enddo
          s30 = s3a + d4_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s102i) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-s102i) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s106i) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-s106i) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s103i + s104i * r(i, 3) + s105i * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * s105i
          d7_b = d12_v * d9_v + d8_b * s104i
          d14_b = s101i * d2_p
          do g_i_ = 1, g_p_
            g_s3ic(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d
     *7_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          s3ic = s101i * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s1ialp * g_rnags(g_i_)
          enddo
          d1_w = s1ialp * rnags
          d3_v = s3ic - s3a
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d6_v = 1.0d0 - d5_v
          d7_b = (-d3_v) * d1_p
          d2_b = 1.0d0 + (-d6_v)
          do g_i_ = 1, g_p_
            g_s3i(g_i_) = d7_b * g_d1_w(g_i_) + d6_v * g_s3ic(g_i_) + d2
     *_b * g_s3a(g_i_)
          enddo
          s3i = s3a + d3_v * d6_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s190c2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-s190c2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s190c6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-s190c6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d5_v = r(i, 3) * s190c3
          d7_v = exp(d2_w)
          d1_p =  d7_v
          d6_b = d5_v * d1_p
          d7_b = d7_v * s190c3
          d9_b = s190c1 * d2_p
          do g_i_ = 1, g_p_
            g_s390c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r(g_i_, i, 3)
     * + d9_b * g_d1_w(g_i_)
          enddo
          s390c = s190c1 * d2_v + d5_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s190alp * g_r(g_i_, i, 2)
          enddo
          d1_w = s190alp * (r(i, 2) - s190rho)
          d4_v = (s390c - s3a) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.0d0 - d6_v
          d7_b = (-d4_v) * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s390(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s390c(g_i_) + 
     *d2_b * g_s3a(g_i_)
          enddo
          s390 = s3a + d4_v * d7_v
C--------
C
          d3_v = s30 - s390
          d8_v = 0.5d0 + 0.5d0 * cosa
          d13_v = 0.5d0 - 0.5d0 * cosa
          d7_b = (-s3i) * 0.5d0 + s390 * 0.5d0
          d10_b = d8_v + (-coscos)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d13_v * g_s3i(g_i_) + d7_b * g_cosa(g_i_) + d3_
     *v * g_coscos(g_i_) + d10_b * g_s390(g_i_) + coscos * g_s30(g_i_)
          enddo
          s3 = d3_v * coscos + s390 * d8_v + s3i * d13_v
C--------
C
C
C The second NaH triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t190c2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-t190c2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t190c6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-t190c6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t190c3 + t190c4 * r(i, 3) + t190c5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * t190c5
          d7_b = d12_v * d9_v + d8_b * t190c4
          d14_b = t190c1 * d2_p
          do g_i_ = 1, g_p_
            g_t390(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d
     *7_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          t390 = t190c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t10c2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-t10c2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t10c6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-t10c6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t10c3 + t10c4 * r(i, 3) + t10c5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * t10c5
          d7_b = d12_v * d9_v + d8_b * t10c4
          d14_b = t10c1 * d2_p
          do g_i_ = 1, g_p_
            g_t30(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7
     *_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          t30 = t10c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t10ac2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-t10ac2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t10ac6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-t10ac6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t10ac3 + t10ac4 * r(i, 3) + t10ac5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * t10ac5
          d7_b = d12_v * d9_v + d8_b * t10ac4
          d14_b = t10ac1 * d2_p
          do g_i_ = 1, g_p_
            g_t30alt(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) +
     * d7_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          t30alt = t10ac1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = t10alp * g_r(g_i_, i, 2)
          enddo
          d1_w = t10alp * (r(i, 2) - t10rho)
          d3_v = t30alt - t30
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d7_v = 0.5d0 + 0.5d0 * d5_v
          d8_b = d3_v * 0.5d0 * d1_p
          d2_b = 1.0d0 + (-d7_v)
          do g_i_ = 1, g_p_
            g_t30(g_i_) = d8_b * g_d1_w(g_i_) + d7_v * g_t30alt(g_i_) + 
     *d2_b * g_t30(g_i_)
          enddo
          t30 = t30 + d3_v * d7_v
C--------
C
          d3_v = t30 - t390
          d2_b = 1.0d0 + (-cstsq)
          do g_i_ = 1, g_p_
            g_t3(g_i_) = d3_v * g_cstsq(g_i_) + cstsq * g_t30(g_i_) + d2
     *_b * g_t390(g_i_)
          enddo
          t3 = t390 + d3_v * cstsq
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t1ic2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-t1ic2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t1ic6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-t1ic6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t1ic3 + t1ic4 * r(i, 3) + t1ic5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * t1ic5
          d7_b = d12_v * d9_v + d8_b * t1ic4
          d14_b = t1ic1 * d2_p
          do g_i_ = 1, g_p_
            g_t3i(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7
     *_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          t3i = t1ic1 * d2_v + d10_v * d12_v
C--------
C
          d4_v = 0.5d0 * (1.0d0 - cosa)
          d6_v = t3i - t3
          d2_b = 1.0d0 + (-d4_v)
          d8_b = -(d6_v * 0.5d0)
          do g_i_ = 1, g_p_
            g_t3(g_i_) = d4_v * g_t3i(g_i_) + d8_b * g_cosa(g_i_) + d2_b
     * * g_t3(g_i_)
          enddo
          t3 = t3 + d4_v * d6_v
C--------
C
          d3_v = t30 - t390
          d8_v = 0.5d0 + 0.5d0 * cosa
          d13_v = 0.5d0 - 0.5d0 * cosa
          d7_b = (-t3i) * 0.5d0 + t390 * 0.5d0
          d10_b = d8_v + (-coscos)
          do g_i_ = 1, g_p_
            g_t3(g_i_) = d13_v * g_t3i(g_i_) + d7_b * g_cosa(g_i_) + d3_
     *v * g_coscos(g_i_) + d10_b * g_t390(g_i_) + coscos * g_t30(g_i_)
          enddo
          t3 = d3_v * coscos + t390 * d8_v + t3i * d13_v
C--------
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
C
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
          d6_b = (-cevau) * rac2 * d1_p
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = d6_b * g_d1_w(g_i_) + cevau * g_coul3(g_i_) +
     * cevau * g_coul2(g_i_) + cevau * g_coul1(g_i_)
          enddo
          e(i) = (coul1 + coul2 + coul3 + deh2 - rac2 * d8_v) * cevau
C--------
C
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_elec(g_i_) + g_e(g_i_, i)
          enddo
          e(i) = e(i) + elec
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
C   This file was generated on 07/24/98 by the version of
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
        double precision r, e, amp, erl, ers, eps, rmx, rmod, ralp, frho
     *, thalp, thrho, phalp, phrho, r1s, r2s, r3s, cevau, pi, f1, f3, rm
     *odswt1, rmodswt3, f1p, f3p, f1mod, f3mod, thetaswt, phi, cosa, rna
     *g, rnagsq, chi, falp, sinchi, alpha1, alpha3, dtanh
C

        integer nt, i, j, k
        logical leven
        dimension r(nt, 3), e(nt)
C
        save
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d2_p, d12_b, d7_b, d8_b, d9_v, d10_b, d10_v, d2
     *_v, d9_b, d2_b
        double precision d1_p, d3_v, d4_v, d5_v, d6_v, d3_b, d4_b, d5_b,
     * d6_b, d1_w
        double precision d7_v, d8_v, g_r1s(g_pmax_), g_r(ldg_r, nt, 3), 
     *g_r2s(g_pmax_), g_r3s(g_pmax_), g_alpha1(g_pmax_), g_alpha3(g_pmax
     *_), g_d1_w(g_pmax_), g_f1(g_pmax_)
        double precision g_f3(g_pmax_), g_rnagsq(g_pmax_), g_rnag(g_pmax
     *_), g_chi(g_pmax_), g_cosa(g_pmax_), g_sinchi(g_pmax_), g_rmodswt1
     *(g_pmax_), g_rmodswt3(g_pmax_), g_f1p(g_pmax_), g_f3p(g_pmax_)
        double precision g_f1mod(g_pmax_), g_f3mod(g_pmax_), g_thetaswt(
     *g_pmax_), g_phi(g_pmax_), g_e(ldg_e, nt)
        integer g_ehfid
        save g_phi
        save g_chi, g_cosa, g_sinchi, g_rmodswt1, g_rmodswt3, g_f1p, g_f
     *3p, g_f1mod, g_f3mod, g_thetaswt
        save g_r1s, g_r2s, g_r3s, g_alpha1, g_alpha3, g_d1_w, g_f1, g_f3
     *, g_rnagsq, g_rnag
        data amp /1.5150d0/, erl /0.08877d0/, ers /11.630d0/, eps /1.0d-
     *4/, rmx /2.622d0/, rmod /0.5409d0/, falp /0.7583d0/, frho /5.1858d
     *0/, thalp /2.4882d0/, thrho /2.0d0/, phalp /0.40787d0/, phrho /7.3
     *9370d0/
C
        data g_ehfid /0/
C
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        cevau = 1.d0 / 27.2113961d0
        pi = dacos(-1.0d0)

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
C here we do the function evaluations . . . 
C
          d4_v = eps + r(i, 1) * r1s
          d5_v = ers / d4_v
          d4_b = (-d5_v) / d4_v
          d5_b = d4_b * r1s
          d6_b = d4_b * r(i, 1)
          do g_i_ = 1, g_p_
            g_alpha1(g_i_) = d6_b * g_r1s(g_i_) + d5_b * g_r(g_i_, i, 1)
          enddo
          alpha1 = erl + d5_v
C--------
          d4_v = eps + r(i, 3) * r3s
          d5_v = ers / d4_v
          d4_b = (-d5_v) / d4_v
          d5_b = d4_b * r3s
          d6_b = d4_b * r(i, 3)
          do g_i_ = 1, g_p_
            g_alpha3(g_i_) = d6_b * g_r3s(g_i_) + d5_b * g_r(g_i_, i, 3)
          enddo
          alpha3 = erl + d5_v
C--------
C
          d5_v = (r(i, 1) - rmx) * (r(i, 1) - rmx)
          d1_p = 2.0d0 * (r(i, 1) - rmx)
          d5_b = (-alpha1) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_r(g_i_, i, 1) + (-d5_v) * g_alpha1(g
     *_i_)
          enddo
          d1_w = (-alpha1) * d5_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d3_b = amp * d1_p
          do g_i_ = 1, g_p_
            g_f1(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          f1 = amp * d2_v
C--------
          d5_v = (r(i, 3) - rmx) * (r(i, 3) - rmx)
          d1_p = 2.0d0 * (r(i, 3) - rmx)
          d5_b = (-alpha3) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_r(g_i_, i, 3) + (-d5_v) * g_alpha3(g
     *_i_)
          enddo
          d1_w = (-alpha3) * d5_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d3_b = amp * d1_p
          do g_i_ = 1, g_p_
            g_f3(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          f3 = amp * d2_v
C--------
C
C here follows some checking procedures
C
          do g_i_ = 1, g_p_
            g_rnagsq(g_i_) = (-0.25d0) * g_r2s(g_i_) + 0.5d0 * g_r3s(g_i
     *_) + 0.5d0 * g_r1s(g_i_)
          enddo
          rnagsq = 0.5d0 * r1s + 0.5d0 * r3s - 0.25d0 * r2s
C--------
          if (rnagsq .gt. 0.0d0) then

            rnag = dsqrt(rnagsq)
            if (rnag .gt. 0.0d0) then
              g_rnag(1)=0.5d0*r(i,1)/rnag
              g_rnag(2)=-0.25d0*r(i,2)/rnag
              g_rnag(3)=0.5d0*r(i,3)/rnag
            else
              g_rnag(1)=0.0d0
              g_rnag(2)=0.0d0
              g_rnag(3)=0.0d0
            endif

C--------
            d5_v = 2.0d0 * rnag
            d7_v = d5_v * r(i, 2)
            d8_v = (r3s - r1s) / d7_v
            d9_v = d8_v * d8_v
            d1_p = 2.0d0 * d8_v
            d3_b = d1_p * (1.0d0 / d7_v)
            d4_b = d1_p * ((-d8_v) / d7_v)
            d6_b = d4_b * d5_v
            d7_b = d4_b * r(i, 2) * 2.0d0
            do g_i_ = 1, g_p_
              g_chi(g_i_) = d6_b * g_r(g_i_, i, 2) + d7_b * g_rnag(g_i_)
     * + (-d3_b) * g_r1s(g_i_) + d3_b * g_r3s(g_i_)
            enddo
            chi = d9_v
C--------
          else
            do g_i_ = 1, g_p_
              g_rnag(g_i_) = 0.0d0
            enddo
            rnag = 0.0d0
C--------
            do g_i_ = 1, g_p_
              g_chi(g_i_) = 0.0d0
            enddo
            chi = 0.0d0
C--------
          endif
C
          d7_v = 2.0d0 * r(i, 1)
          d9_v = d7_v * r(i, 3)
          d10_v = (r1s + r3s - r2s) / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = (-d10_v) / d9_v
          d5_b = d3_b * d7_v
          d6_b = d3_b * r(i, 3) * 2.0d0
          do g_i_ = 1, g_p_
            g_cosa(g_i_) = d5_b * g_r(g_i_, i, 3) + d6_b * g_r(g_i_, i, 
     *1) + (-d2_b) * g_r2s(g_i_) + d2_b * g_r3s(g_i_) + d2_b * g_r1s(g_i
     *_)
          enddo
          cosa = d10_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-10.0d0) * g_rnagsq(g_i_)
          enddo
          d1_w = (-10.0d0) * rnagsq
          d2_v = 1.0d0 - chi
          d5_v = 0.5d0 * cosa + 0.5d0
          d6_v = d2_v * d5_v
          d8_v = exp(d1_w)
          d1_p =  d8_v
          d9_v = 1.0d0 - d8_v
          d5_b = (-d6_v) * d1_p
          d9_b = d9_v * d2_v * 0.5d0
          d10_b = -(d9_v * d5_v)
          do g_i_ = 1, g_p_
            g_sinchi(g_i_) = d5_b * g_d1_w(g_i_) + d9_b * g_cosa(g_i_) +
     * d10_b * g_chi(g_i_)
          enddo
          sinchi = d6_v * d9_v
C--------
C
C here ends the procedure
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = falp * g_r(g_i_, i, 3)
          enddo
          d1_w = falp * (r(i, 3) - frho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_v = 0.5d0 - 0.5d0 * d2_v
          d7_v = 0.5d0 + 0.5d0 * cosa
          d5_b = d4_v * 0.5d0
          d8_b = (-d7_v) * 0.5d0 * d1_p
          do g_i_ = 1, g_p_
            g_rmodswt1(g_i_) = d5_b * g_cosa(g_i_) + d8_b * g_d1_w(g_i_)
          enddo
          rmodswt1 = d4_v * d7_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = falp * g_r(g_i_, i, 1)
          enddo
          d1_w = falp * (r(i, 1) - frho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_v = 0.5d0 - 0.5d0 * d2_v
          d7_v = 0.5d0 + 0.5d0 * cosa
          d5_b = d4_v * 0.5d0
          d8_b = (-d7_v) * 0.5d0 * d1_p
          do g_i_ = 1, g_p_
            g_rmodswt3(g_i_) = d5_b * g_cosa(g_i_) + d8_b * g_d1_w(g_i_)
          enddo
          rmodswt3 = d4_v * d7_v
C--------
C
          d5_v = (r(i, 1) - rmx * rmod) * (r(i, 1) - rmx * rmod)
          d1_p = 2.0d0 * (r(i, 1) - rmx * rmod)
          d5_b = (-alpha1) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_r(g_i_, i, 1) + (-d5_v) * g_alpha1(g
     *_i_)
          enddo
          d1_w = (-alpha1) * d5_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d3_b = amp * d1_p
          do g_i_ = 1, g_p_
            g_f1p(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          f1p = amp * d2_v
C--------
C
          d5_v = (r(i, 3) - rmx * rmod) * (r(i, 3) - rmx * rmod)
          d1_p = 2.0d0 * (r(i, 3) - rmx * rmod)
          d5_b = (-alpha3) * d1_p
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d5_b * g_r(g_i_, i, 3) + (-d5_v) * g_alpha3(g
     *_i_)
          enddo
          d1_w = (-alpha3) * d5_v
          d2_v = exp(d1_w)
          d1_p =  d2_v
          d3_b = amp * d1_p
          do g_i_ = 1, g_p_
            g_f3p(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          f3p = amp * d2_v
C--------
C
          d4_v = f1p - f1
          d3_b = 1.0d0 + (-rmodswt1)
          do g_i_ = 1, g_p_
            g_f1mod(g_i_) = d3_b * g_f1(g_i_) + rmodswt1 * g_f1p(g_i_) +
     * d4_v * g_rmodswt1(g_i_)
          enddo
          f1mod = rmodswt1 * d4_v + f1
C--------
C
          d4_v = f3p - f3
          d3_b = 1.0d0 + (-rmodswt3)
          do g_i_ = 1, g_p_
            g_f3mod(g_i_) = d3_b * g_f3(g_i_) + rmodswt3 * g_f3p(g_i_) +
     * d4_v * g_rmodswt3(g_i_)
          enddo
          f3mod = rmodswt3 * d4_v + f3
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = thalp * g_r(g_i_, i, 2)
          enddo
          d1_w = thalp * (r(i, 2) - thrho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_b = (-0.5d0) * d1_p
          do g_i_ = 1, g_p_
            g_thetaswt(g_i_) = d4_b * g_d1_w(g_i_)
          enddo
          thetaswt = 0.5d0 - 0.5d0 * d2_v
C--------
C
          d3_v = f1mod - f3mod
          d6_v = sinchi - 1.0d0
          d8_v = thetaswt * d6_v + 1.0d0
          d5_b = d3_v * d6_v
          d7_b = d3_v * thetaswt
          do g_i_ = 1, g_p_
            g_phi(g_i_) = d7_b * g_sinchi(g_i_) + d5_b * g_thetaswt(g_i_
     *) + (-d8_v) * g_f3mod(g_i_) + d8_v * g_f1mod(g_i_)
          enddo
          phi = d3_v * d8_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = phalp * g_r(g_i_, i, 2)
          enddo
          d1_w = phalp * (r(i, 2) - phrho)
          d2_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d2_v *  d2_v)
          d4_v = 0.5d0 * (1.d0 - d2_v)
          d7_v = 0.5d0 - 0.5d0 * cosa
          d9_v = 1.0d0 - d4_v * d7_v
          d4_b = cevau * d9_v
          d5_b = -(cevau * phi)
          d9_b = (-(d5_b * d4_v)) * 0.5d0
          d12_b = (-(d5_b * d7_v * 0.5d0)) * d1_p
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = d4_b * g_phi(g_i_) + d9_b * g_cosa(g_i_) + d1
     *2_b * g_d1_w(g_i_)
          enddo
          e(i) = d9_v * phi * cevau
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
C   This file was generated on 07/24/98 by the version of
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
C the upper diabatic surface U22
C
      subroutine g_prepot22(g_p_)
        implicit double precision (a-h, o-z)

        dimension r(nt, 3), e(nt)
        integer i,j,k,nsurf,nt

        save
C
C the NaH singlet data blocks
        integer g_pmax_
        parameter (g_pmax_ = 3)
        integer g_i_, g_p_, ldg_r, ldg_e
        double precision d3_p, d24_b, d2_p, d15_v, d20_v, d16_v, d22_v, 
     *d2_v, d16_b, d2_b
        double precision d1_w, d2_w, d3_v, d4_v, d5_v, d6_v, d3_b, d4_b,
     * d5_b, d6_b
        double precision d1_p, d7_v, d8_v, d9_v, d7_b, d8_b, d9_b, d10_v
     *, d10_b, d11_v
        double precision d12_v, d13_v, d19_v, d11_b, d12_b, d13_b, d14_b
     *, g_r1s(g_pmax_), g_r(ldg_r, nt, 3), g_r2s(g_pmax_)
        double precision g_r3s(g_pmax_), g_d2_w(g_pmax_), g_d1_w(g_pmax_
     *), g_rnag(g_pmax_), g_rnags(g_pmax_), g_rnagc(g_pmax_), g_cstsq(g_
     *pmax_), g_cosa(g_pmax_), g_coscos(g_pmax_), g_rnag1(g_pmax_)
        double precision g_rnag3(g_pmax_), g_rnag1c(g_pmax_), g_rnag3c(g
     *_pmax_), g_cstsq1(g_pmax_), g_cstsq3(g_pmax_), g_s1a(g_pmax_), g_s
     *1ic(g_pmax_), g_s1i(g_pmax_), g_s1(g_pmax_), g_t10(g_pmax_)
        double precision g_t190(g_pmax_), g_t1(g_pmax_), g_t2a(g_pmax_),
     * g_t290c(g_pmax_), g_t290(g_pmax_), g_t20c(g_pmax_), g_t20(g_pmax_
     *), g_t2(g_pmax_), g_s2a(g_pmax_), g_s2mod(g_pmax_)
        double precision g_s290c(g_pmax_), g_s290(g_pmax_), g_s2(g_pmax_
     *), g_s3a(g_pmax_), g_s3ic(g_pmax_), g_s3i(g_pmax_), g_s3(g_pmax_),
     * g_t30(g_pmax_), g_t390(g_pmax_), g_t3(g_pmax_)
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
        save g_c6swt3, g_coul1, g_coul2, g_coul3, g_exch1, g_exch2, g_ex
     *ch3, g_w, g_cplg2
        save g_rmuimu3, g_c6, g_c61, g_dexp1, g_c63, g_dexp3, g_cind1, g
     *_cind3, g_elec, g_c6swt1
        save g_rnahpol3, g_rnahpol3b, g_qqint, g_qq, g_dexp2, g_qhint, g
     *_qh, g_div1, g_rmuimu1, g_div3
        save g_rnahie1, g_rnahie3, g_rnahmu1, g_rnahmu3, g_rnahpar1, g_r
     *nahpar3, g_rnahperp1, g_rnahperp3, g_rnahpol1, g_rnahpol1b
        save g_t390, g_t3, g_h2ie, g_al1, g_al2, g_alt, g_alb, g_alhh, g
     *_h2hex, g_h2quad
        save g_s2a, g_s2mod, g_s290c, g_s290, g_s2, g_s3a, g_s3ic, g_s3i
     *, g_s3, g_t30
        save g_s1, g_t10, g_t190, g_t1, g_t2a, g_t290c, g_t290, g_t20c, 
     *g_t20, g_t2
        save g_coscos, g_rnag1, g_rnag3, g_rnag1c, g_rnag3c, g_cstsq1, g
     *_cstsq3, g_s1a, g_s1ic, g_s1i
        save g_r1s, g_r2s, g_r3s, g_d2_w, g_d1_w, g_rnag, g_rnags, g_rna
     *gc, g_cstsq, g_cosa
        intrinsic dble
        data s1c1 /170.769d0/, s1c2 /1.514d0/, s1c3 /-87.110d0/, s1c4 /5
     *7.480d0/, s1c5 /-14.228d0/, s1c6 /1.303d0/
        data s1ic1 /460.09449d0/, s1ic2 /2.45039d0/, s1ic3 /-29.81890d0/
     *, s1ic4 /42.67717d0/, s1ic5 /-13.11811d0/, s1ic6 /1.303d0/, s1ialp
     * /0.28780d0/
C
C the NaH triplet data blocks
        data t10c1 /385.59843d0/, t10c2 /2.83307d0/, t10c3 /67.06299d0/,
     * t10c4 /-37.70079d0/, t10c5 /6.59055d0/, t10c6 /1.303d0/
        data t190c1 /366.53543d0/, t190c2 /1.74567d0/, t190c3 /-50.0d0/,
     * t190c4 /4.20472d0/, t190c5 /1.22047d0/, t190c6 /1.303d0/
C
C The H2 triplet
        data t2c1 /147.214d0/, t2c2 /3.915d0/, t2c3 /41.275d0/, t2c4 /10
     *.505d0/, t2c5 /4.408d0/, t2c6 /2.032d0/
        data t290c1 /82.88976d0/, t290c2 /3.38898d0/, t290c3 /10.87402d0
     */, t290c4 /-7.88189d0/, t290c5 /-2.35433d0/, t290c6 /2.032d0/, t29
     *0alp /0.09543d0/, t290rho /40.81102d0/
        data t20c1 /159.34646d0/, t20c2 /3.51260d0/, t20c3 /121.65354d0/
     *, t20c4 /9.41732d0/, t20c5 /-11.03937d0/, t20c6 /2.032/, t20alp /0
     *.11260d0/, t20rho /43.22799504d0/
C
C The H2 singlet
C
        data s2c1 /139.7160d0/, s2c2 /-123.8978d0/, s2c3 /3.4031d0/, s2c
     *4 /-6.8725d0/, s2c5 /-23.0440d0/, s2c6 /2.032d0/
        data s290c1 /133.63780d0/, s290c2 /422.23622d0/, s290c3 /5.56890
     *d0/, s290c4 /-4.22441d0/, s290c5 /-9.79213d0/, s290c6 /2.032d0/  s
     *290alp /0.05906d0/, s290rho /12.72449d0/
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
C misc constants
        data c2a /2.18504d0/, c2b /0.10945d0/, c2c /0.15039d0/, deh2 /4.
     *74772265651191638d0/, una2p /2.1037/
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
          rnag = dsqrt(dabs(0.5d0*(r1s+r3s-0.5d0*r2s)))

          if (rnag .gt. 0.0d0) then
            g_rnag(1)=0.5d0*r(i,1)/rnag
            g_rnag(2)=-0.25d0*r(i,2)/rnag
            g_rnag(3)=0.5d0*r(i,3)/rnag
          else
            g_rnag(1)=0.0d0
            g_rnag(2)=0.0d0
            g_rnag(3)=0.0d0
          endif

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
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-10.0d0) * g_rnags(g_i_)
          enddo
          d1_w = (-10.0d0) * rnags
          d4_v = 0.5d0 + 0.5d0 * cosa
          d5_v = cstsq * d4_v
          d7_v = exp(d1_w)
          d1_p =  d7_v
          d8_v = 1.0d0 - d7_v
          d5_b = (-d5_v) * d1_p
          d6_b = d8_v * d4_v
          d9_b = d8_v * cstsq * 0.5d0
          do g_i_ = 1, g_p_
            g_coscos(g_i_) = d5_b * g_d1_w(g_i_) + d9_b * g_cosa(g_i_) +
     * d6_b * g_cstsq(g_i_)
          enddo
          coscos = d5_v * d8_v
C--------
C
          frac = 0.958d0

           rnag1= dsqrt(dabs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))

           if (rnag1 .gt. 0.0d0) then
             g_rnag1(1) = r(i,1)*frac*(frac-1.0d0)/rnag1
             g_rnag1(2) = r(i,2)*(1.0d0-frac)/rnag1
             g_rnag1(3) = r(i,3)*frac/rnag1
           else
             g_rnag1(1) = 0.0d0
             g_rnag1(2) = 0.0d0
             g_rnag1(3) = 0.0d0
           endif

C--------

           if (rnag3 .gt. 0.0d0) then
             g_rnag3(1) = r(i,1)*frac/rnag3
             g_rnag3(2) = r(i,2)*(1.0d0-frac)/rnag3
             g_rnag3(3) = r(i,3)*frac*(frac-1.0d0)/rnag3
           else
             g_rnag3(1) = 0.0d0
             g_rnag3(2) = 0.0d0
             g_rnag3(3) = 0.0d0
           endif

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
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s1c2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-s1c2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s1c6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-s1c6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s1c3 + s1c4 * r(i, 1) + s1c5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * s1c5
          d7_b = d12_v * d9_v + d8_b * s1c4
          d14_b = s1c1 * d2_p
          do g_i_ = 1, g_p_
            g_s1a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7
     *_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          s1a = s1c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s1ic2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-s1ic2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s1ic6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-s1ic6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s1ic3 + s1ic4 * r(i, 1) + s1ic5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * s1ic5
          d7_b = d12_v * d9_v + d8_b * s1ic4
          d14_b = s1ic1 * d2_p
          do g_i_ = 1, g_p_
            g_s1ic(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d
     *7_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          s1ic = s1ic1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s1ialp * g_rnags(g_i_)
          enddo
          d1_w = s1ialp * rnags
          d3_v = s1ic - s1a
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d6_v = 1.0d0 - d5_v
          d7_b = (-d3_v) * d1_p
          d2_b = 1.0d0 + (-d6_v)
          do g_i_ = 1, g_p_
            g_s1i(g_i_) = d7_b * g_d1_w(g_i_) + d6_v * g_s1ic(g_i_) + d2
     *_b * g_s1a(g_i_)
          enddo
          s1i = s1a + d3_v * d6_v
C--------
C
          d4_v = 0.5d0 * (1.0d0 - cosa)
          d6_v = s1i - s1a
          d2_b = 1.0d0 + (-d4_v)
          d8_b = -(d6_v * 0.5d0)
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d4_v * g_s1i(g_i_) + d8_b * g_cosa(g_i_) + d2_b
     * * g_s1a(g_i_)
          enddo
          s1 = s1a + d4_v * d6_v
C--------
C
C
C The first NaH triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t10c2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-t10c2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t10c6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-t10c6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t10c3 + t10c4 * r(i, 1) + t10c5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * t10c5
          d7_b = d12_v * d9_v + d8_b * t10c4
          d14_b = t10c1 * d2_p
          do g_i_ = 1, g_p_
            g_t10(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d7
     *_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          t10 = t10c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t190c2) * g_r(g_i_, i, 1)
          enddo
          d1_w = (-t190c2) * r(i, 1)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t190c6) * g_r(g_i_, i, 1)
          enddo
          d2_w = (-t190c6) * r(i, 1)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t190c3 + t190c4 * r(i, 1) + t190c5 * r1s
          d10_v = r(i, 1) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 1)
          d11_b = d8_b * t190c5
          d7_b = d12_v * d9_v + d8_b * t190c4
          d14_b = t190c1 * d2_p
          do g_i_ = 1, g_p_
            g_t190(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r1s(g_i_) + d
     *7_b * g_r(g_i_, i, 1) + d14_b * g_d1_w(g_i_)
          enddo
          t190 = t190c1 * d2_v + d10_v * d12_v
C--------
C
          d4_v = 0.5d0 - 0.5d0 * cosa
          d6_v = t10 - t190
          d9_v = t10 - t190
          d6_b = coscos + d4_v
          d7_b = -coscos + 1.0d0 + (-d4_v)
          d12_b = (-d6_v) * 0.5d0
          do g_i_ = 1, g_p_
            g_t1(g_i_) = d9_v * g_coscos(g_i_) + d6_b * g_t10(g_i_) + d1
     *2_b * g_cosa(g_i_) + d7_b * g_t190(g_i_)
          enddo
          t1 = t190 + d4_v * d6_v + d9_v * coscos
C--------
C
C The H2 triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t2c2) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-t2c2) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t2c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-t2c6) * r(i, 2)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t2c3 + t2c4 * r(i, 2) + t2c5 * r2s
          d10_v = r(i, 2) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 2)
          d11_b = d8_b * t2c5
          d7_b = d12_v * d9_v + d8_b * t2c4
          d14_b = t2c1 * d2_p
          do g_i_ = 1, g_p_
            g_t2a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d7
     *_b * g_r(g_i_, i, 2) + d14_b * g_d1_w(g_i_)
          enddo
          t2a = t2c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t290c2) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-t290c2) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t290c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-t290c6) * r(i, 2)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t290c3 + r(i, 2) * t290c4 + t290c5 * r2s
          d10_v = r(i, 2) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 2)
          d11_b = d8_b * t290c5
          d7_b = d12_v * d9_v + d8_b * t290c4
          d14_b = t290c1 * d2_p
          do g_i_ = 1, g_p_
            g_t290c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + 
     *d7_b * g_r(g_i_, i, 2) + d14_b * g_d1_w(g_i_)
          enddo
          t290c = t290c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = t290alp * g_rnags(g_i_)
          enddo
          d1_w = t290alp * (rnags - t290rho)
          d3_v = t290c - t2a
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d6_v = 1.0d0 - d5_v
          d5_b = 0.5d0 * d6_v
          d8_b = (-(0.5d0 * d3_v)) * d1_p
          d2_b = 1.0d0 + (-d5_b)
          do g_i_ = 1, g_p_
            g_t290(g_i_) = d8_b * g_d1_w(g_i_) + d5_b * g_t290c(g_i_) + 
     *d2_b * g_t2a(g_i_)
          enddo
          t290 = t2a + d3_v * d6_v * 0.5d0
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t20c2) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-t20c2) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t20c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-t20c6) * r(i, 2)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t20c3 + t20c4 * r(i, 2) + t20c5 * r2s
          d10_v = r(i, 2) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 2)
          d11_b = d8_b * t20c5
          d7_b = d12_v * d9_v + d8_b * t20c4
          d14_b = t20c1 * d2_p
          do g_i_ = 1, g_p_
            g_t20c(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r2s(g_i_) + d
     *7_b * g_r(g_i_, i, 2) + d14_b * g_d1_w(g_i_)
          enddo
          t20c = t20c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = t20alp * g_rnags(g_i_)
          enddo
          d1_w = t20alp * (rnags - t20rho)
          d3_v = t20c - t2a
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d6_v = 1.0d0 - d5_v
          d5_b = 0.5d0 * d6_v
          d8_b = (-(0.5d0 * d3_v)) * d1_p
          d2_b = 1.0d0 + (-d5_b)
          do g_i_ = 1, g_p_
            g_t20(g_i_) = d8_b * g_d1_w(g_i_) + d5_b * g_t20c(g_i_) + d2
     *_b * g_t2a(g_i_)
          enddo
          t20 = t2a + d3_v * d6_v * 0.5d0
C--------
C
          d3_v = t20 - t290
          d2_b = 1.0d0 + (-coscos)
          do g_i_ = 1, g_p_
            g_t2(g_i_) = d3_v * g_coscos(g_i_) + coscos * g_t20(g_i_) + 
     *d2_b * g_t290(g_i_)
          enddo
          t2 = t290 + d3_v * coscos
C--------
C
C The H2 singlet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s2c3) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-s2c3) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s2c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-s2c6) * r(i, 2)
          d3_v = s2c1 + s2c2 * r(i, 2)
          d5_v = exp(d1_w)
          d2_p =  d5_v
          d9_v = s2c4 + s2c5 * r(i, 2)
          d10_v = r2s * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d7_b = d10_v * d1_p
          d8_b = d12_v * d9_v
          d14_b = d3_v * d2_p
          d11_b = d12_v * r2s * s2c5 + d5_v * s2c2
          do g_i_ = 1, g_p_
            g_s2a(g_i_) = d7_b * g_d2_w(g_i_) + d8_b * g_r2s(g_i_) + d14
     *_b * g_d1_w(g_i_) + d11_b * g_r(g_i_, i, 2)
          enddo
          s2a = d3_v * d5_v + d10_v * d12_v + una2p
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
            g_s2mod(g_i_) = d9_b * g_d1_w(g_i_) + d8_b * g_t2(g_i_) + d7
     *_b * g_s2a(g_i_)
          enddo
          s2mod = 0.5d0 * (s2a + t2 - d5_v * d6_v)
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s290c3) * g_r(g_i_, i, 2)
          enddo
          d1_w = (-s290c3) * r(i, 2)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s290c6) * g_r(g_i_, i, 2)
          enddo
          d2_w = (-s290c6) * r(i, 2)
          d3_v = s290c1 + s290c2 * r(i, 2)
          d5_v = exp(d1_w)
          d2_p =  d5_v
          d9_v = s290c4 + s290c5 * r(i, 2)
          d10_v = r2s * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d7_b = d12_v * d9_v
          d13_b = d3_v * d2_p
          d10_b = d12_v * r2s * s290c5 + d5_v * s290c2
          do g_i_ = 1, g_p_
            g_s290c(g_i_) = d6_b * g_d2_w(g_i_) + d7_b * g_r2s(g_i_) + d
     *13_b * g_d1_w(g_i_) + d10_b * g_r(g_i_, i, 2)
          enddo
          s290c = d3_v * d5_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s290alp * g_rnags(g_i_)
          enddo
          d1_w = s290alp * (rnags - s290rho)
          d4_v = (s290c - s2mod) * 0.5d0
          d6_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d6_v *  d6_v)
          d7_v = 1.0d0 - d6_v
          d7_b = (-d4_v) * d1_p
          d8_b = d7_v * 0.5d0
          d2_b = 1.0d0 + (-d8_b)
          do g_i_ = 1, g_p_
            g_s290(g_i_) = d7_b * g_d1_w(g_i_) + d8_b * g_s290c(g_i_) + 
     *d2_b * g_s2mod(g_i_)
          enddo
          s290 = s2mod + d4_v * d7_v
C--------
C
          d3_v = s2mod - s290
          d2_b = 1.0d0 + (-coscos)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d3_v * g_coscos(g_i_) + coscos * g_s2mod(g_i_) 
     *+ d2_b * g_s290(g_i_)
          enddo
          s2 = s290 + d3_v * coscos
C--------
C
C The second NaH singlet
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s1c2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-s1c2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s1c6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-s1c6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s1c3 + s1c4 * r(i, 3) + s1c5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * s1c5
          d7_b = d12_v * d9_v + d8_b * s1c4
          d14_b = s1c1 * d2_p
          do g_i_ = 1, g_p_
            g_s3a(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7
     *_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          s3a = s1c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-s1ic2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-s1ic2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-s1ic6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-s1ic6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = s1ic3 + s1ic4 * r(i, 3) + s1ic5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * s1ic5
          d7_b = d12_v * d9_v + d8_b * s1ic4
          d14_b = s1ic1 * d2_p
          do g_i_ = 1, g_p_
            g_s3ic(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d
     *7_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          s3ic = s1ic1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = s1ialp * g_rnags(g_i_)
          enddo
          d1_w = s1ialp * rnags
          d3_v = s3ic - s3a
          d5_v = tanh (d1_w)
          d1_p = 1.0d0 - ( d5_v *  d5_v)
          d6_v = 1.0d0 - d5_v
          d7_b = (-d3_v) * d1_p
          d2_b = 1.0d0 + (-d6_v)
          do g_i_ = 1, g_p_
            g_s3i(g_i_) = d7_b * g_d1_w(g_i_) + d6_v * g_s3ic(g_i_) + d2
     *_b * g_s3a(g_i_)
          enddo
          s3i = s3a + d3_v * d6_v
C--------
C
          d4_v = 0.5d0 * (1.0d0 - cosa)
          d6_v = s3i - s3a
          d2_b = 1.0d0 + (-d4_v)
          d8_b = -(d6_v * 0.5d0)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d4_v * g_s3i(g_i_) + d8_b * g_cosa(g_i_) + d2_b
     * * g_s3a(g_i_)
          enddo
          s3 = s3a + d4_v * d6_v
C--------
C
C
C The second NaH triplet
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t10c2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-t10c2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t10c6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-t10c6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t10c3 + t10c4 * r(i, 3) + t10c5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * t10c5
          d7_b = d12_v * d9_v + d8_b * t10c4
          d14_b = t10c1 * d2_p
          do g_i_ = 1, g_p_
            g_t30(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d7
     *_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          t30 = t10c1 * d2_v + d10_v * d12_v
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = (-t190c2) * g_r(g_i_, i, 3)
          enddo
          d1_w = (-t190c2) * r(i, 3)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = (-t190c6) * g_r(g_i_, i, 3)
          enddo
          d2_w = (-t190c6) * r(i, 3)
          d2_v = exp(d1_w)
          d2_p =  d2_v
          d9_v = t190c3 + t190c4 * r(i, 3) + t190c5 * r3s
          d10_v = r(i, 3) * d9_v
          d12_v = exp(d2_w)
          d1_p =  d12_v
          d6_b = d10_v * d1_p
          d8_b = d12_v * r(i, 3)
          d11_b = d8_b * t190c5
          d7_b = d12_v * d9_v + d8_b * t190c4
          d14_b = t190c1 * d2_p
          do g_i_ = 1, g_p_
            g_t390(g_i_) = d6_b * g_d2_w(g_i_) + d11_b * g_r3s(g_i_) + d
     *7_b * g_r(g_i_, i, 3) + d14_b * g_d1_w(g_i_)
          enddo
          t390 = t190c1 * d2_v + d10_v * d12_v
C--------
C
          d4_v = 0.5d0 - 0.5d0 * cosa
          d6_v = t30 - t390
          d9_v = t30 - t390
          d6_b = coscos + d4_v
          d7_b = -coscos + 1.0d0 + (-d4_v)
          d12_b = (-d6_v) * 0.5d0
          do g_i_ = 1, g_p_
            g_t3(g_i_) = d9_v * g_coscos(g_i_) + d6_b * g_t30(g_i_) + d1
     *2_b * g_cosa(g_i_) + d7_b * g_t390(g_i_)
          enddo
          t3 = t390 + d4_v * d6_v + d9_v * coscos
C--------
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
          do g_i_ = 1, g_p_
            g_e(g_i_, i) = g_elec(g_i_) + cevau * g_e(g_i_, i)
          enddo
          e(i) = e(i) * cevau + elec
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
