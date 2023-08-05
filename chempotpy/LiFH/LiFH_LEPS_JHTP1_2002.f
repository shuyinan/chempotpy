      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: nt, r(1,3), v(1)
      double precision :: u11(1), u12(1), u22(1)
      double precision :: u(nstates,nstates)
      logical, save :: first_time_data=.true.
      integer :: i, j, k, l

      !initialize 
      u=0.d0
      p=0.d0
      g=0.d0
      d=0.d0

      nt=1
      ! input cartesian is ClHH
      r(1,1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      r(1,2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      r(1,3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211

      if(first_time_data) then
      call prepot
      first_time_data=.false.
      endif

      if (igrad==0) then 
        call pot(r,u11,1,1)
        call pot(r,u12,1,2)
        call pot(r,u22,1,3)
        u(1,1)=u11(1)*27.211386
        u(1,2)=u12(1)*27.211386
        u(2,1)=u12(1)*27.211386
        u(2,2)=u22(1)*27.211386
        call diagonalize(nstates,u,p,t)
      else
        write(*,*) "Only energy is available"
      endif

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
C   Common name:     LiFH surface fit JS
C   Number of derivatives: 0
C   Number of electronic surfaces: 2
C   Interface: 3-1V

C   References:      A. W. Jasper, M. D. Hack, D. G. Truhlar, and P. Piecuch, 
C                    J. Chem. Phys., Vol. 116, 8353 (2002).

c   Notes:    This fit contains cut-off long-range forces.
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
c                   LiFHJS POTENTIAL ENERGY MATRIX
c
c  Created by A. Jasper on August 22, 2001.  
c  Based on the LiFHJ potential energy matrix, which was
c  based on the LiFHH potential energy matrix created by M. Hack,  
c  which was based on the ab initio calculations of P. Piecuch.
c
c  This surface has the long-range forces cut off.
c
c  Reference: A. W. Jasper, M. D. Hack, D. G. Truhlar, and P. Piecuch, 
c  J. Chem. Phys., in preparation.
c
c ********************************************************************
c ********************************************************************

      subroutine prepot
        integer g_p_
        parameter(g_p_=3)
        common /com_para/ g_myrank,g_nprocs
        integer g_myrank, g_nprocs

        call prepot11
        call prepot12
        call prepot22

      return

      end

c ********************************************************************
c ********************************************************************

      subroutine pot(r,e,nt,nsurf)

        implicit none
        double precision r, e
        integer i,nt,nsurf
        dimension r(nt,3),e(nt)

        if (nsurf.eq.1) call pot11(r,e,nt)
        if (nsurf.eq.2) call pot12(r,e,nt)
        if (nsurf.eq.3) call pot22(r,e,nt)
        if (nsurf.gt.3.or.nsurf.lt.1) then
           write(6,*)' WARNING:  nsurf = ',nsurf,
     &           ' in LiFHJS potential routine!'
        endif

      return
      end

c ********************************************************************
c ********************************************************************

c                            U11 SURFACE

      subroutine prepot11
      implicit none

c MISC VARIABLES
      integer i,nt,j
      double precision cevau
      double precision r(nt,3),e(nt)
      double precision r1s,r2s,r3s

c INTERACTION REGION VARIABLES 
c     LiH
      double precision s1,relih,bes1,s1c2,s1c3,s1a,xelih,bes1b,s1c,
     &  rcut2del,rcut2rho,rcut2,r1r3,r1r3al

c     HF
      double precision s2,behf,yhf,rehf,dehf,xehf,s2m,behfm,dehfm,
     & rehfm,xehfm,s2rho,s2alp,dehfc,dehfcc,yhf_r,behf_0,behf_1,
     & behf_2

c     LiF
      double precision s3,res3,bes3,des3,res3c,bes3c,s3c,s3a,xes3,des3c,
     & rcutdel,rcutrho,rcut,r1r3al2,bes3i,bes3z,bes3r

c LONG RANGE FORCES VARIABLES
c     HF physical properties
      double precision hfdip_y,hfdip_a,hfdip_re,hfdip_m(11),hfdip,
     & hfpol,hfpolperp,hfpolpar,hfie,hfpolavg

c     Li physical properties
      double precision li2spol,rhf0did,rhf0disp,li2sie

c     misc long-range variables
      double precision hf_did,hf_disp,e_LR_HF
      double precision mass(3),cauang,caudeb,costheta,temp,bigR,
     & coschi,r2LRcut_a,r2LRcut_r

c MISC DATA
      data cevau /0.036749309d0/

c INTERACTION REGION DATA
c     LiH
      data relih /1.3432062d0/, bes1 /1.8058651d0/, s1c2 /12.937438d0/
      data s1c3 /11.0913d0/, bes1b /1.2111436d0/
      data rcut2rho /1.07143d0/, rcut2del /0.6d0/, r1r3al /0.4d0/

c     HF
      data rehf /1.733d0/, dehf /6.122d0/, rehfm/ 1.6739d0/,
     & dehfc/ 0.344574d0/,dehfcc/0.8660801d0/,s2rho/1.63d0/,
     & s2alp/2.9941d0/,
     & behfm/0.74633431d0/,yhf_r/2.1042d0/,behf_0/1.1622d0/,
     & behf_1/0.025647d0/,behf_2/0.059062d0/

c     LiF
      data res3/3.6d0/,bes3/1.49061584d0/,des3/0.17426687d0/
      data res3c /2.48571d0/, bes3c /1.92857d0/, des3c /0.25d0/
      data rcutrho /0.4d0/, rcutdel /0.89333d0/, r1r3al2 /0.7333d0/,
     & bes3z/0.d0/,bes3r/7.4799608d0/

c LONG RANGE FORCES DATA
c     HF dipole moment (angstroms vs au)
      data hfdip_a /1.13d0/, hfdip_re /0.924212d0/
      data hfdip_m /0.703661d0,0.516815d0,0.240628d0,-0.430194d0,
     & -2.213088d0,5.535609d0,12.872106d0,-42.060172d0,-17.398065d0,
     & 84.207629d0,-41.961465d0/

c     HF physical properties (polarization in au, ie in eV)
      data hfpolperp/4.59d0/,hfpolpar/5.10d0/,hfie/16.044d0/

c     Li(2s) physical properties (polarization in au, ie in eV)
      data li2spol /165.d0/,li2sie/5.392d0/

c     misc long-range parameters (bohr)
      data rhf0did/6.0d0/,rhf0disp/6.0d0/,r2LRcut_a/2.0d0/,
     & r2LRcut_r/3.d0/
      data cauang /0.52917706d0/,caudeb/2.54177d0/
      data mass /7.016003d0,1.00783d0,18.9984d0/

      save

      return

c ---------------------------------------------------------------------- 

      entry pot11(r,e,nt)

      do i=1,nt

c PRECOMPUTE THINGS FOR USE LATER...
c     cosine of HF bond angle (i.e., Li-F-H angle)
      costheta = (r(i,3)**2-r(i,2)**2-r(i,1)**2)/(-2.d0*r(i,1)*r(i,2))
      temp = mass(3)/(mass(2)+mass(3))*r(i,2)
c     Li-FH Jacobi distance
      bigR = dsqrt(r(i,1)**2+temp**2-2.d0*r(i,1)*temp*costheta)
c     cosine of HF Jacobi angle (i.e., Li--HF angle)
      coschi = -(r(i,1)**2-bigR**2-temp**2)/(-2.d0*bigR*temp)
      coschi = max(coschi,-1.d0)
      coschi = min(coschi,1.d0)
c     squared!
      r1s = r(i,1)**2
      r2s = r(i,2)**2
      r3s = r(i,3)**2

c INTERACTION REGION
c     LiH
      xelih = dexp(-bes1*(r(i,1)-relih))
      s1a =  s1c2*xelih 
      s1c = s1c3*dexp(-bes1b*(r(i,1)-relih))

      r1r3 = r(i,1) - r(i,3) + r1r3al*r(i,2)
      rcut2 = 0.5d0*(1.0d0 - dtanh((r1r3-rcut2rho)/rcut2del))

      s1 =  s1a  + (s1c-s1a)*rcut2

c     HF
      yhf = r(i,2) - yhf_r
      behf = behf_0-behf_1*yhf+behf_2*yhf**2
      xehf = dexp(-behf*(r(i,2)-rehf))
      s2 = dehf*xehf*(xehf-2.0d0)

      dehfm = dehf - (dehfcc+dehfc*0.5d0*(1.d0-coschi))
      xehfm = dexp(-behfm*(r(i,2)-rehfm))
      s2m = dehfm*xehfm*(xehfm-2.0d0)
      s2 = s2m + (s2 - s2m)*0.5d0*(1.d0+dtanh((r(i,3)-s2rho)/s2alp))

c     LiF 
      xes3 = dexp(-bes3*(r(i,3)-res3))
      s3 = des3*xes3*(xes3-2.0d0)

c     Sum of diatomics
      e(i) = s1 + s2 + s3 + dehf  

c EVALUATE THE LONG-RANGE FORCES
c HF DIPOLE MOMENT
      hfdip_y = 1.d0 - dexp(-hfdip_a*(r(i,2)*cauang-hfdip_re))
      hfdip = 0.d0
      do j=1,11
      hfdip = hfdip + hfdip_m(j)*hfdip_y**dble(j-1)
      enddo
      hfdip = hfdip*0.5d0*(1.d0 + dtanh(20.d0*(r(i,2)-1.2d0)))

c HF POLARIZATION
      hfpol = hfpolperp*(1.d0 - coschi**2) + hfpolpar*coschi**2
      hfpolavg = 2.d0*hfpolperp/3.d0 + hfpolpar/3.d0

c HF + Li(2s) arrangement
c     DIPOLE-INDUCED-DIPOLE
      hf_did = -0.5d0*(3.d0*coschi**2 + 1.d0)*li2spol*hfdip**2
      hf_did = hf_did/(bigR**7 + rhf0did**7)/cevau*bigR
c     cut off angular dependence in interaction region
      hf_did = hf_did*(1.d0 + (2.5d0/(3.d0*coschi**2 + 1.d0) - 1.d0)
     &               *dexp(-(bigR/rhf0did)**4))

c     DISPERSION
      hf_disp = -1.5d0*li2sie*hfie*li2spol*hfpol/(li2sie+hfie)
      hf_disp = hf_disp/(bigR**7 + rhf0disp**7)*bigR
c     cut off angular dependence in interaction region
      hf_disp = hf_disp*(1.d0 + (hfpolavg/hfpol - 1.d0)
     &                 *dexp(-(bigR/rhf0disp)**4))

c     SUM FORCES
      e_LR_HF = hf_disp + hf_did
c     cut off LR forces for large r[hf] i.e., r(i,2)
      e_LR_HF=e_LR_HF*0.5d0*(1.d0-dtanh(r2LRcut_a*(r(i,2)-r2LRcut_r)))

c LiF + H
c     NO LONG-RANGE FORCES

c LiH + F
c     NO LONG-RANGE FORCES

c CUT OFF LONG-RANGE FORCES
      if (bigR.gt.15.d0) then
        e_LR_HF = 0.d0
      else
        e_LR_HF = e_LR_HF * dexp( .2d0 / (bigR - 15.d0) )
      endif

c COMBINE INTERACTION REGION ENERGY AND LONG-RANGE FORCES
      e(i) = e(i) + e_LR_HF
      e(i) = e(i)*cevau

      enddo
      return

      end

c ********************************************************************
c ********************************************************************

c                      U12 COUPLING SURFACE

      subroutine prepot12

      implicit none 

c MISC VARIABLES
      integer i,nt
      double precision r(nt,3),e(nt)

c LiH VARIABLES
      integer ng1
      double precision g1,g1r,g1a,h1,h1del,h1rho

c HF VARIABLES
      double precision j2p,j2,j2pdel,j2rho,j2prho,j2del

c LiF VARIABLES
      integer ng3
      double precision g3,g3r,g3a,h3,h3del,h3rho

c CUTOFF VARIABLES
      double precision gcrho,cevau,gcdel
      double precision cut2,gc,cut3

c MISC DATA
      data cevau /0.036749309d0/

c LiH DATA
      data g1a /1.27742d0/, g1r /3.0048876d0/, ng1 /6/
      data h1rho /4.9843597d0/, h1del /2.0899315d0/

c LiF DATA
      data g3a /0.48d0/, g3r /3.4975562d0/, ng3 /8/
      data h3rho /2.000978d0/, h3del /0.90615835d0/

c HF DATA
      data  j2rho /1.15484d0/, j2del /1.75806d0/
      data  j2prho /1.45161d0/, j2pdel /0.98387d0/

c CUTOFF DATA
      data  gcrho /4.87097d0/, gcdel /2.d0/

      save
      return

c ---------------------------------------------------------------------- 

      entry pot12(r,e,nt)

      do i=1,nt
  
        g1  = g1a*(r(i,1)/g1r)**ng1
     &             *dexp(-dble(ng1)*(r(i,1)/g1r-1.0d0))

        g3  = g3a*(r(i,3)/g3r)**ng3
     &             *dexp(-dble(ng3)*(r(i,3)/g3r-1.0d0))

        j2  = 0.5d0*(1.0d0 + dtanh( (r(i,2) - j2rho)/j2del ) )
        j2p = 0.5d0*(1.0d0 + dtanh( (r(i,2) - j2prho)/j2pdel ) )

        h1  = 0.5d0*(1.0d0 - dtanh( (r(i,1) - h1rho)/h1del ) )
        h3  = 0.5d0*(1.0d0 - dtanh( (r(i,3) - h3rho)/h3del ) )

        e(i) = g1*j2*(1.0d0-h3) + g3*j2p*(1.0d0-h1)

        gc = r(i,2)-gcrho
        cut2 = 0.5d0 - 0.5d0*dtanh(gc/gcdel)

        cut3 = 0.5d0 - 0.5d0*dtanh((r(i,2) - 5.5d0)/.2d0)

        e(i) = e(i)*cut2*cut3

        e(i) = e(i)*cevau

      enddo
      return
      end

c ********************************************************************
c ********************************************************************

c                            U22 SURFACE

      subroutine prepot22
      implicit none

c MISC VARIABLES
      integer i,j,k,nt
      double precision r(nt,3),e(nt),r1s,r2s,r3s,cs,cevau,uli2p
      double precision c2a,c2b,c2c,w,cplg2

c INTERACTION REGION VARIABLES 
c     LiH
      double precision r1r3,xelih,s1c,bes1,rcut2del,rcut2rho,s1c3
      double precision s1c1,belihc,bet1,relih,s1c2,t1c1,t1c2
      double precision r1r3al,rcut2
      double precision coul1,exch1,s1,t1,reliht

c     HF
      double precision behf,xehfa,dehfm,yhf,xehfa_180,t2a3,t2a2,bet2
      double precision t2a1,t2a4,dehf,rehf
      double precision s2alp,rehft,bet2_180,s2rho,dehfc
      double precision coul2,exch2,t2,s2,dehfcc,s2cut,xehfm,s2m,
     & s2del,s2cc,s2cc_d,s2cc_r,s2cc_b

c     LiF
      integer nb
      double precision bet3,bet3_180,t3a2,t3a3,t3a4,belihc1
      double precision belifff,gam,belifc2,delifc1,delifc2,belifm
      double precision s3c,xelihf_180,xelif,s3a,delifm,gc2,cut2
      double precision del2_0,del2_180,rho2_0,rho2_180,theta2_0,grho
      double precision theta2_180,beliff,belifmod,belif,xelif_180,relif
      double precision belifc1,t3a1,gswt,delif,rho2,del2,theta2,gdel
      double precision coul3,exch3,s3,t3,s3cut,swit,swit_r,swit_a,swit_d

c LONG RANGE FORCES VARIABLES
c     HF physical properties
      double precision hfdip_y,hfdip_a,hfdip_re,hfdip_m(11),hfdip,
     & hfpol,hfpolperp,hfpolpar,hfie,hfpolavg
      double precision hfqua_y,hfqua_a,hfqua_re,hfqua_m(10),hfqua

c     Li(2p) physical properties
      double precision li2ppol,li2pie,li2pqua,li2ppolpar,
     & li2ppolperp,li2ppolavg

c     LiF physical properties
      double precision lifpol,lifie
      double precision lifdip_y,lifdip_re,lifdip_1,lifdip_2,
     & lifdip_3,lifdip,lifdip_4,lifdip_beta,lifdip_d

C     H physical properties
      double precision hie,hpol

c     misc long-range variables
      double precision rhf0did,rhf0disp,hf_disp,e_LR_HF,e_LR_LiF,
     & rlif0did,rlif0disp,lif_disp,lif_did,hf_qq,hf_dq,rhf0qq,rhf0dq,
     & hf_did
      double precision mass(3),cauang,caudeb,costhetaHF,costhetaLiF,
     & bigRHF,bigRLiF,coschiHF,coschiLiF,temp,r2LRcut_a,r2LRcut_r,
     & r3LRcut_a,r3LRcut_r

c MISC DATA
      data c2a /3.5d0/, c2b /0.27362d0/, c2c /0.15d0/
      data uli2p/1.848d0/, cevau /0.036749309d0/

c INTERACTION REGION DATA
c     LiH
      data s1c1 /4.32258d0/, s1c2 /7.06452d0/, t1c1 /1.1935483d0/
      data t1c2 /13.548387d0/, belihc /1.36d0/, bet1 /2.41319648d0/
      data relih /1.2d0/, r1r3al /1.0d0/, rcut2rho /0.72d0/
      data rcut2del /0.5d0/, s1c3 /14.74194d0/, bes1 /0.90667d0/
      data reliht /1.203323d0/

c     HF
      data rehf /1.733d0/, dehf /6.122d0/
      data t2a1 /1.4242424d0/, t2a2 /14.203323d0/, t2a3 /6.15835777d0/
      data t2a4/4.56304985d0/,bet2/2.4105572d0/,bet2_180/1.4046921d0/
      data rehft /1.61329423d0/, s2rho /1.00293d0/
      data s2alp /4.0899d0/, dehfc/0.0486803519d0/,dehfcc/0.72629521d0/
      data s2cc_b/1.799609d0/,s2cc_r/1.60215d0/,s2cc_d/2.3841642d0/

c     LiF
      data belifc1/0.819648d0/,belifc2/0.7795699d0/,delifc1/2.352884d0/
      data delifc2 /5.5904203d0/, belifff /0.13225806d0/, nb /8/ 
      data gam /4.6304985d0/, grho /3.382209d0/, gdel /0.49472140d0/
      data t3a1 /-.0552298d0/, t3a2 /1.8729228d0/, t3a3 /0.94662756d0/
      data t3a4/0.2994134d0/,bet3/0.9519061d0/,bet3_180/0.4531769d0/
      data rho2_0 /0.2017595d0/, rho2_180 /2.d0/, del2_0 /0.51710655d0/
      data del2_180 /0.43695014d0/, theta2_0 /0.12463343d0/
      data theta2_180 /0.14907135d0/, delif /5.947d0/, relif /2.9553d0/
      data swit_r/3.340762d0/,swit_a/0.56353861d0/,swit_d/1.0243998d0/

c LONG RANGE FORCES DATA
c     HF dipole moment
      data hfdip_a /1.13d0/, hfdip_re /0.924212d0/
      data hfdip_m /0.703661d0,0.516815d0,0.240628d0,-0.430194d0,
     & -2.213088d0,5.535609d0,12.872106d0,-42.060172d0,-17.398065d0,
     & 84.207629d0,-41.961465d0/

c     HF quadrupole moment
      data hfqua_a /1.37d0/, hfqua_re /0.9168d0/
      data hfqua_m /1.6527d0,2.072d0,2.0521d0,3.545d0,
     & -9.8031d0,-19.3808d0,76.4423d0,0.d0,-161.5896d0,105.7460d0/

c     HF pol (au) and ie (eV)
      data hfpolperp/4.59d0/,hfpolpar/5.10d0/,hfie/16.044d0/

c     Li(2p) pol (au) and ie (eV) and quadrupole moment (au)
      data li2ppolpar/131.d0/,li2ppolperp/129.d0/,li2pie/3.544d0/,
     & li2pqua/11.1d0/

c     LiF dipole moment (distance and moment in au)
      data lifdip_re/10.4994d0/,lifdip_1/0.02435d0/,
     & lifdip_2/9.9355d-8/,lifdip_3/0.0015999d0/,lifdip_4/4.471959d0/,
     & lifdip_d/9.30039d0/

c     LiF pol (au) and ie (eV)
      data lifpol/72.9d0/,lifie/11.3d0/

c     H pol (au) and ie (eV)
      data hie/13.598d0/,hpol/4.4997d0/

c     misc long-range
      data rhf0did/6.d0/,rhf0disp/6.d0/,rhf0qq/6.d0/,rhf0dq/6.d0/
c      data rlif0did/7.5d0/,rlif0disp/7.5d0/,r2LRcut_a/2.0d0/,
      data rlif0did/7.d0/,rlif0disp/7.d0/,r2LRcut_a/2.0d0/,
     & r2LRcut_r/3.d0/,r3LRcut_a/2.d0/,r3LRcut_r/4.5d0/
      data cauang /0.52917706d0/,caudeb/2.54177d0/
      data mass /7.016003d0,1.00783d0,18.9984d0/

      save 
      return

c ---------------------------------------------------------------------- 

      entry pot22(r,e,nt)

      do 10 i=1,nt

c PRECOMPUTE THINGS FOR USE LATER
c     cosine of HF bond angle (i.e., Li-F-H angle)
      costhetaHF = (r(i,3)**2-r(i,2)**2-r(i,1)**2)/(-2.d0*r(i,1)*r(i,2))
      temp = mass(3)/(mass(2)+mass(3))*r(i,2)
c     Li-FH Jacobi distance
      bigRHF = dsqrt(r(i,1)**2+temp**2-2.d0*r(i,1)*temp*costhetaHF)
c     cosine of HF Jacobi angle (i.e., Li--HF angle)
      coschiHF = -(r(i,1)**2-bigRHF**2-temp**2)/(-2.d0*bigRHF*temp)
      coschiHF = max(-1.d0,coschiHF)
      coschiHF = min(1.d0,coschiHF)

c     cosine of LiF bond angle (i.e., Li-F-H angle)
      costhetaLiF=(r(i,1)**2-r(i,2)**2-r(i,3)**2)/(-2.d0*r(i,3)*r(i,2))
      temp = mass(1)/(mass(1)+mass(3))*r(i,3)
c     H-LiF Jacobi distance
      bigRLiF = dsqrt(r(i,2)**2+temp**2-2.d0*r(i,2)*temp*costhetaLiF)
c     cosine of LiF Jacobi angle (i.e., H--LiF angle)
      coschiLiF = -(r(i,2)**2-bigRLiF**2-temp**2)/(-2.d0*bigRLiF*temp)
      coschiLiF = min(1.d0,coschiLiF)
      coschiLiF = max(-1.d0,coschiLiF)
c     squared!
      r1s = r(i,1)**2
      r2s = r(i,2)**2
      r3s = r(i,3)**2
c     Li-F-H bondangle
      cs = r3s+r2s-r1s
      cs = cs/(2.0d0*r(i,2)*r(i,3))
      cs = min(cs, 1.0d0)
      cs = max(cs,-1.0d0)

c INTERACTION REGION
c     The LiH singlet
      xelih = dexp(-belihc*(r(i,1)-relih))
      s1 =   s1c1*xelih**2  + s1c2*xelih 
      s1c = s1c3*dexp(-bes1*(r(i,1)-relih))
      r1r3 = r(i,1) - r(i,3) + r1r3al*r(i,2)
      rcut2 = 0.5d0*(1.0d0 - dtanh((r1r3-rcut2rho)/rcut2del))
      s1 = s1 + rcut2 * s1c

c     The LiH triplet
      xelih = dexp(-bet1*(r(i,1)-reliht))
      t1 =   t1c1*xelih**2  + t1c2*xelih 
    
c     The HF singlet
      yhf = r(i,2) - 2.1042d0
      behf = 1.1622d0 -0.025647d0*yhf+0.059062d0*yhf**2
      xehfa = dexp(-behf*(r(i,2)-rehf))

      s2 = dehf*xehfa*(xehfa-2.0d0) + uli2p
      s2cut = (0.5d0 - 0.5d0*tanh((r(i,2)-3.d0)/0.5d0))
      s2cc = s2cc_d*(1.d0-dexp(-s2cc_b*(r(i,2)-s2cc_r)))**2
     & - s2cc_d
      s2 = s2cc + (s2 - s2cc)*s2cut

      dehfm = dehf - (dehfcc + dehfc*0.5d0*(1.0d0-cs))
      xehfm = dexp(-behf*(r(i,2)-rehf))

      s2m = dehfm*xehfm*(xehfm-2.0d0) + uli2p
      s2cut = (0.5d0 - 0.5d0*tanh((r(i,2)-2.2d0)/0.5d0))
      s2m = s2m*s2cut

      s2 = s2m + (s2 - s2m)*0.5d0*(1.d0+dtanh((r(i,3)-s2rho)/s2alp))

c     The HF triplet
      xehfa = dexp(-bet2*(r(i,2)-rehft))
      xehfa_180 = dexp(-bet2_180*(r(i,2)-rehft))

      t2  = (t2a1*xehfa**2 + t2a2*xehfa)*0.5d0*(1.0d0+cs)
     &    + (t2a3*xehfa_180**2 + t2a4*xehfa_180)*0.5d0*(1.0d0-cs)
 
c     The LiF singlet
      belif = r(i,3)*103.57d0/(4.6498d0+7.0489d0*r(i,3))**2+0.064076d0

      gswt = 0.5d0-0.5d0*dtanh((r(i,2)-grho)/gdel)
      gswt = gswt * 0.5d0*(1.0d0+cs)
      beliff = belif + gswt*(belifff - belif)

      belifmod = beliff*(belif +(r(i,3)/gam)**nb)/
     &                  (beliff+(r(i,3)/gam)**nb)

      xelif = dexp(-belifmod*(r(i,3)-relif))
      s3a = (delif+uli2p)*xelif*(xelif-2.0d0)+uli2p
      s3cut = 0.5d0 - 0.5d0*tanh((r(i,3)-13.d0)/0.5d0)
      s3a = s3a*s3cut

      swit = swit_d  + (1.d0 - swit_d)*0.5d0
     &  *(1.d0 + dtanh((r(i,2)-swit_r)/swit_a))
      s3a = s3a*swit

      delifm = delifc1 * 0.5d0*(1.0d0-cs)
     &  + delifc2 * 0.5d0*(1.0d0+cs)

      belifm = belifc1 * 0.5d0*(1.0d0-cs)
     &  + belifc2 * 0.5d0*(1.0d0+cs)

      xelif = dexp(-belifm*(r(i,3)-relif))
      s3c = (delifm+uli2p)*xelif*(xelif-2.0d0)

      rho2 = rho2_0 + (rho2_180-rho2_0)*0.5d0*(1.0d0-cs)
      theta2 = theta2_0 + (theta2_180-theta2_0)*0.5d0*(1.0d0-cs)
      del2 = del2_0 + (del2_180-del2_0)*0.5d0*(1.0d0-cs)
 
      gc2 = dcos(theta2)*(rehf-r(i,2)) - dsin(theta2)*(rho2-r(i,3))
      cut2 = 0.5d0 - 0.5d0*dtanh(gc2/del2)

      s3 = s3c + (s3a-s3c)*cut2

c     The LiF triplet
      xelif = dexp(-bet3*(r(i,3)-relif))
      xelif_180 = dexp(-bet3_180*(r(i,3)-relif))
      t3  = (t3a1*xelif**2 + t3a2*xelif)*0.5d0*(1.0d0+cs)
     &    + (t3a3*xelif_180**2 + t3a4*xelif_180)*0.5d0*(1.0d0-cs)

c     LEPS
      coul1=0.5d0*(s1+t1)
      coul2=0.5d0*(s2+t2)
      coul3=0.5d0*(s3+t3)
      exch1=0.5d0*(s1-t1)
      exch2=0.5d0*(s2-t2)
      exch3=0.5d0*(s3-t3)
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      cplg2=c2a*dexp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
      e(i)=coul1+coul2+coul3-dsqrt(w+cplg2*cplg2)/dsqrt(2.0d0) + 
     &   dehf

c LONG-RANGE FORCES * * * * * * * * * * * * * * * * * * * * * * * * * * *
c HF DIPOLE
      hfdip_y = 1.d0 - dexp(-hfdip_a*(r(i,2)*cauang-hfdip_re))
      hfdip = 0.d0
      do j=1,11
      hfdip = hfdip + hfdip_m(j)*hfdip_y**dble(j-1)
      enddo
      hfdip = hfdip*0.5d0*(1.d0 + dtanh(20.d0*(r(i,2)-1.2d0)))

c HF POLARIZATION
      hfpol = hfpolperp*(1.d0 - coschiHF**2) + hfpolpar*coschiHF**2
      hfpolavg = 2.d0*hfpolperp/3.d0 + hfpolpar/3.d0

c HF QUADRUPOLE
      hfqua_y = 1.d0 - dexp(-hfqua_a*(r(i,2)*cauang-hfqua_re))
      hfqua = 0.d0
      do j=1,10
      hfqua = hfqua + hfqua_m(j)*hfqua_y**dble(j-1)
      enddo
      hfqua = hfqua*0.5d0*(1.d0 + dtanh(20.d0*(r(i,2)-1.2d0)))

c LiF DIPOLE
      lifdip_y = r(i,3) - lifdip_4
      lifdip_beta = lifdip_1 + lifdip_2*lifdip_y**8 + lifdip_3*lifdip_y
      lifdip = lifdip_d*dexp(-lifdip_beta*(r(i,3) - lifdip_re)**2)

c Li(2p) POLARIZATION
      li2ppol = li2ppolperp*(1.d0-coschiHF**2)+li2ppolpar*coschiHF**2
      li2ppolavg = 2.d0*li2ppolperp/3.d0 + li2ppolpar/3.d0

c HF + Li(2p)
c     DIPOLE-INDUCED-DIPOLE
      hf_did = -0.5d0*(3.d0*coschiHF**2 + 1.d0)*li2ppol*hfdip**2
      hf_did = hf_did/(bigRHF**7 + rhf0did**7)/cevau*bigRHF
c     cut off angular dependence in interaction region
      hf_did = hf_did*(1.d0 + (2.5d0*li2ppolavg/
     &        (3.d0*coschiHF**2 + 1.d0)/li2ppol - 1.d0)
     &               *dexp(-(bigRHF/rhf0did)**4))

c     DISPERSION
      hf_disp = -1.5d0*li2pie*hfie*li2ppol*hfpol/(li2pie+hfie)
      hf_disp = hf_disp/(bigRHF**7 + rhf0disp**7)*bigRHF
c     cut off angular dependence in interaction region
      hf_disp = hf_disp*(1.d0+(li2ppolavg*hfpolavg/hfpol/li2ppol-1.d0)
     &                 *dexp(-(bigRHF/rhf0disp)**4))

c     QUADRUPOLE-QUADRUPOLE
      hf_qq = 3.d0/4.d0*hfqua*li2pqua*(3.d0 - 7.d0*coschiHF**2)
      hf_qq = hf_qq/(bigRHF**6 + rhf0qq**6)/cevau*bigRHF
c     cut off angular dependence in interaction region
      hf_qq = hf_qq*(1.d0 + (-.5d0/(3.d0-7.d0*coschiHF**2) - 1.d0)
     &          *dexp(-(bigRHF/rhf0qq)**4))

c     DIPOLE-QUADRUPOLE
      hf_dq = 3.d0/4.d0*hfdip*li2pqua*coschiHF
      hf_dq = hf_dq/(bigRHF**5 + rhf0dq**5)/cevau*bigRHF
c     cut off angular dependence in interaction region
      hf_qq = hf_qq*dexp(-(bigRHF/rhf0dq)**4)

c     SUM FORCES
      e_LR_HF = hf_dq + hf_qq + hf_disp + hf_did
c     cut off at large HF distances
      e_LR_HF=e_LR_HF*0.5d0*(1.d0-dtanh(r2LRcut_a*(r(i,2)-r2LRcut_r)))

c LiF + H
c     DIPOLE-INDUCED-DIPOLE
      lif_did = -0.5d0*(3.d0*coschiLiF**2+1.d0)*hpol*lifdip**2
      lif_did = lif_did/(bigRLiF**6 + rlif0did**6)/cevau
c     cut off angular dependence in interaction region
      lif_did = lif_did*(1.d0+(2.5d0/(3.d0*coschiLiF**2+1.d0)-1.d0)
     &           *dexp(-(bigRLiF/rlif0did)**4))

c     DISPERSION
      lif_disp = -1.5d0*lifie*hie*lifpol*hpol/(lifie+hie)
      lif_disp = lif_disp/(bigRLiF**6 + rlif0disp**6)
c     no angular dependence

c     SUM FORCES
      e_LR_LiF = lif_disp + lif_did
c     cut off at large LiF distances
      e_LR_LiF=e_LR_LiF*0.5d0*(1.d0-dtanh(r3LRcut_a*(r(i,3)-r3LRcut_r)))

c LiH + F
c no long-range forces included

c CUT OFF LONG-RANGE FORCES
      if (bigRHF.gt.15.d0) then
        e_LR_HF = 0.d0
      else
        e_LR_HF = e_LR_HF * dexp( .2d0 / (bigRHF - 15.d0) )
      endif

      if (bigRLiF.gt.15.d0) then
        e_LR_LiF = 0.d0
      else
        e_LR_LiF = e_LR_LiF * dexp( .2d0 / (bigRLiF - 15.d0) )
      endif

c COMBINE INTERACTION REGION ENERGY AND LONG-RANGE FORCES
      e(i)=e(i) + e_LR_LiF + e_LR_HF
      e(i)=e(i)*cevau

 10   continue
      end

c***************************************************************************
c**************************************************************************
