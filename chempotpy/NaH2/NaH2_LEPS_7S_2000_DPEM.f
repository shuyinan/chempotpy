      subroutine dpem(x,igrad,u,ug)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: u(nstates,nstates)
      double precision, intent(out) :: ug(nstates,nstates,natoms,3)

      double precision :: nt, r(1,3), v(1)
      double precision :: u11(1), u12(1), u22(1)
      logical, save :: first_time_data=.true.
      integer :: i, j, k, l

      !initialize 
      u=0.d0
      ug=0.d0

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
      else
        write(*,*) "Only DPEM is available"
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

c   System:             NaH2
C   Functional form:    Modified Extended LEPS (2x2 diabatic fit)
c   Common name:        NaH2 surface set 7s
c   Interface:          3-1V
c   Number of electronic surfaces: 2
c   Number of derivatives: 0

c   Notes:               This version does not contain long range forces,
c                        and has no derivatives.  It is based on the surface
c                        set 6 fit.

C   References:         Unpublished, M. D. Hack and D. G. Truhlar (2000).

C   Protocol:
C      PREPOT - initializes the potentials variables
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
C   and the potential coupling.
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
C                            R(2) = R(H-H)
C                            R(3) = R(Na-H)

c
c This is NaH2 potential matrix 7s, finalized on 5/17/99 
c This version does not contain any long range forces.

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
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank, g_nprocs

        call prepot11
        call prepot12
        call prepot22

      return

      end
c
c ********************************************************************
c ********************************************************************
      subroutine pot(r,e,nt,nsurf)

        implicit double precision (a-h,o-z)
        dimension r(nt,3),e(nt)

        if (nsurf .eq. 1) then
          call pot11(r,e,nt)
        else if (nsurf .eq. 2) then
          call pot12(r,e,nt)
        else if (nsurf .eq. 3) then
          call pot22(r,e,nt)
        else
          do i=1,nt
            call pot11(r(i,1),v11,1)
            call pot12(r(i,1),v12,1)
            call pot22(r(i,1),v22,1)
            if (nsurf .eq. 4) then
              e(i)=0.5d0*(v11+v22-dsqrt((v11-v22)**2+4.d0*v12**2))
            else
              e(i)=0.5d0*(v11+v22+dsqrt((v11-v22)**2+4.d0*v12**2))
            end if
          enddo
        end if

      return

      end

c ********************************************************************
c ********************************************************************
c  U11 surface - lower diabatic surface
c
      subroutine prepot11
      implicit double precision (a-h,o-z)

      dimension r(nt,3),e(nt)

      save
c NaH singlet
      data des1/0.55118d0/, bes1 /1.05118d0/, res1 /3.38583d0/

c H2 singlet
      data reh2 /1.401d0/, bf /1.627d0/, b0 /1.194d0/, 
     &     gamh2a /2.323d0/, gamh2b /3.126d0/, deh2/4.74806d0/
      data des20 /0.996d0/, s2alp /3.0d0/, s2rho /7.0d0/

c H-Na-H term
      data dect /2.24409d0/, bect /0.57874d0/, 
     &     rect /3.46457d0/, bect2 /0.87143d0/

c  misc constants

      rac2=1.d0/dsqrt(2.d0)
      cevau=1.d0/27.2113961d0
      onethd=1.0d0/3.0d0
      twothd=2.0d0/3.0d0

      return

      entry pot11(r,e,nt)

      do 10 i=1,nt


c ============================================================
c useful quantities calculated first

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)

      rmid = 0.5d0*(r(i,1)+r(i,3))
      rave = 0.5d0*(r(i,1)-r(i,3))

c cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up

      cos1 = (r1s+r2s-r3s)/(2.0d0*r(i,1)*r(i,2))
      cos2 = (r3s+r2s-r1s)/(2.0d0*r(i,3)*r(i,2))
      cosg = 0.5d0*(cos1-cos2)

      cosa=(r1s+r3s-r2s)/(2*r(i,1)*r(i,3))

c==============================================================k
cThis section contains the diatomic curves

c The first NaH singlet

      xes1 = dexp(-bes1*(r(i,1)-res1))
      s1 = des1*onethd*xes1*(xes1+2.0d0)

c The H2 singlet

      reh2 = 1.401d0
      bf = 1.627d0
      b0 = 1.194d0
      gamh2a = 2.323d0
      gamh2b = 3.126d0

      beh2 = bf*(b0-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2)
     &   /(bf-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2)

      sswt = 0.5d0*(1.0d0-dtanh((rmid-s2rho)/s2alp))*cosg**2
      deh2m = deh2 + deh2*(des20 - 1.0d0)*sswt

      xeh2 = dexp(-beh2*(r(i,2)-reh2))
      s2 = deh2m*xeh2*(xeh2-2.0d0)

c The second NaH singlet

      xes3 = dexp(-bes1*(r(i,3)-res1))
      s3 = des1*onethd*xes3*(xes3+2.0d0)

c =========================================================
c here is where it all gets put together . . .

      e(i) = s1 + s2 + s3 + deh2

c here we will try a post-LEPS correction term
      xect = dexp(-bect*(rmid-rect))
      vmid = dect*xect*(xect-2.0d0)
     &    *dexp(-bect2*rave**2)*(0.5d0*(1.0d0-cosa))**2

      e(i)=e(i)+vmid

      e(i) = e(i)*cevau

10    continue
      return

      end


c ********************************************************************
c ********************************************************************

c U12 coupling surface . . .

      subroutine prepot12

      implicit none

      double precision ap,alp,bet,rnagsq,rnag,cosg,rd0,r20
      double precision r1s,r2s,r3s,cos1,cos2,theta
      double precision gc, gc2

      double precision cevau

      integer nt,i,j,k
      logical leven
      double precision r(nt,3),e(nt)

      data  ap/0.835d0/,alp/4.0d0/,rd0/4.0d0/,bet/0.5d0/,r20/-1.3d0/

      save

      theta = 15.6d0*dacos(-1.0d0)/180.0d0


      cevau=1.d0/27.2113961d0

c set leven true for an even U12 surface, and false for an odd 
c surface
      leven=.true.
      return

      entry pot12(r,e,nt)

      do i=1,nt

      r1s=r(i,1)**2
      r2s=r(i,2)**2
      r3s=r(i,3)**2

c here we do the function evaluations . . . 

      rnagsq=0.5d0*r1s+0.5d0*r3s-0.25d0*r2s
      if (rnagsq .gt. 0.0d0) then
        rnag=dsqrt(rnagsq)
      else
        rnag=0.0d0
      endif

c cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up

      cos1 = (r1s+r2s-r3s)/(2.0d0*r(i,1)*r(i,2))
      cos2 = (r3s+r2s-r1s)/(2.0d0*r(i,3)*r(i,2))
      cosg = 0.5d0*(cos1-cos2)
 
      gc  = dcos(theta)*rnag + dsin(theta)*r(i,2)
      gc2 = dsin(theta)*rnag - dcos(theta)*r(i,2)
      e(i) = ap*cosg
     &   *dexp(-(1.0d0/alp)**4*(gc - rd0)**4)
     &   *dexp(-(1.0d0/bet)**2*(gc2 - r20)**2)

      e(i) = e(i)*cevau 

c set leven equal to true for an even surface, false for an odd one
      if (leven) then
        e(i)=dsqrt(e(i)**4/(e(i)**2+1.0d-8))
      endif

      enddo
      return
      end

c ********************************************************************
c ********************************************************************

c the upper diabatic surface  U22

      subroutine prepot22
      implicit double precision (a-h,o-z)

      dimension r(nt,3),e(nt)

      save

c NaH triplets
      data ret10 /3.22835d0/, bet10 /0.57874d0/,
     & bet10b /2.44094d0/, f10 /9.44882d0/
      data det190 /4.32441d0/, bet190 /2.43701d0/, ret190 /1.95480d0/,
     &  t1alp /0.30778d0/, t1rho /3.21850d0/
      data det1c /4.92126d0/, bet1c /0.65748d0/, ret1c /2.71654d0/

c H2 triplets
      data det20 /3.29921d0/, bet20 /2.64567d0/, ret20 /2.44882d0/
      data det290 /3.97638d0/, bet290 /3.73937d0/, ret290 /0.96024d0/
      data t2alp /1.0d0/, t2rho /9.0d0/

c H2 singlets
      data des290 /0.91433d0/, s2alp /0.66778d0/, s2rho /5.30472d0/
      data deh2 /4.74806d0/
      data s20alp /0.5d0/, s20rho /6.0d0/

c misc constants

      data c2a /1.38661d0/, c2b /0.27362d0/, c2c /0.03764d0/, 
     &   una2p/2.1037d0/


      cevau=1.d0/27.2113961d0
      onethd=1.0d0/3.0d0
      twothd=2.0d0/3.0d0

      return

      entry pot22(r,e,nt)

      do 10 i=1,nt

c ============================================================
c useful quantities calculated first

      rmid = 0.5d0*(r(i,1)+r(i,3))
      rave = 0.5d0*(r(i,1)-r(i,3))

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)


      cosa=(r1s+r3s-r2s)/(2*r(i,1)*r(i,3))

c cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up

      cos1 = (r1s+r2s-r3s)/(2.0d0*r(i,1)*r(i,2))
      cos2 = (r3s+r2s-r1s)/(2.0d0*r(i,3)*r(i,2))
      cosg = 0.5d0*(cos1-cos2)

c =========================================================
c The diatomic curves . . .

c The first NaH singlet

      denah = 1.97134878d0
      renah = 3.566044d0
      bf = 0.864d0
      b0 = 0.594d0
      gamnah = 7.376d0

      beta = bf*(b0+(r(i,1)/gamnah)**8)/(bf+(r(i,1)/gamnah)**8)
      xenah = dexp(-beta*(r(i,1)-renah))

      s1 = denah*xenah*(xenah-2.0d0)

c The first NaH triplet

      g10 = r(i,1)-ret10
      t10 = f10*dexp(-bet10*g10)+dexp(-bet10b*g10)
      t10 = t10/(1.0d0+f10)

      xet190 = dexp(-bet190*(r(i,1)-ret190))
      t190 = det190*onethd*xet190*(xet190 + 2.0d0)

      xet1c = dexp(-bet1c*(r(i,1)-ret1c))
      t1c = det1c*onethd*xet1c*(xet1c + 2.0d0)

      t1m = t190 +(t1c-t190)*
     &       0.5d0*(1.0d0+dtanh(t1alp*(r(i,2)-t1rho)))

      t1 = t1m + (t10-t1m)*cosg**2


c The H2 triplet

      xet20 = dexp(-bet20*(r(i,2)-ret20))
      t20 = det20*onethd*xet20*(xet20 + 2.0d0)

      xet290 = dexp(-bet290*(r(i,2)-ret290))
      t290 = det290*onethd*xet290*(xet290 + 2.0d0)

      t2 = t20 + (t290-t20)*(1.0d0-cosg**2)*
     &   0.5d0*(1.0d0-dtanh(t2alp*(rmid-t2rho)))

c The H2 singlet
c new singlet, designed to minic lower adiabatic behavior
      reh2 = 1.401d0
      bf = 1.627d0
      b0 = 1.194d0
      gamh2a = 2.323d0
      gamh2b = 3.126d0

      beh2 = bf*(b0-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2)
     &   /(bf-(r(i,2)/gamh2a)+(r(i,2)/gamh2b)**2)

      sswt = 0.5d0*(1.0d0-dtanh(s2alp*(rmid-s2rho)))
      deh2m = deh2 + deh2*(des290 - 1.0d0)*sswt

      xeh2 = dexp(-beh2*(r(i,2)-reh2))
      s2a = deh2m*xeh2*(xeh2-2.0d0) + una2p

      s2=0.5d0*(s2a+t2-dtanh((s2a-t2)/0.1d0)*(s2a-t2))

      xeh2 = dexp(-1.34d0*beh2*(r(i,2)-reh2))
      s20b = (deh2-una2p)*xeh2*(xeh2-2.0d0)

      s2 = s2 + (s20b - s2)*cosg**2*
     &  0.5d0*(1.0d0-dtanh(s20alp*(rmid-s20rho)))

c The second NaH singlet
      denah = 1.97134878d0
      renah = 3.566044d0
      bf = 0.864d0
      b0 = 0.594d0
      gamnah = 7.376d0

      beta = bf*(b0+(r(i,3)/gamnah)**8)/(bf+(r(i,3)/gamnah)**8)
      xenah = dexp(-beta*(r(i,3)-renah))

      s3 = denah*xenah*(xenah-2.0d0)

c The second NaH triplet

      xet3 = dexp(-bet1*(r(i,3)-ret1))
      t30 = det1*onethd*xet3*(xet3 + 2.0d0)

      g30 = r(i,3)-ret10
      t30 = f10*dexp(-bet10*g30)+dexp(-bet10b*g30)
      t30 = t30/(1.0d0+f10)

      xet390 = dexp(-bet190*(r(i,3)-ret190))
      t390 = det190*onethd*xet390*(xet390 + 2.0d0)

      xet3c = dexp(-bet1c*(r(i,3)-ret1c))
      t3c = det1c*onethd*xet3c*(xet3c + 2.0d0)

      t3m = t390 + (t3c - t390)*
     &       0.5d0*(1.0d0+dtanh(t1alp*(r(i,2)-t1rho)))

      t3 = t3m + (t30 - t3m)*cosg**2


c ===============================================================
c The final LEPS form

      coul1=0.5d0*(s1+t1)
      coul2=0.5d0*(s2+t2)
      coul3=0.5d0*(s3+t3)
      exch1=0.5d0*(s1-t1)
      exch2=0.5d0*(s2-t2)
      exch3=0.5d0*(s3-t3)
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      cplg2=c2a*dexp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
      e(i)=coul1+coul2+coul3+deh2-dsqrt(w+cplg2*cplg2)/dsqrt(2.0d0)

      e(i)=e(i)*cevau

 10   continue
      end



c***************************************************************************
c**************************************************************************
