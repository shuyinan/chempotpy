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

      subroutine prepot
C
C
C   System:          NaH2
C   Functional form:    Modified Extended LEPS (2x2 diabatic fit)
C   Common name:     NaH2 surface set 5F
C   Interface:       3-1V
C   Number of electronic surfaces: 2
C   Number of derivatives: 0

C   Reference:       P. Halvick and D. G. Truhlar, J. Chem. Phys. 96, 
C                    2895 (1992); E 100, 4718 (1994)
C
C   Note:            This is a vectorized version of the 5F surface
C                    without derivatives.
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
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
C      The zero of energy occurs on the lower adiabatic surface when 
C      the Na atom is "infinitely" far from the H2 diatom and the
C      H2 diatom is at its equilibrium geometry.  On this surface this occurs
C      at R(H-H)=1.402721754031838 bohr with an energy of 1.842768357662727E-4
C      hartrees (=5.014430 meV)
C
C
C
        call prepot11
        call prepot12
        call prepot22

      return

      end
c
      subroutine pot(r,e,nt,nsurf)

        implicit double precision (a-h,o-z)

        dimension r(nt,3),e(nt)

        if (nsurf .eq. 1) then
          call pot11(r,e,nt)
        else if (nsurf .eq. 2) then
          call pot12(r,e,nt)
        else
          call pot22(r,e,nt)
        end if

      return

      end
c
c  U11 surface - lower diabatic surface
c
      subroutine prepot11
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt)

      save
      data dt2/2.885656d-2/,ret2/1.237624d0/,betat2/1.038036d+1/, 
     &     y90c1/6.279847d0/,y90c2/3.849003d-3/,y90c3/5.546257d0/,
     &     y90c4/5.573967d-1/,y0c1/1.060311d+1/,y0c2/6.058605d-3/,
     &     y0c3/4.573378d0/,y0c4/5.182708d-1/
      data ads1/0.08d0/,ars1/5.015440d0/,abs1/7.536929d-1/,
     &     cds1/1.323042d0/,crs1/4.00d0/,cbs1/5.4070d-1/
      data swa/0.7d0/,swr0/6.5d0/,swa2/1.1d0/,swr20/1.9d0/,
     &     c2a/1.5d0/,c2b/1.0d0/,c2c/0.15d0/,deh2/4.7469d0/

      rac2=1.d0/dsqrt(2.d0)
      cevau=1.d0/27.211611d0

      return
c
      entry pot11(r,e,nt)

      do 10 i=1,nt

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)
      rnag=dsqrt(0.5d0*dabs(r1s+r3s-0.5d0*r2s))
      if (rnag.lt.1.d-10) stop 'RNAG < 1.D-10'
      if (r(i,2).ne.0.d0) then
        cstsq=(0.5d0*(r1s-r3s)/(r(i,2)*rnag))**2
      else
        cstsq=0.d0
      end if
      c1=y90c1+(y0c1-y90c1)*cstsq
      c2=y90c2+(y0c2-y90c2)*cstsq
      c3=y90c3+(y0c3-y90c3)*cstsq
      c4=y90c4+(y0c4-y90c4)*cstsq

      swt=0.5d0*(1.d0-dtanh(swa*(r(i,1)-swr0)))
     &    *0.5d0*(1.d0-dtanh(swa*(r(i,3)-swr0)))
     &    *0.5d0*(1.d0+dtanh(swa2*(r(i,2)-swr20)))

      axs1=dexp(-abs1*(r(i,1)-ars1))
      as1=ads1*axs1*(axs1-2.d0)
      cxs1=dexp(-cbs1*(r(i,1)-crs1))
      cs1=cds1*cxs1*(cxs1-2.d0)
      s1=as1+(cs1-as1)*swt

      s2=(149.480d0-59.6557d0*r(i,2))*dexp(-4.13792d0*r(i,2))
     &   +r2s*(-23.7299d0+3.91747d0*r(i,2))*dexp(-1.41350d0*r(i,2))

      axs3=dexp(-abs1*(r(i,3)-ars1))
      as3=ads1*axs3*(axs3-2.d0)
      cxs3=dexp(-cbs1*(r(i,3)-crs1))
      cs3=cds1*cxs3*(cxs3-2.d0)
      s3=as3+(cs3-as3)*swt

      t1=c1*dexp(-c2*r(i,1)**6)+c3*dexp(-c4*r(i,1))

      xt2=dexp(-betat2*(r(i,2)-ret2))
      t2=dt2*xt2*(xt2-2.d0)

      t3=c1*dexp(-c2*r(i,3)**6)+c3*dexp(-c4*r(i,3))

      coul1=0.5d0*(s1+t1)
      coul2=0.5d0*(s2+t2)
      coul3=0.5d0*(s3+t3)
      exch1=0.5d0*(s1-t1)
      exch2=0.5d0*(s2-t2)
      exch3=0.5d0*(s3-t3)
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      cplg2=c2a*dexp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
      e(i)=(coul1+coul2+coul3+deh2-rac2*dsqrt(w+cplg2*cplg2))*cevau

10    continue

      return

100   format(/,24('='),' NaH2 - U11 potential ',24('='),
     &  /'   dt2 =',1pe14.6,2x,'  ret2 =',e14.6,2x,'betat2 =',e14.6,
     &  /' y90c1 =',e14.6,2x,' y90c2 =',e14.6,2x,' y90c3 =',e14.6,
     &  /' y90c4 =',e14.6,2x,'  y0c1 =',e14.6,2x,'  y0c2 =',e14.6,
     &  /'  y0c3 =',e14.6,2x,'  y0c4 =',e14.6,2x,'  ads1 =',e14.6,
     &  /'  ars1 =',e14.6,2x,'  abs1 =',e14.6,2x,'  cds1 =',e14.6,
     &  /'  crs1 =',e14.6,2x,'  cbs1 =',e14.6,2x,'   swa =',e14.6,
     &  /'  swr0 =',e14.6,2x,'  swa2 =',e14.6,2x,' swr20 =',e14.6,
     &  /'   c2a =',e14.6,2x,'   c2b =',e14.6,2x,'   c2c =',e14.6,
     &  /70('='))

      end
c
c  U12 surface - coupling surface
c
      subroutine prepot12
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt)

      save
      data eps/1.d-4/
      data ampn,ampp,ampe,ampr,erln,erlp,erle,erlr,
     &     ersn,ersp,erse,ersr,rmxn,rmxp,rmxe,rmxr /
     &  2.35737419d-1,1.51495000d+0,5.02447648d-1, 4.05991868d+0,
     &  2.30065420d+2,7.10000000d-2,3.46985415d+0, 8.95373168d-1,
     &  2.00000000d+0,1.16830000d+1,4.54595148d+1,-2.38192975d+0,
     &  2.33484458d+0,2.62222000d+0,2.04523701d+1, 1.98941949d+0 /       1G25T93

      cevau=1.d0/27.211611d0

      return
c
      entry pot12(r,e,nt)

      do 10 i=1,nt

      amp=ampn+(ampp-ampn)*0.5d0*(dtanh(ampe*(r(i,2)-ampr))+1.d0)
      erl=erln+(erlp-erln)*0.5d0*(dtanh(erle*(r(i,2)-erlr))+1.d0)
      ers=ersn+(ersp-ersn)*0.5d0*(dtanh(erse*(r(i,2)-ersr))+1.d0)
      rmx=rmxn+(rmxp-rmxn)*0.5d0*(dtanh(rmxe*(r(i,2)-rmxr))+1.d0)
      f1=amp*dexp(-(erl+ers/(eps+r(i,1)**3))*(r(i,1)-rmx)**2)
      f3=amp*dexp(-(erl+ers/(eps+r(i,3)**3))*(r(i,3)-rmx)**2)
      ang=((r(i,1)-r(i,3))/r(i,2))**2

      e(i)=ang*(f1+f3)*cevau

10    continue

      return

100   format(/,28('='),' NaH2 - U12 coupling ',29('='),
     &  /'ampn=',1pe13.5,'  ampp=',e13.5,'  ampe=',e13.5,'  ampr=',e13.5,
     &  /'erln=',e13.5,'  erlp=',e13.5,'  erle=',e13.5,'  erlr=',e13.5,
     &  /'ersn=',e13.5,'  ersp=',e13.5,'  erse=',e13.5,'  ersr=',e13.5,
     &  /'rmxn=',e13.5,'  rmxp=',e13.5,'  rmxe=',e13.5,'  rmxr=',e13.5,
     &  /78('='))                                                        1G25T93

      end
c
c  U22 surface - upper diabatic surface
c
      subroutine prepot22
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt)

      save
      data c1s1/2.614452d+2/,c2s1/-3.155166d-4/,c3s1/-1.785467d+0/,
     &     c4s1/-3.548074d+0/,c5s1/3.356766d-1/,c6s1/-7.621511d-1/
      data deh2/4.7469d0/,una2p/2.1037d0/,c2a/0.8d0/,c2b/0.8d0/,
     &     c2c/0.15d0/,cprime/1.0d-1/
      data y0d1/2.417800d-1/,y0r1/3.708055d0/,y0b1/6.960127d-1/,
     &     y90d1/1.372065d-2/,y90r1/7.233552d0/,y90b1/5.127734d-1/,
     &     y0d2/8.870172d-1/,y0r2/2.300265d0/,y0b2/8.497374d-1/,
     &     y90d2/1.106187d0/,y90r2/2.002516d0/,y90b2/9.514887d-1/,
     &     y0swa/0.614393d0/,y0swr0/9.006323d0/,
     &     y90swa/0.872116d0/,y90swr0/8.194467d0/ 
      data a0g1/-8.210245d-2/,b0g1/4.658698d0/,r0g1/4.552860d0/,
     &     a0g2/3.546427d-2/,b0g2/6.456787d-1/,r0g2/6.332225d0/,
     &     a0g3/-1.562327d-2/,b0g3/1.170537d-1/,r0g3/1.097310d+1/
      data a90g1/-9.929297d-2/,b90g1/4.367220d0/,r90g1/4.402447d0/,
     &     a90g2/4.696181d-2/,b90g2/5.470459d-1/,r90g2/6.289526d0/,
     &     a90g3/-1.219900d-2/,b90g3/1.002291d-2/,r90g3/1.118725d+1/

      rac2=1.d0/dsqrt(2.d0)
      cevau=1.d0/27.211611d0

      return
c
      entry pot22(r,e,nt)

      do 10 i=1,nt

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)
      rnag=dsqrt(0.5d0*dabs(r1s+r3s-0.5d0*r2s))
      if (rnag.lt.1.d-10) stop 'RNAG < 1.D-10'
      cstsq=(0.5d0*(r1s-r3s)/(r(i,2)*rnag))**2

      swa=y90swa+(y0swa-y90swa)*cstsq
      swr0=y90swr0+(y0swr0-y90swr0)*cstsq

      s1=(c1s1+c2s1*r(i,1))*dexp(c3s1*r(i,1))
     &   +r1s*(c4s1+c5s1*r(i,1))*dexp(c6s1*r(i,1))

      x90t2=dexp(-y90b2*(r(i,2)-y90r2))
      t902a=y90d2*x90t2*(x90t2-2.d0)
      x0t2=dexp(-y0b2*(r(i,2)-y0r2))
      t02a=y0d2*x0t2*(x0t2-2.d0)
      t2a=t902a+(t02a-t902a)*cstsq
      swt=0.5d0*(1.d0-dtanh( swa*(0.5d0*(r(i,1)+r(i,3))-swr0)))
      s2a=(149.480d0-59.6557d0*r(i,2))*dexp(-4.13792d0*r(i,2))+r2s
     &    *(-23.7299d0+3.91747d0*r(i,2))*dexp(-1.41350d0*r(i,2))+una2p
      s2b=144.893d0*dexp(-3.85716d0*r(i,2))+r(i,2)*(37.5919d0
     &    -4.32985d0*r(i,2)-0.003807d0*r2s)*dexp(-1.52496d0*r(i,2))
      t2=swt*(t2a-s2b)+s2b
      s2=0.5d0*(s2a+t2-dtanh((s2a-t2)/cprime)*(s2a-t2))

      s3=(c1s1+c2s1*r(i,3))*dexp(c3s1*r(i,3))
     &   +r3s*(c4s1+c5s1*r(i,3))*dexp(c6s1*r(i,3))

      x90t1=dexp(-y90b1*(r(i,1)-y90r1))
      t901=y90d1*x90t1*(x90t1+2.d0)
      x0t1=dexp(-y0b1*(r(i,1)-y0r1))
      t01=y0d1*x0t1*(x0t1+2.d0)
      t1=t901+(t01-t901)*cstsq

      x90t3=dexp(-y90b1*(r(i,3)-y90r1))
      t903=y90d1*x90t3*(x90t3+2.d0)
      x0t3=dexp(-y0b1*(r(i,3)-y0r1))
      t03=y0d1*x0t3*(x0t3+2.d0)
      t3=t903+(t03-t903)*cstsq

      coul1=0.5d0*(s1+t1)
      coul2=0.5d0*(s2+t2)
      coul3=0.5d0*(s3+t3)
      exch1=0.5d0*(s1-t1)
      exch2=0.5d0*(s2-t2)
      exch3=0.5d0*(s3-t3)
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      cplg2=c2a*dexp(-c2b*w-c2c*(r(i,1)+r(i,2)+r(i,3)))
      e(i)=coul1+coul2+coul3+deh2-rac2*dsqrt(w+cplg2*cplg2)

      g1 = a0g1*dexp(-b0g1*(rnag-r0g1)**2)
      g2 = a0g2*dexp(-b0g2*(rnag-r0g2)**2)
      g3 = a0g3*dexp(-b0g3*(rnag-r0g3)**2)
      g4 = a90g1*dexp(-b90g1*(rnag-r90g1)**2)
      g5 = a90g2*dexp(-b90g2*(rnag-r90g2)**2)
      g6 = a90g3*dexp(-b90g3*(rnag-r90g3)**2)
      vdisp = g4+g5+g6+(g1+g2+g3-g4-g5-g6)*cstsq
      vdisp = vdisp*dexp(-0.5d0*(r(i,2)/2.d0)**7)

      e(i)=(e(i)+vdisp)*cevau

10    continue

      return

100   format(/,24('='),' NaH2 - U22 potential ',25('='),
     &    /'  y90d1 =',1pe14.6,'    y90r1 =',e14.6,'  y90b1 =',e14.6,
     &    /'  y90d2 =',e14.6,  '    y90r2 =',e14.6,'  y90b2 =',e14.6,
     &    /' y90swa =',e14.6,  '  y90swr0 =',e14.6,
     &    /'   y0d1 =',1pe14.6,'     y0r1 =',e14.6,'   y0b1 =',e14.6,
     &    /'   y0d2 =',e14.6,  '     y0r2 =',e14.6,'   y0b2 =',e14.6,
     &    /'  y0swa =',e14.6,  '   y0swr0 =',e14.6,
     &    /'   c1s1 =',1pe14.6,'     c2s1 =',e14.6,'   c3s1 =',e14.6,
     &    /'   c4s1 =',1pe14.6,'     c5s1 =',e14.6,'   c6s1 =',e14.6,
     &    /'    c2a =',1pe14.6,'      c2b =',e14.6,'    c2c =',e14.6,
     &    /'   a0g1 =',1pe14.6,'     b0g1 =',e14.6,'   r0g1 =',e14.6,
     &    /'   a0g2 =',1pe14.6,'     b0g2 =',e14.6,'   r0g2 =',e14.6,
     &    /'   a0g3 =',1pe14.6,'     b0g3 =',e14.6,'   r0g3 =',e14.6,
     &    /'  a90g1 =',1pe14.6,'    b90g1 =',e14.6,'  r90g1 =',e14.6,
     &    /'  a90g2 =',1pe14.6,'    b90g2 =',e14.6,'  r90g2 =',e14.6,
     &    /'  a90g3 =',1pe14.6,'    b90g3 =',e14.6,'  r90g3 =',e14.6,
     &    /71('='),/)

      end
