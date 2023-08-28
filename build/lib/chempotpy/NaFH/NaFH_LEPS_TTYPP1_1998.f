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
C   System:          NaFH
C   Common name:     N3
C   Functional form: Modified extended LEPS
C   Number of derivatives: 0
C   Number of electronic surfaces: 2
C   Interface: 3-1V

C   Reference:  Unpublished, M. S. Topaler and D. G. Truhlar

C   Notes:  Development version of NaFH potential.
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
C      The classical potential energy is set equal to zero for the Na
C      infinitely far from the HF diatomic and R(HF) set equal to the
C      HF equilibrium diatomic value.
C
C   Parameters:
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Na-H)
C                            R(2) = R(H-F)
C                            R(3) = R(Na-F)

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
      data scale/1.5d0/
      data cprime/4.0d-1/
      data dt2/2.885656d-2/,ret2/1.237624d0/,betat2/1.038036d+1/,
     &     c1/1255969.4d0/,c2/4.7d0/,c3/5.546257d0/,
     &     c4/5.573967d-1/
      data ads1/0.08d0/,ars1/4.0d0/,abs1/1.55d0/,
     &     cds1/2.40d0/,crs1/4.15d0/,cbs1/1.20d0/
      data swa/0.7d0/,swr0/6.5d0/,swa2/1.1d0/,swr20/1.9d0/,
     &     c2a/1.5d0/,c2b/1.0d0/,c2c/0.15d0/,deh2/3.39687d0/
      data cfm/0.0d0/

      rac2=1.d0/dsqrt(2.d0)
      cevau=1.d0/27.2113961d0

      return
c
      entry pot11(r,e,nt)

      do 10 i=1,nt

C***********************************************************************
      r2s=r(i,2)*r(i,2)
      rr1 = r(i,1)*scale
      rr3 = r(i,3)*scale
      rr1s = rr1**2
      rr3s = rr3**2
C***********************************************************************

C***********************************************************************
C*****SWITCHING FUNCTION FOR THE SINGLETS BELOW, AND ITS DERIVATIVES:
      hlp1 = 1.d0-dtanh(swa*(rr1-swr0))
      hlp2 = 1.d0+dtanh(swa2*(r(i,2)-swr20))
      hlp3 = 1.d0-dtanh(swa*(rr3-swr0))
      swt=0.125d0*hlp1*hlp2*hlp3
C***********************************************************************

C***********************************************************************
C*****SINGLETS AND THEIR DERIVATIVES:
C-----------------------------------------------------------------------
C*****MH SINGLET:
      axs1=dexp(-0.5d0*abs1*(rr1-ars1))
      as1=0.5d0*ads1*axs1*(axs1-2.d0)

      cxs1=dexp(-2.80d0*abs1*(rr1-3.500d0))
      cs1=0.35d0*cxs1*(cxs1+2.d0)

      s1=as1+(cs1-as1)*swt
C-----------------------------------------------------------------------
C*****H2 SINGLET, AND ITS DERIVATIVES:
C     ORIGINAL
      s2a=(109.384d0+463.033d0*r(i,2))*dexp(-6.240d0*r(i,2))
     & -r2s*(13.654d0+20.717d0*r(i,2))*dexp(-2.032d0*r(i,2))+1.34357d0

      hlp1 = dexp(-8.d0*(r(i,2)-1.411655d0))
      hlp2 = 0.2d0*dexp(-2.032d0*r(i,2))
      s2h = hlp1 - hlp2

      hlp=dtanh((s2a-s2h)/cprime)
c      hlp1=(s2a-s2h)/(cprime*dcosh((s2a-s2h)/cprime)**2)
      s2a=0.5d0*(s2a+s2h-hlp*(s2a-s2h))

      a = 1.00075d0
      dn = 3.39687d0
      hlp = dexp(-a*(r(i,2)-1.42d0))
      s2i = dn*hlp*(hlp-2.d0)

      hlp = 0.5d0*(r(i,1)+r(i,3)-8.d0)
      hlp = 0.5d0*(1.d0-dtanh(hlp))
      s2 = s2a+(s2i-s2a)*hlp

C-----------------------------------------------------------------------
C*****MF SINGLET:
      axs3=dexp(-1.3d0*abs1*(rr3-ars1))
      as3=ads1*axs3*(axs3-2.d0)

      cxs3=dexp(-cbs1*(rr3-crs1))
      cs3=cds1*cxs3*(cxs3-2.d0)

      s3=as3+(cs3-as3)*swt
C***********************************************************************

C***********************************************************************
C*****TRIPLETS AND THEIR FIRST DERIVATIVES:
C-----------------------------------------------------------------------
C     MH TRIPLET
      hlp1=dexp(-1.10d0*c2*(rr1-0.6d0))
      hlp2=dexp(-1.65d0*c4*(rr1-0.6d0))
      t1=c1*hlp1+c3*hlp2

C-----------------------------------------------------------------------
C     HH TRIPLET
      xt2=dexp(-betat2*(r(i,2)-ret2))
      t2=dt2*xt2*(xt2-2.d0)

C-----------------------------------------------------------------------
C     MF TRIPLET
      hlp1=dexp(-c2*(rr3-0.55d0))
      hlp2=dexp(-2.5d0*c4*(rr3-0.55d0))
      t3=c1*hlp1+c3*hlp2
C***********************************************************************

C***********************************************************************
C-----------------------------------------------------------------------
C*****COULOMB INTEGRALS AND THEIR FIRST DERIVATIVES:

      coul1=0.5d0*(s1+t1)
      coul2=0.5d0*(s2+t2)
      coul3=0.5d0*(s3+t3)

C-----------------------------------------------------------------------
C*****EXCHANGE INTEGRALS AND THEIR FIRST DERIVATIVES:

      exch1=0.5d0*(s1-t1)
      exch2=0.5d0*(s2-t2)
      exch3=0.5d0*(s3-t3)

C***********************************************************************

C***********************************************************************
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2

      cplg2=c2a*dexp(-c2b*w-c2c*(rr1+r(i,2)+rr3))

      rootw = rac2*dsqrt(w+cplg2*cplg2)

      e(i)=((coul1+coul2+coul3+deh2-rootw)*cevau-cfm)

C***********************************************************************

10    continue

      return

100   format(/24('='),' MHF - U11 potential ',24('='),
     &  /'   dt2 =',1pe14.6,2x,'  ret2 =',e14.6,2x,'betat2 =',e14.6,
     &  /' c1 =',e14.6,2x,' c2 =',e14.6,2x,' c3 =',e14.6,
     &  /' c4 =',e14.6,2x,
     &  /'  ads1 =',e14.6,
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

c      BLOCK DATA UONETWO
      double precision ahf,anaf,anah,r21,r23,r13
      double precision ccMF(4),ccMH(4)
c      common/ucoupl/ahf,anaf,anah,r21,r23,r13
c      common/ucoupl1/ccMF(4),ccMH(4)

      save

      data ccMF/3.9960D0,5.4730D-01,-1.1428D0,2.0440D-01/
      data ccMH/1.00d0,0.80d0,-2.67d0,0.456d0/

      data ahf/3.0000d0/,
     &    anaf/18.0000d-1/,
     &    anah/16.0000d-1/
      data r21/0.7000d0/,
     &     r23/-3.3006d-1/,
     &     r13/1.2431d0/
c      END

c      save /ucoupl/,/ucoupl1/

c      common/ucoupl/ahf,anaf,anah,r21,r23,r13
c      common/ucoupl1/ccMF(4),ccMH(4)

      cevau=1.d0/27.2113961d0

      return
c
      entry pot12(r,e,nt)

      do 10 i=1,nt


      hlp = ccMF(3)/(1.d0+ccMF(4)*r(i,3))
      unaf = dexp(-ccMF(1)-ccMF(2)*r(i,3)- hlp*r(i,3))
      unaf = unaf*27.2113961d0

      hlp = ccMH(3)/(1.d0+ccMH(4)*r(i,1))
      unah = dexp(-ccMH(1)-ccMH(2)*r(i,1)- hlp*r(i,1))
C-----------------------------------------------------------------------

      hlp = anah*r(i,1)-anaf*r(i,3)-r13
      sw1 = 0.5d0*(1.d0+dtanh(hlp))

      hlp = ahf*r(i,2)-anah*r(i,1)-r21
      sw21 = 0.5d0*(1.d0+dtanh(hlp))

      hlp = ahf*r(i,2)-anaf*r(i,3)-r23
      sw23 = 0.5d0*(1.d0+dtanh(hlp))

      e(i) = unaf*sw1*sw23+unah*(1.d0-sw1)*sw21

      e(i) = e(i)*cevau

C***********************************************************************
C     MODIFICATION: coupling cut-off at large H-H separations
      gmc = 0.3d0
      amc = 0.40d0
      rmc = 2.7d0
      hlp= gmc*(1.d0-dtanh(amc*(r(i,2)-rmc)))
      e(i) = e(i)*hlp
C***********************************************************************

10    continue

      return

  100 format(/25('='),' MFH - U12 coupling ',25('='))
  101 format(2x,'ccMF = ',5x,4e14.6/10x,2e14.6)
  102 format(2x,'ccMH = ',4x,4e14.6/11x,2e14.6)
  106 format(2x,'ahf = ',7x,3e14.6)
  107 format(2x,'anaf = ',6x,3e14.6)
  108 format(2x,'anah = ',6x,3e14.6)
  109 format(2x,'r21 = ',7x,3e14.6)
  110 format(2x,'r23 = ',7x,3e14.6)
  111 format(2x,'r13 = ',7x,3e14.6
     &    /71('=')/)

      end

c
c  U22 surface - upper diabatic surface
c
      subroutine prepot22
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt)

      save
      data scale/1.5d0/
      data c1s1/203.128d0/,c2s1/1.6167d0/,c3s1/-100.124d0/,
     &     c4s1/68.078d0/,c5s1/-15.832d0/,c6s1/1.303d0/
      data dnaht1/3.35000d0/,dnaht2/10.000d0/,dnaht3/54.6393d0/,
     &     dnaht4/-17.63343d0/,dnaht5/3.0085d0/,dnaht6/1.450d0/,
     &     dnaht7/-2.87977d0/,dnaht61/1.300d0/,dnaht11/2.85000d0/
      data deh2/4.7404422d0/,una2p/2.1037d0/,c2a/0.8d0/,c2b/0.8d0/,
     &     c2c/0.15d0/,cprime/8.0d-1/
      data y0d1/2.417800d-1/,y0r1/3.708055d0/,y0b1/6.960127d-1/,
     &     y90d1/1.372065d-2/,y90r1/7.233552d0/,y90b1/5.127734d-1/,
     &     y0d2/8.870172d-1/,y0r2/2.300265d0/,y0b2/8.497374d-1/,
     &     y90d2/1.106187d0/,y90r2/2.002516d0/,y90b2/9.514887d-1/,
     &     swa/0.78495d0/,swr0/7.47654d0/
      data a0g1/-8.210245d-2/,b0g1/4.658698d0/,r0g1/4.552860d0/,
     &     a0g2/3.546427d-2/,b0g2/6.456787d-1/,r0g2/6.332225d0/,
     &     a0g3/-1.562327d-2/,b0g3/1.170537d-1/,r0g3/1.097310d+1/
      data a90g1/-9.929297d-2/,b90g1/4.367220d0/,r90g1/4.402447d0/,
     &     a90g2/4.696181d-2/,b90g2/5.470459d-1/,r90g2/6.289526d0/,
     &     a90g3/-1.219900d-2/,b90g3/1.002291d-2/,r90g3/1.118725d+1/
      data h2t1/1.5d0/,h2t2/4.47009d0/,h2t3/4.48534d0/,
     &     h2t4/1.77908d0/,h2t5/-6.46276d0/,h2t6/3.15396d0/,
     &     h2t7/-1.59912d0/

      rac2=1.d0/dsqrt(2.d0)
      cevau=1.d0/27.2113961d0

      return
c
      entry pot22(r,e,nt)

      do 10 i=1,nt

C***********************************************************************
      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)
      rr1 = r(i,1)*scale
      rr3 = r(i,3)*scale
      rr1s = rr1**2
      rr3s = rr3**2
C***********************************************************************
C*****SWITCHING FUNCTION FOR THE H2 TRIPLET, AND ITS DERIVATIVES:

      swt=0.5d0*(1.d0-dtanh( swa*(0.5d0*(rr1+rr3)-swr0)))
C***********************************************************************

C***********************************************************************
C*****SINGLETS AND THEIR DERIVATIVES:

C-----------------------------------------------------------------------
C*****MH SINGLET:
      raaa = 2.700d0
      aaa = 2.40d0
      aaa1 = 2.60d0
      daa = 1.80d0
      hlpNAH1 = dexp(-aaa1*(r(i,1)-raaa))
      hlpNAH11 = dexp(-1.65d0*aaa*(r(i,1)-raaa))
      s1 = 0.30d0*daa*(hlpNAH11-2.d0*hlpNAH1)
C-----------------------------------------------------------------------
C*****H2 SINGLET
      s2a=(109.384d0+463.033d0*r(i,2))*dexp(-6.240d0*r(i,2))
     &   -r2s*(13.654d0+20.717d0*r(i,2))*dexp(-2.032d0*r(i,2))+una2p

      b = 1.4d0
      a = 2.d0*b
      r00 = 1.411655d0
      dn = -4.74044d0/(1.d0-a/b)
      hlp1 = dexp(-a*(r(i,2)-r00))
      hlp2 = dexp(-b*(r(i,2)-r00))
      s2i = dn*(hlp1-(a/b)*hlp2)+una2p

      hlp = 0.5d0*(r(i,1)+r(i,3)-7.0d0)
      hlp = 0.5d0*(1.d0-dtanh(hlp))
      s2a = s2a+(s2i-s2a)*hlp

C-----------------------------------------------------------------------
C*****MF SINGLET:
      hlpNAH3 = dexp(-aaa*(r(i,3)-raaa))
      hlpNAH33 = dexp(-1.65d0*aaa*(r(i,3)-raaa))
      s3 = daa*(hlpNAH33-2.d0*hlpNAH3)
C***********************************************************************
C*****TRIPLETS AND THEIR DERIVATIVES:

C-----------------------------------------------------------------------
C*****H2 TRIPLET (interaction region):

      t2a=0.5d0*(dexp(-h2t2*(r(i,2)-h2t1))+(h2t3+h2t5*r(i,2)
     &+h2t6*r2s+h2t7*r2s*r(i,2))*dexp(-h2t4*r(i,2)))

C*****H2 TRIPLET (asymptotic)
      s2b=147.214d0*dexp(-3.915d0*r(i,2))+r(i,2)*(41.275d0+10.505d0
     &   *r(i,2)+4.408d0*r2s)*dexp(-2.032d0*r(i,2))

C*****H2 TRIPLET: TOTAL (switched between the asymptotic and the ineraction one)
      t2=swt*(t2a-s2b)+s2b
C-----------------------------------------------------------------------
C*****TRIPLETS AND THEIR DERIVATIVES:
C!!!!!MH TRIPLET:
      t1=3.d0*(dexp(-dnaht2*(rr1-dnaht11))+
     & (dnaht7+dnaht3*rr1+dnaht4*rr1s+dnaht5*rr1*rr1s)
     &  *dexp(-dnaht61*rr1))
C-----------------------------------------------------------------------
C!!!!!MF TRIPLET: (should be the same as the first one, t1)

      t3=3.d0*(dexp(-dnaht2*(rr3-dnaht1))+
     & (dnaht7+dnaht3*rr3+dnaht4*rr3s+dnaht5*rr3*rr3s)
     &  *dexp(-dnaht6*rr3))
C-----------------------------------------------------------------------
C***********************************************************************
C*****TOTAL H2 SINGLET:
C!!!!!H2 SINGLET: total (switched between the true singlet and a triplet
C                        in the asymptotic region at large HH distances)
      hlp=dtanh((s2a-t2)/cprime)
      s2=0.5d0*(s2a+t2-hlp*(s2a-t2))
C***********************************************************************
C***********************************************************************
C-----------------------------------------------------------------------
C*****COULOMB INTEGRALS AND THEIR FIRST DERIVATIVES:

      coul1=0.5d0*(s1+t1)
      coul2=0.5d0*(s2+t2)
      coul3=0.5d0*(s3+t3)
C-----------------------------------------------------------------------
C*****EXCHANGE INTEGRALS AND THEIR FIRST DERIVATIVES:

      exch1=0.5d0*(s1-t1)
      exch2=0.5d0*(s2-t2)
      exch3=0.5d0*(s3-t3)
C***********************************************************************

      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      cplg2=c2a*dexp(-c2b*w-c2c*(rr1+r(i,2)+rr3))
      rootw = rac2*dsqrt(w+cplg2*cplg2)
      e(i)=coul1+coul2+coul3+deh2-rootw

c***********************************************************************
      rf1 = 3.00d0
      af1 = 0.60d0
      fff = 1.d0 - 0.5d0*(1.d0+dtanh(af1*(r(i,3)-rf1)))
      e(i)=((e(i)-1.35d0)/(1.00d0+1.80d0*fff))*cevau
C***********************************************************************
C***********************************************************************
C***********************************************************************

10    continue

      return

100   format(/24('='),' MHF - U22 potential ',25('='),
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
     &    /71('='))

      end


      BLOCK DATA UONETWO
      double precision ahf,anaf,anah,r21,r23,r13
      double precision ccMF,ccMH
      common/ucoupl/ahf,anaf,anah,r21,r23,r13
      common/ucoupl1/ccMF(4),ccMH(4)
      data ccMF/3.9960D0,5.4730D-01,-1.1428D0,2.0440D-01/
      data ccMH/1.00d0,0.80d0,-2.67d0,0.456d0/

      data ahf/3.0000d0/,
     &    anaf/18.0000d-1/,
     &    anah/16.0000d-1/
      data r21/0.7000d0/,
     &     r23/-3.3006d-1/,
     &     r13/1.2431d0/
      END


