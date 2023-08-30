

      subroutine pes(x,igrad,path,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=4
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: CARTNU(natoms,3), R(6), v
      double precision :: rcm, r1, r2
      logical, save :: first_time_data=.true.

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      CARTNU=0.d0
      do iatom=1,natoms
      do idir=1,3
        CARTNU(iatom,idir)=x(iatom,idir)/0.529177211
      enddo
      enddo

      if(first_time_data) then
      call bjinit(path)
      first_time_data=.false.
      endif


      call CARTTOR(CARTNU,R)
      rcm=R(1)                                                                  
      r1=R(2)                                                                   
      r2=R(3)                                                                   

      call trbunk(rcm,r1,r2,v,4,0,0)

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=v*27.211386
        enddo
      else 
        write (*,*) 'Only energy is available'
      endif

      endsubroutine

      SUBROUTINE CARTTOR(CARTNU,R)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 :: CARTNU(4,3), R(6)
      FLM=18.99840D0
      HYM=1.007825D0
      XCM1=(HYM*CARTNU(1,1)+FLM*CARTNU(2,1))/(FLM+HYM)
      YCM1=(HYM*CARTNU(1,2)+FLM*CARTNU(2,2))/(FLM+HYM)
      ZCM1=(HYM*CARTNU(1,3)+FLM*CARTNU(2,3))/(FLM+HYM)
      XCM2=(HYM*CARTNU(3,1)+FLM*CARTNU(4,1))/(FLM+HYM)
      YCM2=(HYM*CARTNU(3,2)+FLM*CARTNU(4,2))/(FLM+HYM)
      ZCM2=(HYM*CARTNU(3,3)+FLM*CARTNU(4,3))/(FLM+HYM)
      XCM3=XCM2-XCM1
      YCM3=YCM2-YCM1
      ZCM3=ZCM2-ZCM1
      XRM1=CARTNU(1,1)-XCM1
      YRM1=CARTNU(1,2)-YCM1
      ZRM1=CARTNU(1,3)-ZCM1
      THETA1=(XRM1*XCM3+YRM1*YCM3+ZRM1*ZCM3)
      THETA1=THETA1/(SQRT(XRM1**2+YRM1**2+ZRM1**2))
      THETA1=THETA1/(SQRT(XCM3**2+YCM3**2+ZCM3**2))
      IF(THETA1.GT.1.0D0)THETA1=1.0D0
      IF(THETA1.LT.-1.0D0)THETA1=-1.0D0
      THETA1=ACOS(THETA1)
      XRM2=CARTNU(3,1)-XCM2
      YRM2=CARTNU(3,2)-YCM2
      ZRM2=CARTNU(3,3)-ZCM2
      THETA2=(XRM2*(-XCM3)+YRM2*(-YCM3)+ZRM2*(-ZCM3))
      THETA2=THETA2/(SQRT(XRM2**2+YRM2**2+ZRM2**2))
      THETA2=THETA2/(SQRT(XCM3**2+YCM3**2+ZCM3**2))
      IF(THETA2.GT.1.0D0)THETA2=1.0D0
      IF(THETA2.LT.-1.0D0)THETA2=-1.0D0
      THETA2=ACOS(THETA2)
      PI=ACOS(-1.0D0)
      THETA2=PI-THETA2
      Q1=SQRT(XRM1**2+YRM1**2+ZRM1**2)
      Q2=SQRT(XRM2**2+YRM2**2+ZRM2**2)
      CMM=(XCM3**2+YCM3**2+ZCM3**2)
      CMM=SQRT(CMM)
      HHD=(CARTNU(1,1)-CARTNU(3,1))**2 +
     +    (CARTNU(1,2)-CARTNU(3,2))**2 +
     +    (CARTNU(1,3)-CARTNU(3,3))**2
      HHD=SQRT(HHD)
      Q=CMM-Q1*COS(THETA1)+Q2*COS(THETA2)
      Q3=SQRT(ABS(HHD**2-Q**2))
      Q1=Q1*SIN(THETA1)
      Q2=Q2*SIN(THETA2)
      CPHI=(Q1**2+Q2**2-Q3**2)/(2.*Q1*Q2)
      IF(CPHI.LT.-1.0D0)CPHI=-1.0D0
      IF(CPHI.GT.1.0D0)CPHI=1.0D0
      PHI=ACOS(CPHI)
      R(1)=SQRT(XCM3**2+YCM3**2+ZCM3**2)
      R(2)=(SQRT(XRM1**2+YRM1**2+ZRM1**2))*(FLM+HYM)/FLM
      R(3)=(SQRT(XRM2**2+YRM2**2+ZRM2**2))*(FLM+HYM)/FLM
      R(4)=THETA1
      R(5)=THETA2
      R(6)=PHI
      ENDSUBROUTINE CARTTOR 

c bjinit reads in parameters for BJKKL potential from file bjkkl.par
      subroutine bjinit(path)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character(len=1024), intent(in) :: path
      character(len=1024) :: file_path1
      DIMENSION C(359)
      COMMON/PARAMS/C


      file_path1 = trim(path)//"/F2H2/bjkkl.par"
      open(unit=11,file=file_path1)
      READ(11,10)(C(I),I=1,359)
      close(11)
 10   FORMAT(5F14.10)

      return
      end
c
c
c
c  bjinit should be called before the first call to trbunk
c  
      subroutine trbunk(r,r1,r2,v,iq1,iq2,imu)

c     subroutine trbunk transforms the surface of Bunker, Jensen, Kofranek,
c     Lischka and Karpfen(the angular part only; no two-body terms) to the
c     representation of Launay et al.

c     references

c     Bunker potential: P.R. Bunker, P. Jensen, A. Karpfen, M. Kofranek, and
c     H. Lischka, J. Chem. Phys. 92, 7432 (1990).
c     Launay expansion: J.M. Launay, J. Phys. B 10, 3665 (1977).

      implicit double precision (a-h,o-z)
      save
      parameter (lmax=4,lsummax=9,ipara=1)
      parameter (zero=0.0d0,one=1.d0,two=2.d0,four=4.d0)
      parameter (izero=0)
c     dimension a(0:4,0:4,0:8),v(ipara)
c     dimension iq1(ipara),iq2(ipara),imu(ipara)
      dimension a(0:4,0:4,0:8)
      dimension c(359)
      common/acoef/a
      common/PARAMS/c
      common/consts/pi,bo,ha,vl,pc

      bo=0.529177249d0
      ha=219474.63067d0
      pi=3.1415926535898d0
      vl=2.99792458d0
      pc=6.6260755d0

c     next get the a's from subroutine vcalc

      call vcalc(r1,r2,r)

      xl1=float(iq1)
      xl2=float(iq2)
      xm=float(imu)
      lm=min(iq1,iq2)
      ll=abs(iq1-iq2)
      lu=iq1+iq2
      v=0.0d0

      do 20 l=ll,lu,2
        xl=float(l)
        phas=(-1)**(iq1-iq2)
        const1=one
        if(imu.eq.0)const1=two
        const1=one/dsqrt(const1)
        const2=two*xl+one
        cdenom=(four*pi)**3
        const3=dsqrt(two/cdenom)
        const=phas*const1*const2*const3
        v=v+a(iq1,iq2,l)*const*threej(xl1,xl2,xl,xm,-xm,0.d0)
 20   continue

c     divide next two terms by three because there is a bug in Bunker's
c     function subprogram "P" which computes the associated Legendre
c     function. (this is a very small correction).

      if(iq1.eq.4.and.iq2.eq.2.and.imu.eq.2)v=v/3.d0
      if(iq1.eq.2.and.iq2.eq.4.and.imu.eq.2)v=v/3.d0

c     vpot=zero
c     kount=kount-1
c     do 50 j=1,kount
c     write(65,*)v(j)
c50   continue

      end
      FUNCTION THREEJ(A,B,C,AL,BE,GA)                                   HF206430
C                                                                       HF206440
C     CALCULATES WIGNER 3-J SYMBOLS USING FORMULA OF P34 OF BRINK &     HF206450
C     SATCHLER IN "ANGULAT MOMENTUM"                                    HF206460
C     NOTE FAKT(A)=FACT(A)/10**A                                        HF206470
C     A=J1,B=J2,C=J,AL=M1,BE=M2,GA=-M                                   HF206480
C                                                                       HF206490
C                                                                       HF206500
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               HF206510
      SAVE
C                                                                       HF206520
      THREEJ=0.D0                                                       HF206530
C                                                                       HF206540
C     (J1+J2).GE.J AND J.GE.ABS(A-B)    -M=M1+M2    J1,J2,J.GE.0        HF206550
C     ABS(M1).LE.J1    ABS(M2).LE.J2   ABS(M).LE.J                      HF206560
C                                                                       HF206570
      IF(C.GT.A+B)RETURN                                                HF206580
      IF(C.LT.DABS(A-B))RETURN                                          HF206590
      IF(A.LT.0.D0.OR.B.LT.0.D0.OR.C.LT.0.D0)RETURN                     HF206600
      IF(A.LT.DABS(AL).OR.B.LT.DABS(BE).OR.C.LT.DABS(GA))RETURN         HF206610
      IF(-1.D0*GA.NE.AL+BE)RETURN                                       HF206620
C                                                                       HF206630
C                                                                       HF206640
C     COMPUTE DELTA(ABC)                                                HF206650
C                                                                       HF206660
      DELTA=DSQRT(FAKT(A+B-C)*FAKT(A+C-B)*FAKT(B+C-A)/FAKT(A+B+C+1.D0)) HF206670
C                                                                       HF206680
C                                                                       HF206690
      TERM1=FAKT(A+AL)*FAKT(A-AL)                                       HF206700
      TERM2=FAKT(B-BE)*FAKT(B+BE)                                       HF206710
      TERM3=FAKT(C+GA)*FAKT(C-GA)                                       HF206720
      TERM=DSQRT((2.D0*C+1.D0)*TERM1*TERM2*TERM3)                       HF206730
C                                                                       HF206740
C                                                                       HF206750
C     NOW COMPUTE SUMMATION TERM                                        HF206760
C                                                                       HF206770
C     SUM TO GET SUMMATION IN EQ(2.34) OF BRINK AND SATCHLER.  SUM UNTILHF206780
C     A TERM INSIDE FACTORIAL GOES NEGATIVE.  NEW IS INDEX FOR SUMMATIONHF206790
C     .  NOW FIND WHAT THE RANGE OF NEW IS.                             HF206800
C                                                                       HF206810
C                                                                       HF206820
      NEWMIN=IDNINT(DMAX1((A+BE-C),(B-C-AL),0.D0))                      HF206830
      NEWMAX=IDNINT(DMIN1((A-AL),(B+BE),(A+B-C)))                       HF206840
C                                                                       HF206850
C                                                                       HF206860
      SUMM=0.D0                                                         HF206870
C                                                                       HF206880
C                                                                       HF206890
             DO 10 NEW=NEWMIN,NEWMAX                                    HF206900
             DNEW=DFLOAT(NEW)                                           HF206910
             TERM4=FAKT(A-AL-DNEW)*FAKT(C-B+AL+DNEW)                    HF206920
             TERM5=FAKT(B+BE-DNEW)*FAKT(C-A-BE+DNEW)                    HF206930
             TERM6=FAKT(DNEW)*FAKT(A+B-C-DNEW)                          HF206940
             SUMM=SUMM+(-1.D0)**NEW/(TERM4*TERM5*TERM6)                 HF206950
10           CONTINUE                                                   HF206960
C                                                                       HF206970
C     SO CLEBSCH-GORDEN <J1J2M1M2LJM> IS CLEBSH                         HF206980
C                                                                       HF206990
      CLEBSH=DELTA*TERM*SUMM/DSQRT(10.D0)                               HF207000
C                                                                       HF207010
C     CONVERT CLEBSCH-GORDEN TO THREEJ                                  HF207020
C                                                                       HF207030
      IPHASE=IDNINT(A-B-GA)                                             HF207040
      THREEJ=(-1.D0)**(IPHASE)*CLEBSH/DSQRT(2.D0*C+1.D0)                HF207050
C                                                                       HF207060
C                                                                       HF207070
      RETURN                                                            HF207080
      END                                                               HF207090
      SUBROUTINE VCALC(R1,R2,RAB)                                       HF200670
C                                                                       HF200680
C  CALCULATES V AT POINT (R1,R2,RAB,THETA1,THETA2,TAU)                  HF200690
C  USING THE CURRENT VALUES OF THE PARAMETERS.                          HF200700
C    V CALCULATED USING GREEN'S EXPRESSION; SEE JCP 62,2271 (1975).     HF200710
C       G IS THE SCRIPT I USED BY GREEN IN HIS EQ.(A9)                  HF200720
C                                                                       HF200730
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                HF200740
      SAVE
      INTEGER OUT                                                       HF200750
      DIMENSION C(359),G(0:4,0:4,0:8),A(0:4,0:4,0:8),REQ(6)             HF200760
      COMMON/PARAMS/C                                                   HF200770
      COMMON/CONSTS/PI,BO,HA,VL,PC                                      HF200780
      COMMON/ACOEF/A
C                                                                       HF200790
      OUT=6                                                             HF200800
C                                                                       HF200810
c     R1=REQ(1)                                                         HF200820
c     R2=REQ(2)                                                         HF200830
c     RAB=REQ(3)                                                        HF200840
c     T1DEG=REQ(4)                                                      HF200850
c     T2DEG=REQ(5)                                                      HF200860
c     TAUDEG=REQ(6)                                                     HF200870
c     T1=T1DEG*PI/180.D0                                                HF200880
c     T2=T2DEG*PI/180.D0                                                HF200890
c     TAU=TAUDEG*PI/180.D0                                              HF200900
C                                                                       HF200910
C
C     NEXT CALL SUBROUTINE HF2INT WHICH CALCULATES THE RELEVANT ATOM-ATOM
C     DISTANCES NEEDED TO COMPUTE THE TWO-BODY TERMS. THE ANGLES MUST
C     BE IN DEGREES SO CONVERT THEM HERE.                           WCN 10/14/91
      T1DEG=T1*180.D0/PI
      T2DEG=T2*180.D0/PI
      TAUDEG=TAU*180.D0/PI
c     CALL HF2INT(R1,RAB,R2,T1DEG,T2DEG,TAUDEG,RHH,RHF1,RHF2,RFF)       HF200920
C                                                                       HF200930
c     DO 100 L1=0,4                                                     HF200940
c     DO 100 L2=0,4                                                     HF200950
c     LL=IABS(L1-L2)                                                    HF200960
c     LU=L1+L2                                                          HF200970
c     LM=L2                                                             HF200980
c     IF(L1.LT.L2)LM=L1                                                 HF200990
c     DO 100 L=LL,LU,2                                                  HF201000
c     XL=DFLOAT(L)                                                      HF201010
c     ZZ=((2.D0*XL+1.D0)/2.D0)/DSQRT(PI)                                HF201020
c     XL1=DFLOAT(L1)                                                    HF201030
c     XL2=DFLOAT(L2)                                                    HF201040
c     COEFF=((-1)**(L1-L2))*ZZ                                          HF201050
c     G(L1,L2,L)=THREEJ(XL1,XL2,XL,0.D0,0.D0,0.D0)*P(L1,0,T1)*          HF201060
c    X           P(L2,0,T2)                                             HF201070
c     IF(LM.EQ.0)GO TO 100                                              HF201080
c       DO 90 M=1,LM                                                    HF201090
c       XM=DFLOAT(M)                                                    HF201100
c 90    G(L1,L2,L)=G(L1,L2,L)+                                          HF201110
c    X             ((-1)**M)*2.0D0*THREEJ(XL1,XL2,XL,XM,-XM,0.D0)*      HF201120
c    X             P(L1,M,T1)*P(L2,M,T2)*DCOS(XM*TAU)                   HF201130
c100  G(L1,L2,L)=COEFF*G(L1,L2,L)                                       HF201140
C                                                                       HF201150
      DO 20 I=0,4                                                       HF201160
      DO 20 J=0,4                                                       HF201170

      DO 20 K=0,8                                                       HF201180
 20   A(I,J,K)=0.0D0                                                    HF201190
C                                                                       HF201200
      D1=R1-C(1)                                                        HF201210
      D2=R2-C(1)                                                        HF201220
      D3=RAB-C(2)                                                       HF201230
      Y1=1.D0-DEXP(-C(3)*D1)                                            HF201240
      Y2=1.D0-DEXP(-C(3)*D2)                                            HF201250
      Y3=1.D0-DEXP(-D3)                                                 HF201260
      DE1=DEXP(-D3)                                                     HF201270
      DE2=DEXP(-2.0D0*D3)                                               HF201280
      DE3=DEXP(-3.0D0*D3)                                               HF201290
      DE4=DEXP(-4.0D0*D3)                                               HF201300
      R6=RAB**6                                                         HF201310
C                                                                       HF201320
C  DAMPING FACTORS AROUND RIN AND ROUT                                  HF201330
C                                                                       HF201340
      ROUT10=0.5D0*(1.D0-DTANH(0.5D0*(RAB-7.5D0)))                      HF201350
      RIN10=0.5D0*(1.D0-DTANH(2.D0*(RAB-4.0D0)))                        HF201360
      RIN01=0.5D0*(1.D0-DTANH(2.D0*(4.0D0-RAB)))                        HF201370
      RMIDDLE=RIN01*ROUT10                                              HF201380
C                                                                       HF201390
      Y1Y22=Y1**2+Y2**2                                                 HF201400
      Y1Y23=Y1**3+Y2**3                                                 HF201410
      Y1Y24=Y1**4+Y2**4                                                 HF201420
      Y1PY2=Y1+Y2                                                       HF201430
      Y1Y2=Y1*Y2                                                        HF201440
C
C    CD, CQ AND CMQ ARE CONVERSION FACTORS FOR MULTIPOLE MOMENTS        HF201460
C                                                                       HF201470
      CD=(BO**3)*PC*VL*HA/100000.D0                                     HF201480
      CQ=CD*BO*BO                                                       HF201490
      CDQ=CD*BO                                                         HF201500
C                                                                       HF201510
C    C(4) TO C(9), C(91) AND C(111) DESCRIBE THE DISPERSION INTERACTION HF201520
C    C(9) IS IN A(2,2,4).                                               HF201530
C                                                                       HF201540
C    C(41) IS THE HF DIPOLE MOMENT IN DEBYE                             HF201550
C    C(101) IS THE HF QUADRUPOLE MOMENT IN DEBYE.ANGSTROM               HF201560
C    C(131) IS THE HF OCTUPOLE MOMENT IN DEBYE.ANGSTROM**2              HF201570
C                                                                       HF201580
C
c     the zero of energy is chosen as the two monomers infinitely separated
c     at r1=r2=1.744 a.u. thus c(11) has been omitted from a(0,0,0).
C

      A(0,0,0)=C(357)*DE2*RIN10                                         HF201590
     X        +C(20)*Y1Y2*ROUT10                                        HF201610
     X           +(C(12)*DE1+C(13)*DE2+C(14)*DE3+C(15)*DE4)*RMIDDLE     HF201620
     X             +(C(4)/R6 + C(5)/(RAB**7) + C(6)/(RAB**8)            HF201630
     X             +C(7)/(RAB**9)+C(8)/(RAB**10))*RIN01                 HF201640
     X             +(C(157)*(Y1+Y2)*Y3                                  HF201650
     X             +C(158)*Y1*Y2*Y3                                     HF201660
     X             +C(166)*Y1*Y2*Y3**2                                  HF201670
     X             +C(167)*Y1*Y2*Y3**3)*RMIDDLE                         HF201680

      A(1,0,1)=(C(22)*DE1+C(23)*DE2                                     HF201700
     X        +C(24)*DE3+C(25)*DE4                                      HF201710
     X         +C(26)*Y1                                                HF201720
     X         +C(27)*Y2                                                HF201730
     X         +C(21)*Y1*Y3                                             HF201740
     X         +C(156)*Y2*Y3                                            HF201750
     X         +C(28)*Y1**2                                             HF201760
     X         +C(29)*Y2**2                                             HF201770
     X         +C(30)*Y1*Y2                                             HF201780
     X         +C(116)*Y1**3                                            HF201790
     X         +C(117)*Y2**3                                            HF201800
     X         +C(118)*Y1*Y2**2                                         HF201810
     X         +C(119)*Y2*Y1**2                                         HF201820
     X          +C(159)*Y1*Y2*Y3                                        HF201830
     X          +C(169)*Y1*Y2*Y3**2                                     HF201840
     X          +C(179)*Y1*Y2**2*Y3)*RMIDDLE                            HF201850
      A(0,1,1)=-(C(22)*DE1+C(23)*DE2                                    HF201860
     X         +C(24)*DE3+C(25)*DE4                                     HF201870
     X         +C(26)*Y2                                                HF201880
     X         +C(27)*Y1                                                HF201890
     X         +C(21)*Y2*Y3                                             HF201900
     X         +C(156)*Y1*Y3                                            HF201910
     X         +C(28)*Y2**2                                             HF201920
     X         +C(29)*Y1**2                                             HF201930
     X         +C(30)*Y1*Y2                                             HF201940
     X         +C(116)*Y2**3                                            HF201950
     X         +C(117)*Y1**3                                            HF201960
     X         +C(118)*Y2*Y1**2                                         HF201970
     X         +C(119)*Y1*Y2**2                                         HF201980
     X          +C(159)*Y1*Y2*Y3                                        HF201990
     X          +C(169)*Y1*Y2*Y3**2                                     HF202000
     X          +C(179)*Y1**2*Y1*Y3)*RMIDDLE                            HF202010
C                                                                       HF202020
      A(1,1,0)=(C(32)*DE1+C(33)*DE2                                     HF202030
     X        +C(34)*DE3+C(35)*DE4                                      HF202040
     X         +C(36)*Y1Y22                                             HF202050
     X         +C(37)*Y1Y23                                             HF202060
     X         +C(38)*Y1PY2                                             HF202070
     X         +C(39)*Y1PY2*Y1Y2                                        HF202080
     X         +C(40)*Y1Y2)*RMIDDLE                                     HF202090
      A(1,1,2)=(C(42)*DE1+C(43)*DE2                                     HF202100
     X        +C(44)*DE3+C(45)*DE4)*RMIDDLE                             HF202110
     X-((8.D0*PI**1.5D0*DSQRT(2.D0/15.D0))*(C(41)**2)/(CD*RAB**3))*RIN01HF202120
     X         +(C(46)*Y1Y22                                            HF202130
     X         +C(47)*Y1Y23                                             HF202140
     X         +C(48)*Y1PY2                                             HF202150
     X         +C(49)*Y1PY2*Y1Y2                                        HF202160
     X         +C(50)*Y1Y2)*RMIDDLE                                     HF202170
C                                                                       HF202180
      A(2,0,2)=C(356)/R6+(C(52)*DE1+C(53)*DE2                           HF202190
     X        +C(54)*DE3+C(55)*DE4                                      HF202200
     X         +C(56)*Y1                                                HF202210
     X         +C(57)*Y2                                                HF202220
     X         +C(58)*Y1**2                                             HF202230
     X         +C(59)*Y2**2                                             HF202240
     X         +C(60)*Y1*Y2                                             HF202250
     X         +C(126)*Y1**3                                            HF202260
     X         +C(127)*Y2**3                                            HF202270
     X         +C(128)*Y1*Y2**2                                         HF202280
     X         +C(129)*Y2*Y1**2)*RMIDDLE                                HF202290
      A(0,2,2)=C(356)/R6+(C(52)*DE1+C(53)*DE2                           HF202300
     X        +C(54)*DE3+C(55)*DE4                                      HF202310
     X         +C(56)*Y2                                                HF202320
     X         +C(57)*Y1                                                HF202330
     X         +C(58)*Y2**2                                             HF202340
     X         +C(59)*Y1**2                                             HF202350
     X         +C(60)*Y1*Y2                                             HF202360
     X         +C(126)*Y2**3                                            HF202370
     X         +C(127)*Y1**3                                            HF202380
     X         +C(128)*Y2*Y1**2                                         HF202390
     X         +C(129)*Y1*Y2**2)*RMIDDLE                                HF202400
C                                                                       HF202410
      A(2,1,1)=(C(62)*DE1+C(63)*DE2                                     HF202420
     X        +C(64)*DE3+C(65)*DE4                                      HF202430
     X         +C(66)*Y1                                                HF202440
     X         +C(67)*Y2                                                HF202450
     X         +C(68)*Y1**2                                             HF202460
     X         +C(69)*Y2**2                                             HF202470
     X         +C(70)*Y1*Y2                                             HF202480
     X         +C(136)*Y1**3                                            HF202490
     X         +C(137)*Y2**3                                            HF202500
     X         +C(138)*Y1*Y2**2                                         HF202510
     X         +C(139)*Y2*Y1**2)*RMIDDLE                                HF202520
      A(1,2,1)=-(C(62)*DE1+C(63)*DE2                                    HF202530
     X         +C(64)*DE3+C(65)*DE4                                     HF202540
     X         +C(66)*Y2                                                HF202550
     X         +C(67)*Y1                                                HF202560
     X         +C(68)*Y2**2                                             HF202570
     X         +C(69)*Y1**2                                             HF202580
     X         +C(70)*Y1*Y2                                             HF202590
     X         +C(136)*Y2**3                                            HF202600
     X         +C(137)*Y1**3                                            HF202610
     X         +C(138)*Y2*Y1**2                                         HF202620
     X         +C(139)*Y1*Y2**2)*RMIDDLE                                HF202630
      A(2,1,3)=(C(72)*DE1+C(73)*DE2                                     HF202640
     X        +C(74)*DE3+C(75)*DE4)*RMIDDLE                             HF202650
     X    -RIN01*(8.D0*PI**1.5D0/DSQRT(7.D0))*C(41)*C(101)/(CDQ*RAB**4) HF202660
     X         +(C(76)*Y1                                               HF202670
     X         +C(77)*Y2                                                HF202680
     X         +C(78)*Y1**2                                             HF202690
     X         +C(79)*Y2**2                                             HF202700
     X         +C(80)*Y1*Y2                                             HF202710
     X         +C(146)*Y1**3                                            HF202720
     X         +C(147)*Y2**3                                            HF202730
     X         +C(148)*Y1*Y2**2                                         HF202740
     X         +C(149)*Y2*Y1**2)*RMIDDLE                                HF202750
      A(1,2,3)=-(C(72)*DE1+C(73)*DE2                                    HF202760
     X         +C(74)*DE3+C(75)*DE4)*RMIDDLE                            HF202770
     X    +RIN01*(8.D0*PI**1.5D0/DSQRT(7.D0))*C(41)*C(101)/(CDQ*RAB**4) HF202780
     X         -(C(76)*Y2                                               HF202790
     X         +C(77)*Y1                                                HF202800
     X         +C(78)*Y2**2                                             HF202810
     X         +C(79)*Y1**2                                             HF202820
     X         +C(80)*Y1*Y2                                             HF202830
     X         +C(146)*Y2**3                                            HF202840
     X         +C(147)*Y1**3                                            HF202850
     X         +C(148)*Y2*Y1**2                                         HF202860
     X         +C(149)*Y1*Y2**2)*RMIDDLE                                HF202870
C                                                                       HF202880
      A(2,2,0)=(C(82)*DE1+C(83)*DE2                                     HF202890
     X        +C(84)*DE3+C(85)*DE4                                      HF202900
     X        +C(86)*Y1Y22                                              HF202910
     X         +C(87)*Y1Y23                                             HF202920
     X         +C(88)*Y1PY2                                             HF202930
     X         +C(89)*Y1PY2*Y1Y2                                        HF202940
     X         +C(90)*Y1Y2)*RMIDDLE                                     HF202950
      A(2,2,2)=C(91)/R6                                                 HF202960
     X        +(C(92)*DE1+C(93)*DE2                                     HF202970
     X        +C(94)*DE3+C(95)*DE4                                      HF202980
     X         +C(96)*Y1Y22                                             HF202990
     X         +C(97)*Y1Y23                                             HF203000
     X         +C(98)*Y1PY2                                             HF203010
     X         +C(99)*Y1PY2*Y1Y2                                        HF203020
     X         +C(100)*Y1Y2)*RMIDDLE                                    HF203030
      A(2,2,4)=(C(102)*DE1+C(103)*DE2                                   HF203040
     X        +C(104)*DE3+C(105)*DE4)*RMIDDLE                           HF203050
     X +(8.D0*PI**1.5D0/3.D0)*DSQRT(14.D0/5.D0)*RIN01*                  HF203060
     X                                          (C(101)**2)/(CQ*RAB**5)
     X        +(C(106)*Y1Y22                                            HF203070
     X         +C(107)*Y1Y23                                            HF203080
     X         +C(108)*Y1PY2                                            HF203090
     X         +C(109)*Y1PY2*Y1Y2                                       HF203100
     X         +C(110)*Y1Y2)*RMIDDLE                                    HF203110
C                                                                       HF203120
      A(3,0,3)=C(111)/(RAB**7)                                          HF203130
     X        +(C(112)*DE1+C(113)*DE2                                   HF203140
     X        +C(114)*DE3+C(115)*DE4)*RMIDDLE                           HF203150
      A(0,3,3)=-A(3,0,3)                                                HF203160
      A(3,1,2)=(C(122)*DE1+C(123)*DE2                                   HF203170
     X        +C(124)*DE3+C(125)*DE4)*RMIDDLE                           HF203180
      A(1,3,2)=A(3,1,2)                                                 HF203190
      A(3,1,4)=(C(132)*DE1+C(133)*DE2                                   HF203200
     X        +C(134)*DE3+C(135)*DE4)*RMIDDLE                           HF203210
     X   -RIN01*(16.D0*((PI/3.D0)**1.5D0))*C(41)*C(131)/(CQ*RAB**5)     HF203220
      A(1,3,4)=A(3,1,4)                                                 HF203230
      A(3,2,1)=(C(142)*DE1+C(143)*DE2                                   HF203240
     X        +C(144)*DE3+C(145)*DE4)*RMIDDLE                           HF203250
      A(2,3,1)=-A(3,2,1)                                                HF203260
      A(3,2,3)=(C(152)*DE1+C(153)*DE2                                   HF203270
     X        +C(154)*DE3+C(155)*DE4)*RMIDDLE                           HF203280
      A(2,3,3)=-A(3,2,3)                                                HF203290
      A(3,2,5)=(C(162)*DEXP(-0.5D0*D3)+C(163)*DEXP(-0.75D0*D3)          HF203300
     X        +C(164)*DE1+C(165)*DE4)*RMIDDLE                           HF203310
      A(2,3,5)=-A(3,2,5)                                                HF203320
      A(3,3,0)=(C(172)*DE1+C(173)*DE2                                   HF203330
     X        +C(174)*DE3+C(175)*DE4                                    HF203340
     X        +C(176)*Y1Y22+C(177)*Y1Y2)*RMIDDLE                        HF203350
      A(3,3,2)=(C(182)*DE1+C(183)*DE2                                   HF203360
     X        +C(184)*DE3+C(185)*DE4                                    HF203370
     X        +C(186)*Y1Y22+C(187)*Y1Y2)*RMIDDLE                        HF203380
      A(3,3,4)=(C(192)*DE1+C(193)*DE2                                   HF203390
     X        +C(194)*DE3+C(195)*DE4                                    HF203400
     X        +C(196)*Y1Y22+C(197)*Y1Y2)*RMIDDLE                        HF203410
      A(3,3,6)=(C(202)*DE1+C(203)*DE2                                   HF203420
     X        +C(204)*DE3+C(205)*DE4                                    HF203430
     X        +C(206)*Y1Y22+C(207)*Y1Y2)*RMIDDLE                        HF203440
C                                                                       HF203450
      A(4,0,4)=C(212)*DE1+C(213)*DE2                                    HF203460
     X        +C(214)*DE3+C(215)*DE4                                    HF203470
      A(0,4,4)=A(4,0,4)                                                 HF203480
      A(4,1,3)=C(222)*DE1+C(223)*DE2                                    HF203490
     X        +C(224)*DE3+C(225)*DE4                                    HF203500
      A(1,4,3)=-A(4,1,3)                                                HF203510
      A(4,1,5)=C(232)*DEXP(-0.5D0*D3)+C(233)*DEXP(-0.75D0*D3)           HF203520
     X        +C(234)*DE1+C(235)*DE2                                    HF203530
      A(1,4,5)=-A(4,1,5)                                                HF203540
      A(4,2,2)=C(242)*DE1+C(243)*DE2                                    HF203550
     X        +C(244)*DE3+C(245)*DE4                                    HF203560
      A(2,4,2)=A(4,2,2)                                                 HF203570
      A(4,2,4)=C(252)*DE1+C(253)*DE2                                    HF203580
     X        +C(254)*DE3+C(255)*DE4                                    HF203590
      A(2,4,4)=A(4,2,4)                                                 HF203600
      A(4,2,6)=C(262)*DEXP(-0.5D0*D3)+C(263)*DEXP(-0.75D0*D3)           HF203610
     X        +C(264)*DE1+C(265)*DE2                                    HF203620
      A(2,4,6)=A(4,2,6)                                                 HF203630
      A(4,3,1)=C(272)*DE1+C(273)*DE2                                    HF203640
     X        +C(274)*DE3+C(275)*DE4                                    HF203650
      A(3,4,1)=-A(4,3,1)                                                HF203660
      A(4,3,3)=C(282)*DE1+C(283)*DE2                                    HF203670
     X        +C(284)*DE3+C(285)*DE4                                    HF203680
      A(3,4,3)=-A(4,3,3)                                                HF203690
      A(4,3,5)=C(292)*DE1+C(293)*DE2                                    HF203700
     X        +C(294)*DE3+C(295)*DE4                                    HF203710
      A(3,4,5)=-A(4,3,5)                                                HF203720
      A(4,3,7)=C(302)*DE1+C(303)*DE2                                    HF203730
     X        +C(304)*DE3+C(305)*DE4                                    HF203740
      A(3,4,7)=-A(4,3,7)                                                HF203750
      A(4,4,0)=C(312)*DE1+C(313)*DE2                                    HF203760
     X        +C(314)*DE3+C(315)*DE4                                    HF203770
     X        +(C(316)*Y1Y22+C(317)*Y1Y2)                               HF203780
      A(4,4,2)=C(322)*DE1+C(323)*DE2                                    HF203790
     X        +C(324)*DE3+C(325)*DE4                                    HF203800
      A(4,4,4)=C(332)*DE1+C(333)*DE2                                    HF203810
     X        +C(334)*DE3+C(335)*DE4                                    HF203820
      A(4,4,6)=C(342)*DE1+C(343)*DE2                                    HF203830
     X        +C(344)*DE3+C(345)*DE4                                    HF203840
      A(4,4,8)=C(352)*DE1+C(353)*DE2                                    HF203850
     X        +C(354)*DE3+C(355)*DE4                                    HF203860
C                                                                       HF203870
      DO 992 I=0,4                                                      HF203900
      DO 992 J=0,4                                                      HF203910
      IL=IABS(I-J)                                                      HF203920
      IU=I+J                                                            HF203930
      DO 992 K=IL,IU,2                                                  HF203940
      IF(I.NE.4.AND.J.NE.4)GO TO 992                                    HF203950
      A(I,J,K)=A(I,J,K)*RMIDDLE                                         HF203960
 992  CONTINUE                                                          HF203970
C                                                                       HF203980
c     V = 0.0D0                                                         HF203990
c     DO 2 I=0,4                                                        HF204000
c     DO 2 J=0,4                                                        HF204010
c     IL=IABS(I-J)                                                      HF204020
c     IU=I+J                                                            HF204030
c     DO 2 K=IL,IU,2                                                    HF204040
c2    V = V + A(I,J,K)*G(I,J,K)                                         HF204050
C                                                                       HF204060
C  HH, FF, HF1 AND HF2 2-BODY TERMS COMBINED IN TERM VSRG.              HF204070
C  VSRG IS DAMPED TO ZERO IF RHF.GT.3.0                                 HF204080
C  IN ORDER TO PREVENT SPURIOUS 'HOLES' IN V.                           HF204090
C  DAMPED TO ZERO IF RHH.GT.2.3, RHF1/2.GT.2.3, AND RFF.GT.2.3          HF204080
C  SO THAT BOUND STATES CAN BE CALCULATED WITH THE NEGLECT              HF204100
C  OF VSRG AND CLOSE COUPLING PROGRAM 'BOUND' USED WITH V.              HF204100
C  eliminate two-body terms completely in this version wcn 3/10/92
C                                                                       HF204100
c     YRHF=0.25D0*(1.D0-DTANH(7.D0*(R1-3.0D0)))                         HF204110
c    X           *(1.D0-DTANH(7.D0*(R2-3.0D0)))                         HF204120
c     YRHH=0.50D0*(1.D0-DTANH(7.D0*(RHH-2.3D0)))                        HF204130
c     YRFF=0.50D0*(1.D0-DTANH(7.D0*(RFF-2.3D0)))                        HF204130
c     YRHF1=0.50D0*(1.D0-DTANH(7.D0*(RHF1-2.3D0)))                      HF204140
c     YRHF2=0.50D0*(1.D0-DTANH(7.D0*(RHF2-2.3D0)))                      HF204150
c                                                                       HF204160
c     VHH=C(161)*DEXP(-C(171)*RHH)+C(181)*DEXP(-C(191)*RHH)             HF204170
c     VHF1=C(201)*DEXP(-C(211)*RHF1)+C(221)*DEXP(-C(231)*RHF1)          HF204170
c     VHF2=C(201)*DEXP(-C(211)*RHF2)+C(221)*DEXP(-C(231)*RHF2)          HF204170
c     VFF=C(241)*DEXP(-C(251)*RFF)+C(261)*DEXP(-C(271)*RFF)             HF204170
c
c     VSRG=(VHH*YRHH+VHF1*YRHF1+VHF2*YRHF2+VFF*YRFF)*YRHF               HF204200
c
c     V=V+VSRG                                                          HF204210
c                                                                       HF204220
      RETURN                                                            HF204230
      END                                                               HF204240
      FUNCTION FAKT(A)                                                  HF207100
C     THIS FUNCTION COMPUTES N!/(10**N).  STOP 101 IMPLIES TRIED        HF207110
C     TO TAKE FACTORIAL OF A NEGATIVE NO.                               HF207120
C                                                                       HF207130
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               HF207140
      SAVE
C                                                                       HF207150
      AX=A                                                              HF207160
      FAKT=1.D0                                                         HF207170
      IF(AX.EQ.0.D0)RETURN                                              HF207180
      FAKT=.1D0                                                         HF207190
      IF(AX.LT.0.D0)STOP 101                                            HF207200
C                                                                       HF207210
      IC=IDNINT(AX)                                                     HF207220
      AX=AX/10.D0                                                       HF207230
      FAKT=AX                                                           HF207240
             DO 10  I=1,IC-1                                            HF207250
             FAKT=FAKT*(AX-DFLOAT(I)*.1D0)                              HF207260
10           CONTINUE                                                   HF207270
      RETURN                                                            HF207280
      END                                                               HF207290

