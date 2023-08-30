  module shell
  implicit none
  real,parameter::total=8
  integer,parameter::coff=1561
  !real,parameter::alpha(4)=(/2.0,2.0,2.0,2.0/)
  real*8::cof(coff),power(coff,6)
  end module shell

  subroutine pes_init(path)
   use shell
   character(len=1024), intent(in) :: path
   character(len=1024) :: file_path1, file_path2
   integer::i

   file_path1 = trim(path)//"/H2CO/coeff8"
   file_path2 = trim(path)//"/H2CO/index8"
   open(112,status='old',file=file_path1)
   open(113,status='old',file=file_path2)

   read(112,*)
   do i=1,coff
     read(112,*) cof(i)
   end do
   do i=1,coff
     read(113,*) power(i,1),power(i,2),power(i,3),power(i,4),&
       &power(i,5),power(i,6)
   end do
   close(113)
   close(112)

   return
 end subroutine pes_init

 subroutine getpot(cood,v,path)
   use shell
   implicit none
   character(len=1024), intent(in) :: path
   integer::i
   real*8::v,f
   real*8::cood(3,4),com(3)  
   real*8::y(6),r(6)
   real*8::bas(coff)
   real*8 smallCH, smallOH
   real*8 S, x,z,x3,x2,z1,z2
!       This probram expects cood to be in the order C,O,H,H
    
   call edis(r(1),y(1),cood(:,3),cood(:,4))
   call edis(r(2),y(2),cood(:,2),cood(:,3))
   call edis(r(3),y(3),cood(:,1),cood(:,3))
   call edis(r(4),y(4),cood(:,1),cood(:,4))
   call edis(r(5),y(5),cood(:,2),cood(:,4))
   call edis(r(6),y(6),cood(:,1),cood(:,2))


   call basis(bas,y)

    if (r(3)< r(4)) then
       smallCH=r(3)     
       smallOH=r(2)
   else
       smallCH=r(4)
       smallOH=r(5)
    end if

! Switching function has been used here |
       x=r(1)
       z=x-8
       z1=10-8
       z2=z/z1
      x2=z2*z2
      x3=x2*z2

     if (x.ge. 10 ) then
         s=1.0
     else if (x.le. 8 ) then
         s=0.0d0
     else
         s=10*x3-15*x3*z2+6*x3*x2
     end if

     f=0
     do i=1,coff
        f=f+cof(i)*bas(i)
     end do
     f=f+114.332958863-1.059892782251382E-004

    if (x .ge. 8) then
       call hcopot(V,smallCH,r(6),smallOH,path)
       f= f*(1-s)+V*s
     end if 
 
    v=f
 end subroutine getpot

 subroutine edis(dis,edist,atomi,atomj)
   use shell
   real*8::atomi(3),atomj(3)
   real*8::dis,alfa
   real*8,intent(out)::edist
   alfa=2.0
   dis=sqrt((atomi(1)-atomj(1))**2+(atomi(2)-atomj(2))**2+(atomi(3)-atomj(3))**2)
   edist=exp(-dis/alfa)
   return
 end subroutine edis

 subroutine basis(bas,y)
   use shell
   integer m,n,p,q,i,j,k
   real*8::bas(coff)
   real*8::y(6)
   do k=1,coff
     bas(k)=(y(1)**power(k,1))*(y(6)**power(k,6))*((y(2)**power(k,2))*&
       &(y(3)**power(k,3))*(y(4)**power(k,4))*(y(5)**power(k,5))+(y(2)**power(k,5))*&
       &(y(3)**power(k,4))*(y(4)**power(k,3))*(y(5)**power(k,2)))
   end do
   return
 end subroutine

      subroutine hcopot(f,RHC,RCO,ROH,path)
      IMPLICIT NONE
      character(len=1024), intent(in) :: path
      integer j,loop,iv
      real*8 rmin(3),rmax(3),z1,z2,z3,m(3),h(3)
      real*8 pot,potmin,ma,mb,mc,au2ev,a0,pi
      parameter (au2ev = 27.2116d0)
      parameter (a0 = 0.52918d0)
      parameter (pi = 3.1415926536d0)
      real*8 pot3,c3,POTHCO,f
      REAL*8 V, RHC,RCO, CGAMMA, ROH
      integer istate


      ma =  1.0    !  H
      mb = 12.0    !  C
      mc = 16.0    !  O

      istate = 1      !  X-Flaeche
      call pot3in(path)
      CGAMMA=(RHC*RHC+RCO*RCO-ROH*ROH)/(2*RHC*RCO)

      POT=0.0d0
      POT= POTHCO(RHC,RCO,CGAMMA,ISTATE)
      pot=pot+.18140000
      f=pot
      end subroutine hcopot

       SUBROUTINE POT3IN(path)
       IMPLICIT NONE
       character(len=1024), intent(in) :: path
       character(len=1024) :: file_path1
       REAL*8 V, R/1D0/, RKL/1D0/, G/0D0/, POTHCO
       EXTERNAL POTHCO

       file_path1 = trim(path)//"/H2CO/pothco_hmk.para"
       OPEN(10,FILE=file_path1,status='old')
       V= POTHCO(0D0,0D0,0D0,0)
       CLOSE(10)
       RETURN
       END

      FUNCTION POTHCO(RHC,RCO,RCOS, ISTATE)
      IMPLICIT NONE

      INTEGER ISTATE                    ! Initialiserung falls <>0, INPUT
      REAL*8 :: RHC, RCO, RCOS, POTHCO
      ! INPUT: R_HC in a.u., R_CO in a.u.,  COS(gamma_HCO),
      ! OUTPUT: Potentialwert in a.u.
      REAL*8 AU2ANG, AU2EV, EV2CM, KC2EV, AU2CM
      PARAMETER(AU2ANG=0.529177D0)      ! a_0 / Angstroem
      PARAMETER(AU2EV=27.2116D0)        ! hartree / eV
      PARAMETER(EV2CM=8065.44D0)        ! eV / cm^-1
      PARAMETER(KC2EV=23.061D0)         ! {kcal/mol} / eV
      PARAMETER(AU2CM=AU2EV*EV2CM)      ! hartree / cm^-1

      INTEGER, PARAMETER :: NHC=12, NCO=9, NWIN=12

      INTEGER GZHC, GZCO, GZWIN, MLIN(2)   ! GRAD IM ZAEHLER
      INTEGER GNHC, GNCO, GNWIN, MTAN(2)   ! GRAD IM NENNER
      INTEGER IXHC, IXCO, IXOH          ! FESTE VORFAKTOREN IM NENNER.
      INTEGER GWHC, GWCO, GWWIN         ! POLYNOME VOR WURZEL (KONISCHE D.)
      INTEGER GZ2CO, GN2CO, IX2CO       ! DIATOM-POTENTIAL v. C-O.
      INTEGER GZ2HC, GN2HC, IX2HC       ! DIATOM-POTENTIAL v. H-C.
      INTEGER GZ2OH, GN2OH, IX2OH       ! DIATOM-POTENTIAL v. O-H.
      INTEGER MAXLIN, MAXTAN
      PARAMETER(MAXLIN=500, MAXTAN=200) ! Maximalwerte fuer MLIN und MTAN.
      REAL*8 WURADD, NENADD, EASYM
      REAL*8, PARAMETER :: NUHC(12)=(/2.3D0,2.5D0, 2.1D0,3.0D0, 1.9D0,3.5D0,&
               &1.7D0,4.0D0, 1.5D0,5.0D0, 1.3D0,7.0D0/)
      REAL*8, PARAMETER :: NUCO(9)=(/2.2d0,1.9d0, 2.5d0,1.7d0, 3.0d0,2.05d0,&
               &2.35d0,3.5d0, 4.0d0/)
      REAL*8, PARAMETER :: NUWIN(12)=(/100d0,120d0, 80d0,140d0, 60d0,160d0,&    ! NEUE NULLST.
               &40d0,170d0, 20d0,180d0, 10d0,0d0/)
      REAL*8 :: XCOS(NWIN), NUCOS(NWIN), KK(MAXLIN,2),AA(MAXTAN,2),V2CH, V2OH
      EXTERNAL V2CH, V2OH


      INTEGER I, J,K,L, LL,JJ, J1,J2,J3,KS,I1,I2
      LOGICAL QINIT(2)/.false.,.false./ ! zeigt an ob Parameter initialisiert.
      REAL*8 :: PQ,  A, B, C, NENN, RHCSIN,RHCCOS, KONUS
      REAL*8 :: PHC, PCO, PCOS, POH, RHCI, RCOI, RCOSI, ROHI
      REAL*8 :: THY, ACO, AHC, VPLUS, PI
      COMPLEX*16 :: QQ1, QQ2, QQ3, QQ4, WURZ

      IF(ISTATE.LE.0) THEN

        IF(ISTATE.EQ.0)THEN
          I1= 1                         ! 1. Flaeche in Datei (Unit 10).
          I2= 2                         ! letzte Flaeche     "           .
        ELSE IF(ISTATE.EQ.-1) THEN
          I1= 1                         ! 1. Flaeche in Datei (Unit 10).
          I2= 1                         ! letzte Flaeche        "        .
        ELSE IF(ISTATE.EQ.-2) THEN
          I1= 2                         ! 1. Flaeche in Datei (Unit 10).
          I2= 2                         ! letzte Flaeche        "        .
        END IF

        PI= ACOS(-1D0)

        DO J1=1,NHC
        END DO
        DO J2=1,NCO
        END DO
        DO J3=1,NWIN
          NUCOS(J3)= COS(PI/180D0*NUWIN(J3)) ! Nullstellen d. 1d-Polynome.
        END DO

        READ(10,*)
        READ(10,*)
        READ(10,*)   GZHC, GZCO, GZWIN
        READ(10,*)
        READ(10,*)   GNHC, GNCO, GNWIN
        READ(10,*)
        READ(10,*)   IXHC, IXCO, IXOH
        READ(10,*)
        READ(10,*)   GWHC, GWCO, GWWIN
        READ(10,*)
        READ(10,*)   GZ2HC,  GN2HC,  IX2HC
        READ(10,*)
        READ(10,*)   GZ2CO,  GN2CO,  IX2CO
        READ(10,*)
        READ(10,*)   GZ2OH,  GN2OH,  IX2OH
        READ(10,*)
        READ(10,*)   WURADD, NENADD, EASYM

        DO KS=I1,I2

          IF(QINIT(KS)) THEN             ! Test ob schon mal initialisiert.
          ELSE
            QINIT(KS)=.true.
          END IF
          MLIN(KS)=(GZHC+1)*(GZCO+1)*(GZWIN+1)+(GZ2CO+1)
          MTAN(KS)= 2*(GNHC+1)*(GNCO+1)*(GNWIN+1)-2+2*GN2CO 
          IF(KS.EQ.1) THEN               ! Mit Konus:
            MLIN(KS)= MLIN(KS) + (GWHC+1)*(GWCO+1)*(GWWIN+1)
            MTAN(KS)= MTAN(KS) + 4
          END IF
          READ(10,*)
          DO I=1,MLIN(KS)
            READ(10,*) J, KK(I,KS)
            IF(I.NE.J) STOP ' POTHCO: INDEX FALSCH '
          END DO
          READ(10,*)
          DO I=1,MTAN(KS)
            READ(10,*) J, AA(I,KS)
            IF(I.NE.J) STOP ' POTHCO: INDEX FALSCH '
          END DO

        END DO

        POTHCO= 0d0
      ELSE
        RHCI= RHC                                      ! interne Kopien.
        RCOI= RCO
        RCOSI= RCOS
        VPLUS= 0D0                                     ! Potential-Offset.

        THY= TANH(1.1d0*(RHCI-4.4d0))                  ! innen -1, aussen 1.
        A= 0.75d0 - 0.25d0 * THY                       ! innen 1, aussen 1/2.
        AHC= 1d0 + (-0.002d0 + 0.014d0 * RCOSI ) * A   ! Skalierung in R_HC.
        ACO= 1d0 +         0.012d0           * A       ! Skalierung in R_CO.
        RHCI= AHC * (RHCI-2.11d0) + 2.11d0
        RCOI= ACO * (RCOI-2.23d0) + 2.23d0
        VPLUS= ( 0.5d0 * THY - 0.5d0 ) * 50d-3 / 27.2116d0 ! innen -50 meV.
        RHCI= MAX(RHCI,1.3d0)            ! R_HC > 1.3 a.u.
        RCOI= MAX(RCOI,1.5d0)            ! R_CO > 1.5 a.u.
        ROHI= SQRT( RHCI**2 + RCOI**2 - 2d0*RHCI*RCOI*RCOSI )
        ROHI= MAX(ROHI,0.01d0)          ! zur numerischen Stabilitaet.

        KS= ISTATE
        JJ= 1                           ! Index in AA(JJ,KS) (Nenner).

        QQ2=(1d0,0d0)                   ! *** QQ_CO(R_CO)
        PCO= 1d0
        DO J=1,GN2CO
          PCO= PCO * ( RCOI - NUCO(J) ) ! P_{CO,J}(R_CO)
          QQ2= QQ2 + DCMPLX( AA(JJ,KS)*PCO, AA(JJ+1,KS)*PCO )
          JJ= JJ + 2
        END DO

        QQ4=(1d0,0d0)                   ! *** QQ_HCO(R_HC,R_CO,COS())
        PHC= 1d0
        DO J=0,GNHC
          PCO= PHC
          DO K=0,GNCO
            PCOS= PCO
            IF(J.GT.0.OR.K.GT.0) THEN
              QQ4= QQ4 + DCMPLX( AA(JJ,KS)*PCO, AA(JJ+1,KS)*PCO )
              JJ= JJ + 2
            END IF
            DO L=1,GNWIN
              PCOS= PCOS * ( RCOSI - NUCOS(L) )
              QQ4= QQ4 + DCMPLX( AA(JJ,KS)*PCOS, AA(JJ+1,KS)*PCOS )

              JJ= JJ + 2
            END DO
            PCO= PCO * ( RCOI - NUCO(K+1) )
          END DO
          PHC= PHC * ( RHCI - NUHC(J+1) )
        END DO

        IF(KS.EQ.1) THEN
          RHCCOS= RHCI * RCOSI              ! R_HC * cos(gamma)
          RHCSIN= RHCI * SQRT(1d0-RCOSI**2) ! R_HC * sin(gamma)
          A= 1D0 + AA(JJ+1,KS)*RHCCOS&                 ! 1 + a2*RHCCOS +
            &+( AA(JJ+2,KS) + AA(JJ+3,KS)*RCOI ) * RCOI ! + a3*RCOI+a4*RCOI**2
          WURZ= SQRT( A**2 + (AA(JJ,KS)*RHCSIN)**2 + WURADD ) ! WURADD ist
          JJ= JJ + 4
        END IF


        IF(JJ.NE.MTAN(KS)+1) STOP 'POTHCO: MTAN FALSCH !'


        LL= 1

        PQ= EASYM
        PQ= PQ + V2CH(RHCI)


        PCO= RCOI**IX2CO / ABS(QQ2)**2  ! *** PP_CO(R_CO)
        DO J=0,GZ2CO
          PQ= PQ + KK(LL,KS) * PCO
          LL= LL + 1
          PCO= PCO * ( RCOI - NUCO(J+1) ) ! P_{CO,J}(R_CO)
        END DO
        PQ= PQ + V2OH(ROHI)


        NENN= RHCI**(IXHC-IXOH) * RCOI**(IXCO-IXOH) * ROHI**IXOH&
              &/ ( ABS(QQ4)**2 + NENADD ) ! *** PP_HCO(R_HC,R_CO,COS())
        KONUS= 0d0                      ! Konische Durchschneidung:
        IF(KS.EQ.1) THEN                ! (nur auf 1. Flaeche)
          PHC= WURZ * NENN
          DO J=0,GWHC
            PCO= PHC
            DO K=0,GWCO
              PCOS= PCO
              DO L=0,GWWIN
                KONUS= KONUS + KK(LL,KS) * PCOS
                LL= LL + 1
                PCOS= PCOS * ( RCOSI - NUCOS(L+1) )
              END DO
              PCO= PCO * ( RCOI - NUCO(K+1) )
            END DO
            PHC= PHC * ( RHCI - NUHC(J+1) )
          END DO
          PQ= PQ + KONUS
        END IF

        PHC= NENN
        DO J=0,GZHC                     ! normale Polynome:
          PCO= PHC
          DO K=0,GZCO
            PCOS= PCO
            DO L=0,GZWIN
              PQ= PQ + KK(LL,KS) * PCOS
              LL= LL + 1
              PCOS= PCOS * ( RCOSI - NUCOS(L+1) )
            END DO
            PCO= PCO * ( RCOI - NUCO(K+1) )
          END DO
          PHC= PHC * ( RHCI - NUHC(J+1) )
        END DO

        IF(LL.NE.MLIN(KS)+1) STOP ' POTHCO: MLIN FALSCH '


        PQ= PQ + 0.055d0/AU2EV + VPLUS         ! OFFSET(Dez.1994)

        IF (ROHI.LT.1.1d0.OR.RHCI.LT.1.3d0.OR.&
          &RCOI.LT.1.5d0.OR.RCOI+ROHI.LT.2.9d0) THEN
          PQ= MAX( 6d0/AU2EV, PQ )      ! >= 6 eV.
        END IF
        POTHCO= MIN( 20d0/AU2EV, PQ )   ! nach oben begrenzen (20eV)

      END IF

 999  RETURN
      END
      FUNCTION V2OH(R)
      IMPLICIT NONE
      REAL*8 :: V2OH,  R

      REAL*8 AU2ANG, AU2EV
      PARAMETER(AU2ANG=0.529177D0)      ! a_0 / ANGSTROEM
      PARAMETER(AU2EV=27.2116D0)        ! HARTREE / EV

      REAL*8 A1, A2, A3, DE, RE, RHO

      A1= 4.507                    ! Huxley etc. nach Daten aus Herzberg.
      A2= 4.884
      A3= 3.795
      DE= 4.621 / AU2EV
      RE= 0.9696

      RHO= ( R*AU2ANG - RE )

      V2OH=-DE*(1+A1*RHO+A2*RHO**2+A3*RHO**3)*EXP(-A1*RHO)

      RETURN
      END
      FUNCTION V2CH(R)
      IMPLICIT NONE
      REAL*8 :: V2CH, R

      REAL*8 AU2ANG, AU2EV
      PARAMETER(AU2ANG=0.529177D0)      ! a_0 / ANGSTROEM
      PARAMETER(AU2EV=27.2116D0)        ! HARTREE / EV

      REAL*8 A1, A2, A3, DE, RE, RHO

      A1= 3.836
      A2= 3.511
      A3= 2.268
      DE= 3.631 / AU2EV
      RE= 1.1199

      RHO= ( R*AU2ANG - RE )

      V2CH=-DE*(1+A1*RHO+A2*RHO**2+A3*RHO**3 )*EXP(-A1*RHO)

      RETURN
      END

      subroutine pes(x,igrad,path,p,gradient_out,d)
      use shell
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=4
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates)
      double precision, intent(out) :: gradient_out(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      real*8 :: tx(3,natoms), v
      integer :: iatom, idir, j, istate
      logical, save :: first_time_data=.true.
      !initialize 
      gradient_out=0.d0
      d=0.d0

      do iatom=1,natoms
      do idir=1,3
        tx(idir,iatom)=x(iatom,idir)/0.529177211
      enddo 
      enddo

      if (igrad==0) then
        if(first_time_data) then
          call pes_init(path)
          first_time_data=.false.
        endif
        call getpot(tx,v,path)
        v=v*27.211386
        do istate=1,nstates
          p(istate)=v
        enddo
      else
        write (*,*) 'Only energy is available'
      endif

      endsubroutine
