      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad

      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: r(1,3), e(1), v
      integer :: iatom, idir, i, j
      logical, save :: first_time_data=.true.

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      ! input cartesian is HHH
      r(1,1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2&
     &          +(x(1,3)-x(2,3))**2)/0.529177211
      r(1,2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2&
     &          +(x(2,3)-x(3,3))**2)/0.529177211
      r(1,3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2&
     &          +(x(1,3)-x(3,3))**2)/0.529177211

    
      call prepot

      call pot(r,e,1,1)

      v=e(1)*27.211386

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=v
        enddo
      else
        write (*,*) 'Only energy is available'
      endif

      endsubroutine

      SUBROUTINE PREPOT
      implicit double precision (a-h,o-z)
      parameter (n2=5000,nmax=800,mmax=5000,htkcal=627.5096d0)
      save
      dimension rvp(nt, 3), evp(nt),xxx(3),vp(3)


      return 

      ENTRY POT (rvp, evp, nt, nsurf)

         eshift= -0.174495770894d0    
	do 101 ivp = 1, nt
	   do 102 ivp1 = 1, 3
                  xxx(ivp1) = rvp(ivp, ivp1)
102        continue
      call bkmp2(xxx,eee,vp,-1)
      evp(ivp) = eee-eshift
  101 continue
      RETURN
      END
      subroutine surfid(h3id)
      implicit double precision (a-h,o-z)
	h3id = 950621.0
      return
      end
      subroutine bkmp2(r,Vtot,dVtot,id)
      implicit double precision (a-h,o-z)
      dimension r(3),dVtot(3)
      dimension dVlon(3),dVas(3),dVbnda(3),dVbndb(3),&
     &          dCal(3), dCas(3),dCbnda(3),dCbndb(3)
      dimension vb(2,25),cb(2,25)
      dimension dT(3),T(3)
      ipr = 0

 6999 format('invalid geometry, v,dv=',4f8.2)
      Vlon  = 0.d0
      Vas   = 0.d0
      Vbnda = 0.d0
      Vbndb = 0.d0
      Cal   = 0.d0
      Cas   = 0.d0
      cbnda = 0.d0
      cbndb = 0.d0
      do i=1,3
	 dVlon(i) = 0.d0
	 dVas(i)  = 0.d0
	 dVbnda(i)= 0.d0
	 dVbndb(i)= 0.d0
	 dCal(i)  = 0.d0
	 dCas(i)  = 0.d0
	 dCbnda(i)= 0.d0
	 dCbndb(i)= 0.d0
      end do  
      do i=1,25
	 vb(1,i) = 0.d0
	 vb(2,i) = 0.d0
	 cb(1,i) = 0.d0
	 cb(2,i) = 0.d0
      end do  
      call h3lond95( r, Vlon, dVlon )
      call vascal95( r, Vas, dVas )
      call compac95(r,icompc,T,dT)
      if(icompc.ge.1)then
         call csym95  ( r, cal, dCal )
         call casym95 ( r, Cas, dCas,  1 ,T,dT)
      end if
      call vbcb95(r,icompc,T,dT,ipr,&
     &          Vbnda,Vbndb,dVbnda,dVbndb,&
     &          Cbnda,Cbndb,dCbnda,dCbndb)
      Vtot = Vlon + Vas + Vbnda + Vbndb + Cal + Cas + Cbnda + Cbndb
      dVtot(1) =   dVlon(1) + dVas(1) + dVbnda(1) + dVbndb(1)& 
     &           +  dCal(1) + dCas(1) + dCbnda(1) + dCbndb(1)
      dVtot(2) =   dVlon(2) + dVas(2) + dVbnda(2) + dVbndb(2)& 
     &           +  dCal(2) + dCas(2) + dCbnda(2) + dCbndb(2)
      dVtot(3) =   dVlon(3) + dVas(3) + dVbnda(3) + dVbndb(3)& 
     &           +  dCal(3) + dCas(3) + dCbnda(3) + dCbndb(3)

      if(ipr.gt.0) then
	  write(7,*)
	  write(7,*) 'h3tot enter --------------------------------'
          write(7,7000) r,icompc,vtot,&
     &             vlon,vas,cal,cas,vbnda,vbndb,cbnda,cbndb
	  if(ipr.gt.1)then
		   write(7,7100) dVtot,dVlon,dVas,dVbnda,dVbndb,&
     &                           dCal,dCas,dCbnda,dCbndb
	  end if
	  write(7,7400) (vb(1,i),i=1,5),(vb(2,i),i=1,5)
	  write(7,7500) (cb(1,i),i=1,9),(cb(2,i),i=1,9)
          write(7,7999) Vbnda,Vbndb,Cbnda,Cbndb,  &
     &                  dVbnda,dvbndb,dCbnda,dCbndb
	  write(7,*) 'exiting subr.h3tot'
	  write(7,*) 'h3tot exit ---------------------------------'
      end if
      return
 7400 format('Vba values: ',5(1x,g12.6),/,'Vbb values: ',5(1x,g12.6))
 7500 format('Cba values: ',3(1x,f16.8),/,&
     &       '            ',3(1x,f16.8),/,&
     &       '            ',3(1x,f16.8),/,&
     &       'Cbb values: ',3(1x,f16.8),/,&
     &       '            ',3(1x,f16.8),/,&
     &       '            ',3(1x,f16.8))  
 7000 format(5x,'    r =',3(1x,f16.10),' icompac = ',i1,/,&
     &       5x,'Vtot  = ',f18.12,/,&
     &       5x,'Vlon  = ',f18.12,'       Vas   = ',g18.12,/,&
     &       5x,'Cal   = ',g18.12,'       Cas   = ',g18.12,/,&
     &       5x,'Vbnda = ',g18.12,'       Vbndb = ',g18.12,/,&
     &       5x,'Cbnda = ',g18.12,'       Cbndb = ',g18.12)
 7100 format(' dVtot   = ',3(1x,g18.12),/,&
     &       ' dVlon   = ',3(1x,g18.12),/,&
     &       ' dVas    = ',3(1x,g18.12),/,&
     &       ' dVbnda  = ',3(1x,g18.12),/,&
     &       ' dVbndb  = ',3(1x,g18.12),/,&
     &       ' dCal    = ',3(1x,g18.12),/,&
     &       ' dCas    = ',3(1x,g18.12),/,&
     &       ' dCbnda  = ',3(1x,g18.12),/,&
     &       ' dCbndb  = ',3(1x,g18.12))
 7750 format(3(1x,g16.10))
 7751 format(/,3(1x,g16.10))
 7999 format('   Vbnda        Vbndb        Cbnda        Cbndb ',/,&
     &        4(e12.6,2x),/,&
     &       'dVbnda: ',3(e16.10,1x),/,&
     &       'dVbndb: ',3(e16.10,1x),/,&
     &       'dCbnda: ',3(e16.10,1x),/,&
     &       'dCbndb: ',3(e16.10,1x))
      end
      subroutine triplet95(r,E3,id)
      implicit double precision (a-h,o-z)
      dimension E3(3),E1(3)
      parameter( rl=0.95d0, rr=1.15d0 )
      parameter( z1=1.d0, z2=2.d0, z4=4.d0)
      parameter(&
     &  a1= -0.0298546962,a2=-23.9604445036,a3=-42.5185569474,&
     &  a4=  2.0382390988,a5=-11.5214861455,a6=  1.5309487826,&
     &  C1= -0.4106358351531854,C2= -0.0770355790707090,&
     &  C3=  0.4303193846943223) !jun21 fit
      ipr = 0
      E3(3) = 0.d0
      if(r.ge.rr )then
         exdr = dexp( -a4*r )
         rsq  = r*r
         ra6  = r**(-a6)
         E3(1) = a1* ( a2 + r + a3*rsq + a5*ra6 )*exdr 
         ra61 = r**(-a6-z1)
         E3(2) = a1*exdr*&
     &   ( z1 -a2*a4 +(z2*a3-a4)*r -a3*a4*rsq -a5*a6*ra61 -a4*a5*ra6 )
      end if
      if(r.lt.rr )then
          dr = r- rl  
          call vh2opt95(r,e1,2)
	  if( r.le.rl ) then
            E3(1) = E1(1) + c2*dr + c3
	    E3(2) = E1(2) + c2
          else                 
            E3(1) = E1(1) + c1*dr*dr*dr + c2*dr + c3
	    E3(2) = E1(2) + 3.d0*c1*dr*dr + c2
	  end if
      end if
      if(ipr.gt.0) write(7,7100) e3
 7100 format('  triplet:  E3 = ',3(1x,f12.8))
      return
      end
      subroutine vH2opt95(r,E,ideriv)  
      implicit double precision (a-h,o-z)
      dimension E(3)
      parameter(a0=0.03537359271649620,    a1= 2.013977588700072    ,&
     &  a2= -2.827452449964767        ,    a3= 2.713257715593500    ,&
     &  a4= -2.792039234205731        ,    a5= 2.166542078766724    ,&
     &  a6= -1.272679684173909        ,    a7= 0.5630423099212294   ,&
     &  a8= -0.1879397372273814       ,    a9= 0.04719891893374140  ,&
     & a10= -0.008851622656489644     ,   a11= 0.001224998776243630 ,&
     & a12= -1.227820520228028d-04    ,   a13= 8.638783190083473d-06,&
     & a14= -4.036967926499151d-07    ,   a15= 1.123286608335365d-08,&
     & a16= -1.406619156782167d-10 )
      parameter( r0=3.5284882d0,  DD=0.160979391d0,&
     &           c6=6.499027d0,   c8=124.3991d0,    c10=3285.828d0)
      E(1)  = -999.d0
      E(2)  = -999.d0
      E(3)  = -999.d0
      r2  = r  *r
      r3  = r2 *r
      r4  = r3 *r
      r5  = r4 *r
      r6  = r5 *r
      r7  = r6 *r
      r8  = r7 *r
      r9  = r8 *r
      r10 = r9 *r
      r11 = r10*r
      r12 = r11*r
      r13 = r12*r
      r14 = r13*r
      r15 = r14*r
      r02 = r0*r0
      r04 = r02*r02
      r06 = r04*r02
      rr2  = r2 + r02
      rr4  = r4 + r04
      rr6  = r6 + r06
      rr25 = rr2*rr2*rr2*rr2*rr2
      alphaR =  a0/r + a1&
     &            + a2 *r   + a3 *r2  + a4 *r3  + a5 *r4  + a6 *r5&
     &            + a7 *r6  + a8 *r7  + a9 *r8  + a10*r9  + a11*r10&
     &            + a12*r11 + a13*r12 + a14*r13 + a15*r14 + a16*r15
      exalph = dexp(alphaR)
      vsr = DD*(exalph-1.d0)*(exalph-1.d0) - DD
      vlr =  -C6/rr6   -C8/(rr4*rr4)   -C10/rr25
      E(1)=  vsr +  vlr 
      if(ideriv.ge.1)then        
         r3 = r2*r
         r5 = r4*r
         rr26 = rr25*rr2
         rr43 = rr4*rr4*rr4
	 dalphR = -a0/r2 + a2 + 2.d0*a3*r&
     &          +  3.d0 *a4 *r2  +  4.d0 *a5 *r3  +  5.d0 *a6 *r4 &
     &          +  6.d0 *a7 *r5  +  7.d0 *a8 *r6  +  8.d0 *a9 *r7 &
     &          +  9.d0 *a10*r8  + 10.d0 *a11*r9  + 11.d0 *a12*r10&
     &          + 12.d0 *a13*r11 + 13.d0 *a14*r12 + 14.d0 *a15*r13&
     &          + 15.d0 *a16*r14
	 dvsr = 2.d0 *DD *(exalph-1.d0) *exalph *dalphR
	 dvlr =    6.d0*C6*r5 / (rr6*rr6)&
     &          +  8.d0*C8*r3 / rr43  + 10.d0*C10*r / rr26
	 E(2) = dvsr + dvlr
      end if
      if(ideriv.ge.2)then
         r10  = r6*r4
         rr27 = rr26*rr2
         rr44 = rr43*rr4
         rr62 = rr6*rr6
         rr63 = rr62*rr6
	 ddalph = 2.d0*a0/r3 +  2.d0*a3      +  6.d0*a4 *r&
     &     + 12.d0*a5 *r2    + 20.d0*a6 *r3  + 30.d0*a7 *r4&
     &     + 42.d0*a8 *r5    + 56.d0*a9 *r6  + 72.d0*a10*r7&
     &     + 90.d0*a11*r8    +110.d0*a12*r9  +132.d0*a13*r10&
     &     +156.d0*a14*r11   +182.d0*a15*r12 +210.d0*a16*r13
	 ddvsr = 2.d0 *DD *exalph &
     &   *( (2.d0*exalph-1.d0)*dalphR*dalphR + (exalph-1.d0)*ddalph )
	 ddvlr =- 72.d0* C6 *r10 / rr63 -96.d0 *C8 *r6 / rr44  &
     &		-120.d0*C10 *r2  / rr27 +30.d0 *C6 *r4 / rr62&
     &		+ 24.d0* C8 *r2  / rr43 +10.d0 *C10    / rr26
	 E(3) = ddvsr + ddvlr
      end if
      return
      end
      subroutine h3lond95(r,Vlon,dVlon)
      implicit double precision (a-h,o-z)
      real*8 Q(3),J(3),Jt
      dimension r(3),e1(3),e3(3),esing(3),etrip(3),dVlon(3)
      dimension de1(3),de3(3)
      parameter( half=0.5d0, two=2.d0 , eps2=1.d-12 )
      ipr = 0
      do i=1,3
         call vh2opt95(r(i),esing,2)
	  e1(i) = esing(1)
	 de1(i) = esing(2)
	 call triplet95(r(i),etrip,2)
	 if(ipr.gt.0) write(7,7400) i,esing,etrip
	  e3(i) = etrip(1)
	 de3(i) = etrip(2)
         q(i)  = half*(E1(i) + E3(i))  
         j(i)  = half*(E1(i) - E3(i)) 
      end do   
      sumQ  =   q(1) + q(2) + q(3)
      sumJ  =   dabs( j(2)-j(1) )**2 &
     &        + dabs( j(3)-j(2) )**2 &
     &        + dabs( j(3)-j(1) )**2 
      Jt     = half*sumJ + eps2
      rootJt = dsqrt(Jt)
      Vlon   = sumQ - rootJt   
      if(ipr.gt.0) then
	 write(7,7410) sumq,sumj
         write(7,7420) vlon,rootJt
      end if
      dVlon(1) = half*(dE1(1)+de3(1))&
     &         - 0.25d0*(two*j(1)-j(2)-j(3))*(dE1(1)-de3(1))/rootJt
      dVlon(2) = half*(dE1(2)+de3(2))&
     &         - 0.25d0*(two*j(2)-j(3)-j(1))*(dE1(2)-de3(2))/rootJt
      dVlon(3) = half*(dE1(3)+de3(3))&
     &         - 0.25d0*(two*j(3)-j(1)-j(2))*(dE1(3)-de3(3))/rootJt
      if(ipr.gt.0) then
         write(7,7000) r,e1,e3,vlon
	 write(7,7100) q,j
	 write(7,7200) dVlon
      end if
 7000 format('             r = ',3(1x,f12.6),/,&
     &       '            E1 = ',3(1x,f12.8),/,&
     &       '            E3 = ',3(1x,f12.8),/,&
     &       '          vlon = ',1x,f12.8)
 7100 format(13x,'q = ',3(1x,e12.6),/,13x,'j = ',3(1x,e12.6))
 7200 format('         dVlon = ',3(1x,g12.6))
 7400 format('from subr.london: ',/,&
     &       '  using r',i1,':   Esinglet=',3(1x,f12.8),/,&
     &       '         ',1x,'    Etriplet=',3(1x,f12.8))
 7410 format('         sumQ = ',g12.6,'        sumJ = ',g12.6)
 7420 format('         Vlon = ',f12.8,'      rootJt = ',g12.6)
      return
      end
      subroutine Vascal95(rpass,Vas,dVas)
      implicit double precision (a-h,o-z)
      dimension dVas(3),dA(3),dS(3),rpass(3)
      parameter( &
     & aa1=0.3788951192E-02, aa2=0.1478100901E-02,&
     & aa3=-.1848513849E-03, aa4=0.9230803609E-05,&
     & aa5=-.1293180255E-06, aa6=0.5237179303E+00,&
     & aa7=-.1112326215E-02) !jun21 fit
      ipr = 0
      r1 = rpass(1)
      r2 = rpass(2)
      r3 = rpass(3)
      R  = r1 + r2 + r3
      Rsq = R*R
      Rcu = Rsq*R
      call acalc95(r1,r2,r3,A,dA)
        A2 = A *A
        A3 = A2*A
        A4 = A3*A
        A5 = A4*A
        exp1 = dexp(-aa1*Rcu)
	exp6 = dexp(-aa6*R)
        S = aa2*A2 + aa3*A3 + aa4*A4 + aa5*A5
      Vas = S*exp1 +  aa7*A2 *exp6 / R
         dS(1)  = ( 2.d0*aa2*A  +3.d0*aa3*A2 &
     &             +4.d0*aa4*A3 +5.d0*aa5*A4) * dA(1)
         dVas(1) =  -3.d0*aa1*Rsq*S*exp1 + dS(1)*exp1&
     &             -aa7*A2*exp6/Rsq  + 2.d0*aa7*A*dA(1)*exp6/R&
     &             -aa6*aa7*A2*exp6/R
         dS(2)  = ( 2.d0*aa2*A  +3.d0*aa3*A2 &
     &             +4.d0*aa4*A3 +5.d0*aa5*A4) * dA(2)
         dVas(2) =  -3.d0*aa1*Rsq*S*exp1 + dS(2)*exp1&
     &             -aa7*A2*exp6/Rsq  + 2.d0*aa7*A*dA(2)*exp6/R&
     &             -aa6*aa7*A2*exp6/R
         dS(3)  = ( 2.d0*aa2*A  +3.d0*aa3*A2 &
     &             +4.d0*aa4*A3 +5.d0*aa5*A4) * dA(3)
         dVas(3) =  -3.d0*aa1*Rsq*S*exp1 + dS(3)*exp1&
     &             -aa7*A2*exp6/Rsq  + 2.d0*aa7*A*dA(3)*exp6/R&
     &             -aa6*aa7*A2*exp6/R
      if(ipr.gt.0) write(7,7100) Vas,dS,dVas
 7100 format('    from subr.vascal:   Vas  = ',  1x,g12.6, /,&
     &       '                        dS   = ',3(1x,g12.6),/,&
     &       '                        dVas = ',3(1x,g12.6))
      return
      end
      subroutine acalc95(r1,r2,r3,A,dA)
      implicit double precision (a-h,o-z)
      dimension dA(3)
      ipr = 0
      A = (r1-r2)*(r2-r3)*(r3-r1)
      dA(1) = ( -2.d0*r1 + r2 + r3 )*(r2-r3)
      dA(2) = ( -2.d0*r2 + r3 + r1 )*(r3-r1)
      dA(3) = ( -2.d0*r3 + r1 + r2 )*(r1-r2)
      if(A.lt.0.d0)then
         A = -A
         dA(1) = -dA(1)
         dA(2) = -dA(2)
         dA(3) = -dA(3)
      end if
      if(ipr.gt.0) write(7,7100) A,dA
 7100 format('Acalc enter ---------------------------------------',/,&
     &       '     A = ',f15.10,'   dA = ',3(1x,g12.6),/,&
     &       'Acalc exit ----------------------------------------')
      return
      end 
      subroutine compac95(r,icompc,T,dT)
      implicit double precision (a-h,o-z)
      dimension r(3),T(3),dT(3)
      parameter( rr = 1.15d0, rp = 1.25d0 )
      do i=1,3
	 T(i) = 0.d0
	 dT(i)= 0.d0
      end do   
      icompc = 0
      if(r(1).lt.rr) icompc = icompc + 1
      if(r(2).lt.rr) icompc = icompc + 1
      if(r(3).lt.rr) icompc = icompc + 1
      if(icompc.eq.0) return
      do i=1,3
      if(r(i).lt.rr) then
	 top = rr-r(i)
	 bot = rp-r(i)
	 top2 = top  * top
	 top3 = top2 * top
	 bot2 = bot  * bot
         T(i)  = top3/bot
         dT(i) = -3.d0*top2/bot + top3/bot2
      end if
      end do   
      return
      end
      subroutine casym95( r,Cas,dCas,id,T,dT )
      implicit double precision (a-h,o-z)
      dimension r(3),T(3)
      dimension dPr(3),dT(3),dsumT(3),dterm1(3),determ(3),dseries(3),&
     &          dCas(3),dA(3)
      parameter( &
     & u1=0.2210243144E+00, u2=0.4367417579E+00, u3=0.6994985432E-02,&
     & u4=0.1491096501E+01, u5=0.1602896673E+01, u6=-.2821747323E+01,&
     & u7=0.4948310833E+00, u8=-.3540394679E-01, u9=-.3305809954E+01,&
     &u10=0.3644382172E+01,u11=-.9997570970E+00,u12=0.7989919534E-01,&
     &u13=-.1075807322E-02) !jun21 fit
      ipr = 0
      cas = 0.d0
      dCas(1) = 0.d0
      dCas(2) = 0.d0
      dCas(3) = 0.d0
      call acalc95(r(1),r(2),r(3),A,dA)
      A2 = A*A
      sumT = T(1) + T(2) + T(3)
      Sr   = r(1) + r(2) + r(3)
      Pr   = r(1) * r(2) * r(3)
      Sr2  = Sr*Sr
      Sr3  = Sr2*Sr
      Pr2  = Pr*Pr
      Pr3  = Pr2*Pr
      series = 1.d0 + u4/Pr2 +  u5/Pr + u6 +   u7*Pr +  u8*Pr2&
     &           + A*(u9/Pr2 + u10/Pr + u11 + u12*Pr + u13*Pr2)
      term1  = u1/Pr**u2
      eterm  = dexp(-u3*Sr3)
      Cas    = sumt * a2 * term1 * series * eterm
      if(ipr.gt.0) then
	    write(7,*)
            write(7,*) 'Casym: ------------------'
	    write(7,*) ' T(i) = ',T
	    write(7,*) 'dT(i) = ',dT
	    write(7,*) ' id = ',id
	    write(7,7400) series,term1,eterm,Cas
	    write(7,*) ' Cas = ',Cas
      end if
      if(id.gt.0)then 
	    dPr(1) = r(2)*r(3)
	    dPr(2) = r(3)*r(1)
	    dPr(3) = r(1)*r(2)
	 do i=1,3
	    dsumt(i) = dt(i)
	    dterm1(i)=  -1.d0*u1*u2*Pr**(-u2-1.d0)*dPr(i)
	    determ(i)= eterm *(-3.d0*u3*Sr2)
	    dseries(i)=&
     &        dPr(i)*( -2.d0*u4/Pr3  -u5/Pr2 +u7  +2.d0*u8*Pr    )&
     &        +dA(i)*( u9/Pr2 +u10/Pr +u11 +u12*Pr +u13*Pr2    )&
     &        +A*dPr(i)*( -2.d0*u9/Pr3 -u10/Pr2 +u12 + 2.d0*u13*Pr )
	    dCas(i) = &
     &        dsumt(i)      * a2    * term1 * series * eterm&
     &      + 2.d0*a*dA(i)  * sumt  * term1 * series * eterm&
     &      + dterm1(i)     * sumt  * a2    * series * eterm&
     &      + dseries(i)    * sumt  * a2    * term1  * eterm&
     &      + determ(i)     * sumt  * a2    * term1  * series
         end do   
	 if(ipr.gt.0) then
	    write(7,7500) dPr,dt,dterm1,determ,dseries,dCas
	    write(7,*) ' dCas = ',dCas
	 end if
      end if
      return
 7400 format('subr.Cas:   series = ',g12.6,'     term1 = ',g12.6,/,&
     &       '             eterm = ',g12.6,'      Cas  = ',g12.6)
 7500 format('    derivatives: dPr = ',3(1x,g12.6),/,&
     &       '                 dT  = ',3(1x,g12.6),/,&
     &       '              dterm1 = ',3(1x,g12.6),/,&
     &       '              determ = ',3(1x,g12.6),/,&
     &       '             dseries = ',3(1x,g12.6),/,&
     &       '              dCas   = ',3(1x,g12.6))  
      end
      subroutine csym95(r,Cal,dCal)
      implicit double precision (a-h,o-z)
      dimension r(3),dCal(3),G(3)
      dimension dG(3),sumv(3)
      parameter( rr=1.15d0, rp=1.25d0 )
      parameter( &
     & v1=-.2071708868E+00, v2=-.5672350377E+00, v3=0.9058780367E-02) !jun21 fit
      cal   = 0.d0
      ipr   = 0
      Sr    = r(1)+r(2)+r(3)
      Sr2   = Sr*Sr
      Sr3   = Sr*Sr2  
      exp3  = dexp( -v3*Sr3 )
      dexp3 = -3.d0*v3*Sr2*exp3
      do i=1,3
         ri      = r(i)
	 rrri    = rr-ri
	 rrri2   = rrri*rrri
	 rrri3   = rrri*rrri2
	 rpri    = rp-ri
	 rpri2   = rpri*rpri
         G(i)    = 0.d0
	 dG(i)   = 0.d0
	 sumv(i) = v1+v1*v2*ri
         if(ri.lt.rr) then
	    G(i)  = (rrri3/rpri) *sumv(i)
  	    dG(i) = (rrri3/rpri2)*sumv(i)& 
     &                -3.d0*(rrri2/rpri)*sumv(i)&
     &                + (rrri3/rpri)*v1*v2
	 end if
      end do   
      sumG = G(1) + G(2) + G(3)
      cal  = sumG*exp3
      dCal(1) = dG(1)*exp3 + sumG*dexp3
      dCal(2) = dG(2)*exp3 + sumG*dexp3
      dCal(3) = dG(3)*exp3 + sumG*dexp3
      if(ipr.gt.0)then
	 write(7,*)
         write(7,7100) Cal
	 write(7,7000) G,dG,sumv,dCal
	 write(7,*) 'exiting subr.csym'
      end if
 7100 format('Csym enter ----------------------------------------',/,&
     &       '    Cal = ',1x,g20.14)
 7000 format('      G = ',3(1x,g12.6),/,&
     &       '     dG = ',3(1x,g12.6),/,&
     &       '   sumv = ',3(1x,g12.6),/,&
     &       '   dCal = ',3(1x,g16.10),/,&
     &       'Csym exit -----------------------------------------')
      return
      end
      subroutine vbcb95(rpass,icompc,T,dT,ipr,&
     &                  Vbnda,Vbndb,dVbnda,dVbndb,&
     &                  Cbnda,Cbndb,dCbnda,dCbndb)
      implicit double precision (a-h,o-z)
      dimension rpass(3),dB1a(3),dB1b(3),dB2(3),dB3(6)
      dimension dVbnda(3),dVbndb(3)
      dimension dVb1(3),dvb2(3),dvb3(3),dvb4(3),dvb5(3)
      dimension Vbnd(2)
      dimension vb(2,25)   
      dimension Cbnds(2),dCbnds(2,3)
      dimension dCb1(3),dCb2(3),dCb3(3),dCb4(3),dCb5(3),dCb6(3),&
     &          dCb7(3),dCb8(3)
      dimension T(3),dP(3)
      dimension dT(3),dCbnda(3),dCbndb(3)
      dimension cb(2,25)   
      parameter( z58=0.625d0, z38=0.375d0 )  
      parameter(beta1=0.52d0, beta2=0.052d0, beta3=0.79d0 ) 
      parameter( &
     &a11=-.1838073394E+03,a12=0.1334593242E+02,a13=-.2358129537E+00,&
     &a21=-.4668193478E+01,a22=0.7197506670E+01,a23=0.2162004275E+02,&
     &a24=0.2106294028E+02,a31=0.4242962586E+01,a32=0.4453505045E+01,&
     &a41=-.1456918088E+00,a42=-.1692657366E-01,a43=0.1279520698E+01,&
     &a44=-.4898940075E+00,a51=0.1742295219E+03,a52=0.3142175348E+02,&
     &a53=0.5152903406E+01) !jun21 fit
      parameter(& 
     &g11=-.4765732725E+02,g12=0.3648933563E+01,g13=-.7141145244E-01,&
     &g21=0.1002349176E-01,g22=0.9989856329E-02,g23=-.4161953634E-02,&
     &g24=0.9075807910E-03,g31=-.2693628729E+00,g32=-.1399065763E-01,&
     &g41=-.1417634346E-01,g42=-.4870024792E-03,g43=0.1312231847E+00,&
     &g44=-.4409850519E-01,g51=0.5382970863E+02,g52=0.4587102824E+01,&
     &g53=0.1768550515E+01) !jun21 fit
      parameter( &
     &c11=0.1860299931E+04,c12=-.6134458037E+03,c13=0.7337207161E+02,&
     &c14=-.2676717625E+04,c15=0.1344099415E+04,c21=0.1538913137E+03,&
     &c22=0.4348007369E+02,c23=0.1719720677E+03,c24=0.2115963042E+03,&
     &c31=-.7026089414E+02,c32=-.1300938992E+03,c41=0.1310273564E+01,&
     &c42=-.6175149574E+00,c43=-.2679089358E+02,c44=0.5577477171E+01,&
     &c51=-.3543353539E+04,c52=-.3740709591E+03,c53=0.7979303144E+02,&
     &c61=-.1104230585E+04,c62=0.4603572025E+04,c63=-.5593496634E+04,&
     &c71=-.1069406434E+02,c72=0.1021807153E+01,c73=0.6669828341E-01,&
     &c74=0.4168542348E+02,c75=0.1751608567E+02,c81=0.9486883238E+02,&
     &c82=-.1519334221E+02,c83=0.4024697252E+04,c84=-.2225159395E+02) !jun21
      parameter( &
     &d11=0.4203543357E+03,d12=-.4922474096E+02,d13=0.3362942544E+00,&
     &d14=-.3827423082E+03,d15=0.1746726001E+03,d21=0.1699995737E-01,&
     &d22=0.1513036778E-01,d23=0.2659119354E-01,d24=-.5760387483E-02,&
     &d31=0.1020622621E+02,d32=0.1050536271E-01,d41=0.6836172780E+00,&
     &d42=-.1627858240E+00,d43=-.6925485045E+01,d44=0.1632567385E+01,&
     &d51=0.1083595009E+04,d52=0.4641431791E+01,d53=-.8233144461E+00,&
     &d61=-.6157225942E+02,d62=0.3094361471E+03,d63=-.3299631143E+03,&
     &d71=0.8866227120E+01,d72=-.1382126854E+01,d73=0.7620770145E-01,&
     &d74=-.5145757859E+02,d75=0.2046097265E+01,d81=0.2540775558E+01,&
     &d82=-.4889246569E+00,d83=-.1127439280E+04,d84=-.2269932295E+01) !jun21
      cx1 = c51 + c83
      dx1 = d51 + d83
      r1 = rpass(1)
      r2 = rpass(2)
      r3 = rpass(3)
      t1 = r1*r1 -r2*r2 -r3*r3 
      t2 = r2*r2 -r3*r3 -r1*r1 
      t3 = r3*r3 -r1*r1 -r2*r2 
      c1 = t1 / (-2.d0*r2*r3)
      c2 = t2 / (-2.d0*r3*r1)
      c3 = t3 / (-2.d0*r1*r2)
      sum = c1 + c2 + c3
      B1a = 1.d0 - sum         
      c1cube = c1*c1*c1
      c2cube = c2*c2*c2
      c3cube = c3*c3*c3
      cos3t1 = 4.d0*c1cube - 3.d0*c1
      cos3t2 = 4.d0*c2cube - 3.d0*c2
      cos3t3 = 4.d0*c3cube - 3.d0*c3
      sumb   = cos3t1 + cos3t2 + cos3t3 
      B1b    = 1.d0 - ( z58*sumb + z38*sum )
         dc1dr1 = -r1/(r2*r3)
         dc2dr2 = -r2/(r1*r3)
         dc3dr3 = -r3/(r1*r2)
         dc1dr2 = ( t1/(r2*r2) + 2.d0 )/(2.d0*r3)
         dc1dr3 = ( t1/(r3*r3) + 2.d0 )/(2.d0*r2)
         dc2dr1 = ( t2/(r1*r1) + 2.d0 )/(2.d0*r3)
         dc2dr3 = ( t2/(r3*r3) + 2.d0 )/(2.d0*r1)
         dc3dr1 = ( t3/(r1*r1) + 2.d0 )/(2.d0*r2)
         dc3dr2 = ( t3/(r2*r2) + 2.d0 )/(2.d0*r1)
         db1a(1) = -1.d0*( dc1dr1 + dc2dr1 + dc3dr1 )
         db1a(2) = -1.d0*( dc1dr2 + dc2dr2 + dc3dr2 )
         db1a(3) = -1.d0*( dc1dr3 + dc2dr3 + dc3dr3 )
         d1 = 12.d0*c1*c1 - 3.d0
         d2 = 12.d0*c2*c2 - 3.d0
         d3 = 12.d0*c3*c3 - 3.d0
         db1b(1) = -z58*( d1*dc1dr1 + d2*dc2dr1 + d3*dc3dr1 )&
     &             -z38*(    dc1dr1 +    dc2dr1 +    dc3dr1 )
         db1b(2) = -z58*( d1*dc1dr2 + d2*dc2dr2 + d3*dc3dr2 )&
     &             -z38*(    dc1dr2 +    dc2dr2 +    dc3dr2 )
         db1b(3) = -z58*( d1*dc1dr3 + d2*dc2dr3 + d3*dc3dr3 )&
     &             -z38*(    dc1dr3 +    dc2dr3 +    dc3dr3 )
      R    = r1 + r2 + r3
      Rsq  = R*R
      B2   = 1.d0/r1 + 1.d0/r2 + 1.d0/r3
      B3   = (r2-r1)*(r2-r1) + (r3-r2)*(r3-r2) + (r1-r3)*(r1-r3)
      eps2 = 1.d-12
      B3b  = dsqrt( B3 + eps2 )
	 dB2(1) = -1.d0/(r1*r1)
	 dB2(2) = -1.d0/(r2*r2)
	 dB2(3) = -1.d0/(r3*r3)
	 dB3(1) = 4.d0*r1 -2.d0*r2 -2.d0*r3
	 dB3(2) = 4.d0*r2 -2.d0*r3 -2.d0*r1
	 dB3(3) = 4.d0*r3 -2.d0*r1 -2.d0*r2
	 dB3(4) = 0.5d0*dB3(1)/B3b
	 dB3(5) = 0.5d0*dB3(2)/B3b
	 dB3(6) = 0.5d0*dB3(3)/B3b
         exp1   = dexp(-beta1*R)
         exp2   = dexp(-beta2*Rsq)
         exp7   = dexp(-beta3*R)
	 dexp1  = -beta1*exp1
	 dexp2  = -2.d0*beta2*R*exp2
	 dexp7  = -beta3*exp7
         B1    = B1a
         B12   = B1*B1
         B13   = B12*B1
         B14   = B13*B1
         B15   = B14*B1
         asum  = a11 +  a12*R +  a13*Rsq 
         bsum  = a21*B12 +  a22*B13 +  a23*B14 +  a24*B15
	 csum  = a31*B1*exp1 +  a32*B12*exp2
         dsum1 = a41*exp1 +  a42*exp2
         dsum2 = a43*exp1 +  a44*exp2
         fsum  = a51 +  a52*R +  a53*Rsq 
         vb(1,1) = B1*asum*exp1 
         vb(1,2) = bsum * exp2
         vb(1,3) = B2*csum
         vb(1,4) = B1*B3*dsum1 + B1*B3b*dsum2
         vb(1,5) = B1* fsum  *exp7 
         Vbnd(1) = vb(1,1)+vb(1,2)+vb(1,3)+vb(1,4)+vb(1,5)
	 dasum   = a12+2.d0*a13*R
	 dbsum   =  2.d0* a21*B1  + 3.d0* a22*B12&
     &            + 4.d0* a23*B13 + 5.d0* a24*B14
	 ddsum1  =  a41*dexp1 +  a42*dexp2
	 ddsum2  =  a43*dexp1 +  a44*dexp2
	 dfsum   =  a52 +2.d0*a53*R
         dVb1(1) = dB1a(1)*asum*exp1 + B1*dasum*exp1 +B1*asum*dexp1
         dVb1(2) = dB1a(2)*asum*exp1 + B1*dasum*exp1 +B1*asum*dexp1
         dVb1(3) = dB1a(3)*asum*exp1 + B1*dasum*exp1 +B1*asum*dexp1
         dVb2(1) = dbsum*dB1a(1)*exp2 + bsum*dexp2
         dVb2(2) = dbsum*dB1a(2)*exp2 + bsum*dexp2
         dVb2(3) = dbsum*dB1a(3)*exp2 + bsum*dexp2
         dVb3(1)=dB2(1)*csum + B2*(  a31*dB1a(1)*exp1 +  a31*B1*dexp1&
     &                  +2.d0* a32*B1*dB1a(1)*exp2 +  a32*B12*dexp2)
         dVb3(2)=dB2(2)*csum + B2*(  a31*dB1a(2)*exp1 +  a31*B1*dexp1&
     &                  +2.d0* a32*B1*dB1a(2)*exp2 +  a32*B12*dexp2)
         dVb3(3)=dB2(3)*csum + B2*(  a31*dB1a(3)*exp1 +  a31*B1*dexp1&
     &                  +2.d0* a32*B1*dB1a(3)*exp2 +  a32*B12*dexp2)
         dVb4(1)= dB1a(1)*B3 *dsum1+  B1*dB3(1)*dsum1 + B1*B3 *ddsum1&
     &          + dB1a(1)*B3b*dsum2 + B1*dB3(4)*dsum2 + B1*B3b*ddsum2
         dVb4(2)= dB1a(2)*B3 *dsum1+  B1*dB3(2)*dsum1 + B1*B3 *ddsum1&
     &          + dB1a(2)*B3b*dsum2 + B1*dB3(5)*dsum2 + B1*B3b*ddsum2
         dVb4(3)= dB1a(3)*B3 *dsum1+  B1*dB3(3)*dsum1 + B1*B3 *ddsum1&
     &          + dB1a(3)*B3b*dsum2 + B1*dB3(6)*dsum2 + B1*B3b*ddsum2
	 dVb5(1) = dB1a(1)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
	 dVb5(2) = dB1a(2)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
	 dVb5(3) = dB1a(3)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
         dVbnda(1) = dvb1(1)+dvb2(1)+dvb3(1)+dvb4(1)+dvb5(1)
         dVbnda(2) = dvb1(2)+dvb2(2)+dvb3(2)+dvb4(2)+dvb5(2)
         dVbnda(3) = dvb1(3)+dvb2(3)+dvb3(3)+dvb4(3)+dvb5(3)
	 Vbnda = Vbnd(1)
         B1    = B1b
         B12   = B1*B1
         B13   = B12*B1
         B14   = B13*B1
         B15   = B14*B1
         asum  = g11 +  g12*R +  g13*Rsq 
         bsum  = g21*B12 +  g22*B13 +  g23*B14 +  g24*B15
	 csum  = g31*B1*exp1 +  g32*B12*exp2
         dsum1 = g41*exp1 +  g42*exp2
         dsum2 = g43*exp1 +  g44*exp2
         fsum  = g51 +  g52*R +  g53*Rsq 
         vb(2,1) = B1*asum*exp1 
         vb(2,2) = bsum * exp2
         vb(2,3) = B2*csum
         vb(2,4) = B1*B3*dsum1 + B1*B3b*dsum2
         vb(2,5) = B1* fsum  *exp7 
         Vbnd(2) = vb(2,1)+vb(2,2)+vb(2,3)+vb(2,4)+vb(2,5)
	 dasum   = g12+2.d0*g13*R
	 dbsum   =  2.d0* g21*B1  + 3.d0* g22*B12&
     &             +4.d0* g23*B13 + 5.d0* g24*B14
	 ddsum1  =  g41*dexp1 +  g42*dexp2
	 ddsum2  =  g43*dexp1 +  g44*dexp2
	 dfsum   =  g52 +2.d0*g53*R
         dVb1(1) = dB1b(1)*asum*exp1 + B1*dasum*exp1 +B1*asum*dexp1
         dVb1(2) = dB1b(2)*asum*exp1 + B1*dasum*exp1 +B1*asum*dexp1
         dVb1(3) = dB1b(3)*asum*exp1 + B1*dasum*exp1 +B1*asum*dexp1
         dVb2(1) = dbsum*dB1b(1)*exp2 + bsum*dexp2
         dVb2(2) = dbsum*dB1b(2)*exp2 + bsum*dexp2
         dVb2(3) = dbsum*dB1b(3)*exp2 + bsum*dexp2
         dVb3(1)=dB2(1)*csum + B2*(  g31*dB1b(1)*exp1 +  g31*B1*dexp1&
     &                  +2.d0* g32*B1*dB1b(1)*exp2 +  g32*B12*dexp2)
         dVb3(2)=dB2(2)*csum + B2*(  g31*dB1b(2)*exp1 +  g31*B1*dexp1&
     &                  +2.d0* g32*B1*dB1b(2)*exp2 +  g32*B12*dexp2)
         dVb3(3)=dB2(3)*csum + B2*(  g31*dB1b(3)*exp1 +  g31*B1*dexp1&
     &                  +2.d0* g32*B1*dB1b(3)*exp2 +  g32*B12*dexp2)
         dVb4(1)= dB1b(1)*B3 *dsum1+  B1*dB3(1)*dsum1 + B1*B3 *ddsum1&
     &          + dB1b(1)*B3b*dsum2 + B1*dB3(4)*dsum2 + B1*B3b*ddsum2
         dVb4(2)= dB1b(2)*B3 *dsum1+  B1*dB3(2)*dsum1 + B1*B3 *ddsum1&
     &          + dB1b(2)*B3b*dsum2 + B1*dB3(5)*dsum2 + B1*B3b*ddsum2
         dVb4(3)= dB1b(3)*B3 *dsum1+  B1*dB3(3)*dsum1 + B1*B3 *ddsum1&
     &          + dB1b(3)*B3b*dsum2 + B1*dB3(6)*dsum2 + B1*B3b*ddsum2
	 dVb5(1) = dB1b(1)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
	 dVb5(2) = dB1b(2)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
	 dVb5(3) = dB1b(3)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
         dVbndb(1) = dvb1(1)+dvb2(1)+dvb3(1)+dvb4(1)+dvb5(1)
         dVbndb(2) = dvb1(2)+dvb2(2)+dvb3(2)+dvb4(2)+dvb5(2)
         dVbndb(3) = dvb1(3)+dvb2(3)+dvb3(3)+dvb4(3)+dvb5(3)
	 Vbndb = Vbnd(2)
      if(icompc.eq.0) return
      sumT  = T(1) + T(2) + T(3)
      Rcu   = Rsq*R
      P     = r1*r2*r3
      Psq   = P*P
      Pcu   = Psq*P
      dP(1) = r2*r3
      dP(2) = r3*r1
      dP(3) = r1*r2
      Cbnds(1)  = 0.d0
      Cbnds(2)  = 0.d0
      dCbnda(1) = 0.d0
      dCbnda(2) = 0.d0
      dCbnda(3) = 0.d0
      dCbndb(1) = 0.d0
      dCbndb(2) = 0.d0
      dCbndb(3) = 0.d0
            exp7  = dexp( -beta2*Rcu ) 
	    dexp7 = -3.d0*beta2*Rsq*exp7
       B1    = B1a
       B12   = B1*B1
       B13   = B12*B1
       B14   = B13*B1
       B15   = B14*B1
       asum  = c11 + c12*R + c13*Rsq + c14/R + c15/Rsq 
       bsum  = c21*B12 + c22*B13 + c23*B14 + c24*B15
       csum  = c31*B1 *exp1  + c32*B12*exp2  
       dsum1 = c41*exp1 + c42*exp2
       dsum2 = c43*exp1 + c44*exp2
       fsum  = cx1 + c52*R + c53*Rsq 
       gsum  = c61 + c62/R + c63/Rsq 
       aasum = c71 + c72*P + c73*Psq + c74/P + c75/Psq 
       ffsum =       c81*P + c82*Psq + c84/Psq 
       dasum = c12 + 2.d0*c13*R -c14/Rsq -2.d0*c15/Rcu
       dbsum = 2.d0*c21*B1 +3.d0*c22*B12 +4.d0*c23*B13 +5.d0*c24*B14
       ddsum1= c41*dexp1 + c42*dexp2
       ddsum2= c43*dexp1 + c44*dexp2
       dfsum = c52 + 2.d0*c53*R
       dgsum = -c62/Rsq - 2.d0*c63/Rcu
       daasum = c72 +2.d0*c73*P -c74/Psq -2.d0*c75/Pcu
       dffsum = c81 + 2.d0*c82*P -2.d0*c84/Pcu
         cb(1,1)  = B1*asum * exp1  /P
         cb(1,2)  = bsum * exp2  
         cb(1,3)  = B2*csum 
         cb(1,4)  = B1*B3*dsum1 + B1*B3b*dsum2
         cb(1,5)  = B1 * fsum * exp7 /P
         cb(1,6)  = B1 * gsum * exp7  
         cb(1,7)  = B1 * aasum * exp2 
         cb(1,8)  = B1 * ffsum * exp7 
         Cbnds(1) = cb(1,1)+cb(1,2)+cb(1,3)+cb(1,4)+cb(1,5)&
     &                     +cb(1,6)+cb(1,7)+cb(1,8)
	 dCb1(1) = dB1a(1)*asum*exp1/P + B1*dasum*exp1/P&
     &                +B1*asum*dexp1/P - B1*asum*exp1*dP(1)/Psq
	 dCb1(2) = dB1a(2)*asum*exp1/P + B1*dasum*exp1/P&
     &                +B1*asum*dexp1/P - B1*asum*exp1*dP(2)/Psq
	 dCb1(3) = dB1a(3)*asum*exp1/P + B1*dasum*exp1/P&
     &                +B1*asum*dexp1/P - B1*asum*exp1*dP(3)/Psq
	 dCb2(1) = dbsum*dB1a(1)*exp2 + bsum*dexp2
	 dCb2(2) = dbsum*dB1a(2)*exp2 + bsum*dexp2
	 dCb2(3) = dbsum*dB1a(3)*exp2 + bsum*dexp2
	 dCb3(1) = dB2(1)*csum + B2*( c31*dB1a(1)*exp1 +c31*B1*dexp1&
     &            +c32 *2.d0*B1*dB1a(1)*exp2 + c32 *B12*dexp2 )
	 dCb3(2) = dB2(2)*csum + B2*( c31 *dB1a(2)*exp1 +c31 *B1*dexp1&
     &            +c32 *2.d0*B1*dB1a(2)*exp2 + c32 *B12*dexp2 )
	 dCb3(3) = dB2(3)*csum + B2*( c31 *dB1a(3)*exp1 +c31 *B1*dexp1&
     &            +c32 *2.d0*B1*dB1a(3)*exp2 + c32 *B12*dexp2 )
	 dCb4(1) = dB1a(1)*B3 *dsum1 + B1*dB3(1)*dsum1 + B1*B3 *ddsum1&
     &           + dB1a(1)*B3b*dsum2 + B1*dB3(4)*dsum2 + B1*B3b*ddsum2
	 dCb4(2) = dB1a(2)*B3 *dsum1 + B1*dB3(2)*dsum1 + B1*B3 *ddsum1&
     &           + dB1a(2)*B3b*dsum2 + B1*dB3(5)*dsum2 + B1*B3b*ddsum2
	 dCb4(3) = dB1a(3)*B3 *dsum1 + B1*dB3(3)*dsum1 + B1*B3 *ddsum1&
     &           + dB1a(3)*B3b*dsum2 + B1*dB3(6)*dsum2 + B1*B3b*ddsum2
	 dCb5(1) = dB1a(1)*fsum*exp7/P + B1*dfsum*exp7/P&
     &               + B1*fsum*dexp7/P - B1*fsum*exp7*dP(1)/Psq
	 dCb5(2) = dB1a(2)*fsum*exp7/P + B1*dfsum*exp7/P&
     &               + B1*fsum*dexp7/P - B1*Fsum*exp7*dP(2)/Psq
	 dCb5(3) = dB1a(3)*fsum*exp7/P + B1*dfsum*exp7/P&
     &               + B1*fsum*dexp7/P - B1*Fsum*exp7*dP(3)/Psq
	 dCb6(1) = dB1a(1)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
	 dCb6(2) = dB1a(2)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
	 dCb6(3) = dB1a(3)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
	 dCb7(1) = dB1a(1)*aasum*exp2 + B1*daasum*dP(1)*exp2&
     &                               + B1* aasum*dexp2
	 dCb7(2) = dB1a(2)*aasum*exp2 + B1*daasum*dP(2)*exp2&
     &                               + B1* aasum*dexp2
	 dCb7(3) = dB1a(3)*aasum*exp2 + B1*daasum*dP(3)*exp2&
     &                               + B1* aasum*dexp2
	 dCb8(1) = dB1a(1)*ffsum*exp7  +B1*dffsum*dP(1)*exp7 &
     &                                +B1* ffsum*dexp7 
	 dCb8(2) = dB1a(2)*ffsum*exp7  +B1*dffsum*dP(2)*exp7 &
     &                                +B1* ffsum*dexp7 
	 dCb8(3) = dB1a(3)*ffsum*exp7  +B1*dffsum*dP(3)*exp7 &
     &                                +B1* ffsum*dexp7 
           dCbnds(1,1) = dCb1(1) +dCb2(1) +dCb3(1) +dCb4(1) +dCb5(1)&
     &                  +dCb6(1) +dCb7(1) +dCb8(1) 
           dCbnds(1,2) = dCb1(2) +dCb2(2) +dCb3(2) +dCb4(2) +dCb5(2)&
     &                  +dCb6(2) +dCb7(2) +dCb8(2) 
           dCbnds(1,3) = dCb1(3) +dCb2(3) +dCb3(3) +dCb4(3) +dCb5(3)&
     &                  +dCb6(3) +dCb7(3) +dCb8(3) 
       B1    = B1b
       B12   = B1*B1
       B13   = B12*B1
       B14   = B13*B1
       B15   = B14*B1
       asum  = d11 + d12*R + d13*Rsq + d14/R + d15/Rsq 
       bsum  = d21*B12 + d22*B13 + d23*B14 + d24*B15
       csum  = d31*B1 *exp1  + d32*B12*exp2  
       dsum1 = d41*exp1 + d42*exp2
       dsum2 = d43*exp1 + d44*exp2
       fsum  = dx1 + d52*R + d53*Rsq 
       gsum  = d61 + d62/R + d63/Rsq 
       aasum = d71 + d72*P + d73*Psq + d74/P + d75/Psq 
       ffsum =       d81*P + d82*Psq + d84/Psq 
       dasum = d12 + 2.d0*d13*R -d14/Rsq -2.d0*d15/Rcu
       dbsum = 2.d0*d21*B1 +3.d0*d22*B12 +4.d0*d23*B13 +5.d0*d24*B14
       ddsum1= d41*dexp1 + d42*dexp2
       ddsum2= d43*dexp1 + d44*dexp2
       dfsum = d52 + 2.d0*d53*R
       dgsum = -d62/Rsq - 2.d0*d63/Rcu
       daasum = d72 +2.d0*d73*P -d74/Psq -2.d0*d75/Pcu
       dffsum = d81 + 2.d0*d82*P -2.d0*d84/Pcu
         cb(2,1)  = B1*asum * exp1  /P
         cb(2,2)  = bsum * exp2  
         cb(2,3)  = B2*csum 
         cb(2,4)  = B1*B3*dsum1 + B1*B3b*dsum2
         cb(2,5)  = B1 * fsum * exp7 /P
         cb(2,6)  = B1 * gsum * exp7  
         cb(2,7)  = B1 * aasum * exp2 
         cb(2,8)  = B1 * ffsum * exp7  
         Cbnds(2) = cb(2,1)+cb(2,2)+cb(2,3)+cb(2,4)+cb(2,5)&
     &                     +cb(2,6)+cb(2,7)+cb(2,8)
	 dCb1(1) = dB1b(1)*asum*exp1/P + B1*dasum*exp1/P&
     &                +B1*asum*dexp1/P - B1*asum*exp1*dP(1)/Psq
	 dCb1(2) = dB1b(2)*asum*exp1/P + B1*dasum*exp1/P&
     &                +B1*asum*dexp1/P - B1*asum*exp1*dP(2)/Psq
	 dCb1(3) = dB1b(3)*asum*exp1/P + B1*dasum*exp1/P&
     &                +B1*asum*dexp1/P - B1*asum*exp1*dP(3)/Psq
	 dCb2(1) = dbsum*dB1b(1)*exp2 + bsum*dexp2
	 dCb2(2) = dbsum*dB1b(2)*exp2 + bsum*dexp2
	 dCb2(3) = dbsum*dB1b(3)*exp2 + bsum*dexp2
	 dCb3(1)= dB2(1)*csum + B2*( d31 *dB1b(1)*exp1 +d31 *B1*dexp1&
     &           +d32 *2.d0*B1*dB1b(1)*exp2 + d32 *B12*dexp2 )
	 dCb3(2)= dB2(2)*csum + B2*( d31 *dB1b(2)*exp1 +d31 *B1*dexp1&
     &           +d32 *2.d0*B1*dB1b(2)*exp2 + d32 *B12*dexp2 )
	 dCb3(3)= dB2(3)*csum + B2*( d31 *dB1b(3)*exp1 +d31 *B1*dexp1&
     &           +d32 *2.d0*B1*dB1b(3)*exp2 + d32 *B12*dexp2 )
	 dCb4(1) = dB1b(1)*B3 *dsum1 + B1*dB3(1)*dsum1 + B1*B3 *ddsum1&
     &           + dB1b(1)*B3b*dsum2 + B1*dB3(4)*dsum2 + B1*B3b*ddsum2
	 dCb4(2) = dB1b(2)*B3 *dsum1 + B1*dB3(2)*dsum1 + B1*B3 *ddsum1&
     &           + dB1b(2)*B3b*dsum2 + B1*dB3(5)*dsum2 + B1*B3b*ddsum2
	 dCb4(3) = dB1b(3)*B3 *dsum1 + B1*dB3(3)*dsum1 + B1*B3 *ddsum1&
     &           + dB1b(3)*B3b*dsum2 + B1*dB3(6)*dsum2 + B1*B3b*ddsum2
	 dCb5(1) = dB1b(1)*fsum*exp7/P + B1*dfsum*exp7/P&
     &               + B1*fsum*dexp7/P - B1*fsum*exp7*dP(1)/Psq
	 dCb5(2) = dB1b(2)*fsum*exp7/P + B1*dfsum*exp7/P&
     &               + B1*fsum*dexp7/P - B1*fsum*exp7*dP(2)/Psq
	 dCb5(3) = dB1b(3)*fsum*exp7/P + B1*dfsum*exp7/P&
     &               + B1*fsum*dexp7/P - B1*fsum*exp7*dP(3)/Psq
	 dCb6(1) = dB1b(1)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
	 dCb6(2) = dB1b(2)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
	 dCb6(3) = dB1b(3)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
	 dCb7(1) = dB1b(1)*aasum*exp2 + B1*daasum*dP(1)*exp2&
     &                               + B1* aasum*dexp2
	 dCb7(2) = dB1b(2)*aasum*exp2 + B1*daasum*dP(2)*exp2&
     &                               + B1* aasum*dexp2
	 dCb7(3) = dB1b(3)*aasum*exp2 + B1*daasum*dP(3)*exp2&
     &                               + B1* aasum*dexp2
	 dCb8(1) = dB1b(1)*ffsum*exp7  +B1*dffsum*dP(1)*exp7 &
     &                                +B1* ffsum*dexp7 
	 dCb8(2) = dB1b(2)*ffsum*exp7  +B1*dffsum*dP(2)*exp7 &
     &                                +B1* ffsum*dexp7 
	 dCb8(3) = dB1b(3)*ffsum*exp7  +B1*dffsum*dP(3)*exp7 &
     &                                +B1* ffsum*dexp7 
           dCbnds(2,1) = dCb1(1) +dCb2(1) +dCb3(1) +dCb4(1) +dCb5(1)&
     &                  +dCb6(1) +dCb7(1) +dCb8(1) 
           dCbnds(2,2) = dCb1(2) +dCb2(2) +dCb3(2) +dCb4(2) +dCb5(2)&
     &                  +dCb6(2) +dCb7(2) +dCb8(2) 
           dCbnds(2,3) = dCb1(3) +dCb2(3) +dCb3(3) +dCb4(3) +dCb5(3)&
     &                  +dCb6(3) +dCb7(3) +dCb8(3) 
         dCbnda(1) = dT(1) * Cbnds(1)    + sumT * dCbnds(1,1)
         dCbndb(1) = dT(1) * Cbnds(2)    + sumT * dCbnds(2,1) 
         dCbnda(2) = dT(2) * Cbnds(1)    + sumT * dCbnds(1,2)
         dCbndb(2) = dT(2) * Cbnds(2)    + sumT * dCbnds(2,2) 
         dCbnda(3) = dT(3) * Cbnds(1)    + sumT * dCbnds(1,3)
         dCbndb(3) = dT(3) * Cbnds(2)    + sumT * dCbnds(2,3) 
      Cbnda = sumT*Cbnds(1)
      Cbndb = sumT*Cbnds(2)
      return
 2000 format(5x,'r1, r2, r3 = ',3(1x,f12.8))
 2010 format(5x,'B1a, B1b, B2, B3 = ',4(1x,g12.6))
 7100 format('    Vbnda = ',g12.6,'     Vbndb = ',g12.6)
 7650 format('    B2 =    ',g12.6,'       B3 =  ',g12.6)
 7700 format('    dB2 =   ',3(1x,g12.6),/,&
     &       '    dB3 =   ',6(1x,g12.6))
 7150 format(' k=',i1,'     derivatives: ',/,&
     &       ' dVb1  = ',3(1x,g16.10),/,&
     &       ' dVb2  = ',3(1x,g16.10),/,&
     &       ' dVb3  = ',3(1x,g16.10),/,&
     &       ' dVb4  = ',3(1x,g16.10),/,&
     &       ' dVb5  = ',3(1x,g16.10),/,&
     &       ' dVbnd = ',3(1x,g16.10))
 7155 format('  Vb1 = ',g16.10,'  Vb2 = ',g16.10,/,&
     &       '  Vb3 = ',g16.10,'  Vb4 = ',g16.10,/,&
     &       '  Vb5 = ',g16.10)
 2020 format(5x,'k=',i1,' cb1,cb2,cb3= ',3(1x,g16.10),&
     &     /,8x,' cb4,cb5,cb6= ',3(1x,g16.10),&
     &     /,8x,' cb7,cb8,cb9= ',3(1x,g16.10))
 2033 format(5x,'sumT = ',g12.6 ,'  cbnds(k) = ',g12.6)
 7152 format('Subr.Cbend:   k=',i1,'     derivatives: ',/,&
     &       '  dCb1  = ',3(1x,g18.12),/,&
     &       '  dCb2  = ',3(1x,g18.12),/,&
     &       '  dCb3  = ',3(1x,g18.12),/,&
     &       '  dCb4  = ',3(1x,g18.12),/,&
     &       '  dCb5  = ',3(1x,g18.12),/,&
     &       '  dCb6  = ',3(1x,g18.12),/,&
     &       '  dCb7  = ',3(1x,g18.12),/,&
     &       '  dCb8  = ',3(1x,g18.12),/,&
     &       '  dCb9  = ',3(1x,g18.12),/,&
     &       '  dCbnd = ',3(1x,g18.12))
 7200 format('  dT    = ',3(1x,g18.12),/,&
     &       '  dCbnda= ',3(1x,g18.12),/,&
     &       '  dCbndb= ',3(1x,g18.12))
 6000 format(' from subr.b1ab:  r=',3(1x,f9.6))
 6010 format('    b1a = ',f16.10,'     b1b = ',f16.10)
 6020 format('    db1a: ',3(1x,e12.6),/,'    db1b: ',3(1x,e12.6))
 7000 format('    dc(j)/dr(i) derivatives: ',/,&
     &    8x,3(1x,e12.6),/,8x,3(1x,e12.6),/,8x,3(1x,e12.6))
 7010 format('    ri2-rj2-rk2: ',3(1x,e12.6),/,&
     &       '       cosines : ',3(1x,e12.6),/,&
     &       '    12c**2 - 3 : ',3(1x,e12.6))
      end
      subroutine chgeom(r,ivalid)
      implicit double precision (a-h,o-z)
      parameter(derr = 1.d-5)
      parameter( istop1 = 1 , istop2 = 1 )
      dimension r(3)
      ivalid = 1
      rhi = dmax1(r(1),r(2))
      if(r(3).ge.rhi)then
         rmid = rhi
         rhi = r(3)
      else
         rmid = dmax1( r(3) , dmin1(r(1),r(2)) )
      endif
      rlo = dmin1(r(1),r(2),r(3))
      if(rlo+rmid+derr.lt.rhi)then
         write(6,*) 'Warning: invalid geometry ',r
         if(istop1.gt.0) stop ' STOP -- geometry error '
      end if
      if(rlo.lt.0.2d0)then
         write(6,*) 'Warning: invalid geometry ',r
         if(istop2.gt.0) stop ' STOP -- geometry error '
      end if
      return
      end
