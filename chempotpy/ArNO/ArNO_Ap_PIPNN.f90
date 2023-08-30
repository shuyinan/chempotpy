!======================================================================
!       1A' and 1A'' PESs of ArNO at the MRCI+Q level
!======================================================================

      subroutine pes(x,igrad,path,p,g,d)

      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      character(len=1024), intent(in) :: path
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: vtm, vim
      double precision :: r_NO, r_ArN, r_ArO, cthx
      integer :: iatom, idir, j, istate
      !initialize 
      vtm=0.0
      p=0.d0
      g=0.d0
      d=0.d0

      ! x: Ar, N, O order
      r_NO=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2+&
               &(x(2,3)-x(3,3))**2)
      r_NO=r_NO/0.529177211
      r_ArN=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2+&
               &(x(1,3)-x(2,3))**2)
      r_ArN=r_ArN/0.529177211
      r_ArO=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2+&
               &(x(1,3)-x(3,3))**2)
      r_ArO=r_ArO/0.529177211
      cthx=(r_ArN**2+r_NO**2-r_ArO**2)/2d0/r_ArN/r_NO

      if (igrad==0) then 
        call potread(path)
        call ArNOpes(r_NO,r_ArN,cthx,vtm,vim,1)
      else
        write (*,*) 'Only energy is available'
      endif

      vtm=vtm/8065.7112013

      do istate=1,nstates
        p(istate)=vtm
      enddo


      endsubroutine

        subroutine ArNOpes(r_NO,r_ArN,cthx,vtot,vint,istat)
        implicit real*8(a-h,o-z)
        call ArNOpes_sh(r_NO,r_ArN,cthx,v1,istat)
        call longpes(r_NO,r_ArN,cthx,v2,istat)
        S=0.5d0*(1d0-dtanh(0.8d0*(r_ArN-17d0)))
        vtot=S*v1+(1d0-S)*v2
        call pes1D(r_NO,v,istat)
        vint=vtot-v
        return
        end

        subroutine ArNOpes_sh(r10,r20,cthi,va,istate)
        implicit real*8(a-h,o-z)
        parameter (rbohr=0.5291771)
        parameter (pi=3.141592653589793d0)
        common /bvcut/vcut
! : r1=r(N-O), r2=r(Ar-N)
        r1=r10
        r2=r20
        cth=cthi
        if(cth.gt.1d0) cth=1d0
        if(cth.lt.-1d0) cth=-1d0
        th=dacos(cth)*180.0d0/pi
        if (r1.lt.1.60d0) r1=1.60d0
        if (r1.gt.4.0d0) r1=4.0d0
        if (r2.lt.3.50d0) r2=3.50d0
        if (r2.gt.40.d0) r2=40.d0
        if(istate.ne.1.and.istate.ne.2)istate=1
        CALL SPl3(r2,th,r1,va,istate)
!       set zero point at minimum of Ar + NO on A' state 
        va=va*219474.63067d0/27.2114D0-31.4930895583566D0
1002  format(1x,3f12.4,f20.8)
      end

        subroutine potread(path)
        implicit real*8(a-h,o-z)
        character(len=1024), intent(in) :: path
        character(len=1024) :: file_path1
        parameter (pi=3.141592653589793d0)    
        parameter (nroh1=17,nth=13,nroh2=70,ip=2)
        parameter (m=100,n=100,l=100)
        dimension xt(m),y(m),y2(l)  
        data vmin1/-656.852768990000d0/
       data dy1,dyn/1.0d30,1.0d30/
        common /bvcut/vcut
        common /pesa/pesmin,roh2(nroh1,nth,nroh2),roh1(nroh1),&
               &thth(nth),ve(ip,nroh1,nth,nroh2),ind(nroh1,nth)

        file_path1 = trim(path)//"/ArNO/spline.dat"
        open(98,file=file_path1,status='old')
!        open(66,file='check.txt')
        vcut=12d0 ! in eV
        pesmin=10000
!        read(98,*)
        do i=1,nroh1
          do j=1,nth
           ind(i,j)=0
            do k=1,nroh2
            read(98,*)ra1,th1,ra2,ve1,ve2
           if (ra1+ra2.lt.0.1) goto 7
            if (k.eq.1) then
            rra1=ra1
            rra2=ra2
            tth1=th1
            else
            dal=abs(ra1-rra1)+abs(th1-tth1)
            if (dal.gt.0.01) then
            write(*,*)' error  ',ii,ra1,rra1,ra2,rra2
            stop
            endif
            endif
            if (ve1.lt.pesmin) then
            pesmin=ve1
            r1m=ra1
            r2m=ra2
            thm=th1
            endif
           ind(i,j)=ind(i,j)+1
            roh2(i,j,k)=ra2
            roh1(i)=ra1
            thth(j)=th1
            ve(1,i,j,k)=(ve1-vmin1)*27.2114d0
            ve(2,i,j,k)=(ve2-vmin1)*27.2114d0
            enddo
7          continue
          enddo  
        enddo
!        write(*,*)r1m,r2m,thm,pesmin
!        write(*,*)' th1=',(thth(i),i=1,nth)
!        write(*,*)' roh1=',(roh1(i),i=1,nroh1)
          do i=1,nroh1
          do j=1,nth
           do k=1,ind(i,j)
          if (ve(1,i,j,k).gt.vcut) ve(1,i,j,k)=vcut
          if (ve(2,i,j,k).gt.vcut) ve(2,i,j,k)=vcut
!          write(66,1001)roh1(i),thth(j),roh2(i,j,k)
!     &,ve(1,i,j,k),ve(2,i,j,k)
 1001      format (3f8.3,6f12.8)
           enddo
         enddo
       enddo
        close(98)
        return
       end

       subroutine spl3(r1,th,r2,v,istate)
       implicit real*8(a-h,o-z)
       parameter (nroh1=17,nth=13,nroh2=70,ip=2)
        parameter (m=100,n=100,l=100)
        dimension xt(m)
       dimension dty(l),ddty(l),s1(l),ds1(l),dds1(l),h1(l)
       dimension dny(n),ddny(n),s2(l),ds2(l),dds2(l),h2(n)
       dimension dhy(m),ddhy(m),s3(l),ds3(l),dds3(l),h3(m)
       dimension y(m),ss(m),sss(m),y2(l)
       data dy1,dyn/1.0d30,1.0d30/
        common /bvcut/vcut
        common /pesa/pesmin,roh2(nroh1,nth,nroh2),roh1(nroh1),&
               &thth(nth),ve(ip,nroh1,nth,nroh2),ind(nroh1,nth)
                     
!        do iop=istate,istate
         iop=istate
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&  PRIVATE (i,j,nh,k,xt,y,y2,y3,ss,nthth,yw2)
!$OMP DO
        do 20 i=1,nroh1
       do 10 j=1,nth
        nh=0
       do 2 k=1,ind(i,j)
        nh=nh+1
        xt(nh)=roh2(i,j,k)
        y(nh)=ve(iop,i,j,k)
   2   continue
        r1a=r1
        call spline(xt,y,nh,dy1,dyn,y2)
        call splint(xt,y,y2,nh,r1a,y3)
 22     continue
        if (y3.gt.vcut) y3=vcut
        ss(j)=y3
   10   continue
        nthth=0
       do 5 j=1,nth
        nthth=nthth+1
        xt(nthth)=thth(j)
        y(nthth)=ss(j)
   5   continue
        call spline(xt,y,nthth,0.0d0,0.0d0,y2)
        call splint(xt,y,y2,nthth,th,yw2)
 33     continue
        if(yw2.gt.vcut) yw2=vcut
       sss(i)=yw2
   20    continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        nr1=nroh1
!       write(*,*)' vrhh=',(sss(i),i=1,noh)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&  PRIVATE (i)
!$OMP DO
        do i=1,nroh1
!       nr1=nr1+1
        xt(i)=roh1(i)
        y(i)=sss(i)
        enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        call spline(xt,y,nr1,dy1,dyn,y2)
        call splint(xt,y,y2,nr1,r2,yw2)
   44   continue
        if (yw2.gt.vcut) yw2=vcut
        v=yw2
!       enddo
       end

!##################################################################
!# SPLINE ROUTINES
!#            Numerical recipes in fortran
!#            Cambrige University Press
!#            York, 2nd edition, 1992.
!##################################################################
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      implicit double precision  (a-h,o-z)
      DIMENSION xa(n),y2a(n),ya(n)
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.0d0) write(6,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*&
        &(h**2)/6.0d0
      return
      END
!##############################################################################
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit double precision  (a-h,o-z)
      DIMENSION x(n),y(n),y2(n)
      PARAMETER (NMAX=100)
      DIMENSION u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0d0
        y2(i)=(sig-1.0d0)/p
        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+&
          &1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-&
          &sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.0d0
        un=0.0d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

        subroutine longpes(r_NO,r_ArN,cthx,v,istat)
        implicit real*8(a-h,o-z)
        dimension x(17),a1(4,17),a2(4,17),y(17),y2(17),de(17)
        data vmin1/-656.852768990000d0/
!       R in angstrom        
        amN=14.0030740d0
        amO=15.9949146d0
        r2=r_NO
        r1=r_ArN
        t=amO/(amO+amN)*r2
        s=dsqrt(r1**2+t**2-2d0*r1*t*cthx) ! in bohr
        R=s*0.5291772083d0
        cth=(t**2+s**2-r1**2)/2d0/t/s
        cth=min(1d0,cth);cth=max(-1d0,cth)

        tocm=219474.63067d0
        dy1=1d30
        dyn=1d30
        if(r2.lt.1.6d0)r2=1.6d0
        if(r2.gt.4d0)r2=4d0
        x(1)=1.60000000000000d0
        x(2)=1.80000000000000d0
        x(3)=1.90000000000000d0
        x(4)=1.95000000000000d0
        x(5)=2.00000000000000d0
        x(6)=2.05000000000000d0
        x(7)=2.10000000000000d0
        x(8)=2.15000000000000d0
        x(9)=2.20000000000000d0
        x(10)=2.2500000000000d0 
        x(11)=2.3000000000000d0 
        x(12)=2.4000000000000d0 
        x(13)=2.6000000000000d0 
        x(14)=2.8000000000000d0 
        x(15)=3.0000000000000d0 
        x(16)=3.5000000000000d0 
        x(17)=4.0000000000000d0 
        a1(1,1)=-267943.38590546D0   
        a1(2,1)=-38055.810825624D0   
        a1(3,1)=-22967.743334597D0   
        a1(4,1)=11194.869520824D0    
        a1(1,2)=-251809.17227766D0   
        a1(2,2)=-27650.194725372D0   
        a1(3,2)=-36839.299113681D0   
        a1(4,2)=8549.8153161124D0    
        a1(1,3)=-251879.04153352D0   
        a1(2,3)=-23032.822374801D0   
        a1(3,3)=-41744.334275868D0   
        a1(4,3)=6985.8690896270D0    
        a1(1,4)=-250283.76005970D0   
        a1(2,4)=-20018.193520338D0   
        a1(3,4)=-44240.758806033D0   
        a1(4,4)=5962.8283874165D0    
        a1(1,5)=-254138.28642784D0   
        a1(2,5)=-19723.912012097D0   
        a1(3,5)=-46256.683875458D0   
        a1(4,5)=6261.8163891098D0    
        a1(1,6)=-255157.48758001D0   
        a1(2,6)=-18048.919914254D0   
        a1(3,6)=-49265.444654343D0   
        a1(4,6)=5298.6955381574D0    
        a1(1,7)=-255804.28482087D0   
        a1(2,7)=-15646.998402989D0   
        a1(3,7)=-50782.220878825D0   
        a1(4,7)=5032.7896739693D0    
        a1(1,8)=-257992.33687205D0   
        a1(2,8)=-14507.636299835D0   
        a1(3,8)=-53263.086294840D0   
        a1(4,8)=5668.8453773915D0    
        a1(1,9)=-262083.75815139D0   
        a1(2,9)=-14037.602789474D0   
        a1(3,9)=-55936.840259583D0   
        a1(4,9)=5142.3562977947D0    
        a1(1,10)=-264567.17397801D0  
        a1(2,10)=-13316.118074956D0  
        a1(3,10)=-58104.973869090D0  
        a1(4,10)=5459.9310221469D0   
        a1(1,11)=-267843.93076864D0  
        a1(2,11)=-12807.022293995D0  
        a1(3,11)=-60824.874069809D0  
        a1(4,11)=5830.2661402447D0   
        a1(1,12)=-276751.23192033D0  
        a1(2,12)=-13497.471138696D0  
        a1(3,12)=-66418.984554633D0  
        a1(4,12)=5439.6524419885D0   
        a1(1,13)=-290143.29086576D0  
        a1(2,13)=-11329.693163825D0  
        a1(3,13)=-77691.580833384D0  
        a1(4,13)=5266.1638029908D0   
        a1(1,14)=-302025.84562829D0  
        a1(2,14)=-9660.7649810443D0  
        a1(3,14)=-88948.228567055D0  
        a1(4,14)=6062.8633344077D0   
        a1(1,15)=-311769.98048147D0  
        a1(2,15)=-7392.2131839275D0  
        a1(3,15)=-101015.44505332D0  
        a1(4,15)=5729.7740086266D0   
        a1(1,16)=-322755.44544985D0  
        a1(2,16)=-7292.9293756397D0  
        a1(3,16)=-115151.78247626D0  
        a1(4,16)=9316.8888253736D0   
        a1(1,17)=-314058.98722717D0  
        a1(2,17)=-9845.6172327698D0  
        a1(3,17)=-108095.25443221D0  
        a1(4,17)=12298.028717238D0   
        a2(1,1)=-240431.54496463D0   
        a2(2,1)=-24089.327542574D0   
        a2(3,1)=-47953.116074105D0   
        a2(4,1)=-1821.4704857423D0   
        a2(1,2)=-241713.62149381D0   
        a2(2,2)=-19568.897780005D0   
        a2(3,2)=-45838.247416406D0   
        a2(4,2)=1062.8908146373D0    
        a2(1,3)=-243933.09613888D0   
        a2(2,3)=-17425.641520976D0   
        a2(3,3)=-47462.042744971D0   
        a2(4,3)=1884.9096687051D0    
        a2(1,4)=-244670.78972127D0   
        a2(2,4)=-15757.181068769D0   
        a2(3,4)=-48668.574025652D0   
        a2(4,4)=2025.3375037173D0    
        a2(1,5)=-249752.59229822D0   
        a2(2,5)=-16345.228476962D0   
        a2(3,5)=-49790.778143198D0   
        a2(4,5)=3270.3444683698D0    
        a2(1,6)=-251769.89335023D0   
        a2(2,6)=-15313.617119895D0   
        a2(3,6)=-52458.348815466D0   
        a2(4,6)=2325.1240929174D0    
        a2(1,7)=-251678.55830509D0   
        a2(2,7)=-12958.861408430D0   
        a2(3,7)=-53659.913259469D0   
        a2(4,7)=2575.1831916109D0    
        a2(1,8)=-253437.72857935D0   
        a2(2,8)=-12258.648248398D0   
        a2(3,8)=-55795.170116181D0   
        a2(4,8)=3732.0847356143D0    
        a2(1,9)=-258471.67365255D0   
        a2(2,9)=-12799.769627073D0   
        a2(3,9)=-57508.633013489D0   
        a2(4,9)=4303.6465924945D0    
        a2(1,10)=-262187.31147592D0  
        a2(2,10)=-12307.330100049D0  
        a2(3,10)=-59805.051173756D0  
        a2(4,10)=4235.9341267133D0   
        a2(1,11)=-265266.31321490D0  
        a2(2,11)=-11337.914954709D0  
        a2(3,11)=-62140.255879242D0  
        a2(4,11)=5747.4771478058D0   
        a2(1,12)=-277136.94531031D0  
        a2(2,12)=-13950.852000624D0  
        a2(3,12)=-67638.158750528D0  
        a2(4,12)=5348.7496751312D0   
        a2(1,13)=-287623.13431309D0  
        a2(2,13)=-10926.585393021D0  
        a2(3,13)=-78513.457663825D0  
        a2(4,13)=5665.0625303400D0   
        a2(1,14)=-300470.77056637D0  
        a2(2,14)=-10421.152677201D0  
        a2(3,14)=-89919.062773954D0  
        a2(4,14)=7178.0270447738D0   
        a2(1,15)=-312097.32976855D0  
        a2(2,15)=-10237.054620313D0  
        a2(3,15)=-102418.68339132D0  
        a2(4,15)=7772.4692422368D0   
        a2(1,16)=-322700.43808882D0  
        a2(2,16)=-13197.018758343D0  
        a2(3,16)=-118420.01108172D0  
        a2(4,16)=13070.129322234D0   
        a2(1,17)=-307994.98282431D0  
        a2(2,17)=-12733.943399831D0  
        a2(3,17)=-111176.93493823D0  
        a2(4,17)=16996.093766663D0   

        if(istat.ne.1.and.istat.ne.2)istat=1
        if(istat.eq.1)then
        de(1)= -656.45779293d0 
        de(2)= -656.72825267d0 
        de(3)= -656.79453152d0 
        de(4)= -656.81626228d0 
        de(5)= -656.83198474d0 
        de(6)= -656.84267992d0 
        de(7)= -656.84918525d0 
        de(8)= -656.85221531d0 
        de(9)= -656.85237944d0 
        de(10)= -656.85019678d0
        de(11)= -656.84610901d0
        de(12)= -656.83366196d0
        de(13)= -656.79902012d0
        de(14)= -656.76032202d0
        de(15)= -656.72354209d0
        de(16)= -656.65642703d0
        de(17)= -656.62721422d0
        do i=1,17
        call VFUNCS(cth,r,a1(:,i),v)
        y(i)=v+(de(i)-vmin1)*tocm
        enddo
        else
        de(1)= -656.45778154d0
        de(2)= -656.72823764d0
        de(3)= -656.79451533d0
        de(4)= -656.81624562d0
        de(5)= -656.83196768d0
        de(6)= -656.84266252d0
        de(7)= -656.84916757d0
        de(8)= -656.85219739d0
        de(9)= -656.85236132d0
        de(10)=-656.85017849d0
        de(11)=-656.84609058d0
        de(12)=-656.83364329d0
        de(13)=-656.79900106d0
        de(14)=-656.76030267d0
        de(15)=-656.72352380d0
        de(16)=-656.65639961d0
        de(17)=-656.62718628d0
        do i=1,17
        call VFUNCS(cth,r,a2(:,i),v)
        y(i)=v+(de(i)-vmin1)*tocm
        enddo
        endif
        call spline(x,y,17,dy1,dyn,y2)
        call splint(x,y,y2,17,r2,v)
!       set zero point at minimum of Ar + NO on A' state 
        v=v-31.4930895583566D0
        return
        end
!***********************************************************************
       subroutine VFUNCS(theta,r,c,v)
      implicit real*8(a-h,o-z)
      dimension afunc(4),c(4)

        cth=theta
        call Fval(R,cth,F)
        afunc(1)=F/R**6
        afunc(2)=cth*F/R**6
        afunc(3)=0.5d0*(3d0*cth**2-1d0)*F/R**6
        afunc(4)=0.5d0*(5d0*cth**3-3d0*cth)*F/R**6
        v=0d0
        do i=1,4
        v=v+afunc(i)*c(i)
        enddo
        return
      end

        subroutine Bval(cth,B)
        implicit none
        real*8 cth,B,b0,b1,b2
        b0=-3.0318d0
        b1=0.1299d0
        b2=0.0564d0
        B=b0+b1*cth+b2*0.5d0*(3d0*cth**2-1d0)
        return
        end
        subroutine Fval(R,cth,F)
        implicit none
        real*8 R,cth,F,B

        call Bval(cth,B)
        B=abs(B)*R
        F=1d0-dexp(-B)*(1d0+B+B**2/2d0+B**3/6d0+B**4/24d0+B**5/120d0+&
          &B**6/720d0)
        return
        end

        subroutine pes1D(r_NO,v,istat)
        implicit real*8(a-h,o-z)
        dimension x(17),a1(4,17),a2(4,17),y(17),y2(17),de(17)
        data vmin1/-656.852768990000d0/
        r2=r_NO
        tocm=219474.63067d0
        dy1=1d30
        dyn=1d30
        if(r2.lt.1.6d0)r2=1.6d0
        if(r2.gt.4d0)r2=4d0
        x(1)=1.60000000000000d0
        x(2)=1.80000000000000d0
        x(3)=1.90000000000000d0
        x(4)=1.95000000000000d0
        x(5)=2.00000000000000d0
        x(6)=2.05000000000000d0
        x(7)=2.10000000000000d0
        x(8)=2.15000000000000d0
        x(9)=2.20000000000000d0
        x(10)=2.2500000000000d0 
        x(11)=2.3000000000000d0 
        x(12)=2.4000000000000d0 
        x(13)=2.6000000000000d0 
        x(14)=2.8000000000000d0 
        x(15)=3.0000000000000d0 
        x(16)=3.5000000000000d0 
        x(17)=4.0000000000000d0 
        if(istat.ne.1.and.istat.ne.2)istat=1
        if(istat.eq.1)then
        de(1)= -656.45779293d0 
        de(2)= -656.72825267d0 
        de(3)= -656.79453152d0 
        de(4)= -656.81626228d0 
        de(5)= -656.83198474d0 
        de(6)= -656.84267992d0 
        de(7)= -656.84918525d0 
        de(8)= -656.85221531d0 
        de(9)= -656.85237944d0 
        de(10)= -656.85019678d0
        de(11)= -656.84610901d0
        de(12)= -656.83366196d0
        de(13)= -656.79902012d0
        de(14)= -656.76032202d0
        de(15)= -656.72354209d0
        de(16)= -656.65642703d0
        de(17)= -656.62721422d0
        do i=1,17
        y(i)=(de(i)-vmin1)*tocm
        enddo
        else
        de(1)= -656.45778154d0
        de(2)= -656.72823764d0
        de(3)= -656.79451533d0
        de(4)= -656.81624562d0
        de(5)= -656.83196768d0
        de(6)= -656.84266252d0
        de(7)= -656.84916757d0
        de(8)= -656.85219739d0
        de(9)= -656.85236132d0
        de(10)=-656.85017849d0
        de(11)=-656.84609058d0
        de(12)=-656.83364329d0
        de(13)=-656.79900106d0
        de(14)=-656.76030267d0
        de(15)=-656.72352380d0
        de(16)=-656.65639961d0
        de(17)=-656.62718628d0
        do i=1,17
        y(i)=(de(i)-vmin1)*tocm
        enddo
        endif
        call spline(x,y,17,dy1,dyn,y2)
        call splint(x,y,y2,17,r2,v)
!       set zero point at minimum of Ar + NO on A' state 
        v=v-31.4930895583566D0
        return
        end
