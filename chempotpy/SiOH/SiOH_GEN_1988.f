      subroutine pes(xyz,nSi,nO,nH,natoms,igrad,p,g,d)
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, intent(in) :: nSi, nO, nH
      integer, intent(in) :: natoms
      double precision, intent(in) :: xyz(10000,3)
      integer, intent(in) :: igrad
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: cx(10000), cy(10000), cz(10000)
      double precision :: dcx(10000), dcy(10000), dcz(10000)
      double precision :: q(10000)
      integer :: indexlist(10000)
      double precision :: v
      integer :: istate, iatom

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      q=0.d0
      do iatom=1,natoms
         cx(iatom)=xyz(iatom,1)/0.529177211
         cy(iatom)=xyz(iatom,2)/0.529177211
         cz(iatom)=xyz(iatom,3)/0.529177211
      enddo
      do iatom=1,nSi
        indexlist(iatom)=14
      enddo
      do iatom=1,nO
        indexlist(iatom+nSi)=8
      enddo
      do iatom=1,nH
        indexlist(iatom+nSi+nO)=1
      enddo 

      call siohpot(cx,cy,cz,q,indexlist,v,dcx,dcy,dcz,natoms,10000,0)

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=v*27.211386
        enddo
      else if (igrad==1) then
        do istate=1,nstates
          p(istate)=v*27.211386
        enddo
        do istate=1,nstates
        do iatom=1,natoms
          g(istate,iatom,1)=dcx(iatom)*51.422067
          g(istate,iatom,2)=dcy(iatom)*51.422067
          g(istate,iatom,3)=dcz(iatom)*51.422067
        enddo
        enddo
      else
        write (*,*) 'Only energy and gradient are available'
      endif

      endsubroutine

       subroutine siohpot(x,y,z,q,indmm,v,dvdx,dvdy,dvdz,natom,
     &  maxatom,ien)

c   System:                        SiOH
c   Functional form:               Two and three body interaction plus Couloumb interaction
c   Common name:                   SiOH
c   Number of derivatives:         1
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                             
c
c   Notes:  Many-body silica  potential energy function.  The functional 
c           form is from Ref. 1.  
c
c  References: (1) JApplPhys89(11) 2001 6013
c              (2) JChemPhys 89(9) 1988 5818
c
c  Units:
c       energies    - hartrees
c       coordinates - bohrs
c
c  v       --- energy of the system (output)
c  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
c              of the system (input)
c  dvdx,dvdy,dvdz   --- one-dimensional arrays of the gradients with respect
c                       to the Cartesian coordinates of the system (output)
c  natom   --- number of atoms in the system (input)
c  maxatom --- maximum number of atoms (input)
c  ien     --- 1: just take coulombic interaction energy
C          --- 0: take all energy

      implicit double precision(a-h,o-z)
c     include 'param.f'
      parameter(rn=10.393496649d0)
      parameter(pi=3.1415926536d0)
      parameter(pener=4.3597482d-11)
      parameter(pdis=0.52917706d0)
      parameter(pang=57.2957795d0)
      dimension alpha(14,14),beta(14,14),rho(14,14),q(maxatom)
       dimension aa(14,14,3),bb(14,14,3),cc(14,14,3)
      dimension alambda(14,14,14),gamma(22,36),r0(22,36),
     &          theta(14,14,14)
      dimension dx(maxatom,maxatom),rij(maxatom,maxatom),
     &   dy(maxatom,maxatom),dz(maxatom,maxatom),indmm(maxatom)
      double precision dvdy(maxatom),dvdz(maxatom),dvdx(maxatom),
     &   x(maxatom),y(maxatom),z(maxatom)

c Nanoparticle-parameterized Parameters in two-body term(Ref. 1)
        alpha(14,14) =   0.18770d-8/pener
        alpha(14,8)  =   0.29620d-8/pener
        alpha(8,14)  =   alpha(14,8)
        alpha(14,1)  =   0.00690d-8/pener
        alpha(1,14)  =   alpha(14,1)
        alpha(8,8)   =   0.07250d-8/pener
        alpha(8,1)   =   0.03984d-8/pener
        alpha(1,8)   =   alpha(8,1)
        alpha(1,1)   =   0.00340d-8/pener
        beta(14,14)  =   2.29d0/pdis
        beta(14,8)   =   2.34d0/pdis
        beta(8,14)   =   beta(14,8)
        beta(14,1)   =   2.31d0/pdis
        beta(1,14)   =   beta(14,1)
        beta(8,8)    =   2.34d0/pdis
        beta(8,1)    =   2.26d0/pdis
        beta(1,8)    =   beta(8,1)
        beta(1,1)    =   2.10d0/pdis
        rho(14,14)   =   0.29d0/pdis
        rho(14,8)    =   0.29d0/pdis
        rho(8,14)    =   rho(14,8)
        rho(14,1)    =   0.29d0/pdis
        rho(1,14)    =   rho(14,1)
        rho(8,8)     =   0.29d0/pdis
        rho(8,1)     =   0.29d0/pdis
        rho(1,8)     =   rho(8,1)
        rho(1,1)     =   0.35d0/pdis
C        q(14)    =   4.0d0
C        q(8)     =  -2.0d0
C        q(1)     =   1.0d0
        aa(14,1,1)   =  -4.6542d-12/pener
        aa(1,14,1)   =  aa(14,1,1)
        aa(14,1,2)   =   0.d0
        aa(1,14,2)   =  aa(14,1,2)
        aa(14,1,3)   =   0.d0
        aa(1,14,3)   =  aa(14,1,3)
        aa(1,1,1)    =  -5.2793d-12/pener
        aa(1,1,2)    =   0.3473d-12/pener
        aa(1,1,3)    =   0.d0
        aa(8,1,1)    =  -2.0840d-12/pener
        aa(8,1,2)    =   7.6412d-12/pener
        aa(8,1,3)    =  -0.8336d-12/pener
        aa(1,8,1)    =  aa(8,1,1)
        aa(1,8,2)    =  aa(8,1,2)
        aa(1,8,3)    =  aa(8,1,3)
        bb(14,1,1)   =   6.0d0*pdis
        bb(14,1,2)   =   0.d0
        bb(14,1,3)   =   0.d0
        bb(1,14,1)   =  bb(14,1,1)
        bb(1,14,2)   =  bb(14,1,2)
        bb(1,14,3)   =  bb(14,1,3)
        bb(1,1,1)    =   6.0d0*pdis
        bb(1,1,2)    =   2.0d0*pdis
        bb(1,1,3)    =   0.d0
        bb(8,1,1)    =  15.0d0*pdis
        bb(8,1,2)    =   3.2d0*pdis
        bb(8,1,3)    =   5.0d0*pdis
        bb(1,8,1)    =  bb(8,1,1)
        bb(1,8,2)    =  bb(8,1,2)
        bb(1,8,3)    =  bb(8,1,3)
        cc(14,1,1)   =   2.20d0/pdis
        cc(14,1,2)   =   0.d0
        cc(14,1,3)   =   0.d0
        cc(1,14,1)   =  cc(14,1,1)
        cc(1,14,2)   =  cc(14,1,2)
        cc(1,14,3)   =  cc(14,1,3)
        cc(1,1,1)    =   1.51d0/pdis
        cc(1,1,2)    =   2.42d0/pdis
        cc(1,1,3)    =   0.d0
        cc(8,1,1)    =   1.05d0/pdis
        cc(8,1,2)    =   1.50d0/pdis
        cc(8,1,3)    =   2.00d0/pdis
        cc(1,8,1)    =  cc(8,1,1)
        cc(1,8,2)    =  cc(8,1,2)
        cc(1,8,3)    =  cc(8,1,3)
c parameters in the three-body term
        alambda(8,14,8)   =  19.0d-11/pener
        alambda(14,8,14)  =   0.3d-11/pener
        alambda(14,8,1)   =   5.0d-11/pener
        alambda(1,8,14)   =   alambda(14,8,1)
        alambda(1,8,1)    =  35.0d-11/pener
        gamma(22,30)      =   2.8d0/pdis
        gamma(22,36)      =   2.0d0/pdis
        gamma(22,23)      =   2.0d0/pdis
        gamma(9,23)       =   1.2d0/pdis
        gamma(9,10)       =   1.3d0/pdis
        r0(22,30)         =   3.0d0/pdis
        r0(22,36)         =   2.6d0/pdis
        r0(22,23)         =   2.6d0/pdis
        r0(9,10)          =   1.6d0/pdis
        r0(9,23)          =   1.5d0/pdis
        theta(8,14,8)     = 109.5d0/pang
        theta(14,8,14)    = 109.5d0/pang
        theta(14,8,1)     = 109.5d0/pang
        theta(1,8,14)     = theta(14,8,1)
        theta(1,8,1)      = 104.5d0/pang

 	  	  	  	  	   
c Initialize
C         do i=1,natom
C           write(6,100) indmm(i),q(i),x(i)*pdis,y(i)*pdis,z(i)*pdis
C         enddo

      v1 = 0.d0
      v2 = 0.d0
      v3 = 0.d0
      v  = 0.d0
      do i=1,natom
      dvdx(i)  = 0.d0
      dvdy(i)  = 0.d0
      dvdz(i)  = 0.d0
      do j=i+1,natom
      dx(i,j)=x(i)-x(j)
      dy(i,j)=y(i)-y(j)
      dz(i,j)=z(i)-z(j)
      rij(i,j)=dsqrt(dx(i,j)*dx(i,j)+dy(i,j)*dy(i,j)+dz(i,j)*dz(i,j))
      rij(j,i)=rij(i,j)
      dx(j,i)=-dx(i,j)
      dy(j,i)=-dy(i,j)
      dz(j,i)=-dz(i,j)
      enddo
      enddo

c Double loop (main loop)
      do 10 i=1,natom
      do 20 j=1,natom

c Cutoff function
      if (i.eq.j) goto 20
      if (rij(i,j).le.rn) then
      if (i.gt.j) goto 111 
         v1=v1+alpha(indmm(i),indmm(j))*dexp(-rij(i,j)/rho(indmm(i),
     &           indmm(j)))
     &	   +(q(i)*q(j)/rij(i,j))*Erfc(rij(i,j)/
     &      beta(indmm(i),indmm(j)))
          e=-alpha(indmm(i),indmm(j))*dexp(-rij(i,j)/rho(indmm(i),
     &      indmm(j)))/rho(indmm(i),indmm(j))
     & -q(i)*q(j)*Erfc(rij(i,j)/beta(indmm(i),indmm(j)))/
     &     rij(i,j)**2
     &    -2*q(i)*q(j)*dexp(-(rij(i,j)/beta(indmm(i),
     &     indmm(j)))**2)/(dsqrt(pi)*rij(i,j)*beta(indmm(i),indmm(j))) 
          tmp=e/rij(i,j)
          dvdx(i) = dvdx(i)+tmp*dx(i,j)
          dvdy(i) = dvdy(i)+tmp*dy(i,j)
          dvdz(i) = dvdz(i)+tmp*dz(i,j)
          dvdx(j) = dvdx(j)-tmp*dx(i,j)
          dvdy(j) = dvdy(j)-tmp*dy(i,j)
          dvdz(j) = dvdz(j)-tmp*dz(i,j)
          e1=dexp(bb(indmm(i),indmm(j),1)*(rij(i,j)-cc(indmm(i),
     &       indmm(j),1)))
          e2=dexp(bb(indmm(i),indmm(j),2)*(rij(i,j)-cc(indmm(i),
     &       indmm(j),2)))
          e3=dexp(bb(indmm(i),indmm(j),3)*(rij(i,j)-cc(indmm(i),
     &       indmm(j),3)))
         v2=v2+aa(indmm(i),indmm(j),1)/(1+e1)+aa(indmm(i),indmm(j),2)
     &       /(1+e2)+aa(indmm(i),indmm(j),3)/(1+e3)
         f=-aa(indmm(i),indmm(j),1)*bb(indmm(i),indmm(j),1)*e1/(1+e1)**2
     &     -aa(indmm(i),indmm(j),2)*bb(indmm(i),indmm(j),2)*e2/(1+e2)**2
     &     -aa(indmm(i),indmm(j),3)*bb(indmm(i),indmm(j),3)*e3/(1+e3)**2
          tmp=f/rij(i,j)
         dvdx(i) = dvdx(i)+tmp*dx(i,j)
         dvdy(i) = dvdy(i)+tmp*dy(i,j)
         dvdz(i) = dvdz(i)+tmp*dz(i,j)
         dvdx(j) = dvdx(j)-tmp*dx(i,j)
         dvdy(j) = dvdy(j)-tmp*dy(i,j)
         dvdz(j) = dvdz(j)-tmp*dz(i,j)
 111     continue        
            v3=0.d0
       do 30 k=j+1,natom
        if (k.eq.i .or. k.eq.j) goto 30
        if (indmm(i).eq.1) goto 30
        if (indmm(i).eq.8 .and. (indmm(j).eq.8 .or. indmm(k).eq.8)) 
     &     goto 30
        if (indmm(i).eq.14 .and. (indmm(j).ne.8 .or. indmm(k).ne.8))
     &      goto 30
        if (rij(i,j).ge.r0(indmm(i)+indmm(j),indmm(i)+indmm(j)+indmm(k))
     &       .or. rij(i,k).ge.r0(indmm(i)+indmm(k),indmm(i)+indmm(k)+
     &       indmm(j))) then
            v3=0.d0
        else
         g=dexp(gamma(indmm(i)+indmm(j),indmm(i)+indmm(j)+indmm(k))/
     &     (rij(i,j)-r0(indmm(i)+indmm(j),indmm(i)+indmm(j)+indmm(k)))
     &     +gamma(indmm(i)+indmm(k),indmm(i)+indmm(k)+indmm(j))/
     &     (rij(i,k)-r0(indmm(i)+indmm(k),indmm(i)+indmm(k)+indmm(j))))
         h  =(rij(i,j)**2+rij(i,k)**2-rij(j,k)**2)/(2*rij(i,j)*rij(i,k))
         p  =h-dcos(theta(indmm(j),indmm(i),indmm(k)))
         v3 =alambda(indmm(j),indmm(i),indmm(k))*g*p**2
         f3 =-alambda(indmm(j),indmm(i),indmm(k))
     &   *g*gamma(indmm(i)+indmm(j),indmm(i)+indmm(j)+indmm(k))*p**2
     &   /(rij(i,j)-r0(indmm(i)+indmm(j),indmm(i)+indmm(j)+indmm(k)))**2
     &   +2*alambda(indmm(j),indmm(i),indmm(k))*g*p*(rij(i,j)**2
     &           -rij(i,k)**2+rij(j,k)**2)/(2*rij(i,k)*rij(i,j)**2)
         f4=-alambda(indmm(j),indmm(i),indmm(k))
     &   *g*gamma(indmm(i)+indmm(k),indmm(i)+indmm(k)+indmm(j))*p**2
     &   /(rij(i,k)-r0(indmm(i)+indmm(k),indmm(i)+indmm(k)+indmm(j)))**2
     &    +2*alambda(indmm(j),indmm(i),indmm(k))*g*p*(rij(i,k)**2
     &    -rij(i,j)**2+rij(j,k)**2)/(2*rij(i,j)*rij(i,k)**2)
         f5=-2*alambda(indmm(j),indmm(i),indmm(k))*g*p*rij(j,k)
     &          /(rij(i,j)*rij(i,k))
c
c
          tmp=f3/rij(i,j)
         dvdx(i)=dvdx(i)+tmp*dx(i,j)
         dvdy(i)=dvdy(i)+tmp*dy(i,j)
         dvdz(i)=dvdz(i)+tmp*dz(i,j)
         dvdx(j)=dvdx(j)-tmp*dx(i,j)
         dvdy(j)=dvdy(j)-tmp*dy(i,j)
         dvdz(j)=dvdz(j)-tmp*dz(i,j)
          tmp=f4/rij(i,k)
         dvdx(i)=dvdx(i)+tmp*dx(i,k)
         dvdy(i)=dvdy(i)+tmp*dy(i,k)
         dvdz(i)=dvdz(i)+tmp*dz(i,k)
         dvdx(k)=dvdx(k)-tmp*dx(i,k)
         dvdy(k)=dvdy(k)-tmp*dy(i,k)
         dvdz(k)=dvdz(k)-tmp*dz(i,k)
          tmp=f5/rij(j,k)
         dvdx(j)=dvdx(j)+tmp*dx(j,k)
         dvdy(j)=dvdy(j)+tmp*dy(j,k)
         dvdz(j)=dvdz(j)+tmp*dz(j,k)
         dvdx(k)=dvdx(k)-tmp*dx(j,k)
         dvdy(k)=dvdy(k)-tmp*dy(j,k)
         dvdz(k)=dvdz(k)-tmp*dz(j,k)
        endif
	v=v+v3
 30   enddo
      endif
C        write(6,*) "i,j,v2"
C        write(6,*) i,j,v2*autoev
 20   enddo
 10   enddo
      if (ien.eq.1) then
        v = v1
      else
        v=v+v1+v2
      endif
C
 100    format(1x,i5,f10.2,3f10.5)
C      do i=1,natom
C      write(6,'(1x,i5,3f12.6)') i,dvdx(i),dvdy(i),dvdz(i)
C      enddo
       end 
