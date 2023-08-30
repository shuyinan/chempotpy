!   System:                     CH3NH2
!   Functional form:            Anchor-Points Reactive Potential (APRP) 
!                               (2x2 diabatic potential energy matrix)
!   Common name:                Methylamine(Coupled PESs)
!   Number of derivatives:      1
!   Number of bodies:           7
!   Number of electronic surfaces: 2
!   Interface: Section-2
!
!   Reference: Kelsey A. Parker,and Donald G. Truhlar 
!              J. Chem. Phys. 152, 244309 (2020).
!
!   Notes:
!    1. This routine calculates two coupled singlet PESs for
!       nonadiabatic NH photodissociation of CH3NH2 to produce
!       CH3NH + H. It calculates the diabatic potential energy matrix
!       and its gradient and the adiabatic potential energy surfaces
!       and their gradients, and it calculates the space-frame 
!       nonadiabatic couling vector.
!    2. Only the first H is allowed to dissociate.
!    3. The final version of this potential was added to POTLIB 
!       in March, 2020.
!
!   Numbering of atoms:
!       H5         H1
!         \       /
!       H6--C4---N3
!         /       \
!       H7         H2
!
!   Units:      input:   angstrom for coordinates
!               output:  hartree for energy
!                        hartree/bohr for gradients
!                        1/bohr for nonadiabatic coupling vector
!
!   Input:  
!   igrad = 0           Energy only calculation
!         = 1           Energy + Analytic gradient
!   repflag             flag to indicate wheter to use diabatic or 
!                       adiabatic representation (used by ANT program)
!   xx(1,i)             X coordinate of atom i
!   xx(2,i)             Y coordinate of atom i
!   xx(3,i)             Z coordinate of atom i
!
!   Output: 
!   uu(j,j)             diabatic potential of diabatic state j
!   uu(j,k)             diabatic coupling between state j and k 
!   guu(1,i,j,j)        X component of the gradients of diabatic 
!                       potential of state j at atom i
!   guu(2,i,j,j)        Y component of the gradients of diabatic 
!                       potential of state j at atom i
!   guu(3,i,j,j)        Z component of the gradients of diabatic
!                       potential of state j at atom i
!   guu(1,i,j,k)        X component of the gradients of diabatic         
!                       coupling between state j and k at atom i
!   guu(2,i,j,k)        Y component of the gradients of diabatic         
!                       coupling between state j and k at atom i
!   guu(3,i,j,k)        Z component of the gradients of diabatic
!                       coupling between state j and k at atom i
!
!   vv(j)               adiabatic potential of adiabatic state j
!   gvv(1,i,j)          X component of the gradients of adiabatic 
!                       potential of state j at atom i
!   gvv(2,i,j)          Y component of the gradients of adiabatic
!                       potential of state j at atom i
!   gvv(3,i,j)          Z component of the gradients of adiabatic
!                       potential of state j at atom i
!   dvec(1,i,j,k)       X component of nonadiabatic coupling between
!                       state j and k at atom i
!   dvec(2,i,j,k)       Y component of nonadiabatic coupling between
!                       state j and k at atom i
!   dvec(3,i,j,k)       Z component of nonadiabatic coupling between
!                       state j and k at atom i
!   eigvec(2,2)         2*2 orthogonal matrix diagonalizing UU matrix
!                       (CC^T*UU*CC yield diagonal matrix with adiabatic
!                       energies as diagonal elements)              
!
!   Note:       LAPACK library is needed to diagonize diabatic potential
!               matrix (UU) to yield adiabatic potenitials VV(2)
!
!***********************************************************************

      subroutine pes(x,igrad,p,g,d)

      implicit none
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=7
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: xx(3,natoms), uu(nstates,nstates)
      double precision :: guu(3,natoms,nstates,nstates)
      double precision :: vv(nstates), gvv(3,natoms,nstates)
      double precision :: dvec(3,natoms,nstates,nstates)
      double precision :: cc(nstates,nstates)
      integer :: iatom, idir, j, istate, jstate
      integer :: jgrad
      !initialize 
      p=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
        do idir=1,3
          xx(idir,iatom)=x(iatom,idir)  !/0.529177211
        enddo
      enddo

      if (igrad==0) then
        jgrad=0
      else
        jgrad=1
      endif

      call pot(jgrad,xx,uu,guu,vv,gvv,dvec,cc,0)

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=vv(istate)*27.211386
        enddo
      elseif (igrad==1) then
        do istate=1,nstates
          p(istate)=vv(istate)*27.211386
        enddo
        do istate=1,nstates
        do iatom=1,natoms
        do idir=1,3
          g(istate,iatom,idir)=gvv(idir,iatom,istate)*51.422067
        enddo
        enddo
        enddo
      elseif (igrad==2) then
        do istate=1,nstates
          p(istate)=vv(istate)*27.211386
        enddo
        do istate=1,nstates
        do iatom=1,natoms
        do idir=1,3
          g(istate,iatom,idir)=gvv(idir,iatom,istate)*51.422067
        enddo
        enddo
        enddo
        do istate=1,nstates
        do jstate=1,nstates
        do iatom=1,natoms
        do idir=1,3
          d(istate,jstate,iatom,idir)=dvec(idir,iatom,istate,jstate)/0.529177211
        enddo
        enddo
        enddo
        enddo
      endif

      endsubroutine


      subroutine pot(igrad,xx,uu,guu,vv,gvv,dvec,eigvec,repflag)
      implicit none
      double precision xx(3,7)
      double precision X(7),Y(7),Z(7)
      double precision scp(3)
      integer i,j
      double precision e1
      integer nrap,ntap
      double precision dx1(7),dy1(7),dz1(7)
      double precision dx(7),dy(7),dz(7)
      integer igrad
      double precision dUx(7,3),dUy(7,3),dUz(7,3)
      double precision uu(2,2),vv(2)
      double precision tt(1)
      double precision guu(3,7,2,2),gtt(3,7,1,1)
      double precision gvv(3,7,2)
      double precision dvec(3,7,2,2)
      double precision eigvec(2,2)
      integer lwork,info
      integer repflag
      double precision, allocatable :: work(:)
      double precision tmpmat22(2,2)
      integer ist,jst
      double precision, parameter :: evtohartree=27.21139d0
      double precision, parameter :: bohrtoang=0.52917721067d0
      double precision v2,v4

      scp=0d0
      dUx=0d0
      dUy=0d0
      dUz=0d0

      do i=1,7
      X(i)=xx(1,i)
      Y(i)=xx(2,i)
      Z(i)=xx(3,i)
      enddo

! get tertiary then primary contributions
      do j=1,3
       nrap=4
       ntap=3
       call tert(igrad,j,X,Y,Z,nrap,ntap,e1,dx1,dy1,dz1)

        scp(j)=e1
        do i=1,7
         dUx(i,j)=dx1(i)
         dUy(i,j)=dy1(i)
         dUz(i,j)=dz1(i)
        enddo

       call KPev2gm2(igrad,X,Y,Z,v2,j,dx,dy,dz)
       scp(j)=scp(j)+v2
        do i=1,7
         dUx(i,j)=dUx(i,j)+dx(i)
         dUy(i,j)=dUy(i,j)+dy(i)
         dUz(i,j)=dUz(i,j)+dz(i)
        enddo


      enddo

! change units from eV to Hartree and Angstrom to Bohr

       uu(1,1)=scp(1)/evtohartree
       uu(2,2)=scp(2)/evtohartree
       uu(1,2)=scp(3)/evtohartree
       uu(2,1)=scp(3)/evtohartree

       eigvec=uu

      lwork = -1
      allocate (work (1) )
      call dsyev('V','L',2,eigvec,2,vv,work,lwork,info)
      lwork=int(work(1))
      deallocate (work)
      allocate( work(lwork) )
      call dsyev('V','L',2,eigvec,2,vv,work,lwork,info)
      deallocate (work)
      if(info.ne.0) then
       write(6,*) 'Error in diagolization of U matrix'
       stop
      endif

       do i=1,7
       guu(1,i,1,1)=(dUx(i,1)/evtohartree)*bohrtoang
       guu(2,i,1,1)=(dUy(i,1)/evtohartree)*bohrtoang
       guu(3,i,1,1)=(dUz(i,1)/evtohartree)*bohrtoang
!
       guu(1,i,2,2)=(dUx(i,2)/evtohartree)*bohrtoang
       guu(2,i,2,2)=(dUy(i,2)/evtohartree)*bohrtoang
       guu(3,i,2,2)=(dUz(i,2)/evtohartree)*bohrtoang
!
       guu(1,i,1,2)=(dUx(i,3)/evtohartree)*bohrtoang
       guu(2,i,1,2)=(dUy(i,3)/evtohartree)*bohrtoang
       guu(3,i,1,2)=(dUz(i,3)/evtohartree)*bohrtoang
!
       guu(1,i,2,1)=(dUx(i,3)/evtohartree)*bohrtoang
       guu(2,i,2,1)=(dUy(i,3)/evtohartree)*bohrtoang
       guu(3,i,2,1)=(dUz(i,3)/evtohartree)*bohrtoang

       enddo

! convert diabatic potential energy matrix and its gradients to
! adiabatic potentials and gradients

       do i=1,3
        do j=1,7
         tmpmat22(:,:) = guu(i,j,:,:)
         tmpmat22 = matmul(transpose(eigvec), 
     &                     matmul(tmpmat22,eigvec))
         do ist=1,2
          gvv(i,j,ist)=tmpmat22(ist,ist)
         enddo
!
         ist=1
         jst=2
         dvec(i,j,ist,jst)=tmpmat22(ist,jst)/
     &                    (vv(jst)-vv(ist))
         dvec(i,j,jst,ist)= -dvec(i,j,ist,jst)
        enddo
       enddo

      end


      double precision function dihedr1(i,j,k,l,X,Y,Z)
!************************************************************************
! Evaluate the dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
! I created a modified version of this code to avoid issues when the
! amine is flat in a plane
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   
      rij=bndlen(i,j,X,Y,Z)
      rjk=bndlen(j,k,X,Y,Z)
      rkl=bndlen(k,l,X,Y,Z)
      aijk=bndang(i,j,k,X,Y,Z)
      ajkl=bndang(j,k,l,X,Y,Z)

! Check whether i-j-k or j-k-l are collinear
      if ((abs(aijk) .le. tiny) .or. (abs(aijk-pi) .lt. tiny)) then 
        write(*,9999) i,j,k
      endif
      if ((abs(ajkl) .le. tiny) .or. (abs(ajkl-pi) .lt. tiny)) then
        write(*,9999) j,k,l
      endif
 9999 Format('Torsion angle not defined since atoms',3I3,' collinear')

! Calculate eij
      eij(1)=(X(j)-X(i))/rij
      eij(2)=(Y(j)-Y(i))/rij
      eij(3)=(Z(j)-Z(i))/rij
! Calculate ejk
      ejk(1)=(X(k)-X(j))/rjk
      ejk(2)=(Y(k)-Y(j))/rjk
      ejk(3)=(Z(k)-Z(j))/rjk
! Calculate ekl
      ekl(1)=(X(l)-X(k))/rkl
      ekl(2)=(Y(l)-Y(k))/rkl
      ekl(3)=(Z(l)-Z(k))/rkl

!
! u = eij*ejk normal vector in plane i-j-k
! v = ejk*ekl Normal vector in plane j-k-l
! w = u*v used to determine the sign of dihedral angle
!    
      call xprod1(eij,ejk,u)
      call xprod1(ejk,ekl,v)
      call xprod1(u,v,w)

      call sprod(3,u,v,a)

      a=a/(dsin(aijk)*dsin(ajkl))
      a=max(a,-1.0d0)
      a=min(a,1.0d0)
      dihedr1=dacos(a)

! Determine the sign of dihedr based on u,v,ejk
! dihedr >= 0, u,v,ejk right-hand (w*ejk >=0)
! dihedr < 0,  u,v,ejk left-hand (w*ejk < 0)

      wejk=w(1)*ejk(1)+w(2)*ejk(2)+w(3)*ejk(3)
      if (wejk .lt. 0.0d0) then
        dihedr1=-dihedr1
      endif

! Setting the range of dihedral angle between -90.0 degree to 270 degree
! to avoid discontinuity around 180 degrees
!      if (dihedr .lt. -0.5d0*pi) then
!        dihedr=2.0d0*pi+dihedr
!      endif

      end


      double precision function bndlen(i,j,X,Y,Z)
!************************************************************************
!* Function to evaluate the bond length between atom i and j
!* Input: i,j (indices of 2 atoms)
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny
      parameter(tiny=1.0d-13)

      double precision X(*),Y(*),Z(*)

      bndlen=sqrt( (x(i)-x(j))*(x(i)-x(j)) 
     &    +      (y(i)-y(j))*(y(i)-y(j)) 
     &    +      (z(i)-z(j))*(z(i)-z(j)))

      if (bndlen .le. tiny) then
        write(*,9999) i,j
!CDB  KRY 20140718
        write(*,'(I3,3F18.12)') i,x(i),y(i),z(i)
        write(*,'(I3,3F18.12)') j,x(j),y(j),z(j)
        stop
      endif

 9999 Format(' Atom',I3,' and',I3,' are too close to each other!')

      end

      double precision function bndang(i,j,k,X,Y,Z)
!************************************************************************
!* Evaluate the bond angle between atoms i, j, and k
!* Input: i,j,k (indices of 3 atoms, i-j and j-k bonds)
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision X(*),Y(*),Z(*)
      double precision eji(3),ejk(3),rji,rjk

      rji=bndlen(i,j,X,Y,Z)
      rjk=bndlen(j,k,X,Y,Z)

      eji(1)=(X(i)-X(j))/rji
      eji(2)=(Y(i)-Y(j))/rji
      eji(3)=(Z(i)-Z(j))/rji

      ejk(1)=(X(k)-X(j))/rjk
      ejk(2)=(Y(k)-Y(j))/rjk
      ejk(3)=(Z(k)-Z(j))/rjk

      call sprod(3,eji,ejk,a)
      bndang=acos(a)

      end

      subroutine xprod1(x1,x2,x3)
!************************************************************************
!* Calculate cross product of x1 and x2 and return x3
!* Input: X1,X2
!* Output: X3
!************************************************************************

       implicit double precision (a-h,o-z)

       double precision x1(3),x2(3),x3(3)

       x3(1)=x1(2)*x2(3)-x1(3)*x2(2)
       x3(2)=x1(3)*x2(1)-x1(1)*x2(3)
       x3(3)=x1(1)*x2(2)-x1(2)*x2(1)

       end
      subroutine sprod(n,a,b,c)
!************************************************************************
!* Calculate scalar product of a and b and return c
!* Input: a(n),b(n)
!* Output: c
!************************************************************************

       implicit double precision (a-h,o-z)

       integer n
       double precision a(*),b(*),c

       c=0.0d0
       do i=1,n
         c=c+a(i)*b(i)
       enddo

       end

      subroutine convertit(X,Y,Z,dk,dx1,dy1,dz1)
! this subroutine deals with the tertiary derivatives and converts these
! from the various internal coords to XYZ
      implicit none
      integer kk,nat
      integer i,j,k
      double precision :: xyz(4,3)
      double precision beta
      double precision X(7),Y(7),Z(7)
      double precision dk(22)
      double precision r1,t2
      double precision qc(20)
      double precision dx1(7),dy1(7),dz1(7)
      double precision bndlen
      double precision dngX(7),dngY(7),dngZ(7)
      double precision dihedrIX,dihedrIY,dihedrIZ
      double precision dihedrJX,dihedrJY,dihedrJZ
      double precision dihedrKX,dihedrKY,dihedrKZ
      double precision dihedrLX,dihedrLY,dihedrLZ
      double precision dihedr2
      double precision torH1,torH2,rH1HX,rH2HX
      double precision ang1,ang2,ang3
      double precision len1,len2,len3
      integer bonds(2,8),bends(3,8)
      data bonds   /2,3,
     &              3,4,
     &              4,5,
     &              4,6,
     &              4,7,
!
     &              5,6,
     &              6,7,
     &              7,5/
      data bends   /1,3,4,
     &              2,3,4,
     &              3,4,5,
     &              3,4,6,
     &              3,4,7,
!
     &              5,4,6,
     &              6,4,7,
     &              7,4,5/


      r1=bndlen(1,3,X,Y,Z)
      t2=dihedr2(1,3,4,2,X,Y,Z)

       call calcgeom(X,Y,Z,qc)

        torH1=qc(11)
        torH2=qc(12)

        rH1HX=qc(13)
        rH2HX=qc(14)

        ang1=qc(15)
        ang2=qc(16)
        ang3=qc(17)

        len1=qc(18)
        len2=qc(19)
        len3=qc(20)



      dx1=0d0
      dy1=0d0
      dz1=0d0

      do kk=1,5
      i=bonds(1,kk)
      j=bonds(2,kk)
       call disder(X,Y,Z,i,j,dngX,dngY,dngZ)
       do nat=1,7
        dx1(nat)=dx1(nat)+dk(kk)*dngX(nat)
        dy1(nat)=dy1(nat)+dk(kk)*dngY(nat)
        dz1(nat)=dz1(nat)+dk(kk)*dngZ(nat)
       enddo
      enddo

      do kk=1,5
      i=bends(1,kk)
      j=bends(2,kk)
      k=bends(3,kk)
       call angder(X,Y,Z,i,j,k,dngX,dngY,dngZ) 
       do nat=1,7
        dx1(nat)=dx1(nat)+dk(kk+5)*dngX(nat)
        dy1(nat)=dy1(nat)+dk(kk+5)*dngY(nat)
        dz1(nat)=dz1(nat)+dk(kk+5)*dngZ(nat)
       enddo
      enddo



! 11th/12th dk COORD

      dx1(1)=dx1(1)+dk(11)*((X(1)-X(5))/rH1HX)
      dy1(1)=dy1(1)+dk(11)*((Y(1)-Y(5))/rH1HX)
      dz1(1)=dz1(1)+dk(11)*((Z(1)-Z(5))/rH1HX)

      dx1(5)=dx1(5)+dk(11)*((X(5)-X(1))/rH1HX)
      dy1(5)=dy1(5)+dk(11)*((Y(5)-Y(1))/rH1HX)
      dz1(5)=dz1(5)+dk(11)*((Z(5)-Z(1))/rH1HX)
!!!!!!!!

      dx1(2)=dx1(2)+dk(12)*((X(2)-X(5))/rH2HX)
      dy1(2)=dy1(2)+dk(12)*((Y(2)-Y(5))/rH2HX)
      dz1(2)=dz1(2)+dk(12)*((Z(2)-Z(5))/rH2HX)

      dx1(5)=dx1(5)+dk(12)*((X(5)-X(2))/rH2HX)
      dy1(5)=dy1(5)+dk(12)*((Y(5)-Y(2))/rH2HX)
      dz1(5)=dz1(5)+dk(12)*((Z(5)-Z(2))/rH2HX)
!!!!

! 13th/14th dk COORD

!!!!
! H1
      dx1(1)=dx1(1)+dk(13)*dihedrIX(1,3,4,5,X,Y,Z)
      dy1(1)=dy1(1)+dk(13)*dihedrIY(1,3,4,5,X,Y,Z)
      dz1(1)=dz1(1)+dk(13)*dihedrIZ(1,3,4,5,X,Y,Z)
! N
      dx1(3)=dx1(3)+dk(13)*dihedrJX(1,3,4,5,X,Y,Z)
      dy1(3)=dy1(3)+dk(13)*dihedrJY(1,3,4,5,X,Y,Z)
      dz1(3)=dz1(3)+dk(13)*dihedrJZ(1,3,4,5,X,Y,Z)
! C
      dx1(4)=dx1(4)+dk(13)*dihedrKX(1,3,4,5,X,Y,Z)
      dy1(4)=dy1(4)+dk(13)*dihedrKY(1,3,4,5,X,Y,Z)
      dz1(4)=dz1(4)+dk(13)*dihedrKZ(1,3,4,5,X,Y,Z)
! HX
      dx1(5)=dx1(5)+dk(13)*dihedrLX(1,3,4,5,X,Y,Z)
      dy1(5)=dy1(5)+dk(13)*dihedrLY(1,3,4,5,X,Y,Z)
      dz1(5)=dz1(5)+dk(13)*dihedrLZ(1,3,4,5,X,Y,Z)

!!!!
! H2
      dx1(2)=dx1(2)+dk(14)*dihedrIX(2,3,4,5,X,Y,Z)
      dy1(2)=dy1(2)+dk(14)*dihedrIY(2,3,4,5,X,Y,Z)
      dz1(2)=dz1(2)+dk(14)*dihedrIZ(2,3,4,5,X,Y,Z)
! N
      dx1(3)=dx1(3)+dk(14)*dihedrJX(2,3,4,5,X,Y,Z)
      dy1(3)=dy1(3)+dk(14)*dihedrJY(2,3,4,5,X,Y,Z)
      dz1(3)=dz1(3)+dk(14)*dihedrJZ(2,3,4,5,X,Y,Z)
! C
      dx1(4)=dx1(4)+dk(14)*dihedrKX(2,3,4,5,X,Y,Z)
      dy1(4)=dy1(4)+dk(14)*dihedrKY(2,3,4,5,X,Y,Z)
      dz1(4)=dz1(4)+dk(14)*dihedrKZ(2,3,4,5,X,Y,Z)
! HX
      dx1(5)=dx1(5)+dk(14)*dihedrLX(2,3,4,5,X,Y,Z)
      dy1(5)=dy1(5)+dk(14)*dihedrLY(2,3,4,5,X,Y,Z)
      dz1(5)=dz1(5)+dk(14)*dihedrLZ(2,3,4,5,X,Y,Z)


      do kk=6,8
      i=bends(1,kk)
      j=bends(2,kk)
      k=bends(3,kk)
       call angder(X,Y,Z,i,j,k,dngX,dngY,dngZ)
       do nat=1,7
        dx1(nat)=dx1(nat)+dk(9+kk)*dngX(nat)
        dy1(nat)=dy1(nat)+dk(9+kk)*dngY(nat)
        dz1(nat)=dz1(nat)+dk(9+kk)*dngZ(nat)
       enddo
      enddo

      do kk=6,8
      i=bonds(1,kk)
      j=bonds(2,kk)
       call disder(X,Y,Z,i,j,dngX,dngY,dngZ)
       do nat=1,7
        dx1(nat)=dx1(nat)+dk(12+kk)*dngX(nat)
        dy1(nat)=dy1(nat)+dk(12+kk)*dngY(nat)
        dz1(nat)=dz1(nat)+dk(12+kk)*dngZ(nat)
       enddo
      enddo


! 21th dk COORD

      dx1(1)=dx1(1)+dk(21)*(X(1)-X(3))/r1
      dy1(1)=dy1(1)+dk(21)*(Y(1)-Y(3))/r1
      dz1(1)=dz1(1)+dk(21)*(Z(1)-Z(3))/r1

      dx1(3)=dx1(3)+dk(21)*(X(3)-X(1))/r1
      dy1(3)=dy1(3)+dk(21)*(Y(3)-Y(1))/r1
      dz1(3)=dz1(3)+dk(21)*(Z(3)-Z(1))/r1

! 22th dk COORD

      dx1(1)=dx1(1)+dk(22)*dihedrIX(1,3,4,2,X,Y,Z)
      dy1(1)=dy1(1)+dk(22)*dihedrIY(1,3,4,2,X,Y,Z)
      dz1(1)=dz1(1)+dk(22)*dihedrIZ(1,3,4,2,X,Y,Z)

      dx1(2)=dx1(2)+dk(22)*dihedrLX(1,3,4,2,X,Y,Z) 
      dy1(2)=dy1(2)+dk(22)*dihedrLY(1,3,4,2,X,Y,Z) 
      dz1(2)=dz1(2)+dk(22)*dihedrLZ(1,3,4,2,X,Y,Z) 
! N
      dx1(3)=dx1(3)+dk(22)*dihedrJX(1,3,4,2,X,Y,Z)
      dy1(3)=dy1(3)+dk(22)*dihedrJY(1,3,4,2,X,Y,Z) 
      dz1(3)=dz1(3)+dk(22)*dihedrJZ(1,3,4,2,X,Y,Z) 
! C
      dx1(4)=dx1(4)+dk(22)*dihedrKX(1,3,4,2,X,Y,Z) 
      dy1(4)=dy1(4)+dk(22)*dihedrKY(1,3,4,2,X,Y,Z) 
      dz1(4)=dz1(4)+dk(22)*dihedrKZ(1,3,4,2,X,Y,Z) 


      end

      subroutine angder(X,Y,Z,i,j,k,dngX,dngY,dngZ)
! This subroutine converts a bond angle derivative to XYZ
      implicit none
      integer i,j,k
      double precision X(7),Y(7),Z(7)
      double precision ra,rb
      double precision to,be,go,cgo
      double precision dngX(7),dngY(7),dngZ(7)
      double precision bndlen

      dngX=0d0
      dngY=0d0
      dngZ=0d0

      ra=bndlen(j,i,X,Y,Z)
      rb=bndlen(j,k,X,Y,Z)

      to=(X(j)-X(i))*(X(j)-X(k))
     &  +(Y(j)-Y(i))*(Y(j)-Y(k))
     &  +(Z(j)-Z(i))*(Z(j)-Z(k))


      be=((((X(j)-X(i))**2)+((Y(j)-Y(i))**2)+((Z(j)-Z(i))**2))**0.5d0)
     &  *((((X(j)-X(k))**2)+((Y(j)-Y(k))**2)+((Z(j)-Z(k))**2))**0.5d0)

      go=to/be
      cgo=(-1d0/((1d0-go**2d0)**(1d0/2d0)))

!--------- j -----------!

      dngX(j)=cgo*((2d0*X(j)-X(i)-X(k))*be
     &  -( (((X(j)-X(i))/ra)*rb) + (((X(j)-X(k))/rb)*ra) )*to)/(be**2d0)

      dngY(j)=cgo*((2d0*Y(j)-Y(i)-Y(k))*be
     &  -( (((Y(j)-Y(i))/ra)*rb) + (((Y(j)-Y(k))/rb)*ra) )*to)/(be**2d0)

      dngZ(j)=cgo*((2d0*Z(j)-Z(i)-Z(k))*be
     &  -( (((Z(j)-Z(i))/ra)*rb) + (((Z(j)-Z(k))/rb)*ra) )*to)/(be**2d0)


!--------- k -----------!
      dngX(k)=cgo*( -1d0*(X(j)-X(i))*be
     &       -(((X(k)-X(j))/rb)*ra*to))/(be**2d0)

      dngY(k)=cgo*( -1d0*(Y(j)-Y(i))*be
     &       -(((Y(k)-Y(j))/rb)*ra*to))/(be**2d0)

      dngZ(k)=cgo*( -1d0*(Z(j)-Z(i))*be
     &       -(((Z(k)-Z(j))/rb)*ra*to))/(be**2d0)



!--------- i -----------!
      dngX(i)=cgo*( -(X(j)-X(k))*be
     &       -(((X(i)-X(j))/ra)*rb*to))/(be**2d0)

      dngY(i)=cgo*( -(Y(j)-Y(k))*be
     &       -(((Y(i)-Y(j))/ra)*rb*to))/(be**2d0)

      dngZ(i)=cgo*( -(Z(j)-Z(k))*be
     &       -(((Z(i)-Z(j))/ra)*rb*to))/(be**2d0)


      end


      subroutine disder(X,Y,Z,i,j,dngX,dngY,dngZ)
! This subroutine converts a distance derivative to XYZ
      implicit none
      integer i,j,k
      double precision X(7),Y(7),Z(7)
      double precision ra,rb
      double precision to,be,go,cgo
      double precision dngX(7),dngY(7),dngZ(7)
      double precision bndlen

      dngX=0d0
      dngY=0d0
      dngZ=0d0

      dngX(i)=(X(i)-X(j))/bndlen(i,j,X,Y,Z)
      dngY(i)=(Y(i)-Y(j))/bndlen(i,j,X,Y,Z)
      dngZ(i)=(Z(i)-Z(j))/bndlen(i,j,X,Y,Z)

      dngX(j)=(X(j)-X(i))/bndlen(i,j,X,Y,Z)
      dngY(j)=(Y(j)-Y(i))/bndlen(i,j,X,Y,Z)
      dngZ(j)=(Z(j)-Z(i))/bndlen(i,j,X,Y,Z)

      end


      subroutine solver(igrad,ra1,ra2,qc,ii,nrap,ntap,ut,dut)
!***********************************************************************
! Subroutine to evaluate the tertiary contributions to the 
! diabatic potentials and gradients at different anchor points
! ii:   flag to indicate the diabatic states
!       = 1, U11
!       = 2, U22
!       = 3, U12
!***********************************************************************

      implicit none 

      integer nrap,ntap,ra1,ra2
      integer nqs,nqb
      parameter (nqs=5,nqb=5)
      integer igrad,ii
      integer ntot,nqtc
      double precision r,ut,dutdr,ra

C Local variables
      integer i,j,k
      double precision u0
      double precision vs,vb,vbs
      double precision qs(nqs)
      double precision qb(nqb)
      double precision sqb(nqb),cqb(nqb)
      double precision qc(20),dut(20)
      double precision v0(nrap,ntap)
      double precision PI
      double precision dnc1,dnc2
      double precision dch1,dch2
      double precision KNHB1(nrap,ntap),KNHB2(nrap,ntap)
      double precision KNC1(nrap,ntap),KNC2(nrap,ntap)
      double precision KCH1(nrap,ntap),KCH2(nrap,ntap)
      double precision KCH3(nrap,ntap)
      double precision KCNHA1(nrap,ntap),KCNHA2(nrap,ntap)
      double precision KCNHB1(nrap,ntap),KCNHB2(nrap,ntap)
      double precision KNCH1(nrap,ntap),KNCH2(nrap,ntap)
!
      double precision JNHB1(nrap,ntap),JNHB2(nrap,ntap)
      double precision JNHB3(nrap,ntap)
      double precision JNC1(nrap,ntap),JNC2(nrap,ntap)
      double precision JNC3(nrap,ntap)
      double precision JCNHA1(nrap,ntap),JCNHA2(nrap,ntap)
      double precision JCNHA3(nrap,ntap)
      double precision JCNHB1(nrap,ntap),JCNHB2(nrap,ntap)
      double precision JCNHB3(nrap,ntap)
!
      double precision GNHB1(nrap,ntap),GNHB2(nrap,ntap)
      double precision GNHB3(nrap,ntap)
!
      double precision GNC1(nrap,ntap),GNC2(nrap,ntap)
      double precision GNC3(nrap,ntap)
!
      double precision v3,v4
      double precision torH1,torH2,rH1HX,rH2HX
      double precision ang1,ang2,ang3
      double precision len1,len2,len3
      double precision DD1,DD2,RR,ft2,fa2
      double precision ft1(nrap,ntap)
      double precision fa1(nrap,ntap)

      PI=4.D0*DATAN(1.D0)

          ut=0d0
          vs=0d0
          vb=0d0
          dut=0d0


          do i=1,nqs
          qs(i)= qc(i)
          enddo

          do i=1,nqb
          qb(i)= qc(i+nqs) 
          sqb(i)=sin(qb(i))
          cqb(i)=cos(qb(i))
          enddo


        torH1=qc(11)
        torH2=qc(12)

        rH1HX=qc(13)
        rH2HX=qc(14)

        ang1=qc(15)
        ang2=qc(16)
        ang3=qc(17)

        len1=qc(18)
        len2=qc(19)
        len3=qc(20)


       call getpara(ii,nrap,ntap,
     & KNHB1,KNHB2,
     & KNC1,KNC2,
     & KCH1,KCH2,
     & KCH3,
     & KCNHA1,KCNHA2,
     & KCNHB1,KCNHB2,
     & KNCH1,KNCH2,
     & JNHB1,JNHB2,
     & JNHB3,
     & JNC1,JNC2,JNC3,JCNHA1,JCNHA2,
     & JCNHA3,JCNHB1,JCNHB2,JCNHB3,
     & GNHB1,GNHB2,GNHB3,
     & GNC1,GNC2,GNC3,
     & ft1,fa1)


        call getV0(ii,nrap,ntap,v0)

C Calculate potential energy
! CONSTANT
            u0=v0(ra1,ra2)
       if ((ii .eq. 1) .or. (ii .eq. 2)) then
! STRETCHES
! STRETCH 1 (NHB)
!       if (ii .eq. 1) then
        vs=GNHB1(ra1,ra2)
     &       *(1d0-exp(GNHB2(ra1,ra2)*(qc(1)-GNHB3(ra1,ra2))))**2d0
        vs=vs+KNHB1(ra1,ra2)*(qc(1)-KNHB2(ra1,ra2))**4d0
! STRETCH 2 (NC)
        vs=vs+GNC1(ra1,ra2)
     &       *(1d0-exp(GNC2(ra1,ra2)*(qc(2)-GNC3(ra1,ra2))))**2d0
        vs=vs+KNC1(ra1,ra2)*(qc(2)-KNC2(ra1,ra2))**4d0



! STRETCHES 3-5 (CH)
        do i=3,nqs

        vs=vs+KCH1(ra1,ra2)
     &       *(1d0-exp(KCH2(ra1,ra2)*(qc(i)-KCH3(ra1,ra2))))**2d0
        enddo

! BENDS
! BEND 1 (CNHA)
        vb=KCNHA1(ra1,ra2)*(cqb(1)-KCNHA2(ra1,ra2))**2d0
! BEND 2 (CNHB)
        vb=vb+KCNHB1(ra1,ra2)*(cqb(2)-KCNHB2(ra1,ra2))**2d0
! BEND 3-5 (NCH)
        do i=3,nqb
        vb=vb+KNCH1(ra1,ra2)*(qb(i)-KNCH2(ra1,ra2))**2d0
        enddo

      ft2=-3.0d0
      DD1=2.2d0

      v3=
     &   ft1(ra1,ra2)*(1d0-(4d0*(torH2**3d0)-3d0*torH2))
     & * exp(((rH2HX-DD1)**2d0)/ft2)


      DD2=1.881d0
      RR=1.761d0
      fa2=-2.0d0

      v4=(fa1(ra1,ra2)*(ang1-DD2)**2d0)
     & *exp(((len1-RR)**2d0)/fa2)
     &  +(fa1(ra1,ra2)*(ang2-DD2)**2d0)
     & *exp(((len2-RR)**2d0)/fa2)
     &  +(fa1(ra1,ra2)*(ang3-DD2)**2d0)
     & *exp(((len3-RR)**2d0)/fa2)


       elseif (ii .eq. 3) then
! STRETCHES
! STRETCH 1 (NHB)
        vs=JNHB1(ra1,ra2)*exp(JNHB2(ra1,ra2)
     &    *(qc(1)-JNHB3(ra1,ra2))**2d0)
! STRETCH 2 (NC)
        vs=vs+JNC1(ra1,ra2)*exp(JNC2(ra1,ra2)
     &    *(qc(2)-JNC3(ra1,ra2))**2d0)
! STRETCHES 3-5 (CH)





! BENDS
! BEND 1 (CNHA)
        vb=JCNHA1(ra1,ra2)*exp(JCNHA2(ra1,ra2)
     &    *(cqb(1)-JCNHA3(ra1,ra2))**2d0)
! BEND 2 (CNHB)
        vb=vb+JCNHB1(ra1,ra2)*exp(JCNHB2(ra1,ra2)
     &    *(cqb(2)-JCNHB3(ra1,ra2))**2d0)
! BEND 3-5 (NCH)




       endif

! COMBINE
        ut=vs+vb+u0+v3+v4



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! GRADIENT
      if (igrad .eq. 1) then
       if ((ii .eq. 1) .or. (ii .eq. 2)) then
! STRETCHES
! STRETCH 1 (NHB)
!       if (ii .eq. 1) then
        dut(1)=-2d0*GNHB2(ra1,ra2)*GNHB1(ra1,ra2)*
     &    exp(GNHB2(ra1,ra2)*(qc(1)-GNHB3(ra1,ra2)))*
     &    (1d0-exp(GNHB2(ra1,ra2)*(qc(1)-GNHB3(ra1,ra2))))
        dut(1)=dut(1)+4d0*KNHB1(ra1,ra2)*(qc(1)-KNHB2(ra1,ra2))**3d0





! STRETCH 2 (NC)

        dut(2)=-2d0*GNC2(ra1,ra2)*GNC1(ra1,ra2)*
     &    exp(GNC2(ra1,ra2)*(qc(2)-GNC3(ra1,ra2)))*
     &    (1d0-exp(GNC2(ra1,ra2)*(qc(2)-GNC3(ra1,ra2))))
        dut(2)=dut(2)+4d0*KNC1(ra1,ra2)*(qc(2)-KNC2(ra1,ra2))**3d0




! STRETCHES 3-5 (CH)
        do i=3,nqs



       dut(i)=-2d0*KCH2(ra1,ra2)*KCH1(ra1,ra2)*
     &    exp(KCH2(ra1,ra2)*(qs(i)-KCH3(ra1,ra2)))*
     &    (1d0-exp(KCH2(ra1,ra2)*(qs(i)-KCH3(ra1,ra2))))
        enddo

! BENDS
! BEND 1 (CNHA)
        dut(6)=KCNHA1(ra1,ra2)*(cqb(1)-KCNHA2(ra1,ra2))
     &        *2d0*(-sqb(1))
! BEND 2 (CNHB)
        dut(7)=KCNHB1(ra1,ra2)*(cqb(2)-KCNHB2(ra1,ra2))
     &        *2d0*(-sqb(2))
! BEND 3-5 (NCH)
        do i=3,nqb
        dut(nqs+i)=KNCH1(ra1,ra2)*(qb(i)-KNCH2(ra1,ra2))*2d0
        enddo


      dut(11)=2d0*((rH1HX-DD1)/ft2)*exp(((rH1HX-DD1)**2d0)/ft2)
     &       *ft1(ra1,ra2)*(1d0-(4d0*(torH1**3d0)-3d0*torH1))
      dut(12)=2d0*((rH2HX-DD1)/ft2)*exp(((rH2HX-DD1)**2d0)/ft2)
     &       *ft1(ra1,ra2)*(1d0-(4d0*(torH2**3d0)-3d0*torH2))

      dut(13)=ft1(ra1,ra2)*(-12d0*(torH1**2d0)+3d0)
     & * exp(((rH1HX-DD1)**2d0)/ft2)
      dut(14)=ft1(ra1,ra2)*(-12d0*(torH2**2d0)+3d0)
     & * exp(((rH2HX-DD1)**2d0)/ft2)


      dut(15)=2d0*fa1(ra1,ra2)*(ang1-DD2)
     & *exp(((len1-RR)**2d0)/fa2)
      dut(16)=2d0*fa1(ra1,ra2)*(ang2-DD2)
     & *exp(((len2-RR)**2d0)/fa2)
      dut(17)=2d0*fa1(ra1,ra2)*(ang3-DD2)
     & *exp(((len3-RR)**2d0)/fa2)

      dut(18)=(fa1(ra1,ra2)*(ang1-DD2)**2d0)
     &  *2d0*((len1-RR)/fa2)*exp(((len1-RR)**2d0)/fa2)
!
      dut(19)=(fa1(ra1,ra2)*(ang2-DD2)**2d0)
     &  *2d0*((len2-RR)/fa2)*exp(((len2-RR)**2d0)/fa2)
!
      dut(20)=(fa1(ra1,ra2)*(ang3-DD2)**2d0)
     &  *2d0*((len3-RR)/fa2)*exp(((len3-RR)**2d0)/fa2)

      dut(11)=0d0
      dut(13)=0d0

       elseif (ii .eq. 3) then
! STRETCH 1 (NHB)
        dut(1)=JNHB2(ra1,ra2)*(qc(1)-JNHB3(ra1,ra2))*2d0*
     &    JNHB1(ra1,ra2)*exp(JNHB2(ra1,ra2)
     &    *(qc(1)-JNHB3(ra1,ra2))**2d0)
! STRETCH 2 (NC)
        dut(2)=JNC2(ra1,ra2)*(qc(2)-JNC3(ra1,ra2))*2d0*
     &    JNC1(ra1,ra2)*exp(JNC2(ra1,ra2)
     &    *(qc(2)-JNC3(ra1,ra2))**2d0)

! BENDS
! BEND 1 (CNHA)
        dut(6)=JCNHA2(ra1,ra2)*(cqb(1)-JCNHA3(ra1,ra2))*2d0*
     &    JCNHA1(ra1,ra2)*exp(JCNHA2(ra1,ra2)
     &    *(cqb(1)-JCNHA3(ra1,ra2))**2d0)
     &    *(-sqb(1))
! BEND 2 (CNHB)
        dut(7)=JCNHB2(ra1,ra2)*(cqb(2)-JCNHB3(ra1,ra2))*2d0
     &    *JCNHB1(ra1,ra2)*exp(JCNHB2(ra1,ra2)
     &    *(cqb(2)-JCNHB3(ra1,ra2))**2d0)
     &    *(-sqb(2))
       endif
      endif

      end


       subroutine getV0(ii,nrap,ntap,v0)
! subrountine gives constants in tertiary eqns
       implicit none
       integer ii
       integer nrap,ntap
       double precision v0(nrap,ntap)

        v0=0d0

        if (ii .eq. 1) then
        v0(1,1)=-0.2734d0
        v0(1,2)=-0.0359d0
        v0(1,3)=-0.2495d0

        v0(2,1)=-0.2243d0
        v0(2,2)=-0.0108d0
        v0(2,3)=-0.3548d0
!!!!!!!!!!!!!
        v0(3,1)=-0.1393d0
        v0(3,2)=-0.1971d0
        v0(3,3)=-0.4208d0
!!!!!!!!!!!!
        v0(4,1)=-0.5648d0
        v0(4,2)=-0.5648d0
        v0(4,3)=-0.5648d0

        else if (ii .eq. 2) then

        v0(1,1)=-0.1853d0
        v0(1,2)=-0.1269d0
        v0(1,3)=-0.3474d0

        v0(2,1)=-0.1265d0
        v0(2,2)=-0.1550d0
        v0(2,3)=-0.2869d0
!!!!!!!!!!!
        v0(3,1)=-0.1910d0
        v0(3,2)=-0.0865d0
        v0(3,3)=-0.0470d0
!!!!!!!!!!!
        v0(4,1)=-0.0481d0
        v0(4,2)=-0.0481d0
        v0(4,3)=-0.0481d0

        else if (ii .eq. 3) then

        v0(1,1)=-0.4000d0
        v0(1,2)=-0.2125d0
        v0(1,3)=-0.0093d0
        
        v0(2,1)=-0.1812d0
        v0(2,2)=-0.0606d0
        v0(2,3)=-0.0090d0
!!!!!!!!!!!!
        v0(3,1)=-5.0661d0
        v0(3,2)=-3.5468d0
        v0(3,3)=-0.0369d0
!!!!!!!!!!!
        v0(4,1)=-0.0000d0
        v0(4,2)=-0.0000d0
        v0(4,3)=-0.0000d0

        endif

       end

       subroutine getpara(ii,nrap,ntap,
     & KNHB1,KNHB2,
     & KNC1,KNC2,
     & KCH1,KCH2,
     & KCH3,
     & KCNHA1,KCNHA2,
     & KCNHB1,KCNHB2,
     & KNCH1,KNCH2,
     & JNHB1,JNHB2,
     & JNHB3,
     & JNC1,JNC2,JNC3,JCNHA1,JCNHA2,
     & JCNHA3,JCNHB1,JCNHB2,JCNHB3,
     & GNHB1,GNHB2,GNHB3,
     & GNC1,GNC2,GNC3,
     & ft1,fa1)
! subroutine gives parameters for the tertiary contribution
       implicit none
       integer ii
       integer nrap,i,j,ntap
      double precision KNHB1(nrap,ntap),KNHB2(nrap,ntap)
      double precision KNC1(nrap,ntap),KNC2(nrap,ntap)
      double precision KCH1(nrap,ntap),KCH2(nrap,ntap)
      double precision KCH3(nrap,ntap)
      double precision KCNHA1(nrap,ntap),KCNHA2(nrap,ntap)
      double precision KCNHB1(nrap,ntap),KCNHB2(nrap,ntap)
      double precision KNCH1(nrap,ntap),KNCH2(nrap,ntap)
!
      double precision JNHB1(nrap,ntap),JNHB2(nrap,ntap)
      double precision JNHB3(nrap,ntap)
      double precision JNC1(nrap,ntap),JNC2(nrap,ntap)
      double precision JNC3(nrap,ntap)
      double precision JCNHA1(nrap,ntap),JCNHA2(nrap,ntap)
      double precision JCNHA3(nrap,ntap)
      double precision JCNHB1(nrap,ntap),JCNHB2(nrap,ntap)
      double precision JCNHB3(nrap,ntap)
!
      double precision GNHB1(nrap,ntap),GNHB2(nrap,ntap)
      double precision GNHB3(nrap,ntap)
!
      double precision GNC1(nrap,ntap),GNC2(nrap,ntap)
      double precision GNC3(nrap,ntap)
      double precision ft1(nrap,ntap),fa1(nrap,ntap)

      KNHB1=0d0
      KNHB2=0d0
      KNC1=0d0
      KNC2=0d0
      KCH1=0d0
      KCH2=0d0
      KCH3=0d0
      KCNHA1=0d0
      KCNHA2=0d0
      KCNHB1=0d0
      KCNHB2=0d0
      KNCH1=0d0
      KNCH2=0d0
      JNHB1=0d0
      JNHB2=0d0
      JNHB3=0d0
      JNC1=0d0
      JNC2=0d0
      JNC3=0d0
      JCNHA1=0d0
      JCNHA2=0d0
      JCNHA3=0d0
      JCNHB1=0d0
      JCNHB2=0d0
      JCNHB3=0d0



       if (ii .eq. 1) then
! STRETCH 1 (NHB)
           GNHB1(1,2)=6.420621d0
           GNHB1(2,2)=6.420621d0
           GNHB1(3,2)=6.420621d0
           GNHB1(4,2)=6.420621d0
           GNHB1(1,3)=6.420621d0
           GNHB1(2,3)=6.420621d0
           GNHB1(3,3)=6.420621d0
           GNHB1(4,3)=GNHB1(4,2)
           GNHB1(1,1)=4.3d0
           GNHB1(2,1)=4.4d0
           GNHB1(3,1)=8d0
           GNHB1(4,1)=GNHB1(4,2)

           GNHB2(1,2)=-1.85d0
           GNHB2(2,2)=-1.85d0
           GNHB2(3,2)=-1.85d0
           GNHB2(4,2)=-1.96d0
           GNHB2(1,3)=-1.85d0
           GNHB2(2,3)=-1.85d0
           GNHB2(3,3)=-1.85d0
           GNHB2(4,3)=GNHB2(4,2)
           GNHB2(1,1)=-1.85d0
           GNHB2(2,1)=-1.85d0
           GNHB2(3,1)=-2.0d0
           GNHB2(4,1)=GNHB2(4,2)

           GNHB3(1,2)=1.01574847d0
           GNHB3(2,2)=1.01574847d0
           GNHB3(3,2)=1.02d0
           GNHB3(4,2)=1.00368d0
           GNHB3(1,3)=1.01574847d0
           GNHB3(2,3)=1.01574847d0
           GNHB3(3,3)=1.02d0
           GNHB3(4,3)=GNHB3(4,2)
           GNHB3(1,1)=1.08d0
           GNHB3(2,1)=1.07d0
           GNHB3(3,1)=1.0d0
           GNHB3(4,1)=GNHB3(4,2)

           KNHB1(1,2)=1.0d0
           KNHB1(2,2)=1.0d0
           KNHB1(3,2)=0.5d0
           KNHB1(4,2)=0.5d0
           KNHB1(1,3)=6.0d0
           KNHB1(2,3)=3.0d0
           KNHB1(3,3)=1.0d0
           KNHB1(4,3)=KNHB1(4,2)
           KNHB1(1,1)=1.0d0
           KNHB1(2,1)=1.0d0
           KNHB1(3,1)=0.5d0
           KNHB1(4,1)=KNHB1(4,2)

           KNHB2(1,2)=1.01574847d0
           KNHB2(2,2)=1.01574847d0
           KNHB2(3,2)=1.02d0
           KNHB2(4,2)=1.00368d0
           KNHB2(1,3)=1.01574847d0
           KNHB2(2,3)=1.01574847d0
           KNHB2(3,3)=1.02d0
           KNHB2(4,3)=KNHB2(4,2)
           KNHB2(1,1)=1.08d0
           KNHB2(2,1)=1.07d0
           KNHB2(3,1)=1.0d0
           KNHB2(4,1)=KNHB2(4,2)

! STRETCH 2 (NC)
           GNC1(1,2)=4.584032444d0
           GNC1(2,2)=4.543215344d0
           GNC1(3,2)=4.0d0
           GNC1(4,2)=2.7d0
           GNC1(1,3)=4.43273706d0
           GNC1(2,3)=4.412600624d0
           GNC1(3,3)=3.1d0
           GNC1(4,3)=GNC1(4,2)
           GNC1(1,1)=4.90621542d0
           GNC1(2,1)=4.870024258d0
           GNC1(3,1)=4.0d0
           GNC1(4,1)=GNC1(4,2)

           GNC2(1,2)=-2.0d0
           GNC2(2,2)=-2.0d0
           GNC2(3,2)=-1.9d0
           GNC2(4,2)=-2.6d0
           GNC2(1,3)=-2.0d0
           GNC2(2,3)=-2.0d0
           GNC2(3,3)=-2.3d0
           GNC2(4,3)=GNC2(4,2)
           GNC2(1,1)=-1.9d0
           GNC2(2,1)=-1.9d0
           GNC2(3,1)=-2.3d0
           GNC2(4,1)=GNC2(4,2)

           GNC3(1,2)=1.4597436d0
           GNC3(2,2)=1.4597436d0
           GNC3(3,2)=1.4597436d0
           GNC3(4,2)=1.40443d0
           GNC3(1,3)=1.4597436d0
           GNC3(2,3)=1.4597436d0
           GNC3(3,3)=1.44d0
           GNC3(4,3)=GNC3(4,2)
           GNC3(1,1)=1.4597436d0
           GNC3(2,1)=1.4597436d0
           GNC3(3,1)=1.45d0
           GNC3(4,1)=GNC3(4,2)

           KNC1(1,2)=3.0d0
           KNC1(2,2)=5.0d0
           KNC1(3,2)=8.0d0
           KNC1(4,2)=1.7d0
           KNC1(1,3)=2.0d0
           KNC1(2,3)=3.0d0
           KNC1(3,3)=4.5d0
           KNC1(4,3)=KNC1(4,2)
           KNC1(1,1)=4.0d0
           KNC1(2,1)=5.0d0
           KNC1(3,1)=4.0d0
           KNC1(4,1)=KNC1(4,2)

           KNC2(1,2)=1.4597436d0
           KNC2(2,2)=1.4597436d0
           KNC2(3,2)=1.4597436d0
           KNC2(4,2)=1.40443d0
           KNC2(1,3)=1.4597436d0
           KNC2(2,3)=1.4597436d0
           KNC2(3,3)=1.44d0
           KNC2(4,3)=KNC2(4,2)
           KNC2(1,1)=1.4597436d0
           KNC2(2,1)=1.4597436d0
           KNC2(3,1)=1.45d0
           KNC2(4,1)=KNC2(4,2)

! STRETCHES 3-5 (CH)
           KCH1(1,2)=4.680905028d0
           KCH1(2,2)=4.660046401d0
           KCH1(3,2)=3.103188056d0
           KCH1(4,2)=2.783998334d0
           KCH1(1,3)=4.778593954d0
           KCH1(2,3)=4.763899798d0
           KCH1(3,3)=3.552720384d0
           KCH1(4,3)=KCH1(4,2)
           KCH1(1,1)=4.655054198d0
           KCH1(2,1)=4.633285078d0
           KCH1(3,1)=2.861006596d0
           KCH1(4,1)=KCH1(4,2)

           KCH2(1,2)=-2.0d0
           KCH2(2,2)=-1.9d0
           KCH2(3,2)=-2.3d0
           KCH2(4,2)=-2.3d0
           KCH2(1,3)=-2.0d0
           KCH2(2,3)=-2.0d0
           KCH2(3,3)=-2.2d0
           KCH2(4,3)=KCH2(4,2)
           KCH2(1,1)=-2.0d0
           KCH2(2,1)=-2.0d0
           KCH2(3,1)=-2.4d0
           KCH2(4,1)=KCH2(4,2)

           KCH3(1,2)=1.0894707d0
           KCH3(2,2)=1.0894707d0
           KCH3(3,2)=1.0894707d0
           KCH3(4,2)=1.094d0
           KCH3(1,3)=1.0894707d0
           KCH3(2,3)=1.0894707d0
           KCH3(3,3)=1.0894707d0
           KCH3(4,3)=KCH3(4,2)
           KCH3(1,1)=1.0894707d0
           KCH3(2,1)=1.0894707d0
           KCH3(3,1)=1.0894707d0
           KCH3(4,1)=KCH3(4,2)

! BEND 1 (CNHA)
           KCNHA1(1,2)=3.7d0
           KCNHA1(2,2)=3.4d0
           KCNHA1(3,2)=1.0d0
           KCNHA1(4,2)=0.0d0
           KCNHA1(1,3)=3.7d0
           KCNHA1(2,3)=3.8d0
           KCNHA1(3,3)=1.5d0
           KCNHA1(4,3)=KCNHA1(4,2)
           KCNHA1(1,1)=3.0d0
           KCNHA1(2,1)=3.3d0
           KCNHA1(3,1)=2.0d0
           KCNHA1(4,1)=KCNHA1(4,2)

           KCNHA2(1,2)=-0.354096319d0
           KCNHA2(2,2)=-0.354096319d0
           KCNHA2(3,2)=-0.354096319d0
           KCNHA2(4,2)=0.0d0
           KCNHA2(1,3)=-0.544639035d0
           KCNHA2(2,3)=-0.559192903d0
           KCNHA2(3,3)=-0.544639035d0
           KCNHA2(4,3)=KCNHA2(4,2)
           KCNHA2(1,1)=-0.207911691d0
           KCNHA2(2,1)=-0.207911691d0
           KCNHA2(3,1)=-0.258819045d0
           KCNHA2(4,1)=KCNHA2(4,2)

! BEND 2 (CNHB)
           KCNHB1(1,2)=3.9d0
           KCNHB1(2,2)=3.4d0
           KCNHB1(3,2)=3.3d0
           KCNHB1(4,2)=2.6d0
           KCNHB1(1,3)=3.7d0
           KCNHB1(2,3)=3.3d0
           KCNHB1(3,3)=4.0d0
           KCNHB1(4,3)=KCNHB1(4,2)
           KCNHB1(1,1)=3.4d0
           KCNHB1(2,1)=3.5d0
           KCNHB1(3,1)=4.0d0
           KCNHB1(4,1)=KCNHB1(4,2)

           KCNHB2(1,2)=-0.325568154d0
           KCNHB2(2,2)=-0.355262388d0
           KCNHB2(3,2)=-0.573576436d0
           KCNHB2(4,2)=-0.788976538d0
           KCNHB2(1,3)=-0.5d0
           KCNHB2(2,3)=-0.573576436d0
           KCNHB2(3,3)=-0.64278761d0
           KCNHB2(4,3)=KCNHB2(4,2)
           KCNHB2(1,1)=-0.173648178d0
           KCNHB2(2,1)=-0.207911691d0
           KCNHB2(3,1)=-0.5d0
           KCNHB2(4,1)=KCNHB2(4,2)

! BEND 3-5 (NCH)
           KNCH1(1,2)=3.7d0
           KNCH1(2,2)=3.4d0
           KNCH1(3,2)=3.7d0
           KNCH1(4,2)=3.0d0
           KNCH1(1,3)=3.5d0
           KNCH1(2,3)=3.5d0
           KNCH1(3,3)=3.6d0
           KNCH1(4,3)=KNCH1(4,2)
           KNCH1(1,1)=3.5d0
           KNCH1(2,1)=3.6d0
           KNCH1(3,1)=3.5d0
           KNCH1(4,1)=KNCH1(4,2)

           KNCH2(1,2)=1.92d0
           KNCH2(2,2)=1.92d0
           KNCH2(3,2)=1.97d0
           KNCH2(4,2)=1.93d0
           KNCH2(1,3)=1.97d0
           KNCH2(2,3)=1.97d0
           KNCH2(3,3)=1.95d0
           KNCH2(4,3)=KNCH2(4,2)
           KNCH2(1,1)=1.95d0
           KNCH2(2,1)=1.93d0
           KNCH2(3,1)=1.93d0
           KNCH2(4,1)=KNCH2(4,2)

! Methylrot and Methylang
           ft1(1,2)=0.0305d0
           ft1(2,2)=0.0305d0
           ft1(3,2)=0.0035d0
           ft1(4,2)=0.0035d0
           ft1(1,3)=0.0305d0
           ft1(2,3)=0.0305d0
           ft1(3,3)=ft1(3,2)
           ft1(4,3)=ft1(4,2)
           ft1(1,1)=0.0305d0
           ft1(2,1)=0.0305d0
           ft1(3,1)=ft1(3,2)
           ft1(4,1)=ft1(4,2)

           fa1(1,2)=0.5d0
           fa1(2,2)=0.5d0
           fa1(3,2)=1.5d0
           fa1(4,2)=1.5d0
           fa1(1,3)=0.5d0
           fa1(2,3)=0.5d0
           fa1(3,3)=fa1(3,2)
           fa1(4,3)=fa1(4,2)
           fa1(1,1)=0.5d0
           fa1(2,1)=0.5d0
           fa1(3,1)=fa1(3,2)
           fa1(4,1)=fa1(4,2)


       else if (ii .eq. 2) then
! STRETCH 1 (NHB)

           GNHB1(1,2)=1.0d0
           GNHB1(2,2)=1.0d0
           GNHB1(3,2)=1.5d0
           GNHB1(4,2)=3.2d0
           GNHB1(1,3)=1.0d0
           GNHB1(2,3)=1.0d0
           GNHB1(3,3)=1.8d0
           GNHB1(4,3)=GNHB1(4,2)
           GNHB1(1,1)=1.0d0
           GNHB1(2,1)=1.0d0
           GNHB1(3,1)=0.8d0
           GNHB1(4,1)=GNHB1(4,2)

           GNHB2(1,2)=-3.8d0
           GNHB2(2,2)=-3.8d0
           GNHB2(3,2)=-3.8d0
           GNHB2(4,2)=-2.5d0
           GNHB2(1,3)=-3.5d0
           GNHB2(2,3)=-4.2d0
           GNHB2(3,3)=-3.6d0
           GNHB2(4,3)=GNHB2(4,2)
           GNHB2(1,1)=-3.1d0
           GNHB2(2,1)=-3.3d0
           GNHB2(3,1)=-3.3d0
           GNHB2(4,1)=GNHB2(4,2)

           GNHB3(1,2)=1.040695848d0
           GNHB3(2,2)=1.040695848d0
           GNHB3(3,2)=1.0d0
           GNHB3(4,2)=1.030034236d0
           GNHB3(1,3)=1.040695848d0
           GNHB3(2,3)=1.055d0
           GNHB3(3,3)=1.01d0
           GNHB3(4,3)=GNHB3(4,2)
           GNHB3(1,1)=1.11d0
           GNHB3(2,1)=1.08d0
           GNHB3(3,1)=1.08d0
           GNHB3(4,1)=GNHB3(4,2)

           KNHB1(1,2)=2.0d0
           KNHB1(2,2)=2.0d0
           KNHB1(3,2)=3.0d0
           KNHB1(4,2)=11.0d0
           KNHB1(1,3)=2.0d0
           KNHB1(2,3)=2.0d0
           KNHB1(3,3)=6.0d0
           KNHB1(4,3)=KNHB1(4,2)
           KNHB1(1,1)=2.0d0
           KNHB1(2,1)=2.0d0
           KNHB1(3,1)=3.0d0
           KNHB1(4,1)=KNHB1(4,2)

           KNHB2(1,2)=1.040695848d0
           KNHB2(2,2)=1.040695848d0
           KNHB2(3,2)=1.0d0
           KNHB2(4,2)=1.030034236d0
           KNHB2(1,3)=1.040695848d0
           KNHB2(2,3)=1.055d0
           KNHB2(3,3)=1.01d0
           KNHB2(4,3)=KNHB2(4,2)
           KNHB2(1,1)=1.11d0
           KNHB2(2,1)=1.08d0
           KNHB2(3,1)=1.08d0
           KNHB2(4,1)=KNHB2(4,2)

! STRETCH 2 (NC)

           GNC1(1,2)=3.0d0
           GNC1(2,2)=3.0d0
           GNC1(3,2)=3.0d0
           GNC1(4,2)=3.34292049d0
           GNC1(1,3)=3.0d0
           GNC1(2,3)=3.0d0
           GNC1(3,3)=3.8d0
           GNC1(4,3)=GNC1(4,2)
           GNC1(1,1)=3.0d0
           GNC1(2,1)=3.0d0
           GNC1(3,1)=4.0d0
           GNC1(4,1)=GNC1(4,2)
           
           GNC2(1,2)=-2.3d0
           GNC2(2,2)=-2.3d0
           GNC2(3,2)=-2.3d0
           GNC2(4,2)=-2.1d0
           GNC2(1,3)=-2.3d0
           GNC2(2,3)=-2.4d0
           GNC2(3,3)=-2.0d0
           GNC2(4,3)=GNC2(4,2)
           GNC2(1,1)=-2.3d0
           GNC2(2,1)=-2.3d0
           GNC2(3,1)=-1.9d0
           GNC2(4,1)=GNC2(4,2)
           
           GNC3(1,2)=1.43081d0
           GNC3(2,2)=1.43081d0
           GNC3(3,2)=1.43081d0
           GNC3(4,2)=1.44623d0
           GNC3(1,3)=1.43081d0
           GNC3(2,3)=1.43081d0
           GNC3(3,3)=1.445d0
           GNC3(4,3)=GNC3(4,2)
           GNC3(1,1)=1.43081d0
           GNC3(2,1)=1.43081d0
           GNC3(3,1)=1.43081d0
           GNC3(4,1)=GNC3(4,2)

           KNC1(1,2)=2.0d0
           KNC1(2,2)=2.0d0
           KNC1(3,2)=2.0d0
           KNC1(4,2)=2.0d0
           KNC1(1,3)=2.0d0
           KNC1(2,3)=2.0d0
           KNC1(3,3)=2.5d0
           KNC1(4,3)=KNC1(4,2)
           KNC1(1,1)=2.0d0
           KNC1(2,1)=2.0d0
           KNC1(3,1)=1.0d0
           KNC1(4,1)=KNC1(4,2)

           KNC2(1,2)=1.43081d0
           KNC2(2,2)=1.43081d0
           KNC2(3,2)=1.43081d0
           KNC2(4,2)=1.44623d0
           KNC2(1,3)=1.43081d0
           KNC2(2,3)=1.43081d0
           KNC2(3,3)=1.445d0
           KNC2(4,3)=KNC2(4,2)
           KNC2(1,1)=1.43081d0
           KNC2(2,1)=1.43081d0
           KNC2(3,1)=1.43081d0
           KNC2(4,1)=KNC2(4,2)

! STRETCHES 3-5 (CH)
           KCH1(1,2)=2.811209734d0
           KCH1(2,2)=3.037064354d0
           KCH1(3,2)=4.641992726d0
           KCH1(4,2)=4.605257336d0
           KCH1(1,3)=2.164122642d0
           KCH1(2,3)=2.476781628d0
           KCH1(3,3)=4.53750095d0
           KCH1(4,3)=KCH1(4,2)
           KCH1(1,1)=3.084140076d0
           KCH1(2,1)=3.237068144d0
           KCH1(3,1)=4.646074436d0
           KCH1(4,1)=KCH1(4,2)

           KCH2(1,2)=-2.4d0
           KCH2(2,2)=-2.4d0
           KCH2(3,2)=-2.0d0
           KCH2(4,2)=-1.9d0
           KCH2(1,3)=-2.6d0
           KCH2(2,3)=-2.5d0
           KCH2(3,3)=-2.0d0
           KCH2(4,3)=KCH2(4,2)
           KCH2(1,1)=-2.3d0
           KCH2(2,1)=-2.3d0
           KCH2(3,1)=-2.0d0
           KCH2(4,1)=KCH2(4,2)

           KCH3(1,2)=1.0894707d0
           KCH3(2,2)=1.0894707d0
           KCH3(3,2)=1.0894707d0
           KCH3(4,2)=1.0894707d0
           KCH3(1,3)=1.0894707d0
           KCH3(2,3)=1.0894707d0
           KCH3(3,3)=1.0894707d0
           KCH3(4,3)=KCH3(4,2)
           KCH3(1,1)=1.0894707d0
           KCH3(2,1)=1.0894707d0
           KCH3(3,1)=1.0894707d0
           KCH3(4,1)=KCH3(4,2)

! BEND 1 (CNHA)
           KCNHA1(1,2)=3.2d0
           KCNHA1(2,2)=3.0d0
           KCNHA1(3,2)=1.0d0
           KCNHA1(4,2)=0.0d0
           KCNHA1(1,3)=3.6d0
           KCNHA1(2,3)=4.0d0
           KCNHA1(3,3)=0.5d0
           KCNHA1(4,3)=KCNHA1(4,2)
           KCNHA1(1,1)=1.9d0
           KCNHA1(2,1)=2.1d0
           KCNHA1(3,1)=1.0d0
           KCNHA1(4,1)=KCNHA1(4,2)

           KCNHA2(1,2)=-0.4539905d0
           KCNHA2(2,2)=-0.48480962d0
           KCNHA2(3,2)=-0.48480962d0
           KCNHA2(4,2)=0d0
           KCNHA2(1,3)=-0.559192903d0
           KCNHA2(2,3)=-0.505340076d0
           KCNHA2(3,3)=-0.5d0
           KCNHA2(4,3)=KCNHA2(4,2)
           KCNHA2(1,1)=-0.438371147d0
           KCNHA2(2,1)=-0.4539905d0
           KCNHA2(3,1)=-0.64278761d0
           KCNHA2(4,1)=KCNHA2(4,2)

! BEND 2 (CNHB)
           KCNHB1(1,2)=3.7d0
           KCNHB1(2,2)=3.3d0
           KCNHB1(3,2)=2.8d0
           KCNHB1(4,2)=3.3d0
           KCNHB1(1,3)=3.6d0
           KCNHB1(2,3)=3.3d0
           KCNHB1(3,3)=2.5d0
           KCNHB1(4,3)=KCNHB1(4,2)
           KCNHB1(1,1)=2.3d0
           KCNHB1(2,1)=2.1d0
           KCNHB1(3,1)=3.0d0
           KCNHB1(4,1)=KCNHB1(4,2)

           KCNHB2(1,2)=-0.406736643d0
           KCNHB2(2,2)=-0.422618262d0
           KCNHB2(3,2)=-0.406736643d0
           KCNHB2(4,2)=-0.286195791d0
           KCNHB2(1,3)=-0.529919264d0
           KCNHB2(2,3)=-0.528778795d0
           KCNHB2(3,3)=-0.342020143d0
           KCNHB2(4,3)=KCNHB2(4,2)
           KCNHB2(1,1)=-0.342020143d0
           KCNHB2(2,1)=-0.374606593d0
           KCNHB2(3,1)=-0.422618262d0
           KCNHB2(4,1)=KCNHB2(4,2)

! BEND 3-5 (NCH)
           KNCH1(1,2)=3.5d0
           KNCH1(2,2)=3.7d0
           KNCH1(3,2)=3.0d0
           KNCH1(4,2)=2.9d0
           KNCH1(1,3)=3.0d0
           KNCH1(2,3)=2.8d0
           KNCH1(3,3)=2.6d0
           KNCH1(4,3)=KNCH1(4,2)
           KNCH1(1,1)=3.6d0
           KNCH1(2,1)=3.6d0
           KNCH1(3,1)=3.1d0
           KNCH1(4,1)=KNCH1(4,2)

           KNCH2(1,2)=1.88d0
           KNCH2(2,2)=1.88d0
           KNCH2(3,2)=1.89d0
           KNCH2(4,2)=1.92d0
           KNCH2(1,3)=1.87d0
           KNCH2(2,3)=1.88d0
           KNCH2(3,3)=1.9d0
           KNCH2(4,3)=KNCH2(4,2)
           KNCH2(1,1)=1.90d0
           KNCH2(2,1)=1.91d0
           KNCH2(3,1)=1.9d0
           KNCH2(4,1)=KNCH2(4,2)


! Methylrot and Methylang
           ft1(1,2)=0.012d0
           ft1(2,2)=0.012d0
           ft1(3,2)=0.012d0
           ft1(4,2)=0.012d0
           ft1(1,3)=0.012d0
           ft1(2,3)=0.012d0
           ft1(3,3)=ft1(3,2)
           ft1(4,3)=ft1(4,2)
           ft1(1,1)=0.012d0
           ft1(2,1)=0.012d0
           ft1(3,1)=ft1(3,2)
           ft1(4,1)=ft1(4,2)

           fa1(1,2)=0.5d0
           fa1(2,2)=0.5d0
           fa1(3,2)=0.8d0
           fa1(4,2)=0.8d0
           fa1(1,3)=0.5d0
           fa1(2,3)=0.5d0
           fa1(3,3)=fa1(3,2)
           fa1(4,3)=fa1(4,2)
           fa1(1,1)=0.5d0
           fa1(2,1)=0.5d0
           fa1(3,1)=fa1(3,2)
           fa1(4,1)=fa1(4,2)



       else if (ii .eq. 3) then
! STRETCH 1 (NHB)
           JNHB1(1,2)=1.4d0
           JNHB1(2,2)=1.4d0
           JNHB1(3,2)=0.9d0
           JNHB1(4,2)=0d0
           JNHB1(1,3)=0.28d0
           JNHB1(2,3)=0.28d0
           JNHB1(3,3)=1.15d0
           JNHB1(4,3)=JNHB1(4,2)
           JNHB1(1,1)=1.4d0
           JNHB1(2,1)=1.44d0
           JNHB1(3,1)=1.28d0
           JNHB1(4,1)=JNHB1(4,2)

           JNHB2(1,2)=-8.5d0
           JNHB2(2,2)=-10.0d0
           JNHB2(3,2)=-0.4d0
           JNHB2(4,2)=-1d0
           JNHB2(1,3)=-22d0
           JNHB2(2,3)=-22d0
           JNHB2(3,3)=-35d0
           JNHB2(4,3)=JNHB2(4,2)
           JNHB2(1,1)=-6.5d0
           JNHB2(2,1)=-6.5d0
           JNHB2(3,1)=-0.5d0
           JNHB2(4,1)=JNHB2(4,2)

           JNHB3(1,2)=1.8145815d0
           JNHB3(2,2)=1.8145815d0
           JNHB3(3,2)=1.2d0
           JNHB3(4,2)=1d0
           JNHB3(1,3)=1.7145815d0
           JNHB3(2,3)=1.7145815d0
           JNHB3(3,3)=2d0
           JNHB3(4,3)=JNHB3(4,2)
           JNHB3(1,1)=2.0145815d0
           JNHB3(2,1)=2.0145815d0
           JNHB3(3,1)=1.1d0
           JNHB3(4,1)=JNHB3(4,2)

! STRETCH 2 (NC)
           JNC1(1,2)=1.65d0
           JNC1(2,2)=1.6d0
           JNC1(3,2)=0.95
           JNC1(4,2)=0d0
           JNC1(1,3)=0.32d0
           JNC1(2,3)=0.4d0
           JNC1(3,3)=0.34d0
           JNC1(4,3)=JNC1(4,2)
           JNC1(1,1)=2.4d0
           JNC1(2,1)=2.5d0
           JNC1(3,1)=1.3d0
           JNC1(4,1)=JNC1(4,2)

           JNC2(1,2)=-25.0d0
           JNC2(2,2)=-25.0d0
           JNC2(3,2)=-0.7d0
           JNC2(4,2)=-1
           JNC2(1,3)=-30.0d0
           JNC2(2,3)=-30.0d0
           JNC2(3,3)=-30.0d0
           JNC2(4,3)=JNC2(4,2)
           JNC2(1,1)=-30.0d0
           JNC2(2,1)=-30.0d0
           JNC2(3,1)=-1.2d0
           JNC2(4,1)=JNC2(4,2)

           JNC3(1,2)=1.9597436d0
           JNC3(2,2)=2.0597436d0
           JNC3(3,2)=1.7d0
           JNC3(4,2)=1d0
           JNC3(1,3)=2.3597436d0
           JNC3(2,3)=2.3597436d0
           JNC3(3,3)=2.6d0
           JNC3(4,3)=JNC3(4,2)
           JNC3(1,1)=1.9597436d0
           JNC3(2,1)=1.9597436d0
           JNC3(3,1)=1.55d0
           JNC3(4,1)=JNC3(4,2)

! BEND 1 (CNHA)
           JCNHA1(1,2)=0.1d0
           JCNHA1(2,2)=0.06d0
           JCNHA1(3,2)=0.93d0
           JCNHA1(4,2)=0d0
           JCNHA1(1,3)=0.007d0
           JCNHA1(2,3)=0.005d0
           JCNHA1(3,3)=0.048d0
           JCNHA1(4,3)=JCNHA1(4,2)
           JCNHA1(1,1)=0.2d0
           JCNHA1(2,1)=0.18d0
           JCNHA1(3,1)=1.3d0
           JCNHA1(4,1)=JCNHA1(4,2)

           JCNHA2(1,2)=-2.5d0
           JCNHA2(2,2)=-3.0d0
           JCNHA2(3,2)=-1.5d0
           JCNHA2(4,2)=-1d0
           JCNHA2(1,3)=-1.2d0
           JCNHA2(2,3)=-0.5d0
           JCNHA2(3,3)=-3.0d0
           JCNHA2(4,3)=JCNHA2(4,2)
           JCNHA2(1,1)=-0.7d0
           JCNHA2(2,1)=-2.3d0
           JCNHA2(3,1)=-2.0d0
           JCNHA2(4,1)=JCNHA2(4,2)

           JCNHA3(1,2)=-0.354096319d0
           JCNHA3(2,2)=0.160948072d0
           JCNHA3(3,2)=-0.173648178
           JCNHA3(4,2)=1
           JCNHA3(1,3)=0.160948072d0
           JCNHA3(2,3)=-0.354096319d0
           JCNHA3(3,3)=0.329887214
           JCNHA3(4,3)=JCNHA3(4,2)
           JCNHA3(1,1)=-0.186319468d0
           JCNHA3(2,1)=0.160948072d0
           JCNHA3(3,1)=-0.258819045
           JCNHA3(4,1)=JCNHA3(4,2)

! BEND 2 (CNHB)
           JCNHB1(1,2)=0.15d0
           JCNHB1(2,2)=0.06d0
           JCNHB1(3,2)=1.4d0
           JCNHB1(4,2)=0
           JCNHB1(1,3)=0.006d0
           JCNHB1(2,3)=0.005d0
           JCNHB1(3,3)=0.04d0
           JCNHB1(4,3)=JCNHB1(4,2)
           JCNHB1(1,1)=0.32d0
           JCNHB1(2,1)=0.14d0
           JCNHB1(3,1)=1.9d0
           JCNHB1(4,1)=JCNHB1(4,2)

           JCNHB2(1,2)=-0.8d0
           JCNHB2(2,2)=-2.5d0
           JCNHB2(3,2)=-1d0
           JCNHB2(4,2)=-1
           JCNHB2(1,3)=-3.0d0
           JCNHB2(2,3)=-0.5d0
           JCNHB2(3,3)=-1.0d0
           JCNHB2(4,3)=JCNHB2(4,2)
           JCNHB2(1,1)=-1.0d0
           JCNHB2(2,1)=-1.2d0
           JCNHB2(3,1)=-0.9d0
           JCNHB2(4,1)=JCNHB2(4,2)

           JCNHB3(1,2)=0.328709624d0
           JCNHB3(2,2)=0.159717057d0
           JCNHB3(3,2)=0.342020143d0
           JCNHB3(4,2)=1d0
           JCNHB3(1,3)=-0.014128431d0
           JCNHB3(2,3)=0.328709624d0
           JCNHB3(3,3)=0.328709624d0
           JCNHB3(4,3)=JCNHB3(4,2)
           JCNHB3(1,1)=0.328709624d0
           JCNHB3(2,1)=0.328709624d0
           JCNHB3(3,1)=0.342020143d0
           JCNHB3(4,1)=JCNHB3(4,2)

       endif






       end


      subroutine tert(igrad,ii,X,Y,Z,nrap,ntap,EN,dx1,dy1,dz1)
!***********************************************************************
! Subroutine to calculate the tertiary contribution
!***********************************************************************
      implicit none
      integer tk
      integer nrap,ntap
      integer igrad
      double precision r1,t2
      double precision X(7),Y(7),Z(7)
      double precision xyz(4,3)
      double precision dx(22)
      double precision rap(nrap,3),tap(ntap,3)
      integer i,ii,j
      integer fg(2)
      double precision u,v,f(4)
      double precision bndlen,bndang
      double precision H(4),EN
      double precision dH(4,20),dHe(20)
      double precision df(4,2)
      double precision dTr(2)
      double precision dx1(7),dy1(7),dz1(7)
      double precision beta
      double precision dbeta(12)
      double precision dihedr2


! The anchor points:

      do i=1,3
      rap(1,i)=0.7157485d0
      rap(2,i)=1.0157485d0
      rap(3,i)=2.2157485d0
      rap(4,i)=3.5157485d0
      enddo
      do i=1,3
      tap(1,i)=-0.193287305d0
      tap(2,i)=0.482599384d0
      tap(3,i)=1.0d0
      enddo

      r1=bndlen(1,3,X,Y,Z)
      t2=dihedr2(1,3,4,2,X,Y,Z)

! calls the 2D tent
      call bigtent(r1,t2,f,rap,nrap,tap,ntap,fg,ii,tk,df)

! calcs the contribution at the anchor points
       H=0d0
       if (tk.eq.1) then
       call calcHess1(igrad,nrap,ntap,fg,X,Y,Z,ii,H,dH)
       else if (tk.eq.2) then
       call calcHess2(igrad,nrap,ntap,fg,X,Y,Z,ii,H,dH)
       else if (tk.eq.3) then
       call calcHess3(igrad,nrap,ntap,fg,X,Y,Z,ii,H,dH)
       endif

! Combines the tent functions and contribution from the anchor points
       EN=0d0
       do i=1,4
       EN=EN+H(i)*f(i)
       enddo

       dx1=0d0
       dy1=0d0
       dz1=0d0
!!!!!!!!!!! GET GRADIENT !!!!!!!!!
       if (igrad .eq. 1) then
       dHe=0d0
       do j=1,20
        do i=1,4
         dHe(j)=dHe(j)+dH(i,j)*f(i)
        enddo
       enddo

      dTr=0d0
      do j=1,2
       do i=1,4
        dTr(j)=dTr(j)+df(i,j)*H(i)
       enddo
      enddo

       do i=1,20
       dx(i)=dHe(i)
       enddo
       dx(21)=dTr(1)
       dx(22)=dTr(2)

       call convertit(X,Y,Z,dx,dx1,dy1,dz1)
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       end
!-------------------------------------------------!
!-------------------------------------------------!
!-------------------------------------------------!
!-------------------------------------------------!


       subroutine calcHess1(igrad,nrap,ntap,fg,X,Y,Z,ii,H,dH)
! calcs the tertiary contribution for case 1: corners of the 2D tent
! i.e. only need 1 anchor point
      implicit none
      integer igrad
      integer ra1,ra2
      integer nrap,ntap
      double precision r1,t2
      double precision X(7),Y(7),Z(7)
      double precision qc(20)
      double precision e
      double precision dx

      integer i,ii
      integer fg(2)
      double precision u,v,f(4)
      double precision H(4),Hs
      double precision dHs(20),dH(4,20)
      double precision PI

      PI=4.D0*DATAN(1.D0)
      Hs=0d0
      dHs=0d0
      dH=0d0

       call calcgeom(X,Y,Z,qc)
!****** CORNERS
       if( (fg(1)-1.eq.0) .and. (fg(2)-1.eq.0)) then
       ra1=1
       ra2=1
        call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(4)=Hs
         do i=1,20
          dH(4,i)=dHs(i)
         enddo
       else if( (fg(1).eq.nrap+1) .and. (fg(2).eq.ntap+1)) then
       ra1=nrap
       ra2=ntap
       call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(1)=Hs
         do i=1,20
          dH(1,i)=dHs(i)
         enddo
       else if( (fg(1)-1.eq.0) .and. (fg(2).eq.ntap+1)) then
       ra1=1
       ra2=ntap
       call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(3)=Hs
         do i=1,20
          dH(3,i)=dHs(i)
         enddo
       else if( (fg(1).eq.nrap+1) .and. (fg(2)-1.eq.0)) then
       ra1=nrap
       ra2=1
       call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(2)=Hs
         do i=1,20
          dH(2,i)=dHs(i)
         enddo
       endif


       end

      subroutine calcHess2(igrad,nrap,ntap,fg,X,Y,Z,ii,H,dH)
! calcs the tertiary contribution for case 2: sides of the 2D tent
! i.e. between 2 anchor points
      implicit none
      integer igrad
      integer ra1,ra2
      integer nrap,ntap  
      double precision r1,t2
      double precision X(7),Y(7),Z(7)
      double precision qc(20)
      double precision e
      double precision dx

      integer i,ii
      integer fg(2)
      double precision u,v,f(4)
      double precision H(4),Hs
      double precision dHs(20),dH(4,20)
      double precision PI

      PI=4.D0*DATAN(1.D0)
      Hs=0d0
      dHs=0d0
      dH=0d0

       call calcgeom(X,Y,Z,qc)
!****** SIDES
       if( (fg(1)-1.eq.0) ) then
       ra1=1
       ra2=fg(2)
       call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(3)=Hs
         do i=1,20
          dH(3,i)=dHs(i)
         enddo
       call solver(igrad,ra1,ra2-1,qc,ii,nrap,ntap,Hs,dHs)
       H(4)=Hs
         do i=1,20
          dH(4,i)=dHs(i)
         enddo

       else if( (fg(1).eq.nrap+1) ) then
       ra1=nrap
       ra2=fg(2)
       call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(1)=Hs
         do i=1,20
          dH(1,i)=dHs(i)
         enddo 
       call solver(igrad,ra1,ra2-1,qc,ii,nrap,ntap,Hs,dHs)
       H(2)=Hs
         do i=1,20
          dH(2,i)=dHs(i)
         enddo 

       else if( (fg(2)-1.eq.0) ) then
       ra2=1
       ra1=fg(1)
       call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(2)=Hs
         do i=1,20
          dH(2,i)=dHs(i)
         enddo 
       call solver(igrad,ra1-1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(4)=Hs
         do i=1,20
          dH(4,i)=dHs(i)
         enddo 

       else if( (fg(2).eq.ntap+1) ) then
       ra2=ntap
       ra1=fg(1)
       call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(1)=Hs
         do i=1,20
          dH(1,i)=dHs(i)
         enddo 
       call solver(igrad,ra1-1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(3)=Hs
         do i=1,20
          dH(3,i)=dHs(i)
         enddo 

       endif

       end


       subroutine calcHess3(igrad,nrap,ntap,fg,X,Y,Z,ii,H,dH)
! calcs the tertiary contribution for case 3: between 4 anchor points
      implicit none
      integer igrad
      integer ra1,ra2
      integer nrap,ntap 
      double precision r1,t2
      double precision X(7),Y(7),Z(7)
      double precision qc(20)
      double precision e
      double precision dx

      integer i,ii
      integer fg(2)
      double precision u,v,f(4)
      double precision H(4),Hs
      double precision dHs(20),dH(4,20)
      double precision PI

      PI=4.D0*DATAN(1.D0)
      Hs=0d0
      dHs=0d0
      dH=0d0

       ra1=fg(1)
       ra2=fg(2)

       call calcgeom(X,Y,Z,qc)
       call solver(igrad,ra1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(1)=Hs
         do i=1,20
          dH(1,i)=dHs(i)
         enddo
       call solver(igrad,ra1,ra2-1,qc,ii,nrap,ntap,Hs,dHs)
       H(2)=Hs 
         do i=1,20
          dH(2,i)=dHs(i)
         enddo
       call solver(igrad,ra1-1,ra2,qc,ii,nrap,ntap,Hs,dHs)
       H(3)=Hs
         do i=1,20
          dH(3,i)=dHs(i)
         enddo 
       call solver(igrad,ra1-1,ra2-1,qc,ii,nrap,ntap,Hs,dHs)
       H(4)=Hs
         do i=1,20
          dH(4,i)=dHs(i)
         enddo 

       end



      subroutine calcgeom(X,Y,Z,qc)
! get the internal coords for tertiary calcs
      implicit none
      integer i
      double precision X(7),Y(7),Z(7)
      double precision qc(20)
      double precision bndlen,bndang
      double precision dihedr2

      qc(1)=bndlen(2,3,X,Y,Z)
      qc(2)=bndlen(3,4,X,Y,Z)

      qc(3)=bndlen(4,5,X,Y,Z)
      qc(4)=bndlen(4,6,X,Y,Z)
      qc(5)=bndlen(4,7,X,Y,Z)

      qc(6)=bndang(1,3,4,X,Y,Z)
      qc(7)=bndang(2,3,4,X,Y,Z)

      qc(8)=bndang(3,4,5,X,Y,Z)
      qc(9)=bndang(3,4,6,X,Y,Z)
      qc(10)=bndang(3,4,7,X,Y,Z)

        qc(11)=dihedr2(1,3,4,5,X,Y,Z)
        qc(12)=dihedr2(2,3,4,5,X,Y,Z)

        qc(13)=bndlen(1,5,X,Y,Z)
        qc(14)=bndlen(2,5,X,Y,Z)

        qc(15)=bndang(5,4,6,X,Y,Z)
        qc(16)=bndang(6,4,7,X,Y,Z)
        qc(17)=bndang(7,4,5,X,Y,Z)

        qc(18)=bndlen(5,6,X,Y,Z)
        qc(19)=bndlen(6,7,X,Y,Z)
        qc(20)=bndlen(7,5,X,Y,Z)

      end


      subroutine bigtent(r1,t2,f,rap,nrap,tap,ntap,fg,si,tk,df)
!***********************************************************************
! Subroutine to evaluate 2D tent functions
!***********************************************************************
      implicit none
      integer i,si,s,tk
      integer nrap  !number of anchor points
      integer ntap  !number of anchor points
      double precision r1,t2

      double precision rap(nrap,3)
      double precision tap(ntap,3)
      double precision f(4)
      double precision un,up,vn,vp,dnom
      double precision df(4,2)
      double precision dun,dup,dvn,dvp,ddnom(2)
      integer fg(2)

      f=0d0
      fg=0
      


! Check whether the point are end points
      if (r1 .lt. rap(1,si)) then
        fg(1)=1
      else if (r1 .ge. rap(nrap,si)) then
        fg(1)=nrap+1
      endif

      do i=2,nrap
       if ((rap(i-1,si) .le. r1) .and. (r1 .lt. rap(i,si))) then
       fg(1)=i
       endif
      enddo

! Check whether the point are end points
      if (t2 .lt. tap(1,si)) then
        fg(2)=1
      else if (t2 .ge. tap(ntap,si)) then
        fg(2)=ntap+1
      endif

      do i=2,ntap
       if ((tap(i-1,si) .le. t2) .and. (t2 .lt. tap(i,si))) then
       fg(2)=i
       endif
      enddo


        un=0d0
        up=0d0
        vn=0d0
        vp=0d0
!****** CORNERS
       if( (fg(1)-1.eq.0) .and. (fg(2)-1.eq.0)) then
       vp=1d0
       up=1d0
       tk=1
       else if( (fg(1).eq.nrap+1) .and. (fg(2).eq.ntap+1)) then
       un=1d0
       vn=1d0
       tk=1
       else if( (fg(1)-1.eq.0) .and. (fg(2).eq.ntap+1)) then
       up=1d0
       vn=1d0
       tk=1
       else if( (fg(1).eq.nrap+1) .and. (fg(2)-1.eq.0)) then
       un=1d0
       vp=1d0
       tk=1
!****** SIDES
       else if( (fg(1)-1.eq.0) ) then
       up=(r1-rap(fg(1),si))**4d0
       vn=(t2-tap(fg(2)-1,si))**4d0
       vp=(t2-tap(fg(2),si))**4d0
       tk=2
       else if( (fg(1).eq.nrap+1) ) then
       un=(r1-rap((fg(1)-1),si))**4d0
       vn=(t2-tap(fg(2)-1,si))**4d0
       vp=(t2-tap(fg(2),si))**4d0
! Connection for r1=rap
          if (r1 .eq. rap((fg(1)-1),si)) then
           un=1d0
          endif
       tk=2
       else if( (fg(2)-1.eq.0) ) then
       un=(r1-rap((fg(1)-1),si))**4d0
       up=(r1-rap(fg(1),si))**4d0
       vp=(t2-tap(fg(2),si))**4d0
       tk=2
       else if( (fg(2).eq.ntap+1) ) then
       un=(r1-rap((fg(1)-1),si))**4d0
       up=(r1-rap(fg(1),si))**4d0
       vn=(t2-tap(fg(2)-1,si))**4d0
! Connection for t2=rap 
          if (t2 .eq. tap(fg(2)-1,si)) then
           vn=1d0
          endif
       tk=2
!****** MIDDLE
       else
       un=(r1-rap((fg(1)-1),si))**4d0
       up=(r1-rap(fg(1),si))**4d0
       vn=(t2-tap(fg(2)-1,si))**4d0
       vp=(t2-tap(fg(2),si))**4d0
       tk=3
       endif


       dnom = un*vn + un*vp + up*vn + up*vp 

        f(1)=un*vn/dnom
        f(2)=un*vp/dnom
        f(3)=up*vn/dnom
        f(4)=up*vp/dnom





        df=0d0
       if( (fg(1)-1.eq.0) .and. (fg(2)-1.eq.0)) then
       else if( (fg(1).eq.nrap+1) .and. (fg(2).eq.ntap+1)) then
       else if( (fg(1)-1.eq.0) .and. (fg(2).eq.ntap+1)) then
       else if( (fg(1).eq.nrap+1) .and. (fg(2)-1.eq.0)) then
! THEN nothing for corners
!****** CORNERS
! all derivatives equal zero
!****** SIDES
!      if( (fg(1)-1.eq.0) ) then
       else if( (fg(1)-1.eq.0) ) then
       dup=4d0*(r1-rap(fg(1),si))**3d0
       dvn=4d0*(t2-tap(fg(2)-1,si))**3d0
       dvp=4d0*(t2-tap(fg(2),si))**3d0

       else if( (fg(1).eq.nrap+1) ) then
       dun=4d0*(r1-rap((fg(1)-1),si))**3d0
       dvn=4d0*(t2-tap(fg(2)-1,si))**3d0
       dvp=4d0*(t2-tap(fg(2),si))**3d0

       else if( (fg(2)-1.eq.0) ) then
       dun=4d0*(r1-rap((fg(1)-1),si))**3d0
       dup=4d0*(r1-rap(fg(1),si))**3d0
       dvp=4d0*(t2-tap(fg(2),si))**3d0

       else if( (fg(2).eq.nrap+1) ) then
       dun=4d0*(r1-rap((fg(1)-1),si))**3d0
       dup=4d0*(r1-rap(fg(1),si))**3d0
       dvn=4d0*(t2-tap(fg(2)-1,si))**3d0

!****** MIDDLE
       else
       dun=4d0*(r1-rap((fg(1)-1),si))**3d0
       dup=4d0*(r1-rap(fg(1),si))**3d0
       dvn=4d0*(t2-tap(fg(2)-1,si))**3d0
       dvp=4d0*(t2-tap(fg(2),si))**3d0

       endif



       ddnom(1) = dun*vn + dun*vp + dup*vn + dup*vp
       ddnom(2) = un*dvn + un*dvp + up*dvn + up*dvp

       df(1,1)=((dun*vn)*dnom-ddnom(1)*un*vn)/(dnom**2d0)
       df(2,1)=((dun*vp)*dnom-ddnom(1)*un*vp)/(dnom**2d0)
       df(3,1)=((dup*vn)*dnom-ddnom(1)*up*vn)/(dnom**2d0)
       df(4,1)=((dup*vp)*dnom-ddnom(1)*up*vp)/(dnom**2d0)

       df(1,2)=((un*dvn)*dnom-ddnom(2)*un*vn)/(dnom**2d0)
       df(2,2)=((un*dvp)*dnom-ddnom(2)*un*vp)/(dnom**2d0)
       df(3,2)=((up*dvn)*dnom-ddnom(2)*up*vn)/(dnom**2d0)
       df(4,2)=((up*dvp)*dnom-ddnom(2)*up*vp)/(dnom**2d0)

      end

 

      subroutine KPev2gm2(igrad,X,Y,Z,v3,imol,dx,dy,dz) 
! gives the primary coordinate contribution

      implicit none
      integer imol





      double precision a12(4)
      double precision const(4)
      double precision ga1(4),ga2(4),ga3(4)
      double precision gb1(4),gb2(4),gb3(4)
      double precision ka1(4),ka2(4)
      double precision c(0:10)
      double precision re,de,rmin,pi,rx
      double precision fy,dfdr,u,dydr,lj,dljdr
      double precision rr,tt,dv2dr,dv2dt,v3
      integer ri,i,igrad
      double precision X(7),Y(7),Z(7)
      double precision dx(7),dy(7),dz(7)
      double precision ae,e1,e2,e3,f1,f2,f3
      double precision g1,g2,g3,g4,i2,i3,i4
      double precision h2,h3,h4,j2,j3,j4
      double precision n1,n2,n3,n4
      double precision bndlen
      double precision dihedr2
      double precision ak,tt0,gr1,gr2,gr3
      double precision fr1(5),fr2(5),fr3(5)
      double precision ft2(5),ft3(5)
      double precision dihedrIX,dihedrIY,dihedrIZ
      double precision dihedrJX,dihedrJY,dihedrJZ
      double precision dihedrKX,dihedrKY,dihedrKZ
      double precision dihedrLX,dihedrLY,dihedrLZ

        dv2dr=0d0
        dv2dt=0d0
        v3=0d0

        dx=0d0
        dy=0d0
        dz=0d0

        rr=bndlen(1,3,X,Y,Z)
        tt=dihedr2(1,3,4,2,X,Y,Z)


       if (imol .eq. 1) then ! GROUND STATE

      e1=6.420621d0
      e2=-1.85d0
      e3=1.01574847d0
      f1=0.6d0
      f2=-4.0d0
      f3=2.23d0

      v3=e1*(1.0d0-exp(e2*((rr-e3))))**2d0
     &       +f1*exp(f2*((rr-f3)**2d0))

      if (igrad .eq. 1) then
          dv2dr=-2.0d0*e1*e2*exp(e2*(rr-e3))
     &        *(1-exp(e2*(rr-e3)))
     &      + 2d0*f1*f2*(rr-f3)*exp(f2*(rr-f3)**2d0)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ak=0.945760042d0
      tt0=0.482599385d0
      gr1=1.472234651d0
      gr2=-5.656037749d0
      gr3=0.647196191d0

        v3=v3+
     & ak*(tt-tt0)**2d0*
     & gr1*exp(gr2*(rr-gr3)**2d0)


      if (igrad .eq. 1) then
        dv2dr=dv2dr+
     & ak*(tt-tt0)**2d0*
     & 2d0*gr1*gr2*(rr-gr3)*exp(gr2*(rr-gr3)**2d0)

        dv2dt=
     & ak*(tt-tt0)*2d0*
     & gr1*exp(gr2*(rr-gr3)**2d0)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       fr1(1)=3.067752289d0
       fr1(2)=1076.860912d0
       fr1(3)=0.792499074d0
       fr1(4)=-1.373902904d0

       fr2(1)=-0.845025076d0
       fr2(2)=-5.474618306d0
       fr2(3)=-5.214819223d0
       fr2(4)=-7.271506597d0

       fr3(1)=0.98801149d0
       fr3(2)=0.89839555d0
       fr3(3)=1.484837694d0
       fr3(4)=1.89501344d0

       ft2(1)=-3.340905161d0
       ft2(2)=-14.55654461d0
       ft2(3)=-12.70986496d0
       ft2(4)=-11.06393037d0

       ft3(1)=-0.917131468d0
       ft3(2)=-1.489784616d0
       ft3(3)=1.002868947d0
       ft3(4)=-0.334959847d0

        do i=1,4
        v3=v3+
     & fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))*
     & exp(ft2(i)*((tt-ft3(i))**2d0))
        enddo


      if (igrad .eq. 1) then
        do i=1,4
        dv2dr=dv2dr+
     & 2d0*fr2(i)*(rr-fr3(i))*
     & fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))*
     & exp(ft2(i)*((tt-ft3(i))**2d0))
        enddo

        do i=1,4
        dv2dt=dv2dt+
     & 2d0*ft2(i)*(tt-ft3(i))*
     & fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))*
     & exp(ft2(i)*((tt-ft3(i))**2d0))
        enddo
      endif



       else if (imol .eq. 2) then ! EXCITED STATE


      e1=4.509566d0
      n4=1.650445177d0
      e2=-1.460739955d0
      e3=0.582944229d0
      f1=0.234487181d0
      f2=-2.282322532d0
      f3=0.30831821d0
      n1=0.0d0
      n2=300d0
      n3=1.025d0

      v3=e1+n4*exp(e2*(rr-e3))
     &       +(n1+n2*(rr-n3)**2d0)
     &       *f1*exp(f2*((rr-f3)**2d0))

      if (igrad .eq. 1) then
          dv2dr=e2*n4*exp(e2*(rr-e3))
     &      + (n2*(rr-n3)*2d0)
     &      *f1*exp(f2*((rr-f3)**2d0))
     &       +(n1+n2*(rr-n3)**2d0)
     &       *2d0*f1*f2*(rr-f3)*exp(f2*(rr-f3)**2d0)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ak=1.233643227d0
      tt0=1.0d0
      gr1=0.978308804d0
      gr2=-0.981638325d0
      gr3=0.901675465d0

        v3=v3+
     & ak*(tt-tt0)**2d0*
     & gr1*exp(gr2*(rr-gr3)**2d0)


      if (igrad .eq. 1) then
        dv2dr=dv2dr+
     & ak*(tt-tt0)**2d0*
     & 2d0*gr1*gr2*(rr-gr3)*exp(gr2*(rr-gr3)**2d0)

        dv2dt=
     & ak*(tt-tt0)*2d0*
     & gr1*exp(gr2*(rr-gr3)**2d0)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       fr1(1)=5.402431748d0
       fr1(2)=-5.542582536d0
       fr1(3)=-3.460058164d0
       fr1(4)=-0.508881823d0
       fr1(5)=3.692631232d0

       fr2(1)=-5.951923443d0
       fr2(2)=-6.031843482d0
       fr2(3)=-9.434511421d0
       fr2(4)=-2.878092095d0
       fr2(5)=-9.180177476d0

       fr3(1)=1.867955914d0
       fr3(2)=1.870582014d0
       fr3(3)=1.787229743d0
       fr3(4)=2.302428489d0
       fr3(5)=1.746230571d0

       ft2(1)=-9.559209372d0
       ft2(2)=-9.696891789d0
       ft2(3)=-5.670799519d0
       ft2(4)=-8.901043282d0
       ft2(5)=-5.346117363d0

       ft3(1)=-0.590397385d0
       ft3(2)=-0.675159631d0
       ft3(3)=0.465695686d0
       ft3(4)=-0.5580416d0
       ft3(5)=0.427418036d0

        do i=1,5
        v3=v3+
     & fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))*
     & exp(ft2(i)*((tt-ft3(i))**2d0))
        enddo


      if (igrad .eq. 1) then
        do i=1,5
        dv2dr=dv2dr+
     & 2d0*fr2(i)*(rr-fr3(i))*
     & fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))*
     & exp(ft2(i)*((tt-ft3(i))**2d0))
        enddo

        do i=1,5
        dv2dt=dv2dt+
     & 2d0*ft2(i)*(tt-ft3(i))*
     & fr1(i)*exp(fr2(i)*((rr-fr3(i))**2d0))*
     & exp(ft2(i)*((tt-ft3(i))**2d0))
        enddo
      endif



       else if (imol .eq. 3) then ! COUPLING


       a12(1)=-0.16290134d0
       ga1(1)=-0.231270633d0
       ga2(1)=-42.62845833d0
       ga3(1)=-0.788781464d0
       gb1(1)=-0.157578396d0
       gb2(1)=-4.728738508d0
       gb3(1)=-0.492785959d0
       const(1)=0.050222306d0

       a12(2)=1.514877043d0
       ga1(2)=0.38885665d0
       ga2(2)=-52.43604225d0
       ga3(2)=-0.596932233d0
       gb1(2)=2.340167593d0
       gb2(2)=-0.158293136d0
       gb3(2)=0.03631615d0
       const(2)=-2.497398102d0

       a12(3)=0.96486818d0
       ga1(3)=-0.567995439d0
       ga2(3)=-14.16424012d0
       ga3(3)=-0.534406038d0
       gb1(3)=1.682033257d0
       gb2(3)=-0.103900297d0
       gb3(3)=0.059650329d0
       const(3)=-1.859104948d0

       a12(4)=0.730431712d0
       ga1(4)=-1.39153824d0
       ga2(4)=-0.081123542d0
       ga3(4)=-1.342253144d0
       gb1(4)=4.257325153d0
       gb2(4)=-1.963596235d0
       gb3(4)=-2.057979708d0
       const(4)=0.669000566d0


        ka1(1)=-2.904320563d0
        ka1(2)=-5.709140977d0
        ka1(3)=-5.612288257d0
        ka1(4)=-2.876200866d0


        ka2(1)=0.382538122d0
        ka2(2)=1.81713229d0
        ka2(3)=2.40922499d0
        ka2(4)=3.270326934d0


       do i=1,4
       v3=v3
     &    + (a12(i)*sqrt(1.1d0-(tt**2d0))
     &        +const(i)
     & + ga1(i)*exp(ga2(i)*((tt-ga3(i))**2d0))
     & + gb1(i)*exp(gb2(i)*((tt-gb3(i))**2d0)))
     & *(exp(ka1(i)*(rr-ka2(i))**2d0))
       enddo



       if (igrad .eq. 1) then

       do i=1,4
        dv2dt=dv2dt
     & + (-1.0d0*a12(i)*tt
     & *(1.0d0/sqrt(1.1d0-(tt**2d0)))
     & + ga1(i)*exp(ga2(i)*((tt-ga3(i))**2d0))
     & *ga2(i)*((tt-ga3(i))*2d0)
     & + gb1(i)*exp(gb2(i)*((tt-gb3(i))**2d0))
     & *gb2(i)*((tt-gb3(i))*2d0))
     & *(exp(ka1(i)*(rr-ka2(i))**2d0))
        enddo

        do i=1,4
        dv2dr=dv2dr
     &    + (a12(i)*sqrt(1.1d0-(tt**2d0))
     &        +const(i)
     & + ga1(i)*exp(ga2(i)*((tt-ga3(i))**2d0))
     & + gb1(i)*exp(gb2(i)*((tt-gb3(i))**2d0)))
     & *(exp(ka1(i)*(rr-ka2(i))**2d0))
     & *(ka1(i)*(rr-ka2(i))*2d0)
        enddo
       endif



       endif


      dx(1)=dv2dr*((X(1)-X(3))/rr)
      dy(1)=dv2dr*((Y(1)-Y(3))/rr)
      dz(1)=dv2dr*((Z(1)-Z(3))/rr)

      dx(3)=dv2dr*((X(3)-X(1))/rr)
      dy(3)=dv2dr*((Y(3)-Y(1))/rr)
      dz(3)=dv2dr*((Z(3)-Z(1))/rr)


      dx(1)=dx(1)+dv2dt*dihedrIX(1,3,4,2,X,Y,Z)
      dy(1)=dy(1)+dv2dt*dihedrIY(1,3,4,2,X,Y,Z)
      dz(1)=dz(1)+dv2dt*dihedrIZ(1,3,4,2,X,Y,Z)

      dx(2)=dv2dt*dihedrLX(1,3,4,2,X,Y,Z)
      dy(2)=dv2dt*dihedrLY(1,3,4,2,X,Y,Z)
      dz(2)=dv2dt*dihedrLZ(1,3,4,2,X,Y,Z)

      dx(3)=dx(3)+dv2dt*dihedrJX(1,3,4,2,X,Y,Z)
      dy(3)=dy(3)+dv2dt*dihedrJY(1,3,4,2,X,Y,Z)
      dz(3)=dz(3)+dv2dt*dihedrJZ(1,3,4,2,X,Y,Z)

      dx(4)=dv2dt*dihedrKX(1,3,4,2,X,Y,Z)
      dy(4)=dv2dt*dihedrKY(1,3,4,2,X,Y,Z)
      dz(4)=dv2dt*dihedrKZ(1,3,4,2,X,Y,Z)

      end



      double precision function dihedr2(i,j,k,l,X,Y,Z)
!************************************************************************
! Evaluate the dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
! modified version of dihedr1 to work through flat amine geoms
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision kkX,kkY,kkZ,yyX,yyY,yyZ
      double precision llX,llY,llZ,eeX,eeY,eeZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       kkX=X(i)-X(j)
       kkY=Y(i)-Y(j)
       kkZ=Z(i)-Z(j)

       yyX=X(k)-X(j)
       yyY=Y(k)-Y(j)
       yyZ=Z(k)-Z(j)

       llX=X(l)-X(k)
       llY=Y(l)-Y(k)
       llZ=Z(l)-Z(k)

       eeX=X(j)-X(k)
       eeY=Y(j)-Y(k)
       eeZ=Z(j)-Z(k)


       oo1(1)=kkY*yyZ-kkZ*yyY
       oo1(2)=kkZ*yyX-kkX*yyZ
       oo1(3)=kkX*yyY-kkY*yyX

       oo2(1)=llY*eeZ-llZ*eeY
       oo2(2)=llZ*eeX-llX*eeZ
       oo2(3)=llX*eeY-llY*eeX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       ra=sab
!       dihedr2=acos(ra)
       dihedr2=ra


!--------------------------------------------------!

       end


      double precision function dihedrIX(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)
!       dihedr2=ra




!------------------Xi---------------------------
!top (sab deriv for Xi)

       td=-1d0*kjZ*lkZ*jkX
     &          +kjZ*lkX*jkZ
     &          +kjY*lkX*jkY
     &          -kjY*lkY*jkX

!bottom (p1p2 deriv for Xi)
! deriv of sso1 and bd for Xi

       ld=2d0*(X(i)-X(j))*(kjZ**2d0+kjY**2d0)
     &     -2*(kjZ*ijZ*kjX
     &        +kjY*ijY*kjX)

        bd=(sso2*ld)/(2d0*p1p2)

!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrIX=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrIX=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Xi---------------------------


       end




      double precision function dihedrIY(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)


!------------------Yi---------------------------

!top (sab deriv for Yi)

       td=-1d0*kjZ*lkZ*jkY
     &          +kjZ*lkY*jkZ
     &          +kjX*lkY*jkX
     &          -kjX*lkX*jkY

! deriv of sso1 and bd for Yi

       ld=2d0*(Y(i)-Y(j))*(kjZ**2+kjX**2)
     &     -2*(kjZ*ijZ*kjY
     &        +kjX*ijX*kjY)

        bd=(sso2*ld)/(2d0*p1p2)



!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrIY=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrIY=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Yi---------------------------

       end





      double precision function dihedrIZ(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)


!------------------Zi---------------------------

!top (sab deriv for Zi)

       td=-1d0*kjY*lkY*jkZ
     &          +kjY*lkZ*jkY
     &          +kjX*lkZ*jkX
     &          -kjX*lkX*jkZ

! deriv of sso1 and bd for Zi

       ld=2d0*(Z(i)-Z(j))*(kjX**2+kjY**2)
     &     -2*(kjX*ijX*kjZ
     &        +kjY*ijY*kjZ)

        bd=(sso2*ld)/(2d0*p1p2)


!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrIZ=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrIZ=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Zi---------------------------


       end




      double precision function dihedrLX(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)


!------------------Xl---------------------------


!top (sab deriv for Xl)

       td=-1d0*ijZ*kjX*jkZ
     &          +ijX*kjZ*jkZ
     &          +ijX*kjY*jkY
     &          -ijY*kjX*jkY


! deriv of sso2 for Xl

       rd=2d0*(X(l)-X(k))*(jkY**2+jkZ**2)
     &     -2*(jkY*lkY*jkX
     &        +jkZ*lkZ*jkX)

        bd=(sso1*rd)/(2d0*p1p2)


!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrLX=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrLX=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Xl---------------------------

       end
 

      double precision function dihedrLY(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)


!------------------Yl---------------------------

!top (sab deriv for Yl)

       td=-1d0*ijZ*kjY*jkZ
     &          +ijY*kjZ*jkZ
     &          +ijY*kjX*jkX
     &          -ijX*kjY*jkX


! deriv of sso2 for Yl

       rd=2d0*(Y(l)-Y(k))*(jkX**2+jkZ**2)
     &     -2*(jkX*lkX*jkY
     &        +jkZ*lkZ*jkY)

        bd=(sso1*rd)/(2d0*p1p2)



!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrLY=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrLY=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Yl---------------------------

       end



      double precision function dihedrLZ(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)


!------------------Zl---------------------------

!top (sab deriv for Zl)

       td=-1d0*ijY*kjZ*jkY
     &          +ijZ*kjY*jkY
     &          +ijZ*kjX*jkX
     &          -ijX*kjZ*jkX


! deriv of sso2 and bd for Zl

       rd=2d0*(Z(l)-Z(k))*(jkY**2+jkX**2)
     &     -2*(jkY*lkY*jkZ
     &        +jkX*lkX*jkZ)


        bd=(sso1*rd)/(2d0*p1p2)


!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrLZ=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrLZ=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Zl---------------------------
       end


      double precision function dihedrJX(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)


!------------------Xj---------------------------

!top (sab deriv for Xj)


      td=(-2d0*X(j)+2d0*X(k))*ijZ*lkZ
     &     +ijZ*lkX*jkZ
     &     -(-2d0*X(j)+X(k)+X(i))*kjZ*lkZ
     &     -kjZ*lkX*jkZ
     &     -kjY*lkX*jkY
     &     -(-2d0*X(j)+X(k)+X(i))*kjY*lkY
     &     +ijY*lkX*jkY
     &     +(-2d0*X(j)+2d0*X(k))*ijY*lkY

! deriv of sso1 for Xj

       ld=2d0*(X(j)-X(k))*(ijZ**2+ijY**2)
     &     +2d0*(X(j)-X(i))*(kjZ**2+kjY**2)
     &     -2*(2d0*X(j)-X(k)-X(i))*(ijZ*kjZ
     &                             +ijY*kjY)


! deriv of sso2 for Xj

       rd=2d0*(X(j)-X(k))*(lkZ**2+lkY**2)
     &     -2*(lkZ*lkX*jkZ
     &        +lkY*lkX*jkY)


        bd=((ld*sso2)+(sso1*rd))/(2d0*p1p2)



!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrJX=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrJX=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Xj---------------------------
       end



      double precision function dihedrJY(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)


!------------------Yj---------------------------

!top (sab deriv for Yj)


      td=-kjZ*lkY*jkZ
     &     -(-2d0*Y(j)+Y(k)+Y(i))*kjZ*lkZ
     &     +ijZ*lkY*jkZ
     &     +(-2d0*Y(j)+2d0*Y(k))*ijZ*lkZ
     &     +(-2d0*Y(j)+2d0*Y(k))*ijX*lkX
     &     +ijX*lkY*jkX
     &     -(-2d0*Y(j)+Y(k)+Y(i))*kjX*lkX
     &     -kjX*lkY*jkX


! deriv of sso1 for Yj

       ld=2d0*(Y(j)-Y(k))*(ijZ**2+ijX**2)
     &     +2d0*(Y(j)-Y(i))*(kjZ**2+kjX**2)
     &     -2*(2d0*Y(j)-Y(k)-Y(i))*(ijZ*kjZ
     &                             +ijX*kjX)



! deriv of sso2 for Yj

       rd=2d0*(Y(j)-Y(k))*(lkZ**2+lkX**2)
     &     -2*(lkZ*lkY*jkZ
     &        +lkX*lkY*jkX)


        bd=((ld*sso2)+(sso1*rd))/(2d0*p1p2)



!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrJY=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrJY=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Yj---------------------------
       end



      double precision function dihedrJZ(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)


!------------------Zj---------------------------

!top (sab deriv for Zj)

      td=(-2d0*Z(j)+2d0*Z(k))*ijY*lkY
     &     +ijY*lkZ*jkY
     &     -(-2d0*Z(j)+Z(k)+Z(i))*kjY*lkY
     &     -kjY*lkZ*jkY
     &     -kjX*lkZ*jkX
     &     -(-2d0*Z(j)+Z(k)+Z(i))*kjX*lkX
     &     +ijX*lkZ*jkX
     &     +(-2d0*Z(j)+2d0*Z(k))*ijX*lkX

! deriv of sso1 for Zj

       ld=2d0*(Z(j)-Z(k))*(ijY**2+ijX**2)
     &     +2d0*(Z(j)-Z(i))*(kjY**2+kjX**2)
     &     -2*(2d0*Z(j)-Z(k)-Z(i))*(ijY*kjY
     &                             +ijX*kjX)


! deriv of sso2 for Zj

       rd=2d0*(Z(j)-Z(k))*(lkY**2+lkX**2)
     &     -2*(lkY*lkZ*jkY
     &        +lkX*lkZ*jkX)


        bd=((ld*sso2)+(sso1*rd))/(2d0*p1p2)



!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrJZ=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrJZ=id!/(1d0-ra**2d0)**(0.5d0)

!------------------Zj---------------------------
       end



      double precision function dihedrKX(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)









!------------------Xk---------------------------

!top (sab deriv for Xk)


      td=(-2d0*X(k)+2d0*X(j))*ijZ*lkZ
     &     -(-2d0*X(k)+X(j)+X(l))*ijZ*jkZ
     &     +ijX*kjZ*lkZ
     &     -ijX*kjZ*jkZ
     &     -ijX*kjY*jkY
     &     +ijX*kjY*lkY
     &     -(-2d0*X(k)+X(j)+X(l))*ijY*jkY
     &     +(-2d0*X(k)+2d0*X(j))*ijY*lkY


! deriv of sso1 for Xk

       ld=2d0*(X(k)-X(j))*(ijZ**2+ijY**2)
     &     -2*(ijX*ijZ*kjZ
     &        +ijX*ijY*kjY)


! deriv of sso2 for Xk


       rd=2d0*(X(k)-X(j))*(lkZ**2+lkY**2)
     &   +2d0*(X(k)-X(l))*(jkZ**2+jkY**2)
     &   -2*(2d0*X(k)-X(j)-X(l))*(lkZ*jkZ
     &                           +lkY*jkY)


        bd=((ld*sso2)+(sso1*rd))/(2d0*p1p2)


!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrKX=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrKX=id!(-id)/(1d0-ra**2d0)**(0.5d0)

!------------------Xk---------------------------



       end



      double precision function dihedrKY(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!        i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)





!------------------Yk---------------------------

!top (sab deriv for Yk)


      td=(-2d0*Y(k)+2d0*Y(j))*ijZ*lkZ
     &     -(-2d0*Y(k)+Y(j)+Y(l))*ijZ*jkZ
     &     +ijY*kjZ*lkZ
     &     -ijY*kjZ*jkZ
     &     -ijY*kjX*jkX
     &     +ijY*kjX*lkX
     &     -(-2d0*Y(k)+Y(j)+Y(l))*ijX*jkX
     &     +(-2d0*Y(k)+2d0*Y(j))*ijX*lkX


! deriv of sso1 for Yk

       ld=2d0*(Y(k)-Y(j))*(ijZ**2+ijX**2)
     &     -2*(ijY*ijZ*kjZ
     &        +ijY*ijX*kjX)


! deriv of sso2 for Yk


       rd=2d0*(Y(k)-Y(j))*(lkZ**2+lkX**2)
     &   +2d0*(Y(k)-Y(l))*(jkZ**2+jkX**2)
     &   -2*(2d0*Y(k)-Y(j)-Y(l))*(lkZ*jkZ
     &                           +lkX*jkX)


        bd=((ld*sso2)+(sso1*rd))/(2d0*p1p2)


!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrKY=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrKY=id!(-id)/(1d0-ra**2d0)**(0.5d0)

!------------------Yk---------------------------


       end


      double precision function dihedrKZ(i,j,k,l,X,Y,Z)
!************************************************************************
! Converts derivative from internal to XYZ for
! dihedral angle between atoms i, j, k, and l
! Input: i,j,k,l (indices of 4 atoms)
!         i--j--k--l
!************************************************************************

      implicit double precision(a-h,o-z)

      double precision tiny,pi
      parameter(tiny=1.0d-13,pi=3.1415926535898d0)

      double precision X(*),Y(*),Z(*)
      double precision eij(3),ejk(3),ekl(3),u(3),v(3),w(3)
      double precision rij,rjk,rkl,aijk,ajkl,wejk
      double precision ijX,ijY,ijZ,kjX,kjY,kjZ
      double precision lkX,lkY,lkZ,jkX,jkY,jkZ
      double precision oo1(3),oo2(3),so1(3),so2(3)
      double precision sso1,sso2,rso1,rso2
      double precision ab(3)
      double precision sab,p1p2,ra
      double precision ld,rd,td,bd,id,fd

!   eij
! i-----j
!   aijk \ ejk
!         \ 
!          \ ajkl
!           k----l
!             kl   


       ijX=X(i)-X(j)
       ijY=Y(i)-Y(j)
       ijZ=Z(i)-Z(j)

       kjX=X(k)-X(j)
       kjY=Y(k)-Y(j)
       kjZ=Z(k)-Z(j)

       lkX=X(l)-X(k)
       lkY=Y(l)-Y(k)
       lkZ=Z(l)-Z(k)

       jkX=X(j)-X(k)
       jkY=Y(j)-Y(k)
       jkZ=Z(j)-Z(k)


       oo1(1)=ijY*kjZ-ijZ*kjY
       oo1(2)=ijZ*kjX-ijX*kjZ
       oo1(3)=ijX*kjY-ijY*kjX

       oo2(1)=lkY*jkZ-lkZ*jkY
       oo2(2)=lkZ*jkX-lkX*jkZ
       oo2(3)=lkX*jkY-lkY*jkX


       so1(1)=oo1(1)**2d0
       so1(2)=oo1(2)**2d0
       so1(3)=oo1(3)**2d0
       sso1=so1(1)+so1(2)+so1(3)
       rso1=(sso1)**(0.5d0)

       so2(1)=oo2(1)**2d0
       so2(2)=oo2(2)**2d0
       so2(3)=oo2(3)**2d0
       sso2=so2(1)+so2(2)+so2(3)
       rso2=(sso2)**(0.5d0)

       ab(1)=oo1(1)*oo2(1)
       ab(2)=oo1(2)*oo2(2)
       ab(3)=oo1(3)*oo2(3)
       sab=ab(1)+ab(2)+ab(3)

       p1p2=rso1*rso2

       ra=sab/p1p2
!       dihedr2=acos(ra)

!------------------Zk---------------------------

!top (sab deriv for Zk)


      td=(-2d0*Z(k)+2d0*Z(j))*ijY*lkY
     &     -(-2d0*Z(k)+Z(j)+Z(l))*ijY*jkY
     &     +ijZ*kjY*lkY
     &     -ijZ*kjY*jkY
     &     -ijZ*kjX*jkX
     &     +ijZ*kjX*lkX
     &     -(-2d0*Z(k)+Z(j)+Z(l))*ijX*jkX
     &     +(-2d0*Z(k)+2d0*Z(j))*ijX*lkX


! deriv of sso1 for Zk

       ld=2d0*(Z(k)-Z(j))*(ijY**2+ijX**2)
     &     -2*(ijZ*ijY*kjY
     &        +ijZ*ijX*kjX)


! deriv of sso2 for Zk


       rd=2d0*(Z(k)-Z(j))*(lkY**2+lkX**2)
     &   +2d0*(Z(k)-Z(l))*(jkY**2+jkX**2)
     &   -2*(2d0*Z(k)-Z(j)-Z(l))*(lkY*jkY
     &                           +lkX*jkX)


        bd=((ld*sso2)+(sso1*rd))/(2d0*p1p2)


!ra=sab/p1p2
       id=((td*p1p2)-(bd*sab))/(p1p2**2d0)
!       dihedrKZ=(-id)/(1d0-ra**2d0)**(0.5d0)
       dihedrKZ=id!(-id)/(1d0-ra**2d0)**(0.5d0)

!------------------Zk---------------------------


       end

