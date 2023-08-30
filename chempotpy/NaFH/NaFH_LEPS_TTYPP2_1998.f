      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=2
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: nt, r(1,3), r2(3), v(1)
      double precision :: u11(1), u12(1), u22(1)
      double precision :: de_u11(3,1), de_u12(3,1), de_u22(3,1)
      double precision :: u(nstates,nstates), dudr(nstates,nstates,3)
      double precision :: t(nstates,nstates)
      double precision :: tmpmat(nstates,nstates)
      double precision :: hr(nstates,nstates,3), gr(nstates,3)
      double precision :: tx(9), drdx(3,9)
      double precision :: hx(nstates,nstates,9), gx(nstates,9)
      logical, save :: first_time_data=.true.
      integer :: i, j, k, l

      !initialize 
      u=0.d0
      p=0.d0
      g=0.d0
      d=0.d0

      nt=1
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        tx(j)=x(iatom,idir)
      enddo
      enddo
      ! input cartesian is ClHH
      r(1,1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      r(1,2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      r(1,3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211

      NDER=igrad

      if(first_time_data) then
      call prepot
      first_time_data=.false.
      endif

      call pot(r,u11,de_u11,1,1)
      call pot(r,u12,de_u12,1,2)
      call pot(r,u22,de_u22,1,3)
      u(1,1)=u11(1)*27.211386
      u(1,2)=u12(1)*27.211386
      u(2,1)=u12(1)*27.211386
      u(2,2)=u22(1)*27.211386
      dudr(1,1,:)=de_u11(:,1)*51.422067
      dudr(1,2,:)=de_u12(:,1)*51.422067
      dudr(2,1,:)=de_u12(:,1)*51.422067
      dudr(2,2,:)=de_u22(:,1)*51.422067
      call diagonalize(nstates,u,p,t)

      do i=1,3
        tmpmat(:,:)=dudr(:,:,i)
        tmpmat=matmul(transpose(t), matmul(tmpmat,t))
        do j=1,nstates
          do k=j+1,nstates
            if ((p(k)-p(j)) > 1.d-8) then
              hr(j,k,i)=tmpmat(j,k)/(p(k)-p(j))
            else
              hr(j,k,i)=tmpmat(j,k)/1.d-8
            endif
            hr(k,j,i)=-hr(j,k,i)
          enddo 
        enddo
        do j=1,nstates
          gr(j,i)=tmpmat(j,j)
        enddo
      enddo
      
      r2=r(1,:)*0.529177211
      call evdrdx(tx, r2, drdx)
      hx=0.d0
      gx=0.d0
      do i=1,9
      do j=1,3
        hx(:,:,i)=hx(:,:,i)+hr(:,:,j)*drdx(j,i)
        gx(:,i)=gx(:,i)+gr(:,j)*drdx(j,i)
      enddo
      enddo
  
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        g(:,iatom,idir)=gx(:,j)
        d(:,:,iatom,idir)=hx(:,:,j)
      enddo
      enddo


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
      subroutine prepot
C
C
C   System:          NaFH
C   Common name:     NaFH surface fit A
C   Number of derivatives: 1
C   Number of electronic surfaces: 2
C   Interface: 3-2V

C   References:  M.S. Topaler, D.G. Truhlar, X.Y. Chang, P. Piecuch, 
C                and J.C. Polanyi, J. Chem. Phys. 108, 5349 (1998).

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
C      
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank,g_nprocs

        call prepot11
        call prepot12
        call prepot22

      return

      end
c
      subroutine pot(r,e,De,nt,nsurf)
C
C         R     is an array with dimensions R(NT,3)
C               this array contains the geometry(ies) to be computed  (input)
C               Geometries are specified in terms of interatomic distances so
C               that
C               R(i,1) should be the distance between atoms 1 and 2 (A-B),
C               R(i,2) should be the distance between atoms 2 and 3 (B-C),
C               R(i,3) should be the distance between atoms 1 and 3 (A-C),
C               where i is an index in the range 1 <= i <= NT.
C         E     is an array with dimensions E(NT)
C               this array contains the diabatic potential energy(ies)  (output)
C         DE    is an array of the derivatives of the diabatic potential
C               energy(ies) with dimensions DE(3,NT) where the second
C               index labels the derivative with respect to the specific
C               component of R (interatomic distance)  (output)
C         NT    is the number of geometries to be computed  (input)
C         NSURF indicates the surface for which the energy should be
C               computed  (input)
C               The parameter NSURF may have the values 1, 2, or 3.
C               These values select a particular diabatic surface:
C               NSURF = 1  diabatic surface which is the lower one
C                          in the reactant channel
C               NSURF = 2  diabatic coupling surface
C               NSURF = 3  diabatic surface which is the upper one
C                          in the reactant channel
C
        implicit double precision (a-h,o-z)

        dimension r(nt,3),e(nt)
        dimension De(3,nt)

        if (nsurf .eq. 1) then
          call pot11(r,e,De,nt)
        else if (nsurf .eq. 2) then
          call pot12(r,e,De,nt)
        else
          call pot22(r,e,De,nt)
        end if

      return

      end
c
c
c  U11 surface - lower diabatic surface
c
      subroutine prepot11
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt)
      dimension De(3,nt)
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank,g_nprocs

      save /onetwo/,/uone/,/angular1/,rac2,af1,af2

      common/onetwo/CNF(6),CNF1(5),CHF(6),bc(0:5),
     & ads1,ars1,abs1,amh,amf,amna,amu,una2p,asmall,cevau,
     & DeNF,reNF,DeHF,reHF
      common/uone/y(4,3,4),c2a(4),c2b(4),c2c(4),dip1(4),
     &            adip1(4),rdip1(4),acut1(4),rcut1(4),
     &            dipn(4),al2(4),aln0(4),
     &            dipnn(4),al2n(4),aln0n1(4),aln0n2(4),rn0n(4)
      common/angular1/ac(4)
C-----------------------------------------------------------------------

      rac2=1.d0/dsqrt(2.d0)

      af1 = amh/(amf+amh)
      af2 = amf/(amf+amh)

      return
c
      entry pot11(r,e,De,nt)

      do 10 i=1,nt

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)
      rnag2 = af1*r1s+af2*r3s-af1*af2*r2s
      Drnag2D1 = 2.d0*af1*r(i,1)
      Drnag2D2 = -2.d0*af1*af2*r(i,2)
      Drnag2D3 = 2.d0*af2*r(i,3)

C     Cosine of Na-F-H bond angle:
      hlp = 1.d0/(r(i,2)*r(i,3))
      csb = 0.5d0*(r2s+r3s-r1s)*hlp
      DcsbD1 = -2.d0*r(i,1)*hlp
      DcsbD2 = 2.d0*(1.d0/r(i,3)-csb/r(i,2))
      DcsbD3 = 2.d0*(1.d0/r(i,2)-csb/r(i,3))
      csb = 2.d0*(1.d0+csb)

      if (csb.ge.4.d0) csb = 3.99999d0
      l1 = idint(csb)+1
      l2 = l1+1
      if (l1.eq.4) l2=l1
      if (csb.le.bc(l1)) csb = bc(l1)+1.d-13
      if (csb.ge.bc(l1+1)) csb = bc(l1+1)-1.d-13

      axs1=dexp(-abs1*(r(i,1)-ars1))
      s1=ads1*axs1*(axs1-2.d0)
      Ds1D1 = 2.d0*ads1*abs1*axs1*(1.d0-axs1)
C-----------------------------------------------------------------------

      y2 = r(i,2)-reHF
      hlp1 = cHF(1)*dexp(-cHF(2)*y2)
      hlp2 = dexp(-cHF(6)*y2)
      hlp22 = (-1.d0-cHF(1)+cHF(3)*y2+cHF(4)*y2**2+cHF(5)*y2**3)*hlp2
      s2 = DeHF*(hlp1+hlp22)
      Ds2D2 = DeHF*(-cHF(2)*hlp1 - cHF(6)*hlp22
     &   + (cHF(3)+2.d0*cHF(4)*y2+3.d0*cHF(5)*y2**2)*hlp2)
C-----------------------------------------------------------------------

      y3 = r(i,3)-reNF
      hlp1 = dexp(-cNF1(3)*y3)
      hlp11 = (cNF1(1)+cNF1(2)*y3)*hlp1
      hlp2 = cNF1(4)*dexp(-cNF1(5)*y3)
      s3 = hlp11+hlp2
      Ds3D3 = -cNF1(3)*hlp11+cNF1(2)*hlp1-cNF1(5)*hlp2
C-----------------------------------------------------------------------

      e(i) = 0.d0
      De(1,i) = 0.d0
      De(2,i) = 0.d0
      De(3,i) = 0.d0
      anorm = 0.d0
      DanormD1 = 0.d0
      DanormD2 = 0.d0
      DanormD3 = 0.d0
      do 20 l = l1, l2

      hlp1 = y(1,1,l)*dexp(-y(2,1,l)*r(i,1))
      hlp2 = y(3,1,l)*dexp(-y(4,1,l)*r(i,1))
      t1=hlp1+hlp2+s1
      Dt1D1 = -y(2,1,l)*hlp1-y(4,1,l)*hlp2+Ds1D1

      hlp1 = y(1,2,l)*dexp(-y(2,2,l)*r(i,2))
      hlp2 = y(3,2,l)*dexp(-y(4,2,l)*r(i,2))
      t2=hlp1+hlp2+s2
      Dt2D2=-y(2,2,l)*hlp1-y(4,2,l)*hlp2+Ds2D2
    
      hlp1 = y(1,3,l)*dexp(-y(2,3,l)*r(i,3))
      hlp2 = y(3,3,l)*dexp(-y(4,3,l)*r(i,3))
      t3 = hlp1+hlp2+s3
      Dt3D3=-y(2,3,l)*hlp1-y(4,3,l)*hlp2+Ds3D3

      coul1=0.5d0*(s1+t1)
      Dcl1D1=0.5d0*(Ds1D1+Dt1D1)
      coul2=0.5d0*(s2+t2)
      Dcl2D2=0.5d0*(Ds2D2+Dt2D2)
      coul3=0.5d0*(s3+t3)
      Dcl3D3=0.5d0*(Ds3D3+Dt3D3)
      exch1=0.5d0*(s1-t1)
      Dex1D1=0.5d0*(Ds1D1-Dt1D1)
      exch2=0.5d0*(s2-t2)
      Dex2D2=0.5d0*(Ds2D2-Dt2D2)
      exch3=0.5d0*(s3-t3)
      Dex3D3=0.5d0*(Ds3D3-Dt3D3)

      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      DwD1 = 2.d0*(2.d0*exch1-exch2-exch3)*Dex1D1
      DwD2 = 2.d0*(2.d0*exch2-exch1-exch3)*Dex2D2
      DwD3 = 2.d0*(2.d0*exch3-exch1-exch2)*Dex3D3

      cplg2=c2a(l)*dexp(-c2b(l)*w-c2c(l)*(r(i,1)+r(i,2)+r(i,3)))
      DcD1 = (-c2b(l)*DwD1 - c2c(l))*cplg2
      DcD2 = (-c2b(l)*DwD2 - c2c(l))*cplg2
      DcD3 = (-c2b(l)*DwD3 - c2c(l))*cplg2

      rootw = rac2*dsqrt(w+cplg2*cplg2)
      surf=(coul1+coul2+coul3+DeHF-rootw)
c
c check for zero by zero division (a case when all three atoms are asymptotically
c separated)
      if( (DwD1+2.d0*cplg2*DcD1) .eq. 0.0d0 )then
         DsurfD1=Dcl1D1
      else
         DsurfD1=(Dcl1D1-0.25d0*(DwD1+2.d0*cplg2*DcD1)/rootw)
      endif
      if( (DwD2+2.d0*cplg2*DcD2) .eq. 0.0d0 )then
         DsurfD2=Dcl2D2
      else
         DsurfD2=(Dcl2D2-0.25d0*(DwD2+2.d0*cplg2*DcD2)/rootw)
      endif
      if( (DwD3+2.d0*cplg2*DcD3) .eq. 0.0d0 )then
         DsurfD3=Dcl3D3
      else
         DsurfD3=(Dcl3D3-0.25d0*(DwD3+2.d0*cplg2*DcD3)/rootw)
      endif

      hlp1 = adip1(l)*(rnag2-rdip1(l))
      hlp2 = acut1(l)*(r(i,2)-rcut1(l))
      hlp11 = 0.5d0/(dcosh(hlp1)**2)
      hlp22 = 0.5d0/(dcosh(hlp2)**2)
      hlp1 = 0.5d0*(dtanh(hlp1)-1.d0)
      hlp2 = 0.5d0*(1.d0-dtanh(hlp2))
      surf = surf+dip1(l)*hlp1*hlp2
      DsurfD1=DsurfD1+dip1(l)*adip1(l)*Drnag2D1*hlp11*hlp2
      DsurfD2=DsurfD2+dip1(l)*adip1(l)*Drnag2D2*hlp11*hlp2
     &               -dip1(l)*hlp1*acut1(l)*hlp22
      DsurfD3=DsurfD3+dip1(l)*adip1(l)*Drnag2D3*hlp11*hlp2

      hlp1 = al2(l)*(r(i,2)-reHF)
      hlp2 = aln0(l)*rnag2
      hlp11 = dtanh(hlp1)
      hlp22 = dtanh(hlp2)
      hlp1 = 1.d0/dcosh(hlp1)
      hlp2 = 1.d0/dcosh(hlp2)
      surf = surf - dipn(l)*hlp1*hlp2
      DsurfD1=DsurfD1 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*Drnag2D1
      DsurfD2=DsurfD2 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*Drnag2D2
     &                + dipn(l)*hlp2*hlp1*hlp11*al2(l)
      DsurfD3=DsurfD3 + dipn(l)*hlp1*hlp2*hlp22*aln0(l)*Drnag2D3

      hlp1 = al2n(l)*(r(i,2)-reHF)
      hlp2 = aln0n1(l)*(rnag2-rn0n(l)**2)
c try to prevent an overflow
      if( hlp2 .lt. 709.0d0 ) then
         hlp2 = dexp(hlp2)
      else
         hlp2 = 1.0D308
      endif
      hlp3 = dexp(-aln0n2(l)*(rnag2-rn0n(l)**2))
      hlp23 = 1.d0/(hlp2+hlp3)**2
      hlp11 = dtanh(hlp1)
      hlp1 = 1.d0/dcosh(hlp1)
      surf = surf - 2.d0*dipnn(l)*hlp1/(hlp2+hlp3)
      DsurfD1=DsurfD1 + 2.d0*dipnn(l)*hlp1
     &       *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*Drnag2D1
      DsurfD2=DsurfD2 + 2.d0*dipnn(l)*hlp1
     &       *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*Drnag2D2
     &       + 2.d0*dipnn(l)*hlp1*hlp11*al2n(l)/(hlp2+hlp3)
      DsurfD3=DsurfD3 + 2.d0*dipnn(l)*hlp1
     &       *(aln0n1(l)*hlp2-aln0n2(l)*hlp3)*hlp23*Drnag2D3

      hlp = 1.d0/(
     &              (1.d0-(csb-bc(l))/(bc(l-1)-bc(l)))
     &            * (1.d0-(csb-bc(l))/(bc(l+1)-bc(l)))
     &                                                 )
      hlp1 = (csb-bc(l))**2*hlp
      cog = dexp(-ac(l)*hlp1)
      Dcog = -cog*ac(l)*(2.d0*(csb-bc(l))*hlp
     & +hlp1*(1.d0/(bc(l-1)-csb)+1.d0/(bc(l+1)-csb)))
C-----------------------------------------------------------------------
      DcogD1 = Dcog*DcsbD1
      DcogD2 = Dcog*DcsbD2
      DcogD3 = Dcog*DcsbD3

      e(i) = e(i) + surf*cog
      De(1,i) = De(1,i) + DsurfD1*cog+surf*DcogD1
      De(2,i) = De(2,i) + DsurfD2*cog+surf*DcogD2
      De(3,i) = De(3,i) + DsurfD3*cog+surf*DcogD3
      anorm = anorm + cog
      DanormD1 = DanormD1 + DcogD1
      DanormD2 = DanormD2 + DcogD2
      DanormD3 = DanormD3 + DcogD3
   20 continue
c     if (csb.gt.3.99d0) write(6,77)l1,l2,csb,r(i,3),
c    &                         r(i,2),e(i),anorm
      hlp = 1.d0/anorm
      e(i) = e(i)*hlp
      De(1,i) = (De(1,i) - e(i)*DanormD1)*hlp
      De(2,i) = (De(2,i) - e(i)*DanormD2)*hlp
      De(3,i) = (De(3,i) - e(i)*DanormD3)*hlp
c  77 format(2(2x,i3),5(2x,e16.8))

      e(i) = e(i)*cevau
      De(1,i) = De(1,i)*cevau
      De(2,i) = De(2,i)*cevau
      De(3,i) = De(3,i)*cevau
10    continue

C-----------------------------------------------------------------------
  100 format(/26('='),' NaFH - U11 diabat ',26('='),/
     &2x,'mH/(mF+mH) = ',f13.5,2x,'mF/(mF+mH) = ',f13.5)
  102 format(2x,'cNaF1 = ',9x,3e14.6/19x,2e14.6)
  103 format(2x,'DeHF = ',10x,e14.6/2x,'reHF = ',10x,e14.6
     & /2x,'cHF = ',11x,3e14.6/19x,3e14.6)
  104 format(2x,'ads1,ars1,abs1 = ',3e14.6)
  105 format(2x,'una2p,asmall = ',2x,2e14.6)
  106 format(2x,'y(1,1,i) = ',2x,4e13.5)
  107 format(2x,'y(2,1,i) = ',2x,4e13.5)
  108 format(2x,'y(3,1,i) = ',2x,4e13.5)
  109 format(2x,'y(4,1,i) = ',2x,4e13.5)
  110 format(2x,'y(1,2,i) = ',2x,4e13.5)
  111 format(2x,'y(2,2,i) = ',2x,4e13.5)
  112 format(2x,'y(3,2,i) = ',2x,4e13.5)
  113 format(2x,'y(4,2,i) = ',2x,4e13.5)
  114 format(2x,'y(1,3,i) = ',2x,4e13.5)
  115 format(2x,'y(2,3,i) = ',2x,4e13.5)
  116 format(2x,'y(3,3,i) = ',2x,4e13.5)
  117 format(2x,'y(4,3,i) = ',2x,4e13.5)
  118 format(2x,'c2a = ',7x,4e13.5)
  119 format(2x,'c2b = ',7x,4e13.5)
  120 format(2x,'c2c = ',7x,4e13.5)
  121 format(2x,'dip1 = ',6x,4e13.5)
  122 format(2x,'adip1 = ',5x,4e13.5)
  123 format(2x,'rdip1 = ',5x,4e13.5)
  124 format(2x,'acut1 = ',5x,4e13.5)
  125 format(2x,'rcut1 = ',5x,4e13.5)
  126 format(2x,'dipn = ',6x,4e13.5)
  127 format(2x,'al2 = ',7x,4e13.5)
  128 format(2x,'aln0 = ',6x,4e13.5)
  129 format(2x,'dipnn = ',5x,4e13.5)
  130 format(2x,'al2n = ',6x,4e13.5)
  131 format(2x,'aln0n1 = ',4x,4e13.5)
  132 format(2x,'aln0n2 = ',4x,4e13.5)
  133 format(2x,'rn0n = ',6x,4e13.5
     &    /71('=')/)

      return

      end
c
c  U12 surface - coupling surface
c
      subroutine prepot12
      implicit double precision (a-h,o-z)
      dimension r(nt,3),e(nt)
      dimension De(3,nt)
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank,g_nprocs

      save /onetwo/,/ucoupl/,/ucoupl1/,/u12nah/,af1,af2

      common/onetwo/CNF(6),CNF1(5),CHF(6),bc(0:5),
     & ads1,ars1,abs1,amh,amf,amna,amu,una2p,asmall,cevau,
     & DeNF,reNF,DeHF,reHF
      common/ucoupl/ahf(3),anaf(3),anah(3),r21(3),r23(3),r13(3)
      common/ucoupl1/ccNF(4),ccNHa(4),ccNHi(4),ccNHi1(4)
      common/u12nah/xpNH,xpNH1
      dimension ccNH(4)
      dimension DccNHD2(4),DccNHD3(4)

      af1 = amh/(amf+amh)
      af2 = amf/(amf+amh)

      return
c
      entry pot12(r,e,De,nt)

      do 10 i=1,nt

      hlp1 = dexp(-xpNH*r(i,2)**2)
      hlp2 = dexp(-xpNH1*r(i,3)**2)
      ccNH(1) = ccNHa(1)+ccNHi(1)*hlp1+ccNHi1(1)*hlp2
      ccNH(2) = ccNHa(2)+ccNHi(2)*hlp1+ccNHi1(2)*hlp2
      ccNH(3) = ccNHa(3)+ccNHi(3)*hlp1+ccNHi1(3)*hlp2
      ccNH(4) = ccNHa(4)+ccNHi(4)*hlp1+ccNHi1(4)*hlp2

      hlp1 = - 2.d0*xpNH*r(i,2)*hlp1
      hlp2 = - 2.d0*xpNH1*r(i,3)*hlp2
      DccNHD2(1) = ccNHi(1)*hlp1
      DccNHD3(1) = ccNHi1(1)*hlp2
      DccNHD2(2) = ccNHi(2)*hlp1
      DccNHD3(2) = ccNHi1(2)*hlp2
      DccNHD2(3) = ccNHi(3)*hlp1
      DccNHD3(3) = ccNHi1(3)*hlp2
      DccNHD2(4) = ccNHi(4)*hlp1
      DccNHD3(4) = ccNHi1(4)*hlp2

      hlp = ccNF(3)/(1.d0+ccNF(4)*r(i,3))
      unaf = dexp(-ccNF(1)-ccNF(2)*r(i,3)- hlp*r(i,3))
      unaf = unaf*27.2113961d0
      DunafD3 = -unaf*(ccNF(2)+hlp-ccNF(4)*hlp*hlp*r(i,3)/ccNF(3))


      hlp = 1.d0/(1.d0+ccNH(4)*r(i,1))
      unah = dexp(-ccNH(1)-ccNH(2)*r(i,1)- ccNH(3)*hlp*r(i,1))
      DunahD1 = -unah*(ccNH(2)+ccNH(3)*hlp
     &        -ccNH(4)*hlp*hlp*ccNH(3)*r(i,1))
      DunahD2 = -unah*(DccNHD2(1)+DccNHD2(2)*r(i,1)
     &+ DccNHD2(3)*r(i,1)*hlp-ccNH(3)*r(i,1)*r(i,1)*dccNHD2(4)*hlp*hlp)
      DunahD3 = -unah*(DccNHD3(1)+DccNHD3(2)*r(i,1)
     &+ DccNHD3(3)*r(i,1)*hlp-ccNH(3)*r(i,1)*r(i,1)*dccNHD3(4)*hlp*hlp)
C-----------------------------------------------------------------------

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)
      rnag2=af1*r1s+af2*r3s-af1*af2*r2s
      Drnag2D1 = 2.d0*af1*r(i,1)
      Drnag2D2 = -2.d0*af1*af2*r(i,2)
      Drnag2D3 = 2.d0*af2*r(i,3)
      hlp2 = rnag2+asmall
      hlp = dsqrt(hlp2)
      cstsq=0.5d0*(r3s-r1s+(af2-af1)*r2s)
     &               /(r(i,2)*hlp)
      DcstsqD1 = -2.d0*r(i,1)/(r(i,2)*hlp)
     &         - cstsq*Drnag2D1/hlp2
      DcstsqD2 = 2.d0*((af2-af1)/hlp-cstsq/r(i,2))
     &         - cstsq*Drnag2D2/hlp2
      DcstsqD3 = 2.d0*r(i,3)/(r(i,2)*hlp)
     &         - cstsq*Drnag2D3/hlp2
      cstsq = 2.0d0*(1.d0+cstsq)

            hlp = r21(2)-r21(1)-r21(3)
            r21t = r21(1)+hlp*cstsq
     &                   +r21(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*r21(3)*cstsq
            Dr21tD1 = hlp1*DcstsqD1
            Dr21tD2 = hlp1*DcstsqD2
            Dr21tD3 = hlp1*DcstsqD3

            hlp = r23(2)-r23(1)-r23(3)
            r23t = r23(1)+hlp*cstsq
     &                   +r23(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*r23(3)*cstsq
            Dr23tD1 = hlp1*DcstsqD1
            Dr23tD2 = hlp1*DcstsqD2
            Dr23tD3 = hlp1*DcstsqD3

            hlp = r13(2)-r13(1)-r13(3)
            r13t = r13(1)+hlp*cstsq
     &                   +r13(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*r13(3)*cstsq
            Dr13tD1 = hlp1*DcstsqD1
            Dr13tD2 = hlp1*DcstsqD2
            Dr13tD3 = hlp1*DcstsqD3

            hlp = ahf(2)-ahf(1)-ahf(3)
            ahft = ahf(1)+hlp*cstsq
     &                   +ahf(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*ahf(3)*cstsq
            DahftD1 = hlp1*DcstsqD1
            DahftD2 = hlp1*DcstsqD2
            DahftD3 = hlp1*DcstsqD3

            hlp = anaf(2)-anaf(1)-anaf(3)
            anaft = anaf(1)+hlp*cstsq
     &                   +anaf(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*anaf(3)*cstsq
            DanaftD1 = hlp1*DcstsqD1
            DanaftD2 = hlp1*DcstsqD2
            DanaftD3 = hlp1*DcstsqD3

            hlp = anah(2)-anah(1)-anah(3)
            anaht = anah(1)+hlp*cstsq
     &                   +anah(3)*cstsq*cstsq
            hlp1 = hlp+2.d0*anah(3)*cstsq
            DanahtD1 = hlp1*DcstsqD1
            DanahtD2 = hlp1*DcstsqD2
            DanahtD3 = hlp1*DcstsqD3
 
      hlp = anaht*r(i,1)-anaft*r(i,3)-r13t
      sw1 = 0.5d0*(1.d0+dtanh(hlp))
      hlp1 = 0.5d0/(dcosh(hlp)**2)
      Dsw1D1 = hlp1*(anaht+DanahtD1*r(i,1)-DanaftD1*r(i,3)-Dr13tD1)
      Dsw1D2 = hlp1*(DanahtD2*r(i,1)-DanaftD2*r(i,3)-Dr13tD2)
      Dsw1D3 = hlp1*(DanahtD3*r(i,1)-anaft-DanaftD3*r(i,3)-Dr13tD3)

      hlp = ahft*r(i,2)-anaht*r(i,1)-r21t
      sw21 = 0.5d0*(1.d0+dtanh(hlp))
      hlp1 = 0.5d0/(dcosh(hlp)**2)
      Dsw21D1 = hlp1*(DahftD1*r(i,2)-anaht-DanahtD1*r(i,1)-Dr21tD1)
      Dsw21D2 = hlp1*(ahft+DahftD2*r(i,2)-DanahtD2*r(i,1)-Dr21tD2)
      Dsw21D3 = hlp1*(DahftD3*r(i,2)-DanahtD3*r(i,1)-Dr21tD3)

      hlp = ahft*r(i,2)-anaft*r(i,3)-r23t
      sw23 = 0.5d0*(1.d0+dtanh(hlp))
      hlp1 = 0.5d0/(dcosh(hlp)**2)
      Dsw23D1 = hlp1*(DahftD1*r(i,2)-DanaftD1*r(i,3)-Dr23tD1)
      Dsw23D2 = hlp1*(ahft+DahftD2*r(i,2)-DanaftD2*r(i,3)-Dr23tD2)
      Dsw23D3 = hlp1*(DahftD3*r(i,2)-anaft-DanaftD3*r(i,3)-Dr23tD3)

      e(i) = unaf*sw1*sw23+unah*(1.d0-sw1)*sw21
      De(1,i) = unaf*Dsw1D1*sw23+unaf*sw1*Dsw23D1
     &        + DunahD1*(1.d0-sw1)*sw21-unah*Dsw1D1*sw21
     &        + unah*(1.d0-sw1)*Dsw21D1
      De(2,i) = unaf*Dsw1D2*sw23+unaf*sw1*Dsw23D2
     &        + DunahD2*(1.d0-sw1)*sw21-unah*Dsw1D2*sw21
     &        + unah*(1.d0-sw1)*Dsw21D2
      De(3,i) = DunafD3*sw1*sw23+unaf*Dsw1D3*sw23+unaf*sw1*Dsw23D3
     &        + DunahD3*(1.d0-sw1)*sw21-unah*Dsw1D3*sw21
     &        + unah*(1.d0-sw1)*Dsw21D3

      e(i) = e(i)*cevau
      De(1,i) = De(1,i)*cevau
      De(2,i) = De(2,i)*cevau
      De(3,i) = De(3,i)*cevau
10    continue

      return

  100 format(/25('='),' NaFH - U12 coupling ',25('='),/
     &2x,'mF/(mF+mH) = ',f13.5,2x,'mH/(mF+mH) = ',f13.5)
  101 format(2x,'ccNaF = ',5x,4e14.6/10x,2e14.6)
  102 format(2x,'ccNaHa = ',4x,4e14.6/11x,2e14.6)
  103 format(2x,'ccNaHi = ',4x,4e14.6/11x,2e14.6)
  104 format(2x,'ccNaHi1 = ',3x,4e14.6/12x,2e14.6)
  105 format(2x,'xpNH,xpNH1 = ',2e14.6)
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
      dimension De(3,nt)
      common /com_para/ g_myrank,g_nprocs
      integer g_myrank,g_nprocs

      save /onetwo/,/utwo/,/angular2/,rac2,af1,af2

      common/onetwo/CNF(6),CNF1(5),CHF(6),bc(0:5),
     & ads1,ars1,abs1,amh,amf,amna,amu,una2p,asmall,cevau,
     & DeNF,reNF,DeHF,reHF
      common/utwo/y(4,3,4),c2a(4),c2b(4),c2c(4),r20i(4),
     &      r20a,al20i(4),al20i1(4),al20a,psi,rsi,al2a(4),
     &      rr0a(4),dip2(4),adip2(4),rdip2(4),acut2(4),
     &      rcut2(4),ry(4),aly(4),dy(4),dy1(4),dy2(4),
     &      dip2n(4),adip2n(4),rdip2n(4),acut2n(4),rcut2n(4)
      common/angular2/ac(4)
C-----------------------------------------------------------------------

      rac2=1.d0/dsqrt(2.d0)

      af1 = amh/(amf+amh)
      af2 = amf/(amf+amh)
C
      return
c
      entry pot22(r,e,De,nt)
C
      do 10 i=1,nt

      r1s=r(i,1)*r(i,1)
      r2s=r(i,2)*r(i,2)
      r3s=r(i,3)*r(i,3)
      rnag2=af1*r1s+af2*r3s-af1*af2*r2s
      Drnag2D1 = 2.d0*af1*r(i,1)
      Drnag2D2 = -2.d0*af1*af2*r(i,2)
      Drnag2D3 = 2.d0*af2*r(i,3)

C     Cosine of Na-F-H bond angle:
      hlp = 1.d0/(r(i,2)*r(i,3))
      csb = 0.5d0*(r2s+r3s-r1s)*hlp
      DcsbD1 = -2.d0*r(i,1)*hlp
      DcsbD2 = 2.d0*(1.d0/r(i,3)-csb/r(i,2))
      DcsbD3 = 2.d0*(1.d0/r(i,2)-csb/r(i,3))
      csb = 2.d0*(1.d0+csb)

      if (csb.ge.4.d0) csb = 3.99999d0
      l1 = idint(csb)+1
      l2 = l1+1
      if (l1.eq.4) l2=l1
      if (csb.le.bc(l1)) csb = bc(l1)+1.d-13
      if (csb.ge.bc(l1+1)) csb = bc(l1+1)-1.d-13

C-----------------------------------------------------------------------
      s1=0.d0
      Ds1D1 = 0.d0
C-----------------------------------------------------------------------
      y2 = r(i,2)-reHF
      hlp1 = cHF(1)*dexp(-cHF(2)*y2)
      hlp2 = dexp(-cHF(6)*y2)
      hlp22 = (-1.d0-cHF(1)+cHF(3)*y2+cHF(4)*y2**2+cHF(5)*y2**3)*hlp2 
      s2a = DeHF*(hlp1+hlp22) + una2p
      Ds2aD2 = DeHF*(-cHF(2)*hlp1-cHF(6)*hlp22
     &       + (cHF(3)+2.d0*cHF(4)*y2+3.d0*cHF(5)*y2**2)*hlp2)
C-----------------------------------------------------------------------
      y3 = r(i,3)-reNF
      hlp1 = cNF(1)*dexp(-cNF(2)*y3)
      hlp2 = dexp(-cNF(6)*y3)
      hlp22 = (-1.d0-cNF(1)+cNF(3)*y3+cNF(4)*y3*y3+cNF(5)*y3**3)*hlp2
      s3 = DeNF*(hlp1+hlp22)
      Ds3D3 = DeNF*(-cNF(2)*hlp1-cNF(6)*hlp22
     &       + (cNF(3)+2.d0*cNF(4)*y3+3.d0*cNF(5)*y3**2)*hlp2)
C-----------------------------------------------------------------------

      e(i) = 0.d0
      De(1,i) = 0.d0
      De(2,i) = 0.d0
      De(3,i) = 0.d0
      anorm = 0.d0
      DanormD1 = 0.d0
      DanormD2 = 0.d0
      DanormD3 = 0.d0

      do 20 l = l1, l2

      hlp1 = y(1,1,l)*dexp(-y(2,1,l)*r(i,1))
      hlp2 = y(3,1,l)*dexp(-y(4,1,l)*r(i,1))
      t1=hlp1+hlp2+s1
      Dt1D1 = -y(2,1,l)*hlp1-y(4,1,l)*hlp2+Ds1D1

      hlp = al2a(l)*(r(i,3)+r(i,1)-rr0a(l))
      hlp1 = al2a(l)/(dcosh(hlp)**2)
      hlp = 0.5d0*(1.d0+dtanh(hlp))
      hlp2 = r(i,3)+r(i,1)-rsi
      if (hlp2.gt.0.d0) then
c     write(6,*)'R1+R3=',(r(i,1)+r(i,3))
      ctanh = dexp(-psi/hlp2)
      dctanh = psi*ctanh/(r(i,3)+r(i,1)-rsi)**2
      else
      ctanh = 0.d0
      dctanh = 0.d0
      end if
      r20 = r20i(l)+(r20a-r20i(l))*hlp
      Dr20 = 0.5d0*(r20a-r20i(l))*hlp1
      al20 = al20i(l)+(al20i1(l)-al20i(l))*hlp
     &     + (al20a-al20i1(l))*ctanh
      Dal20 = 0.5d0*(al20i1(l)-al20i(l))*hlp1
     &      + (al20a-al20i1(l))*dctanh
   
      hlp1 = y(1,2,l)*dexp(-y(2,2,l)*r(i,2))
      hlp2 =  dexp(-(y(4,2,l)+dy2(l)*ctanh)*r(i,2))
      hlp3 = (y(3,2,l)+dy1(l)*ctanh)*hlp2
      t2=hlp1+hlp3
      Dt2D1=dy1(l)*dctanh*hlp2 - dy2(l)*dctanh*r(i,2)*hlp3
      Dt2D2=-y(2,2,l)*hlp1-(y(4,2,l)+dy2(l)*ctanh)*hlp3
      Dt2D3=Dt2D1

      hlp = 0.5d0*(1.d0+dtanh(al20*(r(i,2)-r20)))
      s2=s2a+(t2-s2a)*hlp
      hlp1 = 0.5d0/(dcosh(al20*(r(i,2)-r20))**2)
      Ds2D1 = (t2-s2a)*hlp1*(Dal20*(r(i,2)-r20)-al20*Dr20)
      Ds2D2 = Ds2aD2+(Dt2D2-Ds2aD2)*hlp
     &      + al20*hlp1*(t2-s2a)
      Ds2D3 = Ds2D1

      hlp = aly(l)*(r(i,3)-ry(l))
      cy=y(2,3,l)+0.5d0*dy(l)*(1.d0+dtanh(hlp))
      dcyD3 = 0.5d0*dy(l)*aly(l)/(dcosh(hlp)**2)
      hlp1 = y(1,3,l)*dexp(-cy*r(i,3))
      hlp2 = y(3,3,l)*dexp(-y(4,3,l)*r(i,3))
      t3 = hlp1+hlp2
      Dt3D3=-(cy+dcyD3*r(i,3))*hlp1-y(4,3,l)*hlp2

      coul1=0.5d0*(s1+t1)
      Dcl1D1=0.5d0*(Ds1D1+Dt1D1)
      coul2=0.5d0*(s2+t2)
      Dcl2D1=0.5d0*(Ds2D1+Dt2D1)
      Dcl2D2=0.5d0*(Ds2D2+Dt2D2)
      Dcl2D3=0.5d0*(Ds2D3+Dt2D3)
      coul3=0.5d0*(s3+t3)
      Dcl3D3=0.5d0*(Ds3D3+Dt3D3)

      exch1=0.5d0*(s1-t1)
      Dex1D1=0.5d0*(Ds1D1-Dt1D1)
      exch2=0.5d0*(s2-t2)
      Dex2D1=0.5d0*(Ds2D1-Dt2D1)
      Dex2D2=0.5d0*(Ds2D2-Dt2D2)
      Dex2D3=0.5d0*(Ds2D3-Dt2D3)
      exch3=0.5d0*(s3-t3)
      Dex3D3=0.5d0*(Ds3D3-Dt3D3)

      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2
      DwD1 = 2.d0*(exch1-exch2)*(Dex1D1-Dex2D1)
     &     + 2.d0*(exch2-exch3)*Dex2D1
     &     - 2.d0*(exch3-exch1)*Dex1D1
      DwD2 = - 2.d0*(exch1-exch2)*Dex2D2
     &     + 2.d0*(exch2-exch3)*Dex2D2
      DwD3 = - 2.d0*(exch1-exch2)*Dex2D3
     &     + 2.d0*(exch2-exch3)*(Dex2D3-Dex3D3)
     &     + 2.d0*(exch3-exch1)*Dex3D3

      cplg2=c2a(l)*dexp(-c2b(l)*w-c2c(l)*(r(i,1)+r(i,2)+r(i,3)))
      DcD1 = (-c2b(l)*DwD1 - c2c(l))*cplg2
      DcD2 = (-c2b(l)*DwD2 - c2c(l))*cplg2
      DcD3 = (-c2b(l)*DwD3 - c2c(l))*cplg2

      rootw = rac2*dsqrt(w+cplg2*cplg2)
      surf=(coul1+coul2+coul3+DeHF-rootw)
      DsurfD1=(Dcl1D1+Dcl2D1-0.25d0*(DwD1+2.d0*cplg2*DcD1)/rootw)
      DsurfD2=(Dcl2D2-0.25d0*(DwD2+2.d0*cplg2*DcD2)/rootw)
      DsurfD3=(Dcl2D3+Dcl3D3-0.25d0*(DwD3+2.d0*cplg2*DcD3)/rootw)

      hlp1 = acut2(l)*(r(i,2)-rcut2(l))
      hlp2 = adip2(l)*(rnag2-rdip2(l)**2)
      hlp11 = dtanh(hlp1)
      hlp22 = dtanh(hlp2)
      hlp1 = 1.d0/dcosh(hlp1)
      hlp2 = 1.d0/dcosh(hlp2)
      surf = surf - dip2(l)*hlp1*hlp2
      DsurfD1=DsurfD1 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*Drnag2D1
      DsurfD2=DsurfD2 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*Drnag2D2
     &                + dip2(l)*hlp2*hlp1*hlp11*acut2(l)
      DsurfD3=DsurfD3 + dip2(l)*hlp1*hlp2*hlp22*adip2(l)*Drnag2D3

      hlp1 = acut2n(l)*(r(i,2)-rcut2n(l))
      hlp2 = adip2n(l)*(rnag2-rdip2n(l)**2)
      hlp11 = dtanh(hlp1)
      hlp22 = dtanh(hlp2)
      hlp1 = 1.d0/dcosh(hlp1)
      hlp2 = 1.d0/dcosh(hlp2)
      surf = surf - dip2n(l)*hlp1*hlp2
      DsurfD1=DsurfD1 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*Drnag2D1
      DsurfD2=DsurfD2 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*Drnag2D2
     &                + dip2n(l)*hlp2*hlp1*hlp11*acut2n(l)
      DsurfD3=DsurfD3 + dip2n(l)*hlp1*hlp2*hlp22*adip2n(l)*Drnag2D3

      hlp = 1.d0/(
     &              (1.d0-(csb-bc(l))/(bc(l-1)-bc(l)))
     &            * (1.d0-(csb-bc(l))/(bc(l+1)-bc(l)))
     &                                                 )
      hlp1 = (csb-bc(l))**2*hlp
      cog = dexp(-ac(l)*hlp1)
      Dcog = -cog*ac(l)*(2.d0*(csb-bc(l))*hlp
     & +hlp1*(1.d0/(bc(l-1)-csb)+1.d0/(bc(l+1)-csb)))
      DcogD1 = Dcog*DcsbD1
      DcogD2 = Dcog*DcsbD2
      DcogD3 = Dcog*DcsbD3

      e(i) = e(i) + surf*cog
      De(1,i) = De(1,i) + DsurfD1*cog+surf*DcogD1
      De(2,i) = De(2,i) + DsurfD2*cog+surf*DcogD2
      De(3,i) = De(3,i) + DsurfD3*cog+surf*DcogD3
      anorm = anorm + cog
      DanormD1 = DanormD1 + DcogD1
      DanormD2 = DanormD2 + DcogD2
      DanormD3 = DanormD3 + DcogD3
   20 continue
      hlp = 1.d0/anorm
      e(i) = e(i)*hlp
      De(1,i) = (De(1,i) - e(i)*DanormD1)*hlp
      De(2,i) = (De(2,i) - e(i)*DanormD2)*hlp
      De(3,i) = (De(3,i) - e(i)*DanormD3)*hlp
      e(i) = e(i)*cevau
      De(1,i) = De(1,i)*cevau
      De(2,i) = De(2,i)*cevau
      De(3,i) = De(3,i)*cevau

10    continue

      return

  100 format(/26('='),' NaFH - U22 diabat ',26('='),
     &      /2x,'mF/(mF+mH) = ',f13.5,2x,'mH/(mF+mH) = ',f13.5)
  101 format(2x,'DeNaF = ',9x,e14.6/2x,'reNaF = ',9x,e14.6,
     & /2x,'cNaF = ',10x,3e14.6/19x,3e14.6)
  103 format(2x,'DeHF = ',10x,e14.6/2x,'reHF = ',10x,e14.6
     & /2x,'cHF = ',11x,3e14.6/19x,3e14.6)
  105 format(2x,'una2p,asmall = ',2x,2e14.6)
  106 format(2x,'y(1,1,i) = ',2x,4e13.5)
  107 format(2x,'y(2,1,i) = ',2x,4e13.5)
  108 format(2x,'y(3,1,i) = ',2x,4e13.5)
  109 format(2x,'y(4,1,i) = ',2x,4e13.5)
  110 format(2x,'y(1,2,i) = ',2x,4e13.5)
  111 format(2x,'y(2,2,i) = ',2x,4e13.5)
  112 format(2x,'y(3,2,i) = ',2x,4e13.5)
  113 format(2x,'y(4,2,i) = ',2x,4e13.5)
  114 format(2x,'y(1,3,i) = ',2x,4e13.5)
  115 format(2x,'y(2,3,i) = ',2x,4e13.5)
  116 format(2x,'y(3,3,i) = ',2x,4e13.5)
  117 format(2x,'y(4,3,i) = ',2x,4e13.5)
  118 format(2x,'c2a = ',7x,4e13.5)
  119 format(2x,'c2b = ',7x,4e13.5)
  120 format(2x,'c2c = ',7x,4e13.5)
  121 format(2x,'r20i = ',6x,4e13.5)
  122 format(2x,'r20a = ',6x,4e13.5)
  123 format(2x,'al20i = ',5x,4e13.5)
  124 format(2x,'al20i1 = ',4x,4e13.5)
  125 format(2x,'al20a = ',5x,4e13.5)
  126 format(2x,'al2a = ',6x,4e13.5)
  127 format(2x,'rr0a = ',6x,4e13.5)
  128 format(2x,'psi = ',7x,4e13.5)
  129 format(2x,'rsi = ',7x,4e13.5)
  130 format(2x,'dip2 = ',6x,4e13.5)
  131 format(2x,'adip2 = ',5x,4e13.5)
  132 format(2x,'rdip2 = ',5x,4e13.5)
  133 format(2x,'acut2 = ',5x,4e13.5)
  134 format(2x,'rcut2 = ',5x,4e13.5)
  135 format(2x,'ry = ',8x,4e13.5)
  136 format(2x,'aly = ',7x,4e13.5)
  137 format(2x,'dy = ',8x,4e13.5)
  138 format(2x,'dip2n = ',5x,4e13.5)
  139 format(2x,'adip2n = ',4x,4e13.5)
  140 format(2x,'rdip2n = ',4x,4e13.5)
  141 format(2x,'acut2n = ',4x,4e13.5)
  142 format(2x,'rcut2n = ',4x,4e13.5
     &    /71('=')/)

      end

      BLOCK DATA GENERAL
      double precision bc
      double precision CNF,CNF1,CHF
      double precision amh,amf,amna,amu
      double precision DeNF,reNF,DeHF,reHF
      double precision ads1,ars1,abs1
      double precision una2p,asmall,cevau

      common/onetwo/CNF(6),CNF1(5),CHF(6),bc(0:5),
     & ads1,ars1,abs1,amh,amf,amna,amu,una2p,asmall,cevau,
     & DeNF,reNF,DeHF,reHF
      data amh/1.007825d0/,amf/18.998403d0/,amna/22.989767d0/,
     &     amu/1822.8885d0/
      data DeNF/4.3830525d0/,reNF/3.6395d0/,
     & CNF/0.168364d0,2.29333d0,-0.848620d-1,-0.109707d-1,
     &   0.246582d-2,0.377012d0/
      data cNF1/0.620115d0,0.694459d-1,0.757570d0,0.283055d0,2.69998d0/
      data DeHF/5.77096d0/,reHF/1.7328d0/,cHF/0.441258d0,3.88056d0,
     & -1.27110d0,-1.37766d0,0.168186d0,2.07230d0/
      data bc/-1.d0,0.d0,1.d0,2.d0,3.0,5.d0/
      data ads1/2.91087d-3/,ars1/7.37707d0/,abs1/7.13223d-1/
      data una2p/2.0973375d0/
      data asmall/0.09d0/
      data cevau/0.036749308867692d0/

      END

      BLOCK DATA UONEONE
      double precision ac
      double precision y,dip1,adip1,rdip1,acut1,rcut1,
     &            c2a,c2b,c2c,dipn,al2,aln0,
     &            dipnn,al2n,aln0n1,aln0n2,rn0n
      common/uone/y(4,3,4),c2a(4),c2b(4),c2c(4),dip1(4),
     &            adip1(4),rdip1(4),acut1(4),rcut1(4),
     &            dipn(4),al2(4),aln0(4),
     &            dipnn(4),al2n(4),aln0n1(4),aln0n2(4),rn0n(4)
      common/angular1/ac(4)
      data y(1,1,1)/0.00000d0/,y(2,1,1)/0.00000d0/,
     &     y(1,1,2)/0.00000d0/,y(2,1,2)/0.00000d0/,
     &     y(1,1,3)/0.00000d0/,y(2,1,3)/0.00000d0/,
     &     y(1,1,4)/0.00000d0/,y(2,1,4)/0.00000d0/,
     &     y(3,1,1)/6.6796d0/,y(4,1,1)/2.7839d0/,
     &     y(3,1,2)/6.6796d0/,y(4,1,2)/2.7839d0/,
     &     y(3,1,3)/7.1360d0/,y(4,1,3)/1.0376d0/,
     &     y(3,1,4)/8.8737d0/,y(4,1,4)/1.9627d0/
      data y(1,2,1)/6.0459d+1/,y(2,2,1)/3.7503d0/,
     &     y(1,2,2)/6.0459d+1/,y(2,2,2)/3.7503d0/,
     &     y(1,2,3)/6.0459d+1/,y(2,2,3)/3.7503d0/,
     &     y(1,2,4)/8.5430d0/,y(2,2,4)/9.00856d-1/,
     &     y(3,2,1)/8.4874d0/,y(4,2,1)/1.0916d+1/,
     &     y(3,2,2)/8.4874d0/,y(4,2,2)/9.9850d0/,
     &     y(3,2,3)/8.4874d0/,y(4,2,3)/5.5898d0/,
     &     y(3,2,4)/9.5899d0/,y(4,2,4)/1.1615d+1/
      data y(1,3,1)/7.3251d0/,y(2,3,1)/2.8774d0/,
     &     y(1,3,2)/7.3251d0/,y(2,3,2)/3.6094d0/,
     &     y(1,3,3)/7.3251d0/,y(2,3,3)/3.9288d0/,
     &     y(1,3,4)/7.5811d0/,y(2,3,4)/4.4272d0/,
     &     y(3,3,1)/1.9300d+3/,y(4,3,1)/2.7091d0/,
     &     y(3,3,2)/3.8969d+3/,y(4,3,2)/2.7591d0/,
     &     y(3,3,3)/5.9130d+3/,y(4,3,3)/2.6881d0/,
     &     y(3,3,4)/5.9007d+3/,y(4,3,4)/3.8929d0/
      data c2a(1)/8.5742d0/,c2b(1)/1.5953d0/,c2c(1)/3.3745d-1/
      data c2a(2)/9.8509d0/,c2b(2)/6.6359d0/,c2c(2)/3.0396d-1/
      data c2a(3)/9.5937d0/,c2b(3)/1.3441d+1/,c2c(3)/2.8812d-1/
      data c2a(4)/1.1386d0/,c2b(4)/8.1112d0/,c2c(4)/2.6061d0/
      data dip1(1)/1.9651d0/,rdip1(1)/8.3109d-1/,
     &     dip1(2)/1.9651d0/,rdip1(2)/8.3109d-1/,
     &     dip1(3)/2.1803d0/,rdip1(3)/7.1647d-1/,
     &     dip1(4)/2.2483d0/,rdip1(4)/9.7128d-1/,
     &     adip1(1)/4.8674d-2/,adip1(2)/5.0603d-2/,
     &     adip1(3)/4.8499d-2/,adip1(4)/8.5924d-2/
      data acut1(1)/6.3761d0/,rcut1(1)/5.2007d0/,
     &     acut1(2)/6.3761d0/,rcut1(2)/3.6585d0/,
     &     acut1(3)/6.0812d0/,rcut1(3)/3.4684d0/,
     &     acut1(4)/5.0516d0/,rcut1(4)/1.6434d0/
      data dipn(1)/1.5754d-1/,al2(1)/8.4809d-2/,
     &     dipn(2)/1.2274d-1/,al2(2)/8.4809d-2/,
     &     dipn(3)/4.1823d-2/,al2(3)/9.4867d-2/,
     &     dipn(4)/3.4350d-1/,al2(4)/6.5542d-2/,
     &      aln0(1)/3.1971d-2/,
     &      aln0(2)/3.5817d-2/,
     &      aln0(3)/1.6686d-1/,
     &      aln0(4)/5.9961d-2/
      data dipnn(1)/0.0000d0/,al2n(1)/0.0000d-2/,
     &     dipnn(2)/3.7501d-3/,al2n(2)/9.2376d0/,
     &     dipnn(3)/7.6323d-3/,al2n(3)/9.3515d-2/,
     &     dipnn(4)/0.0000d0/,al2n(4)/.0000d-2/,
     &      rn0n(1)/0.0d0/,
     &      rn0n(2)/6.6001d0/,
     &      rn0n(3)/5.9949d0/,
     &      rn0n(4)/0.0d0/,
     &      aln0n1(1)/0.0000d-2/,aln0n2(1)/0.0000d-2/,
     &      aln0n1(2)/3.2545d-2/,aln0n2(2)/4.9001d-1/,
     &      aln0n1(3)/8.5914d-2/,aln0n2(3)/2.1382d0/,
     &      aln0n1(4)/0.0000d-2/,aln0n2(4)/0.0000d-2/
      data ac/4.6913d0,5.4684d-1,7.5932d-1,8.0000d-1/

      END

      BLOCK DATA UTWOTWO
      double precision ac
      double precision y,r20i,r20a,al20i,al20a,al2a,
     &            rr0a,ry,aly,dy,c2a,c2b,c2c,al20i1,
     &            dip2,adip2,rdip2,acut2,rcut2,psi,rsi,
     &            dip2n,adip2n,rdip2n,acut2n,rcut2n,
     &            dy1,dy2
      common/utwo/y(4,3,4),c2a(4),c2b(4),c2c(4),r20i(4),
     &      r20a,al20i(4),al20i1(4),al20a,psi,rsi,al2a(4),
     &      rr0a(4),dip2(4),adip2(4),rdip2(4),acut2(4),
     &      rcut2(4),ry(4),aly(4),dy(4),dy1(4),dy2(4),
     &      dip2n(4),adip2n(4),rdip2n(4),acut2n(4),rcut2n(4)
      common/angular2/ac(4)

      data y(1,1,1)/1.0407d+3/,y(2,1,1)/4.5911d0/,
     &     y(1,1,2)/1.1608d+3/,y(2,1,2)/3.8886d0/,
     &     y(1,1,3)/4.9492d+3/,y(2,1,3)/2.6147d0/,
     &     y(1,1,4)/8.3093d+3/,y(2,1,4)/3.2431d0/,
     &     y(3,1,1)/1.4941d0/,y(4,1,1)/4.5525d0/,
     &     y(3,1,2)/1.4941d0/,y(4,1,2)/4.5525d0/,
     &     y(3,1,3)/1.4941d0/,y(4,1,3)/4.5525d0/,
     &     y(3,1,4)/2.8361d0/,y(4,1,4)/4.4686d0/
      data y(1,2,1)/0.0000d0/,y(2,2,1)/0.0000d0/,
     &     y(1,2,2)/0.0000d0/,y(2,2,2)/0.0000d0/,
     &     y(1,2,3)/0.0000d0/,y(2,2,3)/0.0000d0/,
     &     y(1,2,4)/0.0000d0/,y(2,2,4)/0.0000d0/,
     &     y(3,2,1)/1.8076d+3/,y(4,2,1)/2.8815d0/
     &     y(3,2,2)/1.8076d+3/,y(4,2,2)/2.8815d0/,
     &     y(3,2,3)/1.8076d+3/,y(4,2,3)/2.8815d0/,
     &     y(3,2,4)/2.3109d+3/,y(4,2,4)/3.1315d0/
      data y(1,3,1)/5.7330d0/,y(2,3,1)/3.2053d-1/,
     &     y(1,3,2)/5.7330d0/,y(2,3,2)/3.2053d-1/,
     &     y(1,3,3)/5.7330d0/,y(2,3,3)/3.2053d-1/,
     &     y(1,3,4)/5.0393d0/,y(2,3,4)/2.7567d-1/,
     &     y(3,3,1)/2.5774d+5/,y(4,3,1)/4.2267d0/,
     &     y(3,3,2)/2.5774d+5/,y(4,3,2)/4.2267d0/,
     &     y(3,3,3)/2.5774d+5/,y(4,3,3)/4.2267d0/,
     &     y(3,3,4)/1.6506d+5/,y(4,3,4)/3.9246d0/

      data c2a(1)/5.6406d0/,c2b(1)/7.7152d-2/
      data c2a(2)/8.2390d0/,c2b(2)/1.7133d-1/
      data c2a(3)/7.1504d0/,c2b(3)/1.1011d-1/
      data c2a(4)/8.0196d0/,c2b(4)/1.0017d-1/,
     &                         c2c(1)/8.4595d-2/,
     &                         c2c(2)/8.4595d-2/,
     &                         c2c(3)/8.4595d-2/,
     &                         c2c(4)/8.4411d-2/
      data r20i(1)/2.9013d0/,al20i(1)/1.7119d0/
      data r20i(2)/3.1196d0/,al20i(2)/9.4849d-1/
      data r20i(3)/3.0995d0/,al20i(3)/0.0d0/
      data r20i(4)/2.6572d0/,al20i(4)/2.3214d0/,
     &     r20a/3.1305d0/,al20i1(1)/3.5428d0/,
     &     al20a/3.024d+1/,al20i1(2)/6.7257d0/,
     &     psi/0.7d0/,al20i1(3)/3.0240d+1/,
     &     rsi/22.2d0/,al20i1(4)/9.6978d0/,
     &     al2a(1)/5.8088d-1/,rr0a(1)/9.5543d0/,
     &     al2a(2)/5.8088d-1/,rr0a(2)/9.5543d0/,
     &     al2a(3)/5.8088d-1/,rr0a(3)/9.5543d0/,
     &     al2a(4)/8.1525d-1/,rr0a(4)/8.2500d0/
      data dip2(1)/2.6475d-1/
      data dip2(2)/2.6475d-1/
      data dip2(3)/2.6475d-1/
      data dip2(4)/5.0336d-1/
      data adip2(1)/3.2317d-1/,rdip2(1)/3.7260d0/
      data adip2(2)/1.8513d-1/,rdip2(2)/4.5515d0/
      data adip2(3)/4.3999d-2/,rdip2(3)/2.1191d-1/
      data adip2(4)/2.3003d-1/,rdip2(4)/5.0160d0/
      data acut2(1)/2.0799d0/,rcut2(1)/1.2496d0/
      data acut2(2)/1.4204d0/,rcut2(2)/1.0817d0/
      data acut2(3)/7.9086d-1/,rcut2(3)/1.0575d0/
      data acut2(4)/2.8447d0/,rcut2(4)/3.0015d0/
      data dip2n(1)/2.2091d-1/,adip2n(1)/1.0288d-1/
      data dip2n(2)/2.2091d-1/,adip2n(2)/1.0288d-1/
      data dip2n(3)/2.2091d-1/,adip2n(3)/1.0288d-1/
      data dip2n(4)/4.5301d-1/,adip2n(4)/1.1009d-1/
      data rdip2n(1)/5.1554d0/
      data rdip2n(2)/5.1554d0/
      data rdip2n(3)/5.1554d0/
      data rdip2n(4)/6.7003d0/
      data acut2n(1)/3.5607d0/,rcut2n(1)/3.0661d0/
      data acut2n(2)/4.3956d0/,rcut2n(2)/3.0661d0/
      data acut2n(3)/4.1924d0/,rcut2n(3)/3.0661d0/
      data acut2n(4)/3.0765d0/,rcut2n(4)/2.8845d0/
      data ry(1)/1.4328d+1/,aly(1)/4.9044d-1/
      data ry(2)/1.2504d+1/,aly(2)/6.4718d-1/
      data ry(3)/1.2102d+1/,aly(3)/2.5833d-1/
      data ry(4)/7.9381d0/,aly(4)/6.4853d-1/
      data dy(1)/7.5126d-1/,dy1(1)/0.d0/,dy2(1)/0.d0/
      data dy(2)/2.0210d-1/,dy1(2)/0.d0/,dy2(2)/0.d0/
      data dy(3)/1.1749d-1/,dy1(3)/0.d0/,dy2(3)/0.d0/
      data dy(4)/8.1959d-2/,
     &               dy1(4)/-5.03d+2/,dy2(4)/-0.25d0/
      data ac/2.0018d0,4.1908d-1,5.4988d0,3.0178d-1/

      END

      BLOCK DATA UONETWO
      double precision ahf,anaf,anah,r21,r23,r13
      double precision ccNF,ccNHa,ccNHi,ccNHi1,xpNH,xpNH1
      common/ucoupl/ahf(3),anaf(3),anah(3),r21(3),r23(3),r13(3)
      common/ucoupl1/ccNF(4),ccNHa(4),ccNHi(4),ccNHi1(4)
      common/u12nah/xpNH,xpNH1
      data ccNF/3.9960D0,5.4730D-01,-1.1428D0,2.0440D-01/
      data ccNHa/1.00d0,0.80d0,-2.67d0,0.456d0/
      data ccNHi/1.0897d0,-6.5205d-1,7.7968d-1,-1.2933d-2/
      data ccNHi1/-6.1073d-1,-1.7425d-1,-8.3877d-1,
     & -3.5361d-1/
      data xpNH/5.6906d-2/,xpNH1/2.3185d-1/

      data ahf/1.7177d0,2.2707d0,1.2248d-1/,
     &    anaf/6.7786d-1,9.8809d-1,-1.0541d-1/,
     &    anah/6.5680d-1,1.0141d0,-5.1056d-2/
      data r21/3.5114d0,1.3540d0,3.2267d-1/,
     &     r23/-3.3006d-1,2.0804d-1,2.4499d-1/,
     &     r13/1.2431d0,9.9188d-1,1.1051d-1/
      END

