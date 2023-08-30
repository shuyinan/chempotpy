  !********************************************************************** 
  !   System:                     OHBr_3App
  !   Functional form:            permutationally invariant polynomials
  !   Common name:                OHBr (3A")
  !   Number of derivatives:      1
  !   Number of bodies:           3
  !   Number of electronic surfaces: 1
  !   Interface: Section-2
  !**********************************************************************
  !
  !   References: A. G. S. de Oliveira-Filho, F. R. Ornellas and K. A. 
  !   Peterson  "Accurate ab initio potential energy surfaces for the
  !   3A" and 3A' electronic states of the O(3P)+HBr system" 
  !   J. Chem. Phys. 136, 174316 (2012).
  !
  !   Notes:
  !      
  !      r1 = R(O-H)
  !      r2 = R(H-Br)
  !      r3 = R(Br-O)
  !
  !   Input: r1,r2,r3                     in units of bohr
  !   Output: value                       in units of hartree
  !   Output: dr1,dr2,dr3                 in units of hartree/bohr
  !**********************************************************************

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

      double precision :: v
      double precision :: r_OH, r_HBr, r_OBr
      double precision :: dr_OH, dr_HBr, dr_OBr
      double precision :: r(3), tx(9)
      double precision :: dr(3), drdx(3,9), dx(9)
      integer :: iatom, idir, i, j, istate
      logical, save :: first_time_data=.true.
      !initialize
      v=0.d0
      g=0.d0
      d=0.d0

      ! x: H, O, Br order
      r_OBr=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2+&
               &(x(2,3)-x(3,3))**2)
!      r_OBr=r_OBr/0.529177211

      r_OH=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2+&
               &(x(1,3)-x(2,3))**2)
!      r_OH=r_OH/0.529177211

      r_HBr=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2+&
               &(x(1,3)-x(3,3))**2)
!      r_HBr=r_HBr/0.529177211

      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        tx(j)=x(iatom,idir)
      enddo
      enddo 
      r(1)=r_OH
      r(2)=r_HBr
      r(3)=r_OBr

      if(first_time_data) then
      call myprepot(path)
      first_time_data=.false.
      endif

      call mypot(r_OH,r_HBr,r_OBr,v,dr_OH,dr_HBr,dr_OBr)

      v=v*27.211386
      dr(1)=dr_OH*27.211386
      dr(2)=dr_HBr*27.211386
      dr(3)=dr_OBr*27.211386
!      dr(1)=dr_OH*51.422067
!      dr(2)=dr_HBr*51.422067
!      dr(3)=dr_OBr*51.422067

      call evdrdx(tx, r, drdx)
 
      dx=0.d0
      do i=1,9
      do j=1,3
        dx(i)=dx(i)+dr(j)*drdx(j,i)
      enddo
      enddo

      if (igrad==0) then 
        do istate=1,nstates
          p(istate)=v
        enddo
      else if (igrad==1) then 
        do istate=1,nstates
          p(istate)=v
        enddo
        do iatom=1, natoms
        do idir=1,3
          j=(iatom-1)*3+idir
          g(1,iatom,idir)=dx(j)
        enddo
        enddo
      else if (igrad==2) then
        write (*,*) 'Only energy and gradient are available'
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

      dRdX(2,1)=(x(1)-x(7))/r(2)
      dRdX(2,2)=(x(2)-x(8))/r(2)
      dRdX(2,3)=(x(3)-x(9))/r(2)
      dRdX(2,7)=-dRdX(2,1)
      dRdX(2,8)=-dRdX(2,2)
      dRdX(2,9)=-dRdX(2,3)

      dRdX(3,4)=(x(4)-x(7))/r(3)
      dRdX(3,5)=(x(5)-x(8))/r(3)
      dRdX(3,6)=(x(6)-x(9))/r(3)
      dRdX(3,7)=-dRdX(3,4)
      dRdX(3,8)=-dRdX(3,5)
      dRdX(3,9)=-dRdX(3,6)

      endsubroutine

  subroutine myprepot(path)
  !-----------------------------------------------------------------------------
  !..
  !.. The O(3P) + HBr --> OH + Br MRCI+Q/CBS(aug-cc-pVnZ(-PP); n = Q,5)+SO 3A"
  !.. potential surface of A. G. S. de Oliveira-Filho, F. R. Ornellas
  !.. and K. A. Peterson
  !..
  !.. Reference:
  !.. Antonio G. S. de Oliveira-Filho, Fernando R. Ornellas and Kirk A. Peterson
  !.. J. Chem. Phys. 136, 174316 (2012). 
  !..
  !.. The Reproducing Kernel Hilbert Space method of Ho and Rabitz is used
  !.. to interpolate the surface consisting of 1110 geometries spanning
  !.. O-H-Br angles of 60-180 deg.  This surface does not contain the
  !.. H + BrO arrangement and no claims of accuracy are made for O-H-Br
  !.. angles smaller than 60 deg. 
  !-----------------------------------------------------------------------------
  !  USAGE:
  !
  !  On Input:
  !      r1 = R(O-H)
  !      r2 = R(H-Br)
  !      r3 = R(Br-O).
  !
  ! On Output:
  !      value = Energy, relative to the asymptotic reactants 
  !              valley (O + HBr(r_e)), in hartree 
  !      dr1 = Derivative of the potential with respect to r1 (R(O-H))
  !            in hartree/bohr 
  !      dr2 = Derivative of the potential with respect to r2 (R(H-Br))
  !            in hartree/bohr 
  !      dr3 = Derivative of the potential with respect to r3 (R(O-Br))
  !            in hartree/bohr 
  !
  !  NOTE:
  !    Before any actual potential energy calculations are made, a single
  !    call to myprepot must be made:
  !      call myprepot
  !
  !    Later, the potential energy is computed by calling mypot:
  !      call mypot(r1,r2,r3,value,dr1,dr2,dr3)
  !
  !   The parameters are read from the o_hbr_a_pp.par file 
  !
  !-----------------------------------------------------------------------------
  !
  !                              Saddle point
  !
  !        r1                          r2                        r3
  !   2.6289099999999999        2.8486750000000001        5.1064215904894388
  !        energy
  !   7.9866812204093529E-003
  !        dr1                         dr2                       dr3
  !   1.5320730578638475E-007   1.5109533511581397E-006  -1.5577335148869720E-007
  !
  !                             Reactants side vdW well
  !
  !        r1                          r2                        r3
  !   4.3746320000000001        2.6855400000000000        7.0601719999996142     
  !        energy
  !  -2.5958471320488563E-003
  !        dr1                         dr2                       dr3
  !   1.4650519388234480E-003   1.4648796113941348E-003  -1.4650322069311608E-003
  !
  !                             Products  side vdW well
  !
  !        r1                          r2                        r3
  !   1.8372599999999999        4.8222319999999996        4.4678939513580840     
  !        energy
  !  -3.3549004949264938E-002
  !        dr1                         dr2                       dr3
  !   6.5041529057150577E-007  -4.8876198921465885E-010   6.6615884752874166E-009
  !
  ! ----------------------------------------------------------------------------
  implicit none
  character(len=1024), intent(in) :: path
  character(len=1024) :: file_path1
  integer, parameter:: dp=kind(0.d0)                   ! double precision
  real(dp), dimension(14) :: phi_oh, phi_hbr 
  real(dp), dimension(1110) :: x_1,x_2,y_1,amn
  real(dp) :: r1,r2,r3,value,dr1,dr2,dr3
  real(dp) :: r_oh,r_hbr,r_bro
  real(dp) :: sum1, sum2, dzdr, dfdr 
  real(dp) :: dbetadz, dedr_oh, dedr_hbr
  real(dp) :: v_oh, v_hbr, v_bro, v3
  real(dp) :: th_y,theta
  real(dp) :: z,re,de,phi_inf,alpha_s,r_s,z_i,phi_z,f_s
  real(dp) :: x_1_max,x_2_max,y_1_max
  real(dp) :: x_1_min,x_2_min,y_1_min
  real(dp) :: q_1,q_2,q_3,dq_1,dq_2,dq_3
  integer :: i, ncall
  data ncall/1/
  save ncall
  save phi_oh,phi_hbr,x_1,x_2,y_1,amn
  
  file_path1 = trim(path)//"/HOBr/o_hbr_a_pp.par"

  if(ncall.eq.1) then 
      open(unit=45,file=file_path1,status="old")
      do i=1,14
          read(45,*) phi_oh(i)
      end do
      do i=1,14
          read(45,*) phi_hbr(i)
      end do
      do i=1,1110
          read(45,*) x_1(i), x_2(i), y_1(i), amn(i)
      end do
      close(45)
      ncall=2
      return
  end if
   
  entry mypot(r1,r2,r3,value,dr1,dr2,dr3)
  
  r_oh=r1
  r_hbr=r2
  r_bro=r3
  
  re=phi_oh(1)
  de=phi_oh(2)
  phi_inf=phi_oh(3)
  alpha_s=phi_oh(4)
  r_s=phi_oh(5)
  z=(r_oh-re)/(r_oh+re)
  phi_z=-phi_inf
  z_i=1.0_dp
  sum2=0.0_dp
  do i=6,14,1
      phi_z=phi_z+phi_oh(i)*z_i
      sum2=sum2+(i-6)*phi_oh(i)*z_i/z
      z_i=z_i*z
  end do
  sum1=phi_z
  f_s=1.0_dp/(exp(alpha_s*(r_oh-r_s))+1.0_dp)
  phi_z=f_s*phi_z+phi_inf
  v_oh=de*((1.0_dp-(re/r_oh)**6*exp(-phi_z*z))**2-1.0_dp)
  
  dzdr=2.0_dp*re/((r_oh+re)**2)
  dfdr=-alpha_s*exp(alpha_s*(r_oh-r_s))*f_s*f_s
  
  dbetadz=phi_z+((1/dzdr)*dfdr*sum1+f_s*sum2)*z
  dedr_oh=2.0_dp*de*(1.0_dp-(re/r_oh)**6*exp(-phi_z*z))*((re/r_oh)**6)*exp(-phi_z*z)*(dzdr*dbetadz+6.0_dp/r_oh)
  
  re=phi_hbr(1)
  de=phi_hbr(2)
  phi_inf=phi_hbr(3)
  alpha_s=phi_hbr(4)
  r_s=phi_hbr(5)
  z=(r_hbr-re)/(r_hbr+re)
  phi_z=-phi_inf
  z_i=1.0_dp
  sum2=0.0_dp
  do i=6,14,1
      phi_z=phi_z+phi_hbr(i)*z_i
      sum2=sum2+(i-6)*phi_hbr(i)*z_i/z
      z_i=z_i*z
  end do
  sum1=phi_z
  f_s=1.0_dp/(exp(alpha_s*(r_hbr-r_s))+1.0_dp)
  phi_z=f_s*phi_z+phi_inf
  v_hbr=de*((1.0_dp-(re/r_hbr)**6*exp(-phi_z*z))**2-1.0_dp)
  
  dzdr=2.0_dp*re/((r_hbr+re)**2)
  dfdr=-alpha_s*exp(alpha_s*(r_hbr-r_s))*f_s*f_s
  
  dbetadz=phi_z+((1/dzdr)*dfdr*sum1+f_s*sum2)*z
  dedr_hbr=2.0_dp*de*(1.0_dp-(re/r_hbr)**6*exp(-phi_z*z))*((re/r_hbr)**6)*exp(-phi_z*z)*(dzdr*dbetadz+6.0_dp/r_hbr)
  
  v_bro=0.0868406_dp*exp(-1.91_dp*(r_bro-3.2509_dp))
  
  th_y=(1.0_dp-((r_oh*r_oh+r_hbr*r_hbr-r_bro*r_bro)/(2.0_dp*r_oh*r_hbr)))/2.0_dp
  
  v3=0.0_dp
  dr1=0.0_dp
  dr2=0.0_dp
  dr3=0.0_dp
  do i=1,1110,1
  
      if (r_oh.ge.x_1(i)) then
          q_1=(1.0_dp/14.0_dp)*(r_oh**(-7))*(1.0_dp-(7.0_dp/9.0_dp)*(x_1(i)/r_oh))
          dq_1=0.5_dp*(r_oh**(-8))*((8.0_dp/9.0_dp)*(x_1(i)/r_oh)-1.0_dp)
      else
          q_1=(1.0_dp/14.0_dp)*(x_1(i)**(-7))*(1.0_dp-(7.0_dp/9.0_dp)*(r_oh/x_1(i)))
          dq_1=-x_1(i)**(-8)/18.0_dp
      end if
  
      if (r_hbr.ge.x_2(i)) then
          q_2=(1.0_dp/14.0_dp)*(r_hbr**(-7))*(1.0_dp-(7.0_dp/9.0_dp)*(x_2(i)/r_hbr))
          dq_2=0.5_dp*(r_hbr**(-8))*((8.0_dp/9.0_dp)*(x_2(i)/r_hbr)-1.0_dp)
      else
          q_2=(1.0_dp/14.0_dp)*(x_2(i)**(-7))*(1.0_dp-(7.0_dp/9.0_dp)*(r_hbr/x_2(i)))
          dq_2=-x_2(i)**(-8)/18.0_dp
      end if
  
      if (th_y.ge.y_1(i)) then
          q_3=1.0_dp+th_y*y_1(i)+2.0_dp*y_1(i)**2*th_y*(1.0_dp-(1.0_dp/3.0_dp)*(y_1(i)/th_y))
          dq_3=y_1(i)+2.0_dp*y_1(i)**2
      else
          q_3=1.0_dp+th_y*y_1(i)+2.0_dp*th_y**2*y_1(i)*(1.0_dp-(1.0_dp/3.0_dp)*(th_y/y_1(i)))
          dq_3=y_1(i)+4.0_dp*y_1(i)*th_y-2.0_dp*th_y**2
      end if
  
      v3=v3+amn(i)*q_1*q_2*q_3
      dr1=dr1+amn(i)*dq_1*q_2*q_3
      dr2=dr2+amn(i)*q_1*dq_2*q_3
      dr3=dr3+amn(i)*q_1*q_2*dq_3
  end do
  value=0.97714_dp*v3+v_oh+v_hbr+v_bro+0.14242006_dp
  
  dr3=0.97714_dp*dr3*r3/(2.0_dp*r1*r2)
  dr1=0.97714_dp*dr1-dr3*((r1-r2*((r_oh**2+r_hbr**2-r_bro**2)/(2*r_oh*r_hbr)))/r3)
  dr2=0.97714_dp*dr2-dr3*((r2-r1*((r_oh**2+r_hbr**2-r_bro**2)/(2*r_oh*r_hbr)))/r3)
  
  dr1=dedr_oh+dr1
  dr2=dedr_hbr+dr2
  dr3=dr3-1.91_dp*0.0868406_dp*exp(-1.91_dp*(r_bro-3.2509_dp))
  end subroutine
