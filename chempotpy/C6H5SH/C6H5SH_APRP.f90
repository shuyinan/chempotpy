!*****************************************************************
!                Linyao Zhang, July, 2019
!*****************************************************************
!   Version:               v1
!   System:                C6H5SH
!   Functional form:       Anchor-Points Reactive Potential (APRP)
!                          (3x3 diabatic fit)
!   Number of derivatives:          1
!   Number of electronic surfaces:  3
!
!   Reference:
!      "Full-dimensional three-state potential energy surfaces and
!       state couplings for photodissociation of thiophenol"
!      by Linyao Zhang, Donald G. Truhlar and Shaozeng Sun, J. Chem. Phys.,
!      submitted July 2019.
!   Notes:
!    1. This PES subroutine has been tested with
!       the following compiler and software:
!       * Intel ifort with MKL
!       * GCC gfortran with LAPACK
!       * Polyrate
!       * Gaussian
!       * ANT
!    2. This routine calculates three coupled singlet PESs for
!       nonadiabatic photodissociation of PhSH to produce
!       PhS + H.
!    3. LAPACK library is needed to diagonize diabatic potential
!       matrix (UU) to yield adiabatic potenitials VV(3)
!    4. Variable 'debug' should be set to .false. for production runs
!    5. Internal working units (for all subroutines except pot):
!       Angstrom, radian, Hartree. Return energy in hartree.
!
!   Internal numbering of atoms:
!        H11      H12
!         \      /
!          C5---C6
!         /      \
!   H10--C4       C1---S7
!         \      /      \
!          C3---C2       H13
!         /      \
!        H9       H8
!
!*****************************************************************

!*****************************************************************
!   *def pot
!      Main PES subroutine for calculating potential energies and
!      gradients of thiophenol. Designed to interface with ANT.
!   Input:
!   igrad          dummy in current potential
!         = 0      Energy only calculation
!         = 1      Energy + Analytic gradient
!   repflag        flag to indicate wheter to use diabatic or
!                  adiabatic representation (used by ANT program)
!   xx(1,i)        X coordinate of atom i
!   xx(2,i)        Y coordinate of atom i
!   xx(3,i)        Z coordinate of atom i
!
!   Output:
!   uu(j,j)        diabatic potential of diabatic state j
!   uu(j,k)        diabatic coupling between state j and k
!   guu(1,i,j,j)   X component of the gradients of diabatical
!                  potential of state j at atom i
!   guu(2,i,j,j)   Y component of the gradients of diabatical
!                  potential of state j at atom i
!   guu(3,i,j,j)   Z component of the gradients of diabatical
!                  potential of state j at atom i
!   guu(1,i,j,k)   X component of the gradients of diabatical
!                  coupling between state j and k at atom i
!   guu(2,i,j,k)   Y component of the gradients of diabatical
!                  coupling between state j and k at atom i
!   guu(3,i,j,k)   Z component of the gradients of diabatical
!                  coupling between state j and k at atom i
!
!   vv(j)          adiabatic potential of adiabatic state j
!   gvv(1,i,j)     X component of the gradients of adiabatic
!                  potential of state j at atom i
!   gvv(2,i,j)     Y component of the gradients of adiabatic
!                  potential of state j at atom i
!   gvv(3,i,j)     Z component of the gradients of adiabatic
!                  potential of state j at atom i
!   dvec(1,i,j,k)  X component of nonadiabatic coupling between
!                  state j and k at atom i
!   dvec(2,i,j,k)  Y component of nonadiabatic coupling between
!                  state j and k at atom i
!   dvec(3,i,j,k)  Z component of nonadiabatic coupling between
!                  state j and k at atom i
!   cc(3,3)        3*3 orghtonal matrix diagonalizing UU matrix
!                  (CC^T*UU*CC yield diagonal matrix with adiabatic
!                  energies as diagonal elements)
!*****************************************************************

      subroutine pes(x,igrad,p,g,d)

      implicit none
      ! number of electronic state
      integer, parameter :: nstates=3
      integer, parameter :: natoms=13
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
      !initialize 
      p=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
        do idir=1,3
          xx(idir,iatom)=x(iatom,idir)/0.529177211
        enddo
      enddo

      call pot(igrad,xx,uu,guu,vv,gvv,dvec,cc,0)
      !output is in hartree/bohr

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

  subroutine pot(igrad,xx,uu,guu,vv,gvv,dvec,cc,repflag)
  !### energis in Hatree
  !### gradients in Hatree/bohr
  !### zero of energy at V1 equilibrium of XQDPT energy
  implicit none
  !--- {'global' constants
  integer, parameter :: natom=13, nstate=3, ntermmax=200 &
       , ntermtype=4, napmax=10, nqtpset=2, nlcmax=6
  double precision, parameter :: pi=3.14159265358979d0
  double precision, parameter :: eV_hartree=3.674931d-2
  double precision, parameter :: ang_bohr=0.529177249d0
  integer, parameter :: n_ps=2 !number of primary & secondary coordinates
  !--- }'global' constants

  !--- {arguments
  integer :: igrad, repflag
  double precision :: xx(3,natom), uu(nstate,nstate), vv(nstate) &
       , guu(3,natom,nstate,nstate), gvv(3,natom,nstate) &
       , dvec(3,natom,nstate,nstate), cc(nstate,nstate)
  !--- }arguments

  integer :: i, j, k, l, ist, jst
  double precision :: x(3,natom), vt(nstate,nstate) &
       , gvt(3,natom,nstate,nstate), tmpmat33(nstate,nstate)

  !--- internals for prim.&sec. pot, tert. pot., tert. coupl.
  !--- qtp has two sets, for bonded & dissoc. S-H
  double precision :: qps(n_ps),qtp(ntermmax,nqtpset), qtc(ntermmax) &
       , bmatqtc(ntermmax,natom*3) !--- attributes

  !--- number/type of terms and atom index lists for qtp
  integer :: nqtp(nqtpset), itermqtptype(ntermmax,nqtpset) &
       , iqtplist(4,ntermmax,nqtpset) &
       , nqtc, qtcsym(ntermmax), qtctyp(ntermmax)

  !--- for diagonalization of prim+sec U
  double precision :: DTAmat(nstate,nstate), Vtmp(nstate)
  !--- zero of energy
  double precision, parameter :: Ezero=-629.4436077422d0
  !--- debug
  logical :: debug_notert, debug

  !--- debug: using tertiary or not
  debug_notert = .false.
  !--- debug: input cartesians in Angstroms; otherwise in bohr
  debug = .false.

  !--- convert Cartesians from bohr to Angstrom
  if(.not.debug) then
    x(:,:) = xx(:,:)*ang_bohr
  else
    x(:,:) = xx(:,:)
  end if

  !@01 Calculate internal coordinates
  !--- primary and secondary internals
  call calcqps(qps, x)
  !--- redundant tertiary internals for tert. diab. pot.
  call calcqtp(qtp, nqtp, itermqtptype, iqtplist, x)
  !--- nonredundant tertiary internals for diab. coupl.
  call calcqtc(qtc, nqtc, qtcsym, qtctyp, bmatqtc, x)

  uu = 0d0
  guu = 0d0
  !@02 calculate prim.&sec. U elements and their grad. w.r.t. x
  call  calcu1_ps(uu, guu, x, qps)
  call  calcu2_ps(uu, guu, x, qps)
  call  calcu3_ps(uu, guu, x, qps)
  call calcu12_ps(uu, guu, x, qps)
  call calcu13_ps(uu, guu, x, qps)
  call calcu23_ps(uu, guu, x, qps)

  !@03 tertiary elements and their grad. w.r.t. x
  if(.not.debug_notert) then
    !--- calculate tertiary potential and its grad. w.r.t. x
    call calcuii_t(vt, gvt, x, qps, qtp, nqtp, itermqtptype, iqtplist)
    !=== {assuming tert. pot. are diabatic
    uu = uu+vt
    guu = guu+gvt
    !=== }assuming tert. pot. are diabatic

    !--- calculate tert. diab. coupl.
    call calcuij_t(uu, guu, x, qps, qtc, nqtc, qtcsym, qtctyp, bmatqtc)

  end if

  !--- for dynamics the zero of energy is
  !--- setted as that of V1 at equilibrium by XQDPT
  if(.not.debug) then
    uu = uu
  else
    do j = 1, 3
      uu(j,j) = uu(j,j) + Ezero
    end do
  end if

  !--- cc: eigenvectors
  call diagonalize(nstate, uu, vv, cc)
  !call diagonalize(vv, cc, uu, nstate)


  !--- get the nonadiabatic coupling elements
  dvec = 0d0
  do i = 1, 3
    do j = 1, natom
      !--- C^T*dU*C
      tmpmat33(:,:) = guu(i,j,:,:)
      tmpmat33 = matmul(transpose(cc), matmul(tmpmat33, cc))
      !--- dV_i = (C^T*dU*C)_ii
      do ist = 1, nstate
        gvv(i,j,ist) = tmpmat33(ist,ist)
      end do
      !--- F_ij = (C^T*dU*C)_ij / (V_j-V_i)
      do ist = 1, nstate-1
        do jst = ist+1, nstate
           dvec(i,j,ist,jst) = tmpmat33(ist,jst)/(vv(jst)-vv(ist))
           dvec(i,j,jst,ist) = -dvec(i,j,ist,jst)
        end do
      end do
    end do
  end do

  !--- convert gradients from Ang-1 to bohr-1
  if(.not.debug) then
    guu = guu*ang_bohr
    gvv = gvv*ang_bohr
    dvec = dvec*ang_bohr
  end if

  contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      Calculate primary & secondary diabatic potential
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!=================================================================
! *def calcqps
!    Calculate primary & secondary internal coordinates.
!  Input:
!    x: Cartesian coordinates in Angstroms
!  Output:
!    qps: values of primary & secondary internal coordinates
!=================================================================
  subroutine calcqps(qps, x)
  implicit none
  double precision :: qps(:), x(3,natom)

  qps(1) = evalBL(x, 7, 13)
  !--- the CSCC torsion is not a good prim. coordinate;
  !--- instead, use the "special" coordinate phi
  !@@@ using new definition for primary torsion
  qps(2) = 0.5*(abs(evalSP1(x, 2, 1, 6, 7, 13))+ &
                abs(evalSP1(x, 3, 4, 5, 7, 13)))

  end subroutine calcqps

!=================================================================
! *def calu1_ps
!    Calculate prim.+sec. U1 and its grad. w.r.t. Cartesians.
!  Input:
!    uu: primary+secondary U matrix in hartree
!    guu: gradient array (in hartree/Angstrom)
!    x: Cartesian coordinates in Angstrom
!    qps: values of primary & secondary internal coordinates
!  Output:
!    uu: primary+secondary U w/ calculated U1_PS added
!    guu: grad. array w/ calculated dU1/dx added to guu(:,:,1,1)
!=================================================================
  subroutine calcu1_ps(uu, guu, x, qps)
  implicit none
  double precision :: uu(nstate,nstate), guu(3,natom,nstate,nstate) &
                      , x(3,natom), qps(:)
  double precision :: R, phi, gradu1(2)
  double precision :: B(1:6), D(0:6), Ra(1:6)
  double precision :: bmatrow(3*natom), bmat(n_ps,3*natom)
  integer :: n, iatomlist(10)

   B(1:6) = (/2.526, 0.9756, 0.0, 3.061, 0.0, 8.456/)
   D(0:6) = (/3.797, 0.7851, 2.581, 0.0, 1.034, 0.0, 0.2404/)
  Ra(1:6) = (/1.651, 1.132, 0.0, 1.022, 0.0, 1.11/)

  R = qps(1); phi = qps(2)

  !--- U11
  uu(1,1) = D(0)-D(1)*(1+B(1)*(R-Ra(1)))*exp(-B(1)*(R-Ra(1)))
  do n = 1, 3
    uu(1,1) = uu(1,1)-D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2)*cos(2*n*phi)
  end do

  !--- (d U11)/(d R)
  gradu1(1) = 0d0
  gradu1(1) = B(1)**2*D(1)*exp(-B(1)*(R-Ra(1)))*(R-Ra(1))
  do n = 1, 3
    gradu1(1) = gradu1(1)+2*B(2*n)*D(2*n)*(R-Ra(2*n)) &
                *exp(-B(2*n)*(R-Ra(2*n))**2)*cos(2*n*phi)
  end do

  !--- (d U11)/(d phi)
  gradu1(2) = 0d0
  do n = 1, 3
    gradu1(2) = gradu1(2)+2*n*D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2)*sin(2*n*phi)
  end do

  !convert eV to Hartree
  uu(1,1) = uu(1,1)*eV_hartree
  gradu1(1) = gradu1(1)*eV_hartree
  gradu1(2) = gradu1(2)*eV_hartree

  !--- use B matrix to convert dU1/dq to dU1/dx
  iatomlist(1:4) = (/7,13,0,0/)
  call calcbmat(bmatrow, x, iatomlist, 1)
  bmat(1,:) = bmatrow(:)
  !--- special coordinate in place of CCSH torsion
  !>>> modify case(6) of calcbmat
  iatomlist(1:5) = (/2,1,6,7,13/)
  call calcbmat(bmatrow, x, iatomlist, 6)
  bmat(2,:) = bmatrow(:)

  guu(:,:,1,1) = guu(:,:,1,1)+reshape(matmul(gradu1,bmat),(/3,natom/))
  end subroutine calcu1_ps

!=================================================================
! *def calu2_ps
!    Calculate prim.+sec. U2 and its grad. w.r.t. Cartesians.
!  Input:
!    uu: primary+secondary U matrix in hartree
!    guu: gradient array (in hartree/Angstrom)
!    x: Cartesian coordinates in Angstrom
!    qps: values of primary & secondary internal coordinates
!  Output:
!    uu: primary+secondary U w/ calculated U2_PS added
!    guu: grad. array w/ calculated dU2/dx added to guu(:,:,2,2)
!=================================================================
  subroutine calcu2_ps(uu, guu, x, qps)
  implicit none
  double precision :: uu(nstate,nstate), guu(3,natom,nstate,nstate) &
                      , x(3,natom), qps(:)
  double precision :: R, phi, gradu2(2)
  double precision :: B(1:4), D(0:4), Ra(1:4)
  double precision :: bmatrow(3*natom), bmat(n_ps,3*natom)
  integer :: n, iatomlist(10)

   B(1:4) = (/3.034, 8.449, 0.0, 6.72/)
   D(0:4) = (/4.724, 1.628, 0.1262, 0.0, -0.09122/)
  Ra(1:4) = (/1.309, 1.307, 0.0, 2.148/)

  R = qps(1); phi = qps(2)

  !--- U22
  uu(2,2) = D(0)+D(1)*(1-exp(-B(1)*(R-Ra(1))))**2
  do n = 1, 2
    uu(2,2) = uu(2,2)-D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2)*cos(2*n*phi)
  end do

  !--- (d U22)/(d R)
  gradu2(1) = 0d0
  gradu2(1) = 2*B(1)*D(1)*exp(-B(1)*(R-Ra(1)))*(1-exp(-B(1)*(R-Ra(1))))
  do n = 1, 2
    gradu2(1) = gradu2(1)+2*B(2*n)*D(2*n)*(R-Ra(2*n)) &
                *exp(-B(2*n)*(R-Ra(2*n))**2)*cos(2*n*phi)
  end do

  !--- (d U22)/(d phi)
  gradu2(2) = 0d0
  do n = 1, 2
    gradu2(2) = gradu2(2)+2*n*D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2)*sin(2*n*phi)
  end do

  !convert eV to Hartree
  uu(2,2) = uu(2,2)*eV_hartree
  gradu2(1) = gradu2(1)*eV_hartree
  gradu2(2) = gradu2(2)*eV_hartree

  !--- use B matrix to convert dU2/dq to dU2/dx
  iatomlist(1:4) = (/7,13,0,0/)
  call calcbmat(bmatrow, x, iatomlist, 1)
  bmat(1,:) = bmatrow(:)
  !--- special coordinate in place of CCSH torsion
  !>>> modify case(6) of calcbmat
  iatomlist(1:5) = (/2,1,6,7,13/)
  call calcbmat(bmatrow, x, iatomlist, 6)
  bmat(2,:) = bmatrow(:)

  guu(:,:,2,2) = guu(:,:,2,2)+reshape(matmul(gradu2, bmat), (/3,natom/))

  end subroutine calcu2_ps

!=================================================================
! *def calu3_ps
!    Calculate prim.+sec. U3 and its grad. w.r.t. Cartesians.
!  Input:
!    uu: primary+secondary U matrix in hartree
!    guu: gradient array (in hartree/Angstrom)
!    x: Cartesian coordinates in Angstrom
!    qps: values of primary & secondary internal coordinates
!  Output:
!    uu: primary+secondary U w/ calculated U3_PS added
!    guu: grad. array w/ calculated dU3/dx added to guu(:,:,3,3)
!=================================================================
  subroutine calcu3_ps(uu, guu, x, qps)
  implicit none
  double precision :: uu(nstate,nstate), guu(3,natom,nstate,nstate) &
                      , x(3,natom), qps(:)
  double precision :: R, phi, gradu3(2)
  double precision :: B(1:6), D(0:6), Ra(1:6)
  double precision :: bmatrow(3*natom), bmat(n_ps,3*natom)
  integer :: n, iatomlist(10)

   B(1:6) = (/3.087, 0.9997, 0.0, 2.641, 0.0, 4.85/)
   D(0:6) = (/3.537, 1.713, -2.544, 0.0, -0.9263, 0.0, -0.3608/)
  Ra(1:6) = (/1.366, 1.066, 0.0, 1.045, 0.0, 0.913/)

  R = qps(1); phi = qps(2)

  !--- U33
  uu(3,3) = D(0)-D(1)*(1+B(1)*(R-Ra(1)))*exp(-B(1)*(R-Ra(1)))
  do n = 1, 3
    uu(3,3) = uu(3,3)-D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2)*cos(2*n*phi)
  end do

  !--- (d U33)/(d R)
  gradu3(1) = 0d0
  gradu3(1) = B(1)**2*D(1)*exp(-B(1)*(R-Ra(1)))*(R-Ra(1))
  do n = 1, 3
    gradu3(1) = gradu3(1)+2*B(2*n)*D(2*n)*(R-Ra(2*n)) &
                *exp(-B(2*n)*(R-Ra(2*n))**2)*cos(2*n*phi)
  end do

  !--- (d U33)/(d phi)
  gradu3(2) = 0d0
  do n = 1, 3
    gradu3(2) = gradu3(2)+2*n*D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2)*sin(2*n*phi)
  end do

  !convert eV to Hartree
  uu(3,3) = uu(3,3)*eV_hartree
  gradu3(1) = gradu3(1)*eV_hartree
  gradu3(2) = gradu3(2)*eV_hartree

  !--- use B matrix to convert dU3/dq to dU3/dx
  iatomlist(1:4) = (/7,13,0,0/)
  call calcbmat(bmatrow, x, iatomlist, 1)
  bmat(1,:) = bmatrow(:)
  !--- special coordinate in place of CCSH torsion
  !>>> modify case(6) of calcbmat
  iatomlist(1:5) = (/2,1,6,7,13/)
  call calcbmat(bmatrow, x, iatomlist, 6)
  bmat(2,:) = bmatrow(:)

  guu(:,:,3,3) = guu(:,:,3,3)+reshape(matmul(gradu3, bmat), (/3,natom/))

  end subroutine calcu3_ps

!=================================================================
! *def calu12_ps
!    Calculate prim.+sec. U12 and its grad. w.r.t. Cartesians.
!  Input:
!    uu: primary+secondary U matrix in hartree
!    guu: gradient array (in hartree/Angstrom)
!    x: Cartesian coordinates in Angstrom
!    qps: values of primary & secondary internal coordinates
!  Output:
!    uu: primary+secondary U w/ calculated U12_PS added
!    guu: grad. array w/ calculated dU12/dx added to guu(:,:,1,2)
!=================================================================
  subroutine calcu12_ps(uu, guu, x, qps)
  implicit none
  double precision :: uu(nstate,nstate), guu(3,natom,nstate,nstate) &
                      , x(3,natom), qps(:)
  double precision :: R, phi, gradu12(2)
  double precision :: B(1:2), D(1:2), Ra(1:2)
  double precision :: bmatrow(3*natom), bmat(n_ps,3*natom)
  integer :: n, iatomlist(10)

   B(1:2) = (/1.825, 12.01/)
   D(1:2) = (/0.4797, -0.3624/)
  Ra(1:2) = (/2.032, 1.854/)

  R = qps(1); phi = qps(2)

  !--- U12
  do n = 1, 2
    uu(1,2) = uu(1,2)+D(n)*exp(-B(n)*(R-Ra(n))**2)*(sin(2*phi))**2
  end do

  !--- (d U12)/(d R)
  gradu12(1) = 0d0
  do n = 1, 2
    gradu12(1) = gradu12(1)-2*B(n)*D(n)*(R-Ra(n)) &
                 *exp(-B(n)*(R-Ra(n))**2)*(sin(2*phi))**2
  end do

  !--- (d U12)/(d phi)
  gradu12(2) = 0d0
  do n = 1, 2
    gradu12(2) = gradu12(2)+4*D(n)*exp(-B(n)*(R-Ra(n))**2) &
                 *sin(2*phi)*cos(2*phi)
  end do

  !convert eV to Hartree
  uu(1,2) = uu(1,2)*eV_hartree
  uu(2,1) = uu(1,2)
  gradu12(1) = gradu12(1)*eV_hartree
  gradu12(2) = gradu12(2)*eV_hartree

  !--- use B matrix to convert dU12/dq to dU12/dx
  iatomlist(1:4) = (/7,13,0,0/)
  call calcbmat(bmatrow, x, iatomlist, 1)
  bmat(1,:) = bmatrow(:)
  !--- special coordinate in place of CCSH torsion
  !>>> modify case(6) of calcbmat
  iatomlist(1:5) = (/2,1,6,7,13/)
  call calcbmat(bmatrow, x, iatomlist, 6)
  bmat(2,:) = bmatrow(:)

  guu(:,:,1,2) = guu(:,:,1,2)+reshape(matmul(gradu12, bmat), (/3,natom/))
  guu(:,:,2,1) = guu(:,:,1,2)

  end subroutine calcu12_ps

!=================================================================
! *def calu13_ps
!    Calculate prim.+sec. U13 and its grad. w.r.t. Cartesians.
!  Input:
!    uu: primary+secondary U matrix in hartree
!    guu: gradient array (in hatree/Angstrom)
!    x: Cartesian coordinates in Angstrom
!    qps: values of primary & secondary internal coordinates
!  Output:
!    uu: primary+secondary U w/ calculated U13_PS added
!    guu: grad. array w/ calculated dU13/dx added to guu(:,:,1,3)
!=================================================================
  subroutine calcu13_ps(uu, guu, x, qps)
  implicit none
  double precision :: uu(nstate,nstate), guu(3,natom,nstate,nstate) &
                      , x(3,natom), qps(:)
  double precision :: R, phi, gradu13(2)
  double precision :: B(2:6), D(2:6), Ra(2:6)
  double precision :: bmatrow(3*natom), bmat(n_ps,3*natom)
  integer :: n, iatomlist(10)

   B(2:6) = (/0.8257, 0.0, 2.01, 0.0, 6.049/)
   D(2:6) = (/2.554, 0.0, 1.126, 0.0, 0.2316 /)
  Ra(2:6) = (/0.615, 0.0, 0.7497, 0.0, 0.9754/)

  R = qps(1); phi = qps(2)

  !--- U13
  do n = 1, 3
    uu(1,3) = uu(1,3)+D(2*n)*R*exp(-B(2*n)*(R-Ra(2*n))**2)*sin(2*n*phi)
  end do

  !--- (d U13)/(d R)
  gradu13(1) = 0d0
  do n = 1, 3
    gradu13(1) = gradu13(1)+(D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2)-2*B(2*n) &
                 *D(2*n)*R*(R-Ra(2*n))*exp(-B(2*n)*(R-Ra(2*n))**2)) &
                 *sin(2*n*phi)
  end do

  !--- (d U13)/(d phi)
  gradu13(2) = 0d0
  do n = 1, 3
    gradu13(2) = gradu13(2)+2*n*D(2*n)*R*exp(-B(2*n)*(R-Ra(2*n))**2) &
                 *cos(2*n*phi)
  end do

  !convert eV to Hartree
  uu(1,3) = uu(1,3)*eV_hartree
  uu(3,1) = uu(1,3)
  gradu13(1) = gradu13(1)*eV_hartree
  gradu13(2) = gradu13(2)*eV_hartree

  !--- use B matrix to convert dU13/dq to dU13/dx
  iatomlist(1:4) = (/7,13,0,0/)
  call calcbmat(bmatrow, x, iatomlist, 1)
  bmat(1,:) = bmatrow(:)
  !--- special coordinate in place of CCSH torsion
  !>>> modify case(6) of calcbmat
  iatomlist(1:5) = (/2,1,6,7,13/)
  call calcbmat(bmatrow, x, iatomlist, 6)
  bmat(2,:) = bmatrow(:)

  guu(:,:,1,3) = guu(:,:,1,3)+reshape(matmul(gradu13, bmat), (/3,natom/))
  guu(:,:,3,1) = guu(:,:,1,3)

  end subroutine calcu13_ps

!=================================================================
! *def calu23_ps
!    Calculate prim.+sec. U23 and its grad. w.r.t. Cartesians.
!  Input:
!    uu: primary+secondary U matrix in hartree
!    guu: gradient array (in hartree/Angstrom)
!    x: Cartesian coordinates in Angstrom
!    qps: values of primary & secondary internal coordinates
!  Output:
!    uu: primary+secondary U w/ calculated U23_PS added
!    guu: grad. array w/ calculated dU23/dx added to guu(:,:,2,3)
!=================================================================
  subroutine calcu23_ps(uu, guu, x, qps)
  implicit none
  double precision :: uu(nstate,nstate), guu(3,natom,nstate,nstate) &
                      , x(3,natom), qps(:)
  double precision :: R, phi, gradu23(2)
  double precision :: B(0:4), D(0:4), Ra(0:4)
  double precision :: bmatrow(3*natom), bmat(n_ps,3*natom)
  integer :: n, iatomlist(10)

   B(0:4) = (/2.575, 0.0, 4.073, 0.0, 19.18/)
   D(0:4) = (/0.05761, 0.0, -0.3011, 0.0, -0.1517/)
  Ra(0:4) = (/1.44, 0.0, 2.784, 0.0, 2.172/)

  R = qps(1); phi = qps(2)

  !--- U23
  do n = 0, 2
    uu(2,3) = uu(2,3)+D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2)*sin(2**n*phi)
  end do

  !--- 2018.08.27-modify an error in derivative
  !--- (d U23)/(d R)
  gradu23(1) = 0d0
  do n = 0, 2
    gradu23(1) = gradu23(1)-2*B(2*n)*D(2*n)*(R-Ra(2*n)) &
                 *exp(-B(2*n)*(R-Ra(2*n))**2) &
                 *sin(2**n*phi)
  end do

  !--- 2018.08.27-modify an error in derivative
  !--- (d U23)/(d phi)
  gradu23(2) = 0d0
  do n = 0, 2
    gradu23(2) = gradu23(2)+2**n*D(2*n)*exp(-B(2*n)*(R-Ra(2*n))**2) &
                 *cos(2**n*phi)
  end do

  !convert eV to Hartree
  uu(2,3) = uu(2,3)*eV_hartree
  uu(3,2) = uu(2,3)
  gradu23(1) = gradu23(1)*eV_hartree
  gradu23(2) = gradu23(2)*eV_hartree

  !--- use B matrix to convert dU23/dq to dU23/dx
  iatomlist(1:4) = (/7,13,0,0/)
  call calcbmat(bmatrow, x, iatomlist, 1)
  bmat(1,:) = bmatrow(:)
  !--- special coordinate in place of CCSH torsion
  !>>> modify case(6) of calcbmat
  iatomlist(1:5) = (/2,1,6,7,13/)
  call calcbmat(bmatrow, x, iatomlist, 6)
  bmat(2,:) = bmatrow(:)

  guu(:,:,2,3) = guu(:,:,2,3)+reshape(matmul(gradu23, bmat), (/3,natom/))
  guu(:,:,3,2) = guu(:,:,2,3)

  end subroutine calcu23_ps

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!             Calculate tertiary diabatic potential
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!=================================================================
! *def calcuii_t
!    Calculate tert. diab. pot. and its grad. w.r.t. Cartesians.
!  Input:
!    x: Cartesian coordinates in Angstroms
!    qps: values of primary & secondary internal coordinates
!    qtp: values of tert. redund. internals
!    nterm: number of terms in each qtp set
!    itypeqtp: type of each term (and coord.)
!    iqtplist: atom index list array
!  Output:
!    vt: tertiary Uii in hartree
!    gvt: grad. array w/ calculated dU/dx
!=================================================================
  subroutine calcuii_t(vt, gvt, x, qps, qtp, nterm, itypeqtp, iqtplist)
  implicit none
  double precision :: vt(nstate,nstate), gvt(3,natom,nstate,nstate) &
      , x(3,natom), qps(:), qtp(:,:)
  integer :: nterm(nqtpset), itypeqtp(ntermmax,nqtpset) &
      , iqtplist(4,ntermmax,nqtpset)
  integer :: i, j, imin
  double precision :: vmin, vtmp, gvtmp(3,natom)

  vt = 0d0; gvt = 0d0
  !--- calculate vt and gvt in diabatic order
  !--- uii_t is the actual routine
  do i = 1, nstate
     call uii_t(vt, gvt, x, qps, qtp, nterm, itypeqtp, &
          iqtplist, i)
  end do

  end subroutine calcuii_t

  !--- Actual core routine for calcuii_t
  !--- Parameters have the same meaning as in calcuii_t
  subroutine uii_t(vt, gvt, x, qps, qtp, nterm, itypeqtp, &
       iqtplist, istate)
    implicit none
    double precision :: vt(nstate,nstate) &
         , gvt(3,natom,nstate,nstate), x(3,natom), qps(:), qtp(:,:)
    integer :: nterm(nqtpset), itypeqtp(ntermmax,nqtpset) &
         , iqtplist(4,ntermmax,nqtpset)
    integer :: istate
    double precision :: k(ntermmax,napmax) &
         , q0(ntermmax,napmax) &
         , R0(napmax), phi0(napmax), theta0(napmax) &
         , eterm, gterm, energy(0:napmax)
    integer :: iap, iterm, nap, n(ntermmax,napmax), ntermap &
         , iset, iatomlist(10), iqtpset(napmax)
    double precision :: bmatrow(3*natom), gradc(3,natom,0:napmax)
    double precision,allocatable :: bmat(:,:), gg(:)
    double precision :: tent(napmax), dtent(0:napmax), R

    call assignparamuii_t(k, q0, n, istate)
    !--- set up state-dependent variables
    !--- Here 'istate' is diabatic state at phi = 0
    energy = 0d0
    select case(istate)
    case(1)
       !--- # of anchor points
       nap = 7
       !--- anchor points 1-6 use internals set 1
       !--- anchor point  7   uses          set 2
       iqtpset(1:6) = 1
       iqtpset(7) = 2
       !--- prim&sec coord. of each anchor point
       R0(1:7) = (/1.2d0,1.35d0,1.6d0,2.1d0,2.5d0,3.4d0,5.0d0/)
       phi0(1:7) = 0d0
       !--- const. part of tert. potential
       !--- which is relaxed minus unrelaxed energy at APs
       energy(1:7) = (/&
            -0.0024405948,&
            -0.0019666022,&
            -0.0020489942,&
            -0.0023584233,&
            -0.0017178094,&
            -0.0008860257,&
            -0.0007388447/)
    case(2)
       !--- # of anchor points
       nap = 3
       !--- anchor points 1-3 use internals set 1
       iqtpset(1:3) = 1
       !--- prim&sec coord. of each anchor point
       R0(1:3) = (/1.25d0,1.35d0,1.6d0/)
       phi0(1:3) = 0d0
       !--- const. part of tert. potential
       !---  which is relaxed minus unrelaxed energy at APs
       energy(1:3) = (/&
            -0.0064003322,&
            -0.0064962846,&
            -0.0066046951/)
    case(3)
       !--- # of anchor points
       nap = 7
       !--- anchor points 1-3 uses internals set 1
       !--- anchor point 4    uses           set 2
       iqtpset(1:6) = 1
       iqtpset(7) = 2
       !--- prim&sec coord. of each anchor point
       R0(1:7) = (/1.25d0,1.35d0,1.6d0,2.1d0,2.5d0,3.4d0,5.0d0/)
       phi0(1:7) = 0d0
       !--- const. part of tert. potential
       !---  which is relaxed minus unrelaxed energy at APs
       energy(1:7) = (/&
            -0.0086052169,&
            -0.0067687440,&
            -0.0054494440,&
            -0.0042860713,&
            -0.0073215638,&
            -0.0022082657,&
            -0.0019861529/)
    end select
    !energy(:)=0d0

    !--- loop over anchor points
    do iap = 1, nap
       ntermap = nterm(iqtpset(iap))
       iset = iqtpset(iap)
       allocate(bmat(ntermap,3*natom),gg(ntermap))
       do iterm = 1, ntermap
          !--- energy and grad. w.r.t. internals
          call calcFFterm(eterm, gterm, k(iterm,iap), &
               q0(iterm,iap), n(iterm,iap), &
               qtp(iterm,iset), &
               itypeqtp(iterm, iset))
          energy(iap) = energy(iap)+eterm
          gg(iterm) = gterm
          !--- calc. B matrix row
          iatomlist(1:4) = iqtplist(:,iterm,iqtpset(iap))
          call calcbmat(bmatrow, x, iatomlist, &
               itypeqtp(iterm, iset))
          bmat(iterm,:) = bmatrow(:)
       end do
       !--- use B matrix to convert to grad. w.r.t. Cartesians
       gradc(:,:,iap) = reshape(matmul(gg, bmat), (/3,natom/))
       deallocate(bmat,gg)
    end do
    !--- calculate tent func. and dtent/dR
    R = qps(1)
    do iap = 1, nap
       call calctent1(tent(iap), dtent(iap), R, R0, &
            iap-1, iap, iap+1, nap)
    end do

    energy(0) = 0d0; gradc(:,:,0) = 0d0; dtent(0) = 0d0
    do iap = 1, nap
       !--- interpolate anchor point values using tent func.
       energy(0) = energy(0)+energy(iap)*tent(iap)
       gradc(:,:,0) = gradc(:,:,0)+gradc(:,:,iap)*tent(iap)
       !--- interpolate dtent/dR to get dV/dR
       dtent(0) = dtent(0)+dtent(iap)*energy(iap)
    end do
    !--- use B matrix to convert dV/dR to dV/dx
    iatomlist(1:4) = (/7,13,0,0/)
    call calcbmat(bmatrow, x, iatomlist, 1)
    gradc(:,:,0) = gradc(:,:,0)+reshape(dtent(0)*bmatrow, (/3,natom/))

    !--- collect results and convert energy unit to eV
    vt(istate,istate) = energy(0)
    gvt(:,:,istate,istate) = gradc(:,:,0)
  end subroutine uii_t

!=================================================================
! *def calcFFterm
!    Calculate one FF term and its gradient.
!  Input:
!    k: force constant
!    q0: rest value
!    n: multiplicity (for torsion) or dummy (for others)
!    q: internal coordinate
!    itype: type of term
!  Output:
!    ee: energy
!    gg: gradient w.r.t. the internal coordinate
!=================================================================
  subroutine calcFFterm(ee, gg, k, q0, n, q, itype)
  implicit none
  double precision :: ee, gg, k, q0, q
  integer :: n, itype
  logical, parameter :: debug=.false.

  select case (itype)
  !--- bond length
  case(1)
     ee = 0.5d0*k*(1d0-q0/q)**2
     gg = k*(1d0-q0/q)*q0/q**2
  !--- bond angle
  case(2)
     ee = 0.5d0*k*(cos(q)-cos(q0))**2
     gg = -k*(cos(q)-cos(q0))*sin(q)
  !--- torsion
  case(3)
     ee = 0.5d0*k*(1d0-cos(n*(q-q0)))
     gg = 0.5d0*n*k*sin(n*(q-q0))
  !--- oop distance
  case(5)
     ee = 0.5d0*k*(q-q0)**2
     gg = k*(q-q0)
  end select
10   FORMAT(A10, 4f12.6)
  end subroutine calcFFterm

!=================================================================
! *def assignqtp
!    Assign attributes of internal coordinates for tertiary potentials
!  Input:
!    NA
!  Output:
!    nterm: number of terms in each qtp set
!    itypeqtp: type of each term (and coord.)
!    iqtplist: atom index list array
!  ~200 lines of params
!=================================================================
  subroutine assignqtp(nterm, itypeqtp, iqtplist)
  implicit none
  integer :: nterm(nqtpset), itypeqtp(ntermmax,nqtpset) &
       , iqtplist(4,ntermmax,nqtpset)

  nterm(1:2) = (/61,60/)
  !--- type of each FF term of each set
  !--- 1=BL(BondLen), 2=BA(BondAng), 3=TO(Torsion), 5=OD(OOPDist)
  !----- set 1, for short R, corresp. to atom connectivity of Ph-S-H
  itypeqtp( 1:12,1) = 1
  itypeqtp(13:31,1) = 2
  itypeqtp(32:55,1) = 3
  itypeqtp(56:61,1) = 5
  !----- set 2, for long R, corresp. to atom connectivity of Ph-S + H
  itypeqtp( 1:12,2) = 1
  itypeqtp(13:30,2) = 2
  itypeqtp(31:54,2) = 3
  itypeqtp(55:60,2) = 5

  iqtplist = 0
  !--- first set (for short R)
  !--- BEGIN generated by 'python gen_qtp.py v1_ap1.txt'
  iqtplist(1:2,  1,1) = (/ 1, 2/)
  iqtplist(1:2,  2,1) = (/ 1, 6/)
  iqtplist(1:2,  3,1) = (/ 1, 7/)
  iqtplist(1:2,  4,1) = (/ 2, 3/)
  iqtplist(1:2,  5,1) = (/ 2, 8/)
  iqtplist(1:2,  6,1) = (/ 3, 4/)
  iqtplist(1:2,  7,1) = (/ 3, 9/)
  iqtplist(1:2,  8,1) = (/ 4, 5/)
  iqtplist(1:2,  9,1) = (/ 4,10/)
  iqtplist(1:2, 10,1) = (/ 5, 6/)
  iqtplist(1:2, 11,1) = (/ 5,11/)
  iqtplist(1:2, 12,1) = (/ 6,12/)
  iqtplist(1:3, 13,1) = (/ 1, 2, 3/)
  iqtplist(1:3, 14,1) = (/ 1, 2, 8/)
  iqtplist(1:3, 15,1) = (/ 1, 6, 5/)
  iqtplist(1:3, 16,1) = (/ 1, 6,12/)
  iqtplist(1:3, 17,1) = (/ 1, 7,13/)
  iqtplist(1:3, 18,1) = (/ 2, 1, 6/)
  iqtplist(1:3, 19,1) = (/ 2, 1, 7/)
  iqtplist(1:3, 20,1) = (/ 2, 3, 4/)
  iqtplist(1:3, 21,1) = (/ 2, 3, 9/)
  iqtplist(1:3, 22,1) = (/ 3, 2, 8/)
  iqtplist(1:3, 23,1) = (/ 3, 4, 5/)
  iqtplist(1:3, 24,1) = (/ 3, 4,10/)
  iqtplist(1:3, 25,1) = (/ 4, 3, 9/)
  iqtplist(1:3, 26,1) = (/ 4, 5, 6/)
  iqtplist(1:3, 27,1) = (/ 4, 5,11/)
  iqtplist(1:3, 28,1) = (/ 5, 4,10/)
  iqtplist(1:3, 29,1) = (/ 5, 6,12/)
  iqtplist(1:3, 30,1) = (/ 6, 1, 7/)
  iqtplist(1:3, 31,1) = (/ 6, 5,11/)
  iqtplist(1:4, 32,1) = (/ 1, 2, 3, 4/)
  iqtplist(1:4, 33,1) = (/ 1, 2, 3, 9/)
  iqtplist(1:4, 34,1) = (/ 1, 6, 5, 4/)
  iqtplist(1:4, 35,1) = (/ 1, 6, 5,11/)
  iqtplist(1:4, 36,1) = (/ 2, 1, 6, 5/)
  iqtplist(1:4, 37,1) = (/ 2, 1, 6,12/)
  iqtplist(1:4, 38,1) = (/ 2, 3, 4, 5/)
  iqtplist(1:4, 39,1) = (/ 2, 3, 4,10/)
  iqtplist(1:4, 40,1) = (/ 3, 2, 1, 6/)
  iqtplist(1:4, 41,1) = (/ 3, 2, 1, 7/)
  iqtplist(1:4, 42,1) = (/ 3, 4, 5, 6/)
  iqtplist(1:4, 43,1) = (/ 3, 4, 5,11/)
  iqtplist(1:4, 44,1) = (/ 4, 3, 2, 8/)
  iqtplist(1:4, 45,1) = (/ 4, 5, 6,12/)
  iqtplist(1:4, 46,1) = (/ 5, 4, 3, 9/)
  iqtplist(1:4, 47,1) = (/ 5, 6, 1, 7/)
  iqtplist(1:4, 48,1) = (/ 6, 1, 2, 8/)
  iqtplist(1:4, 49,1) = (/ 6, 5, 4,10/)
  iqtplist(1:4, 50,1) = (/11, 5, 4,10/)
  iqtplist(1:4, 51,1) = (/11, 5, 6,12/)
  iqtplist(1:4, 52,1) = (/12, 6, 1, 7/)
  iqtplist(1:4, 53,1) = (/ 8, 2, 1, 7/)
  iqtplist(1:4, 54,1) = (/ 8, 2, 3, 9/)
  iqtplist(1:4, 55,1) = (/ 9, 3, 4,10/)
  iqtplist(1:4, 56,1) = (/ 1, 3, 8, 2/)
  iqtplist(1:4, 57,1) = (/ 1, 5,12, 6/)
  iqtplist(1:4, 58,1) = (/ 2, 4, 9, 3/)
  iqtplist(1:4, 59,1) = (/ 2, 6, 7, 1/)
  iqtplist(1:4, 60,1) = (/ 3, 5,10, 4/)
  iqtplist(1:4, 61,1) = (/ 4, 6,11, 5/)
  !--- END generated by 'python gen_qtp.py v1_ap1.txt'

  !--- second set (for long R)
  !--- BEGIN generated by 'python gen_qtp.py v1_ap4.txt'
  iqtplist(1:2,  1,2) = (/ 1, 2/)
  iqtplist(1:2,  2,2) = (/ 1, 6/)
  iqtplist(1:2,  3,2) = (/ 1, 7/)
  iqtplist(1:2,  4,2) = (/ 2, 3/)
  iqtplist(1:2,  5,2) = (/ 2, 8/)
  iqtplist(1:2,  6,2) = (/ 3, 4/)
  iqtplist(1:2,  7,2) = (/ 3, 9/)
  iqtplist(1:2,  8,2) = (/ 4, 5/)
  iqtplist(1:2,  9,2) = (/ 4,10/)
  iqtplist(1:2, 10,2) = (/ 5, 6/)
  iqtplist(1:2, 11,2) = (/ 5,11/)
  iqtplist(1:2, 12,2) = (/ 6,12/)
  iqtplist(1:3, 13,2) = (/ 1, 2, 3/)
  iqtplist(1:3, 14,2) = (/ 1, 2, 8/)
  iqtplist(1:3, 15,2) = (/ 1, 6, 5/)
  iqtplist(1:3, 16,2) = (/ 1, 6,12/)
  iqtplist(1:3, 17,2) = (/ 2, 1, 6/)
  iqtplist(1:3, 18,2) = (/ 2, 1, 7/)
  iqtplist(1:3, 19,2) = (/ 2, 3, 4/)
  iqtplist(1:3, 20,2) = (/ 2, 3, 9/)
  iqtplist(1:3, 21,2) = (/ 3, 2, 8/)
  iqtplist(1:3, 22,2) = (/ 3, 4, 5/)
  iqtplist(1:3, 23,2) = (/ 3, 4,10/)
  iqtplist(1:3, 24,2) = (/ 4, 3, 9/)
  iqtplist(1:3, 25,2) = (/ 4, 5, 6/)
  iqtplist(1:3, 26,2) = (/ 4, 5,11/)
  iqtplist(1:3, 27,2) = (/ 5, 4,10/)
  iqtplist(1:3, 28,2) = (/ 5, 6,12/)
  iqtplist(1:3, 29,2) = (/ 6, 1, 7/)
  iqtplist(1:3, 30,2) = (/ 6, 5,11/)
  iqtplist(1:4, 31,2) = (/ 1, 2, 3, 4/)
  iqtplist(1:4, 32,2) = (/ 1, 2, 3, 9/)
  iqtplist(1:4, 33,2) = (/ 1, 6, 5, 4/)
  iqtplist(1:4, 34,2) = (/ 1, 6, 5,11/)
  iqtplist(1:4, 35,2) = (/ 2, 1, 6, 5/)
  iqtplist(1:4, 36,2) = (/ 2, 1, 6,12/)
  iqtplist(1:4, 37,2) = (/ 2, 3, 4, 5/)
  iqtplist(1:4, 38,2) = (/ 2, 3, 4,10/)
  iqtplist(1:4, 39,2) = (/ 3, 2, 1, 6/)
  iqtplist(1:4, 40,2) = (/ 3, 2, 1, 7/)
  iqtplist(1:4, 41,2) = (/ 3, 4, 5, 6/)
  iqtplist(1:4, 42,2) = (/ 3, 4, 5,11/)
  iqtplist(1:4, 43,2) = (/ 4, 3, 2, 8/)
  iqtplist(1:4, 44,2) = (/ 4, 5, 6,12/)
  iqtplist(1:4, 45,2) = (/ 5, 4, 3, 9/)
  iqtplist(1:4, 46,2) = (/ 5, 6, 1, 7/)
  iqtplist(1:4, 47,2) = (/ 6, 1, 2, 8/)
  iqtplist(1:4, 48,2) = (/ 6, 5, 4,10/)
  iqtplist(1:4, 49,2) = (/11, 5, 4,10/)
  iqtplist(1:4, 50,2) = (/11, 5, 6,12/)
  iqtplist(1:4, 51,2) = (/12, 6, 1, 7/)
  iqtplist(1:4, 52,2) = (/ 8, 2, 1, 7/)
  iqtplist(1:4, 53,2) = (/ 8, 2, 3, 9/)
  iqtplist(1:4, 54,2) = (/ 9, 3, 4,10/)
  iqtplist(1:4, 55,2) = (/ 1, 3, 8, 2/)
  iqtplist(1:4, 56,2) = (/ 1, 5,12, 6/)
  iqtplist(1:4, 57,2) = (/ 2, 4, 9, 3/)
  iqtplist(1:4, 58,2) = (/ 2, 6, 7, 1/)
  iqtplist(1:4, 59,2) = (/ 3, 5,10, 4/)
  iqtplist(1:4, 60,2) = (/ 4, 6,11, 5/)
  !--- END generated by 'python gen_qtp.py v1_ap4.txt'
  end subroutine assignqtp

!=================================================================
! *def calcqtp
!    Calculate tertiary internal coordinates as defined
!    by QuickFF for diabatic potentials.
!  Input:
!    x: Cartesian coordinates in Angstroms
!  Output:
!    qtp: values of tertiary internal coordinates for potential
!    nterm: number of terms in each qtp set
!    itypeqtp: type of each term (and coord.)
!    iqtplist: atom index list array
!=================================================================
  subroutine calcqtp(qtp, nterm, itypeqtp, iqtplist, x)
  implicit none
  double precision :: qtp(:,:), x(3,natom)
  integer :: nterm(nqtpset), itypeqtp(ntermmax,nqtpset) &
       , iqtplist(4,ntermmax,nqtpset)
  integer :: iterm, iset, ii, jj, kk, ll

  !--- assign attributes of tert. int. for potentials
  call assignqtp(nterm, itypeqtp, iqtplist)

  do iset = 1, nqtpset
     do iterm = 1, nterm(iset)
        ii = iqtplist(1,iterm,iset)
        jj = iqtplist(2,iterm,iset)
        kk = iqtplist(3,iterm,iset)
        ll = iqtplist(4,iterm,iset)
        select case(itypeqtp(iterm,iset))
        !--- bond length
        case(1)
           qtp(iterm,iset) = evalBL(x, ii, jj)
        !--- bond angle
        case(2)
           qtp(iterm,iset) = evalBA(x, ii, jj, kk)
        !--- torsion
        case(3)
           qtp(iterm,iset) = evalTO(x, ii, jj, kk, ll)
        !--- oop distance
        case(5)
           qtp(iterm,iset) = evalOD(x, ii, jj, kk, ll)
        end select
     end do
  end do

  end subroutine calcqtp

!=================================================================
! *def assignparamuii_t
!    Assign tertiary FF parameters to param. arrays
!  Input:
!    istate: electronic state of interest
!  Output:
!    k: force constants in Hatree
!    q0: rest values
!    n: multiplicities (for torsion) or dummy (for others)
!
!  ~2200 lines of params
!=================================================================
  subroutine assignparamuii_t(k, q0, n, istate)
  implicit none
  double precision :: k(ntermmax,napmax), q0(ntermmax,napmax)
  integer :: n(ntermmax,napmax), istate

  !--- unit of parameters have been coverted to hartree, Angstrom, rad
  select case(istate)
  case(1)
    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 1, anchor point 1)
    k (  1, 1) =  2.347640e+00
    q0(  1, 1) =  1.392339e+00
    n (  1, 1) =  0
    k (  2, 1) =  2.314035e+00
    q0(  2, 1) =  1.394857e+00
    n (  2, 1) =  0
    k (  3, 1) =  2.186652e+00
    q0(  3, 1) =  1.769756e+00
    n (  3, 1) =  0
    k (  4, 1) =  2.499818e+00
    q0(  4, 1) =  1.387612e+00
    n (  4, 1) =  0
    k (  5, 1) =  1.522459e+00
    q0(  5, 1) =  1.083546e+00
    n (  5, 1) =  0
    k (  6, 1) =  2.437438e+00
    q0(  6, 1) =  1.390686e+00
    n (  6, 1) =  0
    k (  7, 1) =  1.530598e+00
    q0(  7, 1) =  1.082797e+00
    n (  7, 1) =  0
    k (  8, 1) =  2.444289e+00
    q0(  8, 1) =  1.390135e+00
    n (  8, 1) =  0
    k (  9, 1) =  1.534995e+00
    q0(  9, 1) =  1.082358e+00
    n (  9, 1) =  0
    k ( 10, 1) =  2.474006e+00
    q0( 10, 1) =  1.389680e+00
    n ( 10, 1) =  0
    k ( 11, 1) =  1.527239e+00
    q0( 11, 1) =  1.083327e+00
    n ( 11, 1) =  0
    k ( 12, 1) =  1.515627e+00
    q0( 12, 1) =  1.084342e+00
    n ( 12, 1) =  0
    k ( 13, 1) =  2.284764e-01
    q0( 13, 1) =  2.097222e+00
    n ( 13, 1) =  0
    k ( 14, 1) =  1.533959e-01
    q0( 14, 1) =  2.083699e+00
    n ( 14, 1) =  0
    k ( 15, 1) =  2.325593e-01
    q0( 15, 1) =  2.091301e+00
    n ( 15, 1) =  0
    k ( 16, 1) =  1.502361e-01
    q0( 16, 1) =  2.103077e+00
    n ( 16, 1) =  0
    k ( 17, 1) =  2.115477E-01
    q0( 17, 1) =  1.673750E+00
    n ( 17, 1) =  0
    k ( 18, 1) =  1.158105e-01
    q0( 18, 1) =  2.089139e+00
    n ( 18, 1) =  0
    k ( 19, 1) =  3.286435e-01
    q0( 19, 1) =  2.044293e+00
    n ( 19, 1) =  0
    k ( 20, 1) =  2.304577e-01
    q0( 20, 1) =  2.103734e+00
    n ( 20, 1) =  0
    k ( 21, 1) =  1.667380e-01
    q0( 21, 1) =  2.081506e+00
    n ( 21, 1) =  0
    k ( 22, 1) =  1.816923e-01
    q0( 22, 1) =  2.100937e+00
    n ( 22, 1) =  0
    k ( 23, 1) =  2.399995e-01
    q0( 23, 1) =  2.083836e+00
    n ( 23, 1) =  0
    k ( 24, 1) =  1.697894e-01
    q0( 24, 1) =  2.099039e+00
    n ( 24, 1) =  0
    k ( 25, 1) =  1.645611e-01
    q0( 25, 1) =  2.097651e+00
    n ( 25, 1) =  0
    k ( 26, 1) =  2.280088e-01
    q0( 26, 1) =  2.102177e+00
    n ( 26, 1) =  0
    k ( 27, 1) =  1.652397e-01
    q0( 27, 1) =  2.096157e+00
    n ( 27, 1) =  0
    k ( 28, 1) =  1.699711e-01
    q0( 28, 1) =  2.100769e+00
    n ( 28, 1) =  0
    k ( 29, 1) =  1.798326e-01
    q0( 29, 1) =  2.091346e+00
    n ( 29, 1) =  0
    k ( 30, 1) =  3.314910e-01
    q0( 30, 1) =  2.148903e+00
    n ( 30, 1) =  0
    k ( 31, 1) =  1.660580e-01
    q0( 31, 1) =  2.086352e+00
    n ( 31, 1) =  0
    k ( 32, 1) =  7.291328e-03
    q0( 32, 1) =  0.000000e+00
    n ( 32, 1) =  2
    k ( 33, 1) =  1.485228e-02
    q0( 33, 1) =  0.000000e+00
    n ( 33, 1) =  2
    k ( 34, 1) =  5.852054e-03
    q0( 34, 1) =  0.000000e+00
    n ( 34, 1) =  2
    k ( 35, 1) =  1.423946e-02
    q0( 35, 1) =  0.000000e+00
    n ( 35, 1) =  2
    k ( 36, 1) =  1.335404e-03
    q0( 36, 1) =  0.000000e+00
    n ( 36, 1) =  2
    k ( 37, 1) =  1.443554e-02
    q0( 37, 1) =  0.000000e+00
    n ( 37, 1) =  2
    k ( 38, 1) =  4.167540e-03
    q0( 38, 1) =  0.000000e+00
    n ( 38, 1) =  2
    k ( 39, 1) =  1.130031e-02
    q0( 39, 1) =  0.000000e+00
    n ( 39, 1) =  2
    k ( 40, 1) =  7.946128e-04
    q0( 40, 1) =  0.000000e+00
    n ( 40, 1) =  2
    k ( 41, 1) =  2.883125e-03
    q0( 41, 1) =  0.000000e+00
    n ( 41, 1) =  2
    k ( 42, 1) =  5.849757e-03
    q0( 42, 1) =  0.000000e+00
    n ( 42, 1) =  2
    k ( 43, 1) =  1.137841e-02
    q0( 43, 1) =  0.000000e+00
    n ( 43, 1) =  2
    k ( 44, 1) =  1.212178e-02
    q0( 44, 1) =  0.000000e+00
    n ( 44, 1) =  2
    k ( 45, 1) =  1.197676e-02
    q0( 45, 1) =  0.000000e+00
    n ( 45, 1) =  2
    k ( 46, 1) =  1.060003e-02
    q0( 46, 1) =  0.000000e+00
    n ( 46, 1) =  2
    k ( 47, 1) =  7.246360e-03
    q0( 47, 1) =  0.000000e+00
    n ( 47, 1) =  2
    k ( 48, 1) =  1.259226e-02
    q0( 48, 1) =  0.000000e+00
    n ( 48, 1) =  2
    k ( 49, 1) =  1.234373e-02
    q0( 49, 1) =  0.000000e+00
    n ( 49, 1) =  2
    k ( 50, 1) =  1.486691e-03
    q0( 50, 1) =  0.000000e+00
    n ( 50, 1) =  2
    k ( 51, 1) =  2.225082e-03
    q0( 51, 1) =  0.000000e+00
    n ( 51, 1) =  2
    k ( 52, 1) =  5.211216e-18
    q0( 52, 1) =  0.000000e+00
    n ( 52, 1) =  2
    k ( 53, 1) =  3.648800e-18
    q0( 53, 1) =  0.000000e+00
    n ( 53, 1) =  2
    k ( 54, 1) =  2.415010e-03
    q0( 54, 1) =  0.000000e+00
    n ( 54, 1) =  2
    k ( 55, 1) =  8.554418e-04
    q0( 55, 1) =  0.000000e+00
    n ( 55, 1) =  2
    k ( 56, 1) =  1.780684e-01
    q0( 56, 1) =  4.665495e-05
    n ( 56, 1) =  0
    k ( 57, 1) =  1.474116e-01
    q0( 57, 1) =  1.344938e-05
    n ( 57, 1) =  0
    k ( 58, 1) =  1.946621e-01
    q0( 58, 1) =  5.444235e-05
    n ( 58, 1) =  0
    k ( 59, 1) =  1.965373e-01
    q0( 59, 1) =  1.732530e-06
    n ( 59, 1) =  0
    k ( 60, 1) =  1.837186e-01
    q0( 60, 1) =  6.018712e-06
    n ( 60, 1) =  0
    k ( 61, 1) =  1.755806e-01
    q0( 61, 1) =  3.538403e-05
    n ( 61, 1) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 1, anchor point 2)
    k (  1, 2) =  2.326234e+00
    q0(  1, 2) =  1.394732e+00
    n (  1, 2) =  0
    k (  2, 2) =  2.323994e+00
    q0(  2, 2) =  1.393427e+00
    n (  2, 2) =  0
    k (  3, 2) =  2.178704e+00
    q0(  3, 2) =  1.770804e+00
    n (  3, 2) =  0
    k (  4, 2) =  2.494225e+00
    q0(  4, 2) =  1.388303e+00
    n (  4, 2) =  0
    k (  5, 2) =  1.519302e+00
    q0(  5, 2) =  1.083863e+00
    n (  5, 2) =  0
    k (  6, 2) =  2.436564e+00
    q0(  6, 2) =  1.390784e+00
    n (  6, 2) =  0
    k (  7, 2) =  1.529534e+00
    q0(  7, 2) =  1.083067e+00
    n (  7, 2) =  0
    k (  8, 2) =  2.445907e+00
    q0(  8, 2) =  1.389932e+00
    n (  8, 2) =  0
    k (  9, 2) =  1.534695e+00
    q0(  9, 2) =  1.082414e+00
    n (  9, 2) =  0
    k ( 10, 2) =  2.476604e+00
    q0( 10, 2) =  1.389305e+00
    n ( 10, 2) =  0
    k ( 11, 2) =  1.528320e+00
    q0( 11, 2) =  1.083132e+00
    n ( 11, 2) =  0
    k ( 12, 2) =  1.516863e+00
    q0( 12, 2) =  1.084033e+00
    n ( 12, 2) =  0
    k ( 13, 2) =  2.399468e-01
    q0( 13, 2) =  2.095415e+00
    n ( 13, 2) =  0
    k ( 14, 2) =  1.570554e-01
    q0( 14, 2) =  2.090223e+00
    n ( 14, 2) =  0
    k ( 15, 2) =  2.340478e-01
    q0( 15, 2) =  2.094894e+00
    n ( 15, 2) =  0
    k ( 16, 2) =  1.513838e-01
    q0( 16, 2) =  2.096252e+00
    n ( 16, 2) =  0
    k ( 17, 2) =  2.179719E-01
    q0( 17, 2) =  1.649270E+00
    n ( 17, 2) =  0
    k ( 18, 2) =  1.170628e-01
    q0( 18, 2) =  2.086876e+00
    n ( 18, 2) =  0
    k ( 19, 2) =  2.947384e-01
    q0( 19, 2) =  2.056428e+00
    n ( 19, 2) =  0
    k ( 20, 2) =  2.270511e-01
    q0( 20, 2) =  2.102919e+00
    n ( 20, 2) =  0
    k ( 21, 2) =  1.653498e-01
    q0( 21, 2) =  2.084323e+00
    n ( 21, 2) =  0
    k ( 22, 2) =  1.791418e-01
    q0( 22, 2) =  2.098111e+00
    n ( 22, 2) =  0
    k ( 23, 2) =  2.404640e-01
    q0( 23, 2) =  2.083401e+00
    n ( 23, 2) =  0
    k ( 24, 2) =  1.703805e-01
    q0( 24, 2) =  2.100336e+00
    n ( 24, 2) =  0
    k ( 25, 2) =  1.651287e-01
    q0( 25, 2) =  2.096419e+00
    n ( 25, 2) =  0
    k ( 26, 2) =  2.273249e-01
    q0( 26, 2) =  2.103832e+00
    n ( 26, 2) =  0
    k ( 27, 2) =  1.653546e-01
    q0( 27, 2) =  2.096489e+00
    n ( 27, 2) =  0
    k ( 28, 2) =  1.696972e-01
    q0( 28, 2) =  2.099908e+00
    n ( 28, 2) =  0
    k ( 29, 2) =  1.788575e-01
    q0( 29, 2) =  2.092387e+00
    n ( 29, 2) =  0
    k ( 30, 2) =  3.188216e-01
    q0( 30, 2) =  2.140728e+00
    n ( 30, 2) =  0
    k ( 31, 2) =  1.657247e-01
    q0( 31, 2) =  2.083215e+00
    n ( 31, 2) =  0
    k ( 32, 2) =  6.098455e-03
    q0( 32, 2) =  0.000000e+00
    n ( 32, 2) =  2
    k ( 33, 2) =  1.235282e-02
    q0( 33, 2) =  0.000000e+00
    n ( 33, 2) =  2
    k ( 34, 2) =  6.171266e-03
    q0( 34, 2) =  0.000000e+00
    n ( 34, 2) =  2
    k ( 35, 2) =  1.245596e-02
    q0( 35, 2) =  0.000000e+00
    n ( 35, 2) =  2
    k ( 36, 2) =  3.148000e-19
    q0( 36, 2) =  0.000000e+00
    n ( 36, 2) =  2
    k ( 37, 2) =  1.202348e-02
    q0( 37, 2) =  0.000000e+00
    n ( 37, 2) =  2
    k ( 38, 2) =  5.387721e-03
    q0( 38, 2) =  0.000000e+00
    n ( 38, 2) =  2
    k ( 39, 2) =  1.122949e-02
    q0( 39, 2) =  0.000000e+00
    n ( 39, 2) =  2
    k ( 40, 2) =  8.192029e-19
    q0( 40, 2) =  0.000000e+00
    n ( 40, 2) =  2
    k ( 41, 2) =  1.035050e-02
    q0( 41, 2) =  0.000000e+00
    n ( 41, 2) =  2
    k ( 42, 2) =  5.409366e-03
    q0( 42, 2) =  0.000000e+00
    n ( 42, 2) =  2
    k ( 43, 2) =  1.258348e-02
    q0( 43, 2) =  0.000000e+00
    n ( 43, 2) =  2
    k ( 44, 2) =  1.252389e-02
    q0( 44, 2) =  0.000000e+00
    n ( 44, 2) =  2
    k ( 45, 2) =  1.283397e-02
    q0( 45, 2) =  0.000000e+00
    n ( 45, 2) =  2
    k ( 46, 2) =  1.196208e-02
    q0( 46, 2) =  0.000000e+00
    n ( 46, 2) =  2
    k ( 47, 2) =  1.404531e-02
    q0( 47, 2) =  0.000000e+00
    n ( 47, 2) =  2
    k ( 48, 2) =  1.191647e-02
    q0( 48, 2) =  0.000000e+00
    n ( 48, 2) =  2
    k ( 49, 2) =  1.180136e-02
    q0( 49, 2) =  0.000000e+00
    n ( 49, 2) =  2
    k ( 50, 2) =  2.400516e-03
    q0( 50, 2) =  0.000000e+00
    n ( 50, 2) =  2
    k ( 51, 2) =  1.889423e-03
    q0( 51, 2) =  0.000000e+00
    n ( 51, 2) =  2
    k ( 52, 2) = -9.211373e-20
    q0( 52, 2) =  0.000000e+00
    n ( 52, 2) =  2
    k ( 53, 2) =  2.457347e-04
    q0( 53, 2) =  0.000000e+00
    n ( 53, 2) =  2
    k ( 54, 2) =  2.107653e-03
    q0( 54, 2) =  0.000000e+00
    n ( 54, 2) =  2
    k ( 55, 2) =  1.923286e-03
    q0( 55, 2) =  0.000000e+00
    n ( 55, 2) =  2
    k ( 56, 2) =  1.757168e-01
    q0( 56, 2) =  1.884472e-06
    n ( 56, 2) =  0
    k ( 57, 2) =  1.557202e-01
    q0( 57, 2) =  1.496890e-06
    n ( 57, 2) =  0
    k ( 58, 2) =  1.829432e-01
    q0( 58, 2) =  4.885969e-06
    n ( 58, 2) =  0
    k ( 59, 2) =  2.202718e-01
    q0( 59, 2) =  3.337738e-07
    n ( 59, 2) =  0
    k ( 60, 2) =  1.638524e-01
    q0( 60, 2) =  5.918230e-07
    n ( 60, 2) =  0
    k ( 61, 2) =  1.629773e-01
    q0( 61, 2) =  3.310474e-07
    n ( 61, 2) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 1, anchor point 3)
    k (  1, 3) =  2.279548e+00
    q0(  1, 3) =  1.397059e+00
    n (  1, 3) =  0
    k (  2, 3) =  2.317374e+00
    q0(  2, 3) =  1.392386e+00
    n (  2, 3) =  0
    k (  3, 3) =  2.212747e+00
    q0(  3, 3) =  1.770000e+00
    n (  3, 3) =  0
    k (  4, 3) =  2.492692e+00
    q0(  4, 3) =  1.388585e+00
    n (  4, 3) =  0
    k (  5, 3) =  1.516551e+00
    q0(  5, 3) =  1.084130e+00
    n (  5, 3) =  0
    k (  6, 3) =  2.435536e+00
    q0(  6, 3) =  1.390816e+00
    n (  6, 3) =  0
    k (  7, 3) =  1.529162e+00
    q0(  7, 3) =  1.083228e+00
    n (  7, 3) =  0
    k (  8, 3) =  2.447353e+00
    q0(  8, 3) =  1.389702e+00
    n (  8, 3) =  0
    k (  9, 3) =  1.534832e+00
    q0(  9, 3) =  1.082431e+00
    n (  9, 3) =  0
    k ( 10, 3) =  2.478740e+00
    q0( 10, 3) =  1.389136e+00
    n ( 10, 3) =  0
    k ( 11, 3) =  1.528954e+00
    q0( 11, 3) =  1.083011e+00
    n ( 11, 3) =  0
    k ( 12, 3) =  1.519137e+00
    q0( 12, 3) =  1.083758e+00
    n ( 12, 3) =  0
    k ( 13, 3) =  2.425219e-01
    q0( 13, 3) =  2.094562e+00
    n ( 13, 3) =  0
    k ( 14, 3) =  1.554735e-01
    q0( 14, 3) =  2.094952e+00
    n ( 14, 3) =  0
    k ( 15, 3) =  2.347723e-01
    q0( 15, 3) =  2.097260e+00
    n ( 15, 3) =  0
    k ( 16, 3) =  1.496334e-01
    q0( 16, 3) =  2.090187e+00
    n ( 16, 3) =  0
    k ( 17, 3) =  1.917046E-01
    q0( 17, 3) =  1.669610E+00
    n ( 17, 3) =  0
    k ( 18, 3) =  1.089022e-01
    q0( 18, 3) =  2.084751e+00
    n ( 18, 3) =  0
    k ( 19, 3) =  2.921064e-01
    q0( 19, 3) =  2.051540e+00
    n ( 19, 3) =  0
    k ( 20, 3) =  2.240999e-01
    q0( 20, 3) =  2.102656e+00
    n ( 20, 3) =  0
    k ( 21, 3) =  1.642951e-01
    q0( 21, 3) =  2.085428e+00
    n ( 21, 3) =  0
    k ( 22, 3) =  1.810727e-01
    q0( 22, 3) =  2.095241e+00
    n ( 22, 3) =  0
    k ( 23, 3) =  2.406570e-01
    q0( 23, 3) =  2.082864e+00
    n ( 23, 3) =  0
    k ( 24, 3) =  1.708084e-01
    q0( 24, 3) =  2.101158e+00
    n ( 24, 3) =  0
    k ( 25, 3) =  1.651373e-01
    q0( 25, 3) =  2.096003e+00
    n ( 25, 3) =  0
    k ( 26, 3) =  2.268501e-01
    q0( 26, 3) =  2.105195e+00
    n ( 26, 3) =  0
    k ( 27, 3) =  1.652394e-01
    q0( 27, 3) =  2.096628e+00
    n ( 27, 3) =  0
    k ( 28, 3) =  1.695634e-01
    q0( 28, 3) =  2.099655e+00
    n ( 28, 3) =  0
    k ( 29, 3) =  1.799576e-01
    q0( 29, 3) =  2.094771e+00
    n ( 29, 3) =  0
    k ( 30, 3) =  3.067507e-01
    q0( 30, 3) =  2.148098e+00
    n ( 30, 3) =  0
    k ( 31, 3) =  1.653963e-01
    q0( 31, 3) =  2.081003e+00
    n ( 31, 3) =  0
    k ( 32, 3) =  6.874204e-03
    q0( 32, 3) =  0.000000e+00
    n ( 32, 3) =  2
    k ( 33, 3) =  1.211864e-02
    q0( 33, 3) =  0.000000e+00
    n ( 33, 3) =  2
    k ( 34, 3) =  6.742790e-03
    q0( 34, 3) =  0.000000e+00
    n ( 34, 3) =  2
    k ( 35, 3) =  1.209450e-02
    q0( 35, 3) =  0.000000e+00
    n ( 35, 3) =  2
    k ( 36, 3) =  6.223053e-04
    q0( 36, 3) =  0.000000e+00
    n ( 36, 3) =  2
    k ( 37, 3) =  8.099623e-03
    q0( 37, 3) =  0.000000e+00
    n ( 37, 3) =  2
    k ( 38, 3) =  6.832013e-03
    q0( 38, 3) =  0.000000e+00
    n ( 38, 3) =  2
    k ( 39, 3) =  1.322854e-02
    q0( 39, 3) =  0.000000e+00
    n ( 39, 3) =  2
    k ( 40, 3) =  9.009229e-20
    q0( 40, 3) =  0.000000e+00
    n ( 40, 3) =  2
    k ( 41, 3) =  6.542268e-03
    q0( 41, 3) =  0.000000e+00
    n ( 41, 3) =  2
    k ( 42, 3) =  4.036241e-03
    q0( 42, 3) =  0.000000e+00
    n ( 42, 3) =  2
    k ( 43, 3) =  1.284694e-02
    q0( 43, 3) =  0.000000e+00
    n ( 43, 3) =  2
    k ( 44, 3) =  1.202531e-02
    q0( 44, 3) =  0.000000e+00
    n ( 44, 3) =  2
    k ( 45, 3) =  1.457981e-02
    q0( 45, 3) =  0.000000e+00
    n ( 45, 3) =  2
    k ( 46, 3) =  1.171550e-02
    q0( 46, 3) =  0.000000e+00
    n ( 46, 3) =  2
    k ( 47, 3) =  1.970903e-02
    q0( 47, 3) =  0.000000e+00
    n ( 47, 3) =  2
    k ( 48, 3) =  1.200196e-02
    q0( 48, 3) =  0.000000e+00
    n ( 48, 3) =  2
    k ( 49, 3) =  9.834332e-03
    q0( 49, 3) =  0.000000e+00
    n ( 49, 3) =  2
    k ( 50, 3) =  2.437633e-03
    q0( 50, 3) =  0.000000e+00
    n ( 50, 3) =  2
    k ( 51, 3) =  2.415636e-03
    q0( 51, 3) =  0.000000e+00
    n ( 51, 3) =  2
    k ( 52, 3) =  1.581733e-18
    q0( 52, 3) =  0.000000e+00
    n ( 52, 3) =  2
    k ( 53, 3) =  4.720208e-03
    q0( 53, 3) =  0.000000e+00
    n ( 53, 3) =  2
    k ( 54, 3) =  3.713412e-03
    q0( 54, 3) =  0.000000e+00
    n ( 54, 3) =  2
    k ( 55, 3) =  2.428604e-03
    q0( 55, 3) =  0.000000e+00
    n ( 55, 3) =  2
    k ( 56, 3) =  1.128810e-01
    q0( 56, 3) =  3.497804e-05
    n ( 56, 3) =  0
    k ( 57, 3) =  1.442682e-01
    q0( 57, 3) =  3.732224e-05
    n ( 57, 3) =  0
    k ( 58, 3) =  1.649825e-01
    q0( 58, 3) =  3.340977e-05
    n ( 58, 3) =  0
    k ( 59, 3) =  7.012445e-02
    q0( 59, 3) =  1.107913e-05
    n ( 59, 3) =  0
    k ( 60, 3) =  1.571235e-01
    q0( 60, 3) =  3.445407e-06
    n ( 60, 3) =  0
    k ( 61, 3) =  1.552393e-01
    q0( 61, 3) =  1.812268e-05
    n ( 61, 3) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 1, anchor point 4)
    k (  1, 4) =  2.231258e+00
    q0(  1, 4) =  1.399828e+00
    n (  1, 4) =  0
    k (  2, 4) =  2.310947e+00
    q0(  2, 4) =  1.392583e+00
    n (  2, 4) =  0
    k (  3, 4) =  2.179581e+00
    q0(  3, 4) =  1.763637e+00
    n (  3, 4) =  0
    k (  4, 4) =  2.486675e+00
    q0(  4, 4) =  1.388517e+00
    n (  4, 4) =  0
    k (  5, 4) =  1.513275e+00
    q0(  5, 4) =  1.084104e+00
    n (  5, 4) =  0
    k (  6, 4) =  2.437557e+00
    q0(  6, 4) =  1.390523e+00
    n (  6, 4) =  0
    k (  7, 4) =  1.531134e+00
    q0(  7, 4) =  1.083223e+00
    n (  7, 4) =  0
    k (  8, 4) =  2.447727e+00
    q0(  8, 4) =  1.389509e+00
    n (  8, 4) =  0
    k (  9, 4) =  1.535626e+00
    q0(  9, 4) =  1.082373e+00
    n (  9, 4) =  0
    k ( 10, 4) =  2.480480e+00
    q0( 10, 4) =  1.388825e+00
    n ( 10, 4) =  0
    k ( 11, 4) =  1.529698e+00
    q0( 11, 4) =  1.082983e+00
    n ( 11, 4) =  0
    k ( 12, 4) =  1.519050e+00
    q0( 12, 4) =  1.083605e+00
    n ( 12, 4) =  0
    k ( 13, 4) =  2.560005e-01
    q0( 13, 4) =  2.093942e+00
    n ( 13, 4) =  0
    k ( 14, 4) =  1.566703e-01
    q0( 14, 4) =  2.098527e+00
    n ( 14, 4) =  0
    k ( 15, 4) =  2.413345e-01
    q0( 15, 4) =  2.098190e+00
    n ( 15, 4) =  0
    k ( 16, 4) =  1.479298e-01
    q0( 16, 4) =  2.082111e+00
    n ( 16, 4) =  0
    k ( 17, 4) =  1.464658E-01
    q0( 17, 4) =  1.705890E+00
    n ( 17, 4) =  0
    k ( 18, 4) =  1.019613e-01
    q0( 18, 4) =  2.082829e+00
    n ( 18, 4) =  0
    k ( 19, 4) =  2.501192e-01
    q0( 19, 4) =  2.037807e+00
    n ( 19, 4) =  0
    k ( 20, 4) =  2.145793e-01
    q0( 20, 4) =  2.104071e+00
    n ( 20, 4) =  0
    k ( 21, 4) =  1.603081e-01
    q0( 21, 4) =  2.083135e+00
    n ( 21, 4) =  0
    k ( 22, 4) =  1.769633e-01
    q0( 22, 4) =  2.092175e+00
    n ( 22, 4) =  0
    k ( 23, 4) =  2.422186e-01
    q0( 23, 4) =  2.082103e+00
    n ( 23, 4) =  0
    k ( 24, 4) =  1.717612e-01
    q0( 24, 4) =  2.101302e+00
    n ( 24, 4) =  0
    k ( 25, 4) =  1.665253e-01
    q0( 25, 4) =  2.096890e+00
    n ( 25, 4) =  0
    k ( 26, 4) =  2.240627e-01
    q0( 26, 4) =  2.106130e+00
    n ( 26, 4) =  0
    k ( 27, 4) =  1.655003e-01
    q0( 27, 4) =  2.096842e+00
    n ( 27, 4) =  0
    k ( 28, 4) =  1.693364e-01
    q0( 28, 4) =  2.100268e+00
    n ( 28, 4) =  0
    k ( 29, 4) =  1.776005e-01
    q0( 29, 4) =  2.101989e+00
    n ( 29, 4) =  0
    k ( 30, 4) =  2.610728e-01
    q0( 30, 4) =  2.163951e+00
    n ( 30, 4) =  0
    k ( 31, 4) =  1.640269e-01
    q0( 31, 4) =  2.079842e+00
    n ( 31, 4) =  0
    k ( 32, 4) =  5.458822e-03
    q0( 32, 4) =  0.000000e+00
    n ( 32, 4) =  2
    k ( 33, 4) =  1.193116e-02
    q0( 33, 4) =  0.000000e+00
    n ( 33, 4) =  2
    k ( 34, 4) =  6.519191e-03
    q0( 34, 4) =  0.000000e+00
    n ( 34, 4) =  2
    k ( 35, 4) =  1.154422e-02
    q0( 35, 4) =  0.000000e+00
    n ( 35, 4) =  2
    k ( 36, 4) = -8.287131e-19
    q0( 36, 4) =  0.000000e+00
    n ( 36, 4) =  2
    k ( 37, 4) =  7.951998e-03
    q0( 37, 4) =  0.000000e+00
    n ( 37, 4) =  2
    k ( 38, 4) =  6.361505e-03
    q0( 38, 4) =  0.000000e+00
    n ( 38, 4) =  2
    k ( 39, 4) =  1.219586e-02
    q0( 39, 4) =  0.000000e+00
    n ( 39, 4) =  2
    k ( 40, 4) =  1.007202e-18
    q0( 40, 4) =  0.000000e+00
    n ( 40, 4) =  2
    k ( 41, 4) =  6.021719e-03
    q0( 41, 4) =  0.000000e+00
    n ( 41, 4) =  2
    k ( 42, 4) =  3.978436e-03
    q0( 42, 4) =  0.000000e+00
    n ( 42, 4) =  2
    k ( 43, 4) =  1.328427e-02
    q0( 43, 4) =  0.000000e+00
    n ( 43, 4) =  2
    k ( 44, 4) =  1.223164e-02
    q0( 44, 4) =  0.000000e+00
    n ( 44, 4) =  2
    k ( 45, 4) =  1.438962e-02
    q0( 45, 4) =  0.000000e+00
    n ( 45, 4) =  2
    k ( 46, 4) =  1.165162e-02
    q0( 46, 4) =  0.000000e+00
    n ( 46, 4) =  2
    k ( 47, 4) =  1.876940e-02
    q0( 47, 4) =  0.000000e+00
    n ( 47, 4) =  2
    k ( 48, 4) =  1.207975e-02
    q0( 48, 4) =  0.000000e+00
    n ( 48, 4) =  2
    k ( 49, 4) =  1.064373e-02
    q0( 49, 4) =  0.000000e+00
    n ( 49, 4) =  2
    k ( 50, 4) =  2.562274e-03
    q0( 50, 4) =  0.000000e+00
    n ( 50, 4) =  2
    k ( 51, 4) =  2.268551e-03
    q0( 51, 4) =  0.000000e+00
    n ( 51, 4) =  2
    k ( 52, 4) =  2.839749e-19
    q0( 52, 4) =  0.000000e+00
    n ( 52, 4) =  2
    k ( 53, 4) =  3.300967e-03
    q0( 53, 4) =  0.000000e+00
    n ( 53, 4) =  2
    k ( 54, 4) =  3.382732e-03
    q0( 54, 4) =  0.000000e+00
    n ( 54, 4) =  2
    k ( 55, 4) =  2.054074e-03
    q0( 55, 4) =  0.000000e+00
    n ( 55, 4) =  2
    k ( 56, 4) =  1.247386e-01
    q0( 56, 4) =  5.726373e-06
    n ( 56, 4) =  0
    k ( 57, 4) =  1.446992e-01
    q0( 57, 4) =  3.580872e-05
    n ( 57, 4) =  0
    k ( 58, 4) =  1.731829e-01
    q0( 58, 4) =  3.806675e-05
    n ( 58, 4) =  0
    k ( 59, 4) =  9.069293e-02
    q0( 59, 4) =  8.717881e-06
    n ( 59, 4) =  0
    k ( 60, 4) =  1.603129e-01
    q0( 60, 4) =  5.820413e-07
    n ( 60, 4) =  0
    k ( 61, 4) =  1.573913e-01
    q0( 61, 4) =  3.002728e-05
    n ( 61, 4) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 1, anchor point 5)
    k (  1, 5) =  2.224225e+00
    q0(  1, 5) =  1.398904e+00
    n (  1, 5) =  0
    k (  2, 5) =  2.277618e+00
    q0(  2, 5) =  1.395246e+00
    n (  2, 5) =  0
    k (  3, 5) =  2.146771e+00
    q0(  3, 5) =  1.758248e+00
    n (  3, 5) =  0
    k (  4, 5) =  2.486770e+00
    q0(  4, 5) =  1.388296e+00
    n (  4, 5) =  0
    k (  5, 5) =  1.517587e+00
    q0(  5, 5) =  1.083700e+00
    n (  5, 5) =  0
    k (  6, 5) =  2.440282e+00
    q0(  6, 5) =  1.390368e+00
    n (  6, 5) =  0
    k (  7, 5) =  1.531026e+00
    q0(  7, 5) =  1.083155e+00
    n (  7, 5) =  0
    k (  8, 5) =  2.446709e+00
    q0(  8, 5) =  1.389623e+00
    n (  8, 5) =  0
    k (  9, 5) =  1.535490e+00
    q0(  9, 5) =  1.082318e+00
    n (  9, 5) =  0
    k ( 10, 5) =  2.479237e+00
    q0( 10, 5) =  1.388883e+00
    n ( 10, 5) =  0
    k ( 11, 5) =  1.529385e+00
    q0( 11, 5) =  1.083102e+00
    n ( 11, 5) =  0
    k ( 12, 5) =  1.522250e+00
    q0( 12, 5) =  1.083305e+00
    n ( 12, 5) =  0
    k ( 13, 5) =  2.527379e-01
    q0( 13, 5) =  2.094879e+00
    n ( 13, 5) =  0
    k ( 14, 5) =  1.504106e-01
    q0( 14, 5) =  2.094911e+00
    n ( 14, 5) =  0
    k ( 15, 5) =  2.423258e-01
    q0( 15, 5) =  2.095669e+00
    n ( 15, 5) =  0
    k ( 16, 5) =  1.447724e-01
    q0( 16, 5) =  2.090052e+00
    n ( 16, 5) =  0
    k ( 17, 5) =  1.028653E-01
    q0( 17, 5) =  1.751410E+00
    n ( 17, 5) =  0
    k ( 18, 5) =  1.151130e-01
    q0( 18, 5) =  2.083297e+00
    n ( 18, 5) =  0
    k ( 19, 5) =  2.394527e-01
    q0( 19, 5) =  2.058934e+00
    n ( 19, 5) =  0
    k ( 20, 5) =  2.142067e-01
    q0( 20, 5) =  2.105212e+00
    n ( 20, 5) =  0
    k ( 21, 5) =  1.600837e-01
    q0( 21, 5) =  2.081235e+00
    n ( 21, 5) =  0
    k ( 22, 5) =  1.790453e-01
    q0( 22, 5) =  2.094100e+00
    n ( 22, 5) =  0
    k ( 23, 5) =  2.430724e-01
    q0( 23, 5) =  2.081628e+00
    n ( 23, 5) =  0
    k ( 24, 5) =  1.718860e-01
    q0( 24, 5) =  2.101061e+00
    n ( 24, 5) =  0
    k ( 25, 5) =  1.666574e-01
    q0( 25, 5) =  2.097283e+00
    n ( 25, 5) =  0
    k ( 26, 5) =  2.211735e-01
    q0( 26, 5) =  2.106648e+00
    n ( 26, 5) =  0
    k ( 27, 5) =  1.652871e-01
    q0( 27, 5) =  2.096470e+00
    n ( 27, 5) =  0
    k ( 28, 5) =  1.698874e-01
    q0( 28, 5) =  2.100954e+00
    n ( 28, 5) =  0
    k ( 29, 5) =  1.796925e-01
    q0( 29, 5) =  2.097590e+00
    n ( 29, 5) =  0
    k ( 30, 5) =  2.571206e-01
    q0( 30, 5) =  2.142129e+00
    n ( 30, 5) =  0
    k ( 31, 5) =  1.627543e-01
    q0( 31, 5) =  2.080262e+00
    n ( 31, 5) =  0
    k ( 32, 5) =  3.427996e-03
    q0( 32, 5) =  0.000000e+00
    n ( 32, 5) =  2
    k ( 33, 5) =  1.358425e-02
    q0( 33, 5) =  0.000000e+00
    n ( 33, 5) =  2
    k ( 34, 5) =  7.545969e-04
    q0( 34, 5) =  0.000000e+00
    n ( 34, 5) =  2
    k ( 35, 5) =  1.035871e-02
    q0( 35, 5) =  0.000000e+00
    n ( 35, 5) =  2
    k ( 36, 5) =  4.414514e-03
    q0( 36, 5) =  0.000000e+00
    n ( 36, 5) =  2
    k ( 37, 5) =  1.015838e-02
    q0( 37, 5) =  0.000000e+00
    n ( 37, 5) =  2
    k ( 38, 5) =  5.766217e-03
    q0( 38, 5) =  0.000000e+00
    n ( 38, 5) =  2
    k ( 39, 5) =  9.624610e-03
    q0( 39, 5) =  0.000000e+00
    n ( 39, 5) =  2
    k ( 40, 5) =  8.017042e-03
    q0( 40, 5) =  0.000000e+00
    n ( 40, 5) =  2
    k ( 41, 5) =  1.793989e-18
    q0( 41, 5) =  0.000000e+00
    n ( 41, 5) =  2
    k ( 42, 5) =  1.093724e-02
    q0( 42, 5) =  0.000000e+00
    n ( 42, 5) =  2
    k ( 43, 5) =  1.562397e-02
    q0( 43, 5) =  0.000000e+00
    n ( 43, 5) =  2
    k ( 44, 5) =  1.032460e-02
    q0( 44, 5) =  0.000000e+00
    n ( 44, 5) =  2
    k ( 45, 5) =  1.469065e-02
    q0( 45, 5) =  0.000000e+00
    n ( 45, 5) =  2
    k ( 46, 5) =  1.171390e-02
    q0( 46, 5) =  0.000000e+00
    n ( 46, 5) =  2
    k ( 47, 5) =  2.131533e-02
    q0( 47, 5) =  0.000000e+00
    n ( 47, 5) =  2
    k ( 48, 5) =  1.667150e-02
    q0( 48, 5) =  0.000000e+00
    n ( 48, 5) =  2
    k ( 49, 5) =  1.240539e-02
    q0( 49, 5) =  0.000000e+00
    n ( 49, 5) =  2
    k ( 50, 5) =  4.233082e-03
    q0( 50, 5) =  0.000000e+00
    n ( 50, 5) =  2
    k ( 51, 5) =  3.376412e-03
    q0( 51, 5) =  0.000000e+00
    n ( 51, 5) =  2
    k ( 52, 5) =  9.709543e-19
    q0( 52, 5) =  0.000000e+00
    n ( 52, 5) =  2
    k ( 53, 5) =  4.904397e-03
    q0( 53, 5) =  0.000000e+00
    n ( 53, 5) =  2
    k ( 54, 5) =  5.036038e-03
    q0( 54, 5) =  0.000000e+00
    n ( 54, 5) =  2
    k ( 55, 5) =  3.263352e-03
    q0( 55, 5) =  0.000000e+00
    n ( 55, 5) =  2
    k ( 56, 5) =  8.349949e-02
    q0( 56, 5) =  2.016920e-07
    n ( 56, 5) =  0
    k ( 57, 5) =  1.075575e-01
    q0( 57, 5) =  2.719359e-07
    n ( 57, 5) =  0
    k ( 58, 5) =  1.368283e-01
    q0( 58, 5) =  6.776049e-06
    n ( 58, 5) =  0
    k ( 59, 5) =  1.347852e-01
    q0( 59, 5) =  2.906716e-06
    n ( 59, 5) =  0
    k ( 60, 5) =  1.247791e-01
    q0( 60, 5) =  2.405811e-06
    n ( 60, 5) =  0
    k ( 61, 5) =  1.196619e-01
    q0( 61, 5) =  1.847443e-06
    n ( 61, 5) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 1, anchor point 6)
    k (  1, 6) =  1.030521e+00
    q0(  1, 6) =  1.406720e+00
    n (  1, 6) =  0
    k (  2, 6) =  9.493466e-01
    q0(  2, 6) =  1.406245e+00
    n (  2, 6) =  0
    k (  3, 6) =  5.410092e-01
    q0(  3, 6) =  1.759236e+00
    n (  3, 6) =  0
    k (  4, 6) =  1.016899e+00
    q0(  4, 6) =  1.389937e+00
    n (  4, 6) =  0
    k (  5, 6) =  4.036298e-01
    q0(  5, 6) =  1.085046e+00
    n (  5, 6) =  0
    k (  6, 6) =  1.067908e+00
    q0(  6, 6) =  1.394944e+00
    n (  6, 6) =  0
    k (  7, 6) =  4.051993e-01
    q0(  7, 6) =  1.085787e+00
    n (  7, 6) =  0
    k (  8, 6) =  1.063125e+00
    q0(  8, 6) =  1.394871e+00
    n (  8, 6) =  0
    k (  9, 6) =  4.043513e-01
    q0(  9, 6) =  1.084677e+00
    n (  9, 6) =  0
    k ( 10, 6) =  1.034509e+00
    q0( 10, 6) =  1.390205e+00
    n ( 10, 6) =  0
    k ( 11, 6) =  4.039146e-01
    q0( 11, 6) =  1.085753e+00
    n ( 11, 6) =  0
    k ( 12, 6) =  3.926463e-01
    q0( 12, 6) =  1.084651e+00
    n ( 12, 6) =  0
    k ( 13, 6) =  4.420621e-01
    q0( 13, 6) =  2.093925e+00
    n ( 13, 6) =  0
    k ( 14, 6) =  5.948993e-01
    q0( 14, 6) =  2.096981e+00
    n ( 14, 6) =  0
    k ( 15, 6) =  3.780707e-01
    q0( 15, 6) =  2.093822e+00
    n ( 15, 6) =  0
    k ( 16, 6) =  6.755354e-01
    q0( 16, 6) =  2.096793e+00
    n ( 16, 6) =  0
    k ( 17, 6) =  1.961200E-02
    q0( 17, 6) =  1.839430E+00
    n ( 17, 6) =  0
    k ( 18, 6) =  3.823444e-01
    q0( 18, 6) =  2.082663e+00
    n ( 18, 6) =  0
    k ( 19, 6) =  3.668207e-01
    q0( 19, 6) =  2.094545e+00
    n ( 19, 6) =  0
    k ( 20, 6) =  3.767371e-01
    q0( 20, 6) =  2.108826e+00
    n ( 20, 6) =  0
    k ( 21, 6) =  5.824497e-01
    q0( 21, 6) =  2.078964e+00
    n ( 21, 6) =  0
    k ( 22, 6) =  6.857721e-01
    q0( 22, 6) =  2.092893e+00
    n ( 22, 6) =  0
    k ( 23, 6) =  3.593871e-01
    q0( 23, 6) =  2.081971e+00
    n ( 23, 6) =  0
    k ( 24, 6) =  6.044578e-01
    q0( 24, 6) =  2.101155e+00
    n ( 24, 6) =  0
    k ( 25, 6) =  6.870095e-01
    q0( 25, 6) =  2.096106e+00
    n ( 25, 6) =  0
    k ( 26, 6) =  3.928344e-01
    q0( 26, 6) =  2.109174e+00
    n ( 26, 6) =  0
    k ( 27, 6) =  6.073915e-01
    q0( 27, 6) =  2.096279e+00
    n ( 27, 6) =  0
    k ( 28, 6) =  6.824865e-01
    q0( 28, 6) =  2.100765e+00
    n ( 28, 6) =  0
    k ( 29, 6) =  6.304806e-01
    q0( 29, 6) =  2.093294e+00
    n ( 29, 6) =  0
    k ( 30, 6) =  4.529314e-01
    q0( 30, 6) =  2.106876e+00
    n ( 30, 6) =  0
    k ( 31, 6) =  6.596825e-01
    q0( 31, 6) =  2.078452e+00
    n ( 31, 6) =  0
    k ( 32, 6) =  7.018616e-03
    q0( 32, 6) =  0.000000e+00
    n ( 32, 6) =  2
    k ( 33, 6) =  1.134519e-02
    q0( 33, 6) =  0.000000e+00
    n ( 33, 6) =  2
    k ( 34, 6) =  7.809210e-03
    q0( 34, 6) =  0.000000e+00
    n ( 34, 6) =  2
    k ( 35, 6) =  1.234957e-02
    q0( 35, 6) =  0.000000e+00
    n ( 35, 6) =  2
    k ( 36, 6) = -1.109128e-18
    q0( 36, 6) =  0.000000e+00
    n ( 36, 6) =  2
    k ( 37, 6) =  1.039982e-02
    q0( 37, 6) =  0.000000e+00
    n ( 37, 6) =  2
    k ( 38, 6) =  5.313863e-03
    q0( 38, 6) =  0.000000e+00
    n ( 38, 6) =  2
    k ( 39, 6) =  1.157038e-02
    q0( 39, 6) =  0.000000e+00
    n ( 39, 6) =  2
    k ( 40, 6) =  5.990963e-19
    q0( 40, 6) =  0.000000e+00
    n ( 40, 6) =  2
    k ( 41, 6) =  1.446517e-02
    q0( 41, 6) =  0.000000e+00
    n ( 41, 6) =  2
    k ( 42, 6) =  4.434246e-03
    q0( 42, 6) =  0.000000e+00
    n ( 42, 6) =  2
    k ( 43, 6) =  1.254120e-02
    q0( 43, 6) =  0.000000e+00
    n ( 43, 6) =  2
    k ( 44, 6) =  1.358333e-02
    q0( 44, 6) =  0.000000e+00
    n ( 44, 6) =  2
    k ( 45, 6) =  1.388110e-02
    q0( 45, 6) =  0.000000e+00
    n ( 45, 6) =  2
    k ( 46, 6) =  1.280230e-02
    q0( 46, 6) =  0.000000e+00
    n ( 46, 6) =  2
    k ( 47, 6) =  1.212896e-02
    q0( 47, 6) =  0.000000e+00
    n ( 47, 6) =  2
    k ( 48, 6) =  1.163289e-02
    q0( 48, 6) =  0.000000e+00
    n ( 48, 6) =  2
    k ( 49, 6) =  1.127802e-02
    q0( 49, 6) =  0.000000e+00
    n ( 49, 6) =  2
    k ( 50, 6) =  1.767183e-03
    q0( 50, 6) =  0.000000e+00
    n ( 50, 6) =  2
    k ( 51, 6) =  1.954895e-03
    q0( 51, 6) =  0.000000e+00
    n ( 51, 6) =  2
    k ( 52, 6) =  5.728171e-05
    q0( 52, 6) =  0.000000e+00
    n ( 52, 6) =  2
    k ( 53, 6) =  8.194179e-04
    q0( 53, 6) =  0.000000e+00
    n ( 53, 6) =  2
    k ( 54, 6) =  2.032833e-03
    q0( 54, 6) =  0.000000e+00
    n ( 54, 6) =  2
    k ( 55, 6) =  1.942172e-03
    q0( 55, 6) =  0.000000e+00
    n ( 55, 6) =  2
    k ( 56, 6) =  1.222414e-01
    q0( 56, 6) =  7.660981e-07
    n ( 56, 6) =  0
    k ( 57, 6) =  1.373986e-01
    q0( 57, 6) =  5.838531e-07
    n ( 57, 6) =  0
    k ( 58, 6) =  1.536470e-01
    q0( 58, 6) =  1.299623e-06
    n ( 58, 6) =  0
    k ( 59, 6) =  2.100033e-01
    q0( 59, 6) =  3.272801e-07
    n ( 59, 6) =  0
    k ( 60, 6) =  1.366956e-01
    q0( 60, 6) =  2.057717e-07
    n ( 60, 6) =  0
    k ( 61, 6) =  1.527704e-01
    q0( 61, 6) =  1.846245e-07
    n ( 61, 6) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 1, anchor point 6)
    k (  1, 7) =  2.091985e+00
    q0(  1, 7) =  1.406498e+00
    n (  1, 7) =  0
    k (  2, 7) =  2.091984e+00
    q0(  2, 7) =  1.406498e+00
    n (  2, 7) =  0
    k (  3, 7) =  1.958110e+00
    q0(  3, 7) =  1.758607e+00
    n (  3, 7) =  0
    k (  4, 7) =  2.419031e+00
    q0(  4, 7) =  1.390021e+00
    n (  4, 7) =  0
    k (  5, 7) =  1.482304e+00
    q0(  5, 7) =  1.084994e+00
    n (  5, 7) =  0
    k (  6, 7) =  2.329915e+00
    q0(  6, 7) =  1.394996e+00
    n (  6, 7) =  0
    k (  7, 7) =  1.484484e+00
    q0(  7, 7) =  1.085765e+00
    n (  7, 7) =  0
    k (  8, 7) =  2.329915e+00
    q0(  8, 7) =  1.394996e+00
    n (  8, 7) =  0
    k (  9, 7) =  1.490964e+00
    q0(  9, 7) =  1.084701e+00
    n (  9, 7) =  0
    k ( 10, 7) =  2.419034e+00
    q0( 10, 7) =  1.390021e+00
    n ( 10, 7) =  0
    k ( 11, 7) =  1.484484e+00
    q0( 11, 7) =  1.085765e+00
    n ( 11, 7) =  0
    k ( 12, 7) =  1.482304e+00
    q0( 12, 7) =  1.084994e+00
    n ( 12, 7) =  0
    k ( 13, 7) =  2.707859e-01
    q0( 13, 7) =  2.092975e+00
    n ( 13, 7) =  0
    k ( 14, 7) =  1.420778e-01
    q0( 14, 7) =  2.097110e+00
    n ( 14, 7) =  0
    k ( 15, 7) =  2.707854e-01
    q0( 15, 7) =  2.092975e+00
    n ( 15, 7) =  0
    k ( 16, 7) =  1.420776e-01
    q0( 16, 7) =  2.097114e+00
    n ( 16, 7) =  0
    k ( 17, 7) =  1.160228e-01
    q0( 17, 7) =  2.083147e+00
    n ( 17, 7) =  0
    k ( 18, 7) =  2.120183e-01
    q0( 18, 7) =  2.100583e+00
    n ( 18, 7) =  0
    k ( 19, 7) =  2.171192e-01
    q0( 19, 7) =  2.108362e+00
    n ( 19, 7) =  0
    k ( 20, 7) =  1.531627e-01
    q0( 20, 7) =  2.078968e+00
    n ( 20, 7) =  0
    k ( 21, 7) =  1.748085e-01
    q0( 21, 7) =  2.093643e+00
    n ( 21, 7) =  0
    k ( 22, 7) =  2.520321e-01
    q0( 22, 7) =  2.081950e+00
    n ( 22, 7) =  0
    k ( 23, 7) =  1.688261e-01
    q0( 23, 7) =  2.100868e+00
    n ( 23, 7) =  0
    k ( 24, 7) =  1.616118e-01
    q0( 24, 7) =  2.096324e+00
    n ( 24, 7) =  0
    k ( 25, 7) =  2.171190e-01
    q0( 25, 7) =  2.108362e+00
    n ( 25, 7) =  0
    k ( 26, 7) =  1.616120e-01
    q0( 26, 7) =  2.096327e+00
    n ( 26, 7) =  0
    k ( 27, 7) =  1.688260e-01
    q0( 27, 7) =  2.100869e+00
    n ( 27, 7) =  0
    k ( 28, 7) =  1.748085e-01
    q0( 28, 7) =  2.093646e+00
    n ( 28, 7) =  0
    k ( 29, 7) =  2.120188e-01
    q0( 29, 7) =  2.100580e+00
    n ( 29, 7) =  0
    k ( 30, 7) =  1.531624e-01
    q0( 30, 7) =  2.078966e+00
    n ( 30, 7) =  0
    k ( 31, 7) =  7.404584e-03
    q0( 31, 7) =  0.000000e+00
    n ( 31, 7) =  2
    k ( 32, 7) =  1.230381e-02
    q0( 32, 7) =  0.000000e+00
    n ( 32, 7) =  2
    k ( 33, 7) =  7.401157e-03
    q0( 33, 7) =  0.000000e+00
    n ( 33, 7) =  2
    k ( 34, 7) =  1.230155e-02
    q0( 34, 7) =  0.000000e+00
    n ( 34, 7) =  2
    k ( 35, 7) = -2.666588e-19
    q0( 35, 7) =  0.000000e+00
    n ( 35, 7) =  2
    k ( 36, 7) =  1.100249e-02
    q0( 36, 7) =  0.000000e+00
    n ( 36, 7) =  2
    k ( 37, 7) =  5.014385e-03
    q0( 37, 7) =  0.000000e+00
    n ( 37, 7) =  2
    k ( 38, 7) =  1.139501e-02
    q0( 38, 7) =  0.000000e+00
    n ( 38, 7) =  2
    k ( 39, 7) =  3.523411e-19
    q0( 39, 7) =  0.000000e+00
    n ( 39, 7) =  2
    k ( 40, 7) =  1.238713e-02
    q0( 40, 7) =  0.000000e+00
    n ( 40, 7) =  2
    k ( 41, 7) =  5.017983e-03
    q0( 41, 7) =  0.000000e+00
    n ( 41, 7) =  2
    k ( 42, 7) =  1.253947e-02
    q0( 42, 7) =  0.000000e+00
    n ( 42, 7) =  2
    k ( 43, 7) =  1.368269e-02
    q0( 43, 7) =  0.000000e+00
    n ( 43, 7) =  2
    k ( 44, 7) =  1.368366e-02
    q0( 44, 7) =  0.000000e+00
    n ( 44, 7) =  2
    k ( 45, 7) =  1.253903e-02
    q0( 45, 7) =  0.000000e+00
    n ( 45, 7) =  2
    k ( 46, 7) =  1.238893e-02
    q0( 46, 7) =  0.000000e+00
    n ( 46, 7) =  2
    k ( 47, 7) =  1.100671e-02
    q0( 47, 7) =  0.000000e+00
    n ( 47, 7) =  2
    k ( 48, 7) =  1.139524e-02
    q0( 48, 7) =  0.000000e+00
    n ( 48, 7) =  2
    k ( 49, 7) =  1.849697e-03
    q0( 49, 7) =  0.000000e+00
    n ( 49, 7) =  2
    k ( 50, 7) =  1.926491e-03
    q0( 50, 7) =  0.000000e+00
    n ( 50, 7) =  2
    k ( 51, 7) =  5.560487e-20
    q0( 51, 7) =  0.000000e+00
    n ( 51, 7) =  2
    k ( 52, 7) =  2.298999e-19
    q0( 52, 7) =  0.000000e+00
    n ( 52, 7) =  2
    k ( 53, 7) =  1.924750e-03
    q0( 53, 7) =  0.000000e+00
    n ( 53, 7) =  2
    k ( 54, 7) =  1.848317e-03
    q0( 54, 7) =  0.000000e+00
    n ( 54, 7) =  2
    k ( 55, 7) =  1.386778e-01
    q0( 55, 7) =  6.857471e-07
    n ( 55, 7) =  0
    k ( 56, 7) =  1.386744e-01
    q0( 56, 7) =  6.858292e-07
    n ( 56, 7) =  0
    k ( 57, 7) =  1.526267e-01
    q0( 57, 7) =  5.433153e-07
    n ( 57, 7) =  0
    k ( 58, 7) =  2.237168e-01
    q0( 58, 7) =  2.230362e-07
    n ( 58, 7) =  0
    k ( 59, 7) =  1.369988e-01
    q0( 59, 7) =  1.208118e-07
    n ( 59, 7) =  0
    k ( 60, 7) =  1.525977e-01
    q0( 60, 7) =  5.434899e-07
    n ( 60, 7) =  0
    !--- END generated by 'gen_FFparam.py'

  case(2)
    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 2, anchor point 1)
    k (  1, 1) =  2.079465e+00
    q0(  1, 1) =  1.410030e+00
    n (  1, 1) =  0
    k (  2, 1) =  1.961431e+00
    q0(  2, 1) =  1.415558e+00
    n (  2, 1) =  0
    k (  3, 1) =  2.502119e+00
    q0(  3, 1) =  1.723831e+00
    n (  3, 1) =  0
    k (  4, 1) =  2.463908e+00
    q0(  4, 1) =  1.425356e+00
    n (  4, 1) =  0
    k (  5, 1) =  1.548588e+00
    q0(  5, 1) =  1.081761e+00
    n (  5, 1) =  0
    k (  6, 1) =  2.280753e+00
    q0(  6, 1) =  1.407401e+00
    n (  6, 1) =  0
    k (  7, 1) =  1.550704e+00
    q0(  7, 1) =  1.079899e+00
    n (  7, 1) =  0
    k (  8, 1) =  2.307415e+00
    q0(  8, 1) =  1.407449e+00
    n (  8, 1) =  0
    k (  9, 1) =  1.505291e+00
    q0(  9, 1) =  1.084130e+00
    n (  9, 1) =  0
    k ( 10, 1) =  2.476788e+00
    q0( 10, 1) =  1.425103e+00
    n ( 10, 1) =  0
    k ( 11, 1) =  1.545158e+00
    q0( 11, 1) =  1.080544e+00
    n ( 11, 1) =  0
    k ( 12, 1) =  1.546923e+00
    q0( 12, 1) =  1.082036e+00
    n ( 12, 1) =  0
    k ( 13, 1) =  7.998842e-02
    q0( 13, 1) =  2.056193e+00
    n ( 13, 1) =  0
    k ( 14, 1) =  1.374451e-01
    q0( 14, 1) =  2.096172e+00
    n ( 14, 1) =  0
    k ( 15, 1) =  8.172076e-02
    q0( 15, 1) =  2.045544e+00
    n ( 15, 1) =  0
    k ( 16, 1) =  1.321691e-01
    q0( 16, 1) =  2.113155e+00
    n ( 16, 1) =  0
    k ( 17, 1) =  2.012576E-01
    q0( 17, 1) =  1.667300E+00
    n ( 17, 1) =  0
    k ( 18, 1) =  2.604370e-01
    q0( 18, 1) =  2.167751e+00
    n ( 18, 1) =  0
    k ( 19, 1) =  1.232835e-01
    q0( 19, 1) =  2.014283e+00
    n ( 19, 1) =  0
    k ( 20, 1) =  1.170009e-01
    q0( 20, 1) =  2.071240e+00
    n ( 20, 1) =  0
    k ( 21, 1) =  1.535793e-01
    q0( 21, 1) =  2.106718e+00
    n ( 21, 1) =  0
    k ( 22, 1) =  1.681110e-01
    q0( 22, 1) =  2.130551e+00
    n ( 22, 1) =  0
    k ( 23, 1) =  3.526958e-01
    q0( 23, 1) =  2.150170e+00
    n ( 23, 1) =  0
    k ( 24, 1) =  1.275305e-01
    q0( 24, 1) =  2.067309e+00
    n ( 24, 1) =  0
    k ( 25, 1) =  1.250607e-01
    q0( 25, 1) =  2.105310e+00
    n ( 25, 1) =  0
    k ( 26, 1) =  1.034276e-01
    q0( 26, 1) =  2.076792e+00
    n ( 26, 1) =  0
    k ( 27, 1) =  1.231069e-01
    q0( 27, 1) =  2.099850e+00
    n ( 27, 1) =  0
    k ( 28, 1) =  1.264939e-01
    q0( 28, 1) =  2.066174e+00
    n ( 28, 1) =  0
    k ( 29, 1) =  1.619420e-01
    q0( 29, 1) =  2.125946e+00
    n ( 29, 1) =  0
    k ( 30, 1) =  1.380703e-01
    q0( 30, 1) =  2.101279e+00
    n ( 30, 1) =  0
    k ( 31, 1) =  1.605235e-01
    q0( 31, 1) =  2.107561e+00
    n ( 31, 1) =  0
    k ( 32, 1) =  1.862040e-19
    q0( 32, 1) =  0.000000e+00
    n ( 32, 1) =  2
    k ( 33, 1) =  7.271756e-03
    q0( 33, 1) =  0.000000e+00
    n ( 33, 1) =  2
    k ( 34, 1) = -8.500540e-19
    q0( 34, 1) =  0.000000e+00
    n ( 34, 1) =  2
    k ( 35, 1) =  8.084170e-03
    q0( 35, 1) =  0.000000e+00
    n ( 35, 1) =  2
    k ( 36, 1) =  1.691922e-03
    q0( 36, 1) =  0.000000e+00
    n ( 36, 1) =  2
    k ( 37, 1) =  7.561672e-03
    q0( 37, 1) =  0.000000e+00
    n ( 37, 1) =  2
    k ( 38, 1) =  3.613234e-03
    q0( 38, 1) =  0.000000e+00
    n ( 38, 1) =  2
    k ( 39, 1) =  4.612244e-03
    q0( 39, 1) =  0.000000e+00
    n ( 39, 1) =  2
    k ( 40, 1) =  4.233461e-03
    q0( 40, 1) =  0.000000e+00
    n ( 40, 1) =  2
    k ( 41, 1) =  3.613224e-03
    q0( 41, 1) =  0.000000e+00
    n ( 41, 1) =  2
    k ( 42, 1) =  3.910194e-03
    q0( 42, 1) =  0.000000e+00
    n ( 42, 1) =  2
    k ( 43, 1) =  3.682482e-03
    q0( 43, 1) =  0.000000e+00
    n ( 43, 1) =  2
    k ( 44, 1) =  4.532228e-03
    q0( 44, 1) =  0.000000e+00
    n ( 44, 1) =  2
    k ( 45, 1) =  3.770345e-03
    q0( 45, 1) =  0.000000e+00
    n ( 45, 1) =  2
    k ( 46, 1) =  4.247783e-03
    q0( 46, 1) =  0.000000e+00
    n ( 46, 1) =  2
    k ( 47, 1) =  2.210487e-03
    q0( 47, 1) =  0.000000e+00
    n ( 47, 1) =  2
    k ( 48, 1) =  5.982627e-03
    q0( 48, 1) =  0.000000e+00
    n ( 48, 1) =  2
    k ( 49, 1) =  6.057412e-03
    q0( 49, 1) =  0.000000e+00
    n ( 49, 1) =  2
    k ( 50, 1) =  1.996855e-03
    q0( 50, 1) =  0.000000e+00
    n ( 50, 1) =  2
    k ( 51, 1) =  2.405180e-03
    q0( 51, 1) =  0.000000e+00
    n ( 51, 1) =  2
    k ( 52, 1) = -3.979061e-19
    q0( 52, 1) =  0.000000e+00
    n ( 52, 1) =  2
    k ( 53, 1) = -1.601790e-18
    q0( 53, 1) =  0.000000e+00
    n ( 53, 1) =  2
    k ( 54, 1) =  2.358012e-03
    q0( 54, 1) =  0.000000e+00
    n ( 54, 1) =  2
    k ( 55, 1) =  1.951237e-03
    q0( 55, 1) =  0.000000e+00
    n ( 55, 1) =  2
    k ( 56, 1) =  5.223792e-02
    q0( 56, 1) =  1.683841e-04
    n ( 56, 1) =  0
    k ( 57, 1) = -1.468285e-18
    q0( 57, 1) =  5.465064e-05
    n ( 57, 1) =  0
    k ( 58, 1) =  5.166171e-02
    q0( 58, 1) =  5.024368e-05
    n ( 58, 1) =  0
    k ( 59, 1) =  9.321288e-02
    q0( 59, 1) =  6.326185e-07
    n ( 59, 1) =  0
    k ( 60, 1) =  2.774174e-01
    q0( 60, 1) =  4.179489e-06
    n ( 60, 1) =  0
    k ( 61, 1) =  5.400849e-02
    q0( 61, 1) =  6.254581e-05
    n ( 61, 1) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 2, anchor point 2)
    k (  1, 2) =  2.042088e+00
    q0(  1, 2) =  1.412006e+00
    n (  1, 2) =  0
    k (  2, 2) =  1.968470e+00
    q0(  2, 2) =  1.413828e+00
    n (  2, 2) =  0
    k (  3, 2) =  2.506132e+00
    q0(  3, 2) =  1.723497e+00
    n (  3, 2) =  0
    k (  4, 2) =  2.455485e+00
    q0(  4, 2) =  1.426054e+00
    n (  4, 2) =  0
    k (  5, 2) =  1.547886e+00
    q0(  5, 2) =  1.081912e+00
    n (  5, 2) =  0
    k (  6, 2) =  2.280188e+00
    q0(  6, 2) =  1.406961e+00
    n (  6, 2) =  0
    k (  7, 2) =  1.550155e+00
    q0(  7, 2) =  1.080089e+00
    n (  7, 2) =  0
    k (  8, 2) =  2.298740e+00
    q0(  8, 2) =  1.407431e+00
    n (  8, 2) =  0
    k (  9, 2) =  1.505138e+00
    q0(  9, 2) =  1.084198e+00
    n (  9, 2) =  0
    k ( 10, 2) =  2.478304e+00
    q0( 10, 2) =  1.424958e+00
    n ( 10, 2) =  0
    k ( 11, 2) =  1.546508e+00
    q0( 11, 2) =  1.080395e+00
    n ( 11, 2) =  0
    k ( 12, 2) =  1.547404e+00
    q0( 12, 2) =  1.081936e+00
    n ( 12, 2) =  0
    k ( 13, 2) =  7.957448e-02
    q0( 13, 2) =  2.053887e+00
    n ( 13, 2) =  0
    k ( 14, 2) =  1.382335e-01
    q0( 14, 2) =  2.100706e+00
    n ( 14, 2) =  0
    k ( 15, 2) =  8.060930e-02
    q0( 15, 2) =  2.048386e+00
    n ( 15, 2) =  0
    k ( 16, 2) =  1.324924e-01
    q0( 16, 2) =  2.108602e+00
    n ( 16, 2) =  0
    k ( 17, 2) =  2.091043E-01
    q0( 17, 2) =  1.642840E+00
    n ( 17, 2) =  0
    k ( 18, 2) =  2.631325e-01
    q0( 18, 2) =  2.166850e+00
    n ( 18, 2) =  0
    k ( 19, 2) =  1.185920e-01
    q0( 19, 2) =  2.017152e+00
    n ( 19, 2) =  0
    k ( 20, 2) =  1.134174e-01
    q0( 20, 2) =  2.071194e+00
    n ( 20, 2) =  0
    k ( 21, 2) =  1.547549e-01
    q0( 21, 2) =  2.108216e+00
    n ( 21, 2) =  0
    k ( 22, 2) =  1.675222e-01
    q0( 22, 2) =  2.129155e+00
    n ( 22, 2) =  0
    k ( 23, 2) =  3.544694e-01
    q0( 23, 2) =  2.150396e+00
    n ( 23, 2) =  0
    k ( 24, 2) =  1.276069e-01
    q0( 24, 2) =  2.067759e+00
    n ( 24, 2) =  0
    k ( 25, 2) =  1.244972e-01
    q0( 25, 2) =  2.104266e+00
    n ( 25, 2) =  0
    k ( 26, 2) =  1.031344e-01
    q0( 26, 2) =  2.077022e+00
    n ( 26, 2) =  0
    k ( 27, 2) =  1.239739e-01
    q0( 27, 2) =  2.100530e+00
    n ( 27, 2) =  0
    k ( 28, 2) =  1.264439e-01
    q0( 28, 2) =  2.065518e+00
    n ( 28, 2) =  0
    k ( 29, 2) =  1.628226e-01
    q0( 29, 2) =  2.126635e+00
    n ( 29, 2) =  0
    k ( 30, 2) =  1.239716e-01
    q0( 30, 2) =  2.100265e+00
    n ( 30, 2) =  0
    k ( 31, 2) =  1.595697e-01
    q0( 31, 2) =  2.106058e+00
    n ( 31, 2) =  0
    k ( 32, 2) = -3.220431e-19
    q0( 32, 2) =  0.000000e+00
    n ( 32, 2) =  2
    k ( 33, 2) =  7.676458e-03
    q0( 33, 2) =  0.000000e+00
    n ( 33, 2) =  2
    k ( 34, 2) =  5.788259e-19
    q0( 34, 2) =  0.000000e+00
    n ( 34, 2) =  2
    k ( 35, 2) =  8.520357e-03
    q0( 35, 2) =  0.000000e+00
    n ( 35, 2) =  2
    k ( 36, 2) =  2.083841e-03
    q0( 36, 2) =  0.000000e+00
    n ( 36, 2) =  2
    k ( 37, 2) =  7.715170e-03
    q0( 37, 2) =  0.000000e+00
    n ( 37, 2) =  2
    k ( 38, 2) =  3.919692e-03
    q0( 38, 2) =  0.000000e+00
    n ( 38, 2) =  2
    k ( 39, 2) =  5.130564e-03
    q0( 39, 2) =  0.000000e+00
    n ( 39, 2) =  2
    k ( 40, 2) =  4.489144e-03
    q0( 40, 2) =  0.000000e+00
    n ( 40, 2) =  2
    k ( 41, 2) =  7.043624e-20
    q0( 41, 2) =  0.000000e+00
    n ( 41, 2) =  2
    k ( 42, 2) =  3.622416e-03
    q0( 42, 2) =  0.000000e+00
    n ( 42, 2) =  2
    k ( 43, 2) =  3.049532e-03
    q0( 43, 2) =  0.000000e+00
    n ( 43, 2) =  2
    k ( 44, 2) =  4.394818e-03
    q0( 44, 2) =  0.000000e+00
    n ( 44, 2) =  2
    k ( 45, 2) =  3.650972e-03
    q0( 45, 2) =  0.000000e+00
    n ( 45, 2) =  2
    k ( 46, 2) =  4.087373e-03
    q0( 46, 2) =  0.000000e+00
    n ( 46, 2) =  2
    k ( 47, 2) =  1.215218e-03
    q0( 47, 2) =  0.000000e+00
    n ( 47, 2) =  2
    k ( 48, 2) =  6.402733e-03
    q0( 48, 2) =  0.000000e+00
    n ( 48, 2) =  2
    k ( 49, 2) =  5.903727e-03
    q0( 49, 2) =  0.000000e+00
    n ( 49, 2) =  2
    k ( 50, 2) =  1.689358e-03
    q0( 50, 2) =  0.000000e+00
    n ( 50, 2) =  2
    k ( 51, 2) =  2.452769e-03
    q0( 51, 2) =  0.000000e+00
    n ( 51, 2) =  2
    k ( 52, 2) =  5.644965e-19
    q0( 52, 2) =  0.000000e+00
    n ( 52, 2) =  2
    k ( 53, 2) = -1.030843e-19
    q0( 53, 2) =  0.000000e+00
    n ( 53, 2) =  2
    k ( 54, 2) =  2.515619e-03
    q0( 54, 2) =  0.000000e+00
    n ( 54, 2) =  2
    k ( 55, 2) =  1.853507e-03
    q0( 55, 2) =  0.000000e+00
    n ( 55, 2) =  2
    k ( 56, 2) =  4.663277e-02
    q0( 56, 2) =  3.817474e-06
    n ( 56, 2) =  0
    k ( 57, 2) = -6.006364e-18
    q0( 57, 2) =  3.623691e-06
    n ( 57, 2) =  0
    k ( 58, 2) =  5.302189e-02
    q0( 58, 2) =  7.383191e-06
    n ( 58, 2) =  0
    k ( 59, 2) =  8.761828e-02
    q0( 59, 2) =  1.257897e-06
    n ( 59, 2) =  0
    k ( 60, 2) =  2.813269e-01
    q0( 60, 2) =  3.947519e-07
    n ( 60, 2) =  0
    k ( 61, 2) =  5.675662e-02
    q0( 61, 2) =  1.870055e-06
    n ( 61, 2) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 2, anchor point 3)
    k (  1, 3) =  1.931324e+00
    q0(  1, 3) =  1.416256e+00
    n (  1, 3) =  0
    k (  2, 3) =  1.971143e+00
    q0(  2, 3) =  1.410711e+00
    n (  2, 3) =  0
    k (  3, 3) =  2.550453e+00
    q0(  3, 3) =  1.720360e+00
    n (  3, 3) =  0
    k (  4, 3) =  2.442456e+00
    q0(  4, 3) =  1.426900e+00
    n (  4, 3) =  0
    k (  5, 3) =  1.548070e+00
    q0(  5, 3) =  1.082050e+00
    n (  5, 3) =  0
    k (  6, 3) =  2.283913e+00
    q0(  6, 3) =  1.405824e+00
    n (  6, 3) =  0
    k (  7, 3) =  1.549850e+00
    q0(  7, 3) =  1.080348e+00
    n (  7, 3) =  0
    k (  8, 3) =  2.272927e+00
    q0(  8, 3) =  1.407717e+00
    n (  8, 3) =  0
    k (  9, 3) =  1.505403e+00
    q0(  9, 3) =  1.084278e+00
    n (  9, 3) =  0
    k ( 10, 3) =  2.482900e+00
    q0( 10, 3) =  1.424760e+00
    n ( 10, 3) =  0
    k ( 11, 3) =  1.549028e+00
    q0( 11, 3) =  1.080165e+00
    n ( 11, 3) =  0
    k ( 12, 3) =  1.546858e+00
    q0( 12, 3) =  1.081961e+00
    n ( 12, 3) =  0
    k ( 13, 3) =  7.678447e-02
    q0( 13, 3) =  2.049429e+00
    n ( 13, 3) =  0
    k ( 14, 3) =  1.382654e-01
    q0( 14, 3) =  2.107266e+00
    n ( 14, 3) =  0
    k ( 15, 3) =  7.523201e-02
    q0( 15, 3) =  2.053327e+00
    n ( 15, 3) =  0
    k ( 16, 3) =  1.309652e-01
    q0( 16, 3) =  2.100633e+00
    n ( 16, 3) =  0
    k ( 17, 3) =  1.805942E-01
    q0( 17, 3) =  1.660180E+00
    n ( 17, 3) =  0
    k ( 18, 3) =  2.676916e-01
    q0( 18, 3) =  2.165462e+00
    n ( 18, 3) =  0
    k ( 19, 3) =  1.101378e-01
    q0( 19, 3) =  2.006201e+00
    n ( 19, 3) =  0
    k ( 20, 3) =  1.027079e-01
    q0( 20, 3) =  2.072904e+00
    n ( 20, 3) =  0
    k ( 21, 3) =  1.571628e-01
    q0( 21, 3) =  2.108466e+00
    n ( 21, 3) =  0
    k ( 22, 3) =  1.677307e-01
    q0( 22, 3) =  2.127816e+00
    n ( 22, 3) =  0
    k ( 23, 3) =  3.580636e-01
    q0( 23, 3) =  2.150844e+00
    n ( 23, 3) =  0
    k ( 24, 3) =  1.277508e-01
    q0( 24, 3) =  2.067846e+00
    n ( 24, 3) =  0
    k ( 25, 3) =  1.234885e-01
    q0( 25, 3) =  2.102713e+00
    n ( 25, 3) =  0
    k ( 26, 3) =  1.045249e-01
    q0( 26, 3) =  2.075778e+00
    n ( 26, 3) =  0
    k ( 27, 3) =  1.256801e-01
    q0( 27, 3) =  2.102155e+00
    n ( 27, 3) =  0
    k ( 28, 3) =  1.262500e-01
    q0( 28, 3) =  2.065019e+00
    n ( 28, 3) =  0
    k ( 29, 3) =  1.679536e-01
    q0( 29, 3) =  2.128586e+00
    n ( 29, 3) =  0
    k ( 30, 3) =  9.409439e-02
    q0( 30, 3) =  2.113073e+00
    n ( 30, 3) =  0
    k ( 31, 3) =  1.564671e-01
    q0( 31, 3) =  2.105033e+00
    n ( 31, 3) =  0
    k ( 32, 3) = -3.035246e-18
    q0( 32, 3) =  0.000000e+00
    n ( 32, 3) =  2
    k ( 33, 3) =  4.725086e-03
    q0( 33, 3) =  0.000000e+00
    n ( 33, 3) =  2
    k ( 34, 3) =  1.229260e-18
    q0( 34, 3) =  0.000000e+00
    n ( 34, 3) =  2
    k ( 35, 3) =  4.599719e-03
    q0( 35, 3) =  0.000000e+00
    n ( 35, 3) =  2
    k ( 36, 3) =  1.271177e-03
    q0( 36, 3) =  0.000000e+00
    n ( 36, 3) =  2
    k ( 37, 3) =  1.470479e-03
    q0( 37, 3) =  0.000000e+00
    n ( 37, 3) =  2
    k ( 38, 3) =  6.526457e-03
    q0( 38, 3) =  0.000000e+00
    n ( 38, 3) =  2
    k ( 39, 3) =  6.796678e-03
    q0( 39, 3) =  0.000000e+00
    n ( 39, 3) =  2
    k ( 40, 3) =  8.765119e-19
    q0( 40, 3) =  0.000000e+00
    n ( 40, 3) =  2
    k ( 41, 3) =  3.214762e-03
    q0( 41, 3) =  0.000000e+00
    n ( 41, 3) =  2
    k ( 42, 3) =  3.964372e-03
    q0( 42, 3) =  0.000000e+00
    n ( 42, 3) =  2
    k ( 43, 3) =  6.345482e-03
    q0( 43, 3) =  0.000000e+00
    n ( 43, 3) =  2
    k ( 44, 3) =  5.503521e-03
    q0( 44, 3) =  0.000000e+00
    n ( 44, 3) =  2
    k ( 45, 3) =  6.470799e-03
    q0( 45, 3) =  0.000000e+00
    n ( 45, 3) =  2
    k ( 46, 3) =  5.575892e-03
    q0( 46, 3) =  0.000000e+00
    n ( 46, 3) =  2
    k ( 47, 3) =  1.383535e-02
    q0( 47, 3) =  0.000000e+00
    n ( 47, 3) =  2
    k ( 48, 3) =  2.657937e-03
    q0( 48, 3) =  0.000000e+00
    n ( 48, 3) =  2
    k ( 49, 3) =  4.659333e-03
    q0( 49, 3) =  0.000000e+00
    n ( 49, 3) =  2
    k ( 50, 3) =  3.732391e-03
    q0( 50, 3) =  0.000000e+00
    n ( 50, 3) =  2
    k ( 51, 3) =  2.800938e-03
    q0( 51, 3) =  0.000000e+00
    n ( 51, 3) =  2
    k ( 52, 3) = -6.148079e-19
    q0( 52, 3) =  0.000000e+00
    n ( 52, 3) =  2
    k ( 53, 3) =  5.182819e-19
    q0( 53, 3) =  0.000000e+00
    n ( 53, 3) =  2
    k ( 54, 3) =  3.413332e-03
    q0( 54, 3) =  0.000000e+00
    n ( 54, 3) =  2
    k ( 55, 3) =  3.493926e-03
    q0( 55, 3) =  0.000000e+00
    n ( 55, 3) =  2
    k ( 56, 3) =  2.827930e-02
    q0( 56, 3) =  9.137295e-06
    n ( 56, 3) =  0
    k ( 57, 3) =  1.161135e-03
    q0( 57, 3) =  1.290043e-04
    n ( 57, 3) =  0
    k ( 58, 3) =  3.059245e-02
    q0( 58, 3) =  3.146994e-05
    n ( 58, 3) =  0
    k ( 59, 3) = -2.251376e-18
    q0( 59, 3) =  8.978456e-06
    n ( 59, 3) =  0
    k ( 60, 3) =  2.334060e-01
    q0( 60, 3) =  4.608660e-06
    n ( 60, 3) =  0
    k ( 61, 3) =  2.072143e-02
    q0( 61, 3) =  1.773691e-05
    n ( 61, 3) =  0
    !--- END generated by 'gen_FFparam.py'

  case(3)
    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 3, anchor point 1)
    k (  1, 1) =  2.310558e+00
    q0(  1, 1) =  1.400127e+00
    n (  1, 1) =  0
    k (  2, 1) =  2.083109e+00
    q0(  2, 1) =  1.402304e+00
    n (  2, 1) =  0
    k (  3, 1) =  4.181930e-01
    q0(  3, 1) =  1.952841e+00
    n (  3, 1) =  0
    k (  4, 1) =  2.545152e+00
    q0(  4, 1) =  1.380005e+00
    n (  4, 1) =  0
    k (  5, 1) =  1.416220e+00
    q0(  5, 1) =  1.092025e+00
    n (  5, 1) =  0
    k (  6, 1) =  2.279895e+00
    q0(  6, 1) =  1.398043e+00
    n (  6, 1) =  0
    k (  7, 1) =  1.523799e+00
    q0(  7, 1) =  1.083322e+00
    n (  7, 1) =  0
    k (  8, 1) =  2.271562e+00
    q0(  8, 1) =  1.399392e+00
    n (  8, 1) =  0
    k (  9, 1) =  1.528383e+00
    q0(  9, 1) =  1.083964e+00
    n (  9, 1) =  0
    k ( 10, 1) =  2.540767e+00
    q0( 10, 1) =  1.379427e+00
    n ( 10, 1) =  0
    k ( 11, 1) =  1.506135e+00
    q0( 11, 1) =  1.085019e+00
    n ( 11, 1) =  0
    k ( 12, 1) =  1.351096e+00
    q0( 12, 1) =  1.096719e+00
    n ( 12, 1) =  0
    k ( 13, 1) =  1.932464e-01
    q0( 13, 1) =  2.096564e+00
    n ( 13, 1) =  0
    k ( 14, 1) =  7.345143e-02
    q0( 14, 1) =  2.031994e+00
    n ( 14, 1) =  0
    k ( 15, 1) =  2.847566e-01
    q0( 15, 1) =  2.094612e+00
    n ( 15, 1) =  0
    k ( 16, 1) =  1.149791e-01
    q0( 16, 1) =  2.054076e+00
    n ( 16, 1) =  0
    k ( 17, 1) =  1.551681E-01
    q0( 17, 1) =  1.512550E+00
    n ( 17, 1) =  0
    k ( 18, 1) =  1.274360e-01
    q0( 18, 1) =  2.097687e+00
    n ( 18, 1) =  0
    k ( 19, 1) =  3.704415e-01
    q0( 19, 1) =  2.027580e+00
    n ( 19, 1) =  0
    k ( 20, 1) =  2.325325e-01
    q0( 20, 1) =  2.086893e+00
    n ( 20, 1) =  0
    k ( 21, 1) =  1.542787e-01
    q0( 21, 1) =  2.102575e+00
    n ( 21, 1) =  0
    k ( 22, 1) =  2.209864e-01
    q0( 22, 1) =  2.154009e+00
    n ( 22, 1) =  0
    k ( 23, 1) =  2.440567e-01
    q0( 23, 1) =  2.108652e+00
    n ( 23, 1) =  0
    k ( 24, 1) =  1.751891e-01
    q0( 24, 1) =  2.087313e+00
    n ( 24, 1) =  0
    k ( 25, 1) =  1.578772e-01
    q0( 25, 1) =  2.093774e+00
    n ( 25, 1) =  0
    k ( 26, 1) =  1.971060e-01
    q0( 26, 1) =  2.083343e+00
    n ( 26, 1) =  0
    k ( 27, 1) =  1.647239e-01
    q0( 27, 1) =  2.089744e+00
    n ( 27, 1) =  0
    k ( 28, 1) =  1.838054e-01
    q0( 28, 1) =  2.087561e+00
    n ( 28, 1) =  0
    k ( 29, 1) =  2.086744e-01
    q0( 29, 1) =  2.137323e+00
    n ( 29, 1) =  0
    k ( 30, 1) =  1.987688e-01
    q0( 30, 1) =  2.154926e+00
    n ( 30, 1) =  0
    k ( 31, 1) =  1.384887e-01
    q0( 31, 1) =  2.112004e+00
    n ( 31, 1) =  0
    k ( 32, 1) =  8.086088e-03
    q0( 32, 1) =  0.000000e+00
    n ( 32, 1) =  2
    k ( 33, 1) =  1.403272e-02
    q0( 33, 1) =  0.000000e+00
    n ( 33, 1) =  2
    k ( 34, 1) =  7.987513e-03
    q0( 34, 1) =  0.000000e+00
    n ( 34, 1) =  2
    k ( 35, 1) =  1.962870e-02
    q0( 35, 1) =  0.000000e+00
    n ( 35, 1) =  2
    k ( 36, 1) =  1.016497e-03
    q0( 36, 1) =  0.000000e+00
    n ( 36, 1) =  2
    k ( 37, 1) =  2.379387e-02
    q0( 37, 1) =  0.000000e+00
    n ( 37, 1) =  2
    k ( 38, 1) =  6.685627e-18
    q0( 38, 1) =  0.000000e+00
    n ( 38, 1) =  2
    k ( 39, 1) =  9.657400e-03
    q0( 39, 1) =  0.000000e+00
    n ( 39, 1) =  2
    k ( 40, 1) =  1.206445e-03
    q0( 40, 1) =  0.000000e+00
    n ( 40, 1) =  2
    k ( 41, 1) =  1.293066e-17
    q0( 41, 1) =  0.000000e+00
    n ( 41, 1) =  2
    k ( 42, 1) = -6.840716e-18
    q0( 42, 1) =  0.000000e+00
    n ( 42, 1) =  2
    k ( 43, 1) =  5.414778e-03
    q0( 43, 1) =  0.000000e+00
    n ( 43, 1) =  2
    k ( 44, 1) =  1.380757e-02
    q0( 44, 1) =  0.000000e+00
    n ( 44, 1) =  2
    k ( 45, 1) =  2.958534e-03
    q0( 45, 1) =  0.000000e+00
    n ( 45, 1) =  2
    k ( 46, 1) =  1.171033e-02
    q0( 46, 1) =  0.000000e+00
    n ( 46, 1) =  2
    k ( 47, 1) =  7.194867e-18
    q0( 47, 1) =  0.000000e+00
    n ( 47, 1) =  2
    k ( 48, 1) =  4.590395e-03
    q0( 48, 1) =  0.000000e+00
    n ( 48, 1) =  2
    k ( 49, 1) =  1.340472e-02
    q0( 49, 1) =  0.000000e+00
    n ( 49, 1) =  2
    k ( 50, 1) = -1.825963e-18
    q0( 50, 1) =  0.000000e+00
    n ( 50, 1) =  2
    k ( 51, 1) =  8.785665e-18
    q0( 51, 1) =  0.000000e+00
    n ( 51, 1) =  2
    k ( 52, 1) = -4.023449e-17
    q0( 52, 1) =  0.000000e+00
    n ( 52, 1) =  2
    k ( 53, 1) =  3.837148e-19
    q0( 53, 1) =  0.000000e+00
    n ( 53, 1) =  2
    k ( 54, 1) = -2.353001e-18
    q0( 54, 1) =  0.000000e+00
    n ( 54, 1) =  2
    k ( 55, 1) =  2.590805e-18
    q0( 55, 1) =  0.000000e+00
    n ( 55, 1) =  2
    k ( 56, 1) =  2.140563e-01
    q0( 56, 1) =  6.057128e-05
    n ( 56, 1) =  0
    k ( 57, 1) =  1.501968e-01
    q0( 57, 1) =  4.623978e-05
    n ( 57, 1) =  0
    k ( 58, 1) =  2.290009e-01
    q0( 58, 1) =  1.957332e-05
    n ( 58, 1) =  0
    k ( 59, 1) = -6.318513e-17
    q0( 59, 1) =  3.780271e-05
    n ( 59, 1) =  0
    k ( 60, 1) =  2.936814e-01
    q0( 60, 1) =  3.417183e-06
    n ( 60, 1) =  0
    k ( 61, 1) =  2.272500e-01
    q0( 61, 1) =  2.176340e-05
    n ( 61, 1) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 3, anchor point 2)
    k (  1, 2) =  2.281158e+00
    q0(  1, 2) =  1.403447e+00
    n (  1, 2) =  0
    k (  2, 2) =  2.056345e+00
    q0(  2, 2) =  1.402839e+00
    n (  2, 2) =  0
    k (  3, 2) =  5.701203e-01
    q0(  3, 2) =  1.843035e+00
    n (  3, 2) =  0
    k (  4, 2) =  2.542655e+00
    q0(  4, 2) =  1.380568e+00
    n (  4, 2) =  0
    k (  5, 2) =  1.461359e+00
    q0(  5, 2) =  1.088796e+00
    n (  5, 2) =  0
    k (  6, 2) =  2.283847e+00
    q0(  6, 2) =  1.397692e+00
    n (  6, 2) =  0
    k (  7, 2) =  1.531471e+00
    q0(  7, 2) =  1.082901e+00
    n (  7, 2) =  0
    k (  8, 2) =  2.282535e+00
    q0(  8, 2) =  1.398211e+00
    n (  8, 2) =  0
    k (  9, 2) =  1.531152e+00
    q0(  9, 2) =  1.083797e+00
    n (  9, 2) =  0
    k ( 10, 2) =  2.535499e+00
    q0( 10, 2) =  1.380150e+00
    n ( 10, 2) =  0
    k ( 11, 2) =  1.516598e+00
    q0( 11, 2) =  1.084325e+00
    n ( 11, 2) =  0
    k ( 12, 2) =  1.378970e+00
    q0( 12, 2) =  1.094329e+00
    n ( 12, 2) =  0
    k ( 13, 2) =  1.966876e-01
    q0( 13, 2) =  2.094748e+00
    n ( 13, 2) =  0
    k ( 14, 2) =  7.096395e-02
    q0( 14, 2) =  2.041652e+00
    n ( 14, 2) =  0
    k ( 15, 2) =  2.899411e-01
    q0( 15, 2) =  2.093385e+00
    n ( 15, 2) =  0
    k ( 16, 2) =  1.109109e-01
    q0( 16, 2) =  2.049448e+00
    n ( 16, 2) =  0
    k ( 17, 2) =  1.359533E-01
    q0( 17, 2) =  1.497760E+00
    n ( 17, 2) =  0
    k ( 18, 2) =  1.159916e-01
    q0( 18, 2) =  2.095989e+00
    n ( 18, 2) =  0
    k ( 19, 2) =  3.885487e-01
    q0( 19, 2) =  2.032920e+00
    n ( 19, 2) =  0
    k ( 20, 2) =  2.350427e-01
    q0( 20, 2) =  2.088762e+00
    n ( 20, 2) =  0
    k ( 21, 2) =  1.552872e-01
    q0( 21, 2) =  2.100367e+00
    n ( 21, 2) =  0
    k ( 22, 2) =  2.214173e-01
    q0( 22, 2) =  2.146484e+00
    n ( 22, 2) =  0
    k ( 23, 2) =  2.440212e-01
    q0( 23, 2) =  2.106902e+00
    n ( 23, 2) =  0
    k ( 24, 2) =  1.750206e-01
    q0( 24, 2) =  2.087989e+00
    n ( 24, 2) =  0
    k ( 25, 2) =  1.587225e-01
    q0( 25, 2) =  2.094174e+00
    n ( 25, 2) =  0
    k ( 26, 2) =  1.988587e-01
    q0( 26, 2) =  2.087780e+00
    n ( 26, 2) =  0
    k ( 27, 2) =  1.663275e-01
    q0( 27, 2) =  2.089934e+00
    n ( 27, 2) =  0
    k ( 28, 2) =  1.842158e-01
    q0( 28, 2) =  2.088676e+00
    n ( 28, 2) =  0
    k ( 29, 2) =  2.157657e-01
    q0( 29, 2) =  2.142272e+00
    n ( 29, 2) =  0
    k ( 30, 2) =  2.275259e-01
    q0( 30, 2) =  2.153873e+00
    n ( 30, 2) =  0
    k ( 31, 2) =  1.349261e-01
    q0( 31, 2) =  2.106674e+00
    n ( 31, 2) =  0
    k ( 32, 2) =  6.957902e-03
    q0( 32, 2) =  0.000000e+00
    n ( 32, 2) =  2
    k ( 33, 2) =  1.545732e-02
    q0( 33, 2) =  0.000000e+00
    n ( 33, 2) =  2
    k ( 34, 2) =  5.356156e-03
    q0( 34, 2) =  0.000000e+00
    n ( 34, 2) =  2
    k ( 35, 2) =  1.537744e-02
    q0( 35, 2) =  0.000000e+00
    n ( 35, 2) =  2
    k ( 36, 2) =  9.594041e-19
    q0( 36, 2) =  0.000000e+00
    n ( 36, 2) =  2
    k ( 37, 2) =  1.843538e-02
    q0( 37, 2) =  0.000000e+00
    n ( 37, 2) =  2
    k ( 38, 2) = -7.566669e-19
    q0( 38, 2) =  0.000000e+00
    n ( 38, 2) =  2
    k ( 39, 2) =  1.053734e-02
    q0( 39, 2) =  0.000000e+00
    n ( 39, 2) =  2
    k ( 40, 2) =  9.560561e-19
    q0( 40, 2) =  0.000000e+00
    n ( 40, 2) =  2
    k ( 41, 2) =  1.750327e-03
    q0( 41, 2) =  0.000000e+00
    n ( 41, 2) =  2
    k ( 42, 2) =  7.374537e-19
    q0( 42, 2) =  0.000000e+00
    n ( 42, 2) =  2
    k ( 43, 2) =  9.235224e-03
    q0( 43, 2) =  0.000000e+00
    n ( 43, 2) =  2
    k ( 44, 2) =  1.208511e-02
    q0( 44, 2) =  0.000000e+00
    n ( 44, 2) =  2
    k ( 45, 2) =  1.127996e-02
    q0( 45, 2) =  0.000000e+00
    n ( 45, 2) =  2
    k ( 46, 2) =  9.574813e-03
    q0( 46, 2) =  0.000000e+00
    n ( 46, 2) =  2
    k ( 47, 2) =  8.644190e-03
    q0( 47, 2) =  0.000000e+00
    n ( 47, 2) =  2
    k ( 48, 2) =  1.525221e-02
    q0( 48, 2) =  0.000000e+00
    n ( 48, 2) =  2
    k ( 49, 2) =  9.375076e-03
    q0( 49, 2) =  0.000000e+00
    n ( 49, 2) =  2
    k ( 50, 2) =  1.839942e-19
    q0( 50, 2) =  0.000000e+00
    n ( 50, 2) =  2
    k ( 51, 2) =  2.224012e-03
    q0( 51, 2) =  0.000000e+00
    n ( 51, 2) =  2
    k ( 52, 2) =  2.428539e-18
    q0( 52, 2) =  0.000000e+00
    n ( 52, 2) =  2
    k ( 53, 2) =  1.017065e-18
    q0( 53, 2) =  0.000000e+00
    n ( 53, 2) =  2
    k ( 54, 2) =  2.479856e-03
    q0( 54, 2) =  0.000000e+00
    n ( 54, 2) =  2
    k ( 55, 2) =  1.181751e-19
    q0( 55, 2) =  0.000000e+00
    n ( 55, 2) =  2
    k ( 56, 2) =  1.928861e-01
    q0( 56, 2) =  3.634476e-05
    n ( 56, 2) =  0
    k ( 57, 2) =  1.658853e-01
    q0( 57, 2) =  3.668086e-05
    n ( 57, 2) =  0
    k ( 58, 2) =  2.145445e-01
    q0( 58, 2) =  1.682311e-05
    n ( 58, 2) =  0
    k ( 59, 2) =  8.533899e-02
    q0( 59, 2) =  1.144445e-05
    n ( 59, 2) =  0
    k ( 60, 2) =  2.993355e-01
    q0( 60, 2) =  1.374154e-05
    n ( 60, 2) =  0
    k ( 61, 2) =  2.140906e-01
    q0( 61, 2) =  1.272600e-05
    n ( 61, 2) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 3, anchor point 3)
    k (  1, 3) =  2.164063e+00
    q0(  1, 3) =  1.414867e+00
    n (  1, 3) =  0
    k (  2, 3) =  2.025165e+00
    q0(  2, 3) =  1.413132e+00
    n (  2, 3) =  0
    k (  3, 3) =  1.883090e+00
    q0(  3, 3) =  1.725492e+00
    n (  3, 3) =  0
    k (  4, 3) =  2.565976e+00
    q0(  4, 3) =  1.378936e+00
    n (  4, 3) =  0
    k (  5, 3) =  1.538252e+00
    q0(  5, 3) =  1.083572e+00
    n (  5, 3) =  0
    k (  6, 3) =  2.282997e+00
    q0(  6, 3) =  1.397586e+00
    n (  6, 3) =  0
    k (  7, 3) =  1.539572e+00
    q0(  7, 3) =  1.082296e+00
    n (  7, 3) =  0
    k (  8, 3) =  2.304633e+00
    q0(  8, 3) =  1.395793e+00
    n (  8, 3) =  0
    k (  9, 3) =  1.538583e+00
    q0(  9, 3) =  1.083234e+00
    n (  9, 3) =  0
    k ( 10, 3) =  2.531058e+00
    q0( 10, 3) =  1.381013e+00
    n ( 10, 3) =  0
    k ( 11, 3) =  1.532520e+00
    q0( 11, 3) =  1.083067e+00
    n ( 11, 3) =  0
    k ( 12, 3) =  1.491946e+00
    q0( 12, 3) =  1.085545e+00
    n ( 12, 3) =  0
    k ( 13, 3) =  2.323831e-01
    q0( 13, 3) =  2.103615e+00
    n ( 13, 3) =  0
    k ( 14, 3) =  9.102823e-02
    q0( 14, 3) =  2.055045e+00
    n ( 14, 3) =  0
    k ( 15, 3) =  2.644463e-01
    q0( 15, 3) =  2.095355e+00
    n ( 15, 3) =  0
    k ( 16, 3) =  1.004454e-01
    q0( 16, 3) =  2.058990e+00
    n ( 16, 3) =  0
    k ( 17, 3) =  6.450642E-02
    q0( 17, 3) =  1.505040E+00
    n ( 17, 3) =  0
    k ( 18, 3) =  8.235398e-02
    q0( 18, 3) =  2.076498e+00
    n ( 18, 3) =  0
    k ( 19, 3) =  3.144890e-01
    q0( 19, 3) =  2.043754e+00
    n ( 19, 3) =  0
    k ( 20, 3) =  2.353738e-01
    q0( 20, 3) =  2.091722e+00
    n ( 20, 3) =  0
    k ( 21, 3) =  1.554502e-01
    q0( 21, 3) =  2.095545e+00
    n ( 21, 3) =  0
    k ( 22, 3) =  2.086111e-01
    q0( 22, 3) =  2.124364e+00
    n ( 22, 3) =  0
    k ( 23, 3) =  2.397573e-01
    q0( 23, 3) =  2.101920e+00
    n ( 23, 3) =  0
    k ( 24, 3) =  1.774260e-01
    q0( 24, 3) =  2.090055e+00
    n ( 24, 3) =  0
    k ( 25, 3) =  1.612598e-01
    q0( 25, 3) =  2.096011e+00
    n ( 25, 3) =  0
    k ( 26, 3) =  2.205481e-01
    q0( 26, 3) =  2.098305e+00
    n ( 26, 3) =  0
    k ( 27, 3) =  1.642046e-01
    q0( 27, 3) =  2.091015e+00
    n ( 27, 3) =  0
    k ( 28, 3) =  1.825681e-01
    q0( 28, 3) =  2.091589e+00
    n ( 28, 3) =  0
    k ( 29, 3) =  2.165453e-01
    q0( 29, 3) =  2.130047e+00
    n ( 29, 3) =  0
    k ( 30, 3) =  3.411760e-01
    q0( 30, 3) =  2.163817e+00
    n ( 30, 3) =  0
    k ( 31, 3) =  1.427006e-01
    q0( 31, 3) =  2.094519e+00
    n ( 31, 3) =  0
    k ( 32, 3) =  5.253782e-03
    q0( 32, 3) =  0.000000e+00
    n ( 32, 3) =  2
    k ( 33, 3) =  1.408246e-02
    q0( 33, 3) =  0.000000e+00
    n ( 33, 3) =  2
    k ( 34, 3) =  4.321567e-03
    q0( 34, 3) =  0.000000e+00
    n ( 34, 3) =  2
    k ( 35, 3) =  1.441817e-02
    q0( 35, 3) =  0.000000e+00
    n ( 35, 3) =  2
    k ( 36, 3) = -4.325454e-18
    q0( 36, 3) =  0.000000e+00
    n ( 36, 3) =  2
    k ( 37, 3) =  1.675533e-02
    q0( 37, 3) =  0.000000e+00
    n ( 37, 3) =  2
    k ( 38, 3) = -5.702422e-19
    q0( 38, 3) =  0.000000e+00
    n ( 38, 3) =  2
    k ( 39, 3) =  9.966848e-03
    q0( 39, 3) =  0.000000e+00
    n ( 39, 3) =  2
    k ( 40, 3) =  1.852221e-18
    q0( 40, 3) =  0.000000e+00
    n ( 40, 3) =  2
    k ( 41, 3) =  1.278432e-02
    q0( 41, 3) =  0.000000e+00
    n ( 41, 3) =  2
    k ( 42, 3) =  3.667105e-19
    q0( 42, 3) =  0.000000e+00
    n ( 42, 3) =  2
    k ( 43, 3) =  1.062771e-02
    q0( 43, 3) =  0.000000e+00
    n ( 43, 3) =  2
    k ( 44, 3) =  1.245578e-02
    q0( 44, 3) =  0.000000e+00
    n ( 44, 3) =  2
    k ( 45, 3) =  1.196926e-02
    q0( 45, 3) =  0.000000e+00
    n ( 45, 3) =  2
    k ( 46, 3) =  1.093318e-02
    q0( 46, 3) =  0.000000e+00
    n ( 46, 3) =  2
    k ( 47, 3) =  1.043458e-02
    q0( 47, 3) =  0.000000e+00
    n ( 47, 3) =  2
    k ( 48, 3) =  1.521237e-02
    q0( 48, 3) =  0.000000e+00
    n ( 48, 3) =  2
    k ( 49, 3) =  9.266782e-03
    q0( 49, 3) =  0.000000e+00
    n ( 49, 3) =  2
    k ( 50, 3) =  3.392227e-04
    q0( 50, 3) =  0.000000e+00
    n ( 50, 3) =  2
    k ( 51, 3) =  3.009962e-03
    q0( 51, 3) =  0.000000e+00
    n ( 51, 3) =  2
    k ( 52, 3) =  7.394745e-04
    q0( 52, 3) =  0.000000e+00
    n ( 52, 3) =  2
    k ( 53, 3) =  5.924151e-04
    q0( 53, 3) =  0.000000e+00
    n ( 53, 3) =  2
    k ( 54, 3) =  3.000830e-03
    q0( 54, 3) =  0.000000e+00
    n ( 54, 3) =  2
    k ( 55, 3) =  6.894458e-04
    q0( 55, 3) =  0.000000e+00
    n ( 55, 3) =  2
    k ( 56, 3) =  1.702255e-01
    q0( 56, 3) =  1.127324e-05
    n ( 56, 3) =  0
    k ( 57, 3) =  1.846629e-01
    q0( 57, 3) =  3.660479e-06
    n ( 57, 3) =  0
    k ( 58, 3) =  1.886835e-01
    q0( 58, 3) =  1.322188e-05
    n ( 58, 3) =  0
    k ( 59, 3) =  1.721867e-01
    q0( 59, 3) =  4.538248e-06
    n ( 59, 3) =  0
    k ( 60, 3) =  2.860467e-01
    q0( 60, 3) =  1.179678e-06
    n ( 60, 3) =  0
    k ( 61, 3) =  2.011796e-01
    q0( 61, 3) =  1.254435e-05
    n ( 61, 3) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 3, anchor point 4)
    k (  1, 4) =  2.130042e+00
    q0(  1, 4) =  1.419388e+00
    n (  1, 4) =  0
    k (  2, 4) =  2.119768e+00
    q0(  2, 4) =  1.416988e+00
    n (  2, 4) =  0
    k (  3, 4) =  2.462817e+00
    q0(  3, 4) =  1.704584e+00
    n (  3, 4) =  0
    k (  4, 4) =  2.591727e+00
    q0(  4, 4) =  1.377194e+00
    n (  4, 4) =  0
    k (  5, 4) =  1.550691e+00
    q0(  5, 4) =  1.082830e+00
    n (  5, 4) =  0
    k (  6, 4) =  2.260290e+00
    q0(  6, 4) =  1.399474e+00
    n (  6, 4) =  0
    k (  7, 4) =  1.538672e+00
    q0(  7, 4) =  1.082388e+00
    n (  7, 4) =  0
    k (  8, 4) =  2.310472e+00
    q0(  8, 4) =  1.394863e+00
    n (  8, 4) =  0
    k (  9, 4) =  1.542121e+00
    q0(  9, 4) =  1.082920e+00
    n (  9, 4) =  0
    k ( 10, 4) =  2.522015e+00
    q0( 10, 4) =  1.381699e+00
    n ( 10, 4) =  0
    k ( 11, 4) =  1.539575e+00
    q0( 11, 4) =  1.082446e+00
    n ( 11, 4) =  0
    k ( 12, 4) =  1.548884e+00
    q0( 12, 4) =  1.081952e+00
    n ( 12, 4) =  0
    k ( 13, 4) =  2.505476e-01
    q0( 13, 4) =  2.106023e+00
    n ( 13, 4) =  0
    k ( 14, 4) =  1.088551e-01
    q0( 14, 4) =  2.062002e+00
    n ( 14, 4) =  0
    k ( 15, 4) =  2.722284e-01
    q0( 15, 4) =  2.098895e+00
    n ( 15, 4) =  0
    k ( 16, 4) =  1.149034e-01
    q0( 16, 4) =  2.065255e+00
    n ( 16, 4) =  0
    k ( 17, 4) =  4.942193E-02
    q0( 17, 4) =  1.755600E+00
    n ( 17, 4) =  0
    k ( 18, 4) =  6.662221e-02
    q0( 18, 4) =  2.069944e+00
    n ( 18, 4) =  0
    k ( 19, 4) =  2.957028e-01
    q0( 19, 4) =  2.077300e+00
    n ( 19, 4) =  0
    k ( 20, 4) =  2.358957e-01
    q0( 20, 4) =  2.091314e+00
    n ( 20, 4) =  0
    k ( 21, 4) =  1.555981e-01
    q0( 21, 4) =  2.098046e+00
    n ( 21, 4) =  0
    k ( 22, 4) =  2.002178e-01
    q0( 22, 4) =  2.115642e+00
    n ( 22, 4) =  0
    k ( 23, 4) =  2.377661e-01
    q0( 23, 4) =  2.103321e+00
    n ( 23, 4) =  0
    k ( 24, 4) =  1.790463e-01
    q0( 24, 4) =  2.089139e+00
    n ( 24, 4) =  0
    k ( 25, 4) =  1.617343e-01
    q0( 25, 4) =  2.094195e+00
    n ( 25, 4) =  0
    k ( 26, 4) =  2.304359e-01
    q0( 26, 4) =  2.097884e+00
    n ( 26, 4) =  0
    k ( 27, 4) =  1.635948e-01
    q0( 27, 4) =  2.092924e+00
    n ( 27, 4) =  0
    k ( 28, 4) =  1.807829e-01
    q0( 28, 4) =  2.091104e+00
    n ( 28, 4) =  0
    k ( 29, 4) =  1.991843e-01
    q0( 29, 4) =  2.119433e+00
    n ( 29, 4) =  0
    k ( 30, 4) =  3.073020e-01
    q0( 30, 4) =  2.136934e+00
    n ( 30, 4) =  0
    k ( 31, 4) =  1.473926e-01
    q0( 31, 4) =  2.092692e+00
    n ( 31, 4) =  0
    k ( 32, 4) =  6.114448e-03
    q0( 32, 4) =  0.000000e+00
    n ( 32, 4) =  2
    k ( 33, 4) =  1.416597e-02
    q0( 33, 4) =  0.000000e+00
    n ( 33, 4) =  2
    k ( 34, 4) =  5.361552e-03
    q0( 34, 4) =  0.000000e+00
    n ( 34, 4) =  2
    k ( 35, 4) =  1.459075e-02
    q0( 35, 4) =  0.000000e+00
    n ( 35, 4) =  2
    k ( 36, 4) =  9.322934e-20
    q0( 36, 4) =  0.000000e+00
    n ( 36, 4) =  2
    k ( 37, 4) =  1.761570e-02
    q0( 37, 4) =  0.000000e+00
    n ( 37, 4) =  2
    k ( 38, 4) = -9.841576e-19
    q0( 38, 4) =  0.000000e+00
    n ( 38, 4) =  2
    k ( 39, 4) =  1.043376e-02
    q0( 39, 4) =  0.000000e+00
    n ( 39, 4) =  2
    k ( 40, 4) =  2.261850e-18
    q0( 40, 4) =  0.000000e+00
    n ( 40, 4) =  2
    k ( 41, 4) =  1.707624e-02
    q0( 41, 4) =  0.000000e+00
    n ( 41, 4) =  2
    k ( 42, 4) =  1.539639e-18
    q0( 42, 4) =  0.000000e+00
    n ( 42, 4) =  2
    k ( 43, 4) =  1.126973e-02
    q0( 43, 4) =  0.000000e+00
    n ( 43, 4) =  2
    k ( 44, 4) =  1.294985e-02
    q0( 44, 4) =  0.000000e+00
    n ( 44, 4) =  2
    k ( 45, 4) =  1.146681e-02
    q0( 45, 4) =  0.000000e+00
    n ( 45, 4) =  2
    k ( 46, 4) =  1.144805e-02
    q0( 46, 4) =  0.000000e+00
    n ( 46, 4) =  2
    k ( 47, 4) =  8.236247e-03
    q0( 47, 4) =  0.000000e+00
    n ( 47, 4) =  2
    k ( 48, 4) =  1.417272e-02
    q0( 48, 4) =  0.000000e+00
    n ( 48, 4) =  2
    k ( 49, 4) =  9.781039e-03
    q0( 49, 4) =  0.000000e+00
    n ( 49, 4) =  2
    k ( 50, 4) =  6.783696e-04
    q0( 50, 4) =  0.000000e+00
    n ( 50, 4) =  2
    k ( 51, 4) =  2.947367e-03
    q0( 51, 4) =  0.000000e+00
    n ( 51, 4) =  2
    k ( 52, 4) =  1.767069e-03
    q0( 52, 4) =  0.000000e+00
    n ( 52, 4) =  2
    k ( 53, 4) =  8.142323e-04
    q0( 53, 4) =  0.000000e+00
    n ( 53, 4) =  2
    k ( 54, 4) =  3.114179e-03
    q0( 54, 4) =  0.000000e+00
    n ( 54, 4) =  2
    k ( 55, 4) =  1.090924e-03
    q0( 55, 4) =  0.000000e+00
    n ( 55, 4) =  2
    k ( 56, 4) =  1.628889e-01
    q0( 56, 4) =  6.885659e-07
    n ( 56, 4) =  0
    k ( 57, 4) =  1.782787e-01
    q0( 57, 4) =  8.705742e-07
    n ( 57, 4) =  0
    k ( 58, 4) =  1.727729e-01
    q0( 58, 4) =  2.370261e-07
    n ( 58, 4) =  0
    k ( 59, 4) =  8.864287e-02
    q0( 59, 4) =  1.965272e-08
    n ( 59, 4) =  0
    k ( 60, 4) =  2.714496e-01
    q0( 60, 4) =  2.797603e-07
    n ( 60, 4) =  0
    k ( 61, 4) =  1.898144e-01
    q0( 61, 4) =  7.306042e-07
    n ( 61, 4) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 3, anchor point 5)
    k (  1, 5) =  2.350003e+00
    q0(  1, 5) =  1.420972e+00
    n (  1, 5) =  0
    k (  2, 5) =  2.694068e+00
    q0(  2, 5) =  1.421237e+00
    n (  2, 5) =  0
    k (  3, 5) =  4.475616e+00
    q0(  3, 5) =  1.694030e+00
    n (  3, 5) =  0
    k (  4, 5) =  2.561999e+00
    q0(  4, 5) =  1.378988e+00
    n (  4, 5) =  0
    k (  5, 5) =  1.520409e+00
    q0(  5, 5) =  1.082426e+00
    n (  5, 5) =  0
    k (  6, 5) =  2.319626e+00
    q0(  6, 5) =  1.398069e+00
    n (  6, 5) =  0
    k (  7, 5) =  1.544240e+00
    q0(  7, 5) =  1.082376e+00
    n (  7, 5) =  0
    k (  8, 5) =  2.313492e+00
    q0(  8, 5) =  1.397505e+00
    n (  8, 5) =  0
    k (  9, 5) =  1.537893e+00
    q0(  9, 5) =  1.083394e+00
    n (  9, 5) =  0
    k ( 10, 5) =  2.530460e+00
    q0( 10, 5) =  1.379689e+00
    n ( 10, 5) =  0
    k ( 11, 5) =  1.544508e+00
    q0( 11, 5) =  1.082436e+00
    n ( 11, 5) =  0
    k ( 12, 5) =  1.557013e+00
    q0( 12, 5) =  1.082434e+00
    n ( 12, 5) =  0
    k ( 13, 5) =  3.817523e-01
    q0( 13, 5) =  2.110381e+00
    n ( 13, 5) =  0
    k ( 14, 5) =  3.278262e-01
    q0( 14, 5) =  2.045842e+00
    n ( 14, 5) =  0
    k ( 15, 5) =  1.619503e-01
    q0( 15, 5) =  2.107670e+00
    n ( 15, 5) =  0
    k ( 16, 5) =  1.605017e-01
    q0( 16, 5) =  2.055565e+00
    n ( 16, 5) =  0
    k ( 17, 5) =  3.986381E-02
    q0( 17, 5) =  1.867290E+00
    n ( 17, 5) =  0
    k ( 18, 5) =  2.077577e-01
    q0( 18, 5) =  2.062934e+00
    n ( 18, 5) =  0
    k ( 19, 5) =  6.811215e-18
    q0( 19, 5) =  2.109265e+00
    n ( 19, 5) =  0
    k ( 20, 5) =  2.302737e-01
    q0( 20, 5) =  2.085202e+00
    n ( 20, 5) =  0
    k ( 21, 5) =  1.546808e-01
    q0( 21, 5) =  2.104534e+00
    n ( 21, 5) =  0
    k ( 22, 5) =  1.038466e-01
    q0( 22, 5) =  2.127224e+00
    n ( 22, 5) =  0
    k ( 23, 5) =  2.335461e-01
    q0( 23, 5) =  2.114309e+00
    n ( 23, 5) =  0
    k ( 24, 5) =  1.818217e-01
    q0( 24, 5) =  2.084691e+00
    n ( 24, 5) =  0
    k ( 25, 5) =  1.658523e-01
    q0( 25, 5) =  2.093728e+00
    n ( 25, 5) =  0
    k ( 26, 5) =  2.765262e-01
    q0( 26, 5) =  2.087077e+00
    n ( 26, 5) =  0
    k ( 27, 5) =  1.537724e-01
    q0( 27, 5) =  2.092978e+00
    n ( 27, 5) =  0
    k ( 28, 5) =  1.700442e-01
    q0( 28, 5) =  2.084804e+00
    n ( 28, 5) =  0
    k ( 29, 5) =  1.355961e-01
    q0( 29, 5) =  2.120541e+00
    n ( 29, 5) =  0
    k ( 30, 5) =  4.026281e-01
    q0( 30, 5) =  2.111764e+00
    n ( 30, 5) =  0
    k ( 31, 5) =  1.803650e-01
    q0( 31, 5) =  2.103575e+00
    n ( 31, 5) =  0
    k ( 32, 5) =  1.385288e-02
    q0( 32, 5) =  0.000000e+00
    n ( 32, 5) =  2
    k ( 33, 5) =  1.568087e-03
    q0( 33, 5) =  0.000000e+00
    n ( 33, 5) =  2
    k ( 34, 5) =  7.301002e-03
    q0( 34, 5) =  0.000000e+00
    n ( 34, 5) =  2
    k ( 35, 5) =  6.019684e-03
    q0( 35, 5) =  0.000000e+00
    n ( 35, 5) =  2
    k ( 36, 5) =  4.985267e-17
    q0( 36, 5) =  0.000000e+00
    n ( 36, 5) =  2
    k ( 37, 5) =  3.684564e-02
    q0( 37, 5) =  0.000000e+00
    n ( 37, 5) =  2
    k ( 38, 5) = -1.317018e-17
    q0( 38, 5) =  0.000000e+00
    n ( 38, 5) =  2
    k ( 39, 5) =  2.035660e-02
    q0( 39, 5) =  0.000000e+00
    n ( 39, 5) =  2
    k ( 40, 5) = -8.348360e-17
    q0( 40, 5) =  0.000000e+00
    n ( 40, 5) =  2
    k ( 41, 5) =  7.617600e-02
    q0( 41, 5) =  0.000000e+00
    n ( 41, 5) =  2
    k ( 42, 5) =  1.520926e-18
    q0( 42, 5) =  0.000000e+00
    n ( 42, 5) =  2
    k ( 43, 5) =  1.231268e-02
    q0( 43, 5) =  0.000000e+00
    n ( 43, 5) =  2
    k ( 44, 5) =  1.085541e-02
    q0( 44, 5) =  0.000000e+00
    n ( 44, 5) =  2
    k ( 45, 5) =  1.171598e-03
    q0( 45, 5) =  0.000000e+00
    n ( 45, 5) =  2
    k ( 46, 5) =  1.786666e-02
    q0( 46, 5) =  0.000000e+00
    n ( 46, 5) =  2
    k ( 47, 5) =  3.524048e-17
    q0( 47, 5) =  0.000000e+00
    n ( 47, 5) =  2
    k ( 48, 5) =  1.960997e-02
    q0( 48, 5) =  0.000000e+00
    n ( 48, 5) =  2
    k ( 49, 5) =  7.352789e-03
    q0( 49, 5) =  0.000000e+00
    n ( 49, 5) =  2
    k ( 50, 5) =  6.747023e-04
    q0( 50, 5) =  0.000000e+00
    n ( 50, 5) =  2
    k ( 51, 5) =  6.364265e-18
    q0( 51, 5) =  0.000000e+00
    n ( 51, 5) =  2
    k ( 52, 5) =  2.573680e-17
    q0( 52, 5) =  0.000000e+00
    n ( 52, 5) =  2
    k ( 53, 5) =  1.094557e-17
    q0( 53, 5) =  0.000000e+00
    n ( 53, 5) =  2
    k ( 54, 5) =  2.602755e-18
    q0( 54, 5) =  0.000000e+00
    n ( 54, 5) =  2
    k ( 55, 5) =  4.980248e-03
    q0( 55, 5) =  0.000000e+00
    n ( 55, 5) =  2
    k ( 56, 5) =  2.237955e-01
    q0( 56, 5) =  7.723842e-10
    n ( 56, 5) =  0
    k ( 57, 5) =  3.284633e-01
    q0( 57, 5) =  4.830951e-06
    n ( 57, 5) =  0
    k ( 58, 5) =  1.460862e-01
    q0( 58, 5) =  7.530821e-07
    n ( 58, 5) =  0
    k ( 59, 5) =  1.048746e-16
    q0( 59, 5) =  1.811819e-06
    n ( 59, 5) =  0
    k ( 60, 5) =  2.402438e-01
    q0( 60, 5) =  8.323602e-08
    n ( 60, 5) =  0
    k ( 61, 5) =  2.542964e-01
    q0( 61, 5) =  3.699948e-06
    n ( 61, 5) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 3, anchor point 6)
    k (  1, 6) =  2.041288e+00
    q0(  1, 6) =  1.417458e+00
    n (  1, 6) =  0
    k (  2, 6) =  2.041289e+00
    q0(  2, 6) =  1.417458e+00
    n (  2, 6) =  0
    k (  3, 6) =  2.033774e+00
    q0(  3, 6) =  1.732475e+00
    n (  3, 6) =  0
    k (  4, 6) =  2.421895e+00
    q0(  4, 6) =  1.386383e+00
    n (  4, 6) =  0
    k (  5, 6) =  1.503032e+00
    q0(  5, 6) =  1.084746e+00
    n (  5, 6) =  0
    k (  6, 6) =  2.237141e+00
    q0(  6, 6) =  1.399117e+00
    n (  6, 6) =  0
    k (  7, 6) =  1.485305e+00
    q0(  7, 6) =  1.085510e+00
    n (  7, 6) =  0
    k (  8, 6) =  2.237140e+00
    q0(  8, 6) =  1.399118e+00
    n (  8, 6) =  0
    k (  9, 6) =  1.489029e+00
    q0(  9, 6) =  1.085581e+00
    n (  9, 6) =  0
    k ( 10, 6) =  2.421895e+00
    q0( 10, 6) =  1.386383e+00
    n ( 10, 6) =  0
    k ( 11, 6) =  1.485305e+00
    q0( 11, 6) =  1.085510e+00
    n ( 11, 6) =  0
    k ( 12, 6) =  1.503032e+00
    q0( 12, 6) =  1.084746e+00
    n ( 12, 6) =  0
    k ( 13, 6) =  2.766702e-01
    q0( 13, 6) =  2.103863e+00
    n ( 13, 6) =  0
    k ( 14, 6) =  1.149595e-01
    q0( 14, 6) =  2.068768e+00
    n ( 14, 6) =  0
    k ( 15, 6) =  2.766702e-01
    q0( 15, 6) =  2.103864e+00
    n ( 15, 6) =  0
    k ( 16, 6) =  1.149593e-01
    q0( 16, 6) =  2.068768e+00
    n ( 16, 6) =  0
    k ( 17, 6) =  1.253335E-02
    q0( 17, 6) =  1.883350E+00
    n ( 17, 6) =  0
    k ( 18, 6) =  5.867528e-02
    q0( 18, 6) =  2.068959e+00
    n ( 18, 6) =  0
    k ( 19, 6) =  2.910010e-01
    q0( 19, 6) =  2.107630e+00
    n ( 19, 6) =  0
    k ( 20, 6) =  2.345380e-01
    q0( 20, 6) =  2.096166e+00
    n ( 20, 6) =  0
    k ( 21, 6) =  1.527723e-01
    q0( 21, 6) =  2.093581e+00
    n ( 21, 6) =  0
    k ( 22, 6) =  1.894821e-01
    q0( 22, 6) =  2.111081e+00
    n ( 22, 6) =  0
    k ( 23, 6) =  2.461660e-01
    q0( 23, 6) =  2.098801e+00
    n ( 23, 6) =  0
    k ( 24, 6) =  1.737555e-01
    q0( 24, 6) =  2.092422e+00
    n ( 24, 6) =  0
    k ( 25, 6) =  1.594971e-01
    q0( 25, 6) =  2.093903e+00
    n ( 25, 6) =  0
    k ( 26, 6) =  2.345383e-01
    q0( 26, 6) =  2.096165e+00
    n ( 26, 6) =  0
    k ( 27, 6) =  1.594975e-01
    q0( 27, 6) =  2.093904e+00
    n ( 27, 6) =  0
    k ( 28, 6) =  1.737551e-01
    q0( 28, 6) =  2.092420e+00
    n ( 28, 6) =  0
    k ( 29, 6) =  1.894825e-01
    q0( 29, 6) =  2.111081e+00
    n ( 29, 6) =  0
    k ( 30, 6) =  2.909994e-01
    q0( 30, 6) =  2.107630e+00
    n ( 30, 6) =  0
    k ( 31, 6) =  1.527721e-01
    q0( 31, 6) =  2.093580e+00
    n ( 31, 6) =  0
    k ( 32, 6) =  1.046807e-02
    q0( 32, 6) =  0.000000e+00
    n ( 32, 6) =  2
    k ( 33, 6) =  1.410852e-02
    q0( 33, 6) =  0.000000e+00
    n ( 33, 6) =  2
    k ( 34, 6) =  7.961954e-03
    q0( 34, 6) =  0.000000e+00
    n ( 34, 6) =  2
    k ( 35, 6) =  1.410728e-02
    q0( 35, 6) =  0.000000e+00
    n ( 35, 6) =  2
    k ( 36, 6) =  2.918819e-03
    q0( 36, 6) =  0.000000e+00
    n ( 36, 6) =  2
    k ( 37, 6) =  1.076466e-02
    q0( 37, 6) =  0.000000e+00
    n ( 37, 6) =  2
    k ( 38, 6) =  3.488921e-03
    q0( 38, 6) =  0.000000e+00
    n ( 38, 6) =  2
    k ( 39, 6) =  1.070845e-02
    q0( 39, 6) =  0.000000e+00
    n ( 39, 6) =  2
    k ( 40, 6) =  3.693467e-04
    q0( 40, 6) =  0.000000e+00
    n ( 40, 6) =  2
    k ( 41, 6) =  6.093330e-03
    q0( 41, 6) =  0.000000e+00
    n ( 41, 6) =  2
    k ( 42, 6) =  5.971783e-03
    q0( 42, 6) =  0.000000e+00
    n ( 42, 6) =  2
    k ( 43, 6) =  1.108861e-02
    q0( 43, 6) =  0.000000e+00
    n ( 43, 6) =  2
    k ( 44, 6) =  1.371040e-02
    q0( 44, 6) =  0.000000e+00
    n ( 44, 6) =  2
    k ( 45, 6) =  1.370711e-02
    q0( 45, 6) =  0.000000e+00
    n ( 45, 6) =  2
    k ( 46, 6) =  1.109123e-02
    q0( 46, 6) =  0.000000e+00
    n ( 46, 6) =  2
    k ( 47, 6) =  6.093332e-03
    q0( 47, 6) =  0.000000e+00
    n ( 47, 6) =  2
    k ( 48, 6) =  1.076531e-02
    q0( 48, 6) =  0.000000e+00
    n ( 48, 6) =  2
    k ( 49, 6) =  1.070534e-02
    q0( 49, 6) =  0.000000e+00
    n ( 49, 6) =  2
    k ( 50, 6) =  2.385851e-03
    q0( 50, 6) =  0.000000e+00
    n ( 50, 6) =  2
    k ( 51, 6) =  4.196294e-03
    q0( 51, 6) =  0.000000e+00
    n ( 51, 6) =  2
    k ( 52, 6) =  8.092265e-20
    q0( 52, 6) =  0.000000e+00
    n ( 52, 6) =  2
    k ( 53, 6) = -6.069924e-20
    q0( 53, 6) =  0.000000e+00
    n ( 53, 6) =  2
    k ( 54, 6) =  4.199041e-03
    q0( 54, 6) =  0.000000e+00
    n ( 54, 6) =  2
    k ( 55, 6) =  2.389473e-03
    q0( 55, 6) =  0.000000e+00
    n ( 55, 6) =  2
    k ( 56, 6) =  1.608127e-01
    q0( 56, 6) =  9.050319e-07
    n ( 56, 6) =  0
    k ( 57, 6) =  1.608599e-01
    q0( 57, 6) =  9.038175e-07
    n ( 57, 6) =  0
    k ( 58, 6) =  1.218387e-01
    q0( 58, 6) =  1.668569e-07
    n ( 58, 6) =  0
    k ( 59, 6) =  2.364301e-01
    q0( 59, 6) =  2.152578e-07
    n ( 59, 6) =  0
    k ( 60, 6) =  1.805203e-01
    q0( 60, 6) =  1.068481e-06
    n ( 60, 6) =  0
    k ( 61, 6) =  1.219424e-01
    q0( 61, 6) =  1.656789e-07
    n ( 61, 6) =  0
    !--- END generated by 'gen_FFparam.py'

    !--- BEGIN generated by 'gen_FFparam.py'
    !--- (state 3, anchor point 7)
    k (  1, 7) =  2.040746e+00
    q0(  1, 7) =  1.417435e+00
    n (  1, 7) =  0
    k (  2, 7) =  2.040746e+00
    q0(  2, 7) =  1.417435e+00
    n (  2, 7) =  0
    k (  3, 7) =  2.077015e+00
    q0(  3, 7) =  1.732362e+00
    n (  3, 7) =  0
    k (  4, 7) =  2.422221e+00
    q0(  4, 7) =  1.386393e+00
    n (  4, 7) =  0
    k (  5, 7) =  1.503421e+00
    q0(  5, 7) =  1.084729e+00
    n (  5, 7) =  0
    k (  6, 7) =  2.237308e+00
    q0(  6, 7) =  1.399112e+00
    n (  6, 7) =  0
    k (  7, 7) =  1.485368e+00
    q0(  7, 7) =  1.085506e+00
    n (  7, 7) =  0
    k (  8, 7) =  2.237308e+00
    q0(  8, 7) =  1.399112e+00
    n (  8, 7) =  0
    k (  9, 7) =  1.489048e+00
    q0(  9, 7) =  1.085583e+00
    n (  9, 7) =  0
    k ( 10, 7) =  2.422222e+00
    q0( 10, 7) =  1.386393e+00
    n ( 10, 7) =  0
    k ( 11, 7) =  1.485368e+00
    q0( 11, 7) =  1.085506e+00
    n ( 11, 7) =  0
    k ( 12, 7) =  1.503422e+00
    q0( 12, 7) =  1.084729e+00
    n ( 12, 7) =  0
    k ( 13, 7) =  2.781837e-01
    q0( 13, 7) =  2.103711e+00
    n ( 13, 7) =  0
    k ( 14, 7) =  1.157935e-01
    q0( 14, 7) =  2.068870e+00
    n ( 14, 7) =  0
    k ( 15, 7) =  2.781840e-01
    q0( 15, 7) =  2.103712e+00
    n ( 15, 7) =  0
    k ( 16, 7) =  1.157934e-01
    q0( 16, 7) =  2.068869e+00
    n ( 16, 7) =  0
    k ( 17, 7) =  5.302041e-02
    q0( 17, 7) =  2.069146e+00
    n ( 17, 7) =  0
    k ( 18, 7) =  2.887871e-01
    q0( 18, 7) =  2.107536e+00
    n ( 18, 7) =  0
    k ( 19, 7) =  2.341827e-01
    q0( 19, 7) =  2.096197e+00
    n ( 19, 7) =  0
    k ( 20, 7) =  1.525717e-01
    q0( 20, 7) =  2.093544e+00
    n ( 20, 7) =  0
    k ( 21, 7) =  1.891481e-01
    q0( 21, 7) =  2.111137e+00
    n ( 21, 7) =  0
    k ( 22, 7) =  2.464018e-01
    q0( 22, 7) =  2.098832e+00
    n ( 22, 7) =  0
    k ( 23, 7) =  1.737700e-01
    q0( 23, 7) =  2.092410e+00
    n ( 23, 7) =  0
    k ( 24, 7) =  1.595498e-01
    q0( 24, 7) =  2.093885e+00
    n ( 24, 7) =  0
    k ( 25, 7) =  2.341826e-01
    q0( 25, 7) =  2.096197e+00
    n ( 25, 7) =  0
    k ( 26, 7) =  1.595498e-01
    q0( 26, 7) =  2.093889e+00
    n ( 26, 7) =  0
    k ( 27, 7) =  1.737694e-01
    q0( 27, 7) =  2.092410e+00
    n ( 27, 7) =  0
    k ( 28, 7) =  1.891483e-01
    q0( 28, 7) =  2.111134e+00
    n ( 28, 7) =  0
    k ( 29, 7) =  2.887868e-01
    q0( 29, 7) =  2.107538e+00
    n ( 29, 7) =  0
    k ( 30, 7) =  1.525714e-01
    q0( 30, 7) =  2.093544e+00
    n ( 30, 7) =  0
    k ( 31, 7) =  9.608718e-03
    q0( 31, 7) =  0.000000e+00
    n ( 31, 7) =  2
    k ( 32, 7) =  1.409565e-02
    q0( 32, 7) =  0.000000e+00
    n ( 32, 7) =  2
    k ( 33, 7) =  8.780380e-03
    q0( 33, 7) =  0.000000e+00
    n ( 33, 7) =  2
    k ( 34, 7) =  1.410145e-02
    q0( 34, 7) =  0.000000e+00
    n ( 34, 7) =  2
    k ( 35, 7) =  2.030424e-03
    q0( 35, 7) =  0.000000e+00
    n ( 35, 7) =  2
    k ( 36, 7) =  1.075349e-02
    q0( 36, 7) =  0.000000e+00
    n ( 36, 7) =  2
    k ( 37, 7) =  4.326097e-03
    q0( 37, 7) =  0.000000e+00
    n ( 37, 7) =  2
    k ( 38, 7) =  1.072261e-02
    q0( 38, 7) =  0.000000e+00
    n ( 38, 7) =  2
    k ( 39, 7) =  1.176453e-03
    q0( 39, 7) =  0.000000e+00
    n ( 39, 7) =  2
    k ( 40, 7) =  6.102503e-03
    q0( 40, 7) =  0.000000e+00
    n ( 40, 7) =  2
    k ( 41, 7) =  5.143853e-03
    q0( 41, 7) =  0.000000e+00
    n ( 41, 7) =  2
    k ( 42, 7) =  1.108646e-02
    q0( 42, 7) =  0.000000e+00
    n ( 42, 7) =  2
    k ( 43, 7) =  1.369349e-02
    q0( 43, 7) =  0.000000e+00
    n ( 43, 7) =  2
    k ( 44, 7) =  1.369764e-02
    q0( 44, 7) =  0.000000e+00
    n ( 44, 7) =  2
    k ( 45, 7) =  1.109463e-02
    q0( 45, 7) =  0.000000e+00
    n ( 45, 7) =  2
    k ( 46, 7) =  6.107321e-03
    q0( 46, 7) =  0.000000e+00
    n ( 46, 7) =  2
    k ( 47, 7) =  1.075222e-02
    q0( 47, 7) =  0.000000e+00
    n ( 47, 7) =  2
    k ( 48, 7) =  1.071361e-02
    q0( 48, 7) =  0.000000e+00
    n ( 48, 7) =  2
    k ( 49, 7) =  2.379739e-03
    q0( 49, 7) =  0.000000e+00
    n ( 49, 7) =  2
    k ( 50, 7) =  4.195018e-03
    q0( 50, 7) =  0.000000e+00
    n ( 50, 7) =  2
    k ( 51, 7) = -1.771924e-20
    q0( 51, 7) =  0.000000e+00
    n ( 51, 7) =  2
    k ( 52, 7) = -1.148656e-19
    q0( 52, 7) =  0.000000e+00
    n ( 52, 7) =  2
    k ( 53, 7) =  4.191420e-03
    q0( 53, 7) =  0.000000e+00
    n ( 53, 7) =  2
    k ( 54, 7) =  2.385780e-03
    q0( 54, 7) =  0.000000e+00
    n ( 54, 7) =  2
    k ( 55, 7) =  1.612998e-01
    q0( 55, 7) =  5.262817e-07
    n ( 55, 7) =  0
    k ( 56, 7) =  1.612068e-01
    q0( 56, 7) =  5.260694e-07
    n ( 56, 7) =  0
    k ( 57, 7) =  1.221269e-01
    q0( 57, 7) =  3.144837e-07
    n ( 57, 7) =  0
    k ( 58, 7) =  2.372395e-01
    q0( 58, 7) =  5.205649e-07
    n ( 58, 7) =  0
    k ( 59, 7) =  1.806178e-01
    q0( 59, 7) =  3.922352e-07
    n ( 59, 7) =  0
    k ( 60, 7) =  1.221685e-01
    q0( 60, 7) =  3.144473e-07
    n ( 60, 7) =  0
    !--- END generated by 'gen_FFparam.py'
  end select

  end subroutine assignparamuii_t

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!             Calculate tertiary diabatic couplings
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!=================================================================
! *def calcuij_t
!    Calculate tert. diab. coupl. and its grad. w.r.t. Cart.
!  Input:
!    uu: diab. matrix (in hartree)
!    guu: gradient array (in hartree/Angstrom)
!    x: Cartesian coordinates in Angstroms
!    qps: values of primary & secondary internal coordinates
!    qtc: values of tert. internals for couplings
!    nqtc: number of qtc
!    qtcsym: symmetry (a' or a") of qtc
!    qtctyp: type of qtc (CSC bend or other)
!    bmatqtc: B matrix, B(i, j) = dqtc(i)/dx(j)
!  Output:
!    uu: diab. matrix w/ calculated uij added
!    guu: grad. array w/ calculated duij/dx added
!=================================================================
  subroutine calcuij_t(uu, guu, x, qps, qtc, nqtc, &
     qtcsym, qtctyp, bmatqtc)
  implicit none
  double precision :: uu(nstate,nstate), qps(:), qtc(:) &
       , guu(3,natom,nstate,nstate), x(3,natom) &
       , bmatqtc(ntermmax,natom*3)
  integer :: nqtc, qtcsym(ntermmax), qtctyp(ntermmax)
  double precision :: uij, duij(3,natom)

  call uij_t(uij, duij, x, qps, qtc, &
       nqtc, qtcsym, qtctyp, bmatqtc, 12)
  uu(1,2) = uu(1,2)+uij; uu(2,1) = uu(1,2)
  guu(:,:,1,2) = guu(:,:,1,2)+duij(:,:)
  guu(:,:,2,1) = guu(:,:,1,2)
  call uij_t(uij, duij, x, qps, qtc, &
       nqtc, qtcsym, qtctyp, bmatqtc, 13)
  uu(1,3) = uu(1,3)+uij; uu(3,1) = uu(1,3)
  guu(:,:,1,3) = guu(:,:,1,3)+duij(:,:)
  guu(:,:,3,1) = guu(:,:,1,3)
  call uij_t(uij, duij, x, qps, qtc, &
       nqtc, qtcsym, qtctyp, bmatqtc, 23)
  uu(2,3) = uu(2,3)+uij; uu(3,2) = uu(2,3)
  guu(:,:,2,3) = guu(:,:,2,3)+duij(:,:)
  guu(:,:,3,2) = guu(:,:,2,3)

  end subroutine calcuij_t
  !--- Actual core routine of calcuij_t
  !--- Parameters have the same meaning as in calcuij_t
  subroutine uij_t(uij, duij, x, qps, qtc, nqtc, &
       qtcsym, qtctyp, bmatqtc, istate)
    implicit none
    double precision :: uij, duij(3,natom), x(3,natom), qps(:) &
         , qtc(:), bmatqtc(ntermmax,natom*3)
    integer :: nqtc, qtcsym(ntermmax), qtctyp(ntermmax), istate
    integer,parameter :: nR0=3, nphi0=5, nk=2
    integer :: iR0, iphi0, iqtc, iatomlist(10)
    double precision, allocatable :: gg(:), bmat(:,:)
    double precision :: ee(nR0,nphi0), e, g &
         , k(nk,nR0,nphi0,nqtc) &
         , gradc(3,natom,nR0,nphi0), R, phi, qtc0(nqtc) &
         , tentR(0:nR0), dtentR(0:nR0), R0(nR0), phi0(nphi0) &
         , tentphi(0:nphi0), dtentphi(0:nphi0), cos2phi0(nphi0) &
         , qtcv(nqtc), bmatrow(3*natom) &
         , sij(nphi0), dsij(nphi0)

    !--- initialization
    !--- prim&sec coord.
    R = qps(1); phi = qps(2)
    qtcv(1:nqtc) = qtc(1:nqtc)
    !--- R and phi values at anchor points
    R0(:) = (/1.45d0,2.9d0,4.0d0/)
    phi0(:) = (/0d0,30d0,50d0,70d0,90d0/)
    phi0(:) = phi0(:)/180d0*pi
    do iphi0 = 1, nphi0
       cos2phi0(iphi0) = -cos(2*phi0(iphi0))
    end do

    !--- state-dependent params
    call assignparamuij_t(k, qtc0, nqtc, nk, nR0, nphi0, istate)

    ee = 0d0
    allocate(gg(nqtc),bmat(nqtc,3*natom))
    bmat(1:nqtc,:) = bmatqtc(1:nqtc,:)
    !--- loop over anchor points (each with 2 subscripts iR0,iphi0)
    do iR0 = 1, nR0
       do iphi0 = 1, nphi0
          do iqtc = 1, nqtc
             !--- Calculate uij and duij/dqtc at anchor points:
             !---   the working unit for calcTDCterm is degree;
             !---   this is an exception due to the parameterization
             call calcTDCterm(e, g, &
                  k(:,iR0,iphi0,iqtc), nk, qtc0(iqtc), &
                  qtcv(iqtc)/pi*180d0, qtctyp(iqtc))
             if(iqtc.eq.13 .or. iqtc.eq.23 .or. iqtc.eq.27) then
                e = 0d0; g = 0d0
             end if
             !--- Zero torsion and 90 deg. torsion are set to zero
             if (iphi0.eq.1 .or. iphi0.eq.nphi0) then
                e = 0d0; g = 0d0
             end if

             e = e*eV_hartree
             g = g*eV_hartree

             !--- convert gg from deg^-1 to rad^-1
             g = g*180d0/pi
             ee(iR0,iphi0) = ee(iR0,iphi0)+e
             gg(iqtc) = g
          end do
          !--- use B matrix to convert to grad. w.r.t. Cartesians
          gradc(:,:,iR0,iphi0) = reshape(matmul(gg, bmat), (/3,natom/))
       end do
    end do
    deallocate(gg,bmat)
    !--- calculate tent func., sij,  and dtentR/dR & dtentphi/dphi
    !--- & dsij/dphi
    do iR0 = 1, nR0
       call calctent1(tentR(iR0), dtentR(iR0), R, R0, &
            iR0-1, iR0, iR0+1, nR0)
    end do
    do iphi0 = 1, nphi0
       call calctent1(tentphi(iphi0), dtentphi(iphi0), &
            -cos(2*phi), cos2phi0, &
            iphi0-1, iphi0, iphi0+1, nphi0)
       !--- argument of tent func. is actually -cos(2*phi)
       !--- dtentphi/dphi = dtentphi/d-cos2phi * d-cos2phi/dphi
       dtentphi(iphi0) = dtentphi(iphi0)*(2*sin(2*phi))
       if (istate.eq.12) then
          sij(iphi0) = 1d0
          dsij(iphi0) = 0d0
       else
          sij(iphi0) = sign(1d0, sin(2*phi))
          dsij(iphi0) = 0d0
       end if
    end do

    uij = 0d0; duij(:,:) = 0d0
    dtentR(0) = 0d0; dtentphi(0) = 0d0
    do iR0 = 1, nR0
       do iphi0 = 1, nphi0
          !--- interpolate anchor point values using tent func.
          !--- to get uij and duij/dqtc
          uij = uij+ &
               ee(iR0,iphi0)*tentR(iR0)*tentphi(iphi0)&
               *sij(iphi0)
          duij(:,:) = duij(:,:) &
               +gradc(:,:,iR0,iphi0)*tentR(iR0)*tentphi(iphi0)&
                *sij(iphi0)
          !--- interpolate dtentq/dq to get duij/dq (q=R, phi)
          !--- (dtentq(0) is duij/dq)
          dtentR(0) = dtentR(0) &
               +ee(iR0,iphi0)*dtentR(iR0)*tentphi(iphi0)&
                *sij(iphi0)
          dtentphi(0) = dtentphi(0) &
               +ee(iR0,iphi0)*tentR(iR0)&
                *(dtentphi(iphi0)*sij(iphi0)&
                  +tentphi(iphi0)*dsij(iphi0))
       end do
    end do
    !--- use B matrix to convert duij/dq to duij/dx
    iatomlist(1:4) = (/7,13,0,0/)
    call calcbmat(bmatrow, x, iatomlist, 1)
    duij(:,:) = duij(:,:) &
         +reshape(dtentR(0)*bmatrow, (/3,natom/))
    !--- special coordinate in place of CCSH torsion
    !>>> modify case(6) of calcbmat
    iatomlist(1:5) = (/2,1,6,7,13/)
    call calcbmat(bmatrow, x, iatomlist, 6)
    duij(:,:) = duij(:,:) &
         +reshape(dtentphi(0)*bmatrow, (/3,natom/))

  end subroutine uij_t

!=================================================================
! *def calcTDCterm
!    Calculate tertiary diabatic coupling term and its grad.,
!      which is a nk-degree polynomial (w/ intercept=0).
!    (Unlike most internal routines, the input and output unit
!     here for angles is degree due to parametrization.)
!  Input:
!    k: coefficients of polynomial
!    nk: number of degree of polynomial
!    q0: center of expansion
!    q: value of internal coordinate
!    ityp: type of TDC term; see below
!  Output:
!    ee: diab. coupling
!    gg: gradient w.r.t. the internal coordinate
!=================================================================
  subroutine calcTDCterm(ee, gg, k, nk, q0, q, ityp)
  implicit none
  integer :: nk, ityp
  double precision :: ee, gg, k(nk), q0, q
  integer :: i
  double precision :: dq, cc, expo, qq, qq0

  select case(ityp)
  !--- ityp 1: use the coordinate as variable
  case(1)
     dq = q-q0
     ee = 0d0; gg = 0d0
     do i = 1, nk
        ee = ee+k(i)*dq**i
        gg = gg+i*k(i)*dq**(i-1)
     end do

     !--- damp the couplings at large dq
     cc = 100d0
     expo = exp(-dq**2/cc**2)
     gg = gg*expo+ee*(-2*dq/cc**2)*expo
     ee = ee*expo
  !--- Note that q needs deg->rad first,
  !---   and gg needs rad-1 -> deg-1 last
  case(2)
     qq = q/180d0*pi; qq0 = q0/180d0*pi
     dq = cos(qq)-cos(qq0)
     ee = 0d0; gg = 0d0
     do i = 1, nk
        ee = ee+k(i)*dq**i
        gg = gg+i*k(i)*dq**(i-1)
     end do
     gg = -gg*sin(qq)/180d0*pi
  end select

  end subroutine calcTDCterm

!=================================================================
! *def calcqtc
!    Calculate nonredundant tert. int. coord. for diab. coupl.
!  Input:
!    x: Cartesian coordinates in Angstroms
!  Output:
!    qtc: values of tert. internals
!    nqtc: number of qtc
!    qtcsym: symmetry (a' or a") of qtc
!    qtctyp: type of qtc (CSC bend or other)
!    bmatqtc: B matrix, B(i, j) = dqtc(i)/dx(j)
!=================================================================
  subroutine calcqtc(qtc, nqtc, qtcsym, qtctyp, bmatqtc, x)
  implicit none
  integer :: qtcsym(ntermmax), nqtc
  double precision :: qtc(ntermmax), x(3,natom) &
       , bmatqtc(ntermmax,natom*3)
  double precision :: qr(ntermmax), lccoef(nlcmax,ntermmax) &
       , lcqr(nlcmax), bmatrow(natom*3), bmatr(ntermmax,natom*3)
  integer :: nqr, itypeqr(ntermmax), iqrlist(10,ntermmax) &
       , nqnr, lcindex(nlcmax,ntermmax), nlc(ntermmax) &
       , qnrindex(ntermmax), qtctyp(ntermmax) &
       , iqtc, iqnr, iqr, ilc

  !--- calculate redundant q
  call calcqr(qr, x, nqr, itypeqr, iqrlist)
  !--- attributes of nonredundant q (qnr),
  !---  esp. their relation to qr
  call assignqnr(nqnr, lcindex, lccoef, nlc)
  !--- attributes of qtc,
  !---  esp. their relation to qnr
  call assignqtc(nqtc, qtcsym, qnrindex, qtctyp)

  !--- value of qtc
  !--- loop over qtc index
  do iqtc = 1, nqtc
    !--- corresponding qnr index
    iqnr = qnrindex(iqtc)

    !@WRONG
    !"QR has value 0 which is less than the lower bound of 1"
    !lcqr(1:nlcmax) = qr(lcindex(1:nlcmax,iqnr))

    !--- value of qr involved in qtc(i)
    lcqr(1:nlc(iqnr)) = qr(lcindex(1:nlc(iqnr),iqnr))
    !--- take lin.comb. of qr to get qtc
    qtc(iqtc) = dot_product(lccoef(:,iqnr),lcqr)
  end do

  !--- B matrix for qr
  !--- loop over qr index
  do iqr = 1, nqr
    call calcbmat(bmatrow, x, iqrlist(:,iqr), itypeqr(iqr))
    bmatr(iqr, :) = bmatrow(:)
  end do

  !--- B matrix for qtc from lin.comb. of bmatr
  bmatqtc = 0d0
  !--- loop over qtc
  do iqtc = 1, nqtc
    !--- corresponding qnr index
    iqnr = qnrindex(iqtc)
    !--- loop over qr involved in the lin.comb.
    do ilc = 1, nlc(iqnr)
       !--- lin.comb. of bmatr rows
       bmatqtc(iqtc,:) = bmatqtc(iqtc,:) &
            +lccoef(ilc,iqnr)*bmatr(lcindex(ilc,iqnr),:)
    end do
  end do

  end subroutine calcqtc

!=================================================================
! *def calcqr
!    Calculate redundant q (qr) for later use of qtc
!  Input:
!    x: Cartesian coordinates in Angstroms
!    nqr: number of qr
!    itypeqr: type of each qr
!    iqrlist: atom index list array
!  Output:
!    qr: values of qr
!=================================================================
  subroutine calcqr(qr, x, nqr, itypeqr, iqrlist)
  implicit none
  double precision :: qr(ntermmax), x(3,natom)
  integer :: nqr, itypeqr(ntermmax), iqrlist(10,ntermmax)
  integer :: iterm, ii, jj, kk, ll

  !--- attributes of redundant q (qr)
  call assignqr(nqr, itypeqr, iqrlist)

  do iterm = 1, nqr
    ii = iqrlist(1,iterm)
    jj = iqrlist(2,iterm)
    kk = iqrlist(3,iterm)
    ll = iqrlist(4,iterm)
    select case(itypeqr(iterm))
       !--- bond length
      case(1)
        qr(iterm) = evalBL(x, ii, jj)
      !--- bond angle
      case(2)
        qr(iterm) = evalBA(x, ii, jj, kk)
      !--- torsion
      case(3)
        qr(iterm) = evalTO(x, ii, jj, kk, ll)
      !--- oop bend
      case(4)
        qr(iterm) = evalOB(x, ii, jj, kk, ll)
    end select
  end do
  end subroutine calcqr

!=================================================================
! *def assignqr
!    Assign attributes of redundant q (qr) for later use of qtc
!  Input:
!    NA
!  Output:
!    nqr: number of qr
!    itypeqr: type of each qr
!    iqrlist: atom index list array
!=================================================================
  subroutine assignqr(nqr, itypeqr, iqrlist)
  implicit none
  integer :: nqr, itypeqr(ntermmax), iqrlist(10,ntermmax)

  nqr = 45
  !--- type of qr
  !--- 1=BL(BondLen), 2=BA(BondAng), 3=TO(Torsion), 4=OB(OOPBend)
  itypeqr( 1:13) = 1
  itypeqr(14:32) = 2
  itypeqr(33:39) = 3
  itypeqr(40:45) = 4

  !--- BEGIN generated by 'python gen_qr_forUij-t.py rdef0.txt'
  !--- bond length
  iqrlist(1:2, 1) = (/ 1, 2/)
  iqrlist(1:2, 2) = (/ 2, 3/)
  iqrlist(1:2, 3) = (/ 3, 4/)
  iqrlist(1:2, 4) = (/ 4, 5/)
  iqrlist(1:2, 5) = (/ 5, 6/)
  iqrlist(1:2, 6) = (/ 6, 1/)
  iqrlist(1:2, 7) = (/ 1, 7/)
  iqrlist(1:2, 8) = (/ 2, 8/)
  iqrlist(1:2, 9) = (/ 3, 9/)
  iqrlist(1:2,10) = (/ 4,10/)
  iqrlist(1:2,11) = (/ 5,11/)
  iqrlist(1:2,12) = (/ 6,12/)
  iqrlist(1:2,13) = (/ 7,13/)
  !--- bend
  iqrlist(1:3,14) = (/ 6, 1, 2/)
  iqrlist(1:3,15) = (/ 1, 2, 3/)
  iqrlist(1:3,16) = (/ 2, 3, 4/)
  iqrlist(1:3,17) = (/ 3, 4, 5/)
  iqrlist(1:3,18) = (/ 4, 5, 6/)
  iqrlist(1:3,19) = (/ 5, 6, 1/)
  iqrlist(1:3,20) = (/ 7, 1, 2/)
  iqrlist(1:3,21) = (/ 7, 1, 6/)
  iqrlist(1:3,22) = (/ 8, 2, 1/)
  iqrlist(1:3,23) = (/ 8, 2, 3/)
  iqrlist(1:3,24) = (/ 9, 3, 2/)
  iqrlist(1:3,25) = (/ 9, 3, 4/)
  iqrlist(1:3,26) = (/10, 4, 3/)
  iqrlist(1:3,27) = (/10, 4, 5/)
  iqrlist(1:3,28) = (/11, 5, 4/)
  iqrlist(1:3,29) = (/11, 5, 6/)
  iqrlist(1:3,30) = (/12, 6, 5/)
  iqrlist(1:3,31) = (/12, 6, 1/)
  iqrlist(1:3,32) = (/13, 7, 1/)
  !--- torsion
  iqrlist(1:4,33) = (/ 6, 1, 2, 3/)
  iqrlist(1:4,34) = (/ 1, 2, 3, 4/)
  iqrlist(1:4,35) = (/ 2, 3, 4, 5/)
  iqrlist(1:4,36) = (/ 3, 4, 5, 6/)
  iqrlist(1:4,37) = (/ 4, 5, 6, 1/)
  iqrlist(1:4,38) = (/ 5, 6, 1, 2/)
  iqrlist(1:4,39) = (/13, 7, 1, 2/)
  !--- oop bend
  iqrlist(1:4,40) = (/ 7, 2, 6, 1/)
  iqrlist(1:4,41) = (/ 8, 3, 1, 2/)
  iqrlist(1:4,42) = (/ 9, 4, 2, 3/)
  iqrlist(1:4,43) = (/10, 5, 3, 4/)
  iqrlist(1:4,44) = (/11, 6, 4, 5/)
  iqrlist(1:4,45) = (/12, 1, 5, 6/)
  !--- END generated by 'python gen_qr_forUij-t.py rdef0.txt'

  end subroutine assignqr

!=================================================================
! *def assignqnr
!    Assign attributes of nonredund. q (qnr) for later use of qtc
!  Input:
!    NA
!  Output:
!    nqnr: number of qr
!    lcindex: index of qr involved in lin.comb. to get each qnr
!    lccoef: coefficients of lin.comb. to get each qnr
!    nlc: number of terms in lin.comb. to get each qnr
!  ~200 lines of params
!=================================================================
  subroutine assignqnr(nqnr, lcindex, lccoef, nlc)
  implicit none
  integer :: nqnr, lcindex(nlcmax,ntermmax), nlc(ntermmax)
  double precision :: lccoef(nlcmax,ntermmax)

  nqnr = 33
  !--- BEGIN generated by python gen_qnr_forUij-t.py intnred.txt
  nlc    (           1) = 1
  lcindex(1:nlc( 1), 1) = (/ 1/)
  lccoef (1:nlc( 1), 1) = (/ 1.00000/)
  nlc    (           2) = 1
  lcindex(1:nlc( 2), 2) = (/ 2/)
  lccoef (1:nlc( 2), 2) = (/ 1.00000/)
  nlc    (           3) = 1
  lcindex(1:nlc( 3), 3) = (/ 3/)
  lccoef (1:nlc( 3), 3) = (/ 1.00000/)
  nlc    (           4) = 1
  lcindex(1:nlc( 4), 4) = (/ 4/)
  lccoef (1:nlc( 4), 4) = (/ 1.00000/)
  nlc    (           5) = 1
  lcindex(1:nlc( 5), 5) = (/ 5/)
  lccoef (1:nlc( 5), 5) = (/ 1.00000/)
  nlc    (           6) = 1
  lcindex(1:nlc( 6), 6) = (/ 6/)
  lccoef (1:nlc( 6), 6) = (/ 1.00000/)
  nlc    (           7) = 1
  lcindex(1:nlc( 7), 7) = (/ 7/)
  lccoef (1:nlc( 7), 7) = (/ 1.00000/)
  nlc    (           8) = 1
  lcindex(1:nlc( 8), 8) = (/ 8/)
  lccoef (1:nlc( 8), 8) = (/ 1.00000/)
  nlc    (           9) = 1
  lcindex(1:nlc( 9), 9) = (/ 9/)
  lccoef (1:nlc( 9), 9) = (/ 1.00000/)
  nlc    (          10) = 1
  lcindex(1:nlc(10),10) = (/10/)
  lccoef (1:nlc(10),10) = (/ 1.00000/)
  nlc    (          11) = 1
  lcindex(1:nlc(11),11) = (/11/)
  lccoef (1:nlc(11),11) = (/ 1.00000/)
  nlc    (          12) = 1
  lcindex(1:nlc(12),12) = (/12/)
  lccoef (1:nlc(12),12) = (/ 1.00000/)
  nlc    (          13) = 1
  lcindex(1:nlc(13),13) = (/13/)
  lccoef (1:nlc(13),13) = (/ 1.00000/)
  nlc    (          14) = 6
  lcindex(1:nlc(14),14) = (/14, 15, 16, 17, 18, 19/)
  lccoef (1:nlc(14),14) = (/ 0.40825, -0.40825,  0.40825, -0.40825,  0.40825, -0.40825/)
  nlc    (          15) = 6
  lcindex(1:nlc(15),15) = (/14, 15, 16, 17, 18, 19/)
  lccoef (1:nlc(15),15) = (/ 0.57735, -0.28868, -0.28868,  0.57735, -0.28868, -0.28868/)
  nlc    (          16) = 6
  lcindex(1:nlc(16),16) = (/14, 15, 16, 17, 18, 19/)
  lccoef (1:nlc(16),16) = (/ 0.00000,  0.50000, -0.50000,  0.00000,  0.50000, -0.50000/)
  nlc    (          17) = 2
  lcindex(1:nlc(17),17) = (/20, 21/)
  lccoef (1:nlc(17),17) = (/ 0.70711, -0.70711/)
  nlc    (          18) = 2
  lcindex(1:nlc(18),18) = (/22, 23/)
  lccoef (1:nlc(18),18) = (/ 0.70711, -0.70711/)
  nlc    (          19) = 2
  lcindex(1:nlc(19),19) = (/24, 25/)
  lccoef (1:nlc(19),19) = (/ 0.70711, -0.70711/)
  nlc    (          20) = 2
  lcindex(1:nlc(20),20) = (/26, 27/)
  lccoef (1:nlc(20),20) = (/ 0.70711, -0.70711/)
  nlc    (          21) = 2
  lcindex(1:nlc(21),21) = (/28, 29/)
  lccoef (1:nlc(21),21) = (/ 0.70711, -0.70711/)
  nlc    (          22) = 2
  lcindex(1:nlc(22),22) = (/30, 31/)
  lccoef (1:nlc(22),22) = (/ 0.70711, -0.70711/)
  nlc    (          23) = 1
  lcindex(1:nlc(23),23) = (/32/)
  lccoef (1:nlc(23),23) = (/ 1.00000/)
  nlc    (          24) = 6
  lcindex(1:nlc(24),24) = (/33, 34, 35, 36, 37, 38/)
  lccoef (1:nlc(24),24) = (/ 0.40825, -0.40825,  0.40825, -0.40825,  0.40825, -0.40825/)
  nlc    (          25) = 6
  lcindex(1:nlc(25),25) = (/33, 34, 35, 36, 37, 38/)
  lccoef (1:nlc(25),25) = (/-0.28868,  0.57735, -0.28868, -0.28868,  0.57735, -0.28868/)
  nlc    (          26) = 6
  lcindex(1:nlc(26),26) = (/33, 34, 35, 36, 37, 38/)
  lccoef (1:nlc(26),26) = (/ 0.50000,  0.00000, -0.50000,  0.50000,  0.00000, -0.50000/)
  nlc    (          27) = 1
  lcindex(1:nlc(27),27) = (/39/)
  lccoef (1:nlc(27),27) = (/ 1.00000/)
  nlc    (          28) = 1
  lcindex(1:nlc(28),28) = (/40/)
  lccoef (1:nlc(28),28) = (/ 1.00000/)
  nlc    (          29) = 1
  lcindex(1:nlc(29),29) = (/41/)
  lccoef (1:nlc(29),29) = (/ 1.00000/)
  nlc    (          30) = 1
  lcindex(1:nlc(30),30) = (/42/)
  lccoef (1:nlc(30),30) = (/ 1.00000/)
  nlc    (          31) = 1
  lcindex(1:nlc(31),31) = (/43/)
  lccoef (1:nlc(31),31) = (/ 1.00000/)
  nlc    (          32) = 1
  lcindex(1:nlc(32),32) = (/44/)
  lccoef (1:nlc(32),32) = (/ 1.00000/)
  nlc    (          33) = 1
  lcindex(1:nlc(33),33) = (/45/)
  lccoef (1:nlc(33),33) = (/ 1.00000/)
  !--- END generated by python gen_qnr_forUij-t.py intnred.txt

  end subroutine assignqnr

!=================================================================
! *def assignqtc
!    Assign attributes of internal coord. for tert. coupl. (qtc)
!  Input:
!    NA
!  Output:
!    nqtc: number of qtc
!    qtcsym: symmetry (a' or a") of qtc
!    qnrindex: qnr index to which each qtc corresponds
!    qtctyp: type of qtc; used for choosing form of Uij[3] term
!=================================================================
  subroutine assignqtc(nqtc, qtcsym, qnrindex, qtctyp)
  implicit none
  integer :: nqtc, qtcsym(ntermmax), qnrindex(ntermmax) &
       , qtctyp(ntermmax)

  nqtc = 11
  !--- correspondence of qtc to qnr
  qnrindex(1:nqtc) = (/17,23,24,25, 26,28,29,30,31,32,33/)
  !--- type of qtc: 1=normal; 2=angle involving a breaking bond
  qtctyp  (1:nqtc) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
  !--- symm. of qtc at Cs geom. 1=a', 2=a"
  qtcsym(1:nqtc) = (/1,1,1,1, 2,2,2,2,2,2,2/)
  end subroutine assignqtc

!=================================================================
! *def assignparamuij_t
!    Assign param. for tert. couplings
!  Input:
!    nqtc: number of tert. coord. for coupling
!    nk: number of parameters for each anchor point
!    nR0: number of R-type anchor points
!    nphi0: number of phi-type anchor points
!    istate: the 'ij' value of uij of interest
!  Output:
!    k: force coefficients in eV
!    qtc0: rest values
!
!  ~400 lines of params
!=================================================================
  subroutine assignparamuij_t(k, qtc0, nqtc, nk, nR0, nphi0, istate)
  implicit none
  integer :: nqtc, nk, nR0, nphi0, istate
  double precision :: k(nk,nR0,nphi0,nqtc), qtc0(nqtc)
  !--- rest values
  !--- BEGIN generated by "python gen_uijparam.py u12-t_param.txt 2"
  qtc0( 1) =  3.684273e+00
  qtc0( 2) =  9.671927e+01
  qtc0( 3) =  1.443754e-07
  qtc0( 4) = -4.506356e-07
  qtc0( 5) = -4.268868e-07
  qtc0( 6) =  0.000000e+00
  qtc0( 7) =  0.000000e+00
  qtc0( 8) =  0.000000e+00
  qtc0( 9) =  0.000000e+00
  qtc0(10) =  0.000000e+00
  qtc0(11) =  0.000000e+00
  !--- END generated by "python gen_uijparam.py u12-t_param.txt 2"

  !--- add the last AP
  k(:, 3, :, :) =  0.000000e+00
  select case(istate)
  !--- k(K1orK2, iR0, iphi0, iAP)
  !--- BEGIN generated by "python gen_uijk.py tert_k_para.txt"
  case(12)
    k(:, 1, 1, 1) = (/-1.733313e-03, 6.140927e-05/)
    k(:, 1, 1, 2) = (/-1.123851e-03, 1.584981e-05/)
    k(:, 1, 1, 3) = (/ 7.432695e-05,-1.145820e-06/)
    k(:, 1, 1, 4) = (/ 1.934644e-05, 3.490189e-06/)
    k(:, 1, 1, 5) = (/-1.833795e-04,-6.853786e-06/)
    k(:, 1, 1, 6) = (/ 4.604703e-05,-3.963427e-06/)
    k(:, 1, 1, 7) = (/ 4.036499e-05, 4.480713e-06/)
    k(:, 1, 1, 8) = (/-8.452112e-05, 2.014610e-06/)
    k(:, 1, 1, 9) = (/-2.069535e-04,-2.635311e-06/)
    k(:, 1, 1,10) = (/-3.934874e-04, 7.413461e-06/)
    k(:, 1, 1,11) = (/ 1.555130e-04, 1.000857e-05/)
    k(:, 2, 1, 1) = (/-9.200494e-06, 5.684526e-06/)
    k(:, 2, 1, 2) = (/ 2.598292e-07, 4.572700e-06/)
    k(:, 2, 1, 3) = (/-1.064928e-04, 2.805651e-05/)
    k(:, 2, 1, 4) = (/-1.835741e-04, 1.326380e-04/)
    k(:, 2, 1, 5) = (/-1.115714e-05,-1.523882e-05/)
    k(:, 2, 1, 6) = (/ 4.165478e-04, 1.185261e-05/)
    k(:, 2, 1, 7) = (/ 3.905996e-05,-5.531793e-06/)
    k(:, 2, 1, 8) = (/-1.767802e-04, 9.281913e-05/)
    k(:, 2, 1, 9) = (/ 1.237090e-04, 3.900445e-05/)
    k(:, 2, 1,10) = (/-2.318209e-04, 8.995037e-05/)
    k(:, 2, 1,11) = (/ 1.074156e-04,-8.892633e-05/)
    k(:, 1, 2, 1) = (/ 3.550205e-03, 5.415768e-05/)
    k(:, 1, 2, 2) = (/ 7.169084e-03,-1.098199e-04/)
    k(:, 1, 2, 3) = (/-1.931807e-03, 3.703723e-05/)
    k(:, 1, 2, 4) = (/-1.394182e-03, 4.300702e-05/)
    k(:, 1, 2, 5) = (/-9.278868e-04, 4.873925e-05/)
    k(:, 1, 2, 6) = (/ 2.569027e-03, 6.953337e-05/)
    k(:, 1, 2, 7) = (/ 5.827668e-04, 4.437021e-05/)
    k(:, 1, 2, 8) = (/ 1.416050e-04,-2.296609e-06/)
    k(:, 1, 2, 9) = (/ 4.308071e-05, 8.358967e-06/)
    k(:, 1, 2,10) = (/-1.979060e-03,-8.668516e-07/)
    k(:, 1, 2,11) = (/ 2.138308e-04, 1.613456e-05/)
    k(:, 2, 2, 1) = (/-3.925206e-04, 6.034605e-05/)
    k(:, 2, 2, 2) = (/-2.759939e-03, 1.957102e-04/)
    k(:, 2, 2, 3) = (/ 5.102108e-04, 2.403553e-06/)
    k(:, 2, 2, 4) = (/-3.271674e-03, 1.793564e-06/)
    k(:, 2, 2, 5) = (/ 6.092980e-05,-4.068323e-06/)
    k(:, 2, 2, 6) = (/-1.249070e-03, 4.021620e-05/)
    k(:, 2, 2, 7) = (/-7.927348e-05, 2.155305e-05/)
    k(:, 2, 2, 8) = (/ 2.170419e-03,-1.068113e-05/)
    k(:, 2, 2, 9) = (/-6.602038e-04, 6.769210e-06/)
    k(:, 2, 2,10) = (/-1.829907e-03, 8.533659e-06/)
    k(:, 2, 2,11) = (/-1.280150e-03,-1.127389e-05/)
    k(:, 1, 3, 1) = (/ 5.530290e-03,-1.334283e-04/)
    k(:, 1, 3, 2) = (/ 1.079832e-02,-3.475413e-04/)
    k(:, 1, 3, 3) = (/-4.310390e-04,-1.618155e-04/)
    k(:, 1, 3, 4) = (/-3.437741e-04,-1.583457e-04/)
    k(:, 1, 3, 5) = (/ 1.103053e-03,-1.071814e-04/)
    k(:, 1, 3, 6) = (/-5.266252e-04,-6.957921e-05/)
    k(:, 1, 3, 7) = (/-3.278192e-04,-1.663693e-04/)
    k(:, 1, 3, 8) = (/-1.916320e-04,-1.330304e-04/)
    k(:, 1, 3, 9) = (/ 6.961883e-04,-1.099526e-04/)
    k(:, 1, 3,10) = (/-8.202996e-04,-1.699004e-04/)
    k(:, 1, 3,11) = (/-6.685734e-04,-1.460806e-04/)
    k(:, 2, 3, 1) = (/-2.452063e-04, 4.256961e-05/)
    k(:, 2, 3, 2) = (/ 2.540847e-03,-2.777013e-05/)
    k(:, 2, 3, 3) = (/ 6.369245e-04, 2.459653e-05/)
    k(:, 2, 3, 4) = (/-3.206825e-03, 1.696150e-05/)
    k(:, 2, 3, 5) = (/-1.641192e-04, 3.005646e-05/)
    k(:, 2, 3, 6) = (/-2.458476e-03, 6.255257e-05/)
    k(:, 2, 3, 7) = (/-3.014940e-04, 3.662379e-05/)
    k(:, 2, 3, 8) = (/ 2.001268e-03, 2.419227e-05/)
    k(:, 2, 3, 9) = (/-2.960397e-04, 2.428507e-05/)
    k(:, 2, 3,10) = (/-1.527550e-03, 3.101250e-05/)
    k(:, 2, 3,11) = (/-1.551237e-04, 2.426056e-05/)
    k(:, 1, 4, 1) = (/ 2.425588e-03,-1.764648e-05/)
    k(:, 1, 4, 2) = (/ 2.235716e-03,-4.466839e-05/)
    k(:, 1, 4, 3) = (/ 8.355758e-04, 1.917921e-05/)
    k(:, 1, 4, 4) = (/-1.456925e-03,-2.716395e-06/)
    k(:, 1, 4, 5) = (/ 1.309175e-03, 5.699794e-06/)
    k(:, 1, 4, 6) = (/ 3.301770e-03,-8.334269e-05/)
    k(:, 1, 4, 7) = (/-7.397992e-04, 1.332061e-05/)
    k(:, 1, 4, 8) = (/-9.299288e-04, 7.076232e-06/)
    k(:, 1, 4, 9) = (/-2.293984e-04, 4.539787e-06/)
    k(:, 1, 4,10) = (/-1.481681e-03,-2.656398e-06/)
    k(:, 1, 4,11) = (/-3.081512e-04,-2.758696e-05/)
    k(:, 2, 4, 1) = (/-5.662650e-05,-1.311049e-05/)
    k(:, 2, 4, 2) = (/ 1.119997e-03,-3.671780e-05/)
    k(:, 2, 4, 3) = (/ 3.495025e-04,-2.619119e-05/)
    k(:, 2, 4, 4) = (/-3.478137e-04, 9.504889e-05/)
    k(:, 2, 4, 5) = (/ 3.326808e-04,-2.142887e-05/)
    k(:, 2, 4, 6) = (/-2.809782e-03,-9.661015e-06/)
    k(:, 2, 4, 7) = (/-4.692720e-05,-1.711932e-05/)
    k(:, 2, 4, 8) = (/ 1.966424e-03,-1.665941e-05/)
    k(:, 2, 4, 9) = (/-1.768390e-04,-2.929108e-05/)
    k(:, 2, 4,10) = (/-1.649492e-03,-3.313789e-05/)
    k(:, 2, 4,11) = (/ 3.223081e-04,-2.652242e-05/)
    k(:, 1, 5, 1) = (/-7.519468e-05, 3.104521e-05/)
    k(:, 1, 5, 2) = (/ 1.186169e-03, 8.814958e-05/)
    k(:, 1, 5, 3) = (/ 1.944601e-04, 2.037035e-05/)
    k(:, 1, 5, 4) = (/ 4.086014e-04, 9.424137e-06/)
    k(:, 1, 5, 5) = (/ 1.442327e-04, 5.376781e-05/)
    k(:, 1, 5, 6) = (/-4.510180e-04, 1.086084e-04/)
    k(:, 1, 5, 7) = (/ 2.233551e-05, 3.460483e-05/)
    k(:, 1, 5, 8) = (/ 8.393966e-05, 3.275722e-05/)
    k(:, 1, 5, 9) = (/ 6.375066e-05, 2.086896e-05/)
    k(:, 1, 5,10) = (/ 1.023903e-04, 2.758298e-05/)
    k(:, 1, 5,11) = (/-1.322786e-04, 1.868933e-05/)
    k(:, 2, 5, 1) = (/-2.739543e-04,-2.591254e-05/)
    k(:, 2, 5, 2) = (/ 1.644175e-04,-8.241038e-06/)
    k(:, 2, 5, 3) = (/-7.112248e-05,-1.676814e-05/)
    k(:, 2, 5, 4) = (/-1.890528e-04,-1.494142e-04/)
    k(:, 2, 5, 5) = (/ 3.549912e-04,-1.716656e-05/)
    k(:, 2, 5, 6) = (/ 4.193161e-04, 6.918794e-05/)
    k(:, 2, 5, 7) = (/-2.461462e-05, 3.203342e-06/)
    k(:, 2, 5, 8) = (/ 8.062714e-06,-7.998681e-05/)
    k(:, 2, 5, 9) = (/-1.388941e-04,-1.971303e-05/)
    k(:, 2, 5,10) = (/-4.549589e-04,-6.265402e-05/)
    k(:, 2, 5,11) = (/-1.716187e-04,-1.111957e-05/)
  case(13)
    k(:, 1, 1, 1) = (/ 4.233611e-05,-1.440031e-06/)
    k(:, 1, 1, 2) = (/-1.147258e-05, 3.175519e-07/)
    k(:, 1, 1, 3) = (/ 1.223195e-04,-2.188840e-04/)
    k(:, 1, 1, 4) = (/-1.393881e-04,-6.757574e-04/)
    k(:, 1, 1, 5) = (/-6.040445e-05,-2.014182e-04/)
    k(:, 1, 1, 6) = (/-3.002734e-04,-1.968271e-03/)
    k(:, 1, 1, 7) = (/ 1.130252e-04,-5.310176e-04/)
    k(:, 1, 1, 8) = (/ 1.607426e-04,-4.283343e-04/)
    k(:, 1, 1, 9) = (/-1.112961e-05,-5.752479e-05/)
    k(:, 1, 1,10) = (/ 5.549380e-05,-4.514205e-04/)
    k(:, 1, 1,11) = (/-9.626649e-05,-2.946564e-04/)
    k(:, 2, 1, 1) = (/-1.760119e-06, 3.961948e-08/)
    k(:, 2, 1, 2) = (/-1.038482e-05, 6.424664e-07/)
    k(:, 2, 1, 3) = (/-7.166214e-05, 4.045936e-05/)
    k(:, 2, 1, 4) = (/ 3.565803e-05, 7.708932e-05/)
    k(:, 2, 1, 5) = (/ 3.438323e-05, 2.835983e-05/)
    k(:, 2, 1, 6) = (/ 1.130593e-05, 1.475977e-05/)
    k(:, 2, 1, 7) = (/ 1.180181e-04, 5.709408e-05/)
    k(:, 2, 1, 8) = (/ 2.097644e-04, 4.520880e-06/)
    k(:, 2, 1, 9) = (/ 6.315090e-06, 9.911447e-06/)
    k(:, 2, 1,10) = (/ 5.185873e-05, 4.830811e-06/)
    k(:, 2, 1,11) = (/ 6.341090e-05,-1.725945e-05/)
    k(:, 1, 2, 1) = (/-3.037932e-03,-2.059012e-04/)
    k(:, 1, 2, 2) = (/ 5.132848e-03,-3.592990e-04/)
    k(:, 1, 2, 3) = (/ 3.610291e-04,-3.651568e-05/)
    k(:, 1, 2, 4) = (/-6.054996e-04,-4.498317e-05/)
    k(:, 1, 2, 5) = (/-1.686852e-04,-2.751556e-05/)
    k(:, 1, 2, 6) = (/-1.826274e-03,-6.476976e-04/)
    k(:, 1, 2, 7) = (/-1.470543e-03,-2.496202e-05/)
    k(:, 1, 2, 8) = (/ 1.268344e-04,-1.898746e-05/)
    k(:, 1, 2, 9) = (/-1.327626e-04,-1.402519e-05/)
    k(:, 1, 2,10) = (/ 5.597512e-04,-3.403968e-05/)
    k(:, 1, 2,11) = (/ 6.718388e-04,-2.120539e-05/)
    k(:, 2, 2, 1) = (/ 1.039921e-03,-1.462828e-06/)
    k(:, 2, 2, 2) = (/-1.439244e-03, 3.747561e-05/)
    k(:, 2, 2, 3) = (/-5.992639e-04,-1.173022e-05/)
    k(:, 2, 2, 4) = (/ 2.353771e-03,-1.874110e-05/)
    k(:, 2, 2, 5) = (/-7.848634e-04,-4.693980e-06/)
    k(:, 2, 2, 6) = (/ 1.892128e-04,-1.011777e-05/)
    k(:, 2, 2, 7) = (/ 1.479365e-03,-2.278629e-05/)
    k(:, 2, 2, 8) = (/-5.875646e-04,-1.083229e-05/)
    k(:, 2, 2, 9) = (/ 3.636021e-04,-1.266371e-05/)
    k(:, 2, 2,10) = (/ 5.768094e-04,-8.569099e-06/)
    k(:, 2, 2,11) = (/-4.876912e-04,-2.076882e-05/)
    k(:, 1, 3, 1) = (/-1.007583e-02, 1.150186e-04/)
    k(:, 1, 3, 2) = (/ 5.088219e-03, 8.660688e-05/)
    k(:, 1, 3, 3) = (/ 1.052535e-04,-3.004025e-05/)
    k(:, 1, 3, 4) = (/ 4.685570e-03,-2.208504e-05/)
    k(:, 1, 3, 5) = (/ 9.116534e-04,-6.609631e-05/)
    k(:, 1, 3, 6) = (/-2.475618e-02, 2.598621e-04/)
    k(:, 1, 3, 7) = (/-3.682167e-03,-2.373912e-05/)
    k(:, 1, 3, 8) = (/-4.645296e-03,-8.436536e-05/)
    k(:, 1, 3, 9) = (/ 2.361530e-03,-6.735597e-05/)
    k(:, 1, 3,10) = (/ 2.431437e-03,-3.147485e-06/)
    k(:, 1, 3,11) = (/ 3.397998e-03,-3.939173e-05/)
    k(:, 2, 3, 1) = (/ 8.694413e-04,-1.931188e-05/)
    k(:, 2, 3, 2) = (/-1.504769e-03, 5.152748e-05/)
    k(:, 2, 3, 3) = (/-5.132036e-04,-3.239069e-05/)
    k(:, 2, 3, 4) = (/ 3.090883e-03,-2.944202e-05/)
    k(:, 2, 3, 5) = (/-3.289352e-04,-2.679156e-05/)
    k(:, 2, 3, 6) = (/-1.076673e-03,-2.017112e-05/)
    k(:, 2, 3, 7) = (/ 1.225747e-03,-2.805845e-05/)
    k(:, 2, 3, 8) = (/-4.588306e-04,-2.708493e-05/)
    k(:, 2, 3, 9) = (/ 3.819308e-04,-2.248080e-05/)
    k(:, 2, 3,10) = (/ 3.312654e-04,-3.498314e-05/)
    k(:, 2, 3,11) = (/-3.984610e-04,-1.402845e-05/)
    k(:, 1, 4, 1) = (/-1.182391e-02, 4.912747e-05/)
    k(:, 1, 4, 2) = (/ 4.284753e-03,-1.404770e-05/)
    k(:, 1, 4, 3) = (/-7.029479e-05,-6.147952e-05/)
    k(:, 1, 4, 4) = (/ 5.090030e-03,-7.177128e-05/)
    k(:, 1, 4, 5) = (/ 4.229768e-04,-6.127896e-05/)
    k(:, 1, 4, 6) = (/-2.273869e-02, 4.952074e-05/)
    k(:, 1, 4, 7) = (/-2.928830e-03,-9.363366e-05/)
    k(:, 1, 4, 8) = (/-2.320434e-03,-8.314641e-05/)
    k(:, 1, 4, 9) = (/ 8.605732e-06,-6.477149e-05/)
    k(:, 1, 4,10) = (/ 2.422198e-03,-6.245140e-05/)
    k(:, 1, 4,11) = (/ 2.745523e-03,-4.687781e-05/)
    k(:, 2, 4, 1) = (/ 6.521794e-04,-7.603284e-05/)
    k(:, 2, 4, 2) = (/-7.214905e-04,-4.368251e-05/)
    k(:, 2, 4, 3) = (/-6.660444e-04,-8.813874e-05/)
    k(:, 2, 4, 4) = (/ 2.619028e-03,-9.436314e-05/)
    k(:, 2, 4, 5) = (/-4.988018e-04,-7.519505e-05/)
    k(:, 2, 4, 6) = (/-1.964510e-03,-7.596191e-05/)
    k(:, 2, 4, 7) = (/ 1.102339e-03,-8.151012e-05/)
    k(:, 2, 4, 8) = (/-9.307223e-04,-7.123780e-05/)
    k(:, 2, 4, 9) = (/ 2.841984e-04,-9.019627e-05/)
    k(:, 2, 4,10) = (/ 4.839706e-04,-8.099391e-05/)
    k(:, 2, 4,11) = (/-3.007515e-04,-1.048443e-04/)
    k(:, 1, 5, 1) = (/-1.218093e-03,-4.362820e-04/)
    k(:, 1, 5, 2) = (/-4.973274e-04,-2.576429e-06/)
    k(:, 1, 5, 3) = (/-3.224718e-04,-1.618338e-06/)
    k(:, 1, 5, 4) = (/ 5.412665e-04,-1.704953e-04/)
    k(:, 1, 5, 5) = (/-5.517984e-04,-1.078339e-05/)
    k(:, 1, 5, 6) = (/ 3.260079e-03,-8.525028e-04/)
    k(:, 1, 5, 7) = (/-5.883795e-04,-4.080349e-05/)
    k(:, 1, 5, 8) = (/-9.265579e-04,-8.151791e-05/)
    k(:, 1, 5, 9) = (/-4.641611e-05,-1.669791e-05/)
    k(:, 1, 5,10) = (/ 6.471915e-04,-7.850791e-05/)
    k(:, 1, 5,11) = (/ 4.535030e-04,-8.465086e-05/)
    k(:, 2, 5, 1) = (/ 1.260931e-04,-2.870881e-05/)
    k(:, 2, 5, 2) = (/-1.481947e-04, 1.690554e-05/)
    k(:, 2, 5, 3) = (/-3.573728e-04,-2.109928e-05/)
    k(:, 2, 5, 4) = (/ 1.365223e-04,-1.117911e-04/)
    k(:, 2, 5, 5) = (/-2.848717e-05, 2.472967e-06/)
    k(:, 2, 5, 6) = (/ 9.888985e-05,-2.209138e-05/)
    k(:, 2, 5, 7) = (/ 7.557246e-05,-1.097679e-05/)
    k(:, 2, 5, 8) = (/ 3.560006e-05,-1.641905e-05/)
    k(:, 2, 5, 9) = (/-1.228012e-04, 9.355295e-06/)
    k(:, 2, 5,10) = (/ 5.434784e-05,-1.744800e-05/)
    k(:, 2, 5,11) = (/ 4.837745e-05,-4.883250e-07/)
  case(23)
    k(:, 1, 1, 1) = (/-5.618929e-07, 3.648557e-08/)
    k(:, 1, 1, 2) = (/ 9.248460e-08,-4.373233e-09/)
    k(:, 1, 1, 3) = (/-1.563725e-04,-6.110504e-05/)
    k(:, 1, 1, 4) = (/-2.240980e-04,-3.831003e-05/)
    k(:, 1, 1, 5) = (/ 7.337798e-05,-4.587000e-05/)
    k(:, 1, 1, 6) = (/ 2.441892e-04,-1.093111e-05/)
    k(:, 1, 1, 7) = (/ 2.508612e-04,-2.316993e-05/)
    k(:, 1, 1, 8) = (/ 3.360223e-04, 1.761043e-06/)
    k(:, 1, 1, 9) = (/-6.221440e-05, 1.285138e-05/)
    k(:, 1, 1,10) = (/ 9.663138e-06,-7.700481e-05/)
    k(:, 1, 1,11) = (/-6.132585e-05,-1.502945e-06/)
    k(:, 2, 1, 1) = (/ 1.133452e-03, 5.852317e-07/)
    k(:, 2, 1, 2) = (/-5.005573e-04, 5.347888e-05/)
    k(:, 2, 1, 3) = (/-3.944040e-05,-1.179655e-05/)
    k(:, 2, 1, 4) = (/ 2.145334e-05, 8.283277e-05/)
    k(:, 2, 1, 5) = (/-2.139553e-04,-1.566626e-05/)
    k(:, 2, 1, 6) = (/ 1.927716e-04, 5.484677e-05/)
    k(:, 2, 1, 7) = (/ 3.585232e-04, 2.448634e-05/)
    k(:, 2, 1, 8) = (/-7.676862e-05, 4.959070e-05/)
    k(:, 2, 1, 9) = (/ 1.189968e-04, 3.423867e-05/)
    k(:, 2, 1,10) = (/-1.861311e-05, 4.716546e-05/)
    k(:, 2, 1,11) = (/-1.628233e-04, 3.810242e-05/)
    k(:, 1, 2, 1) = (/-2.072190e-03,-5.199952e-05/)
    k(:, 1, 2, 2) = (/-1.034962e-04, 1.262660e-04/)
    k(:, 1, 2, 3) = (/-9.623645e-04,-1.696602e-05/)
    k(:, 1, 2, 4) = (/-4.786655e-04,-2.472436e-05/)
    k(:, 1, 2, 5) = (/-7.515358e-04,-2.344332e-05/)
    k(:, 1, 2, 6) = (/-2.244887e-03,-6.052707e-05/)
    k(:, 1, 2, 7) = (/ 1.562763e-04,-1.294151e-05/)
    k(:, 1, 2, 8) = (/-6.664864e-04,-2.484866e-05/)
    k(:, 1, 2, 9) = (/ 2.565785e-05,-1.045594e-05/)
    k(:, 1, 2,10) = (/-3.121890e-04,-9.000814e-06/)
    k(:, 1, 2,11) = (/ 6.585018e-05,-4.316256e-05/)
    k(:, 2, 2, 1) = (/ 8.104219e-04,-2.411988e-05/)
    k(:, 2, 2, 2) = (/ 2.665220e-05, 6.163214e-06/)
    k(:, 2, 2, 3) = (/ 8.722384e-04, 8.834701e-05/)
    k(:, 2, 2, 4) = (/ 1.360684e-03,-6.597433e-06/)
    k(:, 2, 2, 5) = (/ 8.228824e-05, 3.409907e-05/)
    k(:, 2, 2, 6) = (/ 1.346192e-03,-2.470937e-05/)
    k(:, 2, 2, 7) = (/-1.159925e-03, 5.913071e-05/)
    k(:, 2, 2, 8) = (/ 2.777742e-05,-4.529754e-05/)
    k(:, 2, 2, 9) = (/ 2.685689e-04,-2.182174e-05/)
    k(:, 2, 2,10) = (/ 3.884088e-04,-1.786362e-05/)
    k(:, 2, 2,11) = (/ 6.607354e-04, 4.594803e-05/)
    k(:, 1, 3, 1) = (/ 1.898945e-03,-4.114615e-05/)
    k(:, 1, 3, 2) = (/ 2.124939e-03,-1.062008e-04/)
    k(:, 1, 3, 3) = (/-2.136585e-04,-3.471270e-05/)
    k(:, 1, 3, 4) = (/ 3.304956e-04,-1.038134e-05/)
    k(:, 1, 3, 5) = (/ 3.302591e-04,-3.584463e-05/)
    k(:, 1, 3, 6) = (/ 1.872065e-03,-7.902102e-05/)
    k(:, 1, 3, 7) = (/-2.890172e-04,-7.647184e-05/)
    k(:, 1, 3, 8) = (/ 1.526834e-04,-2.790347e-05/)
    k(:, 1, 3, 9) = (/ 7.562355e-04,-1.988505e-05/)
    k(:, 1, 3,10) = (/-7.584423e-04,-3.636800e-05/)
    k(:, 1, 3,11) = (/ 1.776274e-04,-3.672316e-05/)
    k(:, 2, 3, 1) = (/ 2.197848e-03,-3.207625e-05/)
    k(:, 2, 3, 2) = (/ 7.075983e-05, 5.283304e-06/)
    k(:, 2, 3, 3) = (/ 1.940724e-04, 7.353277e-05/)
    k(:, 2, 3, 4) = (/ 1.104931e-03,-5.173726e-06/)
    k(:, 2, 3, 5) = (/ 3.425042e-04, 4.981673e-05/)
    k(:, 2, 3, 6) = (/ 3.407230e-04, 4.376735e-06/)
    k(:, 2, 3, 7) = (/-2.256775e-04, 2.938070e-05/)
    k(:, 2, 3, 8) = (/-1.471331e-04,-3.947064e-05/)
    k(:, 2, 3, 9) = (/-3.516193e-04,-2.389501e-05/)
    k(:, 2, 3,10) = (/ 2.010860e-04, 2.274867e-06/)
    k(:, 2, 3,11) = (/ 1.323799e-04, 4.627243e-06/)
    k(:, 1, 4, 1) = (/-5.433620e-04,-3.498279e-05/)
    k(:, 1, 4, 2) = (/ 1.370427e-03,-8.051840e-05/)
    k(:, 1, 4, 3) = (/ 8.402759e-04, 2.069944e-05/)
    k(:, 1, 4, 4) = (/ 2.342223e-04, 1.203957e-05/)
    k(:, 1, 4, 5) = (/-4.406597e-04,-5.395087e-05/)
    k(:, 1, 4, 6) = (/ 2.929828e-03, 3.333899e-05/)
    k(:, 1, 4, 7) = (/-4.453296e-04,-1.548554e-05/)
    k(:, 1, 4, 8) = (/ 6.919792e-05, 3.127016e-05/)
    k(:, 1, 4, 9) = (/ 3.283214e-06, 2.286256e-05/)
    k(:, 1, 4,10) = (/ 2.432282e-04,-1.301495e-05/)
    k(:, 1, 4,11) = (/-9.170427e-05, 4.179726e-05/)
    k(:, 2, 4, 1) = (/ 4.142901e-04, 1.322291e-05/)
    k(:, 2, 4, 2) = (/ 1.086775e-04, 7.772806e-06/)
    k(:, 2, 4, 3) = (/ 2.677530e-04, 4.101547e-05/)
    k(:, 2, 4, 4) = (/ 4.568884e-04,-3.686458e-06/)
    k(:, 2, 4, 5) = (/-1.292690e-04, 2.361311e-05/)
    k(:, 2, 4, 6) = (/ 4.734334e-04, 1.039946e-05/)
    k(:, 2, 4, 7) = (/-1.983331e-04, 1.311048e-05/)
    k(:, 2, 4, 8) = (/ 2.328728e-04,-3.049184e-05/)
    k(:, 2, 4, 9) = (/-2.109337e-04,-1.693929e-05/)
    k(:, 2, 4,10) = (/ 1.117041e-04,-3.389942e-05/)
    k(:, 2, 4,11) = (/-3.681104e-04,-1.271709e-05/)
    k(:, 1, 5, 1) = (/ 1.000661e-04, 4.864518e-06/)
    k(:, 1, 5, 2) = (/ 9.362372e-05, 4.056004e-06/)
    k(:, 1, 5, 3) = (/ 1.027685e-04, 2.978778e-06/)
    k(:, 1, 5, 4) = (/ 6.354179e-06,-1.034004e-05/)
    k(:, 1, 5, 5) = (/ 9.528239e-05, 2.052767e-05/)
    k(:, 1, 5, 6) = (/-6.196094e-04,-1.227168e-04/)
    k(:, 1, 5, 7) = (/ 1.595500e-04,-2.537411e-05/)
    k(:, 1, 5, 8) = (/-1.253366e-04, 1.240045e-05/)
    k(:, 1, 5, 9) = (/-2.212461e-04, 1.553848e-05/)
    k(:, 1, 5,10) = (/-2.613636e-05, 1.216166e-05/)
    k(:, 1, 5,11) = (/-1.447876e-04,-8.188116e-06/)
    k(:, 2, 5, 1) = (/ 4.572047e-04,-2.947590e-06/)
    k(:, 2, 5, 2) = (/-4.324040e-05,-4.776077e-05/)
    k(:, 2, 5, 3) = (/-3.215957e-05,-5.019978e-05/)
    k(:, 2, 5, 4) = (/-1.605995e-04, 2.240962e-06/)
    k(:, 2, 5, 5) = (/ 3.184796e-05,-3.594095e-05/)
    k(:, 2, 5, 6) = (/ 3.854867e-04,-6.698726e-06/)
    k(:, 2, 5, 7) = (/ 9.213132e-05,-5.696104e-06/)
    k(:, 2, 5, 8) = (/ 5.087784e-04, 2.928541e-05/)
    k(:, 2, 5, 9) = (/ 3.860606e-05, 1.301820e-05/)
    k(:, 2, 5,10) = (/ 2.928029e-04, 1.929021e-05/)
    k(:, 2, 5,11) = (/-6.358518e-05,-5.136726e-06/)
  !--- END generated by "python gen_uijparam.py tert_k_para.txt"
  end select

  end subroutine assignparamuij_t

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!             Useful subroutines used before
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!=================================================================
! *def evalBL
!    Get bond length of atoms i and j.
!  Input:
!    x: Cartesian coordinates in Angstroms
!    i, j: atom indices
!  Return:
!    bond length in Angstroms
!=================================================================
  double precision function evalBL(x, i, j)
  implicit none
  integer :: i, j
  double precision :: x(3,natom)
  double precision :: r(3)

  r(:) = x(:,i)-x(:,j)
  evalBL = sqrt(dot_product(r, r))
  end function evalBL

!=================================================================
! *def evalBA
!    Get bond angle of atoms i ,j, k.
!  Input:
!    x: Cartesian coordinates in Angstroms
!    i, j, k: atom indices
!  Return:
!    bond angle in radians
!=================================================================
  double precision function evalBA(x, i, j, k)
  implicit none
  integer :: i, j, k
  double precision :: x(3,natom)
  double precision :: eji(3), ejk(3), cosBA

  eji(:) = (x(:,i)-x(:,j))/evalBL(x, i, j)
  ejk(:) = (x(:,k)-x(:,j))/evalBL(x, j, k)
  cosBA = dot_product(eji, ejk)
  !--- avoid numerical problem
  if(cosBA>1) then
     evalBA = 0d0
  else if(cosBA<-1) then
     evalBA = pi
  else
     evalBA = acos(cosBA)
  end if
  end function evalBA

!=================================================================
! *def evalTO
!    Get torsion of atoms i ,j, k, l.
!  Input:
!    x: Cartesian coordinates in Angstroms
!    i, j, k, l: atom indices
!  Return:
!    torsion in radians
!=================================================================
  double precision function evalTO(x, i, j, k, l)
  implicit none
  integer :: i, j, k, l
  double precision :: x(3,natom)
  integer :: sign
  double precision :: eij(3), ejk(3), ekl(3) &
       , sinijk, sinjkl, cosTO &
       , cross1(3), cross2(3), cross3(3)

  eij(:) = (x(:,j)-x(:,i))/evalBL(x, i, j)
  ejk(:) = (x(:,k)-x(:,j))/evalBL(x, j, k)
  ekl(:) = (x(:,l)-x(:,k))/evalBL(x, k, l)
  sinijk = sin(evalBA(x, i, j, k))
  sinjkl = sin(evalBA(x, j, k, l))
  call xprod(cross1, eij, ejk)
  call xprod(cross2, ejk, ekl)
  cosTO  = dot_product(cross1, cross2)/sinijk/sinjkl
  !--- determine sign of torsion
  call xprod(cross1, eij, ejk)
  call xprod(cross2, ejk, ekl)
  call xprod(cross3, cross1, cross2)
  if(dot_product(cross3, ejk)>0) then
     sign = 1
  else
     sign = -1
  end if
  !--- avoid numerical problem
  if(cosTO>1) then
     evalTO = 0d0
  else if(cosTO<-1) then
     evalTO = pi
  else
     evalTO = sign*acos(cosTO)
  end if
  end function evalTO

!=================================================================
! *def evalOB
!    Get out-of-plane bend of atoms i ,j, k, l.
!  Input:
!    x: Cartesian coordinates in Angstroms
!    i, j, k, l: atom indices
!  Return:
!    oop bend in radians
!=================================================================
  double precision function evalOB(x, i, j, k, l)
  implicit none
  integer :: i, j, k, l
  double precision :: x(3,natom)
  double precision :: eli(3), elj(3), elk(3), sinjlk, sinOB &
       , cross1(3)

  eli(:) = (x(:,i)-x(:,l))/evalBL(x, l, i)
  elj(:) = (x(:,j)-x(:,l))/evalBL(x, l, j)
  elk(:) = (x(:,k)-x(:,l))/evalBL(x, l, k)
  sinjlk = sin(evalBA(x, j, l, k))
  call xprod(cross1, elj, elk)
  sinOB  = dot_product(cross1, eli)/sinjlk
  !--- avoid numerical problem
  if(sinOB>1) then
     evalOB = pi/2
  else if(sinOB<-1) then
     evalOB = -pi/2
  else
     evalOB = asin(sinOB)
  end if
  end function evalOB

!=================================================================
! *def evalOD
!    Get out-of-plane distance of atoms i ,j, k, l.
!  Input:
!    x: Cartesian coordinates in Angstroms
!    i, j, k, l: atom indices
!  Return:
!    oop distance in Angstroms
!=================================================================
  double precision function evalOD(x, i, j, k, l)
  implicit none
  integer :: i, j, k, l
  double precision :: x(3,natom)
  double precision :: eij(3), eik(3), eijxeik(3), ril(3), sinjik

  ril(:) = x(:,l)-x(:,i)
  eij(:) = (x(:,j)-x(:,i))/evalBL(x, i, j)
  eik(:) = (x(:,k)-x(:,i))/evalBL(x, i, k)
  call xprod(eijxeik, eij, eik)
  sinjik = sin(evalBA(x, j, i, k))
  evalOD = abs(dot_product(ril, eijxeik)/sinjik)
  end function evalOD

!=================================================================
! *def evalSP1
!    Get a special coordinate to replace torsion CCSH
!  Input:
!    x: Cartesian coordinates in Angstroms
!    i, j, k, l, m: atom indices
!  Return:
!    The coordinate
!=================================================================
  double precision function evalSP1(x, i, j, k, l, m)
  implicit none
  integer :: i, j, k, l, m
  double precision :: x(3,natom)
  double precision :: e12(3), e23(3), e13(3), eZ(3), eX(3), ep(3) &
       , e45(3), e24(3), etmp(3), etmp2(3), etmp3(3) &
       , sin123, sin13Z, cos13Z, e12xe23(3), cosSP1, sgn

  e12(:) = (x(:,j)-x(:,i))/evalBL(x, i, j)
  e23(:) = (x(:,k)-x(:,j))/evalBL(x, j, k)
  e13(:) = (x(:,k)-x(:,i))/evalBL(x, i, k)
  e45(:) = (x(:,m)-x(:,l))/evalBL(x, l, m)
  e24(:) = (x(:,l)-x(:,j))/evalBL(x, j, l)
  call xprod(etmp, e12, e23)
  eZ(:) = etmp(:)/sqrt(dot_product(etmp, etmp))
  call xprod(etmp, e13, eZ)
  eX(:) = etmp(:)/sqrt(dot_product(etmp, etmp))
  etmp(:) = e45(:)-eX(:)*dot_product(e45, eX)
  ep(:) = etmp(:)/sqrt(dot_product(etmp, etmp))
  cosSP1 = dot_product(ep, e13)

  !--- determine sign
  call xprod(etmp, e13, ep)
  if(dot_product(etmp, eX)>0) then
     sgn = 1
  else
     sgn = -1
  end if
  !--- avoid numerical problem
  if(cosSP1>1) then
     evalSP1 = 0d0
  else if(cosSP1<-1) then
     evalSP1 = pi
  else
     evalSP1 = sgn*acos(cosSP1)
  end if
  end function evalSP1

!=================================================================
! *def xprod
!    Cross product of two 3D vectors.
!  Input:
!    a(3), b(3): two vectors
!  Return:
!    axb(3): = a x b
!=================================================================
  subroutine xprod(axb, a, b)
  implicit none
  double precision :: axb(3), a(3), b(3)

  axb(1) = a(2)*b(3)-a(3)*b(2)
  axb(2) = a(3)*b(1)-a(1)*b(3)
  axb(3) = a(1)*b(2)-a(2)*b(1)
  end subroutine xprod

!=================================================================
! *def calctent1
!    Calculate tent func. and its grad. (for interpolating R)
!  Input:
!    q: internal coordinate as argument
!    q0: ref. q values at anchor points
!    i, j, k: anchor point indexes
!    nq0: number of anchor points
!  Output:
!    t: tent function value (unitless)
!    dtdq: dt/dq
!=================================================================
  subroutine calctent1(t, dtdq, q, q0, i, j, k, nq0)
  implicit none
  double precision :: t, dtdq, q, q0(:)
  integer :: i, j, k, nq0

  !--- if j is the 1st anchor point
  if(i.eq.0) then
     if(q<q0(j)) then
        t = 1d0; dtdq = 0d0
     else if(q>q0(k)) then
        t = 0d0; dtdq = 0d0
     else
        call tent1(t, dtdq, q, q0(k), q0(j))
     end if
  !--- if j is the last anchor point
  else if(k>nq0) then
     if(q>q0(j)) then
        t = 1d0; dtdq = 0d0
     else if(q<q0(i)) then
        t = 0d0; dtdq = 0d0
     else
        call tent1(t, dtdq, q, q0(i), q0(j))
     end if
  else
     if(q<q0(i).or.q>q0(k)) then
        t = 0d0; dtdq = 0d0
     else if(q<q0(j)) then
        call tent1(t, dtdq, q, q0(i), q0(j))
     else
        call tent1(t, dtdq, q, q0(k), q0(j))
     end if
  end if

  end subroutine calctent1
  !--- called by calctent1
  !--- tent(q1) = 0, tent(q2) = 1
  subroutine tent1(t, dtdq, q, q1, q2)
    implicit none
    double precision :: t, dtdq, q, q1, q2
    double precision :: dq1, dq2
    dq1 = q-q1
    dq2 = q-q2
    t = dq1**4/(dq1**4+dq2**4)
    dtdq = ( 4*dq1**3*(dq1**4+dq2**4)-4*dq1**4*(dq1**3+dq2**3) ) &
         / (dq1**4+dq2**4)**2
  end subroutine tent1

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

!=================================================================
! *def calcbmat
!    Calculate a row of B matrix for an internal coordinate.
!  Input:
!    x: Cartesian coordinates in Angstroms
!    iatomlist: indexes of atoms involved in the internal coord.
!    inttype: type of the internal coordinate
!  Output:
!    bmatrow: row of B matrix
!=================================================================
  subroutine calcbmat(bmatrow, x, iatomlist, inttype)
  implicit none
  double precision :: bmatrow(natom*3), x(3,natom)
  integer :: inttype, iatomlist(10)
  double precision :: rij, rij1, rij2, rij3, rjk1, rjk2, rjk3, rjk, pijk
  double precision :: cijk, sijk, pxmg, pymg, pzmg, sjkl, sijkl
  double precision :: rkl, rkl1, rkl2, rkl3, cjkl, cijkl, pjkl
  integer :: i, j, ijk, im, k, l, m, ijkl, ilisttmp(4), ii, iii, jj
  integer :: nn_iatomlist(10)
  double precision :: A243,S243,C243,E41(3),E42(3),E43(3) &
       ,R41,R42,R43,OOPB,COOPB,TOOPB,TMPV(3),BMV1(3),BMV2(3),BMV3(3) &
       ,BMV4(3),e12(3),e13(3),v14(3),r12,r13,theta,rtmp &
       ,tmpv2(3),tmpv3(3),rowtmp(3*natom),xtmp(3,natom) &
       ,coord1,coord2,step

  bmatrow(:) = 0d0
  select case(inttype)
  !--- bond length; adapted from polyrate
  case(1)
     I=iatomlist(1)
     J=iatomlist(2)
     RIJ1=X(1,J)-X(1,I)
     RIJ2=x(2,J)-x(2,I)
     RIJ3=x(3,J)-x(3,I)

     RIJ=evalBL(x, i, j)

     BMATROW(3*I-2)=-RIJ1/RIJ
     BMATROW(3*I-1)=-RIJ2/RIJ
     BMATROW(3*I)=-RIJ3/RIJ
     BMATROW(3*J-2)=RIJ1/RIJ
     BMATROW(3*J-1)=RIJ2/RIJ
     BMATROW(3*J)=RIJ3/RIJ
  !--- bending; adapted from polyrate
  case(2)
     I=iatomlist(1)
     J=iatomlist(2)
     K=iatomlist(3)

     RIJ1=X(1,J)-X(1,I)
     RIJ2=X(2,J)-X(2,I)
     RIJ3=x(3,J)-x(3,I)
     RJK1=X(1,K)-X(1,J)
     RJK2=X(2,K)-X(2,J)
     RJK3=x(3,K)-x(3,J)
     RIJ=evalBL(x, I, J)
     RJK=evalBL(x, J, K)
     PIJK=evalBA(x, I, J, K)
     CIJK=cos(PIJK)
     SIJK=sin(PIJK)

     BMATROW(3*I-2)=(-CIJK*RIJ1*RJK-RIJ*RJK1)/(SIJK*RIJ**2*RJK)
     BMATROW(3*I-1)=(-CIJK*RIJ2*RJK-RIJ*RJK2)/(SIJK*RIJ**2*RJK)
     BMATROW(3*I)=(-CIJK*RIJ3*RJK-RIJ*RJK3)/(SIJK*RIJ**2*RJK)
     BMATROW(3*J-2)= (-RIJ*RIJ1*RJK+CIJK*RIJ1*RJK**2-CIJK*RIJ**2*RJK1+RIJ*RJK*RJK1) &
          /(SIJK*RIJ**2*RJK**2)
     BMATROW(3*J-1)= (-RIJ*RIJ2*RJK+CIJK*RIJ2*RJK**2-CIJK*RIJ**2*RJK2+RIJ*RJK*RJK2) &
          /(SIJK*RIJ**2*RJK**2)
     BMATROW(3*J)= (-RIJ*RIJ3*RJK+CIJK*RIJ3*RJK**2-CIJK*RIJ**2*RJK3+RIJ*RJK*RJK3)  &
          /(SIJK*RIJ**2*RJK**2)
     BMATROW(3*K-2)=(RIJ1*RJK+CIJK*RIJ*RJK1)/(SIJK*RIJ*RJK**2)
     BMATROW(3*K-1)=(RIJ2*RJK+CIJK*RIJ*RJK2)/(SIJK*RIJ*RJK**2)
     BMATROW(3*K)=(RIJ3*RJK+CIJK*RIJ*RJK3)/(SIJK*RIJ*RJK**2)

  !--- torsion; adpated from polyrate
  case(3)
     I=IATOMLIST(1)
     J=IATOMLIST(2)
     K=IATOMLIST(3)
     L=IATOMLIST(4)

     RIJ1=X(1,J)-X(1,I)
     RIJ2=X(2,J)-X(2,I)
     RIJ3=X(3,J)-X(3,I)
     RJK1=X(1,K)-X(1,J)
     RJK2=X(2,K)-X(2,J)
     RJK3=X(3,K)-X(3,J)
     RKL1=X(1,L)-X(1,K)
     RKL2=X(2,L)-X(2,K)
     RKL3=X(3,L)-X(3,K)
     RIJ=evalBL(x, I,J)
     RJK=evalBL(x, J,K)
     RKL=evalBL(x, K,L)

     PIJK=evalBA(x, I,J,K)
     CIJK=COS(PIJK)
     SIJK=SIN(PIJK)
     PJKL=evalBA(x, J,K,L)
     CJKL=COS(PJKL)
     SJKL=SIN(PJKL)
     CIJKL=((-RIJ2*RJK1+RIJ1*RJK2)*(-RJK2*RKL1+RJK1*RKL2)+ &
          (RIJ3*RJK1-RIJ1*RJK3)*(RJK3*RKL1-RJK1*RKL3)+     &
          (-RIJ3*RJK2+RIJ2*RJK3)*(-RJK3*RKL2+RJK2*RKL3))/  &
          (SIJK*SJKL*RIJ*RJK*RJK*RKL)
     SIJKL=((-RIJ3*RJK2+RIJ2*RJK3)*RKL1+(RIJ3*RJK1-RIJ1*RJK3)*RKL2+  &
          (-(RIJ2*RJK1)+RIJ1*RJK2)*RKL3)/(RIJ*RJK*RKL*SIJK*SJKL)
     BMATROW(3*I-2)=SIJKL*RIJ1/(CIJKL*SIJK**2*RIJ**2)+ &
          CIJK*SIJKL*RJK1/(CIJKL*SIJK**2*RIJ*RJK)+     &
          (RJK3*RKL2-RJK2*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
     BMATROW(3*I-1)=SIJKL*RIJ2/(CIJKL*SIJK**2*RIJ**2)+ &
          CIJK*SIJKL*RJK2/(CIJKL*SIJK**2*RIJ*RJK)+     &
          (-RJK3*RKL1+RJK1*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
     BMATROW(3*I)=SIJKL*RIJ3/(CIJKL*SIJK**2*RIJ**2)+ &
          CIJK*SIJKL*RJK3/(CIJKL*SIJK**2*RIJ*RJK)+   &
          (RJK2*RKL1-RJK1*RKL2)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
     BMATROW(3*J-2)=-(SIJKL*RIJ1/(CIJKL*SIJK**2*RIJ**2))+  &
          CIJK*SIJKL*(RIJ1-RJK1)/(CIJKL*SIJK**2*RIJ*RJK)-  &
          SIJKL*RJK1/(CIJKL*RJK**2)+SIJKL*RJK1/(CIJKL*SIJK**2*RJK**2)+  &
          SIJKL*RJK1/(CIJKL*SJKL**2*RJK**2)+  &
          CJKL*SIJKL*RKL1/(CIJKL*SJKL**2*RJK*RKL)+  &
          (-RIJ3*RKL2-RJK3*RKL2+RIJ2*RKL3+RJK2*RKL3)/  &
          (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
     BMATROW(3*J-1)=-SIJKL*RIJ2/(CIJKL*SIJK**2*RIJ**2)+  &
          CIJK*SIJKL*(RIJ2-RJK2)/(CIJKL*SIJK**2*RIJ*RJK)- &
          SIJKL*RJK2/(CIJKL*RJK**2)+SIJKL*RJK2/(CIJKL*SIJK**2*RJK**2)+  &
          SIJKL*RJK2/(CIJKL*SJKL**2*RJK**2)+ &
          CJKL*SIJKL*RKL2/(CIJKL*SJKL**2*RJK*RKL)+  &
          (RIJ3*RKL1+RJK3*RKL1-RIJ1*RKL3-RJK1*RKL3)/  &
          (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
     BMATROW(3*J)=-SIJKL*RIJ3/(CIJKL*SIJK**2*RIJ**2)+   &
          CIJK*SIJKL*(RIJ3-RJK3)/(CIJKL*SIJK**2*RIJ*RJK)-  &
          SIJKL*RJK3/(CIJKL*RJK**2)+SIJKL*RJK3/(CIJKL*SIJK**2*RJK**2)+  &
          SIJKL*RJK3/(CIJKL*SJKL**2*RJK**2)+   &
          (-RIJ2*RKL1-RJK2*RKL1+RIJ1*RKL2+RJK1*RKL2)/  &
          (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)+   &
          CJKL*SIJKL*RKL3/(CIJKL*SJKL**2*RJK*RKL)
     BMATROW(3*K-2)=-CIJK*SIJKL*RIJ1/(CIJKL*SIJK**2*RIJ*RJK)+  &
          SIJKL*RJK1/(CIJKL*RJK**2)-SIJKL*RJK1/(CIJKL*SIJK**2*RJK**2)-  &
          SIJKL*RJK1/(CIJKL*SJKL**2*RJK**2)+  &
          CJKL*SIJKL*(RJK1-RKL1)/(CIJKL*SJKL**2*RJK*RKL)+  &
          SIJKL*RKL1/(CIJKL*SJKL**2*RKL**2)+  &
          (RIJ3*RJK2-RIJ2*RJK3+RIJ3*RKL2-RIJ2*RKL3)/  &
          (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
     BMATROW(3*K-1)=-CIJK*SIJKL*RIJ2/(CIJKL*SIJK**2*RIJ*RJK)+  &
          SIJKL*RJK2/(CIJKL*RJK**2)-SIJKL*RJK2/(CIJKL*SIJK**2*RJK**2)-  &
          SIJKL*RJK2/(CIJKL*SJKL**2*RJK**2)+  &
          CJKL*SIJKL*(RJK2-RKL2)/(CIJKL*SJKL**2*RJK*RKL)+  &
          SIJKL*RKL2/(CIJKL*SJKL**2*RKL**2)+  &
          (-RIJ3*RJK1+RIJ1*RJK3-RIJ3*RKL1+RIJ1*RKL3)/  &
          (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)
     BMATROW(3*K)=-CIJK*SIJKL*RIJ3/(CIJKL*SIJK**2*RIJ*RJK)+  &
          SIJKL*RJK3/(CIJKL*RJK**2)-SIJKL*RJK3/(CIJKL*SIJK**2*RJK**2)-  &
          SIJKL*RJK3/(CIJKL*SJKL**2*RJK**2)+  &
          (RIJ2*RJK1-RIJ1*RJK2+RIJ2*RKL1-RIJ1*RKL2)/  &
          (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)+  &
          CJKL*SIJKL*(RJK3-RKL3)/(CIJKL*SJKL**2*RJK*RKL)+  &
          SIJKL*RKL3/(CIJKL*SJKL**2*RKL**2)
     BMATROW(3*L-2)=-CJKL*SIJKL*RJK1/(CIJKL*SJKL**2*RJK*RKL)+ &
          (-RIJ3*RJK2+RIJ2*RJK3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)- &
          SIJKL*RKL1/(CIJKL*SJKL**2*RKL**2)
     BMATROW(3*L-1)=-CJKL*SIJKL*RJK2/(CIJKL*SJKL**2*RJK*RKL)+  &
          (RIJ3*RJK1-RIJ1*RJK3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)- &
          SIJKL*RKL2/(CIJKL*SJKL**2*RKL**2)
     BMATROW(3*L)=(-RIJ2*RJK1+RIJ1*RJK2)/  &
          (CIJKL*SIJK*SJKL*RIJ*RJK*RKL)-   &
          CJKL*SIJKL*RJK3/(CIJKL*SJKL**2*RJK*RKL)-  &
          SIJKL*RKL3/(CIJKL*SJKL**2*RKL**2)

  !--- out-of-plane bend
  case(4)
     I=IATOMLIST(1)
     J=IATOMLIST(2)
     K=IATOMLIST(3)
     L=IATOMLIST(4)
     e41(1)=X(1,I)-X(1,L)
     e41(2)=X(2,I)-X(2,L)
     e41(3)=X(3,I)-X(3,L)
     e42(1)=X(1,j)-X(1,L)
     e42(2)=X(2,j)-X(2,L)
     e42(3)=X(3,j)-X(3,L)
     e43(1)=X(1,k)-X(1,L)
     e43(2)=X(2,k)-X(2,L)
     e43(3)=X(3,k)-X(3,L)
     R41=evalBL(x, I,L)
     R42=evalBL(x, J,L)
     R43=evalBL(x, K,L)
     e41=e41/R41
     e42=e42/R42
     e43=e43/r43
     a243=evalBA(x, j,l,k)
     c243=cos(a243)
     s243=sin(a243)
     call xprod(tmpv,e42,e43)
     oopb=asin(dot_product(tmpv,e41)/s243)
     coopb=cos(oopb)
     toopb=tan(oopb)

     call xprod(tmpv,e42,e43)
     bmv1=(tmpv/coopb/s243-toopb*e41)/r41
     bmatrow(3*I-2)=bmv1(1)
     bmatrow(3*i-1)=bmv1(2)
     bmatrow(3*i)=bmv1(3)
     call xprod(tmpv,e43,e41)
     bmv2=(tmpv/coopb/s243-toopb/s243**2*(e42-e43*c243))/r42
     bmatrow(3*j-2)=bmv2(1)
     bmatrow(3*j-1)=bmv2(2)
     bmatrow(3*j)=bmv2(3)
     call xprod(tmpv,e41,e42)
     bmv3=(tmpv/coopb/s243-toopb/s243**2*(e43-e42*c243))/r43
     bmatrow(3*k-2)=bmv3(1)
     bmatrow(3*k-1)=bmv3(2)
     bmatrow(3*k)=bmv3(3)
     bmv4=-bmv1-bmv2-bmv3
     bmatrow(3*l-2)=bmv4(1)
     bmatrow(3*l-1)=bmv4(2)
     bmatrow(3*l)=bmv4(3)

  !--- out-of-plane distance (aka pyramid height)
  case(5)
     I=IATOMLIST(1)
     J=IATOMLIST(2)
     K=IATOMLIST(3)
     L=IATOMLIST(4)
     e12(1)=X(1,j)-X(1,i)
     e12(2)=X(2,j)-X(2,i)
     e12(3)=X(3,j)-X(3,i)
     e13(1)=X(1,k)-X(1,i)
     e13(2)=X(2,k)-X(2,i)
     e13(3)=X(3,k)-X(3,i)
     v14(1)=X(1,l)-X(1,i)
     v14(2)=X(2,l)-X(2,i)
     v14(3)=X(3,l)-X(3,i)
     r12=evalBL(x, i,j)
     r13=evalBL(x, i,k)
     e12(:)=e12(:)/r12
     e13(:)=e13(:)/r13
     theta=evalBA(x, j, i, k)
     call xprod(tmpv2, e12, e13)
     rtmp=dot_product(v14, tmpv2)

     call xprod(tmpv, e13, v14)
     tmpv3(:)=(cos(theta)*e12(:)-e13(:))/r12/sin(theta)
     bmv2 = (tmpv-rtmp*e12)/sin(theta)/r12 &
          - cos(theta)*rtmp/sin(theta)**2*tmpv3

     call xprod(tmpv, v14, e12)
     tmpv3(:)=(cos(theta)*e13(:)-e12(:))/r13/sin(theta)
     bmv3 = (tmpv-rtmp*e13)/sin(theta)/r13 &
          - cos(theta)*rtmp/sin(theta)**2*tmpv3

     call xprod(tmpv, e12, e13)
     bmv4 = tmpv/sin(theta)

     bmv1 = -bmv2-bmv3-bmv4

     if(rtmp/sin(theta)<0) then
        bmv1 = -bmv1
        bmv2 = -bmv2
        bmv3 = -bmv3
        bmv4 = -bmv4
     end if

     bmatrow(3*I-2)=bmv1(1)
     bmatrow(3*i-1)=bmv1(2)
     bmatrow(3*i)=bmv1(3)
     bmatrow(3*j-2)=bmv2(1)
     bmatrow(3*j-1)=bmv2(2)
     bmatrow(3*j)=bmv2(3)
     bmatrow(3*k-2)=bmv3(1)
     bmatrow(3*k-1)=bmv3(2)
     bmatrow(3*k)=bmv3(3)
     bmatrow(3*l-2)=bmv4(1)
     bmatrow(3*l-1)=bmv4(2)
     bmatrow(3*l)=bmv4(3)

  !--- special coordinate replacing torsion CCSH
  !--- numerical gradient
  case(6)
     step = 1d-5
     i = iatomlist(1)
     j = iatomlist(2)
     k = iatomlist(3)
     l = iatomlist(4)
     m = iatomlist(5)
     nn_iatomlist=[i, j, k, l, m, 3, 4, 5, 7, 13]

     !--- loop over atoms
     do ii = 1, 10
        iii = nn_iatomlist(ii)
        !--- loop over x,y,z
        do jj = 1, 3
           xtmp = x
           xtmp(jj, iii) = x(jj, iii)+step
           coord1 = 0.5*(abs(evalSP1(xtmp, i, j, k, l, m)) &
                        +abs(evalSP1(xtmp, 3, 4, 5, 7, 13)))
           xtmp(jj, iii) = x(jj, iii)-step
           coord2 = 0.5*(abs(evalSP1(xtmp, i, j, k, l, m)) &
                        +abs(evalSP1(xtmp, 3, 4, 5, 7, 13)))
           bmatrow(3*iii-3+jj) = (coord1-coord2)/step/2
        end do
     end do
  !--- special coordinate replacing torsion CCSH
  !--- numerical gradient
  case(7)
     step = 1d-5
     i = iatomlist(1)
     j = iatomlist(2)
     k = iatomlist(3)
     l = iatomlist(4)
     m = iatomlist(5)
     !--- loop over atoms
     do ii = 1, 5
        iii = iatomlist(ii)
        !--- loop over x,y,z
        do jj = 1, 3
           xtmp = x
           xtmp(jj, iii) = x(jj, iii)+step
           coord1 = evalSP1(xtmp, i, j, k, l, m)
           xtmp(jj, iii) = x(jj, iii)-step
           coord2 = evalSP1(xtmp, i, j, k, l, m)
           bmatrow(3*iii-3+jj) = (coord1-coord2)/step/2
        end do
     end do
  end select
  end subroutine calcbmat

  end subroutine pot
