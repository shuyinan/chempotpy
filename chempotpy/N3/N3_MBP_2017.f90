!**********************************************************************
!   System:                     N3
!   Functional form:
!   Common name:		N3
!   Number of derivatives:      0
!   Number of bodies:           3
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References: doi 10.1063/1.4983813
!               Author:  Tapan K. Mankodi, Upendra V. Bhandarkar,
!               Bhalchandra P. Puranik
!
!   Notes:
!
!***********************************************************************

!***********************************************************************!
! Global Potential energy surface for N3 system using CASSCF+CASPT2  	!
! Permutational Invariant Global Least Square fit                       !
! Input Variable: Interatomic distance R1-R2-R3 in Angstrom             !
! Output Variable: Energy in kcal/mol                                   !	
! Author of the Code: Tapan K. Mankodi                                  !
! 		      Indian Institute of Technology Bombay             !
! 		      Mumbai, India                                     !
! Ref: http://aip.scitation.org/doi/abs/10.1063/1.4983813               !
! Author:  Tapan K. Mankodi, Upendra V. Bhandarkar,                     !
!    	   Bhalchandra P. Puranik                                       ! 
! Corresponding Address: tapan.mankodi@iitb.ac.in                       !
!						 bhandarkar@iitb.ac.in  !					
!***********************************************************************!


subroutine pes(x,igrad,p,g,d)

      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: v, tx(3), sum_d
      integer :: ibond, iatom, idir, j, istate
      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      if (igrad==0) then
        tx(1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2+(x(1,3)-x(2,3))**2)
        tx(2)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2+(x(1,3)-x(3,3))**2)
        tx(3)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2+(x(2,3)-x(3,3))**2)
        call n3(v,tx(1),tx(2),tx(3))
        v=v/23.0609
        do istate=1,nstates
          p(istate)=v
        enddo
      else
        write (*,*) 'Only energy is available'
      endif

      endsubroutine


subroutine n3(v,rx,ry,rz)

	implicit none
        double precision, intent(in) :: rx, ry, rz
        double precision, intent(out) :: v

	integer							::	mpoints
	double precision				::	coeff(52)
	double precision				::	indices(52,3)
	double precision				::	ae,re
        double precision :: xrx,xry,xrz

	integer							::	i,j,k
	double precision				::	a
        integer                 ::      ii

	re = 1.098d0
	ae = 0.9d0


        xrx = exp(-(rx-re)/ae)
        xry = exp(-(ry-re)/ae)
        xrz = exp(-(rz-re)/ae)

        call assign(coeff, indices)
        
        v=0.d0

        mpoints=52

        do ii=1,mpoints
          v = v + coeff(ii)*(xrx**(indices(ii,1))*xry**(indices(ii,2))*xrz**(indices(ii,3)) + &
          xrx**(indices(ii,1))*xry**(indices(ii,3))*xrz**(indices(ii,2)) + &
          xrx**(indices(ii,2))*xry**(indices(ii,3))*xrz**(indices(ii,1)) + &
          xrx**(indices(ii,2))*xry**(indices(ii,1))*xrz**(indices(ii,3)) + &
          xrx**(indices(ii,3))*xry**(indices(ii,1))*xrz**(indices(ii,2)) + &
          xrx**(indices(ii,3))*xry**(indices(ii,2))*xrz**(indices(ii,1)))
        enddo
        
endsubroutine

subroutine assign(coeff, indices)
  double precision,intent(inout) :: coeff(52), indices(52,3)

indices(1,1)=0.0
indices(1,2)=0.0
indices(1,3)=0.0
indices(2,1)=0.0
indices(2,2)=0.0
indices(2,3)=1.0
indices(3,1)=0.0
indices(3,2)=0.0
indices(3,3)=2.0
indices(4,1)=0.0
indices(4,2)=0.0
indices(4,3)=3.0
indices(5,1)=0.0
indices(5,2)=0.0
indices(5,3)=4.0
indices(6,1)=0.0
indices(6,2)=0.0
indices(6,3)=5.0
indices(7,1)=0.0
indices(7,2)=0.0
indices(7,3)=6.0
indices(8,1)=0.0
indices(8,2)=0.0
indices(8,3)=7.0
indices(9,1)=0.0
indices(9,2)=0.0
indices(9,3)=8.0
indices(10,1)=0.0
indices(10,2)=1.0
indices(10,3)=1.0
indices(11,1)=0.0
indices(11,2)=1.0
indices(11,3)=2.0
indices(12,1)=0.0
indices(12,2)=1.0
indices(12,3)=3.0
indices(13,1)=0.0
indices(13,2)=1.0
indices(13,3)=4.0
indices(14,1)=0.0
indices(14,2)=1.0
indices(14,3)=5.0
indices(15,1)=0.0
indices(15,2)=1.0
indices(15,3)=6.0
indices(16,1)=0.0
indices(16,2)=1.0
indices(16,3)=7.0
indices(17,1)=0.0
indices(17,2)=1.0
indices(17,3)=8.0
indices(18,1)=0.0
indices(18,2)=2.0
indices(18,3)=2.0
indices(19,1)=0.0
indices(19,2)=2.0
indices(19,3)=3.0
indices(20,1)=0.0
indices(20,2)=2.0
indices(20,3)=4.0
indices(21,1)=0.0
indices(21,2)=2.0
indices(21,3)=5.0
indices(22,1)=0.0
indices(22,2)=2.0
indices(22,3)=6.0
indices(23,1)=0.0
indices(23,2)=2.0
indices(23,3)=7.0
indices(24,1)=0.0
indices(24,2)=3.0
indices(24,3)=3.0
indices(25,1)=0.0
indices(25,2)=3.0
indices(25,3)=4.0
indices(26,1)=0.0
indices(26,2)=3.0
indices(26,3)=5.0
indices(27,1)=0.0
indices(27,2)=3.0
indices(27,3)=6.0
indices(28,1)=0.0
indices(28,2)=4.0
indices(28,3)=4.0
indices(29,1)=0.0
indices(29,2)=4.0
indices(29,3)=5.0
indices(30,1)=1.0
indices(30,2)=1.0
indices(30,3)=1.0
indices(31,1)=1.0
indices(31,2)=1.0
indices(31,3)=2.0
indices(32,1)=1.0
indices(32,2)=1.0
indices(32,3)=3.0
indices(33,1)=1.0
indices(33,2)=1.0
indices(33,3)=4.0
indices(34,1)=1.0
indices(34,2)=1.0
indices(34,3)=5.0
indices(35,1)=1.0
indices(35,2)=1.0
indices(35,3)=6.0
indices(36,1)=1.0
indices(36,2)=1.0
indices(36,3)=7.0
indices(37,1)=1.0
indices(37,2)=2.0
indices(37,3)=2.0
indices(38,1)=1.0
indices(38,2)=2.0
indices(38,3)=3.0
indices(39,1)=1.0
indices(39,2)=2.0
indices(39,3)=4.0
indices(40,1)=1.0
indices(40,2)=2.0
indices(40,3)=5.0
indices(41,1)=1.0
indices(41,2)=2.0
indices(41,3)=6.0
indices(42,1)=1.0
indices(42,2)=3.0
indices(42,3)=3.0
indices(43,1)=1.0
indices(43,2)=3.0
indices(43,3)=4.0
indices(44,1)=1.0
indices(44,2)=3.0
indices(44,3)=5.0
indices(45,1)=1.0
indices(45,2)=4.0
indices(45,3)=4.0
indices(46,1)=2.0
indices(46,2)=2.0
indices(46,3)=2.0
indices(47,1)=2.0
indices(47,2)=2.0
indices(47,3)=3.0
indices(48,1)=2.0
indices(48,2)=2.0
indices(48,3)=4.0
indices(49,1)=2.0
indices(49,2)=2.0
indices(49,3)=5.0
indices(50,1)=2.0
indices(50,2)=3.0
indices(50,3)=3.0
indices(51,1)=2.0
indices(51,2)=3.0
indices(51,3)=4.0
indices(52,1)=3.0
indices(52,2)=3.0
indices(52,3)=3.0

coeff(1)=76.2027492386
coeff(2)=-29.2187456109
coeff(3)=478.216454362
coeff(4)=-2147.2354857
coeff(5)=2855.65549848
coeff(6)=-1866.34740294
coeff(7)=749.316219703
coeff(8)=-185.211752171
coeff(9)=29.2339861286
coeff(10)=481.859972236
coeff(11)=-3544.42763148
coeff(12)=3643.93394468
coeff(13)=2327.69020666
coeff(14)=-7937.75105597
coeff(15)=6936.84223471
coeff(16)=-2809.62542808
coeff(17)=460.279175648
coeff(18)=3281.89897882
coeff(19)=-1449.68238423
coeff(20)=-855.567473508
coeff(21)=-12359.5497875
coeff(22)=14201.5436577
coeff(23)=-3935.55084811
coeff(24)=3345.76391264
coeff(25)=5781.39695959
coeff(26)=-12081.0026827
coeff(27)=3039.16420473
coeff(28)=1079.60292321
coeff(29)=1320.20183514
coeff(30)=1005.42318879
coeff(31)=-2637.2433705
coeff(32)=-13723.319821
coeff(33)=14091.6436177
coeff(34)=14070.1822273
coeff(35)=-19986.7084944
coeff(36)=5578.36737242
coeff(37)=17011.6685761
coeff(38)=-49075.98983
coeff(39)=9702.54716608
coeff(40)=11609.1404971
coeff(41)=-3108.62534529
coeff(42)=16494.8907779
coeff(43)=1436.38756694
coeff(44)=-6583.14188182
coeff(45)=-2255.72399868
coeff(46)=18893.0012481
coeff(47)=-13735.259834
coeff(48)=-18895.4441394
coeff(49)=5348.42098676
coeff(50)=10801.560473
coeff(51)=10627.4574663
coeff(52)=-7482.93770079

endsubroutine 
