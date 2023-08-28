!********************************************************************** 
!   System:                     NO2
!   Name of this surface:       NO2_6Ap_MB-PIP-MEG
!   Functional form:            permutationally invariant polynomials
!   Common name:                NO2(adiabatic sextet ground state)
!   Number of derivatives:      1
!   Number of bodies:           3
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References:: Z. Varga, Y. Liu, J. Li, Y. Paukku, H. Guo, and
!                D. G. Truhlar, in preparation
!
!   Notes:    -Sextet A' (6Ap) surface 
!             -Mixed-exponential-Gaussian (MEG)
!              variables are applied
!             -12-th order polynomial
!
!     O1--O2
!
!     N3
!
!   Input: X(3),Y(3),Z(3)               in units of bohr
!   Output: E                           in units of hartree
!   Output: dEdX(3),dEdY(3),dEdZ(3)     in units of hartree/bohr
!**********************************************************************
      module NO2_6Ap_ZV_par

! Conversion factors
! Cconv: bohr to Angstrom 
!        1 bohr = 0.52917721092 angstrom
! Econv: kcal/mol to hartree 
!        1 kcal/mol = 0.159360144 * 10^-2 hartree
! Gconv: kcal/(mol*Angstrom) to hartree/bohr
!        1 kcal mol^-1 angstrom^-1 = 0.843297564 * 10^-3 hartree/bohr
      double precision,parameter :: Cconv = 0.52917721092d0
      double precision,parameter :: Econv = 0.159360144d-2
      double precision,parameter :: Gconv = 0.843297564d-3
! Common variables
! R(3):          Interatomic bond distance
! rMs(3)             Array to store base terms
! rM(0:4)            Array to store monomials
! P(0:251)           Array to store polynomials
! B(1:227)       Array to store basis functions
! dMsdR(3,3):    The derivative of base terms w.r.t. R
! dMdR(3,0:3):   The derivative of monomials w.r.t. R
! dPdR(3,0:251): The derivative of basis functions w.r.t. R
! dVdR(3):       The derivative of V w.r.t. R
! dRdX(3,9):     The derivative of R w.r.t. X
! dBdR(3,227)    The derivative of B w.r.t. R 
  
      double precision :: R(3)
      double precision :: rMs(3),rM(0:4),P(0:251),B(227)
      double precision :: dMsdR(3,3),dMdR(3,0:4),dPdR(3,0:251)
      double precision :: dVdR(3),dRdX(3,9),dBdR(3,227)
  
! Nonlinear parameters:
! a(in Ang)
! ab (in Ang^2)
! re (in Ang)
! order of data: OO,NO,NO,NO,NO,NN
      double precision,parameter :: a(3)  = (/ 1.25d0,0.94d0,0.94d0 /)
      double precision,parameter :: ab(3) = (/ 2.10d0,3.90d0,3.90d0 /)
      double precision,parameter :: re(3) = (/ 1.208d0,1.4219d0,&
      1.4219d0 /) 
! Reference energy of infinitely separated four atoms in hartree
      double precision,parameter :: Eref = -204.78312431d0
! Initialized v to be De(N2) + De(O2) = 228.4 + 120.143 kcal/mol
      double precision,parameter :: totdiss = 348.643d0
! Linear parameters optimized by the weighted-least square fitting
      double precision,parameter :: C(227) = (/ &
         0.197851973467D+03 &
      , -0.255484057397D+03 &
      ,  0.102451508727D+04 &
      , -0.126574980991D+04 &
      ,  0.275229975480D+04 &
      ,  0.733673569320D+03 &
      , -0.518870736871D+04 &
      ,  0.201817758606D+04 &
      ,  0.891703858476D+03 &
      ,  0.507316828077D+04 &
      , -0.101579474125D+05 &
      , -0.350104191179D+04 &
      , -0.337155016986D+04 &
      ,  0.971442829339D+04 &
      ,  0.168379520427D+05 &
      , -0.105186569326D+05 &
      ,  0.123602811150D+04 &
      , -0.259860701557D+04 &
      , -0.583209633828D+04 &
      ,  0.152118875687D+05 &
      ,  0.991348255848D+04 &
      ,  0.600977702834D+04 &
      ,  0.711983606150D+04 &
      , -0.112666103912D+05 &
      , -0.196576443468D+05 &
      , -0.152049051200D+05 &
      ,  0.707407150137D+04 &
      ,  0.195350176198D+05 &
      , -0.445877017302D+04 &
      , -0.109170577599D+05 &
      ,  0.136475419992D+05 &
      , -0.713740263918D+04 &
      , -0.602423893153D+04 &
      , -0.155469012866D+05 &
      ,  0.842567027256D+04 &
      , -0.116600624284D+05 &
      , -0.568317643645D+04 &
      , -0.573210695082D+04 &
      ,  0.162954212098D+05 &
      , -0.330644370138D+04 &
      , -0.679270263122D+04 &
      ,  0.286451564814D+05 &
      ,  0.145220800878D+05 &
      , -0.400234667228D+05 &
      ,  0.851746181247D+04 &
      , -0.625017097990D+04 &
      ,  0.459881481068D+04 &
      , -0.118928671299D+05 &
      ,  0.154696846468D+05 &
      , -0.759688748088D+04 &
      ,  0.552526417473D+04 &
      , -0.121643944936D+05 &
      , -0.726093907242D+04 &
      ,  0.445906853203D+03 &
      , -0.142216891785D+04 &
      , -0.598237155146D+04 &
      ,  0.155914829358D+05 &
      ,  0.622685213824D+04 &
      ,  0.769234646079D+04 &
      ,  0.741994254446D+04 &
      , -0.652683939907D+04 &
      ,  0.117787810809D+05 &
      , -0.739058011120D+04 &
      ,  0.185247130315D+04 &
      , -0.265371525365D+04 &
      , -0.982702489649D+04 &
      ,  0.261478366815D+04 &
      ,  0.316387605270D+04 &
      , -0.138684003793D+04 &
      , -0.417191707457D+04 &
      ,  0.126919676090D+04 &
      ,  0.706451507056D+04 &
      ,  0.322749421172D+04 &
      ,  0.105085457929D+05 &
      ,  0.760100817180D+04 &
      , -0.154452941021D+03 &
      ,  0.116591690913D+05 &
      ,  0.545158062256D+04 &
      , -0.147009768999D+05 &
      ,  0.837232861928D+04 &
      , -0.215839979930D+05 &
      ,  0.259296687114D+04 &
      , -0.200350447569D+05 &
      ,  0.715220756699D+04 &
      , -0.287718965267D+05 &
      ,  0.465283734010D+05 &
      , -0.293538119662D+05 &
      ,  0.360801302323D+04 &
      , -0.711375038801D+04 &
      ,  0.207165303287D+05 &
      , -0.429871435991D+04 &
      ,  0.146086049742D+05 &
      , -0.716645835750D+04 &
      ,  0.621995624878D+04 &
      , -0.473847386335D+04 &
      ,  0.105333303930D+05 &
      ,  0.263169711349D+04 &
      , -0.172818407904D+05 &
      ,  0.463716566887D+04 &
      , -0.281408500373D+04 &
      , -0.894014094886D+04 &
      ,  0.200179311844D+05 &
      , -0.441979470138D+03 &
      ,  0.323489831002D+04 &
      , -0.907414122311D+03 &
      ,  0.508257104647D+02 &
      ,  0.936674512833D+04 &
      ,  0.129397387217D+05 &
      , -0.208555663847D+05 &
      ,  0.641186135074D+04 &
      , -0.116098450869D+05 &
      , -0.997845013200D+04 &
      ,  0.382079016880D+05 &
      ,  0.100290883426D+05 &
      ,  0.248234518634D+05 &
      , -0.245829719823D+05 &
      ,  0.308607637422D+05 &
      ,  0.234677776723D+05 &
      , -0.375273450900D+05 &
      , -0.109881222021D+05 &
      , -0.571478320286D+04 &
      ,  0.814295997139D+04 &
      ,  0.284182465748D+04 &
      ,  0.104086164722D+05 &
      , -0.410536834285D+03 &
      ,  0.354029040396D+03 &
      , -0.843965401937D+04 &
      ,  0.906493399432D+04 &
      , -0.172590472400D+05 &
      ,  0.952079534524D+04 &
      ,  0.126554024450D+05 &
      , -0.981471615810D+04 &
      ,  0.479612615911D+04 &
      , -0.516333206145D+04 &
      , -0.168329719551D+05 &
      , -0.848459247106D+04 &
      , -0.126788858769D+05 &
      , -0.145757850273D+04 &
      , -0.126305473974D+05 &
      , -0.610797260940D+04 &
      ,  0.844685790027D+04 &
      , -0.142421471160D+05 &
      , -0.715744509654D+04 &
      ,  0.246162603687D+05 &
      , -0.824628990267D+04 &
      ,  0.117111101687D+05 &
      ,  0.351150401387D+03 &
      , -0.292494493414D+05 &
      ,  0.140305020196D+04 &
      , -0.422729151490D+05 &
      ,  0.216155783998D+05 &
      ,  0.818521936322D+04 &
      , -0.132237293517D+05 &
      , -0.258008632244D+05 &
      ,  0.280840425206D+05 &
      , -0.111869656942D+05 &
      ,  0.186793437556D+05 &
      ,  0.111486007952D+05 &
      , -0.145059853978D+04 &
      ,  0.892769732663D+03 &
      , -0.540868185268D+04 &
      ,  0.109063146952D+04 &
      , -0.146293827803D+05 &
      ,  0.423685045475D+04 &
      , -0.184042887679D+04 &
      ,  0.378425075516D+04 &
      , -0.629363370015D+04 &
      ,  0.101283066202D+05 &
      , -0.815033042720D+04 &
      , -0.342032862346D+04 &
      ,  0.528707825716D+04 &
      , -0.471057640315D+04 &
      ,  0.129718064645D+05 &
      , -0.784552516827D+04 &
      ,  0.152080863164D+05 &
      ,  0.123734944846D+05 &
      , -0.909484795350D+04 &
      ,  0.899359101011D+04 &
      ,  0.837285273537D+04 &
      ,  0.450266986362D+04 &
      , -0.408792307120D+03 &
      ,  0.307028846770D+04 &
      ,  0.104772975173D+04 &
      , -0.591207739289D+04 &
      ,  0.207157128033D+04 &
      ,  0.926540004974D+04 &
      , -0.109545683769D+05 &
      ,  0.535526969523D+04 &
      ,  0.355079207629D+04 &
      ,  0.342973320658D+03 &
      ,  0.176505414329D+04 &
      , -0.349594085218D+04 &
      ,  0.250628788190D+04 &
      , -0.842520291547D+03 &
      ,  0.138136890343D+05 &
      , -0.276137469472D+03 &
      , -0.119356198512D+05 &
      ,  0.651433421422D+04 &
      , -0.452855594330D+04 &
      ,  0.453039184733D+04 &
      ,  0.290856136590D+04 &
      , -0.136535894560D+04 &
      ,  0.159837882572D+04 &
      , -0.926400783938D+04 &
      , -0.434731735828D+04 &
      ,  0.111788331706D+05 &
      ,  0.395860487370D+04 &
      , -0.267210631811D+04 &
      ,  0.543433196451D+03 &
      , -0.641876015513D+03 &
      ,  0.415457036349D+04 &
      , -0.427222709067D+04 &
      ,  0.307291280121D+04 &
      , -0.245813507346D+04 &
      ,  0.187946788742D+04 &
      ,  0.278640613016D+03 &
      , -0.967554521946D+03 &
      ,  0.142923552437D+04 &
      , -0.584461294728D+04 &
      ,  0.119942280466D+05 &
      , -0.243417599475D+05 &
      ,  0.292835698311D+05 &
      , -0.304103983818D+05 &
      ,  0.189934641723D+05 &
      , -0.749215811846D+04 &
      , -0.123338225832D+04 &
      , -0.978741526281D+03 &
      /)

      end module NO2_6Ap_ZV_par

      subroutine pes(x,igrad,potential,gradient,dvec)

      use NO2_6Ap_ZV_par
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: potential(nstates)
      double precision, intent(out) :: gradient(nstates,natoms,3)
      double precision, intent(out) :: dvec(nstates,nstates,natoms,3)

      double precision :: v, tx(9), tg(9)
      integer :: iatom, idir, j, istate
      !initialize 
      potential=0.d0
      gradient=0.d0
      dvec=0.d0

      do iatom=1,natoms
        do idir=1,3
          j=3*(iatom-1)+idir
          tx(j)=x(iatom, idir)
        enddo
      enddo
      tg=0.d0

      call no2pes(tx,v,tg,igrad)

      v=v/23.0609 !+Eref*27.211386d0
      tg=tg/23.0609

      do istate=1,nstates
        potential(istate)=v
      enddo

      do istate=1,nstates
        do iatom=1,natoms
          do idir=1,3
            j=3*(iatom-1)+idir
            gradient(istate,iatom,idir)=tg(j)
          enddo
        enddo
      enddo

      dvec=0.d0

      endsubroutine

!********************************************************************** 
! Local variables used in the N2O2 PES subroutine
! input coordinate matrix: X(9)in Ang
!                          O1: X(1),X(2),X(3)
!                          O2: X(4),X(5),X(6)
!                          N3: X(7),X(8),X(9)
! input flag: igrad     igrad=0 energy calculation only
!                       igrad=1 calculation of energy and gradient
! output potential energy:      v    in kcal/mol
! output gradient:              dVdX in kcal/(mol*Ang)
!********************************************************************** 
      subroutine no2pes(X,v,dVdX,igrad)
      use NO2_6Ap_ZV_par
!**********************************************************************
! Subroutine to calculate the potential energy V and gradient dVdX
! for given Cartesian coordinates X(9)  
! R:            Interatomic bond distances (3)
! V:            Calculated potential energy
! dVdX:         The derivative of V w.r.t. X, dim(9)
! dVdR:         The derivative of V w.r.t. R, dim(3) 
! dPdR:         The derivative of basis functions w.r.t. R
!               dim(3*161)
! dMdR:         The derivative of monomials w.r.t. R
!               dim(3*3)
! dRdX:         The derivative of R w.r.t. X, dim(3*9)
!**********************************************************************
      integer i,igrad,j,nob,k
      double precision V
      double precision dVdX(9),X(9) 
! Obtain Cartesian coordinates from input file
      call coord_convt(X)
      if (igrad .le. 1) then
! Call subroutine Evv to evaluate potential energy V
        call evv(V)
        if (igrad .eq. 1) then
! Call subroutine EvdVdX to evaluate the derivatives of V w.r.t. X
          call evdvdx(X,dVdX)
        endif
      else
        write (*,*) 'Only igrad = 0, 1 is allowed!'
      endif
      end subroutine no2pes
      subroutine coord_convt(X)
      use NO2_6Ap_ZV_par
!**********************************************************************
!  This subroutine calculates the three interatomic distances 
!  from the XYZ coordinates
!
!  r1 = r(O1O2)
!  r2 = r(O1N3)
!  r3 = r(O2N3)
!
!**********************************************************************
      integer i
      double precision X(9)
      R(1)=Sqrt((X(4)-X(1))**2  + (X(5)-X(2))**2  + (X(6)-X(3))**2)
      R(2)=Sqrt((X(7)-X(1))**2  + (X(8)-X(2))**2  + (X(9)-X(3))**2)
      R(3)=Sqrt((X(4)-X(7))**2  + (X(5)-X(8))**2  + (X(6)-X(9))**2)
 
      return
      end subroutine coord_convt
      subroutine EvV(V)
      use NO2_6Ap_ZV_par
!**********************************************************************
! This subroutine evaluates V for given R 
! V(R) = C*B
! C:            Coefficients 
! P:            Basis functions evaluated for given R
! rMs:          rMs(3), three mixed exponential Gaussian (MEG) terms
! a:            Nonlinear parameters in MEG terms (Angstrom)
! ab:           Nonlinear parameters in MEG terms (Angstrom^2)
! re:           Equilibrium bond length (Angstrom)
! nop:          number of points
! nom:          number of monomials
! nob:          number of basis functions (polynomials)
! rM(0:4):      Array to store monomials
! P(0:251):     Array to store polynomials
! B(1:227): Array to store basis functions
!**********************************************************************
      integer i,j,k
      double precision dist,dv2dr,V,V21,V22,V23
      double precision Vpf,dVpfdR(3)
! Calculate the six MEG terms for each point
      call evmeg
! Calculate the monomials for each point by using six MEG terms
      call evmono
! Calculate the polynomials (basis functions) by using monomials
      call evpoly 
! Calculate the basis functions by removing unconnected and 2-body terms
      call evbas
! Initialized v to be totdiss
      v=totdiss 
! Evaluate 2-body interactions
 
        dist=r(1)
        call ev2gm2(dist,v21,dv2dr,4,0)         ! OO
        dist=r(2)
        call ev2gm2(dist,v22,dv2dr,3,0)         ! NO
        dist=r(3)
        call ev2gm2(dist,v23,dv2dr,3,0)         ! NO
        v= v + v21 + v22 + v23
! Evaluate V by taken the product of C and Basis function arrays
      do i=1,227
         v=v + c(i)*b(i)
      enddo
! Add a patch function to get the barrier between 2A1 and 2B2 minima
!       call patchfunc(Vpf,dVpfdR,0)
!       v=v + Vpf
!      Write(*,9999) V 
! 9999 Format('The potential energy is ',F20.14,' kcal/mol')
      return
      end subroutine EvV
      subroutine EvdVdX(X,dVdX)
      use NO2_6Ap_ZV_par
!**********************************************************************
! This subroutine evaluates dRdX for given R and X 
! R:            R(3), 3 bond lengths
! X:            X(9), 12 Cartesian coordinates
! rM(0:4):    Array to store monomials
! P(0:251):   Array to store polynomials
! dVdX:         dVdX(9), derivatives of V w.r.t. Cartesian coordinates 
! dVdR:         dVdR(3), derivatives of V w.r.t. 3 bond lengths
! dRdX:         dRdX(3,9), derivatives of R(3) w.r.t. 9  
!               Cartesian coordinates
!**********************************************************************
      integer i,j
      double precision dVdX(9),X(9)
! Initialize dVdX
      do i=1,9
        dVdX(i)=0.0d0
      enddo
! Call EvdVdR to evaluate dVdR(3)
      Call evdvdr
! Call EvdRdX to evaluate dRdX(3,9)
      Call evdrdx(X)  
! Calculate dVdX by using chain rule: dV/dXi=(dV/dRj)*(dRj/dXi),
! j=1 to 3
      do i=1,9
        do j=1,3
          dVdX(i)=dVdX(i) + dVdR(j)*dRdX(j,i)
        enddo
      enddo
      return
      end subroutine EvdVdX
      subroutine EvMEG
      use NO2_6Ap_ZV_par
!**********************************************************************
!  This subroutine calculates the MEG function
!  MEG term rms = exp(-(r-re)/a-(r-re)^2/ab)
!  re:   equilibrium bond length (Angstrom)
!  a:    nonlinear parameter (Angstrom)
!  ab:   nonlinear parameter (Angstrom^2)
!**********************************************************************
      integer i
      do i=1,3
         rms(i)=Exp(-(r(i)-re(i))/a(i)-(r(i)-re(i))**2.0d0/ab(i))
      enddo
      end subroutine EvMEG
      subroutine EvMono
      use NO2_6Ap_ZV_par
!**********************************************************************
!  This subroutine takes six MEG variables(X) and calculates the
!  monomials(RM) and that do not have usable decomposition.
!  For degree 10, the number of monomials is nom.
!**********************************************************************

      integer i

      rm(0) = 1.0d0
      rm(1) = rms(3)
      rm(2) = rms(2)
      rm(3) = rms(1)
      rm(4) = rm(1)*rm(2)

      return
      end subroutine EvMono
      subroutine EvPoly
      use NO2_6Ap_ZV_par
!**********************************************************************
!  This subroutine takes monomials(RM) and calculates the
!  permutationally invariant polynomials(p)
!  For degree 10, the number of polynomials is nob.
!**********************************************************************
      p(0) = rm(0)
      p(1) = rm(1) + rm(2)
      p(2) = rm(3)
      p(3) = rm(4)
      p(4) = p(2)*p(1)
      p(5) = p(1)*p(1) - p(3) - p(3)
      p(6) = p(2)*p(2)
      p(7) = p(2)*p(3)
      p(8) = p(1)*p(3)
      p(9) = p(2)*p(5)
      p(10) = p(2)*p(4)
      p(11) = p(1)*p(5) - p(8)
      p(12) = p(2)*p(6)
      p(13) = p(2)*p(8)
      p(14) = p(2)*p(7)
      p(15) = p(3)*p(3)
      p(16) = p(3)*p(5)
      p(17) = p(2)*p(11)
      p(18) = p(2)*p(9)
      p(19) = p(2)*p(10)
      p(20) = p(1)*p(11) - p(16)
      p(21) = p(2)*p(12)
      p(22) = p(2)*p(15)
      p(23) = p(2)*p(16)
      p(24) = p(2)*p(13)
      p(25) = p(2)*p(14)
      p(26) = p(1)*p(15)
      p(27) = p(3)*p(11)
      p(28) = p(2)*p(20)
      p(29) = p(2)*p(17)
      p(30) = p(2)*p(18)
      p(31) = p(2)*p(19)
      p(32) = p(1)*p(20) - p(27)
      p(33) = p(2)*p(21)
      p(34) = p(2)*p(26)
      p(35) = p(2)*p(27)
      p(36) = p(2)*p(22)
      p(37) = p(2)*p(23)
      p(38) = p(2)*p(24)
      p(39) = p(2)*p(25)
      p(40) = p(3)*p(15)
      p(41) = p(3)*p(16)
      p(42) = p(3)*p(20)
      p(43) = p(2)*p(32)
      p(44) = p(2)*p(28)
      p(45) = p(2)*p(29)
      p(46) = p(2)*p(30)
      p(47) = p(2)*p(31)
      p(48) = p(1)*p(32) - p(42)
      p(49) = p(2)*p(33)
      p(50) = p(2)*p(40)
      p(51) = p(2)*p(41)
      p(52) = p(2)*p(42)
      p(53) = p(2)*p(34)
      p(54) = p(2)*p(35)
      p(55) = p(2)*p(36)
      p(56) = p(2)*p(37)
      p(57) = p(2)*p(38)
      p(58) = p(2)*p(39)
      p(59) = p(1)*p(40)
      p(60) = p(3)*p(27)
      p(61) = p(3)*p(32)
      p(62) = p(2)*p(48)
      p(63) = p(2)*p(43)
      p(64) = p(2)*p(44)
      p(65) = p(2)*p(45)
      p(66) = p(2)*p(46)
      p(67) = p(2)*p(47)
      p(68) = p(1)*p(48) - p(61)
      p(69) = p(2)*p(49)
      p(70) = p(2)*p(59)
      p(71) = p(2)*p(60)
      p(72) = p(2)*p(61)
      p(73) = p(2)*p(50)
      p(74) = p(2)*p(51)
      p(75) = p(2)*p(52)
      p(76) = p(2)*p(53)
      p(77) = p(2)*p(54)
      p(78) = p(2)*p(55)
      p(79) = p(2)*p(56)
      p(80) = p(2)*p(57)
      p(81) = p(2)*p(58)
      p(82) = p(3)*p(40)
      p(83) = p(3)*p(41)
      p(84) = p(3)*p(42)
      p(85) = p(3)*p(48)
      p(86) = p(2)*p(68)
      p(87) = p(2)*p(62)
      p(88) = p(2)*p(63)
      p(89) = p(2)*p(64)
      p(90) = p(2)*p(65)
      p(91) = p(2)*p(66)
      p(92) = p(2)*p(67)
      p(93) = p(1)*p(68) - p(85)
      p(94) = p(2)*p(69)
      p(95) = p(2)*p(82)
      p(96) = p(2)*p(83)
      p(97) = p(2)*p(84)
      p(98) = p(2)*p(85)
      p(99) = p(2)*p(70)
      p(100) = p(2)*p(71)
      p(101) = p(2)*p(72)
      p(102) = p(2)*p(73)
      p(103) = p(2)*p(74)
      p(104) = p(2)*p(75)
      p(105) = p(2)*p(76)
      p(106) = p(2)*p(77)
      p(107) = p(2)*p(78)
      p(108) = p(2)*p(79)
      p(109) = p(2)*p(80)
      p(110) = p(2)*p(81)
      p(111) = p(1)*p(82)
      p(112) = p(3)*p(60)
      p(113) = p(3)*p(61)
      p(114) = p(3)*p(68)
      p(115) = p(2)*p(93)
      p(116) = p(2)*p(86)
      p(117) = p(2)*p(87)
      p(118) = p(2)*p(88)
      p(119) = p(2)*p(89)
      p(120) = p(2)*p(90)
      p(121) = p(2)*p(91)
      p(122) = p(2)*p(92)
      p(123) = p(1)*p(93) - p(114)
      p(124) = p(2)*p(94)
      p(125) = p(2)*p(111)
      p(126) = p(2)*p(112)
      p(127) = p(2)*p(113)
      p(128) = p(2)*p(114)
      p(129) = p(2)*p(95)
      p(130) = p(2)*p(96)
      p(131) = p(2)*p(97)
      p(132) = p(2)*p(98)
      p(133) = p(2)*p(99)
      p(134) = p(2)*p(100)
      p(135) = p(2)*p(101)
      p(136) = p(2)*p(102)
      p(137) = p(2)*p(103)
      p(138) = p(2)*p(104)
      p(139) = p(2)*p(105)
      p(140) = p(2)*p(106)
      p(141) = p(2)*p(107)
      p(142) = p(2)*p(108)
      p(143) = p(2)*p(109)
      p(144) = p(2)*p(110)
      p(145) = p(3)*p(82)
      p(146) = p(3)*p(83)
      p(147) = p(3)*p(84)
      p(148) = p(3)*p(85)
      p(149) = p(3)*p(93)
      p(150) = p(2)*p(123)
      p(151) = p(2)*p(115)
      p(152) = p(2)*p(116)
      p(153) = p(2)*p(117)
      p(154) = p(2)*p(118)
      p(155) = p(2)*p(119)
      p(156) = p(2)*p(120)
      p(157) = p(2)*p(121)
      p(158) = p(2)*p(122)
      p(159) = p(1)*p(123) - p(149)
      p(160) = p(2)*p(124)
      p(161) = p(2)*p(145)
      p(162) = p(2)*p(146)
      p(163) = p(2)*p(147)
      p(164) = p(2)*p(148)
      p(165) = p(2)*p(149)
      p(166) = p(2)*p(125)
      p(167) = p(2)*p(126)
      p(168) = p(2)*p(127)
      p(169) = p(2)*p(128)
      p(170) = p(2)*p(129)
      p(171) = p(2)*p(130)
      p(172) = p(2)*p(131)
      p(173) = p(2)*p(132)
      p(174) = p(2)*p(133)
      p(175) = p(2)*p(134)
      p(176) = p(2)*p(135)
      p(177) = p(2)*p(136)
      p(178) = p(2)*p(137)
      p(179) = p(2)*p(138)
      p(180) = p(2)*p(139)
      p(181) = p(2)*p(140)
      p(182) = p(2)*p(141)
      p(183) = p(2)*p(142)
      p(184) = p(2)*p(143)
      p(185) = p(2)*p(144)
      p(186) = p(1)*p(145)
      p(187) = p(3)*p(112)
      p(188) = p(3)*p(113)
      p(189) = p(3)*p(114)
      p(190) = p(3)*p(123)
      p(191) = p(2)*p(159)
      p(192) = p(2)*p(150)
      p(193) = p(2)*p(151)
      p(194) = p(2)*p(152)
      p(195) = p(2)*p(153)
      p(196) = p(2)*p(154)
      p(197) = p(2)*p(155)
      p(198) = p(2)*p(156)
      p(199) = p(2)*p(157)
      p(200) = p(2)*p(158)
      p(201) = p(1)*p(159) - p(190)
      p(202) = p(2)*p(160)
      p(203) = p(2)*p(186)
      p(204) = p(2)*p(187)
      p(205) = p(2)*p(188)
      p(206) = p(2)*p(189)
      p(207) = p(2)*p(190)
      p(208) = p(2)*p(161)
      p(209) = p(2)*p(162)
      p(210) = p(2)*p(163)
      p(211) = p(2)*p(164)
      p(212) = p(2)*p(165)
      p(213) = p(2)*p(166)
      p(214) = p(2)*p(167)
      p(215) = p(2)*p(168)
      p(216) = p(2)*p(169)
      p(217) = p(2)*p(170)
      p(218) = p(2)*p(171)
      p(219) = p(2)*p(172)
      p(220) = p(2)*p(173)
      p(221) = p(2)*p(174)
      p(222) = p(2)*p(175)
      p(223) = p(2)*p(176)
      p(224) = p(2)*p(177)
      p(225) = p(2)*p(178)
      p(226) = p(2)*p(179)
      p(227) = p(2)*p(180)
      p(228) = p(2)*p(181)
      p(229) = p(2)*p(182)
      p(230) = p(2)*p(183)
      p(231) = p(2)*p(184)
      p(232) = p(2)*p(185)
      p(233) = p(3)*p(145)
      p(234) = p(3)*p(146)
      p(235) = p(3)*p(147)
      p(236) = p(3)*p(148)
      p(237) = p(3)*p(149)
      p(238) = p(3)*p(159)
      p(239) = p(2)*p(201)
      p(240) = p(2)*p(191)
      p(241) = p(2)*p(192)
      p(242) = p(2)*p(193)
      p(243) = p(2)*p(194)
      p(244) = p(2)*p(195)
      p(245) = p(2)*p(196)
      p(246) = p(2)*p(197)
      p(247) = p(2)*p(198)
      p(248) = p(2)*p(199)
      p(249) = p(2)*p(200)
      p(250) = p(1)*p(201) - p(238)
      p(251) = p(2)*p(202)

      return
      end subroutine EvPoly
      subroutine evbas
      use NO2_6Ap_ZV_par
!**********************************************************************
!  This subroutine eliminates the 2-body terms in the permutationally
!  invariant formulation
!**********************************************************************
      
      integer i
      double precision b1(252) 
! Pass P(0:251) to BM1(1:252)
      do i=1,252
        b1(i)=p(i-1)
      enddo
! Remove unconnected terms and 2-body terms and pass to B(1:227)
      do i=1,2
        b(i)=b1(i+3)
      enddo
        do i=3,6
        b(i)=b1(i+5)
      enddo
      do i=7,13  
        b(i)=b1(i+7)
      enddo
      do i=14,23
        b(i)=b1(i+9)
      enddo
      do i=24,37
        b(i)=b1(i+11)  
      enddo
      do i=38,55
        b(i)=b1(i+13)
      enddo
      do i=56,78
        b(i)=b1(i+15)
      enddo
      do i=79,106
        b(i)=b1(i+17)     
      enddo
      do i=107,140
        b(i)=b1(i+19)
      enddo
      do i=141,180
        b(i)=b1(i+21)
      enddo
      do i=181,227
        b(i)=b1(i+23)
      enddo

      return
      end subroutine evbas
      subroutine ev2gm2(r,v,grad,imol,igrad) 
!**********************************************************************
!
! Compute the diatomic potential of ground-state triplet O2
!
! References: J. Chem. Phys. 132, 074307 (2010)
!
! Input:  r      interatomic distance in Angstrom
! Output: V      potential in kcal/mol
!         grad   gradient (kcal/mol)/Angstrom
!
!**********************************************************************
      implicit none
      integer,intent(in) :: imol, igrad
      double precision,intent(in)  :: r
      double precision,intent(out) :: v, grad
! Parameters of analytical even-tempered Gaussian expansions for the
! ground state potential energy curve of O2 CBS+SR+SO+CV.
! Units: alpha in Angstrom^-2, beta=dimensionless, as_k in milihartree.
      double precision :: alpha,beta,as(0:7)
      integer :: k 
! NO potential variables
      double precision :: re,de,c(12),u,dfdr
! Dispersion variables
      double precision :: dist(1),disp,dispdr(1)
! O2 pairwise
      if (imol .eq. 4) then
! Original parameters
!      alpha = 0.785d0
!      beta = 1.307d0
!      as(0) = -2388.5641690d0
!      as(1) = 18086.977116d0
!      as(2) = -71760.197585d0
!      as(3) = 154738.09175d0
!      as(4) = -215074.85646d0
!      as(5) = 214799.54567d0
!      as(6) = -148395.42850d0
!      as(7) = 73310.781453d0
! Modified parameters for D3(BJ)
       alpha = 9.439784362354936d-1
       beta =  1.262242998506810d0
       as(0) = -1.488979427684798d3
       as(1) =  1.881435846488955d4
       as(2) = -1.053475425838226d5
       as(3) =  2.755135591229064d5
       as(4) = -4.277588997761775d5
       as(5) =  4.404104009614092d5
       as(6) = -2.946204062950765d5
       as(7) =  1.176861219078620d5
      v=0.d0
      do k=0,7
       v= v + as(k)*dexp(-alpha*beta**k*r**2)
      enddo
! From milihartree to kcal/mol
      v=v*627.509523475149d-3
! Add D3 dispersion correction
!        write(*,*) "voo= ",v
        dist(1)=r
        call d3disp(dist,disp,dispdr,0,imol)
        v=v+disp
!        write(*,*) "dispoo= ",disp
  
! Compute the gradient if needed
      if (igrad.eq.1) then
       grad=0.d0
       do k=0,7
        grad=grad-2.d0*as(k)*alpha*beta**k*r*dexp(-alpha*beta**k*r**2)
       enddo
! Convert from milihartree/A to (kcal/mol)/A
         grad=grad*627.509523475149d-3
! Add analytical gradient of D3 dispersion correction
        call d3disp(dist,disp,dispdr,1,imol)
        grad= grad + dispdr(1)
      endif
! NO pairwise for quartet NO
      else if (imol .eq. 3) then
!!! quartet NO only for sextet NO2 surface !!!											  
! Short range NO parameters for NO dissociation
! The long range term comes from dispersion
! The difference between the points of the 8th order MEG NO diatomic
! potential and the dispersion corrections was refitted by a 10th
! order MEG polynomial
        re=   1.4219d0
        de=   42.7d0   !152.6d0
        c(1)  =  -0.8037362024081990d+00
        c(2)  =  12.8033640295562000d+00
        c(3)  = -38.6920219767776000d+00
        c(4)  =  69.7090314448948000d+00
        c(5)  = -77.2484600687765000d+00
        c(6)  =  49.2765652793876000d+00
        c(7)  = -15.2482809053106000d+00
        c(8)  =   0.3234724591257560d+00
        c(9)  =   0.9904058013501440d+00
        c(10) =  -0.1746291297566160d+00
        c(11) =   0.5122786262498600d+00 ! delta b
        c(12) =   1.1224168097377800d+00 ! delta c
!  for NO
! MEG variable
       u=exp(-(r-re)/c(11)-(r-re)**2.0d0/c(12))
! Pair-wise potential
       v=0.d0
       v= (-de)*(c(1)*u**1.0d0 + c(2)*u**2.0d0 + c(3)*u**3.0d0 &
               + c(4)*u**4.0d0 + c(5)*u**5.0d0 + c(6)*u**6.0d0 &
               + c(7)*u**7.0d0 + c(8)*u**8.0d0 + c(9)*u**9.0d0 &
               + c(10)*u**10.0d0 )
! Add D3 dispersion correction
!        write(*,*) "vno= ",v
        dist(1)=r
        call d3disp(dist,disp,dispdr,0,imol)
        v=v+disp
!        write(*,*) "dispno= ",disp
! 1st derivate by r
      if (igrad .eq. 1) then
! multiplier part of the MEG derivate
       dfdr=(-2.0d0*(r-re)/c(12)-1.0d0/c(11))
       grad =  - de*(c(1)*1.0d0*dfdr*u**1.0d0 &
                   + c(2)*2.0d0*dfdr*u**2.0d0 &
                   + c(3)*3.0d0*dfdr*u**3.0d0 &
                   + c(4)*4.0d0*dfdr*u**4.0d0 &
                   + c(5)*5.0d0*dfdr*u**5.0d0 &
                   + c(6)*6.0d0*dfdr*u**6.0d0 &
                   + c(7)*7.0d0*dfdr*u**7.0d0 &
                   + c(8)*8.0d0*dfdr*u**8.0d0 &
                   + c(9)*9.0d0*dfdr*u**9.0d0 &
                   + c(10)*10.0d0*dfdr*u**10.0d0 )
! Add analytical gradient of D3 dispersion correction
        call d3disp(dist,disp,dispdr,1,imol)
        grad= grad + dispdr(1)
      endif
      endif
      return
      end subroutine ev2gm2

       subroutine patchfunc(Vpf,dVpfdR,igrad)
      use NO2_6Ap_ZV_par
!**********************************************************************
!      A patch function to correct the barrier between 2A1 and 2B1
!      minima of the doublet surface
!**********************************************************************
       integer igrad
       double precision Vpf
       double precision tmp,tmp2,tmp3
       double precision dVpfdR(3)
       double precision rNOe,aOOe,a0,a1,a2

       a0 = 6.5d0  ! in kcal/mol
       a1 = 0.009d0 ! in Ang**2
       a2 = 0.002d0 ! in rad**2
       rNOe = 1.25d0  ! in Ang
! aOOe corresponds to cos(107 deg) 
       aOOe = cos(1.8675023d0)  ! in rad

       tmp = -(( R(2) - rNOe )**2)/a1
       tmp = tmp -(( R(3) - rNOe )**2)/a1
       tmp2 = ( R(2)**2 + R(3)**2 - R(1)**2 )/(2*R(2)*R(3))
       tmp = tmp -(( tmp2 - aOOe ))**2/a2
       Vpf = a0 * Exp(tmp)

       if (igrad .eq. 1) then       
! R(1) is O1O2
       tmp3 = 2.0d0* R(1)*(tmp2 - aOOe)/(a2*R(2)*R(3))
       dVpfdR(1) = Vpf * tmp3
! R(2) is NO1
       tmp3 = -2.0d0*(R(2)-rNOe)/a1-2.0d0*(1/R(3)-tmp2)*(tmp2-aOOe)/a2
       dVpfdR(2) = Vpf * tmp3
! R(3) is NO2
       tmp3 = -2.0d0*(R(3)-rNOe)/a1-2.0d0*(1/R(2)-tmp2)*(tmp2-aOOe)/a2
       dVpfdR(3) = Vpf * tmp3

       endif

       end subroutine patchfunc

      subroutine EvdVdR
      use NO2_6Ap_ZV_par
!**********************************************************************
! This subroutine evaluates dVdR for given R 
! dVdR = dV2dR + C*dBdR
! C:            Coefficients, stored in 'dim.inc' 
! P:            Basis functions evaluated for given R
! M:            Monomials evaluated for given R
! dV2dR:        Gradient of 2-body interactions
! dMsdR:        dMsdR(3,3), 3 MEG terms w.r.t. 3 bond lengths
! dMdR:         dMdR(3,nom), nom monomials w.r.t. 3 bond lengths
! dPdR:         dPdR(3,nob), nop polynomial basis functions 
!               w.r.t. 3 bond lengthss
! nom:          number of monomials
! nob:          number of basis functions (polynomials)
! M(nom):       Array to store monomials
! P(nob):       Array to store polynomials
!**********************************************************************
      
      integer i,j
      double precision dist,V,V21,V22,V23
      double precision dv2dr1,dv2dr2,dv2dr3
      double precision Vpf,dVpfdR(3)
! Initialize dVdR(6)
      do i=1,3
        dVdR(i)=0.0d0
      enddo
! Add dV2dR(i) to dVdR
        dist=R(1)
        call ev2gm2(dist,v21,dv2dr1,4,1)  ! OO
        dVdR(1)=dv2dr1
        dist=R(2)
        call ev2gm2(dist,v22,dv2dr2,3,1)  ! NO
        dVdR(2)=dv2dr2
        dist=R(3)
        call ev2gm2(dist,v23,dv2dr3,3,1)  ! NO
        dVdR(3)=dv2dr3

! Calculate dMEG/dr(3,3) for given R(3)
      call evdmsdr
! Calculate the monomials for each point by using 3 MEG terms
      call evdmdr
! Calculate the polynomials by using monomials
      call evdpdr 
! Remove 2-body interactions and unconnected terms from polynomials
      call evdbdr
! Evaluate dVdR(3) by taking the product of C(j) and dPdR(i,j)
      do i=1,3      
        do j=1,227
         dVdR(i)=dVdR(i) + c(j)*dBdR(i,j)
        enddo
      enddo

! Add the derivative of patch function to get the barrier
! between 2A1 and 2B2 minima
!      call patchfunc(Vpf,dVpfdR,1)

!      do i=1,3 
!      dVdR(i) = dVdR(i) + dVpfdR(i)
!      enddo

      return
      end subroutine EvdVdR
      subroutine EvdRdX(X)
      use NO2_6Ap_ZV_par
!**********************************************************************
! Subroutine to evaluate dRdX for given R and X 
! R:            R(3), 3 bond lengths
! X:            X(9), 9 Cartesian coordinates
! 
! dMdR:         dMdR(3,nom), nom monomials w.r.t.3 bond length
! dPdR:         dPdR(3,nob), nop polynomial basis functions 
!               w.r.t. 3 bond lengths
! M(nom):       Array to store monomials
! P(nob):       Array to store polynomials
!**********************************************************************
      integer i,j
      double precision X(9)
! Initialize dRdX(3,9)
      do i=1,3
        do j=1,9
          dRdX(i,j)=0.0d0
        enddo
      enddo
! Start to calculate the non-zero dRdX
! dr1dx
      dRdX(1,1)=(x(1)-x(4))/r(1)
      dRdX(1,2)=(x(2)-x(5))/r(1)
      dRdX(1,3)=(x(3)-x(6))/r(1)
      dRdX(1,4)=-dRdX(1,1)
      dRdX(1,5)=-dRdX(1,2)
      dRdX(1,6)=-dRdX(1,3)
! dr2dx
      dRdX(2,1)=(x(1)-x(7))/r(2)
      dRdX(2,2)=(x(2)-x(8))/r(2)
      dRdX(2,3)=(x(3)-x(9))/r(2)
      dRdX(2,7)=-dRdX(2,1)
      dRdX(2,8)=-dRdX(2,2)
      dRdX(2,9)=-dRdX(2,3)
! dr3dx
      dRdX(3,4)=(x(4)-x(7))/r(3)
      dRdX(3,5)=(x(5)-x(8))/r(3)
      dRdX(3,6)=(x(6)-x(9))/r(3)
      dRdX(3,7)=-dRdX(3,4)
      dRdX(3,8)=-dRdX(3,5)
      dRdX(3,9)=-dRdX(3,6)
! Finish the calculation of non-zero dRdX
      return
      end subroutine EvdRdX
      subroutine EvdMsdR
      use NO2_6Ap_ZV_par
!**********************************************************************
! This subroutine evaluates the derivatives of MEG term X
! w.r.t. interatomic distances R(3)
! dmsdR:        Local variables, dirm(3,3)
! a:            Nonlinear parameter(Angstrom)
! ab:           Nonlinear parameter(Angstrom^2)
! re:           equilibrium bond length(Angstrom)
!**********************************************************************
      integer i,j
! Initialize dmsdr
      do i=1,3
        do j=1,3
          dmsdr(i,j)=0.0d0
        enddo
      enddo
! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
! dmsdr(i,j)=0  i!=j
! dmsdr(i,i)= (-1/a-2*(r-re)/ab)*Exp(-(r-re)/a-(r-re)^2/ab)
      do i=1,3
       dmsdr(i,i)=(-2.0d0*(r(i)-re(i))/ab(i)-1.0d0/a(i)) &
       * Exp(-(r(i)-re(i))/a(i)-((r(i)-re(i))**2.0d0)/ab(i))
      enddo
      return
      end subroutine EvdMsdR
      subroutine EvdMdR
      use NO2_6Ap_ZV_par
!**********************************************************************
!  This subroutine takes M(nom) and dMSdR(3,3) and calculates the
!  dMdR(6,nom) that do not have usable decomposition.
!  For degree 10, the number of monomials is nom.
!**********************************************************************
      integer i
      do i=1,3
      dmdr(i,0) = 1.0d0
      dmdr(i,1) = dmsdr(i,3)
      dmdr(i,2) = dmsdr(i,2)
      dmdr(i,3) = dmsdr(i,1)
      dmdr(i,4) = dmdr(i,1)*rm(2) + rm(1)*dmdr(i,2)
      enddo
      return
      end subroutine EvdMdR
      subroutine EvdPdr
      use NO2_6Ap_ZV_par
!**********************************************************************
!  This subroutine takes monomials(m) and calculates the
!  permutationally invariant polynomials(p)
!  For degree 10, the number of polynomials is nob.
!**********************************************************************
      integer i
      do i=1,3
      dpdr(i,0) = dmdr(i,0)
      dpdr(i,1) = dmdr(i,1) + dmdr(i,2)
      dpdr(i,2) = dmdr(i,3)
      dpdr(i,3) = dmdr(i,4)
      dpdr(i,4) = dpdr(i,2)*p(1) + p(2)*dpdr(i,1)
      dpdr(i,5) = dpdr(i,1)*p(1) + p(1)*dpdr(i,1) -dpdr(i,3) -dpdr(i,3)
      dpdr(i,6) = dpdr(i,2)*p(2) + p(2)*dpdr(i,2)
      dpdr(i,7) = dpdr(i,2)*p(3) + p(2)*dpdr(i,3)
      dpdr(i,8) = dpdr(i,1)*p(3) + p(1)*dpdr(i,3)
      dpdr(i,9) = dpdr(i,2)*p(5) + p(2)*dpdr(i,5)
      dpdr(i,10) = dpdr(i,2)*p(4)  + p(2)*dpdr(i,4)
      dpdr(i,11) = dpdr(i,1)*p(5)  + p(1)*dpdr(i,5) - dpdr(i,8)
      dpdr(i,12) = dpdr(i,2)*p(6)  + p(2)*dpdr(i,6)
      dpdr(i,13) = dpdr(i,2)*p(8)  + p(2)*dpdr(i,8)
      dpdr(i,14) = dpdr(i,2)*p(7)  + p(2)*dpdr(i,7)
      dpdr(i,15) = dpdr(i,3)*p(3)  + p(3)*dpdr(i,3)
      dpdr(i,16) = dpdr(i,3)*p(5)  + p(3)*dpdr(i,5)
      dpdr(i,17) = dpdr(i,2)*p(11) + p(2)*dpdr(i,11)
      dpdr(i,18) = dpdr(i,2)*p(9)  + p(2)*dpdr(i,9)
      dpdr(i,19) = dpdr(i,2)*p(10) + p(2)*dpdr(i,10)
      dpdr(i,20) = dpdr(i,1)*p(11) + p(1)*dpdr(i,11) - dpdr(i,16)
      dpdr(i,21) = dpdr(i,2)*p(12) + p(2)*dpdr(i,12) 
      dpdr(i,22) = dpdr(i,2)*p(15) + p(2)*dpdr(i,15) 
      dpdr(i,23) = dpdr(i,2)*p(16) + p(2)*dpdr(i,16) 
      dpdr(i,24) = dpdr(i,2)*p(13) + p(2)*dpdr(i,13) 
      dpdr(i,25) = dpdr(i,2)*p(14) + p(2)*dpdr(i,14) 
      dpdr(i,26) = dpdr(i,1)*p(15) + p(1)*dpdr(i,15) 
      dpdr(i,27) = dpdr(i,3)*p(11) + p(3)*dpdr(i,11) 
      dpdr(i,28) = dpdr(i,2)*p(20) + p(2)*dpdr(i,20) 
      dpdr(i,29) = dpdr(i,2)*p(17) + p(2)*dpdr(i,17) 
      dpdr(i,30) = dpdr(i,2)*p(18) + p(2)*dpdr(i,18) 
      dpdr(i,31) = dpdr(i,2)*p(19) + p(2)*dpdr(i,19) 
      dpdr(i,32) = dpdr(i,1)*p(20) + p(1)*dpdr(i,20) - dpdr(i,27)
      dpdr(i,33) = dpdr(i,2)*p(21) + p(2)*dpdr(i,21) 
      dpdr(i,34) = dpdr(i,2)*p(26) + p(2)*dpdr(i,26) 
      dpdr(i,35) = dpdr(i,2)*p(27) + p(2)*dpdr(i,27) 
      dpdr(i,36) = dpdr(i,2)*p(22) + p(2)*dpdr(i,22) 
      dpdr(i,37) = dpdr(i,2)*p(23) + p(2)*dpdr(i,23) 
      dpdr(i,38) = dpdr(i,2)*p(24) + p(2)*dpdr(i,24) 
      dpdr(i,39) = dpdr(i,2)*p(25) + p(2)*dpdr(i,25) 
      dpdr(i,40) = dpdr(i,3)*p(15) + p(3)*dpdr(i,15) 
      dpdr(i,41) = dpdr(i,3)*p(16) + p(3)*dpdr(i,16) 
      dpdr(i,42) = dpdr(i,3)*p(20) + p(3)*dpdr(i,20) 
      dpdr(i,43) = dpdr(i,2)*p(32) + p(2)*dpdr(i,32) 
      dpdr(i,44) = dpdr(i,2)*p(28) + p(2)*dpdr(i,28) 
      dpdr(i,45) = dpdr(i,2)*p(29) + p(2)*dpdr(i,29) 
      dpdr(i,46) = dpdr(i,2)*p(30) + p(2)*dpdr(i,30) 
      dpdr(i,47) = dpdr(i,2)*p(31) + p(2)*dpdr(i,31) 
      dpdr(i,48) = dpdr(i,1)*p(32) + p(1)*dpdr(i,32) - dpdr(i,42)
      dpdr(i,49) = dpdr(i,2)*p(33) + p(2)*dpdr(i,33)
      dpdr(i,50) = dpdr(i,2)*p(40) + p(2)*dpdr(i,40)
      dpdr(i,51) = dpdr(i,2)*p(41) + p(2)*dpdr(i,41)
      dpdr(i,52) = dpdr(i,2)*p(42) + p(2)*dpdr(i,42)
      dpdr(i,53) = dpdr(i,2)*p(34) + p(2)*dpdr(i,34)
      dpdr(i,54) = dpdr(i,2)*p(35) + p(2)*dpdr(i,35)
      dpdr(i,55) = dpdr(i,2)*p(36) + p(2)*dpdr(i,36)
      dpdr(i,56) = dpdr(i,2)*p(37) + p(2)*dpdr(i,37)
      dpdr(i,57) = dpdr(i,2)*p(38) + p(2)*dpdr(i,38)
      dpdr(i,58) = dpdr(i,2)*p(39) + p(2)*dpdr(i,39)
      dpdr(i,59) = dpdr(i,1)*p(40) + p(1)*dpdr(i,40)
      dpdr(i,60) = dpdr(i,3)*p(27) + p(3)*dpdr(i,27)
      dpdr(i,61) = dpdr(i,3)*p(32) + p(3)*dpdr(i,32)
      dpdr(i,62) = dpdr(i,2)*p(48) + p(2)*dpdr(i,48)
      dpdr(i,63) = dpdr(i,2)*p(43) + p(2)*dpdr(i,43)
      dpdr(i,64) = dpdr(i,2)*p(44) + p(2)*dpdr(i,44)
      dpdr(i,65) = dpdr(i,2)*p(45) + p(2)*dpdr(i,45)
      dpdr(i,66) = dpdr(i,2)*p(46) + p(2)*dpdr(i,46)
      dpdr(i,67) = dpdr(i,2)*p(47) + p(2)*dpdr(i,47)
      dpdr(i,68) = dpdr(i,1)*p(48) + p(1)*dpdr(i,48) - dpdr(i,61)
      dpdr(i,69) = dpdr(i,2)*p(49) + p(2)*dpdr(i,49) 
      dpdr(i,70) = dpdr(i,2)*p(59) + p(2)*dpdr(i,59) 
      dpdr(i,71) = dpdr(i,2)*p(60) + p(2)*dpdr(i,60) 
      dpdr(i,72) = dpdr(i,2)*p(61) + p(2)*dpdr(i,61) 
      dpdr(i,73) = dpdr(i,2)*p(50) + p(2)*dpdr(i,50) 
      dpdr(i,74) = dpdr(i,2)*p(51) + p(2)*dpdr(i,51) 
      dpdr(i,75) = dpdr(i,2)*p(52) + p(2)*dpdr(i,52) 
      dpdr(i,76) = dpdr(i,2)*p(53) + p(2)*dpdr(i,53) 
      dpdr(i,77) = dpdr(i,2)*p(54) + p(2)*dpdr(i,54) 
      dpdr(i,78) = dpdr(i,2)*p(55) + p(2)*dpdr(i,55) 
      dpdr(i,79) = dpdr(i,2)*p(56) + p(2)*dpdr(i,56) 
      dpdr(i,80) = dpdr(i,2)*p(57) + p(2)*dpdr(i,57) 
      dpdr(i,81) = dpdr(i,2)*p(58) + p(2)*dpdr(i,58) 
      dpdr(i,82) = dpdr(i,3)*p(40) + p(3)*dpdr(i,40) 
      dpdr(i,83) = dpdr(i,3)*p(41) + p(3)*dpdr(i,41) 
      dpdr(i,84) = dpdr(i,3)*p(42) + p(3)*dpdr(i,42) 
      dpdr(i,85) = dpdr(i,3)*p(48) + p(3)*dpdr(i,48) 
      dpdr(i,86) = dpdr(i,2)*p(68) + p(2)*dpdr(i,68) 
      dpdr(i,87) = dpdr(i,2)*p(62) + p(2)*dpdr(i,62) 
      dpdr(i,88) = dpdr(i,2)*p(63) + p(2)*dpdr(i,63) 
      dpdr(i,89) = dpdr(i,2)*p(64) + p(2)*dpdr(i,64) 
      dpdr(i,90) = dpdr(i,2)*p(65) + p(2)*dpdr(i,65) 
      dpdr(i,91) = dpdr(i,2)*p(66) + p(2)*dpdr(i,66) 
      dpdr(i,92) = dpdr(i,2)*p(67) + p(2)*dpdr(i,67) 
      dpdr(i,93) = dpdr(i,1)*p(68) + p(1)*dpdr(i,68) - dpdr(i,85)
      dpdr(i,94) = dpdr(i,2)*p(69) + p(2)*dpdr(i,69)     
      dpdr(i,95) = dpdr(i,2)*p(82) + p(2)*dpdr(i,82)     
      dpdr(i,96) = dpdr(i,2)*p(83) + p(2)*dpdr(i,83)     
      dpdr(i,97) = dpdr(i,2)*p(84) + p(2)*dpdr(i,84)     
      dpdr(i,98) = dpdr(i,2)*p(85) + p(2)*dpdr(i,85)     
      dpdr(i,99) = dpdr(i,2)*p(70) + p(2)*dpdr(i,70)     
      dpdr(i,100) = dpdr(i,2)*p(71)  + p(2)*dpdr(i,71)     
      dpdr(i,101) = dpdr(i,2)*p(72)  + p(2)*dpdr(i,72)     
      dpdr(i,102) = dpdr(i,2)*p(73)  + p(2)*dpdr(i,73)     
      dpdr(i,103) = dpdr(i,2)*p(74)  + p(2)*dpdr(i,74)     
      dpdr(i,104) = dpdr(i,2)*p(75)  + p(2)*dpdr(i,75)     
      dpdr(i,105) = dpdr(i,2)*p(76)  + p(2)*dpdr(i,76)     
      dpdr(i,106) = dpdr(i,2)*p(77)  + p(2)*dpdr(i,77)     
      dpdr(i,107) = dpdr(i,2)*p(78)  + p(2)*dpdr(i,78)     
      dpdr(i,108) = dpdr(i,2)*p(79)  + p(2)*dpdr(i,79)     
      dpdr(i,109) = dpdr(i,2)*p(80)  + p(2)*dpdr(i,80)     
      dpdr(i,110) = dpdr(i,2)*p(81)  + p(2)*dpdr(i,81)     
      dpdr(i,111) = dpdr(i,1)*p(82)  + p(1)*dpdr(i,82)     
      dpdr(i,112) = dpdr(i,3)*p(60)  + p(3)*dpdr(i,60)     
      dpdr(i,113) = dpdr(i,3)*p(61)  + p(3)*dpdr(i,61)     
      dpdr(i,114) = dpdr(i,3)*p(68)  + p(3)*dpdr(i,68)     
      dpdr(i,115) = dpdr(i,2)*p(93)  + p(2)*dpdr(i,93)     
      dpdr(i,116) = dpdr(i,2)*p(86)  + p(2)*dpdr(i,86)     
      dpdr(i,117) = dpdr(i,2)*p(87)  + p(2)*dpdr(i,87)     
      dpdr(i,118) = dpdr(i,2)*p(88)  + p(2)*dpdr(i,88)     
      dpdr(i,119) = dpdr(i,2)*p(89)  + p(2)*dpdr(i,89)     
      dpdr(i,120) = dpdr(i,2)*p(90)  + p(2)*dpdr(i,90)     
      dpdr(i,121) = dpdr(i,2)*p(91)  + p(2)*dpdr(i,91)     
      dpdr(i,122) = dpdr(i,2)*p(92)  + p(2)*dpdr(i,92)     
      dpdr(i,123) = dpdr(i,1)*p(93)  + p(1)*dpdr(i,93) - dpdr(i,114)
      dpdr(i,124) = dpdr(i,2)*p(94)  + p(2)*dpdr(i,94)
      dpdr(i,125) = dpdr(i,2)*p(111) + p(2)*dpdr(i,111)
      dpdr(i,126) = dpdr(i,2)*p(112) + p(2)*dpdr(i,112)
      dpdr(i,127) = dpdr(i,2)*p(113) + p(2)*dpdr(i,113)
      dpdr(i,128) = dpdr(i,2)*p(114) + p(2)*dpdr(i,114)
      dpdr(i,129) = dpdr(i,2)*p(95)  + p(2)*dpdr(i,95)
      dpdr(i,130) = dpdr(i,2)*p(96)  + p(2)*dpdr(i,96)
      dpdr(i,131) = dpdr(i,2)*p(97)  + p(2)*dpdr(i,97)
      dpdr(i,132) = dpdr(i,2)*p(98)  + p(2)*dpdr(i,98)
      dpdr(i,133) = dpdr(i,2)*p(99)  + p(2)*dpdr(i,99)
      dpdr(i,134) = dpdr(i,2)*p(100) + p(2)*dpdr(i,100)
      dpdr(i,135) = dpdr(i,2)*p(101) + p(2)*dpdr(i,101)
      dpdr(i,136) = dpdr(i,2)*p(102) + p(2)*dpdr(i,102)
      dpdr(i,137) = dpdr(i,2)*p(103) + p(2)*dpdr(i,103)
      dpdr(i,138) = dpdr(i,2)*p(104) + p(2)*dpdr(i,104)
      dpdr(i,139) = dpdr(i,2)*p(105) + p(2)*dpdr(i,105)
      dpdr(i,140) = dpdr(i,2)*p(106) + p(2)*dpdr(i,106)
      dpdr(i,141) = dpdr(i,2)*p(107) + p(2)*dpdr(i,107)
      dpdr(i,142) = dpdr(i,2)*p(108) + p(2)*dpdr(i,108)
      dpdr(i,143) = dpdr(i,2)*p(109) + p(2)*dpdr(i,109)
      dpdr(i,144) = dpdr(i,2)*p(110) + p(2)*dpdr(i,110)
      dpdr(i,145) = dpdr(i,3)*p(82)  + p(3)*dpdr(i,82)
      dpdr(i,146) = dpdr(i,3)*p(83)  + p(3)*dpdr(i,83)
      dpdr(i,147) = dpdr(i,3)*p(84)  + p(3)*dpdr(i,84)
      dpdr(i,148) = dpdr(i,3)*p(85)  + p(3)*dpdr(i,85)
      dpdr(i,149) = dpdr(i,3)*p(93)  + p(3)*dpdr(i,93)
      dpdr(i,150) = dpdr(i,2)*p(123) + p(2)*dpdr(i,123)
      dpdr(i,151) = dpdr(i,2)*p(115) + p(2)*dpdr(i,115)
      dpdr(i,152) = dpdr(i,2)*p(116) + p(2)*dpdr(i,116)
      dpdr(i,153) = dpdr(i,2)*p(117) + p(2)*dpdr(i,117)
      dpdr(i,154) = dpdr(i,2)*p(118) + p(2)*dpdr(i,118)
      dpdr(i,155) = dpdr(i,2)*p(119) + p(2)*dpdr(i,119)
      dpdr(i,156) = dpdr(i,2)*p(120) + p(2)*dpdr(i,120)
      dpdr(i,157) = dpdr(i,2)*p(121) + p(2)*dpdr(i,121)
      dpdr(i,158) = dpdr(i,2)*p(122) + p(2)*dpdr(i,122)
      dpdr(i,159) = dpdr(i,1)*p(123) + p(1)*dpdr(i,123) - dpdr(i,149)
      dpdr(i,160) = dpdr(i,2)*p(124) + p(2)*dpdr(i,124)
      dpdr(i,161) = dpdr(i,2)*p(145) + p(2)*dpdr(i,145)
      dpdr(i,162) = dpdr(i,2)*p(146) + p(2)*dpdr(i,146)
      dpdr(i,163) = dpdr(i,2)*p(147) + p(2)*dpdr(i,147)
      dpdr(i,164) = dpdr(i,2)*p(148) + p(2)*dpdr(i,148)
      dpdr(i,165) = dpdr(i,2)*p(149) + p(2)*dpdr(i,149)
      dpdr(i,166) = dpdr(i,2)*p(125) + p(2)*dpdr(i,125)
      dpdr(i,167) = dpdr(i,2)*p(126) + p(2)*dpdr(i,126)
      dpdr(i,168) = dpdr(i,2)*p(127) + p(2)*dpdr(i,127)
      dpdr(i,169) = dpdr(i,2)*p(128) + p(2)*dpdr(i,128)
      dpdr(i,170) = dpdr(i,2)*p(129) + p(2)*dpdr(i,129)
      dpdr(i,171) = dpdr(i,2)*p(130) + p(2)*dpdr(i,130)
      dpdr(i,172) = dpdr(i,2)*p(131) + p(2)*dpdr(i,131)
      dpdr(i,173) = dpdr(i,2)*p(132) + p(2)*dpdr(i,132)
      dpdr(i,174) = dpdr(i,2)*p(133) + p(2)*dpdr(i,133)
      dpdr(i,175) = dpdr(i,2)*p(134) + p(2)*dpdr(i,134)
      dpdr(i,176) = dpdr(i,2)*p(135) + p(2)*dpdr(i,135)
      dpdr(i,177) = dpdr(i,2)*p(136) + p(2)*dpdr(i,136)
      dpdr(i,178) = dpdr(i,2)*p(137) + p(2)*dpdr(i,137)
      dpdr(i,179) = dpdr(i,2)*p(138) + p(2)*dpdr(i,138)
      dpdr(i,180) = dpdr(i,2)*p(139) + p(2)*dpdr(i,139)
      dpdr(i,181) = dpdr(i,2)*p(140) + p(2)*dpdr(i,140)
      dpdr(i,182) = dpdr(i,2)*p(141) + p(2)*dpdr(i,141)
      dpdr(i,183) = dpdr(i,2)*p(142) + p(2)*dpdr(i,142)
      dpdr(i,184) = dpdr(i,2)*p(143) + p(2)*dpdr(i,143)
      dpdr(i,185) = dpdr(i,2)*p(144) + p(2)*dpdr(i,144)
      dpdr(i,186) = dpdr(i,1)*p(145) + p(1)*dpdr(i,145)
      dpdr(i,187) = dpdr(i,3)*p(112) + p(3)*dpdr(i,112)
      dpdr(i,188) = dpdr(i,3)*p(113) + p(3)*dpdr(i,113)
      dpdr(i,189) = dpdr(i,3)*p(114) + p(3)*dpdr(i,114)
      dpdr(i,190) = dpdr(i,3)*p(123) + p(3)*dpdr(i,123)
      dpdr(i,191) = dpdr(i,2)*p(159) + p(2)*dpdr(i,159)
      dpdr(i,192) = dpdr(i,2)*p(150) + p(2)*dpdr(i,150)
      dpdr(i,193) = dpdr(i,2)*p(151) + p(2)*dpdr(i,151)
      dpdr(i,194) = dpdr(i,2)*p(152) + p(2)*dpdr(i,152)
      dpdr(i,195) = dpdr(i,2)*p(153) + p(2)*dpdr(i,153)
      dpdr(i,196) = dpdr(i,2)*p(154) + p(2)*dpdr(i,154)
      dpdr(i,197) = dpdr(i,2)*p(155) + p(2)*dpdr(i,155)
      dpdr(i,198) = dpdr(i,2)*p(156) + p(2)*dpdr(i,156)
      dpdr(i,199) = dpdr(i,2)*p(157) + p(2)*dpdr(i,157)
      dpdr(i,200) = dpdr(i,2)*p(158) + p(2)*dpdr(i,158)
      dpdr(i,201) = dpdr(i,1)*p(159) + p(1)*dpdr(i,159) - dpdr(i,190)
      dpdr(i,202) = dpdr(i,2)*p(160) + p(2)*dpdr(i,160)
      dpdr(i,203) = dpdr(i,2)*p(186) + p(2)*dpdr(i,186)
      dpdr(i,204) = dpdr(i,2)*p(187) + p(2)*dpdr(i,187)
      dpdr(i,205) = dpdr(i,2)*p(188) + p(2)*dpdr(i,188)
      dpdr(i,206) = dpdr(i,2)*p(189) + p(2)*dpdr(i,189)
      dpdr(i,207) = dpdr(i,2)*p(190) + p(2)*dpdr(i,190)
      dpdr(i,208) = dpdr(i,2)*p(161) + p(2)*dpdr(i,161)
      dpdr(i,209) = dpdr(i,2)*p(162) + p(2)*dpdr(i,162)
      dpdr(i,210) = dpdr(i,2)*p(163) + p(2)*dpdr(i,163)
      dpdr(i,211) = dpdr(i,2)*p(164) + p(2)*dpdr(i,164)
      dpdr(i,212) = dpdr(i,2)*p(165) + p(2)*dpdr(i,165)
      dpdr(i,213) = dpdr(i,2)*p(166) + p(2)*dpdr(i,166)
      dpdr(i,214) = dpdr(i,2)*p(167) + p(2)*dpdr(i,167)
      dpdr(i,215) = dpdr(i,2)*p(168) + p(2)*dpdr(i,168)
      dpdr(i,216) = dpdr(i,2)*p(169) + p(2)*dpdr(i,169)
      dpdr(i,217) = dpdr(i,2)*p(170) + p(2)*dpdr(i,170)
      dpdr(i,218) = dpdr(i,2)*p(171) + p(2)*dpdr(i,171)
      dpdr(i,219) = dpdr(i,2)*p(172) + p(2)*dpdr(i,172)
      dpdr(i,220) = dpdr(i,2)*p(173) + p(2)*dpdr(i,173)
      dpdr(i,221) = dpdr(i,2)*p(174) + p(2)*dpdr(i,174)
      dpdr(i,222) = dpdr(i,2)*p(175) + p(2)*dpdr(i,175)
      dpdr(i,223) = dpdr(i,2)*p(176) + p(2)*dpdr(i,176)
      dpdr(i,224) = dpdr(i,2)*p(177) + p(2)*dpdr(i,177)
      dpdr(i,225) = dpdr(i,2)*p(178) + p(2)*dpdr(i,178)
      dpdr(i,226) = dpdr(i,2)*p(179) + p(2)*dpdr(i,179)
      dpdr(i,227) = dpdr(i,2)*p(180) + p(2)*dpdr(i,180)
      dpdr(i,228) = dpdr(i,2)*p(181) + p(2)*dpdr(i,181)
      dpdr(i,229) = dpdr(i,2)*p(182) + p(2)*dpdr(i,182)
      dpdr(i,230) = dpdr(i,2)*p(183) + p(2)*dpdr(i,183)
      dpdr(i,231) = dpdr(i,2)*p(184) + p(2)*dpdr(i,184)
      dpdr(i,232) = dpdr(i,2)*p(185) + p(2)*dpdr(i,185)
      dpdr(i,233) = dpdr(i,3)*p(145) + p(3)*dpdr(i,145)
      dpdr(i,234) = dpdr(i,3)*p(146) + p(3)*dpdr(i,146)
      dpdr(i,235) = dpdr(i,3)*p(147) + p(3)*dpdr(i,147)
      dpdr(i,236) = dpdr(i,3)*p(148) + p(3)*dpdr(i,148)
      dpdr(i,237) = dpdr(i,3)*p(149) + p(3)*dpdr(i,149)
      dpdr(i,238) = dpdr(i,3)*p(159) + p(3)*dpdr(i,159)
      dpdr(i,239) = dpdr(i,2)*p(201) + p(2)*dpdr(i,201)
      dpdr(i,240) = dpdr(i,2)*p(191) + p(2)*dpdr(i,191)
      dpdr(i,241) = dpdr(i,2)*p(192) + p(2)*dpdr(i,192)
      dpdr(i,242) = dpdr(i,2)*p(193) + p(2)*dpdr(i,193)
      dpdr(i,243) = dpdr(i,2)*p(194) + p(2)*dpdr(i,194)
      dpdr(i,244) = dpdr(i,2)*p(195) + p(2)*dpdr(i,195)
      dpdr(i,245) = dpdr(i,2)*p(196) + p(2)*dpdr(i,196)
      dpdr(i,246) = dpdr(i,2)*p(197) + p(2)*dpdr(i,197)
      dpdr(i,247) = dpdr(i,2)*p(198) + p(2)*dpdr(i,198)
      dpdr(i,248) = dpdr(i,2)*p(199) + p(2)*dpdr(i,199)
      dpdr(i,249) = dpdr(i,2)*p(200) + p(2)*dpdr(i,200)
      dpdr(i,250) = dpdr(i,1)*p(201) + p(1)*dpdr(i,201) - dpdr(i,238)
      dpdr(i,251) = dpdr(i,2)*p(202) + p(2)*dpdr(i,202)

      enddo
      return
      end subroutine EvdPdr
      subroutine evdbdr
      use NO2_6Ap_ZV_par
!**********************************************************************
!  This subroutine eliminates the 2-body terms in the permutationally
!  invariant formulation.
!**********************************************************************
      
      integer i,j
      double precision db1dr(3,252) 
! Pass P(0:251) to BM1(1:252)
      do j=1,3
      do i=1,252
        db1dr(j,i)=dpdr(j,i-1)
      enddo
      enddo
! Remove unconnected terms and 2-body terms and pass to B(1:227)
      do j=1,3 
      do i=1,2
        dbdr(j,i)=db1dr(j,i+3)
      enddo
      do i=3,6
        dbdr(j,i)=db1dr(j,i+5)
      enddo
      do i=7,13  
        dbdr(j,i)=db1dr(j,i+7)
      enddo
      do i=14,23
        dbdr(j,i)=db1dr(j,i+9)
      enddo
      do i=24,37
        dbdr(j,i)=db1dr(j,i+11)  
      enddo
      do i=38,55
        dbdr(j,i)=db1dr(j,i+13)
      enddo
      do i=56,78
        dbdr(j,i)=db1dr(j,i+15)
      enddo
      do i=79,106
        dbdr(j,i)=db1dr(j,i+17)     
      enddo
      do i=107,140
        dbdr(j,i)=db1dr(j,i+19)
      enddo
      do i=141,180
        dbdr(j,i)=db1dr(j,i+21)
      enddo
      do i=181,227
        dbdr(j,i)=db1dr(j,i+23)
      enddo
      enddo
  
      return
      end subroutine evdbdr
      subroutine d3disp(dist,disp,dispdr,igrad,imol)
!**********************************************************************
! Dispersion correction based on Grimme's D3(BJ) calculation for
! diatomic pairs
!
! Several subroutines of DFTD3 V3.1 Rev 1 by Grimme were merged into 
! subroutine edisp and they have been heavily modified to calculate
! only dispersion energy corrections that are needed.
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
! and 
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011),
! 1456-1465
!
! The C6 values are fixed.
!
!**********************************************************************
      double precision cn(2),s6,s8,rs6,rs8
      double precision dist(1), e6(1), e8(1), disp, dispdr(1), c6(2)
      double precision e6dr(1),e8dr(1)
      integer iz(2), mxc(94), i, j, igrad
      double precision c6ab(94,94,5,5,3)
      double precision r2r4(94)
      double precision autoang,autokcal
      integer imol
      autoang =0.52917726d0
      autokcal=627.509541d0
! Generalized parameters for BJ damping from P. Verma, B. Wang, 
! L. E. Fernandez, and D. G. Truhlar, J. Phys. Chem. A 121, 2855 (2017)
      s6= 1.0d0
      s8= 2.0d0
      rs6= 0.5299d0
      rs8= 2.20d0
      do i=1,1
      dist(i)=dist(i)/autoang
      enddo
      if (imol.eq.4) then
! iz for O2 system
      iz(1)=8
      iz(2)=8
! C6 for O2 system
      c6(1)=12.8d0
      c6(2)=12.8d0
      else if (imol.eq.3) then
! iz for NO system
      iz(1)=8
      iz(2)=7
      c6(1)=16.7d0  !C6 for NO doubly-coordinated
      c6(2)=16.7d0  !C6 for NO doubly-coordinated
      else
! currently imol = 4 is used
      stop
      endif
! Calculate dispersion correction
      call edisp(94,5,2,dist,iz,mxc, &
           rs6,rs8,e6,e8,e6dr,e8dr,c6,0)
      disp = 0.0d0
      do i=1,1
      disp =disp + (-s6*e6(i)-s8*e8(i))*autokcal
      enddo
      if (igrad .eq. 1) then
      call edisp(94,5,2,dist,iz,mxc, &
           rs6,rs8,e6,e8,e6dr,e8dr,c6,1)
      dispdr(:) = 0.0d0
      do i=1,1
      dispdr(i) =dispdr(i) + (-s6*e6dr(i)-s8*e8dr(i))*autokcal/autoang
      enddo
      endif
      do i=1,1
      dist(i)=dist(i)*autoang
      enddo
      end subroutine d3disp
!**********************************************************************
! compute energy
!**********************************************************************
      subroutine edisp(max_elem,maxc,n,dist,iz,mxc, &
                 rs6,rs8,e6,e8,e6dr,e8dr,c6a,igrad)
      integer n,iz(2),max_elem,maxc,mxc(max_elem) 
      double precision dist(1),r2r4(max_elem),r0ab(max_elem,max_elem)
      double precision rs6,rs8,rcov(max_elem)
      double precision c6ab(max_elem,max_elem,maxc,maxc,3)
      double precision e6(1), e8(1), c6a(2), e6dr(1), e8dr(1)
       
      integer iat,jat,igrad
      double precision r,tmp,c6,c8,a1,a2
      double precision damp6,damp8
      double precision cn(n)
      double precision r2ab(n*n),cc6ab(n*n),dmp(n*n)
      integer step
      e6(:) =0.0d0
      e8(:) =0.0d0
      e6dr(:) =0.0d0
      e8dr(:) =0.0d0
      a1=rs6
      a2=rs8 
!  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
!  the large number of digits is just to keep the results consistent
!  with older versions. They should not imply any higher accuracy than
!  the old values
      r2r4(1:94)=(/ &
      2.00734898d0,  1.56637132d0,  5.01986934d0,  3.85379032d0, &
      3.64446594d0,  3.10492822d0,  2.71175247d0,  2.59361680d0, &
      2.38825250d0,  2.21522516d0,  6.58585536d0,  5.46295967d0, &
      5.65216669d0,  4.88284902d0,  4.29727576d0,  4.04108902d0, &
      3.72932356d0,  3.44677275d0,  7.97762753d0,  7.07623947d0, &
      6.60844053d0,  6.28791364d0,  6.07728703d0,  5.54643096d0, &
      5.80491167d0,  5.58415602d0,  5.41374528d0,  5.28497229d0, &
      5.22592821d0,  5.09817141d0,  6.12149689d0,  5.54083734d0, &
      5.06696878d0,  4.87005108d0,  4.59089647d0,  4.31176304d0, &
      9.55461698d0,  8.67396077d0,  7.97210197d0,  7.43439917d0, &
      6.58711862d0,  6.19536215d0,  6.01517290d0,  5.81623410d0, &
      5.65710424d0,  5.52640661d0,  5.44263305d0,  5.58285373d0, &
      7.02081898d0,  6.46815523d0,  5.98089120d0,  5.81686657d0, &
      5.53321815d0,  5.25477007d0, 11.02204549d0,  0.15679528d0, &
      9.35167836d0,  9.06926079d0,  8.97241155d0,  8.90092807d0, &
      8.85984840d0,  8.81736827d0,  8.79317710d0,  7.89969626d0, &
      8.80588454d0,  8.42439218d0,  8.54289262d0,  8.47583370d0, &
      8.45090888d0,  8.47339339d0,  7.83525634d0,  8.20702843d0, &
      7.70559063d0,  7.32755997d0,  7.03887381d0,  6.68978720d0, &
      6.05450052d0,  5.88752022d0,  5.70661499d0,  5.78450695d0, &
      7.79780729d0,  7.26443867d0,  6.78151984d0,  6.67883169d0, &
      6.39024318d0,  6.09527958d0, 11.79156076d0, 11.10997644d0, &
      9.51377795d0,  8.67197068d0,  8.77140725d0,  8.65402716d0, &
      8.53923501d0,  8.85024712d0 /)
! these new data are scaled with k2=4./3. and converted to a_0 via
! autoang=0.52917726d0
      rcov(1:94)=(/ &
       0.80628308d0, 1.15903197d0, 3.02356173d0, 2.36845659d0, &
       1.94011865d0, 1.88972601d0, 1.78894056d0, 1.58736983d0, &
       1.61256616d0, 1.68815527d0, 3.52748848d0, 3.14954334d0, &
       2.84718717d0, 2.62041997d0, 2.77159820d0, 2.57002732d0, &
       2.49443835d0, 2.41884923d0, 4.43455700d0, 3.88023730d0, &
       3.35111422d0, 3.07395437d0, 3.04875805d0, 2.77159820d0, &
       2.69600923d0, 2.62041997d0, 2.51963467d0, 2.49443835d0, &
       2.54483100d0, 2.74640188d0, 2.82199085d0, 2.74640188d0, &
       2.89757982d0, 2.77159820d0, 2.87238349d0, 2.94797246d0, &
       4.76210950d0, 4.20778980d0, 3.70386304d0, 3.50229216d0, &
       3.32591790d0, 3.12434702d0, 2.89757982d0, 2.84718717d0, &
       2.84718717d0, 2.72120556d0, 2.89757982d0, 3.09915070d0, &
       3.22513231d0, 3.17473967d0, 3.17473967d0, 3.09915070d0, &
       3.32591790d0, 3.30072128d0, 5.26603625d0, 4.43455700d0, &
       4.08180818d0, 3.70386304d0, 3.98102289d0, 3.95582657d0, &
       3.93062995d0, 3.90543362d0, 3.80464833d0, 3.82984466d0, &
       3.80464833d0, 3.77945201d0, 3.75425569d0, 3.75425569d0, &
       3.72905937d0, 3.85504098d0, 3.67866672d0, 3.45189952d0, &
       3.30072128d0, 3.09915070d0, 2.97316878d0, 2.92277614d0, &
       2.79679452d0, 2.82199085d0, 2.84718717d0, 3.32591790d0, &
       3.27552496d0, 3.27552496d0, 3.42670319d0, 3.30072128d0, &
       3.47709584d0, 3.57788113d0, 5.06446567d0, 4.56053862d0, &
       4.20778980d0, 3.98102289d0, 3.82984466d0, 3.85504098d0, &
       3.88023730d0, 3.90543362d0 /)
! DFT-D3
      step=0
      do iat=1,n-1
         do jat=iat+1,n
         step=step+1
         r=dist(step)
         c6=c6a(step)
! r2r4 stored in main as sqrt
         c8 =3.0d0*c6*r2r4(iz(iat))*r2r4(iz(jat))
! energy for BJ damping
          tmp=sqrt(c8/c6)
          e6(step)= c6/(r**6+(a1*tmp+a2)**6)
          e8(step)= c8/(r**8+(a1*tmp+a2)**8)
! calculate gradients
         if (igrad .eq. 1) then
! grad for BJ damping
          e6dr(step)=c6*(-6*r**5)/(r**6+(a1*tmp+a2)**6)**2
          e8dr(step)=c8*(-8*r**7)/(r**8+(a1*tmp+a2)**8)**2
         endif
         enddo
      enddo
      end subroutine edisp
  

