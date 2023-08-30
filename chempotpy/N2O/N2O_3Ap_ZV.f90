!********************************************************************** 
!   System:                     N2O (lowest 3A' state)
!   Name of this surface:       N2O_3Ap_MB-PIP-MEG2 
!   Functional form:            permutationally invariant polynomials
!   Common name:                N2O-3Ap (adiabatic triplet ground state)
!   Number of derivatives:      1
!   Number of bodies:           3
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References: W. Lin, Z. Varga, G. Song, Y. Paukku, and D. G. Truhlar, Global Triplet 
!               Potential Energy Surfaces for the N2(X 1Σ) + O(3P) → NO(X 2Π) + N(4S) 
!               Reaction, J. Chem. Phys. 2016, 144, 024309.
!               Z. Varga, W. Lin, G. Song, Y. Paukku, and D. G. Truhlar, 
!               MEG2 surfaces for N2O,
!               in: POTLIB: An Online Library of Potential Energy Surfaces,
!               comp.chem.umn.edu/potlib/ 
!
!   Notes:    -PES of N2O expressed as sum of pairwise and many-body potentials 
!             -Mixed-exponential-Gaussian (MEG) variables are applied
!             -This is a modified version of a potential originally published in 2016. 
!              The difference between this surface and the published surface is that 
!              the pairwise potentials of NO and N2 are replaced with new ones that
!              include dispersion terms, After this replacement, the least squares 
!              fits of the three-body potential were redetermined. With these changes
!              to the pairwise potentials, we have a unified set of potentials for 
!              systems containing N2, O2, and NO that all use the same diatomic 
!              potentials:
!                 N4_1A_MB-PIP-MEG3              
!                 N2O2_3A_MB-PIP-MEG2
!                 N2O_3Ap_MB-PIP-MEG2
!                 N2O_3App_MB-PIP-MEG2
!                 PES_O3_1_1Ap_umn_v1
!                 PES_O3_1_1App_umn_v1
!                 PES_O3_1_3Ap_umn_v1
!                 PES_O3_1_3App_umn_v1
!                 PES_O3_1_5Ap_umn_v1
!                 PES_O3_1_5App_umn_v1
!                 PES_O3_2_1Ap_umn_v1
!                 PES_O3_2_3Ap_umn_v1
!                 PES_O3_2_5Ap_umn_v1
!                 PES_O4_quintet_umn_v1
!                 PES_O4_singlet_umn_v1
!                 PES_O4_triplet_umn_v2
!                 NO2_2Ap_MB-PIP-MEG
!                 NO2_4Ap_MB-PIP-MEG
!
!     O1
!
!     N2--N3
!
!   Input: X(3),Y(3),Z(3)               in units of bohr
!   Output: E                           in units of hartree
!   Output: dEdX(3),dEdY(3),dEdZ(3)     hartree/bohr
!**********************************************************************

      module N2O_3Ap_ZV_par

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
! R(6):         Interatomic bond distance
! rMs(6)                Array to store base terms
! rM(0:17)              Array to store monomials
! P(0:2253)     Array to store polynomials
! B(1:256)      Array to store basis functions
! B2(1:2153)    Array to store basis functions
! dMsdR(6,6):   The derivative of base terms w.r.t. R
! dMdR(6,18):   The derivative of monomials w.r.t. R
! dPdR(6,2254): The derivative of basis functions w.r.t. R
! dVdR(6):      The derivative of V w.r.t. R
! dRdX(6,12):   The derivative of R w.r.t. X
! dB2dR(6,2153) The derivative of B w.r.t. R (four-body) 
! dBdR(6,256)   The derivative of B w.r.t. R (three-body)

      double precision :: R(6)
      double precision :: rMs(6),rM(0:17),P(0:2253),B(256),B2(2153)
      double precision :: dMsdR(6,6),dMdR(6,0:17),dPdR(6,0:2253)
      double precision :: dVdR(6),dRdX(6,12)
      double precision :: dB2dR(6,2153),dBdR(6,256)
  
! Nonlinear parameters:
! a(in Ang)
! ab (in Ang^2)
! re (in Ang)
! order of data: OO,NO,NO,NO,NO,NN

      double precision,parameter :: a(6)  = (/ 1.30d0,1.20d0,1.20d0,&
      1.20d0,1.20d0,0.90d0 /)
      double precision,parameter :: ab(6) = (/ 1.60d0,2.35d0,2.35d0,&
      2.35d0,2.35d0,1.90d0 /)
      double precision,parameter :: re(6) = (/ 1.208d0,1.1508d0,&
      1.1508d0,1.1508d0,1.1508d0,1.098d0 /)

! Reference energy of N2 + O in hartree for N2O system
      double precision,parameter :: Eref = -0.191619504727d0

! De(N2) + De(O2) = 228.4 + 120.143 kcal/mol 
      double precision,parameter :: totdiss = 348.643d0

! Linear parameters optimized by the weighted-least square fitting
      double precision,parameter :: C(256) = (/ &
        -0.215637524806D+01 &
      ,  0.170177729501D+04 &
      , -0.334328137244D+03 &
      , -0.132128990130D+03 &
      ,  0.762400190155D+03 &
      ,  0.562983055088D+03 &
      , -0.407913183643D+04 &
      , -0.805679763552D+02 &
      ,  0.349030883065D+03 &
      , -0.832830465631D+03 &
      ,  0.652844051811D+02 &
      , -0.231072300689D+03 &
      ,  0.668519029794D+03 &
      ,  0.596879158350D+04 &
      ,  0.232079769136D+03 &
      ,  0.339861527808D+04 &
      , -0.125942963093D+03 &
      , -0.306475569790D+03 &
      ,  0.543964293001D+03 &
      , -0.294212929184D+03 &
      ,  0.164025139190D+03 &
      , -0.455728518169D+02 &
      ,  0.257707957452D+03 &
      ,  0.359373572165D+02 &
      ,  0.165231417288D+03 &
      , -0.134846966931D+03 &
      ,  0.221227416474D+05 &
      , -0.216907037120D+02 &
      ,  0.693255697527D+02 &
      ,  0.371936995632D+02 &
      ,  0.124720175332D+05 &
      , -0.438040415304D+03 &
      , -0.147957885915D+03 &
      ,  0.266292554005D-02 &
      ,  0.109083502366D-03 &
      , -0.746672792502D-04 &
      ,  0.780773620832D-03 &
      , -0.218205400597D-02 &
      , -0.374398534802D-03 &
      , -0.208234417329D-02 &
      ,  0.873746200796D-04 &
      ,  0.284916291776D-03 &
      , -0.136042459027D+05 &
      , -0.604390680752D+04 &
      ,  0.245669797419D-02 &
      , -0.206602274421D+05 &
      , -0.427000824601D+05 &
      , -0.235225165528D-03 &
      , -0.120830073612D+05 &
      , -0.107708406233D-03 &
      , -0.664409491606D-03 &
      ,  0.171717312992D-02 &
      ,  0.222027608993D-02 &
      ,  0.142345706763D-03 &
      , -0.844662096512D-03 &
      , -0.533704859890D-03 &
      , -0.762619392390D-03 &
      ,  0.996582230869D-03 &
      ,  0.590031504544D-05 &
      , -0.558084852173D-03 &
      ,  0.275761571629D-03 &
      ,  0.305085927721D-04 &
      , -0.873676802939D-04 &
      ,  0.814120783616D-03 &
      ,  0.458260993184D-03 &
      , -0.393196796722D-05 &
      , -0.936248534344D-05 &
      , -0.229253646467D-05 &
      , -0.203262061405D-06 &
      , -0.413652651332D-05 &
      ,  0.412392409999D-05 &
      ,  0.841946530272D-06 &
      , -0.586336682318D-06 &
      , -0.106134007183D-05 &
      , -0.277322851616D-05 &
      , -0.111854358965D+06 &
      ,  0.190947048395D-05 &
      , -0.259601620201D-05 &
      , -0.219575622395D-05 &
      ,  0.212981649383D-05 &
      ,  0.381943415340D+04 &
      , -0.136990092869D-05 &
      ,  0.183395968634D-05 &
      ,  0.426426322513D+05 &
      , -0.123176505440D-05 &
      ,  0.258766885963D-05 &
      , -0.236030246015D-06 &
      ,  0.399367490900D-05 &
      , -0.917758084689D+05 &
      ,  0.159470801009D-07 &
      ,  0.150499545271D-06 &
      ,  0.122480196296D-05 &
      , -0.185886165127D-06 &
      , -0.464075128548D-06 &
      ,  0.258754880633D-06 &
      , -0.179657945409D-06 &
      ,  0.131958950078D-05 &
      ,  0.142645149026D-07 &
      ,  0.165728124557D-07 &
      , -0.120926415548D-07 &
      , -0.593536242377D-07 &
      ,  0.574800651520D-08 &
      ,  0.222826201934D-08 &
      , -0.672298483551D-08 &
      , -0.257559804595D-07 &
      ,  0.122963683680D-07 &
      , -0.149884726852D-07 &
      ,  0.100408215076D-08 &
      ,  0.839054337121D-08 &
      ,  0.642761229033D+05 &
      ,  0.635120720311D+04 &
      ,  0.589787785857D+05 &
      , -0.630006979918D-08 &
      ,  0.520613207365D+05 &
      ,  0.131371713319D+06 &
      ,  0.873114913702D-08 &
      ,  0.317435223145D+05 &
      , -0.261934474111D-08 &
      , -0.236832420342D-08 &
      ,  0.109139364213D-08 &
      , -0.419095158577D-08 &
      ,  0.634463503957D-08 &
      , -0.165891833603D-08 &
      ,  0.478758011013D-08 &
      ,  0.107684172690D-08 &
      , -0.186264514923D-08 &
      , -0.547152012587D-08 &
      ,  0.576255843043D-08 &
      , -0.125728547573D-07 &
      ,  0.909494701773D-10 &
      , -0.785803422332D-09 &
      ,  0.582076609135D-10 &
      , -0.349245965481D-09 &
      ,  0.116415321827D-09 &
      ,  0.582076609135D-10 &
      , -0.303771230392D-09 &
      ,  0.232830643654D-09 &
      , -0.162981450558D-08 &
      , -0.145519152284D-10 &
      ,  0.316805724908D+06 &
      , -0.360698009971D+05 &
      ,  0.613935636317D+05 &
      , -0.126152816824D+05 &
      , -0.642014938204D+05 &
      ,  0.234475692736D+06 &
      , -0.146237994102D+06 &
      , -0.520314656888D+05 &
      ,  0.202398756086D+05 &
      , -0.159181343091D+06 &
      , -0.189387423144D+06 &
      , -0.259528146035D+05 &
      , -0.206763994133D+06 &
      , -0.812958974635D+05 &
      , -0.487597367512D+06 &
      ,  0.400336496538D+05 &
      , -0.253960651064D+06 &
      ,  0.201276664540D+06 &
      ,  0.316031421509D+05 &
      , -0.623711285647D+05 &
      , -0.104298666604D+06 &
      ,  0.104439580609D+06 &
      , -0.293089804249D+06 &
      ,  0.186566144966D+06 &
      ,  0.111313353132D+06 &
      , -0.569178251024D+05 &
      ,  0.168099315973D+05 &
      ,  0.176580373226D+06 &
      ,  0.282441991193D+06 &
      , -0.694599483482D+05 &
      ,  0.171624526072D+06 &
      ,  0.147109150846D+06 &
      ,  0.418501059972D+06 &
      , -0.339378337678D+05 &
      ,  0.226664473506D+06 &
      , -0.150564246133D+06 &
      ,  0.317189187012D+05 &
      , -0.374156508526D+06 &
      , -0.930035749069D+05 &
      ,  0.258693105922D+06 &
      ,  0.694852928398D+05 &
      , -0.174574143041D+05 &
      , -0.148280576710D+06 &
      ,  0.252923430828D+06 &
      , -0.135748760155D+06 &
      , -0.108716943728D+06 &
      ,  0.339476060153D+05 &
      ,  0.400675616645D+05 &
      , -0.193631675600D+06 &
      , -0.388747604464D+05 &
      , -0.615304152132D+05 &
      , -0.209303970256D+06 &
      ,  0.105030369332D+06 &
      , -0.721858974389D+05 &
      , -0.152213264036D+06 &
      , -0.188758134292D+06 &
      ,  0.226697141272D+05 &
      , -0.610670025761D+05 &
      ,  0.477026125770D+04 &
      , -0.435020516049D+05 &
      , -0.225860198325D+06 &
      ,  0.128759234888D+06 &
      ,  0.179038356226D+06 &
      ,  0.225731682725D+06 &
      , -0.245425281837D+06 &
      , -0.975225583205D+05 &
      , -0.416541023501D+06 &
      ,  0.107644827458D+06 &
      ,  0.210237752431D+06 &
      ,  0.917124604548D+05 &
      , -0.158252856477D+06 &
      ,  0.525924526149D+05 &
      ,  0.512679911890D+05 &
      , -0.268930739542D+04 &
      , -0.300267708686D+05 &
      ,  0.478820084516D+04 &
      ,  0.742695904006D+05 &
      ,  0.126679196149D+05 &
      , -0.126907834795D+05 &
      ,  0.691941630835D+05 &
      , -0.533603567569D+05 &
      ,  0.140360785573D+05 &
      ,  0.809317128910D+05 &
      ,  0.345834633437D+05 &
      , -0.567982415302D+04 &
      , -0.760873197940D+04 &
      ,  0.183763229251D+05 &
      ,  0.747743834306D+04 &
      ,  0.435252002638D+05 &
      , -0.341671324734D+05 &
      , -0.431898598057D+04 &
      ,  0.101209763749D+05 &
      , -0.310960177307D+05 &
      , -0.669851242060D+04 &
      , -0.607679514479D+05 &
      , -0.436626546076D+05 &
      ,  0.441578657660D+05 &
      ,  0.413403948250D+05 &
      ,  0.160573840836D+06 &
      , -0.573082659981D+05 &
      , -0.103848529745D+06 &
      , -0.202475208030D+05 &
      ,  0.458642503038D+05 &
      , -0.843491821560D+04 &
      , -0.928375739302D+04 &
      , -0.224256375528D+04 &
      ,  0.261372647101D+04 &
      ,  0.108511856405D+05 &
      , -0.405988157452D+05 &
      , -0.122420031626D+05 &
      ,  0.970440619411D+04 &
      ,  0.516951274408D+04 &
      ,  0.231547028433D+04 &
      , -0.700715903307D+04 &
      ,  0.102938699829D+05 &
      , -0.127297896327D+04 &
      , -0.172828826899D+05 &
      /)

      end module N2O_3Ap_ZV_par

      subroutine pes(x,igrad,potential,gradient,dvec)

      use N2O_3Ap_ZV_par
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad
      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: potential(nstates)
      double precision, intent(out) :: gradient(nstates,natoms,3)
      double precision, intent(out) :: dvec(nstates,nstates,natoms,3)


      double precision :: v, tx(12), tg(12)
      integer :: iatom, idir, j, istate
      !initialize 
      potential=0.d0
      gradient=0.d0
      dvec=0.d0

      ! Notice that this N2O surface is weird, it is using N2O2 surface
      ! and simply make the first three coordinates as zero. 
      do idir=1,3
        tx(idir)=1.0d10
      enddo
      do iatom=1,3
        do idir=1,3
          j=3*(iatom-1)+idir+3
          tx(j)=x(iatom, idir)
        enddo
      enddo
      tg=0.d0

      call n2o2pes(tx, v, tg, igrad)

      ! output v is in kcal/mol, but there is a reference energy that is
      ! in hartree, so convert everything in eV.
      v=v/23.0609 !+ Eref*27.211386d0
      
      ! output g is in kcal/mol/ang
      tg=tg/23.0609

      do istate=1,nstates
        potential(istate)=v
      enddo

      do istate=1,nstates
        do iatom=1,natoms
          do idir=1,3
            j=3*(iatom-1)+idir+3
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
!                          N4: X(10),X(11),X(12)
! input flag: igrad     igrad=0 energy calculation only
!                       igrad=1 calculation of energy and gradient
! output potential energy:      v    in kcal/mol
! output gradient:              dVdX in kcal/(mol*Ang)
!********************************************************************** 

      subroutine n2o2pes(X,v,dVdX,igrad)
      use N2O_3Ap_ZV_par
!**********************************************************************
! Subroutine to calculate the potential energy V and gradient dVdX
! for given Cartesian coordinates X(12)  
! R:            Interatomic bond distances (6)
! V:            Calculated potential energy
! dVdX:         The derivative of V w.r.t. X, dim(12)
! dVdR:         The derivative of V w.r.t. R, dim(6) 
! dPdR:         The derivative of basis functions w.r.t. R
!               dim(6*2254)
! dMdR:         The derivative of monomials w.r.t. R
!               dim(6*18)
! dRdX:         The derivative of R w.r.t. X, dim(6*12)
!**********************************************************************

      integer i,igrad,j,nob,k
      double precision V
      double precision dVdX(12),X(12) 

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
        write (*,*) 'Only energy and gradient are available'
      endif
      end subroutine n2o2pes

      subroutine coord_convt(X)
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine calculates the six interatomic distances 
!  from the XYZ coordinates
!
!  r1 = r(O1N2)    r2 = r(O1N3)
!  r3 = r(O1N4)    r4 = r(O2N3)
!  r5 = r(O2N4)    r6 = r(N3N4)     
!**********************************************************************

      integer i
      double precision X(12)

      R(1)=Sqrt((X(4)-X(1))**2  + (X(5)-X(2))**2  + (X(6)-X(3))**2)
      R(2)=Sqrt((X(7)-X(1))**2  + (X(8)-X(2))**2  + (X(9)-X(3))**2)
      R(3)=Sqrt((X(10)-X(1))**2 + (X(11)-X(2))**2 + (X(12)-X(3))**2)
      R(4)=Sqrt((X(4)-X(7))**2  + (X(5)-X(8))**2  + (X(6)-X(9))**2)
      R(5)=Sqrt((X(4)-X(10))**2 + (X(5)-X(11))**2 + (X(6)-X(12))**2)
      R(6)=Sqrt((X(7)-X(10))**2 + (X(8)-X(11))**2 + (X(9)-X(12))**2)
 
      return
      end subroutine coord_convt

      subroutine EvV(V)
      use N2O_3Ap_ZV_par
!**********************************************************************
! This subroutine evaluates V for given R 
! V(R) = C*P
! C:            Coefficients, stored in 'dim.inc' 
! P:            Basis functions evaluated for given R
! rMs:          rMs(6), six mixed exponential Gaussian (MEG) terms
! a:            Nonlinear parameters in MEG terms (Angstrom)
! ab:           Nonlinear parameters in MEG terms (Angstrom^2)
! re:           Equilibrium bond length (Angstrom)
! nop:          number of points
! nom:          number of monomials
! nob:          number of basis functions (polynomials)
! rM(0:17):     Array to store monomials
! P(0:2253):    Array to store polynomials
! B(1:2153):    Array to store basis functions
!**********************************************************************

      integer i,j,k
      double precision dist,dv2dr,V,V21,V22,V23,V24,V25,V26

! Calculate the six MEG terms for each point
      call evmeg

! Calculate the monomials for each point by using six MEG terms
      call evmono

! Calculate the polynomials (basis functions) by using monomials
      call evpoly 

! Calculate the basis functions by removing unconnected and 2-body terms
      call evbas

! Eliminate the unused four-body terms
      call evbas2

! Initialized v to be De(N2) + De(O2) = 228.4 + 120.143 kcal/mol
      v=totdiss
! Evaluate 2-body interactions
        dist=r(1)
        call ev2gm2(dist,v21,dv2dr,4,0)         ! OO
        dist=r(2)
        call ev2gm2(dist,v22,dv2dr,3,0)         ! NO
        dist=r(3)
        call ev2gm2(dist,v23,dv2dr,3,0)         ! NO
        dist=r(4)
        call ev2gm2(dist,v24,dv2dr,3,0)         ! NO
        dist=r(5)
        call ev2gm2(dist,v25,dv2dr,3,0)         ! NO
        dist=r(6)
        call ev2gm2(dist,v26,dv2dr,2,0)         ! NN

        v= v + v21 + v22 + v23 + v24 + v25 + v26

! Evaluate V by taken the product of C and Basis function arrays
      do i=1,256
        v=v + c(i)*b(i)
      enddo

!      Write(*,9999) V 
! 9999 Format('The potential energy is ',F20.14,' kcal/mol')

      return
      end subroutine EvV

      subroutine EvdVdX(X,dVdX)
      use N2O_3Ap_ZV_par
!**********************************************************************
! This subroutine evaluates dRdX for given R and X 
! R:            R(6), 6 bond lengths
! X:            X(12), 12 Cartesian coordinates
! rM(0:17):    Array to store monomials
! P(0:2253):   Array to store polynomials
! dVdX:         dVdX(12), derivatives of V w.r.t. Cartesian coordinates 
! dVdR:         dVdR(6), derivatives of V w.r.t.6 bond length
! dRdX:         dRdX(6,12), derivatives of R(6) w.r.t. 12  
!               Cartesian coordinates
!**********************************************************************

      integer i,j
      double precision dVdX(12),X(12)

! Initialize dVdX
      do i=1,12
        dVdX(i)=0.0d0
      enddo

! Call EvdVdR to evaluate dVdR(6)
      Call evdvdr

! Call EvdRdX to evaluate dRdX(6,12)
      Call evdrdx(X)  

! Calculate dVdX by using chain rule: dV/dXi=(dV/dRj)*(dRj/dXi),
! j=1 to 6
      do i=1,12
        do j=1,6
          dVdX(i)=dVdX(i) + dVdR(j)*dRdX(j,i)
        enddo
      enddo

      return
      end subroutine EvdVdX

      subroutine EvMEG
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine calculates the MEG function
!  MEG term rms = exp(-(r-re)/a-(r-re)^2/ab)
!  re:   equilibrium bond length (Angstrom)
!  a:    nonlinear parameter (Angstrom)
!  ab:   nonlinear parameter (Angstrom^2)
!**********************************************************************

      integer i

      do i=1,6
         rms(i)=Exp(-(r(i)-re(i))/a(i)-(r(i)-re(i))**2.0d0/ab(i))
      enddo
      end subroutine EvMEG

      subroutine EvMono
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine takes six MEG variables(X) and calculates the
!  monomials(RM) and that do not have usable decomposition.
!  For degree 10, the number of monomials is nom.
!**********************************************************************

      rm(0) = 1.0d0
      rm(1) = rms(6)
      rm(2) = rms(5)
      rm(3) = rms(4)
      rm(4) = rms(3)
      rm(5) = rms(2)
      rm(6) = rms(1)
      rm(7) = rm(3)*rm(4)
      rm(8) = rm(2)*rm(5)
      rm(9) = rm(2)*rm(4)
      rm(10) = rm(3)*rm(5)
      rm(11) = rm(2)*rm(3)
      rm(12) = rm(4)*rm(5)
      rm(13) = rm(2)*rm(7)
      rm(14) = rm(2)*rm(10)
      rm(15) = rm(2)*rm(12)
      rm(16) = rm(3)*rm(12)
      rm(17) = rm(2)*rm(16)

      return

      end subroutine EvMono

      subroutine EvPoly
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine takes monomials(RM) and calculates the
!  permutationally invariant polynomials(p)
!  For degree 10, the number of polynomials is nob.
!**********************************************************************

      p(0) = rm(0)
      p(1) = rm(1)
      p(2) = rm(2) + rm(3) + rm(4) + rm(5)
      p(3) = rm(6)
      p(4) = p(1)*p(2)
      p(5) = rm(7) + rm(8)
      p(6) = rm(9) + rm(10)
      p(7) = rm(11) + rm(12)
      p(8) = p(1)*p(3)
      p(9) = p(3)*p(2)
      p(10) = p(1)*p(1)
      p(11) = p(2)*p(2) - p(7) - p(6) - p(5) - p(7) &
            - p(6) - p(5)
      p(12) = p(3)*p(3)
      p(13) = p(1)*p(5)
      p(14) = p(1)*p(6)
      p(15) = p(1)*p(7)
      p(16) = rm(13) + rm(14) + rm(15) + rm(16)
      p(17) = p(1)*p(9)
      p(18) = p(3)*p(5)
      p(19) = p(3)*p(6)
      p(20) = p(3)*p(7)
      p(21) = p(1)*p(4)
      p(22) = p(1)*p(11)
      p(23) = p(2)*p(5) - p(16)
      p(24) = p(2)*p(6) - p(16)
      p(25) = p(2)*p(7) - p(16)
      p(26) = p(1)*p(8)
      p(27) = p(3)*p(11)
      p(28) = p(1)*p(12)
      p(29) = p(3)*p(9)
      p(30) = p(1)*p(10)
      p(31) = p(2)*p(11) - p(25) - p(24) - p(23)
      p(32) = p(3)*p(12)
      p(33) = p(1)*p(16)
      p(34) = rm(17)
      p(35) = p(1)*p(18)
      p(36) = p(1)*p(19)
      p(37) = p(1)*p(20)
      p(38) = p(3)*p(16)
      p(39) = p(1)*p(13)
      p(40) = p(1)*p(14)
      p(41) = p(1)*p(15)
      p(42) = p(1)*p(23)
      p(43) = p(1)*p(24)
      p(44) = p(5)*p(6)
      p(45) = p(1)*p(25)
      p(46) = p(5)*p(7)
      p(47) = p(6)*p(7)
      p(48) = p(1)*p(17)
      p(49) = p(1)*p(27)
      p(50) = p(3)*p(23)
      p(51) = p(3)*p(24)
      p(52) = p(3)*p(25)
      p(53) = p(1)*p(29)
      p(54) = p(3)*p(18)
      p(55) = p(3)*p(19)
      p(56) = p(3)*p(20)
      p(57) = p(1)*p(21)
      p(58) = p(1)*p(22)
      p(59) = p(5)*p(5) - p(34) - p(34)
      p(60) = p(6)*p(6) - p(34) - p(34)
      p(61) = p(7)*p(7) - p(34) - p(34)
      p(62) = p(1)*p(31)
      p(63) = p(5)*p(11) - p(47)
      p(64) = p(6)*p(11) - p(46)
      p(65) = p(7)*p(11) - p(44)
      p(66) = p(1)*p(26)
      p(67) = p(3)*p(31)
      p(68) = p(1)*p(28)
      p(69) = p(3)*p(27)
      p(70) = p(1)*p(32)
      p(71) = p(3)*p(29)
      p(72) = p(1)*p(30)
      p(73) = p(2)*p(31) - p(65) - p(64) - p(63)
      p(74) = p(3)*p(32)
      p(75) = p(1)*p(34)
      p(76) = p(1)*p(38)
      p(77) = p(3)*p(34)
      p(78) = p(1)*p(33)
      p(79) = p(1)*p(44)
      p(80) = p(1)*p(46)
      p(81) = p(1)*p(47)
      p(82) = p(34)*p(2)
      p(83) = p(1)*p(35)
      p(84) = p(1)*p(36)
      p(85) = p(1)*p(37)
      p(86) = p(1)*p(50)
      p(87) = p(1)*p(51)
      p(88) = p(3)*p(44)
      p(89) = p(1)*p(52)
      p(90) = p(3)*p(46)
      p(91) = p(3)*p(47)
      p(92) = p(1)*p(54)
      p(93) = p(1)*p(55)
      p(94) = p(1)*p(56)
      p(95) = p(3)*p(38)
      p(96) = p(1)*p(39)
      p(97) = p(1)*p(40)
      p(98) = p(1)*p(41)
      p(99) = p(1)*p(42)
      p(100) = p(1)*p(59)
      p(101) = p(1)*p(43)
      p(102) = p(1)*p(60)
      p(103) = p(1)*p(45)
      p(104) = p(5)*p(16) - p(82)
      p(105) = p(6)*p(16) - p(82)
      p(106) = p(1)*p(61)
      p(107) = p(7)*p(16) - p(82)
      p(108) = p(1)*p(63)
      p(109) = p(1)*p(64)
      p(110) = p(5)*p(24) - p(105)
      p(111) = p(1)*p(65)
      p(112) = p(5)*p(25) - p(107)
      p(113) = p(6)*p(25) - p(107)
      p(114) = p(1)*p(48)
      p(115) = p(1)*p(49)
      p(116) = p(3)*p(59)
      p(117) = p(3)*p(60)
      p(118) = p(3)*p(61)
      p(119) = p(1)*p(67)
      p(120) = p(3)*p(63)
      p(121) = p(3)*p(64)
      p(122) = p(3)*p(65)
      p(123) = p(1)*p(53)
      p(124) = p(1)*p(69)
      p(125) = p(3)*p(50)
      p(126) = p(3)*p(51)
      p(127) = p(3)*p(52)
      p(128) = p(1)*p(71)
      p(129) = p(3)*p(54)
      p(130) = p(3)*p(55)
      p(131) = p(3)*p(56)
      p(132) = p(1)*p(57)
      p(133) = p(1)*p(58)
      p(134) = p(1)*p(62)
      p(135) = p(2)*p(59) - p(104)
      p(136) = p(2)*p(60) - p(105)
      p(137) = p(2)*p(61) - p(107)
      p(138) = p(1)*p(73)
      p(139) = p(5)*p(31) - p(113)
      p(140) = p(6)*p(31) - p(112)
      p(141) = p(7)*p(31) - p(110)
      p(142) = p(1)*p(66)
      p(143) = p(3)*p(73)
      p(144) = p(1)*p(68)
      p(145) = p(3)*p(67)
      p(146) = p(1)*p(70)
      p(147) = p(3)*p(69)
      p(148) = p(1)*p(74)
      p(149) = p(3)*p(71)
      p(150) = p(1)*p(72)
      p(151) = p(2)*p(73) - p(141) - p(140) - p(139)
      p(152) = p(3)*p(74)
      p(153) = p(1)*p(77)
      p(154) = p(1)*p(75)
      p(155) = p(1)*p(82)
      p(156) = p(1)*p(76)
      p(157) = p(1)*p(88)
      p(158) = p(1)*p(90)
      p(159) = p(1)*p(91)
      p(160) = p(3)*p(82)
      p(161) = p(1)*p(95)
      p(162) = p(3)*p(77)
      p(163) = p(1)*p(78)
      p(164) = p(1)*p(79)
      p(165) = p(1)*p(80)
      p(166) = p(1)*p(104)
      p(167) = p(1)*p(81)
      p(168) = p(34)*p(5)
      p(169) = p(1)*p(105)
      p(170) = p(34)*p(6)
      p(171) = p(1)*p(107)
      p(172) = p(34)*p(7)
      p(173) = p(1)*p(110)
      p(174) = p(1)*p(112)
      p(175) = p(1)*p(113)
      p(176) = p(34)*p(11)
      p(177) = p(1)*p(83)
      p(178) = p(1)*p(84)
      p(179) = p(1)*p(85)
      p(180) = p(1)*p(86)
      p(181) = p(1)*p(116)
      p(182) = p(1)*p(87)
      p(183) = p(1)*p(117)
      p(184) = p(1)*p(89)
      p(185) = p(3)*p(104)
      p(186) = p(3)*p(105)
      p(187) = p(1)*p(118)
      p(188) = p(3)*p(107)
      p(189) = p(1)*p(120)
      p(190) = p(1)*p(121)
      p(191) = p(3)*p(110)
      p(192) = p(1)*p(122)
      p(193) = p(3)*p(112)
      p(194) = p(3)*p(113)
      p(195) = p(1)*p(92)
      p(196) = p(1)*p(93)
      p(197) = p(1)*p(94)
      p(198) = p(1)*p(125)
      p(199) = p(1)*p(126)
      p(200) = p(3)*p(88)
      p(201) = p(1)*p(127)
      p(202) = p(3)*p(90)
      p(203) = p(3)*p(91)
      p(204) = p(1)*p(129)
      p(205) = p(1)*p(130)
      p(206) = p(1)*p(131)
      p(207) = p(3)*p(95)
      p(208) = p(1)*p(96)
      p(209) = p(1)*p(97)
      p(210) = p(1)*p(98)
      p(211) = p(1)*p(99)
      p(212) = p(1)*p(100)
      p(213) = p(1)*p(101)
      p(214) = p(1)*p(102)
      p(215) = p(1)*p(103)
      p(216) = p(1)*p(106)
      p(217) = p(5)*p(47) - p(176)
      p(218) = p(1)*p(108)
      p(219) = p(1)*p(135)
      p(220) = p(1)*p(109)
      p(221) = p(6)*p(59)
      p(222) = p(1)*p(136)
      p(223) = p(5)*p(60)
      p(224) = p(1)*p(111)
      p(225) = p(7)*p(59)
      p(226) = p(7)*p(60)
      p(227) = p(1)*p(137)
      p(228) = p(5)*p(61)
      p(229) = p(6)*p(61)
      p(230) = p(1)*p(139)
      p(231) = p(1)*p(140)
      p(232) = p(5)*p(64) - p(226)
      p(233) = p(1)*p(141)
      p(234) = p(5)*p(65) - p(229)
      p(235) = p(6)*p(65) - p(228)
      p(236) = p(1)*p(114)
      p(237) = p(1)*p(115)
      p(238) = p(1)*p(119)
      p(239) = p(3)*p(135)
      p(240) = p(3)*p(136)
      p(241) = p(3)*p(137)
      p(242) = p(1)*p(143)
      p(243) = p(3)*p(139)
      p(244) = p(3)*p(140)
      p(245) = p(3)*p(141)
      p(246) = p(1)*p(123)
      p(247) = p(1)*p(124)
      p(248) = p(3)*p(116)
      p(249) = p(3)*p(117)
      p(250) = p(3)*p(118)
      p(251) = p(1)*p(145)
      p(252) = p(3)*p(120)
      p(253) = p(3)*p(121)
      p(254) = p(3)*p(122)
      p(255) = p(1)*p(128)
      p(256) = p(1)*p(147)
      p(257) = p(3)*p(125)
      p(258) = p(3)*p(126)
      p(259) = p(3)*p(127)
      p(260) = p(1)*p(149)
      p(261) = p(3)*p(129)
      p(262) = p(3)*p(130)
      p(263) = p(3)*p(131)
      p(264) = p(1)*p(132)
      p(265) = p(1)*p(133)
      p(266) = p(1)*p(134)
      p(267) = p(5)*p(59) - p(168)
      p(268) = p(6)*p(60) - p(170)
      p(269) = p(7)*p(61) - p(172)
      p(270) = p(1)*p(138)
      p(271) = p(5)*p(63) - p(176)
      p(272) = p(6)*p(64) - p(176)
      p(273) = p(7)*p(65) - p(176)
      p(274) = p(1)*p(151)
      p(275) = p(5)*p(73) - p(235)
      p(276) = p(6)*p(73) - p(234)
      p(277) = p(7)*p(73) - p(232)
      p(278) = p(1)*p(142)
      p(279) = p(3)*p(151)
      p(280) = p(1)*p(144)
      p(281) = p(3)*p(143)
      p(282) = p(1)*p(146)
      p(283) = p(3)*p(145)
      p(284) = p(1)*p(148)
      p(285) = p(3)*p(147)
      p(286) = p(1)*p(152)
      p(287) = p(3)*p(149)
      p(288) = p(1)*p(150)
      p(289) = p(2)*p(151) - p(277) - p(276) - p(275)
      p(290) = p(3)*p(152)
      p(291) = p(1)*p(153)
      p(292) = p(1)*p(160)
      p(293) = p(1)*p(162)
      p(294) = p(1)*p(154)
      p(295) = p(1)*p(155)
      p(296) = p(1)*p(168)
      p(297) = p(1)*p(170)
      p(298) = p(1)*p(172)
      p(299) = p(1)*p(176)
      p(300) = p(1)*p(156)
      p(301) = p(1)*p(157)
      p(302) = p(1)*p(158)
      p(303) = p(1)*p(185)
      p(304) = p(1)*p(159)
      p(305) = p(3)*p(168)
      p(306) = p(1)*p(186)
      p(307) = p(3)*p(170)
      p(308) = p(1)*p(188)
      p(309) = p(3)*p(172)
      p(310) = p(1)*p(191)
      p(311) = p(1)*p(193)
      p(312) = p(1)*p(194)
      p(313) = p(3)*p(176)
      p(314) = p(1)*p(161)
      p(315) = p(1)*p(200)
      p(316) = p(1)*p(202)
      p(317) = p(1)*p(203)
      p(318) = p(3)*p(160)
      p(319) = p(1)*p(207)
      p(320) = p(3)*p(162)
      p(321) = p(1)*p(163)
      p(322) = p(1)*p(164)
      p(323) = p(1)*p(165)
      p(324) = p(1)*p(166)
      p(325) = p(1)*p(167)
      p(326) = p(1)*p(169)
      p(327) = p(1)*p(171)
      p(328) = p(1)*p(217)
      p(329) = p(34)*p(16)
      p(330) = p(1)*p(173)
      p(331) = p(1)*p(221)
      p(332) = p(1)*p(223)
      p(333) = p(1)*p(174)
      p(334) = p(1)*p(225)
      p(335) = p(1)*p(175)
      p(336) = p(34)*p(23)
      p(337) = p(1)*p(226)
      p(338) = p(34)*p(24)
      p(339) = p(1)*p(228)
      p(340) = p(1)*p(229)
      p(341) = p(34)*p(25)
      p(342) = p(1)*p(232)
      p(343) = p(1)*p(234)
      p(344) = p(1)*p(235)
      p(345) = p(34)*p(31)
      p(346) = p(1)*p(177)
      p(347) = p(1)*p(178)
      p(348) = p(1)*p(179)
      p(349) = p(1)*p(180)
      p(350) = p(1)*p(181)
      p(351) = p(1)*p(182)
      p(352) = p(1)*p(183)
      p(353) = p(1)*p(184)
      p(354) = p(1)*p(187)
      p(355) = p(3)*p(217)
      p(356) = p(1)*p(189)
      p(357) = p(1)*p(239)
      p(358) = p(1)*p(190)
      p(359) = p(3)*p(221)
      p(360) = p(1)*p(240)
      p(361) = p(3)*p(223)
      p(362) = p(1)*p(192)
      p(363) = p(3)*p(225)
      p(364) = p(3)*p(226)
      p(365) = p(1)*p(241)
      p(366) = p(3)*p(228)
      p(367) = p(3)*p(229)
      p(368) = p(1)*p(243)
      p(369) = p(1)*p(244)
      p(370) = p(3)*p(232)
      p(371) = p(1)*p(245)
      p(372) = p(3)*p(234)
      p(373) = p(3)*p(235)
      p(374) = p(1)*p(195)
      p(375) = p(1)*p(196)
      p(376) = p(1)*p(197)
      p(377) = p(1)*p(198)
      p(378) = p(1)*p(248)
      p(379) = p(1)*p(199)
      p(380) = p(1)*p(249)
      p(381) = p(1)*p(201)
      p(382) = p(3)*p(185)
      p(383) = p(3)*p(186)
      p(384) = p(1)*p(250)
      p(385) = p(3)*p(188)
      p(386) = p(1)*p(252)
      p(387) = p(1)*p(253)
      p(388) = p(3)*p(191)
      p(389) = p(1)*p(254)
      p(390) = p(3)*p(193)
      p(391) = p(3)*p(194)
      p(392) = p(1)*p(204)
      p(393) = p(1)*p(205)
      p(394) = p(1)*p(206)
      p(395) = p(1)*p(257)
      p(396) = p(1)*p(258)
      p(397) = p(3)*p(200)
      p(398) = p(1)*p(259)
      p(399) = p(3)*p(202)
      p(400) = p(3)*p(203)
      p(401) = p(1)*p(261)
      p(402) = p(1)*p(262)
      p(403) = p(1)*p(263)
      p(404) = p(3)*p(207)
      p(405) = p(1)*p(208)
      p(406) = p(1)*p(209)
      p(407) = p(1)*p(210)
      p(408) = p(1)*p(211)
      p(409) = p(1)*p(212)
      p(410) = p(1)*p(213)
      p(411) = p(1)*p(214)
      p(412) = p(1)*p(215)
      p(413) = p(1)*p(216)
      p(414) = p(1)*p(218)
      p(415) = p(1)*p(219)
      p(416) = p(1)*p(267)
      p(417) = p(1)*p(220)
      p(418) = p(1)*p(222)
      p(419) = p(5)*p(105) - p(338)
      p(420) = p(1)*p(268)
      p(421) = p(1)*p(224)
      p(422) = p(5)*p(104) - p(329)
      p(423) = p(6)*p(105) - p(329)
      p(424) = p(1)*p(227)
      p(425) = p(5)*p(107) - p(341)
      p(426) = p(5)*p(113) - p(345)
      p(427) = p(1)*p(269)
      p(428) = p(7)*p(107) - p(329)
      p(429) = p(1)*p(230)
      p(430) = p(1)*p(271)
      p(431) = p(1)*p(231)
      p(432) = p(5)*p(110) - p(338)
      p(433) = p(1)*p(272)
      p(434) = p(5)*p(136) - p(423)
      p(435) = p(1)*p(233)
      p(436) = p(5)*p(112) - p(341)
      p(437) = p(6)*p(113) - p(341)
      p(438) = p(1)*p(273)
      p(439) = p(5)*p(137) - p(428)
      p(440) = p(6)*p(137) - p(428)
      p(441) = p(1)*p(275)
      p(442) = p(1)*p(276)
      p(443) = p(5)*p(140) - p(437)
      p(444) = p(1)*p(277)
      p(445) = p(5)*p(141) - p(440)
      p(446) = p(6)*p(141) - p(439)
      p(447) = p(1)*p(236)
      p(448) = p(1)*p(237)
      p(449) = p(1)*p(238)
      p(450) = p(3)*p(267)
      p(451) = p(3)*p(268)
      p(452) = p(3)*p(269)
      p(453) = p(1)*p(242)
      p(454) = p(3)*p(271)
      p(455) = p(3)*p(272)
      p(456) = p(3)*p(273)
      p(457) = p(1)*p(279)
      p(458) = p(3)*p(275)
      p(459) = p(3)*p(276)
      p(460) = p(3)*p(277)
      p(461) = p(1)*p(246)
      p(462) = p(1)*p(247)
      p(463) = p(1)*p(251)
      p(464) = p(3)*p(239)
      p(465) = p(3)*p(240)
      p(466) = p(3)*p(241)
      p(467) = p(1)*p(281)
      p(468) = p(3)*p(243)
      p(469) = p(3)*p(244)
      p(470) = p(3)*p(245)
      p(471) = p(1)*p(255)
      p(472) = p(1)*p(256)
      p(473) = p(3)*p(248)
      p(474) = p(3)*p(249)
      p(475) = p(3)*p(250)
      p(476) = p(1)*p(283)
      p(477) = p(3)*p(252)
      p(478) = p(3)*p(253)
      p(479) = p(3)*p(254)
      p(480) = p(1)*p(260)
      p(481) = p(1)*p(285)
      p(482) = p(3)*p(257)
      p(483) = p(3)*p(258)
      p(484) = p(3)*p(259)
      p(485) = p(1)*p(287)
      p(486) = p(3)*p(261)
      p(487) = p(3)*p(262)
      p(488) = p(3)*p(263)
      p(489) = p(1)*p(264)
      p(490) = p(1)*p(265)
      p(491) = p(1)*p(266)
      p(492) = p(1)*p(270)
      p(493) = p(2)*p(267) - p(422)
      p(494) = p(2)*p(268) - p(423)
      p(495) = p(2)*p(269) - p(428)
      p(496) = p(1)*p(274)
      p(497) = p(5)*p(139) - p(345)
      p(498) = p(6)*p(140) - p(345)
      p(499) = p(7)*p(141) - p(345)
      p(500) = p(1)*p(289)
      p(501) = p(5)*p(151) - p(446)
      p(502) = p(6)*p(151) - p(445)
      p(503) = p(7)*p(151) - p(443)
      p(504) = p(1)*p(278)
      p(505) = p(3)*p(289)
      p(506) = p(1)*p(280)
      p(507) = p(3)*p(279)
      p(508) = p(1)*p(282)
      p(509) = p(3)*p(281)
      p(510) = p(1)*p(284)
      p(511) = p(3)*p(283)
      p(512) = p(1)*p(286)
      p(513) = p(3)*p(285)
      p(514) = p(1)*p(290)
      p(515) = p(3)*p(287)
      p(516) = p(1)*p(288)
      p(517) = p(2)*p(289) - p(503) - p(502) - p(501)
      p(518) = p(3)*p(290)
      p(519) = p(1)*p(291)
      p(520) = p(1)*p(292)
      p(521) = p(1)*p(305)
      p(522) = p(1)*p(307)
      p(523) = p(1)*p(309)
      p(524) = p(1)*p(313)
      p(525) = p(1)*p(293)
      p(526) = p(1)*p(318)
      p(527) = p(1)*p(320)
      p(528) = p(1)*p(294)
      p(529) = p(1)*p(295)
      p(530) = p(1)*p(296)
      p(531) = p(1)*p(297)
      p(532) = p(1)*p(298)
      p(533) = p(1)*p(329)
      p(534) = p(1)*p(299)
      p(535) = p(1)*p(336)
      p(536) = p(1)*p(338)
      p(537) = p(1)*p(341)
      p(538) = p(1)*p(345)
      p(539) = p(1)*p(300)
      p(540) = p(1)*p(301)
      p(541) = p(1)*p(302)
      p(542) = p(1)*p(303)
      p(543) = p(1)*p(304)
      p(544) = p(1)*p(306)
      p(545) = p(1)*p(308)
      p(546) = p(1)*p(355)
      p(547) = p(3)*p(329)
      p(548) = p(1)*p(310)
      p(549) = p(1)*p(359)
      p(550) = p(1)*p(361)
      p(551) = p(1)*p(311)
      p(552) = p(1)*p(363)
      p(553) = p(1)*p(312)
      p(554) = p(3)*p(336)
      p(555) = p(1)*p(364)
      p(556) = p(3)*p(338)
      p(557) = p(1)*p(366)
      p(558) = p(1)*p(367)
      p(559) = p(3)*p(341)
      p(560) = p(1)*p(370)
      p(561) = p(1)*p(372)
      p(562) = p(1)*p(373)
      p(563) = p(3)*p(345)
      p(564) = p(1)*p(314)
      p(565) = p(1)*p(315)
      p(566) = p(1)*p(316)
      p(567) = p(1)*p(382)
      p(568) = p(1)*p(317)
      p(569) = p(3)*p(305)
      p(570) = p(1)*p(383)
      p(571) = p(3)*p(307)
      p(572) = p(1)*p(385)
      p(573) = p(3)*p(309)
      p(574) = p(1)*p(388)
      p(575) = p(1)*p(390)
      p(576) = p(1)*p(391)
      p(577) = p(3)*p(313)
      p(578) = p(1)*p(319)
      p(579) = p(1)*p(397)
      p(580) = p(1)*p(399)
      p(581) = p(1)*p(400)
      p(582) = p(3)*p(318)
      p(583) = p(1)*p(404)
      p(584) = p(3)*p(320)
      p(585) = p(1)*p(321)
      p(586) = p(1)*p(322)
      p(587) = p(1)*p(323)
      p(588) = p(1)*p(324)
      p(589) = p(1)*p(325)
      p(590) = p(1)*p(326)
      p(591) = p(1)*p(327)
      p(592) = p(1)*p(328)
      p(593) = p(34)*p(34)
      p(594) = p(1)*p(330)
      p(595) = p(1)*p(331)
      p(596) = p(1)*p(332)
      p(597) = p(1)*p(419)
      p(598) = p(1)*p(333)
      p(599) = p(1)*p(334)
      p(600) = p(1)*p(422)
      p(601) = p(1)*p(335)
      p(602) = p(34)*p(59)
      p(603) = p(1)*p(337)
      p(604) = p(34)*p(44)
      p(605) = p(1)*p(423)
      p(606) = p(34)*p(60)
      p(607) = p(1)*p(339)
      p(608) = p(1)*p(425)
      p(609) = p(1)*p(340)
      p(610) = p(34)*p(46)
      p(611) = p(1)*p(426)
      p(612) = p(34)*p(47)
      p(613) = p(1)*p(428)
      p(614) = p(34)*p(61)
      p(615) = p(1)*p(342)
      p(616) = p(1)*p(432)
      p(617) = p(1)*p(434)
      p(618) = p(1)*p(343)
      p(619) = p(1)*p(436)
      p(620) = p(1)*p(344)
      p(621) = p(34)*p(63)
      p(622) = p(1)*p(437)
      p(623) = p(34)*p(64)
      p(624) = p(1)*p(439)
      p(625) = p(1)*p(440)
      p(626) = p(34)*p(65)
      p(627) = p(1)*p(443)
      p(628) = p(1)*p(445)
      p(629) = p(1)*p(446)
      p(630) = p(34)*p(73)
      p(631) = p(1)*p(346)
      p(632) = p(1)*p(347)
      p(633) = p(1)*p(348)
      p(634) = p(1)*p(349)
      p(635) = p(1)*p(350)
      p(636) = p(1)*p(351)
      p(637) = p(1)*p(352)
      p(638) = p(1)*p(353)
      p(639) = p(1)*p(354)
      p(640) = p(1)*p(356)
      p(641) = p(1)*p(357)
      p(642) = p(1)*p(450)
      p(643) = p(1)*p(358)
      p(644) = p(1)*p(360)
      p(645) = p(3)*p(419)
      p(646) = p(1)*p(451)
      p(647) = p(1)*p(362)
      p(648) = p(3)*p(422)
      p(649) = p(3)*p(423)
      p(650) = p(1)*p(365)
      p(651) = p(3)*p(425)
      p(652) = p(3)*p(426)
      p(653) = p(1)*p(452)
      p(654) = p(3)*p(428)
      p(655) = p(1)*p(368)
      p(656) = p(1)*p(454)
      p(657) = p(1)*p(369)
      p(658) = p(3)*p(432)
      p(659) = p(1)*p(455)
      p(660) = p(3)*p(434)
      p(661) = p(1)*p(371)
      p(662) = p(3)*p(436)
      p(663) = p(3)*p(437)
      p(664) = p(1)*p(456)
      p(665) = p(3)*p(439)
      p(666) = p(3)*p(440)
      p(667) = p(1)*p(458)
      p(668) = p(1)*p(459)
      p(669) = p(3)*p(443)
      p(670) = p(1)*p(460)
      p(671) = p(3)*p(445)
      p(672) = p(3)*p(446)
      p(673) = p(1)*p(374)
      p(674) = p(1)*p(375)
      p(675) = p(1)*p(376)
      p(676) = p(1)*p(377)
      p(677) = p(1)*p(378)
      p(678) = p(1)*p(379)
      p(679) = p(1)*p(380)
      p(680) = p(1)*p(381)
      p(681) = p(1)*p(384)
      p(682) = p(3)*p(355)
      p(683) = p(1)*p(386)
      p(684) = p(1)*p(464)
      p(685) = p(1)*p(387)
      p(686) = p(3)*p(359)
      p(687) = p(1)*p(465)
      p(688) = p(3)*p(361)
      p(689) = p(1)*p(389)
      p(690) = p(3)*p(363)
      p(691) = p(3)*p(364)
      p(692) = p(1)*p(466)
      p(693) = p(3)*p(366)
      p(694) = p(3)*p(367)
      p(695) = p(1)*p(468)
      p(696) = p(1)*p(469)
      p(697) = p(3)*p(370)
      p(698) = p(1)*p(470)
      p(699) = p(3)*p(372)
      p(700) = p(3)*p(373)
      p(701) = p(1)*p(392)
      p(702) = p(1)*p(393)
      p(703) = p(1)*p(394)
      p(704) = p(1)*p(395)
      p(705) = p(1)*p(473)
      p(706) = p(1)*p(396)
      p(707) = p(1)*p(474)
      p(708) = p(1)*p(398)
      p(709) = p(3)*p(382)
      p(710) = p(3)*p(383)
      p(711) = p(1)*p(475)
      p(712) = p(3)*p(385)
      p(713) = p(1)*p(477)
      p(714) = p(1)*p(478)
      p(715) = p(3)*p(388)
      p(716) = p(1)*p(479)
      p(717) = p(3)*p(390)
      p(718) = p(3)*p(391)
      p(719) = p(1)*p(401)
      p(720) = p(1)*p(402)
      p(721) = p(1)*p(403)
      p(722) = p(1)*p(482)
      p(723) = p(1)*p(483)
      p(724) = p(3)*p(397)
      p(725) = p(1)*p(484)
      p(726) = p(3)*p(399)
      p(727) = p(3)*p(400)
      p(728) = p(1)*p(486)
      p(729) = p(1)*p(487)
      p(730) = p(1)*p(488)
      p(731) = p(3)*p(404)
      p(732) = p(1)*p(405)
      p(733) = p(1)*p(406)
      p(734) = p(1)*p(407)
      p(735) = p(1)*p(408)
      p(736) = p(1)*p(409)
      p(737) = p(1)*p(410)
      p(738) = p(1)*p(411)
      p(739) = p(1)*p(412)
      p(740) = p(1)*p(413)
      p(741) = p(1)*p(414)
      p(742) = p(1)*p(415)
      p(743) = p(1)*p(416)
      p(744) = p(1)*p(417)
      p(745) = p(1)*p(418)
      p(746) = p(1)*p(420)
      p(747) = p(1)*p(421)
      p(748) = p(1)*p(424)
      p(749) = p(5)*p(217) - p(612)
      p(750) = p(5)*p(226) - p(623)
      p(751) = p(1)*p(427)
      p(752) = p(5)*p(229) - p(626)
      p(753) = p(1)*p(429)
      p(754) = p(1)*p(430)
      p(755) = p(1)*p(493)
      p(756) = p(1)*p(431)
      p(757) = p(6)*p(267)
      p(758) = p(1)*p(433)
      p(759) = p(59)*p(60)
      p(760) = p(1)*p(494)
      p(761) = p(5)*p(268)
      p(762) = p(1)*p(435)
      p(763) = p(7)*p(267)
      p(764) = p(7)*p(268)
      p(765) = p(1)*p(438)
      p(766) = p(59)*p(61)
      p(767) = p(60)*p(61)
      p(768) = p(1)*p(495)
      p(769) = p(5)*p(269)
      p(770) = p(6)*p(269)
      p(771) = p(1)*p(441)
      p(772) = p(1)*p(497)
      p(773) = p(1)*p(442)
      p(774) = p(5)*p(232) - p(623)
      p(775) = p(1)*p(498)
      p(776) = p(5)*p(272) - p(764)
      p(777) = p(1)*p(444)
      p(778) = p(5)*p(234) - p(626)
      p(779) = p(6)*p(235) - p(626)
      p(780) = p(1)*p(499)
      p(781) = p(5)*p(273) - p(770)
      p(782) = p(6)*p(273) - p(769)
      p(783) = p(1)*p(501)
      p(784) = p(1)*p(502)
      p(785) = p(5)*p(276) - p(779)
      p(786) = p(1)*p(503)
      p(787) = p(5)*p(277) - p(782)
      p(788) = p(6)*p(277) - p(781)
      p(789) = p(1)*p(447)
      p(790) = p(1)*p(448)
      p(791) = p(1)*p(449)
      p(792) = p(1)*p(453)
      p(793) = p(3)*p(493)
      p(794) = p(3)*p(494)
      p(795) = p(3)*p(495)
      p(796) = p(1)*p(457)
      p(797) = p(3)*p(497)
      p(798) = p(3)*p(498)
      p(799) = p(3)*p(499)
      p(800) = p(1)*p(505)
      p(801) = p(3)*p(501)
      p(802) = p(3)*p(502)
      p(803) = p(3)*p(503)
      p(804) = p(1)*p(461)
      p(805) = p(1)*p(462)
      p(806) = p(1)*p(463)
      p(807) = p(3)*p(450)
      p(808) = p(3)*p(451)
      p(809) = p(3)*p(452)
      p(810) = p(1)*p(467)
      p(811) = p(3)*p(454)
      p(812) = p(3)*p(455)
      p(813) = p(3)*p(456)
      p(814) = p(1)*p(507)
      p(815) = p(3)*p(458)
      p(816) = p(3)*p(459)
      p(817) = p(3)*p(460)
      p(818) = p(1)*p(471)
      p(819) = p(1)*p(472)
      p(820) = p(1)*p(476)
      p(821) = p(3)*p(464)
      p(822) = p(3)*p(465)
      p(823) = p(3)*p(466)
      p(824) = p(1)*p(509)
      p(825) = p(3)*p(468)
      p(826) = p(3)*p(469)
      p(827) = p(3)*p(470)
      p(828) = p(1)*p(480)
      p(829) = p(1)*p(481)
      p(830) = p(3)*p(473)
      p(831) = p(3)*p(474)
      p(832) = p(3)*p(475)
      p(833) = p(1)*p(511)
      p(834) = p(3)*p(477)
      p(835) = p(3)*p(478)
      p(836) = p(3)*p(479)
      p(837) = p(1)*p(485)
      p(838) = p(1)*p(513)
      p(839) = p(3)*p(482)
      p(840) = p(3)*p(483)
      p(841) = p(3)*p(484)
      p(842) = p(1)*p(515)
      p(843) = p(3)*p(486)
      p(844) = p(3)*p(487)
      p(845) = p(3)*p(488)
      p(846) = p(1)*p(489)
      p(847) = p(1)*p(490)
      p(848) = p(1)*p(491)
      p(849) = p(1)*p(492)
      p(850) = p(5)*p(267) - p(602)
      p(851) = p(6)*p(268) - p(606)
      p(852) = p(7)*p(269) - p(614)
      p(853) = p(1)*p(496)
      p(854) = p(5)*p(271) - p(621)
      p(855) = p(6)*p(272) - p(623)
      p(856) = p(7)*p(273) - p(626)
      p(857) = p(1)*p(500)
      p(858) = p(5)*p(275) - p(630)
      p(859) = p(6)*p(276) - p(630)
      p(860) = p(7)*p(277) - p(630)
      p(861) = p(1)*p(517)
      p(862) = p(5)*p(289) - p(788)
      p(863) = p(6)*p(289) - p(787)
      p(864) = p(7)*p(289) - p(785)
      p(865) = p(1)*p(504)
      p(866) = p(3)*p(517)
      p(867) = p(1)*p(506)
      p(868) = p(3)*p(505)
      p(869) = p(1)*p(508)
      p(870) = p(3)*p(507)
      p(871) = p(1)*p(510)
      p(872) = p(3)*p(509)
      p(873) = p(1)*p(512)
      p(874) = p(3)*p(511)
      p(875) = p(1)*p(514)
      p(876) = p(3)*p(513)
      p(877) = p(1)*p(518)
      p(878) = p(3)*p(515)
      p(879) = p(1)*p(516)
      p(880) = p(2)*p(517) - p(864) - p(863) - p(862)
      p(881) = p(3)*p(518)
      p(882) = p(1)*p(519)
      p(883) = p(1)*p(520)
      p(884) = p(1)*p(521)
      p(885) = p(1)*p(522)
      p(886) = p(1)*p(523)
      p(887) = p(1)*p(547)
      p(888) = p(1)*p(524)
      p(889) = p(1)*p(554)
      p(890) = p(1)*p(556)
      p(891) = p(1)*p(559)
      p(892) = p(1)*p(563)
      p(893) = p(1)*p(525)
      p(894) = p(1)*p(526)
      p(895) = p(1)*p(569)
      p(896) = p(1)*p(571)
      p(897) = p(1)*p(573)
      p(898) = p(1)*p(577)
      p(899) = p(1)*p(527)
      p(900) = p(1)*p(582)
      p(901) = p(1)*p(584)
      p(902) = p(1)*p(528)
      p(903) = p(1)*p(529)
      p(904) = p(1)*p(530)
      p(905) = p(1)*p(531)
      p(906) = p(1)*p(532)
      p(907) = p(1)*p(533)
      p(908) = p(1)*p(593)
      p(909) = p(1)*p(534)
      p(910) = p(1)*p(535)
      p(911) = p(1)*p(602)
      p(912) = p(1)*p(536)
      p(913) = p(1)*p(604)
      p(914) = p(1)*p(606)
      p(915) = p(1)*p(537)
      p(916) = p(1)*p(610)
      p(917) = p(1)*p(612)
      p(918) = p(1)*p(614)
      p(919) = p(1)*p(538)
      p(920) = p(1)*p(621)
      p(921) = p(1)*p(623)
      p(922) = p(1)*p(626)
      p(923) = p(1)*p(630)
      p(924) = p(1)*p(539)
      p(925) = p(1)*p(540)
      p(926) = p(1)*p(541)
      p(927) = p(1)*p(542)
      p(928) = p(1)*p(543)
      p(929) = p(1)*p(544)
      p(930) = p(1)*p(545)
      p(931) = p(1)*p(546)
      p(932) = p(3)*p(593)
      p(933) = p(1)*p(548)
      p(934) = p(1)*p(549)
      p(935) = p(1)*p(550)
      p(936) = p(1)*p(645)
      p(937) = p(1)*p(551)
      p(938) = p(1)*p(552)
      p(939) = p(1)*p(648)
      p(940) = p(1)*p(553)
      p(941) = p(3)*p(602)
      p(942) = p(1)*p(555)
      p(943) = p(3)*p(604)
      p(944) = p(1)*p(649)
      p(945) = p(3)*p(606)
      p(946) = p(1)*p(557)
      p(947) = p(1)*p(651)
      p(948) = p(1)*p(558)
      p(949) = p(3)*p(610)
      p(950) = p(1)*p(652)
      p(951) = p(3)*p(612)
      p(952) = p(1)*p(654)
      p(953) = p(3)*p(614)
      p(954) = p(1)*p(560)
      p(955) = p(1)*p(658)
      p(956) = p(1)*p(660)
      p(957) = p(1)*p(561)
      p(958) = p(1)*p(662)
      p(959) = p(1)*p(562)
      p(960) = p(3)*p(621)
      p(961) = p(1)*p(663)
      p(962) = p(3)*p(623)
      p(963) = p(1)*p(665)
      p(964) = p(1)*p(666)
      p(965) = p(3)*p(626)
      p(966) = p(1)*p(669)
      p(967) = p(1)*p(671)
      p(968) = p(1)*p(672)
      p(969) = p(3)*p(630)
      p(970) = p(1)*p(564)
      p(971) = p(1)*p(565)
      p(972) = p(1)*p(566)
      p(973) = p(1)*p(567)
      p(974) = p(1)*p(568)
      p(975) = p(1)*p(570)
      p(976) = p(1)*p(572)
      p(977) = p(1)*p(682)
      p(978) = p(3)*p(547)
      p(979) = p(1)*p(574)
      p(980) = p(1)*p(686)
      p(981) = p(1)*p(688)
      p(982) = p(1)*p(575)
      p(983) = p(1)*p(690)
      p(984) = p(1)*p(576)
      p(985) = p(3)*p(554)
      p(986) = p(1)*p(691)
      p(987) = p(3)*p(556)
      p(988) = p(1)*p(693)
      p(989) = p(1)*p(694)
      p(990) = p(3)*p(559)
      p(991) = p(1)*p(697)
      p(992) = p(1)*p(699)
      p(993) = p(1)*p(700)
      p(994) = p(3)*p(563)
      p(995) = p(1)*p(578)
      p(996) = p(1)*p(579)
      p(997) = p(1)*p(580)
      p(998) = p(1)*p(709)
      p(999) = p(1)*p(581)
      p(1000) = p(3)*p(569)
      p(1001) = p(1)*p(710)
      p(1002) = p(3)*p(571)
      p(1003) = p(1)*p(712)
      p(1004) = p(3)*p(573)
      p(1005) = p(1)*p(715)
      p(1006) = p(1)*p(717)
      p(1007) = p(1)*p(718)
      p(1008) = p(3)*p(577)
      p(1009) = p(1)*p(583)
      p(1010) = p(1)*p(724)
      p(1011) = p(1)*p(726)
      p(1012) = p(1)*p(727)
      p(1013) = p(3)*p(582)
      p(1014) = p(1)*p(731)
      p(1015) = p(3)*p(584)
      p(1016) = p(1)*p(585)
      p(1017) = p(1)*p(586)
      p(1018) = p(1)*p(587)
      p(1019) = p(1)*p(588)
      p(1020) = p(1)*p(589)
      p(1021) = p(1)*p(590)
      p(1022) = p(1)*p(591)
      p(1023) = p(1)*p(592)
      p(1024) = p(1)*p(594)
      p(1025) = p(1)*p(595)
      p(1026) = p(1)*p(596)
      p(1027) = p(1)*p(597)
      p(1028) = p(1)*p(598)
      p(1029) = p(1)*p(599)
      p(1030) = p(1)*p(600)
      p(1031) = p(1)*p(601)
      p(1032) = p(1)*p(603)
      p(1033) = p(1)*p(605)
      p(1034) = p(1)*p(607)
      p(1035) = p(1)*p(608)
      p(1036) = p(1)*p(749)
      p(1037) = p(1)*p(609)
      p(1038) = p(34)*p(104)
      p(1039) = p(1)*p(611)
      p(1040) = p(34)*p(82)
      p(1041) = p(1)*p(750)
      p(1042) = p(34)*p(105)
      p(1043) = p(1)*p(613)
      p(1044) = p(1)*p(752)
      p(1045) = p(34)*p(107)
      p(1046) = p(1)*p(615)
      p(1047) = p(1)*p(616)
      p(1048) = p(1)*p(757)
      p(1049) = p(1)*p(617)
      p(1050) = p(1)*p(759)
      p(1051) = p(1)*p(761)
      p(1052) = p(1)*p(618)
      p(1053) = p(1)*p(619)
      p(1054) = p(1)*p(763)
      p(1055) = p(1)*p(620)
      p(1056) = p(34)*p(135)
      p(1057) = p(1)*p(622)
      p(1058) = p(34)*p(110)
      p(1059) = p(1)*p(764)
      p(1060) = p(34)*p(136)
      p(1061) = p(1)*p(624)
      p(1062) = p(1)*p(766)
      p(1063) = p(1)*p(625)
      p(1064) = p(34)*p(112)
      p(1065) = p(1)*p(767)
      p(1066) = p(34)*p(113)
      p(1067) = p(1)*p(769)
      p(1068) = p(1)*p(770)
      p(1069) = p(34)*p(137)
      p(1070) = p(1)*p(627)
      p(1071) = p(1)*p(774)
      p(1072) = p(1)*p(776)
      p(1073) = p(1)*p(628)
      p(1074) = p(1)*p(778)
      p(1075) = p(1)*p(629)
      p(1076) = p(34)*p(139)
      p(1077) = p(1)*p(779)
      p(1078) = p(34)*p(140)
      p(1079) = p(1)*p(781)
      p(1080) = p(1)*p(782)
      p(1081) = p(34)*p(141)
      p(1082) = p(1)*p(785)
      p(1083) = p(1)*p(787)
      p(1084) = p(1)*p(788)
      p(1085) = p(34)*p(151)
      p(1086) = p(1)*p(631)
      p(1087) = p(1)*p(632)
      p(1088) = p(1)*p(633)
      p(1089) = p(1)*p(634)
      p(1090) = p(1)*p(635)
      p(1091) = p(1)*p(636)
      p(1092) = p(1)*p(637)
      p(1093) = p(1)*p(638)
      p(1094) = p(1)*p(639)
      p(1095) = p(1)*p(640)
      p(1096) = p(1)*p(641)
      p(1097) = p(1)*p(642)
      p(1098) = p(1)*p(643)
      p(1099) = p(1)*p(644)
      p(1100) = p(1)*p(646)
      p(1101) = p(1)*p(647)
      p(1102) = p(1)*p(650)
      p(1103) = p(3)*p(749)
      p(1104) = p(3)*p(750)
      p(1105) = p(1)*p(653)
      p(1106) = p(3)*p(752)
      p(1107) = p(1)*p(655)
      p(1108) = p(1)*p(656)
      p(1109) = p(1)*p(793)
      p(1110) = p(1)*p(657)
      p(1111) = p(3)*p(757)
      p(1112) = p(1)*p(659)
      p(1113) = p(3)*p(759)
      p(1114) = p(1)*p(794)
      p(1115) = p(3)*p(761)
      p(1116) = p(1)*p(661)
      p(1117) = p(3)*p(763)
      p(1118) = p(3)*p(764)
      p(1119) = p(1)*p(664)
      p(1120) = p(3)*p(766)
      p(1121) = p(3)*p(767)
      p(1122) = p(1)*p(795)
      p(1123) = p(3)*p(769)
      p(1124) = p(3)*p(770)
      p(1125) = p(1)*p(667)
      p(1126) = p(1)*p(797)
      p(1127) = p(1)*p(668)
      p(1128) = p(3)*p(774)
      p(1129) = p(1)*p(798)
      p(1130) = p(3)*p(776)
      p(1131) = p(1)*p(670)
      p(1132) = p(3)*p(778)
      p(1133) = p(3)*p(779)
      p(1134) = p(1)*p(799)
      p(1135) = p(3)*p(781)
      p(1136) = p(3)*p(782)
      p(1137) = p(1)*p(801)
      p(1138) = p(1)*p(802)
      p(1139) = p(3)*p(785)
      p(1140) = p(1)*p(803)
      p(1141) = p(3)*p(787)
      p(1142) = p(3)*p(788)
      p(1143) = p(1)*p(673)
      p(1144) = p(1)*p(674)
      p(1145) = p(1)*p(675)
      p(1146) = p(1)*p(676)
      p(1147) = p(1)*p(677)
      p(1148) = p(1)*p(678)
      p(1149) = p(1)*p(679)
      p(1150) = p(1)*p(680)
      p(1151) = p(1)*p(681)
      p(1152) = p(1)*p(683)
      p(1153) = p(1)*p(684)
      p(1154) = p(1)*p(807)
      p(1155) = p(1)*p(685)
      p(1156) = p(1)*p(687)
      p(1157) = p(3)*p(645)
      p(1158) = p(1)*p(808)
      p(1159) = p(1)*p(689)
      p(1160) = p(3)*p(648)
      p(1161) = p(3)*p(649)
      p(1162) = p(1)*p(692)
      p(1163) = p(3)*p(651)
      p(1164) = p(3)*p(652)
      p(1165) = p(1)*p(809)
      p(1166) = p(3)*p(654)
      p(1167) = p(1)*p(695)
      p(1168) = p(1)*p(811)
      p(1169) = p(1)*p(696)
      p(1170) = p(3)*p(658)
      p(1171) = p(1)*p(812)
      p(1172) = p(3)*p(660)
      p(1173) = p(1)*p(698)
      p(1174) = p(3)*p(662)
      p(1175) = p(3)*p(663)
      p(1176) = p(1)*p(813)
      p(1177) = p(3)*p(665)
      p(1178) = p(3)*p(666)
      p(1179) = p(1)*p(815)
      p(1180) = p(1)*p(816)
      p(1181) = p(3)*p(669)
      p(1182) = p(1)*p(817)
      p(1183) = p(3)*p(671)
      p(1184) = p(3)*p(672)
      p(1185) = p(1)*p(701)
      p(1186) = p(1)*p(702)
      p(1187) = p(1)*p(703)
      p(1188) = p(1)*p(704)
      p(1189) = p(1)*p(705)
      p(1190) = p(1)*p(706)
      p(1191) = p(1)*p(707)
      p(1192) = p(1)*p(708)
      p(1193) = p(1)*p(711)
      p(1194) = p(3)*p(682)
      p(1195) = p(1)*p(713)
      p(1196) = p(1)*p(821)
      p(1197) = p(1)*p(714)
      p(1198) = p(3)*p(686)
      p(1199) = p(1)*p(822)
      p(1200) = p(3)*p(688)
      p(1201) = p(1)*p(716)
      p(1202) = p(3)*p(690)
      p(1203) = p(3)*p(691)
      p(1204) = p(1)*p(823)
      p(1205) = p(3)*p(693)
      p(1206) = p(3)*p(694)
      p(1207) = p(1)*p(825)
      p(1208) = p(1)*p(826)
      p(1209) = p(3)*p(697)
      p(1210) = p(1)*p(827)
      p(1211) = p(3)*p(699)
      p(1212) = p(3)*p(700)
      p(1213) = p(1)*p(719)
      p(1214) = p(1)*p(720)
      p(1215) = p(1)*p(721)
      p(1216) = p(1)*p(722)
      p(1217) = p(1)*p(830)
      p(1218) = p(1)*p(723)
      p(1219) = p(1)*p(831)
      p(1220) = p(1)*p(725)
      p(1221) = p(3)*p(709)
      p(1222) = p(3)*p(710)
      p(1223) = p(1)*p(832)
      p(1224) = p(3)*p(712)
      p(1225) = p(1)*p(834)
      p(1226) = p(1)*p(835)
      p(1227) = p(3)*p(715)
      p(1228) = p(1)*p(836)
      p(1229) = p(3)*p(717)
      p(1230) = p(3)*p(718)
      p(1231) = p(1)*p(728)
      p(1232) = p(1)*p(729)
      p(1233) = p(1)*p(730)
      p(1234) = p(1)*p(839)
      p(1235) = p(1)*p(840)
      p(1236) = p(3)*p(724)
      p(1237) = p(1)*p(841)
      p(1238) = p(3)*p(726)
      p(1239) = p(3)*p(727)
      p(1240) = p(1)*p(843)
      p(1241) = p(1)*p(844)
      p(1242) = p(1)*p(845)
      p(1243) = p(3)*p(731)
      p(1244) = p(1)*p(732)
      p(1245) = p(1)*p(733)
      p(1246) = p(1)*p(734)
      p(1247) = p(1)*p(735)
      p(1248) = p(1)*p(736)
      p(1249) = p(1)*p(737)
      p(1250) = p(1)*p(738)
      p(1251) = p(1)*p(739)
      p(1252) = p(1)*p(740)
      p(1253) = p(1)*p(741)
      p(1254) = p(1)*p(742)
      p(1255) = p(1)*p(743)
      p(1256) = p(1)*p(744)
      p(1257) = p(1)*p(745)
      p(1258) = p(1)*p(746)
      p(1259) = p(1)*p(747)
      p(1260) = p(1)*p(748)
      p(1261) = p(1)*p(751)
      p(1262) = p(5)*p(426) - p(1066)
      p(1263) = p(1)*p(753)
      p(1264) = p(1)*p(754)
      p(1265) = p(1)*p(755)
      p(1266) = p(1)*p(850)
      p(1267) = p(1)*p(756)
      p(1268) = p(1)*p(758)
      p(1269) = p(5)*p(419) - p(1042)
      p(1270) = p(1)*p(760)
      p(1271) = p(5)*p(423) - p(1060)
      p(1272) = p(1)*p(851)
      p(1273) = p(1)*p(762)
      p(1274) = p(5)*p(422) - p(1038)
      p(1275) = p(6)*p(423) - p(1042)
      p(1276) = p(1)*p(765)
      p(1277) = p(5)*p(425) - p(1045)
      p(1278) = p(5)*p(437) - p(1078)
      p(1279) = p(1)*p(768)
      p(1280) = p(5)*p(428) - p(1069)
      p(1281) = p(5)*p(440) - p(1081)
      p(1282) = p(1)*p(852)
      p(1283) = p(7)*p(428) - p(1045)
      p(1284) = p(1)*p(771)
      p(1285) = p(1)*p(772)
      p(1286) = p(1)*p(854)
      p(1287) = p(1)*p(773)
      p(1288) = p(5)*p(432) - p(1058)
      p(1289) = p(1)*p(775)
      p(1290) = p(5)*p(434) - p(1060)
      p(1291) = p(1)*p(855)
      p(1292) = p(5)*p(494) - p(1275)
      p(1293) = p(1)*p(777)
      p(1294) = p(5)*p(436) - p(1064)
      p(1295) = p(6)*p(437) - p(1066)
      p(1296) = p(1)*p(780)
      p(1297) = p(5)*p(439) - p(1069)
      p(1298) = p(5)*p(446) - p(1085)
      p(1299) = p(1)*p(856)
      p(1300) = p(5)*p(495) - p(1283)
      p(1301) = p(6)*p(495) - p(1283)
      p(1302) = p(1)*p(783)
      p(1303) = p(1)*p(858)
      p(1304) = p(1)*p(784)
      p(1305) = p(5)*p(443) - p(1078)
      p(1306) = p(1)*p(859)
      p(1307) = p(5)*p(498) - p(1295)
      p(1308) = p(1)*p(786)
      p(1309) = p(5)*p(445) - p(1081)
      p(1310) = p(6)*p(446) - p(1081)
      p(1311) = p(1)*p(860)
      p(1312) = p(5)*p(499) - p(1301)
      p(1313) = p(6)*p(499) - p(1300)
      p(1314) = p(1)*p(862)
      p(1315) = p(1)*p(863)
      p(1316) = p(5)*p(502) - p(1310)
      p(1317) = p(1)*p(864)
      p(1318) = p(5)*p(503) - p(1313)
      p(1319) = p(6)*p(503) - p(1312)
      p(1320) = p(1)*p(789)
      p(1321) = p(1)*p(790)
      p(1322) = p(1)*p(791)
      p(1323) = p(1)*p(792)
      p(1324) = p(3)*p(850)
      p(1325) = p(3)*p(851)
      p(1326) = p(3)*p(852)
      p(1327) = p(1)*p(796)
      p(1328) = p(3)*p(854)
      p(1329) = p(3)*p(855)
      p(1330) = p(3)*p(856)
      p(1331) = p(1)*p(800)
      p(1332) = p(3)*p(858)
      p(1333) = p(3)*p(859)
      p(1334) = p(3)*p(860)
      p(1335) = p(1)*p(866)
      p(1336) = p(3)*p(862)
      p(1337) = p(3)*p(863)
      p(1338) = p(3)*p(864)
      p(1339) = p(1)*p(804)
      p(1340) = p(1)*p(805)
      p(1341) = p(1)*p(806)
      p(1342) = p(1)*p(810)
      p(1343) = p(3)*p(793)
      p(1344) = p(3)*p(794)
      p(1345) = p(3)*p(795)
      p(1346) = p(1)*p(814)
      p(1347) = p(3)*p(797)
      p(1348) = p(3)*p(798)
      p(1349) = p(3)*p(799)
      p(1350) = p(1)*p(868)
      p(1351) = p(3)*p(801)
      p(1352) = p(3)*p(802)
      p(1353) = p(3)*p(803)
      p(1354) = p(1)*p(818)
      p(1355) = p(1)*p(819)
      p(1356) = p(1)*p(820)
      p(1357) = p(3)*p(807)
      p(1358) = p(3)*p(808)
      p(1359) = p(3)*p(809)
      p(1360) = p(1)*p(824)
      p(1361) = p(3)*p(811)
      p(1362) = p(3)*p(812)
      p(1363) = p(3)*p(813)
      p(1364) = p(1)*p(870)
      p(1365) = p(3)*p(815)
      p(1366) = p(3)*p(816)
      p(1367) = p(3)*p(817)
      p(1368) = p(1)*p(828)
      p(1369) = p(1)*p(829)
      p(1370) = p(1)*p(833)
      p(1371) = p(3)*p(821)
      p(1372) = p(3)*p(822)
      p(1373) = p(3)*p(823)
      p(1374) = p(1)*p(872)
      p(1375) = p(3)*p(825)
      p(1376) = p(3)*p(826)
      p(1377) = p(3)*p(827)
      p(1378) = p(1)*p(837)
      p(1379) = p(1)*p(838)
      p(1380) = p(3)*p(830)
      p(1381) = p(3)*p(831)
      p(1382) = p(3)*p(832)
      p(1383) = p(1)*p(874)
      p(1384) = p(3)*p(834)
      p(1385) = p(3)*p(835)
      p(1386) = p(3)*p(836)
      p(1387) = p(1)*p(842)
      p(1388) = p(1)*p(876)
      p(1389) = p(3)*p(839)
      p(1390) = p(3)*p(840)
      p(1391) = p(3)*p(841)
      p(1392) = p(1)*p(878)
      p(1393) = p(3)*p(843)
      p(1394) = p(3)*p(844)
      p(1395) = p(3)*p(845)
      p(1396) = p(1)*p(846)
      p(1397) = p(1)*p(847)
      p(1398) = p(1)*p(848)
      p(1399) = p(1)*p(849)
      p(1400) = p(1)*p(853)
      p(1401) = p(2)*p(850) - p(1274)
      p(1402) = p(2)*p(851) - p(1275)
      p(1403) = p(2)*p(852) - p(1283)
      p(1404) = p(1)*p(857)
      p(1405) = p(5)*p(497) - p(1076)
      p(1406) = p(6)*p(498) - p(1078)
      p(1407) = p(7)*p(499) - p(1081)
      p(1408) = p(1)*p(861)
      p(1409) = p(5)*p(501) - p(1085)
      p(1410) = p(6)*p(502) - p(1085)
      p(1411) = p(7)*p(503) - p(1085)
      p(1412) = p(1)*p(880)
      p(1413) = p(5)*p(517) - p(1319)
      p(1414) = p(6)*p(517) - p(1318)
      p(1415) = p(7)*p(517) - p(1316)
      p(1416) = p(1)*p(865)
      p(1417) = p(3)*p(880)
      p(1418) = p(1)*p(867)
      p(1419) = p(3)*p(866)
      p(1420) = p(1)*p(869)
      p(1421) = p(3)*p(868)
      p(1422) = p(1)*p(871)
      p(1423) = p(3)*p(870)
      p(1424) = p(1)*p(873)
      p(1425) = p(3)*p(872)
      p(1426) = p(1)*p(875)
      p(1427) = p(3)*p(874)
      p(1428) = p(1)*p(877)
      p(1429) = p(3)*p(876)
      p(1430) = p(1)*p(881)
      p(1431) = p(3)*p(878)
      p(1432) = p(1)*p(879)
      p(1433) = p(2)*p(880) - p(1415) - p(1414) - p(1413)
      p(1434) = p(3)*p(881)
      p(1435) = p(1)*p(882)
      p(1436) = p(1)*p(883)
      p(1437) = p(1)*p(884)
      p(1438) = p(1)*p(885)
      p(1439) = p(1)*p(886)
      p(1440) = p(1)*p(887)
      p(1441) = p(1)*p(932)
      p(1442) = p(1)*p(888)
      p(1443) = p(1)*p(889)
      p(1444) = p(1)*p(941)
      p(1445) = p(1)*p(890)
      p(1446) = p(1)*p(943)
      p(1447) = p(1)*p(945)
      p(1448) = p(1)*p(891)
      p(1449) = p(1)*p(949)
      p(1450) = p(1)*p(951)
      p(1451) = p(1)*p(953)
      p(1452) = p(1)*p(892)
      p(1453) = p(1)*p(960)
      p(1454) = p(1)*p(962)
      p(1455) = p(1)*p(965)
      p(1456) = p(1)*p(969)
      p(1457) = p(1)*p(893)
      p(1458) = p(1)*p(894)
      p(1459) = p(1)*p(895)
      p(1460) = p(1)*p(896)
      p(1461) = p(1)*p(897)
      p(1462) = p(1)*p(978)
      p(1463) = p(1)*p(898)
      p(1464) = p(1)*p(985)
      p(1465) = p(1)*p(987)
      p(1466) = p(1)*p(990)
      p(1467) = p(1)*p(994)
      p(1468) = p(1)*p(899)
      p(1469) = p(1)*p(900)
      p(1470) = p(1)*p(1000)
      p(1471) = p(1)*p(1002)
      p(1472) = p(1)*p(1004)
      p(1473) = p(1)*p(1008)
      p(1474) = p(1)*p(901)
      p(1475) = p(1)*p(1013)
      p(1476) = p(1)*p(1015)
      p(1477) = p(1)*p(902)
      p(1478) = p(1)*p(903)
      p(1479) = p(1)*p(904)
      p(1480) = p(1)*p(905)
      p(1481) = p(1)*p(906)
      p(1482) = p(1)*p(907)
      p(1483) = p(1)*p(908)
      p(1484) = p(1)*p(909)
      p(1485) = p(1)*p(910)
      p(1486) = p(1)*p(911)
      p(1487) = p(1)*p(912)
      p(1488) = p(1)*p(913)
      p(1489) = p(1)*p(914)
      p(1490) = p(1)*p(915)
      p(1491) = p(1)*p(916)
      p(1492) = p(1)*p(1038)
      p(1493) = p(1)*p(917)
      p(1494) = p(1)*p(1040)
      p(1495) = p(1)*p(1042)
      p(1496) = p(1)*p(918)
      p(1497) = p(1)*p(1045)
      p(1498) = p(1)*p(919)
      p(1499) = p(1)*p(920)
      p(1500) = p(1)*p(1056)
      p(1501) = p(1)*p(921)
      p(1502) = p(1)*p(1058)
      p(1503) = p(1)*p(1060)
      p(1504) = p(1)*p(922)
      p(1505) = p(1)*p(1064)
      p(1506) = p(1)*p(1066)
      p(1507) = p(1)*p(1069)
      p(1508) = p(1)*p(923)
      p(1509) = p(1)*p(1076)
      p(1510) = p(1)*p(1078)
      p(1511) = p(1)*p(1081)
      p(1512) = p(1)*p(1085)
      p(1513) = p(1)*p(924)
      p(1514) = p(1)*p(925)
      p(1515) = p(1)*p(926)
      p(1516) = p(1)*p(927)
      p(1517) = p(1)*p(928)
      p(1518) = p(1)*p(929)
      p(1519) = p(1)*p(930)
      p(1520) = p(1)*p(931)
      p(1521) = p(1)*p(933)
      p(1522) = p(1)*p(934)
      p(1523) = p(1)*p(935)
      p(1524) = p(1)*p(936)
      p(1525) = p(1)*p(937)
      p(1526) = p(1)*p(938)
      p(1527) = p(1)*p(939)
      p(1528) = p(1)*p(940)
      p(1529) = p(1)*p(942)
      p(1530) = p(1)*p(944)
      p(1531) = p(1)*p(946)
      p(1532) = p(1)*p(947)
      p(1533) = p(1)*p(1103)
      p(1534) = p(1)*p(948)
      p(1535) = p(3)*p(1038)
      p(1536) = p(1)*p(950)
      p(1537) = p(3)*p(1040)
      p(1538) = p(1)*p(1104)
      p(1539) = p(3)*p(1042)
      p(1540) = p(1)*p(952)
      p(1541) = p(1)*p(1106)
      p(1542) = p(3)*p(1045)
      p(1543) = p(1)*p(954)
      p(1544) = p(1)*p(955)
      p(1545) = p(1)*p(1111)
      p(1546) = p(1)*p(956)
      p(1547) = p(1)*p(1113)
      p(1548) = p(1)*p(1115)
      p(1549) = p(1)*p(957)
      p(1550) = p(1)*p(958)
      p(1551) = p(1)*p(1117)
      p(1552) = p(1)*p(959)
      p(1553) = p(3)*p(1056)
      p(1554) = p(1)*p(961)
      p(1555) = p(3)*p(1058)
      p(1556) = p(1)*p(1118)
      p(1557) = p(3)*p(1060)
      p(1558) = p(1)*p(963)
      p(1559) = p(1)*p(1120)
      p(1560) = p(1)*p(964)
      p(1561) = p(3)*p(1064)
      p(1562) = p(1)*p(1121)
      p(1563) = p(3)*p(1066)
      p(1564) = p(1)*p(1123)
      p(1565) = p(1)*p(1124)
      p(1566) = p(3)*p(1069)
      p(1567) = p(1)*p(966)
      p(1568) = p(1)*p(1128)
      p(1569) = p(1)*p(1130)
      p(1570) = p(1)*p(967)
      p(1571) = p(1)*p(1132)
      p(1572) = p(1)*p(968)
      p(1573) = p(3)*p(1076)
      p(1574) = p(1)*p(1133)
      p(1575) = p(3)*p(1078)
      p(1576) = p(1)*p(1135)
      p(1577) = p(1)*p(1136)
      p(1578) = p(3)*p(1081)
      p(1579) = p(1)*p(1139)
      p(1580) = p(1)*p(1141)
      p(1581) = p(1)*p(1142)
      p(1582) = p(3)*p(1085)
      p(1583) = p(1)*p(970)
      p(1584) = p(1)*p(971)
      p(1585) = p(1)*p(972)
      p(1586) = p(1)*p(973)
      p(1587) = p(1)*p(974)
      p(1588) = p(1)*p(975)
      p(1589) = p(1)*p(976)
      p(1590) = p(1)*p(977)
      p(1591) = p(3)*p(932)
      p(1592) = p(1)*p(979)
      p(1593) = p(1)*p(980)
      p(1594) = p(1)*p(981)
      p(1595) = p(1)*p(1157)
      p(1596) = p(1)*p(982)
      p(1597) = p(1)*p(983)
      p(1598) = p(1)*p(1160)
      p(1599) = p(1)*p(984)
      p(1600) = p(3)*p(941)
      p(1601) = p(1)*p(986)
      p(1602) = p(3)*p(943)
      p(1603) = p(1)*p(1161)
      p(1604) = p(3)*p(945)
      p(1605) = p(1)*p(988)
      p(1606) = p(1)*p(1163)
      p(1607) = p(1)*p(989)
      p(1608) = p(3)*p(949)
      p(1609) = p(1)*p(1164)
      p(1610) = p(3)*p(951)
      p(1611) = p(1)*p(1166)
      p(1612) = p(3)*p(953)
      p(1613) = p(1)*p(991)
      p(1614) = p(1)*p(1170)
      p(1615) = p(1)*p(1172)
      p(1616) = p(1)*p(992)
      p(1617) = p(1)*p(1174)
      p(1618) = p(1)*p(993)
      p(1619) = p(3)*p(960)
      p(1620) = p(1)*p(1175)
      p(1621) = p(3)*p(962)
      p(1622) = p(1)*p(1177)
      p(1623) = p(1)*p(1178)
      p(1624) = p(3)*p(965)
      p(1625) = p(1)*p(1181)
      p(1626) = p(1)*p(1183)
      p(1627) = p(1)*p(1184)
      p(1628) = p(3)*p(969)
      p(1629) = p(1)*p(995)
      p(1630) = p(1)*p(996)
      p(1631) = p(1)*p(997)
      p(1632) = p(1)*p(998)
      p(1633) = p(1)*p(999)
      p(1634) = p(1)*p(1001)
      p(1635) = p(1)*p(1003)
      p(1636) = p(1)*p(1194)
      p(1637) = p(3)*p(978)
      p(1638) = p(1)*p(1005)
      p(1639) = p(1)*p(1198)
      p(1640) = p(1)*p(1200)
      p(1641) = p(1)*p(1006)
      p(1642) = p(1)*p(1202)
      p(1643) = p(1)*p(1007)
      p(1644) = p(3)*p(985)
      p(1645) = p(1)*p(1203)
      p(1646) = p(3)*p(987)
      p(1647) = p(1)*p(1205)
      p(1648) = p(1)*p(1206)
      p(1649) = p(3)*p(990)
      p(1650) = p(1)*p(1209)
      p(1651) = p(1)*p(1211)
      p(1652) = p(1)*p(1212)
      p(1653) = p(3)*p(994)
      p(1654) = p(1)*p(1009)
      p(1655) = p(1)*p(1010)
      p(1656) = p(1)*p(1011)
      p(1657) = p(1)*p(1221)
      p(1658) = p(1)*p(1012)
      p(1659) = p(3)*p(1000)
      p(1660) = p(1)*p(1222)
      p(1661) = p(3)*p(1002)
      p(1662) = p(1)*p(1224)
      p(1663) = p(3)*p(1004)
      p(1664) = p(1)*p(1227)
      p(1665) = p(1)*p(1229)
      p(1666) = p(1)*p(1230)
      p(1667) = p(3)*p(1008)
      p(1668) = p(1)*p(1014)
      p(1669) = p(1)*p(1236)
      p(1670) = p(1)*p(1238)
      p(1671) = p(1)*p(1239)
      p(1672) = p(3)*p(1013)
      p(1673) = p(1)*p(1243)
      p(1674) = p(3)*p(1015)
      p(1675) = p(1)*p(1016)
      p(1676) = p(1)*p(1017)
      p(1677) = p(1)*p(1018)
      p(1678) = p(1)*p(1019)
      p(1679) = p(1)*p(1020)
      p(1680) = p(1)*p(1021)
      p(1681) = p(1)*p(1022)
      p(1682) = p(1)*p(1023)
      p(1683) = p(1)*p(1024)
      p(1684) = p(1)*p(1025)
      p(1685) = p(1)*p(1026)
      p(1686) = p(1)*p(1027)
      p(1687) = p(1)*p(1028)
      p(1688) = p(1)*p(1029)
      p(1689) = p(1)*p(1030)
      p(1690) = p(1)*p(1031)
      p(1691) = p(1)*p(1032)
      p(1692) = p(1)*p(1033)
      p(1693) = p(1)*p(1034)
      p(1694) = p(1)*p(1035)
      p(1695) = p(1)*p(1036)
      p(1696) = p(1)*p(1037)
      p(1697) = p(1)*p(1039)
      p(1698) = p(34)*p(168)
      p(1699) = p(1)*p(1041)
      p(1700) = p(34)*p(170)
      p(1701) = p(1)*p(1043)
      p(1702) = p(1)*p(1044)
      p(1703) = p(34)*p(172)
      p(1704) = p(1)*p(1262)
      p(1705) = p(34)*p(217)
      p(1706) = p(1)*p(1046)
      p(1707) = p(1)*p(1047)
      p(1708) = p(1)*p(1048)
      p(1709) = p(1)*p(1049)
      p(1710) = p(1)*p(1050)
      p(1711) = p(1)*p(1269)
      p(1712) = p(1)*p(1051)
      p(1713) = p(1)*p(1271)
      p(1714) = p(1)*p(1052)
      p(1715) = p(1)*p(1053)
      p(1716) = p(1)*p(1054)
      p(1717) = p(1)*p(1274)
      p(1718) = p(1)*p(1055)
      p(1719) = p(34)*p(267)
      p(1720) = p(1)*p(1057)
      p(1721) = p(34)*p(221)
      p(1722) = p(1)*p(1059)
      p(1723) = p(34)*p(223)
      p(1724) = p(1)*p(1275)
      p(1725) = p(34)*p(268)
      p(1726) = p(1)*p(1061)
      p(1727) = p(1)*p(1062)
      p(1728) = p(1)*p(1277)
      p(1729) = p(1)*p(1063)
      p(1730) = p(34)*p(225)
      p(1731) = p(1)*p(1065)
      p(1732) = p(34)*p(176)
      p(1733) = p(1)*p(1278)
      p(1734) = p(34)*p(226)
      p(1735) = p(1)*p(1067)
      p(1736) = p(1)*p(1280)
      p(1737) = p(1)*p(1068)
      p(1738) = p(34)*p(228)
      p(1739) = p(1)*p(1281)
      p(1740) = p(34)*p(229)
      p(1741) = p(1)*p(1283)
      p(1742) = p(34)*p(269)
      p(1743) = p(1)*p(1070)
      p(1744) = p(1)*p(1071)
      p(1745) = p(1)*p(1288)
      p(1746) = p(1)*p(1072)
      p(1747) = p(1)*p(1290)
      p(1748) = p(1)*p(1292)
      p(1749) = p(1)*p(1073)
      p(1750) = p(1)*p(1074)
      p(1751) = p(1)*p(1294)
      p(1752) = p(1)*p(1075)
      p(1753) = p(34)*p(271)
      p(1754) = p(1)*p(1077)
      p(1755) = p(34)*p(232)
      p(1756) = p(1)*p(1295)
      p(1757) = p(34)*p(272)
      p(1758) = p(1)*p(1079)
      p(1759) = p(1)*p(1297)
      p(1760) = p(1)*p(1080)
      p(1761) = p(34)*p(234)
      p(1762) = p(1)*p(1298)
      p(1763) = p(34)*p(235)
      p(1764) = p(1)*p(1300)
      p(1765) = p(1)*p(1301)
      p(1766) = p(34)*p(273)
      p(1767) = p(1)*p(1082)
      p(1768) = p(1)*p(1305)
      p(1769) = p(1)*p(1307)
      p(1770) = p(1)*p(1083)
      p(1771) = p(1)*p(1309)
      p(1772) = p(1)*p(1084)
      p(1773) = p(34)*p(275)
      p(1774) = p(1)*p(1310)
      p(1775) = p(34)*p(276)
      p(1776) = p(1)*p(1312)
      p(1777) = p(1)*p(1313)
      p(1778) = p(34)*p(277)
      p(1779) = p(1)*p(1316)
      p(1780) = p(1)*p(1318)
      p(1781) = p(1)*p(1319)
      p(1782) = p(34)*p(289)
      p(1783) = p(1)*p(1086)
      p(1784) = p(1)*p(1087)
      p(1785) = p(1)*p(1088)
      p(1786) = p(1)*p(1089)
      p(1787) = p(1)*p(1090)
      p(1788) = p(1)*p(1091)
      p(1789) = p(1)*p(1092)
      p(1790) = p(1)*p(1093)
      p(1791) = p(1)*p(1094)
      p(1792) = p(1)*p(1095)
      p(1793) = p(1)*p(1096)
      p(1794) = p(1)*p(1097)
      p(1795) = p(1)*p(1098)
      p(1796) = p(1)*p(1099)
      p(1797) = p(1)*p(1100)
      p(1798) = p(1)*p(1101)
      p(1799) = p(1)*p(1102)
      p(1800) = p(1)*p(1105)
      p(1801) = p(3)*p(1262)
      p(1802) = p(1)*p(1107)
      p(1803) = p(1)*p(1108)
      p(1804) = p(1)*p(1109)
      p(1805) = p(1)*p(1324)
      p(1806) = p(1)*p(1110)
      p(1807) = p(1)*p(1112)
      p(1808) = p(3)*p(1269)
      p(1809) = p(1)*p(1114)
      p(1810) = p(3)*p(1271)
      p(1811) = p(1)*p(1325)
      p(1812) = p(1)*p(1116)
      p(1813) = p(3)*p(1274)
      p(1814) = p(3)*p(1275)
      p(1815) = p(1)*p(1119)
      p(1816) = p(3)*p(1277)
      p(1817) = p(3)*p(1278)
      p(1818) = p(1)*p(1122)
      p(1819) = p(3)*p(1280)
      p(1820) = p(3)*p(1281)
      p(1821) = p(1)*p(1326)
      p(1822) = p(3)*p(1283)
      p(1823) = p(1)*p(1125)
      p(1824) = p(1)*p(1126)
      p(1825) = p(1)*p(1328)
      p(1826) = p(1)*p(1127)
      p(1827) = p(3)*p(1288)
      p(1828) = p(1)*p(1129)
      p(1829) = p(3)*p(1290)
      p(1830) = p(1)*p(1329)
      p(1831) = p(3)*p(1292)
      p(1832) = p(1)*p(1131)
      p(1833) = p(3)*p(1294)
      p(1834) = p(3)*p(1295)
      p(1835) = p(1)*p(1134)
      p(1836) = p(3)*p(1297)
      p(1837) = p(3)*p(1298)
      p(1838) = p(1)*p(1330)
      p(1839) = p(3)*p(1300)
      p(1840) = p(3)*p(1301)
      p(1841) = p(1)*p(1137)
      p(1842) = p(1)*p(1332)
      p(1843) = p(1)*p(1138)
      p(1844) = p(3)*p(1305)
      p(1845) = p(1)*p(1333)
      p(1846) = p(3)*p(1307)
      p(1847) = p(1)*p(1140)
      p(1848) = p(3)*p(1309)
      p(1849) = p(3)*p(1310)
      p(1850) = p(1)*p(1334)
      p(1851) = p(3)*p(1312)
      p(1852) = p(3)*p(1313)
      p(1853) = p(1)*p(1336)
      p(1854) = p(1)*p(1337)
      p(1855) = p(3)*p(1316)
      p(1856) = p(1)*p(1338)
      p(1857) = p(3)*p(1318)
      p(1858) = p(3)*p(1319)
      p(1859) = p(1)*p(1143)
      p(1860) = p(1)*p(1144)
      p(1861) = p(1)*p(1145)
      p(1862) = p(1)*p(1146)
      p(1863) = p(1)*p(1147)
      p(1864) = p(1)*p(1148)
      p(1865) = p(1)*p(1149)
      p(1866) = p(1)*p(1150)
      p(1867) = p(1)*p(1151)
      p(1868) = p(1)*p(1152)
      p(1869) = p(1)*p(1153)
      p(1870) = p(1)*p(1154)
      p(1871) = p(1)*p(1155)
      p(1872) = p(1)*p(1156)
      p(1873) = p(1)*p(1158)
      p(1874) = p(1)*p(1159)
      p(1875) = p(1)*p(1162)
      p(1876) = p(3)*p(1103)
      p(1877) = p(3)*p(1104)
      p(1878) = p(1)*p(1165)
      p(1879) = p(3)*p(1106)
      p(1880) = p(1)*p(1167)
      p(1881) = p(1)*p(1168)
      p(1882) = p(1)*p(1343)
      p(1883) = p(1)*p(1169)
      p(1884) = p(3)*p(1111)
      p(1885) = p(1)*p(1171)
      p(1886) = p(3)*p(1113)
      p(1887) = p(1)*p(1344)
      p(1888) = p(3)*p(1115)
      p(1889) = p(1)*p(1173)
      p(1890) = p(3)*p(1117)
      p(1891) = p(3)*p(1118)
      p(1892) = p(1)*p(1176)
      p(1893) = p(3)*p(1120)
      p(1894) = p(3)*p(1121)
      p(1895) = p(1)*p(1345)
      p(1896) = p(3)*p(1123)
      p(1897) = p(3)*p(1124)
      p(1898) = p(1)*p(1179)
      p(1899) = p(1)*p(1347)
      p(1900) = p(1)*p(1180)
      p(1901) = p(3)*p(1128)
      p(1902) = p(1)*p(1348)
      p(1903) = p(3)*p(1130)
      p(1904) = p(1)*p(1182)
      p(1905) = p(3)*p(1132)
      p(1906) = p(3)*p(1133)
      p(1907) = p(1)*p(1349)
      p(1908) = p(3)*p(1135)
      p(1909) = p(3)*p(1136)
      p(1910) = p(1)*p(1351)
      p(1911) = p(1)*p(1352)
      p(1912) = p(3)*p(1139)
      p(1913) = p(1)*p(1353)
      p(1914) = p(3)*p(1141)
      p(1915) = p(3)*p(1142)
      p(1916) = p(1)*p(1185)
      p(1917) = p(1)*p(1186)
      p(1918) = p(1)*p(1187)
      p(1919) = p(1)*p(1188)
      p(1920) = p(1)*p(1189)
      p(1921) = p(1)*p(1190)
      p(1922) = p(1)*p(1191)
      p(1923) = p(1)*p(1192)
      p(1924) = p(1)*p(1193)
      p(1925) = p(1)*p(1195)
      p(1926) = p(1)*p(1196)
      p(1927) = p(1)*p(1357)
      p(1928) = p(1)*p(1197)
      p(1929) = p(1)*p(1199)
      p(1930) = p(3)*p(1157)
      p(1931) = p(1)*p(1358)
      p(1932) = p(1)*p(1201)
      p(1933) = p(3)*p(1160)
      p(1934) = p(3)*p(1161)
      p(1935) = p(1)*p(1204)
      p(1936) = p(3)*p(1163)
      p(1937) = p(3)*p(1164)
      p(1938) = p(1)*p(1359)
      p(1939) = p(3)*p(1166)
      p(1940) = p(1)*p(1207)
      p(1941) = p(1)*p(1361)
      p(1942) = p(1)*p(1208)
      p(1943) = p(3)*p(1170)
      p(1944) = p(1)*p(1362)
      p(1945) = p(3)*p(1172)
      p(1946) = p(1)*p(1210)
      p(1947) = p(3)*p(1174)
      p(1948) = p(3)*p(1175)
      p(1949) = p(1)*p(1363)
      p(1950) = p(3)*p(1177)
      p(1951) = p(3)*p(1178)
      p(1952) = p(1)*p(1365)
      p(1953) = p(1)*p(1366)
      p(1954) = p(3)*p(1181)
      p(1955) = p(1)*p(1367)
      p(1956) = p(3)*p(1183)
      p(1957) = p(3)*p(1184)
      p(1958) = p(1)*p(1213)
      p(1959) = p(1)*p(1214)
      p(1960) = p(1)*p(1215)
      p(1961) = p(1)*p(1216)
      p(1962) = p(1)*p(1217)
      p(1963) = p(1)*p(1218)
      p(1964) = p(1)*p(1219)
      p(1965) = p(1)*p(1220)
      p(1966) = p(1)*p(1223)
      p(1967) = p(3)*p(1194)
      p(1968) = p(1)*p(1225)
      p(1969) = p(1)*p(1371)
      p(1970) = p(1)*p(1226)
      p(1971) = p(3)*p(1198)
      p(1972) = p(1)*p(1372)
      p(1973) = p(3)*p(1200)
      p(1974) = p(1)*p(1228)
      p(1975) = p(3)*p(1202)
      p(1976) = p(3)*p(1203)
      p(1977) = p(1)*p(1373)
      p(1978) = p(3)*p(1205)
      p(1979) = p(3)*p(1206)
      p(1980) = p(1)*p(1375)
      p(1981) = p(1)*p(1376)
      p(1982) = p(3)*p(1209)
      p(1983) = p(1)*p(1377)
      p(1984) = p(3)*p(1211)
      p(1985) = p(3)*p(1212)
      p(1986) = p(1)*p(1231)
      p(1987) = p(1)*p(1232)
      p(1988) = p(1)*p(1233)
      p(1989) = p(1)*p(1234)
      p(1990) = p(1)*p(1380)
      p(1991) = p(1)*p(1235)
      p(1992) = p(1)*p(1381)
      p(1993) = p(1)*p(1237)
      p(1994) = p(3)*p(1221)
      p(1995) = p(3)*p(1222)
      p(1996) = p(1)*p(1382)
      p(1997) = p(3)*p(1224)
      p(1998) = p(1)*p(1384)
      p(1999) = p(1)*p(1385)
      p(2000) = p(3)*p(1227)
      p(2001) = p(1)*p(1386)
      p(2002) = p(3)*p(1229)
      p(2003) = p(3)*p(1230)
      p(2004) = p(1)*p(1240)
      p(2005) = p(1)*p(1241)
      p(2006) = p(1)*p(1242)
      p(2007) = p(1)*p(1389)
      p(2008) = p(1)*p(1390)
      p(2009) = p(3)*p(1236)
      p(2010) = p(1)*p(1391)
      p(2011) = p(3)*p(1238)
      p(2012) = p(3)*p(1239)
      p(2013) = p(1)*p(1393)
      p(2014) = p(1)*p(1394)
      p(2015) = p(1)*p(1395)
      p(2016) = p(3)*p(1243)
      p(2017) = p(1)*p(1244)
      p(2018) = p(1)*p(1245)
      p(2019) = p(1)*p(1246)
      p(2020) = p(1)*p(1247)
      p(2021) = p(1)*p(1248)
      p(2022) = p(1)*p(1249)
      p(2023) = p(1)*p(1250)
      p(2024) = p(1)*p(1251)
      p(2025) = p(1)*p(1252)
      p(2026) = p(1)*p(1253)
      p(2027) = p(1)*p(1254)
      p(2028) = p(1)*p(1255)
      p(2029) = p(1)*p(1256)
      p(2030) = p(1)*p(1257)
      p(2031) = p(1)*p(1258)
      p(2032) = p(1)*p(1259)
      p(2033) = p(1)*p(1260)
      p(2034) = p(1)*p(1261)
      p(2035) = p(1)*p(1263)
      p(2036) = p(1)*p(1264)
      p(2037) = p(1)*p(1265)
      p(2038) = p(1)*p(1266)
      p(2039) = p(1)*p(1267)
      p(2040) = p(1)*p(1268)
      p(2041) = p(1)*p(1270)
      p(2042) = p(5)*p(750) - p(1734)
      p(2043) = p(1)*p(1272)
      p(2044) = p(1)*p(1273)
      p(2045) = p(1)*p(1276)
      p(2046) = p(5)*p(749) - p(1705)
      p(2047) = p(5)*p(764) - p(1757)
      p(2048) = p(1)*p(1279)
      p(2049) = p(5)*p(752) - p(1740)
      p(2050) = p(5)*p(767) - p(1763)
      p(2051) = p(1)*p(1282)
      p(2052) = p(5)*p(770) - p(1766)
      p(2053) = p(1)*p(1284)
      p(2054) = p(1)*p(1285)
      p(2055) = p(1)*p(1286)
      p(2056) = p(1)*p(1401)
      p(2057) = p(1)*p(1287)
      p(2058) = p(6)*p(850)
      p(2059) = p(1)*p(1289)
      p(2060) = p(60)*p(267)
      p(2061) = p(1)*p(1291)
      p(2062) = p(59)*p(268)
      p(2063) = p(1)*p(1402)
      p(2064) = p(5)*p(851)
      p(2065) = p(1)*p(1293)
      p(2066) = p(7)*p(850)
      p(2067) = p(7)*p(851)
      p(2068) = p(1)*p(1296)
      p(2069) = p(61)*p(267)
      p(2070) = p(61)*p(268)
      p(2071) = p(1)*p(1299)
      p(2072) = p(59)*p(269)
      p(2073) = p(60)*p(269)
      p(2074) = p(1)*p(1403)
      p(2075) = p(5)*p(852)
      p(2076) = p(6)*p(852)
      p(2077) = p(1)*p(1302)
      p(2078) = p(1)*p(1303)
      p(2079) = p(1)*p(1405)
      p(2080) = p(1)*p(1304)
      p(2081) = p(5)*p(774) - p(1755)
      p(2082) = p(1)*p(1306)
      p(2083) = p(5)*p(776) - p(1757)
      p(2084) = p(1)*p(1406)
      p(2085) = p(5)*p(855) - p(2067)
      p(2086) = p(1)*p(1308)
      p(2087) = p(5)*p(778) - p(1761)
      p(2088) = p(6)*p(779) - p(1763)
      p(2089) = p(1)*p(1311)
      p(2090) = p(5)*p(781) - p(1766)
      p(2091) = p(5)*p(788) - p(1782)
      p(2092) = p(1)*p(1407)
      p(2093) = p(5)*p(856) - p(2076)
      p(2094) = p(6)*p(856) - p(2075)
      p(2095) = p(1)*p(1314)
      p(2096) = p(1)*p(1409)
      p(2097) = p(1)*p(1315)
      p(2098) = p(5)*p(785) - p(1775)
      p(2099) = p(1)*p(1410)
      p(2100) = p(5)*p(859) - p(2088)
      p(2101) = p(1)*p(1317)
      p(2102) = p(5)*p(787) - p(1778)
      p(2103) = p(6)*p(788) - p(1778)
      p(2104) = p(1)*p(1411)
      p(2105) = p(5)*p(860) - p(2094)
      p(2106) = p(6)*p(860) - p(2093)
      p(2107) = p(1)*p(1413)
      p(2108) = p(1)*p(1414)
      p(2109) = p(5)*p(863) - p(2103)
      p(2110) = p(1)*p(1415)
      p(2111) = p(5)*p(864) - p(2106)
      p(2112) = p(6)*p(864) - p(2105)
      p(2113) = p(1)*p(1320)
      p(2114) = p(1)*p(1321)
      p(2115) = p(1)*p(1322)
      p(2116) = p(1)*p(1323)
      p(2117) = p(1)*p(1327)
      p(2118) = p(3)*p(1401)
      p(2119) = p(3)*p(1402)
      p(2120) = p(3)*p(1403)
      p(2121) = p(1)*p(1331)
      p(2122) = p(3)*p(1405)
      p(2123) = p(3)*p(1406)
      p(2124) = p(3)*p(1407)
      p(2125) = p(1)*p(1335)
      p(2126) = p(3)*p(1409)
      p(2127) = p(3)*p(1410)
      p(2128) = p(3)*p(1411)
      p(2129) = p(1)*p(1417)
      p(2130) = p(3)*p(1413)
      p(2131) = p(3)*p(1414)
      p(2132) = p(3)*p(1415)
      p(2133) = p(1)*p(1339)
      p(2134) = p(1)*p(1340)
      p(2135) = p(1)*p(1341)
      p(2136) = p(1)*p(1342)
      p(2137) = p(3)*p(1324)
      p(2138) = p(3)*p(1325)
      p(2139) = p(3)*p(1326)
      p(2140) = p(1)*p(1346)
      p(2141) = p(3)*p(1328)
      p(2142) = p(3)*p(1329)
      p(2143) = p(3)*p(1330)
      p(2144) = p(1)*p(1350)
      p(2145) = p(3)*p(1332)
      p(2146) = p(3)*p(1333)
      p(2147) = p(3)*p(1334)
      p(2148) = p(1)*p(1419)
      p(2149) = p(3)*p(1336)
      p(2150) = p(3)*p(1337)
      p(2151) = p(3)*p(1338)
      p(2152) = p(1)*p(1354)
      p(2153) = p(1)*p(1355)
      p(2154) = p(1)*p(1356)
      p(2155) = p(1)*p(1360)
      p(2156) = p(3)*p(1343)
      p(2157) = p(3)*p(1344)
      p(2158) = p(3)*p(1345)
      p(2159) = p(1)*p(1364)
      p(2160) = p(3)*p(1347)
      p(2161) = p(3)*p(1348)
      p(2162) = p(3)*p(1349)
      p(2163) = p(1)*p(1421)
      p(2164) = p(3)*p(1351)
      p(2165) = p(3)*p(1352)
      p(2166) = p(3)*p(1353)
      p(2167) = p(1)*p(1368)
      p(2168) = p(1)*p(1369)
      p(2169) = p(1)*p(1370)
      p(2170) = p(3)*p(1357)
      p(2171) = p(3)*p(1358)
      p(2172) = p(3)*p(1359)
      p(2173) = p(1)*p(1374)
      p(2174) = p(3)*p(1361)
      p(2175) = p(3)*p(1362)
      p(2176) = p(3)*p(1363)
      p(2177) = p(1)*p(1423)
      p(2178) = p(3)*p(1365)
      p(2179) = p(3)*p(1366)
      p(2180) = p(3)*p(1367)
      p(2181) = p(1)*p(1378)
      p(2182) = p(1)*p(1379)
      p(2183) = p(1)*p(1383)
      p(2184) = p(3)*p(1371)
      p(2185) = p(3)*p(1372)
      p(2186) = p(3)*p(1373)
      p(2187) = p(1)*p(1425)
      p(2188) = p(3)*p(1375)
      p(2189) = p(3)*p(1376)
      p(2190) = p(3)*p(1377)
      p(2191) = p(1)*p(1387)
      p(2192) = p(1)*p(1388)
      p(2193) = p(3)*p(1380)
      p(2194) = p(3)*p(1381)
      p(2195) = p(3)*p(1382)
      p(2196) = p(1)*p(1427)
      p(2197) = p(3)*p(1384)
      p(2198) = p(3)*p(1385)
      p(2199) = p(3)*p(1386)
      p(2200) = p(1)*p(1392)
      p(2201) = p(1)*p(1429)
      p(2202) = p(3)*p(1389)
      p(2203) = p(3)*p(1390)
      p(2204) = p(3)*p(1391)
      p(2205) = p(1)*p(1431)
      p(2206) = p(3)*p(1393)
      p(2207) = p(3)*p(1394)
      p(2208) = p(3)*p(1395)
      p(2209) = p(1)*p(1396)
      p(2210) = p(1)*p(1397)
      p(2211) = p(1)*p(1398)
      p(2212) = p(1)*p(1399)
      p(2213) = p(1)*p(1400)
      p(2214) = p(5)*p(850) - p(1719)
      p(2215) = p(6)*p(851) - p(1725)
      p(2216) = p(7)*p(852) - p(1742)
      p(2217) = p(1)*p(1404)
      p(2218) = p(5)*p(854) - p(1753)
      p(2219) = p(6)*p(855) - p(1757)
      p(2220) = p(7)*p(856) - p(1766)
      p(2221) = p(1)*p(1408)
      p(2222) = p(5)*p(858) - p(1773)
      p(2223) = p(6)*p(859) - p(1775)
      p(2224) = p(7)*p(860) - p(1778)
      p(2225) = p(1)*p(1412)
      p(2226) = p(5)*p(862) - p(1782)
      p(2227) = p(6)*p(863) - p(1782)
      p(2228) = p(7)*p(864) - p(1782)
      p(2229) = p(1)*p(1433)
      p(2230) = p(5)*p(880) - p(2112)
      p(2231) = p(6)*p(880) - p(2111)
      p(2232) = p(7)*p(880) - p(2109)
      p(2233) = p(1)*p(1416)
      p(2234) = p(3)*p(1433)
      p(2235) = p(1)*p(1418)
      p(2236) = p(3)*p(1417)
      p(2237) = p(1)*p(1420)
      p(2238) = p(3)*p(1419)
      p(2239) = p(1)*p(1422)
      p(2240) = p(3)*p(1421)
      p(2241) = p(1)*p(1424)
      p(2242) = p(3)*p(1423)
      p(2243) = p(1)*p(1426)
      p(2244) = p(3)*p(1425)
      p(2245) = p(1)*p(1428)
      p(2246) = p(3)*p(1427)
      p(2247) = p(1)*p(1430)
      p(2248) = p(3)*p(1429)
      p(2249) = p(1)*p(1434)
      p(2250) = p(3)*p(1431)
      p(2251) = p(1)*p(1432)
      p(2252) = p(2)*p(1433) - p(2232) - p(2231) - p(2230)
      p(2253) = p(3)*p(1434)

      return
      end subroutine EvPoly

      subroutine evbas
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine eliminates the 2-body terms in the permutationally
!  invariant formulation
!**********************************************************************
      
      integer i
      double precision b1(2254) 

! Pass P(0:2253) to BM1(1:2254)
      do i=1,2254
        b1(i)=p(i-1)
      enddo

! Remove unconnected terms and 2-body terms and pass to B(1:2153)

      b2(1)=b1(5)
 
      do i=2,3
        b2(i)=b1(i+5)
      enddo
 
      b2(4)=b1(10)
 
      do i=5,14
        b2(i)=b1(i+9)
      enddo

      do i=15,16  
        b2(i)=b1(i+10)
      enddo

      b2(17)=b1(28) 

      b2(18)=b1(30)
     
      do i=19,44
        b2(i)=b1(i+15)
      enddo

      do i=45,47
        b2(i)=b1(i+16)  
      enddo

      do i=48,49
        b2(i)=b1(i+17)
      enddo

      b2(50)=b1(68)

      b2(51)=b1(70)

      b2(52)=b1(72)

      do i=53,112
        b2(i)=b1(i+23)
      enddo

      do i=113,115
        b2(i)=b1(i+24)     
      enddo

      do i=116,117
        b2(i)=b1(i+25)
      enddo

      b2(118)=b1(144)

      b2(119)=b1(146)

      b2(120)=b1(148)

      b2(121)=b1(150)

      do i=122,235
        b2(i)=b1(i+32)
      enddo

      do i=236,238
        b2(i)=b1(i+33)
      enddo

      do i=239,241
        b2(i)=b1(i+34)
      enddo

      do i=242,243
        b2(i)=b1(i+35)
      enddo

      b2(244)=b1(280)

      b2(245)=b1(282)

      b2(246)=b1(284)

      b2(247)=b1(286)

      b2(248)=b1(288)

      do i=249,450
        b2(i)=b1(i+43)
      enddo

      do i=451,453
        b2(i)=b1(i+44)
      enddo

      do i=454,456
        b2(i)=b1(i+45)
      enddo

      do i=457,458
        b2(i)=b1(i+46)
      enddo

      b2(459)=b1(506)

      b2(460)=b1(508)

      b2(461)=b1(510)

      b2(462)=b1(512)

      b2(463)=b1(514)

      b2(464)=b1(516)

      do i=465,795
        b2(i)=b1(i+55)
      enddo

      do i=796,798
        b2(i)=b1(i+56)
      enddo

      do i=799,801
        b2(i)=b1(i+57)
      enddo

      do i=802,804
        b2(i)=b1(i+58)
      enddo

      do i=805,806
        b2(i)=b1(i+59)
      enddo

      b2(807)=b1(867)

      b2(808)=b1(869)

      b2(809)=b1(871)

      b2(810)=b1(873)

      b2(811)=b1(875)

      b2(812)=b1(877)

      b2(813)=b1(879)

      do i=814,1332
        b2(i)=b1(i+69)
      enddo

      do i=1333,1335
        b2(i)=b1(i+70)
      enddo

      do i=1336,1338
        b2(i)=b1(i+71)
      enddo

      do i=1339,1341
        b2(i)=b1(i+72)
      enddo

      do i=1342,1343
        b2(i)=b1(i+73)
      enddo

      b2(1344)=b1(1418)

      b2(1345)=b1(1420)

      b2(1346)=b1(1422)
      
      b2(1347)=b1(1424)

      b2(1348)=b1(1426)

      b2(1349)=b1(1428)

      b2(1350)=b1(1430)

      b2(1351)=b1(1432)

      do i=1352,2130
        b2(i)=b1(i+84)
      enddo

      do i=2131,2133
        b2(i)=b1(i+85)
      enddo

      do i=2134,2136
        b2(i)=b1(i+86)
      enddo

      do i=2137,2139
        b2(i)=b1(i+87)
      enddo

      do i=2140,2142
        b2(i)=b1(i+88)
      enddo

      do i=2143,2144
        b2(i)=b1(i+89)
      enddo

      b2(2145)=b1(2235)

      b2(2146)=b1(2237)

      b2(2147)=b1(2239)

      b2(2148)=b1(2241)

      b2(2149)=b1(2243)

      b2(2150)=b1(2245)

      b2(2151)=b1(2247)

      b2(2152)=b1(2249)

      b2(2153)=b1(2251)

      return
      end subroutine evbas

      subroutine evbas2
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine eliminates unused four-body terms
!**********************************************************************

       b(  1) = b2(   1)
       b(  2) = b2(   2)
       b(  3) = b2(   3)
       b(  4) = b2(   4)
       b(  5) = b2(   5)
       b(  6) = b2(   6)
       b(  7) = b2(   7)
       b(  8) = b2(   8)
       b(  9) = b2(   9)
       b( 10) = b2(  10)
       b( 11) = b2(  11)
       b( 12) = b2(  12)
       b( 13) = b2(  13)
       b( 14) = b2(  14)
       b( 15) = b2(  15)
       b( 16) = b2(  16)
       b( 17) = b2(  17)
       b( 18) = b2(  18)
       b( 19) = b2(  19)
       b( 20) = b2(  20)
       b( 21) = b2(  21)
       b( 22) = b2(  22)
       b( 23) = b2(  23)
       b( 24) = b2(  24)
       b( 25) = b2(  25)
       b( 26) = b2(  26)
       b( 27) = b2(  27)
       b( 28) = b2(  28)
       b( 29) = b2(  29)
       b( 30) = b2(  30)
       b( 31) = b2(  31)
       b( 32) = b2(  32)
       b( 33) = b2(  33)
       b( 34) = b2(  34)
       b( 35) = b2(  35)
       b( 36) = b2(  36)
       b( 37) = b2(  37)
       b( 38) = b2(  38)
       b( 39) = b2(  39)
       b( 40) = b2(  40)
       b( 41) = b2(  41)
       b( 42) = b2(  42)
       b( 43) = b2(  43)
       b( 44) = b2(  44)
       b( 45) = b2(  45)
       b( 46) = b2(  46)
       b( 47) = b2(  47)
       b( 48) = b2(  48)
       b( 49) = b2(  49)
       b( 50) = b2(  50)
       b( 51) = b2(  51)
       b( 52) = b2(  52)
       b( 53) = b2(  53)
       b( 54) = b2(  54)
       b( 55) = b2(  55)
       b( 56) = b2(  56)
       b( 57) = b2(  57)
       b( 58) = b2(  58)
       b( 59) = b2(  59)
       b( 60) = b2(  60)
       b( 61) = b2(  61)
       b( 62) = b2(  62)
       b( 63) = b2(  63)
       b( 64) = b2(  64)
       b( 65) = b2(  65)
       b( 66) = b2(  66)
       b( 67) = b2(  67)
       b( 68) = b2(  68)
       b( 69) = b2(  69)
       b( 70) = b2(  70)
       b( 71) = b2(  71)
       b( 72) = b2(  72)
       b( 73) = b2(  73)
       b( 74) = b2(  74)
       b( 75) = b2(  75)
       b( 76) = b2(  76)
       b( 77) = b2(  77)
       b( 78) = b2(  78)
       b( 79) = b2(  79)
       b( 80) = b2(  80)
       b( 81) = b2(  81)
       b( 82) = b2(  82)
       b( 83) = b2(  83)
       b( 84) = b2(  84)
       b( 85) = b2(  85)
       b( 86) = b2(  86)
       b( 87) = b2(  87)
       b( 88) = b2(  88)
       b( 89) = b2(  89)
       b( 90) = b2(  90)
       b( 91) = b2(  91)
       b( 92) = b2(  92)
       b( 93) = b2(  93)
       b( 94) = b2(  94)
       b( 95) = b2(  95)
       b( 96) = b2(  96)
       b( 97) = b2(  97)
       b( 98) = b2(  98)
       b( 99) = b2(  99)
       b(100) = b2( 100)
       b(101) = b2( 101)
       b(102) = b2( 102)
       b(103) = b2( 103)
       b(104) = b2( 104)
       b(105) = b2( 105)
       b(106) = b2( 106)
       b(107) = b2( 107)
       b(108) = b2( 108)
       b(109) = b2( 109)
       b(110) = b2( 110)
       b(111) = b2( 111)
       b(112) = b2( 112)
       b(113) = b2( 113)
       b(114) = b2( 114)
       b(115) = b2( 115)
       b(116) = b2( 116)
       b(117) = b2( 117)
       b(118) = b2( 118)
       b(119) = b2( 119)
       b(120) = b2( 120)
       b(121) = b2( 121)
       b(122) = b2( 122)
       b(123) = b2( 123)
       b(124) = b2( 124)
       b(125) = b2( 125)
       b(126) = b2( 126)
       b(127) = b2( 127)
       b(128) = b2( 128)
       b(129) = b2( 129)
       b(130) = b2( 130)
       b(131) = b2( 131)
       b(132) = b2( 132)
       b(133) = b2( 133)
       b(134) = b2( 134)
       b(135) = b2( 135)
       b(136) = b2( 136)
       b(137) = b2( 137)
       b(138) = b2( 138)
       b(139) = b2( 140)
       b(140) = b2( 179)
       b(141) = b2( 184)
       b(142) = b2( 185)
       b(143) = b2( 193)
       b(144) = b2( 196)
       b(145) = b2( 202)
       b(146) = b2( 233)
       b(147) = b2( 234)
       b(148) = b2( 235)
       b(149) = b2( 237)
       b(150) = b2( 238)
       b(151) = b2( 240)
       b(152) = b2( 241)
       b(153) = b2( 243)
       b(154) = b2( 365)
       b(155) = b2( 370)
       b(156) = b2( 371)
       b(157) = b2( 379)
       b(158) = b2( 382)
       b(159) = b2( 385)
       b(160) = b2( 393)
       b(161) = b2( 396)
       b(162) = b2( 402)
       b(163) = b2( 447)
       b(164) = b2( 448)
       b(165) = b2( 449)
       b(166) = b2( 450)
       b(167) = b2( 452)
       b(168) = b2( 453)
       b(169) = b2( 455)
       b(170) = b2( 456)
       b(171) = b2( 458)
       b(172) = b2( 680)
       b(173) = b2( 685)
       b(174) = b2( 686)
       b(175) = b2( 693)
       b(176) = b2( 694)
       b(177) = b2( 697)
       b(178) = b2( 708)
       b(179) = b2( 711)
       b(180) = b2( 714)
       b(181) = b2( 723)
       b(182) = b2( 726)
       b(183) = b2( 732)
       b(184) = b2( 792)
       b(185) = b2( 793)
       b(186) = b2( 794)
       b(187) = b2( 795)
       b(188) = b2( 797)
       b(189) = b2( 798)
       b(190) = b2( 800)
       b(191) = b2( 801)
       b(192) = b2( 803)
       b(193) = b2( 804)
       b(194) = b2( 806)
       b(195) = b2(1178)
       b(196) = b2(1183)
       b(197) = b2(1184)
       b(198) = b2(1191)
       b(199) = b2(1192)
       b(200) = b2(1193)
       b(201) = b2(1205)
       b(202) = b2(1208)
       b(203) = b2(1211)
       b(204) = b2(1214)
       b(205) = b2(1225)
       b(206) = b2(1228)
       b(207) = b2(1231)
       b(208) = b2(1240)
       b(209) = b2(1243)
       b(210) = b2(1249)
       b(211) = b2(1328)
       b(212) = b2(1329)
       b(213) = b2(1330)
       b(214) = b2(1331)
       b(215) = b2(1332)
       b(216) = b2(1334)
       b(217) = b2(1335)
       b(218) = b2(1337)
       b(219) = b2(1338)
       b(220) = b2(1340)
       b(221) = b2(1341)
       b(222) = b2(1343)
       b(223) = b2(1936)
       b(224) = b2(1941)
       b(225) = b2(1942)
       b(226) = b2(1949)
       b(227) = b2(1950)
       b(228) = b2(1951)
       b(229) = b2(1961)
       b(230) = b2(1962)
       b(231) = b2(1965)
       b(232) = b2(1968)
       b(233) = b2(1982)
       b(234) = b2(1985)
       b(235) = b2(1988)
       b(236) = b2(1991)
       b(237) = b2(2003)
       b(238) = b2(2006)
       b(239) = b2(2009)
       b(240) = b2(2018)
       b(241) = b2(2021)
       b(242) = b2(2027)
       b(243) = b2(2126)
       b(244) = b2(2127)
       b(245) = b2(2128)
       b(246) = b2(2129)
       b(247) = b2(2130)
       b(248) = b2(2132)
       b(249) = b2(2133)
       b(250) = b2(2135)
       b(251) = b2(2136)
       b(252) = b2(2138)
       b(253) = b2(2139)
       b(254) = b2(2141)
       b(255) = b2(2142)
       b(256) = b2(2144)

      end subroutine evbas2

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
! additional N2 potential variables
      double precision :: y, fy
      double precision :: dfdy,rr,dydr

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
! NO pairwise
      else if (imol .eq. 3) then
! Short range NO parameters for NO dissociation
! The long range term comes from dispersion
! The difference between the points of the 8th order MEG NO diatomic
! potential and the dispersion corrections was refitted by a 10th
! order MEG polynomial
        re=   1.1508d0
        de=   149.47844d0   !152.6d0
        c(1)  = -1.385343053807088d-01
        c(2)  =  1.889990874379912d+00
        c(3)  = -4.297653558650669d+00
        c(4)  =  2.143053955670876d+01
        c(5)  = -4.434782706902637d+01
        c(6)  =  4.107240828842039d+01
        c(7)  = -1.090996257595447d+01
        c(8)  = -1.068726081590712d+01
        c(9)  =  9.189735457956482d+00
        c(10) = -2.201965266358976d+00
        c(11) =  8.966018393955543d-01 ! delta b
        c(12) =  2.069542710330378d+00 ! delta c
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
! N2 pairwise
      else if (imol .eq. 2) then
! Parameter for N2 dissociation with SEC (for N2O2)
! This is for generalized Morse equation
! Short range N2 parameters
! Modified parameters for D3(BJ)
        re=1.098d0
        de=  224.9157d0 ! 228.4d0 
        as(0) =  2.7599278840949d0  ! 2.71405774451d0
        as(1) =  0.2318898277373d0  ! 1.32757649829d-1
        as(2) =  0.1908422945648d0  ! 2.66756890408d-1
        as(3) = -0.2727504034613d0  ! 1.95350725241d-1
        as(4) = -0.5345112219335d0  !-4.08663480982d-1
        as(5) =  1.0857331617073d0  ! 3.92451705557d-1
        as(6) =  1.6339897930305d0  ! 1.13006674877d0

      v=0.d0

      y=(r*r*r*r - re*re*re*re)/(r*r*r*r + re*re*re*re)

      fy = as(0) + as(1)*y + as(2)*y*y + as(3)*y*y*y + as(4)*y*y*y*y &
         + as(5)*y*y*y*y*y + as(6)*y*y*y*y*y*y

      u=exp(-fy*(r-re))

      v=de*(1.0d0-u)*(1.0d0-u)-de


! Add D3 dispersion correction
        dist(1)=r
        call d3disp(dist,disp,dispdr,0,imol)
        v=v+disp

! Compute the gradient if needed
      if (igrad.eq.1) then
       grad=0.d0

       dfdy = as(1) + 2.0d0*as(2)*y + 3.0d0*as(3)*y*y &
            + 4.0d0*as(4)*y*y*y + 5.0d0*as(5)*y*y*y*y &
            + 6.0d0*as(6)*y*y*y*y*y
        rr = r*r*r*r + re*re*re*re
        dydr = 8.0d0*r*r*r*re*re*re*re/(rr*rr)
        dfdr = dfdy*dydr
        grad = 2.0d0*de*(1-u)*u*(dfdr*(r-re)+fy)

! Add analytical gradient of D3 dispersion correction
        call d3disp(dist,disp,dispdr,1,imol)
        grad= grad + dispdr(1)

      endif

      endif
      return
      end subroutine ev2gm2 

      subroutine EvdVdR
      use N2O_3Ap_ZV_par
!**********************************************************************
! This subroutine evaluates dVdR for given R 
! dVdR = dV2dR + C*dBdR
! C:            Coefficients, stored in 'dim.inc' 
! P:            Basis functions evaluated for given R
! M:            Monomials evaluated for given R
! dV2dR:        Gradient of 2-body interactions
! dMsdR:        dMsdR(6,6), 6 MEG terms w.r.t. 6 bond length
! dMdR:         dMdR(6,nom), nom monomials w.r.t.6 bond length
! dPdR:         dPdR(6,nob), nop polynomial basis functions 
!               w.r.t. 6 bond lengths
! nom:          number of monomials
! nob:          number of basis functions (polynomials)
! M(nom):       Array to store monomials
! P(nob):       Array to store polynomials
!**********************************************************************
      
      integer i,j
      double precision dist,V,V21,V22,V23,V24,V25,V26
      double precision dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6

! Initialize dVdR(6)
      do i=1,6
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

        dist=R(4)
        call ev2gm2(dist,v24,dv2dr4,3,1)  ! NO
        dVdR(4)=dv2dr4

        dist=R(5)
        call ev2gm2(dist,v25,dv2dr5,3,1)  ! NO
        dVdR(5)=dv2dr5

        dist=R(6)
        call ev2gm2(dist,v26,dv2dr6,2,1)  ! NN
        dVdR(6)=dv2dr6

! Calculate dMEG/dr(6,6) for given R(6)
      call evdmsdr

! Calculate the monomials for each point by using six MEG terms
      call evdmdr

! Calculate the polynomials by using monomials
      call evdpdr 

! Remove 2-body interactions and unconnected terms from polynomials
      call evdbdr

! Eliminate the unused four-body terms
      call evdbdr2

! Evaluate dVdR(6) by taking the product of C(j) and dPdR(i,j)
      do i=1,6      
        do j=1,256
         dVdR(i)=dVdR(i) + c(j)*dBdR(i,j)
        enddo
      enddo

      return
      end subroutine EvdVdR

      subroutine EvdRdX(X)
      use N2O_3Ap_ZV_par
!**********************************************************************
! This subroutine evaluates dRdX for given R and X 
! R:            R(6), 6 bond lengths
! X:            X(12), 12 Cartesian coordinates
! 
! dMdR:         dMdR(6,nom), nom monomials w.r.t.6 bond lengths
! dPdR:         dPdR(6,nob), nop polynomial basis functions 
!               w.r.t. 6 bond lengths
! M(nom):       Array to store monomials
! P(nob):       Array to store polynomials
!**********************************************************************

      integer i,j
      double precision X(12)

! Initialize dRdX(6,12)
      do i=1,6
        do j=1,12
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
      dRdX(3,1)=(x(1)-x(10))/r(3)
      dRdX(3,2)=(x(2)-x(11))/r(3)
      dRdX(3,3)=(x(3)-x(12))/r(3)
      dRdX(3,10)=-dRdX(3,1)
      dRdX(3,11)=-dRdX(3,2)
      dRdX(3,12)=-dRdX(3,3)

! dr4dx
      dRdX(4,4)=(x(4)-x(7))/r(4)
      dRdX(4,5)=(x(5)-x(8))/r(4)
      dRdX(4,6)=(x(6)-x(9))/r(4)
      dRdX(4,7)=-dRdX(4,4)
      dRdX(4,8)=-dRdX(4,5)
      dRdX(4,9)=-dRdX(4,6)

! dr5dx
      dRdX(5,4)=(x(4)-x(10))/r(5)
      dRdX(5,5)=(x(5)-x(11))/r(5)
      dRdX(5,6)=(x(6)-x(12))/r(5)
      dRdX(5,10)=-dRdX(5,4)
      dRdX(5,11)=-dRdX(5,5)
      dRdX(5,12)=-dRdX(5,6)

! dr6dx
      dRdX(6,7)=(x(7)-x(10))/r(6)
      dRdX(6,8)=(x(8)-x(11))/r(6)
      dRdX(6,9)=(x(9)-x(12))/r(6)
      dRdX(6,10)=-dRdX(6,7)
      dRdX(6,11)=-dRdX(6,8)
      dRdX(6,12)=-dRdX(6,9)
! Finish the calculation of non-zero dRdX

      return
      end subroutine EvdRdX

      subroutine EvdMsdR
      use N2O_3Ap_ZV_par
!**********************************************************************
! This subroutine evaluates the derivatives of MEG term X
! w.r.t. interatomic distances R(6)
! dmsdR:        Local variables, dirm(6,6)
! a:            Nonlinear parameter(Angstrom)
! ab:           Nonlinear parameter(Angstrom^2)
! re:           equilibrium bond length(Angstrom)
!**********************************************************************

      integer i,j

! Initialize dmsdr
      do i=1,6
        do j=1,6
          dmsdr(i,j)=0.0d0
        enddo
      enddo

! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
! dmsdr(i,j)=0  i!=j
! dmsdr(i,i)= (-1/a-2*(r-re)/ab)*Exp(-(r-re)/a-(r-re)^2/ab)

      do i=1,6
       dmsdr(i,i)=(-2.0d0*(r(i)-re(i))/ab(i)-1.0d0/a(i)) &
       * Exp(-(r(i)-re(i))/a(i)-((r(i)-re(i))**2.0d0)/ab(i))
      enddo

      return
      end subroutine EvdMsdR

      subroutine EvdMdR
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine takes M(nom) and dMSdR(6,6) and calculates the
!  dMdR(6,nom) that do not have usable decomposition.
!  For degree 10, the number of monomials is nom.
!**********************************************************************

      integer i

      do i=1,6
        dmdr(i,0)  = 0.0d0
        dmdr(i,1)  = dmsdr(i,6)
        dmdr(i,2)  = dmsdr(i,5)
        dmdr(i,3)  = dmsdr(i,4)
        dmdr(i,4)  = dmsdr(i,3)
        dmdr(i,5)  = dmsdr(i,2)
        dmdr(i,6)  = dmsdr(i,1)
        dmdr(i,7)  = dmdr(i,3)*rm(4)  +  rm(3)*dmdr(i,4)
        dmdr(i,8)  = dmdr(i,2)*rm(5)  +  rm(2)*dmdr(i,5)
        dmdr(i,9)  = dmdr(i,2)*rm(4)  +  rm(2)*dmdr(i,4)
        dmdr(i,10) = dmdr(i,3)*rm(5)  +  rm(3)*dmdr(i,5) 
        dmdr(i,11) = dmdr(i,2)*rm(3)  +  rm(2)*dmdr(i,3) 
        dmdr(i,12) = dmdr(i,4)*rm(5)  +  rm(4)*dmdr(i,5) 
        dmdr(i,13) = dmdr(i,2)*rm(7)  +  rm(2)*dmdr(i,7) 
        dmdr(i,14) = dmdr(i,2)*rm(10) +  rm(2)*dmdr(i,10)
        dmdr(i,15) = dmdr(i,2)*rm(12) +  rm(2)*dmdr(i,12)
        dmdr(i,16) = dmdr(i,3)*rm(12) +  rm(3)*dmdr(i,12)
        dmdr(i,17) = dmdr(i,2)*rm(16) +  rm(2)*dmdr(i,16)
      enddo

      return
      end subroutine EvdMdR

      subroutine EvdPdr
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine takes monomials(m) and calculates the
!  permutationally invariant polynomials(p)
!  For degree 10, the number of polynomials is nob.
!**********************************************************************

      integer i

      do i=1,6
      dpdr(i,0) = dmdr(i,0)
      dpdr(i,1) = dmdr(i,1)
      dpdr(i,2) = dmdr(i,2) + dmdr(i,3) + dmdr(i,4) + dmdr(i,5)
      dpdr(i,3) = dmdr(i,6)
      dpdr(i,4) = dpdr(i,1)*p(2) + p(1)*dpdr(i,2)
      dpdr(i,5) = dmdr(i,7) + dmdr(i,8)
      dpdr(i,6) = dmdr(i,9) + dmdr(i,10)
      dpdr(i,7) = dmdr(i,11) + dmdr(i,12)
      dpdr(i,8) = dpdr(i,1)*p(3) + p(1)*dpdr(i,3)
      dpdr(i,9) = dpdr(i,3)*p(2) + p(3)*dpdr(i,2)
      dpdr(i,10) = dpdr(i,1)*p(1) + p(1)*dpdr(i,1)
      dpdr(i,11) = dpdr(i,2)*p(2) + p(2)*dpdr(i,2) & 
                 - dpdr(i,7) - dpdr(i,6) - dpdr(i,5) &
                 - dpdr(i,7) - dpdr(i,6) - dpdr(i,5)
      dpdr(i,12) = dpdr(i,3)*p(3) + p(3)*dpdr(i,3)
      dpdr(i,13) = dpdr(i,1)*p(5) + p(1)*dpdr(i,5)
      dpdr(i,14) = dpdr(i,1)*p(6) + p(1)*dpdr(i,6)
      dpdr(i,15) = dpdr(i,1)*p(7) + p(1)*dpdr(i,7)
      dpdr(i,16) = dmdr(i,13) + dmdr(i,14) + dmdr(i,15) + dmdr(i,16)
      dpdr(i,17) = dpdr(i,1)*p(9) + p(1)*dpdr(i,9) 
      dpdr(i,18) = dpdr(i,3)*p(5) + p(3)*dpdr(i,5) 
      dpdr(i,19) = dpdr(i,3)*p(6) + p(3)*dpdr(i,6) 
      dpdr(i,20) = dpdr(i,3)*p(7) + p(3)*dpdr(i,7) 
      dpdr(i,21) = dpdr(i,1)*p(4) + p(1)*dpdr(i,4) 
      dpdr(i,22) = dpdr(i,1)*p(11) + p(1)*dpdr(i,11)
      dpdr(i,23) = dpdr(i,2)*p(5) + p(2)*dpdr(i,5) - dpdr(i,16)
      dpdr(i,24) = dpdr(i,2)*p(6) + p(2)*dpdr(i,6) - dpdr(i,16)
      dpdr(i,25) = dpdr(i,2)*p(7) + p(2)*dpdr(i,7) - dpdr(i,16)
      dpdr(i,26) = dpdr(i,1)*p(8) + p(1)*dpdr(i,8)
      dpdr(i,27) = dpdr(i,3)*p(11) + p(3)*dpdr(i,11)
      dpdr(i,28) = dpdr(i,1)*p(12) + p(1)*dpdr(i,12)
      dpdr(i,29) = dpdr(i,3)*p(9) + p(3)*dpdr(i,9)
      dpdr(i,30) = dpdr(i,1)*p(10) + p(1)*dpdr(i,10)
      dpdr(i,31) = dpdr(i,2)*p(11) + p(2)*dpdr(i,11) &
                 - dpdr(i,25) - dpdr(i,24) - dpdr(i,23)
      dpdr(i,32) = dpdr(i,3)*p(12) + p(3)*dpdr(i,12)
      dpdr(i,33) = dpdr(i,1)*p(16) + p(1)*dpdr(i,16)
      dpdr(i,34) = dmdr(i,17)
      dpdr(i,35) = dpdr(i,1)*p(18) + p(1)*dpdr(i,18)    
      dpdr(i,36) = dpdr(i,1)*p(19) + p(1)*dpdr(i,19)    
      dpdr(i,37) = dpdr(i,1)*p(20) + p(1)*dpdr(i,20)    
      dpdr(i,38) = dpdr(i,3)*p(16) + p(3)*dpdr(i,16)    
      dpdr(i,39) = dpdr(i,1)*p(13) + p(1)*dpdr(i,13)    
      dpdr(i,40) = dpdr(i,1)*p(14) + p(1)*dpdr(i,14)    
      dpdr(i,41) = dpdr(i,1)*p(15) + p(1)*dpdr(i,15)    
      dpdr(i,42) = dpdr(i,1)*p(23) + p(1)*dpdr(i,23)    
      dpdr(i,43) = dpdr(i,1)*p(24) + p(1)*dpdr(i,24)    
      dpdr(i,44) = dpdr(i,5)*p(6) + p(5)*dpdr(i,6)     
      dpdr(i,45) = dpdr(i,1)*p(25) + p(1)*dpdr(i,25)    
      dpdr(i,46) = dpdr(i,5)*p(7) + p(5)*dpdr(i,7)     
      dpdr(i,47) = dpdr(i,6)*p(7) + p(6)*dpdr(i,7)     
      dpdr(i,48) = dpdr(i,1)*p(17) + p(1)*dpdr(i,17)    
      dpdr(i,49) = dpdr(i,1)*p(27) + p(1)*dpdr(i,27)    
      dpdr(i,50) = dpdr(i,3)*p(23) + p(3)*dpdr(i,23)    
      dpdr(i,51) = dpdr(i,3)*p(24) + p(3)*dpdr(i,24)    
      dpdr(i,52) = dpdr(i,3)*p(25) + p(3)*dpdr(i,25)    
      dpdr(i,53) = dpdr(i,1)*p(29) + p(1)*dpdr(i,29)    
      dpdr(i,54) = dpdr(i,3)*p(18) + p(3)*dpdr(i,18)    
      dpdr(i,55) = dpdr(i,3)*p(19) + p(3)*dpdr(i,19)    
      dpdr(i,56) = dpdr(i,3)*p(20) + p(3)*dpdr(i,20)    
      dpdr(i,57) = dpdr(i,1)*p(21) + p(1)*dpdr(i,21)    
      dpdr(i,58) = dpdr(i,1)*p(22) + p(1)*dpdr(i,22)    
      dpdr(i,59) = dpdr(i,5)*p(5) + p(5)*dpdr(i,5) &
                 - dpdr(i,34) - dpdr(i,34)  
      dpdr(i,60) = dpdr(i,6)*p(6) + p(6)*dpdr(i,6) &
                 - dpdr(i,34) - dpdr(i,34)  
      dpdr(i,61) = dpdr(i,7)*p(7) + p(7)*dpdr(i,7) &
                 - dpdr(i,34) - dpdr(i,34)  
      dpdr(i,62) = dpdr(i,1)*p(31) + p(1)*dpdr(i,31) 
      dpdr(i,63) = dpdr(i,5)*p(11) + p(5)*dpdr(i,11) - dpdr(i,47)  
      dpdr(i,64) = dpdr(i,6)*p(11) + p(6)*dpdr(i,11) - dpdr(i,46)  
      dpdr(i,65) = dpdr(i,7)*p(11) + p(7)*dpdr(i,11) - dpdr(i,44)  
      dpdr(i,66) = dpdr(i,1)*p(26) + p(1)*dpdr(i,26)    
      dpdr(i,67) = dpdr(i,3)*p(31) + p(3)*dpdr(i,31)    
      dpdr(i,68) = dpdr(i,1)*p(28) + p(1)*dpdr(i,28)    
      dpdr(i,69) = dpdr(i,3)*p(27) + p(3)*dpdr(i,27)    
      dpdr(i,70) = dpdr(i,1)*p(32) + p(1)*dpdr(i,32)    
      dpdr(i,71) = dpdr(i,3)*p(29) + p(3)*dpdr(i,29)    
      dpdr(i,72) = dpdr(i,1)*p(30) + p(1)*dpdr(i,30)    
      dpdr(i,73) = dpdr(i,2)*p(31) + p(2)*dpdr(i,31) &
                 - dpdr(i,65) - dpdr(i,64) - dpdr(i,63)  
      dpdr(i,74) = dpdr(i,3)*p(32) + p(3)*dpdr(i,32)    
      dpdr(i,75) = dpdr(i,1)*p(34) + p(1)*dpdr(i,34)    
      dpdr(i,76) = dpdr(i,1)*p(38) + p(1)*dpdr(i,38)    
      dpdr(i,77) = dpdr(i,3)*p(34) + p(3)*dpdr(i,34)    
      dpdr(i,78) = dpdr(i,1)*p(33) + p(1)*dpdr(i,33)    
      dpdr(i,79) = dpdr(i,1)*p(44) + p(1)*dpdr(i,44)    
      dpdr(i,80) = dpdr(i,1)*p(46) + p(1)*dpdr(i,46)    
      dpdr(i,81) = dpdr(i,1)*p(47) + p(1)*dpdr(i,47)    
      dpdr(i,82) = dpdr(i,34)*p(2) + p(34)*dpdr(i,2)    
      dpdr(i,83) = dpdr(i,1)*p(35) + p(1)*dpdr(i,35)    
      dpdr(i,84) = dpdr(i,1)*p(36) + p(1)*dpdr(i,36)    
      dpdr(i,85) = dpdr(i,1)*p(37) + p(1)*dpdr(i,37)    
      dpdr(i,86) = dpdr(i,1)*p(50) + p(1)*dpdr(i,50)    
      dpdr(i,87) = dpdr(i,1)*p(51) + p(1)*dpdr(i,51)    
      dpdr(i,88) = dpdr(i,3)*p(44) + p(3)*dpdr(i,44)    
      dpdr(i,89) = dpdr(i,1)*p(52) + p(1)*dpdr(i,52)    
      dpdr(i,90) = dpdr(i,3)*p(46) + p(3)*dpdr(i,46)    
      dpdr(i,91) = dpdr(i,3)*p(47) + p(3)*dpdr(i,47)    
      dpdr(i,92) = dpdr(i,1)*p(54) + p(1)*dpdr(i,54)    
      dpdr(i,93) = dpdr(i,1)*p(55) + p(1)*dpdr(i,55)    
      dpdr(i,94) = dpdr(i,1)*p(56) + p(1)*dpdr(i,56)    
      dpdr(i,95) = dpdr(i,3)*p(38) + p(3)*dpdr(i,38)    
      dpdr(i,96) = dpdr(i,1)*p(39) + p(1)*dpdr(i,39)    
      dpdr(i,97) = dpdr(i,1)*p(40) + p(1)*dpdr(i,40)    
      dpdr(i,98) = dpdr(i,1)*p(41) + p(1)*dpdr(i,41)    
      dpdr(i,99) = dpdr(i,1)*p(42) + p(1)*dpdr(i,42)    
      dpdr(i,100) = dpdr(i,1)*p(59) + p(1)*dpdr(i,59)    
      dpdr(i,101) = dpdr(i,1)*p(43) + p(1)*dpdr(i,43)    
      dpdr(i,102) = dpdr(i,1)*p(60) + p(1)*dpdr(i,60)    
      dpdr(i,103) = dpdr(i,1)*p(45) + p(1)*dpdr(i,45)    
      dpdr(i,104) = dpdr(i,5)*p(16) + p(5)*dpdr(i,16) - dpdr(i,82)  
      dpdr(i,105) = dpdr(i,6)*p(16) + p(6)*dpdr(i,16) - dpdr(i,82)  
      dpdr(i,106) = dpdr(i,1)*p(61) + p(1)*dpdr(i,61)   
      dpdr(i,107) = dpdr(i,7)*p(16) + p(7)*dpdr(i,16) - dpdr(i,82)  
      dpdr(i,108) = dpdr(i,1)*p(63) + p(1)*dpdr(i,63)   
      dpdr(i,109) = dpdr(i,1)*p(64) + p(1)*dpdr(i,64)   
      dpdr(i,110) = dpdr(i,5)*p(24) + p(5)*dpdr(i,24) - dpdr(i,105)  
      dpdr(i,111) = dpdr(i,1)*p(65) + p(1)*dpdr(i,65)   
      dpdr(i,112) = dpdr(i,5)*p(25) + p(5)*dpdr(i,25) - dpdr(i,107)  
      dpdr(i,113) = dpdr(i,6)*p(25) + p(6)*dpdr(i,25) - dpdr(i,107)  
      dpdr(i,114) = dpdr(i,1)*p(48) + p(1)*dpdr(i,48)    
      dpdr(i,115) = dpdr(i,1)*p(49) + p(1)*dpdr(i,49)    
      dpdr(i,116) = dpdr(i,3)*p(59) + p(3)*dpdr(i,59)    
      dpdr(i,117) = dpdr(i,3)*p(60) + p(3)*dpdr(i,60)    
      dpdr(i,118) = dpdr(i,3)*p(61) + p(3)*dpdr(i,61)    
      dpdr(i,119) = dpdr(i,1)*p(67) + p(1)*dpdr(i,67)    
      dpdr(i,120) = dpdr(i,3)*p(63) + p(3)*dpdr(i,63)    
      dpdr(i,121) = dpdr(i,3)*p(64) + p(3)*dpdr(i,64)    
      dpdr(i,122) = dpdr(i,3)*p(65) + p(3)*dpdr(i,65)    
      dpdr(i,123) = dpdr(i,1)*p(53) + p(1)*dpdr(i,53)    
      dpdr(i,124) = dpdr(i,1)*p(69) + p(1)*dpdr(i,69)    
      dpdr(i,125) = dpdr(i,3)*p(50) + p(3)*dpdr(i,50)    
      dpdr(i,126) = dpdr(i,3)*p(51) + p(3)*dpdr(i,51)    
      dpdr(i,127) = dpdr(i,3)*p(52) + p(3)*dpdr(i,52)    
      dpdr(i,128) = dpdr(i,1)*p(71) + p(1)*dpdr(i,71)    
      dpdr(i,129) = dpdr(i,3)*p(54) + p(3)*dpdr(i,54)    
      dpdr(i,130) = dpdr(i,3)*p(55) + p(3)*dpdr(i,55)    
      dpdr(i,131) = dpdr(i,3)*p(56) + p(3)*dpdr(i,56)    
      dpdr(i,132) = dpdr(i,1)*p(57) + p(1)*dpdr(i,57)    
      dpdr(i,133) = dpdr(i,1)*p(58) + p(1)*dpdr(i,58)    
      dpdr(i,134) = dpdr(i,1)*p(62) + p(1)*dpdr(i,62)    
      dpdr(i,135) = dpdr(i,2)*p(59) + p(2)*dpdr(i,59) - dpdr(i,104)  
      dpdr(i,136) = dpdr(i,2)*p(60) + p(2)*dpdr(i,60) - dpdr(i,105)  
      dpdr(i,137) = dpdr(i,2)*p(61) + p(2)*dpdr(i,61) - dpdr(i,107)  
      dpdr(i,138) = dpdr(i,1)*p(73) + p(1)*dpdr(i,73)   
      dpdr(i,139) = dpdr(i,5)*p(31) + p(5)*dpdr(i,31) - dpdr(i,113)  
      dpdr(i,140) = dpdr(i,6)*p(31) + p(6)*dpdr(i,31) - dpdr(i,112)  
      dpdr(i,141) = dpdr(i,7)*p(31) + p(7)*dpdr(i,31) - dpdr(i,110)  
      dpdr(i,142) = dpdr(i,1)*p(66) + p(1)*dpdr(i,66)    
      dpdr(i,143) = dpdr(i,3)*p(73) + p(3)*dpdr(i,73)    
      dpdr(i,144) = dpdr(i,1)*p(68) + p(1)*dpdr(i,68)    
      dpdr(i,145) = dpdr(i,3)*p(67) + p(3)*dpdr(i,67)    
      dpdr(i,146) = dpdr(i,1)*p(70) + p(1)*dpdr(i,70)    
      dpdr(i,147) = dpdr(i,3)*p(69) + p(3)*dpdr(i,69)    
      dpdr(i,148) = dpdr(i,1)*p(74) + p(1)*dpdr(i,74)    
      dpdr(i,149) = dpdr(i,3)*p(71) + p(3)*dpdr(i,71)    
      dpdr(i,150) = dpdr(i,1)*p(72) + p(1)*dpdr(i,72)    
      dpdr(i,151) = dpdr(i,2)*p(73) + p(2)*dpdr(i,73) &
                  - dpdr(i,141) - dpdr(i,140) - dpdr(i,139)  
      dpdr(i,152) = dpdr(i,3)*p(74) + p(3)*dpdr(i,74)    
      dpdr(i,153) = dpdr(i,1)*p(77) + p(1)*dpdr(i,77)    
      dpdr(i,154) = dpdr(i,1)*p(75) + p(1)*dpdr(i,75)    
      dpdr(i,155) = dpdr(i,1)*p(82) + p(1)*dpdr(i,82)    
      dpdr(i,156) = dpdr(i,1)*p(76) + p(1)*dpdr(i,76)    
      dpdr(i,157) = dpdr(i,1)*p(88) + p(1)*dpdr(i,88)    
      dpdr(i,158) = dpdr(i,1)*p(90) + p(1)*dpdr(i,90)    
      dpdr(i,159) = dpdr(i,1)*p(91) + p(1)*dpdr(i,91)    
      dpdr(i,160) = dpdr(i,3)*p(82) + p(3)*dpdr(i,82)    
      dpdr(i,161) = dpdr(i,1)*p(95) + p(1)*dpdr(i,95)    
      dpdr(i,162) = dpdr(i,3)*p(77) + p(3)*dpdr(i,77)    
      dpdr(i,163) = dpdr(i,1)*p(78) + p(1)*dpdr(i,78)    
      dpdr(i,164) = dpdr(i,1)*p(79) + p(1)*dpdr(i,79)    
      dpdr(i,165) = dpdr(i,1)*p(80) + p(1)*dpdr(i,80)    
      dpdr(i,166) = dpdr(i,1)*p(104) + p(1)*dpdr(i,104)   
      dpdr(i,167) = dpdr(i,1)*p(81) + p(1)*dpdr(i,81)    
      dpdr(i,168) = dpdr(i,34)*p(5) + p(34)*dpdr(i,5)    
      dpdr(i,169) = dpdr(i,1)*p(105) + p(1)*dpdr(i,105)   
      dpdr(i,170) = dpdr(i,34)*p(6) + p(34)*dpdr(i,6)    
      dpdr(i,171) = dpdr(i,1)*p(107) + p(1)*dpdr(i,107)   
      dpdr(i,172) = dpdr(i,34)*p(7) + p(34)*dpdr(i,7)    
      dpdr(i,173) = dpdr(i,1)*p(110) + p(1)*dpdr(i,110)   
      dpdr(i,174) = dpdr(i,1)*p(112) + p(1)*dpdr(i,112)   
      dpdr(i,175) = dpdr(i,1)*p(113) + p(1)*dpdr(i,113)   
      dpdr(i,176) = dpdr(i,34)*p(11) + p(34)*dpdr(i,11)   
      dpdr(i,177) = dpdr(i,1)*p(83) + p(1)*dpdr(i,83)    
      dpdr(i,178) = dpdr(i,1)*p(84) + p(1)*dpdr(i,84)    
      dpdr(i,179) = dpdr(i,1)*p(85) + p(1)*dpdr(i,85)    
      dpdr(i,180) = dpdr(i,1)*p(86) + p(1)*dpdr(i,86)    
      dpdr(i,181) = dpdr(i,1)*p(116) + p(1)*dpdr(i,116)   
      dpdr(i,182) = dpdr(i,1)*p(87) + p(1)*dpdr(i,87)    
      dpdr(i,183) = dpdr(i,1)*p(117) + p(1)*dpdr(i,117)   
      dpdr(i,184) = dpdr(i,1)*p(89) + p(1)*dpdr(i,89)    
      dpdr(i,185) = dpdr(i,3)*p(104) + p(3)*dpdr(i,104)   
      dpdr(i,186) = dpdr(i,3)*p(105) + p(3)*dpdr(i,105)   
      dpdr(i,187) = dpdr(i,1)*p(118) + p(1)*dpdr(i,118)   
      dpdr(i,188) = dpdr(i,3)*p(107) + p(3)*dpdr(i,107)   
      dpdr(i,189) = dpdr(i,1)*p(120) + p(1)*dpdr(i,120)   
      dpdr(i,190) = dpdr(i,1)*p(121) + p(1)*dpdr(i,121)   
      dpdr(i,191) = dpdr(i,3)*p(110) + p(3)*dpdr(i,110)   
      dpdr(i,192) = dpdr(i,1)*p(122) + p(1)*dpdr(i,122)   
      dpdr(i,193) = dpdr(i,3)*p(112) + p(3)*dpdr(i,112)   
      dpdr(i,194) = dpdr(i,3)*p(113) + p(3)*dpdr(i,113)   
      dpdr(i,195) = dpdr(i,1)*p(92) + p(1)*dpdr(i,92)    
      dpdr(i,196) = dpdr(i,1)*p(93) + p(1)*dpdr(i,93)    
      dpdr(i,197) = dpdr(i,1)*p(94) + p(1)*dpdr(i,94)    
      dpdr(i,198) = dpdr(i,1)*p(125) + p(1)*dpdr(i,125)   
      dpdr(i,199) = dpdr(i,1)*p(126) + p(1)*dpdr(i,126)   
      dpdr(i,200) = dpdr(i,3)*p(88) + p(3)*dpdr(i,88)    
      dpdr(i,201) = dpdr(i,1)*p(127) + p(1)*dpdr(i,127)   
      dpdr(i,202) = dpdr(i,3)*p(90) + p(3)*dpdr(i,90)    
      dpdr(i,203) = dpdr(i,3)*p(91) + p(3)*dpdr(i,91)    
      dpdr(i,204) = dpdr(i,1)*p(129) + p(1)*dpdr(i,129)   
      dpdr(i,205) = dpdr(i,1)*p(130) + p(1)*dpdr(i,130)   
      dpdr(i,206) = dpdr(i,1)*p(131) + p(1)*dpdr(i,131)   
      dpdr(i,207) = dpdr(i,3)*p(95) + p(3)*dpdr(i,95)    
      dpdr(i,208) = dpdr(i,1)*p(96) + p(1)*dpdr(i,96)    
      dpdr(i,209) = dpdr(i,1)*p(97) + p(1)*dpdr(i,97)    
      dpdr(i,210) = dpdr(i,1)*p(98) + p(1)*dpdr(i,98)    
      dpdr(i,211) = dpdr(i,1)*p(99) + p(1)*dpdr(i,99)    
      dpdr(i,212) = dpdr(i,1)*p(100) + p(1)*dpdr(i,100)   
      dpdr(i,213) = dpdr(i,1)*p(101) + p(1)*dpdr(i,101)   
      dpdr(i,214) = dpdr(i,1)*p(102) + p(1)*dpdr(i,102)   
      dpdr(i,215) = dpdr(i,1)*p(103) + p(1)*dpdr(i,103)   
      dpdr(i,216) = dpdr(i,1)*p(106) + p(1)*dpdr(i,106)   
      dpdr(i,217) = dpdr(i,5)*p(47) + p(5)*dpdr(i,47) - dpdr(i,176)  
      dpdr(i,218) = dpdr(i,1)*p(108) + p(1)*dpdr(i,108)   
      dpdr(i,219) = dpdr(i,1)*p(135) + p(1)*dpdr(i,135)   
      dpdr(i,220) = dpdr(i,1)*p(109) + p(1)*dpdr(i,109)   
      dpdr(i,221) = dpdr(i,6)*p(59) + p(6)*dpdr(i,59)    
      dpdr(i,222) = dpdr(i,1)*p(136) + p(1)*dpdr(i,136)   
      dpdr(i,223) = dpdr(i,5)*p(60) + p(5)*dpdr(i,60)    
      dpdr(i,224) = dpdr(i,1)*p(111) + p(1)*dpdr(i,111)   
      dpdr(i,225) = dpdr(i,7)*p(59) + p(7)*dpdr(i,59)    
      dpdr(i,226) = dpdr(i,7)*p(60) + p(7)*dpdr(i,60)    
      dpdr(i,227) = dpdr(i,1)*p(137) + p(1)*dpdr(i,137)   
      dpdr(i,228) = dpdr(i,5)*p(61) + p(5)*dpdr(i,61)    
      dpdr(i,229) = dpdr(i,6)*p(61) + p(6)*dpdr(i,61)    
      dpdr(i,230) = dpdr(i,1)*p(139) + p(1)*dpdr(i,139)   
      dpdr(i,231) = dpdr(i,1)*p(140) + p(1)*dpdr(i,140)   
      dpdr(i,232) = dpdr(i,5)*p(64) + p(5)*dpdr(i,64) - dpdr(i,226)  
      dpdr(i,233) = dpdr(i,1)*p(141) + p(1)*dpdr(i,141) 
      dpdr(i,234) = dpdr(i,5)*p(65) + p(5)*dpdr(i,65) - dpdr(i,229)  
      dpdr(i,235) = dpdr(i,6)*p(65) + p(6)*dpdr(i,65) - dpdr(i,228)  
      dpdr(i,236) = dpdr(i,1)*p(114) + p(1)*dpdr(i,114)   
      dpdr(i,237) = dpdr(i,1)*p(115) + p(1)*dpdr(i,115)   
      dpdr(i,238) = dpdr(i,1)*p(119) + p(1)*dpdr(i,119)   
      dpdr(i,239) = dpdr(i,3)*p(135) + p(3)*dpdr(i,135)   
      dpdr(i,240) = dpdr(i,3)*p(136) + p(3)*dpdr(i,136)   
      dpdr(i,241) = dpdr(i,3)*p(137) + p(3)*dpdr(i,137)   
      dpdr(i,242) = dpdr(i,1)*p(143) + p(1)*dpdr(i,143)   
      dpdr(i,243) = dpdr(i,3)*p(139) + p(3)*dpdr(i,139)   
      dpdr(i,244) = dpdr(i,3)*p(140) + p(3)*dpdr(i,140)   
      dpdr(i,245) = dpdr(i,3)*p(141) + p(3)*dpdr(i,141)   
      dpdr(i,246) = dpdr(i,1)*p(123) + p(1)*dpdr(i,123)   
      dpdr(i,247) = dpdr(i,1)*p(124) + p(1)*dpdr(i,124)   
      dpdr(i,248) = dpdr(i,3)*p(116) + p(3)*dpdr(i,116)   
      dpdr(i,249) = dpdr(i,3)*p(117) + p(3)*dpdr(i,117)   
      dpdr(i,250) = dpdr(i,3)*p(118) + p(3)*dpdr(i,118)   
      dpdr(i,251) = dpdr(i,1)*p(145) + p(1)*dpdr(i,145)   
      dpdr(i,252) = dpdr(i,3)*p(120) + p(3)*dpdr(i,120)   
      dpdr(i,253) = dpdr(i,3)*p(121) + p(3)*dpdr(i,121)   
      dpdr(i,254) = dpdr(i,3)*p(122) + p(3)*dpdr(i,122)   
      dpdr(i,255) = dpdr(i,1)*p(128) + p(1)*dpdr(i,128)   
      dpdr(i,256) = dpdr(i,1)*p(147) + p(1)*dpdr(i,147)   
      dpdr(i,257) = dpdr(i,3)*p(125) + p(3)*dpdr(i,125)   
      dpdr(i,258) = dpdr(i,3)*p(126) + p(3)*dpdr(i,126)   
      dpdr(i,259) = dpdr(i,3)*p(127) + p(3)*dpdr(i,127)   
      dpdr(i,260) = dpdr(i,1)*p(149) + p(1)*dpdr(i,149)   
      dpdr(i,261) = dpdr(i,3)*p(129) + p(3)*dpdr(i,129)   
      dpdr(i,262) = dpdr(i,3)*p(130) + p(3)*dpdr(i,130)   
      dpdr(i,263) = dpdr(i,3)*p(131) + p(3)*dpdr(i,131)   
      dpdr(i,264) = dpdr(i,1)*p(132) + p(1)*dpdr(i,132)   
      dpdr(i,265) = dpdr(i,1)*p(133) + p(1)*dpdr(i,133)   
      dpdr(i,266) = dpdr(i,1)*p(134) + p(1)*dpdr(i,134)   
      dpdr(i,267) = dpdr(i,5)*p(59) + p(5)*dpdr(i,59) - dpdr(i,168)  
      dpdr(i,268) = dpdr(i,6)*p(60) + p(6)*dpdr(i,60) - dpdr(i,170)  
      dpdr(i,269) = dpdr(i,7)*p(61) + p(7)*dpdr(i,61) - dpdr(i,172)  
      dpdr(i,270) = dpdr(i,1)*p(138) + p(1)*dpdr(i,138)    
      dpdr(i,271) = dpdr(i,5)*p(63) + p(5)*dpdr(i,63) - dpdr(i,176)  
      dpdr(i,272) = dpdr(i,6)*p(64) + p(6)*dpdr(i,64) - dpdr(i,176)  
      dpdr(i,273) = dpdr(i,7)*p(65) + p(7)*dpdr(i,65) - dpdr(i,176)  
      dpdr(i,274) = dpdr(i,1)*p(151) + p(1)*dpdr(i,151)    
      dpdr(i,275) = dpdr(i,5)*p(73) + p(5)*dpdr(i,73) - dpdr(i,235)  
      dpdr(i,276) = dpdr(i,6)*p(73) + p(6)*dpdr(i,73) - dpdr(i,234)  
      dpdr(i,277) = dpdr(i,7)*p(73) + p(7)*dpdr(i,73) - dpdr(i,232)  
      dpdr(i,278) = dpdr(i,1)*p(142) + p(1)*dpdr(i,142)   
      dpdr(i,279) = dpdr(i,3)*p(151) + p(3)*dpdr(i,151)   
      dpdr(i,280) = dpdr(i,1)*p(144) + p(1)*dpdr(i,144)   
      dpdr(i,281) = dpdr(i,3)*p(143) + p(3)*dpdr(i,143)   
      dpdr(i,282) = dpdr(i,1)*p(146) + p(1)*dpdr(i,146)   
      dpdr(i,283) = dpdr(i,3)*p(145) + p(3)*dpdr(i,145)   
      dpdr(i,284) = dpdr(i,1)*p(148) + p(1)*dpdr(i,148)   
      dpdr(i,285) = dpdr(i,3)*p(147) + p(3)*dpdr(i,147)   
      dpdr(i,286) = dpdr(i,1)*p(152) + p(1)*dpdr(i,152)   
      dpdr(i,287) = dpdr(i,3)*p(149) + p(3)*dpdr(i,149)   
      dpdr(i,288) = dpdr(i,1)*p(150) + p(1)*dpdr(i,150)   
      dpdr(i,289) = dpdr(i,2)*p(151) + p(2)*dpdr(i,151) &
                  - dpdr(i,277) - dpdr(i,276) - dpdr(i,275)  
      dpdr(i,290) = dpdr(i,3)*p(152) + p(3)*dpdr(i,152)   
      dpdr(i,291) = dpdr(i,1)*p(153) + p(1)*dpdr(i,153)   
      dpdr(i,292) = dpdr(i,1)*p(160) + p(1)*dpdr(i,160)   
      dpdr(i,293) = dpdr(i,1)*p(162) + p(1)*dpdr(i,162)   
      dpdr(i,294) = dpdr(i,1)*p(154) + p(1)*dpdr(i,154)   
      dpdr(i,295) = dpdr(i,1)*p(155) + p(1)*dpdr(i,155)   
      dpdr(i,296) = dpdr(i,1)*p(168) + p(1)*dpdr(i,168)   
      dpdr(i,297) = dpdr(i,1)*p(170) + p(1)*dpdr(i,170)   
      dpdr(i,298) = dpdr(i,1)*p(172) + p(1)*dpdr(i,172)   
      dpdr(i,299) = dpdr(i,1)*p(176) + p(1)*dpdr(i,176)   
      dpdr(i,300) = dpdr(i,1)*p(156) + p(1)*dpdr(i,156)   
      dpdr(i,301) = dpdr(i,1)*p(157) + p(1)*dpdr(i,157)   
      dpdr(i,302) = dpdr(i,1)*p(158) + p(1)*dpdr(i,158)   
      dpdr(i,303) = dpdr(i,1)*p(185) + p(1)*dpdr(i,185)   
      dpdr(i,304) = dpdr(i,1)*p(159) + p(1)*dpdr(i,159)   
      dpdr(i,305) = dpdr(i,3)*p(168) + p(3)*dpdr(i,168)   
      dpdr(i,306) = dpdr(i,1)*p(186) + p(1)*dpdr(i,186)   
      dpdr(i,307) = dpdr(i,3)*p(170) + p(3)*dpdr(i,170)   
      dpdr(i,308) = dpdr(i,1)*p(188) + p(1)*dpdr(i,188)   
      dpdr(i,309) = dpdr(i,3)*p(172) + p(3)*dpdr(i,172)   
      dpdr(i,310) = dpdr(i,1)*p(191) + p(1)*dpdr(i,191)   
      dpdr(i,311) = dpdr(i,1)*p(193) + p(1)*dpdr(i,193)   
      dpdr(i,312) = dpdr(i,1)*p(194) + p(1)*dpdr(i,194)   
      dpdr(i,313) = dpdr(i,3)*p(176) + p(3)*dpdr(i,176)   
      dpdr(i,314) = dpdr(i,1)*p(161) + p(1)*dpdr(i,161)   
      dpdr(i,315) = dpdr(i,1)*p(200) + p(1)*dpdr(i,200)   
      dpdr(i,316) = dpdr(i,1)*p(202) + p(1)*dpdr(i,202)   
      dpdr(i,317) = dpdr(i,1)*p(203) + p(1)*dpdr(i,203)   
      dpdr(i,318) = dpdr(i,3)*p(160) + p(3)*dpdr(i,160)   
      dpdr(i,319) = dpdr(i,1)*p(207) + p(1)*dpdr(i,207)   
      dpdr(i,320) = dpdr(i,3)*p(162) + p(3)*dpdr(i,162)   
      dpdr(i,321) = dpdr(i,1)*p(163) + p(1)*dpdr(i,163)   
      dpdr(i,322) = dpdr(i,1)*p(164) + p(1)*dpdr(i,164)   
      dpdr(i,323) = dpdr(i,1)*p(165) + p(1)*dpdr(i,165)   
      dpdr(i,324) = dpdr(i,1)*p(166) + p(1)*dpdr(i,166)   
      dpdr(i,325) = dpdr(i,1)*p(167) + p(1)*dpdr(i,167)   
      dpdr(i,326) = dpdr(i,1)*p(169) + p(1)*dpdr(i,169)   
      dpdr(i,327) = dpdr(i,1)*p(171) + p(1)*dpdr(i,171)   
      dpdr(i,328) = dpdr(i,1)*p(217) + p(1)*dpdr(i,217)   
      dpdr(i,329) = dpdr(i,34)*p(16) + p(34)*dpdr(i,16)   
      dpdr(i,330) = dpdr(i,1)*p(173) + p(1)*dpdr(i,173)   
      dpdr(i,331) = dpdr(i,1)*p(221) + p(1)*dpdr(i,221)   
      dpdr(i,332) = dpdr(i,1)*p(223) + p(1)*dpdr(i,223)   
      dpdr(i,333) = dpdr(i,1)*p(174) + p(1)*dpdr(i,174)   
      dpdr(i,334) = dpdr(i,1)*p(225) + p(1)*dpdr(i,225)   
      dpdr(i,335) = dpdr(i,1)*p(175) + p(1)*dpdr(i,175)   
      dpdr(i,336) = dpdr(i,34)*p(23) + p(34)*dpdr(i,23)   
      dpdr(i,337) = dpdr(i,1)*p(226) + p(1)*dpdr(i,226)   
      dpdr(i,338) = dpdr(i,34)*p(24) + p(34)*dpdr(i,24)   
      dpdr(i,339) = dpdr(i,1)*p(228) + p(1)*dpdr(i,228)   
      dpdr(i,340) = dpdr(i,1)*p(229) + p(1)*dpdr(i,229)   
      dpdr(i,341) = dpdr(i,34)*p(25) + p(34)*dpdr(i,25)   
      dpdr(i,342) = dpdr(i,1)*p(232) + p(1)*dpdr(i,232)   
      dpdr(i,343) = dpdr(i,1)*p(234) + p(1)*dpdr(i,234)   
      dpdr(i,344) = dpdr(i,1)*p(235) + p(1)*dpdr(i,235)   
      dpdr(i,345) = dpdr(i,34)*p(31) + p(34)*dpdr(i,31)   
      dpdr(i,346) = dpdr(i,1)*p(177) + p(1)*dpdr(i,177)   
      dpdr(i,347) = dpdr(i,1)*p(178) + p(1)*dpdr(i,178)   
      dpdr(i,348) = dpdr(i,1)*p(179) + p(1)*dpdr(i,179)   
      dpdr(i,349) = dpdr(i,1)*p(180) + p(1)*dpdr(i,180)   
      dpdr(i,350) = dpdr(i,1)*p(181) + p(1)*dpdr(i,181)   
      dpdr(i,351) = dpdr(i,1)*p(182) + p(1)*dpdr(i,182)   
      dpdr(i,352) = dpdr(i,1)*p(183) + p(1)*dpdr(i,183)   
      dpdr(i,353) = dpdr(i,1)*p(184) + p(1)*dpdr(i,184)   
      dpdr(i,354) = dpdr(i,1)*p(187) + p(1)*dpdr(i,187)   
      dpdr(i,355) = dpdr(i,3)*p(217) + p(3)*dpdr(i,217)   
      dpdr(i,356) = dpdr(i,1)*p(189) + p(1)*dpdr(i,189)   
      dpdr(i,357) = dpdr(i,1)*p(239) + p(1)*dpdr(i,239)   
      dpdr(i,358) = dpdr(i,1)*p(190) + p(1)*dpdr(i,190)   
      dpdr(i,359) = dpdr(i,3)*p(221) + p(3)*dpdr(i,221)   
      dpdr(i,360) = dpdr(i,1)*p(240) + p(1)*dpdr(i,240)   
      dpdr(i,361) = dpdr(i,3)*p(223) + p(3)*dpdr(i,223)   
      dpdr(i,362) = dpdr(i,1)*p(192) + p(1)*dpdr(i,192)   
      dpdr(i,363) = dpdr(i,3)*p(225) + p(3)*dpdr(i,225)   
      dpdr(i,364) = dpdr(i,3)*p(226) + p(3)*dpdr(i,226)   
      dpdr(i,365) = dpdr(i,1)*p(241) + p(1)*dpdr(i,241)   
      dpdr(i,366) = dpdr(i,3)*p(228) + p(3)*dpdr(i,228)   
      dpdr(i,367) = dpdr(i,3)*p(229) + p(3)*dpdr(i,229)   
      dpdr(i,368) = dpdr(i,1)*p(243) + p(1)*dpdr(i,243)   
      dpdr(i,369) = dpdr(i,1)*p(244) + p(1)*dpdr(i,244)   
      dpdr(i,370) = dpdr(i,3)*p(232) + p(3)*dpdr(i,232)   
      dpdr(i,371) = dpdr(i,1)*p(245) + p(1)*dpdr(i,245)   
      dpdr(i,372) = dpdr(i,3)*p(234) + p(3)*dpdr(i,234)   
      dpdr(i,373) = dpdr(i,3)*p(235) + p(3)*dpdr(i,235)   
      dpdr(i,374) = dpdr(i,1)*p(195) + p(1)*dpdr(i,195)   
      dpdr(i,375) = dpdr(i,1)*p(196) + p(1)*dpdr(i,196)   
      dpdr(i,376) = dpdr(i,1)*p(197) + p(1)*dpdr(i,197)   
      dpdr(i,377) = dpdr(i,1)*p(198) + p(1)*dpdr(i,198)   
      dpdr(i,378) = dpdr(i,1)*p(248) + p(1)*dpdr(i,248)   
      dpdr(i,379) = dpdr(i,1)*p(199) + p(1)*dpdr(i,199)   
      dpdr(i,380) = dpdr(i,1)*p(249) + p(1)*dpdr(i,249)   
      dpdr(i,381) = dpdr(i,1)*p(201) + p(1)*dpdr(i,201)   
      dpdr(i,382) = dpdr(i,3)*p(185) + p(3)*dpdr(i,185)   
      dpdr(i,383) = dpdr(i,3)*p(186) + p(3)*dpdr(i,186)   
      dpdr(i,384) = dpdr(i,1)*p(250) + p(1)*dpdr(i,250)   
      dpdr(i,385) = dpdr(i,3)*p(188) + p(3)*dpdr(i,188)   
      dpdr(i,386) = dpdr(i,1)*p(252) + p(1)*dpdr(i,252)   
      dpdr(i,387) = dpdr(i,1)*p(253) + p(1)*dpdr(i,253)   
      dpdr(i,388) = dpdr(i,3)*p(191) + p(3)*dpdr(i,191)   
      dpdr(i,389) = dpdr(i,1)*p(254) + p(1)*dpdr(i,254)   
      dpdr(i,390) = dpdr(i,3)*p(193) + p(3)*dpdr(i,193)   
      dpdr(i,391) = dpdr(i,3)*p(194) + p(3)*dpdr(i,194)   
      dpdr(i,392) = dpdr(i,1)*p(204) + p(1)*dpdr(i,204)   
      dpdr(i,393) = dpdr(i,1)*p(205) + p(1)*dpdr(i,205)   
      dpdr(i,394) = dpdr(i,1)*p(206) + p(1)*dpdr(i,206)   
      dpdr(i,395) = dpdr(i,1)*p(257) + p(1)*dpdr(i,257)   
      dpdr(i,396) = dpdr(i,1)*p(258) + p(1)*dpdr(i,258)   
      dpdr(i,397) = dpdr(i,3)*p(200) + p(3)*dpdr(i,200)   
      dpdr(i,398) = dpdr(i,1)*p(259) + p(1)*dpdr(i,259)   
      dpdr(i,399) = dpdr(i,3)*p(202) + p(3)*dpdr(i,202)   
      dpdr(i,400) = dpdr(i,3)*p(203) + p(3)*dpdr(i,203)   
      dpdr(i,401) = dpdr(i,1)*p(261) + p(1)*dpdr(i,261)   
      dpdr(i,402) = dpdr(i,1)*p(262) + p(1)*dpdr(i,262)   
      dpdr(i,403) = dpdr(i,1)*p(263) + p(1)*dpdr(i,263)   
      dpdr(i,404) = dpdr(i,3)*p(207) + p(3)*dpdr(i,207)   
      dpdr(i,405) = dpdr(i,1)*p(208) + p(1)*dpdr(i,208)   
      dpdr(i,406) = dpdr(i,1)*p(209) + p(1)*dpdr(i,209)   
      dpdr(i,407) = dpdr(i,1)*p(210) + p(1)*dpdr(i,210)   
      dpdr(i,408) = dpdr(i,1)*p(211) + p(1)*dpdr(i,211)   
      dpdr(i,409) = dpdr(i,1)*p(212) + p(1)*dpdr(i,212)   
      dpdr(i,410) = dpdr(i,1)*p(213) + p(1)*dpdr(i,213)   
      dpdr(i,411) = dpdr(i,1)*p(214) + p(1)*dpdr(i,214)   
      dpdr(i,412) = dpdr(i,1)*p(215) + p(1)*dpdr(i,215)   
      dpdr(i,413) = dpdr(i,1)*p(216) + p(1)*dpdr(i,216)   
      dpdr(i,414) = dpdr(i,1)*p(218) + p(1)*dpdr(i,218)   
      dpdr(i,415) = dpdr(i,1)*p(219) + p(1)*dpdr(i,219)   
      dpdr(i,416) = dpdr(i,1)*p(267) + p(1)*dpdr(i,267)   
      dpdr(i,417) = dpdr(i,1)*p(220) + p(1)*dpdr(i,220)   
      dpdr(i,418) = dpdr(i,1)*p(222) + p(1)*dpdr(i,222)   
      dpdr(i,419) = dpdr(i,5)*p(105) + p(5)*dpdr(i,105) - dpdr(i,338)  
      dpdr(i,420) = dpdr(i,1)*p(268) + p(1)*dpdr(i,268)  
      dpdr(i,421) = dpdr(i,1)*p(224) + p(1)*dpdr(i,224)  
      dpdr(i,422) = dpdr(i,5)*p(104) + p(5)*dpdr(i,104) - dpdr(i,329)  
      dpdr(i,423) = dpdr(i,6)*p(105) + p(6)*dpdr(i,105) - dpdr(i,329)  
      dpdr(i,424) = dpdr(i,1)*p(227) + p(1)*dpdr(i,227)  
      dpdr(i,425) = dpdr(i,5)*p(107) + p(5)*dpdr(i,107) - dpdr(i,341)  
      dpdr(i,426) = dpdr(i,5)*p(113) + p(5)*dpdr(i,113) - dpdr(i,345)  
      dpdr(i,427) = dpdr(i,1)*p(269) + p(1)*dpdr(i,269)  
      dpdr(i,428) = dpdr(i,7)*p(107) + p(7)*dpdr(i,107) - dpdr(i,329)  
      dpdr(i,429) = dpdr(i,1)*p(230) + p(1)*dpdr(i,230)  
      dpdr(i,430) = dpdr(i,1)*p(271) + p(1)*dpdr(i,271)  
      dpdr(i,431) = dpdr(i,1)*p(231) + p(1)*dpdr(i,231)  
      dpdr(i,432) = dpdr(i,5)*p(110) + p(5)*dpdr(i,110) - dpdr(i,338)  
      dpdr(i,433) = dpdr(i,1)*p(272) + p(1)*dpdr(i,272)  
      dpdr(i,434) = dpdr(i,5)*p(136) + p(5)*dpdr(i,136) - dpdr(i,423)  
      dpdr(i,435) = dpdr(i,1)*p(233) + p(1)*dpdr(i,233)  
      dpdr(i,436) = dpdr(i,5)*p(112) + p(5)*dpdr(i,112) - dpdr(i,341)  
      dpdr(i,437) = dpdr(i,6)*p(113) + p(6)*dpdr(i,113) - dpdr(i,341)  
      dpdr(i,438) = dpdr(i,1)*p(273) + p(1)*dpdr(i,273)  
      dpdr(i,439) = dpdr(i,5)*p(137) + p(5)*dpdr(i,137) - dpdr(i,428)  
      dpdr(i,440) = dpdr(i,6)*p(137) + p(6)*dpdr(i,137) - dpdr(i,428)  
      dpdr(i,441) = dpdr(i,1)*p(275) + p(1)*dpdr(i,275)  
      dpdr(i,442) = dpdr(i,1)*p(276) + p(1)*dpdr(i,276)  
      dpdr(i,443) = dpdr(i,5)*p(140) + p(5)*dpdr(i,140) - dpdr(i,437)  
      dpdr(i,444) = dpdr(i,1)*p(277) + p(1)*dpdr(i,277)  
      dpdr(i,445) = dpdr(i,5)*p(141) + p(5)*dpdr(i,141) - dpdr(i,440)  
      dpdr(i,446) = dpdr(i,6)*p(141) + p(6)*dpdr(i,141) - dpdr(i,439)  
      dpdr(i,447) = dpdr(i,1)*p(236) + p(1)*dpdr(i,236)   
      dpdr(i,448) = dpdr(i,1)*p(237) + p(1)*dpdr(i,237)   
      dpdr(i,449) = dpdr(i,1)*p(238) + p(1)*dpdr(i,238)   
      dpdr(i,450) = dpdr(i,3)*p(267) + p(3)*dpdr(i,267)   
      dpdr(i,451) = dpdr(i,3)*p(268) + p(3)*dpdr(i,268)   
      dpdr(i,452) = dpdr(i,3)*p(269) + p(3)*dpdr(i,269)   
      dpdr(i,453) = dpdr(i,1)*p(242) + p(1)*dpdr(i,242)   
      dpdr(i,454) = dpdr(i,3)*p(271) + p(3)*dpdr(i,271)   
      dpdr(i,455) = dpdr(i,3)*p(272) + p(3)*dpdr(i,272)   
      dpdr(i,456) = dpdr(i,3)*p(273) + p(3)*dpdr(i,273)   
      dpdr(i,457) = dpdr(i,1)*p(279) + p(1)*dpdr(i,279)   
      dpdr(i,458) = dpdr(i,3)*p(275) + p(3)*dpdr(i,275)   
      dpdr(i,459) = dpdr(i,3)*p(276) + p(3)*dpdr(i,276)   
      dpdr(i,460) = dpdr(i,3)*p(277) + p(3)*dpdr(i,277)   
      dpdr(i,461) = dpdr(i,1)*p(246) + p(1)*dpdr(i,246)   
      dpdr(i,462) = dpdr(i,1)*p(247) + p(1)*dpdr(i,247)   
      dpdr(i,463) = dpdr(i,1)*p(251) + p(1)*dpdr(i,251)   
      dpdr(i,464) = dpdr(i,3)*p(239) + p(3)*dpdr(i,239)   
      dpdr(i,465) = dpdr(i,3)*p(240) + p(3)*dpdr(i,240)   
      dpdr(i,466) = dpdr(i,3)*p(241) + p(3)*dpdr(i,241)   
      dpdr(i,467) = dpdr(i,1)*p(281) + p(1)*dpdr(i,281)   
      dpdr(i,468) = dpdr(i,3)*p(243) + p(3)*dpdr(i,243)   
      dpdr(i,469) = dpdr(i,3)*p(244) + p(3)*dpdr(i,244)   
      dpdr(i,470) = dpdr(i,3)*p(245) + p(3)*dpdr(i,245)   
      dpdr(i,471) = dpdr(i,1)*p(255) + p(1)*dpdr(i,255)   
      dpdr(i,472) = dpdr(i,1)*p(256) + p(1)*dpdr(i,256)   
      dpdr(i,473) = dpdr(i,3)*p(248) + p(3)*dpdr(i,248)   
      dpdr(i,474) = dpdr(i,3)*p(249) + p(3)*dpdr(i,249)   
      dpdr(i,475) = dpdr(i,3)*p(250) + p(3)*dpdr(i,250)   
      dpdr(i,476) = dpdr(i,1)*p(283) + p(1)*dpdr(i,283)   
      dpdr(i,477) = dpdr(i,3)*p(252) + p(3)*dpdr(i,252)   
      dpdr(i,478) = dpdr(i,3)*p(253) + p(3)*dpdr(i,253)   
      dpdr(i,479) = dpdr(i,3)*p(254) + p(3)*dpdr(i,254)   
      dpdr(i,480) = dpdr(i,1)*p(260) + p(1)*dpdr(i,260)   
      dpdr(i,481) = dpdr(i,1)*p(285) + p(1)*dpdr(i,285)   
      dpdr(i,482) = dpdr(i,3)*p(257) + p(3)*dpdr(i,257)   
      dpdr(i,483) = dpdr(i,3)*p(258) + p(3)*dpdr(i,258)   
      dpdr(i,484) = dpdr(i,3)*p(259) + p(3)*dpdr(i,259)   
      dpdr(i,485) = dpdr(i,1)*p(287) + p(1)*dpdr(i,287)   
      dpdr(i,486) = dpdr(i,3)*p(261) + p(3)*dpdr(i,261)   
      dpdr(i,487) = dpdr(i,3)*p(262) + p(3)*dpdr(i,262)   
      dpdr(i,488) = dpdr(i,3)*p(263) + p(3)*dpdr(i,263)   
      dpdr(i,489) = dpdr(i,1)*p(264) + p(1)*dpdr(i,264)   
      dpdr(i,490) = dpdr(i,1)*p(265) + p(1)*dpdr(i,265)   
      dpdr(i,491) = dpdr(i,1)*p(266) + p(1)*dpdr(i,266)   
      dpdr(i,492) = dpdr(i,1)*p(270) + p(1)*dpdr(i,270)   
      dpdr(i,493) = dpdr(i,2)*p(267) + p(2)*dpdr(i,267) - dpdr(i,422)  
      dpdr(i,494) = dpdr(i,2)*p(268) + p(2)*dpdr(i,268) - dpdr(i,423)  
      dpdr(i,495) = dpdr(i,2)*p(269) + p(2)*dpdr(i,269) - dpdr(i,428)  
      dpdr(i,496) = dpdr(i,1)*p(274) + p(1)*dpdr(i,274)  
      dpdr(i,497) = dpdr(i,5)*p(139) + p(5)*dpdr(i,139) - dpdr(i,345)  
      dpdr(i,498) = dpdr(i,6)*p(140) + p(6)*dpdr(i,140) - dpdr(i,345)  
      dpdr(i,499) = dpdr(i,7)*p(141) + p(7)*dpdr(i,141) - dpdr(i,345)  
      dpdr(i,500) = dpdr(i,1)*p(289) + p(1)*dpdr(i,289)  
      dpdr(i,501) = dpdr(i,5)*p(151) + p(5)*dpdr(i,151) - dpdr(i,446)  
      dpdr(i,502) = dpdr(i,6)*p(151) + p(6)*dpdr(i,151) - dpdr(i,445)  
      dpdr(i,503) = dpdr(i,7)*p(151) + p(7)*dpdr(i,151) - dpdr(i,443)  
      dpdr(i,504) = dpdr(i,1)*p(278) + p(1)*dpdr(i,278)   
      dpdr(i,505) = dpdr(i,3)*p(289) + p(3)*dpdr(i,289)   
      dpdr(i,506) = dpdr(i,1)*p(280) + p(1)*dpdr(i,280)   
      dpdr(i,507) = dpdr(i,3)*p(279) + p(3)*dpdr(i,279)   
      dpdr(i,508) = dpdr(i,1)*p(282) + p(1)*dpdr(i,282)   
      dpdr(i,509) = dpdr(i,3)*p(281) + p(3)*dpdr(i,281)   
      dpdr(i,510) = dpdr(i,1)*p(284) + p(1)*dpdr(i,284)   
      dpdr(i,511) = dpdr(i,3)*p(283) + p(3)*dpdr(i,283)   
      dpdr(i,512) = dpdr(i,1)*p(286) + p(1)*dpdr(i,286)   
      dpdr(i,513) = dpdr(i,3)*p(285) + p(3)*dpdr(i,285)   
      dpdr(i,514) = dpdr(i,1)*p(290) + p(1)*dpdr(i,290)   
      dpdr(i,515) = dpdr(i,3)*p(287) + p(3)*dpdr(i,287)   
      dpdr(i,516) = dpdr(i,1)*p(288) + p(1)*dpdr(i,288)   
      dpdr(i,517) = dpdr(i,2)*p(289) + p(2)*dpdr(i,289) &
                  - dpdr(i,503) - dpdr(i,502) - dpdr(i,501)  
      dpdr(i,518) = dpdr(i,3)*p(290) + p(3)*dpdr(i,290)   
      dpdr(i,519) = dpdr(i,1)*p(291) + p(1)*dpdr(i,291)   
      dpdr(i,520) = dpdr(i,1)*p(292) + p(1)*dpdr(i,292)   
      dpdr(i,521) = dpdr(i,1)*p(305) + p(1)*dpdr(i,305)   
      dpdr(i,522) = dpdr(i,1)*p(307) + p(1)*dpdr(i,307)   
      dpdr(i,523) = dpdr(i,1)*p(309) + p(1)*dpdr(i,309)   
      dpdr(i,524) = dpdr(i,1)*p(313) + p(1)*dpdr(i,313)   
      dpdr(i,525) = dpdr(i,1)*p(293) + p(1)*dpdr(i,293)   
      dpdr(i,526) = dpdr(i,1)*p(318) + p(1)*dpdr(i,318)   
      dpdr(i,527) = dpdr(i,1)*p(320) + p(1)*dpdr(i,320)   
      dpdr(i,528) = dpdr(i,1)*p(294) + p(1)*dpdr(i,294)   
      dpdr(i,529) = dpdr(i,1)*p(295) + p(1)*dpdr(i,295)   
      dpdr(i,530) = dpdr(i,1)*p(296) + p(1)*dpdr(i,296)   
      dpdr(i,531) = dpdr(i,1)*p(297) + p(1)*dpdr(i,297)   
      dpdr(i,532) = dpdr(i,1)*p(298) + p(1)*dpdr(i,298)   
      dpdr(i,533) = dpdr(i,1)*p(329) + p(1)*dpdr(i,329)   
      dpdr(i,534) = dpdr(i,1)*p(299) + p(1)*dpdr(i,299)   
      dpdr(i,535) = dpdr(i,1)*p(336) + p(1)*dpdr(i,336)   
      dpdr(i,536) = dpdr(i,1)*p(338) + p(1)*dpdr(i,338)   
      dpdr(i,537) = dpdr(i,1)*p(341) + p(1)*dpdr(i,341)   
      dpdr(i,538) = dpdr(i,1)*p(345) + p(1)*dpdr(i,345)   
      dpdr(i,539) = dpdr(i,1)*p(300) + p(1)*dpdr(i,300)   
      dpdr(i,540) = dpdr(i,1)*p(301) + p(1)*dpdr(i,301)   
      dpdr(i,541) = dpdr(i,1)*p(302) + p(1)*dpdr(i,302)   
      dpdr(i,542) = dpdr(i,1)*p(303) + p(1)*dpdr(i,303)   
      dpdr(i,543) = dpdr(i,1)*p(304) + p(1)*dpdr(i,304)   
      dpdr(i,544) = dpdr(i,1)*p(306) + p(1)*dpdr(i,306)   
      dpdr(i,545) = dpdr(i,1)*p(308) + p(1)*dpdr(i,308)   
      dpdr(i,546) = dpdr(i,1)*p(355) + p(1)*dpdr(i,355)   
      dpdr(i,547) = dpdr(i,3)*p(329) + p(3)*dpdr(i,329)   
      dpdr(i,548) = dpdr(i,1)*p(310) + p(1)*dpdr(i,310)   
      dpdr(i,549) = dpdr(i,1)*p(359) + p(1)*dpdr(i,359)   
      dpdr(i,550) = dpdr(i,1)*p(361) + p(1)*dpdr(i,361)   
      dpdr(i,551) = dpdr(i,1)*p(311) + p(1)*dpdr(i,311)   
      dpdr(i,552) = dpdr(i,1)*p(363) + p(1)*dpdr(i,363)   
      dpdr(i,553) = dpdr(i,1)*p(312) + p(1)*dpdr(i,312)   
      dpdr(i,554) = dpdr(i,3)*p(336) + p(3)*dpdr(i,336)   
      dpdr(i,555) = dpdr(i,1)*p(364) + p(1)*dpdr(i,364)   
      dpdr(i,556) = dpdr(i,3)*p(338) + p(3)*dpdr(i,338)   
      dpdr(i,557) = dpdr(i,1)*p(366) + p(1)*dpdr(i,366)   
      dpdr(i,558) = dpdr(i,1)*p(367) + p(1)*dpdr(i,367)   
      dpdr(i,559) = dpdr(i,3)*p(341) + p(3)*dpdr(i,341)   
      dpdr(i,560) = dpdr(i,1)*p(370) + p(1)*dpdr(i,370)   
      dpdr(i,561) = dpdr(i,1)*p(372) + p(1)*dpdr(i,372)   
      dpdr(i,562) = dpdr(i,1)*p(373) + p(1)*dpdr(i,373)   
      dpdr(i,563) = dpdr(i,3)*p(345) + p(3)*dpdr(i,345)   
      dpdr(i,564) = dpdr(i,1)*p(314) + p(1)*dpdr(i,314)   
      dpdr(i,565) = dpdr(i,1)*p(315) + p(1)*dpdr(i,315)   
      dpdr(i,566) = dpdr(i,1)*p(316) + p(1)*dpdr(i,316)   
      dpdr(i,567) = dpdr(i,1)*p(382) + p(1)*dpdr(i,382)   
      dpdr(i,568) = dpdr(i,1)*p(317) + p(1)*dpdr(i,317)   
      dpdr(i,569) = dpdr(i,3)*p(305) + p(3)*dpdr(i,305)   
      dpdr(i,570) = dpdr(i,1)*p(383) + p(1)*dpdr(i,383)   
      dpdr(i,571) = dpdr(i,3)*p(307) + p(3)*dpdr(i,307)   
      dpdr(i,572) = dpdr(i,1)*p(385) + p(1)*dpdr(i,385)   
      dpdr(i,573) = dpdr(i,3)*p(309) + p(3)*dpdr(i,309)   
      dpdr(i,574) = dpdr(i,1)*p(388) + p(1)*dpdr(i,388)   
      dpdr(i,575) = dpdr(i,1)*p(390) + p(1)*dpdr(i,390)   
      dpdr(i,576) = dpdr(i,1)*p(391) + p(1)*dpdr(i,391)   
      dpdr(i,577) = dpdr(i,3)*p(313) + p(3)*dpdr(i,313)   
      dpdr(i,578) = dpdr(i,1)*p(319) + p(1)*dpdr(i,319)   
      dpdr(i,579) = dpdr(i,1)*p(397) + p(1)*dpdr(i,397)   
      dpdr(i,580) = dpdr(i,1)*p(399) + p(1)*dpdr(i,399)   
      dpdr(i,581) = dpdr(i,1)*p(400) + p(1)*dpdr(i,400)   
      dpdr(i,582) = dpdr(i,3)*p(318) + p(3)*dpdr(i,318)   
      dpdr(i,583) = dpdr(i,1)*p(404) + p(1)*dpdr(i,404)   
      dpdr(i,584) = dpdr(i,3)*p(320) + p(3)*dpdr(i,320)   
      dpdr(i,585) = dpdr(i,1)*p(321) + p(1)*dpdr(i,321)   
      dpdr(i,586) = dpdr(i,1)*p(322) + p(1)*dpdr(i,322)   
      dpdr(i,587) = dpdr(i,1)*p(323) + p(1)*dpdr(i,323)   
      dpdr(i,588) = dpdr(i,1)*p(324) + p(1)*dpdr(i,324)   
      dpdr(i,589) = dpdr(i,1)*p(325) + p(1)*dpdr(i,325)   
      dpdr(i,590) = dpdr(i,1)*p(326) + p(1)*dpdr(i,326)   
      dpdr(i,591) = dpdr(i,1)*p(327) + p(1)*dpdr(i,327)   
      dpdr(i,592) = dpdr(i,1)*p(328) + p(1)*dpdr(i,328)   
      dpdr(i,593) = dpdr(i,34)*p(34) + p(34)*dpdr(i,34)   
      dpdr(i,594) = dpdr(i,1)*p(330) + p(1)*dpdr(i,330)   
      dpdr(i,595) = dpdr(i,1)*p(331) + p(1)*dpdr(i,331)   
      dpdr(i,596) = dpdr(i,1)*p(332) + p(1)*dpdr(i,332)   
      dpdr(i,597) = dpdr(i,1)*p(419) + p(1)*dpdr(i,419)   
      dpdr(i,598) = dpdr(i,1)*p(333) + p(1)*dpdr(i,333)   
      dpdr(i,599) = dpdr(i,1)*p(334) + p(1)*dpdr(i,334)   
      dpdr(i,600) = dpdr(i,1)*p(422) + p(1)*dpdr(i,422)   
      dpdr(i,601) = dpdr(i,1)*p(335) + p(1)*dpdr(i,335)   
      dpdr(i,602) = dpdr(i,34)*p(59) + p(34)*dpdr(i,59)   
      dpdr(i,603) = dpdr(i,1)*p(337) + p(1)*dpdr(i,337)   
      dpdr(i,604) = dpdr(i,34)*p(44) + p(34)*dpdr(i,44)   
      dpdr(i,605) = dpdr(i,1)*p(423) + p(1)*dpdr(i,423)   
      dpdr(i,606) = dpdr(i,34)*p(60) + p(34)*dpdr(i,60)   
      dpdr(i,607) = dpdr(i,1)*p(339) + p(1)*dpdr(i,339)   
      dpdr(i,608) = dpdr(i,1)*p(425) + p(1)*dpdr(i,425)   
      dpdr(i,609) = dpdr(i,1)*p(340) + p(1)*dpdr(i,340)   
      dpdr(i,610) = dpdr(i,34)*p(46) + p(34)*dpdr(i,46)   
      dpdr(i,611) = dpdr(i,1)*p(426) + p(1)*dpdr(i,426)   
      dpdr(i,612) = dpdr(i,34)*p(47) + p(34)*dpdr(i,47)   
      dpdr(i,613) = dpdr(i,1)*p(428) + p(1)*dpdr(i,428)   
      dpdr(i,614) = dpdr(i,34)*p(61) + p(34)*dpdr(i,61)   
      dpdr(i,615) = dpdr(i,1)*p(342) + p(1)*dpdr(i,342)   
      dpdr(i,616) = dpdr(i,1)*p(432) + p(1)*dpdr(i,432)   
      dpdr(i,617) = dpdr(i,1)*p(434) + p(1)*dpdr(i,434)   
      dpdr(i,618) = dpdr(i,1)*p(343) + p(1)*dpdr(i,343)   
      dpdr(i,619) = dpdr(i,1)*p(436) + p(1)*dpdr(i,436)   
      dpdr(i,620) = dpdr(i,1)*p(344) + p(1)*dpdr(i,344)   
      dpdr(i,621) = dpdr(i,34)*p(63) + p(34)*dpdr(i,63)   
      dpdr(i,622) = dpdr(i,1)*p(437) + p(1)*dpdr(i,437)   
      dpdr(i,623) = dpdr(i,34)*p(64) + p(34)*dpdr(i,64)   
      dpdr(i,624) = dpdr(i,1)*p(439) + p(1)*dpdr(i,439)   
      dpdr(i,625) = dpdr(i,1)*p(440) + p(1)*dpdr(i,440)   
      dpdr(i,626) = dpdr(i,34)*p(65) + p(34)*dpdr(i,65)   
      dpdr(i,627) = dpdr(i,1)*p(443) + p(1)*dpdr(i,443)   
      dpdr(i,628) = dpdr(i,1)*p(445) + p(1)*dpdr(i,445)   
      dpdr(i,629) = dpdr(i,1)*p(446) + p(1)*dpdr(i,446)   
      dpdr(i,630) = dpdr(i,34)*p(73) + p(34)*dpdr(i,73)   
      dpdr(i,631) = dpdr(i,1)*p(346) + p(1)*dpdr(i,346)   
      dpdr(i,632) = dpdr(i,1)*p(347) + p(1)*dpdr(i,347)   
      dpdr(i,633) = dpdr(i,1)*p(348) + p(1)*dpdr(i,348)   
      dpdr(i,634) = dpdr(i,1)*p(349) + p(1)*dpdr(i,349)   
      dpdr(i,635) = dpdr(i,1)*p(350) + p(1)*dpdr(i,350)   
      dpdr(i,636) = dpdr(i,1)*p(351) + p(1)*dpdr(i,351)   
      dpdr(i,637) = dpdr(i,1)*p(352) + p(1)*dpdr(i,352)   
      dpdr(i,638) = dpdr(i,1)*p(353) + p(1)*dpdr(i,353)   
      dpdr(i,639) = dpdr(i,1)*p(354) + p(1)*dpdr(i,354)   
      dpdr(i,640) = dpdr(i,1)*p(356) + p(1)*dpdr(i,356)   
      dpdr(i,641) = dpdr(i,1)*p(357) + p(1)*dpdr(i,357)   
      dpdr(i,642) = dpdr(i,1)*p(450) + p(1)*dpdr(i,450)   
      dpdr(i,643) = dpdr(i,1)*p(358) + p(1)*dpdr(i,358)   
      dpdr(i,644) = dpdr(i,1)*p(360) + p(1)*dpdr(i,360)   
      dpdr(i,645) = dpdr(i,3)*p(419) + p(3)*dpdr(i,419)   
      dpdr(i,646) = dpdr(i,1)*p(451) + p(1)*dpdr(i,451)   
      dpdr(i,647) = dpdr(i,1)*p(362) + p(1)*dpdr(i,362)   
      dpdr(i,648) = dpdr(i,3)*p(422) + p(3)*dpdr(i,422)   
      dpdr(i,649) = dpdr(i,3)*p(423) + p(3)*dpdr(i,423)   
      dpdr(i,650) = dpdr(i,1)*p(365) + p(1)*dpdr(i,365)   
      dpdr(i,651) = dpdr(i,3)*p(425) + p(3)*dpdr(i,425)   
      dpdr(i,652) = dpdr(i,3)*p(426) + p(3)*dpdr(i,426)   
      dpdr(i,653) = dpdr(i,1)*p(452) + p(1)*dpdr(i,452)   
      dpdr(i,654) = dpdr(i,3)*p(428) + p(3)*dpdr(i,428)   
      dpdr(i,655) = dpdr(i,1)*p(368) + p(1)*dpdr(i,368)   
      dpdr(i,656) = dpdr(i,1)*p(454) + p(1)*dpdr(i,454)   
      dpdr(i,657) = dpdr(i,1)*p(369) + p(1)*dpdr(i,369)   
      dpdr(i,658) = dpdr(i,3)*p(432) + p(3)*dpdr(i,432)   
      dpdr(i,659) = dpdr(i,1)*p(455) + p(1)*dpdr(i,455)   
      dpdr(i,660) = dpdr(i,3)*p(434) + p(3)*dpdr(i,434)   
      dpdr(i,661) = dpdr(i,1)*p(371) + p(1)*dpdr(i,371)   
      dpdr(i,662) = dpdr(i,3)*p(436) + p(3)*dpdr(i,436)   
      dpdr(i,663) = dpdr(i,3)*p(437) + p(3)*dpdr(i,437)   
      dpdr(i,664) = dpdr(i,1)*p(456) + p(1)*dpdr(i,456)   
      dpdr(i,665) = dpdr(i,3)*p(439) + p(3)*dpdr(i,439)   
      dpdr(i,666) = dpdr(i,3)*p(440) + p(3)*dpdr(i,440)   
      dpdr(i,667) = dpdr(i,1)*p(458) + p(1)*dpdr(i,458)   
      dpdr(i,668) = dpdr(i,1)*p(459) + p(1)*dpdr(i,459)   
      dpdr(i,669) = dpdr(i,3)*p(443) + p(3)*dpdr(i,443)   
      dpdr(i,670) = dpdr(i,1)*p(460) + p(1)*dpdr(i,460)   
      dpdr(i,671) = dpdr(i,3)*p(445) + p(3)*dpdr(i,445)   
      dpdr(i,672) = dpdr(i,3)*p(446) + p(3)*dpdr(i,446)   
      dpdr(i,673) = dpdr(i,1)*p(374) + p(1)*dpdr(i,374)   
      dpdr(i,674) = dpdr(i,1)*p(375) + p(1)*dpdr(i,375)   
      dpdr(i,675) = dpdr(i,1)*p(376) + p(1)*dpdr(i,376)   
      dpdr(i,676) = dpdr(i,1)*p(377) + p(1)*dpdr(i,377)   
      dpdr(i,677) = dpdr(i,1)*p(378) + p(1)*dpdr(i,378)   
      dpdr(i,678) = dpdr(i,1)*p(379) + p(1)*dpdr(i,379)   
      dpdr(i,679) = dpdr(i,1)*p(380) + p(1)*dpdr(i,380)   
      dpdr(i,680) = dpdr(i,1)*p(381) + p(1)*dpdr(i,381)   
      dpdr(i,681) = dpdr(i,1)*p(384) + p(1)*dpdr(i,384)   
      dpdr(i,682) = dpdr(i,3)*p(355) + p(3)*dpdr(i,355)   
      dpdr(i,683) = dpdr(i,1)*p(386) + p(1)*dpdr(i,386)   
      dpdr(i,684) = dpdr(i,1)*p(464) + p(1)*dpdr(i,464)   
      dpdr(i,685) = dpdr(i,1)*p(387) + p(1)*dpdr(i,387)   
      dpdr(i,686) = dpdr(i,3)*p(359) + p(3)*dpdr(i,359)   
      dpdr(i,687) = dpdr(i,1)*p(465) + p(1)*dpdr(i,465)   
      dpdr(i,688) = dpdr(i,3)*p(361) + p(3)*dpdr(i,361)   
      dpdr(i,689) = dpdr(i,1)*p(389) + p(1)*dpdr(i,389)   
      dpdr(i,690) = dpdr(i,3)*p(363) + p(3)*dpdr(i,363)   
      dpdr(i,691) = dpdr(i,3)*p(364) + p(3)*dpdr(i,364)   
      dpdr(i,692) = dpdr(i,1)*p(466) + p(1)*dpdr(i,466)   
      dpdr(i,693) = dpdr(i,3)*p(366) + p(3)*dpdr(i,366)   
      dpdr(i,694) = dpdr(i,3)*p(367) + p(3)*dpdr(i,367)   
      dpdr(i,695) = dpdr(i,1)*p(468) + p(1)*dpdr(i,468)   
      dpdr(i,696) = dpdr(i,1)*p(469) + p(1)*dpdr(i,469)   
      dpdr(i,697) = dpdr(i,3)*p(370) + p(3)*dpdr(i,370)   
      dpdr(i,698) = dpdr(i,1)*p(470) + p(1)*dpdr(i,470)   
      dpdr(i,699) = dpdr(i,3)*p(372) + p(3)*dpdr(i,372)   
      dpdr(i,700) = dpdr(i,3)*p(373) + p(3)*dpdr(i,373)   
      dpdr(i,701) = dpdr(i,1)*p(392) + p(1)*dpdr(i,392)   
      dpdr(i,702) = dpdr(i,1)*p(393) + p(1)*dpdr(i,393)   
      dpdr(i,703) = dpdr(i,1)*p(394) + p(1)*dpdr(i,394)   
      dpdr(i,704) = dpdr(i,1)*p(395) + p(1)*dpdr(i,395)   
      dpdr(i,705) = dpdr(i,1)*p(473) + p(1)*dpdr(i,473)   
      dpdr(i,706) = dpdr(i,1)*p(396) + p(1)*dpdr(i,396)   
      dpdr(i,707) = dpdr(i,1)*p(474) + p(1)*dpdr(i,474)   
      dpdr(i,708) = dpdr(i,1)*p(398) + p(1)*dpdr(i,398)   
      dpdr(i,709) = dpdr(i,3)*p(382) + p(3)*dpdr(i,382)   
      dpdr(i,710) = dpdr(i,3)*p(383) + p(3)*dpdr(i,383)   
      dpdr(i,711) = dpdr(i,1)*p(475) + p(1)*dpdr(i,475)   
      dpdr(i,712) = dpdr(i,3)*p(385) + p(3)*dpdr(i,385)   
      dpdr(i,713) = dpdr(i,1)*p(477) + p(1)*dpdr(i,477)   
      dpdr(i,714) = dpdr(i,1)*p(478) + p(1)*dpdr(i,478)   
      dpdr(i,715) = dpdr(i,3)*p(388) + p(3)*dpdr(i,388)   
      dpdr(i,716) = dpdr(i,1)*p(479) + p(1)*dpdr(i,479)   
      dpdr(i,717) = dpdr(i,3)*p(390) + p(3)*dpdr(i,390)   
      dpdr(i,718) = dpdr(i,3)*p(391) + p(3)*dpdr(i,391)   
      dpdr(i,719) = dpdr(i,1)*p(401) + p(1)*dpdr(i,401)   
      dpdr(i,720) = dpdr(i,1)*p(402) + p(1)*dpdr(i,402)   
      dpdr(i,721) = dpdr(i,1)*p(403) + p(1)*dpdr(i,403)   
      dpdr(i,722) = dpdr(i,1)*p(482) + p(1)*dpdr(i,482)   
      dpdr(i,723) = dpdr(i,1)*p(483) + p(1)*dpdr(i,483)   
      dpdr(i,724) = dpdr(i,3)*p(397) + p(3)*dpdr(i,397)   
      dpdr(i,725) = dpdr(i,1)*p(484) + p(1)*dpdr(i,484)   
      dpdr(i,726) = dpdr(i,3)*p(399) + p(3)*dpdr(i,399)   
      dpdr(i,727) = dpdr(i,3)*p(400) + p(3)*dpdr(i,400)   
      dpdr(i,728) = dpdr(i,1)*p(486) + p(1)*dpdr(i,486)   
      dpdr(i,729) = dpdr(i,1)*p(487) + p(1)*dpdr(i,487)   
      dpdr(i,730) = dpdr(i,1)*p(488) + p(1)*dpdr(i,488)   
      dpdr(i,731) = dpdr(i,3)*p(404) + p(3)*dpdr(i,404)   
      dpdr(i,732) = dpdr(i,1)*p(405) + p(1)*dpdr(i,405)   
      dpdr(i,733) = dpdr(i,1)*p(406) + p(1)*dpdr(i,406)   
      dpdr(i,734) = dpdr(i,1)*p(407) + p(1)*dpdr(i,407)   
      dpdr(i,735) = dpdr(i,1)*p(408) + p(1)*dpdr(i,408)   
      dpdr(i,736) = dpdr(i,1)*p(409) + p(1)*dpdr(i,409)   
      dpdr(i,737) = dpdr(i,1)*p(410) + p(1)*dpdr(i,410)   
      dpdr(i,738) = dpdr(i,1)*p(411) + p(1)*dpdr(i,411)   
      dpdr(i,739) = dpdr(i,1)*p(412) + p(1)*dpdr(i,412)   
      dpdr(i,740) = dpdr(i,1)*p(413) + p(1)*dpdr(i,413)   
      dpdr(i,741) = dpdr(i,1)*p(414) + p(1)*dpdr(i,414)   
      dpdr(i,742) = dpdr(i,1)*p(415) + p(1)*dpdr(i,415)   
      dpdr(i,743) = dpdr(i,1)*p(416) + p(1)*dpdr(i,416)   
      dpdr(i,744) = dpdr(i,1)*p(417) + p(1)*dpdr(i,417)   
      dpdr(i,745) = dpdr(i,1)*p(418) + p(1)*dpdr(i,418)   
      dpdr(i,746) = dpdr(i,1)*p(420) + p(1)*dpdr(i,420)   
      dpdr(i,747) = dpdr(i,1)*p(421) + p(1)*dpdr(i,421)   
      dpdr(i,748) = dpdr(i,1)*p(424) + p(1)*dpdr(i,424)   
      dpdr(i,749) = dpdr(i,5)*p(217) + p(5)*dpdr(i,217) - dpdr(i,612)  
      dpdr(i,750) = dpdr(i,5)*p(226) + p(5)*dpdr(i,226) - dpdr(i,623)  
      dpdr(i,751) = dpdr(i,1)*p(427) + p(1)*dpdr(i,427)  
      dpdr(i,752) = dpdr(i,5)*p(229) + p(5)*dpdr(i,229) - dpdr(i,626)  
      dpdr(i,753) = dpdr(i,1)*p(429) + p(1)*dpdr(i,429)   
      dpdr(i,754) = dpdr(i,1)*p(430) + p(1)*dpdr(i,430)   
      dpdr(i,755) = dpdr(i,1)*p(493) + p(1)*dpdr(i,493)   
      dpdr(i,756) = dpdr(i,1)*p(431) + p(1)*dpdr(i,431)   
      dpdr(i,757) = dpdr(i,6)*p(267) + p(6)*dpdr(i,267)   
      dpdr(i,758) = dpdr(i,1)*p(433) + p(1)*dpdr(i,433)   
      dpdr(i,759) = dpdr(i,59)*p(60) + p(59)*dpdr(i,60)   
      dpdr(i,760) = dpdr(i,1)*p(494) + p(1)*dpdr(i,494)   
      dpdr(i,761) = dpdr(i,5)*p(268) + p(5)*dpdr(i,268)   
      dpdr(i,762) = dpdr(i,1)*p(435) + p(1)*dpdr(i,435)   
      dpdr(i,763) = dpdr(i,7)*p(267) + p(7)*dpdr(i,267)   
      dpdr(i,764) = dpdr(i,7)*p(268) + p(7)*dpdr(i,268)   
      dpdr(i,765) = dpdr(i,1)*p(438) + p(1)*dpdr(i,438)   
      dpdr(i,766) = dpdr(i,59)*p(61) + p(59)*dpdr(i,61)   
      dpdr(i,767) = dpdr(i,60)*p(61) + p(60)*dpdr(i,61)   
      dpdr(i,768) = dpdr(i,1)*p(495) + p(1)*dpdr(i,495)   
      dpdr(i,769) = dpdr(i,5)*p(269) + p(5)*dpdr(i,269)   
      dpdr(i,770) = dpdr(i,6)*p(269) + p(6)*dpdr(i,269)   
      dpdr(i,771) = dpdr(i,1)*p(441) + p(1)*dpdr(i,441)   
      dpdr(i,772) = dpdr(i,1)*p(497) + p(1)*dpdr(i,497)   
      dpdr(i,773) = dpdr(i,1)*p(442) + p(1)*dpdr(i,442)   
      dpdr(i,774) = dpdr(i,5)*p(232) + p(5)*dpdr(i,232) - dpdr(i,623)  
      dpdr(i,775) = dpdr(i,1)*p(498) + p(1)*dpdr(i,498)  
      dpdr(i,776) = dpdr(i,5)*p(272) + p(5)*dpdr(i,272) - dpdr(i,764)  
      dpdr(i,777) = dpdr(i,1)*p(444) + p(1)*dpdr(i,444)  
      dpdr(i,778) = dpdr(i,5)*p(234) + p(5)*dpdr(i,234) - dpdr(i,626)  
      dpdr(i,779) = dpdr(i,6)*p(235) + p(6)*dpdr(i,235) - dpdr(i,626)  
      dpdr(i,780) = dpdr(i,1)*p(499) + p(1)*dpdr(i,499)  
      dpdr(i,781) = dpdr(i,5)*p(273) + p(5)*dpdr(i,273) - dpdr(i,770)  
      dpdr(i,782) = dpdr(i,6)*p(273) + p(6)*dpdr(i,273) - dpdr(i,769)  
      dpdr(i,783) = dpdr(i,1)*p(501) + p(1)*dpdr(i,501)  
      dpdr(i,784) = dpdr(i,1)*p(502) + p(1)*dpdr(i,502)  
      dpdr(i,785) = dpdr(i,5)*p(276) + p(5)*dpdr(i,276) - dpdr(i,779)  
      dpdr(i,786) = dpdr(i,1)*p(503) + p(1)*dpdr(i,503)  
      dpdr(i,787) = dpdr(i,5)*p(277) + p(5)*dpdr(i,277) - dpdr(i,782)  
      dpdr(i,788) = dpdr(i,6)*p(277) + p(6)*dpdr(i,277) - dpdr(i,781)  
      dpdr(i,789) = dpdr(i,1)*p(447) + p(1)*dpdr(i,447)   
      dpdr(i,790) = dpdr(i,1)*p(448) + p(1)*dpdr(i,448)   
      dpdr(i,791) = dpdr(i,1)*p(449) + p(1)*dpdr(i,449)   
      dpdr(i,792) = dpdr(i,1)*p(453) + p(1)*dpdr(i,453)   
      dpdr(i,793) = dpdr(i,3)*p(493) + p(3)*dpdr(i,493)   
      dpdr(i,794) = dpdr(i,3)*p(494) + p(3)*dpdr(i,494)   
      dpdr(i,795) = dpdr(i,3)*p(495) + p(3)*dpdr(i,495)   
      dpdr(i,796) = dpdr(i,1)*p(457) + p(1)*dpdr(i,457)   
      dpdr(i,797) = dpdr(i,3)*p(497) + p(3)*dpdr(i,497)   
      dpdr(i,798) = dpdr(i,3)*p(498) + p(3)*dpdr(i,498)   
      dpdr(i,799) = dpdr(i,3)*p(499) + p(3)*dpdr(i,499)   
      dpdr(i,800) = dpdr(i,1)*p(505) + p(1)*dpdr(i,505)   
      dpdr(i,801) = dpdr(i,3)*p(501) + p(3)*dpdr(i,501)   
      dpdr(i,802) = dpdr(i,3)*p(502) + p(3)*dpdr(i,502)   
      dpdr(i,803) = dpdr(i,3)*p(503) + p(3)*dpdr(i,503)   
      dpdr(i,804) = dpdr(i,1)*p(461) + p(1)*dpdr(i,461)   
      dpdr(i,805) = dpdr(i,1)*p(462) + p(1)*dpdr(i,462)   
      dpdr(i,806) = dpdr(i,1)*p(463) + p(1)*dpdr(i,463)   
      dpdr(i,807) = dpdr(i,3)*p(450) + p(3)*dpdr(i,450)   
      dpdr(i,808) = dpdr(i,3)*p(451) + p(3)*dpdr(i,451)   
      dpdr(i,809) = dpdr(i,3)*p(452) + p(3)*dpdr(i,452)   
      dpdr(i,810) = dpdr(i,1)*p(467) + p(1)*dpdr(i,467)   
      dpdr(i,811) = dpdr(i,3)*p(454) + p(3)*dpdr(i,454)   
      dpdr(i,812) = dpdr(i,3)*p(455) + p(3)*dpdr(i,455)   
      dpdr(i,813) = dpdr(i,3)*p(456) + p(3)*dpdr(i,456)   
      dpdr(i,814) = dpdr(i,1)*p(507) + p(1)*dpdr(i,507)   
      dpdr(i,815) = dpdr(i,3)*p(458) + p(3)*dpdr(i,458)   
      dpdr(i,816) = dpdr(i,3)*p(459) + p(3)*dpdr(i,459)   
      dpdr(i,817) = dpdr(i,3)*p(460) + p(3)*dpdr(i,460)   
      dpdr(i,818) = dpdr(i,1)*p(471) + p(1)*dpdr(i,471)   
      dpdr(i,819) = dpdr(i,1)*p(472) + p(1)*dpdr(i,472)   
      dpdr(i,820) = dpdr(i,1)*p(476) + p(1)*dpdr(i,476)   
      dpdr(i,821) = dpdr(i,3)*p(464) + p(3)*dpdr(i,464)   
      dpdr(i,822) = dpdr(i,3)*p(465) + p(3)*dpdr(i,465)   
      dpdr(i,823) = dpdr(i,3)*p(466) + p(3)*dpdr(i,466)   
      dpdr(i,824) = dpdr(i,1)*p(509) + p(1)*dpdr(i,509)   
      dpdr(i,825) = dpdr(i,3)*p(468) + p(3)*dpdr(i,468)   
      dpdr(i,826) = dpdr(i,3)*p(469) + p(3)*dpdr(i,469)   
      dpdr(i,827) = dpdr(i,3)*p(470) + p(3)*dpdr(i,470)   
      dpdr(i,828) = dpdr(i,1)*p(480) + p(1)*dpdr(i,480)   
      dpdr(i,829) = dpdr(i,1)*p(481) + p(1)*dpdr(i,481)   
      dpdr(i,830) = dpdr(i,3)*p(473) + p(3)*dpdr(i,473)   
      dpdr(i,831) = dpdr(i,3)*p(474) + p(3)*dpdr(i,474)   
      dpdr(i,832) = dpdr(i,3)*p(475) + p(3)*dpdr(i,475)   
      dpdr(i,833) = dpdr(i,1)*p(511) + p(1)*dpdr(i,511)   
      dpdr(i,834) = dpdr(i,3)*p(477) + p(3)*dpdr(i,477)   
      dpdr(i,835) = dpdr(i,3)*p(478) + p(3)*dpdr(i,478)   
      dpdr(i,836) = dpdr(i,3)*p(479) + p(3)*dpdr(i,479)   
      dpdr(i,837) = dpdr(i,1)*p(485) + p(1)*dpdr(i,485)   
      dpdr(i,838) = dpdr(i,1)*p(513) + p(1)*dpdr(i,513)   
      dpdr(i,839) = dpdr(i,3)*p(482) + p(3)*dpdr(i,482)   
      dpdr(i,840) = dpdr(i,3)*p(483) + p(3)*dpdr(i,483)   
      dpdr(i,841) = dpdr(i,3)*p(484) + p(3)*dpdr(i,484)   
      dpdr(i,842) = dpdr(i,1)*p(515) + p(1)*dpdr(i,515)   
      dpdr(i,843) = dpdr(i,3)*p(486) + p(3)*dpdr(i,486)   
      dpdr(i,844) = dpdr(i,3)*p(487) + p(3)*dpdr(i,487)   
      dpdr(i,845) = dpdr(i,3)*p(488) + p(3)*dpdr(i,488)   
      dpdr(i,846) = dpdr(i,1)*p(489) + p(1)*dpdr(i,489)   
      dpdr(i,847) = dpdr(i,1)*p(490) + p(1)*dpdr(i,490)   
      dpdr(i,848) = dpdr(i,1)*p(491) + p(1)*dpdr(i,491)   
      dpdr(i,849) = dpdr(i,1)*p(492) + p(1)*dpdr(i,492)   
      dpdr(i,850) = dpdr(i,5)*p(267) + p(5)*dpdr(i,267) - dpdr(i,602)  
      dpdr(i,851) = dpdr(i,6)*p(268) + p(6)*dpdr(i,268) - dpdr(i,606)  
      dpdr(i,852) = dpdr(i,7)*p(269) + p(7)*dpdr(i,269) - dpdr(i,614)  
      dpdr(i,853) = dpdr(i,1)*p(496) + p(1)*dpdr(i,496)  
      dpdr(i,854) = dpdr(i,5)*p(271) + p(5)*dpdr(i,271) - dpdr(i,621)  
      dpdr(i,855) = dpdr(i,6)*p(272) + p(6)*dpdr(i,272) - dpdr(i,623)  
      dpdr(i,856) = dpdr(i,7)*p(273) + p(7)*dpdr(i,273) - dpdr(i,626)  
      dpdr(i,857) = dpdr(i,1)*p(500) + p(1)*dpdr(i,500)  
      dpdr(i,858) = dpdr(i,5)*p(275) + p(5)*dpdr(i,275) - dpdr(i,630)  
      dpdr(i,859) = dpdr(i,6)*p(276) + p(6)*dpdr(i,276) - dpdr(i,630)  
      dpdr(i,860) = dpdr(i,7)*p(277) + p(7)*dpdr(i,277) - dpdr(i,630)  
      dpdr(i,861) = dpdr(i,1)*p(517) + p(1)*dpdr(i,517)  
      dpdr(i,862) = dpdr(i,5)*p(289) + p(5)*dpdr(i,289) - dpdr(i,788)  
      dpdr(i,863) = dpdr(i,6)*p(289) + p(6)*dpdr(i,289) - dpdr(i,787)  
      dpdr(i,864) = dpdr(i,7)*p(289) + p(7)*dpdr(i,289) - dpdr(i,785)  
      dpdr(i,865) = dpdr(i,1)*p(504) + p(1)*dpdr(i,504)   
      dpdr(i,866) = dpdr(i,3)*p(517) + p(3)*dpdr(i,517)   
      dpdr(i,867) = dpdr(i,1)*p(506) + p(1)*dpdr(i,506)   
      dpdr(i,868) = dpdr(i,3)*p(505) + p(3)*dpdr(i,505)   
      dpdr(i,869) = dpdr(i,1)*p(508) + p(1)*dpdr(i,508)   
      dpdr(i,870) = dpdr(i,3)*p(507) + p(3)*dpdr(i,507)   
      dpdr(i,871) = dpdr(i,1)*p(510) + p(1)*dpdr(i,510)   
      dpdr(i,872) = dpdr(i,3)*p(509) + p(3)*dpdr(i,509)   
      dpdr(i,873) = dpdr(i,1)*p(512) + p(1)*dpdr(i,512)   
      dpdr(i,874) = dpdr(i,3)*p(511) + p(3)*dpdr(i,511)   
      dpdr(i,875) = dpdr(i,1)*p(514) + p(1)*dpdr(i,514)   
      dpdr(i,876) = dpdr(i,3)*p(513) + p(3)*dpdr(i,513)   
      dpdr(i,877) = dpdr(i,1)*p(518) + p(1)*dpdr(i,518)   
      dpdr(i,878) = dpdr(i,3)*p(515) + p(3)*dpdr(i,515)   
      dpdr(i,879) = dpdr(i,1)*p(516) + p(1)*dpdr(i,516)   
      dpdr(i,880) = dpdr(i,2)*p(517) + p(2)*dpdr(i,517) &
                  - dpdr(i,864) - dpdr(i,863) - dpdr(i,862)  
      dpdr(i,881) = dpdr(i,3)*p(518) + p(3)*dpdr(i,518)   
      dpdr(i,882) = dpdr(i,1)*p(519) + p(1)*dpdr(i,519)   
      dpdr(i,883) = dpdr(i,1)*p(520) + p(1)*dpdr(i,520)   
      dpdr(i,884) = dpdr(i,1)*p(521) + p(1)*dpdr(i,521)   
      dpdr(i,885) = dpdr(i,1)*p(522) + p(1)*dpdr(i,522)   
      dpdr(i,886) = dpdr(i,1)*p(523) + p(1)*dpdr(i,523)   
      dpdr(i,887) = dpdr(i,1)*p(547) + p(1)*dpdr(i,547)   
      dpdr(i,888) = dpdr(i,1)*p(524) + p(1)*dpdr(i,524)   
      dpdr(i,889) = dpdr(i,1)*p(554) + p(1)*dpdr(i,554)   
      dpdr(i,890) = dpdr(i,1)*p(556) + p(1)*dpdr(i,556)   
      dpdr(i,891) = dpdr(i,1)*p(559) + p(1)*dpdr(i,559)   
      dpdr(i,892) = dpdr(i,1)*p(563) + p(1)*dpdr(i,563)   
      dpdr(i,893) = dpdr(i,1)*p(525) + p(1)*dpdr(i,525)   
      dpdr(i,894) = dpdr(i,1)*p(526) + p(1)*dpdr(i,526)   
      dpdr(i,895) = dpdr(i,1)*p(569) + p(1)*dpdr(i,569)   
      dpdr(i,896) = dpdr(i,1)*p(571) + p(1)*dpdr(i,571)   
      dpdr(i,897) = dpdr(i,1)*p(573) + p(1)*dpdr(i,573)   
      dpdr(i,898) = dpdr(i,1)*p(577) + p(1)*dpdr(i,577)   
      dpdr(i,899) = dpdr(i,1)*p(527) + p(1)*dpdr(i,527)   
      dpdr(i,900) = dpdr(i,1)*p(582) + p(1)*dpdr(i,582)   
      dpdr(i,901) = dpdr(i,1)*p(584) + p(1)*dpdr(i,584)   
      dpdr(i,902) = dpdr(i,1)*p(528) + p(1)*dpdr(i,528)   
      dpdr(i,903) = dpdr(i,1)*p(529) + p(1)*dpdr(i,529)   
      dpdr(i,904) = dpdr(i,1)*p(530) + p(1)*dpdr(i,530)   
      dpdr(i,905) = dpdr(i,1)*p(531) + p(1)*dpdr(i,531)   
      dpdr(i,906) = dpdr(i,1)*p(532) + p(1)*dpdr(i,532)   
      dpdr(i,907) = dpdr(i,1)*p(533) + p(1)*dpdr(i,533)   
      dpdr(i,908) = dpdr(i,1)*p(593) + p(1)*dpdr(i,593)   
      dpdr(i,909) = dpdr(i,1)*p(534) + p(1)*dpdr(i,534)   
      dpdr(i,910) = dpdr(i,1)*p(535) + p(1)*dpdr(i,535)   
      dpdr(i,911) = dpdr(i,1)*p(602) + p(1)*dpdr(i,602)   
      dpdr(i,912) = dpdr(i,1)*p(536) + p(1)*dpdr(i,536)   
      dpdr(i,913) = dpdr(i,1)*p(604) + p(1)*dpdr(i,604)   
      dpdr(i,914) = dpdr(i,1)*p(606) + p(1)*dpdr(i,606)   
      dpdr(i,915) = dpdr(i,1)*p(537) + p(1)*dpdr(i,537)   
      dpdr(i,916) = dpdr(i,1)*p(610) + p(1)*dpdr(i,610)   
      dpdr(i,917) = dpdr(i,1)*p(612) + p(1)*dpdr(i,612)   
      dpdr(i,918) = dpdr(i,1)*p(614) + p(1)*dpdr(i,614)   
      dpdr(i,919) = dpdr(i,1)*p(538) + p(1)*dpdr(i,538)   
      dpdr(i,920) = dpdr(i,1)*p(621) + p(1)*dpdr(i,621)   
      dpdr(i,921) = dpdr(i,1)*p(623) + p(1)*dpdr(i,623)   
      dpdr(i,922) = dpdr(i,1)*p(626) + p(1)*dpdr(i,626)   
      dpdr(i,923) = dpdr(i,1)*p(630) + p(1)*dpdr(i,630)   
      dpdr(i,924) = dpdr(i,1)*p(539) + p(1)*dpdr(i,539)   
      dpdr(i,925) = dpdr(i,1)*p(540) + p(1)*dpdr(i,540)   
      dpdr(i,926) = dpdr(i,1)*p(541) + p(1)*dpdr(i,541)   
      dpdr(i,927) = dpdr(i,1)*p(542) + p(1)*dpdr(i,542)   
      dpdr(i,928) = dpdr(i,1)*p(543) + p(1)*dpdr(i,543)   
      dpdr(i,929) = dpdr(i,1)*p(544) + p(1)*dpdr(i,544)   
      dpdr(i,930) = dpdr(i,1)*p(545) + p(1)*dpdr(i,545)   
      dpdr(i,931) = dpdr(i,1)*p(546) + p(1)*dpdr(i,546)   
      dpdr(i,932) = dpdr(i,3)*p(593) + p(3)*dpdr(i,593)   
      dpdr(i,933) = dpdr(i,1)*p(548) + p(1)*dpdr(i,548)   
      dpdr(i,934) = dpdr(i,1)*p(549) + p(1)*dpdr(i,549)   
      dpdr(i,935) = dpdr(i,1)*p(550) + p(1)*dpdr(i,550)   
      dpdr(i,936) = dpdr(i,1)*p(645) + p(1)*dpdr(i,645)   
      dpdr(i,937) = dpdr(i,1)*p(551) + p(1)*dpdr(i,551)   
      dpdr(i,938) = dpdr(i,1)*p(552) + p(1)*dpdr(i,552)   
      dpdr(i,939) = dpdr(i,1)*p(648) + p(1)*dpdr(i,648)   
      dpdr(i,940) = dpdr(i,1)*p(553) + p(1)*dpdr(i,553)   
      dpdr(i,941) = dpdr(i,3)*p(602) + p(3)*dpdr(i,602)   
      dpdr(i,942) = dpdr(i,1)*p(555) + p(1)*dpdr(i,555)   
      dpdr(i,943) = dpdr(i,3)*p(604) + p(3)*dpdr(i,604)   
      dpdr(i,944) = dpdr(i,1)*p(649) + p(1)*dpdr(i,649)   
      dpdr(i,945) = dpdr(i,3)*p(606) + p(3)*dpdr(i,606)   
      dpdr(i,946) = dpdr(i,1)*p(557) + p(1)*dpdr(i,557)   
      dpdr(i,947) = dpdr(i,1)*p(651) + p(1)*dpdr(i,651)   
      dpdr(i,948) = dpdr(i,1)*p(558) + p(1)*dpdr(i,558)   
      dpdr(i,949) = dpdr(i,3)*p(610) + p(3)*dpdr(i,610)   
      dpdr(i,950) = dpdr(i,1)*p(652) + p(1)*dpdr(i,652)   
      dpdr(i,951) = dpdr(i,3)*p(612) + p(3)*dpdr(i,612)   
      dpdr(i,952) = dpdr(i,1)*p(654) + p(1)*dpdr(i,654)   
      dpdr(i,953) = dpdr(i,3)*p(614) + p(3)*dpdr(i,614)   
      dpdr(i,954) = dpdr(i,1)*p(560) + p(1)*dpdr(i,560)   
      dpdr(i,955) = dpdr(i,1)*p(658) + p(1)*dpdr(i,658)   
      dpdr(i,956) = dpdr(i,1)*p(660) + p(1)*dpdr(i,660)   
      dpdr(i,957) = dpdr(i,1)*p(561) + p(1)*dpdr(i,561)   
      dpdr(i,958) = dpdr(i,1)*p(662) + p(1)*dpdr(i,662)   
      dpdr(i,959) = dpdr(i,1)*p(562) + p(1)*dpdr(i,562)   
      dpdr(i,960) = dpdr(i,3)*p(621) + p(3)*dpdr(i,621)   
      dpdr(i,961) = dpdr(i,1)*p(663) + p(1)*dpdr(i,663)   
      dpdr(i,962) = dpdr(i,3)*p(623) + p(3)*dpdr(i,623)   
      dpdr(i,963) = dpdr(i,1)*p(665) + p(1)*dpdr(i,665)   
      dpdr(i,964) = dpdr(i,1)*p(666) + p(1)*dpdr(i,666)   
      dpdr(i,965) = dpdr(i,3)*p(626) + p(3)*dpdr(i,626)   
      dpdr(i,966) = dpdr(i,1)*p(669) + p(1)*dpdr(i,669)   
      dpdr(i,967) = dpdr(i,1)*p(671) + p(1)*dpdr(i,671)   
      dpdr(i,968) = dpdr(i,1)*p(672) + p(1)*dpdr(i,672)   
      dpdr(i,969) = dpdr(i,3)*p(630) + p(3)*dpdr(i,630)   
      dpdr(i,970) = dpdr(i,1)*p(564) + p(1)*dpdr(i,564)   
      dpdr(i,971) = dpdr(i,1)*p(565) + p(1)*dpdr(i,565)   
      dpdr(i,972) = dpdr(i,1)*p(566) + p(1)*dpdr(i,566)   
      dpdr(i,973) = dpdr(i,1)*p(567) + p(1)*dpdr(i,567)   
      dpdr(i,974) = dpdr(i,1)*p(568) + p(1)*dpdr(i,568)   
      dpdr(i,975) = dpdr(i,1)*p(570) + p(1)*dpdr(i,570)   
      dpdr(i,976) = dpdr(i,1)*p(572) + p(1)*dpdr(i,572)   
      dpdr(i,977) = dpdr(i,1)*p(682) + p(1)*dpdr(i,682)   
      dpdr(i,978) = dpdr(i,3)*p(547) + p(3)*dpdr(i,547)   
      dpdr(i,979) = dpdr(i,1)*p(574) + p(1)*dpdr(i,574)   
      dpdr(i,980) = dpdr(i,1)*p(686) + p(1)*dpdr(i,686)   
      dpdr(i,981) = dpdr(i,1)*p(688) + p(1)*dpdr(i,688)   
      dpdr(i,982) = dpdr(i,1)*p(575) + p(1)*dpdr(i,575)   
      dpdr(i,983) = dpdr(i,1)*p(690) + p(1)*dpdr(i,690)   
      dpdr(i,984) = dpdr(i,1)*p(576) + p(1)*dpdr(i,576)   
      dpdr(i,985) = dpdr(i,3)*p(554) + p(3)*dpdr(i,554)   
      dpdr(i,986) = dpdr(i,1)*p(691) + p(1)*dpdr(i,691)   
      dpdr(i,987) = dpdr(i,3)*p(556) + p(3)*dpdr(i,556)   
      dpdr(i,988) = dpdr(i,1)*p(693) + p(1)*dpdr(i,693)   
      dpdr(i,989) = dpdr(i,1)*p(694) + p(1)*dpdr(i,694)   
      dpdr(i,990) = dpdr(i,3)*p(559) + p(3)*dpdr(i,559)   
      dpdr(i,991) = dpdr(i,1)*p(697) + p(1)*dpdr(i,697)   
      dpdr(i,992) = dpdr(i,1)*p(699) + p(1)*dpdr(i,699)   
      dpdr(i,993) = dpdr(i,1)*p(700) + p(1)*dpdr(i,700)   
      dpdr(i,994) = dpdr(i,3)*p(563) + p(3)*dpdr(i,563)   
      dpdr(i,995) = dpdr(i,1)*p(578) + p(1)*dpdr(i,578)   
      dpdr(i,996) = dpdr(i,1)*p(579) + p(1)*dpdr(i,579)   
      dpdr(i,997) = dpdr(i,1)*p(580) + p(1)*dpdr(i,580)   
      dpdr(i,998) = dpdr(i,1)*p(709) + p(1)*dpdr(i,709)   
      dpdr(i,999) = dpdr(i,1)*p(581) + p(1)*dpdr(i,581)   
      dpdr(i,1000) = dpdr(i,3)*p(569) + p(3)*dpdr(i,569)   
      dpdr(i,1001) = dpdr(i,1)*p(710) + p(1)*dpdr(i,710)   
      dpdr(i,1002) = dpdr(i,3)*p(571) + p(3)*dpdr(i,571)   
      dpdr(i,1003) = dpdr(i,1)*p(712) + p(1)*dpdr(i,712)   
      dpdr(i,1004) = dpdr(i,3)*p(573) + p(3)*dpdr(i,573)   
      dpdr(i,1005) = dpdr(i,1)*p(715) + p(1)*dpdr(i,715)   
      dpdr(i,1006) = dpdr(i,1)*p(717) + p(1)*dpdr(i,717)   
      dpdr(i,1007) = dpdr(i,1)*p(718) + p(1)*dpdr(i,718)   
      dpdr(i,1008) = dpdr(i,3)*p(577) + p(3)*dpdr(i,577)   
      dpdr(i,1009) = dpdr(i,1)*p(583) + p(1)*dpdr(i,583)   
      dpdr(i,1010) = dpdr(i,1)*p(724) + p(1)*dpdr(i,724)   
      dpdr(i,1011) = dpdr(i,1)*p(726) + p(1)*dpdr(i,726)   
      dpdr(i,1012) = dpdr(i,1)*p(727) + p(1)*dpdr(i,727)   
      dpdr(i,1013) = dpdr(i,3)*p(582) + p(3)*dpdr(i,582)   
      dpdr(i,1014) = dpdr(i,1)*p(731) + p(1)*dpdr(i,731)   
      dpdr(i,1015) = dpdr(i,3)*p(584) + p(3)*dpdr(i,584)   
      dpdr(i,1016) = dpdr(i,1)*p(585) + p(1)*dpdr(i,585)   
      dpdr(i,1017) = dpdr(i,1)*p(586) + p(1)*dpdr(i,586)   
      dpdr(i,1018) = dpdr(i,1)*p(587) + p(1)*dpdr(i,587)   
      dpdr(i,1019) = dpdr(i,1)*p(588) + p(1)*dpdr(i,588)   
      dpdr(i,1020) = dpdr(i,1)*p(589) + p(1)*dpdr(i,589)   
      dpdr(i,1021) = dpdr(i,1)*p(590) + p(1)*dpdr(i,590)   
      dpdr(i,1022) = dpdr(i,1)*p(591) + p(1)*dpdr(i,591)   
      dpdr(i,1023) = dpdr(i,1)*p(592) + p(1)*dpdr(i,592)   
      dpdr(i,1024) = dpdr(i,1)*p(594) + p(1)*dpdr(i,594)   
      dpdr(i,1025) = dpdr(i,1)*p(595) + p(1)*dpdr(i,595)   
      dpdr(i,1026) = dpdr(i,1)*p(596) + p(1)*dpdr(i,596)   
      dpdr(i,1027) = dpdr(i,1)*p(597) + p(1)*dpdr(i,597)   
      dpdr(i,1028) = dpdr(i,1)*p(598) + p(1)*dpdr(i,598)   
      dpdr(i,1029) = dpdr(i,1)*p(599) + p(1)*dpdr(i,599)   
      dpdr(i,1030) = dpdr(i,1)*p(600) + p(1)*dpdr(i,600)   
      dpdr(i,1031) = dpdr(i,1)*p(601) + p(1)*dpdr(i,601)   
      dpdr(i,1032) = dpdr(i,1)*p(603) + p(1)*dpdr(i,603)   
      dpdr(i,1033) = dpdr(i,1)*p(605) + p(1)*dpdr(i,605)   
      dpdr(i,1034) = dpdr(i,1)*p(607) + p(1)*dpdr(i,607)   
      dpdr(i,1035) = dpdr(i,1)*p(608) + p(1)*dpdr(i,608)   
      dpdr(i,1036) = dpdr(i,1)*p(749) + p(1)*dpdr(i,749)   
      dpdr(i,1037) = dpdr(i,1)*p(609) + p(1)*dpdr(i,609)   
      dpdr(i,1038) = dpdr(i,34)*p(104) + p(34)*dpdr(i,104)  
      dpdr(i,1039) = dpdr(i,1)*p(611) + p(1)*dpdr(i,611)   
      dpdr(i,1040) = dpdr(i,34)*p(82) + p(34)*dpdr(i,82)   
      dpdr(i,1041) = dpdr(i,1)*p(750) + p(1)*dpdr(i,750)   
      dpdr(i,1042) = dpdr(i,34)*p(105) + p(34)*dpdr(i,105)  
      dpdr(i,1043) = dpdr(i,1)*p(613) + p(1)*dpdr(i,613)   
      dpdr(i,1044) = dpdr(i,1)*p(752) + p(1)*dpdr(i,752)   
      dpdr(i,1045) = dpdr(i,34)*p(107) + p(34)*dpdr(i,107)  
      dpdr(i,1046) = dpdr(i,1)*p(615) + p(1)*dpdr(i,615)   
      dpdr(i,1047) = dpdr(i,1)*p(616) + p(1)*dpdr(i,616)   
      dpdr(i,1048) = dpdr(i,1)*p(757) + p(1)*dpdr(i,757)   
      dpdr(i,1049) = dpdr(i,1)*p(617) + p(1)*dpdr(i,617)   
      dpdr(i,1050) = dpdr(i,1)*p(759) + p(1)*dpdr(i,759)   
      dpdr(i,1051) = dpdr(i,1)*p(761) + p(1)*dpdr(i,761)   
      dpdr(i,1052) = dpdr(i,1)*p(618) + p(1)*dpdr(i,618)   
      dpdr(i,1053) = dpdr(i,1)*p(619) + p(1)*dpdr(i,619)   
      dpdr(i,1054) = dpdr(i,1)*p(763) + p(1)*dpdr(i,763)   
      dpdr(i,1055) = dpdr(i,1)*p(620) + p(1)*dpdr(i,620)   
      dpdr(i,1056) = dpdr(i,34)*p(135) + p(34)*dpdr(i,135)  
      dpdr(i,1057) = dpdr(i,1)*p(622) + p(1)*dpdr(i,622)   
      dpdr(i,1058) = dpdr(i,34)*p(110) + p(34)*dpdr(i,110)  
      dpdr(i,1059) = dpdr(i,1)*p(764) + p(1)*dpdr(i,764)   
      dpdr(i,1060) = dpdr(i,34)*p(136) + p(34)*dpdr(i,136)  
      dpdr(i,1061) = dpdr(i,1)*p(624) + p(1)*dpdr(i,624)   
      dpdr(i,1062) = dpdr(i,1)*p(766) + p(1)*dpdr(i,766)   
      dpdr(i,1063) = dpdr(i,1)*p(625) + p(1)*dpdr(i,625)   
      dpdr(i,1064) = dpdr(i,34)*p(112) + p(34)*dpdr(i,112)  
      dpdr(i,1065) = dpdr(i,1)*p(767) + p(1)*dpdr(i,767)   
      dpdr(i,1066) = dpdr(i,34)*p(113) + p(34)*dpdr(i,113)  
      dpdr(i,1067) = dpdr(i,1)*p(769) + p(1)*dpdr(i,769)   
      dpdr(i,1068) = dpdr(i,1)*p(770) + p(1)*dpdr(i,770)   
      dpdr(i,1069) = dpdr(i,34)*p(137) + p(34)*dpdr(i,137)  
      dpdr(i,1070) = dpdr(i,1)*p(627) + p(1)*dpdr(i,627)   
      dpdr(i,1071) = dpdr(i,1)*p(774) + p(1)*dpdr(i,774)   
      dpdr(i,1072) = dpdr(i,1)*p(776) + p(1)*dpdr(i,776)   
      dpdr(i,1073) = dpdr(i,1)*p(628) + p(1)*dpdr(i,628)   
      dpdr(i,1074) = dpdr(i,1)*p(778) + p(1)*dpdr(i,778)   
      dpdr(i,1075) = dpdr(i,1)*p(629) + p(1)*dpdr(i,629)   
      dpdr(i,1076) = dpdr(i,34)*p(139) + p(34)*dpdr(i,139)  
      dpdr(i,1077) = dpdr(i,1)*p(779) + p(1)*dpdr(i,779)   
      dpdr(i,1078) = dpdr(i,34)*p(140) + p(34)*dpdr(i,140)  
      dpdr(i,1079) = dpdr(i,1)*p(781) + p(1)*dpdr(i,781)   
      dpdr(i,1080) = dpdr(i,1)*p(782) + p(1)*dpdr(i,782)   
      dpdr(i,1081) = dpdr(i,34)*p(141) + p(34)*dpdr(i,141)  
      dpdr(i,1082) = dpdr(i,1)*p(785) + p(1)*dpdr(i,785)   
      dpdr(i,1083) = dpdr(i,1)*p(787) + p(1)*dpdr(i,787)   
      dpdr(i,1084) = dpdr(i,1)*p(788) + p(1)*dpdr(i,788)   
      dpdr(i,1085) = dpdr(i,34)*p(151) + p(34)*dpdr(i,151)  
      dpdr(i,1086) = dpdr(i,1)*p(631) + p(1)*dpdr(i,631)   
      dpdr(i,1087) = dpdr(i,1)*p(632) + p(1)*dpdr(i,632)   
      dpdr(i,1088) = dpdr(i,1)*p(633) + p(1)*dpdr(i,633)   
      dpdr(i,1089) = dpdr(i,1)*p(634) + p(1)*dpdr(i,634)   
      dpdr(i,1090) = dpdr(i,1)*p(635) + p(1)*dpdr(i,635)   
      dpdr(i,1091) = dpdr(i,1)*p(636) + p(1)*dpdr(i,636)   
      dpdr(i,1092) = dpdr(i,1)*p(637) + p(1)*dpdr(i,637)   
      dpdr(i,1093) = dpdr(i,1)*p(638) + p(1)*dpdr(i,638)   
      dpdr(i,1094) = dpdr(i,1)*p(639) + p(1)*dpdr(i,639)   
      dpdr(i,1095) = dpdr(i,1)*p(640) + p(1)*dpdr(i,640)   
      dpdr(i,1096) = dpdr(i,1)*p(641) + p(1)*dpdr(i,641)   
      dpdr(i,1097) = dpdr(i,1)*p(642) + p(1)*dpdr(i,642)   
      dpdr(i,1098) = dpdr(i,1)*p(643) + p(1)*dpdr(i,643)   
      dpdr(i,1099) = dpdr(i,1)*p(644) + p(1)*dpdr(i,644)   
      dpdr(i,1100) = dpdr(i,1)*p(646) + p(1)*dpdr(i,646)   
      dpdr(i,1101) = dpdr(i,1)*p(647) + p(1)*dpdr(i,647)   
      dpdr(i,1102) = dpdr(i,1)*p(650) + p(1)*dpdr(i,650)   
      dpdr(i,1103) = dpdr(i,3)*p(749) + p(3)*dpdr(i,749)   
      dpdr(i,1104) = dpdr(i,3)*p(750) + p(3)*dpdr(i,750)   
      dpdr(i,1105) = dpdr(i,1)*p(653) + p(1)*dpdr(i,653)   
      dpdr(i,1106) = dpdr(i,3)*p(752) + p(3)*dpdr(i,752)   
      dpdr(i,1107) = dpdr(i,1)*p(655) + p(1)*dpdr(i,655)   
      dpdr(i,1108) = dpdr(i,1)*p(656) + p(1)*dpdr(i,656)   
      dpdr(i,1109) = dpdr(i,1)*p(793) + p(1)*dpdr(i,793)   
      dpdr(i,1110) = dpdr(i,1)*p(657) + p(1)*dpdr(i,657)   
      dpdr(i,1111) = dpdr(i,3)*p(757) + p(3)*dpdr(i,757)   
      dpdr(i,1112) = dpdr(i,1)*p(659) + p(1)*dpdr(i,659)   
      dpdr(i,1113) = dpdr(i,3)*p(759) + p(3)*dpdr(i,759)   
      dpdr(i,1114) = dpdr(i,1)*p(794) + p(1)*dpdr(i,794)   
      dpdr(i,1115) = dpdr(i,3)*p(761) + p(3)*dpdr(i,761)   
      dpdr(i,1116) = dpdr(i,1)*p(661) + p(1)*dpdr(i,661)   
      dpdr(i,1117) = dpdr(i,3)*p(763) + p(3)*dpdr(i,763)   
      dpdr(i,1118) = dpdr(i,3)*p(764) + p(3)*dpdr(i,764)   
      dpdr(i,1119) = dpdr(i,1)*p(664) + p(1)*dpdr(i,664)   
      dpdr(i,1120) = dpdr(i,3)*p(766) + p(3)*dpdr(i,766)   
      dpdr(i,1121) = dpdr(i,3)*p(767) + p(3)*dpdr(i,767)   
      dpdr(i,1122) = dpdr(i,1)*p(795) + p(1)*dpdr(i,795)   
      dpdr(i,1123) = dpdr(i,3)*p(769) + p(3)*dpdr(i,769)   
      dpdr(i,1124) = dpdr(i,3)*p(770) + p(3)*dpdr(i,770)   
      dpdr(i,1125) = dpdr(i,1)*p(667) + p(1)*dpdr(i,667)   
      dpdr(i,1126) = dpdr(i,1)*p(797) + p(1)*dpdr(i,797)   
      dpdr(i,1127) = dpdr(i,1)*p(668) + p(1)*dpdr(i,668)   
      dpdr(i,1128) = dpdr(i,3)*p(774) + p(3)*dpdr(i,774)   
      dpdr(i,1129) = dpdr(i,1)*p(798) + p(1)*dpdr(i,798)   
      dpdr(i,1130) = dpdr(i,3)*p(776) + p(3)*dpdr(i,776)   
      dpdr(i,1131) = dpdr(i,1)*p(670) + p(1)*dpdr(i,670)   
      dpdr(i,1132) = dpdr(i,3)*p(778) + p(3)*dpdr(i,778)   
      dpdr(i,1133) = dpdr(i,3)*p(779) + p(3)*dpdr(i,779)   
      dpdr(i,1134) = dpdr(i,1)*p(799) + p(1)*dpdr(i,799)   
      dpdr(i,1135) = dpdr(i,3)*p(781) + p(3)*dpdr(i,781)   
      dpdr(i,1136) = dpdr(i,3)*p(782) + p(3)*dpdr(i,782)   
      dpdr(i,1137) = dpdr(i,1)*p(801) + p(1)*dpdr(i,801)   
      dpdr(i,1138) = dpdr(i,1)*p(802) + p(1)*dpdr(i,802)   
      dpdr(i,1139) = dpdr(i,3)*p(785) + p(3)*dpdr(i,785)   
      dpdr(i,1140) = dpdr(i,1)*p(803) + p(1)*dpdr(i,803)   
      dpdr(i,1141) = dpdr(i,3)*p(787) + p(3)*dpdr(i,787)   
      dpdr(i,1142) = dpdr(i,3)*p(788) + p(3)*dpdr(i,788)   
      dpdr(i,1143) = dpdr(i,1)*p(673) + p(1)*dpdr(i,673)   
      dpdr(i,1144) = dpdr(i,1)*p(674) + p(1)*dpdr(i,674)   
      dpdr(i,1145) = dpdr(i,1)*p(675) + p(1)*dpdr(i,675)   
      dpdr(i,1146) = dpdr(i,1)*p(676) + p(1)*dpdr(i,676)   
      dpdr(i,1147) = dpdr(i,1)*p(677) + p(1)*dpdr(i,677)   
      dpdr(i,1148) = dpdr(i,1)*p(678) + p(1)*dpdr(i,678)   
      dpdr(i,1149) = dpdr(i,1)*p(679) + p(1)*dpdr(i,679)   
      dpdr(i,1150) = dpdr(i,1)*p(680) + p(1)*dpdr(i,680)   
      dpdr(i,1151) = dpdr(i,1)*p(681) + p(1)*dpdr(i,681)   
      dpdr(i,1152) = dpdr(i,1)*p(683) + p(1)*dpdr(i,683)   
      dpdr(i,1153) = dpdr(i,1)*p(684) + p(1)*dpdr(i,684)   
      dpdr(i,1154) = dpdr(i,1)*p(807) + p(1)*dpdr(i,807)   
      dpdr(i,1155) = dpdr(i,1)*p(685) + p(1)*dpdr(i,685)   
      dpdr(i,1156) = dpdr(i,1)*p(687) + p(1)*dpdr(i,687)   
      dpdr(i,1157) = dpdr(i,3)*p(645) + p(3)*dpdr(i,645)   
      dpdr(i,1158) = dpdr(i,1)*p(808) + p(1)*dpdr(i,808)   
      dpdr(i,1159) = dpdr(i,1)*p(689) + p(1)*dpdr(i,689)   
      dpdr(i,1160) = dpdr(i,3)*p(648) + p(3)*dpdr(i,648)   
      dpdr(i,1161) = dpdr(i,3)*p(649) + p(3)*dpdr(i,649)   
      dpdr(i,1162) = dpdr(i,1)*p(692) + p(1)*dpdr(i,692)   
      dpdr(i,1163) = dpdr(i,3)*p(651) + p(3)*dpdr(i,651)   
      dpdr(i,1164) = dpdr(i,3)*p(652) + p(3)*dpdr(i,652)   
      dpdr(i,1165) = dpdr(i,1)*p(809) + p(1)*dpdr(i,809)   
      dpdr(i,1166) = dpdr(i,3)*p(654) + p(3)*dpdr(i,654)   
      dpdr(i,1167) = dpdr(i,1)*p(695) + p(1)*dpdr(i,695)   
      dpdr(i,1168) = dpdr(i,1)*p(811) + p(1)*dpdr(i,811)   
      dpdr(i,1169) = dpdr(i,1)*p(696) + p(1)*dpdr(i,696)   
      dpdr(i,1170) = dpdr(i,3)*p(658) + p(3)*dpdr(i,658)   
      dpdr(i,1171) = dpdr(i,1)*p(812) + p(1)*dpdr(i,812)   
      dpdr(i,1172) = dpdr(i,3)*p(660) + p(3)*dpdr(i,660)   
      dpdr(i,1173) = dpdr(i,1)*p(698) + p(1)*dpdr(i,698)   
      dpdr(i,1174) = dpdr(i,3)*p(662) + p(3)*dpdr(i,662)   
      dpdr(i,1175) = dpdr(i,3)*p(663) + p(3)*dpdr(i,663)   
      dpdr(i,1176) = dpdr(i,1)*p(813) + p(1)*dpdr(i,813)   
      dpdr(i,1177) = dpdr(i,3)*p(665) + p(3)*dpdr(i,665)   
      dpdr(i,1178) = dpdr(i,3)*p(666) + p(3)*dpdr(i,666)   
      dpdr(i,1179) = dpdr(i,1)*p(815) + p(1)*dpdr(i,815)   
      dpdr(i,1180) = dpdr(i,1)*p(816) + p(1)*dpdr(i,816)   
      dpdr(i,1181) = dpdr(i,3)*p(669) + p(3)*dpdr(i,669)   
      dpdr(i,1182) = dpdr(i,1)*p(817) + p(1)*dpdr(i,817)   
      dpdr(i,1183) = dpdr(i,3)*p(671) + p(3)*dpdr(i,671)   
      dpdr(i,1184) = dpdr(i,3)*p(672) + p(3)*dpdr(i,672)   
      dpdr(i,1185) = dpdr(i,1)*p(701) + p(1)*dpdr(i,701)   
      dpdr(i,1186) = dpdr(i,1)*p(702) + p(1)*dpdr(i,702)   
      dpdr(i,1187) = dpdr(i,1)*p(703) + p(1)*dpdr(i,703)   
      dpdr(i,1188) = dpdr(i,1)*p(704) + p(1)*dpdr(i,704)   
      dpdr(i,1189) = dpdr(i,1)*p(705) + p(1)*dpdr(i,705)   
      dpdr(i,1190) = dpdr(i,1)*p(706) + p(1)*dpdr(i,706)   
      dpdr(i,1191) = dpdr(i,1)*p(707) + p(1)*dpdr(i,707)   
      dpdr(i,1192) = dpdr(i,1)*p(708) + p(1)*dpdr(i,708)   
      dpdr(i,1193) = dpdr(i,1)*p(711) + p(1)*dpdr(i,711)   
      dpdr(i,1194) = dpdr(i,3)*p(682) + p(3)*dpdr(i,682)   
      dpdr(i,1195) = dpdr(i,1)*p(713) + p(1)*dpdr(i,713)   
      dpdr(i,1196) = dpdr(i,1)*p(821) + p(1)*dpdr(i,821)   
      dpdr(i,1197) = dpdr(i,1)*p(714) + p(1)*dpdr(i,714)   
      dpdr(i,1198) = dpdr(i,3)*p(686) + p(3)*dpdr(i,686)   
      dpdr(i,1199) = dpdr(i,1)*p(822) + p(1)*dpdr(i,822)   
      dpdr(i,1200) = dpdr(i,3)*p(688) + p(3)*dpdr(i,688)   
      dpdr(i,1201) = dpdr(i,1)*p(716) + p(1)*dpdr(i,716)   
      dpdr(i,1202) = dpdr(i,3)*p(690) + p(3)*dpdr(i,690)   
      dpdr(i,1203) = dpdr(i,3)*p(691) + p(3)*dpdr(i,691)   
      dpdr(i,1204) = dpdr(i,1)*p(823) + p(1)*dpdr(i,823)   
      dpdr(i,1205) = dpdr(i,3)*p(693) + p(3)*dpdr(i,693)   
      dpdr(i,1206) = dpdr(i,3)*p(694) + p(3)*dpdr(i,694)   
      dpdr(i,1207) = dpdr(i,1)*p(825) + p(1)*dpdr(i,825)   
      dpdr(i,1208) = dpdr(i,1)*p(826) + p(1)*dpdr(i,826)   
      dpdr(i,1209) = dpdr(i,3)*p(697) + p(3)*dpdr(i,697)   
      dpdr(i,1210) = dpdr(i,1)*p(827) + p(1)*dpdr(i,827)   
      dpdr(i,1211) = dpdr(i,3)*p(699) + p(3)*dpdr(i,699)   
      dpdr(i,1212) = dpdr(i,3)*p(700) + p(3)*dpdr(i,700)   
      dpdr(i,1213) = dpdr(i,1)*p(719) + p(1)*dpdr(i,719)   
      dpdr(i,1214) = dpdr(i,1)*p(720) + p(1)*dpdr(i,720)   
      dpdr(i,1215) = dpdr(i,1)*p(721) + p(1)*dpdr(i,721)   
      dpdr(i,1216) = dpdr(i,1)*p(722) + p(1)*dpdr(i,722)   
      dpdr(i,1217) = dpdr(i,1)*p(830) + p(1)*dpdr(i,830)   
      dpdr(i,1218) = dpdr(i,1)*p(723) + p(1)*dpdr(i,723)   
      dpdr(i,1219) = dpdr(i,1)*p(831) + p(1)*dpdr(i,831)   
      dpdr(i,1220) = dpdr(i,1)*p(725) + p(1)*dpdr(i,725)   
      dpdr(i,1221) = dpdr(i,3)*p(709) + p(3)*dpdr(i,709)   
      dpdr(i,1222) = dpdr(i,3)*p(710) + p(3)*dpdr(i,710)   
      dpdr(i,1223) = dpdr(i,1)*p(832) + p(1)*dpdr(i,832)   
      dpdr(i,1224) = dpdr(i,3)*p(712) + p(3)*dpdr(i,712)   
      dpdr(i,1225) = dpdr(i,1)*p(834) + p(1)*dpdr(i,834)   
      dpdr(i,1226) = dpdr(i,1)*p(835) + p(1)*dpdr(i,835)   
      dpdr(i,1227) = dpdr(i,3)*p(715) + p(3)*dpdr(i,715)   
      dpdr(i,1228) = dpdr(i,1)*p(836) + p(1)*dpdr(i,836)   
      dpdr(i,1229) = dpdr(i,3)*p(717) + p(3)*dpdr(i,717)   
      dpdr(i,1230) = dpdr(i,3)*p(718) + p(3)*dpdr(i,718)   
      dpdr(i,1231) = dpdr(i,1)*p(728) + p(1)*dpdr(i,728)   
      dpdr(i,1232) = dpdr(i,1)*p(729) + p(1)*dpdr(i,729)   
      dpdr(i,1233) = dpdr(i,1)*p(730) + p(1)*dpdr(i,730)   
      dpdr(i,1234) = dpdr(i,1)*p(839) + p(1)*dpdr(i,839)   
      dpdr(i,1235) = dpdr(i,1)*p(840) + p(1)*dpdr(i,840)   
      dpdr(i,1236) = dpdr(i,3)*p(724) + p(3)*dpdr(i,724)   
      dpdr(i,1237) = dpdr(i,1)*p(841) + p(1)*dpdr(i,841)   
      dpdr(i,1238) = dpdr(i,3)*p(726) + p(3)*dpdr(i,726)   
      dpdr(i,1239) = dpdr(i,3)*p(727) + p(3)*dpdr(i,727)   
      dpdr(i,1240) = dpdr(i,1)*p(843) + p(1)*dpdr(i,843)   
      dpdr(i,1241) = dpdr(i,1)*p(844) + p(1)*dpdr(i,844)   
      dpdr(i,1242) = dpdr(i,1)*p(845) + p(1)*dpdr(i,845)   
      dpdr(i,1243) = dpdr(i,3)*p(731) + p(3)*dpdr(i,731)   
      dpdr(i,1244) = dpdr(i,1)*p(732) + p(1)*dpdr(i,732)   
      dpdr(i,1245) = dpdr(i,1)*p(733) + p(1)*dpdr(i,733)   
      dpdr(i,1246) = dpdr(i,1)*p(734) + p(1)*dpdr(i,734)   
      dpdr(i,1247) = dpdr(i,1)*p(735) + p(1)*dpdr(i,735)   
      dpdr(i,1248) = dpdr(i,1)*p(736) + p(1)*dpdr(i,736)   
      dpdr(i,1249) = dpdr(i,1)*p(737) + p(1)*dpdr(i,737)   
      dpdr(i,1250) = dpdr(i,1)*p(738) + p(1)*dpdr(i,738)   
      dpdr(i,1251) = dpdr(i,1)*p(739) + p(1)*dpdr(i,739)   
      dpdr(i,1252) = dpdr(i,1)*p(740) + p(1)*dpdr(i,740)   
      dpdr(i,1253) = dpdr(i,1)*p(741) + p(1)*dpdr(i,741)   
      dpdr(i,1254) = dpdr(i,1)*p(742) + p(1)*dpdr(i,742)   
      dpdr(i,1255) = dpdr(i,1)*p(743) + p(1)*dpdr(i,743)   
      dpdr(i,1256) = dpdr(i,1)*p(744) + p(1)*dpdr(i,744)   
      dpdr(i,1257) = dpdr(i,1)*p(745) + p(1)*dpdr(i,745)   
      dpdr(i,1258) = dpdr(i,1)*p(746) + p(1)*dpdr(i,746)   
      dpdr(i,1259) = dpdr(i,1)*p(747) + p(1)*dpdr(i,747)   
      dpdr(i,1260) = dpdr(i,1)*p(748) + p(1)*dpdr(i,748)   
      dpdr(i,1261) = dpdr(i,1)*p(751) + p(1)*dpdr(i,751)   
      dpdr(i,1262) = dpdr(i,5)*p(426) + p(5)*dpdr(i,426) - dpdr(i,1066) 
      dpdr(i,1263) = dpdr(i,1)*p(753) + p(1)*dpdr(i,753)  
      dpdr(i,1264) = dpdr(i,1)*p(754) + p(1)*dpdr(i,754)  
      dpdr(i,1265) = dpdr(i,1)*p(755) + p(1)*dpdr(i,755)  
      dpdr(i,1266) = dpdr(i,1)*p(850) + p(1)*dpdr(i,850)  
      dpdr(i,1267) = dpdr(i,1)*p(756) + p(1)*dpdr(i,756)  
      dpdr(i,1268) = dpdr(i,1)*p(758) + p(1)*dpdr(i,758)  
      dpdr(i,1269) = dpdr(i,5)*p(419) + p(5)*dpdr(i,419) - dpdr(i,1042) 
      dpdr(i,1270) = dpdr(i,1)*p(760) + p(1)*dpdr(i,760)  
      dpdr(i,1271) = dpdr(i,5)*p(423) + p(5)*dpdr(i,423) - dpdr(i,1060) 
      dpdr(i,1272) = dpdr(i,1)*p(851) + p(1)*dpdr(i,851)  
      dpdr(i,1273) = dpdr(i,1)*p(762) + p(1)*dpdr(i,762)  
      dpdr(i,1274) = dpdr(i,5)*p(422) + p(5)*dpdr(i,422) - dpdr(i,1038) 
      dpdr(i,1275) = dpdr(i,6)*p(423) + p(6)*dpdr(i,423) - dpdr(i,1042) 
      dpdr(i,1276) = dpdr(i,1)*p(765) + p(1)*dpdr(i,765)  
      dpdr(i,1277) = dpdr(i,5)*p(425) + p(5)*dpdr(i,425) - dpdr(i,1045) 
      dpdr(i,1278) = dpdr(i,5)*p(437) + p(5)*dpdr(i,437) - dpdr(i,1078) 
      dpdr(i,1279) = dpdr(i,1)*p(768) + p(1)*dpdr(i,768)  
      dpdr(i,1280) = dpdr(i,5)*p(428) + p(5)*dpdr(i,428) - dpdr(i,1069) 
      dpdr(i,1281) = dpdr(i,5)*p(440) + p(5)*dpdr(i,440) - dpdr(i,1081) 
      dpdr(i,1282) = dpdr(i,1)*p(852) + p(1)*dpdr(i,852)  
      dpdr(i,1283) = dpdr(i,7)*p(428) + p(7)*dpdr(i,428) - dpdr(i,1045) 
      dpdr(i,1284) = dpdr(i,1)*p(771) + p(1)*dpdr(i,771)               
      dpdr(i,1285) = dpdr(i,1)*p(772) + p(1)*dpdr(i,772)               
      dpdr(i,1286) = dpdr(i,1)*p(854) + p(1)*dpdr(i,854)               
      dpdr(i,1287) = dpdr(i,1)*p(773) + p(1)*dpdr(i,773)               
      dpdr(i,1288) = dpdr(i,5)*p(432) + p(5)*dpdr(i,432) - dpdr(i,1058) 
      dpdr(i,1289) = dpdr(i,1)*p(775) + p(1)*dpdr(i,775)  
      dpdr(i,1290) = dpdr(i,5)*p(434) + p(5)*dpdr(i,434) - dpdr(i,1060) 
      dpdr(i,1291) = dpdr(i,1)*p(855) + p(1)*dpdr(i,855)  
      dpdr(i,1292) = dpdr(i,5)*p(494) + p(5)*dpdr(i,494) - dpdr(i,1275) 
      dpdr(i,1293) = dpdr(i,1)*p(777) + p(1)*dpdr(i,777)  
      dpdr(i,1294) = dpdr(i,5)*p(436) + p(5)*dpdr(i,436) - dpdr(i,1064) 
      dpdr(i,1295) = dpdr(i,6)*p(437) + p(6)*dpdr(i,437) - dpdr(i,1066) 
      dpdr(i,1296) = dpdr(i,1)*p(780) + p(1)*dpdr(i,780)  
      dpdr(i,1297) = dpdr(i,5)*p(439) + p(5)*dpdr(i,439) - dpdr(i,1069) 
      dpdr(i,1298) = dpdr(i,5)*p(446) + p(5)*dpdr(i,446) - dpdr(i,1085) 
      dpdr(i,1299) = dpdr(i,1)*p(856) + p(1)*dpdr(i,856)
      dpdr(i,1300) = dpdr(i,5)*p(495) + p(5)*dpdr(i,495)  - dpdr(i,1283)
      dpdr(i,1301) = dpdr(i,6)*p(495) + p(6)*dpdr(i,495)  - dpdr(i,1283)
      dpdr(i,1302) = dpdr(i,1)*p(783) + p(1)*dpdr(i,783)                
      dpdr(i,1303) = dpdr(i,1)*p(858) + p(1)*dpdr(i,858)                
      dpdr(i,1304) = dpdr(i,1)*p(784) + p(1)*dpdr(i,784)                
      dpdr(i,1305) = dpdr(i,5)*p(443) + p(5)*dpdr(i,443) - dpdr(i,1078) 
      dpdr(i,1306) = dpdr(i,1)*p(859) + p(1)*dpdr(i,859)  
      dpdr(i,1307) = dpdr(i,5)*p(498) + p(5)*dpdr(i,498) - dpdr(i,1295) 
      dpdr(i,1308) = dpdr(i,1)*p(786) + p(1)*dpdr(i,786)  
      dpdr(i,1309) = dpdr(i,5)*p(445) + p(5)*dpdr(i,445) - dpdr(i,1081) 
      dpdr(i,1310) = dpdr(i,6)*p(446) + p(6)*dpdr(i,446) - dpdr(i,1081) 
      dpdr(i,1311) = dpdr(i,1)*p(860) + p(1)*dpdr(i,860)  
      dpdr(i,1312) = dpdr(i,5)*p(499) + p(5)*dpdr(i,499) - dpdr(i,1301) 
      dpdr(i,1313) = dpdr(i,6)*p(499) + p(6)*dpdr(i,499) - dpdr(i,1300) 
      dpdr(i,1314) = dpdr(i,1)*p(862) + p(1)*dpdr(i,862)  
      dpdr(i,1315) = dpdr(i,1)*p(863) + p(1)*dpdr(i,863)  
      dpdr(i,1316) = dpdr(i,5)*p(502) + p(5)*dpdr(i,502) - dpdr(i,1310) 
      dpdr(i,1317) = dpdr(i,1)*p(864) + p(1)*dpdr(i,864)  
      dpdr(i,1318) = dpdr(i,5)*p(503) + p(5)*dpdr(i,503) - dpdr(i,1313) 
      dpdr(i,1319) = dpdr(i,6)*p(503) + p(6)*dpdr(i,503) - dpdr(i,1312) 
      dpdr(i,1320) = dpdr(i,1)*p(789) + p(1)*dpdr(i,789)   
      dpdr(i,1321) = dpdr(i,1)*p(790) + p(1)*dpdr(i,790)   
      dpdr(i,1322) = dpdr(i,1)*p(791) + p(1)*dpdr(i,791)   
      dpdr(i,1323) = dpdr(i,1)*p(792) + p(1)*dpdr(i,792)   
      dpdr(i,1324) = dpdr(i,3)*p(850) + p(3)*dpdr(i,850)   
      dpdr(i,1325) = dpdr(i,3)*p(851) + p(3)*dpdr(i,851)   
      dpdr(i,1326) = dpdr(i,3)*p(852) + p(3)*dpdr(i,852)   
      dpdr(i,1327) = dpdr(i,1)*p(796) + p(1)*dpdr(i,796)   
      dpdr(i,1328) = dpdr(i,3)*p(854) + p(3)*dpdr(i,854)   
      dpdr(i,1329) = dpdr(i,3)*p(855) + p(3)*dpdr(i,855)   
      dpdr(i,1330) = dpdr(i,3)*p(856) + p(3)*dpdr(i,856)   
      dpdr(i,1331) = dpdr(i,1)*p(800) + p(1)*dpdr(i,800)   
      dpdr(i,1332) = dpdr(i,3)*p(858) + p(3)*dpdr(i,858)   
      dpdr(i,1333) = dpdr(i,3)*p(859) + p(3)*dpdr(i,859)   
      dpdr(i,1334) = dpdr(i,3)*p(860) + p(3)*dpdr(i,860)   
      dpdr(i,1335) = dpdr(i,1)*p(866) + p(1)*dpdr(i,866)   
      dpdr(i,1336) = dpdr(i,3)*p(862) + p(3)*dpdr(i,862)   
      dpdr(i,1337) = dpdr(i,3)*p(863) + p(3)*dpdr(i,863)   
      dpdr(i,1338) = dpdr(i,3)*p(864) + p(3)*dpdr(i,864)   
      dpdr(i,1339) = dpdr(i,1)*p(804) + p(1)*dpdr(i,804)   
      dpdr(i,1340) = dpdr(i,1)*p(805) + p(1)*dpdr(i,805)   
      dpdr(i,1341) = dpdr(i,1)*p(806) + p(1)*dpdr(i,806)   
      dpdr(i,1342) = dpdr(i,1)*p(810) + p(1)*dpdr(i,810)   
      dpdr(i,1343) = dpdr(i,3)*p(793) + p(3)*dpdr(i,793)   
      dpdr(i,1344) = dpdr(i,3)*p(794) + p(3)*dpdr(i,794)   
      dpdr(i,1345) = dpdr(i,3)*p(795) + p(3)*dpdr(i,795)   
      dpdr(i,1346) = dpdr(i,1)*p(814) + p(1)*dpdr(i,814)   
      dpdr(i,1347) = dpdr(i,3)*p(797) + p(3)*dpdr(i,797)   
      dpdr(i,1348) = dpdr(i,3)*p(798) + p(3)*dpdr(i,798)   
      dpdr(i,1349) = dpdr(i,3)*p(799) + p(3)*dpdr(i,799)   
      dpdr(i,1350) = dpdr(i,1)*p(868) + p(1)*dpdr(i,868)   
      dpdr(i,1351) = dpdr(i,3)*p(801) + p(3)*dpdr(i,801)   
      dpdr(i,1352) = dpdr(i,3)*p(802) + p(3)*dpdr(i,802)   
      dpdr(i,1353) = dpdr(i,3)*p(803) + p(3)*dpdr(i,803)   
      dpdr(i,1354) = dpdr(i,1)*p(818) + p(1)*dpdr(i,818)   
      dpdr(i,1355) = dpdr(i,1)*p(819) + p(1)*dpdr(i,819)   
      dpdr(i,1356) = dpdr(i,1)*p(820) + p(1)*dpdr(i,820)   
      dpdr(i,1357) = dpdr(i,3)*p(807) + p(3)*dpdr(i,807)   
      dpdr(i,1358) = dpdr(i,3)*p(808) + p(3)*dpdr(i,808)   
      dpdr(i,1359) = dpdr(i,3)*p(809) + p(3)*dpdr(i,809)   
      dpdr(i,1360) = dpdr(i,1)*p(824) + p(1)*dpdr(i,824)   
      dpdr(i,1361) = dpdr(i,3)*p(811) + p(3)*dpdr(i,811)   
      dpdr(i,1362) = dpdr(i,3)*p(812) + p(3)*dpdr(i,812)   
      dpdr(i,1363) = dpdr(i,3)*p(813) + p(3)*dpdr(i,813)   
      dpdr(i,1364) = dpdr(i,1)*p(870) + p(1)*dpdr(i,870)   
      dpdr(i,1365) = dpdr(i,3)*p(815) + p(3)*dpdr(i,815)   
      dpdr(i,1366) = dpdr(i,3)*p(816) + p(3)*dpdr(i,816)   
      dpdr(i,1367) = dpdr(i,3)*p(817) + p(3)*dpdr(i,817)   
      dpdr(i,1368) = dpdr(i,1)*p(828) + p(1)*dpdr(i,828)   
      dpdr(i,1369) = dpdr(i,1)*p(829) + p(1)*dpdr(i,829)   
      dpdr(i,1370) = dpdr(i,1)*p(833) + p(1)*dpdr(i,833)   
      dpdr(i,1371) = dpdr(i,3)*p(821) + p(3)*dpdr(i,821)   
      dpdr(i,1372) = dpdr(i,3)*p(822) + p(3)*dpdr(i,822)   
      dpdr(i,1373) = dpdr(i,3)*p(823) + p(3)*dpdr(i,823)   
      dpdr(i,1374) = dpdr(i,1)*p(872) + p(1)*dpdr(i,872)   
      dpdr(i,1375) = dpdr(i,3)*p(825) + p(3)*dpdr(i,825)   
      dpdr(i,1376) = dpdr(i,3)*p(826) + p(3)*dpdr(i,826)   
      dpdr(i,1377) = dpdr(i,3)*p(827) + p(3)*dpdr(i,827)   
      dpdr(i,1378) = dpdr(i,1)*p(837) + p(1)*dpdr(i,837)   
      dpdr(i,1379) = dpdr(i,1)*p(838) + p(1)*dpdr(i,838)   
      dpdr(i,1380) = dpdr(i,3)*p(830) + p(3)*dpdr(i,830)   
      dpdr(i,1381) = dpdr(i,3)*p(831) + p(3)*dpdr(i,831)   
      dpdr(i,1382) = dpdr(i,3)*p(832) + p(3)*dpdr(i,832)   
      dpdr(i,1383) = dpdr(i,1)*p(874) + p(1)*dpdr(i,874)   
      dpdr(i,1384) = dpdr(i,3)*p(834) + p(3)*dpdr(i,834)   
      dpdr(i,1385) = dpdr(i,3)*p(835) + p(3)*dpdr(i,835)   
      dpdr(i,1386) = dpdr(i,3)*p(836) + p(3)*dpdr(i,836)   
      dpdr(i,1387) = dpdr(i,1)*p(842) + p(1)*dpdr(i,842)   
      dpdr(i,1388) = dpdr(i,1)*p(876) + p(1)*dpdr(i,876)   
      dpdr(i,1389) = dpdr(i,3)*p(839) + p(3)*dpdr(i,839)   
      dpdr(i,1390) = dpdr(i,3)*p(840) + p(3)*dpdr(i,840)   
      dpdr(i,1391) = dpdr(i,3)*p(841) + p(3)*dpdr(i,841)   
      dpdr(i,1392) = dpdr(i,1)*p(878) + p(1)*dpdr(i,878)   
      dpdr(i,1393) = dpdr(i,3)*p(843) + p(3)*dpdr(i,843)   
      dpdr(i,1394) = dpdr(i,3)*p(844) + p(3)*dpdr(i,844)   
      dpdr(i,1395) = dpdr(i,3)*p(845) + p(3)*dpdr(i,845)   
      dpdr(i,1396) = dpdr(i,1)*p(846) + p(1)*dpdr(i,846)   
      dpdr(i,1397) = dpdr(i,1)*p(847) + p(1)*dpdr(i,847)   
      dpdr(i,1398) = dpdr(i,1)*p(848) + p(1)*dpdr(i,848)   
      dpdr(i,1399) = dpdr(i,1)*p(849) + p(1)*dpdr(i,849)   
      dpdr(i,1400) = dpdr(i,1)*p(853) + p(1)*dpdr(i,853)   
      dpdr(i,1401) = dpdr(i,2)*p(850) + p(2)*dpdr(i,850) - dpdr(i,1274) 
      dpdr(i,1402) = dpdr(i,2)*p(851) + p(2)*dpdr(i,851) - dpdr(i,1275) 
      dpdr(i,1403) = dpdr(i,2)*p(852) + p(2)*dpdr(i,852) - dpdr(i,1283) 
      dpdr(i,1404) = dpdr(i,1)*p(857) + p(1)*dpdr(i,857)   
      dpdr(i,1405) = dpdr(i,5)*p(497) + p(5)*dpdr(i,497) - dpdr(i,1076) 
      dpdr(i,1406) = dpdr(i,6)*p(498) + p(6)*dpdr(i,498) - dpdr(i,1078) 
      dpdr(i,1407) = dpdr(i,7)*p(499) + p(7)*dpdr(i,499) - dpdr(i,1081) 
      dpdr(i,1408) = dpdr(i,1)*p(861) + p(1)*dpdr(i,861)   
      dpdr(i,1409) = dpdr(i,5)*p(501) + p(5)*dpdr(i,501) - dpdr(i,1085) 
      dpdr(i,1410) = dpdr(i,6)*p(502) + p(6)*dpdr(i,502) - dpdr(i,1085) 
      dpdr(i,1411) = dpdr(i,7)*p(503) + p(7)*dpdr(i,503) - dpdr(i,1085) 
      dpdr(i,1412) = dpdr(i,1)*p(880) + p(1)*dpdr(i,880)   
      dpdr(i,1413) = dpdr(i,5)*p(517) + p(5)*dpdr(i,517) - dpdr(i,1319) 
      dpdr(i,1414) = dpdr(i,6)*p(517) + p(6)*dpdr(i,517) - dpdr(i,1318) 
      dpdr(i,1415) = dpdr(i,7)*p(517) + p(7)*dpdr(i,517) - dpdr(i,1316) 
      dpdr(i,1416) = dpdr(i,1)*p(865) + p(1)*dpdr(i,865)                
      dpdr(i,1417) = dpdr(i,3)*p(880) + p(3)*dpdr(i,880)   
      dpdr(i,1418) = dpdr(i,1)*p(867) + p(1)*dpdr(i,867)   
      dpdr(i,1419) = dpdr(i,3)*p(866) + p(3)*dpdr(i,866)   
      dpdr(i,1420) = dpdr(i,1)*p(869) + p(1)*dpdr(i,869)   
      dpdr(i,1421) = dpdr(i,3)*p(868) + p(3)*dpdr(i,868)   
      dpdr(i,1422) = dpdr(i,1)*p(871) + p(1)*dpdr(i,871)   
      dpdr(i,1423) = dpdr(i,3)*p(870) + p(3)*dpdr(i,870)   
      dpdr(i,1424) = dpdr(i,1)*p(873) + p(1)*dpdr(i,873)   
      dpdr(i,1425) = dpdr(i,3)*p(872) + p(3)*dpdr(i,872)   
      dpdr(i,1426) = dpdr(i,1)*p(875) + p(1)*dpdr(i,875)   
      dpdr(i,1427) = dpdr(i,3)*p(874) + p(3)*dpdr(i,874)   
      dpdr(i,1428) = dpdr(i,1)*p(877) + p(1)*dpdr(i,877)   
      dpdr(i,1429) = dpdr(i,3)*p(876) + p(3)*dpdr(i,876)   
      dpdr(i,1430) = dpdr(i,1)*p(881) + p(1)*dpdr(i,881)   
      dpdr(i,1431) = dpdr(i,3)*p(878) + p(3)*dpdr(i,878)   
      dpdr(i,1432) = dpdr(i,1)*p(879) + p(1)*dpdr(i,879)   
      dpdr(i,1433) = dpdr(i,2)*p(880) + p(2)*dpdr(i,880) &
                   - dpdr(i,1415) - dpdr(i,1414) - dpdr(i,1413)  
      dpdr(i,1434) = dpdr(i,3)*p(881) + p(3)*dpdr(i,881)   
      dpdr(i,1435) = dpdr(i,1)*p(882) + p(1)*dpdr(i,882)   
      dpdr(i,1436) = dpdr(i,1)*p(883) + p(1)*dpdr(i,883)   
      dpdr(i,1437) = dpdr(i,1)*p(884) + p(1)*dpdr(i,884)   
      dpdr(i,1438) = dpdr(i,1)*p(885) + p(1)*dpdr(i,885)   
      dpdr(i,1439) = dpdr(i,1)*p(886) + p(1)*dpdr(i,886)   
      dpdr(i,1440) = dpdr(i,1)*p(887) + p(1)*dpdr(i,887)   
      dpdr(i,1441) = dpdr(i,1)*p(932) + p(1)*dpdr(i,932)   
      dpdr(i,1442) = dpdr(i,1)*p(888) + p(1)*dpdr(i,888)   
      dpdr(i,1443) = dpdr(i,1)*p(889) + p(1)*dpdr(i,889)   
      dpdr(i,1444) = dpdr(i,1)*p(941) + p(1)*dpdr(i,941)   
      dpdr(i,1445) = dpdr(i,1)*p(890) + p(1)*dpdr(i,890)   
      dpdr(i,1446) = dpdr(i,1)*p(943) + p(1)*dpdr(i,943)   
      dpdr(i,1447) = dpdr(i,1)*p(945) + p(1)*dpdr(i,945)   
      dpdr(i,1448) = dpdr(i,1)*p(891) + p(1)*dpdr(i,891)   
      dpdr(i,1449) = dpdr(i,1)*p(949) + p(1)*dpdr(i,949)   
      dpdr(i,1450) = dpdr(i,1)*p(951) + p(1)*dpdr(i,951)   
      dpdr(i,1451) = dpdr(i,1)*p(953) + p(1)*dpdr(i,953)   
      dpdr(i,1452) = dpdr(i,1)*p(892) + p(1)*dpdr(i,892)   
      dpdr(i,1453) = dpdr(i,1)*p(960) + p(1)*dpdr(i,960)   
      dpdr(i,1454) = dpdr(i,1)*p(962) + p(1)*dpdr(i,962)   
      dpdr(i,1455) = dpdr(i,1)*p(965) + p(1)*dpdr(i,965)   
      dpdr(i,1456) = dpdr(i,1)*p(969) + p(1)*dpdr(i,969)   
      dpdr(i,1457) = dpdr(i,1)*p(893) + p(1)*dpdr(i,893)   
      dpdr(i,1458) = dpdr(i,1)*p(894) + p(1)*dpdr(i,894)   
      dpdr(i,1459) = dpdr(i,1)*p(895) + p(1)*dpdr(i,895)   
      dpdr(i,1460) = dpdr(i,1)*p(896) + p(1)*dpdr(i,896)   
      dpdr(i,1461) = dpdr(i,1)*p(897) + p(1)*dpdr(i,897)   
      dpdr(i,1462) = dpdr(i,1)*p(978) + p(1)*dpdr(i,978)   
      dpdr(i,1463) = dpdr(i,1)*p(898) + p(1)*dpdr(i,898)   
      dpdr(i,1464) = dpdr(i,1)*p(985) + p(1)*dpdr(i,985)   
      dpdr(i,1465) = dpdr(i,1)*p(987) + p(1)*dpdr(i,987)   
      dpdr(i,1466) = dpdr(i,1)*p(990) + p(1)*dpdr(i,990)   
      dpdr(i,1467) = dpdr(i,1)*p(994) + p(1)*dpdr(i,994)   
      dpdr(i,1468) = dpdr(i,1)*p(899) + p(1)*dpdr(i,899)   
      dpdr(i,1469) = dpdr(i,1)*p(900) + p(1)*dpdr(i,900)   
      dpdr(i,1470) = dpdr(i,1)*p(1000) + p(1)*dpdr(i,1000)  
      dpdr(i,1471) = dpdr(i,1)*p(1002) + p(1)*dpdr(i,1002)  
      dpdr(i,1472) = dpdr(i,1)*p(1004) + p(1)*dpdr(i,1004)  
      dpdr(i,1473) = dpdr(i,1)*p(1008) + p(1)*dpdr(i,1008)  
      dpdr(i,1474) = dpdr(i,1)*p(901) + p(1)*dpdr(i,901)   
      dpdr(i,1475) = dpdr(i,1)*p(1013) + p(1)*dpdr(i,1013)  
      dpdr(i,1476) = dpdr(i,1)*p(1015) + p(1)*dpdr(i,1015)  
      dpdr(i,1477) = dpdr(i,1)*p(902) + p(1)*dpdr(i,902)   
      dpdr(i,1478) = dpdr(i,1)*p(903) + p(1)*dpdr(i,903)   
      dpdr(i,1479) = dpdr(i,1)*p(904) + p(1)*dpdr(i,904)   
      dpdr(i,1480) = dpdr(i,1)*p(905) + p(1)*dpdr(i,905)   
      dpdr(i,1481) = dpdr(i,1)*p(906) + p(1)*dpdr(i,906)   
      dpdr(i,1482) = dpdr(i,1)*p(907) + p(1)*dpdr(i,907)   
      dpdr(i,1483) = dpdr(i,1)*p(908) + p(1)*dpdr(i,908)   
      dpdr(i,1484) = dpdr(i,1)*p(909) + p(1)*dpdr(i,909)   
      dpdr(i,1485) = dpdr(i,1)*p(910) + p(1)*dpdr(i,910)   
      dpdr(i,1486) = dpdr(i,1)*p(911) + p(1)*dpdr(i,911)   
      dpdr(i,1487) = dpdr(i,1)*p(912) + p(1)*dpdr(i,912)   
      dpdr(i,1488) = dpdr(i,1)*p(913) + p(1)*dpdr(i,913)   
      dpdr(i,1489) = dpdr(i,1)*p(914) + p(1)*dpdr(i,914)   
      dpdr(i,1490) = dpdr(i,1)*p(915) + p(1)*dpdr(i,915)   
      dpdr(i,1491) = dpdr(i,1)*p(916) + p(1)*dpdr(i,916)   
      dpdr(i,1492) = dpdr(i,1)*p(1038) + p(1)*dpdr(i,1038)  
      dpdr(i,1493) = dpdr(i,1)*p(917) + p(1)*dpdr(i,917)   
      dpdr(i,1494) = dpdr(i,1)*p(1040) + p(1)*dpdr(i,1040)  
      dpdr(i,1495) = dpdr(i,1)*p(1042) + p(1)*dpdr(i,1042)  
      dpdr(i,1496) = dpdr(i,1)*p(918) + p(1)*dpdr(i,918)   
      dpdr(i,1497) = dpdr(i,1)*p(1045) + p(1)*dpdr(i,1045)  
      dpdr(i,1498) = dpdr(i,1)*p(919) + p(1)*dpdr(i,919)   
      dpdr(i,1499) = dpdr(i,1)*p(920) + p(1)*dpdr(i,920)   
      dpdr(i,1500) = dpdr(i,1)*p(1056) + p(1)*dpdr(i,1056)  
      dpdr(i,1501) = dpdr(i,1)*p(921) + p(1)*dpdr(i,921)   
      dpdr(i,1502) = dpdr(i,1)*p(1058) + p(1)*dpdr(i,1058)  
      dpdr(i,1503) = dpdr(i,1)*p(1060) + p(1)*dpdr(i,1060)  
      dpdr(i,1504) = dpdr(i,1)*p(922) + p(1)*dpdr(i,922)   
      dpdr(i,1505) = dpdr(i,1)*p(1064) + p(1)*dpdr(i,1064)  
      dpdr(i,1506) = dpdr(i,1)*p(1066) + p(1)*dpdr(i,1066)  
      dpdr(i,1507) = dpdr(i,1)*p(1069) + p(1)*dpdr(i,1069)  
      dpdr(i,1508) = dpdr(i,1)*p(923) + p(1)*dpdr(i,923)   
      dpdr(i,1509) = dpdr(i,1)*p(1076) + p(1)*dpdr(i,1076)  
      dpdr(i,1510) = dpdr(i,1)*p(1078) + p(1)*dpdr(i,1078)  
      dpdr(i,1511) = dpdr(i,1)*p(1081) + p(1)*dpdr(i,1081)  
      dpdr(i,1512) = dpdr(i,1)*p(1085) + p(1)*dpdr(i,1085)  
      dpdr(i,1513) = dpdr(i,1)*p(924) + p(1)*dpdr(i,924)   
      dpdr(i,1514) = dpdr(i,1)*p(925) + p(1)*dpdr(i,925)   
      dpdr(i,1515) = dpdr(i,1)*p(926) + p(1)*dpdr(i,926)   
      dpdr(i,1516) = dpdr(i,1)*p(927) + p(1)*dpdr(i,927)   
      dpdr(i,1517) = dpdr(i,1)*p(928) + p(1)*dpdr(i,928)   
      dpdr(i,1518) = dpdr(i,1)*p(929) + p(1)*dpdr(i,929)   
      dpdr(i,1519) = dpdr(i,1)*p(930) + p(1)*dpdr(i,930)   
      dpdr(i,1520) = dpdr(i,1)*p(931) + p(1)*dpdr(i,931)   
      dpdr(i,1521) = dpdr(i,1)*p(933) + p(1)*dpdr(i,933)   
      dpdr(i,1522) = dpdr(i,1)*p(934) + p(1)*dpdr(i,934)   
      dpdr(i,1523) = dpdr(i,1)*p(935) + p(1)*dpdr(i,935)   
      dpdr(i,1524) = dpdr(i,1)*p(936) + p(1)*dpdr(i,936)   
      dpdr(i,1525) = dpdr(i,1)*p(937) + p(1)*dpdr(i,937)   
      dpdr(i,1526) = dpdr(i,1)*p(938) + p(1)*dpdr(i,938)   
      dpdr(i,1527) = dpdr(i,1)*p(939) + p(1)*dpdr(i,939)   
      dpdr(i,1528) = dpdr(i,1)*p(940) + p(1)*dpdr(i,940)   
      dpdr(i,1529) = dpdr(i,1)*p(942) + p(1)*dpdr(i,942)   
      dpdr(i,1530) = dpdr(i,1)*p(944) + p(1)*dpdr(i,944)   
      dpdr(i,1531) = dpdr(i,1)*p(946) + p(1)*dpdr(i,946)   
      dpdr(i,1532) = dpdr(i,1)*p(947) + p(1)*dpdr(i,947)   
      dpdr(i,1533) = dpdr(i,1)*p(1103) + p(1)*dpdr(i,1103)  
      dpdr(i,1534) = dpdr(i,1)*p(948) + p(1)*dpdr(i,948)   
      dpdr(i,1535) = dpdr(i,3)*p(1038) + p(3)*dpdr(i,1038)  
      dpdr(i,1536) = dpdr(i,1)*p(950) + p(1)*dpdr(i,950)   
      dpdr(i,1537) = dpdr(i,3)*p(1040) + p(3)*dpdr(i,1040)  
      dpdr(i,1538) = dpdr(i,1)*p(1104) + p(1)*dpdr(i,1104)  
      dpdr(i,1539) = dpdr(i,3)*p(1042) + p(3)*dpdr(i,1042)  
      dpdr(i,1540) = dpdr(i,1)*p(952) + p(1)*dpdr(i,952)   
      dpdr(i,1541) = dpdr(i,1)*p(1106) + p(1)*dpdr(i,1106)  
      dpdr(i,1542) = dpdr(i,3)*p(1045) + p(3)*dpdr(i,1045)  
      dpdr(i,1543) = dpdr(i,1)*p(954) + p(1)*dpdr(i,954)   
      dpdr(i,1544) = dpdr(i,1)*p(955) + p(1)*dpdr(i,955)   
      dpdr(i,1545) = dpdr(i,1)*p(1111) + p(1)*dpdr(i,1111)  
      dpdr(i,1546) = dpdr(i,1)*p(956) + p(1)*dpdr(i,956)   
      dpdr(i,1547) = dpdr(i,1)*p(1113) + p(1)*dpdr(i,1113)  
      dpdr(i,1548) = dpdr(i,1)*p(1115) + p(1)*dpdr(i,1115)  
      dpdr(i,1549) = dpdr(i,1)*p(957) + p(1)*dpdr(i,957)   
      dpdr(i,1550) = dpdr(i,1)*p(958) + p(1)*dpdr(i,958)   
      dpdr(i,1551) = dpdr(i,1)*p(1117) + p(1)*dpdr(i,1117)  
      dpdr(i,1552) = dpdr(i,1)*p(959) + p(1)*dpdr(i,959)   
      dpdr(i,1553) = dpdr(i,3)*p(1056) + p(3)*dpdr(i,1056)  
      dpdr(i,1554) = dpdr(i,1)*p(961) + p(1)*dpdr(i,961)   
      dpdr(i,1555) = dpdr(i,3)*p(1058) + p(3)*dpdr(i,1058)  
      dpdr(i,1556) = dpdr(i,1)*p(1118) + p(1)*dpdr(i,1118)  
      dpdr(i,1557) = dpdr(i,3)*p(1060) + p(3)*dpdr(i,1060)  
      dpdr(i,1558) = dpdr(i,1)*p(963) + p(1)*dpdr(i,963)   
      dpdr(i,1559) = dpdr(i,1)*p(1120) + p(1)*dpdr(i,1120)  
      dpdr(i,1560) = dpdr(i,1)*p(964) + p(1)*dpdr(i,964)   
      dpdr(i,1561) = dpdr(i,3)*p(1064) + p(3)*dpdr(i,1064)  
      dpdr(i,1562) = dpdr(i,1)*p(1121) + p(1)*dpdr(i,1121)  
      dpdr(i,1563) = dpdr(i,3)*p(1066) + p(3)*dpdr(i,1066)  
      dpdr(i,1564) = dpdr(i,1)*p(1123) + p(1)*dpdr(i,1123)  
      dpdr(i,1565) = dpdr(i,1)*p(1124) + p(1)*dpdr(i,1124)  
      dpdr(i,1566) = dpdr(i,3)*p(1069) + p(3)*dpdr(i,1069)  
      dpdr(i,1567) = dpdr(i,1)*p(966) + p(1)*dpdr(i,966)   
      dpdr(i,1568) = dpdr(i,1)*p(1128) + p(1)*dpdr(i,1128)  
      dpdr(i,1569) = dpdr(i,1)*p(1130) + p(1)*dpdr(i,1130)  
      dpdr(i,1570) = dpdr(i,1)*p(967) + p(1)*dpdr(i,967)   
      dpdr(i,1571) = dpdr(i,1)*p(1132) + p(1)*dpdr(i,1132)  
      dpdr(i,1572) = dpdr(i,1)*p(968) + p(1)*dpdr(i,968)   
      dpdr(i,1573) = dpdr(i,3)*p(1076) + p(3)*dpdr(i,1076)  
      dpdr(i,1574) = dpdr(i,1)*p(1133) + p(1)*dpdr(i,1133)  
      dpdr(i,1575) = dpdr(i,3)*p(1078) + p(3)*dpdr(i,1078)  
      dpdr(i,1576) = dpdr(i,1)*p(1135) + p(1)*dpdr(i,1135)  
      dpdr(i,1577) = dpdr(i,1)*p(1136) + p(1)*dpdr(i,1136)  
      dpdr(i,1578) = dpdr(i,3)*p(1081) + p(3)*dpdr(i,1081)  
      dpdr(i,1579) = dpdr(i,1)*p(1139) + p(1)*dpdr(i,1139)  
      dpdr(i,1580) = dpdr(i,1)*p(1141) + p(1)*dpdr(i,1141)  
      dpdr(i,1581) = dpdr(i,1)*p(1142) + p(1)*dpdr(i,1142)  
      dpdr(i,1582) = dpdr(i,3)*p(1085) + p(3)*dpdr(i,1085)  
      dpdr(i,1583) = dpdr(i,1)*p(970) + p(1)*dpdr(i,970)   
      dpdr(i,1584) = dpdr(i,1)*p(971) + p(1)*dpdr(i,971)   
      dpdr(i,1585) = dpdr(i,1)*p(972) + p(1)*dpdr(i,972)   
      dpdr(i,1586) = dpdr(i,1)*p(973) + p(1)*dpdr(i,973)   
      dpdr(i,1587) = dpdr(i,1)*p(974) + p(1)*dpdr(i,974)   
      dpdr(i,1588) = dpdr(i,1)*p(975) + p(1)*dpdr(i,975)   
      dpdr(i,1589) = dpdr(i,1)*p(976) + p(1)*dpdr(i,976)   
      dpdr(i,1590) = dpdr(i,1)*p(977) + p(1)*dpdr(i,977)   
      dpdr(i,1591) = dpdr(i,3)*p(932) + p(3)*dpdr(i,932)   
      dpdr(i,1592) = dpdr(i,1)*p(979) + p(1)*dpdr(i,979)   
      dpdr(i,1593) = dpdr(i,1)*p(980) + p(1)*dpdr(i,980)   
      dpdr(i,1594) = dpdr(i,1)*p(981) + p(1)*dpdr(i,981)   
      dpdr(i,1595) = dpdr(i,1)*p(1157) + p(1)*dpdr(i,1157)  
      dpdr(i,1596) = dpdr(i,1)*p(982) + p(1)*dpdr(i,982)   
      dpdr(i,1597) = dpdr(i,1)*p(983) + p(1)*dpdr(i,983)   
      dpdr(i,1598) = dpdr(i,1)*p(1160) + p(1)*dpdr(i,1160)  
      dpdr(i,1599) = dpdr(i,1)*p(984) + p(1)*dpdr(i,984)   
      dpdr(i,1600) = dpdr(i,3)*p(941) + p(3)*dpdr(i,941)   
      dpdr(i,1601) = dpdr(i,1)*p(986) + p(1)*dpdr(i,986)   
      dpdr(i,1602) = dpdr(i,3)*p(943) + p(3)*dpdr(i,943)   
      dpdr(i,1603) = dpdr(i,1)*p(1161) + p(1)*dpdr(i,1161)  
      dpdr(i,1604) = dpdr(i,3)*p(945) + p(3)*dpdr(i,945)   
      dpdr(i,1605) = dpdr(i,1)*p(988) + p(1)*dpdr(i,988)   
      dpdr(i,1606) = dpdr(i,1)*p(1163) + p(1)*dpdr(i,1163)  
      dpdr(i,1607) = dpdr(i,1)*p(989) + p(1)*dpdr(i,989)   
      dpdr(i,1608) = dpdr(i,3)*p(949) + p(3)*dpdr(i,949)   
      dpdr(i,1609) = dpdr(i,1)*p(1164) + p(1)*dpdr(i,1164)  
      dpdr(i,1610) = dpdr(i,3)*p(951) + p(3)*dpdr(i,951)   
      dpdr(i,1611) = dpdr(i,1)*p(1166) + p(1)*dpdr(i,1166)  
      dpdr(i,1612) = dpdr(i,3)*p(953) + p(3)*dpdr(i,953)   
      dpdr(i,1613) = dpdr(i,1)*p(991) + p(1)*dpdr(i,991)   
      dpdr(i,1614) = dpdr(i,1)*p(1170) + p(1)*dpdr(i,1170)  
      dpdr(i,1615) = dpdr(i,1)*p(1172) + p(1)*dpdr(i,1172)  
      dpdr(i,1616) = dpdr(i,1)*p(992) + p(1)*dpdr(i,992)   
      dpdr(i,1617) = dpdr(i,1)*p(1174) + p(1)*dpdr(i,1174)  
      dpdr(i,1618) = dpdr(i,1)*p(993) + p(1)*dpdr(i,993)   
      dpdr(i,1619) = dpdr(i,3)*p(960) + p(3)*dpdr(i,960)   
      dpdr(i,1620) = dpdr(i,1)*p(1175) + p(1)*dpdr(i,1175)  
      dpdr(i,1621) = dpdr(i,3)*p(962) + p(3)*dpdr(i,962)   
      dpdr(i,1622) = dpdr(i,1)*p(1177) + p(1)*dpdr(i,1177)  
      dpdr(i,1623) = dpdr(i,1)*p(1178) + p(1)*dpdr(i,1178)  
      dpdr(i,1624) = dpdr(i,3)*p(965) + p(3)*dpdr(i,965)   
      dpdr(i,1625) = dpdr(i,1)*p(1181) + p(1)*dpdr(i,1181)  
      dpdr(i,1626) = dpdr(i,1)*p(1183) + p(1)*dpdr(i,1183)  
      dpdr(i,1627) = dpdr(i,1)*p(1184) + p(1)*dpdr(i,1184)  
      dpdr(i,1628) = dpdr(i,3)*p(969) + p(3)*dpdr(i,969)   
      dpdr(i,1629) = dpdr(i,1)*p(995) + p(1)*dpdr(i,995)   
      dpdr(i,1630) = dpdr(i,1)*p(996) + p(1)*dpdr(i,996)   
      dpdr(i,1631) = dpdr(i,1)*p(997) + p(1)*dpdr(i,997)   
      dpdr(i,1632) = dpdr(i,1)*p(998) + p(1)*dpdr(i,998)   
      dpdr(i,1633) = dpdr(i,1)*p(999) + p(1)*dpdr(i,999)   
      dpdr(i,1634) = dpdr(i,1)*p(1001) + p(1)*dpdr(i,1001)  
      dpdr(i,1635) = dpdr(i,1)*p(1003) + p(1)*dpdr(i,1003)  
      dpdr(i,1636) = dpdr(i,1)*p(1194) + p(1)*dpdr(i,1194)  
      dpdr(i,1637) = dpdr(i,3)*p(978) + p(3)*dpdr(i,978)   
      dpdr(i,1638) = dpdr(i,1)*p(1005) + p(1)*dpdr(i,1005)  
      dpdr(i,1639) = dpdr(i,1)*p(1198) + p(1)*dpdr(i,1198)  
      dpdr(i,1640) = dpdr(i,1)*p(1200) + p(1)*dpdr(i,1200)  
      dpdr(i,1641) = dpdr(i,1)*p(1006) + p(1)*dpdr(i,1006)  
      dpdr(i,1642) = dpdr(i,1)*p(1202) + p(1)*dpdr(i,1202)  
      dpdr(i,1643) = dpdr(i,1)*p(1007) + p(1)*dpdr(i,1007)  
      dpdr(i,1644) = dpdr(i,3)*p(985) + p(3)*dpdr(i,985)   
      dpdr(i,1645) = dpdr(i,1)*p(1203) + p(1)*dpdr(i,1203)  
      dpdr(i,1646) = dpdr(i,3)*p(987) + p(3)*dpdr(i,987)   
      dpdr(i,1647) = dpdr(i,1)*p(1205) + p(1)*dpdr(i,1205)  
      dpdr(i,1648) = dpdr(i,1)*p(1206) + p(1)*dpdr(i,1206)  
      dpdr(i,1649) = dpdr(i,3)*p(990) + p(3)*dpdr(i,990)   
      dpdr(i,1650) = dpdr(i,1)*p(1209) + p(1)*dpdr(i,1209)  
      dpdr(i,1651) = dpdr(i,1)*p(1211) + p(1)*dpdr(i,1211)  
      dpdr(i,1652) = dpdr(i,1)*p(1212) + p(1)*dpdr(i,1212)  
      dpdr(i,1653) = dpdr(i,3)*p(994) + p(3)*dpdr(i,994)   
      dpdr(i,1654) = dpdr(i,1)*p(1009) + p(1)*dpdr(i,1009)  
      dpdr(i,1655) = dpdr(i,1)*p(1010) + p(1)*dpdr(i,1010)  
      dpdr(i,1656) = dpdr(i,1)*p(1011) + p(1)*dpdr(i,1011)  
      dpdr(i,1657) = dpdr(i,1)*p(1221) + p(1)*dpdr(i,1221)  
      dpdr(i,1658) = dpdr(i,1)*p(1012) + p(1)*dpdr(i,1012)  
      dpdr(i,1659) = dpdr(i,3)*p(1000) + p(3)*dpdr(i,1000)  
      dpdr(i,1660) = dpdr(i,1)*p(1222) + p(1)*dpdr(i,1222)  
      dpdr(i,1661) = dpdr(i,3)*p(1002) + p(3)*dpdr(i,1002)  
      dpdr(i,1662) = dpdr(i,1)*p(1224) + p(1)*dpdr(i,1224)  
      dpdr(i,1663) = dpdr(i,3)*p(1004) + p(3)*dpdr(i,1004)  
      dpdr(i,1664) = dpdr(i,1)*p(1227) + p(1)*dpdr(i,1227)  
      dpdr(i,1665) = dpdr(i,1)*p(1229) + p(1)*dpdr(i,1229)  
      dpdr(i,1666) = dpdr(i,1)*p(1230) + p(1)*dpdr(i,1230)  
      dpdr(i,1667) = dpdr(i,3)*p(1008) + p(3)*dpdr(i,1008)  
      dpdr(i,1668) = dpdr(i,1)*p(1014) + p(1)*dpdr(i,1014)  
      dpdr(i,1669) = dpdr(i,1)*p(1236) + p(1)*dpdr(i,1236)  
      dpdr(i,1670) = dpdr(i,1)*p(1238) + p(1)*dpdr(i,1238)  
      dpdr(i,1671) = dpdr(i,1)*p(1239) + p(1)*dpdr(i,1239)  
      dpdr(i,1672) = dpdr(i,3)*p(1013) + p(3)*dpdr(i,1013)  
      dpdr(i,1673) = dpdr(i,1)*p(1243) + p(1)*dpdr(i,1243)  
      dpdr(i,1674) = dpdr(i,3)*p(1015) + p(3)*dpdr(i,1015)  
      dpdr(i,1675) = dpdr(i,1)*p(1016) + p(1)*dpdr(i,1016)  
      dpdr(i,1676) = dpdr(i,1)*p(1017) + p(1)*dpdr(i,1017)  
      dpdr(i,1677) = dpdr(i,1)*p(1018) + p(1)*dpdr(i,1018)  
      dpdr(i,1678) = dpdr(i,1)*p(1019) + p(1)*dpdr(i,1019)  
      dpdr(i,1679) = dpdr(i,1)*p(1020) + p(1)*dpdr(i,1020)  
      dpdr(i,1680) = dpdr(i,1)*p(1021) + p(1)*dpdr(i,1021)  
      dpdr(i,1681) = dpdr(i,1)*p(1022) + p(1)*dpdr(i,1022)  
      dpdr(i,1682) = dpdr(i,1)*p(1023) + p(1)*dpdr(i,1023)  
      dpdr(i,1683) = dpdr(i,1)*p(1024) + p(1)*dpdr(i,1024)  
      dpdr(i,1684) = dpdr(i,1)*p(1025) + p(1)*dpdr(i,1025)  
      dpdr(i,1685) = dpdr(i,1)*p(1026) + p(1)*dpdr(i,1026)  
      dpdr(i,1686) = dpdr(i,1)*p(1027) + p(1)*dpdr(i,1027)  
      dpdr(i,1687) = dpdr(i,1)*p(1028) + p(1)*dpdr(i,1028)  
      dpdr(i,1688) = dpdr(i,1)*p(1029) + p(1)*dpdr(i,1029)  
      dpdr(i,1689) = dpdr(i,1)*p(1030) + p(1)*dpdr(i,1030)  
      dpdr(i,1690) = dpdr(i,1)*p(1031) + p(1)*dpdr(i,1031)  
      dpdr(i,1691) = dpdr(i,1)*p(1032) + p(1)*dpdr(i,1032)  
      dpdr(i,1692) = dpdr(i,1)*p(1033) + p(1)*dpdr(i,1033)  
      dpdr(i,1693) = dpdr(i,1)*p(1034) + p(1)*dpdr(i,1034)  
      dpdr(i,1694) = dpdr(i,1)*p(1035) + p(1)*dpdr(i,1035)  
      dpdr(i,1695) = dpdr(i,1)*p(1036) + p(1)*dpdr(i,1036)  
      dpdr(i,1696) = dpdr(i,1)*p(1037) + p(1)*dpdr(i,1037)  
      dpdr(i,1697) = dpdr(i,1)*p(1039) + p(1)*dpdr(i,1039)  
      dpdr(i,1698) = dpdr(i,34)*p(168) + p(34)*dpdr(i,168)  
      dpdr(i,1699) = dpdr(i,1)*p(1041) + p(1)*dpdr(i,1041)  
      dpdr(i,1700) = dpdr(i,34)*p(170) + p(34)*dpdr(i,170)  
      dpdr(i,1701) = dpdr(i,1)*p(1043) + p(1)*dpdr(i,1043)  
      dpdr(i,1702) = dpdr(i,1)*p(1044) + p(1)*dpdr(i,1044)  
      dpdr(i,1703) = dpdr(i,34)*p(172) + p(34)*dpdr(i,172)  
      dpdr(i,1704) = dpdr(i,1)*p(1262) + p(1)*dpdr(i,1262)  
      dpdr(i,1705) = dpdr(i,34)*p(217) + p(34)*dpdr(i,217)  
      dpdr(i,1706) = dpdr(i,1)*p(1046) + p(1)*dpdr(i,1046)  
      dpdr(i,1707) = dpdr(i,1)*p(1047) + p(1)*dpdr(i,1047)  
      dpdr(i,1708) = dpdr(i,1)*p(1048) + p(1)*dpdr(i,1048)  
      dpdr(i,1709) = dpdr(i,1)*p(1049) + p(1)*dpdr(i,1049)  
      dpdr(i,1710) = dpdr(i,1)*p(1050) + p(1)*dpdr(i,1050)  
      dpdr(i,1711) = dpdr(i,1)*p(1269) + p(1)*dpdr(i,1269)  
      dpdr(i,1712) = dpdr(i,1)*p(1051) + p(1)*dpdr(i,1051)  
      dpdr(i,1713) = dpdr(i,1)*p(1271) + p(1)*dpdr(i,1271)  
      dpdr(i,1714) = dpdr(i,1)*p(1052) + p(1)*dpdr(i,1052)  
      dpdr(i,1715) = dpdr(i,1)*p(1053) + p(1)*dpdr(i,1053)  
      dpdr(i,1716) = dpdr(i,1)*p(1054) + p(1)*dpdr(i,1054)  
      dpdr(i,1717) = dpdr(i,1)*p(1274) + p(1)*dpdr(i,1274)  
      dpdr(i,1718) = dpdr(i,1)*p(1055) + p(1)*dpdr(i,1055)  
      dpdr(i,1719) = dpdr(i,34)*p(267) + p(34)*dpdr(i,267)  
      dpdr(i,1720) = dpdr(i,1)*p(1057) + p(1)*dpdr(i,1057)  
      dpdr(i,1721) = dpdr(i,34)*p(221) + p(34)*dpdr(i,221)  
      dpdr(i,1722) = dpdr(i,1)*p(1059) + p(1)*dpdr(i,1059)  
      dpdr(i,1723) = dpdr(i,34)*p(223) + p(34)*dpdr(i,223)  
      dpdr(i,1724) = dpdr(i,1)*p(1275) + p(1)*dpdr(i,1275)  
      dpdr(i,1725) = dpdr(i,34)*p(268) + p(34)*dpdr(i,268)  
      dpdr(i,1726) = dpdr(i,1)*p(1061) + p(1)*dpdr(i,1061)  
      dpdr(i,1727) = dpdr(i,1)*p(1062) + p(1)*dpdr(i,1062)  
      dpdr(i,1728) = dpdr(i,1)*p(1277) + p(1)*dpdr(i,1277)  
      dpdr(i,1729) = dpdr(i,1)*p(1063) + p(1)*dpdr(i,1063)  
      dpdr(i,1730) = dpdr(i,34)*p(225) + p(34)*dpdr(i,225)  
      dpdr(i,1731) = dpdr(i,1)*p(1065) + p(1)*dpdr(i,1065)  
      dpdr(i,1732) = dpdr(i,34)*p(176) + p(34)*dpdr(i,176)  
      dpdr(i,1733) = dpdr(i,1)*p(1278) + p(1)*dpdr(i,1278)  
      dpdr(i,1734) = dpdr(i,34)*p(226) + p(34)*dpdr(i,226)  
      dpdr(i,1735) = dpdr(i,1)*p(1067) + p(1)*dpdr(i,1067)  
      dpdr(i,1736) = dpdr(i,1)*p(1280) + p(1)*dpdr(i,1280)  
      dpdr(i,1737) = dpdr(i,1)*p(1068) + p(1)*dpdr(i,1068)  
      dpdr(i,1738) = dpdr(i,34)*p(228) + p(34)*dpdr(i,228)  
      dpdr(i,1739) = dpdr(i,1)*p(1281) + p(1)*dpdr(i,1281)  
      dpdr(i,1740) = dpdr(i,34)*p(229) + p(34)*dpdr(i,229)  
      dpdr(i,1741) = dpdr(i,1)*p(1283) + p(1)*dpdr(i,1283)  
      dpdr(i,1742) = dpdr(i,34)*p(269) + p(34)*dpdr(i,269)  
      dpdr(i,1743) = dpdr(i,1)*p(1070) + p(1)*dpdr(i,1070)  
      dpdr(i,1744) = dpdr(i,1)*p(1071) + p(1)*dpdr(i,1071)  
      dpdr(i,1745) = dpdr(i,1)*p(1288) + p(1)*dpdr(i,1288)  
      dpdr(i,1746) = dpdr(i,1)*p(1072) + p(1)*dpdr(i,1072)  
      dpdr(i,1747) = dpdr(i,1)*p(1290) + p(1)*dpdr(i,1290)  
      dpdr(i,1748) = dpdr(i,1)*p(1292) + p(1)*dpdr(i,1292)  
      dpdr(i,1749) = dpdr(i,1)*p(1073) + p(1)*dpdr(i,1073)  
      dpdr(i,1750) = dpdr(i,1)*p(1074) + p(1)*dpdr(i,1074)  
      dpdr(i,1751) = dpdr(i,1)*p(1294) + p(1)*dpdr(i,1294)  
      dpdr(i,1752) = dpdr(i,1)*p(1075) + p(1)*dpdr(i,1075)  
      dpdr(i,1753) = dpdr(i,34)*p(271) + p(34)*dpdr(i,271)  
      dpdr(i,1754) = dpdr(i,1)*p(1077) + p(1)*dpdr(i,1077)  
      dpdr(i,1755) = dpdr(i,34)*p(232) + p(34)*dpdr(i,232)  
      dpdr(i,1756) = dpdr(i,1)*p(1295) + p(1)*dpdr(i,1295)  
      dpdr(i,1757) = dpdr(i,34)*p(272) + p(34)*dpdr(i,272)  
      dpdr(i,1758) = dpdr(i,1)*p(1079) + p(1)*dpdr(i,1079)  
      dpdr(i,1759) = dpdr(i,1)*p(1297) + p(1)*dpdr(i,1297)  
      dpdr(i,1760) = dpdr(i,1)*p(1080) + p(1)*dpdr(i,1080)  
      dpdr(i,1761) = dpdr(i,34)*p(234) + p(34)*dpdr(i,234)  
      dpdr(i,1762) = dpdr(i,1)*p(1298) + p(1)*dpdr(i,1298)  
      dpdr(i,1763) = dpdr(i,34)*p(235) + p(34)*dpdr(i,235)  
      dpdr(i,1764) = dpdr(i,1)*p(1300) + p(1)*dpdr(i,1300)  
      dpdr(i,1765) = dpdr(i,1)*p(1301) + p(1)*dpdr(i,1301)  
      dpdr(i,1766) = dpdr(i,34)*p(273) + p(34)*dpdr(i,273)  
      dpdr(i,1767) = dpdr(i,1)*p(1082) + p(1)*dpdr(i,1082)  
      dpdr(i,1768) = dpdr(i,1)*p(1305) + p(1)*dpdr(i,1305)  
      dpdr(i,1769) = dpdr(i,1)*p(1307) + p(1)*dpdr(i,1307)  
      dpdr(i,1770) = dpdr(i,1)*p(1083) + p(1)*dpdr(i,1083)  
      dpdr(i,1771) = dpdr(i,1)*p(1309) + p(1)*dpdr(i,1309)  
      dpdr(i,1772) = dpdr(i,1)*p(1084) + p(1)*dpdr(i,1084)  
      dpdr(i,1773) = dpdr(i,34)*p(275) + p(34)*dpdr(i,275)  
      dpdr(i,1774) = dpdr(i,1)*p(1310) + p(1)*dpdr(i,1310)  
      dpdr(i,1775) = dpdr(i,34)*p(276) + p(34)*dpdr(i,276)  
      dpdr(i,1776) = dpdr(i,1)*p(1312) + p(1)*dpdr(i,1312)  
      dpdr(i,1777) = dpdr(i,1)*p(1313) + p(1)*dpdr(i,1313)  
      dpdr(i,1778) = dpdr(i,34)*p(277) + p(34)*dpdr(i,277)  
      dpdr(i,1779) = dpdr(i,1)*p(1316) + p(1)*dpdr(i,1316)  
      dpdr(i,1780) = dpdr(i,1)*p(1318) + p(1)*dpdr(i,1318)  
      dpdr(i,1781) = dpdr(i,1)*p(1319) + p(1)*dpdr(i,1319)  
      dpdr(i,1782) = dpdr(i,34)*p(289) + p(34)*dpdr(i,289)  
      dpdr(i,1783) = dpdr(i,1)*p(1086) + p(1)*dpdr(i,1086)  
      dpdr(i,1784) = dpdr(i,1)*p(1087) + p(1)*dpdr(i,1087)  
      dpdr(i,1785) = dpdr(i,1)*p(1088) + p(1)*dpdr(i,1088)  
      dpdr(i,1786) = dpdr(i,1)*p(1089) + p(1)*dpdr(i,1089)  
      dpdr(i,1787) = dpdr(i,1)*p(1090) + p(1)*dpdr(i,1090)  
      dpdr(i,1788) = dpdr(i,1)*p(1091) + p(1)*dpdr(i,1091)  
      dpdr(i,1789) = dpdr(i,1)*p(1092) + p(1)*dpdr(i,1092)  
      dpdr(i,1790) = dpdr(i,1)*p(1093) + p(1)*dpdr(i,1093)  
      dpdr(i,1791) = dpdr(i,1)*p(1094) + p(1)*dpdr(i,1094)  
      dpdr(i,1792) = dpdr(i,1)*p(1095) + p(1)*dpdr(i,1095)  
      dpdr(i,1793) = dpdr(i,1)*p(1096) + p(1)*dpdr(i,1096)  
      dpdr(i,1794) = dpdr(i,1)*p(1097) + p(1)*dpdr(i,1097)  
      dpdr(i,1795) = dpdr(i,1)*p(1098) + p(1)*dpdr(i,1098)  
      dpdr(i,1796) = dpdr(i,1)*p(1099) + p(1)*dpdr(i,1099)  
      dpdr(i,1797) = dpdr(i,1)*p(1100) + p(1)*dpdr(i,1100)  
      dpdr(i,1798) = dpdr(i,1)*p(1101) + p(1)*dpdr(i,1101)  
      dpdr(i,1799) = dpdr(i,1)*p(1102) + p(1)*dpdr(i,1102)  
      dpdr(i,1800) = dpdr(i,1)*p(1105) + p(1)*dpdr(i,1105)  
      dpdr(i,1801) = dpdr(i,3)*p(1262) + p(3)*dpdr(i,1262)  
      dpdr(i,1802) = dpdr(i,1)*p(1107) + p(1)*dpdr(i,1107)  
      dpdr(i,1803) = dpdr(i,1)*p(1108) + p(1)*dpdr(i,1108)  
      dpdr(i,1804) = dpdr(i,1)*p(1109) + p(1)*dpdr(i,1109)  
      dpdr(i,1805) = dpdr(i,1)*p(1324) + p(1)*dpdr(i,1324)  
      dpdr(i,1806) = dpdr(i,1)*p(1110) + p(1)*dpdr(i,1110)  
      dpdr(i,1807) = dpdr(i,1)*p(1112) + p(1)*dpdr(i,1112)  
      dpdr(i,1808) = dpdr(i,3)*p(1269) + p(3)*dpdr(i,1269)  
      dpdr(i,1809) = dpdr(i,1)*p(1114) + p(1)*dpdr(i,1114)  
      dpdr(i,1810) = dpdr(i,3)*p(1271) + p(3)*dpdr(i,1271)  
      dpdr(i,1811) = dpdr(i,1)*p(1325) + p(1)*dpdr(i,1325)  
      dpdr(i,1812) = dpdr(i,1)*p(1116) + p(1)*dpdr(i,1116)  
      dpdr(i,1813) = dpdr(i,3)*p(1274) + p(3)*dpdr(i,1274)  
      dpdr(i,1814) = dpdr(i,3)*p(1275) + p(3)*dpdr(i,1275)  
      dpdr(i,1815) = dpdr(i,1)*p(1119) + p(1)*dpdr(i,1119)  
      dpdr(i,1816) = dpdr(i,3)*p(1277) + p(3)*dpdr(i,1277)  
      dpdr(i,1817) = dpdr(i,3)*p(1278) + p(3)*dpdr(i,1278)  
      dpdr(i,1818) = dpdr(i,1)*p(1122) + p(1)*dpdr(i,1122)  
      dpdr(i,1819) = dpdr(i,3)*p(1280) + p(3)*dpdr(i,1280)  
      dpdr(i,1820) = dpdr(i,3)*p(1281) + p(3)*dpdr(i,1281)  
      dpdr(i,1821) = dpdr(i,1)*p(1326) + p(1)*dpdr(i,1326)  
      dpdr(i,1822) = dpdr(i,3)*p(1283) + p(3)*dpdr(i,1283)  
      dpdr(i,1823) = dpdr(i,1)*p(1125) + p(1)*dpdr(i,1125)  
      dpdr(i,1824) = dpdr(i,1)*p(1126) + p(1)*dpdr(i,1126)  
      dpdr(i,1825) = dpdr(i,1)*p(1328) + p(1)*dpdr(i,1328)  
      dpdr(i,1826) = dpdr(i,1)*p(1127) + p(1)*dpdr(i,1127)  
      dpdr(i,1827) = dpdr(i,3)*p(1288) + p(3)*dpdr(i,1288)  
      dpdr(i,1828) = dpdr(i,1)*p(1129) + p(1)*dpdr(i,1129)  
      dpdr(i,1829) = dpdr(i,3)*p(1290) + p(3)*dpdr(i,1290)  
      dpdr(i,1830) = dpdr(i,1)*p(1329) + p(1)*dpdr(i,1329)  
      dpdr(i,1831) = dpdr(i,3)*p(1292) + p(3)*dpdr(i,1292)  
      dpdr(i,1832) = dpdr(i,1)*p(1131) + p(1)*dpdr(i,1131)  
      dpdr(i,1833) = dpdr(i,3)*p(1294) + p(3)*dpdr(i,1294)  
      dpdr(i,1834) = dpdr(i,3)*p(1295) + p(3)*dpdr(i,1295)  
      dpdr(i,1835) = dpdr(i,1)*p(1134) + p(1)*dpdr(i,1134)  
      dpdr(i,1836) = dpdr(i,3)*p(1297) + p(3)*dpdr(i,1297)  
      dpdr(i,1837) = dpdr(i,3)*p(1298) + p(3)*dpdr(i,1298)  
      dpdr(i,1838) = dpdr(i,1)*p(1330) + p(1)*dpdr(i,1330)  
      dpdr(i,1839) = dpdr(i,3)*p(1300) + p(3)*dpdr(i,1300)  
      dpdr(i,1840) = dpdr(i,3)*p(1301) + p(3)*dpdr(i,1301)  
      dpdr(i,1841) = dpdr(i,1)*p(1137) + p(1)*dpdr(i,1137)  
      dpdr(i,1842) = dpdr(i,1)*p(1332) + p(1)*dpdr(i,1332)  
      dpdr(i,1843) = dpdr(i,1)*p(1138) + p(1)*dpdr(i,1138)  
      dpdr(i,1844) = dpdr(i,3)*p(1305) + p(3)*dpdr(i,1305)  
      dpdr(i,1845) = dpdr(i,1)*p(1333) + p(1)*dpdr(i,1333)  
      dpdr(i,1846) = dpdr(i,3)*p(1307) + p(3)*dpdr(i,1307)  
      dpdr(i,1847) = dpdr(i,1)*p(1140) + p(1)*dpdr(i,1140)  
      dpdr(i,1848) = dpdr(i,3)*p(1309) + p(3)*dpdr(i,1309)  
      dpdr(i,1849) = dpdr(i,3)*p(1310) + p(3)*dpdr(i,1310)  
      dpdr(i,1850) = dpdr(i,1)*p(1334) + p(1)*dpdr(i,1334)  
      dpdr(i,1851) = dpdr(i,3)*p(1312) + p(3)*dpdr(i,1312)  
      dpdr(i,1852) = dpdr(i,3)*p(1313) + p(3)*dpdr(i,1313)  
      dpdr(i,1853) = dpdr(i,1)*p(1336) + p(1)*dpdr(i,1336)  
      dpdr(i,1854) = dpdr(i,1)*p(1337) + p(1)*dpdr(i,1337)  
      dpdr(i,1855) = dpdr(i,3)*p(1316) + p(3)*dpdr(i,1316)  
      dpdr(i,1856) = dpdr(i,1)*p(1338) + p(1)*dpdr(i,1338)  
      dpdr(i,1857) = dpdr(i,3)*p(1318) + p(3)*dpdr(i,1318)  
      dpdr(i,1858) = dpdr(i,3)*p(1319) + p(3)*dpdr(i,1319)  
      dpdr(i,1859) = dpdr(i,1)*p(1143) + p(1)*dpdr(i,1143)  
      dpdr(i,1860) = dpdr(i,1)*p(1144) + p(1)*dpdr(i,1144)  
      dpdr(i,1861) = dpdr(i,1)*p(1145) + p(1)*dpdr(i,1145)  
      dpdr(i,1862) = dpdr(i,1)*p(1146) + p(1)*dpdr(i,1146)  
      dpdr(i,1863) = dpdr(i,1)*p(1147) + p(1)*dpdr(i,1147)  
      dpdr(i,1864) = dpdr(i,1)*p(1148) + p(1)*dpdr(i,1148)  
      dpdr(i,1865) = dpdr(i,1)*p(1149) + p(1)*dpdr(i,1149)  
      dpdr(i,1866) = dpdr(i,1)*p(1150) + p(1)*dpdr(i,1150)  
      dpdr(i,1867) = dpdr(i,1)*p(1151) + p(1)*dpdr(i,1151)  
      dpdr(i,1868) = dpdr(i,1)*p(1152) + p(1)*dpdr(i,1152)  
      dpdr(i,1869) = dpdr(i,1)*p(1153) + p(1)*dpdr(i,1153)  
      dpdr(i,1870) = dpdr(i,1)*p(1154) + p(1)*dpdr(i,1154)  
      dpdr(i,1871) = dpdr(i,1)*p(1155) + p(1)*dpdr(i,1155)  
      dpdr(i,1872) = dpdr(i,1)*p(1156) + p(1)*dpdr(i,1156)  
      dpdr(i,1873) = dpdr(i,1)*p(1158) + p(1)*dpdr(i,1158)  
      dpdr(i,1874) = dpdr(i,1)*p(1159) + p(1)*dpdr(i,1159)  
      dpdr(i,1875) = dpdr(i,1)*p(1162) + p(1)*dpdr(i,1162)  
      dpdr(i,1876) = dpdr(i,3)*p(1103) + p(3)*dpdr(i,1103)  
      dpdr(i,1877) = dpdr(i,3)*p(1104) + p(3)*dpdr(i,1104)  
      dpdr(i,1878) = dpdr(i,1)*p(1165) + p(1)*dpdr(i,1165)  
      dpdr(i,1879) = dpdr(i,3)*p(1106) + p(3)*dpdr(i,1106)  
      dpdr(i,1880) = dpdr(i,1)*p(1167) + p(1)*dpdr(i,1167)  
      dpdr(i,1881) = dpdr(i,1)*p(1168) + p(1)*dpdr(i,1168)  
      dpdr(i,1882) = dpdr(i,1)*p(1343) + p(1)*dpdr(i,1343)  
      dpdr(i,1883) = dpdr(i,1)*p(1169) + p(1)*dpdr(i,1169)  
      dpdr(i,1884) = dpdr(i,3)*p(1111) + p(3)*dpdr(i,1111)  
      dpdr(i,1885) = dpdr(i,1)*p(1171) + p(1)*dpdr(i,1171)  
      dpdr(i,1886) = dpdr(i,3)*p(1113) + p(3)*dpdr(i,1113)  
      dpdr(i,1887) = dpdr(i,1)*p(1344) + p(1)*dpdr(i,1344)  
      dpdr(i,1888) = dpdr(i,3)*p(1115) + p(3)*dpdr(i,1115)  
      dpdr(i,1889) = dpdr(i,1)*p(1173) + p(1)*dpdr(i,1173)  
      dpdr(i,1890) = dpdr(i,3)*p(1117) + p(3)*dpdr(i,1117)  
      dpdr(i,1891) = dpdr(i,3)*p(1118) + p(3)*dpdr(i,1118)  
      dpdr(i,1892) = dpdr(i,1)*p(1176) + p(1)*dpdr(i,1176)  
      dpdr(i,1893) = dpdr(i,3)*p(1120) + p(3)*dpdr(i,1120)  
      dpdr(i,1894) = dpdr(i,3)*p(1121) + p(3)*dpdr(i,1121)  
      dpdr(i,1895) = dpdr(i,1)*p(1345) + p(1)*dpdr(i,1345)  
      dpdr(i,1896) = dpdr(i,3)*p(1123) + p(3)*dpdr(i,1123)  
      dpdr(i,1897) = dpdr(i,3)*p(1124) + p(3)*dpdr(i,1124)  
      dpdr(i,1898) = dpdr(i,1)*p(1179) + p(1)*dpdr(i,1179)  
      dpdr(i,1899) = dpdr(i,1)*p(1347) + p(1)*dpdr(i,1347)  
      dpdr(i,1900) = dpdr(i,1)*p(1180) + p(1)*dpdr(i,1180)  
      dpdr(i,1901) = dpdr(i,3)*p(1128) + p(3)*dpdr(i,1128)  
      dpdr(i,1902) = dpdr(i,1)*p(1348) + p(1)*dpdr(i,1348)  
      dpdr(i,1903) = dpdr(i,3)*p(1130) + p(3)*dpdr(i,1130)  
      dpdr(i,1904) = dpdr(i,1)*p(1182) + p(1)*dpdr(i,1182)  
      dpdr(i,1905) = dpdr(i,3)*p(1132) + p(3)*dpdr(i,1132)  
      dpdr(i,1906) = dpdr(i,3)*p(1133) + p(3)*dpdr(i,1133)  
      dpdr(i,1907) = dpdr(i,1)*p(1349) + p(1)*dpdr(i,1349)  
      dpdr(i,1908) = dpdr(i,3)*p(1135) + p(3)*dpdr(i,1135)  
      dpdr(i,1909) = dpdr(i,3)*p(1136) + p(3)*dpdr(i,1136)  
      dpdr(i,1910) = dpdr(i,1)*p(1351) + p(1)*dpdr(i,1351)  
      dpdr(i,1911) = dpdr(i,1)*p(1352) + p(1)*dpdr(i,1352)  
      dpdr(i,1912) = dpdr(i,3)*p(1139) + p(3)*dpdr(i,1139)  
      dpdr(i,1913) = dpdr(i,1)*p(1353) + p(1)*dpdr(i,1353)  
      dpdr(i,1914) = dpdr(i,3)*p(1141) + p(3)*dpdr(i,1141)  
      dpdr(i,1915) = dpdr(i,3)*p(1142) + p(3)*dpdr(i,1142)  
      dpdr(i,1916) = dpdr(i,1)*p(1185) + p(1)*dpdr(i,1185)  
      dpdr(i,1917) = dpdr(i,1)*p(1186) + p(1)*dpdr(i,1186)  
      dpdr(i,1918) = dpdr(i,1)*p(1187) + p(1)*dpdr(i,1187)  
      dpdr(i,1919) = dpdr(i,1)*p(1188) + p(1)*dpdr(i,1188)  
      dpdr(i,1920) = dpdr(i,1)*p(1189) + p(1)*dpdr(i,1189)  
      dpdr(i,1921) = dpdr(i,1)*p(1190) + p(1)*dpdr(i,1190)  
      dpdr(i,1922) = dpdr(i,1)*p(1191) + p(1)*dpdr(i,1191)  
      dpdr(i,1923) = dpdr(i,1)*p(1192) + p(1)*dpdr(i,1192)  
      dpdr(i,1924) = dpdr(i,1)*p(1193) + p(1)*dpdr(i,1193)  
      dpdr(i,1925) = dpdr(i,1)*p(1195) + p(1)*dpdr(i,1195)  
      dpdr(i,1926) = dpdr(i,1)*p(1196) + p(1)*dpdr(i,1196)  
      dpdr(i,1927) = dpdr(i,1)*p(1357) + p(1)*dpdr(i,1357)  
      dpdr(i,1928) = dpdr(i,1)*p(1197) + p(1)*dpdr(i,1197)  
      dpdr(i,1929) = dpdr(i,1)*p(1199) + p(1)*dpdr(i,1199)  
      dpdr(i,1930) = dpdr(i,3)*p(1157) + p(3)*dpdr(i,1157)  
      dpdr(i,1931) = dpdr(i,1)*p(1358) + p(1)*dpdr(i,1358)  
      dpdr(i,1932) = dpdr(i,1)*p(1201) + p(1)*dpdr(i,1201)  
      dpdr(i,1933) = dpdr(i,3)*p(1160) + p(3)*dpdr(i,1160)  
      dpdr(i,1934) = dpdr(i,3)*p(1161) + p(3)*dpdr(i,1161)  
      dpdr(i,1935) = dpdr(i,1)*p(1204) + p(1)*dpdr(i,1204)  
      dpdr(i,1936) = dpdr(i,3)*p(1163) + p(3)*dpdr(i,1163)  
      dpdr(i,1937) = dpdr(i,3)*p(1164) + p(3)*dpdr(i,1164)  
      dpdr(i,1938) = dpdr(i,1)*p(1359) + p(1)*dpdr(i,1359)  
      dpdr(i,1939) = dpdr(i,3)*p(1166) + p(3)*dpdr(i,1166)  
      dpdr(i,1940) = dpdr(i,1)*p(1207) + p(1)*dpdr(i,1207)  
      dpdr(i,1941) = dpdr(i,1)*p(1361) + p(1)*dpdr(i,1361)  
      dpdr(i,1942) = dpdr(i,1)*p(1208) + p(1)*dpdr(i,1208)  
      dpdr(i,1943) = dpdr(i,3)*p(1170) + p(3)*dpdr(i,1170)  
      dpdr(i,1944) = dpdr(i,1)*p(1362) + p(1)*dpdr(i,1362)  
      dpdr(i,1945) = dpdr(i,3)*p(1172) + p(3)*dpdr(i,1172)  
      dpdr(i,1946) = dpdr(i,1)*p(1210) + p(1)*dpdr(i,1210)  
      dpdr(i,1947) = dpdr(i,3)*p(1174) + p(3)*dpdr(i,1174)  
      dpdr(i,1948) = dpdr(i,3)*p(1175) + p(3)*dpdr(i,1175)  
      dpdr(i,1949) = dpdr(i,1)*p(1363) + p(1)*dpdr(i,1363)  
      dpdr(i,1950) = dpdr(i,3)*p(1177) + p(3)*dpdr(i,1177)  
      dpdr(i,1951) = dpdr(i,3)*p(1178) + p(3)*dpdr(i,1178)  
      dpdr(i,1952) = dpdr(i,1)*p(1365) + p(1)*dpdr(i,1365)  
      dpdr(i,1953) = dpdr(i,1)*p(1366) + p(1)*dpdr(i,1366)  
      dpdr(i,1954) = dpdr(i,3)*p(1181) + p(3)*dpdr(i,1181)  
      dpdr(i,1955) = dpdr(i,1)*p(1367) + p(1)*dpdr(i,1367)  
      dpdr(i,1956) = dpdr(i,3)*p(1183) + p(3)*dpdr(i,1183)  
      dpdr(i,1957) = dpdr(i,3)*p(1184) + p(3)*dpdr(i,1184)  
      dpdr(i,1958) = dpdr(i,1)*p(1213) + p(1)*dpdr(i,1213)  
      dpdr(i,1959) = dpdr(i,1)*p(1214) + p(1)*dpdr(i,1214)  
      dpdr(i,1960) = dpdr(i,1)*p(1215) + p(1)*dpdr(i,1215)  
      dpdr(i,1961) = dpdr(i,1)*p(1216) + p(1)*dpdr(i,1216)  
      dpdr(i,1962) = dpdr(i,1)*p(1217) + p(1)*dpdr(i,1217)  
      dpdr(i,1963) = dpdr(i,1)*p(1218) + p(1)*dpdr(i,1218)  
      dpdr(i,1964) = dpdr(i,1)*p(1219) + p(1)*dpdr(i,1219)  
      dpdr(i,1965) = dpdr(i,1)*p(1220) + p(1)*dpdr(i,1220)  
      dpdr(i,1966) = dpdr(i,1)*p(1223) + p(1)*dpdr(i,1223)  
      dpdr(i,1967) = dpdr(i,3)*p(1194) + p(3)*dpdr(i,1194)  
      dpdr(i,1968) = dpdr(i,1)*p(1225) + p(1)*dpdr(i,1225)  
      dpdr(i,1969) = dpdr(i,1)*p(1371) + p(1)*dpdr(i,1371)  
      dpdr(i,1970) = dpdr(i,1)*p(1226) + p(1)*dpdr(i,1226)  
      dpdr(i,1971) = dpdr(i,3)*p(1198) + p(3)*dpdr(i,1198)  
      dpdr(i,1972) = dpdr(i,1)*p(1372) + p(1)*dpdr(i,1372)  
      dpdr(i,1973) = dpdr(i,3)*p(1200) + p(3)*dpdr(i,1200)  
      dpdr(i,1974) = dpdr(i,1)*p(1228) + p(1)*dpdr(i,1228)  
      dpdr(i,1975) = dpdr(i,3)*p(1202) + p(3)*dpdr(i,1202)  
      dpdr(i,1976) = dpdr(i,3)*p(1203) + p(3)*dpdr(i,1203)  
      dpdr(i,1977) = dpdr(i,1)*p(1373) + p(1)*dpdr(i,1373)  
      dpdr(i,1978) = dpdr(i,3)*p(1205) + p(3)*dpdr(i,1205)  
      dpdr(i,1979) = dpdr(i,3)*p(1206) + p(3)*dpdr(i,1206)  
      dpdr(i,1980) = dpdr(i,1)*p(1375) + p(1)*dpdr(i,1375)  
      dpdr(i,1981) = dpdr(i,1)*p(1376) + p(1)*dpdr(i,1376)  
      dpdr(i,1982) = dpdr(i,3)*p(1209) + p(3)*dpdr(i,1209)  
      dpdr(i,1983) = dpdr(i,1)*p(1377) + p(1)*dpdr(i,1377)  
      dpdr(i,1984) = dpdr(i,3)*p(1211) + p(3)*dpdr(i,1211)  
      dpdr(i,1985) = dpdr(i,3)*p(1212) + p(3)*dpdr(i,1212)  
      dpdr(i,1986) = dpdr(i,1)*p(1231) + p(1)*dpdr(i,1231)  
      dpdr(i,1987) = dpdr(i,1)*p(1232) + p(1)*dpdr(i,1232)  
      dpdr(i,1988) = dpdr(i,1)*p(1233) + p(1)*dpdr(i,1233)  
      dpdr(i,1989) = dpdr(i,1)*p(1234) + p(1)*dpdr(i,1234)  
      dpdr(i,1990) = dpdr(i,1)*p(1380) + p(1)*dpdr(i,1380)  
      dpdr(i,1991) = dpdr(i,1)*p(1235) + p(1)*dpdr(i,1235)  
      dpdr(i,1992) = dpdr(i,1)*p(1381) + p(1)*dpdr(i,1381)  
      dpdr(i,1993) = dpdr(i,1)*p(1237) + p(1)*dpdr(i,1237)  
      dpdr(i,1994) = dpdr(i,3)*p(1221) + p(3)*dpdr(i,1221)  
      dpdr(i,1995) = dpdr(i,3)*p(1222) + p(3)*dpdr(i,1222)  
      dpdr(i,1996) = dpdr(i,1)*p(1382) + p(1)*dpdr(i,1382)  
      dpdr(i,1997) = dpdr(i,3)*p(1224) + p(3)*dpdr(i,1224)  
      dpdr(i,1998) = dpdr(i,1)*p(1384) + p(1)*dpdr(i,1384)  
      dpdr(i,1999) = dpdr(i,1)*p(1385) + p(1)*dpdr(i,1385)  
      dpdr(i,2000) = dpdr(i,3)*p(1227) + p(3)*dpdr(i,1227)  
      dpdr(i,2001) = dpdr(i,1)*p(1386) + p(1)*dpdr(i,1386)  
      dpdr(i,2002) = dpdr(i,3)*p(1229) + p(3)*dpdr(i,1229)  
      dpdr(i,2003) = dpdr(i,3)*p(1230) + p(3)*dpdr(i,1230)  
      dpdr(i,2004) = dpdr(i,1)*p(1240) + p(1)*dpdr(i,1240)  
      dpdr(i,2005) = dpdr(i,1)*p(1241) + p(1)*dpdr(i,1241)  
      dpdr(i,2006) = dpdr(i,1)*p(1242) + p(1)*dpdr(i,1242)  
      dpdr(i,2007) = dpdr(i,1)*p(1389) + p(1)*dpdr(i,1389)  
      dpdr(i,2008) = dpdr(i,1)*p(1390) + p(1)*dpdr(i,1390)  
      dpdr(i,2009) = dpdr(i,3)*p(1236) + p(3)*dpdr(i,1236)  
      dpdr(i,2010) = dpdr(i,1)*p(1391) + p(1)*dpdr(i,1391)  
      dpdr(i,2011) = dpdr(i,3)*p(1238) + p(3)*dpdr(i,1238)  
      dpdr(i,2012) = dpdr(i,3)*p(1239) + p(3)*dpdr(i,1239)  
      dpdr(i,2013) = dpdr(i,1)*p(1393) + p(1)*dpdr(i,1393)  
      dpdr(i,2014) = dpdr(i,1)*p(1394) + p(1)*dpdr(i,1394)  
      dpdr(i,2015) = dpdr(i,1)*p(1395) + p(1)*dpdr(i,1395)  
      dpdr(i,2016) = dpdr(i,3)*p(1243) + p(3)*dpdr(i,1243)  
      dpdr(i,2017) = dpdr(i,1)*p(1244) + p(1)*dpdr(i,1244)  
      dpdr(i,2018) = dpdr(i,1)*p(1245) + p(1)*dpdr(i,1245)  
      dpdr(i,2019) = dpdr(i,1)*p(1246) + p(1)*dpdr(i,1246)  
      dpdr(i,2020) = dpdr(i,1)*p(1247) + p(1)*dpdr(i,1247)  
      dpdr(i,2021) = dpdr(i,1)*p(1248) + p(1)*dpdr(i,1248)  
      dpdr(i,2022) = dpdr(i,1)*p(1249) + p(1)*dpdr(i,1249)  
      dpdr(i,2023) = dpdr(i,1)*p(1250) + p(1)*dpdr(i,1250)  
      dpdr(i,2024) = dpdr(i,1)*p(1251) + p(1)*dpdr(i,1251)  
      dpdr(i,2025) = dpdr(i,1)*p(1252) + p(1)*dpdr(i,1252)  
      dpdr(i,2026) = dpdr(i,1)*p(1253) + p(1)*dpdr(i,1253)  
      dpdr(i,2027) = dpdr(i,1)*p(1254) + p(1)*dpdr(i,1254)  
      dpdr(i,2028) = dpdr(i,1)*p(1255) + p(1)*dpdr(i,1255)  
      dpdr(i,2029) = dpdr(i,1)*p(1256) + p(1)*dpdr(i,1256)  
      dpdr(i,2030) = dpdr(i,1)*p(1257) + p(1)*dpdr(i,1257)  
      dpdr(i,2031) = dpdr(i,1)*p(1258) + p(1)*dpdr(i,1258)  
      dpdr(i,2032) = dpdr(i,1)*p(1259) + p(1)*dpdr(i,1259)  
      dpdr(i,2033) = dpdr(i,1)*p(1260) + p(1)*dpdr(i,1260)  
      dpdr(i,2034) = dpdr(i,1)*p(1261) + p(1)*dpdr(i,1261)  
      dpdr(i,2035) = dpdr(i,1)*p(1263) + p(1)*dpdr(i,1263)  
      dpdr(i,2036) = dpdr(i,1)*p(1264) + p(1)*dpdr(i,1264)  
      dpdr(i,2037) = dpdr(i,1)*p(1265) + p(1)*dpdr(i,1265)  
      dpdr(i,2038) = dpdr(i,1)*p(1266) + p(1)*dpdr(i,1266)  
      dpdr(i,2039) = dpdr(i,1)*p(1267) + p(1)*dpdr(i,1267)  
      dpdr(i,2040) = dpdr(i,1)*p(1268) + p(1)*dpdr(i,1268)  
      dpdr(i,2041) = dpdr(i,1)*p(1270) + p(1)*dpdr(i,1270)  
      dpdr(i,2042) = dpdr(i,5)*p(750) + p(5)*dpdr(i,750) - dpdr(i,1734) 
      dpdr(i,2043) = dpdr(i,1)*p(1272) + p(1)*dpdr(i,1272)   
      dpdr(i,2044) = dpdr(i,1)*p(1273) + p(1)*dpdr(i,1273)   
      dpdr(i,2045) = dpdr(i,1)*p(1276) + p(1)*dpdr(i,1276)   
      dpdr(i,2046) = dpdr(i,5)*p(749) + p(5)*dpdr(i,749) - dpdr(i,1705) 
      dpdr(i,2047) = dpdr(i,5)*p(764) + p(5)*dpdr(i,764) - dpdr(i,1757) 
      dpdr(i,2048) = dpdr(i,1)*p(1279) + p(1)*dpdr(i,1279)   
      dpdr(i,2049) = dpdr(i,5)*p(752) + p(5)*dpdr(i,752) - dpdr(i,1740) 
      dpdr(i,2050) = dpdr(i,5)*p(767) + p(5)*dpdr(i,767) - dpdr(i,1763) 
      dpdr(i,2051) = dpdr(i,1)*p(1282) + p(1)*dpdr(i,1282)   
      dpdr(i,2052) = dpdr(i,5)*p(770) + p(5)*dpdr(i,770) - dpdr(i,1766) 
      dpdr(i,2053) = dpdr(i,1)*p(1284) + p(1)*dpdr(i,1284)  
      dpdr(i,2054) = dpdr(i,1)*p(1285) + p(1)*dpdr(i,1285)  
      dpdr(i,2055) = dpdr(i,1)*p(1286) + p(1)*dpdr(i,1286)  
      dpdr(i,2056) = dpdr(i,1)*p(1401) + p(1)*dpdr(i,1401)  
      dpdr(i,2057) = dpdr(i,1)*p(1287) + p(1)*dpdr(i,1287)  
      dpdr(i,2058) = dpdr(i,6)*p(850) + p(6)*dpdr(i,850)   
      dpdr(i,2059) = dpdr(i,1)*p(1289) + p(1)*dpdr(i,1289)  
      dpdr(i,2060) = dpdr(i,60)*p(267) + p(60)*dpdr(i,267)  
      dpdr(i,2061) = dpdr(i,1)*p(1291) + p(1)*dpdr(i,1291)  
      dpdr(i,2062) = dpdr(i,59)*p(268) + p(59)*dpdr(i,268)  
      dpdr(i,2063) = dpdr(i,1)*p(1402) + p(1)*dpdr(i,1402)  
      dpdr(i,2064) = dpdr(i,5)*p(851) + p(5)*dpdr(i,851)   
      dpdr(i,2065) = dpdr(i,1)*p(1293) + p(1)*dpdr(i,1293)  
      dpdr(i,2066) = dpdr(i,7)*p(850) + p(7)*dpdr(i,850)   
      dpdr(i,2067) = dpdr(i,7)*p(851) + p(7)*dpdr(i,851)   
      dpdr(i,2068) = dpdr(i,1)*p(1296) + p(1)*dpdr(i,1296)  
      dpdr(i,2069) = dpdr(i,61)*p(267) + p(61)*dpdr(i,267)  
      dpdr(i,2070) = dpdr(i,61)*p(268) + p(61)*dpdr(i,268)  
      dpdr(i,2071) = dpdr(i,1)*p(1299) + p(1)*dpdr(i,1299)  
      dpdr(i,2072) = dpdr(i,59)*p(269) + p(59)*dpdr(i,269)  
      dpdr(i,2073) = dpdr(i,60)*p(269) + p(60)*dpdr(i,269)  
      dpdr(i,2074) = dpdr(i,1)*p(1403) + p(1)*dpdr(i,1403)  
      dpdr(i,2075) = dpdr(i,5)*p(852) + p(5)*dpdr(i,852)   
      dpdr(i,2076) = dpdr(i,6)*p(852) + p(6)*dpdr(i,852)   
      dpdr(i,2077) = dpdr(i,1)*p(1302) + p(1)*dpdr(i,1302)  
      dpdr(i,2078) = dpdr(i,1)*p(1303) + p(1)*dpdr(i,1303)  
      dpdr(i,2079) = dpdr(i,1)*p(1405) + p(1)*dpdr(i,1405)  
      dpdr(i,2080) = dpdr(i,1)*p(1304) + p(1)*dpdr(i,1304)  
      dpdr(i,2081) = dpdr(i,5)*p(774) + p(5)*dpdr(i,774) - dpdr(i,1755) 
      dpdr(i,2082) = dpdr(i,1)*p(1306) + p(1)*dpdr(i,1306)   
      dpdr(i,2083) = dpdr(i,5)*p(776) + p(5)*dpdr(i,776) - dpdr(i,1757) 
      dpdr(i,2084) = dpdr(i,1)*p(1406) + p(1)*dpdr(i,1406)   
      dpdr(i,2085) = dpdr(i,5)*p(855) + p(5)*dpdr(i,855) - dpdr(i,2067) 
      dpdr(i,2086) = dpdr(i,1)*p(1308) + p(1)*dpdr(i,1308)   
      dpdr(i,2087) = dpdr(i,5)*p(778) + p(5)*dpdr(i,778) - dpdr(i,1761) 
      dpdr(i,2088) = dpdr(i,6)*p(779) + p(6)*dpdr(i,779) - dpdr(i,1763) 
      dpdr(i,2089) = dpdr(i,1)*p(1311) + p(1)*dpdr(i,1311)   
      dpdr(i,2090) = dpdr(i,5)*p(781) + p(5)*dpdr(i,781) - dpdr(i,1766) 
      dpdr(i,2091) = dpdr(i,5)*p(788) + p(5)*dpdr(i,788) - dpdr(i,1782) 
      dpdr(i,2092) = dpdr(i,1)*p(1407) + p(1)*dpdr(i,1407)   
      dpdr(i,2093) = dpdr(i,5)*p(856) + p(5)*dpdr(i,856) - dpdr(i,2076) 
      dpdr(i,2094) = dpdr(i,6)*p(856) + p(6)*dpdr(i,856) - dpdr(i,2075) 
      dpdr(i,2095) = dpdr(i,1)*p(1314) + p(1)*dpdr(i,1314)   
      dpdr(i,2096) = dpdr(i,1)*p(1409) + p(1)*dpdr(i,1409)   
      dpdr(i,2097) = dpdr(i,1)*p(1315) + p(1)*dpdr(i,1315)   
      dpdr(i,2098) = dpdr(i,5)*p(785) + p(5)*dpdr(i,785) - dpdr(i,1775) 
      dpdr(i,2099) = dpdr(i,1)*p(1410) + p(1)*dpdr(i,1410)   
      dpdr(i,2100) = dpdr(i,5)*p(859) + p(5)*dpdr(i,859) - dpdr(i,2088) 
      dpdr(i,2101) = dpdr(i,1)*p(1317) + p(1)*dpdr(i,1317)   
      dpdr(i,2102) = dpdr(i,5)*p(787) + p(5)*dpdr(i,787) - dpdr(i,1778) 
      dpdr(i,2103) = dpdr(i,6)*p(788) + p(6)*dpdr(i,788) - dpdr(i,1778) 
      dpdr(i,2104) = dpdr(i,1)*p(1411) + p(1)*dpdr(i,1411)   
      dpdr(i,2105) = dpdr(i,5)*p(860) + p(5)*dpdr(i,860) - dpdr(i,2094) 
      dpdr(i,2106) = dpdr(i,6)*p(860) + p(6)*dpdr(i,860) - dpdr(i,2093) 
      dpdr(i,2107) = dpdr(i,1)*p(1413) + p(1)*dpdr(i,1413)   
      dpdr(i,2108) = dpdr(i,1)*p(1414) + p(1)*dpdr(i,1414)   
      dpdr(i,2109) = dpdr(i,5)*p(863) + p(5)*dpdr(i,863) - dpdr(i,2103) 
      dpdr(i,2110) = dpdr(i,1)*p(1415) + p(1)*dpdr(i,1415)   
      dpdr(i,2111) = dpdr(i,5)*p(864) + p(5)*dpdr(i,864) - dpdr(i,2106) 
      dpdr(i,2112) = dpdr(i,6)*p(864) + p(6)*dpdr(i,864) - dpdr(i,2105) 
      dpdr(i,2113) = dpdr(i,1)*p(1320) + p(1)*dpdr(i,1320)  
      dpdr(i,2114) = dpdr(i,1)*p(1321) + p(1)*dpdr(i,1321)  
      dpdr(i,2115) = dpdr(i,1)*p(1322) + p(1)*dpdr(i,1322)  
      dpdr(i,2116) = dpdr(i,1)*p(1323) + p(1)*dpdr(i,1323)  
      dpdr(i,2117) = dpdr(i,1)*p(1327) + p(1)*dpdr(i,1327)  
      dpdr(i,2118) = dpdr(i,3)*p(1401) + p(3)*dpdr(i,1401)  
      dpdr(i,2119) = dpdr(i,3)*p(1402) + p(3)*dpdr(i,1402)  
      dpdr(i,2120) = dpdr(i,3)*p(1403) + p(3)*dpdr(i,1403)  
      dpdr(i,2121) = dpdr(i,1)*p(1331) + p(1)*dpdr(i,1331)  
      dpdr(i,2122) = dpdr(i,3)*p(1405) + p(3)*dpdr(i,1405)  
      dpdr(i,2123) = dpdr(i,3)*p(1406) + p(3)*dpdr(i,1406)  
      dpdr(i,2124) = dpdr(i,3)*p(1407) + p(3)*dpdr(i,1407)  
      dpdr(i,2125) = dpdr(i,1)*p(1335) + p(1)*dpdr(i,1335)  
      dpdr(i,2126) = dpdr(i,3)*p(1409) + p(3)*dpdr(i,1409)  
      dpdr(i,2127) = dpdr(i,3)*p(1410) + p(3)*dpdr(i,1410)  
      dpdr(i,2128) = dpdr(i,3)*p(1411) + p(3)*dpdr(i,1411)  
      dpdr(i,2129) = dpdr(i,1)*p(1417) + p(1)*dpdr(i,1417)  
      dpdr(i,2130) = dpdr(i,3)*p(1413) + p(3)*dpdr(i,1413)  
      dpdr(i,2131) = dpdr(i,3)*p(1414) + p(3)*dpdr(i,1414)  
      dpdr(i,2132) = dpdr(i,3)*p(1415) + p(3)*dpdr(i,1415)  
      dpdr(i,2133) = dpdr(i,1)*p(1339) + p(1)*dpdr(i,1339)  
      dpdr(i,2134) = dpdr(i,1)*p(1340) + p(1)*dpdr(i,1340)  
      dpdr(i,2135) = dpdr(i,1)*p(1341) + p(1)*dpdr(i,1341)  
      dpdr(i,2136) = dpdr(i,1)*p(1342) + p(1)*dpdr(i,1342)  
      dpdr(i,2137) = dpdr(i,3)*p(1324) + p(3)*dpdr(i,1324)  
      dpdr(i,2138) = dpdr(i,3)*p(1325) + p(3)*dpdr(i,1325)  
      dpdr(i,2139) = dpdr(i,3)*p(1326) + p(3)*dpdr(i,1326)  
      dpdr(i,2140) = dpdr(i,1)*p(1346) + p(1)*dpdr(i,1346)  
      dpdr(i,2141) = dpdr(i,3)*p(1328) + p(3)*dpdr(i,1328)  
      dpdr(i,2142) = dpdr(i,3)*p(1329) + p(3)*dpdr(i,1329)  
      dpdr(i,2143) = dpdr(i,3)*p(1330) + p(3)*dpdr(i,1330)  
      dpdr(i,2144) = dpdr(i,1)*p(1350) + p(1)*dpdr(i,1350)  
      dpdr(i,2145) = dpdr(i,3)*p(1332) + p(3)*dpdr(i,1332)  
      dpdr(i,2146) = dpdr(i,3)*p(1333) + p(3)*dpdr(i,1333)  
      dpdr(i,2147) = dpdr(i,3)*p(1334) + p(3)*dpdr(i,1334)  
      dpdr(i,2148) = dpdr(i,1)*p(1419) + p(1)*dpdr(i,1419)  
      dpdr(i,2149) = dpdr(i,3)*p(1336) + p(3)*dpdr(i,1336)  
      dpdr(i,2150) = dpdr(i,3)*p(1337) + p(3)*dpdr(i,1337)  
      dpdr(i,2151) = dpdr(i,3)*p(1338) + p(3)*dpdr(i,1338)  
      dpdr(i,2152) = dpdr(i,1)*p(1354) + p(1)*dpdr(i,1354)  
      dpdr(i,2153) = dpdr(i,1)*p(1355) + p(1)*dpdr(i,1355)  
      dpdr(i,2154) = dpdr(i,1)*p(1356) + p(1)*dpdr(i,1356)  
      dpdr(i,2155) = dpdr(i,1)*p(1360) + p(1)*dpdr(i,1360)  
      dpdr(i,2156) = dpdr(i,3)*p(1343) + p(3)*dpdr(i,1343)  
      dpdr(i,2157) = dpdr(i,3)*p(1344) + p(3)*dpdr(i,1344)  
      dpdr(i,2158) = dpdr(i,3)*p(1345) + p(3)*dpdr(i,1345)  
      dpdr(i,2159) = dpdr(i,1)*p(1364) + p(1)*dpdr(i,1364)  
      dpdr(i,2160) = dpdr(i,3)*p(1347) + p(3)*dpdr(i,1347)  
      dpdr(i,2161) = dpdr(i,3)*p(1348) + p(3)*dpdr(i,1348)  
      dpdr(i,2162) = dpdr(i,3)*p(1349) + p(3)*dpdr(i,1349)  
      dpdr(i,2163) = dpdr(i,1)*p(1421) + p(1)*dpdr(i,1421)  
      dpdr(i,2164) = dpdr(i,3)*p(1351) + p(3)*dpdr(i,1351)  
      dpdr(i,2165) = dpdr(i,3)*p(1352) + p(3)*dpdr(i,1352)  
      dpdr(i,2166) = dpdr(i,3)*p(1353) + p(3)*dpdr(i,1353)  
      dpdr(i,2167) = dpdr(i,1)*p(1368) + p(1)*dpdr(i,1368)  
      dpdr(i,2168) = dpdr(i,1)*p(1369) + p(1)*dpdr(i,1369)  
      dpdr(i,2169) = dpdr(i,1)*p(1370) + p(1)*dpdr(i,1370)  
      dpdr(i,2170) = dpdr(i,3)*p(1357) + p(3)*dpdr(i,1357)  
      dpdr(i,2171) = dpdr(i,3)*p(1358) + p(3)*dpdr(i,1358)  
      dpdr(i,2172) = dpdr(i,3)*p(1359) + p(3)*dpdr(i,1359)  
      dpdr(i,2173) = dpdr(i,1)*p(1374) + p(1)*dpdr(i,1374)  
      dpdr(i,2174) = dpdr(i,3)*p(1361) + p(3)*dpdr(i,1361)  
      dpdr(i,2175) = dpdr(i,3)*p(1362) + p(3)*dpdr(i,1362)  
      dpdr(i,2176) = dpdr(i,3)*p(1363) + p(3)*dpdr(i,1363)  
      dpdr(i,2177) = dpdr(i,1)*p(1423) + p(1)*dpdr(i,1423)  
      dpdr(i,2178) = dpdr(i,3)*p(1365) + p(3)*dpdr(i,1365)  
      dpdr(i,2179) = dpdr(i,3)*p(1366) + p(3)*dpdr(i,1366)  
      dpdr(i,2180) = dpdr(i,3)*p(1367) + p(3)*dpdr(i,1367)  
      dpdr(i,2181) = dpdr(i,1)*p(1378) + p(1)*dpdr(i,1378)  
      dpdr(i,2182) = dpdr(i,1)*p(1379) + p(1)*dpdr(i,1379)  
      dpdr(i,2183) = dpdr(i,1)*p(1383) + p(1)*dpdr(i,1383)  
      dpdr(i,2184) = dpdr(i,3)*p(1371) + p(3)*dpdr(i,1371)  
      dpdr(i,2185) = dpdr(i,3)*p(1372) + p(3)*dpdr(i,1372)  
      dpdr(i,2186) = dpdr(i,3)*p(1373) + p(3)*dpdr(i,1373)  
      dpdr(i,2187) = dpdr(i,1)*p(1425) + p(1)*dpdr(i,1425)  
      dpdr(i,2188) = dpdr(i,3)*p(1375) + p(3)*dpdr(i,1375)  
      dpdr(i,2189) = dpdr(i,3)*p(1376) + p(3)*dpdr(i,1376)  
      dpdr(i,2190) = dpdr(i,3)*p(1377) + p(3)*dpdr(i,1377)  
      dpdr(i,2191) = dpdr(i,1)*p(1387) + p(1)*dpdr(i,1387)  
      dpdr(i,2192) = dpdr(i,1)*p(1388) + p(1)*dpdr(i,1388)  
      dpdr(i,2193) = dpdr(i,3)*p(1380) + p(3)*dpdr(i,1380)  
      dpdr(i,2194) = dpdr(i,3)*p(1381) + p(3)*dpdr(i,1381)  
      dpdr(i,2195) = dpdr(i,3)*p(1382) + p(3)*dpdr(i,1382)  
      dpdr(i,2196) = dpdr(i,1)*p(1427) + p(1)*dpdr(i,1427)  
      dpdr(i,2197) = dpdr(i,3)*p(1384) + p(3)*dpdr(i,1384)  
      dpdr(i,2198) = dpdr(i,3)*p(1385) + p(3)*dpdr(i,1385)  
      dpdr(i,2199) = dpdr(i,3)*p(1386) + p(3)*dpdr(i,1386)  
      dpdr(i,2200) = dpdr(i,1)*p(1392) + p(1)*dpdr(i,1392)  
      dpdr(i,2201) = dpdr(i,1)*p(1429) + p(1)*dpdr(i,1429)  
      dpdr(i,2202) = dpdr(i,3)*p(1389) + p(3)*dpdr(i,1389)  
      dpdr(i,2203) = dpdr(i,3)*p(1390) + p(3)*dpdr(i,1390)  
      dpdr(i,2204) = dpdr(i,3)*p(1391) + p(3)*dpdr(i,1391)  
      dpdr(i,2205) = dpdr(i,1)*p(1431) + p(1)*dpdr(i,1431)  
      dpdr(i,2206) = dpdr(i,3)*p(1393) + p(3)*dpdr(i,1393)  
      dpdr(i,2207) = dpdr(i,3)*p(1394) + p(3)*dpdr(i,1394)  
      dpdr(i,2208) = dpdr(i,3)*p(1395) + p(3)*dpdr(i,1395)  
      dpdr(i,2209) = dpdr(i,1)*p(1396) + p(1)*dpdr(i,1396)  
      dpdr(i,2210) = dpdr(i,1)*p(1397) + p(1)*dpdr(i,1397)  
      dpdr(i,2211) = dpdr(i,1)*p(1398) + p(1)*dpdr(i,1398)  
      dpdr(i,2212) = dpdr(i,1)*p(1399) + p(1)*dpdr(i,1399)  
      dpdr(i,2213) = dpdr(i,1)*p(1400) + p(1)*dpdr(i,1400)  
      dpdr(i,2214) = dpdr(i,5)*p(850) + p(5)*dpdr(i,850) - dpdr(i,1719) 
      dpdr(i,2215) = dpdr(i,6)*p(851) + p(6)*dpdr(i,851) - dpdr(i,1725) 
      dpdr(i,2216) = dpdr(i,7)*p(852) + p(7)*dpdr(i,852) - dpdr(i,1742) 
      dpdr(i,2217) = dpdr(i,1)*p(1404) + p(1)*dpdr(i,1404)    
      dpdr(i,2218) = dpdr(i,5)*p(854) + p(5)*dpdr(i,854) - dpdr(i,1753) 
      dpdr(i,2219) = dpdr(i,6)*p(855) + p(6)*dpdr(i,855) - dpdr(i,1757) 
      dpdr(i,2220) = dpdr(i,7)*p(856) + p(7)*dpdr(i,856) - dpdr(i,1766) 
      dpdr(i,2221) = dpdr(i,1)*p(1408) + p(1)*dpdr(i,1408)    
      dpdr(i,2222) = dpdr(i,5)*p(858) + p(5)*dpdr(i,858) - dpdr(i,1773) 
      dpdr(i,2223) = dpdr(i,6)*p(859) + p(6)*dpdr(i,859) - dpdr(i,1775) 
      dpdr(i,2224) = dpdr(i,7)*p(860) + p(7)*dpdr(i,860) - dpdr(i,1778) 
      dpdr(i,2225) = dpdr(i,1)*p(1412) + p(1)*dpdr(i,1412)    
      dpdr(i,2226) = dpdr(i,5)*p(862) + p(5)*dpdr(i,862) - dpdr(i,1782) 
      dpdr(i,2227) = dpdr(i,6)*p(863) + p(6)*dpdr(i,863) - dpdr(i,1782) 
      dpdr(i,2228) = dpdr(i,7)*p(864) + p(7)*dpdr(i,864) - dpdr(i,1782) 
      dpdr(i,2229) = dpdr(i,1)*p(1433) + p(1)*dpdr(i,1433)    
      dpdr(i,2230) = dpdr(i,5)*p(880) + p(5)*dpdr(i,880) - dpdr(i,2112) 
      dpdr(i,2231) = dpdr(i,6)*p(880) + p(6)*dpdr(i,880) - dpdr(i,2111) 
      dpdr(i,2232) = dpdr(i,7)*p(880) + p(7)*dpdr(i,880) - dpdr(i,2109) 
      dpdr(i,2233) = dpdr(i,1)*p(1416) + p(1)*dpdr(i,1416)  
      dpdr(i,2234) = dpdr(i,3)*p(1433) + p(3)*dpdr(i,1433)  
      dpdr(i,2235) = dpdr(i,1)*p(1418) + p(1)*dpdr(i,1418)  
      dpdr(i,2236) = dpdr(i,3)*p(1417) + p(3)*dpdr(i,1417)  
      dpdr(i,2237) = dpdr(i,1)*p(1420) + p(1)*dpdr(i,1420)  
      dpdr(i,2238) = dpdr(i,3)*p(1419) + p(3)*dpdr(i,1419)  
      dpdr(i,2239) = dpdr(i,1)*p(1422) + p(1)*dpdr(i,1422)  
      dpdr(i,2240) = dpdr(i,3)*p(1421) + p(3)*dpdr(i,1421)  
      dpdr(i,2241) = dpdr(i,1)*p(1424) + p(1)*dpdr(i,1424)  
      dpdr(i,2242) = dpdr(i,3)*p(1423) + p(3)*dpdr(i,1423)  
      dpdr(i,2243) = dpdr(i,1)*p(1426) + p(1)*dpdr(i,1426)  
      dpdr(i,2244) = dpdr(i,3)*p(1425) + p(3)*dpdr(i,1425)  
      dpdr(i,2245) = dpdr(i,1)*p(1428) + p(1)*dpdr(i,1428)  
      dpdr(i,2246) = dpdr(i,3)*p(1427) + p(3)*dpdr(i,1427)  
      dpdr(i,2247) = dpdr(i,1)*p(1430) + p(1)*dpdr(i,1430)  
      dpdr(i,2248) = dpdr(i,3)*p(1429) + p(3)*dpdr(i,1429)  
      dpdr(i,2249) = dpdr(i,1)*p(1434) + p(1)*dpdr(i,1434)  
      dpdr(i,2250) = dpdr(i,3)*p(1431) + p(3)*dpdr(i,1431)  
      dpdr(i,2251) = dpdr(i,1)*p(1432) + p(1)*dpdr(i,1432)  
      dpdr(i,2252) = dpdr(i,2)*p(1433) + p(2)*dpdr(i,1433) &
                   - dpdr(i,2232) - dpdr(i,2231) - dpdr(i,2230)  
      dpdr(i,2253) = dpdr(i,3)*p(1434) + p(3)*dpdr(i,1434)
      enddo

      return
      end subroutine EvdPdr

      subroutine evdbdr
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine eliminates the 2-body terms in the permutationally
!  invariant formulation.
!**********************************************************************
      
      integer i,j
      double precision db1dr(6,2254) 

! Pass P(0:2253) to BM1(1:2254)
      do j=1,6
      do i=1,2254
        db1dr(j,i)=dpdr(j,i-1)
      enddo
      enddo

! Remove unconnected terms and 2-body terms and pass to B(1:2153)
      do j=1,6 

      db2dr(j,1)=db1dr(j,5)
 
      do i=2,3
        db2dr(j,i)=db1dr(j,i+5)
      enddo
 
      db2dr(j,4)=db1dr(j,10)
 
      do i=5,14
        db2dr(j,i)=db1dr(j,i+9)
      enddo

      do i=15,16  
        db2dr(j,i)=db1dr(j,i+10)
      enddo

      db2dr(j,17)=db1dr(j,28) 

      db2dr(j,18)=db1dr(j,30)
     
      do i=19,44
        db2dr(j,i)=db1dr(j,i+15)
      enddo

      do i=45,47
        db2dr(j,i)=db1dr(j,i+16)  
      enddo

      do i=48,49
        db2dr(j,i)=db1dr(j,i+17)
      enddo

      db2dr(j,50)=db1dr(j,68)

      db2dr(j,51)=db1dr(j,70)

      db2dr(j,52)=db1dr(j,72)

      do i=53,112
        db2dr(j,i)=db1dr(j,i+23)
      enddo

      do i=113,115
        db2dr(j,i)=db1dr(j,i+24)     
      enddo

      do i=116,117
        db2dr(j,i)=db1dr(j,i+25)
      enddo

      db2dr(j,118)=db1dr(j,144)

      db2dr(j,119)=db1dr(j,146)

      db2dr(j,120)=db1dr(j,148)

      db2dr(j,121)=db1dr(j,150)

      do i=122,235
        db2dr(j,i)=db1dr(j,i+32)
      enddo

      do i=236,238
        db2dr(j,i)=db1dr(j,i+33)
      enddo

      do i=239,241
        db2dr(j,i)=db1dr(j,i+34)
      enddo

      do i=242,243
        db2dr(j,i)=db1dr(j,i+35)
      enddo

      db2dr(j,244)=db1dr(j,280)

      db2dr(j,245)=db1dr(j,282)

      db2dr(j,246)=db1dr(j,284)

      db2dr(j,247)=db1dr(j,286)

      db2dr(j,248)=db1dr(j,288)

      do i=249,450
        db2dr(j,i)=db1dr(j,i+43)
      enddo

      do i=451,453
        db2dr(j,i)=db1dr(j,i+44)
      enddo

      do i=454,456
        db2dr(j,i)=db1dr(j,i+45)
      enddo

      do i=457,458
        db2dr(j,i)=db1dr(j,i+46)
      enddo

      db2dr(j,459)=db1dr(j,506)

      db2dr(j,460)=db1dr(j,508)

      db2dr(j,461)=db1dr(j,510)

      db2dr(j,462)=db1dr(j,512)

      db2dr(j,463)=db1dr(j,514)

      db2dr(j,464)=db1dr(j,516)

      do i=465,795
        db2dr(j,i)=db1dr(j,i+55)
      enddo

      do i=796,798
        db2dr(j,i)=db1dr(j,i+56)
      enddo

      do i=799,801
        db2dr(j,i)=db1dr(j,i+57)
      enddo

      do i=802,804
        db2dr(j,i)=db1dr(j,i+58)
      enddo

      do i=805,806
        db2dr(j,i)=db1dr(j,i+59)
      enddo

      db2dr(j,807)=db1dr(j,867)

      db2dr(j,808)=db1dr(j,869)

      db2dr(j,809)=db1dr(j,871)

      db2dr(j,810)=db1dr(j,873)

      db2dr(j,811)=db1dr(j,875)

      db2dr(j,812)=db1dr(j,877)

      db2dr(j,813)=db1dr(j,879)

      do i=814,1332
        db2dr(j,i)=db1dr(j,i+69)
      enddo

      do i=1333,1335
        db2dr(j,i)=db1dr(j,i+70)
      enddo

      do i=1336,1338
        db2dr(j,i)=db1dr(j,i+71)
      enddo

      do i=1339,1341
        db2dr(j,i)=db1dr(j,i+72)
      enddo

      do i=1342,1343
        db2dr(j,i)=db1dr(j,i+73)
      enddo

      db2dr(j,1344)=db1dr(j,1418)

      db2dr(j,1345)=db1dr(j,1420)

      db2dr(j,1346)=db1dr(j,1422)
      
      db2dr(j,1347)=db1dr(j,1424)

      db2dr(j,1348)=db1dr(j,1426)

      db2dr(j,1349)=db1dr(j,1428)

      db2dr(j,1350)=db1dr(j,1430)

      db2dr(j,1351)=db1dr(j,1432)

      do i=1352,2130
        db2dr(j,i)=db1dr(j,i+84)
      enddo

      do i=2131,2133
        db2dr(j,i)=db1dr(j,i+85)
      enddo

      do i=2134,2136
        db2dr(j,i)=db1dr(j,i+86)
      enddo

      do i=2137,2139
        db2dr(j,i)=db1dr(j,i+87)
      enddo

      do i=2140,2142
        db2dr(j,i)=db1dr(j,i+88)
      enddo

      do i=2143,2144
        db2dr(j,i)=db1dr(j,i+89)
      enddo

      db2dr(j,2145)=db1dr(j,2235)

      db2dr(j,2146)=db1dr(j,2237)

      db2dr(j,2147)=db1dr(j,2239)

      db2dr(j,2148)=db1dr(j,2241)

      db2dr(j,2149)=db1dr(j,2243)

      db2dr(j,2150)=db1dr(j,2245)

      db2dr(j,2151)=db1dr(j,2247)

      db2dr(j,2152)=db1dr(j,2249)

      db2dr(j,2153)=db1dr(j,2251)

      enddo
          
      return
      end subroutine evdbdr

      subroutine evdbdr2
      use N2O_3Ap_ZV_par
!**********************************************************************
!  This subroutine eliminates the unused four-body terms
!**********************************************************************
      
      integer j

      do j=1,6
      dbdr(j,  1) = db2dr(j,   1)
      dbdr(j,  2) = db2dr(j,   2)
      dbdr(j,  3) = db2dr(j,   3)
      dbdr(j,  4) = db2dr(j,   4)
      dbdr(j,  5) = db2dr(j,   5)
      dbdr(j,  6) = db2dr(j,   6)
      dbdr(j,  7) = db2dr(j,   7)
      dbdr(j,  8) = db2dr(j,   8)
      dbdr(j,  9) = db2dr(j,   9)
      dbdr(j, 10) = db2dr(j,  10)
      dbdr(j, 11) = db2dr(j,  11)
      dbdr(j, 12) = db2dr(j,  12)
      dbdr(j, 13) = db2dr(j,  13)
      dbdr(j, 14) = db2dr(j,  14)
      dbdr(j, 15) = db2dr(j,  15)
      dbdr(j, 16) = db2dr(j,  16)
      dbdr(j, 17) = db2dr(j,  17)
      dbdr(j, 18) = db2dr(j,  18)
      dbdr(j, 19) = db2dr(j,  19)
      dbdr(j, 20) = db2dr(j,  20)
      dbdr(j, 21) = db2dr(j,  21)
      dbdr(j, 22) = db2dr(j,  22)
      dbdr(j, 23) = db2dr(j,  23)
      dbdr(j, 24) = db2dr(j,  24)
      dbdr(j, 25) = db2dr(j,  25)
      dbdr(j, 26) = db2dr(j,  26)
      dbdr(j, 27) = db2dr(j,  27)
      dbdr(j, 28) = db2dr(j,  28)
      dbdr(j, 29) = db2dr(j,  29)
      dbdr(j, 30) = db2dr(j,  30)
      dbdr(j, 31) = db2dr(j,  31)
      dbdr(j, 32) = db2dr(j,  32)
      dbdr(j, 33) = db2dr(j,  33)
      dbdr(j, 34) = db2dr(j,  34)
      dbdr(j, 35) = db2dr(j,  35)
      dbdr(j, 36) = db2dr(j,  36)
      dbdr(j, 37) = db2dr(j,  37)
      dbdr(j, 38) = db2dr(j,  38)
      dbdr(j, 39) = db2dr(j,  39)
      dbdr(j, 40) = db2dr(j,  40)
      dbdr(j, 41) = db2dr(j,  41)
      dbdr(j, 42) = db2dr(j,  42)
      dbdr(j, 43) = db2dr(j,  43)
      dbdr(j, 44) = db2dr(j,  44)
      dbdr(j, 45) = db2dr(j,  45)
      dbdr(j, 46) = db2dr(j,  46)
      dbdr(j, 47) = db2dr(j,  47)
      dbdr(j, 48) = db2dr(j,  48)
      dbdr(j, 49) = db2dr(j,  49)
      dbdr(j, 50) = db2dr(j,  50)
      dbdr(j, 51) = db2dr(j,  51)
      dbdr(j, 52) = db2dr(j,  52)
      dbdr(j, 53) = db2dr(j,  53)
      dbdr(j, 54) = db2dr(j,  54)
      dbdr(j, 55) = db2dr(j,  55)
      dbdr(j, 56) = db2dr(j,  56)
      dbdr(j, 57) = db2dr(j,  57)
      dbdr(j, 58) = db2dr(j,  58)
      dbdr(j, 59) = db2dr(j,  59)
      dbdr(j, 60) = db2dr(j,  60)
      dbdr(j, 61) = db2dr(j,  61)
      dbdr(j, 62) = db2dr(j,  62)
      dbdr(j, 63) = db2dr(j,  63)
      dbdr(j, 64) = db2dr(j,  64)
      dbdr(j, 65) = db2dr(j,  65)
      dbdr(j, 66) = db2dr(j,  66)
      dbdr(j, 67) = db2dr(j,  67)
      dbdr(j, 68) = db2dr(j,  68)
      dbdr(j, 69) = db2dr(j,  69)
      dbdr(j, 70) = db2dr(j,  70)
      dbdr(j, 71) = db2dr(j,  71)
      dbdr(j, 72) = db2dr(j,  72)
      dbdr(j, 73) = db2dr(j,  73)
      dbdr(j, 74) = db2dr(j,  74)
      dbdr(j, 75) = db2dr(j,  75)
      dbdr(j, 76) = db2dr(j,  76)
      dbdr(j, 77) = db2dr(j,  77)
      dbdr(j, 78) = db2dr(j,  78)
      dbdr(j, 79) = db2dr(j,  79)
      dbdr(j, 80) = db2dr(j,  80)
      dbdr(j, 81) = db2dr(j,  81)
      dbdr(j, 82) = db2dr(j,  82)
      dbdr(j, 83) = db2dr(j,  83)
      dbdr(j, 84) = db2dr(j,  84)
      dbdr(j, 85) = db2dr(j,  85)
      dbdr(j, 86) = db2dr(j,  86)
      dbdr(j, 87) = db2dr(j,  87)
      dbdr(j, 88) = db2dr(j,  88)
      dbdr(j, 89) = db2dr(j,  89)
      dbdr(j, 90) = db2dr(j,  90)
      dbdr(j, 91) = db2dr(j,  91)
      dbdr(j, 92) = db2dr(j,  92)
      dbdr(j, 93) = db2dr(j,  93)
      dbdr(j, 94) = db2dr(j,  94)
      dbdr(j, 95) = db2dr(j,  95)
      dbdr(j, 96) = db2dr(j,  96)
      dbdr(j, 97) = db2dr(j,  97)
      dbdr(j, 98) = db2dr(j,  98)
      dbdr(j, 99) = db2dr(j,  99)
      dbdr(j,100) = db2dr(j, 100)
      dbdr(j,101) = db2dr(j, 101)
      dbdr(j,102) = db2dr(j, 102)
      dbdr(j,103) = db2dr(j, 103)
      dbdr(j,104) = db2dr(j, 104)
      dbdr(j,105) = db2dr(j, 105)
      dbdr(j,106) = db2dr(j, 106)
      dbdr(j,107) = db2dr(j, 107)
      dbdr(j,108) = db2dr(j, 108)
      dbdr(j,109) = db2dr(j, 109)
      dbdr(j,110) = db2dr(j, 110)
      dbdr(j,111) = db2dr(j, 111)
      dbdr(j,112) = db2dr(j, 112)
      dbdr(j,113) = db2dr(j, 113)
      dbdr(j,114) = db2dr(j, 114)
      dbdr(j,115) = db2dr(j, 115)
      dbdr(j,116) = db2dr(j, 116)
      dbdr(j,117) = db2dr(j, 117)
      dbdr(j,118) = db2dr(j, 118)
      dbdr(j,119) = db2dr(j, 119)
      dbdr(j,120) = db2dr(j, 120)
      dbdr(j,121) = db2dr(j, 121)
      dbdr(j,122) = db2dr(j, 122)
      dbdr(j,123) = db2dr(j, 123)
      dbdr(j,124) = db2dr(j, 124)
      dbdr(j,125) = db2dr(j, 125)
      dbdr(j,126) = db2dr(j, 126)
      dbdr(j,127) = db2dr(j, 127)
      dbdr(j,128) = db2dr(j, 128)
      dbdr(j,129) = db2dr(j, 129)
      dbdr(j,130) = db2dr(j, 130)
      dbdr(j,131) = db2dr(j, 131)
      dbdr(j,132) = db2dr(j, 132)
      dbdr(j,133) = db2dr(j, 133)
      dbdr(j,134) = db2dr(j, 134)
      dbdr(j,135) = db2dr(j, 135)
      dbdr(j,136) = db2dr(j, 136)
      dbdr(j,137) = db2dr(j, 137)
      dbdr(j,138) = db2dr(j, 138)
      dbdr(j,139) = db2dr(j, 140)
      dbdr(j,140) = db2dr(j, 179)
      dbdr(j,141) = db2dr(j, 184)
      dbdr(j,142) = db2dr(j, 185)
      dbdr(j,143) = db2dr(j, 193)
      dbdr(j,144) = db2dr(j, 196)
      dbdr(j,145) = db2dr(j, 202)
      dbdr(j,146) = db2dr(j, 233)
      dbdr(j,147) = db2dr(j, 234)
      dbdr(j,148) = db2dr(j, 235)
      dbdr(j,149) = db2dr(j, 237)
      dbdr(j,150) = db2dr(j, 238)
      dbdr(j,151) = db2dr(j, 240)
      dbdr(j,152) = db2dr(j, 241)
      dbdr(j,153) = db2dr(j, 243)
      dbdr(j,154) = db2dr(j, 365)
      dbdr(j,155) = db2dr(j, 370)
      dbdr(j,156) = db2dr(j, 371)
      dbdr(j,157) = db2dr(j, 379)
      dbdr(j,158) = db2dr(j, 382)
      dbdr(j,159) = db2dr(j, 385)
      dbdr(j,160) = db2dr(j, 393)
      dbdr(j,161) = db2dr(j, 396)
      dbdr(j,162) = db2dr(j, 402)
      dbdr(j,163) = db2dr(j, 447)
      dbdr(j,164) = db2dr(j, 448)
      dbdr(j,165) = db2dr(j, 449)
      dbdr(j,166) = db2dr(j, 450)
      dbdr(j,167) = db2dr(j, 452)
      dbdr(j,168) = db2dr(j, 453)
      dbdr(j,169) = db2dr(j, 455)
      dbdr(j,170) = db2dr(j, 456)
      dbdr(j,171) = db2dr(j, 458)
      dbdr(j,172) = db2dr(j, 680)
      dbdr(j,173) = db2dr(j, 685)
      dbdr(j,174) = db2dr(j, 686)
      dbdr(j,175) = db2dr(j, 693)
      dbdr(j,176) = db2dr(j, 694)
      dbdr(j,177) = db2dr(j, 697)
      dbdr(j,178) = db2dr(j, 708)
      dbdr(j,179) = db2dr(j, 711)
      dbdr(j,180) = db2dr(j, 714)
      dbdr(j,181) = db2dr(j, 723)
      dbdr(j,182) = db2dr(j, 726)
      dbdr(j,183) = db2dr(j, 732)
      dbdr(j,184) = db2dr(j, 792)
      dbdr(j,185) = db2dr(j, 793)
      dbdr(j,186) = db2dr(j, 794)
      dbdr(j,187) = db2dr(j, 795)
      dbdr(j,188) = db2dr(j, 797)
      dbdr(j,189) = db2dr(j, 798)
      dbdr(j,190) = db2dr(j, 800)
      dbdr(j,191) = db2dr(j, 801)
      dbdr(j,192) = db2dr(j, 803)
      dbdr(j,193) = db2dr(j, 804)
      dbdr(j,194) = db2dr(j, 806)
      dbdr(j,195) = db2dr(j,1178)
      dbdr(j,196) = db2dr(j,1183)
      dbdr(j,197) = db2dr(j,1184)
      dbdr(j,198) = db2dr(j,1191)
      dbdr(j,199) = db2dr(j,1192)
      dbdr(j,200) = db2dr(j,1193)
      dbdr(j,201) = db2dr(j,1205)
      dbdr(j,202) = db2dr(j,1208)
      dbdr(j,203) = db2dr(j,1211)
      dbdr(j,204) = db2dr(j,1214)
      dbdr(j,205) = db2dr(j,1225)
      dbdr(j,206) = db2dr(j,1228)
      dbdr(j,207) = db2dr(j,1231)
      dbdr(j,208) = db2dr(j,1240)
      dbdr(j,209) = db2dr(j,1243)
      dbdr(j,210) = db2dr(j,1249)
      dbdr(j,211) = db2dr(j,1328)
      dbdr(j,212) = db2dr(j,1329)
      dbdr(j,213) = db2dr(j,1330)
      dbdr(j,214) = db2dr(j,1331)
      dbdr(j,215) = db2dr(j,1332)
      dbdr(j,216) = db2dr(j,1334)
      dbdr(j,217) = db2dr(j,1335)
      dbdr(j,218) = db2dr(j,1337)
      dbdr(j,219) = db2dr(j,1338)
      dbdr(j,220) = db2dr(j,1340)
      dbdr(j,221) = db2dr(j,1341)
      dbdr(j,222) = db2dr(j,1343)
      dbdr(j,223) = db2dr(j,1936)
      dbdr(j,224) = db2dr(j,1941)
      dbdr(j,225) = db2dr(j,1942)
      dbdr(j,226) = db2dr(j,1949)
      dbdr(j,227) = db2dr(j,1950)
      dbdr(j,228) = db2dr(j,1951)
      dbdr(j,229) = db2dr(j,1961)
      dbdr(j,230) = db2dr(j,1962)
      dbdr(j,231) = db2dr(j,1965)
      dbdr(j,232) = db2dr(j,1968)
      dbdr(j,233) = db2dr(j,1982)
      dbdr(j,234) = db2dr(j,1985)
      dbdr(j,235) = db2dr(j,1988)
      dbdr(j,236) = db2dr(j,1991)
      dbdr(j,237) = db2dr(j,2003)
      dbdr(j,238) = db2dr(j,2006)
      dbdr(j,239) = db2dr(j,2009)
      dbdr(j,240) = db2dr(j,2018)
      dbdr(j,241) = db2dr(j,2021)
      dbdr(j,242) = db2dr(j,2027)
      dbdr(j,243) = db2dr(j,2126)
      dbdr(j,244) = db2dr(j,2127)
      dbdr(j,245) = db2dr(j,2128)
      dbdr(j,246) = db2dr(j,2129)
      dbdr(j,247) = db2dr(j,2130)
      dbdr(j,248) = db2dr(j,2132)
      dbdr(j,249) = db2dr(j,2133)
      dbdr(j,250) = db2dr(j,2135)
      dbdr(j,251) = db2dr(j,2136)
      dbdr(j,252) = db2dr(j,2138)
      dbdr(j,253) = db2dr(j,2139)
      dbdr(j,254) = db2dr(j,2141)
      dbdr(j,255) = db2dr(j,2142)
      dbdr(j,256) = db2dr(j,2144)

      enddo
      end subroutine evdbdr2

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
      else if (imol.eq.2) then
! iz for N2 system
      iz(1)=7
      iz(2)=7
! C6 for N2 system
      c6(1)=19.7d0
      c6(2)=19.7d0
      else
! currently imol = 4,3,2 are used
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

