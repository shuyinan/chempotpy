!**********************************************************************
!   System:                     N4
!   Name of this surface:       N4_1A_MB-PIP-MEG3
!   Functional form:            permutation-invariant polynomials
!   Common name:                N4(adiabatic ground state)
!   Number of derivatives:      1
!   Number of bodies:           4
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References:: J. Li, Z. Varga, D. G. Truhlar, and H. Guo
!                J. Chem. Theory Comput. 16, 4822 (2020).
!                Z. Varga, J. Li, D. G. Truhlar, and H. Guo, 
!                MEG3 surface for N4,
!                in: POTLIB: An Online Library of Potential Energy Surfaces,
!                comp.chem.umn.edu/potlib/ 
!
!   Notes:    -PES of N4 with special emphasize for N2 + N2 --> N2 + N + N             
!             -Instead of Morse variables mixed-exponential- gaussian (MEG)
!              variable is applied to describe the  long-range interactions
!              in a better way  
!             -The N4 data set was exteded to 21,406 points to make this fit
!             -This is a modified version of a potential originally published in 2020. 
!              The difference between this surface and the published surface is that 
!              the pairwise potentials of N2 are replaced with new ones that has a 
!              slightly smaller dissociation energy. After this replacement, the least
!              squares fits of the four-body potential were redetermined. With these
!              changes to the pairwise potentials, we have a unified set of potentials
!              for systems containing N2, O2, and NO that all use the same diatomic 
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
!
!     N1--N2
!
!     N3--N4
!
!
!
!   Input: X(4),Y(4),Z(4)               in units of bohr
!   Output: E                           in units of hartree
!   Output: dEdX(4),dEdY(4),dEdZ(4)     hartree/bohr
!**********************************************************************

      module N4_1A_ZV_par
!**********************************************************************
! This code is based on n4pes-gpip-meg.f (FORTRAN77 version) where some  
! keywords were repleced by more modern ones.
! The surface itself was not changed.
!**********************************************************************

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
! rMs(6):		Array to store base terms
! rM(0:111): 	Array to store monomials
! P(0:305): 	Array to store polynomials
! B(1:276):     Array to store basis functions
! dMsdR(6,6):   The derivative of base terms w.r.t. R
! dMdR(6,112):  The derivative of monomials w.r.t. R
! dPdR(6,306):  The derivative of basis functions w.r.t. R
! dVdR(6):      The derivative of V w.r.t. R
! dRdX(6,12):   The derivative of R w.r.t. X
! dBdR(6,276)	The derivative of B w.r.t. R 

      double precision :: R(6)
      double precision :: rMs(6),rM(0:111),P(0:305),B(276)
      double precision :: dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305)
      double precision :: dVdR(6),dRdX(6,12),dBdR(6,276)

! Nonlinear parameters:
! a(in Ang)
! ab (in Ang^2)
! re (in Ang)

      double precision,parameter :: a  = 1.0d0
      double precision,parameter :: ab = 1.5d0
      double precision,parameter :: re = 1.098d0

! Reference energy of infinitely separated N2 + N2 in hartree
      double precision,parameter :: Eref = &
        -218.40801323d0

! For N2 + N2 framework total diss. energy is 2De 2*228.4 kcal/mol 
      double precision,parameter :: totdiss = 456.8d0 

! Linear parameters optimized by the weighted-least square fitting
      double precision,parameter :: C(276)=(/ &
            0.188192931784D+03 &
          ,-0.306237562152D+03 &
          , 0.130602784443D+04 &
          ,-0.381429493182D+04 &
          , 0.366480199317D+03 &
          ,-0.115822947537D+05 &
          , 0.244890548213D+04 &
          ,-0.300501870278D+04 &
          , 0.349528727923D+04 &
          ,-0.226812742458D+04 &
          , 0.672933155942D+04 &
          , 0.738663493191D+04 &
          ,-0.553741947184D+04 &
          ,-0.116174474121D+05 &
          , 0.236514245770D+05 &
          ,-0.115827973083D+05 &
          ,-0.890607430768D+04 &
          ,-0.227070793998D+04 &
          , 0.584075896525D+04 &
          ,-0.923762797949D+04 &
          , 0.343366252227D+05 &
          ,-0.610081280668D+04 &
          , 0.155216809468D+05 &
          , 0.182955608656D+04 &
          ,-0.270762627366D+05 &
          ,-0.757005817323D+04 &
          ,-0.136267424389D+05 &
          , 0.171525707803D+05 &
          , 0.260001927913D+05 &
          , 0.488983563764D+04 &
          , 0.227593300424D+05 &
          ,-0.214507302673D+05 &
          , 0.132403556339D+05 &
          ,-0.162121812885D+05 &
          ,-0.234906947078D+05 &
          , 0.414001698466D+04 &
          , 0.934536680201D+04 &
          ,-0.321875656495D+05 &
          , 0.537108656402D+05 &
          , 0.125569495382D+05 &
          , 0.125827565028D+05 &
          , 0.617812384268D+04 &
          , 0.358105825296D+05 &
          , 0.689133913438D+04 &
          ,-0.112232440747D+05 &
          , 0.332322419445D+04 &
          , 0.143717862953D+05 &
          ,-0.514017008889D+05 &
          ,-0.167553499259D+04 &
          ,-0.298220045120D+05 &
          ,-0.143141726979D+05 &
          , 0.661871758185D+05 &
          , 0.932425833315D+04 &
          , 0.399181782780D+05 &
          , 0.653692283061D+04 &
          ,-0.239847242048D+05 &
          ,-0.242736693889D+05 &
          , 0.368668622270D+05 &
          , 0.149135204508D+05 &
          ,-0.120093927984D+05 &
          ,-0.213610885096D+05 &
          ,-0.982881587530D+04 &
          ,-0.123308425270D+05 &
          , 0.974725903137D+03 &
          , 0.180542233317D+05 &
          ,-0.243734199176D+05 &
          , 0.130789283408D+05 &
          , 0.145398844898D+05 &
          ,-0.309279778236D+05 &
          ,-0.298658246368D+05 &
          , 0.118355970237D+05 &
          , 0.259383063820D+05 &
          , 0.167714868005D+05 &
          ,-0.402760450928D+04 &
          ,-0.965094223701D+04 &
          ,-0.296261890206D+04 &
          , 0.204539028952D+05 &
          ,-0.684273581720D+05 &
          , 0.504938428655D+04 &
          ,-0.143733240001D+05 &
          , 0.912688917187D+04 &
          , 0.702365765217D+04 &
          ,-0.254638104362D+05 &
          ,-0.121043177737D+05 &
          ,-0.365263529440D+04 &
          ,-0.199658466702D+05 &
          , 0.597785266248D+05 &
          , 0.167249919316D+04 &
          , 0.104177096820D+05 &
          , 0.412717831787D+04 &
          , 0.458866378532D+04 &
          , 0.237236109775D+05 &
          , 0.496616378553D+04 &
          , 0.290083041659D+05 &
          , 0.105300678893D+05 &
          ,-0.544052249590D+05 &
          ,-0.682620876459D+04 &
          ,-0.297676028478D+05 &
          , 0.437948636203D+03 &
          , 0.168176871755D+05 &
          , 0.264412544651D+06 &
          ,-0.703298221236D+05 &
          , 0.437308434252D+05 &
          ,-0.210459792761D+05 &
          , 0.299011474876D+05 &
          ,-0.367449649294D+04 &
          , 0.207478298677D+05 &
          ,-0.450437768899D+05 &
          , 0.213163833075D+04 &
          , 0.153180223085D+05 &
          , 0.125152194411D+05 &
          , 0.987763003709D+04 &
          , 0.129531375662D+05 &
          ,-0.994456197027D+04 &
          ,-0.315107896702D+04 &
          , 0.112270895535D+05 &
          , 0.105227153267D+05 &
          , 0.185363730619D+05 &
          , 0.935339038161D+04 &
          ,-0.681313218429D+04 &
          ,-0.848428099047D+04 &
          , 0.104494779241D+05 &
          ,-0.244583475433D+04 &
          ,-0.804144763527D+04 &
          , 0.672613140154D+04 &
          , 0.234831157828D+04 &
          ,-0.205217461468D+05 &
          , 0.378050024658D+04 &
          ,-0.202752985250D+05 &
          , 0.487120262882D+03 &
          , 0.488799070169D+04 &
          ,-0.118843029077D+05 &
          , 0.210109485120D+05 &
          , 0.212394070709D+05 &
          ,-0.134081259330D+05 &
          ,-0.727298857900D+04 &
          ,-0.620346197265D+04 &
          ,-0.194921767050D+04 &
          , 0.895025904845D+04 &
          , 0.477160999708D+04 &
          ,-0.405499296058D+03 &
          , 0.339158632088D+05 &
          ,-0.143319790589D+05 &
          , 0.732934722586D+04 &
          ,-0.644007876470D+04 &
          , 0.754309910372D+04 &
          , 0.165120801084D+05 &
          , 0.566176764522D+04 &
          ,-0.430880244492D+04 &
          ,-0.721294942700D+04 &
          , 0.122802412901D+05 &
          ,-0.163773374737D+05 &
          , 0.315794156980D+04 &
          , 0.991397882845D+04 &
          , 0.840696528434D+04 &
          ,-0.265511935399D+05 &
          ,-0.478281541053D+04 &
          ,-0.449165333270D+04 &
          ,-0.362395225677D+04 &
          ,-0.108208967519D+05 &
          , 0.392595872067D+04 &
          ,-0.571760479408D+03 &
          ,-0.142046228936D+05 &
          , 0.360929747248D+03 &
          , 0.147850894009D+05 &
          , 0.345708299732D+03 &
          , 0.139692253303D+05 &
          , 0.113781489127D+05 &
          ,-0.107722426704D+04 &
          ,-0.545020740489D+04 &
          , 0.246535214545D+05 &
          ,-0.379050930238D+05 &
          ,-0.588836773928D+05 &
          ,-0.863518153891D+05 &
          , 0.254345105592D+05 &
          ,-0.117631679105D+05 &
          , 0.561354740383D+04 &
          , 0.192192418038D+04 &
          , 0.209853420114D+05 &
          , 0.917078313937D+04 &
          , 0.124138676580D+04 &
          ,-0.576787184887D+04 &
          ,-0.961772675450D+04 &
          , 0.292815384727D+04 &
          ,-0.154997225182D+04 &
          , 0.982738772033D+04 &
          , 0.804178317921D+04 &
          ,-0.204517484586D+05 &
          ,-0.986947881930D+02 &
          ,-0.730473836049D+04 &
          ,-0.451070530070D+04 &
          , 0.108242251735D+05 &
          ,-0.158087291728D+04 &
          ,-0.451726511855D+04 &
          ,-0.465630672307D+04 &
          ,-0.325586250630D+04 &
          ,-0.471732267072D+04 &
          , 0.711509671306D+04 &
          ,-0.483140408246D+04 &
          ,-0.898292555492D+03 &
          , 0.356353376348D+04 &
          ,-0.128133505742D+04 &
          , 0.236521591534D+04 &
          ,-0.606941261756D+04 &
          ,-0.934365415793D+03 &
          , 0.403623542444D+04 &
          ,-0.150415373610D+05 &
          , 0.218593690649D+04 &
          , 0.239425912292D+05 &
          , 0.118751815439D+04 &
          , 0.707746248672D+04 &
          ,-0.400959849220D+04 &
          ,-0.471221871177D+04 &
          ,-0.247453655730D+04 &
          ,-0.272060731517D+03 &
          , 0.359524281907D+04 &
          ,-0.767471874167D+03 &
          , 0.104562533220D+05 &
          ,-0.113530544437D+04 &
          ,-0.627005552570D+03 &
          ,-0.528294776271D+04 &
          ,-0.365881325444D+04 &
          , 0.802753001935D+04 &
          ,-0.755810620415D+03 &
          , 0.564906105226D+04 &
          , 0.841082956480D+04 &
          , 0.229038895725D+04 &
          ,-0.212225711775D+04 &
          ,-0.703508176770D+04 &
          , 0.355220351539D+04 &
          ,-0.393611409989D+04 &
          ,-0.571895247651D+04 &
          , 0.418207698558D+04 &
          ,-0.215998142449D+04 &
          , 0.108434080968D+04 &
          , 0.207284825950D+04 &
          ,-0.848695782933D+03 &
          , 0.144345809753D+03 &
          ,-0.255910695059D+04 &
          ,-0.521423247976D+04 &
          , 0.548206863634D+04 &
          ,-0.202149024490D+04 &
          , 0.403417606532D+03 &
          , 0.367549694781D+04 &
          ,-0.280190411808D+04 &
          , 0.656052307886D+03 &
          , 0.136449764528D+04 &
          ,-0.345055721086D+03 &
          ,-0.247697661346D+04 &
          ,-0.350809399685D+04 &
          ,-0.979188799848D+03 &
          ,-0.206935671572D+04 &
          , 0.472254638903D+04 &
          , 0.257532596476D+03 &
          , 0.864582177363D+03 &
          , 0.185994286891D+04 &
          ,-0.187970402618D+04 &
          , 0.590223360521D+04 &
          ,-0.315680401716D+04 &
          ,-0.242196459504D+04 &
          ,-0.192118722423D+04 &
          , 0.568224413848D+04 &
          , 0.341351291630D+04 &
          , 0.783527450538D+03 &
          , 0.990336060360D+03 &
          , 0.386766649184D+04 &
          ,-0.383330346343D+04 &
          ,-0.107913792772D+04 &
          , 0.273834922602D+04 &
          ,-0.165595102315D+04 &
          , 0.148212652224D+03 &
          , 0.108781989583D+04 &
          ,-0.248989912421D+04 &
          ,-0.205426012347D+04 &
          , 0.222220909309D+03 &
          , 0.553262644100D+03 &
      /)

      end module N4_1A_ZV_par

      subroutine pes(x,igrad,potential,gradient,dvec)

      use N4_1A_ZV_par
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=4
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
      do iatom=1,natoms
        do idir=1,3
          j=3*(iatom-1)+idir
          tx(j)=x(iatom, idir)
        enddo
      enddo
      tg=0.d0

      call n4pes(tx, v, tg, igrad)

      ! output v is in kcal/mol, but there is a reference energy that is
      ! in hartree, so convert everything in eV.
      v=v/23.0609 ! + Eref*27.211386d0

      ! output g is in kcal/mol/ang
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

      subroutine n4pes(X,v,dVdX,igrad)
      use N4_1A_ZV_par
!**********************************************************************
! Subroutine to calculate the potential energy V and gradient dVdX
! for given Cartesian coordinates X(12)  
! R:		Interatomic bond distance (6)
! V:		Calculated potential energy
! dVdX:		The derivative of V w.r.t. X, dim(12)
! dVdR:		The derivative of V w.r.t. R, dim(6) 
! dPdR:		The derivative of basis functions w.r.t. R
!		dim(6*306)
! dMdR:		The derivative of monomials w.r.t. R
!		dim(6*112)
! dRdX:		The derivative of R w.r.t. X, dim(6*12)
!**********************************************************************

      integer i,igrad,j,nob,k
      double precision V
      double precision dVdX(12),X(12) 

! Read Cartesian coordinate from input file
      call coord_convt(X)

      if (igrad .le. 1) then
! Call subroutine Evv to evaluate potential energy V
        call evv(V)

        if (igrad .eq. 1) then
! Call EvdVdX to evaluate the derivatives of V w.r.t. X
          call evdvdx(X,dVdX)
        endif
      else
        write (*,*) 'Only igrad = 0, 1 is allowed!'
      endif

      end subroutine n4pes

      subroutine coord_convt(X)
      use N4_1A_ZV_par
!**********************************************************************
!  Program to calculate the six interatomic distance 
!  by reading XYZ coordinate
!**********************************************************************

      integer i
      double precision X(12)
      
!**********************************************************************
!  Now, calculate the inter-atomic distance
!  r1 = r(N1N2)    r2 = r(N1N3)
!  r3 = r(N1N4)    r4 = r(N2N3)
!  r5 = r(N2N4)    r6 = r(N3N4)       
!**********************************************************************

      R(1)=Sqrt((X(4)-X(1))**2 + (X(5)-X(2))**2 + (X(6)-X(3))**2)
      R(2)=Sqrt((X(7)-X(1))**2 + (X(8)-X(2))**2 + (X(9)-X(3))**2)
      R(3)=Sqrt((X(10)-X(1))**2 + (X(11)-X(2))**2 + (X(12)-X(3))**2)
      R(4)=Sqrt((X(4)-X(7))**2 + (X(5)-X(8))**2 + (X(6)-X(9))**2)
      R(5)=Sqrt((X(4)-X(10))**2 + (X(5)-X(11))**2 + (X(6)-X(12))**2)
      R(6)=Sqrt((X(7)-X(10))**2 + (X(8)-X(11))**2 + (X(9)-X(12))**2)
 
      return

      end subroutine coord_convt

      subroutine EvV(V)
      use N4_1A_ZV_par
!**********************************************************************
! Subroutine to evaluate V for given R 
! V(R) = C*P
! C:		Coefficients, stored in 'dim.inc' 
! P:		Basis functions evaluated for given R
! rMs:		rMs(6), six mixed exponential gaussian terms (MEG)
! a:		Nonlinear parameters in Morse terms(Angstrom)
! ab:		Nonlinear parameters in Gauss terms(Angstrom^2)
! re:		Equilibrium bond length(Angstrom)
! nop:		number of points
! nom:  	number of monomials
! nob:  	number of basis functions(polynomials)
! rM(0:111):	Array to store monomials
! P(0:305):	Array to store polynomials
! B(1:276):     Array to store basis functions
!**********************************************************************

      integer i,j,k
      double precision dist,dv2dr,V,V2

! Calculate the six MEG terms for each point
      call evmorse

! Calculate the monomials for each point by using six MEG terms
      call evmono

! Calculate the polynomials (basis functions) by using monomials
      call evpoly 

! Calculate the basis functions by removing unconnected and 2-body terms
      call evbas

! Initialized v to be totdiss
      v = totdiss
! Evaluate 2-body interactions
      do i=1,6
        dist=r(i)
        call ev2gm2(dist,v2,dv2dr,1,0)
        v=v+v2
      enddo

! Evaluate V by taken the product of C and Basis function array
      do i=1,276
        v=v + c(i)*b(i)
      enddo

!      Write(*,9999) V 
! 9999 Format('The potential energy is ',F20.14,' kcal/mol')

      return

      end subroutine EvV

      subroutine EvdVdX(X,dVdX)
      use N4_1A_ZV_par
!**********************************************************************
! Subroutine to evaluate dRdX for given R and X 
! R:		R(6), 6 bond lengths
! X:		X(12), 12 Cartesian coordinates
! rM(0:111):    Array to store monomials
! P(0:305):     Array to store polynomials
! dVdX:		dVdX(12), derivatives of V w.r.t. Cartesian coordinates 
! dVdR:		dVdR(6), derivatives of V w.r.t.6 bond length
! dRdX:		dRdX(6,12), derivatives of R(6) w.r.t. 12  
!		Cartesian coordinates
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

! Calculate dVdX by using chain rule: dV/dXi=(dV/dRj)*(dRj/dXi), j=1 to 6
      do i=1,12
        do j=1,6
          dVdX(i)=dVdX(i) + dVdR(j)*dRdX(j,i)
        enddo
      enddo

!      write(*,*) 'The 12 dVdX are:'
!      Write(*,9999) (dVdX(i),i=1,12) 
! 9999 Format(1x,3F15.8)

      return
      end subroutine EvdVdX

      subroutine EvMorse
      use N4_1A_ZV_par
!**********************************************************************
! mixed exponential gaussian term rms = exp(-(r-re)/a-(r-re)^2/ab)
! re:	equlibrium bond length
! a: 	nonlinear parameter, unit Ang
! ab:	nonlinear parameter, unit Ang^2   
!**********************************************************************

      integer i
      
      do i=1,6
         rms(i)=Exp(-(r(i)-re)/a-((r(i)-re)**2.0d0)/ab)
      enddo

      end subroutine EvMorse

      subroutine EvMono
      use N4_1A_ZV_par
!**********************************************************************
!  The subroutine reads six MEG variables(X) and calculates the
!  monomials(M) that do not have usable decomposition.
!  For A4 with max. degree 9, the number of monomials is nom.
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
      rm(9) = rm(1)*rm(6)
      rm(10) = rm(1)*rm(2)
      rm(11) = rm(1)*rm(3)
      rm(12) = rm(2)*rm(3)
      rm(13) = rm(1)*rm(4)
      rm(14) = rm(2)*rm(4)
      rm(15) = rm(1)*rm(5)
      rm(16) = rm(3)*rm(5)
      rm(17) = rm(4)*rm(5)
      rm(18) = rm(2)*rm(6)
      rm(19) = rm(3)*rm(6)
      rm(20) = rm(4)*rm(6)
      rm(21) = rm(5)*rm(6)
      rm(22) = rm(1)*rm(7)
      rm(23) = rm(2)*rm(7)
      rm(24) = rm(1)*rm(8)
      rm(25) = rm(2)*rm(16)
      rm(26) = rm(2)*rm(17)
      rm(27) = rm(3)*rm(17)
      rm(28) = rm(1)*rm(18)
      rm(29) = rm(1)*rm(19)
      rm(30) = rm(1)*rm(20)
      rm(31) = rm(3)*rm(20)
      rm(32) = rm(1)*rm(21)
      rm(33) = rm(2)*rm(21)
      rm(34) = rm(1)*rm(12)
      rm(35) = rm(1)*rm(17)
      rm(36) = rm(2)*rm(20)
      rm(37) = rm(3)*rm(21)
      rm(38) = rm(1)*rm(14)
      rm(39) = rm(1)*rm(16)
      rm(40) = rm(2)*rm(19)
      rm(41) = rm(4)*rm(21)
      rm(42) = rm(2)*rm(27)
      rm(43) = rm(1)*rm(31)
      rm(44) = rm(1)*rm(33)
      rm(45) = rm(1)*rm(23)
      rm(46) = rm(1)*rm(25)
      rm(47) = rm(1)*rm(26)
      rm(48) = rm(1)*rm(27)
      rm(49) = rm(1)*rm(40)
      rm(50) = rm(1)*rm(36)
      rm(51) = rm(2)*rm(31)
      rm(52) = rm(1)*rm(37)
      rm(53) = rm(2)*rm(37)
      rm(54) = rm(1)*rm(41)
      rm(55) = rm(2)*rm(41)
      rm(56) = rm(3)*rm(41)
      rm(57) = rm(1)*rm(42)
      rm(58) = rm(1)*rm(51)
      rm(59) = rm(1)*rm(53)
      rm(60) = rm(1)*rm(55)
      rm(61) = rm(1)*rm(56)
      rm(62) = rm(2)*rm(56)
      rm(63) = rm(1)*rm(62)
      rm(64) = rm(2)*rm(57)
      rm(65) = rm(3)*rm(57)
      rm(66) = rm(4)*rm(57)
      rm(67) = rm(5)*rm(57)
      rm(68) = rm(1)*rm(58)
      rm(69) = rm(3)*rm(58)
      rm(70) = rm(4)*rm(58)
      rm(71) = rm(1)*rm(59)
      rm(72) = rm(2)*rm(59)
      rm(73) = rm(1)*rm(60)
      rm(74) = rm(2)*rm(60)
      rm(75) = rm(1)*rm(61)
      rm(76) = rm(2)*rm(62)
      rm(77) = rm(3)*rm(61)
      rm(78) = rm(3)*rm(62)
      rm(79) = rm(4)*rm(61)
      rm(80) = rm(4)*rm(62)
      rm(81) = rm(5)*rm(59)
      rm(82) = rm(5)*rm(60)
      rm(83) = rm(5)*rm(62)
      rm(84) = rm(6)*rm(58)
      rm(85) = rm(6)*rm(59)
      rm(86) = rm(6)*rm(60)
      rm(87) = rm(6)*rm(61)
      rm(88) = rm(2)*rm(64)
      rm(89) = rm(3)*rm(65)
      rm(90) = rm(4)*rm(66)
      rm(91) = rm(5)*rm(67)
      rm(92) = rm(1)*rm(68)
      rm(93) = rm(3)*rm(69)
      rm(94) = rm(4)*rm(70)
      rm(95) = rm(1)*rm(71)
      rm(96) = rm(2)*rm(72)
      rm(97) = rm(1)*rm(73)
      rm(98) = rm(2)*rm(74)
      rm(99) = rm(1)*rm(75)
      rm(100) = rm(2)*rm(76)
      rm(101) = rm(3)*rm(77)
      rm(102) = rm(3)*rm(78)
      rm(103) = rm(4)*rm(79)
      rm(104) = rm(4)*rm(80)
      rm(105) = rm(5)*rm(81)
      rm(106) = rm(5)*rm(82)
      rm(107) = rm(5)*rm(83)
      rm(108) = rm(6)*rm(84)
      rm(109) = rm(6)*rm(85)
      rm(110) = rm(6)*rm(86)
      rm(111) = rm(6)*rm(87)

      return

      end subroutine EvMono

      subroutine EvPoly
      use N4_1A_ZV_par
!**********************************************************************
!  The subroutine reads monomials(m) and calculates the
!  permutation-invariant polynomials(p)
!  For A4 with max. degree 9, the number of polynomials is nob.
!**********************************************************************

      p(0) = rm(0)
      p(1) = rm(1) + rm(2) + rm(3) + rm(4) + rm(5) + rm(6)
      p(2) = rm(7) + rm(8) + rm(9)
      p(3) = rm(10) + rm(11) + rm(12) + rm(13) + rm(14) + rm(15) &
           + rm(16) + rm(17) + rm(18) + rm(19) + rm(20) + rm(21)
      p(4) = p(1)*p(1) - p(3) - p(2) - p(3) - p(2)
      p(5) = rm(22) + rm(23) + rm(24) + rm(25) + rm(26) + rm(27) &
           + rm(28) + rm(29) + rm(30) + rm(31) + rm(32) + rm(33)
      p(6) = rm(34) + rm(35) + rm(36) + rm(37)
      p(7) = rm(38) + rm(39) + rm(40) + rm(41)
      p(8) = p(1)*p(2) - p(5)
      p(9) = p(1)*p(3) - p(6) - p(7) - p(5) - p(6) &
            - p(7) - p(5) - p(6) - p(7)
      p(10) = p(1)*p(4) - p(9) - p(8)
      p(11) = rm(42) + rm(43) + rm(44)
      p(12) = rm(45) + rm(46) + rm(47) + rm(48) + rm(49) + rm(50) &
            + rm(51) + rm(52) + rm(53) + rm(54) + rm(55) + rm(56)
      p(13) = p(2)*p(3) - p(12)
      p(14) = p(1)*p(5) - p(12) - p(11) - p(13) - p(12) &
            - p(11) - p(11) - p(11)
      p(15) = p(1)*p(6) - p(12)
      p(16) = p(1)*p(7) - p(12)
      p(17) = p(2)*p(2) - p(11) - p(11)
      p(18) = p(3)*p(3) - p(12) - p(11) - p(15) - p(16) &
            - p(14) - p(12) - p(11) - p(15) - p(16) - p(14) &
            - p(12) - p(11) - p(12) - p(11)
      p(19) = p(2)*p(4) - p(14)
      p(20) = p(3)*p(4) - p(15) - p(16) - p(13)
      p(21) = p(1)*p(10) - p(20) - p(19)
      p(22) = rm(57) + rm(58) + rm(59) + rm(60) + rm(61) + rm(62)
      p(23) = p(1)*p(11) - p(22)
      p(24) = p(2)*p(6)
      p(25) = p(2)*p(7)
      p(26) = p(1)*p(12) - p(22) - p(24) - p(25) - p(22) &
            - p(22) - p(22)
      p(27) = p(2)*p(5) - p(22) - p(23) - p(22)
      p(28) = p(3)*p(5) - p(22) - p(26) - p(24) - p(25) &
            - p(23) - p(22) - p(24) - p(25) - p(23) - p(22) &
            - p(22)
      p(29) = p(3)*p(6) - p(22) - p(26) - p(22)
      p(30) = p(3)*p(7) - p(22) - p(26) - p(22)
      p(31) = p(2)*p(9) - p(26) - p(28)
      p(32) = p(1)*p(14) - p(26) - p(23) - p(28)
      p(33) = p(4)*p(6) - p(25)
      p(34) = p(4)*p(7) - p(24)
      p(35) = p(1)*p(17) - p(27)
      p(36) = p(1)*p(18) - p(29) - p(30) - p(28)
      p(37) = p(2)*p(10) - p(32)
      p(38) = p(3)*p(10) - p(33) - p(34) - p(31)
      p(39) = p(1)*p(21) - p(38) - p(37)
      p(40) = rm(63)
      p(41) = rm(64) + rm(65) + rm(66) + rm(67) + rm(68) + rm(69) &
            + rm(70) + rm(71) + rm(72) + rm(73) + rm(74) + rm(75) &
            + rm(76) + rm(77) + rm(78) + rm(79) + rm(80) + rm(81) &
            + rm(82) + rm(83) + rm(84) + rm(85) + rm(86) + rm(87)
      p(42) = p(1)*p(22) - p(40) - p(41) - p(40) - p(40) &
            - p(40) - p(40) - p(40)
      p(43) = p(2)*p(11) - p(40) - p(40) - p(40)
      p(44) = p(2)*p(12) - p(41)
      p(45) = p(3)*p(11) - p(41)
      p(46) = p(5)*p(6) - p(41)
      p(47) = p(5)*p(7) - p(41)
      p(48) = p(6)*p(7) - p(40) - p(40) - p(40) - p(40)
      p(49) = p(4)*p(11) - p(42)
      p(50) = p(2)*p(15) - p(46)
      p(51) = p(2)*p(16) - p(47)
      p(52) = p(4)*p(12) - p(41) - p(50) - p(51)
      p(53) = p(2)*p(14) - p(42) - p(49) - p(42)
      p(54) = p(6)*p(6) - p(42) - p(42)
      p(55) = p(7)*p(7) - p(42) - p(42)
      p(56) = p(3)*p(17) - p(44)
      p(57) = p(2)*p(18) - p(48)
      p(58) = p(3)*p(14) - p(41) - p(52) - p(46) - p(47) &
            - p(45) - p(45)
      p(59) = p(6)*p(9) - p(41) - p(52) - p(47)
      p(60) = p(7)*p(9) - p(41) - p(52) - p(46)
      p(61) = p(2)*p(20) - p(52) - p(58)
      p(62) = p(1)*p(32) - p(52) - p(49) - p(58)
      p(63) = p(6)*p(10) - p(51)
      p(64) = p(7)*p(10) - p(50)
      p(65) = p(2)*p(17) - p(43)
      p(66) = p(3)*p(18) - p(46) - p(47) - p(45) - p(59) &
            - p(60) - p(58)
      p(67) = p(2)*p(19) - p(49)
      p(68) = p(1)*p(36) - p(59) - p(60) - p(58) - p(57) &
            - p(66) - p(66)
      p(69) = p(2)*p(21) - p(62)
      p(70) = p(3)*p(21) - p(63) - p(64) - p(61)
      p(71) = p(1)*p(39) - p(70) - p(69)
      p(72) = p(40)*p(1)
      p(73) = p(2)*p(22) - p(72)
      p(74) = p(6)*p(11)
      p(75) = p(7)*p(11)
      p(76) = p(3)*p(22) - p(72) - p(74) - p(75) - p(72) &
            - p(72) - p(72)
      p(77) = rm(88) + rm(89) + rm(90) + rm(91) + rm(92) + rm(93) &
            + rm(94) + rm(95) + rm(96) + rm(97) + rm(98) + rm(99) &
            + rm(100) + rm(101) + rm(102) + rm(103) + rm(104) &
            + rm(105) + rm(106) + rm(107) + rm(108) + rm(109) &
            + rm(110) + rm(111)
      p(78) = p(1)*p(42) - p(72) - p(76)
      p(79) = p(5)*p(11) - p(72) - p(73) - p(72)
      p(80) = p(2)*p(26) - p(76) - p(77)
      p(81) = p(6)*p(12) - p(72) - p(76) - p(72)
      p(82) = p(7)*p(12) - p(72) - p(76) - p(72)
      p(83) = p(8)*p(11) - p(72)
      p(84) = p(6)*p(17)
      p(85) = p(7)*p(17)
      p(86) = p(9)*p(11) - p(76) - p(77)
      p(87) = p(2)*p(29) - p(81)
      p(88) = p(6)*p(14) - p(75) - p(75)
      p(89) = p(2)*p(30) - p(82)
      p(90) = p(7)*p(14) - p(74) - p(74)
      p(91) = p(1)*p(48) - p(76) - p(81) - p(82)
      p(92) = p(10)*p(11) - p(78)
      p(93) = p(2)*p(33) - p(88)
      p(94) = p(2)*p(34) - p(90)
      p(95) = p(10)*p(12) - p(77) - p(93) - p(94)
      p(96) = p(2)*p(27) - p(73) - p(79)
      p(97) = p(2)*p(28) - p(76) - p(86)
      p(98) = p(1)*p(53) - p(80) - p(79) - p(97)
      p(99) = p(1)*p(54) - p(81)
      p(100) = p(1)*p(55) - p(82)
      p(101) = p(5)*p(18) - p(76) - p(91) - p(87) - p(89) &
             - p(86)
      p(102) = p(6)*p(18) - p(74) - p(90)
      p(103) = p(7)*p(18) - p(75) - p(88)
      p(104) = p(2)*p(31) - p(77) - p(86)
      p(105) = p(2)*p(36) - p(91) - p(101)
      p(106) = p(3)*p(32) - p(77) - p(95) - p(88) - p(90) &
             - p(86)
      p(107) = p(4)*p(29) - p(82) - p(80) - p(99)
      p(108) = p(4)*p(30) - p(81) - p(80) - p(100)
      p(109) = p(2)*p(38) - p(95) - p(106)
      p(110) = p(1)*p(62) - p(95) - p(92) - p(106)
      p(111) = p(6)*p(21) - p(94)
      p(112) = p(7)*p(21) - p(93)
      p(113) = p(1)*p(65) - p(96)
      p(114) = p(1)*p(66) - p(102) - p(103) - p(101)
      p(115) = p(2)*p(37) - p(92)
      p(116) = p(10)*p(18) - p(99) - p(100) - p(97)
      p(117) = p(2)*p(39) - p(110)
      p(118) = p(3)*p(39) - p(111) - p(112) - p(109)
      p(119) = p(1)*p(71) - p(118) - p(117)
      p(120) = p(40)*p(2)
      p(121) = p(40)*p(3)
      p(122) = p(40)*p(4)
      p(123) = p(11)*p(12) - p(121)
      p(124) = p(2)*p(42) - p(122)
      p(125) = p(6)*p(22) - p(121)
      p(126) = p(7)*p(22) - p(121)
      p(127) = p(2)*p(41) - p(121) - p(123) - p(121)
      p(128) = p(6)*p(23) - p(123)
      p(129) = p(7)*p(23) - p(123)
      p(130) = p(3)*p(41) - p(122) - p(121) - p(120) - p(125) &
             - p(128) - p(126) - p(129) - p(124) - p(123) - p(122) &
             - p(121) - p(120) - p(125) - p(126) - p(124) - p(123) &
             - p(122) - p(121) - p(120) - p(122) - p(121) - p(120) &
             - p(120) - p(120) - p(120) - p(120)
      p(131) = p(3)*p(42) - p(121) - p(125) - p(126) - p(121)
      p(132) = p(4)*p(41) - p(121) - p(131) - p(128) - p(129) &
             - p(127) - p(121)
      p(133) = p(1)*p(78) - p(122) - p(131)
      p(134) = p(11)*p(11) - p(120) - p(120)
      p(135) = p(2)*p(48) - p(130)
      p(136) = p(11)*p(17) - p(120)
      p(137) = p(2)*p(44) - p(123)
      p(138) = p(2)*p(45) - p(121)
      p(139) = p(11)*p(14) - p(122) - p(124) - p(122)
      p(140) = p(6)*p(27) - p(127)
      p(141) = p(2)*p(54)
      p(142) = p(7)*p(27) - p(127)
      p(143) = p(2)*p(52) - p(131) - p(132)
      p(144) = p(1)*p(81) - p(125) - p(141) - p(135) - p(125)
      p(145) = p(2)*p(55)
      p(146) = p(1)*p(82) - p(126) - p(135) - p(145) - p(126)
      p(147) = p(11)*p(18) - p(130)
      p(148) = p(6)*p(28) - p(129) - p(123) - p(143)
      p(149) = p(7)*p(28) - p(128) - p(123) - p(143)
      p(150) = p(6)*p(30) - p(121) - p(146)
      p(151) = p(11)*p(19) - p(122)
      p(152) = p(2)*p(50) - p(128)
      p(153) = p(2)*p(51) - p(129)
      p(154) = p(11)*p(20) - p(131) - p(132)
      p(155) = p(2)*p(59) - p(144) - p(148)
      p(156) = p(6)*p(32) - p(129)
      p(157) = p(2)*p(60) - p(146) - p(149)
      p(158) = p(7)*p(32) - p(128)
      p(159) = p(6)*p(34) - p(122) - p(145) - p(122)
      p(160) = p(11)*p(21) - p(133)
      p(161) = p(2)*p(63) - p(156)
      p(162) = p(2)*p(64) - p(158)
      p(163) = p(12)*p(21) - p(132) - p(161) - p(162)
      p(164) = p(2)*p(53) - p(124) - p(139)
      p(165) = p(2)*p(58) - p(131) - p(154)
      p(166) = p(3)*p(54) - p(125) - p(144)
      p(167) = p(3)*p(55) - p(126) - p(146)
      p(168) = p(3)*p(65) - p(137)
      p(169) = p(17)*p(18) - p(135)
      p(170) = p(1)*p(98) - p(143) - p(139) - p(165)
      p(171) = p(4)*p(54) - p(135)
      p(172) = p(4)*p(55) - p(135)
      p(173) = p(2)*p(66) - p(150)
      p(174) = p(1)*p(101) - p(148) - p(149) - p(147) - p(173) &
            - p(165) - p(147)
      p(175) = p(1)*p(102) - p(150) - p(148) - p(166)
      p(176) = p(1)*p(103) - p(150) - p(149) - p(167)
      p(177) = p(2)*p(61) - p(132) - p(154)
      p(178) = p(2)*p(68) - p(159) - p(174)
      p(179) = p(3)*p(62) - p(132) - p(163) - p(156) - p(158) &
             - p(154)
      p(180) = p(6)*p(38) - p(132) - p(163) - p(157)
      p(181) = p(7)*p(38) - p(132) - p(163) - p(155)
      p(182) = p(2)*p(70) - p(163) - p(179)
      p(183) = p(1)*p(110) - p(163) - p(160) - p(179)
      p(184) = p(6)*p(39) - p(162)
      p(185) = p(7)*p(39) - p(161)
      p(186) = p(2)*p(65) - p(136)
      p(187) = p(3)*p(66) - p(148) - p(149) - p(147) - p(175) &
             - p(176) - p(174)
      p(188) = p(2)*p(67) - p(151)
      p(189) = p(4)*p(66) - p(166) - p(167) - p(165)
      p(190) = p(2)*p(69) - p(160)
      p(191) = p(18)*p(21) - p(171) - p(172) - p(169)
      p(192) = p(2)*p(71) - p(183)
      p(193) = p(3)*p(71) - p(184) - p(185) - p(182)
      p(194) = p(1)*p(119) - p(193) - p(192)
      p(195) = p(40)*p(5)
      p(196) = p(40)*p(6)
      p(197) = p(40)*p(7)
      p(198) = p(40)*p(8)
      p(199) = p(40)*p(9)
      p(200) = p(40)*p(10)
      p(201) = p(11)*p(22) - p(195)
      p(202) = p(12)*p(22) - p(196) - p(197) - p(195) - p(196) &
             - p(197) - p(195) - p(196) - p(197)
      p(203) = p(17)*p(22) - p(198)
      p(204) = p(6)*p(43)
      p(205) = p(7)*p(43)
      p(206) = p(11)*p(26) - p(199) - p(202)
      p(207) = p(2)*p(76) - p(199) - p(202)
      p(208) = p(2)*p(78) - p(200)
      p(209) = p(6)*p(41) - p(199) - p(195) - p(202) - p(195)
      p(210) = p(6)*p(42) - p(197) - p(197) - p(197)
      p(211) = p(7)*p(41) - p(199) - p(195) - p(202) - p(195)
      p(212) = p(7)*p(42) - p(196) - p(196) - p(196)
      p(213) = p(11)*p(29) - p(209)
      p(214) = p(11)*p(30) - p(211)
      p(215) = p(18)*p(22) - p(199) - p(213) - p(214)
      p(216) = p(2)*p(77) - p(199) - p(206)
      p(217) = p(6)*p(49) - p(205)
      p(218) = p(7)*p(49) - p(204)
      p(219) = p(3)*p(77) - p(200) - p(199) - p(198) - p(209) &
             - p(217) - p(211) - p(218) - p(207) - p(204) - p(205) &
             - p(200) - p(199) - p(198) - p(200) - p(198) - p(200) &
             - p(198)
      p(220) = p(3)*p(78) - p(199) - p(210) - p(212)
      p(221) = p(10)*p(41) - p(199) - p(220) - p(217) - p(218) &
             - p(216)
      p(222) = p(1)*p(133) - p(200) - p(220)
      p(223) = p(11)*p(27) - p(195) - p(203)
      p(224) = p(1)*p(134) - p(201)
      p(225) = p(2)*p(81) - p(209)
      p(226) = p(2)*p(80) - p(202) - p(206)
      p(227) = p(2)*p(82) - p(211)
      p(228) = p(2)*p(91) - p(215) - p(219)
      p(229) = p(11)*p(28) - p(199) - p(207)
      p(230) = p(6)*p(53) - p(205)
      p(231) = p(5)*p(54) - p(209)
      p(232) = p(7)*p(53) - p(204)
      p(233) = p(7)*p(54) - p(196)
      p(234) = p(5)*p(55) - p(211)
      p(235) = p(6)*p(55) - p(197)
      p(236) = p(11)*p(35) - p(198)
      p(237) = p(6)*p(65)
      p(238) = p(7)*p(65)
      p(239) = p(2)*p(86) - p(199) - p(229)
      p(240) = p(1)*p(139) - p(206) - p(229) - p(224)
      p(241) = p(17)*p(29) - p(225)
      p(242) = p(2)*p(99) - p(231)
      p(243) = p(17)*p(30) - p(227)
      p(244) = p(2)*p(95) - p(220) - p(221)
      p(245) = p(4)*p(81) - p(202) - p(242) - p(227)
      p(246) = p(2)*p(100) - p(234)
      p(247) = p(4)*p(82) - p(202) - p(225) - p(246)
      p(248) = p(11)*p(36) - p(215) - p(219)
      p(249) = p(2)*p(102) - p(233)
      p(250) = p(6)*p(58) - p(206) - p(214) - p(244) - p(214)
      p(251) = p(2)*p(103) - p(235)
      p(252) = p(7)*p(58) - p(213) - p(206) - p(244) - p(213)
      p(253) = p(1)*p(150) - p(215) - p(233) - p(235)
      p(254) = p(11)*p(37) - p(200)
      p(255) = p(2)*p(93) - p(217)
      p(256) = p(2)*p(94) - p(218)
      p(257) = p(11)*p(38) - p(220) - p(221)
      p(258) = p(2)*p(107) - p(245) - p(250)
      p(259) = p(6)*p(62) - p(218)
      p(260) = p(2)*p(108) - p(247) - p(252)
      p(261) = p(7)*p(62) - p(217)
      p(262) = p(6)*p(64) - p(200) - p(246) - p(200)
      p(263) = p(11)*p(39) - p(222)
      p(264) = p(2)*p(111) - p(259)
      p(265) = p(2)*p(112) - p(261)
      p(266) = p(12)*p(39) - p(221) - p(264) - p(265)
      p(267) = p(2)*p(98) - p(208) - p(240)
      p(268) = p(6)*p(54) - p(210)
      p(269) = p(7)*p(55) - p(212)
      p(270) = p(2)*p(96) - p(203) - p(223)
      p(271) = p(2)*p(97) - p(207) - p(229)
      p(272) = p(2)*p(101) - p(215) - p(248)
      p(273) = p(2)*p(106) - p(220) - p(257)
      p(274) = p(6)*p(59) - p(220) - p(215) - p(209)
      p(275) = p(7)*p(60) - p(220) - p(215) - p(211)
      p(276) = p(5)*p(66) - p(215) - p(253) - p(249) - p(251) &
             - p(248)
      p(277) = p(6)*p(66) - p(213) - p(252)
      p(278) = p(7)*p(66) - p(214) - p(250)
      p(279) = p(2)*p(104) - p(216) - p(239)
      p(280) = p(2)*p(105) - p(219) - p(248)
      p(281) = p(1)*p(170) - p(244) - p(240) - p(273)
      p(282) = p(10)*p(54) - p(227)
      p(283) = p(10)*p(55) - p(225)
      p(284) = p(2)*p(114) - p(253) - p(276)
      p(285) = p(1)*p(174) - p(250) - p(252) - p(248) - p(276) &
             - p(273)
      p(286) = p(6)*p(68) - p(217) - p(261) - p(251)
      p(287) = p(7)*p(68) - p(218) - p(259) - p(249)
      p(288) = p(2)*p(109) - p(221) - p(257)
      p(289) = p(2)*p(116) - p(262) - p(285)
      p(290) = p(3)*p(110) - p(221) - p(266) - p(259) - p(261) &
             - p(257)
      p(291) = p(6)*p(70) - p(221) - p(266) - p(260)
      p(292) = p(7)*p(70) - p(221) - p(266) - p(258)
      p(293) = p(2)*p(118) - p(266) - p(290)
      p(294) = p(1)*p(183) - p(266) - p(263) - p(290)
      p(295) = p(6)*p(71) - p(265)
      p(296) = p(7)*p(71) - p(264)
      p(297) = p(1)*p(186) - p(270)
      p(298) = p(1)*p(187) - p(277) - p(278) - p(276)
      p(299) = p(2)*p(115) - p(254)
      p(300) = p(1)*p(189) - p(286) - p(287) - p(285) - p(284) &
             - p(298)
      p(301) = p(2)*p(117) - p(263)
      p(302) = p(18)*p(39) - p(282) - p(283) - p(280)
      p(303) = p(2)*p(119) - p(294)
      p(304) = p(3)*p(119) - p(295) - p(296) - p(293)
      p(305) = p(1)*p(194) - p(304) - p(303)

      return

      end subroutine EvPoly

      subroutine evbas
      use N4_1A_ZV_par
!**********************************************************************
!  The subroutine eliminates the 2-body terms in Bowman's approach
!**********************************************************************

      integer i
      double precision b1(306) 

!	  Pass P(0:305) to BM1(1:306)
      do i=1,306
        b1(i)=p(i-1)
      enddo


! Remove unconnected terms and 2-body terms and pass to B(1:276)
      b(1)=b1(4)

      do i=2,4
        b(i)=b1(i+4)
      enddo

      b(5)=b1(10)

      do i=6,11
        b(i)=b1(i+6)
      enddo

      b(12)=b1(19)
      b(13)=b1(21)

      do i=14,26
        b(i)=b1(i+9)
      enddo

      b(27)=b1(37)
      b(28)=b1(39)

      do i=29,53
        b(i)=b1(i+12)
      enddo

      b(54)=b1(67)
      b(55)=b1(69)
      b(56)=b1(71)

      do i=57,97
        b(i)=b1(i+16)
      enddo

      b(98)=b1(115)
      b(99)=b1(117)
      b(100)=b1(119)

      do i=101,166
        b(i)=b1(i+20)
      enddo

      b(167)=b1(188)
      b(168)=b1(190)
      b(169)=b1(192)
      b(170)=b1(194)

      do i=171,272
        b(i)=b1(i+25)
      enddo

      b(273)=b1(299)
      b(274)=b1(301)
      b(275)=b1(303)
      b(276)=b1(305)

      return

      end subroutine evbas

      subroutine ev2gm2(r,v,grad,imol,igrad)
!**********************************************************************
!
! Compute the diatomic potential of ground-state N2 (singlet)
!
! Input:  r      interatomic distance in Angstrom
!         imol   1 - stands for N2 for N4 PES
!         igrad  1 - gradient is calculated
!                0 - gradient is not calculated
! Output: V      potential in kcal/mol
!         grad   gradient (kcal/mol)/Angstrom
!
!**********************************************************************
      implicit none

      integer,intent(in) :: imol, igrad
      double precision :: r
      double precision :: v, grad
      double precision :: re, de, as(0:6)
      double precision :: y, fy, u
      double precision :: dfdy,rr,dydr,dfdr
      integer :: k
! Dispersion variables
      double precision :: dist(1),disp,dispdr(1)

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
      return
      end subroutine ev2gm2

      subroutine EvdVdR
      use N4_1A_ZV_par
!**********************************************************************
! Subroutine to evaluate dVdR for given R 
! dVdR = dV2dR + C*dBdR
! C:		Coefficients, stored in 'dim.inc' 
! P:		Basis functions evaluated for given R
! M:		Monomials evaluated for given R
! dV2dR:        Gradient of 2-body interactions
! dMsdR:	dMsdR(6,6), 6 MEG terms w.r.t. 6 bond length
! dMdR:		dMdR(6,nom), nom monomials w.r.t.6 bond length
! dPdR:		dPdR(6,nob), nop polynomial basis functions 
!		w.r.t. 6 bond length
! nom:  	number of monomials
! nob:  	number of basis functions(polynomials)
! M(nom):	Array to store monomials
! P(nob):	Array to store polynomials
!**********************************************************************
      
      integer i,j
      double precision dist,v2,dv2dr

! Initialize dVdR(6)
      do i=1,6
        dVdR(i)=0.0d0
      enddo

! Add dV2dR(i) to dVdR
      do i=1,6
        dist=R(i)
        call ev2gm2(dist,v2,dv2dr,1,1)
        dVdR(i)=dv2dr
      enddo

! Calculate dMEG/dr(6,6) for given R(6)
      call evdmsdr

! Calculate the monomials for each point by using six MEG terms
      call evdmdr

! Calculate the polynomials by using monomials
      call evdpdr 

! Remove 2-body interactions and unconnected terms from polynomials
      call evdbdr

! Evaluate dVdR(6) by taken the product of C(j) and dPdR(i,j)
      do i=1,6      
        do j=1,276
         dVdR(i)=dVdR(i) + c(j)*dBdR(i,j)
        enddo
      enddo

      return
      end subroutine EvdVdR

      subroutine EvdRdX(X)
      use N4_1A_ZV_par
!**********************************************************************
! Subroutine to evaluate dRdX for given R and X 
! R:		R(6), 6 bond lengths
! X:		X(12), 12 Cartesian coordinates
! 
! dMdR:		dMdR(6,nom), nom monomials w.r.t.6 bond length
! dPdR:		dPdR(6,nob), nop polynomial basis functions 
!		w.r.t. 6 bond length
! M(nom):	Array to store monomials
! P(nob):	Array to store polynomials
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
      use N4_1A_ZV_par
!**********************************************************************
! Subroutine to evaluate the derivatives of MEG term X
! w.r.t. interatomic distance R(6)
! dmsdR:	Local variables, dirm(6,6)
! a:		Nonlinear parameter(Angstrom)
! ab:       Nonlinear parameter(Angstrom^2)
! re:		equilibrium bond length(Angstrom)
!**********************************************************************

      integer i,j

! Initialize dmsdr
      do i=1,6
        do j=1,6
          dmsdr(i,j)=0.0d0
        enddo
      enddo

! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
! dmsdr(i,j)=0	i!=j

      do i=1,6
         dmsdr(i,i)=(-2.0d0*(r(i)-re)/ab-1/a) &
      * Exp(-(r(i)-re)/a-((r(i)-re)**2.0d0)/ab)
 
      enddo 

      return

      end subroutine EvdMsdR

      subroutine EvdMdR
      use N4_1A_ZV_par
!**********************************************************************
!  The subroutine reads M(nom) and dMSdR(6,6) and calculates the
!  dMdR(6,nom) that do not have usable decomposition.
!  For A4 with max. degree 9, the number of monomials is nom.
!**********************************************************************

      integer i

      do i=1,6
        dmdr(i,0) = 0.0d0
        dmdr(i,1) = dmsdr(i,6)
        dmdr(i,2) = dmsdr(i,5)
        dmdr(i,3) = dmsdr(i,4)
        dmdr(i,4) = dmsdr(i,3)
        dmdr(i,5) = dmsdr(i,2)
        dmdr(i,6) = dmsdr(i,1)
        dmdr(i,7) = dmdr(i,3)*rm(4) + rm(3)*dmdr(i,4)
        dmdr(i,8) = dmdr(i,2)*rm(5) + rm(2)*dmdr(i,5)
        dmdr(i,9) = dmdr(i,1)*rm(6) + rm(1)*dmdr(i,6)
        dmdr(i,10) = dmdr(i,1)*rm(2) + rm(1)*dmdr(i,2)
        dmdr(i,11) = dmdr(i,1)*rm(3) + rm(1)*dmdr(i,3)
        dmdr(i,12) = dmdr(i,2)*rm(3) + rm(2)*dmdr(i,3)
        dmdr(i,13) = dmdr(i,1)*rm(4) + rm(1)*dmdr(i,4)
        dmdr(i,14) = dmdr(i,2)*rm(4) + rm(2)*dmdr(i,4)
        dmdr(i,15) = dmdr(i,1)*rm(5) + rm(1)*dmdr(i,5)
        dmdr(i,16) = dmdr(i,3)*rm(5) + rm(3)*dmdr(i,5)
        dmdr(i,17) = dmdr(i,4)*rm(5) + rm(4)*dmdr(i,5)
        dmdr(i,18) = dmdr(i,2)*rm(6) + rm(2)*dmdr(i,6)
        dmdr(i,19) = dmdr(i,3)*rm(6) + rm(3)*dmdr(i,6)
        dmdr(i,20) = dmdr(i,4)*rm(6) + rm(4)*dmdr(i,6)
        dmdr(i,21) = dmdr(i,5)*rm(6) + rm(5)*dmdr(i,6)
        dmdr(i,22) = dmdr(i,1)*rm(7) + rm(1)*dmdr(i,7)
        dmdr(i,23) = dmdr(i,2)*rm(7) + rm(2)*dmdr(i,7)
        dmdr(i,24) = dmdr(i,1)*rm(8) + rm(1)*dmdr(i,8)
        dmdr(i,25) = dmdr(i,2)*rm(16) + rm(2)*dmdr(i,16)
        dmdr(i,26) = dmdr(i,2)*rm(17) + rm(2)*dmdr(i,17)
        dmdr(i,27) = dmdr(i,3)*rm(17) + rm(3)*dmdr(i,17)
        dmdr(i,28) = dmdr(i,1)*rm(18) + rm(1)*dmdr(i,18)
        dmdr(i,29) = dmdr(i,1)*rm(19) + rm(1)*dmdr(i,19)
        dmdr(i,30) = dmdr(i,1)*rm(20) + rm(1)*dmdr(i,20)
        dmdr(i,31) = dmdr(i,3)*rm(20) + rm(3)*dmdr(i,20)
        dmdr(i,32) = dmdr(i,1)*rm(21) + rm(1)*dmdr(i,21)
        dmdr(i,33) = dmdr(i,2)*rm(21) + rm(2)*dmdr(i,21)
        dmdr(i,34) = dmdr(i,1)*rm(12) + rm(1)*dmdr(i,12)
        dmdr(i,35) = dmdr(i,1)*rm(17) + rm(1)*dmdr(i,17)
        dmdr(i,36) = dmdr(i,2)*rm(20) + rm(2)*dmdr(i,20)
        dmdr(i,37) = dmdr(i,3)*rm(21) + rm(3)*dmdr(i,21)
        dmdr(i,38) = dmdr(i,1)*rm(14) + rm(1)*dmdr(i,14)
        dmdr(i,39) = dmdr(i,1)*rm(16) + rm(1)*dmdr(i,16)
        dmdr(i,40) = dmdr(i,2)*rm(19) + rm(2)*dmdr(i,19)
        dmdr(i,41) = dmdr(i,4)*rm(21) + rm(4)*dmdr(i,21)
        dmdr(i,42) = dmdr(i,2)*rm(27) + rm(2)*dmdr(i,27)
        dmdr(i,43) = dmdr(i,1)*rm(31) + rm(1)*dmdr(i,31)
        dmdr(i,44) = dmdr(i,1)*rm(33) + rm(1)*dmdr(i,33)
        dmdr(i,45) = dmdr(i,1)*rm(23) + rm(1)*dmdr(i,23)
        dmdr(i,46) = dmdr(i,1)*rm(25) + rm(1)*dmdr(i,25)
        dmdr(i,47) = dmdr(i,1)*rm(26) + rm(1)*dmdr(i,26)
        dmdr(i,48) = dmdr(i,1)*rm(27) + rm(1)*dmdr(i,27)
        dmdr(i,49) = dmdr(i,1)*rm(40) + rm(1)*dmdr(i,40)
        dmdr(i,50) = dmdr(i,1)*rm(36) + rm(1)*dmdr(i,36)
        dmdr(i,51) = dmdr(i,2)*rm(31) + rm(2)*dmdr(i,31)
        dmdr(i,52) = dmdr(i,1)*rm(37) + rm(1)*dmdr(i,37)
        dmdr(i,53) = dmdr(i,2)*rm(37) + rm(2)*dmdr(i,37)
        dmdr(i,54) = dmdr(i,1)*rm(41) + rm(1)*dmdr(i,41)
        dmdr(i,55) = dmdr(i,2)*rm(41) + rm(2)*dmdr(i,41)
        dmdr(i,56) = dmdr(i,3)*rm(41) + rm(3)*dmdr(i,41)
        dmdr(i,57) = dmdr(i,1)*rm(42) + rm(1)*dmdr(i,42)
        dmdr(i,58) = dmdr(i,1)*rm(51) + rm(1)*dmdr(i,51)
        dmdr(i,59) = dmdr(i,1)*rm(53) + rm(1)*dmdr(i,53)
        dmdr(i,60) = dmdr(i,1)*rm(55) + rm(1)*dmdr(i,55)
        dmdr(i,61) = dmdr(i,1)*rm(56) + rm(1)*dmdr(i,56)
        dmdr(i,62) = dmdr(i,2)*rm(56) + rm(2)*dmdr(i,56)
        dmdr(i,63) = dmdr(i,1)*rm(62) + rm(1)*dmdr(i,62)
        dmdr(i,64) = dmdr(i,2)*rm(57) + rm(2)*dmdr(i,57)
        dmdr(i,65) = dmdr(i,3)*rm(57) + rm(3)*dmdr(i,57)
        dmdr(i,66) = dmdr(i,4)*rm(57) + rm(4)*dmdr(i,57)
        dmdr(i,67) = dmdr(i,5)*rm(57) + rm(5)*dmdr(i,57)
        dmdr(i,68) = dmdr(i,1)*rm(58) + rm(1)*dmdr(i,58)
        dmdr(i,69) = dmdr(i,3)*rm(58) + rm(3)*dmdr(i,58)
        dmdr(i,70) = dmdr(i,4)*rm(58) + rm(4)*dmdr(i,58)
        dmdr(i,71) = dmdr(i,1)*rm(59) + rm(1)*dmdr(i,59)
        dmdr(i,72) = dmdr(i,2)*rm(59) + rm(2)*dmdr(i,59)
        dmdr(i,73) = dmdr(i,1)*rm(60) + rm(1)*dmdr(i,60)
        dmdr(i,74) = dmdr(i,2)*rm(60) + rm(2)*dmdr(i,60)
        dmdr(i,75) = dmdr(i,1)*rm(61) + rm(1)*dmdr(i,61)
        dmdr(i,76) = dmdr(i,2)*rm(62) + rm(2)*dmdr(i,62)
        dmdr(i,77) = dmdr(i,3)*rm(61) + rm(3)*dmdr(i,61)
        dmdr(i,78) = dmdr(i,3)*rm(62) + rm(3)*dmdr(i,62)
        dmdr(i,79) = dmdr(i,4)*rm(61) + rm(4)*dmdr(i,61)
        dmdr(i,80) = dmdr(i,4)*rm(62) + rm(4)*dmdr(i,62)
        dmdr(i,81) = dmdr(i,5)*rm(59) + rm(5)*dmdr(i,59)
        dmdr(i,82) = dmdr(i,5)*rm(60) + rm(5)*dmdr(i,60)
        dmdr(i,83) = dmdr(i,5)*rm(62) + rm(5)*dmdr(i,62)
        dmdr(i,84) = dmdr(i,6)*rm(58) + rm(6)*dmdr(i,58)
        dmdr(i,85) = dmdr(i,6)*rm(59) + rm(6)*dmdr(i,59)
        dmdr(i,86) = dmdr(i,6)*rm(60) + rm(6)*dmdr(i,60)
        dmdr(i,87) = dmdr(i,6)*rm(61) + rm(6)*dmdr(i,61)
        dmdr(i,88) = dmdr(i,2)*rm(64) + rm(2)*dmdr(i,64)
        dmdr(i,89) = dmdr(i,3)*rm(65) + rm(3)*dmdr(i,65)
        dmdr(i,90) = dmdr(i,4)*rm(66) + rm(4)*dmdr(i,66)
        dmdr(i,91) = dmdr(i,5)*rm(67) + rm(5)*dmdr(i,67)
        dmdr(i,92) = dmdr(i,1)*rm(68) + rm(1)*dmdr(i,68)
        dmdr(i,93) = dmdr(i,3)*rm(69) + rm(3)*dmdr(i,69)
        dmdr(i,94) = dmdr(i,4)*rm(70) + rm(4)*dmdr(i,70)
        dmdr(i,95) = dmdr(i,1)*rm(71) + rm(1)*dmdr(i,71)
        dmdr(i,96) = dmdr(i,2)*rm(72) + rm(2)*dmdr(i,72)
        dmdr(i,97) = dmdr(i,1)*rm(73) + rm(1)*dmdr(i,73)
        dmdr(i,98) = dmdr(i,2)*rm(74) + rm(2)*dmdr(i,74)
        dmdr(i,99) = dmdr(i,1)*rm(75) + rm(1)*dmdr(i,75)
        dmdr(i,100) = dmdr(i,2)*rm(76) + rm(2)*dmdr(i,76)
        dmdr(i,101) = dmdr(i,3)*rm(77) + rm(3)*dmdr(i,77)
        dmdr(i,102) = dmdr(i,3)*rm(78) + rm(3)*dmdr(i,78)
        dmdr(i,103) = dmdr(i,4)*rm(79) + rm(4)*dmdr(i,79)
        dmdr(i,104) = dmdr(i,4)*rm(80) + rm(4)*dmdr(i,80)
        dmdr(i,105) = dmdr(i,5)*rm(81) + rm(5)*dmdr(i,81)
        dmdr(i,106) = dmdr(i,5)*rm(82) + rm(5)*dmdr(i,82)
        dmdr(i,107) = dmdr(i,5)*rm(83) + rm(5)*dmdr(i,83)
        dmdr(i,108) = dmdr(i,6)*rm(84) + rm(6)*dmdr(i,84)
        dmdr(i,109) = dmdr(i,6)*rm(85) + rm(6)*dmdr(i,85)
        dmdr(i,110) = dmdr(i,6)*rm(86) + rm(6)*dmdr(i,86)
        dmdr(i,111) = dmdr(i,6)*rm(87) + rm(6)*dmdr(i,87)
      enddo

      return

      end subroutine EvdMdR

      subroutine EvdPdr
      use N4_1A_ZV_par
!**********************************************************************
!  The subroutine reads monomials(m) and calculates the
!  permutation-invariant polynomials(p)
!  For A4 with max. degree 9, the number of polynomials is nob.
!**********************************************************************

      integer i

      do i=1,6
       dpdr(i,0) = dmdr(i,0)
       dpdr(i,1) = dmdr(i,1) + dmdr(i,2) + dmdr(i,3) + dmdr(i,4) &
             + dmdr(i,5) + dmdr(i,6)
       dpdr(i,2) = dmdr(i,7) + dmdr(i,8) + dmdr(i,9)
       dpdr(i,3) = dmdr(i,10) + dmdr(i,11) + dmdr(i,12) + dmdr(i,13) &
            + dmdr(i,14) + dmdr(i,15) + dmdr(i,16) + dmdr(i,17) &
            + dmdr(i,18) + dmdr(i,19) + dmdr(i,20) + dmdr(i,21)
       dpdr(i,4) = dpdr(i,1)*p(1) + p(1)*dpdr(i,1) &
             - dpdr(i,3) - dpdr(i,2) - dpdr(i,3) - dpdr(i,2)
       dpdr(i,5) = dmdr(i,22) + dmdr(i,23) + dmdr(i,24) + dmdr(i,25) &
            + dmdr(i,26) + dmdr(i,27) + dmdr(i,28) + dmdr(i,29) &
            + dmdr(i,30) + dmdr(i,31) + dmdr(i,32) + dmdr(i,33)
       dpdr(i,6) = dmdr(i,34) + dmdr(i,35) + dmdr(i,36) + dmdr(i,37)
       dpdr(i,7) = dmdr(i,38) + dmdr(i,39) + dmdr(i,40) + dmdr(i,41)
       dpdr(i,8) = dpdr(i,1)*p(2) + p(1)*dpdr(i,2) &
             - dpdr(i,5)
       dpdr(i,9) = dpdr(i,1)*p(3) + p(1)*dpdr(i,3) &
             - dpdr(i,6) - dpdr(i,7) - dpdr(i,5) - dpdr(i,6) &
             - dpdr(i,7) - dpdr(i,5) - dpdr(i,6) - dpdr(i,7)
       dpdr(i,10) = dpdr(i,1)*p(4) + p(1)*dpdr(i,4) &
             - dpdr(i,9) - dpdr(i,8)
       dpdr(i,11) = dmdr(i,42) + dmdr(i,43) + dmdr(i,44)
       dpdr(i,12) = dmdr(i,45) + dmdr(i,46) + dmdr(i,47) + dmdr(i,48) &
            + dmdr(i,49) + dmdr(i,50) + dmdr(i,51) + dmdr(i,52) &
            + dmdr(i,53) + dmdr(i,54) + dmdr(i,55) + dmdr(i,56)
       dpdr(i,13) = dpdr(i,2)*p(3) + p(2)*dpdr(i,3) &
             - dpdr(i,12)
       dpdr(i,14) = dpdr(i,1)*p(5) + p(1)*dpdr(i,5) &
             - dpdr(i,12) - dpdr(i,11) - dpdr(i,13) - dpdr(i,12) &
             - dpdr(i,11) - dpdr(i,11) - dpdr(i,11)
       dpdr(i,15) = dpdr(i,1)*p(6) + p(1)*dpdr(i,6) &
             - dpdr(i,12)
       dpdr(i,16) = dpdr(i,1)*p(7) + p(1)*dpdr(i,7) &
             - dpdr(i,12)
       dpdr(i,17) = dpdr(i,2)*p(2) + p(2)*dpdr(i,2) &
             - dpdr(i,11) - dpdr(i,11)
       dpdr(i,18) = dpdr(i,3)*p(3) + p(3)*dpdr(i,3) &
             - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16) &
             - dpdr(i,14) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) &
             - dpdr(i,16) - dpdr(i,14) - dpdr(i,12) - dpdr(i,11) &
             - dpdr(i,12) - dpdr(i,11)
       dpdr(i,19) = dpdr(i,2)*p(4) + p(2)*dpdr(i,4) &
             - dpdr(i,14)
       dpdr(i,20) = dpdr(i,3)*p(4) + p(3)*dpdr(i,4) &
             - dpdr(i,15) - dpdr(i,16) - dpdr(i,13)
       dpdr(i,21) = dpdr(i,1)*p(10) + p(1)*dpdr(i,10) &
             - dpdr(i,20) - dpdr(i,19)
       dpdr(i,22) = dmdr(i,57) + dmdr(i,58) + dmdr(i,59) + dmdr(i,60) &
            + dmdr(i,61) + dmdr(i,62)
       dpdr(i,23) = dpdr(i,1)*p(11) + p(1)*dpdr(i,11) &
             - dpdr(i,22)
       dpdr(i,24) = dpdr(i,2)*p(6) + p(2)*dpdr(i,6)
       dpdr(i,25) = dpdr(i,2)*p(7) + p(2)*dpdr(i,7)
       dpdr(i,26) = dpdr(i,1)*p(12) + p(1)*dpdr(i,12) &
             - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,22) &
             - dpdr(i,22) - dpdr(i,22)
       dpdr(i,27) = dpdr(i,2)*p(5) + p(2)*dpdr(i,5) &
             - dpdr(i,22) - dpdr(i,23) - dpdr(i,22)
       dpdr(i,28) = dpdr(i,3)*p(5) + p(3)*dpdr(i,5) &
             - dpdr(i,22) - dpdr(i,26) - dpdr(i,24) - dpdr(i,25) &
             - dpdr(i,23) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) &
             - dpdr(i,23) - dpdr(i,22) - dpdr(i,22)
       dpdr(i,29) = dpdr(i,3)*p(6) + p(3)*dpdr(i,6) &
             - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
       dpdr(i,30) = dpdr(i,3)*p(7) + p(3)*dpdr(i,7) &
             - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
       dpdr(i,31) = dpdr(i,2)*p(9) + p(2)*dpdr(i,9) &
             - dpdr(i,26) - dpdr(i,28)
       dpdr(i,32) = dpdr(i,1)*p(14) + p(1)*dpdr(i,14) &
             - dpdr(i,26) - dpdr(i,23) - dpdr(i,28)
       dpdr(i,33) = dpdr(i,4)*p(6) + p(4)*dpdr(i,6) &
             - dpdr(i,25)
       dpdr(i,34) = dpdr(i,4)*p(7) + p(4)*dpdr(i,7) &
             - dpdr(i,24)
       dpdr(i,35) = dpdr(i,1)*p(17) + p(1)*dpdr(i,17) &
             - dpdr(i,27)
       dpdr(i,36) = dpdr(i,1)*p(18) + p(1)*dpdr(i,18) &
             - dpdr(i,29) - dpdr(i,30) - dpdr(i,28)
       dpdr(i,37) = dpdr(i,2)*p(10) + p(2)*dpdr(i,10) &
             - dpdr(i,32)
       dpdr(i,38) = dpdr(i,3)*p(10) + p(3)*dpdr(i,10) &
             - dpdr(i,33) - dpdr(i,34) - dpdr(i,31)
       dpdr(i,39) = dpdr(i,1)*p(21) + p(1)*dpdr(i,21) &
             - dpdr(i,38) - dpdr(i,37)
       dpdr(i,40) = dmdr(i,63)
       dpdr(i,41) = dmdr(i,64) + dmdr(i,65) + dmdr(i,66) + dmdr(i,67) &
            + dmdr(i,68) + dmdr(i,69) + dmdr(i,70) + dmdr(i,71) &
            + dmdr(i,72) + dmdr(i,73) + dmdr(i,74) + dmdr(i,75) &
            + dmdr(i,76) + dmdr(i,77) + dmdr(i,78) + dmdr(i,79) &
            + dmdr(i,80) + dmdr(i,81) + dmdr(i,82) + dmdr(i,83) &
            + dmdr(i,84) + dmdr(i,85) + dmdr(i,86) + dmdr(i,87)
       dpdr(i,42) = dpdr(i,1)*p(22) + p(1)*dpdr(i,22) &
             - dpdr(i,40) - dpdr(i,41) - dpdr(i,40) - dpdr(i,40) &
             - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
       dpdr(i,43) = dpdr(i,2)*p(11) + p(2)*dpdr(i,11) &
             - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
       dpdr(i,44) = dpdr(i,2)*p(12) + p(2)*dpdr(i,12) &
             - dpdr(i,41)
       dpdr(i,45) = dpdr(i,3)*p(11) + p(3)*dpdr(i,11) &
             - dpdr(i,41)
       dpdr(i,46) = dpdr(i,5)*p(6) + p(5)*dpdr(i,6) &
             - dpdr(i,41)
       dpdr(i,47) = dpdr(i,5)*p(7) + p(5)*dpdr(i,7) &
             - dpdr(i,41)
       dpdr(i,48) = dpdr(i,6)*p(7) + p(6)*dpdr(i,7) &
             - dpdr(i,40) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
       dpdr(i,49) = dpdr(i,4)*p(11) + p(4)*dpdr(i,11) &
             - dpdr(i,42)
       dpdr(i,50) = dpdr(i,2)*p(15) + p(2)*dpdr(i,15) &
             - dpdr(i,46)
       dpdr(i,51) = dpdr(i,2)*p(16) + p(2)*dpdr(i,16) &
             - dpdr(i,47)
       dpdr(i,52) = dpdr(i,4)*p(12) + p(4)*dpdr(i,12) &
             - dpdr(i,41) - dpdr(i,50) - dpdr(i,51)
       dpdr(i,53) = dpdr(i,2)*p(14) + p(2)*dpdr(i,14) &
             - dpdr(i,42) - dpdr(i,49) - dpdr(i,42)
       dpdr(i,54) = dpdr(i,6)*p(6) + p(6)*dpdr(i,6) &
             - dpdr(i,42) - dpdr(i,42)
       dpdr(i,55) = dpdr(i,7)*p(7) + p(7)*dpdr(i,7) &
             - dpdr(i,42) - dpdr(i,42)
       dpdr(i,56) = dpdr(i,3)*p(17) + p(3)*dpdr(i,17) &
             - dpdr(i,44)
       dpdr(i,57) = dpdr(i,2)*p(18) + p(2)*dpdr(i,18) &
             - dpdr(i,48)
       dpdr(i,58) = dpdr(i,3)*p(14) + p(3)*dpdr(i,14) &
             - dpdr(i,41) - dpdr(i,52) - dpdr(i,46) - dpdr(i,47) &
             - dpdr(i,45) - dpdr(i,45)
       dpdr(i,59) = dpdr(i,6)*p(9) + p(6)*dpdr(i,9) &
             - dpdr(i,41) - dpdr(i,52) - dpdr(i,47)
       dpdr(i,60) = dpdr(i,7)*p(9) + p(7)*dpdr(i,9) &
             - dpdr(i,41) - dpdr(i,52) - dpdr(i,46)
       dpdr(i,61) = dpdr(i,2)*p(20) + p(2)*dpdr(i,20) &
             - dpdr(i,52) - dpdr(i,58)
       dpdr(i,62) = dpdr(i,1)*p(32) + p(1)*dpdr(i,32) &
             - dpdr(i,52) - dpdr(i,49) - dpdr(i,58)
       dpdr(i,63) = dpdr(i,6)*p(10) + p(6)*dpdr(i,10) &
             - dpdr(i,51)
       dpdr(i,64) = dpdr(i,7)*p(10) + p(7)*dpdr(i,10) &
             - dpdr(i,50)
       dpdr(i,65) = dpdr(i,2)*p(17) + p(2)*dpdr(i,17) &
             - dpdr(i,43)
       dpdr(i,66) = dpdr(i,3)*p(18) + p(3)*dpdr(i,18) &
             - dpdr(i,46) - dpdr(i,47) - dpdr(i,45) - dpdr(i,59) &
             - dpdr(i,60) - dpdr(i,58)
       dpdr(i,67) = dpdr(i,2)*p(19) + p(2)*dpdr(i,19) &
             - dpdr(i,49)
       dpdr(i,68) = dpdr(i,1)*p(36) + p(1)*dpdr(i,36) &
             - dpdr(i,59) - dpdr(i,60) - dpdr(i,58) - dpdr(i,57) &
             - dpdr(i,66) - dpdr(i,66)
       dpdr(i,69) = dpdr(i,2)*p(21) + p(2)*dpdr(i,21) &
             - dpdr(i,62)
       dpdr(i,70) = dpdr(i,3)*p(21) + p(3)*dpdr(i,21) &
             - dpdr(i,63) - dpdr(i,64) - dpdr(i,61)
       dpdr(i,71) = dpdr(i,1)*p(39) + p(1)*dpdr(i,39) &
             - dpdr(i,70) - dpdr(i,69)
       dpdr(i,72) = dpdr(i,40)*p(1) + p(40)*dpdr(i,1)
       dpdr(i,73) = dpdr(i,2)*p(22) + p(2)*dpdr(i,22) &
             - dpdr(i,72)
       dpdr(i,74) = dpdr(i,6)*p(11) + p(6)*dpdr(i,11)
       dpdr(i,75) = dpdr(i,7)*p(11) + p(7)*dpdr(i,11)
       dpdr(i,76) = dpdr(i,3)*p(22) + p(3)*dpdr(i,22) &
             - dpdr(i,72) - dpdr(i,74) - dpdr(i,75) - dpdr(i,72) &
             - dpdr(i,72) - dpdr(i,72)
       dpdr(i,77) = dmdr(i,88) + dmdr(i,89) + dmdr(i,90) + dmdr(i,91) &
             + dmdr(i,92) + dmdr(i,93) + dmdr(i,94) + dmdr(i,95) &
             + dmdr(i,96) + dmdr(i,97) + dmdr(i,98) + dmdr(i,99) &
             + dmdr(i,100) + dmdr(i,101) + dmdr(i,102) + dmdr(i,103) &
             + dmdr(i,104) + dmdr(i,105) + dmdr(i,106) + dmdr(i,107) &
             + dmdr(i,108) + dmdr(i,109) + dmdr(i,110) + dmdr(i,111)
       dpdr(i,78) = dpdr(i,1)*p(42) + p(1)*dpdr(i,42) &
             - dpdr(i,72) - dpdr(i,76)
       dpdr(i,79) = dpdr(i,5)*p(11) + p(5)*dpdr(i,11) &
             - dpdr(i,72) - dpdr(i,73) - dpdr(i,72)
       dpdr(i,80) = dpdr(i,2)*p(26) + p(2)*dpdr(i,26) &
             - dpdr(i,76) - dpdr(i,77)
       dpdr(i,81) = dpdr(i,6)*p(12) + p(6)*dpdr(i,12) &
             - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
       dpdr(i,82) = dpdr(i,7)*p(12) + p(7)*dpdr(i,12) &
             - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
       dpdr(i,83) = dpdr(i,8)*p(11) + p(8)*dpdr(i,11) &
             - dpdr(i,72)
       dpdr(i,84) = dpdr(i,6)*p(17) + p(6)*dpdr(i,17)
       dpdr(i,85) = dpdr(i,7)*p(17) + p(7)*dpdr(i,17)
       dpdr(i,86) = dpdr(i,9)*p(11) + p(9)*dpdr(i,11) &
             - dpdr(i,76) - dpdr(i,77)
       dpdr(i,87) = dpdr(i,2)*p(29) + p(2)*dpdr(i,29) &
             - dpdr(i,81)
       dpdr(i,88) = dpdr(i,6)*p(14) + p(6)*dpdr(i,14) &
             - dpdr(i,75) - dpdr(i,75)
       dpdr(i,89) = dpdr(i,2)*p(30) + p(2)*dpdr(i,30) &
             - dpdr(i,82)
       dpdr(i,90) = dpdr(i,7)*p(14) + p(7)*dpdr(i,14) &
             - dpdr(i,74) - dpdr(i,74)
       dpdr(i,91) = dpdr(i,1)*p(48) + p(1)*dpdr(i,48) &
             - dpdr(i,76) - dpdr(i,81) - dpdr(i,82)
       dpdr(i,92) = dpdr(i,10)*p(11) + p(10)*dpdr(i,11) &
             - dpdr(i,78)
       dpdr(i,93) = dpdr(i,2)*p(33) + p(2)*dpdr(i,33) &
             - dpdr(i,88)
       dpdr(i,94) = dpdr(i,2)*p(34) + p(2)*dpdr(i,34) &
             - dpdr(i,90)
       dpdr(i,95) = dpdr(i,10)*p(12) + p(10)*dpdr(i,12) &
             - dpdr(i,77) - dpdr(i,93) - dpdr(i,94)
       dpdr(i,96) = dpdr(i,2)*p(27) + p(2)*dpdr(i,27) &
             - dpdr(i,73) - dpdr(i,79)
       dpdr(i,97) = dpdr(i,2)*p(28) + p(2)*dpdr(i,28) &
             - dpdr(i,76) - dpdr(i,86)
       dpdr(i,98) = dpdr(i,1)*p(53) + p(1)*dpdr(i,53) &
             - dpdr(i,80) - dpdr(i,79) - dpdr(i,97)
       dpdr(i,99) = dpdr(i,1)*p(54) + p(1)*dpdr(i,54) &
             - dpdr(i,81)
       dpdr(i,100) = dpdr(i,1)*p(55) + p(1)*dpdr(i,55) &
             - dpdr(i,82)
       dpdr(i,101) = dpdr(i,5)*p(18) + p(5)*dpdr(i,18) &
             - dpdr(i,76) - dpdr(i,91) - dpdr(i,87) - dpdr(i,89) &
             - dpdr(i,86)
       dpdr(i,102) = dpdr(i,6)*p(18) + p(6)*dpdr(i,18) &
             - dpdr(i,74) - dpdr(i,90)
       dpdr(i,103) = dpdr(i,7)*p(18) + p(7)*dpdr(i,18) &
             - dpdr(i,75) - dpdr(i,88)
       dpdr(i,104) = dpdr(i,2)*p(31) + p(2)*dpdr(i,31) &
             - dpdr(i,77) - dpdr(i,86)
       dpdr(i,105) = dpdr(i,2)*p(36) + p(2)*dpdr(i,36) &
             - dpdr(i,91) - dpdr(i,101)
       dpdr(i,106) = dpdr(i,3)*p(32) + p(3)*dpdr(i,32) &
             - dpdr(i,77) - dpdr(i,95) - dpdr(i,88) - dpdr(i,90) &
             - dpdr(i,86)
       dpdr(i,107) = dpdr(i,4)*p(29) + p(4)*dpdr(i,29) &
             - dpdr(i,82) - dpdr(i,80) - dpdr(i,99)
       dpdr(i,108) = dpdr(i,4)*p(30) + p(4)*dpdr(i,30) &
             - dpdr(i,81) - dpdr(i,80) - dpdr(i,100)
       dpdr(i,109) = dpdr(i,2)*p(38) + p(2)*dpdr(i,38) &
             - dpdr(i,95) - dpdr(i,106)
       dpdr(i,110) = dpdr(i,1)*p(62) + p(1)*dpdr(i,62) &
             - dpdr(i,95) - dpdr(i,92) - dpdr(i,106)
       dpdr(i,111) = dpdr(i,6)*p(21) + p(6)*dpdr(i,21) &
             - dpdr(i,94)
       dpdr(i,112) = dpdr(i,7)*p(21) + p(7)*dpdr(i,21) &
             - dpdr(i,93)
       dpdr(i,113) = dpdr(i,1)*p(65) + p(1)*dpdr(i,65) &
             - dpdr(i,96)
       dpdr(i,114) = dpdr(i,1)*p(66) + p(1)*dpdr(i,66) &
             - dpdr(i,102) - dpdr(i,103) - dpdr(i,101)
       dpdr(i,115) = dpdr(i,2)*p(37) + p(2)*dpdr(i,37) &
             - dpdr(i,92)
       dpdr(i,116) = dpdr(i,10)*p(18) + p(10)*dpdr(i,18) &
             - dpdr(i,99) - dpdr(i,100) - dpdr(i,97)
       dpdr(i,117) = dpdr(i,2)*p(39) + p(2)*dpdr(i,39) &
             - dpdr(i,110)
       dpdr(i,118) = dpdr(i,3)*p(39) + p(3)*dpdr(i,39) &
             - dpdr(i,111) - dpdr(i,112) - dpdr(i,109)
       dpdr(i,119) = dpdr(i,1)*p(71) + p(1)*dpdr(i,71) &
             - dpdr(i,118) - dpdr(i,117)
       dpdr(i,120) = dpdr(i,40)*p(2) + p(40)*dpdr(i,2)
       dpdr(i,121) = dpdr(i,40)*p(3) + p(40)*dpdr(i,3)
       dpdr(i,122) = dpdr(i,40)*p(4) + p(40)*dpdr(i,4)
       dpdr(i,123) = dpdr(i,11)*p(12) + p(11)*dpdr(i,12) &
             - dpdr(i,121)
       dpdr(i,124) = dpdr(i,2)*p(42) + p(2)*dpdr(i,42) &
             - dpdr(i,122)
       dpdr(i,125) = dpdr(i,6)*p(22) + p(6)*dpdr(i,22) &
             - dpdr(i,121)
       dpdr(i,126) = dpdr(i,7)*p(22) + p(7)*dpdr(i,22) &
             - dpdr(i,121)
       dpdr(i,127) = dpdr(i,2)*p(41) + p(2)*dpdr(i,41) &
             - dpdr(i,121) - dpdr(i,123) - dpdr(i,121)
       dpdr(i,128) = dpdr(i,6)*p(23) + p(6)*dpdr(i,23) &
             - dpdr(i,123)
       dpdr(i,129) = dpdr(i,7)*p(23) + p(7)*dpdr(i,23) &
             - dpdr(i,123)
       dpdr(i,130) = dpdr(i,3)*p(41) + p(3)*dpdr(i,41) &
             - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,125) &
             - dpdr(i,128) - dpdr(i,126) - dpdr(i,129) - dpdr(i,124) &
             - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) &
             - dpdr(i,125) - dpdr(i,126) - dpdr(i,124) - dpdr(i,123) &
             - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,122) &
             - dpdr(i,121) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120) &
             - dpdr(i,120) - dpdr(i,120)
       dpdr(i,131) = dpdr(i,3)*p(42) + p(3)*dpdr(i,42) &
             - dpdr(i,121) - dpdr(i,125) - dpdr(i,126) - dpdr(i,121)
       dpdr(i,132) = dpdr(i,4)*p(41) + p(4)*dpdr(i,41) &
             - dpdr(i,121) - dpdr(i,131) - dpdr(i,128) - dpdr(i,129) &
             - dpdr(i,127) - dpdr(i,121)
       dpdr(i,133) = dpdr(i,1)*p(78) + p(1)*dpdr(i,78) &
             - dpdr(i,122) - dpdr(i,131)
       dpdr(i,134) = dpdr(i,11)*p(11) + p(11)*dpdr(i,11) &
             - dpdr(i,120) - dpdr(i,120)
       dpdr(i,135) = dpdr(i,2)*p(48) + p(2)*dpdr(i,48) &
             - dpdr(i,130)
       dpdr(i,136) = dpdr(i,11)*p(17) + p(11)*dpdr(i,17) &
             - dpdr(i,120)
       dpdr(i,137) = dpdr(i,2)*p(44) + p(2)*dpdr(i,44) &
             - dpdr(i,123)
       dpdr(i,138) = dpdr(i,2)*p(45) + p(2)*dpdr(i,45) &
             - dpdr(i,121)
       dpdr(i,139) = dpdr(i,11)*p(14) + p(11)*dpdr(i,14) &
             - dpdr(i,122) - dpdr(i,124) - dpdr(i,122)
       dpdr(i,140) = dpdr(i,6)*p(27) + p(6)*dpdr(i,27) &
             - dpdr(i,127)
       dpdr(i,141) = dpdr(i,2)*p(54) + p(2)*dpdr(i,54)
       dpdr(i,142) = dpdr(i,7)*p(27) + p(7)*dpdr(i,27) &
             - dpdr(i,127)
       dpdr(i,143) = dpdr(i,2)*p(52) + p(2)*dpdr(i,52) &
             - dpdr(i,131) - dpdr(i,132)
       dpdr(i,144) = dpdr(i,1)*p(81) + p(1)*dpdr(i,81) &
             - dpdr(i,125) - dpdr(i,141) - dpdr(i,135) - dpdr(i,125)
       dpdr(i,145) = dpdr(i,2)*p(55) + p(2)*dpdr(i,55)
       dpdr(i,146) = dpdr(i,1)*p(82) + p(1)*dpdr(i,82) &
             - dpdr(i,126) - dpdr(i,135) - dpdr(i,145) - dpdr(i,126)
       dpdr(i,147) = dpdr(i,11)*p(18) + p(11)*dpdr(i,18) &
             - dpdr(i,130)
       dpdr(i,148) = dpdr(i,6)*p(28) + p(6)*dpdr(i,28) &
             - dpdr(i,129) - dpdr(i,123) - dpdr(i,143)
       dpdr(i,149) = dpdr(i,7)*p(28) + p(7)*dpdr(i,28) &
             - dpdr(i,128) - dpdr(i,123) - dpdr(i,143)
       dpdr(i,150) = dpdr(i,6)*p(30) + p(6)*dpdr(i,30) &
             - dpdr(i,121) - dpdr(i,146)
       dpdr(i,151) = dpdr(i,11)*p(19) + p(11)*dpdr(i,19) &
             - dpdr(i,122)
       dpdr(i,152) = dpdr(i,2)*p(50) + p(2)*dpdr(i,50) &
             - dpdr(i,128)
       dpdr(i,153) = dpdr(i,2)*p(51) + p(2)*dpdr(i,51) &
             - dpdr(i,129)
       dpdr(i,154) = dpdr(i,11)*p(20) + p(11)*dpdr(i,20) &
             - dpdr(i,131) - dpdr(i,132)
       dpdr(i,155) = dpdr(i,2)*p(59) + p(2)*dpdr(i,59) &
             - dpdr(i,144) - dpdr(i,148)
       dpdr(i,156) = dpdr(i,6)*p(32) + p(6)*dpdr(i,32) &
             - dpdr(i,129)
       dpdr(i,157) = dpdr(i,2)*p(60) + p(2)*dpdr(i,60) &
             - dpdr(i,146) - dpdr(i,149)
       dpdr(i,158) = dpdr(i,7)*p(32) + p(7)*dpdr(i,32) &
             - dpdr(i,128)
       dpdr(i,159) = dpdr(i,6)*p(34) + p(6)*dpdr(i,34) &
             - dpdr(i,122) - dpdr(i,145) - dpdr(i,122)
       dpdr(i,160) = dpdr(i,11)*p(21) + p(11)*dpdr(i,21) &
             - dpdr(i,133)
       dpdr(i,161) = dpdr(i,2)*p(63) + p(2)*dpdr(i,63) &
             - dpdr(i,156)
       dpdr(i,162) = dpdr(i,2)*p(64) + p(2)*dpdr(i,64) &
             - dpdr(i,158)
       dpdr(i,163) = dpdr(i,12)*p(21) + p(12)*dpdr(i,21) &
             - dpdr(i,132) - dpdr(i,161) - dpdr(i,162)
       dpdr(i,164) = dpdr(i,2)*p(53) + p(2)*dpdr(i,53) &
             - dpdr(i,124) - dpdr(i,139)
       dpdr(i,165) = dpdr(i,2)*p(58) + p(2)*dpdr(i,58) &
             - dpdr(i,131) - dpdr(i,154)
       dpdr(i,166) = dpdr(i,3)*p(54) + p(3)*dpdr(i,54) &
             - dpdr(i,125) - dpdr(i,144)
       dpdr(i,167) = dpdr(i,3)*p(55) + p(3)*dpdr(i,55) &
             - dpdr(i,126) - dpdr(i,146)
       dpdr(i,168) = dpdr(i,3)*p(65) + p(3)*dpdr(i,65) &
             - dpdr(i,137)
       dpdr(i,169) = dpdr(i,17)*p(18) + p(17)*dpdr(i,18) &
             - dpdr(i,135)
       dpdr(i,170) = dpdr(i,1)*p(98) + p(1)*dpdr(i,98) &
             - dpdr(i,143) - dpdr(i,139) - dpdr(i,165)
       dpdr(i,171) = dpdr(i,4)*p(54) + p(4)*dpdr(i,54) &
             - dpdr(i,135)
       dpdr(i,172) = dpdr(i,4)*p(55) + p(4)*dpdr(i,55) &
             - dpdr(i,135)
       dpdr(i,173) = dpdr(i,2)*p(66) + p(2)*dpdr(i,66) &
             - dpdr(i,150)
       dpdr(i,174) = dpdr(i,1)*p(101) + p(1)*dpdr(i,101) &
             - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,173) &
             - dpdr(i,165) - dpdr(i,147)
       dpdr(i,175) = dpdr(i,1)*p(102) + p(1)*dpdr(i,102) &
             - dpdr(i,150) - dpdr(i,148) - dpdr(i,166)
       dpdr(i,176) = dpdr(i,1)*p(103) + p(1)*dpdr(i,103) &
             - dpdr(i,150) - dpdr(i,149) - dpdr(i,167)
       dpdr(i,177) = dpdr(i,2)*p(61) + p(2)*dpdr(i,61) &
             - dpdr(i,132) - dpdr(i,154)
       dpdr(i,178) = dpdr(i,2)*p(68) + p(2)*dpdr(i,68) &
             - dpdr(i,159) - dpdr(i,174)
       dpdr(i,179) = dpdr(i,3)*p(62) + p(3)*dpdr(i,62) &
             - dpdr(i,132) - dpdr(i,163) - dpdr(i,156) - dpdr(i,158) &
             - dpdr(i,154)
       dpdr(i,180) = dpdr(i,6)*p(38) + p(6)*dpdr(i,38) &
             - dpdr(i,132) - dpdr(i,163) - dpdr(i,157)
       dpdr(i,181) = dpdr(i,7)*p(38) + p(7)*dpdr(i,38) &
             - dpdr(i,132) - dpdr(i,163) - dpdr(i,155)
       dpdr(i,182) = dpdr(i,2)*p(70) + p(2)*dpdr(i,70) &
             - dpdr(i,163) - dpdr(i,179)
       dpdr(i,183) = dpdr(i,1)*p(110) + p(1)*dpdr(i,110) &
             - dpdr(i,163) - dpdr(i,160) - dpdr(i,179)
       dpdr(i,184) = dpdr(i,6)*p(39) + p(6)*dpdr(i,39) &
             - dpdr(i,162)
       dpdr(i,185) = dpdr(i,7)*p(39) + p(7)*dpdr(i,39) &
             - dpdr(i,161)
       dpdr(i,186) = dpdr(i,2)*p(65) + p(2)*dpdr(i,65) &
             - dpdr(i,136)
       dpdr(i,187) = dpdr(i,3)*p(66) + p(3)*dpdr(i,66) &
             - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,175) &
             - dpdr(i,176) - dpdr(i,174)
       dpdr(i,188) = dpdr(i,2)*p(67) + p(2)*dpdr(i,67) &
             - dpdr(i,151)
       dpdr(i,189) = dpdr(i,4)*p(66) + p(4)*dpdr(i,66) &
             - dpdr(i,166) - dpdr(i,167) - dpdr(i,165)
       dpdr(i,190) = dpdr(i,2)*p(69) + p(2)*dpdr(i,69) &
             - dpdr(i,160)
       dpdr(i,191) = dpdr(i,18)*p(21) + p(18)*dpdr(i,21) &
             - dpdr(i,171) - dpdr(i,172) - dpdr(i,169)
       dpdr(i,192) = dpdr(i,2)*p(71) + p(2)*dpdr(i,71) &
             - dpdr(i,183)
       dpdr(i,193) = dpdr(i,3)*p(71) + p(3)*dpdr(i,71) &
             - dpdr(i,184) - dpdr(i,185) - dpdr(i,182)
       dpdr(i,194) = dpdr(i,1)*p(119) + p(1)*dpdr(i,119) &
             - dpdr(i,193) - dpdr(i,192)
       dpdr(i,195) = dpdr(i,40)*p(5) + p(40)*dpdr(i,5)
       dpdr(i,196) = dpdr(i,40)*p(6) + p(40)*dpdr(i,6)
       dpdr(i,197) = dpdr(i,40)*p(7) + p(40)*dpdr(i,7)
       dpdr(i,198) = dpdr(i,40)*p(8) + p(40)*dpdr(i,8)
       dpdr(i,199) = dpdr(i,40)*p(9) + p(40)*dpdr(i,9)
       dpdr(i,200) = dpdr(i,40)*p(10) + p(40)*dpdr(i,10)
       dpdr(i,201) = dpdr(i,11)*p(22) + p(11)*dpdr(i,22) &
             - dpdr(i,195)
       dpdr(i,202) = dpdr(i,12)*p(22) + p(12)*dpdr(i,22) &
             - dpdr(i,196) - dpdr(i,197) - dpdr(i,195) - dpdr(i,196) &
             - dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - dpdr(i,197)
       dpdr(i,203) = dpdr(i,17)*p(22) + p(17)*dpdr(i,22) &
             - dpdr(i,198)
       dpdr(i,204) = dpdr(i,6)*p(43) + p(6)*dpdr(i,43)
       dpdr(i,205) = dpdr(i,7)*p(43) + p(7)*dpdr(i,43)
       dpdr(i,206) = dpdr(i,11)*p(26) + p(11)*dpdr(i,26) &
             - dpdr(i,199) - dpdr(i,202)
       dpdr(i,207) = dpdr(i,2)*p(76) + p(2)*dpdr(i,76) &
             - dpdr(i,199) - dpdr(i,202)
       dpdr(i,208) = dpdr(i,2)*p(78) + p(2)*dpdr(i,78) &
             - dpdr(i,200)
       dpdr(i,209) = dpdr(i,6)*p(41) + p(6)*dpdr(i,41) &
             - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
       dpdr(i,210) = dpdr(i,6)*p(42) + p(6)*dpdr(i,42) &
             - dpdr(i,197) - dpdr(i,197) - dpdr(i,197)
       dpdr(i,211) = dpdr(i,7)*p(41) + p(7)*dpdr(i,41) &
             - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
       dpdr(i,212) = dpdr(i,7)*p(42) + p(7)*dpdr(i,42) &
             - dpdr(i,196) - dpdr(i,196) - dpdr(i,196)
       dpdr(i,213) = dpdr(i,11)*p(29) + p(11)*dpdr(i,29) &
             - dpdr(i,209)
       dpdr(i,214) = dpdr(i,11)*p(30) + p(11)*dpdr(i,30) &
             - dpdr(i,211)
       dpdr(i,215) = dpdr(i,18)*p(22) + p(18)*dpdr(i,22) &
             - dpdr(i,199) - dpdr(i,213) - dpdr(i,214)
       dpdr(i,216) = dpdr(i,2)*p(77) + p(2)*dpdr(i,77) &
             - dpdr(i,199) - dpdr(i,206)
       dpdr(i,217) = dpdr(i,6)*p(49) + p(6)*dpdr(i,49) &
             - dpdr(i,205)
       dpdr(i,218) = dpdr(i,7)*p(49) + p(7)*dpdr(i,49) &
             - dpdr(i,204)
       dpdr(i,219) = dpdr(i,3)*p(77) + p(3)*dpdr(i,77) &
             - dpdr(i,200) - dpdr(i,199) - dpdr(i,198) - dpdr(i,209) &
             - dpdr(i,217) - dpdr(i,211) - dpdr(i,218) - dpdr(i,207) &
             - dpdr(i,204) - dpdr(i,205) - dpdr(i,200) - dpdr(i,199) &
             - dpdr(i,198) - dpdr(i,200) - dpdr(i,198) - dpdr(i,200) &
             - dpdr(i,198)
       dpdr(i,220) = dpdr(i,3)*p(78) + p(3)*dpdr(i,78) &
             - dpdr(i,199) - dpdr(i,210) - dpdr(i,212)
       dpdr(i,221) = dpdr(i,10)*p(41) + p(10)*dpdr(i,41) &
             - dpdr(i,199) - dpdr(i,220) - dpdr(i,217) - dpdr(i,218) &
             - dpdr(i,216)
       dpdr(i,222) = dpdr(i,1)*p(133) + p(1)*dpdr(i,133) &
             - dpdr(i,200) - dpdr(i,220)
       dpdr(i,223) = dpdr(i,11)*p(27) + p(11)*dpdr(i,27) &
             - dpdr(i,195) - dpdr(i,203)
       dpdr(i,224) = dpdr(i,1)*p(134) + p(1)*dpdr(i,134) &
             - dpdr(i,201)
       dpdr(i,225) = dpdr(i,2)*p(81) + p(2)*dpdr(i,81) &
             - dpdr(i,209)
       dpdr(i,226) = dpdr(i,2)*p(80) + p(2)*dpdr(i,80) &
             - dpdr(i,202) - dpdr(i,206)
       dpdr(i,227) = dpdr(i,2)*p(82) + p(2)*dpdr(i,82) &
             - dpdr(i,211)
       dpdr(i,228) = dpdr(i,2)*p(91) + p(2)*dpdr(i,91) &
             - dpdr(i,215) - dpdr(i,219)
       dpdr(i,229) = dpdr(i,11)*p(28) + p(11)*dpdr(i,28) &
             - dpdr(i,199) - dpdr(i,207)
       dpdr(i,230) = dpdr(i,6)*p(53) + p(6)*dpdr(i,53) &
             - dpdr(i,205)
       dpdr(i,231) = dpdr(i,5)*p(54) + p(5)*dpdr(i,54) &
             - dpdr(i,209)
       dpdr(i,232) = dpdr(i,7)*p(53) + p(7)*dpdr(i,53) &
             - dpdr(i,204)
       dpdr(i,233) = dpdr(i,7)*p(54) + p(7)*dpdr(i,54) &
             - dpdr(i,196)
       dpdr(i,234) = dpdr(i,5)*p(55) + p(5)*dpdr(i,55) &
             - dpdr(i,211)
       dpdr(i,235) = dpdr(i,6)*p(55) + p(6)*dpdr(i,55) &
             - dpdr(i,197)
       dpdr(i,236) = dpdr(i,11)*p(35) + p(11)*dpdr(i,35) &
             - dpdr(i,198)
       dpdr(i,237) = dpdr(i,6)*p(65) + p(6)*dpdr(i,65)
       dpdr(i,238) = dpdr(i,7)*p(65) + p(7)*dpdr(i,65)
       dpdr(i,239) = dpdr(i,2)*p(86) + p(2)*dpdr(i,86) &
             - dpdr(i,199) - dpdr(i,229)
       dpdr(i,240) = dpdr(i,1)*p(139) + p(1)*dpdr(i,139) &
             - dpdr(i,206) - dpdr(i,229) - dpdr(i,224)
       dpdr(i,241) = dpdr(i,17)*p(29) + p(17)*dpdr(i,29) &
             - dpdr(i,225)
       dpdr(i,242) = dpdr(i,2)*p(99) + p(2)*dpdr(i,99) &
             - dpdr(i,231)
       dpdr(i,243) = dpdr(i,17)*p(30) + p(17)*dpdr(i,30) &
             - dpdr(i,227)
       dpdr(i,244) = dpdr(i,2)*p(95) + p(2)*dpdr(i,95) &
             - dpdr(i,220) - dpdr(i,221)
       dpdr(i,245) = dpdr(i,4)*p(81) + p(4)*dpdr(i,81) &
             - dpdr(i,202) - dpdr(i,242) - dpdr(i,227)
       dpdr(i,246) = dpdr(i,2)*p(100) + p(2)*dpdr(i,100) &
             - dpdr(i,234)
       dpdr(i,247) = dpdr(i,4)*p(82) + p(4)*dpdr(i,82) &
             - dpdr(i,202) - dpdr(i,225) - dpdr(i,246)
       dpdr(i,248) = dpdr(i,11)*p(36) + p(11)*dpdr(i,36) &
             - dpdr(i,215) - dpdr(i,219)
       dpdr(i,249) = dpdr(i,2)*p(102) + p(2)*dpdr(i,102) &
             - dpdr(i,233)
       dpdr(i,250) = dpdr(i,6)*p(58) + p(6)*dpdr(i,58) &
             - dpdr(i,206) - dpdr(i,214) - dpdr(i,244) - dpdr(i,214)
       dpdr(i,251) = dpdr(i,2)*p(103) + p(2)*dpdr(i,103) &
             - dpdr(i,235)
       dpdr(i,252) = dpdr(i,7)*p(58) + p(7)*dpdr(i,58) &
             - dpdr(i,213) - dpdr(i,206) - dpdr(i,244) - dpdr(i,213)
       dpdr(i,253) = dpdr(i,1)*p(150) + p(1)*dpdr(i,150) &
             - dpdr(i,215) - dpdr(i,233) - dpdr(i,235)
       dpdr(i,254) = dpdr(i,11)*p(37) + p(11)*dpdr(i,37) &
             - dpdr(i,200)
       dpdr(i,255) = dpdr(i,2)*p(93) + p(2)*dpdr(i,93) &
             - dpdr(i,217)
       dpdr(i,256) = dpdr(i,2)*p(94) + p(2)*dpdr(i,94) &
             - dpdr(i,218)
       dpdr(i,257) = dpdr(i,11)*p(38) + p(11)*dpdr(i,38) &
             - dpdr(i,220) - dpdr(i,221)
       dpdr(i,258) = dpdr(i,2)*p(107) + p(2)*dpdr(i,107) &
             - dpdr(i,245) - dpdr(i,250)
       dpdr(i,259) = dpdr(i,6)*p(62) + p(6)*dpdr(i,62) &
             - dpdr(i,218)
       dpdr(i,260) = dpdr(i,2)*p(108) + p(2)*dpdr(i,108) &
             - dpdr(i,247) - dpdr(i,252)
       dpdr(i,261) = dpdr(i,7)*p(62) + p(7)*dpdr(i,62) &
             - dpdr(i,217)
       dpdr(i,262) = dpdr(i,6)*p(64) + p(6)*dpdr(i,64) &
             - dpdr(i,200) - dpdr(i,246) - dpdr(i,200)
       dpdr(i,263) = dpdr(i,11)*p(39) + p(11)*dpdr(i,39) &
             - dpdr(i,222)
       dpdr(i,264) = dpdr(i,2)*p(111) + p(2)*dpdr(i,111) &
             - dpdr(i,259)
       dpdr(i,265) = dpdr(i,2)*p(112) + p(2)*dpdr(i,112) &
             - dpdr(i,261)
       dpdr(i,266) = dpdr(i,12)*p(39) + p(12)*dpdr(i,39) &
             - dpdr(i,221) - dpdr(i,264) - dpdr(i,265)
       dpdr(i,267) = dpdr(i,2)*p(98) + p(2)*dpdr(i,98) &
             - dpdr(i,208) - dpdr(i,240)
       dpdr(i,268) = dpdr(i,6)*p(54) + p(6)*dpdr(i,54) &
             - dpdr(i,210)
       dpdr(i,269) = dpdr(i,7)*p(55) + p(7)*dpdr(i,55) &
             - dpdr(i,212)
       dpdr(i,270) = dpdr(i,2)*p(96) + p(2)*dpdr(i,96) &
             - dpdr(i,203) - dpdr(i,223)
       dpdr(i,271) = dpdr(i,2)*p(97) + p(2)*dpdr(i,97) &
             - dpdr(i,207) - dpdr(i,229)
       dpdr(i,272) = dpdr(i,2)*p(101) + p(2)*dpdr(i,101) &
             - dpdr(i,215) - dpdr(i,248)
       dpdr(i,273) = dpdr(i,2)*p(106) + p(2)*dpdr(i,106) &
             - dpdr(i,220) - dpdr(i,257)
       dpdr(i,274) = dpdr(i,6)*p(59) + p(6)*dpdr(i,59) &
             - dpdr(i,220) - dpdr(i,215) - dpdr(i,209)
       dpdr(i,275) = dpdr(i,7)*p(60) + p(7)*dpdr(i,60) &
             - dpdr(i,220) - dpdr(i,215) - dpdr(i,211)
       dpdr(i,276) = dpdr(i,5)*p(66) + p(5)*dpdr(i,66) &
             - dpdr(i,215) - dpdr(i,253) - dpdr(i,249) - dpdr(i,251) &
             - dpdr(i,248)
       dpdr(i,277) = dpdr(i,6)*p(66) + p(6)*dpdr(i,66) &
             - dpdr(i,213) - dpdr(i,252)
       dpdr(i,278) = dpdr(i,7)*p(66) + p(7)*dpdr(i,66) &
             - dpdr(i,214) - dpdr(i,250)
       dpdr(i,279) = dpdr(i,2)*p(104) + p(2)*dpdr(i,104) &
             - dpdr(i,216) - dpdr(i,239)
       dpdr(i,280) = dpdr(i,2)*p(105) + p(2)*dpdr(i,105) &
             - dpdr(i,219) - dpdr(i,248)
       dpdr(i,281) = dpdr(i,1)*p(170) + p(1)*dpdr(i,170) &
             - dpdr(i,244) - dpdr(i,240) - dpdr(i,273)
       dpdr(i,282) = dpdr(i,10)*p(54) + p(10)*dpdr(i,54) &
             - dpdr(i,227)
       dpdr(i,283) = dpdr(i,10)*p(55) + p(10)*dpdr(i,55) &
             - dpdr(i,225)
       dpdr(i,284) = dpdr(i,2)*p(114) + p(2)*dpdr(i,114) &
             - dpdr(i,253) - dpdr(i,276)
       dpdr(i,285) = dpdr(i,1)*p(174) + p(1)*dpdr(i,174) &
             - dpdr(i,250) - dpdr(i,252) - dpdr(i,248) - dpdr(i,276) &
             - dpdr(i,273)
       dpdr(i,286) = dpdr(i,6)*p(68) + p(6)*dpdr(i,68) &
             - dpdr(i,217) - dpdr(i,261) - dpdr(i,251)
       dpdr(i,287) = dpdr(i,7)*p(68) + p(7)*dpdr(i,68) &
             - dpdr(i,218) - dpdr(i,259) - dpdr(i,249)
       dpdr(i,288) = dpdr(i,2)*p(109) + p(2)*dpdr(i,109) &
             - dpdr(i,221) - dpdr(i,257)
       dpdr(i,289) = dpdr(i,2)*p(116) + p(2)*dpdr(i,116) &
             - dpdr(i,262) - dpdr(i,285)
       dpdr(i,290) = dpdr(i,3)*p(110) + p(3)*dpdr(i,110) &
             - dpdr(i,221) - dpdr(i,266) - dpdr(i,259) - dpdr(i,261) &
             - dpdr(i,257)
       dpdr(i,291) = dpdr(i,6)*p(70) + p(6)*dpdr(i,70) &
             - dpdr(i,221) - dpdr(i,266) - dpdr(i,260)
       dpdr(i,292) = dpdr(i,7)*p(70) + p(7)*dpdr(i,70) &
             - dpdr(i,221) - dpdr(i,266) - dpdr(i,258)
       dpdr(i,293) = dpdr(i,2)*p(118) + p(2)*dpdr(i,118) &
             - dpdr(i,266) - dpdr(i,290)
       dpdr(i,294) = dpdr(i,1)*p(183) + p(1)*dpdr(i,183) &
             - dpdr(i,266) - dpdr(i,263) - dpdr(i,290)
       dpdr(i,295) = dpdr(i,6)*p(71) + p(6)*dpdr(i,71) &
             - dpdr(i,265)
       dpdr(i,296) = dpdr(i,7)*p(71) + p(7)*dpdr(i,71) &
             - dpdr(i,264)
       dpdr(i,297) = dpdr(i,1)*p(186) + p(1)*dpdr(i,186) &
             - dpdr(i,270)
       dpdr(i,298) = dpdr(i,1)*p(187) + p(1)*dpdr(i,187) &
             - dpdr(i,277) - dpdr(i,278) - dpdr(i,276)
       dpdr(i,299) = dpdr(i,2)*p(115) + p(2)*dpdr(i,115) &
             - dpdr(i,254)
       dpdr(i,300) = dpdr(i,1)*p(189) + p(1)*dpdr(i,189) &
             - dpdr(i,286) - dpdr(i,287) - dpdr(i,285) - dpdr(i,284) &
             - dpdr(i,298)
       dpdr(i,301) = dpdr(i,2)*p(117) + p(2)*dpdr(i,117) &
             - dpdr(i,263)
       dpdr(i,302) = dpdr(i,18)*p(39) + p(18)*dpdr(i,39) &
             - dpdr(i,282) - dpdr(i,283) - dpdr(i,280)
       dpdr(i,303) = dpdr(i,2)*p(119) + p(2)*dpdr(i,119) &
             - dpdr(i,294)
       dpdr(i,304) = dpdr(i,3)*p(119) + p(3)*dpdr(i,119) &
             - dpdr(i,295) - dpdr(i,296) - dpdr(i,293)
       dpdr(i,305) = dpdr(i,1)*p(194) + p(1)*dpdr(i,194) &
             - dpdr(i,304) - dpdr(i,303)
      enddo

      return

      end subroutine EvdPdR

      subroutine evdbdr
      use N4_1A_ZV_par
!**********************************************************************
!  The subroutine eliminates the 2-body terms in Bowman's approach
!**********************************************************************
      
      integer i,j
      double precision db1dr(6,306) 

! Pass P(0:305) to BM1(1:306)
      do j=1,6
      do i=1,306
        db1dr(j,i)=dpdr(j,i-1)
      enddo
      enddo

! Remove unconnected terms and 2-body terms and pass to B(1:276)
      do j=1,6
      dbdr(j,1)=db1dr(j,4)

      do i=2,4
        dbdr(j,i)=db1dr(j,i+4)
      enddo

      dbdr(j,5)=db1dr(j,10)

      do i=6,11
        dbdr(j,i)=db1dr(j,i+6)
      enddo

      dbdr(j,12)=db1dr(j,19)
      dbdr(j,13)=db1dr(j,21)

      do i=14,26
        dbdr(j,i)=db1dr(j,i+9)
      enddo

      dbdr(j,27)=db1dr(j,37)
      dbdr(j,28)=db1dr(j,39)

      do i=29,53
        dbdr(j,i)=db1dr(j,i+12)
      enddo

      dbdr(j,54)=db1dr(j,67)
      dbdr(j,55)=db1dr(j,69)
      dbdr(j,56)=db1dr(j,71)

      do i=57,97
        dbdr(j,i)=db1dr(j,i+16)
      enddo

      dbdr(j,98)=db1dr(j,115)
      dbdr(j,99)=db1dr(j,117)
      dbdr(j,100)=db1dr(j,119)

      do i=101,166
        dbdr(j,i)=db1dr(j,i+20)
      enddo

      dbdr(j,167)=db1dr(j,188)
      dbdr(j,168)=db1dr(j,190)
      dbdr(j,169)=db1dr(j,192)
      dbdr(j,170)=db1dr(j,194)

      do i=171,272
        dbdr(j,i)=db1dr(j,i+25)
      enddo

      dbdr(j,273)=db1dr(j,299)
      dbdr(j,274)=db1dr(j,301)
      dbdr(j,275)=db1dr(j,303)
      dbdr(j,276)=db1dr(j,305)
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

      if (imol.eq.1) then
! iz for N2 system
      iz(1)=7
      iz(2)=7

! C6 for N2 system
      c6(1)=19.7d0
      c6(2)=19.7d0
      else
! currently imol = 1 is used
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

