C   System:          H3
C   Functional form:
C   Common name:     A3
C   Number of derivatives: 0
C   Number of bodies: 3
C   Number of electronic surfaces: 1
C   Interface: 3-1V
C   Data file:
C
C   References:      SL Mielke, BC Garrett, and KA Peterson, J. Chem. Phys. 116 (2002) 4142.
C
C   Notes:    3 H-H distances passed in bohr. 
C             Energy returned in Hatrees relative to the classical minimum of H+H2.
C
      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad

      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: r(1,3), e(1), v
      integer :: iatom, idir, i, j
      logical, save :: first_time_data=.true.

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      ! input cartesian is HHH
      r(1,1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      r(1,2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      r(1,3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211

      call pot(r,e,1)

      v=e(1)
      v=v*27.211386

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=v
        enddo
      else
        write (*,*) 'Only energy is available'
      endif

      endsubroutine

      module a3data
      save
      double precision ::  eshift=0.1729928728793815d0
      double precision ::  beta=0.9291969140819419d0
      double precision ::  rnought=2.d0                         
      double precision, dimension(89) ::  xlin1=(/-0.1523291433408970d0,   
     &     0.8814508465229537d0,         6.537053395327828d0,   
     &     -11.34209491672260d0,         2.504271304238934d0,   
     &     4.326668451461330d0,         -2.239606610799320d0,   
     &     -0.4128984024042961d0,         0.5186887309355593d0,   
     &     -0.1196607910913041d0,         8.6003994509478154d-3,
     &     41.15414291845340d0,         86.45559795459729d0,   
     &     -116.3713616382037d0,         12.58972597209381d0,   
     &     33.81929163465899d0,         -17.15593990627275d0,   
     &     4.196478522411963d0,         -0.3411919636902630d0,   
     &     -3.4206236345359939d-3,       -127.6850556787161d0,   
     &     12.44205218121803d0,         -46.65469830727286d0,   
     &     2.132190766110241d0,         5.029968337261709d0,   
     &     -0.6386010902854538d0,         5.0081780443521871d-3,
     &     32.66461240039222d0,         -12.69768482032793d0,   
     &     0.8275601413171687d0,         3.9468614677319450d-2,
     &     1.5477194647200569d-2,       1.089948589115662d0,   
     &     0.7206200046869276d0,         -7.8425404380267669d-2,
     &     -0.1603878592389938d0,         -9.724970811924141d0,   
     &     135.5222874536157d0,         270.0021723964958d0,   
     &     -175.4825227217210d0,         -58.58906289532513d0,   
     &     100.5430564127615d0,         -37.00315174873428d0,   
     &     6.447449239215345d0,         8.2570662069006429d-3,
     &     -7.1098331653255367d-2,       715.3552761780156d0,   
     &     703.9102525501766d0,         -1861.078911265242d0,   
     &     -125.9587491921284d0,         225.8731414579132d0,   
     &     -46.07417109551250d0,         1.490246977236569d0,   
     &     0.2982638988343370d0,         -999.4775422429434d0,   
     &     151.7862715652945d0,         -133.0715723394372d0,   
     &     -12.30680794373364d0,         5.825782052452050d0,   
     &     -0.2439712905709224d0,         85.16524762327094d0,   
     &     16.46204547051492d0,         -7.133858657001830d0,   
     &     0.4030617759680488d0,         0.6909270689548164d0,   
     &     0.3097014450084296d0,         2556.785909933932d0,   
     &     5103.315112906826d0,         -2912.860957206383d0,   
     &     -160.8805563304034d0,         173.7542607534209d0,   
     &     -7.478371849288498d0,         -1.555852990003698d0,   
     &     2468.715016248928d0,         529.4607240124748d0,   
     &     -162.6956551745768d0,         -4.113114408714173d0,   
     &     1.611041727204823d0,         56.28848656173622d0,   
     &     9.735269264701914d0,         -3.377129269530819d0,   
     &     2.600888740672099d0,         -3246.118792587758d0,   
     &     161.4154798439835d0,         14.65512118402112d0,   
     &     1.684352485791373d0,         -19.56031027772914d0,   
     &     -8.9132714864690765d-2,       -0.6919443797821230d0/)  
      double precision, dimension(89) :: xlin2=(/2.6159085145117750d-2,
     &     1.2241854330557397d-2,     -0.6763682764316550d0,       
     &     -0.1562836764561804d0,       1.996594730253928d0,       
     &     -1.707513510844818d0,       0.2797466754569802d0,       
     &     0.2576622535966453d0,       -0.1404344725580116d0,       
     &     2.6052658461384465d-2,     -1.6693444722953859d-3,
     &     -0.4833014947701142d0,       -11.17866172577437d0,       
     &     7.607092641774579d0,       9.519929742191135d0,       
     &     -12.37037283407236d0,       4.707623579732490d0,       
     &     -0.5730444600084309d0,       -3.4653598431171197d-2,
     &     7.8187537452263211d-3,     38.72674320125338d0,       
     &     -7.269232169491375d0,       -9.879886268080352d0,       
     &     8.517458481793820d0,       -0.8199496374616673d0,       
     &     -0.1559314546076156d0,       1.5401566524384906d-2,
     &     9.554678972364389d0,       3.201632506258832d0,       
     &     -2.321941631684996d0,       3.6930680231572893d-2,
     &     3.2742627057071150d-2,     1.900429427415868d0,       
     &     -0.3327189910649105d0,       -7.7188825506846991d-3,
     &     9.3760580872625787d-2,     -1.014102607634918d0,       
     &     22.02809848224747d0,       -24.76473768346871d0,       
     &     -37.02983718150585d0,       58.80030375084864d0,       
     &     -28.47495308394478d0,       5.757365917529928d0,       
     &     4.3246175682466526d-2,     -0.1676517786924019d0,       
     &     1.5097351831581923d-2,     135.0679946550904d0,       
     &     -213.7000109608329d0,       95.89118538821953d0,       
     &     25.40918010829430d0,       9.782035588385490d0,       
     &     -4.168396476529826d0,       0.6883577777053669d0,       
     &     -6.9367683522399748d-2,     -262.5899913042082d0,       
     &     -53.22977186145765d0,       6.121343202849875d0,       
     &     7.823177388578925d0,       -0.2711421320419219d0,       
     &     -0.1043327957891995d0,       53.77598190726489d0,       
     &     -12.01741605017266d0,       -1.4784757119912589d-2,
     &     0.1488112444684154d0,       2.732678494168803d0,       
     &     -0.1970131259131940d0,       -467.4973230906062d0,       
     &     -444.5887032753010d0,       -315.4834996683269d0,       
     &     48.67913126874383d0,       24.84508509899326d0,       
     &     -4.729140122994481d0,       0.2069558855324623d0,       
     &     781.9447144772034d0,       -143.3043797424338d0,       
     &     -28.68693796003896d0,       2.344077954027283d0,       
     &     0.4918333516655770d0,       72.15250134833730d0,       
     &     -2.105811641390979d0,       -0.8842156350073892d0,       
     &     0.7824630993546668d0,       303.3073086352659d0,       
     &     -32.34761963775958d0,       -2.157141352189826d0,       
     &     -0.8000872307506535d0,       5.610906605940581d0,       
     &     1.199667001896809d0,       -4.678892380269374d0/)  
      double precision, dimension(71) ::  xlin3=(/5.3353477682006475d-2,      
     &    -0.1599144960966421d0,         -4.312643718679910d0,         
     &    6.814315854917867d0,         -2.205456828188728d0,         
     &    -0.7143749272312716d0,         -0.4114395291526620d0,         
     &    0.9520707372222905d0,         -0.3745853003354750d0,         
     &    4.2811012132827980d-2,       -21.36048998715569d0,         
     &    -74.47880934502300d0,         59.54491453184337d0,         
     &    26.38093015799090d0,         -40.49367250732830d0,         
     &    15.47625249844092d0,         -3.376506055551386d0,         
     &    0.1062027161738525d0,         18.81967136271011d0,         
     &    71.16187431300146d0,         14.71335233143087d0,         
     &    6.721378533743240d0,         -3.339779478411798d0,         
     &    -0.2339343234862702d0,         -72.80454829701864d0,         
     &    30.21027407834825d0,         -1.905603230726993d0,         
     &    -0.3996844993912013d0,         2.076285530599899d0,         
     &    -1.070380536340887d0,         8.598507152662780d0,         
     &    -73.61703213544804d0,         -219.8778318595011d0,         
     &    48.95327228138328d0,         109.3741906195230d0,         
     &    -104.0415972823556d0,         35.11211365107425d0,         
     &    -6.534860269104033d0,         -8.0422603679405902d-3,      
     &    -600.3417941867114d0,         -878.4454407436597d0,         
     &    1243.686176699445d0,         591.3052045918021d0,         
     &    -243.0727894093502d0,         39.64018613563970d0,         
     &    -0.8090321241234495d0,         753.6354030176653d0,         
     &    592.5340756904949d0,         109.1519077002292d0,         
     &    3.933746884407120d0,         3.506148188290761d0,         
     &    -189.3502887433559d0,         13.94268601946511d0,         
     &    0.6424807674312493d0,         1.419771258623844d0,         
     &    -2050.644716844225d0,         -5080.398907091004d0,         
     &    2822.341530025133d0,         638.1501839778731d0,         
     &    -217.6432578060511d0,         9.608610580227475d0,         
     &    -4953.965720751209d0,         -460.4461493597798d0,         
     &    40.61534110405616d0,         -18.26218983170062d0,         
     &    -113.4926963135811d0,         3.382375857465263d0,         
     &    2169.169337193152d0,         154.5462714936868d0,         
     &    84.83291904037159d0,         -71.12679286953023d0/)        
      double precision, dimension(71) :: xlin4=(/-9.4876856817530143d-3,      
     &    -3.6483362893877616d-2,       0.4205707974558235d0,         
     &    -0.1059572058101154d0,         -0.3707179832723506d0,         
     &    -0.1445341683023923d0,         0.5958052755316040d0,         
     &    -0.3854655745032128d0,         0.1010458670628552d0,         
     &    -9.4255692193045623d-3,       0.1509885666028168d0,         
     &    7.682780170568845d0,         -0.6840277749251158d0,         
     &    -11.98139522620621d0,         9.946678161114002d0,         
     &    -2.590064247380303d0,         7.6410732889741748d-2,      
     &    4.3207950587928659d-2,       -15.82021851395566d0,         
     &    -11.55029253337579d0,         9.357730866212668d0,         
     &    -1.817290225301895d0,         -2.051176763692164d0,         
     &    0.3422211063996889d0,         5.935560909793717d0,         
     &    -5.899083969710368d0,         -0.3644702094335688d0,         
     &    0.5301830356858073d0,         -2.119504022080795d0,         
     &    0.4144161928859211d0,         0.2636535169916772d0,         
     &    -9.256727433362473d0,         1.010605951519258d0,         
     &    48.54597955205001d0,         -48.87059486081539d0,         
     &    17.98526173595513d0,         -2.237557599107648d0,         
     &    -0.3789995719107901d0,         8.6650123871130774d-2,      
     &    -119.5839332405538d0,         133.3734002478926d0,         
     &    -28.89035273972550d0,         -56.60471823462028d0,         
     &    -9.1459472557995802d-2,       -2.754447481273484d0,         
     &    9.6332311722992170d-2,       253.4638121163468d0,         
     &    112.6549501199819d0,         -2.630187545871118d0,         
     &    -11.07273872864548d0,         0.3646509060305644d0,         
     &    -68.65504615698455d0,         10.86627678597302d0,         
     &    -0.6358917122975185d0,         -2.352160870433293d0,         
     &    84.30979852066358d0,         797.2622595018731d0,         
     &    214.5994942620052d0,         30.59962513407079d0,         
     &    -31.32236829922999d0,         2.833941401102188d0,         
     &    -476.2563698124825d0,         122.6747471786113d0,         
     &    17.34552327162847d0,         -0.8126973285879644d0,         
     &    -38.30796451594332d0,         1.661856665041447d0,         
     &    -1019.748561415903d0,         64.84858716655161d0,         
     &    5.856777175810681d0,         -8.373832093263392d0/)        
      end module
      subroutine prepot
      implicit none
      integer          :: ivp,nt
      double precision :: rvp(nt,3),evp(nt),r1x,r2x,r3x
      call prepota3
      return

      entry pot(rvp, evp, nt)
c
      do ivp = 1, nt
         r1x=rvp(ivp,1)  
         r2x=rvp(ivp,2)  
         r3x=rvp(ivp,3)  
         call pota3(r1x,r2x,r3x,evp(ivp))
      enddo    
c
      return
      end

      subroutine prepota3 
      use a3data 
      implicit double precision (a-h,o-z)
      save
      parameter (n2=5000)
      common /indy/iind(n2),jind(n2),kind(n2),ipar(n2),maxi,maxi2
      double precision :: rho1t(0:12),rho2t(0:12),rho3t(0:12)

      write(6,*)' Using A3 H3 PES of 8/15/01'
      write(6,*)' S. L. Mielke, B. C. Garrett, and K. A. Peterson, J. Ch
     &em. Phys., 116 (2002) 4142'

      epslon=1.d-12
      epscem=1.d-12
      ald=0.02d0
      betad=0.72d0

      maxi2=0
      maxi=0
      call indexa3(12,12)     
      maxi2a=maxi2

      maxi=0
      call indexa3(11,11)        
      maxi2b=maxi2-maxi2a
      return

      entry pota3(r1x,r2x,r3x,eee)
      call trip(r1x,etrip1)
      call trip(r2x,etrip2)
      call trip(r3x,etrip3)
      call singlet(r1x,es1)
      call singlet(r2x,es2)
      call singlet(r3x,es3)
      q1=0.5d0*(es1+etrip1)
      q2=0.5d0*(es2+etrip2)
      q3=0.5d0*(es3+etrip3)
      xj1=0.5d0*(es1-etrip1)
      xj2=0.5d0*(es2-etrip2)
      xj3=0.5d0*(es3-etrip3)
      xj=sqrt(epslon+0.5d0*((xj3-xj1)**2+(xj2-xj1)**2+(xj3-xj2)**2))    
      xjs=0.25d0*(xj1+xj2+xj3)
      vlondon=q1+q2+q3-xj

      rho1=exp(-beta*(r1x-rnought))
      rho2=exp(-beta*(r2x-rnought))
      rho3=exp(-beta*(r3x-rnought))

      if((r1x.le.betad).or.(r2x.le.betad).or.(r3x.le.betad))then
        damper=0.d0
      else
        rdamp1=ald/(betad-r1x)
        rdamp2=ald/(betad-r2x)
        rdamp3=ald/(betad-r3x)
        damper=exp(rdamp1+rdamp2+rdamp3)
      endif

      bb=sqrt(epscem+(r1x-r3x)**2+(r1x-r2x)**2+(r2x-r3x)**2)
      warp=1.d0/(1.d0/r1x+1.d0/r2x+1.d0/r3x)

      v=eshift+vlondon

      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      sum4=0.d0
      
      rho1t(0)=1.d0
      rho2t(0)=1.d0
      rho3t(0)=1.d0

      do it=1,12
        rho1t(it)=rho1t(it-1)*rho1
        rho2t(it)=rho2t(it-1)*rho2
        rho3t(it)=rho3t(it-1)*rho3
      enddo    

      do ii=1,maxi2a
        rprod=rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
        sum1=sum1+xlin1(ipar(ii))*rprod
        sum2=sum2+xlin2(ipar(ii))*rprod
      enddo     

      do ii=maxi2a+1,maxi2a+maxi2b
        rprod=rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
        sum3=sum3+xlin3(ipar(ii))*rprod
        sum4=sum4+xlin4(ipar(ii))*rprod
      enddo     

      eee=v+damper*(sum1+sum2*bb+sum3*warp+sum4*bb*warp)

      return
      end
      
      subroutine indexa3(mbig,mtop)            
      parameter (n2=5000)
c     calculate the fitting indices for A3 symmetry
      common /indy/iind(n2),jind(n2),kind(n2),ipar(n2),maxi,maxi2

      l=maxi2
      l2=maxi

      do  i=0,min(mtop,mbig)
       do  j=i,min(mtop,mbig)
        do  k=j,min(mtop,mbig)
         isum=i+j+k
          if(isum.le.mbig)then
           if(((isum-i).gt.0).and.((isum-j).gt.0).and.
     &        ((isum-k).gt.0))then     
            l2=l2+1
            l=l+1
            ipar(l)=l2
            iind(l)=i
            jind(l)=j
            kind(l)=k
            if(i.ne.j)then
             if(j.ne.k)then
               l=l+1
               ipar(l)=l2
               iind(l)=i
               jind(l)=k
               kind(l)=j
               l=l+1
               ipar(l)=l2
               iind(l)=j
               jind(l)=i
               kind(l)=k
               l=l+1
               ipar(l)=l2
               iind(l)=j
               jind(l)=k
               kind(l)=i
               l=l+1
               ipar(l)=l2
               iind(l)=k
               jind(l)=i
               kind(l)=j
               l=l+1
               ipar(l)=l2
               iind(l)=k
               jind(l)=j
               kind(l)=i
              else
               l=l+1
               ipar(l)=l2
               iind(l)=j
               jind(l)=i
               kind(l)=k
               l=l+1
               ipar(l)=l2
               iind(l)=j
               jind(l)=k
               kind(l)=i
             endif
            elseif (j.ne.k)then
               l=l+1
               ipar(l)=l2
               iind(l)=j
               jind(l)=k
               kind(l)=i
               l=l+1
               ipar(l)=l2
               iind(l)=k
               jind(l)=j
               kind(l)=i
            endif
           endif
          endif
        enddo     
       enddo    
      enddo    
      maxi=l2
      maxi2=l
      return
      end

      subroutine trip(r,v)
c     H2 triplet curve for FCI/aug-cc-pVTZ
      implicit none
      integer          :: j
      double precision :: r,v,damp,xd,prefac
      double precision :: c6=-6.499027d0
      double precision :: c8=-124.3991d0
      double precision :: c10=-3285.828d0
      double precision :: beta=0.2d0
      double precision :: re=1.401d0

      double precision :: alpha=3.8716418156250960d0
      double precision,dimension(17) :: xlp=(/2.1948631629921860d-1,
     &        5.8892583097268760d-1, 8.6235600154107670d-1,
     &        8.6194702003339070d-1, 6.8604820043881530d-1,
     &        4.5436692359636230d-1, 2.3371369821182240d-1,
     &        6.9843424078195320d-2, -2.6494129220450960d-3,
     &        -8.5475462878786820d-3, 3.0498779190429310d-3,
     &        3.0448279356991980d-3, -5.6431515849655190d-4,
     &        -4.0765667872606770d-4, 1.9838486248362640d-4,
     &        -3.2129748801006750d-5, 1.9377307867602390d-6/)

      damp=1.d0-exp(-beta*r**2)
      xd=damp/r  
      v=c6*xd**6+c8*xd**8+c10*xd**10

      prefac=exp(- alpha*(r-re))
      if((r-re).ne.0.d0)then
        do j=1,17
          v=v+xlp(j)*(r-re)**(j-1)*prefac
        enddo    
      else
        v=v+xlp(1)*prefac
      endif

      return
      end

      subroutine singlet(r,v)
c     H2 singlet potential for FCI/aug-cc-pVTZ
      implicit none
      integer          :: j
      double precision :: r,v,damp,xd,prefac
      double precision :: c6=-6.499027d0
      double precision :: c8=-124.3991d0
      double precision :: c10=-3285.828d0
      double precision :: beta=0.2d0
      double precision :: re=1.401d0
      double precision :: alpha=  3.9259981290585980d0
      double precision, dimension(17) :: xlp= (/-1.6948148539493060d-1,
     &        -6.5341157509085330d-1, -1.0568246834574850d0,
     &        -1.0380175566115990d0, -6.7657643586978360d-1,
     &        -3.0434043205949370d-1, -1.0305186076403050d-1,
     &        -4.9128168199975780d-2, -3.7985765587970510d-2,
     &        -1.8881471279090120d-2, 4.4201285958012090d-4,
     &        3.2421117771421200d-3, -3.9290175622769260d-4,
     &        -5.7935611043841150d-4, 2.0970556978282260d-4,
     &        -2.8399222117118230d-5, 1.2735256346272810d-6/)


      damp=1.d0-exp(-beta*r**2)
      xd=damp/r 
      prefac=exp(- alpha*(r-re))
      v=c6*xd**6+c8*xd**8+c10*xd**10

      if((r-re).ne.0.d0)then
        do j=1,17
          v=v+xlp(j)*(r-re)**(j-1)*prefac
        enddo    
      else
        v=v+xlp(1)*prefac
      endif

      return
      end
