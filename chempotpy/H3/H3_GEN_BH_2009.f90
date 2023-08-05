C   System:          H3
C   Functional form:
C   Common name:     Born-Huang
C   Number of derivatives: 1
C   Number of bodies: 3
C   Number of electronic surfaces: 1
C   Interface: Section-2
C   Data file:
C
C   References:  Steven L. Mielke, David W. Schwenke, George C. Schatz, Bruce C. Garrett, and Kirk A. Peterson, J. Phys. Chem. A 113, 4479 (2009).
C                Steven L. Mielke, Bruce C. Garrett, and Kirk A. Peterson, J. Chem. Phys. 116, 4142 (2002).
C
C   Notes:  Version with analytical gradients. This is the Born-Oppenheimer diagonal correction PES published in
C           Steven L. Mielke, David W. Schwenke, George C. Schatz, Bruce C. Garrett, and Kirk A. Peterson, J. Phys. 
C           Chem. A, 2009, 113, 4479 combined with the CCI PES published in Steven L. Mielke, Bruce C. Garrett, 
C           and Kirk A. Peterson, J. Chem. Phys. 116, 4142 (2002).  This version is set up for the HHH mass combination. 
C           The Born-Huang PES is mass dependent and to select other mass combinations one should reset the parameters 
C           xmass1, xmass2, and xmass3 with the nuclear masses of the desired isotopes.  
C           Input and output is via the common block "common /potcm/ rvp,evp,dvp" where rvp(3) is a 3-vector, 
C           (R1, R2, R3) of HH distances, evp is the energy, and dvp(3) is the gradient, all in atomic units.  
C           R1 is the distance between mass1 and mass2, R2 is the distance between mass2 and mass3, and R3 is the 
C           distance between mass1 and mass3.
C
C
      module bodcdata
      save
      double precision :: xlin1(180)=  (/
     &      2.451839720665176D+02,   -3.023514674602468D+02,
     &     -3.643179260300581D+02,    1.033283164173005D+03,
     &     -9.649202796347731D+02,    4.723249503697386D+02,
     &     -1.312024315948558D+02,    2.053242829600808D+01,
     &     -1.663116010000327D+00,    5.317955346689944D-02,
     &      3.693835631800750D+02,    9.643897029819937D+01,
     &     -3.794338103064540D+02,    2.184383291329870D+02,
     &     -3.787026565563110D+01,   -3.960932914729867D+00,
     &      1.733261249407424D+00,   -1.233088325428533D-01,
     &      7.634113468160084D+02,   -6.089390321694615D+02,
     &      2.145121698941533D+02,   -4.303659182480859D+01,
     &      4.793100434454963D+00,   -2.362075732031712D-01,
     &      5.161649762181137D+02,   -1.379218427454547D+02,
     &      1.222214446630326D+01,   -1.700714347704069D-01,
     &      2.590291319801722D+01,   -1.187599837006957D+00,
     &      2.949665817120014D+02,   -4.423967394072385D+02,
     &     -1.499691959357717D+02,    8.615959233030945D+02,
     &     -9.063992109731515D+02,    4.800278837250469D+02,
     &     -1.444131347943261D+02,    2.489091872075205D+01,
     &     -2.288354420765698D+00,    8.706064188372947D-02,
     &     -1.561653042358118D+03,    2.795630390969673D+03,
     &     -4.442497231862677D+03,    2.653340911293203D+03,
     &     -7.274300003583212D+02,    1.223835539827257D+02,
     &     -1.767422177777780D+01,    1.080509007725006D+00,
     &      8.007412890755861D-02,   -6.866852114392939D+02,
     &      1.965474846931427D+03,   -8.327550428277026D+02,
     &      2.921969416962440D+00,    5.934336236873139D+01,
     &     -1.157242040507434D+01,    6.777877552427620D-01,
     &      6.809701703798451D+02,   -6.460103807549953D+02,
     &      1.080711387406044D+02,   -2.581355973993933D+00,
     &     -2.702457089997727D-01,    3.373040929271790D+02,
     &     -3.319320528532528D+01,   -8.252443297418132D-01,
     &      2.242813487013438D+00,   -4.786033544785251D+02,
     &      7.379465395733814D+02,   -2.997606193372937D+02,
     &     -1.444004106661323D+02,    1.160834932238466D+02,
     &     -1.055041157003008D+01,   -6.084935176995607D+00,
     &      1.255348914553237D+00,   -4.967414449352265D-02,
     &      2.869805514805477D+03,   -2.901158537312932D+02,
     &      2.593678864031728D+03,   -1.474577262776355D+03,
     &      1.739661050571814D+02,    4.108033331334562D+01,
     &     -9.031479084332030D+00,    3.265756035116496D-01,
     &     -9.394389739950026D+03,    1.281273536289386D+03,
     &      9.186863730169819D+02,   -2.224206878555418D+02,
     &      2.852507226985430D+01,   -2.975868406356395D+00,
     &     -2.887753159216629D+02,   -6.232442825355130D+01,
     &     -2.080038525889382D+01,    5.846990317948459D+00,
     &      1.494611359198648D+01,    1.430793238496208D+00,
     &     -1.161743802353405D+02,   -4.435211957384280D+02,
     &      7.949179810694009D+02,   -4.621878870840682D+02,
     &      1.829448335674127D+02,   -4.859832146779939D+01,
     &      6.467586280687250D+00,   -3.014121429251370D-01,
     &     -4.727165211936207D+03,    1.681945178807448D+03,
     &      6.001879815346814D+02,   -5.284816509082527D+02,
     &      9.509901063900256D+01,   -8.355776044216967D+00,
     &      4.667856406455721D-01,    5.648107714717939D+02,
     &     -2.902912088460810D+02,   -9.407789952750028D+01,
     &     -1.631433254869286D+00,    5.242491229250340D+00,
     &      3.214464002696617D+01,    2.913754508725327D+01,
     &     -1.568292873005413D+01,    5.462484623657479D+00,
     &      9.083635945770520D+02,    9.799890386730289D+01,
     &     -6.411746299039472D+02,    4.298716219171475D+02,
     &     -1.205657322203302D+02,    1.315814978109771D+01,
     &     -3.975978827444066D-01,    3.180124260914723D+03,
     &     -5.872907359889355D+02,   -6.259492461216178D+02,
     &      3.434767513553659D+02,   -3.640741128738960D+01,
     &     -9.236985549909328D-01,    1.155215128322806D+03,
     &     -9.784160598115349D+01,    1.777745091732398D+01,
     &      2.485226250947354D-01,    1.521755613790718D+01,
     &      7.156889276881836D+00,   -1.016142747224836D+03,
     &     -5.200008131651686D+01,    2.946853566528521D+02,
     &     -1.316509239684804D+02,    2.366837195107599D+01,
     &     -1.174923797943324D+00,   -1.159057124509925D+03,
     &     -1.521893512764090D+02,    1.104768477989274D+02,
     &     -3.656291398178829D+01,    3.488333934608939D+00,
     &     -2.286294041696426D+02,   -7.953044862181258D+00,
     &     -2.999788001060466D-02,   -1.546533722850765D+01,
     &      5.670669717065498D+02,    5.041037222513963D+01,
     &     -7.690800748197464D+01,    1.401968534663515D+01,
     &     -1.189136309799305D+00,    2.944715172811239D+02,
     &      1.011155288939896D+02,   -5.420733436869213D+00,
     &     -8.900491551445843D-01,    2.785072259806540D+01,
     &      6.092864636803840D+00,   -1.799957321893311D+02,
     &     -1.892894466520360D+01,    1.036160665441249D+01,
     &     -3.761882545122726D-01,   -5.436396122931865D+01,
     &     -1.702038609057893D+01,    4.969068890501028D-03,
     &     -3.519006024232350D+00,    3.284752951627190D+01,
     &      2.783127622573391D+00,   -5.624072554917517D-01,
     &      5.370237229208703D+00,    9.768755790120242D-01,
     &     -3.206502231938864D+00,   -1.267702449737564D-01,
     &     -1.430186406716300D-01,    1.296229968602821D-01/)
      double precision :: xlin2(180)=  (/
     &      8.145631479853698D-01,    1.770730113974052D+01,
     &     -2.173959339189124D+01,    1.961131724453491D+01,
     &     -4.104999264119586D+01,    2.337793724047610D+01,
     &     -5.449082199443552D+00,    6.201541014865459D-01,
     &     -4.728365816655994D-02,    3.321586948535380D-03,
     &     -6.652255446067896D+01,    5.962375314289732D+01,
     &      1.469998744381738D+01,    1.153758314494926D+00,
     &     -7.591897754372836D+00,    2.103874186533575D+00,
     &     -1.702965824754904D-01,    1.346878024842888D-04,
     &     -1.178260275080853D+02,    3.489556535478315D+01,
     &     -1.251243111695901D+00,    2.471519774632120D-01,
     &     -2.474793708434775D-01,    2.687242988983517D-02,
     &     -2.026760191030008D+01,    2.756556273693077D+00,
     &      8.357078625422904D-02,   -1.640242746249622D-02,
     &     -6.876043147213573D-01,    2.609337260717558D-02,
     &      3.002603577505805D+00,    1.498828705957069D+01,
     &     -2.086306249380632D+01,    2.081499551542078D+01,
     &     -3.171758822478321D+01,    1.723622638287100D+01,
     &     -4.028169514216751D+00,    4.671596528203490D-01,
     &     -3.670331530469922D-02,    2.629265006741015D-03,
     &     -8.927371067901906D+01,    8.867171901574056D+01,
     &     -1.133027069526800D+02,    1.602734361579233D+02,
     &     -1.023650687832928D+02,    2.471829619377287D+01,
     &     -2.300433949989794D+00,    2.592872235150836D-01,
     &     -4.379218716825992D-02,   -1.231186961261672D+01,
     &      1.420549968400411D+01,    1.232687494983029D+00,
     &      1.841045801282512D+01,   -6.141054432632874D+00,
     &      2.034847145705890D-01,    8.860867509074010D-02,
     &     -9.644503856306309D+01,    2.392130829529128D+01,
     &     -3.146354866365471D+00,    9.562642394369082D-01,
     &     -2.269962411748502D-01,   -7.627358272973419D+00,
     &      4.570619067451160D-01,    2.237747868045461D-01,
     &     -3.169396242810660D-01,    1.231481888044332D+01,
     &     -5.861438736123795D+01,    5.294643266448595D+01,
     &     -4.088837194494829D+00,    5.495540512301918D+00,
     &     -6.300138454926564D+00,    1.765098522094756D+00,
     &     -1.780074643720451D-01,    4.114464551103923D-03,
     &      8.917318947290573D+01,   -1.770180683611587D+01,
     &     -5.602678714811219D+00,    2.973755372611144D+01,
     &      5.605823953149696D+00,   -4.949227356210042D+00,
     &      2.420977258547833D-01,    8.167453315236323D-02,
     &     -1.135834657778073D+02,    8.290397544044251D+01,
     &     -7.296380654367582D+01,    1.347009356803447D+01,
     &     -2.500870639262219D-01,   -1.272556171177672D-01,
     &      5.669201615035288D+00,    4.449130216981585D+00,
     &     -1.272014612117903D+00,    2.139168528456700D-01,
     &     -4.777213854750650D-02,   -1.275863250945433D-02,
     &     -1.642934979019225D+01,    5.594030291218376D+01,
     &     -8.305010271115899D+01,    3.140219381847346D+01,
     &     -2.534163029719606D+00,    1.142479699516450D-01,
     &     -1.263578178350560D-01,    1.629095589290311D-02,
     &     -1.136985820244078D+02,    3.262096657063820D+01,
     &     -1.087121813266270D+02,    2.673092202326776D+01,
     &     -5.444160084123269D-01,    6.138628773543448D-01,
     &     -2.227712501792724D-01,    5.512545979043862D+01,
     &      2.559974311865491D+01,    5.216798906963263D-01,
     &     -1.496009262806457D+00,    2.638941938405173D-01,
     &     -9.723650088603947D+00,    6.868905539663074D-01,
     &     -4.177757286042882D-01,    2.781285313472535D-01,
     &      1.563554076363040D+01,   -6.542137033344432D+00,
     &      2.193791357437970D+01,   -1.576158545476312D+01,
     &      2.587945022189958D+00,   -3.893301430803128D-02,
     &     -5.190707866429820D-03,    1.367205185988990D+02,
     &      1.188103897312893D+01,    1.606164459906229D+01,
     &     -9.262867579726493D+00,    6.271894032241021D-01,
     &      2.285179285451767D-01,   -7.103397905922931D+01,
     &      3.121061890947882D+00,    6.354928967335751D-01,
     &     -4.163368478300122D-02,    1.020366559611126D+00,
     &      1.837284088182978D-01,   -2.832630010141147D+01,
     &      1.275654589067907D+01,   -4.603547842975937D-01,
     &      1.995827205765568D+00,   -5.232664167379651D-01,
     &      1.710301952383997D-02,   -8.590072407323460D+01,
     &      1.674416687913584D+01,    1.953210020027019D-02,
     &      2.520910415921704D-01,   -3.423339521867145D-01,
     &      1.169816725323726D+01,   -1.761561723165250D+00,
     &      3.013931319831423D-02,   -4.165165498393316D-01,
     &      1.386201564606759D+01,   -9.091983474080589D+00,
     &      2.768047770222391D-01,    5.763259032772141D-02,
     &      3.029683243309628D-02,    1.740587394078216D+01,
     &     -6.199982916108434D+00,    9.489522270750911D-01,
     &      2.535531600295191D-01,    1.698179626605738D-01,
     &      2.305666871236848D-01,   -2.509462487790768D+00,
     &      2.075944787757823D+00,   -2.275822797369400D-01,
     &     -2.065529414930606D-02,   -1.074236086914455D+00,
     &      6.803905517919883D-02,   -2.615499139877834D-01,
     &     -1.335426493176628D-01,    1.840487213722926D-01,
     &     -1.569973762062688D-01,    3.050419120227986D-02,
     &      2.653719569186909D-01,    1.094965343816732D-01,
     &     -2.187339773412837D-02,   -1.873849846674096D-03,
     &     -5.603548924512692D-02,    3.558364313199398D-03/)
      double precision :: xlin3(180)=  (/
     &     -7.638557823809794D+01,    1.703707427310439D+02,
     &     -1.025601533591727D+02,   -5.296578652461962D+01,
     &      1.093595183012879D+02,   -6.588531176536077D+01,
     &      2.040793884905717D+01,   -3.479631418111812D+00,
     &      3.088169472711390D-01,   -1.111590540608995D-02,
     &     -3.361962954982916D+02,    2.446590342879648D+02,
     &     -6.622902002288096D+01,    1.191674304065106D+01,
     &     -1.201463499030519D+01,    6.122631636838099D+00,
     &     -1.187879235377073D+00,    7.926806077855969D-02,
     &     -5.725326950390917D+01,    2.868589802572544D+01,
     &     -2.971320475494586D+01,    1.008864240170261D+01,
     &     -1.365325434437133D+00,    6.554733610556563D-02,
     &     -4.045186273038085D+01,    1.874648094047193D+01,
     &     -2.454684686634876D+00,    6.765838416051127D-02,
     &     -4.059257274704468D+00,    2.165639317481252D-01,
     &     -7.534329078823575D+01,    1.796639204536162D+02,
     &     -1.564733447866279D+02,    3.775541937573816D+01,
     &      3.375720237731807D+01,   -3.005181622792514D+01,
     &      1.029486036217851D+01,   -1.791612509905809D+00,
     &      1.538953705090221D-01,   -5.019936328398713D-03,
     &     -4.774634443417375D+01,    3.863906310035196D+02,
     &      1.909693601494610D+00,    7.896215873407198D+01,
     &     -1.295233067421028D+02,    4.285454822181649D+01,
     &     -3.947563931965661D+00,   -2.453031697183805D-02,
     &     -4.072260116080702D-03,   -1.512501065486055D+03,
     &      9.514665008404671D+02,   -6.253456734186879D+02,
     &      2.584487207114875D+02,   -5.147932463240105D+01,
     &      4.978444306783625D+00,   -2.242912909135229D-01,
     &     -5.814167788685131D+02,    2.286668120955779D+02,
     &     -6.032832571577611D+01,    8.065674990175529D+00,
     &     -3.404150604733470D-01,   -3.061436485306763D+01,
     &      2.200637753633254D+00,    1.057897733012469D-01,
     &      7.289136325949297D-03,    1.896105176076998D+02,
     &     -3.648470785265107D+02,    2.595841400715867D+02,
     &     -6.319152600268254D+01,    8.494896397888191D+00,
     &     -1.053255418482307D+01,    5.467310061893341D+00,
     &     -1.031401599859863D+00,    6.598583555227435D-02,
     &      4.919698991437962D+02,   -1.523587257667702D+03,
     &      1.014167387686189D+03,   -6.508692392742124D+02,
     &      2.713378987371789D+02,   -5.735977155603892D+01,
     &      5.830582807234884D+00,   -2.440949592570337D-01,
     &     -6.735007851702401D+02,    7.901400404656945D+02,
     &     -1.229867349843578D+02,   -3.009484917235691D+01,
     &      2.526743872042880D+00,    5.970110295569673D-01,
     &     -3.532469213851834D+02,    4.005928931048262D+01,
     &      5.976702173098687D+00,   -1.614616921246194D+00,
     &     -9.871235251782917D+00,    3.290059681876524D-01,
     &     -1.813749415305320D+02,    3.109697012254384D+02,
     &     -3.205683056506578D+01,   -3.657410087464648D+00,
     &     -2.624482850677591D+01,    1.173325195979073D+01,
     &     -1.590994311096534D+00,    5.486195466674193D-02,
     &     -2.440497562681578D+02,    1.038109037697877D+03,
     &     -7.152807104622796D+02,    2.980486621854296D+02,
     &     -7.381873725420134D+01,    9.536084102100837D+00,
     &     -4.034487399636460D-01,    6.749359685299933D+02,
     &     -3.958160095292775D+02,    7.043339305226182D+01,
     &      4.280031318910807D+00,   -1.799431596685355D+00,
     &      1.572684765203245D+02,   -1.497852043495127D+01,
     &      3.450483851677588D+00,   -1.651601456234870D+00,
     &      6.267522188135509D+01,   -1.490721159969828D+02,
     &      7.203117438959657D+01,   -4.910294960766687D+01,
     &      2.126436044051092D+01,   -3.286263182633637D+00,
     &      1.432920523463851D-01,    3.482161123451701D+02,
     &     -6.634624347379950D+02,    2.978996288910582D+02,
     &     -6.102181144139994D+01,    4.913896432859977D+00,
     &      9.220877926952298D-02,   -9.040330715822043D+01,
     &      3.763878231734758D+01,   -1.466110603121045D+01,
     &      4.286621953679466D-01,   -1.171184819157715D+01,
     &     -8.282936408674097D-01,    2.329172506772613D+01,
     &      7.205432276781195D+01,   -7.134506264990834D+01,
     &      2.778750878256592D+01,   -4.744233867644683D+00,
     &      2.489849788505943D-01,   -3.035952094511171D+02,
     &      2.774672648114976D+02,   -7.050154798680367D+01,
     &      7.641437548457593D+00,   -1.880492973015346D-01,
     &     -1.984086397967394D+01,    7.604978254357780D+00,
     &      2.014923057356265D-01,    3.308108128880802D+00,
     &     -2.969993736353019D+01,   -3.682715000434093D+01,
     &      2.210331765988686D+01,   -4.301393898599036D+00,
     &      2.671370386619945D-01,    1.124784775141184D+02,
     &     -6.041199428997342D+01,    6.988933983482634D+00,
     &     -1.362655996129701D-01,   -2.997626634973072D-02,
     &     -1.846855278627695D+00,    1.154319604248443D+01,
     &      1.192127579658776D+01,   -2.699670781911932D+00,
     &      1.920744268281359D-01,   -1.990490891239322D+01,
     &      7.011127199328610D+00,   -9.502550138734989D-02,
     &      7.777433674389070D-01,   -2.235484930636709D+00,
     &     -1.909433565086661D+00,    1.029497435932162D-01,
     &      1.784299385800847D+00,   -3.937636905107279D-01,
     &      2.159590396303356D-01,    1.161199705458321D-01,
     &     -7.827610505544675D-02,   -8.180052967228596D-03/)
      end module

      module zmasses
        implicit none
        double precision :: zmass1, zmass2
      end module

!      subroutine prepot(newmasses)
      subroutine pot(rvp,evp,dvp)
      use zmasses
      implicit none
      double precision, intent(in) :: rvp(3)
      double precision, intent(out) :: evp, dvp(3)
      integer          :: icase, info, i
      double precision :: r1x,r2x,r3x, ecci
      double precision :: dcci(3),d1(3),d2(3),d3(3)
      double precision :: xmh, xmd, xmminv(3,3), amat(3), t
      double precision :: e1,e2,e3,emat(3), ezero, hatom_bodc,xmin
      double precision :: twobzero=-17.5988842162493135d0  ! lowest H2 2body bodc
      double precision :: ccishift= 0.1744759989649372d0,ezero2
      double precision :: towave=219474.7d0  ! hartree to cm-1 conversion

!   CASE for  H + H2
      double precision :: xmass1=1836.15264d0, xmass2=1836.15264d0, 
     &                    xmass3=1836.15264d0
      save xmh, xmd, xmass1, xmass2, xmass3, icase, xmminv, ezero, 
     &     towave

!     write(6,*)present(newmasses)
!     if(present(newmasses))then
!       xmass1=newmasses(1)
!       xmass2=newmasses(2)
!       xmass3=newmasses(3)
!     endif

!    We need to zero the PES with respect to one arrangement
!    If only one mass is distinct, choose the target mass as that, otherwise choose 
!    the lightest particle.  If another choice is desired, add appropriate assignment
!    statements here

      if(xmass1.eq.xmass2)then
        zmass1=xmass1
        zmass2=xmass2
      elseif(xmass1.eq.xmass3)then
        zmass1=xmass1
        zmass2=xmass3
      elseif(xmass2.eq.xmass3)then
        zmass1=xmass2
        zmass2=xmass3
      else
        zmass1=max(xmass1,xmass2,xmass3)
        if(xmass1.lt.zmass1)then
          zmass2=max(xmass1,min(xmass2,xmass3))
        else
          zmass2=max(xmass2,xmass3)
        endif
      endif

      xmh=1836.15264d0   ! H nuclear mass (after substracting 1 for electron mass)
      xmd=3670.48293d0   ! D nuclear mass (after substracting 1 for electron mass)
      hatom_bodc=59.6818d0  ! the BODC in cm-1 of H with an aug-cc-pVTZ basis set
 
      ezero=hatom_bodc*xmh*(1d0/xmass1+1d0/xmass2+1d0/xmass3)
      call golden(1.3d0,1.4d0,1.5d0,1.d-14,xmin,ezero2)

      ezero=ezero/towave+ezero2+ccishift

      if( (xmass1.eq.xmd) .and. (xmass2.eq.xmh) 
     &      .and. (xmass3.eq.xmh) )then
        icase=1  ! We are solving for D + H2 so we don't need to solve for amat
      else
        icase=-1 ! We are not solving for D + H2 so we need 3 times as many evaluations so we can determine amat
        t=(2.0d0*xmd+xmh)*(xmd-xmh)
        xmminv(:,:)=xmh*xmd*xmd/t
        t=-xmh*xmd*(xmd+xmh)/t
        xmminv(1,1)=t
        xmminv(2,2)=t
        xmminv(3,3)=t
      endif
 
      call prepotcbs
      call prepot_leps
      call prepotbodcdh2

! R1 is the distance between mass1 and mass2, R2 is the distance between mass2 and mass3,
! and R3 is the distance between mass1 and mass3.  Note that the D+H2 code requires the homonuclear
! diatom distance to be R3, so before calling all the various subroutines we need to switch R2 and R3
!

      r1x=rvp(1)
      r3x=rvp(2)
      r2x=rvp(3)
      call potcbs(r1x,r2x,r3x,ecci,dcci)

      call potbodcdh2(r1x,r2x,r3x,e1,d1)

      if(icase.eq.1)then
         evp=e1 - ezero
         dvp(:)=d1(:)
      else
         call potbodcdh2(r1x,r3x,r2x,e2,d2)
         t = d2(2)
         d2(2) = d2(3)
         d2(3) = t
         call potbodcdh2(r3x,r2x,r1x,e3,d3)
         t = d3(1)
         d3(1) = d3(3)
         d3(3) = t
         emat(1)=e1
         emat(2)=e2
         emat(3)=e3
         amat(:) = matmul(xmminv,emat)
         evp=amat(1)/xmass1+amat(2)/xmass2+amat(3)/xmass3 - ezero
         do i =1,3
            emat(1) = d1(i)
            emat(2) = d2(i)
            emat(3) = d3(i)
            amat(:) = matmul(xmminv,emat)
            dvp(i) = amat(1)/xmass1+amat(2)/xmass2+amat(3)/xmass3
         enddo
      endif

      evp=evp/towave+ecci-ezero
      dvp(1)=dvp(1)/towave+dcci(1)
      t =dvp(2)/towave+dcci(2)
      dvp(2)=dvp(3)/towave+dcci(3)
      dvp(3) = t

      end


      subroutine prepotbodcdh2
      use bodcdata
      implicit none
      save

      integer,parameter :: n2 = 5000
      common /indy2/iind(n2),jind(n2),kind(n2),ipar(n2),maxi,maxi2
      integer :: iind,jind,kind,ipar  ! arrays of length n2
      integer :: ii, it, maxi, maxi2, maxi2a, i
      double precision :: r1x, r2x, r3x, rho1, rho2, rho3, eee, dv(1:3)
      double precision :: drho1, drho2, drho3
      double precision :: e1body, xmd, xmh, hatom_bodc, bodc2body
      double precision :: round, ald, betad, beta, deldamp, rhodamp
      double precision :: bb2, bb3, bb4, bb5, v, rnought, towave
      double precision :: xma, xmb, xmc, v1, v2, gtot, rsum,
     &                    betadsum, rdampsum, v3c
      double precision :: sum1, sum2, sum3, sum4, sum5, rprod
      double precision :: bigq, slit, slitsq, damper, rpass(3),
     &                    sdamp, rdamp1, rdamp2, rdamp3
      double precision :: slit1, slit2, slit3, ddamper(1:3), dsdamp
      double precision :: dgtot(1:3),dbb2(1:3),dbb3(1:3)
      double precision :: dv1, dv2, dv3
      double precision :: dsum1(1:3), dsum2(1:3), dsum3(1:3)
      double precision :: rho1t(0:13),rho2t(0:13),rho3t(0:13)
      double precision :: drho1t(0:13),drho2t(0:13),drho3t(0:13)
      double precision :: drprod1, drprod2, drprod3


!     distances are in Bohr, energies are in cm-1

      towave=219474.7d0  ! hartree to cm-1 conversion
      xmh=1836.15264d0   ! H nuclear mass (after substracting 1 for electron mass)
      xmd=3670.48293d0   ! D nuclear mass (after substracting 1 for electron mass)
      hatom_bodc=59.6818d0  ! the BODC in cm-1 of H with an aug-cc-pVTZ basis set
      xma=xmd
      xmb=xmh
      xmc=xmh

      round=1.d-2  
      ald=0.02d0
      betad=0.72d0
      betadsum=3.3d0
      deldamp=0.25d0
      rhodamp=0.01d0
      beta=  1.414254132358239d0 
      rnought=2.0d0
      maxi2=0
      maxi=0
      call indexa2(11,11)     
      maxi2a=maxi2
      return


      entry potbodcdh2(r1x,r2x,r3x,eee,dv)

      rsum=r1x+r2x+r3x
      if((r1x.le.betad).or.(r2x.le.betad).or.(r3x.le.betad)
     &    .or.(rsum.le.betadsum) )then
        damper=0.0d0
        ddamper(:)=0.0d0
      else
        bigq=r1x*r1x+r2x*r2x+r3x*r3x
        slit1=2.0d0*r1x**2-r2x**2-r3x**2
        slit2=r2x**2-r3x**2
        slit3=slit1*slit1+3.0d0*slit2*slit2
        slit=sqrt(slit3)/bigq
        if(slit.le.deldamp)then
           damper=0.0d0
           ddamper(:)=0.0d0
        else
           sdamp=rhodamp/(deldamp-slit)
           rdampsum=ald/(betadsum-rsum)
           rdamp1=ald/(betad-r1x)
           rdamp2=ald/(betad-r2x)
           rdamp3=ald/(betad-r3x)
           damper=exp(rdampsum+rdamp1+rdamp2+rdamp3+sdamp)

           dsdamp=sdamp/(deldamp-slit)
           ddamper(1)= damper*(dsdamp*slit*r1x*
     &          (4.0d0*slit1/slit3 - 2.0d0/bigq)
     &          +rdampsum/(betadsum-rsum)+rdamp1/(betad-r1x))
           ddamper(2)=damper*(-dsdamp*slit*r2x*
     &          (4.0d0*(r1x*r1x-2.0d0*r2x*r2x+r3x*r3x)/slit3+2.0d0/bigq)
     &          +rdampsum/(betadsum-rsum)+rdamp2/(betad-r2x))
           ddamper(3)=damper*(-dsdamp*slit*r3x*
     &          (4.0d0*(r1x*r1x+r2x*r2x-2.0d0*r3x*r3x)/slit3+2.0d0/bigq)
     &          +rdampsum/(betadsum-rsum)+rdamp3/(betad-r3x))
        endif
      endif

      rpass(1)=r1x
      rpass(3)=r2x
      rpass(2)=r3x
      Call Leps(rpass,V1,V2,xma,xmb,xmc,gtot,dgtot)
      sum1 = dgtot(2)
      dgtot(2) = dgtot(3)
      dgtot(3) = sum1


      bb2=towave*gtot
      do i = 1,3
         dbb2(i) = towave*dgtot(i)
      enddo
      bb3=sqrt(round+(r1x-r3x)**2+(r1x-r2x)**2+(r2x-r3x)**2)
      dbb3(1)=(2.0d0*r1x-r2x-r3x)/bb3
      dbb3(2)=(2.0d0*r2x-r1x-r3x)/bb3
      dbb3(3)=(2.0d0*r3x-r1x-r2x)/bb3

      e1body = (2.d0 + xmh/xmd) * hatom_bodc 
 
      v = e1body + ((xmd+xmh)/(2.d0*xmd))*(bodc2body(r1x,dv1)
     &     +bodc2body(r2x,dv2)) + bodc2body(r3x,dv3)

      dv(1) = ((xmd+xmh)/(2.d0*xmd))*dv1
      dv(2) = ((xmd+xmh)/(2.d0*xmd))*dv2
      dv(3) = dv3

      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      dsum1(:)=0.0d0
      dsum2(:)=0.0d0
      dsum3(:)=0.0d0
      
      rho1=exp(-beta*(r1x-rnought))
      drho1=(1.0-beta*r1x)*rho1
      rho1=r1x*rho1
      rho2=exp(-beta*(r2x-rnought))
      drho2=(1.0-beta*r2x)*rho2
      rho2=r2x*rho2
      rho3=exp(-beta*(r3x-rnought))
      drho3=(1.0-beta*r3x)*rho3
      rho3=r3x*rho3

      rho1t(0)=1.d0
      rho2t(0)=1.d0
      rho3t(0)=1.d0

      do it=1,13
        rho1t(it)=rho1t(it-1)*rho1
        rho2t(it)=rho2t(it-1)*rho2
        rho3t(it)=rho3t(it-1)*rho3

        drho1t(it)=it*rho1t(it-1)*drho1
        drho2t(it)=it*rho2t(it-1)*drho2
        drho3t(it)=it*rho3t(it-1)*drho3
      enddo    

      do ii=1,maxi2a
        rprod=rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
        sum1=sum1+xlin1(ipar(ii))*rprod
        sum2=sum2+xlin2(ipar(ii))*rprod
        sum3=sum3+xlin3(ipar(ii))*rprod

        drprod1=drho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
        drprod2=rho1t(kind(ii))*drho2t(jind(ii))*rho3t(iind(ii))
        drprod3=rho1t(kind(ii))*rho2t(jind(ii))*drho3t(iind(ii))
        dsum1(1)=dsum1(1)+xlin1(ipar(ii))*drprod1
        dsum1(2)=dsum1(2)+xlin1(ipar(ii))*drprod2
        dsum1(3)=dsum1(3)+xlin1(ipar(ii))*drprod3
        dsum2(1)=dsum2(1)+xlin2(ipar(ii))*drprod1
        dsum2(2)=dsum2(2)+xlin2(ipar(ii))*drprod2
        dsum2(3)=dsum2(3)+xlin2(ipar(ii))*drprod3
        dsum3(1)=dsum3(1)+xlin3(ipar(ii))*drprod1
        dsum3(2)=dsum3(2)+xlin3(ipar(ii))*drprod2
        dsum3(3)=dsum3(3)+xlin3(ipar(ii))*drprod3
      enddo     


      v3c=sum1+sum2*bb2+sum3*bb3
      eee=v+damper*v3c

      do i = 1,3
         dv(i)=dv(i) +ddamper(i)*v3c +damper*(sum2*dbb2(i)+sum3*dbb3(i)
     &        +dsum1(i)+dsum2(i)*bb2+dsum3(i)*bb3)
      enddo


      return
      end


      double precision function bodc2body(r1x,dv)
      implicit none
      double precision :: xlp(14), beta, re, v, r1x, dv, t
      integer :: j, indx
!
!     distances are in Bohr, energies are in cm-1
!     The global minimum is at: R=2.3178250879779552  E= -17.5988842162493135 cm-1
!
      beta= 3.4059169d0           
      xlp( 1) = -4.942622126d0 
      xlp( 2) = -4.742197464d1 
      xlp( 3) = -1.115075510d2      
      xlp( 4) = -1.419033104d2  
      xlp( 5) = -1.165969889d2 
      xlp( 6) = -6.279822284d1 
      xlp( 7) = -2.626186059d1      
      xlp( 8) = -1.795750569d1        
      xlp( 9) = -6.826792412d0       
      xlp(10) =  7.563692060d0       
      xlp(11) =  2.630227612d0        
      xlp(12) = -2.595885106d0       
      xlp(13) =  5.647782847d-1       
      xlp(14) = -4.042171086d-2        
      re=1.401d0

      v=xlp(14)
      dv=-beta*xlp(14)
      indx=14
      do j=2,14
         indx=indx-1
         v=v*(r1x-re)+xlp(indx)
         dv=dv*(r1x-re)+indx*xlp(indx+1)-beta*xlp(indx)
      enddo
      t=exp(-beta*(r1x-re))
      v=v*t
      dv=dv*t
      bodc2body=v

      end

      subroutine indexa2(mbig,mtop)            
      implicit none
      integer :: mbig,mtop,n2,iind,jind,kind,ipar,maxi,maxi2
      integer :: l,l2,i,j,k,isum
      parameter (n2=5000)
c     calculate the fitting indices for AB2 symmetry
      common /indy2/iind(n2),jind(n2),kind(n2),ipar(n2),maxi,maxi2
      l=0
      l2=0
      do 10 i=0,min(mtop,mbig)
       do 20 j=0,min(mtop,mbig)
        do 30 k=j,min(mtop,mbig)
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
            if(j.ne.k)then
               l=l+1
               ipar(l)=l2
               iind(l)=i
               jind(l)=k
               kind(l)=j
            endif
           endif
          endif
30      continue
20     continue
10    continue
      maxi=l2
      maxi2=l
      return
      end


      SUBROUTINE prepot_leps
!
!       PREPOT must be called once before any calls to LEPS
!       potential parameters are read in PREPOT
!
      implicit none
      integer          :: i
      double precision :: d(3),re(3),beta(3),sato(3),rt3
      common /lepscm/  d,re,beta,sato,rt3
      save

      d(:)=109.4583d0          
      re(:)=.7412871d0           
      beta(:)=1.973536d0             
      sato(:)=.132d0
      RT3 = SQRT(3d0)

! Convert to atomic units

      Do I = 1,3
         D(I)=D(I)/627.509552d0
         Re(I) = Re(I)/0.52917725d0
         Beta(I) = Beta(I)*0.52917725d0
      End do

!     WRITE (6,600) (D(I)*627.509552,I=1,3),Re,Beta,Sato
      RETURN
600   FORMAT (//,10x,' LEPS Potential Energy Parameters',
     & //,10x,' Bond',27x,'AB',8x,'BC',8x,'AC',/,10x,
     & ' Diss. Energies (kcal/mol)',T40,3F10.5,/,10x,' Req (au)',
     & T40,3F10.5,/,10x,' Morse Beta (au)',T40,3F10.5,/,10x,
     & ' Sato Parameters',T40,3F10.5)
      END


      subroutine leps (r,v1,v2,xma,xmb,xmc,gtot,dgtot)
!
!   r1 is the distance between A and B
!   r2 is the distance between B and C
!   r3 is the distance between A and C
!   v1 is the lower energy root
!   v2 is the upper energy root
!   xma, xmb, xmc are the atomic masses
!   ga, gb, gc are the contributions to the nonadiabatic corrections
!      from each atom
!
      implicit none
      integer          :: i,j
      double precision :: r(3),v1,v2,xma,xmb,xmc,gtot,dgtot(3)
      double precision :: ex,vs(3),vt(3),dvs(3),dvt(3)
      double precision :: h11,h12,h22,hsum,hdif,dhdif(3),rad,delh
      double precision :: f12(3),cos1,cos2,cos3,ga,gb,gc
      double precision :: d2vs(3),d2vt(3),dh11(3),dh12(3),dh22(3)
      double precision :: d2h12(3),d2hdif(3),drad(3),df12(3,3)
      double precision :: dcos1(3),dcos2(3),dcos3(3)

      double precision :: d(3),re(3),beta(3),sato(3),rt3
      common /lepscm/  d,re,beta,sato,rt3
!
! Compute single and triplet diatomic potentials
      do i = 1,3
         ex = exp(-beta(i)*(r(i)-re(i)))
         vs(i) = d(i)*(ex - 2.0d0)*ex
         dvs(i) = -2.0d0*d(i)*beta(i)*(ex-1.0d0)*ex
         d2vs(i) = 2.0d0*d(i)*beta(i)*beta(i)*(2.d0*ex-1.0d0)*ex
         vt(i) = 0.5d0*d(i)*(ex + 2.0)*ex*(1.0-sato(i))/(1.0+sato(i))
         dvt(i) = -d(i)*beta(i)*(ex+1.0)*ex*(1.0-sato(i))/(1.0+sato(i))
         d2vt(i) = d(i)*beta(i)*beta(i)*(2.0d0*ex+1.0)*ex
     &        *(1.0-sato(i))/(1.0+sato(i))
      end do
!
! Construct 2x2 DIM Hamiltonian
      h11 = d(2) + vs(2) + 0.25d0*(vs(1)+vs(3)) 
     &   + 0.75d0*(vt(1)+vt(3))
      h12 = 0.25d0*rt3*(vs(1) - vs(3) - vt(1) + vt(3))
      h22 = d(2) + vt(2) + 0.75*(vs(1)+vs(3)) 
     &   + 0.25d0*(vt(1)+vt(3))

!
! Calculate eigenvalues (V1 and V2) and quantities needed for matrix elements
!   of d/dRi
!
      hsum = h11 + h22
      hdif = h22 - h11
      dhdif(1) = 0.5d0*(dvs(1)-dvt(1))
      dhdif(2) = -(dvs(2) - dvt(2))
      dhdif(3) = 0.5d0*(dvs(3)-dvt(3))
      rad = 0.25d0*hdif*hdif + h12*h12
      delh = sqrt(rad)
      v1 = 0.5d0*hsum - delh
      v2 = 0.5d0*hsum + delh

! derivatives
!
      dh11(1) = 0.25d0*dvs(1)+0.75d0*dvt(1)
      dh11(2) = dvs(2)
      dh11(3) = 0.25d0*dvs(3)+0.75d0*dvt(3)

      dh12(1) = 0.25d0*rt3*(dvs(1)-dvt(1))
      dh12(2) = 0.0d0
      dh12(3) =-0.25d0*rt3*(dvs(3)-dvt(3))

      dh22(1) = 0.75d0*dvs(1)+0.25d0*dvt(1)
      dh22(2) = dvt(2)
      dh22(3) = 0.75d0*dvs(3)+0.25d0*dvt(3)

      d2h12(1) = 0.25d0*rt3*(d2vs(1)-d2vt(1))
      d2h12(2) = 0.0d0
      d2h12(3) =-0.25d0*rt3*(d2vs(3)-d2vt(3))

      d2hdif(1) = 0.5d0*(d2vs(1) - d2vt(1))
      d2hdif(2) = -(d2vs(2) - d2vt(2))
      d2hdif(3) = 0.5d0*(d2vs(3) - d2vt(3))

      drad(1) = 0.5d0*hdif*dhdif(1) +2.0d0*h12*dh12(1)
      drad(2) = 0.5d0*hdif*dhdif(2) +2.0d0*h12*dh12(2)
      drad(3) = 0.5d0*hdif*dhdif(3) +2.0d0*h12*dh12(3)


! Calculate matrix elements of d/dRi; note that at conical intersection f12
!   diverges, but is set to zero in the code.
      do i = 1,3
         f12(i) = 0.25d0*(h12*dhdif(i) - hdif*dh12(i))
         if (rad.ne.0.0d0) then
            f12(i) = f12(i)/rad
            do j = 1,3
               df12(i,j) = (0.25d0*(dh12(j)*dhdif(i)-dhdif(j)*dh12(i))
     &          -f12(i)*drad(j))/rad
            enddo
            df12(i,i) = df12(i,i) 
     &       + 0.25d0*(h12*d2hdif(i)-hdif*d2h12(i))/rad
         endif
      enddo


!
! Calculate diagonal matrix element of kinetic energy operator
!
      cos1 = (r(2)*r(2) + r(3)*r(3) - r(1)*r(1))/(2.0*r(2)*r(3))
      cos2 = (r(1)*r(1) + r(3)*r(3) - r(2)*r(2))/(2.0*r(1)*r(3))
      cos3 = (r(1)*r(1) + r(2)*r(2) - r(3)*r(3))/(2.0*r(1)*r(2))
      ga=(f12(1)*(f12(1) + 2.0*Cos2*f12(3)) + f12(3)*f12(3))/(2.d0*xma)
      gb=(f12(1)*(f12(1) + 2.0*Cos3*f12(2)) + f12(2)*f12(2))/(2.d0*xmb)
      gc=(f12(2)*(f12(2) + 2.0*Cos1*f12(3)) + f12(3)*f12(3))/(2.d0*xmc)
      gtot=ga+gb+gc

      dcos1(1) = -r(1)/(r(2)*r(3))
      dcos1(2) = 1.0d0/r(3) - cos1/r(2)
      dcos1(3) = 1.0d0/r(2) - cos1/r(3)
      dcos2(1) = 1.0d0/r(3) - cos2/r(1)
      dcos2(2) = -r(2)/(r(1)*r(3))
      dcos2(3) = 1.0d0/r(1) - cos2/r(3)
      dcos3(1) = 1.0d0/r(2) - cos3/r(1)
      dcos3(2) = 1.0d0/r(1) - cos3/r(2)
      dcos3(3) = -r(3)/(r(1)*r(2))


      do i = 1,3
! contribution from ga
         dgtot(i) = 
     &       (f12(1)*(df12(1,i) + dcos2(i)*f12(3) + cos2*df12(3,i))
     &        + cos2*df12(1,i)*f12(3) + f12(3)*df12(3,i))/xma
! contribution from ga
         dgtot(i) = dgtot(i) +
     &       (f12(1)*(df12(1,i) + dcos3(i)*f12(2) + cos3*df12(2,i))
     &        + cos3*df12(1,i)*f12(2) + f12(2)*df12(2,i))/xmb
! contribution from ga
         dgtot(i) = dgtot(i) +
     &       (f12(2)*(df12(2,i) + dcos1(i)*f12(3) + cos1*df12(3,i))
     &        + cos1*df12(2,i)*f12(3) + f12(3)*df12(3,i))/xmc
      enddo

      return



      return
      end


      module ccidata
      save
      double precision ::  eshift= 0.1744759989649372d0
      double precision ::  beta=0.9246250191375499d0    
      double precision ::  rnought=2.d0                         
      double precision, dimension(89) ::  xlin1=(/-0.1682626669809615d0,        
     &     1.697952225005269d0,      3.959367441888326d0,   
     &    -7.159854944676235d0,     -2.462488887191634d0,  
     &     7.998434264390882d0,     -3.665712084901418d0, 
     &    -0.2155463656483861d0,     0.5620418518154215d0, 
     &    -0.1369812453279634d0,     1.0063614579308723d-2,
     &    41.07358514471853d0,      67.02217912830267d0, 
     &  -100.0014758720550d0,        2.481875721304641d0, 
     &    32.46244740624601d0,     -15.32460432587821d0, 
     &     3.814631583515083d0,     -0.3155985331636701d0, 
     &    -4.2198516063816239d-3, -157.7203098415152d0, 
     &    17.11141015275416d0,     -61.48638359701442d0, 
     &     4.244841908917134d0,      5.300393723900312d0, 
     &    -0.6170558358973661d0,    -6.2143674437616274d-3,
     &     5.785766830955630d0,     -0.1133784523542150d0, 
     &    -0.2539145873400511d0,    -4.3406559735048131d-2,
     &     2.4548851428702549d-2,   -1.092841685240530d0, 
     &     0.8131178021738225d0,    -7.7933784638615652d-2,
     &    -0.1484455564689179d0,    -7.351078313940695d0, 
     &   128.7480567271992d0,      199.6058760392980d0, 
     &   -50.18444648711906d0,    -178.2184990028168d0, 
     &   129.6832433574917d0,      -37.57036725437851d0, 
     &     5.747838042162251d0,      0.1434298979682676d0, 
     &    -8.3110343245055174d-2,  428.7635440340811d0, 
     &   495.9564113906064d0,    -1742.123837168904d0, 
     &  -259.0564024322437d0,      250.2352551214438d0, 
     &   -44.15232947384926d0,       0.9963955068997190d0, 
     &     0.3136096547683901d0, -1032.921631265462d0, 
     &   199.9708129167335d0,     -135.4683186684037d0, 
     &   -14.50534672808297d0,       5.650612567917722d0, 
     &    -0.1473535040913008d0,     1.351710421967908d0, 
     &    29.60972395839821d0,      -7.224575522867537d0, 
     &     0.3582134126835411d0,     0.2589923362572861d0, 
     &     0.1819402887127863d0,  1312.633253770123d0, 
     &  5189.281073844522d0,     -2727.020597875299d0, 
     &  -122.2193440638810d0,      152.4859537357477d0, 
     &    -5.562127541501408d0,     -1.549418235204051d0, 
     &  2647.710764337011d0,       502.1991025406292d0, 
     &  -152.3069488960497d0,       -5.740041181528469d0, 
     &     1.362209926084833d0,     17.96504193139269d0, 
     &    10.13930464010770d0,      -2.287463377365149d0, 
     &     1.932492687974671d0,  -3337.181754475553d0, 
     &   182.9671187483376d0,       20.29606065582956d0, 
     &     1.192500349837268d0,    -22.69079850990068d0, 
     &    -0.5967747962499963d0,     1.091111370927579d0/)
      double precision, dimension(89) :: xlin2=(/3.2289532152874986d-2,   
     &    -8.6686506736620789d-2,   -0.3266024449171150d0,   
     &    -0.7972798579742405d0,     2.812854025206053d0,   
     &    -2.374871328104482d0,      0.6089164699935086d0,   
     &     0.1635658322754857d0,    -0.1263080841790984d0,   
     &     2.5221027611065153d-2,   -1.6677760633386995d-3,
     &    -0.3372069966080046d0,    -8.377537389356892d0, 
     &     3.694668699246826d0,     11.36868357221171d0, 
     &   -11.88889428297186d0,       4.270226804928180d0, 
     &    -0.4526863369836274d0,    -5.2518975083660345d-2,
     &     9.0223742463199624d-3,   38.03127084604298d0, 
     &   -12.82065963469141d0,      -5.663705994784566d0, 
     &     7.300080908068341d0,     -0.5157232521089076d0, 
     &    -0.2200954472805637d0,     2.0686337667851093d-2,
     &    14.19614073502447d0,       1.726255538125508d0, 
     &    -2.526391344985463d0,      6.0598667264186397d-2,
     &     3.8560006141770319d-2,    2.306725001955670d0, 
     &    -0.3197660708712612d0,    -1.9982819587147187d-2,
     &     0.1036787525350017d0,    -2.030231243513758d0, 
     &    25.16474735933574d0,     -14.13898387128746d0, 
     &   -62.01219930765056d0,      75.56471352787402d0, 
     &   -30.25154027581496d0,       4.257663341071418d0, 
     &     0.6142213624561507d0,    -0.2458345005137628d0, 
     &     1.9135101060142316d-2,  168.2482649287792d0, 
     &  -270.3587962365308d0,       72.46467150857366d0, 
     &    32.12278113625301d0,      11.89435258275814d0, 
     &    -5.521268211450904d0,      0.8934653854937415d0, 
     &    -7.6805320500971247d-2, -422.2672942154591d0, 
     &   -35.92489169679082d0,      14.74665733902546d0, 
     &     6.020706922914545d0,      2.4538599994799737d-3,
     &    -0.1370624616096560d0,    44.47654636843973d0,
     &   -11.18166310447657d0,       3.6110526441458540d-2,
     &     0.1539535582849851d0,     2.065981296572923d0, 
     &    -0.1624439939014451d0,  -658.3723053145102d0, 
     &  -883.4947620258348d0,     -263.8479489759550d0, 
     &    44.68359598895833d0,      29.88970883430300d0, 
     &    -5.663152586850906d0,      0.2413124485342995d0, 
     &   656.9813540494517d0,      -91.56710742448996d0, 
     &   -30.38175110815412d0,       2.066013906655882d0, 
     &     0.5314632792073760d0,    54.99422880821747d0, 
     &    -0.7887249891070589d0,    -0.9164068632403487d0, 
     &     0.6576346543435602d0,   292.2665706493461d0, 
     &   -32.04152511166628d0,      -2.104943871500433d0, 
     &   -0.7251129002870176d0,      4.709122797258667d0,   
     &    1.119971732384046d0,      -4.066519840824071d0/)
      double precision, dimension(71) ::  xlin3=(/6.3998636332917941d-2,    
     &    -0.5826083150212633d0,    -3.294968443823759d0, 
     &     5.751817019955082d0,     -1.486176533347650d0, 
     &    -0.6845411439568367d0,    -0.9083627566591157d0, 
     &     1.305678536333854d0,     -0.4821780067143284d0, 
     &     5.5182580699196751d-2,  -22.57213192201081d0, 
     &   -61.37436407473290d0,      52.64511100643983d0, 
     &    29.99521825023724d0,     -38.31341936937510d0, 
     &    14.86681566870566d0,      -3.460352146284503d0, 
     &     0.1343449644795501d0,    53.21459427503272d0, 
     &    64.60699757442232d0,      30.89466488573955d0, 
     &     8.412732049211240d0,     -4.274560634850636d0, 
     &    -0.1788031871862021d0,   -46.47154313325355d0, 
     &    24.88995914944731d0,      -2.757433396373881d0, 
     &    -0.3635530438220282d0,     1.049428472748707d0, 
     &    -0.9909460260393127d0,     7.623381850055868d0, 
     &   -75.76387696233419d0,    -162.8879929453359d0, 
     &   -29.98317728892977d0,     168.7972614497797d0, 
     &  -102.0864921135087d0,       31.38363768116054d0, 
     &    -6.084136544334489d0,     -2.9400115678577254d-2,
     &  -417.4808742201730d0,     -601.3246381855070d0,
     &  1111.017825017578d0,       749.0346863822200d0, 
     &  -254.9798923160399d0,       37.02896288679148d0, 
     &    -0.7010161055223685d0,   928.6100310886904d0, 
     &   513.3554862860957d0,      138.4248255502364d0, 
     &    -0.9347938923678116d0,     3.867048618649862d0, 
     &  -109.6794304656507d0,       10.30110710850976d0, 
     &     0.5118243051520034d0,     0.2069543381912151d0, 
     &  -796.4010440692600d0,    -4693.456970734772d0, 
     &  2597.246850944603d0,       630.4845628327919d0, 
     &  -221.0729177223342d0,        9.839396652075672d0, 
     & -5227.790162633692d0,      -529.4879950032479d0, 
     &    30.53701394166368d0,     -14.22864486129778d0, 
     &   -35.14731691232372d0,       2.140914121987492d0,  
     &  2250.367423678786d0,       138.5994085512685d0,  
     &    70.30416282155014d0,     -62.74844809925895d0/) 
      double precision, dimension(71) :: xlin4=(/-1.2612782515309912d-2,
     &     1.5724446750473936d-2,     0.2653887372495336d0, 
     &     0.1216019078797281d0,     -0.6102844930670878d0, 
     &     4.6420446436398755d-3,     0.5563651571664677d0, 
     &    -0.3890198684837317d0,      0.1052569138798340d0, 
     &    -1.0080347996402828d-2,     0.2409741384095758d0, 
     &     5.745106082351191d0,       1.047295525278513d0, 
     &   -12.03507114814484d0,        9.035186991905373d0, 
     &    -2.267198360038944d0,       2.0836099694222385d-2,
     &     4.7445938282087222d-2,   -18.26537244396756d0, 
     &    -4.809068965727743d0,       5.658647992475390d0, 
     &    -1.284093508478602d0,      -2.164733987255563d0, 
     &     0.3630073042947094d0,      3.760490137844927d0, 
     &    -5.927286374087972d0,       0.1278839922444694d0, 
     &     0.5381095462111573d0,     -2.179370013756247d0, 
     &     0.4127566904703107d0,      0.7168269469196163d0, 
     &    -9.502943939736065d0,      -9.597638311349975d0, 
     &    65.60860814240448d0,      -56.04003254284753d0, 
     &    15.38423202570607d0,       -0.3167914651468724d0, 
     &    -0.7643266917502949d0,      0.1110705439245336d0, 
     &  -146.5638710976931d0,       151.8079548350792d0, 
     &    15.86486188082157d0,      -65.56160562330304d0, 
     &    -1.126379669058077d0,      -1.759058937005427d0, 
     &    -3.1582037336166588d-4,   410.9371301349821d0, 
     &   142.2842073695382d0,       -15.14809406881419d0, 
     &    -9.858662824266736d0,       0.2622049408092991d0, 
     &   -72.67984057864001d0,       11.28350883179339d0, 
     &    -0.5323962789000201d0,     -2.447540583751350d0, 
     &   148.9575982135905d0,      1249.939675860206d0, 
     &   268.6478693121586d0,        33.67689052007242d0, 
     &   -35.65384689269180d0,        3.464021415311834d0, 
     &  -103.4190791219546d0,        76.76961997409407d0, 
     &    11.59182442022548d0,       -0.5254501536114633d0, 
     &   -32.07157614897808d0,        1.102243506822152d0, 
     &  -952.9757997893194d0,        52.11461967121677d0, 
     &     6.929202177070909d0,      -6.322557284152793d0/)
      end module


      subroutine prepotcbs
      use ccidata
      implicit none
      integer          :: n2,iind,jind,kind,ipar,maxi,maxi2
      integer          :: maxi2a,maxi2b,ii,i,it
      save
      parameter (n2=5000)
      common /indy/iind(n2),jind(n2),kind(n2),ipar(n2),maxi,maxi2

      double precision :: epslon,epscem,ald,betad
      double precision :: r1x,r2x,r3x,eee,dv(3)
      double precision :: etrip1,etrip2,etrip3,detrip1,detrip2,detrip3
      double precision :: es1,es2,es3,des1,des2,des3
      double precision :: q1,q2,q3,xj1,xj2,xj3,xj,vlondon
      double precision :: dq1,dq2,dq3,dxj1,dxj2,dxj3,dvlon(3)
      double precision :: damper,ddamp(3),rdamp1,rdamp2,rdamp3
      double precision :: bb,dbb(3),warp,dwarp(3),rho1,rho2,rho3
      double precision :: rho1t(0:12),rho2t(0:12),rho3t(0:12)
      double precision :: sum1,sum2,sum3,sum4
      double precision :: dsum1(3),dsum2(3),dsum3(3),dsum4(3)
      double precision :: rprod,t1,v3c


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

      entry potcbs(r1x,r2x,r3x,eee,dv)

      call trip(r1x,etrip1,detrip1)
      call trip(r2x,etrip2,detrip2)
      call trip(r3x,etrip3,detrip3)
      call singlet(r1x,es1,des1)
      call singlet(r2x,es2,des2)
      call singlet(r3x,es3,des3)

C      write (99,*) ' detrip,des=',detrip1,detrip2,detrip3,des1,des2,des3

      q1=0.5d0*(es1+etrip1)
      q2=0.5d0*(es2+etrip2)
      q3=0.5d0*(es3+etrip3)
      xj1=0.5d0*(es1-etrip1)
      xj2=0.5d0*(es2-etrip2)
      xj3=0.5d0*(es3-etrip3)
      xj=sqrt(epslon+0.5d0*((xj3-xj1)**2+(xj2-xj1)**2+(xj3-xj2)**2))    
      vlondon=q1+q2+q3-xj
      dq1=0.5d0*(des1+detrip1)
      dq2=0.5d0*(des2+detrip2)
      dq3=0.5d0*(des3+detrip3)
      dxj1=0.5d0*(des1-detrip1)
      dxj2=0.5d0*(des2-detrip2)
      dxj3=0.5d0*(des3-detrip3)
      dvlon(1)=dq1-0.5d0*dxj1*(2.0d0*xj1-xj2-xj3)/xj
      dvlon(2)=dq2-0.5d0*dxj2*(2.0d0*xj2-xj1-xj3)/xj
      dvlon(3)=dq3-0.5d0*dxj3*(2.0d0*xj3-xj1-xj2)/xj

      if((r1x.le.betad).or.(r2x.le.betad).or.(r3x.le.betad))then
        damper=0.d0
        ddamp(1)=0.d0
        ddamp(2)=0.d0
        ddamp(3)=0.d0
      else
        rdamp1=ald/(betad-r1x)
        rdamp2=ald/(betad-r2x)
        rdamp3=ald/(betad-r3x)
        damper=exp(rdamp1+rdamp2+rdamp3)
        ddamp(1)=rdamp1/(betad-r1x)
        ddamp(2)=rdamp2/(betad-r2x)
        ddamp(3)=rdamp3/(betad-r3x)
      endif

      bb=sqrt(epscem+(r1x-r3x)**2+(r1x-r2x)**2+(r2x-r3x)**2)
      dbb(1) = (2.d0*r1x-r2x-r3x)/bb
      dbb(2) = (2.d0*r2x-r1x-r3x)/bb
      dbb(3) = (2.d0*r3x-r1x-r2x)/bb

      warp=1.d0/(1.d0/r1x+1.d0/r2x+1.d0/r3x)
      dwarp(1)=(warp/r1x)**2
      dwarp(2)=(warp/r2x)**2
      dwarp(3)=(warp/r3x)**2

      rho1=exp(-beta*(r1x-rnought))
      rho2=exp(-beta*(r2x-rnought))
      rho3=exp(-beta*(r3x-rnought))
      rho1t(0)=1.d0
      rho2t(0)=1.d0
      rho3t(0)=1.d0
      do it=1,12
        rho1t(it)=rho1t(it-1)*rho1
        rho2t(it)=rho2t(it-1)*rho2
        rho3t(it)=rho3t(it-1)*rho3
      enddo    

      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      sum4=0.d0
      do i = 1,3
         dsum1(i)=0.d0
         dsum2(i)=0.d0
         dsum3(i)=0.d0
         dsum4(i)=0.d0
      enddo

      do ii=1,maxi2a
        rprod=rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
        t1=xlin1(ipar(ii))*rprod
        sum1=sum1+t1
        dsum1(1)=dsum1(1)-beta*kind(ii)*t1
        dsum1(2)=dsum1(2)-beta*jind(ii)*t1
        dsum1(3)=dsum1(3)-beta*iind(ii)*t1
        t1=xlin2(ipar(ii))*rprod
        sum2=sum2+t1
        dsum2(1)=dsum2(1)-beta*kind(ii)*t1
        dsum2(2)=dsum2(2)-beta*jind(ii)*t1
        dsum2(3)=dsum2(3)-beta*iind(ii)*t1
      enddo     

      do ii=maxi2a+1,maxi2a+maxi2b
        rprod=rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))
        t1=xlin3(ipar(ii))*rprod
        sum3=sum3+t1
        dsum3(1)=dsum3(1)-beta*kind(ii)*t1
        dsum3(2)=dsum3(2)-beta*jind(ii)*t1
        dsum3(3)=dsum3(3)-beta*iind(ii)*t1
        t1=xlin4(ipar(ii))*rprod
        sum4=sum4+t1
        dsum4(1)=dsum4(1)-beta*kind(ii)*t1
        dsum4(2)=dsum4(2)-beta*jind(ii)*t1
        dsum4(3)=dsum4(3)-beta*iind(ii)*t1
      enddo     

      v3c=damper*(sum1+sum2*bb+warp*(sum3+sum4*bb))
      eee=vlondon+v3c+eshift

      do i = 1,3
         dv(i)=dvlon(i)
     *        +ddamp(i)*v3c
     *        +damper*(dsum1(i)+dbb(i)*(sum2+warp*sum4)+bb*dsum2(i)
     *        +dwarp(i)*(sum3+bb*sum4)+warp*(dsum3(i)+bb*dsum4(i)))
      enddo

      return
      end
      

      subroutine indexa3(mbig,mtop)            
      implicit none
      integer :: mbig,mtop
      integer :: n2,iind,jind,kind,ipar,maxi,maxi2
      integer :: l,l2,i,j,k,isum
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


      subroutine trip(r,v,dv)
c     H2 triplet curve for FCI/CBS       
      implicit none
      integer          :: j
      double precision :: r,v,dv,t0,t1,t2,vlr,dvlr,vsr,dvsr
      double precision :: damp,xd,xd2,prefac
      double precision :: c6=-6.499027d0
      double precision :: c8=-124.3991d0
      double precision :: c10=-3285.828d0
      double precision :: beta=0.2d0
      double precision :: re=1.401d0

      double precision :: alpha=4.342902497495857d0    
      double precision,dimension(17) :: xlp=(/2.190254292077946d-1,
     &           6.903970886899540d-1, 1.163805164813647d+0,
     &           1.337274224278977d+0, 1.176131501286896d+0,
     &           9.093025316067311d-1, 5.908490740317524d-1,
     &           1.467271248985162d-1, 6.218667627370255d-2,
     &           1.105747186790804d-1, -7.112594221275334d-2,
     &          -8.866929977503379d-3, 3.823523348718375d-2,
     &          -2.055915301253834d-2, 5.419433399532862d-3,
     &          -7.263483303294653d-4, 4.154351972583585d-5/)

      t0=r-re
      prefac=exp(-alpha*t0)
      if((r-re).ne.0.d0)then
         vsr = 0.0
         dvsr = 0.0
         t1=1.0d0
         t2=0.0d0
         do j=1,17
C            v=v+xlp(j)*(r-re)**(j-1)*prefac
            vsr=vsr+xlp(j)*t1*prefac
              dvsr=dvsr+(j-1)*xlp(j)*t2*prefac
            t2=t1
            t1=t1*t0
         enddo   
      else
C         v=v+xlp(1)*prefac
         vsr=xlp(1)*prefac
         dvsr=xlp(2)*prefac
      endif
      dvsr=dvsr-alpha*vsr

C      damp=1.d0-exp(-beta*r**2)
      t1 = r*r
      t2 = exp(-beta*t1)
      damp=1.d0-t2
      xd=damp/r  
      xd2=xd*xd
C      v=c6*xd**6+c8*xd**8+c10*xd**10
      vlr=(c6+xd2*(c8+xd2*c10))*xd2**3
      dvlr=-(6.0d0*c6+xd2*(8.0d0*c8+xd2*10.0d0*c10))*
     *   ((1/r)-2.0d0*beta*t2/xd)*xd2**3

      v=vlr+vsr
      dv=dvlr+dvsr

      return
      end

      subroutine singlet(r,v,dv)
c     H2 singlet potential for FCI/CBS
      implicit none
      integer          :: j
      double precision :: r,v,dv,t0,t1,t2,vsr,damp,xd,xd2,prefac
      double precision :: c6=-6.499027d0
      double precision :: c8=-124.3991d0
      double precision :: c10=-3285.828d0
      double precision :: beta=0.2d0
      double precision :: re=1.401d0
      double precision :: alpha=3.980917850296971d0    
      double precision, dimension(17) :: xlp= (/ -1.709663318869429d-1,
     &          -6.675286769482015d-1, -1.101876055072129d+0, 
     &          -1.106460095658282d+0, -7.414724284525231d-1,
     &          -3.487882731923523d-1, -1.276736255756598d-1,
     &          -5.875965709151867d-2, -4.030128017933840d-2,
     &          -2.038653237221007d-2, -7.639198558462706d-4,
     &           2.912954920483885d-3, -2.628273116815280d-4,
     &          -4.622088855684211d-4, 1.278468948126147d-4,
     &          -1.157434070240206d-5, -2.609691840882097d-12/)
C      damp=1.d0-exp(-beta*r**2)
      t1 = r*r
      t2 = exp(-beta*t1)
      damp=1.d0-t2
      xd=damp/r  
      xd2=xd*xd
C      v=c6*xd**6+c8*xd**8+c10*xd**10
      v=(c6+xd2*(c8+xd2*c10))*xd2**3
      dv=-(6.0d0*c6+xd2*(8.0d0*c8+xd2*10.0d0*c10))*
     *   ((1/r)-2.0d0*beta*t2/xd)*xd2**3

      t0=r-re
      prefac=exp(-alpha*t0)
      vsr = 0.0
      if((r-re).ne.0.d0)then
         t1=1.0d0
         t2=0.0d0
         do j=1,17
C            v=v+xlp(j)*(r-re)**(j-1)*prefac
            vsr=vsr+xlp(j)*t1*prefac
            dv=dv+(j-1)*xlp(j)*t2*prefac
            t2=t1
            t1=t1*t0
         enddo   
      else
C         v=v+xlp(1)*prefac
         vsr=xlp(1)*prefac
         dv=dv+xlp(2)*prefac
      endif
      v=v+vsr
      dv=dv-alpha*vsr

      return
      end


      subroutine golden(ax,bx,cx,tol,xmin,fval)   
       implicit none
       external fcn

       double precision, parameter :: R=0.61803399d0, C = 1d0-R
       double precision :: f1, f2, x0, x1, x2, x3, tol
       double precision :: fval, xmin, bx, ax, cx, fcn
       double precision :: fa, fb, fc
       integer            :: icount
       integer, parameter :: maxcount=5000

       fa=fcn(ax)
       fb=fcn(bx)
       fc=fcn(cx)
       if((fb.ge.fa).or.(fb.ge.fc)) then
         write(6,*)'bracketing error in golden'
         write(6,*)ax,fa
         write(6,*)bx,fb
         write(6,*)cx,fc
         stop 
       endif

       x0=ax
       x3=cx

       if(abs(cx-bx) .gt. abs(bx-ax))then
         x1=bx
         f1=fb
         x2=bx+c*(cx-bx)
         f2=fcn(x2)
       else
         x2=bx
         f2=fb
         x1=bx-c*(bx-ax)
         f1=fcn(x1)
       endif
         
       icount = 0

       do
         if(abs(x3-x0) .le. tol*(abs(x1)+abs(x2))) exit

         icount=icount+1
         if(icount.ge.maxcount)then
           write(6,*)'too many evaluations in golden'
           stop
         endif

         if(f2.lt.f1)then
           x0=x1
           x1=x2
           x2=R*x1+C*x3
           f1=f2
           f2=fcn(x2)
         else
           x3=x2
           x2=x1
           x1=R*x2+C*x0
           f2=f1
           f1=fcn(x1)
         endif
       enddo

       if(f1.lt.f2)then
         fval=f1
         xmin=x1
        else
         fval=f2
         xmin=x2
       endif

       end


       double precision function fcn(x)
       use zmasses
       implicit none
       double precision :: x, e, bodc2body, v, dv, towave, xmh, massfact
!      Evaluates the Born-Huang energy of an isomer of H2 with masses of zmass1 and zmass2
       xmh=1836.15264d0   ! H nuclear mass (after substracting 1 for electron mass)
       massfact=0.5d0*xmh*(zmass1+zmass2)/(zmass1*zmass2)
       towave=219474.7d0 
       call singlet(x,v,dv)
       fcn=massfact*bodc2body(x,dv)/towave+v
       end    


      subroutine pes(x,igrad,p,g,d)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, parameter :: natoms=3
      integer, intent(in) :: igrad

      double precision, intent(in) :: x(natoms,3)
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: r(3), dvp(3), v
      double precision :: tx(9), r2(3), drdx(3,9)
      double precision :: gx(nstates,9)
      integer :: iatom, idir, i, j
      logical, save :: first_time_data=.true.


      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        tx(j)=x(iatom,idir)
      enddo
      enddo
      ! input cartesian is HHH
      r(1)=sqrt((x(1,1)-x(2,1))**2+(x(1,2)-x(2,2))**2
     *          +(x(1,3)-x(2,3))**2)/0.529177211
      r(2)=sqrt((x(2,1)-x(3,1))**2+(x(2,2)-x(3,2))**2
     *          +(x(2,3)-x(3,3))**2)/0.529177211
      r(3)=sqrt((x(1,1)-x(3,1))**2+(x(1,2)-x(3,2))**2
     *          +(x(1,3)-x(3,3))**2)/0.529177211

      r2=r*0.529177211

      call pot(r,v,dvp)

      p=v*27.211386
      dvp=dvp*51.422067

      call evdrdx(tx, r2, drdx)
      gx=0.d0
      do i=1,9
      do j=1,3
        gx(:,i)=gx(:,i)+dvp(j)*drdx(j,i)
      enddo
      enddo
  
      do iatom=1, natoms
      do idir=1,3
        j=(iatom-1)*3+idir
        g(:,iatom,idir)=gx(:,j)
      enddo
      enddo

      if (igrad==2) then
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
