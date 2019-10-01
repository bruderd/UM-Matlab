function [Y,Xf,Af] = nnet_2019_09_30_17_11(X,~,~)
%NNET_2019_09_30_17_11 neural network simulation function.
%
% Auto-generated by MATLAB, 30-Sep-2019 17:11:34.
% 
% [Y] = nnet_2019_09_30_17_11(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 120xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.210683480440972;-0.278877667618903;-1;-1;-1;-1;-1;-0.310441377421948;-0.476120590899149;-0.402413531056711;-0.101630002142874;-0.278792025847369;-2.63775279396583e-10;-0.212473807754962;-0.13122529556891;-1.45591740486478e-10;-1.98181946547971e-10;-0.162577409599759;-0.296454749608797;-0.155397845598898;-0.291220882895356;-0.359558048455584;-0.199419005175685;-0.379496173330994;-0.165360437120805;-0.370789484231199;-0.328540683230395;-0.335711531722636;-0.222391304958488;-0.451892016961284;-0.195826019319836;-0.379338598051715;-0.276333664168532;-0.343725926132977;-0.337091538715417;-0.400215464156118;-0.479220902824089;-0.401173798197997;-0.42164397673007;-0.20540007185747;-0.297766831420243;-0.234732975533792;-0.0888038744031764;-0.0886111614883502;-7.46339740009816e-10;-0.250599156461505;-0.23343434008084;-0.197491161533254;-0.0794531117304958;-0.069401882923825;-9.68153446301301e-10;-0.119788378056081;-0.0930286834757229;-1.29962085537727e-09;-7.03182307520978e-10;-0.405455152751332;-0.289376072489118;-0.415480902912099;-0.239470088048253;-0.128754324607593;-0.22903152303676;-0.115494733248293;-0.273432583225153;-0.0686319839999293;-0.266476298088247;-0.367723629486796;-0.380275198228709;-0.35571454159205;-0.158748555150473;-1.12315337672989;-0.344948090742903;-0.249354457406686;-0.452406562177771;-0.285498967699957;-0.124436117168525;-0.268925003620865;-0.116548916258759;-0.292307554823078;-0.0734396322975566;-0.308891087525364;-0.249400804771191;-0.309641993258685;-0.238135256197765;-0.145843288611547;-0.56531267484451;-0.203286618564999;-0.265900728070923;-0.241489466959267;-0.237946623846318;-0.663144384019225;-0.955322729843654;-0.295886082232451;-0.242465148442536;-0.413401683299765;-0.276501941340049;-0.0916550312185518;-0.358780815133582;-0.108022489823466;-0.341364911978146;-0.0926439001049119;-0.356377365089337;-0.272500855036183;-0.153036696749193;-0.332612334417867;-0.172526876517726;-0.759462565223943;-0.179665354049678;-0.166627337068796;-0.214598513873774;-0.212466732799794;-0.260380905743512;-0.649217044069121;-0.207148232760036;-0.140576690126271;-0.279745345274449;-0.164777038112655;-0.709481515210769;-0.646201656120608;-0.935471142904879;-1.0298428776423e-11];
x1_step1.gain = [4.75187287884821;4.2972540100953;1;1;1;1;1;2.68837434099474;2.22231338789035;2.80928935221552;8.57778787101221;4.3198394271407;3671051121.2174;4.46417948165;8.58719253499922;6651676104.25069;5112773441.4997;5.43390734620574;3.56973498246104;7.1764162163727;3.54207184044604;3.05338872874895;5.45829318862608;2.72974166671106;6.60584738400952;2.83956914073657;3.39812307186642;3.24947653524747;4.99219133623163;2.34208909816724;5.11177638176143;2.53732200405894;3.67800176331527;3.32223111453657;3.11329069625556;2.69741047842681;2.30024333058731;2.34602788766431;2.22890338287927;4.7781919141861;3.63651069534184;4.74581665590726;11.897087326475;9.34749653905283;1415509445.87335;4.62412444608399;4.20434038816018;4.06706890969944;14.4680398291659;14.862192474245;1046912265.74188;9.31892842455341;10.8286892694689;810818709.289397;1372057853.94951;2.85550201512499;4.06069397368422;2.75971481833753;4.11593753642651;7.63487053366397;4.17875956407209;8.7184315772816;3.72487766545889;14.00450135311;3.76421701824463;2.75630337380801;3.05616279468727;2.85260620922826;5.17366160762839;0.920109079538482;2.94365803758406;4.86033285117087;2.68825698367066;3.91333636717379;9.29350447546835;3.89166774432037;9.7232018652409;3.52626721894329;14.0414150000766;3.31605654613932;3.55818041048953;3.35937884171942;4.47136679339308;7.40393167549623;1.64724051672711;4.62840259345926;4.11993170567037;3.97919550487344;4.84624814604425;1.49015975157078;1.00626892999905;2.98736742319123;4.09639197019122;2.48789215348642;3.67957868266444;9.18733279633873;2.92530561481105;9.43786592079649;3.04080260746139;11.2344093115627;2.8047737991279;3.41370779623295;6.67574820405393;3.01458548196527;5.82202531287992;1.46238972510508;5.65366622211361;6.91473325460636;4.75287893701736;5.7004169640132;3.14738343453435;1.61031677670911;4.31072428558303;6.14279899967542;3.52469183902657;5.65778913511327;1.37749708234755;1.47238435846855;1.08036320721073;93726384270.1012];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.2766141668812984733;-0.43555778871300820931;0.43639049127375495463;-0.68087286268884228413;0.23412741052343086223;0.17042293184227835789;-1.4197323856608126125;0.42719995882149847111;-0.63723680887915223892;1.1447652997780322082];
IW1_1 = [0.46960490812657268433 0.23874964585923594274 127465920249656.01563 342425638068581.375 -2.3680737142635135406 -0.12357248646146212767 -0.13153024166255458671 3.5533560154970809108 1.7218693335954786683 2.1730060052859836439 -1.8731853914556206853 -1.7572640841751230312 -0.42221969870229647537 1.1567285988017894649 -1.1208038867294300989 -1.4098378000906146568 0.43046817028721451326 -0.41205503410410670995 0.40300847259884070573 0.94569942708358833983 1.2328797558038264714 0.033987758803398437912 0.058045179449492327606 0.35148215517554898257 0.82083100112932139147 1.2817834522490396676 -0.087033265765660172697 0.58468935493521501101 -0.27716499169921954282 -0.0431088954662315893 0.50683823308417086384 1.0160281315122707824 0.73286402213940804451 0.39562949381324802323 0.50883357641799253912 -4.0489980490300538918 -0.78727717096471150349 -2.4892121227882122447 0.78261315379224927202 0.12972173987201507384 -2.1387343437123469414 -2.2354841448621294298 1.5652853860883066961 2.0524012867469014765 -1.2445572210330599905 0.27879838449301441727 1.8362290241570180527 -4.3037441949211254766 -0.16034213588084370028 1.8971106954436829373 -1.2170248271524970463 -0.89250846842650710666 2.0412595995949849303 0.45496934080852924343 2.0715625204255192493 1.2548479688387030873 -2.0409962764493285725 -1.8894438484170592574 0.25011123665265411198 0.49282991630717881693 -2.5691943824940266516 0.95640648794058880089 0.36056353928552353771 -0.076753444218340599603 -0.25179073476656443908 -0.59593050259323820939 0.88310115643037156996 2.6266931943643347935 -0.66802686405339062325 0.89385612531696923 0.029496474639242237226 -0.6477026570109817305 -1.9039767663491418137 -1.722745818320088107 0.25227962310999119122 0.68777243176999802721 0.1119946112027444679 -0.4657738239249987755 0.07567092185057180953 1.2234088671091096767 0.92032669322645166421 -0.94086927889597238295 1.5018214548950858944 0.13701638719736772276 0.88586774362193754229 0.145659749566046548 -1.1941545135652380516 0.32250465959558671791 0.76979590854140655942 1.2670037041682626899 0.41690373198327074222 -0.30768590687716274035 -0.99936485849590994679 0.76589209404071634246 -0.38762060151051863643 -0.11457073512216214839 1.3871008129959569999 0.64253093367040925443 0.64048286729157299924 -0.30678544533030360775 0.0086234912567814819839 1.9645123329870641804 0.82831786932547546431 -2.3926480210654497149 -0.093154373564907716476 1.4132125804161446148 -0.78982867149462532286 -0.25877097330742432346 1.1147473654422723754 0.85685797564450449926 0.06406305282118572586 -0.7123393725548948785 -0.95701339748752700665 -0.28562106170250678439 0.43483772054987412936 -0.055831721034898536271 0.55299828540790141318 -0.80700982347915839465 0.58070130606407555085 -0.18654798584492934821;-0.084179172959579717794 0.26000090786804819221 16121501542820.853516 191905620442671.5625 0.18737197286843704225 -0.36989360625414369199 -0.1498654954872498668 -0.51759013987798807577 1.0753537470982001967 0.1836934064541071443 0.14711161639837649906 -0.60631243337559614925 0.16240400952629696629 -0.22075920976107077376 -0.041958592625949475696 0.055041100848297358494 0.047922867223353572685 0.15457392552518905982 0.03466292467953865919 -0.0233523361647305569 -0.099480982583400073826 -0.18795471293052842809 0.054208342431866854416 0.047084368971532425396 -0.01483181485946823662 -0.14660464924435490097 0.050892237799304421153 -0.42482458432640207713 0.14307230377021296075 -0.082704852234996725779 -0.17514170820701663311 0.40237574476146464875 -0.094122905646579713146 0.31943473966488628024 -0.10072743990236046607 0.76985812937447772875 -1.268896664742042546 -0.10713295947163706257 -1.3574692120415920993 -0.79962572448844071449 0.12428778038776169823 -0.51614387284460183647 -0.15120993499763091261 0.23929470163842578301 0.12315152040525942367 0.015243430165763380973 1.1630189492924034589 0.43557725654649720459 0.17972774395474788722 0.095561349160919253909 -0.022579410438991170951 -0.68926493145158107545 0.20413994938510682253 -0.0017007292188823989673 0.2479431246047417714 0.038745271059951742754 0.88168080310354313145 -0.084082240632963753391 0.33224396536143119274 -0.69015133941728434852 -0.19739587279628553751 -0.29346965120021273332 0.1857099645575188418 0.41919270378181267001 -0.11667153172876884581 0.17063644083878964919 0.26676141708834893995 -0.16745929280644414328 -0.21951637525588707511 0.056859981498975004088 -0.10701986685839769142 0.12841095625881873699 -0.72723377633087582517 0.37087308729059298429 0.0018679968407598317864 -0.37489563073460130393 -0.013103165473026458232 1.0207824597080594753 0.11220179235234871429 -0.43824717217271957903 0.40733079903981078695 0.33641335444062175908 -0.3312651022924354427 -0.22943339982264313037 -0.10589765855283317375 0.16213660676893160151 -0.065399139051781088328 -0.16587776044877580728 -0.079698227632895177197 -0.028771391689040699824 0.17885877177779255232 0.44104250981060066916 0.0041914479864771300721 -0.061949282015470556739 -0.91303788230238192014 -0.12994237458981500177 0.64523128103633042851 0.10037286291165450902 -0.81425099213065221804 0.19418276346658178144 0.97120968701209076102 0.15738214966937677519 0.1383728951059513268 -0.0043271578341588998184 -0.35133221023908789427 -0.048156100222269862532 -0.01357784893500204107 -0.11465339822655556956 0.074392417481103187415 0.086845967252331485176 0.013127161710156706048 0.0095247740885759518692 0.13191635429282677627 -0.0072454628034767896982 -0.03222826706399233021 -0.058684312238862036681 -0.22706671338677311689 0.15086951013549554368 -0.059025705193341063137 -0.1833097923639009641;-0.13944930537467287168 -0.27788817917414410097 15651767456462.472656 -327560859469778.9375 -0.50853166456163922948 -0.51428094072687935423 0.27426755803566643666 0.88381178355072853581 1.7421700460122486742 0.37751168366382986674 -0.99451773771910745303 -0.8020412232490526927 0.2841114029560539711 -0.40448828762921401614 -0.29516065816174469338 0.058522449644255229717 0.18072133131224399349 0.047741683830836864622 -0.031246019126249031939 0.034298148031768629918 -0.25029514820673753883 0.19896367466710049299 -0.24742838922795612477 0.24851515434270551874 0.04262566531328330216 -0.4815571988198190656 -0.21476625334465698769 -0.10895433413363825315 -0.23970215136594807959 -0.31526893334216005149 -0.20669947795198714835 0.30472975227061893655 -0.37907681930539993687 0.42960195946666290823 0.014106292507382489632 0.31773429328850905096 0.14849013543757075273 0.074146340361765719629 0.57946425536998491523 -0.15823343158942132503 2.3976023302698745709 -0.34014427944811975735 -0.39032613800503329626 -0.88142612688056509729 -0.05212106851038937616 0.22224517941894117445 0.4160434654135726773 0.27690150801199819863 -0.70126780687684930893 0.0050484506265423669327 0.25569922928515770311 -0.52497646446790002805 -0.37065473732756598624 -0.10214786458018230986 -0.28203653263319977684 1.3702187258360456479 -0.031105685908359850117 0.33679419007591082647 -0.50346782729511341259 0.15495359773179015472 -0.18922019031707371095 0.053736506973307993085 0.11164686364108294958 0.18304781685846013839 -0.28478356465161430711 0.14597270693189742796 -0.040595534556305060137 -0.49721814928682750434 -0.26238769660973038089 -0.35677208109157915183 0.30803459880885281086 0.028448184692292540626 0.55693460932817651088 0.58657011746004172892 0.32026039058938127191 -0.7423414909849839205 -0.13966995100982207423 -0.42703312470968363135 0.073266090504067085676 0.27268664633858374913 0.57154595030621868634 0.37254276641300643202 -0.82016981267999256033 -0.28050513916533637016 0.36573631542156059826 0.33635419487118745074 0.17125007100060168908 -0.58643013708643176862 -0.1122087863111787126 -0.020863660405280336219 0.48144920920030381817 -0.64892270951444674942 -0.6433347500181451073 -0.01244367989639124511 1.165087223246795034 0.21332619420238835284 -1.0856545296052213523 0.37059933598912853103 -0.23350871385861682317 -0.086299368024406922228 0.012547962769528355279 0.40168166002161853623 0.53063903975440929184 -0.26299959482558782353 -0.63689126010839614001 0.34339544985295689106 0.00094005112499825601094 0.3121418961010541393 0.10192006348448845454 0.010843854197500849001 -0.051592894233830059736 -0.31801220298436855272 0.17083424635614902209 0.31602291027967616133 -0.31596168457126155849 -0.23242555755616567392 0.60745618968581149844 0.34924635673854620421 -0.17790744978463060577 -0.28669866079029743267;0.20012493532564520904 -0.19748254742423235197 44653179758308.851563 63733799291987.53125 -2.3961373225939248677 -1.8558553037273533004 -0.21053806942094352439 0.53667098272530988634 -1.3868537844936317249 1.8205052375240495799 -0.21270593358006925633 0.29093698244208487935 -0.071739647271484643931 1.2479441853636015214 -1.2332587419175307097 -0.18569057881294429024 -0.32090557074836034701 -0.44919401198533842612 0.043953301632591881831 -0.30944050811306006432 1.0003006177497220097 -0.28293019500385235832 0.043922475370912102133 0.0078028425720636914939 -0.44672450909643041106 0.77875971431100399212 -0.27356103929563840849 -0.21236740186078789661 -0.46521700151810668755 -0.11107160264921722703 -0.079321316913044073726 0.37071413052093316054 0.011620543883138576696 0.030582395512080175204 0.088908910989840328476 1.2325970630064115863 -1.7559757080781823735 -1.4143263553256317788 -2.1553592790832292181 -3.29114897234478887 2.3211627509831136429 0.76408448086104763242 1.8243937109166403765 -0.87564157526560426614 -0.51430811286872535337 0.2549762268567983825 1.4208205983303929809 4.6582670026847692313 -0.0074823948668166634368 0.082623946179985141569 -0.55589054480749255926 -0.70905415537074123478 -2.3115203923233047512 0.38342004126036560718 0.45786154009389845587 0.19449462007961362819 0.36149556726639731075 0.35721176886687677765 -0.14799839171554426964 -0.29203000879445273341 -0.9310090872067746659 -0.46666073199691399109 -1.8111613420230066041 0.045356196753543319422 1.366789915173544756 -0.43987373625188169735 0.67094485139840265298 1.9888087936716727366 -0.55820897996920770812 0.72329111648770139098 2.2047863749670262301 -0.19046580518547753802 1.69511620876044522 -2.5114724851115277637 0.30020533403085969404 1.0062628366592150098 -0.548720833446910361 -2.2827003499744926707 0.13094172722027772049 1.4273524169089113478 -0.093322597011887206686 0.105240653092926989 0.68371391016432703225 -0.14456121015139006269 0.86921951452795953585 -0.12394713513395311244 0.47599503253108671474 0.26734123739041315559 -0.53782378605947345207 0.9387140021638283871 -0.13202179412430892214 -0.9420760006061607017 0.1541971459999467442 -1.0433390529200481378 0.35801746367157949447 0.11213892502031135023 0.093845002816947159507 -0.60129565199881840343 2.0502767666880838959 0.40099412366389536322 -0.63244980488149293851 0.34994854418795917628 0.23909234022889208626 -0.1132413456333916657 -0.38236487918156303323 0.53396752773102740264 -0.71380418242372800641 -0.18016208749240125719 0.769646374563000224 0.24772986341656125231 0.2776643969202577944 -0.19627649993804963113 -0.40790651879453510009 -0.036663718404101000581 0.38895971256774847236 -0.1335539111018314018 0.3230284227567741806 0.067994669761601470026 0.61886186216279159833 0.085279653107431310777;-0.046463237846831526023 0.043195024858424496683 -11230160784144.380859 84043629160563.578125 0.12325403311584025301 0.0062925015006530821832 -0.16915010524677523418 -0.42511340711256578873 0.62998012016178239758 -0.23450901745435143164 0.13621127270046287649 -0.37020995975659393151 0.06509376251502412003 -0.061184386470486996912 0.15315973182908690364 -0.089543772357727424716 0.14110258195351785271 0.09066090144728151945 -0.05667095592120231784 -0.074528253007458175605 0.0033999735962634307015 -0.085729728412170969065 0.0091438318798152569411 0.020166689135665140375 -0.035628255176268462667 -0.064435288932947920593 -0.0051067830408751818672 -0.068168302155663135222 0.080693708523506560426 -0.0029835127653249004717 -0.13979000654345433152 0.1205247697729898565 -0.017243985149316295308 0.040378051243170301476 -0.073851966517716444161 0.076304840138659874027 -0.43347626930586508065 0.20436951687927143451 -0.37113898844299153534 0.050625762342914507053 0.25610320851065521142 -0.42922159724601205433 -0.24747994604585155098 0.053663478824006530243 -0.32074357999585095813 0.073603445821182000852 0.39726300176785617202 0.07828040190379789709 -0.019492802826100476266 0.081583107170473970204 -0.21036780774353325274 -0.30633929806464316181 0.10530264647414071255 -0.16133539483209816123 -0.17763619214234210864 0.17815839644437747591 0.2028026737732122331 -0.25665907379527308496 -0.014938621691633865649 -0.1795751863411089444 -0.005104663037938135875 0.0013594183670300053878 0.35916522137340872733 0.13067766950693568417 -0.23492998202448106548 0.05977546195234327564 -0.075018122246243143514 -0.21207653991569522245 -0.0017793479787326610224 0.054187088978049011145 -0.25001956691171434155 0.046792993250967335195 -0.514868255875805092 0.50877193410695920051 -0.013978642226719129038 -0.42203893641916906398 0.03489330634938563741 0.76225749090795857921 0.035489267573894665853 -0.48271763505991766019 0.20623696715856995598 0.27256666755541930636 -0.18397103424698724172 -0.15646455290973795837 -0.0046825309243320283048 0.11594084004657707887 -0.046567427252269345817 -0.18608891058092436044 -0.004628185695892007151 -0.13172222018457691028 0.26654645946588384486 -0.049937992364504212028 -0.0024962141608602093391 -0.089467630310049131581 -0.050909045183725132322 -0.025615138827603217697 0.2320419918592538655 0.047871080988125251643 -0.33093524061817247706 0.033850363103431861667 0.48830031687965996801 0.059222300366834540464 0.053218406973369249136 -0.095877630702771957094 -0.071137677555667896967 -0.076666615466750395735 0.045237973635772192604 0.0048212702716204820785 -0.01078360052768366148 -0.033660775393469802297 0.014123975148286404407 -0.16688743789686252827 0.13396490567622121315 0.060829466469576473631 -0.16284684119418385562 -0.093044977929007274819 -0.072612710510125258456 0.14308143871561343485 -0.073270390060861792514 -0.14769992897543052224;0.14460696070506484645 0.47482793083143215451 68535271277421.640625 101610365164100.42188 -0.34970177767865984197 -0.31160138082880789456 -0.33562662745252219709 0.097618262586479223608 -1.7375456329909879649 0.66402537818784745927 0.33564697374571489386 0.61891495632370574498 0.17927537685620439922 0.70577880710988694091 -0.14032340861170877755 -0.063409040296426077976 0.14146795085035307471 -0.1564026554888220677 0.27183769366578064419 0.14797537176839484219 0.32266754533343888767 0.01453290343499975823 0.17225226569172960578 -0.030990856210077584454 -0.064917005583322542894 0.32934784906337644106 -0.22701139670217504452 -0.034164650738868165836 0.16334625974980368501 0.18052859412234353864 0.014328565897995464151 -0.039600430802951835929 -0.096621810758062864899 -0.080302610323754700095 -0.0023448473925667470685 1.3325659032026226924 -0.72844564282393775567 -1.3260236741960058637 -0.88621453499846725155 -3.0176843087727300841 -2.2833073175500033791 0.26040499591463833751 1.3731596857342407425 0.76829455222221265664 -0.43508808104730711808 -1.1042251337459214611 1.4116095368486154538 -0.45519697595440883342 1.4325216044375395086 0.0855044951311761392 -0.30923542829008909827 -0.42310543934853916959 0.46778880637157799205 -0.18459906663167824359 0.017602566691235645158 -1.1557759769838544806 1.2295144170705669673 -0.03655772620050012911 1.1019079779418154974 -1.3324034228985079942 -0.94955620783017535835 -0.5466200916905330498 -0.88357517186846812329 0.346728928995462371 0.95720839265059720624 -0.40305698871996320509 0.80064498708441811381 1.2895446722425165653 -0.20535752340400126315 -0.13784808853198632361 1.0238624008325252213 -0.60996449744017333661 -0.32912779968254451957 -2.14558824362981726 0.17417965429472920258 0.99027990167588275661 0.13844232586978225119 0.41748291170580387188 -0.026941577581779431461 0.0040625300937764223519 -0.002915863566518137398 -0.44867049225262140011 0.66714761291509472318 0.23131730701762934221 0.10080551838452628255 -0.19902565305213790192 0.053478155641082117389 0.50097952442883808555 -0.086977823394606593177 -0.25473373561068990423 0.40956863081969646645 0.24221285709304901657 0.085286265051875961851 -0.91955847875069351982 -1.7135398128782775018 -0.4315389014545126245 1.7517695109116040442 -0.19437282555432319753 0.9957972017753240479 0.66011969163991890408 0.45474564260513361535 0.35629706771528818487 -0.020036063367046184902 -0.23389719571673256437 -0.26732360732260679326 0.31036854459765012004 -0.27579073900691619681 -0.479465890902889047 0.18704489701680640512 0.27135838453405453086 0.3002521200055135453 -0.11783035683529591298 -0.079427426971000650302 -0.24341355688885765729 0.33713492672170042308 0.032704573406117640189 -0.74693878643334310485 -0.0778531890774444707 0.26437697292064182841 -0.34444208883160698509;0.19000872343500274653 0.43555056757607413243 110025548183640.5 -380807616534591.5 1.0646056235188237071 0.52448054704274726934 0.72874691997808493937 0.86672679948601905053 -1.2235182844318055473 0.88539259018936977697 -0.33488311052986852179 1.5233047411051032771 0.47718332862420898 -0.68227347235888213817 -0.29477620615501298751 0.19060564664604870422 0.038957440954880129258 0.012131393911988704037 0.27902186702112169447 0.54462886681734068439 -0.89872696895665327155 0.13199337885442202767 0.085885575555988288854 -0.17958453769677434897 0.43073241157677427671 -0.33845730676449597985 -0.097117815531606102586 -0.1232867500925901999 0.23011690216160196609 -0.25141690648400022434 0.51451999084198096757 0.17627992984919491071 -0.1196274531438537847 0.12312818056849615667 0.27693203646364816795 1.2313192555783176729 1.1477082641356959325 -2.2777759257735508136 0.58321301076844855871 -1.8955455380890582351 -2.7501207863308860802 1.8986309101388081366 0.93707350052021898268 0.59101886660359248538 0.66192758575241295027 -1.4841645003511760859 0.51786614292134292015 -1.9684045788874706684 1.0386844996484607684 -0.58791238840114723452 0.52317467774765535893 0.347852987596000085 0.8156663858137375156 -0.28902785881349807307 -0.52924392161448241545 -1.0866676965745136751 1.3725323165768572586 1.2380108942851482112 0.79252741973984319745 -0.93886119106570298243 0.76713662178943597958 -0.88787624840824808103 -1.1567010518786380313 0.64911460198258752463 1.0236612968758878051 0.34796588002905032289 1.2516066794616058644 -0.32461455548881967825 -0.60020694427783571889 -0.27301090356373558521 1.6131111682367904603 0.056281313960220749615 1.7535382751878096919 -2.3333284906511835288 -0.044416661119607611652 1.6573657272609951807 0.087860640774100384753 -2.0311527440855061677 0.020569932846846522567 1.2135597835771199104 -0.35001861746773288431 -0.81496640016849197785 -0.18300042657614903896 0.53055445743044338247 -0.49631870262441701858 -0.42632125811060145315 0.11579827129340278247 0.96006252041416229037 0.01764630376869860337 -0.42196548663320437855 0.42262995031453665096 2.6286731013133679546 -0.86191487055656945682 1.0400460233455586234 -3.2253835863270090556 -0.15262951272769495481 0.65086984114118084932 0.81242570565042193653 -1.4362875368394498743 -0.055082782760309312309 -0.062510121009965946559 -0.14730848997771472941 -0.058982040366332033621 1.0253163434480934146 -0.44046376771342848899 -0.28540874965266765928 -0.060834044047857464832 -0.1827993037208368865 0.28253231955895624639 0.093257955426045896496 0.043469982252508858334 0.34413616124815815667 -0.22507606955977663143 -0.34457291236893439201 0.88686517663879838036 0.46068861746617445485 -0.41721268374073183782 0.27273761887625508349 -0.080453800446083870157 0.33995766419234896283;0.0069258501385753602764 0.073224586694325077829 -43495669889343.382813 310301169931845.125 0.36036863740153823876 0.1510479454274318567 -0.67616948810884680565 -1.2917676804485700437 -0.02099383738147152581 -0.98969394566043800143 0.82977080506283451555 -0.34112069923213261813 -0.16406111643469931116 0.39993358499906572767 0.61103978913907919246 -0.10495541386556503372 -0.1012993869972057398 0.091539897106727982634 -0.21723113393447349151 -0.19796115974113942526 0.38893591353347872763 -0.05417520889073858259 0.038978033321485874851 -0.20496885785337765062 -0.15417704278210955793 0.29234390772736679187 0.087561972425043738122 0.024002029928606373466 0.19389029250370573565 0.26602329825919529549 -0.27863550694440747391 -0.19757868907995809749 0.0062483813642738427613 -0.20964763870826122472 -0.27055981690381886162 -0.10507764037878442043 -0.59878364599120792455 0.87524313437264722459 -0.41401830504773151631 0.46977564124246351529 -0.23898937132680447837 -0.9391625675078845159 -0.36647692738817011504 0.30425080895272621495 -0.16495372940033842246 0.14542554211364785988 0.4234096254543672222 -0.021260466084936178316 0.17507578973715362536 0.22701626254706069874 0.040230222999572531895 -0.30634372650799535087 0.23245846745636772335 0.24844113973328785239 0.067852021217716079748 -0.4008483356523370178 -0.015087721457941762754 -0.64643361827364342886 0.15013325816971220816 -0.1734217333561417218 -0.083808807298619072967 0.20987695416168344886 0.36550598089664887436 -0.14134173736955141232 -0.17847703031582137778 -0.1464528549760671372 -0.49239041187106746822 -0.053556711652908588206 0.47504689275783534841 -0.048243755598772218585 -0.89404938936826749973 -0.12514081426046266832 -1.7628813004477545956 0.63816495385748672753 -0.16209352678604088682 -0.28224865300111573463 0.14653985346565223669 2.0050120766771999392 0.01316966164441347141 -0.99987737198224657664 0.1344506761366324743 0.4297113982109107555 0.27483054332716383428 -0.21556081135597685017 -0.10447610950126144658 0.010920165764889040033 -0.042963436793121242974 -0.096991591868861720682 0.010328694243166210501 -0.23415324557726635568 0.10663684764985986719 -0.56954861665294109763 0.71418027885534540733 -0.64402699692926457242 -0.011744156080283917054 -0.20549916936129580813 1.1469079764210365902 -0.49486126508934297208 0.051029050627098414272 0.30868835707109337774 1.0386224222182500743 -0.054713365811863569343 -0.2089786401825535489 -0.14076573538455783297 0.35219242000537132053 0.031803884519592137692 0.086666641041122885469 -0.10607151126173545996 -0.15723602900961211604 -0.077481453321106869003 -0.014456773228049644292 0.059582602402064997749 0.19904184480931649825 0.034502314750494263162 -0.25244080529479279162 -0.16382304317232948554 -0.30132441907986029506 0.12756662301779056468 0.20286941099241967423 -0.23066681866237007115;0.069554538449853500759 -0.092270531614694187295 -18895659506398.515625 138801261595891.5625 -0.48385833743564676812 -0.014376322169769773299 0.086401576126414608714 -0.42564331522895648696 1.5928190387729184163 0.18293963815925040617 -0.0043855329008062431118 -0.64817654098822885533 0.057526062555879622917 -0.46742332115265677883 -0.12544101638354568373 0.062393673884920008987 -0.062424315040055966541 0.043094053613413514647 -0.11082029289603488598 -0.20479045076019466998 -0.34205335596060160164 -0.5090314312390634699 0.15457562915915015744 0.091852216364207942001 -0.13970573212550471687 -0.16175916436742393545 0.037152093682905312777 -0.23991205250926034287 0.042369157839254115161 -0.36110295733378938721 -0.076833752117311721053 0.67960471872234495816 0.28249819653884850057 0.083680462130322125214 0.10543114608322171066 -0.84556153770028374961 -0.36683975645257249631 -0.43868278011358530222 -0.92042562186380938449 0.72339106735662461745 1.0380726532326953748 -0.07038151528659929157 -0.38407652238932965538 -0.2434252454705464197 -0.11403524051484338619 0.27343868981230989856 0.30468972906174657256 1.2481540843780294381 -0.43013890817082434115 0.01242518562470935721 0.17511699641318109744 -0.29960186542774158713 -0.48132453729540269904 -0.0090028575151864763393 -0.016202948637531855436 0.75804842857660315047 0.22279221676383056949 -0.8055854663806888416 -0.026874451857353782308 -0.058444424169043288075 0.088610433208866329347 -0.1479219583438270702 1.5480769441521462504 0.39614335768597075393 -0.71736758265677047497 -0.018425013182088223473 0.25957490831114343743 -0.53658883993552886071 -0.59295084192037783577 -0.27881215768604228344 -0.064663682270683198183 0.25098349910745787472 0.68623545383689577015 1.0282728893843207096 0.11291245925256578053 -0.86691792142973111979 -0.058530753303321546299 -0.51819011129102443558 0.0025737001844764397125 -0.093677780183671358061 0.69789358210015450812 0.38206385293639527623 -0.99961273000061712235 -0.23227980601855746645 0.25802359211276593154 0.18882872703005998405 -0.054769901621454467056 -0.32046872162152062469 -0.12071121666626351276 -0.16561022358208329908 0.1443874359885304326 0.30233155895898755627 -0.47878896718106450647 0.78753244740185013306 0.044790400763741883283 0.31731102002884387092 -0.28528083073678439652 0.44446449342096877588 -1.5571635143528417267 -0.41476920840751996256 0.84610276930327021905 0.15530535686472221601 0.059768060709657302998 -0.13977402005844163924 -0.035561878818403572233 0.12060743004462229355 -0.20582250594354370432 0.1388494514660539858 0.34494498499897657817 -0.15496021343370547374 -0.14189578137499878063 -0.31719891077763534559 -0.028052469893577049176 0.071345782700688781253 -0.038309310079142754046 -0.011231532752296227717 0.3267610619878101863 0.32784341499223657035 -0.21156475140664060541 -0.16641993578975355916;0.19096921970872843266 0.24882786497633863076 54908601853047.867188 -379118909851175.5 -0.19206041545947549753 1.2158538339822257957 -0.40279748199992820856 -0.89941460766486447387 0.27600705467616382105 -0.78442273146921959182 0.38910596086803950211 0.2087019022883147279 -0.031434096474282718914 0.045623401723982888456 0.66367613625157406698 -0.28382177350334336641 -0.12861117694418749968 0.17641328703094819041 -0.1551156622356651138 -0.058078430352277179982 0.067019993671192726326 -0.22399955567435020365 -0.045927921408750987942 0.41748340310995918356 0.19162803397193556298 -0.50430240190336028672 -0.66458891940505171458 0.71117755326429965557 0.39208621460851555973 -0.31957759599266488681 -0.2518307681701453804 0.43469762839735098803 0.08528459711665377363 -0.62681166734239945892 0.085631351062842656074 -0.58552101938412359772 -0.11843189519760344286 -1.7555083909766227013 0.21954669599690088133 -1.0881309648824513125 -0.88958476699087030148 -0.13955651907172175741 0.47120459568885741941 0.45569545077114609333 0.13242559743764170666 -0.50271813841734946937 1.3580632729644241952 -1.1052163233473182302 0.71494082838298378402 0.083316625285389786515 0.18740205470055593118 -0.62251012490250712261 0.64586667500210759307 -0.29374008662512385293 -0.0026153990463093235858 1.2260458208821385373 0.20189082537972954934 -0.48642024682096507426 -0.38944532870762821331 -0.57860032342171208164 -0.01095496695264141368 0.085088586186609374562 1.1341927704082015449 0.38884612842597493509 -0.55040214816396404451 -0.21458078315922782253 0.13372205565993658394 -0.099774021502974818976 -0.14593553373121048722 0.33800679213178003168 0.43855598388622596406 -0.025218370696787766744 -0.68133254877248172754 0.17878240918986942898 -0.19632526460691163139 -1.0028477185897022306 0.36588480766517550435 1.8343163007128111541 0.029909680333932871926 -1.7910259759553397441 0.36350711424228582258 0.56450572632516715643 -0.40315466036063740685 -0.27242438766250876014 0.056931707782956850195 0.20763485569956269261 -0.28302389127770111177 -0.38732018331544582068 0.36881297734854401327 0.28713078590405050061 -0.6682986013969339778 0.2446857059491623354 -1.3545161505483414022 0.015928433554440533937 -0.25958591520114893125 0.061070607070844114272 0.54689687290836386246 0.98562484841252662182 -0.36067124012980655845 -0.20551840599364593065 0.83383532843515617561 0.32577236397540049184 0.28080579310175057817 -0.56917005681904075942 -0.40237885744748247419 -0.090948543245318411166 -0.27287046328976327647 -0.0084352015524842144378 0.40746162768994348413 -0.15098918009255737904 -0.11122573554981962651 0.11402790838590494416 0.35658634448303394571 0.13832459495807622463 -0.44649993293692807894 -0.24280823849819296867 0.038736572556607576368 0.21859238426441332526 -0.2390502694375639392 -0.0075993080873288654242];

% Layer 2
b2 = [0.00099093398950099442768;0.072763924658430442771;0.40265068307933138403];
LW2_1 = [0.12003036927897615949 0.20294296585042745362 -0.34735070232770209886 0.3959325056489455541 0.97223391175641060702 -0.57733756827151594404 0.25746765075102562292 0.023071743600870123214 -0.69375051653019625952 0.19428493825122100658;0.17282670715790771077 -1.0217876673823851341 0.23623648787455034381 0.30383512902185705329 1.5163100719426079177 -0.16842850981161178181 0.53750573731739248462 0.66251128986880403193 0.083740048942831712964 -0.41470828442589796614;0.1341454756805689208 0.29936568479782316166 0.44169051445212792473 0.12731717435316924902 -1.6913096448124660398 -0.27737250944311647727 0.32913057013212670965 1.1086913307653760707 0.35665257916149190054 0.14401518104045452073];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = [1;1;1];
y1_step1.xoffset = [-1;-1;-1];

% ===== SIMULATION ========

% Format Input Arguments
isCellX = iscell(X);
if ~isCellX
  X = {X};
end

% Dimensions
TS = size(X,2); % timesteps
if ~isempty(X)
  Q = size(X{1},2); % samples/series
else
  Q = 0;
end

% Allocate Outputs
Y = cell(1,TS);

% Time loop
for ts=1:TS

    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = repmat(b2,1,Q) + LW2_1*a1;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a2,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

% Format Output Arguments
if ~isCellX
  Y = cell2mat(Y);
end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end