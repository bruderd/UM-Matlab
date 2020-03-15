function [Y,Xf,Af] = nnet_2020_03_15_16_52(X,~,~)
%NNET_2020_03_15_16_52 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Mar-2020 16:52:10.
% 
% [Y] = nnet_2020_03_15_16_52(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 17xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-2.44476982008442;-2.31300497436095;-2.1603244616845;-0.967066313401745;-0.802125813047507;-0.866868869931027;-0.780848949366417;-0.720433274148633;-0.69232963651746;-1;-1;-1;-1;-1;-1;-1;-1];
x1_step1.gain = [0.393329218500476;0.620771070077988;0.541995878218238;0.835568589577178;0.704385023583641;0.986318377317454;0.864693114042257;1.01808715517802;1.34097346910526;1;1;1;1;1;1;1;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.5916049122491715551;1.8320286928983573471;1.479892402130263962;1.3428318371708816148;1.3494051159439015564;-0.7064290713220007234;-1.115699372692974034;-1.0811923222495176677;0.62490047676284987688;0.87142758795195063293;1.1111940134303326388;0.28287623138531459155;0.49036140053899163194;0.57433455429264901326;-0.24688081296379610352;0.18452734001596349245;0.16146148745894195464;-0.94203025363046111185;-0.041297050694037842899;0.72452494153170621249;-0.81673500567651402537;1.0629053913100219297;1.1566955967216014933;-0.81841738715100875723;1.2035633216870993412;-1.0116528370324058272;1.7479642393496253039;-1.6496781965034243989;-1.8473317843748933154;1.6395178344130691972];
IW1_1 = [0.2197674573090398209 -0.32263192397906731435 -0.041261133803309486823 0.81857373885941020486 -0.013788452578126527207 0.85536747056835726877 -0.27686107551435745755 -0.39796275140959691541 0.016646041341137531727 0.51629403073923840761 0.44260798893760411321 -0.81597418777732944051 -0.062903345983226291449 -0.089534411138066205549 -0.53735261127496269928 0.64567693943128379175 0.12685787112222371431;-0.52579172667276052877 -0.22971257122249194582 -0.67915500934406425326 -0.22437158153230338242 -0.12257934671397065218 0.32222969006569712747 -0.7214506964159729252 -0.30624208407087244455 0.3936040712816469278 0.35107309398169822545 0.48618179297115643456 0.22923239494652944415 0.43739508737355520251 0.030893060352053150486 0.15705403411713061979 -0.26440688318933963652 0.47415107961657326774;-0.44554878672740888534 0.25958967341800670559 -0.35820018587156360335 0.097708900315323490848 0.58883371419962893256 -0.063116526236777631409 -0.013314817400603876085 0.18467900954859278739 0.018661189152414180137 0.18036217030778178838 0.0061316018899978388598 0.51559432341382083198 -0.61283174814146201825 -0.36603961784355687481 0.45520970324400950036 1.0071868432508412905 -0.024739316355096876293;-0.28559724548855647264 -0.030320454762992356951 -0.029401683453618699532 -0.88316951276108623947 -0.19901985444868777742 0.039947004278075667894 0.16312768585072218164 0.027210506530323353563 -0.035327757667363693383 0.23045435733310060122 -0.43876091176399267813 -0.61570689042598358842 0.71666178482574749431 0.10582914703492851682 -0.48423823745727689882 0.39500406697049800231 0.25949680316703938043;0.39847095220705180951 -0.1303211170313995082 -0.24120415902417854936 0.28989811493140632548 0.076618500733666619618 0.24585168638500801142 0.30178265611611543306 -0.0041603822286856275217 0.43662776251072621392 0.50168552693600632075 0.58541030473790744981 0.35021751894833746954 -0.98215389916918671709 0.62718954197798870709 -0.3217202607824035665 0.37027937976305985712 0.29262538978535485024;0.19055919581265609586 0.13345770336693257851 0.39761648140066052637 -0.21356226938923805747 -0.56202255005103840535 0.12430260818478074125 -0.089417333567401877104 0.10604962446621002403 -0.6761036862559497429 -0.91804502140396204712 -0.6971095855052541479 -0.1652264117016942202 0.68989375188988200094 -0.29167091072894341108 -0.044288418652863779934 0.34033631523742630121 0.29860968461514186822;-0.6080317640158854342 0.8277688447087332424 -0.7631089301934714797 0.33043870979699341284 0.052092329501549253234 -0.38444699416722249996 0.08455948355685166673 -0.38002904361843936609 -0.13469790076205725993 -0.18171002587445830856 -0.44498427007766555574 -0.48901317729352111741 -0.05319708692786197235 -0.82986769271770832734 -0.37696155999950370363 -0.053189942676566591417 -0.35195768718938569286;0.26098894952896839472 -0.37586588079904620363 -1.1924600687089967899 -0.051455849708997483671 0.030806072780445538373 -0.060806664887879009751 0.072151687742781217549 -0.26131045739207320322 -0.29294513644943842579 -0.40558520448998797159 0.073090691079641958994 -0.10482453824272043585 0.57723717607543800767 0.40470307351366729964 -0.37360182610887943122 -0.35410343907387714646 -0.38698890324343937497;-0.13800756171412967288 0.47231374062762665478 -0.23445021395763945837 -0.35344426628298719617 0.4202233474184778883 0.084532266826140151084 0.015385577162176689719 0.47808024252271624777 0.1428641597807486785 0.5525044098261405745 -0.30507177821291042719 0.11041304788175067275 -0.088338301103379035517 0.52971639673637593848 0.67245549615227007134 -0.31351326478541585718 0.3202508508977700008;0.042994418180346737068 0.96922133332822024698 -0.16981010734276971719 -0.16522097458662787783 0.32770235621619381305 0.20920743957099813271 0.20939833509466104777 -0.10165210805065928723 -0.39318256087875974591 1.2335106575542749141 -1.1123563416325690234 -0.22518406159270859535 0.15876520054289580686 -0.53187661383590512099 0.3765205343532050386 0.65932782977160642179 -0.47141535624842390018;-0.20366035428108988925 0.41119085170371133531 0.62807703661416192098 -0.075768059021777411921 0.18931169208307510599 0.24220463797497737701 0.13814712299905915205 -0.088079452097140950584 0.32834442632513571869 0.42716387944274908373 0.36443711842720477589 0.055030230200090006487 0.47867465909667694657 -0.42165520355468849578 0.61274957942409302714 -0.18957980347653538566 0.94699999703724335731;-0.47106166414769873318 -0.20691872049881496198 -0.29672322609906276458 -0.15625947997139214007 0.18985989883539475254 -0.23025903575025899461 0.31237949654610069894 -0.20486754906363777162 -0.0025815651701540558507 -0.14562682727298525975 -0.16001006784026142182 -0.65824071486795354691 0.32771178729880318103 -0.35155166548810667937 -0.33612630946536042886 -0.037225838380660319915 -0.30440045196961934559;-0.56908792095894911167 0.020611638361359096994 0.0095217373740981922636 -0.094195492359024302065 -0.38093581002633208898 -0.43914840481689726426 0.11611469665691986908 -0.039605642708364821591 0.033855254366922438169 0.32513152654114158002 0.4086293790743967258 -0.039949875186928669735 0.20355016863009808836 -0.53874934425447751352 -0.42981768776534617915 0.52108074016176508803 0.42505092604288141622;-0.12687691746129101111 0.0433047285684709099 0.2745892884250483168 -0.030876938259970777151 -0.13426675721735184532 0.22202819831412162754 -0.17977582078338058547 0.2504155003674503055 0.25352482820209421011 0.23459298580181575655 0.097693076932696507053 0.45597371735213620969 0.17712096103854019224 0.057452156963019224145 0.17281881898240247786 0.10526944968344041642 -0.024710927617774849224;-0.23681361116397595112 0.4019664520079254677 -0.80215241494942657319 -0.40751608095547120492 -0.0378403747714291111 -0.1542581654836611138 0.13279688375489606678 0.25057914268567393901 -0.17713998935405178559 -1.0341701310314519446 -0.53544940145123320896 -0.066380706146597023842 0.66785395385471857388 0.55743730143910896757 0.20984001817234759768 -0.57511287237124064653 0.33338186455082352566;0.13666752982032145947 -0.14003164642404133944 0.19564002805803162865 -0.096145793250038080746 0.059427830039461296274 0.70331001176535423713 0.26276104424825624806 0.15989475942512074158 -0.081403753451674956687 0.018902543494068883845 -0.47123242992249075867 -0.050888253714574969988 0.82696619981020580425 0.47895963878117481238 -0.27432297217343248885 0.092265854847138464856 0.31115322091456054654;0.15743604405950492597 -0.4210977060087565671 0.31989287944202671188 0.30538931006459679729 0.26369121964261182933 0.078893896477883904139 -0.15451688610446562278 -0.078874180574975541469 0.10023237195949373513 -0.33067164036168106289 0.40696789481263678789 0.68321934495561853495 -0.16776746381018831089 0.45313342295304770424 -0.32520805524866625014 -0.32550045645443653752 -0.5235341840402730762;-0.22480886904688671124 -0.018465631551592344428 -0.30841701479067784586 0.17222005611712251638 -0.41428921667381490845 -0.11254563776703703504 -0.068374403357391630442 -0.036789528049015363109 0.57442537061591336656 0.40026554818461235197 0.4221570178355553149 -0.10420911224314388654 -0.64243420941945716951 0.04230437440718795572 -0.45688869258981523735 0.27308617803992413231 0.089091437236030926461;0.329046949657936183 0.31704269024575071345 -1.0867510630007535255 0.31907397622819222427 -0.17069970517758734085 -0.60249196183629982748 0.33873511009212814438 -0.064903036137163722108 0.22218583727624030555 -0.88769765198913497262 0.10916528110477093094 -0.016911615508185411788 0.18506980840904446728 0.6035442808082429611 -0.33667496206557429028 0.46306845366412763276 0.2251378787183154484;-0.2227516505620165399 0.23924389837606416243 -0.60481871898746009109 0.24187347699240452359 0.19719904701169735373 -0.50477692619167180421 -0.15306968598249581781 0.20063161515221758502 0.31723648547211874549 0.49452248007575327904 0.57367467454091924939 0.53550121469504630145 -0.34568452545484740224 -0.25610353926003187608 -0.77914878364447981429 -0.35017592561064769807 0.31123209012735220114;-0.04566743290471104999 -0.83537371848498909177 0.81304412778589618682 -0.13331442043944102815 -0.3333029985512118154 0.35870504697583061038 -0.40571617277033811266 0.17458552532330295914 -0.11999893641066988748 0.36180009180081174058 0.32631391900061901268 -0.60163381022282269583 0.49103230350462961917 -0.084402035672473313399 -0.79140175924206757063 0.15544456187456709428 0.38128814369927999106;0.059703809549412213686 -0.10729913809080171561 -0.77888642824426379541 -0.35409792202939172068 -0.13121668329223018512 0.024954309947474821729 0.16157967056361197122 -0.14701928176815193261 -0.14997682573314549681 -0.68330117986299854227 -0.28800132428147146024 0.2510702565714151846 0.20942888780613327926 0.6577345683775853713 -0.11996644303341122417 -0.47020269273494647022 0.45023187502593059461;0.3370178233518652533 -0.28662829330099692804 1.2544380992408328268 0.078884120634951052087 0.29774035250582336909 -0.39370282768828568365 -0.49310196039433423776 0.2165997441122188838 -0.039406929623186377765 0.18944920987401245305 -0.52336924533068041043 -0.38937273795660570785 1.0343685951116179389 -0.39670964551268456111 -0.48457327096390806487 0.18297687351677038858 -0.66809353530767223717;0.20398180451690745008 0.0044320117065257376779 0.62219371222520769571 -0.68081026926410370681 0.40119650140405804795 0.50677081279953639648 0.4242603766271798027 0.092408413072540177002 -0.15817377793334116975 1.2306143332160561776 -0.4851350032159733594 0.096603242150769361163 -0.072262425389457748959 -0.16197081473816021813 0.18347905613211665243 -0.25003234139273938785 0.52721743291278211352;0.20147722435048232015 -0.38633275451502185893 -0.81274897795013267832 0.33583390805787360911 0.19191853790536986257 0.70887276898837736638 0.2658868776176451787 0.22434276263915128857 0.016927638685806906116 -0.2529500687715203755 0.024595029055444635885 0.4046965704988631507 0.36978408419463110191 0.64970841436483528497 -0.61762189574997505837 0.20804371437925420207 -0.59275302162011933582;-0.19623983879517775675 0.34040635706215716105 -0.3195636696595711701 0.051953076971950784557 0.71764377481195906761 0.44173132138661785895 0.39628692657937408805 -0.54025952871872973127 -0.19486953704512413599 0.16591866927642417995 -0.35215719840223247417 -0.52339521603719885157 -0.078199580100276486205 0.41155393795752343022 -0.20645516210639880406 -0.63280151070812928893 0.18054564615641149894;0.091647127502678379929 0.1868294918637844193 0.10455130778606709641 0.075275350060466858504 0.79200726545761390707 0.193648301029772818 0.24087095256646331998 0.048077573400090591049 -0.30170763651140386497 -0.64921847649444841544 -0.18037544085631798674 0.073932661376071531767 -0.014450857185031758256 -0.26830391234110168019 0.38900240088983867315 0.69246791965339948849 -0.80383191088580907557;-0.39316716951417912673 0.93731471667595189778 -0.34337818223454874556 -0.012106020893569213737 -0.68495647684266702804 0.19002040040572126811 0.07209555183043787896 0.14250753526846707264 0.43383649766419035698 0.39636032169162588312 -0.22037386430276792448 0.32513691362543734042 -0.44233516937058009466 0.71000528447383981501 -0.04368081322470628719 0.11449672581704395391 0.017740478342072985019;-0.35591281314954525961 -0.054979647642291920584 -0.66690775748646757837 0.14985120513047656132 -0.098764151117188753171 -0.070910455205695263614 0.099368706005402415715 -0.43462192298238361277 -0.33683403590720062404 -0.046868862631782311212 -0.31370366579186009792 -0.016423533919571520201 0.39377603950470518868 0.083547148849190219777 0.82361577119943729652 0.48735648049261126591 0.81526220637914681433;0.24516422565602116745 0.36594598587002358236 0.47392254763055824407 -0.026463237887750946237 -0.097278105484463014974 -0.04264963994248946455 0.28669857840409668226 -0.3361958587422481104 -0.2584539105699693029 0.79217690151462538672 0.065505950305898189012 0.048634946441165893827 -0.95707389116339880442 0.18074222647348395232 -0.015072924787825583981 0.066455587788524297932 0.20757273603040188314];

% Layer 2
b2 = [-0.7783965055184880466;0.44906799689778115203;-1.0195996494571928359];
LW2_1 = [-0.28127002630483566303 0.020803476978158888522 -0.20611827604821283955 -0.74520429455345182923 0.11947417919537113573 0.36203103800719177441 0.47405327902500687953 0.23509562977330270006 0.53525442001265810532 0.273557447700492673 -0.058098194891525692385 0.2184094167539074316 0.79067258937469331315 0.34719449641639105186 -0.36225715266157232364 0.61167739308414281751 0.30851396174699352581 0.48283074396766911818 -0.25033457340577436323 0.19716517322376517218 0.12663735894218466904 0.43431269530829275105 -0.056120587457643991702 -0.21374072143059141826 -0.039803371229251850727 -0.13031788784076578369 -0.21470263378882145644 -0.49327167111341785333 -0.0035262796716401681799 0.42236014959027512461;-0.39143785885143800307 -0.11793695950034296938 0.089141955960608421083 -0.27399129587604931224 -0.51324872907484586637 -0.38061387577381361469 0.68964084337904618671 0.21952425359175867237 -0.43919467835487890905 -0.70848618066158053352 -0.22199844671764923376 -0.85031231884905766893 -0.45236686051093433214 1.2571025711291334837 0.37138465376259133244 -0.17269060232583749448 0.25278519206152361143 0.62319089131669358483 0.52059929277541339143 -0.78001497116942053456 -0.46108113042221243738 0.30461227944198737694 -0.14075324588329093944 0.7122069832666925171 0.063491256945230722941 -0.36623107821382550053 0.50878226252014280284 -0.22385775086459688832 -0.21577500278302608194 0.59938000904514043832;-0.82841866216126869293 -0.31106995414865445948 1.1472217309053380774 -1.6324245414385132324 0.65751299034732346716 -0.23006280492578334851 -0.31369401501691768219 0.84351392781137679577 -1.0091408761881708323 0.57352766849541225724 0.16273145563299748484 0.70949833583881349863 -0.46546173018148290224 0.49578963350301930024 -0.45609011393938964085 0.34875116952126750114 -1.3699036518275453655 -0.25414371295270854478 0.14149882239516547999 0.4929969056565921548 1.1018976715875874461 1.7338895331179442483 0.64771157116785127439 -0.50400822546207502128 -0.63410472962069330816 0.21296302145422779661 0.6230168080631273142 0.8221228133846590902 -0.55303383588150556438 -0.52012870563787827471];

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