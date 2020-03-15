function [Y,Xf,Af] = nnet_2020_03_05_18_31(X,~,~)
%NNET_2020_03_05_18_31 neural network simulation function.
%
% Auto-generated by MATLAB, 05-Mar-2020 18:31:29.
% 
% [Y] = nnet_2020_03_05_18_31(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 12xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-2.27996354434799;-2.57484418426004;-0.729448704229556;-0.593846911759953;-0.859839997460717;-1.2622649669569;-1;-1;-1;-1;-1;-1];
x1_step1.gain = [0.63164371502551;0.452882988170205;0.987881852399362;0.957496537241149;1.21564920211708;1.10567597683335;1;1;1;1;1;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.77312781392354223975;-2.4456979261558213778;-2.5079133784110889671;-1.2554216219902312179;1.4050560476815219246;-0.8784492560817519502;0.26972320282281031512;-0.24725323763723938697;0.055943472165048704903;0.032352598604036539309;1.0051163414158807985;0.27231108664719128054;0.3826694869608569527;1.4547365596513104347;-2.0413402714768422719;-0.071997148918296463571;0.075914175821131940625;-0.25973484091134413365;0.52082232159237507041;1.0343618001981866872;-0.3378493834367637394;-0.18031154684799149956;-1.020873329789039996;-0.78937394201716282183;-0.24970606977378043134;-0.0095680082974377855809;2.4109281741416226374;-1.1457864683925156335;1.3468830278223662411;-2.5468620149403928821];
IW1_1 = [0.42039166946558975679 0.16385545010870358751 0.82275180205254228039 -0.27608454996856274244 0.10156324682423083228 -0.0071294721458220657734 1.2160868393670256626 1.6890358696258007676 -0.24238484055533654482 0.30834112348134978809 1.1120048339548920247 -0.2243567737050905786;0.011796841897514682984 0.14965334165859886517 0.089102329380473896392 0.29999056443581534337 -0.31922708305140007345 0.22939100876199866597 0.83461707221076919527 -1.3810560760101444622 -3.210408322264221237 1.5229681258596059479 2.7158894804277049495 0.87807992655009636795;0.18281009166301870339 0.87309146402945314591 -0.26470440544796169924 0.53859218786026519243 -0.54003707309539672465 0.28287680149372362859 0.81390143759481714536 0.31721494564641172786 -0.87962656180366627368 1.3929518270113427914 1.2334819532277336851 0.3694278030959552428;-0.87982063582799596713 0.38571072449866183884 -0.18385712427348968601 0.18998099289606562579 0.084659858104176310967 -0.014322739217512826457 1.6047168945799283879 2.0213676895917251919 0.40488436743769601334 0.03999217658286705751 -2.0232866494904144794 -1.1114730431564678881;1.103912046429321725 0.42698959839127686999 -0.18387717643452383687 -0.11961901048606453368 -0.30616294019654965997 0.032388553219514219783 -0.090040097304024846459 -0.6959252747671390571 -1.0570574474320133529 -1.1847828623022984562 0.37682639546442942002 0.11564819695816841438;0.25475154673449851295 0.063156279041205617975 0.69927579692504993059 -0.67180069180998780443 0.16001280195273703821 -0.036039438348166003001 1.0800253078203534418 2.2961544720154494215 0.69492056717593375659 -0.096099769152315578369 1.037674968373945017 0.29282144303593732326;1.6129927474925322972 0.63313373238889669725 1.0149944943336131065 0.77301877280013298055 0.03128047013111764918 0.039373772625475651765 -1.0099551082973949878 -4.3415774909647728563 -2.6854819596721455355 -0.99978390819197093808 -1.1239757086713406675 -1.2453457626687525206;0.20653870512465224984 -0.27192001890877420811 0.039145558212113885921 -0.23352282150277506312 -0.00092581941975973329993 0.13852111582349341257 0.50918431932947505469 2.3625429400503290545 1.3311116286052460111 -0.39122659528650616245 -1.2789602851335228895 -1.8443021508360912808;0.36070567987643831209 0.05483687895352465369 0.195353093606845718 0.36157552441118512299 0.40499495781024386076 0.11107860039339201297 0.57094478334099107997 1.4416598720329951711 0.73679772051759528928 -0.3068989081548157416 -0.30638220988134357015 -0.48553251388139251343;-0.022658168333638994268 -0.027831297993173281097 -0.047887499899261894154 0.38909618768168729996 0.11820652778416211592 -0.15162504430261297683 -0.39858659853051248145 1.7606692885588479669 2.4997941311036404599 0.45274346577835805627 -0.012718484236598570983 -0.80903064747555308855;-0.64455452472753982374 0.10027132773815269906 -0.51951660162207824456 -0.21684380770045502684 0.42016479169022280526 -0.055479532572284544212 -0.99452868237790326145 -2.806828361044369835 -1.2086499924644509285 -0.38996054726190515494 -0.47664997111986873168 -0.19620581325868668987;0.11648337789459091429 0.079819062097035237136 0.12931038539846789748 0.10947934518549708172 0.24116097739424438573 0.015108616527468830826 0.64815210679901080759 -0.20894065873214373585 -0.32842547628625379241 0.83825077272638459647 1.9125967152003799931 1.7058533585928350451;0.38774379960185850758 -0.11828380996309560691 -0.01669882289835351169 -0.38288448186858708988 -0.37857379516976136236 0.063778443218641275858 -0.24791193201742814378 -1.8740379318119388952 -1.0115240273066561372 0.10927455102897021266 1.5359654157900186622 1.5859611726640321372;-0.53497378255452343598 -0.07444539832031921156 -0.22558936368783469306 -0.033405368776049056134 0.11493878185861523822 -0.026954027476756878962 1.4352822676637266852 1.3850496196562178675 0.060687129775182487024 -0.20686569348855993189 0.64756046331042371733 1.4267208256750598672;0.15657663043035038863 -0.1601097743298347531 0.077914363649933177092 -0.021513257573289246843 -0.15989112549801848262 0.0067537167079153043875 -1.3465478340474725805 0.30279543863132152293 1.9689811019931313663 0.47535912940655489223 2.048286341047864223 1.3639508774328021889;-0.24784450871269114725 0.087573280701895961364 -0.10889538337852831817 -0.15473231608841347962 0.11331398990641386448 0.21708688739468789963 0.69749045294412204132 1.964118375931593885 1.3628326549386196032 0.22028365127187549133 1.2231093022814050197 1.7077915713975919676;0.10428072558541369441 -0.2009342679072110005 -0.16008008495689490802 -0.47743783007061646462 -0.35483317865878338804 0.051231321202934325398 0.96592409657725886163 1.1380693741418999032 0.44964666467873004585 -0.35927558670637399496 0.31591855849475386497 0.63014588001478810675;0.015876807572285468739 0.99691354029129264358 0.30586825822523722174 0.52816795363917612427 0.51789822089950743322 -0.40471384822566408568 0.34100679376126780129 -2.813521421077381568 -0.018084573415010241071 1.1325385975210668921 1.8870969531126746688 1.4769770861710040943;-1.1905071206524295846 0.20900863776045358611 0.26686730035449968135 0.29151843073589628341 0.55387770482527343141 0.021236289725985259963 -0.582609741563399397 -0.13505559378883369437 -0.44887631502692476593 1.3763548749477254596 -0.38528897983932497429 -0.16770652417141351798;0.86795295884423639965 0.49034644541335603574 0.62810813184516267604 -0.11924229470496981731 -0.51803015660222651029 0.25168537486728825847 -1.6057598052047843407 0.73209432472381408274 0.73320026336067256612 -1.8099712759535464546 -2.7959462390089742101 -1.4693263052660063828;-0.72423642004244070947 0.82359760371566848036 -0.67191513108637379581 0.43580245401431061403 0.51039115338799967425 -0.0095376718748429456518 -0.053796156999457576953 -0.30292306997448747996 -0.0053983061477575131071 0.25133106543744820893 0.93161311060890694691 0.44674445344271057623;-0.1287287588136337424 -0.075513394987018422899 0.30435529483984957722 0.3391322990034660867 0.300131194169212967 0.0079507530175516524795 0.99848312861759269676 1.9521871281826455036 -0.1603971864142615944 0.63335925588071884107 0.58941915127610444713 -0.24693908058409513862;0.38562026650855713017 -0.42983749360119771765 -0.028271623914844359188 -0.53075225913991219695 -0.080668063154411856486 -0.10674289437048874174 0.75654778834017577527 -0.80286370372400051032 -0.55171809132710236678 -0.15888866231436180043 -0.21302891707534332699 -0.40923996758287289888;0.055602480433368732526 0.42735547887702773906 0.53128656160520038032 0.36443608318445525063 0.1113860117952559653 0.20696513472943456335 -0.98162019269216616024 -0.8573277367134923832 -1.0175526375285448566 0.25888860433563704033 0.78059365641063871433 1.0762451071119996104;-0.76248517099753609916 0.44191168528561552975 -0.9857914988048807059 0.58574749892225741466 0.75380634762968035112 0.22159523057297136894 0.36579205891511551441 0.65741210444550091108 0.57265473752443607491 -0.79201530567972211472 -0.20111257364271045711 0.1062102066268583167;0.32419772001634022951 -0.1077905672201095455 0.090217897156188664454 0.42602472264917379263 0.19800581253343269017 -0.052754582195805636846 -0.64642911794059476627 2.3706148725302567115 3.5722004672910889767 0.00667669029271112692 0.015953723694849281795 -0.45672861714472451355;-2.7945422562954473911 -0.15401307813351966525 -0.49172109213619502599 0.37845428483306753353 0.48756601028278756527 -0.14514499958839502392 0.53336378117030325274 2.9907622687370150061 0.17504581976854516334 2.1462277160289797706 -1.7320531709075976945 0.14104945971387597425;-0.70034822967892185197 -0.1028343587589442415 -0.11937889098371716468 -0.47427023060960149525 0.10401705447711387387 -0.26447540836510002871 -0.26851016829710894207 -3.3234368110586780531 -3.6170897465532698689 0.80533839468976153864 0.8032216632659251454 0.093123852316413568975;-0.18958651706397769421 0.22619993009885228097 -0.30817745720851702496 0.37288783508695338975 0.568600397380094158 -0.13898126411508859857 1.6192981761631604076 0.92043504315554292017 0.0072950032709942803921 0.2240709054391408217 0.029380460950413581772 0.096113628341031820401;-1.367298404556494873 -2.0056571432681216116 -1.5208717757279468508 -2.2174407768108714833 1.2034032362938298366 0.7142729887600192118 0.089527244251564622179 1.2839173834449448819 -1.496174578992538251 0.16338657058275332656 2.1629055060591482196 2.4177744849872713218];

% Layer 2
b2 = [-0.69054673659989085355;2.9448841488135801825;0.00026137958297961332021];
LW2_1 = [-1.158839649260195781 0.61632857215481462454 -0.80700322271714441946 1.8314278368686027765 1.2932233402309838421 1.5031499963229040251 0.58384517756885190298 1.7568800855738631039 -3.0275121300731702156 0.17880588727886506661 -0.022689893091930320168 0.70027743554924326563 -0.098945153211437208984 -0.45282914189570605545 -0.24209799025620282187 0.97125496359972895988 -1.3627179234240287542 0.24964806264651218104 -1.1725661990684033142 -0.87107614319539605763 -1.6030178590900303615 1.1157055375824098675 -0.91398822305868610005 0.19488313369656146135 1.1024852028028027107 -0.33885647448503530832 0.34532967667066616357 -0.3691305665494141186 1.5958902031064401417 -0.29412725419518925829;-0.037269546560671606983 0.58074349562721261364 -0.79460392452696726462 3.6045809851967280935 0.6739464646620004773 -0.013434813830179678784 0.015625550182986069209 0.79940350182262431833 -1.0096485866604658632 -1.504871603566663163 -1.7141768173885596482 -1.1254747988696176542 2.0296725107253426046 2.2255571626149301956 0.057808508342710465677 1.6433329908335403413 -4.463755924794363672 0.16716846800345161461 -0.51860428178180717484 0.49674004604242555372 1.5820718261280577455 1.4165988182205959678 2.1586241080952035531 -1.5822243661080159782 -0.5934713223659312531 1.2627940884304720282 0.12446944628301223612 -0.62684172521725189942 -1.2207243536986238652 -0.17693250061482315716;2.441187130271576855 -3.6711006942658324093 3.1512433563111605928 -0.69726548578258285005 -0.57592283602970428813 -1.7785010331440898046 -0.95149705927238392711 2.1101537644233574476 -1.7951409376106350901 -4.6634897081813067032 0.71199323143532555758 2.2122724106974813729 -2.4043750262797178863 1.7771086062899421609 5.1118925729599276053 -3.1994394275284920859 0.67471330515394978278 0.23889600375550429368 0.016273863981428691311 -0.38335535761454997417 0.25220662832589346714 -1.4353914338789524052 -1.6782509037956996334 0.60472953151322195176 -0.36903263631793170596 2.9815123077425247367 -0.029120270963761882982 -2.1760304788831899003 1.2872883717769165735 0.11880680817555057394];

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