function [Y,Xf,Af] = nnet_2020_03_12_16_01(X,~,~)
%NNET_2020_03_12_16_01 neural network simulation function.
%
% Auto-generated by MATLAB, 12-Mar-2020 16:01:48.
% 
% [Y] = nnet_2020_03_12_16_01(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-4.66722932348959;-1.87069250390506;-2.76693312377929;-1.96463294853913;-1.50418417126776;-1.51821272217858;-1.24104321265912;-0.943822029337829;-1.38976317019856;-1;-1;-1;-1;-1;-1;-1;-1];
x1_step1.gain = [0.247209017235578;0.258712630681784;0.444337646118018;0.700193484405396;0.422911992480351;0.654351841141191;0.692711071625961;1.18765044002801;0.879204187193534;1;1;1;1;1;1;1;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.1649320303642842767;2.1058913313592362293;1.233964803546707012;1.4606510246789981;0.99589578115162913363;-0.5916232117611087693;1.1858682259709276696;1.1025587030842753933;-1.1547096900743747394;-0.054214619911355892146;-0.76916794891117856725;0.99420857772933746155;0.17943154594433130544;-0.50539962807873828332;0.49427238291465658282;0.051641877202113932455;-0.19080695619170623401;-0.50345299334151538684;0.96525979240499282241;0.35433004480494989519;-0.06095447038011797547;0.65233493184072322002;1.1515833040670049314;-1.2356150197541209135;1.4693949874405660694;1.0497707747197835282;1.3932783292858073043;0.99756007880351071915;-1.779912287849030994;-1.9099626797440916182];
IW1_1 = [-0.70375171594143448583 -0.45759795785673867963 -0.28753664108383214515 -0.10701000085569599418 0.24082842889620564852 0.76424711684588120253 0.013526870129761656147 -0.29549480551600065814 0.15975278877950885303 -0.34313590873150107141 -0.87898630529763677188 0.067471031276407938959 0.088773614193786876125 -0.72469386963042914385 -0.28851281530427741462 0.15810281380877144741 -0.064344320003008367848;-0.95254406495229571039 0.59183118493255915116 0.10654123616041914768 0.84311686328703216997 0.37917513599833374549 0.18871778189902360667 0.093658175961130338272 0.18826025479536748208 0.051583097302299069087 1.0989080657124226104 -0.47648427629876133338 -1.2462294079620965803 -0.56125940805212282125 -0.32964391730139536163 -0.62662724756918330993 0.60025150748743294482 0.0061822265634078384794;-1.0120437998798574863 -0.55788430381389442214 0.73202235864511211449 0.28298593634011909304 -0.078716067105847262186 0.57624660048739595464 0.012625563691037748962 -0.51006313883475096205 0.039398039262342275657 0.40773153844212495223 -0.38677902112388279354 -0.26873877139591345609 -0.22943119212985482736 0.70118399693641952908 -0.29285973380795288223 0.43302278098604474188 -0.45678813874970275455;-0.26360985050706753352 -0.031890494442224914806 -0.14690787261503132632 -0.08040556085969322031 0.50981749626256867014 -0.33705534240463252971 -0.78307176494474606265 0.23065160106874499357 -0.22621686902315410461 -0.3658362435680563296 -0.46891479478931546998 -0.011093542133958488038 0.25979871354007694206 0.37887729801108771532 0.17892686382339484386 0.64005544164672578678 0.53152797228465797907;0.49596226499989515268 -0.41917464107122637262 0.33784611056009339336 -0.30329690312211704528 -0.37896440317485402538 -0.54951237837330446645 -0.14629227175595235777 0.25680996954531148946 -0.014965394259362959847 0.75055559641435731777 -0.092045301745402541282 0.53769280588362944329 -0.35507343697332471866 0.8436545089783530349 0.60923254169536811276 -0.29007745691412634148 -0.42675523872092901989;-0.46394958439957384444 -0.53967845809404091373 -0.13766866609248021303 -0.37780626202336803354 0.77888581876123941239 0.24195412605356245161 -0.35566396066611222215 0.25335049542457688254 0.038149397171378704452 -1.2778152387265799828 0.31478366700563403802 -0.057427865142827896894 -0.39696682150907364184 -0.42487812105803673957 0.54083679452168742507 0.2746673606025344383 0.087112395364384614949;-0.14092160220846028373 1.3996544389950602572 -0.62182419882700334934 0.40678836866936846084 -0.055042437245004453839 0.603829840575166088 -0.7762062764654465985 -0.14397578974494765203 0.16440891363035500805 -0.017386893520042051575 -0.40323652493916356532 0.57248894216046564143 -0.53587693453479134931 0.19479578293484398066 0.84200896842960082367 0.16338554109206548937 0.12945417286676610957;-0.92116182673936120029 -0.48640294372650866306 -0.5348644997760779729 -0.078618298433185726526 0.13817812143486629672 -0.35235893731585038458 -0.74136936653856988499 0.20067206166725845562 0.24814587471333235391 -0.88243845486671268485 0.30092763025065633054 -0.49533084945835442747 0.41956551672089181793 0.026120482261406920838 -0.59612392291039850711 0.027838562086114706368 -0.28432201089121572624;-0.21920545976408170596 0.60662372360375726998 0.4813199845035833202 0.32597385370945503036 0.35578700717181388491 0.4834697665455531812 0.056967669004893778384 -0.39915702503733729856 0.20468344384451125917 -0.21619734649413127614 -0.2713680926340510946 -0.74832915224185936331 0.24158725327208743772 -0.30878251256764660537 -0.25679075095104492155 -0.34328895683540261441 -0.64466286876999034927;-0.98187396088077538714 0.6801537955076362163 0.12536977426700399163 -0.087973675987799829223 0.71203736660375693823 0.47552809203011164874 -0.10141831788008175141 0.054805860915453474824 0.46157062412311467003 -0.034995504864067862594 -0.51054777024887487702 -0.6166007625402553094 0.16247953011930604816 -0.64669112905642134326 -0.15146053024260555664 -0.44278504188895578819 -0.29199442800960040412;0.47422022088849730137 1.1540315580804747508 0.17622280604589643516 0.044526982073708996912 -0.17490635036659676893 -0.77512092074223182259 0.39546001357554877442 0.12720781308946041421 -0.17672013491717022959 0.058910374289615617671 -0.13573654120007638824 -0.31541867432189241161 1.2703148139311084552 -0.55009576910828639207 -0.13668589890210300353 0.07955359320161535519 -0.037620701326944769316;-0.20198993603221929161 0.20299606732641942153 0.010324576075384986687 0.11250979719147798097 -0.036304326360099462012 1.0006130239875750121 0.036568026817350260393 -0.32489825648453535267 0.2097547492321381879 -0.30457291847578044086 -0.17680989034264754656 0.42391175955278836973 -0.11442761090832405391 0.35121835323456834788 0.01582224478530067785 0.96969864002519012836 -0.062677217414725233668;-0.57275443797898062392 0.40693798999706437547 0.34343731231072172072 -0.43305846944713699553 -0.43779809934883656686 -0.075975542548876812132 -0.086773679886534171857 0.05630719698024019837 -0.12625668889070199019 0.025077091008968866587 -0.2347153250777728184 -0.81239494154107527724 0.18380727792648418606 -0.68386456961912012176 -0.05153079367859045995 -0.37919329092712938944 -0.5990111135147418997;-0.77588608869129827017 -0.52659172101581608327 0.032953930868546990607 0.071482667836242325121 0.49918782262895083912 0.19943832994200952347 -0.64725027894783671023 0.20564561854607746327 -0.047572890976011646125 -1.1005896921596449722 0.2742441554509180146 0.30847954404624589042 -0.54621101473745026222 -0.24980787240412472006 0.028404723976054278844 0.18543168465854215077 0.22607514364735686874;1.1956639056623907802 -0.055093506395912157647 -0.8557672553563304696 -0.16473108520491180906 -0.46010800864386230957 0.064942766543200761142 -0.0034361779890249646852 -0.28632970608920876865 -0.343866785177013079 -0.24149587837269789747 0.032419714943151260733 0.12669764926628737656 0.27983529510013566144 -0.30517011205552091457 0.31130272931550162507 -0.616386791565036396 -0.59792219785233824414;0.46765154455052998506 -0.27393658811048954282 -0.37551406119696095764 -0.57850143630850059395 -0.15305737736132338234 0.33479508938391805417 0.23536512403932341453 0.070678728205947394247 -0.0023183048904501539142 0.058995452836667947039 0.050460653715478920311 0.3612972390582132598 -0.43729725736578506368 0.34259829285691950629 0.10920581680986843343 0.025040081708513572956 -0.27698655640398367472;-0.096962910572962926481 -0.413213371250519379 -0.0014177103044049816673 -0.094439757632488421923 0.2087508093112062646 -0.1146236062413325818 -0.0049262172604837631273 0.034520906252022524885 0.041024899443628561058 0.44024525136320930008 0.16173717206536264435 -0.29062730078587867277 -0.64318567774675883264 0.28537398768728111653 0.20235232989121579306 0.52039164591035669272 0.096548398726398934611;-0.028378801317363876872 -0.44636251222451794884 0.13665140993620178422 -0.27689410704458200518 0.041422251994515817131 -0.082676781746876440438 0.21391167089568791093 0.10252617997517476933 -0.042127768268460971668 -0.02281385980893010923 0.045636536957344986387 -0.88586631029014262229 0.080031580045963515713 0.052926130948711874735 0.021229098333605631033 -0.4243907072227670052 -0.42038704714059638157;-0.4742113367675517277 -0.43179611684785362824 -0.2410312880398237656 0.35574028211492114693 0.099057442560962119527 0.62311499680226145603 0.54357494764929892295 0.25552350869661866195 0.30579663828910425005 0.0059401091558845941359 -0.36896888054275539526 -0.30526173184457605636 0.082729237339686956254 -0.85105172677612672683 -0.25989667745888572759 -0.26406487021976421392 0.31315687523201002351;0.6539406745249283226 -1.4497884078390372853 -0.35284339064830105848 0.11031776140811304365 0.31241441280627857413 -0.32255659862388508241 0.40205968738082603631 0.39417010971825955368 -0.20207307179233333083 0.16451462927414159387 0.46227251042568523465 0.15925075692798174032 -0.71584169781342210648 0.38464882852225895649 -0.62554420826773970976 0.20769217987803528791 0.33579425769303244964;-0.54817521101935362537 -0.76134471980135165481 0.38123485640813215003 0.35192190343082091486 -0.46075846390261365748 -0.15107444658607779853 -0.15381569900984010113 -0.092087713859998407417 -0.13535532838387237531 -0.058688945068143559547 -0.22897726125911715966 -0.081860628425018483467 1.0063457508446094213 -0.25177595035880839802 -0.42086494596152140835 0.23010110569498196109 0.24291450949371021228;1.2056760376066950613 0.31524866611546309425 -0.714272897005509666 0.23185033925365361451 0.60507316604409500105 -0.73091956730295515321 -0.52050549186796501999 0.63187893255112237423 -0.19165481792911659431 -0.45796822445961937742 -0.01898942634158877138 -0.82027394273988629703 0.78843777381528257209 -0.33938476894755698066 0.20651867384578578846 -0.27805910622376200925 -0.89413089368723774886;0.4059989134984076653 0.40616607596143255998 -0.50529725436484218548 0.046670766040954055065 -0.24886373655005505912 -0.31535465959686220794 0.40166656605694650573 0.11379412735548019575 -0.28765663981254174786 0.23710067778825463791 0.53029365569497788346 -0.31262666704944552798 -0.074030229446855846942 -0.284942237789375441 0.071318714227592755472 -0.030211307558307767818 0.66307927884516959516;-0.91720525646875816772 -0.14961220465823399373 -0.13754972870331677592 0.58185058347708018545 0.38664722149638847126 -0.076100459378731569182 -0.21955606212045614134 0.10805486947768919159 -0.16398095320388905716 0.12727425760731392312 -0.55267484753847995194 -0.16384817511622437602 -0.20780139546384063243 -0.65866432895944904136 -0.082748463596745827631 0.33060012694818846635 -0.11008457395762939746;0.978413382478664162 0.24282812158790473278 0.76920842554999757645 -0.61812628038855788049 -0.31580925434026457799 -0.061420078723563101275 0.081183019736746991901 -0.075455060301051105065 -0.24340203053305636827 -0.13734962566741293344 -0.7690021676450065069 -0.43222985669738378522 0.0634172297658177897 -0.16923856364172962241 -0.034087238641493711189 0.032946388878245101706 0.037257345749280502967;0.37434992418557438976 0.083463935542994929784 -0.16007608202506667938 0.13951718402746388081 -0.64975229533870282328 0.16718157618956888677 0.090241080674105508819 -0.16636453963303876602 -0.19045168846769652826 -0.37254820370739522017 -0.5094367756053713725 -0.11665649696269374258 0.64276200729019716817 0.47681388644710887004 -0.25535944551815942249 0.99098960293560822699 0.8503104660919071911;0.73307970345778505905 -0.50904967904765430209 -0.039755571519340268671 0.70568400652906271286 0.2885909159845196359 0.29240703235205328658 -0.31916181373417201739 0.17321949955380028241 0.27721586194528746727 0.38472830685164366438 0.0021602274183693198975 -0.11620060852182483024 0.12730227409596236998 0.1543248980067223608 -0.013967212545670676729 -0.28858407392029022942 -0.027756725553026405373;0.4278209390600236639 -0.6061664688050512062 0.16186123458028450761 0.45967292471793419484 0.10713833583202080069 0.97228987438185576764 -0.39028484694128684485 -0.62360823830955858238 0.52221658484654709387 -0.48225400524285166881 -0.083347974102069885549 0.3087786797032524011 0.12785857689127705372 0.81117359064798566504 -0.2785774166216854475 1.5617796542342763377 -0.44178037330903829893;-0.73392428950361798812 0.082674918453657431083 -0.41514378752266350991 0.18241202948493859259 0.55738722305112553546 0.22219698981991925502 -0.16868970396091353958 -0.13443662176450113521 0.47553777016462667193 -0.38895361529520577903 0.57649575831684796157 0.69083862942411045083 -1.0626511831276188591 0.75244329717289915038 -0.45885180633194305733 -0.3331320377322542492 1.1145816803233545844;-0.0032585725422322095977 0.08963712617150305928 -0.47810181262465278884 0.52121609963921244901 -0.14275072866602356303 -0.57481269637043119758 -0.17715484440373632791 -0.14672866682567359264 -0.31661282676002244774 -0.71567619340840316067 -0.20039536480742797808 -0.81811429744146568765 0.26258553588090222108 -0.38358749087479043771 0.16388951572564949832 0.58847989127203303994 0.037466241230018924868];

% Layer 2
b2 = [-0.18455249558913677799;1.3169334752277590539;-0.32953020818421846494];
LW2_1 = [-0.74578571777848601254 0.033659568046085337756 0.0572551381731746592 0.44237952508052252876 -0.94544223676864580508 -1.0635352204001040644 -0.3152703735965085019 0.62050589011557666375 -0.17745363867601796071 0.22612018219867330293 0.3024251786571047762 -0.44285807015244316176 0.023399243514747890271 0.41563819114492256412 0.12684797513006035308 0.80615906792904201694 0.21074858936895501937 -1.394615524454467792 0.5511596267823769546 -0.034173479095667647076 0.81519159716814693084 -0.08571728532269985712 -0.0041932533772147898704 0.46030113798616173471 0.92447396152647942635 -0.55944893150227248402 -0.051659373838190997263 0.10475617898222311375 0.41506132337677614785 -0.34830078187669472056;1.4601965570614321432 -0.80448483805398096624 -0.28313210296013680001 -0.25516905415355023434 0.64484744537208782411 0.30773219983051502968 0.13058993039327948527 -0.75809644924332098004 0.88972850450581775217 0.7101775289357248333 0.76722473606259877776 0.833690220212823907 0.33795541846676457887 0.66535070734020484995 0.17157155248751657095 -0.20153819304193631989 -0.39973423810191444083 -1.0339081643854453052 -0.7685016374611690626 0.83387918256551241125 0.73038511206345535509 -0.24231119255190067263 -0.4328286787724610063 -0.94966670277971643177 -0.76977610330475310096 -0.26112567572625355883 0.1148539869993119994 -0.41701462899142160978 0.40048272522781275251 0.29295109311546657649;-0.10668289572191876069 0.83838855030909553712 -1.3657158493268097832 -0.30233982279900811774 0.57057998354752037518 0.44979708528863976413 -1.3148711850604277718 -0.52464294483599260843 0.19467004765888942797 -0.32173600384126382901 -1.4119092316800505582 1.9169381588854070841 0.89235994920809835751 -0.61784623745489120061 0.88173203663183863377 -1.037654676133353604 2.4303803981176841376 -1.4752659341200096144 0.1517736797412552463 -1.0977420493436897964 -0.28089299201197087674 0.55575631077365073018 -0.61539063637255642103 -0.093618379132615170413 -0.5153633728974974515 -0.11134526217779747159 1.0388169103632838297 -0.81942939178792129074 -0.87904812850800251312 0.85978364454207067968];

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
