function [Y,Xf,Af] = nnet_2020_03_28_17_11(X,~,~)
%NNET_2020_03_28_17_11 neural network simulation function.
%
% Auto-generated by MATLAB, 28-Mar-2020 17:11:29.
% 
% [Y] = nnet_2020_03_28_17_11(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 11xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 2xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-1.8712794308374;-0.657253600993518;-1.14536857211879;-0.621477272756634;-0.60395073081581;-1;-1;-1;-1;-1;-1];
x1_step1.gain = [0.797951306493166;0.698510287303379;0.863112388626723;1.31326373837786;1.54290319083816;1;1;1;1;1;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.0837862166306828637;-2.9638487218609634866;-1.4921659552596981957;-1.1080802234450428667;3.0914578489244046722;1.5950111901502306289;-1.1016406664240776525;-1.357567106240459287;-2.4618214787073511296;2.0655691601954755399;-0.11282933856075887114;-0.039505985256699363561;1.5123158974733241777;-0.66775403953665590784;-0.5575997051427636908;1.342827932519495393;-0.85750461815714251745;-0.93042732812237216589;-0.22560078381721326846;0.87298543503372239005;-1.7832179304096118511;0.84655281457713948967;2.1589844830593007963;-0.1828036363569434386;-0.54847058424833106116;-1.0838783186045468643;-3.1046043692423310745;-2.6612250554123932922;1.8807242371696109018;2.2384143357341810265];
IW1_1 = [0.49934491315335199912 -0.36691501895417710788 0.19790073658014151192 -0.6537807981608823038 -0.11647197345891248388 0.057310084673772569708 0.48044145754459732789 0.16347313828831472327 0.52883922720358156333 -1.2394668289790125026 0.30146413439049085659;0.56522557645802229409 1.1619642646730217184 0.9837009797148704715 0.075228503068367050255 -0.084323984932720502661 -1.0971876261368753713 1.4036681479805330763 -0.5614517153087714707 0.16185696620048109495 0.93850983937034193616 0.17185101008840966696;0.45462504009701149021 0.27512815661310335846 0.24212076585249081573 -0.56488231312700731834 0.11070044436734026438 1.2611920578625006595 0.16480041018551003229 0.86093916019440031029 0.12567190123133470392 -0.36584340221612937638 0.39565909548044680966;-0.020970618017154132179 -0.38955256377597008433 -0.70399874510012538931 0.30351351396729087107 -0.1423851446463705317 -1.0045548111693840543 0.098191025426928910114 0.072372978365046178317 0.072938891528023311328 -0.38550909953455547718 -0.94706382508663500275;-1.7461354388634462431 0.66292214024814577478 -0.68985811139667962966 0.022797517049480221785 0.042314757458627609565 -0.28620344278641396407 -1.0213800392849596843 0.55088618131045197224 -0.27024030242683766412 0.28486031261085303923 0.63367155146075915262;0.76963214106102917356 -0.2653206873008854827 0.74344405504823929398 -0.092337208396009742839 -0.2903110275617597491 1.2433369756295784558 -0.32620571431244471672 0.40311563756839752237 0.14039715313922082607 -0.65695707405950698465 -0.42621209325888759212;-0.046488746654167995009 0.28296074713744995277 0.71802085823196648562 -0.40111726840697414787 -0.28885761693749029133 -0.62122879480857373657 -0.2801630366448109033 0.31106564874391817721 0.31221178722777043202 0.66041590078868217528 0.88192932810049473602;0.31963721149752849504 -0.21905953768343303101 0.58134803713007787707 0.018840365556937500868 0.23133946149532663705 -0.27281418530494899732 -0.011387018907904684212 -0.13882964983034110085 1.338450407160435951 -0.13187617033941392841 -0.12460855724731104144;-0.0353615999977405851 -1.4272059229709606054 -1.3737150909912041463 -0.14610101208136458406 0.43743585348119762291 -1.5548420876333934881 0.75570801903131012978 0.028946134559620237403 -0.68999101019173136873 -0.26850334941882431083 0.97624396229953280812;0.062665421992026931752 0.58014234924799479298 -0.43579096378601006601 -0.18745163116201493669 0.12831512126393146422 -1.0299125415401109152 -0.28483223414484415414 -0.50134719799546934027 -0.48887496535713265144 0.46479391119945240307 0.20306915042549139594;0.37933565417378128926 -0.3959372788405426391 0.56630213921313432568 0.14406032940211688231 -0.096426233939352679059 0.35095289928852396732 -0.030682680640097732133 -0.81255139236729601038 0.23075158397850950576 0.24168050808778743233 0.83680386945913176611;1.4949056553512907186 0.85406637771369142698 1.9782791937124546422 -0.80744898239576434751 0.10097091275087476114 1.4579070887995397854 -1.0772643769033836136 -0.067400972727103172777 0.0013667764909157855552 -0.30825221214541737025 0.60924412693360907589;-1.091016417980200881 1.3582837114840622394 -0.12892664686155752563 -0.92176160986803212793 -0.007772651521104925916 -0.84941034122076952695 -0.11834995866243841478 1.4727244271926402952 0.045036679234393875437 1.0067263307540399442 -0.23900697193262201989;-0.14565863743175591627 -0.66538383143947321674 -0.18734744486220639148 -0.92181284899366777719 1.5833981275130435407 0.65143569910201326589 -0.54770023473340601061 0.64393493855593930686 -0.010355467764220872492 1.2546225949970040148 0.6878657111682001668;0.65387118441218849973 -0.91688565595243121109 -0.99611652832692076931 0.0052057885685833022241 0.43976483299854612907 -0.56041977073643489682 -0.33527868346695682833 0.094428478168185198127 0.10091766398213328215 0.23436570348964674837 -0.17441480530964623141;-0.26116401608632766607 0.49307921957478423325 -0.40854354104478468601 0.67427720579345085117 -0.23195454677306634461 -1.6895718643159918226 -0.46605363371199093336 0.26602491965931374107 -0.60842933781794106274 0.17344255806772340689 0.075956539145362336751;0.30253731046638904223 -0.6431995322595004394 0.33725112920698502306 0.19685365713808605781 0.050369685419083376487 -0.13126755635577644354 0.32888516711909254742 -1.1332465036711691919 -0.11468301180397112715 0.56642638542136436453 0.39148232563435547693;0.16028518098466812414 0.86447109257605281307 0.52079656365080329028 -0.79678892138558854441 -1.4413923522291804957 -0.0979420467585607607 0.00260064816082489908 -0.78345444205580205832 0.11368983510601457509 -0.48528510748362646554 0.77167649288100226279;-0.58165552841570511422 0.37299254018053534798 -0.54036102181141121292 -0.068228303052692354513 0.35553985536854165739 0.13805625578346922278 0.086507519447859962392 0.25638640669520090798 0.50754802075219940249 0.77494408751367116484 -0.4727905287761444697;-0.35746329350557237703 0.47607114472096789815 0.20658263376537316924 -0.51263936006100851994 -0.51050756924781870438 -0.76820214588380231202 -0.34560682975179990439 0.73662688635657669334 0.56658184941868894313 1.1452040761295811944 -0.76549578031786080778;-0.81495070186162443804 -0.77405993075444989859 -0.98675882417980942396 -0.29768338174191011358 0.82402202434283600141 0.69614135909323560103 -0.1511395944369605393 1.2632713306024598854 0.79472807571453940589 0.45201012857996630823 1.0725816300458019814;0.050659237992211510193 -1.1468438678489452087 -0.89874913193877803419 0.85214051400353540977 1.1740200675591334623 0.24785947467633584473 -0.029804860478162713278 0.24159377345935303949 -0.1396128570873085839 0.14700043381480953908 -1.2800051257230977697;1.9315951240223796503 0.96187859260663799432 -0.28293385837997159227 -0.0080701437766438812588 0.60692531342607747291 -0.47163026335268454226 -1.9753630249323399504 -0.78089486246325656715 -0.8283039114177013218 0.093353070912685492955 0.8403964401768060366;-0.059032852640153607449 -1.1677443653331884832 0.25406528268418804295 0.20589478683552167593 0.3444963624125665369 0.23230674062716874517 0.033513858482336573996 0.11777590674113104507 -0.16233237953360971084 -0.15873701008302421034 0.50763190447940786942;0.47237502505465528424 -0.44246803726066569196 0.7137991789078529381 -0.35537751715684823273 0.18929844770930923104 0.71727003807383560297 0.013505917808277541348 0.50199630814271900547 -0.12353017175984835896 0.53166336838560257494 -0.17664822300217644657;-0.59141697159552031327 0.58397792620057986657 -0.58184266742557566587 -0.39792111344083541669 -0.043093519488639160653 -0.33030909553301107673 0.31240896635084941524 0.042555873261774831984 -1.0497061844434008027 1.1488542453544889899 0.47221106619632952084;0.40996914098256442838 -1.3031446102237367324 -1.5405585584859831094 -0.22182841652837431701 0.034549080368430712396 -1.6109565463747306246 0.18456870973646141709 -0.95755928784493704597 1.0720398003149045341 0.43659471218141282511 -0.27919065874004156402;-1.3687215589039183339 0.067504725650464378339 0.41398504052486440052 -0.023388955202338034683 -0.76508698499831906403 0.94762825129635375188 1.2425970497895639788 0.57824823764221433553 1.1285516508656368018 0.14006251097174021725 -0.048798039959611447247;0.27117203772224945935 -0.93189732677828851504 -0.11016481440487697896 -0.36429746008864910545 0.79833255747339793018 -1.1282845657655358185 -1.3893988480189953805 -0.47588101996067899702 0.30169277264543625794 -0.67803187643812412588 0.5005706018427125148;-0.18390516668247289256 0.14953638738606936376 1.0409357701034773758 -0.95083252808197449291 0.28610637179549897047 1.6604475859395071424 -0.48785164561398536298 1.2302460876483505281 0.098854417826363769062 -0.55785714842544853465 0.83218572650213618047];

% Layer 2
b2 = [1.0197664126899270709;-0.16462081698140035302];
LW2_1 = [0.031025667631827812121 0.23048064433051140831 -1.1559178946015848943 1.2769480742104444282 0.12638999800077296842 -0.48810953736236967204 -0.26729647910053627724 1.1255745003420782879 -1.1744042298183237083 0.37087360787264656015 0.8308740166377909242 0.05374853685361923733 -0.22013109823681742405 -0.078111787782649030887 -0.19667129129278146382 -0.6179991918738338974 -1.6441286222675164375 0.17695040070167480928 -0.19481221424475969606 0.47665872396000885658 0.090437113481980205276 0.24079364188535851143 0.46177833117113231687 0.97936536642229365945 0.38382890503773303692 1.0793587806075017888 1.0714950484427800959 0.60396411843404929076 -0.50600169670896655294 0.46870148379713150177;0.23147083750879599018 -1.0358323913739255762 1.173296650372122274 -0.96505922221869255839 -0.40805134597816333786 1.6642434931176228741 1.1469921125261539618 -0.6218080469740118188 1.2970601266759893555 -1.549760688100011663 -0.98504578691163402482 -0.16659132848141292027 0.51637066855904223406 0.27964280370845268076 0.71581274177934572922 0.8659844314272550081 1.1076871027383139445 -0.64892656316293928498 1.4398610348505829215 -1.0687841045566541798 -0.40357764645739735432 -0.77478026868442018138 -0.43786676438859989391 -0.42799771336677672195 -0.76329642325220214971 -0.53871994430405800713 -1.4043931625049743683 -1.0338000287258297671 0.17534775762164792057 -0.91775183820186534422];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = [1;1];
y1_step1.xoffset = [-1;-1];

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
