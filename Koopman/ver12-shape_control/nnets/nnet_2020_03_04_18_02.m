function [Y,Xf,Af] = nnet_2020_03_04_18_02(X,~,~)
%NNET_2020_03_04_18_02 neural network simulation function.
%
% Auto-generated by MATLAB, 04-Mar-2020 18:02:01.
% 
% [Y] = nnet_2020_03_04_18_02(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 7xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 2xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.584137318713385;-1;-1;-1;-1;-1;-1];
x1_step1.gain = [1.42467250875656;1;1;1;1;1;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.3479000219907635305;-2.3200796144749338268;1.9166209320235818581;-2.4725969854309286511;1.9113417289542313782;1.9792112492970106441;1.8984950710086980763;-0.47892268398033394083;-0.9765034707762422217;-0.75021281295592656679;0.049660659806735327204;0.65451422883207210912;-1.3744329081475581411;0.062845310191560213586;0.45910662081637765297;-0.15973950418238128113;1.0770507598009082528;-0.59012500528694589175;0.97474189509518827101;0.68654172304016525263;1.0973964096359554166;0.9381389839330485092;1.3389700381498541581;1.0105605563500557054;1.8169064802951315762;2.3097141690693394622;2.0840012444666498048;-1.8555308160166261722;1.9606392971190245156;2.3087803957896997353];
IW1_1 = [-1.1998787513963764528 1.1859131628800956282 -0.71119036303323501791 -0.61516319692162091393 -0.3129570294968053723 0.70421927757534175107 0.50572717389306165447;1.0485539523103972748 -0.31655182233655931867 -0.26319423284578569167 -1.0297743054630272841 -1.1614402990732528576 1.0202436505294307789 -0.04423633624726507757;-0.11585511478127913376 0.73693197956479017385 -0.41330616019478022549 -0.76581340796490682266 1.2984565044941389633 -1.1243063505369765842 0.78022815878878082341;0.72400048508796654989 0.61261015961037657007 -0.38010305603440724953 -0.92729544653185458358 -0.33077254883028378574 -0.994812378961460686 -1.2012086751537816021;-0.897715009857942281 0.33479145771236773488 -0.55342646284890795183 -0.81842507298061983878 -0.50139973378495239498 0.64427110545176435874 -1.228067911890305286;-0.72705896917510626221 0.68862964911531077128 0.4939366600048483158 1.1616639357318072623 0.27936769282117901447 -0.30376280143365369524 0.75152319208037021703;-0.96700774754944718303 -0.27441453511585961467 -0.6580812724928511992 0.35856366610220991031 0.60914693612927206523 -0.97952914430267090662 1.0485103992802544681;1.7086061372385443757 -1.6468217607871606933 0.54456609912258102657 0.85442431994953604857 -0.41511104132337134365 0.19807928167954685916 -0.2411838012318807678;0.56013735157151034461 -0.297589817004446755 -0.81466277363698902381 1.1820127499291392947 -0.21597377778013154837 -1.2163094688002242805 -0.11136889807630703797;0.047434118730901635164 -0.160356847856332696 0.15673670840972078744 -0.81097035524132921047 -1.330224596309333096 -0.2084545792243353568 -0.89254147091266533565;-2.2753302402754602696 0.31707307020095759631 0.057702315406782625495 1.4558121269219439498 0.72449986816689493363 -0.48703130401743627287 -0.40012826761776937312;0.4889911128888731473 0.53841270208208746872 -0.4395882316852587901 0.10670881368765831876 -0.01032123226968110867 -0.63106377997692375548 0.47572410897135275487;2.4961276843475070919 -0.28199732064783711305 0.64552685414896127103 -0.25175438320766052591 -0.26345294327361529207 -0.60558746032415478044 -1.3313597412962239197;-2.4242662246639574697 1.1839894725756310301 -0.13209574703105780857 0.63345578969668592251 0.43030026449109420872 -0.11276701514210646371 -0.55449454346503945779;0.95225555218285495851 -0.67751398170212084882 -0.56976847905378968928 1.2504163865523523214 0.14892873111936041908 0.29249278353913582551 1.0220439697567931248;-1.2808067018904498724 0.0799636011829139709 0.097282847383159135046 -1.1296254480971801293 0.73346028832887277815 1.1953831687897626956 0.1939935310956802661;0.87879589560998927489 -0.38887652473219802562 -0.13018601614180527748 1.3928500545850475145 0.66168922403789987019 0.22650067971171411463 0.43580542919532683221;-1.1246881720930923532 -0.088317738765331221806 0.34983745323509929381 0.78659720416963563316 -0.88999634985125497533 -1.170176461842632909 0.23819423483504348349;0.54730903220772275652 -0.072669532723319371637 -0.94748568113884901187 0.85238600656727747662 0.75577185505787702891 0.906553761687804327 0.4070390091354527673;0.9061965784630575671 -0.22183783516133365188 -0.23841914706033870419 0.97797780333901074279 1.127864045979281471 0.68025196554384992353 0.8618547262833886613;1.5091503252900138055 -0.74721177934699301026 -0.93687670573563319731 -0.056229761012078577354 0.12844540747677543613 -1.3173874644819743018 0.070986822063910812175;0.80821084412400945318 -1.4518015510820343472 -0.14506555625450515135 -0.25867782100934583189 0.78979811243855368996 -0.24938419451739896049 1.4713806295831259874;1.9114887227668124758 -0.84813337707212521632 0.92855169680398663701 0.49629806200429643637 -0.76572885758512954091 -0.1738600984516739667 -0.34569760390377413106;0.35394154487857321101 -0.39530718331045000546 -0.74368583074507510489 -1.1576583347220446019 0.99310699444009065129 0.70154062276681317112 -1.2040867360835847411;-0.44907452904149169903 -1.6665925043030909158 -0.13400667663832929732 0.5151950843105660427 -0.22445895406096663671 -1.0791303388726951606 1.5285206330443428246;1.8643559019552713441 -0.33606312461964205074 -0.35693084569653782401 -1.0890864794969323803 0.71692699744043597576 0.52037083492316182909 -1.1301087109726375424;1.4576826767585824296 -0.073154898967446924107 -0.33763406734461282044 0.52749560049559762742 0.84650978205110105623 0.56583657518390884089 0.64850381572504267602;-1.4304154906946866532 0.56244217363656268116 -0.29746351674494730943 1.0124729102919070023 -0.8860405640852698772 -1.1266849041208473192 -0.026819788204254076114;0.1786538874004485189 0.79548660348078803484 0.41871411635611360813 -0.13063277531707320755 1.2843005651062224803 0.27025846830597938242 -0.97344698641929749972;0.88809053839374119121 -1.5894355404897899664 -0.30361993383284702208 1.0253376166771639433 0.8423565645551132075 0.6210966571722644769 -0.89379179863620195068];

% Layer 2
b2 = [-0.022181362850700370865;-1.149188765242758814];
LW2_1 = [-1.1869933649541954956 0.40259826383190011612 0.7010108943460012787 -0.3256006881312787371 0.072481985335741822007 0.70342432958532119525 -0.73874656857431364454 0.40919840622126912866 0.22880126289905580728 0.020190013830941672818 -0.18064276750747670919 0.56792906861628089654 -0.45805043566949865674 0.40707820870183830664 -0.74667704679048518823 -0.31514474993214808274 -0.6343527498241990914 -0.77539293925869368085 0.46272523203843302086 0.32680775699239061893 0.12540810151020923002 -0.34125979040144888144 0.29703438443368002098 -0.18820473732314851878 0.58729530471090662935 -0.45049732604382325407 0.43405770195505977149 0.9953244331621273755 0.50779076954954283085 0.65371734559247918206;0.39617332637492669312 0.58465856291604434425 0.2786548875874603759 0.80183015936074464314 -0.45948642832874514452 0.3329149502303300201 -0.63182931612802173404 0.053124522492258866024 -0.97218350722639645944 0.41444846135510376994 -1.3685213222635759234 -0.04743043879549575792 -0.55219835128273808689 1.1439902644513983976 -0.41584689944757247115 0.5773013514485501041 0.61840913844453559989 -0.44229922856691916699 0.23051170966746137303 -0.15612556956770379246 -0.23671081450847150651 0.32163489049175220114 -0.2842827346767252128 1.2698481706293123228 -0.27874312188158462122 1.4960446251673098494 0.81218236622904560473 0.57133291432013511013 -0.2162761940404157035 -0.45637269106764410953];

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