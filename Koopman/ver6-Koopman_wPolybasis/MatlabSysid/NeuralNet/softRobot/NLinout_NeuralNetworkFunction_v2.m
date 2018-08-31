function [y1,xf1] = NLinout_NeuralNetworkFunction_v2(x1,xi1)
%MYNEURALNETWORKFUNCTION neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 24-Aug-2018 11:53:05.
%
% [y1,xf1] = myNeuralNetworkFunction(x1,xi1) takes these arguments:
%   x1 = 3xTS matrix, input #1
%   xi1 = 3x1 matrix, initial 1 delay states for input #1.
% and returns:
%   y1 = 6xTS matrix, output #1
%   xf1 = 3x1 matrix, final 1 delay states for input #1.
% where TS is the number of timesteps.

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [0;0;0];
x1_step1.gain = [0.20053699350483;0.200544041106296;0.200280392549569];
x1_step1.ymin = -1;

% Layer 1
b1 = [4.7019649867488038;1.4088342045752278;8.6371245836944599;3.4449651157643211;-2.0108992210885468;-10.827313702142458;0.97249870811438721;1.2762057142621106;-5.022042986198282;-2.8286742381894801;-2.0509719112595697;-1.1213536898607095;0.1458411481234311;0.92291697014705265;0.38485739123281498;-0.25119899710888449;-0.52881646850474184;0.10197828671629132;4.8905802843188884;1.1822197472801905;1.4186502901931097;0.48574318013790946;0.99512849006982085;4.1801543190475199;-2.2190037104054108;-2.3563244128844039;-0.88452966786197584;7.6631826906085383;-4.5169901543093491;12.458132941976002];
IW1_1 = [-0.56786571115486828 -0.22629574747586625 -5.2753637004878096;-0.62319311336116567 -0.32017309272734973 0.07781856979932017;-3.5316892616043867 -8.0708174474065633 7.177467349284032;1.538700617934859 3.3819056259280087 1.1016702906600062;3.8861286877558014 -2.5723819588124202 1.4783622894513511;10.84760147359286 -3.7962461896264115 -5.3007261852940584;-1.1209469238516658 0.27686191663542364 0.2638097195956835;-3.1496310358324835 2.105271299051692 -1.2411526365362704;7.5118344866878388 6.421812872330289 -0.48410967553825052;2.7307322707741917 4.8599917418150671 4.7093718632246429;3.2547004221116689 0.3427253790783602 -3.6554985644707139;3.3717768633014042 -2.313341130818936 1.3432720303532608;0.80506039789804895 0.045853696815741606 0.87686237483753049;-0.27065477138718153 1.3688771071001951 -3.4458201318283068;-2.38359876644419 -4.1364388149219451 5.9695302264988372;2.8480178870418982 0.21462040502438776 0.45949393475894124;-0.39958478073876846 0.69492805915240996 0.049384572150602088;5.1227631937554907 2.3934056481395056 -1.8868713102939345;11.173902507318683 3.0143496903170357 -6.2695056872514181;-0.078991255875876534 -0.90766173736920408 1.5917657685186803;-0.030879562228806921 -0.95755820336056707 2.1915177430033315;1.058547334973815 0.23830694628187513 0.86423931569217061;-0.10046810773952161 -0.88090133558442052 0.97914153431756623;3.4909513281724793 0.90936123083031972 -0.70147143148084845;-0.28268285840692348 -1.7711915794404764 1.9108602953161695;-1.1844759471657327 -0.12935637275409756 0.2304439163758;0.15867664504136414 0.54885479023376682 -0.20715923244322237;1.6179132066763324 -8.093118617913003 1.3155837008336466;-0.40158547304465197 1.5623865792633949 4.3664684418387463;6.640810070675677 10.55517660849374 7.6930020704632307];

% Layer 2
b2 = [6.453591620792146;-3.1277589328502549;0.37629316093845544;-2.4342438155063695;1.3052526961295847;-1.0475365387549398];
LW2_1 = [-0.30852967354951283 12.800010934071096 -0.032110256620824403 0.022611183879493228 0.57084587127338371 -0.060020174817951528 -1.2738472118360451 1.6018271289527726 -0.08302848806882715 -0.038998314049940197 -0.14731705388080069 1.0131215354255778 1.8193452833901931 0.011402208356967383 -0.042100212063110515 -0.19322670059223848 -2.3232756655703626 -0.1699989530852605 -0.017817757626137323 -8.7382103951034011 3.1044976772718478 -1.2239318217534845 9.5439023838999741 0.94909220536526229 -0.21621110292163886 10.349219936365426 17.23296732118235 0.19862766709582899 -0.41792818495679279 0.039734177498948094;-0.025671481172788762 -5.5030413511410039 0.074151581060427227 0.203314701288348 1.6494593923664882 0.057434633344927433 0.65609116329950046 4.8043868905004405 0.011395049538714367 0.043518476382909975 0.052406554305436892 2.9245924877650022 -0.37278630692249598 -0.10020634774086913 0.039166923229749316 0.19855364476789519 0.99294409913614712 0.021218758873120549 -0.044544135349907955 6.9218427120099637 -2.3630753678463017 0.17556385477107461 -7.1816489686104088 -0.66656610482339806 0.091421738389973153 -4.3995675656736735 -8.9977567665791582 -0.060680188588462573 0.2384553561788387 -0.077452564204740837;0.05465873455265554 -0.46366632232746019 0.040342608357948749 -0.10482436369749554 0.47782096869163027 -0.0094456179701638797 0.62112348330140366 1.2511259487537232 -0.0098831519303474863 -0.025520812303936023 -0.014206975927527758 0.74932340564864797 0.15301135162185114 0.023777405333777084 0.00077545285880807049 0.10134129718426042 -0.28406212809762232 -0.00018487331496332989 0.0019619431612591144 -2.4496497617466333 1.1234152151827153 -0.24596469302264667 1.7587654904544165 0.29410105991366758 -0.21733949019132398 0.77089027873290239 0.57540096578488031 0.0790660817544694 -0.050979602641560977 0.042267601997120707;0.14345162982656551 -5.3529699253695036 0.0028896298394720087 0.0013747029202416826 -0.23479799196808238 0.028642672875898622 0.78652997932072044 -0.63173037676120081 0.043896622668867263 0.01197340109583602 0.045766094795530683 -0.4026586593534901 -0.66145916239394364 0.041282037602363826 0.008521616351752654 0.078773329808748863 1.4805051248879477 0.064936517438246938 0.0043065065766095898 2.7570722346794541 -1.0502340773190633 0.45475740644492246 -3.1057410381720345 -0.40793513131942488 0.088979725518969169 -4.2207737417030549 -6.78323126142344 -0.1020615501405998 0.13833721350241018 -0.016568242678502952;-0.059039187116380211 1.5871017326595425 -0.037801548137792178 -0.099392407117914652 -0.65428435060214152 -0.01292243267705925 0.028052214276732754 -1.9322941952837982 -0.0024615743608958338 -0.0023940662603184606 -0.0010870398260699398 -1.176306478410883 0.20713814321967522 0.0028481811294063923 -0.0012006034077365476 -0.10953847264312712 -0.26350685360936993 -0.0071906109227770675 0.010211788259166911 -1.5728010894088491 0.53885119925492442 0.040812442258196267 1.6791621612727503 0.23142609058610439 -0.010916961801449977 2.2845302702590105 1.8873091832561693 0.037078836305708754 -0.093373657483009007 0.046443714576734858;0.060148657921401547 -1.2539769744291904 0.00019177999590042283 0.003119396582456491 -0.12205119706871129 -0.01577596337762939 0.42794485287544626 -0.32808245569969879 0.0096103441860493841 -0.00083260092478646503 -0.0021574880278895331 -0.18964721684058403 -0.18268908328193417 0.018761983871162406 0.00057721212053831415 0.023673376024073253 0.038661367194451859 0.010632865159931967 -0.00032236800824938714 0.082172994410011163 -0.09687662707974963 0.018170113674803538 0.17164746090137831 -0.11649867714645405 -0.05454779480560102 -1.3547654521066952 -0.61903479183868604 0.047950833679119065 -0.0002780983918529178 0.0024881824086025284];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = [6.02925470600584;5.80239714539318;19.5174170890666;7.66459323490276;7.02419411493605;19.1826794012276];
y1_step1.xoffset = [-0.156274843142796;-0.148885540376267;-0.106342153265596;-0.133510159664189;-0.139378367975141;-0.0515426322556757];

% ===== SIMULATION ========

% Dimensions
TS = size(x1,2); % timesteps

% Input 1 Delay States
xd1 = mapminmax_apply(xi1,x1_step1);
xd1 = [xd1 zeros(3,1)];

% Allocate Outputs
y1 = zeros(6,TS);

% Time loop
for ts=1:TS
    
    % Rotating delay state position
    xdts = mod(ts+0,2)+1;
    
    % Input 1
    xd1(:,xdts) = mapminmax_apply(x1(:,ts),x1_step1);
    
    % Layer 1
    tapdelay1 = reshape(xd1(:,mod(xdts-1-1,2)+1),3,1);
    a1 = tansig_apply(b1 + IW1_1*tapdelay1);
    
    % Layer 2
    a2 = b2 + LW2_1*a1;
    
    % Output 1
    y1(:,ts) = mapminmax_reverse(a2,y1_step1);
end

% Final delay states
finalxts = TS+(1: 1);
xits = finalxts(finalxts<=1);
xts = finalxts(finalxts>1)-1;
xf1 = [xi1(:,xits) x1(:,xts)];
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
