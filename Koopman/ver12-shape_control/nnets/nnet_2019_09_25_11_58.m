function [Y,Xf,Af] = nnet_2019_09_25_11_58(X,~,~)
%NNET_2019_09_25_11_58 neural network simulation function.
%
% Auto-generated by MATLAB, 25-Sep-2019 11:58:44.
% 
% [Y] = nnet_2019_09_25_11_58(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 10xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.404691787807594;-0.380631614568998;-0.386224446278155;-0.203587195769292;-0.178867175199713;-0.212205404694113;-0.173319344073958;-0.180874614606994;-0.198749268335464;-1];
x1_step1.gain = [2.53283559356745;2.58858814749415;2.56215682071682;4.54288063098757;4.98693803687895;4.60844757321881;5.68156258045277;5.420361940415;4.73921520260683;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.9987962307783426752;2.4852365397990965334;2.0545861381544949786;3.1348648285048219542;2.3665039486420438308;4.1294794684604276824;1.2942653402917885241;0.57566228942825126147;-4.6476482316072651457;6.7124137189743100507];
IW1_1 = [0.79733844991977931294 0.18064960420694803345 -0.14398881953199560435 -0.21456189786997262314 0.24398517929404342497 0.51208164553012003317 0.41093755915926188926 0.24257700032508225374 0.042700121086593101349 115361668635393.53125;0.28493048925096048363 -0.23715453823410573286 0.59184468300019077436 -0.50599462034725017912 1.1544732007575322363 0.36187967886691879693 -1.2162995286813347295 -0.28102981350208938727 -0.049571389058234605063 -111656039571579.89063;-1.6646755502585026676 -1.6806813325834590866 -0.18290261947269817466 -0.35043108305582676421 -1.2077714789854352428 -1.4673455740594549823 0.16109406049760399693 -0.12009030195052264323 -0.50585583496844077622 -84591786582335.796875;-0.1471027579393652529 0.19759643267317825166 -0.06055779311200397419 0.00075124716600563949068 0.20457711999488481514 0.14498327382710568534 0.15155102460579272794 0.068846883176702317053 -0.15436379923203447095 -107792525097543.70313;-0.22073595881809773833 -0.64864495799030763568 -0.013647315499675623934 -0.022655810727164354956 0.35228471494876029135 0.22179935555056462526 0.11369248959380087727 0.13772332875609155645 -0.29277953844206111578 -122058547464097.3125;-0.2660482435980445759 -0.47595068829029940627 0.055012259030566022788 0.10768533884671550871 0.11445359940566479651 0.045699122603512012941 0.031833414955283738379 0.095034635416236037941 -0.31113322565252971019 -134931444103078.25;-0.83406426781290754491 0.08961086212873548007 0.068703628196658941829 0.019042609862221952094 0.09875828527303102955 0.087032924097462396884 -0.041879448129928213695 0.055065330387590351313 -0.11643905938441302672 -47151775430887.289063;-0.014802107709531440657 0.10229524224429759061 0.33823911752663510244 0.0011634305055793583462 0.029878072743731914906 0.013820038351593037437 0.042081577006865628077 0.0092959992329991331106 -0.026520046270495010315 -20053986279568.945313;0.19078280072541717316 0.19251662971934890844 -0.025133051313732258925 0.023064114540107434614 -0.32783904612022463088 -0.26470348811651922549 -0.14445828074223726611 -0.075274511904031768394 0.15270980000459180581 145946356324355.8125;-2.1625175999546057959 -1.38511060203423364 -0.47561867735962676251 0.35886097726809890185 0.42191650475996006486 0.52613880084955944039 -0.60889007673712325364 -0.20192680830481840593 0.69373772920251908225 -166545399191602];

% Layer 2
b2 = [0.52212820046605945556;0.35102300873108505108;0.15836561804566656741];
LW2_1 = [0.55952244645566606707 0.15031556465604475892 -0.14203313123858055111 0.99596758996431800082 -0.031529922979280304574 1.1783503841621034702 -1.3928363838452422208 0.17382675185344925728 1.7121520757062389517 0.17168094118316742924;0.52848690689140442878 0.13988207511776334924 -0.10895290968461147918 2.3672029115086830053 -0.93392473958168253922 0.37192206747715289561 0.19521825631396874656 0.28020141203837956034 2.4829924023507294883 0.20312193470716427557;0.10502487876407466671 0.017695831934936248619 -0.016314205588377359013 -0.95780365986366289288 0.18353524415075728182 0.16892398707599109997 0.061160379260873476426 2.8368455562335554987 -0.093434355296356727538 0.038459875594539540133];

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