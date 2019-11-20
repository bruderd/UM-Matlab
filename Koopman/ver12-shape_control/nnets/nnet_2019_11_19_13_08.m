function [Y,Xf,Af] = nnet_2019_11_19_13_08(X,~,~)
%NNET_2019_11_19_13_08 neural network simulation function.
%
% Auto-generated by MATLAB, 19-Nov-2019 13:08:38.
% 
% [Y] = nnet_2019_11_19_13_08(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 3xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.204085338122521;-0.0755507633573603;-1];
x1_step1.gain = [3.81267506288467;6.18782479625196;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.3296079145594608129;-4.4946558847449677998;-3.2251133479256677283;7.5086645833018499374;-3.1199150527809833555;1.2161682352284131081;-1.4836010260682090944;4.3971014300569857269;-2.109232793484722368;-0.21190308491688730763;-0.46807600898020573776;-1.180168714247733508;0.35172213074966424129;0.042693004733138362194;1.090620011437330561;-0.090655766589954125956;0.82613494489642746998;0.7915538321492285867;-0.69474994248304988709;2.8915200877635327537;-2.1477698072976840749;-1.0363248694236195657;-6.7106265275167471529;-4.0732665714642424248;-4.3215854341065895738;-4.0943930420126086744;-4.0895059640772633003;-5.4446074678378755252;-5.2346529173562617032;-4.0544320013588892593];
IW1_1 = [-4.0086143088363668241 -2.1716923058767596544 0.89313769118086827614;5.4081483277236417351 0.25246716568927407565 1.0914007654868866926;0.728359462556553372 0.10718520295417828525 2.6577353145166084936;0.86982758923303216125 -0.30194590138921062472 -8.3475659466549227972;1.7519481867392121544 1.4425424783095810444 2.8534921379271276365;-2.2420971339833188907 -4.2825617895358085718 -2.3581273938235152166;1.3506650485503377546 0.68541443271030744189 -1.211606416911218842;-1.0409062901048307115 3.4068568172044311204 -0.52313257248213473893;-0.036316860908548839748 -4.2772642873757478199 -3.7171149954165416673;9.2958616922036796382 -0.61063340297707613402 6.04732726117908026;3.8258093605962604222 2.0664650451404487086 6.2758355777020380373;4.5802886187980824673 -5.4810090749503652674 -1.5003071063438793153;-1.2791337157999100871 -4.9442124392011761458 -1.8011585480827076999;1.2911064467247133081 1.9979702949220887742 18.522226832034807131;-6.565396510864336399 0.23955760456333624608 2.2938324849210922629;-0.7086401045104361307 -4.0117555673011189299 -0.72087487560560647992;1.1751184985689397955 3.1396367024093021492 3.0678263709475053211;0.12504746515103390414 -5.813332437871204128 5.6214305457957980394;-0.42363327272584738736 1.1864756999464620701 -5.5109268895530139076;1.6441947341491649492 2.176387994377633639 -4.6188047998240211456;-2.6200449798806744184 -2.7065720270457966912 -7.2491504827937163213;-0.55839576682260516627 -3.3923018589453479699 -2.6017905022984137986;-1.0106430685113174661 0.41570240723282492201 7.2842380208889379034;-1.4266656351033228933 0.62862533557318256427 3.9290713158475658595;-2.0585251346645341286 -4.1736748166658470538 1.5257362774933307392;-3.1317098533975000763 -0.088959230185886417619 1.8655479301518198554;0.15323125789672431551 -2.5502982084348801983 3.0687011158784587295;0.16213962167318465846 2.8911965881635670605 -5.7592032714721463904;-0.95657360700320026581 1.2927200332066610677 -2.3027378485934706376;-5.2699423580456645766 0.045543001399366067616 -0.27235188324039821817];

% Layer 2
b2 = 1.5146648698520568388;
LW2_1 = [-1.1796368882754380003 0.1082011103024838744 4.4765062541642093308 3.3818095539916654957 -1.7133188283752354142 -0.26374696136960484871 -2.1145702728798827508 0.90229700781491195549 -0.48014439818859516729 0.22376060653245333221 -0.36101783792574249077 0.087384799250390388092 0.67607573252654884488 0.1796370141439667667 -0.12961721403299128474 -0.81419169087852905076 2.1422279902594434375 -0.10899789617353432269 0.18090458712958601417 0.16707139642525831169 0.28886674401852274752 1.9371142007749622138 4.9363441264551086718 -2.8213723816776372644 0.098865357070070603029 2.8172327317879486408 -0.5321991975320237156 0.82281123003256728676 -0.59274787807466988721 -1.4709300728614522757];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 1;
y1_step1.xoffset = -1;

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
