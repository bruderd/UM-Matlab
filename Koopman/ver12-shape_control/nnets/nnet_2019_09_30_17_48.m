function [Y,Xf,Af] = nnet_2019_09_30_17_48(X,~,~)
%NNET_2019_09_30_17_48 neural network simulation function.
%
% Auto-generated by MATLAB, 30-Sep-2019 17:48:10.
% 
% [Y] = nnet_2019_09_30_17_48(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 3xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.878104971371236;-0.718512783764323;-5.91431892442893e-15];
x1_step1.gain = [1.18608476967992;1.24865548709588;38763343309333.1];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.0685440600639850572;1.9672133162391394201;2.4387679014543999223;2.1198828242256078092;0.73987140893958935894;-0.61653173493829105389;0.37525858153229535441;1.5401133224626435414;-2.8519566312818986908;2.7654998066401579315];
IW1_1 = [1.6593285671791391245 -1.5919417160460540561 2.091291591230877156;-1.6825034492559043375 4.2302065576323455431 -4.2159654503759620425;-9.2628416061145131266 0.56511122242689992845 2.323945553882441839;-0.30672444781868946073 0.097852342451611756813 -3.0162690010006527608;0.42446199187033900602 2.1980580569178060912 -2.547344337325161856;-0.777146045668330121 -3.9191186277576699482 4.5742256574491682031;1.9784954667330754141 0.35052937537836648119 -0.50419500027457153646;2.1938136731047217509 1.3706828545570952649 -1.8593425257094331471;-1.986830276534108819 -1.4457976864493209401 -3.6192866275243416041;2.3917181555630762091 3.2067933659415546188 -3.3300031899574840111];

% Layer 2
b2 = [1.6814716304642693956;0.013799693451433718561;-0.38743058336219010629];
LW2_1 = [2.2235317964285505887 -0.95544661693149768311 0.31675075634756383325 1.0711025215584100145 4.2614370301838739863 1.4683674884502713187 0.83058590264872855791 -4.7036124452230385629 0.265233744383290726 2.3756981768285814205;0.061017218062510786281 -0.35547242569183268301 -0.027160107758901919778 -0.18549504790332829729 0.87135553589155123611 0.25588522129504986946 -0.68919347394964391995 1.5631149419565202319 -0.029139886169291272233 -1.0021308584696668742;-0.28614125200462681153 -0.047938528259469294357 -0.14201026862825744113 -0.1124545253508750231 -0.68110893965625873836 -0.42068656380861035871 -0.72030161758988398013 1.5456348695392190251 -0.078444903864325496445 -0.57118909600720801922];

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
