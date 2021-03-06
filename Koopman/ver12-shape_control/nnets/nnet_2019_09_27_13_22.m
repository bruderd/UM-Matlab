function [Y,Xf,Af] = nnet_2019_09_27_13_22(X,~,~)
%NNET_2019_09_27_13_22 neural network simulation function.
%
% Auto-generated by MATLAB, 27-Sep-2019 13:22:07.
% 
% [Y] = nnet_2019_09_27_13_22(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 1xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = -1.5391474617351;
x1_step1.gain = 0.650402606650308;
x1_step1.ymin = -1;

% Layer 1
b1 = [12.453360117575924093;-9.1331768862749935778;-6.2891985590758014979;-1.5370375156555724505;-0.24971468996154139175;1.5784358216952636411;-3.7566599532577256149;3.6521986394598999759;-9.2439551513472952848;-12.854209362338380629];
IW1_1 = [-14.035636602936934381;12.028443869281428391;10.273142799958943527;4.0735290414976681461;5.8161454568084662853;9.2173273837150837551;-11.047928593168283129;6.4781192265272160569;-11.989280605225461684;-14.220541042207214488];

% Layer 2
b2 = -0.0063918655462897287778;
LW2_1 = [-0.072816351392919878571 0.05637243682308407805 0.047972793405459263072 0.2301354198554612096 0.15156218124073861375 0.080342062883215742408 -0.061684180999661109013 0.15629691906958667613 -0.067382552468055828543 -0.06675473687492444741];

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
