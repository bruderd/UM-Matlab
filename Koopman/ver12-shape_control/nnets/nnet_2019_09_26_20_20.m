function [Y,Xf,Af] = nnet_2019_09_26_20_20(X,~,~)
%NNET_2019_09_26_20_20 neural network simulation function.
%
% Auto-generated by MATLAB, 26-Sep-2019 20:20:36.
% 
% [Y] = nnet_2019_09_26_20_20(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1.61601525396146;-1.40859477467407;-1.50579991371487];
x1_step1.gain = [0.630965482559027;0.726399427713015;0.690404926467999];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.3305208366817038446;-2.3682320945464336859;1.045228666791541583;-0.43364777218436156847;0.39311480273898308235;0.026494847091875786682;-1.0535389927070903937;1.6884969107992284787;-3.042844601731540255;-2.5659098103732027596];
IW1_1 = [0.10555680852248433554 -2.3501043402437353791 -1.3698255414554871834;1.2137340671954259808 -1.7807280220633443513 -1.6865400545192601456;-0.41649995379297383158 -1.1000576860028061787 0.43852338755345893206;-0.23455551202002461464 -0.65851094934692666616 0.28121503563639671519;-0.46630784976072242598 -0.021758249172718003434 -0.55185026769181910744;0.27096701930000205438 -0.33983390858079193242 -0.33536955862374617787;-1.0038648579828004515 -0.10584500057045635735 -1.1199964864236513495;2.2223457594144782057 1.1837674777220315026 0.72792456454267839838;-2.2247621123348233496 0.21737539450103901539 1.3287351508107809739;-1.4074057441638430088 1.3750272372476195493 1.4761626101634746089];

% Layer 2
b2 = [0.30954133189056598496;0.080680546682275186554;-0.12490462915888936313];
LW2_1 = [0.059997298223119777116 -0.0064680194658854528109 0.014935656866206513047 -0.036928168566246637183 -1.9035191789604544166 0.25052760429405973586 -0.48851980716834902951 -0.018965090168088654526 0.016634977467358310482 -0.017107107107404696117;0.049981365261387293575 -0.117864951843726945 0.15405899943763104454 0.44200295744557566158 0.28660770159189463691 -2.6107902690458084649 0.059143189742290110111 -0.011633210869901126558 -0.052545444390656746714 0.1505588638761468645;0.01033906936802676986 -0.010669268789159528504 -0.60710126771639061349 -1.594602099455989519 -0.18174987771930498082 -0.011003607200782736242 -0.027166751893902163495 0.014801156689046329751 0.062889235086791228646 -0.02305662292981600775];

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