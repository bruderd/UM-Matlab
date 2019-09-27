function [Y,Xf,Af] = nnet_2019_09_26_19_20(X,~,~)
%NNET_2019_09_26_19_20 neural network simulation function.
%
% Auto-generated by MATLAB, 26-Sep-2019 19:20:17.
% 
% [Y] = nnet_2019_09_26_19_20(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1.57328270982114;-1.43269065472291;-1.38738492688586];
x1_step1.gain = [0.62233024158695;0.690336905052307;0.7434131339883];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.6598714160355991787;-2.5756488996383475332;2.9411543613177046907;0.27388400418112285317;0.066995098463803040567;-0.016246235591523255221;-0.062617840480421016469;-1.3468829093776635286;-2.4370794487334745959;3.5117888757534210242];
IW1_1 = [-0.55528079638468696988 -0.13436471298373703354 1.809459016327986669;1.4732777983323437354 1.3924289179370374914 -1.2269511561681356593;-2.5730711027042700856 0.26025212553597126819 1.6439732595131189008;-0.10797877896563423672 0.55504240104447977178 -0.022060640010737493644;0.040065407387329193822 0.083288807394165553788 -0.25786396281772416827;-0.22721148477911212726 -0.27541731597425928024 -0.045009678187717228059;0.15492435468344978311 0.34937320287496537174 -0.14504704969691478711;-0.25364961160417009989 1.5288836065662032571 -0.015670027344258367907;1.284955578905748963 0.52880058497124771133 -1.0733780876935301585;0.29702064491261059853 -2.1589808483459567157 -1.5574189816593011493];

% Layer 2
b2 = [0.77130812631404221946;0.1972447654012381979;-0.56731054349866960607];
LW2_1 = [-0.083311338511421792474 -0.085665856717361613115 0.039290628256104839777 -3.2663620182000814296 0.83132390436414194124 -1.8843166955094952097 2.1062712578260245344 -0.47693846131479600281 0.20415182396128522813 -0.083789948377279072078;-0.78787503224563537607 0.37058712805899907128 -0.22744152426310110249 0.22514884178485444766 -0.33081718499836920699 -1.3887696338033881371 2.1387218946104278317 -0.033897147402491410728 -1.2878295474602310389 -0.0097000443725538398887;0.92121749025237120723 -0.039159479529654594998 -0.012552468793669949859 -0.76804443651073772692 -4.4993228337681285822 -2.3119400328935575395 0.40606802756589893688 -0.13152662653288657579 0.0045800184013265604621 -0.022506087540313968187];

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
