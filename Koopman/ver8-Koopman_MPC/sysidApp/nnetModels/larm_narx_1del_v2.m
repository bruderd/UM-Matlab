function [y1,xf1,xf2] = larm_narx_1del_v2(x1,x2,xi1,xi2)
%MYNEURALNETWORKFUNCTION neural network simulation function.
% GENERATE IN SAME WAY AS V1... probably not needed
% Generated by Neural Network Toolbox function genFunction, 22-Jan-2019 15:50:58.
%
% [y1,xf1,xf2] = myNeuralNetworkFunction(x1,x2,xi1,xi2) takes these arguments:
%   x1 = 3xTS matrix, input #1
%   x2 = 2xTS matrix, input #2
%   xi1 = 3x1 matrix, initial 1 delay states for input #1.
%   xi2 = 2x1 matrix, initial 1 delay states for input #2.
% and returns:
%   y1 = 2xTS matrix, output #1
%   xf1 = 3x1 matrix, final 1 delay states for input #1.
%   xf2 = 2x1 matrix, final 1 delay states for input #2.
% where TS is the number of timesteps.

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [0;0;0];
x1_step1.gain = [2.22222222222222;2.22222222222222;2.22222222222222];
x1_step1.ymin = -1;

% Input 2
x2_step1.xoffset = [-0.897930439478848;-0.9];
x2_step1.gain = [1.11239008811694;1.16699751868495];
x2_step1.ymin = -1;

% Layer 1
b1 = [0.26889002327086536;2.0194568533060546;-0.068307040406867123;1.2482425818189693;1.2080296810632731;4.0644123277222191;1.470373628807891;-0.26787172090462275;-0.071293590864110218;9.9868186051529193];
IW1_1 = [-0.0051817071319579175 -0.046354814831364034 -0.019741099886148122;0.020039056997330579 0.19217177746898806 0.040911808212126491;-0.00076132241760973264 -0.04020271475710506 -0.027432486334340377;-0.12799712771370722 0.25579588624838751 0.33604962177117625;-0.098752673702310875 -0.086433925397371783 0.29515692561250861;-0.40464282345390629 -0.8539087512113076 2.7078313333540751;0.75692288741875935 -0.15623428125017449 0.85196055485458222;-0.025010967955894189 -0.024282056032403881 0.05696442306625029;-0.039191970919062173 -0.023672930278425292 0.1010979814021961;3.194075566959818 2.1984473410042984 0.11054768430318644];
IW1_2 = [0.10872620686372532 -0.13732907267519043;-0.29080398118316464 0.38793059884856146;-0.035586830175471582 0.13832739483067213;-0.83144334337995784 0.78057643217266492;0.46732382658919447 0.53043029309893852;-2.0896568892985936 -3.0105228768271441;-2.2472238868582859 1.8834967566602108;0.13218891845273972 0.12670212162061825;-0.2022814674688416 -0.056215752109212216;-0.70755494492272775 4.8604048369475148];

% Layer 2
b2 = [-0.37142218882414829;0.08932710798623747];
LW2_1 = [1.4300530665321358 0.43391221999531387 -2.1301363961024142 -0.05909098611857496 0.16580596377559853 -0.0053866130631109914 -0.0035784578706588061 2.8901388797719334 -1.9374798764431385 -0.028685614571501577;-2.6522894746688292 1.168117791000542 2.6436675969627692 -0.014272743471308694 0.098639751881498844 -0.0010059353710143411 0.001604257900585172 1.7221346058796683 -0.79012464808830363 -0.015066165477250129];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = [1.11239008811694;1.16699751868495];
y1_step1.xoffset = [-0.897930439478848;-0.9];

% ===== SIMULATION ========

% Dimensions
TS = size(x1,2); % timesteps

% Input 1 Delay States
xd1 = mapminmax_apply(xi1,x1_step1);
xd1 = [xd1 zeros(3,1)];

% Input 2 Delay States
xd2 = mapminmax_apply(xi2,x2_step1);
xd2 = [xd2 zeros(2,1)];

% Allocate Outputs
y1 = zeros(2,TS);

% Time loop
for ts=1:TS
    
    % Rotating delay state position
    xdts = mod(ts+0,2)+1;
    
    % Input 1
    xd1(:,xdts) = mapminmax_apply(x1(:,ts),x1_step1);
    
    % Input 2
    xd2(:,xdts) = mapminmax_apply(x2(:,ts),x2_step1);
    
    % Layer 1
    tapdelay1 = reshape(xd1(:,mod(xdts-1-1,2)+1),3,1);
    tapdelay2 = reshape(xd2(:,mod(xdts-1-1,2)+1),2,1);
    a1 = tansig_apply(b1 + IW1_1*tapdelay1 + IW1_2*tapdelay2);
    
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
xf2 = [xi2(:,xits) x2(:,xts)];
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
