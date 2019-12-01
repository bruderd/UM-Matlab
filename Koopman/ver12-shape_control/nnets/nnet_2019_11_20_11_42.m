function [Y,Xf,Af] = nnet_2019_11_20_11_42(X,~,~)
%NNET_2019_11_20_11_42 neural network simulation function.
%
% Auto-generated by MATLAB, 20-Nov-2019 11:42:59.
% 
% [Y] = nnet_2019_11_20_11_42(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 7xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-1.56771344656424;-0.29596390680745;-0.168266734365819;-0.506463872079437;-0.469372496884669;-1;-1];
x1_step1.gain = [1.03915875162499;1.7292838457274;3.36016357505227;2.07784881685257;2.56216104668591;1;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.1289905823643864302;1.8817630377695959698;1.9155483826256005653;-2.7949208907418738868;1.8720020840841609111;-0.046450471300387409512;2.0025184613915580911;0.38903768537189065135;-0.2252636945693590742;-0.19850792117209789533;-2.084426715465462987;-0.10121893458857564307;2.1539009042753995971;-0.34511199641662304582;0.89617826284226487221;1.1452580243632894952;-2.1009737569538411428;1.6250359240932144012;-2.5013633316812802221;-1.3179624627826751482;-0.8031732629369948695;-0.11229862572359228245;1.0605313343724385167;-3.355806821539738749;2.5570035611174724011;2.2975350941268488647;2.3899559389399271936;-2.0990378788348262873;-3.0621728225618993768;1.8475881208069184769];
IW1_1 = [-0.039931899015356521476 -0.23021055084439245864 3.0042639124998151701 -0.13431610728604387317 -0.044846466097507253767 2.0549265694971854401 1.2409452725279845886;-0.26954247574855133562 0.52196559497855510301 0.24229447177775409461 0.33886090594992296365 -0.081180547149613396796 2.8501983011544518298 1.9959997699809934879;-1.2874696946117003549 -0.60310591334150154275 -0.0047747134720649009321 0.99180194972127333486 0.58352467849166356206 0.8252924028214593255 1.2676833389633310745;0.76159928744093696729 0.3699509575878723866 1.4002125310031647842 0.7349250742566632022 0.75448713797846300722 -0.74497364328280002432 0.066360789621519061798;-0.88843220764402019984 -0.032171983106192195145 1.1853676898949028207 -0.01148115740145356535 -0.12262264883923734704 2.3019589180895696501 1.263684510380199022;0.073578420184635437429 -0.72071054035441683983 -0.04280718304745628755 0.18110383553451589456 0.21149114548244501988 -3.0967839237204035641 -0.25693381106750229348;0.26745544737993048878 0.1961926350703638855 1.3299868751129164668 -0.94613051831889072307 0.83241093041586677881 5.547429497870989934 -4.9311274084714495913;-0.011867805397241004348 -0.48334130914941930701 -0.16679593254450328366 0.20395464460754708647 0.16358091911283859088 -3.134033248316975051 0.29271187895637429399;1.2150892063380231978 0.44911849805627135357 0.96680919357335437514 -0.78736526213011548414 0.40088319612362138811 2.3457812668870068329 -2.6477478650178620789;-0.1262330782954158015 0.91378677599316160407 0.36858131339700961338 -0.34176480892358529484 -0.29636645824109064673 2.4542504179016280119 -1.1571208701069153335;-0.20426988295190659506 0.33970947952165597705 -0.70200948997660383988 0.35369641336371165918 -0.12384262324323491722 -2.137372758553013341 5.0598657032145126777;-1.0243262884933150847 -1.0984799869228181191 -3.0631302730612377516 0.51150588748222691748 -1.1833132352414506006 -0.5691144668684946506 -0.23952631329668783167;-4.4083076840376200778 -0.6392075288842243852 -0.016650990671096699974 -0.022553324204756482196 0.32392184149826191986 1.2522232195774003483 -3.2739007019278272281;0.99579265460952848255 0.61008090468471931622 1.9752777760776114579 -0.21657578267484503542 0.23967470966457876669 1.5570269533063421097 0.58407105580506102349;-0.74426217493993218621 -0.080317445051893171515 0.48822383541017233544 0.82382948586659277002 -0.97005328921871059045 0.21845844276527867645 -2.3668256881900853195;0.50358628129943761031 0.83985099743735125433 -1.6957990894863077536 -0.049350275941495470866 0.15529607900032163514 0.93754093009205641263 -0.6259380808213708125;-2.2377504583344443034 2.2424474752262226751 -2.782488881448829332 0.72652493269073625815 -0.52443632946628682134 -4.4233336635241542467 7.3787431595259560524;0.041332454342630539634 -1.220706849163163854 -0.40616146216186010953 0.2629822161344118614 -0.66188777987833680161 -2.065031632604967804 -3.4309172446779583332;-0.073238122514383413586 0.44251561726215449255 -0.71566308225124131148 0.41090756328350513948 -0.14097277661891557954 -0.81952720269645495943 4.8324176627623156577;0.50397348606934389004 -0.53148569629685127325 -0.59373531465343598601 -0.9967281357142701026 0.62528281329003765343 0.70798250511857829803 2.4442470381393492929;0.097869957718077715825 -1.8754736957876338987 0.32327619187032397008 -0.32344655475174249881 0.95242314328651811461 -1.0486530017093662615 0.40327428806405907213;0.029724667092471169366 -1.7395689337303568145 -0.68129127575175041365 0.22221380989284980623 0.43935877744764451069 0.19198760765047470644 -0.14878781057877790706;0.43500803827399225909 -0.013237360928704387797 0.86851323925176648899 -0.18363064810699222806 0.40951719323561902586 1.2650340232461882994 -0.12330357047409001847;-1.147828278683788561 0.53360816163600066631 -1.3092013150520169784 0.40579575171536835798 -0.81504508349193971473 2.7581554486831612039 4.4573039681429627024;0.29065764746157524145 0.9678820520704239927 -0.15290855800560287925 0.54101054568290007918 0.35275121969564054991 2.0642698566313986142 1.3565932223365086973;-0.4751945581277773889 -1.3572147910799430282 -0.3142180351135411831 0.51802201195255825894 -0.89958555435125875377 -2.6230693810501377605 -3.0028922518351901516;0.23083248161372479923 1.3429688251941214183 0.047553784752305240657 1.1060458880890915445 -0.033751993297040906727 1.0590413614629641259 2.1090934357353647854;-1.0624524391725462813 -1.3312323795930047776 -2.2695342978235304443 0.61478324016154517473 -0.091345496352471897139 2.595187012429313711 1.8907113675893516547;-0.68356966223944126693 -0.69195406436219775159 1.9636143796026397368 0.31667749325727756471 -0.13492371175839845687 0.50782626843327183241 2.7820658779630220714;0.99615011726171154294 -1.0045225210852768427 -0.84240615602514179461 0.58061377791832036177 -0.89970133357960702103 -1.341096185845320532 -0.41308525428619502495];

% Layer 2
b2 = 0.64245832525783541644;
LW2_1 = [0.91271525366357275644 0.92864761212760138509 0.12839177851224481408 -0.61070605836520452936 -1.9107178068682266758 -2.2435373091041084237 0.44898299099723382177 2.7033975315709222187 0.70048041814277484107 1.3900920917183325898 -2.5347421358549975245 0.2964364756388828126 -0.18489189441845643724 0.67306244106940249505 -0.9997279654355610834 -1.6073927436017423354 -0.25922093118660133326 -0.82892504024685387254 3.4852890536215910267 -1.154911287924601293 0.56419895034022471414 -0.27583912067504429144 -3.5572815527773369304 -1.1781033805701386363 2.1529522609571847269 0.67255522104802056838 -0.70838357770192550422 0.48773852513493282101 -0.987134831443958638 0.04456152189348590309];

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