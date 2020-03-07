function [Y,Xf,Af] = nnet_2020_03_06_18_56(X,~,~)
%NNET_2020_03_06_18_56 neural network simulation function.
%
% Auto-generated by MATLAB, 06-Mar-2020 18:56:17.
% 
% [Y] = nnet_2020_03_06_18_56(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 5xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.874231897569386;-0.739350650421629;-0.429952038078063;-1;-1];
x1_step1.gain = [1.50523225798103;1.26648817538866;2.59657386551262;1;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.7068999289761044302;1.2824706157747967072;-4.8314317104515884793;-0.24741664035479579709;-6.8649599591173435797;1.723090824500531637;-0.2822181618918063184;-2.031397827860204508;-15.095716095675500767;0.29566461782599146035;1.2828175030636226328;-2.5684412361856194806;-2.638458510897274234;-0.43931166650818120889;2.8588491804739430968;-0.30768934258337454768;0.91706484488358308571;-5.0836984775975304629;6.6953435781156125373;-2.2329360419584136643;0.1669578946790739038;6.0111959145704085827;7.2342669044220668795;6.2970992248321282503;3.743635602118637884;-5.0371304937465559703;-1.2194310460905484206;3.909804393064432837;4.4941240637525856627;9.6239590352397961226];
IW1_1 = [3.1721085249440679199 2.0858325997697022558 0.48019853154436537013 3.5854572247676741448 -1.5891221404161903763;-1.0636388312970881387 -0.65291277907062705843 -0.25741064054013174101 -0.18583636639669731672 0.91903374479551491749;-1.0021226637002653703 2.2238121373337271791 0.68147729836682324756 -4.3129783808991577487 6.7203473311194148465;0.59221965325676306779 0.33120969355477847262 0.84878201227813654395 -14.706845169980127253 -1.3506829229811996296;4.9316817972903335132 0.99279679148480670658 0.73935554795243063619 -13.357007406030657393 3.1647922298994224377;-0.062554229094781282905 1.937993949087307044 1.5432594666668411065 1.2432164278127333379 2.4224614122376868508;-0.52340782205190372611 0.14694357638464936033 -0.59651326745544863694 -0.61084189837513158672 -0.96886001870308591233;-0.46093534512414030546 -0.57947342991372208232 -0.52474328807913173112 -1.5166142220788416672 -4.8442135799913215166;8.988075505652096453 6.178282062910992245 -4.9038302978428633949 10.178137955420593741 19.402795494309671653;1.4540550404672754059 -2.0183726531742474108 -0.40733804205865381531 0.48077837285117269728 -2.7133787665704942604;1.2006011132496612337 -0.0055050634258096212137 0.73557137512223536024 -13.474891576530598769 -3.9959005270232816898;-9.3178239170939818337 5.8938164804014343545 1.0107269878284963482 0.94580998941292826476 10.76327771641198261;1.0780530896658238049 0.55522062621551104833 0.91113698475148197087 -15.492857112212110593 1.042599191283435589;0.26212140033181080545 -0.0019140429145257029098 -0.12025898298963849653 -0.59480488312558366104 1.0227946799339033834;5.2856148739415980131 3.7069935777314384495 1.8738278922440585461 1.1008265023312264752 -9.7970974174551503921;-1.3662829196028243306 1.4641932603750122954 0.35068305408572181259 -0.28206033715654132976 2.0730803286286922926;7.1334250662039124791 -7.9427796295396069226 3.6973125776324966552 18.096255801641316197 -24.302141072145346357;-0.80508771669744572108 -1.9243425310166735187 -1.0499221604129502783 29.907827361583429138 11.167325420896231591;1.3956412760583774002 -2.739119734123167671 -0.95533224803999616803 6.2528089338043564283 -9.3719110677824257749;1.9505318236104312035 -0.79598055324734595661 -0.95767755151587030671 23.01609409609181256 -18.575347935746965078;-0.31769390374054801152 -0.24692646231353898578 -0.007928658059038203032 -1.784869477821952044 1.5640239001090632254;-2.5520330097048948836 -0.85265982552811048212 -0.84382557074566177935 12.108360336019012848 -2.8459947086378845071;-1.1405526791753171079 -3.3849453259563428986 -0.2688537796247488787 -11.456671015645069289 -1.6920897639005385571;-2.6294535125978448775 3.1524312892596109137 0.84983436274910328923 -6.5790639895875342447 8.7023909012850282352;5.5000878445583021303 3.2962157416073933902 2.2436208623849909216 5.964554285287865909 0.92509659008479305786;2.7871422293635590428 0.63312801343105673535 0.69809583050170098684 -9.9829370815768889713 2.3125269154816932016;0.60010455799717044023 0.47270997280354942749 0.016255720263038071771 2.9758464712576446765 -0.94294718239353958467;4.7639495284725743574 2.2455998378293884343 2.434168856004774284 5.9816592950298286269 0.69085416399996701919;2.1760703442092017568 1.7812799741753615468 0.41850347149246219791 -3.0815974663328660199 -2.2660773483635123782;-4.7581934454860492778 2.4583856734987326043 -0.60112119064249691647 -11.33849218708925477 2.9107192828138104268];

% Layer 2
b2 = [2.2043702281617068017;-2.0699244625204622494;-1.4335255635870665092];
LW2_1 = [-0.9920904909029018004 -2.817668284894208508 -0.94164447143506280202 -4.6221844850923252679 0.56429507247203014586 -0.45845629623965389321 -0.61636436407735939458 -0.12269760137332813565 0.011836266270773372275 0.84553650100674260237 2.6390371451295830951 0.0048897800724336613204 3.3612144733394666574 0.10852737254955209523 0.2006898643750322786 1.0683772689026829195 -0.1134628032674271092 0.75651235027206897854 -0.70833667333424543688 0.19854441519733712074 2.9355700722358410992 -1.8128305933826858176 0.2169478648181516478 -0.23532507061204999133 -0.50664855471025871392 -3.0292498953947970541 3.2033769643360763268 0.54180685727858146628 -0.3971824520293691374 0.23774527220619348644;0.95773692692480882904 2.6665960246938653633 2.8710680092707123734 4.0540936466824444295 -1.249500210158131619 -1.09257863434649094 -1.1133953294402896272 -0.36864518929444956008 0.073770564959677570327 -2.3184887318226743425 -2.5558359677513977637 0.12276996442128715326 -2.9236319926965976812 -3.1492194903816543849 -0.25905736771159643128 -3.2677350557953963595 0.087066537354387213421 -0.8188448009286324325 1.9933502307540120757 -0.37279932109634478232 -1.4225325207569707864 2.8208197000651242803 -0.24201606212192905421 0.14416661946685152551 2.393719494528558478 4.5665691309917031404 -2.5475118253904707188 -2.6079646230720010536 0.28776975359627149365 0.050294544620127293322;0.56833296407421629581 2.4657119256233919202 -2.2652694652294478317 5.1732957399002748033 0.094463817850830369971 0.58770221202532213489 0.12600552566091521856 0.68416948784595110755 -0.15496550634822844894 3.3616850617990721162 -2.3454585420798497175 -0.2887625809150113354 -4.35089482850928988 3.0202509786146780613 -0.079489731479042452511 4.6343938720211914628 0.078448096575371856209 -0.47457084142312683239 -1.4323358968685961923 0.063090489651105932656 -3.6006198351673592839 1.0425108916263172709 0.27708691878721009028 -0.46390610020035094596 -2.0692552624339253242 2.1994710184312125101 -1.866810649028339153 2.2068763669105635472 0.45340524260650805211 0.10584430026144830972];

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
