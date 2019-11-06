function [Y,Xf,Af] = nnet_2019_11_05_19_05(X,~,~)
%NNET_2019_11_05_19_05 neural network simulation function.
%
% Auto-generated by MATLAB, 05-Nov-2019 19:05:08.
% 
% [Y] = nnet_2019_11_05_19_05(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 6xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.0414680022291759;-0.0262478593634704;-0.0194075445527305;-0.0150152300885408;-1;-1];
x1_step1.gain = [21.4031105149736;33.788291631072;50.4516621593694;70.8605720879517;1;1];
x1_step1.ymin = -1;

% Layer 1
b1 = [16.229700756550975882;54.37152878163755787;2.9209573525290069007;-9.5361872593902763384;-11.704828224616882437;-0.35502004525281694036;-107.35035438857022427;-0.043336144268783856337;0.23956418539108645382;-4.6690233806485466062;-141.82543841295932907;4.1364337011496923324;321.72200146361149109;-6.9311740288709025037;35.470674214284208858;-45.058692441054795097;61.605607602044656801;-427.47420842844991284;-3.9299495459911311102;-26.087080454419584896;6.0886653996202504402;562.92551450158464377;84.293358391523838691;306.94333783803187998;6.2542944588465188005;13.009291561990181663;-400.3706610828633643;-499.9833289400738181;-6.5366684963888932458;13.022958369420353009];
IW1_1 = [-3.0888325369484435079 -7.7691162921074123915 -18.055958038985352232 -6.9747730395931020198 27.906470400852352753 5.7014802673118447629;-22.368284398334768781 -6.7703606868980861222 -34.327242010274403583 -31.485669212560591035 90.754287390674278413 -2.5724265842143405614;-1.6664212960665616414 0.61796020850219679499 -0.78846654624847856141 -2.0690716588494861838 2.9094750695746722791 -0.67182992871775759625;13.463806769135244323 -4.7149686573700959613 0.80099432792684510662 -1.2429403679096771018 -2.4363102893778516567 0.79923424873841408811;5.3824997943277974599 8.3644365131811841962 -7.0812175055402315493 -3.1455213827287957784 -21.756390777405112402 -1.9345445222378965244;-0.23892352772472910671 0.19798491259876174753 -0.27442350678926302754 -0.90985333131452883482 0.24460624288803034521 -0.010273555004087829221;4.18152054630033998 55.242992521902216652 -82.740839588588130482 -43.615329551591351276 -186.90148800122528883 -9.600634126257800105;-0.20150385561947029012 0.21340548072211648623 -0.31775825902309834392 -1.0084756731614694658 -0.10443952704294622025 -0.1239439908866679374;-0.37158854198785190581 0.27356902981694591892 -0.22609400087265338342 -0.86467454052551262489 -0.42688472623375844917 -0.22209406530233247179;-1.4783668269570944265 -3.3740491074406291361 -1.0734662552143561243 0.54838779670955128953 -5.8421510694201401748 0.23240387465598369432;70.53750389169353241 27.275006446066907273 -97.046209479351603022 -49.368173023194145799 -249.49782916456666726 -8.7403825851888754528;-0.85758795834802659108 3.6895594674272933489 1.4570790036382232913 -0.96543764605482595975 6.0427042609847108068 -0.21364836289495081623;0.091196841676578654923 -39.917866159591788744 -59.061578534296423015 -33.481357758092443078 523.99210981153009925 3.706500513153471843;-28.590633552710428944 -9.1886983842548684009 -0.41082495291244242042 2.0926671095386866206 -2.9392002330415314404 -4.6406254404783870982;-3.388311281903093608 3.0613832325972483162 24.514662395653800786 4.2526241519938503544 57.233559815367911483 2.6283795670113683229;11.95489210168076788 2.9515305124080097876 -37.614280402107965529 -3.0399918446279596829 -74.677057484142224553 -3.5081940174465322535;27.984523405337213831 -2.1739778977894466117 -14.544979338061267526 -30.172769330822813316 92.259547945876775543 -0.51893986557215743272;-33.283347045870648628 27.27665237352557881 150.05069337765064574 71.942420732232434943 -685.55489134496008319 -1.3275983169122393956;0.68899295147199179912 10.038950954276231187 1.4200773931084844648 2.7374831042499305767 -9.8827633319416001711 -0.033827977486208860025;1.8797085165194622469 36.025587143726717443 -20.110439764922986683 26.830656856479148331 -47.488682593794848685 -2.5134123307861506724;2.172217638856030586 4.7933686754283542086 0.90192105243213938071 -0.67293061922733654612 7.0475809430547942114 -0.33420879640763168084;-72.820792396420472414 -67.97815107536865753 30.485571393408463337 90.190455593989099725 942.64874351793844198 12.365555647989635801;-55.513937406344581404 -22.897405732951686019 51.622868776038181693 23.906445681891391075 151.6470849139267898 5.5469808069898673963;-13.686923178605757911 -64.442117563258520363 -81.034229909933685576 -28.550189420441604682 507.94907754712522774 3.0752165453528679251;2.1193038145464009325 -7.3663784743156712764 0.64126195818300923968 -2.5985746944042764639 11.080516002049623836 0.025995996406612556506;3.2046493656831382957 5.0149822677799189563 -14.881487977671472578 -8.2797316200633481742 11.944666713068734509 -22.337876218632530367;-3.274024368036236865 60.531779995499469749 73.395347948476114652 25.528340176386993932 -655.85124619804776103 -4.9743686121934063138;50.344360662901053161 63.614599637347374994 17.274919663981101792 -36.248952321196505011 -831.81609801223214617 -8.7334238751231954012;4.2685475280265148257 6.763517381046206367 0.1126866398402178171 4.8329789659085538034 -12.594660061107894933 -0.03077399405986189565;-8.5913018842214885495 -4.8179237333421101042 4.8923527625723233214 4.3985567328827004374 23.723511071674710848 1.5891413326850316601];

% Layer 2
b2 = 2.6957477389382598076;
LW2_1 = [-0.69376338457824959782 -27.446163138908801216 -4.5357026044750075044 1.0854134502540564444 10.369261808462260532 -8.4291900729458646424 -38.113888992423056834 15.51145711949120809 -8.305264883667449638 -50.127220071620158137 15.356754222391458597 -17.772705100051318539 177.16000519461803719 -0.020564061122601957898 -9.9251988661220718058 -12.855939124428127585 27.901701335754705013 -46.247222467412569813 -1.2096718998945961143 -1.1006693313604014595 -24.952400683604118115 -78.941421246040761162 -29.244014673844187513 -77.402918235544007075 -26.000606296591435296 0.020858433568522839058 236.49598860611280315 -170.93235626736080235 -20.518567493076165675 11.284993925839142292];

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
