% setParams.m
%   Defines FREE parameters and measured values
%   Should run this before plots_FBstress and elastomerSysID

clear
params = struct;

%% CHOOSE THE RELAXED FIBER ANGLE
%   1 = 20deg, 2 = 40deg, 3 = 50deg, 4 = 60deg
choice = 2;

%% Set values of relaxed parameters (except gamma, which you chose above)
P_rest = 0;
r_rest = 0.1875;    % 3/16 in
L_rest = 3.4;     % medium one
phi_rest = 0;
T_rest = 0;     % fiber tension
t_rest = 1/16;   % tube wall thickness

params.t_rest = t_rest;   % tube wall thickness

%% Set pressure range, loads, elastic modulus and other parameters
params.Pmin = 0;
params.Pmax = 15;        % maximum pressure tested

% Set the values of the external loads
params.load = [0, 0];        % [F_load, M_load];

% Set values of moduli of elastomer (if not being set specifically for different angles below)
% params.modulus = [500, 200]; % Young's and shear modulus of elastomer: [E, G] (constant modulus)
% params.modulus = [270, 200];
% params.modulus = [3117, 218.5; -440.4, 129.3]; % Linear coefficients for moduli: [c1, c2; c3, c4] --> E = c1*dLnorm + c2, G = c3*dphinorm + c4 (varying modulus)
% params.modulus = [3117, 218.5; -150/0.2, 200]; % Linear coefficients for moduli: [c1, c2; c3, c4] --> E = c1*dLnorm + c2, G = c3*dphinorm + c4 (varying modulus)


%% Measured data points at each fiber angle (no load)

% 20 deg FREE
P_20 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 10, 8.5, 7.5, 5, 6.5, 4, 2.5, 6.5];
dL_20 = P_20*(-1e-3);
dphi_20 = [0.1419999879,...
           0.2387610417, 0.3443185548, 0.469982261, 0.6283185307, 0.7879114375, 1.03044239, 1.043008761, 0.8168140899,...
           0.6710441908, 0.3933274002, 0.5152211952, 0.266407057, 0.1445132621, 0.4611858015];

% 40 deg FREE
P_40 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 11, 8.5, 6.5, 4, 2.5];
dL_40 = [0, 0, 0, 0, 0, 0.004, 0.02, 0.028, 0.04, 0.06, 0.012, 0.012, 0.012, 0.012, 0.012, 0.016, 0.02, 0.028, 0.044, 0.06, 0.06, 0.06, 0.048, 0.036, 0.032, 0, 0, 0];
dphi_40 = [0.3254689989, 0.511451284, 0.6810972873, 0.9060353213, 1.129716718, 1.323238826, 1.672583929, 1.954070631,...
           2.248123703, 2.561026331, 0.3619114737, 0.5943893301, 0.7715751557, 1.005309649, 1.23150432, 1.425026428,...
           1.744212241, 1.964123727, 2.307185645, 2.626371458, 2.096070618, 1.543150311, 1.084477784, 0.6785840132, 0.3782477555];

       
% 60 deg FREE
P_60 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 10, 12, 8.5, 5, 2.5];
dL_60 = [0, 0, 0.012, 0.032, 0.052, 0.064, 0.084, 0.108, 0.132, 0.164, 0, 0, 0.012, 0.024, 0.04, 0.056, 0.084, 0.108, 0.136, 0.144, 0.136, 0.144, 0.104, 0.064, 0.024, 0, 0, 0];
dphi_60 = [0.3154159024, 0.5767964112, 0.8004778081, 1.105840614, 1.306902544, 1.5054512, 1.851026391, 2.084760885,...
           2.481858196, 2.8651325, 0.4033804967, 0.6094689748, 0.879645943, 1.080707873, 1.290566262, 1.577079512,...
           1.888725503, 2.222990962, 2.618831636, 2.84376967, 2.168955568, 2.648990926, 1.734159145, 1.024159205, 0.5026548246];

       
% 50 deg FREE
P_50 = [2.5, 4, 5, 6, 7.5, 8.5, 10, 11, 12.5, 13, 12.5, 13.5, 11, 12, 10, 11, 8.5, 10, 7.5, 9, 6.5, 8, 5.4, 6.5,...
        4, 5, 3.4, 0];
dL_50 = [0, 0, 0.004, 0.008, 0.024, 0.048, 0.06, 0.084, 0.104, 0.148, 0.148, 0.164, 0.14, 0.148, 0.116,...
         0.124, 0.096, 0.1, 0.096, 0.096, 0.076, 0.076, 0.064, 0.064, 0.056, 0.056, 0.052, 0.032];
dphi_50 = [0.3506017401, 0.5981592412, 0.7841415263, 0.981433545, 1.303132633, 1.559486593, 1.901291874, 2.280796267,...
           2.816123655, 3.204424507, 3.183061677, 3.507274038, 2.790990913, 3.204424507, 2.454212181,...
           2.729415697, 2.038265314, 2.408973247, 1.807044094, 2.065911329, 1.501681288, 1.771858257, 1.221451224,...
           1.476548547, 0.9437344331, 1.109610525, 0.7690618816, 0.1445132621];
       
% 40 deg FREE with sandpaper
P_40v2 = [3, 4.3, 5.4, 6.3, 8, 8.8, 10.1, 11.1, 12.2, 13.7, 11, 12.3, 10, 11.2, 8.9, 10.1, 7.8, 9, 7,...
          8, 5.2, 6.6, 4, 5, 3, 0];
dL_40v2 = [0, 0.004, 0.008, 0.012, 0.024, 0.036, 0.052, 0.064, 0.084, 0.112, 0.096, 0.104, 0.084, 0.084, 0.08, 0.08, 0.064, 0.064, 0.056, 0.056, 0.048,...
           0.048, 0.044, 0.044, 0.036, 0.024];
dphi_40v2 = [0.3581415625, 0.5504070329, 0.736389318, 0.9022654101, 1.198831757, 1.394867138, 1.751752064, 2.00559275, 2.333575023, 2.788477639,...
             2.141309553, 2.469291826, 1.922654704, 2.173982116, 1.632371543, 1.928937889, 1.394867138, 1.642424639, 1.21893795, 1.422513154,...
             0.8984954989, 1.112123799, 0.6723008279, 0.8432034682, 0.5215043805, 0.08168140899];
       
       
%% Choose which measured values to use: 1 = 20 deg FREE, 2 = 40 deg FREE, 3 = 60 deg FREE, 4 = 50 deg FREE
%   notes:
%       rcoeff: coefficients of relationship between r and P. r = c1 + c2*P.
%       modulus: Young's (E) and shear (G) moduli: [E, G]

if choice == 1
    gama_rest = deg2rad(20);
    P = P_20;
    dL = dL_20;
    dphi = dphi_20;
    params.rcoeff = [r_rest, 0.0118/2];     % 20 degree FREE
    params.modulus = [1200, 700];
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('data/load20deg_v2.csv',2,46);
elseif choice == 2
    gama_rest = deg2rad(40);
    P = P_40v2;
    dL = dL_40v2;
    dphi = dphi_40v2;
    params.rcoeff = [r_rest, 0.0038/2];     % 40 degree FREE
    params.modulus = [1200, 700];
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('data/load40deg_v3.csv',2,44);
elseif choice == 3
    gama_rest = deg2rad(60);
    P = P_60;
    dL = dL_60;
    dphi = dphi_60;
    params.rcoeff = [r_rest, 0.0038/2];     % 60 degree FREE
    params.modulus = [1200, 700];
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('data/load60deg_v2.csv',2,50);
elseif choice == 4
    gama_rest = deg2rad(50);
    P = P_50;
    dL = dL_50;
    dphi = dphi_50;
    params.rcoeff = [r_rest, 0.0038/2];     % 50 degree FREE
    params.modulus = [1200, 700];
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('data/load50deg_v2.csv',2,48);
end

%% Measured data at each fiber angle, with loads. Separate data by load

index1 = 0;
index2 = 0;

for index = 1:length(loadbolts)
    if loadbolts(index) == 1
        index1 = index1 + 1;
        Pvphi_1bolt(index1,:) = [loadP(index), loadphi(index)];
    elseif loadbolts(index) == 2
        index2 = index2 + 1;
        Pvphi_2bolt(index2,:) = [loadP(index), loadphi(index)];
    end
end

%% Compile resting parameters into a single vector, x_rest

params.theta_rest = -L_rest*tan(gama_rest)/(r_rest);
params.x_rest = [P_rest, gama_rest, r_rest, L_rest, phi_rest, T_rest];

