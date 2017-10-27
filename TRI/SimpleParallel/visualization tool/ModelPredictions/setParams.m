% setParams.m
%   Set resting values (when P = 0) for state (xrest), and other FREE parameters.
clear
clc

params = struct;

%% USER DEFINED FREE PARAMETERS

params.GamaDeg = 40;  % relaxed fiber angle (deg)
params.R = 6.25e-3;    % relaxed FREE radius (m)
params.L = 97.4852e-3; % relaxed FREE length (m)

% params.kelast = [-4e1, -4e-6];   % spring constants for the elastomer
params.kelast = [-1e-1, -1e-3];   % spring constants for the elastomer

params.load = [0, 0];   % loads on FREE [Fload, Mload] (N)
% params.load = [20, -0.1];   % loads on FREE [Fload, Mload] (N)

params.Ptest = [21, 35, 61]; % pressure(s) at which plots generated (kPa)

% How many points to use for s and w
params.res = 20;

%% Optional parameters

params.srange = [-5e-3 ,5e-3];     % [smin, smax] (m)
params.wrange = [-120 ,120];     % [wmin, wmax] (deg)

params.sinc = 1e-3;       % s increment size (m)
params.winc = 10;       % w increment size (deg)

%% Dependent parameters

params.GamaRad = deg2rad(params.GamaDeg);
params.B = params.L/cos(params.GamaRad);   % fiber length
params.N = -params.L/(2*pi*params.R) * tan(params.GamaRad); % total fiber windings in revolutions (when relaxed)




