% /VSA/setParams.m
%   Set resting values (when P = 0) for state (xrest), and other McKibbon parameters.
clear
clc

params = struct;

%% USER DEFINED FREE PARAMETERS

params.Gama = 0.75;  % relaxed fiber angle (rad)
params.R = 1.5e-2;    % relaxed FREE radius (m)
params.L = 10e-2; % relaxed FREE length (m)

params.kelast = -4e3;   % spring constants for the elastomer

params.Ptest = 69e3; % pressure at which plots generated (Pa)

params.I = 1;   % moment of inertia of rotating joint

params.rho = 2.5e-2;  % radius of rotating joint (m)

%% User defined test parameters

% range of pressures to iterate over
params.Prange1 = [0, 200e3];
params.Prange2 = [0, 200e3];

% number of steps to take along each FREEs pressure;
params.steps1 = 50;
params.steps2 = 50;

%% User defined control parameters (for simple control demos)

params.theta_des = pi/4;    % desired rotation (rad)
params.thetadot_des = 0;    % desired angular velocity (rad/s)

params.K_des = -20;     % desired stiffness (N/m)

params.C = [1,0];  % proportional gain constants


%% Dependent parameters

params.B = params.L/cos(params.Gama);   % fiber length
params.N = -params.L/(2*pi*params.R) * tan(params.Gama); % total fiber windings in revolutions (when relaxed)

params.x_des = [params.theta_des; params.thetadot_des]; % desired state