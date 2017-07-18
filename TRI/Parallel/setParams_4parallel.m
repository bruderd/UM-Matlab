% /SimpleParallel/setParams_4parallel.m
%   Set resting values (when p = 0) for state (xrest), and other FREE parameters.
clear
clc

params = struct;

%% USER DEFINED FREE PARAMETERS

params.numFREEs = 4;    % total number of FREEs
params.Lspine = 20e-2;  % length of central spine

% Parameters for the first FREE, designated by postscript "1"
params.Gama1 = 0.75;  % relaxed fiber angle (rad)
params.R1 = 0.5e-2;    % relaxed FREE radius (m)
params.L1 = 20e-2; % relaxed FREE length (m)
params.kelast1 = [-4e3, -4e-2];   % spring constants for the elastomer
% params.kelast1 = [0, 0];   % spring constants for the elastomer
params.Ptest1 = 70e3; % test pressure (Pa). Not needed for most functions.
params.attach1 = [2e-2, 2e-2]'; %coordinates of attachment point of FREE

% Parameters for the second FREE, designated by postscript "2"
params.Gama2 = -0.75;  % relaxed fiber angle (rad)
params.R2 = 0.5e-2;    % relaxed FREE radius (m)
params.L2 = 20e-2; % relaxed FREE length (m)
params.kelast2 = [-4e3, -4e-2];   % spring constants for the elastomer
% params.kelast2 = [0, 0];   % spring constants for the elastomer
params.Ptest2 = 35e3; % test pressure (Pa). Not needed for most functions.
params.attach2 = [2e-2, -2e-2]'; %coordinates of attachment point of FREE

% Parameters for the third FREE, designated by postscript "3"
params.Gama3 = 0.75;  % relaxed fiber angle (rad)
params.R3 = 0.5e-2;    % relaxed FREE radius (m)
params.L3 = 20e-2; % relaxed FREE length (m)
params.kelast3 = [-4e3, -4e-2];   % spring constants for the elastomer
% params.kelast3 = [0, 0];   % spring constants for the elastomer
params.Ptest3 = 35e3; % test pressure (Pa). Not needed for most functions.
params.attach3 = [-2e-2, -2e-2]'; %coordinates of attachment point of FREE

% Parameters for the fourth FREE, designated by postscript "4"
params.Gama4 = -0.75;  % relaxed fiber angle (rad)
params.R4 = 0.5e-2;    % relaxed FREE radius (m)
params.L4 = 20e-2; % relaxed FREE length (m)
params.kelast4 = [-4e3, -4e-2];   % spring constants for the elastomer
% params.kelast4 = [0, 0];   % spring constants for the elastomer
params.Ptest4 = 35e3; % test pressure (Pa). Not needed for most functions.
params.attach4 = [-2e-2, 2e-2]'; %coordinates of attachment point of FREE

% Parameters shared by all FREEs
params.load = [0; 0];   % loads on FREE [Fload, Mload] (N)

%% USER DEFINED TEST PARAMETERS

% range of pressures to iterate over
params.Prange1 = [0, 200e3];
params.Prange2 = [0, 200e3];
params.Prange3 = [0, 200e3];
params.Prange4 = [0, 200e3];

% number of steps to take along each FREEs pressure;
params.steps1 = 20;
params.steps2 = 20;
params.steps3 = 20;
params.steps4 = 20;

%% USER DEFINED PLOTTING PARAMETERS 

params.thickness = 0.01;    % thickness of the top and bottom blocks
params.width = 0.04;        % width of the top and bottom blocks


%% Dependent parameters (do not edit below this line)

% FREE 1
params.B1 = abs(params.L1/cos(params.Gama1));   % fiber length (must be positive)
params.N1 = -params.L1/(2*pi*params.R1) * tan(params.Gama1); % total fiber windings in revolutions (when relaxed)

% FREE 2
params.B2 = abs(params.L2/cos(params.Gama2));   % fiber length (must be positive)
params.N2 = -params.L2/(2*pi*params.R2) * tan(params.Gama2); % total fiber windings in revolutions (when relaxed)

% FREE 3
params.B3 = abs(params.L3/cos(params.Gama3));   % fiber length (must be positive)
params.N3 = -params.L3/(2*pi*params.R3) * tan(params.Gama3); % total fiber windings in revolutions (when relaxed)

% FREE 4
params.B4 = abs(params.L4/cos(params.Gama4));   % fiber length (must be positive)
params.N4 = -params.L4/(2*pi*params.R4) * tan(params.Gama4); % total fiber windings in revolutions (when relaxed)


