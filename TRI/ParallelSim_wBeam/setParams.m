%setParams.m
%   Set resting values (when P = 0) for state (xrest), and other FREE parameters.
clear
clc
params = struct;

%% USER DEFINED FREE PARAMETERS

params.numFREEs = 4;    % total number of FREEs
params.Lspine = 27.94e-2;  % length of central spine

params.Gama = [0.75, -0.75, 0.75, -0.75]';        % relaxed fiber angle (rad)
% params.Gama = deg2rad([42.9, -42.9, 42.9, -42.9])';        % relaxed fiber angle (rad)
params.R = [0.5e-2, 0.5e-2, 0.5e-2, 0.5e-2]';    % relaxed FREE radius (m)
params.L = [params.Lspine, params.Lspine, params.Lspine, params.Lspine]'; % relaxed FREE lenght (m)
params.kelast_s = 1e-2 * [-4e3, -4e3, -4e3, -4e3]';     % spring constants for the elastomer in extension
params.kelast_w = 1e0 * [-4e-2, -4e-2, -4e-2, -4e-2]'; % spring constants for the elastomer in twist
params.xattach = 1*[2e-2, -2e-2, -2e-2, 2e-2]';      % x-coordinate of attachment point of FREE
params.yattach = 1*[2e-2, 2e-2, -2e-2, -2e-2]';      % y-coordinate of attachment point of FREE

%% USER DEFINED END EFFECTOR PARAMETERS
params.m = 0.01;         % end effector mass (kg)
params.effdim = [4e-2, 4e-2, 2e-2];     % [lenght, width, height] (m) ... assuming end effector is rectangular prism. 
params.I = params.m/12 * [params.effdim(1)^2 + params.effdim(3)^2; params.effdim(2)^2 + params.effdim(3)^2; params.effdim(1)^2 + params.effdim(2)^2]; % end effector rotational moments of inertia [Mx, My, Mz]'
% params.I = [0.1; 0.1; 0.01]; % end effector rotational moments of inertia [Mx, My, Mz]'
params.g = 9.81;        % acceleration due to gravity (m/s^2)
params.damp = 0.5*[1.3e0; 1.3e0; 1.25*8e-2];    % damping in each direction [damppsi, damptheta, dampphi]'

%% USER DEFINED SPINE/BEAM PARAMETERS
params.dbeam = 0.5e-3;    % beam diameter (m)
params.Ebeam = 200e9;  % beam Young's modulus (Pa)
params.Ibeam = pi*params.dbeam^4 / 64;

%% USER DEFINED TEST PARAMETERS

% Initial conditions (could be made more generic in the future)
params.x0 = [0, 0, 0, 0, 0, 0]';
params.xdot0 = [params.x0(4), params.x0(5), params.x0(6), 0, 0, 0]';
params.u0 = [1000, 0, 1000, 0]';

params.T = 1;   %final time
params.N = 50; %number of steps
params.dt = params.T/params.N;    %size of one time step
params.n = length(params.x0);  %dimension of state vector x
params.m = length(params.u0);  %dimension of state vector u

% Maximum pressure for the FREEs
params.Pmax = 400e3;    % (Pa)

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

params.thickness = 0.02;    % thickness of the top and bottom blocks
params.width = 0.04;        % width of the top and bottom blocks

%% Dependent parameters (do not edit below this line)

params.B = abs(params.L ./ cos(params.Gama));   % fiber length (must be positive))
params.Nf = -params.L ./ (2*pi*params.R) .* tan(params.Gama); % total fiber windings in revolutions (when relaxed)
