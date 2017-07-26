%setParams.m
%   Set resting values (when P = 0) for state (xrest), and other FREE parameters.
clear
clc
params = struct;

%% USER DEFINED FREE PARAMETERS

params.numFREEs = 4;    % total number of FREEs
params.Lspine = 20e-2;  % length of central spine

params.Gama = [0.75, -0.75, 0.75, -0.75]';        % relaxed fiber angle (rad)
params.R = [0.5e-2, 0.5e-2, 0.5e-2, 0.5e-2]';    % relaxed FREE radius (m)
params.L = [params.Lspine, params.Lspine, params.Lspine, params.Lspine]'; % relaxed FREE lenght (m)
params.kelast_s = [-4e3, -4e3, -4e3, -4e3]';     % spring constants for the elastomer in extension
params.kelast_w = [-4e-2, -4e-2, -4e-2, -4e-2]'; % spring constants for the elastomer in twist
params.xattach = [2e-2, -2e-2, -2e-2, 2e-2]';      % x-coordinate of attachment point of FREE
params.yattach = [2e-2, 2e-2, -2e-2, -2e-2]';      % y-coordinate of attachment point of FREE

%% USER DEFINED END EFFECTOR PARAMETERS
params.m = 0.1;         % end effector mass (kg)
params.I = [0.1; 0.1; 0.01]; % end effector rotational moments of inertia [Mx, My, Mz]'
params.g = 9.81;        % acceleration due to gravity (m/s^2)
params.damp = [5e-1; 5e-1; 5e-1];    % damping in each direction [damppsi, damptheta, dampphi]'

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

params.B = abs(params.L ./ cos(params.Gama));   % fiber length (must be positive))
params.Nf = -params.L ./ (2*pi*params.R) .* tan(params.Gama); % total fiber windings in revolutions (when relaxed)
