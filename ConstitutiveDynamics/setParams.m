% setParams.m
%   Set initial values for state (x0) and input (u0), and other solver
%   parameters.
clear
clc

params = struct;

params.x0 = [0, 3/16, 0.1, 0.012, 0]';    %state vector initial condition (any size)
params.u0 = [0]';        %input vector initial condition (any size)

params.T = 1;   %final time
params.N = 10; %number of steps
params.dt = params.T/params.N;    %size of one time step
params.n = length(params.x0);  %dimension of state vector x
params.m = length(params.u0);  %dimension of state vector u

%% OTHER USER DEFINED CONSTANTS
params.C = [1,1];   % material constants for neo-hookean model [C1, C2]
params.Gama = deg2rad(40);  % relaxed fiber angle (rad)
params.R = [0.01, 0.012];    % relaxed FREE radius [Ri, Ro] (m)
params.L = 0.1; % relaxed FREE length (m)
params.load = [0, 0];   % loads on FREE [Fload, Mload]


%% 
params.vf = @(x, u) vf(x, u);

% creates functions to evaluate the dynamics at any point
dynamics_symbolic(params);

% tests the gradients of the dynamics at x0, u0.
[grad, num_grad] = run_test_gradient(params);




