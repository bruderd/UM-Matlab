% setParams.m
%   Set initial values for state (x0) and input (u0), and other solver
%   parameters.
clear
clc

params = struct;

params.x0 = [1 2]';    %state vector initial condition (any size)
params.u0 = [3]';        %input vector initial condition (any size)

params.T = 20;   %final time
params.N = 100; %number of steps
params.dt = params.T/params.N;    %size of one time step
params.n = length(params.x0);  %dimension of state vector x
params.m = length(params.u0);  %dimension of state vector u

params.vf = @(x, u) vf(x, u);

% creates functions to evaluate the dynamics at any point
dynamics_symbolic(params);




