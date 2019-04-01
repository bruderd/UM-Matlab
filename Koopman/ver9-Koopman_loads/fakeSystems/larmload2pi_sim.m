% armloadv2_sim.m
%   Main script for simulating a 1dof arm with load attached system

%% set double pendulum parameters
params = struct;

% system parameters
params.name = 'armloadv2_100s_load0';
params.n = 3;   % number of states (load is the third state)
params.p = 1;   % number of inputs

% physical parameters
params.M = 1;
params.K = 20;
params.L = 1;
params.D = -10;
params.g = 9.81;
params.load = 8;    % LOAD ATTACHED TO END EFFECTOR

% simulation parameters
params.Ts = 0.02;
params.tf = 100;
params.x0 = [ 0, 0, params.load]';

% step/ramp input parameters
params.steplen = 1;  % duration of each step
params.steps = pi * rand(2000, params.p);   % random list of 1000 step inputs between 0 and pi
% params.steps = pi/2 * ones(2000,params.p);
params.tconstant = 2;

% other params
params.saveon = false;

%% simulate double pendulum
data = gen_fakeData( params.name, @armload_vf, @armload_input, params );

%% save function that evaluates dynamics with these system parameters
x = sym('x' , [params.n , 1]);
u = sym('u' , [params.p , 1]);
symbolic_dynamics = armload_vf(x, u, params);