% arm_sim.m
%   Main script for simulating a 1dof arm with load attached system

%% set double pendulum parameters
params = struct;

% system parameters
params.name = 'arm_1000s_load2';
params.n = 2;   % number of states
params.p = 1;   % number of inputs

% physical parameters
params.M = 1;
params.K = 20;
params.L = 1;
params.D = -10;
params.g = 9.81;
params.load = 2;    % load attached to the end effector

% simulation parameters
params.Ts = 0.02;
params.tf = 1000;
params.x0 = [ 0, 0]';

% step/ramp input parameters
params.steplen = 1;  % duration of each step
params.steps = pi * rand(2000, params.p);   % random list of 1000 step inputs between -10 and 10
% params.steps = pi/2 * ones(2000,params.p);
params.tconstant = 2;

% other params
params.saveon = true;

%% simulate double pendulum
data = gen_fakeData( params.name, @arm_vf, @arm_input, params );

%% save function that evaluates dynamics with these system parameters
x = sym('x' , [params.n , 1]);
u = sym('u' , [params.p , 1]);
symbolic_dynamics = arm_vf(x, u, params);