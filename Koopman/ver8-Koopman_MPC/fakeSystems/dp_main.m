% dp_main.m
%   Main script for simulating a double pendulum system

%% set double pendulum parameters
params = struct;

% system parameters
params.name = 'doublePendulum_100s';
params.n = 4;   % number of states
params.p = 1;   % number of inputs

% physical parameters
params.phi1                = 0; % (pi/2)*rand - pi/4;
params.dtphi1              = 0;
params.phi2                = pi/3; % (pi/2)*rand - pi/4;
params.dtphi2              = 0;
params.g                   = 9.81; 
params.m1                  = 1; 
params.m2                  = 1; 
params.l1                  = 1; 
params.l2                  = 1;

% simulation parameters
params.Ts = 0.02;
params.tf = 100;
params.x0 = [params.phi1, params.phi2, params.dtphi1, params.dtphi2]';

% step/ramp input parameters
params.steplen = 1;  % duration of each step
params.steps = 16 * rand(2000, params.p) - 8;   % random list of 1000 step inputs between -10 and 10
params.tconstant = 1.5;

%% simulate double pendulum
data = gen_fakeData( params.name, @dp_vf, @dp_input, params );

%% save function that evaluates dynamics with these system parameters
x = sym('x' , [params.n , 1]);
u = sym('u' , [params.p , 1]);
symbolic_dynamics = dp_vf(x, u, params);