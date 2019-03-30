% SMD_main.m
%   Main script for simulating a DC motor system. This is pulled from the
%   Korda and Mezic Koopman MPC paper.

%% set double system parameters
params = struct;

% system parameters
params.name = 'SMD_10s';
params.n = 2;   % number of states
params.p = 1;   % number of inputs

% physical parameters
params.m = 1;
params.b = 1;
params.k = 1;


% simulation parameters
params.Ts = 0.01;
params.tf = 10;
params.x0 = 1 * (2*rand([params.n,1]) - 1);

% step/ramp input parameters
params.steplen = 0.5;  % duration of each step
params.steps = 1 * (2*rand(2000, params.p) - 1);   % random list of 1000 step inputs between -1 and 1
params.tconstant = 1.5;

% other params
params.saveon = true;

%% simulate system
data = gen_fakeData( params.name, @SMD_vf, @SMD_input, params );

%% save function that evaluates dynamics with these system parameters
x = sym('x' , [params.n , 1]);
u = sym('u' , [params.p , 1]);
symbolic_dynamics = SMD_vf(x, u, params);