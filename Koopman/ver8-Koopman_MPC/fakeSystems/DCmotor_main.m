% DCmotor_main.m
%   Main script for simulating a DC motor system. This is pulled from the
%   Korda and Mezic Koopman MPC paper.

%% set double system parameters
params = struct;

% system parameters
params.name = 'DCmotor_100s';
params.n = 2;   % number of states
params.p = 1;   % number of inputs

% physical parameters
params.La      = 0.314; 
params.Ra      = 12.345;
params.km      = 0.253;
params.J       = 0.00441;
params.B       = 0.00732; 
params.tau1    = 1.47; 
params.ua      = 60; 


% simulation parameters
params.Ts = 0.01;
params.tf = 100;
params.x0 = 1 * (2*rand([2,1]) - 1);

% step/ramp input parameters
params.steplen = 0.5;  % duration of each step
params.steps = 4 * (2*rand(2000, params.p) - 1);   % random list of 1000 step inputs between -10 and 10
params.tconstant = 1.5;

%% simulate system
data = gen_fakeData( params.name, @DCmotor_vf, @DCmotor_input, params );