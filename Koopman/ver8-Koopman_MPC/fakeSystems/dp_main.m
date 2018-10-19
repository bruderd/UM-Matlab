% dp_main.m
%   Main script for simulating a double pendulum system

%% set double pendulum parameters
params = struct;

% physical parameters
params.phi1                = pi/6; % (pi/2)*rand - pi/4;
params.dtphi1              = 0;
params.phi2                = -pi/4.5; % (pi/2)*rand - pi/4;
params.dtphi2              = 0;
params.g                   = 9.81; 
params.m1                  = 1; 
params.m2                  = 1; 
params.l1                  = 1; 
params.l2                  = 1;

% simulation parameters
params.Ts = 0.02;
params.tf = 100;
params.x0 = [params.phi1, phi2, dtphi1, dtphi2]';

%% simulate double pendulum