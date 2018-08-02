% Simply call 
%
%       >> double_pendulum_init
%
% to run the double pendulum simulation with the below parameters. This
% script calls double_pendulum.
%
%   ---------------------------------------------------------------------

params = struct;

params.phi1                = pi;
params.dtphi1              = 0;
params.phi2                = pi;
params.dtphi2              = 0.1;
params.g                   = 9.81; 
params.m1                  = 1; 
params.m2                  = 1; 
params.l1                  = 1; 
params.l2                  = 1;
params.duration            = 10;
params.fps                 = 30;
params.movie               = true;

clc; figure;

% interval=[0, duration];
% ivp=[phi1; dtphi1; phi2; dtphi2; g; m1; m2; l1; l2];

double_pendulum(params);