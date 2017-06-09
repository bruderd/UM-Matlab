% /SimpleParallel/setParams.m
%   Set resting values (when p = 0) for state (xrest), and other FREE parameters.
clear
clc

params = struct;

%% USER DEFINED FREE PARAMETERS
params.Gama = deg2rad(40);  % relaxed fiber angle (rad)
params.R = 3/16;    % relaxed FREE radius (in)
params.L = 3.4; % relaxed FREE length (in)
params.B = params.L/cos(params.Gama);   % fiber length (in)
params.N = -params.L/(2*pi*params.R) * tan(params.Gama); % total fiber windings in revolutions (when relaxed)

% params.kelast = [-4373.98, -128.76];   % spring constants for the elastomer
params.kelast = [-117.4588, -1.3208];   % spring constants for the elastomer

params.load = [0, 0];   % loads on FREE [Fload, Mload]

params.Ptest = 10; % pressure at which plots generated (psi)



