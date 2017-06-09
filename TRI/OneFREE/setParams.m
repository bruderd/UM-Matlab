% /SimpleParallel/setParams.m
%   Set resting values (when p = 0) for state (xrest), and other FREE parameters.
clear
clc

params = struct;

%% USER DEFINED FREE PARAMETERS

params.Gama = 0.75;  % relaxed fiber angle (rad)
params.R = 0.5e-2;    % relaxed FREE radius (m)
params.L = 10e-2; % relaxed FREE length (m)
params.B = params.L/cos(params.Gama);   % fiber length
params.N = -params.L/(2*pi*params.R) * tan(params.Gama); % total fiber windings in revolutions (when relaxed)

params.kelast = [-4e3, -4e-2];   % spring constants for the elastomer

params.load = [0, 0];   % loads on FREE [Fload, Mload]

params.Ptest = 69e3; % pressure at which plots generated (Pa)



%% Parameters for an actual FREE (roughly)
% params.Gama = deg2rad(45);  % relaxed fiber angle (rad)
% params.R = 4.7625e-3;    % relaxed FREE radius (m)
% params.L = 86.36e-3; % relaxed FREE length (m)
% params.B = params.L/cos(params.Gama);   % fiber length
% params.N = -params.L/(2*pi*params.R) * tan(params.Gama); % total fiber windings in revolutions (when relaxed)
% 
% % params.kelast = [-4373.98, -128.76];   % spring constants for the elastomer
% params.kelast = [-4000, -0.05];   % spring constants for the elastomer
% 
% params.load = [0, 0];   % loads on FREE [Fload, Mload]
% 
% params.Ptest = 69e3; % pressure at which plots generated (Pa)



