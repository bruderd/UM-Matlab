% /OneFREE/setParams.m
%   Set resting values (when P = 0) for state (xrest), and other FREE parameters.
clear
clc

params = struct;

%% USER DEFINED FREE PARAMETERS

params.Gama = 0.75;  % relaxed fiber angle (rad)
params.R = 0.5e-2;    % relaxed FREE radius (m)
params.L = 10e-2; % relaxed FREE length (m)

params.kelast = [-4e3, -4e-2];   % spring constants for the elastomer

% params.load = [0, 0];   % loads on FREE [Fload, Mload] (N)
params.load = [20, -0.1];   % loads on FREE [Fload, Mload] (N)

params.Ptest = 69e3; % pressure at which plots generated (Pa)

%% Dependent parameters

params.B = params.L/cos(params.Gama);   % fiber length
params.N = -params.L/(2*pi*params.R) * tan(params.Gama); % total fiber windings in revolutions (when relaxed)




