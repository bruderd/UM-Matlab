% call_plotPos_array
%   This is just a script to call plotPos_array.m so that you don't have to
%   do it from the command line.
%
%   Must run setParams_4parallel.m before this script

% Set the value of the pressures you want in each FREE for each of the 4
% plots
P1 = [1e6, 1e6, 0, 0]';
P2 = [0, 0, 1e6, 1e6]';
P3 = [1e6, 0, 1e6, 0]';
P4 = [1e4, 1e6, 1e5, 0]';

plotPos_array(P1, P2, P3, P4, params);

