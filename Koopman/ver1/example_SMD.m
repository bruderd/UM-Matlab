% example: Spring-mass-damper system
clear

%% set parameters
n = 2;
maxDegree = 3;

m = 1;
b = 1;
k = 1;

p = def_polyLift(n, maxDegree);

Ts = 0.1;
numTrials = 3;
x0max = [1, 0]';
tf = 10;

%% Simulate and find Koopman operator from "measurements"

[x,y] = get_snapshotPairs_SMD(Ts, numTrials, x0max, tf, m, b, k);

U = get_Koopman(x,y);