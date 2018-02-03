function f = staticBalance(x, u, params)
%staticBalance: force balance function of manipulator
%   Detailed explanation goes here

%% Define local names of global parameters

p = params.p;       % total number of modules in manipulator
n = params.n;       % number of actuators in each module (a vector)
m = params.m;       % masses of the blocks
I = params.I;       % moment of inertia matrices of blocks
L = params.L;       % lengths of the modules

%% Caclulate the manipulator forces (zeta0)

P = u;
% P = u.^2;      % input is vector of pressures

x0 = x_orient2x0_orient(x, params);
qdot = zeros(2*sum(n),1);
xdot = zeros(3*p,1);
x0dot = zeros(3*p,1);

[Z, Vdot] = actuatorBlock(P, qdot, x, params);
[zeta, ~] = moduleBlock(Z, xdot, x, params);
[zeta0, ~] = manipulatorBlock(zeta, x0dot, x0, params);


%% Force balance equation
f = zeta;


end

