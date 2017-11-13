function [f, dfdX, dfdu, dfdXdot] = vf_v1(X0, u, X0dot, params)
%vf_v1: Dynamics soft robotic manipulator
%   Detailed explanation goes here

%% Define local names of global parameters

p = params.p;       % total number of modules in manipulator
m = params.m;       % masses of the blocks
I = params.I;       % moment of inertia matrices of blocks

%% Caclulate the manipulator forces (zeta0)
x0 = X0(1:6*p);
x0dot = X0(6*p+1 : 2*(6*p));
x = x02x(x0, params);   % convert x0 to local coordinates x
P = u;      % input is vector of pressures

[~, xdot] = manipulatorBlock(0, x0dot, x0, params);
[~, qdot] = moduleBlock(0, xdot, x, params);
[Z, Vdot] = actuatorBlock(P, qdot, x, params);
[zeta, qdot] = moduleBlock(Z, xdot, x, params);
[zeta0, xdot] = manipulatorBlock(zeta, x0dot, x0, params);


%% Equations of Motion, f(xddot, xdot, x, u) = 0
f(1:6*p, 1) = X0dot(1:6*p) - X0(6*p+1 : 2*(6*p));
f(6*p+1 : 2*(6*p), 1) = EOM(X0, X0dot, zeta0, m, I);

%% Gradients (not needed for now)
dfdX = [0];
dfdu = [0];
dfdXdot = [0];

end

