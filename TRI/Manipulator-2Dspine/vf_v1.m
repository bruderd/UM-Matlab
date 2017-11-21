function [f, dfdX, dfdu, dfdXdot] = vf_v1(X0, u, X0dot, params)
%vf_v1: Dynamics soft robotic manipulator
%   Detailed explanation goes here

%% Define local names of global parameters

p = params.p;       % total number of modules in manipulator

%% Caclulate the manipulator forces (zeta0)
alpha = X0(1:p);
alphadot = X0(p+1 : 2*p);
P = u;      % input is vector of pressures
gama = zeros(3*p,1);  % no external loads

[~, xdot] = manipulatorBlock(gama, alphadot, alpha, params);
[~, qdot] = moduleBlock(0, alphadot, alpha, params);
[Z, Vdot] = actuatorBlock(P, qdot, alpha, params);
[zeta, qdot] = moduleBlock(Z, alphadot, alpha, params);
[tau, xdot] = manipulatorBlock(gama, alphadot, alpha, params);


%% Equations of Motion, f(xddot, xdot, x, u) = 0
f(1:p, 1) = X0dot(1:p) - X0(p+1 : 2*p);
f(p+1 : 2*p, 1) = EOM(X0, X0dot, zeta, gama) + 1e-4 * ones(p,1);

%% Gradients (not needed for now)
dfdX = [0];
dfdu = [0];
dfdXdot = [0];

end

