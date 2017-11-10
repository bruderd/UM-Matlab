function [f, dfdx, dfdu, dfdxdot] = vf_v1(x, u, xdot, params)
%vf_v1: Dynamics soft robotic manipulator
%   Detailed explanation goes here

%% Define local names of global parameters

p = params.p;       % total number of modules in manipulator





%% Equations of Motion, f(xddot, xdot, x, u) = 0
f(1:6*p, 1) = xdot(1:6*p) - poop;
f(6*p+1 : 2*(6*p), 1) = EOM(x, xdot, zeta0, m, I);



end

