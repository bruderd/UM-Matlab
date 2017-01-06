%% ODE15i Simulation
%   Use this to check of the dynamics are broken, or if it's something else

tspan = [0, 20];
x0 = [0.001, deg2rad(40), 3/16, 5.68]';       % resting point
% x0 = [9.99, 0.67183, 0.21476, 5.8307]';         % ss value taken from ICRA
xdot0 = [0 0 0 0]';
u = 10;

[t, y] = ode15i(@(t, x, dxdt)vf(x,u,dxdt,params), tspan, x0, xdot0);

