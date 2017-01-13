%% ODE15i Simulation
%   Use this to check of the dynamics are broken, or if it's something else

tspan = [0, 10];
x_rest = [0.001, deg2rad(40), 3/16, 5.68]';     % resting point
x0 = [0.001, deg2rad(40), 3/16, 5.68]';       % initial point
% x0 = [9.99, 0.67183, 0.21476, 5.8307]';         % ss value taken from ICRA
% x0 = [0.06	0.71892	0.18821	5.6084]';
% xdot0 = [3	-0.017	0.007	0.05]';             % ss value calculated from ICRA 
xdot0 = [0 0 0 0]';
u = 10;

options = odeset('abstol', 1e-2, 'reltol', 1e-1);

fixed_x0 = [1 1 1 1];
fixed_xdot0 = [0 0 0 0];

[x0_new,xdot0_new] = decic(@(t, x, dxdt)vf(x,u,dxdt,params),0,x0,fixed_x0,xdot0,fixed_xdot0);

[t, y] = ode15i(@(t, x, dxdt)vf(x,u,dxdt,params), tspan, x0_new, xdot0_new);

