%% ODE15i Simulation - 2 fibers
%   Use this to check if the dynamics are broken, or if it's something else

tspan = [0, 100];
x_rest = params.x_rest;     % resting point
x0 = params.x0;       % initial point 
xdot0 = [0 0 0 0 0 0]';
u = 10;

% options = odeset('abstol', 1e-6, 'reltol', 1e-6, 'NonNegative', 1);

fixed_x0 = [1 1 1 1 1 1];
fixed_xdot0 = [0 0 0 0 0 0];

[x0_new,xdot0_new] = decic(@(t, x, dxdt)vf2(x,u,dxdt,params),0,x0,fixed_x0,xdot0,fixed_xdot0);

[t, y] = ode15i(@(t, x, dxdt)vf2(x,u,dxdt,params), tspan, x0_new, xdot0_new);