%% ODE15i Simulation

tspan = [0, 20];
x0 = [0.001, deg2rad(40), 3/16, 5.68]';
xdot0 = [0 0 0 0]';
u = 10;

[t, y] = ode15i(@(t, x, dxdt)vf(x,u,dxdt,params), tspan, x0, xdot0);

