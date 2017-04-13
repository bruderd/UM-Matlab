%% ODE15i Simulation - 2 fibers
%   Use this to check if the dynamics are broken, or if it's something else
%   Version 2: Used for state vector [P, gama, betta, r, L, phi, T_gama, T_betta]'

tspan = [0, 10];
x_rest = params.x_rest;     % resting point
x0 = params.x0;       % initial point 
x0 = [-0.000233112795968   0.610821992664272  -0.349036299887671   0.187489686900412   5.000081964340222   0.000000909290820  -0.000107952274380   -0.000233337373754];
xdot0 = [0 0 0 0 0 0 0 0]';
u = 10;

% options = odeset('abstol', 1e-6, 'reltol', 1e-6, 'NonNegative', 1);

fixed_x0 = [1 1 1 1 1 1 0 0];
fixed_xdot0 = [0 0 0 0 0 0 0 0];

[x0_new,xdot0_new] = decic(@(t, x, dxdt)vf(x,u,dxdt,params),0,x0,fixed_x0,xdot0,fixed_xdot0);

[t, y] = ode15i(@(t, x, dxdt)vf(x,u,dxdt,params), tspan, x0_new, xdot0_new);


% Plot the results
figure
set(gcf,'numbertitle','off','name','ode15i Results') % See the help for GCF

subplot(2,3,1)       
plot(t, y(:,1))
title('Pressure (psi)')

subplot(2,3,2)       
plot(t, rad2deg(y(:,2)))
title('gamma (deg)')

subplot(2,3,3)       
plot(t, rad2deg(y(:,3)))
title('beta (deg)')

subplot(2,3,4)       
plot(t, y(:,4))
title('radius (in)')

subplot(2,3,5)       
plot(t, y(:,5))
title('Length (in)')

subplot(2,3,6)       
plot(t, rad2deg(y(:,6)))
title('phi (deg)')


