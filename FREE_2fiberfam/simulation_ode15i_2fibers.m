%% ODE15i Simulation - 2 fibers
%   Use this to check if the dynamics are broken, or if it's something else

tspan = [0, 10];
x_rest = params.x_rest;     % resting point
x0 = params.x0;       % initial point 
x0 = [1.6667    0.9731    0.5344    0.1895    5.5099    0.0866];    % inital point. Taken from optimization, rather than resting condition
xdot0 = [0 0 0 0 0 0]';
u = 10;

% options = odeset('abstol', 1e-6, 'reltol', 1e-6, 'NonNegative', 1);

fixed_x0 = [1 1 1 1 1 1];
fixed_xdot0 = [0 0 0 0 0 0];

[x0_new,xdot0_new] = decic(@(t, x, dxdt)vf2(x,u,dxdt,params),0,x0,fixed_x0,xdot0,fixed_xdot0);

[t, y] = ode15i(@(t, x, dxdt)vf2(x,u,dxdt,params), tspan, x0_new, xdot0_new);


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


