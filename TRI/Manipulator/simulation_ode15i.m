% simulation.m

clear
clc

%% Set system parameter values

p = 1;      % total number of modules
n = [4];     % numer of actuators in each module (a vector)

% module parameters: module = [L, block density, EI of spine], size(module) = p x 3
L = [0.1]';    % length of each module (m)
density = [1e3]';      % density of each module block (kg/m^3)
EI = [0]';       % effective stiffness/inertia of spine material (currently is not factored in since spine is considered forceless)
module = [L, density, EI];

% actuator paramters: free = [Gama, R, xattach, yattach], size(free) = sum(n) x 4
Gama = deg2rad([20, -20, 20, -20]');  % fiber angle (rad)
R = 1 * ones(sum(n), 1) * 1e-2;     % internal radius (m)
xattach = [2, -2, -2, 2]' * 1e-2;
yattach = [2, 2, -2, -2]' * 1e-2;
free = [Gama, R, xattach, yattach];

% simulation parameters (might want to put some simulation parameters here at some point...)
% tspan = [0, 10];
% x0_0 = zeros(6*p, 1);
% x0dot_0 = zeros(6*p, 1);
sim = 0;

% create params struct to be passed to many functions
params = setParams(p, n , module, free, sim);
setJacobians(params);   % derive the appropriate Jacobian matrices
setEOM(params);     % derive the equations of motion

%% ODE15i Simulation

tspan = [0, 10];    % initial and final time to simulate
x0 = 1e-6*[0 0 0 0 0 0 0 0 0 0 0 0]';       % initial point
xdot0 = 0*1e-6*[1 1 1 1 1 1 1 1 1 1 1 1]';

% options = odeset('abstol', 1e-6, 'reltol', 1e-6, 'NonNegative', 1);

% % Check that initial conditions make sense
% fixed_x0 = [1 1 1 0 0 0];
% fixed_xdot0 = [0 0 0 0 0 0];
% [x0_new,xdot0_new] = decic(@(t, x, xdot)vf(x,setInput(t, params),xdot,params),0,x0,fixed_x0,xdot0,fixed_xdot0);

% Simulate system response
[t, y] = ode15i(@(t, x, xdot)vf(x, setInput(t, params), xdot, params), tspan, x0, xdot0);

%% Create vector of inputs u(t), for plotting
u = zeros(length(t), sum(params.n));
for i = 1:length(t)
    u(i,:) = setInput(t(i), params);
end

%% Plot the results
figure
plot(t,y(:,1:3))
legend('psi','theta','phi')

figure
plot(t,y(:,4:6))
legend('psidot','thetadot','phidot')

figure
plot(t, u*10^(-3))
title('Control Inputs')
legend('FREE-1', 'FREE-2', 'FREE-3', 'FREE-4', 'location', 'southeast')
xlabel('time (s)')
ylabel('Pressure (kPa)')