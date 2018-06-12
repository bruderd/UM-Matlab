%% ODE15i Simulation for MWRW2018

%% ODE15i Simulation
%   Use this to check if the dynamics are broken, or if it's something else

% tspan = [0, params.tfinal];
tspan = [0:1e-3:params.tfinal];
x0 = 1e-6*[0 0 0 0 0 0]';       % initial point
xdot0 = 0*1e-6*[1 1 1 1 1 1]';
% u = [1000, 0, 0, 1000]';

% options = odeset('abstol', 1e-6, 'reltol', 1e-6, 'NonNegative', 1);

% Check that initial conditions make sense
fixed_x0 = [1 1 1 0 0 0];
fixed_xdot0 = [0 0 0 0 0 0];
[x0_new,xdot0_new] = decic(@(t, x, xdot)vf(x,calc_u(t, params),xdot,params),0,x0,fixed_x0,xdot0,fixed_xdot0);

% global errP; % global error variable declaration
global errI; % global integral error variable declaration (needed for the integral controller)
global time; % global variable that keeps track of time
% errP = zeros(6,1);
errI = zeros(6,1);
time = 0;

% Simulate system response
% options = odeset('OutputFcn',@(t,x) calc_u_MWRW(x, t, params));
[t, y] = ode15i(@(t, x, xdot) vf(x, calc_u_MWRW(x, t, params), xdot, params), tspan, x0_new, xdot0_new);

%% Create vector of inputs u(t)
errI = zeros(6,1);  % reinitialize global variables
time = 0;
u = zeros(length(t), params.m);
error = zeros(6, length(t));
yeul_des = zeros(3, length(t));
for i = 1:length(t)
    [u(i,:), error(:,i), yeul_des(:,i)] = calc_u_MWRW(y(i,:)', t(i), params);
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

figure
plot(t,error(:,:))
title('Error')
legend('x', 'y', 'z', 'psi','theta','phi')



