% simulation.m

clear
clc

%% Set system parameter values

params = struct;

p = 6;
n = [2]';
params.p = p;      % total number of spine segments
params.n = n;     % numer of actuators in each module (a vector)

% spine parameters
params.L = [0.01]';%[0.3048]';    % length of each module (m)
params.dl = params.L / params.p;
params.K = -0.001 * eye(p);
params.D = -0.1 * eye(p);
params.M = 0.005 * eye(3*p);   % mass of each vertebrae (kg)
params.M(3*p-2 : 3*p, 3*p-2 : 3*p) = eye(3) * 0.030;  % mass of end block (kg)
params.I = 0.0001 * eye(p);   % inertia of each vertebrae
params.I(p,p) = 0.0030;  % inertia of end block

% actuator paramters: free = [Gama, R, xattach, yattach], size(free) = sum(n) x 4
params.Lfree = ones(1,n) * params.L;
params.Gama = deg2rad([40, 40]');  % fiber angle (rad)
params.R = 1 * ones(sum(n), 1) * 1e-2;     % internal radius (m)
params.attach = [2, -2]' * 1e-2;   % x-coordinates of attachment points
params.B = abs(params.Lfree ./ cos(params.Gama));   % fiber length (must be positive))
params.Nf = -params.Lfree ./ (2*pi*params.R) .* tan(params.Gama); % total fiber windings in revolutions (when relaxed)

disp('Initialized parameters.')

%% Create various functions

setJacobians(params);   % derive the appropriate Jacobian matrices
disp('Symbolically derived Jacobian matrices.')

setEOM(params);     % derive the equations of motion
disp('Symbolically derived Euler-Lagrange equations of motion.')

%% ODE15i Simulation
disp('Simulating dynamics...')

tspan = [0, 10];    % initial and final time to simulate
% x0 = 1e-6*[0 0 0 0 0 0 0 0 0 0 0 0]';       % initial point
% xdot0 = 0*1e-6*[1 1 1 1 1 1 1 1 1 1 1 1]';
x0 = zeros(2*p,1);
xdot0 = zeros(2*p,1);

% options = odeset('abstol', 1e-6, 'reltol', 1e-6, 'NonNegative', 1);

% Check that initial conditions make sense
fixed_x0 = ones(size(x0'));
fixed_xdot0 = zeros(size(xdot0'));
[x0_new,xdot0_new] = decic(@(t, x, xdot)vf(x,setInput(t, params),xdot,params),0,x0,fixed_x0,xdot0,fixed_xdot0);

% Simulate system response
[t, y] = ode15i(@(t, x, xdot)vf(x, setInput(t, params), xdot, params), tspan, x0_new, xdot0_new);

%% Create vector of inputs u(t), for plotting
u = zeros(length(t), sum(params.n));
for i = 1:length(t)
    u(i,:) = setInput(t(i), params);
end

%% Plot final position in xy-coordinates
X_final = alpha2x(y(end, :), params);
x_final = X_final(1:3:end);
y_final = X_final(2:3:end);

figure
plot(x_final, y_final)

%% Convert x0_orient to x0
% Y = zeros(length(y), 6*p);
% for i = 1 : length(t)
%     Y(i, 1:6*p) = x0_orient2x0(y(i, 1:3*p)', params)';
% end

%% Plot the results
% figure
% plot(t,y(:,1:3))
% legend('psi','theta','phi')
% 
% figure
% plot(t,y(:,4:6))
% legend('psidot','thetadot','phidot')

% figure
% plot(t, u*10^(-3))
% title('Control Inputs')
% legend('FREE-1', 'FREE-2', 'FREE-3', 'FREE-4', 'location', 'southeast')
% xlabel('time (s)')
% ylabel('Pressure (kPa)')