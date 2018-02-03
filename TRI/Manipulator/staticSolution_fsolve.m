% simulation.m

clear
clc

%% Set system parameter values

p = 3;      % total number of modules
n = [4, 4, 4]';     % numer of actuators in each module (a vector)

% module parameters: module = [L, block density, EI of spine], size(module) = p x 3
L = [0.1, 0.1, 0.1]';    % length of each module (m)
density = [1e3, 1e3, 1e3]';      % density of each module block (kg/m^3)
dim = [5, 5, 5; 5, 5, 5; 5, 5, 5] * 1e-2;     % length (along x-axis), width (along y-axis), height (along z-axis) of end blocks, assuming they are rectangular prisms (m)
EIspine = [0 0 0]';       % effective stiffness/inertia of spine material (currently is not factored in since spine is considered forceless)
module = [L, density, dim, EIspine];

% actuator paramters: free = [Gama, R, xattach, yattach], size(free) = sum(n) x 4
Gama = deg2rad([40, -40, 40, -40, 40, -40, 40, -40, 40, -40, 40, -40]');  % fiber angle (rad)
R = 1 * ones(sum(n), 1) * 1e-2;     % internal radius (m)
xattach = [2, -2, -2, 2, 2, -2, -2, 2, 2, -2, -2, 2]' * 1e-2;   % x-coordinates of attachment points
yattach = [2, 2, -2, -2, 2, 2, -2, -2, 2, 2, -2, -2]' * 1e-2;   % y-coordinates of attachment points
free = [Gama, R, xattach, yattach];

% simulation parameters (might want to put some simulation parameters here at some point...)
% tspan = [0, 10];
% x0_0 = zeros(6*p, 1);
% x0dot_0 = zeros(6*p, 1);
sim = 0;

% create params struct to be passed to many functions
params = setParams(p, n , module, free, sim);
disp('Initialized parameters.')

%% Create various functions

setJacobians_v2(params);   % derive the appropriate Jacobian matrices
disp('Symbolically derived Jacobian matrices.')

%% Solve for state of the system under control inputs defined

% set the input pressure at each time step
% t = linspace(0,10,100);
% for i = 1:length(t)
%    u(:,i) = setInput(t(i), params);
%    x_star(:,i) = fsolve(@(x) staticBalance(x, u(:,i), params), zeros(3*p,1));
% end

% State given Input, not a bunch in "time"
t = 0;  % placeholder, does nothing
u = setInput(t, params);
x_star = fsolve(@(x) staticBalance(x, u, params), zeros(3*p,1));

% % Input Given State (should return the same thing as setInput-- so far this does not work...)
% u_star = fsolve(@(u) staticBalance(x_star, u, params), zeros(sum(n),1));

%% Solve for u instead of x

% x_star = [0 0 0.5 0 0 0 0 0 0]';
% u_star = fsolve(@(u) staticBalance(x_star, u, params), zeros(4*p,1));
% u_star = u_star.^2;


%% Convert from local orientation to global orientation and position

x0_star = x_orient2x0_orient(x_star, params);
X0_star = x0_orient2x0(x0_star, params);


% %% Create vector of inputs u(t), for plotting
% u = zeros(length(t), sum(params.n));
% for i = 1:length(t)
%     u(i,:) = setInput(t(i), params);
% end
% 
% %% Convert x0_orient to x0
% Y = zeros(length(y), 6*p);
% for i = 1 : length(t)
%     Y(i, 1:6*p) = x0_orient2x0(y(i, 1:3*p)', params)';
% end
% 
% %% Plot the results
% figure
% plot(t,y(:,1:3))
% legend('psi','theta','phi')
% 
% figure
% plot(t,y(:,4:6))
% legend('psidot','thetadot','phidot')
% 
% % figure
% % plot(t, u*10^(-3))
% % title('Control Inputs')
% % legend('FREE-1', 'FREE-2', 'FREE-3', 'FREE-4', 'location', 'southeast')
% % xlabel('time (s)')
% % ylabel('Pressure (kPa)')