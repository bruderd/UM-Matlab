%% LQR Test

clear
clc

params = struct;

params.m = 1; %kg
params.k = 1; %N/m
params.b = 0.5; %Ns/m (damping constant)
params.F = 0; %N (input force) %unforced for now

params.A = [0 1; -params.k/params.m -params.b/params.m];
params.B = [0 1/params.m]';
params.C = [1 0];
params.D = [0];

params.x0 = [0.2 0]';    %state vector initial condition
params.u0 = [0]';        %input vector initial condition

params.T = 20;   %final time
params.N = 200; %number of steps
params.dt = params.T/params.N;    %size of one time step
params.n = length(params.x0);  %dimension of state vector x
params.m = length(params.u0);  %dimension of state vector u

s0 = state_encode0(params);

% Dummy inequality constraints. Defined as A*s <= B
A = eye(length(s0));
B = ones(length(s0),1)*100;

% % defining options for fmincon
options = optimoptions('fmincon','Algorithm','sqp',...
   'GradObj', 'on', 'GradConstr', 'on', 'Display', 'iter',  ...
   'MaxIter', 1000, 'ScaleProblem', 'obj-and-constr', 'Hessian', 'lbfgs'  );

% Lower and upper bounds of states/inputs
lb = -10*ones(length(s0),1);
ub = 10*ones(length(s0),1);

% NL program to solve optimization problem
s_star = fmincon(@(s)obj(s, params), s0, [], [],[],[], lb, ub, @(s)dynamics_LQR_v2(s,params), options);
% s_star = fmincon(@(s)obj(s, params), s0, [], [],[],[], lb, ub, [], options); %no constraints

% convert results of optimization to states and inputs
[x_star, u_star] = state_decode(s_star, params);
