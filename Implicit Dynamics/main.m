clear
clc

params = struct;

params.x0 = [1]';    %state vector initial condition
params.u0 = [5]';        %input vector initial condition

params.T = 20;   %final time
params.N = 100; %number of steps
params.dt = params.T/params.N;    %size of one time step
params.n = length(params.x0);  %dimension of state vector x
params.m = length(params.u0);  %dimension of state vector u

params.vf = @( x, u ) vf( x, u );

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
s_star = fmincon(@(s)obj(s, params), s0, [], [],[],[], lb, ub, @(s)implicit_Dynamics(s,params), options);
% s_star = fmincon(@(s)obj(s, params), s0, [], [],[],[], lb, ub, [], options); %no constraints

% convert results of optimization to states and inputs
[x_star, u_star] = state_decode(s_star, params);
