% main.m
%   Runs optimization to find optimal solution
%   Note: Must run setParams.m before this script

ready = exist('params', 'var');

% Check that param struct has been created
if ready
    
    s0 = state_encode0(params);
    
    % Dummy inequality constraints. Defined as A*s <= B
    A = eye(length(s0));
    B = ones(length(s0),1)*100;
    
    % % defining options for fmincon
    options = optimoptions('fmincon','Algorithm','interior-point',...
        'GradObj', 'on', 'GradConstr', 'on', 'Display', 'iter',  ...
        'MaxIter', 1000, 'ScaleProblem', 'obj-and-constr', 'Hessian', 'lbfgs'  );
    
    % Lower and upper bounds of states/inputs
    lb = -10000*ones(length(s0),1);
    ub = 10000*ones(length(s0),1);
    
    % NL program to solve optimization problem
    s_star = fmincon(@(s)obj(s, params), s0, [], [],[],[], lb, ub, @(s)constraints(s,params), options);
    
    % convert results of optimization to states and inputs
    [x_star, u_star] = state_decode(s_star, params);

else
    'You must set initial conditions and solvler parameters by running setParams.m first bro'
end
    
