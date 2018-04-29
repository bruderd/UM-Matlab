function results = gurobiPressure( x, params, fload)
%calcPressure: Calculates the input pressure needed to acheive some
%end effector position x by solving for the equilibrium point.
%   FYI the last element of psol is the constraint tolerance at that point

% check for optional argument fload
if ~exist('fload','var')
     % fload input does not exist, so default it to 0
      fload = zeros(6,1);
end

num = params.num;


%% calculate quadratic cost function matrices
[H, f, A, b, Aeq, beq] = gurobiConst(x, params, fload, params.penalty);
% A = zeros(num, num);
% b = zeros(num, 1);
% Aeq = zeros(num, num);
% beq = zeros(num, 1);

%% solve QP using gurobi
names = {'p1', 'p2', 'p3', 'tol1', 'tol2', 'tol3', 'tol4', 'tol5', 'tol6'};
model.varnames = names;
model.lb = [params.pmin, zeros(1,6)];
model.ub = [params.pmax, params.tol];

model.Q = sparse(H);
model.A = sparse(A);
model.obj = f;
model.rhs = b;
model.sense = '<';

results = gurobi(model);

% %% solve for the pressure at equilibrium point
% options = optimoptions('quadprog', 'Display', 'iter');  % 'ConstraintTolerance', 2.5e-2);
% tolmin = zeros(1,6); % the minimum equality constraint tolerance
% tolmax = params.tol;  % the maximum equality constraint tolerance






end






