function [psol, exitflag] = calcPressure( x, params, fload)
%calcPressure: Calculates the input pressure needed to acheive some
%end effector position x by solving for the equilibrium point.
%   FYI the last element of psol is the constraint tolerance at that point

% check for optional argument fload
if ~exist('fload','var')
     % fload input does not exist, so default it to 0
      fload = zeros(6,1);
end

num = params.num;

%% Check that the desired end effector position is even physically possible given physical constraints of system
% If it is possible, replace the desired x with the desired x written in
% terms of the independent components.

% [xgoal, feasability] = checkfeas_xdes( x, params );
% if ~feasability
%     psol = [];
%     exitflag = -2;
%     return;
% else
%     x = xgoal';
% end


%% calculate quadratic cost function matrices
[H, f, A, b, Aeq, beq] = quadCost(x, params, fload, params.penalty);
% A = zeros(num, num);
% b = zeros(num, 1);
% Aeq = zeros(num, num);
% beq = zeros(num, 1);

%% solve for the pressure at equilibrium point
options = optimoptions('quadprog', 'Display', 'iter');  % 'ConstraintTolerance', 2.5e-2);
tolmin = zeros(1,6); % the minimum equality constraint tolerance
tolmax = params.tol;  % the maximum equality constraint tolerance

% pmax = (6894.76) * [12 12 22]; % unused

[psol, ~, exitflag] = quadprog(H, f, A, b, [], [], [params.pmin, tolmin], [params.pmax, tolmax], [], options);


%% older version where tolerance is not optimized over
% [psol, ~, exitflag] = quadprog(H, f, [], [], Aeq, beq, params.pmin, params.pmax, [], options);


%% solution using lsqlin for comparision. It has no penalty term so solutions differ from qp
% plin = lsqlin(Aeq,beq,[],[],Aeq,beq,params.pmin,params.pmax);


%% solve using lsqnonlin (cannot handle casese with non-unique solutions)
% fun = @(p) calcf(x, p, params) + fload;
% options = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective');
% plsq = lsqnonlin(fun, p0, params.pmin', params.pmax', options);


end

