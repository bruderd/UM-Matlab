function [psol, exitflag] = calcPressure( x, params, fload)
%calcPressure: Calculates the input pressure needed to acheive some
%end effector position x by solving for the equilibrium point.
%   Detailed explanation goes here

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
[H, f, Aeq, beq] = quadCost(x, params, fload, params.penalty);
A = zeros(num, num);
b = zeros(num, 1);
% Aeq = zeros(num, num);
% beq = zeros(num, 1);

%% solve for the pressure at equilibrium point
[psol, ~, exitflag] = quadprog(H,f,[],[],Aeq,beq,params.pmin,params.pmax);


%% solution using lsqlin for comparision. It has no penalty term so solutions differ from qp
% plin = lsqlin(Aeq,beq,[],[],Aeq,beq,params.pmin,params.pmax);


%% solve using lsqnonlin (cannot handle casese with non-unique solutions)
% fun = @(p) calcf(x, p, params) + fload;
% options = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective');
% plsq = lsqnonlin(fun, p0, params.pmin', params.pmax', options);


end

