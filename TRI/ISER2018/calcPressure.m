function [psol, plsq] = calcPressure( x, params, p0, fload)
%calcPressure: Calculates the input pressure needed to acheive some
%end effector position x by solving for the equilibrium point.
%   Detailed explanation goes here

% check for optional argument fload
if ~exist('fload','var')
     % fload input does not exist, so default it to 0
      fload = zeros(6,1);
end

num = params.num;

% calculate quadratic cost function matrices
[H,f] = quadCost(x, params, fload, params.penalty);
A = zeros(num, num);
b = zeros(num, 1);
Aeq = zeros(num, num);
beq = zeros(num, 1);

% solve for the pressure at equilibrium point
psol = quadprog(H,f,A,b,Aeq,beq,params.pmin,params.pmax);


%% solve using lsqnonlin (cannot handle casese with non-unique solutions)
fun = @(p) calcf(x, p, params) + fload;
options = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective');
plsq = lsqnonlin(fun, p0, params.pmin', params.pmax', options);


end

