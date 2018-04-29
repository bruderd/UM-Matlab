function [ H, f , A, b, Aeq, beq] = gurobiConst( x, params, fload, penalty )
%quadCost: Calculates the H and f matrices for the quadratic cost function
%whose minimum lies has at the equilibrium point.
%   See Matlab documentation on 'quadprog' for more info.

% check for optional arguments
if ~exist('fload','var')
     % fload input does not exist, so default it to 0
      fload = zeros(6,1);
end
if ~exist('penalty','var')
     % penalty input does not exist, so default it to 1e-5
      penalty = 1e-5;
end

num = params.num;
[C, D] = deal(params.C, params.D);

q = x2q(x);
Jq = calcJq(q, params);

% calculate the elastomer force
felast = calcFelast( x, params );

%% Put the tolerance value into the cost and optimize over it along with p, don't even minimize pressure

H = blkdiag( zeros(num,num), eye(6) );
f = zeros(num + 6,1);

% need more slack so will use inequality constraints with tolerance
A = [D*Jq', -eye(6);...
    -D*Jq', -eye(6)];
b = [-(felast + fload);...
    (felast + fload)];

Aeq = [];
beq = [];