function [ H, f ] = quadCost( x, params, fload, penalty )
%quadCost: Calculates the H and f matrices for the quadratic cost function
%whose minimum lies has at the equilibrium point.
%   See Matlab documentation on 'quadprog' for more info.

% check for optional arguments
if ~exist('fload','var')
     % fload input does not exist, so default it to 0
      fload = zeros(6,1);
end
if ~exist('penalty','var')
     % penalty input does not exist, so default it to 1
      penalty = 1;
end

num = params.num;
[C, D] = deal(params.C, params.D);

q = x2q(x);
Jq = calcJq(q, params);

% define output matrices
H = 2 * Jq*(D'*D)*Jq' + penalty*eye(num);
f = 2 * Jq*D'*(C*q + fload);

end

