function params = setParams( Gama, L, R)
%setParams: Defines resting FREE parameters
%   Also defines cost function used for optimization that finds the maximal volume shape of a FREE.

params.Gama = Gama;
params.L = L;
params.R = R;
params.Theta = L*tan(Gama)/R;
params.S = L/cos(params.Gama);
params.slices = 1e6;    % number of slices in fiber constraint

syms l phi x a0 a1 a2 a3 a4 a5 a6

r = poly2sym([a6 a5 a4 a3 a2 a1 a0],x);

cost = -int(r^2,x,0,l);
gradcost = jacobian(cost,[l phi a0 a1 a2 a3 a4 a5 a6]);

drdx = diff(r,x);

% fiber inextensibility constraint (unused)
fibconst = int(sqrt(drdx^2 + (((params.Theta + phi)/l)^2 - drdx^2) * r^2 + 1),x,0,l) - L/cos(Gama);
% fibconst = int(sqrt(drdx^2 + 1),x,0,l) - L; % test with straight fiber (0 deg)


% create matlab functions for the cost and fiber constraints
matlabFunction(cost, 'File', 'cost', 'Vars', {[l, phi, a0, a1, a2, a3, a4, a5, a6]});
matlabFunction(gradcost, 'File', 'gradcost', 'Vars', {[l, phi, a0, a1, a2, a3, a4, a5, a6]});
matlabFunction(fibconst, 'File', 'fibconst', 'Vars', {[l, phi, a0, a1, a2, a3, a4, a5, a6]}, 'Optimize',false);

end

