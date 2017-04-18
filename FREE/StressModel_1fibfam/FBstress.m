% FBhoop.m
%   Force balance equations that account for the effect of the elastomer
%   bearing hoop stress
%
%   x = [P, gama, r, L, phi, T]
%   u = P_desired
%   params includes:
%       params.x_rest = "states of FREE when P=0 and no loads applied"
%       params.load = [Fload, Mload]
%       params.t_rest = "thickness of tube when P=0 and no loads applied"
%
%   NOTE: This is for a FREE with 1 fiber family

function [F, J] = FBstress(x, u, params)

[P, gama, r, L, phi, T] = deal(x(1), x(2), x(3), x(4), x(5), x(6));
[P0, gama0, r0, L0, phi0, T0] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4), params.x_rest(5), params.x_rest(6));

%% Definition of modulus of elasticity equations
dL_norm = (L - L0)/L0;
dphi_norm = atan(r*phi)/L0;

% E = params.modulus(1);    % constant modulus
% G = params.modulus(2);    % constant modulus
E = params.modulus(1,1)*dL_norm + params.modulus(1,2);  % varying modulus
G = params.modulus(2,1)*dphi_norm + params.modulus(2,2);    % varying modulus

%% Evaluate FB equations and gradient at current point
F = Feval(x, params.x_rest, params.t_rest, E, G, params.load, u);
F = F*1e4;

if nargout > 1  % Two output arguments
    J = Jeval(x, params.x_rest, params.t_rest, E, G, params.load, u);     % Jacobian of F evaluated at x
    J = J*1e4;
end

end