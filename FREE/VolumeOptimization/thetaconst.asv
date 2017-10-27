function [ err ] = thetaconst( x, params )
%thetaconst: Constraint due to the inextensibility of the fiber
%   Uses a discrete sum to approximate the total fiber windings in radians (theta).
%   Then compares it to what theta should be (Theta + phi)

[l, phi, a0, a1, a2, a3, a4, a5, a6] = deal(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9));
[Gama, L, R, Theta] = deal(params.Gama, params.L, params.R, params.Theta);
n = params.slices;

p = [a6, a5, a4, a3, a2, a1, a0];
sk = L/cos(Gama) *  (1/n);
tk = l/n;

num = sqrt(sk^2 - tk^2);
xk = linspace(0,l,n);
denom = polyval(p, xk);

thetak = num ./ denom;

theta = sum(thetak);

err = theta - (Theta + phi);

% stop for debugging if returns complex value
if ~isreal(err)
    ;
end

end

