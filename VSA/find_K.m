function K = find_K(theta, u, params)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here

[B, N, L, I, rho] = deal(params.B, params.N, params.L, params.I, params.rho);

K = (-3*rho^2)/(2*pi*N^2) * ( (L + rho*theta)*u(1) +  (L - rho*theta)*u(2) );

end

