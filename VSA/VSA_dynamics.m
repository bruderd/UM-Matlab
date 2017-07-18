function xdot = VSA_dynamics(t, x, u, tau, params)
%Differential equation describing the dynamics of the VSA
%   Detailed explanation goes here

[B, N, L, I, rho] = deal(params.B, params.N, params.L, params.I, params.rho);
[u1, u2] = deal(u(1), u(2));
[x1, x2] = deal(x(1), x(2));

damp = -1;

xdot(1,1) = x2;
xdot(2,1) = (rho/(4*pi*I*N^2)) * ( (B^2 - 3*(L+rho*x1)^2)*u1 + (-B^2 + 3*(L-rho*x1)^2)*u2 ) + tau/I + damp*x2;

end

