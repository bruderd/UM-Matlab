function u = controlLaw_p(x, params)
%Proportional control law
%   Detailed explanation goes here

[B, N, L, I, rho] = deal(params.B, params.N, params.L, params.I, params.rho);
K = params.K_des;
[theta, thetadot] = deal(x(1), x(2));

error = -x + params.x_des;

M = params.C * error;

u(1,1) = pi*N^2*(B^2*K - 3*(L-rho*theta)*(2*M*rho + K*(L-rho*theta))) / (3*L*rho^2 * (-B^2 + 3*L^2 - 3*rho^2*theta^2));
u(1,2) = pi*N^2*(B^2*K - 3*(L+rho*theta)*(-2*M*rho + K*(L+rho*theta))) / (3*L*rho^2 * (-B^2 + 3*L^2 - 3*rho^2*theta^2));


end

