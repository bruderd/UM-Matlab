% solve_FBhoop just solves for deformation of 1fiber FREE at given pressure

function [P, gama, r, L, phi, T] = solve_FBhoop(P_test, x_rest, load)

[P0, gama0, r0, L0, phi0, T0] = deal(x_rest(1), x_rest(2), x_rest(3), x_rest(4), x_rest(5), x_rest(6));

dx = 1e-3;

% Don't fix anything but the input pressure
LB = [P_test-dx, -pi/2, r0, L0*0.25, -Inf, 0]; 
UB = [P_test+dx, pi/2, r0*4, L0*2, Inf, Inf];

% solver variable
x = lsqnonlin( @(x) FBhoop(x, P_test, x_rest, x_rest, load), x_rest, LB, UB );

P = x(1);
gama = x(2);
r = x(3);
L = x(4);
phi = x(5);
T = x(6);

end