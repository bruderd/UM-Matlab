% solve_FBstress just solves for deformation of 1fiber FREE at given pressure

function [P, gama, r, L, phi, T] = solve_FBstress(u, params)

[P_rest, gama_rest, r_rest, L_rest, phi_rest, T_rest] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4), params.x_rest(5), params.x_rest(6));

% Don't fix anything but the input pressure
dx = 1e-3;
LB = [u-dx, -pi/2, r_rest, L_rest*0.25, -Inf, 0]; 
UB = [u+dx, pi/2, r_rest*4, L_rest*2, Inf, Inf];

% solver variable
x = lsqnonlin( @(x) FBstress(x, u, params), params.x_rest, LB, UB );

P = x(1);
gama = x(2);
r = x(3);
L = x(4);
phi = x(5);
T = x(6);

end