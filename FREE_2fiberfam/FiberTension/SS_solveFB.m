% SS_solve_FB - solve for the full state given a steady state pressure, u

function [Tgama, Tbetta, P, gama, betta, r, L, phi] = SS_solveFB(u, x_rest)


[P0, gama0, betta0, r0, L0, phi0] = deal(x_rest(1), x_rest(2), x_rest(3), x_rest(4), x_rest(5), x_rest(6));

x_test = [u, gama0, betta0, r0, L0, phi0];

dx = 1e-3;

% LB = [0, 0, u-dx, gama0-dx, betta0-dx, r0-dx, L0-dx, phi0-dx]; 
% UB = [inf, inf, u+dx, gama0+dx, betta0+dx, r0-dx, L0+dx, phi0+dx];
LB = [0, 0, u-dx, -pi/2, -pi/2, r0*1.001, L0*0.6, -Inf]; 
UB = [inf, inf, u+dx, pi/2, pi/2, r0*3, Inf, Inf];

% solver variable
Tx = lsqnonlin(@(Tx) FB(u, Tx(1:2), Tx(3:8), x_rest), [0, 0, x_rest], LB, UB);

Tgama = Tx(1);
Tbetta = Tx(2);
P = Tx(3);
gama = Tx(4);
betta = Tx(5);
r = Tx(6);
L = Tx(7);
phi = Tx(8);

end