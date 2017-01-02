clear all;

% syms r0 L0 gama0 P gama r L dP dgama dr dL c1 c2 c3 c4 c5 c6
syms r0 L0 gama0 dP dgama dr dL c1 c2 c3 c4 c5 c6

assume(r0, 'real')
assumeAlso(L0, 'real')
assumeAlso(gama0, 'real')
assumeAlso(dP, 'real')
assumeAlso(dgama, 'real')
assumeAlso(dr, 'real')
assumeAlso(dL, 'real')
assumeAlso(c1, 'real')
assumeAlso(c2, 'real')
assumeAlso(c3, 'real')
assumeAlso(c4, 'real')
assumeAlso(c5, 'real')
assumeAlso(c6, 'real')


syms t

P = sym('P(t)');
gama = sym('gama(t)');
r = sym('r(t)');
L = sym('L(t)');
phi = (-tan(gama)*L/r + tan(gama0)*L0/r0);

assumeAlso(P, 'real')
assumeAlso(gama, 'real')
assumeAlso(r, 'real')
assumeAlso(L, 'real')

% Definition of elastomer spring force functions, constants from sys id
% experimental data
F_elast = [c1 c2 c3] * [L^2, L, 1]';
M_elast = [c4 c5 c6] * [phi^2, phi, 1]';


% Differentiates the system of equations wrt time
force_balance = diff(0 == -P*pi*r^2 + 2*P*pi*r^2*cot(gama)^2 + F_elast,t);
torque_balance = diff(0 == -2*P*pi*r^3*cot(gama) + M_elast,t);
geometry_constraint = diff(0 == -cos(gama) + (L/L0)*cos(gama0),t);

system = [force_balance; torque_balance; geometry_constraint];

% Here we take 'system' into a text editor to make the following
% substitutions:
%	P(t) ?> P
%	r(t) ?> r
%	gama(t) ?> gama
%	L(t) ?> L
%	diff(P(t), t) ?> dP
%	diff(r(t), t) ?> dr
%	diff(gama(t), t) ?> dgama
%	diff(L(t), t) ?> dL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Take 'system' and find the derivatves with respect to:
%       x = [P, gama, r, L, phi]'
%       u = u
%       xdot = [dP, dgama, dr, dL, dphi]'


syms r0 L0 gama0 dP dgama dr dL c1 c2 c3 c4 c5 c6 phi phi0 x u dxdt dphi P gama r L

assume(r0, 'real')
assumeAlso(L0, 'real')
assumeAlso(gama0, 'real')
assumeAlso(phi0, 'real')
assumeAlso(dP, 'real')
assumeAlso(dgama, 'real')
assumeAlso(dr, 'real')
assumeAlso(dL, 'real')
assumeAlso(dphi, 'real')
assumeAlso(c1, 'real')
assumeAlso(c2, 'real')
assumeAlso(c3, 'real')
assumeAlso(c4, 'real')
assumeAlso(c5, 'real')
assumeAlso(c6, 'real')
assumeAlso(P, 'real')
assumeAlso(gama, 'real')
assumeAlso(r, 'real')
assumeAlso(L, 'real')
assumeAlso(phi, 'real')


x = [P, gama, r, L, phi]';
xdot = [dP, dgama, dr, dL, dphi]';

f = [-dP + (u - P);...
         c2*dL + 2*c1*L*dL - pi*r^2*dP - 2*pi*P*r*dr + 2*pi*cot(gama)^2*r^2*dP + 4*pi*cot(gama)^2*P*r*dr - 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1);...
         2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP - c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 6*pi*cot(gama)*P*r^2*dr;...
         sin(gama)*dgama + (cos(gama0)*dL)/L0;...
         -phi + (-tan(gama)*L/r + tan(gama0)*L0/r0)];

dfdx = jacobian(f, x);

dfdu = jacobian(f, u);

dfdxdot = jacobian(f, xdot);



