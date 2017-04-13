% function [f] = (params

%% Check dynamics constraints at a single point

% Point at which we want to evaluate the feasability
x = [10.02, 0.67163, 0.21484, 5.8316]';
xdot = [0 0 0 0]';  % asssume steady state

n = params.n;
m = params.m;
    
[c1, c2, c3] = deal(params.Felast_consts(1), params.Felast_consts(2), params.Felast_consts(3));
[c4, c5, c6] = deal(params.Melast_consts(1), params.Melast_consts(2), params.Melast_consts(3));
    
[P0, gama0, r0, L0] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4));
[P, gama, r, L] = deal(x(1), x(2), x(3), x(4));
u = P;  % we assume it is at steady state
[dP, dgama, dr, dL] = deal(xdot(1), xdot(2), xdot(3), xdot(4));

f = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary
     c2*dL + 2*c1*L*dL - pi*r^2*dP - 2*pi*P*r*dr + 2*pi*cot(gama)^2*r^2*dP + 4*pi*cot(gama)^2*P*r*dr - 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1);...
     c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr;...
     sin(gama)*dgama + (cos(gama0)*dL)/L0];
 
f_v2 = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary
         c2*dL + 2*c1*L*dL + pi*r^2*dP + 2*pi*P*r*dr - 2*pi*cot(gama)^2*r^2*dP - 4*pi*cot(gama)^2*P*r*dr + 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1);...
         c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr;...
         sin(gama)*dgama + (cos(gama0)*dL)/L0];