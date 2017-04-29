% elastomerModulus.m 

function [Emod, Gmod] = elastomerModulus(x, params)

[P, gama, r, L, phi] = deal(x(1), x(2), x(3), x(4), x(5));
[P0, gama0, r0, L0, phi0, T0] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4), params.x_rest(5), params.x_rest(6));
t0 = params.t_rest;
[Fload, Mload] = deal(params.load(1), params.load(2));

syms E G T

%% Stress/strain equations
sig_z = E * (L - L0)/L0;
sig_theta = E * (r - r0)/r0;        % removed factor of 2*pi
tau_ztheta = G * atan(r*phi/L);     % added atan since x ~= tan(x) only for small x

%% Tube wall thickness equation
% t = -r + sqrt(r^2 + 2*r0*t0 +t0^2); 
t = t0;
%% Force Balance Equations
hoop_balance = 2*pi*P*r^2*cot(gama) - 2*T*sin(gama) - 2*sig_theta*(pi*r*cot(gama))*t;
force_balance = P*pi*r^2 - 2*T*cos(gama) - pi*(2*r*t + t^2)*sig_z + Fload;
torque_balance = 2*(r+t)*T*sin(gama) - pi*(2*r*t + t^2)*(r + t/2)*tau_ztheta + Mload;

modulus = solve([hoop_balance; force_balance; torque_balance], [E, G, T]);

Emod = modulus.E;
Gmod = modulus.G;

end