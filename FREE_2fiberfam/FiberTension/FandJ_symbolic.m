% Defines system of force balance equations and its jacobian
clear

%% Define symbolic variables
syms r0 L0 gama0 dP dgama dr dL c1 c2 c3 c4 c5 c6 betta0 nrat T_gama T_betta

assume(r0, 'real')
assumeAlso(L0, 'real')
assumeAlso(gama0, 'real')
assumeAlso(betta0, 'real')
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
assumeAlso(nrat, 'real')

% syms t 
syms P gama betta r L phi T_gama T_betta

% P = sym('P(t)');
% gama = sym('gama(t)');
% betta = sym('betta(t)');
% r = sym('r(t)');
% L = sym('L(t)');
% phi = sym('phi(t)');
% T_gama = sym('T_gama(t)');
% T_betta = sym('T_betta(t)');

assumeAlso(P, 'real')
assumeAlso(gama, 'real')
assumeAlso(betta, 'real')
assumeAlso(r, 'real')
assumeAlso(L, 'real')
assumeAlso(phi, 'real')
assumeAlso(T_gama, 'real')
assumeAlso(T_betta, 'real')

%%

% Define state vector
x = [P, gama, betta, r, L, phi, T_gama, T_betta];

% Difine some intermediary parameters
theta_gama0 = -tan(gama0)*L0/r0;     
theta_betta0 = -tan(betta0)*L0/r0;   
theta_gama = -tan(gama)*L/r; 
theta_betta = -tan(betta)*L/r; 

F_elast = c1*(-1)*(L-L0);
M_elast = c4*(-1)*phi;
% F_elast = 0;
% M_elast = 0;

force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)           
geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
geometry_constraint3 = (theta_gama - theta_gama0) - phi; 
geometry_constraint4 = (theta_betta - theta_betta0) - phi;
extra_constraint1 = P*r - (T_gama*sin(abs(gama)) + T_betta*sin(abs(betta)));
% extra_constraint1 = P*r - (T_gama*sin(gama) + T_betta*sin(betta));  % no absolute value
extra_constraint2 = P*r - (T_gama*sin(abs(-gama)) + T_betta*sin(abs(-betta)));
% extra_constraint2 = T_gama - T_betta;       % shot in the dark, what if the tensions are equal?-DOESN'T SEEM TO WORK

% Force Balance system
F = [force_balance;...
      torque_balance;...
      geometry_constraint1;...
      geometry_constraint2;...
      geometry_constraint3;...
      geometry_constraint4;...
      extra_constraint1;...
      extra_constraint2];
  
% Jacobian, or gradient of force balance system  
J = jacobian(F,x);