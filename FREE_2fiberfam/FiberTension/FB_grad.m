% FB_grad.m
%   Force balance equations, with the gradient as well

function [F,J] = FB_grad(P_test, T, x, x_rest)

P = P_test;
[T_gama, T_betta] = deal(T(1), T(2));
[P, gama, betta, r, L, phi] = deal(x(1), x(2), x(3), x(4), x(5), x(6));
[P0, gama0, betta0, r0, L0, phi0] = deal(x_rest(1), x_rest(2), x_rest(3), x_rest(4), x_rest(5), x_rest(6));

theta_gama0 = -tan(gama0)*L0/r0;     
theta_betta0 = -tan(betta0)*L0/r0;   
theta_gama = -tan(gama)*L/r; 
theta_betta = -tan(betta)*L/r; 

% Elastomer forces
c1 = 7;
c4 = 7;
F_elast = c1*(-1)*(L-L0);
M_elast = c4*(-1)*phi;
% F_elast = 0;
% M_elast = 0;

% Force balance system
force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)           
geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
geometry_constraint3 = (theta_gama - theta_gama0) - phi; 
geometry_constraint4 = (theta_betta - theta_betta0) - phi;
extra_constraint1 = P*r - (T_gama*sin(abs(gama)) + T_betta*sin(abs(betta)));
extra_constraint2 = P*r - (T_gama*sin(abs(-gama)) + T_betta*sin(abs(-betta)));
% extra_constraint2 = T_gama - T_betta;       % shot in the dark, what if the tensions are equal?-DOESN'T SEEM TO WORK

F = [force_balance;...
      torque_balance;...
      geometry_constraint1;...
      geometry_constraint2;...
      geometry_constraint3;...
      geometry_constraint4;...
      extra_constraint1;...
      extra_constraint2];

% Jacobian of force balance system
  if nargout > 1   % Two output arguments
      
      J = [ pi*r^2,                                                               2*T_gama*sin(gama),                                                                  2*T_betta*sin(betta),                                  2*pi*P*r,           -c1,            0,    -2*cos(gama),    -2*cos(betta);...
          0,                                                             2*T_gama*r*cos(gama),                                                                2*T_betta*r*cos(betta), 2*T_betta*sin(betta) + 2*T_gama*sin(gama),             0,          -c4,   2*r*sin(gama),   2*r*sin(betta);...
          0, (L*sin(gama))/cos(gama)^2 - (r*cos(gama)*(phi - (L0*tan(gama0))/r0))/sin(gama)^2,                                                                                     0,      (phi - (L0*tan(gama0))/r0)/sin(gama),   1/cos(gama),  r/sin(gama),               0,                0;...
          0,                                                                                0, (L*sin(betta))/cos(betta)^2 - (r*cos(betta)*(phi - (L0*tan(betta0))/r0))/sin(betta)^2,    (phi - (L0*tan(betta0))/r0)/sin(betta),  1/cos(betta), r/sin(betta),               0,                0;...
          0,                                                         -(L*(tan(gama)^2 + 1))/r,                                                                                     0,                         (L*tan(gama))/r^2,  -tan(gama)/r,           -1,               0,                0;...
          0,                                                                                0,                                                             -(L*(tan(betta)^2 + 1))/r,                        (L*tan(betta))/r^2, -tan(betta)/r,           -1,               0,                0;...
          r,                                                -T_gama*cos(abs(gama))*sign(gama),                                                  -T_betta*cos(abs(betta))*sign(betta),                                         P,             0,            0, -sin(abs(gama)), -sin(abs(betta));...
          r,                                                -T_gama*cos(abs(gama))*sign(gama),                                                  -T_betta*cos(abs(betta))*sign(betta),                                         P,             0,            0, -sin(abs(gama)), -sin(abs(betta))];
  
  end
  

end