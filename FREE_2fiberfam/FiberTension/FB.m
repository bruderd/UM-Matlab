% FB.m
%   Force balance equations

function FB = FB(P_test, T, x, x_rest)

P = P_test;
[T_gama, T_betta] = deal(T(1), T(2));
[P, gama, betta, r, L, phi] = deal(x(1), x(2), x(3), x(4), x(5), x(6));
[P0, gama0, betta0, r0, L0, phi0] = deal(x_rest(1), x_rest(2), x_rest(3), x_rest(4), x_rest(5), x_rest(6));

theta_gama0 = -tan(gama0)*L0/r0;     
theta_betta0 = -tan(betta0)*L0/r0;   
theta_gama = -tan(gama)*L/r; 
theta_betta = -tan(betta)*L/r; 

F_elast = 0;
M_elast = 0;

force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)           
geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
geometry_constraint3 = (theta_gama - theta_gama0) - phi; 
geometry_constraint4 = (theta_betta - theta_betta0) - phi;
extra_constraint1 = P*r - (T_gama*sin(abs(gama)) + T_betta*sin(abs(betta)));
extra_constraint2 = P*r - (T_gama*sin(abs(-gama)) + T_betta*sin(abs(-betta)));

FB = [force_balance;...
      torque_balance;...
      geometry_constraint1;...
      geometry_constraint2;...
      geometry_constraint3;...
      geometry_constraint4;...
      extra_constraint1;...
      extra_constraint2];
    

end