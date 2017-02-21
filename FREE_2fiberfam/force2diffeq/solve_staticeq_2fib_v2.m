%% Solve system of static equations
%    New states added: T_gama, T_betta
%    Uses 8 dimensional set of force balance equations
%    Uses nonlinear solver

clear


%% Solve steady state equation

% Set desired parameter values
u = 5;     % steady state pressure (input)
P0 = 0;
gama0 = deg2rad(40);
betta0 = deg2rad(-40); 
r0 = 3/16;
L0 = 5.62;
phi0 = 0;
T_gama0 = 0;
T_betta0 = 0;
x0 = [P0, gama0, betta0, r0, L0, phi0, T_gama0, T_betta0];

fun = @(x)staticeq_2fib_v2(x, u, x0);
LB = [u-1, -pi/2, -pi/2, r0*1.0, L0*0.4, -Inf, 0, 0];
UB = [u+1, pi/2, pi/2, r0*5, L0*2, Inf, Inf, Inf];

sol = lsqnonlin(fun, x0, LB, UB);

P = sol(1);
gama = sol(2);
betta = sol(3);
r = sol(4);
L = sol(5);
phi = sol(6);





