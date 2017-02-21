%% Solve system of static equations
%    New states added: T_gama, T_betta
%    Uses 8 dimensional set of force balance equations
%    Uses nonlinear solver

clear


%% Solve steady state equation

% Set desired parameter values
u = 10;     % steady state pressure (input)
P0 = 0;
gama0 = deg2rad(40);
betta0 = deg2rad(-30); 
r0 = 3/16;
L0 = 5;
phi0 = 0;
T_gama0 = 0;
T_betta0 = 0;
x0 = [P0, gama0, betta0, r0, L0, phi0, T_gama0, T_betta0];


% Solve for the tension forces with all other parameters fixed
tensioneq = @(x)tensioneq_2fib_v3(x, u, x0);
dx = 0.01;
LB_tens = [u-dx, gama0-dx, betta0-dx, r0, L0-dx, phi0-dx, 0, 0];
UB_tens = [u+dx, gama0+dx, betta0+dx, r0+dx, L0+dx, phi0+dx, Inf, Inf];
tens = lsqnonlin(tensioneq, x0, LB_tens, UB_tens);


% Solve for steady state behavior with the tension forces
fun = @(x)staticeq_2fib_v2(x, u, x0);
T_gama = tens(7);
T_betta = tens(8);
dtens = 0.001;
LB = [u-1, -pi/2, -pi/2, r0*1.01, L0*0.4, -Inf, T_gama-dtens, T_betta-dtens];
UB = [u+1, pi/2, pi/2, r0*5, L0*2, Inf, T_gama+dtens, T_betta+dtens];

sol = lsqnonlin(fun, x0, LB, UB);

P = sol(1);
gama = sol(2);
betta = sol(3);
r = sol(4);
L = sol(5);
phi = sol(6);





