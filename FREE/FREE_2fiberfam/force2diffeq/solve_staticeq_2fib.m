%% Solve system of static equations
clear


%% Solve steady state equation
%    New state added: phi = "twist"
%    Uses nonlinear solver

P = 10;     % steady state pressure
P0 = 0;
gama0 = deg2rad(40);
betta0 = deg2rad(30); 
r0 = 3/16;
L0 = 5.62;
phi0 = 0;
x0 = [0,gama0,betta0,r0,L0, phi0];

fun = @(x)staticeq_2fib(x,P,x0);
LB = [9,gama0*0.5,betta0*0.5,r0*1.1,L0*0.5,-Inf];
UB = [11,gama0*1.5,betta0*1.5,r0*5,L0*1.5,Inf];

sol = lsqnonlin(fun, x0, LB, UB);

P = sol(1);
gama = sol(2);
betta = sol(3);
r = sol(4);
L = sol(5);
phi = sol(6);





