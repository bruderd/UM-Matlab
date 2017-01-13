%% solve_ICRA

clear

P = 10;     % steady state pressure

fun = @(x)static_balance_ICRA(P,x);

gama0 = 0.7037;    
r0 = 3/16;    
L0 = 5.68;

x0 = [gama0,r0,L0];
x_sol = lsqnonlin(fun,x0,[-pi/2,r0,1],[pi/2,1,5]);


%% Plug solution back into equations to check:

F = static_balance_ICRA(P,x_sol);
