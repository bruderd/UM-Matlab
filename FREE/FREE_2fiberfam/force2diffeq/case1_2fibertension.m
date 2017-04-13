%% Case 1: gama,betta>0, gama>betta

clear
clc

syms P gama betta r L T_gama T_betta nrat

Fup_gama = 4*pi*P*r^2*cot(gama);
Fup_betta = 4*pi*P*r^2*cot(betta);
eq1 = Fup_betta - 2*T_betta*sin(betta) - (nrat+1)*T_gama*sin(gama) + T_gama*sin(gama)*cos(pi*(tan(gama)/tan(betta) - nrat));
eq2 = Fup_gama - 2*T_gama*sin(gama) - T_betta*sin(betta) + T_betta*sin(betta)*cos(pi*tan(betta)/tan(gama));

solution = solve(0 == eq1, 0 == eq2, T_gama, T_betta);