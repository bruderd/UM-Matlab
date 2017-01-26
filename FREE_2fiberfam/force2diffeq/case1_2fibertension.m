clear
clc

syms P gama betta r L T_gama T_betta

F_up = 4*pi*P*r^2*cot(gama);
eq1 = F_up - 2*T_gama*sin(gama) - 2*T_betta*sin(betta);
eq2 = F_up - 2*T_gama*sin(gama) + 2*T_betta*sin(betta)*cos(pi*tan(betta)/tan(gama));

solution = solve(eq1, eq2, T_gama, T_betta);