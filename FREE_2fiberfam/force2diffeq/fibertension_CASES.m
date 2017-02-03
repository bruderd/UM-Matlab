%% Case 1: gama,betta>0

clear
clc

syms P gama betta r L T_gama T_betta nrat

%% Case 1: gama,betta>0
Fup_gama1 = 4*pi*P*r^2*cot(gama);
Fup_betta1 = 4*pi*P*r^2*cot(betta);
psi1 = pi*(tan(gama)/tan(betta) - nrat);
eq_lgama1 = Fup_gama1 - 2*T_gama*sin(gama) - T_betta*sin(betta)*(1 - cos(pi*tan(betta)/tan(gama)));
eq_lbetta1 = Fup_betta1 - 2*T_betta*sin(betta) - T_gama*sin(gama)*((nrat+1) - cos(psi1));

case1 = solve(0 == eq_lbetta1, 0 == eq_lgama1, T_gama, T_betta);

%% Case 2: gama,betta<0
Fup_gama2 = 4*pi*P*r^2*cot(-gama);
Fup_betta2 = 4*pi*P*r^2*cot(-betta);
psi2 = pi*(tan(-gama)/tan(-betta) - nrat);
eq_lgama2 = Fup_gama2 - 2*T_gama*sin(-gama) - T_betta*sin(-betta)*(1 - cos(pi*tan(-betta)/tan(-gama)));
eq_lbetta2 = Fup_betta2 - 2*T_betta*sin(-betta) - T_gama*sin(-gama)*((nrat+1) - cos(psi1));

case2 = solve(0 == eq_lbetta2, 0 == eq_lgama2, T_gama, T_betta);

%% Case 3: gama>0, betta<0
Fup_gama3 = 4*pi*P*r^2*cot(gama);
Fup_betta3 = 4*pi*P*r^2*cot(-betta);

eq_lgama3 = Fup_gama3 - 2*T_gama*sin(gama);
eq_lbetta3 = Fup_betta3 - 2*T_betta*sin(-betta) - T_gama*sin(gama)*(nrat);

case3 = solve(0 == eq_lbetta3, 0 == eq_lgama3, T_gama, T_betta);

%% Case 4: gama<0, betta>0
Fup_gama4 = 4*pi*P*r^2*cot(-gama);
Fup_betta4 = 4*pi*P*r^2*cot(betta);

eq_lgama4 = Fup_gama4 - 2*T_gama*sin(-gama);
eq_lbetta4 = Fup_betta4 - 2*T_betta*sin(betta) - T_gama*sin(-gama)*(nrat);

case4 = solve(0 == eq_lbetta4, 0 == eq_lgama4, T_gama, T_betta);

%% Case 5: gama,betta>0
Fup_gama5 = 4*pi*P*r^2*cot(gama);
Fup_betta5 = 4*pi*P*r^2*cot(betta);

eq_lgama5 = Fup_gama5 - 2*T_gama*sin(gama) - T_betta*sin(betta)*(1-cos(pi*tan(betta)/tan(gama)));
eq_lbetta5 = Fup_betta5 - 2*T_betta*sin(betta) - T_gama*sin(gama)*(nrat + 1);

case5 = solve(0 == eq_lbetta5, 0 == eq_lgama5, T_gama, T_betta);

%% Case 6: gama,betta<0
Fup_gama6 = 4*pi*P*r^2*cot(-gama);
Fup_betta6 = 4*pi*P*r^2*cot(-betta);

eq_lgama6 = Fup_gama6 - 2*T_gama*sin(-gama) - T_betta*sin(-betta)*(1-cos(pi*tan(-betta)/tan(-gama)));
eq_lbetta6 = Fup_betta6 - 2*T_betta*sin(-betta) - T_gama*sin(-gama)*(nrat + 1);

case6 = solve(0 == eq_lbetta6, 0 == eq_lgama6, T_gama, T_betta);

%% Case 7: gama>0, betta<0
Fup_gama7 = 4*pi*P*r^2*cot(gama);
Fup_betta7 = 4*pi*P*r^2*cot(-betta);
psi7 = pi*(tan(gama)/tan(-betta) - nrat);
eq_lgama7 = Fup_gama7 - 2*T_gama*sin(gama);
eq_lbetta7 = Fup_betta7 - 2*T_betta*sin(-betta) - T_gama*sin(gama)*(nrat - cos(psi7));

case7 = solve(0 == eq_lbetta7, 0 == eq_lgama7, T_gama, T_betta);

%% Case 8: gama<0, betta>0
Fup_gama8 = 4*pi*P*r^2*cot(-gama);
Fup_betta8 = 4*pi*P*r^2*cot(betta);
psi8 = pi*(tan(-gama)/tan(betta) - nrat);
eq_lgama8 = Fup_gama8 - 2*T_gama*sin(-gama);
eq_lbetta8 = Fup_betta8 - 2*T_betta*sin(betta) - T_gama*sin(-gama)*(nrat - cos(psi8));

case8 = solve(0 == eq_lbetta8, 0 == eq_lgama8, T_gama, T_betta);


