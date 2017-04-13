clear

syms P gama betta r L T_gama T_betta nrat


%% ONLY 2 CASES 
% gama, betta have same sign
Fup_gama = 4*pi*P*r^2*cot(gama);
Fup_betta = 4*pi*P*r^2*cot(betta);
psi = pi*(tan(gama)/tan(betta) - nrat);
eq_lgama = Fup_gama - 2*T_gama*sin(gama) - T_betta*sin(betta)*(1 - cos(pi*tan(betta)/tan(gama)));
eq_lbetta = Fup_betta - 2*T_betta*sin(betta) - T_gama*sin(gama)*((2*nrat+1) - cos(psi));
system_same = [eq_lgama; eq_lbetta];

% tension_same = solve(0 == eq_lbetta, 0 == eq_lgama, T_gama, T_betta);

% gama, betta have different sign
Fup_gama2 = 4*pi*P*r^2*cot(gama);
Fup_betta2 = 4*pi*P*r^2*cot(-betta);
psi2 = pi*(tan(gama)/tan(-betta) - nrat);
eq_lgama2 = Fup_gama2 - 2*T_gama*sin(gama) - T_betta*sin(-betta)*(1 - cos(pi*tan(-betta)/tan(gama)));
eq_lbetta2 = Fup_betta2 - 2*T_betta*sin(-betta) - T_gama*sin(gama)*((2*nrat+1) - cos(psi2));
system_dif = [eq_lgama2; eq_lbetta2];

% tension_dif = solve(0 == eq_lbetta, 0 == eq_lgama, T_gama, T_betta);


%% Plug values in for gama, betta

Pnew = 1;
rnew = 1;
Lnew = 4;
gam = deg2rad(40);
bet = deg2rad(-38);
nr = floor(abs(tan(gam)/tan(bet)));

poop_same = subs(system_same, [P, gama, betta, r, L, nrat], [Pnew, gam, bet, rnew, Lnew, nr]);
poop_dif = subs(system_dif, [P, gama, betta, r, L, nrat], [Pnew, gam, bet, rnew, Lnew, nr]);

