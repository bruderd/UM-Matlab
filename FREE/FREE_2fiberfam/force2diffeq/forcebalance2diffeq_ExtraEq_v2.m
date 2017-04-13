clear all;

syms r0 L0 gama0 dP dgama dr dL c1 c2 c3 c4 c5 c6 betta0 nrat T_gama T_betta

assume(r0, 'real')
assumeAlso(L0, 'real')
assumeAlso(gama0, 'real')
assumeAlso(betta0, 'real')
assumeAlso(dP, 'real')
assumeAlso(dgama, 'real')
assumeAlso(dr, 'real')
assumeAlso(dL, 'real')
assumeAlso(c1, 'real')
assumeAlso(c2, 'real')
assumeAlso(c3, 'real')
assumeAlso(c4, 'real')
assumeAlso(c5, 'real')
assumeAlso(c6, 'real')
assumeAlso(nrat, 'real')


syms t

P = sym('P(t)');
gama = sym('gama(t)');
betta = sym('betta(t)');
r = sym('r(t)');
L = sym('L(t)');
phi = sym('phi(t)');
T_gama = sym('T_gama(t)');
T_betta = sym('T_betta(t)');

assumeAlso(P, 'real')
assumeAlso(gama, 'real')
assumeAlso(betta, 'real')
assumeAlso(r, 'real')
assumeAlso(L, 'real')
assumeAlso(phi, 'real')
assumeAlso(T_gama, 'real')
assumeAlso(T_betta, 'real')

%% Definition of elastomer spring force functions, constants from sys id experimental data
% F_elast = [c1 c2 c3] * [L^2, L, 1]';
% M_elast = [c4 c5 c6] * [theta^2, theta, 1]';

%simpler version
F_elast = c1*(L0-L);    
M_elast = c4 * (-1) * phi;    % removed (-1) becasue I changed phi up there (1/26/2017). Put it back (1/30/2017)

% %simplest version
% F_elast = 0;   
% M_elast = 0;  


%% Take time derivatives of force balance and geometry equations

theta_gama0 = -tan(gama0)*L0/r0;     % (-) fixes sign convention (1/28/2017)
theta_betta0 = -tan(betta0)*L0/r0;   % (-) fixes sign convention (1/28/2017)
theta_gama = -tan(gama)*L/r;       % (-) fixes sign convention (1/28/2017)
theta_betta = -tan(betta)*L/r;       % (-) fixes sign convention (1/28/2017)

% Case 1: gama, betta > 0
force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)           
geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
geometry_constraint3 = (theta_gama - theta_gama0) - phi; 
geometry_constraint4 = (theta_betta - theta_betta0) - phi;
extra_constraint1 = 2*pi*P*r^2 - (tan(gama)*T_gama*sin(gama) + tan(betta)*T_betta*sin(betta));

f_case1 = diff([force_balance; torque_balance; geometry_constraint1; geometry_constraint2; geometry_constraint3; geometry_constraint4; extra_constraint1], t);

% Case 2: gama, betta < 0
force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)           
geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
geometry_constraint3 = (theta_gama - theta_gama0) - phi; 
geometry_constraint4 = (theta_betta - theta_betta0) - phi;
extra_constraint1 = 2*pi*P*r^2 - (tan(-gama)*T_gama*sin(-gama) + tan(-betta)*T_betta*sin(-betta));

f_case2 = diff([force_balance; torque_balance; geometry_constraint1; geometry_constraint2; geometry_constraint3; geometry_constraint4; extra_constraint1], t);

% Case 3: gama > 0, betta < 0
force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)           
geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
geometry_constraint3 = (theta_gama - theta_gama0) - phi; 
geometry_constraint4 = (theta_betta - theta_betta0) - phi;
extra_constraint1 = 2*pi*P*r^2 - (tan(gama)*T_gama*sin(gama) + tan(-betta)*T_betta*sin(-betta));

f_case3 = diff([force_balance; torque_balance; geometry_constraint1; geometry_constraint2; geometry_constraint3; geometry_constraint4; extra_constraint1], t);

% Case 4: gama < 0, betta > 0
force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)           
geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
geometry_constraint3 = (theta_gama - theta_gama0) - phi; 
geometry_constraint4 = (theta_betta - theta_betta0) - phi;
extra_constraint1 = 2*pi*P*r^2 - (tan(-gama)*T_gama*sin(-gama) + tan(betta)*T_betta*sin(betta));

f_case4 = diff([force_balance; torque_balance; geometry_constraint1; geometry_constraint2; geometry_constraint3; geometry_constraint4; extra_constraint1], t);


