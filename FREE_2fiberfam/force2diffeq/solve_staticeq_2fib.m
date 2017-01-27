%% Solve system of static equations
clear

% syms r0 L0 gama0 dP dgama dr dL c1 c2 c3 c4 c5 c6 betta0 P gama betta r L phi
% 
% assume(r0, 'real')
% assumeAlso(L0, 'real')
% assumeAlso(gama0, 'real')
% assumeAlso(betta0, 'real')
% assumeAlso(dP, 'real')
% assumeAlso(dgama, 'real')
% assumeAlso(dr, 'real')
% assumeAlso(dL, 'real')
% assumeAlso(c1, 'real')
% assumeAlso(c2, 'real')
% assumeAlso(c3, 'real')
% assumeAlso(c4, 'real')
% assumeAlso(c5, 'real')
% assumeAlso(c6, 'real')
% 
% assumeAlso(P, 'real')
% assumeAlso(gama, 'real')
% assumeAlso(betta, 'real')
% assumeAlso(r, 'real')
% assumeAlso(L, 'real')
% assumeAlso(phi, 'real')
% 
% %simpler version
% F_elast = c1*(L0-L);    
% M_elast = c4 * phi;    % removed (-1) becasue I changed phi up there (1/26/2017)

%% Case 1: gama,betta>0, gama>betta. Differentiates the system of equations wrt time (2 fibers)
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

sol = lsqnonlin(fun,x0,[9,gama0*0.5,betta0*0.5,r0*1.1,L0*0.5,-Inf],[11,gama0*1.5,betta0*1.5,r0*5,L0*1.5,Inf]);



%% Case 1: gama,betta>0, gama>betta. Differentiates the system of equations wrt time (2 fibers)
%    New state added: phi = "twist"

% P = 0;
% gama0 = deg2rad(40);
% betta0 = deg2rad(30);
% r0 = 3/16;
% L0 = 5.62;
%     
% T_gama = (2*(pi*P*r^2*cot(gama) + P*r^2*pi*cos((pi*tan(betta))/tan(gama))*cot(betta)))/(sin(gama)*(cos((pi*tan(betta))/tan(gama)) + 1));
% T_betta = (2*(pi*P*r^2*cot(betta) - pi*P*r^2*cot(gama)))/(sin(betta) + cos((pi*tan(betta))/tan(gama))*sin(betta));
% theta_gama0 = tan(gama0)*L0/r0;
% theta_betta0 = tan(betta0)*L0/r0;
% 
% force_balance = P*pi*r^2 - (T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
% torque_balance = (T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;               
% geometry_constraint1 = L/cos(gama) - r*(theta_gama0 + phi)/sin(gama);
% geometry_constraint2 = L/cos(betta) - r*(theta_betta0 + phi)/sin(betta);
% geometry_constraint3 = (L/r)*tan(gama) - (L0/r0)*tan(gama0) - phi;
% 
% sol = solve(0 == force_balance, 0 == torque_balance, 0 == geometry_constraint1, 0 == geometry_constraint2, 0 == geometry_constraint3, gama, betta, r, L, phi);

%% Case 1: gama,betta>0, gama>betta. Differentiates the system of equations wrt time (2 fibers)
%     % New state added: phi = "twist"
% 
% P = 0;
% gama0 = deg2rad(40);
% betta0 = deg2rad(30);
% r0 = 3/16;
% L0 = 5.62;
% c1 = 1;
% c2 = 1;
% 
% T_gama = (2*(pi*P*r^2*cot(gama) + P*r^2*pi*cos((pi*tan(betta))/tan(gama))*cot(betta)))/(sin(gama)*(cos((pi*tan(betta))/tan(gama)) + 1));
% T_betta = (2*(pi*P*r^2*cot(betta) - pi*P*r^2*cot(gama)))/(sin(betta) + cos((pi*tan(betta))/tan(gama))*sin(betta));
% theta_gama = tan(gama)*L/r;
% theta_betta = tan(betta)*L/r;
% theta_gama0 = tan(gama0)*L0/r0;
% theta_betta0 = tan(betta0)*L0/r0;
% 
% force_balance = P*pi*r^2 - (T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
% torque_balance = (T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;               
% geometry_constraint1 = L/cos(gama) - L0/cos(gama0);
% geometry_constraint2 = L/cos(betta) - L0/cos(betta0);
% % geometry_constraint3 = (L/r)*tan(gama) - (L0/r0)*tan(gama0) + phi;
% 
% sol = solve(0 == force_balance, 0 == torque_balance, 0 == geometry_constraint1, 0 == geometry_constraint2, gama, betta, r, L);
% 

%% Simpler. Doesn't use complicated tension expressions
% 
% % P = 0;
% % gama0 = deg2rad(40);
% % betta0 = deg2rad(-40);
% % r0 = 3/16;
% % L0 = 5.62;
% % c1 = 7;
% % c2 = 1;
% 
% force_balance = P*pi*r^2 - P*pi*r^2*(cot(gama)^2 + cot(betta)^2) + F_elast;   
% torque_balance = P*pi*r^3*(cot(gama) + cot(betta)) + M_elast;               
% geometry_constraint1 = -cos(gama) + (L/L0)*cos(gama0);
% geometry_constraint2 = -cos(betta) + (L/L0)*cos(betta0);
% 
% sol = solve(force_balance, torque_balance, geometry_constraint1, geometry_constraint2, gama, betta, r, L);



