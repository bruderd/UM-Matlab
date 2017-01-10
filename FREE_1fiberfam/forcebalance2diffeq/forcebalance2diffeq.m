clear all;

% syms r0 L0 gama0 P gama r L dP dgama dr dL c1 c2 c3 c4 c5 c6
syms r0 L0 gama0 dP dgama dr dL c1 c2 c3 c4 c5 c6

assume(r0, 'real')
assumeAlso(L0, 'real')
assumeAlso(gama0, 'real')
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


syms t

P = sym('P(t)');
gama = sym('gama(t)');
r = sym('r(t)');
L = sym('L(t)');
phi = (-tan(gama)*L/r + tan(gama0)*L0/r0);      

assumeAlso(P, 'real')
assumeAlso(gama, 'real')
assumeAlso(r, 'real')
assumeAlso(L, 'real')

% Definition of elastomer spring force functions, constants from sys id
% experimental data
F_elast = [c1 c2 c3] * [L^2, L, 1]';
M_elast = [c4 c5 c6] * [phi^2, phi, 1]';


% Differentiates the system of equations wrt time
force_balance = diff( 0 == P*pi*r^2 - 2*P*pi*r^2*cot(gama)^2 + F_elast ,  t);   % changed signs of first 2 terms to match ICRA
torque_balance = diff( 0 == -2*P*pi*r^3*cot(gama) + M_elast ,  t);               % changed signs back to match ICRA
geometry_constraint = diff( 0 == -cos(gama) + (L/L0)*cos(gama0) ,  t);

System = [force_balance; torque_balance; geometry_constraint];

%%
% % I've taken that system and replaced all the F(t) with just F's. That is
% % the one written below...
% % Now solve for each of the time derivatives
% S = solve([...
%     pi*r^2*dP + 2*pi*P*r*dr == 2*pi*cot(gama)^2*r^2*dP - ((c2 - c1*r)*dL)/L0 - (c1*dr*(L0 - L))/L0 + 4*pi*cot(gama)^2*P*r*dr - 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1),...
%     2*pi*cot(gama)*r^3*dP - 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) + 6*pi*cot(gama)*P*r^2*dr == (c3*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*dr)/L0 - ((c4 + c3*r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2))/L0,...
%     -sin(gama)*dgama == (cos(gama0)*dL)/L0],...
%     [dgama, dr, dL]);

% %%
% % Substitutes in the values of c1, c2, c3, c4 derived from experimental
% % data
% constants = [248.19, -17.146, 5.6038, 0.2483];
% dgama_dt = subs(S.dgama, [c1, c2, c3, c4], constants);
% dr_dt = subs(S.dr, [c1, c2, c3, c4], constants);
% dL_dt = subs(S.dL, [c1, c2, c3, c4], constants);
% 
% % converts differential equations into a matlab function
% f(dgama, dr, dL, dP) = [dgama_dt, dr_dt, dL_dt];
% g = matlabFunction(f);
% 
% 

