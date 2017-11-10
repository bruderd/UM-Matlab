% Pendulum_demo  This demo illustrates the usage of the EulerLagrange tool
%                from beginning to end. The system Lagrangian is defined,
%                the tool creates the corresponding MATLAB function, which 
%                is integrated using ode45. The plots show the state vector 
%                and the total energy of the system.
%
% Copyright 2015-2016 The MathWorks, Inc.

% Derive DE for pendulum
syms th thd m l g
T   = m*l^2*thd.^2/2;
V   = m*g*l*(1 - cos(th));
L   = T - V;
E   = T + V;
X   = {th thd};
Q_i = {0}; Q_e = {0};
R   = 0;
par = {m g l};
% Solve Lagrange equations and save DE as .m file
VF  = EulerLagrange(L,X,Q_i,Q_e,R,par,'m','Pendulum_sys');

% Solve DE numerically using ode45
m_num   = 1;
g_num   = 9.81;
l_num   = 1;
tspan   = [0 10];
Y0      = [0 30]*pi/180;
options = odeset('RelTol',1e-8);
% Use created .m file to solve DE 
[t, Y]  = ode45(@Pendulum_sys,tspan,Y0,options,m_num,g_num,l_num);

% Plot state vector and total energy
subplot(2,1,1)
plot(t,Y)
title('Pendulum')
xlabel('t')
ylabel('{\theta}(t), d{\theta}/dt(t)')
subplot(2,1,2)
E_num = double(subs(E,{th,thd,m,g,l},{Y(:,1),Y(:,2),m_num,g_num,l_num}));
plot(t,E_num)
xlabel('t')
ylabel('E = T + V')
% End of script