% VSAsim_ode45_test.m
%   Must run setParams.m before this script

[B, N, L, I, rho] = deal(params.B, params.N, params.L, params.I, params.rho);

%% Simulate dynamics using ode45

u = [1000; 2000]; % set the value of the input u = [P1, P2]'
tau = 0;    % set value external torque being applied to joint

tspan = [0 20];
y0 = [0, 0];
[t,y] = ode45(@(t,y) VSA_dynamics(t,y,u,tau,params), tspan, y0);

%% Calculate the stiffness, K
theta = y(:,1);

K = abs((-3*rho^2)/(2*pi*N^2) * ( (L + rho*theta)*u(1) +  (L - rho*theta)*u(2) ));

%% Caluculate the analytic steady state solution

[u1, u2] = deal(u(1), u(2));

% % With extraaneous solution
% theta_ss1 = ( 3*L*(u2+u1) + sqrt(3)*sqrt( B^2*(u2-u1)^2 + 12*L^2*u2*u1 ) ) / (3*rho*(u2-u1));
% theta_ss2 = ( 3*L*(u2+u1) - sqrt(3)*sqrt( B^2*(u2-u1)^2 + 12*L^2*u2*u1 ) ) / (3*rho*(u2-u1));
% xss = [ones(size(t)) * theta_ss1, ones(size(t)) * theta_ss2];

% Without extraaneous solution
if u1 ~= u2
    theta_ss = ( 3*L*(u2+u1) - sqrt(3)*sqrt( B^2*(u2-u1)^2 + 12*L^2*u2*u1 ) ) / (3*rho*(u2-u1));
else
    theta_ss = 0;
end
xss = ones(size(t)) * theta_ss;

%% plot the results

figure
hold on
plot(t,y(:,1),'-o',t,y(:,2),'-o')
plot(t,xss)
plot(t,K)
hold off
title('Solution of VSA dynamics with ODE45');
xlabel('Time, t');
ylabel('State, x');
h = legend('$\theta$','$\dot{\theta}$','$\theta_{ss}$','$|K|$');
set(h,'Interpreter','latex');

