% fbControl_p.m
%   Simulates the system dynamics under proportional feedback control

[B, N, L, I, rho] = deal(params.B, params.N, params.L, params.I, params.rho);

%% simulate dynamics using ode45

tau = 0.5;    % set value external torque being applied to joint

tspan = [0 20];
y0 = [0, 0];
[t,y] = ode45(@(t,y) VSA_dynamics_fbControl_p(t,y,tau,params), tspan, y0);

%% Calculate the stiffness, K
theta = y(:,1);

% calculate the control input at each time step, then compute the stiffness
for i = 1:length(theta)
    
    u(:,i) = controlLaw_p(y(i,:)', params);
    
    K(i) = (-3*rho^2)/(2*pi*N^2) * ( (L + rho*theta(i))*u(1,i) +  (L - rho*theta(i))*u(2,i) );
    
end
    


%% Make vectors from constant values so that they can be plotted

theta_des = ones(size(t)) * params.theta_des;
thetadot_des = ones(size(t)) * params.thetadot_des;
K_des = ones(size(t)) * params.K_des;

x_des = [theta_des'; thetadot_des'];

%% plot the results

figure
hold on
plot(t,y(:,1),'-o',t,y(:,2),'-o')
plot(t,x_des)
hold off
title('Solution of VSA dynamics using proportional control with ODE45');
xlabel('Time, t');
ylabel('State, x');
h = legend('$\theta$','$\dot{\theta}$','$\theta_{des}$','$\dot{\theta}_{des}$');
set(h,'Interpreter','latex');

% % plot the desired stiffness with the actual stiffness
% figure
% hold on
% plot(t,K)
% plot(t,K_des)
% hold off
% h2 = legend('$K_{des}$','$K$');
% set(h2,'Interpreter','latex');
