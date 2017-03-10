% sim_fsolve.m


clear

%% Set desired parameter values
Pss = 10;     % steady state pressure (input)

% Relaxed parameters (when P = 0)
P_rest = 0;
gama_rest = deg2rad(40);
betta_rest = deg2rad(-20); 
r_rest = 3/16;
L_rest = 5;
phi_rest = 0;
T_gama_rest = 0;
T_betta_rest = 0;
x_rest = [P_rest, gama_rest, betta_rest, r_rest, L_rest, phi_rest, T_gama_rest, T_betta_rest];

N = 1000;        % number of pressure steps to reach steady state value
dP = (Pss/N);   % pressure step size


%% Iteratively solve for state at each pressure step

x = zeros(N+1,8);
error = zeros(1,length(x(:,1)));
x(1,:) = x_rest;
dx = 1;
dtens = 0.01;
delta = 1.0;

for k = 1:N
    u = dP*k;
    
%     options = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective','Display', 'iter', 'TolFun', 1e-10, 'TolX', 1e-10);
    [P0, gama0, betta0, r0, L0, phi0, T_gama0, T_betta0] = deal(x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), x(k,7), x(k,8));
    x0 = [P0, gama0, betta0, r0, L0, phi0, T_gama0, T_betta0];
    x_prev = x(k,:);
    
    
%     % Solve for the tension forces with all other parameters fixed
%     tensioneq = @(x)tensioneq_2fib_v2(x, u, x_rest);
%     LB_tens = [u-dx, gama0-dx, betta0-dx, r0, L0-dx, phi0-dx, 0, 0];
%     UB_tens = [u+dx, gama0+dx, betta0+dx, r0+dx, L0+dx, phi0+dx, Inf, Inf];
% %    UB_tens = [u+dx, gama0+dx, betta0+dx, r0+dx, L0+dx, phi0+dx, T_gama0+delta, T_betta0+delta];
%     tens = lsqnonlin(tensioneq, x0, LB_tens, UB_tens);
%     
    % Solve for steady state behavior with the tension forces
%    options = optimoptions('lsqnonlin','Jacobian','on', 'Display', 'iter' );
%     fun = @(x)staticeq_2fib_v2(x, u, x_rest);
    fun = @(x)FB_v2(x, u, x_prev, x_rest);
%     T_gama = tens(7);
%     T_betta = tens(8);
%     LB = [u-dx, -pi/2, -pi/2, r0*1.0, L0*0.4, -Inf, T_gama-dtens, T_betta-dtens];
%     UB = [u+dx, pi/2, pi/2, r0*5, L0*2, Inf, T_gama+dtens, T_betta+dtens];

    % Try bounds that impose small changes between time/pressure steps
%     LB = [u-dx, gama0-delta, betta0-delta, r0*1.0, L0-delta, phi0-delta, T_gama-dtens, T_betta-dtens];
%     UB = [u+dx, gama0+delta, betta0+delta, r0*(1+delta), L0+delta, phi0+delta, T_gama+dtens, T_betta+dtens];
% %    UB = [u+dx, gama0+delta, betta0+delta, r_rest*3, L0+delta, phi0+delta, T_gama+dtens, T_betta+dtens];

    options = optimoptions('fsolve', 'Display', 'iter' );
    x(k+1,:) = fsolve(fun, x0, options);
    
    % Calculate error
    Fval = FB_v2(x(k+1,:), u, x_prev, x_rest);
    error(k) = norm(Fval.^2, Inf);
    
end

maxerror = max(error);

%% Check magnitude of the force balance equations at each step
% error = zeros(1,length(x));
% for j = 1:N
%     u = dP*j;
%     
%     Fval = FB_grad_v2(x(j+1,:), u, x_rest);
%     error(j) = norm(Fval.^2, Inf);
% end
% 
% maxerror = max(error);


%% Plot results

figure
set(gcf,'numbertitle','off','name','Steady State Iterative Simulation Results') % See the help for GCF

subplot(2,3,1)       
plot(x(:,1))
title('Pressure (psi)')

subplot(2,3,2)       
plot(rad2deg(x(:,2)))
title('gamma (deg)')

subplot(2,3,3)       
plot(rad2deg(x(:,3)))
title('beta (deg)')

subplot(2,3,4)       
plot(x(:,4))
title('radius (in)')

subplot(2,3,5)       
plot(x(:,5))
title('Length (in)')

subplot(2,3,6)       
plot(rad2deg(x(:,6)))
title('phi (deg)')

