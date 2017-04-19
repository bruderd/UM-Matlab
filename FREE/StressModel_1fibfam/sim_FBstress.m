% sim_FBstress.m
%
% NOTE: Must run setParams.m before using this script


% Set desired input values
Pss = 5;     % steady state pressure (input)

N = 500;        % number of pressure steps to reach steady state value
dP = (Pss/N);   % pressure step size

%% Iteratively solve for state at each pressure step

x = zeros(N,length(params.x_rest));
x(1,:) = params.x_rest;

for k = 2:N+1
    u = dP*(k-1);
    x0 = x(k-1,:);  % x at last time step, can be used as IC of next step
    x(k,:) = solve_FBstress(u, params, x0);  % remove x0 if you want IC to be x_rest
    if( norm( FBstress(x(k,:), u, params) ) > 1e-3 )
        norm( FBstress(x(k,:), u, params) );
    end
end

%% Pull out parts of x
P = x(:,1);
gama = x(:,2);
r = x(:,3);
L = x(:,4);
phi = x(:,5);
T = x(:,6);

%% Plot results

figure
set(gcf,'numbertitle','off','name','Steady State Iterative Simulation Results') % See the help for GCF

subplot(2,3,1)       
plot(P)
title('Pressure (psi)')

subplot(2,3,2)       
plot(rad2deg(gama))
title('gamma (deg)')

subplot(2,3,3)       
plot(r)
title('radius (in)')

subplot(2,3,4)       
plot(L)
title('Length (in)')

subplot(2,3,5)       
plot(rad2deg(phi))
title('phi (deg)')

subplot(2,3,6)       
plot(T)
title('Tension (lbf)')