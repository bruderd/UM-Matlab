% sim_FBstress.m
%
% NOTE: Must run setParams.m before using this script


% Set desired input values
Pss = 10;     % steady state pressure (input)

N = 100;        % number of pressure steps to reach steady state value
dP = (Pss/N);   % pressure step size

%% Iteratively solve for state at each pressure step

for k = 1:N
    u = dP*k;
    [P(k), gama(k), r(k), L(k), phi(k), T(k)] = solve_FBstress(u, params);
end


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