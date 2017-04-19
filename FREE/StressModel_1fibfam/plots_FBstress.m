% plots_FBstress.m
%   Steady State Iterative Simulation
%   Define FREE parameters
%   Iterate through pressures and calculate states at each one
%
%   NOTE: Must run setParams.m before using this script

%% Number of points at which to evaluate the steady state behavior
res = 500; 

%% Simulates steady state behavior of FREE for pressures from Pmin-Pmax psi
x(1,:) = params.x_rest + [0, 0, 1e-6, 0, 0, 0];    % initial simulation point (slightly more than zero so solver doesn't complain)
x1bolt(1,:) = params.x_rest + [0, 0, 1e-6, 0, 0, 0];    % initial simulation point (slightly more than zero so solver doesn't complain)
x2bolt(1,:) = params.x_rest + [0, 0, 1e-6, 0, 0, 0];    % initial simulation point (slightly more than zero so solver doesn't complain) 


% Mload = Fload = 0
params.load = [0, 0];        % [F_load, M_load];
for i = 1:res
    P_des(i+1) = params.Pmin + ((params.Pmax-params.Pmin)/res) * i;
    x0 = x(i,:);    %set initial guess of solution as state at last P_des
    x(i+1, :)  = solve_FBstress(P_des(i+1), params, x0);
end

% Mload = 1 bolt, Fload = 0
params.load = [0, -(0.117933525 + 0.0)];        % [F_load, M_load];
for j = 1:res
    P_des(j+1) = params.Pmin + ((params.Pmax-params.Pmin)/res) * j;
    x0 = x1bolt(j,:);    %set initial guess of solution as state at last P_des
    x1bolt(j+1,:) = solve_FBstress(P_des(j+1), params, x0);
end

% Mload = 2 bolts, Fload = 0
params.load = [0, -(0.23586705 + 0.0)];        % [F_load, M_load];
for k = 1:res
    P_des(k+1) = params.Pmin + ((params.Pmax-params.Pmin)/res) * k;
    x0 = x2bolt(k,:);    %set initial guess of solution as state at last P_des
    x2bolt(k+1,:) = solve_FBstress(P_des(k+1), params, x0);
end


%% Plots of steady state values at different pressures
% figure
% plot(x(:,1), rad2deg(x(:,2)))
% title('SS fiber angle vs pressure')
% 
% figure
% plot(x(:,1), x(:,3))
% title('SS radius vs pressure')
% 
% figure
% hold on
% plot(x(:,1), x(:,4) - L_rest)
% title('SS length change vs pressure')
% plot(P, dL,'*')
% hold off
% 
% figure
% hold on
% plot(x(:,1), rad2deg(x(:,5)))
% title('SS FREE twist vs pressure')
% plot(P, rad2deg(dphi),'*')
% hold off

%% Generates Rotation vs. Pressure Plots Suitable for publication

% Convert psi to kPa
Psim_kPa = x(:,1) * 6.89476;
Psim_kPa_1bolt = x1bolt(:,1) * 6.89476;
Psim_kPa_2bolt = x2bolt(:,1) * 6.89476;

Pemp_kPa = P * 6.89476;
Pemp_kPa_1bolt = Pvphi_1bolt(:,1) * 6.89476;
Pemp_kPa_2bolt = Pvphi_2bolt(:,1) * 6.89476;


% plot it
figure
set(gca,'fontsize', 32);
%axis([5, 12*6, -100, 150]);

hold on
plot(Psim_kPa, rad2deg(x(:,5)),'k')
plot(Psim_kPa_1bolt, rad2deg(x1bolt(:,5)),'r')
plot(Psim_kPa_2bolt, rad2deg(x2bolt(:,5)),'b')

plot(Pemp_kPa, rad2deg(dphi),'*k')
plot(Pemp_kPa_1bolt, rad2deg(Pvphi_1bolt(:,2)),'*r')
plot(Pemp_kPa_2bolt, rad2deg(Pvphi_2bolt(:,2)),'*b')


%legend('No Load','Load = 9.73e-3 Nm','Load = 19.46e-3 Nm','Load = 29.61e-3 Nm','Location','southeast');
xlabel('Pressure, $P$ (kPa)', 'fontsize', 32, 'Interpreter', 'LaTex');
ylabel('Rotation, $\Phi$ (degrees)', 'fontsize', 32, 'Interpreter', 'LaTex');
hold off

%% Caluclate average absolute error for each loading case
fudge = 0.01;   % tolerance for equality comparisons

% No load error
for j = 1:length(P)
    for k = 1:length(x(:,1))
        if P(j) < x(k,1)+fudge && P(j) > x(k,1)-fudge
            error_noload(j) = abs(rad2deg(dphi(j) - x(k,5)));
        end
    end
end  
abserror_noload = sum(error_noload) / length(P);


% 1 bolt load error
for j = 1:length(Pvphi_1bolt(:,1))
    for k = 1:length(x1bolt(:,1))
        if Pvphi_1bolt(j,1) < x1bolt(k,1)+fudge && Pvphi_1bolt(j,1) > x1bolt(k,1)-fudge
            error_1bolt(j) = abs(rad2deg(Pvphi_1bolt(j,2) - x1bolt(k,5)));
        end
    end
end  
abserror_1bolt = sum(error_1bolt) / length(Pvphi_1bolt(:,1));


% 2 bolt load error
for j = 1:length(Pvphi_2bolt(:,1))
    for k = 1:length(x2bolt(:,1))
        if Pvphi_2bolt(j,1) < x2bolt(k,1)+fudge && Pvphi_2bolt(j,1) > x2bolt(k,1)-fudge
            error_2bolt(j) = abs(rad2deg(Pvphi_2bolt(j,2) - x2bolt(k,5)));
        end
    end
end  
abserror_2bolt = sum(error_2bolt) / length(Pvphi_2bolt(:,1));

% Put all errors into a matrix so it's easy to read later
error = [0, abserror_noload;...
         1, abserror_1bolt;...
         2, abserror_2bolt];
     
maxerror = [0, max(error_noload);...
            1, max(error_1bolt);...
            2, max(error_2bolt)];

