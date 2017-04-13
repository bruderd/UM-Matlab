% Steady State Iterative Simulation
%   Define FREE parameters
%   Iterate through pressures and calculate states at each one

clear

params = struct;

% Set values of relaxed parameters (except gamma, that comes in sections 3)
P_rest = 0;
r_rest = 0.1875;    % 3/16 in
L_rest = 3.4;     % medium one
phi_rest = 0;

%% Measured data points at each fiber angle (no load)

% 20 deg FREE
P_20 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 10, 8.5, 7.5, 5, 6.5, 4, 2.5, 6.5, 0, 0, 0];
dL_20 = P_20*(-1e-3);
% dphi_20 = [0.5378406623, 0.9500176184, 1.352141478, 1.88998214, 2.392636965, 2.920424531, 4.272566009, 0.5679999518, 0.9550441667,...
%         1.377274219, 1.879929044, 2.513274123, 3.15164575, 4.121769562, 4.172035044, 3.26725636, 2.684176763, 1.573309601,...
%         2.060884781, 1.065628228, 0.5780530483, 1.844743206, 0, 0, 0];
dphi_20 = [0.1419999879,...
           0.2387610417, 0.3443185548, 0.469982261, 0.6283185307, 0.7879114375, 1.03044239, 1.043008761, 0.8168140899,...
           0.6710441908, 0.3933274002, 0.5152211952, 0.266407057, 0.1445132621, 0.4611858015, 0, 0, 0];

% 40 deg FREE
P_40 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 11, 8.5, 6.5, 4, 2.5, 0, 0, 0];
dL_40 = [0, 0, 0, 0, 0, 0.004, 0.02, 0.028, 0.04, 0.06, 0.012, 0.012, 0.012, 0.012, 0.012, 0.016, 0.02, 0.028, 0.044, 0.06, 0.06, 0.06, 0.048, 0.036, 0.032, 0, 0, 0];
% dphi_40 = [1.301875996, 2.045805136, 2.724389149, 3.624141285, 4.518866873, 5.292955303, 6.690335715, 7.816282522, 8.992494812, 10.24410532, 1.447645895,...
%            2.37755732, 3.086300623, 4.021238597, 4.926017281, 5.700105711, 6.976848965, 7.856494908, 9.228742579, 10.50548583, 8.384282474, 6.172601246,...
%            4.337911136, 2.714336053, 1.512991022, 0, 0, 0];
dphi_40 = [0.3254689989, 0.511451284, 0.6810972873, 0.9060353213, 1.129716718, 1.323238826, 1.672583929, 1.954070631,...
           2.248123703, 2.561026331, 0.3619114737, 0.5943893301, 0.7715751557, 1.005309649, 1.23150432, 1.425026428,...
           1.744212241, 1.964123727, 2.307185645, 2.626371458, 2.096070618, 1.543150311, 1.084477784, 0.6785840132, 0.3782477555, 0, 0, 0];

       
% 60 deg FREE
P_60 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 10, 12, 8.5, 5, 2.5, 0, 0, 0];
dL_60 = [0, 0, 0.012, 0.032, 0.052, 0.064, 0.084, 0.108, 0.132, 0.164, 0, 0, 0.012, 0.024, 0.04, 0.056, 0.084, 0.108, 0.136, 0.144, 0.136, 0.144, 0.104, 0.064, 0.024, 0, 0, 0];
% dphi_60 = [1.26166361, 2.307185645, 3.201911233, 4.423362456, 5.227610176, 6.021804798, 7.404105566, 8.33904354, 9.927432785, 11.46053, 1.613521987, 2.437875899,... 
%            3.518583772, 4.322831491, 5.162265048, 6.308318048, 7.554902013, 8.891963847, 10.47532654, 11.37507868, 8.675822272, 10.5959637, 6.936636579, 4.09663682, 2.010619298, 0, 0, 0];
dphi_60 = [0.3154159024, 0.5767964112, 0.8004778081, 1.105840614, 1.306902544, 1.5054512, 1.851026391, 2.084760885,...
           2.481858196, 2.8651325, 0.4033804967, 0.6094689748, 0.879645943, 1.080707873, 1.290566262, 1.577079512,...
           1.888725503, 2.222990962, 2.618831636, 2.84376967, 2.168955568, 2.648990926, 1.734159145, 1.024159205, 0.5026548246, 0, 0, 0];

       
% 50 deg FREE
P_50 = [0, 2.5, 4, 5, 6, 7.5, 8.5, 10, 11, 12.5, 13, 12.5, 13.5, 11, 12, 10, 11, 8.5, 10, 7.5, 9, 6.5, 8, 5.4, 6.5,...
        4, 5, 3.4, 0];
dL_50 = [0, 0, 0, 0.004, 0.008, 0.024, 0.048, 0.06, 0.084, 0.104, 0.148, 0.148, 0.164, 0.14, 0.148, 0.116,...
         0.124, 0.096, 0.1, 0.096, 0.096, 0.076, 0.076, 0.064, 0.064, 0.056, 0.056, 0.052, 0.032];
dphi_50 = [0, 0.3506017401, 0.5981592412, 0.7841415263, 0.981433545, 1.303132633, 1.559486593, 1.901291874, 2.280796267,...
           2.816123655, 3.204424507, 3.183061677, 3.507274038, 2.790990913, 3.204424507, 2.454212181,...
           2.729415697, 2.038265314, 2.408973247, 1.807044094, 2.065911329, 1.501681288, 1.771858257, 1.221451224,...
           1.476548547, 0.9437344331, 1.109610525, 0.7690618816, 0.1445132621];
       
% 40 deg FREE with sandpaper
P_40v2 = [0, 3, 4.3, 5.4, 6.3, 8, 8.8, 10.1, 11.1, 12.2, 13.7, 11, 12.3, 10, 11.2, 8.9, 10.1, 7.8, 9, 7,...
          8, 5.2, 6.6, 4, 5, 3, 0];
dL_40v2 = [0, 0, 0.004, 0.008, 0.012, 0.024, 0.036, 0.052, 0.064, 0.084, 0.112, 0.096, 0.104, 0.084, 0.084, 0.08, 0.08, 0.064, 0.064, 0.056, 0.056, 0.048,...
           0.048, 0.044, 0.044, 0.036, 0.024];
dphi_40v2 = [0, 0.3581415625, 0.5504070329, 0.736389318, 0.9022654101, 1.198831757, 1.394867138, 1.751752064, 2.00559275, 2.333575023, 2.788477639,...
             2.141309553, 2.469291826, 1.922654704, 2.173982116, 1.632371543, 1.928937889, 1.394867138, 1.642424639, 1.21893795, 1.422513154,...
             0.8984954989, 1.112123799, 0.6723008279, 0.8432034682, 0.5215043805, 0.08168140899];
       
       
%% Choose which measured values to use: 1 = 20 deg FREE, 2 = 40 deg FREE, 3 = 60 deg FREE, 4 = 50 deg FREE
%   notes:
%       rcoeff: coefficients of relationship between r and P. r = c1 + c2*P.
%       kelast: coefficients of elastomer functions [1, dL_norm, dl_norm^2,...; 1, dphi_norm, dphi_norm^2,...]

choice = 5; % MAKE YOUR CHOICE HERE

if choice == 1
    gama_rest = deg2rad(20);
    P = P_20;
    dL = dL_20;
    dphi = dphi_20;
    params.rcoeff = [r_rest, 0.0118/2];     % 20 degree FREE
    params.kelast = [-0.0010e3, -7.0179e3, 0, 0, 0; -0.0217, -7.5372, 0, 0, 0];      % 20 deg FREE
%     params.kelast = [-0.0010e3, -7.0179e3, 0, 0, 0; 0, 0, 0, 0, 0];      % no elastomer (rotation)
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('load20deg_v2.csv',2,46);
elseif choice == 2
    gama_rest = deg2rad(40);
    P = P_40;
    dL = dL_40;
    dphi = dphi_40;
    params.rcoeff = [r_rest, 0.0038/2];     % 40 degree FREE
    params.kelast = [0.9120, 147.5603, 0, 0, 0; 0.0001, -1.3240, 0, 0, 0];      % 40 deg FREE
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('load40deg_v2.csv',2,48);
elseif choice == 3
    gama_rest = deg2rad(60);
    P = P_60;
    dL = dL_60;
    dphi = dphi_60;
    params.rcoeff = [r_rest, 0.0038/2];     % 60 degree FREE
    params.kelast = [-0.1282, -10.2866, 0, 0, 0; -0.0034, -0.5421, 0, 0, 0];      % 60 deg FREE
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('load60deg_v2.csv',2,50);
elseif choice == 4
    gama_rest = deg2rad(50);
    P = P_50;
    dL = dL_50;
    dphi = dphi_50;
    params.rcoeff = [r_rest, 0.0038/2];     % 50 degree FREE
    params.kelast = [0.0146146, 24.9761, 0, 0, 0; 0.0072126, -0.73525, 0, 0, 0];      % 50 deg FREE
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('load50deg_v2.csv',2,48);
elseif choice == 5
    gama_rest = deg2rad(40);
    P = P_40v2;
    dL = dL_40v2;
    dphi = dphi_40v2;
    params.rcoeff = [r_rest, 0.0038/2];     % 40 degree FREE
    params.kelast = [0.3142976, 117.4588, 0, 0, 0; 0.011105, -1.3208, 0, 0, 0];      % 40 deg FREE v2 (with sandpaper)
    [loadbolts,loadNm,loadP,loadphi] = importfile_loads('load40deg_v3.csv',2,44);
end


% % rcoeff: coefficients of relationship between r and P. r = c1 + c2*P.
% params.rcoeff = [r_rest, 0.0118/2];     % 20 degree FREE
% params.rcoeff = [r_rest, 0.0038/2];     % 40 degree FREE
% params.rcoeff = [r_rest, 0.0038/2];     % 60 degree FREE
% 
% % kelast: coefficients of elastomer functions [1, dL_norm, dl_norm^2,...; 1, dphi_norm, dphi_norm^2,...]
% params.kelast = [-0.0010e3, -7.0179e3, 0, 0, 0; -0.0217, -7.5372, 0, 0, 0];      % 20 deg FREE
% params.kelast = [0.9120, 147.5603, 0, 0, 0; 0.0001, -1.3240, 0, 0, 0];      % 40 deg FREE
% params.kelast = [-0.1282, -10.2866, 0, 0, 0; -0.0034, -0.5421, 0, 0, 0];      % 60 deg FREE

% params.kelast = [0 0 0 0 0; 0 0 0 0 0];      % no elastomer

%% Measured data at each fiber angle, with loads. Separate data by load

index1 = 0;
index2 = 0;

for index = 1:length(loadbolts)
    if loadbolts(index) == 1
        index1 = index1 + 1;
        Pvphi_1bolt(index1,:) = [loadP(index), loadphi(index)];
    elseif loadbolts(index) == 2
        index2 = index2 + 1;
        Pvphi_2bolt(index2,:) = [loadP(index), loadphi(index)];
    end
end

%% Set Loads and other Parameter values

params.theta_rest = -L_rest*tan(gama_rest)/(r_rest);
params.x_rest = [P_rest, gama_rest, r_rest, L_rest, phi_rest, 0];
params.xmin = [-pi/2, r_rest, 1, -Inf];
params.xmax = [pi/2, 0.75, 10, Inf];
params.Pmin = 3;
params.Pmax = 15;        % maximum pressure tested

% % Set the values of the external loads
% params.load = [0, 0];        % [F_load, M_load];

x(1,:) = params.x_rest + [0, 0, 1e-6, 0, 0, 0];    % initial simulation point (slightly more than zero so solver doesn't complain)
x1bolt(1,:) = params.x_rest + [0, 0, 1e-6, 0, 0, 0];    % initial simulation point (slightly more than zero so solver doesn't complain)
x2bolt(1,:) = params.x_rest + [0, 0, 1e-6, 0, 0, 0];    % initial simulation point (slightly more than zero so solver doesn't complain)

% Number of points at which to evaluate the steady state behavior
res = 50;  


%% Simulates steady state behavior of FREE for pressures from 0-15 psi
% Mload = Fload = 0
params.load = [0, 0];        % [F_load, M_load];
for i = 1:res
    P_des(i+1) = params.Pmin + (params.Pmax/res) * i;
    x0 = x(i,:);    %set initial guess of solution as state at last P_des
    [x(i+1,1), x(i+1,2), x(i+1,3), x(i+1,4), x(i+1,5), x(i+1,6)]  = solve_FBhoop(P_des(i+1), params.x_rest, params.load);
end

% Mload = 1 bolt, Fload = 0
params.load = [0, -(0.117933525 + 0.0)];        % [F_load, M_load];
% params.load = [0, -(0.117933525 + 0.1)];        % [F_load, M_load];
for j = 1:res
    P_des(j+1) = params.Pmin + (params.Pmax/res) * j;
    x0 = x1bolt(j,:);    %set initial guess of solution as state at last P_des
    [x1bolt(j+1,1), x1bolt(j+1,2), x1bolt(j+1,3), x1bolt(j+1,4), x1bolt(j+1,5), x1bolt(j+1,6)] = solve_FBhoop(P_des(j+1), params.x_rest, params.load);
end

% Mload = 2 bolts, Fload = 0
params.load = [0, -(0.23586705 + 0.0)];        % [F_load, M_load];
% params.load = [0, -(0.23586705 + 0.1)];        % [F_load, M_load];
for k = 1:res
    P_des(k+1) = params.Pmin + (params.Pmax/res) * k;
    x0 = x2bolt(k,:);    %set initial guess of solution as state at last P_des
    [x2bolt(k+1,1), x2bolt(k+1,2), x2bolt(k+1,3), x2bolt(k+1,4), x2bolt(k+1,5), x2bolt(k+1,6)] = solve_FBhoop(P_des(k+1), params.x_rest, params.load);
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

%% Generates Rotation vs. Pressure Plot Suitable for publication

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

%% Caluclate average relative error for each loading case
% fudge = 0.01;   % tolerance for equality comparisons
% 
% % No load error
% for j = 1:length(P)
%     for k = 1:length(x(:,1))
%         if P(j) < x(k,1)+fudge && P(j) > x(k,1)-fudge && abs(dphi(j)) > fudge
%             error_noload(j) = abs(rad2deg(dphi(j) - x(k,5)))/abs(rad2deg(dphi(j)));
%         end
%     end
% end  
% relerror_noload = sum(error_noload) / length(P);
% 
% 
% % 1 bolt load error
% for j = 1:length(Pvphi_1bolt(:,1))
%     for k = 1:length(x1bolt(:,1))
%         if Pvphi_1bolt(j,1) < x1bolt(k,1)+fudge && Pvphi_1bolt(j,1) > x1bolt(k,1)-fudge
%             error_1bolt(j) = abs(rad2deg(Pvphi_1bolt(j,2) - x1bolt(k,5)))/abs(rad2deg(Pvphi_1bolt(j,2)));
%         end
%     end
% end  
% relerror_1bolt = sum(error_1bolt) / length(Pvphi_1bolt(:,1));
% 
% 
% % 2 bolt load error
% for j = 1:length(Pvphi_2bolt(:,1))
%     for k = 1:length(x2bolt(:,1))
%         if Pvphi_2bolt(j,1) < x2bolt(k,1)+fudge && Pvphi_2bolt(j,1) > x2bolt(k,1)-fudge
%             error_2bolt(j) = abs(rad2deg(Pvphi_2bolt(j,2) - x2bolt(k,5)))/abs(rad2deg(Pvphi_2bolt(j,2)));
%         end
%     end
% end  
% relerror_2bolt = sum(error_2bolt) / length(Pvphi_2bolt(:,1));
% 
% % Put all errors into a matrix so it's easy to read later
% relerror = [0, relerror_noload;...
%          1, relerror_1bolt;...
%          2, relerror_2bolt];