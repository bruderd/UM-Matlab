% Uses polynomial fits of P vs L and P vs theta to determine the spring
% constant functions of the system.
clear
clc

%% Measured values

P_20 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 10, 8.5, 7.5, 5, 6.5, 4, 2.5, 6.5];
dL_20 = P_20*(-1e-3);
dphi_20 = [0.1419999879,...
           0.2387610417, 0.3443185548, 0.469982261, 0.6283185307, 0.7879114375, 1.03044239, 1.043008761, 0.8168140899,...
           0.6710441908, 0.3933274002, 0.5152211952, 0.266407057, 0.1445132621, 0.4611858015];       

% 40 deg FREE
P_40 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 11, 8.5, 6.5, 4, 2.5];
dL_40 = [0, 0, 0, 0, 0, 0.004, 0.02, 0.028, 0.04, 0.06, 0.012, 0.012, 0.012, 0.012, 0.012, 0.016, 0.02, 0.028, 0.044, 0.06, 0.06, 0.06, 0.048, 0.036, 0.032];
% dphi_40 = [1.301875996, 2.045805136, 2.724389149, 3.624141285, 4.518866873, 5.292955303, 6.690335715, 7.816282522, 8.992494812, 10.24410532, 1.447645895,...
%            2.37755732, 3.086300623, 4.021238597, 4.926017281, 5.700105711, 6.976848965, 7.856494908, 9.228742579, 10.50548583, 8.384282474, 6.172601246,...
%            4.337911136, 2.714336053, 1.512991022, 0, 0, 0];
dphi_40 = [0.3254689989, 0.511451284, 0.6810972873, 0.9060353213, 1.129716718, 1.323238826, 1.672583929, 1.954070631,...
           2.248123703, 2.561026331, 0.3619114737, 0.5943893301, 0.7715751557, 1.005309649, 1.23150432, 1.425026428,...
           1.744212241, 1.964123727, 2.307185645, 2.626371458, 2.096070618, 1.543150311, 1.084477784, 0.6785840132, 0.3782477555];

       
% 60 deg FREE
P_60 = [2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 2.5, 4, 5, 6.5, 7.5, 8.5, 10, 11, 12, 13, 10, 12, 8.5, 5, 2.5];
dL_60 = [0, 0, 0.012, 0.032, 0.052, 0.064, 0.084, 0.108, 0.132, 0.164, 0, 0, 0.012, 0.024, 0.04, 0.056, 0.084, 0.108, 0.136, 0.144, 0.136, 0.144, 0.104, 0.064, 0.024];
% dphi_60 = [1.26166361, 2.307185645, 3.201911233, 4.423362456, 5.227610176, 6.021804798, 7.404105566, 8.33904354, 9.927432785, 11.46053, 1.613521987, 2.437875899,... 
%            3.518583772, 4.322831491, 5.162265048, 6.308318048, 7.554902013, 8.891963847, 10.47532654, 11.37507868, 8.675822272, 10.5959637, 6.936636579, 4.09663682, 2.010619298, 0, 0, 0];
dphi_60 = [0.3154159024, 0.5767964112, 0.8004778081, 1.105840614, 1.306902544, 1.5054512, 1.851026391, 2.084760885,...
           2.481858196, 2.8651325, 0.4033804967, 0.6094689748, 0.879645943, 1.080707873, 1.290566262, 1.577079512,...
           1.888725503, 2.222990962, 2.618831636, 2.84376967, 2.168955568, 2.648990926, 1.734159145, 1.024159205, 0.5026548246];


% 50 deg FREE
P_50 = [2.5, 4, 5, 6, 7.5, 8.5, 10, 11, 12.5, 13, 14.5, 12.5, 13.5, 11, 12, 10, 11, 8.5, 10, 7.5, 9, 6.5, 8, 5.4, 6.5,...
        4, 5, 3.4, 0];
dL_50 = [0, 0, 0.004, 0.008, 0.024, 0.048, 0.06, 0.084, 0.104, 0.148, 0.18, 0.148, 0.164, 0.14, 0.148, 0.116,...
         0.124, 0.096, 0.1, 0.096, 0.096, 0.076, 0.076, 0.064, 0.064, 0.056, 0.056, 0.052, 0.032];
dphi_50 = [0.3506017401, 0.5981592412, 0.7841415263, 0.981433545, 1.303132633, 1.559486593, 1.901291874, 2.280796267,...
           2.816123655, 3.204424507, 3.751061628, 3.183061677, 3.507274038, 2.790990913, 3.204424507, 2.454212181,...
           2.729415697, 2.038265314, 2.408973247, 1.807044094, 2.065911329, 1.501681288, 1.771858257, 1.221451224,...
           1.476548547, 0.9437344331, 1.109610525, 0.7690618816, 0.1445132621];


% 40 deg FREE with sandpaper
P_40v2 = [3, 4.3, 5.4, 6.3, 8, 8.8, 10.1, 11.1, 12.2, 13.7, 11, 12.3, 10, 11.2, 8.9, 10.1, 7.8, 9, 7,...
          8, 5.2, 6.6, 4, 5, 3, 0];
dL_40v2 = [0, 0.004, 0.008, 0.012, 0.024, 0.036, 0.052, 0.064, 0.084, 0.112, 0.096, 0.104, 0.084, 0.084, 0.08, 0.08, 0.064, 0.064, 0.056, 0.056, 0.048,...
           0.048, 0.044, 0.044, 0.036, 0.024];
dphi_40v2 = [0.3581415625, 0.5504070329, 0.736389318, 0.9022654101, 1.198831757, 1.394867138, 1.751752064, 2.00559275, 2.333575023, 2.788477639,...
             2.141309553, 2.469291826, 1.922654704, 2.173982116, 1.632371543, 1.928937889, 1.394867138, 1.642424639, 1.21893795, 1.422513154,...
             0.8984954989, 1.112123799, 0.6723008279, 0.8432034682, 0.5215043805, 0.08168140899];
         
% 40 deg FREE with sandpaper
P_40v2 = [3, 4.3, 5.4, 6.3, 8, 8.8, 10.1, 11.1, 12.2, 13.7, 11, 12.3, 10, 11.2, 8.9, 10.1, 7.8, 9, 7,...
          8, 5.2, 6.6, 4, 5, 3];
dL_40v2 = [0, 0.004, 0.008, 0.012, 0.024, 0.036, 0.052, 0.064, 0.084, 0.112, 0.096, 0.104, 0.084, 0.084, 0.08, 0.08, 0.064, 0.064, 0.056, 0.056, 0.048,...
           0.048, 0.044, 0.044, 0.036];
dphi_40v2 = [0.3581415625, 0.5504070329, 0.736389318, 0.9022654101, 1.198831757, 1.394867138, 1.751752064, 2.00559275, 2.333575023, 2.788477639,...
             2.141309553, 2.469291826, 1.922654704, 2.173982116, 1.632371543, 1.928937889, 1.394867138, 1.642424639, 1.21893795, 1.422513154,...
             0.8984954989, 1.112123799, 0.6723008279, 0.8432034682, 0.5215043805];
       
%% Choose which measured values to use: 1 = 20 deg FREE, 2 = 40 deg FREE, 3 = 60 deg FREE

choice = 2; % MAKE YOUR CHOICE HERE

if choice == 1
    gama_rest = deg2rad(20);
    P = P_20;
    dL = dL_20;
    dphi = dphi_20;
    dr = 0.0118*P/2;    % for 20 deg FREE
elseif choice == 2
    gama_rest = deg2rad(40);
    P = P_40v2;
    dL = dL_40v2;
    dphi = dphi_40v2;
    dr = 0.0038*P/2;    % for 40 deg FREE
elseif choice == 3
    gama_rest = deg2rad(60);
    P = P_60;
    dL = dL_60;
    dphi = dphi_60;
    dr = 0.0038*P/2;    % for 60 deg FREE
elseif choice == 4
    gama_rest = deg2rad(50);
    P = P_50;
    dL = dL_50;
    dphi = dphi_50;
    dr = 0.0038*P/2;    % for 60 deg FREE
end

%% Define other initial parameters
P_rest = 0;   
r_rest = 0.1875;    
L_rest = 3.4;
phi_rest = 0;
theta_rest = -tan(gama_rest)*L_rest/r_rest;

%% Generate elastomer force/torques (Fs and Ms) at each pressure

% Generates r at each pressure in P from linear equation fit
r = dr + r_rest*ones(size(P)); 

% Generates L at each pressure in P
L = dL + L_rest*ones(size(dL));

% Generates theta at each pressure in P
theta = dphi + theta_rest*ones(size(dphi));

% Generates phi at each pressure in P
phi = dphi + 0;

% Generates gama at each pressure in P
gama = atan(r.*(-theta)./L);

% Normalized displacements (divided by L0)
dL_norm = dL/L_rest;
dphi_norm = (dphi.*r)/L_rest;




% Generates T, E, and G at each time step
t = 1/32 * ones(size(P)); % + P*(1/32)/15;   % thickness of tube (made up)
for k = 1:length(P)
    
%     % with tube thickness inclued (assumed constant)
%     E(k) = (L_rest*P(k)*r(k)*r_rest*(sin(gama(k)) - 2*cos(gama(k))*cot(gama(k))))/(2*t(k)*(L(k)*r_rest*sin(gama(k)) - L_rest*r_rest*sin(gama(k)) - L_rest*r(k)*cos(gama(k))*cot(gama(k)) + L_rest*r_rest*cos(gama(k))*cot(gama(k))));
%     G(k) = -(L(k)*P(k)*cot(gama(k))*(L_rest*r(k)*sin(gama(k)) - 2*L(k)*r_rest*sin(gama(k)) + L_rest*r_rest*sin(gama(k))))/(2*phi(k)*t(k)*(L(k)*r_rest*sin(gama(k)) - L_rest*r_rest*sin(gama(k)) - L_rest*r(k)*cos(gama(k))*cot(gama(k)) + L_rest*r_rest*cos(gama(k))*cot(gama(k))));
    
   % without tube thickness (2D planar stresses)
    E(k) = (L_rest*P(k)*r(k)*r_rest*(sin(gama(k)) - 2*cos(gama(k))*cot(gama(k))))/(2*(L(k)*r_rest*sin(gama(k)) - L_rest*r_rest*sin(gama(k)) - L_rest*r(k)*cos(gama(k))*cot(gama(k)) + L_rest*r_rest*cos(gama(k))*cot(gama(k))));
    G(k) = -(L(k)*P(k)*cot(gama(k))*(L_rest*r(k)*sin(gama(k)) - 2*L(k)*r_rest*sin(gama(k)) + L_rest*r_rest*sin(gama(k))))/(2*phi(k)*(L(k)*r_rest*sin(gama(k)) - L_rest*r_rest*sin(gama(k)) - L_rest*r(k)*cos(gama(k))*cot(gama(k)) + L_rest*r_rest*cos(gama(k))*cot(gama(k))));
    
end
    


% % Generates tension in fiber, T, at each pressure in P
% T = 2.*P.*pi.*r.^2.*cot(gama)./sin(gama);
% 
% % Generates tension induced axial force at each pressure in P
% Ft = T.*cos(gama);
% 
% % Generates tension induced torque at each pressure in P
% Mt = r.*T.*sin(gama);
% 
% % Generates outward pressure force on endcaps at each pressure in P
% Fp = pi.*P.*r.^2;
% 
% % Generates the spring counter forces/torques at each pressure in P
% Fs = (Ft - Fp);    % Fp - Ft + Fs = 0
% Ms = -Mt;          % Mt + Ms = 0


%% Fits a polynomial to each spring FORCE/MOMENT as a functions of P and phi normalized
reg_E = MultiPolyRegress(dL_norm',E',1);
reg_G = MultiPolyRegress(dphi_norm',G',1);

E_poly = polyfit(dL_norm',E',1);
G_poly = polyfit(dphi_norm',G',1);

%% Generate Plots of interest
figure

% F_elast vs. P
subplot(2,2,1)
set(gca,'fontsize', 32);
%axis([5, 12*6, -100, 150]);
xlabel('$P$ (psi)', 'fontsize', 32, 'Interpreter', 'LaTex');
ylabel('$E$ (lb/in)', 'fontsize', 32, 'Interpreter', 'LaTex');
hold on
plot(P,E,'*')
hold off

% M_elast vs. P
subplot(2,2,2)
set(gca,'fontsize', 32);
%axis([5, 12*6, -100, 150]);
xlabel('$P$ (psi)', 'fontsize', 32, 'Interpreter', 'LaTex');
ylabel('$G$ (lb/in)', 'fontsize', 32, 'Interpreter', 'LaTex');
hold on
plot(P,G,'*')
hold off

% F_elast vs. (L-L0/L0)
subplot(2,2,3)
set(gca,'fontsize', 32);
xlabel('$\frac{l-L}{L}$ (\%)', 'fontsize', 32, 'Interpreter', 'LaTex');
ylabel('$F_{elast}$ (lb)', 'fontsize', 32, 'Interpreter', 'LaTex');
hold on
plot(dL_norm,E,'*')
plot(linspace(0,0.05,100), polyval(E_poly, linspace(0,0.05,100)))
hold off


% M_elast vs. (phi/L0)
subplot(2,2,4)
set(gca,'fontsize', 32);
xlabel('$\frac{\Phi r}{L}$ (rad/in)', 'fontsize', 32, 'Interpreter', 'LaTex');
ylabel('$M_{elast}$ (lb-in)', 'fontsize', 32, 'Interpreter', 'LaTex');
hold on
plot(dphi_norm,G,'*')
plot(linspace(0,0.2,100), polyval(G_poly, linspace(0,0.2,100)))
hold off




