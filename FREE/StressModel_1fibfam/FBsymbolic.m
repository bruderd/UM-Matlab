function [Fsymbolic, Jsymbolic] = FBsymbolic(params)

syms P gama r L phi T P0 gama0 r0 L0 phi0 T0 t0 Fload Mload u

x = [P, gama, r, L, phi, T];

theta_gama0 = -tan(gama0)*L0/r0;       
theta_gama = -tan(gama)*L/r; 

%% Tube wall thickness equation
t = -r + sqrt(r^2 + 2*r0*t0 +t0^2); 

%% Stress/strain equations
dL_norm = (L - L0)/L0;
dphi_norm = atan(r*phi)/L0;

E = params.modulus(1,1)*dL_norm + params.modulus(1,2);  % varying modulus
G = params.modulus(2,1)*dphi_norm + params.modulus(2,2);    % varying modulus

sig_z = E * (L - L0)/L0;
sig_theta = E * (r - r0)/r0;        % removed factor of 2*pi
tau_ztheta = G * atan((r+t/2)*phi/L);     % added atan since x ~= tan(x) only for small x

%% Force Balance Equations
inputeq = u - P;    % not sure if needed or can just set value of P directly
hoop_balance = 2*pi*P*r^2*cot(gama) - 2*T*sin(gama) - 2*sig_theta*(pi*r*cot(gama))*t;
force_balance = P*pi*r^2 - 2*T*cos(gama) - pi*(2*r*t + t^2)*sig_z + Fload;
torque_balance = 2*(r+t)*T*sin(gama) - pi*(2*r*t + t^2)*(r + t/2)*tau_ztheta + Mload;

%% Geometric Constraints
gama_inextensible = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
twist_gama = (theta_gama - theta_gama0) - phi; 


% %% Volumetric transduction constraints
% % Decide whether to used dV/dL or -dV/dL expression
% if sign(gama_prev) == sign(betta_prev)
%     alpha = betta;
%     if abs(gama_prev) < deg2rad(54.7)
%         vl = -1;
%     else
%         vl = 1;
%     end
% else
%     if gama_prev^2 + betta_prev^2 < (pi/2)^2
%         vl = -1;
%     else
%         vl = 1;
%     end
% end
% 
% % Decide whether to used dV/dphi or -dV/dphi expression
% if gama_prev >= 0 && betta_prev >= 0      % quadrant I triangle      
%     vt = 1;
% elseif gama_prev < 0 && betta_prev < 0      % quadrant III triangle    
%     vt = -1;
% elseif gama_prev >= 0 && betta_prev < 0
%     if gama_prev^2 + betta_prev^2 < (pi/2)^2    % green pie slice
%         vt = 1;
%         alpha = betta;
%     elseif gama_prev^2 + betta_prev^2 >= (pi/2)^2    %orange sliver
%         vt = -1;
%         alpha = gama;
%     end
% elseif gama_prev < 0 && betta_prev >= 0
%     if gama_prev^2 + betta_prev^2 < (pi/2)^2    % red pie slice
%         vt = -1;
%         alpha = gama;
%     elseif gama_prev^2 + betta_prev^2 >= (pi/2)^2    % blue sliver
%         vt = 1;
%         alpha = betta;
%     end
% end
% 
% % Volumetric transduction equations
% vol_length = P*pi*(L*r^2 - L_prev*r_prev^2) - vl * pi*r0^2*(1 - 2*cot(gama)^2)*(L - L_prev);
% vol_twist = P*pi*(L*r^2 - L_prev*r_prev^2) + vt * 2*pi*r0^3*cot(alpha)*(phi - phi_prev);
% 

%% Create system of equations
Fsym = [hoop_balance;...
      force_balance;...
      torque_balance;...
      inputeq;...
      gama_inextensible;...
      twist_gama;...
%      vol_length;...
%      vol_twist;...
      ];

%% The Jacobian of Fsym, or dFsym/dx where x = [P, gama, r, L, phi, T]
Jsym = jacobian(Fsym, x);


%% Create Matlab Functions to calculate Fsym and Jsym

Feval = matlabFunction(Fsym, 'File', 'Feval', 'Vars', {[P, gama, r, L, phi, T], [P0, gama0, r0, L0, phi0, T0], t0, [Fload, Mload], u});
Jeval = matlabFunction(Jsym, 'File', 'Jeval', 'Vars', {[P, gama, r, L, phi, T], [P0, gama0, r0, L0, phi0, T0], t0, [Fload, Mload], u});

end