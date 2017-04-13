% FBhoop.m
%   Force balance equations that account for the effect of the elastomer
%   bearing hoop stress
%
%   NOTE: This is for a FREE with 1 fiber family

function FB = FBhoop(x, u, x_prev, x_rest, load)

[Fload, Mload] = deal(load(1), load(2));

T = x(6);
[P, gama, r, L, phi] = deal(x(1), x(2), x(3), x(4), x(5));
% [P_prev, gama_prev, r_prev, L_prev, phi_prev] = deal(x_prev(1), x_prev(2), x_prev(3), x_prev(4), x_prev(5));
[P0, gama0, r0, L0, phi0] = deal(x_rest(1), x_rest(2), x_rest(3), x_rest(4), x_rest(5));

theta_gama0 = -tan(gama0)*L0/r0;       
theta_gama = -tan(gama)*L/r; 

% F_elast = -1*(L0-L);
% M_elast = 0.08*(-1)*phi;

%% Stress/strain equations
% E = 13.16174.*1 + 142.2554.*((L - L0)/L0);  % Young's modulus of latex (find via sysid)
% G = 7.02151.*1 - 16.2889.*(r*phi/L);  % Shear modululs of latex (find via sysid)
E = 20;
G = 9;

sig_x = E * (L - L0)/L0;
sig_r = E * (r - r0)/r0;        % removed factor of 2*pi
tau_rx = G * atan(r*phi/L);     % added atan since x ~= tan(x) only for small x

%% Force Balance Equations
inputeq = u - P;    % not sure if needed or can just set value of P directly
hoop_balance = 2*pi*P*r^2*cot(gama) - 2*T*sin(gama) - 2*sig_r*(pi*r*cot(gama));
force_balance = P*pi*r^2 - 2*T*cos(gama) - 2*pi*r*sig_x + Fload;
torque_balance = 2*r*T*sin(gama) - 2*pi*r^2*tau_rx + Mload;

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
FB = [hoop_balance;...
      force_balance;...
      torque_balance;...
      inputeq;...
      gama_inextensible;...
      twist_gama;...
%      vol_length;...
%      vol_twist;...
      ];



end