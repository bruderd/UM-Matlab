% FBhoop.m
%   Force balance equations that account for the effect of the elastomer
%   bearing hoop stress
%
%   x = [P, gama, r, L, phi, T]
%   u = P_desired
%   params includes:
%       params.x_rest = "states of FREE when P=0 and no loads applied"
%       params.load = [Fload, Mload]
%       params.t_rest = "thickness of tube when P=0 and no loads applied"
%
%   NOTE: This is for a FREE with 1 fiber family

function FB = FBstress(x, u, params)

[Fload, Mload] = deal(params.load(1), params.load(2));

[P, gama, r, L, phi, T] = deal(x(1), x(2), x(3), x(4), x(5), x(6));
% [P_prev, gama_prev, r_prev, L_prev, phi_prev] = deal(x_prev(1), x_prev(2), x_prev(3), x_prev(4), x_prev(5));
[P0, gama0, r0, L0, phi0, T0] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4), params.x_rest(5), params.x_rest(6));
t0 = params.t_rest;

theta_gama0 = -tan(gama0)*L0/r0;       
theta_gama = -tan(gama)*L/r; 


%% Stress/strain equations
E = params.modulus(1);
G = params.modulus(2);

sig_z = E * (L - L0)/L0;
sig_theta = E * (r - r0)/r0;        % removed factor of 2*pi
tau_ztheta = G * atan(r*phi/L);     % added atan since x ~= tan(x) only for small x

%% Tube wall thickness equation
t = -r + sqrt(r^2 + 2*r0*t0 +t0^2); 

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