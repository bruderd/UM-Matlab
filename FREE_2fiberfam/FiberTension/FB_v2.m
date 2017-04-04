% FB_v2.m
%   Force balance equations

function FB = FB_v2(x, u, x_prev, x_rest)

P = u;
[T_gama, T_betta] = deal(x(7), x(8));
[P, gama, betta, r, L, phi] = deal(x(1), x(2), x(3), x(4), x(5), x(6));
[P_prev, gama_prev, betta_prev, r_prev, L_prev, phi_prev] = deal(x_prev(1), x_prev(2), x_prev(3), x_prev(4), x_prev(5), x_prev(6));
[P0, gama0, betta0, r0, L0, phi0] = deal(x_rest(1), x_rest(2), x_rest(3), x_rest(4), x_rest(5), x_rest(6));

theta_gama0 = -tan(gama0)*L0/r0;     
theta_betta0 = -tan(betta0)*L0/r0;   
theta_gama = -tan(gama)*L/r; 
theta_betta = -tan(betta)*L/r; 

F_elast = 1*(L0-L);
M_elast = 0.08*(-1)*phi;
% F_elast = 0;
% M_elast = 0;

inputeq = u - P;
force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)           
gama_inextensible = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
betta_inextensible = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
twist_gama = (theta_gama - theta_gama0) - phi; 
twist_betta = (theta_betta - theta_betta0) - phi;

eta = sqrt(L0^2 + (r0*theta_gama0)^2)/sqrt(L0^2 + (r0*theta_betta0)^2);
fibers_coupled = gama - acos(cos(betta)/eta);


% Let's try some new things with the extra constraints
% %version 1
% extra_constraint1 = 2*pi*P*r^2 - (tan(abs(gama))*T_gama*sin(abs(gama)) + tan(abs(betta))*T_betta*sin(abs(betta)));  
% extra_constraint2 = 2*pi*P*r^2 - (tan(abs(-gama))*T_gama*sin(abs(-gama)) + tan(abs(-betta))*T_betta*sin(abs(-betta))); 

%version 2
ngama = floor(L*tan(abs(gama))/(r*pi));
nbetta = floor(L*tan(abs(betta))/(r*pi));
psigama = L*tan(abs(gama))/r - ngama*pi;
psibetta = L*tan(abs(betta))/r - nbetta*pi;
extra_constraint1 = 4*P*r*L - ( T_gama*sin(abs(gama)) * (2*ngama + 1 + cos(psigama)) + T_betta*sin(abs(betta)) * (2*nbetta + 1 + cos(psibetta)) );  
extra_constraint2 = 4*P*r*L - ( T_gama*sin(abs(-gama)) * (2*ngama + 1 + cos(psigama)) + T_betta*sin(abs(-betta)) * (2*nbetta + 1 + cos(psibetta)) );

% virtual work constraint PdV = FdL + Mdphi (DOESN'T WORK)
virtualwork = P*pi*(L*r^2 - L_prev*r_prev^2) - (force_balance - F_elast)*(L - L_prev) - (torque_balance - M_elast)*(phi - phi_prev);

% Volumetric transduction constraints
% Decide whether to used dV/dL or -dV/dL expression
if sign(gama_prev) == sign(betta_prev)
    alpha = betta;
    if abs(gama_prev) < deg2rad(54.7)
        vl = -1;
    else
        vl = 1;
    end
else
    if gama_prev^2 + betta_prev^2 < (pi/2)^2
        vl = -1;
    else
        vl = 1;
    end
end

% Decide whether to used dV/dphi or -dV/dphi expression
if gama_prev >= 0 && betta_prev >= 0      % quadrant I triangle      
    vt = 1;
elseif gama_prev < 0 && betta_prev < 0      % quadrant III triangle    
    vt = -1;
elseif gama_prev >= 0 && betta_prev < 0
    if gama_prev^2 + betta_prev^2 < (pi/2)^2    % green pie slice
        vt = 1;
        alpha = betta;
    elseif gama_prev^2 + betta_prev^2 >= (pi/2)^2    %orange sliver
        vt = -1;
        alpha = gama;
    end
elseif gama_prev < 0 && betta_prev >= 0
    if gama_prev^2 + betta_prev^2 < (pi/2)^2    % red pie slice
        vt = -1;
        alpha = gama;
    elseif gama_prev^2 + betta_prev^2 >= (pi/2)^2    % blue sliver
        vt = 1;
        alpha = betta;
    end
end

% Volumetric transduction equations
vol_length = P*pi*(L*r^2 - L_prev*r_prev^2) - vl * pi*r0^2*(1 - 2*cot(gama)^2)*(L - L_prev);
vol_twist = P*pi*(L*r^2 - L_prev*r_prev^2) + vt * 2*pi*r0^3*cot(alpha)*(phi - phi_prev);

FB = [force_balance;...
      torque_balance;...
%      virtualwork;...
      inputeq;...
      gama_inextensible;...
      betta_inextensible;...
      twist_gama;...
      twist_betta;...
%      fibers_coupled;...
      extra_constraint1;...
%       extra_constraint2;...
      vol_length;...
      vol_twist;...
      ];



end