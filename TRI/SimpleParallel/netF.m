function netF = netF( q, P, params )
% netF.m
%   Determine the net force and torque on the central element 

[L_l, B_l, N_l] = deal(params.L_l, params.B_l, params.N_l); 
[L_r, B_r, N_r] = deal(params.L_r, params.B_r, params.N_r);
K = [params.kelast_l(1) + params.kelast_r(1), 0; 0, params.kelast_l(2) + params.kelast_r(2)];
load = params.load;
[s, w] = deal(q(1), q(2));


% Jacobian varies based on sign of fiber angle due to the absolute value in the radius expression
J_l = [pi*(B_l^2 - 3*(L_l+s)^2) / (2*pi*N_l+w)^2,...
    2*pi*(L_l+s)*((L_l+s)^2 - B_l^2) / (2*pi*N_l+w)^3];
    
J_r = [-pi*(B_r^2 - 3*(L_r-s)^2) / (2*pi*N_r-w)^2,...
   -2*pi*(L_r-s)*((L_r-s)^2 - B_r^2) / (2*pi*N_r-w)^3];


% assemble full jacobian
J = [J_l; J_r];

netF = J'*P + K*q + load;


end

