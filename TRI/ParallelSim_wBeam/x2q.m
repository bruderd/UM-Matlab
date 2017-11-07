function q = x2q( x, params )
%x2q: Converts module state to actuator states
%   Detailed explanation goes here

p = params.p;   % number of modules
n = params.n;   % number of actuators in each module (a vector)

s_ki = -L_k + sqrt( (L_k - sin(theta_k)*a_ki + cos(theta_k)*sin(psi_k)*b_ki)^2 + (a_ki^2 + b_ki^2)*phi_k^2 );
w_ki = phi_k;


end

