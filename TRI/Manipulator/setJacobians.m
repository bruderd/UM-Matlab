function setJacobians(params)
%setJacobians: Defines the Jacobians of the system symbolically then
%creates Matlab functions to evaluate them.
%   Detailed explanation goes here


%% Module Block Jacobian

syms a_ki b_ki L_k x_k y_k z_k psi_k theta_k phi_k real

% The state of the ith actuator of the kth module (q_ki) in terms of the state of
% the kth module in local coordinate frame (X_k).
s_ki = -L_k + sqrt( (L_k - sin(theta_k)*a_ki + cos(theta_k)*sin(psi_k)*b_ki)^2 + (a_ki^2 + b_ki^2)*phi_k^2 );
w_ki = phi_k;

q_ki = [s_ki, w_ki]';
X_k = [x_k, y_k, z_k, psi_k, theta_k, phi_k];

% define the module Jacobian
J_ki = jacobian(q_ki, X_k);

% Create Matlab function for evaluating J_ki
% matlabFunction(J_ki, 'File', 'J_ki', 'Vars', {[x_k, y_k, z_k, psi_k, theta_k, phi_k], [a_ki, b_ki], L_k});
matlabFunction(J_ki, 'File', 'J_ki', 'Vars', {X_k, [a_ki, b_ki], L_k});



%% Manipulator Block Jacobian

% p = params.p;   % total number of modules in the manipulator
p = 2;

% selection matrix to isolate position component of module state
Spos = [eye(3), zeros(3,3); zeros(3,3), zeros(3,3)];
Seul = [zeros(3,3), zeros(3,3); zeros(3,3), eye(3)];

% state vectors (in local and global coordinates)
x = sym('x', [6*p,1]);      % local coordinates
x0 = sym('x0', [6*p,1]);    % global coordinates
T = zeros(6,6,p);       % Rotation transformation matrix from local to global
T = sym(T);

% define x0 in terms of x
x0(1:6,1) = x(1:6,1);
T(:,:,1) = eye(6);
for i = 2:p
    xi = x(1+6*(i-1) : 6*i, 1);
    Roti = Reuler(xi(4:6,1));
    T(:,:,i) = T(:,:, i-1) * [Roti, zeros(3,3); zeros(3,3), Roti];
    
    x0(1+6*(i-1) : 6*i, 1) = Spos * x0(1+6*(i-2) : 6*(i-1), 1) + T(:,:,i) * xi;
end

% define the manipulator Jacobian
Jx = jacobian(x0, x);
Jxinv = inv(Jx);

% Create Matlab function for evaluating Jx
matlabFunction(Jx, 'File', 'Jx', 'Vars', {x});
matlabFunction(Jxinv, 'File', 'Jxinv', 'Vars', {x});

end

