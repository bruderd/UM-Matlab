function setJacobians_v2(params)
%setJacobians: Defines the Jacobians of the system symbolically then
%creates Matlab functions to evaluate them. Here the state, x, is just
%Euler angles
%   Detailed explanation goes here


%% Module Block Jacobian

syms a_ki b_ki L_k psi_k theta_k phi_k real

% The state of the ith actuator of the kth module (q_ki) in terms of the state of
% the kth module in local coordinate frame (X_k).
s_ki = -L_k + sqrt( (L_k - sin(theta_k)*a_ki + cos(theta_k)*sin(psi_k)*b_ki)^2 + (a_ki^2 + b_ki^2)*phi_k^2 );
w_ki = phi_k;

q_ki = [s_ki, w_ki]';
X_k = [psi_k, theta_k, phi_k]';

% define the module Jacobian
J_ki = jacobian(q_ki, X_k);

% Create Matlab function for evaluating J_ki
% matlabFunction(J_ki, 'File', 'J_ki', 'Vars', {[x_k, y_k, z_k, psi_k, theta_k, phi_k], [a_ki, b_ki], L_k});
matlabFunction(J_ki, 'File', 'J_ki', 'Vars', {X_k, [a_ki, b_ki], L_k});



%% Manipulator Block Jacobian (gets hung up when tries to take inverse)
% 
% % p = params.p;   % total number of modules in the manipulator
% p = 2;
% 
% % selection matrix to isolate position component of module state
% Spos = [eye(3), zeros(3,3); zeros(3,3), zeros(3,3)];
% Seul = [zeros(3,3), zeros(3,3); zeros(3,3), eye(3)];
% 
% % state vectors (in local and global coordinates)
% x = sym('x', [6*p,1], 'real');      % local coordinates
% x0 = sym('x0', [6*p,1], 'real');    % global coordinates
% T = zeros(6,6,p);       % Rotation transformation matrix from local to global
% T = sym(T);
% 
% % define x0 in terms of x
% x0(1:6,1) = x(1:6,1);
% 
% R0 = zeros(3*p,3);
% R0 = sym(R0);
% R0(1:3, 1:3) = Reuler(x(4:6, 1));
% for i = 2:p
%     xi = x(1+6*(i-1) : 6*i, 1);     % local coordinates of ith module
%     x0im1 = x0(1+6*(i-2) : 6*(i-1), 1);     % global coordinates of (i-1)th module
%     
%     Rim1_i = Reuler(xi(4:6));   % rotation matrix from ith to (i-1)th frame
%     R0_im1 = R0(1+3*(i-2) : 3*(i-1), 1:3); % rotation matrix from (i-1)th frame to 0 (global) frame
%     R0_i = R0_im1 * Rim1_i; % rotation matrix from ith frame to 0 (global) frame
%     
%     % perform coordinate transformation of xi to x0i
%     x0i(1:3, 1) = x0im1(1:3,1) + R0_im1 * xi(1:3, 1);
%     x0i(4:6, 1) = rot2euler_sym(R0_i);
%     
%     x0(1+6*(i-1) : 6*i, 1) = x0i;       % concatenate x0's of each module
%     
%     R0(1+3*(i-1) : 3*i, 1:3) = R0_i;    % storing the local to global rotations in the R0 matrix
% 
% end
% 
% % old version of for loop (doesn't work)
% % T(:,:,1) = eye(6);
% % 
% % for i = 2:p
% %     xi = x(1+6*(i-1) : 6*i, 1);
% %     Roti = Reuler(xi(4:6,1));
% %     T(:,:,i) = T(:,:, i-1) * [Roti, zeros(3,3); zeros(3,3), Roti];
% %     
% %     x0(1+6*(i-1) : 6*i, 1) = Spos * x0(1+6*(i-2) : 6*(i-1), 1) + T(:,:,i) * xi;
% % end
% 
% % define the manipulator Jacobian
% Jx = jacobian(x0, x);
% Jxinv = pinv(Jx);
% 
% % Create Matlab function for evaluating Jx
% matlabFunction(Jx, 'File', 'Jx', 'Vars', {x});
% matlabFunction(Jxinv, 'File', 'Jxinv', 'Vars', {x});


%% New version of Manipulator Jacobian that is a function of x0 not x

p = params.p;   % total number of modules in manipulator

% state vectors (in local and global coordinates)
x = sym('x', [3*p,1], 'real');      % local coordinates
x0 = sym('x0', [3*p,1], 'real');    % global coordinates

x = x0_orient2x_orient_sym(x0, params);   % define x in terms of x0

% define the manipulator Jacobian
Jx = jacobian(x, x0);

% Create Matlab function for evaluating Jx
matlabFunction(Jx, 'File', 'Jx', 'Vars', {x0}, 'Optimize', false);

%% Another version of Manipulator Jacobian that is a function of x not x0

% p = params.p;   % total number of modules in manipulator
% 
% % state vectors (in local and global coordinates)
% x = sym('x', [3*p,1], 'real');      % local coordinates
% x0 = sym('x0', [3*p,1], 'real');    % global coordinates
% 
% x0 = x_orient2x0_orient_sym(x, params);   % define x in terms of x0
% 
% % define the manipulator Jacobian
% Jx0 = jacobian(x0, x);
% 
% % Create Matlab function for evaluating Jx
% matlabFunction(Jx0, 'File', 'Jx0', 'Vars', {x}, 'Optimize', false);


%% Jacobian to transfer between orientation and position velocities

syms psi_ theta_ phi_ x y z L_ real

x_orient = [psi_, theta_, phi_]';
x_pos = euler2cart(x_orient, L_);

% define the Jacobian
Jorient2pos =  jacobian(x_pos, x_orient);

% Create Matlab function for evaluating Jx
matlabFunction(Jorient2pos, 'File', 'Jorient2pos', 'Vars', {x_orient, L_}, 'Optimize', false);

end

