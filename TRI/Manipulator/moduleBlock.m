function [zeta, qdot] = moduleBlock( Z, xdot, x, params)
%moduleBlock: Relates actuator torques to module torques.
%   Detailed explanation goes here

p = params.p;   % number of modules in the manipulator
n = params.n;   % number of actuators in each module (a vector)


%% for x state vector that includes only orientation
% Assemble the Jacobian out of Jacobians of each actuator
Jki = zeros(2, 3, max(n));     % Jacobian for all frees combined into single 3d matrix
J = zeros(2*sum(n), 3*p);
for i = 1:p
    xi = x(1+3*(i-1) : 3*i, 1);
      
    Jk = zeros(2*n(i), 3);  % ensure the module Jacobian is the right size.
    for j = 1:n(i)
       Jki(:,:,j) = J_ki(xi, params.attach(j,:), params.L(i));      % Jacobian for each actuator in ith module
       Jk(2*(j-1)+1 : 2*j, :) = Jki(:,:,j);   % Jacobian for the ith module. Stacks actuator Jacobians vertically.
    end
    
    if i == 1
        J(1 : 2*n(i), 1 : 3*i) = Jk(:,:);
    else
        J(2*(n(i-1))+1 : 2*(n(i-1))+2*n(i), 3*(i-1)+1 : 3*i) = Jk(:,:);      % Jacobian for the full manipulator. Stacks module Jacobians diagonally.
    end



%% for x state vector that includes position and orientation
% % Assemble the Jacobian out of Jacobians of each actuator
% Jki = zeros(2, 6, max(n));     % Jacobian for all frees combined into single 3d matrix
% J = zeros(2*sum(n), 6*p);
% for i = 1:p
%     xi = x(1+6*(i-1) : 6*i, 1);
%       
%     Jk = zeros(2*n(i), 6);  % ensure the module Jacobian is the right size.
%     for j = 1:n(i)
%        Jki(:,:,j) = J_ki(xi, params.attach(i,:), params.L(i));      % Jacobian for each actuator in ith module
%        Jk(2*(j-1)+1 : 2*j, :) = Jki(:,:,j);   % Jacobian for the ith module. Stacks actuator Jacobians vertically.
%     end
%     
%     if i == 1
%         J(1 : 2*n(i), 1 : 6*i) = Jk(:,:);
%     else
%         J(2*(n(i-1))+1 : 2*(n(i-1))+2*n(i), 6*(i-1)+1 : 6*i) = Jk(:,:);      % Jacobian for the full manipulator. Stacks module Jacobians diagonally.
%     end
% end

%% set output values
qdot = J * xdot;
zeta = J' * Z;

end

