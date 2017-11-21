function [tau, xdot] = manipulatorBlock(gama, alphadot, alpha , params)
%manipulatorBlock: Relates module torques to manipulator torques in global
%frame
%   Detailed explanation goes here

Jax = J_ax(alpha);   % define coordinate transformation Jacobian to go from global to local coordinates

%% set output values
xdot = Jax * alphadot;
tau = Jax' * gama;

end