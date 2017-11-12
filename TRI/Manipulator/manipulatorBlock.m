function [zeta0, xdot] = manipulatorBlock( zeta, xdot0, x0 , params)
%manipulatorBlock: Relates module torques to manipulator torques in global
%frame
%   Detailed explanation goes here

J_x = Jx(x0);   % invert coordinate transformation Jacobian to go from global to local coordinates

%% set output values
xdot = J_x * xdot0;
zeta0 = J_x' * zeta;

end