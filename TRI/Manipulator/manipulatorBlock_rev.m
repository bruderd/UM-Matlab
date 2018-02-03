function [zeta, xdot0] = manipulatorBlock_rev( zeta0, xdot, x , params)
%manipulatorBlock: Relates module torques to manipulator torques in global
%frame
%   Detailed explanation goes here

J_x0 = Jx0(x);   % define coordinate transformation Jacobian to go from global to local coordinates

%% set output values
xdot0 = J_x0 * xdot;
zeta = J_x0' * zeta0;

end