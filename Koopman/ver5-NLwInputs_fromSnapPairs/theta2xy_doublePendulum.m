function [ x_xy ] = theta2xy_doublePendulum( x_theta, params )
%theta2xy_doublePendulum: Converts state in terms of theta to state in
%terms of xy coordinates
%   Detailed explanation goes here

theta1 = x_theta(1);
theta2 = x_theta(2);
dtheta1 = x_theta(3);
dtheta2 = x_theta(4);
x1 = params.l1 * sin(theta1);
y1 = -params.l1 * cos(theta1);
x2 = x1 + params.l2 * sin(theta2);
y2 = y1 - params.l2 * cos(theta2);
x1dot = params.l1 * cos(theta1) .* dtheta1;
y1dot = params.l1 * sin(theta1) .* dtheta1;
x2dot = x1dot + params.l2 * cos(theta2) .* dtheta2;
y2dot = y1dot + params.l2 * sin(theta2) .* dtheta2;

x_xy = [x1, y1, x2, y2, x1dot, y1dot, x2dot, y2dot]';    % full state in xy coordinates

end

