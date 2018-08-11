function xdot = myode_2d( t, x, u )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

xdot(1) = x(2);
xdot(2) = x(1) - x(1)^3 - 0.2*x(2) + 0.2*x(1)^2*u;

end

