function [ t,y ] = sim_SpringMassDamper( x0, tf, m, b, k )
%sim_SpringMassDamper: Simulate the behavior of a spring mass damper system
%of the form: mx'' + bx' + kx = 0
%   x0: initial state
%   tf: final time of simulation
%   m: mass
%   b: damping coefficent
%   k: spring constant


tspan = [0, tf];
[t, y] = ode45(@(t,x) dynamics(x, m, b, k), tspan, x0);


end

function xdot = dynamics(x, m, b, k)
    xdot = [0, 1; -k/m, -b/m] * x;
end