function [zeta, qdot] = moduleBlock( Z, alphadot, alpha, params)
%moduleBlock: Relates actuator torques to module torques.
%   Detailed explanation goes here

p = params.p;   % number of modules in the manipulator
n = params.n;   % number of actuators in each module (a vector)
K = params.K;   % spine stiffness
D = params.D;   % spine damping

Jaq = J_aq(alpha);

%% set output values
qdot = Jaq * alphadot;
zeta = Jaq' * Z;% + K * alpha + D * alphadot;

end

