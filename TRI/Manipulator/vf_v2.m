function [f, dfdX, dfdu, dfdXdot] = vf_v2(X0_orient, u, X0dot_orient, params)
%vf_v1: Dynamics soft robotic manipulator
%   Detailed explanation goes here

%% Define local names of global parameters

p = params.p;       % total number of modules in manipulator
m = params.m;       % masses of the blocks
I = params.I;       % moment of inertia matrices of blocks
L = params.L;       % lengths of the modules

%% Caclulate the manipulator forces (zeta0)

% % isolate just the euler angles in the input vectors
% X0_orient = kron(eye(p), [zeros(3,3), eye(3)]) * X0;
% X0dot_orient = kron(eye(p), [zeros(3,3), eye(3)]) * X0dot;

% x0 = X0(1:6*p);
% for i = 1:p
%     x0i_orient = X0_orient(3*(i-1)+1 : 3*i);
%     x0(6*(i-1)+1, 6*(i-1)+3) = euler2cart
%     x0(6*(i-1)+4, 6*i) = x0i_orient;
% end

% x0_orient = X0_orient(1:3*p, 1);
% x0dot_orient = X0dot_orient(3*p+1 : 2*3*p, 1);
% 
% % derive x0dot that includes position components from just the euler angles
% x0dot = zeros(6*p, 1);
% for i = 1:p
%     x0i_orient = X0_orient(3*(i-1)+1 : 3*i);
%     x0doti_orient = X0dot_orient(3*(i-1)+1 : 3*i);
%     x0dot(6*(i-1)+1 : 6*(i-1)+3, 1) = Jorient2pos(x0i_orient, L(i)) * x0doti_orient;
%     x0dot(6*(i-1)+4 : 6*i, 1) = x0doti_orient;
% end
% 
% % x = x02x(x0, params);   % convert x0 to local coordinates x
% x = x0_orient2x(x0_orient, params);      % derive x from just x0_orient, i.e. euler angles in global frame
% x0 = x0_orient2x0(x0_orient, params);      % derive x0 from just x0_orient, i.e. euler angles in global frame

P = u;      % input is vector of pressures
x0 = X0_orient(1 : 3*p, 1);
x0dot = X0_orient(3*p+1 : 2*3*p, 1);
x = x0_orient2x_orient(x0, params);

[~, xdot] = manipulatorBlock(0, x0dot, x0, params);
[~, qdot] = moduleBlock(0, xdot, x, params);
[Z, Vdot] = actuatorBlock(P, qdot, x, params);
[zeta, qdot] = moduleBlock(Z, xdot, x, params);
[zeta0, xdot] = manipulatorBlock(zeta, x0dot, x0, params);

% zeta0_orient = kron(eye(p), [zeros(3,3), eye(3)]) * zeta0;  % isolate just the torques


%% Equations of Motion, f(xddot, xdot, x, u) = 0
f(1:3*p, 1) = X0dot_orient(1:3*p) - X0_orient(3*p+1 : 2*(3*p));
f(3*p+1 : 2*(3*p), 1) = EOM(X0_orient, X0dot_orient, zeta0, m, I) + X0dot_orient(1:3*p, 1)*(0.01);     % damping term haphazardly added

%% Gradients (not needed for now)
dfdX = [0];
dfdu = [0];
dfdXdot = [0];

end

