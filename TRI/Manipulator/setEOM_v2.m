function setEOM_v2(params)
%setDynamics: Symbolically derives Euler-Lagrange EOM, then creates
%function to evaluate them.
%   Assumptions: 
%       -all mass is concentrated at the module connecting blocks (i.e. FREEs and spine are massless).
%       -spine does not contribute forces.

% suppress inconsequential warning
warn_id = 'symbolic:sym:sym:DeprecateExpressions';
warning('off', warn_id);

%% Define symbolic parameters

p = params.p;   % total number of modules
n = params.n;   % number of actuators in each module (a vector)
L = params.L;   % length of each module

g = 9.81;   % acceleration due to gravity, 9.81 m/s^s

m = sym('m', [p,1], 'real');    % masses of the module blocks
I = sym('I', [3*p,3], 'real');      % rotational moment of inertia matrices of the module blocks expressed in the local coordinate frame, vertically concatenated


x0_orient = sym('x0_orient', [3*p,1], 'real');
% x0 = sym('x0', [6*p,1], 'real');

% for i = 1:p
%     x0i_orient = x0_orient(1+3*(i-1) : 3*i, 1);
%     x0i(1:3, 1) = euler2cart(x0i_orient, L(i));
%     x0i(4:6, 1) = x0i_orient;
%     x0(1+6*(i-1) : 6*i, 1) = x0i;
% end
x0 = x0_orient2x0_sym(x0_orient, params);

dx0dx0_orient = jacobian(x0, x0_orient);   % deriviative of the x0 vector with respect to the last 3 components of each module, which describe orientation.
% dx0_orientdx0 = kron(eye(p), [zeros(3,3), eye(3)]);     % deriviative of the x0_orientation vector with respect to x0. Basically eliminates everything except
% dx0_orientdx0 = pinv(dx0dx0_orient);

x0dot_orient = sym('x0dot_orient', [3*p, 1], 'real');
% x0dot = sym('x0dot', [6*p,1], 'real');
x0dot = dx0dx0_orient * x0dot_orient;   % chain rule

x0ddot = sym('x0ddot', [3*p,1], 'real');
zeta0_orient = sym('zeta0_orient', [3*p,1], 'real');

%% Definition of matrices to simplify later expressions

I0 = zeros(size(I));    % rotational moment of inertia matrices of the module blocks expressed in the global coordinate frame, vertically concatenated
M = zeros(6*p, 6*p);  % inertia matrix

I0 = sym(I0);   % convert to symbolic
M = sym(M);     % convert to symbolic
for i = 1:p
   x0i = x0(1+6*(i-1) : 6*i, 1);      % state of the ith module
   Roti = Reuler(x0i(4:6,1));        % rotation between frame zero and ith frame
   
   I0(3*(i-1)+1 : 3*i, :) = Roti * I(3*(i-1)+1 : 3*i, :) * Roti';
   M(6*(i-1)+1 : 6*i,  6*(i-1)+1 : 6*i) = [m(i)*eye(3), zeros(3,3); zeros(3,3), I0(3*(i-1)+1 : 3*i, :)]; 
end

% Selection matrix to isolate z-component of position vector
sz = [0 0 1 0 0 0];
Sz = kron(eye(p), sz);

%% Define the Lagrangian

% kinetic energy
KE = (1/2) * x0dot' * M * x0dot;

% potential energy
% x = x02x_sym(x0, params);
Lpe = zeros(p,1);
for i = 1:p
%     xi = x(1+6*(i-1) : 6*i, 1);      % state of the ith module
    Lpe(i,1) = sum(L(1:i)); 
%     h(i,1) = xi(3) - L(i);
%     z0(i,1) = L(i)/2 * (cos(xi(5)) * cos(xi(4)) - 1);
end
% PE = -g * m' * h;
PE = g * m' * (Lpe - Sz*x0);
% PE = -g * m' * Sz*x0;
% PE = 0;

% lagrangian
Lagrangian = KE - PE;

%% Euler Lagrange Equations of Motion (using external Matlab function) DID NOT WORK
% 
% % combination of x0 and x0dot of the form [x1, x1dot, x2, x2dot, ...]'
% X0 = zeros(2*length(x0),1);
% X0 = sym(X0);       % convert to symbolic
% X0(1:2:end) = x0;
% X0(2:2:end) = x0dot;
% X0 = sym2cell(X0);
% 
% Q_i = zeta0;      % internal forces, control forces
% Q_i = sym2cell(Q_i);
% Q_e = zeros(size(x0));      % external forces, loads
% Q_e = num2cell(Q_e);
% R = 0;  % friction term
% par = {m(1), I(1,1), I(1,2), I(1,3), I(2,1), I(2,2), I(2,3), I(3,1), I(3,2), I(3,3)};    % system parameters
% 
% VF = EulerLagrange(Lagrangian, X0, Q_i, Q_e, R, par, 'm', 'EOM');


%% Euler Lagrange Equations of Motion (using my own code)

% % rename stuff to make it work
% x0_ = sym('x0_', [3*p,1], 'real');          
% x0dot_ = sym('x0dot_', [3*p,1], 'real');
% 
% % dLdx0dot = (jacobian(Lagrangian, x0dot_orient) * dx0_orientdx0)';
% % dLdx0 = (jacobian(Lagrangian, x0_orient) * dx0_orientdx0)';
dLdx0dot = jacobian(Lagrangian, x0dot_orient)';
dLdx0 = jacobian(Lagrangian, x0_orient)';
% for i = 1 : p
%     dLdx0dot = subs(dLdx0dot, x0dot_orient(3*(i-1)+1 : 3*i) , x0dot_(6*(i-1)+4 : 6*i));
%     dLdx0 = subs(dLdx0, x0dot_orient(3*(i-1)+1 : 3*i) , x0dot_(6*(i-1)+4 : 6*i));
%     
%     dLdx0dot = subs(dLdx0dot, x0_orient(3*(i-1)+1 : 3*i) , x0_(6*(i-1)+4 : 6*i));
%     dLdx0 = subs(dLdx0, x0_orient(3*(i-1)+1 : 3*i) , x0_(6*(i-1)+4 : 6*i));
% end


% define x0 and x0dot as functions of time
syms t
x0t = zeros(3*p,1);
x0t = sym(x0t);
x0tdot = zeros(3*p,1);
x0tdot = sym(x0tdot);
for j = 1 : 3*p
   jstr = num2str(j);
   
   x0t(j) = sym(strcat('x0t', jstr, '(t)'));
   x0tdot(j) = sym(strcat('x0tdot', jstr, '(t)'));
end

dLdx0dot_t = subs(dLdx0dot, [x0_orient, x0dot_orient], [x0t, x0tdot]);
dLdx0_t = subs(dLdx0, [x0_orient, x0dot_orient], [x0t, x0tdot]);

% Euler Lagrange Equations of motion
EOM_raw = diff(dLdx0dot_t, t) - dLdx0 - zeta0_orient;    % assuming no load on the system 

% Character substitutions to get rid of all the 'diff(x(t), t)' stuff in EOM_raw
Dx0t = sym( zeros(3*p,1) );        % x0dot written in gross way, e.g. x0dot = diff(x0(t), t)
Dx0tdot = sym( zeros(3*p,1) );     % x0ddot written in gross way, e.g. x0ddot = diff(x0dot(t), t)
for i = 1:3*p
   istr = num2str(i);
   Dx0t(i,1) = sym(strcat( 'diff(x0t', istr, '(t), t)' )); 
   Dx0tdot(i,1) = sym(strcat( 'diff(x0tdot', istr, '(t), t)' )); 
end

EOM = subs(EOM_raw, [x0t; x0tdot; Dx0t; Dx0tdot], [x0_orient; x0dot_orient; x0dot_orient; x0ddot]);      % replace all instances of 't'

%% Creates Matlab function for evaluating the Equations of Motion
X0 = [x0_orient; x0dot_orient];       % dynamics state vector, x0 and x0dot vertically concatenated
X0dot = [x0dot_orient; x0ddot];
matlabFunction(EOM, 'File', 'EOM', 'Vars', {X0, X0dot, zeta0_orient, m, I}, 'Optimize', false);


%% Reactivate previously suppressed warning
warning('on', warn_id);
end










