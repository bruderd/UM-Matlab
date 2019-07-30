function EOM = set_EOM(params)
%setEOM: Find symbolic expression of the equations of motion
%   Also saves a function for evaluating the equations of motion

%% define symbolic variables
alpha = sym('alpha', [params.Nlinks,1], 'real');
alphadot = sym('alphadot', [params.Nlinks,1], 'real');
alphaddot = sym('alphaddot', [params.Nlinks,1], 'real');

theta = alpha2theta(alpha);
thetadot = alpha2theta(alphadot);

[ x , xcm ]= alpha2x(alpha, params);

%% define Jacobians

J_theta_alpha = jacobian( theta , alpha );

J_xcm_alpha = jacobian( xcm , alpha );
xcmdot = J_xcm_alpha * alphadot;    % velocity of link COMs 

%% define useful matrices

% mass matrix
M = eye(params.Nlinks) * params.m;
I = eye(params.Nlinks) * params.i;
K = ones(1 , params.Nlinks) * params.k;
D = eye(params.Nlinks) * params.d;

%% define Lagrangian (L = KE - PE)

% mass matrix
Dq = params.m * J_xcm_alpha' * J_xcm_alpha + params.i * J_theta_alpha' * J_theta_alpha;

% kinetic energy
KE = (1/2) * alphadot' * Dq * alphadot;

% potential energy (needs minus sign since "down" is positive)
PE = - params.m * params.g * ones(1 , length(xcm)/2) * xcm(2:2:end) + ...
     (1/2) * alpha' * params.k * alpha;

Lagrangian = KE - PE;

%% derive equations of motion

% save mass matrix as a function
EOM.fcns.get_massMatrix = matlabFunction(Dq, 'Vars', { alpha }, 'Optimize', false);

% derive non-inertial part of dynamics
% creata a variable alpha that is a function of t
syms t
alpha_t = zeros( params.Nlinks , 1 );
alpha_t = sym(alpha_t); 
for i = 1 : params.Nlinks
    istr = num2str(i);
    alpha_t(i) = str2sym(strcat('alpha_t', istr, '(t)'));
end

% write mass matrix as a function of t
Dq_t = subs( Dq , alpha , alpha_t );

% differentiate mass matrix wrt t
Dq_dt = diff( Dq_t , t );

% character substitutions to get rid of all the 'diff(x(t), t)' stuff
alpha_dt = zeros( params.Nlinks , 1 );
alpha_dt = sym(alpha_dt); 
for i = 1 : params.Nlinks
    istr = num2str(i);
    alpha_dt(i) = str2sym(strcat( 'diff(alpha_t', istr, '(t), t)' )); 
end
Dq_dt = subs( Dq_dt , [ alpha_t , alpha_dt ] , [ alpha , alphadot ] ); % replace all t's

dLdalpha = jacobian(Lagrangian, alpha)';

% include damping and input terms
% damping
damp = params.d * alphadot;
EOM.fcns.get_damp = matlabFunction(damp, 'Vars', { alphadot }, 'Optimize', false);

% input
u = sym('u', [params.Nmods,1], 'real'); % input. Desired joint angle for all joints in each module
input = -params.ku * ( kron( u , ones( params.nlinks , 1) ) - alpha );   % vector of all joint torques
EOM.fcns.get_input = matlabFunction(input, 'Vars', { alpha , u }, 'Optimize', false);

% save damping and input as a function
dampNinput = damp + input;
EOM.fcns.get_dampNinput = matlabFunction(dampNinput, 'Vars', { alpha , alphadot , u }, 'Optimize', false);

% save non-inertial part of dynamics as a function
nonInert = Dq_dt * alphadot - dLdalpha + damp + input;
EOM.fcns.get_nonInert = matlabFunction(nonInert, 'Vars', { alpha , alphadot , u }, 'Optimize', false);


%% set output
EOM.massMatrix = Dq;
EOM.nonInert = nonInert;
EOM.damp = damp;
EOM.input = input;

%% save the equations of motion for this system

%% Symbolically derive Euler-Lagrange equations of motion (takes too long)
% COMMENT EVERYTHING BELOW HERE FOR THE NUMERIC VERSION

% EOM = struct;   % equations of motion struct
% 
% dLdalphadot = jacobian(Lagrangian, alphadot)';
% dLdalpha = jacobian(Lagrangian, alpha)';
% 
% % define x0 and x0dot as functions of time
% syms t
% x0t = zeros(params.Nlinks,1);
% x0t = sym(x0t);
% x0tdot = zeros(params.Nlinks,1);
% x0tdot = sym(x0tdot);
% for j = 1 : params.Nlinks
%    jstr = num2str(j);
%    
%    x0t(j) = str2sym(strcat('x0t', jstr, '(t)'));
%    x0tdot(j) = str2sym(strcat('x0tdot', jstr, '(t)'));
% end
% 
% dLdx0dot_t = subs(dLdalphadot, [alpha, alphadot], [x0t, x0tdot]);
% dLdx0_t = subs(dLdalpha, [alpha, alphadot], [x0t, x0tdot]);
% 
% % Euler Lagrange Equations of motion
% EOM_raw = diff(dLdx0dot_t, t) - dLdalpha;    % assuming no load on the system 
% 
% % Character substitutions to get rid of all the 'diff(x(t), t)' stuff in EOM_raw
% Dx0t = sym( zeros(params.Nlinks,1) );        % x0dot written in gross way, e.g. x0dot = diff(x0(t), t)
% Dx0tdot = sym( zeros(params.Nlinks,1) );     % x0ddot written in gross way, e.g. x0ddot = diff(x0dot(t), t)
% for i = 1:params.Nlinks
%    istr = num2str(i);
%    Dx0t(i,1) = str2sym(strcat( 'diff(x0t', istr, '(t), t)' )); 
%    Dx0tdot(i,1) = str2sym(strcat( 'diff(x0tdot', istr, '(t), t)' )); 
% end
% 
% EOM.auton.implicit = subs(EOM_raw, [x0t; x0tdot; Dx0t; Dx0tdot], [alpha; alphadot; alphadot; alphaddot]);      % replace all instances of 't'
% 
% addot = solve( EOM.auton.implicit , alphaddot );   % convert EOM to explicit form, i.e. xddot = f(x,xdot)
% 
% addot_cell = struct2cell(addot);
% vf_auton = [ alphadot ; zeros( params.Nlinks , 1 ) ];   % converts to vector field form, i.e. Xdot = f(X), where X = [x,xdot]'. No inputs or damping
% for i = 1 : params.Nlinks
%     vf_auton(params.Nlinks + i) = addot_cell{i};
% end
% EOM.auton.vf = vf_auton;
% 
% %% Add damping and input into the dynamics
% 
% % damping
% damp = params.d * alphadot;
% 
% % input
% u = sym('u', [params.Nmods,1], 'real'); % input. Applies a torque at all joints in each module
% input = kron( u , ones( params.nlinks , 1) );   % vector of all joint torques
% 
% EOM.wInput.vf = EOM.auton.vf + [ zeros(size(alphadot)) ; damp] + [ zeros(size(alphadot)) ; input ];
% 
% 
% %% Creates Matlab function for evaluating the Equations of Motion
% 
% X0 = [alpha; alphadot];       % dynamics state vector, x0 and x0dot vertically concatenated
% X0dot = [alphadot; alphaddot];
% matlabFunction(EOM.wInput.vf, 'File', 'vf_sym', 'Vars', { t , X0 , u }, 'Optimize', false);


end

