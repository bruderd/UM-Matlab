function setJacobians(params)
%setJacobians: Defines the Jacobians of the system symbolically then
%creates Matlab functions to evaluate them.
%   Detailed explanation goes here

p = params.p;

%% Connect Block Jacobian

alpha = sym('alpha', [p,1], 'real');
q = alpha2q_sym(alpha, params);

J_aq = jacobian(q, alpha);

% Create Matlab function for evaluating J_aq
matlabFunction(J_aq, 'File', 'J_aq', 'Vars', {alpha});


%% Spine Block Jacobian

x = alpha2x_sym_exact(alpha, params);

J_ax = jacobian(x, alpha);

% Create Matlab function for evaluating J_ax
matlabFunction(J_ax, 'File', 'J_ax', 'Vars', {alpha});


end

