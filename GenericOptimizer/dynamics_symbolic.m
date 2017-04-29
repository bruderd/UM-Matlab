% dynamics_symbolic
%   User defines their system dynamics sybolically here, then matlab
%   creates functions for them to be called by vf.m


function dynamics_symbolic(params)
%% USER SPECIFIED SECTION.

% defining symbolic variables
syms a b c dadt dbdt

% define state vector (x) and its time derivative (xdot)
x = [a, b];
xdot = [dadt, dbdt];

% defing input (u)
u = c;

% define (implicit) dynamics of the form f = 0 (f should be nx1)
f_sym = [a + dadt;
    b + c];


%% DO NOT EDIT BELOW THIS LINE---------------------------------------------
% uses symbolic differentiation to find gradients
dfdx_sym = jacobian(f_sym, x);
dfdu_sym = jacobian(f_sym, u);
dfdxdot_sym = jacobian(f_sym, xdot);

% creates matlab functions for evaluating dynamics
matlabFunction(f_sym, 'File', 'f_eval', 'Vars', {x, u, xdot});
matlabFunction(dfdx_sym, 'File', 'dfdx_eval', 'Vars', {x, u, xdot});
matlabFunction(dfdu_sym, 'File', 'dfdu_eval', 'Vars', {x, u, xdot});
matlabFunction(dfdxdot_sym, 'File', 'dfdxdot_eval', 'Vars', {x, u, xdot});

end