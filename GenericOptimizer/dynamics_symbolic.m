% dynamics_symbolic
%   User defines their system dynamics sybolically here, then matlab
%   creates functions for them to be called by vf.m


function dynamics_symbolic(params)
%% USER SPECIFIED SECTION.

% define symbolic variables
syms V turn xs ys thetas dxs dys dthetas

% define state vector (x), its time derivative (xdot), and the input (u)
x = [xs, ys, thetas];
xdot = [dxs, dys, dthetas];
u = [V, turn];

% define (implicit) dynamics of the form f = 0 (f should be nx1)
f_sym = [V*cos(thetas) - dxs;...
         V*sin(thetas) - dys;...
         turn - dthetas];


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