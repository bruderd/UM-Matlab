% main: Run the Koopman sysid process on a user defined nonlinear system

%% 


%% USER EDIT SECTION

% define parameters
params.struct

params.n = 2;   % dimension of state space
params.maxDegree = 4;   % maximum degree of monomial basis
params = def_polyLift(params);  % creates the lifting function, polyLift


% define dynamics
x = sym('x', [params.n, 1]);

vf = [x(2);...
     -x(1)];
 
matlabFunction(vf, 'File', 'vf_real', 'Vars', {x}); % create plant dynamics func.

