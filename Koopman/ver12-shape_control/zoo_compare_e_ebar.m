% zoo_compare_e_ebar
%
% Solves for nu over one time step. Then solves for u given the structure
% of e is known. Then reconstructs e using that value of u, which we call
% ebar. Then compares e to ebar. Hopefully they are very similar, because
% this implies that mpc chooses close to a "valid" nu over one time step.

function [ e , ebar ] = zoo_compare_e_ebar( x_now , u_last , y_des , sys )

% lift the state
zx_now = sys.lift.zx( x_now );

% solve for nu over one timestep
nu = zoo_check_nu1step_2( zx_now , y_des , sys );
e = sys.model.Beta_pinv * nu + sys.traindata.emean';

% solve for x one timestep in the future given nu
zx_plus = sys.model.A * zx_now + e;
x_plus = sys.model.C * zx_plus;

% solve for u based on nu and structure of e
u_now = zoo_solve4u( nu , x_now , u_last , sys );

Febar = matlabFunction( sys.model.B * sys.basis.zu , 'Vars' , { sys.params.x , sys.params.u } );
ebar = Febar( x_now , u_now );

end