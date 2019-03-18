function [xdot] = armload_vf(x,u,params)
% arm_vf: Dynamic equations for a 1 dof arm with load attached at the end
% effector

% if x is symbolic, make xdot symbolic
if isa(x,'sym')
    xdot = sym( 'xdot' , [params.n,1] );
else
    xdot=zeros(params.n,1);
end


xdot(1,:) = x(2);
xdot(2,:) = ( 1 / ( ( (1/3)*params.M + x(3) ) * params.L^2 ) ) * ( x(3)*params.L*params.g*cos(x(1)) + params.K*(u-x(1)) + params.D*x(2) );
xdot(3,:) = 0;


% if input is symbolic, create matlab function that evaluates these
% dynamics with the given parameters.
if isa( x , 'sym' )
    name = [ 'simDynamics' , filesep , params.name , '_dynamics.m' ];
    matlabFunction( xdot, 'File', name, 'Vars', {x , u} );
end

end