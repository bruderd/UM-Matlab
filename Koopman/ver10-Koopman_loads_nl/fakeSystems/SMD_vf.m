function xdot = SMD_vf(x,u,params)
% DCmotor_vf: Ordinary differential equations for double pendulum.
%
%   author:  Alexander Erlich (alexander.erlich@gmail.com)
%
%   parameters:
%
%   t       Column vector of time points 
%   xdot    Solution array. Each row in xdot corresponds to the solution at a
%           time returned in the corresponding row of t.
%
%
%   ---------------------------------------------------------------------

m = params.m;
b = params.b;
k = params.k;

if isa(x,'sym')
    xdot = sym( 'xdot' , [params.n,1] );
else
    xdot = zeros(params.n,1);
end
xdot = [0, 1; -k/m, -b/m] * x + [0, 1]' * u;


% if input is symbolic, create matlab function that evaluates these
% dynamics with the given parameters.
if isa( x , 'sym' )
    name = [ 'simDynamics' , filesep , params.name , '_dynamics.m' ];
    matlabFunction( xdot, 'File', name, 'Vars', {x , u} );
end

end