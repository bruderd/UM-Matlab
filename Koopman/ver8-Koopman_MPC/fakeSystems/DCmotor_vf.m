function xdot = DCmotor_vf(x,u,params)
% vf_doublePendulum: Ordinary differential equations for double pendulum.
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

La = params.La;
Ra = params.Ra;
km = params.km;
J = params.J;
B = params.B;
tau1 = params.tau1;
ua = params.ua;

xdot = zeros(2,1);
xdot(1) = -(Ra/La) * x(1) - (km/La) * x(2) * u + ua/La;
xdot(2) = -(B/J) * x(2) + (km/J) * x(1) * u - tau1/J;

end