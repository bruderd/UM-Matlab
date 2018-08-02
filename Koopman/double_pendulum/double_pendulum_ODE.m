function xdot = double_pendulum_ODE(t,x,u,params)
% DOUBLE_PENDULUM_ODE Ordinary differential equations for double pendulum.
%
%   author:  Alexander Erlich (alexander.erlich@gmail.com)
%
%   parameters:
%
%   t       Column vector of time points 
%   xdot    Solution array. Each row in xdot corresponds to the solution at a
%           time returned in the corresponding row of t.
%
%   This function calls is called by double_pendulum.
%
%   ---------------------------------------------------------------------

g=params.g; m1=params.m1; m2=params.m2; l1=params.l1; l2=params.l2;

xdot=zeros(4,1);


xdot(1)=x(3);

xdot(2)=x(4);

xdot(3)=-((g*(2*m1+m2)*sin(x(1))+m2*(g*sin(x(1)-2*x(2))+2*(l2*x(4)^2+...
    l1*x(3)^2*cos(x(1)-x(2)))*sin(x(1)-x(2))))/...
    (2*l1*(m1+m2-m2*cos(x(1)-x(2))^2))) + ...
    -1*x(3) + u(1);   % damping and input torque

xdot(4)=(((m1+m2)*(l1*x(3)^2+g*cos(x(1)))+l2*m2*x(4)^2*cos(x(1)-x(2)))*...
    sin(x(1)-x(2)))/(l2*(m1+m2-m2*cos(x(1)-x(2))^2)) + ...
    -1*(x(4) - x(3));   % damping







% xdot(1)=x(2);
% 
% xdot(2)=-((g*(2*m1+m2)*sin(x(1))+m2*(g*sin(x(1)-2*x(3))+2*(l2*x(4)^2+...
%     l1*x(2)^2*cos(x(1)-x(3)))*sin(x(1)-x(3))))/...
%     (2*l1*(m1+m2-m2*cos(x(1)-x(3))^2))) + ...
%     -1*x(2) + u(1);   % damping and input torque
% 
% xdot(3)=x(4);
% 
% xdot(4)=(((m1+m2)*(l1*x(2)^2+g*cos(x(1)))+l2*m2*x(4)^2*cos(x(1)-x(3)))*...
%     sin(x(1)-x(3)))/(l2*(m1+m2-m2*cos(x(1)-x(3))^2)) + ...
%     -1*(x(4) - x(2));   % damping