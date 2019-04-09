function xdot = larmload_100s_load10_dynamics(in1,u1)
%LARMLOAD_100S_LOAD10_DYNAMICS
%    XDOT = LARMLOAD_100S_LOAD10_DYNAMICS(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    01-Apr-2019 17:44:57

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
xdot = [x2;(u1.*2.0e1-x1.*2.0e1-x2.*1.0e1+x3.*cos(x1).*(9.81e2./5.0e1))./(x3.*4.0+4.0./3.0);0.0];
