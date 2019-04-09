function xdot = armloadv2_10s_load7p5_dynamics(in1,u1)
%ARMLOADV2_10S_LOAD7P5_DYNAMICS
%    XDOT = ARMLOADV2_10S_LOAD7P5_DYNAMICS(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    01-Apr-2019 19:12:13

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
xdot = [x2;-(u1.*-2.0e1+x1.*2.0e1+x2.*1.0e1+x3.*cos(x1).*(9.81e2./1.0e2))./(x3+1.0./3.0);0.0];
