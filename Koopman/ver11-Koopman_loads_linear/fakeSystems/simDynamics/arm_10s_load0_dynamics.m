function xdot = arm_10s_load1_dynamics(in1,u1)
%ARM_10S_LOAD1_DYNAMICS
%    XDOT = ARM_10S_LOAD1_DYNAMICS(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    11-Mar-2019 16:15:02

x1 = in1(1,:);
x2 = in1(2,:);
xdot = [x2;u1.*6.0e1-x1.*6.0e1-x2.*3.0e1];
