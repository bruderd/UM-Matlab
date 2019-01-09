function xdot = manip_jang_stiff_2000s_dynamics(in1,in2)
%MANIP_JANG_STIFF_2000S_DYNAMICS
%    XDOT = MANIP_JANG_STIFF_2000S_DYNAMICS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    04-Jan-2019 09:30:55

u1 = in2(1,:);
u2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t3 = x1-x2;
t2 = cos(t3);
t4 = sin(t3);
t5 = t2.^2;
t6 = x4.^2;
t7 = x3.^2;
xdot = [x3;x4;u1.*5.0e3-x1.*5.0e3-x3.*1.76e2+(sin(x1-x2.*2.0).*(9.81e2./1.0e2)+sin(x1).*2.943e1+t4.*(t6.*2.0+t2.*t7.*2.0))./(t5.*2.0-4.0);u2.*5.0e3+x1.*5.0e3-x2.*5.0e3+x3.*1.75e2-x4.*1.75e2-(t4.*(t7.*2.0+cos(x1).*(9.81e2./5.0e1)+t2.*t6))./(t5-2.0)];
