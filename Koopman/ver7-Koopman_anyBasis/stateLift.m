function polyBasis = stateLift(in1,u1)
%STATELIFT
%    POLYBASIS = STATELIFT(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    06-Nov-2019 15:20:06

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t2 = u1.^2;
t3 = x1.^2;
t4 = x2.^2;
t5 = x3.^2;
polyBasis = [1.0;x1;x2;x3;u1;t3;x1.*x2;t4;x1.*x3;x2.*x3;t5;u1.*x1;u1.*x2;u1.*x3;t2;x1.^3;t3.*x2;t4.*x1;x2.^3;t3.*x3;x1.*x2.*x3;t4.*x3;t5.*x1;t5.*x2;x3.^3;t3.*u1;u1.*x1.*x2;t4.*u1;u1.*x1.*x3;u1.*x2.*x3;t5.*u1;t2.*x1;t2.*x2;t2.*x3;u1.^3];
