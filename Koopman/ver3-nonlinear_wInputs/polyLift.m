function polyBasis = polyLift(in1,u1)
%POLYLIFT
%    POLYBASIS = POLYLIFT(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    31-Jul-2018 16:47:27

x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1.^2;
t3 = x2.^2;
t4 = u1.^2;
polyBasis = [1.0;x1;x2;u1;t2;x1.*x2;t3;u1.*x1;u1.*x2;t4;t2.*x1;t2.*x2;t3.*x1;t3.*x2;t2.*u1;u1.*x1.*x2;t3.*u1;t4.*x1;t4.*x2;t4.*u1];
