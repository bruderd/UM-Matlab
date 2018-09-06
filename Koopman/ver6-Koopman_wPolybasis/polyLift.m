function polyBasis = polyLift(in1,in2)
%POLYLIFT
%    POLYBASIS = POLYLIFT(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    23-Aug-2018 16:01:54

u1 = in2(1,:);
u2 = in2(2,:);
u3 = in2(3,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
t2 = x1.^2;
t3 = x2.^2;
t4 = x3.^2;
t5 = x4.^2;
t6 = x5.^2;
t7 = x6.^2;
t8 = u1.^2;
t9 = u2.^2;
t10 = u3.^2;
polyBasis = [1.0;x1;x2;x3;x4;x5;x6;u1;u2;u3;t2;x1.*x2;t3;x1.*x3;x2.*x3;t4;x1.*x4;x2.*x4;x3.*x4;t5;x1.*x5;x2.*x5;x3.*x5;x4.*x5;t6;x1.*x6;x2.*x6;x3.*x6;x4.*x6;x5.*x6;t7;u1.*x1;u1.*x2;u1.*x3;u1.*x4;u1.*x5;u1.*x6;t8;u2.*x1;u2.*x2;u2.*x3;u2.*x4;u2.*x5;u2.*x6;u1.*u2;t9;u3.*x1;u3.*x2;u3.*x3;u3.*x4;u3.*x5;u3.*x6;u1.*u3;u2.*u3;t10;t2.*x1;t2.*x2;t3.*x1;t3.*x2;t2.*x3;x1.*x2.*x3;t3.*x3;t4.*x1;t4.*x2;t4.*x3;t2.*x4;x1.*x2.*x4;t3.*x4;x1.*x3.*x4;x2.*x3.*x4;t4.*x4;t5.*x1;t5.*x2;t5.*x3;t5.*x4;t2.*x5;x1.*x2.*x5;t3.*x5;x1.*x3.*x5;x2.*x3.*x5;t4.*x5;x1.*x4.*x5;x2.*x4.*x5;x3.*x4.*x5;t5.*x5;t6.*x1;t6.*x2;t6.*x3;t6.*x4;t6.*x5;t2.*x6;x1.*x2.*x6;t3.*x6;x1.*x3.*x6;x2.*x3.*x6;t4.*x6;x1.*x4.*x6;x2.*x4.*x6;x3.*x4.*x6;t5.*x6;x1.*x5.*x6;x2.*x5.*x6;x3.*x5.*x6;x4.*x5.*x6;t6.*x6;t7.*x1;t7.*x2;t7.*x3;t7.*x4;t7.*x5;t7.*x6;t2.*u1;u1.*x1.*x2;t3.*u1;u1.*x1.*x3;u1.*x2.*x3;t4.*u1;u1.*x1.*x4;u1.*x2.*x4;u1.*x3.*x4;t5.*u1;u1.*x1.*x5;u1.*x2.*x5;u1.*x3.*x5;u1.*x4.*x5;t6.*u1;u1.*x1.*x6;u1.*x2.*x6;u1.*x3.*x6;u1.*x4.*x6;u1.*x5.*x6;t7.*u1;t8.*x1;t8.*x2;t8.*x3;t8.*x4;t8.*x5;t8.*x6;t8.*u1;t2.*u2;u2.*x1.*x2;t3.*u2;u2.*x1.*x3;u2.*x2.*x3;t4.*u2;u2.*x1.*x4;u2.*x2.*x4;u2.*x3.*x4;t5.*u2;u2.*x1.*x5;u2.*x2.*x5;u2.*x3.*x5;u2.*x4.*x5;t6.*u2;u2.*x1.*x6;u2.*x2.*x6;u2.*x3.*x6;u2.*x4.*x6;u2.*x5.*x6;t7.*u2;u1.*u2.*x1;u1.*u2.*x2;u1.*u2.*x3;u1.*u2.*x4;u1.*u2.*x5;u1.*u2.*x6;t8.*u2;t9.*x1;t9.*x2;t9.*x3;t9.*x4;t9.*x5;t9.*x6;t9.*u1;t9.*u2;t2.*u3;u3.*x1.*x2;t3.*u3;u3.*x1.*x3;u3.*x2.*x3;t4.*u3;u3.*x1.*x4;u3.*x2.*x4;u3.*x3.*x4;t5.*u3;u3.*x1.*x5;u3.*x2.*x5;u3.*x3.*x5;u3.*x4.*x5;t6.*u3;u3.*x1.*x6;u3.*x2.*x6;u3.*x3.*x6;u3.*x4.*x6;u3.*x5.*x6;t7.*u3;u1.*u3.*x1;u1.*u3.*x2;u1.*u3.*x3;u1.*u3.*x4;u1.*u3.*x5;u1.*u3.*x6;t8.*u3;u2.*u3.*x1;u2.*u3.*x2;u2.*u3.*x3;u2.*u3.*x4;u2.*u3.*x5;u2.*u3.*x6;u1.*u2.*u3;t9.*u3;t10.*x1;t10.*x2;t10.*x3;t10.*x4;t10.*x5;t10.*x6;t10.*u1;t10.*u2;t10.*u3];