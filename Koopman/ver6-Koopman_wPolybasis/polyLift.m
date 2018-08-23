function polyBasis = polyLift(in1,in2)
%POLYLIFT
%    POLYBASIS = POLYLIFT(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    22-Aug-2018 15:54:45

u1 = in2(1,:);
u2 = in2(2,:);
u3 = in2(3,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
polyBasis = [1.0;x1;x2;x3;x4;x5;x6;u1;u2;u3];
