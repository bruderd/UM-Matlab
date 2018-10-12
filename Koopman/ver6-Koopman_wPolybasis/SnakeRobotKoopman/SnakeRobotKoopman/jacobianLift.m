function dlift = jacobianLift(in1,u1)
%JACOBIANLIFT
%    DLIFT = JACOBIANLIFT(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    12-Oct-2018 19:22:22

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t2 = x1.^2;
t3 = x1.*x2.*2.0;
t4 = x2.^2;
t5 = x1.*x3.*2.0;
t6 = x3.^2;
t7 = x2.*x3.*2.0;
t8 = u1.*x1;
t9 = u1.*x3;
t10 = u1.*x2;
t11 = u1.^2;
dlift = reshape([0.0,1.0,0.0,0.0,0.0,x1.*2.0,x2,0.0,x3,0.0,0.0,u1,0.0,0.0,0.0,t2.*3.0,t3,t4,0.0,t5,x2.*x3,0.0,t6,0.0,0.0,u1.*x1.*2.0,t10,0.0,t9,0.0,0.0,t11,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,x1,x2.*2.0,0.0,x3,0.0,0.0,u1,0.0,0.0,0.0,t2,t3,t4.*3.0,0.0,x1.*x3,t7,0.0,t6,0.0,0.0,t8,u1.*x2.*2.0,0.0,t9,0.0,0.0,t11,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,x1,x2,x3.*2.0,0.0,0.0,u1,0.0,0.0,0.0,0.0,0.0,t2,x1.*x2,t4,t5,t7,t6.*3.0,0.0,0.0,0.0,t8,t10,u1.*x3.*2.0,0.0,0.0,t11,0.0],[35,3]);
