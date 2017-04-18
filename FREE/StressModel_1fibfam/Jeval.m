function Jsym = Jeval(in1,in2,t0,E,G,in6,u)
%JEVAL
%    JSYM = JEVAL(IN1,IN2,T0,E,G,IN6,U)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    18-Apr-2017 17:40:17

L = in1(:,4);
L0 = in2(:,4);
P = in1(:,1);
T = in1(:,6);
gama = in1(:,2);
gama0 = in2(:,2);
phi = in1(:,5);
r = in1(:,3);
r0 = in2(:,3);
t2 = r.^2;
t3 = cot(gama);
t4 = t3.^2;
t5 = t4+1.0;
t6 = 1.0./r0;
t7 = r-r0;
t8 = r0.*t0.*2.0;
t9 = t0.^2;
t10 = t2+t8+t9;
t11 = sqrt(t10);
t12 = r-t11;
t13 = sin(gama);
t14 = 1.0./sqrt(t10);
t15 = r.*t14;
t16 = t15-1.0;
t17 = 1.0./L0;
t18 = cos(gama);
t19 = t12.^2;
t30 = r.*t12.*2.0;
t20 = t19-t30;
t21 = 1.0./L;
t22 = phi.*r.*t21;
t23 = atan(t22);
t24 = r.*2.0;
t25 = t11.*2.0;
t26 = t24-t25-r.*t16.*2.0+t12.*t16.*2.0;
t27 = r.*(1.0./2.0);
t28 = t11.*(1.0./2.0);
t29 = t27+t28;
t31 = 1.0./L.^2;
t32 = phi.^2;
t33 = t2.*t31.*t32;
t34 = t33+1.0;
t35 = 1.0./t34;
t36 = tan(gama0);
t37 = phi-L0.*t6.*t36;
t38 = 1.0./t13;
t39 = tan(gama);
t40 = 1.0./r;
Jsym = reshape([t2.*t3.*pi.*2.0,t2.*pi,0.0,-1.0,0.0,0.0,T.*t18.*-2.0-P.*t2.*t5.*pi.*2.0-E.*r.*t5.*t6.*t7.*t12.*pi.*2.0,T.*t13.*2.0,T.*t11.*t18.*2.0,0.0,L.*t13.*1.0./t18.^2-r.*1.0./t13.^2.*t18.*t37,-L.*t40.*(t39.^2+1.0),P.*r.*t3.*pi.*4.0+E.*r.*t3.*t6.*t12.*pi.*2.0+E.*t3.*t6.*t7.*t12.*pi.*2.0-E.*r.*t3.*t6.*t7.*t16.*pi.*2.0,P.*r.*pi.*2.0+E.*t17.*t26.*pi.*(L-L0),T.*r.*t13.*t14.*2.0+G.*t23.*t26.*t29.*pi-G.*t20.*t23.*pi.*(r.*t14.*(1.0./2.0)+1.0./2.0)-G.*phi.*t20.*t21.*t29.*t35.*pi,0.0,t37.*t38,L.*1.0./r.^2.*t39,0.0,-E.*t17.*t20.*pi,G.*phi.*r.*t20.*t29.*t31.*t35.*pi,0.0,1.0./t18,-t39.*t40,0.0,0.0,-G.*r.*t20.*t21.*t29.*t35.*pi,0.0,r.*t38,-1.0,t13.*-2.0,t18.*-2.0,t11.*t13.*2.0,0.0,0.0,0.0],[6,6]);
