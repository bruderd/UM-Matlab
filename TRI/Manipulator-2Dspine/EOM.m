function EOM = EOM(in1,in2,in3,in4)
%EOM
%    EOM = EOM(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    20-Nov-2017 21:07:20

alpha2 = in1(2,:);
alpha3 = in1(3,:);
alpha4 = in1(4,:);
alpha5 = in1(5,:);
alpha6 = in1(6,:);
alphadot2 = in1(8,:);
alphadot3 = in1(9,:);
alphadot4 = in1(10,:);
alphadot5 = in1(11,:);
alphadot6 = in1(12,:);
alphaddot1 = in2(7,:);
alphaddot2 = in2(8,:);
alphaddot3 = in2(9,:);
alphaddot4 = in2(10,:);
alphaddot5 = in2(11,:);
alphaddot6 = in2(12,:);
gama1 = in4(1,:);
gama2 = in4(2,:);
gama3 = in4(3,:);
gama4 = in4(4,:);
gama5 = in4(5,:);
gama6 = in4(6,:);
zeta1 = in3(1,:);
zeta2 = in3(2,:);
zeta3 = in3(3,:);
zeta4 = in3(4,:);
zeta5 = in3(5,:);
zeta6 = in3(6,:);
EOM = [alphaddot1.*1.0e-4-gama1-zeta1;alpha2.*1.22625e-3+alpha3.*8.175e-4+alpha4.*4.905e-4+alpha5.*2.4525e-4+alpha6.*8.175e-5+alphaddot2.*5.01025e-2+alphaddot3.*4.500194444444444e-2+alphaddot4.*4.000140277777778e-2+alphaddot5.*3.500088888888889e-2+alphaddot6.*3.000041666666667e-2-gama2-zeta2+(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*8.333333333333333e-6+alpha3.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)))+(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot4.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot5.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot4.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot5.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*1.666666666666667e-5+alpha3.*1.25e-5+alpha4.*8.333333333333333e-6+alpha5.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)))+alpha2.*alphadot2.^2.*1.388888888888889e-8+alpha2.^2.*alphaddot2.*1.388888888888889e-8+(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot4.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot4.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*1.25e-5+alpha3.*8.333333333333333e-6+alpha4.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)))+(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2))+(alpha2.*1.25e-4+alpha3.*1.0e-4+alpha4.*7.5e-5+alpha5.*5.0e-5+alpha6.*2.5e-5).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)));alpha2.*8.175e-4+alpha3.*8.175e-4+alpha4.*4.905e-4+alpha5.*2.4525e-4+alpha6.*8.175e-5+alphaddot2.*4.500194444444444e-2+alphaddot3.*4.510152777777778e-2+alphaddot4.*4.000111111111111e-2+alphaddot5.*3.500070833333333e-2+alphaddot6.*3.000033333333333e-2-gama3-zeta3+(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*4.166666666666667e-6+alpha3.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)))+(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot4.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot5.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot4.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot5.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*1.25e-5+alpha3.*1.25e-5+alpha4.*8.333333333333333e-6+alpha5.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)))+(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot4.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot4.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*8.333333333333333e-6+alpha3.*8.333333333333333e-6+alpha4.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)))+(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2))+(alpha2.*1.0e-4+alpha3.*1.0e-4+alpha4.*7.5e-5+alpha5.*5.0e-5+alpha6.*2.5e-5).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)));alpha2.*4.905e-4+alpha3.*4.905e-4+alpha4.*4.905e-4+alpha5.*2.4525e-4+alpha6.*8.175e-5+alphaddot2.*4.000140277777778e-2+alphaddot3.*4.000111111111111e-2+alphaddot4.*4.010081944444444e-2+alphaddot5.*3.500052777777778e-2+alphaddot6.*3.000025e-2-gama4-zeta4+(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot4.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot5.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot4.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot5.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*8.333333333333333e-6+alpha3.*8.333333333333333e-6+alpha4.*8.333333333333333e-6+alpha5.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)))+(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot4.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot4.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*4.166666666666667e-6+alpha3.*4.166666666666667e-6+alpha4.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)))+(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2))+(alpha2.*7.5e-5+alpha3.*7.5e-5+alpha4.*7.5e-5+alpha5.*5.0e-5+alpha6.*2.5e-5).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)));alpha2.*2.4525e-4+alpha3.*2.4525e-4+alpha4.*2.4525e-4+alpha5.*2.4525e-4+alpha6.*8.175e-5+alphaddot2.*3.500088888888889e-2+alphaddot3.*3.500070833333333e-2+alphaddot4.*3.500052777777778e-2+alphaddot5.*3.510034722222222e-2+alphaddot6.*3.000016666666667e-2-gama5-zeta5+(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot3.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot4.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphaddot5.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot2.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot3.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot4.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2)+alphadot5.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)).*(1.0./4.0e2))+(alpha2.*4.166666666666667e-6+alpha3.*4.166666666666667e-6+alpha4.*4.166666666666667e-6+alpha5.*4.166666666666667e-6).*(alphaddot2.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)))+(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2))+(alpha2.*5.0e-5+alpha3.*5.0e-5+alpha4.*5.0e-5+alpha5.*5.0e-5+alpha6.*2.5e-5).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)));alpha2.*8.175e-5+alpha3.*8.175e-5+alpha4.*8.175e-5+alpha5.*8.175e-5+alpha6.*8.175e-5+alphaddot2.*3.000041666666667e-2+alphaddot3.*3.000033333333333e-2+alphaddot4.*3.000025e-2+alphaddot5.*3.000016666666667e-2+alphaddot6.*3.300008333333333e-2-gama6-zeta6+(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2)).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2)+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)).*(3.0./2.0e2))+(alpha2.*2.5e-5+alpha3.*2.5e-5+alpha4.*2.5e-5+alpha5.*2.5e-5+alpha6.*2.5e-5).*(alphaddot2.*(alpha2.*(1.0./1.2e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot3.*(alpha2.*(1.0./1.5e2)+alpha3.*(1.0./1.5e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot4.*(alpha2.*(1.0./2.0e2)+alpha3.*(1.0./2.0e2)+alpha4.*(1.0./2.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphaddot5.*(alpha2.*(1.0./3.0e2)+alpha3.*(1.0./3.0e2)+alpha4.*(1.0./3.0e2)+alpha5.*(1.0./3.0e2)+alpha6.*(1.0./6.0e2))+alphadot2.*(alphadot2.*(1.0./1.2e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot3.*(alphadot2.*(1.0./1.5e2)+alphadot3.*(1.0./1.5e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphaddot6.*(alpha2.*(1.0./6.0e2)+alpha3.*(1.0./6.0e2)+alpha4.*(1.0./6.0e2)+alpha5.*(1.0./6.0e2)+alpha6.*(1.0./6.0e2))+alphadot4.*(alphadot2.*(1.0./2.0e2)+alphadot3.*(1.0./2.0e2)+alphadot4.*(1.0./2.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot5.*(alphadot2.*(1.0./3.0e2)+alphadot3.*(1.0./3.0e2)+alphadot4.*(1.0./3.0e2)+alphadot5.*(1.0./3.0e2)+alphadot6.*(1.0./6.0e2))+alphadot6.*(alphadot2.*(1.0./6.0e2)+alphadot3.*(1.0./6.0e2)+alphadot4.*(1.0./6.0e2)+alphadot5.*(1.0./6.0e2)+alphadot6.*(1.0./6.0e2)))];