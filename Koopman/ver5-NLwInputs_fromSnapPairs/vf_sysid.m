function vf2 = vf_sysid(in1,u1)
%VF_SYSID
%    VF2 = VF_SYSID(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    03-Aug-2018 14:25:01

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
vf2 = [x1.*(-1.33527615982945e-5)+x2.*8.904365768491434e-6+x3.*1.001465621845637-x4.*7.853103238297686e-4+8.33861698200869e-6;x1.*8.21256123991616e-6-x2.*2.2183552268823e-5-x3.*1.254475319906842e-3+x4.*1.001617715967966-5.548605522981327e-6;x1.*(-5.694321062669435e-2)+x2.*1.735328105567626e-2-x3.*2.979183823817346e-1-x4.*2.570818753734801e-1-1.279194409752649e-2;x1.*2.598091294969637e-3-x2.*8.222366100266959e-2+x3.*8.631213752058974e-1-x4.*2.192901236168455e-1+2.490131271077826e-2];
